// Copyright 2017 10x Genomics

//! Containers for path-compressed De Bruijn graphs

use bit_set::BitSet;
use log::{debug, trace};
use serde_derive::{Deserialize, Serialize};
use smallvec::SmallVec;
use std::borrow::Borrow;
use std::cmp::min;
use std::collections::HashSet;
use std::collections::VecDeque;
use std::f32;
use std::fmt::{self, Debug};
use std::fs::File;
use std::hash::Hash;
use std::io::Error;
use std::io::Write;
use std::iter::FromIterator;
use std::marker::PhantomData;
use std::path::Path;

use boomphf::hashmap::BoomHashMap;

use serde_json;
use serde_json::Value;

type SmallVec4<T> = SmallVec<[T; 4]>;
type SmallVec8<T> = SmallVec<[T; 8]>;

use crate::compression::CompressionSpec;
use crate::dna_string::{DnaString, DnaStringSlice, PackedDnaStringSet};
use crate::Dir;
use crate::Exts;
use crate::Kmer;
use crate::Mer;
use crate::Vmer;

/// A compressed DeBruijn graph carrying auxiliary data on each node of type `D`.
/// This type does not carry the sorted index arrays the allow the graph
/// to be walked efficiently. The `DeBruijnGraph` type wraps this type and add those
/// vectors.
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct BaseGraph<K, D> {
    pub sequences: PackedDnaStringSet,
    pub exts: Vec<Exts>,
    pub data: Vec<D>,
    pub stranded: bool,
    phantom: PhantomData<K>,
}

impl<K, D> BaseGraph<K, D> {
    pub fn new(stranded: bool) -> Self {
        BaseGraph {
            sequences: PackedDnaStringSet::new(),
            exts: Vec::new(),
            data: Vec::new(),
            phantom: PhantomData,
            stranded,
        }
    }

    pub fn len(&self) -> usize {
        self.sequences.len()
    }

    pub fn combine<I: Iterator<Item = BaseGraph<K, D>>>(graphs: I) -> Self {
        let mut sequences = PackedDnaStringSet::new();
        let mut exts = Vec::new();
        let mut data = Vec::new();
        let mut stranded = Vec::new();

        for g in graphs {
            for s in 0..g.sequences.len() {
                sequences.add(&g.sequences.get(s));
            }

            exts.extend(g.exts);
            data.extend(g.data);
            stranded.push(g.stranded);
        }

        if !stranded.iter().all(|x| *x) && !stranded.iter().all(|x| !*x) {
            panic!("attempted to combine stranded and unstranded graphs");
        }

        let out_stranded = stranded.iter().all(|x| *x);

        BaseGraph {
            sequences,
            stranded: out_stranded,
            exts,
            data,
            phantom: PhantomData,
        }
    }
}

impl<K: Kmer, D> BaseGraph<K, D> {
    pub fn add<'b, R: Borrow<u8>, S: IntoIterator<Item = R>>(
        &mut self,
        sequence: S,
        exts: Exts,
        data: D,
    ) {
        self.sequences.add(sequence);
        self.exts.push(exts);
        self.data.push(data);
    }
}

impl<K: Kmer + Send + Sync, D> BaseGraph<K, D> {
    pub fn finish(self) -> DebruijnGraph<K, D> {
        let indices: Vec<u32> = (0..self.len() as u32).collect();

        let left_order = {
            let mut kmers: Vec<K> = Vec::with_capacity(self.len());
            for idx in &indices {
                kmers.push(self.sequences.get(*idx as usize).first_kmer());
            }
            BoomHashMap::new_parallel(kmers, indices.clone())
        };

        let right_order = {
            let mut kmers: Vec<K> = Vec::with_capacity(self.len());
            for idx in &indices {
                kmers.push(self.sequences.get(*idx as usize).last_kmer());
            }
            BoomHashMap::new_parallel(kmers, indices)
        };

        DebruijnGraph {
            base: self,
            left_order,
            right_order,
        }
    }
}

impl<K: Kmer, D> BaseGraph<K, D> {
    pub fn finish_serial(self) -> DebruijnGraph<K, D> {
        let indices: Vec<u32> = (0..self.len() as u32).collect();

        let left_order = {
            let mut kmers: Vec<K> = Vec::with_capacity(self.len());
            for idx in &indices {
                kmers.push(self.sequences.get(*idx as usize).first_kmer());
            }
            BoomHashMap::new(kmers, indices.clone())
        };

        let right_order = {
            let mut kmers: Vec<K> = Vec::with_capacity(self.len());
            for idx in &indices {
                kmers.push(self.sequences.get(*idx as usize).last_kmer());
            }
            BoomHashMap::new(kmers, indices)
        };

        DebruijnGraph {
            base: self,
            left_order,
            right_order,
        }
    }
}

/// A compressed DeBruijn graph carrying auxiliary data on each node of type `D`.
/// The struct carries sorted index arrays the allow the graph
/// to be walked efficiently.
#[derive(Serialize, Deserialize, Debug)]
pub struct DebruijnGraph<K: Hash, D> {
    pub base: BaseGraph<K, D>,
    left_order: BoomHashMap<K, u32>,
    right_order: BoomHashMap<K, u32>,
}

impl<K: Kmer, D: Debug> DebruijnGraph<K, D> {
    /// Total number of nodes in the DeBruijn graph
    pub fn len(&self) -> usize {
        self.base.len()
    }

    /// Get a node given it's `node_id`
    pub fn get_node<'a>(&'a self, node_id: usize) -> Node<'a, K, D> {
        Node {
            node_id: node_id,
            graph: self,
        }
    }

    /// Get a node given it's `node_id`
    pub fn get_node_kmer<'a>(&'a self, node_id: usize) -> NodeKmer<'a, K, D> {
        let node = self.get_node(node_id);
        let node_seq = node.sequence();

        NodeKmer {
            node_id: node_id,
            node_seq_slice: node_seq,
            phantom_d: PhantomData,
            phantom_k: PhantomData,
        }
    }

    /// Return an iterator over all nodes in the graph
    pub fn iter_nodes<'a>(&'a self) -> NodeIter<'a, K, D> {
        NodeIter {
            graph: self,
            node_id: 0,
        }
    }

    /// Find the edges leaving node `node_id` in direction `Dir`. Should generally be
    /// accessed via a Node wrapper object
    fn find_edges(&self, node_id: usize, dir: Dir) -> SmallVec4<(usize, Dir, bool)> {
        let exts = self.base.exts[node_id];
        let sequence = self.base.sequences.get(node_id);
        let kmer: K = sequence.term_kmer(dir);
        let mut edges = SmallVec4::new();

        for i in 0..4 {
            if exts.has_ext(dir, i) {
                let link = self.find_link(kmer.extend(i, dir), dir); //.expect("missing link");
                match link {
                    Some(l) => edges.push(l),
                    // This edge doesn't exist within this shard, so ignore it.
                    // NOTE: this should be allowed in a 'complete' DBG
                    None => (),
                }
            }
        }

        edges
    }

    /// Seach for the kmer `kmer`, appearing at the given `side` of a node sequence.
    fn search_kmer(&self, kmer: K, side: Dir) -> Option<usize> {
        match side {
            Dir::Left => match self.left_order.get(&kmer) {
                Some(pos) => Some(*pos as usize),
                _ => None,
            },
            Dir::Right => match self.right_order.get(&kmer) {
                Some(pos) => Some(*pos as usize),
                _ => None,
            },
        }
    }

    /// Find a link in the graph, possibly handling a RC switch.
    pub fn find_link(&self, kmer: K, dir: Dir) -> Option<(usize, Dir, bool)> {
        // Only test self-consistent paths through
        // the edges
        // Avoids problems due to single kmer edges
        // (dir, next_side_incoming, rc)
        // (Dir::Left, Dir::Right, false) => true,
        // (Dir::Left, Dir::Left,  true) => true,
        // (Dir::Right, Dir::Left, false) => true,
        // (Dir::Right, Dir::Right, true) => true,

        let rc = kmer.rc();

        match dir {
            Dir::Left => {
                match self.search_kmer(kmer, Dir::Right) {
                    Some(idx) => return Some((idx, Dir::Right, false)),
                    _ => (),
                }

                if !self.base.stranded {
                    match self.search_kmer(rc, Dir::Left) {
                        Some(idx) => return Some((idx, Dir::Left, true)),
                        _ => (),
                    }
                }
            }

            Dir::Right => {
                match self.search_kmer(kmer, Dir::Left) {
                    Some(idx) => return Some((idx, Dir::Left, false)),
                    _ => (),
                }

                if !self.base.stranded {
                    match self.search_kmer(rc, Dir::Right) {
                        Some(idx) => return Some((idx, Dir::Right, true)),
                        _ => (),
                    }
                }
            }
        }

        return None;
    }

    /// Check whether the graph is fully compressed. Return `None` if it's compressed,
    /// otherwise return `Some(node1, node2)` representing a pair of node that could
    /// be collapsed. Probably only useful for testing.
    pub fn is_compressed<S: CompressionSpec<D>>(&self, spec: &S) -> Option<(usize, usize)> {
        for i in 0..self.len() {
            let n = self.get_node(i);

            for dir in vec![Dir::Left, Dir::Right] {
                let dir_edges = n.edges(dir);
                if dir_edges.len() == 1 {
                    let (next_id, return_dir, _) = dir_edges[0];
                    let next = self.get_node(next_id);

                    let ret_edges = next.edges(return_dir);
                    if ret_edges.len() == 1 {
                        // Test for us being a palindrome: this makes it OK
                        if n.len() == K::k() && n.sequence().first_kmer::<K>().is_palindrome() {
                            continue;
                        }

                        // Test for a neighbor being a palindrome: this makes it OK
                        if next.len() == K::k() && next.sequence().first_kmer::<K>().is_palindrome()
                        {
                            continue;
                        }

                        // Test for this edge representing a smooth circle (biting it's own tail):
                        if n.node_id == next_id {
                            continue;
                        }

                        if spec.join_test(n.data(), next.data()) {
                            // Found a unbranched edge that should have been eliminated
                            return Some((i, next_id));
                        }
                    }
                }
            }
        }

        None
    }

    /// Remove non-existent extensions that may be created due to filtered kmers
    pub fn fix_exts(&mut self, valid_nodes: Option<&BitSet>) {
        for i in 0..self.len() {
            let valid_exts = self.get_valid_exts(i, valid_nodes);
            self.base.exts[i] = valid_exts;
        }
    }

    pub fn get_valid_exts(&self, node_id: usize, valid_nodes: Option<&BitSet>) -> Exts {
        let mut new_exts = Exts::empty();
        let node = self.get_node(node_id);
        let exts = node.exts();
        let l_kmer: K = node.sequence().first_kmer();
        let r_kmer: K = node.sequence().last_kmer();

        let check_node = |id| match valid_nodes {
            Some(ref bs) => bs.contains(id),
            None => true,
        };

        for i in 0..4 {
            if exts.has_ext(Dir::Left, i) {
                match self.find_link(l_kmer.extend_left(i), Dir::Left) {
                    Some((target, _, _)) if check_node(target) => {
                        new_exts = new_exts.set(Dir::Left, i)
                    }
                    _ => (),
                }
            }

            if exts.has_ext(Dir::Right, i) {
                match self.find_link(r_kmer.extend_right(i), Dir::Right) {
                    Some((target, _, _)) if check_node(target) => {
                        new_exts = new_exts.set(Dir::Right, i)
                    }
                    _ => (),
                }
            }
        }

        new_exts
    }

    /// Find the highest-scoring, unambiguous path in the graph. Each node get a score
    /// given by `score`. Any node where `solid_path(node) == True` are valid paths -
    /// paths will be terminated if there are multiple valid paths emanating from a node.
    pub fn max_path<F, F2>(&self, score: F, solid_path: F2) -> Vec<(usize, Dir)>
    where
        F: Fn(&D) -> f32,
        F2: Fn(&D) -> bool,
    {
        if self.len() == 0 {
            return vec![];
        }

        let mut best_node = 0;
        let mut best_score = f32::MIN;
        for i in 0..self.len() {
            let node = self.get_node(i);
            let node_score = score(node.data());

            if node_score > best_score {
                best_node = i;
                best_score = node_score;
            }
        }

        let oscore = |state| match state {
            None => 0.0,
            Some((id, _)) => score(&self.get_node(id).data()),
        };

        let osolid_path = |state| match state {
            None => false,
            Some((id, _)) => solid_path(&self.get_node(id).data()),
        };

        // Now expand in each direction, greedily taking the best path. Stop if we hit a node we've
        // already put into the path
        let mut used_nodes = HashSet::new();
        let mut path = VecDeque::new();

        // Start w/ initial state
        used_nodes.insert(best_node);
        path.push_front((best_node, Dir::Left));

        for init in [(best_node, Dir::Left, false), (best_node, Dir::Right, true)].iter() {
            let &(start_node, dir, do_flip) = init;
            let mut current = (start_node, dir);
            debug!("start: {:?}", current);

            loop {
                let mut next = None;
                let (cur_id, incoming_dir) = current;
                let node = self.get_node(cur_id);
                let edges = node.edges(incoming_dir.flip());
                debug!("{:?}", node);

                let mut solid_paths = 0;
                for (id, dir, _) in edges {
                    let cand = Some((id, dir));
                    if osolid_path(cand) {
                        solid_paths += 1;
                    }

                    if oscore(cand) > oscore(next) {
                        next = cand;
                    }
                }

                if solid_paths > 1 {
                    break;
                }

                match next {
                    Some((next_id, next_incoming)) if !used_nodes.contains(&next_id) => {
                        if do_flip {
                            path.push_front((next_id, next_incoming.flip()));
                        } else {
                            path.push_back((next_id, next_incoming));
                        }

                        used_nodes.insert(next_id);
                        current = (next_id, next_incoming);
                    }
                    _ => break,
                }
            }
        }

        debug!("path:{:?}", path);
        Vec::from_iter(path)
    }

    /// Get the sequence of a path through the graph. The path is given as a sequence of node_id integers
    pub fn sequence_of_path<'a, I: 'a + Iterator<Item = &'a (usize, Dir)>>(
        &self,
        path: I,
    ) -> DnaString {
        let mut seq = DnaString::new();

        for (idx, &(node_id, dir)) in path.enumerate() {
            let start = if idx == 0 { 0 } else { K::k() - 1 };

            let node_seq = match dir {
                Dir::Left => self.get_node(node_id).sequence(),
                Dir::Right => self.get_node(node_id).sequence().rc(),
            };

            for p in start..node_seq.len() {
                seq.push(node_seq.get(p))
            }
        }

        seq
    }

    fn node_to_dot<F: Fn(&D) -> String>(
        &self,
        node: &Node<'_, K, D>,
        node_label: &F,
        f: &mut dyn Write,
    ) {
        let label = node_label(node.data());
        writeln!(
            f,
            "n{} [label=\"id:{} len:{}  {}\",style=filled]",
            node.node_id,
            node.node_id,
            node.sequence().len(),
            label
        )
        .unwrap();

        for (id, incoming_dir, _) in node.l_edges() {
            let color = match incoming_dir {
                Dir::Left => "blue",
                Dir::Right => "red",
            };
            writeln!(f, "n{} -> n{} [color={}]", id, node.node_id, color).unwrap();
        }

        for (id, incoming_dir, _) in node.r_edges() {
            let color = match incoming_dir {
                Dir::Left => "blue",
                Dir::Right => "red",
            };
            writeln!(f, "n{} -> n{} [color={}]", node.node_id, id, color).unwrap();
        }
    }

    /// Write the graph to a dot file.
    pub fn to_dot<P: AsRef<Path>, F: Fn(&D) -> String>(&self, path: P, node_label: &F) {
        let mut f = File::create(path).expect("couldn't open file");

        writeln!(&mut f, "digraph {{").unwrap();
        for i in 0..self.len() {
            self.node_to_dot(&self.get_node(i), node_label, &mut f);
        }
        writeln!(&mut f, "}}").unwrap();
    }

    fn node_to_gfa<F: Fn(&Node<'_, K, D>) -> String>(
        &self,
        node: &Node<'_, K, D>,
        w: &mut dyn Write,
        tag_func: Option<&F>,
    ) -> Result<(), Error> {
        match tag_func {
            Some(f) => {
                let tags = (f)(node);
                writeln!(
                    w,
                    "S\t{}\t{}\t{}",
                    node.node_id,
                    node.sequence().to_dna_string(),
                    tags
                )?;
            }
            _ => writeln!(
                w,
                "S\t{}\t{}",
                node.node_id,
                node.sequence().to_dna_string()
            )?,
        }

        for (target, dir, _) in node.l_edges() {
            if target >= node.node_id as usize {
                let to_dir = match dir {
                    Dir::Left => "+",
                    Dir::Right => "-",
                };
                writeln!(
                    w,
                    "L\t{}\t{}\t{}\t{}\t{}M",
                    node.node_id,
                    "-",
                    target,
                    to_dir,
                    K::k() - 1
                )?;
            }
        }

        for (target, dir, _) in node.r_edges() {
            if target > node.node_id as usize {
                let to_dir = match dir {
                    Dir::Left => "+",
                    Dir::Right => "-",
                };
                writeln!(
                    w,
                    "L\t{}\t{}\t{}\t{}\t{}M",
                    node.node_id,
                    "+",
                    target,
                    to_dir,
                    K::k() - 1
                )?;
            }
        }

        Ok(())
    }

    /// Write the graph to GFA format
    pub fn to_gfa<P: AsRef<Path>>(&self, gfa_out: P) -> Result<(), Error> {
        let wtr = File::create(gfa_out)?;
        self.write_gfa(&mut std::io::BufWriter::new(wtr))
    }

    pub fn write_gfa(&self, wtr: &mut impl Write) -> Result<(), Error> {
        writeln!(wtr, "H\tVN:Z:debruijn-rs")?;

        // Hack to generate a None value with the right type.
        let dummy_func = |_n: &Node<'_, K, D>| "".to_string();
        let mut dummy_opt = Some(&dummy_func);
        let _ = dummy_opt.take();

        for i in 0..self.len() {
            let n = self.get_node(i);
            self.node_to_gfa(&n, wtr, dummy_opt)?;
        }

        Ok(())
    }

    /// Write the graph to GFA format
    pub fn to_gfa_with_tags<P: AsRef<Path>, F: Fn(&Node<'_, K, D>) -> String>(
        &self,
        gfa_out: P,
        tag_func: F,
    ) -> Result<(), Error> {
        let mut wtr = File::create(gfa_out)?;
        writeln!(wtr, "H\tVN:Z:debruijn-rs")?;

        for i in 0..self.len() {
            let n = self.get_node(i);
            self.node_to_gfa(&n, &mut wtr, Some(&tag_func))?;
        }

        Ok(())
    }

    pub fn to_json_rest<W: Write, F: Fn(&D) -> Value>(
        &self,
        fmt_func: F,
        mut writer: &mut W,
        rest: Option<Value>,
    ) {
        writeln!(writer, "{{\n\"nodes\": [").unwrap();
        for i in 0..self.len() {
            let node = self.get_node(i);
            node.to_json(&fmt_func, writer);
            if i == self.len() - 1 {
                write!(writer, "\n").unwrap();
            } else {
                write!(writer, ",\n").unwrap();
            }
        }
        writeln!(writer, "],").unwrap();

        writeln!(writer, "\"links\": [").unwrap();
        for i in 0..self.len() {
            let node = self.get_node(i);
            match node.edges_to_json(writer) {
                true => {
                    if i == self.len() - 1 {
                        write!(writer, "\n").unwrap();
                    } else {
                        write!(writer, ",\n").unwrap();
                    }
                }
                _ => continue,
            }
        }
        writeln!(writer, "]").unwrap();

        match rest {
            Some(Value::Object(v)) => {
                for (k, v) in v.iter() {
                    writeln!(writer, ",").expect("io error");
                    write!(writer, "\"{}\": ", k).expect("io error");
                    serde_json::to_writer(&mut writer, v).expect("io error");
                    writeln!(writer, "").expect("io error");
                }
            }
            _ => {
                writeln!(writer, "").expect("io error");
            }
        }

        writeln!(writer, "}}").expect("io error");
    }

    /// Write the graph to JSON
    pub fn to_json<W: Write, F: Fn(&D) -> Value, RF: Fn(&mut W) -> ()>(
        &self,
        fmt_func: F,
        writer: &mut W,
    ) {
        self.to_json_rest(fmt_func, writer, None);
    }

    /// Print a text representation of the graph.
    pub fn print(&self) {
        println!("DebruijnGraph {{ len: {}, K: {} }} :", self.len(), K::k());
        for node in self.iter_nodes() {
            println!("{:?}", node);
        }
    }

    pub fn print_with_data(&self) {
        println!("DebruijnGraph {{ len: {}, K: {} }} :", self.len(), K::k());
        for node in self.iter_nodes() {
            println!("{:?} ({:?})", node, node.data());
        }
    }

    pub fn max_path_beam<F, F2>(&self, beam: usize, score: F, _solid_path: F2) -> Vec<(usize, Dir)>
    where
        F: Fn(&D) -> f32,
        F2: Fn(&D) -> bool,
    {
        if self.len() == 0 {
            return vec![];
        }

        let mut states = Vec::new();

        for i in 0..self.len() {
            let node = self.get_node(i);

            // Initialize beam search on terminal nodes
            if node.exts().num_exts_l() == 0 || node.exts().num_exts_r() == 0 {
                let dir = if node.exts().num_exts_l() > 0 {
                    Dir::Right
                } else {
                    Dir::Left
                };

                let status = if node.exts().num_exts_l() == 0 && node.exts().num_exts_r() == 0 {
                    Status::End
                } else {
                    Status::Active
                };

                let mut path = SmallVec8::new();
                path.push((i as u32, dir));

                let s = State {
                    path,
                    status,
                    score: score(node.data()),
                };
                states.push(s);
            }
        }

        // No end nodes -- just start on the first node
        if states.len() == 0 {
            // Make a start
            let node = self.get_node(0);
            let mut path = SmallVec8::new();
            path.push((0, Dir::Left));
            states.push(State {
                path: path,
                status: Status::Active,
                score: score(node.data()),
            });
        }

        // Beam search until we can't find any more expansions
        let mut active = true;
        while active {
            let mut new_states = Vec::new();
            active = false;

            for s in states {
                if s.status == Status::Active {
                    active = true;
                    let expanded = self.expand_state(&s, &score);
                    new_states.extend(expanded);
                } else {
                    new_states.push(s)
                }
            }

            // workaround to sort by descending score - will panic if there are NaN scores
            new_states.sort_by(|a, b| (-(a.score)).partial_cmp(&-(b.score)).unwrap());
            new_states.truncate(beam);
            states = new_states;
        }

        for i in 0..min(5, states.len()) {
            trace!("i:{}  -- {:?}", i, states[i]);
        }

        // convert back to using usize for node_id
        states[0]
            .path
            .iter()
            .map(|&(node, dir)| (node as usize, dir))
            .collect()
    }

    fn expand_state<F>(&self, state: &State, score: &F) -> SmallVec4<State>
    where
        F: Fn(&D) -> f32,
    {
        if state.status != Status::Active {
            panic!("only attempt to expand active states")
        }

        let (node_id, dir) = state.path[state.path.len() - 1];
        let node = self.get_node(node_id as usize);
        let mut new_states = SmallVec4::new();

        for (next_node_id, incoming_dir, _) in node.edges(dir.flip()) {
            let next_node = self.get_node(next_node_id);
            let new_score = state.score + score(next_node.data());

            let cycle = state
                .path
                .iter()
                .any(|&(prev_node, _)| prev_node == (next_node_id as u32));

            let status = if cycle {
                Status::Cycle
            } else if next_node.edges(incoming_dir.flip()).len() == 0 {
                Status::End
            } else {
                Status::Active
            };

            let mut new_path = state.path.clone();
            new_path.push((next_node_id as u32, incoming_dir));

            let next_state = State {
                path: new_path,
                score: new_score,
                status: status,
            };

            new_states.push(next_state);
        }

        new_states
    }
}

#[derive(Debug, Eq, PartialEq)]
enum Status {
    Active,
    End,
    Cycle,
}

#[derive(Debug)]
struct State {
    path: SmallVec8<(u32, Dir)>,
    score: f32,
    status: Status,
}

impl State {}

/// Iterator over nodes in a `DeBruijnGraph`
pub struct NodeIter<'a, K: Kmer + 'a, D: Debug + 'a> {
    graph: &'a DebruijnGraph<K, D>,
    node_id: usize,
}

impl<'a, K: Kmer + 'a, D: Debug + 'a> Iterator for NodeIter<'a, K, D> {
    type Item = Node<'a, K, D>;

    fn next(&mut self) -> Option<Node<'a, K, D>> {
        if self.node_id < self.graph.len() {
            let node = self.graph.get_node(self.node_id);
            self.node_id += 1;
            Some(node)
        } else {
            None
        }
    }
}

impl<'a, K: Kmer + 'a, D: Debug + 'a> IntoIterator for &'a DebruijnGraph<K, D> {
    type Item = NodeKmer<'a, K, D>;
    type IntoIter = NodeIntoIter<'a, K, D>;

    fn into_iter(self) -> Self::IntoIter {
        NodeIntoIter {
            graph: self,
            node_id: 0,
        }
    }
}

/// Iterator over nodes in a `DeBruijnGraph`
pub struct NodeIntoIter<'a, K: Kmer + 'a, D: Debug + 'a> {
    graph: &'a DebruijnGraph<K, D>,
    node_id: usize,
}

impl<'a, K: Kmer + 'a, D: Debug + 'a> Iterator for NodeIntoIter<'a, K, D> {
    type Item = NodeKmer<'a, K, D>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.node_id < self.graph.len() {
            let node_id = self.node_id;
            let node = self.graph.get_node(node_id);
            let node_seq = node.sequence();

            self.node_id += 1;
            Some(NodeKmer {
                node_id: node_id,
                node_seq_slice: node_seq,
                phantom_d: PhantomData,
                phantom_k: PhantomData,
            })
        } else {
            None
        }
    }
}

/// A `DebruijnGraph` node with a reference to the sequence of the node.
#[derive(Clone)]
pub struct NodeKmer<'a, K: Kmer + 'a, D: Debug + 'a> {
    pub node_id: usize,
    node_seq_slice: DnaStringSlice<'a>,
    phantom_k: PhantomData<K>,
    phantom_d: PhantomData<D>,
}

/// An iterator over the kmers in a `DeBruijn graph node`
pub struct NodeKmerIter<'a, K: Kmer + 'a, D: Debug + 'a> {
    kmer_id: usize,
    kmer: K,
    num_kmers: usize,
    node_seq_slice: DnaStringSlice<'a>,
    phantom_k: PhantomData<K>,
    phantom_d: PhantomData<D>,
}

impl<'a, K: Kmer + 'a, D: Debug + 'a> IntoIterator for NodeKmer<'a, K, D> {
    type Item = K;
    type IntoIter = NodeKmerIter<'a, K, D>;

    fn into_iter(self) -> Self::IntoIter {
        let num_kmers = self.node_seq_slice.len() - K::k() + 1;

        let kmer = if num_kmers > 0 {
            self.node_seq_slice.get_kmer::<K>(0)
        } else {
            K::empty()
        };

        NodeKmerIter {
            kmer_id: 0,
            kmer,
            num_kmers,
            node_seq_slice: self.node_seq_slice,
            phantom_k: PhantomData,
            phantom_d: PhantomData,
        }
    }
}

impl<'a, K: Kmer + 'a, D: Debug + 'a> Iterator for NodeKmerIter<'a, K, D> {
    type Item = K;

    fn next(&mut self) -> Option<Self::Item> {
        if self.num_kmers == self.kmer_id {
            None
        } else {
            let current_kmer = self.kmer;

            self.kmer_id += 1;
            if self.kmer_id < self.num_kmers {
                let next_base = self.node_seq_slice.get(self.kmer_id + K::k() - 1);
                let new_kmer = self.kmer.extend_right(next_base);
                self.kmer = new_kmer;
            }

            Some(current_kmer)
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        return (self.num_kmers, Some(self.num_kmers));
    }

    /// Provide a 'fast-forward' capability for this iterator
    /// MPHF will use this to reduce the number of kmers that
    /// need to be produced.
    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        if n <= 4 {
            // for small skips forward, shift one base at a time
            for _ in 0..n {
                self.next();
            }
        } else {
            self.kmer_id += n;
            self.kmer = self.node_seq_slice.get_kmer::<K>(self.kmer_id);
        }

        self.next()
    }
}

/// Marker signifying that NodeKmerIter has a known size.
impl<'a, K: Kmer + 'a, D: Debug + 'a> ExactSizeIterator for NodeKmerIter<'a, K, D> {}

/// Unbranched sequence in the DeBruijn graph
pub struct Node<'a, K: Kmer + 'a, D: 'a> {
    pub node_id: usize,
    pub graph: &'a DebruijnGraph<K, D>,
}

impl<'a, K: Kmer, D: Debug> Node<'a, K, D> {
    /// Length of the sequence of this node
    pub fn len(&self) -> usize {
        self.graph.base.sequences.get(self.node_id).len()
    }

    /// Sequence of the node
    pub fn sequence(&self) -> DnaStringSlice<'a> {
        self.graph.base.sequences.get(self.node_id)
    }

    /// Reference to auxiliarly data associated with the node
    pub fn data(&self) -> &'a D {
        &self.graph.base.data[self.node_id]
    }

    /// Extension bases from this node
    pub fn exts(&self) -> Exts {
        self.graph.base.exts[self.node_id]
    }

    /// Edges leaving the left side of the node in the format
    //// (target_node id, incoming side of target node, whether target node has is flipped)
    pub fn l_edges(&self) -> SmallVec4<(usize, Dir, bool)> {
        self.graph.find_edges(self.node_id, Dir::Left)
    }

    /// Edges leaving the right side of the node in the format
    //// (target_node id, incoming side of target node, whether target node has is flipped)
    pub fn r_edges(&self) -> SmallVec4<(usize, Dir, bool)> {
        self.graph.find_edges(self.node_id, Dir::Right)
    }

    /// Edges leaving the 'dir' side of the node in the format
    //// (target_node id, incoming side of target node, whether target node has is flipped)
    pub fn edges(&self, dir: Dir) -> SmallVec4<(usize, Dir, bool)> {
        self.graph.find_edges(self.node_id, dir)
    }

    fn to_json<F: Fn(&D) -> Value>(&self, func: &F, f: &mut dyn Write) {
        write!(
            f,
            "{{\"id\":\"{}\",\"L\":{},\"D\":{},\"Se\":\"{:?}\"}}",
            self.node_id,
            self.sequence().len(),
            (func)(self.data()),
            self.sequence(),
        )
        .unwrap();
    }

    fn edges_to_json(&self, f: &mut dyn Write) -> bool {
        let mut wrote = false;
        let edges = self.r_edges();
        for (idx, &(id, incoming_dir, _)) in edges.iter().enumerate() {
            write!(
                f,
                "{{\"source\":\"{}\",\"target\":\"{}\",\"D\":\"{}\"}}",
                self.node_id,
                id,
                match incoming_dir {
                    Dir::Left => "L",
                    Dir::Right => "R",
                }
            )
            .unwrap();

            if idx < edges.len() - 1 {
                write!(f, ",").unwrap();
            }

            wrote = true;
        }
        wrote
    }
}

/*
impl<'a, K: Kmer, D> fmt::Debug for Node<'a, K, D> where D:Debug {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Node: {}, Seq: {}, Exts:{:?}, Data: {:?}", self.node_id, self.sequence().to_string(), self.exts(), self.data())
    }
}
*/

impl<'a, K: Kmer, D> fmt::Debug for Node<'a, K, D>
where
    D: Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Node {{ id:{}, Exts: {:?}, L:{:?} R:{:?}, Seq: {:?}, Data: {:?} }}",
            self.node_id,
            self.exts(),
            self.l_edges(),
            self.r_edges(),
            self.sequence().len(),
            self.data()
        )
    }
}
