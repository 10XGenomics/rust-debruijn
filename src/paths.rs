//! Compute path-compressed De Bruijn graphs from kmers or intermediate sized De Bruijn graph fragments

use std::marker::PhantomData;
use fx::{FxHashMap, FxLMap, FxHasher};
use std::collections::VecDeque;
use bit_set::BitSet;
use smallvec::SmallVec;
use std::iter::FromIterator;
use std::collections::HashSet;
use std::io::Write;
use std::fs::File;
use std::path::Path;
use std::fmt::{self, Debug};

use std::hash::BuildHasherDefault;
use std::ops::Index;
use std::f32;

use serde_json;
use serde_json::Value;

type SmallVec4<T> = SmallVec<[T; 4]>;

use Mer;
use Kmer;
use Vmer;
use Dir;
use Exts;
use dna_string::{DnaString, DnaStringSlice};

#[derive(Serialize, Deserialize)]
pub struct PackedDnaStringSet {
    pub sequence: DnaString,
    pub start: Vec<usize>,
    pub length: Vec<u32>,
}

impl<'a> PackedDnaStringSet {
    fn new() -> Self {
        PackedDnaStringSet {
            sequence: DnaString::new(),
            start: Vec::new(),
            length: Vec::new(),
        }
    }

    pub fn get(&'a self, i: usize) -> DnaStringSlice<'a> {
        DnaStringSlice {
            dna_string: &self.sequence,
            start: self.start[i],
            length: self.length[i] as usize,
            is_rc: false,
        }
    }

    pub fn len(&self) -> usize {
        self.start.len()
    }

    fn add<'b, S: IntoIterator<Item = &'b u8>>(&mut self, sequence: S) {
        let start = self.sequence.len();
        self.start.push(start);

        let mut length = 0;
        for b in sequence {
            self.sequence.push(*b);
            length += 1;
        }
        self.length.push(length as u32);
    }
}

#[derive(Serialize, Deserialize)]
pub struct BaseGraph<K, D> {
    pub sequences: PackedDnaStringSet,
    pub exts: Vec<Exts>,
    pub data: Vec<D>,
    phantom: PhantomData<K>,
}

impl<K, D> BaseGraph<K, D> {
    pub fn new() -> Self {
        BaseGraph {
            sequences: PackedDnaStringSet::new(),
            exts: Vec::new(),
            data: Vec::new(),
            phantom: PhantomData,
        }
    }

    pub fn len(&self) -> usize {
        self.sequences.len()
    }
}

impl<K: Kmer, D> BaseGraph<K, D> {
    pub fn add<'a, S: IntoIterator<Item = &'a u8>>(&mut self, sequence: S, exts: Exts, data: D) {
        self.sequences.add(sequence);
        self.exts.push(exts);
        self.data.push(data);
    }

    pub fn finish(self) -> DebruijnGraph<K, D> {

        let mut left_sort: Vec<u32> = Vec::with_capacity(self.len());
        let mut right_sort: Vec<u32> = Vec::with_capacity(self.len());
        for i in 0..self.len() {
            left_sort.push(i as u32);
            right_sort.push(i as u32);
        }

        left_sort.sort_by_key(|idx| -> K { self.sequences.get(*idx as usize).first_kmer() });
        right_sort.sort_by_key(|idx| -> K { self.sequences.get(*idx as usize).last_kmer() });

        DebruijnGraph {
            base: self,
            left_order: left_sort,
            right_order: right_sort,
        }
    }
}

#[derive(Serialize, Deserialize)]
pub struct DebruijnGraph<K, D> {
    base: BaseGraph<K, D>,
    left_order: Vec<u32>,
    right_order: Vec<u32>,
}

impl<K: Kmer, D:Debug> DebruijnGraph<K, D> {
    pub fn len(&self) -> usize {
        self.base.len()
    }

    pub fn get_node<'a>(&'a self, node_id: usize) -> Node<'a, K, D> {
        Node {
            node_id: node_id,
            graph: self,
            //l_edges: self.find_edges(node_id, Dir::Left),
            //r_edges: self.find_edges(node_id, Dir::Right),
        }
    }

    pub fn find_edges(&self, node_id: usize, dir: Dir) -> SmallVec4<(usize, Dir, bool)> {

        let exts = self.base.exts[node_id];
        let sequence = self.base.sequences.get(node_id);
        let kmer: K = sequence.term_kmer(dir);
        let mut edges = SmallVec4::new();

        for i in 0..4 {
            if exts.has_ext(dir, i) {
                let link = self.find_link(kmer.extend(i, dir), dir)
                    .expect("missing link");
                edges.push(link);
            }
        }

        edges
    }

    pub fn search_kmer(&self, kmer: K, side: Dir) -> Option<usize> {
        match side {
            Dir::Left => {
                let pos = self.left_order
                    .binary_search_by_key(&kmer, |idx| {
                        self.base.sequences.get(*idx as usize).first_kmer()
                    });

                match pos {
                    Ok(idx) => Some(self.left_order[idx] as usize),
                    _ => None,
                }
            }
            Dir::Right => {
                let pos =
                    self.right_order
                        .binary_search_by_key(&kmer, |idx| {
                            self.base.sequences.get(*idx as usize).last_kmer()
                        });
                match pos {
                    Ok(idx) => Some(self.right_order[idx] as usize),
                    _ => None,
                }
            }
        }
    }

    pub fn find_link(&self, kmer: K, dir: Dir) -> Option<(usize, Dir, bool)> {
        let rc = kmer.rc();

        // FIXME -- should this be constrained for allow_rc???

        // Only test self-consistent paths through
        // the edges
        // Avoids problems due to single kmer edges
        // (dir, next_side_incoming, rc)
        // (Dir::Left, Dir::Right, false) => true,
        // (Dir::Left, Dir::Left,  true) => true,
        // (Dir::Right, Dir::Left, false) => true,
        // (Dir::Right, Dir::Right, true) => true,
        //

        match dir {
            Dir::Left => {
                match self.search_kmer(kmer, Dir::Right) {
                    Some(idx) => return Some((idx, Dir::Right, false)),
                    _ => (),
                }

                match self.search_kmer(rc, Dir::Left) {
                    Some(idx) => return Some((idx, Dir::Left, true)),
                    _ => (),
                }
            }

            Dir::Right => {
                match self.search_kmer(kmer, Dir::Left) {
                    Some(idx) => return Some((idx, Dir::Left, false)),
                    _ => (),
                }

                match self.search_kmer(rc, Dir::Right) {
                    Some(idx) => return Some((idx, Dir::Right, true)),
                    _ => (),
                }
            }
        }

        return None;
    }

    /// Remove non-existent extensions that may be created due to filtered kmers
    pub fn fix_exts(&mut self, valid_nodes: Option<&BitSet>)
    {
        for i in 0 .. self.len() {
            let valid_exts = self.get_valid_exts(i, valid_nodes);
            self.base.exts[i] = valid_exts;
        }
    }

    pub fn get_valid_exts(&self, node_id: usize, valid_nodes: Option<&BitSet>) -> Exts
    {
        let mut new_exts = Exts::empty();
        let node = self.get_node(node_id);
        let exts = node.exts();
        let l_kmer: K = node.sequence().first_kmer();
        let r_kmer: K = node.sequence().last_kmer();

        let check_node = |id| {
            match valid_nodes {
                Some(ref bs) => bs.contains(id),
                None => true,
            }
        };

        for i in 0..4
        {
            if exts.has_ext(Dir::Left, i) {
                match self.find_link(l_kmer.extend_left(i), Dir::Left) {
                    Some((target, _, _)) if check_node(target) => 
                        new_exts = new_exts.set(Dir::Left, i),
                    _ => (),
                }
            }

            if exts.has_ext(Dir::Right, i) {
                match self.find_link(r_kmer.extend_right(i), Dir::Right) {
                    Some((target, _, _)) if check_node(target) =>
                        new_exts = new_exts.set(Dir::Right, i),
                    _ => (),
                }
            }
        }

        new_exts
    }


    pub fn max_path<F, F2>(&self, score: F, solid_path: F2) -> Vec<(usize, Dir)>
        where F: Fn(&D) -> f32, F2: Fn(&D) -> bool {

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

        let oscore = |state| {
            match state  {
                None => 0.0,
                Some((id, _)) => score(&self.get_node(id).data()),
            }
        };

        let osolid_path = |state| {
            match state  {
                None => false,
                Some((id, _)) => solid_path(&self.get_node(id).data()),
            }
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
            println!("start: {:?}", current);

            loop {
                let mut next = None;
                let (cur_id, incoming_dir) = current;
                let node = self.get_node(cur_id);
                let edges = node.edges(incoming_dir.flip());
                println!("{:?}", node);

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
                    },
                    _ => break,
                }
            }
        }

        println!("path:{:?}", path);
        Vec::from_iter(path)
    }


    pub fn sequence_of_path(&self, path: &Vec<(usize, Dir)>) -> DnaString {
        let mut seq = DnaString::new();

        for (idx, &(node_id, dir)) in path.iter().enumerate() {
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


    fn node_to_dot<F: Fn(&D) -> String>(&self, node: &Node<K,D>, node_label: &F, f: &mut Write) {

        let label = node_label(node.data());
        writeln!(f, "n{} [label=\"id:{} len:{}  {}\",style=filled]", node.node_id, node.node_id, node.sequence().len(), label).unwrap();


        for (id, incoming_dir, _) in node.l_edges() {
            let color = match incoming_dir { Dir::Left => "blue", Dir::Right => "red"};
            writeln!(f, "n{} -> n{} [color={}]", id, node.node_id, color).unwrap();
        }

        for (id, incoming_dir, _) in node.r_edges() {
            let color = match incoming_dir { Dir::Left => "blue", Dir::Right => "red"};
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


    fn node_to_gfa(&self, node: &Node<K,D>, w: &mut Write) {
        writeln!(w, "S\t{}\t{}", node.node_id, node.sequence().to_dna_string()).unwrap();

        for (target, dir, _) in node.l_edges() {
            if target > node.node_id as usize {
                let to_dir = match dir { Dir::Left => "+", Dir::Right => "-" };
                writeln!(w, "L\t{}\t{}\t{}\t{}\t{}M", node.node_id, "-", target, to_dir, K::k()-1).unwrap();
            }
        }

        for (target, dir, _) in node.r_edges() {
            if target > node.node_id as usize {
                let to_dir = match dir { Dir::Left => "+", Dir::Right => "-" };
                writeln!(w, "L\t{}\t{}\t{}\t{}\t{}M", node.node_id, "+", target, to_dir, K::k()-1).unwrap();
            }
        }
    }

    /// Write the graph to GFA format
    pub fn to_gfa<P: AsRef<Path>>(&self, gfa_out: P) 
    {
        let mut wtr = File::create(gfa_out).unwrap();
        writeln!(wtr, "H\tVN:Z:debruijn-rs").unwrap();

        for i in 0 .. self.len() {
            let n = self.get_node(i);
            self.node_to_gfa(&n, &mut wtr);
        }
    }


    pub fn to_json_rest<W: Write, F: Fn(&D) -> Value>(&self, fmt_func: F, mut writer: &mut W, rest: Option<Value>) {

        writeln!(writer, "{{\n\"nodes\": [").unwrap();
        for i in 0 .. self.len() {
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
        for i in 0 .. self.len() {
            let node = self.get_node(i);
            match node.edges_to_json(writer) {
                true => {
                    if i == self.len() - 1 {
                        write!(writer, "\n").unwrap();
                    } else {
                        write!(writer, ",\n").unwrap();
                    }
                },
                _ => continue,
            }
        }
        writeln!(writer, "]").unwrap();

        match rest {
            Some(Value::Object(v)) => {
                for (k,v) in v.iter() {
                    writeln!(writer, ",");
                    write!(writer, "\"{}\": ", k);
                    serde_json::to_writer(&mut writer, v);
                    writeln!(writer, "");
                }
            },
            _ => { writeln!(writer, ""); }
        }

        writeln!(writer, "}}");
    }

    
    pub fn to_json<W: Write, F: Fn(&D) -> Value, RF: Fn(&mut W) -> ()>(&self, fmt_func: F, writer: &mut W) {
        self.to_json_rest(fmt_func, writer, None);
    }
}


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

    fn to_json<F: Fn(&D) -> Value>(&self, func: &F, f: &mut Write) {
        write!(f,
               "{{\"id\":\"{}\",\"L\":{},\"D\":{},\"Se\":\"{:?}\"}}",
               self.node_id,
               self.sequence().len(),
               (func)(self.data()),
               self.sequence(),
        ).unwrap();
    }

    fn edges_to_json(&self, f: &mut Write) -> bool {
        let mut wrote = false;
        let edges = self.r_edges();
        for (idx, &(id, incoming_dir, _)) in edges.iter().enumerate() {
            write!(f,
                     "{{\"source\":\"{}\",\"target\":\"{}\",\"D\":\"{}\"}}",
                     self.node_id,
                     id,
                     match incoming_dir {
                         Dir::Left => "L",
                         Dir::Right => "R",
                     }
            ).unwrap();

            if idx < edges.len() - 1 {
                write!(f, ",").unwrap();
            }

            wrote = true;
        }
        wrote
    }
}


impl<'a, K: Kmer, D> fmt::Debug for Node<'a, K, D> where D:Debug {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Node: {}, Seq: {}, Exts:{:?}, Data: {:?}", self.node_id, self.sequence().to_string(), self.exts(), self.data())
    }
}


#[derive(Copy, Clone)]
enum ExtMode<K: Kmer> {
    Unique(K, Dir, Exts),
    Terminal(Exts),
}


#[derive(Copy, Clone)]
enum ExtModeNode {
    Unique(usize, Dir, Exts),
    Terminal(Exts),
}


type Seq = Vec<u8>;

pub trait BuildGraph<K: Kmer, V: Vmer<K>, D> {
    // Basic progression
    fn kmers_from_sequences(seqs: Vec<(Seq, D)>) -> Vec<(K, Exts, D)>;
    fn spaths_from_kmers(kmers: Vec<(K, Exts, D)>) -> Vec<(V, Exts, D)>;
    fn graph_from_spaths(spaths: Vec<(V, Exts, D)>) -> DebruijnGraph<K, D>;

    // Derived methods
    fn spaths_from_kmers_dict<I: Index<K, Output = (Exts, D)>>(kmers: I) -> Vec<(V, Exts, D)>;
    fn build_graph<I: Index<K, Output = Exts>>(kmers: I) -> DebruijnGraph<K, D>;
}


/// Customize the path-compression process. Implementing this trait lets the user
/// control how the per-kmer data (of type `D`) is summarized into a per-path
/// summary data (of type `DS`). It also let's the user introduce new breaks into
/// paths by inspecting in the per-kmer data of a proposed with `join_test_kmer`
/// function.
pub trait CompressionSpec<D> {
    fn reduce(&self, path_object: D, kmer_object: &D) -> D;
    fn join_test(&self, d1: &D, d2: &D) -> bool;
}

/// Simple implementation of `CompressionSpec` that lets you provide that data reduction function as a closure
pub struct SimpleCompress<D, F> {
    func: F,
    d: PhantomData<D>,
}

impl<D, F> SimpleCompress<D, F> {
    pub fn new(func: F) -> SimpleCompress<D, F> {
        SimpleCompress {
            func: func,
            d: PhantomData,
        }
    }
}

impl<D, F> CompressionSpec<D> for SimpleCompress<D,F> 
    where for<'r> F: Fn(D, &'r D) -> D {
    fn reduce(&self, d: D, other: &D) -> D {
        (self.func)(d, other)
    }

    fn join_test(&self, _: &D, _: &D) -> bool {
        true
    }
}


pub struct PathCompression<K: Kmer, V: Vmer<K>, D, S: CompressionSpec<D>> {
    allow_rc: bool,
    k: PhantomData<K>,
    v: PhantomData<V>,
    d: PhantomData<D>,
    spec: S,
}

/// Compression of paths in Debruijn graph
impl<K: Kmer, V: Vmer<K>, D: Clone + Debug, S: CompressionSpec<D>> PathCompression<K, V, D, S> {

    pub fn new(allow_rc: bool, spec: S) -> Self {
        PathCompression {
            allow_rc: allow_rc,
            spec: spec,
            k: PhantomData,
            v: PhantomData,
            d: PhantomData
        }
    }

    fn get_kmer_data<'a>(&self, kmer: &K, kmer_exts: &'a Vec<(K, (Exts, D))>) -> &'a (Exts, D) {
        let pos = kmer_exts.binary_search_by_key(kmer, |x| x.0).expect("couldn't find kmer");
        &kmer_exts[pos].1
    }

    /// Attempt to extend kmer v in direction dir. Return:
    ///  - Unique(nextKmer, nextDir) if a single unique extension
    ///    is possible.  nextDir indicates the direction to extend nextMker
    ///    to preserve the direction of the extension.
    /// - Term(ext) no unique extension possible, indicating the extensions at this end of the line
    fn try_extend_kmer(&self,
                       available_kmers: &FxLMap<K, ()>,
                       kmer_exts: &Vec<(K, (Exts, D))>,
                       kmer: K,
                       dir: Dir)
                       -> ExtMode<K> {

        // metadata of start kmer

        let &(exts, ref kmer_data) = self.get_kmer_data(&kmer, kmer_exts);

        if exts.num_ext_dir(dir) != 1 || kmer == kmer.rc() {
            ExtMode::Terminal(exts.single_dir(dir))
        } else {
            // Get the next kmer
            let ext_base = exts.get_unique_extension(dir).expect("should be unique");

            let mut next_kmer = kmer.extend(ext_base, dir);
            let mut do_flip = false;

            if self.allow_rc {
                let flip_rc = next_kmer.min_rc_flip();
                do_flip = flip_rc.1;
                next_kmer = flip_rc.0;
            }

            //let is_palindrome = next_kmer == next_kmer.rc();
            let is_palindrome = false;

            let next_dir = dir.cond_flip(do_flip);

            // We can include this kmer in the line if:
            // a) it exists in the partition, and is still unused
            // b) the kmer we go to has a unique extension back in our direction
            if !available_kmers.contains_key(&next_kmer) {
                // This kmer isn't in this partition, or we've already used it
                return ExtMode::Terminal(exts.single_dir(dir));
            }

            // Direction we're approaching the new kmer from
            let new_incoming_dir = dir.flip().cond_flip(do_flip);
            let next_kmer_r = self.get_kmer_data(&next_kmer, kmer_exts);
            let &(next_kmer_exts, ref next_kmer_data) = next_kmer_r;
            let incoming_count = next_kmer_exts.num_ext_dir(new_incoming_dir);
            let outgoing_exts = next_kmer_exts.single_dir(new_incoming_dir.flip());

            // Test if the spec let's us combine these into the same path
            let can_join = self.spec.join_test(kmer_data, next_kmer_data);

            if incoming_count == 0 && !is_palindrome {
                panic!("unreachable");
            } else if can_join && incoming_count == 1 && !is_palindrome {
                // We have a unique path to next_kmer -- include it
                ExtMode::Unique(next_kmer, next_dir, outgoing_exts)
            } else {
                // there's more than one path
                // into the target kmer - don't include it
                ExtMode::Terminal(exts.single_dir(dir))
            }
        }
    }


    /// Build the maximal line starting at kmer in direction dir, at most max_dist long.
    /// Also return the extensions at the end of this line.
    /// Sub-lines break if their extensions are not available in this shard
    #[inline(never)]
    fn extend_kmer(&self,
                   kmer_exts: &Vec<(K, (Exts, D))>,
                   available_kmers: &mut FxLMap<K, ()>,
                   kmer: K,
                   start_dir: Dir,
                   max_dist: usize,
                   path: &mut Vec<(K, Dir)>)
                   -> Exts {

        let mut current_dir = start_dir;
        let mut current_kmer = kmer;
        path.clear();

        let mut final_exts: Exts; // must get set below

        if max_dist == 0 {
            let first_exts = self.get_kmer_data(&current_kmer, kmer_exts).0;
            return first_exts.single_dir(start_dir);
        }

        let _ = available_kmers.remove(&kmer);

        loop {
            let ext_result =
                self.try_extend_kmer(available_kmers, kmer_exts, current_kmer, current_dir);

            match ext_result {
                ExtMode::Unique(next_kmer, next_dir, next_ext) => {
                    path.push((next_kmer, next_dir));
                    available_kmers.remove(&next_kmer);
                    current_kmer = next_kmer;
                    current_dir = next_dir;
                    final_exts = next_ext
                }
                ExtMode::Terminal(ext) => {
                    final_exts = ext;
                    break;
                }
            }

            if path.len() >= max_dist {
                break;
            }
        }

        final_exts
    }


    /// Build the edge surrounding a kmer
    #[inline(never)]
    fn build_node(&self,
                  kmer_exts: &Vec<(K, (Exts, D))>,
                  available_kmers: &mut FxLMap<K, ()>,
                  seed: K,
                  path: &mut Vec<(K, Dir)>,
                  edge_seq: &mut VecDeque<u8>)
                  -> (Exts, D) {


        edge_seq.clear();
        for i in 0..K::k() {
            edge_seq.push_back(seed.get(i));
        }

        let mut node_data = self.get_kmer_data(&seed, kmer_exts).1.clone();

        let l_ext = self.extend_kmer(kmer_exts,
                                    available_kmers,
                                    seed,
                                    Dir::Left,
                                    V::max_len() - K::k(),
                                    path);

        // Add on the left path
        for &(next_kmer, dir) in path.iter() {
            let kmer = match dir {
                Dir::Left => next_kmer,
                Dir::Right => next_kmer.rc(),
            };

            edge_seq.push_front(kmer.get(0));

            // Reduce the data object
            let &(_, ref kmer_data) = self.get_kmer_data(&next_kmer, kmer_exts);
            node_data = self.spec.reduce(node_data, kmer_data)
        }

        let left_extend = match path.last() {
            None => l_ext,
            Some(&(_, Dir::Left)) => l_ext,
            Some(&(_, Dir::Right)) => l_ext.complement(),
        };


        let r_ext = self.extend_kmer(kmer_exts,
                             available_kmers,
                             seed,
                             Dir::Right,
                             V::max_len() - edge_seq.len(),
                             path);


        // Add on the right path
        for &(next_kmer, dir) in path.iter() {
            let kmer = match dir {
                Dir::Left => next_kmer.rc(),
                Dir::Right => next_kmer,
            };

            edge_seq.push_back(kmer.get(K::k() - 1));

            let &(_, ref kmer_data) = self.get_kmer_data(&next_kmer, kmer_exts);
            node_data = self.spec.reduce(node_data, kmer_data)
        }

        let right_extend = match path.last() {
            None => r_ext,
            Some(&(_, Dir::Left)) => r_ext.complement(),
            Some(&(_, Dir::Right)) => r_ext,
        };

        let mut edge = V::new(edge_seq.len());
        for (pos, base) in edge_seq.iter().enumerate() {
            edge.set_mut(pos, *base);
        }
        (Exts::from_single_dirs(left_extend, right_extend), node_data)
    }

    /// Build all sedges until all kmers are exhausted
    #[inline(never)]
    pub fn build_nodes(&self, kmer_exts: &Vec<(K, (Exts, D))>) -> BaseGraph<K,D> {

        let h = BuildHasherDefault::<FxHasher>::default();
        let mut available_kmers: FxLMap<K, ()> = FxLMap::with_capacity_and_hasher(kmer_exts.len(), h);
        available_kmers.extend(kmer_exts.iter().map(|&(k, _)| (k.clone(), ())));

        // Path-compressed De Bruijn graph will be created here
        let mut graph = BaseGraph::new();

        // Paths will be get assembled here
        let mut path_buf = Vec::new();

        // Node sequences will get assembled here
        let mut edge_seq_buf = VecDeque::new();

        loop {
            match available_kmers.front() {
                Some((&start_kmer, _)) => {
                    let (node_exts, node_data) = self.build_node(kmer_exts, 
                    &mut available_kmers, 
                    start_kmer,
                    &mut path_buf,
                    &mut edge_seq_buf);

                    graph.add(&edge_seq_buf, node_exts, node_data);
                }
                None => break,
            }
        }

        graph
    }


    /// Try to extend edge 'edge' in direction 'dir', returning:
    /// - Unique(usize, dir, exts) if there a unique extension into a another sedge
    /// - Terminal(exts) if there isn't -- exts denotes the neighboring kmers
    #[inline(never)]
    fn try_extend_node(graph: &DebruijnGraph<K,D>,
                    available_nodes: &BitSet,
                    node_id: usize,
                    dir: Dir)
                    -> ExtModeNode {

        let node = graph.get_node(node_id);
        let exts = node.exts();
        let bases = node.sequence();

        if exts.num_ext_dir(dir) != 1 || (bases.len() == K::k()) { // && bases.is_palindrome()) {
            ExtModeNode::Terminal(exts.single_dir(dir))
        } else {
            // Get the next kmer
            let ext_base = exts.get_unique_extension(dir).expect("should be unique");
            let end_kmer: K = bases.term_kmer(dir);

            let next_kmer = end_kmer.extend(ext_base, dir);
            let (next_node_id, next_side_incoming, rc) = match graph.find_link(next_kmer, dir) {
                Some(e) => e,
                None => {
                    println!("dir: {:?}, Lmer: {:?}, exts: {:?}", dir, bases, exts);
                    println!("end kmer: {:?}", end_kmer);
                    println!("No kmer: {:?}", next_kmer);
                    println!("rc: {:?}", next_kmer.min_rc());
                    panic!(format!("No kmer: {:?}", next_kmer))
                }
            };

            let next_node = graph.get_node(next_node_id);
            let next_exts = node.exts();
            let next_bases = node.sequence();

            let consistent = (next_bases.len() == K::k()) ||
                            match (dir, next_side_incoming, rc) {
                (Dir::Left, Dir::Right, false) => true,
                (Dir::Left, Dir::Left, true) => true,
                (Dir::Right, Dir::Left, false) => true,
                (Dir::Right, Dir::Right, true) => true,
                _ => {
                    println!("dir: {:?}, Lmer: {:?}, exts: {:?}", dir, bases, exts);
                    println!("end kmer: {:?}", end_kmer);
                    println!("next kmer: {:?}", next_kmer);
                    println!("rc: {:?}", next_kmer.min_rc());
                    println!("next bases: {:?}, next_side_incoming: {:?}, rc: {:?}",
                            next_bases,
                            next_side_incoming,
                            rc);
                    false
                }
            };
            assert!(consistent);

            // We can include this kmer in the line if:
            // a) it exists in the partition, and is still unused
            // b) the kmer we go to has a unique extension back in our direction
            // c) the new edge is not of length K and a palindrome

            if !available_nodes.contains(next_node_id) || (next_bases.len() == K::k()) { // && next_bases.is_palindrome()) {
                // This kmer isn't in this partition, or we've already used it
                return ExtModeNode::Terminal(exts.single_dir(dir));
            }

            // orientation of next edge
            let next_side_outgoing = next_side_incoming.flip();

            let incoming_count = next_exts.num_ext_dir(next_side_incoming);
            let outgoing_exts = next_exts.single_dir(next_side_outgoing);

            if incoming_count == 0 {
                println!("dir: {:?}, Lmer: {:?}, exts: {:?}", dir, bases, exts);
                println!("end kmer: {:?}", end_kmer);
                panic!("unreachable");
            } else if incoming_count == 1 {
                // We have a unique path to next_kmer -- include it
                ExtModeNode::Unique(next_node_id, next_side_outgoing, outgoing_exts)
            } else {
                // there's more than one path
                // into the target kmer - don't include it
                ExtModeNode::Terminal(exts.single_dir(dir))
            }
        }
    }
}
