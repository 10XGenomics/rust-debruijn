//! Compute path-compressed De Bruijn graphs from kmers or intermediate sized De Bruijn graph fragments

use std::marker::PhantomData;
use std::collections::VecDeque;
use bit_set::BitSet;
use smallvec::SmallVec;
use std::iter::FromIterator;
use std::collections::HashSet;
use std::io::Write;
use std::fs::File;
use std::path::Path;
use std::fmt::{self, Debug};
use std::borrow::Borrow;

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

    fn add<'b, R: Borrow<u8>, S: IntoIterator<Item = R>>(&mut self, sequence: S) {
        let start = self.sequence.len();
        self.start.push(start);

        let mut length = 0;
        for b in sequence {
            self.sequence.push(b.borrow().clone());
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

    pub fn combine<I: Iterator<Item=BaseGraph<K,D>>>(graphs: I) -> Self {
        let mut sequences = PackedDnaStringSet::new();
        let mut exts = Vec::new();
        let mut data = Vec::new();
        let mut stranded = Vec::new();

        for g in graphs {
            for s in 0 .. g.sequences.len() {
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
    pub fn add<'b, R: Borrow<u8>, S: IntoIterator<Item = R>>(&mut self, sequence: S, exts: Exts, data:D) {
    //pub fn add<'a, S: IntoIterator<Item = u8>>(&mut self, sequence: S, exts: Exts, data: D) {
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
    pub base: BaseGraph<K, D>,
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
        }
    }

    pub fn iter_nodes<'a>(&'a self) -> NodeIter<'a, K, D> {
        NodeIter { graph: self, node_id: 0 }
    }

    pub fn find_edges(&self, node_id: usize, dir: Dir) -> SmallVec4<(usize, Dir, bool)> {

        let exts = self.base.exts[node_id];
        let sequence = self.base.sequences.get(node_id);
        let kmer: K = sequence.term_kmer(dir);
        let mut edges = SmallVec4::new();

        for i in 0..4 {
            if exts.has_ext(dir, i) {
                let link = self.find_link(kmer.extend(i, dir), dir);//.expect("missing link");
                match link {
                    Some(l) => edges.push(l),
                    None => {
                        println!("seq: {:?}, kmer: {:?}, exts: {:?}, next:{:?}, p: {}",
                            sequence, kmer, exts, kmer.extend(i, dir), kmer.is_palindrome() );
                        panic!("no link");
                    },
                }
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

    pub fn is_compressed(&self) -> Option<(usize, usize)>  {
        for i in 0 .. self.len() {
            let n = self.get_node(i);

            for dir in vec![Dir::Left, Dir::Right] {
            
                let dir_edges = n.edges(dir);
                if dir_edges.len() == 1 {
                    let (next_id, return_dir, _) = dir_edges[0];
                    let next = self.get_node(next_id);

                    let ret_edges = next.edges(return_dir);
                    if ret_edges.len() == 1 {

                        if n.len() == K::k() && Vmer::<K>::first_kmer(&n.sequence()).is_palindrome() {
                            return None
                        }

                        if next.len() == K::k() && Vmer::<K>::first_kmer(&next.sequence()).is_palindrome() {
                            return None
                        }

                        // Found a unbranched edge that should have been eliminated
                        return Some((i, next_id));
                    }
                }
            }
        }

        None
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
                    },
                    _ => break,
                }
            }
        }

        debug!("path:{:?}", path);
        Vec::from_iter(path)
    }


    pub fn sequence_of_path<'a, I: 'a + Iterator<Item=&'a (usize, Dir)>>(&self, path: I) -> DnaString {
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


    fn node_to_gfa<F: (Fn(&Node<K,D>) -> String)>(&self, node: &Node<K,D>, w: &mut Write, tag_func: Option<&F>) {
        
        match tag_func {
            Some(f) => {
                let tags = (f)(node);
                writeln!(w, "S\t{}\t{}\t{}", node.node_id, node.sequence().to_dna_string(), tags).unwrap();
            }
            _ => writeln!(w, "S\t{}\t{}", node.node_id, node.sequence().to_dna_string()).unwrap(),
        }

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

        // Hack to generate a None value with the right type.
        let dummy_func = |n: &Node<K,D>| { "".to_string() };
        let mut dummy_opt = Some(&dummy_func);
        let _ = dummy_opt.take();

        for i in 0 .. self.len() {
            let n = self.get_node(i);
            self.node_to_gfa(&n, &mut wtr, dummy_opt);
        }
    }

    /// Write the graph to GFA format
    pub fn to_gfa_with_tags<P: AsRef<Path>, F: (Fn(&Node<K,D>) -> String)>(&self, gfa_out: P, tag_func: F) 
    {
        let mut wtr = File::create(gfa_out).unwrap();
        writeln!(wtr, "H\tVN:Z:debruijn-rs").unwrap();

        for i in 0 .. self.len() {
            let n = self.get_node(i);
            self.node_to_gfa(&n, &mut wtr, Some(&tag_func));
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
                    writeln!(writer, ",").expect("io error");
                    write!(writer, "\"{}\": ", k).expect("io error");
                    serde_json::to_writer(&mut writer, v).expect("io error");
                    writeln!(writer, "").expect("io error");
                }
            },
            _ => { writeln!(writer, "").expect("io error"); }
        }

        writeln!(writer, "}}").expect("io error");
    }

    
    pub fn to_json<W: Write, F: Fn(&D) -> Value, RF: Fn(&mut W) -> ()>(&self, fmt_func: F, writer: &mut W) {
        self.to_json_rest(fmt_func, writer, None);
    }

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
}

pub struct NodeIter<'a, K: Kmer + 'a, D: Debug + 'a> {
    graph: &'a DebruijnGraph<K,D>,
    node_id: usize
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

/*
impl<'a, K: Kmer, D> fmt::Debug for Node<'a, K, D> where D:Debug {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Node: {}, Seq: {}, Exts:{:?}, Data: {:?}", self.node_id, self.sequence().to_string(), self.exts(), self.data())
    }
}
*/

impl<'a, K: Kmer, D> fmt::Debug for Node<'a, K, D> where D: Debug {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Node {{ id:{}, Exts: {:?}, L:{:?} R:{:?}, Seq: {:?} }}", self.node_id, self.exts(), self.l_edges(), self.r_edges(), self.sequence())
    }
}
