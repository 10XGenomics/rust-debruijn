use std::marker::PhantomData;
use fx::{FxLMap, FxHasher};
use std::collections::VecDeque;
use bit_set::BitSet;
use smallvec::SmallVec;
use std::fmt::Debug;

use std::hash::BuildHasherDefault;

type SmallVec4<T> = SmallVec<[T; 4]>;

use Kmer;
use Vmer;
use Dir;
use Exts;
use paths::{DebruijnGraph, BaseGraph};
use dna_string::DnaString;


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


pub struct PathCompression<K: Kmer, D, S: CompressionSpec<D>> {
    stranded: bool,
    k: PhantomData<K>,
    d: PhantomData<D>,
    spec: S,
}

/// Compression of paths in Debruijn graph
impl<K: Kmer, D: Clone + Debug, S: CompressionSpec<D>> PathCompression<K, D, S> {

    pub fn new(stranded: bool, spec: S) -> Self {
        PathCompression {
            stranded: stranded,
            spec: spec,
            k: PhantomData,
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

        if exts.num_ext_dir(dir) != 1 || (!self.stranded && kmer.is_palindrome()) {
            ExtMode::Terminal(exts.single_dir(dir))
        } else {
            // Get the next kmer
            let ext_base = exts.get_unique_extension(dir).expect("should be unique");

            let mut next_kmer = kmer.extend(ext_base, dir);
            
            let mut do_flip = false;

            if !self.stranded {
                let flip_rc = next_kmer.min_rc_flip();
                do_flip = flip_rc.1;
                next_kmer = flip_rc.0;
            }

            let next_dir = dir.cond_flip(do_flip);
            let is_palindrome = !self.stranded && next_kmer.is_palindrome();

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
                   path: &mut Vec<(K, Dir)>)
                   -> Exts {

        let mut current_dir = start_dir;
        let mut current_kmer = kmer;
        path.clear();

        let mut final_exts: Exts; // must get set below

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

        (Exts::from_single_dirs(left_extend, right_extend), node_data)
    }

    /// Build all sedges until all kmers are exhausted
    #[inline(never)]
    pub fn build_nodes(&self, kmer_exts: &Vec<(K, (Exts, D))>) -> BaseGraph<K,D> {

        let h = BuildHasherDefault::<FxHasher>::default();
        let mut available_kmers: FxLMap<K, ()> = FxLMap::with_capacity_and_hasher(kmer_exts.len(), h);
        available_kmers.extend(kmer_exts.iter().map(|&(k, _)| (k.clone(), ())));

        // Path-compressed De Bruijn graph will be created here
        let mut graph = BaseGraph::new(self.stranded);

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
}


pub struct Simplify<'a, K: 'a + Kmer, D:'a, S: CompressionSpec<D>> {
    stranded: bool,
    d: PhantomData<D>,
    spec: S,
    available_nodes: BitSet,
    graph: &'a DebruijnGraph<K, D>,
}

impl<'a, K:Kmer, D: Debug + Clone, S: CompressionSpec<D>> Simplify<'a, K, D, S> {

    #[inline(never)]
    fn try_extend_node(&mut self, node: usize, dir: Dir) -> ExtModeNode {

        let node = self.graph.get_node(node);
        let bases = node.sequence();
        let exts = node.exts();

        if exts.num_ext_dir(dir) != 1 || (!self.stranded && node.len() == K::k() && Vmer::<K>::get_kmer(&bases,0).is_palindrome()) {
            ExtModeNode::Terminal(exts.single_dir(dir))
        } else {
            // Get the next kmer
            let ext_base = exts.get_unique_extension(dir).expect("should be unique");
            let end_kmer: K = bases.term_kmer(dir);

            let next_kmer = end_kmer.extend(ext_base, dir);
            let (next_node_id, next_side_incoming, rc) = 
                match self.graph.find_link(next_kmer, dir) {
                    Some(e) => e,
                    None => {
                        println!("dir: {:?}, Lmer: {:?}, exts: {:?}", dir, bases, exts);
                        println!("end kmer: {:?}", end_kmer);
                        println!("No kmer: {:?}", next_kmer);
                        println!("rc: {:?}", next_kmer.min_rc());
                        panic!(format!("No kmer: {:?}", next_kmer))
                    }
                };

            let next_node = self.graph.get_node(next_node_id);
            let next_exts = next_node.exts();

            let consistent = (next_node.len() == K::k()) ||
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
                            next_node.sequence(),
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

            if !self.available_nodes.contains(next_node_id) || (!self.stranded && next_kmer.is_palindrome())  {
                // Next kmer isn't in this partition,
                // or we've already used it,
                // or it's palindrom and we are not stranded
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

    /// Generate complete unbranched edges
    fn extend_node(&mut self, start_node: usize, start_dir: Dir) -> (Vec<(usize, Dir)>, Exts) {

        let mut current_dir = start_dir;
        let mut current_node = start_node;
        let mut path = Vec::new();
        let final_exts: Exts; // must get set below

        self.available_nodes.remove(start_node);

        loop {
            let ext_result = self.try_extend_node(current_node, current_dir);

            match ext_result {
                ExtModeNode::Unique(next_node, next_dir_outgoing, _) => {
                    let next_dir_incoming = next_dir_outgoing.flip();
                    path.push((next_node, next_dir_incoming));
                    self.available_nodes.remove(next_node);
                    current_node = next_node;
                    current_dir = next_dir_outgoing;
                }
                ExtModeNode::Terminal(ext) => {
                    final_exts = ext;
                    break;
                }
            }
        }

        (path, final_exts)
    }


    // Determine the sequence and extensions of the maximal unbranched
    // edge, centered around the given edge number
    #[inline(never)]
    fn build_node(&mut self, seed_node: usize)
                -> (DnaString, Exts, VecDeque<(usize, Dir)>, D) {

        let (l_path, l_ext) = self.extend_node(seed_node, Dir::Left);
        let (r_path, r_ext) = self.extend_node(seed_node, Dir::Right);

        // Stick together edge chunks to get full edge sequence
        let mut node_path = VecDeque::new();

        let mut node_data: D = self.graph.get_node(seed_node).data().clone();
        node_path.push_back((seed_node, Dir::Left));

        // Add on the left path
        for &(next_node, incoming_dir) in l_path.iter() {
            node_path.push_front((next_node, incoming_dir.flip()));
            node_data = self.spec.reduce(node_data, self.graph.get_node(next_node).data());
        }

        // Add on the right path
        for &(next_node, incoming_dir) in r_path.iter() {
            node_path.push_back((next_node, incoming_dir));
            node_data = self.spec.reduce(node_data, self.graph.get_node(next_node).data());
        }

        let left_extend = match l_path.last() {
            None => l_ext,
            Some(&(_, Dir::Left)) => l_ext.complement(),
            Some(&(_, Dir::Right)) => l_ext,
        };

        let right_extend = match r_path.last() {
            None => r_ext,
            Some(&(_, Dir::Left)) => r_ext,
            Some(&(_, Dir::Right)) => r_ext.complement(),
        };

        let path_seq = self.graph.sequence_of_path(node_path.iter());

        // return sequence and extensions
        (path_seq, Exts::from_single_dirs(left_extend, right_extend), node_path, node_data)
    }


    /// Simplify a compressed Debruijn graph by merging adjacent unbranched nodes, and optionally
    /// censoring some nodes
    pub fn simplify(mut old_graph: DebruijnGraph<K,D>, censor_nodes: Option<Vec<usize>>, stranded: bool, compression: S) -> DebruijnGraph<K, D> {


        let n_nodes = old_graph.len();
        

        let mut available_nodes = BitSet::with_capacity(n_nodes);
        for i in 0 .. n_nodes {
            available_nodes.insert(i);
        }

        match censor_nodes {
            Some(c) => {
                for censor in c {
                    available_nodes.remove(censor);
                }
            },
            None => (),
        }

        old_graph.fix_exts(Some(&available_nodes));

        let mut simplify = 
            Simplify {
                spec: compression,
                stranded: stranded,
                graph: &old_graph,
                available_nodes: available_nodes,
                d: PhantomData,
            };

        // FIXME -- clarify requirements around state of extensions
        let mut graph = BaseGraph::new(stranded);

        for node_counter in 0 .. n_nodes {

            if simplify.available_nodes.contains(node_counter) {
                let (seq, exts, _, data) = simplify.build_node(node_counter);
                graph.add(&seq, exts, data);
            }
        }

        // We will have some hanging exts due to
        let mut dbg = graph.finish();
        dbg.fix_exts(None);
        //assert!(dbg.is_compressed() == None);
        dbg

        
    }
}
