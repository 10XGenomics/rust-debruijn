use std::marker::PhantomData;
use fx::{FxHashMap, FxLMap, FxHasher};
use std::collections::VecDeque;
//use bit_set::BitSet;

use std::hash::BuildHasherDefault;
use std::ops::Index;

use Mer;
use Kmer;
use dna_string::{DnaString, DnaStringSlice};
use vmer::Vmer;
use Dir;
use Exts;


#[derive(Copy, Clone)]
pub enum ExtMode<K: Kmer> {
    Unique(K, Dir, Exts),
    Terminal(Exts),
}


struct PackedDnaStringSet {
    pub sequence: DnaString,
    pub start: Vec<usize>,
    pub length: Vec<u32>,
}

impl<'a> PackedDnaStringSet {
    fn get(&'a self, i: usize) -> DnaStringSlice<'a> {
        DnaStringSlice {
            dna_string: &self.sequence,
            start: self.start[i],
            length: self.length[i] as usize,
        }
    }

    fn len(&self) -> usize {
        self.start.len()
    }

    fn add<S: IntoIterator<Item=u8>>(&mut self, sequence: S) {
        let start = self.sequence.len();
        self.start.push(start);

        let mut length = 0;
        for b in sequence {
            self.sequence.push(b);
            length += 1;
        }
        self.length.push(length as u32);
    }
}


struct BaseGraph<K, D> {
    pub sequences: PackedDnaStringSet,
    pub exts: Vec<Exts>,
    pub data: Vec<D>,
    phantom: PhantomData<K>,
}

impl<K, D> BaseGraph<K, D> {
    pub fn len(&self) -> usize {
        self.sequences.len()
    }
}

impl<K: Kmer, D> BaseGraph<K,D> {
    pub fn add<S: IntoIterator<Item=u8>>(&mut self, sequence: S, exts: Exts, data: D) {
        self.sequences.add(sequence);
        self.exts.push(exts);
        self.data.push(data);
    }

    pub fn finish(self) -> DebruijnGraph<K,D> {

        let mut left_sort: Vec<u32> = Vec::with_capacity(self.len());
        let mut right_sort: Vec<u32> = Vec::with_capacity(self.len());
        for i in 0 .. self.len()
        {
            left_sort.push(i as u32);
            right_sort.push(i as u32);
        }

        left_sort.sort_by_key( |idx| -> K { self.sequences.get(*idx as usize).first_kmer() });
        right_sort.sort_by_key(|idx| -> K { self.sequences.get(*idx as usize).last_kmer() });

        DebruijnGraph {
            base: self,
            left_order: left_sort,
            right_order: right_sort,
        }
    }
}

pub struct DebruijnGraph<K, D> {
    base: BaseGraph<K, D>,
    left_order: Vec<u32>,
    right_order: Vec<u32>,
}

impl<K: Kmer, D> DebruijnGraph<K, D> {
    pub fn get_node<'a>(&'a self, node_id: usize) -> Node<'a, K, D> {
        Node {
            node_id: node_id,
            graph: self,
            l_edges: self.find_edges(node_id, Dir::Left),
            r_edges: self.find_edges(node_id, Dir::Right),
        }
    }

    pub fn find_edges(&self, node_id: usize, dir: Dir) -> Vec<(usize, Dir, bool)> {

        let exts = self.base.exts[node_id];
        let sequence = self.base.sequences.get(node_id);
        let kmer: K = sequence.term_kmer(dir);
        let mut edges = Vec::new();

        for i in 0..4
        {
            if exts.has_ext(dir, i) {
                let link = self.find_link(kmer.extend(i, dir), dir).expect("missing link");
                edges.push(link);
            }
        }

        edges
    }

    pub fn search_kmer(&self, kmer: K, side: Dir) -> Option<usize> {
        match side {
            Dir::Left => {
                let pos = self.left_order.binary_search_by_key(&kmer,
                    |idx| self.base.sequences.get(*idx as usize).first_kmer());

                match pos {
                    Ok(idx) => Some(self.left_order[idx] as usize),
                    _ => None,
                }
            },
            Dir::Right => {
                let pos = self.right_order.binary_search_by_key(&kmer,
                    |idx| self.base.sequences.get(*idx as usize).last_kmer());
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
}


// Unbranched edge in the DeBruijn graph
pub struct Node<'a, K: Kmer + 'a, D: 'a> {
    pub node_id: usize,
    pub graph: &'a DebruijnGraph<K,D>,
    pub l_edges: Vec<(usize, Dir, bool)>,
    pub r_edges: Vec<(usize, Dir, bool)>,
}

impl<'a, K: Kmer, D> Node<'a, K, D> {
    pub fn sequence(&self) -> DnaStringSlice<'a> {
        self.graph.base.sequences.get(self.node_id)
    }

    pub fn data(&self) -> &'a D {
        &self.graph.base.data[self.node_id]
    }

    pub fn exts(&self) -> Exts {
        self.graph.base.exts[self.node_id]
    }
}



pub struct PathCompression<K: Kmer, V:Vmer<K>, D, B: Fn(D,D) -> bool, R: Fn(D,D)->D> {
    allow_rc: bool,
    k: PhantomData<K>,
    v: PhantomData<V>,
    d: PhantomData<D>,
    break_fn: B,
    reduce: R,
}

type Seq = Vec<u8>;

pub trait BuildGraph<K: Kmer, V:Vmer<K>, D> {
    // Basic progression
    fn kmers_from_sequences(seqs: Vec<(Seq, D)>) -> Vec<(K, Exts, D)>;
    fn spaths_from_kmers(kmers: Vec<(K, Exts, D)>) -> Vec<(V, Exts, D)>;
    fn graph_from_spaths(spaths: Vec<(V, Exts, D)>) -> DebruijnGraph<K, D>;

    // Derived methods
    fn spaths_from_kmers_dict<I: Index<K, Output=(Exts, D)>>(kmers: I) -> Vec<(V, Exts, D)>;
    fn build_graph<I: Index<K, Output=Exts>>(kmers: I) -> DebruijnGraph<K, D>;
}

/// Compression of paths in Debruijn graph
impl<K: Kmer, V:Vmer<K>, D, B: Fn(D,D) -> bool, R: Fn(D,D)->D> PathCompression<K, V, D, B, R> {

    /// Attempt to extend kmer v in direction dir. Return:
    ///  - Unique(nextKmer, nextDir) if a single unique extension
    ///    is possible.  nextDir indicates the direction to extend nextMker
    ///    to preserve the direction of the extension.
    /// - Term(ext) no unique extension possible, indicating the extensions at this end of the line
    fn try_extend_kmer(&self, available_kmers: &FxLMap<K, ()>,
                    kmer_exts: &FxHashMap<K, Exts>,
                    v: K,
                    dir: Dir)
                    -> ExtMode<K> {
        let exts = kmer_exts.get(&v).expect("didn't have kmer");
        if exts.num_ext_dir(dir) != 1 || v == v.rc() {
            ExtMode::Terminal(exts.single_dir(dir))
        } else {
            // Get the next kmer
            let ext_base = exts.get_unique_extension(dir).expect("should be unique");

            let mut next_kmer = v.extend(ext_base, dir);
            let mut do_flip = false;

            if self.allow_rc {
                let flip_rc = next_kmer.min_rc_flip();
                do_flip = flip_rc.1;
                next_kmer = flip_rc.0;
            }

            let is_palindrome = next_kmer == next_kmer.rc();

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
            let next_kmer_exts = kmer_exts.get(&next_kmer).expect("must have kmer");
            let incoming_count = next_kmer_exts.num_ext_dir(new_incoming_dir);
            let outgoing_exts = next_kmer_exts.single_dir(new_incoming_dir.flip());

            if incoming_count == 0 && !is_palindrome {
                panic!("unreachable");
            } else if incoming_count == 1 && !is_palindrome {
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
    ///
    /// Also return the extensions at the end of this line.
    /// Sub-lines break if their extensions are not available in this shard
    #[inline(never)]
    fn extend_kmer(&self, kmer_exts: &FxHashMap<K, Exts>,
                    available_kmers: &mut FxLMap<K, ()>,
                    kmer: K,
                    start_dir: Dir,
                    max_dist: usize)
                    -> (Vec<(K, Dir)>, Exts) {

        let mut current_dir = start_dir;
        let mut current_kmer = kmer;
        let mut path = Vec::new();

        let mut final_exts: Exts; // must get set below

        if max_dist == 0 {
            let first_exts = kmer_exts.get(&current_kmer).expect("didn't have kmer");
            return (path, first_exts.single_dir(start_dir));
        }

        let _ = available_kmers.remove(&kmer);

        loop {
            let ext_result = self.try_extend_kmer(available_kmers, kmer_exts, current_kmer, current_dir);

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

        (path, final_exts)
    }


    /// Build the edge surrounding a kmer
    #[inline(never)]
    fn build_node(&self, kmer_exts: &FxHashMap<K, Exts>,
                available_kmers: &mut FxLMap<K, ()>,
                seed: K)
                -> (V, Exts) {

        let (l_path, l_ext) = self.extend_kmer(kmer_exts, available_kmers, seed, Dir::Left, V::max_len() - K::k());
        let (r_path, r_ext) = self.extend_kmer(kmer_exts,
                                        available_kmers,
                                        seed,
                                        Dir::Right,
                                        V::max_len() - K::k() - l_path.len());

        let mut edge_seq = VecDeque::new();
        for i in 0..K::k() {
            edge_seq.push_back(seed.get(i));
        }

        // Add on the left path
        for &(next_kmer, dir) in l_path.iter() {
            let kmer = match dir {
                Dir::Left => next_kmer,
                Dir::Right => next_kmer.rc(),
            };

            edge_seq.push_front(kmer.get(0))
        }

        // Add on the right path
        for &(next_kmer, dir) in r_path.iter() {
            let kmer = match dir {
                Dir::Left => next_kmer.rc(),
                Dir::Right => next_kmer,
            };

            edge_seq.push_back(kmer.get(K::k() - 1))
        }

        let left_extend = match l_path.last() {
            None => l_ext,
            Some(&(_, Dir::Left)) => l_ext,
            Some(&(_, Dir::Right)) => l_ext.complement(),
        };

        let right_extend = match r_path.last() {
            None => r_ext,
            Some(&(_, Dir::Left)) => r_ext.complement(),
            Some(&(_, Dir::Right)) => r_ext,
        };

        let mut edge = V::new(edge_seq.len());
        for (pos, base) in edge_seq.iter().enumerate() {
            edge.set_mut(pos, *base);
        }
        (edge, Exts::from_single_dirs(left_extend, right_extend))
    }

    /// Build all sedges until all kmers are exhausted
    #[inline(never)]
    pub fn build_sedges(&self, kmer_exts: &FxHashMap<K, Exts>) -> Vec<(V, Exts)> {
        let mut edges: Vec<(V, Exts)> = Vec::new();

        let h = BuildHasherDefault::<FxHasher>::default();
        let mut available_kmers: FxLMap<K, ()> = FxLMap::with_capacity_and_hasher(kmer_exts.len(), h);
        available_kmers.extend(kmer_exts.keys().map(|k| (k.clone(), ())));

        loop {
            match available_kmers.front() {
                Some((&k, _)) => {
                    let sedge: (V, Exts) = self.build_sedge(kmer_exts, &mut available_kmers, k);
                    edges.push(sedge)
                }
                None => break,
            }
        }

        edges
    }


    /// Build the edge surrounding a kmer
    #[inline(never)]
    fn build_sedge(&self, kmer_exts: &FxHashMap<K, Exts>,
                available_kmers: &mut FxLMap<K, ()>,
                seed: K)
                -> (V, Exts) {
        let (l_path, l_ext) = self.extend_sedge(kmer_exts, available_kmers, seed, Dir::Left, V::max_len() - K::k());
        let (r_path, r_ext) = self.extend_sedge(kmer_exts,
                                        available_kmers,
                                        seed,
                                        Dir::Right,
                                        V::max_len() - K::k() - l_path.len());

        let mut edge_seq = VecDeque::new();
        for i in 0..K::k() {
            edge_seq.push_back(seed.get(i));
        }

        // Add on the left path
        for &(next_kmer, dir) in l_path.iter() {
            let kmer = match dir {
                Dir::Left => next_kmer,
                Dir::Right => next_kmer.rc(),
            };

            edge_seq.push_front(kmer.get(0))
        }

        // Add on the right path
        for &(next_kmer, dir) in r_path.iter() {
            let kmer = match dir {
                Dir::Left => next_kmer.rc(),
                Dir::Right => next_kmer,
            };

            edge_seq.push_back(kmer.get(K::k() - 1))
        }

        let left_extend = match l_path.last() {
            None => l_ext,
            Some(&(_, Dir::Left)) => l_ext,
            Some(&(_, Dir::Right)) => l_ext.complement(),
        };

        let right_extend = match r_path.last() {
            None => r_ext,
            Some(&(_, Dir::Left)) => r_ext.complement(),
            Some(&(_, Dir::Right)) => r_ext,
        };

        let mut edge = V::new(edge_seq.len());
        for (idx, base) in edge_seq.iter().enumerate() {
            edge.set_mut(idx, *base);
        }
        (edge, Exts::from_single_dirs(left_extend, right_extend))
    }

    /// Build the maximal line starting at kmer in direction dir, at most max_dist long.
    ///
    /// Also return the extensions at the end of this line.
    /// Sub-lines break if their extensions are not available in this shard
    #[inline(never)]
    fn extend_sedge(&self, kmer_exts: &FxHashMap<K, Exts>,
                    available_kmers: &mut FxLMap<K, ()>,
                    kmer: K,
                    start_dir: Dir,
                    max_dist: usize)
                    -> (Vec<(K, Dir)>, Exts) {
        let mut current_dir = start_dir;
        let mut current_kmer = kmer;
        let mut path = Vec::new();
        let mut final_exts: Exts; // must get set below

        if max_dist == 0 {
            let first_exts = kmer_exts.get(&current_kmer).expect("didn't have kmer");
            return (path, first_exts.single_dir(start_dir));
        }

        let _ = available_kmers.remove(&kmer);

        loop {
            let ext_result = self.try_extend_kmer(available_kmers, kmer_exts, current_kmer, current_dir);

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

        (path, final_exts)
    }
}
