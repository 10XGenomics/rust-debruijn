// Copyright 2017 10x Genomics

//! # debruijn: a De Bruijn graph library for DNA seqeunces in Rust.
//! This library provides tools for efficient construction DeBruijn graphs (dBG)
//! from DNA sequences, tracking arbitrary metadata associated with kmers in the
//! graph, and performing path-compression of unbranched graph paths to improve
//! speed and reduce memory consumption.
//!
//! Most applications of `debruijn` will follow this general workflow:
//! 1. You generate a set of sequences to make a dBG from.
//! 2. You pass those sequences to the `filter_kmers` function, which converts the sequences into kmers, while tracking 'metadata' about each kmer in a very customizable way. The metadata could be read count, a set of colors, a set of read counts split by haplotype, a UMI count, etc.
//! 3. The the library will convert the kmers to a compressed dBG. You can also customize the rules for how to compress the dBG and how to 'combine' the per-kmer metadata.
//!
//! Then you can use the final compressed dBG how you like. There are some methods for simplifying and re-building the  graph, but those could be developed more.
//!
//! ## Examples
//! - [Local phased SV assembly tool in our Long Ranger package](https://github.com/10XGenomics/longranger/blob/master/lib/pvc/src/asm_caller.rs#L205)
//! - [Single-cell VDJ assember](https://github.com/10XGenomics/cellranger/blob/master/lib/rust/vdj_asm/src/asm.rs#L191)
//! - [Build a colored, compressed dBG of a transcriptome reference](https://github.com/10XGenomics/rust-pseudoaligner/blob/master/src/build_index.rs#L40)
//!
//! All the data structures in debruijn-rs are specialized to the 4 base DNA alphabet,
//! and use 2-bit packed encoding of base-pairs into integer types, and efficient methods for
//! reverse complement, enumerating kmers from longer sequences, and transfering data between
//! sequences.
//!
//! ## Encodings
//! Most methods for ingesting sequence data into the library have a form named 'bytes',
//! which expects bases encoded as the integers 0,1,2,3, and a separate form names 'ascii',
//! which expects bases encoded as the ASCII letters A,C,G,T.

use serde_derive::{Deserialize, Serialize};
use std::fmt;
use std::hash::Hash;

pub mod clean_graph;
pub mod compression;
pub mod dna_string;
pub mod filter;
pub mod graph;
pub mod kmer;
pub mod msp;
pub mod vmer;

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
mod bitops_avx2;

#[cfg(test)]
pub mod test;

/// Convert a 2-bit representation of a base to a char
#[inline]
pub fn bits_to_ascii(c: u8) -> u8 {
    match c {
        0u8 => 'A' as u8,
        1u8 => 'C' as u8,
        2u8 => 'G' as u8,
        3u8 => 'T' as u8,
        _ => 'X' as u8,
    }
}

/// Convert an ASCII-encoded DNA base to a 2-bit representation
#[inline]
pub fn base_to_bits(c: u8) -> u8 {
    match c {
        b'A' | b'a' => 0u8,
        b'C' | b'c' => 1u8,
        b'G' | b'g' => 2u8,
        b'T' | b't' => 3u8,
        _ => 0u8,
    }
}

#[inline]
pub fn dna_only_base_to_bits(c: u8) -> Option<u8> {
    match c {
        b'A' | b'a' => Some(0u8),
        b'C' | b'c' => Some(1u8),
        b'G' | b'g' => Some(2u8),
        b'T' | b't' => Some(3u8),
        _ => None,
    }
}

/// Convert an ASCII-encoded DNA base to a 2-bit representation
#[inline]
pub fn is_valid_base(c: u8) -> bool {
    match c {
        b'A' | b'C' | b'G' | b'T' => true,
        b'a' | b'c' | b'g' | b't' => true,
        _ => false,
    }
}

/// Convert a 2-bit representation of a base to a char
#[inline]
pub fn bits_to_base(c: u8) -> char {
    match c {
        0u8 => 'A',
        1u8 => 'C',
        2u8 => 'G',
        3u8 => 'T',
        _ => 'X',
    }
}

/// The complement of a 2-bit encoded base
#[inline(always)]
pub fn complement(base: u8) -> u8 {
    (!base) & 0x3u8
}

/// Trait for interacting with DNA sequences
pub trait Mer: Sized + fmt::Debug {
    /// Length of DNA sequence
    fn len(&self) -> usize;

    /// Get 2-bit encoded base at position `pos`
    fn get(&self, pos: usize) -> u8;

    /// Set base at `pos` to 2-bit encoded base `val`
    fn set_mut(&mut self, pos: usize, val: u8);

    /// Set `nbases` positions in the sequence, starting at `pos`.
    /// Values must  be packed into the upper-most bits of `value`.
    fn set_slice_mut(&mut self, pos: usize, nbases: usize, value: u64);

    /// Return a new object containing the reverse complement of the sequence
    fn rc(&self) -> Self;

    /// Iterate over the bases in the sequence
    fn iter<'a>(&'a self) -> MerIter<'a, Self> {
        MerIter {
            sequence: self,
            i: 0,
        }
    }

    /// Count the number of A/T bases in the kmer
    fn at_count(&self) -> u32 {
        let mut count = 0;
        for i in 0..self.len() {
            let base = self.get(i);
            if base == 0 || base == 3 {
                count += 1;
            }
        }
        count
    }

    /// Count the number of G/C bases in the kmer
    fn gc_count(&self) -> u32 {
        let mut count = 0;
        for i in 0..self.len() {
            let base = self.get(i);
            if base == 1 || base == 2 {
                count += 1;
            }
        }
        count
    }
}

/// Iterator over bases of a DNA sequence (bases will be unpacked into bytes).
pub struct MerIter<'a, M: 'a + Mer> {
    sequence: &'a M,
    i: usize,
}

impl<'a, M: 'a + Mer> Iterator for MerIter<'a, M> {
    type Item = u8;

    fn next(&mut self) -> Option<u8> {
        if self.i < self.sequence.len() {
            let value = self.sequence.get(self.i);
            self.i += 1;
            Some(value)
        } else {
            None
        }
    }
}

/// Encapsulates a Kmer sequence with statically known K.
pub trait Kmer: Mer + Sized + Copy + PartialEq + PartialOrd + Eq + Ord + Hash {
    /// Create a Kmer initialized to all A's
    fn empty() -> Self;

    /// K value for this concrete type.
    fn k() -> usize;

    /// Return the rank of this kmer in an lexicographic ordering of all kmers
    /// E.g. 'AAAA' -> 0, 'AAAT' -> 1, etc. This will panic if K > 32.
    fn to_u64(&self) -> u64;

    // Construct a kmer from the given lexicographic rank of the kmer.
    // If K > 32, the leads bases will be A's.
    fn from_u64(value: u64) -> Self;

    /// Add the base `v` to the left side of the sequence, and remove the rightmost base
    fn extend_left(&self, v: u8) -> Self;

    /// Add the base `v` to the right side of the sequence, and remove the leftmost base
    fn extend_right(&self, v: u8) -> Self;

    /// Add the base `v` to the side of sequence given by `dir`, and remove a base at the opposite side
    fn extend(&self, v: u8, dir: Dir) -> Self {
        match dir {
            Dir::Left => self.extend_left(v),
            Dir::Right => self.extend_right(v),
        }
    }

    /// Generate all the extension of this sequence given by `exts` in direction `Dir`
    fn get_extensions(&self, exts: Exts, dir: Dir) -> Vec<Self> {
        let ext_bases = exts.get(dir);
        ext_bases
            .iter()
            .map(|b| self.extend(b.clone(), dir))
            .collect()
    }

    /// Return the minimum of the kmer and it's reverse complement, and a flag indicating if sequence was flipped
    fn min_rc_flip(&self) -> (Self, bool) {
        let rc = self.rc();
        if *self < rc {
            (self.clone(), false)
        } else {
            (rc, true)
        }
    }

    // Return the minimum of the kmer and it's reverse complement
    fn min_rc(&self) -> Self {
        let rc = self.rc();
        if *self < rc {
            self.clone()
        } else {
            rc
        }
    }

    /// Test if this Kmer and it's reverse complement are the same
    fn is_palindrome(&self) -> bool {
        self.len() % 2 == 0 && *self == self.rc()
    }

    /// Create a Kmer from the first K bytes of `bytes`, which must be encoded as the integers 0-4.
    fn from_bytes(bytes: &[u8]) -> Self {
        if bytes.len() < Self::k() {
            panic!("bytes not long enough to form kmer")
        }

        let mut k0 = Self::empty();

        for i in 0..Self::k() {
            k0.set_mut(i, bytes[i])
        }

        k0
    }

    /// Create a Kmer from the first K bytes of `bytes`, which must be encoded as ASCII letters A,C,G, or T.
    fn from_ascii(bytes: &[u8]) -> Self {
        if bytes.len() < Self::k() {
            panic!("bytes not long enough to form kmer")
        }

        let mut k0 = Self::empty();

        for i in 0..Self::k() {
            k0.set_mut(i, base_to_bits(bytes[i]))
        }

        k0
    }

    /// Return String containing Kmer sequence
    fn to_string(&self) -> String {
        let mut s = String::new();
        for pos in 0..self.len() {
            s.push(bits_to_base(self.get(pos)))
        }
        s
    }

    /// Generate vector of all kmers contained in `str` encoded as 0-4.
    fn kmers_from_bytes(str: &[u8]) -> Vec<Self> {
        let mut r = Vec::new();

        if str.len() < Self::k() {
            return r;
        }

        let mut k0 = Self::empty();

        for i in 0..Self::k() {
            k0.set_mut(i, str[i]);
        }

        r.push(k0.clone());

        for i in Self::k()..str.len() {
            k0 = k0.extend_right(str[i]);
            r.push(k0.clone());
        }

        r
    }

    /// Generate vector of all kmers contained in `str`, encoded as ASCII ACGT.
    fn kmers_from_ascii(str: &[u8]) -> Vec<Self> {
        let mut r = Vec::new();

        if str.len() < Self::k() {
            return r;
        }

        let mut k0 = Self::empty();

        for i in 0..Self::k() {
            k0.set_mut(i, base_to_bits(str[i]));
        }

        r.push(k0.clone());

        for i in Self::k()..str.len() {
            k0 = k0.extend_right(base_to_bits(str[i]));
            r.push(k0.clone());
        }

        r
    }
}

/// An immutable interface to a Mer sequence.
pub trait MerImmut: Mer + Clone {
    fn set(&self, pos: usize, val: u8) -> Self {
        let mut new = self.clone();
        new.set_mut(pos, val);
        new
    }

    fn set_slice(&self, pos: usize, nbases: usize, bits: u64) -> Self {
        let mut new = self.clone();
        new.set_slice_mut(pos, nbases, bits);
        new
    }
}

impl<T> MerImmut for T where T: Mer + Clone {}

/// A DNA sequence with run-time variable length, up to a statically known maximum length
pub trait Vmer: Mer + PartialEq + Eq {
    /// Create a new sequence with length `len`, initialized to all A's
    fn new(len: usize) -> Self;

    /// Maximum sequence length that can be stored in this type
    fn max_len() -> usize;

    /// Create a Vmer from a sequence of bytes
    fn from_slice(seq: &[u8]) -> Self {
        let mut vmer = Self::new(seq.len());
        for i in 0..seq.len() {
            vmer.set_mut(i, seq[i]);
        }

        vmer
    }

    /// Efficiently extract a Kmer from the sequence
    fn get_kmer<K: Kmer>(&self, pos: usize) -> K;

    /// Get the first Kmer from the sequence
    fn first_kmer<K: Kmer>(&self) -> K {
        self.get_kmer(0)
    }

    /// Get the last kmer in the sequence
    fn last_kmer<K: Kmer>(&self) -> K {
        self.get_kmer(self.len() - K::k())
    }

    /// Get the terminal kmer of the sequence, on the both side of the sequence
    fn both_term_kmer<K: Kmer>(&self) -> (K, K) {
        (self.first_kmer(), self.last_kmer())
    }

    /// Get the terminal kmer of the sequence, on the side of the sequence given by dir
    fn term_kmer<K: Kmer>(&self, dir: Dir) -> K {
        match dir {
            Dir::Left => self.first_kmer(),
            Dir::Right => self.last_kmer(),
        }
    }

    /// Iterate over the kmers in the sequence
    fn iter_kmers<K: Kmer>(&self) -> KmerIter<'_, K, Self> {
        let kmer = if self.len() >= K::k() {
            self.first_kmer()
        } else {
            // Default kmer, will not be used
            K::empty()
        };

        KmerIter {
            bases: self,
            kmer: kmer,
            pos: K::k(),
        }
    }

    /// Iterate over the kmers and their extensions, given the extensions of the whole sequence
    fn iter_kmer_exts<K: Kmer>(&self, seq_exts: Exts) -> KmerExtsIter<'_, K, Self> {
        let kmer = if self.len() >= K::k() {
            self.first_kmer()
        } else {
            // Default kmer, will not be used
            K::empty()
        };

        KmerExtsIter {
            bases: self,
            exts: seq_exts,
            kmer: kmer,
            pos: K::k(),
        }
    }
}

/// A newtype wrapper around a `Vec<u8>` with implementations
// of the `Mer` and `Vmer` traits.
#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd)]
pub struct DnaBytes(pub Vec<u8>);

impl Mer for DnaBytes {
    fn len(&self) -> usize {
        self.0.len()
    }

    fn get(&self, pos: usize) -> u8 {
        self.0[pos]
    }

    /// Set base at `pos` to 2-bit encoded base `val`
    fn set_mut(&mut self, pos: usize, val: u8) {
        self.0[pos] = val
    }

    /// Set `nbases` positions in the sequence, starting at `pos`.
    /// Values must  be packed into the upper-most bits of `value`.
    fn set_slice_mut(&mut self, _pos: usize, _nbases: usize, _value: u64) {
        unimplemented!();
        //for i in pos .. (pos + nbases) {
        //
        //}
    }

    /// Return a new object containing the reverse complement of the sequence
    fn rc(&self) -> Self {
        unimplemented!();
    }
}

impl Vmer for DnaBytes {
    /// Create a new sequence with length `len`, initialized to all A's
    fn new(len: usize) -> Self {
        DnaBytes(vec![0u8; len])
    }

    /// Maximum sequence length that can be stored in this type
    fn max_len() -> usize {
        1 << 48
    }

    /// Efficiently extract a Kmer from the sequence
    fn get_kmer<K: Kmer>(&self, pos: usize) -> K {
        K::from_bytes(&self.0[pos..pos + K::k()])
    }
}

/// A newtype wrapper around a `&[u8]` with implementations
// of the `Mer` and `Vmer` traits.
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct DnaSlice<'a>(pub &'a [u8]);

impl<'a> Mer for DnaSlice<'a> {
    fn len(&self) -> usize {
        self.0.len()
    }

    fn get(&self, pos: usize) -> u8 {
        self.0[pos]
    }

    /// Set base at `pos` to 2-bit encoded base `val`
    fn set_mut(&mut self, _pos: usize, _val: u8) {
        unimplemented!()
    }

    /// Set `nbases` positions in the sequence, starting at `pos`.
    /// Values must  be packed into the upper-most bits of `value`.
    fn set_slice_mut(&mut self, _pos: usize, _nbases: usize, _value: u64) {
        unimplemented!();
        //for i in pos .. (pos + nbases) {
        //
        //}
    }

    /// Return a new object containing the reverse complement of the sequence
    fn rc(&self) -> Self {
        unimplemented!();
    }
}

impl<'a> Vmer for DnaSlice<'a> {
    /// Create a new sequence with length `len`, initialized to all A's
    fn new(_len: usize) -> Self {
        unimplemented!();
    }

    /// Maximum sequence length that can be stored in this type
    fn max_len() -> usize {
        1 << 48
    }

    /// Efficiently extract a Kmer from the sequence
    fn get_kmer<K: Kmer>(&self, pos: usize) -> K {
        K::from_bytes(&self.0[pos..pos + K::k()])
    }
}

/// Direction of motion in a DeBruijn graph
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub enum Dir {
    Left,
    Right,
}

impl Dir {
    /// Return a fresh Dir with the opposite direction
    pub fn flip(&self) -> Dir {
        match *self {
            Dir::Left => Dir::Right,
            Dir::Right => Dir::Left,
        }
    }

    /// Return a fresh Dir opposite direction if do_flip == True
    pub fn cond_flip(&self, do_flip: bool) -> Dir {
        if do_flip {
            self.flip()
        } else {
            *self
        }
    }

    /// Pick between two alternatives, depending on the direction
    pub fn pick<T>(&self, if_left: T, if_right: T) -> T {
        match self {
            &Dir::Left => if_left,
            &Dir::Right => if_right,
        }
    }
}

/// Store single-base extensions for a DNA Debruijn graph.
///
/// 8 bits, 4 higher order ones represent extensions to the right, 4 lower order ones
/// represent extensions to the left. For each direction the bits (from lower order
/// to higher order) represent whether there exists an extension with each of the
/// letters A, C, G, T. So overall the bits are:
///  right   left
/// T G C A T G C A
#[derive(Eq, PartialEq, Copy, Clone, Ord, PartialOrd, Hash, Serialize, Deserialize)]
pub struct Exts {
    pub val: u8,
}

impl Exts {
    pub fn new(val: u8) -> Self {
        Exts { val: val }
    }

    pub fn empty() -> Exts {
        Exts { val: 0u8 }
    }

    pub fn from_single_dirs(left: Exts, right: Exts) -> Exts {
        Exts {
            val: (right.val << 4) | (left.val & 0xf),
        }
    }

    pub fn merge(left: Exts, right: Exts) -> Exts {
        Exts {
            val: left.val & 0x0f | right.val & 0xf0,
        }
    }

    pub fn add(&self, v: Exts) -> Exts {
        Exts {
            val: self.val | v.val,
        }
    }

    pub fn set(&self, dir: Dir, pos: u8) -> Exts {
        let shift = pos
            + match dir {
                Dir::Right => 4,
                Dir::Left => 0,
            };

        let new_val = self.val | (1u8 << shift);
        Exts { val: new_val }
    }

    #[inline]
    fn dir_bits(&self, dir: Dir) -> u8 {
        match dir {
            Dir::Right => self.val >> 4,
            Dir::Left => self.val & 0xf,
        }
    }

    pub fn get(&self, dir: Dir) -> Vec<u8> {
        let bits = self.dir_bits(dir);
        let mut v = Vec::new();
        for i in 0..4 {
            if bits & (1 << i) > 0 {
                v.push(i);
            }
        }

        v
    }

    pub fn has_ext(&self, dir: Dir, base: u8) -> bool {
        let bits = self.dir_bits(dir);
        (bits & (1 << base)) > 0
    }

    pub fn from_slice_bounds(src: &[u8], start: usize, length: usize) -> Exts {
        let l_extend = if start > 0 {
            1u8 << (src[start - 1])
        } else {
            0u8
        };
        let r_extend = if start + length < src.len() {
            1u8 << src[start + length]
        } else {
            0u8
        };

        Exts {
            val: (r_extend << 4) | l_extend,
        }
    }

    pub fn from_dna_string(src: &dna_string::DnaString, start: usize, length: usize) -> Exts {
        let l_extend = if start > 0 {
            1u8 << (src.get(start - 1))
        } else {
            0u8
        };
        let r_extend = if start + length < src.len() {
            1u8 << src.get(start + length)
        } else {
            0u8
        };

        Exts {
            val: (r_extend << 4) | l_extend,
        }
    }

    pub fn num_exts_l(&self) -> u8 {
        self.num_ext_dir(Dir::Left)
    }

    pub fn num_exts_r(&self) -> u8 {
        self.num_ext_dir(Dir::Right)
    }

    pub fn num_ext_dir(&self, dir: Dir) -> u8 {
        let e = self.dir_bits(dir);
        ((e & 1u8) >> 0) + ((e & 2u8) >> 1) + ((e & 4u8) >> 2) + ((e & 8u8) >> 3)
    }

    pub fn mk_left(base: u8) -> Exts {
        Exts::empty().set(Dir::Left, base)
    }

    pub fn mk_right(base: u8) -> Exts {
        Exts::empty().set(Dir::Right, base)
    }

    pub fn mk(left_base: u8, right_base: u8) -> Exts {
        Exts::merge(Exts::mk_left(left_base), Exts::mk_right(right_base))
    }

    pub fn get_unique_extension(&self, dir: Dir) -> Option<u8> {
        if self.num_ext_dir(dir) != 1 {
            None
        } else {
            let e = self.dir_bits(dir);
            for i in 0..4 {
                if (e & (1 << i)) > 0 {
                    return Some(i);
                }
            }

            None
        }
    }

    pub fn single_dir(&self, dir: Dir) -> Exts {
        match dir {
            Dir::Right => Exts { val: self.val >> 4 },
            Dir::Left => Exts {
                val: self.val & 0xfu8,
            },
        }
    }

    /// Complement the extension bases for each direction
    pub fn complement(&self) -> Exts {
        let v = self.val;

        // swap bits
        let mut r = (v & 0x55u8) << 1 | ((v >> 1) & 0x55u8);

        // swap pairs
        r = (r & 0x33u8) << 2 | ((r >> 2) & 0x33u8);
        Exts { val: r }
    }

    pub fn reverse(&self) -> Exts {
        let v = self.val;
        let r = (v & 0xf) << 4 | (v >> 4);
        Exts { val: r }
    }

    pub fn rc(&self) -> Exts {
        self.reverse().complement()
    }
}

impl fmt::Debug for Exts {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();

        for b in self.get(Dir::Left) {
            s.push(bits_to_base(b));
        }
        s.push('|');

        for b in self.get(Dir::Right) {
            s.push(bits_to_base(b));
        }

        write!(f, "{}", s)
    }
}

/// Iterate over the `Kmer`s of a DNA sequence efficiently
pub struct KmerIter<'a, K: Kmer, D>
where
    D: 'a,
{
    bases: &'a D,
    kmer: K,
    pos: usize,
}

impl<'a, K: Kmer, D: Mer> Iterator for KmerIter<'a, K, D> {
    type Item = K;

    #[inline]
    fn next(&mut self) -> Option<K> {
        if self.pos <= self.bases.len() {
            let retval = self.kmer;

            if self.pos < self.bases.len() {
                self.kmer = self.kmer.extend_right(self.bases.get(self.pos));
            }

            self.pos = self.pos + 1;
            Some(retval)
        } else {
            None
        }
    }
}

/// Iterate over the `(Kmer, Exts)` tuples of a sequence and it's extensions efficiently
pub struct KmerExtsIter<'a, K: Kmer, D>
where
    D: 'a,
{
    bases: &'a D,
    exts: Exts,
    kmer: K,
    pos: usize,
}

impl<'a, K: Kmer, D: Mer> Iterator for KmerExtsIter<'a, K, D> {
    type Item = (K, Exts);

    fn next(&mut self) -> Option<(K, Exts)> {
        if self.pos <= self.bases.len() {
            let next_base = if self.pos < self.bases.len() {
                self.bases.get(self.pos)
            } else {
                0u8
            };

            let cur_left = if self.pos == K::k() {
                self.exts
            } else {
                Exts::mk_left(self.bases.get(self.pos - K::k() - 1))
            };

            let cur_right = if self.pos < self.bases.len() {
                Exts::mk_right(next_base)
            } else {
                self.exts
            };

            let cur_exts = Exts::merge(cur_left, cur_right);

            let retval = self.kmer;
            self.kmer = self.kmer.extend_right(next_base);
            self.pos = self.pos + 1;
            Some((retval, cur_exts))
        } else {
            None
        }
    }
}
