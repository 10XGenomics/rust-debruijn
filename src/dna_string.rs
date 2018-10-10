// Copyright 2014 Johannes Köster and 10x Genomics
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! A 2-bit encoding of arbitrary length DNA sequences.
//!
//! Store arbitrary-length DNA strings in a packed 2-bit encoding. Individual base values are encoded
//! as the integers 0,1,2,3 corresponding to A,C,G,T.
//!
//! # Example
//! ```
//! use debruijn::Kmer;
//! use debruijn::dna_string::*;
//! use debruijn::kmer::Kmer16;
//! use debruijn::Vmer;
//!
//! // Construct a new DNA string
//! let dna_string1 = DnaString::from_dna_string("ACAGCAGCAGCACGTATGACAGATAGTGACAGCAGTTTGTGACCGCAAGAGCAGTAATATGATG");
//!
//! // Get an immutable view into the sequence
//! let slice1 = dna_string1.slice(10, 40);
//!
//! // Get a kmer from the DNA string
//! let first_kmer: Kmer16 = slice1.get_kmer(0);
//! assert_eq!(first_kmer, Kmer16::from_ascii(b"CACGTATGACAGATAG"))


use std::fmt;
use std::borrow::Borrow;

use Kmer;
use bits_to_base;
use base_to_bits;
use bits_to_ascii;
use dna_only_base_to_bits;
use std::cmp::min;
use kmer::IntHelp;
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;

use Mer;
use MerIter;
use Vmer;

const BLOCK_BITS: usize = 64;
const WIDTH: usize = 2;

const MASK: u64 = 0x3;

/// A container for sequence of DNA bases.
#[derive(Ord, PartialOrd, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct DnaString {
    storage: Vec<u64>,
    len: usize,
}

impl Mer for DnaString {
    fn len(&self) -> usize {
        self.len
    }

    /// Get the value at position `i`.
    fn get(&self, i: usize) -> u8 {
        let (block, bit) = self.addr(i);
        self.get_by_addr(block, bit)
    }

    /// Set the value as position `i`.
    fn set_mut(&mut self, i: usize, value: u8) {
        let (block, bit) = self.addr(i);
        self.set_by_addr(block, bit, value);
    }

    fn set_slice_mut(&mut self, _: usize, _: usize, _: u64) {
        unimplemented!()
    }

    fn rc(&self) -> DnaString {
        let mut dna_string = DnaString::new();
        for i in (0..self.len()).rev() {
            let v = 3 - self.get(i);
            dna_string.push(v);
        }
        dna_string
    }
}

impl Vmer for DnaString
{
    fn new(len: usize) -> Self {
        Self::empty(len)
    }

    fn max_len() -> usize {
        <usize>::max_value()
    }

    /// Get the kmer starting at position pos
    fn get_kmer<K: Kmer>(&self, pos: usize) -> K {
        assert!(self.len() - pos >= K::k());

        // Which block has the first base
        let (mut block, _) = self.addr(pos);

        // Where we are in the kmer
        let mut kmer_pos = 0;

        // Where are in the block
        let mut block_pos = pos % 32;

        let mut kmer = K::empty();

        while kmer_pos < K::k() {
            // get relevent bases for current block
            let nb = min(K::k() - kmer_pos, 32 - block_pos);

            let v = self.storage[block].reverse_by_twos();
            let val = v << (2 * block_pos);
            kmer.set_slice_mut(kmer_pos, nb, val);

            // move to next block, move ahead in kmer.
            block += 1;
            kmer_pos += nb;
            // alway start a beginning of next block
            block_pos = 0;
        }

        kmer
    }
}


impl DnaString {
    /// Create an empty DNA string
    pub fn new() -> DnaString {
        DnaString {
            storage: Vec::new(),
            len: 0,
        }
    }

    /// Length of the sequence
    pub fn len(&self) -> usize {
        self.len
    }

    /// Create a new instance with a given capacity.
    pub fn empty(n: usize) -> Self {
        let blocks = (n * WIDTH >> 6) + (if n * WIDTH & 0x3F > 0 { 1 } else { 0 });
        let storage = vec![0; blocks];

        DnaString {
            storage: storage,
            len: n,
        }
    }

    /// Create a DnaString corresponding to an ACGT-encoded str.
    pub fn from_dna_string(dna: &str) -> DnaString {
        let mut dna_string = DnaString {
            storage: Vec::new(),
            len: 0,
        };

        for c in dna.chars() {
            dna_string.push(base_to_bits(c as u8));
        }

        dna_string
    }

    /// Create a DnaString corresponding to an ACGT-encoded str.
    pub fn from_dna_only_string(dna: &str) -> Vec<DnaString> {
        let mut dna_vector: Vec<DnaString> = Vec::new();
        let mut dna_string = DnaString::new();

        for c in dna.chars() {
            match dna_only_base_to_bits(c as u8) {
                Some(bit) => {
                    dna_string.push(bit);
                }
                None => {
                    if dna_string.len() > 0 {
                        dna_vector.push(dna_string);
                        dna_string = DnaString::new();
                    }
                }
            }
        }
        if dna_string.len() > 0 {
            dna_vector.push(dna_string);
        }

        dna_vector
    }


    /// Create a DnaString from an ASCII ACGT-encoded byte slice.
    /// Non ACGT positions will be converted to 'A'
    pub fn from_acgt_bytes(bytes: &[u8]) -> DnaString {
        let mut dna_string = DnaString {
            storage: Vec::new(),
            len: 0,
        };

        for b in bytes.iter() {
            dna_string.push(base_to_bits(*b));
        }

        dna_string
    }

    /// Create a DnaString from an ACGT-encoded byte slice,
    /// Non ACGT positions will be converted to repeatable random base determined
    /// by a hash of the read name and the position within the string.
    pub fn from_acgt_bytes_hashn(bytes: &[u8], read_name: &[u8]) -> DnaString {

        let mut hasher = DefaultHasher::new();
        read_name.hash(&mut hasher);
        
        let mut dna_string = DnaString {
            storage: Vec::new(),
            len: 0,
        };

        for (pos, c) in bytes.iter().enumerate() {

            let v = match c {
                b'A' => 0u8,
                b'C' => 1u8,
                b'G' => 2u8,
                b'T' => 3u8,
                _ => {
                    let mut hasher_clone = hasher.clone();
                    pos.hash(&mut hasher_clone);
                    (hasher_clone.finish() % 4) as u8
                    },
            };

            dna_string.push(v);
        }

        dna_string
    }

    /// Create a DnaString from a 0-4 encoded byte slice
    pub fn from_bytes(bytes: &[u8]) -> DnaString {
        let mut dna_string = DnaString {
            storage: Vec::new(),
            len: 0,
        };

        for b in bytes.iter() {
            dna_string.push(*b)
        }

        dna_string
    }

    /// Convert sequence to a String
    pub fn to_string(&self) -> String {
        let mut dna: String = String::new();
        for v in self.iter() {
            dna.push(bits_to_base(v));
        }
        dna
    }

    /// Convert sequence to a Vector of 0-4 encoded bytes
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();
        bytes.extend(self.iter());
        bytes
    }

    /// Convert sequence to a Vector of ascii-encoded bytes
    pub fn to_ascii_vec(&self) -> Vec<u8> {
        let mut res = Vec::new();
        for v in self.iter() {
            res.push(bits_to_ascii(v));
        }
        res
    }

    /// Append a 0-4 encoded base.
    pub fn push(&mut self, value: u8) {
        let (block, bit) = self.addr(self.len);
        if bit == 0 && block >= self.storage.len() {
            self.storage.push(0);
        }
        self.set_by_addr(block, bit, value);
        self.len += 1;
    }

    /// Push 0-4 encoded bases from a byte array.
    ///
    /// # Arguments
    /// `bytes`: byte array to read values from
    /// `seq_length`: how many values to read from the byte array. Note that this
    /// is number of values not number of elements of the byte array.
    pub fn push_bytes(&mut self, bytes: &Vec<u8>, seq_length: usize) {
        assert!(
            seq_length <= bytes.len() * 8 / WIDTH,
            "Number of elements to push exceeds array length"
        );

        for i in 0..seq_length {
            let byte_index = (i * WIDTH) / 8;
            let byte_slot = (i * WIDTH) % 8;

            let v = bytes[byte_index];
            let bits = (v >> byte_slot) & (MASK as u8);

            self.push(bits);
        }
    }

    /// Iterate over stored values (values will be unpacked into bytes).
    pub fn iter(&self) -> DnaStringIter {
        DnaStringIter {
            dna_string: self,
            i: 0,
        }
    }

    /// Clear the sequence.
    pub fn clear(&mut self) {
        self.storage.clear();
        self.len = 0;
    }

    fn get_by_addr(&self, block: usize, bit: usize) -> u8 {
        ((self.storage[block] >> bit) & MASK) as u8
    }

    fn set_by_addr(&mut self, block: usize, bit: usize, value: u8) {
        let mask = MASK << bit;
        self.storage[block] |= mask;
        self.storage[block] ^= mask;
        self.storage[block] |= (value as u64 & MASK) << bit;
    }

    fn addr(&self, i: usize) -> (usize, usize) {
        let k = i * WIDTH;
        (k / BLOCK_BITS, k % BLOCK_BITS)
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Get the length `k` prefix of the DnaString
    pub fn prefix(&self, k: usize) -> DnaStringSlice {
        assert!(k <= self.len, "Prefix size exceeds number of elements.");
        DnaStringSlice {
            dna_string: self,
            start: 0,
            length: k,
            is_rc: false,
        }
    }

    /// Get the length `k` suffix of the DnaString
    pub fn suffix(&self, k: usize) -> DnaStringSlice {
        assert!(k <= self.len, "Suffix size exceeds number of elements.");

        DnaStringSlice {
            dna_string: self,
            start: self.len() - k,
            length: k,
            is_rc: false,
        }
    }

    /// Get slice containing the interval [`start`, `end`) of `self`
    pub fn slice(&self, start: usize, end: usize) -> DnaStringSlice {
        assert!(start <= self.len, "coordinate exceeds number of elements.");
        assert!(end <= self.len, "coordinate exceeds number of elements.");

        DnaStringSlice {
            dna_string: self,
            start: start,
            length: end - start,
            is_rc: false,
        }
    }

    /// Create a fresh DnaString containing the reverse of `self`
    pub fn reverse(&self) -> DnaString {
        let values: Vec<u8> = self.iter().collect();
        let mut dna_string = DnaString::new();
        for v in values.iter().rev() {
            dna_string.push(*v);
        }
        dna_string
    }

    // pub fn complement(&self) -> DnaString {
    //    assert!(self.width == 2, "Complement only supported for 2bit encodings.");
    //    let values: Vec<u32> = Vec::with_capacity(self.len());
    //    for i, v in self.storage.iter() {
    //        values[i] = v;
    //    }
    //    values[values.len() - 1] =
    // }
}


impl fmt::Debug for DnaString {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s = String::new();
        for pos in 0..self.len() {
            s.push(bits_to_base(self.get(pos)))
        }

        write!(f, "{}", s)
    }
}

/// Iterator over values of a DnaStringoded sequence (values will be unpacked into bytes).
pub struct DnaStringIter<'a> {
    dna_string: &'a DnaString,
    i: usize,
}

impl<'a> Iterator for DnaStringIter<'a> {
    type Item = u8;

    fn next(&mut self) -> Option<u8> {
        if self.i < self.dna_string.len() {
            let value = self.dna_string.get(self.i);
            self.i += 1;
            Some(value)
        } else {
            None
        }

    }
}

impl<'a> IntoIterator for &'a DnaString {
    type Item = u8;
    type IntoIter = DnaStringIter<'a>;

    fn into_iter(self) -> DnaStringIter<'a> {
        self.iter()
    }
}


/// An immutable slice into a DnaString
#[derive(Clone)]
pub struct DnaStringSlice<'a> {
    pub dna_string: &'a DnaString,
    pub start: usize,
    pub length: usize,
    pub is_rc: bool,
}

impl<'a> PartialEq for DnaStringSlice<'a> {
    fn eq( &self, other: &DnaStringSlice ) -> bool {
        let n = self.length;

        if other.length != n { return false; }

        println!( "\nn = {}", n );
        println!( "self.is_rc = {}", self.is_rc );
        println!( "other.is_rc = {}", other.is_rc );
        for i in 0..n {
            println!( "{} vs {}", self.get(i), other.get(i) );
        }

        if self.is_rc == other.is_rc {
            for i in 0..n {
                if self.get( self.start + i ) != other.get( other.start + i ) { 

                    // XXX:
                    println!( "eq returning false for {} and {}, whose actual\
                        equality is {}",
                        self.to_string(), other.to_string(),
                        self.to_string() == other.to_string() );

                    return false; 
                }
            }
        }
        else {
            for i in 0..n {
                if self.get( self.start + i ) 
                        != ::complement( other.get( other.start + n - i - 1 ) ) { 

                    // XXX:
                    println!( "eq-rc returning false for {} and {}, whose actual\
                        equality is {}",
                        self.to_string(), other.to_string(),
                        self.to_string() == other.to_string() );

                    return false; 
                }
            }
        }

        // XXX:
        println!( "eq returning true for {}.{} and {}.{}, \
            whose actual equality is {}",
            self.to_string(), self.is_rc, other.to_string(), other.is_rc,
            self.to_string() == other.to_string() );

        true
    }
}
impl<'a> Eq for DnaStringSlice<'a> { }

impl<'a> Mer for DnaStringSlice<'a> {
    fn len(&self) -> usize {
        self.length
    }

    /// Get the value at position `i`.
    fn get(&self, i: usize) -> u8 {
        if !self.is_rc {
            self.dna_string.get(i + self.start)
        } else {
            ::complement(self.dna_string.get(self.start + self.length - 1 - i))
        }
    }

    /// Set the value as position `i`.
    fn set_mut(&mut self, _: usize, _: u8) {
        unimplemented!()
        //debug_assert!(i < self.length);
        //self.dna_string.set(i + self.start, value);
    }

    fn set_slice_mut(&mut self, _: usize, _: usize, _: u64) {
        unimplemented!();
    }

    fn rc(&self) -> DnaStringSlice<'a> {
        DnaStringSlice {
            dna_string: self.dna_string,
            start: self.start,
            length: self.length,
            is_rc: !self.is_rc,
        }
    }
}

impl<'a> Vmer for DnaStringSlice<'a>
{
    fn new(_: usize) -> Self {
        unimplemented!()
    }

    fn max_len() -> usize {
        <usize>::max_value()
    }

    /// Get the kmer starting at position pos
    fn get_kmer<K: Kmer>(&self, pos: usize) -> K {
        debug_assert!(pos + K::k() <= self.length);
        self.dna_string.get_kmer(self.start + pos)
    }
}



impl<'a> DnaStringSlice<'a> {
    pub fn is_palindrome(&self) -> bool {
        unimplemented!();
    }

    pub fn bytes(&self) -> Vec<u8> {
        let mut v = Vec::new();
        for pos in 0..self.length {
            v.push(self.get(pos));
        }
        v
    }

    pub fn ascii(&self) -> Vec<u8> {
        let mut v = Vec::new();
        for pos in 0..self.length {
            v.push(bits_to_ascii(self.get(pos)));
        }
        v
    }

    pub fn to_dna_string(&self) -> String {
        let mut dna: String = String::new();
        for pos in 0..self.length {
            dna.push(bits_to_base(self.get(pos)));
        }
        dna
    }

    pub fn to_string(&self) -> String {
        String::from_utf8(self.ascii()).unwrap()
    }

    pub fn to_owned(&self) -> DnaString {
        // FIXME make this faster
        let mut be = DnaString::empty(self.length);
        for pos in 0..self.length {
            be.set_mut(pos, self.get(pos));
        }

        be
    }
        /// Get slice containing the interval [`start`, `end`) of `self`
    pub fn slice(&self, start: usize, end: usize) -> DnaStringSlice {
        assert!(start <= self.length, "coordinate exceeds number of elements.");
        assert!(end <= self.length, "coordinate exceeds number of elements.");

        DnaStringSlice {
            dna_string: self.dna_string,
            start: self.start + start,
            length: end - start,
            is_rc: false,
        }
    }

}


impl<'a> fmt::Debug for DnaStringSlice<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s = String::new();
        if self.length < 256 {
            for pos in self.start..(self.start + self.length) {
                s.push(bits_to_base(self.dna_string.get(pos)))
            }
            write!(f, "{}", s)
        } else {
            write!(f, "start: {}, len: {}, is_rc: {}", self.start, self.length, self.is_rc)
        }
    }
}

impl<'a> IntoIterator for &'a DnaStringSlice<'a> {
    type Item = u8;
    type IntoIter = MerIter<'a, DnaStringSlice<'a>>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}


/// Container for many distinct sequences, concatenated into a single DnaString.  Each
/// sequence is accessible by index as a DnaStringSlice.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PackedDnaStringSet {
    pub sequence: DnaString,
    pub start: Vec<usize>,
    pub length: Vec<u32>,
}

impl<'a> PackedDnaStringSet {
    /// Create an empty `PackedDnaStringSet`
    pub fn new() -> Self {
        PackedDnaStringSet {
            sequence: DnaString::new(),
            start: Vec::new(),
            length: Vec::new(),
        }
    }

    /// Get a `DnaStringSlice` containing `i`th sequence in the set
    pub fn get(&'a self, i: usize) -> DnaStringSlice<'a> {
        DnaStringSlice {
            dna_string: &self.sequence,
            start: self.start[i],
            length: self.length[i] as usize,
            is_rc: false,
        }
    }

    /// Get a `DnaStringSlice` containing `i`th sequence in the set
    pub fn slice(&'a self, i: usize, start: usize, end: usize) -> DnaStringSlice<'a> {
        assert!(start <= self.length[i] as usize);
        assert!(end <= self.length[i] as usize);

        DnaStringSlice {
            dna_string: &self.sequence,
            start: self.start[i] + start,
            length: end - start,
            is_rc: false,
        }
    }


    /// Number of sequences in the set
    pub fn len(&self) -> usize {
        self.start.len()
    }

    pub fn add<'b, R: Borrow<u8>, S: IntoIterator<Item = R>>(&mut self, sequence: S) {
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


#[cfg(test)]
mod tests {
    use super::*;
    use kmer::IntKmer;

    #[test]
    fn test_dna_string() {
        let mut dna_string = DnaString::new();
        dna_string.push(0);
        dna_string.push(2);
        dna_string.push(1);
        let mut values: Vec<u8> = dna_string.iter().collect();
        assert_eq!(values, [0, 2, 1]);
        dna_string.set_mut(1, 3);
        values = dna_string.iter().collect();
        assert_eq!(values, [0, 3, 1]);
    }

    #[test]
    fn test_push_bytes() {
        let in_values: Vec<u8> = vec![2, 20];

        let mut dna_string = DnaString::new();
        dna_string.push_bytes(&in_values, 8);
        // Contents should be 00010100 00000010
        let values: Vec<u8> = dna_string.iter().collect();
        assert_eq!(values, [2, 0, 0, 0, 0, 1, 1, 0]);
        assert_eq!(dna_string.storage, [5122]);

        let mut dna_string = DnaString::new();
        dna_string.push_bytes(&in_values, 2);
        // Contents should be 00000010
        let values: Vec<u8> = dna_string.iter().collect();
        assert_eq!(values, [2, 0]);
        assert_eq!(dna_string.storage, [2]);
    }

    #[test]
    fn test_from_dna_string() {
        let dna = "ACGTACGT";
        let dna_string = DnaString::from_dna_string(dna);
        let values: Vec<u8> = dna_string.iter().collect();

        assert_eq!(dna_string.len, 8);
        assert_eq!(values, [0, 1, 2, 3, 0, 1, 2, 3]);

        let dna_cp = dna_string.to_string();
        assert_eq!(dna, dna_cp);
    }

    #[test]
    fn test_prefix() {
        let in_values: Vec<u8> = vec![2, 20];
        let mut dna_string = DnaString::new();
        dna_string.push_bytes(&in_values, 8);
        // Contents should be 00010100 00000010

        let pref_dna_string = dna_string.prefix(0).to_owned();
        assert_eq!(pref_dna_string.len(), 0);

        let pref_dna_string = dna_string.prefix(8).to_owned();
        assert_eq!(pref_dna_string, dna_string);

        let pref_dna_string = dna_string.prefix(4).to_owned();
        assert_eq!(pref_dna_string.storage, [2]);

        let pref_dna_string = dna_string.prefix(6).to_owned();
        assert_eq!(pref_dna_string.storage, [1026]);
        let values: Vec<u8> = pref_dna_string.iter().collect();
        assert_eq!(values, [2, 0, 0, 0, 0, 1]);

        dna_string.push_bytes(&in_values, 8);
        dna_string.push_bytes(&in_values, 8);

        let pref_dna_string = dna_string.prefix(17).to_owned();
        let values: Vec<u8> = pref_dna_string.iter().collect();
        assert_eq!(values, [2, 0, 0, 0, 0, 1, 1, 0, 2, 0, 0, 0, 0, 1, 1, 0, 2]);
    }

    #[test]
    fn test_suffix() {
        let in_values: Vec<u8> = vec![2, 20];
        let mut dna_string = DnaString::new();
        dna_string.push_bytes(&in_values, 8);
        // Contents should be 00010100 00000010

        let suf_dna_string = dna_string.suffix(0).to_owned();
        assert_eq!(suf_dna_string.len(), 0);

        let suf_dna_string = dna_string.suffix(8).to_owned();
        assert_eq!(suf_dna_string, dna_string);

        let suf_dna_string = dna_string.suffix(4).to_owned();
        assert_eq!(suf_dna_string.storage, [20]);

        // 000101000000 64+256
        let suf_dna_string = dna_string.suffix(6).to_owned();
        assert_eq!(suf_dna_string.storage, [320]);
        let values: Vec<u8> = suf_dna_string.iter().collect();
        assert_eq!(values, [0, 0, 0, 1, 1, 0]);

        dna_string.push_bytes(&in_values, 8);
        dna_string.push_bytes(&in_values, 8);

        let suf_dna_string = dna_string.suffix(17).to_owned();
        let values: Vec<u8> = suf_dna_string.iter().collect();
        assert_eq!(values, [0, 2, 0, 0, 0, 0, 1, 1, 0, 2, 0, 0, 0, 0, 1, 1, 0]);
    }

    #[test]
    fn test_reverse() {
        let in_values: Vec<u8> = vec![2, 20];
        let mut dna_string = DnaString::new();
        let rev_dna_string = dna_string.reverse();
        assert_eq!(dna_string, rev_dna_string);

        dna_string.push_bytes(&in_values, 8);
        // Contents should be 00010100 00000010

        let rev_dna_string = dna_string.reverse();
        let values: Vec<u8> = rev_dna_string.iter().collect();
        assert_eq!(values, [0, 1, 1, 0, 0, 0, 0, 2]);
    }

    #[test]
    fn test_kmers() {
        let dna = "TGCATTAGAAAACTCCTTGCCTGTCAGCCCGACAGGTAGAAACTCATTAATCCACACATTGA".to_string() +
            "CTCTATTTCAGGTAAATATGACGTCAACTCCTGCATGTTGAAGGCAGTGAGTGGCTGAAACAGCATCAAGGCGTGAAGGC";
        let dna_string = DnaString::from_dna_string(&dna);

        let kmers: Vec<IntKmer<u64>> = dna_string.iter_kmers().collect();
        kmer_test::<IntKmer<u64>>(&kmers, &dna, &dna_string);
    }

    #[test]
    fn test_kmers_too_short() {
        let dna = "TGCATTAGAA".to_string();
        let dna_string = DnaString::from_dna_string(&dna);

        let kmers: Vec<IntKmer<u64>> = dna_string.iter_kmers().collect();
        assert_eq!(kmers, vec![]);
    }

    fn kmer_test<K: Kmer>(kmers: &Vec<K>, dna: &String, dna_string: &DnaString) {
        for i in 0..(dna.len() - K::k() + 1) {
            assert_eq!(kmers[i].to_string(), &dna[i..(i + K::k())]);
        }

        let last_kmer: K = dna_string.last_kmer();
        assert_eq!(last_kmer.to_string(), &dna[(dna.len() - K::k())..]);

        for (idx, &k) in kmers.iter().enumerate() {
            assert_eq!(k, dna_string.get_kmer(idx));
        }
    }
}
