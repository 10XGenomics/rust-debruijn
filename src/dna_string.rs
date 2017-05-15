// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! A fixed-width bit encoding implementation. This allows to store a sequence of values over
//! a reduced alphabet by packing them bit-encoded into a sequence of bytes.

use std::fmt;
use Kmer;
use bits_to_base;
use base_to_bits;
use bits_to_ascii;
use std::cmp::min;
use IntHelp;

use Mer;
use vmer::Vmer;

const BLOCK_BITS: usize = 64;
const WIDTH: usize = 2;

const MASK: u64 = 0x3;

/// A sequence of DnaStringoded values. This implementation is specialized to 2-bit alphabets
#[derive(Ord, PartialOrd, Clone, PartialEq, Eq, Hash)]
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

    fn set_slice_mut(&mut self, pos: usize, nbases: usize, bits: u64) {
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

    fn extend_left(&self, _: u8) -> Self {
        unimplemented!()
    }

    fn extend_right(&self, _: u8) -> Self {
        unimplemented!()
    }
}

impl<K> Vmer<K> for DnaString
    where K: Kmer
{
    fn new(len: usize) -> Self {
        Self::with_capacity(len)
    }

    fn max_len() -> usize {
        <usize>::max_value()
    }

    /// Get the kmer starting at position pos
    fn get_kmer(&self, pos: usize) -> K {
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
    /// Create a new instance with a given encoding width (e.g. width=2 for using two bits per value).
    pub fn new() -> DnaString {
        DnaString {
            storage: Vec::new(),
            len: 0,
        }
    }

    /// Create a new instance with a given capacity.
    pub fn with_capacity(n: usize) -> Self {
        DnaString {
            storage: Vec::with_capacity(n * WIDTH / 64 + 1),
            len: 0,
        }
    }

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

    pub fn from_bytes(bytes: &Vec<u8>) -> DnaString {
        let mut dna_string = DnaString {
            storage: Vec::new(),
            len: 0,
        };

        for b in bytes.iter() {
            dna_string.push(*b)
        }

        dna_string
    }

    pub fn to_dna_string(&self) -> String {
        let mut dna: String = String::new();
        for v in self.iter() {
            dna.push(bits_to_base(v));
        }
        dna
    }

    pub fn to_ascii_vec(&self) -> Vec<u8> {
        let mut res = Vec::new();
        for v in self.iter() {
            res.push(bits_to_ascii(v));
        }
        res
    }

    /// Append a value.
    pub fn push(&mut self, value: u8) {
        let (block, bit) = self.addr(self.len);
        if bit == 0 {
            self.storage.push(0);
        }
        self.set_by_addr(block, bit, value);
        self.len += 1;
    }

    /// Push values read from a byte array.
    ///
    /// # Arguments
    /// `bytes`: byte array to read values from
    /// `seq_length`: how many values to read from the byte array. Note that this
    /// is number of values not number of elements of the byte array.
    pub fn push_bytes(&mut self, bytes: &Vec<u8>, seq_length: usize) {
        assert!(seq_length <= bytes.len() * 8 / WIDTH,
                "Number of elements to push exceeds array length");

        for i in 0..seq_length {
            let byte_index = (i * WIDTH) / 8;
            let byte_slot = (i * WIDTH) % 8;

            let v = bytes[byte_index];
            let bits = (v >> byte_slot) & (MASK as u8);

            self.push(bits);
        }
    }

    /// Append `n` times the given value.
    pub fn push_values(&mut self, mut n: usize, value: u8) {
        {
            // fill the last block
            let (block, mut bit) = self.addr(self.len);
            if bit > 0 {
                // TODO use step_by once it has been stabilized: for bit in (bit..64).step_by(self.width) {
                while bit <= 64 {
                    self.set_by_addr(block, bit, value);
                    n -= 1;
                    bit = bit + WIDTH
                }
            }
        }

        // pack the value into a block
        let mut value_block = 0;
        {
            let mut v = value as u64;
            for _ in 0..(64 / WIDTH) {
                value_block |= v;
                v <<= WIDTH;
            }
        }

        // push as many value blocks as needed
        let i = self.len + n;
        let (block, bit) = self.addr(i);
        for _ in self.storage.len()..block {
            self.storage.push(value_block);
        }

        if bit > 0 {
            // add the remaining values to a final block
            self.storage.push(value_block >> (64 - bit));
        }

        self.len = i;
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

    pub fn width(&self) -> usize {
        WIDTH
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Get the prefix of a given size.
    pub fn pref(&self, k: usize) -> DnaString {
        assert!(k <= self.len, "Prefix size exceeds number of elements.");
        if k == 0 {
            return DnaString::new();
        }
        // Get position of k-th element (element with index k - 1)
        let (block, bitpos) = self.addr(k - 1);
        // We need to keep the first block-1 elements of the storage vector
        // plus the first bitpos positions of the block-th element.
        let mut new_storage: Vec<u64> = self.storage.iter().take(block + 1).cloned().collect();
        // zero-out all bits after bitpos
        if bitpos + WIDTH < 64 {
            new_storage[block] &= (1 << (bitpos + WIDTH)) - 1;
        }
        DnaString {
            storage: new_storage,
            len: k,
        }
    }

    pub fn suf(&self, k: usize) -> DnaString {
        assert!(k <= self.len, "Suffix size exceeds number of elements.");
        // There's a faster way but involves shifting bits down from "higher" blocks...
        // let (block, bitpos) = self.addr(self.len - k);
        let values: Vec<u8> = self.iter().skip(self.len - k).collect();
        let mut dna_string = DnaString::new();
        for v in values {
            dna_string.push(v);
        }
        dna_string
    }

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


#[derive(Eq, PartialEq, Clone)]
pub struct DnaStringSlice<'a> {
    pub dna_string: &'a DnaString,
    pub start: usize,
    pub length: usize,
}


impl<'a> Mer for DnaStringSlice<'a> {
    fn len(&self) -> usize {
        self.length
    }

    /// Get the value at position `i`.
    fn get(&self, i: usize) -> u8 {
        self.dna_string.get(i + self.start)
    }

    /// Set the value as position `i`.
    fn set_mut(&mut self, _: usize, _: u8) {
        unimplemented!()
        //debug_assert!(i < self.length);
        //self.dna_string.set_mut(i + self.start, value);
    }

    fn set_slice_mut(&mut self, _: usize, _: usize, _: u64) {
        unimplemented!();
    }

    fn rc(&self) -> DnaStringSlice<'a> {
        unimplemented!();
    }

    fn extend_left(&self, _: u8) -> Self {
        unimplemented!()
    }

    fn extend_right(&self, _: u8) -> Self {
        unimplemented!()
    }
}

impl<'a, K> Vmer<K> for DnaStringSlice<'a>
    where K: Kmer
{
    fn new(_: usize) -> Self {
        unimplemented!()
    }

    fn max_len() -> usize {
        <usize>::max_value()
    }

    /// Get the kmer starting at position pos
    fn get_kmer(&self, pos: usize) -> K {
        debug_assert!(self.length + pos + K::k() <= self.length);
        self.dna_string.get_kmer(self.start + pos)
    }
}



impl<'a> DnaStringSlice<'a> {

    pub fn is_palindrome(&self) -> bool {
        // FIXME
        return false;
    }

    pub fn ascii(&self) -> Vec<u8> {
        let mut v = Vec::new();
        for pos in self.start..(self.start + self.length) {
            v.push(bits_to_ascii(self.dna_string.get(pos)));
        }
        v
    }

    pub fn to_dna_string(&self) -> String {
        let mut dna: String = String::new();
        for pos in self.start..(self.start + self.length) {
            dna.push(bits_to_base(self.dna_string.get(pos)));
        }
        dna
    }

    pub fn to_string(&self) -> String {
        String::from_utf8(self.ascii()).unwrap()
    }

    pub fn to_owned(&self) -> DnaString {
        let mut be = DnaString::with_capacity(self.length);
        for pos in self.start..(self.start + self.length) {
            be.push(self.dna_string.get(pos));
        }

        be
    }
}


impl<'a> fmt::Debug for DnaStringSlice<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s = String::new();
        for pos in self.start..(self.start + self.length) {
            s.push(bits_to_base(self.dna_string.get(pos)))
        }

        write!(f, "{}", s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use IntKmer;
    use KmerIter;

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
    fn test_push_values() {
        let mut dna_string = DnaString::new();
        dna_string.push_values(32, 0);
        assert_eq!(dna_string.storage, [0]);
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

        let dna_cp = dna_string.to_dna_string();
        assert_eq!(dna, dna_cp);
    }

    #[test]
    fn test_pref() {
        let in_values: Vec<u8> = vec![2, 20];
        let mut dna_string = DnaString::new();
        dna_string.push_bytes(&in_values, 8);
        // Contents should be 00010100 00000010

        let pref_dna_string = dna_string.pref(0);
        assert_eq!(pref_dna_string.storage.len(), 0);

        let pref_dna_string = dna_string.pref(8);
        assert_eq!(pref_dna_string, dna_string);

        let pref_dna_string = dna_string.pref(4);
        assert_eq!(pref_dna_string.storage, [2]);

        let pref_dna_string = dna_string.pref(6);
        assert_eq!(pref_dna_string.storage, [1026]);
        let values: Vec<u8> = pref_dna_string.iter().collect();
        assert_eq!(values, [2, 0, 0, 0, 0, 1]);

        dna_string.push_bytes(&in_values, 8);
        dna_string.push_bytes(&in_values, 8);

        let pref_dna_string = dna_string.pref(17);
        let values: Vec<u8> = pref_dna_string.iter().collect();
        assert_eq!(values, [2, 0, 0, 0, 0, 1, 1, 0, 2, 0, 0, 0, 0, 1, 1, 0, 2]);
    }

    #[test]
    fn test_suf() {
        let in_values: Vec<u8> = vec![2, 20];
        let mut dna_string = DnaString::new();
        dna_string.push_bytes(&in_values, 8);
        // Contents should be 00010100 00000010

        let suf_dna_string = dna_string.suf(0);
        assert_eq!(suf_dna_string.storage.len(), 0);

        let suf_dna_string = dna_string.suf(8);
        assert_eq!(suf_dna_string, dna_string);

        let suf_dna_string = dna_string.suf(4);
        assert_eq!(suf_dna_string.storage, [20]);

        // 000101000000 64+256
        let suf_dna_string = dna_string.suf(6);
        assert_eq!(suf_dna_string.storage, [320]);
        let values: Vec<u8> = suf_dna_string.iter().collect();
        assert_eq!(values, [0, 0, 0, 1, 1, 0]);

        dna_string.push_bytes(&in_values, 8);
        dna_string.push_bytes(&in_values, 8);

        let suf_dna_string = dna_string.suf(17);
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
