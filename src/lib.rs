#![allow(dead_code)]

extern crate num;
extern crate extprim;
extern crate rand;
extern crate linked_hash_map;
#[macro_use]
extern crate serde_derive;
extern crate serde;
extern crate smallvec;
extern crate bit_set;
extern crate itertools;

use std::hash::Hash;
use std::fmt;
use num::PrimInt;
use num::FromPrimitive;
use extprim::u128::u128;

pub mod dna_string;
pub mod exts;
pub mod dir;
pub mod paths;
pub mod vmer;
pub mod msp;
pub mod filter;
mod fx;
mod test;

use dir::Dir;
use exts::Exts;


/// Convert a 2-bit representation of a base to a char
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
pub fn base_to_bits(c: u8) -> u8 {
    match c {
        b'A' => 0u8,
        b'C' => 1u8,
        b'G' => 2u8,
        b'T' => 3u8,
        _ => 0u8,
    }
}


/// Convert a 2-bit representation of a base to a char
pub fn bits_to_base(c: u8) -> char {
    match c {
        0u8 => 'A',
        1u8 => 'C',
        2u8 => 'G',
        3u8 => 'T',
        _ => 'X',
    }
}

pub trait IntHelp: PrimInt + FromPrimitive {
    fn reverse_by_twos(&self) -> Self;
}

impl IntHelp for u128 {
    #[inline]
    fn reverse_by_twos(&self) -> u128 {
        u128::from_parts(self.low64().reverse_by_twos(),
                         self.high64().reverse_by_twos())
    }
}

impl IntHelp for u64 {
    #[inline]
    fn reverse_by_twos(&self) -> u64 {
        // swap adjacent pairs
        let mut r = ((self & 0x3333333333333333u64) << 2) | ((self >> 2) & 0x3333333333333333u64);

        // swap nibbles
        r = ((r & 0x0F0F0F0F0F0F0F0Fu64) << 4) | ((r >> 4) & 0x0F0F0F0F0F0F0F0Fu64);

        // swap bytes
        r = ((r & 0x00FF00FF00FF00FFu64) << 8) | ((r >> 8) & 0x00FF00FF00FF00FFu64);

        // swap 2 bytes
        r = ((r & 0x0000FFFF0000FFFFu64) << 16) | ((r >> 16) & 0x0000FFFF0000FFFFu64);

        // swap 4 bytes
        r = ((r & 0x00000000FFFFFFFFu64) << 32) | ((r >> 32) & 0x00000000FFFFFFFFu64);

        r
    }
}

impl IntHelp for u32 {
    #[inline]
    fn reverse_by_twos(&self) -> u32 {
        // swap adjacent pairs
        let mut r = ((self & 0x33333333u32) << 2) | ((self >> 2) & 0x33333333u32);

        // swap nibbles
        r = ((r & 0x0F0F0F0Fu32) << 4) | ((r >> 4) & 0x0F0F0F0Fu32);

        // swap bytes
        r = ((r & 0x00FF00FFu32) << 8) | ((r >> 8) & 0x00FF00FFu32);

        // swap 2 bytes
        r = ((r & 0x0000FFFFu32) << 16) | ((r >> 16) & 0x0000FFFFu32);

        r
    }
}

impl IntHelp for u16 {
    #[inline]
    fn reverse_by_twos(&self) -> u16 {
        // swap adjacent pairs
        let mut r = ((self & 0x3333u16) << 2) | ((self >> 2) & 0x3333u16);

        // swap nibbles
        r = ((r & 0x0F0Fu16) << 4) | ((r >> 4) & 0x0F0Fu16);

        // swap bytes
        r = ((r & 0x00FFu16) << 8) | ((r >> 8) & 0x00FFu16);

        r
    }
}


/// The complement of a 2-bit encoded base
pub fn complement(base: u8) -> u8 {
    (!base) & 0x3u8
}


pub trait Mer: Sized + fmt::Debug {
    fn len(&self) -> usize;
    fn get(&self, pos: usize) -> u8;

    fn set_mut(&mut self, pos: usize, val: u8);
    fn set_slice_mut(&mut self, pos: usize, nbases: usize, bits: u64);

    fn rc(&self) -> Self;

    fn extend_left(&self, v: u8) -> Self;
    fn extend_right(&self, v: u8) -> Self;

    fn extend(&self, v: u8, dir: Dir) -> Self {
        match dir {
            Dir::Left => self.extend_left(v),
            Dir::Right => self.extend_right(v),
        }
    }

    fn get_extensions(&self, exts: Exts, dir: Dir) -> Vec<Self> {
        let ext_bases = exts.get(dir);
        ext_bases
            .iter()
            .map(|b| self.extend(b.clone(), dir))
            .collect()
    }
}

pub trait Kmer: Sized + Copy + Mer + PartialEq + PartialOrd + Eq + Ord + Hash {
    fn empty() -> Self;
    fn k() -> usize;

    fn min_rc_flip(&self) -> (Self, bool) {
        let rc = self.rc();
        if *self < rc {
            (self.clone(), false)
        } else {
            (rc, true)
        }
    }

    fn min_rc(&self) -> Self {
        let rc = self.rc();
        if *self < rc { self.clone() } else { rc }
    }

    /// Test if this Lmer is a rc-palindrome
    fn is_palindrome(&self) -> bool {
        self.len() % 2 == 0 && *self == self.rc()
    }

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

    fn to_string(&self) -> String {
        let mut s = String::new();
        for pos in 0..self.len() {
            s.push(bits_to_base(self.get(pos)))
        }
        s
    }

    /// Generate all kmers from string
    fn kmers_from_string(str: &[u8]) -> Vec<Self> {
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
}

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

/// A fixed-length Kmer sequence.
#[derive(Copy, Clone, PartialEq, PartialOrd, Eq, Ord, Hash, Serialize, Deserialize)]
pub struct IntKmer<T: PrimInt + FromPrimitive + IntHelp> {
    pub storage: T,
}

impl<T: PrimInt + FromPrimitive + Hash + IntHelp> IntKmer<T> {
    fn msk() -> T {
        T::one() << 1 | T::one()
    }

    fn to_byte(v: T) -> u8 {
        T::to_u8(&v).unwrap()
    }

    fn from_byte(v: u8) -> T {
        T::from_u8(v).unwrap()
    }

    fn from_u64(v: u64) -> T {
        T::from_u64(v).unwrap()
    }

    pub fn from_bytes(bytes: &Vec<u8>) -> Vec<Self> {
        let mut r = Vec::new();

        let mut k0 = Self::empty();

        for i in 0..Self::k() {
            k0.set_mut(i, bytes[i])
        }

        r.push(k0);

        for i in 1..(bytes.len() - Self::k() + 1) {
            k0 = k0.clone().extend_right(bytes[Self::k() + i - 1]);
            r.push(k0);
        }
        r
    }

    fn addr(&self, pos: usize) -> usize {
        let top_base = Self::k() - 1;
        let bitpos = (top_base - pos) * 2;
        bitpos
    }

    fn _k() -> usize {
        // 4 bases per byte
        std::mem::size_of::<T>() * 4
    }

    fn _bits() -> usize {
        std::mem::size_of::<T>() * 8
    }

    fn top_mask(n_bases: usize) -> T {
        if n_bases > 0 {
            // first pos bases
            let one = T::one();
            ((one << (n_bases * 2)) - one) << (Self::_bits() - n_bases * 2)
        } else {
            T::zero()
        }
    }

    fn bottom_mask(n_bases: usize) -> T {
        if n_bases > 0 {
            // first pos bases
            let one = T::one();
            ((one << (n_bases * 2)) - one)
        } else {
            T::zero()
        }
    }
}


impl<T: PrimInt + FromPrimitive + Hash + IntHelp> Mer for IntKmer<T> {
    fn len(&self) -> usize {
        Self::_k()
    }

    /// Get the letter at the given position.
    fn get(&self, pos: usize) -> u8 {
        let bit = self.addr(pos);
        Self::to_byte((self.storage >> bit & Self::msk()))
    }

    fn set_mut(&mut self, pos: usize, v: u8) {
        let bit = self.addr(pos);
        let mask = !(Self::msk() << bit);

        self.storage = (self.storage & mask) | (Self::from_byte(v) << bit);
    }

    /// Set a slice of bases in the kmer, using the packed representation in value.
    /// Sets n_bases, starting at pos. Bases must always be packed into the upper-most
    /// bits of the value.
    fn set_slice_mut(&mut self, pos: usize, n_bases: usize, value: u64) {
        debug_assert!(pos + n_bases <= Self::k());

        let v_shift = if Self::_bits() < 64 {
            value >> (64 - Self::_bits())
        } else {
            value
        };

        let v = if Self::_bits() > 64 {
            Self::from_u64(v_shift) << (Self::_bits() - 64)
        } else {
            Self::from_u64(v_shift)
        };

        let top_mask = Self::top_mask(pos);
        let bottom_mask = Self::bottom_mask(Self::k() - (pos + n_bases));
        let mask = top_mask | bottom_mask;

        let value_slide = v >> (2 * pos);

        self.storage = (self.storage & mask) | (value_slide & !mask);
    }

    /// Return the reverse complement of this kmer
    fn rc(&self) -> Self {
        // not bits to get complement, then reverse order
        let new = !self.storage.reverse_by_twos();

        // FIXME - deal with case when the kmer doesn't fill the bits
        //let up_shift = 2 * (64 - Self::k());
        //let down_shift = 64 - up_shift;
        //let u2 = new >> down_shift;

        IntKmer { storage: new }
    }

    /// Shift the base v into the left end of the kmer
    fn extend_left(&self, v: u8) -> Self {
        let new = self.storage >> 2 | (Self::from_byte(v) << (Self::k() - 1) * 2);
        IntKmer { storage: new }
    }

    fn extend_right(&self, v: u8) -> Self {
        let new = self.storage << 2;
        let mut kmer = IntKmer { storage: new };
        kmer.set_mut(Self::k() - 1, v);
        kmer
    }
}

impl<T: PrimInt + FromPrimitive + Hash + IntHelp> Kmer for IntKmer<T> {
    fn empty() -> Self {
        IntKmer { storage: T::zero() }
    }

    fn k() -> usize {
        Self::_k()
    }
}

impl<T: PrimInt + FromPrimitive + Hash + IntHelp> fmt::Debug for IntKmer<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s = String::new();
        for pos in 0..Self::k() {
            s.push(bits_to_base(self.get(pos)))
        }

        write!(f, "{}", s)
    }
}

/// Iterate over the Kmers of a sequence efficiently
pub struct KmerIter<'a, K: Kmer, D>
    where D: 'a
{
    bases: &'a D,
    kmer: K,
    pos: usize,
}

impl<'a, K: Kmer, D: Mer> Iterator for KmerIter<'a, K, D> {
    type Item = K;

    fn next(&mut self) -> Option<K> {
        if self.pos <= self.bases.len() {
            let retval = self.kmer;
            self.kmer = self.kmer.extend_right(self.bases.get(self.pos));
            self.pos = self.pos + 1;
            Some(retval)
        } else {
            None
        }
    }
}

/// Iterate over the Kmers of a sequence efficiently
pub struct KmerExtsIter<'a, K: Kmer, D>
    where D: 'a
{
    bases: &'a D,
    exts: Exts,
    kmer: K,
    pos: usize,
}

impl<'a, K: Kmer, D: Mer> Iterator for KmerExtsIter<'a, K, D> {
    type Item = (K, Exts);

    fn next(&mut self) -> Option<(K,Exts)> {
        if self.pos <= self.bases.len() {

            let next_base = self.bases.get(self.pos);

            let cur_left = 
                if self.pos == K::k() {
                    self.exts
                } else {
                    Exts::mk_left(self.bases.get(self.pos - K::k() - 1))
                };

            let cur_right = 
                if self.pos < self.bases.len() {
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



#[cfg(test)]
mod tests {
    use super::*;
    use rand::{self, Rng};
    use extprim::u128::u128;

    use vmer::Lmer;
    use vmer::Vmer;

    use KmerIter;

    fn check_kmer<T: Kmer>() {
        let K = T::k();

        let km = random_kmer::<T>();

        let rc = km.rc();
        let double_rc = rc.rc();
        assert!(km == double_rc);

        for i in 0..K {
            assert!(km.get(i) == (3 - rc.get(K - 1 - i)))
        }


        for i in 0..K {
            // Get and set
            let km2 = km.set(i, 0);
            assert!(km2.get(i) == 0);
        }

        let mut copy_kmer = T::empty();
        for i in 0..K {
            copy_kmer = copy_kmer.set(i, km.get(i));
        }
        assert!(km == copy_kmer);


        // Extend right
        let nb = random_base();
        let ext_r = km.extend_right(nb);
        assert!(ext_r.get(K - 1) == nb);
        assert!(km.get(1) == ext_r.get(0));

        // Extend Left
        let nb = random_base();
        let ext_l = km.extend_left(nb);
        assert!(ext_l.get(0) == nb);
        assert!(ext_l.get(1) == km.get(0));

        // Shift twice
        let l_base = random_base();
        let r_base = random_base();
        let ts = km.set(0, l_base).set(K - 1, r_base);

        let double_shift = ts.extend_left(0).extend_right(r_base);
        assert!(ts == double_shift)
    }


    fn check_vmer<V: Vmer<T>, T: Kmer>() {

        let vm = random_vmer::<V, T>();
        let l = vm.len();

        let rc = vm.rc();

        for i in 0..l {

            if vm.get(i) != (3 - rc.get(l - 1 - i)) {
                println!("km: {:?}, rc: {:?}", vm, rc);
            }

            assert!(vm.get(i) == (3 - rc.get(l - 1 - i)))
        }


        let double_rc = rc.rc();
        assert!(vm == double_rc);

        for i in 0..l {
            // Get and set
            let vm2 = vm.set(i, 0);
            assert!(vm2.get(i) == 0);
        }

        let mut copy_vmer = V::new(l);
        for i in 0..l {
            copy_vmer = copy_vmer.set(i, vm.get(i));
        }
        assert!(vm == copy_vmer);


        let kmers: Vec<T> = vm.iter_kmers().collect();
        assert_eq!(kmers.len(), vm.len() - T::k() + 1);

        for (pos, kmer) in kmers.iter().enumerate() {
            for i in 0..T::k() {
                assert_eq!(kmer.get(i), vm.get(i + pos))
            }
        }

        assert!(vm.last_kmer() == kmers[kmers.len() - 1]);
    }


    pub fn random_vmer<V: Vmer<T>, T: Kmer>() -> V {
        let mut r = rand::thread_rng();
        let len = r.gen_range(T::k(), V::max_len());

        let mut vmer = V::new(len);
        for pos in 0..len {
            let b = (r.next_u64() % 4) as u8;
            vmer = vmer.set(pos, b);
        }
        vmer
    }



    pub fn random_kmer<T: Kmer>() -> T {
        let mut r = rand::thread_rng();
        let mut kmer = T::empty();
        for pos in 0..T::k() {
            let b = (r.next_u64() % 4) as u8;
            kmer = kmer.set(pos, b);
        }
        kmer
    }

    pub fn random_base() -> u8 {
        let mut r = rand::thread_rng();
        (r.next_u64() % 4) as u8
    }


    #[test]
    fn test_lmer_3_kmer_64() {
        for _ in 0..10000 {
            check_vmer::<Lmer<IntKmer<u128>, [u64; 3]>, IntKmer<u128>>();
        }
    }

    #[test]
    fn test_lmer_3_kmer_32() {
        for _ in 0..10000 {
            check_vmer::<Lmer<IntKmer<u64>, [u64; 3]>, IntKmer<u64>>();
        }
    }

    #[test]
    fn test_lmer_2_kmer_32() {
        for _ in 0..10000 {
            check_vmer::<Lmer<IntKmer<u64>, [u64; 3]>, IntKmer<u64>>();
        }
    }

    #[test]
    fn test_lmer_1_kmer_16() {
        for _ in 0..10000 {
            check_vmer::<Lmer<IntKmer<u32>, [u64; 1]>, IntKmer<u32>>();
        }
    }

    #[test]
    fn test_kmer_64() {
        for _ in 0..10000 {
            check_kmer::<IntKmer<u128>>();
        }
    }

    #[test]
    fn test_kmer_32() {
        for _ in 0..10000 {
            check_kmer::<IntKmer<u64>>();
        }
    }

    #[test]
    fn test_kmer_16() {
        for _ in 0..10000 {
            check_kmer::<IntKmer<u32>>();
        }
    }

    #[test]
    fn test_kmer_8() {
        for _ in 0..10000 {
            check_kmer::<IntKmer<u16>>();
        }
    }
}
