// Copyright 2017 10x Genomics

//! Represent kmers with statically know length in compact integer types
//!
//! A kmer is a DNA sequence with statically-known length K. The sequence is stored the smallest possible integer type
//! Efficient methods for reverse complement and shifting new bases into the left or right of a sequence are provided.
//! Kmers implement `Eq` to test if two kmers represent the same string.
//! Kmers implement `Ord`, which corresponds to the natural lexicographic ordering.
//!
//! ```
//! use debruijn::*;
//! use debruijn::kmer::*;
//!
//! let k1 = Kmer16::from_ascii(b"ACGTACGTACGTACGT");
//!
//! // Reverse complement
//! let rc_k1 = k1.rc();
//!
//! // Double reverse complement
//! let k1_copy = rc_k1.rc();
//! assert_eq!(k1, k1_copy);
//!
//! // Push one base onto the left
//! assert_eq!(k1.extend_left(base_to_bits(b'T')), Kmer16::from_ascii(b"TACGTACGTACGTACG"));
//!
//! // Generate a set of kmers from a string, then sort
//! let mut all_kmers = Kmer16::kmers_from_ascii(b"TACGTACGTACGTACGTT");
//! all_kmers.sort();
//! assert_eq!(all_kmers,
//!     vec![
//!         Kmer16::from_ascii(b"ACGTACGTACGTACGT"),
//!         Kmer16::from_ascii(b"CGTACGTACGTACGTT"),
//!         Kmer16::from_ascii(b"TACGTACGTACGTACG")
//!     ]);


use std;
use std::hash::Hash;
use std::fmt;
use num::PrimInt;
use num::FromPrimitive;
use extprim::u128::u128;
use std::marker::PhantomData;

use Mer;
use Kmer;
use bits_to_base;


// Pre-defined kmer types

/// 64-base kmer, backed by a single u128
pub type Kmer64 = IntKmer<u128>;

/// 48-base kmer, backed by a single u128
pub type Kmer48 = VarIntKmer<u128, K48>;

/// 40-base kmer, backed by a single u128
pub type Kmer40 = VarIntKmer<u128, K40>;

/// 32-base kmer, backed by a single u64
pub type Kmer32 = IntKmer<u64>;

/// 24-base kmer, backed by a single u64
pub type Kmer24 = VarIntKmer<u64, K24>;

/// 16-base kmer, backed by a single u32
pub type Kmer16 = IntKmer<u32>;

/// 14-base kmer, backed by a single u32
pub type Kmer14 = VarIntKmer<u32, K14>;

/// 16-base kmer, backed by a single u16
pub type Kmer8 = IntKmer<u16>;



/// Trait for specialized integer operations used in DeBruijn Graph
pub trait IntHelp: PrimInt + FromPrimitive {
    /// Reverse the order of 2-bit units of the integer
    fn reverse_by_twos(&self) -> Self;
}

impl IntHelp for u128 {
    #[inline]
    fn reverse_by_twos(&self) -> u128 {
        u128::from_parts(
            self.low64().reverse_by_twos(),
            self.high64().reverse_by_twos(),
        )
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

impl IntHelp for u8 {
    #[inline]
    fn reverse_by_twos(&self) -> u8 {
        // swap adjacent pairs
        let mut r = ((self & 0x33u8) << 2) | ((self >> 2) & 0x33u8);

        // swap nibbles
        r = ((r & 0x0Fu8) << 4) | ((r >> 4) & 0x0Fu8);

        r
    }
}


/// A Kmer sequence with a statically know K. K will fill the underlying integer type.
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

    fn t_from_byte(v: u8) -> T {
        T::from_u8(v).unwrap()
    }

    fn t_from_u64(v: u64) -> T {
        T::from_u64(v).unwrap()
    }

    #[inline]
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

    pub fn top_mask(n_bases: usize) -> T {
        if n_bases > 0 {
            // first pos bases
            let one = T::one();
            ((one << (n_bases * 2)) - one) << (Self::_bits() - n_bases * 2)
        } else {
            T::zero()
        }
    }

    pub fn bottom_mask(n_bases: usize) -> T {
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
        Self::to_byte(self.storage >> bit & Self::msk())
    }

    fn set_mut(&mut self, pos: usize, v: u8) {
        let bit = self.addr(pos);
        let mask = !(Self::msk() << bit);

        self.storage = (self.storage & mask) | (Self::t_from_byte(v) << bit);
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
            Self::t_from_u64(v_shift) << (Self::_bits() - 64)
        } else {
            Self::t_from_u64(v_shift)
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

        // NOTE: IntKmer always fills the bits, so we don't need to shift here.
        IntKmer { storage: new }
    }

    /// Shift the base v into the left end of the kmer
    fn extend_left(&self, v: u8) -> Self {
        let new = self.storage >> 2 | (Self::t_from_byte(v) << (Self::k() - 1) * 2);
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

    fn from_u64(v: u64) -> IntKmer<T> {
        IntKmer { storage: Self::t_from_u64(v) }
    }

    fn to_u64(&self) -> u64 {
        T::to_u64(&self.storage).unwrap()
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

/// Helper trait for declaring the K value of a Kmer. Will be removed when const generics are available
pub trait KmerSize: Ord + Hash + Copy + fmt::Debug {
    #[allow(non_snake_case)]
    #[inline]
    fn K() -> usize;
}

/// A fixed-length Kmer sequence that may not fill the bits of T
///
/// side:             L           R
/// bases: 0   0   0  A  C  G  T  T
/// bits:             H  ........ L
/// bit :  14  12  10 8  6  4  2  0
///
/// sorting the integer will give a lexicographic sorting of the corresponding string.
///  kmers that don't fill `storage` are always aligned to the least signifcant bits
#[derive(Copy, Clone, PartialEq, PartialOrd, Eq, Ord, Hash, Serialize, Deserialize)]
pub struct VarIntKmer<T: PrimInt + FromPrimitive + IntHelp, KS: KmerSize> {
    pub storage: T,
    pub phantom: PhantomData<KS>,
}


impl<T: PrimInt + FromPrimitive + Hash + IntHelp, KS: KmerSize> Kmer for VarIntKmer<T, KS> {
    fn empty() -> Self {
        VarIntKmer {
            storage: T::zero(),
            phantom: PhantomData,
        }
    }

    #[inline]
    fn k() -> usize {
        Self::_k()
    }

    fn to_u64(&self) -> u64 {
        T::to_u64(&self.storage).unwrap()
    }

    fn from_u64(v: u64) -> Self {
        VarIntKmer { storage: Self::t_from_u64(v), phantom: PhantomData}
    }
}


impl<T: PrimInt + FromPrimitive + Hash + IntHelp, KS: KmerSize> VarIntKmer<T, KS> {
    fn msk() -> T {
        T::one() << 1 | T::one()
    }

    fn to_byte(v: T) -> u8 {
        T::to_u8(&v).unwrap()
    }

    fn t_from_byte(v: u8) -> T {
        T::from_u8(v).unwrap()
    }

    fn t_from_u64(v: u64) -> T {
        T::from_u64(v).unwrap()
    }

    fn addr(&self, pos: usize) -> usize {
        let top_base = Self::k() - 1;
        let bitpos = (top_base - pos) * 2;
        bitpos
    }

    // K of this kmer
    fn _k() -> usize {
        KS::K()
    }

    // Bits used by this kmer
    fn _bits() -> usize {
        Self::_k() * 2
    }

    fn _total_bits() -> usize {
        std::mem::size_of::<T>() * 8
    }

    // mask the unused bits at the top, plus the requested number of bases
    pub fn top_mask(n_bases: usize) -> T {
        let unused_bits = Self::_total_bits() - Self::_bits();
        let mask_bits = n_bases * 2 + unused_bits;

        if mask_bits > 0 {
            let one = T::one();
            ((one << mask_bits) - one) << (Self::_total_bits() - mask_bits)
        } else {
            T::zero()
        }
    }

    pub fn bottom_mask(n_bases: usize) -> T {
        if n_bases > 0 {
            // first pos bases
            let one = T::one();
            ((one << (n_bases * 2)) - one)
        } else {
            T::zero()
        }
    }
}


impl<T: PrimInt + FromPrimitive + Hash + IntHelp, KS: KmerSize> Mer for VarIntKmer<T, KS> {
    fn len(&self) -> usize {
        Self::_k()
    }

    /// Get the letter at the given position.
    fn get(&self, pos: usize) -> u8 {
        let bit = self.addr(pos);
        Self::to_byte(self.storage >> bit & Self::msk())
    }

    fn set_mut(&mut self, pos: usize, v: u8) {
        let bit = self.addr(pos);
        let mask = !(Self::msk() << bit);

        self.storage = (self.storage & mask) | (Self::t_from_byte(v) << bit);
    }

    /// Set a slice of bases in the kmer, using the packed representation in value.
    /// Sets n_bases, starting at pos. Incoming bases must always be packed into the upper-most
    /// bits of the value.
    fn set_slice_mut(&mut self, pos: usize, n_bases: usize, value: u64) {
        debug_assert!(pos + n_bases <= Self::k());

        // Move bases to the top of the smaller type
        let v_shift = if Self::_total_bits() < 64 {
            value >> (64 - Self::_total_bits())
        } else {
            value
        };

        // Move bases up to the top of this type
        let v = if Self::_total_bits() > 64 {
            Self::t_from_u64(v_shift) << (Self::_total_bits() - 64)
        } else {
            Self::t_from_u64(v_shift)
        };

        // Mask for where the bases will sit in the kmer
        let top_mask = Self::top_mask(pos);
        let bottom_mask = Self::bottom_mask(Self::k() - (pos + n_bases));
        let mask = top_mask | bottom_mask;

        // Move the base down to their home: position + unused high-order bits
        let shift = (2 * pos) + (Self::_total_bits() - Self::_bits());
        let value_slide = v >> shift;

        self.storage = (self.storage & mask) | (value_slide & !mask);
    }

    /// Return the reverse complement of this kmer
    fn rc(&self) -> Self {
        // not bits to get complement, then reverse order
        let mut new = !self.storage.reverse_by_twos();

        // deal with case when the kmer doesn't fill the bits
        if Self::k() < std::mem::size_of::<T>() * 4 {
            let up_shift = 2 * (std::mem::size_of::<T>() * 4 - Self::k());
            new = new >> up_shift;
        }

        VarIntKmer {
            storage: new,
            phantom: PhantomData,
        }
    }

    /// Shift the base v into the left end of the kmer
    fn extend_left(&self, v: u8) -> Self {
        let new = self.storage >> 2;
        let mut kmer = VarIntKmer {
            storage: new,
            phantom: PhantomData,
        };
        kmer.set_mut(0, v);
        kmer
    }

    fn extend_right(&self, v: u8) -> Self {
        let new = self.storage << 2 & !Self::top_mask(0);
        let mut kmer = VarIntKmer {
            storage: new,
            phantom: PhantomData,
        };
        kmer.set_mut(Self::k() - 1, v);
        kmer
    }
}

impl<T: PrimInt + FromPrimitive + Hash + IntHelp, KS: KmerSize> fmt::Debug for VarIntKmer<T, KS> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s = String::new();
        for pos in 0..Self::k() {
            s.push(bits_to_base(self.get(pos)))
        }

        write!(f, "{}", s)
    }
}

/// Marker struct for generating K=48 Kmers
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K48;

impl KmerSize for K48 {
    fn K() -> usize {
        48
    }
}

/// Marker trait for generating K=40 Kmers
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K40;

impl KmerSize for K40 {
    fn K() -> usize {
        40
    }
}

/// Marker trait for generating K=24 Kmers
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K24;

impl KmerSize for K24 {
    fn K() -> usize {
        24
    }
}

/// Marker trait for generating K=24 Kmers
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K14;

impl KmerSize for K14 {
    #[inline]
    fn K() -> usize {
        14
    }
}




#[cfg(test)]
mod tests {
    use super::*;
    use rand::{self, Rng};
    use extprim::u128::u128;

    use vmer::Lmer;

    use Vmer;
    use MerImmut;

    // Generate random kmers & test the methods for manipulating them
    fn check_kmer<T: Kmer>() {
        let K = T::k();

        // Random kmer
        let km = random_kmer::<T>();

        // Reverse complementing
        let rc = km.rc();
        let double_rc = rc.rc();
        assert!(km == double_rc);

        for i in 0..K {
            assert!(km.get(i) == (3 - rc.get(K - 1 - i)))
        }

        // Get and set bases
        for i in 0..K {
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
        assert!(ts == double_shift);

        if T::k() <= 32 {
            // Convert to and from u64.
            let u64_1 = km.to_u64();
            let km2   = T::from_u64(u64_1);
            let u64_2 = km2.to_u64();
            assert_eq!(km, km2);
            assert_eq!(u64_1, u64_2);
        }
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

            if vm.get_kmer(pos) != *kmer {
                println!(
                    "didn't get same kmer: i:{}, vm: {:?} kmer_iter: {:?}, kmer_get: {:?}",
                    pos,
                    vm,
                    kmer,
                    vm.get_kmer(pos)
                )
            }
            assert_eq!(vm.get_kmer(pos), *kmer);
        }

        assert!(vm.first_kmer() == kmers[0]);
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
    fn test_lmer_3_kmer_48() {
        for _ in 0..10000 {
            check_vmer::<Lmer<VarIntKmer<u128, K48>, [u64; 3]>, VarIntKmer<u128, K48>>();
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
    fn test_lmer_1_kmer_24() {
        for _ in 0..10000 {
            check_vmer::<Lmer<VarIntKmer<u64, K24>, [u64; 1]>, VarIntKmer<u64, K24>>();
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
    fn test_kmer_48() {
        for _ in 0..10000 {
            check_kmer::<VarIntKmer<u128, K48>>();
        }
    }

    #[test]
    fn test_kmer_32() {
        for _ in 0..10000 {
            check_kmer::<IntKmer<u64>>();
        }
    }

    #[test]
    fn test_kmer_24() {
        for _ in 0..10000 {
            check_kmer::<VarIntKmer<u64, K24>>();
        }
    }

    #[test]
    fn test_kmer_16() {
        for _ in 0..10000 {
            check_kmer::<IntKmer<u32>>();
        }
    }

    #[test]
    fn test_kmer_14() {
        for _ in 0..10000 {
            check_kmer::<VarIntKmer<u32, K14>>();
        }
    }

    #[test]
    fn test_kmer_8() {
        for _ in 0..10000 {
            check_kmer::<IntKmer<u16>>();
        }
    }
}
