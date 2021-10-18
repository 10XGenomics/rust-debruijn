//! Kmers represented by an array of integers.

use num_traits::AsPrimitive;
use num_traits::FromPrimitive;
use num_traits::PrimInt;
use std;
use std::fmt;
use std::fmt::Debug;
use std::hash::Hash;

use crate::bits_to_base;
use crate::kmer::IntHelp;
use crate::Kmer;
use crate::Mer;

/// The [`ArrayKmer`] struct implements a [`Kmer`] backed by an
/// array of integers. The type of the integers, the number of integers in the array and value of `K` are
/// specified with generic and const-generic parameters. The type of the integer is specified as `T`, the `K` value is specified as `KS`, and the number of
/// elements of the backing array is `B`.
///
/// This implementation of [`Kmer`] will generally be slower
/// than [`crate::kmer::IntKmer`] or [`crate::kmer::VarIntKmer`], but can be smaller by having less padding, and can go up to
/// K > 64, which is the current limit for [`crate::kmer::IntKmer`].
///
/// Implementation note: any padding required must always occupy the highest-order bits of the first element of the storage array.
/// The padding bits _must_ be maintained as 0 by all manipulations of the kmer.
/// The bits are ordered such that the natural ordering of the storage array corresponds to the lexicographic ordering of the kmers.
#[derive(Copy, Clone, PartialEq, PartialOrd, Eq, Ord, Hash)]
pub struct ArrayKmer<T: PrimInt + FromPrimitive + IntHelp, const KS: usize, const B: usize> {
    pub storage: [T; B],
}

impl<
        T: 'static + PrimInt + FromPrimitive + Hash + IntHelp + Debug,
        const KS: usize,
        const B: usize,
    > Kmer for ArrayKmer<T, KS, B>
where
    u64: AsPrimitive<T>,
{
    fn empty() -> Self {
        ArrayKmer {
            storage: [T::zero(); B],
        }
    }

    #[inline]
    fn k() -> usize {
        KS
    }

    fn to_u64(&self) -> u64 {
        let mut r = 0u64;
        for i in 0..B {
            let v = self.storage[B - i - 1].to_u64().unwrap();
            r = r | v << (i * Self::slot_bits());
        }

        r
    }

    fn from_u64(v: u64) -> Self {
        let mut new = Self::empty();

        // shift v to the upper bits of the u64, which is what copy_bases_u64 expects
        let vshift = v << (64 - KS * 2);

        new.copy_bases_u64(vshift, KS, 0);
        new
    }

    /// Shift the base v into the left end of the kmer
    fn extend_left(&self, v: u8) -> Self {
        let mut new: [T; B] = [T::zero(); B];

        new[0] = self.storage[0] >> 2;
        for i in 1..B {
            new[i] = self.storage[i] >> 2 | self.storage[i - 1] << (Self::slot_bits() - 2);
        }

        let mut kmer = ArrayKmer { storage: new };

        kmer.set_mut(0, v);
        kmer
    }

    fn extend_right(&self, v: u8) -> Self {
        let mut new = [T::zero(); B];

        for i in 0..B - 1 {
            new[i] = self.storage[i] << 2 | self.storage[i + 1] >> (Self::slot_bits() - 2);

            if i == 0 {
                let lower_mask = Self::mask_below((Self::slot_bases() - Self::padding_bases()) * 2);
                new[i] = new[i] & lower_mask;
            }
        }
        new[B - 1] = self.storage[B - 1] << 2;

        let mut kmer = ArrayKmer { storage: new };
        kmer.set_mut(Self::k() - 1, v);
        kmer
    }

    fn hamming_dist(&self, other: Self) -> u32 {
        let mut total_diffs = 0;
        for i in 0..B {
            let bit_diffs = self.storage[i] ^ other.storage[i];
            let two_bit_diffs = (bit_diffs | bit_diffs >> 1) & IntHelp::lower_of_two();
            total_diffs += two_bit_diffs.count_ones();
        }

        total_diffs
    }
}

impl<
        T: 'static + PrimInt + FromPrimitive + Hash + IntHelp + Debug,
        const KS: usize,
        const B: usize,
    > ArrayKmer<T, KS, B>
where
    u64: AsPrimitive<T>,
{
    #[inline(always)]
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
        // move upper bits of
        let shift = std::cmp::max(64 - Self::slot_bits(), 0);
        let b = v >> shift;
        b.as_()
    }

    #[inline(always)]
    fn slot_bits() -> usize {
        std::mem::size_of::<T>() * 8
    }

    #[inline(always)]
    fn slot_bases() -> usize {
        Self::slot_bits() / 2
    }

    #[inline(always)]
    fn padding_bases() -> usize {
        B * Self::slot_bases() - KS
    }

    /// return (slot, pos_in_slot, bitpos_in_slot) for base posistion `pos`
    #[inline(always)]
    fn addr(pos: usize) -> (usize, usize, usize) {
        let padding = Self::padding_bases();
        let padded_pos = pos + padding;

        let slot = padded_pos / Self::slot_bases();
        let pos_in_slot = padded_pos % Self::slot_bases();

        let top_base = Self::slot_bases() - 1;
        let bitpos = (top_base - pos_in_slot) * 2;

        let r = (slot, pos_in_slot, bitpos);
        r
    }

    /// Copy upper-most `nbases` bases from `src` into self, starting a base position `dest`
    fn copy_bases_u64(&mut self, mut src: u64, mut nbases: usize, mut dest: usize) {
        let (mut slot, mut slot_base, mut bit_pos) = Self::addr(dest);

        while nbases > 0 {
            let cp_bits = std::cmp::min(nbases * 2, bit_pos + 2);

            let top_mask = Self::mask_above(bit_pos + 1);
            let bottom_mask = Self::mask_below(bit_pos + 2 - cp_bits);

            let old_val = self.storage[slot];
            let new_bits = Self::t_from_u64(src >> (slot_base * 2));

            let new_val = (old_val & top_mask)
                | (old_val & bottom_mask)
                | (new_bits & !top_mask & !bottom_mask);

            self.storage[slot] = new_val;

            src = src << cp_bits;
            nbases = nbases - cp_bits / 2;
            dest = dest + cp_bits / 2;
            slot = slot + 1;
            slot_base = 0;
            bit_pos = Self::slot_bits() - 2;
        }
    }

    /// Copy upper-most `nbases` bases from `src` into self, starting a base position `dest`
    fn copy_bases(&mut self, mut src: T, mut nbases: usize, mut dest: usize) {
        let (mut slot, mut slot_base, mut bit_pos) = Self::addr(dest);

        loop {
            let cp_bits = std::cmp::min(nbases * 2, (Self::slot_bases() - slot_base) * 2);

            let top_mask = Self::mask_above(bit_pos + 1);
            let bottom_mask = Self::mask_below(bit_pos + 2 - cp_bits);

            let old_val = self.storage[slot];

            let new_bits = src >> (slot_base * 2);

            let new_val = (old_val & top_mask)
                | (old_val & bottom_mask)
                | (new_bits & !top_mask & !bottom_mask);

            self.storage[slot] = new_val;

            nbases = nbases - cp_bits / 2;

            if nbases == 0 {
                break;
            }

            src = src << cp_bits;
            dest = dest + cp_bits / 2;
            slot = slot + 1;
            slot_base = 0;
            bit_pos = Self::slot_bits() - 2;
        }
    }

    #[inline(always)]
    fn mask_above(bit: usize) -> T {
        if bit + 1 >= Self::slot_bits() {
            return T::zero();
        }

        let one = T::one();
        let start_bit = bit + 1;
        let mask_bits = Self::slot_bits() - start_bit;
        let res = ((one << mask_bits) - one) << start_bit;

        for i in 0..Self::slot_bits() {
            let res_bit = (res >> i) & T::one();
            if i > bit {
                debug_assert_eq!(res_bit, T::one());
            } else {
                debug_assert_eq!(res_bit, T::zero());
            }
        }

        res
    }

    #[inline(always)]
    fn mask_below(bit: usize) -> T {
        // special case to mask everything
        if bit == Self::slot_bits() {
            return !T::zero();
        }

        let one = T::one();
        let res = (one << bit) - one;

        for i in 0..Self::slot_bits() {
            let res_bit = (res >> i) & T::one();
            if i < bit {
                debug_assert_eq!(res_bit, T::one());
            } else {
                debug_assert_eq!(res_bit, T::zero());
            }
        }

        res
    }
}

impl<
        T: 'static + PrimInt + FromPrimitive + Hash + IntHelp + Debug,
        const KS: usize,
        const B: usize,
    > Mer for ArrayKmer<T, KS, B>
where
    u64: AsPrimitive<T>,
{
    #[inline(always)]
    fn len(&self) -> usize {
        KS
    }

    /// Get the letter at the given position.
    fn get(&self, pos: usize) -> u8 {
        let (slot, _, bit) = Self::addr(pos);
        Self::to_byte(self.storage[slot] >> bit & Self::msk())
    }

    fn set_mut(&mut self, pos: usize, v: u8) {
        let (slot, _, bit) = Self::addr(pos);
        let mask = !(Self::msk() << bit);

        self.storage[slot] = (self.storage[slot] & mask) | (Self::t_from_byte(v) << bit);
    }

    /// Set a slice of bases in the kmer, using the packed representation in value.
    /// Sets n_bases, starting at pos. Incoming bases must always be packed into the upper-most
    /// bits of the value.
    #[inline(always)]
    fn set_slice_mut(&mut self, pos: usize, n_bases: usize, value: u64) {
        debug_assert!(pos + n_bases <= Self::k());
        self.copy_bases_u64(value, n_bases, pos);
    }

    /// Return the reverse complement of this kmer
    fn rc(&self) -> Self {
        let mut new = ArrayKmer::empty();
        let mut cur_pos = 0;

        while cur_pos < KS {
            let (slot, _, bitpos) = Self::addr(cur_pos);

            let s = self.storage[slot];
            let s_rc = !s.reverse_by_twos();
            let nbases = (bitpos / 2) + 1;

            let dest = KS - cur_pos - nbases;
            new.copy_bases(s_rc, nbases, dest);
            cur_pos += nbases;
        }

        new
    }

    fn at_count(&self) -> u32 {
        // A's and T's have upper_bit ^ lower_bit == 0
        // count how many of these are present

        let mut count = 0;
        let slot0 = self.storage[0];

        let slot0_bases = KS - Self::slot_bases() * (B - 1);
        let mix_base_bits = !((slot0 >> 1) ^ slot0);
        let at_bases = mix_base_bits & Self::mask_below(slot0_bases * 2) & IntHelp::lower_of_two();
        count += at_bases.count_ones();

        for i in 1..B {
            let mix_base_bits = !((self.storage[i] >> 1) ^ self.storage[i]);
            let at_bases = mix_base_bits & IntHelp::lower_of_two();
            count += at_bases.count_ones();
        }

        count
    }

    fn gc_count(&self) -> u32 {
        (KS as u32) - self.at_count()
    }
}

impl<
        T: 'static + PrimInt + FromPrimitive + Hash + IntHelp + Debug,
        const KS: usize,
        const B: usize,
    > fmt::Debug for ArrayKmer<T, KS, B>
where
    u64: AsPrimitive<T>,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();
        for pos in 0..Self::k() {
            s.push(bits_to_base(self.get(pos)))
        }

        write!(f, "{}", s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer;

    #[test]
    fn simple_arr_test() {
        type K = ArrayKmer<u16, 26, 4>;

        let km = K::from_ascii(b"CCCCCCCCCCCCCCCCCCCCCCCCCC");
        println!("{:?}", km);

        let k1 = km.extend_left(0);
        println!("{:?}", k1);

        let k1 = km.extend_right(0);
        println!("{:?}", k1);

        let km = K::from_ascii(b"ACGTACGTACGTACGTACGTACGTAC");
        println!("{:?}", km);

        let k1 = km.extend_left(0);
        println!("{:?}", k1);

        let k1 = km.extend_right(0);
        println!("{:?}", k1);

        let km = K::from_ascii(b"CCCCCCCCCCCCCCCCCCCCCCCCCC");
        println!("{:?}", km);

        let mut k2 = km.clone();
        k2.copy_bases_u64(0xF << 60, 2, 0);
        assert_eq!(k2, K::from_ascii(b"TTCCCCCCCCCCCCCCCCCCCCCCCC"));

        let mut k2 = km.clone();
        k2.copy_bases_u64(0xF << 60, 2, 1);
        assert_eq!(k2, K::from_ascii(b"CTTCCCCCCCCCCCCCCCCCCCCCCC"));

        let mut k2 = km.clone();
        k2.copy_bases_u64(0xF << 60, 2, 2);
        assert_eq!(k2, K::from_ascii(b"CCTTCCCCCCCCCCCCCCCCCCCCCC"));

        let mut k2 = km.clone();
        k2.copy_bases_u64(0xF << 60, 2, 3);
        assert_eq!(k2, K::from_ascii(b"CCCTTCCCCCCCCCCCCCCCCCCCCC"));

        let mut k2 = km.clone();
        k2.copy_bases_u64(0xF << 60, 2, 7);
        assert_eq!(k2, K::from_ascii(b"CCCCCCCTTCCCCCCCCCCCCCCCCC"));

        let mut k2 = km.clone();
        k2.copy_bases_u64(0xFF << 56, 4, 7);
        assert_eq!(k2, K::from_ascii(b"CCCCCCCTTTTCCCCCCCCCCCCCCC"));

        let km = K::from_ascii(b"CCCCCCCCCCCCCCCCCCCCCCCCCC");
        assert_eq!(km.rc(), K::from_ascii(b"GGGGGGGGGGGGGGGGGGGGGGGGGG"));
    }

    #[test]
    fn test_arr_kmer15() {
        kmer::tests::kmer_test_suite::<ArrayKmer<u16, 15, 2>>();
    }

    #[test]
    fn test_arr_kmer26() {
        kmer::tests::kmer_test_suite::<ArrayKmer<u16, 26, 4>>();
    }

    #[test]
    fn test_arr_kmer23() {
        kmer::tests::kmer_test_suite::<ArrayKmer<u16, 23, 3>>();
    }

    #[test]
    fn test_arr_kmer48() {
        kmer::tests::kmer_test_suite::<ArrayKmer<u32, 48, 3>>();
    }

    #[test]
    fn test_arr_kmer47() {
        kmer::tests::kmer_test_suite::<ArrayKmer<u32, 47, 3>>();
    }

    #[test]
    fn test_arr_kmer65_u16() {
        kmer::tests::kmer_test_suite::<ArrayKmer<u16, 65, 9>>();
    }

    #[test]
    fn test_arr_kmer68_u16() {
        kmer::tests::kmer_test_suite::<ArrayKmer<u16, 68, 9>>();
    }

    #[test]
    fn test_arr_kmer65_u32() {
        kmer::tests::kmer_test_suite::<ArrayKmer<u32, 65, 5>>();
    }
}
