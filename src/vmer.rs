use std::marker::PhantomData;
use std::cmp::{max, min};
use std::fmt;
use std::hash::Hash;

use Kmer;
use KmerIter;
use KmerExtsIter;
use Mer;
use Dir;
use Exts;
use IntKmer;
use IntHelp;
use bits_to_base;
use base_to_bits;

fn block_set(kmer: u64, pos: usize, val: u8) -> u64 {
    let offset = (31 - pos) * 2;
    let mask = !(3 << offset);

    (kmer & mask) | ((val as u64) << offset)
}

fn block_get(kmer: u64, pos: usize) -> u8 {
    let offset = (31 - pos) * 2;
    ((kmer >> offset) & 3) as u8
}

pub trait Vmer<K: Kmer>: Mer + PartialEq + Eq + Clone {
    fn new(len: usize) -> Self;
    fn max_len() -> usize;

    fn get_kmer(&self, pos: usize) -> K;

    fn first_kmer(&self) -> K {
        self.get_kmer(0)
    }

    fn last_kmer(&self) -> K {
        self.get_kmer(self.len() - K::k())
    }

    /// Get the terminal kmer of the sequence, on the side of the sequence given by dir
    fn term_kmer(&self, dir: Dir) -> K {
        match dir {
            Dir::Left => self.first_kmer(),
            Dir::Right => self.last_kmer(),
        }
    }

    fn iter_kmers(&self) -> KmerIter<K, Self> {
        KmerIter {
            bases: self,
            kmer: self.first_kmer(),
            pos: K::k(),
        }
    }

    fn iter_kmer_exts(&self, seq_exts: Exts) -> KmerExtsIter<K, Self> {
        KmerExtsIter {
            bases: self,
            exts: seq_exts,
            kmer: self.first_kmer(),
            pos: K::k(),
        }
    }

    fn from_slice(seq: &[u8]) -> Self {
        let mut vmer = Self::new(seq.len());
        for i in 0 .. seq.len() {
            vmer.set_mut(i, seq[i]);
        }

        vmer
    }
}


/// Store a variable-length DNA sequence in a packed 2-bit encoding, up 92bp in length
/// The length of the sequence is stored in the lower 8 bits of storage
#[derive(Hash, Copy, Clone, PartialEq, PartialOrd, Eq, Ord, Serialize, Deserialize)]
pub struct Lmer<K: Kmer, A: Array> {
    storage: A,
    phantom: PhantomData<K>,
}


impl<K: Kmer, A: Array<Item = u64> + Copy + Eq + Ord + Hash> Mer for Lmer<K, A> {
    /// The length of the DNA string
    fn len(&self) -> usize {
        (self.storage.as_slice()[A::size() - 1] & 0xff) as usize
    }

    /// Get the base at position pos
    fn get(&self, pos: usize) -> u8 {
        let block = pos / 32;
        let offset = pos % 32;
        block_get(self.storage.as_slice()[block], offset)
    }

    /// Return a new Lmer with position pos set to base val
    fn set_mut(&mut self, pos: usize, val: u8) {
        let block = pos / 32;
        let offset = pos % 32;

        let block_val = block_set(self.storage.as_slice()[block], offset, val);
        self.storage.as_mut_slice()[block] = block_val;
    }

    fn set_slice_mut(&mut self, pos: usize, n_bases: usize, value: u64) {
        let slc = self.storage.as_mut_slice();
        let b0 = pos / 32;
        let block_pos = pos % 32;
        let top_mask = IntKmer::<u64>::top_mask(block_pos);
        let mut bottom_mask =
            IntKmer::<u64>::bottom_mask(max(0, 32usize.saturating_sub(block_pos + n_bases)));
        if b0 == A::size() - 1 {
            bottom_mask = bottom_mask | 0xFF
        }
        let mask = top_mask | bottom_mask;

        let nb0 = 32 - block_pos;
        let value_top = value >> (block_pos * 2);
        let v0 = (slc[b0] & mask) | (value_top & !mask);
        slc[b0] = v0;

        if n_bases > nb0 {
            let b1 = b0 + 1;
            let nb1 = n_bases - nb0;
            let bottom_mask = IntKmer::<u64>::bottom_mask(32 - nb1);
            let value_bottom = value << (nb0 * 2);
            let v1 = (slc[b1] & bottom_mask) | (value_bottom & !bottom_mask);
            slc[b1] = v1
        }
    }


    fn rc(&self) -> Self {
        let slc = self.storage.as_slice();

        let mut new_lmer = Self::new(self.len());
        let mut block = 0;
        let mut pos = 0;

        while pos < self.len() {
            let n_bases = min(32, self.len() - pos);

            let mut v = slc[block];
            // Mask the length field packed into the last block
            if block == A::size() - 1 {
                v = v & !0xFF
            }

            let v_rc = !v.reverse_by_twos() << (64 - n_bases * 2);
            new_lmer.set_slice_mut(self.len() - pos - n_bases, n_bases, v_rc);
            block += 1;
            pos += n_bases;
        }

        new_lmer
    }

    fn extend_left(&self, _: u8) -> Self {
        unimplemented!();
    }

    fn extend_right(&self, _: u8) -> Self {
        unimplemented!();
    }
}


impl<K: Kmer, A: Array<Item = u64> + Copy + Eq + Ord + Hash> Vmer<K> for Lmer<K, A> {
    fn max_len() -> usize {
        (A::size() * 64 - 8) / 2
    }

    /// Initialize an blank Lmer of length len.
    /// Will initially represent all A's.
    fn new(len: usize) -> Lmer<K, A> {
        let mut arr = A::new();
        {
            let slc = arr.as_mut_slice();

            // Write the length into the last 8 bits
            slc[A::size() - 1] = (len as u64) & 0xff;
        }
        Lmer {
            storage: arr,
            phantom: PhantomData,
        }
    }

    /// Get the kmer starting at position pos
    fn get_kmer(&self, pos: usize) -> K {
        assert!(self.len() - pos >= K::k());
        let slc = self.storage.as_slice();

        // Which block has the first base
        let mut block = pos / 32;

        // Where we are in the kmer
        let mut kmer_pos = 0;

        // Where are in the block
        let mut block_pos = pos % 32;

        let mut kmer = K::empty();

        while kmer_pos < K::k() {
            // get relevent bases for current block
            let nb = min(K::k() - kmer_pos, 32 - block_pos);

            let val = slc[block] << (2 * block_pos);
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

impl<K: Kmer, A: Array<Item = u64> + Copy + Eq + Ord + Hash> fmt::Debug for Lmer<K, A> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s = String::new();
        for pos in 0..self.len() {
            s.push(bits_to_base(self.get(pos)))
        }

        write!(f, "{}", s)
    }
}

/// Types that can be used as the backing store for a SmallVec
pub trait Array {
    type Item;
    fn new() -> Self;
    fn size() -> usize;

    fn as_slice(&self) -> &[Self::Item];
    fn as_mut_slice(&mut self) -> &mut [Self::Item];
}

macro_rules! impl_array(
    ($($size:expr),+) => {
        $(
            impl<T: Default + Copy + Eq + Ord> Array for [T; $size] {
                type Item = T;
                fn new() -> [T; $size] { [T::default(); $size] }
                fn size() -> usize { $size }
                fn as_slice(&self) -> &[T] { self }
                fn as_mut_slice(&mut self) -> &mut [T] { self }

            }
        )+
    }
);

impl_array!(1, 2, 3, 4, 5, 6);
