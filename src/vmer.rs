use std::marker::PhantomData;
use std::cmp::{max, min};
use std::fmt;
use std::hash::Hash;

use Kmer;
use Mer;
use Dir;
use IntKmer;
use IntHelp;
use bits_to_base;

fn block_set(kmer: u64, pos: usize, val: u8) -> u64 {
    let offset = (31 - pos) * 2;
    let mask = !(3 << offset);

    (kmer & mask) | ((val as u64) << offset)
}

fn block_get(kmer: u64, pos: usize) -> u8 {
    let offset = (31 - pos) * 2;
    ((kmer >> offset) & 3) as u8
}

pub trait Vmer<T: Kmer>: Mer {
    fn new(len: usize) -> Self;
    fn max_len() -> usize;

    fn get_kmer(&self, pos: usize) -> T;
    fn first_kmer(&self) -> T;
    fn last_kmer(&self) -> T;

    /// Get the terminal kmer of the sequence, on the side of the sequence given by dir
    fn term_kmer(&self, dir: Dir) -> T {
        match dir {
            Dir::Left => self.first_kmer(),
            Dir::Right => self.last_kmer(),
        }
    }

    fn iter_kmers(&self) -> KmerIter<T, Self>;
}


/// Store a variable-length DNA sequence in a packed 2-bit encoding, up 92bp in length
/// The length of the sequence is stored in the lower 8 bits of storage
#[derive(Hash, Copy, Clone, PartialEq, PartialOrd, Eq, Ord)]
pub struct Lmer<T: Kmer, A: Array> {
    storage: A,
    phantom:  PhantomData<T>,
}


impl<T: Kmer, A: Array<Item=u64> + Copy + Eq + Ord + Hash> Mer for Lmer<T, A> {

    /// The length of the DNA string
    fn len(&self) -> usize {
        (self.storage.as_slice()[A::size() - 1] & 0xff) as usize
    }

        /// Return a new Lmer with position pos set to base val
    fn set(&self, pos: usize, val: u8) -> Lmer<T, A> {
        let block = pos / 32;
        let offset = pos % 32;

        let block_val = block_set(self.storage.as_slice()[block], offset, val);
        let mut new_bases = self.storage.clone();
        new_bases.as_mut_slice()[block] = block_val;
        Lmer { storage: new_bases, phantom: PhantomData }
    }

    /// Get the base at position pos
    fn get(&self, pos: usize) -> u8 {
        let block = pos / 32;
        let offset = pos % 32;
        block_get(self.storage.as_slice()[block], offset)
    }

    fn set_slice(&self, pos: usize, n_bases: usize, value: u64) -> Self {

        let mut storage = self.storage.clone();
        {
            let slc = storage.as_mut_slice();

            let b0 = pos / 32;
            let block_pos = pos % 32;
            let top_mask = IntKmer::<u64>::top_mask(block_pos);
            let mut bottom_mask = IntKmer::<u64>::bottom_mask(max(0, 32usize.saturating_sub(block_pos + n_bases)));
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

        Lmer { storage: storage, phantom: PhantomData }
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

            let v_rc = !v.reverse_by_twos() << (64 - n_bases*2);
            new_lmer = new_lmer.set_slice(self.len() - pos - n_bases, n_bases, v_rc);
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


impl<T: Kmer, A: Array<Item=u64> + Copy + Eq + Ord + Hash> Vmer<T> for Lmer<T, A> {
    fn max_len() -> usize {
        (A::size() * 64 - 8) / 2
    }

    /// Initialize an blank Lmer of length len.
    /// Will initially represent all A's.
    fn new(len: usize) -> Lmer<T,A> {
        let mut arr = A::new();
        {
            let slc = arr.as_mut_slice();

            // Write the length into the last 8 bits
            slc[A::size() - 1] = (len as u64) & 0xff;
        }
        Lmer { storage: arr, phantom: PhantomData }
    }

    /// Get the kmer starting at position pos
    fn get_kmer(&self, pos: usize) -> T {
        assert!(self.len() - pos >= T::k());
        let slc = self.storage.as_slice();

        // Which block has the first base
        let mut block = pos / 32;

        // Where we are in the kmer
        let mut kmer_pos = 0;

        // Where are in the block
        let mut block_pos = pos % 32;

        let mut kmer = T::empty();

        while kmer_pos < T::k() {
            // get relevent bases for current block
            let nb = min(T::k() - kmer_pos, 32 - block_pos);
            
            let val = slc[block] << (2*block_pos);
            kmer = kmer.set_slice(kmer_pos, nb, val);

            // move to next block, move ahead in kmer. 
            block += 1;
            kmer_pos += nb;
            // alway start a beginning of next block
            block_pos = 0;
        }

        kmer
    }

    fn first_kmer(&self) -> T { self.get_kmer(0) }
    fn last_kmer(&self) -> T { self.get_kmer(self.len() - T::k()) }

    /// Get the terminal kmer of the sequence, on the side of the sequence given by dir
    fn term_kmer(&self, dir: Dir) -> T {
        match dir {
            Dir::Left => self.first_kmer(),
            Dir::Right => self.last_kmer(),
        }
    }

    /// Efficiently iterate over the kmers in the sequence
    fn iter_kmers(&self) -> KmerIter<T,Self>
    {
        KmerIter {
            bases: self,
            kmer: self.first_kmer(),
            pos: T::k(),
        }
    }
}


/// Iterate over the Kmers of an Lmer efficiently
pub struct KmerIter<'a,T: Kmer, V: Vmer<T>> where T: 'a, V: 'a
{
    bases: &'a V,
    kmer: T,
    pos: usize
}

impl<'a,T: Kmer, V: Vmer<T>> Iterator for KmerIter<'a,T, V> {
    type Item=T;

    fn next(&mut self) -> Option<T>
    {
        if self.pos <= self.bases.len()
        {
            let retval = self.kmer;
            self.kmer = self.kmer.extend_right(self.bases.get(self.pos));
            self.pos = self.pos + 1;
            Some(retval)
        }
        else
        {
            None
        }
    }
}


impl<T: Kmer, A: Array<Item=u64> + Copy + Eq + Ord + Hash> fmt::Debug for Lmer<T, A> {
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
