// Copyright 2017 10x Genomics

//! Methods for minimum substring partitioning of a DNA string
//!
//! simple_scan method is based on:
//! Li, Yang. "MSPKmerCounter: a fast and memory efficient approach for k-mer counting." arXiv preprint arXiv:1505.06550 (2015).

use crate::DnaSlice;
use crate::Exts;
use crate::Kmer;
use crate::Vmer;
use std::cmp::min;
use std::cmp::Ordering;
use std::iter::Iterator;
use std::ops::Range;

#[derive(Debug)]
pub struct MspInterval {
    bucket: u16,
    start: u32,
    len: u16,
}

impl MspInterval {
    pub fn new(bucket: u16, start: u32, len: u16) -> MspInterval {
        MspInterval {
            bucket: bucket,
            start: start,
            len: len,
        }
    }

    pub fn start(&self) -> usize {
        self.start as usize
    }

    pub fn len(&self) -> usize {
        self.len as usize
    }

    pub fn end(&self) -> usize {
        (self.len as u32 + self.start) as usize
    }

    pub fn range(&self) -> Range<usize> {
        self.start()..self.start() + self.len()
    }

    pub fn bucket(&self) -> u16 {
        self.bucket
    }
}

/// Determine MSP substrings of seq, for given k and p.
/// Returns a vector of tuples indicating the substrings, and the pmer values:
/// (p-mer value, min p-mer position, start position, end position)
/// permutation is a permutation of the lexicographically-sorted set of all pmers.
/// A permutation of pmers sorted by their inverse frequency in the dataset will give the
/// most even bucketing of MSPs over pmers.
#[deprecated(note = "Please use the `Scanner` type instead ")]
pub fn simple_scan<V: Vmer, P: Kmer>(
    k: usize,
    seq: &V,
    permutation: &[usize],
    rc: bool,
) -> Vec<MspInterval> {
    // Can't partition strings shorter than k
    assert!(seq.len() >= k);
    assert!(P::k() <= 8);
    assert!(seq.len() < 1 << 32);

    let score = |pi: &P| {
        if rc {
            min(
                permutation[pi.to_u64() as usize],
                permutation[pi.rc().to_u64() as usize],
            )
        } else {
            permutation[pi.to_u64() as usize]
        }
    };

    let scanner = Scanner::new(seq, score, k);
    let res = scanner.scan();

    res.into_iter()
        .map(|slc| MspInterval {
            bucket: slc.bucket() as u16,
            start: slc.start,
            len: slc.len,
        })
        .collect()
}

/// Represents a sequence interval composed of
/// successive k-mers that share a
/// common minizer p-mer.
#[derive(Debug)]
pub struct MspIntervalP<P> {
    /// The minimizing p-mer in this interval
    pub minimizer: P,
    /// The start of the sequence interval
    pub start: u32,
    /// The length of the sequence interval
    pub len: u16,
    /// The position of the minimizer
    pub minimizer_pos: u32,
}

impl<P: Kmer> MspIntervalP<P> {
    /// The shard 'bucket' identifer
    /// to put this sequence in.
    /// This is the reverse-complement canonicalized
    /// value of the p-mer.
    pub fn bucket(&self) -> u64 {
        self.minimizer.min_rc().to_u64()
    }
}

#[derive(PartialEq, Eq, Clone, Copy, Debug)]
struct MinPos<P> {
    val: usize,
    pos: usize,
    kmer: P,
}

impl<P: Eq> Ord for MinPos<P> {
    fn cmp(&self, other: &Self) -> Ordering {
        let val_cmp = self.val.cmp(&other.val);
        if val_cmp != Ordering::Equal {
            return val_cmp;
        }

        let pos_cmp = self.pos.cmp(&other.pos);
        match pos_cmp {
            Ordering::Equal => Ordering::Equal,
            Ordering::Less => Ordering::Greater,
            Ordering::Greater => Ordering::Less,
        }
    }
}

impl<P: Eq> PartialOrd for MinPos<P> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let val_cmp = self.val.cmp(&other.val);
        if val_cmp != Ordering::Equal {
            return Some(val_cmp);
        }

        let pos_cmp = self.pos.cmp(&other.pos);
        Some(match pos_cmp {
            Ordering::Equal => Ordering::Equal,
            Ordering::Less => Ordering::Greater,
            Ordering::Greater => Ordering::Less,
        })
    }
}

/// Determine MSP substrings of a sequence, for given k and p.
/// The `scan()` method Returns a vector of tuples indicating the substrings,
/// and the p-mer values as a set of `MspIntervalP<P>` values. A user-supplied
/// score function is used to rank p-mers for the purposes of finding the minimizer.
/// A permutation is a permutation of the lexicographically-sorted set of all pmers.
/// A permutation of pmers sorted by their inverse frequency in the dataset will give the
/// most even bucketing of MSPs over pmers.
pub struct Scanner<'a, V, F, P> {
    seq: &'a V,
    score: F,
    k: usize,
    _mp: MinPos<P>,
}

impl<'a, V, F, P> Scanner<'a, V, F, P>
where
    V: Vmer,
    P: Kmer,
    F: Fn(&P) -> usize,
{
    /// Build a scanner for the `sequence`, using minimizer function `score_func`,
    /// and kmer length `k`. The p-mer length is set by the type `P`.
    pub fn new(sequence: &'a V, score_func: F, k: usize) -> Scanner<'a, V, F, P> {
        Scanner {
            seq: sequence,
            score: score_func,
            k,
            _mp: MinPos {
                val: 0,
                pos: 0,
                kmer: P::empty(),
            },
        }
    }

    fn mp(&self, pos: usize) -> MinPos<P> {
        let kmer = self.seq.get_kmer::<P>(pos);
        let val = (self.score)(&kmer);
        MinPos { pos, val, kmer }
    }

    fn incr(&self, mp: &MinPos<P>) -> MinPos<P> {
        let pos = mp.pos + 1;
        let kmer = mp.kmer.extend_right(self.seq.get(pos + P::k() - 1));
        let val = (self.score)(&kmer);
        MinPos { pos, val, kmer }
    }

    pub fn scan(&self) -> Vec<MspIntervalP<P>> {
        // Can't partition strings shorter than k
        assert!(self.seq.len() >= self.k);
        assert!(self.seq.len() < 1 << 32);

        let seq = self.seq;
        let m = seq.len();

        let k = self.k;
        let p = P::k();

        let find_min = |start, stop| {
            let mut min_pos = self.mp(start);
            let mut current = min_pos.clone();

            while current.pos < stop {
                current = self.incr(&current);
                min_pos = min(min_pos, current);
            }

            min_pos
        };

        let mut min_positions = Vec::with_capacity(16);

        let mut min_pos = find_min(0, k - p);
        let mut end_pos = self.mp(k - p);

        min_positions.push((0, min_pos));

        for i in 1..(m - k + 1) {
            // end_pos always corresponds to i + k - p
            end_pos = self.incr(&end_pos);

            if i > min_pos.pos {
                min_pos = find_min(i, i + k - p);
                min_positions.push((i, min_pos));
            } else {
                if end_pos.val < min_pos.val {
                    min_pos = end_pos;
                    min_positions.push((i, min_pos));
                }
            }
        }

        let mut slices = Vec::with_capacity(min_positions.len());

        // Generate the slices of the final string
        for p in 0..min_positions.len() - 1 {
            let (start_pos, min_pos) = min_positions[p];
            let (next_pos, _) = min_positions[p + 1];

            let interval = MspIntervalP {
                minimizer: min_pos.kmer,
                minimizer_pos: min_pos.pos as u32,
                start: start_pos as u32,
                len: (next_pos + k - 1 - start_pos) as u16,
            };
            slices.push(interval);
        }

        let (last_pos, min_pos) = min_positions[min_positions.len() - 1];
        let last_interval = MspIntervalP {
            minimizer: min_pos.kmer,
            minimizer_pos: min_pos.pos as u32,
            start: last_pos as u32,
            len: (m - last_pos) as u16,
        };
        slices.push(last_interval);

        slices
    }
}

pub fn msp_sequence<P, V>(
    k: usize,
    seq: &[u8],
    permutation: Option<&[usize]>,
    rc: bool,
) -> Vec<(u32, Exts, V)>
where
    P: Kmer,
    V: Vmer,
{
    let p = P::k();

    // Make sure the substrings will fit into the Vmers
    assert!(V::max_len() >= 2 * k - p);

    if seq.len() < k {
        return Vec::new();
    }

    let default_perm: Option<Vec<usize>> = match permutation {
        Some(_) => None,
        None => Some((0..1 << (2 * p)).collect()),
    };

    let perm = permutation.unwrap_or_else(|| default_perm.as_ref().unwrap());

    let score = |pi: &P| {
        if rc {
            min(perm[pi.to_u64() as usize], perm[pi.rc().to_u64() as usize])
        } else {
            perm[pi.to_u64() as usize]
        }
    };

    let dna = DnaSlice(seq);
    let scanner = Scanner::new(&dna, score, k);
    let msp_parts = scanner.scan();

    let mut msps = Vec::new();
    for msp in msp_parts {
        let v = V::from_slice(&seq[(msp.start as usize)..(msp.start as usize + msp.len as usize)]);
        let exts = Exts::from_slice_bounds(seq, msp.start as usize, msp.len as usize);
        msps.push((msp.bucket() as u32, exts, v));
    }

    msps
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna_string::DnaString;
    use crate::kmer::*;
    use crate::test;
    use crate::DnaSlice;
    use std::collections::HashSet;
    use std::iter::FromIterator;

    fn all_kmers<T>(k: usize, seq: &[T]) -> Vec<&[T]> {
        (0..(seq.len() - k + 1)).map(|i| &seq[i..i + k]).collect()
    }

    fn test_all_kmers(k: usize, full_seq: &[u8], slices: Vec<MspInterval>) {
        let start_kmers = HashSet::from_iter(all_kmers(k, full_seq));

        let mut sliced_kmers = HashSet::new();

        for msp in slices {
            let slc = &full_seq[(msp.start as usize)..(msp.start as usize + msp.len as usize)];
            sliced_kmers.extend(all_kmers(k, slc));
        }

        if start_kmers != sliced_kmers {
            println!("start kmers: {:?}", start_kmers);
            println!("sliced kmers: {:?}", sliced_kmers);
            panic!("kmer sets not equal");
        }
    }

    #[test]
    fn test1() {
        let v = vec![1u8, 2u8, 3u8, 4u8, 5u8, 6u8];
        let s = &v[..];
        let ak = all_kmers(2, s);

        if ak[0] < ak[1] {
            println!("sorts!")
        }

        println!("{:?}", ak);

        let v = vec![6u8, 5u8, 4u8, 3u8, 2u8, 1u8, 0u8];
        let s = &v[..];
        let mut ak = all_kmers(2, s);

        println!("{:?}", ak);
        ak.sort();

        if ak[0] < ak[1] {
            println!("sorts!")
        }
        println!("{:?}", ak);
    }

    #[test]
    fn test_slice() {
        let p = 8;
        let permutation: Vec<usize> = (0..(1 << 2 * p)).collect();

        for _ in 0..100 {
            let k = 50usize;
            let dna = test::random_dna(150);
            println!("{:?}", dna);
            let slices = super::simple_scan::<_, Kmer8>(k, &DnaSlice(&dna), &permutation, true);
            println!(
                "Made {} slices from dna of length {:?}",
                slices.len(),
                dna.len()
            );

            println!("slices: {:?}", slices);
            test_all_kmers(k, &dna[..], slices);
        }
    }

    fn check_msp_slices<P, F>(
        k: usize,
        full_seq: &DnaString,
        slices: Vec<MspIntervalP<P>>,
        score: F,
    ) where
        P: Kmer,
        F: Fn(&P) -> usize,
    {
        // Check all the correctness properties of the slices returned by MSP.
        // 1. Each pmer is covered by one and only one slice
        // 2. Slices are at least p and most 2k-p long
        // 3. The selected pmer exists in the slice and is minimal given the scoring function
        // 4. No slices can be extended to the right by covering an additional pmer without violating the above constraints

        let p = P::k();

        /*
        For debugging:
        println!("k: {}, p: {}", k, p);
        println!("seq: {:?}", full_seq);
        println!("slices: {:#?}", slices);
        */

        // 1. each kmer is covered exactly once.
        let mut covered = vec![false; full_seq.len() - k + 1];

        for s in &slices {
            let start = s.start as usize;
            let end = s.start as usize + s.len as usize - k + 1;

            for i in start..end {
                if covered[i] {
                    println!("at {}", i);
                    assert!(false, "base already covered!");
                }

                covered[i] = true;
            }
        }
        assert!(covered.iter().all(|x| *x), "a pmer wasn't covered");

        //2. p <= slice.len <= 2k-p
        for s in &slices {
            assert!(s.len as usize >= p);
            assert!(s.len as usize <= 2 * k - p);
        }

        // 3. pmer exists and is the best possible within the slice
        for s in &slices {
            let slice_score = score(&s.minimizer);

            let start = s.start as usize;
            let end = s.start as usize + s.len as usize - p + 1;

            // check the correct minimizer info is reported
            assert_eq!(s.minimizer, full_seq.get_kmer(s.minimizer_pos as usize));

            for i in start..end {
                let pmer: P = full_seq.get_kmer(i);
                let score = score(&pmer);

                if score < slice_score {
                    assert!(false, "found better pmer within slice")
                }
            }
        }

        // 4. No slices can be extended to the right by covering an additional pmer without violating the above constraints
        for s in &slices[0..(slices.len() - 1)] {
            // the next kmer beyond the end of the slice must either:
            // a. not cover the minimizer
            // b. contain a better minimizer
            let next_kmer_pos = s.start as usize + s.len as usize - k + 1;
            let next_kmer_covers_pmer = next_kmer_pos <= s.minimizer_pos as usize;

            let next_pmer_pos = s.start as usize + s.len as usize - p + 1;
            let next_pmer: P = full_seq.get_kmer(next_pmer_pos);
            let next_score = score(&next_pmer);

            let correct_end = next_score < score(&s.minimizer) || !next_kmer_covers_pmer;

            if !correct_end {
                assert!(
                    false,
                    "detected a sequence slice that ended before it should have."
                )
            }
        }
    }

    fn test_new_slicer<P: Kmer>(k: usize) {
        let p = P::k();
        if p >= k {
            return;
        }

        let score = |p: &P| p.to_u64() as usize;

        for i in 0..(20 * k) {
            let len = 2 * i;
            if len < k {
                continue;
            }

            let dna = test::random_dna(len);
            let dna = DnaString::from_bytes(&dna);

            // use lexicographic ordering
            let scanner = Scanner::new(&dna, score, k);
            let slices = scanner.scan();
            check_msp_slices(k, &dna, slices, score);

            // use AT count ordering
            let at_score = |p: &P| p.at_count() as usize;
            let scanner = Scanner::new(&dna, at_score, k);
            let slices = scanner.scan();
            check_msp_slices(k, &dna, slices, at_score);
        }

        for i in 0..(4 * k) {
            let len = i;
            if len < k {
                continue;
            }

            let dna = DnaString::blank(i);

            let scanner = Scanner::new(&dna, score, k);
            let slices = scanner.scan();
            check_msp_slices(k, &dna, slices, score);
        }
    }

    #[test]
    fn test_msp_scanner() {
        for k in 16..64 {
            test_new_slicer::<Kmer5>(k);
            test_new_slicer::<Kmer8>(k);
            test_new_slicer::<Kmer10>(k);
            test_new_slicer::<Kmer12>(k);
            test_new_slicer::<Kmer14>(k);
            test_new_slicer::<Kmer15>(k);
            test_new_slicer::<Kmer16>(k);
        }
    }

    #[test]
    fn test_sample() {
        // for testing MSP on specific error cases

        // let v1 : Vec<u8> = vec![3, 0, 3, 0, 3, 2, 3, 1, 0, 0, 2, 0, 0, 1, 1, 3, 0, 2, 0, 3, 2, 3, 3, 1, 2, 0, 2, 2, 3, 1, 0, 1, 2, 1, 2, 3, 2, 1, 0, 2, 1, 0, 2, 2, 1, 0, 2, 3, 0, 3, 2, 3, 0, 2, 0, 0, 1, 0, 2, 3, 1, 3, 2, 0, 2, 2, 2, 2, 1, 2, 0, 2, 1, 0, 0, 1, 3, 0, 0, 0, 2, 2, 3, 0, 0, 3, 2, 3, 1, 2, 0, 2, 0, 3, 2, 1, 3, 2, 2, 3, 3, 1, 2, 3, 0, 3, 1, 1, 2, 2, 2, 3, 1, 1, 2, 0, 2, 3, 3, 2, 0, 3, 1, 1, 1, 1, 3, 3, 3, 0, 1, 0, 0, 3, 0, 2, 3, 2, 2, 2, 3, 2, 3, 1, 0, 2, 3, 0, 2, 2, 2, 1, 0, 2, 3, 3, 3, 2, 2, 1, 3, 1, 0, 2, 0, 1, 0, 1, 2, 2, 2, 2, 2, 3, 0, 3, 3, 1, 2, 3, 2, 0, 2, 1, 0, 3, 1, 1, 0, 2, 1, 3, 0, 1, 2, 1, 1, 1, 3, 1, 0, 3, 0, 0, 2, 2, 0, 3, 2, 3, 2, 0, 2, 1, 2, 3, 0, 1, 3, 1, 1, 2, 2, 3, 0, 3, 3, 3, 3, 0, 3, 3, 0, 0, 0, 2, 3, 3, 1, 3, 3, 1, 2, 2, 0, 3, 2, 1, 0, 2];
        // let v2 : Vec<u8> = vec![3, 0, 3, 0, 3, 2, 3, 1, 0, 0, 2, 0, 0, 1, 1, 3, 0, 2, 0, 3, 2, 3, 3, 1, 2, 0, 2, 2, 3, 2, 0, 3, 2, 1, 2, 3, 2, 1, 0, 2, 1, 0, 2, 2, 1, 0, 2, 3, 0, 3, 2, 3, 0, 2, 3, 0, 1, 0, 2, 3, 1, 3, 2, 0, 2, 2, 0, 2, 1, 2, 0, 2, 1, 0, 0, 1, 3, 0, 0, 0, 2, 2, 3, 0, 0, 3, 2, 3, 1, 2, 0, 2, 0, 3, 2, 1, 3, 2, 2, 3, 3, 1, 2, 3, 0, 3, 1, 1, 2, 2, 2, 3, 1, 1, 2, 0, 0, 3, 3, 2, 0, 3, 1, 1, 1, 1, 3, 3, 3, 0, 1, 0, 0, 3, 0, 2, 3, 1, 2, 1, 3, 2, 3, 1, 0, 2, 3, 0, 2, 2, 2, 2, 1, 2, 3, 3, 3, 2, 2, 1, 3, 1, 0, 2, 0, 1, 0, 1, 2, 2, 2, 2, 2, 3, 0, 3, 3, 1, 2, 3, 2, 0, 2, 1, 0, 3, 1, 1, 0, 2, 1, 3, 0, 1, 2, 1, 2, 1, 3, 1, 0, 3, 0, 0, 2, 2, 0, 3, 2, 3, 2, 0, 2, 1, 2, 3, 0, 1, 3, 1, 1, 2, 2, 3, 0, 3, 3, 3, 3, 0, 3, 3, 0, 0, 0, 2, 3, 3, 1, 0, 3, 1, 2, 2, 0, 3, 2, 1, 0, 2];

        let v1: Vec<u8> = vec![
            3, 0, 3, 0, 1, 2, 3, 3, 0, 0, 0, 1, 1, 0, 3, 3, 1, 2, 0, 1, 1, 3, 2, 1, 1, 1, 2, 3, 3,
            2, 1, 2, 2, 1, 2, 3, 2, 1, 0, 2, 1, 1, 2, 1, 0, 1, 2, 3, 0, 3, 2, 3, 0, 1, 3, 0, 1, 0,
            2, 3, 3, 3, 2, 3, 2, 2, 0, 2, 1, 2, 2, 2, 0, 0, 1, 0, 1, 2, 1, 0, 3, 2, 3, 1, 0, 3, 2,
            2, 2, 2, 1, 3, 0, 2, 1, 2, 1, 3, 3, 0, 1, 1, 2, 3, 0, 2, 3, 1, 3, 2, 3, 1, 1, 0, 2, 2,
            1, 1, 1, 2, 0, 0, 1, 2, 1, 0, 3, 3, 3, 0, 1, 1, 0, 3, 0, 2, 2, 1, 2, 1, 3, 2, 3, 1, 0,
            2, 3, 0, 2, 2, 0, 3, 3, 2, 3, 3, 0, 3, 0, 1, 1, 3, 0, 1, 0, 1, 0, 1, 2, 2, 0, 3, 3, 3,
            1, 3, 1, 1, 0, 1, 2, 0, 2, 1, 0, 3, 1, 1, 1, 2, 0, 3, 0, 1, 2, 1, 1, 1, 3, 1, 2, 3, 0,
            3, 2, 2, 0, 3, 2, 3, 3, 0, 2, 3, 2, 0, 0, 1, 2, 0, 3, 2, 2, 1, 2, 2, 3, 3, 3, 2, 3, 0,
            3, 1, 0, 2, 3, 3, 0, 3, 1, 0, 3, 2, 1, 3, 2, 1, 1, 2,
        ];
        let v2: Vec<u8> = vec![
            3, 0, 3, 0, 1, 3, 3, 3, 0, 0, 0, 1, 1, 0, 3, 3, 1, 2, 0, 1, 1, 3, 2, 1, 1, 1, 2, 1, 3,
            2, 1, 2, 2, 1, 2, 3, 2, 1, 0, 2, 1, 1, 2, 1, 1, 1, 2, 3, 0, 3, 2, 3, 0, 1, 3, 0, 1, 0,
            2, 3, 3, 3, 2, 3, 2, 2, 0, 2, 1, 2, 2, 2, 0, 0, 1, 0, 1, 2, 1, 0, 3, 2, 3, 1, 0, 3, 2,
            2, 2, 2, 1, 3, 0, 2, 1, 2, 1, 3, 3, 0, 1, 1, 2, 3, 0, 2, 3, 1, 3, 2, 3, 1, 1, 0, 2, 2,
            1, 1, 0, 2, 0, 0, 1, 2, 1, 0, 3, 3, 3, 0, 1, 1, 0, 3, 0, 2, 2, 1, 2, 1, 3, 2, 3, 1, 0,
            2, 3, 0, 2, 2, 0, 3, 3, 2, 3, 3, 0, 3, 0, 1, 1, 3, 0, 1, 0, 1, 0, 1, 2, 2, 0, 3, 3, 3,
            1, 3, 1, 1, 0, 1, 2, 0, 2, 1, 0, 3, 1, 1, 1, 2, 0, 3, 0, 1, 2, 1, 3, 1, 3, 1, 2, 3, 0,
            3, 2, 2, 0, 3, 2, 3, 3, 0, 2, 3, 2, 0, 0, 1, 2, 0, 3, 2, 2, 1, 2, 2, 3, 3, 3, 2, 3, 1,
            3, 1, 0, 2, 2, 3, 0, 3, 1, 0, 3, 3, 1, 3, 2, 1, 1, 2,
        ];

        let p = 5;
        let permutation: Vec<usize> = (0..(1 << 2 * p)).collect();

        let s1 = super::simple_scan::<_, Kmer5>(35, &DnaSlice(&mut &v1), &permutation, true);
        let s2 = super::simple_scan::<_, Kmer5>(35, &DnaSlice(&mut &v2), &permutation, true);

        println!("{:?}", s1);
        println!("{:?}", s2);
    }
}
