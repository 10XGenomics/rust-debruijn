// Copyright 2017 10x Genomics

//! Methods for minimum substring partitioning of a DNA string
//!
//! simple_scan method is based on:
//! Li, Yang. "MSPKmerCounter: a fast and memory efficient approach for k-mer counting." arXiv preprint arXiv:1505.06550 (2015).

use std::cmp::min;
use std::iter::Iterator;
use std::ops::Range;
use crate::Kmer;
use crate::Vmer;
use crate::Exts;
use crate::DnaSlice;

#[derive(Debug)]
pub struct MspInterval {
    bucket: u16,
    start: u32,
    len: u16,
}

impl MspInterval {
    pub fn new(bucket: u16, start: u32, len: u16)
               -> MspInterval{
        MspInterval{
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
        self.start() .. self.start() + self.len()
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
pub fn simple_scan<V: Vmer, P: Kmer>(
    k: usize,
    seq: &V,
    permutation: &[usize],
    rc: bool,
) -> Vec<MspInterval> {

    // Can't partition strings shorter than k
    assert!(seq.len() >= k);
    assert!(P::k() <= 8);
    assert!(seq.len() < 1<<32);

    let p = P::k();

    let pv = |pi: P| {
        if rc {
            min(
                permutation[pi.to_u64() as usize],
                permutation[pi.rc().to_u64() as usize])
        } else {
            permutation[pi.to_u64() as usize]
        }
    };

    let pval = |i: usize| {
        let pi = seq.get_kmer::<P>(i);
        pv(pi)
    };


    let pmin = |i: usize, j: usize| if pval(i) <= pval(j) { i } else { j };

    let m = seq.len();

    let find_min = |start, stop| {
        let mut pos = start;
        let mut min_pos = start;
        let mut pmer = seq.get_kmer(pos);
        let mut min_val = pv(pmer);

        while pos < stop {
            pos = pos + 1;
            pmer = pmer.extend_right(seq.get(pos + p - 1));
            let val = pv(pmer);

            if val < min_val {
                min_pos = pos;
                min_val = val;
            }
        }

        min_pos
    };

    let mut min_positions = Vec::with_capacity(16);
    let mut min_pos = find_min(0, k - p);
    min_positions.push((0, min_pos));

    for i in 0..(m - k + 1) {
        if i > min_pos {
            min_pos = find_min(i, i + k - p);
            min_positions.push((i, min_pos));
        } else {
            let test_min = pmin(min_pos, i + k - p);

            if test_min != min_pos {
                min_pos = test_min;
                min_positions.push((i, min_pos));
            }
        }
    }

    let mut slices = Vec::with_capacity(min_positions.len());

    // Generate the slices of the final string
    for p in 0..min_positions.len() - 1 {
        let (start_pos, min_pos) = min_positions[p];
        let (next_pos, _) = min_positions[p + 1];

        let interval = MspInterval {
            bucket: pval(min_pos) as u16,
            start: start_pos as u32,
            len: (next_pos + k - 1 - start_pos) as u16
        };
        slices.push(interval);
    }

    let (last_pos, min_pos) = min_positions[min_positions.len() - 1];
    let last_interval = MspInterval {
        bucket: pval(min_pos) as u16,
        start: last_pos as u32,
        len: (m - last_pos) as u16,
    };
    slices.push(last_interval);

    slices
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

    let default_perm: Option<Vec<usize>> =
        match permutation {
            Some(_) => None,
            None => Some( (0..1<<(2*p)).collect()),
        };

    let perm = permutation.unwrap_or_else(|| default_perm.as_ref().unwrap());

    let msp_parts = simple_scan::<_, P>(k, &DnaSlice(seq), perm, rc);

    let mut msps = Vec::new();
    for msp in msp_parts {
        let v = V::from_slice(&seq[(msp.start as usize)..(msp.start as usize + msp.len as usize)]);
        let exts = Exts::from_slice_bounds(seq, msp.start as usize, msp.len as usize);
        msps.push((msp.bucket as u32, exts, v));
    }

    msps
}



#[cfg(test)]
mod tests {
    use crate::test;
    use super::*;
    use std::collections::HashSet;
    use std::iter::FromIterator;
    use crate::kmer::{Kmer8, Kmer5};
    use crate::DnaSlice;

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

    #[test]
    fn test_sample() {
        // for testing MSP on specific error cases

        // let v1 : Vec<u8> = vec![3, 0, 3, 0, 3, 2, 3, 1, 0, 0, 2, 0, 0, 1, 1, 3, 0, 2, 0, 3, 2, 3, 3, 1, 2, 0, 2, 2, 3, 1, 0, 1, 2, 1, 2, 3, 2, 1, 0, 2, 1, 0, 2, 2, 1, 0, 2, 3, 0, 3, 2, 3, 0, 2, 0, 0, 1, 0, 2, 3, 1, 3, 2, 0, 2, 2, 2, 2, 1, 2, 0, 2, 1, 0, 0, 1, 3, 0, 0, 0, 2, 2, 3, 0, 0, 3, 2, 3, 1, 2, 0, 2, 0, 3, 2, 1, 3, 2, 2, 3, 3, 1, 2, 3, 0, 3, 1, 1, 2, 2, 2, 3, 1, 1, 2, 0, 2, 3, 3, 2, 0, 3, 1, 1, 1, 1, 3, 3, 3, 0, 1, 0, 0, 3, 0, 2, 3, 2, 2, 2, 3, 2, 3, 1, 0, 2, 3, 0, 2, 2, 2, 1, 0, 2, 3, 3, 3, 2, 2, 1, 3, 1, 0, 2, 0, 1, 0, 1, 2, 2, 2, 2, 2, 3, 0, 3, 3, 1, 2, 3, 2, 0, 2, 1, 0, 3, 1, 1, 0, 2, 1, 3, 0, 1, 2, 1, 1, 1, 3, 1, 0, 3, 0, 0, 2, 2, 0, 3, 2, 3, 2, 0, 2, 1, 2, 3, 0, 1, 3, 1, 1, 2, 2, 3, 0, 3, 3, 3, 3, 0, 3, 3, 0, 0, 0, 2, 3, 3, 1, 3, 3, 1, 2, 2, 0, 3, 2, 1, 0, 2];
        // let v2 : Vec<u8> = vec![3, 0, 3, 0, 3, 2, 3, 1, 0, 0, 2, 0, 0, 1, 1, 3, 0, 2, 0, 3, 2, 3, 3, 1, 2, 0, 2, 2, 3, 2, 0, 3, 2, 1, 2, 3, 2, 1, 0, 2, 1, 0, 2, 2, 1, 0, 2, 3, 0, 3, 2, 3, 0, 2, 3, 0, 1, 0, 2, 3, 1, 3, 2, 0, 2, 2, 0, 2, 1, 2, 0, 2, 1, 0, 0, 1, 3, 0, 0, 0, 2, 2, 3, 0, 0, 3, 2, 3, 1, 2, 0, 2, 0, 3, 2, 1, 3, 2, 2, 3, 3, 1, 2, 3, 0, 3, 1, 1, 2, 2, 2, 3, 1, 1, 2, 0, 0, 3, 3, 2, 0, 3, 1, 1, 1, 1, 3, 3, 3, 0, 1, 0, 0, 3, 0, 2, 3, 1, 2, 1, 3, 2, 3, 1, 0, 2, 3, 0, 2, 2, 2, 2, 1, 2, 3, 3, 3, 2, 2, 1, 3, 1, 0, 2, 0, 1, 0, 1, 2, 2, 2, 2, 2, 3, 0, 3, 3, 1, 2, 3, 2, 0, 2, 1, 0, 3, 1, 1, 0, 2, 1, 3, 0, 1, 2, 1, 2, 1, 3, 1, 0, 3, 0, 0, 2, 2, 0, 3, 2, 3, 2, 0, 2, 1, 2, 3, 0, 1, 3, 1, 1, 2, 2, 3, 0, 3, 3, 3, 3, 0, 3, 3, 0, 0, 0, 2, 3, 3, 1, 0, 3, 1, 2, 2, 0, 3, 2, 1, 0, 2];

        let v1: Vec<u8> = vec![
            3,
            0,
            3,
            0,
            1,
            2,
            3,
            3,
            0,
            0,
            0,
            1,
            1,
            0,
            3,
            3,
            1,
            2,
            0,
            1,
            1,
            3,
            2,
            1,
            1,
            1,
            2,
            3,
            3,
            2,
            1,
            2,
            2,
            1,
            2,
            3,
            2,
            1,
            0,
            2,
            1,
            1,
            2,
            1,
            0,
            1,
            2,
            3,
            0,
            3,
            2,
            3,
            0,
            1,
            3,
            0,
            1,
            0,
            2,
            3,
            3,
            3,
            2,
            3,
            2,
            2,
            0,
            2,
            1,
            2,
            2,
            2,
            0,
            0,
            1,
            0,
            1,
            2,
            1,
            0,
            3,
            2,
            3,
            1,
            0,
            3,
            2,
            2,
            2,
            2,
            1,
            3,
            0,
            2,
            1,
            2,
            1,
            3,
            3,
            0,
            1,
            1,
            2,
            3,
            0,
            2,
            3,
            1,
            3,
            2,
            3,
            1,
            1,
            0,
            2,
            2,
            1,
            1,
            1,
            2,
            0,
            0,
            1,
            2,
            1,
            0,
            3,
            3,
            3,
            0,
            1,
            1,
            0,
            3,
            0,
            2,
            2,
            1,
            2,
            1,
            3,
            2,
            3,
            1,
            0,
            2,
            3,
            0,
            2,
            2,
            0,
            3,
            3,
            2,
            3,
            3,
            0,
            3,
            0,
            1,
            1,
            3,
            0,
            1,
            0,
            1,
            0,
            1,
            2,
            2,
            0,
            3,
            3,
            3,
            1,
            3,
            1,
            1,
            0,
            1,
            2,
            0,
            2,
            1,
            0,
            3,
            1,
            1,
            1,
            2,
            0,
            3,
            0,
            1,
            2,
            1,
            1,
            1,
            3,
            1,
            2,
            3,
            0,
            3,
            2,
            2,
            0,
            3,
            2,
            3,
            3,
            0,
            2,
            3,
            2,
            0,
            0,
            1,
            2,
            0,
            3,
            2,
            2,
            1,
            2,
            2,
            3,
            3,
            3,
            2,
            3,
            0,
            3,
            1,
            0,
            2,
            3,
            3,
            0,
            3,
            1,
            0,
            3,
            2,
            1,
            3,
            2,
            1,
            1,
            2,
        ];
        let v2: Vec<u8> = vec![
            3,
            0,
            3,
            0,
            1,
            3,
            3,
            3,
            0,
            0,
            0,
            1,
            1,
            0,
            3,
            3,
            1,
            2,
            0,
            1,
            1,
            3,
            2,
            1,
            1,
            1,
            2,
            1,
            3,
            2,
            1,
            2,
            2,
            1,
            2,
            3,
            2,
            1,
            0,
            2,
            1,
            1,
            2,
            1,
            1,
            1,
            2,
            3,
            0,
            3,
            2,
            3,
            0,
            1,
            3,
            0,
            1,
            0,
            2,
            3,
            3,
            3,
            2,
            3,
            2,
            2,
            0,
            2,
            1,
            2,
            2,
            2,
            0,
            0,
            1,
            0,
            1,
            2,
            1,
            0,
            3,
            2,
            3,
            1,
            0,
            3,
            2,
            2,
            2,
            2,
            1,
            3,
            0,
            2,
            1,
            2,
            1,
            3,
            3,
            0,
            1,
            1,
            2,
            3,
            0,
            2,
            3,
            1,
            3,
            2,
            3,
            1,
            1,
            0,
            2,
            2,
            1,
            1,
            0,
            2,
            0,
            0,
            1,
            2,
            1,
            0,
            3,
            3,
            3,
            0,
            1,
            1,
            0,
            3,
            0,
            2,
            2,
            1,
            2,
            1,
            3,
            2,
            3,
            1,
            0,
            2,
            3,
            0,
            2,
            2,
            0,
            3,
            3,
            2,
            3,
            3,
            0,
            3,
            0,
            1,
            1,
            3,
            0,
            1,
            0,
            1,
            0,
            1,
            2,
            2,
            0,
            3,
            3,
            3,
            1,
            3,
            1,
            1,
            0,
            1,
            2,
            0,
            2,
            1,
            0,
            3,
            1,
            1,
            1,
            2,
            0,
            3,
            0,
            1,
            2,
            1,
            3,
            1,
            3,
            1,
            2,
            3,
            0,
            3,
            2,
            2,
            0,
            3,
            2,
            3,
            3,
            0,
            2,
            3,
            2,
            0,
            0,
            1,
            2,
            0,
            3,
            2,
            2,
            1,
            2,
            2,
            3,
            3,
            3,
            2,
            3,
            1,
            3,
            1,
            0,
            2,
            2,
            3,
            0,
            3,
            1,
            0,
            3,
            3,
            1,
            3,
            2,
            1,
            1,
            2,
        ];

        let p = 5;
        let permutation: Vec<usize> = (0..(1 << 2 * p)).collect();

        let s1 = super::simple_scan::<_, Kmer5>(35, &DnaSlice(&mut &v1), &permutation, true);
        let s2 = super::simple_scan::<_, Kmer5>(35, &DnaSlice(&mut &v2), &permutation, true);

        println!("{:?}", s1);
        println!("{:?}", s2);
    }
}
