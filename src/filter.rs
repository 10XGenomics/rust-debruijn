// Copyright 2017 10x Genomics

//! Methods for converting sequences into kmers, filtering observed kmers before De Bruijn graph construction, and summarizing 'color' annotations.
use std::fmt::Debug;
use std::marker::PhantomData;
use std::mem;
use std::ops::Deref;

use boomphf::hashmap::BoomHashMap2;
use itertools::Itertools;
use log::debug;

use crate::Dir;
use crate::Exts;
use crate::Kmer;
use crate::Vmer;

fn bucket<K: Kmer>(kmer: K) -> usize {
    (kmer.get(0) as usize) << 6
        | (kmer.get(1) as usize) << 4
        | (kmer.get(2) as usize) << 2
        | (kmer.get(3) as usize)
}

/// Implement this trait to control how multiple observations of a kmer
/// are carried forward into a DeBruijn graph.
pub trait KmerSummarizer<DI, DO> {
    /// The input `items` is an iterator over kmer observations. Input observation
    /// is a tuple of (kmer, extensions, data). The summarize function inspects the
    /// data and returns a tuple indicating:
    /// * whether this kmer passes the filtering criteria (e.g. is there a sufficient number of observation)
    /// * the accumulated Exts of the kmer
    /// * a summary data object of type `DO` that will be used as a color annotation in the DeBruijn graph.
    fn summarize<K, F: Iterator<Item = (K, Exts, DI)>>(&self, items: F) -> (bool, Exts, DO);
}

/// A simple KmerSummarizer that only accepts kmers that are observed
/// at least a given number of times. The metadata returned about a Kmer
/// is the number of times it was observed, capped at 2^16.
pub struct CountFilter {
    min_kmer_obs: usize,
}

impl CountFilter {
    /// Construct a `CountFilter` KmerSummarizer only accepts kmers that are observed
    /// at least `min_kmer_obs` times.
    pub fn new(min_kmer_obs: usize) -> CountFilter {
        CountFilter {
            min_kmer_obs: min_kmer_obs,
        }
    }
}

impl<D> KmerSummarizer<D, u16> for CountFilter {
    fn summarize<K, F: Iterator<Item = (K, Exts, D)>>(&self, items: F) -> (bool, Exts, u16) {
        let mut all_exts = Exts::empty();
        let mut count = 0u16;
        for (_, exts, _) in items {
            count = count.saturating_add(1);
            all_exts = all_exts.add(exts);
        }

        (count as usize >= self.min_kmer_obs, all_exts, count)
    }
}

/// A simple KmerSummarizer that only accepts kmers that are observed
/// at least a given number of times. The metadata returned about a Kmer
/// is a vector of the unique data values observed for that kmer.
pub struct CountFilterSet<D> {
    min_kmer_obs: usize,
    phantom: PhantomData<D>,
}

impl<D> CountFilterSet<D> {
    /// Construct a `CountFilterSet` KmerSummarizer only accepts kmers that are observed
    /// at least `min_kmer_obs` times.
    pub fn new(min_kmer_obs: usize) -> CountFilterSet<D> {
        CountFilterSet {
            min_kmer_obs: min_kmer_obs,
            phantom: PhantomData,
        }
    }
}

impl<D: Ord> KmerSummarizer<D, Vec<D>> for CountFilterSet<D> {
    fn summarize<K, F: Iterator<Item = (K, Exts, D)>>(&self, items: F) -> (bool, Exts, Vec<D>) {
        let mut all_exts = Exts::empty();

        let mut out_data: Vec<D> = Vec::new();

        let mut nobs = 0;
        for (_, exts, d) in items {
            out_data.push(d);
            all_exts = all_exts.add(exts);
            nobs += 1;
        }

        out_data.sort();
        out_data.dedup();
        (nobs as usize >= self.min_kmer_obs, all_exts, out_data)
    }
}

/// Process DNA sequences into kmers and determine the set of valid kmers,
/// their extensions, and summarize associated label/'color' data. The input
/// sequences are converted to kmers of type `K`, and like kmers are grouped together.
/// All instances of each kmer, along with their label data are passed to
/// `summarizer`, an implementation of the `KmerSummarizer` which decides if
/// the kmer is 'valid' by an arbitrary predicate of the kmer data, and
/// summarizes the the individual label into a single label data structure
/// for the kmer. Care is taken to keep the memory consumption small.
/// Less than 4G of temporary memory should be allocated to hold intermediate kmers.
///
///
/// # Arguments
///
/// * `seqs` a slice of (sequence, extensions, data) tuples. Each tuple
///   represents an input sequence. The input sequence must implement `Vmer<K`> The data slot is an arbitrary data
///   structure labeling the input sequence.
///   If complete sequences are passed in, the extensions entry should be
///   set to `Exts::empty()`.
///   In sharded DBG construction (for example when minimizer-based partitioning
///   of the input strings), the input sequence is a sub-string of the original input string.
///   In this case the extensions of the sub-string in the original string
///   should be passed in the extensions.
/// * `summarizer` is an implementation of `KmerSummarizer<D1,DS>` that decides
///   whether a kmer is valid (e.g. based on the number of observation of the kmer),
///   and summarizes the data about the individual kmer observations. See `CountFilter`
///   and `CountFilterSet` for examples.
/// * `stranded`: if true, preserve the strandedness of the input sequences, effectively
///   assuming they are all in the positive strand. If false, the kmers will be canonicalized
///   to the lexicographic minimum of the kmer and it's reverse complement.
/// * `report_all_kmers`: if true returns the vector of all the observed kmers and performs the
///   kmer based filtering
/// * `memory_size`: gives the size bound on the memory in GB to use and automatically determines
///   the number of passes needed.
/// # Returns
/// BoomHashMap2 Object, check rust-boomphf for details
#[inline(never)]
pub fn filter_kmers<K: Kmer, V: Vmer, D1: Clone, DS, S: KmerSummarizer<D1, DS>>(
    seqs: &[(V, Exts, D1)],
    summarizer: &dyn Deref<Target = S>,
    stranded: bool,
    report_all_kmers: bool,
    memory_size: usize,
) -> (BoomHashMap2<K, Exts, DS>, Vec<K>)
where
    DS: Debug,
{
    let rc_norm = !stranded;

    let mut all_kmers = Vec::new();
    let mut valid_kmers = Vec::new();
    let mut valid_exts = Vec::new();
    let mut valid_data = Vec::new();

    // Estimate memory consumed by Kmer vectors, and set iteration count appropriately
    let input_kmers: usize = seqs
        .iter()
        .map(|&(ref vmer, _, _)| vmer.len().saturating_sub(K::k() - 1))
        .sum();
    let kmer_mem = input_kmers * mem::size_of::<(K, D1)>();
    let max_mem = memory_size * (10 as usize).pow(9);
    let slices = kmer_mem / max_mem + 1;
    let sz = 256 / slices + 1;

    let mut bucket_ranges = Vec::new();
    let mut start = 0;
    while start < 256 {
        bucket_ranges.push(start..start + sz);
        start += sz;
    }
    assert!(bucket_ranges[bucket_ranges.len() - 1].end >= 256);
    let n_buckets = bucket_ranges.len();

    if bucket_ranges.len() > 1 {
        debug!(
            "{} sequences, {} kmers, {} passes",
            seqs.len(),
            input_kmers,
            bucket_ranges.len()
        );
    }

    for (i, bucket_range) in bucket_ranges.into_iter().enumerate() {
        debug!("Processing bucket {} of {}", i, n_buckets);

        let mut kmer_buckets: Vec<Vec<(K, Exts, D1)>> = Vec::new();
        for _ in 0..256 {
            kmer_buckets.push(Vec::new());
        }

        for &(ref seq, seq_exts, ref d) in seqs {
            for (kmer, exts) in seq.iter_kmer_exts::<K>(seq_exts) {
                let (min_kmer, flip_exts) = if rc_norm {
                    let (min_kmer, flip) = kmer.min_rc_flip();
                    let flip_exts = if flip { exts.rc() } else { exts };
                    (min_kmer, flip_exts)
                } else {
                    (kmer, exts)
                };
                let bucket = bucket(min_kmer);

                if bucket >= bucket_range.start && bucket < bucket_range.end {
                    kmer_buckets[bucket].push((min_kmer, flip_exts, d.clone()));
                }
            }
        }

        for mut kmer_vec in kmer_buckets {
            kmer_vec.sort_by_key(|elt| elt.0);

            for (kmer, kmer_obs_iter) in &kmer_vec.into_iter().group_by(|elt| elt.0) {
                let (is_valid, exts, summary_data) = summarizer.summarize(kmer_obs_iter);
                if report_all_kmers {
                    all_kmers.push(kmer);
                }
                if is_valid {
                    valid_kmers.push(kmer);
                    valid_exts.push(exts);
                    valid_data.push(summary_data);
                }
            }
        }
    }

    debug!(
        "Unique kmers: {}, All kmers (if returned): {}",
        valid_kmers.len(),
        all_kmers.len(),
    );
    (
        BoomHashMap2::new(valid_kmers, valid_exts, valid_data),
        all_kmers,
    )
}

/// Remove extensions in valid_kmers that point to censored kmers. A censored kmer
/// exists in all_kmers but not valid_kmers. Since the kmer exists in this partition,
/// but was censored, we know that we can delete extensions to it.
/// In sharded kmer processing, we will have extensions to kmers in other shards. We don't
/// know whether these are censored until later, so we retain these extension.
pub fn remove_censored_exts_sharded<K: Kmer, D>(
    stranded: bool,
    valid_kmers: &mut Vec<(K, (Exts, D))>,
    all_kmers: &Vec<K>,
) {
    for idx in 0..valid_kmers.len() {
        let mut new_exts = Exts::empty();
        let kmer = valid_kmers[idx].0;
        let exts = (valid_kmers[idx].1).0;

        for dir in [Dir::Left, Dir::Right].iter() {
            for i in 0..4 {
                if exts.has_ext(*dir, i) {
                    let _ext_kmer = kmer.extend(i, *dir);

                    let ext_kmer = if stranded {
                        _ext_kmer
                    } else {
                        _ext_kmer.min_rc()
                    };

                    let censored = if valid_kmers.binary_search_by_key(&ext_kmer, |d| d.0).is_ok() {
                        // ext_kmer is valid. not censored
                        false
                    } else {
                        // ext_kmer is not valid. if it was in this shard, then we censor it
                        all_kmers.binary_search(&ext_kmer).is_ok()
                    };

                    if !censored {
                        new_exts = new_exts.set(*dir, i);
                    }
                }
            }
        }

        (valid_kmers[idx].1).0 = new_exts;
    }
}

/// Remove extensions in valid_kmers that point to censored kmers. Use this method in a non-partitioned
/// context when valid_kmers includes _all_ kmers that will ultimately be included in the graph.
pub fn remove_censored_exts<K: Kmer, D>(stranded: bool, valid_kmers: &mut Vec<(K, (Exts, D))>) {
    for idx in 0..valid_kmers.len() {
        let mut new_exts = Exts::empty();
        let kmer = valid_kmers[idx].0;
        let exts = (valid_kmers[idx].1).0;

        for dir in [Dir::Left, Dir::Right].iter() {
            for i in 0..4 {
                if exts.has_ext(*dir, i) {
                    let ext_kmer = if stranded {
                        kmer.extend(i, *dir)
                    } else {
                        kmer.extend(i, *dir).min_rc()
                    };

                    let kmer_valid = valid_kmers.binary_search_by_key(&ext_kmer, |d| d.0).is_ok();

                    if kmer_valid {
                        new_exts = new_exts.set(*dir, i);
                    }
                }
            }
        }

        (valid_kmers[idx].1).0 = new_exts;
    }
}
