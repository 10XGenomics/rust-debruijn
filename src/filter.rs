//! Methods for filtering observed kmers before De Bruijn graph construction, and summarizing 'color' annotations.

use std::mem;
use itertools::Itertools;
use std::marker::PhantomData;

use Dir;
use Kmer;
use Exts;
use Vmer;
use fx::FxHashMap;
use fx::FxHashSet;


fn bucket<K: Kmer>(kmer: K) -> usize {
    // FIXME - make 256 mins    
    kmer.get(0) as usize
}

pub trait KmerSummarizer<DI,DO> {
    fn summarize<K, F: Iterator<Item=(K,Exts,DI)>>(&self, items: F) -> (bool, Exts, DO);
}

struct CountFilter {
    min_kmer_obs: usize
}

impl CountFilter {
    pub fn new(min_kmer_obs: usize) -> CountFilter {
        CountFilter { min_kmer_obs: min_kmer_obs }
    }
}

impl<D> KmerSummarizer<D, u16> for CountFilter {
    fn summarize<K, F: Iterator<Item=(K,Exts,D)>>(&self, items: F) -> (bool, Exts, u16) {
        let mut all_exts = Exts::empty();
        let mut count = 0u16;
        for (_, exts, _) in items {
            count = count.saturating_add(1);
            all_exts = all_exts.add(exts);
        }

        (count as usize >= self.min_kmer_obs, all_exts, count)
    }
}

pub struct CountFilterSet<D> {
    min_kmer_obs: usize,
    phantom: PhantomData<D>,
}

impl<D> CountFilterSet<D> {
    pub fn new(min_kmer_obs: usize) -> CountFilterSet<D> {
        CountFilterSet { min_kmer_obs: min_kmer_obs, phantom: PhantomData }
    }
}


impl<D> KmerSummarizer<D, Vec<D>> for CountFilterSet<D> {
    fn summarize<K, F: Iterator<Item=(K,Exts,D)>>(&self, items: F) -> (bool, Exts, Vec<D>) {
        let mut all_exts = Exts::empty();
        
        let mut out_data: Vec<D> = Vec::new();
        
        for (_, exts, d) in items {
            out_data.push(d);
            all_exts = all_exts.add(exts);
        }

        (out_data.len() as usize >= self.min_kmer_obs, all_exts, out_data)
    }
}



/// Read a shard and determine the valid kmers
/// Low memory implementation that should consume < 4G of temporary memory
/// To reduce memory consumption, set track_bcs to false to forget about BC lists.
#[inline(never)]
pub fn filter_kmers<K:Kmer, V:Vmer<K>, D1: Clone, DS, S: KmerSummarizer<D1,DS>>(seqs: &Vec<(V, Exts, D1)>, summarizer: S) ->
 (Vec<(K, (Exts, DS))>, Vec<K>) {

    let mut all_kmers = Vec::new();
    let mut valid_kmers = Vec::new();

    // Estimate memory consumed by Kmer vectors, and set iteration count appropriately
    let input_kmers: usize = seqs.iter().map(|&(ref vmer, _, _)| vmer.len() - K::k() + 1).sum();
    let kmer_mem = input_kmers * mem::size_of::<(K, D1)>();
    let max_mem = 4 * (10 as usize).pow(9);
    let slices = kmer_mem / max_mem + 1;
    let sz = 256 / slices + 1;
    
    let mut bucket_ranges = Vec::new();
    let mut start = 0;
    while start < 256 {
        bucket_ranges.push(start..start+sz);
        start += sz;
    }

    //if bucket_ranges.len() > 1 {
    //    info!("processing {} BSPs in {} passes. Bucket ranges: {:?}", bsps.len(), bucket_ranges.len(), bucket_ranges);
    //}

    for bucket_range in bucket_ranges {

        let mut kmer_buckets: Vec<Vec<(K, Exts, D1)>> = Vec::new();
        for _ in 0..256 {
            kmer_buckets.push(Vec::new());
        }

        //info!("Enumerating kmers...");
        for &(ref seq, seq_exts, ref d) in seqs {
            for (kmer, exts) in seq.iter_kmer_exts(seq_exts) {
                let (min_kmer, flip) = kmer.min_rc_flip();
                let flip_exts = if flip { exts.rc() } else { exts };
                let bucket = bucket(min_kmer);

                if bucket >= bucket_range.start && bucket < bucket_range.end {
                    kmer_buckets[bucket].push((min_kmer, flip_exts, d.clone()));
                }
            }
        }

        //info!("Validating kmers...");
        for mut kmer_vec in kmer_buckets {

            kmer_vec.sort_by_key(|elt| elt.0);

            for (kmer, kmer_obs_iter) in &kmer_vec.into_iter().group_by(|elt| elt.0) {
                let (is_valid, exts, summary_data) = summarizer.summarize(kmer_obs_iter);
                all_kmers.push(kmer);
                if is_valid {valid_kmers.push((kmer, (exts, summary_data))); }
            }
        }
    }

    valid_kmers.sort_by_key(|x| x.0);
    all_kmers.sort();
    fix_exts(&mut valid_kmers, &all_kmers);

    //info!("Total Sequences: {}, Total kmers observed: {}, Unique Kmers observed: {}. Kmers accepted: {}", seqs.len(), total_kmers, unique_kmers, final_kmers.len());
    (valid_kmers, all_kmers)
}

/// Remove extensions in valid_kmers that point to censored kmers. A censored kmer
/// exists in all_kmers but not valid_kmers. We know that we can delete these extensions.
/// In sharded kmer processing, we will have extensions to kmers in other shards. We don't
/// know whether these are censored until later, so we retain the extension.
pub fn fix_exts<K: Kmer, D>(valid_kmers: &mut Vec<(K, (Exts, D))>, all_kmers: &Vec<K>) {

    for idx in 0 .. valid_kmers.len() {
        let mut new_exts = Exts::empty();
        let kmer = valid_kmers[idx].0;
        let exts = (valid_kmers[idx].1).0;

        for dir in [Dir::Left, Dir::Right].iter() {
            for i in 0..4
            {
                if exts.has_ext(*dir, i) {
                    let ext_kmer = kmer.extend(i, *dir);
                    let kmer_observed = all_kmers.binary_search(&ext_kmer).is_ok();
                    let kmer_valid = valid_kmers.binary_search_by_key(&ext_kmer, |d| d.0).is_ok();

                    let censored = kmer_observed && !kmer_valid;

                    if !censored {
                        new_exts = new_exts.set(*dir, i);
                    }
                }
            }
        }

        (valid_kmers[idx].1).0 = new_exts;
    }
}