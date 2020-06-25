// Copyright 2017 10x Genomics

//! Generate random genomes (with lots of re-used sustrings), reassemble them, and check sanity

use crate::complement;
use crate::Kmer;
use crate::Vmer;

use rand::distributions::{Distribution, Gamma, Range};
use rand::{self, Rng, RngCore};
use std::cmp::{max, min};

/// Generate a uniformly random base
pub fn random_base() -> u8 {
    let mut r = rand::thread_rng();
    (r.next_u64() % 4) as u8
}

/// Generate uniformly random DNA sequences
pub fn random_dna(len: usize) -> Vec<u8> {
    let mut r = rand::thread_rng();
    let mut dna = Vec::new();
    for _ in 0..len {
        let b = (r.next_u64() % 4) as u8;
        dna.push(b);
    }

    dna
}

/// Randomly mutate each base with probability `p`
pub fn edit_dna<R: Rng>(seq: &mut Vec<u8>, p: f64, r: &mut R) {
    for b in seq.iter_mut() {
        if r.gen_range(0.0, 1.0) < p {
            *b = random_base();
        }
    }
}

pub fn random_kmer<K: Kmer>() -> K {
    let mut r = rand::thread_rng();
    let mut kmer = K::empty();
    for pos in 0..K::k() {
        let b = (r.next_u64() % 4) as u8;
        kmer.set_mut(pos, b);
    }
    kmer
}

pub fn random_vmer<K: Kmer, V: Vmer>() -> V {
    let mut r = rand::thread_rng();
    let len = r.gen_range(K::k(), min(200, V::max_len()));
    let mut lmer = V::new(len);

    for pos in 0..len {
        let b = (r.next_u64() % 4) as u8;
        lmer.set_mut(pos, b);
    }
    lmer
}

pub fn simple_random_contigs() -> Vec<Vec<u8>> {
    let p1 = random_dna(40);
    let p2 = random_dna(30);

    let pc = random_dna(100);

    let p3 = random_dna(30);
    let p4 = random_dna(40);

    // Simulate contigs
    let mut c1 = Vec::new();
    c1.extend(p1);
    c1.extend(pc.clone());
    c1.extend(p3);

    let mut c2 = Vec::new();
    c2.extend(p2);
    c2.extend(pc);
    c2.extend(p4);

    let mut c3 = Vec::new();
    c3.extend(random_dna(30));

    // Stick a palindrome in the middle
    let palindrome1 = random_dna(33);
    let mut palindrome2 = palindrome1.clone();
    palindrome2.reverse();
    for i in 0..palindrome2.len() {
        palindrome2[i] = complement(palindrome2[i]);
    }

    c3.extend(palindrome1);
    c3.extend(palindrome2);
    c3.extend(random_dna(50));

    let contigs = vec![c1, c2, c3];
    contigs
}

// Generate random contigs with complicated repeats
pub fn random_contigs() -> Vec<Vec<u8>> {
    // Generate a bunch of sequence chunks
    let mut rng = rand::thread_rng();

    let gamma_dist = Gamma::new(0.6, 25.0);

    let nchunks = max(5, gamma_dist.sample(&mut rng) as u32);
    let chunk_sample = Range::new(0, nchunks);

    let length_dist = Gamma::new(1.5, 200.0);

    let mut chunks: Vec<Vec<u8>> = Vec::new();
    for _ in 0..nchunks {
        let len = max(10, length_dist.sample(&mut rng) as usize);
        let seq = random_dna(len);
        chunks.push(seq);
    }

    // Now make a bunch of chromosomes by pasting together chunks
    let nchrom = max(4, gamma_dist.sample(&mut rng) as u32);

    let mut chroms = Vec::new();
    for _ in 0..nchrom {
        let chrom_chunks = max(4, gamma_dist.sample(&mut rng) as u32);

        let mut chrom_seq = Vec::new();
        for _ in 0..chrom_chunks {
            let chunk_idx = chunk_sample.sample(&mut rng) as usize;
            chrom_seq.extend(chunks[chunk_idx].clone());
        }
        chroms.push(chrom_seq);
    }

    chroms
}

#[cfg(test)]
mod tests {

    use crate::clean_graph::CleanGraph;
    use crate::compression::{compress_graph, compress_kmers_with_hash, SimpleCompress};
    use crate::graph::BaseGraph;
    use crate::DnaBytes;
    use crate::{Dir, Exts, Kmer};
    use boomphf::hashmap::BoomHashMap2;
    use boomphf::Mphf;
    use std::collections::{HashMap, HashSet};
    use std::iter::FromIterator;

    use crate::dna_string::DnaString;
    use crate::filter;
    use crate::kmer::Kmer6;
    use crate::kmer::{IntKmer, VarIntKmer, K31};
    use crate::msp;
    use std::ops::Sub;

    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn simple_kmer_compress() {
        let contigs = simple_random_contigs();
        reassemble_contigs::<IntKmer<u64>, DnaString>(contigs, false);
    }

    #[test]
    fn simple_sharded() {
        let contigs = simple_random_contigs();
        reassemble_sharded::<IntKmer<u64>, DnaString>(contigs, false);
    }

    #[test]
    fn degen_seq_asm() {
        let ctg = "AAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAA";
        let seq: Vec<u8> = ctg
            .as_bytes()
            .iter()
            .cloned()
            .map(crate::base_to_bits)
            .collect();

        reassemble_contigs::<VarIntKmer<u64, K31>, DnaString>(vec![seq.clone(), seq], false);
    }

    #[test]
    fn degen_seq_asm_sharded() {
        let ctg = "AAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAA";
        let seq: Vec<u8> = ctg
            .as_bytes()
            .iter()
            .cloned()
            .map(crate::base_to_bits)
            .collect();

        reassemble_sharded::<VarIntKmer<u64, K31>, DnaString>(vec![seq.clone(), seq], false);
    }

    #[test]
    fn complex_kmer_compress() {
        for _ in 0..10 {
            let contigs = random_contigs();
            reassemble_contigs::<IntKmer<u64>, DnaString>(contigs, false);
        }
    }

    #[test]
    fn complex_sharded() {
        for _ in 0..10 {
            let contigs = random_contigs();
            reassemble_sharded::<IntKmer<u64>, DnaString>(contigs, false);
        }
    }

    #[test]
    fn simple_path_compress() {
        let contigs = simple_random_contigs();
        simplify_from_kmers::<IntKmer<u64>>(contigs, false);
    }

    #[test]
    fn complex_path_compress_k31() {
        for _ in 0..100 {
            let contigs = random_contigs();
            simplify_from_kmers::<VarIntKmer<u64, K31>>(contigs, false);
        }
    }

    #[test]
    fn complex_path_compress() {
        for _ in 0..10 {
            let contigs = random_contigs();
            simplify_from_kmers::<IntKmer<u64>>(contigs, false);
        }
    }

    fn simplify_from_kmers<K: Kmer + Send + Sync>(mut contigs: Vec<Vec<u8>>, stranded: bool) {
        let seqs: Vec<(DnaBytes, Exts, ())> = contigs
            .drain(..)
            .map(|x| (DnaBytes(x), Exts::empty(), ()))
            .collect();
        let (valid_kmers, _): (BoomHashMap2<K, Exts, u16>, _) = filter::filter_kmers(
            &seqs,
            &Box::new(filter::CountFilter::new(1)),
            stranded,
            false,
            4,
        );

        let spec =
            SimpleCompress::new(|d1: u16, d2: &u16| ((d1 as u32 + *d2 as u32) % 65535) as u16);
        let from_kmers = compress_kmers_with_hash(stranded, &spec, &valid_kmers).finish();
        let is_cmp = from_kmers.is_compressed(&spec);
        if is_cmp.is_some() {
            println!("not compressed: nodes: {:?}", is_cmp);
            from_kmers.print();
        }
        assert!(from_kmers.is_compressed(&spec) == None);

        // Create a DBG with one node per input kmer
        let mut base_graph: BaseGraph<K, u16> = BaseGraph::new(stranded);

        for (kmer, exts, _) in valid_kmers.iter() {
            base_graph.add(kmer.iter(), *exts, 1);
        }
        let uncompressed_dbg = base_graph.finish();

        // Canonicalize the graph with
        let spec = SimpleCompress::new(|d1: u16, d2: &u16| d1 + d2);
        let simp_dbg = compress_graph(stranded, &spec, uncompressed_dbg, None);

        let is_cmp = simp_dbg.is_compressed(&spec);
        if is_cmp.is_some() {
            println!("not compressed: nodes: {:?}", is_cmp);
            simp_dbg.print();
        }

        assert!(simp_dbg.is_compressed(&spec) == None);

        let total_kmers = valid_kmers.len();

        // Test the Boomphf DBG indexing machinery
        // Make an MPHF of the kmers in the DBG.
        // Each kmer should hash to a unique slot.
        let mphf = Mphf::from_chunked_iterator_parallel(1.7, &simp_dbg, None, valid_kmers.len(), 2);

        let mut got_slot = vec![false; total_kmers];

        for n in simp_dbg.iter_nodes() {
            for kmer in n.sequence().iter_kmers() {
                let r = mphf.try_hash(&kmer).unwrap() as usize;
                assert_eq!(got_slot[r], false);
                got_slot[r] = true;
            }
        }

        assert!(got_slot.iter().all(|x| *x));
    }

    // Take some input contig, which likely form a complicated graph,
    // and test the kmer, bsp, sedge and edge construction machinery
    fn reassemble_contigs<K: Kmer + Copy, V: Vmer + Clone>(contigs: Vec<Vec<u8>>, stranded: bool) {
        let ctg_lens: Vec<_> = contigs.iter().map(|c| c.len()).collect();
        println!("Reassembling contig sizes: {:?}", ctg_lens);

        // kmer vector
        let mut kmers = Vec::new();
        for c in contigs.iter() {
            let mut _kmers = K::kmers_from_bytes(&c);
            kmers.extend(_kmers.iter().map(|k| k.min_rc()));
        }

        // True kmer set
        let mut kmer_set = HashSet::new();
        kmer_set.extend(kmers.iter());

        // Bsps of kmers
        let p = 6;

        let mut seqs: Vec<(V, Exts, u8)> = Vec::new();
        let permutation: Vec<usize> = (0..1 << (2 * p)).collect();
        for c in contigs.iter() {
            let msps =
                msp::msp_sequence::<Kmer6, V>(K::k(), c.as_slice(), Some(&permutation), true);
            seqs.extend(msps.clone().into_iter().map(|(_, e, v)| (v, e, 0u8)));
            seqs.extend(msps.into_iter().map(|(_, e, v)| (v, e, 1u8)));
        }

        // kmer set from bsps
        let mut msp_kmers = HashSet::new();
        for &(ref v, _, _) in seqs.iter() {
            for k in v.iter_kmers::<K>() {
                msp_kmers.insert(k.min_rc());
            }
        }

        // Kmer extensions from BSP match raw kmers
        if kmer_set != msp_kmers {
            println!("{:?}", kmer_set);
            println!("{:?}", msp_kmers);
        }

        // Raw kmers and BSP kmers match
        assert!(kmer_set == msp_kmers);

        // Check the correctness of the process_kmer_shard kmer filtering function
        let (valid_kmers, _): (BoomHashMap2<K, Exts, u16>, _) = filter::filter_kmers(
            &seqs,
            &Box::new(filter::CountFilter::new(2)),
            stranded,
            false,
            4,
        );
        let mut process_kmer_set: HashSet<K> = HashSet::new();
        for k in valid_kmers.iter().map(|x| x.0) {
            process_kmer_set.insert(k.clone());
        }
        assert_eq!(process_kmer_set, kmer_set);

        // Every kmer should be reachable as the extension of a kmer.
        // No new kmers should be reachable
        let mut extension_kmer_set: HashSet<K> = HashSet::new();
        for (kmer, exts, _) in &valid_kmers {
            for e in kmer.get_extensions(*exts, Dir::Left) {
                extension_kmer_set.insert(e.min_rc());
            }

            for e in kmer.get_extensions(*exts, Dir::Right) {
                extension_kmer_set.insert(e.min_rc());
            }
        }

        // Kmer extensions from BSP match raw kmers
        if kmer_set != extension_kmer_set {
            println!("n:{}, {:?}", kmer_set.len(), kmer_set);
            println!("n:{}, {:?}", extension_kmer_set.len(), extension_kmer_set);

            if extension_kmer_set.is_superset(&kmer_set) {
                let invented = extension_kmer_set.sub(&kmer_set);
                println!("Invented kmers: {:?}", invented);
            }
        }

        assert!(kmer_set == extension_kmer_set);

        let spec = SimpleCompress::new(|d1: u16, d2: &u16| d1.saturating_add(*d2));

        // Generate compress DBG for these kmers
        let graph = compress_kmers_with_hash(stranded, &spec, &valid_kmers);

        // Check that all the lines have valid kmers,
        // and have extensions into other valid kmers
        let mut all_contig_set = HashSet::new();

        for seq_id in 0..graph.len() {
            let seq = graph.sequences.get(seq_id);
            let exts = graph.exts[seq_id];

            let contig_set = HashSet::from_iter(seq.iter_kmers().map(|km: K| km.min_rc()));
            all_contig_set.extend(contig_set.iter());
            assert!(kmer_set.is_superset(&contig_set));

            for l_ext in exts.get(Dir::Left) {
                let f: K = seq.first_kmer();
                let ext_kmer: K = f.extend_left(l_ext);
                assert!(kmer_set.contains(&(ext_kmer.min_rc())));
            }

            for r_ext in exts.get(Dir::Right) {
                let f: K = seq.last_kmer();
                let ext_kmer: K = f.extend_right(r_ext);
                assert!(kmer_set.contains(&(ext_kmer.min_rc())));
            }
        }

        assert_eq!(kmer_set, all_contig_set);
    }

    // Take some input contig, which likely form a complicated graph,
    // and the msp / shard_asm / main_asm loop
    fn reassemble_sharded<K: Kmer + Copy + Sync + Send, V: Vmer + Clone>(
        contigs: Vec<Vec<u8>>,
        stranded: bool,
    ) {
        let ctg_lens: Vec<_> = contigs.iter().map(|c| c.len()).collect();
        println!("Reassembling contig sizes: {:?}", ctg_lens);

        // kmer vector
        let mut kmer_set = HashSet::new();
        for c in contigs.iter() {
            let mut _kmers = K::kmers_from_bytes(&c);
            kmer_set.extend(_kmers.iter().map(|k| k.min_rc()));
        }

        // Bsps of kmers
        let mut shards = HashMap::new();

        for ctg in contigs.iter() {
            let msps = msp::msp_sequence::<Kmer6, V>(K::k(), ctg.as_slice(), None, true);

            for (shard, exts, seq) in msps {
                let shard_vec = shards.entry(shard).or_insert_with(|| Vec::new());
                shard_vec.push((seq.clone(), exts, 0u8));
                shard_vec.push((seq, exts, 1u8));
            }
        }

        let mut shard_asms = Vec::new();

        // Do a subassembly in each shard
        for seqs in shards.values() {
            // Check the correctness of the process_kmer_shard kmer filtering function
            let (valid_kmers, _): (BoomHashMap2<K, Exts, u16>, _) = filter::filter_kmers(
                &seqs,
                &Box::new(filter::CountFilter::new(2)),
                stranded,
                false,
                4,
            );

            // Generate compress DBG for this shard
            let spec = SimpleCompress::new(|d1: u16, d2: &u16| d1.saturating_add(*d2));

            //print!("{:?}", valid_kmers);
            let graph = compress_kmers_with_hash(stranded, &spec, &valid_kmers);
            shard_asms.push(graph.clone());
            //graph.finish().print();
        }

        // Shove the subassemblies into a partially compress base graph
        let combined_graph = BaseGraph::combine(shard_asms.into_iter()).finish();
        let cmp = SimpleCompress::new(|a: u16, b: &u16| max(a, *b));
        let dbg_graph = compress_graph(false, &cmp, combined_graph, None);

        // Switch on for debugging
        //dbg_graph.print();
        //dbg_graph.write_gfa(&mut std::io::stdout().lock());

        let graph = dbg_graph.base;

        // Check that all the lines have valid kmers,
        // and have extensions into other valid kmers
        let mut all_contig_set = HashSet::new();

        for seq_id in 0..graph.len() {
            let seq = graph.sequences.get(seq_id);
            let exts = graph.exts[seq_id];

            let contig_set = HashSet::from_iter(seq.iter_kmers().map(|km: K| km.min_rc()));
            all_contig_set.extend(contig_set.iter());
            assert!(kmer_set.is_superset(&contig_set));

            for l_ext in exts.get(Dir::Left) {
                let f: K = seq.first_kmer();
                let ext_kmer: K = f.extend_left(l_ext);
                assert!(kmer_set.contains(&(ext_kmer.min_rc())));
            }

            for r_ext in exts.get(Dir::Right) {
                let f: K = seq.last_kmer();
                let ext_kmer: K = f.extend_right(r_ext);
                assert!(kmer_set.contains(&(ext_kmer.min_rc())));
            }
        }

        assert_eq!(kmer_set, all_contig_set);
    }

    #[test]
    fn simple_tip_clean() {
        let contigs = vec![random_dna(200), random_dna(200)];
        test_tip_cleaning::<IntKmer<u64>>(contigs, false);
    }

    // Take some input contig, which likely form a complicated graph,
    // and test the kmer, bsp, sedge and edge construction machinery
    fn test_tip_cleaning<K: Kmer + Sync + Send>(contigs: Vec<Vec<u8>>, stranded: bool) {
        let mut clean_seqs = Vec::new();
        let mut all_seqs = Vec::new();

        // Generate 5x coverage of the main sequences & add some tips
        for c in contigs {
            if c.len() < K::k() * 3 {
                continue;
            }

            for _i in 0..5 {
                clean_seqs.push((DnaBytes(c.clone()), Exts::empty(), ()));
                all_seqs.push((DnaBytes(c.clone()), Exts::empty(), ()));
            }

            let junk = random_dna(5);
            let mut err_ctg = c.clone();
            let l = err_ctg.len();
            err_ctg.truncate(l / 2);
            err_ctg.extend(junk);
            all_seqs.push((DnaBytes(err_ctg.clone()), Exts::empty(), ()));
            all_seqs.push((DnaBytes(err_ctg.clone()), Exts::empty(), ()));
        }

        // Assemble w/o tips
        let (valid_kmers_clean, _): (BoomHashMap2<K, Exts, u16>, _) = filter::filter_kmers(
            &clean_seqs,
            &Box::new(filter::CountFilter::new(2)),
            stranded,
            false,
            4,
        );
        let spec = SimpleCompress::new(|d1: u16, d2: &u16| d1 + d2);
        let graph = compress_kmers_with_hash(stranded, &spec, &valid_kmers_clean);
        let graph1 = graph.finish();
        graph1.print();

        // Assemble w/ tips
        let (valid_kmers_errs, _): (BoomHashMap2<K, Exts, u16>, _) = filter::filter_kmers(
            &all_seqs,
            &Box::new(filter::CountFilter::new(2)),
            stranded,
            false,
            4,
        );
        let spec = SimpleCompress::new(|d1: u16, d2: &u16| d1 + d2);
        let graph = compress_kmers_with_hash(stranded, &spec, &valid_kmers_errs);
        let graph2 = graph.finish();
        graph2.print();

        // Now try to clean the tips.
        let cleaner = CleanGraph::new(|node| node.len() < K::k() * 2);
        let nodes_to_censor = cleaner.find_bad_nodes(&graph2);

        println!("censor: {:?}", nodes_to_censor);
        let spec = SimpleCompress::new(|d1: u16, d2: &u16| d1 + d2);
        let fixed = compress_graph(stranded, &spec, graph2, Some(nodes_to_censor));
        fixed.print();
    }
}
