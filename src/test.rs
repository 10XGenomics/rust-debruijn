//! Generate random genomes (with lots of re-used sustrings), reassemble them, and check sanity

use Kmer;
use Vmer;
use complement;

use std::cmp::{min, max};
use rand::{self, Rng};
use rand::distributions::{IndependentSample, Gamma, Range};


pub fn random_base() -> u8 {
    let mut r = rand::thread_rng();
    (r.next_u64() % 4) as u8
}

// Generate uniformly random DNA sequences
pub fn random_dna(len: usize) -> Vec<u8> {
    let mut r = rand::thread_rng();
    let mut dna = Vec::new();
    for _ in 0..len {
        let b = (r.next_u64() % 4) as u8;
        dna.push(b);
    }

    dna
}

pub fn edit_dna<R: Rng>(seq: &mut Vec<u8>, p: f64, r: &mut R) {
    for b in seq.iter_mut() {
        if r.next_f64() < p {
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

pub fn random_vmer<K: Kmer, V: Vmer<K>>() -> V {
    let mut r = rand::thread_rng();
    let len = Range::new(K::k(), min(200, V::max_len())).ind_sample(&mut r);
    let mut lmer = V::new(len);

    for pos in 0..len {
        let b = (r.next_u64() % 4) as u8;
        lmer.set_mut(pos, b);
    }
    lmer
}

fn simple_contigs() -> Vec<Vec<u8>> {
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

    let nchunks = max(5, gamma_dist.ind_sample(&mut rng) as u32);
    let chunk_sample = Range::new(0, nchunks);

    let length_dist = Gamma::new(1.5, 200.0);


    let mut chunks: Vec<Vec<u8>> = Vec::new();
    for _ in 0..nchunks {
        let len = max(10, length_dist.ind_sample(&mut rng) as usize);
        let seq = random_dna(len);
        chunks.push(seq);
    }

    // Now make a bunch of chromosomes by pasting together chunks
    let nchrom = max(4, gamma_dist.ind_sample(&mut rng) as u32);

    let mut chroms = Vec::new();
    for _ in 0..nchrom {
        let chrom_chunks = max(4, gamma_dist.ind_sample(&mut rng) as u32);

        let mut chrom_seq = Vec::new();
        for _ in 0..chrom_chunks {
            let chunk_idx = chunk_sample.ind_sample(&mut rng) as usize;
            chrom_seq.extend(chunks[chunk_idx].clone());
        }
        chroms.push(chrom_seq);
    }

    chroms
}

#[cfg(test)]
mod tests {

    use {Kmer, Dir, Exts};
    use clean_graph::CleanGraph;
    use std::collections::HashSet;
    use paths::{BaseGraph, DebruijnGraph};
    use compression::{SimpleCompress, PathCompression, Simplify};
    use std::iter::FromIterator;
    use std::hash::{Hash, SipHasher, Hasher};
    use std::marker::PhantomData;
    use DnaBytes;

    use std::ops::Sub;
    use msp;
    use kmer::IntKmer;    
    use dna_string::DnaString;
    use filter;

    use super::*;

    #[test]
    fn simple_kmer_compress() {
        let contigs = simple_contigs();
        reassemble_contigs::<IntKmer<u64>, DnaString>(contigs, false);
    }

    #[test]
    fn complex_kmer_compress() {
        for _ in 0..5 {
            let contigs = random_contigs();
            reassemble_contigs::<IntKmer<u64>, DnaString>(contigs, false);
        }
    }

    #[test]
    fn simple_path_compress() {
        let contigs = simple_contigs();
        simplify_from_kmers::<IntKmer<u64>>(contigs, false);
    }

    #[test]
    fn complex_path_compress() {
        for _ in 0..5 {
            let contigs = random_contigs();
            simplify_from_kmers::<IntKmer<u64>>(contigs, false);
        }
    }

    fn simplify_from_kmers<K: Kmer>(mut contigs: Vec<Vec<u8>>, stranded: bool) {
        
        use DnaBytes;
        let seqs = contigs.drain(..).map(|x| (DnaBytes(x), Exts::empty(), ())).collect();
        let (valid_kmers, _): (Vec<(K, (Exts, _))>, _) = 
            filter::filter_kmers(&seqs, filter::CountFilter::new(1), stranded);

        // Create a DBG with one node per input kmer
        let mut base_graph: BaseGraph<K, u16> = BaseGraph::new(stranded);

        for (kmer, (exts, _)) in valid_kmers.clone() {
            base_graph.add(kmer.iter(), exts, 1);
        }
        let uncompressed_dbg = base_graph.finish();

        // Canonicalize the graph with 
        let spec = SimpleCompress::new(|mut d1: u16, d2: &u16| { d1 + d2 });
        let simp_dbg = Simplify::simplify(uncompressed_dbg, None, stranded, spec);

        let is_cmp = simp_dbg.is_compressed();
        if is_cmp.is_some() {
            println!("not compressed: nodes: {:?}", is_cmp);
            simp_dbg.print();
        }

        assert!(simp_dbg.is_compressed() == None);
    }

    // Take some input contig, which likely form a complicated graph,
    // and test the kmer, bsp, sedge and edge construction machinery
    fn reassemble_contigs<K:Kmer + Copy, V:Vmer<K>>(contigs: Vec<Vec<u8>>, stranded: bool) {
        let ctg_lens: Vec<_> = contigs.iter().map(|c| c.len()).collect();
        println!("Reassembling contig sizes: {:?}", ctg_lens);

        // kmer vector
        let mut kmers = Vec::new();
        for c in contigs.iter() {
            let mut _kmers = K::kmers_from_string(c.as_slice());
            kmers.extend(_kmers.iter().map(|k| k.min_rc()));
        }

        // True kmer set
        let mut kmer_set = HashSet::new();
        kmer_set.extend(kmers.iter());

        // Bsps of kmers
        let P = 6;

        let mut seqs: Vec<(V, Exts, u8)> = Vec::new();
        let permutation = (0..1 << (2 * P)).collect();
        for c in contigs.iter() {            
            let msps = msp::msp_sequence::<K, V>(P, c.as_slice(), Some(&permutation));
            seqs.extend(msps.clone().into_iter().map(|(_, e, v)| (v, e, 0u8)));
            seqs.extend(msps.into_iter().map(|(_, e, v)| (v, e, 1u8)));
        }

        // kmer set from bsps
        let mut msp_kmers = HashSet::new();
        for &(ref v, _, _) in seqs.iter() {
            for k in v.iter_kmers() {
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
        let (valid_kmers, _) = filter::filter_kmers(&seqs, filter::CountFilter::new(2), stranded);
        let mut process_kmer_set = HashSet::new();
        for k in valid_kmers.iter().map(|x| x.0) {
            process_kmer_set.insert(k);
        }
        assert_eq!(process_kmer_set, kmer_set);

        // Every kmer should be reachable as the extension of a kmer.
        // No new kmers should be reachable
        let mut extension_kmer_set: HashSet<K> = HashSet::new();
        for &(kmer, (exts, _)) in valid_kmers.iter() {
            for e in kmer.get_extensions(exts, Dir::Left) {
                extension_kmer_set.insert(e.min_rc());
            }

            for e in kmer.get_extensions(exts, Dir::Right) {
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

        let spec = SimpleCompress::new(|mut d1: u16, d2: &u16| { d1 + d2 });
        let pc: PathCompression<K,_,_> = PathCompression::new(stranded, spec);

        // Now generate the lines for these kmers
        let graph = pc.build_nodes(&valid_kmers);

        // Check that all the lines have valid kmers,
        // and have extensions into other valid kmers
        for seq_id in 0 .. graph.len() {
            let seq = graph.sequences.get(seq_id);
            let exts = graph.exts[seq_id];
            
            let seq_set = HashSet::from_iter(seq.iter_kmers().map(|km: K| km.min_rc()));
            assert!(kmer_set.is_superset(&seq_set));

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
    }


    #[test]
    fn simple_tip_clean() {
        let contigs = vec![random_dna(200), random_dna(200)];
        test_tip_cleaning::<IntKmer<u64>>(contigs, false);
    }

    // Take some input contig, which likely form a complicated graph,
    // and test the kmer, bsp, sedge and edge construction machinery
    fn test_tip_cleaning<K:Kmer>(contigs: Vec<Vec<u8>>, stranded: bool) {
        
        let mut clean_seqs = Vec::new();
        let mut all_seqs = Vec::new();

        // Generate 5x coverage of the main sequences & add some tips
        for c in contigs {
            if c.len() < K::k() * 3 { continue }

            for i in 0 .. 5 {
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
        let (valid_kmers_clean, _) = filter::filter_kmers(&clean_seqs, filter::CountFilter::new(2), stranded);
        let spec = SimpleCompress::new(|mut d1: u16, d2: &u16| { d1 + d2 });
        let pc: PathCompression<K,_,_> = PathCompression::new(stranded, spec);
        let graph1 = pc.build_nodes(&valid_kmers_clean).finish();
        graph1.print();

        // Assemble w/ tips
        let (valid_kmers_errs, _) = filter::filter_kmers(&all_seqs, filter::CountFilter::new(2), stranded);
        let spec = SimpleCompress::new(|mut d1: u16, d2: &u16| { d1 + d2 });
        let pc: PathCompression<K,_,_> = PathCompression::new(stranded, spec);
        let graph2 = pc.build_nodes(&valid_kmers_errs).finish();
        graph2.print();

        // Now try to clean the tips.
        let cleaner = CleanGraph::new(|node| node.len() < K::k() * 2);
        let nodes_to_censor = cleaner.find_bad_nodes(&graph2);

        println!("censor: {:?}", nodes_to_censor);
        let spec = SimpleCompress::new(|mut d1: u16, d2: &u16| { d1 + d2 });
        let fixed = Simplify::simplify(graph2, Some(nodes_to_censor), stranded, spec);
        fixed.print();
    }
}