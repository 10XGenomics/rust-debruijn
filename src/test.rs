//! Generate random genomes (with lots of re-used sustrings), reassemble them, and check sanity

use Mer;
use Kmer;
use dna_string::{DnaString, DnaStringSlice};
use vmer::Vmer;
use Dir;
use Exts;

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

/*
pub fn random_pmer() -> Pmer {
    let mut r = rand::thread_rng();
    let p = 8;
    let mut kmer = Pmer::empty(p);
    for pos in 0..p {
        let b = (r.next_u64() % 4) as u8;
        kmer = kmer.set(pos, b);
    }
    kmer
}
*/

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

    use {Kmer, Dir, Bsp, Exts, complement};
    use std::collections::HashSet;
    use path;
    use std::iter::FromIterator;
    use std::hash::{Hash, SipHasher, Hasher};
    use fx::FxHashMap;
    use std::fs::{File, remove_file};
    use std::io::{BufReader, BufWriter};
    use std::ops::Sub;
    //use utils;

    use super::*;

    fn hash<T: Hash>(t: &T) -> u64 {
        let mut s = SipHasher::new();
        t.hash(&mut s);
        s.finish()
    }

    #[test]
    fn test_line_construction() {
        for _ in 0..5 {
            let contigs = random_contigs();
            reassemble_contigs(contigs);
        }
    }

    #[test]
    fn test_sn_graph_round_trip()
    {
        for n2 in vec![100, 200, 500, 1000] {
            check_sn_graph_round_trip(5, n2);
        }
    }


    fn check_sn_graph_round_trip(n1: usize, n2: usize) {
        let mut g = debruijn::TempGraph::new(false);

        println!("hello");
        for i in n1..n2 {
            let s : Vec<u8> = random_dna(i);
            let v: Option<Vec<u32>> = None;
            g.add(s, Exts::empty(), v);
        }

        let test_file = "test_sn_graph.hbv";

        {
            let mut _wtr = File::create(test_file).unwrap();
            let mut wtr = BufWriter::new(_wtr);
            g.write_to_sn_format(&mut wtr);
        }
        
        let round_trip = {
             let _r = File::open(test_file).unwrap();
             let mut r = BufReader::new(_r);
             debruijn::TempGraph::read_from_sn_format(&mut r)
        };

        remove_file(test_file).unwrap();
        assert_eq!(g.start_pos, round_trip.start_pos);
        assert_eq!(g.edge_len, round_trip.edge_len);
        assert_eq!(g.sequence, round_trip.sequence);
    }


    #[test]
    fn simple_line_construction() {
        let p1 = random_dna(20);
        let p2 = random_dna(20);

        let pc = random_dna(80);

        let p3 = random_dna(20);
        let p4 = random_dna(20);



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
        c3.extend(random_dna(40));
        
        // Stick a palindrome in the middle
        let palindrome1 = random_dna(30);
        let mut palindrome2 = palindrome1.clone();
        palindrome2.reverse();
        for i in 0..palindrome2.len() {
            palindrome2[i] = complement(palindrome2[i]);
        }

        c3.extend(palindrome1);
        c3.extend(palindrome2);
        c3.extend(random_dna(40));

        let contigs = vec![c1, c2, c3];
        reassemble_contigs(contigs);
    }

    // Take some input contig, which likely form a complicated graph,
    // and test the kmer, bsp, sedge and edge construction machinery
    fn reassemble_contigs<K:Kmer>(contigs: Vec<Vec<u8>>) {
        let ctg_lens: Vec<_> = contigs.iter().map(|c| c.len()).collect();
        println!("Reassembling contig sizes: {:?}", ctg_lens);

        // kmer vector
        let mut kmers = Vec::new();
        for c in contigs.iter() {
            let mut _kmers = Kmer::kmers_from_string(&c[..]);
            kmers.extend(_kmers.iter().map(|k| k.min_rc()));
        }

        // True kmer set
        let mut kmer_set = HashSet::new();
        kmer_set.extend(kmers);

        // Bsps of kmers
        let P = 6;

        let mut bsps = Vec::new();
        let permutation = (0..1 << (2 * P)).collect();
        for (idx, c) in contigs.iter().enumerate() {            
            let _bsps = Bsp::msp_read(kmer::K, P, (idx*2+1) as u32, 0, &c[..], &permutation);
            bsps.extend(_bsps.into_iter().map(|(_, b)| b));

            let _bsps = Bsp::msp_read(kmer::K, P, (idx*2 + 2) as u32, 0, &c[..], &permutation);
            bsps.extend(_bsps.into_iter().map(|(_, b)| b));
        }

        // kmer set from bsps
        let mut bsp_kmers = HashSet::new();
        for b in bsps.iter() {
            let _kmers = b.kmers();
            for k in _kmers {
                bsp_kmers.insert(k.min_rc());
            }
        }

        // Raw kmers and BSP kmers match
        assert!(kmer_set == bsp_kmers);

        // Check the correctness of the process_kmer_shard kmer filtering function
        let (valid_kmers, all_kmers) = utils::process_kmer_shard(&bsps, 2);
        let mut process_kmer_set = HashSet::new();
        for k in all_kmers.keys() {
            process_kmer_set.insert(*k);
        }
        assert_eq!(process_kmer_set, kmer_set);

        // Full set of barcodes from valid kmers should equal 1..contigs.len()*2+1  (0 barcodes are not counted)
        let mut obs_bcs = HashSet::new();
        for (_, bcs) in valid_kmers {
            for bc in bcs { obs_bcs.insert(bc); }
        }

        let mut true_bcs = HashSet::new();
        for bc in 1 .. contigs.len() * 2 + 1 {
            true_bcs.insert(bc as u32);
        }

        assert_eq!(obs_bcs, true_bcs);


        // Now generate and check the kmer extensions
        let mut kmer_hash_map = FxHashMap();
        for k in kmer_set.iter() {
            kmer_hash_map.insert(*k, 0);
        }

        let kmer_exts = debruijn::kmer_extensions(&kmer_hash_map, &kmer_hash_map, &bsps);

        // Every kmer should be reachable as the extension of a kmer.
        // No new kmers should be reachable
        let mut extension_kmer_set: HashSet<K> = HashSet::new();
        for (kmer, exts) in kmer_exts.iter() {
            for e in kmer.get_extensions(*exts, Dir::Left) {
                extension_kmer_set.insert(e.min_rc());
            }

            for e in kmer.get_extensions(*exts, Dir::Right) {
                extension_kmer_set.insert(e.min_rc());
            }
        }

        // Kmer extensions from BSP match raw kmers
        assert!(kmer_set == extension_kmer_set);

        let mut ext_vec: Vec<(K, Exts)> = kmer_exts.clone().into_iter().collect();
        ext_vec.sort();
        let ext_vec_hash = hash(&ext_vec);
        println!("ext_vec hash: {}", ext_vec_hash);



        // Now generate the lines for these kmers
        let sedges = debruijn::build_sedges(&kmer_exts);

        // println!("--- Sedges ---");
        // let sdb = debruijn::SedgeDb::new(sedges.clone());

        // debruijn::SedgeDb::new(sedges.clone()).print();

        // println!("n edges: {:?}", sedges.len());
        // for (buf, exts) in sedges.clone() {
        //    println!("sedge: {:?}", buf);
        //    println!("exts: {:?}", exts);
        // }

        // Check that all the lines have valid kmers,
        // and have extensions into other valid kmers
        for (seq, exts) in sedges.clone() {
            let kmers = seq.kmers();
            let sedge_set = HashSet::from_iter(kmers.iter().map(|km| km.min_rc()));
            assert!(kmer_set.is_superset(&sedge_set));

            for l_ext in exts.get(Dir::Left) {
                let ext_kmer = kmers.first().expect("kmer").extend_left(l_ext);
                assert!(kmer_set.contains(&(ext_kmer.min_rc())));
            }

            for r_ext in exts.get(Dir::Right) {
                let ext_kmer = kmers.last().expect("kmer").extend_right(r_ext);
                assert!(kmer_set.contains(&(ext_kmer.min_rc())));
            }
        }

        // Now upgrade the sedges into complete edges
        let (temp_graph, _) = debruijn::build_edges(sedges, None);

        // let mut tot_length : usize = 0;
        // for l in temp_graph.edge_len.iter() { tot_length += *l as usize };
        // println!("Graph Edges: {}, Total Edge size: {}", temp_graph.edge_len.len(), tot_length);

        // for e in temp_graph.iter()
        // {
        //    println!("{:?}", e);
        // }

        // Check that all the edges and their extension contain valid kmers
        for e in temp_graph.iter() {
            let kmers = e.sequence.kmers();
            let edge_set = HashSet::from_iter(kmers.iter().map(|km| km.min_rc()));

            if !kmer_set.is_superset(&edge_set) {
                let invented = edge_set.sub(&kmer_set);
                println!("Invented kmers: {:?}", invented);
            }

            assert!(kmer_set.is_superset(&edge_set));

            for l_ext in e.exts.get(Dir::Left) {
                let curr = kmers.first().expect("kmer");
                let ext_kmer = curr.extend_left(l_ext);
                let ext_min = ext_kmer.min_rc();

                if !kmer_set.contains(&(ext_kmer.min_rc())) {
                    println!("edge: {:?}", e);
                    println!("end kmer: {:?}", curr);
                    println!("ext kmer: {:?}", ext_kmer);
                    println!("ext rc: {:?}", ext_min);
                }
                assert!(kmer_set.contains(&(ext_kmer.min_rc())));
            }

            for r_ext in e.exts.get(Dir::Right) {
                let curr = kmers.last().expect("kmer");
                let ext_kmer = curr.extend_right(r_ext);
                let ext_min = ext_kmer.min_rc();

                if !kmer_set.contains(&(ext_kmer.min_rc())) {
                    println!("edge: {:?}", e);
                    println!("end kmer: {:?}", curr);
                    println!("ext kmer: {:?}", ext_kmer);
                    println!("ext rc: {:?}", ext_min);
                }

                assert!(kmer_set.contains(&(ext_kmer.min_rc())));
            }
        }
    }
}