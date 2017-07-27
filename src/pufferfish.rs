use paths;
use paths::DebruijnGraph;

use boomphf::Mphf

struct SuccinctBitVector {

}

impl SuccintBitVector {
    pub fn new(size: usize) -> SuccintBitVector {
        SuccintBitVector {}
    }

    pub fn set(point: usize) {

    }

    /// Return the position of the nth set bit
    pub fn select(sum: usize) -> usize {

    }

    /// Return the number of bits set prior to position pos
    pub fn rank(pos: usize) -> usize {

    }
}

struct PufferIdx<K: Kmer> {
    graph: DebruijnGraph<K>,
    kmer_mphf: Mphf<K>,
    pos_vec: Vec<usize>
    contig_bounds: SuccinctBitVector,
}

impl<K: Kmer> PufferIdx<K> {

    pub fn new(graph: DebruijnGraph<K>) -> PufferIdx<K> {

        // Generate contig end markers rank/select structure
        let mut contig_end_markers = Vec::new();

        for (start, len) in graph.base.sequences.start.iter().zip(graph.base.sequences.len.iter()) {
            contig_end_markers.push(start + len as usize - 1);
        }
        let contig_bounds = SuccinctBitVector::new(contig_end_marker);


        // Generate a MPHF on all kmers
        let mut kmers = Vec::new()

        for (ind, start) in graph.base.sequences.start.iter().enumerate() {
            for (pos, kmer) in graph.base.sequences.get(ind).iter_kmers() {
                kmers.push(kmer)
            }
        }

        let n_kmers = kmers.len();
        let kmer_mphf = Mphf::new(kmers);

        // Creat Lookup table from h to pos
        let pos_vec = vec![0; n_kmers];

        for (ind, start) in graph.base.sequences.start.iter().enumerate() {
            for (ctg_pos, kmer) in graph.base.sequences.get(ind).iter_kmers() {
                let h = mphf.query(k);
                pos_vec[h] = start + pos
            }
        }

        PufferIdx {
            graph: graph,
            kmer_mphf: kmer_mphf,
            pos_vec: pos_vec,
            contig_bounds: contig_bounds
        }
    }

    /// Attempt to find the given kmer in the Debruijn graph.
    /// Return the contig id and contig position if found.
    pub fn find_kmer(&self, k: K) -> Option<(u32, u32)> {

        let h = self.kmer_table.query(k);
        let pos = self.pos_vec[h];

        let cseq = self.graph.base.sequences.sequence;
        if cseq.get_kmer(pos) == k {
            // Which contig are we in?
            let ctg_id = self.contig_bounds.select(pos);
            let ctg_start = self.contig_bounds.rank(ctg_id);
            let ctg_pos = pos - ctg_start;
            Some((ctg_id as u32, ctg_pos as u32))
        }

        None
    }


}