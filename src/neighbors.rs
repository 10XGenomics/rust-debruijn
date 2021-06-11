use crate::{Kmer, MerImmut};

/// Generate all Hamming distance 1 neighbors of a kmer
pub struct KmerOneHammingIter<K>
where
    K: Kmer,
{
    source: K,       // Original kmer from which we need to generate values
    position: usize, // Index into kmer where last base was mutated
    char: u8,        // The last base which was used
}

impl<K> KmerOneHammingIter<K>
where
    K: Kmer,
{
    /// Create an iterator over all Hamming distance=1 neighbors of `kmer`.
    pub fn new(kmer: K) -> Self {
        KmerOneHammingIter {
            source: kmer,
            position: 0,
            char: 0,
        }
    }
}

impl<K> Iterator for KmerOneHammingIter<K>
where
    K: Kmer,
{
    type Item = K;

    fn next(&mut self) -> Option<Self::Item> {
        if self.position >= self.source.len() {
            return None;
        }
        let base_at_pos = self.source.get(self.position);

        if self.char >= 4 {
            self.position += 1;
            self.char = 0;
            self.next()
        } else if base_at_pos == self.char {
            self.char += 1;
            self.next()
        } else {
            let next_sseq = self.source.set(self.position, self.char);
            self.char += 1;
            Some(next_sseq)
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::kmer::Kmer12;
    use crate::Kmer;

    #[test]
    fn test_hd1() {
        let k1 = Kmer12::from_ascii(b"AAAGGGTTTCCC");

        let mut hd1: Vec<_> = KmerOneHammingIter::new(k1).collect();
        assert_eq!(hd1.len(), 12 * 3);

        hd1.sort();
        hd1.dedup();
        assert_eq!(hd1.len(), 12 * 3);

        for k in &hd1 {
            assert_eq!(k1.hamming_dist(*k), 1);
        }
    }
}
