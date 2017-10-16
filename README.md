# rust-debruijn
De Bruijn graph construction & path compression libraries.

[Docs](https://10xgenomics.github.io/rust-debruijn/debruijn/index.html)

## Key features
* 2-bit packed fixed-length (Kmer) and variable-length (DnaString) sequence containers
* Statically compiled code paths for different K values
* Ability to track arbitrary auxiliary data through the graph
* Efficient kmer counting & filtering schemes
* Minimum-substring partitioning to shard kmer counting and cDBG construction
* Configurable for stranded and non-stranded input sequence
* Simple tip-removal code


### Kmer structs
- [x] single-integer kmer class that is generic over the integer types (u8, u16, u32, u64, u128)
- [x] use the num crate to get a u128 implementation
- [x] trait hack to support arbitrarily sized kmers (see VarIntKmer) 
- [ ] when there is support for integer type parameters, use that to get intermediate sized kmers.

### Path compression class
- [x] configurable stranded / non-stranded cDBG construction
- [ ] path compression: pass a reducer in to summarize the per-kmer annotation data & associate with graph?
- [ ] path compression:  provide a predicate on pairs of colors that controls where the path compression can start and stop
- [ ] design required for other operations on graphs -- bubble popping, edge trimming, etc.
