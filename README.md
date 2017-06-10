# rust-debruijn
Low-memory De Bruijn graph construction & path compression libraries.

[Docs](https://10xgenomics.github.io/rust-debruijn/debruijn/index.html)

## Key attributes
* Very few performance compromises
* Statically compiled code paths for different K
* Ability to track arbitrary auxiliary data through the graph


## Todos

### Kmer structs
- [x] try to design a single-element kmer class that is generic over the integer types (u8, u16, u32, u64, u128)
- [x] use the num crate to get a u128 implementation
- [ ] when there is support for integer type parameters, use that to get intermediate sized kmers.

- [x] define a kmer trait with the common kmer operations that can be reused
- [x] update DeBruijn methods to accept the kmer trait
- [ ] build an AnnotatedKmer<K, T> struct that lets you attach generic data to a kmer, and proxies it's kmer implementation through to a kmer
- [ ] AnnotatedKmer idea is probably not workable.


### Path compression class
- [ ] think about stranded / non-stranded analysis: RNA-seq generally knows strand - how to generalize for both cases?
- [ ] path compression: can you pass a reducer in to summarize the per-kmer annotation data & associate with graph?
- [ ] path compression: can you provide a predicate on pairs of colors that controls where the path compression can start and stop
- [ ] design required for other operations on graphs -- bubble popping, edge trimming, etc.
