# debruijn-rs
Low-memory De Bruijn graph construction & path compression libraries.

## Key attributes
* Very few performance compromises
* Statically compiled code paths for different K
* Ability to track arbitrary auxiliary data through the graph


### Kmer structs
- [x] try to design a single-element kmer class that is generic over the integer types (u8, u16, u32, u64, u128)
- [x] use the num crate to get a u128 implementation
- [ ] when there is support for integer type parameters, use that to get intermediate sized kmers.

- [ ] define a kmer trait with the common kmer operations that can be reused
- [x] update DeBruijn methods to accept the kmer trait
- [ ] build an AnnotatedKmer<K, T> struct that lets you attach generic data to a kmer, and proxies it's kmer implementation through to a kmer
- [ ] AnnotatedKmer idea is probably not workable.


### Path compression class
- [ ] think about stranded / non-stranded analysis: RNA-seq generally knows strand - how to generalize for both cases?
- [ ] improve naming & ergonomics: TempGraph, Edge, VEdge, etc
- [ ] convenience methods for creating DBG without worrying about Lmer sharding
- [ ] path compression: can you pass a reducer in to summarize the per-kmer annotation data & associate with graph?
- [ ] path compression: can you provide a predicate on pairs of colors that controls where the path compression can start and stop
- [ ] design required for other operations on graphs -- bubble popping, edge trimming, etc.


### Pathing
Is there a good (potentially configurable) algorithm / heuristics for increasing the K of a graph by utilizing 1) full read lengths, and 2) read pairs.
