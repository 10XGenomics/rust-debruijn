# rust-debruijn
De Bruijn graph construction & path compression libraries.

[Docs](https://docs.rs/debruijn/)

## Key features
* 2-bit packed fixed-length (Kmer) and variable-length (DnaString) sequence containers
* Statically compiled code paths for different K values
* Ability to track arbitrary auxiliary data through the DeBruijn graph
* Customizable kmer counting & filtering schemes supporting a variety of use cases
* DeBruijn graph compression
* Minimum-substring partitioning to shard kmers for memory efficient counting and DeBruijn graph compression
* Configurable for stranded and non-stranded input sequence
* Extensive unit test suite
* In production use in Supernova, Long Ranger, Cell Ranger, and Cell Ranger VDJ pipelines from 10x Genomics.
