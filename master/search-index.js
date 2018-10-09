var N = null;var searchIndex = {};
searchIndex["debruijn"]={"doc":"debruijn-rs: a De Bruijn graph for DNA seqeunces in Rust. This library provides tools for efficient construction DeBruijn graphs from DNA sequences, tracking arbitrary metadata associated with kmers in the graph, and performing path-compression of unbranched graph paths to improve speed and reduce memory consumption. All the data structures in debruijn-rs are specialized to the 4 base DNA alphabet, and use 2-bit packed encoding of base-pairs into integer types, and efficient methods for reverse complement, enumerating kmers from longer sequences, and transfering data between sequences.","items":[[3,"MerIter","debruijn","Iterator over bases of a DNA sequence (bases will be unpacked into bytes).",N,N],[3,"DnaBytes","","A newtype wrapper around a `Vec<u8>` with implementations ",N,N],[12,"0","","",0,N],[3,"DnaSlice","","A newtype wrapper around a `&[u8]` with implementations ",N,N],[12,"0","","",1,N],[3,"Exts","","Store single-base extensions for a DNA Debruijn graph.",N,N],[12,"val","","",2,N],[3,"KmerIter","","Iterate over the `Kmer`s of a DNA sequence efficiently",N,N],[3,"KmerExtsIter","","Iterate over the `(Kmer, Exts)` tuples of a sequence and it's extensions efficiently",N,N],[4,"Dir","","Direction of motion in a DeBruijn graph",N,N],[13,"Left","","",3,N],[13,"Right","","",3,N],[5,"bits_to_ascii","","Convert a 2-bit representation of a base to a char",N,[[["u8"]],["u8"]]],[5,"base_to_bits","","Convert an ASCII-encoded DNA base to a 2-bit representation",N,[[["u8"]],["u8"]]],[5,"dna_only_base_to_bits","","",N,[[["u8"]],["option",["u8"]]]],[5,"is_valid_base","","Convert an ASCII-encoded DNA base to a 2-bit representation",N,[[["u8"]],["bool"]]],[5,"bits_to_base","","Convert a 2-bit representation of a base to a char",N,[[["u8"]],["char"]]],[5,"complement","","The complement of a 2-bit encoded base",N,[[["u8"]],["u8"]]],[0,"kmer","","Represent kmers with statically know length in compact integer types",N,N],[3,"IntKmer","debruijn::kmer","A Kmer sequence with a statically know K. K will fill the underlying integer type.",N,N],[12,"storage","","",4,N],[3,"VarIntKmer","","A fixed-length Kmer sequence that may not fill the bits of T",N,N],[12,"storage","","",5,N],[12,"phantom","","",5,N],[3,"K48","","Marker struct for generating K=48 Kmers",N,N],[3,"K40","","Marker trait for generating K=40 Kmers",N,N],[3,"K30","","Marker trait for generating K=40 Kmers",N,N],[3,"K24","","Marker trait for generating K=24 Kmers",N,N],[3,"K20","","Marker trait for generating K=20 Kmers",N,N],[3,"K14","","Marker trait for generating K=14 Kmers",N,N],[3,"K6","","Marker trait for generating K=6 Kmers",N,N],[3,"K5","","Marker trait for generating K=6 Kmers",N,N],[3,"K4","","Marker trait for generating K=6 Kmers",N,N],[3,"K3","","Marker trait for generating K=6 Kmers",N,N],[3,"K2","","Marker trait for generating K=6 Kmers",N,N],[6,"Kmer64","","64-base kmer, backed by a single u128",N,N],[6,"Kmer48","","48-base kmer, backed by a single u128",N,N],[6,"Kmer40","","40-base kmer, backed by a single u128",N,N],[6,"Kmer32","","32-base kmer, backed by a single u64",N,N],[6,"Kmer30","","30-base kmer, backed by a single u64",N,N],[6,"Kmer24","","24-base kmer, backed by a single u64",N,N],[6,"Kmer20","","20-base kmer, backed by a single u64",N,N],[6,"Kmer16","","16-base kmer, backed by a single u32",N,N],[6,"Kmer14","","14-base kmer, backed by a single u32",N,N],[6,"Kmer8","","16-base kmer, backed by a single u16",N,N],[6,"Kmer6","","",N,N],[6,"Kmer5","","",N,N],[6,"Kmer4","","",N,N],[6,"Kmer3","","",N,N],[6,"Kmer2","","",N,N],[8,"IntHelp","","Trait for specialized integer operations used in DeBruijn Graph",N,N],[10,"reverse_by_twos","","Reverse the order of 2-bit units of the integer",6,[[["self"]],["self"]]],[8,"KmerSize","","Helper trait for declaring the K value of a Kmer. Will be removed when const generics are available",N,N],[10,"K","","",7,[[],["usize"]]],[11,"clone","","",4,[[["self"]],["intkmer"]]],[11,"eq","","",4,[[["self"],["intkmer"]],["bool"]]],[11,"ne","","",4,[[["self"],["intkmer"]],["bool"]]],[11,"partial_cmp","","",4,[[["self"],["intkmer"]],["option",["ordering"]]]],[11,"lt","","",4,[[["self"],["intkmer"]],["bool"]]],[11,"le","","",4,[[["self"],["intkmer"]],["bool"]]],[11,"gt","","",4,[[["self"],["intkmer"]],["bool"]]],[11,"ge","","",4,[[["self"],["intkmer"]],["bool"]]],[11,"cmp","","",4,[[["self"],["intkmer"]],["ordering"]]],[11,"hash","","",4,N],[11,"top_mask","","",4,[[["usize"]],["t"]]],[11,"bottom_mask","","",4,[[["usize"]],["t"]]],[11,"len","","",4,[[["self"]],["usize"]]],[11,"get","","Get the letter at the given position.",4,[[["self"],["usize"]],["u8"]]],[11,"set_mut","","",4,[[["self"],["usize"],["u8"]]]],[11,"set_slice_mut","","Set a slice of bases in the kmer, using the packed representation in value. Sets n_bases, starting at pos. Bases must always be packed into the upper-most bits of the value.",4,[[["self"],["usize"],["usize"],["u64"]]]],[11,"rc","","Return the reverse complement of this kmer",4,[[["self"]],["self"]]],[11,"empty","","",4,[[],["self"]]],[11,"k","","",4,[[],["usize"]]],[11,"from_u64","","",4,[[["u64"]],["intkmer"]]],[11,"to_u64","","",4,[[["self"]],["u64"]]],[11,"extend_left","","Shift the base v into the left end of the kmer",4,[[["self"],["u8"]],["self"]]],[11,"extend_right","","",4,[[["self"],["u8"]],["self"]]],[11,"fmt","","",4,[[["self"],["formatter"]],["result"]]],[11,"clone","","",5,[[["self"]],["varintkmer"]]],[11,"eq","","",5,[[["self"],["varintkmer"]],["bool"]]],[11,"ne","","",5,[[["self"],["varintkmer"]],["bool"]]],[11,"partial_cmp","","",5,[[["self"],["varintkmer"]],["option",["ordering"]]]],[11,"lt","","",5,[[["self"],["varintkmer"]],["bool"]]],[11,"le","","",5,[[["self"],["varintkmer"]],["bool"]]],[11,"gt","","",5,[[["self"],["varintkmer"]],["bool"]]],[11,"ge","","",5,[[["self"],["varintkmer"]],["bool"]]],[11,"cmp","","",5,[[["self"],["varintkmer"]],["ordering"]]],[11,"hash","","",5,N],[11,"empty","","",5,[[],["self"]]],[11,"k","","",5,[[],["usize"]]],[11,"to_u64","","",5,[[["self"]],["u64"]]],[11,"from_u64","","",5,[[["u64"]],["self"]]],[11,"extend_left","","Shift the base v into the left end of the kmer",5,[[["self"],["u8"]],["self"]]],[11,"extend_right","","",5,[[["self"],["u8"]],["self"]]],[11,"top_mask","","",5,[[["usize"]],["t"]]],[11,"bottom_mask","","",5,[[["usize"]],["t"]]],[11,"len","","",5,[[["self"]],["usize"]]],[11,"get","","Get the letter at the given position.",5,[[["self"],["usize"]],["u8"]]],[11,"set_mut","","",5,[[["self"],["usize"],["u8"]]]],[11,"set_slice_mut","","Set a slice of bases in the kmer, using the packed representation in value. Sets n_bases, starting at pos. Incoming bases must always be packed into the upper-most bits of the value.",5,[[["self"],["usize"],["usize"],["u64"]]]],[11,"rc","","Return the reverse complement of this kmer",5,[[["self"]],["self"]]],[11,"fmt","","",5,[[["self"],["formatter"]],["result"]]],[11,"fmt","","",8,[[["self"],["formatter"]],["result"]]],[11,"hash","","",8,N],[11,"clone","","",8,[[["self"]],["k48"]]],[11,"cmp","","",8,[[["self"],["k48"]],["ordering"]]],[11,"partial_cmp","","",8,[[["self"],["k48"]],["option",["ordering"]]]],[11,"eq","","",8,[[["self"],["k48"]],["bool"]]],[11,"K","","",8,[[],["usize"]]],[11,"fmt","","",9,[[["self"],["formatter"]],["result"]]],[11,"hash","","",9,N],[11,"clone","","",9,[[["self"]],["k40"]]],[11,"cmp","","",9,[[["self"],["k40"]],["ordering"]]],[11,"partial_cmp","","",9,[[["self"],["k40"]],["option",["ordering"]]]],[11,"eq","","",9,[[["self"],["k40"]],["bool"]]],[11,"K","","",9,[[],["usize"]]],[11,"fmt","","",10,[[["self"],["formatter"]],["result"]]],[11,"hash","","",10,N],[11,"clone","","",10,[[["self"]],["k30"]]],[11,"cmp","","",10,[[["self"],["k30"]],["ordering"]]],[11,"partial_cmp","","",10,[[["self"],["k30"]],["option",["ordering"]]]],[11,"eq","","",10,[[["self"],["k30"]],["bool"]]],[11,"K","","",10,[[],["usize"]]],[11,"fmt","","",11,[[["self"],["formatter"]],["result"]]],[11,"hash","","",11,N],[11,"clone","","",11,[[["self"]],["k24"]]],[11,"cmp","","",11,[[["self"],["k24"]],["ordering"]]],[11,"partial_cmp","","",11,[[["self"],["k24"]],["option",["ordering"]]]],[11,"eq","","",11,[[["self"],["k24"]],["bool"]]],[11,"K","","",11,[[],["usize"]]],[11,"fmt","","",12,[[["self"],["formatter"]],["result"]]],[11,"hash","","",12,N],[11,"clone","","",12,[[["self"]],["k20"]]],[11,"cmp","","",12,[[["self"],["k20"]],["ordering"]]],[11,"partial_cmp","","",12,[[["self"],["k20"]],["option",["ordering"]]]],[11,"eq","","",12,[[["self"],["k20"]],["bool"]]],[11,"K","","",12,[[],["usize"]]],[11,"fmt","","",13,[[["self"],["formatter"]],["result"]]],[11,"hash","","",13,N],[11,"clone","","",13,[[["self"]],["k14"]]],[11,"cmp","","",13,[[["self"],["k14"]],["ordering"]]],[11,"partial_cmp","","",13,[[["self"],["k14"]],["option",["ordering"]]]],[11,"eq","","",13,[[["self"],["k14"]],["bool"]]],[11,"K","","",13,[[],["usize"]]],[11,"fmt","","",14,[[["self"],["formatter"]],["result"]]],[11,"hash","","",14,N],[11,"clone","","",14,[[["self"]],["k6"]]],[11,"cmp","","",14,[[["self"],["k6"]],["ordering"]]],[11,"partial_cmp","","",14,[[["self"],["k6"]],["option",["ordering"]]]],[11,"eq","","",14,[[["self"],["k6"]],["bool"]]],[11,"K","","",14,[[],["usize"]]],[11,"fmt","","",15,[[["self"],["formatter"]],["result"]]],[11,"hash","","",15,N],[11,"clone","","",15,[[["self"]],["k5"]]],[11,"cmp","","",15,[[["self"],["k5"]],["ordering"]]],[11,"partial_cmp","","",15,[[["self"],["k5"]],["option",["ordering"]]]],[11,"eq","","",15,[[["self"],["k5"]],["bool"]]],[11,"K","","",15,[[],["usize"]]],[11,"fmt","","",16,[[["self"],["formatter"]],["result"]]],[11,"hash","","",16,N],[11,"clone","","",16,[[["self"]],["k4"]]],[11,"cmp","","",16,[[["self"],["k4"]],["ordering"]]],[11,"partial_cmp","","",16,[[["self"],["k4"]],["option",["ordering"]]]],[11,"eq","","",16,[[["self"],["k4"]],["bool"]]],[11,"K","","",16,[[],["usize"]]],[11,"fmt","","",17,[[["self"],["formatter"]],["result"]]],[11,"hash","","",17,N],[11,"clone","","",17,[[["self"]],["k3"]]],[11,"cmp","","",17,[[["self"],["k3"]],["ordering"]]],[11,"partial_cmp","","",17,[[["self"],["k3"]],["option",["ordering"]]]],[11,"eq","","",17,[[["self"],["k3"]],["bool"]]],[11,"K","","",17,[[],["usize"]]],[11,"fmt","","",18,[[["self"],["formatter"]],["result"]]],[11,"hash","","",18,N],[11,"clone","","",18,[[["self"]],["k2"]]],[11,"cmp","","",18,[[["self"],["k2"]],["ordering"]]],[11,"partial_cmp","","",18,[[["self"],["k2"]],["option",["ordering"]]]],[11,"eq","","",18,[[["self"],["k2"]],["bool"]]],[11,"K","","",18,[[],["usize"]]],[0,"dna_string","debruijn","A 2-bit encoding of arbitrary length DNA sequences.",N,N],[3,"DnaString","debruijn::dna_string","A container for sequence of DNA bases.",N,N],[3,"DnaStringIter","","Iterator over values of a DnaStringoded sequence (values will be unpacked into bytes).",N,N],[3,"DnaStringSlice","","An immutable slice into a DnaString",N,N],[12,"dna_string","","",19,N],[12,"start","","",19,N],[12,"length","","",19,N],[12,"is_rc","","",19,N],[3,"PackedDnaStringSet","","Container for many distinct sequences, concatenated into a single DnaString.  Each sequence is accessible by index as a DnaStringSlice.",N,N],[12,"sequence","","",20,N],[12,"start","","",20,N],[12,"length","","",20,N],[11,"cmp","","",21,[[["self"],["dnastring"]],["ordering"]]],[11,"partial_cmp","","",21,[[["self"],["dnastring"]],["option",["ordering"]]]],[11,"lt","","",21,[[["self"],["dnastring"]],["bool"]]],[11,"le","","",21,[[["self"],["dnastring"]],["bool"]]],[11,"gt","","",21,[[["self"],["dnastring"]],["bool"]]],[11,"ge","","",21,[[["self"],["dnastring"]],["bool"]]],[11,"clone","","",21,[[["self"]],["dnastring"]]],[11,"eq","","",21,[[["self"],["dnastring"]],["bool"]]],[11,"ne","","",21,[[["self"],["dnastring"]],["bool"]]],[11,"hash","","",21,N],[11,"len","","",21,[[["self"]],["usize"]]],[11,"get","","Get the value at position `i`.",21,[[["self"],["usize"]],["u8"]]],[11,"set_mut","","Set the value as position `i`.",21,[[["self"],["usize"],["u8"]]]],[11,"set_slice_mut","","",21,[[["self"],["usize"],["usize"],["u64"]]]],[11,"rc","","",21,[[["self"]],["dnastring"]]],[11,"new","","",21,[[["usize"]],["self"]]],[11,"max_len","","",21,[[],["usize"]]],[11,"get_kmer","","Get the kmer starting at position pos",21,[[["self"],["usize"]],["k"]]],[11,"new","","Create an empty DNA string",21,[[],["dnastring"]]],[11,"len","","Length of the sequence",21,[[["self"]],["usize"]]],[11,"empty","","Create a new instance with a given capacity.",21,[[["usize"]],["self"]]],[11,"from_dna_string","","Create a DnaString corresponding to an ACGT-encoded str.",21,[[["str"]],["dnastring"]]],[11,"from_dna_only_string","","Create a DnaString corresponding to an ACGT-encoded str.",21,[[["str"]],["vec",["dnastring"]]]],[11,"from_acgt_bytes","","Create a DnaString from an ASCII ACGT-encoded byte slice. Non ACGT positions will be converted to 'A'",21,N],[11,"from_acgt_bytes_hashn","","Create a DnaString from an ACGT-encoded byte slice, Non ACGT positions will be converted to repeatable random base determined by a hash of the read name and the position within the string.",21,N],[11,"from_bytes","","Create a DnaString from a 0-4 encoded byte slice",21,N],[11,"to_string","","Convert sequence to a String",21,[[["self"]],["string"]]],[11,"to_bytes","","Convert sequence to a Vector of 0-4 encoded bytes",21,[[["self"]],["vec",["u8"]]]],[11,"to_ascii_vec","","Convert sequence to a Vector of ascii-encoded bytes",21,[[["self"]],["vec",["u8"]]]],[11,"push","","Append a 0-4 encoded base.",21,[[["self"],["u8"]]]],[11,"push_bytes","","Push 0-4 encoded bases from a byte array.",21,[[["self"],["vec"],["usize"]]]],[11,"iter","","Iterate over stored values (values will be unpacked into bytes).",21,[[["self"]],["dnastringiter"]]],[11,"clear","","Clear the sequence.",21,[[["self"]]]],[11,"is_empty","","",21,[[["self"]],["bool"]]],[11,"prefix","","Get the length `k` prefix of the DnaString",21,[[["self"],["usize"]],["dnastringslice"]]],[11,"suffix","","Get the length `k` suffix of the DnaString",21,[[["self"],["usize"]],["dnastringslice"]]],[11,"slice","","Get slice containing the interval [`start`, `end`) of `self`",21,[[["self"],["usize"],["usize"]],["dnastringslice"]]],[11,"reverse","","Create a fresh DnaString containing the reverse of `self`",21,[[["self"]],["dnastring"]]],[11,"fmt","","",21,[[["self"],["formatter"]],["result"]]],[11,"next","","",22,[[["self"]],["option",["u8"]]]],[11,"eq","","",19,[[["self"],["dnastringslice"]],["bool"]]],[11,"ne","","",19,[[["self"],["dnastringslice"]],["bool"]]],[11,"clone","","",19,[[["self"]],["dnastringslice"]]],[11,"len","","",19,[[["self"]],["usize"]]],[11,"get","","Get the value at position `i`.",19,[[["self"],["usize"]],["u8"]]],[11,"set_mut","","Set the value as position `i`.",19,[[["self"],["usize"],["u8"]]]],[11,"set_slice_mut","","",19,[[["self"],["usize"],["usize"],["u64"]]]],[11,"rc","","",19,[[["self"]],["dnastringslice"]]],[11,"new","","",19,[[["usize"]],["self"]]],[11,"max_len","","",19,[[],["usize"]]],[11,"get_kmer","","Get the kmer starting at position pos",19,[[["self"],["usize"]],["k"]]],[11,"is_palindrome","","",19,[[["self"]],["bool"]]],[11,"bytes","","",19,[[["self"]],["vec",["u8"]]]],[11,"ascii","","",19,[[["self"]],["vec",["u8"]]]],[11,"to_dna_string","","",19,[[["self"]],["string"]]],[11,"to_string","","",19,[[["self"]],["string"]]],[11,"to_owned","","",19,[[["self"]],["dnastring"]]],[11,"slice","","Get slice containing the interval [`start`, `end`) of `self`",19,[[["self"],["usize"],["usize"]],["dnastringslice"]]],[11,"fmt","","",19,[[["self"],["formatter"]],["result"]]],[11,"fmt","","",20,[[["self"],["formatter"]],["result"]]],[11,"clone","","",20,[[["self"]],["packeddnastringset"]]],[11,"new","","Create an empty `PackedDnaStringSet`",20,[[],["self"]]],[11,"get","","Get a `DnaStringSlice` containing `i`th sequence in the set",20,[[["self"],["usize"]],["dnastringslice"]]],[11,"slice","","Get a `DnaStringSlice` containing `i`th sequence in the set",20,[[["self"],["usize"],["usize"],["usize"]],["dnastringslice"]]],[11,"len","","Number of sequences in the set",20,[[["self"]],["usize"]]],[11,"add","","",20,[[["self"],["s"]]]],[0,"graph","debruijn","Containers for path-compressed De Bruijn graphs",N,N],[3,"BaseGraph","debruijn::graph","A compressed DeBruijn graph carrying auxiliary data on each node of type `D`. This type does not carry the sorted index arrays the allow the graph to be walked efficiently. The `DeBruijnGraph` type wraps this type and add those vectors.",N,N],[12,"sequences","","",23,N],[12,"exts","","",23,N],[12,"data","","",23,N],[12,"stranded","","",23,N],[3,"DebruijnGraph","","A compressed DeBruijn graph carrying auxiliary data on each node of type `D`. The struct carries sorted index arrays the allow the graph to be walked efficiently.",N,N],[12,"base","","",24,N],[3,"NodeIter","","Iterator over nodes in a `DeBruijnGraph`",N,N],[3,"NodeIntoIter","","Iterator over nodes in a `DeBruijnGraph`",N,N],[3,"NodeKmer","","Iterator over nodes in a `DeBruijnGraph`",N,N],[12,"node_id","","",25,N],[3,"NodeKmerIter","","",N,N],[3,"Node","","Unbranched sequence in the DeBruijn graph",N,N],[12,"node_id","","",26,N],[12,"graph","","",26,N],[11,"clone","","",23,[[["self"]],["basegraph"]]],[11,"fmt","","",23,[[["self"],["formatter"]],["result"]]],[11,"new","","",23,[[["bool"]],["self"]]],[11,"len","","",23,[[["self"]],["usize"]]],[11,"combine","","",23,[[["i"]],["self"]]],[11,"add","","",23,[[["self"],["s"],["exts"],["d"]]]],[11,"finish","","",23,[[["self"]],["debruijngraph"]]],[11,"finish_serial","","",23,[[["self"]],["debruijngraph"]]],[11,"fmt","","",24,[[["self"],["formatter"]],["result"]]],[11,"len","","Total number of nodes in the DeBruijn graph",24,[[["self"]],["usize"]]],[11,"get_node","","Get a node given it's `node_id`",24,[[["self"],["usize"]],["node"]]],[11,"iter_nodes","","Return an iterator over all nodes in the graph",24,[[["self"]],["nodeiter"]]],[11,"find_link","","Find a link in the graph, possibly handling a RC switch.",24,[[["self"],["k"],["dir"]],["option"]]],[11,"is_compressed","","Check whether the graph is fully compressed. Return `None` if it's compressed, otherwise return `Some(node1, node2)` representing a pair of node that could be collapsed. Probably only useful for testing.",24,[[["self"]],["option"]]],[11,"fix_exts","","Remove non-existent extensions that may be created due to filtered kmers",24,[[["self"],["option",["bitset"]]]]],[11,"get_valid_exts","","",24,[[["self"],["usize"],["option",["bitset"]]],["exts"]]],[11,"max_path","","Find the highest-scoring, unambiguous path in the graph. Each node get a score given by `score`. Any node where `solid_path(node) == True` are valid paths - paths will be terminated if there are multiple valid paths emanating from a node.",24,[[["self"],["f"],["f2"]],["vec"]]],[11,"sequence_of_path","","Get the sequence of a path through the graph. The path is given as a sequence of node_id integers",24,[[["self"],["i"]],["dnastring"]]],[11,"to_dot","","Write the graph to a dot file.",24,[[["self"],["p"],["f"]]]],[11,"to_gfa","","Write the graph to GFA format",24,[[["self"],["p"]],["result",["error"]]]],[11,"to_gfa_with_tags","","Write the graph to GFA format",24,[[["self"],["p"],["f"]],["result",["error"]]]],[11,"to_json_rest","","",24,[[["self"],["f"],["w"],["option",["value"]]]]],[11,"to_json","","Write the graph to JSON",24,[[["self"],["f"],["w"]]]],[11,"print","","Print a text representation of the graph.",24,[[["self"]]]],[11,"print_with_data","","",24,[[["self"]]]],[11,"to_supernova_bv","","",24,[[["self"],["write"]],["result",["error"]]]],[11,"max_path_beam","","",24,[[["self"],["usize"],["f"],["f2"]],["vec"]]],[11,"next","","",27,[[["self"]],["option",["node"]]]],[11,"next","","",28,[[["self"]],["option"]]],[11,"clone","","",25,[[["self"]],["nodekmer"]]],[11,"into_iter","","",25,N],[11,"next","","",29,[[["self"]],["option"]]],[11,"size_hint","","",29,N],[11,"skip_next","","",29,[[["self"]]]],[11,"len","","Length of the sequence of this node",26,[[["self"]],["usize"]]],[11,"sequence","","Sequence of the node",26,[[["self"]],["dnastringslice"]]],[11,"data","","Reference to auxiliarly data associated with the node",26,[[["self"]],["d"]]],[11,"exts","","Extension bases from this node",26,[[["self"]],["exts"]]],[11,"l_edges","","Edges leaving the left side of the node in the format",26,[[["self"]],["smallvec"]]],[11,"r_edges","","Edges leaving the right side of the node in the format",26,[[["self"]],["smallvec"]]],[11,"edges","","Edges leaving the 'dir' side of the node in the format",26,[[["self"],["dir"]],["smallvec"]]],[11,"fmt","","",26,[[["self"],["formatter"]],["result"]]],[0,"vmer","debruijn","Variable-length DNA strings packed into fixed-size structs.",N,N],[3,"Lmer","debruijn::vmer","Store a variable-length DNA sequence in a packed 2-bit encoding, up 92bp in length The length of the sequence is stored in the lower 8 bits of storage",N,N],[6,"Lmer1","","",N,N],[6,"Lmer2","","",N,N],[6,"Lmer3","","",N,N],[8,"Array","","Types that can be used as the backing store for a SmallVec",N,N],[16,"Item","","",30,N],[10,"new","","",30,[[],["self"]]],[10,"size","","",30,[[],["usize"]]],[10,"as_slice","","",30,N],[10,"as_mut_slice","","",30,N],[11,"hash","","",31,N],[11,"clone","","",31,[[["self"]],["lmer"]]],[11,"eq","","",31,[[["self"],["lmer"]],["bool"]]],[11,"ne","","",31,[[["self"],["lmer"]],["bool"]]],[11,"partial_cmp","","",31,[[["self"],["lmer"]],["option",["ordering"]]]],[11,"lt","","",31,[[["self"],["lmer"]],["bool"]]],[11,"le","","",31,[[["self"],["lmer"]],["bool"]]],[11,"gt","","",31,[[["self"],["lmer"]],["bool"]]],[11,"ge","","",31,[[["self"],["lmer"]],["bool"]]],[11,"cmp","","",31,[[["self"],["lmer"]],["ordering"]]],[11,"len","","The length of the DNA string",31,[[["self"]],["usize"]]],[11,"get","","Get the base at position pos",31,[[["self"],["usize"]],["u8"]]],[11,"set_mut","","Return a new Lmer with position pos set to base val",31,[[["self"],["usize"],["u8"]]]],[11,"set_slice_mut","","",31,[[["self"],["usize"],["usize"],["u64"]]]],[11,"rc","","",31,[[["self"]],["self"]]],[11,"max_len","","",31,[[],["usize"]]],[11,"new","","Initialize an blank Lmer of length len. Will initially represent all A's.",31,[[["usize"]],["lmer"]]],[11,"get_kmer","","Get the kmer starting at position pos",31,[[["self"],["usize"]],["k"]]],[11,"fmt","","",31,[[["self"],["formatter"]],["result"]]],[0,"msp","debruijn","Methods for minimum substring partitioning of a DNA string",N,N],[3,"MspInterval","debruijn::msp","",N,N],[5,"simple_scan","","Determine MSP substrings of seq, for given k and p. Returns a vector of tuples indicating the substrings, and the pmer values: (p-mer value, min p-mer position, start position, end position) permutation is a permutation of the lexicographically-sorted set of all pmers. A permutation of pmers sorted by their inverse frequency in the dataset will give the most even bucketing of MSPs over pmers.",N,N],[5,"msp_sequence","","",N,N],[11,"fmt","","",32,[[["self"],["formatter"]],["result"]]],[11,"new","","",32,[[["u16"],["u32"],["u16"]],["mspinterval"]]],[11,"start","","",32,[[["self"]],["usize"]]],[11,"len","","",32,[[["self"]],["usize"]]],[11,"end","","",32,[[["self"]],["usize"]]],[11,"range","","",32,[[["self"]],["range",["usize"]]]],[11,"bucket","","",32,[[["self"]],["u16"]]],[0,"filter","debruijn","Methods for converting sequences into kmers, filtering observed kmers before De Bruijn graph construction, and summarizing 'color' annotations.",N,N],[3,"CountFilter","debruijn::filter","A simple KmerSummarizer that only accepts kmers that are observed at least a given number of times. The metadata returned about a Kmer is the number of times it was observed, capped at 2^16.",N,N],[3,"CountFilterSet","","A simple KmerSummarizer that only accepts kmers that are observed at least a given number of times. The metadata returned about a Kmer is a vector of the unique data values observed for that kmer.",N,N],[3,"CountFilterEqClass","","",N,N],[5,"filter_kmers","","Process DNA sequences into kmers and determine the set of valid kmers, their extensions, and summarize associated label/'color' data. The input sequences are converted to kmers of type `K`, and like kmers are grouped together. All instances of each kmer, along with their label data are passed to `summarizer`, an implementation of the `KmerSummarizer` which decides if the kmer is 'valid' by an arbitrary predicate of the kmer data, and summarizes the the individual label into a single label data structure for the kmer. Care is taken to keep the memory consumption small. Less than 4G of temporary memory should be allocated to hold intermediate kmers.",N,N],[5,"remove_censored_exts_sharded","","Remove extensions in valid_kmers that point to censored kmers. A censored kmer exists in all_kmers but not valid_kmers. Since the kmer exists in this partition, but was censored, we know that we can delete extensions to it. In sharded kmer processing, we will have extensions to kmers in other shards. We don't know whether these are censored until later, so we retain these extension.",N,[[["bool"],["vec"],["vec"]]]],[5,"remove_censored_exts","","Remove extensions in valid_kmers that point to censored kmers. Use this method in a non-partitioned context when valid_kmers includes all kmers that will ultimately be included in the graph.",N,[[["bool"],["vec"]]]],[6,"EqClassIdType","","",N,N],[8,"KmerSummarizer","","Implement this trait to control how multiple observations of a kmer are carried forward into a DeBruijn graph.",N,N],[10,"summarize","","The input `items` is an iterator over kmer observations. Input observation is a tuple of (kmer, extensions, data). The summarize function inspects the data and returns a tuple indicating: * whether this kmer passes the filtering criteria (e.g. is there a sufficient number of observation) * the accumulated Exts of the kmer * a summary data object of type `DO` that will be used as a color annotation in the DeBruijn graph.",33,N],[11,"new","","Construct a `CountFilter` KmerSummarizer only accepts kmers that are observed at least `min_kmer_obs` times.",34,[[["usize"]],["countfilter"]]],[11,"summarize","","",34,N],[11,"new","","Construct a `CountFilterSet` KmerSummarizer only accepts kmers that are observed at least `min_kmer_obs` times.",35,[[["usize"]],["countfilterset"]]],[11,"summarize","","",35,N],[11,"new","","",36,[[["usize"]],["countfiltereqclass"]]],[11,"get_eq_classes","","",36,[[["self"]],["vec",["vec"]]]],[11,"get_number_of_eq_classes","","",36,[[["self"]],["usize"]]],[11,"fetch_add","","",36,[[["self"]],["usize"]]],[11,"summarize","","",36,N],[0,"compression","debruijn","Create compressed DeBruijn graphs from uncompressed DeBruijn graphs, or a collection of disjoint DeBruijn graphs.",N,N],[3,"SimpleCompress","debruijn::compression","Simple implementation of `CompressionSpec` that lets you provide that data reduction function as a closure",N,N],[3,"ScmapCompress","","",N,N],[5,"compress_graph","","Perform path-compression on a (possibly partially compressed) DeBruijn graph",N,[[["bool"],["s"],["debruijngraph"],["option",["vec"]]],["debruijngraph"]]],[5,"compress_kmers_with_hash","","Take a BoomHash Object and build a compressed DeBruijn graph.",N,[[["bool"],["s"],["boomhashmap2"]],["basegraph"]]],[5,"compress_kmers","","Take a BoomHash Object and build a compressed DeBruijn graph.",N,[[["bool"],["s"],["vec"]],["basegraph"]]],[8,"CompressionSpec","","Customize the path-compression process. Implementing this trait lets the user control how the per-kmer data (of type `D`) is summarized into a per-path summary data (of type `DS`). It also let's the user introduce new breaks into paths by inspecting in the per-kmer data of a proposed with `join_test_kmer` function.",N,N],[10,"reduce","","",37,[[["self"],["d"],["d"]],["d"]]],[10,"join_test","","",37,[[["self"],["d"],["d"]],["bool"]]],[11,"new","","",38,[[["f"]],["simplecompress"]]],[11,"reduce","","",38,[[["self"],["d"],["d"]],["d"]]],[11,"join_test","","",38,[[["self"],["d"],["d"]],["bool"]]],[11,"new","","",39,[[],["scmapcompress"]]],[11,"reduce","","",39,[[["self"],["d"],["d"]],["d"]]],[11,"join_test","","",39,[[["self"],["d"],["d"]],["bool"]]],[0,"clean_graph","debruijn","DeBruijn graph simplification routines. Currently tip-removal is implemented.",N,N],[3,"CleanGraph","debruijn::clean_graph","",N,N],[11,"new","","",40,[[["t1"]],["cleangraph"]]],[11,"find_bad_nodes","","",40,[[["self"],["debruijngraph"]],["vec",["usize"]]]],[0,"test","debruijn","Generate random genomes (with lots of re-used sustrings), reassemble them, and check sanity",N,N],[5,"random_base","debruijn::test","Generate a uniformly random base",N,[[],["u8"]]],[5,"random_dna","","Generate uniformly random DNA sequences",N,[[["usize"]],["vec",["u8"]]]],[5,"edit_dna","","Randomly mutate each base with probability `p`",N,[[["vec"],["f64"],["r"]]]],[5,"random_kmer","","",N,[[],["k"]]],[5,"random_vmer","","",N,[[],["v"]]],[5,"simple_random_contigs","","",N,[[],["vec",["vec"]]]],[5,"random_contigs","","",N,[[],["vec",["vec"]]]],[8,"Mer","debruijn","Trait for interacting with DNA sequences",N,N],[10,"len","","Length of DNA sequence",41,[[["self"]],["usize"]]],[10,"get","","Get 2-bit encoded base at position `pos`",41,[[["self"],["usize"]],["u8"]]],[10,"set_mut","","Set base at `pos` to 2-bit encoded base `val`",41,[[["self"],["usize"],["u8"]]]],[10,"set_slice_mut","","Set `nbases` positions in the sequence, starting at `pos`. Values must  be packed into the upper-most bits of `value`.",41,[[["self"],["usize"],["usize"],["u64"]]]],[10,"rc","","Return a new object containing the reverse complement of the sequence",41,[[["self"]],["self"]]],[11,"iter","","Iterate over the bases in the sequence",41,[[["self"]],["meriter"]]],[8,"Kmer","","Encapsulates a Kmer sequence with statically known K.",N,N],[10,"empty","","Create a Kmer initialized to all A's",42,[[],["self"]]],[10,"k","","K value for this concrete type.",42,[[],["usize"]]],[10,"to_u64","","Return the rank of this kmer in an lexicographic ordering of all kmers E.g. 'AAAA' -> 0, 'AAAT' -> 1, etc. This will panic if K > 32.",42,[[["self"]],["u64"]]],[10,"from_u64","","",42,[[["u64"]],["self"]]],[10,"extend_left","","Add the base `v` to the left side of the sequence, and remove the rightmost base",42,[[["self"],["u8"]],["self"]]],[10,"extend_right","","Add the base `v` to the right side of the sequence, and remove the leftmost base",42,[[["self"],["u8"]],["self"]]],[11,"extend","","Add the base `v` to the side of sequence given by `dir`, and remove a base at the opposite side",42,[[["self"],["u8"],["dir"]],["self"]]],[11,"get_extensions","","Generate all the extension of this sequence given by `exts` in direction `Dir`",42,[[["self"],["exts"],["dir"]],["vec"]]],[11,"min_rc_flip","","Return the minimum of the kmer and it's reverse complement, and a flag indicating if sequence was flipped",42,N],[11,"min_rc","","",42,[[["self"]],["self"]]],[11,"is_palindrome","","Test if this Kmer and it's reverse complement are the same",42,[[["self"]],["bool"]]],[11,"from_bytes","","Create a Kmer from the first K bytes of `bytes`, which must be encoded as the integers 0-4.",42,N],[11,"from_ascii","","Create a Kmer from the first K bytes of `bytes`, which must be encoded as ASCII letters A,C,G, or T.",42,N],[11,"to_string","","Return String containing Kmer sequence",42,[[["self"]],["string"]]],[11,"kmers_from_bytes","","Generate vector of all kmers contained in `str` encoded as 0-4.",42,N],[11,"kmers_from_ascii","","Generate vector of all kmers contained in `str`, encoded as ASCII ACGT.",42,N],[8,"MerImmut","","An immutable interface to a Mer sequence.",N,N],[11,"set","","",43,[[["self"],["usize"],["u8"]],["self"]]],[11,"set_slice","","",43,[[["self"],["usize"],["usize"],["u64"]],["self"]]],[8,"Vmer","","A DNA sequence with run-time variable length, up to a statically known maximum length",N,N],[10,"new","","Create a new sequence with length `len`, initialized to all A's",44,[[["usize"]],["self"]]],[10,"max_len","","Maximum sequence length that can be stored in this type",44,[[],["usize"]]],[11,"from_slice","","Create a Vmer from a sequence of bytes",44,N],[10,"get_kmer","","Efficiently extract a Kmer from the sequence",44,[[["self"],["usize"]],["k"]]],[11,"first_kmer","","Get the first Kmer from the sequence",44,[[["self"]],["k"]]],[11,"last_kmer","","Get the last kmer in the sequence",44,[[["self"]],["k"]]],[11,"both_term_kmer","","Get the terminal kmer of the sequence, on the both side of the sequence",44,N],[11,"term_kmer","","Get the terminal kmer of the sequence, on the side of the sequence given by dir",44,[[["self"],["dir"]],["k"]]],[11,"iter_kmers","","Iterate over the kmers in the sequence",44,[[["self"]],["kmeriter"]]],[11,"iter_kmer_exts","","Iterate over the kmers and their extensions, given the extensions of the whole sequence",44,[[["self"],["exts"]],["kmerextsiter"]]],[11,"next","","",45,[[["self"]],["option",["u8"]]]],[11,"fmt","","",0,[[["self"],["formatter"]],["result"]]],[11,"clone","","",0,[[["self"]],["dnabytes"]]],[11,"eq","","",0,[[["self"],["dnabytes"]],["bool"]]],[11,"ne","","",0,[[["self"],["dnabytes"]],["bool"]]],[11,"cmp","","",0,[[["self"],["dnabytes"]],["ordering"]]],[11,"partial_cmp","","",0,[[["self"],["dnabytes"]],["option",["ordering"]]]],[11,"lt","","",0,[[["self"],["dnabytes"]],["bool"]]],[11,"le","","",0,[[["self"],["dnabytes"]],["bool"]]],[11,"gt","","",0,[[["self"],["dnabytes"]],["bool"]]],[11,"ge","","",0,[[["self"],["dnabytes"]],["bool"]]],[11,"len","","",0,[[["self"]],["usize"]]],[11,"get","","",0,[[["self"],["usize"]],["u8"]]],[11,"set_mut","","Set base at `pos` to 2-bit encoded base `val`",0,[[["self"],["usize"],["u8"]]]],[11,"set_slice_mut","","Set `nbases` positions in the sequence, starting at `pos`. Values must  be packed into the upper-most bits of `value`.",0,[[["self"],["usize"],["usize"],["u64"]]]],[11,"rc","","Return a new object containing the reverse complement of the sequence",0,[[["self"]],["self"]]],[11,"new","","Create a new sequence with length `len`, initialized to all A's",0,[[["usize"]],["self"]]],[11,"max_len","","Maximum sequence length that can be stored in this type",0,[[],["usize"]]],[11,"get_kmer","","Efficiently extract a Kmer from the sequence",0,[[["self"],["usize"]],["k"]]],[11,"fmt","","",1,[[["self"],["formatter"]],["result"]]],[11,"eq","","",1,[[["self"],["dnaslice"]],["bool"]]],[11,"ne","","",1,[[["self"],["dnaslice"]],["bool"]]],[11,"cmp","","",1,[[["self"],["dnaslice"]],["ordering"]]],[11,"partial_cmp","","",1,[[["self"],["dnaslice"]],["option",["ordering"]]]],[11,"lt","","",1,[[["self"],["dnaslice"]],["bool"]]],[11,"le","","",1,[[["self"],["dnaslice"]],["bool"]]],[11,"gt","","",1,[[["self"],["dnaslice"]],["bool"]]],[11,"ge","","",1,[[["self"],["dnaslice"]],["bool"]]],[11,"len","","",1,[[["self"]],["usize"]]],[11,"get","","",1,[[["self"],["usize"]],["u8"]]],[11,"set_mut","","Set base at `pos` to 2-bit encoded base `val`",1,[[["self"],["usize"],["u8"]]]],[11,"set_slice_mut","","Set `nbases` positions in the sequence, starting at `pos`. Values must  be packed into the upper-most bits of `value`.",1,[[["self"],["usize"],["usize"],["u64"]]]],[11,"rc","","Return a new object containing the reverse complement of the sequence",1,[[["self"]],["self"]]],[11,"new","","Create a new sequence with length `len`, initialized to all A's",1,[[["usize"]],["self"]]],[11,"max_len","","Maximum sequence length that can be stored in this type",1,[[],["usize"]]],[11,"get_kmer","","Efficiently extract a Kmer from the sequence",1,[[["self"],["usize"]],["k"]]],[11,"clone","","",3,[[["self"]],["dir"]]],[11,"fmt","","",3,[[["self"],["formatter"]],["result"]]],[11,"flip","","Return a fresh Dir with the opposite direction",3,[[["self"]],["dir"]]],[11,"cond_flip","","Return a fresh Dir opposite direction if do_flip == True",3,[[["self"],["bool"]],["dir"]]],[11,"pick","","Pick between two alternatives, depending on the direction",3,[[["self"],["t"],["t"]],["t"]]],[11,"eq","","",2,[[["self"],["exts"]],["bool"]]],[11,"ne","","",2,[[["self"],["exts"]],["bool"]]],[11,"clone","","",2,[[["self"]],["exts"]]],[11,"cmp","","",2,[[["self"],["exts"]],["ordering"]]],[11,"partial_cmp","","",2,[[["self"],["exts"]],["option",["ordering"]]]],[11,"lt","","",2,[[["self"],["exts"]],["bool"]]],[11,"le","","",2,[[["self"],["exts"]],["bool"]]],[11,"gt","","",2,[[["self"],["exts"]],["bool"]]],[11,"ge","","",2,[[["self"],["exts"]],["bool"]]],[11,"hash","","",2,N],[11,"new","","",2,[[["u8"]],["self"]]],[11,"empty","","",2,[[],["exts"]]],[11,"from_single_dirs","","",2,[[["exts"],["exts"]],["exts"]]],[11,"merge","","",2,[[["exts"],["exts"]],["exts"]]],[11,"add","","",2,[[["self"],["exts"]],["exts"]]],[11,"set","","",2,[[["self"],["dir"],["u8"]],["exts"]]],[11,"get","","",2,[[["self"],["dir"]],["vec",["u8"]]]],[11,"has_ext","","",2,[[["self"],["dir"],["u8"]],["bool"]]],[11,"from_slice_bounds","","",2,N],[11,"from_dna_string","","",2,[[["dnastring"],["usize"],["usize"]],["exts"]]],[11,"num_exts_l","","",2,[[["self"]],["u8"]]],[11,"num_exts_r","","",2,[[["self"]],["u8"]]],[11,"num_ext_dir","","",2,[[["self"],["dir"]],["u8"]]],[11,"mk_left","","",2,[[["u8"]],["exts"]]],[11,"mk_right","","",2,[[["u8"]],["exts"]]],[11,"mk","","",2,[[["u8"],["u8"]],["exts"]]],[11,"get_unique_extension","","",2,[[["self"],["dir"]],["option",["u8"]]]],[11,"single_dir","","",2,[[["self"],["dir"]],["exts"]]],[11,"complement","","Complement the extension bases for each direction",2,[[["self"]],["exts"]]],[11,"reverse","","",2,[[["self"]],["exts"]]],[11,"rc","","",2,[[["self"]],["exts"]]],[11,"fmt","","",2,[[["self"],["formatter"]],["result"]]],[11,"next","","",46,[[["self"]],["option"]]],[11,"next","","",47,[[["self"]],["option"]]]],"paths":[[3,"DnaBytes"],[3,"DnaSlice"],[3,"Exts"],[4,"Dir"],[3,"IntKmer"],[3,"VarIntKmer"],[8,"IntHelp"],[8,"KmerSize"],[3,"K48"],[3,"K40"],[3,"K30"],[3,"K24"],[3,"K20"],[3,"K14"],[3,"K6"],[3,"K5"],[3,"K4"],[3,"K3"],[3,"K2"],[3,"DnaStringSlice"],[3,"PackedDnaStringSet"],[3,"DnaString"],[3,"DnaStringIter"],[3,"BaseGraph"],[3,"DebruijnGraph"],[3,"NodeKmer"],[3,"Node"],[3,"NodeIter"],[3,"NodeIntoIter"],[3,"NodeKmerIter"],[8,"Array"],[3,"Lmer"],[3,"MspInterval"],[8,"KmerSummarizer"],[3,"CountFilter"],[3,"CountFilterSet"],[3,"CountFilterEqClass"],[8,"CompressionSpec"],[3,"SimpleCompress"],[3,"ScmapCompress"],[3,"CleanGraph"],[8,"Mer"],[8,"Kmer"],[8,"MerImmut"],[8,"Vmer"],[3,"MerIter"],[3,"KmerIter"],[3,"KmerExtsIter"]]};
initSearch(searchIndex);
