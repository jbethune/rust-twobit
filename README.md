# twobit

Efficient 2bit file reader, implemented in pure Rust.

[![Build](https://github.com/jbethune/rust-twobit/workflows/CI/badge.svg)](https://github.com/jbethune/rust-twobit/actions?query=branch%3Amaster)
[![Latest Version](https://img.shields.io/crates/v/twobit.svg)](https://crates.io/crates/twobit)
[![Documentation](https://docs.rs/twobit/badge.svg)](https://docs.rs/twobit)
![twobit: rustc 1.51+](https://img.shields.io/badge/twobit-rustc_1.51+-lightblue.svg)
[![MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

The [2bit file format](http://genome.ucsc.edu/FAQ/FAQformat.html#format7) is
used to store genomic sequences on disk. It allows for fast access to specific
parts of the genome.

This crate is inspired by [py2bit](https://github.com/deeptools/py2bit) and tries to
offer somewhat similar functionality with no C-dependency, no external crate dependencies,
and great performance. It follows
[2 bit specification version 0](http://genome.ucsc.edu/FAQ/FAQformat.html#format7).

## Examples

```rust
use twobit::TwoBitFile;

let mut tb = TwoBitFile::open("assets/foo.2bit")?;
assert_eq!(tb.chrom_names(), &["chr1", "chr2"]);
assert_eq!(tb.chrom_sizes(), &[150, 100]);
let expected_seq = "NNACGTACGTACGTAGCTAGCTGATC";
assert_eq!(tb.read_sequence("chr1", 48..74)?, expected_seq);
```

All sequence-related methods expect range argument; one can pass `..` (unbounded range)
in order to query the entire sequence:

```rust
assert_eq!(tb.read_sequence("chr1", ..)?.len(), 150);
```

Files can be fully cached in memory in order to provide fast random access and avoid any
IO operations when decoding:

```rust
let mut tb_mem = TwoBitFile::open_and_read("assets/foo.2bit")?;
let expected_seq = tb.read_sequence("chr1", ..)?;
assert_eq!(tb_mem.read_sequence("chr1", ..)?, expected_seq);
```

2bit files offer two types of masks: N masks (aka hard masks) for unknown or arbitrary
nucleotides, and soft masks for lower-case nucleotides (e.g. "t" instead of "T").

Hard masks are *always enabled*; soft masks are *disabled by default*, but can be enabled
manually:

```rust
let mut tb_soft = tb.enable_softmask(true);
let expected_seq = "NNACGTACGTACGTagctagctGATC";
assert_eq!(tb_soft.read_sequence("chr1", 48..74)?, expected_seq);
```