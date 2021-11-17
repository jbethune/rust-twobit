# twobit

2bit file format library in Rust.

[![Latest Version](https://img.shields.io/crates/v/twobit.svg)](https://crates.io/crates/twobit)
[![Documentation](https://docs.rs/twobit/badge.svg)](https://docs.rs/twobit)
![hdf5: rustc 1.51+](https://img.shields.io/badge/hdf5-rustc_1.51+-lightblue.svg)
[![MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

The [2bit file format](http://genome.ucsc.edu/FAQ/FAQformat.html#format7) is
used to store genomic sequences on disk. It allows for fast access to specific
parts of the genome.

This crate is inspired by [py2bit python package](https://github.com/deeptools/py2bit)
and does not have any external dependencies.
