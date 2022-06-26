//! The `twobit` crate provides an efficient 2bit file reader, implemented in pure Rust.
//!
//! # Brief overview
//!
//! This crate is inspired by [py2bit](https://github.com/deeptools/py2bit) and tries to
//! offer somewhat similar functionality with no C-dependency, no external crate dependencies,
//! and great performance. It follows
//! [2 bit specification version 0](http://genome.ucsc.edu/FAQ/FAQformat.html#format7).
//!
//! The primary type in this crate is [`TwoBitFile`](struct.TwoBitFile.html), a wrapper around
//! a generic IO reader that provides access to 2bit reading routines. Resulting nucleotide
//! sequences are returned as `String`, but can be also handled as bytes since they are
//! guaranteed to be pure ASCII.
//!
//! The set of errors is described by the [`Error`](struct.Error.html) type, with most methods
//! returning results wrapped in [`Result`](type.Result.html) due to possible IO and format errors.
//!
//! # Examples
//!
//! ```
//! use twobit::TwoBitFile;
//!
//! # || -> twobit::Result<()> {
//! let mut tb = TwoBitFile::open("assets/foo.2bit")?;
//! assert_eq!(tb.chrom_names(), &["chr1", "chr2"]);
//! assert_eq!(tb.chrom_sizes(), &[150, 100]);
//! let expected_seq = "NNACGTACGTACGTAGCTAGCTGATC";
//! assert_eq!(tb.read_sequence("chr1", 48..74)?, expected_seq);
//! # Ok(())
//! # }();
//! ```
//!
//! All sequence-related methods expect range argument; one can pass `..` (unbounded range)
//! in order to query the entire sequence:
//!
//! ```
//! # use twobit::TwoBitFile;
//! # || -> twobit::Result<()> {
//! # let mut tb = TwoBitFile::open("assets/foo.2bit")?;
//! assert_eq!(tb.read_sequence("chr1", ..)?.len(), 150);
//! # Ok(())
//! # }();
//! ```
//!
//! Files can be fully cached in memory in order to provide fast random access and avoid any
//! IO operations when decoding:
//!
//! ```
//! # use twobit::TwoBitFile;
//! # || -> twobit::Result<()> {
//! # let mut tb = TwoBitFile::open("assets/foo.2bit")?;
//! let mut tb_mem = TwoBitFile::open_and_read("assets/foo.2bit")?;
//! let expected_seq = tb.read_sequence("chr1", ..)?;
//! assert_eq!(tb_mem.read_sequence("chr1", ..)?, expected_seq);
//! # Ok(())
//! # }();
//! ```
//!
//! 2bit files offer two types of masks: N masks (aka hard masks) for unknown or arbitrary
//! nucleotides, and soft masks for lower-case nucleotides (e.g. "t" instead of "T").
//!
//! Hard masks are *always enabled*; soft masks are *disabled by default*, but can be enabled
//! manually:
//!
//! ```
//! # use twobit::TwoBitFile;
//! # || -> twobit::Result<()> {
//! # let mut tb = TwoBitFile::open("assets/foo.2bit")?;
//! let mut tb_soft = tb.enable_softmask(true);
//! let expected_seq = "NNACGTACGTACGTagctagctGATC";
//! assert_eq!(tb_soft.read_sequence("chr1", 48..74)?, expected_seq);
//! # Ok(())
//! # }();
//! ```

#![warn(clippy::all, clippy::pedantic, clippy::cargo)]
#![allow(
    clippy::missing_errors_doc,
    clippy::module_name_repetitions,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation
)]

mod block;
pub mod convert;
mod counts;
mod error;
pub mod nucleotide;
mod quad;
mod reader;

use std::cmp::min;
use std::default::Default;
use std::fs::File;
use std::io::{BufReader, Cursor, Read, Seek, SeekFrom};
use std::ops::{Bound, Deref, RangeBounds};
use std::path::Path;

use crate::block::{Block, Blocks};
use crate::quad::NUC_QUAD_LOOKUP;
use crate::reader::{Reader, ValueReader};

pub use crate::counts::BaseCounts;
pub use crate::error::{Error, Result};

/// 2bit file reader, a wrapper around [`Read`](std::io::Read) + [`Seek`](std::io::Seek).
pub struct TwoBitFile<R: Read + Seek> {
    reader: ValueReader<R>,
    sequences: SequenceRecords,
    softmask_enabled: bool,
}

/// A type alias for `TwoBitFile` with arbitrary boxed reader.
pub type BoxTwoBitFile = TwoBitFile<Box<dyn Reader>>;

/// A type alias for `TwoBitFile<_>` returned by `open()`.
pub type TwoBitPhysicalFile = TwoBitFile<BufReader<File>>;

/// A type alias for `TwoBitFile<_>` returned by `open_and_read()`.
pub type TwoBitMemoryFile = TwoBitFile<Cursor<Vec<u8>>>;

#[derive(Debug, Clone)]
pub(crate) struct SequenceRecord {
    chr: String,
    offset: u64,
    length: usize,
    blocks_n: Blocks,
    blocks_soft_mask: Blocks,
}

// This wrapper is needed to avoid lifetime problems
#[derive(Debug, Clone)]
pub(crate) struct SequenceRecords(Vec<SequenceRecord>);

impl SequenceRecords {
    #[inline]
    pub fn query(&self, chr: impl AsRef<str>) -> Result<&SequenceRecord> {
        let chr = chr.as_ref();
        self.0
            .iter()
            .find(|seq| seq.chr == chr)
            .ok_or_else(|| Error::MissingName(chr.to_string()))
    }
}

impl Deref for SequenceRecords {
    type Target = [SequenceRecord];

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

/// Information on a particular chromosome or sequence in a 2bit file.
#[derive(Debug, Clone, PartialEq)]
pub struct SequenceInfo {
    /// Chromosome or sequence name
    pub chr: String,
    /// Number of nucleotides in the sequence
    pub length: usize,
    /// Total length of all hard masks
    pub hard_masks_total_length: usize,
    /// Total length of all hard masks
    pub soft_masks_total_length: usize,
}

impl SequenceRecord {
    pub fn info(&self) -> SequenceInfo {
        SequenceInfo {
            chr: self.chr.clone(),
            length: self.length,
            hard_masks_total_length: self.blocks_n.count(),
            soft_masks_total_length: self.blocks_soft_mask.count(),
        }
    }
}

/// Summary information about a 2bit file.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TwoBitFileInfo {
    /// File size in bytes.
    pub file_size: u64,
    /// Number of chromosomes or sequences.
    pub num_chromosomes: usize,
    /// Total number of nucleotides over all sequences.
    pub total_sequence_length: usize,
    /// Total length of all hard masks in nucleotides.
    pub hard_masks_total_length: usize,
    /// Total length of all soft masks in nucleotides.
    pub soft_masks_total_length: usize,
}

impl TwoBitFileInfo {
    pub(crate) const fn new(file_size: u64) -> Self {
        Self {
            file_size,
            num_chromosomes: 0,
            total_sequence_length: 0,
            hard_masks_total_length: 0,
            soft_masks_total_length: 0,
        }
    }

    pub(crate) fn update(&mut self, seq: &SequenceRecord) {
        self.num_chromosomes += 1;
        self.total_sequence_length += seq.length;
        self.hard_masks_total_length += seq.blocks_n.count();
        self.soft_masks_total_length += seq.blocks_soft_mask.count();
    }
}

impl TwoBitPhysicalFile {
    /// Opens a 2bit file from a given file path.
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::from_value_reader(ValueReader::open(path)?)
    }
}

impl<T> TwoBitFile<Cursor<T>>
where
    Cursor<T>: Reader,
{
    /// Opens a 2bit file from a given in-memory buffer.
    pub fn from_buf(buf: T) -> Result<Self> {
        Self::from_value_reader(ValueReader::from_buf(buf)?)
    }
}

impl TwoBitMemoryFile {
    /// Opens a 2bit file from a given file path and reads all of it into memory.
    pub fn open_and_read<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::from_value_reader(ValueReader::open_and_read(path)?)
    }
}

impl<R: Read + Seek> TwoBitFile<R> {
    /// Creates a 2bit file reader from a given generic reader.
    pub fn new(reader: R) -> Result<Self> {
        Self::from_value_reader(ValueReader::new(reader)?)
    }

    /// Returns a file handles with a boxed internal reader.
    ///
    /// This is useful for cases when the exact reader type needs to be erased.
    pub fn boxed(self) -> BoxTwoBitFile
    where
        R: 'static,
    {
        TwoBitFile {
            reader: self.reader.boxed(),
            sequences: self.sequences,
            softmask_enabled: self.softmask_enabled,
        }
    }

    /// Enables or disables soft masks.
    ///
    /// If enabled, lower case nucleotides are returned for soft blocks.
    ///
    /// Note: this option is *disabled by default*.
    #[must_use]
    pub fn enable_softmask(self, softmask_enabled: bool) -> Self {
        Self {
            softmask_enabled,
            ..self
        }
    }

    fn from_value_reader(mut reader: ValueReader<R>) -> Result<Self> {
        let sequences = reader.sequence_records()?;
        Ok(Self {
            reader,
            sequences,
            softmask_enabled: false,
        })
    }

    /// Reads a partial sequence of a chromosome.
    ///
    /// * `chr` – name of the chromosome
    /// * `range` – selected range of the sequence (`..` for full sequence)
    pub fn read_sequence(
        &mut self,
        chr: impl AsRef<str>,
        range: impl RangeBounds<usize>,
    ) -> Result<String> {
        const NUC: &[u8; 4] = b"TCAG";

        let seq = self.sequences.query(&chr)?;
        let reader = &mut self.reader;

        let start = {
            let pos = match range.start_bound() {
                Bound::Included(&v) => v,
                Bound::Excluded(&v) => v + 1,
                Bound::Unbounded => 0,
            };
            min(pos, seq.length)
        };
        let end = {
            let pos = match range.end_bound() {
                Bound::Included(&v) => v + 1,
                Bound::Excluded(&v) => v,
                Bound::Unbounded => seq.length,
            };
            min(pos, seq.length)
        };

        if start >= end {
            return Ok(String::new()); // trivial case, empty return result
        }
        let first_byte = start / 4;
        reader.seek(SeekFrom::Start(seq.offset))?; // beginning of the DNA sequence
        reader.seek(SeekFrom::Current(first_byte as _))?; // position where we want to start reading

        let length = end - start;
        let mut out = vec![0_u8; length];

        let last_byte = (end - 1) / 4; // inclusive (!) index, and here we know that end >= 1
        let skip_start = start % 4; // number of pairs to skip in the first byte

        let parse_byte = |buf: &mut [u8], mut byte: u8| {
            for v in buf {
                *v = NUC[((byte & 192) >> 6) as usize];
                byte <<= 2;
            }
        };

        if first_byte == last_byte {
            parse_byte(&mut out, reader.byte()? << (skip_start * 2)); // special case (single byte)
        } else {
            // most common case where the sequence is spread over 2 or more bytes
            if skip_start != 0 {
                parse_byte(
                    &mut out[..4 - skip_start],
                    reader.byte()? << (skip_start * 2),
                );
            }
            for chunk in out[(4 - skip_start) % 4..].chunks_exact_mut(4) {
                chunk.copy_from_slice(&NUC_QUAD_LOOKUP[reader.byte()? as usize]);
            }
            let include_end = ((end - 1) % 4) + 1; // number of pairs to include in the last byte
            if include_end != 4 {
                parse_byte(&mut out[length - include_end..], reader.byte()?);
            }
        }

        seq.blocks_n.apply_masks::<true>(&mut out, start);
        if self.softmask_enabled {
            seq.blocks_soft_mask.apply_masks::<false>(&mut out, start);
        }

        Ok(unsafe { String::from_utf8_unchecked(out) }) // we know it's ascii so it's ok
    }

    /// Obtain summary information on the 2bit file.
    pub fn info(&mut self) -> Result<TwoBitFileInfo> {
        let mut info = TwoBitFileInfo::new(self.reader.stream_len()?);
        self.sequences.iter().for_each(|seq| info.update(seq));
        Ok(info)
    }

    /// Obtain information on each sequence in the file.
    pub fn sequence_info(&self) -> Vec<SequenceInfo> {
        self.sequences.iter().map(SequenceRecord::info).collect()
    }

    /// Returns names of all chromosomes.
    pub fn chrom_names(&self) -> Vec<String> {
        self.sequences.iter().map(|seq| seq.chr.clone()).collect()
    }

    /// Returns sizes of all chromosomes.
    pub fn chrom_sizes(&self) -> Vec<usize> {
        self.sequences.iter().map(|seq| seq.length).collect()
    }

    /// Count bases in a partial sequence of a chromosome.
    ///
    /// * `chr` – name of the chromosome
    /// * `range` – selected range of the sequence (`..` for full sequence)
    pub fn base_counts(
        &mut self,
        chr: impl AsRef<str>,
        range: impl RangeBounds<usize>,
    ) -> Result<BaseCounts<usize>> {
        let mut counts = BaseCounts::default();
        for &nuc in self.read_sequence(chr, range)?.as_bytes() {
            counts.update(nuc)?;
        }
        Ok(counts)
    }

    /// Get hard blocks (N-blocks) of a region on a chromosome.
    ///
    /// * `chr` – name of the chromosome
    /// * `range` – selected range of the sequence (`..` for full sequence)
    pub fn hard_masked_blocks(
        &mut self,
        chr: impl AsRef<str>,
        range: impl RangeBounds<usize>,
    ) -> Result<Vec<Block>> {
        Ok(self
            .sequences
            .query(chr)?
            .blocks_n
            .overlaps(range)
            .cloned()
            .collect())
    }

    /// Get soft blocks (lower case blocks) of a region of a chromosome.
    ///
    /// * `chr` – name of the chromosome
    /// * `range` – selected range of the sequence (`..` for full sequence)
    pub fn soft_masked_blocks(
        &mut self,
        chr: impl AsRef<str>,
        range: impl RangeBounds<usize>,
    ) -> Result<Vec<Block>> {
        Ok(self
            .sequences
            .query(chr)?
            .blocks_soft_mask
            .overlaps(range)
            .cloned()
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use super::error::Result;
    use super::reader::Reader;
    use super::{BaseCounts, SequenceInfo, TwoBitFile, TwoBitFileInfo};

    const TESTFILE: &str = "assets/foo.2bit";

    fn run_test(softmask_enabled: bool, func: impl Fn(TwoBitFile<Box<dyn Reader>>) -> Result<()>) {
        let mut files = vec![
            TwoBitFile::open(TESTFILE).expect("unit-test").boxed(),
            TwoBitFile::open_and_read(TESTFILE)
                .expect("unit-test")
                .boxed(),
        ];
        for mut tb in files.drain(..) {
            if softmask_enabled {
                tb = tb.enable_softmask(true);
            }
            func(tb).expect("unit-test");
        }
    }

    #[test]
    fn test_chrom() {
        run_test(true, |bit| {
            assert_eq!(bit.chrom_names(), vec!["chr1", "chr2"]);
            assert_eq!(bit.chrom_sizes(), vec![150, 100]);
            Ok(())
        });
    }

    #[test]
    fn test_sequence_info() {
        run_test(true, |bit| {
            let expected_info = vec![
                SequenceInfo {
                    chr: "chr1".to_owned(),
                    length: 150,
                    hard_masks_total_length: 100,
                    soft_masks_total_length: 8,
                },
                SequenceInfo {
                    chr: "chr2".to_owned(),
                    length: 100,
                    hard_masks_total_length: 50,
                    soft_masks_total_length: 0,
                },
            ];
            assert_eq!(bit.sequence_info(), expected_info);
            Ok(())
        });
    }

    #[test]
    fn test_read_sequence() {
        run_test(true, |mut bit| {
            let seq = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTagctagctGATCGATCGTAGCTAGCTAGCTAGCTGATCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
            assert_eq!(bit.read_sequence("chr1", ..)?, seq);
            let n = 12;
            for start in 0..n {
                for end_range in &[start..(start + n), (seq.len() - n)..seq.len()] {
                    for end in end_range.clone() {
                        if end > start {
                            assert_eq!(bit.read_sequence("chr1", start..end)?, &seq[start..end]);
                        }
                        if end > 0 {
                            assert_eq!(bit.read_sequence("chr1", ..end)?, &seq[..end]);
                        }
                        assert_eq!(bit.read_sequence("chr1", start..)?, &seq[start..]);
                        if end < seq.len() {
                            assert_eq!(bit.read_sequence("chr1", start..=end)?, &seq[start..=end]);
                            assert_eq!(bit.read_sequence("chr1", ..=end)?, &seq[..=end]);
                        }
                    }
                }
            }
            Ok(())
        });
    }

    #[test]
    fn test_bases() {
        run_test(true, |mut bit| {
            let full_counts = BaseCounts::new(12, 12, 13, 13, 100);
            let full_percentages = BaseCounts::new(0.08, 0.08, 13. / 150., 13. / 150., 100. / 150.);
            let c = bit.base_counts("chr1", ..)?;
            assert_eq!(c, full_counts);
            assert_eq!(c.percentages(), full_percentages);

            let partial_counts = BaseCounts::new(6, 6, 6, 6, 26);
            let partial_percentages = BaseCounts::new(0.12, 0.12, 0.12, 0.12, 26. / 50.);
            let c = bit.base_counts("chr1", 24..74)?;
            assert_eq!(c, partial_counts);
            assert_eq!(c.percentages(), partial_percentages);

            Ok(())
        });
    }

    #[test]
    fn test_info() {
        run_test(true, |mut bit| {
            let info = TwoBitFileInfo {
                file_size: 161,
                num_chromosomes: 2,
                total_sequence_length: 250,
                hard_masks_total_length: 150,
                soft_masks_total_length: 8,
            };
            assert_eq!(bit.info()?, info);
            Ok(())
        });
    }

    #[test]
    fn test_hard_masked_blocks() {
        run_test(true, |mut bit| {
            assert_eq!(bit.hard_masked_blocks("chr1", ..)?, vec![0..50, 100..150]);
            //TODO hard_masked_blocks()
            Ok(())
        });
    }

    #[test]
    fn test_soft_masked_blocks() {
        run_test(true, |mut bit| {
            assert_eq!(bit.soft_masked_blocks("chr1", ..)?, vec![62..70]);
            //TODO soft_masked_blocks()
            Ok(())
        });
    }

    //TODO IO Errors
}
