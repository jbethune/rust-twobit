//! Read 2bit files in Rust.
//!
//! This crate is inspired by <a href="https://github.com/deeptools/py2bit">py2bit</a> and tries to
//! offer the same functionality. It is written from the ground up with no C-dependency. It follows
//! the <a href="http://genome.ucsc.edu/FAQ/FAQformat.html#format7">2bit specification version
//! 0</a>.
//!
//! Note that most methods perform IO operations and return a `Result` for that reason.
//!
//! The main entry point for users of this crate is the [`TwoBitFile`](struct.TwoBitFile.html)
//! struct. Please see its documentation for details how to use this crate.

#![warn(clippy::all, clippy::pedantic, clippy::nursery, clippy::cargo)]
#![allow(
    clippy::missing_errors_doc,
    clippy::module_name_repetitions,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation
)]

mod block;
mod counts;
mod error;
mod value_reader;

use std::default::Default;
use std::fs::File;
use std::io::{BufReader, Cursor, Read, Seek, SeekFrom};
use std::ops::{Bound, Deref, RangeBounds};
use std::path::Path;

use crate::block::{Block, Blocks};
use crate::value_reader::{Reader, ValueReader};

pub use crate::counts::BaseCounts;
pub use crate::error::{Error, Result};

/// 2bit file reader, a wrapper around [`Read`](std::io::Read) + [`Seek`](std::io::Seek).
///
/// Usage:
///
/// ```
/// use twobit::TwoBitFile;
///
/// // softmask is disabled by default (see below for details)
/// let mut tb = TwoBitFile::open("assets/foo.2bit").unwrap();
/// assert_eq!(tb.chromosomes(), &["chr1", "chr2"]);
/// assert_eq!(tb.chrom_sizes(), &[150, 100]);
/// ```
///
/// 2bit files offer two types of masks: N masks (aka hard masks) for unknown or arbitrary nucleotides
/// and soft masks for lower-case nucleotides (e.g. "t" instead of "T").
///
/// ```
/// # use twobit::TwoBitFile;
/// let mut tb_soft = TwoBitFile::open("assets/foo.2bit").unwrap().enable_softmask(true);
/// let expected_seq = "NNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTagctagctGATC"; // some lower and some upper case
/// assert_eq!(tb_soft.read_sequence("chr1", 24..74).unwrap(), expected_seq);
/// ```
/// It is not possible to disable hard masks but you can disable soft masks:
/// ```
/// # use twobit::TwoBitFile;
/// let mut tb_nosoft = TwoBitFile::open("assets/foo.2bit").unwrap().enable_softmask(false);
/// let expected_seq = "NNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTAGCTAGCTGATC"; // all upper case
/// assert_eq!(tb_nosoft.read_sequence("chr1", 24..74).unwrap(), expected_seq);
/// ```
///
/// Every sequence-related method accepts a range; to access the full sequence, pass `..`.
pub struct TwoBitFile<R: Read + Seek> {
    reader: ValueReader<R>,
    sequences: SequenceRecords,
    softmask_enabled: bool,
}

/// A type alias for `TwoBitFile` with arbitrary boxed reader.
pub type BoxTwoBitFile = TwoBitFile<Box<dyn Reader>>;

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
struct SequenceRecords(Vec<SequenceRecord>);

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
    pub(crate) fn new(file_size: u64) -> Self {
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

impl TwoBitFile<BufReader<File>> {
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

impl TwoBitFile<Cursor<Vec<u8>>> {
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
    /// Note: this option is disabled by default.
    #[must_use]
    pub fn enable_softmask(self, softmask_enabled: bool) -> Self {
        Self {
            softmask_enabled,
            ..self
        }
    }

    fn from_value_reader(mut reader: ValueReader<R>) -> Result<Self> {
        reader.seek_start()?; // rewind to the start of the file, skipping the header

        let mut sequences = vec![];

        let sequence_count = reader.field()?;
        let _reserved = reader.field(); // read the unused reserved field

        for _ in 0..sequence_count {
            let name_size = reader.byte()? as usize;
            let name = reader.string(name_size)?;
            let seq_offset = u64::from(reader.field()?);
            let offset = reader.tell()?;
            reader.seek(SeekFrom::Start(seq_offset))?;
            let seq_record = reader.sequence_record(&name)?;
            reader.seek(SeekFrom::Start(offset))?;
            sequences.push(seq_record);
        }

        Ok(Self {
            reader,
            sequences: SequenceRecords(sequences),
            softmask_enabled: false,
        })
    }

    /// Reads a partial sequence of a chromosome.
    ///
    /// * `chr` - Name of the chromosome
    /// * `range` - Selected range of the sequence (`..` for full sequence)
    pub fn read_sequence(
        &mut self,
        chr: impl AsRef<str>,
        range: impl RangeBounds<usize>,
    ) -> Result<String> {
        const NUC: &[u8; 4] = b"TCAG";

        let seq = self.sequences.query(chr)?;
        let reader = &mut self.reader;

        let start = match range.start_bound() {
            Bound::Included(&v) => v,
            Bound::Excluded(&v) => v + 1,
            Bound::Unbounded => 0,
        };
        let end = match range.end_bound() {
            Bound::Included(&v) => v + 1,
            Bound::Excluded(&v) => v,
            Bound::Unbounded => seq.length,
        };

        let first_byte = start / 4;
        reader.seek(SeekFrom::Start(seq.offset))?; // beginning of the DNA sequence
        reader.seek(SeekFrom::Current(first_byte as _))?; // position where we want to start reading
        if start >= end {
            return Ok(String::new()); // trivial case, empty return result
        }

        let length = end - start;
        let mut out = Vec::with_capacity(length);
        unsafe {
            out.set_len(length);
        }

        let last_byte = (end - 1) / 4; // inclusive (!) index, and here we know that end >= 1
        let skip_start = start % 4; // number of pairs to skip in the first byte
        let partial_start = skip_start != 0;

        let parse_byte =
            |buf: &mut [u8], r: &mut ValueReader<_>, skip, take, shift| -> Result<()> {
                let mut byte = r.byte()? << shift;
                for v in buf.iter_mut().skip(skip).take(take) {
                    *v = NUC[((byte & 192) >> 6) as usize];
                    byte <<= 2;
                }
                Ok(())
            };

        if first_byte == last_byte {
            parse_byte(&mut out, reader, 0, length, skip_start * 2)?; // special case (single byte)
        } else {
            // most common case where the sequence is spread over 2 or more bytes
            let include_end = ((end - 1) % 4) + 1; // number of pairs to include in the last byte
            let partial_end = include_end != 4;

            if partial_start {
                parse_byte(&mut out, reader, 0, 4 - skip_start, skip_start * 2)?;
            }
            let n_mid = (last_byte - first_byte + 1)
                - usize::from(partial_start)
                - usize::from(partial_end);
            let mut pos = (4 - skip_start) % 4;
            for _ in 0..n_mid {
                let byte = reader.byte()? as usize;
                unsafe {
                    *out.get_unchecked_mut(pos) = NUC[(byte & 192) >> 6];
                    *out.get_unchecked_mut(pos + 1) = NUC[(byte & 48) >> 4];
                    *out.get_unchecked_mut(pos + 2) = NUC[(byte & 12) >> 2];
                    *out.get_unchecked_mut(pos + 3) = NUC[byte & 3];
                }
                pos += 4;
            }
            if partial_end {
                parse_byte(&mut out, reader, pos, include_end, 0)?;
            }
        }

        seq.blocks_n.apply_masks::<true>(&mut out, start..end);
        if self.softmask_enabled {
            seq.blocks_soft_mask
                .apply_masks::<false>(&mut out, start..end);
        }

        Ok(unsafe { String::from_utf8_unchecked(out) }) // we know it's ascii so it's ok
    }

    /// Obtain summary information on the 2bit file.
    pub fn info(&mut self) -> Result<TwoBitFileInfo> {
        let mut info = TwoBitFileInfo::new(self.reader.stream_len()?);
        self.sequences.iter().for_each(|seq| info.update(seq));
        Ok(info)
    }

    /// Returns names of all chromosomes.
    pub fn chromosomes(&self) -> Vec<String> {
        self.sequences.iter().map(|seq| seq.chr.clone()).collect()
    }

    /// Returns sizes of all chromosomes.
    pub fn chrom_sizes(&mut self) -> Vec<usize> {
        self.sequences.iter().map(|seq| seq.length).collect()
    }

    /// Count bases in a partial sequence of a chromosome
    ///
    /// * `chr` - Name of the chromosome
    /// * `range` - Selected range of the sequence (`..` for full sequence)
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

    /// Get hard blocks (N-blocks) of a region on a chromosome
    ///
    /// * `chr` - Name of the chromosome
    /// * `range` - Selected range of the sequence (`..` for full sequence)
    pub fn hard_masked_blocks(
        &mut self,
        chr: impl AsRef<str>,
        range: impl RangeBounds<usize>,
    ) -> Result<Vec<Block>> {
        Ok(self
            .sequences
            .query(chr)?
            .blocks_n
            .iter_overlaps(range)
            .cloned()
            .collect())
    }

    /// Get soft blocks (lower case blocks) of a region of a chromosome
    ///
    /// * `chr` - Name of the chromosome
    /// * `range` - Selected range of the sequence (`..` for full sequence)
    pub fn soft_masked_blocks(
        &mut self,
        chr: impl AsRef<str>,
        range: impl RangeBounds<usize>,
    ) -> Result<Vec<Block>> {
        Ok(self
            .sequences
            .query(chr)?
            .blocks_soft_mask
            .iter_overlaps(range)
            .cloned()
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use super::error::Result;
    use super::value_reader::Reader;
    use super::{BaseCounts, TwoBitFile, TwoBitFileInfo};

    const TESTFILE: &str = "assets/foo.2bit";

    fn run_test(softmask_enabled: bool, func: impl Fn(TwoBitFile<Box<dyn Reader>>) -> Result<()>) {
        let mut files = vec![
            TwoBitFile::open(TESTFILE).unwrap().boxed(),
            TwoBitFile::open_and_read(TESTFILE).unwrap().boxed(),
        ];
        for mut tb in files.drain(..) {
            if softmask_enabled {
                tb = tb.enable_softmask(true);
            }
            func(tb).unwrap();
        }
    }

    #[test]
    fn test_chromosomes() {
        run_test(true, |mut bit| {
            assert_eq!(bit.chromosomes(), vec!["chr1", "chr2"]);
            assert_eq!(bit.chrom_sizes(), vec![150, 100]);
            Ok(())
        });
    }

    #[test]
    fn test_sequence() {
        run_test(true, |mut bit| {
            let seq = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTagctagctGATCGATCGTAGCTAGCTAGCTAGCTGATCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN".to_string();
            assert_eq!(seq, bit.read_sequence("chr1", ..)?);
            let seq = "NNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTagctagctGATC";
            assert_eq!(seq, bit.read_sequence("chr1", 24..74)?);
            // let's try different offsets to test positions not divisible by 4
            let seq = "ACGTACGTagctagctGATC";
            assert_eq!(seq, bit.read_sequence("chr1", 54..74)?);
            let seq = "CGTACGTagctagctGATC";
            assert_eq!(seq, bit.read_sequence("chr1", 55..74)?);
            let seq = "GTACGTagctagctGATC";
            assert_eq!(seq, bit.read_sequence("chr1", 56..74)?); // divisible by 4
            let seq = "TACGTagctagctGATC";
            assert_eq!(seq, bit.read_sequence("chr1", 57..74)?);
            Ok(())
        });
    }

    #[test]
    fn test_bases() {
        run_test(true, |mut bit| {
            let full_counts = BaseCounts {
                a: 12,
                c: 12,
                t: 13,
                g: 13,
                n: 100,
            };
            let full_percentages = BaseCounts {
                a: 0.08,
                c: 0.08,
                t: 0.08666666666666667,
                g: 0.08666666666666667,
                n: 100.0 / 150.0,
            };
            let c = bit.base_counts("chr1", ..)?;
            assert_eq!(c, full_counts);
            assert_eq!(c.percentages(), full_percentages);

            let partial_counts = BaseCounts {
                a: 6,
                c: 6,
                t: 6,
                g: 6,
                n: 26,
            };
            let partial_percentages = BaseCounts {
                a: 0.12,
                c: 0.12,
                t: 0.12,
                g: 0.12,
                n: 26.0 / 50.0,
            };
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
