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

pub mod block;
pub mod counts;
pub mod error;
mod value_reader;

use std::collections::HashMap;
use std::default::Default;
use std::fs::File;
use std::io::{BufReader, Cursor, SeekFrom};
use std::ops::{Bound, Deref, RangeBounds};
use std::path::Path;

use crate::block::{Block, Blocks};
use crate::counts::{BaseCounts, BasePercentages};
use crate::error::{Error, Result};
use crate::value_reader::{Reader, ValueReader};

/// Read data from a 2bit file
///
/// Usage:
///
/// ```
/// use twobit::TwoBitFile;
///
/// // softmask is disabled by default (see below for details)
/// let mut tb = TwoBitFile::open("assets/foo.2bit").unwrap();
/// let chromosome_lengths = tb.chroms();
/// assert_eq!(chromosome_lengths["chr1"], 150);
/// assert_eq!(chromosome_lengths["chr2"], 100);
/// ```
///
/// 2bit files offer two types of masks: N masks (aka hard masks) for unknown or arbitrary nucleotides
/// and soft masks for lower-case nucleotides (e.g. "t" instead of "T").
///
/// ```
/// # use twobit::TwoBitFile;
/// let mut tb_soft = TwoBitFile::open("assets/foo.2bit").unwrap().enable_softmask(true);
/// let expected_seq = "NNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTagctagctGATC"; // some lower and some upper case
/// assert_eq!(tb_soft.sequence("chr1", 24..74).unwrap(), expected_seq);
/// ```
/// It is not possible to disable hard masks but you can disable soft masks:
/// ```
/// # use twobit::TwoBitFile;
/// let mut tb_nosoft = TwoBitFile::open("assets/foo.2bit").unwrap().enable_softmask(false);
/// let expected_seq = "NNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTAGCTAGCTGATC"; // all upper case
/// assert_eq!(tb_nosoft.sequence("chr1", 24..74).unwrap(), expected_seq);
/// ```
///
/// Every sequence-related method accepts a range; to access the full sequence, pass `..`.
pub struct TwoBitFile<R: Reader> {
    reader: ValueReader<R>,
    sequences: SequenceRecords,
    softmask_enabled: bool,
}

pub type BoxTwoBitFile = TwoBitFile<Box<dyn Reader>>;

#[derive(Debug, Clone)]
pub(crate) struct SequenceRecord {
    offset: u64,
    length: usize,
    blocks_n: Blocks,
    blocks_soft_mask: Blocks,
}

// This wrapper is needed to avoid lifetime problems
#[derive(Debug, Clone)]
struct SequenceRecords(HashMap<String, SequenceRecord>);

impl SequenceRecords {
    #[inline]
    pub fn query(&self, chr: &str) -> Result<&SequenceRecord> {
        self.0
            .get(chr)
            .ok_or_else(|| Error::MissingName(chr.to_string()))
    }
}

impl Deref for SequenceRecords {
    type Target = HashMap<String, SequenceRecord>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

/// General information about a 2bit file
#[derive(Debug, PartialEq)]
pub struct TwoBitFileInfo {
    /// File size of the 2bit file
    pub file_size: u64,
    /// Number of chromosomes (or sequences) in the file
    pub num_chromosomes: usize,
    /// Total number of nucleotides in the file
    pub total_sequence_length: usize,
    /// Number of hard masks
    pub hard_masks_count: usize,
    /// Number of soft masks
    pub soft_masks_count: usize,
}

impl TwoBitFile<BufReader<File>> {
    /// Open a 2bit file from a given file path.
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::from_value_reader(ValueReader::open(path)?)
    }
}

impl<T> TwoBitFile<Cursor<T>>
where
    Cursor<T>: Reader,
{
    /// Open a 2bit file from a given in-memory buffer.
    pub fn from_buf(buf: T) -> Result<Self> {
        Self::from_value_reader(ValueReader::from_buf(buf)?)
    }
}

impl TwoBitFile<Cursor<Vec<u8>>> {
    /// Open a 2bit file from a given file path and read all of it into memory.
    pub fn open_and_read<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::from_value_reader(ValueReader::open_and_read(path)?)
    }
}

impl<R: Reader> TwoBitFile<R> {
    /// Open a 2bit file from a given generic reader (softmask is disabled by default).
    pub fn new(reader: R) -> Result<Self> {
        Self::from_value_reader(ValueReader::new(reader)?)
    }

    /// Box the reader (useful for type erasure if using multiple reader types).
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

    /// Enable/disable softmask: if enabled, return lower case nucleotides for soft blocks.
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

        let mut sequences = HashMap::new();

        let sequence_count = reader.field()?;
        let _reserved = reader.field(); // read the unused reserved field

        for _ in 0..sequence_count {
            let name_size = reader.byte()? as usize;
            let name = reader.string(name_size)?;
            let seq_offset = u64::from(reader.field()?);
            let offset = reader.tell()?;
            reader.seek(SeekFrom::Start(seq_offset))?;
            let seq_record = reader.sequence_record()?;
            reader.seek(SeekFrom::Start(offset))?;
            sequences.insert(name, seq_record);
        }

        Ok(Self {
            reader,
            sequences: SequenceRecords(sequences),
            softmask_enabled: false,
        })
    }

    /// Get the sizes of chromosomes in a 2bit file as a `HashMap`
    pub fn chroms(&mut self) -> HashMap<String, usize> {
        self.sequences
            .iter()
            .map(|(k, v)| (k.clone(), v.length))
            .collect()
    }

    /// Count bases of a partial chromosome
    ///
    /// * `chr` - Name of the chromosome from the 2bit file
    /// * `range` - Selected range of the sequence (`..` for full sequence)
    pub fn bases(&mut self, chr: &str, range: impl RangeBounds<usize>) -> Result<BaseCounts> {
        let nucs = self.sequence(chr, range)?;
        let mut result = BaseCounts::default();
        for nuc in nucs.chars() {
            match nuc {
                'A' | 'a' => result.a += 1,
                'C' | 'c' => result.c += 1,
                'G' | 'g' => result.g += 1,
                'T' | 't' => result.t += 1,
                'N' => result.n += 1,
                _ => return Err(Error::BadNucleotide(nuc)),
            }
        }
        Ok(result)
    }

    /// Calculates percentages of bases of a partial chromosome
    ///
    /// * `chr` - Name of the chromosome from the 2bit file
    /// * `range` - Selected range of the sequence (`..` for full sequence)
    pub fn bases_percentages(
        &mut self,
        chr: &str,
        range: impl RangeBounds<usize>,
    ) -> Result<BasePercentages> {
        Ok(self.bases(chr, range)?.into())
    }

    /// Obtain general information on the 2bit file you are dealing with
    pub fn info(&mut self) -> Result<TwoBitFileInfo> {
        let (total_length, hard_masks_count, soft_masks_count) =
            self.sequences.values().fold((0, 0, 0), |acc, seq| {
                (
                    acc.0 + seq.length,
                    acc.1 + seq.blocks_n.count(),
                    acc.2 + seq.blocks_soft_mask.count(),
                )
            });
        Ok(TwoBitFileInfo {
            file_size: self.reader.stream_len()?,
            num_chromosomes: self.sequences.len(),
            total_sequence_length: total_length,
            hard_masks_count,
            soft_masks_count,
        })
    }

    /// Get hard blocks (N-blocks) of a region on a chromosome
    ///
    /// * `chr` - Name of the chromosome from the 2bit file
    /// * `range` - Selected range of the sequence (`..` for full sequence)
    pub fn hard_masked_blocks(
        &mut self,
        chr: &str,
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

    /// Get soft blocks (lower case blocks) of a region on a chromosome
    ///
    /// * `chr` - Name of the chromosome from the 2bit file
    /// * `range` - Selected range of the sequence (`..` for full sequence)
    pub fn soft_masked_blocks(
        &mut self,
        chr: &str,
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

    /// Get a partial sequence of a chromosome
    ///
    /// * `chr` - Name of the chromosome from the 2bit file
    /// * `range` - Selected range of the sequence (`..` for full sequence)
    pub fn sequence(&mut self, chr: &str, range: impl RangeBounds<usize>) -> Result<String> {
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
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn test_chroms() {
        run_test(true, |mut bit| {
            let mut chr_count = 0;
            for (chr, len) in bit.chroms() {
                match chr.as_ref() {
                    "chr1" => assert_eq!(len, 150),
                    "chr2" => assert_eq!(len, 100),
                    _ => assert!(false), // unexpected chromosome
                }
                chr_count += 1;
            }
            assert_eq!(chr_count, 2);
            Ok(())
        });
    }

    #[test]
    fn test_sequence() {
        run_test(true, |mut bit| {
            let seq = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTagctagctGATCGATCGTAGCTAGCTAGCTAGCTGATCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN".to_string();
            assert_eq!(seq, bit.sequence("chr1", ..)?);
            let seq = "NNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTagctagctGATC";
            assert_eq!(seq, bit.sequence("chr1", 24..74)?);
            // let's try different offsets to test positions not divisible by 4
            let seq = "ACGTACGTagctagctGATC";
            assert_eq!(seq, bit.sequence("chr1", 54..74)?);
            let seq = "CGTACGTagctagctGATC";
            assert_eq!(seq, bit.sequence("chr1", 55..74)?);
            let seq = "GTACGTagctagctGATC";
            assert_eq!(seq, bit.sequence("chr1", 56..74)?); // divisible by 4
            let seq = "TACGTagctagctGATC";
            assert_eq!(seq, bit.sequence("chr1", 57..74)?);
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
            let full_percentages = BasePercentages {
                a: 0.08,
                c: 0.08,
                t: 0.08666666666666667,
                g: 0.08666666666666667,
                n: 100.0 / 150.0,
            };
            assert_eq!(bit.bases("chr1", ..)?, full_counts);
            assert_eq!(bit.bases_percentages("chr1", ..)?, full_percentages);

            let partial_counts = BaseCounts {
                a: 6,
                c: 6,
                t: 6,
                g: 6,
                n: 26,
            };
            let partial_percentages = BasePercentages {
                a: 0.12,
                c: 0.12,
                t: 0.12,
                g: 0.12,
                n: 26.0 / 50.0,
            };
            assert_eq!(bit.bases("chr1", 24..74)?, partial_counts);
            assert_eq!(bit.bases_percentages("chr1", 24..74)?, partial_percentages);

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
                hard_masks_count: 150,
                soft_masks_count: 8,
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
