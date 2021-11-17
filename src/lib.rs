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
use std::ops::Deref;
use std::path::Path;

use crate::block::Block;
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
/// assert_eq!(tb_soft.sequence("chr1", 24, 74).unwrap(), expected_seq);
/// ```
/// It is not possible to disable hard masks but you can disable soft masks:
/// ```
/// # use twobit::TwoBitFile;
/// let mut tb_nosoft = TwoBitFile::open("assets/foo.2bit").unwrap().enable_softmask(false);
/// let expected_seq = "NNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTAGCTAGCTGATC"; // all upper case
/// assert_eq!(tb_nosoft.sequence("chr1", 24, 74).unwrap(), expected_seq);
/// ```
///
/// There are 2 variants for every sequence-related method:
/// * Those that request information on parts of a chromosome
/// * Those that request information on a complete chromosome
///
/// The partial requests require you to provide a start (0-based, inclusive) and a stop (0-based,
/// exclusive) position as parameters. The full request methods all start with the prefix
/// `full_`.
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
    blocks_n: Vec<Block>,
    blocks_soft_mask: Vec<Block>,
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
    pub chromosomes: usize,
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

    /// Get the full sequence of a chromosome
    ///
    /// * `chr` Name of the chromosome from the 2bit file
    /// * returns the full sequence as a `String` on success
    pub fn full_sequence(&mut self, chr: &str) -> Result<String> {
        let dna_size = self.sequences.query(chr)?.length;
        self.read_sequence(chr, 0, dna_size as _)
    }

    /// Get a partial sequence of a chromosome
    ///
    /// * `chr` Name of the chromosome from the 2bit file
    /// * `start` Start position of the sequence
    /// * `end` Stop position of the sequence
    /// * returns the full sequence as a `String` on success
    pub fn sequence(&mut self, chr: &str, start: usize, end: usize) -> Result<String> {
        self.read_sequence(chr, start, end)
    }

    /// Count bases of a chromosome
    ///
    /// * `chr` Name of the chromosome from the 2bit file
    pub fn full_bases(&mut self, chr: &str) -> Result<BaseCounts> {
        let len = self.chr_length(chr)?;
        self.bases(chr, 0, len)
    }

    /// Calculates percentages of bases of a chromosome
    ///
    /// * `chr` Name of the chromosome from the 2bit file
    pub fn full_bases_percentages(&mut self, chr: &str) -> Result<BasePercentages> {
        let len = self.chr_length(chr)?;
        self.bases_percentages(chr, 0, len)
    }

    /// Count bases of a partial chromosome
    ///
    /// * `chr` Name of the chromosome from the 2bit file
    pub fn bases(&mut self, chr: &str, start: usize, end: usize) -> Result<BaseCounts> {
        let nucs = self.read_sequence(chr, start, end)?;
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
    /// * `chr` Name of the chromosome from the 2bit file
    pub fn bases_percentages(
        &mut self,
        chr: &str,
        start: usize,
        end: usize,
    ) -> Result<BasePercentages> {
        Ok(self.bases(chr, start, end)?.into())
    }

    /// Obtain general information on the 2bit file you are dealing with
    pub fn info(&mut self) -> Result<TwoBitFileInfo> {
        let mut total_length = 0;
        let mut hard_masks_count = 0;
        let mut soft_masks_count = 0;
        for chr in &self.sequences.keys().cloned().collect::<Vec<_>>() {
            total_length += self.chr_length(chr)?;
            let record = self.sequences.query(chr)?;
            hard_masks_count +=
                record.blocks_n.iter().fold(0, |acc, blk| acc + blk.length) as usize;
            soft_masks_count += record
                .blocks_soft_mask
                .iter()
                .fold(0, |acc, blk| acc + blk.length) as usize;
        }
        Ok(TwoBitFileInfo {
            file_size: self.reader.stream_len()?,
            chromosomes: self.sequences.len(),
            total_sequence_length: total_length,
            hard_masks_count,
            soft_masks_count,
        })
    }

    /// Get all hard blocks (N-blocks) of a chromosome
    ///
    /// * `chr` Name of the chromosome from the 2bit file
    pub fn full_hard_masked_blocks(&mut self, chr: &str) -> Result<Vec<Block>> {
        match self.sequences.query(chr) {
            Ok(record) => Ok(record.blocks_n.clone()),
            Err(e) => Err(e),
        }
    }

    /// Get hard blocks (N-blocks) of a region on a chromosome
    ///
    /// * `chr` Name of the chromosome from the 2bit file
    pub fn hard_masked_blocks(
        &mut self,
        chr: &str,
        start: usize,
        end: usize,
    ) -> Result<Vec<Block>> {
        match self.sequences.query(chr) {
            Ok(record) => {
                let mut result = Vec::new();
                for block in &record.blocks_n {
                    let block_end = block.start + block.length;
                    if block_end as usize <= start {
                        continue;
                    }
                    if block.start as usize > end {
                        break;
                    }
                    result.push(*block);
                }
                Ok(result)
            }
            Err(e) => Err(e),
        }
    }

    /// Get all soft blocks (lower case blocks) of a chromosome
    ///
    /// * `chr` Name of the chromosome from the 2bit file
    pub fn full_soft_masked_blocks(&mut self, chr: &str) -> Result<Vec<Block>> {
        match self.sequences.query(chr) {
            Ok(record) => Ok(record.blocks_soft_mask.clone()),
            Err(e) => Err(e),
        }
    }

    /// Get soft blocks (lower case blocks) of a region on a chromosome
    ///
    /// * `chr` Name of the chromosome from the 2bit file
    pub fn soft_masked_blocks(
        &mut self,
        chr: &str,
        start: usize,
        end: usize,
    ) -> Result<Vec<Block>> {
        match self.sequences.query(chr) {
            Ok(record) => {
                let mut result = Vec::new();
                for block in &record.blocks_soft_mask {
                    let block_end = block.start + block.length;
                    if block_end as usize <= start {
                        continue;
                    }
                    if block.start as usize > end {
                        break;
                    }
                    result.push(*block);
                }
                Ok(result)
            }
            Err(e) => Err(e),
        }
    }

    fn chr_length(&mut self, chr: &str) -> Result<usize> {
        match self.chroms().get(chr) {
            Some(v) => Ok(*v as usize),
            None => Err(Error::MissingName(chr.to_string())),
        }
    }

    fn read_sequence(&mut self, chr: &str, start: usize, end: usize) -> Result<String> {
        const NUC: &[u8; 4] = b"TCAG";

        let seq = self.sequences.query(chr)?;
        let reader = &mut self.reader;

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

        let seq_block = Block::new(start as _, length as _);
        replace_blocks::<true>(&mut out, seq_block, &seq.blocks_n);
        if self.softmask_enabled {
            replace_blocks::<false>(&mut out, seq_block, &seq.blocks_soft_mask);
        }

        Ok(unsafe { String::from_utf8_unchecked(out) }) // we know it's ascii so it's ok
    }
}

fn replace_blocks<const HARD: bool>(seq: &mut Vec<u8>, seq_block: Block, blocks: &[Block]) {
    let seq_block_end = seq_block.start + seq_block.length;
    for block in blocks {
        if block.start + seq_block.length <= seq_block.start {
            continue;
        }
        if seq_block_end <= block.start {
            break; // should be the last block assuming ordering is upheld
        }
        let mut range = block
            .overlap(&seq_block)
            .map_or_else(|| unsafe { core::hint::unreachable_unchecked() }, |r| r);
        range.start -= seq_block.start as usize;
        range.end -= seq_block.start as usize;
        for i in range {
            unsafe {
                *seq.get_unchecked_mut(i) = if HARD {
                    b'N'
                } else {
                    seq.get_unchecked(i).to_ascii_lowercase()
                }
            }
        }
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
    fn test_full_sequence() {
        run_test(true, |mut bit| {
            let seq = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTagctagctGATCGATCGTAGCTAGCTAGCTAGCTGATCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN".to_string();
            assert_eq!(seq, bit.full_sequence("chr1")?);
            Ok(())
        });
    }

    #[test]
    fn test_sequence() {
        run_test(true, |mut bit| {
            let seq = "NNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTagctagctGATC";
            assert_eq!(seq, bit.sequence("chr1", 24, 74)?);
            // let's try different offsets to test positions not divisible by 4
            let seq = "ACGTACGTagctagctGATC";
            assert_eq!(seq, bit.sequence("chr1", 54, 74)?);
            let seq = "CGTACGTagctagctGATC";
            assert_eq!(seq, bit.sequence("chr1", 55, 74)?);
            let seq = "GTACGTagctagctGATC";
            assert_eq!(seq, bit.sequence("chr1", 56, 74)?); // divisible by 4
            let seq = "TACGTagctagctGATC";
            assert_eq!(seq, bit.sequence("chr1", 57, 74)?);
            Ok(())
        });
    }

    #[test]
    fn test_full_bases() {
        run_test(true, |mut bit| {
            let percentages = BasePercentages {
                a: 0.08,
                c: 0.08,
                t: 0.08666666666666667,
                g: 0.08666666666666667,
                n: 100.0 / 150.0,
            };
            assert_eq!(percentages, bit.full_bases_percentages("chr1")?);
            Ok(())
        });
    }

    #[test]
    fn test_bases() {
        run_test(true, |mut bit| {
            let counts = BaseCounts {
                a: 6,
                c: 6,
                t: 6,
                g: 6,
                n: 26,
            };
            let percentages = BasePercentages {
                a: 0.12,
                c: 0.12,
                t: 0.12,
                g: 0.12,
                n: 26.0 / 50.0,
            };
            assert_eq!(counts, bit.bases("chr1", 24, 74)?);
            assert_eq!(percentages, bit.bases_percentages("chr1", 24, 74)?);
            Ok(())
        });
    }

    #[test]
    fn test_info() {
        run_test(true, |mut bit| {
            let info = TwoBitFileInfo {
                file_size: 161,
                chromosomes: 2,
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
            let mut i = 0;
            for block in bit.full_hard_masked_blocks("chr1")? {
                match i {
                    0 => {
                        assert_eq!(block.start, 0);
                        assert_eq!(block.length, 50);
                    }
                    1 => {
                        assert_eq!(block.start, 100);
                        assert_eq!(block.length, 50);
                    }
                    _ => assert!(false),
                }
                i += 1
            }
            assert_eq!(2, i);
            //TODO hard_masked_blocks()
            Ok(())
        });
    }

    #[test]
    fn test_soft_masked_blocks() {
        run_test(true, |mut bit| {
            let mut i = 0;
            for block in bit.full_soft_masked_blocks("chr1")? {
                match i {
                    0 => {
                        assert_eq!(block.start, 62);
                        assert_eq!(block.length, 8);
                    }
                    _ => assert!(false),
                }
                i += 1
            }
            assert_eq!(1, i);
            //TODO soft_masked_blocks()
            Ok(())
        });
    }

    //TODO IO Errors
}
