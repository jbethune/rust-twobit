//! Read 2bit files in Rust.
//!
//! This crate is inspired by <a href="https://github.com/deeptools/py2bit">py2bit</a> and tries to
//! offer the same functionality. It is written from the ground up with no C-dependency. It follows
//! the <a href="http://genome.ucsc.edu/FAQ/FAQformat.html#format7">2bit specification version
//! 0</a>.
//!
//! Note that most methods perform IO operations and return a `Result` for that reason.
//!
//! The main entry point for users of this crate is the [TwoBitFile](struct.TwoBitFile.html)
//! struct. Please see its documentation for details how to use this crate.
//!

pub mod block;
pub mod counts;
pub mod error;
pub mod types;
mod value_reader;

use std::collections::HashMap;
use std::default::Default;
use std::fs::metadata;
use std::io::SeekFrom;
use std::path::{Path, PathBuf};

use crate::block::Block;
use crate::counts::{BaseCounts, BasePercentages};
use crate::error::Error;
use crate::types::{Field, FileIndex};
use crate::value_reader::ValueReader;

/// 2bit signature magic number
const SIGNATURE: Field = 0x1A412743;

/// 2bit signature magic number reversed
const REV_SIGNATURE: Field = 0x4327411A;

/// Read data from a 2bit file
///
/// Usage:
///
/// ```
/// extern crate twobit;
/// use twobit::TwoBitFile;
///
/// let enable_softmask=true; //set to false to enforce upper case sequences
/// let tb = TwoBitFile::open("assets/foo.2bit", enable_softmask).unwrap();
/// let chromosome_lengths = tb.chroms();
/// assert_eq!(chromosome_lengths["chr1"], 150);
/// assert_eq!(chromosome_lengths["chr2"], 100);
/// ```
///
/// 2bit files offer two types of masks: N masks (aka hard masks) for unknown or arbitrary nucleotides
/// and soft masks for lower-case nucleotides (e.g. "t" instead of "T").
///
/// ```
/// # extern crate twobit;
/// # use twobit::TwoBitFile;
/// # let enable_softmask=true; //set to false to enforce upper case sequences
/// # let tb = TwoBitFile::open("assets/foo.2bit", enable_softmask).unwrap();
///
/// let expected_seq = "NNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTagctagctGATC"; // some lower and some upper case
/// assert_eq!(tb.sequence("chr1",24,74).unwrap(), expected_seq);
/// ```
/// It is not possible to disable hard masks but you can disable soft masks:
/// ```
/// # extern crate twobit;
/// # use twobit::TwoBitFile;
///
/// let enable_softmask=false;
/// let tb_nosoft = TwoBitFile::open("assets/foo.2bit", enable_softmask).unwrap();
/// let expected_seq = "NNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTAGCTAGCTGATC"; // all upper case
/// assert_eq!(tb_nosoft.sequence("chr1",24,74).unwrap(), expected_seq);
/// ```
///
/// There are 2 variants for every sequence-related method:
/// * Those that request information on parts of a chromosome
/// * Those that request information on a complete chromosome
///
/// The partial requests require you to provide a start (0-based, inclusive) and a stop (0-based,
/// exclusive) position as parameters. The full request methods all start with the prefix
/// `full_`.
pub struct TwoBitFile {
    path: PathBuf,
    softmask_enabled: bool,
    sequences: HashMap<String, FileIndex>,
}

struct SequenceRecord {
    dna_offset: FileIndex,
    dna_size: Field,
    n_blocks: Vec<Block>,
    soft_mask_blocks: Vec<Block>,
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

impl TwoBitFile {
    /// Constructor for a 2bit file.
    ///
    /// * `path` - A path to the 2bit file
    /// * `softmask_enabled` - return lower case nucleotides for soft blocks
    pub fn open<P: AsRef<Path>>(path: P, softmask_enabled: bool) -> Result<TwoBitFile, Error> {
        let mut reader = ValueReader::from_path(&path)?;

        let mut sequences = HashMap::new();

        let sequence_count = reader.field()?;
        let _reserved = reader.field(); // read the unused reserved field

        for _ in 0..sequence_count {
            let name_size = reader.byte()? as usize;
            let name = reader.string(name_size)?;
            let seq_offset = reader.field()?;

            sequences.insert(name, seq_offset);
        }

        let path = path.as_ref().to_path_buf();
        Ok(TwoBitFile {
            path,
            sequences,
            softmask_enabled,
        })
    }

    /// Get the sizes of chromosomes in a 2bit file as a `HashMap`
    pub fn chroms(&self) -> HashMap<String, Field> {
        let mut result = HashMap::new();
        for chr in self.sequences.keys() {
            let seq = self.sequence_record(chr).expect("Chromosome must be valid");
            result.insert((*chr).clone(), seq.dna_size);
        }
        result
    }

    /// Get the full sequence of a chromosome
    ///
    /// * `chr` Name of the chromosome from the 2bit file
    /// * returns the full sequence as a `String` on success
    pub fn full_sequence(&self, chr: &str) -> Result<String, Error> {
        if self.sequences.contains_key(chr) {
            let record = self.sequence_record(chr)?;
            record.sequence(&self.path, 0, record.dna_size as usize)
        } else {
            Err(Error::MissingName(chr.to_string()))
        }
    }

    /// Get a partial sequence of a chromosome
    ///
    /// * `chr` Name of the chromosome from the 2bit file
    /// * `start` Start position of the sequence
    /// * `end` Stop position of the sequence
    /// * returns the full sequence as a `String` on success
    pub fn sequence(&self, chr: &str, start: usize, end: usize) -> Result<String, Error> {
        if self.sequences.contains_key(chr) {
            let record = self.sequence_record(chr)?;
            record.sequence(&self.path, start, end)
        } else {
            Err(Error::MissingName(chr.to_string()))
        }
    }

    /// Count bases of a chromosome
    ///
    /// * `chr` Name of the chromosome from the 2bit file
    pub fn full_bases(&self, chr: &str) -> Result<BaseCounts, Error> {
        self.bases(chr, 0, self.chr_length(chr)?)
    }

    /// Calculates percentages of bases of a chromosome
    ///
    /// * `chr` Name of the chromosome from the 2bit file
    pub fn full_bases_percentages(&self, chr: &str) -> Result<BasePercentages, Error> {
        self.bases_percentages(chr, 0, self.chr_length(chr)?)
    }

    /// Count bases of a partial chromosome
    ///
    /// * `chr` Name of the chromosome from the 2bit file
    pub fn bases(&self, chr: &str, start: usize, end: usize) -> Result<BaseCounts, Error> {
        let nucs = self
            .sequence_record(chr)?
            .sequence(self.path.clone(), start, end)?;
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
        &self,
        chr: &str,
        start: usize,
        end: usize,
    ) -> Result<BasePercentages, Error> {
        Ok(self.bases(chr, start, end)?.into())
    }

    /// Obtain general information on the 2bit file you are dealing with
    pub fn info(&self) -> Result<TwoBitFileInfo, Error> {
        let mut total_length = 0;
        let mut hard_masks_count = 0;
        let mut soft_masks_count = 0;
        for chr in self.sequences.keys() {
            total_length += self.chr_length(chr)?;
            let record = self.sequence_record(chr)?;
            hard_masks_count +=
                record.n_blocks.iter().fold(0, |acc, blk| acc + blk.length) as usize;
            soft_masks_count += record
                .soft_mask_blocks
                .iter()
                .fold(0, |acc, blk| acc + blk.length) as usize;
        }
        Ok(TwoBitFileInfo {
            file_size: metadata(&self.path)?.len(),
            chromosomes: self.sequences.len(),
            total_sequence_length: total_length,
            hard_masks_count,
            soft_masks_count,
        })
    }

    /// Get all hard blocks (N-blocks) of a chromosome
    ///
    /// * `chr` Name of the chromosome from the 2bit file
    pub fn full_hard_masked_blocks(&self, chr: &str) -> Result<Vec<Block>, Error> {
        match self.sequence_record(chr) {
            Ok(record) => Ok(record.n_blocks),
            Err(e) => Err(e),
        }
    }

    /// Get hard blocks (N-blocks) of a region on a chromosome
    ///
    /// * `chr` Name of the chromosome from the 2bit file
    pub fn hard_masked_blocks(
        &self,
        chr: &str,
        start: usize,
        end: usize,
    ) -> Result<Vec<Block>, Error> {
        match self.sequence_record(chr) {
            Ok(record) => {
                let mut result = Vec::new();
                for block in record.n_blocks {
                    let block_end = block.start + block.length;
                    if block_end as usize <= start {
                        continue;
                    } else if block.start as usize > end {
                        break;
                    }
                    result.push(block)
                }
                Ok(result)
            }
            Err(e) => Err(e),
        }
    }

    /// Get all soft blocks (lower case blocks) of a chromosome
    ///
    /// * `chr` Name of the chromosome from the 2bit file
    pub fn full_soft_masked_blocks(&self, chr: &str) -> Result<Vec<Block>, Error> {
        match self.sequence_record(chr) {
            Ok(record) => Ok(record.soft_mask_blocks),
            Err(e) => Err(e),
        }
    }

    /// Get soft blocks (lower case blocks) of a region on a chromosome
    ///
    /// * `chr` Name of the chromosome from the 2bit file
    pub fn soft_masked_blocks(
        &self,
        chr: &str,
        start: usize,
        end: usize,
    ) -> Result<Vec<Block>, Error> {
        match self.sequence_record(chr) {
            Ok(record) => {
                let mut result = Vec::new();
                for block in record.soft_mask_blocks {
                    let block_end = block.start + block.length;
                    if block_end as usize <= start {
                        continue;
                    } else if block.start as usize > end {
                        break;
                    }
                    result.push(block)
                }
                Ok(result)
            }
            Err(e) => Err(e),
        }
    }

    fn sequence_record(&self, chr: &str) -> Result<SequenceRecord, Error> {
        let mut reader = ValueReader::from_path(&self.path)?;
        let offset = self
            .sequences
            .get(chr)
            .ok_or_else(|| Error::MissingName(chr.to_string()))?;
        reader.seek(SeekFrom::Start(*offset as u64))?;
        let dna_size = reader.field()?;
        let n_blocks = reader.blocks()?;
        let soft_mask_blocks = {
            if self.softmask_enabled {
                reader.blocks()?
            } else {
                reader.skip_blocks()?;
                Vec::new()
            }
        };
        let _reserved = reader.field()?;
        let dna_offset = reader.tell()?;

        Ok(SequenceRecord {
            dna_offset,
            dna_size,
            n_blocks,
            soft_mask_blocks,
        })
    }

    fn chr_length(&self, chr: &str) -> Result<usize, Error> {
        match self.chroms().get(chr) {
            Some(v) => Ok(*v as usize),
            None => Err(Error::MissingName(chr.to_string())),
        }
    }
}

impl SequenceRecord {
    fn sequence<P: AsRef<Path>>(&self, path: P, start: usize, end: usize) -> Result<String, Error> {
        let mut reader = ValueReader::from_path(path)?;
        let mut skip = start % 4;
        reader
            .seek(SeekFrom::Start(self.dna_offset as u64))
            .unwrap(); // beginning of the DNA sequence
        reader.seek(SeekFrom::Current(start as i64 / 4))?; // position where we want to start reading

        let length = end - start;

        let mut result = String::with_capacity(length);
        while result.len() < length {
            let mut byte = reader.byte()?;
            for _ in 0..4 {
                if skip > 0 {
                    skip -= 1;
                    byte <<= 2;
                    continue;
                }
                let nuc = match byte & 192 {
                    0 => 'T',
                    64 => 'C',
                    128 => 'A',
                    192 => 'G',
                    _ => panic!("Bit flag magic is broken"),
                };
                byte <<= 2;
                result.push(nuc);
                if result.len() == length {
                    break;
                }
            }
        }

        let seq_block = Block::new(start as Field, length as Field);

        // closure over `seq_block`, `start` and `end`
        let replace_blocks =
            |mut seq: String, blocks: &Vec<Block>, replacement_mode: BlockReplacementMode| {
                for block in blocks {
                    let block_end = (block.start + block.length) as usize;
                    if block_end <= start {
                        continue;
                    } else if end <= block.start as usize {
                        break; // should be the last block assuming ordering is upheld
                    }
                    let mut range = block
                        .overlap(&seq_block)
                        .expect("Block ordering not upheld");

                    // adjust for partial sequence
                    range.start -= start;
                    range.end -= start;

                    let replacement = match replacement_mode {
                        BlockReplacementMode::Hard => "N".repeat(range.end - range.start),
                        BlockReplacementMode::Soft => seq[range.clone()].to_lowercase(),
                    };
                    seq.replace_range(range, &replacement);
                }
                seq
            };
        let result = replace_blocks(result, &self.n_blocks, BlockReplacementMode::Hard);
        let result = replace_blocks(result, &self.soft_mask_blocks, BlockReplacementMode::Soft);

        Ok(result)
    }
}

enum BlockReplacementMode {
    Soft, // convert to lower case
    Hard, // replace with Ns
}

#[cfg(test)]
mod tests {
    use super::*;
    const TESTFILE: &str = "assets/foo.2bit";

    #[test]
    fn test_chroms() {
        let bit = TwoBitFile::open(TESTFILE, true).unwrap();
        let mut chr_count = 0;
        for (chr, len) in bit.chroms() {
            match chr.as_ref() {
                "chr1" => assert_eq!(len, 150),
                "chr2" => assert_eq!(len, 100),
                _ => assert!(false), // unexpected chromosome
            }
            chr_count += 1;
        }
        assert_eq!(chr_count, 2)
    }

    #[test]
    fn test_full_sequence() {
        let seq = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTagctagctGATCGATCGTAGCTAGCTAGCTAGCTGATCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN".to_string();
        let bit = TwoBitFile::open(TESTFILE, true).unwrap();
        assert_eq!(seq, bit.full_sequence("chr1").unwrap());
    }

    #[test]
    fn test_sequence() {
        let seq = "NNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTagctagctGATC";
        let bit = TwoBitFile::open(TESTFILE, true).unwrap();
        assert_eq!(seq, bit.sequence("chr1", 24, 74).unwrap());
        // let's try different offsets to test positions not divisible by 4
        let seq = "ACGTACGTagctagctGATC";
        assert_eq!(seq, bit.sequence("chr1", 54, 74).unwrap());
        let seq = "CGTACGTagctagctGATC";
        assert_eq!(seq, bit.sequence("chr1", 55, 74).unwrap());
        let seq = "GTACGTagctagctGATC";
        assert_eq!(seq, bit.sequence("chr1", 56, 74).unwrap()); // divisible by 4
        let seq = "TACGTagctagctGATC";
        assert_eq!(seq, bit.sequence("chr1", 57, 74).unwrap());
    }

    #[test]
    fn test_full_bases() {
        let percentages = BasePercentages {
            a: 0.08,
            c: 0.08,
            t: 0.08666666666666667,
            g: 0.08666666666666667,
            n: 100.0 / 150.0,
        };
        let bit = TwoBitFile::open(TESTFILE, true).unwrap();
        assert_eq!(percentages, bit.full_bases_percentages("chr1").unwrap());
    }

    #[test]
    fn test_bases() {
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
        let bit = TwoBitFile::open(TESTFILE, true).unwrap();
        assert_eq!(counts, bit.bases("chr1", 24, 74).unwrap());
        assert_eq!(percentages, bit.bases_percentages("chr1", 24, 74).unwrap());
    }

    #[test]
    fn test_info() {
        let info = TwoBitFileInfo {
            file_size: 161,
            chromosomes: 2,
            total_sequence_length: 250,
            hard_masks_count: 150,
            soft_masks_count: 8,
        };
        let bit = TwoBitFile::open(TESTFILE, true).unwrap();
        assert_eq!(bit.info().unwrap(), info);
    }

    #[test]
    fn test_hard_masked_blocks() {
        let bit = TwoBitFile::open(TESTFILE, true).unwrap();
        let mut i = 0;
        for block_ in bit.full_hard_masked_blocks("chr1").unwrap() {
            match i {
                0 => {
                    assert_eq!(block_.start, 0);
                    assert_eq!(block_.length, 50);
                }
                1 => {
                    assert_eq!(block_.start, 100);
                    assert_eq!(block_.length, 50);
                }
                _ => assert!(false),
            }
            i += 1
        }
        assert_eq!(2, i);
        //TODO hard_masked_blocks()
    }

    #[test]
    fn test_soft_masked_blocks() {
        let bit = TwoBitFile::open(TESTFILE, true).unwrap();
        let mut i = 0;
        for block_ in bit.full_soft_masked_blocks("chr1").unwrap() {
            match i {
                0 => {
                    assert_eq!(block_.start, 62);
                    assert_eq!(block_.length, 8);
                }
                _ => assert!(false),
            }
            i += 1
        }
        assert_eq!(1, i);
        //TODO soft_masked_blocks()
    }

    //TODO IO Errors
}
