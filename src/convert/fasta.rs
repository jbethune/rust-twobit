//! Convert twobit files to Fasta and vice-versa
//!
//! This modules assumes Fasta files to be simple: No comments, no UIPAC stuff.
//! Just `ACGTacgtN` and `>header lines`. Soft blocks (lower case nucleotides)
//! will be automatically detected.

use std::collections::HashMap;
use std::convert::TryInto;
use std::fs::File;
use std::io::{BufReader, Cursor, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use crate::block::Block;
use crate::convert::{Nucleotides, SequenceLength, SequenceRead};
use crate::error::Error;
use crate::nucleotide::Nucleotide;
use crate::TwoBitFile;

type FilePos = u64;
trait FileAccess: Read + Seek {}

impl FileAccess for BufReader<File> {}
impl FileAccess for Cursor<Vec<u8>> {}

/// Read a Fasta file and make the extracted information available
pub struct FastaReader {
    file_path: Option<PathBuf>,
    in_memory_data: Option<Cursor<Vec<u8>>>,
    soft_blocks: HashMap<String, Vec<Block>>,
    hard_blocks: HashMap<String, Vec<Block>>,
    sequence_starts: HashMap<String, FilePos>,
    sequence_lengths: Vec<SequenceLength>,
    empty_block_vec: Vec<Block>,
}

impl FastaReader {
    /// Read data from a Fasta file that is stored in a `Vec<u8>`.
    pub fn mem_open(data: Vec<u8>) -> Result<Self, Error> {
        let mut reader = Cursor::new(data);
        let mut result = Self::parse_fasta_file(&mut reader)?;
        result.in_memory_data = Some(reader);
        Ok(result)
    }

    /// Open a Fasta file on disk
    pub fn open<P: AsRef<Path>>(file_path: P) -> Result<Self, Error> {
        let fd = File::open(&file_path)?;
        let mut reader = BufReader::new(fd);
        let mut result = Self::parse_fasta_file(&mut reader)?;
        result.file_path = Some(PathBuf::from(file_path.as_ref()));
        Ok(result)
    }

    /// Extract all the necessary information from the Fasta file
    #[allow(clippy::too_many_lines)]
    fn parse_fasta_file<R: Read>(reader: &mut R) -> Result<Self, Error> {
        let mut in_seq = false;
        let mut seen_newline = true;
        let mut name_buf = String::with_capacity(64);
        let mut chr = String::new();
        let mut hard_blocks = HashMap::<String, Vec<Block>>::new();
        let mut hard_block_start: Option<usize> = None;
        let mut soft_blocks = HashMap::<String, Vec<Block>>::new();
        let mut soft_block_start: Option<usize> = None;
        let mut sequence_starts = HashMap::new();
        let mut sequence_lengths = vec![];
        let mut terminate_any_soft_block = false;
        let mut terminate_any_hard_block = false;
        let mut data = [0_u8; 256];
        let mut pos: FilePos = 0;
        let mut in_seq_pos = 0_usize;
        loop {
            let count = reader.read(&mut data)?;
            if count == 0 {
                break;
            }
            for byte in data.iter().take(count) {
                if in_seq {
                    match byte {
                        b'\n' => {
                            seen_newline = true;
                            pos += 1; // need to count each byte
                            continue;
                        }
                        b'>' => {
                            if seen_newline {
                                in_seq = false;
                                sequence_lengths.push(SequenceLength::new(&chr, in_seq_pos));
                                terminate_any_soft_block = true;
                                terminate_any_hard_block = true;
                            } else {
                                return Err(Error::FileFormat(format!(
                                    "Unexpected \">\" in Fasat file at position {}",
                                    pos
                                )));
                            }
                        }
                        b'A' | b'C' | b'G' | b'T' | b'U' => {
                            if soft_block_start.is_some() || hard_block_start.is_some() {
                                terminate_any_soft_block = true;
                                terminate_any_hard_block = true;
                            }
                        }
                        b'a' | b'c' | b'g' | b't' | b'u' => {
                            if soft_block_start.is_none() {
                                soft_block_start = Some(in_seq_pos);
                            }
                            if hard_block_start.is_none() {
                                terminate_any_hard_block = true;
                            }
                        }
                        b'N' => {
                            if hard_block_start.is_none() {
                                hard_block_start = Some(in_seq_pos);
                            }
                            if soft_block_start.is_some() {
                                terminate_any_soft_block = true;
                            }
                        }
                        _ => {
                            return Err(Error::FileFormat(format!(
                                "Bad symbol in FASTA file: {}",
                                byte
                            )));
                        }
                    }
                    if terminate_any_soft_block {
                        if let Some(start) = soft_block_start {
                            let block = (start as usize)..in_seq_pos;
                            let soft_blocks_entry =
                                soft_blocks.entry(chr.to_string()).or_insert_with(Vec::new);
                            soft_blocks_entry.push(block);
                            soft_block_start = None;
                        }
                        terminate_any_soft_block = false;
                    }
                    if terminate_any_hard_block {
                        if let Some(start) = hard_block_start {
                            let block = (start as usize)..in_seq_pos;
                            let hard_blocks_entry =
                                hard_blocks.entry(chr.to_string()).or_insert_with(Vec::new);
                            hard_blocks_entry.push(block);
                            hard_block_start = None;
                        }
                        terminate_any_hard_block = false;
                    }
                    seen_newline = false;
                    in_seq_pos += 1;
                } else if *byte == b'\n' {
                    in_seq = true;
                    in_seq_pos = 0;
                    sequence_starts.insert(name_buf.clone(), pos + 1);
                    chr = name_buf.clone();
                    name_buf.clear();
                } else if pos > 0 {
                    name_buf.push(*byte as char);
                } else if *byte != b'>' {
                    return Err(Error::FileFormat(
                        "Fasta file does not start with a \">\" symbol".to_string(),
                    ));
                }
                pos += 1;
            }
        }

        //finish off last seq
        if let Some(start) = soft_block_start {
            let block = (start as usize)..in_seq_pos;
            let soft_blocks_entry = soft_blocks.entry(chr.to_string()).or_insert_with(Vec::new);
            soft_blocks_entry.push(block);
        }
        if let Some(start) = hard_block_start {
            let block = (start as usize)..in_seq_pos;
            let hard_blocks_entry = hard_blocks.entry(chr.to_string()).or_insert_with(Vec::new);
            hard_blocks_entry.push(block);
        }
        sequence_lengths.push(SequenceLength::new(&chr, in_seq_pos));
        Ok(Self {
            file_path: None,
            in_memory_data: None,
            soft_blocks,
            hard_blocks,
            sequence_starts,
            sequence_lengths,
            empty_block_vec: vec![],
        })
    }
}

/// Gives access to nucleotide sequences from a Fasta file
struct FastaNucleotides {
    reader: Box<dyn FileAccess>,
    exhausted: bool,
}

impl Nucleotides for FastaNucleotides {
    fn read_chunk(
        &mut self,
        buf: &mut [Nucleotide],
    ) -> std::result::Result<Option<usize>, Box<dyn std::error::Error>> {
        if self.exhausted {
            Ok(None)
        } else {
            let mut newlines = 0;
            let mut values = vec![0_u8; buf.len()];
            let count = self.reader.read(&mut values)?;
            if count == 0 {
                self.exhausted = true;
                return Ok(None);
            }
            for (i, byte) in values.iter().take(count).enumerate() {
                if *byte == b'\n' {
                    newlines += 1;
                } else if *byte == b'>' {
                    // we are running into the next sequence
                    self.exhausted = true;
                    return Ok(Some(i - newlines));
                } else {
                    let nuc: Nucleotide = (*byte).try_into()?;
                    buf[i - newlines] = nuc;
                }
            }
            Ok(Some(count - newlines))
        }
    }
}

impl<'a> SequenceRead<'a> for FastaReader {
    fn sequence_lengths(&'a self) -> Result<&'a [SequenceLength], Box<dyn std::error::Error>> {
        Ok(&self.sequence_lengths)
    }

    fn nucleotides(&self, chr: &str) -> Result<Box<dyn Nucleotides>, Box<dyn std::error::Error>> {
        let mut reader: Box<dyn FileAccess> = {
            if let Some(file_path) = &self.file_path {
                Box::new(BufReader::new(File::open(file_path)?))
            } else if let Some(cursor) = &self.in_memory_data {
                Box::new(cursor.clone())
            } else {
                panic!("Either one of the two constructors must have been used");
            }
        };
        if let Some(seek_to) = self.sequence_starts.get(chr) {
            reader.seek(SeekFrom::Start(*seek_to))?;
            Ok(Box::new(FastaNucleotides {
                reader,
                exhausted: false,
            }))
        } else {
            // no location information about the sequence
            Ok(Box::new(FastaNucleotides {
                reader,
                exhausted: true,
            }))
        }
    }

    fn soft_masked_blocks(&'a self, chr: &str) -> Result<&'a [Block], Box<dyn std::error::Error>> {
        match self.soft_blocks.get(chr) {
            Some(blocks) => Ok(blocks),
            None => Ok(&self.empty_block_vec),
        }
    }

    fn hard_masked_blocks(&'a self, chr: &str) -> Result<&'a [Block], Box<dyn std::error::Error>> {
        match self.hard_blocks.get(chr) {
            Some(blocks) => Ok(blocks),
            None => Ok(&self.empty_block_vec),
        }
    }
}

/// Convert a 2bit file object to a Fasta file
///
/// Remember to call `twobit.enable_softmask()` if you want softmasked
/// nucleotides in your output.
pub fn to_fasta<R>(twobit: &mut TwoBitFile<R>, writer: &mut dyn Write) -> Result<(), Error>
where
    R: Read + Seek,
{
    const LINE_WIDTH: usize = 80;

    for name in twobit.chrom_names() {
        writer.write_all(&[b'>'])?;
        writer.write_all(name.as_bytes())?;
        writer.write_all(&[b'\n'])?;

        let mut start = 0;
        loop {
            let stop = start + LINE_WIDTH;
            let seq = twobit.read_sequence(&name, start..stop)?;
            if seq.is_empty() {
                break;
            }
            writer.write_all(seq.as_bytes())?;
            writer.write_all(&[b'\n'])?;
            start = stop;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    const FASTA_FILE_DATA: &[u8] =
        b">name1\nNNNNAGTCGTcagtcGTCGTAGNNNNNNTCTACGTATgcgtcaNNNN\n>name2\nCATGCA\nACGTACGCAT";

    #[test]
    fn test_parse_fasta_file() {
        let data: Vec<u8> = Vec::from(FASTA_FILE_DATA);
        let fasta_file = FastaReader::mem_open(data).expect("unit-test");
        assert!(fasta_file.file_path.is_none());
        assert!(fasta_file.in_memory_data.is_some());
        assert_eq!(fasta_file.sequence_starts.len(), 2);
        assert_eq!(
            fasta_file.sequence_starts.get("name1").expect("unit-test"),
            &7
        );
        assert_eq!(
            fasta_file.sequence_starts.get("name2").expect("unit-test"),
            &62
        );
    }

    #[test]
    fn test_sequence_lengths() {
        let data: Vec<u8> = Vec::from(FASTA_FILE_DATA);
        let fasta_file = FastaReader::mem_open(data).expect("unit-test");
        let lengths = fasta_file.sequence_lengths().expect("unit-test");
        assert_eq!(lengths.len(), 2);
        assert_eq!(lengths[0].name, "name1");
        assert_eq!(lengths[0].length, 47);
        assert_eq!(lengths[1].name, "name2");
        assert_eq!(lengths[1].length, 16);
    }

    #[test]
    fn test_block_parsing() {
        let data: Vec<u8> = Vec::from(FASTA_FILE_DATA);
        let fasta_file = FastaReader::mem_open(data).expect("unit-test");

        fn check_soft_blocks(blocks: &[Block]) {
            assert_eq!(blocks.len(), 2);
            for (i, block) in blocks.iter().enumerate() {
                match i {
                    0 => assert_eq!(*block, 10usize..15usize),
                    1 => assert_eq!(*block, 37usize..43usize),
                    _ => assert!(false),
                }
            }
        }

        fn check_hard_blocks(blocks: &[Block]) {
            assert_eq!(blocks.len(), 3);
            for (i, block) in blocks.iter().enumerate() {
                match i {
                    0 => assert_eq!(*block, 0usize..4usize),
                    1 => assert_eq!(*block, 22usize..28usize),
                    2 => assert_eq!(*block, 43usize..47usize),
                    _ => assert!(false),
                }
            }
        }

        check_soft_blocks(fasta_file.soft_masked_blocks("name1").expect("unit-test"));
        assert_eq!(
            fasta_file
                .soft_masked_blocks("name2")
                .expect("unit-test")
                .len(),
            0
        );

        check_hard_blocks(fasta_file.hard_masked_blocks("name1").expect("unit-test"));
        assert_eq!(
            fasta_file
                .hard_masked_blocks("name2")
                .expect("unit-test")
                .len(),
            0
        );
    }
}
