//! Generate 2bit files from other file formats
//!
//! Basic usage example:
//!
//! ```
//! use std::fs::File;
//! use std::io::BufWriter;
//! use twobit::convert::fasta::FastaReader;
//! use twobit::convert::to_2bit;
//!
//! let file = File::create("/tmp/my_genome.2bit").expect("doc-test");
//! let mut writer = BufWriter::new(file);
//! let fasta_sequences = FastaReader::open("assets/foo.fasta").expect("doc-test");
//! to_2bit(&mut writer, &fasta_sequences).expect("doc-test");
//! // make sure that writer goes out of scope so that the file will be closed
//!
//! ```
//!
//! The [`FastaReader`](fasta/struct.FastaReader.html) in the above example implements
//! the [`SequenceRead`](trait.SequenceRead.html) trait. If you want to convert other
//! file formats to 2bit, you can use that trait to implement your own converters.

use std::borrow::Cow;
use std::collections::HashMap;
use std::convert::TryInto;
use std::io::Write;
use std::mem::size_of;

use crate::block::Block;
use crate::error::{Error, Result};
use crate::nucleotide::Nucleotide;
use crate::reader::{Field, SIGNATURE};

pub mod fasta;

const FIELD_SIZE: usize = size_of::<Field>();

/// A simple struct, combining name and length of a sequence
#[derive(Clone)]
pub struct SequenceLength {
    pub name: String,
    pub length: usize,
}

impl SequenceLength {
    pub fn new<S: ToString>(name: &S, length: usize) -> Self {
        Self {
            name: name.to_string(),
            length,
        }
    }

    #[must_use]
    pub const fn len(&self) -> usize {
        self.length
    }

    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.length == 0
    }
}

/// Trait for readers of other file formats
///
/// This trait is used in the function [`to_2bit`](fn.to_2bit.html) to
/// retrieve the needed information to produce a 2bit file from another file format.
///
/// The trait has several methods that return a `Cow`. This `Cow` wrapper is intended to give
/// implementors the flexibility of either returning an owned or borrowed result.
/// Depending on the underlying file format, calculating results on the fly might be more effective
/// than returning a reference to a cached result.
pub trait SequenceRead<'a> {
    /// Return the names and lengths of the sequences (ideally in a stable order for reproducibility)
    fn sequence_lengths(&'a self) -> Cow<'a, [SequenceLength]>;
    /// Return pieces of the actual nucleotide sequence
    fn nucleotides(&self, chr: &str) -> Box<dyn Nucleotides>;
    /// Return the soft masked blocks for a sequence
    fn soft_masked_blocks(&'a self, chr: &str) -> Cow<'a, [Block]>;
    /// Return the hard masked blocks for a sequence
    fn hard_masked_blocks(&'a self, chr: &str) -> Cow<'a, [Block]>;
}

/// Trait for providers of nucleotide sequence data.
///
/// The data is typically provided in chunks to avoid allocating too long sequences.
pub trait Nucleotides {
    /// Read the next chunk of nucleotides from the sequence.
    /// The underlying reader will place the nucleotides into the provided
    /// buffer until the buffer is full.
    /// If the sequence is longer than the buffer, then the consumer needs to
    /// call the method again to retrieve more nucleotides.
    ///
    /// This method returns the number of Nucleotides placed into the buffer or
    /// `None` if the underlying reader has reached the end of the sequence.
    fn read_chunk(
        &mut self,
        buf: &mut [Nucleotide],
    ) -> std::result::Result<Option<usize>, Box<dyn std::error::Error>>;
}

/// Convert and write data from a `SequenceRead` into a new 2bit file.
#[allow(clippy::too_many_lines)]
pub fn to_2bit<'a>(
    writer: &mut dyn Write,
    extractor: &'a (dyn SequenceRead<'a> + 'a),
) -> Result<()> {
    let seqs = extractor.sequence_lengths();

    let signature: Field = SIGNATURE;
    let version: Field = 0;
    let sequence_count: Field = seqs.len() as Field;
    let reserved: Field = 0;

    // let's start writing immediately to force IO errors early
    for input in &[&signature, &version, &sequence_count, &reserved] {
        writer.write_all(&input.to_ne_bytes())?;
    }

    // A 2bit file has three sections: 1) header 2) index 3) sequence records.
    // In order to write the index, we first need to calculate the offsets in the sequence
    // records.

    let mut index_size: Field = 0;
    let mut sequence_records_size: Field = 0;
    let mut relative_sequence_offsets = HashMap::<&str, Field>::new(); // chr -> previous_sequence_records_sizes

    for seq_length in seqs.iter() {
        let chr = &seq_length.name;
        if !chr.is_ascii() {
            return Err(Error::FileFormat(format!(
                "Chromosome name is not ASCII: {}",
                &chr
            )));
        }

        let addition_to_index: Field = usize2field(1 + chr.len() + FIELD_SIZE)?;
        index_size += addition_to_index;
        relative_sequence_offsets.insert(chr, sequence_records_size);

        let hard_masked_blocks = extractor.hard_masked_blocks(chr);
        let soft_masked_blocks = extractor.soft_masked_blocks(chr);
        let length: usize = seq_length.length;
        sequence_records_size += {
            let block_bytewidth = size_of::<u32>() * 2; // start & end
            let integer = FIELD_SIZE // dnaSize
                          + FIELD_SIZE // nBlockCount
                          + block_bytewidth * hard_masked_blocks.len()
                          + FIELD_SIZE // maskBlockCount
                          + block_bytewidth * soft_masked_blocks.len()
                          + FIELD_SIZE; // reserved
            let condensed_length = length >> 2; // div by 4, rounded down
            let adjust_partial_quartet = if length % 4 == 0 { 0 } else { 1 };
            usize2field(integer + condensed_length + adjust_partial_quartet)?
        };
    }

    let sequence_records_start: Field = 16 + index_size; // header + index size
    for seq_length in seqs.iter() {
        let chr = &seq_length.name;
        let name = chr.as_bytes();
        let name_size: u8 = name.len().try_into().map_err(|_| {
            Error::FileFormat(format!(
                "chromosome name is longer than 255 characters: {}",
                name.len()
            ))
        })?;

        let mut buf: Vec<u8> = Vec::with_capacity(2 * FIELD_SIZE + name.len());
        buf.push(name_size);
        buf.extend_from_slice(name);
        let sequence_offset: Field = sequence_records_start
            + *relative_sequence_offsets
                .get(chr as &str)
                .expect("previous loop");
        buf.extend_from_slice(&sequence_offset.to_ne_bytes());
        writer.write_all(&buf)?;
    }

    for seq_length in seqs.iter() {
        let chr = &seq_length.name;
        let length = seq_length.length;
        let length_field = (length as Field).to_ne_bytes();
        writer.write_all(&length_field)?; // dnaSize
        write_blocks(writer, &extractor.hard_masked_blocks(chr))?;
        write_blocks(writer, &extractor.soft_masked_blocks(chr))?;
        writer.write_all(&reserved.to_ne_bytes())?; // reserved

        // and finally packedDna
        let mut data_buf = [Nucleotide::N; 80];
        let mut nuc_modulo: u8 = 0;
        let mut byte: u8 = 0;
        let mut buf = [0_u8; 1];
        let mut sequence_length = 0;
        let mut nucleotides = extractor.nucleotides(chr);
        while let Some(num_nucleotides_read) = nucleotides
            .read_chunk(&mut data_buf)
            .map_err(|e| Error::FileFormat(e.to_string()))?
        {
            for (i, nuc) in data_buf.iter().enumerate() {
                if i == num_nucleotides_read {
                    break;
                }
                byte |= nuc.bits();
                nuc_modulo += 1;
                if nuc_modulo == 4 {
                    buf[0] = byte;
                    writer.write_all(&buf)?;
                    byte = 0;
                    nuc_modulo = 0;
                } else {
                    byte <<= 2;
                }
                sequence_length += 1;
            }
        }
        if nuc_modulo > 0 {
            // if there are still 1-3 bits to be written
            byte <<= 2 * (4 - nuc_modulo); // shift all the way to the left
            buf[0] = byte;
            writer.write_all(&buf)?;
        }
        if sequence_length != length {
            let msg = format!(
                "The reported sequence length of {} for sequence {} is wrong (seen {} nucleotides)",
                length, chr, sequence_length,
            );
            return Err(Error::FileFormat(msg));
        }
    }

    Ok(())
}

/// Write `Block`s to an output stream in a 2bit-style.
///
/// First write the number of blocks, then the start positions of each block and then the end
/// positions of each block. Positions refer to the nucleotide sequence.
fn write_blocks(writer: &mut dyn Write, blocks: &[Block]) -> Result<()> {
    let num_blocks: Field = usize2field(blocks.len())?;
    writer.write_all(&num_blocks.to_ne_bytes())?;
    let mut block_lengths = Vec::with_capacity(&blocks.len() * 4); // 32 bit
    for block in blocks {
        let start: Field = usize2field(block.start)?;
        writer.write_all(&start.to_ne_bytes())?;
        let distance: Field = {
            if let Some(diff) = block.end.checked_sub(block.start) {
                usize2field(diff)?
            } else {
                return Err(Error::FileFormat("Block end < Block start".to_string()));
            }
        };
        //let distance: Field = block.end - block.start;
        block_lengths.extend_from_slice(&distance.to_ne_bytes());
    }
    writer.write_all(&block_lengths)?;
    Ok(())
}

fn usize2field(v: usize) -> Result<Field> {
    v.try_into()
        .map_err(|_| Error::FileFormat(format!("Number {} is too big to store in a 2bit Field", v)))
}
