//! Extraction of data from 2bit files

use std::fs::File;
use std::io::{BufReader, Cursor, Read, Seek, SeekFrom};
use std::mem::size_of;
use std::path::Path;

use crate::block::Blocks;
use crate::error::{Error, Result};
use crate::{SequenceRecord, SequenceRecords};

// 2bit field type (uint32)
pub type Field = u32;
/// 2bit field size (4 bytes)
const FIELD_SIZE: usize = size_of::<Field>();
/// 2bit signature magic number
const SIGNATURE: Field = 0x1A41_2743;
/// 2bit signature magic number reversed
const REV_SIGNATURE: Field = 0x4327_411A;

pub trait Reader: Read + Seek {}

impl<T: Read + Seek> Reader for T {}

/// Extract binary data from 2bit files
///
/// This reads all types of fields except the sequences.
pub struct ValueReader<R: Reader> {
    reader: R,
    twobit_version: Field,
    swap_endian: bool,
}

pub type BoxValueReader = ValueReader<Box<dyn Reader>>;

impl ValueReader<BufReader<File>> {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::new(BufReader::new(File::open(path)?))
    }
}

impl<T> ValueReader<Cursor<T>>
where
    Cursor<T>: Read + Seek,
{
    pub fn from_buf(buf: T) -> Result<Self> {
        Self::new(Cursor::new(buf))
    }
}

impl ValueReader<Cursor<Vec<u8>>> {
    pub fn open_and_read<P: AsRef<Path>>(path: P) -> Result<Self> {
        let mut buf = vec![];
        File::open(path)?.read_to_end(&mut buf)?;
        Self::from_buf(buf)
    }
}

impl<R: Reader> ValueReader<R> {
    pub fn new(reader: R) -> Result<Self> {
        let mut result = Self {
            reader,
            twobit_version: 0,
            swap_endian: false,
        };
        let signature = result.field()?;
        if signature != SIGNATURE {
            if signature == REV_SIGNATURE {
                result.swap_endian = true;
            } else {
                return Err(Error::FileFormat(
                    "File does not start with 2bit signature".to_string(),
                ));
            }
        }
        let version = result.field()?;
        if version == 0 {
            result.twobit_version = version;
            Ok(result)
        } else {
            Err(Error::UnsupportedVersion(
                "Versions larger than 0 are not supported".to_string(),
            ))
        }
    }

    pub fn boxed(self) -> BoxValueReader
    where
        R: 'static,
    {
        ValueReader {
            reader: Box::new(self.reader),
            twobit_version: self.twobit_version,
            swap_endian: self.swap_endian,
        }
    }

    #[inline]
    pub fn seek(&mut self, pos: SeekFrom) -> Result<u64> {
        Ok(self.reader.seek(pos)?)
    }

    #[inline]
    pub fn byte(&mut self) -> Result<u8> {
        let mut buf: [u8; 1] = [0; 1];
        self.reader.read_exact(&mut buf)?;
        Ok(buf[0])
    }

    pub fn stream_len(&mut self) -> Result<u64> {
        // borrowed from unstable Seek method in stdlib
        let old_pos = self.reader.stream_position()?;
        let len = self.reader.seek(SeekFrom::End(0))?;
        // Avoid seeking a third time when we were already at the end of the
        // stream. The branch is usually way cheaper than a seek operation.
        if old_pos != len {
            self.reader.seek(SeekFrom::Start(old_pos))?;
        }
        Ok(len)
    }

    fn tell(&mut self) -> Result<u64> {
        Ok(self.reader.seek(SeekFrom::Current(0))?)
    }

    fn field(&mut self) -> Result<Field> {
        let mut field: [u8; FIELD_SIZE] = [0; FIELD_SIZE];
        self.reader.read_exact(&mut field)?;
        Ok(slice_to_field(field, self.swap_endian))
    }

    fn fields(&mut self, n: usize) -> Result<Vec<Field>> {
        (0..n).map(|_| self.field()).collect()
    }

    fn string(&mut self, length: usize) -> Result<String> {
        let mut buf = vec![0_u8; length];
        self.reader.read_exact(&mut buf)?;
        Ok(String::from_utf8(buf)?)
    }

    fn blocks(&mut self) -> Result<Blocks> {
        let n = self.field()? as usize;
        Ok(Blocks(
            self.fields(n)?
                .iter()
                .zip(&self.fields(n)?)
                .map(|(&s, &l)| (s as usize)..((s + l) as usize))
                .collect(),
        ))
    }

    fn sequence_record(&mut self) -> Result<SequenceRecord> {
        let chr_size = self.byte()? as usize;
        let chr = self.string(chr_size)?;
        let seq_offset = u64::from(self.field()?);
        let backup_offset = self.tell()?;
        self.seek(SeekFrom::Start(seq_offset))?;
        let length = self.field()? as usize;
        let blocks_n = self.blocks()?;
        let blocks_soft_mask = self.blocks()?;
        let _reserved = self.field()?;
        let offset = self.tell()?;
        self.seek(SeekFrom::Start(backup_offset))?;

        Ok(SequenceRecord {
            chr,
            offset,
            length,
            blocks_n,
            blocks_soft_mask,
        })
    }

    pub(crate) fn sequence_records(&mut self) -> Result<SequenceRecords> {
        self.seek(SeekFrom::Start(2 * size_of::<Field>() as u64))?; // rewind to start of file
        let n_seq = self.field()?;
        let _reserved = self.field(); // read the unused reserved field

        Ok(SequenceRecords(
            (0..n_seq)
                .map(|_| self.sequence_record())
                .collect::<Result<_>>()?,
        ))
    }
}

fn slice_to_field(slice: [u8; FIELD_SIZE], swap_endian: bool) -> Field {
    let mut result = 0;
    if swap_endian {
        for byte in slice.iter().rev().copied() {
            result <<= 8;
            result += Field::from(byte);
        }
    } else {
        for byte in slice {
            result <<= 8;
            result += Field::from(byte);
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_slice_to_field() {
        let slice: [u8; 4] = [0x1A, 0x41, 0x27, 0x43];
        assert_eq!(slice_to_field(slice, false), 0x1A412743);
        let slice: [u8; 4] = [0x43, 0x27, 0x41, 0x1A];
        assert_eq!(slice_to_field(slice, false), 0x4327411A);

        assert_eq!(
            slice_to_field([0x1A, 0x41, 0x27, 0x43], true),
            slice_to_field([0x43, 0x27, 0x41, 0x1A,], false)
        );

        let slice: [u8; 4] = [0, 2, 3, 4];
        assert_eq!(slice_to_field(slice, false), 2 * 65536 + 3 * 256 + 4);
    }
}
