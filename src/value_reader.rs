//! Extraction of data from 2bit files

use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::mem::size_of;
use std::path::Path;

use crate::block::Block;
use crate::error::Error;
use crate::types::{Field, FileIndex};
use crate::{REV_SIGNATURE, SIGNATURE};

/// Extract binary data from 2bit files
///
/// This reads all types of fields except the sequences.
pub struct ValueReader<R: Read + Seek> {
    reader: R,
    twobit_version: Field,
    swap_endian: bool,
}

impl ValueReader<BufReader<File>> {
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, Error> {
        Self::new(BufReader::new(File::open(path)?))
    }
}

impl<R: Read + Seek> ValueReader<R> {
    pub fn new(reader: R) -> Result<Self, Error> {
        let mut result = ValueReader {
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

    pub fn seek(&mut self, pos: SeekFrom) -> Result<FileIndex, Error> {
        match self.reader.seek(pos) {
            Ok(v) => Ok(v as FileIndex),
            Err(e) => Err(e.into()),
        }
    }

    pub fn seek_start(&mut self) -> Result<(), Error> {
        self.seek(SeekFrom::Start(2 * size_of::<Field>() as u64))?;
        Ok(())
    }

    pub fn tell(&mut self) -> Result<FileIndex, Error> {
        match self.reader.seek(SeekFrom::Current(0)) {
            Ok(v) => Ok(v as FileIndex),
            Err(e) => Err(e.into()),
        }
    }

    pub fn byte(&mut self) -> Result<u8, Error> {
        let mut byte_slice: [u8; 1] = [0; 1];
        self.fill_completely(&mut byte_slice)?;
        Ok(byte_slice[0])
    }

    pub fn field(&mut self) -> Result<Field, Error> {
        const FIELD_SIZE: usize = size_of::<Field>();
        let mut field: [u8; FIELD_SIZE] = [0; FIELD_SIZE];
        self.fill_completely(&mut field)?;
        Ok(slice_to_field(&field, self.swap_endian))
    }

    pub fn string(&mut self, length: usize) -> Result<String, Error> {
        let mut buf = vec![0u8; length as usize];
        self.fill_completely(&mut buf)?;
        Ok(String::from_utf8(buf)?)
    }

    pub fn blocks(&mut self) -> Result<Vec<Block>, Error> {
        let num_blocks = self.field()? as usize;
        let mut result = Vec::with_capacity(num_blocks);
        for _ in 0..num_blocks {
            result.push(Block {
                start: self.field()?,
                length: 0, // will be assigned in the next loop
            });
        }

        for block in &mut result {
            block.length = self.field()?;
        }
        Ok(result)
    }

    pub fn skip_blocks(&mut self) -> Result<(), Error> {
        let num_blocks = self.field()? as usize;
        let skip = num_blocks * 2 * size_of::<Field>();
        self.reader.seek(SeekFrom::Current(skip as i64))?;
        Ok(())
    }

    /// Read bytes from the BufferedReader until the buffer is completely full
    fn fill_completely(&mut self, buf: &mut [u8]) -> Result<(), Error> {
        let n_bytes = buf.len();
        let mut bytes_read = 0;
        while bytes_read < n_bytes {
            // our buffered reader doesn't guarantee that it's always reading enough bytes at once
            bytes_read += self.reader.read(&mut buf[bytes_read..])?;
        }
        assert_eq!(bytes_read, n_bytes);
        Ok(())
    }
}

fn slice_to_field(slice: &[u8; 4], swap_endian: bool) -> Field {
    let mut result = 0;
    if swap_endian {
        for byte in slice.iter().rev() {
            result <<= 8;
            result += *byte as Field;
        }
    } else {
        for byte in slice {
            result <<= 8;
            result += *byte as Field;
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_slice_to_field() {
        let slice: &[u8; 4] = &[0x1A, 0x41, 0x27, 0x43];
        assert_eq!(slice_to_field(slice, false), 0x1A412743);
        let slice: &[u8; 4] = &[0x43, 0x27, 0x41, 0x1A];
        assert_eq!(slice_to_field(slice, false), 0x4327411A);

        assert_eq!(
            slice_to_field(&[0x1A, 0x41, 0x27, 0x43], true),
            slice_to_field(&[0x43, 0x27, 0x41, 0x1A,], false)
        );

        let slice: &[u8; 4] = &[0, 2, 3, 4];
        assert_eq!(slice_to_field(slice, false), 2 * 65536 + 3 * 256 + 4);
    }
}
