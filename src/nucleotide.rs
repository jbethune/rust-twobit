//! Type-safe representations of nucleotide sequences

use std::convert::{TryFrom, TryInto};
use std::fmt;
use std::u32;

use crate::error::Error;

/// All nucleotides that are supported by the 2bit format (case insensitive)
#[derive(Debug, Copy, Clone)]
pub enum Nucleotide {
    /// Thymine (and Uracil for RNA)
    T,
    /// Cytosine
    C,
    /// Adenine
    A,
    /// Guanine
    G,
    /// Any nucleotide/unknown nucleotide
    N,
}

impl Nucleotide {
    /// Encode `Nucleotide`s as bits according to the 2bit specification.
    ///
    /// The nucleotide `N` is encoded as `0`.
    #[must_use]
    pub const fn bits(&self) -> u8 {
        match self {
            Self::T | Self::N => 0,
            Self::C => 1,
            Self::A => 2,
            Self::G => 3,
        }
    }
}

impl TryFrom<u8> for Nucleotide {
    type Error = Error;
    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            84 | 85 | 116 | 117 => Ok(Self::T), // Thymine and Uracil
            67 | 99 => Ok(Self::C),
            65 | 97 => Ok(Self::A),
            71 | 103 => Ok(Self::G),
            78 => Ok(Self::N),
            _ => Err(Error::FileFormat(format!(
                "Bad ascii value for nucleotide: {}",
                value
            ))),
        }
    }
}

impl TryFrom<char> for Nucleotide {
    type Error = Error;
    fn try_from(value: char) -> Result<Self, Self::Error> {
        let numeric: u32 = value.into();
        let numeric: u8 = numeric
            .try_into()
            .map_err(|_| Self::Error::BadNucleotide(value))?;
        numeric.try_into()
    }
}

impl fmt::Display for Nucleotide {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let c = match self {
            Self::T => 'T',
            Self::C => 'C',
            Self::A => 'A',
            Self::G => 'G',
            Self::N => 'N',
        };
        write!(f, "{}", c)
    }
}
