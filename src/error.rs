use std::fmt;
use std::io;
use std::string::FromUtf8Error;

/// An error that may occur while reading and parsing a 2bit file.
#[derive(Debug)]
pub enum Error {
    /// IO problems
    IO(io::Error),
    /// Something is wrong with the file
    FileFormat(String),
    /// 2bit file version is too advanced
    UnsupportedVersion(String),
    /// An entry (usually a chromosome) is not in the 2bit file
    MissingName(String),
    /// An invalid nucleotide is in the file
    BadNucleotide(char),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::IO(io) => io.fmt(f),
            Error::FileFormat(msg) => write!(f, "{}", msg),
            Error::UnsupportedVersion(version) => write!(f, "Unsupported version: {}", version),
            Error::MissingName(name) => write!(f, "Missing name: {}", name),
            Error::BadNucleotide(nuc) => write!(f, "Bad nucleotide: {}", nuc),
        }
    }
}

impl std::error::Error for Error {}

impl From<io::Error> for Error {
    fn from(io: io::Error) -> Self {
        Self::IO(io)
    }
}

impl From<FromUtf8Error> for Error {
    fn from(utf_error: FromUtf8Error) -> Self {
        Self::FileFormat(format!("{}", utf_error))
    }
}

/// A type alias for `Result<T, twobit::Error>`.
pub type Result<T> = std::result::Result<T, Error>;
