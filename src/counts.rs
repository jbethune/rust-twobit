use crate::error::{Error, Result};

/// Number or percentage of bases of each type in a sequence.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BaseCounts<T> {
    pub a: T,
    pub c: T,
    pub g: T,
    pub t: T,
    pub n: T,
}

impl<T> BaseCounts<T> {
    #[allow(dead_code, clippy::many_single_char_names)]
    pub(crate) const fn new(a: T, c: T, g: T, t: T, n: T) -> Self {
        Self { a, c, g, t, n }
    }
}

impl BaseCounts<usize> {
    #[must_use]
    pub const fn sum(&self) -> usize {
        self.a + self.c + self.g + self.t + self.n
    }

    #[must_use]
    pub fn percentages(&self) -> BaseCounts<f64> {
        let mut result = BaseCounts::default();
        let sum = self.sum() as f64;
        result.a = (self.a as f64) / sum;
        result.c = (self.c as f64) / sum;
        result.g = (self.g as f64) / sum;
        result.t = (self.t as f64) / sum;
        result.n = (self.n as f64) / sum;
        result
    }

    #[inline]
    pub(crate) fn update(&mut self, nuc: u8) -> Result<()> {
        // no assumptions about softmasked sequences - "t" is the same as "T"
        match nuc {
            b'A' | b'a' => self.a += 1,
            b'C' | b'c' => self.c += 1,
            b'G' | b'g' => self.g += 1,
            b'T' | b't' => self.t += 1,
            b'N' => self.n += 1,
            _ => return Err(Error::BadNucleotide(char::from(nuc))),
        }
        Ok(())
    }
}

impl From<BaseCounts<usize>> for BaseCounts<f64> {
    fn from(counts: BaseCounts<usize>) -> Self {
        counts.percentages()
    }
}
