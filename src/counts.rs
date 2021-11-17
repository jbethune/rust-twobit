use crate::error::{Error, Result};

/// Number or percentage of bases of each type in a sequence.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BaseCounts {
    pub a: usize,
    pub c: usize,
    pub g: usize,
    pub t: usize,
    pub n: usize,
}

impl BaseCounts {
    #[must_use]
    pub const fn sum(&self) -> usize {
        self.a + self.c + self.g + self.t + self.n
    }

    #[must_use]
    pub fn percentages(&self) -> BasePercentages {
        let mut result = BasePercentages::default();
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

/// Percentages of bases of each type in a sequence
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BasePercentages {
    pub a: f64,
    pub c: f64,
    pub g: f64,
    pub t: f64,
    pub n: f64,
}

impl From<BaseCounts> for BasePercentages {
    fn from(counts: BaseCounts) -> Self {
        counts.percentages()
    }
}
