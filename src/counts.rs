//! Count information on sequences
//!
//! This module makes no assumptions about softmasked sequences. An "t" is the same as a "T".

use std::convert::From;

/// Number of bases of each type in a sequence
#[derive(PartialEq, Debug, Default)]
pub struct BaseCounts {
    pub a: usize,
    pub c: usize,
    pub g: usize,
    pub t: usize,
    pub n: usize,
}

impl BaseCounts {
    pub fn sum(&self) -> usize {
        self.a + self.c + self.g + self.t + self.n
    }
}

/// Percentages of bases of each type in a sequence
#[derive(PartialEq, Debug, Default)]
pub struct BasePercentages {
    pub a: f64,
    pub c: f64,
    pub g: f64,
    pub t: f64,
    pub n: f64,
}

impl From<BaseCounts> for BasePercentages {
    fn from(counts: BaseCounts) -> Self {
        let mut result = BasePercentages::default();
        let sum = counts.sum() as f64;
        result.a = (counts.a as f64) / sum;
        result.c = (counts.c as f64) / sum;
        result.g = (counts.g as f64) / sum;
        result.t = (counts.t as f64) / sum;
        result.n = (counts.n as f64) / sum;
        result
    }
}
