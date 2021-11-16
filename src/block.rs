//! Masking reqions of a sequence

use std::cmp::{max, min};
use std::ops::Range;

use crate::types::Field;

/// A block mask for sequence regions
///
/// Blocks are used to mark regions as soft-masked or hard-masked.
#[derive(PartialEq, Debug, Clone, Copy)]
pub struct Block {
    pub start: Field,
    pub length: Field,
}

impl Block {
    /// Create a new block
    #[must_use]
    pub const fn new(start: Field, length: Field) -> Self {
        Self { start, length }
    }

    /// Determine an overlap with another block
    ///
    /// Convenience function to determine overlaps between sequence regions.
    ///
    /// Note: Blocks of the same type should usually not overlap in practise.
    #[must_use]
    pub fn overlap(&self, other: &Self) -> Option<Range<usize>> {
        let start = max(self.start, other.start);
        let end = min(self.start + self.length, other.start + other.length);
        if start < end {
            Some(start as usize..end as usize)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_overlap() {
        let a = Block::new(5, 10);
        let b = Block::new(10, 10);
        let c = Block::new(15, 10);
        let d = Block::new(100, 25);
        let e = Block::new(90, 40);
        assert_eq!(a.overlap(&b), Some(10..15));
        assert_eq!(a.overlap(&c), None);
        assert_eq!(b.overlap(&c), Some(15..20));
        assert_eq!(a.overlap(&a), Some(5..15));
        assert_eq!(a.overlap(&d), None);
        assert_eq!(d.overlap(&a), None);
        assert_eq!(d.overlap(&e), Some(100..125));
        assert_eq!(e.overlap(&d), Some(100..125));
        assert_eq!(
            Block::new(67920552, 300000).overlap(&Block::new(67859142, 61415)),
            Some(67920552..(67859142 + 61415))
        );
    }
}
