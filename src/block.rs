use std::ops::{Bound, Deref, Range, RangeBounds};

/// Block mask for sequence regions
///
/// Blocks are used to mark regions as soft-masked or hard-masked.
pub type Block = Range<usize>;

/// Sorted collection of block masks
#[derive(Debug, Clone)]
pub struct Blocks(pub Vec<Block>);

impl Deref for Blocks {
    type Target = [Block];

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl Blocks {
    #[inline]
    pub fn overlaps(&self, range: impl RangeBounds<usize>) -> impl Iterator<Item = &Block> {
        let clone_bound = |bound: Bound<&usize>| {
            // stable in Rust 1.55 but not in 1.51
            match bound {
                Bound::Unbounded => Bound::Unbounded,
                Bound::Included(x) => Bound::Included(*x),
                Bound::Excluded(x) => Bound::Excluded(*x),
            }
        };
        let start = clone_bound(range.start_bound());
        let end = clone_bound(range.end_bound());
        self.iter()
            .skip_while(move |block| match start {
                Bound::Included(v) => block.end <= v,
                Bound::Excluded(v) => block.end <= v + 1,
                Bound::Unbounded => false,
            })
            .take_while(move |block| match end {
                Bound::Included(v) => block.start < v + 1,
                Bound::Excluded(v) => block.start < v,
                Bound::Unbounded => true,
            })
    }

    #[must_use]
    pub fn count(&self) -> usize {
        self.iter().fold(0, |acc, block| acc + block.len())
    }

    pub fn apply_masks<const HARD: bool>(&self, seq: &mut [u8], offset: usize) {
        let (start, end) = (offset, offset + seq.len());
        for block in self.overlaps(start..end) {
            let seq_start = start.max(block.start) - start;
            let seq_end = end.min(block.end) - start;
            seq[seq_start..seq_end]
                .iter_mut()
                .for_each(|nuc| *nuc = if HARD { b'N' } else { nuc.to_ascii_lowercase() });
        }
    }
}
