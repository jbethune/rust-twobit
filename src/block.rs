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
    pub fn iter_overlaps(&self, range: impl RangeBounds<usize>) -> impl Iterator<Item = &Block> {
        let start = range.start_bound().cloned();
        let end = range.end_bound().cloned();
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

    pub fn apply_masks<const HARD: bool>(&self, seq: &mut [u8], range: Block) {
        let (start, end) = (range.start, range.end);
        for block in self.iter_overlaps(range) {
            let seq_start = start.max(block.start) - start;
            let seq_end = end.min(block.end) - start;
            for i in seq_start..seq_end {
                unsafe {
                    *seq.get_unchecked_mut(i) = if HARD {
                        b'N'
                    } else {
                        seq.get_unchecked(i).to_ascii_lowercase()
                    }
                }
            }
        }
    }
}
