# Changelog

## 0.2.0

Release date: Nov 25, 2021.

### Added

- Added support for generic readers: not just buffered file readers, but
  any reader object that implements both `Read` and `Seek`. As a consequence,
  `TwoBitFile` is now generic over `R: Read + Seek`.
- Added `boxed()` method that converts the underlying reader objects to a boxed
  trait object in order to help with type erasure if needed.
- Added constructor that reads a file into memory: `open_and_read()`.
- Added type aliases for convenience: `TwoBitPhysicalFile` when reading from a
  buffered file object; `TwoBitMemoryFile` when reading a file into memory first;
  `BoxTwoBitFile` when boxing the internal reader.
- Added constructor that reads from a given buffer: `from_buf()`.
- Added `error::Result` type alias.

### Changed

- Core parsing routines have been rewritten from scratch, yielding an order
  of magnitude performance improvement.
- Remove `softmask_enabled` from constructors. This option is now disabled by
  default and can be enabled via `enable_softmask()` method.
- All methods that involve reading are now `&mut self`.
- `sequence()` method has been renamed to `read_sequence()`.
- Replaced `chroms()` method with `chrom_names()` and `chrom_sizes()` methods.
- Removed all potential panics in reading methods.
- Updated the docs, the code is now clippy-compliant, enabled CI, added benches.
- Removed `full_*()` methods; instead, all reading methods now expect ranges.
  In order to read the entire sequence, `..` range can be passed in.

## 0.1.1

Release date: Sep 28, 2020.

### Fixed

- Fixed off-by-1 error in block overlaps.

## 0.1.0

Release date: Aug 4, 2020.

Initial release.
