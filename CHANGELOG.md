# Change Log

## [0.1.1] - 2024-10-09

### Changed

- `--prng-bits` defaults to the number of bits of output of the generator,
  given by the new constant `Prng::BITS`: previously, testing a 32-bit
  generator without `--prng-bits 32` reported a spurious overflow, whereas now
  `-w 32` is sufficient, and `-w 64` is rejected. A value of `--prng-bits`
  larger than the output size of the generator causes a warning.

### Fixed

- In parallel runs, the additional small-counter arrays no longer exceed
  their memory bound for small _k_.

- `run_test` no longer panics when given `Some(0)` threads.

- The README provides options to replicate the results of the paper for all
  generators (in particular, `xoroshiro128+` was tested with the transitional
  variant).

## [0.1.0] - 2026-10-08

### New

- First release: a Rust port of the C implementation of the test for
  Hamming-weight dependencies described in "A New Test for Hamming-Weight
  Dependencies", by David Blackman and Sebastiano Vigna, with identical
  results.

- Parallel generation (`-P`): the output is identical to that of a
  sequential run.

- Arbitrary jumps (`try_skip`) for the 𝐅₂-linear generators of the paper.
