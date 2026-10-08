# Change Log

## [0.1.0] - 2026-10-08

### New

- First release: a Rust port of the C implementation of the test for
  Hamming-weight dependencies described in "A New Test for Hamming-Weight
  Dependencies", by David Blackman and Sebastiano Vigna, with identical
  results.

- Faithful parallel generation (`-P`): the output is identical to that of a
  sequential run.

- Arbitrary jumps (`try_skip`) for the 𝐅₂-linear generators of the paper.
