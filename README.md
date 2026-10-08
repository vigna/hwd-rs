# A test for Hamming-weight dependencies in pseudorandom number generators

[![crates.io](https://img.shields.io/crates/v/hwd.svg)](https://crates.io/crates/hwd)
[![docs.rs](https://docs.rs/hwd/badge.svg)](https://docs.rs/hwd)
[![rustc](https://img.shields.io/badge/rustc-1.85+-red.svg)](https://rust-lang.github.io/rfcs/2495-min-rust-version.html)
[![CI](https://github.com/vigna/hwd-rs/actions/workflows/rust.yml/badge.svg)](https://github.com/vigna/hwd-rs/actions)
![license](https://img.shields.io/crates/l/hwd)
[![downloads](https://img.shields.io/crates/d/hwd)](https://crates.io/crates/hwd)
[![coveralls](https://coveralls.io/repos/github/vigna/hwd-rs/badge.svg?branch=main)](https://coveralls.io/github/vigna/hwd-rs?branch=main)

This crate implements the test for Hamming-weight dependencies described in
“[A New Test for Hamming-Weight Dependencies]”, by David Blackman and Sebastiano
Vigna. The test finds bias induced by dependencies among the Hamming weights
(the number of ones) of the outputs of a pseudorandom number generator (PRNG), in
particular for generators based on 𝐅₂-linear transformations such as the dSFMT,
`xoroshiro128+`, and WELL512.

The crate is a port of the original C implementation: for the same generator and
parameters, the result of the test is the same. Moreover, the crate can use
multiple cores to generate and examine the output of the generator in parallel,
again without changing the result of the test.

# The test

The test examines a sequence of _w_-bit words extracted from the output of the
generator: usually _w_ is the size of the output, but it is possible to split
each output into two halves, or to glue consecutive outputs into larger words.
Each word is classified by its Hamming weight as _average_ or _extremal_: the
trit (base-3 digit) of a word is 1 if its weight is one of the 2ℓ + 1 central
weights, whose overall probability is close to 1/2 (ℓ = 1 for _w_ = 32, ℓ = 2
for _w_ = 64, ℓ = 3 for _w_ = 128), and otherwise 0 or 2, depending on whether
the weight is smaller or larger.

The trits of _k_ consecutive words form a _signature_, and for each of the 3*ᵏ*
signatures the test computes the average Hamming weight of the words that follow
the signature in the sequence. For a random sequence, these averages, suitably
normalized, look like independent standard normal values. The vector of
normalized averages is then multiplied by the _k_-th Kronecker power of a 3 × 3
orthonormal matrix, in the same vein as the Walsh–Hadamard transform: since the
transform is unitary, for a random sequence the transformed values still look
like independent standard normal values. However, each transformed value
combines the averages of many signatures: a zero trit in its index means “don't
care about that previous word”, a one compares words with low and high weight,
and a two compares extremal and average words. Values with few nonzero trits
thus make simple dependencies emerge, such as “the weight of a word depends on
the weight of the word three steps before”.

The transformed values are divided into _C_ categories by the number of nonzero
trits of their index (all values with at least _C_ nonzero trits are in the last
category; by default, _C_ = ⌊_k_ / 2⌋ + 1). The test computes the _p_-value of
the most extreme value of each category, corrected for the size of the
category, and the final _p_-value is the smallest category _p_-value, corrected
for the number of categories.

For each category, the test reports also the signature of the most extreme
value, written with the most recent word on the right: when a generator fails,
the faulty signature provides insight into its structure. For example,
`xorshift1024` fails with signature `2000000000000001`, as each output depends
on the first and last word of its 16-word state.

There is also a _transitional_ variant of the test (option `-t`), which xors
the sequence of words, seen as a stream of bits, with itself shifted by one bit,
and then examines the result: in this case, the test looks for dependencies
between bit transitions.

# Running the test

The output of the generator is examined in batches, whose size depends on _w_
and _k_ (it can be reduced with `--max-batch-size`). The test stops after a
given number of bytes, or runs indefinitely if no number is specified. With
`--progress` (or when running indefinitely), _p_-values are reported after
about 10⁸, 1.25 · 10⁸, 1.5 · 10⁸, 1.75 · 10⁸, 2 · 10⁸, 2.5 · 10⁸, 3 · 10⁸, 4 ·
10⁸, 5 · 10⁸, 6 · 10⁸, 7 · 10⁸, 8.5 · 10⁸ bytes, and so on, multiplying by ten
each time; the test stops as soon as a _p_-value is below the threshold set with
`--low-pv`. In the paper, tests were performed with _w_ = 32 or _w_ = 64 and _k_
between 8 and 19, stopping after a petabyte of data or at a _p_-value below
10⁻²⁰.

The test needs 28 · 3*ᵏ* bytes of memory, which become about 1.2 GB for _k_ = 16
and 33 GB for _k_ = 19. To improve locality, words are counted in packed 32-bit
counters, which are moved into 64-bit counters at the end of each batch. Batch
sizes are chosen so that the probability that a packed counter overflows is
below 10⁻¹⁰⁰: if this happens, the test reports a _p_-value of 10⁻¹⁰⁰.

# Parallel generation

With option `-P`, the iterations of the test are split into contiguous ranges,
one for each thread: jump-capable generators jump to the start of their range,
whereas the others reach it with a sequential pre-scan. Each thread counts the
words of its range in its own copy of the packed counters, starting a few words
before its range so to have the same signature as a sequential run; at the end
of each batch, the copies are added together, and since packed counters are
updated by addition modulo 2³², the result is exactly the one of a sequential
run, even in case of overflow. Thus, the output of a parallel run is identical
to that of a sequential run.

Each thread needs its own copy of the packed counters, so a parallel run with
_t_ threads uses at least 4 · (_t_ − 1) · 3*ᵏ* additional bytes of memory
(small values of _k_ use a few more copies to make work units longer). The
number of threads is governed by the environment variable `RAYON_NUM_THREADS`
(by default, the number of cores).

All the 𝐅₂-linear generators of the paper provide arbitrary jumps using their
characteristic polynomial (see the [`f2`] module).

# Differences from the C implementation

- Parameters are run-time options rather than compile-time macros: `-w`
  (`HWD_BITS`), `--prng-bits` (`HWD_PRNG_BITS`), `-k` (`HWD_DIM`), `-c`
  (`HWD_NUMCATS`), and `--max-batch-size` (`MAX_BATCH_SIZE`). The generator is
  selected at compile time using Cargo features.

- Generators return 64-bit values: 32-bit generators return their output in the
  upper 32 bits, and must be tested with `--prng-bits 32`. With 64-bit outputs,
  `-w 32` splits each output into two words, upper half first.

- Generators are initialized from a 64-bit seed (option `-S`), and generators
  with more than 64 bits of state fill it using [SplitMix64]: the faulty
  signatures and the amount of data needed to find bias are thus slightly
  different from those reported in the paper.

- _p_-values are printed with full precision, and the reports contain the exact
  number of bytes examined, but not the current date.

# Usage

The generator to test is selected at compilation time using Cargo features.
For example,

```text
cargo run -r -F xoroshiro128plus -- -P --progress --low-pv=1e-20 1e15
Generator: xoroshiro128+ (24, 16, 37)
Seed: 0x0000000000000000
Running a test for Hamming-weight dependencies with k = 8 and 5 categories using 10 parallel generators (jump-ahead) on 64-bit words, the full outputs, analyzing bits (batches of 2311072 words, 0.002 GiB RAM)
Examining 1000000000000000 bytes, stopping if p < 1e-20
mix3 extreme = 2.57485 (sig = 00020000) weight 1 (16), p-value = 0.14893176321045537
mix3 extreme = 2.95929 (sig = 00001020) weight 2 (112), p-value = 0.2924060708331487
mix3 extreme = 2.97046 (sig = 01002200) weight 3 (448), p-value = 0.736610976840799
mix3 extreme = 2.96745 (sig = 21200100) weight 4 (1120), p-value = 0.9655482751525518
mix3 extreme = 4.32354 (sig = 22120101) weight >=5 (4864), p-value = 0.0719628674515154
bits per word = 64 (analyzing bits); min category p-value = 0.0719628674515154

processed 110931456 bytes in 0.057 seconds (1.9524 GB/s, 7.0286 TB/h)

p = 0.3116223400948752
------

[...]

mix3 extreme = 2.21810 (sig = 00002000) weight 1 (16), p-value = 0.3498212859734305
mix3 extreme = 11.28875 (sig = 00000012) weight 2 (112), p-value = 1.6703369397814162e-27
mix3 extreme = 3.43838 (sig = 20000210) weight 3 (448), p-value = 0.23067832880149697
mix3 extreme = 3.21723 (sig = 20221000) weight 4 (1120), p-value = 0.7655706409615056
mix3 extreme = 3.37608 (sig = 21101022) weight >=5 (4864), p-value = 0.9720588981323016
bits per word = 64 (analyzing bits); min category p-value = 1.6703369397814162e-27

processed 10000008544000 bytes in 435.037 seconds (22.9866 GB/s, 82.7517 TB/h)

p = 8.351684698907081e-27
Test completed in 435.04 seconds
```

will test `xoroshiro128+` with _w_ = 64 and _k_ = 8 on an Apple M1 Max with 10
cores, finding bias after 10¹³ bytes in about seven minutes, with the same
faulty signature reported in the paper (`00000012`). The test processes about
23 GB/s; the C implementation, which is sequential, processes about 3.3 GB/s on
the same hardware, and would need about 50 minutes.

The default parameters (_w_ = 64, _k_ = 8) are those used for 64-bit generators
with 128 bits of state; generators with a larger state need a larger _k_. For
example, to replicate the tests of the paper on `xorshift1024` and on the
32-bit generator WELL512a, use

```text
cargo run -r -F xorshift1024 -- -P -k 16 --progress --low-pv=1e-20 1e15
cargo run -r -F well512a -- -P -w 32 --prng-bits 32 -k 16 --progress --low-pv=1e-20 1e15
```

# Adding your own generator

To add a new generator, add a feature in `Cargo.toml` and a corresponding
implementation in the [`prng`] module. If skipping is possible, you can
implement the `try_skip` method, which must succeed for every offset or for
none (see the [`prng`] module documentation): for 𝐅₂-linear generators, you
just need to implement the [`LinearGenerator`] trait and provide the
characteristic polynomial of the transition map, as done, for example, for
`xoroshiro128`.

The [`prng`] module has the same interface as that of the [`coll-birth`] crate.

[A New Test for Hamming-Weight Dependencies]: https://doi.org/10.1145/3527582
[TestU01]: https://doi.org/10.1145/1268776.1268777
[gjrand]: https://gjrand.sourceforge.net/
[SplitMix64]: https://prng.di.unimi.it/splitmix64.c
[`coll-birth`]: https://crates.io/crates/coll-birth
[`prng`]: https://docs.rs/hwd/latest/hwd/prng/index.html
[`f2`]: https://docs.rs/hwd/latest/hwd/f2/index.html
[`LinearGenerator`]: https://docs.rs/hwd/latest/hwd/f2/trait.LinearGenerator.html
