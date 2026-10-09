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
particular for generators based on 𝐅₂-linear transformations such as the [dSFMT],
[`xoroshiro128+`], and [WELL512].

The crate is a port of the [original C implementation]: for the same generator and
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
each time; the test stops as soon as the final _p_-value of a report (the line
`p = …`) is below the threshold set with `--low-pv`. In the paper, tests were
performed with _w_ = 32 or _w_ = 64 and _k_ between 8 and 19, stopping after a
petabyte of data or at a _p_-value below 10⁻²⁰.

The test needs 28 · 3*ᵏ* bytes of memory, which become about 1.2 GB for _k_ = 16
and 33 GB for _k_ = 19. To improve locality, words are counted in packed 32-bit
counters, which are moved into 64-bit counters at the end of each batch. Batch
sizes are chosen so that the probability that a packed counter overflows is
below 10⁻¹⁰⁰: if this happens, the test reports a _p_-value of 10⁻¹⁰⁰.

# Results

The following table, from the [page of the test], reports for some
𝐅₂-linear generators the amount of output, in bytes, after which the test
yields a _p_-value below 10⁻²⁰, and the faulty signature. Ranges (→) appear
when several variants were tested (e.g., several parameters of the Mersenne
Twister): a missing right extreme means that some instances did not fail the
test within a petabyte.

| PRNG                            | _w_ | Period   | _p_ = 10⁻²⁰ @      | Faulty signature                             |
| ------------------------------- | --: | -------- | ------------------ | -------------------------------------------- |
| `xorshift128+`                  |  64 | 2¹²⁸ − 1 | 6 × 10⁹            | `00000012` (transitional)                    |
| `xoroshiro128+`                 |  64 | 2¹²⁸ − 1 | 5 × 10¹²           | `00000012` (transitional)                    |
| Tiny Mersenne Twister (64 bits) |  32 | 2¹²⁷ − 1 | 8 × 10¹³ →         | `00000202`                                   |
| Tiny Mersenne Twister (32 bits) |  32 | 2¹²⁷ − 1 | 4 × 10¹³ →         | `10001021`                                   |
| SFMT (607 bits)                 |  64 | 2⁶⁰⁷ − 1 | 4 × 10⁸            | `001000001000`                               |
| dSFMT (521 bits)                |  32 | 2⁵²¹ − 1 | 6 × 10¹²           | `1001000100100010`                           |
| Mersenne Twister (521 bits)     |  32 | 2⁵²¹ − 1 | 4 × 10¹⁰ →         | `1000000100000000`, `2000000100000000`       |
| Mersenne Twister (607 bits)     |  32 | 2⁶⁰⁷ − 1 | 4 × 10⁸ → 4 × 10¹⁰ | `1000000001000000000`, `2000000001000000000` |
| WELL512a (512 bits)             |  32 | 2⁵¹² − 1 | 3 × 10¹⁵           | `2001002200000000`                           |

Since the generators of this crate are seeded differently (see below), the
amount of output and the faulty signatures might be slightly different.

# Parallel generation

With option `-P`, the iterations of the test are split into contiguous ranges,
one for each thread: jump-capable generators jump to the start of their range,
whereas the others reach it with a sequential pre-scan. Each thread counts the
words of its range in its own copy of the packed counters, starting a few words
before its range so as to have the same signature as a sequential run; at the
end of each batch, the copies are added together, and since packed counters are
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
  upper 32 bits, and `--prng-bits` defaults to the number of bits of output of
  the generator. With 64-bit outputs, `-w 32` splits each output into two words,
  upper half first.

- Generators are initialized from a 64-bit seed (option `-S`), and generators
  with more than 64 bits of state fill it using [SplitMix64]: the faulty
  signatures and the amount of data needed to find bias are thus slightly
  different from those reported in the paper.

- _p_-values are printed with full precision, and the reports contain the exact
  number of bytes examined, but not the current date.

# Usage

The generator to test is selected at compilation time using Cargo features:
besides `incr`, a counter useful only to test the crate, the available
generators are those of the paper, `xorshift128`, `xorshift128plus`,
`xorshift1024`, `xorshift1024plus`, `xoroshiro128`, `xoroshiro128plus`,
`xoroshiro1024`, `xoroshiro1024plus`, and the 32-bit generator `well512a`. For
example,

```text
cargo run -r -F xoroshiro128plus -- -P -t --progress --low-pv=1e-20 1e15
Generator: xoroshiro128+ (24, 16, 37)
Seed: 0x0000000000000000
Running a test for Hamming-weight dependencies with k = 8 and 5 categories using 10 parallel generators (jump-ahead) on 64-bit words, the full outputs, analyzing transitions (batches of 2311072 words, 0.002 GiB RAM)
Examining 1000000000000000 bytes, stopping if p < 1e-20
mix3 extreme = 2.04196 (sig = 00200000) weight 1 (16), p-value = 0.4895265563711558
mix3 extreme = 2.81178 (sig = 10000200) weight 2 (112), p-value = 0.4248712876027823
mix3 extreme = 2.90032 (sig = 10020100) weight 3 (448), p-value = 0.8123541528946184
mix3 extreme = 4.32909 (sig = 02020202) weight 4 (1120), p-value = 0.016630046075989773
mix3 extreme = 4.32152 (sig = 12021102) weight >=5 (4864), p-value = 0.07260107083813899
bits per word = 64 (analyzing transitions); min category p-value = 0.016630046075989773

processed 110931456 bytes in 0.052 seconds (2.1269 GB/s, 7.6568 TB/h)

p = 0.08043025669891243
------

[...]

mix3 extreme = 1.72707 (sig = 01000000) weight 1 (16), p-value = 0.7550051714620413
mix3 extreme = 10.43265 (sig = 00000012) weight 2 (112), p-value = 1.9702516536983596e-23
mix3 extreme = 3.10710 (sig = 10000110) weight 3 (448), p-value = 0.5713964898923892
mix3 extreme = 3.83667 (sig = 11021000) weight 4 (1120), p-value = 0.1303729781298133
mix3 extreme = 3.63088 (sig = 10100222) weight >=5 (4864), p-value = 0.7469225372955031
bits per word = 64 (analyzing transitions); min category p-value = 1.9702516536983596e-23

processed 5000013516288 bytes in 199.500 seconds (25.0627 GB/s, 90.2258 TB/h)

p = 9.851258268491799e-23
Test completed in 199.50 seconds
```

will run the transitional variant of the test on `xoroshiro128+` with _w_ = 64
and _k_ = 8, as in the paper, on an Apple M1 Max with 10 cores, finding bias
after 5 · 10¹² bytes in less than three and a half minutes, with the same amount
of data and the same faulty signature reported in the paper (`00000012`). The
test processes about 25 GB/s; the C implementation, which is sequential,
processes about 3.3 GB/s on the same hardware, and would need about 25 minutes.

The default parameters (_w_ = 64, _k_ = 8) are those used for 64-bit generators
with 128 bits of state; generators with a larger state need a larger _k_. To
replicate the tests of the paper, use

```text
cargo run -r -F FEATURE -- -P OPTIONS --progress --low-pv=1e-20 1e15
```

with the following features and options (the `+` generators were tested with
the transitional variant, and the batch size of the `xorshift1024` generators
was reduced to obtain more frequent reports); the amount of data and the faulty
signature are those reported in the paper:

| Feature             | Options                         | _p_ = 10⁻²⁰ @ | Faulty signature   |
| ------------------- | ------------------------------- | ------------- | ------------------ |
| `xorshift128`       |                                 | 8 × 10⁸       | `00000021`         |
| `xorshift128plus`   | `-t`                            | 6 × 10⁹       | `00000012`         |
| `xorshift1024`      | `-k 16 --max-batch-size 1e7`    | 6 × 10⁸       | `2000000000000001` |
| `xorshift1024plus`  | `-t -k 16 --max-batch-size 1e8` | 9 × 10⁹       | `2000000000000001` |
| `xoroshiro128`      |                                 | 1 × 10¹⁰      | `00000012`         |
| `xoroshiro128plus`  | `-t`                            | 5 × 10¹²      | `00000012`         |
| `xoroshiro1024`     | `-k 16`                         | 5 × 10¹²      | `1100000000000001` |
| `xoroshiro1024plus` | `-t -k 16`                      | 4 × 10¹³      | `1100000000000001` |
| `well512a`          | `-w 32 -k 16`                   | 3 × 10¹⁵      | `2001002200000000` |

WELL512a needs more than a petabyte, so for it you must raise the limit, or omit
it to run indefinitely.

The repository configures Cargo to compile for the native CPU (`-C
target-cpu=native`), as otherwise, for example, population counts on x86-64
would not use the `POPCNT` instruction. Since `cargo install` does not use this
configuration, install the test with, for example,

```text
RUSTFLAGS="-C target-cpu=native" cargo install hwd -F xoroshiro128plus
```

# Adding your own generator

To add a new generator:

1. Add to `Cargo.toml` a feature enabling the internal `_prng` marker (e.g.,
   `mygen = ["_prng"]`): otherwise, the build stops with a “no PRNG selected”
   error.

2. Add to the [`prng`] module, conditionally on the feature, a `Prng` type with
   the constants `NAME` and `BITS` and the methods `new`, `next_u64`, and
   `try_skip` (see the [`prng`] module documentation). The `try_skip` method is
   mandatory, but it can fail for every offset: in this case, parallel runs
   reach the start of their range with a sequential pre-scan.

3. For 𝐅₂-linear generators, implement the [`LinearGenerator`] trait, define
   the characteristic polynomial `CHARPOLY` of the transition map, and
   implement `try_skip` with the `f2_try_skip!` macro, as done, for example, for
   `xoroshiro128`; the feature must be added to the `cfg` lists of the macro and
   of the `charpoly_tests` module. To obtain the characteristic polynomial,
   define `CHARPOLY` as zero and run `cargo test -F mygen test_charpoly`: the
   test will fail, printing the characteristic polynomial, computed as the
   minimal polynomial of the sequence of the lowest bit of the state (the test
   checks that it has full degree).

The [`prng`] module has the same interface as that of the [`coll-birth`] crate,
except for the `BITS` constant.

[A New Test for Hamming-Weight Dependencies]: https://doi.org/10.1145/3527582
[original C implementation]: https://prng.di.unimi.it/hwd.php
[page of the test]: https://prng.di.unimi.it/hwd.php
[SplitMix64]: https://prng.di.unimi.it/splitmix64.c
[`coll-birth`]: https://crates.io/crates/coll-birth
[`prng`]: https://docs.rs/hwd/latest/hwd/prng/index.html
[`f2`]: https://docs.rs/hwd/latest/hwd/f2/index.html
[`LinearGenerator`]: https://docs.rs/hwd/latest/hwd/f2/trait.LinearGenerator.html
[dSFMT]: https://doi.org/10.1007/978-3-642-04107-5_38
[`xoroshiro128+`]: https://doi.org/10.1145/3460772
[WELL512]: https://doi.org/10.1145/1132973.1132974