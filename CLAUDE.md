# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

A Rust port of the C implementation of the test for Hamming-weight dependencies of "A New Test for Hamming-Weight Dependencies" (Blackman and Vigna, ACM TOMACS 2022). The C code lives in a separate repository (`../hwd/c/hwd.c`, with the paper `../hwd/hwd.tex`). **The port must give the same results as the C code**: same batches, same counters, same floating-point operations in the same order, same p-values. Only formatting may differ (p-values are printed with full precision, byte counts are exact, no date). It shares its structure and the interface of its `prng` module with `../coll-birth-rs`, as the two may eventually be merged into a general PRNG tester: keep them similar.

## Build Commands

The PRNG under test is selected at build time via Cargo features. Exactly one feature must be enabled:

```bash
cargo build --release --features <prng_feature>
cargo run --release --features <prng_feature> -- [N] [options]
cargo test --features <prng_feature>
cargo clippy --all-targets --features <prng_feature> -- -D warnings
```

Available PRNG features are the 𝐅₂-linear generators of the paper: `xorshift128`, `xorshift128plus`, `xorshift1024`, `xorshift1024plus`, `xoroshiro128`, `xoroshiro128plus`, `xoroshiro1024`, `xoroshiro1024plus`, `well512a`. CI runs clippy and tests for every feature.

## Checking Equivalence with the C Implementation

`c-check/check.sh` compiles `hwd.c` from `$HWD_C_DIR` (default `../hwd/c`) with `c-check/prngs_hwd.c`, which defines the same generators as `src/prng.rs` (same SplitMix64 seeding), builds the crate with the corresponding feature, runs both, and compares outputs after `c-check/normalize.py` drops timing and rounds Rust's full-precision p-values to C's `%.3g`. `c-check/check-all.sh` runs a configuration matrix (all generators, all word modes, `-t`, categories, small batches, `--low-pv`, sequential and `-P`). Overflows are obtained by examining the 32-bit `well512a` as 64-bit or 128-bit words, without `-t` (C has the output in the lower bits, the crate in the upper bits). Run it after any change to `scan.rs`, `stats.rs`, `common.rs`, or a generator; when adding a generator to the matrix, add it to `prngs_hwd.c` too. Work files go into `target/c-check`.

## Architecture

The crate is both a library (`src/lib.rs`, crate name `hwd`) and a thin binary (`src/main.rs`), as in coll-birth, so that integration tests can call `run_test` directly.

- **`scan.rs`**: everything from the output of the generator to the large counters.
  - `Mode` (`W32P32`, `W32P64`, `W64`, `W128`): the legal combinations of word size *w* and generator output size (`HWD_BITS`/`HWD_PRNG_BITS` in C), with words and calls per main-loop iteration, the central-trit probability *P*, and the batch-size tables copied verbatim from C (indexed by *k*).
  - The main loop `scan::<W, P64, TRANS, COUNT>` (C's `scan_batch`), specialized by `scan_dispatch`. `SigState` holds the signature and the transition carry. Counter updates use `get_unchecked_mut`, made sound by assertions at entry (`cs.len() == 3 · third`, initial signature in range): removing the bounds check is worth ~18% at *w* = 64. For *w* = 128, the two 64-bit words are combined into a `u128` (first output in the low half) so a single population count is used, which also expresses transitions as `w ^ (w << 1) ^ t`.
  - `CountSum` (large counters) and `desat`, which merges any number of copies of the small counters by wrapping addition into the large counters and zeroes them. Wrapping addition makes the merge equal to a single sequential array, even when counters overflow: this is what makes parallel runs faithful.
- **`stats.rs`**: normalization, `mix3` (the transform *T*ₖ; note that, as in C, the third column of the base matrix has the opposite sign with respect to the paper), category extremes (first index wins ties, also in the parallel reduction), and p-values. `erfc` comes from the C math library via `extern "C"`; `ln_1p`/`exp_m1`/`powf`/`sqrt` also map to libm, so results are bit-identical to C on the same platform. Parallelism is used only where it cannot change results (element-wise operations, the recursion of `mix3`, extremes with ordered merge).
- **`common.rs`**: `Buffer` (mmap with transparent huge pages, prefaulted, as `alloc_mmap` in coll-birth's `common.rs`) and the driver `run_test(args, num_cpus, out)`, where `num_cpus` is `None` for a sequential run. `Hwd` holds the state shared by both paths (large counters, position, progress schedule) and implements `next_batch_size`, `end_batch` (desat, overflow check, progress) and `analyze` (C's `analyze`). The parallel path processes *rounds* of batches (enough to give each thread `TARGET_ITERS_PER_THREAD` iterations, within `EXTRA_SMALL_COUNTERS` of memory), splits the round's iterations into contiguous per-thread ranges, cuts them at batch boundaries into `segments`, and gives each segment its own small-counter array. Thread 0 continues from the previous round's final state; the other threads start `warm_up` = ⌈(*k* + 1) / words-per-iteration⌉ iterations before their range (*k* + 1 words because with `-t` the weight of a word depends on the last bit of the previous word), reached by `try_skip` or, for generators without jumps (none at present), a sequential `prescan`. After the parallel phase, batches are ended in order exactly as in a sequential run.
- **`f2.rs`**: arbitrary jumps for 𝐅₂-linear generators, a port of `f2x.c` from prng.di.unimi.it: `x_pow_mod` by square-and-multiply with bit-serial `mul_mod`, and `jump`, the accumulate-and-step loop over the `LinearGenerator` trait (`step`, `to_vector`, `set_vector`, with canonical rotation for generators with a rotating index). Characteristic polynomials are hardcoded in `prng.rs` (`CHARPOLY`); `minimal_polynomial` (Berlekamp–Massey) is test-only and checks them in `prng::charpoly_tests`.
- **`prng.rs`**: the generators of the paper, with the same interface as in coll-birth (one `pub struct Prng` per feature with `NAME`, `new(seed)`, `next_u64`, `try_skip`; 32-bit outputs in the upper bits). Their state is filled with SplitMix64 from the seed; the linear engines return their first state word, as in the C replication code, and the `+` variants the sum of two words.
- **`cli.rs`**: `Args` (clap) and validation. It accepts the C program's invocation (`-t --progress --low-pv=1e-20 1E15`).

Integration tests (`tests/test_parallel.rs`) compare the full reports of parallel runs with 1, 2, 3, and 7 threads against a sequential run, for all word modes, with and without transitions; they run with whatever generator feature is enabled. Since all generators jump, `prescan` is not exercised.

## CLI

```
hwd [N] [-w W] [--prng-bits B] [-k K] [-c C] [-t] [--progress] [--low-pv P] [--max-batch-size M] [-S SEED] [-P]
```

- `N`: bytes to examine (floating-point notation accepted, truncated as in C); if omitted or not positive, run indefinitely with progress.
- `-w`: word size (32, 64, 128; default 64); `--prng-bits`: 32 or 64 (default 64); `-k`: signature length (1–19, default 8); `-c`: categories (default ⌊k/2⌋+1).
- `-t`: transitional variant; `--progress`: report at 10⁸, 1.25·10⁸, … bytes (C's `progsize`); `--low-pv`: stop when a reported p-value is below the threshold (default `f64::MIN_POSITIVE`, as C's `DBL_MIN`); `--max-batch-size`: C's `MAX_BATCH_SIZE` (must be even).
- `-S`: seed; `-P`: parallel generation, using the Rayon pool (`RAYON_NUM_THREADS`).

The report goes to stdout, headers and diagnostics to stderr. The exit status is zero also on early termination, as in C.

## Key Details

- `.cargo/config.toml` sets `target-cpu=native`; release builds use fat LTO, one codegen unit, `panic = "abort"`.
- Memory is 28 · 3ᵏ bytes, plus 4 · 3ᵏ bytes per additional small-counter array in parallel runs.
- Rust never contracts floating-point operations into FMAs; the C code must be compiled without `-march=native`/`-ffp-contract=fast` for comparisons to be exact on x86.
