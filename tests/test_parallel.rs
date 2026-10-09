/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

//! Tests checking that parallel runs give the same results as sequential runs.

use hwd::cli::Args;
use hwd::common::{Outcome, run_test};

fn make_args(n: f64, word_bits: usize, prng_bits: usize, dim: usize) -> Args {
    Args {
        n: Some(n),
        word_bits,
        prng_bits,
        dim,
        categories: None,
        transitions: false,
        progress: true,
        low_pv: None,
        max_batch_size: None,
        seed: 0x0123_4567_89AB_CDEF,
        parallel: false,
    }
}

/// Runs a test and returns its outcome and its report, without timing
/// information.
fn run(args: &Args, num_cpus: Option<usize>) -> anyhow::Result<(Outcome, String)> {
    let mut out = Vec::new();
    let outcome = run_test(args, num_cpus, &mut out)?;
    let report = String::from_utf8(out)?
        .lines()
        .map(|l| {
            if l.starts_with("processed ") {
                l.split(" in ").next().unwrap_or(l)
            } else {
                l
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    Ok((outcome, report))
}

/// Checks that parallel runs with several numbers of threads give the same
/// outcome and report as a sequential run.
fn check(args: &Args) -> anyhow::Result<()> {
    let seq = run(args, None)?;
    // The report must contain at least a p-value.
    assert!(seq.1.contains("p = "), "{}", seq.1);
    for num_cpus in [1, 2, 3, 7] {
        let par = run(args, Some(num_cpus))?;
        assert_eq!(
            par, seq,
            "{num_cpus} parallel generators disagree with a sequential run for {args:?}"
        );
    }
    Ok(())
}

#[test]
fn test_w64() -> anyhow::Result<()> {
    let mut args = make_args(1.3e8, 64, 64, 8);
    check(&args)?;
    args.transitions = true;
    check(&args)
}

#[test]
fn test_w32_p64() -> anyhow::Result<()> {
    // Odd k, so warm-up is not a whole number of iterations.
    let mut args = make_args(1.3e8, 32, 64, 7);
    check(&args)?;
    args.transitions = true;
    check(&args)
}

#[test]
fn test_w32_p32() -> anyhow::Result<()> {
    let mut args = make_args(1.3e8, 32, 32, 6);
    args.categories = Some(2);
    check(&args)?;
    args.transitions = true;
    check(&args)
}

#[test]
fn test_w128() -> anyhow::Result<()> {
    let mut args = make_args(1.3e8, 128, 64, 5);
    check(&args)?;
    args.transitions = true;
    check(&args)
}

// Small batches, so that rounds contain many batches, and progress is reported
// in the middle of rounds.
#[test]
fn test_small_batches() -> anyhow::Result<()> {
    let mut args = make_args(1.2e8, 64, 64, 4);
    args.max_batch_size = Some(1000.0);
    args.transitions = true;
    check(&args)?;
    args.word_bits = 32;
    args.max_batch_size = Some(1002.0);
    check(&args)
}

// A run stopping at the first report (or before, if counters overflow, as it
// happens with 32-bit generators).
#[test]
fn test_low_pv() -> anyhow::Result<()> {
    let mut args = make_args(1.3e8, 64, 64, 8);
    args.low_pv = Some(1.0);
    check(&args)?;
    let (outcome, _) = run(&args, Some(3))?;
    assert_ne!(outcome, Outcome::Completed);
    Ok(())
}

// A run shorter than a batch, which uses a single thread.
#[test]
fn test_short() -> anyhow::Result<()> {
    let mut args = make_args(1e5, 64, 64, 8);
    args.progress = false;
    check(&args)
}

// Zero parallel generators are treated as one.
#[test]
fn test_zero_threads() -> anyhow::Result<()> {
    let args = make_args(1e7, 64, 64, 4);
    assert_eq!(run(&args, Some(0))?, run(&args, None)?);
    Ok(())
}
