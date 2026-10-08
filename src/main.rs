/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

//! Command-line entry point.

use clap::Parser;

use hwd::cli::Args;
use hwd::common::run_test;
use hwd::prng::Prng;

fn main() -> std::io::Result<()> {
    let args = Args::parse();
    args.validate();

    eprintln!("Generator: {}", Prng::NAME);

    // Report the kernel's transparent-huge-page policy once: if it reads "[never]",
    // MADV_HUGEPAGE is ignored system-wide and the large buffers stay base-paged
    // regardless of what the allocation requests (a separate, system-level cause).
    #[cfg(target_os = "linux")]
    if let Ok(thp) = std::fs::read_to_string("/sys/kernel/mm/transparent_hugepage/enabled") {
        eprintln!("Transparent huge pages: {}", thp.trim());
    }

    // As in the C implementation, the exit status is zero also when the test
    // stops early because of a low p-value or an overflow.
    run_test(&args, args.parallel_cpus(), &mut std::io::stdout().lock())?;
    Ok(())
}
