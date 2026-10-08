/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

//! Command-line argument definitions.

use clap::Parser;

use crate::mode::Mode;

#[derive(Parser, Debug, Clone)]
#[command(
    author,
    version,
    about = "Runs a test for Hamming-weight dependencies; use RAYON_NUM_THREADS to customize the number of threads.",
    next_line_help = true,
    max_term_width = 100
)]
pub struct Args {
    /// Number of bytes to examine, possibly in floating-point notation (e.g., 1e15); if omitted or
    /// negative, the test runs indefinitely, reporting progress.​
    #[arg(allow_negative_numbers = true, value_parser = parse_f64)]
    pub n: Option<f64>,

    /// Size w of the examined words: 32, 64, or 128.​
    #[arg(short = 'w', long, default_value_t = 64)]
    pub word_bits: usize,

    /// Number of output bits of the generator: 32 (in the upper bits of the 64-bit output) or 64;
    /// with 64-bit outputs and w = 32 each output provides two words, upper bits first.​
    #[arg(long, default_value_t = 64)]
    pub prng_bits: usize,

    /// Length k of signatures, between 1 and 19; memory usage is proportional to 3ᵏ.​
    #[arg(short = 'k', long, default_value_t = 8)]
    pub dim: usize,

    /// Number of categories of signatures, between 1 and k (default: ⌊k / 2⌋ + 1).​
    #[arg(short = 'c', long)]
    pub categories: Option<usize>,

    /// Test transitions (the words xored with themselves shifted by one bit) instead of bits.​
    #[arg(short = 't', long)]
    pub transitions: bool,

    /// Report p-values periodically, and not only at the end.​
    #[arg(long)]
    pub progress: bool,

    /// Stop as soon as a reported p-value is below this threshold (default: the smallest positive
    /// normal double).​
    #[arg(long, value_parser = parse_f64)]
    pub low_pv: Option<f64>,

    /// Upper bound on the number of words in a batch, possibly in floating-point notation; smaller
    /// batches make progress reports more frequent (MAX_BATCH_SIZE in the C implementation).​
    #[arg(long, value_parser = parse_f64)]
    pub max_batch_size: Option<f64>,

    /// PRNG seed; accepts decimal, or a 0x/0o/0b prefix for hexadecimal, octal, or binary (underscores may separate digits).​
    #[arg(short = 'S', long, default_value_t = 0, value_parser = parse_u64)]
    pub seed: u64,

    /// Generate data in parallel: the orbit is split into contiguous segments, one generated per
    /// thread; jump-capable generators jump to each segment start, others reach it with a
    /// sequential pre-scan; the output is identical to that of a sequential run.​
    #[arg(short = 'P', long)]
    pub parallel: bool,
}

impl Args {
    /// Returns the word mode determined by the word size and the number of
    /// output bits of the generator.
    pub fn mode(&self) -> Mode {
        Mode::new(self.word_bits, self.prng_bits).expect("arguments must be validated")
    }

    /// Returns the number of categories.
    pub fn numcats(&self) -> usize {
        self.categories.unwrap_or(self.dim / 2 + 1)
    }

    /// Returns the number of bytes to examine, or a negative value to run
    /// indefinitely.
    ///
    /// As in the C implementation, the argument is truncated to an integer.
    pub fn bytes(&self) -> i64 {
        self.n.map_or(-1, |n| n as i64)
    }

    /// Returns whether progress must be reported: either it was requested, or
    /// the number of bytes is not positive (as in the C implementation).
    pub fn progress(&self) -> bool {
        self.progress || self.bytes() <= 0
    }

    /// Returns the threshold for early termination.
    pub fn low_pv(&self) -> f64 {
        self.low_pv.unwrap_or(f64::MIN_POSITIVE)
    }

    /// Returns the upper bound on the number of words in a batch, if any.
    pub fn max_batch_size(&self) -> Option<u64> {
        self.max_batch_size.map(|b| b as u64)
    }

    /// Validates argument combinations, reporting inconsistencies through
    /// [`Args::die`] in clap's own error style.
    pub fn validate(&self) {
        if ![32, 64, 128].contains(&self.word_bits) {
            Self::die(&format!(
                "the word size ({}) must be 32, 64, or 128",
                self.word_bits
            ));
        }
        if ![32, 64].contains(&self.prng_bits) {
            Self::die(&format!(
                "the number of output bits ({}) must be 32 or 64",
                self.prng_bits
            ));
        }
        if Mode::new(self.word_bits, self.prng_bits).is_none() {
            Self::die(&format!(
                "{}-bit words are not supported with {}-bit outputs",
                self.word_bits, self.prng_bits
            ));
        }
        if !(1..=19).contains(&self.dim) {
            Self::die(&format!("k ({}) must be between 1 and 19", self.dim));
        }
        if !(1..=self.dim).contains(&self.numcats()) {
            Self::die(&format!(
                "the number of categories ({}) must be between 1 and k ({})",
                self.numcats(),
                self.dim
            ));
        }
        if let Some(b) = self.max_batch_size {
            // An odd batch would split an iteration of the main loop when
            // w = 32 and each output provides two words.
            if !(b >= 2.0 && b < 2f64.powi(63) && b.fract() == 0.0 && b % 2.0 == 0.0) {
                Self::die(&format!(
                    "the maximum batch size ({b}) must be a positive even integer"
                ));
            }
        }
    }

    /// Reports an argument error in clap's own style (message plus usage footer)
    /// and exits with status 2, indistinguishable from clap's native diagnostics.
    pub fn die(msg: &str) -> ! {
        use clap::CommandFactory;
        Args::command()
            .error(clap::error::ErrorKind::ValueValidation, msg)
            .exit()
    }

    /// Resolved parallel generation width: `None` for sequential generation,
    /// `Some(n)` with `n` the Rayon pool size (governed by `RAYON_NUM_THREADS`)
    /// when `--parallel` is set.
    pub fn parallel_cpus(&self) -> Option<usize> {
        if self.parallel {
            Some(crate::util::parallelism())
        } else {
            None
        }
    }
}

/// Parses a floating-point number, accepting also the notation of C's
/// `strtod()` for exponents (e.g., `1E15`).
fn parse_f64(value: &str) -> Result<f64, String> {
    value
        .trim()
        .parse::<f64>()
        .map_err(|e| format!("invalid number {value:?}: {e}"))
}

/// Parses an unsigned 64-bit integer in decimal, or in hexadecimal, octal, or
/// binary when prefixed with `0x`, `0o`, or `0b` (case-insensitive).
/// Underscores are allowed as digit separators (e.g. `0xDEAD_BEEF`).
fn parse_u64(value: &str) -> Result<u64, String> {
    let trimmed = value.trim();
    let (radix, digits) = match trimmed.get(..2) {
        Some("0x") | Some("0X") => (16, &trimmed[2..]),
        Some("0o") | Some("0O") => (8, &trimmed[2..]),
        Some("0b") | Some("0B") => (2, &trimmed[2..]),
        _ => (10, trimmed),
    };
    let digits: String = digits.chars().filter(|&c| c != '_').collect();
    if digits.is_empty() {
        return Err(format!("invalid integer: {value:?}"));
    }
    u64::from_str_radix(&digits, radix).map_err(|e| format!("invalid integer {value:?}: {e}"))
}
