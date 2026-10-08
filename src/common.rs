/*
 * SPDX-FileCopyrightText: 2004-2016 David Blackman
 * SPDX-FileCopyrightText: 2017-2026 David Blackman and Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

//! The test driver, shared by sequential and parallel runs.
//!
//! The output of the generator is examined in batches (see
//! [`Mode::batch_size`]); at the end of each batch the small counters are
//! moved into the large counters, and possibly a *p*-value is computed and
//! reported.
//!
//! A parallel run processes *rounds* of consecutive batches: the iterations of
//! the main loop in a round are split into contiguous ranges, one per thread,
//! and each thread updates its own copy of the small counters for each batch
//! its range overlaps. At the end of the round, for each batch, in order, the
//! copies are merged into the large counters (see [`desat`]) and progress is
//! reported exactly as in a sequential run. A thread starts by examining (but
//! not counting) the words preceding its range, so that its signature is the
//! same of a sequential run. The output of a parallel run is thus identical to
//! that of a sequential run.

use std::io::{self, Write};
use std::marker::PhantomData;
use std::mem::size_of;
use std::time::Instant;

use bytemuck::Pod;
use mmap_rs::{MmapFlags, MmapMut, MmapOptions};
use rayon::prelude::*;

use crate::cli::Args;
use crate::prng::Prng;
use crate::scan::{CountSum, Mode, SigState, desat, scan_dispatch};
use crate::stats::{compute_pvalue, format_p_value};

/// A zero-initialized buffer of `T` allocated with `mmap()`.
pub struct Buffer<T> {
    mmap: MmapMut,
    len: usize,
    _marker: PhantomData<T>,
}

impl<T: Pod> Buffer<T> {
    /// Allocates a zero-initialized buffer of `len` elements, prefaulting it
    /// as transparent huge pages.
    ///
    /// # Implementation Details
    ///
    /// [`MmapFlags::POPULATE`] would prefault the buffer as base (4 KiB)
    /// pages, so we use [`MmapFlags::TRANSPARENT_HUGE_PAGES`] and prefault
    /// the buffer by touching one byte every 2 MiB. If transparent huge pages
    /// are disabled, the buffer is still prefaulted, as base pages.
    ///
    /// [`MmapFlags::POPULATE`]: mmap_rs::MmapFlags::POPULATE
    /// [`MmapFlags::TRANSPARENT_HUGE_PAGES`]: mmap_rs::MmapFlags::TRANSPARENT_HUGE_PAGES
    pub fn new(len: usize) -> Self {
        // A failed allocation is not a bug, so we exit with an explanation
        // rather than panicking.
        fn alloc_error(n: usize, detail: &dyn std::fmt::Display) -> ! {
            eprintln!(
                "\ncannot allocate a buffer of {n} elements: {detail}; \
                 reduce k or the number of threads"
            );
            std::process::exit(1);
        }
        let bytes_len = len
            .checked_mul(size_of::<T>())
            .unwrap_or_else(|| alloc_error(len, &"the size in bytes overflows usize"))
            .max(1);
        let mut mmap = MmapOptions::new(bytes_len)
            .and_then(|options| {
                options
                    .with_flags(MmapFlags::TRANSPARENT_HUGE_PAGES)
                    .map_mut()
            })
            .unwrap_or_else(|e| alloc_error(len, &e));
        const HUGE_PAGE: usize = 2 * 1024 * 1024;
        let bytes: &mut [u8] = &mut mmap;
        bytes
            .par_chunks_mut(HUGE_PAGE)
            .for_each(|chunk| chunk[0] = 0);
        Self {
            mmap,
            len,
            _marker: PhantomData,
        }
    }

    pub fn as_slice(&self) -> &[T] {
        &bytemuck::cast_slice(&self.mmap[..])[..self.len]
    }

    pub fn as_mut_slice(&mut self) -> &mut [T] {
        &mut bytemuck::cast_slice_mut(&mut self.mmap[..])[..self.len]
    }
}

/// How a test ended.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Outcome {
    /// All the requested data was examined.
    Completed,
    /// A reported *p*-value was below the threshold given with `--low-pv`.
    LowPValue,
    /// A small counter overflowed, and a *p*-value of 10⁻¹⁰⁰ was reported.
    Overflow,
}

/// The amounts of data, in bytes, after which progress is reported. When an
/// amount is used, it is multiplied by ten; the zero marks the end of the
/// list, and causes a restart from the beginning.
const PROGRESS_SIZES: [i64; 13] = [
    100000000, 125000000, 150000000, 175000000, 200000000, 250000000, 300000000, 400000000,
    500000000, 600000000, 700000000, 850000000, 0,
];

/// The target number of iterations of the main loop per thread in a round of
/// a parallel run.
const TARGET_ITERS_PER_THREAD: u64 = 1 << 24;

/// The minimum number of iterations of the main loop per thread in a round of
/// a parallel run. It must be at least the number of warm-up iterations.
const MIN_ITERS_PER_THREAD: u64 = 1 << 16;

/// The maximum number of small counters, beyond one array per thread, that a
/// parallel run allocates to make rounds longer.
const EXTRA_SMALL_COUNTERS: usize = 1 << 26;

/// The state of a test, except for the generator and the small counters.
struct Hwd<'a, W: Write> {
    mode: Mode,
    k: usize,
    numcats: usize,
    trans: bool,
    progress: bool,
    /// The number of bytes to examine, or a negative value.
    n: i64,
    low_pv: f64,
    max_batch_size: Option<u64>,
    /// The time at which the test started.
    tstart: Instant,
    count_sum: Buffer<CountSum>,
    norm: Buffer<f64>,
    /// The number of bytes examined so far.
    pos: i64,
    progress_sizes: [i64; 13],
    progress_index: usize,
    next_progress: i64,
    out: &'a mut W,
}

impl<W: Write> Hwd<'_, W> {
    /// Returns the number of words of the next batch if it starts at `pos`
    /// bytes, or zero if the test is over.
    fn next_batch_size(&self, pos: i64) -> u64 {
        if self.n >= 0 && pos >= self.n {
            return 0;
        }
        let mut batch = self.mode.batch_size(self.k) as i64;
        if let Some(max) = self.max_batch_size {
            batch = batch.min(max as i64);
        }
        let bytes = self.mode.word_bytes() as i64;
        if self.n >= 0 && (self.n - pos) / bytes < batch {
            batch = ((self.n - pos) / bytes) & !7;
        }
        batch as u64
    }

    /// Returns the sizes of the next (at most) `max_batches` batches.
    fn plan_round(&self, max_batches: usize) -> Vec<u64> {
        let mut batches = Vec::new();
        let mut pos = self.pos;
        while batches.len() < max_batches {
            let batch = self.next_batch_size(pos);
            if batch == 0 {
                break;
            }
            batches.push(batch);
            pos += (batch * self.mode.word_bytes()) as i64;
        }
        batches
    }

    /// Ends a batch of `batch` words whose small counters are the sum of
    /// `small`, and whose Hamming weights sum to `tot_sums` (only if *w* =
    /// 128).
    ///
    /// The small counters are moved into the large counters; then, progress is
    /// possibly reported. Returns the outcome of the test if it must stop.
    fn end_batch(
        &mut self,
        small: &mut [&mut [u32]],
        batch: u64,
        tot_sums: u64,
    ) -> io::Result<Option<Outcome>> {
        let half_w = (self.mode.word_bits() / 2) as i64;
        let (c, s) = desat(small, self.count_sum.as_mut_slice(), half_w);
        // If w = 128 the sums might overflow, too.
        let overflow = c != batch || (self.mode == Mode::W128 && s != tot_sums);
        if overflow {
            eprintln!(
                "{}",
                if self.mode == Mode::W128 {
                    "Counters or values overflowed. Seriously non-random."
                } else {
                    "Counters overflowed. Seriously non-random."
                }
            );
            writeln!(self.out, "p = {}", format_p_value(1e-100))?;
            return Ok(Some(Outcome::Overflow));
        }

        self.pos += (batch * self.mode.word_bytes()) as i64;

        if self.progress && self.pos >= self.next_progress {
            if self.analyze(false)? {
                return Ok(Some(Outcome::LowPValue));
            }
            self.progress_sizes[self.progress_index] *= 10;
            self.progress_index += 1;
            self.next_progress = self.progress_sizes[self.progress_index];
            if self.next_progress == 0 {
                self.progress_index = 0;
                self.next_progress = self.progress_sizes[0];
            }
        }
        Ok(None)
    }

    /// Computes and reports a *p*-value. Returns true if the *p*-value is
    /// below the threshold.
    fn analyze(&mut self, final_: bool) -> io::Result<bool> {
        let p = self.mode.central_prob();
        if (self.pos as f64) < 2.0 * (2.0 / (1.0 - p)).powf(self.k as f64) {
            writeln!(
                self.out,
                "WARNING: p-values are unreliable, you have to wait (insufficient data for meaningful answer)"
            )?;
        }

        let pvalue = compute_pvalue(
            self.count_sum.as_slice(),
            self.norm.as_mut_slice(),
            self.k,
            self.numcats,
            self.mode.word_bits(),
            self.trans,
            self.out,
        )?;
        let pos = self.pos;
        let elapsed = self.tstart.elapsed().as_secs_f64();
        writeln!(
            self.out,
            "processed {} bytes in {:.3} seconds ({:.4} GB/s, {:.4} TB/h)\n",
            pos,
            elapsed,
            pos as f64 * 1e-9 / elapsed,
            pos as f64 * (3600.0 * 1e-12) / elapsed,
        )?;

        if final_ {
            writeln!(self.out, "final")?;
        }
        writeln!(self.out, "p = {}", format_p_value(pvalue))?;

        if pvalue < self.low_pv {
            return Ok(true);
        }

        if !final_ {
            writeln!(self.out, "------\n")?;
        }
        Ok(false)
    }
}

/// A contiguous range of iterations of a round processed by a thread within a
/// batch.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Segment {
    thread: usize,
    batch: usize,
    len: u64,
}

/// Returns the segments of a round, in order.
///
/// Thread `i` processes iterations starting at `thread_starts[i]` (the first
/// start is zero) up to the next start, or the end of the round; batch `j`
/// contains `batch_lens[j]` iterations.
fn segments(thread_starts: &[u64], batch_lens: &[u64]) -> Vec<Segment> {
    let total: u64 = batch_lens.iter().sum();
    let thread_end = |i: usize| thread_starts.get(i + 1).copied().unwrap_or(total);
    let mut segs = Vec::new();
    let (mut thread, mut batch) = (0, 0);
    let mut batch_end = batch_lens[0];
    let mut pos = 0;
    while pos < total {
        let end = thread_end(thread).min(batch_end);
        segs.push(Segment {
            thread,
            batch,
            len: end - pos,
        });
        pos = end;
        if pos == thread_end(thread) {
            thread += 1;
        }
        if pos == batch_end && batch + 1 < batch_lens.len() {
            batch += 1;
            batch_end += batch_lens[batch];
        }
    }
    segs
}

/// Scans sequentially the generator starting from `start`, and returns its
/// states after each of the numbers of calls in `offsets`, which must be
/// sorted.
fn prescan(start: Prng, offsets: &[u64]) -> Vec<Prng> {
    let mut p = start;
    let mut done = 0;
    offsets
        .iter()
        .map(|&off| {
            for _ in done..off {
                p.next_u64();
            }
            done = off;
            p
        })
        .collect()
}

/// Returns the description of how the output of the generator is produced,
/// for the header.
fn generation_desc(num_cpus: Option<usize>, skip_capable: bool) -> String {
    match num_cpus {
        None => "sequentially".to_string(),
        Some(n) => format!(
            "using {} parallel generator{} ({})",
            n,
            if n == 1 { "" } else { "s" },
            if skip_capable {
                "jump-ahead"
            } else {
                "pre-scan"
            }
        ),
    }
}

/// Runs the test, writing the report to `out`, sequentially if `num_cpus` is
/// `None`, and otherwise using `num_cpus` parallel generators.
///
/// The report is identical to that of the original C implementation, and does
/// not depend on `num_cpus`, except for the timing information.
pub fn run_test(args: &Args, num_cpus: Option<usize>, out: &mut impl Write) -> io::Result<Outcome> {
    let mode = args.mode();
    let k = args.dim;
    let size = 3usize.pow(k as u32);
    let third = (size / 3) as u32;
    let trans = args.transitions;
    let skip_capable = Prng::new(args.seed).try_skip(0).is_ok();

    // The number of iterations of the main loop in a full batch.
    let mut batch = mode.batch_size(k);
    if let Some(max) = args.max_batch_size() {
        batch = batch.min(max);
    }
    // In a parallel run, a round contains enough batches to give each thread
    // enough work, within the memory limit.
    let (max_batches, num_arrays) = match num_cpus {
        None => (1, 1),
        Some(n) => {
            let wanted = (n as u64 * TARGET_ITERS_PER_THREAD).div_ceil(mode.iterations(batch));
            let max_batches = wanted.min((EXTRA_SMALL_COUNTERS / size) as u64).max(1) as usize;
            (max_batches, max_batches + n - 1)
        }
    };

    eprintln!("Seed: {:#018x}", args.seed);
    let gib = (size * (size_of::<CountSum>() + size_of::<f64>() + num_arrays * size_of::<u32>()))
        as f64
        / 2.0f64.powi(30);
    eprintln!(
        "Running a test for Hamming-weight dependencies with k = {} and {} categories {} on {}, analyzing {} (batches of {} words, {:.3} GiB RAM)",
        k,
        args.numcats(),
        generation_desc(num_cpus, skip_capable),
        mode.description(),
        if trans { "transitions" } else { "bits" },
        batch,
        gib
    );
    let bytes = args.bytes();
    eprintln!(
        "Examining {}{}",
        if bytes >= 0 {
            format!("{bytes} bytes")
        } else {
            "the output indefinitely".to_string()
        },
        args.low_pv
            .map_or(String::new(), |p| format!(", stopping if p < {p:e}"))
    );

    let mut hwd = Hwd {
        mode,
        k,
        numcats: args.numcats(),
        trans,
        progress: args.progress(),
        n: bytes,
        low_pv: args.low_pv(),
        max_batch_size: args.max_batch_size(),
        tstart: Instant::now(),
        count_sum: Buffer::new(size),
        norm: Buffer::new(size),
        pos: 0,
        progress_sizes: PROGRESS_SIZES,
        progress_index: 0,
        next_progress: PROGRESS_SIZES[0],
        out,
    };
    let mut small: Vec<Buffer<u32>> = (0..num_arrays).map(|_| Buffer::new(size)).collect();

    // As in the C implementation, timing starts after allocation.
    hwd.tstart = Instant::now();
    let mut prng = Prng::new(args.seed);
    let mut st = SigState::initial(size as u32);

    let outcome = match num_cpus {
        None => loop {
            let batch = hwd.next_batch_size(hwd.pos);
            if batch == 0 {
                break None;
            }
            let cs = small[0].as_mut_slice();
            let tot_sums = scan_dispatch(
                mode,
                trans,
                true,
                &mut prng,
                cs,
                &mut st,
                mode.iterations(batch),
                third,
            );
            if let Some(outcome) = hwd.end_batch(&mut [cs], batch, tot_sums)? {
                break Some(outcome);
            }
        },
        Some(num_cpus) => {
            // Iterations of warm-up necessary to fill the signature: k words,
            // plus one, as when testing transitions the Hamming weight of a
            // word depends on the last bit of the previous word.
            let warm_up = (k as u64 + 1).div_ceil(mode.words_per_iter());
            debug_assert!(warm_up <= MIN_ITERS_PER_THREAD);
            loop {
                let batches = hwd.plan_round(max_batches);
                if batches.is_empty() {
                    break None;
                }
                let batch_lens: Box<[u64]> = batches.iter().map(|&b| mode.iterations(b)).collect();
                let total: u64 = batch_lens.iter().sum();
                let num_threads = num_cpus.min((total / MIN_ITERS_PER_THREAD).max(1) as usize);
                let (base, rem) = (total / num_threads as u64, total % num_threads as u64);
                let starts: Box<[u64]> = (0..num_threads as u64)
                    .map(|i| i * base + i.min(rem))
                    .collect();
                let segs = segments(&starts, &batch_lens);

                // Threads but the first start warm_up iterations before their
                // range; the offsets are in calls to the generator.
                let offsets: Box<[u64]> = starts[1..]
                    .iter()
                    .map(|&s| (s - warm_up) * mode.calls_per_iter())
                    .collect();
                let snapshots = (!skip_capable).then(|| prescan(prng, &offsets));

                let mut per_thread: Vec<Vec<(&mut [u32], u64)>> =
                    (0..num_threads).map(|_| Vec::new()).collect();
                for (seg, buf) in segs.iter().zip(small.iter_mut()) {
                    per_thread[seg.thread].push((buf.as_mut_slice(), seg.len));
                }

                let (round_start, round_st) = (prng, st);
                let results: Vec<(Vec<u64>, Prng, SigState)> = per_thread
                    .into_par_iter()
                    .enumerate()
                    .map(|(i, thread_segs)| {
                        let (mut p, mut s) = if i == 0 {
                            (round_start, round_st)
                        } else {
                            let mut p = match &snapshots {
                                Some(snapshots) => snapshots[i - 1],
                                None => {
                                    let mut p = round_start;
                                    p.try_skip(offsets[i - 1]).expect(
                                        "try_skip must succeed for every offset or for none",
                                    );
                                    p
                                }
                            };
                            // The initial signature is irrelevant, as the
                            // warm-up replaces all its trits.
                            let mut s = SigState::initial(size as u32);
                            scan_dispatch(
                                mode,
                                trans,
                                false,
                                &mut p,
                                &mut [],
                                &mut s,
                                warm_up,
                                third,
                            );
                            (p, s)
                        };
                        let sums = thread_segs
                            .into_iter()
                            .map(|(cs, len)| {
                                scan_dispatch(mode, trans, true, &mut p, cs, &mut s, len, third)
                            })
                            .collect();
                        (sums, p, s)
                    })
                    .collect();

                // The last thread ends where the next round starts.
                let last = results.last().expect("at least one thread");
                (prng, st) = (last.1, last.2);
                let sums: Vec<u64> = results.into_iter().flat_map(|r| r.0).collect();

                // End the batches in order, as a sequential run would.
                let mut first = 0;
                let mut stop = None;
                for (j, &batch) in batches.iter().enumerate() {
                    let end = first + segs[first..].iter().take_while(|s| s.batch == j).count();
                    let mut arrays: Vec<&mut [u32]> = small[first..end]
                        .iter_mut()
                        .map(|b| b.as_mut_slice())
                        .collect();
                    let tot_sums = sums[first..end].iter().sum();
                    if let Some(outcome) = hwd.end_batch(&mut arrays, batch, tot_sums)? {
                        stop = Some(outcome);
                        break;
                    }
                    first = end;
                }
                if stop.is_some() {
                    break stop;
                }
            }
        }
    };

    let outcome = match outcome {
        Some(outcome) => outcome,
        None => {
            if hwd.analyze(true)? {
                Outcome::LowPValue
            } else {
                Outcome::Completed
            }
        }
    };
    eprintln!(
        "Test completed in {:.2} seconds",
        hwd.tstart.elapsed().as_secs_f64()
    );
    Ok(outcome)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_segments() {
        // Two threads, three batches.
        let segs = segments(&[0, 15], &[10, 10, 10]);
        let expected = [(0, 0, 10), (0, 1, 5), (1, 1, 5), (1, 2, 10)];
        assert_eq!(
            segs.iter()
                .map(|s| (s.thread, s.batch, s.len))
                .collect::<Vec<_>>(),
            expected
        );
        // Thread and batch boundaries coinciding.
        let segs = segments(&[0, 10, 20], &[10, 10, 10]);
        let expected = [(0, 0, 10), (1, 1, 10), (2, 2, 10)];
        assert_eq!(
            segs.iter()
                .map(|s| (s.thread, s.batch, s.len))
                .collect::<Vec<_>>(),
            expected
        );
        // One batch, many threads.
        let segs = segments(&[0, 3, 6], &[10]);
        let expected = [(0, 0, 3), (1, 0, 3), (2, 0, 4)];
        assert_eq!(
            segs.iter()
                .map(|s| (s.thread, s.batch, s.len))
                .collect::<Vec<_>>(),
            expected
        );
    }
}
