/*
 * SPDX-FileCopyrightText: 2004-2016 David Blackman
 * SPDX-FileCopyrightText: 2017-2026 David Blackman and Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

//! Small and large counters.
//!
//! To improve locality, the main loop updates *small counters*: a single
//! `u32` packs a 13-bit count (upper bits) and a 19-bit sum of Hamming weights
//! (lower bits), so both fields are updated with a single addition (see
//! [`SUM_BITS`]). At the end of each batch, the small counters are moved into
//! the *large counters* and zeroed, checking that no count overflowed.
//!
//! Since updates are additions modulo 2³², the small counters of a batch can
//! be split among threads, each updating its own copy for a contiguous part of
//! the batch: the sum modulo 2³² of the copies is exactly the value a single
//! sequential scan would have produced, including overflows. This is what
//! makes parallel runs faithful.

use std::marker::PhantomData;
use std::mem::size_of;

use bytemuck::{Pod, Zeroable};
use mmap_rs::{MmapFlags, MmapMut, MmapOptions};
use rayon::prelude::*;

use crate::scan::SUM_BITS;

/// A large counter: the number of words following a signature, and the sum of
/// the differences between their Hamming weights and *w* / 2.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Pod, Zeroable)]
pub struct CountSum {
    pub c: u64,
    pub s: i64,
}

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

/// Moves the content of the small counters into the large counters, and
/// zeroes the small counters (`desat()` in the C implementation).
///
/// The small counters of a batch are the sum modulo 2³² of the arrays in
/// `small`, all of the same length of `large`. `half_w` is *w* / 2.
///
/// Returns the sum of the counts and the sum of the sums of the small
/// counters, which the caller compares with the number of words in the batch
/// and with their total Hamming weight to detect overflows.
pub fn desat(small: &mut [&mut [u32]], large: &mut [CountSum], half_w: i64) -> (u64, u64) {
    const CHUNK: usize = 1 << 16;

    fn desat_chunk(small: &mut [&mut [u32]], large: &mut [CountSum], half_w: i64) -> (u64, u64) {
        let (mut c, mut s) = (0, 0);
        for (i, l) in large.iter_mut().enumerate() {
            let mut st = 0u32;
            for a in small.iter_mut() {
                st = st.wrapping_add(a[i]);
                a[i] = 0;
            }
            let count = st >> SUM_BITS;
            let sum = st & ((1 << SUM_BITS) - 1);
            c += count as u64;
            s += sum as u64;
            l.c += count as u64;
            // In the small counters the Hamming weights are stored as they
            // are; in the large counters, as differences from w / 2.
            l.s += sum as i64 - half_w * count as i64;
        }
        (c, s)
    }

    if large.len() <= CHUNK {
        return desat_chunk(small, large, half_w);
    }

    // Transpose the arrays into lists of chunks, one list per chunk of large
    // counters, so that chunks can be processed in parallel.
    let num_chunks = large.len().div_ceil(CHUNK);
    let mut chunks: Vec<Vec<&mut [u32]>> = (0..num_chunks)
        .map(|_| Vec::with_capacity(small.len()))
        .collect();
    for a in small.iter_mut() {
        for (list, chunk) in chunks.iter_mut().zip(a.chunks_mut(CHUNK)) {
            list.push(chunk);
        }
    }
    chunks
        .into_par_iter()
        .zip(large.par_chunks_mut(CHUNK))
        .map(|(mut list, l)| desat_chunk(&mut list, l, half_w))
        .reduce(|| (0, 0), |(c0, s0), (c1, s1)| (c0 + c1, s0 + s1))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Updates a small counter as the main loop does.
    fn update(cs: &mut [u32], i: usize, bc: u32) {
        cs[i] = cs[i].wrapping_add(bc + (1 << SUM_BITS));
    }

    // Splitting the updates of a batch among several arrays and merging them
    // must give the same result as a single array, also when counters
    // overflow.
    #[test]
    fn test_split_desat_matches_single() {
        for size in [100usize, (1 << 16) + 1000] {
            let mut x = 0x9e3779b97f4a7c15u64;
            let mut updates = Vec::new();
            for j in 0..40_000 {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                // Signature 0 is hit very often, to cause overflows.
                let i = if j % 3 == 0 {
                    0
                } else {
                    (x % size as u64) as usize
                };
                updates.push((i, (x >> 32) as u32 % 65));
            }

            let mut single = vec![0u32; size];
            for &(i, bc) in &updates {
                update(&mut single, i, bc);
            }
            let mut large_single = vec![CountSum::default(); size];
            let res_single = desat(&mut [&mut single[..]], &mut large_single, 32);

            let mut parts = vec![vec![0u32; size]; 3];
            for (j, &(i, bc)) in updates.iter().enumerate() {
                update(&mut parts[j * 3 / updates.len()], i, bc);
            }
            let mut large_split = vec![CountSum::default(); size];
            let mut views: Vec<&mut [u32]> = parts.iter_mut().map(|p| &mut p[..]).collect();
            let res_split = desat(&mut views, &mut large_split, 32);

            assert_eq!(res_split, res_single);
            assert_eq!(large_split, large_single);
            // The overflow must have been detected.
            assert_ne!(res_single.0, updates.len() as u64);
            assert!(parts.iter().flatten().all(|&v| v == 0));
            assert!(single.iter().all(|&v| v == 0));
        }
    }
}
