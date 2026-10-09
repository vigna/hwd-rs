/*
 * SPDX-FileCopyrightText: 2004-2016 David Blackman
 * SPDX-FileCopyrightText: 2017-2026 David Blackman and Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

//! Scanning the output of the generator: the words examined by the test, their
//! signatures, and the counters.
//!
//! The output of the generator is turned into a sequence of *w*-bit words as
//! specified by a [`Mode`]. The main loop, [`scan`], keeps track of the
//! signature of the last *k* words, and accounts for each word in the counter
//! of the signature of the words preceding it.
//!
//! To improve locality, the main loop updates *small counters*: a single
//! `u32` packs a 13-bit count (upper bits) and a 19-bit sum of Hamming weights
//! (lower bits), so both fields are updated with a single addition (see
//! [`SUM_BITS`]). At the end of each batch, [`desat`] moves the small counters
//! into the *large counters* ([`CountSum`]) and zeroes them.
//!
//! Since updates are additions modulo 2³², the small counters of a batch can
//! be split among threads, each updating its own copy for a contiguous part of
//! the batch: the sum modulo 2³² of the copies is exactly the value a single
//! sequential scan would have produced, including overflows.

use bytemuck::{Pod, Zeroable};
use rayon::prelude::*;

use crate::prng::Prng;

/// How the output of the generator is turned into the sequence of *w*-bit
/// words examined by the test.
///
/// The variants correspond to the legal combinations of `HWD_BITS` and
/// `HWD_PRNG_BITS` of the original C implementation. Generators with 32-bit
/// output return it in the upper 32 bits of
/// [`next_u64`](crate::prng::Prng::next_u64).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    /// 32-bit words, one for each output (its upper 32 bits).
    W32P32,
    /// 32-bit words, two for each output (first the upper 32 bits, then the
    /// lower 32 bits).
    W32P64,
    /// 64-bit words, one for each output.
    W64,
    /// 128-bit words, one for each pair of consecutive outputs.
    W128,
}

/// Batch sizes for *w* = 32, indexed by *k* (see [`Mode::batch_size`]).
const BATCH_SIZE_32: [u64; 20] = [
    0,
    16904,
    37848,
    88680,
    213360,
    520784,
    1280664,
    3160976,
    7815952,
    19342248,
    47885112,
    118569000,
    293614056,
    727107408,
    1800643824,
    4459239480,
    11043223056,
    27348419104,
    67728213816,
    167728896072,
];

/// Batch sizes for *w* = 64, indexed by *k* (see [`Mode::batch_size`]).
const BATCH_SIZE_64: [u64; 20] = [
    0, 14744, 28320, 56616, 116264, 242784, 512040, 1086096, 2311072, 4926224, 10510376, 22435504,
    47903280, 102294608, 218459240, 466556056, 996427288, 2128099936, 4545075936, 9707156552,
];

/// Batch sizes for *w* = 128, indexed by *k* (see [`Mode::batch_size`]).
const BATCH_SIZE_128: [u64; 20] = [
    0,
    14856,
    28792,
    58088,
    120392,
    253680,
    539816,
    1155104,
    2479360,
    5330680,
    11471256,
    24696808,
    53183328,
    114541856,
    246706584,
    531387952,
    1144590984,
    2465432776,
    5310537968,
    11438933136,
];

impl Mode {
    /// Returns the mode for the given word size *w* and generator output size,
    /// or `None` if the combination is not supported.
    pub const fn new(word_bits: usize, prng_bits: usize) -> Option<Self> {
        match (word_bits, prng_bits) {
            (32, 32) => Some(Self::W32P32),
            (32, 64) => Some(Self::W32P64),
            (64, 64) => Some(Self::W64),
            (128, 64) => Some(Self::W128),
            _ => None,
        }
    }

    /// Returns the size *w* of the examined words.
    pub const fn word_bits(self) -> usize {
        match self {
            Self::W32P32 | Self::W32P64 => 32,
            Self::W64 => 64,
            Self::W128 => 128,
        }
    }

    /// Returns the number of bytes of a word.
    pub const fn word_bytes(self) -> u64 {
        self.word_bits() as u64 / 8
    }

    /// Returns the number of words examined at each iteration of the main
    /// loop.
    pub const fn words_per_iter(self) -> u64 {
        match self {
            Self::W32P64 => 2,
            _ => 1,
        }
    }

    /// Returns the number of calls to the generator at each iteration of the
    /// main loop.
    pub const fn calls_per_iter(self) -> u64 {
        match self {
            Self::W128 => 2,
            _ => 1,
        }
    }

    /// Returns the number of iterations of the main loop necessary to examine
    /// `words` words (`TEST_ITERATIONS` in the C implementation).
    pub const fn iterations(self, words: u64) -> u64 {
        words / self.words_per_iter()
    }

    /// Returns the probability that a word has a central Hamming weight, that
    /// is, that its trit is one.
    // The literals of the C implementation.
    #[allow(clippy::excessive_precision)]
    pub const fn central_prob(self) -> f64 {
        match self {
            Self::W32P32 | Self::W32P64 => 0.40338510414585471153,
            Self::W64 => 0.46769122397215788544,
            Self::W128 => 0.46373128592889397439,
        }
    }

    /// Returns the batch size, in words, for signatures of length `k`.
    ///
    /// Batch sizes are such that the probability that a small counter
    /// overflows during a batch is negligible (see the paper). They are
    /// even, so that two-word iterations split batches exactly.
    ///
    /// # Panics
    ///
    /// Panics if `k` is not between 1 and 19.
    pub const fn batch_size(self, k: usize) -> u64 {
        assert!(k >= 1 && k <= 19);
        match self {
            Self::W32P32 | Self::W32P64 => BATCH_SIZE_32[k],
            Self::W64 => BATCH_SIZE_64[k],
            Self::W128 => BATCH_SIZE_128[k],
        }
    }

    /// Returns a description of the words examined, for the header.
    pub const fn description(self) -> &'static str {
        match self {
            Self::W32P32 => "32-bit words, the upper halves of the outputs",
            Self::W32P64 => "32-bit words, two per output (upper half first)",
            Self::W64 => "64-bit words, the full outputs",
            Self::W128 => "128-bit words, pairs of consecutive outputs",
        }
    }
}

/// The number of bits of a small counter used for the sum of Hamming weights;
/// the remaining upper 13 bits contain the number of updates.
pub const SUM_BITS: u32 = 19;

/// The state of the scan between words: the current signature and, for the
/// transitional variant of the test, the bit preceding the next word.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SigState {
    /// The signature of the last *k* words, as a base-3 numeral whose most
    /// significant digit is the trit of the most recent word.
    ///
    /// The paper uses the opposite order, but since the transform *T*ₖ is a
    /// Kronecker power, reversing the trits of its indices permutes its inputs
    /// and outputs in the same way, so the results do not change. Printed
    /// signatures have the trit of the oldest word on the left, as in the
    /// paper.
    pub sig: u32,
    /// The last bit of the previous word (always zero if not testing
    /// transitions).
    pub t: u64,
}

impl SigState {
    /// Returns the initial state: the all-ones signature (the most probable)
    /// and a zero bit.
    pub const fn initial(size: u32) -> Self {
        Self {
            sig: (size - 1) / 2,
            t: 0,
        }
    }
}

/// Divides by three using a fixed-point multiplication (exact for the values
/// of signatures up to *k* = 19).
#[inline(always)]
pub const fn div3(x: u32) -> u32 {
    ((x as u64 * 1431655766) >> 32) as u32
}

/// Returns the trit associated with Hamming weight `bc` of a `W`-bit word.
#[inline(always)]
const fn trit<const W: usize>(bc: u32) -> u32 {
    match W {
        32 => (bc >= 15) as u32 + (bc >= 18) as u32,
        64 => (bc >= 30) as u32 + (bc >= 35) as u32,
        _ => (bc >= 61) as u32 + (bc >= 68) as u32,
    }
}

/// Accounts for a word of Hamming weight `bc`: if `COUNT` is true, adds `bc`
/// to the sum of the small counter of the current signature and increments
/// its count (with a single addition); then, updates the signature.
#[inline(always)]
fn account<const W: usize, const COUNT: bool>(cs: &mut [u32], sig: &mut u32, bc: u32, third: u32) {
    if COUNT {
        // SAFETY: scan checks that cs.len() is 3 · third and that the initial
        // signature is smaller. Since div3 is exact on signatures, and trits
        // are at most 2, the update below preserves the invariant.
        let c = unsafe { cs.get_unchecked_mut(*sig as usize) };
        *c = c.wrapping_add(bc + (1 << SUM_BITS));
    }
    *sig = div3(*sig) + trit::<W>(bc) * third;
}

/// Performs `iters` iterations of the main loop starting from state `st`,
/// updating `st` and, if `COUNT` is true, the small counters `cs`.
///
/// `W` is the size of the words, `P64` is true if the generator has 64 bits of
/// output (relevant only if `W` is 32), `TRANS` is true if we are testing
/// transitions, and `third` is 3ᵏ⁻¹.
///
/// Returns the sum of the Hamming weights of the words if `W` is 128 (as in
/// this case it is not guaranteed that the sums in the small counters do not
/// overflow), and zero otherwise.
#[inline(never)]
pub fn scan<const W: usize, const P64: bool, const TRANS: bool, const COUNT: bool>(
    prng: &mut Prng,
    cs: &mut [u32],
    st: &mut SigState,
    iters: u64,
    third: u32,
) -> u64 {
    if COUNT {
        // Guarantees the safety of account().
        assert!(third <= 3u32.pow(18));
        assert_eq!(cs.len(), 3 * third as usize);
        assert!(st.sig < 3 * third);
    }
    let mut sig = st.sig;
    let mut t = st.t;
    let mut tot_sums = 0;

    for _ in 0..iters {
        match W {
            32 => {
                let w64 = prng.next_u64();
                // The upper 32 bits; for 32-bit generators, the whole output.
                let w = (w64 >> 32) as u32;
                let bc = if TRANS {
                    let bc = (w ^ (w << 1) ^ t as u32).count_ones();
                    t = (w >> 31) as u64;
                    bc
                } else {
                    w.count_ones()
                };
                account::<32, COUNT>(cs, &mut sig, bc, third);

                if P64 {
                    let w = w64 as u32;
                    let bc = if TRANS {
                        let bc = (w ^ (w << 1) ^ t as u32).count_ones();
                        t = (w >> 31) as u64;
                        bc
                    } else {
                        w.count_ones()
                    };
                    account::<32, COUNT>(cs, &mut sig, bc, third);
                }
            }
            64 => {
                let w = prng.next_u64();
                let bc = if TRANS {
                    let bc = (w ^ (w << 1) ^ t).count_ones();
                    t = w >> 63;
                    bc
                } else {
                    w.count_ones()
                };
                account::<64, COUNT>(cs, &mut sig, bc, third);
            }
            _ => {
                // The first output provides the lower 64 bits (a single
                // 128-bit population count is faster on some architectures).
                let w0 = prng.next_u64();
                let w1 = prng.next_u64();
                let w = w0 as u128 | (w1 as u128) << 64;
                let bc = if TRANS {
                    let bc = (w ^ (w << 1) ^ t as u128).count_ones();
                    t = (w >> 127) as u64;
                    bc
                } else {
                    w.count_ones()
                };
                tot_sums += bc as u64;
                account::<128, COUNT>(cs, &mut sig, bc, third);
            }
        }
    }

    st.sig = sig;
    st.t = t;
    tot_sums
}

/// Calls the specialization of [`scan`] for the given mode and options.
///
/// If `count` is false, `cs` is not accessed and can be empty.
#[allow(clippy::too_many_arguments)]
pub fn scan_dispatch(
    mode: Mode,
    trans: bool,
    count: bool,
    prng: &mut Prng,
    cs: &mut [u32],
    st: &mut SigState,
    iters: u64,
    third: u32,
) -> u64 {
    macro_rules! go {
        ($w:literal, $p64:literal) => {
            match (trans, count) {
                (false, false) => scan::<$w, $p64, false, false>(prng, cs, st, iters, third),
                (false, true) => scan::<$w, $p64, false, true>(prng, cs, st, iters, third),
                (true, false) => scan::<$w, $p64, true, false>(prng, cs, st, iters, third),
                (true, true) => scan::<$w, $p64, true, true>(prng, cs, st, iters, third),
            }
        };
    }
    match mode {
        Mode::W32P32 => go!(32, false),
        Mode::W32P64 => go!(32, true),
        Mode::W64 => go!(64, true),
        Mode::W128 => go!(128, true),
    }
}

/// A large counter: the number of words following a signature, and the sum of
/// the differences between their Hamming weights and *w* / 2.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Pod, Zeroable)]
pub struct CountSum {
    pub c: u64,
    pub s: i64,
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

    #[test]
    fn test_div3() {
        // Exact on all signatures up to k = 19.
        for x in (0..3u32.pow(19)).step_by(997).chain([3u32.pow(19) - 1]) {
            assert_eq!(div3(x), x / 3, "div3({x})");
        }
    }

    #[test]
    fn test_trit() {
        assert_eq!(
            (0..=32).map(trit::<32>).collect::<Vec<_>>(),
            [vec![0; 15], vec![1; 3], vec![2; 15]].concat()
        );
        assert_eq!(
            (0..=64).map(trit::<64>).collect::<Vec<_>>(),
            [vec![0; 30], vec![1; 5], vec![2; 30]].concat()
        );
        assert_eq!(
            (0..=128).map(trit::<128>).collect::<Vec<_>>(),
            [vec![0; 61], vec![1; 7], vec![2; 61]].concat()
        );
    }

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
