/*
 * SPDX-FileCopyrightText: 2004-2016 David Blackman
 * SPDX-FileCopyrightText: 2017-2026 David Blackman and Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

//! The main loop of the test, which scans the output of the generator and
//! updates the small counters.

use crate::mode::Mode;
use crate::prng::Prng;

/// The number of bits of a small counter used for the sum of Hamming weights;
/// the remaining upper 13 bits contain the number of updates.
pub const SUM_BITS: u32 = 19;

/// The state of the scan between words: the current signature and, for the
/// transitional variant of the test, the bit preceding the next word.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SigState {
    /// The signature of the last *k* words, as a base-3 numeral whose most
    /// significant digit is the trit of the most recent word.
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
}
