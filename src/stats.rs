/*
 * SPDX-FileCopyrightText: 2004-2016 David Blackman
 * SPDX-FileCopyrightText: 2017-2026 David Blackman and Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

//! Analysis of the large counters: normalization, the transform *T*ₖ, and the
//! computation of *p*-values.
//!
//! All floating-point computations are performed in the same order as in the
//! original C implementation, and use the same functions of the C math
//! library, so that the results are bit-for-bit identical. Parallelism is
//! used only where it cannot change the results.

use std::io::{self, Write};

use rayon::prelude::*;

use crate::counters::CountSum;

unsafe extern "C" {
    /// The complementary error function of the C math library (the Rust
    /// standard library does not provide it on stable).
    safe fn erfc(x: f64) -> f64;
}

/// Formats a *p*-value with full precision.
pub fn format_p_value(p: f64) -> String {
    if p == 0.0 {
        "0".to_string()
    } else if p == 1.0 {
        "1".to_string()
    } else {
        format!("{p:?}")
    }
}

/// 1 / √2 (`M_SQRT1_2`).
const SQRT1_2: f64 = std::f64::consts::FRAC_1_SQRT_2;
/// 1 / √3 (the literal of the C implementation).
#[allow(clippy::excessive_precision)]
const CORRECT3: f64 = 0.57735026918962576451;
/// 1 / √6 (the literal of the C implementation).
#[allow(clippy::excessive_precision)]
const CORRECT6: f64 = 0.40824829046386301636;

/// Below this size, transforms and scans are not parallelized.
const PAR_THRESHOLD: usize = 1 << 16;

/// Returns the probability that the smallest of `n` independent uniform
/// values in [0 . . 1) is at most `x`, that is, 1 − (1 − `x`)ⁿ.
pub fn pco_scale(x: f64, n: f64) -> f64 {
    if x >= 1.0 || x <= 0.0 {
        return x;
    }
    // Accurate for small values of x, as 1.0 - (1.0 - x).powf(n) is not.
    -((-x).ln_1p() * n).exp_m1()
}

/// Applies in place the transform *T*ₖ to `v`, whose length must be 3 ·
/// `sig`, where `sig` is a power of three.
///
/// The transform multiplies `v` by the *k*-th Kronecker power of a 3 × 3
/// orthonormal matrix, recursively as in the fast Walsh–Hadamard transform.
pub fn mix3(v: &mut [f64], sig: usize) {
    debug_assert_eq!(v.len(), 3 * sig);
    let (v0, rest) = v.split_at_mut(sig);
    let (p1, p2) = rest.split_at_mut(sig);

    fn butterfly(v0: &mut [f64], p1: &mut [f64], p2: &mut [f64]) {
        for ((x, y), z) in v0.iter_mut().zip(p1.iter_mut()).zip(p2.iter_mut()) {
            let (a, b, c) = (*x, *y, *z);
            *x = (a + b + c) * CORRECT3;
            *y = (a - c) * SQRT1_2;
            *z = (2.0 * b - a - c) * CORRECT6;
        }
    }

    if sig < PAR_THRESHOLD {
        butterfly(v0, p1, p2);
    } else {
        const CHUNK: usize = 1 << 14;
        v0.par_chunks_mut(CHUNK)
            .zip(p1.par_chunks_mut(CHUNK))
            .zip(p2.par_chunks_mut(CHUNK))
            .for_each(|((x, y), z)| butterfly(x, y, z));
    }

    let sig = sig / 3;
    if sig > 0 {
        if 3 * sig < PAR_THRESHOLD {
            mix3(v0, sig);
            mix3(p1, sig);
            mix3(p2, sig);
        } else {
            rayon::join(
                || mix3(v0, sig),
                || rayon::join(|| mix3(p1, sig), || mix3(p2, sig)),
            );
        }
    }
}

/// Writes `sig` in base 3, least significant digit first (so the most recent
/// trit is the rightmost one).
fn sig_string(mut sig: u32, k: usize) -> String {
    (0..k)
        .map(|_| {
            let d = (b'0' + (sig % 3) as u8) as char;
            sig /= 3;
            d
        })
        .collect()
}

/// For each category, the largest absolute value, the first signature at
/// which it is attained, and the size of the category.
#[derive(Clone)]
struct Extremes {
    sigma: Vec<f64>,
    sig: Vec<u32>,
    count: Vec<u32>,
}

impl Extremes {
    fn new(numcats: usize) -> Self {
        Self {
            sigma: vec![f64::MIN_POSITIVE; numcats],
            sig: vec![0; numcats],
            count: vec![0; numcats],
        }
    }

    /// Combines with the extremes of a range of signatures following those
    /// of `self`, preferring earlier signatures on ties, as a sequential scan
    /// does.
    fn merge(mut self, other: Self) -> Self {
        for c in 0..self.sigma.len() {
            if other.sigma[c] > self.sigma[c] {
                self.sigma[c] = other.sigma[c];
                self.sig[c] = other.sig[c];
            }
            self.count[c] += other.count[c];
        }
        self
    }
}

/// Scans the transformed values of signatures in [`start` . . `start` +
/// `norm.len()`), skipping signature zero.
///
/// The category of a signature is the number of its nonzero trits, capped
/// at the number of categories, minus one.
fn scan_extremes(norm: &[f64], start: usize, k: usize, numcats: usize) -> Extremes {
    let mut e = Extremes::new(numcats);
    // Trits of the current signature, least significant first, and the number
    // of nonzero ones.
    let mut trits = vec![0u8; k];
    let mut nonzero = 0;
    let mut s = start;
    for t in trits.iter_mut() {
        *t = (s % 3) as u8;
        nonzero += (*t != 0) as usize;
        s /= 3;
    }
    for (i, &x) in (start..).zip(norm) {
        if i != 0 {
            let c = nonzero.min(numcats) - 1;
            e.count[c] += 1;
            let x = x.abs();
            if x > e.sigma[c] {
                e.sig[c] = i as u32;
                e.sigma[c] = x;
            }
        }
        // Increment the signature.
        for t in trits.iter_mut() {
            if *t == 2 {
                *t = 0;
                nonzero -= 1;
            } else {
                *t += 1;
                nonzero += (*t == 1) as usize;
                break;
            }
        }
    }
    e
}

/// Normalizes the large counters into `norm`, applies the transform, and
/// writes and returns the resulting *p*-value.
///
/// For each category of signatures, we print the most extreme value and its
/// *p*-value, corrected for the number of signatures in the category. The
/// final *p*-value is the smallest category *p*-value, corrected for the
/// number of categories.
pub fn compute_pvalue(
    count_sum: &[CountSum],
    norm: &mut [f64],
    k: usize,
    numcats: usize,
    word_bits: usize,
    trans: bool,
    out: &mut impl Write,
) -> io::Result<f64> {
    let size = norm.len();
    // The variance of the Hamming weight of a w-bit word.
    let db = word_bits as f64 * 0.25;

    let normalize = |(n, cs): (&mut f64, &CountSum)| {
        // We expect mean 0 and standard deviation 1.
        *n = if cs.c == 0 {
            0.0
        } else {
            cs.s as f64 / (cs.c as f64 * db).sqrt()
        };
    };
    if size < PAR_THRESHOLD {
        norm.iter_mut().zip(count_sum).for_each(normalize);
    } else {
        norm.par_iter_mut().zip(count_sum).for_each(normalize);
    }

    // After the transform we still expect mean 0 and standard deviation 1
    // (except for element 0, which we ignore).
    mix3(norm, size / 3);

    const CHUNK: usize = 1 << 16;
    let e = if size < PAR_THRESHOLD {
        scan_extremes(norm, 0, k, numcats)
    } else {
        norm.par_chunks(CHUNK)
            .enumerate()
            .map(|(i, chunk)| scan_extremes(chunk, i * CHUNK, k, numcats))
            .reduce_with(Extremes::merge)
            .expect("at least one chunk")
    };

    let mut overall_pvalue = f64::MAX;
    for i in 0..numcats {
        // Convert the absolute value of an approximately normal value into a
        // p-value, and correct for having picked the smallest one.
        let pvalue = pco_scale(erfc(SQRT1_2 * e.sigma[i]), e.count[i] as f64);
        writeln!(
            out,
            "mix3 extreme = {:.5} (sig = {}) weight {}{} ({}), p-value = {}",
            e.sigma[i],
            sig_string(e.sig[i], k),
            if i == numcats - 1 { ">=" } else { "" },
            i + 1,
            e.count[i],
            format_p_value(pvalue)
        )?;
        if pvalue < overall_pvalue {
            overall_pvalue = pvalue;
        }
    }

    writeln!(
        out,
        "bits per word = {} (analyzing {}); min category p-value = {}\n",
        word_bits,
        if trans { "transitions" } else { "bits" },
        format_p_value(overall_pvalue)
    )?;
    // Again, we picked the smallest of numcats p-values.
    Ok(pco_scale(overall_pvalue, numcats as f64))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Applies the transform by explicit multiplication by the Kronecker
    /// power. The order of operations is different from that of the fast
    /// transform, so results must be compared with a tolerance.
    fn naive_transform(v: &[f64], k: usize) -> Vec<f64> {
        // The base matrix of the paper, but for the sign of the third column,
        // which is negated in the C implementation (only absolute values are
        // used, so the sign is irrelevant).
        let m = [
            [CORRECT3, SQRT1_2, -CORRECT6],
            [CORRECT3, 0.0, 2.0 * CORRECT6],
            [CORRECT3, -SQRT1_2, -CORRECT6],
        ];
        let size = 3usize.pow(k as u32);
        (0..size)
            .map(|j| {
                (0..size)
                    .map(|i| {
                        // Entry (i, j) of the Kronecker power: product over
                        // trit positions.
                        let (mut a, mut b, mut prod) = (i, j, 1.0);
                        for _ in 0..k {
                            prod *= m[a % 3][b % 3];
                            a /= 3;
                            b /= 3;
                        }
                        v[i] * prod
                    })
                    .sum()
            })
            .collect()
    }

    #[test]
    fn test_mix3_is_kronecker_power() {
        for k in 1..=5 {
            let size = 3usize.pow(k as u32);
            let v: Vec<f64> = (0..size)
                .map(|i| ((i * 7919) % 101) as f64 - 50.0)
                .collect();
            let mut w = v.clone();
            mix3(&mut w, size / 3);
            for (x, y) in w.iter().zip(naive_transform(&v, k)) {
                assert!((x - y).abs() < 1e-9, "k = {k}: {x} != {y}");
            }
        }
    }

    #[test]
    fn test_mix3_parallel_is_exact() {
        // Large enough to use the parallel branches; compare with a
        // sequential recursion written as in the C implementation.
        fn seq(v: &mut [f64], sig: usize) {
            for i in 0..sig {
                let (a, b, c) = (v[i], v[i + sig], v[i + 2 * sig]);
                v[i] = (a + b + c) * CORRECT3;
                v[i + sig] = (a - c) * SQRT1_2;
                v[i + 2 * sig] = (2.0 * b - a - c) * CORRECT6;
            }
            let s = sig / 3;
            if s > 0 {
                seq(&mut v[..sig], s);
                seq(&mut v[sig..2 * sig], s);
                seq(&mut v[2 * sig..], s);
            }
        }
        let size = 3usize.pow(12);
        let v: Vec<f64> = (0..size)
            .map(|i| ((i * 7919) % 1009) as f64 / 3.0)
            .collect();
        let (mut a, mut b) = (v.clone(), v);
        mix3(&mut a, size / 3);
        seq(&mut b, size / 3);
        assert_eq!(a, b);
    }

    #[test]
    fn test_scan_extremes_chunks() {
        let k = 9;
        let size = 3usize.pow(k as u32);
        let v: Vec<f64> = (0..size)
            .map(|i| (((i * 7919) % 1013) as f64 - 500.0) / 100.0)
            .collect();
        let whole = scan_extremes(&v, 0, k, 5);
        let merged = v
            .chunks(1000)
            .enumerate()
            .map(|(i, c)| scan_extremes(c, i * 1000, k, 5))
            .reduce(Extremes::merge)
            .unwrap();
        assert_eq!(whole.sigma, merged.sigma);
        assert_eq!(whole.sig, merged.sig);
        assert_eq!(whole.count, merged.count);
        // Category sizes: signatures with j nonzero trits are C(k, j) 2ʲ.
        assert_eq!(whole.count, vec![18, 144, 672, 2016, 16832]);
    }

    #[test]
    fn test_sig_string() {
        assert_eq!(sig_string(0, 4), "0000");
        assert_eq!(sig_string(1, 4), "1000");
        assert_eq!(sig_string(3 * 3 * 3 * 2, 4), "0002");
        assert_eq!(sig_string((81 - 1) / 2, 4), "1111");
    }
}
