/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

//! Arbitrary jumps for 𝐅₂-linear generators.
//!
//! This is a port of `f2x.c` from <https://prng.di.unimi.it/>: calling the
//! next-state function of a linear generator once is equivalent to
//! multiplying the state by *x* in the ring 𝐅₂\[*x*\] / (*p*), where *p* is the
//! characteristic polynomial of the transition map; thus, jumping ahead by *n*
//! steps is equivalent to multiplying by *x*ⁿ mod *p*, which can be computed
//! by square-and-multiply and applied to the state with an
//! accumulate-and-step loop.
//!
//! Polynomials have degree 64 · `N` and are represented by their 64 · `N`
//! low-order coefficients (coefficient *k* is bit *k* % 64 of word ⌊*k* /
//! 64⌋); the leading coefficient of characteristic polynomials is implicit.

/// A linear generator whose state is represented by `N` words.
pub trait LinearGenerator<const N: usize>: Copy {
    /// Advances the state by one step.
    fn step(&mut self);
    /// Returns the state.
    ///
    /// The representation must be canonical: generators using a rotating
    /// index into an array must rotate the array so that the index is zero.
    fn to_vector(&self) -> [u64; N];
    /// Sets the state.
    fn set_vector(&mut self, v: [u64; N]);
}

/// `a` ← `a` · *x* mod `charpoly`.
fn mul_x<const N: usize>(a: &mut [u64; N], charpoly: &[u64; N]) {
    let mut carry = 0;
    for w in a.iter_mut() {
        let next_carry = *w >> 63;
        *w = (*w << 1) | carry;
        carry = next_carry;
    }
    // The coefficient of x^(64N) after the shift.
    if carry != 0 {
        for (w, c) in a.iter_mut().zip(charpoly) {
            *w ^= c;
        }
    }
}

/// Returns `a` · `b` mod `charpoly` (Horner's method over the bits of `a`,
/// from the top).
fn mul_mod<const N: usize>(a: &[u64; N], b: &[u64; N], charpoly: &[u64; N]) -> [u64; N] {
    let mut r = [0; N];
    for k in (0..64 * N).rev() {
        mul_x(&mut r, charpoly);
        if (a[k / 64] >> (k % 64)) & 1 != 0 {
            for (w, c) in r.iter_mut().zip(b) {
                *w ^= c;
            }
        }
    }
    r
}

/// Returns *x*ⁿ mod `charpoly` (square-and-multiply over the bits of *n*,
/// from the most significant down).
pub fn x_pow_mod<const N: usize>(n: u64, charpoly: &[u64; N]) -> [u64; N] {
    let mut r = [0; N];
    r[0] = 1;
    for k in (0..64 - n.leading_zeros()).rev() {
        r = mul_mod(&r, &r, charpoly);
        if (n >> k) & 1 != 0 {
            mul_x(&mut r, charpoly);
        }
    }
    r
}

/// Advances a generator by `n` steps, given the characteristic polynomial of
/// its transition map.
pub fn jump<const N: usize>(g: &mut impl LinearGenerator<N>, n: u64, charpoly: &[u64; N]) {
    if n < 64 * N as u64 {
        // Cheaper than the accumulate-and-step loop.
        for _ in 0..n {
            g.step();
        }
        return;
    }
    let poly = x_pow_mod(n, charpoly);
    let mut acc = [0; N];
    for w in poly {
        for b in 0..64 {
            if (w >> b) & 1 != 0 {
                for (a, s) in acc.iter_mut().zip(g.to_vector()) {
                    *a ^= s;
                }
            }
            g.step();
        }
    }
    g.set_vector(acc);
}

/// Returns the minimal polynomial (the polynomial *m* of smallest degree such
/// that Σᵢ *m*ᵢ *s*ₙ₊ᵢ = 0 for all *n*) of a binary sequence, including the
/// leading coefficient, computed by the Berlekamp–Massey algorithm.
///
/// The result is correct if the sequence is at least twice as long as the
/// degree of the minimal polynomial. It is used to check characteristic
/// polynomials: if the minimal polynomial of the sequence of a bit of the
/// state has full degree, it is the characteristic polynomial.
#[cfg(test)]
pub(crate) fn minimal_polynomial(s: &[bool]) -> Vec<u64> {
    let bit = |p: &[u64], i: usize| (p[i / 64] >> (i % 64)) & 1 != 0;
    let words = s.len() / 64 + 2;
    // c is the connection polynomial, b the previous one.
    let (mut c, mut b) = (vec![0u64; words], vec![0u64; words]);
    c[0] = 1;
    b[0] = 1;
    let (mut l, mut shift) = (0, 1);
    for i in 0..s.len() {
        // Discrepancy between s[i] and the prediction of c.
        let mut d = s[i];
        for j in 1..=l {
            d ^= bit(&c, j) & s[i - j];
        }
        if d {
            let t = c.clone();
            for j in 0..s.len() + 1 - shift {
                if bit(&b, j) {
                    c[(j + shift) / 64] ^= 1 << ((j + shift) % 64);
                }
            }
            if 2 * l <= i {
                l = i + 1 - l;
                b = t;
                shift = 1;
                continue;
            }
        }
        shift += 1;
    }
    // The minimal polynomial is the reciprocal of the connection polynomial.
    let mut m = vec![0u64; l / 64 + 1];
    for j in (0..=l).filter(|&j| bit(&c, j)) {
        m[(l - j) / 64] |= 1 << ((l - j) % 64);
    }
    m
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A 64-bit xorshift generator, whose characteristic polynomial is
    /// primitive.
    #[derive(Clone, Copy)]
    struct Xorshift64(u64);

    impl LinearGenerator<1> for Xorshift64 {
        fn step(&mut self) {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
        }
        fn to_vector(&self) -> [u64; 1] {
            [self.0]
        }
        fn set_vector(&mut self, v: [u64; 1]) {
            self.0 = v[0];
        }
    }

    #[test]
    fn test_minimal_polynomial_fibonacci() {
        // The sequence 0, 1, 1, 0, 1, 1, ... satisfies sₙ₊₂ = sₙ₊₁ + sₙ.
        let s: Vec<bool> = (0..20).map(|i| i % 3 != 0).collect();
        assert_eq!(minimal_polynomial(&s), vec![0b111]);
    }

    #[test]
    fn test_jump_matches_steps() {
        let mut g = Xorshift64(1);
        let bits: Vec<bool> = (0..128)
            .map(|_| {
                g.step();
                g.0 & 1 != 0
            })
            .collect();
        let m = minimal_polynomial(&bits);
        // Full degree: the leading coefficient is in the second word.
        assert_eq!(m[1], 1);
        let charpoly = [m[0]];
        for n in [0u64, 1, 63, 64, 65, 1000, 123_456] {
            let mut a = Xorshift64(0x0123_4567_89ab_cdef);
            let mut b = a;
            jump(&mut a, n, &charpoly);
            for _ in 0..n {
                b.step();
            }
            assert_eq!(a.0, b.0, "jump({n}) disagrees with stepping");
        }
    }
}
