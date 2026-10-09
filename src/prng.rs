/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

//! Pseudorandom number generators selected at build time via Cargo features.
//!
//! Exactly one feature must be enabled when building the crate; each variant
//! exposes a single [`Prng`] type with a `NAME` constant (a description), a
//! `BITS` constant (the number of bits of output, 32 or 64), a `new(seed: u64)
//! -> Self` constructor, a `next_u64(&mut self) -> u64` step function, and a
//! `try_skip(&mut self, n: u64) -> Result<(), ()>` method to skip ahead by `n`
//! steps.
//!
//! `try_skip(n)` must have the same effect as `n` calls to `next_u64`. It must
//! succeed for every `n`, or fail for every `n` (including 0) without changing
//! the state: the parallel runners probe the capability with `try_skip(0)`.
//!
//! Generators outputting less than 64 bits must shift their outputs to the top
//! (e.g., 32-bit generators must return their output shifted to the left by 32).
//! `BITS` is the default of `--prng-bits`.

// When no PRNG feature is selected, the _prng marker (enabled by every PRNG
// feature in Cargo.toml) is off. Emit one clear error AND expose a placeholder Prng
// so the rest of the crate still type-checks.
#[cfg(not(feature = "_prng"))]
compile_error!(
    "no PRNG selected: enable exactly one PRNG feature, as in `--features xoroshiro128plus` \
     (see the [features] table in Cargo.toml)"
);

#[cfg(not(feature = "_prng"))]
mod placeholder {
    #[derive(Clone, Copy)]
    pub struct Prng;
    impl Prng {
        pub const NAME: &str = "(no generator selected)";
        pub const BITS: usize = 64;
        pub fn new(_seed: u64) -> Self {
            Self
        }
        #[inline(always)]
        pub fn next_u64(&mut self) -> u64 {
            0
        }
        pub fn try_skip(&mut self, _n: u64) -> Result<(), ()> {
            Err(())
        }
    }
}

#[cfg(not(feature = "_prng"))]
pub use placeholder::Prng;

// ----- Trivial counters (sanity/baseline) ---------------------------------------------

#[cfg(feature = "incr")]
#[derive(Clone, Copy)]
pub struct Prng {
    x: u64,
}

#[cfg(feature = "incr")]
impl Prng {
    pub const NAME: &str = "incr (counter, x += 1)";
    pub const BITS: usize = 64;
    pub fn new(seed: u64) -> Self {
        Self { x: seed }
    }

    #[inline(always)]
    pub fn next_u64(&mut self) -> u64 {
        self.x = self.x.wrapping_add(1);
        self.x
    }

    pub fn try_skip(&mut self, n: u64) -> Result<(), ()> {
        self.x = self.x.wrapping_add(n);
        Ok(())
    }
}

// ----- 𝐅₂-linear generators -----------------------------------------------------------
//
// The generators of "A New Test for Hamming-Weight Dependencies", by David
// Blackman and Sebastiano Vigna, with the parameters used in the paper. As in
// the replication code of the paper, the linear engines return the word s0 of
// the reference implementation of the corresponding `+` variant, which returns
// s0 + s1. Since the state transition is 𝐅₂-linear, try_skip uses its
// characteristic polynomial (see the f2 module).

/// Returns the next output of SplitMix64 (incrementing first), used to fill the
/// state of generators with more than 64 bits of state. Since the output
/// function of SplitMix64 is a bijection, two consecutive outputs cannot be
/// both zero, so the state of the 64-bit generators is never zero; WELL512a
/// keeps only the upper 32 bits of each output, so its state is zero only with
/// negligible probability.
fn splitmix64(x: &mut u64) -> u64 {
    *x = x.wrapping_add(0x9e3779b97f4a7c15);
    let mut z = *x;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}

// Implements try_skip using the characteristic polynomial CHARPOLY of the
// generator, which must implement LinearGenerator.
#[cfg(any(
    feature = "xorshift128",
    feature = "xorshift128plus",
    feature = "xorshift1024",
    feature = "xorshift1024plus",
    feature = "xoroshiro128",
    feature = "xoroshiro128plus",
    feature = "xoroshiro1024",
    feature = "xoroshiro1024plus",
    feature = "well512a"
))]
macro_rules! f2_try_skip {
    () => {
        pub fn try_skip(&mut self, n: u64) -> Result<(), ()> {
            $crate::f2::jump(self, n, &CHARPOLY);
            Ok(())
        }
    };
}

#[cfg(any(feature = "xorshift128", feature = "xorshift128plus"))]
macro_rules! xorshift128 {
    ($plus:literal) => {
        const A: u32 = 23;
        const B: u32 = 18;
        const C: u32 = 5;
        // The characteristic polynomial of the transition map, without the
        // leading coefficient.
        const CHARPOLY: [u64; 2] = [0x024f06fae9e61daf, 0x2844c5d42caf7db0];

        #[derive(Clone, Copy)]
        pub struct Prng {
            s: [u64; 2],
        }

        impl Prng {
            pub const NAME: &str = if $plus {
                "xorshift128+ (23, 18, 5)"
            } else {
                "xorshift128 (23, 18, 5)"
            };
            pub const BITS: usize = 64;

            pub fn new(seed: u64) -> Self {
                let mut x = seed;
                Self {
                    s: [splitmix64(&mut x), splitmix64(&mut x)],
                }
            }

            #[inline(always)]
            pub fn next_u64(&mut self) -> u64 {
                let mut s1 = self.s[0];
                let s0 = self.s[1];
                let result = if $plus { s0.wrapping_add(s1) } else { s0 };
                self.s[0] = s0;
                s1 ^= s1 << A;
                self.s[1] = s1 ^ s0 ^ (s1 >> B) ^ (s0 >> C);
                result
            }

            f2_try_skip!();
        }

        impl $crate::f2::LinearGenerator<2> for Prng {
            fn step(&mut self) {
                self.next_u64();
            }
            fn to_vector(&self) -> [u64; 2] {
                self.s
            }
            fn set_vector(&mut self, v: [u64; 2]) {
                self.s = v;
            }
        }
    };
}

#[cfg(feature = "xorshift128")]
xorshift128!(false);

#[cfg(feature = "xorshift128plus")]
xorshift128!(true);

#[cfg(any(feature = "xorshift1024", feature = "xorshift1024plus"))]
macro_rules! xorshift1024 {
    ($plus:literal) => {
        const A: u32 = 31;
        const B: u32 = 11;
        const C: u32 = 30;
        // The characteristic polynomial of the transition map, without the
        // leading coefficient.
        const CHARPOLY: [u64; 16] = [
            0x1000000000000001,
            0x2200aa001400f000,
            0x0111e1c02bc18180,
            0x030d535201556130,
            0x4a32d044029b08f7,
            0x34b3216457d7b028,
            0xe860f083d70158c6,
            0xdf6a7cadba32bca9,
            0xbabab341e2554b59,
            0xcd40a7e2537771ea,
            0x0040f0e46e848800,
            0xa1422cb7814f5c68,
            0x53116c08605c805f,
            0x0440024003007b28,
            0x787878786d381540,
            0x0000000000007879,
        ];

        #[derive(Clone, Copy)]
        pub struct Prng {
            s: [u64; 16],
            p: usize,
        }

        impl Prng {
            pub const NAME: &str = if $plus {
                "xorshift1024+ (31, 11, 30)"
            } else {
                "xorshift1024 (31, 11, 30)"
            };
            pub const BITS: usize = 64;

            pub fn new(seed: u64) -> Self {
                let mut x = seed;
                Self {
                    s: ::std::array::from_fn(|_| splitmix64(&mut x)),
                    p: 0,
                }
            }

            #[inline(always)]
            pub fn next_u64(&mut self) -> u64 {
                let s0 = self.s[self.p];
                self.p = (self.p + 1) & 15;
                let mut s1 = self.s[self.p];
                let result = if $plus { s0.wrapping_add(s1) } else { s0 };
                s1 ^= s1 << A;
                self.s[self.p] = s1 ^ s0 ^ (s1 >> B) ^ (s0 >> C);
                result
            }

            f2_try_skip!();
        }

        impl $crate::f2::LinearGenerator<16> for Prng {
            fn step(&mut self) {
                self.next_u64();
            }
            fn to_vector(&self) -> [u64; 16] {
                ::std::array::from_fn(|j| self.s[(self.p + j) & 15])
            }
            fn set_vector(&mut self, v: [u64; 16]) {
                for (j, x) in v.into_iter().enumerate() {
                    self.s[(self.p + j) & 15] = x;
                }
            }
        }
    };
}

#[cfg(feature = "xorshift1024")]
xorshift1024!(false);

#[cfg(feature = "xorshift1024plus")]
xorshift1024!(true);

#[cfg(any(feature = "xoroshiro128", feature = "xoroshiro128plus"))]
macro_rules! xoroshiro128 {
    ($plus:literal) => {
        const A: u32 = 24;
        const B: u32 = 16;
        const C: u32 = 37;
        // The characteristic polynomial of the transition map, without the
        // leading coefficient.
        const CHARPOLY: [u64; 2] = [0x095b8f76579aa001, 0x0008828e513b43d5];

        #[derive(Clone, Copy)]
        pub struct Prng {
            s: [u64; 2],
        }

        impl Prng {
            pub const NAME: &str = if $plus {
                "xoroshiro128+ (24, 16, 37)"
            } else {
                "xoroshiro128 (24, 16, 37)"
            };
            pub const BITS: usize = 64;

            pub fn new(seed: u64) -> Self {
                let mut x = seed;
                Self {
                    s: [splitmix64(&mut x), splitmix64(&mut x)],
                }
            }

            #[inline(always)]
            pub fn next_u64(&mut self) -> u64 {
                let s0 = self.s[0];
                let mut s1 = self.s[1];
                let result = if $plus { s0.wrapping_add(s1) } else { s0 };
                s1 ^= s0;
                self.s[0] = s0.rotate_left(A) ^ s1 ^ (s1 << B);
                self.s[1] = s1.rotate_left(C);
                result
            }

            f2_try_skip!();
        }

        impl $crate::f2::LinearGenerator<2> for Prng {
            fn step(&mut self) {
                self.next_u64();
            }
            fn to_vector(&self) -> [u64; 2] {
                self.s
            }
            fn set_vector(&mut self, v: [u64; 2]) {
                self.s = v;
            }
        }
    };
}

#[cfg(feature = "xoroshiro128")]
xoroshiro128!(false);

#[cfg(feature = "xoroshiro128plus")]
xoroshiro128!(true);

#[cfg(any(feature = "xoroshiro1024", feature = "xoroshiro1024plus"))]
macro_rules! xoroshiro1024 {
    ($plus:literal) => {
        const A: u32 = 25;
        const B: u32 = 27;
        const C: u32 = 36;
        // The characteristic polynomial of the transition map, without the
        // leading coefficient.
        const CHARPOLY: [u64; 16] = [
            0x5cfeb8cc48ddb211,
            0xb73e379d035a06dd,
            0x17d5100a20a0350e,
            0x7550223f68f98cac,
            0x29d373b5c5ed3459,
            0x3689b412ef70de48,
            0xa1d3b6ee079a7cc6,
            0x9bf0b669abd100f8,
            0x955c84e105f60997,
            0x6ca140c61889cddd,
            0xabaf68c5fc3a0e4a,
            0xa46134526b83adc5,
            0x0710704d05683d63,
            0x580d080b44b606a2,
            0x008040a0580158a1,
            0x0000000000800081,
        ];

        #[derive(Clone, Copy)]
        pub struct Prng {
            s: [u64; 16],
            p: usize,
        }

        impl Prng {
            pub const NAME: &str = if $plus {
                "xoroshiro1024+ (25, 27, 36)"
            } else {
                "xoroshiro1024 (25, 27, 36)"
            };
            pub const BITS: usize = 64;

            pub fn new(seed: u64) -> Self {
                let mut x = seed;
                Self {
                    s: ::std::array::from_fn(|_| splitmix64(&mut x)),
                    p: 0,
                }
            }

            #[inline(always)]
            pub fn next_u64(&mut self) -> u64 {
                let q = self.p;
                self.p = (self.p + 1) & 15;
                let s0 = self.s[self.p];
                let mut s15 = self.s[q];
                let result = if $plus { s0.wrapping_add(s15) } else { s0 };
                s15 ^= s0;
                self.s[q] = s0.rotate_left(A) ^ s15 ^ (s15 << B);
                self.s[self.p] = s15.rotate_left(C);
                result
            }

            f2_try_skip!();
        }

        impl $crate::f2::LinearGenerator<16> for Prng {
            fn step(&mut self) {
                self.next_u64();
            }
            fn to_vector(&self) -> [u64; 16] {
                ::std::array::from_fn(|j| self.s[(self.p + j) & 15])
            }
            fn set_vector(&mut self, v: [u64; 16]) {
                for (j, x) in v.into_iter().enumerate() {
                    self.s[(self.p + j) & 15] = x;
                }
            }
        }
    };
}

#[cfg(feature = "xoroshiro1024")]
xoroshiro1024!(false);

#[cfg(feature = "xoroshiro1024plus")]
xoroshiro1024!(true);

#[cfg(feature = "well512a")]
// The characteristic polynomial of the transition map, without the
// leading coefficient.
const CHARPOLY: [u64; 8] = [
    0xe0f4f3e2a7600001,
    0x7d6b79a9cb30e185,
    0x13a524cbf3d46237,
    0xa1381bcb38e3c2d2,
    0x04a72cdaf7ab5f06,
    0xaca072f14e302521,
    0x24aa25c94dd96181,
    0x0000000003c417e7,
];

#[cfg(feature = "well512a")]
#[derive(Clone, Copy)]
pub struct Prng {
    s: [u32; 16],
    i: usize,
}

#[cfg(feature = "well512a")]
impl Prng {
    pub const NAME: &str = "WELL512a";
    pub const BITS: usize = 32;

    pub fn new(seed: u64) -> Self {
        let mut x = seed;
        Self {
            s: std::array::from_fn(|_| (splitmix64(&mut x) >> 32) as u32),
            i: 0,
        }
    }

    #[inline(always)]
    pub fn next_u64(&mut self) -> u64 {
        let i = self.i;
        let z0 = self.s[(i + 15) & 15];
        let v0 = self.s[i];
        let vm1 = self.s[(i + 13) & 15];
        let vm2 = self.s[(i + 9) & 15];
        let z1 = (v0 ^ (v0 << 16)) ^ (vm1 ^ (vm1 << 15));
        let z2 = vm2 ^ (vm2 >> 11);
        let new_v1 = z1 ^ z2;
        self.s[i] = new_v1;
        self.s[(i + 15) & 15] = (z0 ^ (z0 << 2))
            ^ (z1 ^ (z1 << 18))
            ^ (z2 << 28)
            ^ (new_v1 ^ ((new_v1 << 5) & 0xda442d24));
        self.i = (i + 15) & 15;
        (self.s[self.i] as u64) << 32
    }

    f2_try_skip!();
}

#[cfg(feature = "well512a")]
impl crate::f2::LinearGenerator<8> for Prng {
    fn step(&mut self) {
        self.next_u64();
    }
    fn to_vector(&self) -> [u64; 8] {
        std::array::from_fn(|j| {
            self.s[(self.i + 2 * j) & 15] as u64 | (self.s[(self.i + 2 * j + 1) & 15] as u64) << 32
        })
    }
    fn set_vector(&mut self, v: [u64; 8]) {
        for (j, x) in v.into_iter().enumerate() {
            self.s[(self.i + 2 * j) & 15] = x as u32;
            self.s[(self.i + 2 * j + 1) & 15] = (x >> 32) as u32;
        }
    }
}

// Checks the characteristic polynomial of 𝐅₂-linear generators: the minimal
// polynomial of the sequence of the lowest bit of the state must be equal to it
// (if it has full degree, it is the characteristic polynomial).
#[cfg(all(
    test,
    any(
        feature = "xorshift128",
        feature = "xorshift128plus",
        feature = "xorshift1024",
        feature = "xorshift1024plus",
        feature = "xoroshiro128",
        feature = "xoroshiro128plus",
        feature = "xoroshiro1024",
        feature = "xoroshiro1024plus",
        feature = "well512a"
    )
))]
mod charpoly_tests {
    use super::*;
    use crate::f2::{LinearGenerator, minimal_polynomial};

    #[test]
    fn test_charpoly() {
        let n = CHARPOLY.len();
        let mut g = Prng::new(0);
        let bits: Vec<bool> = (0..128 * n)
            .map(|_| {
                let b = g.to_vector()[0] & 1 != 0;
                g.step();
                b
            })
            .collect();
        let m = minimal_polynomial(&bits);
        assert_eq!(
            (m.len(), m[m.len() - 1]),
            (n + 1, 1),
            "the minimal polynomial has not full degree"
        );
        assert_eq!(m[..n], CHARPOLY, "computed: {:#018x?}", &m[..n]);
    }
}

#[cfg(test)]
mod skip_tests {
    use super::*;

    // A single generic test that works for every feature-selected generator.
    // For skip-capable generators it checks that try_skip(n) lands on the same
    // state as n sequential next_u64() calls, and that two skips compose.
    // For generators without jump-ahead (try_skip returns Err), there is
    // nothing to verify and the test passes trivially.
    #[test]
    fn test_skip_matches_repeated_next() -> anyhow::Result<()> {
        let seed = 0x0123_4567_89ab_cdef;
        // Cheap-to-step values (sequential stepping must stay fast).
        for &n in &[0u64, 1, 2, 7, 1000, 100_000] {
            let mut a = Prng::new(seed);
            if a.try_skip(n).is_err() {
                return Ok(()); // generator has no jump-ahead; nothing to check
            }
            let mut b = Prng::new(seed);
            for _ in 0..n {
                b.next_u64();
            }
            for k in 0..64 {
                assert_eq!(
                    a.next_u64(),
                    b.next_u64(),
                    "skip({n}) disagreed with stepping at output {k}"
                );
            }
        }
        // Composition: two successive skips with large offsets must equal one
        // combined skip.
        let (x, y) = (1u64 << 40, (1u64 << 41) + 12_345);
        let skip_err = |()| anyhow::anyhow!("try_skip failed on a skip-capable generator");
        let mut p = Prng::new(seed);
        p.try_skip(x).map_err(skip_err)?;
        p.try_skip(y).map_err(skip_err)?;
        let mut q = Prng::new(seed);
        q.try_skip(x + y).map_err(skip_err)?;
        for k in 0..64 {
            assert_eq!(
                p.next_u64(),
                q.next_u64(),
                "skip composition failed at output {k}"
            );
        }
        Ok(())
    }
}
