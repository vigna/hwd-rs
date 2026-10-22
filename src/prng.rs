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
    feature = "well512a",
    feature = "mt19937"
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

// ----- Mersenne Twisters --------------------------------------------------------------
//
// The original Mersenne Twister and the SFMT with 19937 bits of state, which are
// not in the paper. As in the reference implementations, the state is an array
// of words that is regenerated in place when all of them have been returned, and
// the first output comes from the first regeneration; however, the state is
// filled with SplitMix64 from the seed, rather than by the initialization
// functions of the reference implementations. For jumps, the state is the
// sequence of words that follows the last output, which is a window sliding
// by one word at each step of the recurrence.

/// The original 32-bit Mersenne Twister, a port of `mt19937ar.c` by Makoto
/// Matsumoto and Takuji Nishimura.
#[cfg(feature = "mt19937")]
mod mt19937 {
    use super::splitmix64;

    const N: usize = 624;
    const M: usize = 397;
    const MATRIX_A: u32 = 0x9908b0df;
    const UPPER_MASK: u32 = 0x80000000;
    const LOWER_MASK: u32 = 0x7fffffff;

    // The characteristic polynomial of the transition map on 624 words, without
    // the leading coefficient. As the transition map discards the lower 31 bits
    // of the oldest word, it is x³¹ times the primitive polynomial of degree
    // 19937 of the Mersenne Twister.
    #[rustfmt::skip]
    pub(super) const CHARPOLY: [u64; N / 2] = [
        0x0000000080000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000010,
        0x0000000000000000, 0x0000000000000000, 0x0000008000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000010000, 0x0000000000000400, 0x0000000000000000,
        0x0000000000000000, 0x0000200000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000010000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0010000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000800000,
        0x0000000000000000, 0x0000000000000000, 0x0400000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000020000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000001, 0x0000000000000000,
        0x0000000000000000, 0x0000000800000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000040, 0x0000000000000000, 0x0000000000000000,
        0x0000020000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0001000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x1000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000100000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000100,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0100000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000100000, 0x4000000000000010, 0x0000000000000000,
        0x0000000000000000, 0x0000008000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000400, 0x0000000000000000, 0x0000000000000000,
        0x0000200000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000001, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000001000,
        0x0000000000010000, 0x0000000000000000, 0x0000800000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000001000000000,
        0x0000010000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000080,
        0x0000000000000000, 0x0001000000000000, 0x0010040000000000, 0x0000400000000000,
        0x0000000000000000, 0x0000000001000000, 0x0000000000802000, 0x1000000000000000,
        0x0000000000000000, 0x0400000000000000, 0x0000000000000000, 0x0000001080000000,
        0x0000000000000000, 0x0000000020000000, 0x0000000000000100, 0x1000000000000080,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000100000, 0x0000000000004000, 0x0000000000000100, 0x1080000000000000,
        0x0000000000000000, 0x0000080000000000, 0x0000000100000000, 0x0000000084000000,
        0x0000000000100000, 0x0000000000000000, 0x2000000000000100, 0x0000000000000004,
        0x0000000000000000, 0x0001000000000000, 0x0000002100000000, 0x0000000004000000,
        0x0000000000100000, 0x0000000000000000, 0x0000000000000008, 0x1080000000000000,
        0x0000000000000000, 0x0000004000000000, 0x0000000100000000, 0x0000000080000000,
        0x0000000000100000, 0x0000000000000200, 0x0000000000000000, 0x0000000000000000,
        0x0000100000000000, 0x0000004000000000, 0x0000000100000000, 0x0000000000000000,
        0x0000000000108000, 0x0000000000000000, 0x0000000000000008, 0x0084000000000000,
        0x0000100000000000, 0x0000000000000000, 0x0000000100000000, 0x0000000004200000,
        0x0000000000000000, 0x0000000000000000, 0x2100000000000100, 0x0004000000000000,
        0x0000100000000000, 0x0000000000000000, 0x0000000008000000, 0x0000000000000000,
        0x0000000000108000, 0x4000000000000000, 0x0100000000000000, 0x0080000000000000,
        0x0000100000000000, 0x0000000200000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x4000000000000010, 0x0100000000000000, 0x0000000000000000,
        0x0000108000000000, 0x0000000000000000, 0x0000000008000000, 0x0000000000000000,
        0x0000000000008400, 0x0000000000000010, 0x0100000000000000, 0x0004200000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000210000,
        0x0000000000000400, 0x0000000000000010, 0x0008000000000000, 0x0000000000000000,
        0x0000008000000000, 0x0000000000000000, 0x0000000000400000, 0x0000000000010000,
        0x0000000000000000, 0x0200000000000010, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000010000000, 0x0000000000400000, 0x0000000000010000,
        0x8000000000000000, 0x0000000000000000, 0x0008000000000000, 0x0000000000000000,
        0x0000000400000000, 0x0000000010000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000020, 0x0000000000000000, 0x0000000000000000, 0x0000010000000000,
        0x0000000400000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000800,
        0x0000000000000000, 0x0000000000000000, 0x0000400000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000020000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
    ];

    #[derive(Clone, Copy)]
    pub struct Prng {
        mt: [u32; N],
        /// The index of the next word to return; `N` if the array must be
        /// regenerated.
        i: usize,
    }

    /// Replaces the words of the array with the following ones in the
    /// sequence.
    fn regenerate(mt: &mut [u32; N]) {
        let twist = |u: u32, l: u32, m: u32| {
            let y = (u & UPPER_MASK) | (l & LOWER_MASK);
            m ^ (y >> 1) ^ ((y & 1).wrapping_neg() & MATRIX_A)
        };
        for k in 0..N - M {
            mt[k] = twist(mt[k], mt[k + 1], mt[k + M]);
        }
        for k in N - M..N - 1 {
            mt[k] = twist(mt[k], mt[k + 1], mt[k + M - N]);
        }
        mt[N - 1] = twist(mt[N - 1], mt[0], mt[M - 1]);
    }

    impl Prng {
        pub const NAME: &str = "MT19937";
        pub const BITS: usize = 32;

        pub fn new(seed: u64) -> Self {
            let mut x = seed;
            Self {
                mt: std::array::from_fn(|_| (splitmix64(&mut x) >> 32) as u32),
                i: N,
            }
        }

        #[inline(always)]
        pub fn next_u64(&mut self) -> u64 {
            if self.i >= N {
                regenerate(&mut self.mt);
                self.i = 0;
            }
            let mut y = self.mt[self.i];
            self.i += 1;
            y ^= y >> 11;
            y ^= (y << 7) & 0x9d2c5680;
            y ^= (y << 15) & 0xefc60000;
            y ^= y >> 18;
            (y as u64) << 32
        }

        f2_try_skip!();
    }

    impl crate::f2::LinearGenerator<{ N / 2 }> for Prng {
        fn step(&mut self) {
            self.next_u64();
        }
        // The next N words of the sequence, two per element.
        fn to_vector(&self) -> [u64; N / 2] {
            let mut next = self.mt;
            if self.i > 0 {
                regenerate(&mut next);
            }
            let w = |j: usize| {
                if self.i + j < N {
                    self.mt[self.i + j]
                } else {
                    next[self.i + j - N]
                }
            };
            std::array::from_fn(|j| w(2 * j) as u64 | (w(2 * j + 1) as u64) << 32)
        }
        fn set_vector(&mut self, v: [u64; N / 2]) {
            for (j, x) in v.into_iter().enumerate() {
                self.mt[2 * j] = x as u32;
                self.mt[2 * j + 1] = (x >> 32) as u32;
            }
            self.i = 0;
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        // The reference implementation, initialized by init_genrand(5489)
        // (the default seed), returns 3499211612 first and 4123659995 as
        // 10000th output (the value required by the C++ standard for
        // std::mt19937).
        #[test]
        fn test_reference() {
            let mut mt = [5489u32; N];
            for i in 1..N {
                mt[i] = 1812433253u32
                    .wrapping_mul(mt[i - 1] ^ (mt[i - 1] >> 30))
                    .wrapping_add(i as u32);
            }
            let mut g = Prng { mt, i: N };
            assert_eq!(g.next_u64() >> 32, 3499211612);
            for _ in 1..9999 {
                g.next_u64();
            }
            assert_eq!(g.next_u64() >> 32, 4123659995);
        }
    }
}

#[cfg(feature = "mt19937")]
pub use mt19937::Prng;

#[cfg(all(test, feature = "mt19937"))]
use mt19937::CHARPOLY;

/// The SIMD-oriented Fast Mersenne Twister (SFMT) with 19937 bits of state, a
/// port of the generic (non-SIMD) code of `SFMT.c` by Mutsuo Saito and Makoto
/// Matsumoto, returning 64-bit values as `sfmt_genrand_uint64()`.
#[cfg(feature = "sfmt19937")]
mod sfmt19937 {
    use super::splitmix64;

    // The parameters of SFMT-params19937.h.
    const N: usize = 156;
    const POS1: usize = 122;
    const SL1: u32 = 18;
    const SL2: u32 = 1;
    const SR1: u32 = 11;
    const SR2: u32 = 1;
    const MSK: [u32; 4] = [0xdfffffef, 0xddfecb7f, 0xbffaffff, 0xbffffff6];
    const PARITY: [u32; 4] = [0x00000001, 0x00000000, 0x00000000, 0x13c9e684];

    // The characteristic polynomial of the transition map on 156 128-bit
    // words, without the leading coefficient.
    #[rustfmt::skip]
    pub(super) const CHARPOLY: [u64; 2 * N] = [
        0x0000000000000001, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000020000, 0x0000000000000000, 0x0000280000000000, 0x0000100000010000,
        0x0000000000000000, 0x00000000000000a0, 0x0000000000000140, 0x0000000a00000000,
        0x1100001400000000, 0x8200000000000000, 0x0000000000200000, 0x0000000000540001,
        0x0000800008280000, 0x0011400000000000, 0x00a0800000000000, 0x0000400000000000,
        0x0000000400000040, 0x00000088000a0800, 0x1000004400000500, 0x000080a000000200,
        0x4400001400000020, 0x0000002004400010, 0x5008800108000808, 0x0010000105500001,
        0x000200a0a2200002, 0x0441000011008200, 0x0802804020810400, 0x011100880008000a,
        0x000a040001054100, 0x2020000082a20805, 0x400040500015140a, 0x8820000810002804,
        0x0009415808808808, 0xa00102a400010500, 0xc2800000a0a0202a, 0x000404440002d011,
        0x0a8200020ad444a2, 0x0111040115028080, 0x028c170826050105, 0x1512280000828020,
        0x4220144044141282, 0x8008200044880c20, 0x04d0501029419208, 0x2260a9a201342400,
        0x42b011808200a0a0, 0x6043668044470047, 0x02800a0280028a42, 0x101448840d038108,
        0x8a0249858605200f, 0x0614be1144282080, 0x982846067854480e, 0x005aa00480a040d4,
        0x24915400e1456171, 0x41a8002841a1a172, 0x8236063b70e58020, 0x045a0e4302c2c14e,
        0x1e1281a00c8a020e, 0x363f214518ba8948, 0x0271812261458186, 0x1cc00e4401222930,
        0x2854d800a7263278, 0x68b8aa40a02c9855, 0x0373450904111bdc, 0x2600f3a0602350a0,
        0x644f4c31b0bb80a5, 0xd20a1d4608ef2560, 0x8a0b4b9211312406, 0x82ae7517009f9982,
        0x601756539a200074, 0x640a70122436867c, 0x8020484dc862a715, 0xe2e2e81b48d8b424,
        0x72b122c30548ac38, 0xd2a3870790381012, 0xcb56e0eece1c70b7, 0x3502990347470682,
        0xa8601f8b5c7411e1, 0x3775a12a833a30b5, 0x052143a0016f9a44, 0xce3b6a1212780c22,
        0xc05c5070c11b954e, 0xa6b0b13223bc8d00, 0x26110291d7d998c0, 0x2097e7a161246d50,
        0x8d4d25c4574d475c, 0x152e14187c8b1e6b, 0xda950b3fcb88e537, 0x835944751836d521,
        0x2636a40253002240, 0xffef9c51964912a5, 0x7d4964adc523308c, 0x69f98f32aa726ab9,
        0x47130b37425091dd, 0x401ab0ff24e21061, 0x9453c512e050cd4b, 0x1ac684510d88fa5f,
        0xa16ca218b2933017, 0x5424cd6cea03afba, 0x3df8a93b3b286b75, 0x32873ba3471bc681,
        0x5b798315ecd48145, 0xb45a9468ba2e3b9e, 0xd571d4457ecae4b2, 0xc9d63e3bd3bbfa43,
        0x192bea7fa9441ce2, 0x79b6d1bcc6cfa705, 0xd63fc57efa82ca0b, 0xc839574ca64d7f35,
        0xecf5868d70ee9058, 0x29f4a75568cf95db, 0x67a6382493eac127, 0xd196437f9f4a71cb,
        0x1b3022c27d461c7f, 0xa6a567ce4085d0bc, 0xae311af7b7278a1e, 0xa48c600294a94bfc,
        0xd624ca7a2f95b256, 0x560241615d847f18, 0xc6371879a520d42c, 0xd08d5f07d17e3abd,
        0x3df9d3be7ad73124, 0xc33686612cb4cbfa, 0x2dbe79740e8090c0, 0x30a4a80f6d4c79ec,
        0x5519d7912ce7f435, 0xc764fa909d0b2688, 0x27c655cfecc233f7, 0xe85987a8af20a5f9,
        0xd411bc7314c8d5dc, 0x93899b016b45a3f0, 0x61f5d113c20b0df0, 0xb25da61e4a096903,
        0x0dbe028d6d3567af, 0x9fa2ffe90c694a8b, 0xddbc8fc13fbb001b, 0xd4f0394b007675b1,
        0x82a77db81439b4c5, 0xe3926b17cba15b02, 0x8c9459c774f90065, 0xc96951bd97a7280d,
        0xd05abe912bca7f94, 0x60711d1a815f1c57, 0x042d25ce0d6cfd66, 0xe26807fc63178c4f,
        0x7ce8a197b575c993, 0x40b7cd97348c4e6e, 0x4121abca0b44faf6, 0xe52018057e436e7c,
        0xeee29d71348ff820, 0x5897af73be049411, 0x0a6fdc8a2abfe601, 0x9927489f06e9acb9,
        0x212a9e204d2b3555, 0x726f34b152c7e23b, 0xba18032b9081e787, 0x1e6fd7621f8d4fce,
        0xddc1ca0a680b74f2, 0x0b73fbbb3926fb78, 0x99f11bf5fbcb7c8c, 0xfa95b50d32e55b88,
        0x898481c3f32feb9f, 0x0c5530801a0da142, 0xe8d7a917f97df770, 0x4875f816a8423596,
        0xdbb428b030a50aa9, 0x0e3950a4612c5231, 0xe3e8182323c04d1d, 0x391f65dd70a31feb,
        0xd0037d2ea87036c2, 0x585cb2a68d024115, 0x3ca80652b82e08da, 0x1222a69b8994a108,
        0x4de6d9cdceae67bc, 0xddca8edabd55bf58, 0xf6a0757e4667e48e, 0x9b32d9f9b71a27e7,
        0x40f2769f8f20f8f8, 0x45043e807c88737f, 0xb8ee0dd038f6f4af, 0x1484c5e77d62c435,
        0x8dd2569dfa4d9131, 0x5f523ec999db3861, 0x3418fa6737e8b00d, 0x269f5801674ff9a5,
        0x0cd977b54925f868, 0x0efe2aca2f5aac13, 0x56317da6a2f6b8c4, 0xe534d38250fa24dd,
        0xdfa8dc9afeb39524, 0xf68b95bdbfe9f66f, 0xcd69cc6772132bd7, 0xb5b4dfded98e8544,
        0x0387409dcb87d8d7, 0x8f0023832ffcb147, 0x2765011aafc4140f, 0x83081b652eca2bdd,
        0x4d14a10e4b5b0ac3, 0x7c88af6e819ec2c9, 0x0e191e6f25748090, 0xd6495ebd110a22f4,
        0xdbf1f3cefb3cbcdf, 0x9448bef759c292ca, 0xa5634a3ae4d4acfb, 0x7164a8c8c26ad6a4,
        0x965e5a7cfb55c640, 0xdcf519a0992e424e, 0x8f610efdff342da1, 0xf9242248af2415d8,
        0x10c4b695164603b8, 0x1e87d6082fa1757b, 0x7a57a7a99015387c, 0x286a730fd18197c4,
        0x337303598db3d5d7, 0xfec20b20ffa6cb03, 0x420ebf29112f2932, 0x854a5d8b53939260,
        0xcb1a14d9f27695a2, 0x70d1a3a726ac668e, 0xf1b6da4284c007a7, 0x72a04fdc5cb3134e,
        0x2a3d847fe51d6b08, 0x3b3b804a91cea167, 0xc59263aa363cac3b, 0x034e799408af0885,
        0x006262ed52a6fa26, 0xe0acc024778a11e8, 0xcd4d4ab18447afca, 0x576f160423a6c70c,
        0x10631e8624500040, 0x02221f668cc007fe, 0x4b061c0105120745, 0x2b15ed7d4b520260,
        0x20410d99d63883d1, 0xe3375e48c3b54b20, 0xcc86a05034ecdea6, 0xced1542ae91014a1,
        0x622980024f61246e, 0x08b013659c68f806, 0xf5909002f128b242, 0x67d3234a7a8458be,
        0x201ac293eeaa9176, 0x0cb848026d5fa140, 0x5c02883711114816, 0x1c518a7c4631ec3a,
        0x164ab085407e6130, 0x00609822b1288189, 0x420e03588aad0882, 0xa0558040a144a900,
        0x0054b1a8b0022848, 0x0a974810486c5464, 0x20406990422a4880, 0x04201d5a0c864f08,
        0x00a14580208b518b, 0x2020d0b080740015, 0xc000b3323000a400, 0x13011049400a9948,
        0x8348220c6a884c49, 0x91500a5781080941, 0x16a001b492002140, 0x00a480923051a804,
        0x1b11001460854081, 0x010442001c20810a, 0x001a4d8101a30803, 0x4552001182b32021,
        0x900000c8b61000a0, 0x4831008010402074, 0xa9d1000a00180808, 0x2040020c42038108,
        0x80400040a0a03122, 0x448808048a111020, 0x0e8a1110001440a0, 0x0889100200080804,
        0x2201120805400101, 0x2000000040888030, 0x0450880048841500, 0x0408801100800028,
        0x00a8414002010808, 0x2220010280560201, 0x000000020000a804, 0x20050080000a0050,
        0x01000a0000000000, 0x1100800400000008, 0x0022000000004020, 0x0000000000100080,
        0x0000000000000004, 0x0800000000000000, 0x0010000040000000, 0x0000200000000002,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
    ];

    /// Returns a 128-bit word made of four 32-bit lanes, the first lane being
    /// the least significant one (as in the `w128_t` type on little-endian
    /// architectures).
    const fn lanes(u: [u32; 4]) -> u128 {
        u[0] as u128 | (u[1] as u128) << 32 | (u[2] as u128) << 64 | (u[3] as u128) << 96
    }

    /// `MSK`, without the bits that a lane shifted right by `SR1` receives
    /// from the next lane.
    const MSK_SR1: u128 = lanes(MSK) & lanes([u32::MAX >> SR1; 4]);
    /// Clears the bits that a lane shifted left by `SL1` receives from the
    /// previous lane.
    const MSK_SL1: u128 = lanes([u32::MAX << SL1; 4]);

    /// `do_recursion()` from `SFMT-common.h`, with lane-wise shifts
    /// implemented by masking 128-bit shifts.
    #[inline(always)]
    fn do_recursion(a: u128, b: u128, c: u128, d: u128) -> u128 {
        a ^ (a << (8 * SL2)) ^ ((b >> SR1) & MSK_SR1) ^ (c >> (8 * SR2)) ^ ((d << SL1) & MSK_SL1)
    }

    /// Replaces the words of the array with the following ones in the
    /// sequence (`sfmt_gen_rand_all()` in `SFMT.c`).
    fn regenerate(s: &mut [u128; N]) {
        let (mut r1, mut r2) = (s[N - 2], s[N - 1]);
        for i in 0..N - POS1 {
            s[i] = do_recursion(s[i], s[i + POS1], r1, r2);
            (r1, r2) = (r2, s[i]);
        }
        for i in N - POS1..N {
            s[i] = do_recursion(s[i], s[i + POS1 - N], r1, r2);
            (r1, r2) = (r2, s[i]);
        }
    }

    #[derive(Clone, Copy)]
    pub struct Prng {
        s: [u128; N],
        /// The index of the next 64-bit half to return; 2`N` if the array
        /// must be regenerated.
        i: usize,
    }

    impl Prng {
        pub const NAME: &str = "SFMT19937";
        pub const BITS: usize = 64;

        pub fn new(seed: u64) -> Self {
            let mut x = seed;
            let mut s = std::array::from_fn(|_| {
                let lo = splitmix64(&mut x);
                lo as u128 | (splitmix64(&mut x) as u128) << 64
            });
            period_certification(&mut s);
            Self { s, i: 2 * N }
        }

        #[inline(always)]
        pub fn next_u64(&mut self) -> u64 {
            if self.i >= 2 * N {
                regenerate(&mut self.s);
                self.i = 0;
            }
            let r = (self.s[self.i / 2] >> (64 * (self.i % 2))) as u64;
            self.i += 1;
            r
        }

        pub fn try_skip(&mut self, mut n: u64) -> Result<(), ()> {
            // A step of the transition map yields two outputs, so we jump by
            // pairs of outputs from an even position.
            if n > 0 && self.i % 2 != 0 {
                self.next_u64();
                n -= 1;
            }
            crate::f2::jump(self, n / 2, &CHARPOLY);
            if n % 2 != 0 {
                self.next_u64();
            }
            Ok(())
        }
    }

    /// Ensures that the period is a multiple of 2¹⁹⁹³⁷ − 1
    /// (`period_certification()` in `SFMT.c`).
    fn period_certification(s: &mut [u128; N]) {
        let parity = lanes(PARITY);
        if (s[0] & parity).count_ones() % 2 == 0 {
            s[0] ^= 1 << parity.trailing_zeros();
        }
    }

    // A step yields a 128-bit word, that is, two outputs: the methods of this
    // trait must be called only when the next output is the lower half of a
    // word, as try_skip does.
    impl crate::f2::LinearGenerator<{ 2 * N }> for Prng {
        fn step(&mut self) {
            self.next_u64();
            self.next_u64();
        }
        // The next N words of the sequence, two elements per word, lower half
        // first.
        fn to_vector(&self) -> [u64; 2 * N] {
            debug_assert!(self.i % 2 == 0);
            let j = self.i / 2;
            let mut next = self.s;
            if j > 0 {
                regenerate(&mut next);
            }
            let w = |k: usize| {
                if j + k < N {
                    self.s[j + k]
                } else {
                    next[j + k - N]
                }
            };
            std::array::from_fn(|k| (w(k / 2) >> (64 * (k % 2))) as u64)
        }
        fn set_vector(&mut self, v: [u64; 2 * N]) {
            self.s = std::array::from_fn(|k| v[2 * k] as u128 | (v[2 * k + 1] as u128) << 64);
            self.i = 0;
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        // Outputs of the reference implementation (sfmt_genrand_uint64()),
        // initialized by sfmt_init_gen_rand(1234).
        #[test]
        fn test_reference() {
            let mut u = [1234u32; 4 * N];
            for i in 1..4 * N {
                u[i] = 1812433253u32
                    .wrapping_mul(u[i - 1] ^ (u[i - 1] >> 30))
                    .wrapping_add(i as u32);
            }
            let mut s = std::array::from_fn(|k| lanes(u[4 * k..4 * k + 4].try_into().unwrap()));
            period_certification(&mut s);
            let mut g = Prng { s, i: 2 * N };
            let out: Vec<u64> = (0..10000).map(|_| g.next_u64()).collect();
            assert_eq!(out[0], 0x5d47f5d7cd0d0032);
            assert_eq!(out[1], 0xaea87b245a0afbf6);
            assert_eq!(out[2], 0xe24675a556927984);
            assert_eq!(out[311], 0x993b14e5485e3509);
            assert_eq!(out[312], 0x6162c162e870c52d);
            assert_eq!(out[9999], 0x41e7bc118157ec5b);
        }
    }
}

#[cfg(feature = "sfmt19937")]
pub use sfmt19937::Prng;

#[cfg(all(test, feature = "sfmt19937"))]
use sfmt19937::CHARPOLY;

// Checks the characteristic polynomial of 𝐅₂-linear generators: the minimal
// polynomial of the sequence of the lowest bit of the state must be equal to it
// (if it has full degree, it is the characteristic polynomial). For MT19937,
// the minimal polynomial is the primitive factor of degree 19937, and the
// characteristic polynomial is x³¹ times it.
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
        feature = "well512a",
        feature = "mt19937",
        feature = "sfmt19937"
    )
))]
mod charpoly_tests {
    use super::*;
    use crate::f2::{LinearGenerator, minimal_polynomial};

    #[test]
    fn test_charpoly() {
        let n = CHARPOLY.len();
        // The power of x dividing the characteristic polynomial.
        let x_pow = if cfg!(feature = "mt19937") { 31 } else { 0 };
        let mut g = Prng::new(0);
        let bits: Vec<bool> = (0..128 * n)
            .map(|_| {
                let b = g.to_vector()[0] & 1 != 0;
                g.step();
                b
            })
            .collect();
        let m = minimal_polynomial(&bits);
        let degree = 64 * m.len() - 1 - m[m.len() - 1].leading_zeros() as usize;
        assert_eq!(
            degree,
            64 * n - x_pow,
            "the minimal polynomial has not the expected degree"
        );
        let mut c = vec![0u64; n + 1];
        for k in (0..=degree).filter(|&k| (m[k / 64] >> (k % 64)) & 1 != 0) {
            c[(k + x_pow) / 64] |= 1 << ((k + x_pow) % 64);
        }
        assert_eq!(c[..n], CHARPOLY, "computed: {:#018x?}", &c[..n]);
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
        // Skips start also after some outputs, as generators returning a
        // block of words at a time might be in the middle of a block.
        for pre in [0, 1, 1000] {
            // Cheap-to-step values (sequential stepping must stay fast).
            for &n in &[0u64, 1, 2, 7, 1000, 100_000] {
                let mut a = Prng::new(seed);
                for _ in 0..pre {
                    a.next_u64();
                }
                let mut b = a;
                if a.try_skip(n).is_err() {
                    return Ok(()); // generator has no jump-ahead; nothing to check
                }
                for _ in 0..n {
                    b.next_u64();
                }
                for k in 0..64 {
                    assert_eq!(
                        a.next_u64(),
                        b.next_u64(),
                        "skip({n}) after {pre} outputs disagreed with stepping at output {k}"
                    );
                }
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
