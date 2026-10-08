/*
 * SPDX-FileCopyrightText: 2004-2016 David Blackman
 * SPDX-FileCopyrightText: 2017-2026 David Blackman and Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

//! The ways in which the output of the generator is split into words.

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
