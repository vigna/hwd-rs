/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

/* Generators identical to those of the Rust crate, for compiling the C
   implementation of the test (hwd.c) with check.sh. The state of the
   𝐅₂-linear generators is filled with SplitMix64 from SEED, as in the crate. */
#ifndef SEED
#define SEED 0
#endif
#define ROTL64 static inline uint64_t rotl(const uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

static uint64_t sm_x = SEED;
static uint64_t splitmix64(void) {
	uint64_t z = (sm_x += 0x9e3779b97f4a7c15);
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
	return z ^ (z >> 31);
}

#if defined(INCR)
static uint64_t x = SEED;
static uint64_t inline next() { return ++x; }
static void init(void) {}
#elif defined(XORSHIFT128) || defined(XORSHIFT128PLUS)
#define A 23
#define B 18
#define C 5
static uint64_t s[2];
static uint64_t inline next() {
#include "xorshift128-next.c"
#ifdef XORSHIFT128PLUS
	return result_plus;
#else
	return s0;
#endif
}
static void init(void) { s[0] = splitmix64(); s[1] = splitmix64(); }
#elif defined(XORSHIFT1024) || defined(XORSHIFT1024PLUS)
#define A 31
#define B 11
#define C 30
int p;
static uint64_t s[16];
static uint64_t inline next() {
#include "xorshift1024-next.c"
#ifdef XORSHIFT1024PLUS
	return result_plus;
#else
	return s0;
#endif
}
static void init(void) { for(int i = 0; i < 16; i++) s[i] = splitmix64(); }
#elif defined(XOROSHIRO128) || defined(XOROSHIRO128PLUS)
#define A 24
#define B 16
#define C 37
ROTL64
static uint64_t s[2];
static uint64_t inline next() {
#include "xoroshiro128-next.c"
#ifdef XOROSHIRO128PLUS
	return result_plus;
#else
	return s0;
#endif
}
static void init(void) { s[0] = splitmix64(); s[1] = splitmix64(); }
#elif defined(XOROSHIRO1024) || defined(XOROSHIRO1024PLUS)
#define A 25
#define B 27
#define C 36
ROTL64
int p;
static uint64_t s[16];
static uint64_t inline next() {
#define XOROSHIRO_USE_P
#include "xoroshiro1024-next.c"
#ifdef XOROSHIRO1024PLUS
	return result_plus;
#else
	return s0;
#endif
}
static void init(void) { for(int i = 0; i < 16; i++) s[i] = splitmix64(); }
#elif defined(WELL512A)
#define W 32
#define R 16
#define M1 13
#define M2 9
#define M3 5
#define MAT0POS(t,v) (v^(v>>t))
#define MAT0NEG(t,v) (v^(v<<(-(t))))
#define MAT3NEG(t,v) (v<<(-(t)))
#define MAT4NEG(t,b,v) (v ^ ((v<<(-(t))) & b))
#define V0            STATE[state_i                   ]
#define VM1           STATE[(state_i+M1) & 0x0000000fU]
#define VM2           STATE[(state_i+M2) & 0x0000000fU]
#define VM3           STATE[(state_i+M3) & 0x0000000fU]
#define VRm1          STATE[(state_i+15) & 0x0000000fU]
#define VRm2          STATE[(state_i+14) & 0x0000000fU]
#define newV0         STATE[(state_i+15) & 0x0000000fU]
#define newV1         STATE[state_i                 ]
#define newVRm1       STATE[(state_i+14) & 0x0000000fU]
static unsigned int state_i = 0;
static uint32_t STATE[R];
static uint32_t z0, z1, z2;
static uint32_t well512a(void){
  z0    = VRm1;
  z1    = MAT0NEG (-16,V0)    ^ MAT0NEG (-15, VM1);
  z2    = MAT0POS (11, VM2)  ;
  newV1 = z1                  ^ z2;
  newV0 = MAT0NEG (-2,z0)     ^ MAT0NEG(-18,z1)    ^ MAT3NEG(-28,z2) ^ MAT4NEG(-5,0xda442d24U,newV1) ;
  state_i = (state_i + 15) & 0x0000000fU;
  return STATE[state_i];
}
/* hwd.c declares next() with the type of the output of the generator: with 64
   bits, the output is in the lower bits (in the crate, in the upper bits). */
#if HWD_PRNG_BITS == 64
static inline uint64_t next(void) { return well512a(); }
#else
static inline uint32_t next(void) { return well512a(); }
#endif
static void init(void) { for(int i = 0; i < 16; i++) STATE[i] = splitmix64() >> 32; }
#else
#error "No generator"
#endif

__attribute__((constructor)) static void init_prng(void) { init(); }
