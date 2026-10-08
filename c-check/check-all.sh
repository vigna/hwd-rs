#!/bin/bash
#
# SPDX-FileCopyrightText: 2026 Sebastiano Vigna
#
# SPDX-License-Identifier: Apache-2.0 OR MIT
#
# Runs check.sh on a set of configurations covering all generators, all word
# modes, the transitional variant, categories, small batches, early
# termination, and overflows, both sequentially and in parallel.

set -e
C=$(cd "$(dirname "$0")" && pwd)/check.sh
for RUSTOPT in "" "-P"; do
	export RUSTOPT
	$C XOROSHIRO128PLUS xoroshiro128plus 64 64 8 --progress 3e9
	$C XOROSHIRO128PLUS xoroshiro128plus 64 64 8 -t --progress 2e9
	$C XOROSHIRO128PLUS xoroshiro128plus 32 64 7 -t --progress 2e9
	$C XOROSHIRO128PLUS xoroshiro128plus 32 64 12 --progress 2e9
	$C XOROSHIRO128PLUS xoroshiro128plus 128 64 8 -t --progress 2e9
	NUMCATS=3 $C XOROSHIRO128PLUS xoroshiro128plus 64 64 10 --progress 1e9
	MAXB=1000 $C XOROSHIRO128PLUS xoroshiro128plus 64 64 6 -t --progress 2e8
	$C XORSHIFT128 xorshift128 64 64 8 --progress --low-pv=1e-20 1e10
	$C XORSHIFT128PLUS xorshift128plus 64 64 8 -t --progress --low-pv=1e-20 1e10
	$C XORSHIFT1024 xorshift1024 64 64 16 --progress --low-pv=1e-20 2e9
	$C XORSHIFT1024PLUS xorshift1024plus 32 64 10 -t --progress 2e9
	$C XOROSHIRO128 xoroshiro128 64 64 8 --progress 3e9
	$C XOROSHIRO1024 xoroshiro1024 128 64 4 -t --progress 1e9
	$C XOROSHIRO1024PLUS xoroshiro1024plus 64 64 12 -t --progress 3e9
	$C WELL512A well512a 32 32 8 --progress 1e9
	$C WELL512A well512a 32 32 16 --progress 2e9
	# A 32-bit generator examined as 64-bit or 128-bit words overflows the
	# small counters. C returns the output in the lower bits, the crate in the
	# upper bits, so these runs must not use -t.
	$C WELL512A well512a 64 64 8 --progress 1e9
	$C WELL512A well512a 128 64 4 --progress 1e9
	# Overflows of counters and of sums with transitions.
	$C INCR incr 128 64 4 -t --progress 1e9
done
