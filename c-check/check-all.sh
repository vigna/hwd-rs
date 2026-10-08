#!/bin/bash
#
# SPDX-FileCopyrightText: 2026 Sebastiano Vigna
#
# SPDX-License-Identifier: Apache-2.0 OR MIT
#
# Runs check.sh on a set of configurations covering all word modes, the
# transitional variant, categories, small batches, early termination,
# overflows, and jump-ahead and pre-scan parallel generation.

set -e
C=$(cd "$(dirname "$0")" && pwd)/check.sh
for RUSTOPT in "" "-P"; do
	export RUSTOPT
	$C SPLITMIX splitmix 64 64 8 --progress 3e9
	$C SPLITMIX splitmix 64 64 8 -t --progress 2e9
	$C SPLITMIX splitmix 32 64 7 -t --progress 2e9
	$C SPLITMIX splitmix 32 64 12 --progress 2e9
	$C SPLITMIX splitmix 128 64 8 -t --progress 2e9
	$C LCG32 lcg_32_32_0xec65035 32 32 8 --progress 1e9
	NUMCATS=3 $C SPLITMIX splitmix 64 64 10 --progress 1e9
	MAXB=1000 $C SPLITMIX splitmix 64 64 6 -t --progress 2e8
	$C XORSHIFT128 xorshift128 64 64 8 --progress --low-pv=1e-20 1e10
	$C XORSHIFT128PLUS xorshift128plus 64 64 8 -t --progress --low-pv=1e-20 1e10
	$C XORSHIFT1024 xorshift1024 64 64 16 --progress --low-pv=1e-20 2e9
	$C XOROSHIRO128 xoroshiro128 64 64 8 --progress 3e9
	$C XOROSHIRO1024PLUS xoroshiro1024plus 64 64 12 -t --progress 3e9
	$C WELL512A well512a 32 32 16 --progress 2e9
	$C ROMUTRIO romutrio 32 64 10 -t --progress 2e9
	$C INCR incr 64 64 8 --progress 1e9
	$C INCR incr 128 64 4 -t --progress 1e9
done
