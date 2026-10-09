#!/bin/bash
#
# SPDX-FileCopyrightText: 2026 Sebastiano Vigna
#
# SPDX-License-Identifier: Apache-2.0 OR MIT
#
# Checks that the crate and the C implementation of the test give the same
# results.
#
# Usage: check.sh C_GENERATOR FEATURE W PRNG_BITS K [ARGS...]
#
# C_GENERATOR is the macro selecting the generator in prngs_hwd.c, FEATURE the
# corresponding Cargo feature, W, PRNG_BITS, and K the word size, the number of
# output bits of the generator, and the signature length, and ARGS are passed to
# both implementations (e.g., -t --progress 1e9). The environment variables
# NUMCATS and MAXB set the number of categories and the maximum batch size,
# and RUSTOPT contains additional options for the crate (e.g., -P).
#
# HWD_C_DIR must point to the directory containing hwd.c and the
# *-next.c files of the C implementation (by default, ../hwd/c).
#
# The copy of hwd.c is patched to print p-values with 17 significant digits,
# rather than 3, so that they identify a double. Outputs are normalized by
# normalize.py, which removes timing information and writes p-values in a
# canonical form, so they compare equal only if they are the same double.

set -e
D=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$D/.." && pwd)
HWD_C_DIR=${HWD_C_DIR:-$ROOT/../hwd/c}
W=$ROOT/target/c-check
mkdir -p $W/src $W/bin $W/out
cp $HWD_C_DIR/*-next.c $D/prngs_hwd.c $W/src/
sed -e 's/p-value = %\.3g/p-value = %.17g/' -e 's/"p = %\.3g/"p = %.17g/' $HWD_C_DIR/hwd.c > $W/src/hwd.c
# The p-values of the categories, their minimum, and three final p-values.
if [ $(grep -c '%\.17g' $W/src/hwd.c) -ne 5 ]; then
	echo "Cannot patch the p-value formats of $HWD_C_DIR/hwd.c" >&2
	exit 1
fi

CDEF=$1; FEAT=$2; WB=$3; PB=$4; DIM=$5; shift 5
COPT=""; ROPT="-w $WB --prng-bits $PB -k $DIM $RUSTOPT"
if [ -n "$NUMCATS" ]; then COPT="$COPT -DHWD_NUMCATS=$NUMCATS"; ROPT="$ROPT -c $NUMCATS"; fi
if [ -n "$MAXB" ]; then COPT="$COPT -DMAX_BATCH_SIZE=$MAXB"; ROPT="$ROPT --max-batch-size $MAXB"; fi

CBIN=$W/bin/hwd-$CDEF-$WB-$PB-$DIM$(echo $COPT | tr -c 'a-zA-Z0-9\n' '_')
gcc -O3 -ffp-contract=off -w $W/src/hwd.c -I$W/src -D$CDEF -DHWD_BITS=$WB -DHWD_PRNG_BITS=$PB -DHWD_DIM=$DIM $COPT -o $CBIN -lm
(cd $ROOT && cargo build -q --release --features $FEAT && cp target/release/hwd $W/bin/hwd-$FEAT)

TAG=$CDEF-$WB-$PB-$DIM-$(echo "$@ $COPT $RUSTOPT" | tr -c 'a-zA-Z0-9.\n' '_')
$CBIN "$@" > $W/out/$TAG.c.txt 2>/dev/null
$W/bin/hwd-$FEAT $ROPT "$@" > $W/out/$TAG.rust.txt 2>/dev/null
python3 $D/normalize.py c $W/out/$TAG.c.txt > $W/out/$TAG.c.norm
python3 $D/normalize.py rust $W/out/$TAG.rust.txt > $W/out/$TAG.rust.norm
if cmp -s $W/out/$TAG.c.norm $W/out/$TAG.rust.norm; then
	echo "SAME $TAG ($(grep -c 'p = ' $W/out/$TAG.c.norm) p-values)"
else
	echo "DIFF $TAG"
	diff $W/out/$TAG.c.norm $W/out/$TAG.rust.norm | head -20
	exit 1
fi
