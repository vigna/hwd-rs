# SPDX-FileCopyrightText: 2026 Sebastiano Vigna
#
# SPDX-License-Identifier: Apache-2.0 OR MIT

"""Normalizes HWD output for comparison: drops timing, and rounds the
full-precision numbers of the Rust output to the precision of the C output."""
import re, sys

def g3(x):
    return '%.3g' % x

def norm(lines, rust):
    out = []
    for line in lines:
        line = line.rstrip('\n')
        m = re.match(r'processed (\S+) bytes', line)
        if m:
            b = m.group(1)
            out.append('processed %s bytes' % (g3(float(b)) if rust else b))
            continue
        if rust:
            line = re.sub(r'(p-value = |p = )(\S+)$', lambda m: m.group(1) + g3(float(m.group(2))), line)
        out.append(line)
    return out

rust = sys.argv[1] == 'rust'
print('\n'.join(norm(open(sys.argv[2]).readlines(), rust)))
