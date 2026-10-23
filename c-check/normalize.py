# SPDX-FileCopyrightText: 2026 Sebastiano Vigna
#
# SPDX-License-Identifier: Apache-2.0 OR MIT

"""Normalizes HWD output for comparison: drops timing, keeps only the exact
byte counts (check.sh patches the C code to print them exactly; the Rust code
prints them in parentheses, after an approximation using binary prefixes), and
writes all p-values as the shortest representation of the double they denote,
so that they compare equal only if they are the same double (check.sh patches
the C code to print p-values with 17 significant digits)."""
import re, sys

BYTES = {
    False: r'processed (\d+) bytes',
    True: r'processed \S+ \S+ \((\d+) bytes\)',
}

def norm(lines, rust):
    out = []
    for line in lines:
        line = line.rstrip('\n')
        m = re.match(BYTES[rust], line)
        if m:
            out.append('processed %s bytes' % m.group(1))
            continue
        line = re.sub(r'(p-value = |p = )(\S+)$', lambda m: m.group(1) + repr(float(m.group(2))), line)
        out.append(line)
    return out

rust = sys.argv[1] == 'rust'
print('\n'.join(norm(open(sys.argv[2]).readlines(), rust)))
