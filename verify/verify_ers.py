#!/usr/bin/env python3
"""Proof-based check of the Edel-Rains-Sloane three-layer kissing configurations in dimensions 33-35, 37 and 39.

usage (from the repository root):  python3 verify/verify_ers.py all      or      python3 verify/verify_ers.py 33 39
("all" also checks the constant-weight code certificates/d38/supports_3088.txt, which is not used for a kissing bound)

The configuration in R^n consists of
  dense layer   (-1)^c on coordinates 0..m-1 (zero elsewhere) for every word c of a binary code D, squared norm m;
                n = 33..37: D = Cheng-Sloane [32,17,8] code, m = 32 (certificates/gen32_17.txt);
                n = 39:     D = union of 10 cosets of a [39,15] linear code, m = 39 (certificates/d39/dense_*.txt);
  middle layer  +-2 on a weight-8 support S with an even number of minus signs, squared norm 32, for every S in
                certificates/dNN/supports_N.txt (128 vectors per support);
  root layer    +-4 e_i +-4 e_j (i < j), squared norm 32, 2n(n-1) vectors.
Two directions x, y are compatible (angle >= 60 degrees) iff <x,y> <= |x||y|/2. By pair type:
  dense-dense    m - 2 dist(c,c')            needs dist >= m/4, i.e. 8 for m = 32 and 10 for m = 39   (checked here)
  dense-middle   |<x,y>| <= 2*8 = 16          16^2 <= 32m/4 holds for m >= 32
  dense-root     |<x,y>| <= 2*4 = 8           holds
  middle-middle  same support: two distinct even sign patterns differ in >= 2 places: 4*(8-2*2) = 16 = 32/2
                 different supports S, T: |<x,y>| <= 4|S & T| <= 16 iff |S & T| <= 4                  (checked here)
  middle-root    |<x,y>| <= 2*(2*4) = 16      holds
  root-root      distinct roots: <x,y> <= 16  holds
So validity reduces to the minimum distance of D and the pairwise intersections of the supports; both are checked
exhaustively below. verify_ers_full.c checks every pair of vectors directly instead."""
import os
import sys
from itertools import combinations

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
CERT = os.path.join(ROOT, 'certificates')
SUPPORTS = {33: 'd33/supports_1852.txt', 34: 'd34/supports_2013.txt', 35: 'd35/supports_2202.txt',
            37: 'd37/supports_2860.txt', 39: 'd39/supports_3424.txt'}
PREVIOUS = {33: 364736, 34: 384836, 35: 409676, 37: 498024, 39: 763668}   # Cohn's table, checked 2026-10-22
# A(38,8,8) >= 3088 is a constant-weight code record only: with the 180224-word dense layer of Yang et al. it gives
# 578300 < 605322, so the 38-dimensional kissing bound comes from a different construction (verify/d38_*).
CODE_ONLY = {38: 'd38/supports_3088.txt'}


def binary_rows(path, length):
    rows = [l.strip() for l in open(path) if l.strip() and not l.startswith('#')]
    assert all(len(r) == length and set(r) <= {'0', '1'} for r in rows), f'bad line format in {path}'
    return [int(r[::-1], 2) for r in rows]                     # bit i = coordinate i


def hex_rows(path, length):
    rows = [int(l, 16) for l in (l.strip() for l in open(path)) if l and l[0] not in '#$']
    assert all(r >> length == 0 for r in rows), f'row longer than {length} bits in {path}'
    return rows


def span(gens):
    """All 2^k words of the linear code spanned by gens; fails if the rows are dependent."""
    words = [0]
    for g in gens:
        words += [w ^ g for w in words]
    assert len(set(words)) == len(words), 'generator rows are linearly dependent'
    return words


def dense_layer(n):
    if n == 39:
        kernel = span(hex_rows(os.path.join(CERT, 'd39/dense_kernel.txt'), 39))
        reps = hex_rows(os.path.join(CERT, 'd39/dense_representatives.txt'), 39)
        m, k = 39, len(kernel).bit_length() - 1
        # dist(r_i + a, r_j + b) = wt(r_i ^ r_j ^ a ^ b): minimum over the kernel and over each coset difference
        d_kernel = min(bin(w).count('1') for w in kernel if w)
        d_cross = min(bin(r ^ s ^ w).count('1') for r, s in combinations(reps, 2) for w in kernel)
        d = min(d_kernel, d_cross)
        size = len(reps) * len(kernel)
        desc = f'{len(reps)} cosets of a [39,{k},{d_kernel}] code, {size} words, minimum distance {d}'
    else:
        code = span(binary_rows(os.path.join(CERT, 'gen32_17.txt'), 32))
        m, size = 32, len(code)
        d = min(bin(w).count('1') for w in code if w)
        desc = f'[32,{len(code).bit_length() - 1},{d}] code, {size} words (coordinates 0..31)'
    assert 4 * d >= m, f'dense layer minimum distance {d} too small for squared norm {m}'
    return size, m, desc


def supports(n):
    masks = binary_rows(os.path.join(CERT, SUPPORTS.get(n) or CODE_ONLY[n]), n)
    assert all(bin(s).count('1') == 8 for s in masks), 'support weight != 8'
    assert len(set(masks)) == len(masks), 'duplicate support'
    worst = max(bin(a & b).count('1') for a, b in combinations(masks, 2))
    assert worst <= 4, f'two supports meet in {worst} points'
    return len(masks), worst


def check(n):
    size, m, desc = dense_layer(n)
    N, worst = supports(n)
    tau = size + 128 * N + 2 * n * (n - 1)
    print(f'dim {n}: dense layer: {desc}')
    print(f'        {N} weight-8 supports, max pairwise intersection {worst}  =>  A({n},8,8) >= {N}')
    print(f'        tau_{n} >= {size} + 128*{N} + {2 * n * (n - 1)} = {tau}   (previous best {PREVIOUS[n]}, +{tau - PREVIOUS[n]})')
    return tau


def check_code_only(n):
    N, worst = supports(n)
    print(f'code {n}: {N} weight-8 supports of a {n}-set, max pairwise intersection {worst}  =>  A({n},8,8) >= {N}')


if __name__ == '__main__':
    dims = sorted(SUPPORTS | CODE_ONLY) if sys.argv[1:] in ([], ['all']) else [int(a) for a in sys.argv[1:]]
    for n in dims:
        check(n) if n in SUPPORTS else check_code_only(n)
    print('ALL CHECKS PASSED')
