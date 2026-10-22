#!/usr/bin/env python3
"""Write all 605,322 points of the 38-dimensional configuration as integer rows, for the all-pairs check d38_full.c.

Coordinates follow the Qiushi Engine export (constructions/d38/export_vectors.py there): the isometry
(x, y) -> (x, y/sqrt2, y/sqrt2) from R^24 + R^14 into R^24 + R^14 + R^14 makes every point an integer row of
length 52 (coordinates 38..51 repeat 24..37, so all rows lie in a 38-dimensional subspace).
  caps     6 * 98280 = 589680 rows (+-u, t, t): u a minimal Leech line, t one of the three 14-dim directions of
           its class (squared norm 32 + 8 + 8 = 48)
  equator  2 * 6855  =  13710 rows (+-v, 0, 0): v an added norm-48 Leech vector
  axes     1932 rows (0, Rt, Rt) with R = D * (rational rotation), integer entries up to about 1.6e17
The cap and equator rows go to d38_rows_i8.bin (int8, 64 bytes per row, zero padded); the axis rows go to
d38_axes_i64.bin (little-endian int64, 52 per row). Correctness of the output is not assumed anywhere:
d38_full.c checks every pair of the written rows directly.

usage: python3 verify/d38_export.py certificates/d38/base-data certificates/d38/equator_lines_6855.txt OUTDIR"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from d38_verify import shell  # noqa: E402  (minimal Leech lines in the shipped basis, canonical sorted order)

data, added_path, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
_, lines = shell(data)
labels = np.loadtxt(data / 'class_labels.txt', dtype=np.int64)
triples = np.loadtxt(data / 'triples.txt', dtype=np.int64)
directions = np.loadtxt(data / 'dim14_directions.txt', delimiter=',', dtype=np.int64)
keep = np.loadtxt(data / 'axis_keep.txt', dtype=np.int64)
rotation = np.array([[int(v) for v in row.split()] for row in (data / 'axis_rotation.txt').read_text().splitlines()
                     if row.strip()], dtype=object)
added = np.loadtxt(added_path, dtype=np.int64)

rows = []
for label in range(3):
    tails = directions[triples[labels, label]]
    for sign in (1, -1):
        block = np.zeros((len(lines), 64), dtype=np.int8)
        block[:, :24] = sign * lines
        block[:, 24:38] = tails
        block[:, 38:52] = tails
        rows.append(block)
for sign in (1, -1):
    block = np.zeros((len(added), 64), dtype=np.int8)
    block[:, :24] = sign * added
    rows.append(block)
rows = np.concatenate(rows)

rotated = directions[keep].astype(object) @ rotation.T
assert all(abs(int(v)) < 2**63 for v in rotated.flat), 'axis coordinate exceeds int64'
axes = np.zeros((len(keep), 52), dtype='<i8')
axes[:, 24:38] = rotated.astype(np.int64)
axes[:, 38:52] = rotated.astype(np.int64)

out.mkdir(parents=True, exist_ok=True)
rows.tofile(out / 'd38_rows_i8.bin')
axes.tofile(out / 'd38_axes_i64.bin')
print(f'{len(rows)} cap/equator rows and {len(axes)} axis rows written to {out} '
      f'(total {len(rows) + len(axes)} points)')
