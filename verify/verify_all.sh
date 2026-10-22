#!/bin/sh
# Runs every check in this repository. From the repository root:  sh verify/verify_all.sh
# Needs a C compiler with pthreads, python3 and numpy. Builds into ./build (about 5 minutes on a 15-core laptop).
set -e
C=certificates
B=build
mkdir -p $B
ARCH=""
case "$(uname -m)" in arm64|aarch64) ARCH="-mcpu=native" ;; esac

echo "== SHA-256 of all files"
if command -v sha256sum >/dev/null; then sha256sum -c --quiet SHA256SUMS; else shasum -a 256 -c --quiet SHA256SUMS; fi

echo "== proof-based checks (code distances and support intersections)"
python3 verify/verify_ers.py all

echo "== brute force over all pairs, dimensions 33-35, 37, 39"
cc -O3 -pthread $ARCH -o $B/verify_ers_full verify/verify_ers_full.c
$B/verify_ers_full 33 pad32 $C/gen32_17.txt $C/d33/supports_1852.txt
$B/verify_ers_full 34 pad32 $C/gen32_17.txt $C/d34/supports_2013.txt
$B/verify_ers_full 35 pad32 $C/gen32_17.txt $C/d35/supports_2202.txt
$B/verify_ers_full 37 pad32 $C/gen32_17.txt $C/d37/supports_2860.txt
$B/verify_ers_full 39 coset $C/d39/dense_kernel.txt $C/d39/dense_representatives.txt $C/d39/supports_3424.txt

echo "== dimension 38: equatorial lines, independent C check"
cc -O3 -o $B/d38_check_equator verify/d38_check_equator.c
$B/d38_check_equator $C/d38/base-data/golay_basis.txt $C/d38/equator_lines_6855.txt

echo "== dimension 38: structural check of base configuration plus equator"
python3 verify/d38_verify.py --base-data $C/d38/base-data --lines $C/d38/equator_lines_6855.txt --output $B/d38_report.json > /dev/null
python3 -c "import json; r = json.load(open('$B/d38_report.json')); assert r['passed'] and r['points'] == 605322; print('passed:', r['points'], 'points')"

echo "== dimension 38: brute force over all pairs of the 605322 points"
python3 verify/d38_export.py $C/d38/base-data $C/d38/equator_lines_6855.txt $B
cc -O3 -pthread $ARCH -o $B/d38_full verify/d38_full.c
$B/d38_full $B/d38_rows_i8.bin $B/d38_axes_i64.bin

echo "ALL CHECKS PASSED"
