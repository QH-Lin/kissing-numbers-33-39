# New lower bounds for kissing numbers in dimensions 33–39

Explicit, independently checkable kissing configurations that improve the best known lower bounds in **six dimensions: 33, 34, 35, 37, 38 and 39**, together with new constant-weight codes **A(n, 8, 8)** for n = 33, 34, 35, 37, 38 and 39.

Qihong Lin (linqihong001@qq.com), 2026-10-22

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23554808.svg)](https://doi.org/10.5281/zenodo.23554808) — Zenodo archive, all versions: [doi:10.5281/zenodo.23554808](https://doi.org/10.5281/zenodo.23554808)

Dimension 32 (τ₃₂ ≥ 352320, A(32, 8, 8) ≥ 1713) has its own repository: [QH-Lin/kissing-number-32](https://github.com/QH-Lin/kissing-number-32).

## Results

| n | previous best τₙ | source of previous best | **new τₙ ≥** | gain | construction |
|---:|---:|---|---:|---:|---|
| 33 | 364736 | Pliny the Liberator (A(33,8,8) ≥ 1809) with Edel–Rains–Sloane | **370240** | +5504 | ERS, A(33,8,8) ≥ 1852 |
| 34 | 384836 | Pliny the Liberator (A(34,8,8) ≥ 1965) with Edel–Rains–Sloane | **390980** | +6144 | ERS, A(34,8,8) ≥ 2013 |
| 35 | 409676 | Yang et al., arXiv:2609.35051 | **415308** | +5632 | ERS, A(35,8,8) ≥ 2202 |
| 37 | 498024 | Yang et al., arXiv:2609.35051 | **499816** | +1792 | ERS, A(37,8,8) ≥ 2860 |
| 38 | 591900 | Yang et al., arXiv:2609.35051 | **605322** | +13422 | Kravatskiy's 591612 points + 6855 equatorial lines |
| 39 | 763668 | Yang et al., arXiv:2609.35051 | **768916** | +5248 | ERS with a 327680-word dense layer, A(39,8,8) ≥ 3424 |

"Previous best" is the value in [Cohn's table of kissing numbers](https://cohn.mit.edu/kissing-numbers/), checked on 2026-10-22. Dimension 36 is not improved (our best is a tie with the listed 484760).

New constant-weight codes (binary, length n, weight 8, minimum distance 8, i.e. 8-subsets of an n-set pairwise meeting in at most 4 points), compared with [Brouwer's table](https://aeb.win.tue.nl/codes/Andw.html) on 2026-10-22:

| n | 33 | 34 | 35 | 37 | 38 | 39 |
|---|---:|---:|---:|---:|---:|---:|
| previous A(n,8,8) ≥ | 1809 | 1965 | 2157 | 2846 | 3077 | 3383 |
| **new A(n,8,8) ≥** | **1852** | **2013** | **2202** | **2860** | **3088** | **3424** |

The code for n = 38 is a code record only; it does not give the best kissing configuration in dimension 38.

## Constructions

### Dimensions 33, 34, 35, 37, 39: Edel–Rains–Sloane three layers

Edel, Rains and Sloane, *On kissing numbers in dimensions 32 to 128*, Electron. J. Combin. 5 (1998) R22. In ℝⁿ take

1. **dense layer** — the vectors ((−1)^{c_i}) on coordinates 0, …, m−1 (zero elsewhere) for every word c of a binary code D; squared norm m;
   - n = 33, 34, 35, 37: D is the Cheng–Sloane [32, 17, 8] code (2¹⁷ = 131072 words, m = 32, `certificates/gen32_17.txt`);
   - n = 39: D is the union of 10 cosets of a [39, 15, 12] linear code (327680 words, minimum distance 10, m = 39), the dense layer of Yang et al. (`certificates/d39/dense_*.txt`);
2. **middle layer** — for each weight-8 support S in `certificates/dNN/supports_N.txt`, the 128 vectors equal to ±2 on S with an even number of minus signs; squared norm 32;
3. **root layer** — the 2n(n−1) vectors ±4 eᵢ ± 4 eⱼ; squared norm 32.

Hence τₙ ≥ |D| + 128·N + 2n(n−1):

| n | dense | middle | roots | total |
|---:|---:|---:|---:|---:|
| 33 | 131072 | 128 × 1852 = 237056 | 2112 | **370240** |
| 34 | 131072 | 128 × 2013 = 257664 | 2244 | **390980** |
| 35 | 131072 | 128 × 2202 = 281856 | 2380 | **415308** |
| 37 | 131072 | 128 × 2860 = 366080 | 2664 | **499816** |
| 39 | 327680 | 128 × 3424 = 438272 | 2964 | **768916** |

Two directions x, y are at angle ≥ 60° iff ⟨x, y⟩ ≤ |x||y|/2. By pair type:

- dense–dense: m − 2·dist(c, c′) ≤ m/2 iff dist ≥ m/4, i.e. 8 for m = 32 and 10 for m = 39;
- dense–middle: |⟨x, y⟩| ≤ 2·8 = 16, and 16² ≤ 32m/4 for m ≥ 32; dense–root: at most 8;
- middle–middle on the same support: two distinct even sign patterns differ in at least 2 places, so ⟨x, y⟩ ≤ 4·(8 − 2·2) = 16;
- middle–middle on supports S ≠ T: |⟨x, y⟩| ≤ 4·|S ∩ T| ≤ 16 iff |S ∩ T| ≤ 4;
- middle–root: at most 2·(2·4) = 16; root–root: at most 16.

So validity reduces to the minimum distance of D and to |S ∩ T| ≤ 4 for the supports. The new ingredient in each dimension is the support family, a constant-weight code proving A(n, 8, 8) ≥ N.

### Dimension 38: equatorial lines on Kravatskiy's configuration

The base is the 591612-point configuration of A. Kravatskiy ([KissingNumbers](https://github.com/alexlegeartis/KissingNumbers)): ℝ³⁸ = ℝ²⁴ ⊕ ℝ¹⁴, with 6 × 98280 = 589680 "cap" points (a minimal Leech vector in the head and one of three 14-dimensional directions of its class in the tail) and 1932 "axis" points in the tail; the equator ℝ²⁴ ⊕ 0 is empty. Yang et al. (arXiv:2609.35051, §7) observed that a set of vectors v of squared norm 48 in the same Leech lattice (minimum squared norm 32) can be added to the equator as ±v whenever the added vectors satisfy |v·w| ≤ 24 pairwise: a cap point has head fraction 2/3 of its squared norm, so its cosine with a new point is (v·u)/48 for a minimal vector u, and v·u ≤ 24 holds automatically because |v − u|² ≥ 32. They added 144 lines (591900 points).

We found **6855** such lines (an independent set in the graph on the 8386560 lines of the norm-48 shell, two lines adjacent when |v·w| > 24), giving

591612 + 2 · 6855 = **605322** points.

`certificates/d38/equator_lines_6855.txt` lists the 6855 vectors (24 integers per line) in the coordinates of the Golay basis in `certificates/d38/base-data/`, which is Kravatskiy's base configuration as packaged by Yang et al. (see `certificates/d38/THIRD_PARTY.md`).

## Verification

Run everything with `sh verify/verify_all.sh` (about 5 minutes on a 15-core Apple M5 Pro laptop; needs a C compiler, python3 and numpy).

| program | what it checks | result |
|---|---|---|
| `verify/verify_ers.py` | enumerates the dense codes (minimum distance 8, resp. 10 for the 10-coset code in dimension 39); weight 8, distinctness and pairwise intersections ≤ 4 of every support family | all pass; τ values as above |
| `verify/verify_ers_full.c` | builds every vector of the configurations in dimensions 33, 34, 35, 37, 39 and tests every pair exactly | 0 violations; 6.85 × 10¹⁰ / 7.64 × 10¹⁰ / 8.62 × 10¹⁰ / 1.25 × 10¹¹ / 2.96 × 10¹¹ pairs (20 s / 19 s / 21 s / 29 s / 89 s) |
| `verify/d38_check_equator.c` | written from scratch: Leech membership and squared norm 48 of the 6855 vectors, pairwise \|v·w\| ≤ 24, \|v·u\| ≤ 24 against all 196560 minimal vectors | max 24 / 24, 0 violations |
| `verify/d38_verify.py` | Yang et al.'s structural verifier of the base configuration (Golay code, class partition, 14-dimensional directions, exact rational rotation), extended to any number of equatorial lines | passed, 605322 points |
| `verify/d38_export.py` + `verify/d38_full.c` | writes all 605322 points as integer vectors and tests every pair exactly (int32, and `__int128` for the axis points) | 0 violations; 1.83 × 10¹¹ pairs (2 min) |

Expected last lines: `VALID kissing configuration of 370240 points in R^33` (and likewise for the other dimensions), and `VALID kissing configuration of 605322 points in R^38`.

## Files

- `certificates/gen32_17.txt` — generator matrix of the Cheng–Sloane [32, 17, 8] code (Y. Cheng, N. J. A. Sloane, SIAM J. Discrete Math. 2 (1989) 28–37), the same file as in the dimension-32 repository.
- `certificates/dNN/supports_N.txt` — the N supports for dimension NN, one per line as an NN-character 0/1 string (leftmost character = coordinate 0).
- `certificates/d39/dense_kernel.txt`, `dense_representatives.txt` — the dimension-39 dense layer (hex rows, bit i = coordinate i), from Yang et al.
- `certificates/d38/equator_lines_6855.txt` — the new equatorial lines; `certificates/d38/base-data/` — the base configuration; `certificates/d38/supports_3088.txt` — the code A(38, 8, 8) ≥ 3088.
- `verify/` — the programs above; `SHA256SUMS` — hashes of all files.

## Method (summary)

For the support families: maximum independent set search in the graph on all 8-subsets of {0, …, n−1}, two subsets adjacent when they share at least 5 points (the graph is never stored). The search combines an iterated local search in the style of Andrade–Resende–Werneck with a ruin-and-recreate loop that deletes a large structured part of the current family and rebuilds it; start points were earlier codes and lengthened codes from smaller n. For dimension 38: randomized greedy and local search for an independent set in the conflict graph on the norm-48 Leech lines. The search code will be released with a preprint.

## Citation

Please cite the Zenodo archive: Qihong Lin, *New lower bounds for kissing numbers in dimensions 33 to 39: explicit certificates*, doi:[10.5281/zenodo.23554808](https://doi.org/10.5281/zenodo.23554808). See also `CITATION.cff`.

## License

Files written for this repository are licensed under [CC BY 4.0](LICENSE). Material from other sources (the dimension-38 base data, the dimension-39 dense layer, and the verifier adapted in `verify/d38_verify.py`) keeps its original license; see `certificates/d38/THIRD_PARTY.md`.
