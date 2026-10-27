# C26a: an explicit {0,±1} family for the Boolean Bohnenblust–Hille inequality, 2(1−2⁻ᵈ)^((d+1)/2d) at every degree

**Claim.** For every d ≥ 1 there is an explicit f_d : {±1}^(3d) → {0, ±1} of degree at most d, with exactly
(2^d − 1)·2^d Fourier coefficients, all equal to ±2^(−d). Hence

  BH^{≤d}_{±1} ≥ R_d := 2 (1 − 2^(−d))^((d+1)/(2d)).

For d ≥ 2 this is larger than the degree-d address function's 2^((d−1)/d). R_d < 2 and R_d → 2.

| d | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|
| R_d | 1.611854897735 | 1.829652855011 | 1.920932230356 | 1.962262141550 | 1.981711051583 | 1.991056425313 |
| address 2^((d−1)/d) | 1.414214 | 1.587401 | 1.681793 | 1.741101 | 1.781797 | 1.811447 |

This relates to constant 26a of [teorth/optimizationproblems](https://teorth.github.io/optimizationproblems/constants/26a.html),
C₂₆ₐ = sup_d BH^{≤d}_{±1}. **It does not change the listed lower bound 2**, since every R_d < 2. It improves the
per-degree lower bounds that the page cites, which come from the address functions. The page itself says incremental
lower-bound improvements are not the interesting direction. `WRITEUP.md` §5 proves that no {0,±1}-valued function
can reach 2.

The proof is short and complete in [`WRITEUP.md`](WRITEUP.md) §2–3: sums of codimension-d affine blocks with pairwise
disjoint live sets and pairwise disjoint Fourier cosets, indexed by the proper subsets of ℤ_d. The certificates are
the explicit functions for d ≤ 7, checked exactly.

## Contents

| path | what |
|---|---|
| `WRITEUP.md` | theorem, proof, certificates, checks, literature, what is not claimed, AI-use disclosure |
| `problems/bh_cube/cycle.py` | writes f_d for any d (format `bh_cube/v1`) and scores it exactly |
| `problems/bh_cube/core.py` | search-side exact scorer: integer Walsh–Hadamard sup norm, interval ℓ_p norm (numpy, mpmath) |
| `problems/bh_cube/FORMAT.md` | the certificate format and the exact conditions a verifier must check (frozen spec) |
| `problems/bh_cube/instances/cycle_d{5,6,7}_n*.json` | f_5, f_6, f_7 as written by `cycle.py` |
| `problems/bh_cube/instances/d{2,3,4}_n*.json` | the d ≤ 4 functions first found by search (3d−1 variables, same R_d) |
| `verify/bh_cube.py` | independent verifier, written from `FORMAT.md` alone by a separate agent; exact integers and rationals, uses numpy for the sup norm when n > 16 |
| `verify/test_bh_cube.py` | its tests: hand cases, red controls (tampered/forged/malformed, `python -O`), the six files |
| `requirements.txt` | numpy, mpmath |
| `SHA256SUMS` | hashes of every file above |

## Reproduce

Python 3.10 or later, with `venv`. Run from the repository root:

```sh
sha256sum -c SHA256SUMS
python3 -m venv venv && venv/bin/pip install -q -r requirements.txt
# 1. build f_d for d = 2..7 and score it exactly; d = 5, 6, 7 must be byte-identical to the committed files
mkdir -p out && venv/bin/python -m problems.bh_cube.cycle 2 3 4 5 6 7 --out out
for d in 5 6 7; do cmp out/cycle_d${d}_n$((3*d)).json problems/bh_cube/instances/cycle_d${d}_n$((3*d)).json; done
# 2. search-side exact scorer on the committed certificates: "... R in [lo, hi] addr=... OK"
venv/bin/python -m problems.bh_cube.core problems/bh_cube/instances/*.json
# 3. independent verifier on the committed and the rebuilt files: "VALID [lo, hi] ..." (numpy needed for n > 16)
for f in problems/bh_cube/instances/*.json out/*.json; do venv/bin/python -I verify/bh_cube.py "$f"; done
# 4. the verifier's own tests and red controls (6 tests, OK)
venv/bin/python -I verify/test_bh_cube.py
```

Everything takes well under a minute on one core.

## Credits and AI use

The construction, the proof and the code were found and written with Claude (Anthropic, Claude Opus 5.5), directed
by Dennis Koshta, who takes responsibility for them. The route was numerical search, then structure extraction, then
a brute-force search over sign rules, which produced the closed form. The independent verifier `verify/bh_cube.py`,
its tests, and a reading of the proof were done by a separate AI agent (gemini-4-argon), working from `FORMAT.md` and
the proof text only. **No human has checked the proof or the code**, and the proof has not been formalised.

Cite as: Koshta, Dennis, with Claude (Anthropic, Claude Opus 5.5), *An explicit family for the Boolean
Bohnenblust–Hille inequality*, 2026, https://github.com/DennisKoshta/c26a-bohnenblust-hille-cube-family.

## Licence

MIT; see `LICENSE`. Copyright (c) 2026 Dennis Koshta.
