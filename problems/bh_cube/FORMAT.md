# Certificate format for C26a lower bounds (`bh_cube/v1`)

This file is self-contained: a verifier can be written from it alone.

## The constant

Tao's table, constant 26a (github.com/teorth/optimizationproblems, `constants/26a.md`). For d ≥ 1 put
p_d = 2d/(d+1). BH(d) is the least C such that every f : {±1}^n → ℝ of Fourier degree at most d (any n) satisfies

    ( Σ_{S ⊆ [n], |S| ≤ d} |f̂(S)|^{p_d} )^{1/p_d}  ≤  C · max_{x ∈ {±1}^n} |f(x)|,

where f(x) = Σ_S f̂(S) x^S and x^S = Π_{i∈S} x_i. The table's constant is C26a = sup_{d ≥ 1} BH(d).
Listed lower bound at the time of writing (2026-10-22): **2**, the limit of the address functions, whose ratio at
degree d is 2^{(d−1)/d} < 2. So every single f certifies

    C26a ≥ BH(d) ≥ R(f, d) := ‖f̂‖_{p_d} / ‖f‖_∞      for any d ≥ deg f.

(Taking d larger than deg f is allowed but only lowers R, since ‖·‖_p decreases in p.)

## What a certificate describes

One function f on {±1}^n, n ≥ 1, given by **integer** Fourier coefficients a_S ≠ 0 on a list of distinct subsets
S ⊆ {0, …, n−1} (variables are 0-indexed). R is invariant under scaling f, so integers lose nothing: any rational f is
an integer f times a constant, and the LP vertices our searches produce are rational. With

    F(x) = Σ_S a_S x^S   (an integer for every x ∈ {±1}^n),
    M    = max_{x ∈ {±1}^n} |F(x)|   (an integer, computed by enumerating all 2^n points),
    Σ    = Σ_S |a_S|^{p_d},

the certified ratio is **R = Σ^{1/p_d} / M**.

Conditions the verifier must check:

1. `n` is an integer ≥ 1, `d` an integer ≥ 1.
2. Every term's subset is a strictly increasing list of distinct integers in [0, n−1] with length ≤ d (the empty list
   is the constant term). No subset appears twice.
3. Every coefficient is a nonzero integer, written as a decimal string (optional leading `-`, no leading zeros, no
   `+`).
4. M > 0 (follows from at least one nonzero coefficient, since the characters are linearly independent; check anyway).
5. The claimed bound holds: a rigorous lower enclosure of R is ≥ `claimed_lower_bound`.

R is irrational in general (p_d is not an integer), so condition 5 needs directed rounding or interval arithmetic,
e.g. |a|^{p} = exp(p · ln|a|) in interval arithmetic at ≥ 100 digits. M itself is exact integer arithmetic; a fast
Walsh–Hadamard transform of the coefficient vector (indexed by bitmask, bit i ↔ variable i) gives F at all points in
O(n 2^n). Point x ↔ bitmask m with x_i = −1 iff bit i of m is set, so x^S = (−1)^{popcount(m & mask(S))}.

## File

UTF-8 JSON, a single object, LF line endings, no duplicate keys.

```json
{
  "format": "bh_cube/v1",
  "constant": "C26a",
  "n": 5,
  "d": 2,
  "terms": [[[1], "1"], [[2], "-1"], [[0, 1], "-1"], [[1, 3], "1"]],
  "claimed_lower_bound": "1.6118",
  "claimed_sup": "4",
  "note": "free text"
}
```

(The `terms` above are abbreviated; real files list every nonzero coefficient.)

- `format`: exactly `"bh_cube/v1"`; `constant`: exactly `"C26a"`.
- `n`, `d`: JSON integers.
- `terms`: non-empty list of `[subset, coefficient]` pairs, as in conditions 2–3.
- `claimed_lower_bound`: unsigned decimal string (digits, at most one `.`).
- `claimed_sup` (optional): decimal string; if present the verifier checks it equals M exactly.
- Other keys (`note`, `source`, `value_float`) are informational and ignored.

Scorer in this repo: `problems/bh_cube/core.py` (`python -m problems.bh_cube.core cert.json`). The independent verifier
goes in `verify/bh_cube.py` and must not import from `problems/`.
