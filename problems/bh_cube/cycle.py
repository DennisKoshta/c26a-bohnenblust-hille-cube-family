"""Explicit construction for every d: BH^{<=d} >= 2 (1 - 2^-d)^{(d+1)/(2d)}, with n = 3d variables (no search).

Variables: a cycle group z_0..z_{d-1} and pairs Q_i = (a_i, b_i), i in Z_d. Signs are bits (bit 1 = value -1).
For every proper subset J of Z_d there is one block, with live set A_J given by

    z_i + z_{i+1} = w_J(i)   for i in J      (edge i of the d-cycle on z)
    a_i + b_i     = w_J(i)   for i not in J  (pair Q_i)

where w_J(i) = [i-1 in J] + c_i (mod 2), with c = 0 for d odd and c = e_0 for d even. The block is chi_beta 1_{A_J}.
The coset beta + A_J^perp is the product of {a_i, b_i} over the fixed pairs and, on every path of consecutive edges
of J (m edges, m+1 vertices), the subsets of the path's vertices whose parity is that of m. All of its weights
are <= |J| + (d - |J|) = d.

Proof that f = sum_J chi_beta_J 1_{A_J} is {0,±1}-valued with (2^d - 1) 2^d coefficients of size 2^-d:

* Cosets are disjoint. For J != K, take i in J xor K. Then Q_i is fixed in exactly one of the two blocks, so their
  characters differ on (a_i, b_i), which is {a_i, b_i} in one and empty in the other.
* Live sets are disjoint. Put y = 1_J + 1_K, which is nonzero.
  - If y is not all-ones, the cyclic word y has an i with y_{i-1} = 1 and y_i = 0. Then i lies in both J and K
    (a shared edge) or in neither (a shared pair), and w_J(i) + w_K(i) = y_{i-1} = 1, so the two equations clash.
  - If y is all-ones, then K = Z_d \\ J, no pair is shared, and the only shared nonzero form is
    alpha = sum_{i in J} (z_i + z_{i+1}) = sum_{i not in J} (z_i + z_{i+1}) (the edges of the cycle sum to 0). Its
    value on A_J plus its value on A_K is sum_{i in J} w_J(i) + sum_{i in K} w_K(i)
    = #{i : i-1 and i on the same side} + |c| = d - (an even number of side changes) + |c|, which is odd: a clash.

So ||f||_inf = 1 and ||f^||_p^p = (2^d - 1) 2^d 2^{-dp}, giving R = 2 (1 - 2^-d)^{(d+1)/(2d)}, which tends to 2 from
below. Certificates (exact, via core.py):

    python -m problems.bh_cube.cycle 2 3 4 5 6 [--out DIR]
"""

import argparse
import os
import sys

from problems.bh_cube import core


def blocks(d):
    """Yield (beta, basis, x0) per block: beta a coset representative, basis of A^perp, x0 a live point (bitmasks)."""
    full = (1 << d) - 1
    c = 0 if d % 2 else 1
    z = lambda i: 1 << (i % d)
    qa = lambda i: 1 << (d + 2 * i)
    qb = lambda i: 1 << (d + 2 * i + 1)
    for J in range(full):
        inJ = lambda i: J >> (i % d) & 1
        w = [inJ(i - 1) ^ (c >> i & 1) for i in range(d)]
        beta, basis, x0 = 0, [], 0
        for i in range(d):
            if not inJ(i):
                basis.append(qa(i) | qb(i))
                beta |= qa(i)
                x0 |= qb(i) if w[i] else 0  # a_i = 0, b_i = w_J(i)
        # paths of consecutive edges of J: start at an edge whose predecessor edge is not in J (J is proper)
        for s in range(d):
            if inJ(s) and not inJ(s - 1):
                verts, i, val = [s], s, 0
                while inJ(i):
                    basis.append(z(i) | z(i + 1))
                    val ^= w[i % d]
                    i += 1
                    verts.append(i)
                    x0 |= z(i) if val else 0  # z_s = 0, then propagate z_{i+1} = z_i + w(i)
                if (len(verts) - 1) % 2:  # m odd: tight coset is the odd subsets of the path's vertices
                    beta |= z(s)
        yield beta, basis, x0


def terms(d):
    par = lambda x: bin(x).count("1") & 1
    coef = {}
    for beta, basis, x0 in blocks(d):
        assert len(basis) == d
        V = [0]
        for a in basis:
            V += [v ^ a for v in V]
        assert len(set(V)) == 1 << d
        for al in V:
            m = beta ^ al
            assert m not in coef, "cosets overlap"
            coef[m] = 1 - 2 * par(al & x0)
    return [(tuple(i for i in range(3 * d) if m >> i & 1), a) for m, a in coef.items()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("d", type=int, nargs="+")
    ap.add_argument("--out", default="problems/bh_cube/instances")
    args = ap.parse_args()
    for d in args.d:
        t = terms(d)
        path = os.path.join(args.out, f"cycle_d{d}_n{3 * d}.json")
        lo, M = core.save_cert(path, 3 * d, d, t, note=f"cycle.py: explicit construction, d={d}, n={3 * d}")
        print(f"d={d} n={3 * d} N={len(t)} M={M} R>={lo} formula={2 * (1 - 2.0 ** -d) ** ((d + 1) / (2 * d)):.12f} -> {path}",
              flush=True)
        assert M == 1 << d
    return 0


if __name__ == "__main__":
    sys.exit(main())
