"""Exact scorer for C26a (Bohnenblust–Hille on the Boolean cube); certificate format in FORMAT.md.

R(f, d) = ||f^||_{p_d} / ||f||_inf with p_d = 2d/(d+1). f has integer coefficients a_S; the sup is exact (integer
fast Walsh–Hadamard transform over all 2^n points) and the l_p norm is enclosed with mpmath interval arithmetic.
"""

import argparse
import json
import re
import sys
from fractions import Fraction

import numpy as np
from mpmath import iv
from mpmath.libmp import to_rational

FORMAT = "bh_cube/v1"
IV_DPS = 120
iv.dps = IV_DPS
COEF_RE = re.compile(r"-?[1-9][0-9]*\Z")  # FORMAT.md condition 3
DEC_RE = re.compile(r"[0-9]+(\.[0-9]+)?\Z")  # claimed_lower_bound
NAT_RE = re.compile(r"(0|[1-9][0-9]*)\Z")  # claimed_sup


class CertError(ValueError):
    pass


def _need(cond, msg):
    if not cond:
        raise CertError(msg)


def _is_int(x):
    return isinstance(x, int) and not isinstance(x, bool)


def mask(S):
    m = 0
    for i in S:
        m |= 1 << i
    return m


def fwht_int(vec):
    """In-place integer Walsh–Hadamard transform: out[m] = sum_S vec[S] (-1)^{|m & S|} (numpy int64 or object)."""
    a = vec
    h = 1
    N = a.shape[0]
    while h < N:
        a = a.reshape(-1, 2, h)
        x = a[:, 0, :].copy()
        y = a[:, 1, :]
        a[:, 0, :] += y
        a[:, 1, :] = x - y
        a = a.reshape(N)
        h *= 2
    return a


def sup_exact(n, terms):
    """max_x |sum a_S x^S| over {±1}^n, exact. terms: list of (subset, int)."""
    total = sum(abs(a) for _, a in terms)
    dtype = np.int64 if total < 2**62 else object
    v = np.zeros(1 << n, dtype=dtype)
    for S, a in terms:
        v[mask(S)] += a
    F = fwht_int(v)
    return int(np.abs(F).max())


def lp_norm_iv(coeffs, d):
    """Interval enclosure of (sum |a|^p)^{1/p}, p = 2d/(d+1)."""
    p = iv.mpf(2 * d) / (d + 1)
    s = iv.mpf(0)
    for a in coeffs:
        s += iv.exp(p * iv.log(iv.mpf(abs(int(a)))))
    return iv.exp(iv.log(s) / p)


def ratio_iv(n, d, terms, M=None):
    if M is None:
        M = sup_exact(n, terms)
    _need(M > 0, "sup is 0")
    return lp_norm_iv([a for _, a in terms], d) / iv.mpf(M), M


def ratio_float(coeffs, d, M):
    c = np.abs(np.asarray(coeffs, float))
    p = 2 * d / (d + 1)
    return (c ** p).sum() ** (1 / p) / M


def check_terms(n, d, terms):
    _need(_is_int(n) and n >= 1 and _is_int(d) and d >= 1, f"need integers n >= 1, d >= 1 (got {n!r}, {d!r})")
    _need(len(terms) > 0, "no terms")
    seen = set()
    for S, a in terms:
        S = tuple(S)
        _need(all(_is_int(i) for i in S), f"subset entries must be integers: {S!r}")
        _need(all(x < y for x, y in zip(S, S[1:])) and all(0 <= i < n for i in S), f"bad subset {S!r}")
        _need(len(S) <= d, f"subset {S!r} has size > d")
        _need(S not in seen, f"duplicate subset {S!r}")
        seen.add(S)
        _need(_is_int(a) and a != 0, f"coefficient must be a nonzero integer: {a!r}")


def save_cert(path, n, d, terms, note="", digits=12):
    terms = sorted(((tuple(S), int(a)) for S, a in terms if a != 0), key=lambda t: (len(t[0]), t[0]))
    check_terms(n, d, terms)
    R, M = ratio_iv(n, d, terms)
    lo = iv_lower_decimal(R, digits)
    obj = {"format": FORMAT, "constant": "C26a", "n": n, "d": d,
           "terms": [[list(S), str(a)] for S, a in terms],
           "claimed_lower_bound": lo, "claimed_sup": str(M), "note": note}
    with open(path, "w", newline="\n") as fh:
        json.dump(obj, fh, separators=(",", ":"))
        fh.write("\n")
    return lo, M


def iv_bounds(R):
    """The interval's endpoints as exact Fractions (no conversion through 53-bit floats)."""
    a, b = R._mpi_
    return Fraction(*to_rational(a)), Fraction(*to_rational(b))


def _decimal(x, digits, up):
    q = 10 ** digits
    k = -((-x.numerator * q) // x.denominator) if up else (x.numerator * q) // x.denominator
    s = str(k).rjust(digits + 1, "0")
    return s[:-digits] + "." + s[-digits:]


def iv_lower_decimal(R, digits):
    """A decimal string <= the interval's lower end (rounded down, so safe)."""
    return _decimal(iv_bounds(R)[0], digits, up=False)


def iv_upper_decimal(R, digits):
    """A decimal string >= the interval's upper end."""
    return _decimal(iv_bounds(R)[1], digits, up=True)


def _no_dup_keys(pairs):
    keys = [k for k, _ in pairs]
    _need(len(keys) == len(set(keys)), f"duplicate JSON keys in {keys}")
    return dict(pairs)


def load_cert(path):
    """Parse and validate a bh_cube/v1 certificate (FORMAT.md); raises CertError on any violation."""
    with open(path, "rb") as fh:
        raw = fh.read()
    _need(b"\r" not in raw, "line endings must be LF")
    obj = json.loads(raw.decode("utf-8"), object_pairs_hook=_no_dup_keys)
    _need(isinstance(obj, dict), "top level must be an object")
    _need(obj.get("format") == FORMAT and obj.get("constant") == "C26a", "format/constant mismatch")
    n, d = obj.get("n"), obj.get("d")
    _need(_is_int(n) and n >= 1 and _is_int(d) and d >= 1, f"need JSON integers n >= 1, d >= 1 (got {n!r}, {d!r})")
    T = obj.get("terms")
    _need(isinstance(T, list) and T, "terms must be a non-empty list")
    terms = []
    for t in T:
        _need(isinstance(t, list) and len(t) == 2 and isinstance(t[0], list) and isinstance(t[1], str),
              f"bad term {t!r}")
        _need(COEF_RE.match(t[1]) is not None, f"bad coefficient string {t[1]!r}")
        terms.append((tuple(t[0]), int(t[1])))
    check_terms(n, d, terms)
    c = obj.get("claimed_lower_bound")
    _need(isinstance(c, str) and DEC_RE.match(c) is not None, f"bad claimed_lower_bound {c!r}")
    if "claimed_sup" in obj:
        _need(isinstance(obj["claimed_sup"], str) and NAT_RE.match(obj["claimed_sup"]) is not None,
              f"bad claimed_sup {obj['claimed_sup']!r}")
    return obj, n, d, terms


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cert", nargs="+")
    args = ap.parse_args(argv)
    ok = True
    for path in args.cert:
        try:
            obj, n, d, terms = load_cert(path)
            R, M = ratio_iv(n, d, terms)
        except (CertError, ValueError, OSError) as e:
            ok = False
            print(f"{path}: FAIL ({e})")
            continue
        lo, _ = iv_bounds(R)
        good = lo >= Fraction(obj["claimed_lower_bound"]) and ("claimed_sup" not in obj or int(obj["claimed_sup"]) == M)
        ok &= good
        print(f"{path}: n={n} d={d} N={len(terms)} M={M} R in [{iv_lower_decimal(R, 20)}, {iv_upper_decimal(R, 20)}]"
              f" addr={2 ** ((d - 1) / d):.6f} {'OK' if good else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
