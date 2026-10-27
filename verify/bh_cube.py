#!/usr/bin/env python3
"""Independent exact-WHT and interval verifier for the Bohnenblust-Hille cube constant (C_26a / 26a).

Reference:
  Tao et al. optimizationproblems constant 26a (constants/26a.md):
    C_26a = sup_{d >= 1} BH^{<= d}_{{+-1}};
  Bohnenblust, Hille (1931);
  Blei (2001);
  Defant, Mastylo, Perez (2019);
  Eskenazis, Ivanisvili (2022);
  Arunachalam, Dutt, Escudero Gutierrez, Palazuelos (2025).

File format (`bh_cube/v1`, UTF-8 JSON object, LF line endings, no CR/BOM, no duplicate keys):
  - "format":     exactly "bh_cube/v1"
  - "constant":   exactly "C26a"
  - "n":          JSON integer >= 1 (number of Boolean variables, 0-indexed 0..n-1)
  - "d":          JSON integer >= 1 (claimed degree bound)
  - "terms":      non-empty list of [subset, coeff_str] pairs:
                    - subset: strictly increasing list of ints in [0, n-1] with len(subset) <= d
                    - coeff_str: canonical non-zero base-10 integer string (no '+', no leading 0s)
                    - no duplicate subsets
  - "claimed_lower_bound": unsigned decimal string (digits with at most one '.')
  - "claimed_sup" (optional): canonical positive base-10 integer string equal to ||F||_inf
  - Additional informational keys (such as "note", "source", "value_float") are ignored.

Mathematical definitions (confirmed against constants/26a.md and problems/bh_cube/FORMAT.md):
  Given integer Fourier coefficients a_S != 0 for S subseteq {0, ..., n-1} with |S| <= d:
      F(x) = sum_S a_S x^S,   where x in {+1, -1}^n and x^S = prod_{i in S} x_i.
  Indexing x by bitmask m in {0, ..., 2^n - 1} with x_i = (-1)^{(m >> i) & 1}:
      F(m) = sum_S a_S (-1)^{popcount(m & mask(S))}   (exact integer WHT of a).
  Supremum norm over the Boolean cube:
      M = ||F||_inf = max_{x in {+1, -1}^n} |F(x)| = max_{0 <= m < 2^n} |F(m)|   (exact integer).
  Bohnenblust-Hille exponent and Fourier l_{p_d} norm:
      p_d = 2d / (d + 1),
      ||F_hat||_{p_d} = ( sum_S |a_S|^{p_d} )^{1 / p_d} = ( sum_S |a_S|^{2d/(d+1)} )^{(d+1)/(2d)}.
  Certified lower bound on BH^{<= d}_{{+-1}} <= C_26a:
      R(F, d) = ||F_hat||_{p_d} / ||F||_inf.
"""

import argparse
from collections import Counter
import decimal
from fractions import Fraction
import json
import re
import sys
from typing import Any

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:  # pragma: no cover
    np = None  # type: ignore[assignment]
    HAS_NUMPY = False

try:
    from mpmath import iv as mp_iv

    HAS_MPMATH = True
except ImportError:  # pragma: no cover
    mp_iv = None  # type: ignore[assignment]
    HAS_MPMATH = False

MAX_N = 24
MAX_D = 1000
MAX_N_PURE_PYTHON = 16
PURE_PYTHON_CROSSCHECK_MAX_N = 14
INTERVAL_DPS = 100
OUTPUT_DIGITS = 25
MAX_STR_DIGITS = 2000
DEFAULT_BEAT = "2"

_NONZERO_INT_RE = re.compile(r"\A-?[1-9][0-9]*\Z")
_POS_INT_RE = re.compile(r"\A[1-9][0-9]*\Z")
_CLAIMED_BOUND_RE = re.compile(r"\A(?:0|[1-9][0-9]*)(?:\.[0-9]+)?\Z")
_BEAT_RE = re.compile(r"\A(?:0|[1-9][0-9]*)(?:\.[0-9]+|/[1-9][0-9]*)?\Z")


class _VerifierArgumentParser(argparse.ArgumentParser):
    def error(self, message: str):
        print(f"INVALID: {message}")
        raise SystemExit(1)


class _DuplicateKeyError(ValueError):
    pass


def _strict_object_pairs_hook(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    obj: dict[str, Any] = {}
    for k, v in pairs:
        if k in obj:
            raise _DuplicateKeyError(f"duplicate JSON key {k!r}")
        obj[k] = v
    return obj


def _reject_json_constant(val: str) -> Any:
    raise ValueError(f"non-standard JSON constant {val!r} is not allowed")


def _is_strict_int(val: Any) -> bool:
    return isinstance(val, int) and not isinstance(val, bool)


def _mpf_tuple_to_fraction(t: tuple[int, int, int, int]) -> Fraction:
    """Convert an internal mpmath (sign, man, exp, bc) tuple to an exact Fraction."""
    sign, man, exp, _bc = t
    if man == 0:
        return Fraction(0, 1)
    signed_man = -int(man) if sign else int(man)
    if exp >= 0:
        return Fraction(signed_man << int(exp), 1)
    return Fraction(signed_man, 1 << (-int(exp)))


def wht_pure_python(a_list: list[int], n: int) -> list[int]:
    """Compute exact integer Walsh-Hadamard transform F(m) in pure Python."""
    w = list(a_list)
    size = 1 << n
    if len(w) != size:
        raise ValueError(f"Expected length {size}, got {len(w)}")
    half = 1
    while half < size:
        step = half << 1
        for base in range(0, size, step):
            for j in range(base, base + half):
                u = w[j]
                v = w[j + half]
                w[j] = u + v
                w[j + half] = u - v
        half = step
    return w


def wht_numpy(a_arr: "np.ndarray", n: int) -> "np.ndarray":
    """Compute exact integer Walsh-Hadamard transform F(m) using numpy int64."""
    if not HAS_NUMPY:
        raise RuntimeError("numpy is not available")
    size = 1 << n
    if a_arr.shape != (size,):
        raise ValueError(f"Expected shape ({size},), got {a_arr.shape}")
    w = a_arr.astype(np.int64, copy=True)
    half = 1
    for _ in range(n):
        v = w.reshape(-1, 2, half)
        tmp = v[:, 1, :].copy()
        v[:, 1, :] = v[:, 0, :] - tmp
        v[:, 0, :] += tmp
        half <<= 1
    return w


def _exact_int64_sum(arr: "np.ndarray") -> int:
    """Compute exact Python-int sum of an int64 array (length <= 2^30) without int64 overflow."""
    lo = int((arr & np.int64(0x7FFFFFFF)).sum(dtype=np.int64))
    hi = int((arr >> np.int64(31)).sum(dtype=np.int64))
    return (hi << 31) + lo


def _exact_int64_sum_squares(arr: "np.ndarray", sup_norm: int) -> int:
    """Compute exact Python-int sum of squares of an int64 array with |arr| < 2^62."""
    if sup_norm < (1 << 31):
        return _exact_int64_sum(arr * arr)
    lo = arr & np.int64(0x7FFFFFFF)
    hi = arr >> np.int64(31)
    s_hi2 = _exact_int64_sum(hi * hi)
    s_hilo = _exact_int64_sum(hi * lo)
    s_lo2 = _exact_int64_sum(lo * lo)
    return (s_hi2 << 62) + (s_hilo << 32) + s_lo2


def _enclose_lp_and_ratio_mpmath_iv(
    abs_coeff_counts: list[tuple[int, int]],
    d: int,
    sup_norm: int,
    dps: int = INTERVAL_DPS,
) -> tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]:
    """Enclose ||f_hat||_{p_d} and ||f_hat||_{p_d} / ||f||_inf via mpmath.iv in two formulations."""
    old_dps = mp_iv.dps
    try:
        mp_iv.dps = max(int(dps), 100)
        p_iv = mp_iv.mpf(2 * d) / mp_iv.mpf(d + 1)
        inv_p_iv = mp_iv.mpf(d + 1) / mp_iv.mpf(2 * d)
        sup_iv = mp_iv.mpf(sup_norm)

        # Formulation 1: interval power operator (k ** p_iv) and (sum ** inv_p_iv).
        sum_iv_1 = mp_iv.mpf(0)
        # Formulation 2: interval exp/log composition exp(p_iv * log(k)).
        sum_iv_2 = mp_iv.mpf(0)

        for k, count in abs_coeff_counts:
            c_iv = mp_iv.mpf(count)
            if k == 1:
                sum_iv_1 += c_iv
                sum_iv_2 += c_iv
            else:
                k_iv = mp_iv.mpf(k)
                sum_iv_1 += c_iv * (k_iv ** p_iv)
                sum_iv_2 += c_iv * mp_iv.exp(p_iv * mp_iv.log(k_iv))

        if d == 1:
            lp_iv_1 = sum_iv_1
            lp_iv_2 = sum_iv_2
        else:
            lp_iv_1 = sum_iv_1 ** inv_p_iv
            if abs_coeff_counts == [(1, 1)]:
                lp_iv_2 = mp_iv.mpf(1)
            else:
                lp_iv_2 = mp_iv.exp(inv_p_iv * mp_iv.log(sum_iv_2))

        lp1_lo = _mpf_tuple_to_fraction(lp_iv_1._mpi_[0])
        lp1_hi = _mpf_tuple_to_fraction(lp_iv_1._mpi_[1])
        lp2_lo = _mpf_tuple_to_fraction(lp_iv_2._mpi_[0])
        lp2_hi = _mpf_tuple_to_fraction(lp_iv_2._mpi_[1])

        lp_lo = max(lp1_lo, lp2_lo)
        lp_hi = min(lp1_hi, lp2_hi)
        if lp_lo > lp_hi:
            raise AssertionError("mpmath.iv l_p norm enclosures failed to overlap")

        r_iv_1 = lp_iv_1 / sup_iv
        r_iv_2 = lp_iv_2 / sup_iv
        r1_lo = _mpf_tuple_to_fraction(r_iv_1._mpi_[0])
        r1_hi = _mpf_tuple_to_fraction(r_iv_1._mpi_[1])
        r2_lo = _mpf_tuple_to_fraction(r_iv_2._mpi_[0])
        r2_hi = _mpf_tuple_to_fraction(r_iv_2._mpi_[1])

        r_lo = max(r1_lo, r2_lo, lp_lo / sup_norm)
        r_hi = min(r1_hi, r2_hi, lp_hi / sup_norm)
        if r_lo > r_hi:
            raise AssertionError("mpmath.iv ratio enclosures failed to overlap")

        return (lp_lo, lp_hi), (r_lo, r_hi)
    finally:
        mp_iv.dps = old_dps


def _enclose_lp_and_ratio_decimal_iv(
    abs_coeff_counts: list[tuple[int, int]],
    d: int,
    sup_norm: int,
    dps: int = INTERVAL_DPS,
) -> tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]:
    """Enclose ||f_hat||_{p_d} and ||f_hat||_{p_d} / ||f||_inf using stdlib Decimal intervals."""
    prec = max(int(dps), 100) + 20
    ctx_near = decimal.Context(prec=prec, rounding=decimal.ROUND_HALF_EVEN)
    ctx_lo = decimal.Context(prec=prec, rounding=decimal.ROUND_FLOOR)
    ctx_hi = decimal.Context(prec=prec, rounding=decimal.ROUND_CEILING)

    def ln_enclose(x_lo: decimal.Decimal, x_hi: decimal.Decimal) -> tuple[decimal.Decimal, decimal.Decimal]:
        v_lo = ctx_near.ln(x_lo)
        ulp_lo = v_lo.next_plus(context=ctx_near) - v_lo
        pad_lo = ctx_hi.multiply(ulp_lo, decimal.Decimal(4))

        v_hi = ctx_near.ln(x_hi)
        ulp_hi = v_hi.next_plus(context=ctx_near) - v_hi
        pad_hi = ctx_hi.multiply(ulp_hi, decimal.Decimal(4))
        return ctx_lo.subtract(v_lo, pad_lo), ctx_hi.add(v_hi, pad_hi)

    def exp_enclose(x_lo: decimal.Decimal, x_hi: decimal.Decimal) -> tuple[decimal.Decimal, decimal.Decimal]:
        v_lo = ctx_near.exp(x_lo)
        ulp_lo = v_lo.next_plus(context=ctx_near) - v_lo
        pad_lo = ctx_hi.multiply(ulp_lo, decimal.Decimal(4))

        v_hi = ctx_near.exp(x_hi)
        ulp_hi = v_hi.next_plus(context=ctx_near) - v_hi
        pad_hi = ctx_hi.multiply(ulp_hi, decimal.Decimal(4))
        return ctx_lo.subtract(v_lo, pad_lo), ctx_hi.add(v_hi, pad_hi)

    p_lo = ctx_lo.divide(decimal.Decimal(2 * d), decimal.Decimal(d + 1))
    p_hi = ctx_hi.divide(decimal.Decimal(2 * d), decimal.Decimal(d + 1))
    inv_p_lo = ctx_lo.divide(decimal.Decimal(d + 1), decimal.Decimal(2 * d))
    inv_p_hi = ctx_hi.divide(decimal.Decimal(d + 1), decimal.Decimal(2 * d))

    sum_lo = decimal.Decimal(0)
    sum_hi = decimal.Decimal(0)

    for k, count in abs_coeff_counts:
        c_dec = decimal.Decimal(count)
        if k == 1:
            sum_lo = ctx_lo.add(sum_lo, c_dec)
            sum_hi = ctx_hi.add(sum_hi, c_dec)
        else:
            k_dec = decimal.Decimal(k)
            lnk_lo, lnk_hi = ln_enclose(k_dec, k_dec)
            kp_lo, kp_hi = exp_enclose(
                ctx_lo.multiply(p_lo, lnk_lo),
                ctx_hi.multiply(p_hi, lnk_hi),
            )
            sum_lo = ctx_lo.add(sum_lo, ctx_lo.multiply(c_dec, kp_lo))
            sum_hi = ctx_hi.add(sum_hi, ctx_hi.multiply(c_dec, kp_hi))

    if d == 1:
        lp_lo = Fraction(sum_lo)
        lp_hi = Fraction(sum_hi)
    elif abs_coeff_counts == [(1, 1)]:
        lp_lo = Fraction(1, 1)
        lp_hi = Fraction(1, 1)
    else:
        lns_lo, lns_hi = ln_enclose(sum_lo, sum_hi)
        lp_dec_lo, lp_dec_hi = exp_enclose(
            ctx_lo.multiply(inv_p_lo, lns_lo),
            ctx_hi.multiply(inv_p_hi, lns_hi),
        )
        lp_lo = Fraction(lp_dec_lo)
        lp_hi = Fraction(lp_dec_hi)

    r_lo = lp_lo / sup_norm
    r_hi = lp_hi / sup_norm
    return (lp_lo, lp_hi), (r_lo, r_hi)


def enclose_lp_and_ratio(
    abs_coeff_counts: list[tuple[int, int]],
    d: int,
    sup_norm: int,
    dps: int = INTERVAL_DPS,
) -> tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]:
    """Compute rigorous rational enclosures ((lp_lo, lp_hi), (ratio_lo, ratio_hi))."""
    dec_lp, dec_r = _enclose_lp_and_ratio_decimal_iv(
        abs_coeff_counts, d, sup_norm, dps=dps
    )
    if not HAS_MPMATH:
        return dec_lp, dec_r

    mp_lp, mp_r = _enclose_lp_and_ratio_mpmath_iv(
        abs_coeff_counts, d, sup_norm, dps=dps
    )
    if max(mp_lp[0], dec_lp[0]) > min(mp_lp[1], dec_lp[1]):
        raise AssertionError("mpmath.iv and Decimal l_p norm enclosures disagree")
    if max(mp_r[0], dec_r[0]) > min(mp_r[1], dec_r[1]):
        raise AssertionError("mpmath.iv and Decimal ratio enclosures disagree")
    return mp_lp, mp_r


def format_interval_outward(
    lo: Fraction, hi: Fraction, digits: int = OUTPUT_DIGITS
) -> tuple[str, str, Fraction, Fraction]:
    """Format rational interval [lo, hi] (>= 0) with outward rounding to `digits` decimal places."""
    if lo < 0 or hi < lo:
        raise ValueError(f"Expected 0 <= lo <= hi, got [{lo}, {hi}]")
    scale = 10**digits
    scaled_lo = (lo.numerator * scale) // lo.denominator
    scaled_hi = (hi.numerator * scale + hi.denominator - 1) // hi.denominator

    lo_int, lo_frac_part = divmod(scaled_lo, scale)
    hi_int, hi_frac_part = divmod(scaled_hi, scale)
    lo_str = f"{lo_int}.{lo_frac_part:0{digits}d}"
    hi_str = f"{hi_int}.{hi_frac_part:0{digits}d}"
    lo_dec_frac = Fraction(scaled_lo, scale)
    hi_dec_frac = Fraction(scaled_hi, scale)
    return lo_str, hi_str, lo_dec_frac, hi_dec_frac


def _parse_beat_threshold(beat: str | Fraction) -> tuple[Fraction, str]:
    """Parse a --beat threshold into an exact Fraction and canonical display string."""
    if isinstance(beat, Fraction):
        if beat < 0:
            raise ValueError("beat threshold must be non-negative")
        return beat, str(beat)
    if not isinstance(beat, str):
        raise ValueError("beat threshold must be a string or Fraction")
    if not beat or beat != beat.strip() or any(c.isspace() for c in beat):
        raise ValueError("whitespace junk in --beat value")
    if len(beat) > MAX_STR_DIGITS:
        raise ValueError("--beat value exceeds digit limit")
    if not _BEAT_RE.match(beat):
        raise ValueError(
            "expected canonical unsigned decimal or rational 'p/q' for --beat"
        )
    return Fraction(beat), beat


def verify_bh_cube_obj(
    obj: Any,
    beat: str | Fraction | None = DEFAULT_BEAT,
    max_n: int = MAX_N,
    use_numpy: bool | None = None,
    dps: int = INTERVAL_DPS,
) -> tuple[bool, dict[str, Any] | None, str]:
    """Validate a parsed bh_cube/v1 certificate object and compute rigorous enclosures."""
    if not isinstance(obj, dict):
        return (
            False,
            None,
            f"INVALID: top-level JSON value must be an object, got {type(obj).__name__}",
        )

    for req_key in ("format", "constant", "n", "d", "terms", "claimed_lower_bound"):
        if req_key not in obj:
            return False, None, f"INVALID: missing required key {req_key!r}"

    fmt = obj["format"]
    if not isinstance(fmt, str) or fmt != "bh_cube/v1":
        return False, None, f"INVALID: expected format 'bh_cube/v1', got {fmt!r}"

    const = obj["constant"]
    if not isinstance(const, str) or const != "C26a":
        return False, None, f"INVALID: expected constant 'C26a', got {const!r}"

    n = obj["n"]
    if not _is_strict_int(n):
        return False, None, f"INVALID: n must be a JSON integer, got {n!r}"
    if n < 1:
        return False, None, f"INVALID: require n >= 1, got n={n}"
    if n > max_n:
        return (
            False,
            None,
            f"INVALID: n = {n} exceeds limit ({max_n}); too large for this verifier",
        )

    d = obj["d"]
    if not _is_strict_int(d):
        return False, None, f"INVALID: d must be a JSON integer, got {d!r}"
    if d < 1:
        return False, None, f"INVALID: require d >= 1, got d={d}"
    if d > MAX_D:
        return (
            False,
            None,
            f"INVALID: d = {d} exceeds limit ({MAX_D}); too large for this verifier",
        )

    if use_numpy is None:
        use_numpy_flag = HAS_NUMPY
    else:
        use_numpy_flag = bool(use_numpy)
        if use_numpy_flag and not HAS_NUMPY:
            return False, None, "INVALID: numpy requested but not available"

    if not use_numpy_flag and n > MAX_N_PURE_PYTHON:
        return (
            False,
            None,
            f"INVALID: n = {n} exceeds pure-Python limit ({MAX_N_PURE_PYTHON})",
        )

    claimed_str = obj["claimed_lower_bound"]
    if not isinstance(claimed_str, str):
        return (
            False,
            None,
            f"INVALID: claimed_lower_bound must be a string, got {type(claimed_str).__name__}",
        )
    if len(claimed_str) > MAX_STR_DIGITS or not _CLAIMED_BOUND_RE.match(claimed_str):
        return (
            False,
            None,
            f"INVALID: claimed_lower_bound must be an unsigned decimal string, got {claimed_str!r}",
        )
    claimed_frac = Fraction(claimed_str)

    claimed_sup: int | None = None
    if "claimed_sup" in obj:
        cs_raw = obj["claimed_sup"]
        if not isinstance(cs_raw, str):
            return (
                False,
                None,
                f"INVALID: claimed_sup must be a string, got {type(cs_raw).__name__}",
            )
        if len(cs_raw) > MAX_STR_DIGITS or not _POS_INT_RE.match(cs_raw):
            return (
                False,
                None,
                f"INVALID: claimed_sup must be a canonical positive integer string, got {cs_raw!r}",
            )
        claimed_sup = int(cs_raw)

    if beat is None:
        beat = DEFAULT_BEAT
    try:
        beat_frac, beat_str = _parse_beat_threshold(beat)
    except (ValueError, ZeroDivisionError) as exc:
        return False, None, f"INVALID: invalid --beat threshold {beat!r}: {exc}"

    terms_raw = obj["terms"]
    if not isinstance(terms_raw, list) or len(terms_raw) == 0:
        return False, None, "INVALID: terms must be a non-empty list"
    max_terms = 1 << n
    if len(terms_raw) > max_terms:
        return (
            False,
            None,
            f"INVALID: number of terms ({len(terms_raw)}) exceeds 2^n = {max_terms}",
        )

    seen_masks: set[int] = set()
    parsed_terms: list[tuple[int, int]] = []
    deg = 0
    sum_coeffs = 0
    sum_abs_coeffs = 0
    sum_sq_coeffs = 0
    a_empty = 0

    for term_idx, term in enumerate(terms_raw):
        if not isinstance(term, list) or len(term) != 2:
            return (
                False,
                None,
                f"INVALID: term {term_idx}: expected [subset, coeff] of length 2",
            )
        subset, coeff_raw = term[0], term[1]
        if not isinstance(subset, list):
            return (
                False,
                None,
                f"INVALID: term {term_idx}: subset must be a list of integers, got {type(subset).__name__}",
            )

        mask = 0
        prev_v = -1
        for v_idx, v in enumerate(subset):
            if not _is_strict_int(v):
                return (
                    False,
                    None,
                    f"INVALID: term {term_idx}, element {v_idx}: variable index must be a JSON integer, got {v!r}",
                )
            if v < 0 or v >= n:
                return (
                    False,
                    None,
                    f"INVALID: term {term_idx}: variable index {v} out of range [0, {n - 1}]",
                )
            if v <= prev_v:
                return (
                    False,
                    None,
                    f"INVALID: term {term_idx}: subset {subset} is not strictly increasing",
                )
            prev_v = v
            mask |= 1 << v

        sub_len = len(subset)
        if sub_len > d:
            return (
                False,
                None,
                f"INVALID: term {term_idx}: subset {subset} has degree {sub_len} > d={d} (wrong degree)",
            )
        if sub_len > deg:
            deg = sub_len

        if mask in seen_masks:
            return (
                False,
                None,
                f"INVALID: term {term_idx}: duplicate subset {subset}",
            )
        seen_masks.add(mask)

        if not isinstance(coeff_raw, str):
            return (
                False,
                None,
                f"INVALID: term {term_idx}: coefficient must be a string, got {type(coeff_raw).__name__}",
            )
        if len(coeff_raw) > MAX_STR_DIGITS or not _NONZERO_INT_RE.match(coeff_raw):
            return (
                False,
                None,
                f"INVALID: term {term_idx}: coefficient must be a canonical non-zero integer string, got {coeff_raw!r}",
            )

        coeff = int(coeff_raw)
        if mask == 0:
            a_empty = coeff
        sum_coeffs += coeff
        sum_abs_coeffs += abs(coeff)
        sum_sq_coeffs += coeff * coeff
        parsed_terms.append((mask, coeff))

    size = 1 << n
    expected_cube_sum = size * a_empty
    expected_parseval = size * sum_sq_coeffs

    try:
        if use_numpy_flag and sum_abs_coeffs < (1 << 62):
            a_arr = np.zeros(size, dtype=np.int64)
            masks_np = np.array([m for m, _ in parsed_terms], dtype=np.int64)
            coeffs_np = np.array([c for _, c in parsed_terms], dtype=np.int64)
            a_arr[masks_np] = coeffs_np
            f_arr = wht_numpy(a_arr, n)

            if int(f_arr[0]) != sum_coeffs:
                raise AssertionError(
                    f"WHT evaluation at x=(1,...,1) failed: F(0)={int(f_arr[0])} != {sum_coeffs}"
                )
            cube_sum_np = _exact_int64_sum(f_arr)
            if cube_sum_np != expected_cube_sum:
                raise AssertionError(
                    f"WHT cube sum identity failed: sum F(m)={cube_sum_np} != {expected_cube_sum}"
                )

            sup_norm = int(max(int(f_arr.max()), -int(f_arr.min())))
            parseval_np = _exact_int64_sum_squares(f_arr, sup_norm)
            if parseval_np != expected_parseval:
                raise AssertionError(
                    f"Parseval identity failed: sum F(m)^2={parseval_np} != {expected_parseval}"
                )

            # Pure-Python cross-check for n <= 14.
            if n <= PURE_PYTHON_CROSSCHECK_MAX_N:
                a_py = [0] * size
                for m, c in parsed_terms:
                    a_py[m] = c
                f_py = wht_pure_python(a_py, n)
                if f_arr.tolist() != f_py:
                    raise AssertionError("numpy and pure-Python WHT disagreed")
                sup_py = max(max(f_py), -min(f_py))
                if sup_norm != sup_py:
                    raise AssertionError("numpy and pure-Python ||f||_inf disagreed")
        else:
            if not use_numpy_flag or n <= MAX_N_PURE_PYTHON:
                a_py = [0] * size
                for m, c in parsed_terms:
                    a_py[m] = c
                f_py = wht_pure_python(a_py, n)
                if f_py[0] != sum_coeffs:
                    raise AssertionError(
                        f"WHT evaluation at x=(1,...,1) failed: F(0)={f_py[0]} != {sum_coeffs}"
                    )
                if sum(f_py) != expected_cube_sum:
                    raise AssertionError(
                        f"WHT cube sum identity failed: sum F(m)={sum(f_py)} != {expected_cube_sum}"
                    )
                if sum(v * v for v in f_py) != expected_parseval:
                    raise AssertionError("Parseval identity failed in pure-Python WHT")
                sup_norm = max(max(f_py), -min(f_py))
            else:
                return (
                    False,
                    None,
                    "INVALID: coefficient magnitudes exceed int64 WHT bound",
                )
    except MemoryError as exc:  # pragma: no cover
        return False, None, f"INVALID: out of memory verifying n={n}: {exc}"

    if sup_norm <= 0:
        return (
            False,
            None,
            f"INVALID: degenerate function with ||f||_inf = {sup_norm} <= 0",
        )

    if claimed_sup is not None and sup_norm != claimed_sup:
        return (
            False,
            None,
            f"INVALID: computed ||f||_inf = {sup_norm} does not match claimed_sup = {claimed_sup}",
        )

    abs_counts_map = Counter(abs(c) for _, c in parsed_terms)
    abs_coeff_counts = sorted(abs_counts_map.items())

    eff_dps = max(int(dps), len(claimed_str) + 50, len(beat_str) + 50, 100)
    (lp_lo, lp_hi), (ratio_lo, ratio_hi) = enclose_lp_and_ratio(
        abs_coeff_counts, d, sup_norm, dps=eff_dps
    )

    lo_str, hi_str, lo_dec_frac, hi_dec_frac = format_interval_outward(
        ratio_lo, ratio_hi, digits=OUTPUT_DIGITS
    )
    lp_lo_str, lp_hi_str, _, _ = format_interval_outward(
        lp_lo, lp_hi, digits=OUTPUT_DIGITS
    )

    if ratio_lo < claimed_frac:
        return (
            False,
            None,
            f"INVALID: certified lower bound {lo_str} is below claimed_lower_bound {claimed_str}",
        )

    # Print BEATS only if lo > the --beat value (default 2) exactly.
    beats = (lo_dec_frac > beat_frac) and (ratio_lo > beat_frac)

    msg = (
        f"VALID [{lo_str}, {hi_str}] n={n} d={d} deg={deg} "
        f"terms={len(parsed_terms)} sup={sup_norm} lp=[{lp_lo_str}, {lp_hi_str}]"
    )
    if beats:
        msg += f" BEATS {beat_str}"

    details: dict[str, Any] = {
        "n": n,
        "d": d,
        "deg": deg,
        "num_terms": len(parsed_terms),
        "sup_norm": sup_norm,
        "f_inf": Fraction(sup_norm, 1),
        "p": Fraction(2 * d, d + 1),
        "lp_lo": lp_lo,
        "lp_hi": lp_hi,
        "lp_lo_str": lp_lo_str,
        "lp_hi_str": lp_hi_str,
        "ratio_lo": ratio_lo,
        "ratio_hi": ratio_hi,
        "lo_str": lo_str,
        "hi_str": hi_str,
        "lo_dec": lo_dec_frac,
        "hi_dec": hi_dec_frac,
        "claimed_lower_bound": claimed_frac,
        "claimed_sup": claimed_sup,
        "beats": beats,
        "beat_str": beat_str,
    }
    return True, details, msg


def verify_bh_cube_json(
    raw: str,
    beat: str | Fraction | None = DEFAULT_BEAT,
    max_n: int = MAX_N,
    use_numpy: bool | None = None,
    dps: int = INTERVAL_DPS,
) -> tuple[bool, dict[str, Any] | None, str]:
    """Parse a JSON string in bh_cube/v1 format and verify it."""
    if not isinstance(raw, str) or not raw:
        return False, None, "INVALID: empty file"
    if "\r" in raw:
        return False, None, "INVALID: carriage return whitespace junk in file"
    if raw.startswith("\ufeff"):
        return False, None, "INVALID: UTF-8 BOM is not allowed"
    if not raw.strip():
        return False, None, "INVALID: empty file"

    try:
        obj = json.loads(
            raw,
            object_pairs_hook=_strict_object_pairs_hook,
            parse_constant=_reject_json_constant,
        )
    except (
        _DuplicateKeyError,
        json.JSONDecodeError,
        ValueError,
        RecursionError,
        MemoryError,
    ) as exc:
        return False, None, f"INVALID: JSON parse error: {exc}"

    return verify_bh_cube_obj(
        obj, beat=beat, max_n=max_n, use_numpy=use_numpy, dps=dps
    )


def verify_bh_cube_file(
    filepath: str,
    beat: str | Fraction | None = DEFAULT_BEAT,
    max_n: int = MAX_N,
    use_numpy: bool | None = None,
    dps: int = INTERVAL_DPS,
) -> tuple[bool, dict[str, Any] | None, str]:
    """Read a bh_cube/v1 JSON certificate file and verify it."""
    try:
        with open(filepath, "r", encoding="utf-8", newline="") as f:
            raw = f.read()
    except (OSError, UnicodeDecodeError) as exc:
        return False, None, f"INVALID: cannot read file: {exc}"

    return verify_bh_cube_json(
        raw, beat=beat, max_n=max_n, use_numpy=use_numpy, dps=dps
    )


def main(argv: list[str] | None = None) -> int:
    parser = _VerifierArgumentParser(
        description="Verify a Bohnenblust-Hille cube (constant 26a) bh_cube/v1 JSON certificate."
    )
    parser.add_argument("file", help="Path to the bh_cube/v1 JSON certificate file")
    parser.add_argument(
        "--beat",
        type=str,
        default=DEFAULT_BEAT,
        help="Threshold V (default: '2'); prints BEATS only if lo > V exactly",
    )
    args = parser.parse_args(argv)

    ok, _details, message = verify_bh_cube_file(args.file, beat=args.beat)
    print(message)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
