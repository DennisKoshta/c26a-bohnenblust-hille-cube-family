"""Tests for verify/bh_cube.py: hand cases, red controls (tampered coefficients, wrong degree, forged claims,
malformed JSON, python -O) and the six committed certificates. Extracted unchanged from the project's
verify/test_verify.py (class TestBhCubeVerifier and its fixtures); only the instance loader reads from disk only.

    python -I verify/test_bh_cube.py
"""
import base64
import copy
import decimal
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest
import zlib

VERIFY_DIR = Path(__file__).resolve().parent
if str(VERIFY_DIR) not in sys.path:
    sys.path.insert(0, str(VERIFY_DIR))

import bh_cube  # noqa: E402


BH_CUBE_D2_N5_JSON = (
    '{"format":"bh_cube/v1","constant":"C26a","n":5,"d":2,'
    '"terms":[[[1],"1"],[[2],"-1"],[[3],"1"],[[4],"1"],[[0,1],"-1"],[[0,2],"1"],'
    '[[0,3],"1"],[[0,4],"1"],[[1,3],"1"],[[1,4],"-1"],[[2,3],"1"],[[2,4],"-1"]],'
    '"claimed_lower_bound":"1.611854897735","claimed_sup":"4",'
    '"note":"LP-linearisation ascent (survey probe bh_quick.py, re-found by ascent.py/soft.py); same value for n=5..8"}\n'
)

BH_CUBE_D3_N8_JSON = (
    '{"format":"bh_cube/v1","constant":"C26a","n":8,"d":3,'
    '"terms":[[[],"1"],[[0],"1"],[[1],"-1"],[[2],"1"],[[4],"1"],[[5],"1"],'
    '[[0,1],"1"],[[0,4],"1"],[[0,6],"1"],[[0,7],"-1"],[[1,3],"-1"],[[1,5],"-1"],'
    '[[1,6],"1"],[[3,4],"-1"],[[3,6],"1"],[[3,7],"-1"],[[4,5],"-1"],[[4,6],"1"],'
    '[[5,6],"1"],[[5,7],"-1"],[[6,7],"-1"],[[0,1,2],"-1"],[[0,1,3],"-1"],'
    '[[0,1,7],"-1"],[[0,2,4],"-1"],[[0,2,6],"-1"],[[0,2,7],"1"],[[0,3,4],"1"],'
    '[[0,3,6],"1"],[[0,3,7],"1"],[[0,4,7],"1"],[[0,6,7],"1"],[[1,2,3],"-1"],'
    '[[1,2,5],"1"],[[1,2,6],"1"],[[1,3,5],"1"],[[1,3,6],"1"],[[1,3,7],"1"],'
    '[[1,5,7],"1"],[[1,6,7],"-1"],[[2,3,4],"-1"],[[2,3,6],"1"],[[2,3,7],"-1"],'
    '[[2,4,5],"1"],[[2,4,6],"1"],[[2,5,6],"-1"],[[2,5,7],"1"],[[2,6,7],"-1"],'
    '[[3,4,5],"-1"],[[3,4,6],"-1"],[[3,4,7],"-1"],[[3,5,6],"1"],[[3,5,7],"1"],'
    '[[4,5,7],"-1"],[[4,6,7],"1"],[[5,6,7],"1"]],'
    '"claimed_lower_bound":"1.829652855011","claimed_sup":"8",'
    '"note":"soft.py n=8 d=3 seed=0"}\n'
)

# Expected metadata and SHA-256 digests for all 6 problems/bh_cube/instances/*.json files:
# (filename, git_blob_sha1, sha256, n, d, num_terms, sup_norm, lo_str, hi_str)
BH_CUBE_INSTANCES_META: tuple[
    tuple[str, str, str, int, int, int, int, str, str], ...
] = (
    (
        "d2_n5.json",
        "f90a7cb5659403df64b3f4a6f11ea02faca8d7c0",
        "13691d36694b3bc1fcc2ed8ba62d785f97c3bb0900c0e9f5f2d91f0a7444385d",
        5,
        2,
        12,
        4,
        "1.6118548977353128793300072",
        "1.6118548977353128793300073",
    ),
    (
        "d3_n8.json",
        "71c0d2e956f4ff019e087e7f65b8eb2ef313b83b",
        "af475381723e51c82340c9d618f99cf9769af2d0c0062e1ae4cb3eaa91fb2e53",
        8,
        3,
        56,
        8,
        "1.8296528550114857586190366",
        "1.8296528550114857586190367",
    ),
    (
        "d4_n11.json",
        "9d3a39b1af9e124cef24f41eaa5d0d1a51950e67",
        "4d303d60ce61ddea5435648bec1cf7410c8bea7992127c906727c7f9122277ad",
        11,
        4,
        240,
        16,
        "1.9209322303565740407806854",
        "1.9209322303565740407806855",
    ),
    (
        "cycle_d5_n15.json",
        "3f1f033c58b0ced7a24da2c753980b834483e17b",
        "bbbd7e89949ffa937ad05f2dd57842b6b5c608d82573121e5a6024d6c2dcb8c5",
        15,
        5,
        992,
        32,
        "1.9622621415505083567270316",
        "1.9622621415505083567270317",
    ),
    (
        "cycle_d6_n18.json",
        "158804fb17f5cbe2a1d24f9d5ec20aa90ef9ebaa",
        "f9b01bcf183adb0aabcead0de136f0512d29f470da06f6dc58eda369c9637459",
        18,
        6,
        4032,
        64,
        "1.9817110515836345508216164",
        "1.9817110515836345508216165",
    ),
    (
        "cycle_d7_n21.json",
        "8ebe4f91b46ae3822d1c83524ef8d3d7b70b8706",
        "5b22ba5c5b60dd5f6074ee6547f26f3727b2cf5064dd40998d578eb201a929d3",
        21,
        7,
        16256,
        128,
        "1.9910564253132065783420630",
        "1.9910564253132065783420631",
    ),
)


def _load_bh_cube_instance_raw(filename: str, blob_sha1: str) -> str:
    """Load a bh_cube instance JSON string from disk (embedded copies for d2/d3 as a fallback)."""
    disk_path = VERIFY_DIR.parent / "problems" / "bh_cube" / "instances" / filename
    if disk_path.exists():
        return disk_path.read_text(encoding="utf-8")

    if filename == "d2_n5.json":
        return BH_CUBE_D2_N5_JSON
    if filename == "d3_n8.json":
        return BH_CUBE_D3_N8_JSON
    raise FileNotFoundError(f"Could not locate bh_cube instance {filename}")


def naive_bh_cube_brute_force(
    n: int, d: int, terms: list[list[object]]
) -> tuple[int, int, Fraction, Fraction]:
    """Direct textbook brute-force definition of deg(f), ||f||_inf, ||f_hat||_{p_d}, and R(f, d).

    Evaluates F(x) = sum_S a_S prod_{i in S} x_i directly over all x in {+1, -1}^n
    via itertools.product without WHT.
    """
    parsed: list[tuple[tuple[int, ...], int]] = []
    deg = 0
    for sub, c_str in terms:
        assert isinstance(sub, list)
        assert isinstance(c_str, str)
        sub_t = tuple(int(i) for i in sub)
        if len(sub_t) > deg:
            deg = len(sub_t)
        parsed.append((sub_t, int(c_str)))

    assert deg <= d, f"Expected deg={deg} <= d={d}"

    max_abs = 0
    for x in itertools.product((1, -1), repeat=n):
        val = 0
        for sub_t, coeff in parsed:
            mono = 1
            for i in sub_t:
                mono *= x[i]
            val += coeff * mono
        if abs(val) > max_abs:
            max_abs = abs(val)

    ctx = decimal.Context(prec=70, rounding=decimal.ROUND_HALF_EVEN)
    p_dec = ctx.divide(decimal.Decimal(2 * d), decimal.Decimal(d + 1))
    inv_p_dec = ctx.divide(decimal.Decimal(d + 1), decimal.Decimal(2 * d))
    sum_p = decimal.Decimal(0)
    for _, coeff in parsed:
        k = abs(coeff)
        if k == 1:
            sum_p = ctx.add(sum_p, decimal.Decimal(1))
        else:
            sum_p = ctx.add(
                sum_p, ctx.exp(ctx.multiply(p_dec, ctx.ln(decimal.Decimal(k))))
            )
    lp_dec = (
        sum_p
        if d == 1
        else ctx.exp(ctx.multiply(inv_p_dec, ctx.ln(sum_p)))
    )
    lp_frac = Fraction(lp_dec)
    ratio_frac = lp_frac / max_abs
    return deg, max_abs, lp_frac, ratio_frac


class TestBhCubeVerifier(unittest.TestCase):
    def test_hand_cases_address_dictator_and_parity(self):
        # 1. Address function at d=2 on n=3 variables:
        #    F(x) = (1 + x_0) x_1 + (1 - x_0) x_2 = x_1 + x_2 + x_0 x_1 - x_0 x_2.
        #    4 terms of +-1, ||F||_inf = 2, p_2 = 4/3, ||F_hat||_{4/3} = 4^{3/4} = 2*sqrt(2),
        #    R(F, 2) = sqrt(2) = 1.4142135623730950488016887...
        addr_d2 = {
            "format": "bh_cube/v1",
            "constant": "C26a",
            "n": 3,
            "d": 2,
            "terms": [
                [[1], "1"],
                [[2], "1"],
                [[0, 1], "1"],
                [[0, 2], "-1"],
            ],
            "claimed_lower_bound": "1.4142135623730950488016887",
            "claimed_sup": "2",
        }
        bf_deg, bf_sup, bf_lp, bf_r = naive_bh_cube_brute_force(
            3, 2, addr_d2["terms"]
        )
        self.assertEqual(bf_deg, 2)
        self.assertEqual(bf_sup, 2)
        self.assertAlmostEqual(float(bf_r), math.sqrt(2.0), places=15)

        for use_np in (True, False):
            if use_np and not bh_cube.HAS_NUMPY:
                continue
            ok, det, msg = bh_cube.verify_bh_cube_obj(addr_d2, use_numpy=use_np)
            self.assertTrue(ok, msg)
            assert det is not None
            self.assertEqual(det["n"], 3)
            self.assertEqual(det["d"], 2)
            self.assertEqual(det["deg"], 2)
            self.assertEqual(det["num_terms"], 4)
            self.assertEqual(det["sup_norm"], 2)
            self.assertEqual(det["f_inf"], Fraction(2, 1))
            self.assertEqual(det["p"], Fraction(4, 3))
            # Rigorous interval must enclose sqrt(2) and 2*sqrt(2):
            self.assertLessEqual(det["ratio_lo"] ** 2, Fraction(2, 1))
            self.assertGreaterEqual(det["ratio_hi"] ** 2, Fraction(2, 1))
            self.assertLessEqual(det["lp_lo"] ** 2, Fraction(8, 1))
            self.assertGreaterEqual(det["lp_hi"] ** 2, Fraction(8, 1))
            self.assertEqual(det["lo_str"], "1.4142135623730950488016887")
            self.assertEqual(det["hi_str"], "1.4142135623730950488016888")
            self.assertEqual(det["lp_lo_str"], "2.8284271247461900976033774")
            self.assertEqual(det["lp_hi_str"], "2.8284271247461900976033775")
            self.assertFalse(det["beats"])
            self.assertNotIn("BEATS", msg)

        # Address function at d=3 on n=6 variables:
        # F(x) = sum_{a,b in {0,1}} (1 + (-1)^a x_0)(1 + (-1)^b x_1) x_{2 + 2a + b}.
        # 16 terms of +-1, ||F||_inf = 4, p_3 = 3/2, ||F_hat||_{3/2} = 16^{2/3} = 4 * 2^{2/3},
        # R(F, 3) = 2^{2/3} = 1.5874010519681994747517056...
        addr_d3_terms = []
        for t_idx in range(4):
            a = (t_idx >> 1) & 1
            b = t_idx & 1
            target = 2 + t_idx
            s0 = -1 if a else 1
            s1 = -1 if b else 1
            addr_d3_terms.append([[target], "1"])
            addr_d3_terms.append([[0, target], "1" if s0 == 1 else "-1"])
            addr_d3_terms.append([[1, target], "1" if s1 == 1 else "-1"])
            addr_d3_terms.append([[0, 1, target], "1" if (s0 * s1) == 1 else "-1"])
        addr_d3 = {
            "format": "bh_cube/v1",
            "constant": "C26a",
            "n": 6,
            "d": 3,
            "terms": addr_d3_terms,
            "claimed_lower_bound": "1.5874010519681994747517056",
            "claimed_sup": "4",
        }
        ok3, det3, msg3 = bh_cube.verify_bh_cube_obj(addr_d3)
        self.assertTrue(ok3, msg3)
        assert det3 is not None
        self.assertEqual(det3["deg"], 3)
        self.assertEqual(det3["sup_norm"], 4)
        self.assertLessEqual(det3["ratio_lo"] ** 3, Fraction(4, 1))
        self.assertGreaterEqual(det3["ratio_hi"] ** 3, Fraction(4, 1))
        self.assertEqual(det3["lo_str"], "1.5874010519681994747517056")
        self.assertEqual(det3["hi_str"], "1.5874010519681994747517057")

        # 2. Dictator cases: f(x) = x_0 (deg=1 <= d, ||f||_inf=1, ||f_hat||_p=1, R=1).
        for n_val, d_val in ((1, 1), (2, 2), (4, 3)):
            dict_cert = {
                "format": "bh_cube/v1",
                "constant": "C26a",
                "n": n_val,
                "d": d_val,
                "terms": [[[0], "1"]],
                "claimed_lower_bound": "1",
                "claimed_sup": "1",
            }
            for use_np in (True, False):
                if use_np and not bh_cube.HAS_NUMPY:
                    continue
                ok_d, det_d, msg_d = bh_cube.verify_bh_cube_obj(
                    dict_cert, use_numpy=use_np
                )
                self.assertTrue(ok_d, msg_d)
                assert det_d is not None
                self.assertEqual(det_d["deg"], 1)
                self.assertEqual(det_d["sup_norm"], 1)
                self.assertEqual(det_d["lp_lo"], Fraction(1, 1))
                self.assertEqual(det_d["lp_hi"], Fraction(1, 1))
                self.assertEqual(det_d["ratio_lo"], Fraction(1, 1))
                self.assertEqual(det_d["ratio_hi"], Fraction(1, 1))
                self.assertEqual(det_d["lo_str"], "1.0000000000000000000000000")

        # 3. Parity cases: f(x) = prod_{i=0}^{d-1} x_i (deg=d, ||f||_inf=1, ||f_hat||_p=1, R=1).
        for d_val in (1, 2, 3, 4):
            parity_cert = {
                "format": "bh_cube/v1",
                "constant": "C26a",
                "n": d_val,
                "d": d_val,
                "terms": [[list(range(d_val)), "-1"]],
                "claimed_lower_bound": "1.0",
                "claimed_sup": "1",
            }
            for use_np in (True, False):
                if use_np and not bh_cube.HAS_NUMPY:
                    continue
                ok_p, det_p, msg_p = bh_cube.verify_bh_cube_obj(
                    parity_cert, use_numpy=use_np
                )
                self.assertTrue(ok_p, msg_p)
                assert det_p is not None
                self.assertEqual(det_p["deg"], d_val)
                self.assertEqual(det_p["sup_norm"], 1)
                self.assertEqual(det_p["ratio_lo"], Fraction(1, 1))
                self.assertEqual(det_p["ratio_hi"], Fraction(1, 1))
                self.assertEqual(det_p["lo_str"], "1.0000000000000000000000000")

    def test_all_instances_d2_to_d7(self):
        for (
            fname,
            blob_sha1,
            expected_sha256,
            exp_n,
            exp_d,
            exp_terms,
            exp_sup,
            exp_lo_str,
            exp_hi_str,
        ) in BH_CUBE_INSTANCES_META:
            raw = _load_bh_cube_instance_raw(fname, blob_sha1)
            self.assertEqual(
                hashlib.sha256(raw.encode("utf-8")).hexdigest(),
                expected_sha256,
                f"SHA-256 mismatch for {fname}",
            )

            ok, det, msg = bh_cube.verify_bh_cube_json(raw)
            self.assertTrue(ok, f"{fname}: {msg}")
            assert det is not None
            self.assertEqual(det["n"], exp_n)
            self.assertEqual(det["d"], exp_d)
            self.assertEqual(det["deg"], exp_d)
            self.assertEqual(det["num_terms"], exp_terms)
            self.assertEqual(det["num_terms"], (1 << exp_d) * ((1 << exp_d) - 1))
            self.assertEqual(det["sup_norm"], exp_sup)
            self.assertEqual(det["sup_norm"], 1 << exp_d)
            self.assertEqual(det["lo_str"], exp_lo_str)
            self.assertEqual(det["hi_str"], exp_hi_str)
            # None of the d=2..7 instances exceed 2, so default --beat="2" must NOT trigger BEATS.
            self.assertFalse(det["beats"])
            self.assertNotIn("BEATS", msg)

            # Check pure-Python path explicitly on instances with n <= 16.
            if exp_n <= bh_cube.MAX_N_PURE_PYTHON:
                ok_py, det_py, msg_py = bh_cube.verify_bh_cube_json(
                    raw, use_numpy=False
                )
                self.assertTrue(ok_py, msg_py)
                assert det_py is not None
                self.assertEqual(det_py["sup_norm"], exp_sup)
                self.assertEqual(det_py["lo_str"], exp_lo_str)
                self.assertEqual(det_py["hi_str"], exp_hi_str)

            # Check --beat behavior: beats claimed_lower_bound, does NOT beat exp_lo_str or exp_hi_str.
            obj = json.loads(raw)
            claimed_lb = obj["claimed_lower_bound"]
            ok_b, det_b, msg_b = bh_cube.verify_bh_cube_json(raw, beat=claimed_lb)
            self.assertTrue(ok_b, msg_b)
            assert det_b is not None
            self.assertTrue(det_b["beats"])
            self.assertIn(f"BEATS {claimed_lb}", msg_b)

            ok_eq, det_eq, msg_eq = bh_cube.verify_bh_cube_json(raw, beat=exp_lo_str)
            self.assertTrue(ok_eq, msg_eq)
            assert det_eq is not None
            self.assertFalse(det_eq["beats"])
            self.assertNotIn("BEATS", msg_eq)

    def test_independent_readme_cycle_construction_d2_to_d6(self):
        # Independently build the n=3d cycle construction for d = 2..6 from the README
        # specification alone:
        #   - variables: d-cycle z_0..z_{d-1} (0..d-1) and d pairs (a_i, b_i) (d..2d-1, 2d..3d-1)
        #   - one block per proper subset J of Z_d:
        #       B_J = prod_{i in J} (z_i + (-1)^{w_J(i)} z_{i+1})
        #             * prod_{i not in J} (a_i + (-1)^{w_J(i)} b_i),
        #     where w_J(i) = [i-1 in J] + c_i (mod 2), c = 0 (d odd), c = e_0 (d even).
        expected_enclosures = {
            2: ("1.6118548977353128793300072", "1.6118548977353128793300073"),
            3: ("1.8296528550114857586190366", "1.8296528550114857586190367"),
            4: ("1.9209322303565740407806854", "1.9209322303565740407806855"),
            5: ("1.9622621415505083567270316", "1.9622621415505083567270317"),
            6: ("1.9817110515836345508216164", "1.9817110515836345508216165"),
        }
        for d, (exp_lo, exp_hi) in expected_enclosures.items():
            n = 3 * d
            coeffs: dict[int, int] = {}
            for j_mask in range((1 << d) - 1):
                factors: list[tuple[tuple[int, int], tuple[int, int]]] = []
                for i in range(d):
                    im1_in_j = (j_mask >> ((i - 1) % d)) & 1
                    c_i = 1 if (d % 2 == 0 and i == 0) else 0
                    w = im1_in_j ^ c_i
                    sgn = -1 if w else 1
                    if (j_mask >> i) & 1:
                        u, v = i, (i + 1) % d
                    else:
                        u, v = d + i, 2 * d + i
                    factors.append(((u, 1), (v, sgn)))
                cur: list[tuple[int, int]] = [(0, 1)]
                for f0, f1 in factors:
                    nxt: list[tuple[int, int]] = []
                    for m, s in cur:
                        nxt.append((m ^ (1 << f0[0]), s * f0[1]))
                        nxt.append((m ^ (1 << f1[0]), s * f1[1]))
                    cur = nxt
                for m, s in cur:
                    self.assertNotIn(m, coeffs, f"Coset collision at d={d}")
                    coeffs[m] = s

            terms = [
                [[i for i in range(n) if (m >> i) & 1], str(coeffs[m])]
                for m in sorted(coeffs)
            ]
            cert = {
                "format": "bh_cube/v1",
                "constant": "C26a",
                "n": n,
                "d": d,
                "terms": terms,
                "claimed_lower_bound": exp_lo,
                "claimed_sup": str(1 << d),
            }
            ok, det, msg = bh_cube.verify_bh_cube_obj(cert)
            self.assertTrue(ok, f"d={d}: {msg}")
            assert det is not None
            self.assertEqual(det["n"], 3 * d)
            self.assertEqual(det["d"], d)
            self.assertEqual(det["deg"], d)
            self.assertEqual(det["num_terms"], ((1 << d) - 1) * (1 << d))
            self.assertEqual(det["sup_norm"], 1 << d)
            self.assertEqual(det["lo_str"], exp_lo)
            self.assertEqual(det["hi_str"], exp_hi)

    def test_random_polynomials_vs_brute_force_and_crosscheck(self):
        rng = random.Random(2026102226)
        tol = Fraction(1, 10**55)

        for n in (2, 3, 4, 5):
            for d in range(1, n + 1):
                all_subs = [
                    list(sub)
                    for k in range(d + 1)
                    for sub in itertools.combinations(range(n), k)
                ]
                chosen = rng.sample(all_subs, k=min(len(all_subs), rng.randint(2, 10)))
                terms = []
                for sub in chosen:
                    c = rng.choice([-7, -3, -2, -1, 1, 2, 3, 5])
                    terms.append([sub, str(c)])

                bf_deg, bf_sup, bf_lp, bf_r = naive_bh_cube_brute_force(n, d, terms)
                cert = {
                    "format": "bh_cube/v1",
                    "constant": "C26a",
                    "n": n,
                    "d": d,
                    "terms": terms,
                    "claimed_lower_bound": "0",
                    "claimed_sup": str(bf_sup),
                }

                ok_py, det_py, msg_py = bh_cube.verify_bh_cube_obj(
                    cert, use_numpy=False
                )
                self.assertTrue(ok_py, msg_py)
                assert det_py is not None
                self.assertEqual(det_py["deg"], bf_deg)
                self.assertEqual(det_py["sup_norm"], bf_sup)
                self.assertLessEqual(det_py["lp_lo"] - tol, bf_lp)
                self.assertLessEqual(bf_lp, det_py["lp_hi"] + tol)
                self.assertLessEqual(det_py["ratio_lo"] - tol, bf_r)
                self.assertLessEqual(bf_r, det_py["ratio_hi"] + tol)

                if bh_cube.HAS_NUMPY:
                    ok_np, det_np, msg_np = bh_cube.verify_bh_cube_obj(
                        cert, use_numpy=True
                    )
                    self.assertTrue(ok_np, msg_np)
                    assert det_np is not None
                    self.assertEqual(det_np["deg"], det_py["deg"])
                    self.assertEqual(det_np["sup_norm"], det_py["sup_norm"])
                    self.assertEqual(det_np["lo_str"], det_py["lo_str"])
                    self.assertEqual(det_np["hi_str"], det_py["hi_str"])
                    self.assertEqual(msg_np, msg_py)

        # Cross-check pure-Python vs numpy at n=14 boundary.
        if bh_cube.HAS_NUMPY:
            n14 = 14
            terms14 = [
                [[0, 3], "2"],
                [[1, 5, 13], "-3"],
                [[2, 7, 11], "1"],
                [[4, 8, 12], "5"],
            ]
            cert14 = {
                "format": "bh_cube/v1",
                "constant": "C26a",
                "n": n14,
                "d": 3,
                "terms": terms14,
                "claimed_lower_bound": "0.5",
            }
            ok_np14, det_np14, _ = bh_cube.verify_bh_cube_obj(cert14, use_numpy=True)
            ok_py14, det_py14, _ = bh_cube.verify_bh_cube_obj(cert14, use_numpy=False)
            self.assertTrue(ok_np14)
            self.assertTrue(ok_py14)
            assert det_np14 is not None and det_py14 is not None
            self.assertEqual(det_np14["sup_norm"], det_py14["sup_norm"])
            self.assertEqual(det_np14["lo_str"], det_py14["lo_str"])

            # Verify large coefficients whose cube sum or Parseval sum exceeds 2^63 - 1
            # do not overflow int64 in the numpy path:
            for large_cert, exp_sup in (
                (
                    {
                        "format": "bh_cube/v1",
                        "constant": "C26a",
                        "n": 4,
                        "d": 1,
                        "terms": [[[], str(1 << 61)]],
                        "claimed_lower_bound": "1",
                    },
                    1 << 61,
                ),
                (
                    {
                        "format": "bh_cube/v1",
                        "constant": "C26a",
                        "n": 20,
                        "d": 1,
                        "terms": [[[0], "3000000"], [[1], "1"]],
                        "claimed_lower_bound": "1",
                    },
                    3000001,
                ),
            ):
                ok_lg, det_lg, msg_lg = bh_cube.verify_bh_cube_obj(
                    large_cert, use_numpy=True
                )
                self.assertTrue(ok_lg, msg_lg)
                assert det_lg is not None
                self.assertEqual(det_lg["sup_norm"], exp_sup)

    def test_red_controls_and_malformed_json_rejected(self):
        base_d2 = json.loads(BH_CUBE_D2_N5_JSON)

        # Red control 1: Wrong degree (term subset length > d).
        wrong_deg_1 = copy.deepcopy(base_d2)
        wrong_deg_1["terms"].append([[0, 1, 2], "1"])
        ok, _, msg = bh_cube.verify_bh_cube_obj(wrong_deg_1)
        self.assertFalse(ok)
        self.assertTrue(msg.startswith("INVALID:"))
        self.assertIn("wrong degree", msg)

        wrong_deg_2 = json.loads(BH_CUBE_D3_N8_JSON)
        wrong_deg_2["d"] = 2
        ok, _, msg = bh_cube.verify_bh_cube_obj(wrong_deg_2)
        self.assertFalse(ok)
        self.assertTrue(msg.startswith("INVALID:"))
        self.assertIn("wrong degree", msg)

        # Red control 2: Tampered coefficient.
        #   (a) Changing a coefficient sign in d2_n5 increases ||F||_inf from 4 to 6,
        #       failing both claimed_sup="4" and (even without claimed_sup) claimed_lower_bound.
        tampered_sup = copy.deepcopy(base_d2)
        tampered_sup["terms"][0][1] = "-1"
        ok, _, msg = bh_cube.verify_bh_cube_obj(tampered_sup)
        self.assertFalse(ok)
        self.assertTrue(msg.startswith("INVALID:"))
        self.assertIn("claimed_sup", msg)

        tampered_lb = copy.deepcopy(tampered_sup)
        del tampered_lb["claimed_sup"]
        ok, _, msg = bh_cube.verify_bh_cube_obj(tampered_lb)
        self.assertFalse(ok)
        self.assertTrue(msg.startswith("INVALID:"))
        self.assertIn("below claimed_lower_bound", msg)

        #   (b) Changing a coefficient magnitude from "1" to "2" increases ||F||_inf from 4 to 5.
        tampered_mag = copy.deepcopy(base_d2)
        tampered_mag["terms"][0][1] = "2"
        ok, _, msg = bh_cube.verify_bh_cube_obj(tampered_mag)
        self.assertFalse(ok)
        self.assertTrue(msg.startswith("INVALID:"))
        self.assertIn("claimed_sup", msg)

        #   (c) Non-canonical or zero coefficients ("0", "+1", "01", "-01", "1.0", 1).
        for bad_c in ("0", "-0", "+1", "01", "-01", "1.0", "1 ", " 1", 1, True):
            bad_c_obj = copy.deepcopy(base_d2)
            bad_c_obj["terms"][0][1] = bad_c
            ok, _, msg = bh_cube.verify_bh_cube_obj(bad_c_obj)
            self.assertFalse(ok, f"Coefficient {bad_c!r} should have been rejected")
            self.assertTrue(msg.startswith("INVALID:"))

        # Red control 3: Malformed JSON, duplicate keys, non-standard JSON constants.
        for bad_json in (
            "",
            "   ",
            "{",
            '{"format": "bh_cube/v1",}',
            '{"format": "bh_cube/v1", "format": "bh_cube/v1"}',
            '{"format": "bh_cube/v1", "n": NaN}',
            "[" * 5000 + "]" * 5000,
            "\ufeff" + BH_CUBE_D2_N5_JSON,
            BH_CUBE_D2_N5_JSON.replace("\n", "\r\n"),
        ):
            ok, _, msg = bh_cube.verify_bh_cube_json(bad_json)
            self.assertFalse(ok, f"Bad JSON {bad_json[:40]!r} should have been rejected")
            self.assertTrue(msg.startswith("INVALID:"))

        # Red control 4: Unsorted subset, duplicate variable in subset, duplicate subset, out-of-range variable.
        for bad_sub in ([1, 0], [1, 1], [-1], [5], [True]):
            bad_sub_obj = copy.deepcopy(base_d2)
            bad_sub_obj["terms"][0][0] = bad_sub
            ok, _, msg = bh_cube.verify_bh_cube_obj(bad_sub_obj)
            self.assertFalse(ok, f"Subset {bad_sub!r} should have been rejected")
            self.assertTrue(msg.startswith("INVALID:"))

        dup_sub_obj = copy.deepcopy(base_d2)
        dup_sub_obj["terms"].append([[1], "1"])
        ok, _, msg = bh_cube.verify_bh_cube_obj(dup_sub_obj)
        self.assertFalse(ok)
        self.assertIn("duplicate subset", msg)

        # Red control 5: Doctored claimed_lower_bound above true ratio, including forged
        # claims ~1e-20 above R that round to the exact same 53-bit IEEE-754 float().
        #   (a) d2_n5: R = 1.61185489773531287933000725...
        #       Forged claim 1e-20 above R ("1.61185489773531287934") and 1e-25 above lo_str:
        for forged_lb in (
            "1.61185489773531287934",
            "1.6118548977353128793300073",
        ):
            self.assertEqual(float(forged_lb), float("1.6118548977353128793300072"))
            doctored_lb = copy.deepcopy(base_d2)
            doctored_lb["claimed_lower_bound"] = forged_lb
            ok, _, msg = bh_cube.verify_bh_cube_obj(doctored_lb)
            self.assertFalse(ok, f"Forged bound {forged_lb} must be rejected")
            self.assertTrue(msg.startswith("INVALID:"))
            self.assertIn("below claimed_lower_bound", msg)

            # Also verify --beat at forged_lb does not trigger BEATS:
            ok_nb, det_nb, msg_nb = bh_cube.verify_bh_cube_obj(base_d2, beat=forged_lb)
            self.assertTrue(ok_nb)
            assert det_nb is not None
            self.assertFalse(det_nb["beats"])
            self.assertNotIn("BEATS", msg_nb)

        #   (b) Address function at d=2: R = sqrt(2) = 1.41421356237309504880168872...
        #       Forged claim ~1e-20 above sqrt(2): "1.41421356237309504881"
        forged_sqrt2 = "1.41421356237309504881"
        self.assertEqual(float(forged_sqrt2), math.sqrt(2.0))
        addr_forged = {
            "format": "bh_cube/v1",
            "constant": "C26a",
            "n": 3,
            "d": 2,
            "terms": [[[1], "1"], [[2], "1"], [[0, 1], "1"], [[0, 2], "-1"]],
            "claimed_lower_bound": forged_sqrt2,
            "claimed_sup": "2",
        }
        ok, _, msg = bh_cube.verify_bh_cube_obj(addr_forged)
        self.assertFalse(ok)
        self.assertIn("below claimed_lower_bound", msg)

        #   (c) Dictator function: R = 1, forged claim 1 + 1e-20 ("1.00000000000000000001")
        forged_one = "1.00000000000000000001"
        self.assertEqual(float(forged_one), 1.0)
        dict_forged = {
            "format": "bh_cube/v1",
            "constant": "C26a",
            "n": 2,
            "d": 2,
            "terms": [[[0], "1"]],
            "claimed_lower_bound": forged_one,
            "claimed_sup": "1",
        }
        ok, _, msg = bh_cube.verify_bh_cube_obj(dict_forged)
        self.assertFalse(ok)
        self.assertIn("below claimed_lower_bound", msg)

        # Red control 6: Invalid --beat values.
        for bad_beat in ("1e100", "-1", "+2", "02", "2/0", " 2", "2 "):
            ok_b, _, msg_b = bh_cube.verify_bh_cube_json(
                BH_CUBE_D2_N5_JSON, beat=bad_beat
            )
            self.assertFalse(ok_b, f"--beat {bad_beat!r} should have been rejected")
            self.assertTrue(msg_b.startswith("INVALID:"))

    def test_cli_exit_codes_and_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            valid_file = Path(tmpdir) / "d2_n5.json"
            valid_file.write_text(BH_CUBE_D2_N5_JSON, encoding="utf-8")

            # Default --beat ("2"): exit 0, VALID, no BEATS.
            proc_def = subprocess.run(
                [sys.executable, str(VERIFY_DIR / "bh_cube.py"), str(valid_file)],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc_def.returncode, 0)
            self.assertIn(
                "VALID [1.6118548977353128793300072, 1.6118548977353128793300073]",
                proc_def.stdout,
            )
            self.assertNotIn("BEATS", proc_def.stdout)

            # With --beat below lo_str: exit 0, prints BEATS.
            proc_beat = subprocess.run(
                [
                    sys.executable,
                    str(VERIFY_DIR / "bh_cube.py"),
                    str(valid_file),
                    "--beat",
                    "1.611854897735",
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc_beat.returncode, 0)
            self.assertIn("BEATS 1.611854897735", proc_beat.stdout)

            # With --beat equal to lo_str: exit 0, no BEATS.
            proc_eq = subprocess.run(
                [
                    sys.executable,
                    str(VERIFY_DIR / "bh_cube.py"),
                    str(valid_file),
                    "--beat",
                    "1.6118548977353128793300072",
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc_eq.returncode, 0)
            self.assertNotIn("BEATS", proc_eq.stdout)

            # Invalid file (CRLF): exit 1, INVALID.
            crlf_file = Path(tmpdir) / "crlf.json"
            crlf_file.write_bytes(
                BH_CUBE_D2_N5_JSON.replace("\n", "\r\n").encode("utf-8")
            )
            proc_bad = subprocess.run(
                [sys.executable, str(VERIFY_DIR / "bh_cube.py"), str(crlf_file)],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc_bad.returncode, 1)
            self.assertTrue(proc_bad.stdout.strip().startswith("INVALID:"))

            # Verify red controls are enforced explicitly under `python -O` (with assert disabled).
            forged_obj = json.loads(BH_CUBE_D2_N5_JSON)
            forged_obj["claimed_lower_bound"] = "1.61185489773531287934"
            forged_file = Path(tmpdir) / "forged_1e20.json"
            forged_file.write_text(json.dumps(forged_obj) + "\n", encoding="utf-8")
            proc_opt = subprocess.run(
                [sys.executable, "-O", str(VERIFY_DIR / "bh_cube.py"), str(forged_file)],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc_opt.returncode, 1)
            self.assertTrue(proc_opt.stdout.strip().startswith("INVALID:"))
            self.assertIn("below claimed_lower_bound", proc_opt.stdout)


# ---------------------------------------------------------------------------
# Shannon capacity of C_7 (constant 9a) helpers and tests
# ---------------------------------------------------------------------------

# The 40 extension words from step (v) of Polak-Schrijver (2018), arXiv:1808.07438
# Section 3 & Appendix (the last 40 words of the 367-word code R).
_PS367_STEP_V_40: tuple[str, ...] = (
    "00521", "01005", "02533", "03565", "04052", "04365", "04624", "04660",
    "05046", "05225", "10534", "14246", "15435", "22524", "24615", "24651",
    "32046", "34035", "34043", "36525", "40040", "41246", "42530", "43514",
    "45641", "50531", "51456", "52400", "52563", "53050", "53142", "53320",
    "53412", "56340", "61505", "62425", "64154", "64340", "65105", "66025",
)

# The 8 private pairs from Gao (2026), arXiv:2607.27869 Section 4 (page 5),
# with J_0 = {0, 5, 6} -> 'H' and J_1 = {1, 2, 3, 4, 7} -> 'V'.
_GAO_BASE_PAIRS: tuple[tuple[str, str, str], ...] = (
    ("13446", "23546", "H"),
    ("34035", "24635", "V"),
    ("53134", "53235", "V"),
    ("44616", "54606", "V"),
    ("60645", "61655", "V"),
    ("03565", "63505", "H"),
    ("64340", "64246", "H"),
    ("64532", "65531", "V"),
)

# Explicit 10-element independent set in C_7^2 (alpha(C_7^2) = 10, Tandon Fig. 1).
C7_2_ALPHA_10: tuple[str, ...] = (
    "00", "12", "24", "26", "31", "43", "46", "51", "63", "65",
)

# Explicit 33-element independent set in C_7^3 (alpha(C_7^3) = 33, Baumert et al. 1971).
C7_3_ALPHA_33: tuple[str, ...] = (
    "000", "005", "044", "056", "112", "124", "131", "136", "151", "163",
    "200", "205", "243", "255", "302", "330", "335", "350", "414", "416",
    "432", "452", "464", "501", "544", "546", "566", "613", "620", "625",
    "632", "651", "663",
)




if __name__ == "__main__":
    unittest.main()
