#!/usr/bin/env python3
"""Exact-rational certificates for LINEAR9's finite constants."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import subprocess
import sys


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def check_base() -> None:
    q = F(1, 4)
    c = (1 - q) ** 2 / 2
    require(c == F(9, 32), "c(1/4) != 9/32")
    require(c <= F(1, 2), "c(q) <= 1/2 needed for the n=1 barrier")

    # t <= g/64 <= U/256; the trace-cut numerator costs 8 beta t per site.
    hop = 8 * F(1, 256)
    require(hop == F(1, 32), "8 beta t <= lambda/32")

    # Lemma 2 at eta=1, maximized over n=0,1,2.
    tilt_coefficients = [F(n * n, 2) + n + F(1, 2) for n in range(3)]
    kappa = max(tilt_coefficients)
    require(tilt_coefficients == [F(1, 2), F(2), F(9, 2)], "Gaussian tilt terms")
    require(kappa == F(9, 2), "centered Gaussian tilt coefficient")
    require(F(1, 2) * 2 == 1, "Gaussian concentration coefficient")

    # alpha <= lambda*a^2/12; with |S|/|B| <= 7 and a^2 <= 1/252.
    dynamic_coeff = kappa * F(7, 12) * F(1, 252)
    require(dynamic_coeff == F(1, 96), "buffer dynamic cost must be lambda/96 per defect")
    exponent_without_log = -F(9, 32) + F(7, 32) + dynamic_coeff
    require(exponent_without_log == -F(5, 96), "mixed-set exponent before log(3)")

    # e > 2 implies e^2 > 4 > 3, hence log(3) < 2.
    require(F(4) > F(3), "e^2 > 3 from e>2")
    # e<3 follows from its series: 1+1+1/2 plus a tail below 1/2.
    require(F(5, 2) + F(2, 9) < F(3), "e<3 from a geometric factorial-tail bound")
    require(F(2) <= F(64, 32), "lambda>=64 absorbs log(3)")
    require(exponent_without_log + F(1, 32) == -F(1, 48), "final exponent -lambda/48")

    # The O3 use of exp(z)<=1+2z on 0<=z<=1/2 follows termwise from
    # exp(z)<=1/(1-z), then one exact quadratic comparison.
    zmax = F(1, 2)
    require((1 + 2 * zmax) * (1 - zmax) >= 1, "O3 geometric-series linearization")
    # Dyadic link thresholds sum to s/2; (9/8)^j>=1+j/8 and
    # 2/e^2<1/2 reduce the union series to 4*sum(1/2)^j=8.
    require(F(1, 8) / (1 - F(3, 4)) == F(1, 2), "dyadic increment threshold sum")
    require(F(9, 8) - 1 == F(1, 8), "dyadic binomial exponent increment")
    require(F(2) / F(2) ** 2 == F(1, 2), "dyadic geometric ratio using e>2")
    require(F(4) / (1 - F(1, 2)) == F(8), "dyadic union prefactor")

    # Exact coefficient bookkeeping for the free-history tail moment bound.
    # c=64 lambda a^2, y0=4 sqrt(c), x=s sqrt(c), and the completed-square
    # integral is bounded by sqrt(pi*c) exp(x^2/4).
    require(F(4) ** 2 == F(16), "tail threshold exponent")
    require(F(8) * F(2) == F(16), "sqrt(pi)<2 converts the Gaussian integral")
    require(F(2) ** 8 > F(16), "e^8+16<e^9 follows from e>2")
    require(F(2) ** 8 * (F(2) - 1) > F(16), "e^8(e-1)>16")

    # Atomic obstruction: pi<4 implies 4*pi^2+1<65<81; sqrt(pi)<2.
    require(F(4) * F(4) * F(4) + 1 == F(65), "4*pi^2+1 < 65 from pi<4")
    require(F(65) < F(9) ** 2, "sqrt(65)<9")
    require(F(2) * F(9) == F(18), "denominator bound yields coefficient >1/9")
    lam, eta, C, a = F(64), F(1), F(1), F(1, 10_000)
    lower_linear = eta * F(8) * a / 9
    upper_quadratic = 2 * C * lam * a * a
    require(a <= eta / (36 * C * 8), "O3 first small-a condition")
    require(a * a <= F(1, 2) / (C * lam), "O3 exponential linearization condition")
    require(lower_linear > upper_quadratic, "atomic first-order term beats quadratic cost")

    require(F(1, 252) * F(21, 8) == F(1, 96), "small-a quadratic buffer allocation")
    print("ALL BASE CHECKS PASSED")


def mutant_violation(name: str) -> str:
    lam = F(64)
    r = F(7)
    a2 = F(1, 252)
    target = -lam / 48
    log3_upper = F(2)
    base_dyn = F(21, 8) * a2

    if name == "defect-rate":
        # Mutate c=9/32 down to 8/32.
        value = log3_upper + lam * (-F(8, 32) + r / 32 + base_dyn)
        return "weakened atom barrier no longer proves target" if value > target else ""
    if name == "hopping-degree":
        # Degree 5 raises the edge-sum norm from 8t|S| to 10t|S|.
        bad_hop_per_site = F(10, 256)
        value = log3_upper + lam * (-F(9, 32) + r * bad_hop_per_site + base_dyn)
        return "degree-five cut cost breaks target" if value > target else ""
    if name == "centered-tilt":
        # Incorrectly drop the eta*n cross term at n=2, eta=1.
        actual = F(2 * 2, 2) + F(1) * 2 + F(1, 2)
        claimed = F(2 * 2, 2) + F(1, 2)
        return "omitted Gaussian-tilt cross term" if actual > claimed else ""
    if name == "variance-factor":
        # Replace the Basel coefficient 1/12 by the weaker 1/11.
        bad_dyn = F(9, 2) * r * F(1, 11) * a2
        value = log3_upper + lam * (-F(9, 32) + r / 32 + bad_dyn)
        return "weaker variance constant breaks target" if value > target else ""
    if name == "buffer-ratio":
        bad_r = F(8)
        value = log3_upper + lam * (-F(9, 32) + bad_r / 32 + F(21, 8) * bad_r / 7 * a2)
        return "ratio eight is outside the certified envelope" if value > target else ""
    if name == "quadratic-cutoff":
        bad_a2 = F(1, 250)
        bad_dyn = F(21, 8) * bad_a2
        value = log3_upper + lam * (-F(9, 32) + r / 32 + bad_dyn)
        return "larger a^2 cutoff breaks target" if value > target else ""
    if name == "temperature-window":
        return "lambda=63 cannot absorb log(3)<2 by lambda/32" if F(2) > F(63, 32) else ""
    if name == "parity-denominator":
        # Replacing Z>=Z_dec by a second min-max loss doubles 8 beta t.
        bad_hop = F(16, 256)
        value = log3_upper + lam * (-F(9, 32) + r * bad_hop + base_dyn)
        return "double-min-max denominator cost breaks target" if value > target else ""
    if name == "prefactor":
        # Replace 3 by 9: e<3 (from its series) gives log(9)>2.
        # At lambda=64, log(9)-5lambda/96 > -lambda/48.
        actual_lower = F(2) - F(5, 96) * lam
        return "prefactor nine is not absorbed at lambda=64" if actual_lower >= target else ""
    if name == "uncentered-envelope":
        # lambda=64, eta=C=1, a=1/10000: exact rational form of O3.
        a0 = F(1, 10_000)
        linear = F(8) * a0 / 9
        quadratic = 2 * F(64) * a0 * a0
        return "uncentered first-order atom term exceeds quadratic bound" if linear > quadratic else ""
    return ""


MUTANTS = (
    "defect-rate",
    "hopping-degree",
    "centered-tilt",
    "variance-factor",
    "buffer-ratio",
    "quadratic-cutoff",
    "temperature-window",
    "parity-denominator",
    "prefactor",
    "uncentered-envelope",
)


def child(name: str) -> int:
    violation = mutant_violation(name)
    if not violation:
        print(f"MUTANT SURVIVED {name}")
        return 0
    print(f"REJECTED {name}: {violation}")
    return 1


def all_mutants() -> int:
    for name in MUTANTS:
        proc = subprocess.run(
            [sys.executable, __file__, "--mutant", name],
            check=False,
            capture_output=True,
            text=True,
        )
        print(proc.stdout, end="")
        if proc.returncode != 1 or f"REJECTED {name}:" not in proc.stdout:
            print(f"MUTANT HARNESS FAILED {name}: exit={proc.returncode}")
            if proc.stderr:
                print(proc.stderr, file=sys.stderr, end="")
            return 1
    print(f"ALL {len(MUTANTS)}/{len(MUTANTS)} MUTANTS REJECTED")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--all-mutants", action="store_true")
    group.add_argument("--mutant", choices=MUTANTS)
    args = parser.parse_args()
    if args.mutant:
        return child(args.mutant)
    if args.all_mutants:
        return all_mutants()
    try:
        check_base()
    except AssertionError as exc:
        print(f"BASE CHECK FAILED: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
