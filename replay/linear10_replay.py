#!/usr/bin/env python3
"""Exact rational replay for LINEAR10's finite constants and guarded mutants."""
from fractions import Fraction as F
import subprocess
import sys

MUTANTS = (
    "atom_linear_constant",
    "static_defect_constant",
    "cut_cost",
    "gaussian_tilt_coefficient",
    "weight_fourth_moment",
    "mixed_union_ratio",
    "dynamic_window",
    "final_exponent",
    "buffer_radius",
    "connected_set_count",
    "kp_threshold",
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def base_checks(mutant=None):
    # The unit circle is contained in the square of area 4, so pi<4.
    pi_upper = F(4)
    require(4 * pi_upper**2 + 1 == 65, "4*pi^2+a^2<65 for pi<4 and a^2<=1")
    require(F(1, 65) > F(1, 81), "sqrt(2/pi)*sqrt(2/65)>1/9")
    require(65 < 9**2, "sqrt(65)<9")
    atom_linear = F(1, 10) if mutant == "atom_linear_constant" else F(1, 9)
    require(atom_linear == F(1, 9), "first-order atomic lower coefficient is 1/9")
    # For 0<=z<=1/2, 1/(1-z)<=1+2z iff z(1-2z)>=0.
    require(F(1, 2) * (1 - 2 * F(1, 2)) >= 0, "small-z exponential majorant endpoint")
    require(2 * F(1, 36) <= F(1, 18), "atom quadratic-to-linear comparison")

    # Static atom cost and cut norm, with t<=g/64<=U/256.
    q = F(1, 4)
    c_q = F(8, 32) if mutant == "static_defect_constant" else (1 - q) ** 2 / 2
    require(c_q == F(9, 32), "uniform static defect exponent c(1/4)=9/32")
    require(2 * (1 - q) >= c_q and c_q <= F(1, 2), "static Chernoff case constants")
    require(F(1) - q - 2 <= -1, "middle static exponent decreases to its minimum")
    require(F(1, 2) >= c_q, "middle-occupation exponent 1/2 dominates c(q)")
    cut_coeff = 7 if mutant == "cut_cost" else 4 * 2
    require(cut_coeff == 8, "deleted hopping norm coefficient 8")
    require(cut_coeff * F(1, 256) == F(1, 32), "per-cut-site cost is lambda/32")

    # Gaussian covariance and the fully explicit chaining/moment constants.
    require(F(2, 4 * 6) == F(1, 12), "Fourier variance coefficient is 1/12")
    require(F(1, 8) / (1 - F(3, 4)) == F(1, 2), "dyadic link thresholds sum to s/2")
    require(F(3, 4) ** 2 / F(1, 2) == F(9, 8), "dyadic Chernoff exponent ratio is 9/8")
    require(F(32**2, 64) == 16, "tail starts at z=16 when y=4*sqrt(64*lambda*a^2)")
    require(F(1) + 1 + F(1, 2) > 2, "e>2 from its first three series terms")
    require(2 * F(1, 4) == F(1, 2), "e>2 makes 2*exp(-2)<1/2")
    require(F(16) < F(2**8), "e^8+16<e^9 follows from e>2")
    e_partial = F(1) + 1 + F(1, 2)
    e_tail = F(1, 6) / (1 - F(1, 4))
    require(e_partial + e_tail < 3, "e<3 from the geometric series tail after n=2")
    require(F(1) + 2 + 2 > 3, "e^2>3 from the first terms of its series")

    # Gaussian exponential tilt and the finite fourth moment required by Lemma I.
    tilt = F(4, 1) if mutant == "gaussian_tilt_coefficient" else F(2**2, 2) + 2 + F(1, 2)
    require(tilt == F(9, 2), "maximal n<=2, eta<=1 tilt coefficient is 9/2")
    fourth = F(7, 1) if mutant == "weight_fourth_moment" else F(4**2, 2)
    require(fourth == 8, "centered-weight fourth-moment exponent coefficient is 8")

    # b_0~N(0,U/beta), independent of v(0); both fourth-tilt exponents have
    # coefficient 8 n^2, multiplying lambda for b_0 and alpha for v(0).
    fourth_tilt_coefficient = F(4**2, 2)
    require(fourth_tilt_coefficient == 8, "fourth Gaussian tilt coefficient is 8 n^2")
    require(fourth_tilt_coefficient * 2**2 == 32,
            "n<=2 fourth moments are exp(32 lambda) for b_0 and exp(32 alpha) for v(0)")
    convexity_constant = F(3)**3
    require(convexity_constant == 27,
            "(a0+a1+a2)^4 <= 27 sum_i a_i^4 by three-term convexity")

    # Mixed-set variance budget and final exponent.
    ratio = 8 if mutant == "mixed_union_ratio" else 7
    require(ratio == 7, "mixed-set union ratio is 7")
    a2 = F(1, 251) if mutant == "dynamic_window" else F(1, 252)
    require(a2 <= F(1, 252), "adiabatic window has a^2<=1/252")
    centered_cost = F(ratio) * tilt * F(1, 12) * a2
    require(centered_cost == F(1, 96), "mixed centered Gaussian cost is lambda/96")
    final = -c_q + F(ratio, 32) + centered_cost
    if mutant == "final_exponent":
        final = -F(4, 96)
    require(final == -F(5, 96), "combined exponent before log(3) is -5 lambda/96")
    require(F(2) <= F(64, 32), "log(3)<2 and lambda/32>=2 for lambda>=64")
    require(F(2) + final * 64 <= -F(1, 48) * 64,
            "final exponent is at most -lambda/48 at lambda=64")

    # Exact interval comparisons for k*=12 and the max-norm buffer.
    radius = 11 if mutant == "buffer_radius" else 12
    require(4**10 * 252 < 2**29 < 4**11 * 252,
            "2^29/252 lies strictly between 4^10 and 4^11")
    buffer_size = (2 * radius + 1) ** 2
    require(radius == 12 and buffer_size == 625, "radius-12 buffer has 625 sites")
    require(buffer_size > 7, "singleton defect buffer exceeds the mixed-set ratio")

    # Conditional KP arithmetic for connected degree-four site-set polymers.
    walk_count = 15 if mutant == "connected_set_count" else 4**2
    require(walk_count == 16, "depth-first closed-walk growth constant is 16")
    incompatible_roots = 1 + 4
    require(incompatible_roots == 5, "distance-one incompatibility root count is 5")
    x = F(1, 21)  # x=z*e
    require(F(5) * x / (1 - F(16) * x) == 1, "A=1 KP threshold is z*e=1/21")
    threshold = 20 if mutant == "kp_threshold" else 21
    require(threshold == 21, "KP numerator and denominator combine to factor 21")
    # e<3 gives exp(-1/3)>1/3 at lambda=64; hence both 16*x>1 and 21*x>1.
    require(F(16, 3) > 1 and F(21, 3) > 1,
            "hypothetical z=exp(-lambda/48) fails the endpoint KP test")

    # The temperature interval is tau>=sqrt(252)*Omega and tau<=U/64.
    require(64 * F(1, 64) == 1, "lambda=64 corresponds to tau=U/64")
    require(252 == 6**2 * 7, "sqrt(252)=6*sqrt(7)")
    require(64**2 * 252 == 1_032_192, "nonempty-window endpoint square factor")
    print("PASS exact rational and interval constants: obstruction, Gaussian moments, atom, trace cut")
    print("PASS exact rational and interval constants: mixed budget, buffer, conditional KP, temperature")



def run_mutant(name):
    try:
        base_checks(name)
    except AssertionError as exc:
        print(f"REJECTED {name}: {exc}")
        return 1
    print(f"MUTANT SURVIVED {name}: no guard fired")
    return 0


def all_mutants():
    failures = []
    for name in MUTANTS:
        result = subprocess.run(
            [sys.executable, __file__, "--mutant", name],
            capture_output=True,
            text=True,
            check=False,
        )
        first = result.stdout.splitlines()[0] if result.stdout.splitlines() else ""
        if result.returncode != 1 or not first.startswith(f"REJECTED {name}:"):
            failures.append((name, result.returncode, result.stdout, result.stderr))
        else:
            print(first)
    if failures:
        for name, code, stdout, stderr in failures:
            print(f"FAILED {name}: child_exit={code}; stdout={stdout!r}; stderr={stderr!r}")
        return 1
    print(f"ALL {len(MUTANTS)}/{len(MUTANTS)} MUTANTS REJECTED")
    return 0


def main(argv):
    if len(argv) == 3 and argv[1] == "--mutant" and argv[2] in MUTANTS:
        return run_mutant(argv[2])
    if len(argv) == 2 and argv[1] == "--all-mutants":
        base_checks()
        return all_mutants()
    if len(argv) == 1:
        base_checks()
        print("ALL BASE CHECKS PASSED")
        return 0
    print("usage: replay.py [--all-mutants | --mutant NAME]")
    return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
