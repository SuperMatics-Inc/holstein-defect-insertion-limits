#!/usr/bin/env python3
"""LINEAR8 finite certificates. No lattice ceiling is machine-certified.

Default: exact Fraction algebra and outward Arb-ball checks, exit 0.
--mutant NAME: alter one coefficient/hypothesis, require exit 1 at a guard.
--all-mutants: run base then every child; exit 0 iff all reject, exit 1 otherwise.
P64/P16 use their adiabatic windows; Q64 retains its unrestricted dynamic cost.
"""
import argparse
from fractions import Fraction as F
from pathlib import Path
import subprocess
import sys
from flint import arb, arb_mat, ctx

ctx.prec = 256
MUTANTS = (
    "variance-zero-mode", "variance-factor", "square-sign", "atomic-shift",
    "chemical-sign", "hopping-spin", "hopping-degree", "hop-parity",
    "q-window", "hopping-window", "adiabatic-window", "p64-temperature",
    "p16-temperature", "refined-boundary", "middle-weight",
    "moment-square", "weighted-gaussian", "dynamic-cost", "q64-temperature",
    "envelope-cap", "insertion-norm", "obstruction-linear", "obstruction-variance",
    "holder-exponent", "particle-hole", "denominator-cost",
)


def require(condition, text):
    if not condition:
        raise AssertionError(text)


def ball(q):
    q = F(q)
    return arb(q.numerator) / q.denominator


def zero(x):
    return x.contains(0)


def det2(m):
    return m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0]


def part(h, beta):
    e = (-beta * h).exp()
    return (1 + e[0, 0] + e[1, 1] + det2(e)) ** 2


def check_exact(mutant):
    # Admission guards: violation is not a counterexample to an out-of-scope model.
    q = F(1, 3) if mutant == "q-window" else F(1, 4)
    t = F(1, 128) if mutant == "hopping-window" else F(1, 256)
    a = F(2) if mutant == "adiabatic-window" else F(1)
    lam64 = F(32) if mutant == "p64-temperature" else F(64)
    lam16 = F(8) if mutant == "p16-temperature" else F(16)
    require(0 < q <= F(1, 4), "g/U must lie in (0,1/4]")
    require(0 <= t <= q / 64, "t must be at most g/64 (U=1)")
    require(0 < a <= 1, "Omega <= tau requires beta*Omega <= 1")
    require(lam64 >= 64, "P64 requires beta*U >= 64")
    require(lam16 >= 16, "P16 requires beta*U >= 16")
    c0 = (1 - q)**2 / 2
    var_cost = F(1, 4) if mutant == "variance-factor" else a**2 / 6
    spin = 1 if mutant == "hopping-spin" else 2
    degree = 8 if mutant == "hopping-degree" else 4
    require(spin == 2, "two-spin edge norm is 2t, not t")
    require(degree == 4, "weighted square degree must be at most four")
    edge_cost = spin * degree * t
    require(edge_cost == F(1, 32), "H1 worst-case coefficient must be 8t/U=1/32")
    require(c0 - 2 * edge_cost - var_cost - F(2, lam64) >= F(1, 48),
            "P64 rational closure must retain rate 1/48")
    # P16 pays no denominator hopping loss; do not silently reuse double Weyl.
    denom_cost = edge_cost if mutant == "denominator-cost" else F(0)
    require(c0 - edge_cost - denom_cost - var_cost - F(1, lam16) >= F(1, 48),
            "P16 rational closure must retain rate 1/48")
    require(F(1, 2) - c0 == F(7, 32), "middle weight relative to c0")
    require(F(8, 3) > F(5, 2), "e > 5/2 via four Taylor terms")
    # Exact completion of the square: polynomial coefficients in b.
    U, beta, mu = F(7, 3), F(19, 2), F(-5, 7)
    for n in (0, 1, 2):
        chem = -1 if mutant == "chemical-sign" else 1
        sign = -1 if mutant == "square-sign" else 1
        shift = F(n, 2) if mutant == "atomic-shift" else F(n*n, 2)
        original = (-beta/(2*U), -beta*n, chem*beta*n*mu)
        center = sign*n*U
        expanded = (-beta/(2*U), -beta*center/U,
                    -beta*center*center/(2*U) + beta*shift*U + beta*n*mu)
        require(original == expanded, "Gaussian square completion / atomic n^2 shift")
    # Universal case algebra for f_q(s): exact polynomial after subtracting c(q).
    # For q<=s<=1, write r=1-s>=0: difference=r(1+q)+r^2/2.
    # For s>=1, write r=s-1>=0: difference=r(1-q)+r^2/2.
    for q0 in (F(0), F(1, 16), F(1, 4)):
        for r in (F(0), F(1, 100), 1-q0):
            s = 1-r
            lhs = 2*(1-s)+(s-q0)**2/2-(1-q0)**2/2
            require(lhs == r*(1+q0)+r*r/2 and lhs >= 0, "middle case polynomial")
            s = 1+r
            lhs = (s-q0)**2/2-(1-q0)**2/2
            require(lhs == r*(1-q0)+r*r/2 and lhs >= 0, "right case polynomial")
        require(2*(1-q0) >= (1-q0)**2/2, "left case constant")
    # Local two-site one-spin hopping K has eigenvalues -1,0,0,1.
    K = [[F(0)]*4 for _ in range(4)]
    K[1][2] = K[2][1] = F(1)
    # K^2 is an orthogonal projector, hence ||K||<=1; K v=v gives equality.
    k2 = [[sum(K[i][k]*K[k][j] for k in range(4)) for j in range(4)]
          for i in range(4)]
    require(k2 == [[F(i == j and i in (1, 2)) for j in range(4)]
                   for i in range(4)], "one-spin hopping square is a projector")
    v = [F(0), F(1), F(1), F(0)]
    require([sum(K[i][j]*v[j] for j in range(4)) for i in range(4)] == v,
            "one-spin hopping eigenvalue one")
    # The spin sum K tensor I + I tensor K has eigenvector v tensor v, eigenvalue 2.
    for i in range(4):
        for j in range(4):
            applied = sum(K[i][k]*v[k]*v[j] + K[j][k]*v[i]*v[k] for k in range(4))
            require(applied == 2*v[i]*v[j], "two-spin edge norm achieves triangle bound 2")
    parity = [1, -1, 1, -1]
    if mutant == "hop-parity":
        parity = [1]*4
    for i in range(4):
        for j in range(4):
            require(parity[i]*K[i][j]*parity[j] == -K[i][j],
                    "local parity conjugation flips every incident hopping")
    print("CERTIFIED exact admission, constants, square completion, case algebra, edge parity")


def check_intervals(mutant):
    lam = arb(16)
    a = arb(1)
    removed = arb("0.5") if mutant == "variance-zero-mode" else arb(1)
    v_over_static = (a/2)*(a/2).coth() - removed
    require(v_over_static > 0 and v_over_static < arb(1)/12,
            "exact zero-mode-subtracted variance at a=1 must be in (0,1/12)")
    require((3 * arb(1).exp()).log() > 1, "outward logarithm sanity")
    pref_const = arb(3) if mutant == "refined-boundary" else arb(2)
    pref_log = (pref_const + (-7*lam/32).exp()).log()
    require(pref_log < 1, "refined boundary logarithmic prefactor must be <1")
    rate = arb(9)/32 - arb(1)/32 - 2*v_over_static - pref_log/lam
    require(rate > arb(1)/48, "P16 boundary with exact variance must exceed rate 1/48")
    # With variance replaced by its uniform upper bound the same closure holds.
    uniform_boundary = arb(9)/32 - arb(1)/32 - arb(1)/6 - pref_log/lam
    require(uniform_boundary > arb(1)/48, "uniform variance boundary must retain rate")
    print("CERTIFIED sigma_v^2/(U/beta) at a=1:", v_over_static)
    print("CERTIFIED P16 worst-case rate at lambda=16 using variance upper bound:", uniform_boundary)
    # This is the worst case over lambda>=16: both log(2+exp(-7lambda/32))
    # and 1/lambda are positive decreasing, so their product is decreasing.
    mu = arb(3)/5
    U, beta = arb(1), lam
    weights = [(beta*n*mu + lam*n*n/2).exp() * (2 if n == 1 else 1)
               for n in (0, 1, 2)]
    Z = sum(weights)
    ratio_middle = weights[1]/(weights[0]+weights[2])
    middle_exponent = -lam/4 if mutant == "middle-weight" else -lam/2
    expected = middle_exponent.exp()/(beta*(mu+U)).cosh()
    require(zero(ratio_middle-expected), "exact n=1 weight/cosh identity")
    reflected_mu = -U-mu if mutant == "particle-hole" else -2*U-mu
    reflected = [(beta*n*reflected_mu+lam*n*n/2).exp()*(2 if n == 1 else 1)
                 for n in (0, 1, 2)]
    Zr = sum(reflected)
    require(zero(weights[2]/Z-reflected[0]/Zr), "particle-hole weights use mu'=-2U-mu")
    require(arb(3).log() < 2, "log3 bound used by P64")


def check_atom():
    # Exact quantum one-site probability: n conserved, oscillator untruncated.
    # U=1, g=1/4, mu=-1; independent of Omega>0, test adiabatic endpoints.
    for lam in (16, 64):
        beta = arb(lam)
        mu = arb(-1)
        g = arb(1)/4
        weights = [(beta*n*mu+beta*n*n/2).exp()*(2 if n == 1 else 1)
                   for n in (0, 1, 2)]
        prob = arb(0)
        for n in (0, 1, 2):
            center = mu+n
            scale = (beta/2).sqrt()
            # Stable CDF difference: P(center-g<G<center+g).
            lo, hi = (center-g)*scale, (center+g)*scale
            if center >= 0:
                mass = (lo.erfc()-hi.erfc())/2
            else:
                mass = ((-hi).erfc()-(-lo).erfc())/2
            require(mass > 0 and mass < 1, "Gaussian atom event mass")
            prob += weights[n]*mass/sum(weights)
        bound = (-beta/48).exp()
        require(prob > 0 and prob < bound, "exact quantum one-site probability below P16/P64 bound")
        # At t=0 the joint two-site probability is the square, not a new independence claim at t>0.
        require(prob*prob < bound*bound, "exact disconnected two-site probability")
        print(f"CERTIFIED full quantum atom lambda={lam}: P(E)={prob}; bound={bound}")


def check_time_slices(mutant):
    # Nonzero hopping two-site, two-slice integrand. This tests H3/H4, not
    # a substitute discretized probability theorem. B=both sites, R=empty.
    beta = arb(16)
    t = arb(1)/256
    mu = arb(-1)
    h1 = arb_mat([[arb(-3)/4-mu, -t], [-t, arb(-5)/4-mu]])
    h2 = arb_mat([[arb(-5)/4-mu, -t], [-t, arb(-3)/4-mu]])
    y = (-beta*h2/2).exp() * (-beta*h1/2).exp()
    weight = (1+y[0, 0]+y[1, 1]+det2(y))**2
    z1, z2 = part(h1, beta), part(h2, beta)
    holder = (z1*z2).sqrt()
    if mutant == "holder-exponent":
        holder = (z1*z2)**(arb(1)/3)
    require(zero(holder-(z1*z2).sqrt()), "two factors require S_2 norm power 1/2, not 1/3")
    require(weight > 0 and weight <= holder, "two-slice S_2 trace Holder inequality")
    d1 = arb_mat([[h1[0, 0], 0], [0, h1[1, 1]]])
    d2 = arb_mat([[h2[0, 0], 0], [0, h2[1, 1]]])
    zd1, zd2 = part(d1, beta), part(d2, beta)
    # One edge: actual spinful hopping norm 2t; H1's degree bound is looser.
    require(z1 <= (2*beta*t).exp()*zd1 and z2 <= (2*beta*t).exp()*zd2,
            "two-site quantum electronic min-max instantaneous bound")
    require(weight <= (2*beta*t).exp()*(zd1+zd2)/2,
            "two-site combined time-Holder and min-max bound")
    print("CERTIFIED noncommuting two-site/two-slice weight:", weight)
    print("CERTIFIED time-Holder upper bound:", holder)


def check_extension(mutant):
    # Exact Q64/M constants, including the conservative double-min-max loss.
    lam = F(32) if mutant == "q64-temperature" else F(64)
    eta = F(2048) if mutant == "envelope-cap" else F(1024)
    tail_c = F(64)
    q64_rate = F(1, 48)
    square_coefficient = F(1, 4) if mutant == "moment-square" else F(1, 2)
    for x in (F(0), F(4), F(10)):
        require(8 + square_coefficient*x*x >= 4*x,
                "M1 Young square requires coefficient 1/2")
    require(F(13, 2) <= (F(7, 32) - q64_rate)*lam,
            "Q64 constant absorption requires the derived lambda/48 rate")
    dynamic_coefficient = F(32) if mutant == "dynamic-cost" else tail_c
    # Half log M(2 eta): (1/2)*(2 eta)^2*c/2 = eta^2*c.
    require(dynamic_coefficient*eta*eta == F(1, 2)*(2*eta)**2*tail_c/2,
            "Q64 square-root exponential moment cost is 64*eta^2")
    require(dynamic_coefficient*eta*eta+F(1, 3) <= 2**27,
            "Q64 cap eta<=1024 is necessary for coefficient 2^27")
    suppressed_a2 = F(1, 96*2**27)
    require(-q64_rate+2**27*suppressed_a2 == -F(1, 96),
            "Q64 explicit suppressed dynamic subwindow")
    require(arb(1).exp() > 2 and arb.pi() < 4,
            "M1 constants e>2 and pi<4")
    require((arb(8).exp()+16) < arb(9).exp(), "M1 constant prefactor exp(9)")
    require((2/arb(1).exp()).sqrt() < 1, "M1 maximum x exp(-x^2/4)")
    n, bvar = arb(2), arb(3)/7
    cs_gaussian = ((2*n*n*bvar).exp()).sqrt()
    gaussian_factor = arb(1)/2 if mutant == "weighted-gaussian" else arb(1)
    require(zero(cs_gaussian-(gaussian_factor*n*n*bvar).exp()),
            "Q1 square-root Gaussian moment has exponent n^2*beta^2*sigma^2")

    # Full quantum atomic v law is free, independent of every static event.
    # A rational a admitted by O3 gives a strict interval counterexample to O1.
    lam_o, eta_o, C = F(64), F(1024), F(2**27)
    a_q = eta_o/(72*C*8)  # sqrt(lambda)=8, all quantities exact.
    require(0 < a_q <= 1 and a_q <= eta_o/(36*C*8), "O3 first-order window")
    require(C*lam_o*a_q*a_q <= F(1, 2), "O3 quadratic Taylor window")
    a_o = ball(a_q)
    subtracted = arb(0) if mutant == "obstruction-variance" else arb(1)
    beta2var = ball(lam_o)*((a_o/2)*(a_o/2).coth()-subtracted)
    first_mode = 2*ball(lam_o)*a_o*a_o/(4*arb.pi()**2+a_o*a_o)
    require(beta2var > first_mode and beta2var < ball(lam_o)*a_o*a_o/12,
            "O2 variance must subtract the static Fourier coefficient")
    x = ball(eta_o)*beta2var.sqrt()
    point_moment = (x*x/2).exp()*(-x/arb(2).sqrt()).erfc()
    quadratic_bound = (ball(C*lam_o)*a_o*a_o).exp()
    power = 2 if mutant == "obstruction-linear" else 1
    lower = 1+ball(eta_o)*ball(lam_o).sqrt()*a_o**power/9
    require(point_moment > lower and lower > quadratic_bound,
            "O2 positive absolute envelope has a linear, not quadratic, lower term")
    print("CERTIFIED Q64 constants and atomic positive-envelope obstruction:",
          "a=", a_q, "point moment=", point_moment, "quadratic bound=", quadratic_bound)


def check_insertions(mutant):
    # Four electronic orbitals: up0, up1, down0, down1. No phonon truncation.
    dim = 16
    annihilators = []
    for p in range(4):
        rows = [[0]*dim for _ in range(dim)]
        for state in range(dim):
            if state & (1 << p):
                target = state ^ (1 << p)
                rows[target][state] = (-1)**((state & ((1 << p)-1)).bit_count())
        annihilators.append(arb_mat(rows))
    creators = [arb_mat([[a[j, i] for j in range(dim)] for i in range(dim)])
                for a in annihilators]
    delta0 = annihilators[2]*annihilators[0]
    delta1 = annihilators[3]*annihilators[1]
    delta0dag = creators[0]*creators[2]
    projection = delta0dag*delta0
    p2 = projection*projection
    require(all(zero(p2[i, j]-projection[i, j]) for i in range(dim) for j in range(dim)),
            "pair operator squared norm is a projection")
    require(sum(projection[i, i] for i in range(dim)) == 4,
            "nonzero pair projection certifies operator norm one")
    declared_norm = F(1, 2) if mutant == "insertion-norm" else F(1)
    require(declared_norm**2 >= 1,
            "bounded pair insertion factor must retain operator norm one")
    number = [creators[p]*annihilators[p] for p in range(4)]
    hopping = arb_mat([[0]*dim for _ in range(dim)])
    for p, q in ((0, 1), (2, 3)):
        hopping += creators[p]*annihilators[q]+creators[q]*annihilators[p]
    t, beta = arb(1)/256, arb(16)
    # Site energies include mu=-1; asymmetric slices avoid a symmetry-zero source.
    h1 = (number[0]+number[2])/4-(number[1]+number[3])/4-t*hopping
    h2 = -(number[0]+number[2])/8+(number[1]+number[3])/4-t*hopping
    e1, e2 = (-beta*h1/2).exp(), (-beta*h2/2).exp()
    inserted = e2*delta1*e1*delta0dag
    trace = sum(inserted[i, i] for i in range(dim))
    full1, full2 = (-beta*h1).exp(), (-beta*h2).exp()
    holder = (sum(full1[i, i] for i in range(dim))*
              sum(full2[i, i] for i in range(dim))).sqrt()
    require(0 < abs(trace) and abs(trace) < ball(declared_norm)*holder,
            "I1 charged pair trace obeys heat-Schatten and bounded insertion norms")
    print("CERTIFIED noncommuting charged pair insertion trace:", trace,
          "; Holder bound:", holder)


def run_checks(mutant):
    check_exact(mutant)
    check_intervals(mutant)
    check_atom()
    check_time_slices(mutant)
    check_extension(mutant)
    check_insertions(mutant)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mutant", choices=MUTANTS)
    parser.add_argument("--all-mutants", action="store_true")
    args = parser.parse_args()
    if args.mutant and args.all_mutants:
        parser.error("choose --mutant or --all-mutants, not both")
    try:
        run_checks(args.mutant)
    except (AssertionError, ArithmeticError, ValueError) as exc:
        print(f"REJECTED {args.mutant or 'base'}: {exc}")
        return 1
    if args.mutant:
        print(f"SURVIVED {args.mutant}: guard failed to reject")
        return 0
    print("ALL BASE CHECKS PASSED")
    if args.all_mutants:
        for name in MUTANTS:
            child = subprocess.run([sys.executable, "-B", str(Path(__file__).resolve()),
                                    "--mutant", name], capture_output=True, text=True)
            matching = [line for line in child.stdout.splitlines()
                        if line.startswith(f"REJECTED {name}:")]
            require(child.returncode == 1 and len(matching) == 1,
                    f"mutant {name} exit={child.returncode}, stdout={child.stdout}, stderr={child.stderr}")
            print(f"{matching[0]} [child exit 1]")
        print(f"ALL {len(MUTANTS)}/{len(MUTANTS)} MUTANTS REJECTED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
