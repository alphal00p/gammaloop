#!/usr/bin/env python3
"""End-to-end H->gamma gamma W-loop form factor A_1(tau): reduce -> evaluate -> assemble (see docs/10)."""
import math
import os
import subprocess
import sys
from fractions import Fraction

import sympy as sp

import oneloop_bridge as ob
from wloop_numerator import build

d = sp.Symbol("d")
eps = sp.Symbol("eps")
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
DELTA = Fraction(1, 100000)          # off-shell regulator; A_1 is the delta->0 limit


def run_reduce(s, mwsq, delta):
    a = [s.numerator, s.denominator, mwsq.numerator, mwsq.denominator,
         delta.numerator, delta.denominator]
    out = subprocess.run(
        ["cargo", "run", "--release", "--quiet", "--example", "wloop_reduce",
         "-p", "oneloop", "--", *map(str, a)],
        cwd=REPO, capture_output=True, text=True, check=True,
        env={**os.environ, "CARGO_BUILD_JOBS": "2"},
    ).stdout
    red = {}; cur = None
    for line in out.splitlines():
        t = line.split()
        if not t:
            continue
        if t[0] == "MONO":
            cur = (int(t[1]), int(t[2]), int(t[3])); red[cur] = []
        elif t[0] == "TERM":
            start = line.index("(", line.index("coeff="))
            depth = 0
            for i in range(start, len(line)):
                depth += (line[i] == "(") - (line[i] == ")")
                if depth == 0:
                    end = i; break
            coeff = line[start+1:end].strip()
            master = line[end+1:].split()          # KIND arg arg ...
            red[cur].append((coeff, master[0], master[1:]))
    return red


def master_laurent(kind, args):
    a = [float(sp.Rational(x)) for x in args]
    r = {"A0": ob.one_point, "B0": ob.two_point, "C0": ob.three_point}[kind](*a)
    return (complex(r.epsilon_minus_2), complex(r.epsilon_minus_1), complex(r.epsilon_0))


def dseries(expr):
    ser = sp.series(expr.subs(d, 4 - 2*eps), eps, 0, 3).removeO()
    p = sp.Poly(ser, eps) if ser.has(eps) else sp.Poly(ser + 0*eps, eps)
    return (complex(p.coeff_monomial(1)), complex(p.coeff_monomial(eps)),
            complex(p.coeff_monomial(eps**2)))


def assemble(comb, red):
    totals = {}
    for (a, b, c), dpoly in comb.items():
        p = sum((0 if abs(v) < 1e-7 else v) * d**i for i, v in enumerate(dpoly))
        for coeff_str, kind, args in red[(a, b, c)]:
            key = (kind, tuple(args))
            totals[key] = totals.get(key, sp.Integer(0)) + p * sp.sympify(coeff_str, locals={"d": d})
    em2 = em1 = e0 = 0.0
    for (kind, args), coeff in totals.items():
        M2, M1, M0 = master_laurent(kind, list(args))
        a0, a1, a2 = dseries(coeff)
        em2 += a0*M2; em1 += a0*M1 + a1*M2; e0 += a0*M0 + a1*M1 + a2*M2
    return em2, em1, e0


def analytic_A1(tau):
    f = math.asin(math.sqrt(tau))**2
    return -(2*tau**2 + 3*tau + 3*(2*tau - 1)*f) / tau**2


def main():
    print(f"{'tau':>7} {'poles(1/eps)':>13} {'A_reduced':>12} {'A_analytic':>12} {'rel.err':>10}")
    worst = 0.0
    for tau in (0.6046, 0.3, 0.45, 0.8, 0.2):
        comb, s, fitres = build(tau)
        s_frac = Fraction(s).limit_denominator(10**7)
        red = run_reduce(s_frac, Fraction(1), DELTA)
        em2, em1, e0 = assemble(comb, red)
        A_red = e0.real / s                      # P.M / [(q1.q2)(d-2)]|_{d=4} = P.M / s
        A_an = analytic_A1(tau)
        rel = abs(A_red - A_an) / abs(A_an)
        worst = max(worst, rel)
        print(f"{tau:>7.4f} {abs(em1):>13.2e} {A_red:>12.5f} {A_an:>12.5f} {rel:>10.1e}")
    print(f"\nmax relative error {worst:.2e}  (residual ~ delta={float(DELTA):.0e}; "
          "-> exact as delta->0)")
    print("PASS: reducer reproduces the H->gamma gamma W-loop form factor A_1(tau)"
          if worst < 1e-3 else "FAIL")

    # Physical H->gamma gamma = top (A_1/2, validated in gg->h to 1e-13) + W (A_1, above).
    mH, mt, mW = 125.0, 173.0, 80.379
    tau_t, tau_W = mH**2/(4*mt**2), mH**2/(4*mW**2)
    comb, s, _ = build(tau_W); red = run_reduce(Fraction(s).limit_denominator(10**7),
                                                Fraction(1), DELTA)
    A1_W = assemble(comb, red)[2].real / s                      # W loop from the reducer
    A_half_t = 2/tau_t**2*(tau_t + (tau_t-1)*math.asin(math.sqrt(tau_t))**2)  # = gg->h reducer value
    F = 3*(2/3)**2 * A_half_t + A1_W                            # N_c Q_t^2 A_1/2 + A_1
    alpha, v = 1/137.036, 246.22
    Gam = alpha**2 * mH**3 / (256*math.pi**3 * v**2) * abs(F)**2
    print(f"\nFull H->gamma gamma form factor (top {3*(2/3)**2*A_half_t:+.3f} + W {A1_W:+.3f}) "
          f"= {F:+.4f}   [analytic -6.489]")
    print(f"Gamma(H->gamma gamma) = {Gam*1e6:.3f} keV   (SM LO ~ 9.1 keV)")


if __name__ == "__main__":
    sys.exit(main())
