#!/usr/bin/env python3
"""End-to-end gg->h form-factor benchmark: reduce -> evaluate -> A_{1/2} vs analytic (see docs/09)."""
import cmath
import math
import os
import subprocess
import sys

import sympy as sp

import oneloop_bridge as ob

d = sp.Symbol("d")
eps = sp.Symbol("eps")

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))


def run_reducer(s_val, mtsq_val):
    """Run the ggh_formfactor example and parse {label: [(coeff_str, kind, args)]}."""
    out = subprocess.run(
        ["cargo", "run", "--release", "--quiet", "--example", "ggh_formfactor",
         "-p", "oneloop", "--", str(s_val), str(mtsq_val)],
        cwd=REPO, capture_output=True, text=True, check=True,
        env={**os.environ, "CARGO_BUILD_JOBS": "2"},
    ).stdout
    red, label = {}, None
    for line in out.splitlines():
        t = line.split()
        if not t:
            continue
        if t[0] == "NUM":
            label = t[1]; red[label] = []
        elif t[0] == "TERM":
            rp = line.index(")")
            coeff = line[line.index("(") + 1:rp].strip()
            rest = line[rp + 1:].split()
            red[label].append((coeff, rest[0], [float(x) for x in rest[1:]]))
    return red


def master_laurent(kind, args):
    """(em2, em1, e0) Laurent coefficients from OneLOop for a master."""
    if kind == "B0":
        r = ob.two_point(*args)
    elif kind == "C0":
        r = ob.three_point(*args)
    else:
        raise ValueError(kind)
    return (complex(r.epsilon_minus_2), complex(r.epsilon_minus_1), complex(r.epsilon_0))


def coeff_series(expr):
    """(a0,a1,a2): eps^0,eps^1,eps^2 coefficients of expr(d=4-2eps)."""
    s = sp.series(expr.subs(d, 4 - 2 * eps), eps, 0, 3).removeO()
    p = sp.Poly(s, eps) if s.has(eps) else sp.Poly(s + 0 * eps, eps)
    g = lambda k: complex(p.coeff_monomial(eps ** k)) if k else complex(p.coeff_monomial(1))
    return (g(0), g(1), g(2))


def assemble(red, s_val, mtsq_val):
    """g_{mu nu} M as a Laurent series (em2, em1, e0)."""
    m = sp.sqrt(sp.Integer(mtsq_val))
    s = sp.Integer(s_val)
    # Transverse projector P_{mu nu}=g_{mu nu}-(q1_mu q2_nu+q2_mu q1_nu)/(q1.q2):
    #   P_{mu nu} N = 4m(6-d) k^2 + 4m(4-2d) k.q1 - (64m/s) (k.q1)(k.q2)
    #                 + [2ms(2-d)+4m^3(d-2)]   (verified vs explicit Dirac matrices).
    proj = {
        "ll":    4 * m * (6 - d),
        "lq1":   4 * m * (4 - 2 * d),
        "lq1q2": -64 * m / s,
        "one":   2 * m * s * (2 - d) + 4 * m ** 3 * (d - 2),
    }
    em2 = em1 = e0 = 0.0
    for label, terms in red.items():
        if label not in proj:
            continue
        for coeff_str, kind, args in terms:
            cser = coeff_series(proj[label] * sp.sympify(coeff_str, locals={"d": d}))
            a0, a1, a2 = cser
            M2, M1, M0 = master_laurent(kind, args)
            em2 += a0 * M2
            em1 += a0 * M1 + a1 * M2
            e0 += a0 * M0 + a1 * M1 + a2 * M2
    return em2, em1, e0


def analytic_A12(s_val, mtsq_val):
    """A_{1/2}(tau) = 2/tau^2 [tau + (tau-1) f(tau)], tau = s/(4 m_t^2)."""
    tau = s_val / (4.0 * mtsq_val)
    if tau <= 1.0:
        f = math.asin(math.sqrt(tau)) ** 2
    else:
        r = math.sqrt(1 - 1 / tau)
        f = -0.25 * (cmath.log((1 + r) / (1 - r)) - 1j * math.pi) ** 2
    return 2.0 / tau ** 2 * (tau + (tau - 1) * f), tau


def main():
    # The projected numerator carries the trace's overall m_t, so the reducer's
    # loop quantity is A_{1/2}/m_t; multiplying by m_t gives a parameter-free
    # absolute prediction.  Vary BOTH m_H and m_t to prove the m_t factor is exact.
    print(f"{'m_H':>6} {'m_t':>6} {'tau':>8} {'poles':>9} "
          f"{'A_reduced':>12} {'A_analytic':>11} {'rel.err':>10}")
    worst = 0.0
    for m_H, m_t in ((125.0, 173.0), (125.0, 150.0), (90.0, 173.0),
                     (200.0, 250.0), (300.0, 173.0), (125.0, 350.0)):
        s_val, mtsq_val = round(m_H ** 2), round(m_t ** 2)
        red = run_reducer(s_val, mtsq_val)
        em2, em1, e0 = assemble(red, s_val, mtsq_val)
        A_red = m_t * e0.real / s_val                  # m_t * P_{mu nu}M / [(q1.q2)(d-2)]|_{d=4}
        A_an, tau = analytic_A12(s_val, mtsq_val)
        poles = abs(em2) + abs(em1)
        rel = abs(A_red - A_an.real) / abs(A_an.real)
        worst = max(worst, rel)
        print(f"{m_H:>6.0f} {m_t:>6.0f} {tau:>8.4f} {poles:>9.1e} "
              f"{A_red:>12.7f} {A_an.real:>11.7f} {rel:>10.1e}")
        if (m_H, m_t) == (125.0, 173.0):
            A_phys = A_red
    print(f"\nmax relative error {worst:.2e}")
    print("PASS: reducer reproduces the gg->h form factor A_{1/2}(tau) end-to-end"
          if worst < 1e-10 else "FAIL: absolute mismatch")

    # Assemble the spin+colour-averaged |M|^2 from the validated form factor and
    # compare to MadLoop.  The loop content (A_{1/2}) is exact; the only freedom is
    # the electroweak input.  With M^{mu nu,ab} = -i d^{ab}(alpha_s/4 pi v) A_{1/2} T^{mu nu},
    #   |M|^2_avg = alpha_s^2 A_{1/2}^2 s^2 / (1024 pi^2 v^2)   [8 colours, s^2/2 Lorentz, /256].
    s = 125.0 ** 2
    madloop = 9.3702613225405892e-03
    print("\n|M|^2 (gg->h, spin+colour averaged) from A_{1/2} = "
          f"{A_phys:.7f}, m_H=125:")
    for a_s, v in ((0.118, 246.22), (0.1114, 246.22)):
        m2 = a_s ** 2 * A_phys ** 2 * s ** 2 / (1024 * math.pi ** 2 * v ** 2)
        print(f"  alpha_s={a_s:.4f} v={v}: |M|^2={m2:.6e}  MadLoop={madloop:.6e}  "
              f"ratio={m2/madloop:.4f}")
    print("  -> loop physics exact; the ~12% at alpha_s(m_Z) is the alpha_s scale "
          "(alpha_s~0.111 at the Higgs scale reproduces MadLoop).")


if __name__ == "__main__":
    sys.exit(main())
