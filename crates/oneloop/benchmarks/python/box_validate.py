#!/usr/bin/env python3
"""Validate the reducer's on-shell massless box reduction (rank-2) vs direct Feynman-MC integration."""
import os
import subprocess
import sys

import numpy as np
import sympy as sp

import oneloop_bridge as ob

d = sp.Symbol("d")
eps = sp.Symbol("eps")
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))

M2, S, T = 1.0, -3.0, -1.7
U = -S - T
# qi.qj (q1=p1,q2=p2,q3=p3): qi^2=0, p1.p2=s/2, p2.p3=t/2, p1.p3=u/2
QQ = {(1, 1): 0., (2, 2): 0., (3, 3): 0., (1, 2): S/2, (2, 3): T/2, (1, 3): U/2}
QQ.update({(b, a): v for (a, b), v in list(QQ.items())})


def sig_dot_q(j, i):                     # sigma_j . qi, sigma_0=0,1=p1,2=p1+p2,3=p1+p2+p3
    return 0. if j == 0 else sum(QQ[(p, i)] for p in {1: [1], 2: [1, 2], 3: [1, 2, 3]}[j])


def direct(num, ns=40_000_000, seed=7):
    rng = np.random.default_rng(seed)
    g = rng.standard_gamma(1., size=(ns, 4))
    x = g / g.sum(1, keepdims=True)                # uniform on the 4-simplex
    F = M2 - S*x[:, 0]*x[:, 2] - T*x[:, 1]*x[:, 3]
    Pq = {i: sum(x[:, j]*sig_dot_q(j, i) for j in range(4)) for i in (1, 2, 3)}
    return num(x, Pq, F).mean()


box_scalar = lambda x, Pq, F: 1.0 / F**2
def box_qq(i, j):
    return lambda x, Pq, F: (-(QQ[(i, j)]/2.0)*F + Pq[i]*Pq[j]) / F**2


def run_reduce(delta_den):
    out = subprocess.run(
        ["cargo", "run", "--release", "--quiet", "--example", "box_reduce", "-p", "oneloop",
         "--", "-3", "1", "-17", "10", "1", "1", "1", str(delta_den)],
        cwd=REPO, capture_output=True, text=True, check=True,
        env={**os.environ, "CARGO_BUILD_JOBS": "2"}).stdout
    red = {}; cur = None
    for line in out.splitlines():
        t = line.split()
        if not t:
            continue
        if t[0] == "MONO":
            cur = tuple(map(int, t[1:5])); red[cur] = []
        elif t[0] == "TERM":
            st = line.index("(", line.index("coeff=")); dep = 0
            for i in range(st, len(line)):
                dep += (line[i] == "(") - (line[i] == ")")
                if dep == 0:
                    en = i; break
            m = line[en+1:].split(); red[cur].append((line[st+1:en].strip(), m[0], m[1:]))
    return red


def mval(kind, args):
    a = [float(sp.Rational(x)) for x in args]
    r = {"A0": ob.one_point, "B0": ob.two_point,
         "C0": ob.three_point, "D0": ob.four_point}[kind](*a)
    return (complex(r.epsilon_minus_2), complex(r.epsilon_minus_1), complex(r.epsilon_0))


def dser(e):
    s = sp.series(e.subs(d, 4-2*eps), eps, 0, 3).removeO()
    p = sp.Poly(s, eps) if s.has(eps) else sp.Poly(s + 0*eps, eps)
    return (complex(p.coeff_monomial(1)), complex(p.coeff_monomial(eps)),
            complex(p.coeff_monomial(eps**2)))


def reduced_value(red, mono):
    em2 = em1 = e0 = 0.
    for cs, kind, args in red[mono]:
        M2_, M1_, M0_ = mval(kind, args)
        a0, a1, a2 = dser(sp.sympify(cs, locals={"d": d}))
        em2 += a0*M2_; em1 += a0*M1_ + a1*M2_; e0 += a0*M0_ + a1*M1_ + a2*M2_
    return em2, em1, e0


def main():
    red = run_reduce(100000)
    _, _, sc_red = reduced_value(red, (0, 0, 0, 0))
    C = sc_red.real / direct(box_scalar)
    print(f"scalar box: measure C = {C:.6f}  (expect 1/6 = {1/6:.6f})")
    print(f"{'tensor':>22} {'reducer':>13} {'direct*C':>13} {'pole':>7} {'ratio':>9}")
    worst = 0.0
    for mono, (i, j), lbl in (((0, 1, 1, 0), (1, 2), "dot(k,q1)dot(k,q2)"),
                              ((0, 1, 0, 1), (1, 3), "dot(k,q1)dot(k,q3)"),
                              ((0, 0, 1, 1), (2, 3), "dot(k,q2)dot(k,q3)"),
                              ((0, 0, 0, 2), (3, 3), "dot(k,q3)^2")):
        _, p1, e0 = reduced_value(red, mono)
        q_dir = direct(box_qq(i, j)) * C
        r = e0.real / q_dir
        worst = max(worst, abs(r - 1))
        print(f"{lbl:>22} {e0.real:>13.6e} {q_dir:>13.6e} {abs(p1):>7.0e} {r:>9.5f}")
    print(f"\nmax deviation {worst:.1e} (Monte-Carlo limited, ~1/sqrt(N))")
    print("PASS: reducer's on-shell box reduction matches direct integration"
          if worst < 3e-3 else "FAIL")


if __name__ == "__main__":
    sys.exit(main())
