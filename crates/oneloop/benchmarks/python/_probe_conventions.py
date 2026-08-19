#!/usr/bin/env python3
"""Empirical convention probe (gitignored, throwaway).

Compute the scalar masters A0/B0/C0/D0 at shared kinematic points with THREE
independent engines and print them side by side, so we can read off the exact
normalization + sign relations before trusting any cross-check:

  * OneLOop (avh_olo, van Hameren)  -> oneloop_bridge, Ellis-Zanderighi norm
  * feynalg (Denner 1993 DE closed forms) -> analytic_A0/B0/C0/D0
  * scipy Euclidean Feynman-parameter MC (finite masters only)

Euclidean massive off-shell config (all separations spacelike, all masses ~1):
every C0/D0 here is UV/IR finite, so the scipy MC and the finite parts are
directly comparable.
"""
import numpy as np
from math import factorial
from scipy.special import gamma
import oneloop_bridge as ob
from feynalg.curated.analytic_integrals import (
    analytic_A0, analytic_B0, analytic_C0, analytic_D0,
)

np.random.seed(20260718)
D = 4.0
NS = 4_000_000
TF = ob.TO_FEYNMAN  # -1/(16 pi^2)

# Euclidean 4D offsets (r_1 = 0) + masses; take first N for the N-gon.
OFF = np.array([
    [0.00, 0.00, 0.00, 0.00],
    [0.30, 0.00, 0.00, 0.00],
    [0.10, 0.35, 0.00, 0.00],
    [0.00, 0.20, 0.40, 0.00],
])
MSQ = np.array([1.0, 1.21, 0.81, 1.44])  # m^2 for each line


def dist2(i, j):
    d = OFF[i] - OFF[j]
    return float(d @ d)          # Euclidean squared distance (>=0)


def mink(i, j):
    return -dist2(i, j)          # Minkowski invariant (spacelike, <=0)


def scipy_scalar(idx):
    """Euclidean scalar n-point (Denner-like norm), finite at d=4."""
    n = len(idx)
    x = np.random.dirichlet(np.ones(n), size=NS)
    F = x @ MSQ[list(idx)]
    for a in range(n):
        for b in range(a + 1, n):
            # spacelike point s_ij = -dist2  ->  F = sum x m^2 - sum x x s_ij = ... + dist2
            F = F + x[:, a] * x[:, b] * dist2(idx[a], idx[b])
    integral = np.mean(F ** (D / 2 - n)) / factorial(n - 1)
    return (-1) ** n * gamma(n - D / 2) * integral


print(f"TO_FEYNMAN = {TF:.10e}   (-1/16pi^2 = {-1/(16*np.pi**2):.10e})\n")

# ---- A0(m^2=1.21) : divergent, compare OneLOop vs feynalg -------------------
m2 = 1.21
a_olo = ob.one_point(m2)
print("A0(1.21):")
print(f"  OneLOop   eps0={a_olo.epsilon_0:+.6e}  em1={a_olo.epsilon_minus_1:+.6e}  (xTF eps0={a_olo.epsilon_0*TF:+.6e})")
print(f"  feynalg   duv=0,mu2=1 -> {analytic_A0(m2):+.6e}   (=m^2(1-ln m^2)={m2*(1-np.log(m2)):+.6e})\n")

# ---- B0(p^2 spacelike; m1,m2) : divergent -----------------------------------
p2, m1, m2b = mink(0, 1), 1.0, 1.21
b_olo = ob.two_point(p2, m1, m2b)
print(f"B0(p2={p2:+.4f}; {m1},{m2b}):")
print(f"  OneLOop   eps0={b_olo.epsilon_0:+.6e}  em1={b_olo.epsilon_minus_1:+.6e}  (xTF eps0={b_olo.epsilon_0*TF:+.6e})")
try:
    print(f"  feynalg   -> {analytic_B0(p2, m1, m2b):+.6e}\n")
except Exception as e:
    print(f"  feynalg   -> ERR {e}\n")

# ---- C0 finite triangle (lines 0,1,2): try BOTH invariant signs -------------
c_scipy = scipy_scalar((0, 1, 2))
c_neg = ob.three_point(-dist2(0,1), -dist2(1,2), -dist2(0,2), MSQ[0], MSQ[1], MSQ[2]).epsilon_0
c_pos = ob.three_point(+dist2(0,1), +dist2(1,2), +dist2(0,2), MSQ[0], MSQ[1], MSQ[2]).epsilon_0
print("C0 finite triangle (lines 0,1,2):")
print(f"  scipy MC              {c_scipy:+.8e}")
print(f"  OneLOop legs=-dist2   {c_neg:+.8e}   scipy/this = {c_scipy/c_neg:+.5f}")
print(f"  OneLOop legs=+dist2   {c_pos:+.8e}   scipy/this = {c_scipy/c_pos:+.5f}\n")

# ---- D0 finite box (lines 0,1,2,3): try BOTH signs --------------------------
d_scipy = scipy_scalar((0, 1, 2, 3))
def box(sgn):
    return ob.four_point(sgn*dist2(0,1), sgn*dist2(1,2), sgn*dist2(2,3), sgn*dist2(3,0),
                         sgn*dist2(0,2), sgn*dist2(1,3),
                         MSQ[0], MSQ[1], MSQ[2], MSQ[3]).epsilon_0
print("D0 finite box (lines 0,1,2,3):")
print(f"  scipy MC              {d_scipy:+.8e}")
print(f"  OneLOop legs=-dist2   {box(-1):+.8e}   scipy/this = {d_scipy/box(-1):+.5f}")
print(f"  OneLOop legs=+dist2   {box(+1):+.8e}   scipy/this = {d_scipy/box(+1):+.5f}")
