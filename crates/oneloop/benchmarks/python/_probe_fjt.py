#!/usr/bin/env python3
"""Test the FJT index-lowering IDENTITY independently of oneloop's code.

Dotted triangle [2,1,1], lower k=0.  With nu=[1,1,1], the derived IBP identity is
    I([2,1,1]) = diag(d) I([1,1,1]) - sum_{i!=j} (Y^-1)_{0i} nu_j I(nu; j+, i-)
    diag(d) = sum_i (Y^-1)_{0i} (d - sum nu - nu_i) = (d-4) sum_i (Y^-1)_{0i}
At d=4, diag=0 and every sub-integral is finite, so
    I([2,1,1])|_{d=4} = - sum_{i!=j} (Y^-1)_{0i} I(shift)|_{d=4}
Evaluate both sides with exact Feynman-parameter quadrature.  If they agree, the
formula is correct and any oneloop mismatch is an execution bug.
"""
import numpy as np
from crosscheck import scipy_quad     # exact simplex quadrature (n=2,3,4)

# config A, triangle lines 0,1,2
OFF = [[0, 0, 0, 0], [.30, 0, 0, 0], [.10, .35, 0, 0]]
M = [1.0, 1.21, 0.81]


def dist2(i, j):
    d = np.array(OFF[i]) - np.array(OFF[j])
    return float(d @ d)


# modified Cayley Y_ij = m_i^2 + m_j^2 + dist2  (oneloop stores -dist2, subtracts it)
Y = np.array([[M[i] + M[j] + dist2(i, j) for j in range(3)] for i in range(3)])
Yinv = np.linalg.inv(Y)

I_triangle = scipy_quad(OFF, M, [2, 1, 1])       # true LHS

# RHS: sum over ordered (i,j), i!=j, of -(Y^-1)_{0i} * I(shift j+,i-)
# shift = [1,1,1] with nu_j->2, nu_i->0  ==  dotted bubble on the two lines != i
rhs = 0.0
nu = [1, 1, 1]
for i in range(3):
    for j in range(3):
        if i == j:
            continue
        lines = [x for x in range(3) if x != i]          # the surviving bubble lines
        exps = [2 if x == j else 1 for x in lines]        # j raised to 2
        off_b = [OFF[x] for x in lines]
        m_b = [M[x] for x in lines]
        Ib = scipy_quad(off_b, m_b, exps)
        rhs += -Yinv[0, i] * nu[j] * Ib

print(f"I([2,1,1]) true (quad)          = {I_triangle:+.10e}")
print(f"-sum (Y^-1)_0i I(shift) (quad)  = {rhs:+.10e}")
print(f"identity residual               = {abs(I_triangle - rhs):.2e}  "
      f"(rel {abs(I_triangle-rhs)/abs(I_triangle):.2e})")
print()
print(f"oneloop-reduced (from scorecard) = +1.5662133108e-01")
print(f"  -> identity RHS gives {rhs:+.6e}; true {I_triangle:+.6e}")
