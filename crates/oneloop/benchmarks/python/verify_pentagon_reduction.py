#!/usr/bin/env python3
"""Independent numeric check of the van Neerven-Vermaseren / Melrose / FJT reduction of
N-point one-loop scalars (N >= 5) to (N-1)-point ones -- pentagon -> boxes, hexagon ->
pentagons, and so on down to boxes.

Claim (scalar, exactly d=4, where the (d-4) dimension-shift remainder vanishes because the
(d+2)-dim N-point is finite for N>=5):

    I_N = sum_{i=1}^N c_i I_{N-1}^(i),

    c_i = (YB^{-1})_{i0} / (YB^{-1})_{00},

where I_{N-1}^(i) drops propagator i, and YB is the BORDERED modified Cayley matrix
(indices 0..N): YB_00 = 0, YB_0i = YB_i0 = 1, YB_ij = Y_ij = m_i^2 + m_j^2 - (r_i-r_j)^2.
These c_i are purely kinematic -- no dependence on d.  (Melrose 1965; FJT hep-ph/9907327.)

Scalar n-point at d=4:  I_n = (-1)^n Gamma(n - d/2) INT_simplex F^{d/2-n} dx,
  F = sum_i x_i m_i^2 - sum_{i<j} x_i x_j (r_i-r_j)^2   (U=1, Euclidean, F>0).
INT_simplex g dx = (1/(n-1)!) * mean of g over Dirichlet(1..1) samples.  Gitignored.
"""
import numpy as np
from scipy.special import gamma
from math import factorial

np.random.seed(20260701)
D = 4.0
NSAMP = 16_000_000

# a pool of Euclidean offsets in 4D (r_1 = 0) and masses; take the first N for the N-gon
OFFSETS = np.array([
    [0.00, 0.00, 0.00, 0.00],
    [0.30, 0.00, 0.00, 0.00],
    [0.10, 0.35, 0.00, 0.00],
    [0.00, 0.20, 0.40, 0.00],
    [0.05, 0.00, 0.15, 0.45],
    [0.22, 0.12, 0.00, 0.18],   # r_6: linearly dependent in 4D -- exactly the N>=5 point
])
MASSES = np.array([1.0, 1.1, 0.9, 1.2, 0.8, 1.05])


def scalar_integral(r, msq, idx):
    """I_n for the sub-topology using propagators in `idx` (Feynman MC at d=D)."""
    n = len(idx)
    x = np.random.dirichlet(np.ones(n), size=NSAMP)
    F = x @ msq[list(idx)]
    for a in range(n):
        for bb in range(a + 1, n):
            d = r[idx[a]] - r[idx[bb]]
            F = F - x[:, a] * x[:, bb] * float(d @ d)
    integral = np.mean(F ** (D / 2 - n)) / factorial(n - 1)
    return (-1) ** n * gamma(n - D / 2) * integral


def verify(N):
    r, msq = OFFSETS[:N], MASSES[:N]
    Y = np.array([[msq[i] + msq[j] - float((r[i] - r[j]) @ (r[i] - r[j]))
                   for j in range(N)] for i in range(N)])
    YB = np.zeros((N + 1, N + 1))
    YB[0, 1:] = 1.0
    YB[1:, 0] = 1.0
    YB[1:, 1:] = Y
    YBinv = np.linalg.inv(YB)
    c = YBinv[1:, 0] / YBinv[0, 0]

    I_N = scalar_integral(r, msq, tuple(range(N)))
    subs = [scalar_integral(r, msq, tuple(j for j in range(N) if j != i)) for i in range(N)]
    rhs = sum(c[i] * subs[i] for i in range(N))
    rel = abs(I_N - rhs) / abs(I_N)

    name = {5: "pentagon -> 5 boxes", 6: "hexagon -> 6 pentagons"}[N]
    print(f"\n=== N={N}  ({name}) ===")
    print("  c_i =", np.array2string(c, precision=5), "  sum =", f"{c.sum():.5f}")
    print(f"  I_{N} (direct)          = {I_N:+.8e}")
    print(f"  sum_i c_i I_{N-1}^(i)   = {rhs:+.8e}")
    print(f"  relative difference    = {rel:.2e}   ->  {'PASS' if rel < 2e-3 else 'FAIL'}")
    return rel < 2e-3


ok = all([verify(5), verify(6)])
print("\nALL PASS" if ok else "\nSOME FAILED")
