#!/usr/bin/env python3
"""Step 1 for dotted N>4: verify the FJT/Tarasov index-lowering recurrence generalizes
verbatim from the box (N=4, as in reduce_box) to the pentagon (N=5).

The recurrence (lower the largest power k; a = b - e_k, total = sum a):
    I(b) = 1/((b_k-1) det Y) [ (sum_i adj[k][i] (d - total - a[i])) I(a)
                               + sum_{i!=j} (-a[j] adj[k][i]) I(a + e_j - e_i) ]
with adj = det(Y) * Y^{-1} (adjugate), Y_ij = m_i^2+m_j^2-(r_i-r_j)^2.

It is EXACT in d, so we check LHS vs RHS at a generic d by direct Feynman-parameter
integration of every term (dotted, and pinched when a power hits 0). Euclidean kinematics.
If this holds, the box recurrence generalizes and dotted pentagons reduce (scalar leaves ->
the committed vNV code; pinched boxes -> reduce_box). Gitignored.
"""
import numpy as np
from scipy.special import gamma

np.random.seed(20260707)
D = 3.7
NS = 12_000_000

# same pentagon kinematics style as verify_pentagon_reduction.py
# larger offsets keep Y well away from rank-2 (det Y healthy); big masses keep F > 0
r = np.array([
    [0.00, 0.00, 0.00, 0.00],
    [1.20, 0.00, 0.00, 0.00],
    [0.50, 1.00, 0.00, 0.00],
    [0.00, 0.60, 1.30, 0.00],
    [0.70, 0.20, 0.40, 1.10],
])
msq = np.array([3.0, 3.5, 4.0, 4.5, 5.0])
N = 5


def gram2(i, j):
    dd = r[i] - r[j]
    return float(dd @ dd)


Y = np.array([[msq[i] + msq[j] - gram2(i, j) for j in range(N)] for i in range(N)])
detY = np.linalg.det(Y)
adj = detY * np.linalg.inv(Y)          # adjugate; adj[k][i] matches reduce_box's adj[k][i]


def direct_integral(powers):
    """I_n(powers) over the active lines (power>0), Feynman param via Dirichlet sampling.
       I_n = (-1)^A Gamma(A - d/2)/Gamma(A) * < F^{d/2 - A} >_{Dirichlet(a)},  A = sum a."""
    active = [i for i in range(N) if powers[i] > 0]
    a = np.array([powers[i] for i in active], dtype=float)
    A = int(round(a.sum()))
    x = np.random.dirichlet(a, size=NS)
    F = x @ msq[active]
    for p in range(len(active)):
        for q in range(p + 1, len(active)):
            F = F - x[:, p] * x[:, q] * gram2(active[p], active[q])
    return (-1) ** A * gamma(A - D / 2) / gamma(A) * np.mean(F ** (D / 2 - A))


def tarasov_rhs(b):
    """One index-lowering step, RHS evaluated by direct integration of each term."""
    k = int(np.argmax(b))              # first largest power (matches the Rust loop)
    a = list(b)
    a[k] -= 1
    total = sum(a)
    den = a[k] * detY                  # a[k] is now b[k]-1
    diag = sum(adj[k][i] * (D - total - a[i]) for i in range(N))
    rhs = (diag / den) * direct_integral(a)
    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            ch = list(a)
            ch[j] += 1
            ch[i] -= 1
            rhs += (-a[j] * adj[k][i] / den) * direct_integral(ch)
    return rhs


print(f"det Y = {detY:.4g},  d = {D}\n")
print(f"{'dotted config':18} {'LHS(direct)':>16} {'RHS(recurrence)':>16}  rel.err")
worst = 0.0
for b in ([2, 1, 1, 1, 1], [1, 2, 1, 1, 1], [2, 2, 1, 1, 1], [1, 1, 3, 1, 1]):
    lhs = direct_integral(b)
    rhs = tarasov_rhs(b)
    rel = abs(lhs - rhs) / abs(lhs)
    worst = max(worst, rel)
    print(f"{str(b):18} {lhs:+.9e} {rhs:+.9e}  {rel:.2e}")
print(f"\nworst rel.err = {worst:.2e}  ->  {'PASS' if worst < 3e-3 else 'FAIL'}")
