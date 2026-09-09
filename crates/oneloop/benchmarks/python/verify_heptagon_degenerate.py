#!/usr/bin/env python3
"""Verify the degenerate-Cayley reduction step used for N >= 7 (det(Y) = 0 in 4D).

For a heptagon the modified Cayley matrix Y (7x7) has det(Y) = 0 identically
(rank <= 6), so the standard vNV c_i = cofactor/det(Y) diverges.  The degenerate
reduction instead takes a null vector (c_1..c_7, K) of the BORDERED Cayley
[[Y, 1],[1^T, 0]] with K != 0.  Then sum_i c_i D_i = K is CONSTANT in the loop
momentum, so

    1/(D_1..D_7) = (1/K) sum_i c_i / prod_{j!=i} D_j

is an exact integrand identity (partial-fractioning), not an integrated relation.
We check it to machine precision at random loop momenta -- this validates the exact
c_i and K that the Rust engine's `degenerate_coeffs` computes from the same matrix.
Gitignored.
"""
import numpy as np

np.random.seed(20260711)
G = np.diag([1.0, -1.0, -1.0, -1.0])  # Minkowski


def dot(a, b):
    return float(a @ G @ b)


def run_trial(n):
    r = [np.zeros(4)] + [np.random.randn(4) for _ in range(n - 1)]  # r_1 = 0, n lines
    msq = [np.random.rand() + 0.5 for _ in range(n)]
    Y = np.array([[msq[i] + msq[j] - dot(r[i] - r[j], r[i] - r[j])
                   for j in range(n)] for i in range(n)])
    assert abs(np.linalg.det(Y)) < 1e-6, f"det(Y_{n}) should vanish (N>=7 in 4D)"

    # Null vector (c, K) of the bordered Cayley with K != 0.
    yb = np.zeros((n + 1, n + 1))
    yb[:n, :n] = Y
    yb[:n, n] = 1.0
    yb[n, :n] = 1.0
    _, s, vt = np.linalg.svd(yb)
    null = [vt[i] for i in range(n + 1) if s[i] < 1e-9]
    vec = next((v for v in null if abs(v[n]) > 1e-9), None)
    if vec is None:
        vec = sum(null)
    c = vec[:n] / vec[n]  # c_i / K

    # Check the integrand identity 1/prod D = (1/K) sum_i c_i / prod_{j!=i} D_j.
    worst = 0.0
    for _ in range(50):
        l = np.random.randn(4)
        D = [dot(l + r[i], l + r[i]) - msq[i] for i in range(n)]
        lhs = 1.0 / np.prod(D)
        rhs = sum(c[i] * np.prod([D[j] for j in range(n) if j != i]) ** -1
                  for i in range(n))
        worst = max(worst, abs(lhs - rhs) / (abs(lhs) + 1e-30))
    return worst


for n in (7, 8, 9):
    worst = max(run_trial(n) for _ in range(100))
    print(f"N={n}: worst rel. error of the integrand identity = {worst:.2e}  ->",
          "PASS" if worst < 1e-8 else "FAIL")


def run_dotted(a):
    # Gap 1: dotted N=7 power-lowering  1/prod D_i^{a_i} = sum_l (c_l/K) / (D_l^{a_l-1} prod_{j!=l} D_j^{a_j}).
    n = len(a)
    r = [np.zeros(4)] + [np.random.randn(4) for _ in range(n - 1)]
    msq = [np.random.rand() + 0.5 for _ in range(n)]
    Y = np.array([[msq[i] + msq[j] - dot(r[i] - r[j], r[i] - r[j])
                   for j in range(n)] for i in range(n)])
    yb = np.zeros((n + 1, n + 1))
    yb[:n, :n] = Y
    yb[:n, n] = 1.0
    yb[n, :n] = 1.0
    _, s, vt = np.linalg.svd(yb)
    null = [vt[i] for i in range(n + 1) if s[i] < 1e-9]
    vec = next((v for v in null if abs(v[n]) > 1e-9), None)
    c = vec[:n] / vec[n]
    worst = 0.0
    for _ in range(50):
        l = np.random.randn(4)
        D = [dot(l + r[i], l + r[i]) - msq[i] for i in range(n)]
        lhs = 1.0 / np.prod([D[i] ** a[i] for i in range(n)])
        rhs = 0.0
        for x in range(n):
            b = list(a)
            b[x] -= 1
            rhs += c[x] / np.prod([D[i] ** b[i] for i in range(n)])
        worst = max(worst, abs(lhs - rhs) / (abs(lhs) + 1e-30))
    return worst


for a in ([2, 1, 1, 1, 1, 1, 1], [2, 2, 1, 1, 1, 1, 1], [3, 1, 1, 1, 1, 1, 1]):
    worst = max(run_dotted(a) for _ in range(50))
    print(f"dotted a={a}: worst rel. error of the power-lowering identity = {worst:.2e}  ->",
          "PASS" if worst < 1e-8 else "FAIL")
