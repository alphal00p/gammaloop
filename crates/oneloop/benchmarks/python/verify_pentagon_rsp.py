#!/usr/bin/env python3
"""Derive + verify the RSP substitution rules for a PENTAGON numerator (N=5).

A pentagon has 5 propagators D_i = (l + r_i)^2 - m_i^2 with offsets
    r_1 = 0, r_2 = q1, r_3 = q1+q2, r_4 = q1+q2+q3, r_5 = q1+q2+q3+q4,
so there are 4 external momenta q1..q4 (n_ext = 4). The reducible scalar products
are {l.l, l.q1, l.q2, l.q3, l.q4} -- exactly 5, matching D_1..D_5. So the map is
BIJECTIVE: every external direction is reducible and there is NO irreducible scalar
product at the top pentagon (unlike the box, whose sub-topologies do have ISPs).

Claimed substitution rules (l.l from r_1 = 0; l.q_i from D_{i+1} - D_i):
    l.l   = D_1 + m_1^2
    l.q_i = [ D_{i+1} - D_i - (s_{1,i+1} - s_{1,i}) + m_{i+1}^2 - m_i^2 ] / 2,
where s_{1,j} = (r_1 - r_j)^2 = r_j^2 (since r_1 = 0). For i = 1,2,3 these reproduce
the box rules already in the crate; i = 4 is the new pentagon rule:
    l.q4 = [ D_5 - D_4 - (s_{1,5} - s_{1,4}) + m_5^2 - m_4^2 ] / 2.

These are exact algebraic identities in the loop momentum (no integration needed),
so we check them to machine precision at random Euclidean kinematics. Gitignored.
"""
import numpy as np

np.random.seed(20260710)
DIM = 6  # ambient vector dimension (>= 4 so q1..q4 are generic)


def run_trial():
    l = np.random.randn(DIM)
    q = [np.random.randn(DIM) for _ in range(4)]        # q1..q4
    msq = np.random.rand(5) + 0.5                        # m_1^2 .. m_5^2
    r = [np.zeros(DIM)]                                  # r_1 = 0
    for a in range(4):
        r.append(r[-1] + q[a])                           # r_2..r_5
    D = [float((l + r[i]) @ (l + r[i]) - msq[i]) for i in range(5)]
    s1 = [float((r[0] - r[j]) @ (r[0] - r[j])) for j in range(5)]  # s_{1,j} = r_j^2

    errs = []
    # l.l = D_1 + m_1^2
    errs.append(abs(float(l @ l) - (D[0] + msq[0])))
    # l.q_i = [D_{i+1} - D_i - (s_{1,i+1}-s_{1,i}) + m_{i+1}^2 - m_i^2]/2
    for i in range(1, 5):                                # i = 1..4  (q index i)
        lhs = float(l @ q[i - 1])
        rhs = (D[i] - D[i - 1] - (s1[i] - s1[i - 1]) + msq[i] - msq[i - 1]) / 2
        errs.append(abs(lhs - rhs))

    # Bijectivity: the 5x5 Jacobian d(RSP)/d(D) must be invertible (det != 0).
    # RSP = [l.l, l.q1, l.q2, l.q3, l.q4] is affine-linear in D with this Jacobian:
    J = np.zeros((5, 5))
    J[0, 0] = 1.0                                        # l.l = D1 + const
    for i in range(1, 5):
        J[i, i] = 0.5
        J[i, i - 1] = -0.5
    return max(errs), abs(np.linalg.det(J))


worst, dets = 0.0, []
for _ in range(2000):
    e, det = run_trial()
    worst = max(worst, e)
    dets.append(det)

print(f"worst RSP-identity error over 2000 random pentagons = {worst:.2e}")
print(f"Jacobian det(d RSP / d D) = {dets[0]:.4f}  (nonzero => bijective, no ISP)")
print("->", "PASS" if worst < 1e-9 and min(dets) > 1e-9 else "FAIL")
