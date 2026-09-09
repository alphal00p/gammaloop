#!/usr/bin/env python3
"""Independent tensor one-loop integrator for numerator (dot(k,k))^p  (gitignored).

Convention-free (only uses r0=0): expands (k^2)^p = ((l - P)^2)^p after the Feynman shift
k = l - P, P = sum x_i r_i, and integrates the loop momentum l analytically via the standard
symmetric-moment rules, leaving a numeric Feynman-parameter (Dirichlet-MC) integral.

Euclidean spacelike setup, matching scipy_quad:  Delta = F = sum x_i m_i^2 + sum_{i<j} x_i x_j dist2.
  <l0^n0 l1^n1 l2^n2 l3^n3> / <1> = [prod (n_i-1)!!] / 2^m * Gamma(A-m-d/2)/Gamma(A-d/2) * Delta^m
  (m = sum n_i / 2, zero if any n_i odd), so at d=4:
  I[(k^2)^p] = (-1)^A / Gamma(A) * < sum_m coeff_m(P) * Gamma(A-m-2) * Delta^(2-A+m) >_{Dir(a)}.

Validated below against numerator=1 (scalar) and k^2 (= D0 + m0^2 decomposition).
"""
import numpy as np
import sympy as sp
from math import gamma, prod


def dfac(k):
    """Double factorial with the (-1)!! = 0!! = 1 convention."""
    r = 1
    while k > 1:
        r *= k
        k -= 2
    return r


def moment_buckets(p):
    """Symbolic: return {m: lambdified coeff_m(P0..P3)} for numerator (k^2)^p."""
    l = sp.symbols("l0 l1 l2 l3")
    P = sp.symbols("P0 P1 P2 P3")
    ksq = sum((l[i] - P[i]) ** 2 for i in range(4))
    poly = sp.Poly(sp.expand(ksq ** p), *l)
    buckets = {}
    for monom, coeff in poly.terms():
        if any(n % 2 for n in monom):
            continue
        m = sum(monom) // 2
        fac = 1
        for n in monom:
            fac *= dfac(n - 1)                      # (n-1)!!,  (-1)!! = 1
        buckets[m] = buckets.get(m, 0) + coeff * sp.Rational(fac, 2 ** m)
    return {m: sp.lambdify(P, sp.expand(e), "numpy") for m, e in buckets.items()}


def tensor_direct(offsets, masses, exps, p, ns=4_000_000, seed=20260718):
    """I[(k^2)^p] via Feynman-parameter Dirichlet MC.  Returns (value, mc_err)."""
    rng = np.random.default_rng(seed)
    a = np.array(exps, dtype=float)
    A = int(sum(exps))
    n = len(exps)
    r = np.array([o[:4] for o in offsets])         # (n,4) offsets  (r0 = 0)
    msq = np.array(masses)
    buckets = moment_buckets(p)

    x = rng.dirichlet(a, size=ns)                  # (ns,n)
    P = x @ r                                       # (ns,4)  = sum x_i r_i
    Delta = x @ msq
    for i in range(n):
        for j in range(i + 1, n):
            dij = r[i] - r[j]
            Delta = Delta + x[:, i] * x[:, j] * float(dij @ dij)
    integrand = np.zeros(ns)
    for m, fm in buckets.items():
        cm = np.broadcast_to(fm(P[:, 0], P[:, 1], P[:, 2], P[:, 3]), (ns,)).astype(float)
        integrand = integrand + cm * gamma(A - m - 2.0) * Delta ** (2 - A + m)
    # dot(k,k)_Minkowski = -k^2_Euclidean, so the physical numerator (k^2)^p carries (-1)^p
    pref = (-1) ** A * (-1) ** p / gamma(A)
    return pref * np.mean(integrand), abs(pref) * np.std(integrand) / np.sqrt(ns)


def q_vectors(offsets):
    """oneloop's external momenta q_j = r_j - r_{j-1}  (leg differences)."""
    r = np.array([o[:4] for o in offsets])
    return [r[j] - r[j - 1] for j in range(1, len(offsets))]


def tensor_direct_general(offsets, masses, exps, num_builder, ns=4_000_000, seed=20260718):
    """General numerator (polynomial in dot(k,k) and dot(k,q_j)).  num_builder(l,P,qv,LL,LQ)
    returns a sympy expr in the shifted loop momentum l (dot(k,k)_Mink=-k_E^2 via LL)."""
    l = sp.symbols("l0 l1 l2 l3")
    Psy = sp.symbols("P0 P1 P2 P3")
    qv = q_vectors(offsets)

    def LL():                                             # dot(k,k)_Mink = -(l-P)^2_Euclid
        return -sum((l[i] - Psy[i]) ** 2 for i in range(4))

    def LQ(q, sgn):                                       # dot(k,q)_Mink = sgn*(l-P).q_Euclid
        return sgn * sum((l[i] - Psy[i]) * float(q[i]) for i in range(4))

    poly = sp.Poly(sp.expand(num_builder(qv, LL, LQ)), *l)
    buckets = {}
    for monom, coeff in poly.terms():
        if any(nn % 2 for nn in monom):
            continue
        m = sum(monom) // 2
        fac = 1
        for nn in monom:
            fac *= dfac(nn - 1)
        buckets[m] = buckets.get(m, 0) + coeff * sp.Rational(fac, 2 ** m)
    funcs = {m: sp.lambdify(Psy, sp.expand(e), "numpy") for m, e in buckets.items()}

    rng = np.random.default_rng(seed)
    a = np.array(exps, dtype=float)
    A = int(sum(exps))
    n = len(exps)
    r = np.array([o[:4] for o in offsets])
    msq = np.array(masses)
    x = rng.dirichlet(a, size=ns)
    P = x @ r
    Delta = x @ msq
    for i in range(n):
        for j in range(i + 1, n):
            dij = r[i] - r[j]
            Delta = Delta + x[:, i] * x[:, j] * float(dij @ dij)
    integrand = np.zeros(ns)
    for m, fm in funcs.items():
        cm = np.broadcast_to(fm(P[:, 0], P[:, 1], P[:, 2], P[:, 3]), (ns,)).astype(float)
        integrand = integrand + cm * gamma(A - m - 2.0) * Delta ** (2 - A + m)
    pref = (-1) ** A / gamma(A)
    return pref * np.mean(integrand), abs(pref) * np.std(integrand) / np.sqrt(ns)


def tensor_covariant(offsets, masses, exps, num, ns=4_000_000, seed=20260718):
    """Metric-correct tensor integral via covariant Wick contraction.  The loop-momentum
    contractions use the MINKOWSKI metric: for external vectors G(u,v) = dot(u,v)_Mink =
    -(u.v)_Euclid (from the invariants), and the trace g_mu^mu = d.  Moment factor
    c_m = (-Delta/2)^m Gamma(A-m-d/2)/Gamma(A-d/2).  Handles the specific benchmark numerators."""
    rng = np.random.default_rng(seed)
    a = np.array(exps, dtype=float)
    A = int(sum(exps))
    n = len(exps)
    r = np.array([o[:4] for o in offsets])
    msq = np.array(masses)
    qv = [r[j] - r[j - 1] for j in range(1, n)]        # q_j = r_j - r_{j-1}
    d = 4.0

    x = rng.dirichlet(a, size=ns)
    P = x @ r                                           # (ns,4) Euclidean
    Delta = x @ msq
    for i in range(n):
        for j in range(i + 1, n):
            dij = r[i] - r[j]
            Delta = Delta + x[:, i] * x[:, j] * float(dij @ dij)
    c1 = -Delta / (2 * (A - 1 - d / 2))                 # c_m = (-D/2)^m / prod_{k=1}^m (A-k-d/2)
    c2 = Delta ** 2 / (4 * (A - 1 - d / 2) * (A - 2 - d / 2))

    def GP(u):                                          # G(P,u)_Mink = -(P.u)_Euclid, per x
        return -(P @ np.asarray(u, dtype=float))

    def GE(u, v):                                       # G(u,v)_Mink scalar
        return -float(np.asarray(u) @ np.asarray(v))

    GPP = -np.sum(P * P, axis=1)                        # G(P,P)_Mink

    import re
    if num == "ll":                                     # dot(k,k)
        N = d * c1 + GPP
    elif num == "ll2":                                  # (dot(k,k))^2
        N = d * (d + 2) * c2 + GPP * c1 * (4 + 2 * d) + GPP ** 2
    elif re.fullmatch(r"q\d", num):                     # dot(k, q_i)
        i = int(num[1]) - 1
        N = -GP(qv[i])
    elif re.fullmatch(r"q\dq\d", num):                  # dot(k,q_i) dot(k,q_j)
        i, j = int(num[1]) - 1, int(num[3]) - 1
        N = GE(qv[i], qv[j]) * c1 + GP(qv[i]) * GP(qv[j])
    elif re.fullmatch(r"llq\d", num):                   # dot(k,k) dot(k,q_i)
        i = int(num[3]) - 1
        N = -c1 * GP(qv[i]) * (d + 2) - GPP * GP(qv[i])
    else:
        raise ValueError(num)

    pref = (-1) ** A * gamma(A - 2.0) / gamma(A)
    integrand = Delta ** (2 - A) * N
    return pref * np.mean(integrand), abs(pref) * np.std(integrand) / np.sqrt(ns)


def tensor_laurent_llp(offsets, masses, exps, p, ns=8_000_000, seed=20260719):
    """Full Laurent (pole 1/eps, finite) of the (k^2)^p tensor integral, d=4-2eps, in the
    OneLOop MS-bar scheme.  Works for finite and single-log-divergent cases.  Uses the
    rank-decomposition from moment_buckets(p): I(eps) = (-1)^{A+p}/Gamma(A) sum_m
    Gamma(A-m-2+eps) E[coeff_m(P) Delta^{2-A+m-eps}], expanded to O(eps^0)."""
    rng = np.random.default_rng(seed)
    a = np.array(exps, dtype=float)
    A = int(sum(exps))
    n = len(exps)
    r = np.array([o[:4] for o in offsets])
    msq = np.array(masses)
    buckets = moment_buckets(p)
    x = rng.dirichlet(a, size=ns)
    P = x @ r
    Delta = x @ msq
    for i in range(n):
        for j in range(i + 1, n):
            dij = r[i] - r[j]
            Delta = Delta + x[:, i] * x[:, j] * float(dij @ dij)
    L = np.log(Delta)

    # <(k^2)^p> = sum_s C(p,s) [prod_{l=s}^{p-1}(d+2l)] c_{p-s} w^s,  w = G(P,P) = -(P.P)_E,
    # c_m = (-Delta/2)^m Gamma(A-m-d/2)/Gamma(A-d/2).  d = 4 - 2 eps kept symbolic (the
    # trace factors carry the eps-dependence that a d=4 moment reduction would drop).
    w = -np.sum(P * P, axis=1)
    e = sp.Symbol("eps")
    d = 4 - 2 * e
    from math import comb
    expr = sp.Integer(0)
    for s in range(p + 1):
        trace = sp.Integer(1)
        for l in range(s, p):
            trace *= (d + 2 * l)
        m = p - s
        base = 2 - A + m                                       # power of Delta after folding
        mS = float(np.mean(w ** s * Delta ** base))
        mSL = float(np.mean(w ** s * Delta ** base * L))
        expr += (comb(p, s) * trace * sp.Rational((-1) ** m, 2 ** m)
                 * sp.gamma(A - m - 2 + e) * (mS - e * mSL))
    expr *= sp.Integer((-1) ** A) / sp.gamma(A)
    ser = sp.series(expr, e, 0, 1).removeO()
    pole = complex(ser.coeff(e, -1)).real
    finite = complex(ser.coeff(e, 0)).real
    finite += np.euler_gamma * pole                            # bare Gamma(eps) -> MS-bar
    return pole, finite


# ---- geometry (config A / B), matching emit_reductions.rs -------------------
OFF_A = [[0, 0, 0, 0], [.30, 0, 0, 0], [.10, .35, 0, 0], [0, .20, .40, 0],
         [.05, 0, .15, .45], [.22, .12, 0, .18], [.08, .25, .30, .10]]
MSQ_A = [1.0, 1.21, 0.81, 1.44, 0.80, 1.05, 0.95]


def scalar_quad(offsets, masses, exps):
    """Deterministic scalar reference on <=4 lines (reused from crosscheck idea)."""
    from crosscheck import scipy_quad
    return float(scipy_quad(offsets, masses, exps))


if __name__ == "__main__":
    # --- validation 1: p=0 tensor == scalar (triangle, box) ---
    for n in (3, 4):
        off, m = OFF_A[:n], MSQ_A[:n]
        v0, e0 = tensor_direct(off, m, [1] * n, 0)
        vs = scalar_quad(off, m, [1] * n)
        print(f"p=0 N={n}: tensor {v0:+.6e} vs scalar-quad {vs:+.6e}  rel {abs(v0-vs)/abs(vs):.1e}")

    # --- validation 2: p=1 (k^2) == I(drop line0) + m0^2 I(scalar) [finite: use box] ---
    for n in (4, 5):
        off, m = OFF_A[:n], MSQ_A[:n]
        v1, e1 = tensor_direct(off, m, [1] * n, 1)
        drop = scalar_quad(off[1:], m[1:], [1] * (n - 1)) if n - 1 <= 4 else None
        full = scalar_quad(off, m, [1] * n) if n <= 4 else None
        if drop is not None and full is not None:
            true = drop + m[0] * full
            print(f"p=1 N={n}: tensor {v1:+.6e} vs D0+m0^2 decomp {true:+.6e}  rel {abs(v1-true)/abs(true):.1e}")
        else:
            print(f"p=1 N={n}: tensor {v1:+.6e} (+-{e1:.1e})  [decomp needs 5pt; skipped]")

    # --- (k^2)^2 finite on pentagon (A=5) and hexagon (A=6) ---
    for n in (5, 6):
        off, m = OFF_A[:n], MSQ_A[:n]
        v2, e2 = tensor_direct(off, m, [1] * n, 2)
        print(f"(k^2)^2 N={n}: tensor {v2:+.8e}  (+-{e2:.1e})")
