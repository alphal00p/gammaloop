"""Build the transverse-projected gauge-invariant H->gamma gamma W-loop numerator (see docs/10)."""
import numpy as np

RANK = 6
DVALS = list(range(4, 12))


def monomials():
    out = []
    for a in range(RANK//2 + 1):
        for b in range(RANK - 2*a + 1):
            for c in range(RANK - 2*a - b + 1):
                out.append((a, b, c))
    return out


MONS = monomials()


def _vgamma(kph, k2, k3, g):
    d = len(kph); V = np.zeros((d, d, d))
    A, B, C = k2 - k3, k3 - kph, kph - k2
    for mu in range(d):
        V[mu] += g * A[mu] + np.outer(B, g[:, mu]) + np.outer(g[mu, :], C)
    return V


def _prop(p, g, mW):
    return -(g - np.outer(p, p)/mW**2)


def _triangle(qf, qs, g, mW, k):
    pa, pb, pc = k, k+qf, k+qf+qs
    Da, Db, Dc = _prop(pa, g, mW), _prop(pb, g, mW), _prop(pc, g, mW)
    Vf, Vs = _vgamma(qf, pa, -pb, g), _vgamma(qs, pb, -pc, g)
    Hww = mW**2 * g; d = len(g); N = np.zeros((d, d))
    for i in range(d):
        for j in range(d):
            N[i, j] = np.trace(Da @ g @ Vf[i] @ g @ Db @ g @ Vs[j] @ g @ Dc @ g @ Hww @ g)
    return N


def _seagull(g, mW, k, q1, q2):
    Da, Db, Hww = _prop(k, g, mW), _prop(k+q1+q2, g, mW), mW**2 * g
    d = len(g); N = np.zeros((d, d))
    for mu in range(d):
        for nu in range(d):
            Vq = 2*g[mu, nu]*g - np.outer(g[mu], g[nu]) - np.outer(g[nu], g[mu])
            N[mu, nu] = np.trace(Da @ g @ Hww @ g @ Db @ g @ Vq @ g)
    return N


def _projected(d, k, q1, q2, mW, which):
    g = np.diag([1.0] + [-1.0]*(d-1))
    N = (_triangle(q1, q2, g, mW, k) if which == "N1" else
         _triangle(q2, q1, g, mW, k) if which == "N2" else
         _seagull(g, mW, k, q1, q2))
    q1q2 = q1 @ g @ q2; ql1, ql2 = g @ q1, g @ q2
    P = g.copy()
    for m in range(d):
        for n in range(d):
            P[m, n] -= (ql1[m]*ql2[n] + ql2[m]*ql1[n]) / q1q2
    return np.einsum('mn,mn->', P, N)


def _fit_diagram(d, which, s, mW, nsamp=140, seed=11):
    E = np.sqrt(s)/2.0
    q1 = np.zeros(d); q2 = np.zeros(d)
    q1[0] = E; q1[3] = E; q2[0] = E; q2[3] = -E
    g = np.diag([1.0] + [-1.0]*(d-1))
    rng = np.random.default_rng(seed + d)
    A = np.zeros((nsamp, len(MONS))); b = np.zeros(nsamp)
    for si in range(nsamp):
        k = rng.standard_normal(d) * rng.uniform(0.4, 1.6)
        ll, lq1, lq2 = k @ g @ k, k @ g @ q1, k @ g @ q2
        A[si] = [ll**a * lq1**bb * lq2**cc for (a, bb, cc) in MONS]
        b[si] = _projected(d, k, q1, q2, mW, which)
    coef, *_ = np.linalg.lstsq(A, b, rcond=None)
    return {MONS[i]: coef[i] for i in range(len(MONS))}, np.abs(A @ coef - b).max()


def build(tau, mW=1.0, seed=11):
    """Return (combined {mon: d-poly coeffs}, s, max_fit_residual)."""
    s = 4.0*tau
    per = {w: {m: [] for m in MONS} for w in ("N1", "N2", "seag")}
    maxr = 0.0
    for d in DVALS:
        for w in per:
            c, r = _fit_diagram(d, w, s, mW, seed=seed)
            maxr = max(maxr, r)
            for m in MONS:
                per[w][m].append(c[m])
    V = np.vander(np.array(DVALS, float), 8, increasing=True)
    poly = {w: {m: np.linalg.lstsq(V, np.array(per[w][m]), rcond=None)[0] for m in MONS}
            for w in per}
    from collections import defaultdict
    comb = defaultdict(lambda: np.zeros(8))
    for m in MONS:
        comb[m] = comb[m] + poly["N1"][m]
    for (a, b, c) in MONS:
        comb[(a, c, b)] = comb[(a, c, b)] + poly["N2"][(a, b, c)]      # crossed: q1<->q2
    for (a, b, c) in MONS:                                            # seagull * (ll+2 lq1 - mW^2)
        p = poly["seag"][(a, b, c)]
        comb[(a+1, b, c)] = comb[(a+1, b, c)] + p
        comb[(a, b+1, c)] = comb[(a, b+1, c)] + 2.0*p
        comb[(a, b, c)] = comb[(a, b, c)] - (mW**2)*p
    # keep physical monomials (rank <= RANK); anything above / tiny is fit noise.
    return ({m: v for m, v in comb.items()
             if np.abs(v).max() > 1e-4 and 2*m[0] + m[1] + m[2] <= RANK}, s, maxr)
