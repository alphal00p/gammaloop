"""Finite-eps retarded spectral integrand with a deterministic radial quadrature.

Loop momenta l_i = lam * rho_i * nhat_i.  For every outer point (rho, nhat) the radial
integral over lam in [0, inf) is done with Gauss-Legendre panels that are sinh-mapped
around every root of every threshold function h(lam) = sum_j h_j E_j(lam) = E (the
locations of the eps-broadened poles of R(E+i eps)) and around near-tangent extrema
of h.  The outer integral is left to vegas.
"""

from __future__ import annotations

import numpy as np

# isort: split
from matsubara import Compiled, build_terms, get_graph


def gauss(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return x, w


class Problem:
    def __init__(
        self,
        graph: str,
        energy: float,
        momentum: float,
        temperature: float,
        eps_list,
        ng: int = 10,
        nsub: int = 6,
        scale: float = 1.0,
        grid: int = 400,
        width_factor: float = 0.5,
        extremum_window: float = 0.15,
        fit_weights=None,
    ):
        self.gname = graph
        self.g = get_graph(graph)
        self.C = Compiled(self.g, build_terms(self.g))
        self.E = float(energy)
        self.k = float(momentum)
        self.T = float(temperature)
        self.eps = np.array(eps_list, dtype=float)
        self.z = self.E + 1j * self.eps
        self.ng, self.nsub = ng, nsub
        self.s = scale
        self.L = self.g.n_loops
        self.N = self.g.n_lines
        self.Hm = np.array(self.C.thresholds, dtype=float)  # (nh, N)
        self.lc = np.array(self.g.loop_coeff, dtype=float)  # (N, L)
        self.ec = np.array(self.g.ext_coeff, dtype=float)  # (N,)
        self.m2 = np.array(self.g.masses) ** 2
        self.Cc = self.ec**2 * self.k**2 + self.m2  # (N,)
        self.general = self.k != 0.0
        self.dim = max(3 * self.L - 2 if self.general else 3 * self.L - 4, 1)
        xg = np.linspace(0.0, 0.998, grid)
        self.lam_grid = self.s * xg / (1.0 - xg)
        self.width_factor = width_factor
        self.ext_window = extremum_window
        self.gx, self.gw = gauss(ng)
        # first component: fixed linear combination of the eps ladder (fit at integrand level)
        self.fit_weights = (
            None if fit_weights is None else np.asarray(fit_weights, float)
        )
        self.norm = self.g.kappa / self.E / (2 * np.pi) ** (3 * self.L)

    # ------------------------------------------------------------------ outer map
    def outer(self, u):
        """u (B, dim) -> v (B, N, 3) line direction vectors (lam=1 part), jac (B,)."""
        B = u.shape[0]
        L = self.L
        col = 0
        if L == 1:
            rho = np.ones((B, 1))
            jac = np.ones(B)
        elif L == 2:
            a = 0.5 * np.pi * u[:, 0]
            rho = np.stack([np.cos(a), np.sin(a)], axis=1)
            jac = 0.5 * np.pi * np.ones(B)
            col = 1
        else:
            a1 = 0.5 * np.pi * u[:, 0]
            a2 = 0.5 * np.pi * u[:, 1]
            rho = np.stack(
                [np.cos(a1), np.sin(a1) * np.cos(a2), np.sin(a1) * np.sin(a2)], axis=1
            )
            jac = (0.5 * np.pi) ** 2 * np.sin(a1)
            col = 2
        jac = jac * np.prod(rho**2, axis=1)
        nhat = np.zeros((B, L, 3))
        if self.general:
            ct = 2 * u[:, col] - 1
            col += 1
            st = np.sqrt(np.clip(1 - ct * ct, 0, None))
            nhat[:, 0] = np.stack([st, np.zeros(B), ct], axis=1)
            jac = jac * 2 * 2 * np.pi
            start = 1
        else:
            nhat[:, 0, 2] = 1.0
            jac = jac * 4 * np.pi
            if L >= 2:
                ct = 2 * u[:, col] - 1
                col += 1
                st = np.sqrt(np.clip(1 - ct * ct, 0, None))
                nhat[:, 1] = np.stack([st, np.zeros(B), ct], axis=1)
                jac = jac * 2 * 2 * np.pi
            start = 2
        for i in range(start, L):
            ct = 2 * u[:, col] - 1
            ph = 2 * np.pi * u[:, col + 1]
            col += 2
            st = np.sqrt(np.clip(1 - ct * ct, 0, None))
            nhat[:, i] = np.stack([st * np.cos(ph), st * np.sin(ph), ct], axis=1)
            jac = jac * 2 * 2 * np.pi
        assert col == self.dim or (L == 1 and not self.general and col == 0)
        loopdir = rho[:, :, None] * nhat  # (B, L, 3)
        v = np.einsum("nl,bld->bnd", self.lc, loopdir)  # (B, N, 3)
        return v, jac

    def quad_coeffs(self, v):
        A = np.sum(v * v, axis=2)  # (B, N)
        Bc = self.ec[None, :] * v[:, :, 2] * self.k  # k along z
        return A, Bc

    def energies(self, A, Bc, lam):
        """A, Bc (B, N); lam (B, P) -> E (B, P, N), dE (B, P, N)."""
        q = (
            A[:, None, :] * lam[:, :, None] ** 2
            + 2 * Bc[:, None, :] * lam[:, :, None]
            + self.Cc
        )
        E = np.sqrt(q)
        dE = (A[:, None, :] * lam[:, :, None] + Bc[:, None, :]) / E
        return E, dE

    # ------------------------------------------------------------------ centres
    def _h_point(self, A, Bc, hrow, lam):
        """Vectorised h, h', h'' for bracket lists: A,Bc (M,N), hrow (M,N), lam (M,)."""
        q = A * lam[:, None] ** 2 + 2 * Bc * lam[:, None] + self.Cc
        E = np.sqrt(q)
        dE = (A * lam[:, None] + Bc) / E
        d2E = (A * self.Cc - Bc**2) / E**3
        return (hrow * E).sum(1), (hrow * dE).sum(1), (hrow * d2E).sum(1)

    def centres(self, A, Bc, eps_min):
        B = A.shape[0]
        lg = np.broadcast_to(self.lam_grid, (B, self.lam_grid.size))
        E, dE = self.energies(A, Bc, lg)  # (B, G, N)
        h = E @ self.Hm.T - self.E  # (B, G, nh)
        hp = dE @ self.Hm.T
        out_b, out_lam, out_w = [], [], []
        # ---- simple roots
        sc = np.signbit(h[:, :-1, :]) != np.signbit(h[:, 1:, :])
        bi, gi, hi = np.nonzero(sc)
        if bi.size:
            lo = self.lam_grid[gi].copy()
            hi_l = self.lam_grid[gi + 1].copy()
            Ab, Bb, Hr = A[bi], Bc[bi], self.Hm[hi]
            flo = self._h_point(Ab, Bb, Hr, lo)[0] - self.E
            for _ in range(60):
                mid = 0.5 * (lo + hi_l)
                fm = self._h_point(Ab, Bb, Hr, mid)[0] - self.E
                same = np.signbit(fm) == np.signbit(flo)
                lo = np.where(same, mid, lo)
                flo = np.where(same, fm, flo)
                hi_l = np.where(same, hi_l, mid)
            lr = 0.5 * (lo + hi_l)
            _, d1, d2 = self._h_point(Ab, Bb, Hr, lr)
            w = eps_min / np.maximum(np.abs(d1), 1e-300)
            wt = np.sqrt(2 * eps_min / np.maximum(np.abs(d2), 1e-300))
            w = np.minimum(w, wt) * self.width_factor
            out_b.append(bi)
            out_lam.append(lr)
            out_w.append(w)
        # ---- extrema close to E (near tangencies)
        sc = np.signbit(hp[:, :-1, :]) != np.signbit(hp[:, 1:, :])
        bi, gi, hi = np.nonzero(sc)
        if bi.size:
            lo = self.lam_grid[gi].copy()
            hi_l = self.lam_grid[gi + 1].copy()
            Ab, Bb, Hr = A[bi], Bc[bi], self.Hm[hi]
            flo = self._h_point(Ab, Bb, Hr, lo)[1]
            for _ in range(60):
                mid = 0.5 * (lo + hi_l)
                fm = self._h_point(Ab, Bb, Hr, mid)[1]
                same = np.signbit(fm) == np.signbit(flo)
                lo = np.where(same, mid, lo)
                flo = np.where(same, fm, flo)
                hi_l = np.where(same, hi_l, mid)
            le = 0.5 * (lo + hi_l)
            hv, _, d2 = self._h_point(Ab, Bb, Hr, le)
            delta = np.abs(hv - self.E)
            keep = delta < self.ext_window
            if np.any(keep):
                w = np.sqrt(
                    2 * np.maximum(delta, eps_min) / np.maximum(np.abs(d2), 1e-300)
                )
                w = w * self.width_factor
                out_b.append(bi[keep])
                out_lam.append(le[keep])
                out_w.append(w[keep])
        if out_b:
            b = np.concatenate(out_b)
            lam = np.concatenate(out_lam)
            w = np.concatenate(out_w)
        else:
            b = np.zeros(0, int)
            lam = np.zeros(0)
            w = np.zeros(0)
        # x-space
        x = lam / (self.s + lam)
        wx = w * self.s / (self.s + lam) ** 2
        wx = np.clip(wx, 1e-13, 1.0)
        counts = np.bincount(b, minlength=B)
        nmax = max(int(counts.max()) if counts.size else 0, 1)
        X = np.full((B, nmax), np.nan)
        W = np.full((B, nmax), np.nan)
        order = np.argsort(b, kind="stable")
        b, x, wx = b[order], x[order], wx[order]
        starts = np.concatenate([[0], np.cumsum(counts)[:-1]])
        slot = np.arange(b.size) - starts[b]
        X[b, slot] = x
        W[b, slot] = wx
        # pad: rows with no centre get a broad centre; others duplicate the first centre
        empty = counts == 0
        X[empty, 0] = 0.5
        W[empty, 0] = 1.0
        firstx = X[:, 0:1].copy()
        firstw = W[:, 0:1].copy()
        pad = np.isnan(X)
        X = np.where(pad, firstx, X)
        W = np.where(pad, firstw, W)
        srt = np.argsort(X, axis=1)
        X = np.take_along_axis(X, srt, 1)
        W = np.take_along_axis(W, srt, 1)
        return X, W, counts

    def nodes(self, X, W, xmax=1.0 - 1e-7):
        B = X.shape[0]
        mids = np.concatenate(
            [np.zeros((B, 1)), 0.5 * (X[:, 1:] + X[:, :-1]), np.full((B, 1), xmax)],
            axis=1,
        )
        tlo = np.arcsinh((mids[:, :-1] - X) / W)  # (B, n)
        thi = np.arcsinh((mids[:, 1:] - X) / W)
        # nsub sub-intervals, ng points each
        frac = (
            np.arange(self.nsub)[:, None] + 0.5 * (self.gx[None, :] + 1)
        ) / self.nsub  # (nsub, ng)
        frac = frac.reshape(-1)
        dt = (thi - tlo)[:, :, None]
        t = tlo[:, :, None] + dt * frac[None, None, :]
        wt = np.tile(self.gw, self.nsub)[None, None, :] * dt / (2 * self.nsub)
        x = X[:, :, None] + W[:, :, None] * np.sinh(t)
        wx = wt * W[:, :, None] * np.cosh(t)
        return x.reshape(B, -1), wx.reshape(B, -1)

    # ------------------------------------------------------------------ integrand
    def radial(self, v, x, wx):
        """Radial integral at fixed outer point: returns (B, neps) of int dlam lam^(3L-1) Im R."""
        B, P = x.shape
        A, Bc = self.quad_coeffs(v)
        lam = self.s * x / (1.0 - x)
        jac = wx * self.s / (1.0 - x) ** 2 * lam ** (3 * self.L - 1)
        E, _ = self.energies(A, Bc, lam)
        Ef = E.reshape(B * P, self.N)
        vals = self.evaluate_guarded(Ef, A, Bc, lam)
        vals = vals.reshape(B, P, -1)
        # the lam^(3L-1) weight makes lam < 1e-6 irrelevant; there all momenta vanish
        # together and spurious denominators are exactly degenerate at k=0.
        vals[lam < 1e-6] = 0.0
        jac[lam < 1e-6] = 0.0
        return np.einsum("bp,bpe->be", jac, vals)

    def evaluate_guarded(self, Ef, A, Bc, lam):
        C = self.C
        vals = C.evaluate(Ef, self.z, self.T).imag
        # guard spurious-surface roundoff: tiny z-independent denominators
        nargs = np.abs(Ef @ C.narg_matrix.T).min(axis=1)
        om = np.abs(Ef @ C.r_a.T)
        rmin = (
            np.abs(Ef[:, C.r_line] - om).min(axis=1)
            if C.r_a.size
            else np.full(Ef.shape[0], np.inf)
        )
        bad = (np.minimum(nargs, rmin) < 1e-9) | ~np.all(np.isfinite(vals), axis=1)
        bad &= lam.reshape(-1) >= 1e-6
        if np.any(bad):
            idx = np.nonzero(bad)[0]
            P = lam.shape[1]
            bb, pp = np.divmod(idx, P)
            acc = 0
            for f in (1 + 1e-6, 1 - 1e-6):
                l2 = lam[bb, pp] * f
                q = A[bb] * l2[:, None] ** 2 + 2 * Bc[bb] * l2[:, None] + self.Cc
                acc = acc + C.evaluate(np.sqrt(q), self.z, self.T).imag
            vals[idx] = 0.5 * acc
            self.n_guarded = getattr(self, "n_guarded", 0) + idx.size
        # exactly degenerate points (e.g. a loop momentum identically zero on the
        # boundary rho_i -> 0 of the outer map, where the Jacobian vanishes)
        nonfinite = ~np.all(np.isfinite(vals), axis=1)
        if np.any(nonfinite):
            vals[nonfinite] = 0.0
            self.n_nonfinite = getattr(self, "n_nonfinite", 0) + int(nonfinite.sum())
        return vals

    def __call__(self, u):
        v, jac = self.outer(u)
        A, Bc = self.quad_coeffs(v)
        X, W, _ = self.centres(A, Bc, self.eps.min())
        x, wx = self.nodes(X, W)
        res = self.radial(v, x, wx) * (jac * self.norm)[:, None]
        if self.fit_weights is not None:
            res = np.concatenate([res @ self.fit_weights[:, None], res], axis=1)
        return res
