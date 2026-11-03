"""Exact Matsubara sums for scalar self-energy graphs by iterated residues.

Independent of GammaLoop.  At fixed spatial loop momenta the Euclidean summand is
prod_k 1/(omega_k^2 + E_k^2) = prod_k 1/(E_k^2 - Omega_k^2), with Omega_k = i*omega_k a
linear combination of loop variables z_l = i*omega_l and of the external z_K = i*omega_K.
Each loop sum T sum_n f(i omega_n) = -sum_{poles of f} Res[f(z) n_B(z)] is done in turn.
At a pole Omega_k = sigma*E_k of propagator k (sigma=+-1) the residue gives
sigma*c_kl/(2E_k) * n_B(z_l*) with z_l* the pole position.  Because every remaining
imaginary part of z_l* is a bosonic Matsubara frequency (other loops and the external
one), n_B(z_l*) = n_B(Re z_l*): this is the replacement n_B(x + i omega_ext) = n_B(x)
applied BEFORE any analytic continuation.  The result is a rational function of z_K
whose only z_K dependence sits in the uncut propagators.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction

import numpy as np

H_MASS, L_MASS = 0.6, 0.25


@dataclass
class Graph:
    name: str
    n_loops: int
    # line momentum = sum_l loop[l]*P_l + ext*K, as integer coefficient tuples
    loop_coeff: list[tuple[int, ...]]
    ext_coeff: list[int]
    masses: list[float]
    kappa: float
    labels: list[str] = field(default_factory=list)

    @property
    def n_lines(self):
        return len(self.masses)


def double_triangle() -> Graph:
    # Routing of reference_derivation.typ: p1=p, p2=q, p3=p-K, p4=q-K, p5=p-q
    H, L = H_MASS, L_MASS
    return Graph(
        "dt",
        2,
        loop_coeff=[(1, 0), (0, 1), (1, 0), (0, 1), (1, -1)],
        ext_coeff=[0, 0, -1, -1, 0],
        masses=[H, L, H, L, H],
        kappa=0.5,
        labels=["p1=P (H)", "p2=Q (L)", "p3=P-K (H)", "p4=Q-K (L)", "p5=P-Q (H)"],
    )


def triangle_box_triangle() -> Graph:
    # Edges of triangle_box_triangle.dot with loops P=e0, Q=e3, R=e6 and K in at v7.
    H, L = H_MASS, L_MASS
    return Graph(
        "tbt",
        3,
        loop_coeff=[
            (1, 0, 0),  # e0 v7->v8   P      H
            (-1, 0, 0),  # e1 v7->v12  K-P    H
            (1, -1, 0),  # e2 v8->v12  P-Q    L
            (0, 1, 0),  # e3 v8->v9   Q      H
            (0, -1, 0),  # e4 v12->v11 K-Q    H
            (0, -1, 1),  # e5 v11->v9  R-Q    H
            (0, 0, 1),  # e6 v9->v10  R      L
            (0, 0, -1),  # e7 v11->v10 K-R    L
        ],
        ext_coeff=[0, 1, 0, 0, 1, 0, 0, 1],
        masses=[H, H, L, H, H, H, L, L],
        kappa=0.5,
        labels=[
            "e0 P (H)",
            "e1 K-P (H)",
            "e2 P-Q (L)",
            "e3 Q (H)",
            "e4 K-Q (H)",
            "e5 R-Q (H)",
            "e6 R (L)",
            "e7 K-R (L)",
        ],
    )


def bubble(m1: float, m2: float) -> Graph:
    """One-loop control: lines P (m1) and K-P (m2), distinct daughters, kappa=1."""
    return Graph(
        "bubble",
        1,
        loop_coeff=[(1,), (-1,)],
        ext_coeff=[0, 1],
        masses=[m1, m2],
        kappa=1.0,
    )


GRAPHS = {"dt": double_triangle, "tbt": triangle_box_triangle}


def get_graph(name: str) -> Graph:
    """'dt', 'tbt' or 'bubble:m1:m2'."""
    if name.startswith("bubble"):
        _, m1, m2 = name.split(":")
        return bubble(float(m1), float(m2))
    return GRAPHS[name]()


# ---------------------------------------------------------------------------
# Term algebra.  A propagator is (line k, zc, ec) with Omega = zc . (z_1..z_L, z_K)
# + ec . (E_1..E_N).  A term is (pref, cut lines, nB args, propagators).


@dataclass(frozen=True)
class Term:
    pref: Fraction  # rational prefactor (product of sigma*c signs)
    cuts: tuple[int, ...]  # lines that supply 1/(2E)
    nargs: tuple[tuple[int, ...], ...]  # n_B(sum_j a_j E_j)
    props: tuple[tuple[int, tuple[int, ...], tuple[int, ...]], ...]


def initial_term(g: Graph) -> Term:
    props = []
    for k in range(g.n_lines):
        zc = tuple(g.loop_coeff[k]) + (g.ext_coeff[k],)
        props.append((k, zc, tuple([0] * g.n_lines)))
    return Term(Fraction(1), (), (), tuple(props))


def sum_loop(terms: list[Term], loop: int, n_lines: int) -> list[Term]:
    out = []
    for t in terms:
        for idx, (k, zc, ec) in enumerate(t.props):
            c = zc[loop]
            if c == 0:
                continue
            assert c in (1, -1)
            for sigma in (1, -1):
                # pole: Omega_k = sigma E_k  =>  z_loop = c*(sigma E_k - rest)
                # z part of z_loop: -c*zc[m] for m != loop; energy: c*sigma at k, -c*ec
                zsub = [-c * zc[m] if m != loop else 0 for m in range(len(zc))]
                esub = [-c * ec[j] for j in range(n_lines)]
                esub[k] += c * sigma
                # n_B argument: the energy part only (all z are bosonic Matsubara)
                narg = tuple(esub)
                new_props = []
                for jdx, (k2, zc2, ec2) in enumerate(t.props):
                    if jdx == idx:
                        continue
                    a = zc2[loop]
                    if a == 0:
                        new_props.append((k2, zc2, ec2))
                        continue
                    zc_new = tuple(
                        (zc2[m] + a * zsub[m]) if m != loop else 0
                        for m in range(len(zc2))
                    )
                    ec_new = tuple(ec2[j] + a * esub[j] for j in range(n_lines))
                    new_props.append((k2, zc_new, ec_new))
                out.append(
                    Term(
                        t.pref * sigma * c,
                        t.cuts + (k,),
                        t.nargs + (narg,),
                        tuple(new_props),
                    )
                )
    return out


def build_terms(g: Graph, order=None) -> list[Term]:
    terms = [initial_term(g)]
    for loop in order if order is not None else range(g.n_loops):
        terms = sum_loop(terms, loop, g.n_lines)
    for t in terms:
        for k, zc, ec in t.props:
            assert all(v == 0 for v in zc[:-1])
    return terms


# ---------------------------------------------------------------------------
# Numerical evaluation


def n_bose(x, T):
    """Bose function; T=0 gives the vacuum step -theta(-x)."""
    x = np.asarray(x, dtype=float)
    if T == 0:
        return np.where(x < 0, -1.0, 0.0)
    with np.errstate(over="ignore"):
        return 1.0 / np.expm1(x / T)


class Compiled:
    """Vectorised evaluator for R(z) = T^L sum prod G at fixed spatial momenta."""

    def __init__(self, g: Graph, terms: list[Term]):
        self.g = g
        self.terms = terms
        N = g.n_lines
        nargs = sorted({a for t in terms for a in t.nargs})
        self.narg_index = {a: i for i, a in enumerate(nargs)}
        self.narg_matrix = np.array(nargs, dtype=float).reshape(len(nargs), N)
        zprops, rprops = {}, {}
        for t in terms:
            for k, zc, ec in t.props:
                key = (k, zc[-1], ec)
                (zprops if zc[-1] != 0 else rprops).setdefault(
                    key, len(zprops if zc[-1] != 0 else rprops)
                )
        self.zprop_keys = list(zprops)
        self.rprop_keys = list(rprops)
        self.zprop_index = zprops
        self.rprop_index = rprops
        # z-dependent propagators: 1/(E_k^2 - (c z + a.E)^2)
        self.z_line = np.array([k for k, c, e in self.zprop_keys], dtype=int)
        self.z_c = np.array([c for k, c, e in self.zprop_keys], dtype=float)
        self.z_a = np.array([e for k, c, e in self.zprop_keys], dtype=float).reshape(
            -1, N
        )
        self.r_line = np.array([k for k, c, e in self.rprop_keys], dtype=int)
        self.r_a = np.array([e for k, c, e in self.rprop_keys], dtype=float).reshape(
            -1, N
        )
        nt = len(terms)
        self.pref = np.array([float(t.pref) for t in terms])
        self.cut_mask = np.zeros((nt, N))
        maxn = max(len(t.nargs) for t in terms)
        maxz = max(sum(1 for p in t.props if p[1][-1] != 0) for t in terms)
        maxr = max(max(sum(1 for p in t.props if p[1][-1] == 0) for t in terms), 1)
        self.n_idx = np.zeros((nt, maxn), dtype=int)
        self.z_idx = np.full((nt, maxz), len(self.zprop_keys), dtype=int)  # pad -> 1
        self.r_idx = np.full((nt, maxr), len(self.rprop_keys), dtype=int)
        for i, t in enumerate(terms):
            for k in t.cuts:
                self.cut_mask[i, k] += 1
            for j, a in enumerate(t.nargs):
                self.n_idx[i, j] = self.narg_index[a]
            zl = [
                self.zprop_index[(k, zc[-1], ec)]
                for k, zc, ec in t.props
                if zc[-1] != 0
            ]
            rl = [
                self.rprop_index[(k, zc[-1], ec)]
                for k, zc, ec in t.props
                if zc[-1] == 0
            ]
            self.z_idx[i, : len(zl)] = zl
            self.r_idx[i, : len(rl)] = rl
        self.n_terms = nt
        # group terms by their z-propagator tuple to share products
        self.thresholds = self._thresholds()

    def _thresholds(self):
        """Signed line combinations h with a pole of R at z = h (z-dependent props)."""
        hs = set()
        for k, c, ec in self.zprop_keys:
            for sigma in (1, -1):
                # E_k^2 = (c z + a.E)^2 -> z = c*(sigma E_k - a.E)
                h = [-c * v for v in ec]
                h[k] += c * sigma
                hs.add(tuple(int(v) for v in h))
        return sorted(hs)

    def real_coefficients(self, E, T):
        """Per-term z-independent factor, shape (npts, nterms)."""
        npts = E.shape[0]
        nvals = n_bose(E @ self.narg_matrix.T, T)  # (npts, nargs)
        inv2E = 1.0 / (2.0 * E)
        # product over cut lines of 1/(2E)
        logpart = np.log(inv2E) @ self.cut_mask.T  # (npts, nterms)
        coef = self.pref[None, :] * np.exp(logpart)
        coef = coef * np.prod(nvals[:, self.n_idx], axis=2)
        if len(self.rprop_keys):
            om = E @ self.r_a.T
            Ek = E[:, self.r_line]
            rvals = 1.0 / ((Ek - om) * (Ek + om))
            rvals = np.concatenate([rvals, np.ones((npts, 1))], axis=1)
            coef = coef * np.prod(rvals[:, self.r_idx], axis=2)
        return coef

    def _merge(self):
        """Merge terms with identical z-propagator products (sum their real coefficients)."""
        keys = [tuple(sorted(r)) for r in self.z_idx.tolist()]
        uniq = sorted(set(keys))
        pos = {k: i for i, k in enumerate(uniq)}
        M = np.zeros((self.n_terms, len(uniq)))
        for i, k in enumerate(keys):
            M[i, pos[k]] = 1.0
        self.merge_matrix = M
        self.uz_idx = np.array(uniq, dtype=int).T.copy()  # (maxz, nuniq)

    def evaluate_ref(self, E, z, T, coef=None, chunk=512):
        """R(z) for energies E (npts, N) and complex z array (nz,) -> (npts, nz)."""
        if not hasattr(self, "merge_matrix"):
            self._merge()
        z = np.atleast_1d(np.asarray(z, dtype=complex))
        npts = E.shape[0]
        out = np.empty((npts, z.size), dtype=complex)
        nzp = len(self.zprop_keys)
        for s in range(0, npts, chunk):
            Ec = E[s : s + chunk]
            n = Ec.shape[0]
            cf = (
                self.real_coefficients(Ec, T) if coef is None else coef[s : s + chunk]
            ) @ self.merge_matrix
            om0 = (Ec @ self.z_a.T).T  # (nzp, n)
            Ek = Ec[:, self.z_line].T[:, :, None]
            zv = np.ones((nzp + 1, n, z.size), dtype=complex)
            omz = self.z_c[:, None, None] * z[None, None, :] + om0[:, :, None]
            zv[:nzp] = 1.0 / ((Ek - omz) * (Ek + omz))
            prod = zv[self.uz_idx[0]]
            for j in range(1, self.uz_idx.shape[0]):
                prod *= zv[self.uz_idx[j]]
            out[s : s + chunk] = np.einsum("pt,tpz->pz", cf, prod)
        return out

    # ------------------------------------------------------------------ fast path
    def _prepare_fast(self):
        """Factor-major gathers, reduceat merging of equal z-products, prefix-tree products."""
        N = self.g.n_lines
        nn = len(self.narg_index)
        nr = len(self.rprop_keys)
        one = N + nn + nr  # index of the row of ones in the real factor table
        rows = []
        for i, t in enumerate(self.terms):
            r = [k for k in t.cuts] + [N + self.narg_index[a] for a in t.nargs]
            r += [
                N + nn + self.rprop_index[(k, zc[-1], ec)]
                for k, zc, ec in t.props
                if zc[-1] == 0
            ]
            rows.append(r)
        width = max(len(r) for r in rows)
        ridx = np.full((len(rows), width), one, dtype=int)
        for i, r in enumerate(rows):
            ridx[i, : len(r)] = r
        ztup = [
            tuple(sorted(int(x) for x in row if x != len(self.zprop_keys)))
            for row in self.z_idx.tolist()
        ]
        order = sorted(range(len(ztup)), key=lambda i: ztup[i])
        self.f_order = np.array(order, dtype=int)
        self.f_ridx = ridx[self.f_order].T.copy()  # (width, nt) in merged order
        self.f_pref = self.pref[self.f_order]
        uniq, starts = [], []
        for pos, i in enumerate(order):
            if not uniq or ztup[i] != uniq[-1]:
                uniq.append(ztup[i])
                starts.append(pos)
        self.f_starts = np.array(starts, dtype=int)
        self.f_one = one
        # prefix tree over the sorted factor tuples
        node_of = {(): -1}
        levels = []  # per depth: (parent node ids, factor ids)
        node_count = 0
        maxlen = max(len(u) for u in uniq)
        for d in range(1, maxlen + 1):
            par, fac, ids = [], [], []
            for u in uniq:
                if len(u) >= d and u[:d] not in node_of:
                    node_of[u[:d]] = node_count
                    node_count += 1
                    par.append(node_of[u[: d - 1]])
                    fac.append(u[d - 1])
                    ids.append(node_of[u[:d]])
            levels.append(
                (
                    np.array(par, dtype=int),
                    np.array(fac, dtype=int),
                    np.array(ids, dtype=int),
                )
            )
        self.f_levels = levels
        self.f_nnodes = node_count
        self.f_leaf = np.array([node_of[u] for u in uniq], dtype=int)

    def evaluate(self, E, z, T, chunk=512):
        """R(z) for energies E (npts, N) and complex z array (nz,) -> (npts, nz)."""
        if not hasattr(self, "f_levels"):
            self._prepare_fast()
        z = np.atleast_1d(np.asarray(z, dtype=complex))
        nz = z.size
        npts = E.shape[0]
        out = np.empty((npts, nz), dtype=complex)
        N = self.g.n_lines
        for s in range(0, npts, chunk):
            Et = np.ascontiguousarray(E[s : s + chunk].T)  # (N, n)
            n = Et.shape[1]
            # real factor table
            F = np.empty((self.f_one + 1, n))
            F[:N] = 0.5 / Et
            nn = self.narg_matrix.shape[0]
            F[N : N + nn] = n_bose(self.narg_matrix @ Et, T)
            if self.r_a.size:
                om = self.r_a @ Et
                Ek = Et[self.r_line]
                F[N + nn : self.f_one] = 1.0 / ((Ek - om) * (Ek + om))
            F[self.f_one] = 1.0
            coef = F[self.f_ridx[0]] * self.f_pref[:, None]
            for j in range(1, self.f_ridx.shape[0]):
                coef *= F[self.f_ridx[j]]
            merged = np.add.reduceat(coef, self.f_starts, axis=0)  # (nuniq, n)
            # z-dependent propagators
            om0 = self.z_a @ Et  # (nzp, n)
            Ek = Et[self.z_line][:, :, None]
            omz = self.z_c[:, None, None] * z[None, None, :] + om0[:, :, None]
            zv = 1.0 / ((Ek - omz) * (Ek + omz))  # (nzp, n, nz)
            vals = np.empty((self.f_nnodes, n, nz), dtype=complex)
            for d, (par, fac, ids) in enumerate(self.f_levels):
                if d == 0:
                    vals[ids] = zv[fac]
                else:
                    vals[ids] = vals[par] * zv[fac]
            leaf = vals[self.f_leaf]  # (nuniq, n, nz)
            res = np.einsum(
                "up,upzc->pzc", merged, leaf.view(float).reshape(leaf.shape + (2,))
            )
            out[s : s + chunk] = res[..., 0] + 1j * res[..., 1]
        return out


def line_energies(g: Graph, loops, kvec):
    """loops: (npts, L, 3); kvec: (3,) -> E (npts, N)."""
    lc = np.array(g.loop_coeff, dtype=float)  # (N, L)
    ec = np.array(g.ext_coeff, dtype=float)
    mom = np.einsum("nl,pld->pnd", lc, loops) + ec[None, :, None] * kvec[None, None, :]
    m = np.array(g.masses)
    return np.sqrt(np.sum(mom**2, axis=2) + m[None, :] ** 2)


def brute_force(g: Graph, E, T, m_ext, nmax):
    """Direct truncated Matsubara sum at omega_K = 2 pi T m_ext for one point E (N,)."""
    L = g.n_loops
    w = 2 * np.pi * T * np.arange(-nmax, nmax + 1)
    grids = np.meshgrid(*([w] * L), indexing="ij")
    wK = 2 * np.pi * T * m_ext
    total = np.ones_like(grids[0])
    for k in range(g.n_lines):
        om = sum(g.loop_coeff[k][l] * grids[l] for l in range(L)) + g.ext_coeff[k] * wK
        total = total / (om**2 + E[k] ** 2)
    return T**L * total.sum()
