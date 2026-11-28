"""The report's named monomials, built by hand for the comparison with the generated bases only (gate V4).

Origin: report 03b sec:Qbasis (tab:Qbasis, the 97 monomials Q_alpha), the two relative-sign variants of Q_R and
Q_sigma of the plane-wave engine (``SCRATCHPAD/st/bases.py``: ``R_<f><r>-``, ``sig_<r>-``; with them the 109
monomials of ``st_btable_unit.json``), eq:QYM and eq:closed (the 25 closed directions), tab:Nbasis (the monomials
of N that appear in eq:btable-gauge and eq:btable-higgs) and eq:evanescent, eq:Nhat. Names are the engine names
of ``SCRATCHPAD/st/bases.py`` and ``bases_n.py``. Conventions of the report: bars omitted (four-dimensional
heads ``gbar``, ``gammabar``), ``psi_f = P_c Psi_f``, ``Am^(r)_mu = g_r T^a A^a_mu`` restricted to the multiplet,
``wm^(r) = g_r T^a omega^a``, ``fg = g_r f``; ``rhotilde`` is the shifted gauge source; ``Q_rho`` and ``Q_zeta``
tensors as listed after tab:Qbasis. This module is test data: the production code never uses it.
"""

from __future__ import annotations

from fractions import Fraction

from symbolica import Expression

from hepkit_renorm import Density, LocalFunctional, Occurrence
from hepkit_renorm import tensors as tn
from hepkit_renorm.action import yukawa_terms
from hepkit_renorm.bases import shift_gauge_sources

I = tn.I
FERMIONS = ("Q", "U", "D", "L", "E")
FACTORS_OF = {"Q": "GWB", "U": "GB", "D": "GB", "L": "WB", "E": "B"}


def num(q) -> Expression:
    q = Fraction(q)
    return Expression.num(q.numerator) / q.denominator


class Report:
    """Builders of the report's monomials for a model with the field content of the one-family SM."""

    def __init__(self, model) -> None:
        self.m = model
        self.F = model.fields
        self.c = {r: model.coupling(r) for r in model.factors}

    # ------------------------------------------------------------------ labels and occurrences
    def lor(self):
        return tn.fresh(tn.LORENTZ)

    def adj(self, r):
        fac = self.m.factors[r]
        return None if fac.abelian else tn.fresh(fac.adjoint_space())

    def gauge(self, r, mu, a, ders=()):
        return Occurrence(r, (mu,) if a is None else (mu, a), tuple(ders))

    def ghost(self, r, a, ders=()):
        return Occurrence(f"gh_{r}", () if a is None else (a,), tuple(ders))

    def func(self, *parts) -> LocalFunctional:
        total = Density(self.F)
        for coeff, occs in parts:
            total = total + Density.monomial(self.F, coeff, *occs)
        return shift_gauge_sources(self.m, LocalFunctional.integrate(total))

    # ------------------------------------------------------------------ matter helpers
    def matter_labels(self, name):
        return [tn.fresh(s.space) for s in self.F[name].slots]

    def generator(self, multiplet, r, a, row, col):
        """``g_r (T^a)_{row col}`` on the slot of factor ``r`` of the multiplet (charge for U(1)); deltas on the
        other group slots. ``row``/``col``: label lists of the conjugate and of the field."""
        fac = self.m.factors[r]
        coeff = self.c[r]
        info = self.F[multiplet]
        for k, s in enumerate(info.slots):
            if s.factor is None:
                continue
            if s.factor == r:
                coeff = coeff * tn.T(a, row[k], col[k])
            else:
                coeff = coeff * tn.delta(row[k], col[k])
        if fac.abelian:
            coeff = coeff * num(info.charges[r])
        return coeff

    def identity(self, multiplet, row, col):
        coeff = Expression.num(1)
        for k, s in enumerate(self.F[multiplet].slots):
            if s.factor is not None:
                coeff = coeff * tn.delta(row[k], col[k])
        return coeff

    def spinor(self, f, row, col, *gammas):
        """``Psibar P_cbar gammabar^{mu_1} .. P_c Psi`` between the spinor labels."""
        proj = tn.PL if self.F[f].chirality == "L" else tn.PR
        u = [tn.fresh(tn.SPINOR) for _ in range(len(gammas) + 1)]
        e = proj(u[-1], col[0])
        for k, mu in enumerate(gammas):
            e = e * tn.gammabar(mu, u[k], u[k + 1])
        return _row_proj(self.F[f].chirality)(row[0], u[0]) * e

    # ------------------------------------------------------------------ Q: gauge sector
    def T(self, r):
        mu, nu, rho, sig, a = self.lor(), self.lor(), self.lor(), self.lor(), self.adj(r)
        b = a
        g = tn.gbar
        half = num(Fraction(-1, 2))
        return self.func((half * g(mu, rho) * g(nu, sig), [self.gauge(r, nu, a, [mu]), self.gauge(r, sig, b, [rho])]),
                         (-half * g(mu, sig) * g(nu, rho), [self.gauge(r, nu, a, [mu]), self.gauge(r, sig, b, [rho])]))

    def L(self, r):
        mu, nu, a = self.lor(), self.lor(), self.adj(r)
        return self.func((num(Fraction(-1, 2)), [self.gauge(r, mu, a, [mu]), self.gauge(r, nu, a, [nu])]))

    def M(self, r):
        mu, a = self.lor(), self.adj(r)
        return self.func((num(Fraction(1, 2)), [self.gauge(r, mu, a), self.gauge(r, mu, a)]))

    def f3(self, r):
        mu, nu = self.lor(), self.lor()
        a, b, c = self.adj(r), self.adj(r), self.adj(r)
        return self.func((-self.c[r] * tn.f(a, b, c), [self.gauge(r, nu, a, [mu]), self.gauge(r, mu, b),
                                                      self.gauge(r, nu, c)]))

    def d3(self):
        mu, nu = self.lor(), self.lor()
        a, b, c = self.adj("G"), self.adj("G"), self.adj("G")
        return self.func((self.c["G"] ** 3 * tn.dsym(a, b, c), [self.gauge("G", mu, a, [mu]), self.gauge("G", nu, b),
                                                                self.gauge("G", nu, c)]))

    def B3(self):
        mu, nu = self.lor(), self.lor()
        return self.func((self.c["B"] ** 3, [self.gauge("B", mu, None, [mu]), self.gauge("B", nu, None),
                                             self.gauge("B", nu, None)]))

    def rB(self, r, k):
        mu, nu, a = self.lor(), self.lor(), self.adj(r)
        w = self.c[r] ** 2 * self.c["B"]
        if k == 1:
            return self.func((w, [self.gauge("B", mu, None, [mu]), self.gauge(r, nu, a), self.gauge(r, nu, a)]))
        if k == 2:
            return self.func((w, [self.gauge(r, mu, a, [mu]), self.gauge(r, nu, a), self.gauge("B", nu, None)]))
        if k == 3:
            return self.func((w, [self.gauge(r, mu, a), self.gauge(r, nu, a), self.gauge("B", nu, None, [mu])]))
        rho, sig = self.lor(), self.lor()
        return self.func((w * tn.eps(mu, nu, rho, sig), [self.gauge(r, mu, a), self.gauge(r, rho, a, [nu]),
                                                          self.gauge("B", sig, None)]))

    def AA(self, r, s, mu, a, b):
        """Occurrences of ``A_(r)^a . A_(s)^b`` (one Lorentz label ``mu``)."""
        return [self.gauge(r, mu, a), self.gauge(s, mu, b)]

    def four(self, name):
        g = self.c
        mu, nu = self.lor(), self.lor()
        if name in ("G4a", "G4b", "W4a", "W4b"):
            r = name[0]
            a, b = self.adj(r), self.adj(r)
            if name.endswith("a"):
                occs = self.AA(r, r, mu, a, a) + self.AA(r, r, nu, b, b)
            else:
                occs = self.AA(r, r, mu, a, b) + self.AA(r, r, nu, a, b)
            return self.func((g[r] ** 4, occs))
        if name == "G4c":
            a, b, c, d, e = (self.adj("G") for _ in range(5))
            return self.func((g["G"] ** 4 * tn.dsym(a, b, e) * tn.dsym(c, d, e),
                              self.AA("G", "G", mu, a, b) + self.AA("G", "G", nu, c, d)))
        if name == "B4":
            return self.func((g["B"] ** 4, self.AA("B", "B", mu, None, None) + self.AA("B", "B", nu, None, None)))
        if name in ("GWa", "GWb", "GBa", "GBb", "WBa", "WBb"):
            r, s = name[0], name[1]
            a, b = self.adj(r), self.adj(s)
            w = g[r] ** 2 * g[s] ** 2
            if name.endswith("a"):
                return self.func((w, self.AA(r, r, mu, a, a) + self.AA(s, s, nu, b, b)))
            return self.func((w, self.AA(r, s, mu, a, b) + self.AA(r, s, nu, a, b)))
        if name == "GGGB":
            a, c, d = self.adj("G"), self.adj("G"), self.adj("G")
            return self.func((g["G"] ** 3 * g["B"] * tn.dsym(a, c, d),
                              self.AA("G", "G", mu, a, c) + self.AA("G", "B", nu, d, None)))
        if name in ("GGGBeps", "WWWBeps"):
            r = name[0]
            a, b, c = self.adj(r), self.adj(r), self.adj(r)
            rho, sig = self.lor(), self.lor()
            return self.func((g[r] ** 3 * g["B"] * tn.eps(mu, nu, rho, sig) * tn.f(a, b, c),
                              [self.gauge(r, mu, a), self.gauge(r, nu, b), self.gauge(r, rho, c),
                               self.gauge("B", sig, None)]))
        raise KeyError(name)

    # ------------------------------------------------------------------ Q: fermions and Higgs
    def K(self, f):
        row, col, mu = self.matter_labels(f + "bar"), self.matter_labels(f), self.lor()
        coeff = I * self.spinor(f, row, col, mu) * self.identity(f, row, col)
        return self.func((coeff, [Occurrence(f + "bar", tuple(row)), Occurrence(f, tuple(col), (mu,))]))

    def V(self, f, r):
        row, col, mu, a = self.matter_labels(f + "bar"), self.matter_labels(f), self.lor(), self.adj(r)
        coeff = self.spinor(f, row, col, mu) * self.generator(f, r, a, row, col)
        return self.func((coeff, [Occurrence(f + "bar", tuple(row)), Occurrence(f, tuple(col)), self.gauge(r, mu, a)]))

    def higgs(self, name):
        row, col = self.matter_labels("Hdag"), self.matter_labels("H")
        ident = self.identity("H", row, col)
        def Hd(d=()):
            return Occurrence("Hdag", tuple(row), tuple(d))

        def H(d=()):
            return Occurrence("H", tuple(col), tuple(d))

        mu = self.lor()
        if name == "SS":
            return self.func((ident, [Hd([mu]), H([mu])]))
        if name == "m":
            return self.func((ident, [Hd(), H()]))
        if name == "lam":
            row2, col2 = self.matter_labels("Hdag"), self.matter_labels("H")
            return self.func((ident * self.identity("H", row2, col2),
                              [Hd(), H(), Occurrence("Hdag", tuple(row2)), Occurrence("H", tuple(col2))]))
        if name in ("J_W", "J_B"):
            r = name[-1]
            a = self.adj(r)
            gen = self.generator("H", r, a, row, col)
            return self.func((I * gen, [Hd([mu]), H(), self.gauge(r, mu, a)]),
                             (-I * gen, [Hd(), H([mu]), self.gauge(r, mu, a)]))
        if name in ("D_W", "D_B"):
            r = name[-1]
            a = self.adj(r)
            return self.func((self.generator("H", r, a, row, col), [Hd(), H(), self.gauge(r, mu, a, [mu])]))
        if name in ("HWW", "HBB", "HWB"):
            r, s = {"HWW": "WW", "HBB": "BB", "HWB": "WB"}[name]
            mid = self.matter_labels("H")
            a, b = self.adj(r), self.adj(s)
            coeff = self.generator("H", r, a, row, mid) * self.generator("H", s, b, mid, col)
            if name == "HWB":
                coeff = coeff * 2
            return self.func((coeff, [Hd(), H(), self.gauge(r, mu, a), self.gauge(s, mu, b)]))
        if name == "HGG":
            a = self.adj("G")
            return self.func((self.c["G"] ** 2 * ident, [Hd(), H(), self.gauge("G", mu, a), self.gauge("G", mu, a)]))
        raise KeyError(name)

    def yukawa(self, k, sign):
        """``Q^+_{y_k} = -(O_k + O_k^dag)`` and ``Q^-_{y_k} = -i (O_k - O_k^dag)``, from
        ``L_Yuk = -sum y_k (O_k + h.c.)``."""
        y = self.m.param(f"y_{k}")
        lyuk = yukawa_terms(self.m)
        part = Density.from_monomials(self.F, [mo for mo in lyuk.monomials() if contains_param(mo.coefficient, y)])
        part = part.map_coefficients(lambda c: (c / y).expand())
        lead = {"u": "Qbar", "d": "Qbar", "e": "Lbar"}[k]
        o = Density.from_monomials(self.F, [mo for mo in part.monomials() if lead in mo.fields()])
        odag = part - o
        # part = -(O + O^dag): O = -o, O^dag = -odag
        if sign > 0:
            return LocalFunctional.integrate(o + odag)
        return LocalFunctional.integrate((o - odag).scale(I))

    # ------------------------------------------------------------------ Q: ghosts and sources
    def rho_d(self, r):
        mu, a = self.lor(), self.adj(r)
        src = Occurrence(f"src_{r}", (mu,) if a is None else (mu, a))
        return self.func((Expression.num(1), [src, self.ghost(r, a, [mu])]))

    def rho_t(self, tensor, x, y, z):
        """``int t^{abc} rhotilde^mu_a A^b_mu omega^c`` with ``t`` built by ``tensor(a, b, c)``."""
        mu, a, b, c = self.lor(), self.adj(x), self.adj(y), self.adj(z)
        src = Occurrence(f"src_{x}", (mu,) if a is None else (mu, a))
        return self.func((tensor(a, b, c), [src, self.gauge(y, mu, b), self.ghost(z, c)]))

    def zeta_t(self, tensor, x, y, z):
        a, b, c = self.adj(x), self.adj(y), self.adj(z)
        src = Occurrence(f"src_gh_{x}", () if a is None else (a,))
        return self.func((tensor(a, b, c), [src, self.ghost(y, b), self.ghost(z, c)]))

    def R(self, f, r, sign):
        """``int [i Rbar wm psi +- i omega^a psibar Tg^a R]`` (``src_f`` = Rbar, ``src_fbar`` = R)."""
        a = self.adj(r)
        row, col = self.matter_labels(f"src_{f}"), self.matter_labels(f)
        proj = tn.PL if self.F[f].chirality == "L" else tn.PR
        first = I * self.generator(f, r, a, row, col) * proj(row[0], col[0])
        row2, col2 = self.matter_labels(f + "bar"), self.matter_labels(f"src_{f}bar")
        second = sign * I * self.generator(f, r, a, row2, col2) * _row_proj(self.F[f].chirality)(row2[0], col2[0])
        return self.func((first, [Occurrence(f"src_{f}", tuple(row)), self.ghost(r, a), Occurrence(f, tuple(col))]),
                         (second, [self.ghost(r, a), Occurrence(f + "bar", tuple(row2)),
                                   Occurrence(f"src_{f}bar", tuple(col2))]))

    def sigma(self, r, sign):
        """``int [i sigma^dag wm H -+ i omega^a H^dag Tg^a sigma]`` (``sign = +1``: the report's Q_sigma)."""
        a = self.adj(r)
        row, col = self.matter_labels("src_H"), self.matter_labels("H")
        row2, col2 = self.matter_labels("Hdag"), self.matter_labels("src_Hdag")
        return self.func((I * self.generator("H", r, a, row, col),
                          [Occurrence("src_H", tuple(row)), self.ghost(r, a), Occurrence("H", tuple(col))]),
                         (-sign * I * self.generator("H", r, a, row2, col2),
                          [self.ghost(r, a), Occurrence("Hdag", tuple(row2)), Occurrence("src_Hdag", tuple(col2))]))

    # ------------------------------------------------------------------ the named basis
    def q_basis(self) -> dict[str, LocalFunctional]:
        """The 109 monomials of ``st_btable_unit.json`` (the 97 of tab:Qbasis plus ``R_*-`` and ``sig_*-``)."""
        g = self.c
        out: dict[str, LocalFunctional] = {}
        for r in "GWB":
            out[f"T_{r}"], out[f"L_{r}"], out[f"M_{r}"] = self.T(r), self.L(r), self.M(r)
        for r in "GW":
            out[f"3f_{r}"] = self.f3(r)
        out["3d"], out["3B"] = self.d3(), self.B3()
        for r in "GW":
            for k in (1, 2, 3):
                out[f"{r}B{k}"] = self.rB(r, k)
            out[f"{r}Beps"] = self.rB(r, 4)
        for name in ("G4a", "G4b", "G4c", "W4a", "W4b", "B4", "GWa", "GWb", "GBa", "GBb", "WBa", "WBb", "GGGB",
                     "GGGBeps", "WWWBeps"):
            out[name] = self.four(name)
        for f in FERMIONS:
            out[f"K_{f}"] = self.K(f)
            for r in FACTORS_OF[f]:
                out[f"V_{f}{r}"] = self.V(f, r)
        for name in ("SS", "m", "lam", "J_W", "J_B", "D_W", "D_B", "HWW", "HBB", "HWB", "HGG"):
            out[name] = self.higgs(name)
        for k in "ude":
            out[f"y{k}+"], out[f"y{k}-"] = self.yukawa(k, 1), self.yukawa(k, -1)
        for r in "GWB":
            out[f"rd_{r}"] = self.rho_d(r)
        for r in "GW":
            out[f"rf_{r}"] = self.rho_t(lambda a, b, c, r=r: g[r] * tn.f(a, b, c), r, r, r)
        out["rd"] = self.rho_t(lambda a, b, c: g["G"] * tn.dsym(a, b, c), "G", "G", "G")
        for x, y, z in ("WWB", "BWW", "WBW", "GGB", "BGG", "GBG", "BBB"):
            def kron(a, b, c, x=x, y=y, z=z):
                if x == y == z:
                    return Expression.num(1)
                if x == y:
                    return tn.delta(a, b)
                if y == z:
                    return tn.delta(b, c)
                return tn.delta(a, c)
            out[f"r{x}{y}{z}"] = self.rho_t(kron, x, y, z)
        for r in "GW":
            out[f"zf_{r}"] = self.zeta_t(lambda a, b, c, r=r: num(Fraction(-1, 2)) * g[r] * tn.f(a, b, c), r, r, r)
        for r in "GW":
            out[f"z{r}{r}B"] = self.zeta_t(lambda a, b, c: tn.delta(a, b), r, r, "B")
        for f in FERMIONS:
            for r in FACTORS_OF[f]:
                out[f"R_{f}{r}+"], out[f"R_{f}{r}-"] = self.R(f, r, 1), self.R(f, r, -1)
        for r in "WB":
            out[f"sig_{r}+"], out[f"sig_{r}-"] = self.sigma(r, 1), self.sigma(r, -1)
        return out

    @staticmethod
    def report_names(names) -> list[str]:
        """The 97 monomials of tab:Qbasis (the engine's ``-`` variants dropped)."""
        return [n for n in names if not n.endswith("-") or n.startswith("y")]

    # ------------------------------------------------------------------ N: the monomials of eq:btable-gauge, -higgs
    def n_basis(self) -> dict[str, LocalFunctional]:
        """The monomials of tab:Nbasis that occur in eq:btable-gauge (two-point, 3f and quartic lines) and
        eq:btable-higgs."""
        g, out = self.c, {}
        for r in "GWB":
            mu, nu, a = self.lor(), self.lor(), self.adj(r)
            out[f"N1_{r}"] = self.func((Expression.num(1), [self.ghost(r, a), self.gauge(r, mu, a, [mu])]))
            out[f"N3_{r}"] = self.func((Expression.num(1), [self.ghost(r, a), self.gauge(r, mu, a, [mu, nu, nu])]))
        for r in "GW":
            mu, nu = self.lor(), self.lor()
            a, b, c, e, y = (self.adj(r) for _ in range(5))
            fg = g[r] * tn.f(a, b, c)
            out[f"Nf1_{r}"] = self.func((fg, [self.ghost(r, a), self.gauge(r, mu, b), self.gauge(r, mu, c, [nu, nu])]))
            out[f"Nf2_{r}"] = self.func((fg, [self.ghost(r, a), self.gauge(r, mu, b), self.gauge(r, nu, c, [mu, nu])]))
            ff = g[r] ** 2 * tn.f(e, a, b) * tn.f(e, c, y)
            out[f"Nff_{r}"] = self.func((ff, [self.ghost(r, a, [mu]), self.gauge(r, nu, b), self.gauge(r, mu, c),
                                              self.gauge(r, nu, y)]))
        for name, r in (("NG4a", "G"), ("NG4b", "G"), ("NW4a", "W"), ("NB4", "B")):
            mu, nu = self.lor(), self.lor()
            a, b = self.adj(r), self.adj(r)
            if name.endswith("b"):
                occs = [self.ghost(r, a, [mu]), self.gauge(r, mu, b), self.gauge(r, nu, a), self.gauge(r, nu, b)]
            else:
                occs = [self.ghost(r, a, [mu]), self.gauge(r, mu, a), self.gauge(r, nu, b), self.gauge(r, nu, b)]
            out[name] = self.func((g[r] ** 4, occs))
        for r, s in ("GW", "GB", "WB"):
            mu, nu = self.lor(), self.lor()
            a, b = self.adj(r), self.adj(s)
            w = g[r] ** 2 * g[s] ** 2
            out[f"N{r}{s}a"] = self.func((w, [self.ghost(r, a, [mu]), self.gauge(r, mu, a), self.gauge(s, nu, b),
                                              self.gauge(s, nu, b)]))
            out[f"N{r}{s}b"] = self.func((w, [self.gauge(r, mu, a), self.gauge(r, mu, a), self.ghost(s, b, [nu]),
                                              self.gauge(s, nu, b)]))
            out[f"N{r}{s}c"] = self.func((w, [self.ghost(r, a, [mu]), self.gauge(s, mu, b), self.gauge(r, nu, a),
                                              self.gauge(s, nu, b)]))
            out[f"N{r}{s}d"] = self.func((w, [self.gauge(r, mu, a), self.ghost(s, b, [mu]), self.gauge(r, nu, a),
                                              self.gauge(s, nu, b)]))
        row, col, mid = self.matter_labels("Hdag"), self.matter_labels("H"), self.matter_labels("H")
        def Hd(d=()):
            return Occurrence("Hdag", tuple(row), tuple(d))

        def H(d=()):
            return Occurrence("H", tuple(col), tuple(d))
        for r in "WB":
            mu, nu, a = self.lor(), self.lor(), self.adj(r)
            gen = self.generator("H", r, a, row, col)
            out[f"NS_{r}"] = self.func((I * gen, [Hd([mu]), self.ghost(r, a, [mu]), H()]),
                                       (-I * gen, [Hd(), self.ghost(r, a, [mu]), H([mu])]))
            out[f"NS2_{r}"] = self.func((gen, [Hd(), self.ghost(r, a, [mu, mu]), H()]))
            for s in "WB":
                b = self.adj(s)
                first = self.generator("H", r, a, row, mid) * self.generator("H", s, b, mid, col)
                second = self.generator("H", s, b, row, mid) * self.generator("H", r, a, mid, col)
                occs = [Hd(), self.ghost(r, a, [mu]), self.gauge(s, mu, b), H()]
                out[f"NSV_{r}{s}"] = self.func((first + second, occs))
        mu, a, b = self.lor(), self.adj("W"), self.adj("W")
        comm = self.generator("H", "W", b, row, mid) * self.generator("H", "W", a, mid, col) - \
            self.generator("H", "W", a, row, mid) * self.generator("H", "W", b, mid, col)
        out["NSC"] = self.func((-I * comm, [Hd(), self.ghost("W", a, [mu]), self.gauge("W", mu, b), H()]))
        a = self.adj("G")
        out["NSG"] = self.func((g["G"] ** 2 * self.identity("H", row, col),
                                [Hd(), H(), self.ghost("G", a, [mu]), self.gauge("G", mu, a)]))
        return out

    # ------------------------------------------------------------------ eq:closed
    def closed(self, q: dict[str, LocalFunctional]) -> dict[str, LocalFunctional]:
        """The 25 closed directions of eq:closed (with eq:QYM)."""
        g, m = self.c, self.m

        def comb(*terms):
            total = Density(self.F)
            for c, name in terms:
                total = total + q[name].scale(c if isinstance(c, Expression) else num(c))
            return LocalFunctional.integrate(total)

        one = Expression.num(1)
        ym = {"G": [(-one / (4 * g["G"] ** 2) * num(Fraction(1, 2)), "G4a"), (one / (4 * g["G"] ** 2), "G4b"),
                    (-one / (4 * g["G"] ** 2) * num(Fraction(3, 2)), "G4c")],
              "W": [(-one / (4 * g["W"] ** 2), "W4a"), (one / (4 * g["W"] ** 2), "W4b")]}
        h = {"W": [(-1, "J_W"), (2, "HWW"), (1, "HWB")], "B": [(-1, "J_B"), (2, "HBB"), (1, "HWB")], "G": []}
        out = {}
        for r in "GW":
            vs = [(1, f"V_{f}{r}") for f in FERMIONS if r in FACTORS_OF[f]]
            rs = [(1, f"R_{f}{r}+") for f in FERMIONS if r in FACTORS_OF[f]]
            sig = [(1, "sig_W+")] if r == "W" else []
            ymr = [(c * 4, n) for c, n in ym[r]]
            out[f"O_A_{r}"] = comb((2, f"T_{r}"), (3, f"3f_{r}"), *ymr, *vs, *h[r], (-1, f"rd_{r}"))
            out[f"O_c_{r}"] = comb((1, f"rd_{r}"), (1, f"rf_{r}"), (1, f"zf_{r}"), *rs, *sig)
            out[f"O_g_{r}"] = comb((1, f"3f_{r}"), *[(c * 2, n) for c, n in ym[r]], *vs, *rs, *h[r], *sig,
                                   (1, f"rf_{r}"), (1, f"zf_{r}"))
        out["T_B"] = q["T_B"]
        out["O_A_B"] = comb((2, "T_B"), *[(1, f"V_{f}B") for f in FERMIONS], *h["B"], (-1, "rd_B"))
        out["O_c_B"] = comb((1, "rd_B"), *[(1, f"R_{f}B+") for f in FERMIONS], (1, "sig_B+"))
        yuk_of = {"Q": "ud", "U": "u", "D": "d", "L": "e", "E": "e"}
        for f in FERMIONS:
            out[f"O_psi_{f}"] = comb((2, f"K_{f}"), *[(2, f"V_{f}{r}") for r in FACTORS_OF[f]],
                                     *[(m.param(f"y_{k}"), f"y{k}+") for k in yuk_of[f]])
        out["O_H"] = comb((2, "SS"), (-2, "J_W"), (-2, "J_B"), (2, "HWW"), (2, "HBB"), (2, "HWB"),
                          (-2 * m.param("m2"), "m"), (-4 * m.param("lambda"), "lam"),
                          *[(m.param(f"y_{k}"), f"y{k}+") for k in "ude"])
        for k in "ude":
            out[f"y{k}+"], out[f"y{k}-"] = q[f"y{k}+"], q[f"y{k}-"]
        out["lam"], out["m"] = q["lam"], q["m"]
        hyper = {f: m.fields[f].charges["B"] for f in FERMIONS}
        for qname, charges in (("baryon", {"Q": Fraction(1, 3), "U": Fraction(1, 3), "D": Fraction(1, 3)}),
                               ("lepton", {"L": 1, "E": 1})):
            terms = []
            for f, qf in charges.items():
                terms += [(Fraction(qf) / hyper[f], f"V_{f}B"), (Fraction(qf) / hyper[f], f"R_{f}B+")]
            out[f"O_{qname}"] = comb(*terms)
        return out


def _row_proj(chirality):
    """Projector on the right of ``Psibar`` (``P_cbar``)."""
    return tn.PR if chirality == "L" else tn.PL


def contains_param(expr: Expression, sym: Expression) -> bool:
    return expr.derivative(sym).expand() != Expression.num(0)
