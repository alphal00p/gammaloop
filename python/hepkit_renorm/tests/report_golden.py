"""Golden functionals of the report for the algebraic end of the recipe (gates V7, V8): the singular counterterm
``deltaS^(1)``, the breaking ``X_1`` and the restoring counterterm ``DeltaS^(1)_ST`` as local functionals, built
from ``tests/data/final_counterterms.json`` (a verbatim copy of ``SCRATCHPAD/phase4/final_counterterms.json``) and
the report's named monomials.

Origin of the definitions: report 03b tab:Qbasis, eq:closed (through ``report_bases.Report``), tab:Nbasis and the
plane-wave engine ``SCRATCHPAD/st/bases_n.py`` (the N monomials that ``X_1`` needs beyond ``report_bases``: ``F``,
``FV``, ``VF``, ``RR``, ``RRb``, the Yukawa ghost insertions and ``GGGB1-3``), eq:evanescent (``E_A``, ``E_H``,
``Etilde_F``, ``Etilde_H``) and 04z eq:sts-anom (the anomaly part of ``X_1`` for an arbitrary fermion content).
Units: ``C`` with ``deltaS^(1) = C / (16 pi^2 eps)``; ``X`` and ``DeltaS`` in units of ``1/(16 pi^2)``.
This module is test data: the production code never uses it.
"""

from __future__ import annotations

from fractions import Fraction

from report_bases import FACTORS_OF, FERMIONS, Report, num
from symbolica import Expression, S

from hepkit_renorm import Density, LocalFunctional, Occurrence
from hepkit_renorm import tensors as tn
from hepkit_renorm.action import yukawa_terms
from hepkit_renorm.bases import BasisElement, OperatorBasis, continuation, external_classes, kinetic_invariants

I = tn.I
_NAMES = {"g_s": "g_s", "g": "g", "gp": "gp", "g'": "gp", "y_u": "y_u", "y_d": "y_d", "y_e": "y_e",
          "lambda": "lambda", "m2": "m2"}


def parse(model, text: str) -> Expression:
    """A sympy-style string of final_counterterms.json as an expression in the model's parameter symbols."""
    t = text.replace("g'", "gp").replace("**", "^")
    e = Expression.parse(t)
    for name in ("g_s", "gp", "y_u", "y_d", "y_e", "lambda", "m2", "g"):
        e = e.replace(S(name), model.param(name))
    return e.replace(S("I"), I).expand()


def _occ(name, labels, ders=()):
    return Occurrence(name, tuple(labels), tuple(ders))


class Golden(Report):
    """The report's functionals for a model with the field content of the one-family SM."""

    # ------------------------------------------------------------------ N monomials beyond report_bases
    def NF(self, f, r):
        row, col, mu, a = self.matter_labels(f + "bar"), self.matter_labels(f), self.lor(), self.adj(r)
        coeff = self.spinor(f, row, col, mu) * self.generator(f, r, a, row, col)
        return self.func((coeff, [_occ(f + "bar", row), self.ghost(r, a, [mu]), _occ(f, col)]))

    def FV(self, f, r, s, ghost_first=True):
        row, col, mid = self.matter_labels(f + "bar"), self.matter_labels(f), self.matter_labels(f)
        mu, a, b = self.lor(), self.adj(r), self.adj(s)
        if ghost_first:     # psibar gamma wm^(r) Am^(s) psi
            coeff = self.spinor(f, row, col, mu) * self.generator(f, r, a, row, mid) * self.generator(f, s, b, mid, col)
        else:               # psibar gamma Am^(s) wm^(r) psi
            coeff = self.spinor(f, row, col, mu) * self.generator(f, s, b, row, mid) * self.generator(f, r, a, mid, col)
        return self.func((coeff, [_occ(f + "bar", row), self.ghost(r, a), self.gauge(s, mu, b), _occ(f, col)]))

    def RR(self, f, r, s, bar=False):
        a, b = self.adj(r), self.adj(s)
        if not bar:         # Rbar wm^(r) wm^(s) psi
            row, col, mid = self.matter_labels(f"src_{f}"), self.matter_labels(f), self.matter_labels(f)
            proj = tn.PL if self.F[f].chirality == "L" else tn.PR
            coeff = proj(row[0], col[0]) * self.generator(f, r, a, row, mid) * self.generator(f, s, b, mid, col)
            return self.func((coeff, [_occ(f"src_{f}", row), self.ghost(r, a), self.ghost(s, b), _occ(f, col)]))
        row, col, mid = self.matter_labels(f + "bar"), self.matter_labels(f"src_{f}bar"), self.matter_labels(f)
        proj = tn.PR if self.F[f].chirality == "L" else tn.PL
        coeff = proj(row[0], col[0]) * self.generator(f, r, a, row, mid) * self.generator(f, s, b, mid, col)
        return self.func((coeff, [_occ(f + "bar", row), self.ghost(r, a), self.ghost(s, b), _occ(f"src_{f}bar", col)]))

    def yukawa_O(self, k, dagger=False) -> Density:
        """``O_k`` (or ``O_k^dagger``) with ``L_Yuk = -sum y_k (O_k + O_k^dagger)``."""
        y = self.m.param(f"y_{k}")
        lyuk = yukawa_terms(self.m)
        part = [mo for mo in lyuk.monomials() if mo.coefficient.derivative(y).expand() != 0]
        part = Density.from_monomials(self.F, part).map_coefficients(lambda c: (-c / y).expand())
        lead = {"u": "Qbar", "d": "Qbar", "e": "Lbar"}[k]
        sel = [mo for mo in part.monomials() if (lead in mo.fields()) != dagger]
        return Density.from_monomials(self.F, sel)

    def Y(self, k, dagger, side, r):
        """``i y_k O[psibar -> psibar wm^(r)]`` (side ``"L"``) or ``i y_k O[psi -> wm^(r) psi]`` (``"R"``)."""
        y = self.m.param(f"y_{k}")
        total = Density(self.F)
        for mo in self.yukawa_O(k, dagger).monomials():
            stored = list(mo.occurrences)
            kinds = [self.F[o.field].kind for o in stored]
            ib, ip = kinds.index("fermion_bar"), kinds.index("fermion")
            sign = -1 if ip < ib else 1          # written order: psibar first, psi last (one odd exchange)
            occs = [stored[ib]] + [o for k, o in enumerate(stored) if k not in (ib, ip)] + [stored[ip]]
            pos = 0 if side == "L" else len(occs) - 1
            o = occs[pos]
            mult = self.F[o.field].multiplet
            new = self.matter_labels(o.field)
            new = (o.slots[0],) + tuple(new[1:])           # keep the spinor label, rotate the group labels
            a = self.adj(r)
            if side == "L":
                gen = self.generator(mult, r, a, list(new), list(o.slots))
                occs[pos:pos + 1] = [Occurrence(o.field, new, o.derivatives), self.ghost(r, a)]
            else:
                gen = self.generator(mult, r, a, list(o.slots), list(new))
                occs[pos:pos + 1] = [self.ghost(r, a), Occurrence(o.field, new, o.derivatives)]
            total = total + Density.monomial(self.F, I * y * sign * mo.coefficient * gen, *occs)
        return LocalFunctional.integrate(total)

    def GGGB(self, k):
        g = self.c
        mu, nu = self.lor(), self.lor()
        a, b, c = self.adj("G"), self.adj("G"), self.adj("G")
        w = g["G"] ** 3 * g["B"] * tn.dsym(a, b, c)
        if k == 1:      # (d omega^A . G^B)(G^C . B)
            occs = [self.ghost("G", a, [mu]), self.gauge("G", mu, b), self.gauge("G", nu, c), self.gauge("B", nu, None)]
        elif k == 2:    # (G^A . G^B)(G^C . d omega^B)
            occs = [self.gauge("G", mu, a), self.gauge("G", mu, b), self.gauge("G", nu, c), self.ghost("B", None, [nu])]
        else:           # (G^A . G^B)(d omega^C . B)
            occs = [self.gauge("G", mu, a), self.gauge("G", mu, b), self.ghost("G", c, [nu]), self.gauge("B", nu, None)]
        return self.func((w, occs))

    def eps_anomaly(self, kind, r=None):
        """The epsilon monomials of eq:sts-anom-basis (report tab:Nbasis cVV, cVVV epsilon rows)."""
        g = self.c
        mu, nu, al, be = self.lor(), self.lor(), self.lor(), self.lor()
        eps = tn.eps(mu, nu, al, be)
        if kind == "eps_B":      # g'^3 eps gh^B d B d B
            return self.func((g["B"] ** 3 * eps, [self.ghost("B", None), self.gauge("B", nu, None, [mu]),
                                                  self.gauge("B", be, None, [al])]))
        if kind == "eps_d":      # gs^3 d^ABC eps gh^A d G^B d G^C
            a, b, c = self.adj("G"), self.adj("G"), self.adj("G")
            return self.func((g["G"] ** 3 * tn.dsym(a, b, c) * eps,
                              [self.ghost("G", a), self.gauge("G", nu, b, [mu]), self.gauge("G", be, c, [al])]))
        a = self.adj(r)
        w = g[r] ** 2 * g["B"] * eps
        if kind in ("eps_rBf", "eps_4B"):
            b, c = self.adj(r), self.adj(r)
            if kind == "eps_rBf":    # c_r^2 g' eps fg^{abc} (d_mu omega^a) A^b_nu A^c_al B_be
                return self.func((w * g[r] * tn.f(a, b, c), [self.ghost(r, a, [mu]), self.gauge(r, nu, b),
                                                              self.gauge(r, al, c), self.gauge("B", be, None)]))
            return self.func((w * g[r] * tn.f(a, b, c), [self.gauge(r, mu, a), self.gauge(r, nu, b),
                                                         self.gauge(r, al, c), self.ghost("B", None, [be])]))
        if kind == "eps_Brr":    # c_r^2 g' eps gh^B d A^a d A^a
            return self.func((w, [self.ghost("B", None), self.gauge(r, nu, a, [mu]), self.gauge(r, be, a, [al])]))
        if kind == "eps_rB":     # c_r^2 g' eps gh^a d A^a d B
            return self.func((w, [self.ghost(r, a), self.gauge(r, nu, a, [mu]), self.gauge("B", be, None, [al])]))
        raise KeyError(kind)

    def n_all(self) -> dict[str, LocalFunctional]:
        """Every N monomial with a non-zero coefficient in X_1 (engine names of final_counterterms.json)."""
        base = self.n_basis()
        out = {}
        for r in "GWB":
            out[f"cV3_{r}"] = base[f"N3_{r}"]
        for r in "GW":
            out[f"f1_{r}"], out[f"f2_{r}"], out[f"ff_{r}"] = base[f"Nf1_{r}"], base[f"Nf2_{r}"], base[f"Nff_{r}"]
        for n in ("G4a", "G4b", "W4a", "B4"):
            out[n] = base["N" + n]
        for rs in ("GW", "GB", "WB"):
            for x in "abcd":
                out[f"{rs}{x}"] = base[f"N{rs}{x}"]
        for k in (1, 2, 3):
            out[f"GGGB{k}"] = self.GGGB(k)
        for f in FERMIONS:
            for r in FACTORS_OF[f]:
                out[f"F_{f}{r}"] = self.NF(f, r)
                if r in "GW":
                    out[f"FV_{f}{r}{r}"] = self.FV(f, r, r, True)
                    out[f"VF_{f}{r}{r}"] = self.FV(f, r, r, False)
                    out[f"RR_{f}{r}{r}"] = self.RR(f, r, r, False)
                    out[f"RRb_{f}{r}{r}"] = self.RR(f, r, r, True)
        for r in "WB":
            out[f"S_{r}"] = base[f"NS_{r}"]
            for s in "WB":
                out[f"SV_{r}{s}"] = base[f"NSV_{r}{s}"]
        out["SG"] = base["NSG"]
        for k in "ude":
            out[f"Y{k}L_W"] = self.Y(k, False, "L", "W")
            out[f"Y{k}bR_W"] = self.Y(k, True, "R", "W")
        return out

    # ------------------------------------------------------------------ evanescent operators (eq:evanescent)
    def E_A(self, r):
        mu, nu, rho, sig, a = self.lor(), self.lor(), self.lor(), self.lor(), self.adj(r)
        coeff = num(Fraction(1, 2)) * tn.gbar(mu, nu) * tn.ghat(rho, sig)
        return LocalFunctional.integrate(Density.monomial(self.F, coeff, self.gauge(r, mu, a),
                                                          self.gauge(r, nu, a, [rho, sig])))

    def E_H(self):
        row, col, mu, nu = self.matter_labels("Hdag"), self.matter_labels("H"), self.lor(), self.lor()
        coeff = tn.ghat(mu, nu) * self.identity("H", row, col)
        return LocalFunctional.integrate(Density.monomial(self.F, coeff, _occ("Hdag", row, [mu]),
                                                          _occ("H", col, [nu])))

    def tilde(self, name):
        """``P(K) - K`` for the kinetic invariant ``name`` (``F2[r]`` or ``DS2[H]``): Etilde_F^(r), Etilde_H."""
        simp = self.m.simplifier()
        k = dict(kinetic_invariants(self.m))[name]
        return LocalFunctional.integrate(k.physical(simp) - k)

    # ------------------------------------------------------------------ golden functionals
    def closed_d(self, q):
        """The 25 closed directions continued to d dimensions (report sec:evanescent, 'Continuation')."""
        out = {}
        for name, f in self.closed(q).items():
            e = BasisElement(name, f, (), (), 0, "monomial", name)
            out[name] = continuation(e).functional
        return out

    def counterterm(self, doc, q) -> LocalFunctional:
        """``C`` with ``deltaS^(1) = C / (16 pi^2 eps)``: closed directions (d-dimensional) with the coefficients of
        ``delta_S.closed_direction_coefficients_4d`` plus the evanescent operators (``delta_S.evanescent.SM``)."""
        names = {"O_A^G": "O_A_G", "O_c^G": "O_c_G", "O_g^G": "O_g_G", "O_A^W": "O_A_W", "O_c^W": "O_c_W",
                 "O_g^W": "O_g_W", "Q_T^B": "T_B", "O_A^B": "O_A_B", "O_c^B": "O_c_B", "O_H": "O_H",
                 "Q_lambda": "lam", "Q_m": "m", "O_baryon": "O_baryon", "O_lepton": "O_lepton"}
        for f in FERMIONS:
            names[f"O_psi^{f}"] = f"O_psi_{f}"
        for k in "ude":
            names[f"Q^+_y{k}"], names[f"Q^-_y{k}"] = f"y{k}+", f"y{k}-"
        closed = self.closed_d(q)
        total = Density(self.F)
        for key, value in doc["delta_S"]["closed_direction_coefficients_4d"].items():
            c = parse(self.m, value)
            if c != Expression.num(0):
                total = total + closed[names[key]].scale(c)
        return LocalFunctional.integrate(total + self.evanescent(doc))

    def evanescent(self, doc) -> LocalFunctional:
        ev = doc["delta_S"]["evanescent"]["SM"]
        total = Density(self.F)
        for r in "GWB":
            total = total + self.E_A(r).scale(parse(self.m, ev[f"E_A^({r})"]))
            total = total + self.tilde(f"F2[{r}]").scale(parse(self.m, ev[f"Etilde_F^({r})"]))
        total = total + self.E_H().scale(parse(self.m, ev["E_H"]))
        total = total + self.tilde("DS2[H]").scale(parse(self.m, ev["Etilde_H"]))
        return LocalFunctional.integrate(total)

    def breaking(self, doc) -> LocalFunctional:
        """``X_1 = sum X[N] N`` (``X.engine_names``)."""
        ns = self.n_all()
        total = Density(self.F)
        for name, value in doc["X"]["engine_names"].items():
            total = total + ns[name].scale(parse(self.m, value))
        return LocalFunctional.integrate(total)

    def restoring(self, doc, q, key=("r", "all_97_engine_names")) -> LocalFunctional:
        """``DeltaS^(1)_ST = sum r[Q] Q`` (``r.all_97_engine_names``)."""
        values = doc
        for k in key:
            values = values[k]
        total = Density(self.F)
        for name, value in values.items():
            c = parse(self.m, value)
            if c != Expression.num(0):
                total = total + q[name].scale(c)
        return LocalFunctional.integrate(total)


# ---------------------------------------------------------------------------------------------- report 05a
def sm_invariants(model) -> dict:
    """The group and Yukawa invariants of report 04z ('Symbols used throughout') for the one-family SM."""
    gs, g, gp = model.param("g_s"), model.param("g"), model.param("gp")
    yu, yd, ye, lam = (model.param(n) for n in ("y_u", "y_d", "y_e", "lambda"))
    q = lambda a, b: num(Fraction(a, b))
    return {
        "S": {"G": 2 * gs**2, "W": 2 * g**2, "B": q(10, 3) * gp**2},
        "S_H": {"G": Expression.num(0), "W": g**2 / 2, "B": gp**2 / 2},
        "C_ad": {"G": 3 * gs**2, "W": 2 * g**2, "B": Expression.num(0)},
        "C_H": q(3, 4) * g**2 + q(1, 4) * gp**2,
        "Yt": 3 * yu**2 + 3 * yd**2 + ye**2,
        "Y4": 3 * yu**4 + 3 * yd**4 + ye**4,
        "G4": 3 * g**4 + 2 * g**2 * gp**2 + gp**4,
        "C": {"Q": q(4, 3) * gs**2 + q(3, 4) * g**2 + q(1, 36) * gp**2, "U": q(4, 3) * gs**2 + q(4, 9) * gp**2,
              "D": q(4, 3) * gs**2 + q(1, 9) * gp**2, "L": q(3, 4) * g**2 + q(1, 4) * gp**2, "E": gp**2},
        "Y": {"Q": (yu**2 + yd**2) / 2, "U": yu**2, "D": yd**2, "L": ye**2 / 2, "E": ye**2},
        "nu": {("Q", "G"): yu**2 + yd**2, ("Q", "W"): Expression.num(0), ("Q", "B"): 4 * yu**2 - 2 * yd**2,
               ("U", "G"): 2 * yu**2, ("U", "B"): yu**2 / 2, ("D", "G"): 2 * yd**2, ("D", "B"): -yd**2,
               ("L", "W"): Expression.num(0), ("L", "B"): 2 * ye**2, ("E", "B"): ye**2},
        "yukawa_z": {"y_u": -(q(16, 3) * gs**2 + q(3, 4) * g**2 + q(25, 36) * gp**2 + yd**2),
                     "y_d": -(q(16, 3) * gs**2 + q(3, 4) * g**2 + q(1, 36) * gp**2 + yu**2),
                     "y_e": -(q(3, 4) * g**2 + q(9, 4) * gp**2)},
        "lambda": lam,
    }


class CounterRules(Golden):
    """The counterterm Feynman rules of report 05a as local functionals per leg class (the rule of a class is
    ``i V`` of its functional): cr-tab-inv (``z`` times the tree rule), cr-tab-full and cr-ev-rules (unsplit,
    barred and hatted parts), and the finite counterterm eq:sts-result / eq:c-DS (gauge part in the trace form
    cr-fin-4 and the first line of eq:c-DS)."""

    def __init__(self, model, action) -> None:
        super().__init__(model)
        self.action = action
        self.simp = model.simplifier()
        self.inv = sm_invariants(model)

    def tree(self, cls) -> LocalFunctional:
        key = sorted(cls)
        return LocalFunctional.integrate(Density.from_monomials(
            self.F, [mo for mo in self.action.monomials() if sorted(mo.fields()) == key]))

    def barred(self, f: LocalFunctional) -> LocalFunctional:
        return LocalFunctional.integrate(f.physical(self.simp))

    def _lin(self, *terms) -> LocalFunctional:
        total = Density(self.F)
        for c, f in terms:
            total = total + f.scale(c.expand() if isinstance(c, Expression) else Expression.num(c))
        return LocalFunctional.integrate(total)

    def singular(self) -> dict[tuple[str, ...], tuple[LocalFunctional, bool]]:
        """``{class: (C_class, physical)}`` with ``deltaS^(1) = C/(16 pi^2 eps)``; ``physical``: compare on
        physical external spinors (the fermion kinetic rule contains only the physical components)."""
        inv, order = self.inv, self.F.order
        canon = lambda *c: tuple(sorted(c, key=order))
        out: dict[tuple[str, ...], tuple[LocalFunctional, bool]] = {}
        dZf = {f: -(inv["C"][f] + inv["Y"][f]) for f in FERMIONS}
        for r in "GWB":
            S, SH, Ca = inv["S"][r], inv["S_H"][r], inv["C_ad"][r]
            aa = canon(r, r)
            t = self.tree(aa)
            out[aa] = (self._lin((num(Fraction(5, 3)) * Ca - SH / 3, t), (-2 * S / 3, self.barred(t)),
                                 (-S / 3, self.E_A(r))), False)
            half = Ca / 2
            src, gh, agh = self.F[r].source, f"gh_{r}", f"agh_{r}"
            for cls, z in ((canon(agh, gh), half), (canon(src, gh), half), (canon(agh, r, gh), -half),
                           (canon(src, r, gh), -half), (canon(f"src_{gh}", gh, gh), -half)):
                t = self.tree(cls)
                if not t.is_zero():
                    out[cls] = (self._lin((z, t)), False)
            if r != "B":
                for n, (a, b) in ((3, (num(Fraction(2, 3)) * Ca - SH / 3, -2 * S / 3)),
                                  (4, (-(Ca + SH) / 3, -2 * S / 3))):
                    cls = canon(*([r] * n))
                    t = self.tree(cls)
                    out[cls] = (self._lin((a, t), (b, self.barred(t))), False)
            for f in FERMIONS:
                if r not in FACTORS_OF[f]:
                    continue
                cls = canon(f + "bar", f, r)
                out[cls] = (self._lin((dZf[f] - Ca, self.tree(cls))), False)
                for s_, x in ((f"src_{f}", f), (f + "bar", f"src_{f}bar")):
                    cls = canon(s_, gh, x)
                    out[cls] = (self._lin((-half, self.tree(cls))), False)
            cls = canon("Hdag", "H", r)
            t = self.tree(cls)
            out[cls] = (self._lin((2 * inv["C_H"] - Ca, t), (-inv["Yt"], self.barred(t))), False)
            if r != "G":
                for cls in (canon("src_H", gh, "H"), canon("Hdag", gh, "src_Hdag")):
                    out[cls] = (self._lin((-half, self.tree(cls))), False)
        for r, s in (("W", "W"), ("W", "B"), ("B", "B")):
            cls = canon("Hdag", "H", r, s)
            t = self.tree(cls)
            out[cls] = (self._lin((2 * inv["C_H"] - inv["C_ad"][r] - inv["C_ad"][s], t),
                                  (-inv["Yt"], self.barred(t))), False)
        for f in FERMIONS:
            out[canon(f + "bar", f)] = (self._lin((dZf[f], self.tree(canon(f + "bar", f)))), True)
        hh = self.tree(canon("Hdag", "H"))
        m2 = self.m.param("m2")
        kin = LocalFunctional.integrate(Density.from_monomials(
            self.F, [mo for mo in hh.monomials() if any(o.derivatives for o in mo.occurrences)]))
        lam = inv["lambda"]
        out[canon("Hdag", "H")] = (self._lin((2 * inv["C_H"], kin), (-inv["Yt"], self.barred(kin)),
                                             (-2 * inv["Yt"] / 3, self.E_H()),
                                             (-m2 * (6 * lam - inv["C_H"]), self.higgs("m"))), False)
        cls = canon("Hdag", "H", "Hdag", "H")
        zl = (12 * lam**2 - 2 * lam * inv["C_H"] - inv["Y4"] + num(Fraction(3, 16)) * inv["G4"]) / lam
        out[cls] = (self._lin((zl, self.tree(cls))), False)
        for k, fb, s_, f in (("y_u", "Qbar", "Hdag", "U"), ("y_d", "Qbar", "H", "D"), ("y_e", "Lbar", "H", "E")):
            for cls in (canon(fb, s_, f), canon(self.F[fb].conjugate, self.F[s_].conjugate, self.F[f].conjugate)):
                out[cls] = (self._lin((inv["yukawa_z"][k], self.tree(cls))), False)
        return out

    # ------------------------------------------------------------------ finite counterterm, eq:c-DS
    def trace(self, f, legs):
        """``Tr_f(Am_1 .. Am_n)`` for gauge legs ``[(factor, mu, derivatives), ...]``: the sum over the factors
        acting on ``f`` is done by the caller; returns ``(coefficient, occurrences)``."""
        labels = [self.matter_labels(f) for _ in legs]
        coeff = Expression.num(1)
        occs = []
        for k, (r, mu, ders) in enumerate(legs):
            a = self.adj(r)
            coeff = coeff * self.generator(f, r, a, labels[k], labels[(k + 1) % len(legs)])
            occs.append(self.gauge(r, mu, a, ders))
        return coeff, occs

    def finite_gauge(self) -> LocalFunctional:
        """``sum_f Tr_f[1/6 (d.Am)^2 - i/3 (d_mu Am_nu)[Am^mu, Am^nu] + 1/3 (Am.Am)^2 - 1/4 Am_mu Am_nu Am^mu
        Am^nu]`` (four-dimensional: ``gbar``)."""
        import itertools as it

        parts = []
        for f in FERMIONS:
            facs = FACTORS_OF[f]
            for rs in it.product(facs, repeat=2):
                mu, nu = self.lor(), self.lor()
                c, o = self.trace(f, [(rs[0], mu, [mu]), (rs[1], nu, [nu])])
                parts.append((num(Fraction(1, 6)) * c, o))
            for rs in it.product(facs, repeat=3):
                mu, nu, a1, a2 = self.lor(), self.lor(), self.lor(), self.lor()
                g = tn.gbar(mu, a1) * tn.gbar(nu, a2)
                c1, o1 = self.trace(f, [(rs[0], nu, [mu]), (rs[1], a1, []), (rs[2], a2, [])])
                c2, o2 = self.trace(f, [(rs[0], nu, [mu]), (rs[1], a2, []), (rs[2], a1, [])])
                parts.append((-I / 3 * g * c1, o1))
                parts.append((I / 3 * g * c2, o2))
            for rs in it.product(facs, repeat=4):
                mu, nu = self.lor(), self.lor()
                c, o = self.trace(f, [(rs[0], mu, []), (rs[1], mu, []), (rs[2], nu, []), (rs[3], nu, [])])
                parts.append((num(Fraction(1, 3)) * c, o))
                mu, nu = self.lor(), self.lor()
                c, o = self.trace(f, [(rs[0], mu, []), (rs[1], nu, []), (rs[2], mu, []), (rs[3], nu, [])])
                parts.append((num(Fraction(-1, 4)) * c, o))
        total = Density(self.F)
        for c, o in parts:
            total = total + Density.monomial(self.F, c, *o)
        return LocalFunctional.integrate(total)

    def finite(self) -> LocalFunctional:
        """eq:sts-result: the gauge trace form plus the fermion, source and Higgs terms."""
        inv = self.inv
        yu, yd, ye = (self.m.param(n) for n in ("y_u", "y_d", "y_e"))
        terms = [(Expression.num(1), self.finite_gauge())]
        for f in FERMIONS:
            for r in FACTORS_OF[f]:
                terms.append((-(inv["C"][f] + inv["nu"][(f, r)] / 2), self.V(f, r)))
                if r != "B":
                    terms.append((-inv["C_ad"][r] / 4, self.R(f, r, 1)))
        for r in "WB":
            terms.append((inv["Yt"] / 3, self.higgs(f"J_{r}")))
        terms += [(-inv["Yt"], self.higgs("HWW")), (-(yu**2 + 2 * yd**2), self.higgs("HWB")),
                  (-(num(Fraction(13, 3)) * yu**2 + num(Fraction(7, 3)) * yd**2 + 3 * ye**2), self.higgs("HBB")),
                  (-(yu**2 + yd**2) / 2, self.higgs("HGG"))]
        return self._lin(*terms)


# ---------------------------------------------------------------------------------------------- the report's basis
REMOVED = {
    "N1": ("K_Q", "K_U", "K_D", "K_L", "K_E", "SS", "m", "lam", "yu+", "yu-", "yd+", "yd-", "ye+", "ye-"),
    "N2": ("T_G", "T_W", "T_B"),
    "N3": ("rd_G", "rd_W", "rd_B", "rf_G", "rf_W"),
    "N4": ("rWWB", "rGGB", "rBBB", "zGGB", "zWWB", "R_QB+", "R_UB+", "R_DB+", "R_LB+", "R_EB+", "sig_B+"),
}
"""The monomials removed by (N1)-(N4) in the report's basis (03b sec:normalisation; engine names)."""


def report_basis(model, q) -> OperatorBasis:
    """The 97 monomials of tab:Qbasis as an :class:`OperatorBasis` (no fingerprint layout: for st_solve)."""
    elements = []
    for name in Report.report_names(q):
        f = q[name]
        classes = tuple(external_classes(model, f))
        degrees = {sum(len(o.derivatives) for o in mo.occurrences) for mo in f.monomials()}
        elements.append(BasisElement(name, f, classes[0], classes, degrees.pop() if len(degrees) == 1 else -1,
                                     "monomial", name))
    return OperatorBasis(model, 0, Fraction(4), True, elements, [], (), 4)


def anomaly_breaking(model, golden, D) -> LocalFunctional:
    """eq:sts-anom-basis for anomaly numbers ``D = {"GGG", "BBB", "BGG", "BWW"}`` (the d-tensor term D_GGG is
    omitted: it vanishes for the contents tested here)."""
    q3 = num(Fraction(-1, 3))
    q12 = num(Fraction(1, 12))
    d = {k: num(v) for k, v in D.items()}
    terms = [(q3 * d["BBB"], golden.eps_anomaly("eps_B"))]
    for r in ("G", "W"):
        if D[f"B{r}{r}"]:
            terms += [(q3 * d[f"B{r}{r}"], golden.eps_anomaly("eps_Brr", r)),
                      (q3 * 2 * d[f"B{r}{r}"], golden.eps_anomaly("eps_rB", r)),
                      (q12 * d[f"B{r}{r}"], golden.eps_anomaly("eps_rBf", r)),
                      (-q12 * d[f"B{r}{r}"], golden.eps_anomaly("eps_4B", r))]
    total = Density(model.fields)
    for c, f in terms:
        total = total + f.scale(c)
    return LocalFunctional.integrate(total)
