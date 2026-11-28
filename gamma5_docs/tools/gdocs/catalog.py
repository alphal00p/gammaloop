"""Operator catalogue of the unbroken phase: ``site/data/unbroken.json`` (DATA_FORMAT.md, section 4).

The renormalised one-loop action of the one-family Standard Model in the unbroken phase,
``S_ren = S_0 + deltaS^(1) + DeltaS^(1)``, split into operators with definite field content ("terms"). Everything
is derived from the outputs of the automation (``site/resources/run_sm_one_family``: the order store, the order
data and the run report) and compared with the report's final numbers (``site/resources/final_counterterms.json``).

Method
======

1. Operators. The operators are the report's monomials ``Q_alpha`` of tab:Qbasis (97 monomials, built by the
   package's test helper ``tests/report_bases.py``, which is a hand transcription of tab:Qbasis used read-only),
   restricted to one class of external fields each (``Q^(r)_rho_d``, ``Q_rho[f]``, ``Q^(fr)_R``, ``Q^(r)_sigma`` and
   ``Q^+_(y_k)`` have two field contents and give two terms each); the evanescent operators ``E^(r)_A``, ``E_H``,
   ``Etilde^(r)_F``, ``Etilde_H`` of eq:evanescent, split into their two-, three- and four-point parts; the tree-level
   operators outside tab:Qbasis (gauge fixing ``b d.A`` and ``xi/2 b b``, the sterile kinetic term
   ``i chibar gammabar d chi`` and the evanescent fermion kinetic term ``i Psibar gammahat d Psi``); and the ten parts
   ``Delta_1^(fr)`` of the breaking (eq:bd-Oc), listed for their Feynman rules.
   Each operator has a d-dimensional form (tree level, singular counterterm: unsplit metrics in bosonic monomials,
   the barred vertex in fermion monomials; report sec:evanescent "Continuation") and, for the ``Q_alpha``, the
   four-dimensional form ``P Q_alpha`` (all objects barred) that multiplies the finite coefficient ``r_alpha``.
2. Coefficients. For every class of external fields the order-store functionals ``S_0``, ``16 pi^2 eps deltaS1``,
   ``16 pi^2 DeltaS1`` and ``Delta1`` restricted to that class are decomposed exactly on the operators of the class:
   plane-wave fingerprints (``hepkit_renorm.fingerprint``) in ``D = 6`` over ``GF(p)``, both embeddings of ``i``,
   rational reconstruction over several primes (``hepkit_renorm.linalg``), coefficients Laurent polynomials in the
   couplings. Nothing is typed: the coefficients come from the store.
3. Checks during the build (a mismatch raises): the 39 non-zero ``r_alpha`` equal ``r.engine_names`` of
   final_counterterms.json, the evanescent coefficients equal ``delta_S.evanescent.SM``, the renormalisation
   constants of the run report equal ``Z_factors``, and every singular coefficient of a tree-level operator equals
   its decomposition through the constants (the first-order term of the bare monomial, report eq:c-bare and
   cr-factor). ``tools/tests/test_catalog.py`` repeats the decisive checks from the written file.
4. Feynman rules: every vertex function on the page is ``mathview.vertex`` of (coefficient x operator), i.e.
   ``LocalFunctional.vertex_function`` of the package (the code path of ``hepkit_renorm.rules.vertex_rule``).

Interface
=========

``build_unbroken(site_dir=SITE, *, classes=True, log=print) -> dict``: writes ``data/unbroken.json`` and (with
``classes=True``) fills ``label``, ``result`` and ``summary`` of ``data/classes.json`` (:func:`decorate_classes`).
``python -m gdocs.catalog [--site DIR] [--no-classes]``. Needs ``SYMBOLICA_LICENSE`` in the environment.

``Context(site_dir)`` (model, store, report monomials, run report, final numbers), ``operator_specs(ctx)`` (the
:class:`OpSpec` of every operator, also those with zero coefficients), ``solve_coefficients(ctx, specs)`` and
``catalogue(ctx)`` (the specs with their coefficients, ``[(OpSpec, {part: Expression})]``) are used by the test.

Extra keys of the written file (ignored by the page): ``constants`` (the 20 renormalisation constants, values in
units of ``1/(16 pi^2 eps)``) and, per term, ``data = {"monomial", "report_name", "values"}`` with the report
monomial (engine name of tab:Qbasis / typst name) and the coefficients as Symbolica strings in the names
``g_s g gp y_u y_d y_e lambda m2 xi I`` (``pole``/``evanescent`` in units of ``1/(16 pi^2 eps)``, ``finite`` in
units of ``1/(16 pi^2)``).
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from fractions import Fraction
from functools import cached_property
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
DOCS = HERE.parents[2]
REPO = HERE.parents[3]
SITE = DOCS / "site"
PKG = REPO / "python" / "hepkit_renorm"
RUN = "resources/run_sm_one_family"

for _p in (PKG, PKG / "tests"):
    if str(_p) not in sys.path:
        sys.path.append(str(_p))

from symbolica import AtomType, Expression, S  # noqa: E402

from hepkit_renorm import Density, LocalFunctional, Occurrence, builtin_model  # noqa: E402
from hepkit_renorm import linalg  # noqa: E402
from hepkit_renorm import tensors as tn  # noqa: E402
from hepkit_renorm.bases import BasisElement, continuation  # noqa: E402
from hepkit_renorm.class_results import physical_samples  # noqa: E402
from hepkit_renorm.fingerprint import Fingerprinter  # noqa: E402
from hepkit_renorm.orders import loop_factor  # noqa: E402
from hepkit_renorm.store import OrderStore  # noqa: E402

from . import mathview as mv  # noqa: E402
from . import names as nm  # noqa: E402
from .graphs import class_id  # noqa: E402

__all__ = ["Context", "OpSpec", "build_unbroken", "catalogue", "decorate_classes", "operator_specs",
           "solve_coefficients", "values_string", "parse_value"]

FERMIONS = ("Q", "U", "D", "L", "E")
FACTORS = ("G", "W", "B")
FACTORS_OF = {"Q": "GWB", "U": "GB", "D": "GB", "L": "WB", "E": "B"}
GROUP = {"G": "SU(3)", "W": "SU(2)", "B": "U(1)"}
COUPLING = {"G": "g_s", "W": "g", "B": "gp"}
PARAMS = ("g_s", "g", "gp", "m2", "lambda", "y_u", "y_d", "y_e", "xi")
PARTS = ("tree", "pole", "evanescent", "finite", "breaking")

# ---------------------------------------------------------------------------------------------- report references
# Section numbers of the labelled headings (from the Typst heading structure, checked against the PDF); the pages
# are read from site/resources/report_headings.txt.
SECTIONS = {
    "sec:action": "3.3", "sec:brst": "3.4", "sec:delta": "3.5", "sec:feynman-rules": "3.6", "sec:Qbasis": "3.8.4",
    "sec:closed": "3.8.5", "sec:evanescent": "3.8.6", "sec:Nbasis": "3.8.7", "sec:normalisation": "3.8.9",
    "sec:g0_VV": "4.1.1", "sec:g0_VVV": "4.1.2", "sec:g0_VVVV": "4.1.3", "sec:g0_FF": "4.1.4", "sec:g0_FFV": "4.1.5",
    "sec:g0_SS": "4.1.6", "sec:g0_SSV_SSVV": "4.1.7", "sec:g0_SSSS": "4.1.8", "sec:g0_FFS": "4.1.9",
    "sec:g0_ghost_sources": "4.1.10", "sec:g1": "4.2", "sec:st-solution": "4.3", "sec:st-dS": "4.3.1",
    "sec:st-poles": "4.3.2", "sec:st-X": "4.3.3", "sec:st-solve": "4.3.4", "sec:st-result": "4.3.5",
    "sec:st-lit": "4.3.6", "sec:st-checks": "4.3.7", "sec:Sren1": "5.1", "sec:ct-rules": "5.2", "cr-inv": "5.2.1",
    "cr-ev": "5.2.2", "cr-fin": "5.2.3", "sec:m-oneloop": "2.4", "sec:m-objects": "2.1", "sec:fields": "3.2",
}
# Tables (number, page of the caption in report/one_loop_gamma5.pdf)
TABLES = {"Qbasis": (8, 33), "sts-r": (56, 142), "c-Z": (57, 147), "c-summary": (58, 150), "cr-tab-inv": (59, 152),
          "cr-tab-full": (60, 153), "cr-tab-quartic": (61, 154)}


def _pages(site: Path) -> dict[str, int]:
    out = {}
    for line in (site / "resources" / "report_headings.txt").read_text().splitlines():
        m = re.match(r"\s*page\s+(\d+)\s+level\s+\d+\s+<([^>]+)>", line)
        if m:
            out[m.group(2)] = int(m.group(1))
    return out


# ---------------------------------------------------------------------------------------------- expressions
def _num(q) -> Expression:
    q = Fraction(q)
    return Expression.num(q.numerator) / q.denominator


def _is_zero(e: Expression) -> bool:
    return e.expand() == Expression.num(0)


_PLAIN = {n: S(f"cat_plain::{n}") for n in PARAMS + ("I",)}


def values_string(model, expr: Expression) -> str:
    """``expr`` (package parameter symbols, ``i``) as a Symbolica string in the plain names of :data:`PARAMS`
    and ``I`` (module docstring; :func:`parse_value` inverts it)."""
    e = expr.expand()
    for n in PARAMS:
        e = e.replace(model.param(n), _PLAIN[n])
    e = e.replace(tn.I, _PLAIN["I"])
    return e.format(color_top_level_sum=False, color_builtin_symbols=False, bracket_level_colors=None,
                    max_line_length=None, num_exp_as_superscript=False, multiplication_operator="*",
                    show_namespaces=False)


def parse_value(model, text: str) -> Expression:
    """A coefficient string (``values`` of the written file, or the sympy-style strings of
    final_counterterms.json with ``**``, ``g'``) as an expression in the package's parameter symbols."""
    t = text.replace("g'", "gp").replace("**", "^")
    e = Expression.parse(t)
    for n in PARAMS:
        e = e.replace(S(n), model.param(n))
    return e.replace(S("I"), tn.I).expand()


def restrict(F: LocalFunctional, cls: Sequence[str]) -> LocalFunctional:
    """The part of ``F`` with the field content ``cls``."""
    key = sorted(cls)
    return LocalFunctional.integrate(Density.from_monomials(F.fields, [m for m in F.monomials()
                                                                       if sorted(m.fields()) == key]))


_W = S("cat_wa_", "cat_wb_", "cat_wc_", "cat_we_", "cat_wm_")


def tidy(F: LocalFunctional) -> LocalFunctional:
    """Projector chains: ``P_R gammabar^mu P_L`` written as ``gammabar^mu P_L`` (and L <-> R), ``P P = P`` and
    ``P_L P_R = 0``, repeated to a fixpoint; kept out of the canonical form (which expands the projectors) so that
    the display keeps them. Equal to ``F``."""
    a, b, c, e, m = _W
    ex = F.expr
    while True:
        new = ex
        for p1, p2 in ((tn.PR, tn.PL), (tn.PL, tn.PR)):
            new = new.replace(p1(a, b) * tn.gammabar(m, b, c) * p2(c, e), tn.gammabar(m, a, c) * p2(c, e))
            new = new.replace(p1(a, b) * p1(b, c), p1(a, c))
            new = new.replace(p1(a, b) * p2(b, c), Expression.num(0))
        if new == ex:
            break
        ex = new
    return LocalFunctional(F.fields, ex, _canonical=True)


_TW = S("cat_ta_", "cat_tb_", "cat_ti_", "cat_tj_")


def simplify_groups(F: LocalFunctional) -> LocalFunctional:
    """Traces of two generators of one factor, ``Tr(T^a T^b) = delta^(ab)/2`` (package convention 4), then the
    projector chains of :func:`tidy`. Used for the class results, which are shown as stored by the run."""
    a, b, i, j = _TW
    ex = F.expr
    while True:
        new = ex.replace(tn.T(a, i, j) * tn.T(b, j, i), Expression.num(1) / 2 * tn.delta(a, b))
        if new == ex:
            break
        ex = new
    return tidy(LocalFunctional.integrate(LocalFunctional(F.fields, ex, _canonical=True)))


def _has_spinor(F: LocalFunctional) -> bool:
    fields = F.fields
    return any(s.space.kind == "spinor" for m in F.monomials() for o in m.occurrences for s in fields[o.field].slots)


# ---------------------------------------------------------------------------------------------- context
class Context:
    """Inputs of the catalogue (lazy)."""

    def __init__(self, site_dir: str | Path = SITE, log: Callable[[str], None] = print) -> None:
        self.site = Path(site_dir)
        self.log = log
        self.model = builtin_model("sm_one_family")
        self.fields = self.model.fields
        self.simp = self.model.simplifier()
        run = self.site / RUN
        self.store = OrderStore.from_json(self.model, (run / "sm_one_family_store.json").read_text())
        self.report = json.loads((run / "run_report.json").read_text())
        self.order_data = json.loads((run / "sm_one_family_order_data.json").read_text())["orders"][0]
        self.final = json.loads((self.site / "resources" / "final_counterterms.json").read_text())
        self.pages = _pages(self.site)
        self.lf = loop_factor(1)

    # ------------------------------------------------------------------ store functionals (units without 1/(16 pi^2))
    @cached_property
    def S0(self) -> LocalFunctional:
        return self.store.entry("S0").coefficient(0)

    @cached_property
    def Delta1(self) -> LocalFunctional:
        return self.store.entry("Delta1").coefficient(0)

    @cached_property
    def dS(self) -> LocalFunctional:
        """``16 pi^2 eps deltaS^(1)``."""
        return LocalFunctional.integrate(self.store.entry("deltaS1").coefficient(-1).scale(1 / self.lf))

    @cached_property
    def DS(self) -> LocalFunctional:
        """``16 pi^2 DeltaS^(1)``."""
        return LocalFunctional.integrate(self.store.entry("DeltaS1").coefficient(0).scale(1 / self.lf))

    @cached_property
    def bd_dS(self) -> LocalFunctional:
        """``16 pi^2 eps b_d deltaS^(1)``."""
        return LocalFunctional.integrate(self.store.entry("bd_deltaS1").coefficient(-1).scale(1 / self.lf))

    @cached_property
    def X1(self) -> LocalFunctional:
        """``16 pi^2 X_1`` (order data)."""
        return LocalFunctional(self.fields, Expression.parse(self.order_data["X"]), _canonical=True)

    def role(self, name: str) -> LocalFunctional:
        return {"tree": self.S0, "pole": self.dS, "finite": self.DS, "breaking": self.Delta1}[name]

    # ------------------------------------------------------------------ report monomials
    @cached_property
    def golden(self):
        from report_golden import Golden

        return Golden(self.model)

    @cached_property
    def q(self) -> dict[str, LocalFunctional]:
        """The report's monomials (engine names of tab:Qbasis; the 97 of the report, without the engine's '-'
        variants of Q_R and Q_sigma)."""
        from report_bases import Report

        allq = self.golden.q_basis()
        return {n: allq[n] for n in Report.report_names(allq)}

    def q_d(self, name: str) -> LocalFunctional:
        """d-dimensional continuation (report sec:evanescent): the barred vertex in fermion monomials, unsplit
        metrics in bosonic ones."""
        f = self.q[name]
        if _has_spinor(f):
            return f
        bar = LocalFunctional.integrate(f.physical(self.simp))
        return continuation(BasisElement(name, bar, (), (), 0, "monomial", name)).functional

    def q_4(self, name: str) -> LocalFunctional:
        """The four-dimensional monomial ``P Q_alpha`` (all objects barred)."""
        f = self.q[name]
        if _has_spinor(f):
            return f
        return LocalFunctional.integrate(f.physical(self.simp))

    # ------------------------------------------------------------------ run report numbers
    def constant(self, key: str) -> Expression:
        """A constant of the run report (``Z[X]``, ``d[p]``) in units of ``1/(16 pi^2 eps)``."""
        return parse_value(self.model, self.report["results"]["1"]["constants"][key])

    @cached_property
    def fp6(self) -> Fingerprinter:
        return Fingerprinter(self.model, dimension=6)

    def nonzero_classes(self, F: LocalFunctional, dimension: int = 6) -> list[tuple[str, ...]]:
        fp = self.fp6 if dimension == 6 else Fingerprinter(self.model, dimension=dimension)
        out = []
        for c in F.leg_classes():
            smp = fp.samples(list(c), 3, salt=1)
            if np.any(fp.functional(restrict(F, c), smp) % fp.field.p):
                out.append(tuple(c))
        return out

    @cached_property
    def brst(self):
        from hepkit_renorm.brst import BRST

        return BRST(self.model)

    @cached_property
    def class_ids(self) -> set[str]:
        data = json.loads((self.site / "data" / "classes.json").read_text())
        return {c["id"] for c in data["classes"]}

    def ref(self, label: str, text: str | None = None) -> dict:
        page = self.pages.get(label)
        if page is None:
            raise KeyError(f"no page for {label}")
        return {"label": text or f"Section {SECTIONS[label]}", "page": page}

    @staticmethod
    def table(key: str, what: str = "") -> dict:
        n, page = TABLES[key]
        return {"label": f"Table {n}" + (f" ({what})" if what else ""), "page": page}


# ---------------------------------------------------------------------------------------------- operator specs
@dataclass
class OpSpec:
    """One operator of the catalogue (module docstring)."""

    id: str
    name: str
    cls: tuple[str, ...]
    sector: str
    theory: str
    op: LocalFunctional                  # d-dimensional operator: tree level, singular counterterm, breaking
    op4: LocalFunctional | None = None   # four-dimensional operator of the finite part (report monomials)
    monomial: str | None = None          # engine name of the report monomial (key of report_bases)
    family: str = "Q"                    # "Q", "extra", "evanescent", "breaking"
    renormalised: bool = True            # False: the operator is not rescaled by the bare quantities
    order: tuple = ()
    evname: str | None = None            # name of the evanescent operator ("E_A^(G)", ...)
    partner: str | None = None           # id of the other field-content piece of the same report monomial
    remarks: list[str] = field(default_factory=list)

    @property
    def evanescent(self) -> bool:
        return self.family == "evanescent" or self.id.startswith("evkin-")


# names of the report monomials: (sub, sup) of the symbol Q (and r), and the typst name of final_counterterms.json
def _qname(n: str) -> tuple[str | None, str | None, str]:
    m = re.fullmatch(r"(T|L|M|3f)_([GWB])", n)
    if m:
        sub = {"3f": "3f"}.get(m.group(1), m.group(1))
        return sub, f"({m.group(2)})", f"Q^(({m.group(2)}))_({' '.join(sub)})" if sub == "3f" else \
            f"Q^(({m.group(2)}))_{sub}"
    if n in ("G4a", "G4b", "G4c", "W4a", "W4b", "B4", "GGGB"):
        return n, None, f"Q_({' '.join(n)})"
    m = re.fullmatch(r"(GW|GB|WB)([ab])", n)
    if m:
        return m.group(2), f"({m.group(1)})", f"Q^(({' '.join(m.group(1))}))_{m.group(2)}"
    m = re.fullmatch(r"K_([QUDLE])", n)
    if m:
        return "K", m.group(1), f"Q^{m.group(1)}_K"
    m = re.fullmatch(r"V_([QUDLE])([GWB])", n)
    if m:
        return "V", f"{m.group(1)}{m.group(2)}", f"Q^({m.group(1)} {m.group(2)})_V"
    m = re.fullmatch(r"R_([QUDLE])([GWB])\+", n)
    if m:
        return "R", f"{m.group(1)}{m.group(2)}", f"Q^({m.group(1)} {m.group(2)})_R"
    m = re.fullmatch(r"(J|D)_([WB])", n)
    if m:
        return m.group(1), f"({m.group(2)})", f"Q^(({m.group(2)}))_{m.group(1)}"
    m = re.fullmatch(r"H(WW|BB|WB|GG)", n)
    if m:
        return m.group(1), None, f"Q_({' '.join(m.group(1))})"
    if n in ("SS",):
        return "SS", None, "Q_(S S)"
    if n == "m":
        return "m", None, "Q_m"
    if n == "lam":
        return "λ", None, "Q_lambda"
    m = re.fullmatch(r"y([ude])([+-])", n)
    if m:
        return mv.parse(f"y_{m.group(1)}"), m.group(2), f"Q^{m.group(2)}_(y_{m.group(1)})"
    m = re.fullmatch(r"rd_([GWB])", n)
    if m:
        return "ρ∂", f"({m.group(1)})", f"Q^(({m.group(1)}))_(rho partial)"
    m = re.fullmatch(r"rf_([GW])", n)
    if m:
        return "ρ", None, f"Q_rho [fg_(({m.group(1)}))]"
    m = re.fullmatch(r"zf_([GW])", n)
    if m:
        return "ζ", None, f"Q_zeta [-1/2 fg_(({m.group(1)}))]"
    m = re.fullmatch(r"sig_([WB])\+", n)
    if m:
        return "σ", f"({m.group(1)})", f"Q^(({m.group(1)}))_sigma"
    return n, None, n


def q_symbol(n: str, base: str = "Q") -> Expression:
    sub, sup, _ = _qname(n)
    return mv.symbol(base, sub=sub, sup=sup)


def q_suffix_html(n: str) -> str:
    """Bracket suffix of the ghost-sector monomials Q_rho[t], Q_zeta[t] (plain HTML)."""
    m = re.fullmatch(r"(rf|zf)_([GW])", n)
    if not m:
        return ""
    return f"[𝔣<sub>({m.group(2)})</sub>]" if m.group(1) == "rf" else f"[−½𝔣<sub>({m.group(2)})</sub>]"


# removed by the normalisation conditions (report sec:normalisation; report_golden.REMOVED)
NORMALISED = {
    "N1": ("K_Q", "K_U", "K_D", "K_L", "K_E", "SS", "m", "lam", "yu+", "yu-", "yd+", "yd-", "ye+", "ye-"),
    "N2": ("T_G", "T_W", "T_B"),
    "N3": ("rd_G", "rd_W", "rd_B", "rf_G", "rf_W"),
    "N4": ("rWWB", "rGGB", "rBBB", "zGGB", "zWWB", "R_QB+", "R_UB+", "R_DB+", "R_LB+", "R_EB+", "sig_B+"),
}


def normalisation_of(n: str | None) -> str | None:
    for k, names in NORMALISED.items():
        if n in names:
            return k
    return None


GAUGE_NAME = {"G": "gluon", "W": "SU(2) gauge field", "B": "hypercharge gauge field"}
GAUGE_ADJ = {"G": "gluon", "W": "SU(2)", "B": "hypercharge"}
FERMION_NAME = {"Q": "quark doublet", "U": "up-quark singlet", "D": "down-quark singlet", "L": "lepton doublet",
                "E": "electron singlet"}


def _describe_q(n: str, cls: tuple[str, ...], F) -> tuple[str, str, str, str]:
    """``(id, name, sector, theory)`` of the piece of the report monomial ``n`` with field content ``cls``."""
    m = re.fullmatch(r"(T|L|M|3f)_([GWB])", n)
    if m:
        k, r = m.groups()
        return {"T": (f"kin-{r}", f"Kinetic term of the {GAUGE_NAME[r]}", "gauge", "gauge-kinetic"),
                "L": (f"long-{r}", f"Longitudinal two-point term of the {GAUGE_NAME[r]}", "gauge", "st-gauge"),
                "M": (f"mass-{r}", f"Mass-type two-point term of the {GAUGE_NAME[r]}", "gauge", "st-gauge"),
                "3f": (f"cubic-{r}", {"G": "Three-gluon vertex", "W": "Cubic SU(2) gauge vertex"}.get(r, ""),
                       "gauge", "gauge-self")}[k]
    quartic = {"G4a": ("quartic-G-a", "Four-gluon operator (G·G)²"),
               "G4b": ("quartic-G-b", "Four-gluon operator (Gᴬ·Gᴮ)²"),
               "G4c": ("quartic-G-c", "Four-gluon operator with two d-tensors"),
               "W4a": ("quartic-W-a", "Four-W operator (W·W)²"),
               "W4b": ("quartic-W-b", "Four-W operator (Wᴵ·Wᴶ)²"),
               "B4": ("quartic-B", "Four-B operator (B·B)²"),
               "GWa": ("mixed-GW-a", "Gluon–W quartic operator (G·G)(W·W)"),
               "GWb": ("mixed-GW-b", "Gluon–W quartic operator (Gᴬ·Wᴵ)²"),
               "GBa": ("mixed-GB-a", "Gluon–B quartic operator (G·G)(B·B)"),
               "GBb": ("mixed-GB-b", "Gluon–B quartic operator (Gᴬ·B)²"),
               "WBa": ("mixed-WB-a", "W–B quartic operator (W·W)(B·B)"),
               "WBb": ("mixed-WB-b", "W–B quartic operator (Wᴵ·B)²"),
               "GGGB": ("mixed-GGGB", "Three-gluon–B operator with a d-tensor")}
    if n in quartic:
        i, nm = quartic[n]
        th = "gauge-self" if n in ("G4a", "G4b", "G4c", "W4a", "W4b") else "st-gauge"
        return i, nm, "gauge", th
    m = re.fullmatch(r"K_([QUDLE])", n)
    if m:
        f = m.group(1)
        return f"kin-{f}", f"Kinetic term of the physical {FERMION_NAME[f]}", "fermion", "fermion-kinetic"
    m = re.fullmatch(r"V_([QUDLE])([GWB])", n)
    if m:
        f, r = m.groups()
        return f"vertex-{f}-{r}", f"{FERMION_NAME[f].capitalize()}–{GAUGE_ADJ[r]} vertex", "fermion", "fermion-gauge"
    m = re.fullmatch(r"R_([QUDLE])([GWB])\+", n)
    if m:
        f, r = m.groups()
        th = "st-source" if r != "B" else "brst-sources"
        if f"src_{f}" in cls:
            return f"src-R-{f}-{r}", f"Source term R̄ s{f} ({GAUGE_ADJ[r]}-ghost part)", "source", th
        return f"src-Rbar-{f}-{r}", f"Source term (s{f}̄) R ({GAUGE_ADJ[r]}-ghost part)", "source", th
    m = re.fullmatch(r"(J|D)_([WB])", n)
    if m:
        k, r = m.groups()
        if k == "J":
            return f"higgs-current-{r}", f"Higgs–{GAUGE_ADJ[r]} current term", "scalar", "higgs-kinetic"
        return f"higgs-div-{r}", f"Higgs–{GAUGE_ADJ[r]} divergence term", "scalar", "st-scalar"
    hh = {"HWW": ("higgs-WW", "Higgs–W–W seagull term", "higgs-kinetic"),
          "HBB": ("higgs-BB", "Higgs–B–B seagull term", "higgs-kinetic"),
          "HWB": ("higgs-WB", "Higgs–W–B seagull term", "higgs-kinetic"),
          "HGG": ("higgs-GG", "Higgs–gluon–gluon operator", "st-scalar"),
          "SS": ("higgs-kin", "Higgs kinetic term", "higgs-kinetic"),
          "m": ("higgs-mass", "Higgs mass term", "higgs-potential"),
          "lam": ("higgs-quartic", "Higgs quartic self-coupling", "higgs-potential")}
    if n in hh:
        i, nm, th = hh[n]
        return i, nm, "scalar", th
    m = re.fullmatch(r"y([ude])([+-])", n)
    if m:
        k, sign = m.groups()
        lead = {"u": "Qbar", "d": "Qbar", "e": "Lbar"}[k]
        name = {"u": "Up-quark Yukawa term", "d": "Down-quark Yukawa term", "e": "Electron Yukawa term"}[k]
        suffix = "" if sign == "+" else "-imag"
        if lead in cls:
            return f"yukawa-{k}{suffix}", f"{name} −O_{k} (part of Q⁺_y{k})", "yukawa", "yukawa"
        return f"yukawa-{k}{suffix}-conj", f"{name} −O_{k}† (part of Q⁺_y{k})", "yukawa", "yukawa"
    m = re.fullmatch(r"rd_([GWB])", n)
    if m:
        r = m.group(1)
        if f"agh_{r}" in cls:
            return f"ghost-kin-{r}", f"Ghost kinetic term ({GAUGE_ADJ[r]})", "ghost", "ghost"
        return f"src-rho-{r}", f"Source term ρ ∂ω ({GAUGE_ADJ[r]})", "source", "brst-sources"
    m = re.fullmatch(r"rf_([GW])", n)
    if m:
        r = m.group(1)
        if f"agh_{r}" in cls:
            return f"ghost-gauge-{r}", f"Ghost–{GAUGE_ADJ[r]} vertex", "ghost", "ghost"
        return f"src-rhoA-{r}", f"Source term ρ A ω ({GAUGE_ADJ[r]})", "source", "brst-sources"
    m = re.fullmatch(r"zf_([GW])", n)
    if m:
        r = m.group(1)
        return f"src-zeta-{r}", f"Ghost source term ζ ω ω ({GAUGE_ADJ[r]})", "source", "brst-sources"
    m = re.fullmatch(r"sig_([WB])\+", n)
    if m:
        r = m.group(1)
        if "src_H" in cls:
            return f"src-sigma-{r}", f"Source term σ† sH ({GAUGE_ADJ[r]}-ghost part)", "source", "brst-sources"
        return f"src-sigmabar-{r}", f"Source term (sH†) σ ({GAUGE_ADJ[r]}-ghost part)", "source", "brst-sources"
    # monomials of tab:Qbasis that never occur with a non-zero coefficient
    return f"q-{re.sub(r'[^A-Za-z0-9]', '_', n)}-{'-'.join(cls)}", f"Monomial {n}", "gauge", "st-gauge"


SECTOR_ORDER = ("gauge", "fermion", "scalar", "yukawa", "gauge-fixing", "ghost", "source", "breaking")


_SCALAR_ORDER = ("higgs-kin", "higgs-mass", "ev-EH", "ev-EtH-2", "higgs-quartic", "higgs-current-W", "ev-EtH-3-W",
                 "higgs-div-W", "higgs-current-B", "ev-EtH-3-B", "higgs-div-B", "higgs-WW", "ev-EtH-4-WW", "higgs-WB",
                 "ev-EtH-4-WB", "higgs-BB", "ev-EtH-4-BB", "higgs-GG")


def _order_key(spec_id: str, sector: str, cls) -> tuple:
    """Display order: the sectors of eq:L0, then factor or multiplet, number of legs, the operator before its
    finite-only and evanescent companions."""
    sec = SECTOR_ORDER.index(sector)
    rr = {"G": 0, "W": 1, "B": 2}
    fr = {f: i for i, f in enumerate(FERMIONS)}
    if sector == "gauge":
        m = re.search(r"-([GWB])(?:-|$)", spec_id)
        r = rr[m.group(1)] if m and not spec_id.startswith("mixed") else 3
        fam = 2 if spec_id.startswith("ev-") else (1 if spec_id.startswith(("long", "mass")) else 0)
        return (sec, r, len(cls), fam, spec_id)
    if sector == "fermion":
        m = re.search(r"-([QUDLE])(?:-|$)", spec_id)
        kind = next(i for i, k in enumerate(("kin-", "sterile-", "evkin-", "vertex-")) if spec_id.startswith(k))
        r = rr.get(spec_id[-1], 0) if kind == 3 else 0
        return (sec, fr[m.group(1)], kind, r, spec_id)
    if sector == "scalar":
        return (sec, _SCALAR_ORDER.index(spec_id) if spec_id in _SCALAR_ORDER else 99, spec_id)
    if sector == "yukawa":
        return (sec, "ude".index(spec_id.split("-")[1][0]), spec_id.endswith("conj"), spec_id)
    if sector in ("gauge-fixing", "ghost"):
        return (sec, rr[spec_id[-1]], spec_id)
    if sector == "source":
        r = next((rr[x[-1]] for x in cls if x.startswith("gh_")), 3)
        kind = next(i for i, k in enumerate(("src-rho-", "src-rhoA-", "src-zeta-", "src-R-", "src-Rbar-",
                                             "src-sigma-", "src-sigmabar-")) if spec_id.startswith(k))
        m = re.search(r"-([QUDLE])-", spec_id)
        return (sec, r, kind if kind < 3 else (3 if kind < 5 else 5), fr[m.group(1)] if m else 0, kind, spec_id)
    if sector == "breaking":
        _, f, r = spec_id.split("-")
        return (sec, fr[f], rr[r], spec_id)
    return (sec, spec_id)


def _projectors(chirality: str):
    """``(P_c, P_cbar)`` of a multiplet."""
    return (tn.PL, tn.PR) if chirality == "L" else (tn.PR, tn.PL)


def _extra_ops(ctx: Context) -> list[OpSpec]:
    """Tree-level operators outside tab:Qbasis: gauge fixing, sterile and evanescent kinetic terms."""
    F, g = ctx.fields, ctx.golden
    out = []
    xi = ctx.model.param("xi")
    for r in FACTORS:
        cls = tuple(sorted((r, f"b_{r}"), key=F.order))
        out.append(OpSpec(f"gf-bA-{r}", f"Gauge-fixing term b ∂·A ({GAUGE_ADJ[r]})", cls, "gauge-fixing",
                          "gauge-fixing", restrict(ctx.S0, cls), family="extra", renormalised=False))
        cls = (f"b_{r}", f"b_{r}")
        bb = LocalFunctional.integrate(restrict(ctx.S0, cls).scale(2 / xi))
        out.append(OpSpec(f"gf-bb-{r}", f"Gauge-fixing term ξ/2 b b ({GAUGE_ADJ[r]})", cls, "gauge-fixing",
                          "gauge-fixing", bb, family="extra", renormalised=False))
    for f in FERMIONS:
        pc, pcb = _projectors(F[f].chirality)
        row, col, mu = g.matter_labels(f + "bar"), g.matter_labels(f), g.lor()
        s = tn.fresh(tn.SPINOR)
        ident = g.identity(f, row, col)
        cls = tuple(sorted((f, f + "bar"), key=F.order))
        sterile = Density.monomial(F, tn.I * tn.gammabar(mu, row[0], s) * pcb(s, col[0]) * ident,
                                   Occurrence(f + "bar", tuple(row)), Occurrence(f, tuple(col), (mu,)))
        out.append(OpSpec(f"sterile-{f}", f"Kinetic term of the sterile partner of {f}", cls, "fermion",
                          "fermion-kinetic", LocalFunctional(F, sterile.expr, _canonical=True), family="extra",
                          renormalised=False))
        row, col, mu = g.matter_labels(f + "bar"), g.matter_labels(f), g.lor()
        ident = g.identity(f, row, col)
        ev = Density.monomial(F, tn.I * tn.gammahat(mu, row[0], col[0]) * ident,
                              Occurrence(f + "bar", tuple(row)), Occurrence(f, tuple(col), (mu,)))
        out.append(OpSpec(f"evkin-{f}", f"Evanescent kinetic term of {f}", cls, "fermion", "fermion-kinetic",
                          LocalFunctional.integrate(ev), family="extra", renormalised=False))
    return out


def _evanescent_ops(ctx: Context) -> list[OpSpec]:
    g = ctx.golden
    out = []
    npt = {2: "two-point", 3: "three-point", 4: "four-point"}
    for r in FACTORS:
        E = g.E_A(r)
        for cls in E.leg_classes():
            out.append(OpSpec(f"ev-EA-{r}", f"Evanescent operator E_A of the {GAUGE_NAME[r]}", tuple(cls), "gauge",
                              "ct-evanescent", restrict(E, cls), family="evanescent", evname=f"E_A^({r})"))
        Et = g.tilde(f"F2[{r}]")
        for cls in Et.leg_classes():
            n = len(cls)
            out.append(OpSpec(f"ev-EtF-{r}-{n}", f"Evanescent operator Ẽ_F of the {GAUGE_NAME[r]}, {npt[n]} part",
                              tuple(cls), "gauge", "ct-evanescent", restrict(Et, cls), family="evanescent",
                              evname=f"Etilde_F^({r})"))
    E = g.E_H()
    for cls in E.leg_classes():
        out.append(OpSpec("ev-EH", "Evanescent Higgs operator E_H", tuple(cls), "scalar", "ct-evanescent",
                          restrict(E, cls), family="evanescent", evname="E_H"))
    Et = g.tilde("DS2[H]")
    for cls in Et.leg_classes():
        gauge = [x for x in cls if x in FACTORS]
        if not gauge:
            sid, name = "ev-EtH-2", "Evanescent Higgs operator Ẽ_H, two-point part"
        else:
            sid = f"ev-EtH-{len(cls)}-{''.join(gauge)}"
            name = f"Evanescent Higgs operator Ẽ_H, {npt[len(cls)]} part with {' '.join(gauge)}"
        out.append(OpSpec(sid, name, tuple(cls), "scalar", "ct-evanescent", restrict(Et, cls), family="evanescent",
                          evname="Etilde_H"))
    return out


def _breaking_ops(ctx: Context) -> list[OpSpec]:
    """``Delta_1^(fr)`` (report eq:bd-Oc), built from eq:Delta1 with ``gammahat``; equal to the store's Delta1
    class by class (checked by the decomposition, which must give the coefficient 1)."""
    F, g = ctx.fields, ctx.golden
    out = []
    for f in FERMIONS:
        pc, pcb = _projectors(F[f].chirality)
        for r in FACTORS_OF[f]:
            a = g.adj(r)
            gh = Occurrence(f"gh_{r}", () if a is None else (a,))
            total = Density(F)
            for which in (0, 1):
                row, col, mu = g.matter_labels(f + "bar"), g.matter_labels(f), g.lor()
                s = tn.fresh(tn.SPINOR)
                gen = g.generator(f, r, a, row, col)
                if which == 0:   # omega Psibar T gammahat P_cbar d Psi
                    coeff = gen * tn.gammahat(mu, row[0], s) * pcb(s, col[0])
                    occs = (gh, Occurrence(f + "bar", tuple(row)), Occurrence(f, tuple(col), (mu,)))
                else:            # omega (d Psibar) gammahat P_c T Psi
                    coeff = gen * tn.gammahat(mu, row[0], s) * pc(s, col[0])
                    occs = (gh, Occurrence(f + "bar", tuple(row), (mu,)), Occurrence(f, tuple(col)))
                total = total + Density.monomial(F, coeff, *occs)
            cls = tuple(sorted((f"gh_{r}", f, f + "bar"), key=F.order))
            out.append(OpSpec(f"Delta1-{f}-{r}", f"Breaking Δ₁ of {f}, {GAUGE_ADJ[r]} ghost", cls, "breaking",
                              "breaking", LocalFunctional(F, total.expr, _canonical=True), family="breaking",
                              renormalised=False))
    return out


def operator_specs(ctx: Context) -> list[OpSpec]:
    """Every operator of the catalogue, with zero coefficients included (module docstring, step 1)."""
    specs: list[OpSpec] = []
    for n in ctx.q:
        if re.fullmatch(r"y[ude]-", n):
            # Q^-_(y_k) = -i (O_k - O_k^dag) restricted to one field content is a multiple of the piece of
            # Q^+_(y_k): class by class the Q^+ pieces suffice (the report's coefficients of Q^- are zero)
            continue
        qd, q4 = ctx.q_d(n), ctx.q_4(n)
        classes = [tuple(c) for c in qd.leg_classes()]
        pieces = []
        for cls in classes:
            sid, name, sector, theory = _describe_q(n, cls, ctx.fields)
            pieces.append(OpSpec(sid, name, cls, sector, theory, tidy(restrict(qd, cls)), tidy(restrict(q4, cls)),
                                 monomial=n))
        if len(pieces) == 2:
            pieces[0].partner, pieces[1].partner = pieces[1].id, pieces[0].id
        specs += pieces
    specs += _extra_ops(ctx) + _evanescent_ops(ctx) + _breaking_ops(ctx)
    for s in specs:
        s.order = _order_key(s.id, s.sector, s.cls)
    ids = [s.id for s in specs]
    dup = {i for i in ids if ids.count(i) > 1}
    if dup:
        raise RuntimeError(f"duplicate operator ids {sorted(dup)}")
    return specs


# ---------------------------------------------------------------------------------------------- exact decomposition
def _exponents(key: Expression, syms) -> tuple[int, ...]:
    out = []
    for s in syms:
        k, t = 0, key
        while not _is_zero(t.derivative(s)):
            t = (t / s).expand()
            k += 1
        out.append(k)
    return tuple(out)


def _expand(ctx: Context, F: LocalFunctional) -> dict[tuple[int, ...], LocalFunctional]:
    """``F = sum_e p^e F_e`` with parameter-free ``F_e``."""
    syms = [ctx.model.param(n) for n in PARAMS]
    out: dict[tuple[int, ...], LocalFunctional] = {}
    if F.is_zero():
        return out
    for key, coeff in F.expr.expand().coefficient_list(*syms):
        e = _exponents(key, syms)
        part = LocalFunctional(F.fields, coeff, _canonical=True)
        out[e] = part if e not in out else LocalFunctional.integrate(out[e] + part)
    return out


def _monomial(ctx: Context, e: Sequence[int]) -> Expression:
    t = Expression.num(1)
    for n, k in zip(PARAMS, e):
        t = t * ctx.model.param(n) ** k
    return t


def decompose(ctx: Context, target: LocalFunctional, ops: Sequence[LocalFunctional], cls: Sequence[str], *,
              dimension: int = 6, salt: int = 11) -> list[Expression]:
    """Exact coefficients ``x_i`` (Laurent polynomials in the parameters, complex rationals) with
    ``target = sum_i x_i ops_i`` on the class ``cls`` (module docstring, step 2). Raises ``ValueError`` if the
    operators are dependent or the target is not in their span."""
    oe = []
    for o in ops:
        ex = _expand(ctx, o)
        if len(ex) != 1:
            raise ValueError(f"operator not homogeneous in the parameters on {cls}: {list(ex)}")
        oe.append(next(iter(ex.items())))
    tx = _expand(ctx, target)
    keys = list(tx)
    n = len(ops)
    if not keys:
        return [Expression.num(0)] * n

    def compute(field_, i_sign):
        fp = Fingerprinter(ctx.model, field_, dimension, i_sign=i_sign)
        smp = fp.samples(list(cls), n + 4, salt=salt)
        A = np.stack([fp.functional(f, smp) for _, f in oe], axis=1)
        B = np.stack([fp.functional(tx[k], smp) for k in keys], axis=1)
        if linalg.rank(A, field_) != n:
            raise ValueError(f"operators linearly dependent on {cls}")
        try:
            return linalg.solve(A, B, field_, require_unique=True)
        except linalg.InconsistentSystemError as exc:
            raise ValueError(f"store functional outside the span of the operators on {cls}") from exc

    def both(field_):
        re_, im_ = linalg.split_gaussian(compute(field_, 1), compute(field_, -1), field_)
        return np.stack([re_, im_])

    vals = linalg.reconstruct_rationals(both)
    out = []
    for i in range(n):
        x = Expression.num(0)
        for j, k in enumerate(keys):
            re_, im_ = vals[0][i, j], vals[1][i, j]
            if re_ == 0 and im_ == 0:
                continue
            c = _num(re_) + tn.I * _num(im_)
            x = x + c * _monomial(ctx, tuple(a - b for a, b in zip(k, oe[i][0])))
        out.append(x.expand())
    return out


def solve_coefficients(ctx: Context, specs: Sequence[OpSpec]) -> dict[str, dict[str, Expression]]:
    """``{spec id: {part: coefficient}}`` for every operator (non-zero parts only); raises if a non-zero class of
    a store functional has no operators or is outside their span."""
    by_class: dict[tuple[str, ...], list[OpSpec]] = {}
    for s in specs:
        by_class.setdefault(s.cls, []).append(s)
    targets = {}
    for role in ("tree", "pole", "finite", "breaking"):
        for c in ctx.nonzero_classes(ctx.role(role)):
            targets.setdefault(c, set()).add(role)
    missing = [c for c in targets if c not in by_class]
    if missing:
        raise RuntimeError(f"classes of the store without operators: {missing}")
    out: dict[str, dict[str, Expression]] = {s.id: {} for s in specs}
    deferred = []

    def store(use, xs, role):
        for s, x in zip(use, xs):
            if _is_zero(x):
                continue
            part = "evanescent" if role == "pole" and s.family == "evanescent" else role
            out[s.id][part] = x

    for c, roles in sorted(targets.items(), key=lambda kv: [ctx.fields.order(x) for x in kv[0]]):
        group = by_class[c]
        for role in sorted(roles):
            if role == "finite":
                use = [s for s in group if s.op4 is not None]
                ops = [s.op4 for s in use]
            elif role == "breaking":
                use = [s for s in group if s.family == "breaking"]
                ops = [s.op for s in use]
            elif role == "tree":
                use = [s for s in group if s.family in ("Q", "extra")]
                ops = [s.op for s in use]
            else:
                use = [s for s in group if s.family != "breaking"]
                ops = [s.op for s in use]
            if not use:
                raise RuntimeError(f"no operator for the {role} part of the class {c}")
            try:
                xs = decompose(ctx, restrict(ctx.role(role), c), ops, c)
            except ValueError as exc:
                if role != "pole" or "dependent" not in str(exc):
                    raise
                deferred.append((c, use))
                continue
            store(use, xs, role)
        ctx.log(f"  class {' '.join(c)}: {', '.join(sorted(roles))}")
    # classes in which pieces of different evanescent operators coincide (the two-point parts of E_H and
    # Etilde_H are both proportional to (dhat H)^dag dhat H): an evanescent operator has one coefficient for all
    # its parts, so the coefficient found in its other classes is used there and the rest is solved
    for c, use in deferred:
        known = {}
        for s in specs:
            if s.evname and "evanescent" in out[s.id] and s.cls != c:
                known.setdefault(s.evname, set()).add(out[s.id]["evanescent"].to_canonical_string())
        fixed = [s for s in use if s.evname in known]
        if any(len(known[s.evname]) != 1 for s in fixed):
            raise RuntimeError(f"evanescent operator with different coefficients in its parts: {known}")
        target = restrict(ctx.dS, c)
        for s in fixed:
            coeff = Expression.parse(next(iter(known[s.evname])))
            target = LocalFunctional.integrate(target - s.op.scale(coeff))
            out[s.id]["evanescent"] = coeff
        rest = [s for s in use if s not in fixed]
        store(rest, decompose(ctx, target, [s.op for s in rest], c), "pole")
        ctx.log(f"  class {' '.join(c)}: pole (with the coefficients of {sorted({s.evname for s in fixed})} "
                "from their other parts)")
    return out


# ---------------------------------------------------------------------------------------------- bare-quantity rule
def bare_factor(ctx: Context, spec: OpSpec, tree: Expression) -> tuple[Expression, Expression]:
    """The first-order factor of report cr-factor for ``tree x operator``: ``(site expression, value)``, value in
    units of ``1/(16 pi^2 eps)``. Legs: ``+1/2 dZ`` per field, ``-1/2 dZ`` per antighost, auxiliary field and source;
    couplings: ``dZ_(c_r)`` per power of ``c_r``, ``delta p/p`` per power of ``y_k``, ``lambda``, ``m^2``."""
    F = ctx.fields
    site, val = Expression.num(0), Expression.num(0)
    half = _num(Fraction(1, 2))
    for x in spec.cls:
        info = F[x]
        kind = info.kind
        if kind == "source":
            base = info.source_of
            sign = -1
        elif kind in ("antighost", "auxiliary"):
            base, sign = info.factor, -1
        else:
            base, sign = x, 1
        bi = F[base]
        if bi.kind == "gauge" or kind in ("antighost", "auxiliary"):
            r = bi.factor
            site += sign * half * mv.parse(f"dZ(A, {r})")
            val += sign * half * ctx.constant(f"Z[{r}]")
        elif bi.kind == "ghost":
            r = bi.factor
            site += sign * half * mv.parse(f"dZ(omega, {r})")
            val += sign * half * ctx.constant(f"Z[gh_{r}]")
        elif bi.kind in ("fermion", "fermion_bar"):
            f = bi.multiplet
            site += sign * half * mv.parse(f"dZ({f})")
            val += sign * half * ctx.constant(f"Z[{f}]")
        elif bi.kind in ("scalar", "scalar_conj"):
            site += sign * half * mv.parse("dZ(H)")
            val += sign * half * ctx.constant("Z[H]")
        else:
            raise ValueError(f"no renormalisation constant for {x}")
    total = _expand(ctx, LocalFunctional.integrate(spec.op.scale(tree)))
    if len(total) != 1:
        raise ValueError(f"{spec.id}: tree term not homogeneous")
    expo = dict(zip(PARAMS, next(iter(total))))
    for p, k in expo.items():
        if not k:
            continue
        if p in ("g_s", "g", "gp"):
            site += k * mv.parse(f"dZc({p})")
            val += k * ctx.constant(f"d[{p}]") / ctx.model.param(p)
        elif p in ("y_u", "y_d", "y_e", "lambda", "m2"):
            site += k * mv.parse(f"dpar({p})") / mv.parse(p)
            val += k * ctx.constant(f"d[{p}]") / ctx.model.param(p)
        else:
            raise ValueError(f"{spec.id}: parameter {p} in a renormalised term")
    return site, val.expand()


def decomposition(ctx: Context, spec: OpSpec, tree: Expression) -> tuple[Expression, Expression]:
    """``(site expression, value)`` of the singular coefficient through the constants: ``tree x bare factor``
    (expanded, unless the tree coefficient has a negative power of a coupling, as for the quartic monomials)."""
    site, val = bare_factor(ctx, spec, tree)
    t_site = mv.coefficient(tree)
    if "^(-" in tree.expand().to_canonical_string():
        expr = t_site * mv.parse("par(x)").replace(mv.parse("x"), site)
    else:
        expr = (t_site * site).expand()
    return expr, (tree * val).expand()


# ---------------------------------------------------------------------------------------------- the catalogue
def catalogue(ctx: Context, specs: Sequence[OpSpec] | None = None) -> list[tuple[OpSpec, dict[str, Expression]]]:
    """The operators with a non-zero coefficient and their coefficients, in display order."""
    specs = list(specs) if specs is not None else operator_specs(ctx)
    coeffs = solve_coefficients(ctx, specs)
    out = [(s, coeffs[s.id]) for s in specs if coeffs[s.id]]
    out.sort(key=lambda sc: sc[0].order)
    return out


# ---------------------------------------------------------------------------------------------- checks (report)
def _final_r(ctx: Context) -> dict[str, tuple[str, Expression]]:
    """``{engine name: (typst name, value)}`` of the 39 non-zero r_alpha of final_counterterms.json."""
    eng = ctx.final["r"]["engine_names"]
    typ = ctx.final["r"]["typst_names"]
    if len(eng) != len(typ):
        raise RuntimeError("final_counterterms.json: engine and typst name lists differ in length")
    out = {}
    for (e, v), (t, w) in zip(eng.items(), typ.items()):
        a, b = parse_value(ctx.model, v), parse_value(ctx.model, w)
        if not _is_zero(a - b):
            raise RuntimeError(f"final_counterterms.json: {e} and {t} differ")
        out[e] = (t, a)
    return out


EV_KEYS = {"E_A^(G)": "E_A^(G)", "E_A^(W)": "E_A^(W)", "E_A^(B)": "E_A^(B)", "Etilde_F^(G)": "Etilde_F^(G)",
           "Etilde_F^(W)": "Etilde_F^(W)", "Etilde_F^(B)": "Etilde_F^(B)", "E_H": "E_H", "Etilde_H": "Etilde_H"}

CONSTANTS = [  # (name, site source, run-report key, final_counterterms Z_factors path)
    ("dZ_A^(G)", "dZ(A, G)", "Z[G]", ("dZ_A", "G")), ("dZ_A^(W)", "dZ(A, W)", "Z[W]", ("dZ_A", "W")),
    ("dZ_A^(B)", "dZ(A, B)", "Z[B]", ("dZ_A", "B")),
    ("dZ_omega^(G)", "dZ(omega, G)", "Z[gh_G]", ("dZ_omega", "G")),
    ("dZ_omega^(W)", "dZ(omega, W)", "Z[gh_W]", ("dZ_omega", "W")),
    ("dZ_omega^(B)", "dZ(omega, B)", "Z[gh_B]", ("dZ_omega", "B")),
    ("dZ_g_s", "dZc(g_s)", "d[g_s]", ("dZ_g", "G")), ("dZ_g", "dZc(g)", "d[g]", ("dZ_g", "W")),
    ("dZ_gp", "dZc(gp)", "d[gp]", ("dZ_g", "B")),
    ("dZ_Q", "dZ(Q)", "Z[Q]", ("dZ_f", "Q")), ("dZ_U", "dZ(U)", "Z[U]", ("dZ_f", "U")),
    ("dZ_D", "dZ(D)", "Z[D]", ("dZ_f", "D")), ("dZ_L", "dZ(L)", "Z[L]", ("dZ_f", "L")),
    ("dZ_E", "dZ(E)", "Z[E]", ("dZ_f", "E")), ("dZ_H", "dZ(H)", "Z[H]", ("dZ_H",)),
    ("dy_u", "dpar(y_u)", "d[y_u]", ("dy", "u")), ("dy_d", "dpar(y_d)", "d[y_d]", ("dy", "d")),
    ("dy_e", "dpar(y_e)", "d[y_e]", ("dy", "e")), ("dlambda", "dpar(lambda)", "d[lambda]", ("dlambda",)),
    ("dm2", "dpar(m2)", "d[m2]", ("dm2",)),
]


def constants(ctx: Context) -> list[tuple[str, str, Expression]]:
    """The 20 renormalisation constants ``[(name, site source, value)]`` (units ``1/(16 pi^2 eps)``; the coupling
    constants ``dZ_(c_r) = delta c_r / c_r``), from the run report and checked against final_counterterms.json."""
    z = ctx.final["delta_S"]["Z_factors"]["values_in_kappa_units"]
    out = []
    for name, src, key, path in CONSTANTS:
        val = ctx.constant(key)
        if key.startswith("d[g"):
            val = (val / ctx.model.param(key[2:-1])).expand()
        ref = z
        for k in path:
            ref = ref[k]
        if not _is_zero(val - parse_value(ctx.model, ref)):
            raise RuntimeError(f"renormalisation constant {name}: run report {val} != report {ref}")
        out.append((name, src, val))
    return out


def check_against_report(ctx: Context, cl) -> dict:
    """Build-time checks of the module docstring, step 3 (raises on a mismatch); returns counts."""
    final_r = _final_r(ctx)
    seen = {}
    for s, c in cl:
        if "finite" not in c:
            continue
        if s.monomial not in final_r:
            raise RuntimeError(f"{s.id}: finite coefficient for {s.monomial}, which has r = 0 in the report")
        name, val = final_r[s.monomial]
        if not _is_zero(c["finite"] - val):
            raise RuntimeError(f"{s.id}: r = {c['finite']} differs from the report's {name} = {val}")
        seen[s.monomial] = name
    missing = sorted(set(final_r) - set(seen))
    if missing:
        raise RuntimeError(f"r_alpha of the report without a term: {missing}")
    ev = ctx.final["delta_S"]["evanescent"]["SM"]
    evseen = {}
    for s, c in cl:
        if s.family != "evanescent":
            continue
        val = parse_value(ctx.model, ev[s.evname])
        if not _is_zero(c.get("evanescent", Expression.num(0)) - val):
            raise RuntimeError(f"{s.id}: evanescent coefficient {c.get('evanescent')} != report {val}")
        evseen.setdefault(s.evname, 0)
        evseen[s.evname] += 1
    if set(evseen) != set(ev):
        raise RuntimeError(f"evanescent operators: {sorted(evseen)} != {sorted(ev)}")
    consts = constants(ctx)
    ndec = 0
    for s, c in cl:
        if s.family not in ("Q", "extra") or "tree" not in c:
            continue
        if not s.renormalised:
            if "pole" in c:
                raise RuntimeError(f"{s.id}: pole on an operator that is not renormalised")
            continue
        _, val = decomposition(ctx, s, c["tree"])
        if not _is_zero(val - c.get("pole", Expression.num(0))):
            raise RuntimeError(f"{s.id}: pole {c.get('pole')} != bare-quantity decomposition {val}")
        ndec += 1
    return {"r": len(seen), "evanescent_operators": len(evseen), "evanescent_terms": sum(evseen.values()),
            "constants": len(consts), "decompositions": ndec}


# ---------------------------------------------------------------------------------------------- classes of a term
def _leg_rank(F, x: str) -> int:
    info = F[x]
    k = info.kind
    if k == "antighost":
        return 0
    if k == "source":
        base = F[info.source_of].kind
        if base in ("gauge", "ghost", "fermion", "scalar"):
            return 0          # rho, zeta, Rbar_f (src_f), sigma^dag (src_H): row objects, first
        return 8              # R_f (src_fbar), sigma (src_Hdag): column objects, last
    return {"ghost": 1, "fermion_bar": 2, "fermion": 3, "scalar_conj": 4, "scalar": 5, "gauge": 6,
            "auxiliary": 7}[k]


def rule_legs(ctx: Context, spec: OpSpec) -> tuple[str, ...]:
    """Leg order of the displayed vertex function: the report's (Psibar before Psi, the source or antighost of a
    row first, the ghost before Psibar or H^dag), sec:feynman-rules."""
    return tuple(sorted(spec.cls, key=lambda x: _leg_rank(ctx.fields, x)))


def pole_class(ctx: Context, spec: OpSpec) -> str | None:
    """The ghost-number-0 class whose graphs fix the singular coefficient (antighost -> gauge source: the
    effective action depends on omegabar only through rhotilde, eq:rhotilde)."""
    F = ctx.fields
    legs = []
    for x in spec.cls:
        if F[x].kind == "antighost":
            legs.append(F[F[x].factor].source)
        else:
            legs.append(x)
    cid = class_id(legs, 0)
    return cid if cid in ctx.class_ids else None


def breaking_classes(ctx: Context, monomial: str, cache: dict) -> list[str]:
    """Ghost-number-1 classes on which ``b P Q_alpha`` has components (D = 4, physical spinors)."""
    if monomial in cache:
        return cache[monomial]
    fp = Fingerprinter(ctx.model, dimension=4)
    image = LocalFunctional.integrate(ctx.brst.b(ctx.q_4(monomial)))
    out = []
    for c in image.leg_classes():
        part = restrict(image, c)
        smp = fp.samples(list(c), 3, salt=17)
        if any(ctx.fields[x].kind in ("fermion", "fermion_bar") for x in c):
            smp = physical_samples(fp, smp)
        if np.any(fp.functional(part, smp) % fp.field.p):
            # antighost legs: the class with the gauge source (rhotilde = rho + d omegabar, eq:rhotilde)
            legs = [ctx.fields[ctx.fields[x].factor].source if ctx.fields[x].kind == "antighost" else x for x in c]
            cid = class_id(legs, 1)
            if cid not in ctx.class_ids:
                raise RuntimeError(f"b {monomial} has a component on {c}, which is not a class of classes.json")
            if cid not in out:
                out.append(cid)
    cache[monomial] = out
    return out


# ---------------------------------------------------------------------------------------------- math objects
_LONG_HTML = 2500  # static MathML longer than this (characters) gets an interactive view


def _M(e, *, interactive: bool | None = None, display: bool = True, threshold: int = _LONG_HTML,
       budget: int | None = None) -> dict:
    """``mathview.M``; expressions whose static MathML is longer than ``threshold`` characters get the interactive
    view of the Symbolica pager (``interactive=None``), the others are static."""
    out = mv.M(e, interactive=False, display=display)
    if interactive is None:
        interactive = len(out["html"]) > threshold
    if interactive:
        out = mv.M(e, interactive=True, display=display, budget=budget or mv.VIEW_BUDGET)
    if "⋯" in out["html"]:
        raise RuntimeError(f"Symbolica omitted part of an expression (fragment too large): {out['text'][:120]}")
    return out


LOOP = {"tree": None, "pole": "pole", "evanescent": "pole", "finite": "finite", "breaking": None}


def _negative(e: Expression) -> bool:
    """True if the numerical factor of a product (or the number) is negative (real or imaginary part)."""
    nums = [e] if e.get_type() == AtomType.Num else ([f for f in e if f.get_type() == AtomType.Num]
                                                    if e.get_type() == AtomType.Mul else [])
    return bool(nums) and nums[0].to_canonical_string().lstrip("(").startswith("-")


def _with_loop(e: Expression, part: str) -> Expression:
    """``e`` with the loop factor of ``part`` written first: ``1/(16 pi^2) 1/eps x`` (``lfe``), ``1/(16 pi^2) x``
    (``lf``); an overall minus sign is taken in front and a sum is parenthesised."""
    lfac = {"pole": "lfe", "evanescent": "lfe", "finite": "lf"}.get(part)
    if lfac is None:
        return e
    neg = _negative(e)
    inner = -e if neg else e
    if inner.get_type() == AtomType.Add:
        inner = mv.parse("par(x)").replace(mv.parse("x"), inner)
    out = mv.parse(f"ord({lfac}, x)").replace(mv.parse("x"), inner)
    out = -out if neg else out
    if "⋯" in mv.M(out)["html"]:
        # too large for one custom-printed fragment: Symbolica's own product (loop factor after the number)
        return mv.parse(lfac) * e
    return out


def _coef_site(c: Expression, part: str) -> Expression:
    return _with_loop(mv.coefficient(c), part)


def _link(spec_id: str, text: str) -> str:
    """Prose source of a link to a term; ``text`` is its plain name (symbols become inline MathML)."""
    return f'<a data-term="unbroken:{spec_id}">{nm.name_source(text)}</a>'


def _qlabel(n: str) -> str:
    """Inline name of a report monomial for prose: ``[[Q-symbol]]`` plus the bracket of Q_rho[t], Q_zeta[t]."""
    return f"[[{mv.M(q_symbol(n))['src']}]]" + q_suffix_html(n)


N_TEXT = {
    "N1": "(N1), which excludes the monomials that survive when all gauge fields, ghosts and sources vanish",
    "N2": "(N2), which excludes the transverse gauge two-point monomials",
    "N3": "(N3), which excludes the ghost-sector source monomials of the closed ghost and coupling directions",
    "N4": "(N4), which excludes every monomial with the abelian ghost",
}


def remarks(ctx: Context, spec: OpSpec, coeffs: dict, names: dict[str, str]) -> str | None:
    """One or two sentences specific to the term (rich text)."""
    out: list[tuple[int, str]] = []

    def add(priority: int, sentence: str) -> None:
        out.append((priority, sentence))

    n = spec.monomial
    F = ctx.fields
    if spec.partner and spec.partner in names:
        what = "Q^+ monomial" if n and n.startswith("y") else "monomial"
        add(0, f"Together with {_link(spec.partner, names[spec.partner])} it forms the {what} {_qlabel(n)} of "
                   "Table 8.")
    elif n and spec.family == "Q":
        add(0, f"Monomial {_qlabel(n)} of Table 8.")
    if spec.family == "Q" and "finite" in coeffs and not _has_spinor(spec.op) and ("tree" in coeffs or
                                                                                   "pole" in coeffs):
        add(2, "The tree-level and singular coefficients multiply the d-dimensional operator shown; the finite "
                   "coefficient multiplies its four-dimensional form, in which all fields, derivatives and metrics "
                   "are barred.")
    if n and n.startswith("L_"):
        add(2, "It is the only longitudinal counterterm; it cancels the finite longitudinal part that the closed "
                   "fermion loops give to the renormalised gauge self-energy.")
    if spec.family == "Q" and set(coeffs) == {"finite"}:
        add(2, "It has no tree-level and no singular coefficient: it is generated only by the restoring "
                   "condition, from the fermion loops with one insertion of [[Delta1]].")
    norm = normalisation_of(n)
    if norm and spec.family == "Q":
        add(3, f"Its finite coefficient is zero by the normalisation condition {N_TEXT[norm]}.")
    elif spec.family == "Q" and "finite" not in coeffs and n and not n.startswith("y"):
        add(3, "Its finite coefficient is zero in the unique solution of [[eqs(ord(bop, DeltaS1ST), -X1)]].")
    if any(F[x].kind == "antighost" for x in spec.cls):
        add(2, "The antighost enters only through "
                   "[[eqs(rhotilde(mu, a), rho(mu, a) + pd(up(mu), omegabar(a)))]] (ghost equation), so its "
                   "counterterm is that of the class with [[rho(mu, a)]] in place of [[omegabar(a)]].")
    if spec.family == "Q" and any(x == "gh_B" or F[x].factor == "B" and F[x].kind in ("ghost",) for x in spec.cls) \
            and "pole" not in coeffs:
        add(1, "The abelian ghost is free: every ghost and source counterterm of the hypercharge factor vanishes, "
                   "because [[eqs(dZ(omega, B), dZ(A, B))]] and [[eqs(dZ(omega, B), -2*dZc(gp))]].")
    if spec.id.startswith("gf-"):
        add(2, "The gauge fixing is not renormalised; all loop calculations use the Feynman gauge "
                   "[[eqs(xi, 1)]].")
    if spec.id.startswith("sterile-"):
        add(2, "The sterile partner has no gauge interaction and is not renormalised.")
    if spec.id.startswith("evkin-"):
        add(2, "This evanescent term is the only part of the classical action that is not BRST invariant; its "
                   "variation is the breaking [[Delta1]]. It is not rescaled by the bare quantities.")
    if n and n.startswith("V_"):
        add(2, "Only the physical chirality couples to the gauge field, through [[gammabar(mu)]].")
        if n[-1] in "GW" and "finite" in coeffs:
            add(1, "Its finite coefficient is fixed through "
                   f"{_eqa('sts-d1')}: [[eqs(rc(Qop(V, up(fr))), rc(Qop(R, up(fr))) + field(x, fr))]], so the classes "
                   "[[ord(omega, omega, Rbar, Psi)]] that fix the source coefficient [[rc(Qop(R, up(fr)))]] enter "
                   "as well as [[ord(omega, Psibar, Psi)]].")
    if n and n.startswith("K_"):
        add(2, "The kinetic counterterm contains only the physical chirality; the sterile and evanescent kinetic "
                   "terms are not renormalised.")
    if n == "m":
        add(2, "The two-point function receives [[dpar(m2)]] in addition to [[dZ(H)]].")
    if spec.evname in ("E_A^(G)", "E_A^(W)", "E_A^(B)"):
        r = spec.evname[-2]
        add(2, f"[[Eop(A, {r})]] is not [[bdop]]-closed: "
                   f"[[eqs(ord(bdop, Eop(A, {r})), Nhat(A, {r}) + Nhat(AA, {r}))]]. It comes from the closed fermion "
                   "loops.")
    if spec.evname == "E_H":
        add(2, "[[Eop(H)]] is not [[bdop]]-closed. In the two-point function it coincides with minus the "
                   "two-point part of [[Etil(H)]]; the coefficient of [[Etil(H)]] is fixed by its three- and "
                   "four-point "
                   "parts.")
    if spec.evname and spec.evname.startswith("Etilde_F") and len(spec.cls) > 2:
        add(1, f"[[eqs(Etil(F, {spec.evname[-2]}), -1/4*intd(minus(ord(bar(F(mu, nu, a)), bar(F(up(mu), up(nu), a))), "
               f"ord(F(mu, nu, a), F(up(mu), up(nu), a)))))]] (report Section 5.1); shown is its "
               f"{'three' if len(spec.cls) == 3 else 'four'}-point part. The relative sign of the unbarred and the "
               "barred piece sits in the order of the indices of the structure constants.")
    if spec.evname and spec.evname.startswith("Etilde"):
        add(2, "It is [[bdop]]-closed and has one coefficient for all its parts. It arises because the closed "
                   "fermion loops give barred invariants, which differ from the continued ones by this operator.")
    if spec.family == "breaking":
        f = [x for x in spec.cls if F[x].kind == "fermion"][0]
        gh = [x for x in spec.cls if F[x].kind == "ghost"][0]
        cid = class_id([gh, f, f + "bar"], 1)
        add(2, "Not part of the action: [[eqs(Delta1, ord(sop, S0))]] is the BRST variation of the evanescent "
                   "kinetic term. It enters the graphs of the ghost-number-1 Green's functions as an insertion, for "
                   f'example those of <a data-class="{cid}">the class with the legs of this vertex</a>.')
    if n and n.startswith("y") and spec.family == "Q":
        k = n[1]
        add(2, f"The imaginary combination [[{mv.M(q_symbol(f'y{k}-'))['src']}]] of the same monomials has zero "
                   "coefficients.")
    if not out:
        return None
    chosen = sorted(out, key=lambda ps: ps[0])[:3]
    return mv.prose("<p>" + " ".join(sentence for _, sentence in chosen) + "</p>")


def _eqa(label: str, text: str | None = None) -> str:
    """Rich-text link to a numbered equation of the report (``Eq. (N)`` with its page)."""
    from gdocs import theory

    return theory.eq_anchor(label, text)


def _groups(ctx: Context, spec: OpSpec) -> list[str]:
    F = ctx.fields
    out = set()
    for x in spec.cls:
        info = F[x]
        base = F[info.source_of] if info.kind == "source" else info
        if base.kind in ("gauge", "ghost", "antighost", "auxiliary"):
            out.add(base.factor)
    for r in FACTORS:
        if not _is_zero(spec.op.expr.derivative(ctx.model.param(COUPLING[r]))):
            out.add(r)
    return [GROUP[r] for r in FACTORS if r in out]


def _refs(ctx: Context, spec: OpSpec, coeffs: dict, cls_json: dict) -> list[dict]:
    refs = []
    if "tree" in coeffs:
        refs.append(ctx.ref("sec:brst") if spec.sector == "source" else ctx.ref("sec:action"))
    if spec.family == "Q":
        refs.append(ctx.table("Qbasis", "operator basis"))
    if spec.family == "breaking":
        refs += [ctx.ref("sec:delta"), ctx.ref("sec:feynman-rules")]
    if "pole" in coeffs or "evanescent" in coeffs:
        cid = pole_class(ctx, spec)
        if cid and cid in cls_json:
            refs += cls_json[cid]
        if "evanescent" in coeffs:
            refs += [ctx.ref("sec:evanescent"), ctx.ref("cr-ev")]
        else:
            refs += [ctx.table("c-Z", "renormalisation constants"), ctx.ref("cr-inv")]
    if "finite" in coeffs:
        refs += [ctx.table("sts-r", "finite coefficients"), ctx.ref("cr-fin")]
    out, seen = [], set()
    for r in refs:
        key = (r["label"], r["page"])
        if key not in seen:
            seen.add(key)
            out.append(r)
    return out[:6]


def term_json(ctx: Context, spec: OpSpec, coeffs: dict[str, Expression], *, names: dict[str, str],
              cls_refs: dict[str, list[dict]], bcache: dict) -> dict:
    """The Term object of DATA_FORMAT.md section 4."""
    model = ctx.model
    types = [p for p in PARTS if p in coeffs]
    op_site = mv.operator(spec.op)
    op4_site = mv.operator(spec.op4) if spec.op4 is not None and "finite" in coeffs else None
    differs = op4_site is not None and not _has_spinor(spec.op) and ("tree" in coeffs or "pole" in coeffs)
    # an operator with a finite coefficient only is the four-dimensional monomial P Q_alpha (report Section 5.1:
    # the restoring counterterm enters with every object barred)
    finite_only = op4_site is not None and not _has_spinor(spec.op) and set(types) == {"finite"}
    if finite_only:
        op_site = op4_site
    co = {}
    for p in types:
        c = coeffs[p]
        entry = {"math": _M(_coef_site(c, p), interactive=False)}
        if p == "pole" and spec.renormalised and "tree" in coeffs:
            dec, val = decomposition(ctx, spec, coeffs["tree"])
            if not _is_zero(val - c):
                raise RuntimeError(f"{spec.id}: decomposition does not reproduce the pole")
            entry["decomposition"] = _M(dec, interactive=False)
            entry["note"] = mv.prose("<p>First-order term of the monomial of [[S0]] written in bare fields and "
                                     "couplings: one [[1/2*dZ(X)]] per leg (with a minus sign for antighosts and "
                                     "sources) and [[dZc(c)]] or [[dpar(p)/p]] per power of a coupling. Each "
                                     "constant includes [[lfe]].</p>")
        if p == "evanescent":
            gen = {"E_A": ("-1/3*Ssum({r})", "S^(r)"), "Etilde_F": ("-2/3*Ssum({r})", "S^(r)"),
                   "E_H": ("-2/3*Ytil", None), "Etilde_H": ("-Ytil", None)}
            base, _, r = spec.evname.partition("^")
            src = gen[base][0].format(r=r.strip("()"))
            entry["decomposition"] = _M(mv.parse(f"lfe*({src})"), interactive=False)
            if base.startswith(("E_A", "Etilde_F")):
                rr = r.strip("()")
                val = {"G": "2*g_s^2", "W": "2*g^2", "B": "10/3*gp^2"}[rr]
                entry["note"] = mv.prose(f"<p>General form of the report (Section 5.1): [[Ssum({rr})]] is the sum of "
                                         f"[[Sinv(f, {rr})]] over the fermion multiplets, [[eqs(Ssum({rr}), {val})]] "
                                         "in the Standard Model.</p>")
            else:
                entry["note"] = mv.prose("<p>General form of the report (Section 5.1), with "
                                         "[[eqs(Ytil, 3*y_u^2 + 3*y_d^2 + y_e^2)]].</p>")
        if p == "finite" and spec.monomial:
            entry["symbol"] = _M(q_symbol(spec.monomial, "r"), interactive=False)
            note = "Coefficient [[" + mv.M(q_symbol(spec.monomial, "r"))["src"] + "]] of Table 56; it does not " \
                   "include [[lf]]."
            if differs:
                note += " It multiplies the four-dimensional operator [[" + mv.M(op4_site)["src"] + "]]."
            elif finite_only:
                note += (" It multiplies the four-dimensional operator shown: all fields, indices and derivatives are "
                         "barred (report Section 5.1).")
            entry["note"] = mv.prose("<p>" + note + "</p>")
        if p == "tree" and spec.id.startswith("gf-bb"):
            entry["note"] = mv.prose("<p>Feynman gauge: [[eqs(xi, 1)]] in all loop calculations.</p>")
        co[p] = entry
    # total: (sum of the coefficients, in the order tree, singular, finite) times the operator
    def ordered_sum(items: list[Expression]) -> Expression:
        acc = None
        for e in items:
            neg = _negative(e) or (e.get_type() == AtomType.Mul and any(
                f.get_type() == AtomType.Num and f.to_canonical_string() == "-1" for f in e))
            if acc is None:
                acc = e
            elif neg:
                acc = mv.parse("minus(x, y)").replace(mv.parse("x"), acc).replace(mv.parse("y"), -e)
            else:
                acc = mv.parse("plus(x, y)").replace(mv.parse("x"), acc).replace(mv.parse("y"), e)
        return acc

    def times(c: Expression, op: Expression, several: bool) -> Expression:
        if several:
            c = mv.parse("par(x)").replace(mv.parse("x"), c)
        if _negative(op) or op.get_type() == AtomType.Add:
            op = mv.parse("par(x)").replace(mv.parse("x"), op)
        if several:
            # a sum of several coefficients follows the operator, so that the operator stays in view
            return mv.parse("ord(y, x)").replace(mv.parse("x"), c).replace(mv.parse("y"), op)
        return mv.parse("ord(x, y)").replace(mv.parse("x"), c).replace(mv.parse("y"), op)

    def is_neg(e: Expression) -> bool:
        return _negative(e) or (e.get_type() == AtomType.Mul and any(
            f.get_type() == AtomType.Num and f.to_canonical_string() == "-1" for f in e))

    same = [p for p in types if not (p == "finite" and differs)]
    if len(same) == 1 and _is_zero(coeffs[same[0]] - 1):
        total = op_site
    else:
        total = times(ordered_sum([_coef_site(coeffs[p], p) for p in same]), op_site, len(same) > 1)
    if differs:
        cf = _coef_site(coeffs["finite"], "finite")
        neg = is_neg(cf)
        fin = times(-cf if neg else cf, op4_site, False)
        total = mv.parse(f"{'minus' if neg else 'plus'}(x, y)").replace(mv.parse("x"), total).replace(
            mv.parse("y"), fin)
    # vertex function
    legs = rule_legs(ctx, spec)
    parts = {}
    explicit = set()
    for p in types:
        info = {}
        opp = spec.op4 if p == "finite" and spec.op4 is not None else spec.op
        V = mv.vertex(opp, legs, info=info)
        explicit |= set(info.get("explicit") or ())
        parts[p] = _M(_with_loop(mv.coefficient(coeffs[p]) * V, p))
    leg_labels = mv.legs(legs, explicit=explicit or None)
    note = ("<p>The finite part is the vertex function of the four-dimensional operator: its metrics and momenta "
            "are barred.</p>" if differs else None)
    rule = {"title": "Vertex function",
            "legs": [{"field": f, "label": _M(lab, interactive=False, display=False)}
                     for f, lab in zip(legs, leg_labels)],
            "parts": parts,
            "note": mv.prose(note) if note else None}
    # classes
    classes = {"pole": [], "breaking": []}
    if "pole" in coeffs or "evanescent" in coeffs:
        cid = pole_class(ctx, spec)
        if cid is None:
            raise RuntimeError(f"{spec.id}: no ghost-number-0 class for {spec.cls}")
        classes["pole"] = [cid]
    if spec.family == "Q" and spec.monomial and not normalisation_of(spec.monomial) and \
            ("finite" in coeffs or "pole" in coeffs):
        classes["breaking"] = breaking_classes(ctx, spec.monomial, bcache)
        m = re.fullmatch(r"V_([QUDLE])([GW])", spec.monomial)
        if m and "finite" in coeffs:
            # r[Q_V^fr] = r[Q_R^fr] + x_fr (report eq:sts-d1): the classes that fix r[Q_R^fr] enter as well
            f, r = m.groups()
            for legs in ((f"gh_{r}", f"gh_{r}", f, f"src_{f}"), (f"gh_{r}", f"gh_{r}", f + "bar", f"src_{f}bar")):
                cid = class_id(list(legs), 1)
                if cid not in ctx.class_ids:
                    raise RuntimeError(f"{spec.id}: no class {cid}")
                if cid not in classes["breaking"]:
                    classes["breaking"].append(cid)
    F = ctx.fields
    values = {p: values_string(model, coeffs[p]) for p in types}
    return {
        "id": spec.id,
        "name": spec.name,
        "name_html": nm.name_html(spec.name),
        "operator": _M(op_site),
        "fields": list(spec.cls),
        "legs": len(spec.cls),
        "ghost_number": sum(F[x].ghost_number for x in spec.cls),
        "sector": spec.sector,
        "groups": _groups(ctx, spec),
        "evanescent": spec.evanescent,
        "types": types,
        "coefficients": co,
        "total": _M(total),
        "rules": [rule],
        "classes": classes,
        "theory": spec.theory,
        "remarks": remarks(ctx, spec, coeffs, names),
        "report": _refs(ctx, spec, coeffs, cls_refs),
        "data": {"monomial": spec.monomial, "report_name": _qname(spec.monomial)[2] if spec.monomial else
                 spec.evname, "values": values},
    }


# ---------------------------------------------------------------------------------------------- overview, fields
def _sym(base: str, sub: str | None = None, sup: str | None = None) -> str:
    """Source of a name symbol (registered in this process) for use in site-notation strings."""
    return mv.M(mv.symbol(base, sub=sub, sup=sup))["src"]


def _seq(*terms: str) -> str:
    """Site source of a sum printed in the written order (``plus``/``minus`` of mathview): ``terms`` are sources,
    a leading ``-`` makes a subtracted term."""
    acc = terms[0]
    for t in terms[1:]:
        t = t.strip()
        acc = f"minus({acc}, {t[1:].strip()})" if t.startswith("-") else f"plus({acc}, {t})"
    return acc


def overview(ctx: Context, counts: dict, consts) -> list[dict]:
    """The master formulae above the table (DATA_FORMAT.md section 4, ``overview``)."""
    Ok = _sym("O", sub="k")
    Lgf, Lext = _sym("ℒ", sub="gf"), _sym("ℒ", sub="ext")
    dSinv, dSev = _sym("δS", sub="inv", sup="(1)"), _sym("δS", sub="ev", sup="(1)")
    QTB, Qp = mv.M(q_symbol("T_B"))["src"], mv.M(mv.symbol("Q", sub=mv.parse("y_k"), sup="+"))["src"]
    cr = _sym("c", sub="r")
    wC = _sym("𝛚", sub="C")
    AW, AB = _sym("𝒜", sup="(W)"), _sym("𝒜", sup="(B)")
    AWmu, ABmu = _sym("𝒜", sub="μ", sup="(W)"), _sym("𝒜", sup="(B)μ")
    Ar = _sym("𝒜", sub="μ", sup="(r)")
    nu = _sym("ν", sub="fr")
    rinf = _sym("r∈f")
    out = []
    tot = counts["terms"]
    by_type = counts["by_type"]
    out.append({
        "id": "S-ren", "title": "The renormalised action",
        "math": _M(mv.parse("eqs(Sren, S0 + deltaS1 + DeltaS1ST)"), interactive=False),
        "text": mv.prose(
            f"<p>At one loop the renormalised action is the classical action, the singular counterterm (MS-bar poles "
            f"in [[ep]], with [[eqs(d, 4 - 2*ep)]]) and the finite counterterm that restores the Slavnov–Taylor "
            f"identity (report Section 5.1). The table splits it into {tot} operators with definite "
            f"field content: {by_type.get('tree', 0)} with a tree-level coefficient, {by_type.get('pole', 0)} with a "
            f"singular coefficient on a non-evanescent operator (continued to d dimensions as in report Section 5.1), {by_type.get('evanescent', 0)} parts of the eight "
            f"evanescent operators, {by_type.get('finite', 0)} with a finite coefficient, and the "
            f"{by_type.get('breaking', 0)} parts of the breaking [[Delta1]], which is not part of the action.</p>")})
    ferm = _seq("𝑖*ord(psibar(f), gammabar(mu), Dcov(mu), psi(f))", "𝑖*ord(chibar(f), gammabar(mu), pd(mu, chi(f)))",
                "𝑖*ord(psibar(f), gammahat(mu), pd(mu, chi(f)))", "𝑖*ord(chibar(f), gammahat(mu), pd(mu, psi(f)))")
    s0 = "eqs(S0, intd(" + _seq(
        "-1/4*ord(F(mu, nu, a), F(up(mu), up(nu), a))", f"sumover(f, par({ferm}))",
        "ord(dagger(par(ord(Dcov(mu), H))), Dcov(up(mu)), H)", "-m2*ord(Hdag, H)", "-lambda*par(ord(Hdag, H))^2",
        f"-sumover(k, y_k*par(plus({Ok}, dagger({Ok}))))", Lgf, Lext) + "))"
    out.append({
        "id": "S0", "title": "Classical action",
        "math": _M(mv.parse(s0), interactive=False),
        "text": mv.prose(
            "<p>The regularised classical action in [[eqs(d, 4 - 2*ep)]] dimensions (report Section 3.3): the "
            "d-dimensional gauge term; for each multiplet [[f]] (Q, U, D, L, E) the physical field "
            "[[eqs(psi(f), ord(Pc(f), Psi(f)))]] with the four-dimensional covariant kinetic term, its sterile partner "
            "[[eqs(chi(f), ord(Pcbar(f), Psi(f)))]] with a free four-dimensional kinetic term, and the evanescent "
            "kinetic term that couples the two; the Higgs doublet; the Yukawa terms with "
            "[[eqs(ord(nm_O_u), ord(Qbar, Htilde, PR, U))]].</p>".replace("ord(nm_O_u)", _sym("O", sub="u")))
        + mv.prose("<p>The evanescent kinetic term is the only part that is not BRST invariant (Section 3.5).</p>")})
    Dw = "par(plus(pd(mu, omega(a)), fg(a, b, c)*ord(A(mu, b), omega(c))))"
    gf = f"eqs(plus({Lgf}, {Lext}), " + _seq(
        "ord(b(a), pd(up(mu), A(mu, a)))", "1/2*xi*ord(b(a), b(a))", f"ord(pd(up(mu), omegabar(a)), {Dw})",
        f"ord(rho(mu, a), {Dw})", "-1/2*fg(a, b, c)*ord(zsrc(a), omega(b), omega(c))",
        "sumover(f, par(minus(ord(Rbar(f), sop, Psi(f)), ord(par(ord(sop, Psibar(f))), R(f)))))",
        "ord(sigmadag, sop, H)", "ord(par(ord(sop, Hdag)), sigma)") + ")"
    out.append({
        "id": "S0-brst", "title": "Gauge fixing, ghosts and BRST sources",
        "math": _M(mv.parse(gf), interactive=False),
        "text": mv.prose(
            "<p>Linear covariant gauge fixing with the auxiliary field [[b(a)]], the ghost terms and the sources of "
            "the "
            f"non-linear BRST variations (report {_eqa('Lgf')} and {_eqa('Lext')}), with "
            "[[eqs(ord(sop, Psi(f)), 𝑖*ord(omega(a), Tg(a, f), Pc(f), Psi(f)))]] and "
            "[[eqs(ord(sop, H), 𝑖*ord(omega(a), Tg(a, H), H))]]; the sterile partners do not transform. All loop "
            "calculations use the Feynman gauge [[eqs(xi, 1)]]. The gauge fixing is not renormalised.</p>")})
    d1rhs = ("sumover(f, intd(ord(omega(a), par(plus(ord(Psibar(f), Tg(a, f), gammahat(mu), Pcbar(f), "
             "pd(mu, Psi(f))), ord(pd(mu, Psibar(f)), gammahat(mu), Pc(f), Tg(a, f), Psi(f)))))))")
    out.append({
        "id": "Delta1", "title": "Tree-level breaking",
        "math": _M(mv.parse(f"eqs(Delta1, {d1rhs})"), interactive=False),
        "text": mv.prose(
            f"<p>[[eqs(Delta1, ord(sop, S0))]] (report {_eqa('Delta1')}) couples a ghost, a physical fermion and its sterile "
            "partner through an evanescent matrix; [[eqs(ord(Pphys, Delta1), 0)]]. It is not part of the action. The "
            "Green's functions with one insertion of [[Delta1]] give the finite breaking [[X1]].</p>")})
    inv = f"eqs({dSinv}, " + _seq(
        "sumover(r, par(" + _seq("1/2*ord(dZ(A, r), Oop(A, r))", f"ord(dZc({cr}), Oop(g, r))",
                                 "1/2*ord(dZ(omega, r), Oop(c, r))") + "))",
        f"ord(dZ(A, B), {QTB})", "sumover(f, 1/2*ord(dZ(f), Oop(psi, f)))", "1/2*ord(dZ(H), Oop(H))",
        f"sumover(k, ord(dpar(y_k), {Qp}))", "-ord(dpar(lambda), Qop(lambda))", "-ord(dpar(m2), Qop(m))") + ")"
    items = "".join(f"<li>[[eqs({src}, lfe*({mv.M(mv.coefficient(v))['src']}))]]</li>" for _, src, v in consts)
    out.append({
        "id": "dS-inv", "title": "Singular counterterm: invariant part",
        "math": _M(mv.parse(inv), interactive=False),
        "text": mv.prose(
            "<p>The sum over [[r]] runs over [[G]] and [[W]]. The closed directions [[Oop(A, r)]], [[Oop(g, r)]], "
            "[[Oop(c, r)]], [[Oop(psi, f)]], [[Oop(H)]] of Section 3.8.5 are continued to d dimensions (unsplit "
            "indices "
            "in the bosonic monomials, the barred vertex in the fermion monomials). Equivalently, "
            f"[[{dSinv}]] is the first-order term of [[S0]] written in bare fields and couplings; the gauge fixing, "
            "the sterile partners and the evanescent kinetic term are not renormalised. The 20 constants of the run "
            "(MS-bar, Feynman gauge; each includes [[lfe]], and [[dZc(g_s)]] stands for [[dpar(g_s)/g_s]]) agree "
            "with Table 57 of the report:</p>") + mv.prose("<ul>" + items + "</ul>")})
    ev = (f"eqs({dSev}, -lfe*par(plus(sumover(r, ord(Ssum(r), par(plus(2/3*Etil(F, r), 1/3*Eop(A, r))))), "
          "ord(Ytil, par(plus(Etil(H), 2/3*Eop(H)))))))")
    out.append({
        "id": "dS-ev", "title": "Singular counterterm: evanescent part",
        "math": _M(mv.parse(ev), interactive=False),
        "text": mv.prose(
            "<p>The sum runs over [[r]] = G, W, B. The operators are "
            "[[eqs(Eop(A, r), 1/2*intd(ord(bar(A(mu, a)), pdhat2(bar(A(up(mu), a))))))]], "
            "[[eqs(Eop(H), intd(ord(dagger(par(pdhat(mu, H))), pdhat(up(mu), H))))]], and [[Etil(F, r)]], [[Etil(H)]], "
            f"the barred minus the d-dimensional kinetic invariants (report {_eqa('evanescent')}). In the Standard Model "
            "[[eqs(Ssum(G), 2*g_s^2)]], [[eqs(Ssum(W), 2*g^2)]], [[eqs(Ssum(B), 10/3*gp^2)]] and "
            "[[eqs(Ytil, 3*y_u^2 + 3*y_d^2 + y_e^2)]]. [[Etil(F, r)]] and [[Etil(H)]] are [[bdop]]-closed; "
            "[[Eop(A, r)]] and [[Eop(H)]] are not.</p>")})
    tr = _seq("1/6*par(dot(pd, Am))^2", "-𝑖/3*ord(pd(mu, Am(nu)), comm(Am(up(mu)), Am(up(nu))))",
              "1/3*par(dot(Am, Am))^2", "-1/4*ord(Am(mu), Am(nu), Am(up(mu)), Am(up(nu)))")
    # Symbolica omits a custom-printed fragment above about 4 kB: the formula is split into three equations
    LF, LH = _sym("ℒ", sub="ψ", sup="ST"), _sym("ℒ", sub="H", sup="ST")
    ds = "eqs(DeltaS1ST, lf*int4(" + _seq(f"sumover(f, Tr({tr}))", LF, LH) + "))"
    dsf = f"eqs({LF}, " + _seq(
        f"-sumover(f, ord(psibar(f), gam(mu), par(plus(Cinv(f)*Am(mu), 1/2*sumover({rinf}, ord({nu}, {Ar})))), "
        "psi(f)))",
        f"-1/4*sumover(f, par(minus(𝑖*ord(Rbar(f), {wC}, psi(f)), 𝑖*ord(psibar(f), {wC}, R(f)))))") + ")"
    dsh = f"eqs({LH}, " + _seq(
        "1/3*Ytil*𝑖*par(minus(ord(dagger(par(pd(mu, H))), Am(up(mu)), H), ord(Hdag, Am(up(mu)), pd(mu, H))))",
        f"-Ytil*ord(Hdag, dot({AW}, {AW}), H)", f"-(y_u^2 + 2*y_d^2)*ord(Hdag, acomm({AWmu}, {ABmu}), H)",
        f"-(13/3*y_u^2 + 7/3*y_d^2 + 3*y_e^2)*ord(Hdag, dot({AB}, {AB}), H)",
        "-1/2*(y_u^2 + y_d^2)*g_s^2*ord(Hdag, H, G(mu, A), G(up(mu), A))") + ")"
    out.append({
        "id": "DS", "title": "Finite restoring counterterm",
        "math": _M(mv.parse(ds), interactive=False),
        "text": mv.prose(
            "<p>The unique solution of the restoring condition under the normalisation conditions (N1)–(N4) (report "
            f"{_eqa('sts-result')}), in matrix notation: [[Am(mu)]] and [[wm]] in the representation of the multiplet, "
            "all fields, indices, derivatives and Dirac matrices four-dimensional. The fermion, source and Higgs "
            f"terms [[{LF}]] and [[{LH}]] are written out below. In the basis [[Qop(alpha)]] of Table 8, "
            "[[eqs(DeltaS1ST, lf*sumover(alpha, rc(alpha)*Qop(alpha)))]] with "
            f"{counts['r']} non-zero [[rc(alpha)]]; they are the finite coefficients of the table. There is no ghost "
            "term other than the fermion-source term and no Levi-Civita tensor.</p>")})
    out.append({
        "id": "DS-fermion", "title": "Finite restoring counterterm: fermion and source terms",
        "math": _M(mv.parse(dsf), interactive=False),
        "text": mv.prose(
            f"<p>With [[eqs({wC}, sumover(r, Cad(r)*{_sym('𝛚', sup='(r)')}))]], [[Cinv(f)]] the Casimir of the "
            f"multiplet and [[{nu}]] the Yukawa numbers of Section 4.2.4. Only the physical chirality "
            "[[psi(f)]] enters.</p>")})
    out.append({
        "id": "DS-higgs", "title": "Finite restoring counterterm: Higgs terms",
        "math": _M(mv.parse(dsh), interactive=False),
        "text": mv.prose(
            "<p>In the Higgs terms [[Am(mu)]] is the gauge matrix on the doublet, with the parts "
            f"[[{AW}]] and [[{AB}]] of the two factors, and [[eqs(Ytil, 3*y_u^2 + 3*y_d^2 + y_e^2)]]. There is no "
            "Higgs-source, [[ord(par(ord(Hdag, H))^2)]] or [[m2]] term.</p>")})
    out.append({
        "id": "ST", "title": "The restored Slavnov–Taylor identity",
        "math": _M(mv.parse("eqs(ord(bop, DeltaS1ST), -X1)"), interactive=False),
        "text": mv.prose(
            "<p>[[X1]] is the finite breaking: the four-dimensional projection of the finite part of the one-loop "
            "Green's functions with one insertion of [[Delta1]], a local functional of ghost number 1 "
            f"(report {_eqa('m-X1')}). [[bop]] is the linearised Slavnov–Taylor operator of [[S0]] in four dimensions. With "
            "this condition the one-loop effective action of [[Sren]] satisfies the Slavnov–Taylor identity in four "
            f"dimensions: [[eqs(ord(Pphys, STop(Gamma)), 0)]] at one loop for [[to(ep, 0)]] (report {_eqa('m-restore')}). "
            "[[X1]] has no anomaly component, because "
            "the anomaly coefficients [[dsym(a, b, c)]]-type traces vanish for the fermion content of the Standard "
            "Model.</p>".replace("[[dsym(a, b, c)]]-type traces", "D<sup>abc</sup>"))})
    return out


FIELD_NAMES = {
    "G": "gluon field", "W": "SU(2) gauge field", "B": "hypercharge gauge field",
    "Q": "quark doublet", "U": "up-quark singlet", "D": "down-quark singlet", "L": "lepton doublet",
    "E": "electron singlet", "H": "Higgs doublet", "Hdag": "conjugate Higgs doublet",
}


def chip_symbol(ctx: Context, name: str) -> Expression:
    """Symbol of a field in the field lists, chips and graph captions: ghosts, antighosts, auxiliary fields and the
    sources of gauge fields and ghosts carry their gauge factor as a subscript (ω_G, ω̄_W, b_B, ρ_G, ζ_W), as in the
    graph drawings; in operators a superscript on these symbols is an adjoint index."""
    F = ctx.fields
    info = F[name]
    head = {"ghost": "omega", "antighost": "omegabar", "auxiliary": "b"}.get(info.kind)
    if head:
        return nm.field_chip(head, info.factor)
    if info.kind == "source":
        of = F[info.source_of]
        if of.kind == "gauge":
            return nm.field_chip("rho", of.factor)
        if of.kind == "ghost":
            return nm.field_chip("zsrc", of.factor)
    return mv.field(name)


def fields_json(ctx: Context) -> list[dict]:
    """Every field, sterile partner, ghost, antighost, auxiliary field and source of the model."""
    F = ctx.fields
    kinds = {"gauge": "gauge", "ghost": "ghost", "antighost": "antighost", "auxiliary": "auxiliary",
             "fermion": "fermion", "fermion_bar": "fermion", "scalar": "scalar", "scalar_conj": "scalar",
             "source": "source"}
    out = []
    for name, info in F.items():
        base = F[info.source_of] if info.kind == "source" else info
        r = base.factor if base.kind in ("gauge", "ghost", "antighost", "auxiliary") else None
        if info.kind == "gauge":
            fname = FIELD_NAMES[name]
        elif info.kind in ("ghost", "antighost", "auxiliary"):
            fname = f"{GROUP[r]} {'auxiliary field' if info.kind == 'auxiliary' else info.kind}"
        elif info.kind in ("fermion", "scalar", "scalar_conj"):
            fname = FIELD_NAMES[name]
        elif info.kind == "fermion_bar":
            fname = f"Dirac conjugate of the {FIELD_NAMES[info.conjugate].split(' (')[0]}"
        else:
            src_of = info.source_of
            sb = F[src_of]
            if sb.kind == "gauge":
                fname = f"BRST source ρ of the {FIELD_NAMES[src_of]}"
            elif sb.kind == "ghost":
                fname = f"BRST source ζ of the {GROUP[sb.factor]} ghost"
            elif sb.kind == "fermion":
                fname = f"BRST source R̄ of the {FIELD_NAMES[src_of]}"
            elif sb.kind == "fermion_bar":
                fname = f"BRST source R of the conjugate {FIELD_NAMES[sb.conjugate]}"
            elif sb.kind == "scalar":
                fname = "BRST source σ† paired with H"
            else:
                fname = "BRST source σ paired with H†"
        note = ""
        if info.kind in ("fermion", "fermion_bar"):
            phys = "left" if F[info.multiplet].chirality == "L" else "right"
            note = (f"Dirac field of the multiplet {info.multiplet}: the physical {phys}-handed field and its sterile "
                    "partner of the opposite chirality; terms with the sterile partner are listed with this field.")
        out.append({"id": name, "symbol": _M(chip_symbol(ctx, name), interactive=False, display=False), "name": fname,
                    "kind": kinds[info.kind], "group": GROUP[r] if r else None, "ghost_number": info.ghost_number,
                    "statistics": "fermi" if info.odd else "bose", "note": note})
    for f in FERMIONS:
        chir = "right-handed" if F[f].chirality == "L" else "left-handed"
        for bar in ("", "bar"):
            out.append({"id": f"chi{bar}_{f}", "symbol": _M(mv.parse(f"chi{bar}({f})"), interactive=False,
                                                              display=False),
                        "name": f"sterile partner of {f} ({chir})" + (", Dirac conjugate" if bar else ""),
                        "kind": "sterile", "group": None, "ghost_number": 0, "statistics": "fermi",
                        "note": f"The component {'P_R' if chir == 'right-handed' else 'P_L'} Ψ_{f} of the Dirac field "
                                f"{f}; a gauge singlet that never couples to a gauge field. It has no separate leg: "
                                f"its terms are listed under {f} and {f}bar."})
    return out


FACETS = {
    "sector": [("gauge", "Gauge"), ("fermion", "Fermion"), ("scalar", "Higgs"), ("yukawa", "Yukawa"),
               ("gauge-fixing", "Gauge fixing"), ("ghost", "Ghost"), ("source", "BRST source"),
               ("breaking", "Breaking")],
    "type": [("tree", "Tree"), ("pole", "Singular"), ("evanescent", "Evanescent"), ("finite", "Finite restoring"),
             ("breaking", "Breaking")],
    "group": [("SU(3)", "SU(3)"), ("SU(2)", "SU(2)"), ("U(1)", "U(1)")],
}


def intro(ctx: Context, counts: dict) -> str:
    """The phase text of :func:`gdocs.theory.phase_intro`, followed by what this build derived and checked."""
    from gdocs import theory

    return theory.phase_intro("unbroken") + mv.prose(
        "<p>The operators are the monomials of the report's basis (<a data-report-page=\"33\">Table 8</a>), the "
        "evanescent operators and the tree-level terms outside that basis. Their coefficients are the exact "
        "decomposition of the automation's order store on these operators and agree with the report: the "
        f"{counts['r']} non-zero [[rc(alpha)]], the {counts['constants']} renormalisation constants and the "
        f"{counts['evanescent_operators']} evanescent operators.</p>")


# ---------------------------------------------------------------------------------------------- classes.json
def _class_label(ctx: Context, fields: Sequence[str], ghost_number: int) -> Expression:
    """``Gamma(legs)``; ``Delta_1 . Gamma(legs)`` for the classes with one insertion of the breaking."""
    legs = mv.legs(tuple(sorted(fields, key=lambda x: _leg_rank(ctx.fields, x))))
    args = ", ".join(mv.M(leg[0])["src"] for leg in legs)
    if ghost_number == 0:
        return mv.parse(f"Gamma({args})")
    return mv.parse(f"dot(Delta1, Gamma({args}))")


def decorate_classes(site_dir: str | Path = SITE, ctx: Context | None = None, terms: list[dict] | None = None,
                     log=print) -> dict:
    """Fill ``label``, ``result`` and ``summary`` of ``data/classes.json`` (module docstring)."""
    ctx = ctx or Context(site_dir, log=log)
    path = Path(site_dir) / "data" / "classes.json"
    data = json.loads(path.read_text())
    by_pole, by_break = {}, {}
    for t in terms or []:
        for cid in (t.get("classes") or {}).get("pole") or []:
            by_pole.setdefault(cid, []).append(t)
        for cid in (t.get("classes") or {}).get("breaking") or []:
            by_break.setdefault(cid, []).append(t)
    n_res = 0
    d1ops = {s.cls: s for s in _breaking_ops(ctx)}
    for c in data["classes"]:
        fields = tuple(c["fields"])
        g = c["ghost_number"]
        label = _class_label(ctx, fields, g)
        c["label"] = _M(label, interactive=False)
        cls = tuple(sorted(fields, key=ctx.fields.order))
        result = {}
        pole_note = None
        if g == 0:
            pole = LocalFunctional.integrate(restrict(ctx.dS, cls).scale(-ctx.lf))
            result["pole"] = (_M(mv.operator(simplify_groups(pole), loop="pole"), threshold=6000, budget=20000)
                              if not pole.is_zero() else _M(Expression.num(0)))
        else:
            result["pole"], pole_note = _g1_pole(ctx, cls, d1ops)
        if g == 1:
            x = restrict(ctx.X1, cls)
            result["breaking"] = (_M(mv.operator(simplify_groups(LocalFunctional.integrate(x.scale(ctx.lf))),
                                                 loop="finite"), threshold=6000, budget=20000)
                                  if not x.is_zero() else _M(Expression.num(0)))
        c["result"] = result
        n_res += 1
        parts = []
        if g == 0:
            parts.append("Pole part of the one-loop 1PI function, written as a local operator: "
                         "[[eqs(Kop(Gamma), -deltaS1)]] restricted to these legs.")
        else:
            parts.append("Green's function with one insertion of [[Delta1]]. Pole part: "
                         "[[-ord(bdop, deltaS1)]] restricted to these legs (pole identity). Finite breaking: "
                         "the component of [[X1]] on these legs (the projected coefficient of [[eta/ep]]).")
            if pole_note:
                parts.append(pole_note)
        if c["id"] in by_pole:
            names = ", ".join(f'<a data-term="unbroken:{t["id"]}">{nm.name_source(t["name"])}</a>' for t in by_pole[c["id"]])
            parts.append(f"It fixes the singular coefficients of: {names}.")
        if c["id"] in by_break:
            names = ", ".join(f'<a data-term="unbroken:{t["id"]}">{nm.name_source(t["name"])}</a>' for t in by_break[c["id"]])
            parts.append(f"Its finite breaking enters the restoring condition that fixes the finite coefficients of: "
                         f"{names}.")
        c["summary"] = mv.prose("<p>" + " ".join(parts) + "</p>")
    data.pop("mock", None)
    bad = _omissions(data)
    if bad:
        raise RuntimeError(f"classes.json: expressions with omitted parts: {bad[:10]}")
    path.write_text(json.dumps(data, ensure_ascii=False, separators=(",", ":")))
    log(f"classes.json: label, result and summary of {n_res} classes")
    return {"classes": n_res}


def _g1_pole(ctx: Context, cls: tuple[str, ...], d1ops: dict) -> tuple[dict, str | None]:
    """Pole part of a ghost-number-1 class, ``-b_d deltaS^(1)`` on its legs, in a form that shows that it is
    evanescent (report Section 4.2): on the classes omega Psibar Psi the exact multiple of the breaking
    ``Delta_1^(fr)`` (fingerprint decomposition), elsewhere the functional integrated by parts so that no
    derivative acts on the ghost (``int (d omega) X = - int omega d X``), with the projector chains simplified."""
    raw = restrict(ctx.bd_dS, cls)
    if raw.is_zero():
        return _M(Expression.num(0)), None
    if cls in d1ops:
        try:
            k = decompose(ctx, raw, [d1ops[cls].op], cls)[0]
        except ValueError:
            k = None
        if k is not None:
            spec = d1ops[cls]
            F = LocalFunctional.integrate(spec.op.scale(-ctx.lf * k))
            return (_M(mv.operator(F, loop="pole"), threshold=6000, budget=20000),
                    f"The pole is {_coef_text(-k)} times the breaking [[Delta1]] on these legs, written out; it "
                    "is evanescent.")
    gh = next(x for x in cls if ctx.fields[x].kind == "ghost")
    F = tidy(simplify_groups(LocalFunctional.integrate(raw.integrate_by_parts(gh).scale(-ctx.lf))))
    note = "The pole is written after integration by parts, with no derivative on the ghost."
    M = _M(mv.operator(F, loop="pole"), threshold=6000, budget=20000)
    if not ctx.nonzero_classes(LocalFunctional.integrate(raw.integrate_by_parts(gh)), dimension=4):
        note += " Its four-dimensional part vanishes (checked numerically by the build): the pole is evanescent."
        if "Tr(" in M["text"]:
            note += (" The display does not apply the trace identities of three generators; its terms without "
                     "a hatted object cancel only once these identities are used.")
    return M, note


def _coef_text(c: Expression) -> str:
    """``[[lfe*c]]`` as prose source."""
    return f"[[{mv.M(_with_loop(mv.coefficient(c), 'pole'))['src']}]]"


# ---------------------------------------------------------------------------------------------- build
def _omissions(obj, path: str = "") -> list[str]:
    """Paths of strings (outside interactive ``view`` payloads) that contain Symbolica's omission mark."""
    out = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k != "view":
                out += _omissions(v, f"{path}.{k}")
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            out += _omissions(v, f"{path}[{i}]")
    elif isinstance(obj, str) and "⋯" in obj:
        out.append(path)
    return out


def build_unbroken(site_dir: str | Path = SITE, *, classes: bool = True, log=print) -> dict:
    """Write ``data/unbroken.json`` (and decorate ``data/classes.json``); returns counts."""
    t0 = time.time()
    ctx = Context(site_dir, log=log)
    specs = operator_specs(ctx)
    log(f"{len(specs)} operators; decomposing the store class by class")
    cl = catalogue(ctx, specs)
    checks = check_against_report(ctx, cl)
    log(f"checks against the report passed: {checks}")
    consts = constants(ctx)
    names = {s.id: s.name for s, _ in cl}
    cls_json = json.loads((Path(site_dir) / "data" / "classes.json").read_text())
    groups = {g["key"]: g for g in cls_json["groups"]}
    cls_refs = {c["id"]: groups[c["group"]]["report"] for c in cls_json["classes"]}
    bcache: dict = {}
    terms = []
    for s, c in cl:
        terms.append(term_json(ctx, s, c, names=names, cls_refs=cls_refs, bcache=bcache))
    by_type, by_sector = {}, {}
    for t in terms:
        for ty in t["types"]:
            by_type[ty] = by_type.get(ty, 0) + 1
        by_sector[t["sector"]] = by_sector.get(t["sector"], 0) + 1
    counts = {"terms": len(terms), "by_type": by_type, "by_sector": by_sector, **checks}
    doc = {
        "phase": "unbroken",
        "title": "Unbroken phase",
        "intro": intro(ctx, counts),
        "overview": overview(ctx, counts, consts),
        "fields": fields_json(ctx),
        "facets": {k: [{"id": i, "label": lab} for i, lab in v] for k, v in FACETS.items()},
        "terms": terms,
        "constants": {name: {"symbol": _M(mv.parse(src), interactive=False, display=False),
                             "value": values_string(ctx.model, v)} for name, src, v in consts},
    }
    bad = _omissions(doc)
    if bad:
        raise RuntimeError(f"expressions with omitted parts: {bad[:10]}")
    out = Path(site_dir) / "data" / "unbroken.json"
    text = json.dumps(doc, ensure_ascii=False, separators=(",", ":"))
    out.write_text(text)
    log(f"wrote {out.name}: {len(terms)} terms, {len(text) / 1e6:.2f} MB; types {by_type}; sectors {by_sector}")
    if classes:
        decorate_classes(site_dir, ctx, terms, log=log)
    log(f"done in {time.time() - t0:.0f} s")
    return counts


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m gdocs.catalog", description=__doc__.split("\n\n")[0])
    ap.add_argument("--site", default=str(SITE), help="site directory (default: gamma5_docs/site)")
    ap.add_argument("--no-classes", action="store_true", help="do not decorate data/classes.json")
    args = ap.parse_args(argv)
    build_unbroken(args.site, classes=not args.no_classes)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
