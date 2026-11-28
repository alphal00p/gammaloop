"""Independent review of the end-to-end run (gates V5 to V8): ``renormalise("sm_one_family")`` from the
specification, compared with the report without the driver's comparison code (``tests/driver_golden.py``) and
without the package's evaluators.

* Renormalisation constants and beta functions: the run's exact polynomials against ``final_counterterms.json``
  (``delta_S.Z_factors``, ``beta_functions``), compared as rational polynomials with sympy.
* Singular counterterm ``deltaS^(1)`` (invariant and evanescent part): the vertex functions of the run's ``C``
  (``deltaS^(1) = (1/(16 pi^2)) C / eps``) on every class against the Feynman rules of report 05a, tables
  cr-tab-inv and cr-tab-full with cr-ev-rules, evaluated in ``D = 6`` (two hatted components, so the evanescent
  coefficients of ``E_A^(r)``, ``Etilde_F^(r)``, ``E_H``, ``Etilde_H`` are tested through the hatted and barred
  structures of the rules). Classes not listed there must vanish.
* Finite counterterm ``DeltaS^(1)`` (the ``r_alpha``): its vertex functions on every class against 05a cr-fin-2 to
  cr-fin-higgs and the trace form cr-fin-4 of the quartic gauge rules; no other class.
* ``X_1``: ``V_{X_1} + V_{b DeltaS^(1)} = 0`` on every class (the restoring condition, re-checked with the independent
  evaluator), and ``X_1`` directly in the classes ``gh A`` (``N^((r))_3``) and ``gh psibar psi`` (``N^(f r)_F``)
  against ``final_counterterms.json`` ``X.engine_names``, with the ``N`` monomials typed here from 03b tab:Nbasis.

Where a rule is "z times the tree rule", the tree rule is the package's ``S_0`` (gate V1, reviewed independently in
``test_independent_v1_rules.py``); a "barred" tree rule is the tree vertex function evaluated on four-dimensional
polarisations and momenta. Everything is evaluated with ``tests/independent_support.py`` (Weyl representation, own
generators, 62-bit prime) at random points.

Slow: one end-to-end run (about 25 min cold on 100 processes; a few minutes with a warm per-graph cache:
``HEPKIT_RENORM_CACHE=<dir outside the repository>``; ``HEPKIT_RENORM_PROCESSES`` sets the workers). Run with
``pytest --runslow tests/test_independent_v5_v8_run.py``.
"""

from __future__ import annotations

import itertools
import json
import os
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest
import sympy
from independent_support import IU, Evaluator, P, Sampler, contract, fe
from symbolica import AtomType, Expression

from hepkit_renorm import tensors as tn
from hepkit_renorm.driver import RunOptions, renormalise
from hepkit_renorm.functional import Density, Leg, LocalFunctional, Occurrence

pytestmark = [pytest.mark.slow,
              pytest.mark.filterwarnings("ignore:This process .* is multi-threaded:DeprecationWarning")]

DATA = Path(__file__).parent / "data"
FACTOR_DIMS = {"G": 3, "W": 2}
NAMES = ("g_s", "g", "gp", "y_u", "y_d", "y_e", "lambda", "m2")

# ---------------------------------------------------------------------------------------------- report data (typed)
# 03_setup tab:fields; 05a / tab:invariants (SM values, units of the couplings squared)
HYPERCHARGE = {"Q": Fraction(1, 6), "U": Fraction(2, 3), "D": Fraction(-1, 3), "L": Fraction(-1, 2), "E": Fraction(-1)}
COUPLING = {"G": "g_s", "W": "g", "B": "gp"}


# ---------------------------------------------------------------------------------------------- the run
@pytest.fixture(scope="module")
def run(sm, tmp_path_factory):
    out = tmp_path_factory.mktemp("independent_run")
    cache = os.environ.get("HEPKIT_RENORM_CACHE") or str(out / "graph_cache")
    opts = RunOptions(output=out, cache=cache,
                      processes=int(os.environ.get("HEPKIT_RENORM_PROCESSES", str(os.cpu_count() or 1))))
    result = renormalise(sm, opts)
    assert result.report.passed
    return result, out


@pytest.fixture(scope="module")
def golden():
    with open(DATA / "final_counterterms.json") as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------------------------- polynomials
def _sym(text: str) -> sympy.Expr:
    t = text.replace("g'", "gp").replace("lambda", "lam").replace("^", "**")
    return sympy.sympify(t)


def _poly(ctx, poly) -> sympy.Expr:
    names = [n.replace("lambda", "lam") for n in ctx.space.names]
    syms = [sympy.Symbol(n) for n in names]
    total = sympy.Integer(0)
    for e, c in poly.items():
        term = sympy.Rational(c.numerator, c.denominator) if not isinstance(c, tuple) else \
            sympy.Rational(c[0].numerator, c[0].denominator) + sympy.I * sympy.Rational(c[1].numerator, c[1].denominator)
        for s, n in zip(syms, e):
            term *= s ** n
        total += term
    return total


def test_renormalisation_constants(run, golden):
    """All 20 constants equal Z_factors (report convention: delta Z_Phi; delta p = p delta Z_p for couplings)."""
    result, _ = run
    sing = result.orders[0].singular
    ctx = sing.context
    z = golden["delta_S"]["Z_factors"]["values_in_kappa_units"]
    expect = {}
    for r, c in COUPLING.items():
        expect[f"Z[{r}]"] = z["dZ_A"][r]
        expect[f"Z[gh_{r}]"] = z["dZ_omega"][r]
        expect[f"d[{c}]"] = f"({z['dZ_g'][r]})*{c}"
    for f, v in z["dZ_f"].items():
        expect[f"Z[{f}]"] = v
    expect["Z[H]"] = z["dZ_H"]
    for k, v in z["dy"].items():
        expect[f"d[y_{k}]"] = v
    expect["d[lambda]"] = z["dlambda"]
    expect["d[m2]"] = z["dm2"]
    got = sing.constants
    assert set(got) == set(expect), (set(got) ^ set(expect))
    bad = {k: (str(_poly(ctx, got[k])), v) for k, v in expect.items()
           if sympy.expand(_poly(ctx, got[k]) - _sym(v)) != 0}
    assert not bad, bad


def test_beta_functions(run, golden):
    """beta_p = 2 delta p (one loop) equals beta_functions.beta_times_16pi2 for all 8 parameters."""
    result, _ = run
    beta = result.report.results["1"]["beta"]
    want = golden["delta_S"]["beta_functions"]["beta_times_16pi2"]
    names = {"g_s": "g_s", "g": "g", "gp": "g'", "y_u": "y_u", "y_d": "y_d", "y_e": "y_e", "lambda": "lambda",
             "m2": "m^2"}
    bad = {k: (beta[k], want[v]) for k, v in names.items() if sympy.expand(_sym(beta[k]) - _sym(want[v])) != 0}
    assert not bad, bad


# ---------------------------------------------------------------------------------------------- evaluation helpers
class Point:
    """Random data of one evaluation: parameters, momenta and polarisations of the legs (independent evaluator)."""

    def __init__(self, model, D: int, seed: int) -> None:
        self.model = model
        self.smp = Sampler(seed)
        self.params = {n: self.smp.element() for n in NAMES}
        self.params["xi"] = 1
        self.ev = Evaluator(D, self.params, FACTOR_DIMS)
        self.D = D

    def legs(self, fields):
        out = []
        for k, f in enumerate(fields, start=1):
            labs = tuple(sl.space.label(70_000 + 10 * k + j) for j, sl in enumerate(self.model.fields[f].slots))
            out.append(Leg(f, labs, k))
        return out

    def data(self, legs):
        moms = self.smp.momenta(len(legs), self.D)
        pols = {}
        for leg in legs:
            for lab in leg.labels:
                pols[self.ev.label_key(lab)] = self.smp.vector(self.ev.dim(self.ev.space(lab)))
        return moms, pols

    def par(self, text: str) -> int:
        """A golden coefficient (Python syntax) at this point."""
        return self.ev._scalar(Expression.parse(text.replace("**", "^").replace("g'", "gp")), {})

    def value(self, functional: LocalFunctional, legs, moms, pols, cache: dict | None = None) -> int:
        key = (id(functional), tuple((l.field, l.momentum) for l in legs))
        if cache is not None and key in cache:
            expr = cache[key]
        else:
            expr = functional.vertex_function(legs)
            if cache is not None:
                cache[key] = expr
        return self.ev.value(expr, pols, moms)

    def barred(self, legs, moms, pols):
        """Four-dimensional (barred) Lorentz polarisations and momenta."""
        m = self.ev.bar_mask
        bm = {k: v * m % P for k, v in moms.items()}
        bp = dict(pols)
        for leg in legs:
            for lab in leg.labels:
                if self.ev.space(lab)[0] == "lor":
                    key = self.ev.label_key(lab)
                    bp[key] = pols[key] * m % P
        return bm, bp


def _lor(pt: Point, leg):
    return next(pt.ev.label_key(l) for l in leg.labels if pt.ev.space(l)[0] == "lor")


def _nonlor(pt: Point, leg):
    return [pt.ev.label_key(l) for l in leg.labels if pt.ev.space(l)[0] not in ("lor", "spin")]


def _spin(pt: Point, leg):
    return next(pt.ev.label_key(l) for l in leg.labels if pt.ev.space(l)[0] == "spin")


def _dot(pt: Point, a, b, mask=None) -> int:
    w = pt.ev.metric if mask is None else pt.ev.metric * mask % P
    return int(np.sum(np.asarray(a, dtype=object) * w * np.asarray(b, dtype=object)) % P)


def _plain(a, b, mask=None) -> int:
    a = np.asarray(a, dtype=object)
    b = np.asarray(b, dtype=object)
    return int(np.sum(a * b * (1 if mask is None else mask)) % P)


def _deltas(pt: Point, la, lb, pols) -> int:
    """The Kronecker deltas of the group slots of two legs of the same multiplet (or gauge field)."""
    out = 1
    for ka, kb in zip(_nonlor(pt, la), _nonlor(pt, lb)):
        out = out * _plain(pols[ka], pols[kb]) % P
    return out


# ---------------------------------------------------------------------------------------------- invariants
def _inv(pt: Point):
    v = {n: pt.params[n] for n in NAMES}
    gs2, g2, gp2 = v["g_s"] ** 2 % P, v["g"] ** 2 % P, v["gp"] ** 2 % P
    C_ad = {"G": 3 * gs2 % P, "W": 2 * g2 % P, "B": 0}
    S = {"G": 2 * gs2 % P, "W": 2 * g2 % P, "B": fe(Fraction(10, 3)) * gp2 % P}
    S_H = {"G": 0, "W": fe(Fraction(1, 2)) * g2 % P, "B": fe(Fraction(1, 2)) * gp2 % P}
    C_H = (fe(Fraction(3, 4)) * g2 + fe(Fraction(1, 4)) * gp2) % P
    Yt = (3 * v["y_u"] ** 2 + 3 * v["y_d"] ** 2 + v["y_e"] ** 2) % P
    return v, C_ad, S, S_H, C_H, Yt


def _f(x) -> int:
    return fe(Fraction(x))


# ---------------------------------------------------------------------------------------------- deltaS^(1)
def _expected_singular(pt: Point, model, tree: LocalFunctional, fields, golden, cache) -> int | None:
    """V of C (units 1/(16 pi^2 eps)) from 05a cr-tab-inv / cr-tab-full; None = no rule (must vanish)."""
    fi = [model.fields[f] for f in fields]
    kinds = sorted(x.kind for x in fi)
    v, C_ad, S, S_H, C_H, Yt = _inv(pt)
    z = golden["delta_S"]["Z_factors"]["values_in_kappa_units"]
    legs = pt.legs(fields)
    moms, pols = pt.data(legs)
    ev = pt.ev

    def tree_v(bar=False):
        if bar:
            m, p = pt.barred(legs, moms, pols)
            return pt.value(tree, legs, m, p, cache)
        return pt.value(tree, legs, moms, pols, cache)

    factors = {x.factor for x in fi if x.factor is not None}
    r = next(iter(factors)) if len(factors) == 1 else None
    if kinds == ["gauge", "gauge"] and len(factors) == 1:
        a, b = legs
        k = moms[1]
        e1, e2 = pols[_lor(pt, a)], pols[_lor(pt, b)]
        dl = _deltas(pt, a, b, pols)
        bar, hat = ev.bar_mask, ev.hat_mask
        T = (_dot(pt, k, k) * _dot(pt, e1, e2) - _plain(k, e1) * _plain(k, e2)) % P
        Tb = (_dot(pt, k, k, bar) * _dot(pt, e1, e2, bar) - _plain(k, e1, bar) * _plain(k, e2, bar)) % P
        Th = _dot(pt, k, k, hat) * _dot(pt, e1, e2, bar) % P
        c1 = (_f(Fraction(5, 3)) * C_ad[r] - _f(Fraction(1, 3)) * S_H[r]) % P
        c2, c3 = _f(Fraction(2, 3)) * S[r] % P, _f(Fraction(1, 3)) * S[r] % P
        # calibration of the reading: the tree two-point function of S_0 (gauge fixing through the auxiliary
        # field b) is the transverse -(p^2 g - p p) delta (report eq:propagators before eliminating b)
        assert tree_v() == (-T * dl) % P
        expect = (-(c1 * T - c2 * Tb - c3 * Th) * dl) % P
        return expect, legs, moms, pols
    if kinds == ["gauge"] * 3 and r is not None:
        val = ((_f(Fraction(2, 3)) * C_ad[r] - _f(Fraction(1, 3)) * S_H[r]) * tree_v()
               - _f(Fraction(2, 3)) * S[r] * tree_v(True)) % P
        return val, legs, moms, pols
    if kinds == ["gauge"] * 4 and r is not None:
        val = (-_f(Fraction(1, 3)) * (C_ad[r] + S_H[r]) * tree_v() - _f(Fraction(2, 3)) * S[r] * tree_v(True)) % P
        return val, legs, moms, pols
    if kinds == ["antighost", "ghost"] and r is not None:
        return _f(Fraction(1, 2)) * C_ad[r] * tree_v() % P, legs, moms, pols
    if kinds == ["antighost", "gauge", "ghost"] and r is not None:
        return -_f(Fraction(1, 2)) * C_ad[r] * tree_v() % P, legs, moms, pols
    srcs = [x for x in fi if x.kind == "source"]
    if len(srcs) == 1 and kinds.count("ghost") >= 1:
        src = srcs[0]
        target = model.fields[src.source_of]
        ghost_factor = next(x.factor for x in fi if x.kind == "ghost")
        if target.kind == "gauge" and kinds == ["ghost", "source"]:
            return _f(Fraction(1, 2)) * C_ad[ghost_factor] * tree_v() % P, legs, moms, pols
        if (target.kind == "gauge" and kinds == ["gauge", "ghost", "source"]) or \
                (target.kind == "ghost" and kinds == ["ghost", "ghost", "source"]) or \
                target.kind in ("fermion", "fermion_bar", "scalar", "scalar_conj"):
            return -_f(Fraction(1, 2)) * C_ad[ghost_factor] * tree_v() % P, legs, moms, pols
    if kinds == ["fermion", "fermion_bar"]:
        order = sorted(fields, key=lambda f: 0 if model.fields[f].kind == "fermion_bar" else 1)
        legs = pt.legs(order)
        moms, pols = pt.data(legs)
        lb, lf = legs
        mult = model.fields[lf.field].multiplet
        k = moms[lf.momentum]
        eb, ef = pols[_spin(pt, lb)], pols[_spin(pt, lf)]
        dl = _deltas(pt, lb, lf, pols)
        slash = sum(ev.metric[m] * k[m] * ev.gammas[m] for m in range(pt.D)) % P
        slashbar = sum(ev.metric[m] * k[m] * ev.gammas[m] for m in range(4)) % P
        one = np.eye(ev.S, dtype=object)
        proj = (one + (-1 if model.fields[lf.field].chirality == "L" else 1) * ev.g5) * fe(Fraction(1, 2)) % P
        typed_tree = int(eb.dot(slash % P).dot(ef) % P) * dl % P
        assert pt.value(tree, legs, moms, pols, cache) == typed_tree        # calibration: V_tree = pslash
        zf = pt.par(z["dZ_f"][mult])
        return zf * int(eb.dot(slashbar.dot(proj) % P).dot(ef) % P) % P * dl % P, legs, moms, pols
    if kinds == ["scalar", "scalar_conj"]:
        a, b = legs
        k = moms[1]
        dl = _deltas(pt, a, b, pols)
        p2, pb2, ph2 = (_dot(pt, k, k, m) for m in (None, ev.bar_mask, ev.hat_mask))
        m2, lam = v["m2"], v["lambda"]
        assert tree_v() == (p2 - m2) * dl % P                               # calibration: rule i (p^2 - m^2)
        val = (2 * C_H * p2 - Yt * pb2 - _f(Fraction(2, 3)) * Yt * ph2 - m2 * (6 * lam - C_H)) % P
        return val * dl % P, legs, moms, pols
    if kinds == ["gauge", "scalar", "scalar_conj"]:
        return ((2 * C_H - C_ad[r]) * tree_v() - Yt * tree_v(True)) % P, legs, moms, pols
    if kinds == ["gauge", "gauge", "scalar", "scalar_conj"]:
        rs = [x.factor for x in fi if x.kind == "gauge"]
        return ((2 * C_H - C_ad[rs[0]] - C_ad[rs[1]]) * tree_v() - Yt * tree_v(True)) % P, legs, moms, pols
    if kinds == ["scalar", "scalar", "scalar_conj", "scalar_conj"]:
        zz = (pt.par(z["dlambda"]) * pow(v["lambda"], -1, P) + 2 * pt.par(z["dZ_H"])) % P
        return zz * tree_v() % P, legs, moms, pols
    if kinds == ["fermion", "fermion_bar", "gauge"]:
        mult = next(x.multiplet for x in fi if x.kind == "fermion")
        rr = next(x.factor for x in fi if x.kind == "gauge")
        return (pt.par(z["dZ_f"][mult]) - C_ad[rr]) * tree_v() % P, legs, moms, pols
    if sorted(kinds) in (["fermion", "fermion_bar", "scalar"], ["fermion", "fermion_bar", "scalar_conj"]):
        mults = {x.multiplet for x in fi if x.kind in ("fermion", "fermion_bar")}
        gs2, g2, gp2 = (v[n] ** 2 % P for n in ("g_s", "g", "gp"))
        if "U" in mults:
            zz = -(_f(Fraction(16, 3)) * gs2 + _f(Fraction(3, 4)) * g2 + _f(Fraction(25, 36)) * gp2 + v["y_d"] ** 2)
        elif "D" in mults:
            zz = -(_f(Fraction(16, 3)) * gs2 + _f(Fraction(3, 4)) * g2 + _f(Fraction(1, 36)) * gp2 + v["y_u"] ** 2)
        else:
            zz = -(_f(Fraction(3, 4)) * g2 + _f(Fraction(9, 4)) * gp2)
        return zz % P * tree_v() % P, legs, moms, pols
    return None, legs, moms, pols


def _classes(F: LocalFunctional):
    return F.leg_classes()


def test_singular_counterterm_rules(run, golden, sm):
    """deltaS^(1) of the run equals 05a cr-tab-inv + cr-ev-rules (cr-tab-full) on every class, D = 6."""
    result, _ = run
    C = result.orders[0].singular.functional()
    tree = result.brst.action
    cache: dict = {}
    seen, bad, nonzero_unlisted = 0, [], []
    for seed in (11, 12):
        pt = Point(sm, 6, seed)
        for cls in _classes(C):
            expect, legs, moms, pols = _expected_singular(pt, sm, tree, cls, golden, cache)
            got = pt.value(C, legs, moms, pols, cache)
            if expect is None:
                if got:
                    nonzero_unlisted.append(cls)
                continue
            seen += 1
            if got != expect:
                bad.append(cls)
    assert not nonzero_unlisted, f"classes of deltaS without a rule in 05a: {sorted(set(nonzero_unlisted))}"
    assert not bad, f"{len(set(bad))} classes differ from 05a: {sorted(set(bad))}"
    assert seen >= 2 * 40


def test_singular_counterterm_rules_negative_control(run, golden, sm):
    """The comparison sees a change of one evanescent coefficient: Etilde_F^(G) scaled by 2 in the run's C."""
    result, _ = run
    sing = result.orders[0].singular
    name = next(e.name for e in sing.context.singular.evanescent if sing.coefficients.get(e.name)
                and e.name.endswith("[G G]") and "P(" in e.origin)
    coeffs = dict(sing.coefficients)
    coeffs[name] = {k: 2 * c for k, c in coeffs[name].items()}
    from hepkit_renorm.matching import _combine
    C2 = _combine(sing.context, sing.elements, coeffs)
    pt = Point(sm, 6, 13)
    expect, legs, moms, pols = _expected_singular(pt, sm, result.brst.action, ("G", "G"), golden, {})
    assert pt.value(C2, legs, moms, pols) != expect


# ---------------------------------------------------------------------------------------------- DeltaS^(1)
def _multiplet_generators(pt: Point, model, mult: str, factor: str, coeff_vec) -> np.ndarray:
    """``Tg = g_r c^a T^a`` of a gauge leg (adjoint polarisation ``coeff_vec``; scalar for U(1)) on the product
    space colour x weak of a fermion multiplet (zero matrix if the factor does not act)."""
    fm = model.fermions[mult]
    dims = [3 if fm.reps.get("G") not in (None, "1") else 1, 2 if fm.reps.get("W") not in (None, "1") else 1]
    gr = pt.params[COUPLING[factor]]
    if factor == "B":
        n = dims[0] * dims[1]
        return np.eye(n, dtype=object) * (gr * fe(HYPERCHARGE[mult]) * int(coeff_vec)) % P
    pos = 0 if factor == "G" else 1
    if dims[pos] == 1:
        return np.zeros((dims[0] * dims[1],) * 2, dtype=object)
    t = np.tensordot(np.asarray(coeff_vec, dtype=object), pt.ev.gens[factor], axes=(0, 0)) % P
    t = t * gr % P
    one = [np.eye(dims[0], dtype=object), np.eye(dims[1], dtype=object)]
    mats = [one[0], one[1]]
    mats[pos] = t
    return np.kron(mats[0], mats[1]) % P


def _quartic_trace(pt: Point, model, legs, pols) -> int:
    """05a cr-fin-4: sum over permutations and multiplets of Tr(Tg Tg Tg Tg)[1/3 gb gb - 1/4 gb gb]."""

    bar = pt.ev.bar_mask
    eps = [pols[_lor(pt, l)] for l in legs]
    cvec = []
    for l in legs:
        nl = _nonlor(pt, l)
        cvec.append(pols[nl[0]] if nl else 1)
    total = 0
    for mult in model.fermions:
        gens = [_multiplet_generators(pt, model, mult, model.fields[l.field].factor, c) for l, c in zip(legs, cvec)]
        for perm in itertools.permutations(range(4)):
            m = gens[perm[0]].dot(gens[perm[1]]) % P
            m = m.dot(gens[perm[2]]) % P
            m = m.dot(gens[perm[3]]) % P
            tr = int(np.trace(m) % P)
            if not tr:
                continue
            e = [eps[k] for k in perm]
            lor = (_f(Fraction(1, 3)) * _dot(pt, e[0], e[1], bar) * _dot(pt, e[2], e[3], bar)
                   - _f(Fraction(1, 4)) * _dot(pt, e[0], e[2], bar) * _dot(pt, e[1], e[3], bar)) % P
            total = (total + tr * lor) % P
    return total


def _expected_finite(pt: Point, model, tree: LocalFunctional, fields, cache):
    """V of sum r_alpha Q_alpha (units 1/(16 pi^2)) from 05a cr-fin-2 .. cr-fin-higgs; None: no rule."""
    fi = [model.fields[f] for f in fields]
    kinds = sorted(x.kind for x in fi)
    v, _, S, _, _, Yt = _inv(pt)
    legs = pt.legs(fields)
    moms, pols = pt.data(legs)
    ev = pt.ev
    factors = sorted({x.factor for x in fi if x.factor is not None})
    gs2, g2, gp2 = (v[n] ** 2 % P for n in ("g_s", "g", "gp"))
    yu2, yd2, ye2 = (v[n] ** 2 % P for n in ("y_u", "y_d", "y_e"))

    def tree_v(bar=False):
        if bar:
            m, p = pt.barred(legs, moms, pols)
            return pt.value(tree, legs, m, p, cache)
        return pt.value(tree, legs, moms, pols, cache)

    if kinds == ["gauge", "gauge"] and len(factors) == 1:
        r = factors[0]
        a, b = legs
        k = moms[1]
        e1, e2 = pols[_lor(pt, a)], pols[_lor(pt, b)]
        val = _f(Fraction(1, 3)) * S[r] * _plain(k, e1, ev.bar_mask) % P * _plain(k, e2, ev.bar_mask) % P
        return val * _deltas(pt, a, b, pols) % P, legs, moms, pols
    if kinds == ["gauge"] * 3 and len(factors) == 1 and factors[0] != "B":
        return -_f(Fraction(1, 3)) * S[factors[0]] * tree_v(True) % P, legs, moms, pols
    if kinds == ["gauge"] * 4:
        return _quartic_trace(pt, model, legs, pols), legs, moms, pols
    if kinds == ["fermion", "fermion_bar", "gauge"]:
        mult = next(x.multiplet for x in fi if x.kind == "fermion")
        r = next(x.factor for x in fi if x.kind == "gauge")
        cf = {("Q", "G"): (Fraction(4, 3), Fraction(3, 4), Fraction(1, 36), Fraction(1, 2), Fraction(1, 2), 0),
              ("Q", "W"): (Fraction(4, 3), Fraction(3, 4), Fraction(1, 36), 0, 0, 0),
              ("Q", "B"): (Fraction(4, 3), Fraction(3, 4), Fraction(1, 36), 2, -1, 0),
              ("U", "G"): (Fraction(4, 3), 0, Fraction(4, 9), 1, 0, 0),
              ("U", "B"): (Fraction(4, 3), 0, Fraction(4, 9), Fraction(1, 4), 0, 0),
              ("D", "G"): (Fraction(4, 3), 0, Fraction(1, 9), 0, 1, 0),
              ("D", "B"): (Fraction(4, 3), 0, Fraction(1, 9), 0, Fraction(-1, 2), 0),
              ("L", "W"): (0, Fraction(3, 4), Fraction(1, 4), 0, 0, 0),
              ("L", "B"): (0, Fraction(3, 4), Fraction(1, 4), 0, 0, 1),
              ("E", "B"): (0, 0, 1, 0, 0, Fraction(1, 2))}[(mult, r)]      # 05a cr-tab-ffv
        c = sum(_f(x) * y for x, y in zip(cf, (gs2, g2, gp2, yu2, yd2, ye2))) % P
        return -c * tree_v() % P, legs, moms, pols
    srcs = [x for x in fi if x.kind == "source"]
    if len(srcs) == 1 and kinds.count("ghost") == 1 and \
            model.fields[srcs[0].source_of].kind in ("fermion", "fermion_bar"):
        r = next(x.factor for x in fi if x.kind == "ghost")
        ra = {"G": -_f(Fraction(3, 4)) * gs2, "W": -_f(Fraction(1, 2)) * g2, "B": 0}[r] % P
        return ra * tree_v() % P, legs, moms, pols
    if kinds == ["gauge", "scalar", "scalar_conj"] and factors[0] in ("W", "B"):
        return -Yt * _f(Fraction(1, 3)) * tree_v(True) % P, legs, moms, pols
    if kinds == ["gauge", "gauge", "scalar", "scalar_conj"]:
        if factors == ["G"]:
            a, b = [l for l in legs if model.fields[l.field].kind == "gauge"]
            h = [l for l in legs if model.fields[l.field].kind != "gauge"]
            val = -(yu2 + yd2) * gs2 % P * _deltas(pt, a, b, pols) % P * _deltas(pt, h[0], h[1], pols) % P
            return val * _dot(pt, pols[_lor(pt, a)], pols[_lor(pt, b)], ev.bar_mask) % P, legs, moms, pols
        coeff = {("W",): Yt, ("B", "W"): yu2 + 2 * yd2,
                 ("B",): _f(Fraction(13, 3)) * yu2 + _f(Fraction(7, 3)) * yd2 + 3 * ye2}.get(tuple(factors))
        if coeff is None:
            return None, legs, moms, pols
        return -coeff * tree_v(True) % P, legs, moms, pols
    return None, legs, moms, pols


def test_finite_counterterm_rules(run, sm):
    """DeltaS^(1) of the run equals 05a cr-fin-2 to cr-fin-higgs (and the trace form cr-fin-4) on every class; it
    has no other vertex."""
    result, _ = run
    dS = result.orders[0].restoring.functional()
    tree = result.brst.action
    cache: dict = {}
    seen, bad, unlisted = 0, [], []
    for seed in (21, 22):
        pt = Point(sm, 6, seed)
        for cls in _classes(dS):
            expect, legs, moms, pols = _expected_finite(pt, sm, tree, cls, cache)
            got = pt.value(dS, legs, moms, pols, cache)
            if expect is None:
                if got:
                    unlisted.append(cls)
                continue
            seen += 1
            if got != expect:
                bad.append(cls)
    assert not unlisted, f"classes of DeltaS without a rule in 05a: {sorted(set(unlisted))}"
    assert not bad, f"{len(set(bad))} classes differ from 05a: {sorted(set(bad))}"
    assert seen >= 2 * 30
    assert len(result.orders[0].restoring.nonzero()) == 39


# ---------------------------------------------------------------------------------------------- X_1
@pytest.fixture(scope="module")
def X1(run, sm):
    _, out = run
    doc = json.loads((Path(out) / "sm_one_family_order_data.json").read_text())
    return LocalFunctional(sm.fields, Expression.parse(doc["orders"][0]["X"]), _canonical=True)


def test_breaking_restoring_condition(run, X1, sm):
    """V_{X_1} + V_{b DeltaS^(1)} = 0 on every ghost-number-1 class of either (D = 4, independent evaluator)."""
    result, _ = run
    dS = result.orders[0].restoring.functional()
    bdS = LocalFunctional.integrate(result.brst.b(dS))
    classes = sorted(set(X1.leg_classes()) | set(bdS.leg_classes()))
    pt = Point(sm, 4, 31)
    bad = []
    for cls in classes:
        legs = pt.legs(cls)
        moms, pols = pt.data(legs)
        _physical(pt, sm, legs, pols)                 # X_1 is defined on P N1
        total = (pt.value(X1, legs, moms, pols) + pt.value(bdS, legs, moms, pols)) % P
        if total:
            bad.append(cls)
    assert len(classes) > 30
    assert not bad, bad


def _physical(pt: Point, model, legs, pols) -> None:
    """Project the spinor polarisation of every fermion (and fermion-source) leg on its physical chirality
    (P_c for Psi, P_cbar for Psibar, that of the field for a source; gamma5 is symmetric here)."""
    for leg in legs:
        fi = model.fields[leg.field]
        base = model.fields[fi.source_of] if fi.kind == "source" else fi
        if base.kind not in ("fermion", "fermion_bar"):
            continue
        chir = base.chirality if base.kind == "fermion" else {"L": "R", "R": "L"}[base.chirality]
        proj = (np.eye(pt.ev.S, dtype=object) + (-1 if chir == "L" else 1) * pt.ev.g5) * fe(Fraction(1, 2)) % P
        key = _spin(pt, leg)
        pols[key] = proj.dot(pols[key]) % P


def _n3(model, r: str) -> LocalFunctional:
    """03b tab:Nbasis N^((r))_3 = int gh^a d^2 d.A^a."""
    mu, nu = tn.fresh(tn.LORENTZ), tn.fresh(tn.LORENTZ)
    adj = model.slot_labels(f"gh_{r}")
    return LocalFunctional.integrate(Density.monomial(model.fields, 1, Occurrence(f"gh_{r}", adj),
                                                      Occurrence(r, (mu,) + adj, (nu, nu, mu))))


def _nf(model, f: str, r: str) -> LocalFunctional:
    """03b tab:Nbasis N^(f r)_F = int psibar_f gamma^mu (d_mu omega^(r)) psi_f, omega = gh^a Tg^a, psi = P_c Psi."""
    fm = model.fermions[f]
    mu = tn.fresh(tn.LORENTZ)
    s1, s2, s3 = (tn.fresh(tn.SPINOR) for _ in range(3))
    bar_l, psi_l = list(model.slot_labels(f + "bar")), list(model.slot_labels(f))
    bar_l[0], psi_l[0] = s1, s3
    gl = model.slot_labels(f"gh_{r}")
    proj = tn.PL if fm.chirality == "L" else tn.PR
    coeff = model.coupling(r) * tn.gammabar(mu, s1, s2) * proj(s2, s3)
    slots = model.fields[f].slots
    for j, sl in enumerate(slots):
        if sl.space.kind == "spinor":
            continue
        if sl.space.factor == r:
            coeff = coeff * tn.T(gl[0], bar_l[j], psi_l[j])
        else:
            coeff = coeff * tn.delta(bar_l[j], psi_l[j])
    if r == "B":
        y = HYPERCHARGE[f]
        coeff = coeff * Expression.num(y.numerator) / Expression.num(y.denominator)
    return LocalFunctional.integrate(Density.monomial(
        model.fields, coeff, Occurrence(f + "bar", tuple(bar_l)), Occurrence(f"gh_{r}", gl, (mu,)),
        Occurrence(f, tuple(psi_l))))


def test_breaking_against_report_classes(X1, sm, golden):
    """X_1 in the classes gh A (N^((r))_3) and gh psibar psi (N^(f r)_F) equals X.engine_names (D = 4)."""
    eng = golden["X"]["engine_names"]
    pt = Point(sm, 4, 41)
    checks = []
    for r in ("G", "W", "B"):
        checks.append(((f"gh_{r}", r), _n3(sm, r), eng[f"cV3_{r}"]))
    for f in sm.fermions:
        for r in ("G", "W", "B"):
            key = f"F_{f}{r}"
            if key in eng:
                cls = tuple(sorted((f + "bar", f"gh_{r}", f), key=sm.fields.order))
                checks.append((cls, _nf(sm, f, r), eng[key]))
    assert len(checks) == 13
    bad = []
    for cls, mono, coeff in checks:
        legs = pt.legs(cls)
        moms, pols = pt.data(legs)
        _physical(pt, sm, legs, pols)
        got = pt.value(X1, legs, moms, pols)
        want = pt.par(coeff) * pt.value(mono, legs, moms, pols) % P
        if got != want or not want:
            bad.append((cls, coeff))
    assert not bad, bad


# ---------------------------------------------------------------------------------------------- the written file
def _parse_ufo(text: str) -> Expression:
    return Expression.parse(text.replace("bmhv::", "").replace("1𝑖", "iu").replace("𝑖", "iu"))


def _terms(expr: Expression) -> list[Expression]:
    e = expr.expand()
    return list(e) if e.get_type() == AtomType.Add else [e]


def _ufo_value(pt: Point, text: str, kind: str, legs: dict, moms: dict) -> int:
    """Value of a UFO Lorentz or colour structure of the written model (own reading, report conventions):
    ``Metric``, ``P(mu, k)`` (incoming momentum of leg ``k``), ``Gamma(mu, i, j) = (gamma^mu)_{ij}``,
    ``ProjM``/``ProjP``, ``bar``/``hat`` markers projecting an index; colour ``T(a, i, j) = (T^a)_{j i}`` (row = the
    antiquark index, as validated for the hand model in test_independent_v1_rules.py), ``f``, ``Identity``.
    ``legs[k]``: ``{"lor", "spin", "col", "colrep"}`` vectors of written leg ``k``."""
    ev = pt.ev
    total = 0
    for term in _terms(_parse_ufo(text)):
        factors = list(term) if term.get_type() == AtomType.Mul else [term]
        scalar, tensors, marks = 1, [], []
        for fa in factors:
            if fa.get_type() != AtomType.Fn:
                scalar = scalar * ev._scalar(fa, moms) % P
                continue
            name = fa.get_name().split("::")[-1]
            args = list(fa)

            def lab(a, space):
                head = a.get_name().split("::")[-1]
                nums = [int(x.to_canonical_string()) for x in a]
                key = f"{space}{nums[1]}" if head == "idx" else f"{space}d{nums[0]}"
                return key, (("lor", None) if space == "L" else ("spin", None))

            if kind == "lorentz":
                if name == "Gamma":
                    tensors.append((ev.gamma, [lab(args[0], "L"), lab(args[1], "S"), lab(args[2], "S")]))
                elif name in ("ProjM", "ProjP"):
                    sgn = -1 if name == "ProjM" else 1
                    proj = (np.eye(ev.S, dtype=object) + sgn * ev.g5) * fe(Fraction(1, 2)) % P
                    tensors.append((proj, [lab(args[0], "S"), lab(args[1], "S")]))
                elif name == "Metric":
                    tensors.append((np.diag(ev.metric), [lab(args[0], "L"), lab(args[1], "L")]))
                elif name == "P":
                    tensors.append((moms[legs[int(args[1].to_canonical_string())]["momentum"]],
                                    [lab(args[0], "L")]))
                elif name in ("bar", "hat"):
                    marks.append((name, lab(args[0], "L")[0]))
                else:
                    raise ValueError(name)
            else:
                def clab(a, rep):
                    k = int(a.to_canonical_string())
                    return (f"Cd{-k}", (rep, "G")) if k < 0 else (f"C{k}", (legs[k]["colrep"], "G"))

                if name == "T":
                    t = np.transpose(ev.gens["G"], (0, 2, 1))
                    tensors.append((t, [clab(args[0], "adj"), clab(args[1], "fund"), clab(args[2], "fund")]))
                elif name == "f":
                    tensors.append((ev.fd["G"][0], [clab(a, "adj") for a in args]))
                elif name == "Identity":
                    rep = legs[int(args[0].to_canonical_string())]["colrep"]
                    tensors.append((np.eye(3 if rep == "fund" else 8, dtype=object), [clab(a, rep) for a in args]))
                else:
                    raise ValueError(name)
        for which, key in marks:
            mask = ev.bar_mask if which == "bar" else ev.hat_mask
            for n, (a, labs) in enumerate(tensors):
                pos = [j for j, (k, _) in enumerate(labs) if k == key]
                if pos:
                    shape = [1] * a.ndim
                    shape[pos[0]] = -1
                    tensors[n] = (a * mask.reshape(shape) % P, labs)
                    break
        pol = {}
        for k, data in legs.items():
            for space, key in (("lor", f"L{k}"), ("spin", f"S{k}"), ("col", f"C{k}")):
                if data.get(space) is not None:
                    pol[key] = data[space]
        total = (total + scalar * contract(tensors, pol, ev)) % P
    return total


def _written_value(pt: Point, doc, rules, legs: dict, moms: dict, scale: Expression) -> int:
    lorentz = {x["name"]: x["structure"] for x in doc["lorentz_structures"]}
    couplings = {x["name"]: x["expression"] for x in doc["couplings"]}
    total = 0
    for rule in rules:
        for ci, col in enumerate(rule["color_structures"]):
            cval = _ufo_value(pt, col, "colour", legs, moms)
            for li, lname in enumerate(rule["lorentz_structures"]):
                gc = rule["couplings"][ci][li]
                if gc is None or not cval:
                    continue
                c = (_parse_ufo(couplings[gc]) * scale).expand()
                cpl = pt.ev._scalar(c, moms)
                total = (total + cpl * cval % P * _ufo_value(pt, lorentz[lname], "lorentz", legs, moms)) % P
    return total


def test_written_counterterm_vertices(run, golden, sm):
    """The counterterm vertices of the written renormalised model (UFO strings, read with an own reader), summed
    over the rules of one leg class, equal i (1/(16 pi^2)) (1/eps) times the 05a rules: every deltaS1 and DeltaS1
    vertex whose legs need no SU(2) components (gluons, B, U, D, E, their ghosts and sources)."""
    result, out = run
    doc = json.loads((Path(out) / "sm_one_family_renormalised.json").read_text())
    rmap = doc["hepkit_renorm"]["reverse_map"]
    tree = result.brst.action
    pi, eps = Expression.parse("pi"), Expression.parse("eps")
    groups: dict = {}
    for rule in doc["vertex_rules"]:
        rec = rmap["vertices"][rule["name"]]
        if rule["operator"] in ("deltaS1", "DeltaS1") and all(c == [] for c in rec["components"]):
            groups.setdefault((rule["operator"], tuple(rec["fields"])), []).append(rule)
    assert len(groups) >= 30
    pt = Point(sm, 6, 51)
    pt.ev.params["iu"] = IU
    cache: dict = {}
    bad = []
    for (op, fields), rules in sorted(groups.items()):
        if op == "deltaS1":
            expect, legs, moms, pols = _expected_singular(pt, sm, tree, fields, golden, cache)
            scale = Expression.num(16) * pi ** 2 * eps
        else:
            expect, legs, moms, pols = _expected_finite(pt, sm, tree, fields, cache)
            scale = Expression.num(16) * pi ** 2
        assert expect is not None, (op, fields)
        used: dict = {}
        wl = {}
        order = []                      # position of written leg k among the legs of the expectation
        for k, f in enumerate(fields, start=1):
            j = used.get(f, 0)
            used[f] = j + 1
            leg = [l for l in legs if l.field == f][j]
            order.append(legs.index(leg))
            data = {"momentum": leg.momentum, "colrep": None}
            for lab in leg.labels:
                kind = pt.ev.space(lab)
                key = pt.ev.label_key(lab)
                if kind[0] == "lor":
                    data["lor"] = pols[key]
                elif kind[0] == "spin":
                    data["spin"] = pols[key]
                elif kind[1] == "G":
                    data["col"], data["colrep"] = pols[key], kind[0]
            wl[k] = data
        got = _written_value(pt, doc, rules, wl, moms, scale)
        odd = [o for o, f in zip(order, fields) if sm.fields.odd(f)]
        sign = (-1) ** sum(1 for a, b in itertools.combinations(odd, 2) if a > b)    # written: package leg order
        if got != sign * IU * expect % P:
            bad.append((op, fields))
    assert not bad, bad
