"""Independent review, gate V1: the tree-level Feynman rules of the one-family Standard Model.

The expected values are taken from the report (03_setup, eq:propagators to eq:V-delta, typed below from the
report) and from the hand-written HepKIT model ``tests/data/sm_unbroken_bmhv.json`` (174 rules, verbatim copy of
``SCRATCHPAD/model/sm_unbroken_bmhv.json``). The package's rules are its derived tree-level action
``BRST(model).action`` and breaking ``BRST(model).breaking``, read through ``LocalFunctional.vertex_function``.
Both sides are evaluated with the independent evaluator of ``tests/independent_support.py`` (Weyl representation,
own SU(N) generators, a 62-bit prime field), never with the package's own evaluators.

Conventions: Feynman rule = ``i V``; the hand model lists its legs in the order in which the report writes the
rule, except for the ghost vertex (hand: ``omega, omegabar, A``; report: ``omegabar, A, omega``) and the source
vertices ``Psibar omega R`` and ``H^dag omega sigma`` (report: ghost first). The hand model encodes ghosts as
commuting, so its values are the report's values in the report's order. Colour: the UFO structure ``T(a, i, jb)``
(``i`` the quark leg, ``jb`` the antiquark leg) is the matrix element ``(T^a)_{jb i}``, row = barred leg, as in
``Psibar T Psi``.
"""

from __future__ import annotations

import itertools
import json
import re
from fractions import Fraction
from functools import cache
from pathlib import Path

import numpy as np
import pytest
from independent_support import IU, SQRT3, Evaluator, P, Sampler, contract, fe
from symbolica import AtomType, Expression, S

from hepkit_renorm import builtin_model
from hepkit_renorm.brst import BRST
from hepkit_renorm.functional import Leg
from hepkit_renorm.rules import propagators

DATA = Path(__file__).parent / "data"
D = 6
FACTOR_DIMS = {"G": 3, "W": 2}

# ---------------------------------------------------------------------------------------------- report data (typed)
Y = {"Q": Fraction(1, 6), "U": Fraction(2, 3), "D": Fraction(-1, 3), "L": Fraction(-1, 2), "E": Fraction(-1),
     "H": Fraction(1, 2)}                                     # tab:fields
CHIRALITY = {"Q": "L", "U": "R", "D": "R", "L": "L", "E": "R"}
NONABELIAN = {"Q": ("G", "W"), "U": ("G",), "D": ("G",), "L": ("W",), "E": (), "H": ("W",)}
COUPLING = {"G": "g_s", "W": "g", "B": "gp"}
COLOUR_ROW_IS_BARRED = True
"""Reading of the UFO colour structure T(a, i, jb); the negative control switches it off."""
HAND_PARAM = {"gs": "g_s", "gw": "g", "gp": "gp", "yu": "y_u", "yd": "y_d", "ye": "y_e", "lam": "lambda"}


def sym(name: str) -> Expression:
    return S(f"review::{name}")


@pytest.fixture(scope="module")
def model():
    return builtin_model("sm_one_family")


@pytest.fixture(scope="module")
def brst(model):
    return BRST(model)


@pytest.fixture(scope="module")
def params():
    smp = Sampler(101)
    vals = {n: smp.element() for n in ("g_s", "g", "gp", "y_u", "y_d", "y_e", "lambda", "m2", "xi")}
    vals.update({h: vals[p] for h, p in HAND_PARAM.items()})
    vals["iu"] = IU
    return vals


@pytest.fixture(scope="module")
def ev(params):
    return Evaluator(D, params, FACTOR_DIMS)


# ---------------------------------------------------------------------------------------------- legs and data
def make_leg(model, field: str, k: int) -> Leg:
    labs = tuple(sl.space.label(50_000 + 10 * k + j) for j, sl in enumerate(model.fields[field].slots))
    return Leg(field, labs, k)


def leg_pols(model, ev, leg: Leg, smp: Sampler, components: dict | None = None,
             shared: dict | None = None) -> dict:
    """Polarisations of a leg: random Lorentz / spinor / colour vectors (or ``shared`` ones, keyed by kind),
    a unit vector for a fixed SU(2) component."""
    out = {}
    for lab, sl in zip(leg.labels, model.fields[leg.field].slots):
        kind = sl.space.kind
        n = ev.dim(ev.space(lab))
        if sl.space.factor == "W" and components is not None and "W" in components:
            vec = np.zeros(n, dtype=object)
            vec[components["W"]] = 1
        elif shared is not None and kind in shared:
            vec = shared[kind]
        else:
            vec = smp.vector(n)
        out[ev.label_key(lab)] = vec
    return out


# ---------------------------------------------------------------------------------------------- hand model
_HAND_STEM = re.compile(r"^(Q|L|H|W|ghW|rhoW|zetaW|RQ|RL|sigma)(\d)$")


def hand_field(name: str) -> tuple[str, int | None]:
    """(package field, SU(2) component) of a particle name of the hand-written model."""
    bar = name.endswith("~")
    base = name[:-1] if bar else name
    m = _HAND_STEM.match(base)
    comp = int(m.group(2)) - 1 if m else None
    stem = m.group(1) if m else base
    if stem in ("Q", "U", "D", "L", "E"):
        return (stem + "bar" if bar else stem), comp
    if stem == "H":
        return ("Hdag" if bar else "H"), comp
    if stem in ("G", "W", "B"):
        return stem, comp
    if stem.startswith("gh"):
        return (("agh_" if bar else "gh_") + stem[2]), comp
    if stem.startswith("rho"):
        return "src_" + stem[3], comp
    if stem.startswith("zeta"):
        return "src_gh_" + stem[4], comp
    if stem.startswith("R"):                    # R_f = (Psibar_f)^*, Rbar_f = (Psi_f)^*
        return ("src_" + stem[1] if bar else "src_" + stem[1] + "bar"), comp
    if stem == "sigma":                         # sigma = (H^dag)^*, sigma^dag = H^*
        return ("src_H" if bar else "src_Hdag"), comp
    raise ValueError(name)


def hand_name(field: str, comp: int | None) -> str:
    c = "" if comp is None else str(comp + 1)
    if field in ("Q", "U", "D", "L", "E", "H"):
        return field + c
    if field.endswith("bar") and field[:-3] in ("Q", "U", "D", "L", "E"):
        return field[:-3] + c + "~"
    if field == "Hdag":
        return "H" + c + "~"
    if field in ("G", "W", "B"):
        return field + c
    if field.startswith("gh_"):
        return "gh" + field[3] + c
    if field.startswith("agh_"):
        return "gh" + field[4] + c + "~"
    if field.startswith("src_gh_"):
        return "zeta" + field[7] + c
    if field in ("src_H", "src_Hdag"):
        return "sigma" + c + ("~" if field == "src_H" else "")
    if field in ("src_G", "src_W", "src_B"):
        return "rho" + field[4] + c
    stem = field[4:]
    return ("R" + stem[:-3] + c) if stem.endswith("bar") else ("R" + stem + c + "~")


def report_order(model, fields: list[str]) -> list[int]:
    """Positions of the hand legs in the order in which the report writes the rule."""
    kinds = [model.fields[f].kind for f in fields]
    if kinds == ["ghost", "antighost", "gauge"]:
        return [1, 2, 0]
    if kinds in (["fermion_bar", "ghost", "source"], ["scalar_conj", "ghost", "source"]):
        return [1, 0, 2]
    return list(range(len(fields)))


@cache
def _parse(text: str) -> Expression:
    clean = text.replace("UFO::{}::", "").replace("bmhv::", "").replace("𝑖", "*iu")
    return Expression.parse(clean)


def _ufo_terms(expr: Expression) -> list[Expression]:
    e = expr.expand()
    return list(e) if e.get_type() == AtomType.Add else [e]


def ufo_value(text: str, kind: str, ev: Evaluator, legs: dict, momenta: dict) -> int:
    """Value of a UFO Lorentz or colour structure contracted with the leg polarisations.

    ``legs[k]`` gives, for hand leg ``k``: ``lor``, ``spin``, ``col`` vectors and ``colrep`` (``fund``/``adj``).
    """
    total = 0
    for term in _ufo_terms(_parse(text)):
        factors = list(term) if term.get_type() == AtomType.Mul else [term]
        scalar, tensors, marks = 1, [], []
        for fa in factors:
            if fa.get_type() != AtomType.Fn:
                scalar = scalar * ev._scalar(fa, momenta) % P
                continue
            name = fa.get_name().split("::")[-1]
            args = list(fa)
            if kind == "lorentz":
                def lab(a, space):
                    head = a.get_name().split("::")[-1]
                    nums = [int(x.to_canonical_string()) for x in a]
                    key = f"{space}{nums[1]}" if head == "idx" else f"{space}d{nums[0]}"
                    return key, (("lor", None) if space == "L" else ("spin", None))
                if name == "Gamma":
                    tensors.append((ev.gamma, [lab(args[0], "L"), lab(args[1], "S"), lab(args[2], "S")]))
                elif name in ("ProjM", "ProjP"):
                    sgn = -1 if name == "ProjM" else 1
                    proj = (np.eye(ev.S, dtype=object) + sgn * ev.g5) * fe(Fraction(1, 2)) % P
                    tensors.append((proj, [lab(args[0], "S"), lab(args[1], "S")]))
                elif name == "Metric":
                    tensors.append((np.diag(ev.metric), [lab(args[0], "L"), lab(args[1], "L")]))
                elif name == "P":
                    tensors.append((momenta[int(args[1].to_canonical_string())], [lab(args[0], "L")]))
                elif name in ("bar", "hat"):
                    marks.append((name, lab(args[0], "L")[0]))
                else:
                    raise ValueError(name)
            else:
                def clab(a, adj=None):
                    k = int(a.to_canonical_string())
                    if k < 0:
                        return f"Cd{-k}", ("adj", "G")
                    return f"C{k}", (legs[k]["colrep"], "G")
                if name == "T":                   # (T^a)_{jb i} indexed [a, i, jb]: row = antiquark leg
                    t = np.transpose(ev.gens["G"], (0, 2, 1)) if COLOUR_ROW_IS_BARRED else ev.gens["G"]
                    tensors.append((t, [clab(args[0]), clab(args[1]), clab(args[2])]))
                elif name == "f":
                    tensors.append((ev.fd["G"][0], [clab(a) for a in args]))
                elif name == "Identity":
                    rep = legs[int(args[0].to_canonical_string())]["colrep"]
                    tensors.append((np.eye(3 if rep == "fund" else 8, dtype=object), [clab(a) for a in args]))
                else:
                    raise ValueError(name)
        for which, key in marks:            # bar/hat: project one occurrence of the index
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


def hand_rule_value(rule: dict, lorentz: dict, couplings: dict, ev, legs, momenta) -> int:
    total = 0
    for ci, col in enumerate(rule["color_structures"]):
        cval = ufo_value(col, "colour", ev, legs, momenta)
        for li, lname in enumerate(rule["lorentz_structures"]):
            gc = rule["couplings"][ci][li]
            if gc is None:
                continue
            lval = ufo_value(lorentz[lname], "lorentz", ev, legs, momenta)
            cpl = ev._scalar(_parse(couplings[gc]), momenta)
            total = (total + cpl * cval % P * lval) % P
    return total


@pytest.fixture(scope="module")
def hand_model():
    with open(DATA / "sm_unbroken_bmhv.json") as fh:
        doc = json.load(fh)
    lorentz = {l["name"]: l["structure"] for l in doc["lorentz_structures"]}
    couplings = {c["name"]: c["expression"] for c in doc["couplings"]}
    return doc["vertex_rules"], lorentz, couplings


def _colrep(model, field):
    for sl in model.fields[field].slots:
        if sl.space.factor == "G":
            return "adj" if sl.space.kind == "adjoint" else "fund"
    return None


def compare_hand_rule(model, ev, functional, rule, lorentz, couplings, smp, cache_v) -> tuple[int, int]:
    """(i V of the package, hand value) for one hand rule on one random configuration."""
    names = rule["particles"]
    fc = [hand_field(n) for n in names]
    fields = [f for f, _ in fc]
    n = len(names)
    momenta = smp.momenta(n, D)
    legs_hand, pol = {}, {}
    pkg_legs = {}
    for k, (f, comp) in enumerate(fc, start=1):
        leg = make_leg(model, f, k)
        shared = {"lorentz": smp.vector(D), "spinor": smp.vector(ev.S)}
        rep = _colrep(model, f)
        if rep is not None:
            shared["fund" if rep == "fund" else "adjoint"] = smp.vector(3 if rep == "fund" else 8)
        pol.update(leg_pols(model, ev, leg, smp, {"W": comp} if comp is not None else None, shared))
        legs_hand[k] = {"lor": shared["lorentz"], "spin": shared["spinor"],
                        "col": shared.get("fund", shared.get("adjoint")), "colrep": rep}
        pkg_legs[k] = leg
    order = report_order(model, fields)
    legs = [pkg_legs[j + 1] for j in order]
    key = tuple((l.field, l.momentum) for l in legs)
    if key not in cache_v:
        cache_v[key] = functional.vertex_function(legs)
    v = ev.value(cache_v[key], pol, momenta)
    return IU * v % P, hand_rule_value(rule, lorentz, couplings, ev, legs_hand, momenta)


# ---------------------------------------------------------------------------------------------- tests: hand model
def test_hand_rules_equal_derived_rules(model, brst, ev, hand_model):
    """All 174 hand-written rules equal i V of the derived S_0 (or Delta_1 for the DELTA order), two points each."""
    rules, lorentz, couplings = hand_model
    assert len(rules) == 174
    smp = Sampler(7)
    cache_v: dict = {}
    bad = []
    for rule in rules:
        functional = brst.breaking if rule["name"].startswith("V_DELTA") else brst.action
        for _ in range(2):
            mine, hand = compare_hand_rule(model, ev, functional, rule, lorentz, couplings, smp, cache_v)
            if mine != hand or hand == 0:
                bad.append(rule["name"])
                break
    assert not bad, f"{len(bad)} hand rules disagree: {bad}"


def test_hand_rules_negative_control(model, brst, ev, hand_model):
    """The comparison is sensitive to Grassmann order: with the hand leg order (ghost not moved) exactly the 7 ghost
    vertices and the 23 Psibar-omega-R vertices flip sign (moving the ghost past the even H^dag changes nothing)."""
    rules, lorentz, couplings = hand_model
    smp = Sampler(8)
    global report_order
    saved = report_order
    try:
        report_order = lambda model, fields: list(range(len(fields)))
        flipped = []
        for rule in rules:
            functional = brst.breaking if rule["name"].startswith("V_DELTA") else brst.action
            mine, hand = compare_hand_rule(model, ev, functional, rule, lorentz, couplings, smp, {})
            if mine != hand:
                assert mine == (-hand) % P
                flipped.append(rule["name"])
    finally:
        report_order = saved
    kinds = [tuple(model.fields[hand_field(n)[0]].kind for n in r["particles"]) for r in rules
             if r["name"] in flipped]
    assert len(flipped) == 7 + 23           # ghost vertices and Psibar omega R; Hdag is even, so Hdag omega sigma
    assert all(k in (("ghost", "antighost", "gauge"), ("fermion_bar", "ghost", "source")) for k in kinds)


def test_hand_rules_colour_control(model, brst, ev, hand_model):
    """With the other reading of T(a, i, jb) exactly the rules with a quark and a gluon (or gluon ghost) fail."""
    global COLOUR_ROW_IS_BARRED
    rules, lorentz, couplings = hand_model
    smp = Sampler(9)
    COLOUR_ROW_IS_BARRED = False
    try:
        failing = set()
        for rule in rules:
            functional = brst.breaking if rule["name"].startswith("V_DELTA") else brst.action
            mine, hand = compare_hand_rule(model, ev, functional, rule, lorentz, couplings, smp, {})
            if mine != hand:
                failing.add(rule["name"])
    finally:
        COLOUR_ROW_IS_BARRED = True
    with_t = {r["name"] for r in rules if any("T(" in c for c in r["color_structures"])}
    assert failing == with_t and len(with_t) == 4 + 2 * 4 + 4   # ffv, R-omega-Psi pairs and Delta_1 with gluons


def test_evaluator_representation(ev):
    """The independent representation: Clifford algebra in D = 6, gamma5 = i g0 g1 g2 g3 anticommutes with the barred
    and commutes with the hatted gammas, tr(g^mu g^nu g^rho g^sigma g5) = -4 i eps (eq:tr5, tr 1 = 4),
    [T^a, T^b] = i f^{abc} T^c, Tr T^a T^b = delta/2, f^{123} = 1, f^{458} = sqrt(3)/2 for SU(3)."""
    g, g5, S_ = ev.gammas, ev.g5, ev.S
    one = np.eye(S_, dtype=object)
    for m, n in itertools.product(range(D), repeat=2):
        anti = (g[m].dot(g[n]) + g[n].dot(g[m])) % P
        assert (anti == (2 * ev.metric[m] * one % P if m == n else 0 * one)).all()
        c5 = (g5.dot(g[m]) + (1 if m < 4 else -1) * g[m].dot(g5)) % P
        assert not c5.any()
    assert (g5.dot(g5) % P == one).all()
    norm = fe(Fraction(4, S_))
    for idx in itertools.product(range(4), repeat=4):
        tr = int(np.trace(g[idx[0]].dot(g[idx[1]]).dot(g[idx[2]]).dot(g[idx[3]]).dot(g5)) % P) * norm % P
        assert tr == (-4 * IU * ev.eps[idx]) % P                # upper indices: eps^{mu nu rho sigma}, eps^{0123} = 1
    for F_, n in FACTOR_DIMS.items():
        t = ev.gens[F_]
        f, _ = ev.fd[F_]
        for a, b in itertools.product(range(n * n - 1), repeat=2):
            comm = (t[a].dot(t[b]) - t[b].dot(t[a])) % P
            assert (comm == IU * np.einsum("c,cij->ij", f[a, b], t) % P).all()
            assert int(np.trace(t[a].dot(t[b])) % P) == (fe(Fraction(1, 2)) if a == b else 0)
    f3 = ev.fd["G"][0]
    assert f3[0, 1, 2] == 1 and f3[3, 4, 7] == SQRT3 * fe(Fraction(1, 2)) % P


@pytest.mark.parametrize("seed", [11])
def test_hand_rules_are_complete(model, brst, ev, hand_model, seed):
    """Every SU(2) component of every derived vertex class (>= 3 legs, or a source) either has a hand rule with the
    same value or vanishes; every hand rule is reached."""
    rules, lorentz, couplings = hand_model
    by_key = {}
    for r in rules:
        key = tuple(sorted(r["particles"]))
        assert key not in by_key
        by_key[key] = r
    smp = Sampler(seed)
    reached, cache_v, nonzero_missing = set(), {}, []
    for functional, is_delta in ((brst.action, False), (brst.breaking, True)):
        for cls in functional.leg_classes():
            kinds = [model.fields[f].kind for f in cls]
            if "auxiliary" in kinds or (len(cls) < 3 and "source" not in kinds):
                continue
            comp_ranges = []
            for f in cls:
                w = [sl for sl in model.fields[f].slots if sl.space.factor == "W"]
                comp_ranges.append(range(3 if w[0].space.kind == "adjoint" else 2) if w else [None])
            seen = set()
            for comps in itertools.product(*comp_ranges):
                key = tuple(sorted(hand_name(f, c) for f, c in zip(cls, comps)))
                if key in seen:
                    continue
                seen.add(key)
                if key in by_key:
                    rule = by_key[key]
                    assert rule["name"].startswith("V_DELTA") == is_delta
                    mine, hand = compare_hand_rule(model, ev, functional, rule, lorentz, couplings, smp, cache_v)
                    assert mine == hand, rule["name"]
                    reached.add(key)
                    continue
                legs = [make_leg(model, f, k) for k, f in enumerate(cls, start=1)]
                momenta = smp.momenta(len(cls), D)
                pol = {}
                for leg, c in zip(legs, comps):
                    pol.update(leg_pols(model, ev, leg, smp, {"W": c} if c is not None else None))
                ck = tuple((l.field, l.momentum) for l in legs)
                if ck not in cache_v:
                    cache_v[ck] = functional.vertex_function(legs)
                if ev.value(cache_v[ck], pol, momenta) != 0:
                    nonzero_missing.append(key)
    assert not nonzero_missing, f"derived components without a hand rule: {nonzero_missing}"
    assert reached == set(by_key), f"hand rules never reached: {set(by_key) - reached}"


# ---------------------------------------------------------------------------------------------- tests: report rules
def _generator(model, fbar_legs, f_legs, factor: str, a_label, multiplet: str) -> Expression:
    """(Tg^a_f)_{row col} with row = the barred leg: T for the factor, deltas on the other group slots (typed from
    eq:Tg-explicit: for U(1) the hypercharge of tab:fields)."""
    expr = sym(COUPLING[factor]) if factor != "B" else sym(COUPLING["B"]) * Expression.num(Y[multiplet].numerator) \
        / Y[multiplet].denominator
    for (lb, sl), lf in zip(fbar_legs, f_legs):
        if sl.space.kind == "spinor":
            continue
        if sl.space.factor == factor:
            expr = expr * S("review::T")(a_label, lb, lf)
        else:
            expr = expr * S("review::delta")(lb, lf)
    return expr


def _rv(name):
    return S(f"review::{name}")


def _slots(model, leg):
    return list(zip(leg.labels, model.fields[leg.field].slots))


def _spin(model, leg):
    return next(l for l, sl in _slots(model, leg) if sl.space.kind == "spinor")


def _lor(model, leg):
    return next(l for l, sl in _slots(model, leg) if sl.space.kind == "lorentz")


def _adj(model, leg):
    return next((l for l, sl in _slots(model, leg) if sl.space.kind == "adjoint"), None)


def _check(ev, model, functional, legs, expected_v, smp, n_points=2):
    """V of the package equals the expected V (typed from the report) at random points."""
    v = functional.vertex_function(legs)
    for _ in range(n_points):
        momenta = smp.momenta(len(legs), D)
        pol = {}
        for leg in legs:
            pol.update(leg_pols(model, ev, leg, smp))
        assert ev.value(v, pol, momenta) == ev.value(expected_v, pol, momenta), [l.field for l in legs]
        assert ev.value(expected_v, pol, momenta) != 0


FERMION_FACTORS = [(f, r) for f in ("Q", "U", "D", "L", "E") for r in ("G", "W", "B")
                   if r == "B" or r in NONABELIAN[f]]


def test_fermion_rules_against_report(model, brst, ev):
    """eq:V-ffv, eq:V-src-matter (Rbar omega Psi, omega Psibar R) and eq:V-delta for all 10 (f, r) pairs."""
    assert len(FERMION_FACTORS) == 10
    smp = Sampler(21)
    s0, d1 = brst.action, brst.breaking
    for f, r in FERMION_FACTORS:
        pc, pcb = (_rv("PL"), _rv("PR")) if CHIRALITY[f] == "L" else (_rv("PR"), _rv("PL"))
        x = _rv("lor")(99_001)
        s = _rv("spin")(99_002)
        # Psibar(1) Psi(2) A(3): V = gbar^mu P_c Tg^a
        lb, lf, la = make_leg(model, f + "bar", 1), make_leg(model, f, 2), make_leg(model, r, 3)
        mu, a = _lor(model, la), _adj(model, la)
        tg = _generator(model, _slots(model, lb), lf.labels, r, a, f)
        exp = _rv("gammabar")(mu, _spin(model, lb), s) * pc(s, _spin(model, lf)) * tg
        _check(ev, model, s0, [lb, lf, la], exp, smp)
        # Rbar(1) omega(2) Psi(3): i V = -Tg P_c  ->  V = i Tg P_c
        lr, lg, lf3 = make_leg(model, "src_" + f, 1), make_leg(model, "gh_" + r, 2), make_leg(model, f, 3)
        tg = _generator(model, _slots(model, lr), lf3.labels, r, _adj(model, lg), f)
        _check(ev, model, s0, [lr, lg, lf3], Expression.I * tg * pc(_spin(model, lr), _spin(model, lf3)), smp)
        # omega(1) Psibar(2) R(3): i V = -P_cbar Tg  ->  V = i P_cbar Tg
        lg1, lb2, lr3 = make_leg(model, "gh_" + r, 1), make_leg(model, f + "bar", 2), make_leg(model, "src_" + f + "bar", 3)
        tg = _generator(model, _slots(model, lb2), lr3.labels, r, _adj(model, lg1), f)
        _check(ev, model, s0, [lg1, lb2, lr3], Expression.I * tg * pcb(_spin(model, lb2), _spin(model, lr3)), smp)
        # eq:V-delta: omega(q) Psibar(p' outgoing) Psi(p); V = -i Tg (phat-slash P_cbar - p'hat-slash P_c)
        lg1, lb2, lf3 = make_leg(model, "gh_" + r, 1), make_leg(model, f + "bar", 2), make_leg(model, f, 3)
        tg = _generator(model, _slots(model, lb2), lf3.labels, r, _adj(model, lg1), f)
        p_in, p_out = _rv("p")(3, x), -_rv("p")(2, x)          # p' = - (incoming momentum of leg 2)
        sb, sf = _spin(model, lb2), _spin(model, lf3)
        exp = -Expression.I * tg * (_rv("gammahat")(x, sb, s) * p_in * pcb(s, sf)
                                    - _rv("gammahat")(x, sb, s) * p_out * pc(s, sf))
        _check(ev, model, d1, [lg1, lb2, lf3], exp, smp)


def test_boson_rules_against_report(model, brst, ev):
    """eq:V-vvv, eq:V-vvvv, eq:V-ghost for G and W; eq:V-ssv, eq:V-ssvv for all factor pairs; eq:V-ssss."""
    smp = Sampler(22)
    s0 = brst.action
    g, f, p, delta = _rv("g"), _rv("f"), _rv("p"), _rv("delta")
    for r in ("G", "W"):
        c = sym(COUPLING[r])
        legs = [make_leg(model, r, k) for k in (1, 2, 3)]
        (m1, a1), (m2, a2), (m3, a3) = [(_lor(model, l), _adj(model, l)) for l in legs]
        exp = c * f(a1, a2, a3) * (g(m1, m2) * (p(1, m3) - p(2, m3)) + g(m2, m3) * (p(2, m1) - p(3, m1))
                                   + g(m3, m1) * (p(3, m2) - p(1, m2)))
        _check(ev, model, s0, legs, -Expression.I * exp, smp)          # the rule is fg[...]: V = -i fg[...]
        legs = [make_leg(model, r, k) for k in (1, 2, 3, 4)]
        (m1, a), (m2, b), (m3, cc), (m4, e) = [(_lor(model, l), _adj(model, l)) for l in legs]
        x = _rv("adj")(S(f"hepkit_renorm_factor::{r}"), 99_100)
        def ff(i, j, k, l, c=c, x=x):
            return c * c * f(i, j, x) * f(k, l, x)
        exp = -Expression.I * (ff(a, b, cc, e) * (g(m1, m3) * g(m2, m4) - g(m1, m4) * g(m2, m3))
                               + ff(a, cc, b, e) * (g(m1, m2) * g(m3, m4) - g(m1, m4) * g(m2, m3))
                               + ff(a, e, b, cc) * (g(m1, m2) * g(m3, m4) - g(m1, m3) * g(m2, m4)))
        _check(ev, model, s0, legs, -Expression.I * exp, smp)
        # omegabar^a(p' outgoing, leg 1) A^b_mu(2) omega^c(3): - fg^{abc} p'^mu
        legs = [make_leg(model, "agh_" + r, 1), make_leg(model, r, 2), make_leg(model, "gh_" + r, 3)]
        exp = -c * f(_adj(model, legs[0]), _adj(model, legs[1]), _adj(model, legs[2])) * (-p(1, _lor(model, legs[1])))
        _check(ev, model, s0, legs, -Expression.I * exp, smp)
    # Higgs: Hdag_i(p' outgoing, leg 1) H_j(p, leg 2) A^a_mu: i (p + p')^mu Tg_H
    for r in ("W", "B"):
        legs = [make_leg(model, "Hdag", 1), make_leg(model, "H", 2), make_leg(model, r, 3)]
        tg = _generator(model, _slots(model, legs[0]), legs[1].labels, r, _adj(model, legs[2]), "H")
        mu = _lor(model, legs[2])
        _check(ev, model, s0, legs, (p(2, mu) - p(1, mu)) * tg, smp)
    for r1, r2 in (("W", "W"), ("W", "B"), ("B", "B")):
        legs = [make_leg(model, "Hdag", 1), make_leg(model, "H", 2), make_leg(model, r1, 3), make_leg(model, r2, 4)]
        i, j = legs[0].labels[0], legs[1].labels[0]
        k = _rv("fund")(S("hepkit_renorm_factor::W"), 99_200)
        def tgh(r, a, row, col):
            if r == "W":
                return sym("g") * _rv("T")(a, row, col)
            return sym("gp") * delta(row, col) / 2
        a, b = _adj(model, legs[2]), _adj(model, legs[3])
        anti = tgh(r1, a, i, k) * tgh(r2, b, k, j) + tgh(r2, b, i, k) * tgh(r1, a, k, j)
        _check(ev, model, s0, legs, g(_lor(model, legs[2]), _lor(model, legs[3])) * anti, smp)
    legs = [make_leg(model, n, k) for k, n in enumerate(("Hdag", "H", "Hdag", "H"), start=1)]
    i, j, k, l = (leg.labels[0] for leg in legs)
    _check(ev, model, s0, legs, -2 * sym("lambda") * (delta(i, j) * delta(k, l) + delta(i, l) * delta(k, j)), smp)


def test_yukawa_and_gauge_source_rules_against_report(model, brst, ev):
    """The six rules of eq:V-yuk and eq:V-src-gauge (rho omega, rho A omega, zeta omega omega) for G, W, B."""
    smp = Sampler(23)
    s0 = brst.action
    delta, epsf = _rv("delta"), _rv("epsf")
    pl, pr = _rv("PL"), _rv("PR")
    cases = [("Qbar", "U", "Hdag", "y_u", "eps", pr), ("Ubar", "Q", "H", "y_u", "eps", pl),
             ("Qbar", "D", "H", "y_d", "delta", pr), ("Dbar", "Q", "Hdag", "y_d", "delta", pl),
             ("Lbar", "E", "H", "y_e", "delta", pr), ("Ebar", "L", "Hdag", "y_e", "delta", pl)]
    for fb, fq, h, y, w, proj in cases:
        legs = [make_leg(model, fb, 1), make_leg(model, fq, 2), make_leg(model, h, 3)]
        sl = {l.field: _slots(model, l) for l in legs}
        doublet = next(lab for l in legs[:2] for lab, s in _slots(model, l) if s.space.factor == "W")
        hlab = legs[2].labels[0]
        wexpr = epsf(doublet, hlab) if w == "eps" else delta(doublet, hlab)
        colour = Expression.num(1)
        cols = [lab for l in legs[:2] for lab, s in _slots(model, l) if s.space.factor == "G"]
        if cols:
            colour = delta(cols[0], cols[1])
        exp = -sym(y) * wexpr * colour * proj(_spin(model, legs[0]), _spin(model, legs[1]))
        _check(ev, model, s0, legs, exp, smp)               # i V = -i y ...  ->  V = -y ...
        del sl
    g, f, p = _rv("g"), _rv("f"), _rv("p")
    for r in ("G", "W", "B"):
        legs = [make_leg(model, "src_" + r, 1), make_leg(model, "gh_" + r, 2)]
        mu = _lor(model, legs[0])
        a, b = _adj(model, legs[0]), _adj(model, legs[1])
        dab = delta(a, b) if a is not None else Expression.num(1)
        _check(ev, model, s0, legs, -Expression.I * dab * p(2, mu), smp)    # i V = delta p_mu (p of the ghost)
        if r == "B":
            continue
        c = sym(COUPLING[r])
        legs = [make_leg(model, "src_" + r, 1), make_leg(model, r, 2), make_leg(model, "gh_" + r, 3)]
        exp = c * f(_adj(model, legs[0]), _adj(model, legs[1]), _adj(model, legs[2])) * \
            g(_lor(model, legs[0]), _lor(model, legs[1]))
        _check(ev, model, s0, legs, exp, smp)                               # i V = i fg g
        legs = [make_leg(model, "src_gh_" + r, 1), make_leg(model, "gh_" + r, 2), make_leg(model, "gh_" + r, 3)]
        exp = -c * f(_adj(model, legs[0]), _adj(model, legs[1]), _adj(model, legs[2]))
        _check(ev, model, s0, legs, exp, smp)                               # i V = -i fg
    for src, order in (("src_H", ("src_H", "gh_W", "H")), ("src_Hdag", ("gh_W", "Hdag", "src_Hdag"))):
        legs = [make_leg(model, n, k) for k, n in enumerate(order, start=1)]
        row = legs[0] if src == "src_H" else legs[1]
        col = legs[2]
        gh = legs[1] if src == "src_H" else legs[0]
        tg = sym("g") * _rv("T")(_adj(model, gh), row.labels[0], col.labels[0])
        sign = 1 if src == "src_H" else -1                  # i V = -T_H and +T_H
        _check(ev, model, s0, legs, sign * Expression.I * tg, smp)


def test_two_point_functions_and_propagators(model, brst, ev):
    """Tree two-point vertex functions (02_method sec:m-ii) and eq:propagators: D = i K^-1 in xi gauge."""
    smp = Sampler(24)
    s0 = brst.action
    p, g, delta, gamma = _rv("p"), _rv("g"), _rv("delta"), _rv("gamma")
    sp = _rv("sp")
    for f in ("Q", "U", "D", "L", "E"):
        legs = [make_leg(model, f + "bar", 1), make_leg(model, f, 2)]
        x = _rv("lor")(99_300)
        deltas = Expression.num(1)
        for (lb, sl), lf in zip(_slots(model, legs[0]), legs[1].labels):
            if sl.space.kind != "spinor":
                deltas = deltas * delta(lb, lf)
        _check(ev, model, s0, legs, p(2, x) * gamma(x, _spin(model, legs[0]), _spin(model, legs[1])) * deltas, smp)
    legs = [make_leg(model, "Hdag", 1), make_leg(model, "H", 2)]
    _check(ev, model, s0, legs, (sp(2, 2) - sym("m2")) * delta(legs[0].labels[0], legs[1].labels[0]), smp)
    for r in ("G", "W", "B"):
        legs = [make_leg(model, r, 1), make_leg(model, r, 2)]
        m1, m2 = _lor(model, legs[0]), _lor(model, legs[1])
        a, b = _adj(model, legs[0]), _adj(model, legs[1])
        dab = delta(a, b) if a is not None else Expression.num(1)
        _check(ev, model, s0, legs, -(sp(1, 1) * g(m1, m2) - p(1, m1) * p(1, m2)) * dab, smp)
    checked = set()
    for prop in propagators(s0, model):
        x, y = prop.fields
        lx, ly = prop.labels
        pol = {}
        for labs, name in ((lx, x), (ly, y)):
            for lab, sl in zip(labs, model.fields[name].slots):
                pol[ev.label_key(lab)] = smp.vector(ev.dim(ev.space(lab)))
        k = smp.vector(D)
        momenta = {0: k}
        val = ev.value(prop.value, pol, momenta)
        kinds = (model.fields[x].kind, model.fields[y].kind)
        if kinds not in (("fermion", "fermion_bar"), ("ghost", "antighost"), ("scalar", "scalar_conj"),
                         ("gauge", "gauge")):
            continue                                # <A b>, <b A>, <b b> (b is kept as a field)
        k2 = int(np.sum(k * ev.metric * k) % P)
        inv_k2 = pow(k2, -1, P)
        group = [(lab, sl) for lab, sl in zip(lx, model.fields[x].slots) if sl.space.kind not in ("spinor", "lorentz")]
        group_y = [lab for lab, sl in zip(ly, model.fields[y].slots) if sl.space.kind not in ("spinor", "lorentz")]
        deltas = [(a, b) for (a, _), b in zip(group, group_y)]
        dval = 1
        for a, b in deltas:
            dval = dval * int(np.sum(pol[ev.label_key(a)] * pol[ev.label_key(b)]) % P) % P
        if kinds == ("fermion", "fermion_bar"):
            sx, sy = lx[0], ly[0]
            kslash = sum(ev.gammas[m] * (ev.metric[m] * k[m] % P) for m in range(D)) % P
            mat = (IU * kslash * inv_k2) % P
            exp = int(pol[ev.label_key(sx)].dot(mat).dot(pol[ev.label_key(sy)]) % P) * dval % P
        elif kinds == ("ghost", "antighost"):
            exp = IU * inv_k2 * dval % P
        elif kinds == ("scalar", "scalar_conj"):
            exp = IU * pow((k2 - ev.params["m2"]) % P, -1, P) * dval % P
        elif kinds == ("gauge", "gauge"):
            ex, ey = pol[ev.label_key(lx[0])], pol[ev.label_key(ly[0])]
            gxy = int(np.sum(ex * ev.metric * ey) % P)
            kx, ky = int(np.sum(ex * k) % P), int(np.sum(ey * k) % P)
            xi = ev.params["xi"]
            exp = -IU * inv_k2 * (gxy - (1 - xi) * kx * ky * inv_k2) % P * dval % P
        assert val == exp, prop.fields
        checked.add(kinds)
    assert len(checked) == 4
