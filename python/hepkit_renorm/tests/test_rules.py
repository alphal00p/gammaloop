"""Gate V1: the tree-level Feynman rules of the one-family Standard Model.

The rules are derived (``brst.BRST(model).action`` = S_0 with gauge fixing, ghosts and BRST sources;
``BRST(model).breaking`` = Delta_1 = s S_0) and compared

* with every rule of report 03_setup sec:feynman-rules (eq:propagators to eq:V-delta) and with the tree vertex
  functions of report 02_method sec:m-ii, typed in below from the report;
* with the 174 vertex rules of the hand-written HepKIT model (golden digest in ``tests/data``), by their values at
  random kinematics, polarisations and parameters, component by component in SU(2) (the encoding of that model),
  in both directions (every hand rule is reproduced, every non-zero component of ours has a hand rule).

The sign conventions a user needs are pinned here: leg order of odd fields, outgoing momenta, the propagator
of an odd pair, the matrix-row convention.
"""

from __future__ import annotations

import itertools
import json

import pytest
from component_values import (
    SPINOR_DIM,
    ComponentEvaluator,
    Kinematics,
    hand_rule_expression,
    model_parameters,
    unit_vector,
)
from symbolica import Expression, S

from hepkit_renorm import classical_action
from hepkit_renorm import tensors as tn
from hepkit_renorm.brst import BRST
from hepkit_renorm.rules import PROPAGATOR_MOMENTUM, propagators, rule_table, tree_rules, vertex_rule

I = tn.I
FERMIONS = ("Q", "U", "D", "L", "E")
_W = [S(n) for n in ("ca_", "ci_", "cj_")]


# ---------------------------------------------------------------------------------------------- fixtures
@pytest.fixture(scope="module")
def brst(sm):
    return BRST(sm)


@pytest.fixture(scope="module")
def s0(brst):
    return brst.action


@pytest.fixture(scope="module")
def delta1(brst):
    return brst.breaking


@pytest.fixture(scope="module")
def rules0(sm, s0):
    return tree_rules(s0, sm)


@pytest.fixture(scope="module")
def rules_delta(sm, delta1):
    return tree_rules(delta1, sm)


@pytest.fixture(scope="module")
def ev(sm, field):
    return ComponentEvaluator(sm, field)


@pytest.fixture(scope="module")
def same(sm, ev, field):
    """``same(a, b, legs)``: equal values at three random points (all free labels of ``legs`` polarised)."""
    names = model_parameters(sm)

    def check(a: Expression, b: Expression, legs, fixed=None) -> bool:
        for seed in range(3):
            kin = Kinematics(field, len(legs), names, seed)
            pol = polarisations(sm, ev, kin, legs, fixed)
            if ev.value(a, pol, kin.momenta, kin.params) != ev.value(b, pol, kin.momenta, kin.params):
                return False
        return True

    return check


def polarisations(sm, ev, kin, legs, fixed=None):
    """Random polarisation vectors for every label of every leg; ``fixed[k] = c`` puts the SU(2) slot of leg
    ``k`` (0-based) on the basis vector ``c``."""
    pol = {}
    for k, leg in enumerate(legs):
        for lab, slot in zip(leg.labels, sm.fields[leg.field].slots):
            n = ev.dim(lab)
            c = None if fixed is None else fixed[k]
            pol[lab.to_canonical_string()] = (unit_vector(ev.F, n, c) if slot.factor == "W" and c is not None
                                              else kin.vector(n))
    return pol


# ---------------------------------------------------------------------------------------------- report objects
def group_labels(sm, leg):
    return [lab for lab, s in zip(leg.labels, sm.fields[leg.field].slots) if s.space.kind not in ("lorentz", "spinor")]


def spinor_label(sm, leg):
    return next(lab for lab, s in zip(leg.labels, sm.fields[leg.field].slots) if s.space.kind == "spinor")


def lorentz_label(sm, leg):
    return next(lab for lab, s in zip(leg.labels, sm.fields[leg.field].slots) if s.space.kind == "lorentz")


def adjoint_label(sm, leg):
    labs = group_labels(sm, leg)
    return labs[0] if labs else None


def num(q):
    return Expression.num(q.numerator) / q.denominator


def Tg(sm, multiplet, factor, a, row, col):
    """Report eq:Tg: ``(cal(T)^a_f)_{row col}`` of the multiplet (``f`` or ``H``) for an adjoint label ``a`` of
    ``factor`` (``None`` for U(1)); ``row``/``col`` are the group labels of the barred / unbarred leg."""
    info = sm.fields[multiplet]
    out = sm.coupling(factor)
    if sm.factors[factor].abelian:
        if factor not in info.charges:
            return Expression.num(0)
        out = out * num(info.charges[factor])
    elif not any(s.factor == factor for s in info.slots):
        return Expression.num(0)
    gslots = [s for s in info.slots if s.space.kind not in ("lorentz", "spinor")]
    for slot, r, c in zip(gslots, row, col):
        out = out * (tn.T(a, r, c) if slot.factor == factor else tn.delta(r, c))
    return out


def projector(sm, f, physical=True):
    left = sm.fermions[f].chirality == "L"
    return tn.PL if left == physical else tn.PR


def slash(head, k, s, t):
    mu = tn.fresh(tn.LORENTZ)
    return head(mu, s, t) * tn.p(k, mu)


def factors_of(sm, multiplet):
    info = sm.fields[multiplet]
    return [n for n, fa in sm.factors.items()
            if (n in info.charges) or (not fa.abelian and any(s.factor == n for s in info.slots))]


# ---------------------------------------------------------------------------------------------- the rule set
def test_rule_set_content(sm, rules0, rules_delta, golden_rules):
    """Multiplet-level vertex classes: 59 of S_0, 10 of Delta_1 (one per fermion multiplet and gauge factor),
    and exactly the families of the 174 hand-written rules."""
    assert len(rules0) == 59 and len(rules_delta) == 10
    assert {c for c in rules_delta} == {tuple(sorted(("gh_" + r, f, f + "bar"), key=sm.fields.order))
                                        for f in FERMIONS for r in factors_of(sm, f)}
    hand = {tuple(sorted((hand_particle(p)[0] for p in r["particles"]), key=sm.fields.order))
            for r in golden_rules["vertex_rules"]}
    assert set(rules0) | set(rules_delta) == hand
    for cls, rule in rules0.items():
        assert rule.fields == cls and rule.eliminated == len(cls) and not rule.outgoing


def test_rules_are_deterministic_and_serialisable(sm, s0, rules0):
    again = tree_rules(s0, sm)
    assert [r.value.to_canonical_string() for r in again.values()] == \
        [r.value.to_canonical_string() for r in rules0.values()]
    table = rule_table(rules0)
    assert json.loads(json.dumps(table)) == table
    assert table[0]["legs"][0]["momentum"] == 1


# ---------------------------------------------------------------------------------------------- propagators
def test_propagators(sm, s0):
    """eq:propagators, with xi kept: -i (g - (1 - xi) k k / k^2) delta / k^2 (Feynman gauge: -i g delta / k^2),
    i kslash / k^2 for the Dirac field (sterile partner included, d-dimensional), +i delta / k^2 for omega omegabar
    with the antighost first in V, i delta / (k^2 - m^2) for H H^dagger."""
    simp = sm.simplifier()

    def zero(e):
        return simp(e).together() == 0
    props = {p.fields: p for p in propagators(s0, sm)}
    k = PROPAGATOR_MOMENTUM
    k2, xi = tn.sp(k, k), sm.gauge_fixing.xi
    expected = {("Q", "Qbar"), ("U", "Ubar"), ("D", "Dbar"), ("L", "Lbar"), ("E", "Ebar"), ("H", "Hdag"),
                ("Hdag", "H")}
    for r in sm.factors:
        expected |= {(r, r), (r, "b_" + r), ("b_" + r, r), ("gh_" + r, "agh_" + r)}
    assert set(props) == expected

    def deltas(labs_a, labs_b, internal):
        out = Expression.num(1)
        for x, y in zip(labs_a, labs_b):
            if tn.label_space(x) in internal:
                out = out * tn.delta(x, y)
        return out

    for r in sm.factors:
        p = props[(r, r)]
        (mu, *a), (nu, *b) = p.labels
        d = deltas(a, b, {tn.label_space(x) for x in a})
        ref = -I * (tn.g(mu, nu) - (1 - xi) * tn.p(k, mu) * tn.p(k, nu) / k2) / k2 * d
        assert zero(p.value - ref)
        assert zero(p.value.replace(xi, 1) + I * tn.g(mu, nu) / k2 * d)
        pg = props[("gh_" + r, "agh_" + r)]
        a, b = pg.labels
        assert zero(pg.value - I / k2 * deltas(a, b, {tn.label_space(x) for x in a}))
        # vector-auxiliary entries of the linear gauge: <A b> = -k^mu / k^2, <b A> = k^nu / k^2, <b b> = 0
        (mu, *a), b = props[(r, "b_" + r)].labels
        assert zero(props[(r, "b_" + r)].value + tn.p(k, mu) / k2 * deltas(a, b, {tn.label_space(x) for x in a}))
    for f in FERMIONS:
        p = props[(f, f + "bar")]
        (s, *a), (t, *b) = p.labels
        d = deltas(a, b, {tn.label_space(x) for x in a})
        assert zero(p.value - I * slash(tn.gamma, k, s, t) / k2 * d)
        # eq:phys-prop: P_c (i kslash / k^2) P_cbar = i kslashbar P_cbar / k^2; the sterile part P_c D P_c is hatted
        u, v = tn.fresh(tn.SPINOR), tn.fresh(tn.SPINOR)
        pc, pcb = projector(sm, f), projector(sm, f, physical=False)
        sandwiched = pc(u, s) * p.value * pcb(t, v)
        w = tn.fresh(tn.SPINOR)
        assert zero(sandwiched - I * slash(tn.gammabar, k, u, w) * pcb(w, v) / k2 * d)
        sterile = pc(u, s) * p.value * pc(t, v)
        assert zero(sterile - I * slash(tn.gammahat, k, u, w) * pc(w, v) / k2 * d)
    p = props[("H", "Hdag")]
    (i,), (j,) = p.labels
    assert zero(p.value - I * tn.delta(i, j) / (k2 - sm.param("m2")))


def test_singular_quadratic_form_without_gauge_fixing(sm):
    with pytest.raises(ValueError, match="singular"):
        propagators(classical_action(sm), sm)


# ---------------------------------------------------------------------------------------------- report rules
def test_gauge_and_fermion_vertices(sm, s0, same):
    """eq:V-ffv, eq:V-vvv, eq:V-vvvv (legs all incoming, in the order written)."""
    g = tn.g
    for f in FERMIONS:
        for r in factors_of(sm, f):
            rule = vertex_rule(s0, [f + "bar", f, r], sm.simplifier())
            lb, l, la = rule.legs
            u = tn.fresh(tn.SPINOR)
            ref = I * tn.gammabar(lorentz_label(sm, la), spinor_label(sm, lb), u) * projector(sm, f)(
                u, spinor_label(sm, l)) * Tg(sm, f, r, adjoint_label(sm, la) if not sm.factors[r].abelian else None,
                                             group_labels(sm, lb), group_labels(sm, l))
            assert same(rule.feynman_rule, ref, rule.legs), (f, r)
    for r, fa in sm.factors.items():
        if fa.abelian:
            continue
        c = sm.coupling(r)
        rule = vertex_rule(s0, [r] * 3, sm.simplifier(), eliminate=None)
        (m1, a1), (m2, a2), (m3, a3) = rule.labels()
        p = tn.p
        ref = c * tn.f(a1, a2, a3) * (g(m1, m2) * (p(1, m3) - p(2, m3)) + g(m2, m3) * (p(2, m1) - p(3, m1))
                                      + g(m3, m1) * (p(3, m2) - p(1, m2)))
        assert same(rule.feynman_rule, ref, rule.legs)
        rule = vertex_rule(s0, [r] * 4, sm.simplifier())
        (m1, a), (m2, b), (m3, c_), (m4, e) = rule.labels()
        x = tn.fresh(fa.adjoint_space())

        def ff(i, j, k, l_, c=c, x=x):
            return c**2 * tn.f(i, j, x) * tn.f(k, l_, x)
        ref = -I * (ff(a, b, c_, e) * (g(m1, m3) * g(m2, m4) - g(m1, m4) * g(m2, m3))
                    + ff(a, c_, b, e) * (g(m1, m2) * g(m3, m4) - g(m1, m4) * g(m2, m3))
                    + ff(a, e, b, c_) * (g(m1, m2) * g(m3, m4) - g(m1, m3) * g(m2, m4)))
        assert same(rule.feynman_rule, ref, rule.legs)


def test_ghost_vertex(sm, s0, same):
    """eq:V-ghost: omegabar^a (p' outgoing) A^b_mu omega^c: - fg^{abc} p'^mu, in the order written (the
    antighost first); the canonical order (A, omega, omegabar) differs by the exchange of the two ghosts."""
    for r, fa in sm.factors.items():
        if fa.abelian:
            continue
        rule = vertex_rule(s0, ["agh_" + r, r, "gh_" + r], sm.simplifier(), outgoing=[1], eliminate=3)
        (a,), (mu, b), (c,) = rule.labels()
        ref = -sm.coupling(r) * tn.f(a, b, c) * tn.p(1, mu)
        assert same(rule.feynman_rule, ref, rule.legs)
        assert same(rule.reordered([1, 2, 0]).feynman_rule, -ref, rule.legs)  # (A, omega, omegabar): ghosts swapped


def test_higgs_and_yukawa_vertices(sm, s0, same):
    """eq:V-ssv (H_j incoming p, H^dagger_i outgoing p'), eq:V-ssss, eq:V-yuk."""
    simp = sm.simplifier()
    for r in factors_of(sm, "H"):
        rule = vertex_rule(s0, ["Hdag", "H", r], simp, outgoing=[1], eliminate=3)
        lh, l, la = rule.legs
        a = adjoint_label(sm, la) if not sm.factors[r].abelian else None
        mu = lorentz_label(sm, la)
        ref = I * (tn.p(2, mu) + tn.p(1, mu)) * Tg(sm, "H", r, a, group_labels(sm, lh), group_labels(sm, l))
        assert same(rule.feynman_rule, ref, rule.legs), r
    for r1, r2 in itertools.combinations_with_replacement(factors_of(sm, "H"), 2):
        rule = vertex_rule(s0, ["Hdag", "H", r1, r2], simp)
        lh, l, l1, l2 = rule.legs
        a = adjoint_label(sm, l1) if not sm.factors[r1].abelian else None
        b = adjoint_label(sm, l2) if not sm.factors[r2].abelian else None
        k = [tn.fresh(sm.fields["H"].slots[0].space)]
        row, col = group_labels(sm, lh), group_labels(sm, l)
        anti = (Tg(sm, "H", r1, a, row, k) * Tg(sm, "H", r2, b, k, col)
                + Tg(sm, "H", r2, b, row, k) * Tg(sm, "H", r1, a, k, col))
        ref = I * tn.g(lorentz_label(sm, l1), lorentz_label(sm, l2)) * anti
        assert same(rule.feynman_rule, ref, rule.legs), (r1, r2)
    rule = vertex_rule(s0, ["Hdag", "H", "Hdag", "H"], simp)
    (i,), (j,), (k,), (l_,) = rule.labels()
    ref = -2 * I * sm.param("lambda") * (tn.delta(i, j) * tn.delta(k, l_) + tn.delta(i, l_) * tn.delta(k, j))
    assert same(rule.feynman_rule, ref, rule.legs)
    yu, yd, ye = (sm.param(n) for n in ("y_u", "y_d", "y_e"))
    cases = [(["Qbar", "U", "Hdag"], -I * yu, "eps", tn.PR), (["Ubar", "Q", "H"], -I * yu, "eps", tn.PL),
             (["Qbar", "D", "H"], -I * yd, "delta", tn.PR), (["Dbar", "Q", "Hdag"], -I * yd, "delta", tn.PL),
             (["Lbar", "E", "H"], -I * ye, "delta", tn.PR), (["Ebar", "L", "Hdag"], -I * ye, "delta", tn.PL)]
    for names, c, kind, proj in cases:
        rule = vertex_rule(s0, names, simp)
        l1, l2, l3 = rule.legs
        g1, g2, (k,) = group_labels(sm, l1), group_labels(sm, l2), group_labels(sm, l3)
        doublet = g1[-1] if names[0] in ("Qbar", "Lbar") else g2[-1]
        ref = c * proj(spinor_label(sm, l1), spinor_label(sm, l2)) * (tn.epsf(doublet, k) if kind == "eps"
                                                                     else tn.delta(doublet, k))
        if names[0][0] in "QUD":
            ref = ref * tn.delta(g1[0], g2[0])
        assert same(rule.feynman_rule, ref, rule.legs), names


def test_source_vertices(sm, s0, same):
    """eq:V-src-gauge and eq:V-src-matter, legs in the order written; the source is a leg with momentum and no
    propagator. Row index of a matrix: the barred leg (Rbar_f, Psibar_f, sigma^dagger, H^dagger)."""
    simp = sm.simplifier()
    for r, fa in sm.factors.items():
        rule = vertex_rule(s0, ["src_" + r, "gh_" + r], simp, eliminate=1)  # srcA^mu_a omega^b(p): delta^{ab} p_mu
        (mu, *a), b = rule.labels()
        ref = tn.p(2, mu) * (tn.delta(a[0], b[0]) if a else 1)
        assert same(rule.feynman_rule, ref, rule.legs), r
        if fa.abelian:
            continue
        c = sm.coupling(r)
        rule = vertex_rule(s0, ["src_" + r, r, "gh_" + r], simp)
        (mu, a), (nu, b), (c_,) = rule.labels()
        assert same(rule.feynman_rule, I * c * tn.f(a, b, c_) * tn.g(mu, nu), rule.legs)
        rule = vertex_rule(s0, ["src_gh_" + r, "gh_" + r, "gh_" + r], simp)
        (a,), (b,), (c_,) = rule.labels()
        assert same(rule.feynman_rule, -I * c * tn.f(a, b, c_), rule.legs)
    for f in FERMIONS:
        for r in factors_of(sm, f):
            ab = not sm.factors[r].abelian
            rule = vertex_rule(s0, ["src_" + f, "gh_" + r, f], simp)  # Rbar_f omega^a Psi_f: -cal(T)^a P_c
            lr, lg, l = rule.legs
            ref = -Tg(sm, f, r, adjoint_label(sm, lg) if ab else None, group_labels(sm, lr), group_labels(sm, l)) \
                * projector(sm, f)(spinor_label(sm, lr), spinor_label(sm, l))
            assert same(rule.feynman_rule, ref, rule.legs), (f, r)
            rule = vertex_rule(s0, ["gh_" + r, f + "bar", "src_" + f + "bar"], simp)  # omega Psibar R: -P_cbar T
            lg, lb, lr = rule.legs
            ref = -projector(sm, f, physical=False)(spinor_label(sm, lb), spinor_label(sm, lr)) \
                * Tg(sm, f, r, adjoint_label(sm, lg) if ab else None, group_labels(sm, lb), group_labels(sm, lr))
            assert same(rule.feynman_rule, ref, rule.legs), (f, r)
    for r in factors_of(sm, "H"):
        ab = not sm.factors[r].abelian
        rule = vertex_rule(s0, ["src_H", "gh_" + r, "H"], simp)  # sigma^dagger omega H: -cal(T)_H
        ls, lg, l = rule.legs
        ref = -Tg(sm, "H", r, adjoint_label(sm, lg) if ab else None, group_labels(sm, ls), group_labels(sm, l))
        assert same(rule.feynman_rule, ref, rule.legs), r
        rule = vertex_rule(s0, ["gh_" + r, "Hdag", "src_Hdag"], simp)  # omega H^dagger sigma: +cal(T)_H
        lg, lh, ls = rule.legs
        ref = Tg(sm, "H", r, adjoint_label(sm, lg) if ab else None, group_labels(sm, lh), group_labels(sm, ls))
        assert same(rule.feynman_rule, ref, rule.legs), r


def test_breaking_vertex(sm, delta1, rules_delta, same):
    """eq:V-delta: omega^a (ghost momentum p' - p eliminated), Psibar_f outgoing p', Psi_f incoming p:
    i V = cal(T)^a (phat-slash P_cbar - p'hat-slash P_c), V = -i cal(T)^a (...)."""
    simp = sm.simplifier()
    for f in FERMIONS:
        for r in factors_of(sm, f):
            rule = vertex_rule(delta1, ["gh_" + r, f + "bar", f], simp, outgoing=[2], eliminate=1)
            lg, lb, l = rule.legs
            s, t, u = spinor_label(sm, lb), spinor_label(sm, l), tn.fresh(tn.SPINOR)
            T = Tg(sm, f, r, adjoint_label(sm, lg) if not sm.factors[r].abelian else None, group_labels(sm, lb),
                   group_labels(sm, l))
            bracket = (slash(tn.gammahat, 3, s, u) * projector(sm, f, physical=False)(u, t)
                       - slash(tn.gammahat, 2, s, u) * projector(sm, f)(u, t))
            assert same(rule.feynman_rule, T * bracket, rule.legs), (f, r)
            assert same(rule.value, -I * T * bracket, rule.legs)
            # the rule of the rule set is the same vertex in canonical order (omega, Psi, Psibar), all incoming
            canon = rules_delta[tuple(sorted((lg.field, lb.field, l.field), key=sm.fields.order))]
            assert canon.fields == ("gh_" + r, f, f + "bar")


def test_tree_vertex_functions(sm, s0, same):
    """Report 02_method sec:m-ii, tree level (V, not i V), including Gamma_{rho omega} = -Gamma_{omega rho}."""
    simp = sm.simplifier()
    r = "G"
    rule = vertex_rule(s0, ["src_" + r, "gh_" + r], simp, eliminate=1)  # (-q, q): -i q_mu delta^{ae}
    (mu, e), (a,) = rule.labels()
    assert same(rule.value, -I * tn.p(2, mu) * tn.delta(a, e), rule.legs)
    swapped = vertex_rule(s0, ["gh_" + r, "src_" + r], simp, eliminate=2)  # Gamma_{omega rho} = -Gamma_{rho omega}
    (a,), (mu, e) = swapped.labels()
    assert same(swapped.value, I * tn.p(1, mu) * tn.delta(a, e), swapped.legs)
    rule = vertex_rule(s0, ["src_" + r, "gh_" + r, r], simp)  # fg^{eba} g_mu^nu
    (mu, e), (a,), (nu, b) = rule.labels()
    assert same(rule.value, sm.coupling(r) * tn.f(e, b, a) * tn.g(mu, nu), rule.legs)
    for f in ("Q", "E"):
        for rr in factors_of(sm, f):
            ab = not sm.factors[rr].abelian
            rule = vertex_rule(s0, [rr, f + "bar", f], simp)  # Gamma_{A Psibar Psi} = gbar^mu P_c cal(T)^e
            la, lb, l = rule.legs
            u = tn.fresh(tn.SPINOR)
            ref = tn.gammabar(lorentz_label(sm, la), spinor_label(sm, lb), u) * projector(sm, f)(
                u, spinor_label(sm, l)) * Tg(sm, f, rr, adjoint_label(sm, la) if ab else None, group_labels(sm, lb),
                                             group_labels(sm, l))
            assert same(rule.value, ref, rule.legs)
            rule = vertex_rule(s0, ["src_" + f, "gh_" + rr, f], simp)  # Gamma_{Rbar omega Psi} = i cal(T) P_c
            lr, lg, l = rule.legs
            ref = I * Tg(sm, f, rr, adjoint_label(sm, lg) if ab else None, group_labels(sm, lr),
                         group_labels(sm, l)) * projector(sm, f)(spinor_label(sm, lr), spinor_label(sm, l))
            assert same(rule.value, ref, rule.legs)
    rule = vertex_rule(s0, ["Qbar", "Q"], simp)  # Gamma_{Psibar Psi}(p) = pslash, p the incoming Psi momentum
    lb, l = rule.legs
    ref = slash(tn.gamma, 2, spinor_label(sm, lb), spinor_label(sm, l)) * tn.delta(
        group_labels(sm, lb)[0], group_labels(sm, l)[0]) * tn.delta(group_labels(sm, lb)[1], group_labels(sm, l)[1])
    assert same(rule.value, ref, rule.legs)
    rule = vertex_rule(s0, ["Hdag", "H"], simp)  # Gamma_{H^dag H}(p) = p^2 - m^2
    (i,), (j,) = rule.labels()
    assert same(rule.value, (tn.sp(2, 2) - sm.param("m2")) * tn.delta(i, j), rule.legs)


def test_tree_ward_identity(sm, s0, delta1, same):
    """eq:m-tree-WI: -i q_mu V^mu[A(q) Psibar(p') Psi(p)] = V_Delta1[omega(q) Psibar(p') Psi(p)]
    + i [cal(T) P_cbar V_2(p) - V_2(p') P_c cal(T)], with V_2(p) = pslash; it ties the signs of the gauge vertex,
    the kinetic term, the BRST transformations and the breaking vertex together."""
    simp = sm.simplifier()
    for f, r in (("Q", "G"), ("Q", "W"), ("L", "B"), ("U", "B")):
        ab = not sm.factors[r].abelian
        va = vertex_rule(s0, [r, f + "bar", f], simp, outgoing=[2], eliminate=1)
        vd = vertex_rule(delta1, ["gh_" + r, f + "bar", f], simp, outgoing=[2], eliminate=1)
        la, lb, l = va.legs
        lg = vd.legs[0]
        mu = lorentz_label(sm, la)
        lhs = -I * va.value * (tn.p(2, mu) - tn.p(3, mu))  # q = p' - p
        if ab:
            lhs = lhs.replace(adjoint_label(sm, la), adjoint_label(sm, lg))
        s, t, u = spinor_label(sm, lb), spinor_label(sm, l), tn.fresh(tn.SPINOR)
        T = Tg(sm, f, r, adjoint_label(sm, lg) if ab else None, group_labels(sm, lb), group_labels(sm, l))
        rest = I * T * (projector(sm, f, physical=False)(s, u) * slash(tn.gamma, 3, u, t)
                        - slash(tn.gamma, 2, s, u) * projector(sm, f)(u, t))
        assert same(lhs, vd.value + rest, vd.legs), (f, r)


def test_reordering_signs(sm, s0):
    """Exchanging two odd legs flips the sign; moving an even leg does not."""
    simp = sm.simplifier()
    rule = vertex_rule(s0, ["src_Q", "gh_G", "Q"], simp)
    assert rule.parities == (0, 1, 1)
    assert (rule.reordered([0, 2, 1]).value + rule.value).expand() == 0
    assert (rule.reordered([1, 0, 2]).value - rule.value).expand() == 0
    assert (rule.reordered([2, 1, 0]).value + rule.value).expand() == 0
    with pytest.raises(ValueError):
        rule.reordered([0, 0, 1])


# ---------------------------------------------------------------------------------------------- hand-written model
def hand_particle(name):
    """Particle of the hand-written model -> (field, SU(2) component or None). ``X~`` is the conjugate leg:
    ``ghX~`` the antighost, ``RX~`` = Rbar_f (source of Psi_f), ``sigma~`` = sigma^dagger (source of H)."""
    bar = name.endswith("~")
    base = name[:-1] if bar else name

    def comp(s):
        return int(s) - 1 if s else None
    if base.startswith("gh"):
        return ("agh_" if bar else "gh_") + base[2], comp(base[3:])
    if base.startswith("rho"):
        return "src_" + base[3], comp(base[4:])
    if base.startswith("zeta"):
        return "src_gh_" + base[4], comp(base[5:])
    if base.startswith("sigma"):
        return ("src_H" if bar else "src_Hdag"), comp(base[5:])
    if base[0] == "R" and len(base) > 1:
        return ("src_" + base[1] if bar else "src_" + base[1] + "bar"), comp(base[2:])
    if base[0] == "H":
        return ("Hdag" if bar else "H"), comp(base[1:])
    if base[0] in "GWB":
        return base[0], comp(base[1:])
    return (base[0] + "bar" if bar else base[0]), comp(base[1:])


def written_order(sm, fields):
    """The leg order of the report (in which the hand-written couplings are given) for the leg order of a hand
    rule: the ghost vertex is written omegabar A omega (eq:V-ghost), the conjugate matter source vertices
    omega Psibar R and omega H^dagger sigma (eq:V-src-matter); all other families in the hand order."""
    kinds = tuple(sm.fields[f].kind for f in fields)
    if kinds == ("ghost", "antighost", "gauge"):
        return [1, 2, 0]
    if len(kinds) == 3 and kinds[0] in ("fermion_bar", "scalar_conj") and kinds[1:] == ("ghost", "source"):
        return [1, 0, 2]
    return list(range(len(fields)))


HAND_PARAMS = {"gs": ("coupling", "G"), "gw": ("coupling", "W"), "gp": ("coupling", "B"), "yu": ("param", "y_u"),
               "yd": ("param", "y_d"), "ye": ("param", "y_e"), "lam": ("param", "lambda")}


def compare_hand_rules(sm, s0, delta1, ev, field, golden, *, transpose_colour=False, use_written_order=True):
    """Names of the hand rules whose value differs from ours (empty: all agree)."""
    simp = sm.simplifier()
    params = {k: (sm.coupling(v) if kind == "coupling" else sm.param(v)) for k, (kind, v) in HAND_PARAMS.items()}
    names = model_parameters(sm)
    G = sm.factors["G"]
    spaces = {"lorentz": tn.LORENTZ, "spinor": tn.SPINOR, "colour_adj": G.adjoint_space(),
              "colour_fund": G.rep_space("F")}
    bad = []
    for rule_json in golden["vertex_rules"]:
        pc = [hand_particle(p) for p in rule_json["particles"]]
        fields, comps = [f for f, _ in pc], [c for _, c in pc]
        breaking = any(o[0] == "DELTA" for row in rule_json["couplings"] for c in row if c for o in c["orders"])
        rule = vertex_rule(delta1 if breaking else s0, fields, simp, eliminate=None)
        value = rule.reordered(written_order(sm, fields)).value if use_written_order else rule.value
        legs = rule.legs
        dummies = {}

        def leg_label(k, space, legs=legs):
            leg = legs[k - 1]
            for lab, s in zip(leg.labels, sm.fields[leg.field].slots):
                kind = s.space.kind
                if (space == kind or space == "colour" and s.factor == "G" or
                        space == "colour_adj" and kind == "adjoint" and s.factor == "G" or
                        space == "colour_fund" and kind == "fund" and s.factor == "G"):
                    return lab
            raise KeyError((leg.field, space))

        def dummy_label(key, space, dummies=dummies):
            if key not in dummies:
                dummies[key] = tn.fresh(spaces[space])
            return dummies[key]
        hand = hand_rule_expression(rule_json, leg_label, dummy_label, params)
        if transpose_colour:  # the other reading of UFO's T(a, i, jb)
            hand = hand.replace(tn.T(*_W), tn.T(_W[0], _W[2], _W[1]))
        ours = (I * value).expand()
        for seed in range(2):
            kin = Kinematics(field, len(legs), names, seed)
            pol = polarisations(sm, ev, kin, legs, comps)
            if ev.value(ours, pol, kin.momenta, kin.params) != ev.value(hand, pol, kin.momenta, kin.params):
                bad.append(rule_json["name"])
                break
    return bad


def test_hand_written_rules(sm, s0, delta1, ev, field, golden_rules):
    """All 174 rules of the hand-written model: i V of the derived rule, legs in the report's written order,
    equals the hand rule (coupling x colour x Lorentz) at random kinematics, SU(2) component by component."""
    assert compare_hand_rules(sm, s0, delta1, ev, field, golden_rules) == []


def test_hand_written_rules_negative_controls(sm, s0, delta1, ev, field, golden_rules):
    """The comparison sees the leg-order sign of odd legs and the colour index order."""
    unsigned = compare_hand_rules(sm, s0, delta1, ev, field, golden_rules, use_written_order=False)
    ghost_vertices = [r["name"] for r in golden_rules["vertex_rules"]
                      if [hand_particle(p)[0][:3] for p in r["particles"]][:2] == ["gh_", "agh"]]
    conj_sources = [r["name"] for r in golden_rules["vertex_rules"]
                    if hand_particle(r["particles"][0])[0].endswith("bar")
                    and hand_particle(r["particles"][2])[0].startswith("src_")]
    assert len(ghost_vertices) == 7 and len(conj_sources) == 23
    assert sorted(unsigned) == sorted(ghost_vertices + conj_sources)
    transposed = compare_hand_rules(sm, s0, delta1, ev, field, golden_rules, transpose_colour=True)

    def has_quark(rule):
        return any(s.factor == "G" and s.space.kind == "fund"
                   for p in rule["particles"] for s in sm.fields[hand_particle(p)[0]].slots)
    assert transposed and all(has_quark(r) for r in golden_rules["vertex_rules"] if r["name"] in transposed)


def test_hand_written_rules_are_complete(sm, rules0, rules_delta, ev, field, golden_rules):
    """Every non-zero SU(2) component of every derived vertex has a hand rule, and conversely (unordered legs)."""
    hand = {tuple(sorted(hand_particle(p) for p in r["particles"])) for r in golden_rules["vertex_rules"]}
    assert len(hand) == len(golden_rules["vertex_rules"]) == 174
    names = model_parameters(sm)
    ours = set()
    for rules in (rules0, rules_delta):
        for rule in rules.values():
            ranges = []
            for leg in rule.legs:
                ws = [s for s in sm.fields[leg.field].slots if s.factor == "W"]
                ranges.append(range(ws[0].space.dim) if ws else [None])
            for comps in itertools.product(*ranges):
                values = []
                for seed in (11, 12):
                    kin = Kinematics(field, len(rule.legs), names, seed)
                    values.append(ev.value(rule.value, polarisations(sm, ev, kin, rule.legs, comps), kin.momenta,
                                           kin.params))
                if any(values):
                    ours.add(tuple(sorted(zip(rule.fields, comps), key=lambda x: (x[0], -1 if x[1] is None else x[1]))))
    hand = {tuple(sorted(h, key=lambda x: (x[0], -1 if x[1] is None else x[1]))) for h in hand}
    assert ours == hand


# ---------------------------------------------------------------------------------------------- evaluator
def test_component_representation(ev, field):
    """The explicit BMHV representation of the test evaluator: Clifford algebra in D = 6, gamma5 anticommutes
    with the four-dimensional and commutes with the evanescent gammas, gamma5 = i g^0 g^1 g^2 g^3, eps^{0123} = 1."""
    F = field
    g = ev.gamma
    for mu, nu in itertools.product(range(6), repeat=2):
        anti = F.reduce(F.matmul(g[mu], g[nu]) + F.matmul(g[nu], g[mu]))
        assert (anti == F.reduce(2 * ev.eta[mu, nu] * F.array([[int(i == j) for j in range(SPINOR_DIM)]
                                                                for i in range(SPINOR_DIM)]))).all()
    for mu in range(6):
        c = F.matmul(ev.gamma5, g[mu])
        d = F.matmul(g[mu], ev.gamma5)
        assert (F.reduce(c + d) == 0).all() if mu < 4 else (F.reduce(c - d) == 0).all()
    prod = F.matmul(F.matmul(g[0], g[1]), F.matmul(g[2], g[3]))
    assert (F.reduce(F.i * prod - ev.gamma5) == 0).all()
    assert ev.eps[0, 1, 2, 3] == 1 and ev.eps[1, 0, 2, 3] == F.p - 1
