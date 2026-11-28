"""The one-loop back-end (hepkit_renorm.loops): gates V5/V6 on hand-built graph records.

References (golden data, see tests/data/README.md):

* ``loops_toolkit_golden.json``: N0 (and P N1 of the chiral cases) of the stand-alone toolkit ``bmhvkit.loop.uv_pole``
  for QED-like and chiral integrands (custom vertex functions on the fields ``B``, ``E`` of the one-family model).
* The class results of the report (``SCRATCHPAD/classes/<key>/final_result.json``), typed in below with their origin:
  g0_VV (pole of the gauge two-point function including the evanescent ``phat^2 gbar`` term, G and W factors, every
  topology), g1_cV (pole and N1 of the omega A breaking, one multiplet), g1_cFF (pole proportional to the Delta_1
  vertex and the finite part P N1 of the omega Psibar Psi breaking, multiplet Q, ghost of SU(3), per topology and
  summed), g1_cVVV_cVVVV (no pole and no finite part of omega A A A A). Units: ``i/(16 pi^2 eps)``; ``Tg = g_r T``
  and ``C_ad = g_r^2 N``, ``S_f`` as in tab:invariants.

Every comparison is exact (fingerprints over GF(p)) in D = 4 and, for the d-dimensional structure, in D = 6. All
tests together take about 20 s, so none is marked slow.
"""

from __future__ import annotations

import json
from fractions import Fraction

import numpy as np
import pytest
from conftest import DATA
from symbolica import Expression

from hepkit_renorm import rules, standard_legs
from hepkit_renorm import tensors as tn
from hepkit_renorm.brst import BRST
from hepkit_renorm.fingerprint import Fingerprinter, Representation
from hepkit_renorm.graphs import ExternalEnd, GraphPropagator, GraphRecord, GraphVertex, PropagatorEnd
from hepkit_renorm.loops import (
    Laurent,
    OneLoopBackend,
    _OddRepresentation,
    _with_dimension,
    embed_samples,
    evaluate_graphs,
    local_uv,
)

E, Pe = ExternalEnd, PropagatorEnd
I = tn.I
BACKEND = OneLoopBackend()


def _q(x) -> Expression:
    x = Fraction(x)
    return Expression.num(x.numerator) / x.denominator


@pytest.fixture(scope="module")
def setup(sm):
    brst = BRST(sm)
    action = brst.action
    props = {p.fields: GraphPropagator.from_rule(p, {sm.gauge_fixing.xi: Expression.num(1)})
             for p in rules.propagators(action, sm)}
    return rules.tree_rules(action, sm), rules.tree_rules(brst.breaking, sm), props


def _vertex(setup, cls, ends, delta=False):
    tree, breaking, _ = setup
    return GraphVertex.from_rule((breaking if delta else tree)[tuple(cls)], ends, dimension=5 if delta else 4)


def _values(sm, graph, expected0=None, expected1=None, dims=(4, 6), n=4, n1_dims=(4,)):
    """Compare N0 (all dims) and N1 (``n1_dims``) with expected expressions (``0`` for zero); returns the Laurent
    data per dimension."""
    pg = BACKEND.prepare(graph, sm)
    out = {}
    for D in dims:
        fp = Fingerprinter(sm, dimension=D, seed=7)
        smp = fp.samples([leg.field for leg in graph.legs], n)
        lau = BACKEND.evaluate(pg, fp, smp)
        for exp, got, use in ((expected0, lau.N0, True), (expected1, lau.N1, D in n1_dims)):
            if exp is None or not use:
                continue
            want = np.zeros(n, dtype=np.int64) if isinstance(exp, int) else fp.evaluate(exp, smp)
            assert np.array_equal(got, want), (graph.name, D)
        out[D] = (fp, smp, lau)
    return out


# ---------------------------------------------------------------------------------------------- representation
@pytest.mark.parametrize("dim", [5, 7, 9])
def test_odd_representation(field, dim):
    rep = _OddRepresentation(field, dim)
    p = field.p
    g = rep.gamma.astype(object)
    eta = rep.metric
    for a in range(dim):
        for b in range(dim):
            anti = (g[a].dot(g[b]) + g[b].dot(g[a])) % p
            want = 2 * eta[a] * (a == b) % p
            assert np.array_equal(anti, want * np.eye(rep.spinor_dim, dtype=object) % p)
    g5 = rep.gamma5.astype(object)
    for a in range(dim):
        comm = (g5.dot(g[a]) - g[a].dot(g5)) % p if a >= 4 else (g5.dot(g[a]) + g[a].dot(g5)) % p
        assert not comm.any()
    assert rep.spinor_dim == Representation(field, dim - 1).spinor_dim


def test_embedded_samples_keep_tree_values(sm, setup):
    """A vertex function without internal summations has the same value on embedded samples."""
    rule = setup[0][("W", "Q", "Qbar")]
    fp = Fingerprinter(sm, dimension=6, seed=2)
    smp = fp.samples(list(rule.fields), 3)
    for dim in (7, 8, 9):
        big = _with_dimension(fp, dim)
        assert np.array_equal(fp.evaluate(rule.value, smp), big.evaluate(rule.value, embed_samples(smp, fp, big)))


def test_laurent_operators(field):
    one = np.array([1, 2], dtype=np.int64)
    lau = Laurent(2, field, 2, {(0, 2): one, (1, 2): 3 * one, (2, 2): 5 * one, (1, 1): 7 * one})
    assert set(lau.K()) == {2, 1} and np.array_equal(lau.K()[1], 3 * one)      # eta/eps^2 -> 1/eps
    assert np.array_equal(lau.K()[2], one)
    assert set(lau.Kf()) == {(0, 2), (1, 2), (2, 2), (1, 1)}
    ke = lau.Ke()
    assert np.array_equal(ke[0], 12 * one % field.p) and np.array_equal(ke[1], 3 * one)
    total = lau + lau.scale(-1)
    assert all(not v.any() for v in total.coefficients.values())
    with pytest.raises(ValueError):
        _ = lau.N0


# ---------------------------------------------------------------------------------------------- toolkit cases
def _toolkit_cases(sm, setup):
    props = setup[2]
    pE, pB = props[("E", "Ebar")], props[("B", "B")]
    (b0,), (b1,) = pB.labels
    mB = GraphPropagator(("B", "B"), pB.labels, -I * tn.g(b0, b1), ((sm.param("m2"), 1),))
    mB2 = GraphPropagator(("B", "B"), pB.labels, -I * sm.param("m2") * tn.g(b0, b1), ((sm.param("m2"), 2),))

    def vertex(kind, ends):
        legs = tuple(standard_legs(sm.fields, ["B", "Ebar", "E"]))
        (lb,), (sbar,), (se,) = [leg.labels for leg in legs]
        if kind == "vec":
            val = tn.gamma(lb, sbar, se)
        else:
            val = tn.gammabar(lb, sbar, tn.spin(7)) * tn.PL(tn.spin(7), se)
        return GraphVertex(legs, val, (ends[0], ends[2], ends[1]))       # ends given as (B, E, Ebar)

    def rec(fields, verts, props_):
        return GraphRecord(tuple(standard_legs(sm.fields, fields)), tuple(verts), tuple(props_))

    out = {}
    for name, kind, photon in (("qed_self_energy", "vec", pB), ("chiral_self_energy", "chi", pB),
                               ("qed_self_energy_massive_photon", "vec", mB),
                               ("qed_self_energy_squared_massive_photon", "vec", mB2)):
        out[name] = rec(["Ebar", "E"], [vertex(kind, [Pe(1, 0), E(1), Pe(0, 1)]),
                                        vertex(kind, [Pe(1, 1), Pe(0, 0), E(0)])], [pE, photon])
    for name, kind in (("qed_vacuum_polarisation", "vec"), ("chiral_vacuum_polarisation", "chi")):
        out[name] = rec(["B", "B"], [vertex(kind, [E(0), Pe(0, 0), Pe(1, 1)]), vertex(kind, [E(1), Pe(1, 0), Pe(0, 1)])],
                        [pE, pE])
    for name, kind in (("qed_vertex", "vec"), ("chiral_vertex", "chi")):
        out[name] = rec(["B", "Ebar", "E"], [vertex(kind, [Pe(2, 0), E(2), Pe(0, 1)]),
                                             vertex(kind, [E(0), Pe(0, 0), Pe(1, 1)]),
                                             vertex(kind, [Pe(2, 1), Pe(1, 0), E(1)])], [pE, pE, pB])
    for name, kind in (("qed_box", "vec"), ("chiral_box", "chi")):
        out[name] = rec(["B"] * 4, [vertex(kind, [E(0), Pe(0, 0), Pe(1, 1)]), vertex(kind, [E(1), Pe(1, 0), Pe(2, 1)]),
                                    vertex(kind, [E(2), Pe(2, 0), Pe(3, 1)]), vertex(kind, [E(3), Pe(3, 0), Pe(0, 1)])],
                        [pE] * 4)
    return out


def test_toolkit_cases(sm, setup):
    with open(DATA / "loops_toolkit_golden.json") as fh:
        golden = json.load(fh)["cases"]
    graphs = _toolkit_cases(sm, setup)
    assert set(graphs) == set(golden)
    for name, graph in graphs.items():
        ref = golden[name]
        e0 = Expression.parse(ref["N0"])
        e1 = Expression.parse(ref["PN1"]) if ref["PN1"] else None
        _values(sm, graph, e0, e1)


# ---------------------------------------------------------------------------------------------- g0_VV
def _vv_graphs(sm, setup, X):
    props = setup[2]
    legs = tuple(standard_legs(sm.fields, [X, X]))
    out = []

    def V(cls, ends):
        return _vertex(setup, cls, ends)

    for f in {"G": ["Q", "U", "D"], "W": ["Q", "L"]}[X]:
        fb = f + "bar"
        out.append(GraphRecord(legs, (V((X, f, fb), [E(0), Pe(0, 0), Pe(1, 1)]), V((X, f, fb), [E(1), Pe(1, 0), Pe(0, 1)])),
                               (props[(f, fb)],) * 2, name=f"fermion {f}"))
    out.append(GraphRecord(legs, (V((X, X, X), [E(0), Pe(0, 0), Pe(1, 0)]), V((X, X, X), [E(1), Pe(0, 1), Pe(1, 1)])),
                           (props[(X, X)],) * 2, symmetry_factor=Fraction(1, 2), name="gauge bubble"))
    out.append(GraphRecord(legs, (V((X, X, X, X), [E(0), E(1), Pe(0, 0), Pe(0, 1)]),), (props[(X, X)],),
                           symmetry_factor=Fraction(1, 2), name="gauge tadpole"))
    gh, agh = "gh_" + X, "agh_" + X
    out.append(GraphRecord(legs, (V((X, gh, agh), [E(0), Pe(0, 0), Pe(1, 1)]), V((X, gh, agh), [E(1), Pe(1, 0), Pe(0, 1)])),
                           (props[(gh, agh)],) * 2, name="ghost loop"))
    if X == "W":
        out.append(GraphRecord(legs, (V((X, "H", "Hdag"), [E(0), Pe(0, 0), Pe(1, 1)]),
                                      V((X, "H", "Hdag"), [E(1), Pe(1, 0), Pe(0, 1)])),
                               (props[("H", "Hdag")],) * 2, name="Higgs bubble"))
        out.append(GraphRecord(legs, (V((X, X, "H", "Hdag"), [E(0), E(1), Pe(0, 0), Pe(0, 1)]),), (props[("H", "Hdag")],),
                               name="Higgs seagull"))
    return legs, out


@pytest.mark.parametrize("factor, bosonic", [("G", 5), ("W", Fraction(19, 6))])
def test_g0_vv_class(sm, setup, factor, bosonic):
    """g0_VV final_result.json: pole of Gamma_{A A} = g^2 delta^{ab} [c_d (p^2 g - p p) - 4/3 (pbar^2 gbar - pbar pbar)
    - 2/3 phat^2 gbar], c_d = 5 (G), 19/6 (W)."""
    legs, graphs = _vv_graphs(sm, setup, factor)
    (mu, a1), (nu, a2) = legs[0].labels, legs[1].labels
    x, y = tn.lor(1), tn.lor(2)
    pbmu, pbnu = tn.p(1, x) * tn.gbar(x, mu), tn.p(1, y) * tn.gbar(y, nu)
    expected = sm.coupling(factor) ** 2 * tn.delta(a1, a2) * (
        _q(bosonic) * (tn.sp(1, 1) * tn.g(mu, nu) - tn.p(1, mu) * tn.p(1, nu))
        - _q(Fraction(4, 3)) * (tn.spbar(1, 1) * tn.gbar(mu, nu) - pbmu * pbnu)
        - _q(Fraction(2, 3)) * tn.sphat(1, 1) * tn.gbar(mu, nu))
    prepared = [BACKEND.prepare(g, sm) for g in graphs]
    for D in (4, 6):
        fp = Fingerprinter(sm, dimension=D, seed=5)
        smp = {(factor, factor): fp.samples([factor, factor], 4)}
        values = evaluate_graphs(prepared, fp, smp)
        total = values[0]
        for v in values[1:]:
            total = total + v
        assert np.array_equal(total.N0, fp.evaluate(expected, smp[(factor, factor)]))
        for g, v in zip(graphs, values):
            if g.name.startswith("fermion"):
                assert not v.N1.any()               # four-dimensional fermion numerators: N1 = 0


# ---------------------------------------------------------------------------------------------- g1_cV
def test_g1_cv_omega_a(sm, setup):
    """g1_cV final_result.json (multiplet Q, factor G): N0 = (i/3) S_Q qhat^2 qbar^nu delta^ab, N1 = (i/3) S_Q q^2
    qbar^nu delta^ab, S_Q = g_s^2; the full N1 also agrees in D = 6 for this graph."""
    props = setup[2]
    legs = tuple(standard_legs(sm.fields, ["gh_G", "G"]))
    g = GraphRecord(legs, (_vertex(setup, ("gh_G", "Q", "Qbar"), [E(0), Pe(0, 0), Pe(1, 1)], delta=True),
                           _vertex(setup, ("G", "Q", "Qbar"), [E(1), Pe(1, 0), Pe(0, 1)])),
                    (props[("Q", "Qbar")],) * 2, name="cA Q")
    (a1,), (nu, a2) = legs[0].labels, legs[1].labels
    x = tn.lor(1)
    pre = I / 3 * sm.coupling("G") ** 2 * tn.p(1, x) * tn.gbar(x, nu) * tn.delta(a1, a2)
    _values(sm, g, pre * tn.sphat(1, 1), pre * (tn.spbar(1, 1) + tn.sphat(1, 1)), n1_dims=(4, 6))


# ---------------------------------------------------------------------------------------------- g1_cFF
def _cff(sm, setup):
    props = setup[2]
    legs = tuple(standard_legs(sm.fields, ["gh_G", "Qbar", "Q"]))
    (a,), (sb, ib, jb), (sq, iq, jq) = [leg.labels for leg in legs]
    gs, g2, gp = sm.coupling("G"), sm.coupling("W"), sm.coupling("B")
    yu, yd = sm.param("y_u"), sm.param("y_d")
    x, s0 = tn.lor(1), tn.spin(1)
    tg = gs * tn.T(a, ib, iq) * tn.delta(jb, jq)

    def sl(k, head, proj, sign=1):                 # sign * p(k).head proj, row Qbar, column Q
        return sign * tn.p(k, x) * head(x, sb, s0) * proj(s0, sq)

    q_b = sl(1, tn.gammabar, tn.PL)
    pb, ppb = sl(3, tn.gammabar, tn.PL), sl(2, tn.gammabar, tn.PL, -1)        # p = p3, p' = -p2
    ph_r, pph_l = sl(3, tn.gammahat, tn.PR), sl(2, tn.gammahat, tn.PL, -1)

    def V(cls, ends, delta=False):
        return _vertex(setup, cls, ends, delta)

    cad = 3 * gs ** 2
    graphs = {}
    for X, c in (("G", -gs ** 2 / 6), ("W", _q(Fraction(3, 4)) * g2 ** 2), ("B", gp ** 2 / 36)):
        graphs["a " + X] = (GraphRecord(legs, (V(("gh_G", "Q", "Qbar"), [E(0), Pe(0, 0), Pe(1, 1)], True),
                                               V((X, "Q", "Qbar"), [Pe(2, 0), E(2), Pe(0, 1)]),
                                               V((X, "Q", "Qbar"), [Pe(2, 1), Pe(1, 0), E(1)])),
                                        (props[("Q", "Qbar")],) * 2 + (props[(X, X)],)), 0, -I * c * tg * q_b)
    graphs["b U"] = (GraphRecord(legs, (V(("gh_G", "U", "Ubar"), [E(0), Pe(0, 0), Pe(1, 1)], True),
                                        V(("Q", "Ubar", "H"), [E(2), Pe(0, 1), Pe(2, 0)]),
                                        V(("Qbar", "U", "Hdag"), [E(1), Pe(1, 0), Pe(2, 1)])),
                                 (props[("U", "Ubar")],) * 2 + (props[("H", "Hdag")],)), 0, -I / 2 * yu ** 2 * tg * q_b)
    graphs["b D"] = (GraphRecord(legs, (V(("gh_G", "D", "Dbar"), [E(0), Pe(0, 0), Pe(1, 1)], True),
                                        V(("Q", "Dbar", "Hdag"), [E(2), Pe(0, 1), Pe(2, 1)]),
                                        V(("Qbar", "D", "H"), [E(1), Pe(1, 0), Pe(2, 0)])),
                                 (props[("D", "Dbar")],) * 2 + (props[("H", "Hdag")],)), 0, -I / 2 * yd ** 2 * tg * q_b)
    ghost = (props[("gh_G", "agh_G")], props[("Q", "Qbar")], props[("G", "G")])
    graphs["c1"] = (GraphRecord(legs, (V(("G", "gh_G", "agh_G"), [Pe(2, 0), E(0), Pe(0, 1)]),
                                       V(("gh_G", "Q", "Qbar"), [Pe(0, 0), E(2), Pe(1, 1)], True),
                                       V(("G", "Q", "Qbar"), [Pe(2, 1), Pe(1, 0), E(1)])), ghost),
                    -I / 2 * cad * tg * ph_r, I / 12 * cad * tg * (2 * pb - ppb))
    graphs["c2"] = (GraphRecord(legs, (V(("G", "gh_G", "agh_G"), [Pe(2, 0), E(0), Pe(0, 1)]),
                                       V(("gh_G", "Q", "Qbar"), [Pe(0, 0), Pe(1, 0), E(1)], True),
                                       V(("G", "Q", "Qbar"), [Pe(2, 1), E(2), Pe(1, 1)])), ghost),
                    I / 2 * cad * tg * pph_l, I / 12 * cad * tg * (pb - 2 * ppb))
    pole = -I * _q(Fraction(3, 2)) * gs ** 2 * tg * (ph_r - pph_l)
    finite = I * (27 * g2 ** 2 + 21 * gs ** 2 + gp ** 2 + 18 * yd ** 2 + 18 * yu ** 2) / 36 * tg * (pb - ppb)
    return graphs, pole, finite


def test_g1_cff_omega_psibar_psi(sm, setup):
    """g1_cFF final_result.json (Q, gluon ghost): per topology N0 and P N1, and the class sum: pole (C_ad/2) V_Delta1
    (evanescent, D = 6), finite part i (27 g^2 + 21 g_s^2 + g'^2 + 18 y_d^2 + 18 y_u^2)/36 (pbar - pbar')slash Tg P_L."""
    graphs, pole, finite = _cff(sm, setup)
    total = None
    for g, e0, e1 in graphs.values():
        got = _values(sm, g, e0, e1)
        total = got if total is None else {D: (o[0], o[1], total[D][2] + o[2]) for D, o in got.items()}
    for D, (fp, smp, lau) in total.items():
        assert np.array_equal(lau.N0, fp.evaluate(pole, smp))
        if D == 4:
            assert np.array_equal(lau.N1, fp.evaluate(finite, smp))


@pytest.mark.filterwarnings("ignore:This process .* is multi-threaded:DeprecationWarning")
def test_parallel_and_galois(sm, setup):
    """Forked workers give the same values; i -> -i gives the Galois conjugate (here: real parts are recovered)."""
    graphs, _, finite = _cff(sm, setup)
    prepared = [BACKEND.prepare(g, sm) for g, _, _ in graphs.values()]
    fp = Fingerprinter(sm, dimension=4, seed=9)
    smp = {("gh_G", "Qbar", "Q"): fp.samples(["gh_G", "Qbar", "Q"], 3)}
    serial = evaluate_graphs(prepared, fp, smp)
    parallel = evaluate_graphs(prepared, fp, smp, processes=2)
    assert all(np.array_equal(a.N0, b.N0) and np.array_equal(a.N1, b.N1) for a, b in zip(serial, parallel))
    conj = fp.replace(i_sign=-1)
    smp_c = {k: conj.samples(list(k), 3) for k in smp}
    total = None
    for pg in prepared:
        lau = local_uv(pg.graph, sm, conj, smp_c[("gh_G", "Qbar", "Q")])
        total = lau if total is None else total + lau
    assert np.array_equal(total.N1, conj.evaluate(finite, smp_c[("gh_G", "Qbar", "Q")]))


def test_rejects_two_loops(sm, setup):
    props = setup[2]
    legs = tuple(standard_legs(sm.fields, ["G", "G"]))
    q = props[("Q", "Qbar")]
    V = lambda ends: _vertex(setup, ("G", "Q", "Qbar"), ends)
    g = GraphRecord(legs, (V([E(0), Pe(0, 0), Pe(1, 1)]), V([E(1), Pe(3, 0), Pe(0, 1)]),
                           V([Pe(4, 0), Pe(1, 0), Pe(2, 1)]), V([Pe(4, 1), Pe(2, 0), Pe(3, 1)])),
                    (q, q, q, q, props[("G", "G")]))
    assert not BACKEND.supports(g)
    with pytest.raises(NotImplementedError):
        BACKEND.prepare(g, sm)


def test_five_point_breaking_vanishes(sm, setup):
    """g1_cVVV_cVVVV final_result.json: omega A A A A has no pole and no finite breaking; one E loop with four B."""
    props = setup[2]
    legs = tuple(standard_legs(sm.fields, ["gh_B", "B", "B", "B", "B"]))
    verts = [_vertex(setup, ("gh_B", "E", "Ebar"), [E(0), Pe(4, 0), Pe(0, 1)], delta=True)]
    verts += [_vertex(setup, ("B", "E", "Ebar"), [E(j), Pe(j - 1, 0), Pe(j, 1)]) for j in range(1, 5)]
    g = GraphRecord(legs, tuple(verts), (props[("E", "Ebar")],) * 5)
    _values(sm, g, 0, 0, dims=(4, 6), n=3)
