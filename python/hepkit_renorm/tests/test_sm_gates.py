"""Gates V7 and V8 for the one-family Standard Model on the report's golden data (slow: about 10 minutes; run with
``pytest --runslow tests/test_sm_gates.py``).

Golden data: ``tests/data/final_counterterms.json`` (verbatim copy of ``SCRATCHPAD/phase4/final_counterterms.json``)
turned into local functionals by ``report_golden.py``: the singular counterterm (closed directions continued to d
dimensions plus the evanescent operators), the breaking ``X_1`` (``X.engine_names`` on the N monomials of
tab:Nbasis / ``st/bases_n.py``), the restoring counterterm (``r``) and the counterterm Feynman rules of report
05a. The algebraic modules (matching, st_solve, orders, the writer) take them through the same interface as loop
data (``FunctionalData``); the loop pipeline enters in ``test_sm_pipeline.py``.

V7: renormalisation constants and evanescent coefficients of deltaS^(1); closure of its four-dimensional limit;
the pole identity against eq:sts-bd-dS; b X_1 = 0 and no anomaly; the 39 non-zero r_alpha with (N1)-(N4) (by name
on the report's basis: rank 64 on 64 allowed monomials, 8 existence conditions, 5 of them testable in source
classes and 3 in the Yukawa classes); the same with the model-independent statement of (N1)-(N4) on the report's
and on the generated basis; the automatic fallback (N4') with a gh^B Yukawa component; the kinetic representative
and the literature one differ from the report's by closed directions; with the leptons removed the anomaly is
reported. V8: the vertices of deltaS^(1) and DeltaS^(1) agree with report 05a class by class, and the written
renormalised model loads in HepKIT with rules equal to the store entries.
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np
import pytest
from report_golden import REMOVED, CounterRules, Golden, anomaly_breaking, parse, report_basis
from test_matching import equal_external

from hepkit_renorm import Density, LocalFunctional, model_from_dict
from hepkit_renorm.baction import breaking_by_factor
from hepkit_renorm.brst import BRST
from hepkit_renorm.exact import default_field
from hepkit_renorm.fingerprint import Fingerprinter
from hepkit_renorm.matching import FunctionalData, OrderContext, match_breaking, match_singular, pole_identity
from hepkit_renorm.orders import OrderResult, counterterm_model, loop_factor, update_store
from hepkit_renorm.st_solve import (
    AnomalyError,
    ElementCondition,
    Normalisation,
    solve_restoring,
    standard_normalisation,
)
from hepkit_renorm.store import OrderStore, tree_store

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def brst(sm):
    return BRST(sm)


@pytest.fixture(scope="module")
def golden(sm):
    return Golden(sm)


@pytest.fixture(scope="module")
def q(golden):
    return golden.q_basis()


@pytest.fixture(scope="module")
def C(golden, golden_counterterms, q):
    return golden.counterterm(golden_counterterms, q)


@pytest.fixture(scope="module")
def X(golden, golden_counterterms):
    return golden.breaking(golden_counterterms)


@pytest.fixture(scope="module")
def DeltaS(golden, golden_counterterms, q):
    return golden.restoring(golden_counterterms, q)


@pytest.fixture(scope="module")
def ctx(sm, brst):
    return OrderContext(sm, brst=brst)


@pytest.fixture(scope="module")
def singular(ctx, brst, C):
    return match_singular(ctx, FunctionalData(brst, counterterm=C))


@pytest.fixture(scope="module")
def restoring(ctx, brst, X):
    return solve_restoring(ctx, FunctionalData(brst, breaking=X))


def expr(ctx, poly):
    return ctx.space.expression(poly)


# ---------------------------------------------------------------------------------------------- V7: deltaS
def test_v7_renormalisation_constants(sm, ctx, singular, golden_counterterms):
    z = golden_counterterms["delta_S"]["Z_factors"]["values_in_kappa_units"]
    P = lambda s: parse(sm, s)
    for r in "GWB":
        assert expr(ctx, singular.constants[f"Z[{r}]"]) == P(z["dZ_A"][r])
        assert expr(ctx, singular.constants[f"Z[gh_{r}]"]) == P(z["dZ_omega"][r])
        assert singular.relative(sm.factors[r].coupling) == P(z["dZ_g"][r])
    for f in "QUDLE":
        assert expr(ctx, singular.constants[f"Z[{f}]"]) == P(z["dZ_f"][f])
    assert expr(ctx, singular.constants["Z[H]"]) == P(z["dZ_H"])
    for k in "ude":
        assert expr(ctx, singular.constants[f"d[y_{k}]"]) == P(z["dy"][k])
    assert expr(ctx, singular.constants["d[lambda]"]) == P(z["dlambda"])
    assert expr(ctx, singular.constants["d[m2]"]) == P(z["dm2"])
    assert not singular.unmatched_directions


def test_v7_closure_and_evanescent_part(sm, ctx, singular, golden, golden_counterterms, q):
    # the four-dimensional limit is b-closed: all 84 independent relations of b on the 109 coefficients hold
    # (72 on the report's 97), and it equals the report's closed-direction form
    assert singular.closure_rank == 84
    four = singular.four_dimensional_functional()
    closed = golden.closed(q)
    names = {"O_A^G": "O_A_G", "O_c^G": "O_c_G", "O_g^G": "O_g_G", "O_A^W": "O_A_W", "O_c^W": "O_c_W",
             "O_g^W": "O_g_W", "Q_T^B": "T_B", "O_H": "O_H", "Q_lambda": "lam", "Q_m": "m"}
    total = Density(sm.fields)
    for key, value in golden_counterterms["delta_S"]["closed_direction_coefficients_4d"].items():
        name = names.get(key) or {f"O_psi^{f}": f"O_psi_{f}" for f in "QUDLE"}.get(key) or \
            {f"Q^+_y{k}": f"y{k}+" for k in "ude"}.get(key)
        c = parse(sm, value)
        if name and c.expand() != 0:
            total = total + closed[name].scale(c)
    assert equal_external(sm, four, LocalFunctional.integrate(total), (4,))
    # evanescent part: eight operators (E_A, Etilde_F per factor, E_H, Etilde_H), four b_d-closed
    assert len(singular.evanescent) == 8
    assert sorted(cl for _, cl in singular.evanescent.values()) == [False] * 4 + [True] * 4
    ours = LocalFunctional.integrate(sum((e.functional.scale(expr(ctx, singular.coefficients[e.name]))
                                          for e in ctx.singular.evanescent), Density(sm.fields)))
    assert equal_external(sm, ours, golden.evanescent(golden_counterterms), (6, 8))


def test_v7_pole_identity(sm, ctx, brst, singular, C, golden):
    """The ghost-number-1 poles predicted by eq:sts-bd-dS: b_d deltaS = -kappa[sum_(G,W) C_ad/2 Delta_1^(r) +
    sum_r S/3 b_d E_A^(r) + 2/3 Yt b_d E_H] (D = 6, all 168 classes)."""
    inv = {"C_ad": {"G": 3 * sm.param("g_s") ** 2, "W": 2 * sm.param("g") ** 2},
           "S": {"G": 2 * sm.param("g_s") ** 2, "W": 2 * sm.param("g") ** 2, "B": sm.param("gp") ** 2 * 10 / 3}}
    yt = 3 * sm.param("y_u") ** 2 + 3 * sm.param("y_d") ** 2 + sm.param("y_e") ** 2
    parts = breaking_by_factor(brst)
    total = Density(sm.fields)
    for r in "GW":
        total = total + parts[r].scale(-inv["C_ad"][r] / 2)
    for r in "GWB":
        total = total + LocalFunctional.integrate(brst.b_d(golden.E_A(r))).scale(-inv["S"][r] / 3)
    total = total + LocalFunctional.integrate(brst.b_d(golden.E_H())).scale(-yt * 2 / 3)
    data = FunctionalData(brst, counterterm=C, counterterm_variation=LocalFunctional.integrate(total))
    res = pole_identity(ctx, singular, data)
    assert res.holds and len(res.classes) == 168


# ---------------------------------------------------------------------------------------------- V7: X and DeltaS
def test_v7_breaking_closed_not_anomalous(ctx, brst, X):
    res = match_breaking(ctx, FunctionalData(brst, breaking=X))
    assert res.closed and not res.anomaly.anomalous
    assert res.anomaly.image_rank == 84 and len(res.anomaly.complement) == len(ctx.N) - 84


def test_v7_restoring_generated_basis(sm, restoring, DeltaS):
    res = restoring
    assert res.unique and (res.unknowns, res.allowed, res.rank_allowed, res.rank_all) == (109, 70, 70, 84)
    assert res.existence_total == len(res.existence) == 14
    assert not any(c.relaxed for c in res.conditions)
    assert [c.rows for c in res.conditions] == [14, 18, 14, 18]
    ours = res.functional()
    assert equal_external(sm, ours, DeltaS, (4,))
    # the d-dimensional form is purely barred (report sec:Sren1): equal to its projection P in D = 6
    barred = LocalFunctional.integrate(ours.physical(sm.simplifier()))
    assert equal_external(sm, ours, barred, (6,))
    assert len(res.nonzero()) == 39


def test_v7_report_basis_by_name(sm, ctx, brst, X, q, golden_counterterms):
    """(N1)-(N4) as the report's removal lists on the report's 97 monomials: the 39 r_alpha exactly."""
    basis = report_basis(sm, q)
    norm = Normalisation([ElementCondition(k, k, v) for k, v in REMOVED.items()])
    res = solve_restoring(ctx, FunctionalData(brst, breaking=X), norm, basis=basis)
    assert (res.unknowns, res.allowed, res.rank_allowed, res.rank_all) == (97, 64, 64, 72) and res.unique
    assert res.existence_total == len(res.existence) == 8
    names = [n for n, _ in res.existence]
    sources = {"rWWB", "rGGB", "rBBB", "zGGB", "zWWB"}
    assert sources <= set(names)
    yukawa = [(n, cls) for n, cls in res.existence if n not in sources]
    assert len(yukawa) == 3 and all(n.startswith(("R_", "sig_")) for n, _ in yukawa)
    for _, cls in yukawa:      # testable only in the Yukawa classes gh^B psibar psi H
        assert cls and all(len(c) == 4 and "gh_B" in c and {"H", "Hdag"} & set(c) for c in cls)
    golden = golden_counterterms["r"]["all_97_engine_names"]
    nonzero = 0
    for name, value in golden.items():
        got = expr(ctx, res.coefficients[name])
        assert got == parse(sm, value), name
        nonzero += bool(res.coefficients[name])
    assert nonzero == 39


def test_v7_generic_conditions_on_report_basis(sm, ctx, brst, X, q, golden_counterterms):
    """The model-independent statement of (N1)-(N4) gives the report's removed space on its basis."""
    res = solve_restoring(ctx, FunctionalData(brst, breaking=X), standard_normalisation(sm),
                          basis=report_basis(sm, q))
    assert (res.allowed, res.rank_allowed, res.existence_total) == (64, 64, 8)
    for name, value in golden_counterterms["r"]["all_97_engine_names"].items():
        assert expr(ctx, res.coefficients[name]) == parse(sm, value), name


def test_v7_fallback_n4prime(sm, ctx, brst, X, q, restoring):
    """A breaking with a component along i y_k g' gh^B (O_k - O_k^dag): (N4) admits no solution, its fallback is
    applied automatically (pivots plus every compatible equation) and reported; the result differs from the
    shifted input by closed directions."""
    shift = LocalFunctional.integrate(q["R_UB+"].scale(sm.param("y_u") ** 2 / 3))
    Xp = LocalFunctional.integrate(X - brst.b(shift))
    res = solve_restoring(ctx, FunctionalData(brst, breaking=Xp))
    n4 = {c.name: c for c in res.conditions}["N4"]
    assert n4.relaxed and n4.dropped and res.unique
    diff = LocalFunctional.integrate(res.functional() - restoring.functional() - shift)
    assert Fingerprinter(sm).vanishes(LocalFunctional.integrate(brst.b(diff)), 3)


def test_v7_kinetic_and_literature_representatives(sm, ctx, brst, X, q, restoring, golden, golden_counterterms):
    kin = solve_restoring(ctx, FunctionalData(brst, breaking=X), standard_normalisation(sm, fermion_kinetic=True))
    assert kin.unique and {c.name for c in kin.conditions} >= {"N1", "N1k"}
    kinetic = [e.name for e in ctx.Q.elements if len(e.leg_class) == 2 and
               {sm.fields[f].kind for f in e.leg_class} == {"fermion", "fermion_bar"}]
    assert any(kin.coefficients[n] for n in kinetic)
    fp = Fingerprinter(sm)
    for other in (kin.functional(),
                  golden.restoring(golden_counterterms, q, ("r_literature_representative",
                                                                "BIMS_2004.14398_eq6.46_equals_LP_Stilde",
                                                                "r_engine_names"))):
        diff = LocalFunctional.integrate(other - restoring.functional())
        assert fp.vanishes(LocalFunctional.integrate(brst.b(diff)), 3)


def test_v7_anomaly_without_leptons(sm, golden_counterterms):
    """With the leptons removed the hypercharge anomaly does not cancel: the anomaly numbers (from the model's
    anomaly tensor) are D_BBB = -3/2, D_BWW = 1/2, D_BGG = D_GGG = 0, and the solver reports the anomaly instead
    of a solution."""
    import tomllib
    from importlib import resources

    doc = tomllib.loads(resources.files("hepkit_renorm.models").joinpath("sm_one_family.toml").read_text())
    for f in ("L", "E"):
        doc["fermions"].pop(f)
    doc["yukawa"] = [y for y in doc["yukawa"] if y["coupling"] != "y_e"]
    doc["parameters"].pop("y_e")
    doc["name"] = "sm_no_leptons"
    model = model_from_dict(doc)
    field = default_field()
    D = model.anomaly_tensor(field)
    offsets, k = {}, 0
    for name, fac in model.factors.items():
        offsets[name] = k
        k += fac.dim_adjoint
    b, w, g = offsets["B"], offsets["W"], offsets["G"]
    numbers = {"BBB": D[b, b, b], "BWW": D[b, w, w], "BGG": D[b, g, g]}
    assert {k: field.reconstruct(int(v)) for k, v in numbers.items()} == {"BBB": Fraction(-3, 2),
                                                                           "BWW": Fraction(1, 2), "BGG": 0}
    brst = BRST(model)
    gold = Golden(model)
    Xa = anomaly_breaking(model, gold, {"BBB": Fraction(-3, 2), "BWW": Fraction(1, 2), "BGG": 0, "GGG": 0})
    ctx = OrderContext(model, brst=brst)
    with pytest.raises(AnomalyError) as err:
        solve_restoring(ctx, FunctionalData(brst, breaking=Xa))
    assert err.value.report.anomalous
    # the full Standard Model content: the same formula with its anomaly numbers (all zero) is no anomaly
    Dsm = sm.anomaly_tensor(field)
    assert not np.any(Dsm % field.p)


# ---------------------------------------------------------------------------------------------- V8
@pytest.fixture(scope="module")
def store(sm, brst, ctx, singular, restoring):
    return update_store(tree_store(sm, brst), OrderResult(1, singular, restoring))


def test_v8_singular_rules(sm, brst, store):
    rules = CounterRules(sm, brst.action).singular()
    lf = loop_factor(1)
    ours = LocalFunctional.integrate(store.entry("deltaS1").coefficient(-1).scale(1 / lf))
    classes = {c for c in ours.leg_classes() if not any(sm.fields[n].kind == "auxiliary" for n in c)}
    assert classes <= set(rules), sorted(classes - set(rules))
    for cls, (expected, physical) in rules.items():
        part = LocalFunctional.integrate(Density.from_monomials(
            sm.fields, [m for m in ours.monomials() if sorted(m.fields()) == sorted(cls)]))
        assert equal_external(sm, part, expected, (6, 8), physical_classes=(cls,) if physical else ()), cls


def test_v8_finite_rules(sm, brst, store):
    lf = loop_factor(1)
    ours = LocalFunctional.integrate(store.entry("DeltaS1").coefficient(0).scale(1 / lf))
    expected = CounterRules(sm, brst.action).finite()
    barred = LocalFunctional.integrate(expected.physical(sm.simplifier()))   # all objects four-dimensional
    assert equal_external(sm, ours, barred, (4, 6))


def test_v8_store_and_model_file(sm, store):
    names = [e.name for e in store]
    assert names == ["S0", "Delta1", "deltaS1", "DeltaS1", "bd_deltaS1", "bd_DeltaS1"]
    assert [e.bookkeeping for e in store.entries(roles=["breaking"])] == ["", "A", "A"]
    back = OrderStore.from_json(sm, store.to_json())
    assert back.entry("DeltaS1").coefficient(0) == store.entry("DeltaS1").coefficient(0)
    written = counterterm_model(sm, store)
    hk_model = written.load()
    ops = {}
    for v in written.data["vertex_rules"]:
        ops.setdefault(v["operator"], []).append(v)
    assert set(ops) == {"S0", "Delta1", "deltaS1", "DeltaS1", "bd_deltaS1", "bd_DeltaS1"}
    for v in ops["DeltaS1"] + ops["deltaS1"]:
        orders = dict(written.reverse_map.vertex(v["name"]).orders)
        assert orders["CT"] == 1 and ("CTPOLE" in orders) == (v["operator"] == "deltaS1")
    assert len(hk_model.vertex_rules) == len(written.data["vertex_rules"])


def test_v8_written_rules_equal_entries(sm, store, field):
    """Every written counterterm vertex, read back from its UFO strings (the independent converter of
    test_output.py) and summed over the rules of the same particles (the writer splits a vertex by coupling
    orders), equals i V of the store entry at its components."""
    from component_values import ComponentEvaluator
    from symbolica import Expression
    from test_output import EPS, UfoConverter, by_name, json_rule_expression, same_values, ufo_params

    from hepkit_renorm import tensors as tn
    from hepkit_renorm.store import EVANESCENT_TAG

    written = counterterm_model(sm, store)
    ev = ComponentEvaluator(sm, field)
    d = written.data
    couplings, lorentz = by_name(d["couplings"]), by_name(d["lorentz_structures"])
    params = ufo_params(sm, written.reverse_map)
    groups: dict[tuple, list[dict]] = {}
    for v in d["vertex_rules"]:
        if v["operator"] in ("deltaS1", "DeltaS1"):
            groups.setdefault((v["operator"], tuple(v["particles"])), []).append(v)
    extra = [EPS.get_name(), EVANESCENT_TAG.get_name(), Expression.PI.get_name()]
    for key, rules in groups.items():
        first = rules[0]
        rec = written.reverse_map.vertex(first["name"])
        target = Expression.num(0)
        legs = None
        for k in rec.powers:
            rule = written.package_rule(first["name"], k)
            legs = rule.legs
            target = target + (tn.I * rule.value * EPS ** k).replace(tn.d, 4 - 2 * EPS)
        conv = UfoConverter(sm, legs, params)
        ours = Expression.num(0)
        for v in rules:
            assert written.reverse_map.vertex(v["name"]).components == rec.components
            ours = ours + json_rule_expression(sm, v, couplings, lorentz, conv)
        assert same_values(sm, ev, field, ours.expand(), target.expand(), legs, rec.components, extra), key
    assert len(groups) == 252
