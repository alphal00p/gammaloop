"""The diagram pipeline's results through the interface of the algebraic modules (slow: a few minutes; run with
``pytest --runslow tests/test_sm_pipeline.py``).

Loop data (``LoopData``: :class:`~hepkit_renorm.class_results.ClassComputer` with the one-loop back-end, from the
tree-level store) are values at one parameter point. At the rational point of ``tests/data/class_results_golden.json``
every parameter is fixed, and the matching and the restoring solve give exact numbers there. A selection of classes
(gauge, fermion and Higgs two-point functions at ghost number 0; ghost-gauge, ghost-fermion and ghost-Higgs classes
at ghost number 1) is served by the pipeline, the others by the golden functionals (``CombinedData``): the
renormalisation constants, the pole identity and the restoring counterterm agree with the report at that point. The
complete one-loop run (every class from the pipeline) is the end-to-end step of the next milestone.
"""

from __future__ import annotations

from fractions import Fraction

import pytest
from report_golden import Golden, parse
from symbolica import Expression

from hepkit_renorm.brst import BRST
from hepkit_renorm.class_results import ClassComputer
from hepkit_renorm.matching import CombinedData, FunctionalData, LoopData, OrderContext, match_singular, pole_identity
from hepkit_renorm.st_solve import solve_restoring
from hepkit_renorm.store import tree_store

pytestmark = pytest.mark.slow

POINT = {"g_s": Fraction(2), "g": Fraction(3, 2), "gp": Fraction(1, 2), "y_u": Fraction(5, 4), "y_d": Fraction(3, 4),
         "y_e": Fraction(2), "lambda": Fraction(3, 2), "m2": Fraction(1, 4), "xi": Fraction(1)}
POLES = [("G", "G"), ("W", "W"), ("B", "B"), ("Q", "Qbar"), ("H", "Hdag")]
GHOST = [("G", "gh_G"), ("W", "gh_W"), ("B", "gh_B"), ("gh_G", "Q", "Qbar"), ("gh_W", "L", "Lbar")]


def at_point(sm, e: Expression) -> Fraction:
    for name, v in POINT.items():
        e = e.replace(sm.param(name), Expression.num(v.numerator) / v.denominator)
    e = e.expand()
    return Fraction(e.to_canonical_string())


@pytest.fixture(scope="module")
def setup(sm, golden_counterterms):
    brst = BRST(sm)
    gold = Golden(sm)
    q = gold.q_basis()
    ctx = OrderContext(sm, brst=brst, fixed=POINT)
    loop = LoopData(ClassComputer(sm, tree_store(sm, brst)))
    C, X = gold.counterterm(golden_counterterms, q), gold.breaking(golden_counterterms)
    return brst, ctx, loop, C, X


def test_pipeline_poles(sm, setup, golden_counterterms):
    brst, ctx, loop, C, _ = setup
    data = CombinedData(FunctionalData(brst, counterterm=C), {c: loop for c in POLES + GHOST})
    res = match_singular(ctx, data)
    z = golden_counterterms["delta_S"]["Z_factors"]["values_in_kappa_units"]
    for r in "GWB":
        assert res.constants[f"Z[{r}]"] == {(): at_point(sm, parse(sm, z["dZ_A"][r]))}
    assert res.constants["Z[Q]"] == {(): at_point(sm, parse(sm, z["dZ_f"]["Q"]))}
    assert res.constants["Z[H]"] == {(): at_point(sm, parse(sm, z["dZ_H"]))}
    assert {k[:2] for k in loop.results} >= {("pole", c) for c in POLES}
    ident = pole_identity(ctx, res, data, classes=GHOST)
    assert ident.holds


def test_pipeline_breaking(sm, setup, golden_counterterms):
    brst, ctx, loop, _, X = setup
    data = CombinedData(FunctionalData(brst, breaking=X), {c: loop for c in GHOST})
    res = solve_restoring(ctx, data)
    assert res.unique
    gold = Golden(sm)
    expected = gold.restoring(golden_counterterms, gold.q_basis())
    ours = res.functional()
    from test_matching import equal_external

    pointed = ours
    for name, v in POINT.items():
        sym, val = sm.param(name), Expression.num(v.numerator) / v.denominator
        expected = expected.map_coefficients(lambda c, s=sym, x=val: c.replace(s, x).expand())
        pointed = pointed.map_coefficients(lambda c, s=sym, x=val: c.replace(s, x).expand())
    assert equal_external(sm, pointed, expected, (4,))


def test_pipeline_bookkeeping_b_remainder(sm, setup, golden_counterterms):
    """eq:m-LP: the evanescent part of the eta/eps coefficients of the ghost-number-1 classes, matched onto the
    hatted basis of those classes; it vanishes in four dimensions, and N1 = V_X + V_remainder in D = 6 on fresh
    samples. It enters the store as hatDelta1 (bookkeeping B)."""
    import numpy as np

    from hepkit_renorm.matching import match_remainder
    from hepkit_renorm.orders import OrderResult, order_entries

    brst, ctx, loop, C, X = setup
    gold = Golden(sm)
    DeltaS = gold.restoring(golden_counterterms, gold.q_basis())
    classes = [("G", "gh_G"), ("gh_G", "Q", "Qbar")]
    rem = match_remainder(ctx, loop, DeltaS, classes)
    R = rem.functional()
    assert not R.is_zero()
    from hepkit_renorm.exact import default_field
    from hepkit_renorm.functional import LocalFunctional

    fp4 = ctx.fingerprinter(default_field(), 4)
    for c in classes:
        assert not np.any(fp4.functional(R, fp4.samples(list(c), 3, salt=5)))
    fp6 = ctx.fingerprinter(default_field(), 6)
    Xf = LocalFunctional.integrate(brst.b(DeltaS).scale(-1))
    for c in classes:
        smp = fp6.samples(list(c), 3, salt=99)
        n1 = loop.values("breaking", c, fp6, smp)
        assert np.array_equal(n1 % fp6.field.p, (fp6.functional(Xf, smp) + fp6.functional(R, smp)) % fp6.field.p)
    singular = match_singular(ctx, FunctionalData(brst, counterterm=C))
    restoring = solve_restoring(ctx, FunctionalData(brst, breaking=X))
    entries = order_entries(OrderResult(1, singular, restoring, remainder=R))
    hat = [e for e in entries if e.name == "hatDelta1"]
    assert hat and hat[0].bookkeeping == "B"
