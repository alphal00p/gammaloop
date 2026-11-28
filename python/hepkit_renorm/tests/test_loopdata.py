"""Loop data for the algebraic steps (hepkit_renorm.loopdata) on small models (fast): the derived parameter supports
of graphs, the resolution into parameter monomials against direct evaluation, the batched job evaluation, the dry
run of a consumer, and the local form of the poles of chiral QED (its renormalisation constants are known in closed
form: delta Z_A = -4/3 e^2, delta Z_psi = -e^2 in Feynman gauge, delta e = 2/3 e^3, units 1/(16 pi^2 eps))."""

from __future__ import annotations

import numpy as np
import pytest

from hepkit_renorm import model_from_dict
from hepkit_renorm.brst import BRST
from hepkit_renorm.class_results import ClassComputer
from hepkit_renorm.exact import default_field, second_field
from hepkit_renorm.loopdata import Request, ResolvedLoopData, localise, recorded
from hepkit_renorm.loops import evaluate_jobs
from hepkit_renorm.matching import POLE, FunctionalData, MatchingError, OrderContext, match_singular, pole_identity
from hepkit_renorm.store import tree_store

pytestmark = pytest.mark.filterwarnings("ignore:This process .* is multi-threaded:DeprecationWarning")

QED = {"name": "chiral_qed", "gauge": {"X": {"group": "U(1)", "coupling": "e"}},
       "fermions": {"P": {"chirality": "L", "X": "1"}, "N": {"chirality": "R", "X": "1"}},
       "gauge_fixing": {"kind": "linear_covariant", "xi": "xi"}}
"""A Dirac fermion of charge 1 written as two chiral multiplets (anomaly free)."""

YUKAWA = {"name": "yukawa_u1", "gauge": {"X": {"group": "U(1)", "coupling": "e"}},
          "fermions": {"P": {"chirality": "L", "X": "1"}, "N": {"chirality": "R", "X": "0"}},
          "scalars": {"phi": {"complex": True, "X": "1"}},
          "parameters": {"m2": {}, "lam": {}, "y": {}},
          "potential": [{"coupling": "m2", "bilinears": [["phi", "phi"]]},
                        {"coupling": "lam", "bilinears": [["phi", "phi"], ["phi", "phi"]]}],
          "yukawa": [{"coupling": "y", "fermion_bar": "P", "scalar": "phi", "fermion": "N"}],
          "gauge_fixing": {"kind": "linear_covariant", "xi": "xi"}}
"""A massive complex scalar with a Yukawa coupling (supports with mass monomials)."""


def setup(doc, processes=2):
    model = model_from_dict(doc)
    brst = BRST(model)
    ctx = OrderContext(model, brst=brst)
    data = ResolvedLoopData(ClassComputer(model, tree_store(model, brst)), ctx.space, processes=processes)
    return model, brst, ctx, data


@pytest.fixture(scope="module")
def qed():
    return setup(QED)


@pytest.fixture(scope="module")
def yukawa():
    return setup(YUKAWA)


def test_graph_supports_are_derived(yukawa):
    model, _, ctx, data = yukawa
    names = ctx.space.names
    _, _, _, sups = data.graphs(tuple(sorted(("phi", model.fields["phi"].conjugate), key=model.fields.order)))
    # every graph of the scalar two-point function: one coupling monomial, times m2^0, m2^1 where a propagator is
    # massive (superficial degree 2, m2 of mass dimension 2)
    assert sups and all(1 <= len(s) <= 2 for s in sups)
    assert any(len(s) == 2 for s in sups)
    m2 = names.index("m2")
    for s in sups:
        if len(s) == 2:
            assert {e[m2] for e in s} == {0, 1}


def test_resolved_data_equal_direct_evaluation(yukawa):
    model, _, ctx, data = yukawa
    phi = model.fields["phi"]
    cls = tuple(sorted(("phi", phi.conjugate), key=model.fields.order))
    fp = ctx.fingerprinter(default_field(), 6)
    smp = fp.samples(list(cls), 3, salt=5)
    poly = data.parametric(POLE, cls, fp, ctx.space, smp)
    point = ctx.space.random_points(fp.field, 1, seed=9)[0]
    direct = data.values(POLE, cls, fp.replace(values=ctx.space.values(point)), smp)
    assert np.array_equal(poly.at(point)[:, 0] % fp.field.p, direct % fp.field.p)
    assert np.any(direct)


def test_jobs_in_parallel_equal_serial(qed):
    _, _, ctx, data = qed
    fp = ctx.fingerprinter(default_field(), 4)
    jobs = []
    for cls in [("X", "X"), ("X", "P", "Pbar")]:
        _, _, prepared, _ = data.graphs(cls)
        smp = fp.samples(list(cls), 2, salt=3)
        jobs += [(pg, fp, smp) for pg in prepared]
    serial = evaluate_jobs(jobs, processes=1)
    parallel = evaluate_jobs(jobs, processes=3)
    for a, b in zip(serial, parallel):
        assert set(a.coefficients) == set(b.coefficients)
        assert all(np.array_equal(a.coefficients[k], b.coefficients[k]) for k in a.coefficients)


def test_localised_poles_of_chiral_qed(qed):
    model, brst, ctx, data = qed
    form = localise(ctx, data, POLE, ctx.singular.elements, ctx.classes(0), dimension=6)
    assert form.primes == [default_field().p] and form.confirmed == second_field().p
    from hepkit_renorm import LocalFunctional

    C = LocalFunctional.integrate(form.functional(ctx).scale(-1))
    res = match_singular(ctx, FunctionalData(brst, counterterm=C))
    e = model.param("e")
    expr = lambda k: ctx.space.expression(res.constants[k])
    assert expr("Z[X]") == (-4 * e**2 / 3).expand()
    assert expr("Z[P]") == expr("Z[N]") == (-(e**2)).expand()
    assert expr("d[e]") == (2 * e**3 / 3).expand()
    assert len(res.evanescent) == 2
    # the vacuum polarisation is not in the span of the fermion elements alone
    fermions = [el for el in ctx.singular.elements if "X" not in el.leg_class]
    with pytest.raises(MatchingError) as err:
        localise(ctx, data, POLE, fermions, [("X", "X")], dimension=6)
    assert ("X", "X") in err.value.classes


def test_dry_run_records_the_requests_of_a_consumer(qed):
    _, brst, ctx, data = qed
    form = localise(ctx, data, POLE, ctx.singular.elements, ctx.classes(0), dimension=6)
    from hepkit_renorm import LocalFunctional

    C = LocalFunctional.integrate(form.functional(ctx).scale(-1))
    singular = match_singular(ctx, FunctionalData(brst, counterterm=C))
    classes = [("gh_X", "P", "Pbar")]
    reqs = recorded(lambda d: pole_identity(ctx, singular, d, classes=classes), ctx.space)
    assert [r.fields for r in reqs] == classes and all(isinstance(r, Request) and r.resolved for r in reqs)
    before = data.stats.batches
    data.prefetch(reqs)
    assert data.stats.batches == before + 1
    res = pole_identity(ctx, singular, data, classes=classes)
    assert res.holds and data.stats.batches == before + 1        # served from the batch
