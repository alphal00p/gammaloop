"""The order-l step and the order store it leaves for order l + 1 (hepkit_renorm.orders), on the toy model of
test_matching.py (fast); the Standard-Model store and model file are tested in test_sm_gates.py (slow)."""

import pytest
from test_matching import TOY, equal_external, q

from hepkit_renorm import Density, LocalFunctional, Occurrence, model_from_dict
from hepkit_renorm.bases import shift_gauge_sources
from hepkit_renorm.brst import BRST
from hepkit_renorm.matching import FunctionalData, OrderContext, field_directions
from hepkit_renorm.orders import counterterm_model, loop_factor, order_entries, prune, run_order, update_store
from hepkit_renorm.store import OrderStore, Role, StoreEntry, tree_store


@pytest.fixture(scope="module")
def toy():
    return model_from_dict(TOY)


@pytest.fixture(scope="module")
def brst(toy):
    return BRST(toy)


@pytest.fixture(scope="module")
def result(toy, brst):
    ctx = OrderContext(toy, brst=brst)
    gx, gy = toy.param("gx"), toy.param("gy")
    dirs = {d.name: d.functional for d in field_directions(toy, brst.action)}
    C = LocalFunctional.integrate(dirs["Z[X]"].scale(q(3, 2) * gx**2) + dirs["Z[gh_X]"].scale(gx**2)
                                  + ctx.singular.evanescent[0].functional.scale(-gx**2))
    vertex = next(e for e in ctx.Q.elements if set(e.leg_class) == {"P", "Pbar", "Y"})
    X = LocalFunctional.integrate(brst.b(vertex.functional).scale(-gx**2 * gy))
    data = FunctionalData(brst, counterterm=C, breaking=X)
    return run_order(ctx, data, data), vertex


def test_run_order(result):
    res, vertex = result
    assert res.order == 1 and res.poles.holds and res.breaking.closed and not res.breaking.anomaly.anomalous
    assert res.restoring.unique and set(res.restoring.nonzero()) == {vertex.name}


def test_store_entries_and_bookkeeping(toy, brst, result):
    res, _ = result
    entries = order_entries(res)
    by = {e.name: e for e in entries}
    assert list(by) == ["deltaS1", "DeltaS1", "bd_deltaS1", "bd_DeltaS1"]
    assert by["deltaS1"].role is Role.CT_SINGULAR and by["deltaS1"].powers == (-1,)
    assert by["DeltaS1"].role is Role.CT_FINITE and by["DeltaS1"].powers == (0,)
    assert by["bd_deltaS1"].bookkeeping == by["bd_DeltaS1"].bookkeeping == "A"
    store = update_store(tree_store(toy, brst), res)
    assert [e.name for e in store.entries(bookkeeping="B")] == ["S0", "Delta1", "deltaS1", "DeltaS1"]
    assert len(store.entries(bookkeeping="A")) == 6 and len(store.entries(max_loop_order=0)) == 2
    back = OrderStore.from_json(toy, store.to_json())
    assert [e.bookkeeping for e in back] == [e.bookkeeping for e in store]
    # the loop factor is written out: deltaS1 = (1/(16 pi^2)) C / eps
    lf = loop_factor(1)
    c = LocalFunctional.integrate(by["deltaS1"].coefficient(-1).scale(1 / lf))
    assert equal_external(toy, c, res.singular.functional(), (6,))
    # entries of one bookkeeping only
    entry = StoreEntry("x", 1, Role.BREAKING, {0: by["bd_DeltaS1"].coefficient(0)}, "", "B")
    assert entry.bookkeeping == "B"
    with pytest.raises(ValueError):
        StoreEntry("x", 1, Role.BREAKING, {0: by["bd_DeltaS1"].coefficient(0)}, "", "C")


def test_prune_removes_cancelling_classes(toy, brst):
    """b_d of the shifted gauge source has auxiliary-field terms that cancel (-d b + d b); prune drops the class."""
    labels = toy.slot_labels("gh_X")
    src = toy.fields["X"].source
    mu = toy.slot_labels(src)[0]
    rho = LocalFunctional.integrate(Density.monomial(toy.fields, 1, Occurrence(src, (mu,) + labels),
                                                     Occurrence("gh_X", labels, (mu,))))
    image = LocalFunctional.integrate(brst.b_d(shift_gauge_sources(toy, rho)))
    aux = [c for c in image.leg_classes() if "b_X" in c]
    assert aux
    assert not [c for c in prune(toy, image).leg_classes() if "b_X" in c]


def test_counterterm_model_loads(toy, brst, result):
    res, _ = result
    written = counterterm_model(toy, update_store(tree_store(toy, brst), res))
    ops = {v["operator"] for v in written.data["vertex_rules"]}
    assert {"S0", "Delta1", "deltaS1", "DeltaS1", "bd_deltaS1", "bd_DeltaS1"} <= ops
    assert len(written.load().vertex_rules) == len(written.data["vertex_rules"])
    # squared momenta are written with two dummies and a metric, never as a repeated factor P(d,k)*P(d,k)
    import re

    for lz in written.data["lorentz_structures"]:
        factors = re.findall(r"P\(dummy\(\d+\),\d+\)", lz["structure"])
        assert len(factors) == len(set(factors)), lz["structure"]
