"""Tests of the fixes made in the review of MS2 to MS7 (fast).

* ``diagrams.store_families`` inserts the breaking entries of one bookkeeping only (report sec:m-beyond,
  tab:m-bookkeeping): before, an order-2 store with both ``bd_*`` (A) and ``hatDelta`` (B) entries gave both.
* ``st_solve`` caches the ``b`` images per basis object, never by a bare ``id`` that another basis can reuse.
* ``matching.feynman_gauge`` is the one derivation of the Feynman-gauge substitution (driver, orders, contexts).
"""

from __future__ import annotations

import pytest
from test_loopdata import QED

from hepkit_renorm import model_from_dict
from hepkit_renorm.brst import BRST
from hepkit_renorm.class_results import ClassComputer
from hepkit_renorm.diagrams import store_families
from hepkit_renorm.functional import LocalFunctional
from hepkit_renorm.matching import OrderContext, feynman_gauge
from hepkit_renorm.st_solve import _b_images
from hepkit_renorm.store import Role, StoreEntry, tree_store


@pytest.fixture(scope="module")
def qed():
    model = model_from_dict(QED)
    return model, BRST(model)


def _store_with_both_bookkeepings(model, brst):
    store = tree_store(model, brst)
    delta = brst.breaking
    store.add(StoreEntry("bd_x1", 1, Role.BREAKING, {0: delta}, "bookkeeping A insertion", "A"))
    store.add(StoreEntry("hat_x1", 1, Role.BREAKING, {0: delta}, "bookkeeping B insertion", "B"))
    return store


def test_store_families_respect_the_bookkeeping(qed):
    model, brst = qed
    store = _store_with_both_bookkeepings(model, brst)
    names = {b: {f.name for f in store_families(model, store, 2, b)} for b in ("A", "B")}
    assert names["A"] == {"S0", "Delta1", "bd_x1"}
    assert names["B"] == {"S0", "Delta1", "hat_x1"}
    assert {f.name for f in store_families(model, store, 1)} == {"S0", "Delta1"}      # default A, order 1
    with pytest.raises(ValueError, match="bookkeeping"):
        store_families(model, store, 2, "C")
    cc = ClassComputer(model, store, 2, bookkeeping="B")
    assert {f.name for f in cc.families} == names["B"]


def test_b_images_are_cached_per_basis_object(qed):
    model, brst = qed
    ctx = OrderContext(model, brst=brst)
    q = ctx.Q
    first = _b_images(ctx, q)
    assert _b_images(ctx, q) is first
    # a stale entry under the same id (another basis object) is rebuilt, not returned
    ctx.replace_cached(("b", id(q)), (object(), ["stale"]))
    again = _b_images(ctx, q)
    assert again != ["stale"] and len(again) == len(q)
    assert all(isinstance(f, LocalFunctional) for f in again)


def test_feynman_gauge_substitution(qed):
    model, _ = qed
    assert feynman_gauge(model) == {"xi": 1}
    ctx = OrderContext(model)
    assert dict(ctx.space.fixed) == {"xi": 1}
