"""The restoring condition with normalisation conditions (hepkit_renorm.st_solve) on the toy SU(2) x U(1) model
(fast; the Standard-Model gate V7 is in test_sm_gates.py, slow)."""

import numpy as np
import pytest
from test_matching import TOY, q

from hepkit_renorm import LocalFunctional, model_from_dict
from hepkit_renorm.brst import BRST
from hepkit_renorm.fingerprint import Fingerprinter
from hepkit_renorm.matching import FunctionalData, OrderContext
from hepkit_renorm.st_solve import (
    ClassCondition,
    ElementCondition,
    ExistenceError,
    Normalisation,
    solve_restoring,
    standard_normalisation,
)


@pytest.fixture(scope="module")
def toy():
    return model_from_dict(TOY)


@pytest.fixture(scope="module")
def brst(toy):
    return BRST(toy)


@pytest.fixture(scope="module")
def ctx(toy, brst):
    return OrderContext(toy, brst=brst)


def _element(ctx, fields):
    return next(e for e in ctx.Q.elements if set(e.leg_class) == set(fields) and len(e.leg_class) == len(fields))


@pytest.fixture(scope="module")
def target(ctx):
    """A restoring counterterm in the allowed space: the hypercharge fermion vertex and a longitudinal two-point
    term of the SU(2) gauge field, with polynomial coefficients."""
    gx, gy = ctx.model.param("gx"), ctx.model.param("gy")
    vertex = _element(ctx, ("Pbar", "P", "Y"))
    two = [e for e in ctx.Q.elements if e.leg_class == ("X", "X") and e.degree == 2]
    # the longitudinal combination: the one whose vertex is p^mu p^nu (fixed by the normalisation N2 itself)
    return vertex, two, (q(5, 4) * gx**2 * gy, q(-2, 3) * gx**2)


def _solve(ctx, brst, functional, normalisation=None):
    X = LocalFunctional.integrate(brst.b(functional).scale(-1))
    return solve_restoring(ctx, FunctionalData(brst, breaking=X), normalisation)


def test_standard_normalisation_unique(ctx, brst, target):
    vertex, _, (cv, _) = target
    f = LocalFunctional.integrate(vertex.functional.scale(cv))
    res = _solve(ctx, brst, f)
    assert res.unique and res.kernel == 0 and res.allowed == res.rank_allowed
    assert res.unknowns == len(ctx.Q)
    assert {c.name for c in res.conditions} == {"N1", "N2", "N3", "N4"}
    assert not any(c.relaxed for c in res.conditions)
    assert ctx.space.expression(res.coefficients[vertex.name]) == cv.expand()
    assert set(res.nonzero()) == {vertex.name}
    # existence conditions: rank(b on Q) - rank(b on the allowed space), located with their test classes
    assert res.existence_total == res.rank_all - res.rank_allowed == len(res.existence)
    assert all(classes for _, classes in res.existence)


def test_closed_directions_are_removed(ctx, brst, target):
    """Adding a closed direction (the X field renormalisation) to the breaking's source does not change the
    solution: the normalisation removes it."""
    from hepkit_renorm.matching import field_directions

    vertex, _, (cv, _) = target
    f = LocalFunctional.integrate(vertex.functional.scale(cv))
    closed = next(d for d in field_directions(ctx.model, brst.action) if d.name == "Z[X]").functional
    a = _solve(ctx, brst, f)
    b = _solve(ctx, brst, LocalFunctional.integrate(f + closed.scale(ctx.model.param("gx") ** 2)))
    assert a.coefficients == b.coefficients


def test_minimal_mode_and_fallback(ctx, brst, target):
    vertex, _, (cv, _) = target
    f = LocalFunctional.integrate(vertex.functional.scale(cv))
    minimal = _solve(ctx, brst, f, standard_normalisation(ctx.model, mode="minimal"))
    assert minimal.unique and minimal.allowed >= _solve(ctx, brst, f).allowed
    fp = Fingerprinter(ctx.model)
    diff = LocalFunctional.integrate(minimal.functional() - f)
    assert fp.vanishes(LocalFunctional.integrate(brst.b(diff)))      # same breaking: differ by closed directions
    # a condition that removes a needed monomial: the solution does not exist with all its equations; the fallback
    # keeps the pivots (none: the vertex is not a closed direction) and the equations compatible with X
    extra = ClassCondition("Nx", "no hypercharge fermion vertex", lambda m, c: set(c) == {"Pbar", "P", "Y"})
    norm = standard_normalisation(ctx.model)
    relaxed = _solve(ctx, brst, f, Normalisation(norm.conditions + [extra]))
    rep = {c.name: c for c in relaxed.conditions}["Nx"]
    assert rep.relaxed and rep.dropped and relaxed.unique
    assert ctx.space.expression(relaxed.coefficients[vertex.name]) == cv.expand()
    with pytest.raises(ExistenceError):
        _solve(ctx, brst, f, Normalisation(norm.conditions + [extra], fallback=False))


def test_element_condition_and_kinetic_representative(ctx, brst, target):
    vertex, _, (cv, _) = target
    f = LocalFunctional.integrate(vertex.functional.scale(cv))
    norm = standard_normalisation(ctx.model)
    named = Normalisation(norm.conditions + [ElementCondition("Nn", "by name", (vertex.name,))])
    res = _solve(ctx, brst, f, named)
    assert {c.name: c for c in res.conditions}["Nn"].relaxed
    kin = _solve(ctx, brst, f, standard_normalisation(ctx.model, fermion_kinetic=True))
    assert kin.unique
    names = {c.name for c in kin.conditions}
    assert "N1k" in names
    diff = LocalFunctional.integrate(kin.functional() - f)
    assert Fingerprinter(ctx.model).vanishes(LocalFunctional.integrate(brst.b(diff)))


def test_transverse_condition_fixes_the_gauge_kinetic_direction(ctx, brst, target):
    """N2: the two-point part of the solution is longitudinal; with a transverse input term (a closed direction
    together with its cubic and quartic partners) the solution is unchanged."""
    _, two, _ = target
    norm = standard_normalisation(ctx.model)
    fp = Fingerprinter(ctx.model)
    rows = norm.conditions[1].rows(ctx, ctx.Q, ctx.fingerprinter(fp.field, 4))
    m = rows.at(np.zeros(len(ctx.space), dtype=np.int64) + 3)
    k = [ctx.Q.index(e.name) for e in two]
    from hepkit_renorm import linalg

    assert linalg.rank(m[:, k], fp.field) == 1             # one transverse direction per factor in the degree-2 part
