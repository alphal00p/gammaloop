"""Operator bases and the action of b (modules ``bases`` and ``baction``): the fast part of gate V4.

The full one-family Standard Model gate (97 / 109 monomials Q, 353 N, rank 72 / 84, kernel 25, eq:closed, the
b-table of st_btable_unit.json, the named rows of eq:btable-gauge and eq:btable-higgs) is in
``test_bases_sm.py`` (marked slow, about 1.5 min; ``pytest --runslow tests/test_bases_sm.py``). Here: the generator
pieces, completeness and hermiticity on a toy model and on Standard-Model classes, the shifted gauge source, the
b-table / kernel / derived directions of the toy model, the evanescent operators of report eq:evanescent and
eq:b-evanescent and the d-dimensional statement eq:bd-Oc.
"""

from fractions import Fraction

import numpy as np
import pytest
from symbolica import Expression

from hepkit_renorm import Density, LocalFunctional, Occurrence, linalg, model_from_dict
from hepkit_renorm import tensors as tn
from hepkit_renorm.baction import (
    b_closure,
    b_table,
    bd_singular,
    breaking_by_factor,
    charge_deformation,
    coupling_direction,
    ghost_direction,
    global_symmetries,
    invariant_directions,
)
from hepkit_renorm.bases import (
    IncompleteBasisError,
    _group_tensors,
    candidates,
    continuation,
    counterterm_basis,
    external_classes,
    operator_basis,
    shift_gauge_sources,
    singular_basis,
)
from hepkit_renorm.brst import BRST
from hepkit_renorm.classes import invariant_trace, leg_classes
from hepkit_renorm.fingerprint import Fingerprinter, is_zero
from hepkit_renorm.groups import GroupEvaluator

TOY = {"name": "toy", "gauge": {"X": {"group": "SU(2)", "coupling": "gx"}, "Y": {"group": "U(1)", "coupling": "gy"}},
       "fermions": {"P": {"chirality": "L", "X": "F", "Y": "1/3"}, "N": {"chirality": "R", "Y": "1/3"}},
       "scalars": {"phi": {"complex": False, "X": "A"}},
       "parameters": {"mphi": {}},
       "potential": [{"coupling": "mphi", "bilinears": [["phi", "phi"]]}],
       "gauge_fixing": {"kind": "linear_covariant", "xi": "1"}}


@pytest.fixture(scope="module")
def toy():
    return model_from_dict(TOY)


@pytest.fixture(scope="module")
def toy_brst(toy):
    return BRST(toy)


@pytest.fixture(scope="module")
def toy_q(toy):
    return counterterm_basis(toy)


@pytest.fixture(scope="module")
def toy_n(toy):
    return operator_basis(toy, 1, hermitian=False)


@pytest.fixture(scope="module")
def toy_table(toy_q, toy_n, toy_brst):
    return b_table(toy_q, toy_n, brst=toy_brst)


def zero(model, f) -> bool:
    return is_zero(model, LocalFunctional.integrate(f))


# ---------------------------------------------------------------------------------------------- generator pieces
@pytest.mark.parametrize("factor,reps", [("G", "AAA"), ("G", "AAAA"), ("G", ["F", "Fbar", "A"]), ("G", ["F", "F", "F"]),
                                         ("W", ["F", "F"]), ("W", ["Fbar", "Fbar", "A"]), ("W", "AAAA"),
                                         ("G", ["F", "Fbar", "F", "Fbar"])])
def test_group_tensors_span_the_invariants(sm, field, factor, reps):
    """The trace / chain / epsilon candidates of one factor span the invariant tensors: their rank equals the
    dimension of the null space of the generators (character count of :func:`classes.invariant_trace`)."""
    fac = sm.factors[factor]
    labels = [tn.fresh(fac.adjoint_space() if r == "A" else fac.rep_space("F")) for r in reps]
    tensors = _group_tensors(fac, list(zip(labels, reps)))
    ev = GroupEvaluator(sm.factors, field)
    arrays = [ev.evaluate(t.expand(), labels).reshape(-1) for t in tensors]
    expected = invariant_trace(fac, [(r, 1) for r in reps])
    assert linalg.rank(np.stack(arrays, axis=1) if arrays else np.zeros((1, 0)), field) == expected
    assert expected > 0


def test_candidates_are_deterministic_and_four_dimensional(sm):
    a = [str(m.coefficient) for m in candidates(sm, ("G", "gh_G", "Q", "Qbar"), 1)][:5]
    b = [str(m.coefficient) for m in candidates(sm, ("G", "gh_G", "Q", "Qbar"), 1)][:5]
    assert len(a) == 5
    assert [len(x) for x in a] == [len(x) for x in b]           # same structure (labels are fresh)
    for m in candidates(sm, ("G", "G", "G"), 1):
        assert not tn.is_evanescent(m.coefficient)
        assert not tn.contains_head(m.coefficient, tn.g)


def test_ibp_field_has_no_derivative(sm):
    """Ghost legs of N stay undifferentiated (report tab:Nbasis); derivatives are distributed elsewhere."""
    for m in candidates(sm, ("G", "G", "gh_G"), 2):
        ghost = [o for o in m.occurrences if o.field == "gh_G"]
        assert all(not o.derivatives for o in ghost)


# ---------------------------------------------------------------------------------------------- completeness
def test_toy_bases_are_complete_and_independent(toy, toy_q, toy_n):
    """Every class reaches its exact structure count (classes.structure_counts): sum over the classes; the basis
    matrix has full column rank on its own layout."""
    for basis, g in ((toy_q, 0), (toy_n, 1)):
        expected = sum(c.n_structures for c in leg_classes(toy, g))
        assert len(basis) == expected
        fp = Fingerprinter(toy)
        assert linalg.rank(basis.matrix(fp).values, fp.field) == len(basis)
    assert len(toy_q) == 51 and len(toy_n) == 135


def test_toy_q_is_hermitian(toy, toy_q):
    for e in toy_q:
        assert zero(toy, e.functional - e.functional.conjugate()), e.name
    assert {e.kind for e in toy_q} == {"hermitian+", "hermitian-"}


def test_sm_classes_match_tab_qbasis_rows(sm):
    """Single sectors of tab:Qbasis: 9 gauge two-point, 12 three-point, 5 Q_K, 10 Q_V, 6 Yukawa monomials."""
    vv = [("G", "G"), ("W", "W"), ("B", "B")]
    assert len(operator_basis(sm, 0, classes=vv)) == 9
    vvv = [c.fields for c in leg_classes(sm, 0) if len(c.fields) == 3 and all(sm.fields[n].kind == "gauge"
                                                                               for n in c.fields)]
    assert len(operator_basis(sm, 0, classes=vvv)) == 12
    ff = [c.fields for c in leg_classes(sm, 0) if c.signature == "FF"]
    assert len(operator_basis(sm, 0, classes=ff)) == 5
    ffs = [c.fields for c in leg_classes(sm, 0) if c.signature == "FFS"]
    q = operator_basis(sm, 0, classes=ffs)
    assert len(q) == 6 and sorted(e.kind for e in q) == ["hermitian+"] * 3 + ["hermitian-"] * 3


def test_rejects_unknown_classes(sm):
    with pytest.raises(ValueError):
        operator_basis(sm, 0, classes=[("G", "W")])
    assert issubclass(IncompleteBasisError, RuntimeError)


def test_coordinates_and_span(sm):
    """A sum of basis elements has the expected coordinates; a functional outside the span is rejected."""
    q = operator_basis(sm, 0, classes=[("H", "Hdag"), ("G", "G")])
    fp = Fingerprinter(sm)
    f = LocalFunctional.integrate(q[0].functional.scale(3) + q[2].functional.scale(-2))
    c = q.coordinates([f], fp)[0]
    assert [int(x) for x in c] == [3, 0, fp.field(-2)] + [0] * (len(q) - 3)
    mu, a = tn.fresh(tn.LORENTZ), tn.fresh(sm.factors["W"].adjoint_space())
    outside = LocalFunctional.integrate(Density.monomial(sm.fields, 1, Occurrence("W", (mu, a)),
                                                         Occurrence("W", (mu, a))))
    assert not q.contains([outside], fp)


# ---------------------------------------------------------------------------------------------- shifted source
def test_gauge_sources_enter_as_rhotilde(sm):
    """report eq:rhotilde: the elements with a gauge source satisfy the ghost equation, and b of them has no
    auxiliary-field legs (the b-terms cancel)."""
    brst = BRST(sm)
    q = operator_basis(sm, 0, classes=[("G", "gh_G", "src_G"), ("gh_B", "src_B")])
    assert len(q) == 3                           # Q_rho[fg], Q_rho[g_s d] and Q_rho_partial^(B)
    for e in q:
        for r in ("G", "B"):
            assert brst.ghost_equation(e.functional, r).simplify(sm.simplifier()).is_zero(), e.name
        image = LocalFunctional.integrate(brst.b(e.functional))
        aux = Density.from_monomials(sm.fields, [m for m in image.monomials()
                                                 if any(n.startswith("b_") for n in m.fields())])
        assert zero(sm, aux), e.name
    lab = sm.slot_labels("H")
    plain = LocalFunctional.integrate(Density.monomial(sm.fields, 1, Occurrence("H", lab), Occurrence("Hdag", lab)))
    assert shift_gauge_sources(sm, plain) == plain


# ---------------------------------------------------------------------------------------------- b-table, kernel
def test_toy_b_table(toy, toy_q, toy_n, toy_table):
    """b Q = B N on fresh samples (another salt), rank + kernel = |Q|, the kernel is closed."""
    t = toy_table
    assert t.rank + t.kernel.shape[0] == len(toy_q)
    assert (t.rank, t.kernel.shape[0]) == (39, 12)
    fp = t.fingerprinter
    lhs = toy_n.matrix(fp, t.images, salt=5).values.astype(object)
    rhs = (toy_n.matrix(fp, salt=5).values.astype(object) @ t.matrix.T.astype(object)) % fp.field.p
    assert not np.any((lhs - rhs) % fp.field.p)
    for vec in t.kernel:
        assert t.is_closed(vec)


def test_toy_b_table_is_rational_at_unit_couplings(toy, toy_table):
    """Exact complex-rational table at unit couplings (two embeddings of i, several primes); consistent with the
    table over the prime field at the same values."""
    unit = {"gx": 1, "gy": 1, "mphi": 1}
    re, im = toy_table.rational(unit)
    flat = [x for x in list(re.flat) + list(im.flat)]
    assert all(isinstance(x, Fraction) for x in flat)
    assert max(x.denominator for x in flat) <= 6 and any(x for x in im.flat)
    fp = Fingerprinter(toy, values=unit)
    direct = b_table(toy_table.source, toy_table.target, fingerprinter=fp).matrix
    F = fp.field
    image = np.vectorize(lambda a, b: (F(a) + F(b) * F.i) % F.p, otypes=[object])(re, im)
    assert np.array_equal(image.astype(np.int64), direct)


def test_toy_kernel_is_spanned_by_derived_directions(toy, toy_table, toy_brst):
    """Field renormalisations (b-exact), coupling derivatives, charge deformations and gauge-invariant operators
    span the kernel; the quartic invariant (phi.phi)^2 is closed but not a coupling derivative of the toy action."""
    k = invariant_directions(toy_table, toy_brst)
    assert k.dimension == 12 and k.spanned == 12 and k.remainder.shape[0] == 0
    assert k.rank_by_kind == {"exact": 7, "coupling": 10, "phase": 10, "charge": 11, "invariant": 12}
    assert all(d.closed for d in k.directions if d.kind in ("exact", "coupling", "charge", "invariant"))


def test_global_symmetries_and_charge_deformation(sm):
    """The SM one-family global U(1)s are spanned by hypercharge, baryon and lepton number; the hypercharge
    deformation is g' dS_0/dg'."""
    syms = global_symmetries(sm)
    names = ["Q", "U", "D", "L", "E", "H"]
    vec = lambda d: [Fraction(d.get(n, 0)) for n in names]
    hyper = {n: sm.fields[n].charges["B"] for n in names}
    baryon = {"Q": Fraction(1, 3), "U": Fraction(1, 3), "D": Fraction(1, 3)}
    lepton = {"L": 1, "E": 1}
    from hepkit_renorm.exact import default_field
    F = default_field()
    m = [[F(x) for x in vec(s)] for s in syms]
    assert len(syms) == 3
    assert linalg.rank(np.array(m + [[F(x) for x in vec(d)] for d in (hyper, baryon, lepton)]), F) == 3
    brst = BRST(sm)
    deform = charge_deformation(sm, "B", hyper, brst)
    gp = sm.coupling("B")
    direct = LocalFunctional.integrate(brst.action.map_coefficients(lambda c: (c.derivative(gp) * gp).expand()))
    assert Fingerprinter(sm).equal(deform, direct)


# ---------------------------------------------------------------------------------------------- d dimensions
@pytest.fixture(scope="module")
def sm_singular(sm):
    q = operator_basis(sm, 0, classes=[("G", "G"), ("W", "W"), ("B", "B"), ("H", "Hdag"), ("Q", "Qbar")])
    return q, singular_basis(q)


def _report_evanescent(sm):
    """eq:evanescent and eq:Nhat (test data): E_A^(r), E_H and Nhat_A^(r) + Nhat_AA^(r)."""
    F, out = sm.fields, {}
    for r, fac in sm.factors.items():
        mu, nu, rho, sig = (tn.fresh(tn.LORENTZ) for _ in range(4))
        a = () if fac.abelian else (tn.fresh(fac.adjoint_space()),)
        ea = Density.monomial(F, Expression.num(1) / 2 * tn.gbar(mu, nu) * tn.ghat(rho, sig),
                              Occurrence(r, (mu,) + a), Occurrence(r, (nu,) + a, (rho, sig)))
        na = Density.monomial(F, tn.gbar(mu, nu) * tn.ghat(rho, sig), Occurrence(f"gh_{r}", a, (mu,)),
                              Occurrence(r, (nu,) + a, (rho, sig)))
        if not fac.abelian:
            b, c = tn.fresh(fac.adjoint_space()), tn.fresh(fac.adjoint_space())
            na = na + Density.monomial(F, sm.coupling(r) * tn.f(a[0], b, c) * tn.gbar(mu, nu) * tn.ghat(rho, sig),
                                       Occurrence(r, (mu, b)), Occurrence(f"gh_{r}", (c,)),
                                       Occurrence(r, (nu, a[0]), (rho, sig)))
        out[f"E_A_{r}"] = (LocalFunctional.integrate(ea), LocalFunctional.integrate(na))
    mu, nu = tn.fresh(tn.LORENTZ), tn.fresh(tn.LORENTZ)
    lab = sm.slot_labels("H")
    eh = Density.monomial(F, tn.ghat(mu, nu), Occurrence("Hdag", lab, (mu,)), Occurrence("H", lab, (nu,)))
    out["E_H"] = (LocalFunctional.integrate(eh), None)
    return out


def test_continuation(sm):
    q = operator_basis(sm, 0, classes=[("G", "G"), ("Q", "Qbar")])
    for e in q:
        c = continuation(e).functional
        spinor = any(n in ("Q", "Qbar") for n in e.leg_class)
        assert (c == e.functional) == spinor, e.name
        if not spinor:
            assert not tn.contains_head(c.expr, tn.gbar)


def test_evanescent_operators_of_eq_evanescent(sm, sm_singular):
    """report eq:evanescent: E_A^(r) (3), E_H, Etilde_F^(r) (3), Etilde_H: eight operators, the report's ones in
    their span; fermion bilinears give none (P_cbar ghat P_c = 0)."""
    _, sing = sm_singular
    ev = sing.evanescent
    assert len(ev) == 8
    assert sorted(e.origin for e in ev if e.origin.startswith("P(")) == [
        "P(DS2[H]) - DS2[H]", "P(F2[B]) - F2[B]", "P(F2[G]) - F2[G]", "P(F2[W]) - F2[W]"]
    assert not any("Q" in e.leg_class for e in ev)
    for e in ev:
        assert zero(sm, e.functional - e.functional.conjugate()), e.name
    fps = [Fingerprinter(sm, dimension=D) for D in (6, 8)]
    base = np.concatenate([fp.matrix([e.functional for e in ev], layout=sing.layout, salt=3).values for fp in fps])
    F = fps[0].field
    assert linalg.rank(base, F) == 8
    for name, (op, _) in _report_evanescent(sm).items():
        extra = np.concatenate([fp.matrix([op], layout=sing.layout, salt=3).values for fp in fps])
        assert linalg.rank(np.concatenate([base, extra], axis=1), F) == 8, name
        assert Fingerprinter(sm).vanishes(op), name         # evanescent: zero in four dimensions


def test_bd_on_the_singular_basis(sm, sm_singular):
    """report sec:evanescent: Etilde_F and Etilde_H are b_d-closed, E_A and E_H are not, b_d E_A = Nhat_A +
    Nhat_AA (eq:b-evanescent, eq:Nhat); the continued Q_m (H^dag H) is closed."""
    brst = BRST(sm)
    _, sing = sm_singular
    result = {e.name: closed for e, _, closed in bd_singular(sing, brst)}
    for e in sing.evanescent:
        assert result[e.name] == e.origin.startswith("P("), e.name
    closed_continued = [e.name for e in sing.continued if result[e.name]]
    assert closed_continued == [e.name for e in sing.continued if e.leg_class == ("H", "Hdag") and e.degree == 0]
    for name, (op, nhat) in _report_evanescent(sm).items():
        if nhat is not None:
            assert zero(sm, brst.b_d(op) - nhat), name


@pytest.mark.parametrize("factor", ["G", "W", "B"])
def test_bd_of_coupling_and_ghost_directions_is_the_breaking(sm, factor):
    """report eq:bd-Oc: b_d O_g^(r) = b_d O_c^(r) = Delta_1^(r), with Delta_1 = sum_r Delta_1^(r)."""
    brst = BRST(sm)
    parts = breaking_by_factor(brst)
    total = Density(sm.fields)
    for p in parts.values():
        total = total + p
    assert zero(sm, total - brst.breaking)
    assert not zero(sm, parts[factor])
    assert zero(sm, brst.b_d(coupling_direction(brst, factor)) - parts[factor])
    assert zero(sm, brst.b_d(ghost_direction(brst, factor)) - parts[factor])


def test_external_classes_drop_antighost_and_auxiliary(sm):
    brst = BRST(sm)
    cls = external_classes(sm, brst.action)
    assert cls and not any(n.startswith(("agh_", "b_")) for c in cls for n in c)


def test_hatted_fields_give_the_complete_hatted_basis(sm):
    """Beyond one loop (hatted_fields=True) hatted field components enter: more evanescent operators in the
    abelian two-point class than the one-loop set (E_A and Etilde_F)."""
    q = operator_basis(sm, 0, classes=[("B", "B")])
    one_loop = singular_basis(q)
    full = singular_basis(q, hatted_fields=True)
    assert len(one_loop.evanescent) == 2 + 2 + 1      # E_A^(B), Etilde_F^(G, W, B), Etilde_H
    assert len(full.evanescent) > len(one_loop.evanescent)


@pytest.mark.slow
def test_b_squared_vanishes_on_the_tables(toy, toy_n, toy_table, toy_brst):
    """b^2 = 0 in four dimensions as a matrix identity: (Q -> N) (N -> K) = 0 with the complete ghost-number-2
    basis K of the toy model (235 monomials, about 25 s)."""
    k = operator_basis(toy, 2, hermitian=False)
    second = b_table(toy_n, k, brst=toy_brst, fingerprinter=toy_table.fingerprinter)
    p = toy_table.field.p
    assert len(k) == sum(c.n_structures for c in leg_classes(toy, 2))
    assert not np.any((toy_table.matrix.astype(object) @ second.matrix.astype(object)) % p)
    assert second.rank > 0


def test_b_closure_annihilates_the_image_of_b(toy, toy_n, toy_table, toy_brst):
    """The condition b X = 0 without a ghost-number-2 basis: the fingerprint rows of b N on the ghost-number-2
    classes annihilate every row of the b-table (b^2 = 0 in four dimensions) and not a generic X."""
    _, rows = b_closure(toy_n, brst=toy_brst, fingerprinter=toy_table.fingerprinter)
    p = toy_table.field.p
    vals = rows.values.astype(object)
    assert not np.any((vals @ toy_table.matrix.T.astype(object)) % p)
    generic = np.arange(1, len(toy_n) + 1, dtype=object)
    assert np.any((vals @ generic) % p)
