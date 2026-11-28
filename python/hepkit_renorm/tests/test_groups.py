from fractions import Fraction

import numpy as np
import pytest

from hepkit_renorm import tensors as tn
from hepkit_renorm.groups import (
    GaugeFactor,
    GroupEvaluator,
    casimir_value,
    dynkin_index,
    generator_expr,
    invariant_tensors,
    product_generators,
    rep_generators,
    su_d_tensor,
    su_generators,
    su_structure_constants,
)


@pytest.mark.parametrize("n", [2, 3, 4, 5])
def test_lie_algebra_and_normalisation(field, n):
    F = field
    t = su_generators(n, F)
    f = su_structure_constants(n, F)
    comm = F.reduce(F.einsum("aij,bjk->abik", t, t) - F.einsum("bij,ajk->abik", t, t))
    assert F.is_zero(F.reduce(comm - F.einsum("abc,cik->abik", f, t) * F.i))
    trace = F.einsum("aij,bji->ab", t, t)
    assert F.is_zero(F.reduce(trace - np.eye(n * n - 1, dtype=object) * F(Fraction(1, 2))))
    assert F.is_zero(F.reduce(f + np.transpose(f, (1, 0, 2)))) and F.is_zero(F.reduce(f - np.transpose(f, (1, 2, 0))))
    fac = GaugeFactor("X", f"SU({n})", "gx")
    assert dynkin_index(fac, "F", F) == Fraction(1, 2) and dynkin_index(fac, "A", F) == n
    assert casimir_value(fac, "F", F) == Fraction(n * n - 1, 2 * n) and casimir_value(fac, "A", F) == n
    assert casimir_value(fac, "Fbar", F) == Fraction(n * n - 1, 2 * n)


def test_su3_conventions(field):
    """Gell-Mann ordering: f^{123} = 1, f^{458} = sqrt(3)/2, d^{118} = 1/sqrt(3), f^{147} = 1/2."""
    F = field
    f, d = su_structure_constants(3, F), su_d_tensor(3, F)
    assert f[0, 1, 2] == 1 and F.reconstruct(int(f[0, 3, 6])) == Fraction(1, 2)
    assert (2 * f[3, 4, 7]) ** 2 % F.p == 3
    assert d[0, 0, 7] ** 2 * 3 % F.p == 1
    assert F.is_zero(su_d_tensor(2, F))
    assert int(F.reduce(F.einsum("abc,abc->", d, d))) == F(Fraction(40, 3))  # d^{abc} d^{abc} = (N^2-4)(N^2-1)/N


def test_representations(field):
    F = field
    su3 = GaugeFactor("G", "SU(3)", "g_s")
    t, tb, ta = (rep_generators(su3, r, F) for r in ("F", "Fbar", "A"))
    assert F.is_zero(F.reduce(tb + np.transpose(t, (0, 2, 1))))
    comm = F.reduce(F.einsum("aij,bjk->abik", ta, ta) - F.einsum("bij,ajk->abik", ta, ta))
    assert F.is_zero(F.reduce(comm - F.einsum("abc,cik->abik", su_structure_constants(3, F), ta) * F.i))
    u1 = GaugeFactor("B", "U(1)", "gp")
    assert rep_generators(u1, Fraction(-1, 3), F)[0, 0, 0] == F(Fraction(-1, 3))
    with pytest.raises(ValueError):
        GaugeFactor("X", "SO(3)", "g")


def test_invariant_tensor_spaces(field):
    F = field
    su3, su2 = GaugeFactor("G", "SU(3)", "g_s"), GaugeFactor("W", "SU(2)", "g")
    assert len(invariant_tensors(su3, ["A", "A", "A"], F)) == 2  # f and d
    assert len(invariant_tensors(su2, ["A", "A", "A"], F)) == 1  # epsilon
    assert len(invariant_tensors(su3, ["Fbar", "F"], F)) == 1
    assert len(invariant_tensors(su3, ["F", "F"], F)) == 0
    assert len(invariant_tensors(su3, ["F", "F", "F"], F)) == 1
    (eps2,) = invariant_tensors(su2, ["F", "F"], F)
    assert eps2[0, 0] == 0 and (eps2[0, 1] + eps2[1, 0]) % F.p == 0
    assert len(invariant_tensors(su3, ["Fbar", "A", "F"], F)) == 1  # the generator itself


def test_product_generators_close(field):
    F = field
    facs = [GaugeFactor("G", "SU(3)", "g_s"), GaugeFactor("W", "SU(2)", "g"), GaugeFactor("B", "U(1)", "gp")]
    t = product_generators(facs, ["F", "F", Fraction(1, 6)], F)
    assert t.shape == (12, 6, 6)
    # generators of different factors commute
    assert F.is_zero(F.reduce(F.einsum("ij,jk->ik", t[0], t[8]) - F.einsum("ij,jk->ik", t[8], t[0])))


def test_symbolic_builders_evaluate_to_the_explicit_tensors(field):
    F = field
    su3 = GaugeFactor("G", "SU(3)", "g_s")
    ev = GroupEvaluator({"G": su3}, F)
    G = tn.factor_symbol("G")
    a, b, i, j = tn.adj(G, 1), tn.adj(G, 2), tn.fund(G, 1), tn.fund(G, 2)
    for rep, (x, y) in (("F", (i, j)), ("Fbar", (i, j)), ("A", (a, b))):
        c = tn.adj(G, 9)
        arr = ev.evaluate(generator_expr(su3, rep, c, x, y), [c, x, y])
        assert np.array_equal(arr, F.reduce(rep_generators(su3, rep, F)))
    ff = ev.evaluate(tn.f(a, tn.adj(G, 3), tn.adj(G, 4)) * tn.f(b, tn.adj(G, 3), tn.adj(G, 4)), [a, b])
    assert np.array_equal(ff, F.reduce(np.eye(8, dtype=object) * 3))
    assert generator_expr(GaugeFactor("B", "U(1)", "gp"), Fraction(2, 3), None, None, None).to_canonical_string() == "2/3"
