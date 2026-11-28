"""The single implementations of two concepts used everywhere: exact elimination / reconstruction
(:mod:`hepkit_renorm.linalg`, to which :class:`hepkit_renorm.exact.PrimeField` delegates) and the certified zero test
of a local functional (:func:`hepkit_renorm.fingerprint.is_zero`, behind :func:`hepkit_renorm.brst.vanishes`)."""

from fractions import Fraction

import numpy as np
import pytest

from hepkit_renorm import Density, LocalFunctional, Occurrence, linalg
from hepkit_renorm import tensors as tn
from hepkit_renorm.brst import vanishes
from hepkit_renorm.exact import PrimeField, default_field, find_prime
from hepkit_renorm.fingerprint import is_zero


def _random_rank_deficient(field: PrimeField, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    u = [[int(x) for x in row] for row in rng.integers(0, 1000, size=(7, 3))]
    v = [[int(x) for x in row] for row in rng.integers(0, 1000, size=(3, 9))]
    return np.asarray(u, dtype=object).dot(np.asarray(v, dtype=object)) % field.p


@pytest.mark.parametrize("p", [default_field().p, find_prime(2**61, squares=(-1,))])
def test_prime_field_methods_are_the_linalg_ones(p):
    """PrimeField.rank / nullspace / row_reduce / reconstruct agree with linalg, also beyond 2^31 (object path)."""
    field = PrimeField(p)
    m = _random_rank_deficient(field, 3)
    assert field.rank(m) == linalg.rank(m, field) == 3
    ker = field.nullspace(m)
    assert len(ker) == 6 and all(isinstance(x, int) for x in ker[0])
    assert np.array_equal(np.array(ker, dtype=object) % p, np.asarray(linalg.nullspace(m, field), dtype=object))
    assert all(field.is_zero(field.matmul(m, v)) for v in ker)
    rref, pivots = field.row_reduce(m)
    assert rref.shape == m.shape and pivots == linalg.independent_columns(m, field) and not rref[3:].any()
    for q in (Fraction(-95, 162), Fraction(7, 3), Fraction(0)):
        assert field.reconstruct(field(q)) == linalg.rational_reconstruction(field(q), p) == q
    with pytest.raises(ValueError):
        field.reconstruct(field.sqrt(-1))


def test_large_prime_solve():
    field = PrimeField(find_prime(2**61, squares=(-1,)))
    m = _random_rank_deficient(field, 5)
    x0 = np.array([1, 2, 0, 5, 0, 0, 7, 1, 3], dtype=object)
    b = m.dot(x0) % field.p
    x = linalg.solve(m, b, field)
    assert x.dtype == object and not np.any((m.dot(x) - b) % field.p)


def _qbar_q(sm, q, coeff, derivatives) -> LocalFunctional:
    return LocalFunctional.integrate(Density.monomial(sm.fields, coeff, Occurrence("Qbar", q),
                                                      Occurrence("Q", (_T,) + q[1:], derivatives)))


_T, _U = tn.fresh(tn.SPINOR), tn.fresh(tn.SPINOR)


def test_certified_zero_test_sees_identities_the_simplifier_does_not_apply(sm):
    """{gamma^mu, gamma^nu} = 2 g^{mu nu}: the simplifier does not reorder gamma matrices, the fingerprints decide."""
    simp = sm.simplifier()
    q = sm.slot_labels("Q")
    s = q[0]
    mu, nu = tn.fresh(tn.LORENTZ), tn.fresh(tn.LORENTZ)
    anti = tn.gamma(mu, s, _U) * tn.gamma(nu, _U, _T) + tn.gamma(nu, s, _U) * tn.gamma(mu, _U, _T)
    identity = _qbar_q(sm, q, anti - 2 * tn.g(mu, nu) * tn.delta(s, _T), (mu, nu))
    assert not identity.vanishes(simp)                 # symbolic test: inconclusive
    assert is_zero(sm, identity) and vanishes(sm, identity)
    wrong = _qbar_q(sm, q, anti - tn.g(mu, nu) * tn.delta(s, _T), (mu, nu))
    assert not is_zero(sm, wrong) and not vanishes(sm, wrong)    # certified by a non-zero fingerprint
    # evanescent: zero in D = 4 only; the d-dimensional test sees it
    hat = _qbar_q(sm, q, tn.gammahat(mu, s, _T), (mu,))
    assert is_zero(sm, hat, dimensions=(4,)) and not is_zero(sm, hat)
