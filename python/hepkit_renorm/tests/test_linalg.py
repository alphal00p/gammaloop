"""Exact linear algebra, rational reconstruction and interpolation (hepkit_renorm.linalg)."""

from __future__ import annotations

import time
from fractions import Fraction

import numpy as np
import pytest

from hepkit_renorm import linalg as la
from hepkit_renorm.exact import default_field, second_field


def _random(field, shape, seed):
    return np.random.default_rng(seed).integers(0, field.p, size=shape, dtype=np.int64)


def test_rank_nullspace_solve_small(field):
    m = la.as_field_matrix([[1, 2, 3], [2, 4, 6], [1, 0, 1]], field)
    assert la.rank(m, field) == 2
    assert la.independent_columns(m, field) == [0, 1]
    ker = la.nullspace(m, field)
    assert ker.shape == (1, 3)
    assert not np.any(la.as_field_matrix(m.astype(object).dot(ker[0].astype(object)), field))
    # the reduced kernel vector is (-1, -1, 1): rational reconstruction recovers it
    assert [field.reconstruct(int(x)) for x in ker[0]] == [-1, -1, 1]
    x = la.solve(m, la.as_field_matrix([6, 12, 2], field), field)
    assert np.array_equal(la.as_field_matrix(m.astype(object).dot(x.astype(object)), field), [6, 12, 2])
    with pytest.raises(la.InconsistentSystemError) as err:
        la.solve(m, la.as_field_matrix([[6, 1], [12, 0], [2, 0]], field), field)
    assert err.value.columns == [1]
    assert not la.column_space_contains(m, la.as_field_matrix([1, 0, 0], field), field)
    with pytest.raises(ValueError):
        la.solve(m, la.as_field_matrix([6, 12, 2], field), field, require_unique=True)


def test_random_low_rank(field):
    a = _random(field, (300, 17), 1).astype(object)
    b = _random(field, (17, 40), 2).astype(object)
    m = la.as_field_matrix(a.dot(b), field)
    assert la.rank(m, field) == 17
    ker = la.nullspace(m, field)
    assert ker.shape == (23, 40)
    assert not np.any(la.as_field_matrix(m.astype(object).dot(ker.T.astype(object)), field))
    assert la.rank(ker, field) == 23


def test_performance_hundred_unknowns(field):
    """About a hundred unknowns against a few thousand rows (the size of the Slavnov-Taylor systems)."""
    rng = np.random.default_rng(3)
    a, b = rng.integers(0, 2**15, size=(3000, 90)), rng.integers(0, 2**15, size=(90, 110))
    m = la.as_field_matrix(a @ b, field)                       # rank 90 (entries < 2^37: exact in int64)
    x_true = rng.integers(0, 2**15, size=110)
    rhs = la.as_field_matrix(a @ (b @ x_true), field)
    t = time.perf_counter()
    assert la.rank(m, field) == 90
    x = la.solve(m, rhs, field)
    assert np.array_equal(la.as_field_matrix(m.astype(object).dot(x.astype(object)), field), rhs)
    assert time.perf_counter() - t < 10


def test_rational_reconstruction_and_crt(field):
    for q in (Fraction(0), Fraction(-3, 7), Fraction(8, 3), Fraction(-1201, 1999)):
        assert la.rational_reconstruction(field(q), field.p) == q
    f2 = second_field()
    x, m = la.crt([3, 5], [field.p, f2.p])
    assert x % field.p == 3 and x % f2.p == 5 and m == field.p * f2.p


def test_reconstruct_rationals_multi_prime():
    """Rationals far beyond the single-prime bound sqrt(p/2) need several primes (CRT) and are confirmed."""
    targets = [Fraction(123456789011, 98765432), Fraction(-7, 3), Fraction(0), Fraction(10**12 + 39, 17)]
    calls = []

    def compute(f):
        calls.append(f.p)
        return np.array([f(q) for q in targets], dtype=object)

    out = la.reconstruct_rationals(compute)
    assert list(out) == targets
    assert len(calls) >= 3 and len(set(calls)) == len(calls)
    assert la.reconstruct_rationals(lambda f: f(Fraction(5, 2))) == Fraction(5, 2)
    with pytest.raises(la.ReconstructionError):          # never confirmed by a second prime
        la.reconstruct_rationals(compute, max_primes=1)


def test_split_gaussian(field):
    z = (Fraction(2, 3), Fraction(-5, 7))
    plus = (field(z[0]) + field(z[1]) * field.i) % field.p
    minus = (field(z[0]) - field(z[1]) * field.i) % field.p
    re, im = la.split_gaussian(plus, minus, field)
    assert field.reconstruct(int(re)) == z[0] and field.reconstruct(int(im)) == z[1]


def test_interpolation_and_polynomial_reconstruction():
    # f0 = 3/4 g^2 h - 2 h^3 + 1/3,  f1 = g h / 5
    def evaluate(f, pt):
        g, h = pt
        return [(f(Fraction(3, 4)) * g * g * h - 2 * h**3 + f(Fraction(1, 3))) % f.p, g * h * f(Fraction(1, 5)) % f.p]

    polys = la.reconstruct_polynomials(evaluate, 2, 3)
    assert polys == [{(2, 1): Fraction(3, 4), (0, 3): Fraction(-2), (0, 0): Fraction(1, 3)}, {(1, 1): Fraction(1, 5)}]
    assert la.monomials(2, 2, homogeneous=True) == [(2, 0), (1, 1), (0, 2)]
    with pytest.raises(la.InconsistentSystemError):        # degree bound too small is detected
        la.interpolate(evaluate, 2, 2, default_field())
