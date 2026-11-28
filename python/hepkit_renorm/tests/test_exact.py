from fractions import Fraction

import numpy as np
import pytest

from hepkit_renorm.exact import PrimeField, default_field, find_prime, is_prime, second_field


def test_primality_and_search_are_deterministic():
    assert [n for n in range(30) if is_prime(n)] == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    assert is_prime(2**61 - 1) and not is_prime(2**61 + 1)
    p = find_prime(10**6, squares=(-1, 3))
    assert p == find_prime(10**6, squares=(-1, 3)) and p % 12 == 1


def test_roots_and_reconstruction():
    F = default_field()
    assert F.i * F.i % F.p == F.p - 1
    for n in (2, 3, 5, 6, 7, 10):
        r = F.sqrt(n)
        assert r * r % F.p == n and r <= F.p - r
    for q in (Fraction(-95, 162), Fraction(11, 216), Fraction(0), Fraction(-7)):
        assert F.reconstruct(F(q)) == q
    with pytest.raises(ValueError):
        F.reconstruct(F.sqrt(3))  # sqrt(3) is not a small rational
    assert default_field().p != second_field().p


def test_linear_algebra_over_the_field():
    F = PrimeField(101)
    m = F.array([[1, 2, 3], [2, 4, 6], [1, 0, 1]])
    assert F.rank(m) == 2
    (v,) = F.nullspace(m)
    assert F.is_zero(F.matmul(m, v))
    assert F.is_zero(F.einsum("ij,j->i", m, v))
    assert F.zeros((2, 2)).shape == (2, 2) and isinstance(F.zeros((2,))[0], int)
    assert np.array_equal(F.reduce(np.array([-1, 102], dtype=object)), np.array([100, 1], dtype=object))
