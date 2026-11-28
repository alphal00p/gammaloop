"""Exact arithmetic: prime fields, square roots, rational reconstruction.

Every rank, kernel, invariance or zero test of the package is decided exactly. Numbers that are not rational
(for example the ``sqrt(3)`` of the SU(3) generator ``T^8`` or the imaginary unit) are represented in a prime field
``GF(p)`` in which the required square roots exist; the map ``Q(i, sqrt(n1), ...) -> GF(p)`` is the ring
homomorphism fixed by the choice of roots made by :meth:`PrimeField.sqrt` (the smaller representative). Rational
results are recovered by rational reconstruction and confirmed on a second, independent prime by the caller when
needed.

Elements of ``GF(p)`` are plain Python ``int`` in ``[0, p)``; arrays of elements are ``numpy`` arrays of dtype
``object`` holding Python ints, so no overflow can occur.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from fractions import Fraction
from math import isqrt

import numpy as np

__all__ = ["PrimeField", "default_field", "find_prime", "is_prime", "second_field"]

_MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def is_prime(n: int) -> bool:
    """Deterministic Miller-Rabin test, exact for ``n < 3.3e24``."""
    if n < 2:
        return False
    for q in _MR_BASES:
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in _MR_BASES:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _legendre(a: int, p: int) -> int:
    return pow(a % p, (p - 1) // 2, p)


def find_prime(start: int, squares: Iterable[int] = (-1,)) -> int:
    """Smallest prime ``p >= start`` in which every integer of ``squares`` is a non-zero square.

    Deterministic: the same arguments always give the same prime.
    """
    req = tuple(squares)
    p = start
    while True:
        if is_prime(p) and all(a % p != 0 and _legendre(a, p) == 1 for a in req):
            return p
        p += 1


@dataclass(frozen=True)
class PrimeField:
    """The prime field ``GF(p)``.

    Conventions: :meth:`sqrt` returns the smaller of the two roots, so ``i = sqrt(-1)`` and every other root are
    fixed once per field; all identities of ``Q(i, sqrt(n), ...)`` hold in ``GF(p)`` under this embedding.
    """

    p: int

    def __post_init__(self) -> None:
        if not is_prime(self.p):
            raise ValueError(f"{self.p} is not prime")

    # ------------------------------------------------------------------ scalars
    def __call__(self, x: int | Fraction | str) -> int:
        """Image of an integer or rational number (``str`` is parsed by :class:`fractions.Fraction`)."""
        q = Fraction(x)
        if q.denominator % self.p == 0:
            raise ZeroDivisionError(f"denominator of {q} vanishes modulo {self.p}")
        return q.numerator % self.p * pow(q.denominator, -1, self.p) % self.p

    def inv(self, x: int) -> int:
        """Multiplicative inverse."""
        if x % self.p == 0:
            raise ZeroDivisionError("inverse of 0")
        return pow(x, -1, self.p)

    def sqrt(self, a: int | Fraction) -> int:
        """Square root (Tonelli-Shanks), the smaller representative of the two roots."""
        p = self.p
        a = self(a)
        if a == 0:
            return 0
        if _legendre(a, p) != 1:
            raise ValueError(f"{a} is not a square modulo {p}")
        q, s = p - 1, 0
        while q % 2 == 0:
            q //= 2
            s += 1
        z = 2
        while _legendre(z, p) != p - 1:
            z += 1
        m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
        while t != 1:
            i, t2 = 0, t
            while t2 != 1:
                t2 = t2 * t2 % p
                i += 1
            b = pow(c, 1 << (m - i - 1), p)
            m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
        return min(r, p - r)

    @property
    def i(self) -> int:
        """The imaginary unit ``sqrt(-1)``."""
        return self.sqrt(-1)

    def reconstruct(self, x: int) -> Fraction:
        """Rational reconstruction: the unique ``n/d`` with ``|n|, d <= sqrt(p/2)`` congruent to ``x``.

        Raises ``ValueError`` if no such fraction exists (then ``x`` is not the image of a small rational).
        """
        p = self.p
        bound = isqrt(p // 2)
        r0, r1, s0, s1 = p, x % p, 0, 1
        while r1 > bound:
            q = r0 // r1
            r0, r1, s0, s1 = r1, r0 - q * r1, s1, s0 - q * s1
        if s1 == 0 or abs(s1) > bound:
            raise ValueError(f"{x} has no rational reconstruction modulo {p}")
        f = Fraction(r1, s1)
        if self(f) != x % p:
            raise ValueError(f"{x} has no rational reconstruction modulo {p}")
        return f

    # ------------------------------------------------------------------ arrays
    def array(self, values: Sequence | np.ndarray) -> np.ndarray:
        """Object array of field elements from nested sequences of ints / Fractions."""
        a = np.array(values, dtype=object)
        flat = [self(v) for v in a.flat]
        return np.array(flat, dtype=object).reshape(a.shape)

    def zeros(self, shape: tuple[int, ...]) -> np.ndarray:
        """Array of zeros (Python ints)."""
        return np.zeros(shape, dtype=object)

    def reduce(self, a: np.ndarray) -> np.ndarray:
        """Entrywise reduction modulo ``p`` (for arrays produced by ``+``, ``*``, ``einsum``)."""
        return np.asarray(a, dtype=object) % self.p

    def einsum(self, spec: str, *ops: np.ndarray) -> np.ndarray:
        """``numpy.einsum`` on object arrays followed by reduction."""
        out = np.einsum(spec, *ops, optimize="greedy" if len(ops) > 2 else False)
        return self.reduce(np.asarray(out, dtype=object))

    def matmul(self, a: np.ndarray, b: np.ndarray) -> np.ndarray:
        return self.reduce(np.dot(a, b))

    def is_zero(self, a: np.ndarray | int) -> bool:
        if isinstance(a, np.ndarray):
            return all(int(v) % self.p == 0 for v in a.flat)
        return int(a) % self.p == 0

    def row_reduce(self, m: np.ndarray) -> tuple[np.ndarray, list[int]]:
        """Reduced row-echelon form and the list of pivot columns (Gauss-Jordan over ``GF(p)``)."""
        p = self.p
        a = [[int(v) % p for v in row] for row in np.asarray(m, dtype=object)]
        nrows = len(a)
        ncols = len(a[0]) if nrows else 0
        pivots: list[int] = []
        r = 0
        for c in range(ncols):
            piv = next((i for i in range(r, nrows) if a[i][c]), None)
            if piv is None:
                continue
            a[r], a[piv] = a[piv], a[r]
            inv = pow(a[r][c], -1, p)
            a[r] = [v * inv % p for v in a[r]]
            for i in range(nrows):
                if i != r and a[i][c]:
                    fac = a[i][c]
                    ri = a[i]
                    rr = a[r]
                    a[i] = [(x - fac * y) % p for x, y in zip(ri, rr)]
            pivots.append(c)
            r += 1
            if r == nrows:
                break
        return np.array(a, dtype=object).reshape(nrows, ncols), pivots

    def rank(self, m: np.ndarray) -> int:
        return len(self.row_reduce(m)[1])

    def nullspace(self, m: np.ndarray) -> list[np.ndarray]:
        """Basis of ``{v : m v = 0}``; vector ``k`` has a ``1`` at the ``k``-th free column (deterministic)."""
        m = np.asarray(m, dtype=object)
        ncols = m.shape[1]
        if m.shape[0] == 0:
            return [self.array([1 if j == k else 0 for j in range(ncols)]) for k in range(ncols)]
        rref, pivots = self.row_reduce(m)
        free = [c for c in range(ncols) if c not in pivots]
        basis = []
        for fc in free:
            v = [0] * ncols
            v[fc] = 1
            for row, pc in enumerate(pivots):
                v[pc] = (-int(rref[row][fc])) % self.p
            basis.append(np.array(v, dtype=object))
        return basis


def default_field() -> PrimeField:
    """The package's default field: a 31-bit prime with ``sqrt(-1)``, ``sqrt(2)``, ``sqrt(3)``, ``sqrt(5)``,
    ``sqrt(6)``, ``sqrt(7)`` and ``sqrt(10)`` (enough for the orthonormal generators of SU(N), N <= 5)."""
    return PrimeField(find_prime(2**31 - 2**20, squares=(-1, 2, 3, 5, 6, 7, 10)))


def second_field() -> PrimeField:
    """An independent prime with the same roots, for confirming reconstructed rationals."""
    return PrimeField(find_prime(2**30 + 2**19, squares=(-1, 2, 3, 5, 6, 7, 10)))
