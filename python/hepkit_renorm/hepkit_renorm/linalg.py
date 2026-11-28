"""Exact linear algebra over prime fields, rational reconstruction and interpolation.

Every decision of the Slavnov-Taylor layer (equality, rank, kernel, image, membership, existence and solution of
linear systems) is taken here, exactly, over ``GF(p)`` (:class:`hepkit_renorm.exact.PrimeField`). Rational results
are recovered by rational reconstruction, with the Chinese remainder theorem over several primes when one prime
is not enough, and accepted only once an additional prime confirms them. Dependence on parameters (couplings) is
recovered by dense polynomial interpolation with a degree bound.

Representation: matrices are ``numpy`` arrays of dtype ``int64`` with entries in ``[0, p)``; the prime must satisfy
``p < 2^31`` so that the product of two entries fits into 63 bits (all row operations are vectorised on that basis).
Vectors are 1-d arrays. ``a x = b`` means a matrix acting on column vectors: rows are equations (fingerprint
samples), columns are unknowns (functionals).

Probabilistic statements: the matrices are exact; only their entries come from random evaluation points (see
:mod:`hepkit_renorm.fingerprint`), which is where the Schwartz-Zippel bound applies.
"""

from __future__ import annotations

import itertools
from collections.abc import Callable, Iterator, Sequence
from fractions import Fraction
from math import gcd, isqrt, prod

import numpy as np

from .exact import PrimeField, default_field, find_prime, second_field

__all__ = [
    "InconsistentSystemError",
    "ReconstructionError",
    "as_field_matrix",
    "column_space_contains",
    "crt",
    "independent_columns",
    "interpolate",
    "monomials",
    "nullspace",
    "prime_fields",
    "rank",
    "rational_reconstruction",
    "reconstruct_polynomials",
    "reconstruct_rationals",
    "row_echelon",
    "solve",
    "split_gaussian",
]

_MAX_P = 2**31
_SQUARES = (-1, 2, 3, 5, 6, 7, 10)


class InconsistentSystemError(ValueError):
    """``a x = b`` has no solution; ``columns`` lists the right-hand sides (columns of ``b``) that fail."""

    def __init__(self, columns: Sequence[int]) -> None:
        super().__init__(f"inconsistent linear system (right-hand sides {list(columns)})")
        self.columns = list(columns)


class ReconstructionError(ArithmeticError):
    """A rational or polynomial reconstruction did not stabilise within the allowed number of primes."""


# ---------------------------------------------------------------------------------------------- matrices
def _check(field: PrimeField) -> int:
    if field.p >= _MAX_P:
        raise ValueError(f"linalg needs p < 2^31, got {field.p}")
    return field.p


def as_field_matrix(m: np.ndarray | Sequence, field: PrimeField) -> np.ndarray:
    """``m`` reduced modulo ``p`` as an ``int64`` array (accepts Python-int object arrays and nested lists)."""
    p = _check(field)
    a = np.asarray(m, dtype=object)
    return np.asarray(a % p, dtype=np.int64)


def row_echelon(m: np.ndarray, field: PrimeField) -> tuple[np.ndarray, list[int]]:
    """Reduced row-echelon form of ``m`` over ``GF(p)`` and its pivot columns.

    The returned matrix has ``rank`` rows (zero rows are dropped). Pivots are chosen as the first non-zero entry
    from the top in each column (deterministic): the pivot columns are the lexicographically first maximal set of
    independent columns.
    """
    p = _check(field)
    a = as_field_matrix(m, field).copy()
    if a.ndim != 2:
        raise ValueError("row_echelon needs a matrix")
    nrows, ncols = a.shape
    pivots: list[int] = []
    r = 0
    for c in range(ncols):
        if r == nrows:
            break
        nz = np.flatnonzero(a[r:, c])
        if nz.size == 0:
            continue
        piv = r + int(nz[0])
        if piv != r:
            a[[r, piv]] = a[[piv, r]]
        a[r] = a[r] * pow(int(a[r, c]), -1, p) % p
        col = a[:, c].copy()
        col[r] = 0
        rows = np.flatnonzero(col)
        if rows.size:
            a[rows] = (a[rows] - col[rows, None] * a[r][None, :] % p) % p
        pivots.append(c)
        r += 1
    return a[:r], pivots


def rank(m: np.ndarray, field: PrimeField) -> int:
    """Rank over ``GF(p)``."""
    if np.size(m) == 0:
        return 0
    return len(row_echelon(m, field)[1])


def independent_columns(m: np.ndarray, field: PrimeField) -> list[int]:
    """The lexicographically first maximal set of linearly independent columns (greedy from the left)."""
    if np.size(m) == 0:
        return []
    return row_echelon(m, field)[1]


def nullspace(m: np.ndarray, field: PrimeField) -> np.ndarray:
    """Basis of ``{v : m v = 0}`` as the rows of a ``(k, ncols)`` array. Vector ``j`` has a ``1`` at the ``j``-th
    free (non-pivot) column and zeros at the other free columns (deterministic, reduced form)."""
    p = _check(field)
    m = as_field_matrix(m, field)
    ncols = m.shape[1]
    if m.shape[0] == 0:
        return np.eye(ncols, dtype=np.int64)
    rref, pivots = row_echelon(m, field)
    free = [c for c in range(ncols) if c not in set(pivots)]
    basis = np.zeros((len(free), ncols), dtype=np.int64)
    for k, fc in enumerate(free):
        basis[k, fc] = 1
        for row, pc in enumerate(pivots):
            basis[k, pc] = (-int(rref[row, fc])) % p
    return basis


def solve(a: np.ndarray, b: np.ndarray, field: PrimeField, *, require_unique: bool = False) -> np.ndarray:
    """A solution ``x`` of ``a x = b`` over ``GF(p)`` (``b`` a vector or a matrix of right-hand sides).

    The particular solution returned has all free unknowns set to zero. Raises :class:`InconsistentSystemError`
    if some right-hand side is not in the column space of ``a`` (this is the existence test), and ``ValueError``
    with ``require_unique`` if ``a`` has a non-trivial kernel.
    """
    a = as_field_matrix(a, field)
    vec = np.ndim(b) == 1
    bm = as_field_matrix(b, field).reshape(a.shape[0], -1)
    n = a.shape[1]
    rref, pivots = row_echelon(np.concatenate([a, bm], axis=1), field)
    bad = [c - n for c in pivots if c >= n]
    if bad:
        raise InconsistentSystemError(bad)
    if require_unique and len(pivots) < n:
        raise ValueError(f"solution not unique: kernel of dimension {n - len(pivots)}")
    x = np.zeros((n, bm.shape[1]), dtype=np.int64)
    for row, pc in enumerate(pivots):
        x[pc] = rref[row, n:]
    return x[:, 0] if vec else x


def column_space_contains(a: np.ndarray, b: np.ndarray, field: PrimeField) -> bool:
    """True if every column of ``b`` (or the vector ``b``) lies in the column space of ``a``."""
    try:
        solve(a, b, field)
    except InconsistentSystemError:
        return False
    return True


def split_gaussian(plus: np.ndarray | int, minus: np.ndarray | int, field: PrimeField) -> tuple[np.ndarray, np.ndarray]:
    """Real and imaginary parts from the two embeddings of ``Q(i)``.

    ``plus`` is a value computed with ``i -> field.i`` and ``minus`` the same computation with ``i -> -field.i``
    (all other inputs being images of rationals, e.g. :class:`hepkit_renorm.fingerprint.Fingerprinter` with
    ``i_sign = +1`` and ``-1`` on the same samples). Then ``plus = re + i im`` with ``re, im`` the images of the
    real and imaginary parts: ``re = (plus + minus)/2``, ``im = (plus - minus)/(2 i)``.
    """
    p = _check(field)
    plus, minus = np.asarray(plus, dtype=np.int64) % p, np.asarray(minus, dtype=np.int64) % p
    half = field(Fraction(1, 2))
    re = (plus + minus) % p * half % p
    im = (plus - minus) % p * (half * field.inv(field.i) % p) % p
    return re, im


# ---------------------------------------------------------------------------------------------- reconstruction
def rational_reconstruction(x: int, modulus: int) -> Fraction:
    """The unique ``n/d`` with ``|n|, d <= sqrt(modulus/2)`` and ``n = x d (mod modulus)`` (Wang's algorithm).

    Raises :class:`ReconstructionError` if there is none.
    """
    bound = isqrt(modulus // 2)
    r0, r1, s0, s1 = modulus, x % modulus, 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1, s0, s1 = r1, r0 - q * r1, s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound or gcd(r1, abs(s1)) != 1:
        raise ReconstructionError(f"{x} has no rational reconstruction modulo {modulus}")
    return Fraction(r1, s1)


def crt(residues: Sequence[int], moduli: Sequence[int]) -> tuple[int, int]:
    """Chinese remainder theorem: ``(x, M)`` with ``x = r_k (mod m_k)`` and ``M = prod m_k`` (coprime moduli)."""
    x, m = 0, 1
    for r, mk in zip(residues, moduli):
        t = (int(r) - x) * pow(m, -1, mk) % mk
        x, m = x + m * t, m * mk
    return x % m, m


def prime_fields(squares: Sequence[int] = _SQUARES) -> Iterator[PrimeField]:
    """Deterministic sequence of distinct primes ``p < 2^31`` in which all ``squares`` are squares, starting with
    :func:`~hepkit_renorm.exact.default_field` and :func:`~hepkit_renorm.exact.second_field` (for the default
    ``squares``)."""
    seen: set[int] = set()
    if tuple(squares) == _SQUARES:
        for f in (default_field(), second_field()):
            seen.add(f.p)
            yield f
    start = 2**30 + 2**25
    while True:
        p = find_prime(start, squares)
        if p >= _MAX_P:
            raise ReconstructionError("ran out of primes below 2^31")
        start = p + 1
        if p not in seen:
            seen.add(p)
            yield PrimeField(p)


def _reconstruct_all(values: list[np.ndarray], moduli: list[int]) -> np.ndarray | None:
    shape = values[0].shape
    out = np.empty(shape, dtype=object)
    for idx in np.ndindex(shape):
        x, m = crt([int(v[idx]) for v in values], moduli)
        try:
            out[idx] = rational_reconstruction(x, m)
        except ReconstructionError:
            return None
    return out


def reconstruct_rationals(compute: Callable[[PrimeField], np.ndarray | int], *, max_primes: int = 8,
                          squares: Sequence[int] = _SQUARES) -> np.ndarray:
    """Rational numbers from their images: ``compute(field)`` returns the images (an array or an int) of fixed
    rationals in ``GF(p)``; primes are added (CRT) until the reconstruction from the first ``k`` primes is confirmed
    by prime ``k + 1``. Returns an object array of :class:`fractions.Fraction` (0-d for a scalar).

    Raises :class:`ReconstructionError` after ``max_primes`` primes.
    """
    values: list[np.ndarray] = []
    moduli: list[int] = []
    candidate = None
    for field in itertools.islice(prime_fields(squares), max_primes):
        v = np.array(np.asarray(compute(field), dtype=object) % field.p, dtype=object)
        if candidate is not None:
            check = np.array([field(q) for q in candidate.flat], dtype=object).reshape(candidate.shape)
            if np.array_equal(check, v):
                return candidate
        values.append(v)
        moduli.append(field.p)
        candidate = _reconstruct_all(values, moduli)
    raise ReconstructionError(f"no stable rational reconstruction with {max_primes} primes")


# ---------------------------------------------------------------------------------------------- interpolation
def monomials(n_vars: int, degree: int, homogeneous: bool = False) -> list[tuple[int, ...]]:
    """Exponent vectors of total degree ``<= degree`` (``== degree`` if ``homogeneous``), in graded
    lexicographic order."""
    out = []
    for total in ([degree] if homogeneous else range(degree + 1)):
        for c in itertools.combinations_with_replacement(range(n_vars), total):
            e = [0] * n_vars
            for v in c:
                e[v] += 1
            out.append(tuple(e))
    return sorted(set(out), key=lambda e: (sum(e), [-x for x in e]))


def interpolate(evaluate: Callable[[PrimeField, Sequence[int]], np.ndarray | int], n_vars: int, degree: int,
                field: PrimeField, *, homogeneous: bool = False, seed: int = 0,
                extra_points: int = 3) -> tuple[list[tuple[int, ...]], np.ndarray]:
    """Dense interpolation over ``GF(p)``: ``evaluate(field, point)`` returns a vector of values of polynomials of
    total degree ``<= degree`` in ``n_vars`` variables at ``point``. Returns the exponent list and the coefficient
    matrix ``c[monomial, output]``.

    ``len(monomials) + extra_points`` random points are used (deterministic for a given ``seed``); the extra points
    test the degree bound: if the overdetermined system is inconsistent, :class:`InconsistentSystemError` is
    raised (the values are not polynomials of that degree), and a degenerate sample raises ``ValueError``.
    """
    p = _check(field)
    monos = monomials(n_vars, degree, homogeneous)
    rng = np.random.default_rng([seed, p, n_vars, degree])
    npts = len(monos) + extra_points
    points = rng.integers(1, p, size=(npts, n_vars), dtype=np.int64)
    rows = [np.atleast_1d(np.asarray(evaluate(field, [int(x) for x in pt]), dtype=object) % p) for pt in points]
    b = as_field_matrix(np.array(rows, dtype=object), field)
    a = np.array([[prod(pow(int(x), e, p) for x, e in zip(pt, mono)) % p for mono in monos] for pt in points],
                 dtype=np.int64)
    return monos, solve(a, b, field, require_unique=True)


def reconstruct_polynomials(evaluate: Callable[[PrimeField, Sequence[int]], np.ndarray | int], n_vars: int,
                            degree: int, *, homogeneous: bool = False, seed: int = 0, max_primes: int = 8,
                            squares: Sequence[int] = _SQUARES) -> list[dict[tuple[int, ...], Fraction]]:
    """Polynomials with rational coefficients from their values: :func:`interpolate` in each prime of
    :func:`prime_fields`, coefficients reconstructed by :func:`reconstruct_rationals` (CRT, confirmed by an extra
    prime). Returns one ``{exponents: coefficient}`` dict (non-zero coefficients only) per output."""
    def coefficients(field: PrimeField) -> np.ndarray:
        return interpolate(evaluate, n_vars, degree, field, homogeneous=homogeneous, seed=seed)[1]

    monos = monomials(n_vars, degree, homogeneous)
    coeffs = reconstruct_rationals(coefficients, max_primes=max_primes, squares=squares)
    return [{monos[k]: coeffs[k, j] for k in range(len(monos)) if coeffs[k, j] != 0} for j in range(coeffs.shape[1])]
