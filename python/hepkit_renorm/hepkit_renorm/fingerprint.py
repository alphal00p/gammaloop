"""Exact plane-wave fingerprints of local functionals and vertex functions.

A fingerprint is the value of a vertex function on random external data, computed exactly in ``GF(p)``. Linear
relations between functionals (equality, rank, kernel, image, membership, solutions of linear systems) are then
decided by :mod:`hepkit_renorm.linalg` on fingerprint matrices: rows are samples, grouped by class of external legs,
columns are functionals.

What is evaluated
-----------------
For a leg class ``(Phi_1, ..., Phi_n)`` (an ordered tuple of fields and sources) and a sample (momenta ``p_k`` with
``sum p_k = 0``, polarisations ``e_k`` over the index slots of ``Phi_k``), the fingerprint of ``Y`` is::

    V_Y(e_1, p_1; ...; e_n, p_n) = coefficient of theta_1 ... theta_n in Y[sum_k theta_k e_k exp(-i p_k.x)]

with ``(2 pi)^d delta(sum p)`` stripped, i.e. the vertex function of :meth:`LocalFunctional.vertex_function`
(report 02_method sec:m-objects; Feynman rule ``i V``) contracted with the polarisations. ``theta_k`` is odd for odd
fields; the Grassmann sign of every term is the one of that single code path. Exchanging two odd legs of the class
flips the sign. Any expression in the vocabulary of :mod:`hepkit_renorm.tensors` whose free labels are the labels of
the legs (as returned by :func:`hepkit_renorm.functional.standard_legs`) can be evaluated (:meth:`Fingerprinter.
evaluate`), for example a vertex function produced elsewhere.

Explicit components (binding conventions, those of the stand-alone engines ``st/engine.py`` and
``phase4/pole_consistency/engine_d.py``)
----------------------------------------------------------------------------------------------------------------
* Dimension ``D = 4 + 2k`` (``D = 4, 6, 8, ...``); the symbol ``d`` takes the value ``D``. Lorentz components
  ``0..3`` are barred (four-dimensional), ``4..D-1`` hatted. Metric ``g = diag(+1, -1, ..., -1)``.
* All Lorentz components are contravariant: ``p_k^mu``, ``e_k^mu``, ``gamma^mu``, ``eps^{mu nu rho sigma}``
  (``eps^{0123} = +1``, non-zero only on barred components), ``g^{mu nu}``, ``gbar^{mu nu}``, ``ghat^{mu nu}``. Every
  label that occurs twice is contracted with ``g_{mu nu}``. (``st/engine.py`` stores lower components: its
  ``e_mu``, ``p_mu`` are ``g e^mu``, ``g p^mu`` here.)
* Dirac matrices: the Dirac representation ``gamma^0 = diag(1, 1, -1, -1)``, ``gamma^j = [[0, sigma_j], [-sigma_j,
  0]]``, ``gamma5 = i gamma^0 gamma^1 gamma^2 gamma^3``. For ``D > 4``: ``Gamma^mu = gamma^mu x 1`` (``mu < 4``),
  ``Gamma^(4+j) = gamma5 x i E_j`` with the Euclidean Clifford generators ``E_j`` of dimension ``2^k``
  (Jordan-Wigner; ``E_1 = sigma_1``, ``E_2 = sigma_2`` for ``k = 1``), and ``gamma5 x 1`` as ``gamma5``; hence
  ``{gammabar, gamma5} = 0``, ``[gammahat, gamma5] = 0``. The spinor dimension is ``4 * 2^k``, with components in
  Kronecker order ``s = s4 * 2^k + s'``. Every closed spinor cycle (a trace) is multiplied by
  ``2^-k``, so that ``tr 1 = 4`` as in :class:`hepkit_renorm.tensors.Simplifier`.
* Group components: the bases of :mod:`hepkit_renorm.groups` (generalised Gell-Mann matrices / 2 in Gell-Mann order;
  combined adjoint index of a product group per factor); fundamental indices ``0..N-1``; flavour ``0..N_f-1``.
* Numbers live in ``GF(p)`` (:class:`hepkit_renorm.exact.PrimeField`, ``p < 2^31``) with the embedding
  ``i -> i_sign * field.i`` and the fixed square roots of :mod:`hepkit_renorm.groups`. With ``i_sign = -1`` every
  ``i`` (numbers, Dirac matrices, generators: ``T -> T^T``) is replaced by ``-i``: on the same samples this gives the
  Galois-conjugate value, from which :func:`hepkit_renorm.linalg.split_gaussian` extracts real and imaginary parts.
* Parameters (couplings, masses, ``xi``, flavour matrices) take deterministic random non-zero values unless given;
  ``conj(y)`` of a complex parameter is an independent variable (polynomial identities in ``y, ybar``).

Samples are deterministic: sample ``j`` of a class depends only on ``(seed, class, salt, j, p, D)``, so a shorter
fingerprint is a prefix of a longer one, and the samples do not depend on ``i_sign`` or on the parameter values.
Exactness: every number is exact; a relation that holds on the samples holds identically with probability at least
``1 - deg/p`` per sample (Schwartz-Zippel), with ``deg`` the polynomial degree in the random data.
"""

from __future__ import annotations

import hashlib
import itertools
import weakref
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from functools import cached_property

import numpy as np
from symbolica import AtomType, Expression

from . import tensors as tn
from .exact import PrimeField, default_field
from .functional import Density, Leg, LocalFunctional, Monomial, standard_legs
from .groups import complex_rational, su_d_tensor, su_generators, su_structure_constants
from .spec import Model

__all__ = ["DEFAULT_DIMENSIONS", "FingerprintMatrix", "Fingerprinter", "Representation", "Samples", "is_zero",
           "shared_fingerprinter"]

_BATCH = "z"
_LETTERS = list("abcdefghijklmnopqrstuvwxyABCDEFGHIJKLMNOPQRSTUVWXYZ")   # "z" is the sample axis
_LIMB = 1 << 16
_CONJ = "symbolica::conj"


# ---------------------------------------------------------------------------------------------- modular einsum
def _pair(a: np.ndarray, sa: str, b: np.ndarray, sb: str, keep: set[str], p: int) -> tuple[np.ndarray, str]:
    """``sum`` over the letters of ``sa + sb`` not in ``keep`` of ``a * b``, modulo ``p`` (``b`` split into 16-bit
    limbs so that no intermediate exceeds 63 bits)."""
    out = "".join(dict.fromkeys(c for c in sa + sb if c in keep))
    dims = {c: n for s, x in ((sa, a), (sb, b)) for c, n in zip(s, x.shape)}
    length = 1
    for c in set(sa + sb) - set(out):
        length *= dims[c]
    if length >= _LIMB:
        raise ValueError(f"contraction of length {length} exceeds the exact int64 bound")
    spec = f"{sa},{sb}->{out}"
    lo = np.einsum(spec, a, b & (_LIMB - 1)) % p
    hi = np.einsum(spec, a, b >> 16) % p
    return (lo + hi * _LIMB) % p, out


def _contract(ops: list[np.ndarray], subs: list[str], out: str, p: int) -> np.ndarray:
    """Exact ``einsum(subs -> out)`` over ``GF(p)`` along numpy's greedy contraction path."""
    ops, subs = list(ops), list(subs)
    if len(ops) > 1:
        path = np.einsum_path(",".join(subs) + "->" + out, *ops, optimize="greedy")[0][1:]
        for step in path:
            sel = [(ops.pop(i), subs.pop(i)) for i in sorted(step, reverse=True)]
            a, sa = sel[0]
            for k in range(1, len(sel)):
                b, sb = sel[k]
                keep = set("".join(subs)) | set(out) | set("".join(s for _, s in sel[k + 1:]))
                a, sa = _pair(a, sa, b, sb, keep, p)
            ops.append(a)
            subs.append(sa)
    return np.einsum(f"{subs[0]}->{out}", ops[0]) % p


def _rng(*key: object) -> np.random.Generator:
    digest = hashlib.sha256(repr(key).encode()).digest()
    return np.random.default_rng(int.from_bytes(digest[:8], "little"))


# ---------------------------------------------------------------------------------------------- representation
@dataclass(frozen=True)
class Representation:
    """Explicit Lorentz and Dirac tensors in dimension ``D`` over ``GF(p)`` (module docstring)."""

    field: PrimeField
    dimension: int = 4
    i_sign: int = 1

    def __post_init__(self) -> None:
        if self.dimension < 4 or self.dimension % 2:
            raise ValueError("the explicit dimension must be even and >= 4")
        if self.field.p >= 2**31:
            raise ValueError("fingerprints need p < 2^31")
        if self.i_sign not in (1, -1):
            raise ValueError("i_sign is +1 or -1")

    @property
    def p(self) -> int:
        return self.field.p

    @cached_property
    def i(self) -> int:
        """The image of the imaginary unit."""
        return self.field.i if self.i_sign == 1 else self.p - self.field.i

    def complex(self, re: int | Fraction | str, im: int | Fraction | str = 0) -> int:
        """Image of ``re + i im``."""
        return (self.field(re) + self.field(im) * self.i) % self.p

    def number(self, num: Expression) -> int:
        """Image of a Symbolica complex rational."""
        return self.complex(*complex_rational(num))

    @cached_property
    def metric(self) -> np.ndarray:
        """``diag(g)``: ``[1, -1, ..., -1]`` modulo ``p``."""
        return np.array([1] + [self.p - 1] * (self.dimension - 1), dtype=np.int64)

    @cached_property
    def bar(self) -> np.ndarray:
        return np.array([1] * 4 + [0] * (self.dimension - 4), dtype=np.int64)

    @cached_property
    def hat(self) -> np.ndarray:
        return 1 - self.bar

    @property
    def hat_pairs(self) -> int:
        """``k = (D - 4)/2``."""
        return (self.dimension - 4) // 2

    @property
    def spinor_dim(self) -> int:
        return 4 << self.hat_pairs

    @cached_property
    def trace_norm(self) -> int:
        """``2^-k``: the factor of every closed spinor cycle."""
        return self.field.inv(1 << self.hat_pairs)

    def _mat(self, rows: Sequence[Sequence[int]]) -> np.ndarray:
        return np.array(rows, dtype=object) % self.p

    def _kron(self, *ms: np.ndarray) -> np.ndarray:
        out = np.array([[1]], dtype=object)
        for m in ms:
            out = np.kron(out, m) % self.p
        return out

    @cached_property
    def gamma(self) -> np.ndarray:
        """``gamma^mu`` as ``[mu, s, t]`` (``int64``)."""
        p, i = self.p, self.i
        one, zero = self._mat([[1, 0], [0, 1]]), self._mat([[0, 0], [0, 0]])
        sig = [self._mat([[0, 1], [1, 0]]), self._mat([[0, -i], [i, 0]]), self._mat([[1, 0], [0, -1]])]
        g4 = [np.block([[one, zero], [zero, -one]]) % p] + [np.block([[zero, s], [-s, zero]]) % p for s in sig]
        g5 = (i * (g4[0].dot(g4[1]).dot(g4[2]).dot(g4[3]))) % p
        k = self.hat_pairs
        eye_k = np.eye(1 << k, dtype=object)
        full = [self._kron(g, eye_k) for g in g4]
        for e in _euclidean_clifford(k, sig):
            full.append(self._kron(g5, (i * e) % p))
        return np.array([np.asarray(g, dtype=np.int64) for g in full], dtype=np.int64)

    @cached_property
    def gamma5(self) -> np.ndarray:
        g = [x.astype(object) for x in self.gamma[:4]]
        return np.asarray((self.i * g[0].dot(g[1]).dot(g[2]).dot(g[3])) % self.p, dtype=np.int64)

    @cached_property
    def identity(self) -> np.ndarray:
        return np.eye(self.spinor_dim, dtype=np.int64)

    @cached_property
    def projectors(self) -> tuple[np.ndarray, np.ndarray]:
        """``(PL, PR) = ((1 - gamma5)/2, (1 + gamma5)/2)``."""
        half = self.field(Fraction(1, 2))
        pl = (self.identity - self.gamma5) % self.p * half % self.p
        pr = (self.identity + self.gamma5) % self.p * half % self.p
        return pl, pr

    @cached_property
    def eps(self) -> np.ndarray:
        """``eps^{mu nu rho sigma}``, ``eps^{0123} = +1``."""
        e = np.zeros((self.dimension,) * 4, dtype=np.int64)
        for perm in itertools.permutations(range(4)):
            inv = sum(1 for x, y in itertools.combinations(perm, 2) if x > y)
            e[perm] = 1 if inv % 2 == 0 else self.p - 1
        return e

    def metric_tensor(self, kind: str) -> np.ndarray:
        """``g``, ``gbar`` or ``ghat`` with upper indices (diagonal)."""
        mask = {"g": 1, "gbar": self.bar, "ghat": self.hat}[kind]
        return np.diag(self.metric * mask % self.p).astype(np.int64)


def _euclidean_clifford(k: int, sig: Sequence[np.ndarray]) -> list[np.ndarray]:
    """``2k`` hermitian matrices of dimension ``2^k`` with ``{E_a, E_b} = 2 delta_ab`` (Jordan-Wigner)."""
    out = []
    for m in range(k):
        for s in (sig[0], sig[1]):
            mats = [sig[2]] * m + [s] + [np.eye(2, dtype=object)] * (k - m - 1)
            e = np.array([[1]], dtype=object)
            for x in mats:
                e = np.kron(e, x)
            out.append(e)
    return out


# ---------------------------------------------------------------------------------------------- samples
@dataclass(frozen=True)
class Samples:
    """External data for ``n`` evaluations of one leg class.

    ``legs``: the legs (fields, free labels, momentum numbers); ``momenta[k]``: ``(n, D)`` contravariant components
    of the momentum number ``k``; ``polarisations[j]``: ``(n, *slot dims)`` components of leg ``j`` (slot order of the
    field). All entries are elements of ``GF(p)``.
    """

    legs: tuple[Leg, ...]
    momenta: Mapping[int, np.ndarray]
    polarisations: tuple[np.ndarray, ...]

    @property
    def n(self) -> int:
        return self.polarisations[0].shape[0] if self.polarisations else next(iter(self.momenta.values())).shape[0]


@dataclass(frozen=True)
class FingerprintMatrix:
    """A fingerprint matrix: ``values[row, column]`` with the rows grouped by leg class.

    ``layout`` lists ``(class, number of samples)`` in row order; pass it to :meth:`Fingerprinter.matrix` to
    fingerprint further functionals on exactly the same rows.
    """

    field: PrimeField
    layout: tuple[tuple[tuple[str, ...], int], ...]
    values: np.ndarray

    def block(self, leg_class: Sequence[str]) -> np.ndarray:
        """The rows of one class."""
        start = 0
        for cls, n in self.layout:
            if cls == tuple(leg_class):
                return self.values[start:start + n]
            start += n
        raise KeyError(leg_class)


# ---------------------------------------------------------------------------------------------- compiled terms
@dataclass(frozen=True)
class _Operand:
    kind: str                     # const | pol | mom | sp
    letters: str
    data: object = None           # array (const), leg index (pol), momentum number (mom), (kind, k, l) (sp)


@dataclass(frozen=True)
class _Term:
    coefficient: int
    operands: tuple[_Operand, ...]


class Fingerprinter:
    """Exact evaluation of vertex functions of a model on plane-wave samples (module docstring).

    Parameters
    ----------
    model:
        the :class:`~hepkit_renorm.spec.Model` (field slots, group factors, parameters).
    field:
        the prime field (default :func:`~hepkit_renorm.exact.default_field`).
    dimension:
        explicit dimension ``D`` (4, 6, 8, ...).
    seed:
        seed of the samples and of the random parameter values.
    values:
        parameter values overriding the random ones: ``{name: value}`` with the specification's parameter and
        coupling names, ``"conj:<name>"`` for the conjugate of a complex parameter, and ``numpy`` arrays
        ``(N_f, N_f)`` for flavour matrices; scalars are ints / ``Fraction`` (images are taken).
    i_sign:
        embedding of ``i`` (``+1`` or ``-1``).
    """

    def __init__(self, model: Model, field: PrimeField | None = None, dimension: int = 4, *, seed: int = 0,
                 values: Mapping[str, object] | None = None, i_sign: int = 1) -> None:
        self.model = model
        self.rep = Representation(field or default_field(), dimension, i_sign)
        self.seed = seed
        self._overrides = dict(values or {})
        self._values = self._parameter_values(self._overrides)
        self._cache: dict[tuple, list[_Term]] = {}
        self._splits: dict[int, tuple[weakref.ref, dict[tuple[str, ...], LocalFunctional]]] = {}

    @property
    def field(self) -> PrimeField:
        return self.rep.field

    @property
    def dimension(self) -> int:
        return self.rep.dimension

    def replace(self, *, field: PrimeField | None = None, dimension: int | None = None, seed: int | None = None,
                values: Mapping[str, object] | None = None, i_sign: int | None = None) -> Fingerprinter:
        """A copy with some settings changed (``values`` are merged into the current overrides)."""
        merged = dict(self._overrides)
        merged.update(values or {})
        return Fingerprinter(self.model, field or self.field, dimension or self.dimension,
                             seed=self.seed if seed is None else seed, values=merged,
                             i_sign=self.rep.i_sign if i_sign is None else i_sign)

    # ------------------------------------------------------------------ parameters
    def _parameter_values(self, overrides: Mapping[str, object]) -> dict[str, object]:
        m, F = self.model, self.field
        names = [f.coupling for f in m.factors.values()] + list(m.parameters)
        unknown = {k.removeprefix("conj:") for k in overrides} - set(names)
        if unknown:
            raise KeyError(f"unknown parameters {sorted(unknown)}")
        out: dict[str, object] = {tn.d.get_name(): self.dimension % F.p}
        nf = m.generations
        for name in names:
            par = m.parameters.get(name)
            matrix = par is not None and par.flavour and nf > 1
            keys = [name] + ([] if par is None or par.real else [f"conj:{name}"])
            for key in keys:
                shape = (nf, nf) if matrix else ()
                if key in overrides:
                    v = overrides[key]
                    val = (np.vectorize(lambda x: F(x), otypes=[object])(np.asarray(v, dtype=object)).astype(np.int64)
                           if matrix else F(v))
                else:
                    r = _rng("parameter", self.seed, key, F.p).integers(1, F.p, size=shape, dtype=np.int64)
                    val = r if matrix else int(r)
                sym = m.param(name).get_name()
                out[sym if key == name else f"{_CONJ}:{sym}"] = val
        return out

    # ------------------------------------------------------------------ samples
    def slot_dims(self, field_name: str) -> tuple[int, ...]:
        """Component dimensions of the slots of a field."""
        return tuple(self._space_dim(s.space.kind, s.space.factor) for s in self.model.fields[field_name].slots)

    def _space_dim(self, kind: str, factor: str | None) -> int:
        if kind == "lorentz":
            return self.dimension
        if kind == "spinor":
            return self.rep.spinor_dim
        if kind == "flavour":
            return self.model.generations
        fac = self.model.factors[factor]
        return fac.dim_adjoint if kind == "adjoint" else fac.n

    def legs(self, leg_class: Sequence[str]) -> list[Leg]:
        """The standard legs of a class (:func:`~hepkit_renorm.functional.standard_legs`)."""
        return standard_legs(self.model.fields, list(leg_class))

    def samples(self, leg_class: Sequence[str], n: int, salt: int = 0) -> Samples:
        """``n`` random samples of a class: momenta with ``sum p = 0``, polarisations uniform in ``GF(p)``."""
        legs = self.legs(leg_class)
        p, D = self.field.p, self.dimension
        cls = tuple(leg_class)
        moms = np.zeros((n, len(legs), D), dtype=np.int64)
        pols = [np.zeros((n,) + self.slot_dims(leg.field), dtype=np.int64) for leg in legs]
        for j in range(n):
            rng = _rng("sample", self.seed, cls, salt, j, p, D)
            moms[j, :-1] = rng.integers(0, p, size=(len(legs) - 1, D), dtype=np.int64)
            moms[j, -1] = (-moms[j, :-1].sum(axis=0)) % p
            for k, pol in enumerate(pols):
                pol[j] = rng.integers(0, p, size=pol.shape[1:], dtype=np.int64)
        return Samples(tuple(legs), {leg.momentum: moms[:, k] for k, leg in enumerate(legs)}, tuple(pols))

    def explicit_samples(self, leg_class: Sequence[str], momenta: np.ndarray,
                         polarisations: Sequence[np.ndarray]) -> Samples:
        """Samples from given field elements: ``momenta`` of shape ``(n, legs, D)`` (contravariant, conserved),
        ``polarisations[j]`` of shape ``(n, *slot_dims(leg j))`` (use :meth:`Representation.complex` for the images
        of complex rationals)."""
        legs = self.legs(leg_class)
        p = self.field.p
        moms = np.asarray(momenta, dtype=np.int64) % p
        if moms.shape[1:] != (len(legs), self.dimension):
            raise ValueError(f"momenta must have shape (n, {len(legs)}, {self.dimension})")
        if np.any(moms.sum(axis=1) % p):
            raise ValueError("momenta are not conserved")
        pols = tuple(np.asarray(e, dtype=np.int64) % p for e in polarisations)
        for leg, e in zip(legs, pols):
            if e.shape[1:] != self.slot_dims(leg.field):
                raise ValueError(f"polarisation of {leg.field} must have shape (n, {self.slot_dims(leg.field)})")
        return Samples(tuple(legs), {leg.momentum: moms[:, k] for k, leg in enumerate(legs)}, pols)

    # ------------------------------------------------------------------ evaluation
    def evaluate(self, expr: Expression, samples: Samples) -> np.ndarray:
        """Values ``(n,)`` of an expression whose free labels are the labels of ``samples.legs``, contracted with the
        polarisations (a vertex function, momenta ``p(k, mu)`` with the legs' momentum numbers)."""
        key = (expr.to_canonical_string(), tuple((l.field, tuple(x.to_canonical_string() for x in l.labels),
                                                   l.momentum) for l in samples.legs))
        terms = self._cache.get(key)
        if terms is None:
            terms = [self._compile(t, samples.legs) for t in _terms(expr.expand())]
            self._cache[key] = terms
        p = self.field.p
        total = np.zeros(samples.n, dtype=np.int64)
        for term in terms:
            if term.coefficient:
                total = (total + self._run(term, samples) * term.coefficient) % p
        return total

    def functional(self, functional: LocalFunctional, samples: Samples) -> np.ndarray:
        """Fingerprint ``(n,)`` of a local functional on the samples of one leg class."""
        piece = self._split(functional).get(tuple(sorted(leg.field for leg in samples.legs)))
        return np.zeros(samples.n, dtype=np.int64) if piece is None else self._evaluate_piece(piece, samples)

    def matrix(self, functionals: Sequence[LocalFunctional], classes: Iterable[Sequence[str]] | None = None,
               samples_per_class: int | Mapping[tuple[str, ...], int] | None = None, *,
               layout: Sequence[tuple[Sequence[str], int]] | None = None, salt: int = 0) -> FingerprintMatrix:
        """Fingerprint matrix of ``functionals`` (columns) on the samples of the given classes (rows).

        Default classes: the union of the leg classes of the functionals. Default number of samples per class: the
        number of functionals that contribute to the class plus two (enough for an exact rank with probability
        ``>= 1 - deg * n / p``). ``layout`` (from a previous :class:`FingerprintMatrix`) fixes classes and counts.
        """
        pieces = [self._split(f) for f in functionals]
        if layout is None:
            if classes is None:
                order = self.model.fields.order
                found = {c for pc in pieces for piece in pc.values() for c in piece.leg_classes()}
                classes = sorted(found, key=lambda c: [order(x) for x in c])
            classes = [tuple(c) for c in classes]
            if samples_per_class is None:
                counts = {c: sum(tuple(sorted(c)) in pc for pc in pieces) + 2 for c in classes}
            elif isinstance(samples_per_class, int):
                counts = {c: samples_per_class for c in classes}
            else:
                counts = {c: samples_per_class[c] for c in classes}
            layout = [(c, counts[c]) for c in classes]
        layout = tuple((tuple(c), int(n)) for c, n in layout)
        blocks = []
        for cls, n in layout:
            smp = self.samples(cls, n, salt)
            key = tuple(sorted(cls))
            cols = [self._evaluate_piece(pc[key], smp) if key in pc else np.zeros(n, dtype=np.int64)
                    for pc in pieces]
            blocks.append(np.stack(cols, axis=1) if cols else np.zeros((n, 0), dtype=np.int64))
        values = np.concatenate(blocks, axis=0) if blocks else np.zeros((0, len(functionals)), dtype=np.int64)
        return FingerprintMatrix(self.field, layout, values)

    def _split(self, functional: LocalFunctional) -> dict[tuple[str, ...], LocalFunctional]:
        """:func:`_by_class`, cached per functional object (densities are immutable)."""
        hit = self._splits.get(id(functional))
        if hit is not None and hit[0]() is functional:
            return hit[1]
        split = _by_class(functional)
        self._splits[id(functional)] = (weakref.ref(functional), split)
        return split

    def _evaluate_piece(self, piece: LocalFunctional, samples: Samples) -> np.ndarray:
        v = piece.vertex_function(list(samples.legs), None, momentum_conservation=False)
        return self.evaluate(v, samples)

    def vanishes(self, functional: LocalFunctional, samples_per_class: int = 2, salt: int = 0) -> bool:
        """True if every vertex function of the functional vanishes on the samples (zero up to integration by parts
        and every identity of the explicit algebra)."""
        classes = functional.leg_classes()
        if not classes:
            return True
        m = self.matrix([functional], classes, samples_per_class, salt=salt)
        return not np.any(m.values)

    def equal(self, a: LocalFunctional, b: LocalFunctional, samples_per_class: int = 2) -> bool:
        """``a == b`` as functionals (:meth:`vanishes` of the difference)."""
        return self.vanishes(LocalFunctional.integrate(a - b), samples_per_class)

    # ------------------------------------------------------------------ compilation
    def _compile(self, term: Expression, legs: Sequence[Leg]) -> _Term:
        p = self.field.p
        coeff = 1
        raw: list[tuple[str, list[str], object]] = []
        label_space: dict[str, str] = {}

        def labels_of(args: Sequence[Expression]) -> list[str]:
            out = []
            for a in args:
                k = a.to_canonical_string()
                label_space[k] = tn.label_space(a)
                out.append(k)
            return out

        for fa in (list(term) if term.get_type() == AtomType.Mul else [term]):
            reps = 1
            base = fa
            if fa.get_type() == AtomType.Pow:
                base, ex = list(fa)
                if ex.get_type() != AtomType.Num:
                    raise ValueError(f"non-numeric power {fa}")
                reps = int(ex.to_canonical_string())
            t = base.get_type()
            if t == AtomType.Num:
                coeff = coeff * pow(self.rep.number(base), reps, p) % p
            elif t == AtomType.Var or (t == AtomType.Fn and base.get_name() == _CONJ
                                       and next(iter(base)).get_type() == AtomType.Var):
                coeff = coeff * pow(self._scalar(base), reps, p) % p
            elif t == AtomType.Fn:
                if reps < 1:
                    raise ValueError(f"negative power of a tensor {fa}")
                for _ in range(reps):
                    raw.append(self._head(base, labels_of))
            else:
                raise ValueError(f"cannot evaluate {fa} (expand first)")
        leg_letters: dict[str, int] = {}
        for j, leg in enumerate(legs):
            keys = labels_of(leg.labels)
            raw.append(("pol", keys, j))
            for k in keys:
                leg_letters[k] = j
        counts: dict[str, int] = {}
        for _, keys, _ in raw:
            for k in keys:
                counts[k] = counts.get(k, 0) + 1
        bad = [k for k, c in counts.items() if c != 2]
        if bad:
            raise ValueError(f"labels {bad} do not occur exactly twice in {term} with the legs")
        letters = {k: _LETTERS[n] for n, k in enumerate(counts)} if len(counts) <= len(_LETTERS) else None
        if letters is None:
            raise ValueError("too many labels in one term")
        operands = []
        for kind, keys, data in raw:
            batch = _BATCH if kind in ("pol", "mom", "sp") else ""
            operands.append(_Operand(kind, batch + "".join(letters[k] for k in keys), data))
        for k, sp in label_space.items():
            if sp == "lorentz":
                operands.append(_Operand("const", letters[k], self.rep.metric))
        closed = _closed_spinor_cycles(raw, label_space, set(leg_letters))
        coeff = coeff * pow(self.rep.trace_norm, closed, p) % p
        return _Term(coeff, tuple(operands))

    def _scalar(self, e: Expression) -> int:
        if e.get_type() == AtomType.Fn:
            key = f"{_CONJ}:{next(iter(e)).get_name()}"
        else:
            key = e.get_name()
        if key not in self._values:
            raise KeyError(f"no value for the symbol {e}")
        v = self._values[key]
        if isinstance(v, np.ndarray):
            raise TypeError(f"{e} is a flavour matrix and needs flavour labels")
        return int(v)

    def _head(self, fn: Expression, labels_of) -> tuple[str, list[str], object]:
        name = fn.get_name()
        args = list(fn)
        rep = self.rep
        if name == tn.p.get_name():
            return "mom", labels_of(args[1:]), int(args[0].to_canonical_string())
        if name in (tn.sp.get_name(), tn.spbar.get_name(), tn.sphat.get_name()):
            ks = tuple(int(a.to_canonical_string()) for a in args)
            return "sp", [], (name.split("::")[-1],) + ks
        if name == _CONJ:
            inner = args[0]
            key = f"{_CONJ}:{inner.get_name()}"
            return "const", labels_of(list(inner)), self._matrix_value(key, inner)
        if name in self._values:
            return "const", labels_of(args), self._matrix_value(name, fn)
        keys = labels_of(args)
        if name in (tn.g.get_name(), tn.gbar.get_name(), tn.ghat.get_name()):
            return "const", keys, rep.metric_tensor(name.split("::")[-1])
        if name == tn.eps.get_name():
            return "const", keys, rep.eps
        if name == tn.gamma.get_name():
            return "const", keys, rep.gamma
        if name == tn.gammabar.get_name():
            return "const", keys, rep.gamma * rep.bar[:, None, None]
        if name == tn.gammahat.get_name():
            return "const", keys, rep.gamma * rep.hat[:, None, None]
        if name == tn.gamma5.get_name():
            return "const", keys, rep.gamma5
        if name == tn.PL.get_name():
            return "const", keys, rep.projectors[0]
        if name == tn.PR.get_name():
            return "const", keys, rep.projectors[1]
        if name == tn.delta.get_name():
            kind, _, fac = tn.label_space(args[0]).partition(":")
            return "const", keys, np.eye(self._space_dim(kind, fac or None), dtype=np.int64)
        if name in (tn.T.get_name(), tn.f.get_name(), tn.dsym.get_name(), tn.epsf.get_name()):
            return "const", keys, self._group_array(name, tn.label_space(args[0]).split(":")[1])
        raise ValueError(f"cannot evaluate the head {name}")

    def _matrix_value(self, key: str, fn: Expression) -> np.ndarray:
        if key not in self._values:
            raise KeyError(f"no value for the tensor {fn}")
        v = self._values[key]
        if not isinstance(v, np.ndarray):
            raise TypeError(f"{fn}: parameter without flavour slots")
        return v

    @cached_property
    def _groups(self) -> dict[tuple[str, str], np.ndarray]:
        out = {}
        F = self.field
        for name, fac in self.model.factors.items():
            if fac.abelian:
                continue
            t = np.asarray(su_generators(fac.n, F), dtype=np.int64)
            out[(tn.T.get_name(), name)] = t if self.rep.i_sign == 1 else np.ascontiguousarray(t.transpose(0, 2, 1))
            out[(tn.f.get_name(), name)] = np.asarray(su_structure_constants(fac.n, F), dtype=np.int64)
            out[(tn.dsym.get_name(), name)] = np.asarray(su_d_tensor(fac.n, F), dtype=np.int64)
            e = np.zeros((fac.n,) * fac.n, dtype=np.int64)
            for perm in itertools.permutations(range(fac.n)):
                inv = sum(1 for x, y in itertools.combinations(perm, 2) if x > y)
                e[perm] = 1 if inv % 2 == 0 else F.p - 1
            out[(tn.epsf.get_name(), name)] = e
        return out

    def _group_array(self, head: str, factor: str) -> np.ndarray:
        return self._groups[(head, factor)]

    def _run(self, term: _Term, samples: Samples) -> np.ndarray:
        p = self.field.p
        ops, subs = [], []
        for op in term.operands:
            if op.kind == "const":
                ops.append(op.data)
            elif op.kind == "pol":
                ops.append(samples.polarisations[op.data])
            elif op.kind == "mom":
                ops.append(samples.momenta[op.data])
            else:
                kind, k, l = op.data
                mask = {"sp": 1, "spbar": self.rep.bar, "sphat": self.rep.hat}[kind]
                w = self.rep.metric * mask % p
                ops.append(_contract([samples.momenta[k], samples.momenta[l], w], ["zm", "zm", "m"], "z", p))
            subs.append(op.letters)
        if not any(_BATCH in s for s in subs):
            value = _contract(ops, subs, "", p) if ops else np.array(1, dtype=np.int64)
            return np.full(samples.n, int(value), dtype=np.int64)
        return _contract(ops, subs, _BATCH, p)


# ---------------------------------------------------------------------------------------------- certified zero test
DEFAULT_DIMENSIONS = (4, 6, 8)
_SHARED: dict[int, tuple[weakref.ref, dict[tuple[int, int], Fingerprinter]]] = {}


def shared_fingerprinter(model: Model, dimension: int = 4, field: PrimeField | None = None) -> Fingerprinter:
    """A :class:`Fingerprinter` with the default seed, cached per model, dimension and prime (its compiled terms
    are reused across calls)."""
    field = field or default_field()
    hit = _SHARED.get(id(model))
    if hit is None or hit[0]() is not model:
        hit = (weakref.ref(model), {})
        _SHARED[id(model)] = hit
    key = (dimension, field.p)
    if key not in hit[1]:
        hit[1][key] = Fingerprinter(model, field, dimension)
    return hit[1][key]


def is_zero(model: Model, functional: LocalFunctional, *, simplifier: tn.Simplifier | None = None,
            field: PrimeField | None = None, dimensions: Sequence[int] = DEFAULT_DIMENSIONS,
            samples_per_class: int = 2) -> bool:
    """The package's zero test of a local functional (zero up to integration by parts), class by class.

    1. If the vertex function of a class vanishes after the BMHV simplifier (default ``model.simplifier()``), or
       after explicit evaluation of its group tensors (:meth:`Model.vanishes_modulo_group`), the class is zero:
       both apply identities only, so a zero is a proof.
    2. Otherwise the class is fingerprinted exactly in every explicit dimension of ``dimensions`` (``d = D``,
       cheapest first). A non-zero value certifies that the functional is not zero (``False``). Vanishing on all
       samples in every dimension gives ``True``: an identity in the random data with probability ``>= 1 - deg/p``
       per sample (Schwartz-Zippel), for the dependence on ``d`` that the dimensions resolve. A structure that
       vanishes at ``D = 4, 6, 8`` without vanishing identically in ``d`` (a polynomial with three roots there, e.g.
       five or more antisymmetrised hatted indices) is not detected; pass more dimensions for such structures.
    """
    simp = simplifier or model.simplifier()
    field = field or default_field()
    pending: list[tuple[str, ...]] = []
    for cls in functional.leg_classes():
        legs = standard_legs(model.fields, cls)
        value = functional.vertex_function(legs, simp)
        if value.expand() == Expression.num(0):
            continue
        open_labels = [lab for leg in legs for lab, sl in zip(leg.labels, model.fields[leg.field].slots)
                       if sl.factor is not None]
        if not model.vanishes_modulo_group(value, open_labels, field):
            pending.append(cls)
    if not pending:
        return True
    split = _by_class(functional)
    rest = [split[tuple(sorted(c))] for c in pending]
    for dim in sorted(dimensions):
        fp = shared_fingerprinter(model, dim, field)
        if any(not fp.vanishes(piece, samples_per_class) for piece in rest):
            return False
    return True


def _terms(expr: Expression) -> list[Expression]:
    if expr.get_type() == AtomType.Add:
        return list(expr)
    if expr.get_type() == AtomType.Num and expr.to_canonical_string() == "0":
        return []
    return [expr]


def _by_class(functional: LocalFunctional) -> dict[tuple[str, ...], LocalFunctional]:
    """The functional split by field content: ``{sorted field names: functional of those monomials}``."""
    groups: dict[tuple[str, ...], list[Monomial]] = {}
    for m in functional.monomials():
        groups.setdefault(tuple(sorted(m.fields())), []).append(m)
    return {k: LocalFunctional.integrate(Density.from_monomials(functional.fields, ms)) for k, ms in groups.items()}


def _closed_spinor_cycles(raw: Sequence[tuple[str, list[str], object]], label_space: Mapping[str, str],
                          leg_labels: set[str]) -> int:
    """Number of connected components of spinor labels that touch no leg (closed traces)."""
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        while parent.setdefault(x, x) != x:
            x = parent[x]
        return x

    for _, keys, _ in raw:
        spin = [k for k in keys if label_space[k] == "spinor"]
        for a, b in itertools.pairwise(spin):
            parent[find(a)] = find(b)
    roots = {find(k) for k, sp in label_space.items() if sp == "spinor"}
    open_roots = {find(k) for k in leg_labels if label_space[k] == "spinor"}
    return len(roots - open_roots)
