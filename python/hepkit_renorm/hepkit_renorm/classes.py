"""Classes of external legs of the 1PI functions the recipe needs, derived from the field content of a model.

A *leg class* is a multiset of fields and sources ``(Phi_1, ..., Phi_n)`` (stored in the canonical field order of
:class:`~hepkit_renorm.fields.FieldTable`) of fixed ghost number ``g`` and total dimension ``sum_j dim Phi_j <=
D_g``. It is kept if some local monomial with exactly these fields is non-zero, Lorentz invariant and rigidly
invariant. The classes depend on power counting and symmetry only, not on the loop order: the same sets serve
every order ``l`` of report 02_method sec:m-algorithm (step 2: ghost number 0; step 7: ghost number 1; ghost
number 2 for the consistency condition ``b X_l = 0``). Which classes have a graph at a given order follows from
diagram generation, not from here.

Power counting
--------------
``dim s`` (the dimension of the BRST operator) is read off the source terms: ``dim X^* + dim X + dim s = 4`` for
every pair (:func:`brst_dimension`; 1 for the report's conventions). A local functional of ghost number ``g``
built from ``s`` and the action has dimension ``<= D_g = 4 + g dim s`` (:func:`dimension_bound`: 4, 5, 6 for
``g = 0, 1, 2``). The superficial degree of divergence of a one-loop graph of the class is
``D_g - sum_j dim Phi_j`` (report sec:integrals; ``g = 1``: one vertex of dimension 5).

External legs (:func:`external_fields`)
--------------------------------------
All fields and sources except those of kind ``auxiliary`` (``b``: gauge condition, tree level only) and
``antighost`` (ghost equation: ``omegabar`` enters only with ``d^mu`` through the shifted gauge source, which the
source leg ``src_F`` stands for), report eq:linear-ids. With ``physical=True`` (default) a fermion leg is the
physical chiral component (report ``Pphys``: the sterile partners have no external legs): the spinor slot of
``Psi_f`` carries the Lorentz representation of ``P_c Psi_f``, that of ``Psibar_f`` the one of ``Psibar_f
P_cbar``, and a source carries the representation of its field (``X^*`` contracts with ``s X`` without a Dirac
matrix, report eq:Phistar). ``physical=False`` counts full Dirac spinors (both chiralities).

Selection (exact, no hand rules)
--------------------------------
A local monomial of the class modulo total derivatives is the same as its vertex function: a polynomial in the
independent momenta (``sum p = 0``) with values in the tensor product of the slot spaces, which is

* Lorentz covariant: invariant under ``so(1,3)`` acting on momenta, Lorentz and spinor slots (``epsilon``
  included: the group is the connected one);
* rigidly invariant: annihilated by every generator ``sum_k T^a_(k)`` of every gauge factor acting on the
  product representation (the null space used by :func:`hepkit_renorm.groups.invariant_tensors`);
* graded symmetric: exchanging two identical legs (their slots and momenta) gives the sign ``-1`` if they are
  odd and ``+1`` otherwise; flavour indices are free (a spurion coefficient).

The dimension of that space at momentum degree ``k`` is computed exactly by characters
(:func:`structure_counts`): with ``H`` the group of permutations of identical legs and ``eps(sigma)`` the sign on
the odd legs,

    N_k = 1/|H| sum_{sigma in H} eps(sigma) tr(sigma | L_k^{so(4)}) prod_F tr(sigma | R_F^{G_F}) N_f(sigma),

where ``tr(sigma | V^G) = int_G dg tr(sigma g | V)`` (``sigma`` commutes with the group) is the constant term of
a Laurent polynomial on the maximal torus times the Weyl denominator ``prod_{alpha > 0} (1 - e^{-alpha})``
(Weyl integration formula), and a cycle of length ``l`` of ``sigma`` on identical legs with character ``chi``
contributes ``chi(g^l)``. For ``sigma = 1`` the rigid factor is the dimension of the null space of the generators
on the product representation (tested against :func:`~hepkit_renorm.groups.invariant_tensors`); the character
form also gives the traces of the leg permutations on that null space, which the graded symmetry needs, and it is
fast for products of many adjoints (``SU(3)`` adjoint^6 at ghost number 2), where a dense null space is not.
Tori: ``SU(N)`` with weights ``e_i`` of the fundamental (``F``), ``-e_i`` (``Fbar``) and ``chi_A = chi_F chi_Fbar
- 1``; ``U(1)`` with the charge; ``so(4) = su(2) + su(2)`` with half-unit exponents, vector ``(1/2, 1/2)``, spinors
``(1/2, 0)`` (``P_L``) and ``(0, 1/2)`` (``P_R``). The momentum polynomials enter through
``sum_k t^k tr(sigma g | Sym^k(W)) = 1/det(1 - t sigma g | W)`` with ``W = (C^n / C(1,..,1)) x vector``.

Grouping (:func:`kind_signature`, :func:`sector`)
-------------------------------------------------
The signature of a class is the multiset of field kinds, written with ``c`` (ghost), ``F`` (fermion or
conjugate), ``S`` (scalar or conjugate), ``V`` (gauge field) in this order, followed by ``*`` plus the letter of
the field of each source (``*V`` = gauge source, ``*c`` = ghost source, ``*F``, ``*S``): ``"cFFV"``,
``"SSVV"``, ``"cc*V"``. Its sector is the signature, or ``"sources"`` for every class with a source. The report's
eighteen classes (tab:classes) are unions of sectors at ghost numbers 0 and 1.
"""

from __future__ import annotations

import functools
import itertools
from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from fractions import Fraction
from math import factorial, floor

from .fields import FieldInfo
from .groups import GaugeFactor
from .spec import Model

__all__ = [
    "LegClass",
    "brst_dimension",
    "by_sector",
    "conjugate_fields",
    "dimension_bound",
    "external_fields",
    "invariant_trace",
    "kind_signature",
    "leg_classes",
    "sector",
    "structure_counts",
]

EXCLUDED_KINDS = ("auxiliary", "antighost")
_KIND_LETTER = {"ghost": "c", "fermion": "F", "fermion_bar": "F", "scalar": "S", "scalar_conj": "S", "gauge": "V"}
_LETTER_ORDER = "cFSV"

# ---------------------------------------------------------------------------------------------- Laurent polynomials
# A Laurent polynomial is a dict {exponent tuple: int}; exponents are ints or Fractions.
Laurent = dict


def _mul(a: Laurent, b: Laurent) -> Laurent:
    out: dict = {}
    for ea, ca in a.items():
        for eb, cb in b.items():
            e = tuple(x + y for x, y in zip(ea, eb))
            out[e] = out.get(e, 0) + ca * cb
    return {e: c for e, c in out.items() if c}


def _add(a: Laurent, b: Laurent, scale: int = 1) -> Laurent:
    out = dict(a)
    for e, c in b.items():
        out[e] = out.get(e, 0) + scale * c
    return {e: c for e, c in out.items() if c}


def _power_map(a: Laurent, k: int) -> Laurent:
    """``chi(g) -> chi(g^k)``: every exponent times ``k``."""
    return {tuple(k * x for x in e): c for e, c in a.items()}


def _one(nvars: int) -> Laurent:
    return {(0,) * nvars: 1}


def _from_weights(weights: Iterable[tuple]) -> Laurent:
    out: dict = {}
    for w in weights:
        out[w] = out.get(w, 0) + 1
    return out


def _constant_term(a: Laurent) -> int:
    for e, c in a.items():
        if not any(e):
            return c
    return 0


# ---------------------------------------------------------------------------------------------- tori
class _Torus:
    """Maximal torus of a compact group: number of variables and the Weyl denominator."""

    def __init__(self, nvars: int, weyl: Laurent) -> None:
        self.nvars = nvars
        self.weyl = weyl

    def integrate(self, character: Laurent) -> int:
        """``int_G chi``: the constant term of ``chi`` times the Weyl denominator (``chi`` a class function)."""
        return _constant_term(_mul(character, self.weyl))


def _su_weight(n: int, vec: Sequence[int]) -> tuple[int, ...]:
    """A weight ``sum a_i e_i`` of SU(n) modulo ``e_1 + .. + e_n = 0``, in the coordinates ``a_i - a_n``."""
    return tuple(vec[i] - vec[n - 1] for i in range(n - 1))


@functools.cache
def _su_torus(n: int) -> _Torus:
    weyl = _one(n - 1)
    for i, j in itertools.combinations(range(n), 2):
        root = [0] * n
        root[i], root[j] = -1, 1                      # e^{-alpha}, alpha = e_i - e_j > 0
        weyl = _mul(weyl, {(0,) * (n - 1): 1, _su_weight(n, root): -1})
    return _Torus(n - 1, weyl)


@functools.cache
def _su_character(n: int, rep: str) -> tuple[tuple[tuple[int, ...], int], ...]:
    unit = [[1 if k == i else 0 for k in range(n)] for i in range(n)]
    fund = _from_weights(_su_weight(n, v) for v in unit)
    if rep == "F":
        out = fund
    elif rep == "Fbar":
        out = _power_map(fund, -1)
    elif rep == "A":
        out = _add(_mul(fund, _power_map(fund, -1)), _one(n - 1), -1)
    elif rep == "1":
        out = _one(n - 1)
    else:
        raise ValueError(f"unknown SU({n}) representation {rep!r}")
    return tuple(sorted(out.items()))


_U1_TORUS = _Torus(1, _one(1))

# so(4) = su(2) + su(2): exponents in half units (spin 1/2 has exponents +-1), roots +-2.
_LORENTZ_TORUS = _Torus(2, _mul({(0, 0): 1, (-2, 0): -1}, {(0, 0): 1, (0, -2): -1}))
_VECTOR = {(s, t): 1 for s in (1, -1) for t in (1, -1)}
_SPINOR = {"L": {(1, 0): 1, (-1, 0): 1}, "R": {(0, 1): 1, (0, -1): 1}}
_DIRAC = _add(_SPINOR["L"], _SPINOR["R"])
_FLIP = {"L": "R", "R": "L"}


# ---------------------------------------------------------------------------------------------- per-field data
def brst_dimension(model: Model) -> Fraction:
    """``dim s`` from the source terms: ``4 - dim X^* - dim X``, the same for every pair (checked)."""
    values = {4 - fi.dimension - model.fields[fi.source_of].dimension
              for fi in model.fields.values() if fi.is_source}
    if len(values) != 1:
        raise ValueError(f"inconsistent dimensions of the BRST sources: dim s in {sorted(values)}")
    return Fraction(values.pop())


def dimension_bound(model: Model, ghost_number: int) -> Fraction:
    """``D_g = 4 + g dim s``: the largest dimension of a local functional of ghost number ``g``."""
    return 4 + ghost_number * brst_dimension(model)


def external_fields(model: Model) -> list[str]:
    """The fields and sources that occur as external legs (module docstring), in canonical order."""
    return [n for n, fi in model.fields.items() if fi.kind not in EXCLUDED_KINDS]


def _spinor_chirality(model: Model, fi: FieldInfo) -> str:
    """Physical chirality label of the spinor slot (module docstring): ``P_c`` for ``Psi``, ``P_cbar`` for
    ``Psibar``, that of the field for a source."""
    base = model.fields[fi.source_of] if fi.is_source else fi
    if base.chirality not in _FLIP:
        raise NotImplementedError(f"fermion {base.name!r} without a definite chirality")
    return base.chirality if base.kind == "fermion" else _FLIP[base.chirality]


def _lorentz_character(model: Model, name: str, physical: bool) -> Laurent:
    fi = model.fields[name]
    out = _one(2)
    for slot in fi.slots:
        if slot.space.kind == "lorentz":
            out = _mul(out, _VECTOR)
        elif slot.space.kind == "spinor":
            out = _mul(out, _SPINOR[_spinor_chirality(model, fi)] if physical else _DIRAC)
    return out


def _torus(factor: GaugeFactor) -> _Torus:
    return _U1_TORUS if factor.abelian else _su_torus(factor.n)


def _rep_character(factor: GaugeFactor, rep: str | Fraction) -> Laurent:
    if factor.abelian:
        return {(Fraction(rep),): 1}
    return dict(_su_character(factor.n, str(rep)))


def _group_character(model: Model, name: str, factor: GaugeFactor) -> Laurent:
    fi = model.fields[name]
    if factor.abelian:
        return _rep_character(factor, fi.charges.get(factor.name, 0))
    out = _one(factor.n - 1)
    for slot in fi.slots:
        if slot.factor == factor.name and slot.space.kind in ("adjoint", "fund"):
            out = _mul(out, _rep_character(factor, slot.rep))
    return out


def _cycle_trace(torus: _Torus, cycles: Iterable[tuple[Laurent, int]]) -> int:
    """``tr(sigma | V^G)`` for ``sigma`` with the given cycles ``(character of the legs, length)``."""
    chi = _one(torus.nvars)
    for character, length in cycles:
        chi = _mul(chi, _power_map(character, length))
    return torus.integrate(chi)


def invariant_trace(factor: GaugeFactor, cycles: Sequence[tuple[str | Fraction, int]]) -> int:
    """Trace of a permutation of identical legs on the invariant subspace (the null space of the generators) of a
    product representation of one factor. ``cycles``: ``(representation, length)`` per cycle; ``length``
    identical legs in that representation are permuted cyclically (all lengths 1: the dimension of the null
    space, as :func:`hepkit_renorm.groups.invariant_tensors`). Representations as in :mod:`hepkit_renorm.groups`
    (``"1"``, ``"F"``, ``"Fbar"``, ``"A"``, or a U(1) charge)."""
    return _cycle_trace(_torus(factor), [(_rep_character(factor, rep), length) for rep, length in cycles])


def _flavour_dim(model: Model, name: str) -> int:
    out = 1
    for slot in model.fields[name].slots:
        if slot.space.kind == "flavour":
            out *= model.generations
    return out


def _partitions(m: int, largest: int | None = None) -> Iterable[tuple[int, ...]]:
    largest = m if largest is None else largest
    if m == 0:
        yield ()
        return
    for first in range(min(m, largest), 0, -1):
        for rest in _partitions(m - first, first):
            yield (first,) + rest


def _centraliser(parts: Sequence[int]) -> int:
    """``z_lambda = prod_l l^{m_l} m_l!``: the order of the centraliser of a permutation of cycle type ``parts``."""
    out = 1
    for length, mult in Counter(parts).items():
        out *= length**mult * factorial(mult)
    return out


def _momentum_series(cycles: Sequence[tuple[Laurent, int]], max_degree: int) -> list[Laurent]:
    """Coefficients ``t^0 .. t^max_degree`` of ``det(1 - t g | V) prod_c chi_c(g^l) / det(1 - t^l g^l | V)``."""
    series = [_one(2)] + [{} for _ in range(max_degree)]

    def times(poly: Sequence[tuple[int, Laurent]]) -> None:   # multiply by sum_j t^{deg_j} c_j (truncated)
        nonlocal series
        new = [{} for _ in range(max_degree + 1)]
        for k, s in enumerate(series):
            if not s:
                continue
            for deg, c in poly:
                if k + deg <= max_degree:
                    new[k + deg] = _add(new[k + deg], _mul(s, c))
        series = new

    for w in _VECTOR:
        times([(0, _one(2)), (1, {w: -1})])
    for chi, length in cycles:
        times([(0, _power_map(chi, length))])
        for w in _VECTOR:
            geometric = [(length * j, {tuple(length * j * x for x in w): 1}) for j in range(max_degree // length + 1)]
            times(geometric)
    return series


# ---------------------------------------------------------------------------------------------- counting
def structure_counts(model: Model, fields: Sequence[str], max_degree: int, *, physical: bool = True) -> tuple[int, ...]:
    """``(N_0, .., N_max_degree)``: the number of independent Lorentz-covariant, rigidly invariant, graded-symmetric
    vertex structures of the legs ``fields`` that are homogeneous of degree ``k`` in the momenta (module
    docstring). ``sum_k N_k > 0`` iff a non-vanishing local monomial with these fields and at most
    ``max_degree`` derivatives exists. Exact (integer character sums)."""
    if max_degree < 0:
        return ()
    content = Counter(fields)
    names = sorted(content, key=model.fields.order)
    lorentz = {n: _lorentz_character(model, n, physical) for n in names}
    factors = list(model.factors.values())
    group = {(n, f.name): _group_character(model, n, f) for n in names for f in factors}
    totals = [Fraction(0)] * (max_degree + 1)
    for cycle_types in itertools.product(*[list(_partitions(content[n])) for n in names]):
        weight = Fraction(1)
        sign = 1
        flavour = 1
        for n, parts in zip(names, cycle_types):
            weight /= _centraliser(parts)
            if model.fields[n].odd:
                sign *= (-1) ** sum(length - 1 for length in parts)
            flavour *= _flavour_dim(model, n) ** len(parts)
        rigid = 1
        for f in factors:
            rigid *= _cycle_trace(_torus(f), [(group[(n, f.name)], length)
                                              for n, parts in zip(names, cycle_types) for length in parts])
            if not rigid:
                break
        if not rigid:
            continue
        cycles = [(lorentz[n], length) for n, parts in zip(names, cycle_types) for length in parts]
        for k, coefficient in enumerate(_momentum_series(cycles, max_degree)):
            if coefficient:
                totals[k] += weight * sign * flavour * rigid * _LORENTZ_TORUS.integrate(coefficient)
    if any(t.denominator != 1 for t in totals):
        raise ArithmeticError(f"non-integer structure count {totals} for {list(fields)}")
    return tuple(int(t) for t in totals)


# ---------------------------------------------------------------------------------------------- classes
@dataclass(frozen=True)
class LegClass:
    """One class of external legs.

    ``fields``: the legs in canonical field order (the key used by :mod:`hepkit_renorm.functional` and
    :mod:`hepkit_renorm.fingerprint`). ``ghost_number``; ``dimension`` = ``sum_j dim Phi_j``; ``bound`` = ``D_g``;
    ``structures[k]``: number of independent vertex structures of momentum degree ``k`` (:func:`structure_counts`);
    ``signature`` and ``sector`` (:func:`kind_signature`, :func:`sector`).
    """

    fields: tuple[str, ...]
    ghost_number: int
    dimension: Fraction
    bound: Fraction
    structures: tuple[int, ...]
    signature: str
    sector: str

    @property
    def degree(self) -> Fraction:
        """Superficial degree of divergence ``D_g - sum_j dim Phi_j`` of a one-loop graph of the class."""
        return self.bound - self.dimension

    @property
    def n_structures(self) -> int:
        """Number of independent local monomials of the class (all momentum degrees)."""
        return sum(self.structures)

    @property
    def key(self) -> str:
        """``"g<ghost number>_<sector>"``, e.g. ``"g0_FFV"``, ``"g1_sources"``."""
        return f"g{self.ghost_number}_{self.sector}"


def kind_signature(model: Model, fields: Iterable[str]) -> str:
    """The multiset of field kinds of a class (module docstring), e.g. ``"cFFV"`` or ``"cc*FF"``."""
    plain, sources = [], []
    for name in fields:
        fi = model.fields[name]
        if fi.is_source:
            sources.append(_KIND_LETTER[model.fields[fi.source_of].kind])
        elif fi.kind in _KIND_LETTER:
            plain.append(_KIND_LETTER[fi.kind])
        else:
            raise ValueError(f"field {name!r} of kind {fi.kind!r} is not an external leg")
    order = _LETTER_ORDER.index
    return "".join(sorted(plain, key=order)) + "".join("*" + s for s in sorted(sources, key=order))


def sector(model: Model, fields: Iterable[str]) -> str:
    """``"sources"`` for a class with a BRST source, otherwise its :func:`kind_signature`."""
    fields = list(fields)
    return "sources" if any(model.fields[n].is_source for n in fields) else kind_signature(model, fields)


def conjugate_fields(model: Model, fields: Iterable[str]) -> tuple[str, ...]:
    """The legs of the hermitian-conjugate class, in canonical order."""
    return tuple(sorted((model.fields[n].conjugate for n in fields), key=model.fields.order))


def _multisets(model: Model, names: Sequence[str], bound: Fraction) -> Iterable[tuple[str, ...]]:
    """All non-empty multisets of ``names`` with total dimension ``<= bound`` (canonical order)."""
    dims = [model.fields[n].dimension for n in names]
    if any(d <= 0 for d in dims):
        raise ValueError("every external field must have positive dimension for a finite enumeration")

    def rec(start: int, budget: Fraction, chosen: tuple[str, ...]) -> Iterable[tuple[str, ...]]:
        if chosen:
            yield chosen
        for k in range(start, len(names)):
            if dims[k] <= budget:
                yield from rec(k, budget - dims[k], chosen + (names[k],))

    return rec(0, Fraction(bound), ())


def leg_classes(model: Model, ghost_number: int, *, bound: Fraction | int | None = None, physical: bool = True,
                fields: Sequence[str] | None = None) -> list[LegClass]:
    """All leg classes of ghost number ``ghost_number`` with ``sum dim <= bound`` (default
    :func:`dimension_bound`) that admit a non-vanishing, Lorentz-invariant, rigidly invariant local monomial.

    ``fields`` restricts the legs (default :func:`external_fields`); ``physical`` selects the physical chiral
    components of the fermions (module docstring). Deterministic order: by number of legs, then canonical field
    order.
    """
    bound = dimension_bound(model, ghost_number) if bound is None else Fraction(bound)
    names = sorted(external_fields(model) if fields is None else fields, key=model.fields.order)
    out = []
    for legs in _multisets(model, names, bound):
        if sum(model.fields[n].ghost_number for n in legs) != ghost_number:
            continue
        dimension = sum((model.fields[n].dimension for n in legs), Fraction(0))
        counts = structure_counts(model, legs, floor(bound - dimension), physical=physical)
        if not any(counts):
            continue
        out.append(LegClass(legs, ghost_number, dimension, bound, counts, kind_signature(model, legs),
                            sector(model, legs)))
    order = model.fields.order
    return sorted(out, key=lambda c: (len(c.fields), [order(n) for n in c.fields]))


def by_sector(classes: Iterable[LegClass]) -> dict[str, list[LegClass]]:
    """The classes grouped by :attr:`LegClass.key` (ghost number and sector); keys ordered by ghost number, then
    by the number of legs of their smallest class."""
    groups: dict[str, list[LegClass]] = {}
    for c in sorted(classes, key=lambda c: (c.ghost_number, len(c.fields))):
        groups.setdefault(c.key, []).append(c)
    return groups
