"""Operator bases: independent local monomials per class of external legs, generated from the field content.

The recipe (report 02_method sec:m-algorithm, 03b sec:Qbasis, sec:Nbasis, sec:evanescent) needs, at every order
``l``, the same bases:

* ``Q``: ghost number 0, dimension ``<= 4``, four-dimensional, Lorentz and rigidly invariant, hermitian with real
  coefficients (the restoring counterterm ``Delta S^(l) = sum r_alpha Q_alpha``; :func:`counterterm_basis`);
* ``N``: ghost number 1, dimension ``<= 5``, four-dimensional (the breaking ``X_l`` and ``b Q``;
  :func:`breaking_basis`); ghost number 2 (the consistency condition ``b X_l = 0``) and ghost number -1 (``Q``
  directions that are ``b``-exact) come from the same generator (:func:`operator_basis`);
* the d-dimensional singular basis: the continuation of ``Q`` plus the evanescent operators that can carry poles
  (:func:`singular_basis`).

Nothing is listed by hand: classes, group invariants, Lorentz / Dirac structures and derivative placements are
derived from the model, and every selection is an exact rank of plane-wave fingerprints
(:mod:`hepkit_renorm.fingerprint`, :mod:`hepkit_renorm.linalg`), which takes care of integration by parts,
momentum conservation, Schouten, Chisholm and Fierz identities and Grassmann symmetry at once.

Generation (overcomplete, per class and momentum degree)
--------------------------------------------------------
For a leg class ``(Phi_1, .., Phi_n)`` of :func:`hepkit_renorm.classes.leg_classes` and each momentum degree ``k``
(``k`` partial derivatives) the candidates are products of

* a derivative placement: ``k`` derivatives distributed over the occurrences, none on one occurrence of the
  *integration-by-parts field* (the field of largest ghost number, first in canonical order: complete by
  integration by parts, and the ghost of ``N`` stays undifferentiated as in report tab:Nbasis); identical fields
  carry non-increasing numbers of derivatives (they are symmetrised by the monomial anyway);
* a Lorentz / Dirac structure on the Lorentz labels of the field slots and derivatives: every spinor line
  (a row spinor ``Psibar``, ``Rbar`` paired with a column spinor ``Psi``, ``R``) carries the projectors of its
  physical chiralities and an ordered product of up to four ``gammabar`` (with up to two internal contracted
  pairs between different lines); the remaining labels are contracted by ``gbar``, or four of them by ``eps`` and
  the rest by ``gbar`` (four-dimensional: barred heads only);
* one invariant tensor per gauge factor: every adjoint slot is written as ``T^a_{ij}`` and the fundamental-type
  indices are paired with the antifundamental-type ones by ``delta``, with ``epsf`` on the excess (the first
  fundamental theorem of SU(N): products of traces and open chains of generators span every invariant); U(1)
  factors only require neutrality (guaranteed by the class).

The physical chiral projectors remove the sterile partners; the auxiliary field and the antighost are not legs
(classes). A gauge source enters as the shifted source ``rhotilde = rho + d omegabar`` (report eq:rhotilde;
:func:`shift_gauge_sources`), so that ``b`` acting on a basis element has no ``b``-field terms; fingerprints are
taken on the classes without antighost, where ``rhotilde`` and ``rho`` coincide.

Reduction (exact)
-----------------
Candidates are fingerprinted (dimension 4) on the samples of their class (and of the conjugate class for the
hermitian bases) and kept greedily when they raise the rank (:func:`hepkit_renorm.linalg.rank`). The number of
independent structures of each class and degree is known exactly from :func:`hepkit_renorm.classes.
structure_counts` (characters); generation stops when it is reached and fails loudly if the candidates are
exhausted first, so completeness is certified class by class. For the hermitian basis a candidate ``M`` gives
``M + M^dagger`` and ``i (M - M^dagger)``; over a pair of conjugate classes ``C != Cbar`` the target is
``n(C) + n(Cbar)``, over a self-conjugate class ``n(C)`` (a hermitian family that is independent over the complex
numbers is independent over the reals, and conversely).

Evanescent operators (:func:`singular_basis`)
---------------------------------------------
The continuation of a ``Q`` element to d dimensions (report sec:evanescent, 'Continuation') keeps the barred heads
of monomials with spinor lines (barred fermion vertex, ``psibar i dslash psi = psibar i dbarslash psi`` for a
physical field) and replaces ``gbar`` by the unsplit ``g`` in bosonic monomials. The evanescent candidates that can
carry one-loop poles are derived, not listed:

* *hatted derivatives*: the generator above with every pair of derivative labels contracted by ``ghat`` instead
  of ``gbar`` and ``gammahat`` on derivative labels in spinor lines (external fields keep barred indices: they
  couple to fermion loops only through the barred vertex, and bosonic loops are d-dimensional covariant);
* *barred minus continued invariants*: ``P(K) - K`` for every kinetic invariant ``K`` of the classical action
  (per gauge factor ``-1/4 F^2``, per scalar multiplet ``|D phi|^2``), i.e. the difference between the invariant
  built from barred objects (a closed fermion loop) and its d-dimensional continuation.

Both sets are reduced by the rank of fingerprints in ``D = 6`` and ``8`` jointly (they vanish in four dimensions).
For the one-family Standard Model this gives report eq:evanescent (``E_A``, ``E_H``, ``Etilde_F``, ``Etilde_H``).
Beyond one loop the complete hatted basis (hatted field components included) is generated with
``hatted_fields=True``.

Conventions: those of :mod:`hepkit_renorm.functional` (Grassmann order, vertex functions) and
:mod:`hepkit_renorm.tensors` (heads); hermitian conjugation is :meth:`Density.conjugate`. The generation is
deterministic: candidate order, pivots and names depend only on the model.
"""

from __future__ import annotations

import itertools
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass, field
from fractions import Fraction

import numpy as np
from symbolica import Expression, S

from . import linalg
from . import tensors as tn
from .action import field_strength, scalar_kinetic
from .classes import EXCLUDED_KINDS, LegClass, conjugate_fields, dimension_bound, leg_classes
from .exact import PrimeField
from .fingerprint import Fingerprinter, FingerprintMatrix, Samples
from .functional import Density, Derivation, FieldImage, LocalFunctional, Monomial, Occurrence
from .groups import GaugeFactor
from .spec import Model

__all__ = [
    "BasisElement",
    "IncompleteBasisError",
    "OperatorBasis",
    "SingularBasis",
    "breaking_basis",
    "candidates",
    "continuation",
    "counterterm_basis",
    "operator_basis",
    "shift_gauge_sources",
    "singular_basis",
]

_MAX_GAMMAS = 4
_MARGIN = 4          # extra samples per class beyond the number of structures
_PREFIX = {-1: "Y", 0: "Q", 1: "N", 2: "K"}
_WA, _WB = S("a_"), S("b_")


class IncompleteBasisError(RuntimeError):
    """The candidates of a class were exhausted before the exact number of structures was reached."""


# ---------------------------------------------------------------------------------------------- basis objects
@dataclass(frozen=True)
class BasisElement:
    """One basis monomial.

    ``name``: deterministic, ``"<prefix><k>[<fields>]"`` (prefix ``Q``, ``N``, ``K``, ``Y`` for ghost number 0, 1,
    2, -1; ``E`` for evanescent operators); ``functional``: the local functional (gauge sources shifted to
    ``rhotilde``); ``leg_class``: the class it was generated in; ``classes``: the classes it lives in (with the
    conjugate class for hermitian elements); ``degree``: number of derivatives (``-1`` for an operator with parts
    of several degrees, e.g. ``P(K) - K``); ``kind``: ``"monomial"``,
    ``"hermitian+"`` (``M + M^dagger``), ``"hermitian-"`` (``i (M - M^dagger)``), ``"continued"``,
    ``"evanescent"``; ``origin``: the generating monomial (or the operator it was derived from) in readable form.
    """

    name: str
    functional: LocalFunctional
    leg_class: tuple[str, ...]
    classes: tuple[tuple[str, ...], ...]
    degree: int
    kind: str
    origin: str


@dataclass
class OperatorBasis:
    """An independent set of local functionals of one ghost number, with the fingerprint layout that certifies
    its independence.

    ``classes``: the leg classes (:class:`~hepkit_renorm.classes.LegClass`) it covers; ``layout``: ``(class,
    samples)`` of the rows used to select the elements (dimension ``dimension``); ``hermitian``: real basis of
    hermitian monomials (``True``) or a complex basis (``False``).
    """

    model: Model
    ghost_number: int
    bound: Fraction
    hermitian: bool
    elements: list[BasisElement]
    classes: list[LegClass]
    layout: tuple[tuple[tuple[str, ...], int], ...]
    dimension: int = 4

    def __len__(self) -> int:
        return len(self.elements)

    def __iter__(self) -> Iterator[BasisElement]:
        return iter(self.elements)

    def __getitem__(self, k: int) -> BasisElement:
        return self.elements[k]

    @property
    def names(self) -> list[str]:
        return [e.name for e in self.elements]

    @property
    def functionals(self) -> list[LocalFunctional]:
        return [e.functional for e in self.elements]

    def index(self, name: str) -> int:
        return self.names.index(name)

    def by_class(self) -> dict[tuple[str, ...], list[BasisElement]]:
        """Elements grouped by the class they were generated in."""
        out: dict[tuple[str, ...], list[BasisElement]] = {}
        for e in self.elements:
            out.setdefault(e.leg_class, []).append(e)
        return out

    def counts(self) -> dict[tuple[str, ...], int]:
        """Number of elements living in each class (conjugate classes counted for hermitian elements)."""
        out: dict[tuple[str, ...], int] = {}
        for e in self.elements:
            for c in e.classes:
                out[c] = out.get(c, 0) + 1
        return out

    def matrix(self, fingerprinter: Fingerprinter, functionals: Sequence[LocalFunctional] | None = None,
               salt: int = 0) -> FingerprintMatrix:
        """Fingerprints of the basis (or of other ``functionals``) on the rows of :attr:`layout`."""
        funcs = self.functionals if functionals is None else list(functionals)
        return fingerprinter.matrix(funcs, layout=self.layout, salt=salt)

    def coordinates(self, functionals: Sequence[LocalFunctional], fingerprinter: Fingerprinter | None = None,
                    ) -> np.ndarray:
        """Coordinates ``c`` (shape ``(len(functionals), len(self))``, over the fingerprinter's prime field) with
        ``F_k = sum_alpha c[k, alpha] Q_alpha``, compared on every class of the basis. Classes with an antighost or
        an auxiliary field are not compared (gauge sources enter as ``rhotilde``). Raises
        :class:`~hepkit_renorm.linalg.InconsistentSystemError` if a functional is not in the span, also when it
        has a component in a class that the basis does not cover."""
        fp = fingerprinter or Fingerprinter(self.model, dimension=self.dimension)
        base = self.matrix(fp).values
        rhs = self.matrix(fp, functionals).values
        known = {c for c, _ in self.layout}
        for k, f in enumerate(functionals):
            for cls in external_classes(self.model, f):
                if cls not in known and not fp.vanishes(LocalFunctional.integrate(_restrict(f, cls))):
                    raise linalg.InconsistentSystemError([k])
        return linalg.solve(base, rhs, fp.field).T

    def contains(self, functionals: Sequence[LocalFunctional], fingerprinter: Fingerprinter | None = None) -> bool:
        """True if every functional lies in the span of the basis (:meth:`coordinates` succeeds)."""
        try:
            self.coordinates(functionals, fingerprinter)
        except linalg.InconsistentSystemError:
            return False
        return True


def external_classes(model: Model, functional: LocalFunctional) -> list[tuple[str, ...]]:
    """The leg classes of a functional without antighost and auxiliary legs (the classes that fingerprints
    compare; module docstring)."""
    return [c for c in functional.leg_classes() if not any(model.fields[n].kind in EXCLUDED_KINDS for n in c)]


def _restrict(functional: LocalFunctional, cls: tuple[str, ...]) -> Density:
    key = sorted(cls)
    return Density.from_monomials(functional.fields, [m for m in functional.monomials() if sorted(m.fields()) == key])


# ---------------------------------------------------------------------------------------------- shifted sources
def shift_gauge_sources(model: Model, functional: LocalFunctional) -> LocalFunctional:
    """``F[rho -> rhotilde]`` with ``rhotilde^mu_a = rho^mu_a + d^mu omegabar_a`` (report eq:rhotilde) for the
    gauge source of every factor that has an antighost: ``exp(delta) F`` with the even derivation ``delta rho =
    d omegabar`` (polynomial: the series terminates)."""
    images = {}
    for name in model.factors:
        rho = model.fields[name].source
        agh = f"agh_{name}"
        if rho is None or agh not in model.fields:
            continue
        labels = model.slot_labels(rho)
        mu, adj = labels[0], labels[1:]
        images[rho] = FieldImage(labels, Density.monomial(model.fields, 1, Occurrence(agh, tuple(adj), (mu,))))
    delta = Derivation(0, images)
    total, term, k = functional, functional, 1
    while True:
        term = LocalFunctional.integrate(delta(term).scale(Expression.num(1) / k))
        if term.is_zero():
            return total
        total = LocalFunctional.integrate(total + term)
        k += 1


# ---------------------------------------------------------------------------------------------- group tensors
def _set_partitions(items: Sequence, block: int) -> Iterator[list[tuple]]:
    """Partitions of ``items`` into blocks of size ``block`` (each block in increasing position order)."""
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for others in itertools.combinations(range(len(rest)), block - 1):
        blk = (first,) + tuple(rest[i] for i in others)
        remaining = [x for i, x in enumerate(rest) if i not in others]
        for tail in _set_partitions(remaining, block):
            yield [blk] + tail


def _group_tensors(factor: GaugeFactor, slots: Sequence[tuple[Expression, str]]) -> list[Expression]:
    """Spanning set of the invariant tensors of one SU(N) factor on the given slot labels (module docstring):
    adjoint slots as ``T^a_{ij}``, fundamental-type indices paired with antifundamental-type ones by ``delta``,
    ``epsf`` on the excess. Single traces ``T^a_{ii}`` (zero) are skipped."""
    if factor.abelian:
        return [Expression.num(1)]
    fund = factor.rep_space("F")
    heads = Expression.num(1)
    up: list[tuple[Expression, int | None]] = []      # fundamental-type index, owning generator
    down: list[tuple[Expression, int | None]] = []    # antifundamental-type index, owning generator
    for k, (lab, rep) in enumerate(slots):
        if rep == "A":
            i, j = tn.fresh(fund), tn.fresh(fund)
            heads = heads * tn.T(lab, i, j)
            up.append((i, k))
            down.append((j, k))
        elif rep == "F":
            up.append((lab, None))
        elif rep == "Fbar":
            down.append((lab, None))
        else:
            raise ValueError(f"unsupported representation {rep!r} of {factor.name}")
    excess = len(up) - len(down)
    n = factor.n
    if excess % n:
        return []
    side, other = (up, down) if excess >= 0 else (down, up)
    n_eps = abs(excess) // n
    out = []
    for chosen in itertools.combinations(range(len(side)), n_eps * n):
        eps_part = [side[i] for i in chosen]
        rest = [x for i, x in enumerate(side) if i not in chosen]
        for blocks in _set_partitions(eps_part, n):
            eps = Expression.num(1)
            for blk in blocks:
                eps = eps * tn.epsf(*[x[0] for x in blk])
            for perm in itertools.permutations(other):
                if any(a[1] is not None and a[1] == b[1] for a, b in zip(rest, perm)):
                    continue                                    # Tr T^a = 0
                deltas = Expression.num(1)
                for a, b in zip(rest, perm):
                    deltas = deltas * tn.delta(a[0], b[0])
                out.append(tn.absorb(heads * eps * deltas))
    return out


# ---------------------------------------------------------------------------------------------- Lorentz structures
@dataclass(frozen=True)
class _Line:
    row: Expression           # spinor label of the row spinor
    col: Expression           # spinor label of the column spinor
    row_projector: Expression  # head applied to the right of the row spinor
    col_projector: Expression  # head applied to the left of the column spinor


def _spinor_side(model: Model, name: str) -> tuple[str, Expression] | None:
    """``("row" | "col", projector head)`` of the spinor slot of a field or source (physical chirality)."""
    fi = model.fields[name]
    if not any(s.space.kind == "spinor" for s in fi.slots):
        return None
    base = model.fields[fi.source_of] if fi.is_source else fi
    if base.chirality not in ("L", "R"):
        raise NotImplementedError(f"fermion {base.name!r} without a definite chirality")
    chir = base.chirality if base.kind == "fermion" else {"L": "R", "R": "L"}[base.chirality]
    head = tn.PL if chir == "L" else tn.PR
    is_row = (base.kind == "fermion") if fi.is_source else (fi.kind == "fermion_bar")
    return ("row" if is_row else "col"), head


def _matchings(labels: Sequence) -> Iterator[list[tuple]]:
    if not labels:
        yield []
        return
    first, rest = labels[0], labels[1:]
    for k in range(len(rest)):
        for tail in _matchings(rest[:k] + rest[k + 1:]):
            yield [(first, rest[k])] + tail


def _contractions(rest: Sequence[tuple[Expression, str]], hats: bool) -> Iterator[Expression]:
    """``gbar`` matchings of ``rest``, and ``eps`` on four of them times matchings of the others; with ``hats``
    a pair of two derivative labels may be contracted by ``ghat`` instead."""
    if len(rest) % 2:
        return
    options = [(None, list(rest))]
    if len(rest) >= 4:
        for quad in itertools.combinations(range(len(rest)), 4):
            options.append((tn.eps(*[rest[i][0] for i in quad]), [x for i, x in enumerate(rest) if i not in quad]))
    for eps, others in options:
        for match in _matchings(others):
            choices = []
            for a, b in match:
                metrics = [tn.gbar(a[0], b[0])]
                if hats and a[1] == "der" and b[1] == "der":
                    metrics.append(tn.ghat(a[0], b[0]))
                choices.append(metrics)
            for picked in itertools.product(*choices):
                c = Expression.num(1) if eps is None else eps
                for x in picked:
                    c = c * x
                yield c


def _strings(pool: Sequence[tuple[Expression, str]], n_lines: int) -> Iterator[tuple[tuple, ...]]:
    """Ordered disjoint tuples of pool entries, one per line, each of length ``<= _MAX_GAMMAS``."""
    if n_lines == 0:
        yield ()
        return
    for length in range(min(len(pool), _MAX_GAMMAS) + 1):
        for chosen in itertools.permutations(range(len(pool)), length):
            rest = [x for i, x in enumerate(pool) if i not in chosen]
            for tail in _strings(rest, n_lines - 1):
                yield (tuple(pool[i] for i in chosen),) + tail


def _lorentz_structures(pool: Sequence[tuple[Expression, str]], lines: Sequence[_Line], hats: bool,
                        ) -> Iterator[Expression]:
    """Lorentz / Dirac structures on the pool labels (module docstring)."""
    for n_internal in range(3 if len(lines) >= 2 else 1):
        internal = []
        for _ in range(n_internal):
            x, y = tn.fresh(tn.LORENTZ), tn.fresh(tn.LORENTZ)
            internal.append((x, y))
        ext = list(pool) + [(x, "int") for pair in internal for x in pair]
        for strings in _strings(ext, len(lines)):
            used = {x[0].to_canonical_string(): k for k, s in enumerate(strings) for x in s}
            ok = all(x.to_canonical_string() in used and y.to_canonical_string() in used
                     and used[x.to_canonical_string()] != used[y.to_canonical_string()] for x, y in internal)
            if not ok:
                continue
            for dirac in _dirac_products(strings, lines, hats):
                coeff = dirac
                for x, y in internal:
                    coeff = coeff * tn.gbar(x, y)
                rest = [e for e in pool if e[0].to_canonical_string() not in used]
                for c in _contractions(rest, hats):
                    yield coeff * c


def _dirac_products(strings: Sequence[tuple], lines: Sequence[_Line], hats: bool) -> Iterator[Expression]:
    per_line = []
    for s, line in zip(strings, lines):
        variants = []
        heads = [[tn.gammabar] + ([tn.gammahat] if hats and kind == "der" else []) for _, kind in s]
        for picked in itertools.product(*heads):
            u = [tn.fresh(tn.SPINOR) for _ in range(len(s) + 1)]
            e = line.row_projector(line.row, u[0])
            for k, (head, (lab, _)) in enumerate(zip(picked, s)):
                e = e * head(lab, u[k], u[k + 1])
            variants.append(e * line.col_projector(u[-1], line.col))
        per_line.append(variants)
    for picked in itertools.product(*per_line):
        out = Expression.num(1)
        for x in picked:
            out = out * x
        yield out


# ---------------------------------------------------------------------------------------------- candidates
def _ibp_field(model: Model, legs: Sequence[str]) -> str:
    """The field kept free of derivatives: largest ghost number, first in canonical order."""
    return max(legs, key=lambda n: (model.fields[n].ghost_number, -model.fields.order(n)))


def _derivative_counts(legs: Sequence[str], degree: int, free: str) -> Iterator[tuple[int, ...]]:
    """Derivatives per occurrence: non-increasing within identical fields, one occurrence of ``free`` without."""
    groups: list[tuple[str, int]] = []
    for name in legs:
        if groups and groups[-1][0] == name:
            groups[-1] = (name, groups[-1][1] + 1)
        else:
            groups.append((name, 1))

    def tuples(m: int, total: int, largest: int) -> Iterator[tuple[int, ...]]:
        if m == 0:
            if total == 0:
                yield ()
            return
        for first in range(min(total, largest), -1, -1):
            for rest in tuples(m - 1, total - first, first):
                yield (first,) + rest

    def rec(k: int, budget: int) -> Iterator[tuple[int, ...]]:
        if k == len(groups):
            if budget == 0:
                yield ()
            return
        name, m = groups[k]
        for used in range(budget + 1):
            for t in tuples(m, used, used):
                if name == free and t[-1] != 0:
                    continue
                for rest in rec(k + 1, budget - used):
                    yield t + rest

    yield from rec(0, degree)


def candidates(model: Model, legs: Sequence[str], degree: int, *, hats: bool = False,
               hatted_fields: bool = False) -> Iterator[Monomial]:
    """Overcomplete candidate monomials of the class ``legs`` (canonical order) with ``degree`` derivatives
    (module docstring), four-dimensional unless ``hats``: then derivative pairs may be contracted by ``ghat`` and
    derivative labels may carry ``gammahat``; ``hatted_fields`` also allows ``ghat`` on field labels (complete
    d-dimensional basis). Deterministic order."""
    if model.generations != 1:
        raise NotImplementedError("operator bases with flavour indices (generations > 1)")
    legs = list(legs)
    free = _ibp_field(model, legs)
    for counts in _derivative_counts(legs, degree, free):
        occs, pool, rows, cols = [], [], [], []
        group_slots: dict[str, list[tuple[Expression, str]]] = {f: [] for f in model.factors}
        for name, nder in zip(legs, counts):
            fi = model.fields[name]
            slots = []
            for s in fi.slots:
                lab = tn.fresh(s.space)
                slots.append(lab)
                if s.space.kind == "lorentz":
                    pool.append((lab, "der" if hatted_fields else "field"))
                elif s.space.kind == "spinor":
                    side, head = _spinor_side(model, name)
                    (rows if side == "row" else cols).append((lab, head))
                elif s.factor is not None:
                    group_slots[s.factor].append((lab, s.rep))
            ders = tuple(tn.fresh(tn.LORENTZ) for _ in range(nder))
            pool += [(mu, "der") for mu in ders]
            occs.append(Occurrence(name, tuple(slots), ders))
        if len(rows) != len(cols):
            continue
        groups = [_group_tensors(model.factors[f], group_slots[f]) for f in model.factors
                  if group_slots[f] or not model.factors[f].abelian]
        if any(not g for g in groups):
            continue
        for perm in itertools.permutations(cols):
            lines = [_Line(r[0], c[0], r[1], c[1]) for r, c in zip(rows, perm)]
            for lorentz in _lorentz_structures(pool, lines, hats):
                if hats and not tn.is_evanescent(lorentz):
                    continue
                for combo in itertools.product(*groups):
                    coeff = lorentz
                    for g in combo:
                        coeff = coeff * g
                    yield Monomial(coeff, tuple(occs))


# ---------------------------------------------------------------------------------------------- reduction
class _Reducer:
    """Greedy exact selection of independent columns on fixed rows: ``(fingerprinter, samples)`` per row block;
    a column is kept when it raises the rank (:func:`linalg.rank`)."""

    def __init__(self, blocks: Sequence[tuple[Fingerprinter, Samples]]) -> None:
        self.blocks = list(blocks)
        self.columns: list[np.ndarray] = []
        self.field = self.blocks[0][0].field

    def vector(self, functional: LocalFunctional) -> np.ndarray:
        return np.concatenate([fp.functional(functional, s) for fp, s in self.blocks])

    def offer(self, functional: LocalFunctional) -> bool:
        v = self.vector(functional)
        if not np.any(v):
            return False
        trial = np.stack(self.columns + [v], axis=1)
        if linalg.rank(trial, self.field) == len(self.columns) + 1:
            self.columns.append(v)
            return True
        return False


def _hermitian_pair(f: LocalFunctional) -> list[tuple[str, LocalFunctional]]:
    """``(M + M^dagger)/2`` and ``i (M - M^dagger)/2`` (``M`` itself for a hermitian ``M``)."""
    conj = f.conjugate()
    half = Expression.num(1) / 2
    plus = LocalFunctional.integrate((f + conj).scale(half))
    minus = LocalFunctional.integrate((f - conj).scale(tn.I * half))
    return [("hermitian+", plus), ("hermitian-", minus)]


def operator_basis(model: Model, ghost_number: int, *, bound: Fraction | int | None = None,
                   hermitian: bool | None = None, classes: Iterable[Sequence[str]] | None = None,
                   fingerprinter: Fingerprinter | None = None) -> OperatorBasis:
    """The basis of local monomials of ghost number ``ghost_number`` and dimension ``<= bound`` (default
    :func:`~hepkit_renorm.classes.dimension_bound`), four-dimensional, reduced exactly (module docstring).

    ``hermitian`` (default: ``ghost_number == 0``) gives hermitian monomials with real coefficients; otherwise a
    complex basis. ``classes`` restricts the leg classes (a hermitian basis adds the conjugate classes). The
    fingerprints use ``fingerprinter`` (default: dimension 4, default prime and seed). Raises
    :class:`IncompleteBasisError` if a class cannot be completed (a generator bug)."""
    hermitian = ghost_number == 0 if hermitian is None else hermitian
    bound = dimension_bound(model, ghost_number) if bound is None else Fraction(bound)
    fp = fingerprinter or Fingerprinter(model)
    all_classes = {c.fields: c for c in leg_classes(model, ghost_number, bound=bound)}
    if classes is None:
        wanted = list(all_classes)
    else:
        wanted = [tuple(sorted(c, key=model.fields.order)) for c in classes]
        missing = [c for c in wanted if c not in all_classes]
        if missing:
            raise ValueError(f"not leg classes of ghost number {ghost_number}: {missing}")
    order = model.fields.order
    units: list[tuple[tuple[str, ...], ...]] = []
    seen: set[tuple[str, ...]] = set()
    for c in sorted(wanted, key=lambda c: (len(c), [order(x) for x in c])):
        if c in seen:
            continue
        unit = (c,)
        if hermitian:
            cbar = conjugate_fields(model, c)
            unit = tuple(sorted({c, cbar}, key=lambda x: [order(n) for n in x]))
        seen.update(unit)
        units.append(unit)
    elements: list[BasisElement] = []
    layout: list[tuple[tuple[str, ...], int]] = []
    prefix = _PREFIX.get(ghost_number, f"G{ghost_number}_")
    for unit in units:
        lead = unit[0]
        info = [all_classes[c] for c in unit]
        layout += [(c, info[0].n_structures * len(unit) + _MARGIN) for c in unit]
        for degree in range(len(info[0].structures)):
            target = sum(c.structures[degree] for c in info) if hermitian else info[0].structures[degree]
            if target == 0:
                continue
            n = target + _MARGIN
            reducer = _Reducer([(fp, fp.samples(c, n, salt=degree)) for c in unit])
            found: list[BasisElement] = []
            seen: set[str] = set()
            for mono in candidates(model, lead, degree):
                f = LocalFunctional.integrate(Density.from_monomials(model.fields, [mono]))
                key = f.expr.to_canonical_string()
                if f.is_zero() or key in seen:
                    continue
                seen.update((key, (-f).expr.to_canonical_string()))
                options = _hermitian_pair(f) if hermitian else [("monomial", f)]
                for kind, g in options:
                    if not g.is_zero() and reducer.offer(g):
                        found.append(BasisElement("", g, lead, unit, degree, kind, str(f.expr)))
                if len(found) == target:
                    break
            if len(found) != target:
                raise IncompleteBasisError(f"class {lead}, degree {degree}: {len(found)} of {target} structures")
            elements += found
    named = []
    for k, e in enumerate(elements, start=1):
        func = shift_gauge_sources(model, e.functional)
        named.append(BasisElement(f"{prefix}{k}[{' '.join(e.leg_class)}]", func, e.leg_class, e.classes, e.degree,
                                  e.kind, e.origin))
    return OperatorBasis(model, ghost_number, bound, hermitian, named,
                         [all_classes[c] for u in units for c in u], tuple(layout), fp.dimension)


def counterterm_basis(model: Model, **kwargs) -> OperatorBasis:
    """``Q``: ghost number 0, dimension ``<= 4``, hermitian with real coefficients (report sec:Qbasis)."""
    return operator_basis(model, 0, **kwargs)


def breaking_basis(model: Model, **kwargs) -> OperatorBasis:
    """``N``: ghost number 1, dimension ``<= 5``, complex (report sec:Nbasis; complete in every class)."""
    kwargs.setdefault("hermitian", False)
    return operator_basis(model, 1, **kwargs)


# ---------------------------------------------------------------------------------------------- singular basis
def continuation(element: BasisElement) -> BasisElement:
    """The d-dimensional continuation of a four-dimensional basis element (report sec:evanescent): barred heads
    kept in monomials with spinor lines, ``gbar -> g`` in bosonic monomials."""
    f = element.functional
    fields = f.fields
    total = Density(fields)
    for m in f.monomials():
        spinor = any(s.space.kind == "spinor" for o in m.occurrences for s in fields[o.field].slots)
        coeff = m.coefficient if spinor else m.coefficient.replace(tn.gbar(_WA, _WB), tn.g(_WA, _WB))
        total = total + Density.from_monomials(fields, [Monomial(coeff, m.occurrences)])
    return BasisElement(element.name, LocalFunctional.integrate(total), element.leg_class, element.classes,
                        element.degree, "continued", element.name)


@dataclass
class SingularBasis:
    """The d-dimensional basis of singular counterterms: ``continued`` (the continuation of ``Q``) and
    ``evanescent`` (operators that vanish in four dimensions and can carry poles; module docstring), with the
    fingerprint layout in ``D = 6, 8`` that certifies the independence of the evanescent ones."""

    continued: list[BasisElement]
    evanescent: list[BasisElement]
    layout: tuple[tuple[tuple[str, ...], int], ...] = field(default=())

    @property
    def elements(self) -> list[BasisElement]:
        return self.continued + self.evanescent


def kinetic_invariants(model: Model) -> list[tuple[str, LocalFunctional]]:
    """The kinetic invariants of the classical action: ``-1/4 F^2`` per gauge factor and the scalar kinetic term
    per scalar multiplet (d-dimensional, as in the action)."""
    out = []
    for name, fac in model.factors.items():
        mu, nu = tn.fresh(tn.LORENTZ), tn.fresh(tn.LORENTZ)
        a = None if fac.abelian else tn.fresh(fac.adjoint_space())
        fs = field_strength(model, name, mu, nu, a)
        out.append((f"F2[{name}]", LocalFunctional.integrate((fs * fs).scale(Expression.num(-1) / 4))))
    kin = scalar_kinetic(model)
    for name in model.scalars:
        members = {n for n, fi in model.fields.items() if fi.multiplet == name and not fi.is_source}
        part = Density.from_monomials(model.fields, [m for m in kin.monomials() if members & set(m.fields())])
        if not part.is_zero():
            out.append((f"DS2[{name}]", LocalFunctional.integrate(part)))
    return out


def singular_basis(basis: OperatorBasis, *, hatted_fields: bool = False, dimensions: Sequence[int] = (6, 8),
                   field_: PrimeField | None = None) -> SingularBasis:
    """The d-dimensional singular basis built on a four-dimensional ghost-number-0 basis (module docstring).

    Evanescent candidates: hatted-derivative monomials of every class of ``basis`` (``hatted_fields``: also
    hatted field components), then ``P(K) - K`` for the kinetic invariants; reduced by exact rank of fingerprints
    in ``dimensions`` jointly; names ``E<k>[<class>]``. For a hermitian ``basis`` the hatted candidates enter as
    hermitian combinations (``P(K) - K`` is hermitian)."""
    model = basis.model
    simp = model.simplifier()
    fps = [Fingerprinter(model, field_, D) for D in dimensions]
    conts = [continuation(e) for e in basis.elements]
    cands: list[tuple[tuple[str, ...], int, str, LocalFunctional]] = []
    for cls in basis.classes:
        for degree in range(len(cls.structures)):
            for mono in candidates(model, cls.fields, degree, hats=True, hatted_fields=hatted_fields):
                f = LocalFunctional.integrate(Density.from_monomials(model.fields, [mono]))
                if f.is_zero():
                    continue
                for _, g in (_hermitian_pair(f) if basis.hermitian else [("monomial", f)]):
                    if not g.is_zero():
                        cands.append((cls.fields, degree, str(f.expr), g))
    for name, k in kinetic_invariants(model):
        diff = LocalFunctional.integrate(k.physical(simp) - k)
        if not diff.is_zero():
            lead = min(external_classes(model, diff), key=lambda c: (len(c), [model.fields.order(x) for x in c]))
            cands.append((lead, -1, f"P({name}) - {name}", diff))
    classes = sorted({c for cand in cands for c in external_classes(model, cand[3])},
                     key=lambda c: [model.fields.order(x) for x in c])
    per_class = {c: sum(1 for cand in cands if c in external_classes(model, cand[3])) + _MARGIN for c in classes}
    layout = tuple((c, per_class[c]) for c in classes)
    reducer = _Reducer([(fp, fp.samples(c, n, salt=7)) for fp in fps for c, n in layout])
    found = []
    for cls, degree, origin, f in cands:
        if reducer.offer(f):
            found.append(BasisElement(f"E{len(found) + 1}[{' '.join(cls)}]", f, cls,
                                      tuple(external_classes(model, f)), degree, "evanescent", origin))
    return SingularBasis(conts, found, layout)
