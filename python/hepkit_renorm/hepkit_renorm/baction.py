"""The action of the linearised Slavnov-Taylor operator on the operator bases.

At order ``l`` (report 02_method sec:m-algorithm, step 9) the restoring counterterm solves ``b Delta S^(l) = -X_l``
in the four-dimensional basis ``Q`` (:func:`hepkit_renorm.bases.counterterm_basis`), with ``X_l`` in the
ghost-number-1 basis ``N`` (:func:`hepkit_renorm.bases.breaking_basis`). This module provides the linear algebra
of ``b`` on these bases, exactly over a prime field (:mod:`hepkit_renorm.linalg`), for any model and order (the
bases and ``b`` do not depend on ``l``):

* :func:`b_closure`: the condition ``b X = 0`` on the ghost-number-1 basis, as fingerprint rows on the
  ghost-number-2 classes (no ghost-number-2 basis is needed);
* :func:`b_table`: ``b Q_alpha = sum_beta B[alpha, beta] N_beta`` (report sec:btable), the image (rank), the
  kernel (the invariant, *closed* directions, report sec:closed), with rational reconstruction of the table at
  rational parameter values (:meth:`BTable.rational`);
* :func:`invariant_directions`: a derived description of the kernel by generators that exist for every model
  (report sec:closed): ``b``-exact directions ``b Y`` for the ghost-number ``-1`` basis ``Y`` (field
  renormalisations), coupling derivatives ``p dS_0/dp`` and phase deformations of parameters that multiply
  non-hermitian monomials, charge deformations ``dS_0/d epsilon`` under ``q_U(1) -> q_U(1) + epsilon q`` for every
  global U(1) symmetry ``q`` of the classical action, and the gauge-invariant source-free operators (the kernel
  restricted to the monomials without ghosts and sources); the part of the kernel they do not span is reported;
* :func:`bd_singular`: ``b_d`` on the d-dimensional singular basis (which elements are ``b_d``-closed), and
  :func:`breaking_by_factor`, :func:`coupling_direction`, :func:`ghost_direction` for the statement that ``b_d`` of
  the coupling and ghost directions of a gauge factor is that factor's part of ``Delta_1`` (report eq:bd-Oc).

``b`` is :attr:`hepkit_renorm.brst.BRST.b` (four-dimensional limit of ``b_d``); four-dimensional statements are
fingerprinted in ``D = 4`` (where ``d = 4`` and every component is barred), d-dimensional ones by
:func:`hepkit_renorm.fingerprint.is_zero`. Parameter values are those of the :class:`~hepkit_renorm.fingerprint.
Fingerprinter` (random by default, i.e. generic couplings).
"""

from __future__ import annotations

import dataclasses
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction

import numpy as np
from symbolica import Expression

from . import linalg
from . import tensors as tn
from .bases import BasisElement, OperatorBasis, SingularBasis, external_classes, operator_basis
from .brst import BRST
from .classes import conjugate_fields, leg_classes
from .exact import PrimeField
from .fingerprint import Fingerprinter, FingerprintMatrix, is_zero
from .functional import Density, LocalFunctional, Occurrence
from .spec import Model

__all__ = [
    "BTable",
    "Direction",
    "KernelAnalysis",
    "b_closure",
    "b_table",
    "bd_singular",
    "breaking_by_factor",
    "charge_deformation",
    "coupling_direction",
    "ghost_direction",
    "global_symmetries",
    "invariant_directions",
]


# ---------------------------------------------------------------------------------------------- the b-table
@dataclass
class BTable:
    """``b`` on a ghost-number-0 basis, expressed in a ghost-number-1 basis.

    ``images[alpha] = b Q_alpha``; ``matrix[alpha, beta] = B_{alpha beta}`` over ``fingerprinter.field`` at the
    fingerprinter's parameter values; ``image_values``: fingerprints of the images on the rows of
    ``target.layout`` (columns ``alpha``).
    """

    source: OperatorBasis
    target: OperatorBasis
    images: list[LocalFunctional]
    fingerprinter: Fingerprinter
    matrix: np.ndarray
    image_values: np.ndarray

    @property
    def field(self) -> PrimeField:
        return self.fingerprinter.field

    @property
    def rank(self) -> int:
        """Dimension of ``b(span Q)``."""
        return linalg.rank(self.image_values, self.field)

    @property
    def kernel(self) -> np.ndarray:
        """Basis of the closed directions as rows of ``Q`` coordinates (:func:`linalg.nullspace`)."""
        return linalg.nullspace(self.image_values, self.field)

    def restricted(self, columns: Sequence[int]) -> np.ndarray:
        """Kernel of ``b`` on the span of the given ``Q`` elements, embedded as rows of full ``Q`` coordinates."""
        sub = linalg.nullspace(self.image_values[:, list(columns)], self.field)
        out = np.zeros((sub.shape[0], len(self.source)), dtype=sub.dtype)
        out[:, list(columns)] = sub
        return out

    def is_closed(self, coordinates: np.ndarray) -> bool:
        """True if ``b (sum_alpha c_alpha Q_alpha) = 0`` (``coordinates`` a vector over the field)."""
        v = np.asarray(coordinates, dtype=object) % self.field.p
        prod = (np.asarray(self.image_values, dtype=object) @ v) % self.field.p
        return not np.any(prod)

    def rational(self, values: Mapping[str, object], max_primes: int = 8) -> tuple[np.ndarray, np.ndarray]:
        """The table at rational parameter values (all parameters given, e.g. unit couplings): real and imaginary
        parts as object arrays of :class:`~fractions.Fraction`, from the two embeddings of ``i``
        (:func:`linalg.split_gaussian`) reconstructed over several primes (:func:`linalg.reconstruct_rationals`)."""
        model = self.source.model

        def compute(field: PrimeField) -> np.ndarray:
            parts = []
            for sign in (1, -1):
                fp = Fingerprinter(model, field, self.fingerprinter.dimension, seed=self.fingerprinter.seed,
                                   values=values, i_sign=sign)
                parts.append(_solve_table(self.target, self.images, fp)[0])
            return np.stack(linalg.split_gaussian(parts[0], parts[1], field))

        both = linalg.reconstruct_rationals(compute, max_primes=max_primes)
        return both[0], both[1]


def _solve_table(target: OperatorBasis, images: Sequence[LocalFunctional], fp: Fingerprinter,
                 ) -> tuple[np.ndarray, np.ndarray]:
    base = target.matrix(fp).values
    rhs = target.matrix(fp, images).values
    known = {c for c, _ in target.layout}
    for k, f in enumerate(images):
        for cls in external_classes(target.model, f):
            if cls not in known:
                piece = LocalFunctional.integrate(Density.from_monomials(
                    f.fields, [m for m in f.monomials() if sorted(m.fields()) == sorted(cls)]))
                if not fp.vanishes(piece):
                    raise linalg.InconsistentSystemError([k])
    return linalg.solve(base, rhs, fp.field).T, rhs


def b_table(source: OperatorBasis, target: OperatorBasis, *, brst: BRST | None = None,
            fingerprinter: Fingerprinter | None = None) -> BTable:
    """``B`` with ``b Q_alpha = sum_beta B[alpha, beta] N_beta`` (report sec:btable) at the fingerprinter's parameter
    values (default: random, generic). Raises :class:`~hepkit_renorm.linalg.InconsistentSystemError` if some
    ``b Q_alpha`` is not in the span of ``target`` (impossible for a complete ``target``)."""
    if source.ghost_number + 1 != target.ghost_number:
        raise ValueError("the target basis must have ghost number one more than the source")
    brst = brst or BRST(source.model)
    fp = fingerprinter or Fingerprinter(source.model, dimension=4)
    images = [LocalFunctional.integrate(brst.b(e.functional)) for e in source.elements]
    matrix, values = _solve_table(target, images, fp)
    return BTable(source, target, images, fp, matrix, values)


def b_closure(basis: OperatorBasis, *, brst: BRST | None = None, fingerprinter: Fingerprinter | None = None,
              samples_per_class: int | None = None) -> tuple[list[LocalFunctional], FingerprintMatrix]:
    """``b`` of every element of a basis of ghost number ``g`` and the fingerprints of the images on every leg
    class of ghost number ``g + 1`` (rows). ``x`` satisfies ``b (sum_beta x_beta N_beta) = 0`` iff the matrix
    annihilates ``x``: the consistency condition ``b X_l = 0`` (report sec:Nbasis) stated without a basis of
    ghost number ``g + 1``. Default rows per class: its structure count plus four (the images in one class span
    at most that many directions, :func:`hepkit_renorm.classes.structure_counts`)."""
    model = basis.model
    brst = brst or BRST(model)
    fp = fingerprinter or Fingerprinter(model)
    images = [LocalFunctional.integrate(brst.b(e.functional)) for e in basis.elements]
    rows = leg_classes(model, basis.ghost_number + 1)
    layout = [(c.fields, samples_per_class or c.n_structures + 4) for c in rows]
    return images, fp.matrix(images, layout=layout)


# ---------------------------------------------------------------------------------------------- closed directions
@dataclass(frozen=True)
class Direction:
    """A derived generator of closed directions: ``functional`` (``None`` for directions known only by their
    coordinates), its ``Q`` ``coordinates`` (over the table's field) and whether ``b`` annihilates it."""

    name: str
    kind: str
    functional: LocalFunctional | None
    coordinates: np.ndarray
    closed: bool


@dataclass
class KernelAnalysis:
    """The kernel of ``b`` on ``span Q`` and its derived description.

    ``dimension``: ``dim ker b``; ``directions``: the derived generators (all kinds); ``rank_by_kind``: the rank
    of the closed generators accumulated in the order exact, coupling, phase, charge, invariant; ``spanned``:
    rank of all closed generators; ``remainder``: kernel vectors (rows) completing their span to the kernel.
    """

    dimension: int
    directions: list[Direction]
    rank_by_kind: dict[str, int]
    spanned: int
    remainder: np.ndarray


def global_symmetries(model: Model) -> list[dict[str, Fraction]]:
    """A basis of the U(1) charge assignments ``q`` (one charge per matter multiplet) that leave the Yukawa and
    potential terms invariant: the global U(1) symmetries of the classical action that act by phases on whole
    multiplets (they commute with the gauge group). Exact (rational null space, reconstructed)."""
    names = list(model.fermions) + [n for n, s in model.scalars.items() if s.complex]
    rows = []
    for y in model.yukawas:
        row = dict.fromkeys(names, 0)
        row[y.fermion_bar] -= 1
        row[y.fermion] += 1
        base = model.fields[y.scalar].multiplet
        row[base] += 1 if model.fields[y.scalar].kind == "scalar" else -1
        rows.append([row[n] for n in names])
    for term in model.potential:
        row = dict.fromkeys(names, 0)
        for a, b in term.bilinears:
            if a in row:
                row[a] -= 1
            if b in row:
                row[b] += 1
        rows.append([row[n] for n in names])

    def kernel(field: PrimeField) -> np.ndarray:
        m = np.asarray(rows, dtype=object).reshape(len(rows), len(names)) if rows else np.zeros((0, len(names)))
        return linalg.nullspace(linalg.as_field_matrix(m, field), field)

    basis = linalg.reconstruct_rationals(kernel)
    return [{n: Fraction(x) for n, x in zip(names, vec)} for vec in np.atleast_2d(basis)]


def charge_deformation(model: Model, factor: str, charges: Mapping[str, Fraction], brst: BRST | None = None,
                       ) -> LocalFunctional:
    """``d S_0 / d epsilon`` at ``epsilon = 0`` for the U(1) charges of ``factor`` shifted by ``epsilon q``
    (report eq:closed, ``cal(O)_q``: it gauges the global symmetry ``q`` with the gauge field of ``factor``).
    ``S_0`` is a polynomial of degree two in ``epsilon`` (scalar kinetic terms); the derivative is the exact
    three-point difference ``(-3 S(0) + 4 S(1) - S(2)) / 2``."""
    if not model.factors[factor].abelian:
        raise ValueError(f"{factor} is not abelian")
    actions = []
    for eps in (0, 1, 2):
        fermions = {n: dataclasses.replace(f, reps={**f.reps, factor: Fraction(f.reps.get(factor, 0))
                                                   + eps * Fraction(charges.get(n, 0))})
                    for n, f in model.fermions.items()}
        scalars = {n: dataclasses.replace(s, reps={**s.reps, factor: Fraction(s.reps.get(factor, 0))
                                                  + eps * Fraction(charges.get(n, 0))})
                   for n, s in model.scalars.items()}
        deformed = model if eps == 0 else dataclasses.replace(model, fermions=fermions, scalars=scalars)
        action = (brst if eps == 0 and brst is not None else BRST(deformed)).action
        actions.append(LocalFunctional(model.fields, action.expr))
    a0, a1, a2 = actions
    half = Expression.num(1) / 2
    return LocalFunctional.integrate((a0.scale(-3) + a1.scale(4) - a2).scale(half))


def _parameter_derivative(action: LocalFunctional, symbol: Expression) -> LocalFunctional:
    return LocalFunctional.integrate(action.map_coefficients(lambda c: (c.derivative(symbol) * symbol).expand()))


def _phase_part(model: Model, functional: LocalFunctional) -> LocalFunctional:
    """``i sum_{C < Cbar} (F_C - F_Cbar)`` over the pairs of conjugate classes of ``F`` (zero if every class of
    ``F`` is self-conjugate)."""
    order = model.fields.order
    total = Density(model.fields)
    for cls in functional.leg_classes():
        cbar = conjugate_fields(model, cls)
        if cbar == cls:
            continue
        sign = 1 if [order(n) for n in cls] < [order(n) for n in cbar] else -1
        part = Density.from_monomials(model.fields, [m for m in functional.monomials()
                                                     if sorted(m.fields()) == sorted(cls)])
        total = total + part.scale(tn.I * sign)
    return LocalFunctional.integrate(total)


def _try_coordinates(basis: OperatorBasis, functional: LocalFunctional, fp: Fingerprinter) -> np.ndarray | None:
    try:
        return basis.coordinates([functional], fp)[0]
    except linalg.InconsistentSystemError:
        return None


def invariant_directions(table: BTable, brst: BRST | None = None) -> KernelAnalysis:
    """The kernel of ``b`` on ``span Q`` and its derived generators (module docstring). Generators that do not lie
    in ``span Q`` (e.g. the gauge-fixing parameter, whose derivative lives in auxiliary-field classes) are
    skipped; generators outside the kernel are kept with ``closed = False``."""
    Q, fp = table.source, table.fingerprinter
    model = Q.model
    brst = brst or BRST(model)
    dirs: list[Direction] = []

    def add(name: str, kind: str, functional: LocalFunctional | None, coords: np.ndarray | None) -> None:
        if coords is None or not np.any(coords):
            return
        dirs.append(Direction(name, kind, functional, coords, table.is_closed(coords)))

    exact = operator_basis(model, Q.ghost_number - 1, hermitian=False, fingerprinter=Fingerprinter(model))
    for y in exact.elements:
        image = LocalFunctional.integrate(brst.b(y.functional))
        add(f"b {y.name}", "exact", image, _try_coordinates(Q, image, fp))
    action = brst.action
    params = [f.coupling for f in model.factors.values()] + [n for n, p in model.parameters.items() if p.real]
    for name in params:
        deriv = _parameter_derivative(action, model.param(name))
        add(f"{name} dS0/d{name}", "coupling", deriv, _try_coordinates(Q, deriv, fp))
        phase = _phase_part(model, deriv)
        if not phase.is_zero():
            add(f"phase {name}", "phase", phase, _try_coordinates(Q, phase, fp))
    for factor, fac in model.factors.items():
        if not fac.abelian:
            continue
        for q in global_symmetries(model):
            deform = charge_deformation(model, factor, q, brst)
            label = ", ".join(f"{n}: {v}" for n, v in q.items() if v)
            add(f"charge {factor} [{label}]", "charge", deform, _try_coordinates(Q, deform, fp))
    plain = [k for k, e in enumerate(Q.elements)
             if all(model.fields[n].ghost_number == 0 and not model.fields[n].is_source for c in e.classes for n in c)]
    for k, vec in enumerate(table.restricted(plain)):
        add(f"invariant {k + 1}", "invariant", None, vec)
    kern = table.kernel
    dim = kern.shape[0]
    rank_by_kind: dict[str, int] = {}
    rows: list[np.ndarray] = []
    for kind in ("exact", "coupling", "phase", "charge", "invariant"):
        rows += [d.coordinates for d in dirs if d.kind == kind and d.closed]
        rank_by_kind[kind] = linalg.rank(np.stack(rows), table.field) if rows else 0
    derived = np.stack(rows) if rows else np.zeros((0, len(Q)), dtype=kern.dtype)
    pivots = linalg.independent_columns(np.concatenate([derived, kern]).T, table.field)
    extra = [kern[p - len(rows)] for p in pivots if p >= len(rows)]
    remainder = np.stack(extra) if extra else np.zeros((0, len(Q)), dtype=kern.dtype)
    return KernelAnalysis(dim, dirs, rank_by_kind, rank_by_kind["invariant"], remainder)


# ---------------------------------------------------------------------------------------------- b_d, d dimensions
def bd_singular(singular: SingularBasis, brst: BRST, *, dimensions: Sequence[int] = (4, 6, 8),
                ) -> list[tuple[BasisElement, LocalFunctional, bool]]:
    """``(element, b_d element, closed)`` for every element of the singular basis; ``closed`` is the certified zero
    test :func:`hepkit_renorm.fingerprint.is_zero` in ``dimensions``."""
    model = brst.model
    out = []
    for e in singular.elements:
        image = LocalFunctional.integrate(brst.b_d(e.functional))
        out.append((e, image, is_zero(model, image, dimensions=dimensions)))
    return out


def breaking_by_factor(brst: BRST) -> dict[str, LocalFunctional]:
    """``Delta_1 = sum_r Delta_1^(r)``: the monomials of ``Delta_1`` with a ghost of factor ``r`` (report eq:bd-Oc;
    ``Delta_1`` is linear in the ghosts)."""
    model = brst.model
    out = {}
    for factor in model.factors:
        ghost = f"gh_{factor}"
        part = [m for m in brst.breaking.monomials() if ghost in m.fields()]
        out[factor] = LocalFunctional.integrate(Density.from_monomials(model.fields, part))
    return out


def coupling_direction(brst: BRST, factor: str) -> LocalFunctional:
    """``c_r dS_0/dc_r`` in d dimensions (report ``cal(O)^(r)_g``, continued)."""
    model = brst.model
    return _parameter_derivative(brst.action, model.coupling(factor))


def ghost_direction(brst: BRST, factor: str) -> LocalFunctional:
    """``-b_d int zeta_a omega^a`` of the factor (report ``cal(O)^(r)_c``, d-dimensional)."""
    model = brst.model
    ghost = f"gh_{factor}"
    labels = model.slot_labels(ghost)
    source = model.fields[ghost].source
    y = LocalFunctional.integrate(Density.monomial(model.fields, -1, Occurrence(source, labels),
                                                   Occurrence(ghost, labels)))
    return LocalFunctional.integrate(brst.b_d(y))
