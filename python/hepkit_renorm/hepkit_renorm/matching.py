"""Matching of the summed loop data of the classes onto operator bases (report 02_method sec:m-algorithm, steps 5,
6 and 8 at order ``l``; 04z_st_solution sec:st-dS, sec:st-poles, sec:st-X).

Input data (:class:`ClassData`)
-------------------------------
The per-class results of the diagrams, in the units of :mod:`hepkit_renorm.class_results`: at order ``l`` a class
gives ``T Gamma = (i/(16 pi^2))^l sum c[r, k] eta^r / eps^k``. Two quantities are read (one loop):

* ``"pole"`` (:data:`POLE`): ``N0 = c[0, 1]`` on any samples (``D = 6`` with full spinors here). At ghost number 0
  it equals ``-V_C`` with ``deltaS^(1) = kappa C``, ``kappa = 1/(16 pi^2 eps)`` (eq:m-Z1); at ghost number 1 it is
  ``-V_{b_d C}`` (the pole identity);
* ``"breaking"`` (:data:`BREAKING`): ``P N1 = c[1, 1]`` in ``D = 4`` on physical samples; ``X_1 = sum P N1``
  in units of ``1/(16 pi^2)`` (eq:m-X1).

Providers: :class:`FunctionalData` (local functionals ``C`` and ``X``: golden data, or any symbolic result) and
:class:`LoopData` (the diagram pipeline, :class:`~hepkit_renorm.class_results.ClassComputer`). Both serve values at
the parameter point of a fingerprinter; :class:`FunctionalData` also serves them as polynomials in the free
parameters (:func:`hepkit_renorm.parametric.parametric_rows`), which is what an exact symbolic result needs. Loop data
are values at one point: use a :class:`~hepkit_renorm.parametric.ParameterSpace` with every parameter fixed (to
rational values; the result is then exact at that point).

Outputs
-------
* :func:`match_singular`: the coefficients of ``C`` on the d-dimensional singular basis (continued ``Q`` plus the
  evanescent operators), its four-dimensional limit, the statement that the limit is ``b``-closed (all relations
  that ``b`` imposes on the coefficients hold: the cross-class conditions of sec:st-dS), its decomposition on the
  derived field and parameter directions (the renormalisation constants: bare ``Phi_0 = Z_Phi^(1/2) Phi`` for
  every field multiplet with its source as ``Z_Phi^(-1/2)``, ``p_0 = p + delta p``; for fermions only the physical
  chirality), and the split of the evanescent part into ``b_d``-closed and non-closed operators.
* :func:`pole_identity`: ``N0 + V_{b_d C} = 0`` class by class at ghost number 1.
* :func:`match_breaking`: ``X_l`` in the basis ``N``, the condition ``b X = 0`` on the ghost-number-2 classes, and
  the anomaly analysis: the components of ``X`` outside ``b(span Q)`` (zero exactly when ``X`` is ``b``-exact;
  since ``b X = 0`` they represent its class in the cohomology of ``b``, linear in the anomaly traces), and on
  request the dimension of that cohomology. For the Standard Model it is 13: besides the anomaly directions it
  contains abelian-ghost terms ``omega_B x`` (gauge-invariant operators), which the abelian ghost equation (N4)
  excludes.
* :func:`match_remainder`: the evanescent part ``Evan[sum N_1]`` of the ghost-number-1 classes (eq:m-LP, the
  remainder ``hatDelta^(l)`` of bookkeeping B), on the complete hatted basis of the classes it is asked for.

Exactness: every number is exact over ``GF(p)``; symbolic coefficients are reconstructed with
:func:`hepkit_renorm.parametric.reconstruct` on the graded supports derived by
:class:`~hepkit_renorm.parametric.ParameterSpace`. Results are deterministic (fixed seeds).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from fractions import Fraction
from functools import cached_property
from typing import Protocol, runtime_checkable

import numpy as np
from symbolica import Expression

from . import linalg
from . import tensors as tn
from .bases import BasisElement, OperatorBasis, SingularBasis, breaking_basis, counterterm_basis, singular_basis
from .brst import BRST
from .class_results import ClassComputer, physical_samples
from .classes import _spinor_chirality, leg_classes
from .exact import PrimeField, default_field
from .fingerprint import Fingerprinter, Samples, is_zero
from .functional import Density, Derivation, FieldImage, LocalFunctional, Occurrence
from .parametric import ParameterSpace, ParametricMatrix, Polynomial, evaluate_polynomial, parametric_rows, reconstruct
from .spec import Model

__all__ = [
    "BREAKING",
    "POLE",
    "AnomalyReport",
    "BreakingResult",
    "ClassData",
    "CombinedData",
    "Direction",
    "FunctionalData",
    "LoopData",
    "MatchingError",
    "OrderContext",
    "PoleIdentity",
    "RemainderResult",
    "SingularResult",
    "evanescent_basis",
    "feynman_gauge",
    "field_directions",
    "match_breaking",
    "match_remainder",
    "match_singular",
    "parameter_directions",
    "pole_identity",
]

POLE = "pole"
BREAKING = "breaking"


class MatchingError(ValueError):
    """Loop data that do not lie in the span of the basis (classes listed in ``classes``)."""

    def __init__(self, message: str, classes: Sequence[tuple[str, ...]] = ()) -> None:
        super().__init__(message + (f": {list(classes)}" if classes else ""))
        self.classes = list(classes)


# ---------------------------------------------------------------------------------------------- data providers
@runtime_checkable
class ClassData(Protocol):
    """Summed loop data of the classes (module docstring)."""

    def values(self, quantity: str, fields: tuple[str, ...], fp: Fingerprinter, samples: Samples) -> np.ndarray:
        """Values on ``samples`` at the parameter point of ``fp``."""


class FunctionalData:
    """Class data generated by local functionals: ``counterterm`` ``C`` (``deltaS^(l) = (1/(16 pi^2))^l C/eps``,
    d-dimensional) gives the poles (``-V_C`` at ghost number 0, ``-V_{b_d C}`` at ghost number 1, ``b_d C``
    computed unless ``counterterm_variation`` is given); ``breaking`` ``X`` (``X_l`` in units of
    ``(1/(16 pi^2))^l``) gives the breaking data ``V_X``."""

    def __init__(self, brst: BRST, counterterm: LocalFunctional | None = None,
                 breaking: LocalFunctional | None = None,
                 counterterm_variation: LocalFunctional | None = None) -> None:
        self.brst = brst
        self.counterterm = counterterm
        self.breaking = breaking
        self._variation = counterterm_variation

    @cached_property
    def variation(self) -> LocalFunctional | None:
        if self._variation is not None or self.counterterm is None:
            return self._variation
        return LocalFunctional.integrate(self.brst.b_d(self.counterterm))

    def _functional(self, quantity: str, fields: Sequence[str]) -> tuple[LocalFunctional | None, int]:
        g = sum(self.brst.fields[n].ghost_number for n in fields)
        if quantity == POLE:
            return (self.counterterm, -1) if g == 0 else (self.variation, -1)
        if quantity == BREAKING:
            return self.breaking, 1
        raise ValueError(f"unknown quantity {quantity!r}")

    def values(self, quantity: str, fields: tuple[str, ...], fp: Fingerprinter, samples: Samples) -> np.ndarray:
        f, sign = self._functional(quantity, fields)
        if f is None:
            raise ValueError(f"no functional for {quantity!r} on {fields}")
        return fp.functional(f, samples) * sign % fp.field.p

    def parametric(self, quantity: str, fields: tuple[str, ...], fp: Fingerprinter, space: ParameterSpace,
                   samples: Samples) -> ParametricMatrix:
        f, sign = self._functional(quantity, fields)
        if f is None:
            raise ValueError(f"no functional for {quantity!r} on {fields}")
        return parametric_rows(fp, space, [f], [(fields, samples)]).scale(sign)


class LoopData:
    """Class data from the diagram pipeline (:class:`~hepkit_renorm.class_results.ClassComputer`, order 1):
    ``pole`` is ``N0``, ``breaking`` is ``N1`` (``P N1`` on physical samples in ``D = 4``)."""

    def __init__(self, computer: ClassComputer) -> None:
        if computer.order != 1:
            raise NotImplementedError("the loop readout is implemented at one loop (N0, N1)")
        self.computer = computer
        self.results: dict[tuple, object] = {}

    def values(self, quantity: str, fields: tuple[str, ...], fp: Fingerprinter, samples: Samples) -> np.ndarray:
        res = self.computer.compute(fields, fp, samples)
        self.results[(quantity, fields, fp.field.p, fp.dimension)] = res
        if quantity == POLE:
            return res.N0 % fp.field.p
        if quantity == BREAKING:
            return res.N1 % fp.field.p
        raise ValueError(f"unknown quantity {quantity!r}")


class CombinedData:
    """Class data from several providers: ``providers[class]`` for the listed classes, ``default`` for the others
    (e.g. loop data for a selection of classes and golden data elsewhere)."""

    def __init__(self, default: ClassData, providers: Mapping[tuple[str, ...], ClassData]) -> None:
        self.default = default
        self.providers = dict(providers)

    def _pick(self, fields: tuple[str, ...]) -> ClassData:
        return self.providers.get(tuple(fields), self.default)

    def values(self, quantity: str, fields: tuple[str, ...], fp: Fingerprinter, samples: Samples) -> np.ndarray:
        return self._pick(fields).values(quantity, fields, fp, samples)

    def parametric(self, quantity: str, fields: tuple[str, ...], fp: Fingerprinter, space: ParameterSpace,
                   samples: Samples) -> ParametricMatrix:
        return _data_rows(self._pick(fields), quantity, fp, space, [(tuple(fields), samples)])


def _data_rows(data: ClassData, quantity: str, fp: Fingerprinter, space: ParameterSpace,
               blocks: Sequence[tuple[tuple[str, ...], Samples]]) -> ParametricMatrix:
    parts = []
    for cls, smp in blocks:
        if hasattr(data, "parametric"):
            parts.append(data.parametric(quantity, cls, fp, space, smp).padded(len(space)))
        elif len(space):
            raise NotImplementedError("loop data are values at one parameter point: fix every parameter of the "
                                      "ParameterSpace (an exact result at that point)")
        else:
            parts.append(ParametricMatrix.from_values(fp.field, data.values(quantity, cls, fp, smp)))
    if not parts:
        return ParametricMatrix(fp.field, (0, 1), {})
    return ParametricMatrix.vstack(parts)


# ---------------------------------------------------------------------------------------------- context
def feynman_gauge(model: Model) -> dict[str, int]:
    """``{name: 1}`` for the parameter that is the gauge-fixing parameter (Feynman gauge, as in the diagrams);
    empty if the gauge-fixing parameter is a number."""
    xi = model.gauge_fixing.xi
    for name in model.parameters:
        if model.param(name) == xi:
            return {name: 1}
    return {}


class OrderContext:
    """Everything an order-``l`` step shares: the model, ``BRST``, the parameter space, the bases (built on first
    use), fingerprinters and sample blocks per prime (deterministic: ``seed``).

    ``fixed`` defaults to Feynman gauge (the gauge-fixing parameter set to 1, as in the diagrams); pass all
    parameters at rational values for a computation at one point (loop data).
    """

    def __init__(self, model: Model, order: int = 1, *, brst: BRST | None = None,
                 fixed: Mapping[str, Fraction | int] | None = None, Q: OperatorBasis | None = None,
                 N: OperatorBasis | None = None, singular: SingularBasis | None = None, seed: int = 0,
                 margin: int = 4) -> None:
        if order < 1:
            raise ValueError("the order of a step is l >= 1")
        self.model = model
        self.order = order
        self.brst = brst or BRST(model)
        self.space = ParameterSpace(model, self.brst.action, feynman_gauge(model) if fixed is None else fixed)
        self.seed = seed
        self.margin = margin
        self._Q, self._N, self._singular = Q, N, singular
        self._fps: dict[tuple[int, int, int], Fingerprinter] = {}
        self._cache: dict[tuple, object] = {}

    # ------------------------------------------------------------------ bases
    @property
    def Q(self) -> OperatorBasis:
        if self._Q is None:
            self._Q = counterterm_basis(self.model)
        return self._Q

    @property
    def N(self) -> OperatorBasis:
        if self._N is None:
            self._N = breaking_basis(self.model)
        return self._N

    @property
    def singular(self) -> SingularBasis:
        if self._singular is None:
            self._singular = singular_basis(self.Q)
        return self._singular

    # ------------------------------------------------------------------ evaluation set-up
    def fingerprinter(self, field_: PrimeField, dimension: int, i_sign: int = 1) -> Fingerprinter:
        key = (field_.p, dimension, i_sign)
        if key not in self._fps:
            self._fps[key] = Fingerprinter(self.model, field_, dimension, seed=self.seed, values=self.space.fixed,
                                           i_sign=i_sign)
        return self._fps[key]

    def blocks(self, fp: Fingerprinter, counts: Sequence[tuple[tuple[str, ...], int]], salt: int,
               physical: bool = False) -> list[tuple[tuple[str, ...], Samples]]:
        out = []
        for cls, n in counts:
            smp = fp.samples(list(cls), n, salt)
            out.append((cls, physical_samples(fp, smp) if physical else smp))
        return out

    def cached(self, key: tuple, build):
        if key not in self._cache:
            self._cache[key] = build()
        return self._cache[key]

    def replace_cached(self, key: tuple, value) -> None:
        self._cache[key] = value

    def rows(self, fp: Fingerprinter, functionals: Sequence[LocalFunctional],
             blocks: Sequence[tuple[tuple[str, ...], Samples]]) -> ParametricMatrix:
        return parametric_rows(fp, self.space, functionals, blocks)

    def classes(self, ghost_number: int) -> list[tuple[str, ...]]:
        return [c.fields for c in leg_classes(self.model, ghost_number)]

    def structure_counts(self, ghost_number: int) -> list[tuple[tuple[str, ...], int]]:
        """``(class, number of structures + margin)`` for every leg class of a ghost number: enough rows to
        represent any functional of the class exactly."""
        return self.cached(("counts", ghost_number), lambda: [(c.fields, c.n_structures + self.margin)
                                                              for c in leg_classes(self.model, ghost_number)])

    def support(self, element: LocalFunctional) -> list[tuple[int, ...]]:
        """The support of the coefficient of ``element`` in an order-``l`` functional (grade ``2 l``, excess 0);
        the constant monomial when every parameter is fixed."""
        if not len(self.space):
            return [()]
        grade, excess = self.space.grade(element)
        return self.space.support(2 * self.order - grade, excess)


def _counts(elements: Sequence[BasisElement], classes: Sequence[tuple[str, ...]], margin: int,
            ) -> list[tuple[tuple[str, ...], int]]:
    per = {c: margin for c in classes}
    for e in elements:
        for c in e.classes:
            if c in per:
                per[c] += 1
    return [(c, per[c]) for c in classes]


# ---------------------------------------------------------------------------------------------- directions
@dataclass(frozen=True)
class Direction:
    """A derived generator of closed directions: ``kind`` ``"field"`` (``O_Phi``, coefficient ``delta Z_Phi / 2``)
    or ``"parameter"`` (``dS_0/dp``, coefficient ``delta p``); ``fields`` the fields it rescales (``+1``) and the
    sources (``-1``)."""

    name: str
    kind: str
    functional: LocalFunctional
    fields: tuple[tuple[str, int], ...] = ()


def _chiral_image(model: Model, name: str, labels: tuple[Expression, ...]) -> Density:
    """``P_c Psi`` (``Psi Pbar_c`` for ``Psibar``) on the spinor slot, the field itself otherwise."""
    fi = model.fields[name]
    spin = [k for k, s in enumerate(fi.slots) if s.space.kind == "spinor"]
    if fi.kind not in ("fermion", "fermion_bar") or not spin:
        return Density.monomial(model.fields, 1, Occurrence(name, labels))
    k = spin[0]
    inner = tn.fresh(tn.SPINOR)
    proj = tn.PL if _spinor_chirality(model, fi) == "L" else tn.PR
    slots = labels[:k] + (inner,) + labels[k + 1:]
    coeff = proj(labels[k], inner) if fi.kind == "fermion" else proj(inner, labels[k])
    return Density.monomial(model.fields, coeff, Occurrence(name, slots))


def field_directions(model: Model, action: LocalFunctional) -> list[Direction]:
    """``O_Phi = int (Phi dS_0/dPhi - Phi^* dS_0/dPhi^*)`` for every field multiplet that is renormalised: the
    gauge field and the ghost of each factor, every fermion multiplet (physical chirality, field and conjugate)
    and every scalar multiplet (field and conjugate), each with its BRST source (report eq:sts-bare); the
    antighost and the auxiliary field of a factor scale inversely to its gauge field, like its source, so that
    the shifted source ``rhotilde = rho + d omegabar`` scales as a whole (the gauge-fixing parameter is not
    rescaled: its classes with auxiliary legs are never compared). In ``d`` dimensions a fermion direction
    also contains the evanescent kinetic term ``psibar dhat chi + chibar dhat psi``, which the counterterm does not
    contain (the sterile partner is not rescaled); the directions are compared in four dimensions."""
    groups: dict[str, list[str]] = {}
    partners: dict[str, list[str]] = {}
    for name, fi in model.fields.items():
        if fi.kind in ("gauge", "ghost"):
            groups.setdefault(name, []).append(name)
        elif fi.kind in ("fermion", "fermion_bar", "scalar", "scalar_conj"):
            groups.setdefault(fi.multiplet, []).append(name)
        elif fi.kind in ("antighost", "auxiliary") and fi.factor is not None:
            partners.setdefault(fi.factor, []).append(name)
    out = []
    for key, members in groups.items():
        images, signed = {}, []
        inverse = [model.fields[n].source for n in members if model.fields[n].source is not None]
        if model.fields[members[0]].kind == "gauge":
            inverse += partners.get(model.fields[members[0]].factor, [])   # rhotilde = rho + d omegabar scales as rho
        for name in members:
            labels = model.slot_labels(name)
            images[name] = FieldImage(labels, _chiral_image(model, name, labels))
            signed.append((name, 1))
        for name in inverse:
            sl = model.slot_labels(name)
            images[name] = FieldImage(sl, Density.monomial(model.fields, -1, Occurrence(name, sl)))
            signed.append((name, -1))
        image = LocalFunctional.integrate(Derivation(0, images)(action))
        if not image.is_zero():
            out.append(Direction(f"Z[{key}]", "field", image, tuple(signed)))
    return out


def parameter_directions(model: Model, action: LocalFunctional) -> list[Direction]:
    """``dS_0/dp`` for every coupling and real parameter except the gauge-fixing parameter (the coefficient is
    ``delta p``; the gauge fixing is not renormalised)."""
    out = []
    names = [f.coupling for f in model.factors.values()] + [n for n, p in model.parameters.items() if p.real]
    for name in names:
        if model.param(name) == model.gauge_fixing.xi:
            continue
        sym = model.param(name)
        deriv = LocalFunctional.integrate(action.map_coefficients(lambda c, s=sym: c.derivative(s).expand()))
        if not deriv.is_zero():
            out.append(Direction(f"d[{name}]", "parameter", deriv))
    return out


# ---------------------------------------------------------------------------------------------- the singular counterterm
@dataclass
class SingularResult:
    """``deltaS^(l) = (1/(16 pi^2))^l C / eps`` with ``C = sum_alpha c_alpha S_alpha`` on the singular basis.

    ``coefficients[name]``: polynomial ``c_alpha``; ``four_dimensional[name]``: the coefficients of the ``Q``
    elements (the four-dimensional limit is ``sum c_alpha Q_alpha``); ``constants[name]``: the coefficients on
    the derived directions (``Z[...]``: ``delta Z_Phi``, i.e. twice the coefficient of ``O_Phi``; ``d[p]``:
    ``delta p``), all without the factor ``(1/(16 pi^2))^l / eps``; ``evanescent``: ``{name: (polynomial,
    b_d-closed)}``; ``closure_rank``: number of independent linear conditions that ``b`` imposes on the ``Q``
    coefficients (all hold); ``fit_rows``: number of exact equations matched.
    """

    context: OrderContext
    elements: list[BasisElement]
    coefficients: dict[str, Polynomial]
    four_dimensional: dict[str, Polynomial]
    constants: dict[str, Polynomial]
    evanescent: dict[str, tuple[Polynomial, bool]]
    directions: list[Direction]
    closure_rank: int
    fit_rows: int
    unmatched_directions: list[str] = field(default_factory=list)

    def functional(self) -> LocalFunctional:
        """``C`` as a local functional (d-dimensional, symbolic in the free parameters)."""
        return _combine(self.context, self.elements, self.coefficients)

    def four_dimensional_functional(self) -> LocalFunctional:
        return _combine(self.context, self.context.Q.elements, self.four_dimensional)

    def relative(self, name: str) -> Expression:
        """``delta p / p`` for a parameter (the report's ``delta Z_{c_r}`` for a gauge coupling)."""
        return (self.context.space.expression(self.constants[f"d[{name}]"]) / self.context.model.param(name)).expand()


def _combine(ctx: OrderContext, elements: Sequence[BasisElement], coefficients: Mapping[str, Polynomial],
             ) -> LocalFunctional:
    total = Density(ctx.model.fields)
    for e in elements:
        poly = coefficients.get(e.name)
        if poly:
            total = total + e.functional.scale(ctx.space.expression(poly))
    return LocalFunctional.integrate(total)


def _solve_point(a: np.ndarray, b: np.ndarray, field_: PrimeField, blocks_layout, what: str) -> np.ndarray:
    try:
        return linalg.solve(a, b, field_, require_unique=True)
    except linalg.InconsistentSystemError:
        bad = _failing_classes(a, b, field_, blocks_layout)
        raise MatchingError(f"{what}: the data are not in the span of the basis", bad) from None


def _failing_classes(a, b, field_, layout) -> list[tuple[str, ...]]:
    bad, start = [], 0
    for cls, n in layout:
        try:
            linalg.solve(a[start:start + n], b[start:start + n], field_)
        except linalg.InconsistentSystemError:
            bad.append(cls)
        start += n
    return bad


def match_singular(ctx: OrderContext, data: ClassData, *, dimensions: Sequence[int] = (6,),
                   max_primes: int = 6) -> SingularResult:
    """Match the ghost-number-0 poles onto the singular basis (module docstring)."""
    model, space = ctx.model, ctx.space
    sing = ctx.singular
    elements = sing.elements
    Q = ctx.Q
    nq = len(Q)
    directions = field_directions(model, ctx.brst.action) + parameter_directions(model, ctx.brst.action)
    classes0 = ctx.classes(0)
    counts = _counts(elements, classes0, ctx.margin)
    counts1 = ctx.structure_counts(1)
    layout = [(c, n) for _ in dimensions for c, n in counts]
    info: dict[str, object] = {}

    def prepare(field_: PrimeField, i_sign: int):
        mats, rhs = [], []
        for D in dimensions:
            fp = ctx.fingerprinter(field_, D, i_sign)
            blocks = ctx.blocks(fp, counts, salt=11)
            mats.append(ctx.rows(fp, [e.functional for e in elements], blocks))
            rhs.append(_data_rows(data, POLE, fp, space, blocks))
        A, R = ParametricMatrix.vstack(mats), ParametricMatrix.vstack(rhs)
        fp4 = ctx.fingerprinter(field_, 4, i_sign)
        # b on the four-dimensional limit: images of Q on the ghost-number-1 classes (closure condition)
        images = ctx.cached(("bQ",), lambda: [LocalFunctional.integrate(ctx.brst.b(e.functional)) for e in Q])
        bq = ctx.rows(fp4, images, ctx.blocks(fp4, counts1, salt=13))
        # the derived directions in Q coordinates: Q fingerprints (D = 4) against the directions
        qblocks = ctx.blocks(fp4, list(Q.layout), salt=17)
        qm = ctx.rows(fp4, Q.functionals, qblocks)
        dm = ctx.rows(fp4, [d.functional for d in directions], qblocks)

        def evaluate(point: np.ndarray) -> np.ndarray:
            p = field_.p
            c = _solve_point(A.at(point), (-R.at(point)) % p, field_, layout, "singular counterterm")[:, 0]
            cq = c[:nq]
            if np.any(bq.at(point) @ cq.astype(object) % p):
                raise MatchingError("the four-dimensional limit of the singular counterterm is not b-closed")
            qa = qm.at(point)
            dirs = linalg.solve(qa, dm.at(point), field_)     # Q coordinates of the directions (columns)
            if "closure_rank" not in info:
                info["closure_rank"] = linalg.rank(bq.at(point), field_)
                info["dir_rank"] = linalg.rank(dirs, field_)
                info["independent"] = linalg.independent_columns(dirs, field_)
            ind = info["independent"]
            try:
                z = linalg.solve(dirs[:, ind], cq, field_)
            except linalg.InconsistentSystemError:
                raise MatchingError("the four-dimensional limit is not spanned by the field and parameter "
                                    "directions") from None
            zz = np.zeros(len(directions), dtype=np.int64)
            zz[ind] = z
            for k, d in enumerate(directions):
                if d.kind == "field":
                    zz[k] = 2 * zz[k] % p          # delta Z = 2 x coefficient of O_Phi
            return np.concatenate([c, zz])

        return evaluate

    supports = [ctx.support(e.functional) for e in elements] + [ctx.support(d.functional) for d in directions]
    polys = reconstruct(prepare, supports, space, max_primes=max_primes, seed=ctx.seed)
    coeffs = {e.name: polys[k] for k, e in enumerate(elements)}
    constants = {d.name: polys[len(elements) + k] for k, d in enumerate(directions)}
    four = {e.name: coeffs[e.name] for e in Q.elements}
    ind = set(info["independent"])
    unmatched = [d.name for k, d in enumerate(directions) if k not in ind]
    closed = {}
    for e in sing.evanescent:
        if coeffs[e.name]:
            image = LocalFunctional.integrate(ctx.brst.b_d(e.functional))
            closed[e.name] = (coeffs[e.name], is_zero(model, image, dimensions=(4, 6, 8)))
    return SingularResult(ctx, elements, coeffs, four, constants, closed, directions, int(info["closure_rank"]),
                          sum(n for _, n in layout), unmatched)


# ---------------------------------------------------------------------------------------------- pole identity
@dataclass
class PoleIdentity:
    """``N0 + V_{b_d C} = 0`` on every ghost-number-1 class (``failed``: the classes where it does not hold),
    checked at ``points`` random parameter points in ``dimension``."""

    classes: list[tuple[str, ...]]
    failed: list[tuple[str, ...]]
    dimension: int
    points: int

    @property
    def holds(self) -> bool:
        return not self.failed


def pole_identity(ctx: OrderContext, singular: SingularResult, data: ClassData, *, dimension: int = 6,
                  samples: int = 2, points: int = 2, classes: Sequence[tuple[str, ...]] | None = None
                  ) -> PoleIdentity:
    """The pole identity of sec:m-oneloop: the ghost-number-1 poles equal ``-V_{b_d C}``."""
    field_ = default_field()
    fp = ctx.fingerprinter(field_, dimension)
    cls = list(classes) if classes is not None else ctx.classes(1)
    blocks = ctx.blocks(fp, [(c, samples) for c in cls], salt=19)
    elements = [e for e in singular.elements if singular.coefficients.get(e.name)]
    images = [LocalFunctional.integrate(ctx.brst.b_d(e.functional)) for e in elements]
    V = ctx.rows(fp, images, blocks)
    R = _data_rows(data, POLE, fp, ctx.space, blocks)
    pts = ctx.space.random_points(field_, points, seed=ctx.seed + 101) if len(ctx.space) else np.zeros((1, 0), int)
    failed: set[tuple[str, ...]] = set()
    p = field_.p
    for pt in pts:
        c = np.array([evaluate_polynomial(singular.coefficients[e.name], pt, field_) for e in elements],
                     dtype=object)
        total = (V.at(pt).astype(object) @ c + R.at(pt)[:, 0].astype(object)) % p if elements else R.at(pt)[:, 0]
        start = 0
        for (c_, smp) in blocks:
            if np.any(np.asarray(total[start:start + smp.n], dtype=object) % p):
                failed.add(c_)
            start += smp.n
    return PoleIdentity(cls, sorted(failed, key=lambda c: cls.index(c)), dimension, len(pts))


# ---------------------------------------------------------------------------------------------- the breaking
@dataclass
class AnomalyReport:
    """The part of ``X`` outside ``b(span Q)``, at the analysis point over ``field``.

    ``complement``: the names of the ``N`` monomials that complete ``b(span Q)`` to ``span N`` (pivots in the order
    of ``N``); ``coefficients``: the components of ``X`` along them (all zero iff ``X`` is ``b``-exact; since
    ``b X = 0`` they represent the class of ``X`` in the cohomology of ``b``, i.e. the anomaly, linear in the
    anomaly traces). ``cohomology``: ``dim ker(b on N) - rank b(span Q)`` (the number of independent anomaly
    directions), computed only on request (it needs ``b N`` on enough rows of every ghost-number-2 class);
    ``image_rank``: ``rank b(span Q)``."""

    complement: list[str]
    coefficients: np.ndarray
    field: PrimeField
    point: np.ndarray
    image_rank: int
    cohomology: int | None = None

    @property
    def anomalous(self) -> bool:
        return bool(np.any(self.coefficients % self.field.p))

    def components(self) -> dict[str, int]:
        """The non-zero components ``{N name: value over the field}``."""
        return {n: int(c) for n, c in zip(self.complement, self.coefficients) if c % self.field.p}


@dataclass
class BreakingResult:
    """``X_l = sum_beta x_beta N_beta`` (``(1/(16 pi^2))^l`` units) at the analysis point (``coordinates`` over
    ``field``), optionally reconstructed (``polynomials``); ``closed``: ``b X = 0`` on every ghost-number-2 class
    that ``b N`` reaches (``closure_rows`` rows); ``anomaly``: the part of ``X`` outside ``b(span Q)``."""

    context: OrderContext
    field: PrimeField
    point: np.ndarray
    coordinates: np.ndarray
    closed: bool
    closure_rows: int
    anomaly: AnomalyReport
    polynomials: dict[str, Polynomial] | None = None

    def functional(self) -> LocalFunctional:
        if self.polynomials is None:
            raise ValueError("reconstruct the coordinates first (match_breaking(..., reconstruct_coordinates=True))")
        return _combine(self.context, self.context.N.elements, self.polynomials)


def _breaking_rows(ctx: OrderContext, data: ClassData, fp: Fingerprinter):
    N = ctx.N
    counts = [(c, n) for c, n in N.layout]
    blocks = ctx.blocks(fp, counts, salt=23, physical=True)
    nm = ctx.rows(fp, N.functionals, blocks)
    xr = _data_rows(data, BREAKING, fp, ctx.space, blocks)
    return blocks, nm, xr


def _closure_rows(ctx: OrderContext, fp: Fingerprinter, full: bool):
    images = ctx.cached(("bN",), lambda: [LocalFunctional.integrate(ctx.brst.b(e.functional)) for e in ctx.N])
    touching: dict[tuple[str, ...], int] = {}
    for img in images:
        for c in img.leg_classes():
            touching[c] = touching.get(c, 0) + 1
    counts = [(c.fields, min(c.n_structures, touching[c.fields]) + ctx.margin if full else 2)
              for c in leg_classes(ctx.model, 2) if c.fields in touching]
    return ctx.rows(fp, images, ctx.blocks(fp, counts, salt=29))


def match_breaking(ctx: OrderContext, data: ClassData, *, cohomology: bool = False,
                   reconstruct_coordinates: bool = False) -> BreakingResult:
    """``X_l`` in the basis ``N``, the consistency condition ``b X = 0`` and the anomaly analysis (module
    docstring) at a random analysis point over the default prime.

    ``b X = 0`` is tested on two samples of every ghost-number-2 class that ``b N`` reaches (a non-zero ``b X``
    is detected with probability ``>= 1 - deg/p`` per sample). ``cohomology=True`` also computes the dimension of
    the cohomology of ``b`` on ``N``, with ``b N`` on as many rows per class as its image can have dimensions
    there (slow for the Standard Model: minutes). ``reconstruct_coordinates`` reconstructs the coordinates as
    complex polynomials in the free parameters."""
    field_ = default_field()
    p = field_.p
    fp = ctx.fingerprinter(field_, 4)
    N, Q = ctx.N, ctx.Q
    blocks, nm, xr = _breaking_rows(ctx, data, fp)
    point = ctx.space.random_points(field_, 1, seed=ctx.seed + 202)[0] if len(ctx.space) else np.zeros(0, int)
    layout = [(c, s.n) for c, s in blocks]
    x = _solve_point(nm.at(point), xr.at(point), field_, layout, "breaking")[:, 0]
    bn = _closure_rows(ctx, fp, cohomology).at(point)
    closed = not np.any(bn.astype(object) @ x.astype(object) % p)
    qimages = ctx.cached(("bQ",), lambda: [LocalFunctional.integrate(ctx.brst.b(e.functional)) for e in Q])
    bq = linalg.solve(nm.at(point), ctx.rows(fp, qimages, blocks).at(point), field_).T    # rows: b Q_alpha
    rank_q = linalg.rank(bq, field_)
    eye = np.eye(len(N), dtype=np.int64)
    piv = linalg.independent_columns(np.concatenate([bq, eye]).T, field_)
    comp = [k - bq.shape[0] for k in piv if k >= bq.shape[0]]
    basis = np.concatenate([bq, eye[comp]])
    ya = linalg.solve(basis.T, x, field_)
    a = np.asarray(ya[bq.shape[0]:], dtype=np.int64)
    dim = None
    if cohomology:
        dim = linalg.nullspace(bn, field_).shape[0] - rank_q
    anomaly = AnomalyReport([N.elements[k].name for k in comp], a, field_, point, rank_q, dim)
    polys = None
    if reconstruct_coordinates:
        def prepare(f: PrimeField, i_sign: int):
            fpp = ctx.fingerprinter(f, 4, i_sign)
            bl, nmm, xrr = _breaking_rows(ctx, data, fpp)
            lay = [(c, s.n) for c, s in bl]
            return lambda pt: _solve_point(nmm.at(pt), xrr.at(pt), f, lay, "breaking")[:, 0]
        supports = [ctx.support(e.functional) for e in N.elements]
        res = reconstruct(prepare, supports, ctx.space, complex_values=True, seed=ctx.seed)
        polys = {e.name: res[k] for k, e in enumerate(N.elements)}
    return BreakingResult(ctx, field_, point, x, closed, bn.shape[0], anomaly, polys)


# ---------------------------------------------------------------------------------------------- bookkeeping B
@dataclass
class RemainderResult:
    """``Evan[sum N_1]`` of the ghost-number-1 classes (eq:m-LP, bookkeeping B) on an evanescent basis of those
    classes: ``elements`` (hatted candidates of :func:`hepkit_renorm.bases.candidates`, reduced in ``D = 6, 8``),
    ``coefficients`` (polynomials), in ``(1/(16 pi^2))^l`` units."""

    context: OrderContext
    elements: list[BasisElement]
    coefficients: dict[str, Polynomial]
    classes: list[tuple[str, ...]]

    def functional(self) -> LocalFunctional:
        return _combine(self.context, self.elements, self.coefficients)


def evanescent_basis(ctx: OrderContext, classes: Sequence[tuple[str, ...]], dimensions: Sequence[int] = (6, 8)
                     ) -> list[BasisElement]:
    """Independent evanescent local monomials of the given classes (complete hatted basis: hatted derivatives and
    hatted field components, :func:`hepkit_renorm.bases.candidates` with ``hats`` and ``hatted_fields``), reduced
    by the exact rank of fingerprints in ``dimensions`` jointly."""
    from .bases import _Reducer, candidates

    model = ctx.model
    info = {c.fields: c for g in sorted({sum(model.fields[n].ghost_number for n in cls) for cls in classes})
            for c in leg_classes(model, g)}
    out: list[BasisElement] = []
    for cls in classes:
        lc = info[tuple(cls)]
        cands = []
        for degree in range(len(lc.structures)):
            for mono in candidates(model, cls, degree, hats=True, hatted_fields=True):
                f = LocalFunctional.integrate(Density.from_monomials(model.fields, [mono]))
                if not f.is_zero():
                    cands.append((degree, f))
        n = len(cands) + ctx.margin
        fps = [Fingerprinter(model, default_field(), D, seed=ctx.seed) for D in dimensions]
        reducer = _Reducer([(fp, fp.samples(list(cls), n, salt=41)) for fp in fps])
        for degree, f in cands:
            if reducer.offer(f):
                out.append(BasisElement(f"R{len(out) + 1}[{' '.join(cls)}]", f, tuple(cls), (tuple(cls),), degree,
                                        "evanescent", str(f.expr)))
    return out


def match_remainder(ctx: OrderContext, data: ClassData, restoring_functional: LocalFunctional,
                    classes: Sequence[tuple[str, ...]], *, dimensions: Sequence[int] = (6, 8),
                    max_primes: int = 6) -> RemainderResult:
    """The evanescent part of the ``eta/eps`` coefficients of the ghost-number-1 classes (eq:m-LP):
    ``Evan[sum N_1] = sum N_1 - V_{X_l}`` in ``D > 4`` on full spinors, with ``X_l = -b DeltaS^(l)`` from the
    restoring counterterm (``restoring_functional``, four-dimensional), matched onto :func:`evanescent_basis`. The
    data must provide the quantity ``"breaking"`` in every dimension of ``dimensions``."""
    space = ctx.space
    elements = evanescent_basis(ctx, classes, dimensions)
    X = LocalFunctional.integrate(ctx.brst.b(restoring_functional).scale(-1))
    counts = [(tuple(c), sum(1 for e in elements if e.leg_class == tuple(c)) + ctx.margin) for c in classes]
    layout = [(c, n) for _ in dimensions for c, n in counts]

    def prepare(field_: PrimeField, i_sign: int):
        mats, rhs = [], []
        for D in dimensions:
            fp = ctx.fingerprinter(field_, D, i_sign)
            blocks = ctx.blocks(fp, counts, salt=43)
            mats.append(ctx.rows(fp, [e.functional for e in elements], blocks))
            data_rows = _data_rows(data, BREAKING, fp, space, blocks)
            xr = ctx.rows(fp, [X], blocks)
            rhs.append(ParametricMatrix.hstack([data_rows, xr.scale(-1)]))
        A, R = ParametricMatrix.vstack(mats), ParametricMatrix.vstack(rhs)

        def evaluate(point: np.ndarray) -> np.ndarray:
            r = R.at(point)
            b = (r[:, 0] + r[:, 1]) % field_.p
            return _solve_point(A.at(point), b[:, None], field_, layout, "evanescent remainder")[:, 0]

        return evaluate

    supports = [ctx.support(e.functional) for e in elements]
    polys = reconstruct(prepare, supports, space, max_primes=max_primes, seed=ctx.seed, complex_values=True)
    return RemainderResult(ctx, elements, {e.name: polys[k] for k, e in enumerate(elements)}, list(classes))
