"""The restoring condition ``b DeltaS^(l) = -X_l`` with normalisation conditions (report 02_method
sec:m-algorithm step 9; 03b sec:normalisation; 04z sec:st-solve).

``DeltaS^(l) = (1/(16 pi^2))^l sum_alpha r_alpha Q_alpha`` over a ghost-number-0 basis (default: the generated
hermitian basis :func:`hepkit_renorm.bases.counterterm_basis`, or any :class:`~hepkit_renorm.bases.OperatorBasis`).
The equation is imposed on fingerprint rows of every ghost-number-1 class (``D = 4``, physical external spinors):
``sum_alpha r_alpha V_{b Q_alpha} = -V_{X_l}`` with ``X_l`` from the class data (:mod:`hepkit_renorm.matching`).
The solution is fixed up to the kernel of ``b`` (the closed directions); normalisation conditions remove it.

Normalisation conditions (model independent)
--------------------------------------------
A condition is a set of linear equations on the coefficients ``r`` (constraint rows), stated for any field
content. :func:`standard_normalisation` gives the report's (N1)-(N4) in this form:

* (N1) ``DeltaS`` has no part that survives when all gauge fields, ghosts and sources vanish (no pure-matter
  classes: kinetic, mass, potential and Yukawa monomials);
* (N2) for each gauge factor, the derivative part of the two-point function of ``DeltaS`` is longitudinal
  (``V^{mu nu}(p) e_mu = 0`` for ``e.p = 0``): no transverse kinetic term ``-1/4 (dA - dA)^2``;
* (N3) for each gauge factor, no ghost-source two-point term ``rhotilde.d omega`` (class ``(src, gh)``) and, for a
  non-abelian factor, no part of the ``(rhotilde, A, omega)`` vertex that is antisymmetric in the adjoint indices
  of ``A`` and ``omega`` (the structure-constant term ``rhotilde fg A omega``);
* (N4) no monomial with the ghost of an abelian factor (the abelian ghost equation).

The conditions are applied in this priority order. In the default mode (``"full"``) every equation of a condition
is imposed when the restoring condition can still be solved; otherwise (with ``fallback``) the condition is relaxed
automatically to its *pivot* equations (those needed to fix the remaining closed directions, which never obstruct
existence) plus every further equation that keeps the system solvable, in order, and the relaxation is reported.
For (N4) this is the generalisation of the report's (N4'). Mode ``"minimal"`` imposes the pivot equations only
(the minimal removed set of the driver design). ``fermion_kinetic=True`` keeps the fermion kinetic terms (the
representative with fermion kinetic counterterms, for comparison with the literature): the fermion two-point
classes are exempt from (N1), and the fermion field directions are fixed instead by removing the gauge vertex of
the first gauge factor that acts on each multiplet (N1k).

Reported (:class:`RestoringResult`): the coefficients (exact polynomials in the free parameters), uniqueness as
a rank statement (rank of ``b`` on the allowed space against its dimension), the existence conditions (``rank b
(span Q) - rank b(allowed)``, located class by class, all satisfied by ``X``), the conditions with their applied
and relaxed equations. A breaking outside ``b(span Q)`` raises :class:`AnomalyError` with the anomaly analysis of
:func:`hepkit_renorm.matching.match_breaking` (the components of ``X`` along the cohomology of ``b``).
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field

import numpy as np

from . import linalg
from .bases import OperatorBasis
from .exact import PrimeField
from .fingerprint import Fingerprinter, Samples
from .functional import LocalFunctional
from .matching import BREAKING, AnomalyReport, ClassData, OrderContext, _data_rows, match_breaking
from .parametric import ParametricMatrix, Polynomial, reconstruct
from .spec import Model

__all__ = [
    "AnomalyError",
    "ClassCondition",
    "Condition",
    "ConditionReport",
    "ElementCondition",
    "ExistenceError",
    "Normalisation",
    "ProjectedCondition",
    "RestoringResult",
    "antisymmetric_vertex",
    "solve_restoring",
    "standard_normalisation",
    "transverse_two_point",
]

_MATTER = ("fermion", "fermion_bar", "scalar", "scalar_conj")


class AnomalyError(ValueError):
    """The breaking has a component outside ``b(span Q)``: the gauge anomaly (``report``: the analysis)."""

    def __init__(self, report: AnomalyReport) -> None:
        comps = report.components()
        super().__init__(f"the breaking is not b-exact (gauge anomaly): non-zero components along {sorted(comps)}, "
                         f"N monomials outside b(span Q) (rank {report.image_rank})")
        self.report = report


class ExistenceError(ValueError):
    """A normalisation condition admits no solution (and ``fallback`` is off)."""


# ---------------------------------------------------------------------------------------------- conditions
class Condition:
    """A normalisation condition: constraint rows on the coefficients of a basis (module docstring)."""

    name: str
    description: str

    def rows(self, ctx: OrderContext, basis: OperatorBasis, fp: Fingerprinter) -> ParametricMatrix:
        raise NotImplementedError


@dataclass
class ClassCondition(Condition):
    """``DeltaS`` has no part in the classes selected by ``predicate(model, class)``: one unit row per basis
    element living in such a class."""

    name: str
    description: str
    predicate: Callable[[Model, tuple[str, ...]], bool]

    def rows(self, ctx: OrderContext, basis: OperatorBasis, fp: Fingerprinter) -> ParametricMatrix:
        sel = [k for k, e in enumerate(basis.elements) if any(self.predicate(basis.model, c) for c in e.classes)]
        m = np.zeros((len(sel), len(basis)), dtype=np.int64)
        for row, k in enumerate(sel):
            m[row, k] = 1
        return ParametricMatrix(fp.field, m.shape, {(0,) * len(ctx.space): m} if len(sel) else {})


@dataclass
class ElementCondition(Condition):
    """No component along the named basis elements (a removal set given by hand, e.g. the report's lists of
    sec:normalisation for its own basis)."""

    name: str
    description: str
    elements: tuple[str, ...]

    def rows(self, ctx: OrderContext, basis: OperatorBasis, fp: Fingerprinter) -> ParametricMatrix:
        m = np.zeros((len(self.elements), len(basis)), dtype=np.int64)
        for row, name in enumerate(self.elements):
            m[row, basis.index(name)] = 1
        return ParametricMatrix(fp.field, m.shape, {(0,) * len(ctx.space): m} if len(m) else {})


@dataclass
class ProjectedCondition(Condition):
    """``sum_j sign_j V_{DeltaS}(samples_j) = 0`` on transformed samples of one class (``transform(fp, samples)``
    returns ``[(sign, samples), ...]``), for the basis elements of the class with ``degree >= min_degree``
    (elements of mixed degree included)."""

    name: str
    description: str
    leg_class: tuple[str, ...]
    transform: Callable[[Fingerprinter, Samples], list[tuple[int, Samples]]]
    min_degree: int = 0

    def rows(self, ctx: OrderContext, basis: OperatorBasis, fp: Fingerprinter) -> ParametricMatrix:
        cls = self.leg_class
        sel = [k for k, e in enumerate(basis.elements) if cls in e.classes and (e.degree < 0 or e.degree >=
                                                                                  self.min_degree)]
        n = len(sel) + ctx.margin
        base = fp.samples(list(cls), n, salt=37)
        total = None
        for sign, smp in self.transform(fp, base):
            part = ctx.rows(fp, [basis.elements[k].functional for k in sel], [(cls, smp)]).scale(sign)
            total = part if total is None else _add(total, part)
        out = ParametricMatrix(fp.field, (n, len(basis)), {})
        for e, m in (total.terms.items() if total is not None else []):
            full = np.zeros((n, len(basis)), dtype=np.int64)
            full[:, sel] = m
            out.terms[e] = full
        return out


def _add(a: ParametricMatrix, b: ParametricMatrix) -> ParametricMatrix:
    terms = dict(a.terms)
    for e, m in b.terms.items():
        terms[e] = (terms[e] + m) % a.field.p if e in terms else m
    return ParametricMatrix(a.field, a.shape, terms)


def _dot(fp: Fingerprinter, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """``a.b`` with the metric over the Lorentz axis (last axis), per sample."""
    p = fp.field.p
    g = fp.rep.metric.astype(object)
    return (np.asarray(a, dtype=object) * g * np.asarray(b, dtype=object)).sum(axis=-1) % p


def transverse_two_point(model: Model, factor: str) -> ProjectedCondition:
    """(N2) for one gauge factor: the derivative part of the ``(A, A)`` vertex of ``DeltaS`` vanishes when the
    first leg is polarised transversely to its momentum (``V^{mu nu} propto p^mu p^nu``)."""
    cls = (factor, factor)

    def transform(fp: Fingerprinter, smp: Samples) -> list[tuple[int, Samples]]:
        p = fp.field.p
        leg = smp.legs[0]
        k = np.asarray(smp.momenta[leg.momentum], dtype=object)            # (n, D)
        e = np.asarray(smp.polarisations[0], dtype=object)                  # (n, D, ...)
        k2 = _dot(fp, k, k)
        inv = np.array([pow(int(x), -1, p) for x in k2], dtype=object)
        ek = np.moveaxis(e, 1, -1)                                          # (n, ..., D)
        kk = k.reshape((k.shape[0],) + (1,) * (ek.ndim - 2) + (k.shape[1],))
        proj = (_dot(fp, ek, kk)[..., None] * kk * inv.reshape((-1,) + (1,) * (ek.ndim - 1))) % p
        et = np.moveaxis((ek - proj) % p, -1, 1).astype(np.int64)
        moms = np.stack([smp.momenta[l.momentum] for l in smp.legs], axis=1)
        new = fp.explicit_samples(list(cls), moms, [et] + list(smp.polarisations[1:]))
        return [(1, new)]

    return ProjectedCondition("N2", f"two-point {factor} part of DeltaS longitudinal (no transverse kinetic term)",
                              cls, transform, min_degree=1)


def antisymmetric_vertex(model: Model, factor: str) -> ProjectedCondition:
    """(N3, non-abelian factor): no part of the ``(rhotilde, A, omega)`` vertex antisymmetric in the adjoint indices
    of ``A`` and ``omega`` (the structure-constant term of ``rhotilde.D omega``)."""
    src, ghost = model.fields[factor].source, f"gh_{factor}"
    cls = tuple(sorted((src, factor, ghost), key=model.fields.order))
    ia, ig = cls.index(factor), cls.index(ghost)

    def transform(fp: Fingerprinter, smp: Samples) -> list[tuple[int, Samples]]:
        p = fp.field.p
        e = np.asarray(smp.polarisations[ia], dtype=np.int64)                # (n, D, adj)
        lor = e[:, :, 0]
        c_a = e[:, 0, :]
        c_g = np.asarray(smp.polarisations[ig], dtype=np.int64)
        moms = np.stack([smp.momenta[l.momentum] for l in smp.legs], axis=1)

        def build(ca, cg):
            pols = list(smp.polarisations)
            pols[ia] = (lor[:, :, None].astype(object) * ca[:, None, :].astype(object) % p).astype(np.int64)
            pols[ig] = cg
            return fp.explicit_samples(list(cls), moms, pols)

        return [(1, build(c_a, c_g)), (-1, build(c_g, c_a))]

    return ProjectedCondition("N3", f"no ({src}, {factor}, {ghost}) term antisymmetric in the adjoint indices",
                              cls, transform)


def _pure_matter(model: Model, cls: tuple[str, ...]) -> bool:
    return all(model.fields[n].kind in _MATTER for n in cls)


def _fermion_two_point(model: Model, cls: tuple[str, ...]) -> bool:
    return len(cls) == 2 and {model.fields[n].kind for n in cls} == {"fermion", "fermion_bar"}


@dataclass
class Normalisation:
    """Normalisation conditions in priority order (module docstring)."""

    conditions: list[Condition]
    mode: str = "full"
    fallback: bool = True

    def __post_init__(self) -> None:
        if self.mode not in ("full", "minimal"):
            raise ValueError("mode is 'full' or 'minimal'")


@dataclass
class _Group(Condition):
    """Several conditions reported under one name (e.g. the parts of N3)."""

    name: str
    description: str
    parts: list[Condition]

    def rows(self, ctx: OrderContext, basis: OperatorBasis, fp: Fingerprinter) -> ParametricMatrix:
        mats = [c.rows(ctx, basis, fp).padded(len(ctx.space)) for c in self.parts]
        mats = [m for m in mats if m.shape[0]]
        if not mats:
            return ParametricMatrix(fp.field, (0, len(basis)), {})
        return ParametricMatrix.vstack(mats)


def standard_normalisation(model: Model, *, mode: str = "full", fallback: bool = True,
                           fermion_kinetic: bool = False) -> Normalisation:
    """The report's (N1)-(N4), stated for any model (module docstring)."""
    if fermion_kinetic:
        n1: Condition = ClassCondition("N1", "no pure-matter monomial except the fermion kinetic terms",
                                       lambda m, c: _pure_matter(m, c) and not _fermion_two_point(m, c))
        firsts = {}
        for name, fm in model.fermions.items():
            acting = [f for f in model.factors if fm.reps.get(f) not in (None, "1") and fm.reps.get(f) != 0]
            if acting:
                firsts[name] = acting[0]
        vertex = {tuple(sorted((name, name + "bar", r), key=model.fields.order)) for name, r in firsts.items()}
        n1k = ClassCondition("N1k", "no fermion gauge vertex of the first gauge factor of each multiplet",
                             lambda m, c: tuple(c) in vertex)
        first = [n1, n1k]
    else:
        first = [ClassCondition("N1", "no pure-matter monomial (fields, sources and ghosts vanish)", _pure_matter)]
    two = [transverse_two_point(model, f) for f in model.factors]
    pairs = {tuple(sorted((model.fields[f].source, f"gh_{f}"), key=model.fields.order)) for f in model.factors
             if model.fields[f].source is not None and f"gh_{f}" in model.fields}
    three: list[Condition] = [ClassCondition("N3", "no gauge-source ghost two-point term rhotilde.d omega",
                                             lambda m, c: tuple(c) in pairs)]
    three += [antisymmetric_vertex(model, f) for f, fac in model.factors.items() if not fac.abelian]
    abelian = {f"gh_{f}" for f, fac in model.factors.items() if fac.abelian}
    four = ClassCondition("N4", "no abelian ghost", lambda m, c: bool(abelian & set(c)))
    conds = first + [_Group("N2", "transverse two-point gauge terms removed", two),
                     _Group("N3", "ghost-sector source terms removed", three), four]
    return Normalisation(conds, mode, fallback)


# ---------------------------------------------------------------------------------------------- solution
@dataclass
class ConditionReport:
    """``rows``: equations of the condition; ``applied``: indices of the equations imposed; ``relaxed``: True if
    the condition was relaxed (fallback), with ``dropped`` the equations not imposed."""

    name: str
    description: str
    rows: int
    applied: list[int]
    relaxed: bool = False
    dropped: list[int] = field(default_factory=list)


@dataclass
class RestoringResult:
    """``DeltaS^(l) = (1/(16 pi^2))^l sum r_alpha Q_alpha`` (module docstring)."""

    context: OrderContext
    basis: OperatorBasis
    coefficients: dict[str, Polynomial]
    conditions: list[ConditionReport]
    unknowns: int
    allowed: int
    rank_allowed: int
    rank_all: int
    existence: list[tuple[str, tuple[tuple[str, ...], ...]]]
    existence_total: int
    rows: int

    @property
    def unique(self) -> bool:
        return self.rank_allowed == self.allowed

    @property
    def kernel(self) -> int:
        """Number of closed directions left free (zero for a unique solution)."""
        return self.allowed - self.rank_allowed

    def functional(self) -> LocalFunctional:
        """``sum r_alpha Q_alpha`` (four-dimensional, symbolic in the free parameters)."""
        from .matching import _combine

        return _combine(self.context, self.basis.elements, self.coefficients)

    def nonzero(self) -> dict[str, Polynomial]:
        return {k: v for k, v in self.coefficients.items() if v}


def _consistent(a: np.ndarray, b: np.ndarray, field_: PrimeField) -> bool:
    return linalg.column_space_contains(a, b, field_)


def _stack(a: np.ndarray, c: np.ndarray, x: np.ndarray, field_: PrimeField) -> tuple[np.ndarray, np.ndarray]:
    lhs = np.concatenate([a, c]) if len(c) else a
    rhs = np.concatenate([x, np.zeros(len(c), dtype=np.int64)]) if len(c) else x
    return lhs, rhs % field_.p


def _pivot_rows(a: np.ndarray, current: np.ndarray, rows: np.ndarray, field_: PrimeField) -> list[int]:
    kernel = linalg.nullspace(np.concatenate([a, current]) if len(current) else a, field_)
    if not len(kernel) or not len(rows):
        return []
    m = rows.astype(object) @ kernel.T.astype(object) % field_.p
    return linalg.independent_columns(np.asarray(m, dtype=np.int64).T, field_)


def _b_images(ctx: OrderContext, basis: OperatorBasis) -> list[LocalFunctional]:
    """``b Q_alpha`` of a basis, cached in the context per basis object (the entry keeps the basis alive and is
    checked by identity, so a reused ``id`` of another basis never returns stale images)."""
    key = ("b", id(basis))
    hit = ctx.cached(key, lambda: (basis, [LocalFunctional.integrate(ctx.brst.b(e.functional))
                                           for e in basis.elements]))
    if hit[0] is not basis:
        hit = (basis, [LocalFunctional.integrate(ctx.brst.b(e.functional)) for e in basis.elements])
        ctx.replace_cached(key, hit)
    return hit[1]


def solve_restoring(ctx: OrderContext, data: ClassData, normalisation: Normalisation | None = None, *,
                    basis: OperatorBasis | None = None, max_primes: int = 6) -> RestoringResult:
    """Solve ``b DeltaS^(l) = -X_l`` (module docstring)."""
    basis = basis or ctx.Q
    norm = normalisation or standard_normalisation(ctx.model)
    space = ctx.space
    counts = ctx.structure_counts(1)
    images = _b_images(ctx, basis)
    n = len(basis)
    plan: dict[str, object] = {"names": basis.names}

    def prepare(field_: PrimeField, i_sign: int):
        fp = ctx.fingerprinter(field_, 4, i_sign)
        blocks = ctx.blocks(fp, counts, salt=31, physical=True)
        A = ctx.rows(fp, images, blocks)
        X = _data_rows(data, BREAKING, fp, space, blocks)
        crows = [c.rows(ctx, basis, fp).padded(len(space)) for c in norm.conditions]
        if "selection" not in plan:
            _analyse(ctx, data, norm, A, X, crows, blocks, field_, plan)
        sel = plan["selection"]
        csel = [r.rows(idx) for r, idx in zip(crows, sel) if idx]
        C = ParametricMatrix.vstack(csel) if csel else ParametricMatrix(field_, (0, n), {})
        unique = plan["unique"]

        def evaluate(point: np.ndarray) -> np.ndarray:
            lhs, rhs = _stack(A.at(point), C.at(point), (-X.at(point)[:, 0]) % field_.p, field_)
            return linalg.solve(lhs, rhs, field_, require_unique=unique)

        return evaluate

    supports = [ctx.support(e.functional) for e in basis.elements]
    polys = reconstruct(prepare, supports, space, max_primes=max_primes, seed=ctx.seed)
    reports = []
    for cond, idx, total, dropped in zip(norm.conditions, plan["selection"], plan["sizes"], plan["dropped"]):
        reports.append(ConditionReport(cond.name, cond.description, total, list(idx), bool(dropped), dropped))
    return RestoringResult(ctx, basis, {e.name: polys[k] for k, e in enumerate(basis.elements)}, reports, n,
                           plan["allowed"], plan["rank_allowed"], plan["rank_all"], plan["existence"],
                           plan["existence_total"], plan["rows"])


def _analyse(ctx: OrderContext, data: ClassData, norm: Normalisation, A: ParametricMatrix, X: ParametricMatrix,
             crows: Sequence[ParametricMatrix], blocks, field_: PrimeField, plan: dict) -> None:
    """Choose the imposed equations at a generic point (deterministic) and the rank statements."""
    p = field_.p
    space = ctx.space
    point = space.random_points(field_, 1, seed=ctx.seed + 303)[0] if len(space) else np.zeros(0, int)
    a = A.at(point)
    x = (-X.at(point)[:, 0]) % p
    n = a.shape[1]
    if not _consistent(a, x, field_):
        raise AnomalyError(match_breaking(ctx, data).anomaly)
    current = np.zeros((0, n), dtype=np.int64)
    selection, dropped_all = [], []
    for cond, rows_m in zip(norm.conditions, crows):
        rows = rows_m.at(point) if rows_m.shape[0] else np.zeros((0, n), dtype=np.int64)
        pivots = _pivot_rows(a, current, rows, field_)
        if norm.mode == "minimal":
            chosen = pivots
            dropped: list[int] = []
        else:
            trial = np.concatenate([current, rows])
            if _consistent(*_stack(a, trial, x, field_), field_):
                chosen, dropped = list(range(len(rows))), []
            elif not norm.fallback:
                raise ExistenceError(f"normalisation condition {cond.name} admits no solution")
            else:
                chosen = list(pivots)
                for k in range(len(rows)):
                    if k in chosen:
                        continue
                    attempt = np.concatenate([current, rows[sorted(chosen + [k])]])
                    if _consistent(*_stack(a, attempt, x, field_), field_):
                        chosen.append(k)
                chosen = sorted(chosen)
                dropped = [k for k in range(len(rows)) if k not in chosen]
        chosen = sorted(chosen)
        current = np.concatenate([current, rows[chosen]]) if chosen else current
        selection.append(chosen)
        dropped_all.append(dropped)
    allowed_basis = linalg.nullspace(current, field_) if len(current) else np.eye(n, dtype=np.int64)
    av = (a.astype(object) @ allowed_basis.T.astype(object) % p).astype(np.int64)
    rank_all = linalg.rank(a, field_)
    rank_allowed = linalg.rank(av, field_) if av.size else 0
    existence = _existence_conditions(a, av, current, blocks, plan["names"], field_)
    plan.update(selection=selection, dropped=dropped_all, sizes=[m.shape[0] for m in crows],
                allowed=allowed_basis.shape[0], rank_allowed=rank_allowed, rank_all=rank_all,
                existence=existence, existence_total=rank_all - rank_allowed,
                unique=rank_allowed == allowed_basis.shape[0], rows=a.shape[0])


def _existence_conditions(a: np.ndarray, av: np.ndarray, constraints: np.ndarray, blocks, names: Sequence[str],
                          field_: PrimeField) -> list[tuple[str, tuple[tuple[str, ...], ...]]]:
    """The existence conditions: a complement of ``b(allowed)`` in ``b(span Q)`` spanned by the images of
    removed directions (constraint rows, in order; a unit row is named by its basis element), each with the classes
    in which its image leaves the span of ``b(allowed)``."""
    p = field_.p
    out = []
    span = av.copy()
    for row in constraints:
        img = (a.astype(object) @ row.astype(object) % p).astype(np.int64)
        trial = np.concatenate([span, img[:, None]], axis=1) if span.size else img[:, None]
        if linalg.rank(trial, field_) == (linalg.rank(span, field_) if span.size else 0):
            continue
        nz = np.flatnonzero(row % p)
        name = names[nz[0]] if len(nz) == 1 else "+".join(names[k] for k in nz[:3]) + ("..." if len(nz) > 3 else "")
        where, start = [], 0
        for cls, smp in blocks:
            blv, bi = av[start:start + smp.n], img[start:start + smp.n]
            if np.any(bi) and not linalg.column_space_contains(blv, bi, field_):
                where.append(cls)
            start += smp.n
        out.append((name, tuple(where)))
        span = trial
    return out
