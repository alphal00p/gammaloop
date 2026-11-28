"""Loop data of an order for the algebraic steps: every class evaluated by the diagram pipeline, resolved in the free
parameters, batched over worker processes and cached; and the local form of the data on operator bases.

The algebraic steps (:mod:`hepkit_renorm.matching`, :mod:`hepkit_renorm.st_solve`) read class data through the
:class:`~hepkit_renorm.matching.ClassData` protocol, either as values at the parameter point of a fingerprinter or as
polynomials in the free parameters (``parametric``). This module serves both from the graphs.

Parameter resolution (:class:`ResolvedLoopData`)
------------------------------------------------
A graph's value is a polynomial in the free parameters. Its support is derived from the graph, never guessed: the
products of the parameter monomials of its vertex functions (exactly, by expansion), times the monomials of the
parameters in its propagators (masses) whose mass dimension does not exceed the superficial degree ``omega`` (the
back-end keeps the terms of total small degree ``omega``, in which a propagator parameter counts by its dimension).
The graph is evaluated at ``|support|`` points of ``GF(p)^n`` (deterministic, shared by all graphs) and the
coefficients are obtained by solving the Vandermonde system exactly. Most graphs have one monomial, so the cost is
about one evaluation per graph, and the class data become polynomials (:class:`~hepkit_renorm.parametric.
ParametricMatrix`) exactly like the parametric fingerprints of local functionals.

Batching
--------
Requests ``(class, fingerprinter, samples)`` are collected and evaluated together (:meth:`ResolvedLoopData.prefetch`):
all graphs of all classes and all parameter points go into one pool of forked workers
(:func:`hepkit_renorm.loops.evaluate_jobs`), with the per-graph :class:`~hepkit_renorm.class_results.ResultCache`
on disk. A consumer that requests classes one by one can be served in one batch after a dry run that records its
requests (:class:`RequestRecorder`).

Local form (:func:`localise`)
-----------------------------
The data of a set of classes are matched once onto basis elements (the singular basis for the poles of ghost number
0, the complete breaking basis ``N`` for ``P N1`` of ghost number 1), monomial by monomial in the parameters, at one
prime (both embeddings of ``i`` for complex coefficients), and the rational coefficients are reconstructed. The
reconstruction is confirmed at a further prime: the local functional and the loop data are compared at a random
parameter point on fresh samples of every class (one evaluation per graph); if they differ, a full resolved pass at
that prime is added (CRT) and the confirmation repeated. The result (:class:`LocalForm`) is a symbolic local
functional: every later step (any prime, any samples) is served by
:class:`~hepkit_renorm.matching.FunctionalData` at the cost of fingerprints, while the data's membership in the span
of the basis is decided exactly in the fit.
"""

from __future__ import annotations

import time
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from fractions import Fraction
from itertools import pairwise

import numpy as np
from symbolica import Expression

from . import linalg
from .bases import BasisElement
from .class_results import ClassComputer, ResultCache, _samples_digest, _shift
from .exact import PrimeField
from .fingerprint import Fingerprinter, Samples
from .functional import Density, LocalFunctional
from .graphs import GraphRecord
from .loops import Laurent, PreparedGraph, evaluate_jobs
from .matching import BREAKING, POLE, MatchingError, OrderContext, _data_rows, _failing_classes
from .parametric import (
    Exponent,
    ParameterSpace,
    ParametricMatrix,
    Polynomial,
    _monomial_exponents,
    interpolate_outputs,
)

__all__ = ["LocalForm", "LoopStatistics", "Request", "RequestRecorder", "ResolvedLoopData", "graph_support",
           "localise", "recorded"]

_READOUT = {POLE: (0, 1), BREAKING: (1, 1)}
"""One loop: the Laurent coefficient ``c[r, k]`` behind each quantity (``N0``, ``N1``)."""


@dataclass(frozen=True)
class Request:
    """One class-data request: the class, the fingerprinter (prime, dimension, embedding of ``i``, parameter point)
    and the samples; ``resolved``: as polynomials in the free parameters (else at the point of ``fp``)."""

    fields: tuple[str, ...]
    fp: Fingerprinter
    samples: Samples
    resolved: bool = True


# ---------------------------------------------------------------------------------------------- supports
def _exponents(expr: Expression, space: ParameterSpace) -> set[Exponent]:
    e = space.substitute(expr).expand()
    if not len(space):
        return {()}
    if e == Expression.num(0):
        return set()
    return {_monomial_exponents(key, space.symbols) for key, _ in e.coefficient_list(*space.symbols)}


def graph_support(graph: GraphRecord, space: ParameterSpace, omega: Fraction) -> list[Exponent]:
    """The derived support of a graph's value in the free parameters (module docstring), deterministic order."""
    k = len(space)
    if not k:
        return [()]
    total: set[Exponent] = {(0,) * k}
    for v in graph.vertices:
        ev = _exponents(v.value, space)
        if not ev:
            return []
        total = {tuple(a + b for a, b in zip(x, y)) for x in total for y in ev}
    massive: set[int] = set()
    for prop in graph.propagators:
        for ex in [prop.numerator] + [m for m, _ in prop.denominators]:
            for e in _exponents(ex, space):
                massive |= {j for j, n in enumerate(e) if n}
    shifts: list[Exponent] = [(0,) * k]
    for j in sorted(massive):
        dim = space.dimensions[j]
        if dim <= 0:
            raise ValueError(f"propagator parameter {space.names[j]!r} without positive mass dimension: the "
                             "support of the graph is not bounded")
        shifts = [s[:j] + (s[j] + n,) + s[j + 1:] for s in shifts for n in range(int(omega / dim) + 1)
                  if sum(Fraction(s[i]) * space.dimensions[i] for i in range(k)) + n * dim <= omega]
    out = {tuple(a + b for a, b in zip(x, s)) for x in total for s in shifts}
    return sorted(out, key=lambda e: (sum(e), [-x for x in e]))


def _power(point: Sequence[int], e: Exponent, p: int) -> int:
    c = 1
    for x, n in zip(point, e):
        c = c * pow(int(x), int(n), p) % p
    return c


# ---------------------------------------------------------------------------------------------- loop data
@dataclass
class LoopStatistics:
    """What the loop data cost: graph evaluations by the back-end, cache hits, seconds, requests served."""

    evaluations: int = 0
    cached: int = 0
    seconds: float = 0.0
    requests: int = 0
    batches: int = 0


class ResolvedLoopData:
    """Class data from the graphs of a :class:`~hepkit_renorm.class_results.ClassComputer` (order 1: ``pole`` is
    ``N0``, ``breaking`` is ``N1``; ``P N1`` on physical samples), at the point of a fingerprinter (``values``) or as
    polynomials in the free parameters of ``space`` (``parametric``), see the module docstring.

    ``processes``: forked workers of a batch; ``cache``: optional per-graph cache; ``seed``: of the resolution
    points.
    """

    def __init__(self, computer: ClassComputer, space: ParameterSpace, *, processes: int = 1,
                 cache: ResultCache | None = None, seed: int = 0) -> None:
        if computer.order != 1:
            raise NotImplementedError("the loop readout is implemented at one loop (N0, N1); at order l > 1 the "
                                      "graphs of loop number L < l need their units (i/(16 pi^2))^(l-L) and the "
                                      "forest subtraction")
        self.computer = computer
        self.space = space
        self.processes = processes
        self.cache = cache
        self.seed = seed
        self.stats = LoopStatistics()
        self._records: dict[tuple[str, ...], tuple[object, list[GraphRecord], list[PreparedGraph],
                                                   list[list[Exponent]]]] = {}
        self._memo: dict[tuple, dict[Exponent, Laurent]] = {}
        self._points: dict[int, np.ndarray] = {}
        self._fps: dict[tuple, tuple[Fingerprinter, Fingerprinter]] = {}

    # ------------------------------------------------------------------ graphs
    def graphs(self, fields: Sequence[str]):
        """``(enumeration, records, prepared graphs, supports)`` of a class (cached)."""
        fields = tuple(fields)
        if fields not in self._records:
            enum, recs = self.computer.records(fields)
            prepared = [self.computer.backend.prepare(g, self.computer.model) for g in recs]
            supports = [graph_support(g, self.space, pg.omega) for g, pg in zip(recs, prepared)]
            self._records[fields] = (enum, recs, prepared, supports)
        return self._records[fields]

    def _point(self, field_: PrimeField, j: int) -> np.ndarray:
        if field_.p not in self._points or len(self._points[field_.p]) <= j:
            n = max(j + 1, 2 * len(self._points.get(field_.p, ())), 8)
            self._points[field_.p] = self.space.random_points(field_, n, seed=self.seed + 4099)
        return self._points[field_.p][j]

    def _at_point(self, fp: Fingerprinter, j: int) -> Fingerprinter:
        key = (id(fp), j)              # the entry keeps fp alive, so its id is not reused
        if key not in self._fps:
            self._fps[key] = (fp, fp.replace(values=self.space.values(self._point(fp.field, j))))
        return self._fps[key][1]

    @staticmethod
    def _key(req: Request) -> tuple:
        return (req.fields, req.fp.field.p, req.fp.dimension, req.fp.rep.i_sign, req.resolved,
                _samples_digest(req.fp, req.samples))

    # ------------------------------------------------------------------ evaluation
    def prefetch(self, requests: Sequence[Request]) -> None:
        """Evaluate every request not yet known, all graphs and points in one batch."""
        todo, seen = [], set()
        for req in requests:
            if tuple(leg.field for leg in req.samples.legs) != tuple(req.fields):
                raise ValueError("the samples belong to another class")
            key = self._key(req)
            if key not in self._memo and key not in seen:
                seen.add(key)
                todo.append(req)
        if not todo:
            return
        t0 = time.time()
        comp = self.computer
        plan = []                        # per request: per graph: list of job slots (one per point)
        jobs, slots, keys = [], [], []
        values: dict[int, Laurent] = {}
        for req in todo:
            _, recs, prepared, supports = self.graphs(req.fields)
            per_graph = []
            for g, pg, sup in zip(recs, prepared, supports):
                if not sup:                      # a vertex vanishes at the fixed parameters
                    per_graph.append([])
                    continue
                npts = len(sup) if req.resolved else 1
                fps = [self._at_point(req.fp, j) for j in range(npts)] if req.resolved else [req.fp]
                idx = []
                for fpj in fps:
                    slot = len(slots)
                    slots.append(None)
                    idx.append(slot)
                    key = ResultCache.key(g, comp.backend, fpj, req.samples) if self.cache is not None else None
                    keys.append(key)
                    hit = self.cache.get(key, fpj, req.samples.n) if self.cache is not None else None
                    if hit is not None:
                        values[slot] = hit
                        self.stats.cached += 1
                    else:
                        jobs.append((slot, (pg, fpj, req.samples), g.name))
                per_graph.append(idx)
            plan.append(per_graph)
        if jobs:
            results = evaluate_jobs([j for _, j, _ in jobs], backend=comp.backend, processes=self.processes)
            for (slot, _, name), lau in zip(jobs, results):
                values[slot] = lau
                if self.cache is not None:
                    self.cache.put(keys[slot], lau, name)
            self.stats.evaluations += len(jobs)
        for req, per_graph in zip(todo, plan):
            enum, recs, _, supports = self.graphs(req.fields)
            total: dict[Exponent, Laurent] = {}
            p = req.fp.field.p
            for top, sup, idx in zip(enum.topologies, supports, per_graph):
                if not idx:
                    continue
                vals = [_shift(values[s], top.power) for s in idx]
                if req.resolved:
                    pts = [self._point(req.fp.field, j) for j in range(len(sup))]
                    vander = np.array([[_power(pt, e, p) for e in sup] for pt in pts], dtype=np.int64)
                    rks = sorted({rk for v in vals for rk in v.coefficients})
                    coeffs = {}
                    for rk in rks:
                        rhs = np.stack([v.coefficient(*rk) for v in vals])
                        coeffs[rk] = linalg.solve(vander, rhs, req.fp.field, require_unique=True)
                    parts = {e: Laurent(vals[0].loops, req.fp.field, req.samples.n,
                                        {rk: np.asarray(c[i], dtype=np.int64) % p for rk, c in coeffs.items()})
                             for i, e in enumerate(sup)}
                else:
                    parts = {(): vals[0]}
                for e, lau in parts.items():
                    total[e] = total[e] + lau if e in total else lau
            self._memo[self._key(req)] = total
        self.stats.seconds += time.time() - t0
        self.stats.batches += 1

    def _data(self, req: Request) -> dict[Exponent, Laurent]:
        self.stats.requests += 1
        key = self._key(req)
        if key not in self._memo:
            self.prefetch([req])
        return self._memo[key]

    def laurent(self, fields: Sequence[str], fp: Fingerprinter, samples: Samples) -> dict[Exponent, Laurent]:
        """The summed Laurent data of a class as ``{exponents: data}`` in the free parameters."""
        return self._data(Request(tuple(fields), fp, samples, True))

    # ------------------------------------------------------------------ the ClassData protocol
    def values(self, quantity: str, fields: tuple[str, ...], fp: Fingerprinter, samples: Samples) -> np.ndarray:
        """``quantity`` at the parameter point of ``fp`` (one evaluation per graph)."""
        data = self._data(Request(tuple(fields), fp, samples, False))
        lau = data.get(())
        return lau.coefficient(*_READOUT[quantity]) % fp.field.p if lau is not None else np.zeros(samples.n, np.int64)

    def parametric(self, quantity: str, fields: tuple[str, ...], fp: Fingerprinter, space: ParameterSpace,
                   samples: Samples) -> ParametricMatrix:
        """``quantity`` as a polynomial in the free parameters (a one-column :class:`ParametricMatrix`)."""
        if tuple(space.names) != tuple(self.space.names):
            raise ValueError("the parameter space differs from the one of the loop data")
        data = self.laurent(fields, fp, samples)
        terms = {}
        for e, lau in data.items():
            col = np.asarray(lau.coefficient(*_READOUT[quantity]), dtype=np.int64) % fp.field.p
            if np.any(col):
                terms[tuple(e)] = col[:, None]
        return ParametricMatrix(fp.field, (samples.n, 1), terms)


class RequestRecorder:
    """A dry-run provider: records the requests of a consumer and answers zero (consistent with every basis), so
    that the requests can be evaluated in one batch (:meth:`ResolvedLoopData.prefetch`) before the real run."""

    def __init__(self, space: ParameterSpace) -> None:
        self.space = space
        self.requests: list[Request] = []

    def values(self, quantity: str, fields: tuple[str, ...], fp: Fingerprinter, samples: Samples) -> np.ndarray:
        self.requests.append(Request(tuple(fields), fp, samples, False))
        return np.zeros(samples.n, dtype=np.int64)

    def parametric(self, quantity: str, fields: tuple[str, ...], fp: Fingerprinter, space: ParameterSpace,
                   samples: Samples) -> ParametricMatrix:
        self.requests.append(Request(tuple(fields), fp, samples, True))
        return ParametricMatrix(fp.field, (samples.n, 1), {})


def recorded(consumer: Callable[[object], object], space: ParameterSpace) -> list[Request]:
    """The requests a consumer ``consumer(data)`` makes (dry run with :class:`RequestRecorder`)."""
    rec = RequestRecorder(space)
    consumer(rec)
    return rec.requests


# ---------------------------------------------------------------------------------------------- local form
@dataclass
class LocalForm:
    """Class data as a local functional ``F = sum_alpha c_alpha E_alpha`` with ``V_F`` equal to the data on every
    class (module docstring): ``coefficients`` are exact polynomials in the free parameters on the graded supports
    of an order-``l`` coefficient (``(re, im)`` pairs when complex), ``primes`` the primes of the full passes,
    ``confirmed`` the prime of the confirmation, ``rows`` the number of matched equations per prime and embedding."""

    quantity: str
    elements: list[BasisElement]
    coefficients: dict[str, Polynomial]
    classes: list[tuple[str, ...]]
    dimension: int
    physical: bool
    primes: list[int]
    confirmed: int
    rows: int

    def functional(self, ctx: OrderContext) -> LocalFunctional:
        """The local functional (symbolic in the free parameters)."""
        total = Density(ctx.model.fields)
        for e in self.elements:
            poly = self.coefficients.get(e.name)
            if poly:
                total = total + e.functional.scale(ctx.space.expression(poly))
        return LocalFunctional.integrate(total)

    def nonzero(self) -> dict[str, Polynomial]:
        return {k: v for k, v in self.coefficients.items() if v}


def _components(elements: Sequence[BasisElement], classes: Sequence[tuple[str, ...]]):
    """The independent blocks of the fit: classes joined by elements that span several of them; per block its
    classes and element indices."""
    parent = {c: c for c in classes}

    def find(c):
        while parent[c] != c:
            c = parent[c]
        return c

    touching = []
    for e in elements:
        cs = [tuple(c) for c in e.classes if tuple(c) in parent]
        touching.append(cs)
        for a, b in pairwise(cs):
            parent[find(a)] = find(b)
    blocks: dict[tuple, tuple[list, list]] = {}
    for c in classes:
        blocks.setdefault(find(c), ([], []))[0].append(c)
    for j, cs in enumerate(touching):
        if cs:
            blocks[find(cs[0])][1].append(j)
    return list(blocks.values())


def localise(ctx: OrderContext, data: ResolvedLoopData, quantity: str, elements: Sequence[BasisElement],
             classes: Sequence[tuple[str, ...]], *, dimension: int, physical: bool = False,
             complex_values: bool = False, salt: int = 51, check_samples: int = 2, max_primes: int = 4) -> LocalForm:
    """Match the class data onto ``elements`` (module docstring): ``F`` with ``V_F = data``. Raises
    :class:`~hepkit_renorm.matching.MatchingError` (with the failing classes) if the data of a class are not in the
    span of the elements, :class:`~hepkit_renorm.linalg.InconsistentSystemError` if a coefficient is not a
    polynomial on its graded support, and :class:`~hepkit_renorm.linalg.ReconstructionError` if no confirmed
    reconstruction is found with ``max_primes`` full passes."""
    elements = list(elements)
    names = [e.name for e in elements]
    classes = [tuple(c) for c in classes]
    per = {c: ctx.margin for c in classes}
    for e in elements:
        for c in e.classes:
            if tuple(c) in per:
                per[tuple(c)] += 1
    counts = [(c, per[c]) for c in classes]
    comps = _components(elements, classes)
    supports = [ctx.support(e.functional) for e in elements]
    offsets = np.cumsum([0] + [len(s) for s in supports])
    space = ctx.space
    signs = (1, -1) if complex_values else (1,)

    def full_pass(field_: PrimeField) -> tuple[np.ndarray, int]:
        """The images of the coefficients on their supports (real and imaginary parts stacked when complex)."""
        bls = {}
        for i_sign in signs:
            fp = ctx.fingerprinter(field_, dimension, i_sign)
            bls[i_sign] = ctx.blocks(fp, counts, salt=salt, physical=physical)
            data.prefetch([Request(c, fp, smp) for c, smp in bls[i_sign]])
        per_sign, rows = [], 0
        for i_sign in signs:
            fp = ctx.fingerprinter(field_, dimension, i_sign)
            bl = bls[i_sign]
            A = ctx.rows(fp, [e.functional for e in elements], bl)
            R = _data_rows(data, quantity, fp, space, bl).padded(len(space))
            start, where = 0, {}
            for c, smp in bl:
                where[c] = list(range(start, start + smp.n))
                start += smp.n
            rows += start

            def evaluate(point: np.ndarray, A=A, R=R, where=where) -> np.ndarray:
                a, r = A.at(point), R.at(point)[:, 0]
                x = np.zeros(len(elements), dtype=np.int64)
                for cls, cols in comps:
                    idx = [k for c in cls for k in where[c]]
                    sub, rhs = a[np.ix_(idx, cols)], r[idx]
                    try:
                        if cols:
                            x[cols] = linalg.solve(sub, rhs[:, None], field_, require_unique=True)[:, 0]
                        elif np.any(rhs % field_.p):
                            raise linalg.InconsistentSystemError("no elements")
                    except linalg.InconsistentSystemError:
                        layout = [(c, len(where[c])) for c in cls]
                        bad = _failing_classes(sub, rhs[:, None], field_, layout) if cols else \
                            [c for c in cls if np.any(r[where[c]] % field_.p)]
                        raise MatchingError(f"{quantity} data are not in the span of the basis", bad) from None
                return x

            per_sign.append(interpolate_outputs(evaluate, supports, space, field_, seed=ctx.seed + 4111))
        if not complex_values:
            return np.asarray(per_sign[0], dtype=object), rows
        re, im = linalg.split_gaussian(per_sign[0].astype(np.int64), per_sign[1].astype(np.int64), field_)
        return np.concatenate([re, im]).astype(object), rows

    def polynomials(rat: np.ndarray) -> dict[str, Polynomial]:
        n = int(offsets[-1])
        out: dict[str, Polynomial] = {}
        for j, name in enumerate(names):
            poly: Polynomial = {}
            for k, e in enumerate(supports[j]):
                pos = int(offsets[j]) + k
                c = (Fraction(rat[pos]), Fraction(rat[n + pos])) if complex_values else Fraction(rat[pos])
                if c not in (0, (0, 0)):
                    poly[tuple(e)] = c
            out[name] = poly
        return out

    def confirm(field_: PrimeField, coeffs: dict[str, Polynomial]) -> bool:
        form = LocalForm(quantity, elements, coeffs, classes, dimension, physical, [], 0, 0)
        F = form.functional(ctx)
        fp = ctx.fingerprinter(field_, dimension, 1)
        bl = ctx.blocks(fp, [(c, check_samples) for c in classes], salt=salt + 1, physical=physical)
        data.prefetch([Request(c, fp, smp, False) for c, smp in bl])
        return all(np.array_equal(np.asarray(fp.functional(F, smp)) % field_.p,
                                  data.values(quantity, c, fp, smp) % field_.p) for c, smp in bl)

    values, moduli = [], []
    rows = 0
    primes = linalg.prime_fields()
    f = next(primes)
    for _ in range(max_primes):
        v, rows = full_pass(f)
        values.append(v)
        moduli.append(f.p)
        rat = linalg.reconstruct_all(values, moduli)
        f = next(primes)
        if rat is not None and confirm(f, coeffs := polynomials(rat)):
            return LocalForm(quantity, elements, coeffs, classes, dimension, physical, moduli, f.p, rows)
    raise linalg.ReconstructionError(f"no confirmed reconstruction of the {quantity} data with {max_primes} primes")
