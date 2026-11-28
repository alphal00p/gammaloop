"""Local results per class of external legs: the summed ultraviolet data of all graphs of a class.

For a leg class ``(Phi_1, ..., Phi_n)`` and an order ``l`` this module enumerates the graphs
(:mod:`hepkit_renorm.diagrams`), evaluates each with a loop back-end (:mod:`hepkit_renorm.loops`) on fingerprint
samples, and sums them (report 02_method sec:m-algorithm steps 2-4 and 7-8)::

    T Gamma_{Phi_1..Phi_n}^(l)  =  (i/(16 pi^2))^L  sum_{r, k}  c[r, k]  eta^r / eps^k  + ...

``c[r, k]`` is a value on every sample (:class:`~hepkit_renorm.loops.Laurent`). At one loop the two numbers the
recipe needs are

* ``N0 = c[0, 1]`` (every class): the pole. At ghost number zero it equals ``-16 pi^2 eps V_{deltaS}`` (the vertex
  function of the counterterm with ``1/(16 pi^2 eps)`` written as unity); at ghost number one it equals
  ``-16 pi^2 eps V_{b_d deltaS}`` (the pole identity ``[(Delta . Gamma)^(1)]_div + b_d deltaS^(1) = 0``);
* ``N1 = c[1, 1]`` at ghost number one in ``D = 4`` with physical external spinors: the finite breaking ``P N1``
  (report eq:N0N1; ``X_1`` is its local functional).

The values are exact (``GF(p)``) on the samples of a :class:`~hepkit_renorm.fingerprint.Fingerprinter`, so the
same samples evaluate any local functional (:meth:`ClassResult.compare`, or a fingerprint matrix of a basis on
``ClassResult.samples`` for :mod:`hepkit_renorm.linalg`). Graphs with vertices of non-zero ``eps`` power (counterterm
insertions, ``l >= 2``) shift ``k`` by that power.

Physical external legs
----------------------
:func:`physical_samples` projects the spinor polarisation of every fermion leg on its physical chirality (report
``Pphys``: the sterile partners have no external legs): ``P_c u`` for ``Psi``, ``P_cbar`` for ``Psibar`` and the
chirality of the field for a source, as in :mod:`hepkit_renorm.classes` (``gamma5`` is symmetric in the explicit
representation, so row and column legs are projected alike). The finite breaking ``P N1`` is read on physical
samples in ``D = 4``; the poles are compared on full Dirac spinors in ``D > 4``, where the evanescent poles of
graphs with the breaking vertex on an external fermion (``Delta_1`` couples the physical to the sterile chirality)
are non-zero.

Insertions
----------
The graphs of a class of ghost number ``0`` have no insertion (``Gamma``); those of ghost number ``1`` exactly one
vertex of an insertion family (``Delta . Gamma``, linear in the breaking). Higher ghost numbers have no loop graphs
in the recipe.

Cache (:class:`ResultCache`)
----------------------------
The value of one graph on one set of samples is stored as a JSON file under a directory chosen by the user, named
by a SHA-256 key of everything it depends on: the graph record (vertex functions, ends, dimensions and orders,
propagators, symmetry factor), the back-end, the prime, the dimension, the embedding of ``i``, the parameter values
and the sample arrays. A changed rule or sample gives a new key. The back-end's code enters only through its class
name and :data:`CACHE_VERSION`, which must be bumped when a back-end changes its values.
"""

from __future__ import annotations

import hashlib
import json
import os
import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path

import numpy as np

from .classes import _spinor_chirality
from .diagrams import Enumeration, VertexFamily, enumerate_graphs, graph_propagators, store_families
from .fingerprint import Fingerprinter, Samples
from .graphs import ExternalEnd, GraphPropagator, GraphRecord
from .loops import Laurent, LoopBackend, OneLoopBackend, evaluate_graphs
from .spec import Model
from .store import OrderStore

__all__ = ["ClassComputer", "ClassResult", "GraphResult", "ResultCache", "physical_samples"]

CACHE_VERSION = "1"
"""Version of the cache layout and of the key; bumped whenever a value could change for an unchanged key."""


# ---------------------------------------------------------------------------------------------- samples
def physical_samples(fp: Fingerprinter, samples: Samples) -> Samples:
    """The samples with the spinor polarisation of every fermion leg projected on its physical chirality."""
    table = fp.model.fields
    pl, pr = fp.rep.projectors
    proj = {"L": pl.astype(object), "R": pr.astype(object)}
    p = fp.field.p
    pols = []
    for leg, pol in zip(samples.legs, samples.polarisations):
        fi = table[leg.field]
        out = pol
        for axis, slot in enumerate(fi.slots, start=1):
            if slot.space.kind == "spinor":
                m = proj[_spinor_chirality(fp.model, fi)]
                moved = np.moveaxis(out.astype(object), axis, -1)
                out = np.moveaxis(moved.dot(m.T) % p, -1, axis).astype(np.int64)
        pols.append(out)
    return Samples(samples.legs, samples.momenta, tuple(pols))


# ---------------------------------------------------------------------------------------------- results
@dataclass(frozen=True)
class GraphResult:
    """The data of one graph: its description, symmetry factor, Wick sign, ``eps`` power of its vertices and its
    Laurent coefficients on the samples (power shift applied)."""

    name: str
    symmetry_factor: Fraction
    sign: int
    power: int
    laurent: Laurent


@dataclass
class ClassResult:
    """The summed local data of the graphs of one leg class on one set of samples (module docstring)."""

    fields: tuple[str, ...]
    ghost_number: int
    order: int
    samples: Samples
    fingerprinter: Fingerprinter
    graphs: list[GraphResult] = field(default_factory=list)
    vanishing: list[str] = field(default_factory=list)
    seconds: float = 0.0
    cached: int = 0

    def by_loops(self) -> dict[int, Laurent]:
        """The sums over the graphs of each loop number ``L`` (units ``(i/(16 pi^2))^L``; at order ``l > 1`` the
        graphs with counterterm vertices have ``L < l`` and carry the remaining ``1/(16 pi^2)`` in their vertices)."""
        out: dict[int, Laurent] = {}
        for g in self.graphs:
            lau = g.laurent
            out[lau.loops] = out[lau.loops] + lau if lau.loops in out else lau
        return out

    @property
    def laurent(self) -> Laurent:
        """Sum over the graphs, all of one loop number (zero data of loop number ``order`` without graphs)."""
        parts = self.by_loops()
        if len(parts) > 1:
            raise ValueError(f"graphs with {sorted(parts)} loops: use by_loops()")
        return next(iter(parts.values())) if parts else Laurent(self.order, self.fingerprinter.field,
                                                                self.samples.n, {})

    @property
    def N0(self) -> np.ndarray:
        """One loop: the pole residue ``c[0, 1]`` summed over the graphs."""
        return self.laurent.N0

    @property
    def N1(self) -> np.ndarray:
        """One loop: the ``eta / eps`` coefficient ``c[1, 1]`` (``P N1`` in ``D = 4`` with physical legs)."""
        return self.laurent.N1

    def compare(self, values: np.ndarray, r: int = 0, k: int = 1) -> bool:
        """True if ``c[r, k]`` equals ``values`` on every sample (e.g. ``fingerprinter.functional(F, samples)``)."""
        return bool(np.array_equal(self.laurent.coefficient(r, k) % self.fingerprinter.field.p,
                                   np.asarray(values, dtype=np.int64) % self.fingerprinter.field.p))


def _shift(lau: Laurent, power: int) -> Laurent:
    if not power:
        return lau
    return Laurent(lau.loops, lau.field, lau.n, {(r, k - power): c for (r, k), c in lau.coefficients.items()})


# ---------------------------------------------------------------------------------------------- cache
def _canonical(e) -> str:
    return e.to_canonical_string()


def graph_digest(graph: GraphRecord) -> str:
    """A stable description of a graph record (everything its value depends on)."""
    parts = ["legs:" + ";".join(f"{leg.field}|{','.join(map(_canonical, leg.labels))}|{leg.momentum}"
                                for leg in graph.legs)]
    for v in graph.vertices:
        ends = ",".join(f"e{x.leg}" if isinstance(x, ExternalEnd) else f"p{x.propagator}.{x.side}" for x in v.ends)
        legs = ";".join(f"{leg.field}|{','.join(map(_canonical, leg.labels))}|{leg.momentum}" for leg in v.legs)
        parts.append(f"v:{legs}|{ends}|{sorted(v.outgoing)}|{v.dimension}|{v.order}|{_canonical(v.value)}")
    for pr in graph.propagators:
        labs = "|".join(",".join(map(_canonical, side)) for side in pr.labels)
        dens = ";".join(f"{_canonical(m)}^{n}" for m, n in pr.denominators)
        parts.append(f"p:{pr.fields}|{labs}|{_canonical(pr.numerator)}|{dens}")
    parts.append(f"s:{graph.symmetry_factor}")
    return "\n".join(parts)


def _samples_digest(fp: Fingerprinter, samples: Samples) -> str:
    h = hashlib.sha256()
    h.update(repr((fp.model.name, fp.field.p, fp.dimension, fp.rep.i_sign,
                   sorted((k, np.asarray(v).tolist() if isinstance(v, np.ndarray) else int(v))
                          for k, v in fp._values.items()))).encode())
    for leg in samples.legs:
        h.update(f"{leg.field}|{leg.momentum}".encode())
        h.update(np.ascontiguousarray(samples.momenta[leg.momentum], dtype=np.int64).tobytes())
    for pol in samples.polarisations:
        h.update(repr(pol.shape).encode())
        h.update(np.ascontiguousarray(pol, dtype=np.int64).tobytes())
    return h.hexdigest()


class ResultCache:
    """Per-graph values on disk (module docstring): ``directory/<key[:2]>/<key>.json``."""

    def __init__(self, directory: str | os.PathLike) -> None:
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def key(graph: GraphRecord, backend: LoopBackend, fp: Fingerprinter, samples: Samples) -> str:
        h = hashlib.sha256()
        for part in (CACHE_VERSION, type(backend).__qualname__, graph_digest(graph), _samples_digest(fp, samples)):
            h.update(part.encode())
            h.update(b"\0")
        return h.hexdigest()

    def _path(self, key: str) -> Path:
        return self.directory / key[:2] / f"{key}.json"

    def get(self, key: str, fp: Fingerprinter, n: int) -> Laurent | None:
        path = self._path(key)
        if not path.exists():
            return None
        doc = json.loads(path.read_text())
        if doc.get("key") != key:
            return None
        coeffs = {tuple(int(x) for x in rk.split(",")): np.array(v, dtype=np.int64)
                  for rk, v in doc["coefficients"].items()}
        return Laurent(int(doc["loops"]), fp.field, n, coeffs)

    def put(self, key: str, lau: Laurent, name: str = "") -> None:
        path = self._path(key)
        path.parent.mkdir(parents=True, exist_ok=True)
        doc = {"key": key, "name": name, "loops": lau.loops, "prime": lau.field.p,
               "coefficients": {f"{r},{k}": np.asarray(c, dtype=np.int64).tolist()
                                for (r, k), c in sorted(lau.coefficients.items())}}
        tmp = path.with_suffix(".tmp")
        tmp.write_text(json.dumps(doc))
        tmp.replace(path)


# ---------------------------------------------------------------------------------------------- driver
class ClassComputer:
    """Graphs and summed local data of leg classes at a fixed order, for one model and order store.

    ``order`` is the loop order ``l`` of sec:m-algorithm (families: :func:`~hepkit_renorm.diagrams.store_families`
    of the store, with the insertions of ``bookkeeping`` ``"A"`` or ``"B"``); ``backend`` the loop back-end
    (default: one loop); ``cache`` an optional :class:`ResultCache`; ``processes`` the number of forked workers for
    the graph evaluation.
    """

    def __init__(self, model: Model, store: OrderStore, order: int = 1, *, backend: LoopBackend | None = None,
                 cache: ResultCache | None = None, processes: int = 1,
                 families: Sequence[VertexFamily] | None = None,
                 propagators: Mapping[tuple[str, str], GraphPropagator] | None = None,
                 bookkeeping: str = "A") -> None:
        self.model, self.store, self.order = model, store, order
        self.backend = backend or OneLoopBackend()
        self.cache, self.processes = cache, processes
        self.families = tuple(families) if families is not None else \
            tuple(store_families(model, store, order, bookkeeping))
        self.propagators = dict(propagators) if propagators is not None else graph_propagators(model, store)

    def ghost_number(self, fields: Sequence[str]) -> int:
        return sum(self.model.fields[f].ghost_number for f in fields)

    def insertions(self, fields: Sequence[str]) -> int:
        """Number of insertion vertices of the graphs of a class (module docstring)."""
        g = self.ghost_number(fields)
        if g not in (0, 1):
            raise ValueError(f"classes of ghost number {g} have no loop graphs in the recipe")
        return g

    def enumerate(self, fields: Sequence[str]) -> Enumeration:
        """The graphs of a class (:func:`~hepkit_renorm.diagrams.enumerate_graphs`)."""
        return enumerate_graphs(self.model, fields, self.families, self.propagators, order=self.order,
                                insertions=self.insertions(fields))

    def records(self, fields: Sequence[str]) -> tuple[Enumeration, list[GraphRecord]]:
        """The enumeration of a class and the graph records of its topologies (same order)."""
        enum = self.enumerate(fields)
        return enum, enum.records(self.model, self.propagators)

    def samples(self, fp: Fingerprinter, fields: Sequence[str], n: int, *, salt: int = 0,
                physical: bool = False) -> Samples:
        """``n`` samples of the class; ``physical=True`` projects the external spinors (module docstring)."""
        smp = fp.samples(list(fields), n, salt)
        return physical_samples(fp, smp) if physical else smp

    def compute(self, fields: Sequence[str], fp: Fingerprinter, samples: Samples) -> ClassResult:
        """Enumerate, evaluate (with the cache) and sum the graphs of a class on ``samples``."""
        t0 = time.time()
        fields = tuple(fields)
        if tuple(leg.field for leg in samples.legs) != fields:
            raise ValueError("the samples belong to another class")
        enum, recs = self.records(fields)
        table = self.model.fields
        values: dict[int, Laurent] = {}
        todo: list[int] = []
        keys = {}
        for i, g in enumerate(recs):
            if self.cache is not None:
                keys[i] = ResultCache.key(g, self.backend, fp, samples)
                hit = self.cache.get(keys[i], fp, samples.n)
                if hit is not None:
                    values[i] = hit
                    continue
            todo.append(i)
        if todo:
            prepared = [self.backend.prepare(recs[i], self.model) for i in todo]
            for i, lau in zip(todo, evaluate_graphs(prepared, fp, {fields: samples}, backend=self.backend,
                                                    processes=self.processes)):
                values[i] = lau
                if self.cache is not None:
                    self.cache.put(keys[i], lau, recs[i].name)
        out = ClassResult(fields, self.ghost_number(fields), self.order, samples, fp,
                          vanishing=[t.describe(enum.families) for t in enum.vanishing],
                          cached=len(recs) - len(todo))
        for i, (top, g) in enumerate(zip(enum.topologies, recs)):
            out.graphs.append(GraphResult(g.name, g.symmetry_factor, g.grassmann_sign(table), top.power,
                                          _shift(values[i], top.power)))
        out.seconds = time.time() - t0
        return out
