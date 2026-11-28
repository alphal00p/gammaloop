"""Local ultraviolet data of Feynman graphs: the loop back-end interface and the one-loop back-end.

The interface (:class:`LoopBackend`) takes a :class:`~hepkit_renorm.graphs.GraphRecord` and returns the local
ultraviolet data of the graph as exact values on fingerprint samples (:class:`Laurent`), so that the driver can
compare a sum of graphs with the vertex functions of local operators through :mod:`hepkit_renorm.fingerprint`
and decide every equality and every solve in :mod:`hepkit_renorm.linalg`. Nothing is returned in a symbolic
canonical form. Any back-end (this module's one-loop back-end, later an ``l``-loop back-end on Idenso and the Rust
extraction path) must return the same object.

What is returned (report 03_setup sec:integrals, eq:N0N1, eq:extraction; 02_method sec:m-K)
-------------------------------------------------------------------------------------------
For a graph ``G`` with ``L`` loops (value as in :mod:`hepkit_renorm.graphs`: Feynman rules ``i V``, propagators,
symmetry factor, Grassmann sign; amputated, legs ``Phi_1 .. Phi_n``)::

    T G  =  (i / (16 pi^2))^L  sum_{r >= 0, k >= 1}  c[r, k]  eta^r / eps^k  +  (finite, M_UV-dependent)

``d = 4 - 2 eps``; ``eta`` is the evanescent tag. ``c[r, k]`` is a polynomial in the external momenta and masses,
with the complete ``d``-dimensional tensor structure of the external legs (barred and hatted momenta, metrics and
gamma matrices of the legs are kept). :class:`Laurent` holds ``c[r, k]`` evaluated on samples: at the momenta and
polarisations of :class:`~hepkit_renorm.fingerprint.Samples` (in any explicit dimension ``D = 4, 6, 8, ...``,
contracted with the polarisations exactly like :meth:`~hepkit_renorm.fingerprint.Fingerprinter.evaluate` contracts
a vertex function), with the parameter values, prime and embedding of ``i`` of the
:class:`~hepkit_renorm.fingerprint.Fingerprinter`. At one loop ``N0 = c[0, 1]`` and ``N1 = c[1, 1]``
(eq:N0N1); the pole operators of sec:m-K act on the map (:meth:`Laurent.K`, :meth:`Laurent.Kf`, :meth:`Laurent.Ke`).
So ``-K T G`` is ``Z(G)``, the contribution of ``G`` to the Feynman rule ``i delta S`` of the counterterm
(eq:m-dS-per-graph): compare ``sum_G c[0, 1]`` with ``-16 pi^2 eps V_{delta S}``, i.e. the vertex function of a
local operator with the prefactor ``1/(16 pi^2 eps)`` written as unity.

The tag (sec:m-K, "the tag at any loop order")
-----------------------------------------------
The tensor reduction is done in ``d`` dimensions with the full metric; its coefficients and the vacuum integrals
depend on ``d = 4 - 2 eps`` and are never tagged. A Lorentz index that is summed inside the graph (loop-momentum
pairings, internal propagators, contracted gamma matrices) runs over the four barred directions and the ``d - 4``
hatted ones; every summation over the hatted directions contributes ``d - 4 = -2 eta``. The value is therefore a
polynomial in the internal hatted dimension ``dh``, and ``c[r, k]`` is the coefficient of ``(-2 eta)^r``. This is the
report's rule for closed cycles of hatted metrics (``ghat^mu_mu = -2 eta``, ``khat^2 -> -2 eta k^2 / d``) applied to
every internal summation, also to the unsplit ``d``-dimensional ones (``g^mu_mu = 4 - 2 eta``, ``k^2`` in a numerator
``= kbar^2 + khat^2``). The two uses of the tag at one loop do not depend on how the unsplit objects are treated
(03_setup, after eq:eta): ``N0`` (``eta`` set to zero) and the physical projection ``P N1`` of the graphs with a
ghost-number-one vertex, where every term carries a hatted object. ``N1`` of graphs without such a vertex depends
on the convention (the toolkit's README, item 4, and the report's statement ``N1 = 0`` for them assumes unsplit
bosonic contractions); it is not used at one loop.

Exact evaluation
----------------
External data live in ``D_ext`` explicit dimensions (the samples). The graph is evaluated with every internal index
in ``D_int = D_ext + m`` dimensions, ``m = 0, 1, ..., B`` (external vectors padded with zeros, external spinors
embedded as ``u (x) e_0``, so that ``ubar M u`` is unchanged by the extra directions), which gives the value at
``dh = D_ext - 4 + m``. Odd ``D_int`` use the gammas of ``D_int - 1`` dimensions plus their product as the last one
(:class:`_OddRepresentation`): an extra direction is reached only through internal summations, hence an even number
of times in every word, on which this representation is faithful; the spinor dimension then grows only every second
step. The value is a polynomial in ``dh`` of degree at most ``B``, one factor per internal index line that can run
over the extra directions; per combination of parts ``B`` = half the number of loop-momentum slots that can lie on such
a line (not on an external momentum, ``gbar``, ``gammabar`` or ``eps``), plus the internal propagators with two
Lorentz ends, plus the closed lines inside single vertices. ``B + 1`` values determine the polynomial exactly
(:func:`~hepkit_renorm.linalg.solve` on the Vandermonde system); each combination is evaluated only in the dimensions
its own bound needs. Everything is over ``GF(p)``.

The one-loop back-end (:class:`OneLoopBackend`)
------------------------------------------------
1. Integrand. The vertex functions and propagator numerators are relabelled (vertex legs take the labels of the
   external legs or of the propagator ends they are joined to, dummies become graph-unique), the momenta of the
   routing (:meth:`GraphRecord.routing`) are substituted, and ``d -> 4`` in coefficients (only the pole is needed).
2. ``T^(omega)`` (eq:uv-expansion) with the expansion point ``k^2``: every denominator is expanded,
   ``1/((k+q)^2 - m^2)^n = sum_j binom(-n, j) (2 k.q + q^2 - m^2)^j / (k^2)^(n+j)``. A term contributes to the pole
   iff its vacuum integral is logarithmic, i.e. iff its degree in the small quantities (external momenta and masses,
   counted by mass dimension) equals ``omega(G)`` (:meth:`GraphRecord.superficial_degree`); terms of lower degree
   (power divergences) are proportional to the auxiliary mass ``M_UV`` and cancel in the pole (sec:integrals), so the
   pole is computed from the logarithmic terms alone. Each factor is split into parts of fixed degree in ``k`` (the
   small degree of a part follows from the mass dimension of the factor), and the parts are combined.
3. Tensor reduction (eq:tensor-red): ``k^{mu_1} .. k^{mu_2n} -> (k^2)^n / (d (d+2) .. (d+2n-2)) sum_pairings g g ..``;
   each pairing of the loop-momentum slots is a full metric contraction of the network.
4. Vacuum integral (eq:logpole): ``int (k^2)^n / (k^2)^(n+2) -> i/(16 pi^2) / eps``; the pole coefficient is the
   numerator coefficient at ``d = 4`` (``1 / (4 . 6 ... (2n + 2))`` from the reduction).
5. Numerics. Every part is evaluated once per dimension as a dense tensor in its Lorentz and spinor labels (on the
   samples when it depends on external momenta), times a separate tensor in its group and flavour labels, with the
   head arrays of :class:`~hepkit_renorm.fingerprint.Fingerprinter`; the network of a combination of parts and a
   pairing is contracted exactly; closed spinor cycles are normalised to ``tr 1 = 4`` as in the fingerprints.

Higher loops
------------
:class:`Laurent` keeps every ``(eta^r, eps^-k)`` coefficient and the loop order; :meth:`GraphRecord.loops`, the
vertex ``order`` and the forest of divergent subgraphs (:func:`divergent_subgraphs`) are part of the input, and an
``l``-loop back-end implements the same :class:`LoopBackend` methods. The one-loop back-end refuses graphs with
``L != 1``.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from fractions import Fraction
from functools import cached_property
from typing import Protocol

import numpy as np
from symbolica import AtomType, Expression, Replacement, S

from . import linalg as la
from . import tensors as tn
from .exact import PrimeField
from .fields import FieldTable
from .fingerprint import Fingerprinter, Representation, Samples, _closed_spinor_cycles, _pair
from .graphs import LOOP_MOMENTUM, ExternalEnd, GraphRecord, PropagatorEnd
from .rules import PROPAGATOR_MOMENTUM
from .spec import Model

__all__ = [
    "Laurent",
    "LoopBackend",
    "OneLoopBackend",
    "PreparedGraph",
    "divergent_subgraphs",
    "embed_samples",
    "evaluate_graphs",
    "evaluate_jobs",
    "local_uv",
]

_LETTERS = list("abcdefghijklmnopqrstuvwxyABCDEFGHIJKLMNOPQRSTUVWXYZ")
_BATCH = "#batch"


# ---------------------------------------------------------------------------------------------- results
@dataclass(frozen=True)
class Laurent:
    """Local ultraviolet data on samples: ``T G = (i/(16 pi^2))^loops sum c[r, k] eta^r / eps^k + ...``.

    ``coefficients[(r, k)]`` are arrays of values (one per sample) in ``GF(p)`` (module docstring). Graph values
    add: ``a + b`` is the sum of the two graphs (same samples and loop order).
    """

    loops: int
    field: PrimeField
    n: int
    coefficients: Mapping[tuple[int, int], np.ndarray]

    def coefficient(self, r: int, k: int) -> np.ndarray:
        """``c[r, k]`` (zeros if absent)."""
        c = self.coefficients.get((r, k))
        return c if c is not None else np.zeros(self.n, dtype=np.int64)

    @property
    def N0(self) -> np.ndarray:
        """One loop: ``N0 = c[0, 1]``, the residue of the pole (eq:N0N1)."""
        self._one_loop()
        return self.coefficient(0, 1)

    @property
    def N1(self) -> np.ndarray:
        """One loop: ``N1 = c[1, 1]``, the coefficient of ``eta / eps`` (eq:N0N1)."""
        self._one_loop()
        return self.coefficient(1, 1)

    def _one_loop(self) -> None:
        if self.loops != 1:
            raise ValueError("N0 and N1 are one-loop names; use coefficient(r, k)")

    def K(self) -> dict[int, np.ndarray]:
        """``K``: ``eta = eps`` first, then the pole part: ``{k - r: sum c[r, k]}`` over ``r < k``."""
        out: dict[int, np.ndarray] = {}
        for (r, k), c in self.coefficients.items():
            if r < k:
                out[k - r] = (out.get(k - r, 0) + c) % self.field.p
        return out

    def Kf(self) -> dict[tuple[int, int], np.ndarray]:
        """``K_f``: the pole part at fixed ``eta``, ``eta`` kept: ``{(r, k): c[r, k]}`` for ``k > 0``."""
        return {rk: c for rk, c in self.coefficients.items() if rk[1] > 0}

    def Ke(self) -> dict[int, np.ndarray]:
        """``K_eta = [K_f]_{eta = eps}``: ``{k - r: sum c[r, k]}`` over all ``r`` (keys ``<= 0`` are finite and
        positive powers of ``eps``)."""
        out: dict[int, np.ndarray] = {}
        for (r, k), c in self.Kf().items():
            out[k - r] = (out.get(k - r, 0) + c) % self.field.p
        return out

    def __add__(self, other: Laurent) -> Laurent:
        if other.loops != self.loops or other.field.p != self.field.p or other.n != self.n:
            raise ValueError("incompatible Laurent data")
        keys = set(self.coefficients) | set(other.coefficients)
        p = self.field.p
        return Laurent(self.loops, self.field, self.n,
                       {k: (self.coefficient(*k) + other.coefficient(*k)) % p for k in keys})

    def scale(self, factor: int | Fraction) -> Laurent:
        """The data times a rational number."""
        f = self.field(factor)
        return Laurent(self.loops, self.field, self.n, {k: c * f % self.field.p for k, c in self.coefficients.items()})


@dataclass
class PreparedGraph:
    """The dimension-independent part of a back-end's work on one graph (symbolic integrand, parts, combinations).

    ``data`` is private to the back-end that prepared it; ``degree`` bounds the polynomial degree in the internal
    hatted dimension (module docstring).
    """

    graph: GraphRecord
    loops: int
    omega: Fraction
    degree: int
    data: object = field(repr=False, default=None)


class LoopBackend(Protocol):
    """The interface of a loop back-end (module docstring)."""

    def supports(self, graph: GraphRecord) -> bool:
        """True if the back-end can treat the graph (loop order, forest)."""

    def prepare(self, graph: GraphRecord, model: Model) -> PreparedGraph:
        """Symbolic, dimension-independent preparation (can be reused for any samples and dimension)."""

    def evaluate(self, prepared: PreparedGraph, fingerprinter: Fingerprinter, samples: Samples) -> Laurent:
        """The local ultraviolet data on the samples (in ``fingerprinter.dimension`` dimensions)."""


def local_uv(graph: GraphRecord, model: Model, fingerprinter: Fingerprinter, samples: Samples,
             backend: LoopBackend | None = None) -> Laurent:
    """``backend.evaluate(backend.prepare(graph, model), fingerprinter, samples)`` (default: one loop)."""
    be = backend or OneLoopBackend()
    if not be.supports(graph):
        raise NotImplementedError(f"the back-end does not support the graph {graph.name or ''}")
    return be.evaluate(be.prepare(graph, model), fingerprinter, samples)


def divergent_subgraphs(graph: GraphRecord, table: FieldTable) -> list[frozenset[int]]:
    """Proper 1PI subgraphs with loops and ``omega >= 0`` (the elements of the forests of eq:m-Rbar), as sets of
    propagator indices; each subgraph is the set of propagators together with the vertices they touch."""
    nprop = len(graph.propagators)
    out = []
    ends = [(graph.attachment(PropagatorEnd(e, 0))[0], graph.attachment(PropagatorEnd(e, 1))[0])
            for e in range(nprop)]
    for size in range(1, nprop):
        for sub in itertools.combinations(range(nprop), size):
            verts = sorted({v for e in sub for v in ends[e]})
            loops = len(sub) - len(verts) + 1
            if loops < 1 or not _connected_1pi(sub, ends, verts):
                continue
            omega = Fraction(4 * loops)
            for v in verts:
                vert = graph.vertices[v]
                omega += vert.dimension - sum(table[f].dimension for f in vert.fields)
            for e in sub:
                prop = graph.propagators[e]
                omega += table[prop.fields[0]].dimension + table[prop.fields[1]].dimension - 4
            if omega >= 0:
                out.append(frozenset(sub))
    return out


def _connected_1pi(sub: Sequence[int], ends: Sequence[tuple[int, int]], verts: Sequence[int]) -> bool:
    for cut in [None, *sub]:
        edges = [ends[e] for e in sub if e != cut]
        seen, stack = {verts[0]}, [verts[0]]
        while stack:
            v = stack.pop()
            for a, b in edges:
                for x, y in ((a, b), (b, a)):
                    if x == v and y not in seen:
                        seen.add(y)
                        stack.append(y)
        if len(seen) != len(verts):
            return False
    return True


# ---------------------------------------------------------------------------------------------- samples
def embed_samples(samples: Samples, source: Fingerprinter, target: Fingerprinter) -> Samples:
    """The samples of ``source`` (dimension ``D``) in the larger dimension of ``target``: momenta and Lorentz
    polarisations padded with zeros, spinor polarisations ``u -> u (x) e_0`` (module docstring)."""
    d_src, d_tgt = source.dimension, target.dimension
    if d_tgt == d_src:
        return samples
    if d_tgt < d_src:
        raise ValueError("can only embed into a larger dimension")
    ratio = target.rep.spinor_dim // source.rep.spinor_dim
    moms = {k: np.pad(v, ((0, 0), (0, d_tgt - d_src))) for k, v in samples.momenta.items()}
    pols = []
    for leg, pol in zip(samples.legs, samples.polarisations):
        out = pol
        for axis, slot in enumerate(source.model.fields[leg.field].slots, start=1):
            if slot.space.kind == "lorentz":
                width = [(0, 0)] * out.ndim
                width[axis] = (0, d_tgt - d_src)
                out = np.pad(out, width)
            elif slot.space.kind == "spinor":
                shape = list(out.shape)
                new = np.zeros(shape[:axis] + [shape[axis] * ratio] + shape[axis + 1:], dtype=np.int64)
                index = [slice(None)] * out.ndim
                index[axis] = slice(0, shape[axis] * ratio, ratio)
                new[tuple(index)] = out
                out = new
        pols.append(out)
    return Samples(samples.legs, moms, tuple(pols))


# ---------------------------------------------------------------------------------------------- dense evaluation
def _key(label: Expression) -> str:
    return label.to_canonical_string()


def _factors(term: Expression) -> list[Expression]:
    """The factors of a product term, powers of functions repeated."""
    out = []
    for fa in (list(term) if term.get_type() == AtomType.Mul else [term]):
        if fa.get_type() == AtomType.Pow:
            base, ex = list(fa)
            if ex.get_type() == AtomType.Num and base.get_type() == AtomType.Fn:
                out.extend([base] * int(ex.to_canonical_string()))
                continue
        out.append(fa)
    return out


def _terms(expr: Expression) -> list[Expression]:
    e = expr.expand()
    if e.get_type() == AtomType.Add:
        return list(e)
    if e.get_type() == AtomType.Num and e.to_canonical_string() == "0":
        return []
    return [e]


def _contract(operands: Sequence[tuple[np.ndarray, Sequence[str]]], output: Sequence[str], p: int,
              paths: dict | None = None) -> np.ndarray:
    """Exact contraction of a labelled network over ``GF(p)`` (greedy pairwise order, cached per structure in
    ``paths``; every pairwise step is the limb-split product of :mod:`hepkit_renorm.fingerprint`)."""
    ops = [(a, list(labs)) for a, labs in operands]
    out = list(output)
    signature = (tuple(tuple(labs) for _, labs in ops), tuple(out))
    path = paths.get(signature) if paths is not None else None
    if path is None:
        sizes = {lab: n for a, labs in ops for lab, n in zip(labs, a.shape)}
        path = []
        cur = [list(labs) for _, labs in ops]
        while len(cur) > 1:
            best = None
            for i in range(len(cur)):
                for j in range(i + 1, len(cur)):
                    rest = {lab for k, labs in enumerate(cur) if k not in (i, j) for lab in labs} | set(out)
                    union = list(dict.fromkeys(cur[i] + cur[j]))
                    keep = [lab for lab in union if lab in rest]
                    size = math.prod(sizes[lab] for lab in keep)
                    cost = math.prod(sizes[lab] for lab in union)
                    score = (not (set(cur[i]) & set(cur[j])), cost - math.prod(sizes[lab] for lab in cur[i])
                             - math.prod(sizes[lab] for lab in cur[j]) + size, size)
                    if best is None or score < best[0]:
                        best = (score, i, j, keep)
            _, i, j, keep = best
            path.append((i, j, tuple(keep)))
            cur = [labs for k, labs in enumerate(cur) if k not in (i, j)] + [keep]
        if paths is not None:
            paths[signature] = path
    for i, j, keep in path:
        (a, la_), (b, lb) = ops[i], ops[j]
        ops = [op for k, op in enumerate(ops) if k not in (i, j)] + [(_pairwise(a, la_, b, lb, keep, p), list(keep))]
    a, labs = ops[0]
    letters = {lab: _LETTERS[n] for n, lab in enumerate(dict.fromkeys(labs + out))}
    return np.einsum("".join(letters[x] for x in labs) + "->" + "".join(letters[x] for x in out), a) % p


def _pairwise(a: np.ndarray, la_: Sequence[str], b: np.ndarray, lb: Sequence[str], keep: Sequence[str], p: int
              ) -> np.ndarray:
    """One exact pairwise contraction; summations longer than the int64 bound of :mod:`fingerprint` are split into
    slices of a summed label."""
    summed = [lab for lab in dict.fromkeys(list(la_) + list(lb)) if lab not in keep]
    sizes = {lab: n for labs, x in ((la_, a), (lb, b)) for lab, n in zip(labs, x.shape)}
    if math.prod(sizes[lab] for lab in summed) >= 1 << 16:
        lab = max(summed, key=lambda x: sizes[x])
        total = None
        for idx in range(sizes[lab]):
            sa = a[tuple(idx if x == lab else slice(None) for x in la_)] if lab in la_ else a
            sb = b[tuple(idx if x == lab else slice(None) for x in lb)] if lab in lb else b
            na = [x for x in la_ if x != lab] if lab in la_ else list(la_)
            nb = [x for x in lb if x != lab] if lab in lb else list(lb)
            part = _pairwise(sa, na, sb, nb, keep, p)
            total = part if total is None else (total + part) % p
        return total
    letters = {lab: _LETTERS[n] for n, lab in enumerate(dict.fromkeys(list(la_) + list(lb)))}
    c, _ = _pair(a, "".join(letters[x] for x in la_), b, "".join(letters[x] for x in lb),
                 {letters[x] for x in keep}, p)
    return c


class _Dense:
    """Dense tensors of expressions over their open labels, with the head arrays of a
    :class:`~hepkit_renorm.fingerprint.Fingerprinter` (one dimension, one set of samples)."""

    def __init__(self, fp: Fingerprinter, samples: Samples) -> None:
        self.fp, self.samples = fp, samples
        self.p = fp.field.p
        self.paths: dict = {}

    def dim(self, label: Expression) -> int:
        kind, _, fac = tn.label_space(label).partition(":")
        return self.fp._space_dim(kind, fac or None)

    def tensor(self, expr: Expression, open_labels: Sequence[Expression]) -> tuple[np.ndarray, bool]:
        """``(array, batched)``: the tensor of ``expr`` over ``open_labels`` (contravariant components; labels that
        occur twice are contracted with the metric), with a leading sample axis if ``batched``."""
        keys = [_key(lab) for lab in open_labels]
        dims = tuple(self.dim(lab) for lab in open_labels)
        values: list[tuple[np.ndarray, bool]] = []
        for term in _terms(expr):
            coeff, ops, b = self._term(term, keys)
            if coeff == 0:
                continue
            if not ops:
                values.append((np.full(dims, coeff, dtype=np.int64), False))
                continue
            val = _contract(ops, ([_BATCH] if b else []) + keys, self.p, self.paths)
            values.append((val * coeff % self.p, b))
        batched = any(b for _, b in values)
        shape = ((self.samples.n,) if batched else ()) + dims
        total = np.zeros(shape, dtype=np.int64)
        for val, _ in values:
            total = (total + val) % self.p
        return total, batched

    def _term(self, term: Expression, open_keys: Sequence[str]) -> tuple[int, list, bool]:
        fp, p = self.fp, self.p
        coeff = 1
        raw: list[tuple[str, list[str], object]] = []
        label_space: dict[str, str] = {}

        def labels_of(args: Sequence[Expression]) -> list[str]:
            out = []
            for a in args:
                k = _key(a)
                label_space[k] = tn.label_space(a)
                out.append(k)
            return out

        for fa in (list(term) if term.get_type() == AtomType.Mul else [term]):
            reps, base = 1, fa
            if fa.get_type() == AtomType.Pow:
                base, ex = list(fa)
                reps = int(ex.to_canonical_string())
            t = base.get_type()
            if t == AtomType.Num:
                coeff = coeff * pow(fp.rep.number(base), reps, p) % p
            elif t == AtomType.Var or (t == AtomType.Fn and base.get_name() == "symbolica::conj"
                                       and next(iter(base)).get_type() == AtomType.Var):
                coeff = coeff * pow(fp._scalar(base), reps, p) % p
            elif t == AtomType.Fn:
                for _ in range(reps):
                    raw.append(fp._head(base, labels_of))
            else:
                raise ValueError(f"cannot evaluate {fa}")
        counts: dict[str, int] = {}
        for _, keys, _ in raw:
            for k in keys:
                counts[k] = counts.get(k, 0) + 1
        opened = set(open_keys)
        for k in open_keys:
            if counts.get(k, 0) != 1:
                raise ValueError(f"open label {k} occurs {counts.get(k, 0)} times in {term}")
        bad = [k for k, c in counts.items() if k not in opened and c != 2]
        if bad:
            raise ValueError(f"labels {bad} are neither open nor contracted in {term}")
        ops: list[tuple[np.ndarray, list[str]]] = []
        batched = False
        for kind, keys, data in raw:
            if kind == "const":
                ops.append((data, keys))
            elif kind == "mom":
                ops.append((self.samples.momenta[data], [_BATCH] + keys))
                batched = True
            else:
                name, k, l_ = data
                mask = {"sp": 1, "spbar": fp.rep.bar, "sphat": fp.rep.hat}[name]
                w = fp.rep.metric * mask % p
                ops.append((_contract([(self.samples.momenta[k], [_BATCH, "m"]), (self.samples.momenta[l_],
                                                                                   [_BATCH, "m"]), (w, ["m"])],
                                      [_BATCH], p), [_BATCH]))
                batched = True
        for k, sp in label_space.items():
            if sp == "lorentz" and k not in opened:
                ops.append((fp.rep.metric, [k]))
        closed = _closed_spinor_cycles(raw, label_space, opened)
        coeff = coeff * pow(fp.rep.trace_norm, closed, p) % p
        return coeff, ops, batched


# ---------------------------------------------------------------------------------------------- one-loop integrand
_LABEL_BASE = 50_000_000
_SLOT_BASE = 60_000_000
_x_ = S("x_")


class _OddRepresentation(Representation):
    """Explicit tensors in an odd dimension ``D = 2K + 1``: the gammas of ``D - 1`` dimensions and ``Gamma^(D-1) = c
    Gamma^0 ... Gamma^(D-2)`` with ``c`` fixed by ``(Gamma^(D-1))^2 = -1``; spinor dimension that of ``D - 1``.

    Odd dimensions are exact for the graphs: an extra direction is reached only through internal summations, so it
    occurs an even number of times in every word, and on such words the representation is faithful (the relation
    ``Gamma^0 ... Gamma^(D-1) ~ 1`` of odd dimensions never applies)."""

    def __post_init__(self) -> None:
        if self.dimension < 5 or self.dimension % 2 == 0:
            raise ValueError("odd dimension >= 5 expected")

    @property
    def hat_pairs(self) -> int:
        return (self.dimension - 5) // 2

    @cached_property
    def gamma(self) -> np.ndarray:
        p = self.p
        even = Representation(self.field, self.dimension - 1, self.i_sign).gamma.astype(object)
        prod = np.eye(even.shape[1], dtype=object)
        for gm in even:
            prod = prod.dot(gm) % p
        square = int(prod.dot(prod)[0, 0] % p)
        c = 1 if square == p - 1 else self.i
        last = prod * c % p
        return np.concatenate([even, last[None]]).astype(np.int64)


_FINGERPRINTERS: dict[tuple, tuple[Fingerprinter, Fingerprinter]] = {}


def _with_dimension(fp: Fingerprinter, dimension: int) -> Fingerprinter:
    """``fp`` in another explicit dimension (odd dimensions through :class:`_OddRepresentation`), with the same
    prime, seed, parameter values and embedding of ``i``; cached."""
    key = (id(fp), dimension)
    hit = _FINGERPRINTERS.get(key)
    if hit is None or hit[0] is not fp:
        if dimension % 2 == 0:
            new = fp.replace(dimension=dimension)
        else:
            new = fp.replace(dimension=dimension - 1)
            new.rep = _OddRepresentation(fp.field, dimension, fp.rep.i_sign)
        hit = (fp, new)
        _FINGERPRINTERS[key] = hit
    return hit[1]


def _relabel(label: Expression, n: int) -> Expression:
    args = list(label)
    return S(label.get_name())(*args[:-1], n)


def _num(c: Fraction | int) -> Expression:
    c = Fraction(c)
    return Expression.num(c.numerator) / c.denominator


def _combination(mom: Mapping[int, Fraction], head, *rest) -> Expression:
    out = Expression.num(0)
    for n, c in mom.items():
        out = out + _num(c) * head(n, *rest)
    return out


def _substitute_momenta(expr: Expression, mapping: Mapping[int, Mapping[int, Fraction]]) -> Expression:
    """Simultaneous substitution ``p(m) -> sum c_n p(n)`` (and bilinearly in ``sp``, ``spbar``, ``sphat``)."""
    reps = [(tn.p(m, _x_), _combination(mom, tn.p, _x_)) for m, mom in mapping.items()]
    nums = list(mapping)
    for i, a in enumerate(nums):
        for b in nums[i:]:
            for head in (tn.sp, tn.spbar, tn.sphat):
                val = Expression.num(0)
                for n1, c1 in mapping[a].items():
                    for n2, c2 in mapping[b].items():
                        val = val + _num(c1 * c2) * head(n1, n2)
                reps.append((head(a, b), val))
    return expr.replace_multiple([Replacement(a, b) for a, b in reps]).expand()


def _is_ld(fa: Expression) -> bool:
    """A factor of the Lorentz-Dirac part (numbers, scalars, momenta, Lorentz and spinor heads)."""
    if fa.get_type() != AtomType.Fn:
        return True
    labs = tn.labels_in(fa)
    return all(tn.label_space(lab) in ("lorentz", "spinor") for lab in labs)


@dataclass
class _Piece:
    group: Expression
    group_labels: tuple[Expression, ...]
    ld: Expression
    ld_labels: tuple[Expression, ...]


@dataclass
class _Part:
    small: int                 # degree in the small quantities
    slots: tuple[Expression, ...]
    kpow: int                  # power of 1/k^2 contributed (denominator parts)
    pieces: list[_Piece]
    live: int = 0              # slots that can lie on an index line over the hatted directions
    inner: int = 0             # closed index lines inside one term (over the pieces)

    def __post_init__(self) -> None:
        keys = {_key(s) for s in self.slots}
        live: set[str] = set()
        for piece in self.pieces:
            for term in _terms(piece.ld):
                find, dead = _lorentz_lines(term)
                live |= {k for k in keys if find(k) not in dead}
            self.inner = max(self.inner, _inner_lorentz(piece.ld, piece.ld_labels))
        self.live = len(live)


@dataclass
class _Factor:
    name: str
    labels: tuple[Expression, ...]      # open labels except loop-momentum slots
    lowered: frozenset[str]             # Lorentz labels whose metric this factor carries
    parts: list[_Part]


@dataclass
class _OneLoopData:
    factors: list[_Factor]
    combinations: list[tuple[tuple[int, ...], int, int]]   # (part index per factor, reduction denominator, degree)
    prefactor: Fraction
    spinor_groups: list[list[str]]                   # spinor labels joined inside each factor
    leg_labels: tuple[tuple[Expression, ...], ...]


class OneLoopBackend:
    """The one-loop back-end (module docstring)."""

    max_loops = 1

    def supports(self, graph: GraphRecord) -> bool:
        return graph.loops == 1

    # ------------------------------------------------------------------ preparation
    def prepare(self, graph: GraphRecord, model: Model) -> PreparedGraph:
        table = model.fields
        if graph.loops != 1:
            raise NotImplementedError("the one-loop back-end needs a one-loop graph")
        if not graph.is_one_particle_irreducible():
            raise ValueError("the graph is not one-particle irreducible")
        omega = graph.superficial_degree(table)
        if omega.denominator != 1:
            raise ValueError(f"non-integer superficial degree {omega}")
        routing = graph.routing()
        loop = LOOP_MOMENTUM
        sigma = []
        for e, mom in enumerate(routing):
            s = mom.get(loop, Fraction(0))
            if s not in (1, -1):
                raise ValueError(f"propagator {e} carries no loop momentum")
            sigma.append(int(s))
        counter = itertools.count(_LABEL_BASE)

        def fresh_like(label: Expression) -> Expression:
            return _relabel(label, next(counter))

        # labels of the propagator ends
        end_labels: dict[tuple[int, int], tuple[Expression, ...]] = {}
        for e, prop in enumerate(graph.propagators):
            for side in (0, 1):
                end_labels[(e, side)] = tuple(fresh_like(lab) for lab in prop.labels[side])
        factors: list[_Factor] = []
        dim = {name: table[name].dimension for name in table}
        slot_counter = itertools.count(_SLOT_BASE)
        # vertices
        for v, vert in enumerate(graph.vertices):
            mapping: dict[str, Expression] = {}
            momenta: dict[int, dict[int, Fraction]] = {}
            labels: list[Expression] = []
            for j, (leg, end) in enumerate(zip(vert.legs, vert.ends)):
                target = graph.legs[end.leg].labels if isinstance(end, ExternalEnd) else end_labels[
                    (end.propagator, end.side)]
                for a, b in zip(leg.labels, target):
                    mapping[_key(a)] = b
                labels.extend(target)
                mom = graph.leg_momentum(v, j, routing)
                if leg.momentum in vert.outgoing:
                    mom = {n: -c for n, c in mom.items()}
                momenta[leg.momentum] = mom
            value = _rename_all(tn.I * vert.value, mapping, fresh_like)
            value = _substitute_momenta(value.replace(tn.d, 4), momenta)
            weight = vert.dimension - sum(dim[f] for f in vert.fields)
            factors.append(self._numerator_factor(f"vertex {v}", value, tuple(labels), frozenset(), weight,
                                                  slot_counter, fresh_like))
        # propagators: numerators and denominators
        for e, prop in enumerate(graph.propagators):
            mapping = {}
            labels = []
            for side in (0, 1):
                for a, b in zip(prop.labels[side], end_labels[(e, side)]):
                    mapping[_key(a)] = b
                labels.extend(end_labels[(e, side)])
            k_e = routing[e]
            value = _rename_all(prop.numerator, mapping, fresh_like)
            value = _substitute_momenta(value.replace(tn.d, 4), {PROPAGATOR_MOMENTUM: k_e})
            power = sum(n for _, n in prop.denominators)
            weight = dim[prop.fields[0]] + dim[prop.fields[1]] - 4 + 2 * power
            lowered = frozenset(_key(lab) for lab in labels if tn.label_space(lab) == "lorentz")
            factors.append(self._numerator_factor(f"propagator {e}", value, tuple(labels), lowered, weight,
                                                  slot_counter, fresh_like))
            shift = {n: sigma[e] * c for n, c in k_e.items() if n != loop}
            for mass2, n in prop.denominators:
                factors.append(self._denominator_factor(f"denominator {e}", shift, mass2.replace(tn.d, 4), n,
                                                        int(omega), slot_counter))
        n_vector = sum(1 for prop in graph.propagators
                       if all(any(s.space.kind == "lorentz" for s in table[f].slots) for f in prop.fields))
        combos = self._combinations(factors, int(omega), n_vector)
        sign = graph.grassmann_sign(table)
        prefactor = Fraction(sign) * graph.symmetry_factor
        degree = max((c[2] for c in combos), default=0)
        groups = []
        for fac in factors:
            spin = [_key(lab) for lab in fac.labels if tn.label_space(lab) == "spinor"]
            if len(spin) > 2:
                raise NotImplementedError(f"{fac.name}: more than one spinor line through a vertex")
            groups.append(spin)
        data = _OneLoopData(factors, combos, prefactor, groups, tuple(leg.labels for leg in graph.legs))
        return PreparedGraph(graph, 1, omega, degree, data)

    def _numerator_factor(self, name: str, value: Expression, labels: tuple[Expression, ...],
                          lowered: frozenset[str], weight: Fraction, slot_counter, fresh_like) -> _Factor:
        loop = LOOP_MOMENTUM
        by_r: dict[int, list[Expression]] = {}
        slot_names: dict[int, tuple[Expression, ...]] = {}
        for term in _terms(value):
            new, slots = _extract_slots(term, loop, fresh_like)
            r = len(slots)
            if r not in slot_names:
                slot_names[r] = tuple(tn.lor(next(slot_counter)) for _ in range(r))
            canon = slot_names[r]
            new = tn.rename(new, list(zip(slots, canon)))
            by_r.setdefault(r, []).append(new)
        parts = []
        for r, terms in sorted(by_r.items()):
            small = weight - r
            if small < 0 or small.denominator != 1:
                raise ValueError(f"{name}: a term of k-degree {r} exceeds the mass dimension {weight}")
            expr = Expression.num(0)
            for t in terms:
                expr = expr + t
            pieces = _split_pieces(expr, labels + slot_names[r])
            if pieces:
                parts.append(_Part(int(small), slot_names[r], 0, pieces))
        return _Factor(name, labels, lowered, parts)

    def _denominator_factor(self, name: str, shift: Mapping[int, Fraction], mass2: Expression, power: int,
                            omega: int, slot_counter) -> _Factor:
        parts = []
        q_zero = not shift
        m_zero = mass2.expand() == Expression.num(0)
        q2 = Expression.num(0)
        for n1, c1 in shift.items():
            for n2, c2 in shift.items():
                q2 = q2 + _num(c1 * c2) * tn.sp(n1, n2)
        for j in range(omega + 1):
            for a in range(j + 1):
                for b in range(j - a + 1):
                    c = j - a - b
                    small = a + 2 * b + 2 * c
                    if small > omega or (q_zero and (a or b)) or (m_zero and c):
                        continue
                    coeff = Fraction(math.comb(power + j - 1, j) * (-1) ** j * math.factorial(j),
                                     math.factorial(a) * math.factorial(b) * math.factorial(c)) * 2 ** a * (-1) ** c
                    slots = tuple(tn.lor(next(slot_counter)) for _ in range(a))
                    expr = _num(coeff)
                    if b:
                        expr = expr * q2 ** b
                    if c:
                        expr = expr * mass2 ** c
                    for s in slots:
                        expr = expr * _combination(shift, tn.p, s)
                    pieces = _split_pieces(expr.expand(), slots)
                    if pieces:
                        parts.append(_Part(small, slots, power + j, pieces))
        return _Factor(name, (), frozenset(), parts)

    @staticmethod
    def _combinations(factors: Sequence[_Factor], omega: int, n_vector: int) -> list[tuple[tuple[int, ...], int, int]]:
        """The combinations of parts of total small degree ``omega``, with the reduction denominator ``4 . 6 ...
        (2 + r)`` and the bound on the degree in ``dh`` (module docstring)."""
        out = []

        def rec(f: int, budget: int, chosen: list[int]) -> None:
            if f == len(factors):
                if budget:
                    return
                parts = [factors[k].parts[i] for k, i in enumerate(chosen)]
                r = sum(len(pt.slots) for pt in parts)
                b = sum(pt.kpow for pt in parts)
                if r != 2 * b - 4:
                    raise AssertionError("power counting: the logarithmic terms do not have r = 2 b - 4")
                denom = math.prod(4 + 2 * i for i in range(r // 2))
                degree = sum(pt.live for pt in parts) // 2 + n_vector + sum(pt.inner for pt in parts)
                out.append((tuple(chosen), denom, degree))
                return
            for i, part in enumerate(factors[f].parts):
                if part.small <= budget:
                    rec(f + 1, budget - part.small, chosen + [i])

        rec(0, omega, [])
        return out

    # ------------------------------------------------------------------ evaluation
    def evaluate(self, prepared: PreparedGraph, fingerprinter: Fingerprinter, samples: Samples) -> Laurent:
        data: _OneLoopData = prepared.data
        if [tuple(_key(x) for x in leg.labels) for leg in samples.legs] != [tuple(_key(x) for x in labs)
                                                                         for labs in data.leg_labels]:
            raise ValueError("the samples are not taken on the legs of the graph")
        fp = fingerprinter
        F = fp.field
        p = F.p
        d_ext = fp.dimension
        pref = F(data.prefactor)
        degrees = sorted({c[2] for c in data.combinations})
        values: dict[int, list[np.ndarray]] = {b: [] for b in degrees}
        for m in range(prepared.degree + 1):
            fp_int = fp if m == 0 else _with_dimension(fp, d_ext + m)
            smp = embed_samples(samples, fp, fp_int)
            todo = [b for b in degrees if b >= m]
            for b, val in self._value(data, fp_int, smp, todo).items():
                values[b].append(val)
        total: dict[int, np.ndarray] = {}
        for b, vals in values.items():
            points = [d_ext - 4 + m for m in range(b + 1)]
            vander = np.array([[pow(x, j, p) for j in range(b + 1)] for x in points], dtype=np.int64)
            coeffs = la.solve(vander, np.stack(vals), F, require_unique=True)
            for r in range(b + 1):
                total[r] = (total.get(r, 0) + np.asarray(coeffs[r], dtype=np.int64)) % p
        out = {}
        for r, c in sorted(total.items()):
            c = c * pref % p * pow(-2, r, p) % p
            if np.any(c):
                out[(r, 1)] = c
        return Laurent(1, F, samples.n, out)

    def _value(self, data: _OneLoopData, fp: Fingerprinter, samples: Samples, degrees: Sequence[int]
               ) -> dict[int, np.ndarray]:
        """The pole coefficient (before the prefactor) in the dimension of ``fp``, summed per degree group."""
        p = fp.field.p
        dense = _Dense(fp, samples)
        cache: dict[tuple[int, int, int], tuple[np.ndarray, list[str], np.ndarray, list[str]]] = {}

        def piece_ops(f: int, i: int, k: int):
            key = (f, i, k)
            if key not in cache:
                piece = data.factors[f].parts[i].pieces[k]
                ld, b1 = dense.tensor(piece.ld, piece.ld_labels)
                ld_keys = [_key(lab) for lab in piece.ld_labels]
                lowered = data.factors[f].lowered
                for axis, lab in enumerate(ld_keys, start=int(b1)):
                    if lab in lowered:
                        shape = [1] * ld.ndim
                        shape[axis] = -1
                        ld = ld * fp.rep.metric.reshape(shape) % p
                gr, b2 = dense.tensor(piece.group, piece.group_labels)
                if b2:
                    raise AssertionError("group tensors do not depend on the samples")
                if not piece.group_labels:
                    ld = ld * int(gr) % p
                cache[key] = (ld, ([_BATCH] if b1 else []) + ld_keys, gr, [_key(lab) for lab in piece.group_labels])
            return cache[key]

        pols = []
        for leg, pol in zip(samples.legs, samples.polarisations):
            keys = [_key(lab) for lab in leg.labels]
            out = pol
            for axis, slot in enumerate(fp.model.fields[leg.field].slots, start=1):
                if slot.space.kind == "lorentz":
                    shape = [1] * out.ndim
                    shape[axis] = -1
                    out = out * fp.rep.metric.reshape(shape) % p
            pols.append((out, [_BATCH] + keys))
        # closed spinor cycles of the graph
        parent: dict[str, str] = {}

        def find(x: str) -> str:
            while parent.setdefault(x, x) != x:
                x = parent[x]
            return x

        for spin in data.spinor_groups:
            for a, b in itertools.pairwise(spin):
                parent[find(a)] = find(b)
        roots = {find(k) for spin in data.spinor_groups for k in spin}
        open_roots = {find(_key(lab)) for labs in data.leg_labels for lab in labs
                      if tn.label_space(lab) == "spinor" and _key(lab) in parent}
        norm = pow(fp.rep.trace_norm, len(roots - open_roots), p)
        totals = {b: np.zeros(samples.n, dtype=np.int64) for b in degrees}
        for chosen, denom, degree in data.combinations:
            if degree not in totals:
                continue
            total = totals[degree]
            parts = [data.factors[f].parts[i] for f, i in enumerate(chosen)]
            slots = [_key(s) for pt in parts for s in pt.slots]
            weight = fp.field(Fraction(1, denom)) * norm % p
            for choice in itertools.product(*[range(len(pt.pieces)) for pt in parts]):
                base_ops = []
                for f, (i, k) in enumerate(zip(chosen, choice)):
                    ld, ld_keys, gr, gr_keys = piece_ops(f, i, k)
                    base_ops.append((ld, ld_keys))
                    if gr_keys:
                        base_ops.append((gr, gr_keys))
                for pairing in _pairings(slots):
                    ren = {b: a for a, b in pairing}
                    ops = [(a, [ren.get(x, x) for x in labs]) for a, labs in base_ops]
                    ops += [(fp.rep.metric, [a]) for a, _ in pairing]
                    ops += pols
                    val = _contract(ops, [_BATCH], p, dense.paths)
                    total = (total + val * weight) % p
            totals[degree] = total
        return totals


def _pairings(slots: Sequence[str]) -> Iterable[list[tuple[str, str]]]:
    if not slots:
        yield []
        return
    a, rest = slots[0], slots[1:]
    for i in range(len(rest)):
        for tail in _pairings(rest[:i] + rest[i + 1:]):
            yield [(a, rest[i])] + tail


def _rename_all(expr: Expression, mapping: Mapping[str, Expression], fresh_like) -> Expression:
    """Rename the labels in ``mapping`` and give every other label of ``expr`` a fresh graph-unique name."""
    pairs = dict(mapping)
    for term in _terms(expr):
        for lab in tn.labels_in(term):
            k = _key(lab)
            if k not in pairs:
                pairs[k] = fresh_like(lab)
    return tn.rename(expr, pairs.items()).expand()


def _extract_slots(term: Expression, loop: int, fresh_like) -> tuple[Expression, list[Expression]]:
    """Remove the loop momentum from a product term: every occurrence becomes an open Lorentz slot (``k^mu ->
    g^{mu kappa}``, ``k.q -> q^kappa``, ``kbar.q -> gbar^{kappa lambda} q_lambda``, ``k^2 -> g^{kappa1 kappa2}``, ...).
    Returns the new term and the slot labels in order."""
    out = Expression.num(1)
    slots: list[Expression] = []

    def slot() -> Expression:
        s = tn.lor(10 ** 9 + len(slots))
        slots.append(s)
        return s

    metric = {tn.sp.get_name(): tn.g, tn.spbar.get_name(): tn.gbar, tn.sphat.get_name(): tn.ghat}
    for fa in _factors(term):
        if fa.get_type() != AtomType.Fn:
            out = out * fa
            continue
        name = fa.get_name()
        args = list(fa)
        if name == tn.p.get_name() and _is_number(args[0], loop):
            out = out * tn.g(args[1], slot())
        elif name in metric:
            ks = [a for a in args if _is_number(a, loop)]
            if not ks:
                out = out * fa
            elif len(ks) == 2:
                a, b = slot(), slot()
                out = out * metric[name](a, b)
            else:
                other = next(a for a in args if not _is_number(a, loop))
                s = slot()
                if name == tn.sp.get_name():
                    out = out * tn.p(other, s)
                else:
                    lam = fresh_like(s)
                    out = out * metric[name](s, lam) * tn.p(other, lam)
        else:
            out = out * fa
    return out, slots


def _is_number(e: Expression, n: int) -> bool:
    return e.get_type() == AtomType.Num and e.to_canonical_string() == str(n)


def _split_pieces(expr: Expression, open_labels: Sequence[Expression]) -> list[_Piece]:
    """Split an expression into pieces ``group * (Lorentz-Dirac sum)`` with distinct group structures."""
    group_of: dict[str, tuple[Expression, list[Expression]]] = {}
    for term in _terms(expr):
        ld, gr = Expression.num(1), Expression.num(1)
        for fa in _factors(term):
            if _is_ld(fa):
                ld = ld * fa
            else:
                gr = gr * fa
        gr_open = [lab for lab in open_labels if tn.label_space(lab) not in ("lorentz", "spinor")]
        gr = tn.canonical_relabel(gr, gr_open) if gr.get_type() != AtomType.Num else gr
        key = gr.to_canonical_string()
        group_of.setdefault(key, (gr, []))[1].append(ld)
    pieces = []
    for gr, lds in group_of.values():
        total = Expression.num(0)
        for x in lds:
            total = total + x
        total = total.expand()
        if total == Expression.num(0):
            continue
        ld_labels = tuple(lab for lab in open_labels if tn.label_space(lab) in ("lorentz", "spinor"))
        gr_labels = tuple(lab for lab in open_labels if tn.label_space(lab) not in ("lorentz", "spinor"))
        missing = [lab for lab in gr_labels if not any(_key(x) == _key(lab) for x in tn.labels_in(gr))]
        if missing:
            raise ValueError(f"group labels {missing} do not occur in the group part {gr}")
        pieces.append(_Piece(gr, gr_labels, total, ld_labels))
    return pieces


_DEAD = (tn.p, tn.gbar, tn.gammabar, tn.eps)
_METRICS = (tn.g, tn.ghat)


def _lorentz_lines(term: Expression):
    """Index lines of a product term: ``(find, dead)`` with ``find`` the component of a Lorentz label (labels joined
    through the metrics ``g``, ``ghat``) and ``dead`` the components that touch a head without hatted components
    (``p``, ``gbar``, ``gammabar``, ``eps``)."""
    dead_names = {h.get_name() for h in _DEAD}
    metric_names = {h.get_name() for h in _METRICS}
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        while parent.setdefault(x, x) != x:
            x = parent[x]
        return x

    killed: set[str] = set()
    for fa in _factors(term):
        if fa.get_type() != AtomType.Fn:
            continue
        labs = [_key(lab) for lab in tn.labels_in(fa) if tn.label_space(lab) == "lorentz"]
        for k in labs:
            find(k)
        if fa.get_name() in metric_names and len(labs) == 2:
            parent[find(labs[0])] = find(labs[1])
        if fa.get_name() in dead_names:
            killed.update(labs)
    return find, {find(k) for k in killed}


def _inner_lorentz(expr: Expression, open_labels: Sequence[Expression]) -> int:
    """Maximum number, over the terms, of closed index lines inside one term (:func:`_lorentz_lines` components
    that touch no open label and no dead head)."""
    opened = [_key(lab) for lab in open_labels if tn.label_space(lab) == "lorentz"]
    best = 0
    for term in _terms(expr):
        find, dead = _lorentz_lines(term)
        labels = {_key(lab) for lab in tn.labels_in(term) if tn.label_space(lab) == "lorentz"}
        bad = dead | {find(k) for k in opened}
        best = max(best, len({find(k) for k in labels} - bad))
    return best


# ---------------------------------------------------------------------------------------------- many graphs
_POOL_STATE: dict[str, object] = {}


def _pool_task(index: int) -> tuple[int, dict[tuple[int, int], np.ndarray]]:
    backend, jobs = _POOL_STATE["backend"], _POOL_STATE["jobs"]
    pg, fp, smp = jobs[index]
    return index, dict(backend.evaluate(pg, fp, smp).coefficients)


def evaluate_jobs(jobs: Sequence[tuple[PreparedGraph, Fingerprinter, Samples]], *,
                  backend: LoopBackend | None = None, processes: int = 1) -> list[Laurent]:
    """:meth:`LoopBackend.evaluate` for many ``(prepared graph, fingerprinter, samples)`` jobs (any mix of classes,
    primes, dimensions and parameter points). With ``processes > 1`` the jobs are distributed over forked worker
    processes, the most expensive first (estimated by the number of part combinations and the degree in the hatted
    dimension); the results are identical, since everything is exact and deterministic."""
    be = backend or OneLoopBackend()
    jobs = list(jobs)
    if processes <= 1 or len(jobs) < 2:
        return [be.evaluate(pg, fp, smp) for pg, fp, smp in jobs]
    import multiprocessing

    def cost(k: int) -> int:
        pg, fp, _ = jobs[k]
        combos = len(getattr(pg.data, "combinations", ()) or ())
        return (combos + 1) * (pg.degree + 1) * fp.dimension

    order = sorted(range(len(jobs)), key=lambda k: -cost(k))
    _POOL_STATE.update(backend=be, jobs=jobs)
    try:
        with multiprocessing.get_context("fork").Pool(min(processes, len(jobs))) as pool:
            results = dict(pool.imap_unordered(_pool_task, order))
    finally:
        _POOL_STATE.clear()
    return [Laurent(pg.loops, fp.field, smp.n, results[k]) for k, (pg, fp, smp) in enumerate(jobs)]


def evaluate_graphs(prepared: Sequence[PreparedGraph], fingerprinter: Fingerprinter,
                    samples: Mapping[tuple[str, ...], Samples], *, backend: LoopBackend | None = None,
                    processes: int = 1) -> list[Laurent]:
    """:meth:`LoopBackend.evaluate` for many prepared graphs with one fingerprinter; ``samples`` maps the leg class
    (field names in leg order) to its samples (:func:`evaluate_jobs`)."""
    return evaluate_jobs([(pg, fingerprinter, samples[tuple(leg.field for leg in pg.graph.legs)]) for pg in prepared],
                         backend=backend, processes=processes)
