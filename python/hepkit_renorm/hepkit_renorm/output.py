"""The model writer: an :class:`~hepkit_renorm.store.OrderStore` as a HepKIT / FeynKit JSON model.

The package, not the JSON strings, is the source of truth for values and signs: the JSON carries the combinatorics
(particles, vertex legs, propagators, coupling orders) and values that agree with the package's rules, and the
reverse map (:class:`ReverseMap`, also embedded in the file) sends every particle and vertex rule back to the
package's field, components and vertex function. Diagram signs are never taken from the generator (see below).

Encoding in today's format (schema ``"v1"``, accepted by the installed HepKIT)
------------------------------------------------------------------------------
* **Particles.** One HepKIT particle per component of every written field. Slots of the *native colour* factor
  (the unique ``SU(3)`` factor by default, HepKIT's ``color`` 3 / -3 / 8) stay native; every other group slot
  (``SU(2)``, ...) and the flavour slot is written in components, numbered ``1..n`` in the bases of
  :mod:`hepkit_renorm.groups` (generalised Gell-Mann / 2; the real adjoint basis). Particle and antiparticle are
  the pairs of the quadratic form (``Psi``/``Psibar``, ``omega``/``omegabar``, ``H``/``H^dagger``) and, for
  sources, the sources of such pairs (``R_f``/``Rbar_f``, ``sigma``/``sigma^dagger``); real fields are
  self-conjugate. The particle side is the field of kind ``fermion``, ``ghost``, ``scalar`` or ``gauge``, and for a
  source the source of the conjugate-side field (``R_f`` is the source of ``Psibar_f``). Name: the particle-side
  field name plus the 1-based components, ``~`` appended for the antiparticle. UFO spin ``2s + 1`` from the slots
  (Lorentz 3, spinor 2, else 1) and ``-1`` for ghosts and antighosts. Fields without vertices (``b``, the source
  of an abelian ghost) are not written.
* **Grassmann workaround.** ``ghost_number`` is 0 for every particle (today's generator rejects vertices with
  three odd legs and external ghost/source chains; README of the hand-written model, item 1). Closed ghost loops
  therefore lack their ``-1`` and the ordering signs of ghosts and sources are not computed by HepKIT: signs are
  the package's (one global Wick sign from the leg orders of the vertex rules, scout notes on generation). BRST
  sources have ``propagating: false`` and no propagator; today's generator ignores the flag, so generation must
  exclude internal sources (coupling order ``SRC`` = number of external sources, a filter, a vertex bound).
* **Vertex rules.** Leg ``k`` of a rule is the ``k``-th field of the package's leg class (canonical field order),
  so the rule is ``i V`` of :func:`hepkit_renorm.rules.vertex_rule` for exactly that leg order, with all momenta
  kept (``P(mu, k)`` = incoming momentum of leg ``k``) and the non-native group / flavour labels set to the
  components of the particles. Legs of the same field carry non-decreasing components (one rule per unordered
  set of particles). The value is ``sum_{c, l} coupling[c][l] * colour_c * lorentz_l``; every coupling is ``i``
  times an exact complex rational times parameters (and, for counterterms, powers of the regulator ``eps``).
* **Lorentz structures.** ``Gamma``, ``Gamma5``, ``ProjM`` (= P_L), ``ProjP`` (= P_R), ``Metric``, ``P``,
  ``Identity`` (spinor), ``Epsilon`` (= minus the package's ``eps``, ``eps^{0123} = +1``, assuming UFO's
  ``epsilon_{0123} = +1``; not used at tree level). A four-dimensional (barred) index ``X`` is marked by the factor
  ``bmhv::bar(X)``, an evanescent (hatted) one by ``bmhv::hat(X)``, once per index; both survive today's lowering
  with renamed indices. Leg indices are ``idx(1,k)``, dummies ``dummy(n)``.
* **Colour.** ``T(a, i, jb)`` with ``i`` the quark leg and ``jb`` the antiquark leg (= the package's
  ``T(a, jb, i)``, row = barred leg), ``f``, ``d`` (= ``dsym``), ``Identity``; dummies are negative integers.
* **Propagators** (one per particle and antiparticle): ``D_{Y Ybar}(k) = i K^{-1}`` of
  :func:`hepkit_renorm.rules.propagators` for the particle ``Y`` with the labels of the ``Y`` end as
  ``idx(1,2)``, those of the ``Ybar`` end as ``idx(1,1)`` and ``k`` = ``P`` the edge momentum (along the particle
  flow, into the ``Y`` end). The antiparticle uses the same object with the ends exchanged and ``k -> -k``.
  ``gammabar + gammahat`` is written as the full ``gamma`` (``PSlash``). Internal Kronecker deltas are implicit.
* **Coupling orders.** One order per interaction type of the specification: one per gauge factor
  (``G<factor>``), one for the Yukawa couplings (``YUK``), one for the scalar self-couplings (``POT``, potential
  terms with at least two bilinears; mass terms carry no order); the power of each parameter adds to its order.
  Roles: ``DELTA`` = 1 on breaking vertices, ``SRC`` = number of source legs, ``CT`` = loop order ``L`` of a
  counterterm or breaking vertex of order ``L >= 1``, ``CTPOLE`` = 1 on singular and ``CTFIN`` = 1 on finite
  counterterm vertices. A diagram of order ``l`` has ``loops + CT = l``. Names: :class:`OrderNames`.

Counterterms
------------
A store entry ``sum_k eps^k C_k`` is written with couplings ``sum_k eps^k c_k`` in the expression (``d`` re-expanded
as ``4 - 2 eps``), the regulator printed as ``eps`` and the evanescent tag as ``eta``. Every rule of a
counterterm entry carries its own orders (no mixing of orders inside a rule: terms with different parameter
orders become separate rules with the same particles).

The proposed schema v2 (scout notes on the model format, section 5)
-------------------------------------------------------------------
Every v2 field is written as an additional key, which today's loader ignores: particle ``statistics``,
``field_kind``, ``source_of``, ``representations``, ``multiplet``; vertex ``loop_order``, ``role`` (``tree``,
``source``, ``breaking``, ``ct_singular``, ``ct_finite``), ``operator`` (the store entry), ``grassmann_order``,
``evanescent``; coupling ``laurent`` (``{"variable": "eps", "terms": {power: coefficient}}``); order ``kind``;
top-level ``regularization``, ``gauge_groups``, ``operators``. Values that v2 changes are kept under
``v2_overrides`` (the true ``ghost_number``; the Lorentz ``structure`` with ``Bar(i)`` / ``Hat(i)`` wrappers
instead of markers). ``schema="v2"`` applies the overrides and writes ``schema_version: 2``: the switch to native
v2 fields is this option (it needs the extended loader; today's loader mis-lowers the wrappers).
"""

from __future__ import annotations

import hashlib
import itertools
import json
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from fractions import Fraction
from pathlib import Path

import numpy as np
from symbolica import AtomType, Expression, S

from . import tensors as tn
from .exact import default_field
from .fields import FieldInfo, FieldTable
from .fingerprint import Fingerprinter, Samples
from .groups import complex_rational
from .linalg import reconstruct_rationals, split_gaussian
from .rules import VertexRule, is_vertex_class, propagators, vertex_rule
from .spec import Model
from .store import EVANESCENT_TAG, OrderStore, Role, StoreEntry, tree_store

__all__ = [
    "OrderNames",
    "ParticleRecord",
    "PropagatorRecord",
    "ReverseMap",
    "VertexRecord",
    "WriterOptions",
    "WrittenModel",
    "component_samples",
    "tree_model",
    "write_model",
]

_I = "1\U0001d456"
_E = S("hepkit_renorm_output::eps")          # the regulator while couplings are assembled
_PSLASH = S("hepkit_renorm_output::PSlash")
_w = {n: S(n + "_") for n in ("a", "s", "t")}
_RESERVED = {"P", "PSlash", "Metric", "Gamma", "Gamma5", "ProjM", "ProjP", "Identity", "Epsilon", "T", "f", "d",
             "idx", "dummy", "pi", "ZERO", "Bar", "Hat"}
_COLOUR_CODE = {"F": 3, "Fbar": -3, "A": 8}
_FIELD_KIND = {"gauge": "field", "fermion": "field", "fermion_bar": "field", "scalar": "field",
               "scalar_conj": "field", "ghost": "ghost", "antighost": "antighost", "auxiliary": "auxiliary",
               "source": "brst_source"}
_PARTICLE_KINDS = ("gauge", "fermion", "ghost", "scalar")
_CONJUGATE_KINDS = (("fermion", "fermion_bar"), ("ghost", "antighost"), ("scalar", "scalar_conj"))


# ---------------------------------------------------------------------------------------------- options
@dataclass(frozen=True)
class OrderNames:
    """Names of the coupling orders (module docstring). ``gauge`` is a template with ``{factor}``."""

    gauge: str = "G{factor}"
    yukawa: str = "YUK"
    potential: str = "POT"
    breaking: str = "DELTA"
    source: str = "SRC"
    loop: str = "CT"
    singular: str = "CTPOLE"
    finite: str = "CTFIN"


@dataclass(frozen=True)
class WriterOptions:
    """Options of :func:`write_model`.

    ``schema``: ``"v1"`` (today's format, v2 content as extra keys) or ``"v2"`` (native v2 fields).
    ``native_colour``: the gauge factor written as HepKIT colour (``"auto"``: the unique SU(3) factor, if any;
    ``None``: all factors in components). ``substitutions``: values of parameters by their specification name
    (``{"xi": 1}`` for Feynman gauge). ``max_loop_order``: write only store entries up to this order.
    ``symbols``: further symbols allowed in couplings, ``{Symbolica name: UFO name}``. ``dimensions``,
    ``samples``, ``seed``: the exact zero test of every component vertex (fingerprints in each ``D`` of
    ``dimensions``, ``d = D``; a component is written if it is non-zero in one of them).
    """

    schema: str = "v1"
    name: str | None = None
    native_colour: str | None = "auto"
    substitutions: Mapping[str, int | Fraction] = field(default_factory=dict)
    max_loop_order: int | None = None
    orders: OrderNames = field(default_factory=OrderNames)
    regulator: str = "eps"
    tag: str = "eta"
    symbols: Mapping[str, str] = field(default_factory=dict)
    dimensions: tuple[int, ...] = (6, 8, 4)
    samples: int = 3
    seed: int = 0

    def __post_init__(self) -> None:
        if self.schema not in ("v1", "v2"):
            raise ValueError("schema is 'v1' or 'v2'")


# ---------------------------------------------------------------------------------------------- reverse map
@dataclass(frozen=True)
class ParticleRecord:
    """A HepKIT particle: package ``field`` with ``components`` = ``{slot position: 0-based component}`` of its
    component slots; ``antiname`` is the HepKIT conjugate."""

    name: str
    field: str
    components: tuple[tuple[int, int], ...]
    antiname: str
    pdg_code: int

    def to_dict(self) -> dict:
        return {"field": self.field, "components": [list(c) for c in self.components], "antiname": self.antiname,
                "pdg_code": self.pdg_code}


@dataclass(frozen=True)
class VertexRecord:
    """A vertex rule: legs = ``fields`` (the package leg class, in this order) with ``components`` per leg; the
    rule is ``i V`` of the store entry ``entry`` (``sum_k eps^k`` over ``powers``) for that leg order."""

    name: str
    entry: str
    loop_order: int
    role: str
    fields: tuple[str, ...]
    components: tuple[tuple[tuple[int, int], ...], ...]
    particles: tuple[str, ...]
    orders: tuple[tuple[str, int], ...]
    powers: tuple[int, ...]
    evanescent: bool

    def to_dict(self) -> dict:
        return {"entry": self.entry, "loop_order": self.loop_order, "role": self.role, "fields": list(self.fields),
                "components": [[list(c) for c in leg] for leg in self.components],
                "particles": list(self.particles), "orders": [list(o) for o in self.orders],
                "powers": list(self.powers), "evanescent": self.evanescent}


@dataclass(frozen=True)
class PropagatorRecord:
    """A HepKIT propagator of ``particle``: the package ``D_{fields[0] fields[1]}`` with the ``fields[0]`` end at
    ``idx(1,2)``; ``reversed``: written as the antiparticle (``k -> -k``, ends exchanged)."""

    name: str
    particle: str
    fields: tuple[str, str]
    reversed: bool

    def to_dict(self) -> dict:
        return {"particle": self.particle, "fields": list(self.fields), "reversed": self.reversed}


@dataclass
class ReverseMap:
    """HepKIT names -> package objects. ``parameters``: UFO name -> ``(kind, reference)`` with kind
    ``coupling`` (reference: gauge factor), ``parameter`` (specification name), ``gauge_parameter``,
    ``regulator`` or ``symbol`` (Symbolica name); ``orders``: name -> description."""

    particles: dict[str, ParticleRecord] = field(default_factory=dict)
    vertices: dict[str, VertexRecord] = field(default_factory=dict)
    propagators: dict[str, PropagatorRecord] = field(default_factory=dict)
    parameters: dict[str, tuple[str, str]] = field(default_factory=dict)
    orders: dict[str, dict] = field(default_factory=dict)

    def particle(self, name: str) -> tuple[str, dict[int, int]]:
        """``(field, {slot: component})`` of a HepKIT particle."""
        rec = self.particles[name]
        return rec.field, dict(rec.components)

    def vertex(self, name: str) -> VertexRecord:
        return self.vertices[name]

    def symbol(self, model: Model, ufo_name: str) -> Expression:
        """The package symbol of a UFO parameter name."""
        kind, ref = self.parameters[ufo_name]
        if kind == "coupling":
            return model.coupling(ref)
        if kind == "parameter":
            return model.param(ref)
        if kind == "gauge_parameter":
            return model.gauge_fixing.xi
        if kind == "symbol":
            return S(ref)
        raise KeyError(f"{ufo_name} ({kind}) has no package symbol")

    def to_dict(self) -> dict:
        return {"particles": {k: v.to_dict() for k, v in self.particles.items()},
                "vertices": {k: v.to_dict() for k, v in self.vertices.items()},
                "propagators": {k: v.to_dict() for k, v in self.propagators.items()},
                "parameters": {k: list(v) for k, v in self.parameters.items()},
                "orders": dict(self.orders)}

    @classmethod
    def from_dict(cls, doc: Mapping) -> ReverseMap:
        def comps(c):
            return tuple(tuple(int(x) for x in pair) for pair in c)
        out = cls()
        for k, v in doc["particles"].items():
            out.particles[k] = ParticleRecord(k, v["field"], comps(v["components"]), v["antiname"], v["pdg_code"])
        for k, v in doc["vertices"].items():
            out.vertices[k] = VertexRecord(k, v["entry"], v["loop_order"], v["role"], tuple(v["fields"]),
                                           tuple(comps(c) for c in v["components"]), tuple(v["particles"]),
                                           tuple((o, int(n)) for o, n in v["orders"]), tuple(v["powers"]),
                                           v["evanescent"])
        for k, v in doc["propagators"].items():
            out.propagators[k] = PropagatorRecord(k, v["particle"], tuple(v["fields"]), v["reversed"])
        out.parameters = {k: tuple(v) for k, v in doc["parameters"].items()}
        out.orders = dict(doc["orders"])
        return out


@dataclass
class WrittenModel:
    """The written model: ``data`` (the JSON document), its ``reverse_map``, and the sources (``model``,
    ``store``, ``options``)."""

    data: dict
    reverse_map: ReverseMap
    model: Model
    store: OrderStore
    options: WriterOptions

    def json_string(self, indent: int | None = 1) -> str:
        return json.dumps(self.data, indent=indent, ensure_ascii=False)

    def write(self, path: str | Path) -> Path:
        path = Path(path)
        path.write_text(self.json_string())
        return path

    def load(self):
        """The model loaded by the installed HepKIT (``symbolica.community.hepkit.Model``)."""
        from symbolica.community import hepkit

        return hepkit.Model.from_json(self.json_string(None))

    def package_rule(self, vertex: str, power: int = 0) -> VertexRule:
        """The package's vertex rule (``V``, all momenta kept, legs in the order of the HepKIT rule) behind a written
        vertex, for the coefficient of ``eps^power`` of its store entry; the components of the legs are in
        :attr:`VertexRecord.components` (evaluate with :func:`component_samples`)."""
        rec = self.reverse_map.vertices[vertex]
        entry = self.store.entry(rec.entry)
        rule = vertex_rule(entry.coefficient(power), rec.fields, self.model.simplifier(expand_projectors=False),
                           eliminate=None)
        subs = _substitution_rules(self.model, self.options)
        return replace(rule, value=_substitute(rule.value, subs))


# ---------------------------------------------------------------------------------------------- samples
def component_samples(fingerprinter: Fingerprinter, fields: Sequence[str],
                      components: Sequence[Iterable[tuple[int, int]]], n: int, salt: int = 0) -> Samples:
    """Fingerprint samples of a leg class with the component slots fixed: leg ``j`` has its polarisation
    restricted to the component ``c`` on every slot ``s`` of ``components[j]`` (pairs ``(s, c)``); evaluating a
    vertex function on them evaluates its component ``c`` (times the random polarisation entry)."""
    smp = fingerprinter.samples(list(fields), n, salt)
    pols = []
    for pol, comps in zip(smp.polarisations, components):
        mask = np.ones(pol.shape[1:], dtype=np.int64)
        for slot, c in comps:
            sel = np.zeros(pol.shape[1 + slot], dtype=np.int64)
            sel[c] = 1
            shape = [1] * (pol.ndim - 1)
            shape[slot] = -1
            mask = mask * sel.reshape(shape)
        pols.append(pol * mask[None, ...])
    return Samples(smp.legs, smp.momenta, tuple(pols))


# ---------------------------------------------------------------------------------------------- symbols
@dataclass(frozen=True)
class _Symbol:
    ufo: str
    kind: str
    ref: str
    order: str | None


def _symbol_table(model: Model, options: WriterOptions) -> dict[str, _Symbol]:
    """Symbolica name -> UFO name, kind and order of every symbol allowed in a coupling."""
    on = options.orders
    yukawa = {y.coupling for y in model.yukawas}
    potential = {t.coupling for t in model.potential if len(t.bilinears) >= 2}
    table: dict[str, _Symbol] = {}
    for name, fa in model.factors.items():
        table[model.coupling(name).get_name()] = _Symbol(fa.coupling, "coupling", name, on.gauge.format(factor=name))
    for name in model.parameters:
        order = on.yukawa if name in yukawa else on.potential if name in potential else None
        table[model.param(name).get_name()] = _Symbol(name, "parameter", name, order)
    xi = model.gauge_fixing.xi
    if xi.get_type() == AtomType.Var:
        table[xi.get_name()] = _Symbol(xi.get_name().split("::")[-1], "gauge_parameter", "xi", None)
    table[_E.get_name()] = _Symbol(options.regulator, "regulator", "eps", None)
    table[EVANESCENT_TAG.get_name()] = _Symbol(options.tag, "regulator", "eta", None)
    table[Expression.PI.get_name()] = _Symbol("pi", "constant", "pi", None)
    for sym, ufo in options.symbols.items():
        table[sym] = _Symbol(ufo, "symbol", sym, None)
    for name in options.substitutions:
        if not any(s.ufo == name and s.kind in ("coupling", "parameter", "gauge_parameter") for s in table.values()):
            raise KeyError(f"cannot substitute the unknown parameter {name!r}")
    names = [s.ufo for s in table.values()]
    clash = {n for n in names if names.count(n) > 1} | (set(names) - {"pi"}) & _RESERVED
    if clash:
        raise ValueError(f"UFO names {sorted(clash)} are not unique or collide with UFO heads")
    return table


def _substitution_rules(model: Model, options: WriterOptions) -> list[tuple[Expression, Expression]]:
    out = []
    for sym in _symbol_table(model, options).values():
        if sym.ufo in options.substitutions and sym.kind in ("coupling", "parameter", "gauge_parameter"):
            q = Fraction(options.substitutions[sym.ufo])
            target = (model.gauge_fixing.xi if sym.kind == "gauge_parameter" else
                      model.coupling(sym.ref) if sym.kind == "coupling" else model.param(sym.ref))
            out.append((target, Expression.num(q.numerator) / q.denominator))
    return out


def _substitute(expr: Expression, rules: Sequence[tuple[Expression, Expression]]) -> Expression:
    for a, b in rules:
        expr = expr.replace(a, b)
    return expr.expand()


def _number(num: Expression) -> str:
    re, im = complex_rational(num)
    if im == 0:
        return f"({re})"
    if re == 0:
        return f"({im})*{_I}"
    return f"(({re})+({im})*{_I})"


def _scalar_string(expr: Expression, table: Mapping[str, _Symbol]) -> str:
    """UFO string of a coupling expression (numbers, symbols of ``table``, sums, products, powers)."""
    t = expr.get_type()
    if t == AtomType.Num:
        return _number(expr)
    if t == AtomType.Var:
        name = expr.get_name()
        if name not in table:
            raise ValueError(f"symbol {name} has no UFO name (pass it in WriterOptions.symbols)")
        return table[name].ufo
    if t == AtomType.Add:
        return "(" + "+".join(_scalar_string(a, table) for a in expr) + ")"
    if t == AtomType.Mul:
        return "*".join(_scalar_string(a, table) for a in expr)
    if t == AtomType.Pow:
        base, ex = list(expr)
        if ex.get_type() != AtomType.Num:
            raise NotImplementedError(f"symbolic exponent in {expr}")
        return f"({_scalar_string(base, table)})^({Fraction(ex.to_canonical_string())})"
    raise NotImplementedError(f"cannot write {expr} in a coupling (complex or flavour-matrix parameters)")


def _orders_of(scalar: Expression, table: Mapping[str, _Symbol]) -> dict[str, int]:
    out: dict[str, int] = {}
    for fa in (list(scalar) if scalar.get_type() == AtomType.Mul else [scalar]):
        base, n = fa, 1
        if fa.get_type() == AtomType.Pow:
            base, ex = list(fa)
            n = int(ex.to_canonical_string())
        if base.get_type() == AtomType.Var:
            sym = table.get(base.get_name())
            if sym is not None and sym.order is not None:
                out[sym.order] = out.get(sym.order, 0) + n
    return out


def _laurent(expr: Expression) -> dict[int, Expression]:
    """``{k: c_k}`` of a Laurent polynomial ``sum_k c_k eps^k`` in the assembly regulator."""
    out: dict[int, Expression] = {}
    expr = expr.expand()
    for term in (list(expr) if expr.get_type() == AtomType.Add else [] if expr == Expression.num(0) else [expr]):
        k = 0
        for fa in (list(term) if term.get_type() == AtomType.Mul else [term]):
            if fa == _E:
                k = 1
            elif fa.get_type() == AtomType.Pow and next(iter(fa)) == _E:
                k = int(list(fa)[1].to_canonical_string())
        c = (term / _E**k).expand()
        out[k] = out[k] + c if k in out else c
    return {k: c for k, c in sorted(out.items()) if c.expand() != Expression.num(0)}


# ---------------------------------------------------------------------------------------------- structures
@dataclass
class _Printer:
    """UFO strings of Lorentz / colour / propagator structures. ``index``: label key -> leg index string;
    dummies are numbered by first appearance; barred / hatted indices are collected for the markers."""

    index: Mapping[str, str]
    colour: bool = False
    propagator: bool = False
    dummies: dict[str, str] = field(default_factory=dict)
    bar: dict[str, None] = field(default_factory=dict)
    hat: dict[str, None] = field(default_factory=dict)
    sign: int = 1

    def _dummy(self, key: str | None = None) -> str:
        n = len(self.dummies) + 1
        tok = f"-{n}" if self.colour else f"dummy({n})"
        self.dummies[key if key is not None else f"#{n}"] = tok
        return tok

    def label(self, lab: Expression) -> str:
        key = lab.to_canonical_string()
        if key in self.index:
            return self.index[key]
        if key in self.dummies:
            return self.dummies[key]
        return self._dummy(key)

    def pieces(self, fa: Expression) -> list[tuple[str, list[str], list[str]]]:
        """``[(head, index args, trailing plain args)]`` of one factor."""
        name = fa.get_name()
        args = list(fa)
        H = tn.HEADS
        lab = self.label
        if self.colour:
            if name == H["T"].get_name():
                a, r, c = (lab(x) for x in args)
                return [("T", [a, c, r], [])]
            if name == H["f"].get_name():
                return [("f", [lab(x) for x in args], [])]
            if name == H["dsym"].get_name():
                return [("d", [lab(x) for x in args], [])]
            if name == H["delta"].get_name():
                return [("Identity", [lab(x) for x in args], [])]
            raise NotImplementedError(f"colour head {name} is not written (epsilon tensors of colour)")
        if name == H["p"].get_name():
            k = int(args[0].to_canonical_string())
            mu = lab(args[1])
            return [("P", [mu], [] if self.propagator else [str(k)])]
        if name in (H["sp"].get_name(), H["spbar"].get_name(), H["sphat"].get_name()):
            ks = [int(a.to_canonical_string()) for a in args]
            if ks[0] == ks[1] and not self.propagator:
                # p_k.p_k: two dummies and a metric (P(d,k)*P(d,k) would be read back as the power P(d,k)^2)
                d1, d2 = self._dummy(), self._dummy()
                self._mark(name, [d1, d2])
                return [("Metric", [d1, d2], []), ("P", [d1], [str(ks[0])]), ("P", [d2], [str(ks[1])])]
            d = self._dummy()
            self._mark(name, [d])
            return [("P", [d], [] if self.propagator else [str(k)]) for k in ks]
        if name in (H["g"].get_name(), H["gbar"].get_name(), H["ghat"].get_name()):
            ix = [lab(x) for x in args]
            self._mark(name, ix)
            return [("Metric", ix, [])]
        if name in (H["gamma"].get_name(), H["gammabar"].get_name(), H["gammahat"].get_name()):
            ix = [lab(x) for x in args]
            self._mark(name, ix[:1])
            return [("Gamma", ix, [])]
        if name == _PSLASH.get_name():
            return [("PSlash", [lab(x) for x in args], [])]
        simple = {H["gamma5"].get_name(): "Gamma5", H["PL"].get_name(): "ProjM", H["PR"].get_name(): "ProjP"}
        if name in simple:
            return [(simple[name], [lab(x) for x in args], [])]
        if name == H["delta"].get_name():
            kind = tn.label_space(args[0])
            return [("Identity" if kind == "spinor" else "Metric", [lab(x) for x in args], [])]
        if name == H["eps"].get_name():
            ix = [lab(x) for x in args]
            for x in ix:
                self.bar[x] = None
            self.sign = -self.sign
            return [("Epsilon", ix, [])]
        raise NotImplementedError(f"Lorentz head {name} is not written")

    def _mark(self, head: str, ix: Sequence[str]) -> None:
        short = head.split("::")[-1]
        if short.endswith("bar"):
            for x in ix:
                self.bar[x] = None
        elif short.endswith("hat"):
            for x in ix:
                self.hat[x] = None

    def render(self, pieces: Sequence[tuple[str, list[str], list[str]]], v2: bool) -> str:
        def arg(x: str) -> str:
            if not v2:
                return x
            if x in self.bar:
                return f"Bar({x})"
            if x in self.hat:
                return f"Hat({x})"
            return x
        out = [f"{h}({','.join([arg(x) for x in ix] + extra)})" for h, ix, extra in pieces]
        if not v2:
            out += [f"bmhv::bar({x})" for x in self.bar] + [f"bmhv::hat({x})" for x in self.hat]
        return "*".join(out) if out else "1"


def _factor_list(term: Expression) -> list[Expression]:
    out = []
    for fa in (list(term) if term.get_type() == AtomType.Mul else [term]):
        if fa.get_type() == AtomType.Pow and next(iter(fa)).get_type() == AtomType.Fn:
            base, ex = list(fa)
            n = int(ex.to_canonical_string())
            if n < 1:
                raise NotImplementedError(f"negative power of a tensor head: {fa}")
            out.extend([base] * n)
        else:
            out.append(fa)
    return out


def _strip_deltas(expr: Expression, internal: set[str]) -> Expression:
    """``expr = R * prod(delta)`` with the Kronecker deltas of the ``internal`` labels common to all terms; returns
    ``R`` (component and colour deltas of a propagator are implicit in HepKIT)."""
    out, common = Expression.num(0), None
    for term in (list(expr) if expr.get_type() == AtomType.Add else [expr]):
        keep, deltas = Expression.num(1), []
        for fa in (list(term) if term.get_type() == AtomType.Mul else [term]):
            if fa.get_type() == AtomType.Fn and fa.get_name() == tn.delta.get_name() and \
                    all(x.to_canonical_string() in internal for x in fa):
                deltas.append(fa.to_canonical_string())
            else:
                keep = keep * fa
        if common is None:
            common = sorted(deltas)
        elif sorted(deltas) != common:
            raise NotImplementedError("propagator not diagonal in its internal indices")
        out = out + keep
    return out


_GROUP_HEADS = {tn.T.get_name(), tn.f.get_name(), tn.dsym.get_name(), tn.epsf.get_name()}


@dataclass(frozen=True)
class _Term:
    number: Expression
    scalar: Expression
    group: Expression
    colour: tuple[Expression, ...]
    lorentz: tuple[Expression, ...]


def _split_term(term: Expression, native: str | None) -> _Term:
    """``number * scalar * group(non-native) * colour(native) * lorentz`` of one product term."""
    number, scalar, group = Expression.num(1), Expression.num(1), Expression.num(1)
    colour, lorentz = [], []
    for fa in _factor_list(term):
        t = fa.get_type()
        if t == AtomType.Num:
            number = number * fa
        elif t == AtomType.Var or (t == AtomType.Pow and next(iter(fa)).get_type() == AtomType.Var):
            scalar = scalar * fa
        elif t == AtomType.Fn:
            name = fa.get_name()
            if name in _GROUP_HEADS or name == tn.delta.get_name():
                space = tn.label_space(next(iter(fa)))
                if space in ("lorentz", "spinor"):
                    lorentz.append(fa)
                elif space != "flavour" and space.split(":")[1] == native:
                    colour.append(fa)
                else:
                    group = group * fa
            elif fa.get_name() in {h.get_name() for h in tn.HEADS.values()}:
                lorentz.append(fa)
            else:
                raise NotImplementedError(f"cannot write the factor {fa} (flavour-matrix or complex parameters)")
        else:
            raise ValueError(f"unexpected factor {fa} (expand first)")
    return _Term(number, scalar, group, tuple(colour), tuple(lorentz))


def _exact_group_values(model: Model, groups: Sequence[Expression], open_labels: Sequence[Expression],
                        ) -> tuple[np.ndarray, np.ndarray]:
    """Exact real and imaginary parts ``(n_terms, *component dims)`` of products of non-native group heads and
    numbers at explicit components: values in ``GF(p)`` from the group evaluator in both embeddings of ``i``
    (the Galois conjugate is :func:`tensors.conjugate`), split and rationally reconstructed over several primes."""
    conj = [tn.conjugate(g) for g in groups]

    def compute(fp):
        ev = model.group_evaluator(fp)
        rows = []
        for g, c in zip(groups, conj):
            plus = np.asarray(ev.evaluate(g, open_labels), dtype=np.int64)
            minus = np.asarray(ev.evaluate(c, open_labels), dtype=np.int64)
            rows.append(np.stack(split_gaussian(plus, minus, fp)))
        return np.stack(rows)

    values = reconstruct_rationals(compute)
    return values[:, 0], values[:, 1]


# ---------------------------------------------------------------------------------------------- layout
@dataclass(frozen=True)
class _Layout:
    """Component slots, native colour, conjugation and names of the written fields."""

    native: str | None
    expanded: Mapping[str, tuple[int, ...]]
    dims: Mapping[str, tuple[int, ...]]
    partner: Mapping[str, str]
    particle_side: Mapping[str, bool]

    def components(self, name: str) -> list[tuple[int, ...]]:
        return list(itertools.product(*[range(n) for n in self.dims[name]]))

    def particle_name(self, name: str, comps: Sequence[int]) -> str:
        base = name if self.particle_side[name] else self.partner[name]
        if any(n >= 10 for n in self.dims[name]):
            suffix = "".join(f"_{c + 1}" for c in comps)
        else:
            suffix = "".join(str(c + 1) for c in comps)
        out = base + suffix
        return out if self.particle_side[name] or self.partner[name] == name else out + "~"


def _native_factor(model: Model, choice: str | None) -> str | None:
    if choice != "auto":
        if choice is not None and choice not in model.factors:
            raise KeyError(f"unknown gauge factor {choice!r}")
        if choice is not None and model.factors[choice].group != "SU(3)":
            raise ValueError("HepKIT colour is SU(3)")
        return choice
    su3 = [n for n, f in model.factors.items() if f.group == "SU(3)"]
    return su3[0] if len(su3) == 1 else None


def _is_particle_side(table: FieldTable, info: FieldInfo) -> bool:
    if info.is_source:
        return not _is_particle_side(table, table[info.source_of])
    return info.kind in _PARTICLE_KINDS


def _layout(model: Model, written: set[str], props, native: str | None) -> _Layout:
    table = model.fields
    partner: dict[str, str] = {}
    for prop in props:
        x, y = prop.fields
        if x not in written and y not in written:
            continue
        if x != y and (x not in written or y not in written):
            kinds = {table[x].kind, table[y].kind}
            if not any(kinds == set(pair) for pair in _CONJUGATE_KINDS):
                continue  # mixing with a field without vertices (b): its lines never occur
        for a, b in ((x, y), (y, x)):
            if partner.setdefault(a, b) != b:
                raise NotImplementedError(f"mixed propagators of {a} with {partner[a]} and {b}")
    for name in sorted(written, key=table.order):
        info = table[name]
        if info.is_source:
            p = partner.get(info.source_of, info.source_of)
            partner[name] = table[p].source if table[p].source is not None else name
            partner[partner[name]] = name
        elif name not in partner:
            raise ValueError(f"field {name} has vertices but no propagator")
    expanded, dims = {}, {}
    all_fields = set(written) | set(partner.values())
    for name in all_fields:
        slots = table[name].slots
        pos = tuple(k for k, s in enumerate(slots)
                    if s.space.kind == "flavour" or (s.factor is not None and s.factor != native))
        natives = [k for k, s in enumerate(slots) if s.factor is not None and s.factor == native]
        if len(natives) > 1:
            raise NotImplementedError(f"{name}: more than one native colour slot")
        expanded[name] = pos
        dims[name] = tuple(slots[k].space.dim for k in pos)
    side = {n: _is_particle_side(table, table[n]) for n in all_fields}
    for n in all_fields:
        if partner[n] != n and side[n] == side[partner[n]]:
            raise ValueError(f"cannot tell particle from antiparticle in the pair {n}, {partner[n]}")
    return _Layout(native, expanded, dims, partner, side)


def _spin(info: FieldInfo) -> int:
    if info.kind in ("ghost", "antighost"):
        return -1
    kinds = [s.space.kind for s in info.slots]
    if "lorentz" in kinds and "spinor" in kinds:
        raise NotImplementedError(f"{info.name}: vector-spinor fields")
    return 3 if "lorentz" in kinds else 2 if "spinor" in kinds else 1


def _colour(info: FieldInfo, native: str | None) -> int:
    reps = [s.rep for s in info.slots if s.factor is not None and s.factor == native]
    return _COLOUR_CODE[reps[0]] if reps else 1


# ---------------------------------------------------------------------------------------------- writer
class _Writer:
    """One run of :func:`write_model`: layout, particles, propagators, then the vertex rules entry by entry and
    class by class (component vertices, exact zero test, coupling assembly), then the document."""

    def __init__(self, model: Model, store: OrderStore, options: WriterOptions) -> None:
        self.model, self.store, self.options = model, store, options
        self.table = _symbol_table(model, options)
        self.subs = _substitution_rules(model, options)
        self.readable = model.simplifier(expand_projectors=False)
        self.fps = [Fingerprinter(model, dimension=dim, seed=options.seed) for dim in options.dimensions]
        self.entries = store.entries(options.max_loop_order)
        self.native = _native_factor(model, options.native_colour)
        self.v2 = options.schema == "v2"
        self.lorentz: dict[tuple[str, tuple[int, ...]], dict] = {}
        self.couplings: dict[tuple[str, tuple], dict] = {}
        self.vertices: list[dict] = []
        self.reverse = ReverseMap()
        self.used_orders: dict[str, dict] = {}
        self.symbol_values = []
        for name, sym in self.table.items():
            if sym.kind in ("regulator", "constant", "symbol"):
                digest = hashlib.sha256(repr((name, options.seed)).encode()).digest()
                symbol = Expression.PI if name == Expression.PI.get_name() else Expression.parse(name)
                self.symbol_values.append((symbol, Expression.num(2 + int.from_bytes(digest[:3], "little"))))

    # ------------------------------------------------------------------ classes
    def vertex_classes(self, entry: StoreEntry) -> list[tuple[str, ...]]:
        classes = entry.leg_classes()
        if entry.role is Role.TREE:
            classes = [c for c in classes if is_vertex_class(self.model.fields, c)]
        return classes

    def run(self) -> WrittenModel:
        model = self.model
        tree = self.store.tree_action()
        props = [p for p in propagators(tree, model, self.readable)] if len(tree) else []
        written = {f for e in self.entries for c in self.vertex_classes(e) for f in c}
        self.layout = _layout(model, written, props, self.native)
        particles = self.particles()
        prop_json = self.propagators(props)
        for entry in self.entries:
            for cls in self.vertex_classes(entry):
                self.vertex_class(entry, cls)
        for name, sym in self.table.items():
            if sym.kind in ("coupling", "parameter", "gauge_parameter") and sym.ufo not in self.options.substitutions:
                self.reverse.parameters[sym.ufo] = (sym.kind, sym.ref)
            elif sym.kind in ("regulator", "symbol"):
                self.reverse.parameters[sym.ufo] = (sym.kind, name)
        data = self.document(particles, prop_json)
        return WrittenModel(data, self.reverse, model, self.store, self.options)

    # ------------------------------------------------------------------ particles
    def particles(self) -> list[dict]:
        table, lay = self.model.fields, self.layout
        out = []
        pdg = 9000000
        names = sorted((n for n in lay.partner if lay.particle_side[n] or lay.partner[n] == n), key=table.order)
        for name in names:
            for comps in lay.components(name):
                pdg += 1
                out.append(self.particle(name, comps, pdg))
                if lay.partner[name] != name:
                    out.append(self.particle(lay.partner[name], comps, -pdg))
        return out

    def particle(self, name: str, comps: tuple[int, ...], code: int) -> dict:
        table, lay = self.model.fields, self.layout
        info = table[name]
        pname = lay.particle_name(name, comps)
        anti = lay.particle_name(lay.partner[name], comps)
        slotc = tuple(zip(lay.expanded[name], comps))
        self.reverse.particles[pname] = ParticleRecord(pname, name, slotc, anti, code)
        reps = {}
        for s in info.slots:
            if s.factor is not None:
                reps[s.factor] = s.rep
        for fac, q in info.charges.items():
            reps[fac] = str(q)
        ghost = info.ghost_number
        rec = {
            "pdg_code": code, "name": pname, "antiname": anti, "spin": _spin(info), "color": _colour(info, lay.native),
            "mass": "ZERO", "width": "ZERO", "texname": pname, "antitexname": anti, "charge": 0.0,
            "ghost_number": ghost if self.v2 else 0, "lepton_number": 0, "y_charge": 0,
            "propagating": not info.is_source,
            "statistics": "fermionic" if info.odd else "bosonic",
            "field_kind": _FIELD_KIND[info.kind],
            "representations": reps,
            "multiplet": {"field": name, "multiplet": info.multiplet or info.factor or name,
                          "components": {str(table[name].slots[s].factor or "flavour"): c + 1 for s, c in slotc}},
        }
        if info.is_source:
            rec["source_of"] = lay.particle_name(info.source_of, comps)
        if not self.v2 and ghost != 0:
            rec["v2_overrides"] = {"ghost_number": ghost}
        return rec

    # ------------------------------------------------------------------ propagators
    def propagators(self, props) -> list[dict]:
        lay = self.layout
        by_fields = {p.fields: p for p in props}
        out = []
        for name in sorted(lay.partner, key=self.model.fields.order):
            if not lay.particle_side[name] and lay.partner[name] != name:
                continue
            if self.model.fields[name].is_source:
                continue
            prop = by_fields[(name, lay.partner[name])]
            strings = {False: self.propagator_strings(prop, False)}
            if lay.partner[name] != name:
                strings[True] = self.propagator_strings(prop, True)
            for comps in lay.components(name):
                for rev, (num, den) in strings.items():
                    particle = lay.particle_name(lay.partner[name] if rev else name, comps)
                    pname = f"{particle}_prop"
                    out.append({"name": pname, "particle": particle, "numerator": num, "denominator": den})
                    self.reverse.propagators[pname] = PropagatorRecord(pname, particle, prop.fields, rev)
        return out

    def propagator_strings(self, prop, reverse: bool) -> tuple[str, str]:
        table = self.model.fields
        x, y = prop.fields
        lx, ly = prop.labels
        a, s_, t_ = _w["a"], _w["s"], _w["t"]
        value = _substitute(prop.value, self.subs)
        value = value.replace(tn.gammahat(a, s_, t_), tn.gamma(a, s_, t_) - tn.gammabar(a, s_, t_))
        value = value.replace(tn.sphat(s_, t_), tn.sp(s_, t_) - tn.spbar(s_, t_))
        value = value.replace(tn.ghat(s_, t_), tn.g(s_, t_) - tn.gbar(s_, t_))
        if reverse:
            value = value.replace(tn.p(0, a), -tn.p(0, a))
        internal = {lab.to_canonical_string() for labels, f in ((lx, x), (ly, y))
                    for lab, sl in zip(labels, table[f].slots) if sl.space.kind not in ("lorentz", "spinor")}
        value = _strip_deltas(value.expand(), internal).together()
        value = value.replace(tn.p(0, a) * tn.gamma(a, s_, t_), _PSLASH(s_, t_))
        x_end, y_end = ("idx(1,1)", "idx(1,2)") if reverse else ("idx(1,2)", "idx(1,1)")
        index = {lab.to_canonical_string(): x_end for lab in lx}
        index.update({lab.to_canonical_string(): y_end for lab in ly})
        num, den = Expression.num(1), []
        for fa in (list(value) if value.get_type() == AtomType.Mul else [value]):
            if fa.get_type() == AtomType.Pow and list(fa)[1].get_type() == AtomType.Num and \
                    Fraction(list(fa)[1].to_canonical_string()) < 0:
                den.append((next(iter(fa)), -Fraction(list(fa)[1].to_canonical_string())))
            else:
                num = num * fa
        if len(den) == 1 and den[0][1] == 1:  # write k^2 - m^2, not m^2 - k^2
            c = den[0][0].expand().coefficient(tn.sp(0, 0))
            if c.get_type() == AtomType.Num and complex_rational(c)[0] < 0:
                num, den = -num, [(-den[0][0], 1)]
        return self.numerator_string(num.expand(), index), self.denominator_string(den)

    def numerator_string(self, expr: Expression, index: Mapping[str, str]) -> str:
        """Propagator numerator: a sum of numbers x symbols x structures."""
        if expr == Expression.num(0):
            raise ArithmeticError("vanishing propagator")
        parts = []
        for term in (list(expr) if expr.get_type() == AtomType.Add else [expr]):
            scal, heads = Expression.num(1), []
            for fa in _factor_list(term):
                if fa.get_type() == AtomType.Fn:
                    heads.append(fa)
                else:
                    scal = scal * fa
            pr = _Printer(index, propagator=True)
            pieces = [pc for h in heads for pc in pr.pieces(h)]
            if pr.bar or pr.hat:
                raise NotImplementedError("barred or hatted propagator structures")
            body = pr.render(pieces, False)
            parts.append(_scalar_string(scal * pr.sign, self.table) + ("" if body == "1" else "*" + body))
        return "+".join(parts)

    def denominator_string(self, factors: Sequence[tuple[Expression, Fraction]]) -> str:
        """Propagator denominator: products of powers of polynomials in ``k^2`` (written ``P(idx(1,1))^2``, the
        hand-written model's form) and symbols."""
        if not factors:
            return "1"
        k2 = S("hepkit_renorm_output::k2")
        out = []
        for base, n in factors:
            base = base.replace(tn.sp(0, 0), k2).expand()
            if any(fa.get_type() == AtomType.Fn for t in (list(base) if base.get_type() == AtomType.Add else [base])
                   for fa in _factor_list(t)):
                raise NotImplementedError(f"propagator denominator {base}")
            text = _scalar_string(base, {**self.table, k2.get_name(): _Symbol("(P(idx(1,1)))^2", "", "", None)})
            out.append(f"({text})" + ("" if n == 1 else f"^({n})"))
        return "*".join(out)

    # ------------------------------------------------------------------ vertices
    def written_role(self, entry: StoreEntry, cls: Sequence[str]) -> str:
        if entry.role is Role.TREE and any(self.model.fields[f].is_source for f in cls):
            return "source"
        return entry.role.value

    def role_orders(self, entry: StoreEntry, cls: Sequence[str]) -> dict[str, int]:
        on = self.options.orders
        out: dict[str, int] = {}
        nsrc = sum(self.model.fields[f].is_source for f in cls)
        if nsrc:
            out[on.source] = nsrc
        if entry.role is Role.BREAKING:
            out[on.breaking] = 1
        if entry.loop_order:
            out[on.loop] = entry.loop_order
        if entry.role is Role.CT_SINGULAR:
            out[on.singular] = 1
        if entry.role is Role.CT_FINITE:
            out[on.finite] = 1
        return out

    def vertex_class(self, entry: StoreEntry, cls: tuple[str, ...]) -> None:
        model, lay = self.model, self.layout
        table = model.fields
        rules = {}
        for k in entry.powers:
            rule = vertex_rule(entry.coefficient(k), cls, self.readable, eliminate=None)
            value = _substitute(rule.value, self.subs)
            group_labels = [lab for leg in rule.legs for lab, sl in zip(leg.labels, table[leg.field].slots)
                            if sl.factor is not None]
            # classes that vanish by the Lie algebra (e.g. four-point classes of s S_0) are dropped exactly
            if value != Expression.num(0) and not model.vanishes_modulo_group(value, group_labels, default_field()):
                rules[k] = replace(rule, value=value)
        if not rules:
            return
        legs = next(iter(rules.values())).legs
        open_labels = [leg.labels[s] for leg in legs for s in lay.expanded[leg.field]]
        decomposed = {}
        for k, rule in rules.items():
            terms = [_split_term(t, lay.native) for t in
                     (list(rule.value) if rule.value.get_type() == AtomType.Add else [rule.value])]
            re, im = _exact_group_values(model, [t.number * t.group for t in terms], open_labels)
            decomposed[k] = (terms, re, im)
        leg_index = {}
        colour_index = {}
        for j, leg in enumerate(legs, start=1):
            for lab, s in zip(leg.labels, table[leg.field].slots):
                key = lab.to_canonical_string()
                if s.space.kind in ("lorentz", "spinor"):
                    leg_index[key] = f"idx(1,{j})"
                elif s.factor is not None and s.factor == lay.native:
                    colour_index[key] = str(j)
        spins = tuple(_spin(table[f]) for f in cls)
        role = self.written_role(entry, cls)
        base_orders = self.role_orders(entry, cls)
        per_leg = [lay.components(f) for f in cls]
        for choice in itertools.product(*per_leg):
            if any(cls[i] == cls[j] and choice[i] > choice[j] for i in range(len(cls)) for j in range(i + 1, len(cls))):
                continue
            flat = tuple(c for leg_c in choice for c in leg_c)
            comps = tuple(tuple(zip(lay.expanded[f], c)) for f, c in zip(cls, choice))
            # a component whose exact group factors all vanish is zero; the others are fingerprinted
            alive = [k for k, rule in rules.items()
                     if any(decomposed[k][1][(n,) + flat] or decomposed[k][2][(n,) + flat]
                            for n in range(len(decomposed[k][0]))) and self.nonzero(rule, comps)]
            if not alive:
                continue
            acc: dict[tuple, dict] = {}
            evanescent = False
            for k in alive:
                terms, re, im = decomposed[k]
                for n, term in enumerate(terms):
                    c, ci = Fraction(re[(n,) + flat]), Fraction(im[(n,) + flat])
                    if c == 0 and ci == 0:
                        continue
                    number = Expression.num(c.numerator) / c.denominator + \
                        tn.I * Expression.num(ci.numerator) / ci.denominator
                    orders = dict(base_orders)
                    for o, m in _orders_of(term.scalar, self.table).items():
                        orders[o] = orders.get(o, 0) + m
                    cpr = _Printer(colour_index, colour=True)
                    colour = cpr.render([pc for fa in term.colour for pc in cpr.pieces(fa)], False)
                    lname, ev, sign = self.lorentz_structure(term.lorentz, leg_index, spins)
                    evanescent = evanescent or ev
                    coeff = (tn.I * number * sign * term.scalar * _E**k).replace(tn.d, 4 - 2 * _E)
                    key = (tuple(sorted(orders.items())), colour, lname)
                    acc[key] = (acc[key] + coeff) if key in acc else coeff
            self.emit(entry, cls, choice, comps, role, acc, evanescent)

    def nonzero(self, rule: VertexRule, comps) -> bool:
        """Exact zero test at the components: fingerprints in every ``D`` of ``options.dimensions``; symbols that are
        not model parameters (the tag, ``pi``, user symbols) take deterministic random integer values
        (Schwartz-Zippel)."""
        value = rule.value
        for sym, val in self.symbol_values:
            value = value.replace(sym, val)
        value = value.expand()
        for fp in self.fps:
            smp = component_samples(fp, rule.fields, comps, self.options.samples, salt=1)
            if np.any(fp.evaluate(value, smp)):
                return True
        return False

    def lorentz_structure(self, factors: Sequence[Expression], leg_index: Mapping[str, str],
                          spins: tuple[int, ...]) -> tuple[str, bool, int]:
        pr = _Printer(leg_index)
        pieces = [pc for fa in factors for pc in pr.pieces(fa)]
        v1, v2 = pr.render(pieces, False), pr.render(pieces, True)
        key = (v1, spins)
        if key not in self.lorentz:
            name = f"L{len(self.lorentz) + 1}"
            rec = {"name": name, "spins": list(spins), "structure": v2 if self.v2 else v1}
            if not self.v2 and v2 != v1:
                rec["v2_overrides"] = {"structure": v2}
            self.lorentz[key] = rec
        return self.lorentz[key]["name"], bool(pr.hat), pr.sign

    def coupling(self, expr: Expression, orders: tuple[tuple[str, int], ...], laurent: bool) -> str:
        text = _scalar_string(expr, self.table)
        key = (text, orders)
        if key not in self.couplings:
            rec = {"name": f"GC_{len(self.couplings) + 1}", "expression": text,
                   "orders": [list(o) for o in orders], "value": None}
            if laurent:
                rec["laurent"] = {"variable": self.options.regulator,
                                  "terms": {str(k): _scalar_string(c, self.table)
                                            for k, c in sorted(_laurent(expr).items())}}
            self.couplings[key] = rec
            for o, _ in orders:
                self.used_orders.setdefault(o, {})
        return self.couplings[key]["name"]

    def emit(self, entry, cls, choice, comps, role, acc, evanescent) -> None:
        lay = self.layout
        particles = [lay.particle_name(f, c) for f, c in zip(cls, choice)]
        by_orders: dict[tuple, dict] = {}
        for (orders, colour, lname), expr in acc.items():
            expr = expr.expand()
            if expr == Expression.num(0):
                continue
            by_orders.setdefault(orders, {})[(colour, lname)] = expr
        if not by_orders:
            raise ArithmeticError(f"component vertex {particles} of {entry.name} is non-zero but its terms cancel")
        laurent = entry.loop_order > 0 or entry.powers != (0,)
        base = entry.name + "_" + "_".join(p.replace("~", "bar") for p in particles)
        for n, (orders, entries) in enumerate(sorted(by_orders.items())):
            colours = list(dict.fromkeys(c for c, _ in entries))
            lorentz = list(dict.fromkeys(lz for _, lz in entries))
            matrix = [[self.coupling(entries[(c, lz)], orders, laurent) if (c, lz) in entries else None
                       for lz in lorentz] for c in colours]
            name = base if len(by_orders) == 1 else f"{base}_{n + 1}"
            odd = [k for k, f in enumerate(cls) if self.model.fields[f].odd]
            self.vertices.append({
                "name": name, "particles": particles, "color_structures": colours, "lorentz_structures": lorentz,
                "couplings": matrix, "loop_order": entry.loop_order, "role": role, "operator": entry.name,
                "grassmann_order": odd, "evanescent": evanescent})
            self.reverse.vertices[name] = VertexRecord(name, entry.name, entry.loop_order, role, tuple(cls), comps,
                                                       tuple(particles), orders, entry.powers, evanescent)

    # ------------------------------------------------------------------ document
    def document(self, particles: list[dict], props: list[dict]) -> dict:
        model, opt = self.model, self.options
        on = opt.orders
        kinds = {on.breaking: "insertion", on.source: "source", on.loop: "loop_counterterm",
                 on.singular: "role", on.finite: "role"}
        orders = []
        for k, name in enumerate(sorted(self.used_orders)):
            kind = kinds.get(name, "coupling")
            orders.append({"name": name, "expansion_order": 99, "hierarchy": k + 1, "kind": kind})
            self.reverse.orders[name] = {"kind": kind}
        params = [{"name": "ZERO", "lhablock": None, "lhacode": None, "nature": "internal", "parameter_type": "real",
                   "value": [0.0, 0.0], "expression": None}]
        n = 0
        for sym in self.table.values():
            if sym.kind not in ("coupling", "parameter", "gauge_parameter") or sym.ufo in opt.substitutions:
                continue
            n += 1
            real = sym.kind != "parameter" or model.parameters[sym.ref].real
            params.append({"name": sym.ufo, "lhablock": "HEPKITRENORM", "lhacode": [n], "nature": "external",
                           "parameter_type": "real" if real else "complex", "value": [1.0, 0.0],
                           "expression": None})
        operators = [{"name": e.name, "loop_order": e.loop_order, "role": e.role.value, "description": e.description,
                      "laurent_powers": list(e.powers),
                      "vertex_rules": [v["name"] for v in self.vertices if v["operator"] == e.name]}
                     for e in self.entries]
        data = {
            "name": opt.name or f"{model.name}_bmhv", "restriction": None, "orders": orders, "parameters": params,
            "particles": particles, "propagators": props, "lorentz_structures": list(self.lorentz.values()),
            "couplings": list(self.couplings.values()), "vertex_rules": self.vertices,
            "regularization": {"scheme": model.regularisation.scheme, "epsilon": opt.regulator,
                               "dimension": f"4-2*{opt.regulator}", "evanescent_tag": opt.tag,
                               "renormalization_scheme": "MSbar", "sterile_partners":
                                   model.regularisation.sterile_partners},
            "gauge_groups": [{"name": n_, "kind": "U1" if fa.abelian else "SU", "n": fa.n, "coupling": fa.coupling,
                              "native": "color" if n_ == self.native else None} for n_, fa in model.factors.items()],
            "operators": operators,
            "hepkit_renorm": {"writer": "hepkit_renorm.output", "schema": opt.schema,
                              "conventions": "vertex = i V (package leg order); signs from the package, not from "
                                             "the generator; P(mu,k) incoming momentum of leg k",
                              "reverse_map": self.reverse.to_dict()},
        }
        if self.v2:
            data["schema_version"] = 2
        return data


def write_model(model: Model, store: OrderStore, options: WriterOptions | None = None) -> WrittenModel:
    """Write the entries of ``store`` (up to ``options.max_loop_order``) as a HepKIT JSON model (module
    docstring) with its reverse map."""
    return _Writer(model, store, options or WriterOptions()).run()


def tree_model(model: Model, options: WriterOptions | None = None, brst=None) -> WrittenModel:
    """The tree-level regularised model: ``S_0`` (with gauge fixing, ghosts and sources) and ``Delta_1`` of
    :func:`hepkit_renorm.store.tree_store`, in Feynman gauge (``xi = 1``) unless ``options`` say otherwise."""
    if options is None:
        xi = model.gauge_fixing.xi
        subs = {xi.get_name().split("::")[-1]: 1} if xi.get_type() == AtomType.Var else {}
        options = WriterOptions(substitutions=subs)
    return write_model(model, tree_store(model, brst), options)

