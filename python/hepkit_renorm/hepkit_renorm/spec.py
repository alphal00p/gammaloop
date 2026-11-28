"""Model specification: gauge group, matter, classical Lagrangian, gauge fixing and regularisation.

A specification is a small declarative document (TOML; see ``models/sm_one_family.toml``)::

    name = "..."
    generations = 1                         # flavour multiplicity (one flavour index space)

    [gauge]                                 # factors in order; the factor name is also its gauge field name
    G = { group = "SU(3)", coupling = "g_s" }
    B = { group = "U(1)",  coupling = "gp" }

    [fermions]                              # chiral multiplets: physical chirality, reps per factor, U(1) charges
    Q = { chirality = "L", G = "F", W = "F", B = "1/6" }

    [scalars]
    H = { complex = true, W = "F", B = "1/2" }

    [parameters]                            # Lagrangian parameters besides the gauge couplings (real by default)
    y_u = { flavour = true }                # flavour = true: a matrix in flavour space when generations > 1

    [[potential]]                           # V = sum factor * coupling * prod (phi_i^dagger phi_j); L contains -V
    coupling = "lambda"
    bilinears = [["H", "H"], ["H", "H"]]

    [[yukawa]]                              # L contains -coupling * fbar phi f + h.c.
    coupling = "y_u"
    fermion_bar = "Q"
    scalar = "Hdag"                         # a scalar field or its conjugate
    fermion = "U"

    [gauge_fixing]
    kind = "linear_covariant"
    xi = "xi"                               # a parameter name or a number

    [regularisation]
    scheme = "BMHV"
    sterile_partners = true
    fermion_kinetic = "d_dimensional"
    gauge_vertex = "barred_physical"
    bosonic = "d_dimensional"

Conventions (report 03_setup): ``D_mu = d_mu - i g_r T^a A^a_mu`` on every field; gauge kinetic term
``-1/4 F^a F^a`` with ``F = dA - dA + fg A A``; scalar kinetic term ``(D phi)^dagger D phi``; fermion kinetic term
``Psibar i gamma.d Psi`` (d-dimensional) and gauge vertex ``Psibar gbar^mu g_r T^a P_c Psi A^a_mu``. The group
contraction of a Yukawa term is derived per factor: a pair (F, Fbar) or (A, A) is contracted with ``delta``, a pair
of two F (or two Fbar) of SU(2) with ``epsf`` (``epsf(1,2) = +1``, first-listed index first); the U(1) charges must
add up to zero and the two chiral projectors must agree. Hermitian conjugate terms are generated, never typed.

:class:`Model` is the validated specification together with the derived field registry
(:class:`hepkit_renorm.fields.FieldTable`): gauge fields, ghosts, antighosts, auxiliary fields, matter fields and
their conjugates, and one BRST source for every gauge field, ghost and matter field.
"""

from __future__ import annotations

import tomllib
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
from symbolica import AtomType, Expression, S

from . import tensors as tn
from .exact import PrimeField
from .fields import FieldInfo, FieldTable, Slot
from .functional import Density, Monomial
from .groups import (
    GaugeFactor,
    GroupEvaluator,
    Rep,
    casimir_value,
    conjugate_rep,
    dynkin_index,
    product_generators,
    rep_generators,
)

__all__ = [
    "FermionMultiplet",
    "GaugeFixing",
    "Model",
    "Parameter",
    "PotentialTerm",
    "Regularisation",
    "ScalarMultiplet",
    "YukawaTerm",
    "builtin_model",
    "load_model",
    "model_from_dict",
]

_PARAM_NS = f"{tn.NS}_param::"
_CPARAM_NS = f"{tn.NS}_cparam::"


@dataclass(frozen=True)
class FermionMultiplet:
    """A chiral fermion multiplet ``f``: physical chirality ``"L"``/``"R"``, representation per gauge factor
    (U(1) entries are charges), ``flavoured``: carries the flavour index."""

    name: str
    chirality: str
    reps: Mapping[str, Rep]
    flavoured: bool = True


@dataclass(frozen=True)
class ScalarMultiplet:
    """A scalar multiplet; ``complex`` scalars have a conjugate field ``<name>dag``."""

    name: str
    complex: bool
    reps: Mapping[str, Rep]
    flavoured: bool = False


@dataclass(frozen=True)
class Parameter:
    """A Lagrangian parameter. ``flavour``: a matrix over the flavour indices of the fields it couples (only
    meaningful when ``generations > 1``)."""

    name: str
    real: bool = True
    flavour: bool = False


@dataclass(frozen=True)
class PotentialTerm:
    """``factor * coupling * prod_k (phi_{i_k}^dagger phi_{j_k})`` in ``V`` (``L`` contains ``-V``)."""

    coupling: str
    bilinears: tuple[tuple[str, str], ...]
    factor: Fraction = Fraction(1)


@dataclass(frozen=True)
class YukawaTerm:
    """``-coupling * fermion_bar scalar fermion + h.c.`` in ``L``. ``scalar`` is a scalar field or its
    conjugate (``Hdag``)."""

    coupling: str
    fermion_bar: str
    scalar: str
    fermion: str


@dataclass(frozen=True)
class GaugeFixing:
    """Linear covariant gauge ``L_gf = b^a d.A^a + xi/2 b^a b^a + ...`` (report eq:Lgf); ``xi`` a coefficient."""

    kind: str
    xi: Expression


@dataclass(frozen=True)
class Regularisation:
    """Regularisation options (the report's choices are the only ones implemented; others are rejected)."""

    scheme: str = "BMHV"
    sterile_partners: bool = True
    fermion_kinetic: str = "d_dimensional"
    gauge_vertex: str = "barred_physical"
    bosonic: str = "d_dimensional"

    def __post_init__(self) -> None:
        allowed = {"scheme": "BMHV", "sterile_partners": True, "fermion_kinetic": "d_dimensional",
                   "gauge_vertex": "barred_physical", "bosonic": "d_dimensional"}
        for k, v in allowed.items():
            if getattr(self, k) != v:
                raise NotImplementedError(f"regularisation option {k} = {getattr(self, k)!r} is not implemented")


@dataclass
class Model:
    """A validated model specification and its derived field registry.

    Attributes: ``factors`` (ordered ``{name: GaugeFactor}``), ``fermions``, ``scalars``, ``parameters``,
    ``potential``, ``yukawas``, ``gauge_fixing``, ``regularisation``, ``generations`` and ``fields`` (the
    :class:`FieldTable`).
    """

    name: str
    factors: dict[str, GaugeFactor]
    fermions: dict[str, FermionMultiplet]
    scalars: dict[str, ScalarMultiplet]
    parameters: dict[str, Parameter]
    potential: tuple[PotentialTerm, ...]
    yukawas: tuple[YukawaTerm, ...]
    gauge_fixing: GaugeFixing
    regularisation: Regularisation
    generations: int = 1
    fields: FieldTable = field(init=False)

    def __post_init__(self) -> None:
        self._validate()
        self.fields = FieldTable(self._derive_fields())

    # ------------------------------------------------------------------ parameters
    def param(self, name: str) -> Expression:
        """The Symbolica symbol of a parameter or coupling (real unless declared complex)."""
        if name in self.parameters:
            real = self.parameters[name].real
        elif name in (f.coupling for f in self.factors.values()):
            real = True
        else:
            raise KeyError(f"unknown parameter {name!r}")
        # real and complex parameters live in different namespaces: a Symbolica symbol keeps its attributes
        return S(_PARAM_NS + name, is_real=True) if real else S(_CPARAM_NS + name)

    def coupling(self, factor: str) -> Expression:
        """The coupling symbol ``g_r`` of a gauge factor."""
        return self.param(self.factors[factor].coupling)

    def flavour_tensor(self, name: str, i: Expression | None, j: Expression | None) -> Expression:
        """A flavour-matrix parameter ``y(i, j)``: the plain symbol if the flavour slots are absent."""
        sym = self.param(name)
        if i is None and j is None:
            return sym
        if i is None or j is None:
            raise ValueError("flavour matrix with one missing index")
        return sym(i, j)

    # ------------------------------------------------------------------ tools
    def simplifier(self, expand_projectors: bool = True) -> tn.Simplifier:
        """A :class:`tensors.Simplifier` that knows the dimensions of the model's index spaces."""
        dims = {n: (f.n, f.dim_adjoint) for n, f in self.factors.items() if not f.abelian}
        return tn.Simplifier(dims, self.generations, expand_projectors)

    def group_evaluator(self, field_: PrimeField) -> GroupEvaluator:
        return GroupEvaluator(self.factors, field_, self.generations)

    def slot_labels(self, field_name: str) -> tuple[Expression, ...]:
        """Fresh labels for the slots of a field."""
        return tuple(tn.fresh(s.space) for s in self.fields[field_name].slots)

    def multiplet_reps(self, multiplet: str) -> Mapping[str, Rep]:
        if multiplet in self.fermions:
            return self.fermions[multiplet].reps
        return self.scalars[multiplet].reps

    def gauge_factor_names(self) -> list[str]:
        return list(self.factors)

    def matter_fields(self) -> list[str]:
        """Names of the matter fields (fermions, their conjugates, scalars and conjugates)."""
        return [n for n, fi in self.fields.items() if fi.kind in ("fermion", "fermion_bar", "scalar", "scalar_conj")]

    # ------------------------------------------------------------------ group theory of the content
    def multiplet_generators(self, multiplet: str, field_: PrimeField,
                             couplings: Mapping[str, int] | None = None) -> np.ndarray:
        """Coupling-weighted generators ``cal(T)^a`` of a multiplet on its group components (one generation),
        combined adjoint index over all factors (:func:`groups.product_generators`); ``couplings`` gives numerical
        values of the factor couplings (default 1)."""
        reps = self._full_reps(multiplet)
        facs = list(self.factors.values())
        cs = None if couplings is None else [couplings[f.name] for f in facs]
        return product_generators(facs, [reps[f.name] for f in facs], field_, cs)

    def group_invariants(self, field_: PrimeField) -> dict[str, dict]:
        """The invariants of report eq:invariants as Symbolica expressions in the couplings (one generation):

        * ``"S"``: ``{(multiplet, factor): S}`` with ``Tr(cal(T)^a cal(T)^b) = S delta^{ab}`` (trace over all group
          components of the multiplet, for fermions and scalars);
        * ``"C"``: ``{multiplet: C}`` with ``sum_b cal(T)^b cal(T)^b = C 1``;
        * ``"C_ad"``: ``{factor: C_ad}`` with ``fg^{acd} fg^{bcd} = C_ad delta^{ab}``.

        Values are computed from the explicit generators over ``field_`` and rationally reconstructed.
        """
        out: dict[str, dict] = {"S": {}, "C": {}, "C_ad": {}}
        for name in list(self.fermions) + list(self.scalars):
            reps = self._full_reps(name)
            dims = {n: fac.rep_dim(reps[n]) for n, fac in self.factors.items()}
            cas = Expression.num(0)
            for n, fac in self.factors.items():
                g2 = self.coupling(n) ** 2
                other = 1
                for m_, dm in dims.items():
                    if m_ != n:
                        other *= dm
                out["S"][(name, n)] = (g2 * _num_expr(dynkin_index(fac, reps[n], field_) * other)).expand()
                cas = cas + g2 * _num_expr(casimir_value(fac, reps[n], field_))
            out["C"][name] = cas.expand()
        for n, fac in self.factors.items():
            val = Fraction(0) if fac.abelian else casimir_value(fac, "A", field_)
            out["C_ad"][n] = (self.coupling(n) ** 2 * _num_expr(val)).expand()
        return out

    def anomaly_tensor(self, field_: PrimeField) -> np.ndarray:
        """``D^{abc} = sum_lh Tr(T^a {T^b, T^c}) - sum_rh Tr(T^a {T^b, T^c})`` (report eq:Dabc) at unit couplings,
        summed over generations, as an explicit array over the combined adjoint index."""
        total = None
        for name, fm in self.fermions.items():
            t = self.multiplet_generators(name, field_)
            tt = field_.einsum("bij,cjk->bcik", t, t)
            anti = field_.reduce(tt + np.transpose(tt, (1, 0, 2, 3)))
            dabc = field_.einsum("aki,bcik->abc", t, anti)
            sign = (1 if fm.chirality == "L" else -1) * (self.generations if fm.flavoured else 1)
            total = dabc * sign if total is None else total + dabc * sign
        return field_.reduce(total)

    def rigid_defects(self, density: Density, field_: PrimeField) -> list[Monomial]:
        """Monomials of ``density`` that are not invariant under the global gauge group, monomial by monomial
        (sufficient for invariance): U(1) charges must add up to zero, and the group tensor of the coefficient
        (group heads evaluated over ``field_``, indexed by the group slots of the occurrences) must be annihilated
        by ``sum_k T^a`` acting on the slots ``k`` in their representations, for every generator ``a``."""
        ev = self.group_evaluator(field_)
        bad = []
        for mono in density.monomials():
            sep = mono.separated()
            for n, fac in self.factors.items():
                if fac.abelian:
                    q = sum((self.fields[o.field].charges.get(n, Fraction(0)) for o in sep.occurrences), Fraction(0))
                    if q != 0:
                        bad.append(mono)
                        break
            else:
                if not self._group_tensor_invariant(sep, ev, field_):
                    bad.append(mono)
        return bad

    def vanishes_modulo_group(self, expr: Expression, open_labels: Sequence[Expression], field_: PrimeField) -> bool:
        """True if a (simplified, canonical) expression vanishes once its group tensors are evaluated: the terms
        are collected by their non-group part (everything except group heads and numbers) and the group parts of
        each collection are summed as explicit arrays over ``field_``, indexed by the free group labels
        ``open_labels`` (e.g. the group slots of the legs of a vertex function). Decides identities such as
        ``[T^a, T^b] = i f^{abc} T^c`` that the symbolic simplifier does not know."""
        ev = self.group_evaluator(field_)
        collected: dict[str, np.ndarray] = {}
        terms = list(expr) if expr.get_type() == AtomType.Add else [expr]
        for t in terms:
            if t.get_type() == AtomType.Num and t.to_canonical_string() == "0":
                continue
            group_part, rest = Expression.num(1), Expression.num(1)
            for fa in (list(t) if t.get_type() == AtomType.Mul else [t]):
                if fa.get_type() == AtomType.Num or _is_group_factor(fa):
                    group_part = group_part * fa
                else:
                    rest = rest * fa
            key = rest.to_canonical_string()
            val = ev.evaluate(group_part, open_labels)
            collected[key] = val if key not in collected else field_.reduce(collected[key] + val)
        return all(field_.is_zero(v) for v in collected.values()) if collected else True

    def _group_tensor_invariant(self, mono, ev: GroupEvaluator, field_: PrimeField) -> bool:
        open_labels, slot_info = [], []
        for o in mono.occurrences:
            for lab, slot in zip(o.slots, self.fields[o.field].slots):
                if slot.factor is not None:
                    open_labels.append(lab)
                    slot_info.append(slot)
        group_part = Expression.num(1)
        factors = list(mono.coefficient) if mono.coefficient.get_type() == AtomType.Mul else [mono.coefficient]
        for fa in factors:
            if _is_group_factor(fa):
                group_part = group_part * fa
        t = ev.evaluate(group_part, open_labels)
        for n, fac in self.factors.items():
            if fac.abelian:
                continue
            positions = [k for k, sl in enumerate(slot_info) if sl.factor == n]
            if not positions:
                continue
            for a in range(fac.dim_adjoint):
                acc = field_.zeros(t.shape)
                for k in positions:
                    gen = rep_generators(fac, slot_info[k].rep, field_)[a]
                    moved = np.tensordot(t, gen, axes=([k], [0]))   # index j of T appended last
                    acc = acc + np.moveaxis(moved, -1, k)
                if not field_.is_zero(field_.reduce(acc)):
                    return False
        return True

    # ------------------------------------------------------------------ validation
    def _validate(self) -> None:
        names = list(self.factors) + list(self.fermions) + list(self.scalars)
        if len(set(names)) != len(names):
            raise ValueError("factor, fermion and scalar names must be distinct")
        for f in self.factors.values():
            if f.coupling in self.parameters:
                raise ValueError(f"coupling {f.coupling} is also listed as a parameter")
        for m in list(self.fermions.values()) + list(self.scalars.values()):
            for fac in m.reps:
                if fac not in self.factors:
                    raise ValueError(f"{m.name}: unknown gauge factor {fac}")
        for f in self.fermions.values():
            if f.chirality not in ("L", "R"):
                raise NotImplementedError(f"{f.name}: only chiral (L/R) multiplets are implemented")
        for t in self.potential:
            self._param_or_raise(t.coupling)
            for a, b in t.bilinears:
                if a not in self.scalars or b not in self.scalars:
                    raise ValueError(f"potential bilinear ({a}, {b}) is not made of scalars")
                if self._full_reps(a) != self._full_reps(b):
                    raise ValueError(f"bilinear {a}^dagger {b} is not invariant")
        for y in self.yukawas:
            self._param_or_raise(y.coupling)
            self.yukawa_structure(y)
        if self.generations < 1:
            raise ValueError("generations >= 1")

    def _param_or_raise(self, name: str) -> None:
        if name not in self.parameters:
            raise ValueError(f"parameter {name!r} is not declared")

    def _full_reps(self, multiplet: str) -> dict[str, Rep]:
        reps = self.multiplet_reps(multiplet)
        out: dict[str, Rep] = {}
        for n, fac in self.factors.items():
            r = reps.get(n, Fraction(0) if fac.abelian else "1")
            out[n] = Fraction(r) if fac.abelian else r
        return out

    def yukawa_structure(self, y: YukawaTerm) -> tuple[str, list[tuple[str, str, list[tuple[int, int]]]]]:
        """Derived structure of a Yukawa term: the projector (``"PL"``/``"PR"``) between ``fermion_bar`` and
        ``fermion``, and per non-abelian factor the contraction ``(factor, "delta" | "epsf", [(field, slot)])``
        with field 0 = fermion_bar, 1 = scalar, 2 = fermion and ``slot`` the position of the group slot."""
        if y.fermion_bar not in self.fermions or y.fermion not in self.fermions:
            raise ValueError(f"Yukawa {y.coupling}: unknown fermion")
        fb, ff = self.fermions[y.fermion_bar], self.fermions[y.fermion]
        if fb.chirality == ff.chirality:
            raise ValueError(f"Yukawa {y.coupling}: fbar_{fb.name} P f_{ff.name} vanishes (equal chiralities)")
        proj = "PL" if ff.chirality == "L" else "PR"
        conj_scalar = y.scalar.endswith("dag") and y.scalar[:-3] in self.scalars and self.scalars[y.scalar[:-3]].complex
        sname = y.scalar[:-3] if conj_scalar else y.scalar
        if sname not in self.scalars:
            raise ValueError(f"Yukawa {y.coupling}: unknown scalar {y.scalar}")
        contractions = []
        for n, fac in self.factors.items():
            reps = [conjugate_rep(fac, self._full_reps(fb.name)[n]), self._full_reps(sname)[n],
                    self._full_reps(ff.name)[n]]
            if conj_scalar:
                reps[1] = conjugate_rep(fac, reps[1])
            if fac.abelian:
                if sum(reps, Fraction(0)) != 0:
                    raise ValueError(f"Yukawa {y.coupling}: {n} charges do not add up to zero")
                continue
            slots = [(k, r) for k, r in enumerate(reps) if r != "1"]
            if not slots:
                continue
            if len(slots) != 2:
                raise NotImplementedError(f"Yukawa {y.coupling}: {n} contraction of {len(slots)} indices")
            (k1, r1), (k2, r2) = slots
            if {r1, r2} == {"F", "Fbar"} or (r1, r2) == ("A", "A"):
                kind = "delta"
            elif r1 == r2 and r1 in ("F", "Fbar") and fac.n == 2:
                kind = "epsf"
            else:
                raise ValueError(f"Yukawa {y.coupling}: ({r1}, {r2}) of {fac.group} has no invariant")
            contractions.append((n, kind, [(k1, self._group_slot(n, (fb.name, sname, ff.name)[k1])),
                                           (k2, self._group_slot(n, (fb.name, sname, ff.name)[k2]))]))
        return proj, contractions

    def _group_slot(self, factor: str, multiplet: str) -> int:
        """Position of the slot of ``factor`` in the slots of a multiplet (the same for its conjugate)."""
        return next(i for i, s in enumerate(self._matter_slots(multiplet)) if s.factor == factor)

    # ------------------------------------------------------------------ derived fields
    def _matter_slots(self, multiplet: str, conjugate: bool = False) -> tuple[Slot, ...]:
        reps = self._full_reps(multiplet)
        slots: list[Slot] = []
        if multiplet in self.fermions:
            slots.append(Slot(tn.SPINOR))
        for n, fac in self.factors.items():
            r = reps[n]
            if fac.abelian or r == "1":
                continue
            slots.append(Slot(fac.rep_space(r), n, conjugate_rep(fac, r) if conjugate else r))
        m = self.fermions.get(multiplet) or self.scalars.get(multiplet)
        if m.flavoured and self.generations > 1:
            slots.append(Slot(tn.FLAVOUR(self.generations)))
        return tuple(slots)

    def _charges(self, multiplet: str, sign: int) -> dict[str, Fraction]:
        reps = self._full_reps(multiplet)
        return {n: sign * reps[n] for n, fac in self.factors.items() if fac.abelian and reps[n] != 0}

    def _derive_fields(self) -> list[FieldInfo]:
        infos: list[FieldInfo] = []
        half = Fraction(1, 2)
        for n, fac in self.factors.items():
            adj_slot = () if fac.abelian else (Slot(fac.adjoint_space(), n, "A"),)
            infos += [
                FieldInfo(n, "gauge", False, 0, Fraction(1), (Slot(tn.LORENTZ),) + adj_slot,
                          conjugate=n, conj_sign=1, factor=n, source=f"src_{n}"),
                FieldInfo(f"gh_{n}", "ghost", True, 1, Fraction(1), adj_slot,
                          conjugate=f"gh_{n}", conj_sign=1, factor=n, source=f"src_gh_{n}"),
                FieldInfo(f"agh_{n}", "antighost", True, -1, Fraction(1), adj_slot,
                          conjugate=f"agh_{n}", conj_sign=-1, factor=n),
                FieldInfo(f"b_{n}", "auxiliary", False, 0, Fraction(2), adj_slot,
                          conjugate=f"b_{n}", conj_sign=1, factor=n),
            ]
        for name, fm in self.fermions.items():
            infos += [
                FieldInfo(name, "fermion", True, 0, 3 * half, self._matter_slots(name), self._charges(name, 1),
                          conjugate=f"{name}bar", multiplet=name, chirality=fm.chirality, source=f"src_{name}"),
                FieldInfo(f"{name}bar", "fermion_bar", True, 0, 3 * half, self._matter_slots(name, True),
                          self._charges(name, -1), conjugate=name, multiplet=name, chirality=fm.chirality,
                          source=f"src_{name}bar"),
            ]
        for name, sm in self.scalars.items():
            if sm.complex:
                infos += [
                    FieldInfo(name, "scalar", False, 0, Fraction(1), self._matter_slots(name), self._charges(name, 1),
                              conjugate=f"{name}dag", multiplet=name, source=f"src_{name}"),
                    FieldInfo(f"{name}dag", "scalar_conj", False, 0, Fraction(1), self._matter_slots(name, True),
                              self._charges(name, -1), conjugate=name, multiplet=name, source=f"src_{name}dag"),
                ]
            else:
                if any(r not in ("1", "A", Fraction(0)) for r in self._full_reps(name).values()):
                    raise ValueError(f"real scalar {name} must be in a real representation")
                infos.append(FieldInfo(name, "scalar", False, 0, Fraction(1), self._matter_slots(name), {},
                                       conjugate=name, multiplet=name, source=f"src_{name}"))
        by_name = {fi.name: fi for fi in infos}
        sources = []
        for fi in infos:
            if fi.source is None:
                continue
            partner = by_name[fi.conjugate]
            real = fi.conjugate == fi.name
            slots = tuple(Slot(s.space, s.factor, _conj_rep_str(self.factors, s)) for s in fi.slots)
            sources.append(FieldInfo(
                fi.source, "source", not fi.odd, -1 - fi.ghost_number, 3 - fi.dimension, slots,
                {k: -v for k, v in fi.charges.items()},
                conjugate=fi.source if real else partner.source, conj_sign=-1 if real else 1,
                factor=fi.factor, multiplet=fi.multiplet, chirality=fi.chirality, source_of=fi.name))
        return infos + sources


def _is_group_factor(fa: Expression) -> bool:
    """A head (or power of a head) all of whose labels are adjoint or fundamental group labels."""
    base = next(iter(fa)) if fa.get_type() == AtomType.Pow else fa
    if base.get_type() != AtomType.Fn:
        return False
    labs = tn.labels_in(base)
    return bool(labs) and all(tn.label_space(lab).split(":")[0] in ("adjoint", "fund") for lab in labs)


def _num_expr(q: Fraction | int) -> Expression:
    q = Fraction(q)
    return Expression.num(q.numerator) / q.denominator


def _conj_rep_str(factors: Mapping[str, GaugeFactor], s: Slot) -> str | None:
    if s.factor is None:
        return None
    return conjugate_rep(factors[s.factor], s.rep)


# ---------------------------------------------------------------------------------------------- loading
def model_from_dict(doc: Mapping[str, Any]) -> Model:
    """Build a :class:`Model` from a parsed specification document (the TOML structure of the module
    docstring)."""
    known = {"name", "generations", "gauge", "fermions", "scalars", "parameters", "potential", "yukawa",
             "gauge_fixing", "regularisation"}
    unknown = set(doc) - known
    if unknown:
        raise ValueError(f"unknown specification keys: {sorted(unknown)}")
    factors = {n: GaugeFactor(n, v["group"], v["coupling"]) for n, v in doc["gauge"].items()}

    def reps_of(entry: Mapping[str, Any], extra: set[str]) -> dict[str, Rep]:
        out: dict[str, Rep] = {}
        for k, v in entry.items():
            if k in extra:
                continue
            if k not in factors:
                raise ValueError(f"unknown key or gauge factor {k!r}")
            out[k] = Fraction(str(v)) if factors[k].abelian else str(v)
        return out

    fermions = {n: FermionMultiplet(n, v["chirality"], reps_of(v, {"chirality", "flavoured"}), v.get("flavoured", True))
                for n, v in doc.get("fermions", {}).items()}
    scalars = {n: ScalarMultiplet(n, v.get("complex", True), reps_of(v, {"complex", "flavoured"}),
                                  v.get("flavoured", False))
               for n, v in doc.get("scalars", {}).items()}
    params = {n: Parameter(n, v.get("real", True), v.get("flavour", False)) for n, v in doc.get("parameters", {}).items()}
    gf = doc.get("gauge_fixing", {"kind": "linear_covariant", "xi": "xi"})
    if gf.get("kind") != "linear_covariant":
        raise NotImplementedError("only the linear covariant gauge is implemented")
    xi_raw = gf.get("xi", "xi")
    if isinstance(xi_raw, str) and not _is_number(xi_raw):
        params.setdefault(xi_raw, Parameter(xi_raw))
    potential = tuple(PotentialTerm(t["coupling"], tuple(tuple(b) for b in t["bilinears"]),
                                    Fraction(str(t.get("factor", 1)))) for t in doc.get("potential", []))
    yukawas = tuple(YukawaTerm(t["coupling"], t["fermion_bar"], t["scalar"], t["fermion"]) for t in doc.get("yukawa", []))
    reg = Regularisation(**doc.get("regularisation", {}))
    xi = _num_expr(Fraction(str(xi_raw))) if _is_number(str(xi_raw)) else S(_PARAM_NS + str(xi_raw), is_real=True)
    return Model(doc["name"], factors, fermions, scalars, params, potential, yukawas,
                 GaugeFixing("linear_covariant", xi), reg, int(doc.get("generations", 1)))


def _is_number(s: str) -> bool:
    try:
        Fraction(s)
        return True
    except ValueError:
        return False


def load_model(path: str | Path) -> Model:
    """Load a specification file (TOML)."""
    with open(path, "rb") as fh:
        return model_from_dict(tomllib.load(fh))


def builtin_model(name: str) -> Model:
    """A specification shipped with the package (``hepkit_renorm/models/<name>.toml``), e.g. ``"sm_one_family"``."""
    return load_model(Path(__file__).parent / "models" / f"{name}.toml")
