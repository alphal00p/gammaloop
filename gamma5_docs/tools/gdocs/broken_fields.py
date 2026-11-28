"""Field basis of the broken phase and the substitution ``H -> H_v + phi`` on ``hepkit_renorm`` functionals.

This module is the engine behind :mod:`gdocs.broken`. It defines

* :data:`BT`, a ``hepkit_renorm`` :class:`~hepkit_renorm.fields.FieldTable` of the broken phase in the field basis
  of report Section 5.4.2 (``sec:bl-fields``): the gluon ``G``; the charged and neutral electroweak gauge fields
  ``Wp``, ``Wm``, ``Z``, ``A`` (photon) of eq:bl-rot, and the same rotation for the ghosts, antighosts, auxiliary
  fields and the BRST sources rho and zeta; the scalars ``h``, ``G0``, ``Gp``, ``Gm`` of eq:bp-H; the mass
  eigenstates ``u``, ``d``, ``e``, ``nu`` and the sterile partners ``up``, ``dp``, ``ep`` of eq:bl-fermions (Dirac
  fields, with their conjugates ``...bar``). The report does not rotate the fermion sources ``R_f`` and the Higgs
  source ``sigma``: they are kept in their SU(2) components (``src_Q_1`` is the upper component of ``Rbar_Q``,
  ``src_H_2`` the lower component of ``sigma^dagger``).
* :func:`shift`, the exact substitution of a ``LocalFunctional`` of the unbroken model (package field table) into a
  ``LocalFunctional`` over :data:`BT`: every SU(2) index is written in explicit components (``T^a = tau^a/2``,
  ``f^{abc} = eps^{abc}``, ``eps_{12} = 1``, the conventions of the package), and the components are replaced by

  ``W^1 = (W^+ + W^-)/sqrt2``, ``W^2 = i (W^+ - W^-)/sqrt2``, ``W^3 = c_W Z + s_W A``, ``B = -s_W Z + c_W A``
  (the same for ``gh``, ``agh``, ``b``, ``src`` (rho), ``src_gh`` (zeta)); ``H^1 = G^+``,
  ``H^2 = (v + h + i G^0)/sqrt2`` and the conjugates; ``Psi^1_Q = P_L u + P_R u'``, ``Psi_U = P_R u + P_L u'``,
  ``Psi^2_Q = P_L d + P_R d'``, ``Psi_D = P_R d + P_L d'``, ``Psi^1_L = nu``, ``Psi^2_L = P_L e + P_R e'``,
  ``Psi_E = P_R e + P_L e'``, and ``Psibar = Psi^dagger gamma^0`` accordingly (``ubar P_R + upbar P_L`` for
  ``Psibar^1_Q``).

  ``v``, ``c_W``, ``s_W`` are the site-notation symbols of :mod:`gdocs.mathview`; ``sqrt2`` is the symbol
  :data:`RT2` (``RT2^2`` is reduced to 2). The relations ``c_W^2 + s_W^2 = 1`` and ``g s_W = g' c_W`` are not used
  by the substitution; :func:`reduce_relations` and :func:`canonical_value` apply them to scalar coefficients.
* :class:`BrokenModel`, the minimal model interface that ``hepkit_renorm.fingerprint.Fingerprinter`` needs, so that
  broken-phase functionals can be fingerprinted (exact plane-wave evaluation over GF(p)) like the package's own.

Nothing here is typed from the report's broken-phase formulas: the report's statements are reproduced by
``tools/tests/test_broken.py`` from these definitions.
"""

from __future__ import annotations

import itertools
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction

from symbolica import AtomType, Expression

from hepkit_renorm import Density, LocalFunctional, Monomial, Occurrence, builtin_model
from hepkit_renorm import tensors as tn
from hepkit_renorm.fields import FieldInfo, FieldTable, Slot
from hepkit_renorm.spec import Parameter

from . import mathview as mv

__all__ = ["BT", "BrokenModel", "RT2", "V", "CW", "SW", "shift", "unbroken_model", "reduce_relations",
           "canonical_value", "BROKEN_PARAMS", "simplifier", "display_simplifier", "site_field"]

# ---------------------------------------------------------------------------------------------- parameters
V = mv.parse("v")
CW = mv.parse("c_W")
SW = mv.parse("s_W")
RT2 = mv.sym("sqrt2")
"""The symbol sqrt(2) (a plain site-notation symbol; :func:`for_display` writes it as ``sqrtn(2)``)."""
SQRT2_DISPLAY = mv.parse("sqrtn(2)")


def for_display(e: Expression) -> Expression:
    """A site expression with :data:`RT2` written as ``2/sqrtn(2)`` (printed √2, in the denominator where the
    coefficient allows: ``y_u v/√2``)."""
    return e.replace(RT2, 2 * SQRT2_DISPLAY ** -1)

_MODEL = None


def unbroken_model():
    global _MODEL
    if _MODEL is None:
        _MODEL = builtin_model("sm_one_family")
    return _MODEL


# ---------------------------------------------------------------------------------------------- field table
_LOR = tn.LORENTZ
_SPIN = tn.SPINOR
_ADJ_G = tn.adjoint_space("G", 8)
_FUND_G = tn.fundamental_space("G", 3)


def _fi(name, kind, odd, gh, dim, slots, conj, conj_sign=1, factor=None, multiplet=None, chirality=None,
        source_of=None, source=None):
    return FieldInfo(name, kind, odd, gh, Fraction(dim), tuple(slots), {}, conj, conj_sign, factor, multiplet,
                     chirality, source_of, source)


def _build_table() -> FieldTable:
    L = Slot(_LOR)
    AG = Slot(_ADJ_G, "G", "A")
    SP = Slot(_SPIN)
    CF, CFb = Slot(_FUND_G, "G", "F"), Slot(_FUND_G, "G", "Fbar")
    infos = []
    # gauge sector: gluon, then the electroweak components (+, -, Z, photon)
    infos += [_fi("G", "gauge", False, 0, 1, [L, AG], "G", factor="G", source="src_G"),
              _fi("gh_G", "ghost", True, 1, 1, [AG], "gh_G", factor="G", source="src_gh_G"),
              _fi("agh_G", "antighost", True, -1, 1, [AG], "agh_G", -1, factor="G"),
              _fi("b_G", "auxiliary", False, 0, 2, [AG], "b_G", factor="G")]
    pairs = {"Wp": "Wm", "Wm": "Wp", "Z": "Z", "A": "A"}
    src_of = {"Wp": "Wm", "Wm": "Wp", "Z": "Z", "A": "A"}   # rho_+ multiplies (D omega)^-: it is the source of W^-
    for x in ("Wp", "Wm", "Z", "A"):
        c = pairs[x]
        infos += [_fi(x, "gauge", False, 0, 1, [L], c, factor=x, source=f"src_{src_of[x]}"),
                  _fi(f"gh_{x}", "ghost", True, 1, 1, [], f"gh_{c}", factor=x, source=f"src_gh_{src_of[x]}"),
                  _fi(f"agh_{x}", "antighost", True, -1, 1, [], f"agh_{c}", -1, factor=x),
                  _fi(f"b_{x}", "auxiliary", False, 0, 2, [], f"b_{c}", factor=x)]
    # fermions: mass eigenstates and sterile partners
    for x, coloured in (("u", True), ("up", True), ("d", True), ("dp", True), ("e", False), ("ep", False),
                        ("nu", False)):
        sl = [SP, CF] if coloured else [SP]
        slb = [SP, CFb] if coloured else [SP]
        infos += [_fi(x, "fermion", True, 0, Fraction(3, 2), sl, x + "bar", multiplet=None),
                  _fi(x + "bar", "fermion_bar", True, 0, Fraction(3, 2), slb, x, multiplet=None)]
    # scalars
    infos += [_fi("h", "scalar", False, 0, 1, [], "h"), _fi("G0", "scalar", False, 0, 1, [], "G0"),
              _fi("Gp", "scalar", False, 0, 1, [], "Gm"), _fi("Gm", "scalar_conj", False, 0, 1, [], "Gp")]
    # sources: rho and zeta (rotated), fermion sources and sigma (SU(2) components)
    infos += [_fi("src_G", "source", True, -1, 2, [L, AG], "src_G", -1, factor="G", source_of="G"),
              _fi("src_gh_G", "source", False, -2, 2, [AG], "src_gh_G", -1, factor="G", source_of="gh_G")]
    for x in ("Wp", "Wm", "Z", "A"):
        c = pairs[x]
        infos += [_fi(f"src_{x}", "source", True, -1, 2, [L], f"src_{c}", -1, factor=x, source_of=src_of[x]),
                  _fi(f"src_gh_{x}", "source", False, -2, 2, [], f"src_gh_{c}", -1, factor=x,
                      source_of=f"gh_{src_of[x]}")]
    # fermion sources: src_<f>_<i> = component i of Rbar_f (row object); src_<f>bar_<i> = component of R_f
    comp_field = {("Q", 1): "u", ("Q", 2): "d", ("U", 0): "u", ("D", 0): "d", ("L", 1): "nu", ("L", 2): "e",
                  ("E", 0): "e"}
    for (f, i), phys in comp_field.items():
        nm = f"src_{f}" + (f"_{i}" if i else "")
        nb = f"src_{f}bar" + (f"_{i}" if i else "")
        coloured = f in ("Q", "U", "D")
        infos += [_fi(nm, "source", False, -1, Fraction(3, 2), [SP, CFb] if coloured else [SP], nb, multiplet=f,
                      source_of=phys),
                  _fi(nb, "source", False, -1, Fraction(3, 2), [SP, CF] if coloured else [SP], nm, multiplet=f,
                      source_of=phys + "bar")]
    for i, phys in ((1, "Gp"), (2, "h")):
        infos += [_fi(f"src_H_{i}", "source", True, -1, 2, [], f"src_Hdag_{i}", multiplet="H", source_of=phys),
                  _fi(f"src_Hdag_{i}", "source", True, -1, 2, [], f"src_H_{i}", multiplet="H",
                      source_of="Gm" if i == 1 else "h")]
    return _SiteTable(infos)


class _SiteTable(FieldTable):
    """The broken-phase field table; :meth:`site_field` gives the site-notation symbol of a field (used by
    :func:`gdocs.mathview.field`)."""

    def site_field(self, name: str, idx: Sequence[Expression]) -> Expression | None:
        return site_field(name, idx)


_PLUS = mv.symbol("+")
_MINUS = mv.symbol("−")
_GAM = mv.symbol("γ")
_ZI = mv.sym("Z")
_COMP = {"Wp": _PLUS, "Wm": _MINUS, "Z": _ZI, "A": _GAM}
_FERM = {"u": "uf", "up": "ufprime", "d": "dqf", "dp": "dqfprime", "e": "elf", "ep": "elfprime", "nu": "nuf"}
_FERMB = {"u": "ufbar", "up": "ufprimebar", "d": "dqfbar", "dp": "dqfprimebar", "e": "elfbar", "ep": "elfprimebar",
          "nu": "nufbar"}


def site_field(name: str, idx: Sequence[Expression] = ()) -> Expression | None:
    """Site-notation symbol of a broken-phase field with the given (placed) indices; ``None`` if ``name`` is not
    a broken-phase field."""
    idx = list(idx)
    P = mv.parse
    if name in ("Wp", "Wm", "Z", "A"):
        head = {"Wp": "Wp", "Wm": "Wm", "Z": "Z", "A": "Agam"}[name]
        return P(head)(*idx) if idx else P(head)
    if name == "G":
        return P("G")(*idx)
    for pre, head in (("gh_", "omega"), ("agh_", "omegabar"), ("b_", "b"), ("src_gh_", "zsrc")):
        if name.startswith(pre) and name[len(pre):] in _COMP:
            return P(head)(_COMP[name[len(pre):]])
        if name == pre + "G":
            return P(head)(*(idx or [P("A")]))
    if name.startswith("src_") and name[4:] in _COMP:
        return P("rho")(*(idx + [_COMP[name[4:]]]))
    if name == "src_G":
        return P("rho")(*idx)
    if name in ("h", "G0", "Gp", "Gm"):
        return P(name)
    base = name[:-3] if name.endswith("bar") else name
    if base in _FERM:
        head = _FERMB[base] if name.endswith("bar") else _FERM[base]
        return P(head)(*idx) if idx else P(head)
    if name.startswith("src_H"):
        i = name[-1]
        head = "sigmadag" if name.startswith("src_H_") else "sigma"
        return P(head)(P(i))
    if name.startswith("src_") and name[4] in "QUDLE":
        rest = name[4:]
        f = rest[0]
        bar = "bar" in rest
        comp = rest.split("_")[1] if "_" in rest else None
        args = [P(f)] + ([mv.parse(f"up({comp})")] if comp else []) + idx
        return P("R" if bar else "Rbar")(*args)
    return None


BT = _build_table()

# ---------------------------------------------------------------------------------------------- model for fingerprints
BROKEN_PARAMS = ("v", "c_W", "s_W", "sqrt2")


class BrokenModel:
    """The attributes of :class:`hepkit_renorm.spec.Model` that :class:`~hepkit_renorm.fingerprint.Fingerprinter`
    uses, for the broken-phase field table: one gauge factor with group slots (colour), the parameters of the
    unbroken model plus ``g``, ``gp`` (couplings of factors without slots here) and ``v``, ``c_W``, ``s_W``,
    ``sqrt2``."""

    def __init__(self) -> None:
        um = unbroken_model()
        self.fields = BT
        self.generations = 1
        self.factors = {"G": um.factors["G"]}
        self.parameters = {"g": Parameter("g"), "gp": Parameter("gp"), **dict(um.parameters)}
        for n in BROKEN_PARAMS:
            self.parameters[n] = Parameter(n)
        self._um = um
        self.name = "sm_one_family_broken"

    def param(self, name: str) -> Expression:
        if name in ("v", "c_W", "s_W"):
            return mv.parse(name)
        if name == "sqrt2":
            return RT2
        return self._um.param(name)

    def coupling(self, factor: str) -> Expression:
        return self._um.coupling(factor)

    def simplifier(self, expand_projectors: bool = True) -> tn.Simplifier:
        return tn.Simplifier({"G": (3, 8)}, 1, expand_projectors)


def simplifier() -> tn.Simplifier:
    """Canonical-mode simplifier of the broken table (projectors expanded)."""
    return _SIMP


def display_simplifier() -> tn.Simplifier:
    """Readable-mode simplifier (projectors kept)."""
    return _DSIMP


_SIMP = tn.Simplifier({"G": (3, 8)}, 1, True)
_DSIMP = tn.Simplifier({"G": (3, 8)}, 1, False)


# ---------------------------------------------------------------------------------------------- SU(2) components
_I = tn.I
_HALF = Expression.num(1) / 2


def _n(x) -> Expression:
    x = Fraction(x)
    return Expression.num(x.numerator) / x.denominator


# tau^a / 2, a = 1, 2, 3 (Pauli order), indices i, j = 1, 2
_TAU = {
    1: {(1, 2): _HALF, (2, 1): _HALF},
    2: {(1, 2): -_I * _HALF, (2, 1): _I * _HALF},
    3: {(1, 1): _HALF, (2, 2): -_HALF},
}


def _levi3(a, b, c) -> int:
    if len({a, b, c}) < 3:
        return 0
    perm = (a, b, c)
    inv = sum(1 for x, y in itertools.combinations(perm, 2) if x > y)
    return -1 if inv % 2 else 1


def _su2_value(name: str, vals: Sequence[int]) -> Expression:
    if name == "T":
        return _TAU[vals[0]].get((vals[1], vals[2]), Expression.num(0))
    if name == "f":
        return Expression.num(_levi3(*vals))
    if name == "delta":
        return Expression.num(1 if vals[0] == vals[1] else 0)
    if name == "epsf":
        return Expression.num({(1, 2): 1, (2, 1): -1}.get(tuple(vals), 0))
    raise ValueError(f"SU(2) head {name}")


def _is_su2_label(lab: Expression) -> bool:
    sp = tn.label_space(lab)
    return sp in ("adjoint:W", "fund:W")


_SU2_HEADS = {tn.T.get_name(): "T", tn.f.get_name(): "f", tn.delta.get_name(): "delta", tn.epsf.get_name(): "epsf"}


def _split_su2(term: Expression) -> tuple[list[tuple[str, list[Expression], int]], Expression]:
    """``(SU(2) heads [(name, labels, power)], rest)`` of a product term."""
    factors = list(term) if term.get_type() == AtomType.Mul else [term]
    heads, rest = [], Expression.num(1)
    for fac in factors:
        base, power = fac, 1
        if fac.get_type() == AtomType.Pow:
            b, ex = list(fac)
            if b.get_type() == AtomType.Fn and ex.get_type() == AtomType.Num:
                base, power = b, int(ex.to_canonical_string())
        if base.get_type() == AtomType.Fn and base.get_name() in _SU2_HEADS:
            args = list(base)
            if args and any(a.get_type() == AtomType.Fn and a.get_name() in (tn.adj.get_name(), tn.fund.get_name())
                            and _is_su2_label(a) for a in args):
                heads.append((_SU2_HEADS[base.get_name()], args, power))
                continue
        rest = rest * fac
    return heads, rest


# component maps: unbroken component -> [(coefficient, broken field or None for the constant)]
def _rot_w(comp: int, names: Mapping[str, str]) -> list[tuple[Expression, str]]:
    """``X^a`` of SU(2) in the rotated basis (eq:bl-rot): ``names`` maps "+", "-", "Z", "A" to field names."""
    if comp == 1:
        return [(1 / RT2, names["+"]), (1 / RT2, names["-"])]
    if comp == 2:
        return [(_I / RT2, names["+"]), (-_I / RT2, names["-"])]
    return [(CW, names["Z"]), (SW, names["A"])]


def _rot_b(names: Mapping[str, str]) -> list[tuple[Expression, str]]:
    return [(-SW, names["Z"]), (CW, names["A"])]


def _names(prefix: str) -> dict[str, str]:
    return {"+": prefix + "Wp", "-": prefix + "Wm", "Z": prefix + "Z", "A": prefix + "A"}


# fermion component maps: (multiplet, SU(2) component or 0) -> [(projector head or None, broken field)]
_PSI = {("Q", 1): [("PL", "u"), ("PR", "up")], ("Q", 2): [("PL", "d"), ("PR", "dp")],
        ("U", 0): [("PR", "u"), ("PL", "up")], ("D", 0): [("PR", "d"), ("PL", "dp")],
        ("L", 1): [(None, "nu")], ("L", 2): [("PL", "e"), ("PR", "ep")], ("E", 0): [("PR", "e"), ("PL", "ep")]}
_OTHER = {"PL": "PR", "PR": "PL", None: None}


@dataclass(frozen=True)
class _Opt:
    """One term of the image of an occurrence: ``factor * field[slots] d^derivatives`` (``field`` None: the
    constant part of ``H``, ``factor`` then includes ``v/sqrt2``)."""

    factor: Expression
    field: str | None
    slots: tuple[Expression, ...] = ()
    ders: tuple[Expression, ...] = ()


def _occurrence_images(o: Occurrence, info: FieldInfo, comp: dict[str, int]) -> list[_Opt]:
    name = o.field
    ders = tuple(o.derivatives)
    labs = list(o.slots)
    # SU(2) slot value and the remaining slots
    su2 = None
    keep = []
    for sl, lab in zip(info.slots, labs):
        if sl.space.kind in ("adjoint", "fund") and sl.space.factor == "W":
            su2 = comp[lab.to_canonical_string()]
        else:
            keep.append(lab)
    keep = tuple(keep)
    one = Expression.num(1)
    kind = info.kind
    # gauge sector of G: unchanged
    if name in ("G", "gh_G", "agh_G", "b_G", "src_G", "src_gh_G"):
        return [_Opt(one, name, keep, ders)]
    for pre in ("", "gh_", "agh_", "b_", "src_", "src_gh_"):
        if name == pre + "W":
            return [_Opt(c, f, keep, ders) for c, f in _rot_w(su2, _names(pre))]
        if name == pre + "B":
            return [_Opt(c, f, keep, ders) for c, f in _rot_b(_names(pre))]
    if name in ("H", "Hdag"):
        conj = name == "Hdag"
        if su2 == 1:
            return [_Opt(one, "Gm" if conj else "Gp", (), ders)]
        out = [_Opt(1 / RT2, "h", (), ders), _Opt((-_I if conj else _I) / RT2, "G0", (), ders)]
        if not ders:
            out.insert(0, _Opt(V / RT2, None))
        return out
    if name in ("src_H", "src_Hdag"):
        return [_Opt(one, f"{name}_{su2}", (), ders)]
    if kind in ("fermion", "fermion_bar"):
        mult = info.multiplet
        key = (mult, su2 or 0)
        spin, rest = keep[0], keep[1:]
        out = []
        for proj, f in _PSI[key]:
            if proj is None:
                out.append(_Opt(one, f if kind == "fermion" else f + "bar", (spin,) + rest, ders))
                continue
            t = tn.fresh(tn.SPINOR)
            if kind == "fermion":       # Psi_s = P_{s t} psi_t
                out.append(_Opt(tn.HEADS[proj](spin, t), f, (t,) + rest, ders))
            else:                       # Psibar_s = psibar_t (P')_{t s}, P' = gamma0 P^dagger gamma0 = the other one
                out.append(_Opt(tn.HEADS[_OTHER[proj]](t, spin), f + "bar", (t,) + rest, ders))
        return out
    if kind == "source" and info.multiplet in ("Q", "U", "D", "L", "E"):
        bar = name.endswith("bar")
        base = f"src_{info.multiplet}" + ("bar" if bar else "")
        nm = base + (f"_{su2}" if su2 else "")
        return [_Opt(one, nm, keep, ders)]
    raise ValueError(f"no broken-phase image for the field {name}")


def reduce_sqrt2(e: Expression) -> Expression:
    """``sqrt2^n`` reduced with ``sqrt2^2 = 2`` (expanded)."""
    e = e.expand()
    for _ in range(8):
        new = e.replace(RT2 ** 2, Expression.num(2)).replace(RT2 ** -2, _HALF).expand()
        # higher powers
        new = new.replace(RT2 ** 3, 2 * RT2).replace(RT2 ** 4, Expression.num(4)).replace(RT2 ** -1, RT2 / 2)
        new = new.replace(RT2 ** -3, RT2 / 4).replace(RT2 ** -4, Expression.num(1) / 4).expand()
        if new == e:
            break
        e = new
    return e


def shift(F: LocalFunctional | Density, *, constant: bool = True) -> LocalFunctional:
    """The image of an unbroken-phase functional in the broken-phase basis (module docstring). ``constant=False``
    drops the constant part of ``H`` (``v -> 0``)."""
    src = F.fields
    total = []
    for m in F.monomials():
        coeff = m.coefficient.expand()
        terms = list(coeff) if coeff.get_type() == AtomType.Add else [coeff]
        # SU(2) labels of the monomial (occurrences and coefficient)
        labs: dict[str, Expression] = {}
        for o in m.occurrences:
            for sl, lab in zip(src[o.field].slots, o.slots):
                if sl.space.kind in ("adjoint", "fund") and sl.space.factor == "W":
                    labs[lab.to_canonical_string()] = lab
        split = []
        for t in terms:
            heads, rest = _split_su2(t)
            for _, args, _ in heads:
                for a in args:
                    if _is_su2_label(a):
                        labs.setdefault(a.to_canonical_string(), a)
            split.append((heads, rest))
        keys = list(labs)
        ranges = [(1, 2, 3) if tn.label_space(labs[k]) == "adjoint:W" else (1, 2) for k in keys]
        for vals in itertools.product(*ranges):
            comp = dict(zip(keys, vals))
            c = Expression.num(0)
            for heads, rest in split:
                num = Expression.num(1)
                for name, args, power in heads:
                    hv = _su2_value(name, [comp[a.to_canonical_string()] for a in args])
                    num = num * hv ** power
                    if num == Expression.num(0):
                        break
                if num != Expression.num(0):
                    c = c + num * rest
            c = c.expand()
            if c == Expression.num(0):
                continue
            opts = [_occurrence_images(o, src[o.field], comp) for o in m.occurrences]
            for choice in itertools.product(*opts):
                if not constant and any(op.field is None for op in choice):
                    continue
                cc = c
                occs = []
                for op in choice:
                    cc = cc * op.factor
                    if op.field is not None:
                        occs.append(Occurrence(op.field, op.slots, op.ders))
                total.append(Monomial(cc, tuple(occs)))
    dens = Density.from_monomials(BT, total)
    dens = dens.map_coefficients(reduce_sqrt2)
    return LocalFunctional.integrate(dens)


# ---------------------------------------------------------------------------------------------- scalar relations
def _p(name: str) -> Expression:
    return unbroken_model().param(name)


def _reduce_powers(e: Expression, x: Expression, sq: Expression) -> Expression:
    """Replace ``x^n`` (n >= 2) by ``x^(n mod 2) sq^(n div 2)`` (``sq`` = the value of ``x^2``), expanded."""
    e = e.expand()
    for _ in range(12):
        new = e
        for n in range(12, 1, -1):
            new = new.replace(x ** n, (x ** (n % 2)) * sq ** (n // 2))
        new = new.expand()
        if new == e:
            return e
        e = new
    return e


def vev(e: Expression) -> Expression:
    """The vacuum condition of report eq:bl-vev, ``m^2 = -lambda v^2`` (package symbols)."""
    return e.replace(_p("m2"), -_p("lambda") * V ** 2)


def canonical_value(e: Expression) -> Expression:
    """A canonical form of a scalar coefficient modulo ``sqrt2^2 = 2``, ``c_W^2 + s_W^2 = 1`` and
    ``g s_W = g' c_W`` (report eq:bp-rot): ``g'`` is eliminated, ``s_W`` occurs at most linearly in the numerator,
    the denominator is a monomial in ``c_W``. Two coefficients are equal iff their canonical forms are equal."""
    g, gp = _p("g"), _p("gp")
    x = e.replace(mv.parse("e"), g * SW).replace(gp, g * SW / CW)
    x = _reduce_powers(x, RT2, Expression.num(2))
    x = x.together()
    num, den = _num_den(x)
    num = _reduce_powers(_reduce_powers(num, SW, 1 - CW ** 2), RT2, Expression.num(2))
    den = _reduce_powers(_reduce_powers(den, SW, 1 - CW ** 2), RT2, Expression.num(2))
    if num == Expression.num(0):
        return Expression.num(0)
    return (num / den).together()


def _num_den(x: Expression) -> tuple[Expression, Expression]:
    """Numerator and denominator of a rational expression written over a common denominator."""
    x = x.together()
    if x.get_type() == AtomType.Mul:
        num, den = Expression.num(1), Expression.num(1)
        for f in x:
            if f.get_type() == AtomType.Pow:
                b, ex = list(f)
                if ex.get_type() == AtomType.Num and ex.to_canonical_string().startswith("-"):
                    den = den * b ** (-ex)
                    continue
            num = num * f
        return num.expand(), den.expand()
    if x.get_type() == AtomType.Pow:
        b, ex = list(x)
        if ex.get_type() == AtomType.Num and ex.to_canonical_string().startswith("-"):
            return Expression.num(1), (b ** (-ex)).expand()
    return x.expand(), Expression.num(1)


def is_zero(e: Expression) -> bool:
    return canonical_value(e) == Expression.num(0)


def equal(a: Expression, b: Expression) -> bool:
    return is_zero(a - b)


def to_site(e: Expression) -> Expression:
    """A scalar coefficient in the site notation, not expanded (package parameters -> site symbols); ``sqrt2``
    stays (:func:`for_display` writes it as a root)."""
    for n in _SITE_PARAMS:
        e = e.replace(_p(n), mv.parse(n))
    return e


_SITE_PARAMS = ("g_s", "g", "gp", "y_u", "y_d", "y_e", "lambda", "m2", "xi")


def _length(e: Expression) -> int:
    """Display length: characters of the plain form in site symbols (``e`` counts as long as ``g s_W``, so that
    the electric charge is used only where it shortens the expression)."""
    t = to_site(e).format_plain().replace(" ", "")
    return len(t) + 2 * len(re.findall(r"(?<![A-Za-z_])e(?![A-Za-z_])", t)) + 3 * t.count("gp")


def _forms(e: Expression) -> list[Expression]:
    """Equivalent forms of a coefficient: as given, expanded, factored; with ``g'`` eliminated and ``s_W^2`` or
    ``c_W^2`` reduced; with ``g = e/s_W``, ``g' = e/c_W``."""
    g, gp = _p("g"), _p("gp")
    E_ = mv.parse("e")
    e = e.replace(E_, g * SW)
    base = _reduce_powers(e, RT2, Expression.num(2)) if e.get_type() != AtomType.Mul else e
    out = [e, base]
    exp_ = _reduce_powers(e.expand(), RT2, Expression.num(2))
    out.append(exp_)
    variants = [exp_]
    a = exp_.replace(gp, g * SW / CW)
    variants += [_reduce_powers(a, SW, 1 - CW ** 2), _reduce_powers(a, CW, 1 - SW ** 2)]
    c_ = exp_.replace(g, gp * CW / SW)
    variants += [_reduce_powers(c_, SW, 1 - CW ** 2), _reduce_powers(c_, CW, 1 - SW ** 2)]
    b = exp_.replace(g, E_ / SW).replace(gp, E_ / CW)
    variants += [_reduce_powers(b, SW, 1 - CW ** 2), _reduce_powers(b, CW, 1 - SW ** 2), b.expand()]
    for v in variants:
        v = _reduce_powers(v.expand(), RT2, Expression.num(2))
        out.append(v)
        try:
            t = v.together()
            out.append(t)
            out.append(_reduce_powers(t, RT2, Expression.num(2)))
        except Exception:  # noqa: BLE001 - together() may refuse function arguments; the form is optional
            pass
        try:
            out.append(v.factor())
        except Exception:  # noqa: BLE001
            pass
    return out


def simplest(e: Expression) -> Expression:
    """The shortest of :func:`_forms` that equals ``e`` (checked with :func:`canonical_value`)."""
    ref = canonical_value(e)
    best, blen = e, _length(e)
    for f in _forms(e):
        n = _length(f)
        if n < blen and canonical_value(f) == ref:
            best, blen = f, n
    return best
