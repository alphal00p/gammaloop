"""Operator catalogue of the broken phase: ``site/data/broken.json`` (DATA_FORMAT.md, section 4).

The broken phase is the renormalised action of the unbroken phase in the variable ``phi = H - H_v`` (report
Section 5.4, eq:bl-def), in the same linear gauge at ``xi = 1``:

    S_ren^broken[phi, Phi'] = S_ren[H = H_v + phi, Phi'],   H_v = (0, v/sqrt2),   m^2 = -lambda v^2 (eq:bl-vev).

Nothing is recomputed and nothing is typed: every term is derived from the unbroken catalogue of
:mod:`gdocs.catalog` (its operators, ``catalog.operator_specs``, and its coefficients, ``data.values`` of
``site/data/unbroken.json``).

Method
======

1. The shift. For every unbroken term ``P`` (operator ``O_P``, four-dimensional form ``O4_P`` for the finite part,
   coefficients ``c_P[part]``) the operator is substituted by :func:`gdocs.broken_fields.shift`: every SU(2)
   index in explicit components, ``H = H_v + phi`` (``k`` of the scalar legs of each monomial become the
   constant ``H_v``, for every ``k``), and the field basis of the report (eq:bp-H, eq:bl-rot, eq:bl-fermions):
   ``W^+-``, ``Z``, ``A^gamma`` and the same rotation for ghosts, antighosts, ``b``, ``rho``, ``zeta``; ``h``,
   ``G0``, ``G^+-``; the mass eigenstates ``u, d, e, nu`` and their sterile partners ``u', d', e'``. The fermion
   sources ``R_f`` and the Higgs source ``sigma``, which the report does not rotate, are kept in SU(2) components.
2. Collection. The image of ``O_P`` is split by broken-phase field content. On each content the images of all
   parents are grouped by proportionality (exactly, modulo ``c_W^2 + s_W^2 = 1``, ``g s_W = g' c_W`` and
   ``sqrt2^2 = 2``; Lorentz indices split into barred and hatted parts first): each group is one broken-phase
   operator ``O``, normalised from its first parent, and each parent contributes ``lambda_P O`` with an exact factor
   ``lambda_P`` (a power ``v^k`` times rotation factors). The groups of one content are checked to be linearly
   independent with exact fingerprints in ``D = 6`` (``hepkit_renorm.fingerprint`` on the broken field table);
   where they are not (some quartic gauge contents), the content is written on its elementary monomials instead.
   Two operators of one content whose coefficients are proportional by a number in every part are summed into one
   (for example ``ubar P_R u`` and ``ubar P_L u`` of the two Yukawa monomials give ``ubar u``).
3. Coefficients: ``c[part] = sum_P lambda_P c_P[part]``, then ``m^2 = -lambda v^2``. The finite part multiplies the
   image of ``O4_P`` (checked to be ``lambda_P`` times one four-dimensional operator ``O4``). Coefficients are shown
   in the shortest equivalent form; ``pole.decomposition`` is the same sum written through the renormalisation
   constants (``catalog.decomposition`` of each parent).
4. Feynman rules: ``mathview.vertex`` of coefficient times operator, i.e. ``LocalFunctional.vertex_function`` of the
   package on the broken field table.

Checks: ``tools/tests/test_broken.py`` (exactness of the shift against the order store by an independent route,
the report's statements of Section 5.4, ``v = 0``).

Interface
=========

``build_broken(site_dir=SITE, *, log=print) -> dict`` writes ``data/broken.json`` and returns counts;
``python -m gdocs.broken [--site DIR]``. ``catalogue(ctx) -> list[BTerm]`` is the derived catalogue (used by the
test). Needs ``SYMBOLICA_LICENSE`` in the environment and ``site/data/unbroken.json`` and ``classes.json`` (run
``gdocs.catalog`` first).

Extra keys of the written file (ignored by the page): per term ``data = {"values", "parents", "v_powers"}``:
``values`` the coefficients (Symbolica strings in ``g_s g gp y_u y_d y_e lambda xi v c_W s_W sqrt2 I``; pole and
evanescent in units of ``1/(16 pi^2 eps)``, finite in ``1/(16 pi^2)``); ``parents`` the list of
``{"id", "factor", "v_power"}`` (``factor`` = ``lambda_P`` as a string); top level ``counts``.
"""

from __future__ import annotations

import argparse
import json
import re
import time
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import catalog as cat
from . import mathview as mv

from symbolica import AtomType, Expression, S  # noqa: E402

from hepkit_renorm import Density, LocalFunctional  # noqa: E402
from hepkit_renorm import tensors as tn  # noqa: E402
from hepkit_renorm import linalg  # noqa: E402
from hepkit_renorm.fingerprint import Fingerprinter  # noqa: E402

from . import broken_fields as bf  # noqa: E402
from . import names as nm  # noqa: E402

__all__ = ["BTerm", "Ctx", "build_broken", "catalogue", "fingerprinter", "broken_values"]

SITE = cat.SITE
PARTS = cat.PARTS
LOOP = {"tree": None, "pole": "pole", "evanescent": "pole", "finite": "finite", "breaking": None}
BT = bf.BT
V, CW, SW, RT2 = bf.V, bf.CW, bf.SW, bf.RT2


def _ref(label: str) -> dict:
    from . import theory

    return theory.report_ref(label)


TABLES = {"bl-Z": {"label": "Table 63 (constants in the broken phase)", "page": 160},
          "bl-ST": {"label": "Table 64 (finite counterterm in the broken phase)", "page": 166},
          "bp-twopoint": {"label": "Table 66 (two-point functions)", "page": 172}}


# ---------------------------------------------------------------------------------------------- context
class Ctx:
    """Inputs: the unbroken catalogue context, the operators and the coefficients of unbroken.json."""

    def __init__(self, site_dir: str | Path = SITE, log=print) -> None:
        self.site = Path(site_dir)
        self.log = log
        self.c = cat.Context(self.site, log=lambda s: None)
        self.model = self.c.model
        self.specs = {s.id: s for s in cat.operator_specs(self.c)}
        self.unbroken = json.loads((self.site / "data" / "unbroken.json").read_text())
        if self.unbroken.get("mock"):
            raise RuntimeError("data/unbroken.json is mock data: run gdocs.catalog first")
        self.terms = {t["id"]: t for t in self.unbroken["terms"]}
        self.order = [t["id"] for t in self.unbroken["terms"]]
        self.coeffs = {}
        for t in self.unbroken["terms"]:
            if t["id"] not in self.specs:
                raise RuntimeError(f"unbroken term {t['id']} has no operator spec")
            self.coeffs[t["id"]] = {p: cat.parse_value(self.model, v) for p, v in t["data"]["values"].items()}
        self.bmodel = bf.BrokenModel()


def fingerprinter(dimension: int = 6, seed: int = 0) -> tuple[Fingerprinter, dict]:
    """A fingerprinter of the broken table with parameter values that satisfy ``c_W^2 + s_W^2 = 1``,
    ``g s_W = g' c_W``, ``sqrt2^2 = 2`` and ``m^2 = -lambda v^2``; returns it and the values (GF(p) ints)."""
    from hepkit_renorm.exact import default_field

    field_ = default_field()
    p = field_.p
    vals = broken_values(field_, seed)
    fp = Fingerprinter(bf.BrokenModel(), field_, dimension, values=vals)
    return fp, vals


def broken_values(field_, seed: int = 0) -> dict:
    """Consistent parameter values over GF(p) (see :func:`fingerprinter`)."""
    p = field_.p
    rng = np.random.default_rng(1234 + seed)
    r = lambda: int(rng.integers(2, p - 1))  # noqa: E731
    t = r()
    inv = lambda x: pow(x % p, p - 2, p)  # noqa: E731
    c = (1 - t * t) * inv(1 + t * t) % p
    s = 2 * t * inv(1 + t * t) % p
    g = r()
    lam, v = r(), r()
    vals = {"g": g, "gp": g * s * inv(c) % p, "c_W": c, "s_W": s, "sqrt2": field_.sqrt(2), "v": v,
            "lambda": lam, "m2": (-lam * v * v) % p, "g_s": r(), "y_u": r(), "y_d": r(), "y_e": r(), "xi": r()}
    return vals


# ---------------------------------------------------------------------------------------------- structures
_SCALARS: set[str] = set()


def _scalar_names(model) -> set[str]:
    if not _SCALARS:
        for n in cat.PARAMS:
            _SCALARS.add(model.param(n).get_name())
        for s in (V, CW, SW, RT2, tn.d, mv.parse("e")):
            _SCALARS.add(s.get_name())
    return _SCALARS


def _is_scalar(f: Expression, names: set[str]) -> bool:
    t = f.get_type()
    if t == AtomType.Num:
        return True
    if t == AtomType.Var:
        return f.get_name() in names or f.get_name() in ("symbolica::𝑖",)
    if t == AtomType.Pow:
        b, ex = list(f)
        return ex.get_type() == AtomType.Num and _is_scalar(b, names)
    return False


def _split_term(term: Expression, names: set[str]) -> tuple[Expression, Expression]:
    """``(scalar, structure)`` of a product term of a canonical density expression."""
    factors = list(term) if term.get_type() == AtomType.Mul else [term]
    sc, st = Expression.num(1), Expression.num(1)
    for f in factors:
        if _is_scalar(f, names):
            sc = sc * f
        else:
            st = st * f
    return sc, st


def structures(F: Density, names: set[str]) -> dict[str, tuple[Expression, Expression]]:
    """``{structure key: (structure, scalar)}`` of an (expanded, canonical) density."""
    out: dict[str, list] = {}
    e = F.expr.expand()
    terms = list(e) if e.get_type() == AtomType.Add else ([] if e == Expression.num(0) else [e])
    for t in terms:
        sc, st = _split_term(t, names)
        k = st.to_canonical_string()
        if k in out:
            out[k][1] = out[k][1] + sc
        else:
            out[k] = [st, sc]
    return {k: (v[0], v[1]) for k, v in sorted(out.items())}


def _numeric_part(e: Expression) -> Expression:
    if e.get_type() == AtomType.Num:
        return e
    if e.get_type() == AtomType.Mul:
        for f in e:
            if f.get_type() == AtomType.Num:
                return f
    return Expression.num(1)


# ---------------------------------------------------------------------------------------------- images
@dataclass
class Image:
    """The image of one parent's operator on one broken-phase field content."""

    parent: str
    content: tuple[str, ...]
    k: int                                   # power of v
    read: LocalFunctional                    # readable form (projectors kept)
    canon: dict                              # structure dict of the split canonical form
    read4: LocalFunctional | None = None     # image of the four-dimensional operator (finite part)
    canon4: dict | None = None


def _restrict(F: LocalFunctional, cls) -> LocalFunctional:
    return cat.restrict(F, cls)


def _canon(F: LocalFunctional) -> LocalFunctional:
    return LocalFunctional.integrate(F.split_lorentz(bf.simplifier()))


def _read(F: LocalFunctional) -> LocalFunctional:
    return cat.tidy(LocalFunctional.integrate(F.simplify(bf.display_simplifier())))


def _scalar_legs(spec_cls) -> int:
    return sum(1 for x in spec_cls if x in ("H", "Hdag"))


def parent_images(ctx: Ctx, pid: str) -> dict[tuple[str, ...], Image]:
    spec = ctx.specs[pid]
    coeffs = ctx.coeffs[pid]
    names = _scalar_names(ctx.model)
    img = bf.shift(spec.op)
    canon_all = _canon(img)
    read_all = _read(img)
    img4 = None
    if "finite" in coeffs and spec.op4 is not None:
        img4 = bf.shift(spec.op4)
        canon4_all, read4_all = _canon(img4), _read(img4)
    nh = _scalar_legs(spec.cls)
    out = {}
    for c in read_all.leg_classes():
        c = tuple(c)
        read = _restrict(read_all, c)
        canon = structures(_restrict(canon_all, c), names)
        if not canon:
            continue
        k = nh - sum(1 for x in c if x in ("h", "G0", "Gp", "Gm"))
        im = Image(pid, c, k, read, canon)
        if img4 is not None:
            r4 = _restrict(read4_all, c)
            if not r4.is_zero():
                im.read4 = r4
                im.canon4 = structures(_restrict(canon4_all, c), names)
        out[c] = im
    # contents that only the four-dimensional image reaches (cannot happen for P(O) = O4; checked)
    if img4 is not None:
        for c in read4_all.leg_classes():
            if tuple(c) not in out and structures(_restrict(canon4_all, tuple(c)), names):
                raise RuntimeError(f"{pid}: the four-dimensional operator reaches {c}, the operator does not")
    return out


# ---------------------------------------------------------------------------------------------- proportionality
def ratio(a: dict, b: dict) -> Expression | None:
    """``r`` with ``a = r b`` (structure dicts), or ``None``."""
    if set(a) != set(b):
        return None
    k0 = next(iter(a))
    r = (a[k0][1] / b[k0][1])
    for k in a:
        if not bf.is_zero(a[k][1] - r * b[k][1]):
            return None
    return bf.simplest(r)


def same(a: dict, b: dict) -> bool:
    """Equality of two structure dicts modulo the relations of :func:`gdocs.broken_fields.canonical_value`."""
    for k in set(a) | set(b):
        x = a[k][1] if k in a else Expression.num(0)
        y = b[k][1] if k in b else Expression.num(0)
        if not bf.is_zero(x - y):
            return False
    return True


def _scale_dict(a: dict, r: Expression) -> dict:
    return {k: (st, sc * r) for k, (st, sc) in a.items()}


# ---------------------------------------------------------------------------------------------- broken terms
@dataclass
class BTerm:
    """One operator of the broken-phase Lagrangian."""

    content: tuple[str, ...]
    op: LocalFunctional                     # operator (readable), multiplies tree / pole / evanescent / breaking
    canon: dict                             # its split canonical structure dict
    parents: list[tuple[str, Expression, int]]   # (unbroken id, lambda_P, power of v)
    op4: LocalFunctional | None = None      # four-dimensional operator of the finite part
    canon4: dict | None = None
    coeffs: dict = field(default_factory=dict)    # part -> Expression (after the vacuum condition)
    raw: dict = field(default_factory=dict)       # part -> Expression (before)
    merged: list = field(default_factory=list)    # [(operator, kappa)] of the merged operators (remarks)
    basis: list = field(default_factory=list)     # the parents whose sum gives the coefficients (before merging)
    groups: list = field(default_factory=list)    # merged terms: [(chirality part, its parents)] (remarks)
    elementary: bool = False
    identical: bool = False                       # the unbroken term itself (no SU(2)xU(1) field, no H)
    id: str = ""
    name: str = ""
    kind: str = ""
    sector: str = ""
    theory: str = ""
    name_src: str = ""

    @property
    def parent_ids(self) -> list[str]:
        out = []
        for p, _, _ in self.parents:
            if p not in out:
                out.append(p)
        return out


def _normalise(read: LocalFunctional, names: set[str]) -> tuple[LocalFunctional, Expression]:
    """``(O, lambda0)`` with ``read = lambda0 O``. ``lambda0`` is the scalar of the first structure of the readable
    form without its numerical factor (the whole scalar if it is a sum); a two-point operator with a single
    structure is normalised to the numerical factor ``1/2`` for two equal fields and ``1`` otherwise (mass terms
    ``M^2 (1/2) Z.Z``, ``M_W^2 W+.W-``)."""
    st = structures(read, names)
    s0 = next(iter(st.values()))[1].expand()
    lam = s0 if s0.get_type() == AtomType.Add else (s0 / _numeric_part(s0))
    if len(st) == 1 and s0.get_type() != AtomType.Add and len(read.monomials()[0].occurrences) == 2:
        occ = read.monomials()[0].occurrences
        sym = Expression.num(1) / 2 if occ[0].field == occ[1].field else Expression.num(1)
        lam = s0 / sym
    O = LocalFunctional(BT, (read.expr / lam).expand())
    O = LocalFunctional.integrate(O.map_coefficients(bf.reduce_sqrt2))
    return O, lam


def _fp_rank(fp: Fingerprinter, content, functionals: Sequence[LocalFunctional]) -> int:
    smp = fp.samples(list(content), len(functionals) + 3, salt=7)
    A = np.stack([fp.functional(f, smp) for f in functionals], axis=1)
    return linalg.rank(A, fp.field)


def collect(ctx: Ctx, images: dict[tuple, list[Image]], fp: Fingerprinter) -> list[BTerm]:
    """Steps 2 and 3 of the module docstring (before merging and the vacuum condition)."""
    names = _scalar_names(ctx.model)
    terms: list[BTerm] = []
    for content, ims in images.items():
        groups: list[list[tuple[Image, Expression]]] = []
        for im in ims:
            for gr in groups:
                r = ratio(im.canon, gr[0][0].canon)
                if r is not None:
                    gr.append((im, r))
                    break
            else:
                groups.append([(im, Expression.num(1))])
        reps = [gr[0][0].read for gr in groups]
        independent = len(reps) == 1 or _fp_rank(fp, content, reps) == len(reps)
        if independent:
            for gr in groups:
                im0 = gr[0][0]
                if len(gr) == 1 and tuple(content) == tuple(ctx.specs[im0.parent].cls):
                    # not touched by the shift or the rotation: the operator of the unbroken catalogue
                    t = BTerm(content, im0.read, im0.canon, [(im0.parent, Expression.num(1), im0.k)],
                              identical=True)
                    _attach_op4(t, gr, names)
                    terms.append(t)
                    continue
                O, lam0 = _normalise(gr[0][0].read, names)
                parents = [(im.parent, bf.simplest(r * lam0), im.k) for im, r in gr]
                t = BTerm(content, O, structures(_canon(O), names), parents)
                _attach_op4(t, gr, names)
                terms.append(t)
        else:
            terms += _elementary(ctx, content, ims, names)
    return terms


def _four_dimensional(O: LocalFunctional) -> LocalFunctional:
    """The four-dimensional form ``P O`` (all objects barred) of a bosonic operator; a fermion operator is its own
    four-dimensional form (the barred vertex), as in the unbroken catalogue."""
    if any(BT[o.field].kind in ("fermion", "fermion_bar") for m in O.monomials() for o in m.occurrences):
        return O
    return cat.tidy(LocalFunctional.integrate(O.physical(bf.display_simplifier())))


def _attach_op4(t: BTerm, gr: list[tuple[Image, Expression]], names) -> None:
    """The four-dimensional operator of the finite part, ``O4 = P O``; checked: ``image(O4_P) = lambda_P O4``
    for every parent with a finite part (raises otherwise)."""
    lam = {p: l for p, l, _ in t.parents}
    with4 = [im for im, _ in gr if im.read4 is not None]
    if not with4:
        return
    O4 = _four_dimensional(t.op)
    c4 = structures(_canon(O4), names)
    for im in with4:
        r = ratio(im.canon4, c4)
        if r is None or not bf.equal(r, lam[im.parent]):
            raise RuntimeError(f"{t.content}: the four-dimensional image of {im.parent} is not lambda_P P(O)")
    t.op4, t.canon4 = O4, c4


def _elementary(ctx: Ctx, content, ims: list[Image], names) -> list[BTerm]:
    """Fallback: the content written on the elementary monomials of the (unsplit) canonical images."""
    simp = bf.simplifier()
    keys: dict[str, Expression] = {}
    per: list[tuple[Image, dict]] = []
    for im in ims:
        st = structures(LocalFunctional.integrate(im.read.simplify(simp)), names)
        per.append((im, st))
        for k, (s, _) in st.items():
            keys.setdefault(k, s)
    out = []
    for k, s in keys.items():
        O = LocalFunctional(BT, s)
        parents = [(im.parent, bf.simplest(st[k][1]), im.k) for im, st in per if k in st]
        t = BTerm(content, O, structures(_canon(O), names), parents, elementary=True)
        # finite part: the four-dimensional images on the projected monomials
        with4 = [im for im, _ in per if im.read4 is not None]
        if with4:
            O4 = _four_dimensional(O)
            t.op4, t.canon4 = O4, structures(_canon(O4), names)
        out.append(t)
    # check: every image is reproduced by the monomials (and its four-dimensional image by their projections)
    for im, st in per:
        rebuilt = Expression.num(0)
        rebuilt4 = Expression.num(0)
        for t in out:
            for p, l, _ in t.parents:
                if p == im.parent:
                    rebuilt = rebuilt + l * t.op.expr
                    if t.op4 is not None:
                        rebuilt4 = rebuilt4 + l * t.op4.expr
        if not same(structures(_canon(LocalFunctional(BT, rebuilt)), names), im.canon):
            raise RuntimeError(f"{content}: elementary decomposition of {im.parent} failed")
        if im.read4 is not None and not same(structures(_canon(LocalFunctional(BT, rebuilt4)), names), im.canon4):
            raise RuntimeError(f"{content}: four-dimensional elementary decomposition of {im.parent} failed")
    return out


def _coefficients(ctx: Ctx, t: BTerm) -> None:
    raw: dict[str, Expression] = {}
    for pid, lam, _ in t.parents:
        for part, c in ctx.coeffs[pid].items():
            raw[part] = raw.get(part, Expression.num(0)) + lam * c
    t.raw = raw
    t.coeffs = {}
    for part in PARTS:
        if part not in raw:
            continue
        val = bf.vev(raw[part])
        if bf.is_zero(val):
            continue
        t.coeffs[part] = bf.simplest(val)


def _merge(terms: list[BTerm]) -> list[BTerm]:
    """Sum the operators of one content whose coefficients are proportional by a number in every part."""
    by_content: dict[tuple, list[BTerm]] = {}
    for t in terms:
        by_content.setdefault(t.content, []).append(t)
    out = []
    for content, ts in by_content.items():
        done = [False] * len(ts)
        for i, a in enumerate(ts):
            if done[i]:
                continue
            for j in range(i + 1, len(ts)):
                b = ts[j]
                if done[j] or set(a.coeffs) != set(b.coeffs) or not a.coeffs:
                    continue
                if (a.op4 is None) != (b.op4 is None) or a.elementary or b.elementary:
                    continue
                if _is_ev(a) != _is_ev(b) or not _chirality_partners(a.op, b.op):
                    continue
                kappa = None
                ok = True
                for part in a.coeffs:
                    r = bf.canonical_value(a.coeffs[part] / b.coeffs[part])
                    if r.get_type() != AtomType.Num:
                        ok = False
                        break
                    if kappa is None:
                        kappa = r
                    elif not bf.equal(r, kappa):
                        ok = False
                        break
                if not ok:
                    continue
                # a.coeffs = kappa * b.coeffs: a O_a + b O_b = a (O_a + O_b / kappa)
                if not a.basis:
                    a.basis = list(a.parents)
                if not a.groups:
                    a.groups = [(a.op, list(a.parents))]
                a.groups.append((b.op, [(p, l * kappa, k) for p, l, k in b.parents]))
                a.merged.append((b.op, kappa))
                a.op = LocalFunctional.integrate(a.op + b.op.scale(1 / kappa))
                if a.op4 is not None:
                    a.op4 = LocalFunctional.integrate(a.op4 + b.op4.scale(1 / kappa))
                a.parents = a.parents + [(p, l * kappa, k) for p, l, k in b.parents]
                done[j] = True
            out.append(a)
    return out


def _swap_chirality(F: LocalFunctional) -> LocalFunctional:
    x = tn.HEADS["PL"]
    y = tn.HEADS["PR"]
    tmp = S("gdocs_broken::swap")
    a, b = S("gdocs_broken::a_", "gdocs_broken::b_")
    e = F.expr.replace(x(a, b), tmp(a, b)).replace(y(a, b), x(a, b)).replace(tmp(a, b), y(a, b))
    return LocalFunctional(BT, e)


def _chirality_partners(a: LocalFunctional, b: LocalFunctional) -> bool:
    """``b`` is proportional to ``a`` with ``P_L`` and ``P_R`` exchanged (the two chiralities of one fermion
    bilinear, e.g. the images of a Yukawa monomial and of its conjugate)."""
    if not any(BT[o.field].kind in ("fermion", "fermion_bar") for m in a.monomials() for o in m.occurrences):
        return False
    names = _SCALARS
    sa = structures(_canon(_swap_chirality(a)), names)
    sb = structures(_canon(b), names)
    return ratio(sa, sb) is not None


_EVANESCENT: dict[str, bool] = {}


def _is_ev(t: BTerm) -> bool:
    """The operator is evanescent: all its parents are evanescent operators of the unbroken catalogue."""
    return all(_EVANESCENT.get(p, False) for p in t.parent_ids)


# ---------------------------------------------------------------------------------------------- display form
def from_site(e: Expression) -> Expression:
    """Inverse of :func:`gdocs.broken_fields.to_site` (site parameter symbols -> package symbols)."""
    for n in bf._SITE_PARAMS:
        e = e.replace(mv.parse(n), bf.unbroken_model().param(n))
    return e


def _n_mono(F: LocalFunctional) -> int:
    return len(F.monomials())


def _display_candidates(O: LocalFunctional) -> list[LocalFunctional]:
    """Equivalent forms of an operator: readable (projectors kept), canonical (projectors as (1 -+ gamma5)/2),
    Lorentz indices split into barred and hatted parts (readable)."""
    out = [O]
    try:
        out.append(LocalFunctional.integrate(O.simplify(bf.simplifier())))
    except Exception:  # noqa: BLE001
        pass
    try:
        out.append(cat.tidy(LocalFunctional.integrate(O.split_lorentz(bf.display_simplifier()))))
    except Exception:  # noqa: BLE001
        pass
    return [c for c in out if not c.is_zero()]


def _display_scalar(O: LocalFunctional) -> Expression:
    """Scalar (site notation) of the first displayed term of ``O``."""
    parts = mv.operator_parts(O, table=BT)
    return parts[0][0]


def _num_of(e: Expression) -> tuple[Expression, bool]:
    """``(rational factor, imaginary)`` of a site scalar (product of a number, ``i`` and symbols)."""
    factors = list(e) if e.get_type() == AtomType.Mul else [e]
    num = Expression.num(1)
    imag = False
    for f in factors:
        if f.get_type() == AtomType.Num:
            s = f.to_canonical_string()
            if "𝑖" in s:
                imag = True
                num = num * (f / Expression.I)
            else:
                num = num * f
    return num, imag


def finalise(t: BTerm) -> None:
    """Choose the display form of the operator (fewest monomials) and normalise it: a single displayed term gets
    the numerical factor ``1/2`` (two equal fields) or ``1`` (otherwise), times ``i`` if it is imaginary; several
    terms keep the normalisation of the first parent. The coefficients and the parents' factors are rescaled."""
    if t.identical:
        return
    cands = _display_candidates(t.op)
    best = min(cands, key=_n_mono)
    parts = mv.operator_parts(best, table=BT)
    s_disp = parts[0][0]
    num, imag = _num_of(s_disp)
    if len(parts) == 1:
        occ = best.monomials()[0].occurrences
        half = len(occ) == 2 and occ[0].field == occ[1].field
        target = (Expression.num(1) / 2 if half else Expression.num(1)) * (Expression.I if imag else 1)
    else:
        target = num * (Expression.I if imag else 1)
    mu = bf.simplest(from_site(s_disp / target))
    if bf.is_zero(mu - 1):
        t.op = best
    else:
        t.op = LocalFunctional.integrate(LocalFunctional(BT, (best.expr / mu).expand()).map_coefficients(
            bf.reduce_sqrt2))
        if t.op4 is not None:
            t.op4 = LocalFunctional.integrate(LocalFunctional(BT, (t.op4.expr / mu).expand()).map_coefficients(
                bf.reduce_sqrt2))
        t.coeffs = {p: bf.simplest(c * mu) for p, c in t.coeffs.items()}
        t.raw = {p: (c * mu) for p, c in t.raw.items()}
        t.parents = [(p, bf.simplest(l * mu), k) for p, l, k in t.parents]
        t.basis = [(p, bf.simplest(l * mu), k) for p, l, k in t.basis]
        t.groups = [(op, [(p, bf.simplest(l * mu), k) for p, l, k in pl]) for op, pl in t.groups]
    if t.op4 is not None:
        c4 = _display_candidates(t.op4)
        t.op4 = min(c4, key=_n_mono)


# ---------------------------------------------------------------------------------------------- classification
FIELD_TEXT = {"G": "G", "Wp": "W⁺", "Wm": "W⁻", "Z": "Z", "A": "γ", "h": "h", "G0": "G⁰", "Gp": "G⁺", "Gm": "G⁻",
              "u": "u", "ubar": "ū", "up": "u′", "upbar": "ū′", "d": "d", "dbar": "d̄", "dp": "d′", "dpbar": "d̄′",
              "e": "e", "ebar": "ē", "ep": "e′", "epbar": "ē′", "nu": "ν", "nubar": "ν̄",
              "src_H_1": "σ†₁", "src_H_2": "σ†₂", "src_Hdag_1": "σ₁", "src_Hdag_2": "σ₂"}
_SUP = {"Wp": "⁺", "Wm": "⁻", "Z": "ᶻ", "A": "ᵞ", "G": ""}
for _x in ("Wp", "Wm", "Z", "A", "G"):
    _s = _SUP[_x]
    FIELD_TEXT[f"gh_{_x}"] = "ω" + (_s or "")
    FIELD_TEXT[f"agh_{_x}"] = "ω̄" + (_s or "")
    FIELD_TEXT[f"b_{_x}"] = "b" + (_s or "")
    FIELD_TEXT[f"src_{_x}"] = "ρ" + (_s or "")
    FIELD_TEXT[f"src_gh_{_x}"] = "ζ" + (_s or "")
for _f in "QUDLE":
    for _i in ("", "_1", "_2"):
        FIELD_TEXT[f"src_{_f}{_i}"] = f"R̄_{_f}" + {"": "", "_1": "¹", "_2": "²"}[_i]
        FIELD_TEXT[f"src_{_f}bar{_i}"] = f"R_{_f}" + {"": "", "_1": "¹", "_2": "²"}[_i]

_DISPLAY_RANK = {"antighost": 0, "ghost": 1, "fermion_bar": 2, "fermion": 3, "scalar_conj": 5, "scalar": 4,
                 "gauge": 6, "auxiliary": 7, "source": 8}


def content_source(content: Sequence[str]) -> str:
    """Prose source of a field content: the chip symbols of :func:`chip_symbol` as inline math, in the order of
    :func:`content_text`."""
    return " ".join(f"[[{mv.M(chip_symbol(x))['src']}]]" for x in _content_order(content))


def content_text(content: Sequence[str]) -> str:
    return " ".join(FIELD_TEXT.get(x, x) for x in _content_order(content))


def _content_order(content: Sequence[str]) -> list[str]:
    def rank(x):
        k = BT[x].kind
        if k == "source":
            sk = BT[BT[x].source_of].kind if BT[x].source_of in BT else ""
            return -1 if sk in ("gauge", "ghost") else (2 if sk == "fermion" else (3.5 if sk == "fermion_bar"
                                                                                    else 4.5))
        return _DISPLAY_RANK.get(k, 9)
    return sorted(content, key=rank)


FAMILY = [
    (r"^kin-[GWB]$", "Gauge kinetic term"), (r"^long-", "Longitudinal two-point term"),
    (r"^ev-EA-", "Evanescent operator E_A"), (r"^ev-EtF-", "Evanescent operator Ẽ_F"),
    (r"^cubic-", "Cubic gauge vertex"), (r"^quartic-", "Quartic gauge operator"),
    (r"^mixed-", "Mixed quartic gauge operator"), (r"^kin-[QUDLE]$", "Fermion kinetic term"),
    (r"^sterile-", "Kinetic term of the sterile partner"), (r"^evkin-", "Evanescent fermion kinetic term"),
    (r"^vertex-", "Fermion–gauge vertex"), (r"^higgs-kin$", "Scalar kinetic term"),
    (r"^higgs-mass$", "Scalar mass term"), (r"^higgs-quartic$", "Scalar self-coupling"),
    (r"^higgs-current-", "Scalar current term"), (r"^higgs-(WW|WB|BB)$", "Scalar–gauge seagull"),
    (r"^higgs-GG$", "Scalar–gluon term"), (r"^ev-EH$", "Evanescent operator E_H"),
    (r"^ev-EtH-", "Evanescent operator Ẽ_H"), (r"^yukawa-", "Yukawa term"), (r"^gf-bA-", "Gauge fixing b ∂·A"),
    (r"^gf-bb-", "Gauge fixing ξ/2 b b"), (r"^ghost-kin-", "Ghost kinetic term"),
    (r"^ghost-gauge-", "Ghost–gauge vertex"), (r"^src-rho-", "Source term ρ ∂ω"),
    (r"^src-rhoA-", "Source term ρ A ω"), (r"^src-zeta-", "Source term ζ ω ω"),
    (r"^src-R", "Fermion source term"), (r"^src-sigma", "Higgs source term"), (r"^Delta1-", "Breaking Δ₁"),
]


def family(pid: str) -> str:
    for pat, lab in FAMILY:
        if re.search(pat, pid):
            return lab
    return pid


KINDS = [("tadpole", "Tadpole"), ("mass", "Mass"), ("mixing", "Mixing"), ("kinetic", "Kinetic"),
         ("interaction", "Interaction"), ("source", "Source"), ("evanescent", "Evanescent"),
         ("unchanged", "Unchanged by the shift")]
KIND_NAME = {"tadpole": "Tadpole", "mass": "Mass term", "mixing": "Mixing term", "kinetic": "Kinetic term",
             "interaction": "Vertex", "source": "Source term", "evanescent": "Evanescent term"}
SCALAR_FIELDS = ("h", "G0", "Gp", "Gm")


def _has_derivative(F: LocalFunctional) -> bool:
    return any(o.derivatives for m in F.monomials() for o in m.occurrences)


def _shifted(ctx: Ctx, pid: str) -> bool:
    return any(x in ("H", "Hdag") for x in ctx.specs[pid].cls)


def classify(ctx: Ctx, t: BTerm) -> None:
    """``kind``, ``sector`` and ``theory`` of a broken-phase term."""
    content = t.content
    kinds = {BT[x].kind for x in content}
    shifted = any(_shifted(ctx, p) for p in t.parent_ids)
    n = len(content)
    if not shifted:
        kind = "unchanged"
    elif "source" in kinds:
        kind = "source"
    elif _is_ev(t):
        kind = "evanescent"
    elif n == 1:
        kind = "tadpole"
    elif n == 2:
        a, b = content
        pair = a == b or BT[a].conjugate == b
        kind = ("kinetic" if _has_derivative(t.op) else "mass") if pair else "mixing"
    else:
        kind = "interaction"
    t.kind = kind
    fams = {ctx.specs[p].family for p in t.parent_ids}
    if "source" in kinds:
        sector = "source"
    elif "breaking" in fams:
        sector = "breaking"
    elif kinds & {"ghost", "antighost"}:
        sector = "ghost"
    elif "auxiliary" in kinds:
        sector = "gauge-fixing"
    elif kinds & {"fermion", "fermion_bar"} and kinds & {"scalar", "scalar_conj"}:
        sector = "yukawa"
    elif kinds & {"fermion", "fermion_bar"}:
        sector = "fermion"
    elif kinds & {"scalar", "scalar_conj"}:
        sector = "scalar"
    else:
        sector = "gauge"
    t.sector = sector
    finite_only = set(t.coeffs) == {"finite"}
    if kind == "unchanged":
        theory = ctx.terms[t.parent_ids[0]]["theory"]
    elif kind == "tadpole":
        theory = "broken-tadpole"
    elif kind == "mass":
        theory = "broken-st" if finite_only else "broken-masses"
    elif kind == "mixing":
        theory = "broken-mixing"
    elif kind == "kinetic":
        theory = "broken-shift"
    elif kind == "evanescent":
        if n == 2:
            a, b = content
            pair = a == b or BT[a].conjugate == b
            theory = ("broken-shift" if _has_derivative(t.op) else "broken-masses") if pair else "broken-mixing"
        else:
            theory = "broken-interactions"
    elif kind == "interaction":
        theory = "broken-st" if finite_only else "broken-interactions"
    else:
        theory = "broken-interactions"
    t.theory = theory


def _chirality(t: BTerm) -> str | None:
    """``left-handed`` / ``right-handed`` if the operator contains only ``P_L`` / only ``P_R``."""
    s = t.op.expr.to_canonical_string()
    pl, pr = "::PL(" in s, "::PR(" in s
    if pl and not pr:
        return "P_L"
    if pr and not pl:
        return "P_R"
    return None


def name_terms(ctx: Ctx, terms: list[BTerm]) -> None:
    """Names and ids (unique, stable: field content, kind and the first parent)."""
    for t in terms:
        p0 = t.parent_ids[0]
        ctext = content_text(t.content)
        if t.kind == "unchanged":
            if t.identical:
                t.name, t.id = ctx.terms[p0]["name"], p0
                continue
            base = family(p0)
        else:
            base = "Yukawa vertex" if (t.kind == "interaction" and t.sector == "yukawa") else KIND_NAME[t.kind]
        t.name = f"{base}: {ctext}"
        tag = "_".join(t.content) or "none"
        t.id = f"{p0}--{tag}" if t.kind == "unchanged" else f"{t.kind}--{tag}"
    # disambiguation: chirality, then the family of the first parent, then a structure number
    for step in range(3):
        seen: dict[str, list[BTerm]] = {}
        for t in terms:
            seen.setdefault(t.name, []).append(t)
        for nm, ts in seen.items():
            if len(ts) < 2:
                continue
            for i, t in enumerate(ts):
                if step == 0:
                    ch = _chirality(t)
                    if ch:
                        fam = family(t.parent_ids[0])
                        if fam in ("Fermion kinetic term", "Fermion–gauge vertex"):
                            t.name += f" ({ {'P_L': 'left-handed', 'P_R': 'right-handed'}[ch] })"
                        else:
                            t.name += f" ({ {'P_L': 'with P_L', 'P_R': 'with P_R'}[ch] })"
                elif step == 1:
                    if t.kind != "unchanged":
                        t.name += f" ({family(t.parent_ids[0]).lower()})"
                    else:
                        pname = ctx.terms[t.parent_ids[0]]["name"]
                        struct = pname.split(" operator ", 1)[1] if " operator " in pname else pname.lower()
                        t.name += f" (from {struct})"
                else:
                    t.name += f" (structure {i + 1})"
    seen_ids: dict[str, int] = {}
    for t in terms:
        if t.id in seen_ids:
            seen_ids[t.id] += 1
            alt = f"{t.id}--{t.parent_ids[0]}"
            t.id = alt if alt not in seen_ids else f"{t.id}-{seen_ids[t.id]}"
        seen_ids.setdefault(t.id, 0)
    ids = [t.id for t in terms]
    if len(set(ids)) != len(ids):
        raise RuntimeError("duplicate broken-phase term ids")
    for t in terms:
        t.name_src = name_source(ctx, t)
    for i in ids:
        if not re.fullmatch(r"[A-Za-z0-9_-]+", i):
            raise RuntimeError(f"invalid term id {i}")


def name_source(ctx: Ctx, t: BTerm) -> str:
    """Prose source of the name of a term: the field content as chip symbols, the rest through
    :func:`gdocs.names.name_source`."""
    if t.kind == "unchanged" and t.identical:
        return nm.name_source(t.name)
    base, sep, rest = t.name.partition(": ")
    ctext = content_text(t.content)
    if not sep or not rest.startswith(ctext):
        return nm.name_source(t.name)
    return nm.name_source(base) + ": " + content_source(t.content) + nm.name_source(rest[len(ctext):])


# ---------------------------------------------------------------------------------------------- math objects
def _site(c: Expression) -> Expression:
    return bf.for_display(bf.to_site(c))


_VIEW_THRESHOLD = 6000   # static MathML longer than this gets an interactive view
_VIEW_BUDGET = 12000


def _M(e: Expression, **kw) -> dict:
    kw.setdefault("threshold", _VIEW_THRESHOLD)
    kw.setdefault("budget", _VIEW_BUDGET)
    return cat._M(bf.for_display(e), **kw)


def _coef_math(c: Expression, part: str) -> Expression:
    return cat._with_loop(bf.to_site(c), part)


_DPAR_M2 = mv.parse("dpar(m2)")
_HOLD = mv.sym("gdocs_broken_hold_dm2")


def _site_vev(e: Expression) -> Expression:
    """``m^2 = -lambda v^2`` in a site expression, keeping ``dpar(m2)``."""
    e = e.replace(_DPAR_M2, _HOLD)
    e = e.replace(mv.parse("m2"), -mv.parse("lambda") * V ** 2)
    return e.replace(_HOLD, _DPAR_M2)


_EV_GEN = {"E_A": "-1/3*Ssum({r})", "Etilde_F": "-2/3*Ssum({r})", "E_H": "-2/3*Ytil", "Etilde_H": "-Ytil"}


def _parent_decomposition(ctx: Ctx, pid: str, part: str) -> tuple[Expression, Expression] | None:
    """``(site expression, value)`` of the singular coefficient of a parent through the constants."""
    spec = ctx.specs[pid]
    co = ctx.coeffs[pid]
    if part == "pole":
        if not spec.renormalised or "tree" not in co or "pole" not in co:
            return None
        return cat.decomposition(ctx.c, spec, co["tree"])
    if part == "evanescent" and spec.evname:
        base, _, r = spec.evname.partition("^")
        return mv.parse(_EV_GEN[base].format(r=r.strip("()"))), co["evanescent"]
    return None


def decomposition(ctx: Ctx, t: BTerm, part: str) -> Expression | None:
    """The singular coefficient of ``t`` as the sum over its parents of ``lambda_P`` times their decomposition
    through the renormalisation constants (checked against the coefficient)."""
    basis = t.basis or t.parents
    total, val = Expression.num(0), Expression.num(0)
    n = 0
    for pid, lam, _ in basis:
        if part not in ctx.coeffs[pid]:
            continue
        d = _parent_decomposition(ctx, pid, part)
        if d is None:
            return None
        site, v = d
        lam_site = bf.to_site(lam)
        total = total + lam_site * (mv.parse("par(x)").replace(mv.parse("x"), site)
                                    if site.get_type() == AtomType.Add else site)
        val = val + lam * v
        n += 1
    if not n:
        return None
    if not bf.equal(bf.vev(val), t.coeffs[part]):
        raise RuntimeError(f"{t.id}: the decomposition through the constants does not reproduce the {part} part")
    total = _site_vev(total)
    expanded = _site_vev(total.replace(mv.parse("par(x_)"), mv.parse("x_")).expand())
    if len(expanded.format_plain()) < len(total.format_plain()):
        total = expanded
    return total


def _differs(t: BTerm) -> bool:
    if t.op4 is None or "finite" not in t.coeffs:
        return False
    a = structures(_canon(t.op), _SCALARS)
    b = structures(_canon(t.op4), _SCALARS)
    return not same(a, b)


def _total(t: BTerm, types: list[str], differs: bool) -> Expression:
    op_site = bf.for_display(mv.operator(t.op, table=BT))
    op4_site = bf.for_display(mv.operator(t.op4, table=BT)) if t.op4 is not None and differs else None

    def neg(e):
        return cat._negative(e) or (e.get_type() == AtomType.Mul and any(
            f.get_type() == AtomType.Num and f.to_canonical_string() == "-1" for f in e))

    def ordered_sum(items):
        acc = None
        for e in items:
            if acc is None:
                acc = e
            elif neg(e):
                acc = mv.parse("minus(x, y)").replace(mv.parse("x"), acc).replace(mv.parse("y"), -e)
            else:
                acc = mv.parse("plus(x, y)").replace(mv.parse("x"), acc).replace(mv.parse("y"), e)
        return acc

    def times(c, op, several):
        if several:
            c = mv.parse("par(x)").replace(mv.parse("x"), c)
        if cat._negative(op) or op.get_type() == AtomType.Add:
            op = mv.parse("par(x)").replace(mv.parse("x"), op)
        if several:
            # a sum of several coefficients follows the operator, so that the operator stays in view
            return mv.parse("ord(y, x)").replace(mv.parse("x"), c).replace(mv.parse("y"), op)
        return mv.parse("ord(x, y)").replace(mv.parse("x"), c).replace(mv.parse("y"), op)

    same_ = [p for p in types if not (p == "finite" and differs)]
    if same_:
        cs = [bf.for_display(_coef_math(t.coeffs[p], p)) for p in same_]
        if len(same_) == 1 and bf.is_zero(t.coeffs[same_[0]] - 1):
            total = op_site
        else:
            total = times(ordered_sum(cs), op_site, len(same_) > 1)
    else:
        total = None
    if differs:
        cf = bf.for_display(_coef_math(t.coeffs["finite"], "finite"))
        n_ = neg(cf)
        fin = times(-cf if n_ else cf, op4_site, False)
        if total is None:
            total = -fin if n_ else fin
        else:
            total = mv.parse(f"{'minus' if n_ else 'plus'}(x, y)").replace(mv.parse("x"), total).replace(
                mv.parse("y"), fin)
    if "⋯" in mv.M(total)["html"]:
        # too large for one of Symbolica's custom-printed fragments: one product per part
        total = Expression.num(0)
        for p in types:
            o = op4_site if (p == "finite" and differs) else op_site
            c = bf.for_display(_coef_math(t.coeffs[p], p))
            total = total - times(-c, o, False) if neg(c) else total + times(c, o, False)
        if "⋯" in mv.M(total)["html"]:
            total = Expression.num(0)
            for p in types:
                o = op4_site if (p == "finite" and differs) else op_site
                total = total + bf.for_display(_coef_math(t.coeffs[p], p)) * o
    return total


def _rule_legs(content) -> tuple[str, ...]:
    return tuple(sorted(content, key=lambda x: cat._leg_rank(BT, x)))


# ---------------------------------------------------------------------------------------------- remarks
def _link(pid: str, ctx: Ctx) -> str:
    return f'<a data-term="unbroken:{pid}">{nm.name_source(ctx.terms[pid]["name"])}</a>'


_EQ_CITE = re.compile(r"\beq\. ([A-Za-z][A-Za-z0-9-]*[A-Za-z0-9])")


def _prose(text: str) -> str:
    """``mathview.prose`` after turning the citations ``eq. LABEL`` (report equation labels) into links that show
    the equation number of the PDF."""
    return mv.prose(_EQ_CITE.sub(lambda m: _eqa(m.group(1)), text))


def _eq_plain(text: str) -> str:
    """Plain-text form of the same citations: ``Eq. (N), p. P``."""
    from . import theory

    def rep(m):
        r = theory.eq_ref(m.group(1))
        return f"{r['label']}, p. {r['page']}"
    return _EQ_CITE.sub(rep, text)


def _eqa(label: str) -> str:
    """Rich-text link to a numbered equation of the report (``Eq. (N)`` with its page)."""
    from . import theory

    return theory.eq_anchor(label)


def _factor_src(lam: Expression, k: int) -> str:
    return mv.M(_site(lam))["src"]


_MASSES = None


def _mass_forms():
    """(name, site source, value) of the tree-level masses of eq:bl-masses."""
    global _MASSES
    if _MASSES is None:
        m = bf.unbroken_model()
        g, lam = m.param("g"), m.param("lambda")
        _MASSES = [("M_W^2", g ** 2 * V ** 2 / 4), ("M_Z^2", g ** 2 * V ** 2 / (4 * CW ** 2)),
                   ("M_h^2", 2 * lam * V ** 2), ("M_W", g * V / 2), ("M_Z", g * V / (2 * CW))]
        for k in "ude":
            _MASSES.append((f"m_{k}", m.param(f"y_{k}") * V / RT2))
            _MASSES.append((f"m_{k}^2", m.param(f"y_{k}") ** 2 * V ** 2 / 2))
    return _MASSES


def mass_identification(c: Expression) -> str | None:
    """``r M`` if the tree-level coefficient ``c`` is a rational (or imaginary rational) multiple of a mass of
    eq:bl-masses; the site source of the product."""
    for name, val in _mass_forms():
        r = bf.canonical_value(c / val)
        if r.get_type() == AtomType.Num and r != Expression.num(0):
            return mv.M(mv.parse(name) * r)["src"]
    return None


SPECIAL = {
    ("h",): "The tree-level tadpole [[-v*par(m2 + lambda*v^2)*h]] vanishes by the vacuum condition "
            "[[eqs(m2, -lambda*v^2)]]. The singular tadpole term is [[-v*dpar(T)*h]] with "
            "[[eqs(dpar(T), dpar(m2) + v^2*dpar(lambda) + lambda*v^2*dZ(H))]] (report eq. bl-dT); it cancels the "
            "pole of the one-loop tadpole, whose tadpole graphs are kept. The report checked only the "
            "[[lambda^2]] part of this pole by hand, and it does not compute the finite one-loop tadpole. The "
            "finite restoring counterterm gives no tadpole.",
    ("G0", "G0"): "The Goldstone fields are massless at tree level in this gauge. The mass counterterm equals "
                  "[[dpar(T)]], as required by the relation between the Goldstone two-point function at zero "
                  "momentum and the tadpole once tadpoles are included (report Section 5.4.4).",
    ("Gp", "Gm"): "The Goldstone fields are massless at tree level in this gauge. The mass counterterm equals "
                  "[[dpar(T)]], as required by the relation between the Goldstone two-point function at zero "
                  "momentum and the tadpole (report Section 5.4.4).",
    ("A", "A"): "The finite photon mass-like term is required: at one loop and zero momentum the fermion loops of "
                "the regularised theory give the photon a finite mass term, and this counterterm cancels it "
                "exactly (report Section 5.4.6, eq. bl-massless; checked there also through the Slavnov–Taylor "
                "identity of the ghost–photon legs). There is no singular photon mass term.",
    ("G", "G"): "The finite gluon mass-like term cancels the finite mass term of the quark loops at zero momentum "
                "(report Section 5.4.6, eq. bl-massless).",
    ("Z", "A"): "Z–photon mixing of mass type. The fermion-loop part of the Z–photon two-point function vanishes at "
                "zero momentum once the finite term is included (report Section 5.4.6); the singular term comes "
                "from bosonic loops.",
    ("Z", "G0"): "The tree-level gauge–Goldstone mixing remains because the gauge fixing is linear; only an "
                 "H-dependent gauge fixing such as the R_ξ gauge, which is not used, would cancel it. The finite "
                 "part is [[-Ytil/3]] times the tree-level mixing (report eq. bl-mix).",
    ("Wp", "Gm"): "The tree-level gauge–Goldstone mixing remains because the gauge fixing is linear (report "
                  "eq. bl-L2mix). The finite part is [[-Ytil/3]] times the tree-level mixing (eq. bl-mix).",
    ("Wm", "Gp"): "The tree-level gauge–Goldstone mixing remains because the gauge fixing is linear (report "
                  "eq. bl-L2mix). The finite part is [[-Ytil/3]] times the tree-level mixing (eq. bl-mix).",
    ("A", "G0"): "There is no tree-level and no finite photon–Goldstone term (the two current coefficients "
                 "[[rc(J, up(W))]] and [[rc(J, up(B))]] are equal); the singular term is new in the broken phase "
                 "(report eq. bl-dS-mix).",
}
for _f in ("u", "d", "e"):
    SPECIAL[(_f, _f + "bar")] = None   # set below (mass terms only)


def _projection(op: LocalFunctional) -> str | None:
    """``PL`` / ``PR`` if the operator contains only that projector."""
    s = op.expr.to_canonical_string()
    pl, pr = "::PL(" in s, "::PR(" in s
    return "PL" if pl and not pr else ("PR" if pr and not pl else None)


def remarks(ctx: Ctx, t: BTerm) -> str | None:
    out: list[str] = []
    pids = t.parent_ids
    identical = t.identical
    if identical:
        out.append(f"Not changed by the shift: this is the term {_link(pids[0], ctx)} of the unbroken phase.")
    elif t.merged and t.groups:
        lead = "Not changed by the shift; written in the mass basis. The operator collects the two chirality " \
               "projections of one bilinear. " if t.kind == "unchanged" else \
            "The operator collects the two chirality projections of one bilinear. "
        parts = []
        for op, plist in t.groups:
            pr = _projection(op)
            bits = []
            for pid, lam, k in plist:
                vk = f", order [[v^{k}]]" if k > 1 else (", order [[v]]" if k == 1 else "")
                bits.append(f"{_link(pid, ctx)} (factor [[{_factor_src(lam, k)}]]{vk})")
            what = f"the [[{pr}]] projection" if pr else "one projection"
            parts.append(f"{what} is the image under the shift of " + " and ".join(bits))
        out.append(lead + "; ".join(parts)[0].upper() + "; ".join(parts)[1:] + ". For each projection the sum of "
                   "its parents' coefficients times these factors is the coefficient of the operator shown; the two "
                   "projections have the same coefficient.")
    else:
        bits = []
        for pid, lam, k in t.parents:
            vk = f", order [[v^{k}]]" if k > 1 else (", order [[v]]" if k == 1 else "")
            bits.append(f"{_link(pid, ctx)} (factor [[{_factor_src(lam, k)}]]{vk})")
        lead = "Not changed by the shift; written in the mass basis. It is the image of " if t.kind == "unchanged" \
            else "Image under the shift of "
        out.append(lead + "; ".join(bits) + ": each coefficient is the sum of the parents' coefficients times "
                   "these factors.")
    if any("m2" in str(ctx.coeffs[p].get(x, "")) for p in pids for x in ctx.coeffs[p]) and t.kind != "unchanged":
        out.append("The coefficients use the vacuum condition [[eqs(m2, -lambda*v^2)]].")
    if "tree" in t.coeffs and t.kind in ("mass", "mixing", "interaction"):
        mi = mass_identification(t.coeffs["tree"])
        if mi:
            out.append(f"With the tree-level masses of the report ({_eqa('bl-masses')}) the tree-level coefficient is "
                       f"[[{mi}]].")
    sp = SPECIAL.get(tuple(t.content))
    if sp and t.kind in ("tadpole", "mass", "mixing"):
        out.append(sp)
    if t.kind == "mass" and t.content in (("u", "ubar"), ("d", "dbar"), ("e", "ebar")):
        out.append("Only the physical components couple to the mass; the sterile partner enters the propagator "
                   f"through the evanescent kinetic term (report {_eqa('bp-prop')}). Under (N1) there is no finite "
                   "fermion mass counterterm.")
    if _differs(t) and ("tree" in t.coeffs or "pole" in t.coeffs):
        out.append("The tree-level and singular coefficients multiply the d-dimensional operator shown; the finite "
                   "coefficient multiplies its four-dimensional form, in which all fields, derivatives and metrics "
                   "are barred.")
    if t.elementary:
        out.append("The images of the parents on this field content are not proportional to each other; the "
                   "content is therefore written on its elementary monomials.")
    if t.kind == "source" and len(t.content) == 2:
        out.append("A two-point vertex between a source and a ghost; it has no propagator.")
    if t.kind in ("tadpole", "mass", "mixing", "kinetic", "interaction", "evanescent", "source"):
        out.append("The graphs shown are those of the parents in the unbroken phase; broken-phase graphs are not "
                   "drawn. Pole parts and the finite breaking of the broken phase are the unbroken-phase ones with "
                   "external [[H]] legs set to [[Hv]] at zero momentum (report Section 5.5.1, point 2).")
    return _prose("<p>" + " ".join(out) + "</p>") if out else None


def references(ctx: Ctx, t: BTerm) -> list[dict]:
    refs = []
    if t.kind == "unchanged":
        if not t.identical:
            refs.append(_ref("sec:bl-fields"))
        refs += ctx.terms[t.parent_ids[0]]["report"][:3]
    else:
        labs = {"tadpole": ["sec:bl-fields", "sec:bl-dS", "sec:bp-tadpole"],
                "mass": ["sec:bl-fields", "sec:bl-tree", "sec:bl-dS"],
                "mixing": ["sec:bl-tree", "sec:bl-dS"], "kinetic": ["sec:bl-tree", "sec:bl-dS"],
                "interaction": ["sec:bl-tree", "sec:bl-dS"], "source": ["sec:bl-statement", "sec:bl-dS"],
                "evanescent": ["sec:bl-dS"]}[t.kind]
        if "finite" in t.coeffs:
            labs.append("sec:bl-ST")
        if tuple(t.content) in (("A", "A"), ("G", "G"), ("Z", "A")) and "finite" in t.coeffs:
            labs.append("sec:bl-massless")
        refs += [_ref(x) for x in labs]
        if "pole" in t.coeffs:
            refs.append(TABLES["bl-Z"])
        if "finite" in t.coeffs:
            refs.append(TABLES["bl-ST"])
    out, seen = [], set()
    for r in refs:
        key = (r["label"], r["page"])
        if key not in seen:
            seen.add(key)
            out.append(r)
    return out[:6]


# ---------------------------------------------------------------------------------------------- term json
def _groups(ctx: Ctx, t: BTerm) -> list[str]:
    gs = set()
    for p in t.parent_ids:
        gs |= set(ctx.terms[p].get("groups") or [])
    return [g for g in ("SU(3)", "SU(2)", "U(1)") if g in gs]


def term_json(ctx: Ctx, t: BTerm) -> dict:
    types = [p for p in PARTS if p in t.coeffs]
    differs = _differs(t)
    co = {}
    for p in types:
        entry = {"math": _M(_coef_math(t.coeffs[p], p), interactive=False)}
        if p in ("pole", "evanescent"):
            dec = decomposition(ctx, t, p)
            if dec is not None:
                entry["decomposition"] = _M(dec, interactive=False)
                entry["note"] = _prose(
                    "<p>Sum over the parents of their factor times their coefficient written through the "
                    "renormalisation constants of the unbroken phase (each includes [[lfe]]; "
                    + ("[[Ssum(r)]] and [[Ytil]] as in Section 5.1" if p == "evanescent" else
                       "[[dZc(c)]] stands for [[dpar(c)/c]]") + ").</p>")
        if p == "finite" and differs and set(types) == {"finite"}:
            entry["note"] = _prose("<p>It multiplies the four-dimensional operator shown: all fields, derivatives "
                                   "and metrics are barred (report Section 5.1).</p>")
        elif p == "finite" and differs:
            m4 = mv.M(bf.for_display(mv.operator(t.op4, table=BT)))
            entry["note"] = _prose("<p>It multiplies the four-dimensional form [[" + m4["src"] + "]].</p>") \
                if len(m4["html"]) < 2500 else _prose(
                "<p>It multiplies the four-dimensional form of the operator, in which all fields, derivatives and "
                "metrics are barred.</p>")
        if p == "tree" and t.content and all(x.startswith("b_") for x in t.content):
            entry["note"] = _prose("<p>Feynman gauge: [[eqs(xi, 1)]] in all loop calculations.</p>")
        co[p] = entry
    total = _total(t, types, differs)
    legs = _rule_legs(t.content)
    parts, explicit = {}, set()
    for p in types:
        info = {}
        opp = t.op4 if p == "finite" and t.op4 is not None else t.op
        Vf = mv.vertex(opp, legs, simplifier=bf.display_simplifier(), info=info)
        explicit |= set(info.get("explicit") or ())
        parts[p] = _M(cat._with_loop(bf.to_site(t.coeffs[p]) * Vf, p))
    leg_labels = mv.legs(legs, table=BT, explicit=explicit or None)
    rule = {"title": "Vertex function" if len(legs) > 2 or t.kind in ("source",) else "Two-point vertex function",
            "legs": [{"field": f, "label": _M(lab, interactive=False, display=False)}
                     for f, lab in zip(legs, leg_labels)],
            "parts": parts,
            "note": _prose("<p>The finite part is the vertex function of the four-dimensional operator: its metrics "
                             "and momenta are barred.</p>") if differs else None}
    if len(legs) == 1:
        rule["title"] = "One-point vertex function"
    values = {p: _values_string(t.coeffs[p]) for p in types}
    return {
        "id": t.id,
        "name": t.name,
        "name_html": mv.prose(t.name_src),
        "operator": _M(mv.operator(t.op4 if differs and set(types) == {"finite"} else t.op, table=BT)),
        "fields": list(t.content),
        "legs": len(t.content),
        "ghost_number": sum(BT[x].ghost_number for x in t.content),
        "sector": t.sector,
        "groups": _groups(ctx, t),
        "evanescent": _is_ev(t),
        "types": types,
        "coefficients": co,
        "total": _M(total),
        "rules": [rule],
        "classes": {},
        "theory": t.theory,
        "remarks": remarks(ctx, t),
        "report": references(ctx, t),
        "kind": t.kind,
        "parents": t.parent_ids,
        "data": {"values": values,
                 "parents": [{"id": p, "factor": _values_string(l), "v_power": k} for p, l, k in t.parents],
                 "v_powers": sorted({k for _, _, k in t.parents})},
    }


def _values_string(e: Expression) -> str:
    """A coefficient as a Symbolica string in the plain names of the module docstring."""
    x = bf.to_site(e.expand())
    s = x.format(color_top_level_sum=False, color_builtin_symbols=False, bracket_level_colors=None,
                 max_line_length=None, num_exp_as_superscript=False, multiplication_operator="*",
                 show_namespaces=False)
    return s.replace("𝑖", "I")


def parse_values(text: str) -> Expression:
    """Inverse of :func:`_values_string` (package parameter symbols, ``v``, ``c_W``, ``s_W``, ``sqrt2``)."""
    e = mv.parse(text.replace("I", "𝑖"))
    return from_site(e)


# ---------------------------------------------------------------------------------------------- fields, facets
FIELD_NAMES = {
    "G": "gluon field", "Wp": "W⁺ boson field", "Wm": "W⁻ boson field", "Z": "Z boson field", "A": "photon field",
    "h": "Higgs boson h", "G0": "neutral Goldstone field G⁰", "Gp": "charged Goldstone field G⁺",
    "Gm": "charged Goldstone field G⁻",
    "u": "up quark u (physical components)", "d": "down quark d (physical components)",
    "e": "electron e (physical components)", "nu": "neutrino ν (Dirac field: physical ν_L and sterile P_R ν)",
    "up": "sterile partner u′ of the up quark", "dp": "sterile partner d′ of the down quark",
    "ep": "sterile partner e′ of the electron",
}
_PART_NAME = {"Wp": "W⁺", "Wm": "W⁻", "Z": "Z", "A": "photon", "G": "gluon"}


def fields_json() -> list[dict]:
    kinds = {"gauge": "gauge", "ghost": "ghost", "antighost": "antighost", "auxiliary": "auxiliary",
             "fermion": "fermion", "fermion_bar": "fermion", "scalar": "scalar", "scalar_conj": "scalar",
             "source": "source"}
    out = []
    for name, info in BT.items():
        kind = kinds[info.kind]
        base = name[:-3] if info.kind == "fermion_bar" else name
        if base in ("up", "dp", "ep"):
            kind = "sterile"
        if info.kind in ("fermion", "fermion_bar"):
            fname = FIELD_NAMES[base] if info.kind == "fermion" else f"Dirac conjugate of the {FIELD_NAMES[base]}"
        elif name in FIELD_NAMES:
            fname = FIELD_NAMES[name]
        elif info.kind in ("ghost", "antighost", "auxiliary"):
            f = name.split("_", 1)[1]
            what = {"ghost": "ghost", "antighost": "antighost", "auxiliary": "auxiliary field"}[info.kind]
            fname = f"{what} of the {_PART_NAME[f]}" if f != "G" else f"SU(3) {what}"
        elif name.startswith("src_gh_"):
            f = name[len("src_gh_"):]
            fname = f"BRST source ζ of the {_PART_NAME[f]} ghost component" if f != "G" else \
                "BRST source ζ of the SU(3) ghost"
        elif name.startswith("src_H"):
            i = name[-1]
            fname = f"component {i} of the Higgs source " + ("σ†" if name.startswith("src_H_") else "σ")
        elif name[4] in "QUDLE":
            rest = name[4:]
            f = rest[0]
            comp = rest.split("_")[1] if "_" in rest else None
            what = "R" if "bar" in rest else "R̄"
            fname = f"BRST source {what}_{f}" + (f", SU(2) component {comp}" if comp else "")
        else:
            f = name[4:]
            fname = f"BRST source ρ of the {_PART_NAME[f]}" if f != "G" else "BRST source ρ of the gluon"
        note = ""
        if base in ("u", "d", "e"):
            phys = {"u": ("P_L Ψ¹_Q", "P_R Ψ_U"), "d": ("P_L Ψ²_Q", "P_R Ψ_D"), "e": ("P_L Ψ²_L", "P_R Ψ_E")}[base]
            note = (f"Mass eigenstate built from physical components: {base} = {phys[0]} + {phys[1]} "
                    "(report eq. bl-fermions).")
        elif base in ("up", "dp", "ep"):
            phys = {"up": "P_R Ψ¹_Q + P_L Ψ_U", "dp": "P_R Ψ²_Q + P_L Ψ_D", "ep": "P_R Ψ²_L + P_L Ψ_E"}[base]
            note = f"Built from the sterile components: {FIELD_TEXT[base]} = {phys}; a gauge singlet."
        elif base == "nu":
            note = "ν = Ψ¹_L; P_L ν is physical, P_R ν its sterile partner."
        elif name in ("Wp", "Wm", "Z", "A"):
            note = {"Wp": "W⁺ = (W¹ − i W²)/√2", "Wm": "W⁻ = (W¹ + i W²)/√2", "Z": "Z = c_W W³ − s_W B",
                    "A": "A^γ = s_W W³ + c_W B"}[name] + " (report eq. bl-rot)."
        elif name in ("h", "G0", "Gp", "Gm"):
            note = "Component of φ = H − H_v = (G⁺, (h + i G⁰)/√2) (report eq. bp-H)."
        group = "SU(3)" if (info.factor == "G" or name in ("G",)) else (
            "SU(2)" if info.factor in ("Wp", "Wm") else None)
        out.append({"id": name, "symbol": _M(chip_symbol(name), interactive=False, display=False),
                    "name": fname, "kind": kind, "group": group, "ghost_number": info.ghost_number,
                    "statistics": "fermi" if info.odd else "bose", "note": _eq_plain(note)})
    return out


_CHIP_LABEL = {"G": "G", "Wp": "+", "Wm": "−", "Z": "Z", "A": "γ"}


def chip_symbol(name: str) -> Expression:
    """Symbol of a broken-phase field in the field lists and chips: ghosts, antighosts, auxiliary fields and the
    sources of gauge fields and ghosts carry the field label as a subscript (ω_G, ω_Z, ω̄_+, b_γ, ρ_Z, ζ_−), as in
    the unbroken phase and the graph drawings."""
    for pre, head in (("src_gh_", "zsrc"), ("gh_", "omega"), ("agh_", "omegabar"), ("b_", "b"), ("src_", "rho")):
        rest = name[len(pre):]
        if name.startswith(pre) and rest in _CHIP_LABEL:
            lab = _CHIP_LABEL[rest]
            return nm.field_chip(head, lab if lab in ("G", "Z") else mv.symbol(lab))
    return bf.site_field(name, ())


FACETS = {
    "sector": [("gauge", "Gauge"), ("fermion", "Fermion"), ("scalar", "Higgs"), ("yukawa", "Yukawa"),
               ("gauge-fixing", "Gauge fixing"), ("ghost", "Ghost"), ("source", "BRST source"),
               ("breaking", "Breaking")],
    "type": [("tree", "Tree"), ("pole", "Singular"), ("evanescent", "Evanescent"), ("finite", "Finite restoring"),
             ("breaking", "Breaking")],
    "group": [("SU(3)", "SU(3)"), ("SU(2)", "SU(2)"), ("U(1)", "U(1)")],
    "kind": KINDS,
}


# ---------------------------------------------------------------------------------------------- overview
def _src(e) -> str:
    return mv.M(e)["src"] if isinstance(e, Expression) else e


def overview(ctx: Ctx, counts: dict) -> list[dict]:
    P = mv.parse
    Sb = mv.M(mv.symbol("S", sub="ren", sup="broken"))["src"]
    Phip = mv.M(mv.symbol("Φ", sup="′"))["src"]
    Tk = mv.M(mv.symbol("M", sub="kin"))["src"]
    Sph = mv.M(mv.symbol("S", sub="phys"))["src"]
    up_, ub_ = "ufprime", "ufprimebar"
    out = []
    by_kind = counts["by_kind"]
    out.append({
        "id": "S-broken", "title": "The broken phase as a change of variable",
        "math": _M(P(f"eqs({Sb}(phi, {Phip}), Sren(Hv + phi, {Phip}))"), interactive=False),
        "text": _prose(
            "<p>The renormalised action of the broken phase is the renormalised action of the unbroken phase, "
            "[[eqs(Sren, S0 + deltaS1 + DeltaS1ST)]], written in the variable [[eqs(phi, H - Hv)]], with "
            f"[[{Phip}]] all other fields and sources (report eq. bl-def). This is exact: the shift by the constant "
            "doublet [[Hv]] is a change of variable, the counterterms are local polynomials in [[H]], and the "
            "linear gauge fixing does not depend on [[H]], so the Slavnov–Taylor operator keeps its form and "
            "[[eqs(ord(sop, phi), 𝑖*ord(wm, par(Hv + phi)))]]. No loop calculation beyond those of the unbroken "
            "phase enters, and no input from an [[field(R, xi)]] gauge is used.</p>")
        + _prose(
            f"<p>Every term of the table below was obtained by applying this substitution to the operators of the "
            f"unbroken catalogue and collecting equal operators; its coefficients are sums of the parents' "
            f"coefficients times exact factors. The table has {counts['terms']} terms: "
            + ", ".join(f"{by_kind.get(k, 0)} of kind “{lab.lower()}”" for k, lab in KINDS)
            + ". Field-independent constants (the value of the potential at [[Hv]]) are omitted.</p>")})
    out.append({
        "id": "shift", "title": "Shift and vacuum value",
        "math": _M(P("eqs(H, Hv + phi)"), interactive=False),
        "text": _prose(
            "<p>[[Hv]] has the upper component 0 and the lower component [[v/sqrtn(2)]]; the fluctuation has the "
            "components [[eqs(field(phi, 1), Gp)]] and [[eqs(field(phi, 2), (h + 𝑖*G0)/sqrtn(2))]], with [[h]], "
            "[[G0]] real and [[eqs(Gm, cc(Gp))]] (report eq. bp-H). The vacuum value is defined at the minimum of "
            "the tree-level potential written with the MS-bar parameters, [[eqs(m2, -lambda*v^2)]] with [[m2]] "
            "negative (eq. bl-vev): [[v]] is a function of [[m2]] and [[lambda]], is not renormalised, and the "
            "tadpole graphs are kept. Every coefficient below uses this condition.</p>")})
    out.append({
        "id": "rotation", "title": "Field basis: rotation of the gauge sector",
        "math": _M(P("eqs(Z(mu), minus(c_W*W(mu, 3), s_W*B(mu)))"), interactive=False),
        "text": _prose(
            "<p>With [[eqs(c_W, g/sqrtn(g^2 + gp^2))]], [[eqs(s_W, gp/sqrtn(g^2 + gp^2))]] and "
            "[[eqs(e, eqs(g*s_W, gp*c_W))]] (report eq. bp-rot): [[eqs(Agam(mu), plus(s_W*W(mu, 3), c_W*B(mu)))]], "
            "[[eqs(Wp(mu), par(minus(W(mu, 1), 𝑖*W(mu, 2)))/sqrtn(2))]] and [[eqs(Wm(mu), par(plus(W(mu, 1), 𝑖*W(mu, 2)))/sqrtn(2))]] "
            "(eq. bl-rot). The same rotation, with the renormalised angle, is applied to the ghosts, antighosts, "
            "auxiliary fields and the sources [[rho(mu, a)]] and [[zsrc(a)]]. The report does not rotate the "
            "fermion sources and the Higgs source; they appear in their SU(2) components, for example "
            "[[Rbar(Q, up(1))]] and [[sigmadag(2)]]. The coefficients below are written with [[g]], [[gp]], "
            "[[c_W]], [[s_W]] and [[e]] in whichever equivalent form is shortest.</p>")})
    out.append({
        "id": "fermions", "title": "Fermion mass eigenstates and sterile partners",
        "math": _M(P("eqs(uf, plus(ord(PL, Q(1)), ord(PR, U)))"), interactive=False),
        "text": _prose(
            "<p>The mass eigenstates are Dirac fields built from the physical components, with sterile partners "
            "built from the sterile ones (report eq. bl-fermions): [[eqs(ufprime, plus(ord(PR, Q(1)), ord(PL, U)))]], "
            "[[eqs(dqf, plus(ord(PL, Q(2)), ord(PR, D)))]], [[eqs(dqfprime, plus(ord(PR, Q(2)), ord(PL, D)))]], "
            "[[eqs(elf, plus(ord(PL, L(2)), ord(PR, E)))]], [[eqs(elfprime, plus(ord(PR, L(2)), ord(PL, E)))]] and "
            "[[eqs(nuf, L(1))]], where [[Q(1)]], [[Q(2)]] are the upper and lower SU(2) components of the doublet. "
            "Mass terms, Yukawa and gauge couplings contain physical components only; the sterile partners appear "
            "only in the kinetic terms, through the evanescent terms [[𝑖*ord(ufbar, gammahat(mu), pd(mu, ufprime))]] "
            "and [[𝑖*ord(ufprimebar, gammahat(mu), pd(mu, uf))]].</p>")})
    mass = cat._seq("M_W^2*ord(Wp(mu), Wm(up(mu)))", "1/2*M_Z^2*ord(Z(mu), Z(up(mu)))", "-1/2*M_h^2*ord(h, h)",
                    "-sumover(f, m_f*ord(bar(field(f)), field(f)))")
    out.append({
        "id": "masses", "title": "Tree-level masses",
        "math": _M(P(f"eqs(field(ℒ, mass), {mass})"), interactive=False),
        "text": _prose(
            "<p>The gauge-boson mass terms come from [[ord(dagger(par(ord(Dcov(mu), Hv))), Dcov(up(mu)), Hv)]], "
            "a mass matrix [[eqs(field(ℳ, a*b)^2, ord(dagger(Hv), acomm(Tg(a, H), Tg(b, H)), Hv))]] over the four "
            "SU(2)×U(1) components, with eigenvalues [[M_W^2]] (twice), [[M_Z^2]] and 0 (photon) (report "
            "eq. bp-masses). The tree-level masses are [[eqs(M_W, g*v/2)]], [[eqs(M_Z, M_W/c_W)]], "
            "[[eqs(M_h^2, 2*lambda*v^2)]] and [[eqs(field(m, k), y_k*v/sqrtn(2))]] for k = u, d, e (eq. bl-masses); "
            "the photon, the gluon, the neutrino, the Goldstone fields and all ghosts are massless in this gauge. "
            "They are MS-bar masses at the tree-level minimum, not pole masses. In the table the corresponding "
            "terms are written with [[g]], [[v]], [[y_k]], [[lambda]]; their remarks give the identification.</p>")})
    mix = "plus(M_Z*ord(Z(mu), pd(up(mu), G0)), 𝑖*M_W*par(minus(ord(Wm(mu), pd(up(mu), Gp)), ord(Wp(mu), pd(up(mu), Gm)))))"
    out.append({
        "id": "mixing", "title": "Gauge–Goldstone mixing",
        "math": _M(P(f"eqs(field(ℒ, mix), {mix})"), interactive=False),
        "text": _prose(
            "<p>Because the gauge fixing is linear and does not depend on [[H]], the tree-level mixing between "
            "gauge and Goldstone fields stays in the quadratic action (report eq. bl-L2mix); it enters the "
            "propagators, which have mixed gauge–Goldstone parts and double poles at [[eqs(sq(k), 0)]] "
            "(eq. bl-prop). The singular counterterm adds a photon–Goldstone term, and the finite counterterm "
            "changes the mixing by [[-Ytil/3]] times its tree value (eq. bl-mix).</p>")})
    out.append({
        "id": "tadpole", "title": "Tadpole",
        "math": _M(P("eqs(dpar(T), plus(plus(dpar(m2), v^2*dpar(lambda)), lambda*v^2*dZ(H)))"), interactive=False),
        "text": _prose(
            "<p>The tree-level tadpole [[-v*par(m2 + lambda*v^2)*h]] vanishes by the vacuum condition. The singular "
            "counterterm contains the tadpole term [[-v*dpar(T)*h]] (report eq. bl-dS-scalar, eq. bl-dT), which "
            "cancels the pole of the one-loop tadpole; here "
            "[[eqs(v*dpar(T), lfe*v^3*par(6*lambda^2 - lambda*Cinv(H) - field(Y, 4) + 3/16*field(G, 4)))]] with "
            "[[eqs(Cinv(H), (2*M_W^2 + M_Z^2)/v^2)]], [[eqs(field(Y, 4), 4*(3*m_u^4 + 3*m_d^4 + m_e^4)/v^4)]] and "
            "[[eqs(field(G, 4), 16*(2*M_W^4 + M_Z^4)/v^4)]] (eq. bl-sums). The Goldstone mass counterterm equals "
            "[[dpar(T)]]. Only the [[lambda^2]] part of the tadpole pole was checked in the report; the finite "
            "one-loop tadpole is not computed there. The finite restoring counterterm gives no tadpole.</p>")})
    out.append({
        "id": "counterterms", "title": "How the counterterms follow",
        "math": _M(P("eqs(field(H, 0), sqrtn(field(Z, H))*par(Hv + phi))"), interactive=False),
        "text": _prose(
            "<p>The singular counterterm is the first-order term of [[S0]] in the bare quantities, now at "
            "[[eqs(H, Hv + phi)]]: [[v]] is not renormalised, the rotation uses the renormalised [[g]] and [[gp]], "
            "and the gauge fixing is not renormalised (report Section 5.4.4). The v-dependent mass, mixing and "
            "tadpole counterterms therefore follow from the unbroken-phase constants of Table 57. The finite "
            "restoring counterterm is [[DeltaS1ST]] at [[eqs(H, Hv + phi)]] (Section 5.4.5): its H-dependent "
            "coefficients are purely Yukawa, so after the shift its terms with [[Hv]] are proportional to fermion "
            "masses squared; no coefficient is recomputed. Under (N1) there is no finite counterterm for the "
            "fermion masses, the Yukawa vertices, the scalar potential or the tadpole.</p>")})
    prop = f"eqs({Sph}, (sq(p)*slashbar(p) + m_f*sqbar(p))/(par(sq(p))^2 - m_f^2*sqbar(p)))"
    out.append({
        "id": "propagator", "title": "Fermion propagator with sterile partners",
        "math": _M(P(prop), interactive=False),
        "text": _prose(
            f"<p>For the pair [[uf]], [[{up_}]] the quadratic Lagrangian has the kinetic matrix with the entries "
            "[[slashbar(p) - m_u]] and [[slashbar(p)]] on the diagonal and [[slashhat(p)]] off the diagonal (report "
            f"eq. bp-prop); the propagator is [[𝑖]] times its inverse, whose physical block [[{Sph}]] is shown "
            "(with [[eqs(m_f, m_u)]]). Its denominator depends on [[sq(p)]] and [[sqbar(p)]] separately and is not "
            "d-dimensionally covariant; at vanishing hatted momentum, where external fermions live, it reduces to "
            "the standard [[(slashbar(p) + m_f)/(sqbar(p) - m_f^2)]]. The other blocks are "
            "[[(sq(p)*slashhat(p) + m_f*ord(slashbar(p), slashhat(p)))/(par(sq(p))^2 - m_f^2*sqbar(p))]], "
            "[[(sq(p)*slashhat(p) - m_f*ord(slashbar(p), slashhat(p)))/(par(sq(p))^2 - m_f^2*sqbar(p))]] and "
            "[[((sq(p) - m_f^2)*slashbar(p) + m_f*sqhat(p))/(par(sq(p))^2 - m_f^2*sqbar(p))]]. The neutrino has "
            "the massless propagator. The counterterms are unaffected: the mass can be treated as a two-point "
            "insertion on massless propagators, so poles and the finite breaking are polynomials in the "
            "masses (Section 5.5.2).</p>")})
    return out


def intro(ctx: Ctx, counts: dict) -> str:
    from . import theory

    base = theory.phase_intro("broken")
    return base + _prose(
        f"<p>The {counts['terms']} terms are derived, not typed: the substitution [[to(H, Hv + phi)]] and the "
        "rotation to the mass basis are applied to each operator of the unbroken catalogue (every SU(2) index in "
        "components), equal operators are collected, and each term lists its parents in the unbroken phase with "
        "the factor and the power of [[v]] with which they contribute. Terms that do not involve the Higgs doublet "
        "have the kind “unchanged by the shift”; they are the unbroken terms written in the mass basis.</p>")


# ---------------------------------------------------------------------------------------------- build
def catalogue(ctx: Ctx, log=print) -> list[BTerm]:
    """The derived broken-phase catalogue (module docstring), in display order."""
    for pid in ctx.order:
        _EVANESCENT[pid] = ctx.specs[pid].evanescent
    _scalar_names(ctx.model)
    images: dict[tuple, list[Image]] = {}
    for pid in ctx.order:
        for c, im in parent_images(ctx, pid).items():
            if not c:
                continue          # field-independent constant
            images.setdefault(c, []).append(im)
    log(f"  {sum(len(v) for v in images.values())} images on {len(images)} field contents")
    fp, _ = fingerprinter()
    terms = collect(ctx, images, fp)
    for t in terms:
        _coefficients(ctx, t)
    terms = [t for t in terms if t.coeffs]
    terms = _merge(terms)
    for t in terms:
        finalise(t)
        classify(ctx, t)
    order = {p: i for i, p in enumerate(ctx.order)}
    sec = {s: i for i, s in enumerate(cat.SECTOR_ORDER)}
    kind = {k: i for i, (k, _) in enumerate(KINDS)}
    terms.sort(key=lambda t: (sec[t.sector], kind[t.kind] if t.kind != "unchanged" else -1,
                              min(order[p] for p in t.parent_ids), len(t.content),
                              [BT.order(x) for x in t.content]))
    name_terms(ctx, terms)
    return terms


def check_independence(terms: list[BTerm], fp: Fingerprinter | None = None) -> int:
    """The operators of each field content are linearly independent (exact fingerprints, D = 6); returns the
    number of contents checked."""
    fp = fp or fingerprinter()[0]
    by: dict[tuple, list[BTerm]] = {}
    for t in terms:
        by.setdefault(tuple(t.content), []).append(t)
    n = 0
    for c, ts in by.items():
        if len(ts) < 2:
            continue
        if _fp_rank(fp, c, [t.op for t in ts]) != len(ts):
            raise RuntimeError(f"the operators of {c} are linearly dependent: {[t.id for t in ts]}")
        n += 1
    return n


def build_broken(site_dir: str | Path = SITE, *, log=print) -> dict:
    """Write ``data/broken.json``; returns counts."""
    t0 = time.time()
    ctx = Ctx(site_dir, log=log)
    log(f"{len(ctx.order)} unbroken terms; shifting and collecting")
    terms = catalogue(ctx, log=log)
    nind = check_independence(terms)
    log(f"  {len(terms)} broken-phase terms; independence checked on {nind} contents ({time.time() - t0:.0f} s)")
    out_terms = [term_json(ctx, t) for t in terms]
    by_type, by_kind, by_sector = {}, {}, {}
    for t in out_terms:
        for ty in t["types"]:
            by_type[ty] = by_type.get(ty, 0) + 1
        by_kind[t["kind"]] = by_kind.get(t["kind"], 0) + 1
        by_sector[t["sector"]] = by_sector.get(t["sector"], 0) + 1
    counts = {"terms": len(out_terms), "by_type": by_type, "by_kind": by_kind, "by_sector": by_sector,
              "contents": len({tuple(t["fields"]) for t in out_terms}), "independence_checked": nind}
    doc = {
        "phase": "broken",
        "title": "Broken phase",
        "intro": intro(ctx, counts),
        "overview": overview(ctx, counts),
        "fields": fields_json(),
        "facets": {k: [{"id": i, "label": lab} for i, lab in v] for k, v in FACETS.items()},
        "terms": out_terms,
        "counts": counts,
    }
    bad = cat._omissions(doc)
    if bad:
        raise RuntimeError(f"expressions with omitted parts: {bad[:10]}")
    out = Path(site_dir) / "data" / "broken.json"
    text = json.dumps(doc, ensure_ascii=False, separators=(",", ":"))
    out.write_text(text)
    log(f"wrote {out.name}: {len(out_terms)} terms, {len(text) / 1e6:.2f} MB; kinds {by_kind}; types {by_type}; "
        f"sectors {by_sector} ({time.time() - t0:.0f} s)")
    return counts


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m gdocs.broken", description=__doc__.split("\n\n")[0])
    ap.add_argument("--site", default=str(SITE), help="site directory (default: gamma5_docs/site)")
    args = ap.parse_args(argv)
    build_broken(args.site)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
