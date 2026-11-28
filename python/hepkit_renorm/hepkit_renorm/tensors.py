"""The tensor vocabulary of coefficients, index labels and the BMHV simplifier.

Every coefficient in the package is a Symbolica expression built from the heads defined here (namespace
``hepkit_renorm``). All indices are explicit arguments, so products of heads commute; Grassmann-odd objects never
appear in coefficients (they live in the ordered field occurrences of :mod:`hepkit_renorm.functional`).

Index labels (one head per index space; the space is part of the label)::

    lor(n)        Lorentz index, d-dimensional (abstract; placement is irrelevant: a label that appears twice
                  is contracted with the metric g^{mu nu} = diag(+,-,-,-,-,...))
    spin(n)       Dirac spinor index
    adj(F, n)     adjoint index of the non-abelian gauge factor F
    fund(F, n)    fundamental index of the SU(N) factor F (fields in the antifundamental use the same space)
    flav(n)       flavour (generation) index

Heads (report conventions, sections 03_setup and 03b)::

    g(mu,nu), gbar(mu,nu), ghat(mu,nu)   full, four-dimensional and (d-4)-dimensional metric (projectors);
                                          g.g = g, gbar.gbar = gbar, ghat.ghat = ghat, gbar.ghat = 0,
                                          g(mu,mu) = d, gbar(mu,mu) = 4, ghat(mu,mu) = d - 4
    eps(mu,nu,rho,sigma)                  four-dimensional Levi-Civita tensor, eps^{0123} = +1 (report; Idenso's
                                          symbol is -i times this one); eps . ghat = 0
    p(k, mu)                              incoming momentum of leg k (vertex functions only)
    sp(k,l), spbar(k,l), sphat(k,l)       p_k.p_l with the full, barred and hatted metric
    gamma(mu,s,t)                         (gamma^mu)_{st}, d-dimensional
    gammabar(mu,s,t), gammahat(mu,s,t)    its four-dimensional and (d-4)-dimensional parts, gamma = gammabar +
                                          gammahat; gamma(nu,..) gbar(nu,mu) = gammabar(mu,..), gammabar.ghat = 0
    gamma5(s,t), PL(s,t), PR(s,t)         gamma5 = i g^0 g^1 g^2 g^3, PL = (1 - gamma5)/2, PR = (1 + gamma5)/2;
                                          {gammabar^mu, gamma5} = 0, [gammahat^mu, gamma5] = 0
    delta(i,j)                            Kronecker delta of any space (both labels of the same space)
    T(a,i,j)                              generator (T^a)_{ij} of an SU(N) factor in the fundamental,
                                          Tr(T^a T^b) = delta^{ab}/2, [T^a,T^b] = i f^{abc} T^c
    f(a,b,c), dsym(a,b,c)                 structure constants and d^{abc} = 2 Tr(T^a {T^b,T^c})
    epsf(i,j,...)                         Levi-Civita symbol of N fundamental indices, epsf(1,2,..,N) = +1
    d                                     the dimension d = 4 - 2 epsilon (a real symbol)

The coupling-weighted objects of the report (cal(T)^a = g_r T^a, fg^{abc} = g_r f^{abc}) are products of a coupling
symbol with these heads. The factor of a group head is carried by its labels.

Conjugation (:func:`conjugate`) is complex conjugation of scalars combined with the Dirac adjoint on spinor
heads with transposed spinor labels: ``gamma(mu,s,t) -> gamma(mu,t,s)`` (also for ``gammabar``, ``gammahat``),
``gamma5(s,t) -> -gamma5(t,s)``,
``PL(s,t) -> PR(t,s)``, and ``T(a,i,j) -> T(a,j,i)`` (hermitian generators). Together with the field map
``Psi -> Psibar`` this implements ``(Psibar M Chi)^dagger = Chibar (gamma0 M^dagger gamma0) Psi``.
"""

from __future__ import annotations

import itertools
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass

from symbolica import AtomType, Expression, Replacement, S

__all__ = [
    "FLAVOUR",
    "HAT_HEADS",
    "HEADS",
    "LORENTZ",
    "NS",
    "SPINOR",
    "I",
    "LabelFactory",
    "Simplifier",
    "Space",
    "absorb",
    "adjoint_space",
    "canonical_form",
    "canonical_relabel",
    "conjugate",
    "contains_head",
    "factor_symbol",
    "free_labels",
    "fresh",
    "fundamental_space",
    "is_evanescent",
    "label_space",
    "labels_in",
    "rename",
]

NS = "hepkit_renorm"
I = Expression.I

# ---------------------------------------------------------------------------------------------- symbols
lor = S(f"{NS}::lor", is_real=True)
spin = S(f"{NS}::spin", is_real=True)
adj = S(f"{NS}::adj", is_real=True)
fund = S(f"{NS}::fund", is_real=True)
flav = S(f"{NS}::flav", is_real=True)

g = S(f"{NS}::g", is_symmetric=True, is_real=True)
gbar = S(f"{NS}::gbar", is_symmetric=True, is_real=True)
ghat = S(f"{NS}::ghat", is_symmetric=True, is_real=True)
eps = S(f"{NS}::eps", is_antisymmetric=True, is_real=True)
p = S(f"{NS}::p", is_real=True)
sp = S(f"{NS}::sp", is_symmetric=True, is_real=True)
spbar = S(f"{NS}::spbar", is_symmetric=True, is_real=True)
sphat = S(f"{NS}::sphat", is_symmetric=True, is_real=True)
gamma = S(f"{NS}::gamma", is_real=True)
gammabar = S(f"{NS}::gammabar", is_real=True)
gammahat = S(f"{NS}::gammahat", is_real=True)
gamma5 = S(f"{NS}::gamma5", is_real=True)
PL = S(f"{NS}::PL", is_real=True)
PR = S(f"{NS}::PR", is_real=True)
delta = S(f"{NS}::delta", is_symmetric=True, is_real=True)
T = S(f"{NS}::T", is_real=True)
f = S(f"{NS}::f", is_antisymmetric=True, is_real=True)
dsym = S(f"{NS}::dsym", is_symmetric=True, is_real=True)
epsf = S(f"{NS}::epsf", is_antisymmetric=True, is_real=True)
d = S(f"{NS}::d", is_real=True)

HEADS = {
    "lor": lor, "spin": spin, "adj": adj, "fund": fund, "flav": flav,
    "g": g, "gbar": gbar, "ghat": ghat, "eps": eps, "p": p, "sp": sp, "spbar": spbar, "sphat": sphat,
    "gamma": gamma, "gammabar": gammabar, "gammahat": gammahat, "gamma5": gamma5, "PL": PL, "PR": PR, "delta": delta,
    "T": T, "f": f, "dsym": dsym, "epsf": epsf, "d": d,
}
"""Name -> Symbolica symbol of every head of the vocabulary."""

_LABEL_HEADS = {lor.get_name(), spin.get_name(), adj.get_name(), fund.get_name(), flav.get_name()}

_FACTORS: dict[str, Expression] = {}


def factor_symbol(name: str) -> Expression:
    """The symbol that names the gauge factor ``name`` inside ``adj`` / ``fund`` labels."""
    if name not in _FACTORS:
        _FACTORS[name] = S(f"{NS}_factor::{name}", is_real=True)
    return _FACTORS[name]


# ---------------------------------------------------------------------------------------------- spaces
@dataclass(frozen=True)
class Space:
    """An index space. ``kind`` is one of ``lorentz``, ``spinor``, ``adjoint``, ``fund``, ``flavour``;
    ``factor`` names the gauge factor for ``adjoint`` / ``fund``; ``dim`` is the dimension (``None`` for the
    d-dimensional Lorentz space)."""

    kind: str
    factor: str | None = None
    dim: int | None = None

    def label(self, n: int) -> Expression:
        """The label number ``n`` of this space."""
        if self.kind == "lorentz":
            return lor(n)
        if self.kind == "spinor":
            return spin(n)
        if self.kind == "adjoint":
            return adj(factor_symbol(self.factor), n)
        if self.kind == "fund":
            return fund(factor_symbol(self.factor), n)
        if self.kind == "flavour":
            return flav(n)
        raise ValueError(self.kind)

    @property
    def key(self) -> str:
        """A string that identifies the space (used as canonicalisation group marker)."""
        return self.kind if self.factor is None else f"{self.kind}:{self.factor}"


LORENTZ = Space("lorentz")
SPINOR = Space("spinor", dim=4)


def FLAVOUR(n: int) -> Space:
    """The flavour space with ``n`` generations."""
    return Space("flavour", dim=n)


def adjoint_space(factor: str, dim: int) -> Space:
    return Space("adjoint", factor, dim)


def fundamental_space(factor: str, dim: int) -> Space:
    return Space("fund", factor, dim)


class LabelFactory:
    """Fresh labels. Numbers start far above the canonical pool (0, 1, ...) so that fresh labels never collide
    with canonical ones. Deterministic: the n-th call in a process always returns the same label."""

    def __init__(self, start: int = 1_000_000) -> None:
        self._next = start

    def __call__(self, space: Space) -> Expression:
        n = self._next
        self._next += 1
        return space.label(n)


fresh = LabelFactory()
"""The process-wide label factory used by all modules."""


def label_space(label: Expression) -> str:
    """Group key of a label: ``lorentz``, ``spinor``, ``adjoint:F``, ``fund:F`` or ``flavour``."""
    name = label.get_name()
    if name == lor.get_name():
        return "lorentz"
    if name == spin.get_name():
        return "spinor"
    if name == flav.get_name():
        return "flavour"
    fac = str(next(iter(label))).split("::")[-1]
    return ("adjoint:" if name == adj.get_name() else "fund:") + fac


def _is_label(e: Expression) -> bool:
    return e.get_type() == AtomType.Fn and e.get_name() in _LABEL_HEADS


def _walk_labels(e: Expression) -> Iterator[Expression]:
    t = e.get_type()
    if t == AtomType.Fn:
        if e.get_name() in _LABEL_HEADS:
            yield e
            return
        for a in e:
            yield from _walk_labels(a)
    elif t in (AtomType.Mul, AtomType.Add):
        for a in e:
            yield from _walk_labels(a)
    elif t == AtomType.Pow:
        base, ex = list(e)
        n = int(str(ex)) if ex.get_type() == AtomType.Num else 1
        for _ in range(max(n, 1)):
            yield from _walk_labels(base)


HAT_HEADS = (ghat, gammahat, sphat)
"""Heads that carry a (d-4)-dimensional index or product: a term containing one of them is evanescent."""


def is_evanescent(term: Expression) -> bool:
    """True if a (simplified, split) term contains a hatted object (:data:`HAT_HEADS`)."""
    return any(contains_head(term, h) for h in HAT_HEADS)


def contains_head(expr: Expression, head: Expression) -> bool:
    """True if a function with the symbol ``head`` occurs anywhere in ``expr``."""
    name = head.get_name()
    stack = [expr]
    while stack:
        e = stack.pop()
        t = e.get_type()
        if t == AtomType.Fn:
            if e.get_name() == name:
                return True
            stack.extend(e)
        elif t in (AtomType.Mul, AtomType.Add, AtomType.Pow):
            stack.extend(e)
    return False


def labels_in(term: Expression) -> list[Expression]:
    """All label occurrences in a product term, with multiplicity (``x^2`` counts twice), in a deterministic
    order (Symbolica's canonical argument order)."""
    return list(_walk_labels(term))


def free_labels(term: Expression) -> list[Expression]:
    """Labels that occur exactly once in a product term (its free indices), in first-occurrence order."""
    counts: dict[str, list] = {}
    for lab in labels_in(term):
        k = lab.to_canonical_string()
        counts.setdefault(k, [lab, 0])[1] += 1
    return [v[0] for v in counts.values() if v[1] == 1]


def rename(expr: Expression, mapping: Mapping[str, Expression] | Iterable[tuple[Expression, Expression]]) -> Expression:
    """Simultaneous renaming of labels: ``mapping`` maps old labels (or their canonical strings) to new ones."""
    pairs = mapping.items() if isinstance(mapping, Mapping) else mapping
    reps = []
    for old, new in pairs:
        old_e = Expression.parse(old) if isinstance(old, str) else old
        if old_e.to_canonical_string() != new.to_canonical_string():
            reps.append(Replacement(old_e, new))
    if not reps:
        return expr
    return expr.replace_multiple(reps)


# ---------------------------------------------------------------------------------------------- conjugation
_s, _t, _m, _a = S("s_", "t_", "m_", "a_")
_CONJ_RULES = [
    Replacement(gamma(_m, _s, _t), gamma(_m, _t, _s)),
    Replacement(gammabar(_m, _s, _t), gammabar(_m, _t, _s)),
    Replacement(gammahat(_m, _s, _t), gammahat(_m, _t, _s)),
    Replacement(gamma5(_s, _t), -gamma5(_t, _s)),
    Replacement(PL(_s, _t), PR(_t, _s)),
    Replacement(PR(_s, _t), PL(_t, _s)),
    Replacement(T(_a, _s, _t), T(_a, _t, _s)),
]


def conjugate(expr: Expression) -> Expression:
    """Complex conjugate of a coefficient, with the Dirac adjoint on spinor heads (see the module docstring).

    Real symbols (``is_real=True``) and the real heads of the vocabulary are unchanged; other symbols and
    parameter tensors become ``conj(...)``.
    """
    return expr.replace_multiple(_CONJ_RULES).conj()


# ---------------------------------------------------------------------------------------------- simplifier
_w = {n: S(n + "_") for n in ("a", "b", "c", "k", "l", "m", "n", "F")}
_x, _y, _fn = S("x___", "y___", "f_")


class Simplifier:
    """Rewriting rules of the BMHV tensor algebra (report 03_setup sec:bmhv), applied to fixpoint.

    Parameters
    ----------
    group_dims:
        ``{factor: (dim_fundamental, dim_adjoint)}`` for the non-abelian factors (traces of ``delta``).
    n_flavours:
        dimension of the flavour space.
    expand_projectors:
        canonical mode (default). ``PL``/``PR`` are written as ``(1 -/+ gamma5)/2``, every full ``gamma`` is split
        into ``gammabar + gammahat``, every ``gamma5`` is moved to the right end of its chain, and the remaining full
        metrics between free labels and the full scalar products ``sp`` are split into barred and hatted parts.
        Open chains then read ``Gamma`` or ``Gamma gamma5`` with ``Gamma`` a product of barred and hatted gammas, so
        equal expressions written with the same order of gamma matrices (and both unsplit or both split by
        :meth:`Density.split_lorentz` in their field contractions) have the same form. If false, projectors are kept
        and moved to the right, full gammas are kept (readable form, not canonical).

    Rules in both modes: full metrics and Kronecker deltas are absorbed into any head or field slot;
    ``gbar``/``ghat`` projector algebra and traces (``gbar(mu,mu) = 4``, ``ghat(mu,mu) = d - 4``); ``eps`` is
    four-dimensional; contracted momenta become ``sp``/``spbar``/``sphat``; ``gamma5`` anticommutes with
    ``gammabar`` and commutes with ``gammahat``; ``gamma5 gamma^mu = -gamma^mu gamma5 + 2 gammahat^mu gamma5``. Gamma
    matrices are never reordered among themselves (general equality: use exact numerical fingerprints).
    """

    def __init__(self, group_dims: Mapping[str, tuple[int, int]] | None = None, n_flavours: int = 1,
                 expand_projectors: bool = True) -> None:
        self.expand_projectors = expand_projectors
        a, b, c, k, l_, m = (_w[n] for n in "abcklm")
        n_ = _w["n"]
        x, y, fn = _x, _y, _fn
        zero, four = Expression.num(0), Expression.num(4)
        rules = [
            Replacement(g(a, a), d),
            Replacement(g(a, b) ** 2, d),
            Replacement(g(a, b) * fn(x, b, y), fn(x, a, y)),
            Replacement(delta(spin(n_), spin(n_)), four),
            Replacement(delta(spin(n_), spin(m)) ** 2, four),
            Replacement(delta(flav(n_), flav(n_)), Expression.num(n_flavours)),
            Replacement(delta(flav(n_), flav(m)) ** 2, Expression.num(n_flavours)),
        ]
        for fac, (nf, na) in (group_dims or {}).items():
            F = factor_symbol(fac)
            rules += [
                Replacement(delta(fund(F, n_), fund(F, n_)), Expression.num(nf)),
                Replacement(delta(fund(F, n_), fund(F, m)) ** 2, Expression.num(nf)),
                Replacement(delta(adj(F, n_), adj(F, n_)), Expression.num(na)),
                Replacement(delta(adj(F, n_), adj(F, m)) ** 2, Expression.num(na)),
            ]
        rules += [
            Replacement(delta(a, b) * fn(x, b, y), fn(x, a, y)),
            # projector metrics
            Replacement(gbar(a, a), four),
            Replacement(ghat(a, a), d - 4),
            Replacement(gbar(a, b) ** 2, four),
            Replacement(ghat(a, b) ** 2, d - 4),
            Replacement(gbar(a, b) * ghat(b, c), zero),
            Replacement(gbar(a, b) * ghat(a, b), zero),
            Replacement(gbar(a, b) * gbar(b, c), gbar(a, c)),
            Replacement(ghat(a, b) * ghat(b, c), ghat(a, c)),
            Replacement(eps(x, b, y) * gbar(b, c), eps(x, c, y)),
            Replacement(eps(x, b, y) * ghat(b, c), zero),
            # barred and hatted gammas
            Replacement(gamma(b, k, l_) * gbar(b, c), gammabar(c, k, l_)),
            Replacement(gamma(b, k, l_) * ghat(b, c), gammahat(c, k, l_)),
            Replacement(gammabar(b, k, l_) * gbar(b, c), gammabar(c, k, l_)),
            Replacement(gammahat(b, k, l_) * ghat(b, c), gammahat(c, k, l_)),
            Replacement(gammabar(b, k, l_) * ghat(b, c), zero),
            Replacement(gammahat(b, k, l_) * gbar(b, c), zero),
            Replacement(gammahat(b, k, l_) * eps(x, b, y), zero),
            # momenta
            Replacement(p(k, a) ** 2, sp(k, k)),
            Replacement(p(k, a) * p(l_, a), sp(k, l_)),
            Replacement(p(k, a) * gbar(a, b) * p(l_, b), spbar(k, l_)),
            Replacement(p(k, a) * ghat(a, b) * p(l_, b), sphat(k, l_)),
            # traces of a single matrix
            Replacement(gamma5(a, a), zero),
            Replacement(gamma(m, a, a), zero),
            Replacement(gammabar(m, a, a), zero),
            Replacement(gammahat(m, a, a), zero),
            # gamma5 algebra
            Replacement(gamma5(a, b) * gamma5(b, c), delta(a, c)),
            Replacement(gamma5(a, b) * gammabar(m, b, c), -gammabar(m, a, b) * gamma5(b, c)),
            Replacement(gamma5(a, b) * gammahat(m, b, c), gammahat(m, a, b) * gamma5(b, c)),
            Replacement(gamma5(a, b) * gamma(m, b, c),
                        -gamma(m, a, b) * gamma5(b, c) + 2 * gammahat(m, a, b) * gamma5(b, c)),
        ]
        if expand_projectors:
            rules += [
                Replacement(PL(a, b), (delta(a, b) - gamma5(a, b)) / 2),
                Replacement(PR(a, b), (delta(a, b) + gamma5(a, b)) / 2),
                Replacement(gamma(m, a, b), gammabar(m, a, b) + gammahat(m, a, b)),
            ]
            self._split = [
                Replacement(g(a, b), gbar(a, b) + ghat(a, b)),
                Replacement(sp(k, l_), spbar(k, l_) + sphat(k, l_)),
            ]
        else:
            rules += [
                Replacement(PL(a, a), Expression.num(2)),
                Replacement(PR(a, a), Expression.num(2)),
                Replacement(PL(a, b) * PL(b, c), PL(a, c)),
                Replacement(PR(a, b) * PR(b, c), PR(a, c)),
                Replacement(PL(a, b) * PR(b, c), zero),
                Replacement(PR(a, b) * PL(b, c), zero),
                Replacement(gamma5(a, b) * PL(b, c), -PL(a, c)),
                Replacement(gamma5(a, b) * PR(b, c), PR(a, c)),
                Replacement(PL(a, b) * gamma5(b, c), -PL(a, c)),
                Replacement(PR(a, b) * gamma5(b, c), PR(a, c)),
                Replacement(PL(a, b) * gammabar(m, b, c), gammabar(m, a, b) * PR(b, c)),
                Replacement(PR(a, b) * gammabar(m, b, c), gammabar(m, a, b) * PL(b, c)),
                Replacement(PL(a, b) * gammahat(m, b, c), gammahat(m, a, b) * PL(b, c)),
                Replacement(PR(a, b) * gammahat(m, b, c), gammahat(m, a, b) * PR(b, c)),
                Replacement(PL(a, b) * gamma(m, b, c), gammabar(m, a, b) * PR(b, c) + gammahat(m, a, b) * PL(b, c)),
                Replacement(PR(a, b) * gamma(m, b, c), gammabar(m, a, b) * PL(b, c) + gammahat(m, a, b) * PR(b, c)),
            ]
            self._split = []
        self._rules = rules

    def __call__(self, expr: Expression, max_passes: int = 200) -> Expression:
        """Expand, rewrite to fixpoint and relabel the dummies of every term canonically
        (:func:`canonical_relabel`), so that equal terms merge."""
        e = expr.expand()
        for _ in range(max_passes):
            n = e.replace_multiple(self._rules, repeat=True).expand()
            if n == e and self._split:
                n = e.replace_multiple(self._split).expand()
            if n == e:
                return canonical_form(e)
            e = n
        raise RuntimeError("simplifier did not reach a fixpoint")


_ABSORB = [
    Replacement(g(_w["a"], _w["a"]), d),
    Replacement(g(_w["a"], _w["b"]) * _fn(_x, _w["b"], _y), _fn(_x, _w["a"], _y)),
    Replacement(delta(_w["a"], _w["b"]) * _fn(_x, _w["b"], _y), _fn(_x, _w["a"], _y)),
]


def absorb(expr: Expression) -> Expression:
    """Absorb full metrics ``g`` and Kronecker deltas into the heads (and field slots) they contract with; the
    only rewriting that never changes the form of anything else."""
    return expr.replace_multiple(_ABSORB, repeat=True).expand()


# ---------------------------------------------------------------------------------------------- canonical labels
def canonical_relabel(term: Expression, keep: Iterable[Expression] = ()) -> Expression:
    """Canonical relabelling of the dummy labels of a single product term.

    Dummy labels (those that occur twice) are renamed to the canonical pool ``lor(0), lor(1), ...`` (per
    space, skipping numbers used by ``keep`` or by free labels) and then canonicalised with Symbolica's
    ``canonize_tensors``, which takes the symmetry of ``g``, ``delta``, ``eps``, ``f``, ... into account. Free
    labels and the labels in ``keep`` are never renamed.
    """
    keep_s = {k.to_canonical_string() for k in keep}
    counts: dict[str, list] = {}
    for lab in labels_in(term):
        k = lab.to_canonical_string()
        counts.setdefault(k, [lab, 0])[1] += 1
    dummies = [v[0] for k, v in counts.items() if v[1] >= 2 and k not in keep_s]
    if not dummies:
        return term
    fixed = [v[0] for k, v in counts.items() if v[1] < 2 or k in keep_s]
    used: dict[str, set[int]] = {}
    for lab in fixed:
        used.setdefault(label_space(lab), set()).add(int(str(list(lab)[-1])))
    by_space: dict[str, list[Expression]] = {}
    for lab in dummies:
        by_space.setdefault(label_space(lab), []).append(lab)
    mapping = []
    pool: list[tuple[Expression, Expression]] = []
    for skey, labs in by_space.items():
        taken = used.get(skey, set())
        numbers = (n for n in itertools.count() if n not in taken)
        proto = _prototype(labs[0], skey)
        for lab in labs:
            new = proto(next(numbers))
            mapping.append((lab, new))
            pool.append((new, Expression.num(_space_marker(skey))))
    renamed = rename(term, mapping)
    canon, _, _ = renamed.canonize_tensors(pool)
    return canon


def canonical_form(expr: Expression, keep: Iterable[Expression] = ()) -> Expression:
    """Sum of :func:`canonical_relabel` over the terms of an expanded expression."""
    terms = list(expr) if expr.get_type() == AtomType.Add else [expr]
    out = Expression.num(0)
    for t in terms:
        out = out + canonical_relabel(t, keep)
    return out


_MARKERS: dict[str, int] = {}


def _space_marker(key: str) -> int:
    if key not in _MARKERS:
        _MARKERS[key] = len(_MARKERS) + 1
    return _MARKERS[key]


def _prototype(lab: Expression, skey: str):
    if skey == "lorentz":
        return lambda n: lor(n)
    if skey == "spinor":
        return lambda n: spin(n)
    if skey == "flavour":
        return lambda n: flav(n)
    fac = next(iter(lab))
    head = adj if lab.get_name() == adj.get_name() else fund
    return lambda n: head(fac, n)
