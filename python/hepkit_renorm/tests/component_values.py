"""Test support: exact values of coefficients at explicit components, and the translation of the hand-written
HepKIT rules (golden data) into the package's vocabulary.

Values are computed over the package's default prime field, in ``D = 6`` Lorentz dimensions with metric
``diag(+,-,-,-,-,-)`` (components 0..3 four-dimensional, 4..5 evanescent) and 8-dimensional Dirac matrices that
realise the BMHV algebra exactly::

    Gamma^mu = gamma^mu x 1 (mu < 4),  Gamma^4 = gamma5 x i sigma1,  Gamma^5 = gamma5 x i sigma2,  Gamma5 = gamma5 x 1

(Dirac representation of ``gamma^mu``, ``gamma5 = i gamma^0 gamma^1 gamma^2 gamma^3``), so that ``{Gamma5,
gammabar} = 0`` and ``[Gamma5, gammahat] = 0``. ``eps^{0123} = +1``. A Lorentz label that occurs twice is contracted
with the metric; free labels are contracted with polarisation vectors (upper components, also with the metric).
Group heads use the explicit generators of :mod:`hepkit_renorm.groups`. A rule is fingerprinted by its value for
random polarisations of all free labels, random momenta with ``sum p = 0`` and random parameter values (fixed
seeds): two rules with equal fingerprints at several points are equal with overwhelming probability.

This module only evaluates; it never decides a sign convention.
"""

from __future__ import annotations

import itertools
import random
import re
from collections.abc import Mapping, Sequence

import numpy as np
from symbolica import AtomType, Expression, S

from hepkit_renorm import tensors as tn
from hepkit_renorm.exact import PrimeField
from hepkit_renorm.groups import complex_rational, su_d_tensor, su_generators, su_structure_constants
from hepkit_renorm.spec import Model

DIM = 6
SPINOR_DIM = 8


class ComponentEvaluator:
    """Exact evaluation of coefficients (sums of products of heads) at explicit components over ``field``."""

    def __init__(self, model: Model, field: PrimeField) -> None:
        self.model, self.F = model, field
        F = field
        i = F.i
        self.eta = F.array(np.diag([1] + [-1] * (DIM - 1)).tolist())
        bar = np.diag([1, 1, 1, 1, 0, 0])
        self.etabar = F.reduce(self.eta * bar)
        self.etahat = F.reduce(self.eta * (1 - bar))
        s1 = np.array([[0, 1], [1, 0]], dtype=object)
        s2 = np.array([[0, -i], [i, 0]], dtype=object)
        s3 = np.array([[1, 0], [0, -1]], dtype=object)
        one2, zero2 = np.eye(2, dtype=int).astype(object), np.zeros((2, 2), dtype=object)
        g0 = np.block([[one2, zero2], [zero2, -one2]])
        gs = [np.block([[zero2, s], [-s, zero2]]) for s in (s1, s2, s3)]
        g4d = [F.reduce(g0)] + [F.reduce(g) for g in gs]
        g5 = F.reduce(i * F.matmul(F.matmul(g4d[0], g4d[1]), F.matmul(g4d[2], g4d[3])))
        one = np.eye(2, dtype=int).astype(object)
        gam = [F.reduce(np.kron(g, one)) for g in g4d]
        gam.append(F.reduce(np.kron(g5, F.reduce(i * s1))))
        gam.append(F.reduce(np.kron(g5, F.reduce(i * s2))))
        self.gamma = np.array(gam, dtype=object)
        self.gamma5 = F.reduce(np.kron(g5, one))
        unit = F.array(np.eye(SPINOR_DIM, dtype=int).tolist())
        self.PL = F.reduce((unit - self.gamma5) * F.inv(2))
        self.PR = F.reduce((unit + self.gamma5) * F.inv(2))
        mask = np.array([1, 1, 1, 1, 0, 0], dtype=object).reshape(DIM, 1, 1)
        self.gammabar = F.reduce(self.gamma * mask)
        self.gammahat = F.reduce(self.gamma * (1 - mask))
        eps = np.zeros((DIM,) * 4, dtype=object)
        for perm in itertools.permutations(range(4)):
            inv = sum(1 for x, y in itertools.combinations(perm, 2) if x > y)
            eps[perm] = 1 if inv % 2 == 0 else F.p - 1
        self.eps = eps

    # ------------------------------------------------------------------ label spaces
    def dim(self, label: Expression) -> int:
        key = tn.label_space(label)
        if key == "lorentz":
            return DIM
        if key == "spinor":
            return SPINOR_DIM
        if key == "flavour":
            return self.model.generations
        kind, fac = key.split(":")
        fa = self.model.factors[fac]
        return fa.dim_adjoint if kind == "adjoint" else fa.n

    def number(self, num: Expression) -> int:
        re, im = complex_rational(num)
        return (self.F(re) + self.F(im) * self.F.i) % self.F.p

    # ------------------------------------------------------------------ heads
    def _head(self, fn: Expression, momenta: Mapping[int, np.ndarray]) -> tuple[np.ndarray | int, list[Expression]]:
        F = self.F
        name = fn.get_name()
        args = list(fn)
        H = tn.HEADS
        if name == H["p"].get_name():
            return momenta[int(args[0].to_canonical_string())], [args[1]]
        for sname, metric in (("sp", self.eta), ("spbar", self.etabar), ("sphat", self.etahat)):
            if name == H[sname].get_name():
                a, b = (momenta[int(x.to_canonical_string())] for x in args)
                return int(F.reduce(np.dot(a, np.dot(metric, b)))), []
        table = {"g": self.eta, "gbar": self.etabar, "ghat": self.etahat, "eps": self.eps, "gamma": self.gamma,
                 "gammabar": self.gammabar, "gammahat": self.gammahat, "gamma5": self.gamma5, "PL": self.PL,
                 "PR": self.PR}
        for hname, arr in table.items():
            if name == H[hname].get_name():
                return arr, args
        if name == H["delta"].get_name():
            return F.array(np.eye(self.dim(args[0]), dtype=int).tolist()), args
        fac = self.model.factors[tn.label_space(args[0]).split(":")[1]]
        if name == H["T"].get_name():
            return su_generators(fac.n, F), args
        if name == H["f"].get_name():
            return su_structure_constants(fac.n, F), args
        if name == H["dsym"].get_name():
            return su_d_tensor(fac.n, F), args
        if name == H["epsf"].get_name():
            n = fac.n
            e = np.zeros((n,) * n, dtype=object)
            for perm in itertools.permutations(range(n)):
                inv = sum(1 for x, y in itertools.combinations(perm, 2) if x > y)
                e[perm] = 1 if inv % 2 == 0 else F.p - 1
            return e, args
        raise ValueError(f"cannot evaluate head {fn}")

    def value(self, expr: Expression, polarisations: Mapping[str, np.ndarray], momenta: Mapping[int, np.ndarray],
              params: Mapping[str, int]) -> int:
        """Value of ``expr`` with every free label contracted with ``polarisations[label canonical string]``."""
        F = self.F
        total = 0
        expr = expr.expand()
        terms = list(expr) if expr.get_type() == AtomType.Add else [expr]
        for term in terms:
            if term.get_type() == AtomType.Num and term.to_canonical_string() == "0":
                continue
            total = (total + self._term(term, polarisations, momenta, params)) % F.p
        return total

    def _term(self, term, polarisations, momenta, params) -> int:
        F = self.F
        scalar = 1
        arrays: list[np.ndarray] = []
        labs: list[list[Expression]] = []
        for fa in (list(term) if term.get_type() == AtomType.Mul else [term]):
            t = fa.get_type()
            reps = 1
            if t == AtomType.Pow:
                base, ex = list(fa)
                reps = int(str(ex))
                fa, t = base, base.get_type()
            if t == AtomType.Num:
                scalar = scalar * pow(self.number(fa), reps, F.p) % F.p
            elif t == AtomType.Var:
                v = params[fa.get_name()]
                scalar = scalar * (pow(v, reps, F.p) if reps >= 0 else pow(F.inv(v), -reps, F.p)) % F.p
            elif t == AtomType.Fn:
                arr, args = self._head(fa, momenta)
                if isinstance(arr, int):
                    scalar = scalar * (pow(arr, reps, F.p) if reps >= 0 else pow(F.inv(arr), -reps, F.p)) % F.p
                    continue
                if reps < 0:
                    raise ValueError(f"negative power of a tensor: {fa}")
                for _ in range(reps):
                    arrays.append(arr)
                    labs.append(args)
            else:
                raise ValueError(f"cannot evaluate factor {fa}")
        counts: dict[str, int] = {}
        for ls in labs:
            for lab in ls:
                counts[lab.to_canonical_string()] = counts.get(lab.to_canonical_string(), 0) + 1
        for key, n in counts.items():
            if n == 1:
                if key not in polarisations:
                    raise KeyError(f"no polarisation for the free label {key}")
                arrays.append(polarisations[key])
                labs.append([Expression.parse(key)])
            elif n != 2:
                raise ValueError(f"label {key} occurs {n} times")
        # letters: one per occurrence of a Lorentz label (paired with the metric), one per other label
        letters = iter("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ")
        spec_in: list[str] = []
        seen_lorentz: dict[str, str] = {}
        other: dict[str, str] = {}
        extra: list[tuple[str, str]] = []
        for ls in labs:
            s = ""
            for lab in ls:
                key = lab.to_canonical_string()
                if tn.label_space(lab) == "lorentz":
                    c = next(letters)
                    if key in seen_lorentz:
                        extra.append((seen_lorentz[key], c))
                    else:
                        seen_lorentz[key] = c
                    s += c
                else:
                    if key not in other:
                        other[key] = next(letters)
                    s += other[key]
            spec_in.append(s)
        ops = list(arrays) + [self.eta] * len(extra)
        specs = spec_in + [a + b for a, b in extra]
        if not ops:
            return scalar
        return _contract(F, ops, specs) * scalar % F.p


def _contract(F: PrimeField, ops: list[np.ndarray], specs: list[str]) -> int:
    """Full contraction of a tensor network (every letter occurs exactly twice), pairwise with ``tensordot``."""
    ops = [np.asarray(o, dtype=object) for o in ops]
    specs = list(specs)
    while len(ops) > 1:
        best = None
        for i, j in itertools.combinations(range(len(ops)), 2):
            shared = set(specs[i]) & set(specs[j])
            if shared:
                size = ops[i].size * ops[j].size
                for c in shared:
                    size //= ops[i].shape[specs[i].index(c)]
                if best is None or size < best[0]:
                    best = (size, i, j)
        if best is None:  # disconnected pieces: outer product of the first two
            i, j = 0, 1
        else:
            i, j = best[1], best[2]
        shared = [c for c in specs[i] if c in specs[j]]
        res = np.tensordot(ops[i], ops[j], axes=([specs[i].index(c) for c in shared],
                                                 [specs[j].index(c) for c in shared]))
        spec = "".join(c for c in specs[i] if c not in shared) + "".join(c for c in specs[j] if c not in shared)
        res = F.reduce(np.asarray(res, dtype=object))
        ops = [o for k, o in enumerate(ops) if k not in (i, j)] + [res]
        specs = [s for k, s in enumerate(specs) if k not in (i, j)] + [spec]
    (op,), (spec,) = ops, specs
    for c in sorted(set(spec)):  # a trace within one operand
        a, b = [k for k, x in enumerate(spec) if x == c]
        op = np.trace(op, axis1=a, axis2=b)
        spec = spec.replace(c, "")
    return int(np.asarray(op, dtype=object).reshape(()).item()) % F.p


# ---------------------------------------------------------------------------------------------- fingerprints
class Kinematics:
    """Deterministic random momenta (``sum p = 0``), parameters and polarisations over a prime field."""

    def __init__(self, field: PrimeField, n_legs: int, params: Sequence[str], seed: int) -> None:
        rng = random.Random(seed)
        self.rng, self.F = rng, field
        moms = [field.array([rng.randrange(-50, 50) for _ in range(DIM)]) for _ in range(n_legs - 1)]
        last = field.reduce(-sum(moms, field.zeros(DIM))) if moms else field.zeros(DIM)
        self.momenta = {k + 1: m for k, m in enumerate(moms + [last])}
        self.momenta[0] = field.array([rng.randrange(-50, 50) for _ in range(DIM)])
        self.params = {name: field(rng.randrange(2, 1000)) for name in params}
        self.params[tn.d.get_name()] = field(DIM)

    def vector(self, n: int) -> np.ndarray:
        return self.F.array([self.rng.randrange(-1000, 1000) for _ in range(n)])


def unit_vector(field: PrimeField, n: int, k: int) -> np.ndarray:
    return field.array([int(i == k) for i in range(n)])


def model_parameters(model: Model) -> list[str]:
    names = [model.param(n).get_name() for n in model.parameters]
    names += [model.coupling(f).get_name() for f in model.factors]
    names.append(model.gauge_fixing.xi.get_name())
    return names


# ---------------------------------------------------------------------------------------------- hand-written rules
_W = {n: S(n) for n in ("x_", "s_", "t_")}
_BAR, _HAT = S("hk_bar"), S("hk_hat")
_GAMMA, _GAMMABAR, _GAMMAHAT = S("hk_Gamma"), S("hk_GammaBar"), S("hk_GammaHat")
_IX, _DM = S("hk_IX"), S("hk_DM")
_UFO_LORENTZ = {"Gamma": ("lorentz", "spinor", "spinor"), "ProjM": ("spinor", "spinor"),
                "ProjP": ("spinor", "spinor"), "Metric": ("lorentz", "lorentz"), "P": ("lorentz", None),
                "PSlash": ("spinor", "spinor")}


def parse_ufo(text: str) -> Expression:
    """Parse a UFO / HepKIT string of the hand-written model: namespaces dropped, ``bmhv::bar``/``bmhv::hat``
    markers absorbed into ``Gamma`` (``hk_GammaBar``, ``hk_GammaHat``), indices as ``hk_IX(k)`` / ``hk_DM(n)``."""
    s = text.replace("UFO::{}::", "").replace("bmhv::bar(", "hk_bar(").replace("bmhv::hat(", "hk_hat(")
    s = s.replace("1\U0001d456", "hk_I").replace("idx(1,", "hk_IX(").replace("dummy(", "hk_DM(")
    s = re.sub(r"(?<![\w:])(Gamma|ProjM|ProjP|Metric|PSlash|Identity|P|T|f)\(", r"hk_\1(", s)
    e = Expression.parse(s)
    x, s_, t_ = _W["x_"], _W["s_"], _W["t_"]
    e = e.replace(_BAR(x) * _GAMMA(x, s_, t_), _GAMMABAR(x, s_, t_), repeat=True)
    e = e.replace(_HAT(x) * _GAMMA(x, s_, t_), _GAMMAHAT(x, s_, t_), repeat=True)
    return e.replace(S("hk_I"), tn.I)


def convert_ufo(expr: Expression, leg_label, dummy_label, params: Mapping[str, Expression]) -> Expression:
    """Map a parsed UFO expression to the package vocabulary.

    ``leg_label(k, space)`` returns the label of leg ``k`` (1-based) in the space ``"lorentz"``, ``"spinor"``,
    ``"colour_fund"``, ``"colour_adj"`` or ``"colour"`` (the colour slot of the leg, of either kind); ``dummy_label(key, space)`` a dummy label; ``params`` maps UFO parameter
    names to model symbols. Colour: ``T(a, i, jb)`` (``i`` the quark leg, ``jb`` the antiquark leg) is the matrix
    element ``(T^a)_{jb i}`` in the row convention of the package (row = barred field).
    """
    def lab(arg: Expression, space: str) -> Expression:
        if arg.get_type() == AtomType.Fn and arg.get_name() == _IX.get_name():
            return leg_label(int(next(iter(arg)).to_canonical_string()), space)
        if arg.get_type() == AtomType.Fn and arg.get_name() == _DM.get_name():
            return dummy_label(("dm", next(iter(arg)).to_canonical_string(), space), space)
        n = int(arg.to_canonical_string())
        return leg_label(n, space) if n > 0 else dummy_label(("col", n, space), space)

    def conv(e: Expression) -> Expression:
        t = e.get_type()
        if t in (AtomType.Add, AtomType.Mul):
            out = Expression.num(0 if t == AtomType.Add else 1)
            for a in e:
                out = out + conv(a) if t == AtomType.Add else out * conv(a)
            return out
        if t == AtomType.Pow:
            b, ex = list(e)
            return conv(b) ** ex
        if t == AtomType.Var:
            name = e.get_name().split("::")[-1]
            return params.get(name, e)
        if t == AtomType.Num:
            return e
        name = e.get_name().split("::")[-1]
        args = list(e)
        if name == "hk_Gamma":
            return tn.gamma(lab(args[0], "lorentz"), lab(args[1], "spinor"), lab(args[2], "spinor"))
        if name == "hk_GammaBar":
            return tn.gammabar(lab(args[0], "lorentz"), lab(args[1], "spinor"), lab(args[2], "spinor"))
        if name == "hk_GammaHat":
            return tn.gammahat(lab(args[0], "lorentz"), lab(args[1], "spinor"), lab(args[2], "spinor"))
        if name == "hk_ProjM":
            return tn.PL(lab(args[0], "spinor"), lab(args[1], "spinor"))
        if name == "hk_ProjP":
            return tn.PR(lab(args[0], "spinor"), lab(args[1], "spinor"))
        if name == "hk_Metric":
            return tn.g(lab(args[0], "lorentz"), lab(args[1], "lorentz"))
        if name == "hk_P":
            return tn.p(int(args[1].to_canonical_string()), lab(args[0], "lorentz"))
        if name == "hk_PSlash":
            mu = dummy_label(("pslash", args[0].to_canonical_string()), "lorentz")
            return tn.gamma(mu, lab(args[0], "spinor"), lab(args[1], "spinor")) * tn.p(0, mu)
        if name == "hk_T":
            return tn.T(lab(args[0], "colour_adj"), lab(args[2], "colour_fund"), lab(args[1], "colour_fund"))
        if name == "hk_f":
            return tn.f(*(lab(a, "colour_adj") for a in args))
        if name == "hk_Identity":
            return tn.delta(lab(args[0], "colour"), lab(args[1], "colour"))
        raise ValueError(f"unknown UFO head {name}")

    return conv(expr)


def hand_rule_expression(rule: Mapping, leg_label, dummy_label, params: Mapping[str, Expression]) -> Expression:
    """Sum over (colour i, Lorentz j) of colour_i * lorentz_j * coupling_ij of a digest rule (= Feynman rule)."""
    total = Expression.num(0)
    colours = [parse_ufo(c) for c in rule["color_structures"]]
    lorentz = [parse_ufo(lz) for lz in rule["lorentz_structures"]]
    for i, row in enumerate(rule["couplings"]):
        for j, c in enumerate(row):
            if c is None:
                continue
            term = colours[i] * lorentz[j] * parse_ufo(c["expression"])
            total = total + convert_ufo(term, leg_label, dummy_label, params)
    return total

