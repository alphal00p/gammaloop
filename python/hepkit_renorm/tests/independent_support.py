"""Independent exact evaluator used by the ``test_independent_*`` review tests.

Written for the independent review of milestones MS0/MS1 without the package's evaluators
(:mod:`hepkit_renorm.fingerprint`, ``tests/component_values.py``). It evaluates

* expressions in the package's coefficient vocabulary (vertex functions, report formulas typed in that
  vocabulary), and
* the UFO strings of the hand-written HepKIT model ``tests/data/sm_unbroken_bmhv.json``,

contracted with explicit polarisations, exactly, in a prime field ``GF(P)`` with ``P`` a 62-bit prime in which
``-1`` and ``3`` are squares (``i`` and the ``sqrt(3)`` of the Gell-Mann matrix ``lambda^8``).

Conventions (report 03_setup sec:bmhv): metric ``diag(+,-,-,...)``; components ``0..3`` are four-dimensional
(barred), ``4..D-1`` evanescent (hatted); every tensor is stored with upper Lorentz indices and a contracted
Lorentz pair is contracted with the metric ``g_{mu nu}``; a free label is contracted directly with the
polarisation supplied for it (so a polarisation is a lower-index object). The Dirac matrices are the *Weyl*
representation in four dimensions (deliberately not the Dirac representation used by the package's evaluators),
extended to ``D = 6`` by ``Gamma^4 = gamma5 x i sigma_1``, ``Gamma^5 = gamma5 x i sigma_2``; ``gamma5`` is
computed as ``i gamma^0 gamma^1 gamma^2 gamma^3`` (not typed). A closed spinor trace is normalised to
``tr 1 = 4`` as in the report. ``eps^{0123} = +1`` (upper indices). SU(N) generators are the generalised
Gell-Mann matrices over two, ``f^{abc} = -2 i Tr([T^a, T^b] T^c)``, ``d^{abc} = 2 Tr(T^a {T^b, T^c})``,
``epsf(1, 2, ..) = +1``.
"""

from __future__ import annotations

import itertools
import re
import string
from collections.abc import Mapping, Sequence
from fractions import Fraction

import numpy as np
from symbolica import AtomType, Expression

# ---------------------------------------------------------------------------------------------- the prime field


def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d, s = d // 2, s + 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _find_prime() -> int:
    p = 2**62 + 1
    while not (p % 12 == 1 and _is_prime(p)):
        p += 2
    return p


P = _find_prime()


def _sqrt(a: int) -> int:
    """Tonelli-Shanks square root modulo ``P`` (the smaller root)."""
    a %= P
    assert pow(a, (P - 1) // 2, P) == 1, "not a square"
    q, s = P - 1, 0
    while q % 2 == 0:
        q, s = q // 2, s + 1
    z = 2
    while pow(z, (P - 1) // 2, P) != P - 1:
        z += 1
    m, c, t, r = s, pow(z, q, P), pow(a, q, P), pow(a, (q + 1) // 2, P)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2, i = t2 * t2 % P, i + 1
        b = pow(c, 1 << (m - i - 1), P)
        m, c, t, r = i, b * b % P, t * b * b % P, r * b % P
    return min(r, P - r)


IU = _sqrt(-1)
SQRT3 = _sqrt(3)


def fe(x: int | Fraction) -> int:
    """Image of a rational number."""
    x = Fraction(x)
    return x.numerator % P * pow(x.denominator, -1, P) % P


def arr(values) -> np.ndarray:
    return np.array(values, dtype=object) % P


# ---------------------------------------------------------------------------------------------- explicit algebra
def _kron(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.kron(a, b) % P


def dirac_matrices(D: int) -> tuple[list[np.ndarray], np.ndarray]:
    """``(gammas, gamma5)``: Weyl representation for ``D = 4``, extended to ``D = 6``."""
    i = IU
    one2 = arr([[1, 0], [0, 1]])
    zero2 = arr([[0, 0], [0, 0]])
    sig = [arr([[0, 1], [1, 0]]), arr([[0, -i], [i, 0]]), arr([[1, 0], [0, -1]])]
    g0 = np.block([[zero2, one2], [one2, zero2]]) % P
    gs = [g0] + [np.block([[zero2, s], [(-s) % P, zero2]]) % P for s in sig]
    g5 = (i * (gs[0].dot(gs[1]).dot(gs[2]).dot(gs[3]) % P)) % P
    if D == 4:
        return gs, g5
    if D == 6:
        out = [_kron(g, one2) for g in gs]
        out += [_kron(g5, (i * sig[0]) % P), _kron(g5, (i * sig[1]) % P)]
        return out, _kron(g5, one2)
    raise ValueError("D must be 4 or 6")


def su_generators(n: int) -> np.ndarray:
    """Generalised Gell-Mann matrices over two (``n = 2``: Pauli / 2; ``n = 3``: Gell-Mann / 2)."""
    i = IU
    half = fe(Fraction(1, 2))
    gens = []
    for k in range(1, n):
        for j in range(k):
            s = [[0] * n for _ in range(n)]
            s[j][k] = s[k][j] = half
            a = [[0] * n for _ in range(n)]
            a[j][k], a[k][j] = (-i * half) % P, (i * half) % P
            gens += [s, a]
        h = [[0] * n for _ in range(n)]
        norm = pow({1: 2, 2: 2 * SQRT3}[k], -1, P)  # 1 / sqrt(2 k (k+1)), only n <= 3 is needed
        for r in range(k):
            h[r][r] = norm
        h[k][k] = (-k * norm) % P
        gens.append(h)
    # order (S01, A01, H1, S02, A02, S12, A12, H2) = (lambda_1, ..., lambda_8) / 2; for n = 2 the Pauli matrices
    return arr(gens)


def structure_constants(t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    prod = np.einsum("aij,bjk->abik", t, t) % P
    comm = (prod - np.transpose(prod, (1, 0, 2, 3))) % P
    anti = (prod + np.transpose(prod, (1, 0, 2, 3))) % P
    f = np.einsum("abik,cki->abc", comm, t) * ((-2 * IU) % P) % P
    dd = np.einsum("abik,cki->abc", anti, t) * 2 % P
    return f, dd


def levi_civita(n: int) -> np.ndarray:
    e = np.zeros((n,) * n, dtype=object)
    for perm in itertools.permutations(range(n)):
        inv = sum(1 for a, b in itertools.combinations(range(n), 2) if perm[a] > perm[b])
        e[perm] = 1 if inv % 2 == 0 else P - 1
    return e


# ---------------------------------------------------------------------------------------------- numbers
_IM = re.compile(r"^(?P<num>[+-]?\d*)𝑖(?:/(?P<den>\d+))?$")


def number(text: str) -> int:
    """Image of a Symbolica (complex) rational such as ``-7/3-2𝑖/9``, ``2𝑖/3`` or ``-𝑖``."""
    if "𝑖" not in text:
        return fe(Fraction(text))
    cut = max(text.rfind("+", 1, text.find("𝑖")), text.rfind("-", 1, text.find("𝑖")))
    re_part, im_part = (text[:cut], text[cut:]) if cut > 0 else ("", text)
    m = _IM.match(im_part)
    if not m:
        raise ValueError(f"cannot read the number {text!r}")
    c = {"": 1, "+": 1, "-": -1}.get(m.group("num"))
    c = int(m.group("num")) if c is None else c
    val = fe(Fraction(re_part)) if re_part else 0
    return (val + fe(Fraction(c, int(m.group("den") or 1))) * IU) % P


# ---------------------------------------------------------------------------------------------- the evaluator
class Evaluator:
    """Exact evaluation in ``D`` dimensions (4 or 6) of expressions in the package vocabulary.

    ``params`` maps short parameter names (``g_s``, ``g``, ``gp``, ``y_u``, ``lambda``, ``m2``, ``xi``, ...) to
    field elements; ``momenta`` maps the momentum number of a leg to its (upper-index) components; ``pol`` maps
    the canonical string of a free label to its polarisation vector.
    """

    def __init__(self, D: int, params: Mapping[str, int], factor_dims: Mapping[str, int]) -> None:
        self.D = D
        self.params = dict(params)
        self.gammas, self.g5 = dirac_matrices(D)
        self.S = self.gammas[0].shape[0]
        self.metric = arr([1] + [-1] * (D - 1))
        self.bar_mask = arr([1] * 4 + [0] * (D - 4))
        self.hat_mask = arr([0] * 4 + [1] * (D - 4))
        self.gamma = arr(np.stack(self.gammas))
        self.gens = {F: su_generators(n) for F, n in factor_dims.items()}
        self.fd = {F: structure_constants(t) for F, t in self.gens.items()}
        eps = np.zeros((D,) * 4, dtype=object)
        eps[:4, :4, :4, :4] = levi_civita(4)
        self.eps = eps
        self.factor_dims = dict(factor_dims)

    # -------------------------------------------------------------- label spaces
    @staticmethod
    def label_key(label: Expression) -> str:
        return label.to_canonical_string()

    def space(self, label: Expression) -> tuple[str, str | None]:
        name = label.get_name().split("::")[-1]
        if name == "lor":
            return ("lor", None)
        if name == "spin":
            return ("spin", None)
        if name == "flav":
            return ("flav", None)
        fac = next(iter(label)).get_name().split("::")[-1]
        return (name, fac)

    def dim(self, sp: tuple[str, str | None]) -> int:
        kind, fac = sp
        if kind == "lor":
            return self.D
        if kind == "spin":
            return self.S
        if kind == "flav":
            return 1
        n = self.factor_dims[fac]
        return n if kind == "fund" else n * n - 1

    # -------------------------------------------------------------- public
    def value(self, expr: Expression, pol: Mapping[str, np.ndarray], momenta: Mapping[int, np.ndarray]) -> int:
        e = expr.expand()
        terms = list(e) if e.get_type() == AtomType.Add else [e]
        total = 0
        for t in terms:
            total = (total + self._term(t, pol, momenta)) % P
        return total

    # -------------------------------------------------------------- internals
    def _scalar(self, e: Expression, momenta) -> int:
        """Value of a sub-expression without labels (numbers, parameters, scalar products, sums thereof)."""
        t = e.get_type()
        if t == AtomType.Num:
            return number(e.to_canonical_string())
        if t == AtomType.Var:
            name = e.get_name().split("::")[-1]
            if name == "d":
                return self.D % P
            if name not in self.params:
                raise KeyError(f"no value for the parameter {name!r}")
            return self.params[name]
        if t == AtomType.Add:
            return sum(self._scalar(a, momenta) for a in e) % P
        if t == AtomType.Mul:
            out = 1
            for a in e:
                out = out * self._scalar(a, momenta) % P
            return out
        if t == AtomType.Pow:
            base, ex = list(e)
            n = int(ex.to_canonical_string())
            b = self._scalar(base, momenta)
            return pow(b, n, P) if n >= 0 else pow(pow(b, -1, P), -n, P)
        if t == AtomType.Fn:
            name = e.get_name().split("::")[-1]
            if name in ("sp", "spbar", "sphat"):
                k, l = (int(a.to_canonical_string()) for a in e)
                w = self.metric * {"sp": 1, "spbar": self.bar_mask, "sphat": self.hat_mask}[name] % P
                return int(np.sum(momenta[k] * w * momenta[l]) % P)
        raise ValueError(f"not a scalar: {e}")

    def _tensor(self, e: Expression, momenta) -> tuple[np.ndarray, list[Expression]] | None:
        name = e.get_name().split("::")[-1]
        args = list(e)
        if name in ("sp", "spbar", "sphat"):
            return None
        if name == "p":
            k = int(args[0].to_canonical_string())
            return momenta[k], [args[1]]
        if name in ("g", "gbar", "ghat"):
            mask = {"g": 1, "gbar": self.bar_mask, "ghat": self.hat_mask}[name]
            return np.diag(self.metric * mask % P), args
        if name == "eps":
            return self.eps, args
        if name in ("gamma", "gammabar", "gammahat"):
            if name == "gamma":
                return self.gamma, args
            mask = self.bar_mask if name == "gammabar" else self.hat_mask
            return self.gamma * mask.reshape(-1, 1, 1) % P, args
        if name == "gamma5":
            return self.g5, args
        if name in ("PL", "PR"):
            one = np.eye(self.S, dtype=object)
            sgn = -1 if name == "PL" else 1
            return (one + sgn * self.g5) * fe(Fraction(1, 2)) % P, args
        if name == "delta":
            n = self.dim(self.space(args[0]))
            assert self.space(args[0]) == self.space(args[1]), "delta between different spaces"
            return np.eye(n, dtype=object), args
        if name == "T":
            fac = self.space(args[0])[1]
            return self.gens[fac], args
        if name in ("f", "dsym"):
            fac = self.space(args[0])[1]
            return self.fd[fac][0 if name == "f" else 1], args
        if name == "epsf":
            return levi_civita(len(args)), args
        raise ValueError(f"unknown head {name}")

    def _term(self, term: Expression, pol, momenta) -> int:
        factors = list(term) if term.get_type() == AtomType.Mul else [term]
        scalar = 1
        tensors: list[tuple[np.ndarray, list[Expression]]] = []
        for fa in factors:
            t = fa.get_type()
            if t == AtomType.Fn and self._tensor(fa, momenta) is not None:
                tensors.append(self._tensor(fa, momenta))
            elif t == AtomType.Pow and next(iter(fa)).get_type() == AtomType.Fn and \
                    self._tensor(next(iter(fa)), momenta) is not None:
                base, ex = list(fa)
                n = int(ex.to_canonical_string())
                assert n > 0
                tensors += [self._tensor(base, momenta)] * n
            else:
                scalar = scalar * self._scalar(fa, momenta) % P
        return scalar * contract(tensors, pol, self) % P


def contract(tensors: Sequence[tuple[np.ndarray, Sequence]], pol: Mapping[str, np.ndarray], ev) -> int:
    """Contract tensors whose axes carry labels (Symbolica labels or ``(key, space)`` pairs): Lorentz pairs with the
    metric, other pairs by summation, free labels with ``pol[key]``; closed spinor traces normalised to 4."""
    keyed = []
    spaces: dict[str, tuple] = {}
    for a, labs in tensors:
        ks = []
        for lab in labs:
            if isinstance(lab, tuple):
                k, sp = lab
            else:
                k, sp = ev.label_key(lab), ev.space(lab)
            spaces[k] = sp
            ks.append(k)
        keyed.append((a, ks))
    counts: dict[str, int] = {}
    for _, ks in keyed:
        for k in ks:
            counts[k] = counts.get(k, 0) + 1
    norm = 1
    # closed spinor traces: connected components (through spinor labels) without a free spinor label
    spin_keys = [k for k, sp in spaces.items() if sp[0] == "spin"]
    if spin_keys:
        parent = list(range(len(keyed)))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        owners: dict[str, list[int]] = {}
        for i, (_, ks) in enumerate(keyed):
            for k in ks:
                if spaces[k][0] == "spin":
                    owners.setdefault(k, []).append(i)
        for k, own in owners.items():
            for j in own[1:]:
                parent[find(j)] = find(own[0])
        comps: dict[int, list[str]] = {}
        for k, own in owners.items():
            comps.setdefault(find(own[0]), []).append(k)
        for ks in comps.values():
            if all(counts[k] == 2 for k in ks):
                norm = norm * fe(Fraction(4, ev.S)) % P
    ops = []
    letters: dict[str, str] = {}
    pool = iter(string.ascii_letters)

    def letter(k):
        if k not in letters:
            letters[k] = next(pool)
        return letters[k]
    seen: dict[str, int] = {}
    for a, ks in keyed:
        sub = ""
        for k in ks:
            seen[k] = seen.get(k, 0) + 1
            if counts[k] == 2 and seen[k] == 2 and spaces[k][0] == "lor":
                sub += letter(k + "#2")
            else:
                sub += letter(k)
        ops.append((a, sub))
    for k, c in counts.items():
        if c > 2:
            raise ValueError(f"label {k} occurs {c} times")
        if c == 2 and spaces[k][0] == "lor":
            ops.append((np.diag(ev.metric), letter(k) + letter(k + "#2")))
        if c == 1:
            if k not in pol:
                raise KeyError(f"no polarisation for the free label {k}")
            ops.append((pol[k], letter(k)))
    return _pairwise(ops) * norm % P


def _pairwise(ops: list[tuple[np.ndarray, str]]) -> int:
    """Full contraction of the operands (every letter occurs exactly twice overall), pairwise and greedily by the
    size of the intermediate result, reducing modulo ``P`` after each step (object arrays: no overflow)."""
    if not ops:
        return 1
    ops = [(np.asarray(a, dtype=object), s) for a, s in ops]
    while len(ops) > 1:
        best = None
        for i, j in itertools.combinations(range(len(ops)), 2):
            (a, sa), (b, sb) = ops[i], ops[j]
            rest = "".join(s for n, (_, s) in enumerate(ops) if n not in (i, j))
            out = "".join(dict.fromkeys(c for c in sa + sb if (sa + sb).count(c) == 1 or c in rest))
            out = "".join(c for c in out if not (c in sa and c in sb and c not in rest))
            size = 1
            dims = {c: n for s, x in ((sa, a), (sb, b)) for c, n in zip(s, x.shape)}
            for c in out:
                size *= dims[c]
            if best is None or size < best[0]:
                best = (size, i, j, out)
        _, i, j, out = best
        (a, sa), (b, sb) = ops[i], ops[j]
        res = np.asarray(np.asarray(np.einsum(f"{sa},{sb}->{out}", a, b), dtype=object) % P, dtype=object)
        ops = [op for n, op in enumerate(ops) if n not in (i, j)] + [(res, out)]
    a, s = ops[0]
    if s:  # letters repeated inside one operand (a trace)
        a = np.asarray(np.einsum(f"{s}->", a), dtype=object)
    return int(a) % P


# ---------------------------------------------------------------------------------------------- random data
class Sampler:
    """Deterministic random data in ``GF(P)``."""

    def __init__(self, seed: int) -> None:
        self.rng = np.random.default_rng(seed)

    def element(self) -> int:
        return int(self.rng.integers(1, 2**62)) % P

    def vector(self, n: int) -> np.ndarray:
        return np.array([self.element() for _ in range(n)], dtype=object)

    def momenta(self, n: int, D: int) -> dict[int, np.ndarray]:
        moms = {k: self.vector(D) for k in range(1, n)}
        moms[n] = (-sum(moms.values())) % P if n > 1 else np.zeros(D, dtype=object)
        return moms
