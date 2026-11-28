"""Notation and rendering of every mathematical expression of the gamma5 documentation site.

Every expression on the page is a Symbolica expression in the *site notation* defined here, printed to MathML by
Symbolica's own printer (``Expression.format(mode=PrintMode.MathML)``) with the custom print callbacks registered by
this module, and optionally explored with Symbolica's MathML pager (the viewer of the Marimo notebooks), fed from a
pre-computed payload. This docstring is the reference for the other build modules.

Contents: 1. public interface, 2. the site notation, 3. the converters from ``hepkit_renorm``, 4. the interactive
payload, 5. assets.

1. Public interface (DATA_FORMAT.md, sections 2, 8 and 9)
=========================================================

``M(obj, *, interactive=False, display=True, page_size=None, budget=VIEW_BUDGET) -> dict``
    The math object ``{"html", "text", "src", "view"}``. ``obj`` is a Symbolica ``Expression`` in the site notation,
    a source string in the site notation (parsed with :func:`parse`), or a ``hepkit_renorm`` object that
    :func:`convert` accepts (``Monomial``, ``Density``, ``LocalFunctional``, ``VertexRule``). ``html`` is
    ``<math>...</math>`` from Symbolica's MathML printer; ``display=True`` gives ``displaystyle="true"`` (the
    default; tables and cards), ``display=False`` the inline style (running text). ``text`` is a plain Unicode
    form (search, ``title``/``aria-label``), ``src`` the Symbolica source in the site notation (``parse(src)``
    restores the expression). ``view`` is ``None`` unless ``interactive=True`` (section 4).

``prose(text) -> str``
    Rich text (R of DATA_FORMAT.md): ``[[ source ]]`` spans become inline MathML; the allowed tags ``p em strong a ul
    ol li code sub sup span br`` (and their attributes ``class``, ``title``, ``data-report-page``, ``data-term``,
    ``data-class``, ``data-theory``) are kept, everything else is escaped. A span may contain brackets
    (``[[ comm(T(a), T(b)) ]]``); ``\\[[`` writes a literal ``[[``. A span that does not parse raises ``ValueError``
    naming the span.

``parse(src) -> Expression``
    Parses site-notation source (default namespace ``gdocs``; ``**`` is accepted for ``^``).

``convert(obj, **options) -> Expression``, ``operator``, ``vertex``, ``coefficient``, ``legs``, ``leg``,
``operator_parts``, ``field``: section 3.

``write_assets(site_dir)``: section 5. ``VIEW_BUDGET``: default size budget of a ``view`` payload (section 4).

2. The site notation
====================

Symbols live in the Symbolica namespace ``gdocs``; :func:`parse` uses it as the default namespace, so sources are
written without prefix. Products are printed in the order in which the symbols were registered (numbers first,
then couplings and parameters, loop factors, counterterm names, tensors, momenta, and the ordered field product
``ord(...)`` last), so a product never needs to be written in a particular order. Where order matters
(Grassmann-odd fields, spinor chains, matrices in group space) the expression uses the ordered product
``ord(x1, x2, ...)``, which prints its arguments in the given order.

Index placement: an index argument is printed at the default position of its slot (listed below); ``lo(x)`` forces
a subscript, ``up(x)`` a superscript. ``prime(x)`` prints ``x′``. Index symbols: ``mu nu rho sigma alpha beta
kappa tau`` print as Greek letters; every other symbol prints as its name (``A``, ``C``, ``I``, ``a``, ``i``, ...).
The bare letters of fields (``B``, ``D``, ``E``, ``G``, ``L``, ...) double as index letters: as atoms they print as
the letter, as functions as the field.

Couplings and parameters (symbols)::

    g_s  g  gp (g′)  y_u  y_d  y_e  lambda (λ)  m2 (m²)  xi (ξ)        Y(f) hypercharge Y_f, Nc(f) N^c_f
    v  M_W  M_Z  M_h  m_u  m_d  m_e  m_f  c_W  s_W  e                     broken phase
    ep (ϵ, the regulator, lunate glyph)  eta (η, evanescent tag)  d (dimension)  pi (π, Symbolica's constant)
    lf = 1/(16π²)   lfe = (1/(16π²)) (1/ϵ)                                loop factors as in the report

Counterterm names::

    dZ(X) -> δZ_X;  dZ(X, r) -> δZ^(r)_X  (dZ(A, G) = δZ^(G)_A, dZ(omega, W) = δZ^(W)_ω);  dZc(g_s) -> δZ_{g_s}
    dpar(y_u) -> δy_u (δλ, δm², δg_s, ...);  rc(name) -> r_name;  Qop(name) -> Q_name;  Nop(name) -> N_name
    name symbols with arbitrary sub- and superscripts: symbol(base, sub=..., sup=...) (Python, section 1)
    S0 (S_0)  Delta1 (Δ_1)  deltaS1 (δS^(1))  DeltaS1 (ΔS^(1))  DeltaS1ST (ΔS^(1)_ST)  Gamma (Γ)  Sren (S_ren)
    sop (𝐬, BRST operator)  bop (𝐛)  bdop (𝐛_d)  STop(x) (𝒮(x))  Kop (K)  Pphys (𝒫)  Evan (ℰ)  hc (h.c.)
    method notation (theory cards): Taylor(n) (T^(n))  deg(G)  Kf (K^(η))  Ke (K_η)  Rop (ℛ)  Rbarop (ℛ̄)  Dl(l) (Δ^[l])
    deltaS(l) (δS^(l))  DeltaS(l) (ΔS^(l))  Delta1r(r) (Δ_1^(r))  Ssum(r) (S^(r), sum over the multiplets)
    Eop(A, r) (E^(r)_A)  Etil(F, r) (Ẽ^(r)_F)  Nhat(A, r) (N̂^(r)_A)  Oop(A, r) (𝒪^(r)_A)  fder(F, X) (δF/δX)  sqrtn(x) (√x)
    a second argument of rc/Qop/Nop/Xc/Eop/... written up(x) is a superscript without parentheses:
    Qop(K, up(f)) -> Q^f_K

Fields (functions; the atom prints the field with indices suppressed, as in the report's matrix notation)::

    A(mu, a)  G(mu, A)  W(mu, I)  B(mu)                     gauge fields; Lorentz index lower, adjoint upper
    omega(a)  omegabar(a)  b(a)                             ghost ω^a, antighost ω̄^a, auxiliary field b^a
    rho(mu, a)  rhotilde(mu, a)  zsrc(a)                    BRST sources ρ^μ_a, ρ̃^μ_a, ζ_a (zeta is Symbolica's ζ)
    Q U D L E  Qbar Ubar Dbar Lbar Ebar                     Dirac fields Ψ_f (incl. the sterile partner) and Ψ̄_f
    Q(i) ... Qbar(i)                                        the same with explicit (lower) group indices
    Psi(f) Psibar(f) psi(f) psibar(f) chi(f) chibar(f)      generic multiplet f, physical and sterile components
    R(f)  Rbar(f)                                           fermion sources R_f, R̄_f
    H  Hdag  Htilde  Htildedag  H(i)  Hdag(i)               Higgs doublet, H†, H̃ = ε H*, H̃†
    sigma  sigmadag  sigma(i)  sigmadag(i)                  Higgs sources σ, σ†
    F(mu, nu, a)  Am(mu)  wm                                field strength F^a_{μν}, matrix-valued 𝒜_μ, 𝛚
    Hv  phi  h  G0  Gp  Gm  Z(mu)  Wp(mu)  Wm(mu)  Agam(mu) broken phase: H_v, φ, h, G⁰, G⁺, G⁻, Z_μ, W^±_μ, A^γ_μ
    u d e nu  ubar dbar ebar nubar  uprime dprime eprime    broken-phase Dirac fields (u, ū, u′ = sterile partner)
    field(name, idx...)                                     any other field, printed as its name

Derivatives and accents::

    pd(mu, X) -> ∂_μ X;  pdbar(mu, X), pdhat(mu, X);  pd2(X) -> ∂² X (pdbar2, pdhat2);  Dcov(mu) -> D_μ
    bar(X) -> X with an overline (four-dimensional part);  hat(X) -> X with a hat ((d-4)-dimensional part);
    tilde(X);  dagger(X) -> X†;  par(X) -> (X);  cc(X) -> X* (conj is Symbolica's own conjugation)
    Inside ord(...) a pd(...) that is not the last factor is printed in parentheses, (∂_μ X) Y: a derivative acts
    on the field it is written with only.

Lorentz and Dirac structures::

    g(mu, nu)  gbar(mu, nu)  ghat(mu, nu)                   metrics g^{μν}, ḡ^{μν}, ĝ^{μν} (upper by default)
    levi(mu, nu, rho, sigma)                                Levi-Civita symbol ε^{μνρσ} (curly ε; ϵ is the regulator)
    gam(mu)  gammabar(mu)  gammahat(mu)  gamma5  PL  PR     γ^μ, γ̄^μ, γ̂^μ, γ_5, P_L, P_R (gamma is Symbolica's Γ)
    Pc(f)  Pcbar(f)  sigmaT(mu, nu)                         P_{c_f}, P_{c̄_f}, σ^{μν}
    p(k) -> p_k;  p(k, mu) -> p_k^μ;  q, k, pp (p′) momenta  (any momentum expression m below, e.g. p(1) - p(2))
    vec(m, mu) -> m^μ;  vecbar(m, mu), vechat(m, mu)        e.g. vec(p(1) - p(2), mu) -> (p_1 − p_2)^μ
    slash(m)  slashbar(m)  slashhat(m)                      p̸, its barred part (bar above the slashed letter), hat
    dot(m1, m2)  dotbar  dothat   sq(m)  sqbar(m)  sqhat(m) scalar products and squares (p̄², p̂², ...)
    tr(x1, x2, ...)                                         Dirac trace of the ordered product
    mel(x, i, j) -> (x)_{ij}                                matrix element, e.g. mel(ord(T(I), T(J)), j, i)

Group structures::

    T(a) -> T^a (matrix);  T(a, i, j) -> (T^a)_{ij};  tau(I) -> τ^I;  Tg(a) -> 𝒯^a, Tg(a, f) -> 𝒯^a_f
    f(a, b, c) -> f^{abc};  fg(a, b, c) -> 𝔣^{abc};  dsym(a, b, c) -> d^{abc}
    delta(a, b) -> δ^{ab} (lo/up for δ_{ij}, δ^α_β);  levi(lo(i), lo(j)) -> ε_{ij};  levi(I, J, K) -> ε^{IJK}
    Tr(x1, x2, ...) -> Tr(x1 x2 ...) (trace over gauge indices);  comm(x, y) -> [x, y];  acomm(x, y) -> {x, y}
    Sinv(f, r) -> S_f^(r);  Cinv(f) -> C_f;  Cad(r) -> C_ad^(r);  Ytil -> Ỹ;  Ycal(f) -> 𝒴_f

Ordered products, sums, integrals and legs::

    ord(x1, ..., xn)       non-commutative product, printed in order (fields, spinor chains, group matrices)
    rat(p, q)              the rational p/q as a leading fraction
    rati(p, q)             the imaginary rational (p i)/q
    cf(x, rat(p, q))       a coefficient: the leading rational followed by the scalar product x
                           (M(...) rewrites every rational coefficient into rat/rati/cf, see display_form)
    withhc(x) -> x + h.c.;  intd(x) -> ∫ dᵈx x;  int4(x) -> ∫ d⁴x x;  sumover(i, x) -> Σ_i x
    plus(x, y) -> x + y, minus(x, y) -> x − y: sums printed in the written order (nest them for longer sums)
    leg(X, m) -> X(m), e.g. leg(G(mu, A), p(1)) -> G^A_μ(p_1);  out(m) marks an outgoing momentum: X(m, out)
    eqs(x, y) -> x = y;  defeq(x, y) -> x ≔ y;  to(x, y) -> x → y

Printing rules worth knowing. Symbolica orders a product as: number, symbols (in registration order), sums,
functions (in registration order); this is why a rational coefficient becomes the function ``cf`` and why a short
sum multiplied by a function factor is wrapped in ``par`` (so that ``f^{ACD} (...)`` keeps the report's order).
Symbolica omits a custom-printed fragment larger than about 4 kB (it prints "⋯"); the converters therefore never
wrap long sums (more than 12 terms or more than 3.5 kB of MathML) into ``par``. The MathML is normalised in one
place (``_typeset``): the power of a subscripted symbol is a single ``msubsup`` (g_s², p_1²).

Generic Symbolica structure (sums, products, powers, fractions, numbers) is printed by Symbolica itself. Unknown
symbols of the ``gdocs`` namespace print as their name; a symbol of any other namespace (e.g. an unconverted
``hepkit_renorm`` head) makes :func:`M` raise ``ValueError``.

Examples (source -> printed form)::

    -1/4*ord(F(mu, nu, a), F(up(mu), up(nu), a))   -> −1/4 F^a_{μν} F^{aμν}
    ord(Qbar, gammabar(mu), T(A), PL, Q, G(lo(mu), A)) -> Q̄ γ̄^μ T^A P_L Q G^A_μ
    lfe*11/3*g_s^2                                -> 1/(16π²) 1/ϵ 11/3 g_s²
    g(mu, nu)*vec(p(1) - p(2), rho)               -> g^{μν} (p_1 − p_2)^ρ
    ord(omega(A), Qbar, T(A), gammahat(mu), PR, pd(lo(mu), Q))

3. Converters from hepkit_renorm
================================

All converters return a Symbolica expression in the site notation; anything they cannot represent raises
:class:`ConversionError` (a ``ValueError``) with the reason. They never print package heads.

``operator(obj, *, model=None, strip=False)``
    A Lagrangian-density operator from a ``Monomial``, a ``Density`` or ``LocalFunctional`` (all monomials, summed),
    or a list of monomials. For each monomial: the scalar coefficient (numbers, couplings, ``d``, ``eta``, ``π``),
    the tensors that are not part of a chain (metrics between free labels, ``f``, ``d``, ``ε``, explicit generators),
    and the ordered product of the fields. Conventions of the display:

    * Spinor indices are never shown: every spinor line is a chain ``Ψ̄ Γ Ψ`` (row field, Dirac and group matrices,
      column field), closed lines are Dirac traces ``tr(...)``. Row fields are ``Ψ̄_f`` and ``R̄_f`` (the source
      ``src_f``), column fields ``Ψ_f`` and ``R_f`` (``src_fbar``).
    * Group indices of a fundamental space are suppressed when the space decomposes into chains that read left to
      right in the displayed order (``Q̄ γ̄^μ T^A P_L Q``, ``H† T^I H``, ``Q̄ H P_R D``) or closed traces
      ``Tr(T^A T^C)``; otherwise (e.g. the ``ε_{ij}`` of ``Q̄ H̃ P_R U``) that space is written with explicit lower
      indices ``i, j, k, l`` (SU(2)) and ``α, β, ...`` (colour).
    * Metrics are absorbed: ``ḡ^{μν} X_μ Y_ν`` is displayed as ``X̄_μ Ȳ^μ`` (bar on both objects), ``ĝ`` with hats.
    * The fields are reordered into the display order (antighosts, sources of gauge fields and ghosts, ghosts,
      fermion chains, scalar chains, gauge fields, auxiliary fields, other sources) and the Grassmann sign of the
      reordering is included in the coefficient: the displayed monomial equals the package monomial.
    * Letters: Lorentz μ ν ρ σ α β κ τ (μ ν α β κ τ when a BRST source is present, as in the report), SU(3)
      adjoint A C D E F, SU(2) adjoint I J K L M N, hypercharge components carry the superscript B
      (``ω^B``, ``b^B``), SU(2) structure constants are printed ``ε^{IJK}`` (``f^{IJK} = ε^{IJK}`` for ``T = τ/2``).
    ``strip=True`` removes the overall scalar factor; :func:`operator_parts` returns ``(scalar, operator)``.

``vertex(obj, legs=None, *, model=None, simplifier=None, momentum_conservation=False, outgoing=(), info=None,
loop=None)``
    A vertex function ``V`` (the Feynman rule is ``i V``). ``obj`` is a ``VertexRule`` (its value and legs are used
    as they are), or a ``LocalFunctional`` / ``StoreEntry`` coefficient together with ``legs`` (field names, in the
    Grassmann order of the legs; the value is ``functional.vertex_function(standard_legs(...))``, simplified with
    ``simplifier`` or ``model.simplifier(expand_projectors=False)``). Leg ``k`` has momentum ``p_k`` (incoming;
    ``outgoing`` lists legs whose ``p_k`` is outgoing, as in ``vertex_rule``); its free labels are printed with
    the letters of :func:`legs` (Lorentz μ ν ρ σ in leg order, adjoint letters per factor, fundamental i j k l).
    Matrix-valued rules (a Ψ̄ and a Ψ leg) are written as the matrix between them (``γ̄^μ T^A P_L``); scalar legs
    keep explicit indices (``(T^I)_{ij}``). Terms that differ only by a momentum are collected
    (``f^{ACD} g^{μν} (p_1 − p_2)^ρ``) and a common factor is taken out. ``info={}`` receives ``explicit`` (the
    fundamental spaces written with indices), ``fields`` and ``legs`` (leg labels consistent with the value; pass
    the union of ``explicit`` over the parts of one rule to :func:`legs` when several parts share leg labels).
    ``loop="pole"`` / ``"finite"`` writes the explicit ``1/π²`` of an order-store coefficient as ``lfe`` / ``lf``
    in front (the same option exists for :func:`operator`). Example::

        info = {}
        V = vertex(store.entry("S0").coefficient(0), ("G", "Q", "Qbar"), info=info)   # -g_s γ̄^μ T^A P_L
        info["legs"]                                                                 # G^A_μ(p_1), Q(p_2), Q̄(p_3)
        vertex(store.entry("deltaS1").coefficient(-1), ("G", "Q", "Qbar"), loop="pole")

``coefficient(value, *, loop=None, sign=1)``
    A scalar coefficient: a Symbolica expression (package or site symbols), a string (Symbolica or Python syntax,
    names ``g_s g gp y_u y_d y_e lambda m2 xi``), an ``int`` or a ``Fraction``. ``loop="pole"`` multiplies by
    ``1/(16π²)`` and ``1/ϵ`` written as ``lfe``; ``loop="finite"`` by ``lf``; ``loop="auto"`` recognises an explicit
    ``1/π²`` (as in the order store) and writes it as ``lf`` (``lfe`` when ``pole=True``). Rationals are leading
    fractions (``rat``), as in the report: ``-4/3*g_s^2`` prints ``−4/3 g_s²``.

A field table of another basis (``gdocs.broken_fields.BT``, the broken phase) names its fields itself: if the table
    has a method ``site_field(name, idx)``, :func:`field` uses its result (``None``: the rules above).

``legs(fields, *, table=None, outgoing=(), explicit=None)`` -> list of expressions ``leg(X, p(k))`` such as
    ``G^A_μ(p_1)``, with the letters used by :func:`vertex`; fermion legs carry group indices only for the spaces in
    ``explicit``; ``leg(name, k)`` for one leg. ``field(name, ...)`` gives the site
    symbol of a package field (``gh_G`` -> ``omega(A)``).

4. The interactive payload (``view``)
=====================================

With ``interactive=True`` the math object carries ``view = {"v": 1, "pages": [P0, ...], "previews": [...],
"answers": {key: [...]}, "bytes": n, "complete": bool}``. ``P0`` is ``Expression.mathml_viewer(page_size).page()``;
an answer is ``["p", page index, new state]`` (navigation) or ``["v", preview index, context]`` (popup), the key
is :func:`answer_key` of the state and the request. Without ``page_size`` the page sizes automatic, 25, 12 and 6
are tried and the one with the smallest initial page is used; if that page alone takes more than 70 % of the
budget, ``view`` is ``None`` and the expression is shown statically (the case of the largest counterterm vertex
functions, about 25 kB of MathML with 80 kB of layout alternatives). The answers are Symbolica's own ``handle_message`` responses to the requests a reader can trigger from
the initial state, found by a breadth-first exploration of the pager's protocol (hover and pinned previews of every
omission marker in every responsive layout, inspection and ``back``, pagination inside previews and of inspected
parts, nested previews, expansion and collapse), deduplicated, until the payload reaches ``budget`` characters of
JSON (``VIEW_BUDGET`` = 30000 by default; a long expression worth exploring can be given more). ``complete`` tells
whether the exploration finished within the budget (most interactive views of long expressions are partial: the
answers nearest to the initial view are kept, deeper ones are dropped).
The static adapter of ``assets/viewer/mathview.js`` answers the pager from this table; a request without an answer
is never sent to a server: its control is disabled (marked "not available in the static page") before it can be
used.

5. Assets
=========

``write_assets(site_dir)`` writes ``assets/viewer/symbolica_pager.js`` (exactly ``MathMLViewer.javascript()``) and
``assets/viewer/STIXTwoMath-Regular.woff2`` (the font shipped with the Symbolica wheel). ``mathview.js`` (the
static adapter, ``renderMath`` / ``ready``) is source code kept in the site.
"""

from __future__ import annotations

import html as _html
import itertools
import json
import re
import shutil
from collections.abc import Iterable, Sequence
from fractions import Fraction
from pathlib import Path

from symbolica import AtomType, E, Expression, MathMLViewer, PrintMode, S

__all__ = [
    "ConversionError",
    "M",
    "VIEW_BUDGET",
    "coefficient",
    "convert",
    "field",
    "leg",
    "legs",
    "operator",
    "operator_parts",
    "parse",
    "prose",
    "registry",
    "symbol",
    "vertex",
    "write_assets",
]

NS = "gdocs"
VIEW_BUDGET = 30_000
MATHML_NS = "http://www.w3.org/1998/Math/MathML"


class ConversionError(ValueError):
    """An object that the site notation cannot represent (raised instead of printing garbage)."""


# ================================================================================================ MathML fragments
class _F:
    """A printed fragment: MathML (``ml``) and plain Unicode text (``tx``). ``atom``: needs no parentheses."""

    __slots__ = ("ml", "tx", "atom")

    def __init__(self, ml: str, tx: str, atom: bool = True) -> None:
        self.ml = ml
        self.tx = tx
        self.atom = atom


def _esc(s: str) -> str:
    return _html.escape(s, quote=False)


_UPRIGHT_GREEK = set("ΓΔΘΛΞΠΣΥΦΨΩ")


def _mi(s: str, normal: bool = False) -> _F:
    attr = ' mathvariant="normal"' if (normal or s in _UPRIGHT_GREEK) and len(s) == 1 else ""
    return _F(f"<mi{attr}>{_esc(s)}</mi>", s)


def _mn(s: str) -> _F:
    return _F(f"<mn>{_esc(s)}</mn>", s)


def _mo(s: str, lspace: str | None = None, rspace: str | None = None, stretchy: bool | None = None) -> _F:
    attrs = ""
    if lspace is not None:
        attrs += f' lspace="{lspace}"'
    if rspace is not None:
        attrs += f' rspace="{rspace}"'
    if stretchy is not None:
        attrs += f' stretchy="{str(stretchy).lower()}"'
    return _F(f"<mo{attrs}>{_esc(s)}</mo>", s)


def _row(*parts: _F, sep: str = "") -> _F:
    parts = [p for p in parts if p is not None]
    if len(parts) == 1:
        return parts[0]
    return _F("<mrow>" + "".join(p.ml for p in parts) + "</mrow>", sep.join(p.tx for p in parts), atom=False)


def _join(parts: Sequence[_F]) -> _F:
    """Index list ``μν`` (no separator, text joined)."""
    if len(parts) == 1:
        return parts[0]
    return _F("<mrow>" + "".join(p.ml for p in parts) + "</mrow>", "".join(p.tx for p in parts), atom=False)


def _tx_script(t: str, mark: str) -> str:
    return t if len(t) == 1 else "{" + t + "}"


def _sub(base: _F, sub: _F | None) -> _F:
    if sub is None:
        return base
    return _F(f"<msub>{_wrap(base).ml}{_wrap(sub).ml}</msub>", f"{base.tx}_{_tx_script(sub.tx, '_')}")


def _sup(base: _F, sup: _F | None) -> _F:
    if sup is None:
        return base
    return _F(f"<msup>{_wrap(base).ml}{_wrap(sup).ml}</msup>", f"{base.tx}^{_tx_script(sup.tx, '^')}")


def _subsup(base: _F, sub: _F | None, sup: _F | None) -> _F:
    if sub is None:
        return _sup(base, sup)
    if sup is None:
        return _sub(base, sub)
    return _F(f"<msubsup>{_wrap(base).ml}{_wrap(sub).ml}{_wrap(sup).ml}</msubsup>",
              f"{base.tx}^{_tx_script(sup.tx, '^')}_{_tx_script(sub.tx, '_')}")


def _wrap(f: _F) -> _F:
    """A single MathML element (scripts need exactly one child per position)."""
    ml = f.ml
    if _single_element(ml):
        return f
    return _F(f"<mrow>{ml}</mrow>", f.tx, f.atom)


def _single_element(ml: str) -> bool:
    m = re.match(r"<([a-z]+)[ >]", ml)
    if not m:
        return False
    tag = m.group(1)
    depth = 0
    for t in re.finditer(r"<(/?)([a-z]+)[^>]*?(/?)>", ml):
        if t.group(3):
            if depth == 0 and t.end() != len(ml):
                return False
            continue
        depth += -1 if t.group(1) else 1
        if depth == 0:
            return t.end() == len(ml) and t.group(2) == tag
    return False


_ACCENTS = {"bar": ("‾", "̅"), "hat": ("ˆ", "̂"), "tilde": ("˜", "̃")}


def _accent(base: _F, kind: str | None) -> _F:
    if not kind:
        return base
    mark, comb = _ACCENTS[kind]
    return _F(f'<mover accent="true">{_wrap(base).ml}<mo>{mark}</mo></mover>', f"{base.tx}{comb}")


def _frac(num: _F, den: _F) -> _F:
    return _F(f"<mfrac>{_wrap(num).ml}{_wrap(den).ml}</mfrac>", f"{_ptx(num)}/{_ptx(den)}")


def _ptx(f: _F) -> str:
    return f.tx if f.atom else f"({f.tx})"


def _paren(f: _F, o: str = "(", c: str = ")") -> _F:
    tall = any(t in f.ml for t in ("<mfrac", "<munder", "<msqrt", "<mtable", "∫", "∑"))
    st = "" if tall else ' stretchy="false"'
    return _F(f"<mrow><mo{st}>{_esc(o)}</mo>{f.ml}<mo{st}>{_esc(c)}</mo></mrow>", f"{o}{f.tx}{c}")


def _slashed(base: _F) -> _F:
    """Feynman slash: a solidus with zero advance drawn over the letter (scripts attach after the letter)."""
    return _F(f'<mrow><mpadded width="0" lspace="0.1em"><mo lspace="0" rspace="0">/</mo></mpadded>{_wrap(base).ml}</mrow>',
              f"{base.tx}̸")


def _space(w: str = "0.167em") -> _F:
    return _F(f'<mspace width="{w}"/>', " ")


_PRIME = _F('<mo lspace="0" rspace="0">′</mo>', "′")


# ================================================================================================ registry
_SYMBOLS: dict[str, Expression] = {}
_PRINTERS: dict[str, object] = {}
_DOCS: list[tuple[str, str, str]] = []  # (name, example source, description)

_GREEK = {
    "mu": "μ", "nu": "ν", "rho": "ρ", "sigma": "σ", "alpha": "α", "beta": "β", "kappa": "κ", "tau": "τ",
    "gamma": "γ", "delta": "δ", "iota": "ι", "chi": "χ", "omega": "ω", "phi": "φ", "psi": "ψ", "zeta": "ζ",
    "eta": "η", "xi": "ξ", "theta": "θ", "lambda": "λ", "pi": "π", "upsilon": "υ", "Delta": "Δ", "Gamma": "Γ",
    "Sigma": "Σ", "Omega": "Ω", "Lambda": "Λ", "Phi": "Φ", "Psi": "Ψ", "Theta": "Θ", "epsilon": "ε",
}


def _short(e: Expression) -> str:
    """Name of a symbol or function head without namespace."""
    return e.get_name().split("::")[-1]


def _callback(e: Expression, mode: PrintMode, **kw) -> str | None:
    if mode == PrintMode.MathML:
        return _frag(e).ml
    if mode == PrintMode.Symbolica:
        cpm = kw.get("custom_print_mode") or {}
        if "gdocs_text" in cpm:
            return _frag(e).tx
    return None


def _reg(name: str, printer, doc: str = "", example: str | None = None, **attrs) -> Expression:
    """Register a symbol (or function head) of the site notation with its printer ``printer(e) -> _F``."""
    if name in _SYMBOLS:
        raise RuntimeError(f"duplicate registry entry {name}")
    sym = S(f"{NS}::{name}", print=_callback, **attrs)
    _SYMBOLS[name] = sym
    _PRINTERS[name] = printer
    _DOCS.append((name, example if example is not None else name, doc))
    return sym


def registry() -> list[tuple[str, str, str]]:
    """``[(name, example source, description)]`` of every registered symbol, in registration (display) order."""
    return list(_DOCS)


def _args(e: Expression) -> list[Expression]:
    return list(e) if e.get_type() == AtomType.Fn else []


def _is_fn(e: Expression, name: str | None = None) -> bool:
    if e.get_type() != AtomType.Fn:
        return False
    return name is None or (e.get_name() == f"{NS}::{name}")


def _is_sym(e: Expression, name: str) -> bool:
    return e.get_type() == AtomType.Var and e.get_name() == f"{NS}::{name}"


# ------------------------------------------------------------------------------------------------ generic printing
_STRIP_ATTR = re.compile(r' (?:data-symbolica-[a-z-]+="[^"]*"|style="[^"]*")')
_STRIP_TABLE_OPEN = re.compile(r"<mtable data-symbolica-(?:product|sum)[^>]*><mtr><mtd[^>]*>")
_STRIP_TABLE_CLOSE = re.compile(r"</mtd></mtr></mtable>")


_SUBSUP = re.compile(r"<msup><mrow><msub>((?:(?!</?msub>).)*)</msub></mrow>(<mrow><mn>[0-9]+</mn></mrow>|<mn>[0-9]+"
                     r"</mn>)</msup>")


def _typeset(ml: str) -> str:
    """Typographic normalisation of printed MathML: the power of a subscripted symbol (g_s², p_1²) is one
    ``msubsup`` instead of a superscript after the subscript."""
    return _SUBSUP.sub(r"<msubsup>\1\2</msubsup>", ml)


def _nested_ml(e: Expression) -> str:
    """Symbolica's MathML of a sub-expression, flattened for use inside a custom fragment (no pager metadata,
    no wrapping tables)."""
    ml = _typeset(e.format(mode=PrintMode.MathML))
    if "<mtable data-symbolica" in ml:
        # only single-row tables are produced for short nested expressions; multi-row ones keep their rows
        ml = _STRIP_TABLE_OPEN.sub("", ml)
        ml = _STRIP_TABLE_CLOSE.sub("", ml)
        ml = re.sub(r"</mtd></mtr><mtr><mtd[^>]*>", "", ml)  # row breaks of multi-row layouts (cells may carry a style)
    return _STRIP_ATTR.sub("", ml)


def _text(e: Expression) -> str:
    return e.format(mode=PrintMode.Symbolica, custom_print_mode={"gdocs_text": 1}, max_terms=None,
                    max_line_length=None, color_top_level_sum=False, color_builtin_symbols=False,
                    bracket_level_colors=None, multiplication_operator=" ", hide_namespace=NS,
                    num_exp_as_superscript=True, print_ring=False)


def _source(e: Expression) -> str:
    return e.format(mode=PrintMode.Symbolica, max_terms=None, max_line_length=None, color_top_level_sum=False,
                    color_builtin_symbols=False, bracket_level_colors=None, multiplication_operator="*",
                    hide_namespace=NS, num_exp_as_superscript=False, print_ring=False,
                    double_star_for_exponentiation=False)


def _generic(e: Expression) -> _F:
    """Any expression: registry objects by their printers, the rest by Symbolica."""
    t = e.get_type()
    if t in (AtomType.Var, AtomType.Fn):
        name = e.get_name()
        if name.startswith(NS + "::"):
            short = name[len(NS) + 2:]
            pr = _PRINTERS.get(short)
            if pr is not None:
                return pr(e)
            if t == AtomType.Var:
                return _plain_symbol(short)
            return _function_like(_plain_symbol(short), _args(e))
        if "::" in name and not name.startswith("symbolica::") and not name.startswith("python::"):
            raise ConversionError(f"symbol {name} is not part of the site notation")
        if t == AtomType.Var:
            short = name.split("::")[-1]
            if short in ("𝜋", "π"):
                return _mi("π")
            if short in ("𝑖",):
                return _mi("i")
            return _plain_symbol(short)
    atom = t in (AtomType.Num, AtomType.Var) and not (t == AtomType.Num and ("/" in e.format_plain() or e.format_plain().startswith("-")))
    return _F(_nested_ml(e), _text(e), atom=atom)


def _plain_symbol(short: str) -> _F:
    if short in _GREEK:
        return _mi(_GREEK[short])
    if re.fullmatch(r"[A-Za-z]", short):
        return _mi(short)
    m = re.fullmatch(r"([A-Za-z]+)_([A-Za-z0-9]+)", short)
    if m:
        base = _GREEK.get(m.group(1), m.group(1))
        return _sub(_mi(base), _mi(m.group(2)) if not m.group(2).isdigit() else _mn(m.group(2)))
    return _mi(short)


def _function_like(head: _F, args: Sequence[Expression]) -> _F:
    inner = _sepjoin([_frag(a) for a in args], ",")
    return _row(head, _paren(inner))


def _sepjoin(parts: Sequence[_F], sep: str) -> _F:
    if not parts:
        return _F("<mrow></mrow>", "")
    out = []
    for i, p in enumerate(parts):
        if i:
            out.append(_mo(sep))
        out.append(p)
    return _F("<mrow>" + "".join(p.ml for p in out) + "</mrow>", (sep + " ").join(p.tx for p in parts), atom=False)


def _frag(e: Expression) -> _F:
    return _generic(e)


# ------------------------------------------------------------------------------------------------ indices
def _index(e: Expression) -> tuple[_F, str | None]:
    """An index argument: its fragment and forced position (``"lo"``, ``"up"`` or ``None``)."""
    if _is_fn(e, "lo") or _is_fn(e, "up"):
        inner = _args(e)[0]
        return _index(inner)[0], ("lo" if _is_fn(e, "lo") else "up")
    if _is_fn(e, "prime"):
        f, pos = _index(_args(e)[0])
        return _F(f"<mrow>{f.ml}{_PRIME.ml}</mrow>", f.tx + "′"), pos
    t = e.get_type()
    if t == AtomType.Num:
        return _mn(e.format_plain()), None
    if t == AtomType.Var:
        short = _short(e)
        if e.get_name().startswith(NS + "::"):
            pr = _PRINTERS.get(short)
            if pr is not None and short not in _INDEX_LIKE:
                return pr(e), None
        return _plain_symbol(short), None
    return _frag(e), None


_INDEX_LIKE: set[str] = set()


def _scripted(base: _F, idx: Sequence[Expression], defaults: Sequence[str], accent: str | None = None) -> _F:
    """``base`` with the indices ``idx`` placed at their default positions (``"lo"``/``"up"`` per argument, the
    last entry repeats) unless forced by ``lo``/``up``."""
    lows, ups = [], []
    for k, a in enumerate(idx):
        f, pos = _index(a)
        d = defaults[min(k, len(defaults) - 1)] if defaults else "lo"
        (lows if (pos or d) == "lo" else ups).append(f)
    b = _accent(base, accent)
    return _subsup(b, _join(lows) if lows else None, _join(ups) if ups else None)


# ------------------------------------------------------------------------------------------------ printers
def _p_const(ml: _F):
    def pr(e):
        if e.get_type() == AtomType.Fn:
            return _function_like(ml, _args(e))
        return ml
    return pr


def _p_sub(base: str, sub: str, sub_normal: bool = False, base_normal: bool = False):
    f = _sub(_mi(base, normal=base_normal), _mi(sub, normal=sub_normal) if not sub.isdigit() else _mn(sub))
    return _p_const(f)


def _p_indexed(base: _F, defaults: Sequence[str], accent: str | None = None, atom_base: _F | None = None):
    """A head printed as ``base`` with its arguments as indices (atom: ``atom_base`` or ``base``)."""
    def pr(e):
        if e.get_type() != AtomType.Fn:
            return _accent(atom_base or base, accent)
        return _scripted(base, _args(e), defaults, accent)
    return pr


def _p_matrix_element(base: _F, defaults=("up",)):
    """``T(a)`` -> T^a; ``T(a, i, j)`` -> (T^a)_{ij}."""
    def pr(e):
        if e.get_type() != AtomType.Fn:
            return base
        a = _args(e)
        if len(a) <= 1:
            return _scripted(base, a, defaults)
        head = _scripted(base, a[:1], defaults)
        lows = [_index(x)[0] for x in a[1:]]
        return _sub(_paren(head), _join(lows))
    return pr


def _p_accent_wrapper(kind: str):
    def pr(e):
        a = _args(e)
        if len(a) != 1:
            raise ConversionError(f"{kind}(x) takes one argument")
        inner = a[0]
        if inner.get_type() in (AtomType.Var, AtomType.Fn) and inner.get_name().startswith(NS + "::"):
            short = _short(inner)
            acc = _ACCENTABLE.get(short)
            if acc is not None:
                return acc(inner, kind)
        f = _frag(inner)
        return _accent(f if f.atom else _paren(f), kind)
    return pr


_ACCENTABLE: dict[str, object] = {}


def _field(base: _F, defaults: Sequence[str], conj: str | None = None, deco: _F | None = None):
    """A field head: ``X`` (indices suppressed) or ``X(idx...)``; ``conj``: ``"bar"`` (overline on the letter) or
    ``"dag"`` (dagger superscript, after the indices); ``deco``: a fixed superscript (W^+, A^γ, G^0)."""
    def build(e, accent=None):
        idx = _args(e) if e.get_type() == AtomType.Fn else []
        b = base
        if conj == "bar":
            b = _accent(b, "bar")
        if accent:
            b = _accent(b, accent)
        lows, ups, raised = [], [], []
        for k, a in enumerate(idx):
            f, pos = _index(a)
            d = defaults[min(k, len(defaults) - 1)] if defaults else "lo"
            if (pos or d) == "lo":
                lows.append(f)
            else:
                (ups if d == "up" else raised).append(f)
        # several superscripts: the group index (default up) before raised Lorentz indices (G^{Aμ}, F^{aμν})
        ups += raised
        if deco is not None:
            ups.insert(0, deco)
        if conj == "dag":
            ups.append(_mo("†"))
        return _subsup(b, _join(lows) if lows else None, _join(ups) if ups else None)
    return build


def _reg_field(name: str, base: _F, defaults: Sequence[str], conj: str | None = None, doc: str = "",
               example: str | None = None, deco: _F | None = None) -> Expression:
    build = _field(base, defaults, conj, deco)
    _ACCENTABLE[name] = lambda e, kind: build(e, kind)
    return _reg(name, lambda e: build(e), doc, example)


def _arg_frag(e: Expression) -> _F:
    f = _frag(e)
    return f


def _momentum(e: Expression) -> _F:
    """A momentum argument: ``p(1)`` -> p_1, a symbol, or a combination (parenthesised)."""
    return _frag(e)


def _p_vec(accent: str | None):
    def pr(e):
        a = _args(e)
        if len(a) != 2:
            raise ConversionError("vec(m, mu) takes a momentum and an index")
        src = a[0]
        if accent and not _momentum(src).atom:
            # (p̄_1 − p̄_2)^μ: the accent goes on every momentum of the combination
            src = src.replace(_SYMBOLS["p"](S("x_")), _SYMBOLS["p" + accent](S("x_")))
            return _scripted_existing(_paren(_momentum(src)), a[1:], ("up",))
        m = _momentum(src)
        if accent:
            m = _accent(m, accent) if not _is_fn(src, "p") else _p_momentum_acc(accent)(src)
        elif not m.atom:
            m = _paren(m)
        return _scripted_existing(m, a[1:], ("up",))
    return pr


def _scripted_existing(base: _F, idx: Sequence[Expression], defaults: Sequence[str]) -> _F:
    """Attach indices to an already built fragment; a fragment that is itself scripted (``p_1``) gets the
    superscript on the same base (``p_1^μ``)."""
    lows, ups = [], []
    for k, a in enumerate(idx):
        f, pos = _index(a)
        d = defaults[min(k, len(defaults) - 1)]
        (lows if (pos or d) == "lo" else ups).append(f)
    m = re.fullmatch(r"<msub>(.*)</msub>", base.ml)
    if m and ups and not lows and _split_two(m.group(1)) is not None:
        b, s = _split_two(m.group(1))
        sup = _join(ups)
        tx_base, _, tx_sub = base.tx.partition("_")
        return _F(f"<msubsup>{b}{s}{_wrap(sup).ml}</msubsup>", f"{tx_base}^{_tx_script(sup.tx, '^')}_{tx_sub}")
    return _subsup(base, _join(lows) if lows else None, _join(ups) if ups else None)


def _split_two(ml: str) -> tuple[str, str] | None:
    """Split the content of a two-child element into its children."""
    depth = 0
    for t in re.finditer(r"<(/?)([a-z]+)[^>]*?(/?)>", ml):
        if t.group(3):
            if depth == 0:
                return ml[: t.end()], ml[t.end():]
            continue
        depth += -1 if t.group(1) else 1
        if depth == 0:
            return ml[: t.end()], ml[t.end():]
    return None


def _p_slash(accent: str | None):
    def pr(e):
        a = _args(e)
        if len(a) != 1:
            raise ConversionError("slash(m) takes one momentum")
        m = a[0]
        if _is_fn(m, "p") and len(_args(m)) == 1:
            k = _index(_args(m)[0])[0]
            s = _sub(_accent(_slashed(_mi("p")), accent), k)
            return s
        f = _frag(m)
        if f.atom:
            return _accent(_slashed(f), accent)
        raise ConversionError("slash(m) takes a single momentum (write a sum of slashed momenta instead)")
    return pr


def _p_dot(accent: str | None):
    def pr(e):
        a = _args(e)
        if len(a) != 2:
            raise ConversionError("dot(m1, m2) takes two momenta")
        x, y = (_accented_momentum(v, accent) for v in a)
        return _row(x, _mo("⋅"), y)
    return pr


def _accented_momentum(m: Expression, accent: str | None) -> _F:
    """A momentum with an accent on its letter: p̄_1 (not an overline over p_1)."""
    if _is_fn(m, "p") and len(_args(m)) == 1:
        return _sub(_accent(_mi("p"), accent), _index(_args(m)[0])[0])
    f = _frag(m)
    f = f if f.atom else _paren(f)
    return _accent(f, accent) if accent else f


def _p_sq(accent: str | None):
    def pr(e):
        a = _args(e)
        if len(a) != 1:
            raise ConversionError("sq(m) takes one momentum")
        m = a[0]
        if _is_fn(m, "p") and len(_args(m)) == 1:
            k = _index(_args(m)[0])[0]
            return _subsup(_accent(_mi("p"), accent), k, _mn("2"))
        f = _frag(m)
        f = f if f.atom else _paren(f)
        return _sup(_accent(f, accent), _mn("2"))
    return pr


def _p_momentum(e):
    """``p(k)`` -> p_k, ``p(k, mu)`` -> p_k^μ."""
    if e.get_type() != AtomType.Fn:
        return _mi("p")
    a = _args(e)
    k = _index(a[0])[0]
    if len(a) == 1:
        return _sub(_mi("p"), k)
    return _scripted_existing(_sub(_mi("p"), k), a[1:], ("up",))


def _p_momentum_acc(accent):
    def pr(e):
        a = _args(e)
        if not a:
            return _accent(_mi("p"), accent)
        k = _index(a[0])[0]
        base = _sub(_accent(_mi("p"), accent), k)
        if len(a) == 1:
            return base
        return _scripted_existing(base, a[1:], ("up",))
    return pr


def _p_der(accent: str | None, square: bool = False):
    """``pd(mu, X)`` -> ∂_μ X, ``pd2(X)`` -> ∂² X."""
    def pr(e):
        a = _args(e)
        if e.get_type() != AtomType.Fn:
            return _accent(_mi("∂"), accent)
        d = _accent(_mi("∂"), accent)
        if square:
            if len(a) != 1:
                raise ConversionError("pd2(X) takes one argument")
            op = _sup(d, _mn("2"))
            target = a[0]
        else:
            if len(a) == 1:
                return _scripted(d, a, ("lo",))
            if len(a) != 2:
                raise ConversionError("pd(mu, X) takes an index and a field")
            op = _scripted(d, a[:1], ("lo",))
            target = a[1]
        t = _frag(target)
        if not t.atom and not _is_derivative(target):
            t = _paren(t)
        return _F("<mrow>" + op.ml + t.ml + "</mrow>", op.tx + " " + t.tx, atom=False)
    return pr


_DER_HEADS = {"pd", "pdbar", "pdhat", "pd2", "pdbar2", "pdhat2"}


def _is_derivative(e: Expression) -> bool:
    return e.get_type() == AtomType.Fn and e.get_name().startswith(NS + "::") and _short(e) in _DER_HEADS


def _p_ord(e):
    """Ordered product: arguments in order; a derivative factor followed by others is parenthesised."""
    a = _args(e)
    parts = []
    for k, x in enumerate(a):
        f = _frag(x)
        if _is_derivative(x) and k < len(a) - 1:
            f = _paren(f)
        elif not f.atom and x.get_type() == AtomType.Add:
            f = _paren(f)
        parts.append(f)
    if not parts:
        return _mn("1")
    ml = "<mrow>" + "".join(p.ml for p in parts) + "</mrow>" if len(parts) > 1 else parts[0].ml
    return _F(ml, " ".join(p.tx for p in parts), atom=len(parts) == 1 and parts[0].atom)


def _p_trace(name: str):
    def pr(e):
        a = _args(e)
        inner = _p_ord(e) if a else _mn("1")
        op = _F(f"<mi>{name}</mi>", name)
        return _F("<mrow>" + op.ml + _paren(inner).ml + "</mrow>", f"{name}({inner.tx})", atom=True)
    return pr


def _p_mel(e):
    a = _args(e)
    if len(a) != 3:
        raise ConversionError("mel(x, i, j) takes a matrix and two indices")
    x = _frag(a[0])
    return _sub(_paren(x), _join([_index(a[1])[0], _index(a[2])[0]]))


def _p_brackets(o: str, c: str):
    def pr(e):
        a = _args(e)
        return _paren(_sepjoin([_frag(x) for x in a], ","), o, c)
    return pr


def _p_wrap1(fn):
    def pr(e):
        a = _args(e)
        if len(a) != 1:
            raise ConversionError(f"{_short(e)}(x) takes one argument")
        return fn(_frag(a[0]), a[0])
    return pr


def _p_rat(imag: bool):
    def pr(e):
        a = _args(e)
        if len(a) != 2:
            raise ConversionError("rat(p, q) takes two integers")
        p, q = (x.format_plain() for x in a)
        num = _row(_mn(p), _mi("i")) if imag else _mn(p)
        if imag and p == "1":
            num = _mi("i")
        if q == "1":
            return num
        return _frac(num, _mn(q))
    return pr


def _p_cf(e):
    """``cf(rest, r)``: the leading rational ``r`` followed by the scalar product ``rest``."""
    a = _args(e)
    if len(a) != 2:
        raise ConversionError("cf(x, rat(p, q)) takes two arguments")
    r = _frag(a[1])
    rest = _frag(a[0])
    if a[0].get_type() == AtomType.Add:
        rest = _paren(rest)
    return _F("<mrow>" + r.ml + '<mspace width="0.1667em"/>' + rest.ml + "</mrow>", r.tx + " " + rest.tx,
              atom=False)


def _p_leg(e):
    a = _args(e)
    if len(a) < 2:
        raise ConversionError("leg(X, m) takes a field and a momentum")
    x = _frag(a[0])
    inner = [_frag(m) for m in a[1:]]
    return _row(x, _paren(_sepjoin(inner, ",")))


def _p_int(measure: str):
    def pr(e):
        a = _args(e)
        meas = _row(_mo("∫"), _sup(_mi("d"), _mi(measure) if not measure.isdigit() else _mn(measure)),
                    _mi("x"))
        if not a:
            return meas
        body = _frag(a[0])
        return _row(meas, _space(), body if body.atom else _paren(body))
    return pr


def _p_sumover(e):
    a = _args(e)
    if len(a) != 2:
        raise ConversionError("sumover(i, x) takes an index and a summand")
    idx = _index(a[0])[0]
    body = _frag(a[1])
    op = _F(f"<munder><mo>∑</mo>{_wrap(idx).ml}</munder>", f"Σ_{idx.tx}")
    return _row(op, body if body.atom or a[1].get_type() != AtomType.Add else _paren(body))


def _p_withhc(e):
    a = _args(e)
    body = _frag(a[0]) if a else _mn("0")
    return _row(body, _mo("+"), _F('<mtext>h.c.</mtext>', "h.c."))


def _p_binop(sym: str):
    def pr(e):
        a = _args(e)
        if len(a) != 2:
            raise ConversionError(f"{_short(e)}(x, y) takes two arguments")
        return _row(_frag(a[0]), _mo(sym), _frag(a[1]))
    return pr


def _p_osum(sym: str):
    def pr(e):
        a = _args(e)
        if len(a) != 2:
            raise ConversionError(f"{_short(e)}(x, y) takes two arguments")
        right = _frag(a[1])
        if sym != "+" and a[1].get_type() == AtomType.Add:
            right = _paren(right)
        return _row(_frag(a[0]), _mo(sym), right)
    return pr


def _p_dZ(e):
    """dZ(X) -> δZ_X ; dZ(X, r) -> δZ^(r)_X."""
    a = _args(e)
    base = _row(_mi("δ"), _mi("Z"))
    if not a:
        return base
    sub = _frag(a[0])
    sup = _paren(_frag(a[1])) if len(a) > 1 else None
    return _F(f"<mrow><mi>δ</mi>{_subsup(_mi('Z'), sub, sup).ml}</mrow>",
              "δ" + _subsup(_mi("Z"), sub, sup).tx, atom=True)


def _p_dpar(e):
    a = _args(e)
    if len(a) != 1:
        raise ConversionError("dpar(x) takes one parameter")
    x = _frag(a[0])
    return _F(f"<mrow><mi>δ</mi>{x.ml}</mrow>", "δ" + x.tx, atom=True)


def _p_named(base: str, accent: str | None = None):
    """``rc(x)`` -> r_x etc. (the argument is printed as the subscript). A second argument is the superscript in
    parentheses (``Qop(T, G)`` -> Q^(G)_T), or without them when written ``up(x)`` (``Qop(K, up(f))`` -> Q^f_K).
    ``accent`` puts an accent on the letter only (Ẽ^(r)_F, N̂^(r)_A)."""
    def pr(e):
        a = _args(e)
        b = _accent(_mi(base), accent) if accent else _mi(base)
        if e.get_type() != AtomType.Fn:
            return b
        if len(a) == 1:
            return _sub(b, _index(a[0])[0])
        sup = _index(a[1])[0] if _is_fn(a[1], "up") else _paren(_index(a[1])[0])
        return _subsup(b, _index(a[0])[0], sup)
    return pr


def _p_inv(base: str, sup_paren: bool = True):
    """Group invariants: Sinv(f, r) -> S_f^(r), Cinv(f) -> C_f, Cad(r) -> C_ad^(r)."""
    def pr(e):
        a = _args(e)
        if base == "Cad":
            sub = _mi("ad", normal=True)
            sup = _paren(_index(a[0])[0]) if a else None
            return _subsup(_mi("C"), sub, sup)
        sub = _index(a[0])[0] if a else None
        sup = _paren(_index(a[1])[0]) if len(a) > 1 else None
        return _subsup(_mi(base), sub, sup)
    return pr


# ------------------------------------------------------------------------------------------------ registration
# The order of registration is the display order of factors in a product.
def _register_all() -> None:
    R = _reg
    # leading coefficients
    R("rat", _p_rat(False), "rational p/q printed as a leading fraction", "rat(11, 3)")
    R("rati", _p_rat(True), "imaginary rational (p i)/q", "rati(2, 3)")
    R("cf", _p_cf, "coefficient cf(x, rat(p, q)): prints the leading rational, then x (created by M)",
      "cf(g_s^2, rat(11, 3))")
    # loop factors and regulator
    lf = _frac(_mn("1"), _row(_mn("16"), _sup(_mi("π"), _mn("2"))))
    R("lf", _p_const(_F("<mrow>" + lf.ml + '<mspace width="0.1667em"/></mrow>', lf.tx, atom=True)),
      "one-loop factor 1/(16π²), followed by a thin space")
    # 1/ϵ is a stacked fraction (as in the report, `lf 1/ep`), followed by a thin space: a slash form '1/ϵ g²'
    # would read as 1/(ϵ g²)
    R("lfe", _p_const(_F("<mrow>" + lf.ml + '<mspace width="0.1667em"/><mfrac><mn>1</mn><mi>ϵ</mi></mfrac>'
                         '<mspace width="0.1667em"/></mrow>', "1/(16π²) 1/ϵ", atom=True)),
      "pole factor 1/(16π²) 1/ϵ (as written in the report: two stacked fractions)")
    # couplings and parameters (display order of products)
    R("g_s", _p_sub("g", "s"), "SU(3) coupling g_s")
    R("g", _p_indexed(_mi("g"), ("up",)), "SU(2) coupling g; with two indices the metric g^{μν}", "g")
    R("gp", _p_const(_F("<mrow><mi>g</mi>" + _PRIME.ml + "</mrow>", "g′", atom=True)), "U(1) coupling g′")
    R("e", _p_const(_mi("e")), "unit of electric charge (broken phase); also the index letter e")
    R("y_u", _p_sub("y", "u"), "up-type Yukawa coupling")
    R("y_d", _p_sub("y", "d"), "down-type Yukawa coupling")
    R("y_e", _p_sub("y", "e"), "electron Yukawa coupling")
    R("y_k", _p_sub("y", "k"), "a Yukawa coupling y_k (k = u, d, e)")
    R("lambda", _p_const(_mi("λ")), "quartic coupling λ")
    R("m2", _p_const(_F("<msup><mi>m</mi><mn>2</mn></msup>", "m²")), "mass parameter m² of the potential")
    R("xi", _p_const(_mi("ξ")), "gauge parameter ξ (Feynman gauge ξ = 1)")
    R("v", _p_const(_mi("v")), "vacuum expectation value v")
    R("c_W", _p_sub("c", "W"), "cosine of the weak mixing angle")
    R("s_W", _p_sub("s", "W"), "sine of the weak mixing angle")
    R("M_W", _p_sub("M", "W"), "W mass")
    R("M_Z", _p_sub("M", "Z"), "Z mass")
    R("M_h", _p_sub("M", "h"), "Higgs mass")
    R("m_u", _p_sub("m", "u"), "up-quark mass")
    R("m_d", _p_sub("m", "d"), "down-quark mass")
    R("m_e", _p_sub("m", "e"), "electron mass")
    R("m_f", _p_sub("m", "f"), "fermion mass")
    R("Y", _p_named("Y"), "hypercharge Y_f", "Y(Q)")
    R("Nc", lambda e: _subsup(_mi("N"), _index(_args(e)[0])[0] if _args(e) else _mi("f"), _mi("c")),
      "colour multiplicity N^c_f", "Nc(Q)")
    R("Sinv", _p_inv("S"), "group invariant S_f^(r)", "Sinv(Q, G)")
    R("Cinv", _p_inv("C"), "Casimir C_f", "Cinv(Q)")
    R("Cad", _p_inv("Cad"), "adjoint Casimir C_ad^(r)", "Cad(G)")
    R("Ytil", _p_const(_accent(_mi("Y"), "tilde")), "Ỹ = 3 y_u² + 3 y_d² + y_e²")
    R("Ycal", lambda e: _sub(_mi("𝒴"), _index(_args(e)[0])[0]) if _args(e) else _mi("𝒴"),
      "Yukawa part 𝒴_f of the fermion constants", "Ycal(Q)")
    R("ep", _p_const(_mi("ϵ")), "regulator ϵ (lunate glyph), d = 4 − 2ϵ")
    R("eta", _p_const(_mi("η")), "evanescent tag η")
    R("d", _p_const(_mi("d")), "dimension d = 4 − 2ϵ")
    # counterterm names and operators
    R("dZ", _p_dZ, "renormalisation constant δZ_X or δZ^(r)_X", "dZ(A, G)")
    R("dZc", lambda e: _F("<mrow><mi>δ</mi>" + _sub(_mi("Z"), _frag(_args(e)[0])).ml + "</mrow>",
                          "δZ_" + _frag(_args(e)[0]).tx), "coupling constant δZ_{c}", "dZc(g_s)")
    R("dpar", _p_dpar, "parameter counterterm δx", "dpar(y_u)")
    R("rc", _p_named("r"), "coefficient r_α of the restoring counterterm", "rc(alpha)")
    R("Qop", _p_named("Q"), "basis monomial Q_α", "Qop(alpha)")
    R("Nop", _p_named("N"), "breaking basis monomial N_β", "Nop(beta)")
    R("Xc", _p_named("X"), "coordinate of the breaking X_1", "Xc(beta)")
    R("S0", _p_const(_sub(_mi("S"), _mn("0"))), "tree-level action S_0")
    R("Sren", _p_const(_sub(_mi("S"), _mi("ren", normal=True))), "renormalised action S_ren")
    R("Delta1", _p_const(_sub(_mi("Δ"), _mn("1"))), "tree-level breaking Δ_1 = 𝐬 S_0")
    R("X1", _p_const(_sub(_mi("X"), _mn("1"))), "one-loop breaking X_1")
    R("deltaS1", _p_const(_row(_mi("δ"), _sup(_mi("S"), _paren(_mn("1"))))), "singular counterterm δS^(1)")
    R("DeltaS1", _p_const(_row(_mi("Δ"), _sup(_mi("S"), _paren(_mn("1"))))), "finite counterterm ΔS^(1)")
    R("DeltaS1ST", _p_const(_row(_mi("Δ"), _subsup(_mi("S"), _mi("ST", normal=True), _paren(_mn("1"))))),
      "finite restoring counterterm ΔS^(1)_ST")
    R("Gamma", _p_const(_mi("Γ")), "effective action Γ")
    R("sop", _p_const(_mi("𝐬")), "BRST operator 𝐬")
    R("bop", _p_const(_mi("𝐛")), "linearised Slavnov-Taylor operator 𝐛")
    R("bdop", _p_const(_sub(_mi("𝐛"), _mi("d"))), "𝐛_d (d dimensions)")
    R("STop", lambda e: _row(_mi("𝒮"), _paren(_frag(_args(e)[0]))) if _args(e) else _mi("𝒮"),
      "Slavnov-Taylor operator 𝒮(F)", "STop(S0)")
    R("Kop", _p_const(_mi("K")), "pole-part operator K")
    R("Pphys", _p_const(_mi("𝒫")), "physical projector 𝒫")
    R("Evan", _p_const(_mi("ℰ")), "evanescent projector ℰ = 1 − 𝒫")
    # method notation of the report (theory cards): Taylor operator, degree, pole operators at fixed eta,
    # R-bar, order-l objects, evanescent / closed / breaking operator names with the accent on the letter
    R("Taylor", lambda e: _sup(_mi("T"), _paren(_frag(_args(e)[0]))) if _args(e) else _mi("T"),
      "Taylor (UV-expansion) operator T^(n) of degree n", "Taylor(2)")
    R("deg", lambda e: _row(_mi("deg", normal=True), _paren(_frag(_args(e)[0]))) if _args(e)
      else _mi("deg", normal=True), "superficial degree of divergence deg(G)", "deg(G)")
    R("Kf", _p_const(_sup(_mi("K"), _paren(_mi("η")))), "pole part at fixed η, η kept: K^(η)")
    R("Ke", _p_const(_sub(_mi("K"), _mi("η"))), "pole part at fixed η, then η → ϵ: K_η")
    R("Rop", _p_const(_mi("ℛ")), "R-operation ℛ")
    R("Rbarop", _p_const(_accent(_mi("ℛ"), "bar")), "R-bar operation (subtraction of all subdivergences)")
    R("Dl", lambda e: _sup(_mi("Δ", normal=True), _paren(_frag(_args(e)[0]), "[", "]")),
      "order-l breaking operator Δ^[l]", "Dl(l)")
    R("deltaS", lambda e: _row(_mi("δ"), _sup(_mi("S"), _paren(_frag(_args(e)[0])))),
      "singular counterterm of order l, δS^(l)", "deltaS(l)")
    R("DeltaS", lambda e: _row(_mi("Δ", normal=True), _sup(_mi("S"), _paren(_frag(_args(e)[0])))),
      "finite counterterm of order l, ΔS^(l)", "DeltaS(l)")
    R("Delta1r", lambda e: _subsup(_mi("Δ", normal=True), _mn("1"), _paren(_frag(_args(e)[0]))),
      "part Δ_1^(r) of the breaking with the ghosts of the factor r", "Delta1r(G)")
    R("Ssum", lambda e: _sup(_mi("S"), _paren(_index(_args(e)[0])[0])) if _args(e) else _mi("S"),
      "S^(r) = sum of S_f^(r) over all fermion multiplets", "Ssum(G)")
    R("Eop", _p_named("E"), "evanescent operator E^(r)_X", "Eop(A, G)")
    R("Etil", _p_named("E", accent="tilde"), "evanescent operator Ẽ^(r)_X (tilde on the letter)", "Etil(F, G)")
    R("Nhat", _p_named("N", accent="hat"), "evanescent ghost-number-1 operator N̂^(r)_X", "Nhat(A, G)")
    R("Oop", _p_named("𝒪"), "closed (invariant) direction 𝒪^(r)_X", "Oop(A, G)")
    R("fder", lambda e: _frac(_row(_mi("δ"), _frag(_args(e)[0])), _row(_mi("δ"), _frag(_args(e)[1]))),
      "functional derivative δF/δX", "fder(Gamma, b(a))")
    R("sqrtn", lambda e: _F(f"<msqrt>{_frag(_args(e)[0]).ml}</msqrt>", "√" + _ptx(_frag(_args(e)[0])), atom=True),
      "square root √x kept unevaluated (sqrt(2) would print as a power)", "sqrtn(2)")
    # group invariants of tensors (commutative tensors in products)
    R("gbar", _p_indexed(_mi("g"), ("up",), accent="bar"), "four-dimensional metric ḡ^{μν}", "gbar(mu, nu)")
    R("ghat", _p_indexed(_mi("g"), ("up",), accent="hat"), "(d−4)-dimensional metric ĝ^{μν}", "ghat(mu, nu)")
    R("levi", _p_indexed(_mi("ε"), ("up",)), "Levi-Civita symbol (curly ε)", "levi(mu, nu, rho, sigma)")
    R("delta", _p_indexed(_mi("δ"), ("up",)), "Kronecker delta", "delta(a, b)")
    R("f", _p_indexed(_mi("f"), ("up",)), "structure constants f^{abc}", "f(a, b, c)")
    R("fg", _p_indexed(_mi("𝔣"), ("up",)), "coupling-weighted structure constants 𝔣^{abc}", "fg(a, b, c)")
    R("dsym", _p_indexed(_mi("d"), ("up",)), "SU(3) symmetric tensor d^{abc}", "dsym(A, C, D)")
    R("sq", _p_sq(None), "square of a momentum", "sq(p(1))")
    R("sqbar", _p_sq("bar"), "square of the barred momentum", "sqbar(p(1))")
    R("sqhat", _p_sq("hat"), "square of the hatted momentum", "sqhat(p(1))")
    R("dot", _p_dot(None), "scalar product", "dot(p(1), p(2))")
    R("dotbar", _p_dot("bar"), "four-dimensional scalar product", "dotbar(p(1), p(2))")
    R("dothat", _p_dot("hat"), "(d−4)-dimensional scalar product", "dothat(p(1), p(2))")
    R("vec", _p_vec(None), "component of a momentum combination", "vec(p(1) - p(2), mu)")
    R("vecbar", _p_vec("bar"), "barred component", "vecbar(p(1), mu)")
    R("vechat", _p_vec("hat"), "hatted component", "vechat(p(1), mu)")
    R("pbar", _p_momentum_acc("bar"), "barred momentum p̄_k^μ", "pbar(1, mu)")
    R("phat", _p_momentum_acc("hat"), "hatted momentum p̂_k^μ", "phat(1, mu)")
    R("p", _p_momentum, "momentum p_k (p(k, mu): component)", "p(1, mu)")
    R("q", _p_const(_mi("q")), "momentum q")
    R("k", _p_const(_mi("k")), "momentum k; also the index letter k")
    R("pp", _p_const(_F("<mrow><mi>p</mi>" + _PRIME.ml + "</mrow>", "p′")), "outgoing momentum p′")
    R("T", _p_matrix_element(_mi("T")), "generator T^a (matrix) or (T^a)_{ij}", "T(A)")
    R("tau", _p_matrix_element(_mi("τ")), "Pauli matrix τ^I", "tau(I)")
    R("Tg", lambda e: _scripted(_mi("𝒯"), _args(e)[:1], ("up",)) if len(_args(e)) < 2 else
      _subsup(_mi("𝒯"), _index(_args(e)[1])[0], _index(_args(e)[0])[0]), "coupling-weighted generator 𝒯^a_f",
      "Tg(a, f)")
    R("Tr", _p_trace("Tr"), "trace over gauge indices", "Tr(T(A), T(C))")
    R("tr", _p_trace("tr"), "Dirac trace", "tr(gam(mu), gam(nu))")
    R("mel", _p_mel, "matrix element (x)_{ij} of a product of matrices", "mel(ord(T(I), T(J)), i, j)")
    R("comm", _p_brackets("[", "]"), "commutator", "comm(T(a), T(b))")
    R("acomm", _p_brackets("{", "}"), "anticommutator", "acomm(T(a), T(b))")
    # index helpers
    R("lo", lambda e: _index(e)[0], "force a subscript", "lo(mu)")
    R("up", lambda e: _index(e)[0], "force a superscript", "up(mu)")
    R("prime", lambda e: _index(e)[0], "primed index or symbol x′", "prime(A)")
    # Dirac structures (inside chains)
    R("gam", _p_indexed(_mi("γ"), ("up",)), "Dirac matrix γ^μ (d-dimensional; gamma is Symbolica's Γ function)",
      "gam(mu)")
    R("gammabar", _p_indexed(_mi("γ"), ("up",), accent="bar"), "four-dimensional γ̄^μ", "gammabar(mu)")
    R("gammahat", _p_indexed(_mi("γ"), ("up",), accent="hat"), "(d−4)-dimensional γ̂^μ", "gammahat(mu)")
    R("slash", _p_slash(None), "slashed momentum", "slash(p(1))")
    R("slashbar", _p_slash("bar"), "four-dimensional slashed momentum", "slashbar(p(1))")
    R("slashhat", _p_slash("hat"), "(d−4)-dimensional slashed momentum", "slashhat(p(1))")
    R("sigmaT", _p_indexed(_mi("σ"), ("up",)), "σ^{μν}", "sigmaT(mu, nu)")
    R("gamma5", _p_const(_sub(_mi("γ"), _mn("5"))), "γ_5")
    R("PL", _p_const(_sub(_mi("P"), _mi("L", normal=True))), "projector P_L = (1 − γ_5)/2")
    R("PR", _p_const(_sub(_mi("P"), _mi("R", normal=True))), "projector P_R = (1 + γ_5)/2")
    R("Pc", lambda e: _sub(_mi("P"), _sub(_mi("c"), _index(_args(e)[0])[0])) if _args(e) else
      _sub(_mi("P"), _mi("c")), "projector of the physical chirality P_{c_f}", "Pc(f)")
    R("Pcbar", lambda e: _sub(_mi("P"), _sub(_accent(_mi("c"), "bar"), _index(_args(e)[0])[0])) if _args(e) else
      _sub(_mi("P"), _accent(_mi("c"), "bar")), "the opposite projector P_{c̄_f}", "Pcbar(f)")
    # derivatives
    R("Dcov", _p_indexed(_mi("D"), ("lo",)), "covariant derivative D_μ", "Dcov(mu)")
    R("pd", _p_der(None), "partial derivative acting on one field", "pd(mu, H)")
    R("pdbar", _p_der("bar"), "four-dimensional derivative", "pdbar(mu, H)")
    R("pdhat", _p_der("hat"), "(d−4)-dimensional derivative", "pdhat(mu, H)")
    R("pd2", _p_der(None, square=True), "∂² acting on one field", "pd2(H)")
    R("pdbar2", _p_der("bar", square=True), "∂̄² acting on one field", "pdbar2(H)")
    R("pdhat2", _p_der("hat", square=True), "∂̂² acting on one field", "pdhat2(H)")
    # fields
    F = _reg_field
    F("A", _mi("A"), ("lo", "up"), doc="gauge field A^a_μ (all twelve)", example="A(mu, a)")
    F("G", _mi("G"), ("lo", "up"), doc="gluon G^A_μ", example="G(mu, A)")
    F("W", _mi("W"), ("lo", "up"), doc="W^I_μ", example="W(mu, I)")
    F("B", _mi("B"), ("lo", "up"), doc="hypercharge field B_μ", example="B(mu)")
    F("Z", _mi("Z"), ("lo", "up"), doc="Z_μ (broken phase)", example="Z(mu)")
    F("Wp", _mi("W"), ("lo",), doc="W^+_μ", example="Wp(mu)", deco=_mo("+"))
    F("Wm", _mi("W"), ("lo",), doc="W^−_μ", example="Wm(mu)", deco=_mo("−"))
    F("Agam", _mi("A"), ("lo",), doc="photon A^γ_μ", example="Agam(mu)", deco=_mi("γ"))
    F("F", _mi("F"), ("lo", "lo", "up"), doc="field strength F^a_{μν}", example="F(mu, nu, a)")
    F("Am", _mi("𝒜"), ("lo",), doc="matrix-valued gauge field 𝒜_μ", example="Am(mu)")
    F("omega", _mi("ω"), ("up",), doc="ghost ω^a", example="omega(a)")
    F("omegabar", _mi("ω"), ("up",), conj="bar", doc="antighost ω̄^a", example="omegabar(a)")
    F("wm", _mi("𝛚"), ("up",), doc="matrix-valued ghost 𝛚")
    F("b", _mi("b"), ("up",), doc="auxiliary field b^a; also the index letter b", example="b(a)")
    F("rho", _mi("ρ"), ("up", "lo"), doc="BRST source ρ^μ_a of the gauge field; also the index ρ",
      example="rho(mu, a)")
    F("rhotilde", _accent(_mi("ρ"), "tilde"), ("up", "lo"), doc="ρ̃^μ_a = ρ^μ_a + ∂^μ ω̄^a",
      example="rhotilde(mu, a)")
    F("zsrc", _mi("ζ"), ("lo",), doc="BRST source ζ_a of the ghost (zeta is Symbolica's ζ function)",
      example="zsrc(a)")
    for nm in ("Q", "U", "D", "L", "E"):
        F(nm, _mi(nm), ("lo",), doc=f"Dirac field Ψ_{nm} (physical component and sterile partner)", example=nm)
    for nm in ("Q", "U", "D", "L", "E"):
        F(nm + "bar", _mi(nm), ("lo",), conj="bar", doc=f"Dirac conjugate of {nm}", example=nm + "bar")
    F("Psi", _mi("Ψ"), ("lo",), doc="Dirac field Ψ_f", example="Psi(f)")
    F("Psibar", _mi("Ψ"), ("lo",), conj="bar", doc="Ψ̄_f", example="Psibar(f)")
    F("psi", _mi("ψ"), ("lo",), doc="physical component ψ_f", example="psi(f)")
    F("psibar", _mi("ψ"), ("lo",), conj="bar", doc="ψ̄_f", example="psibar(f)")
    F("chi", _mi("χ"), ("lo",), doc="sterile partner χ_f", example="chi(f)")
    F("chibar", _mi("χ"), ("lo",), conj="bar", doc="χ̄_f", example="chibar(f)")
    F("R", _mi("R"), ("lo",), doc="fermion source R_f", example="R(Q)")
    F("Rbar", _mi("R"), ("lo",), conj="bar", doc="fermion source R̄_f", example="Rbar(Q)")
    F("H", _mi("H"), ("lo",), doc="Higgs doublet H (H(i): component)", example="H(i)")
    F("Hdag", _mi("H"), ("lo",), conj="dag", doc="H†", example="Hdag(i)")
    F("Htilde", _accent(_mi("H"), "tilde"), ("lo",), doc="H̃ = ε H*", example="Htilde")
    F("Htildedag", _accent(_mi("H"), "tilde"), ("lo",), conj="dag", doc="H̃†", example="Htildedag")
    F("sigma", _mi("σ"), ("lo",), doc="Higgs source σ; also the index σ", example="sigma")
    F("sigmadag", _mi("σ"), ("lo",), conj="dag", doc="Higgs source σ†", example="sigmadag")
    F("Hv", _sub(_mi("H"), _mi("v")), ("lo",), doc="vacuum doublet H_v = (0, v/√2)")
    F("phi", _mi("φ"), ("lo",), doc="fluctuation φ = H − H_v")
    F("h", _mi("h"), ("lo",), doc="Higgs boson h")
    F("G0", _mi("G"), ("lo",), doc="neutral Goldstone G⁰", deco=_mn("0"))
    F("Gp", _mi("G"), ("lo",), doc="charged Goldstone G⁺", deco=_mo("+"))
    F("Gm", _mi("G"), ("lo",), doc="charged Goldstone G⁻", deco=_mo("−"))
    for nm, glyph in (("u", "u"), ("dq", "d"), ("el", "e"), ("nu", "ν")):
        F(nm + "f", _mi(glyph), ("lo",), doc=f"broken-phase Dirac field {glyph}", example=nm + "f")
        F(nm + "fbar", _mi(glyph), ("lo",), conj="bar", doc=f"Dirac conjugate of {glyph}", example=nm + "fbar")
    for nm, glyph in (("u", "u"), ("dq", "d"), ("el", "e")):
        F(nm + "fprime", _F(f"<mrow><mi>{glyph}</mi>{_PRIME.ml}</mrow>", glyph + "′"), ("lo",),
          doc=f"sterile partner {glyph}′", example=nm + "fprime")
        F(nm + "fprimebar", _F(f"<mover accent=\"true\"><mrow><mi>{glyph}</mi>{_PRIME.ml}</mrow><mo>‾</mo>"
                               f"</mover>", glyph + "′̅"), ("lo",), doc=f"Dirac conjugate of {glyph}′",
          example=nm + "fprimebar")
    R("field", _p_field_generic, "any other field, printed as its name", "field(X, mu)")
    # accents, conjugation, brackets
    R("bar", _p_accent_wrapper("bar"), "overline: four-dimensional part", "bar(G(mu, A))")
    R("hat", _p_accent_wrapper("hat"), "hat: (d−4)-dimensional part", "hat(p(1))")
    R("tilde", _p_accent_wrapper("tilde"), "tilde", "tilde(H)")
    R("dagger", _p_wrap1(lambda f, x: _sup(f if f.atom else _paren(f), _mo("†"))), "Hermitian conjugate",
      "dagger(pd(mu, H))")
    R("cc", _p_wrap1(lambda f, x: _sup(f if f.atom else _paren(f), _mo("*"))), "complex conjugate x*", "cc(H)")
    R("par", _p_wrap1(lambda f, x: _paren(f)), "explicit parentheses", "par(p(1) - p(2))")
    # ordered products, legs, integrals
    R("ord", _p_ord, "ordered (non-commutative) product", "ord(Qbar, gammabar(mu), PL, Q)")
    R("leg", _p_leg, "leg label X(p)", "leg(G(mu, A), p(1))")
    R("out", _p_const(_mtext("out")), "marks an outgoing momentum in leg(...)")
    R("intd", _p_int("d"), "∫ dᵈx", "intd(ord(H, Hdag))")
    R("int4", _p_int("4"), "∫ d⁴x", "int4(ord(H, Hdag))")
    R("sumover", _p_sumover, "sum over an index", "sumover(alpha, rc(alpha)*Qop(alpha))")
    R("eqs", _p_binop("="), "equation x = y", "eqs(x, y)")
    R("defeq", _p_binop("≔"), "definition x ≔ y", "defeq(x, y)")
    R("to", _p_binop("→"), "x → y", "to(H, Hv + phi)")
    R("withhc", _p_withhc, "x + h.c.", "withhc(ord(Qbar, H, PR, D))")
    R("hc", _p_const(_mtext("h.c.")), "the symbol h.c.")
    # sums that keep the written order of their terms (Symbolica orders the terms of a sum canonically)
    R("plus", _p_osum("+"), "x + y, terms kept in the written order", "plus(x, y)")
    R("minus", _p_osum("−"), "x − y, terms kept in the written order (a sum y is parenthesised)", "minus(x, y)")


def _mtext(s: str) -> _F:
    return _F(f"<mtext>{_esc(s)}</mtext>", s)


def _p_field_generic(e: Expression) -> _F:
    """``field(X, idx...)``: the name ``X`` with its indices (lower by default)."""
    a = _args(e)
    if not a:
        raise ConversionError("field(X, ...) needs a name")
    return _scripted(_plain_symbol(_short(a[0])) if a[0].get_type() == AtomType.Var else _frag(a[0]), a[1:],
                     ("lo",))


_PLAIN: set[str] = set()
_register_all()
# index letters and plain names: registered last (they rarely occur as factors), printed by _plain_symbol
for _n in (["mu", "nu", "alpha", "beta", "kappa", "tau", "iota", "theta", "chi_", "omega_"]
           + [c for c in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"]):
    if _n not in _SYMBOLS and not _n.endswith("_"):
        _reg(_n, (lambda nm: _p_const(_plain_symbol(nm)))(_n), "index letter or plain name", _n)
        _PLAIN.add(_n)
_INDEX_LIKE.update({"b", "e", "k", "d", "f", "g", "q", "h", "u", "v", "D", "E", "L", "Q", "U", "G", "W", "B", "A",
                    "H", "Z", "rho", "sigma", "lambda", "eta", "xi", "F", "R"})


def sym(name: str) -> Expression:
    """The site-notation symbol ``name`` (created plain if not registered)."""
    if name not in _SYMBOLS:
        _SYMBOLS[name] = S(f"{NS}::{name}")
    return _SYMBOLS[name]


_NAMED: dict[tuple, Expression] = {}


def symbol(base: str, sub: str | None = None, sup: str | None = None, *, upright_sub: bool = False,
           text: str | None = None) -> Expression:
    """A name symbol with fixed sub- and superscripts, e.g. ``symbol("Q", sub="T", sup="(G)")`` for Q^{(G)}_T or
    ``symbol("r", sub="α")``. ``base``/``sub``/``sup`` are plain strings (Greek allowed); returns a registered
    symbol usable in any expression (its source name is generated, ``nm_<hash>``)."""
    # ``sub`` may also be an expression of the site notation, printed as such (e.g. parse("y_u") for Q^+_{y_u})
    sub_expr = sub if isinstance(sub, Expression) else None
    if sub_expr is not None:
        sub = "expr:" + sub_expr.to_canonical_string()
    key = (base, sub, sup, upright_sub)
    if key in _NAMED:
        return _NAMED[key]
    import hashlib

    h = hashlib.sha1(repr(key).encode()).hexdigest()[:10]
    name = f"nm_{h}"
    b = _mi(base, normal=False) if len(base) == 1 else _F(f"<mi>{_esc(base)}</mi>", base)
    if sub_expr is not None:
        sb = _frag(sub_expr)
    else:
        sb = None if sub is None else (_mi(sub, normal=upright_sub) if len(sub) == 1 else
                                       _F(f"<mi>{_esc(sub)}</mi>" if not sub.replace(".", "").isdigit()
                                          else f"<mn>{_esc(sub)}</mn>", sub))
    sp = None if sup is None else _F(f"<mi>{_esc(sup)}</mi>" if len(sup) > 1 else f"<mi>{_esc(sup)}</mi>", sup)
    frag = _subsup(b, sb, sp)
    if text:
        frag = _F(frag.ml, text)
    s = S(f"{NS}::{name}", print=_callback)
    _PRINTERS[name] = _p_const(frag)
    _SYMBOLS[name] = s
    _NAMED[key] = s
    return s


# ================================================================================================ parsing
def parse(src: str) -> Expression:
    """Parse site-notation source (default namespace ``gdocs``)."""
    s = src.replace("**", "^")
    try:
        return E(s, default_namespace=NS)
    except Exception as err:  # Symbolica raises a plain exception with the parse position
        raise ValueError(f"cannot parse {src!r}: {err}") from None


# ================================================================================================ math objects
def _check_symbols(e: Expression) -> None:
    for s in e.get_all_symbols(True):
        name = s.get_name()
        if name.startswith(NS + "::") or name.startswith("symbolica::") or name.startswith("python::"):
            continue
        if "::" not in name:
            continue
        raise ConversionError(f"symbol {name} is not part of the site notation (convert it first)")


def _to_expression(obj) -> Expression:
    if isinstance(obj, Expression):
        return obj
    if isinstance(obj, str):
        return parse(obj)
    if isinstance(obj, (int, Fraction)):
        return coefficient(obj)
    return convert(obj)


def _math(inner: str, display: bool, cls: str = "") -> str:
    ds = "true" if display else "false"
    c = f' class="{cls}"' if cls else ""
    return f'<math xmlns="{MATHML_NS}" displaystyle="{ds}"{c}>{inner}</math>'


def M(obj, *, interactive: bool = False, display: bool = True, page_size: int | None = None,
      budget: int = VIEW_BUDGET, lead_rationals: bool = True) -> dict:
    """The math object of DATA_FORMAT.md section 2 (module docstring, section 1). ``lead_rationals``: write
    rational coefficients as leading fractions (``-4/3 g_s²`` instead of ``-(4 g_s²)/3``)."""
    e = _to_expression(obj)
    _check_symbols(e)
    if lead_rationals:
        e = display_form(e)
    ml = _typeset(e.format(mode=PrintMode.MathML, max_terms=None))
    out = {"html": _math(ml, display), "text": _text(e), "src": _source(e), "view": None}
    if interactive:
        out["view"] = _best_view(e, page_size, budget)
    return out


def _best_view(e: Expression, page_size: int | None, budget: int) -> dict | None:
    """The view with the smallest initial page among the page sizes tried; ``None`` if even that page takes more
    than 70 % of the budget (the expression is then shown statically)."""
    sizes = [page_size] if page_size else [None, 25, 12, 6]
    best = None
    for ps in sizes:
        v = e.mathml_viewer(ps) if ps else e.mathml_viewer()
        n = len(_jdump(v.page()))
        if best is None or n < best[0]:
            best = (n, ps)
    if best is None or best[0] > 0.7 * budget:
        return None
    return view_payload(e, page_size=best[1], budget=budget)


# ================================================================================================ prose
_ALLOWED_TAGS = {"p", "em", "strong", "a", "ul", "ol", "li", "code", "sub", "sup", "span", "br"}
_ALLOWED_ATTRS = {"class", "title", "data-report-page", "data-term", "data-class", "data-theory"}
_TAG = re.compile(r"<(/?)([a-zA-Z][a-zA-Z0-9]*)((?:\s+[a-zA-Z_:][-a-zA-Z0-9_:.]*(?:\s*=\s*(?:\"[^\"]*\"|'[^']*'"
                  r"|[^\s\"'=<>`]+))?)*)\s*(/?)>")
_ATTR = re.compile(r"([a-zA-Z_:][-a-zA-Z0-9_:.]*)(?:\s*=\s*(\"[^\"]*\"|'[^']*'|[^\s\"'=<>`]+))?")
_ENTITY = re.compile(r"&(?:[a-zA-Z][a-zA-Z0-9]*|#[0-9]+|#x[0-9a-fA-F]+);")


def _clean_html(text: str) -> str:
    """Keep the allowed tags and attributes, escape everything else."""
    out = []
    pos = 0
    for m in _TAG.finditer(text):
        out.append(_escape_text(text[pos:m.start()]))
        closing, tag, attrs, selfclose = m.group(1), m.group(2).lower(), m.group(3), m.group(4)
        if tag in _ALLOWED_TAGS:
            if closing:
                out.append(f"</{tag}>")
            else:
                kept = []
                for am in _ATTR.finditer(attrs or ""):
                    an = am.group(1).lower()
                    if an in _ALLOWED_ATTRS:
                        v = am.group(2) or '""'
                        if v[0] in "\"'":
                            v = v[1:-1]
                        kept.append(f' {an}="{_html.escape(_html.unescape(v), quote=True)}"')
                out.append(f"<{tag}{''.join(kept)}{' /' if selfclose and tag == 'br' else ''}>")
        else:
            out.append(_html.escape(m.group(0), quote=False))
        pos = m.end()
    out.append(_escape_text(text[pos:]))
    return "".join(out)


def _escape_text(s: str) -> str:
    res = []
    pos = 0
    for m in _ENTITY.finditer(s):
        res.append(_html.escape(s[pos:m.start()], quote=False))
        res.append(m.group(0))
        pos = m.end()
    res.append(_html.escape(s[pos:], quote=False))
    return "".join(res)


def _math_spans(text: str) -> list[tuple[str, str]]:
    """Split into ``("text", s)`` and ``("math", src)`` pieces; ``\\[[`` is a literal ``[[``."""
    pieces: list[tuple[str, str]] = []
    buf = []
    i = 0
    n = len(text)
    while i < n:
        if text.startswith("\\[[", i):
            buf.append("[[")
            i += 3
            continue
        if text.startswith("[[", i):
            depth = 0
            j = i + 2
            while j < n:
                if text.startswith("]]", j) and depth == 0:
                    break
                c = text[j]
                if c == "[":
                    depth += 1
                elif c == "]":
                    depth -= 1
                    if depth < 0:
                        raise ValueError(f"unbalanced ']' in math span starting at {text[i:i + 40]!r}")
                j += 1
            else:
                raise ValueError(f"unterminated math span {text[i:i + 40]!r}")
            pieces.append(("text", "".join(buf)))
            buf = []
            pieces.append(("math", text[i + 2:j].strip()))
            i = j + 2
            continue
        buf.append(text[i])
        i += 1
    pieces.append(("text", "".join(buf)))
    return [p for p in pieces if p[1] or p[0] == "math"]


def prose(text: str) -> str:
    """Rich text with ``[[ ... ]]`` spans converted to inline MathML (module docstring, section 1)."""
    out = []
    for kind, s in _math_spans(text):
        if kind == "text":
            out.append(_clean_html(s))
        else:
            if not s:
                raise ValueError("empty math span [[ ]]")
            try:
                e = display_form(parse(s))
                _check_symbols(e)
                ml = _typeset(e.format(mode=PrintMode.MathML))
            except ValueError as err:
                raise ValueError(f"math span [[ {s} ]]: {err}") from None
            if "<mtable data-symbolica" in ml and "flex" not in ml and ml.count("<mtr") == ml.count("<mtable"):
                # inline math in running text: drop Symbolica's single-cell layout tables (they add padding);
                # expressions laid out on several rows keep their tables
                ml = _STRIP_TABLE_OPEN.sub("", ml)
                ml = _STRIP_TABLE_CLOSE.sub("", ml)
            out.append(_math(_STRIP_ATTR.sub("", ml) if "flex" not in ml else ml, False, "mv-inline"))
    return "".join(out)


# ================================================================================================ display form
_SCALAR_VARS = {"g_s", "g", "gp", "e", "y_u", "y_d", "y_e", "y_k", "lambda", "m2", "xi", "v", "c_W", "s_W", "M_W",
                "M_Z", "M_h", "m_u", "m_d", "m_e", "m_f", "Ytil", "ep", "eta", "d"}
_SCALAR_FNS = {"rat", "rati", "cf", "Y", "Nc", "Sinv", "Cinv", "Cad", "Ycal", "dZ", "dZc", "dpar", "rc", "Xc"}
_LOOP_VARS = {"lf", "lfe"}
_NUM = re.compile(r"(-?)(\d*)(𝑖?)(?:/(\d+))?")


def _is_scalar(f: Expression) -> bool:
    t = f.get_type()
    if t == AtomType.Num:
        return True
    if t == AtomType.Var:
        name = f.get_name()
        if not name.startswith(NS + "::"):
            return True
        short = name[len(NS) + 2:]
        return short in _SCALAR_VARS or short in _PLAIN or (short not in _PRINTERS)
    if t == AtomType.Pow:
        base, ex = list(f)
        return ex.get_type() == AtomType.Num and _is_scalar(base)
    if t in (AtomType.Add, AtomType.Mul):
        return all(_is_scalar(x) for x in f)
    if t == AtomType.Fn:
        name = f.get_name()
        if not name.startswith(NS + "::"):
            return all(_is_scalar(x) for x in f)
        return name[len(NS) + 2:] in _SCALAR_FNS
    return False


def _number_parts(num: Expression) -> tuple[int, int, int, bool] | None:
    """``(sign, p, q, imaginary)`` of a real or purely imaginary rational, ``None`` otherwise."""
    m = _NUM.fullmatch(num.format_plain().replace(" ", ""))
    if not m:
        return None
    sign = -1 if m.group(1) else 1
    p = int(m.group(2)) if m.group(2) else 1
    q = int(m.group(4)) if m.group(4) else 1
    return sign, p, q, bool(m.group(3))


def _prod(xs: Iterable[Expression]) -> Expression:
    out = Expression.num(1)
    for x in xs:
        out = out * x
    return out


def _display_term(term: Expression) -> Expression:
    factors = list(term) if term.get_type() == AtomType.Mul else [term]
    num = None
    scal, loop, other = [], [], []
    for f in factors:
        t = f.get_type()
        if t == AtomType.Num:
            num = f
            continue
        if t == AtomType.Add:
            f = display_form(f)
        if t == AtomType.Var and f.get_name().startswith(NS + "::") and _short(f) in _LOOP_VARS:
            loop.append(f)
        elif _is_scalar(f):
            scal.append(f)
        else:
            other.append(f)
    rebuilt = _prod(([num] if num is not None else []) + loop + scal + other)
    if num is None:
        return rebuilt
    parts = _number_parts(num)
    if parts is None:
        return rebuilt
    sign, p, q, imag = parts
    if q == 1:
        return rebuilt
    long_sums = [o for o in other if o.get_type() == AtomType.Add and
                 (len(list(o)) > _SHORT_SUM or not _fits_fragment(o))]
    if long_sums:
        # a long sum stays a sum for the pager: the rational goes into its terms, the scalars stay in front
        target = long_sums[0]
        inner = display_form((num * target).expand())
        rest = [o for o in other if o is not target]
        return _prod(loop) * _prod(scal) * inner * _prod(rest)
    r = _SYMBOLS["rati" if imag else "rat"](p, q)
    core = _SYMBOLS["cf"](_prod(scal), r) if scal else r
    others = []
    for o in other:
        # a non-scalar sum would be printed before the coefficient (sums sort before functions): parenthesise it
        others.append(_SYMBOLS["par"](o) if o.get_type() == AtomType.Add else o)
    return sign * _prod(loop) * core * _prod(others)


def display_form(e: Expression) -> Expression:
    """``e`` with rational coefficients written as leading fractions (``rat``/``rati``/``cf``); equal in value
    (``rat(p, q) = p/q``)."""
    if e.get_type() == AtomType.Add:
        return _prod([]) * sum((_display_term(t) for t in e), Expression.num(0))
    return _display_term(e)


# ================================================================================================ converters
_MODEL = None


def _default_model():
    global _MODEL
    if _MODEL is None:
        from hepkit_renorm import builtin_model

        _MODEL = builtin_model("sm_one_family")
    return _MODEL


def _tn():
    from hepkit_renorm import tensors

    return tensors


_PARAMS = {"g_s": "g_s", "g": "g", "gp": "gp", "y_u": "y_u", "y_d": "y_d", "y_e": "y_e", "lambda": "lambda",
           "m2": "m2", "xi": "xi"}


def _scalar_to_site(expr: Expression) -> Expression:
    """Replace the package's parameter symbols, ``d``, ``eta`` by the site symbols; raise on anything else."""
    from symbolica import Replacement

    reps = []
    for sy in expr.get_all_symbols(False):
        name = sy.get_name()
        short = name.split("::")[-1]
        if name.startswith(NS + "::") or short in ("𝜋", "𝑖") or name.startswith("symbolica::"):
            continue
        if name.startswith("hepkit_renorm_param::") and short in _PARAMS:
            reps.append(Replacement(sy, _SYMBOLS[_PARAMS[short]]))
        elif name == "hepkit_renorm::d":
            reps.append(Replacement(sy, _SYMBOLS["d"]))
        elif name == "hepkit_renorm::eta":
            reps.append(Replacement(sy, _SYMBOLS["eta"]))
        elif name.startswith("python::") and short in _PARAMS:
            reps.append(Replacement(sy, _SYMBOLS[_PARAMS[short]]))
        else:
            raise ConversionError(f"coefficient contains {name}, which has no site symbol")
    for sy in expr.get_all_symbols(True):
        if sy.get_name().startswith("hepkit_renorm::") and sy.get_name() not in ("hepkit_renorm::d",
                                                                              "hepkit_renorm::eta"):
            raise ConversionError(f"coefficient contains the tensor head {sy.get_name()}")
    return expr.replace_multiple(reps) if reps else expr


def coefficient(value, *, loop: str | None = None, pole: bool = False, sign: int = 1) -> Expression:
    """A scalar coefficient in the site notation (module docstring, section 3)."""
    if isinstance(value, Expression):
        e = value
    elif isinstance(value, bool):
        raise ConversionError("a coefficient cannot be a bool")
    elif isinstance(value, int):
        e = Expression.num(value)
    elif isinstance(value, Fraction):
        e = Expression.num(value.numerator) / value.denominator
    elif isinstance(value, str):
        e = parse(value.replace("gp'", "gp"))
    else:
        raise ConversionError(f"cannot read a coefficient from {type(value).__name__}")
    e = _scalar_to_site(e)
    if sign != 1:
        e = sign * e
    pi = E("𝜋")
    if loop == "auto":
        # an explicit 1/π² (order store) becomes the loop factor
        scaled = (e * 16 * pi ** 2).expand()
        if not scaled.contains(pi):
            e = scaled
            loop = "pole" if pole else "finite"
        else:
            loop = None
    if loop not in (None, "pole", "finite"):
        raise ConversionError(f"loop must be None, 'pole', 'finite' or 'auto', not {loop!r}")
    e = e.expand()
    if loop == "pole":
        return _SYMBOLS["lfe"] * e
    if loop == "finite":
        return _SYMBOLS["lf"] * e
    return e


# ------------------------------------------------------------------------------------------------ term analysis
_GAMMAS = {"gamma": None, "gammabar": "bar", "gammahat": "hat"}
_METRICS = {"g": None, "gbar": "bar", "ghat": "hat"}


class _Obj:
    """One object of a product term: a field occurrence or leg (``kind`` "field"/"leg") or a tensor head
    (``kind`` "head", ``name`` the package head name)."""

    __slots__ = ("kind", "name", "labels", "accent", "derivs", "pos", "momentum", "info")

    def __init__(self, kind, name, labels, pos=None, derivs=None, momentum=None, info=None):
        self.kind = kind
        self.name = name
        self.labels = list(labels)
        self.accent: dict[int, str] = {}
        self.derivs = [[d, None] for d in (derivs or [])]
        self.pos = pos
        self.momentum = momentum
        self.info = info

    def all_labels(self):
        return [(k, lab) for k, lab in enumerate(self.labels) if lab is not None] + \
               [(("d", k), d[0]) for k, d in enumerate(self.derivs)]


def _key(lab: Expression) -> str:
    return lab.to_canonical_string()


def _space_of_label(lab: Expression) -> str:
    return _tn().label_space(lab)


_LABEL_HEADS = {"lor", "spin", "adj", "fund", "flav"}


def _is_label(e: Expression) -> bool:
    return e.get_type() == AtomType.Fn and e.get_name().startswith("hepkit_renorm::") and _short(e) in _LABEL_HEADS


def _head_labels(name: str, args: list[Expression]) -> list:
    """Label positions of a package head (non-labels, e.g. the momentum number of ``p``, become ``None``)."""
    return [a if _is_label(a) else None for a in args]


def _split_term(term: Expression) -> tuple[Expression, list[_Obj]]:
    """Scalar part and tensor heads of a product term of a package coefficient."""
    factors = list(term) if term.get_type() == AtomType.Mul else [term]
    scalar = Expression.num(1)
    heads: list[_Obj] = []
    for f in factors:
        reps = 1
        base = f
        if f.get_type() == AtomType.Pow:
            b, ex = list(f)
            if b.get_type() == AtomType.Fn and b.get_name().startswith("hepkit_renorm::"):
                if ex.get_type() != AtomType.Num or not ex.format_plain().isdigit():
                    raise ConversionError(f"non-integer power of a tensor head: {f}")
                reps, base = int(ex.format_plain()), b
        if base.get_type() == AtomType.Fn and base.get_name().startswith("hepkit_renorm::"):
            name = _short(base)
            if name in ("occ", "dv"):
                raise ConversionError("unexpected field encoding inside a coefficient")
            args = list(base)
            for _ in range(reps):
                heads.append(_Obj("head", name, _head_labels(name, args), info=args))
        elif base.get_type() == AtomType.Add and any(
                s.get_name().startswith("hepkit_renorm::") and _short(s) not in ("d", "eta")
                for s in base.get_all_symbols(True)):
            raise ConversionError("a coefficient term contains an unexpanded sum of tensors")
        else:
            scalar = scalar * f
    return scalar, heads


class _Term:
    """Display analysis of one product term (module docstring, section 3)."""

    def __init__(self, table, scalar: Expression, heads: list[_Obj], fields: list[_Obj], legs: list[_Obj] | None,
                 letters: dict | None = None):
        self.table = table
        self.scalar = scalar
        self.heads = heads
        self.fields = fields          # field occurrences in Grassmann order (operators)
        self.legs = legs or []        # legs (vertex functions)
        self.letters = dict(letters or {})   # label key -> display symbol (free labels of the legs)
        self.dirac_chains: list[tuple] = []   # (row obj, [heads], col obj)
        self.dirac_traces: list[list] = []
        self.group_chains: dict[str, list] = {}
        self.group_traces: dict[str, list] = {}
        self.explicit: set[str] = set()

    # -------------------------------------------------------------------------------------------- labels
    def objects(self) -> list[_Obj]:
        return self.heads + self.fields + self.legs

    def occurrences(self) -> dict[str, list[tuple[_Obj, object]]]:
        occ: dict[str, list] = {}
        for o in self.objects():
            for pos, lab in o.all_labels():
                occ.setdefault(_key(lab), []).append((o, pos))
        return occ

    @staticmethod
    def _set_label(o: _Obj, pos, new: Expression) -> None:
        if isinstance(pos, tuple):
            o.derivs[pos[1]][0] = new
        else:
            o.labels[pos] = new

    @staticmethod
    def _add_accent(o: _Obj, pos, acc: str | None) -> None:
        if acc is None:
            return
        if isinstance(pos, tuple):
            old = o.derivs[pos[1]][1]
            if old and old != acc:
                raise ConversionError("a derivative is both barred and hatted (vanishes)")
            o.derivs[pos[1]][1] = acc
            return
        old = o.accent.get(pos)
        if old and old != acc:
            raise ConversionError(f"{o.name} has a barred and a hatted projection of one index (vanishes)")
        o.accent[pos] = acc

    def absorb_metrics(self) -> None:
        """Absorb metrics and Kronecker deltas (except spinor deltas) into the objects they contract with."""
        changed = True
        while changed:
            changed = False
            for m in list(self.heads):
                if m.name in _METRICS or (m.name == "delta" and _space_of_label(m.labels[0]) != "spinor"):
                    acc = _METRICS.get(m.name) if m.name != "delta" else None
                    if m.name == "g" and m.accent:
                        acc = m.accent.get(0) or m.accent.get(1)
                    l1, l2 = m.labels
                    if _key(l1) == _key(l2):
                        raise ConversionError("unevaluated trace of a metric")
                    occ = self.occurrences()
                    # legs are not displayed: a label shared with a leg is free
                    o1 = [x for x in occ.get(_key(l1), []) if x[0] is not m and x[0].kind != "leg"]
                    o2 = [x for x in occ.get(_key(l2), []) if x[0] is not m and x[0].kind != "leg"]
                    if len(o1) > 1 or len(o2) > 1:
                        raise ConversionError("an index occurs more than twice in a term")
                    if not o1 and not o2:
                        continue
                    self.heads.remove(m)
                    if o1 and o2:
                        (ob2, p2), = o2
                        self._set_label(ob2, p2, l1)
                        self._add_accent(ob2, p2, acc)
                        self._add_accent(o1[0][0], o1[0][1], acc)
                    elif o2:
                        (ob2, p2), = o2
                        self._set_label(ob2, p2, l1)
                        self._add_accent(ob2, p2, acc)
                    else:
                        (ob1, p1), = o1
                        self._set_label(ob1, p1, l2)
                        self._add_accent(ob1, p1, acc)
                    changed = True
                    break
        # accents on remaining metrics, gammas, eps, momenta
        for h in list(self.heads):
            accs = set(h.accent.values())
            if not accs:
                continue
            if len(accs) > 1:
                raise ConversionError(f"{h.name} carries barred and hatted indices (vanishes)")
            acc = accs.pop()
            if h.name in _METRICS:
                cur = _METRICS[h.name]
                if cur and cur != acc:
                    raise ConversionError("product of a barred and a hatted metric (vanishes)")
                h.name = {"bar": "gbar", "hat": "ghat"}[acc]
                h.accent = {}
            elif h.name in _GAMMAS:
                cur = _GAMMAS[h.name]
                if cur and cur != acc:
                    raise ConversionError("barred gamma with a hatted index (vanishes)")
                h.name = {"bar": "gammabar", "hat": "gammahat"}[acc]
                h.accent = {}
            elif h.name == "eps":
                if acc == "hat":
                    raise ConversionError("Levi-Civita symbol with a hatted index (vanishes)")
                h.accent = {}

    # -------------------------------------------------------------------------------------------- chains
    def _spinor_type(self, o: _Obj) -> str | None:
        info = self.table[o.name]
        kind = info.kind
        if kind == "fermion_bar":
            return "row"
        if kind == "fermion":
            return "col"
        if kind == "source" and info.source_of:
            sk = self.table[info.source_of].kind
            if sk == "fermion":
                return "row"
            if sk == "fermion_bar":
                return "col"
        return None

    def build_dirac(self) -> None:
        spin_heads = [h for h in self.heads if h.name in ("gamma", "gammabar", "gammahat", "gamma5", "PL", "PR")
                      or (h.name == "delta" and _space_of_label(h.labels[0]) == "spinor")]
        out: dict[str, _Obj] = {}
        for h in spin_heads:
            r = h.labels[-2]
            if _key(r) in out:
                raise ConversionError("spinor index used twice as a row index")
            out[_key(r)] = h
        ends_row, ends_col = {}, {}
        for o in self.fields + self.legs:
            st = self._spinor_type(o)
            if st is None:
                continue
            info = self.table[o.name]
            k = next((i for i, sl in enumerate(info.slots) if sl.space.kind == "spinor"), None)
            if k is None:
                continue
            (ends_row if st == "row" else ends_col)[_key(o.labels[k])] = o
        for rk, row in ends_row.items():
            seq = []
            cur = rk
            steps = 0
            while cur not in ends_col:
                h = out.pop(cur, None)
                if h is None:
                    raise ConversionError(f"open spinor line starting at {row.name}")
                seq.append(h)
                cur = _key(h.labels[-1])
                steps += 1
                if steps > 64:
                    raise ConversionError("spinor line does not close")
            col = ends_col.pop(cur)
            self.dirac_chains.append((row, [h for h in seq if h.name != "delta"], col))
        if ends_col:
            raise ConversionError("spinor line without a row field")
        while out:
            start_key, h = next(iter(out.items()))
            seq = []
            cur = start_key
            while cur in out:
                h = out.pop(cur)
                seq.append(h)
                cur = _key(h.labels[-1])
            if cur != start_key:
                raise ConversionError("open spinor line without fields")
            self.dirac_traces.append([h for h in seq if h.name != "delta"])
        for h in spin_heads:
            self.heads.remove(h)

    def find_slashes(self) -> None:
        """A gamma contracted with exactly one momentum ``p(k, mu)`` becomes a slashed momentum."""
        occ = self.occurrences()
        moms = {}
        for h in self.heads:
            if h.name == "p":
                moms.setdefault(_key(h.labels[1]), []).append(h)
        for chain in [c[1] for c in self.dirac_chains] + self.dirac_traces:
            for i, g in enumerate(chain):
                if g.name not in _GAMMAS:
                    continue
                lab = g.labels[0]
                k = _key(lab)
                others = [x for x in occ.get(k, []) if x[0] is not g]
                if len(others) == 1 and others[0][0].name == "p" and others[0][0].kind == "head":
                    pm = others[0][0]
                    acc = _GAMMAS[g.name] or pm.accent.get(1)
                    slash = _Obj("head", "slash", [None], info=[pm.info[0]])
                    slash.accent = {0: acc} if acc else {}
                    chain[i] = slash
                    self.heads.remove(pm)

    def _group_spaces(self) -> list[str]:
        spaces = set()
        for o in self.objects():
            for _, lab in o.all_labels():
                sp = _space_of_label(lab)
                if sp.startswith("fund:"):
                    spaces.add(sp)
        return sorted(spaces)

    def _fund_slot(self, o: _Obj, space: str):
        info = self.table[o.name]
        for k, sl in enumerate(info.slots):
            if sl.space.kind == "fund" and f"fund:{sl.space.factor}" == space:
                return k, ("row" if sl.rep == "Fbar" else "col")
        return None

    def build_group(self, space: str) -> bool:
        """Chains of the fundamental space ``space``; ``False`` if it does not decompose (explicit indices)."""
        heads = [h for h in self.heads if h.name in ("T", "delta", "epsf")
                 and any(_space_of_label(l) == space for l in h.labels if l is not None)]
        if any(h.name == "epsf" for h in heads):
            return False
        out = {}
        for h in heads:
            i, j = (h.labels[1], h.labels[2]) if h.name == "T" else (h.labels[0], h.labels[1])
            if _key(i) in out:
                return False
            out[_key(i)] = (h, j)
        for h in heads:
            if h.name == "delta":
                # a Kronecker delta is undirected: it can be entered from either label
                i, j = h.labels[0], h.labels[1]
                if _key(j) not in out:
                    out[_key(j)] = (h, i)
        rows, cols = {}, {}
        for o in self.fields + self.legs:
            fs = self._fund_slot(o, space)
            if fs is None:
                continue
            k, typ = fs
            (rows if typ == "row" else cols)[_key(o.labels[k])] = o
        chains, used = [], set()
        for rk, row in rows.items():
            seq, cur = [], rk
            while cur not in cols:
                nxt = out.get(cur)
                if nxt is None:
                    # a Kronecker delta between two column or two row fields: not a matrix product
                    return False
                h, j = nxt
                if id(h) in used:
                    return False
                used.add(id(h))
                seq.append(h)
                cur = _key(j)
            col = cols.pop(cur)
            chains.append((row, [h for h in seq if h.name == "T"], col))
        if cols:
            return False
        traces = []
        rest = {k: v for k, v in out.items() if id(v[0]) not in used and v[0].name == "T"}
        while rest:
            k0 = next(iter(rest))
            seq, cur = [], k0
            while cur in rest:
                h, j = rest.pop(cur)
                seq.append(h)
                cur = _key(j)
            if cur != k0:
                return False
            if any(h.name != "T" for h in seq):
                return False
            traces.append(seq)
        self.group_chains[space] = chains
        self.group_traces[space] = traces
        return True


# ------------------------------------------------------------------------------------------------ letters
_LOR = ["mu", "nu", "rho", "sigma", "alpha", "beta", "kappa", "tau"]
_LOR_SRC = ["mu", "nu", "alpha", "beta", "kappa", "tau", "iota", "theta"]
_ADJ = {"G": ["A", "C", "D", "E", "F", "K", "M", "N"], "W": ["I", "J", "K", "L", "M", "N", "O", "P"]}
_FUND = {"W": ["i", "j", "k", "l", "m", "n", "r", "s"], "G": ["alpha", "beta", "kappa", "iota", "theta", "tau"]}
_FIELD_LETTERS = {"Q", "U", "D", "L", "E", "H"}


class _Letters:
    """Letters for the labels of one expression: free labels fixed by the legs, dummies per term."""

    def __init__(self, has_source: bool, avoid: Iterable[str] = ()):
        self.pools = {"lorentz": list(_LOR_SRC if has_source else _LOR)}
        for f, pool in _ADJ.items():
            self.pools[f"adjoint:{f}"] = [x for x in pool if x not in avoid]
        for f, pool in _FUND.items():
            self.pools[f"fund:{f}"] = list(pool)
        self.pools["flavour"] = ["m", "n", "r", "s"]
        self.used: dict[str, set] = {}
        self.assigned: dict[str, Expression] = {}

    def take(self, lab: Expression) -> Expression:
        k = _key(lab)
        if k in self.assigned:
            return self.assigned[k]
        sp = _space_of_label(lab)
        pool = self.pools.get(sp)
        if pool is None:
            raise ConversionError(f"no letters for the index space {sp}")
        used = self.used.setdefault(sp, set())
        taken = {str(v) for v in self.assigned.values()}
        for name in pool + [p + "_" + str(n) for n in range(1, 9) for p in pool[:1]]:
            if name in used:
                continue
            if sp.startswith("fund:G") and name in {_short(v) for v in self.assigned.values()}:
                continue
            used.add(name)
            s = sym(name) if "_" not in name else symbol(_GREEK.get(name.split("_")[0], name.split("_")[0]),
                                                         sub=name.split("_")[1])
            self.assigned[k] = s
            return s
        raise ConversionError(f"ran out of letters for {sp}")

    def copy(self) -> _Letters:
        c = _Letters.__new__(_Letters)
        c.pools = self.pools
        c.used = {k: set(v) for k, v in self.used.items()}
        c.assigned = dict(self.assigned)
        return c


# ------------------------------------------------------------------------------------------------ field symbols
_FIELD_HEADS = {"G": "G", "W": "W", "B": "B"}


def field(name: str, idx: Sequence[Expression] = (), table=None, accent: str | None = None) -> Expression:
    """Site expression of the package field ``name`` with the given (already placed) index expressions in slot
    order (Lorentz, spinor (omitted), group, flavour). Examples: ``field("G", [lo(mu), A])`` -> G^A_μ,
    ``field("gh_B")`` -> ω^B, ``field("src_Q")`` -> R̄_Q, ``field("src_H")`` -> σ†."""
    table = table if table is not None else _default_model().fields
    info = table[name]
    kind = info.kind
    idx = list(idx)
    hook = getattr(table, "site_field", None)   # field tables of other bases (gdocs.broken_fields) name their own
    if hook is not None:
        e = hook(name, idx)
        if e is not None:
            return _SYMBOLS[accent](e) if accent else e
    if kind == "gauge":
        head = _SYMBOLS[name] if name in _FIELD_HEADS else None
        if head is None:
            raise ConversionError(f"no site symbol for the gauge field {name}")
        e = head(*idx)
    elif kind in ("ghost", "antighost", "auxiliary"):
        fac = info.factor
        head = {"ghost": "omega", "antighost": "omegabar", "auxiliary": "b"}[kind]
        if not idx:
            idx = [sym(fac)]
        e = _SYMBOLS[head](*idx)
    elif kind in ("fermion", "fermion_bar"):
        mult = info.multiplet
        nm = mult if kind == "fermion" else mult + "bar"
        e = _SYMBOLS[nm](*idx) if idx else _SYMBOLS[nm]
    elif kind in ("scalar", "scalar_conj"):
        nm = "H" if kind == "scalar" else "Hdag"
        if info.multiplet not in (None, "H") and name not in ("H", "Hdag"):
            raise ConversionError(f"no site symbol for the scalar {name}")
        e = _SYMBOLS[nm](*idx) if idx else _SYMBOLS[nm]
    elif kind == "source":
        of = table[info.source_of]
        if of.kind == "gauge":
            fac = of.factor
            if len(idx) == 1:
                idx = idx + [sym(fac)]
            e = _SYMBOLS["rho"](*idx)
        elif of.kind == "ghost":
            if not idx:
                idx = [sym(of.factor)]
            e = _SYMBOLS["zsrc"](*idx)
        elif of.kind == "fermion":
            e = _SYMBOLS["Rbar"](sym(of.multiplet), *idx)
        elif of.kind == "fermion_bar":
            e = _SYMBOLS["R"](sym(of.multiplet), *idx)
        elif of.kind == "scalar":
            e = _SYMBOLS["sigmadag"](*idx) if idx else _SYMBOLS["sigmadag"]
        elif of.kind == "scalar_conj":
            e = _SYMBOLS["sigma"](*idx) if idx else _SYMBOLS["sigma"]
        else:
            raise ConversionError(f"no site symbol for the source {name}")
    else:
        raise ConversionError(f"no site symbol for the field {name} ({kind})")
    if accent:
        e = _SYMBOLS[accent](e)
    return e


_PRIORITY = {"antighost": 0, "source_boson": 1, "ghost": 2, "fermion": 3, "scalar": 4, "gauge": 5, "auxiliary": 6,
             "source": 7}


def _priority(table, name: str) -> int:
    info = table[name]
    if info.kind == "source":
        of = table[info.source_of].kind
        if of in ("gauge", "ghost"):
            return _PRIORITY["source_boson"]
        if of in ("fermion", "fermion_bar"):
            return _PRIORITY["fermion"]
        return _PRIORITY["scalar"]
    if info.kind in ("fermion", "fermion_bar"):
        return _PRIORITY["fermion"]
    if info.kind in ("scalar", "scalar_conj"):
        return _PRIORITY["scalar"]
    return _PRIORITY.get(info.kind, 8)


def _perm_sign(seq: Sequence[int]) -> int:
    inv = sum(1 for i, j in itertools.combinations(range(len(seq)), 2) if seq[i] > seq[j])
    return -1 if inv % 2 else 1


class _Renderer:
    """Turns an analysed :class:`_Term` into a site expression."""

    def __init__(self, term: _Term, letters: _Letters, vertex: bool):
        self.t = term
        self.L = letters
        self.is_vertex = vertex
        self.place: dict[str, list] = {}   # label key -> list of default positions in display order

    # -------------------------------------------------------------------------------------------- indices
    def idx(self, lab: Expression, default: str) -> Expression:
        s = self.L.take(lab)
        sp = _space_of_label(lab)
        if sp == "lorentz":
            k = _key(lab)
            seen = self.place.setdefault(k, [])
            pos = default
            if seen and seen[0] == default:
                pos = "up" if default == "lo" else "lo"
            seen.append(pos)
            return _SYMBOLS[pos](s) if pos != default else s
        if sp.startswith("fund:") or sp == "flavour":
            return _SYMBOLS["lo"](s) if default != "lo" else s
        return s

    # -------------------------------------------------------------------------------------------- heads
    def head(self, h: _Obj) -> Expression:
        n = h.name
        S_ = _SYMBOLS
        if n in _GAMMAS:
            return S_["gam" if n == "gamma" else n](self.idx(h.labels[0], "up"))
        if n == "slash":
            k = h.info[0]
            acc = h.accent.get(0)
            return S_["slash" + (acc or "")](S_["p"](k))
        if n in ("gamma5", "PL", "PR"):
            return S_[n]
        if n in _METRICS:
            return S_[n](self.idx(h.labels[0], "up"), self.idx(h.labels[1], "up"))
        if n == "eps":
            return S_["levi"](*[self.idx(l, "up") for l in h.labels])
        if n == "p":
            k = h.info[0]
            acc = h.accent.get(1)
            nm = {"bar": "pbar", "hat": "phat"}.get(acc, "p")
            return S_[nm](k, self.idx(h.labels[1], "up"))
        if n in ("sp", "spbar", "sphat"):
            k, l = h.info
            suffix = n[2:]
            if k.format_plain() == l.format_plain():
                return S_["sq" + suffix](S_["p"](k))
            return S_["dot" + suffix](S_["p"](k), S_["p"](l))
        if n == "T":
            a = self.idx(h.labels[0], "up")
            sp = _space_of_label(h.labels[1])
            if sp in self.t.explicit or not self._in_chain(h, sp):
                return S_["T"](a, self.idx(h.labels[1], "lo"), self.idx(h.labels[2], "lo"))
            return S_["T"](a)
        if n == "f":
            sp = _space_of_label(h.labels[0])
            args = [self.idx(l, "up") for l in h.labels]
            return S_["levi"](*args) if sp == "adjoint:W" else S_["f"](*args)
        if n == "dsym":
            return S_["dsym"](*[self.idx(l, "up") for l in h.labels])
        if n == "delta":
            sp = _space_of_label(h.labels[0])
            if sp.startswith("fund:") or sp == "flavour":
                return S_["delta"](*[S_["lo"](self.L.take(l)) for l in h.labels])
            return S_["delta"](*[self.idx(l, "up") for l in h.labels])
        if n == "epsf":
            return S_["levi"](*[S_["lo"](self.L.take(l)) for l in h.labels])
        raise ConversionError(f"no site notation for the head {n}")

    def _chained_ids(self) -> set[int]:
        """Heads written inside chains or traces of a suppressed space (T and Kronecker deltas)."""
        ids = set()
        for space, chains in self.t.group_chains.items():
            if space in self.t.explicit:
                continue
            for h in self.t.heads:
                if h.name in ("T", "delta") and any(_space_of_label(l) == space for l in h.labels if l is not None):
                    ids.add(id(h))
        return ids

    def _in_chain(self, h: _Obj, space: str) -> bool:
        for row, ts, col in self.t.group_chains.get(space, []):
            if any(x is h for x in ts):
                return True
        for tr in self.t.group_traces.get(space, []):
            if any(x is h for x in tr):
                return True
        return False

    # -------------------------------------------------------------------------------------------- fields
    def field_expr(self, o: _Obj) -> Expression:
        info = self.t.table[o.name]
        idx = []
        lor_accent = None
        for k, sl in enumerate(info.slots):
            lab = o.labels[k]
            sp = sl.space.kind
            if sp == "spinor":
                continue
            if sp == "lorentz":
                default = "up" if info.kind == "source" else "lo"
                idx.append(self.idx(lab, default))
                lor_accent = o.accent.get(k)
            elif sp == "adjoint":
                idx.append(self.L.take(lab))
            elif sp == "fund":
                space = f"fund:{sl.space.factor}"
                if space in self.t.explicit:
                    idx.append(self.L.take(lab))
            elif sp == "flavour":
                idx.append(self.L.take(lab))
        e = field(o.name, idx, self.t.table, accent=lor_accent)
        # derivatives (innermost first in the occurrence's order)
        ders = list(o.derivs)
        done = set()
        for i in range(len(ders)):
            if i in done:
                continue
            lab, acc = ders[i]
            pair = next((j for j in range(i + 1, len(ders)) if j not in done and _key(ders[j][0]) == _key(lab)),
                        None)
            if pair is not None:
                done.update((i, pair))
                acc2 = acc or ders[pair][1]
                e = _SYMBOLS["pd" + (acc2 or "") + "2"](e)
                continue
            done.add(i)
            e = _SYMBOLS["pd" + (acc or "")](self.idx(lab, "lo"), e)
        return e

    # -------------------------------------------------------------------------------------------- assembly
    def _blocks(self) -> list[list]:
        """Display blocks of an operator: lists of field objects and chain heads (``("T", h)`` for group
        matrices)."""
        t = self.t
        blocks: list[dict] = []
        by_field: dict[int, dict] = {}
        for row, heads, col in t.dirac_chains:
            b = {"kind": "dirac", "row": row, "col": col, "heads": list(heads), "ts": [], "attach": []}
            blocks.append(b)
            by_field[id(row)] = b
            by_field[id(col)] = b
        placed = {id(o) for b in blocks for o in (b["row"], b["col"])}
        for space, chains in t.group_chains.items():
            if space in t.explicit:
                continue
            for row, ts, col in chains:
                br, bc = by_field.get(id(row)), by_field.get(id(col))
                if br is not None and br is bc:
                    br["ts"].extend(ts)
                elif br is not None and br["row"] is row and id(col) not in placed and bc is None:
                    br["attach"].append([row, ts, col, "after_row"])
                    placed.add(id(col))
                    by_field[id(col)] = br
                elif bc is not None and bc["col"] is col and id(row) not in placed and br is None:
                    bc["attach"].append([row, ts, col, "before"])
                    placed.add(id(row))
                    by_field[id(row)] = bc
                elif br is None and bc is None and id(row) not in placed and id(col) not in placed:
                    b = {"kind": "group", "row": row, "col": col, "heads": [], "ts": list(ts), "attach": []}
                    blocks.append(b)
                    by_field[id(row)] = by_field[id(col)] = b
                    placed.update((id(row), id(col)))
                elif br is not None and bc is not None and br["kind"] == "group" and br is not bc:
                    raise _Explicit(space)
                else:
                    raise _Explicit(space)
        for o in t.fields:
            if id(o) not in placed:
                blocks.append({"kind": "single", "row": o, "col": None, "heads": [], "ts": [], "attach": []})
                placed.add(id(o))
        order = {id(o): i for i, o in enumerate(t.fields)}

        def bkey(b):
            first = b["row"]
            return (_priority(t.table, first.name), t.table.order(first.name), order[id(first)])

        blocks.sort(key=bkey)
        seqs = []
        for b in blocks:
            seq = [b["row"]]
            for row, ts, col, how in b["attach"]:
                if how == "after_row":
                    seq += [("T", h) for h in ts] + [col]
                else:
                    seq += [row] + [("T", h) for h in ts]
            heads = list(b["heads"])
            # group matrices of the block go after the last gamma-type head (they commute with Dirac matrices)
            last = max((i for i, h in enumerate(heads) if h.name in _GAMMAS or h.name == "slash"), default=-1)
            seq += heads[: last + 1] + [("T", h) for h in b["ts"]] + heads[last + 1:]
            if b["col"] is not None:
                seq.append(b["col"])
            seqs.append(seq)
        return seqs

    def _check_spaces(self, seq: list) -> None:
        """Each suppressed space must read as consecutive chains in the displayed order."""
        t = self.t
        for space, chains in t.group_chains.items():
            if space in t.explicit:
                continue
            members = []
            for x in seq:
                if isinstance(x, tuple):
                    h = x[1]
                    if any(_space_of_label(l) == space for l in h.labels[1:] if l is not None):
                        members.append(h)
                elif isinstance(x, _Obj) and x.kind == "field" and t._fund_slot(x, space) is not None:
                    members.append(x)
            expect = []
            for row, ts, col in chains:
                expect.append([row] + ts + [col])
            # members must be a concatenation of the chains in some order
            i = 0
            while i < len(members):
                ch = next((c for c in expect if c[0] is members[i]), None)
                if ch is None:
                    raise _Explicit(space)
                seg = members[i:i + len(ch)]
                if len(seg) != len(ch) or any(a is not b for a, b in zip(seg, ch)):
                    raise _Explicit(space)
                expect.remove(ch)
                i += len(ch)

    def operator(self) -> tuple[Expression, Expression]:
        """``(scalar, tensor part)`` of an operator term (the Grassmann sign is in the scalar)."""
        t = self.t
        seqs = self._blocks()
        flat = [x for s in seqs for x in s]
        self._check_spaces(flat)
        fields_in_order = [x for x in flat if isinstance(x, _Obj) and x.kind == "field"]
        if len(fields_in_order) != len(t.fields):
            raise ConversionError("internal error: a field was lost in the display order")
        odd = [t.fields.index(x) for x in fields_in_order if t.table.odd(x.name)]
        sign = _perm_sign(odd)
        # tensors outside the chains first (they are printed first), then the ordered product
        tensors = Expression.num(1)
        chained = self._chained_ids()
        for h in t.heads:
            if id(h) not in chained:
                tensors = tensors * self.head(h)
        for tr in t.dirac_traces:
            tensors = tensors * _SYMBOLS["tr"](*[self.head(h) for h in tr])
        for space, traces in t.group_traces.items():
            if space in t.explicit:
                continue
            for tr in traces:
                tensors = tensors * _SYMBOLS["Tr"](*[self.head(h) for h in tr])
        items = []
        for x in flat:
            if isinstance(x, tuple):
                items.append(self.head(x[1]))
            elif x.kind == "field":
                items.append(self.field_expr(x))
            else:
                items.append(self.head(x))
        prod = _SYMBOLS["ord"](*items) if len(items) > 1 else (items[0] if items else Expression.num(1))
        return sign * _scalar_to_site(t.scalar), tensors * prod

    def vertex(self) -> Expression:
        """The tensor part of a vertex-function term (scalar included)."""
        t = self.t
        out = Expression.num(1)
        for row, heads, col in t.dirac_chains:
            ts = []
            for space, chains in t.group_chains.items():
                if space in t.explicit:
                    continue
                for r, tl, c in chains:
                    if r is row and c is col:
                        ts.extend(tl)
            last = max((i for i, h in enumerate(heads) if h.name in _GAMMAS or h.name == "slash"), default=-1)
            seq = heads[: last + 1] + ts + heads[last + 1:]
            items = [self.head(h) for h in seq]
            if len(items) > 1:
                out = out * _SYMBOLS["ord"](*items)
            elif items:
                out = out * items[0]
        for tr in t.dirac_traces:
            out = out * _SYMBOLS["tr"](*[self.head(h) for h in tr])
        for space, traces in t.group_traces.items():
            if space in t.explicit:
                continue
            for tr in traces:
                out = out * _SYMBOLS["Tr"](*[self.head(h) for h in tr])
        chained = self._chained_ids()
        # explicit spaces: a chain of two or more generators between two legs is one matrix element (T^a T^b)_{ij}
        for space, chains in t.group_chains.items():
            if space not in t.explicit:
                continue
            for row, tl, col in chains:
                if len(tl) < 2:
                    continue
                rk = t._fund_slot(row, space)[0]
                ck = t._fund_slot(col, space)[0]
                inner = _SYMBOLS["ord"](*[_SYMBOLS["T"](self.idx(h.labels[0], "up")) for h in tl])
                out = out * _SYMBOLS["mel"](inner, self.L.take(row.labels[rk]), self.L.take(col.labels[ck]))
                chained.update(id(h) for h in tl)
        for h in t.heads:
            if id(h) not in chained:
                out = out * self.head(h)
        return _scalar_to_site(t.scalar) * out


class _Explicit(Exception):
    """A fundamental space must be written with explicit indices."""

    def __init__(self, space: str):
        super().__init__(space)
        self.space = space


def _analyse(table, scalar, heads, fields, legs, vertex: bool, letters: _Letters):
    """Analyse one term and render it; retries with explicit indices for the spaces that do not chain."""
    explicit: set[str] = set()
    for _attempt in range(4):
        hs = [_Obj(h.kind, h.name, h.labels, info=h.info) for h in heads]
        fs = [_Obj(o.kind, o.name, o.labels, pos=o.pos, derivs=[d[0] for d in o.derivs], momentum=o.momentum)
              for o in fields]
        ls = [_Obj(o.kind, o.name, o.labels, momentum=o.momentum) for o in legs]
        t = _Term(table, scalar, hs, fs, ls)
        t.absorb_metrics()
        t.build_dirac()
        t.find_slashes()
        for sp in t._group_spaces():
            if sp in explicit or not t.build_group(sp):
                explicit.add(sp)
        if vertex:
            # matrix notation only between the two legs of a Dirac chain; scalar legs keep explicit indices
            for sp in list(t.group_chains):
                ok = all(any(r is row and c is col for row, _, col in t.dirac_chains)
                         for r, _, c in t.group_chains[sp])
                if not ok:
                    explicit.add(sp)
        t.explicit = set(explicit)
        r = _Renderer(t, letters.copy(), vertex)
        try:
            if vertex:
                return r.vertex(), set(explicit)
            return r.operator(), set(explicit)
        except _Explicit as ex:
            if ex.space in explicit:
                raise ConversionError(f"cannot write the space {ex.space}") from None
            explicit.add(ex.space)
    raise ConversionError("could not find a display order")


# ------------------------------------------------------------------------------------------------ public converters
def _monomials(obj) -> tuple[object, list]:
    from hepkit_renorm.functional import Density, Monomial

    if isinstance(obj, Monomial):
        raise ConversionError("a bare Monomial carries no field table: pass operator(monomial, table=...)")
    if isinstance(obj, Density):
        return obj.fields, list(obj.monomials())
    raise ConversionError(f"cannot read monomials from {type(obj).__name__}")


def _mono_parts(table, mono) -> list[tuple]:
    """``[(scalar, heads, field objects)]`` for each term of the expanded coefficient of a monomial."""
    coeff = mono.coefficient.expand()
    terms = list(coeff) if coeff.get_type() == AtomType.Add else [coeff]
    out = []
    for term in terms:
        scalar, heads = _split_term(term)
        fields = []
        for pos, o in enumerate(mono.occurrences):
            if o.field not in table:
                raise ConversionError(f"unknown field {o.field}")
            fields.append(_Obj("field", o.field, o.slots, pos=pos, derivs=list(o.derivatives)))
        out.append((scalar, heads, fields))
    return out


def _has_source(table, names: Iterable[str]) -> bool:
    return any(table[n].kind == "source" and table[table[n].source_of].kind in ("gauge", "scalar", "scalar_conj")
               for n in names)


def operator_parts(obj, *, table=None) -> list[tuple[Expression, Expression]]:
    """``[(scalar coefficient, operator)]``, one pair per term of every monomial of ``obj`` (a ``Density`` /
    ``LocalFunctional``, a ``Monomial`` with ``table``, or a list of monomials with ``table``)."""
    from hepkit_renorm.functional import Monomial

    if isinstance(obj, Monomial):
        if table is None:
            table = _default_model().fields
        monos = [obj]
    elif isinstance(obj, (list, tuple)):
        if table is None:
            table = _default_model().fields
        monos = list(obj)
    else:
        table, monos = _monomials(obj)
    out = []
    for mono in monos:
        names = [o.field for o in mono.occurrences]
        avoid = {table[n].multiplet for n in names if table[n].multiplet} | (
            {"H"} if any(table[n].kind in ("scalar", "scalar_conj") for n in names) else set())
        for scalar, heads, fields in _mono_parts(table, mono):
            letters = _Letters(_has_source(table, names), avoid)
            out.append(_analyse(table, scalar, heads, fields, [], False, letters)[0])
    return out


def operator(obj, *, table=None, strip: bool = False, loop: str | None = None) -> Expression:
    """A Lagrangian-density operator in the site notation (module docstring, section 3). ``loop="pole"`` /
    ``"finite"``: the explicit ``1/π²`` of an order-store coefficient is written as ``lfe`` / ``lf`` in front."""
    parts = operator_parts(obj, table=table)
    if loop is not None and not strip:
        total = Expression.num(0)
        for sc, op in parts:
            total = total + sc * op
        return _with_loop(total, loop)
    if not parts:
        return Expression.num(0)
    total = Expression.num(0)
    for sc, op in parts:
        total = total + sc * op
    if strip:
        if len(parts) == 1:
            return parts[0][1]
        common, _rest = _common_factor([sc for sc, _ in parts])
        total = (total / common).expand() if common is not None else total
        return _tidy(total)
    return _tidy(total)


def _tidy(e: Expression) -> Expression:
    """Collect momentum combinations and take out a common factor."""
    if e.get_type() != AtomType.Add:
        return e
    e = _collect_momenta(e)
    if e.get_type() != AtomType.Add:
        return e
    terms = list(e)
    short = len(terms) <= _SHORT_SUM and _fits_fragment(e)
    common, rest = _common_factor(terms, functions=short, rationals=short)
    if common is None:
        return e
    inner = sum(rest, Expression.num(0))
    facs = list(common) if common.get_type() == AtomType.Mul else [common]
    needs_par = any(f.get_type() == AtomType.Fn for f in facs) or any(
        f.get_type() == AtomType.Num and _number_parts(f) is not None and _number_parts(f)[2] != 1 for f in facs)
    if needs_par and inner.get_type() == AtomType.Add:
        # Symbolica prints a sum before the function factors of a product: f^{abc} (...) needs the sum wrapped.
        # Only short sums are wrapped (a wrapped sum is one unit for the pager).
        inner = _SYMBOLS["par"](inner)
    return common * inner


def _factors_of(term: Expression) -> tuple[Expression, dict[str, tuple[Expression, int]]]:
    num = Expression.num(1)
    facs: dict[str, tuple[Expression, int]] = {}
    for f in (list(term) if term.get_type() == AtomType.Mul else [term]):
        if f.get_type() == AtomType.Num:
            num = num * f
            continue
        base, ex = f, 1
        if f.get_type() == AtomType.Pow:
            b, x = list(f)
            if x.get_type() == AtomType.Num and x.format_plain().lstrip("-").isdigit():
                base, ex = b, int(x.format_plain())
        k = base.to_canonical_string()
        facs[k] = (base, facs.get(k, (base, 0))[1] + ex)
    return num, facs


def _rational(num: Expression) -> tuple[Fraction, bool] | None:
    parts = _number_parts(num)
    if parts is None:
        return None
    sign, p, q, imag = parts
    return Fraction(sign * p, q), imag


_SHORT_SUM = 12
_FRAGMENT_LIMIT = 3500   # Symbolica omits a custom-printed fragment above 4 kB ("⋯"); stay well below


def _fits_fragment(e: Expression) -> bool:
    """True if ``par(e)`` prints as one custom fragment that Symbolica keeps (below its size bound)."""
    try:
        return len(_frag(e).ml) < _FRAGMENT_LIMIT
    except ConversionError:
        return False


def _common_factor(terms: Sequence[Expression], functions: bool = True, rationals: bool = True):
    """Common factor of a list of terms: ``(common, [term / common])`` or ``(None, terms)``. ``functions``:
    also take out function factors (f^{abc}, g^{μν}, ...); ``rationals``: also non-integer rationals."""
    if len(terms) < 2:
        return None, list(terms)
    data = [_factors_of(t) for t in terms]
    keys = set(data[0][1])
    for _, f in data[1:]:
        keys &= set(f)
    if not functions:
        keys = {k for k in keys if data[0][1][k][0].get_type() != AtomType.Fn}
    common = Expression.num(1)
    for k in keys:
        pows = [f[k][1] for _, f in data]
        if all(x > 0 for x in pows):
            pw = min(pows)
        elif all(x < 0 for x in pows):
            pw = max(pows)
        else:
            continue
        common = common * data[0][1][k][0] ** pw
    nums = [_rational(n) for n, _ in data]
    if all(n is not None for n in nums):
        imag = all(n[1] for n in nums)
        if imag or not any(n[1] for n in nums):
            vals = [n[0] for n in nums]
            from math import gcd

            g = Fraction(0)
            num_g = 0
            den_l = 1
            for v in vals:
                num_g = gcd(num_g, abs(v.numerator))
                den_l = den_l * v.denominator // gcd(den_l, v.denominator)
            g = Fraction(num_g, den_l) if num_g else Fraction(1)
            if not rationals and g.denominator != 1:
                g = Fraction(1)
            if all(v < 0 for v in vals):
                g = -g
            cnum = Expression.num(g.numerator) / g.denominator
            if imag:
                cnum = cnum * Expression.I
            common = common * cnum
    if common.format_plain() == "1":
        return None, list(terms)
    rest = [(t / common).expand() for t in terms]
    return common, rest


def _momentum_number(term: Expression) -> str:
    """Sort key of a term ``c p(k)``: the momentum ``k``."""
    for f in (list(term) if term.get_type() == AtomType.Mul else [term]):
        if _is_fn(f, "p"):
            return _args(f)[0].format_plain().rjust(4)
    return term.format_plain()


def _collect_momenta(e: Expression) -> Expression:
    """Terms that differ only by the momentum of one component ``p(k, mu)`` are combined into
    ``vec(sum c_k p(k), mu)``."""
    terms = list(e)
    for _ in range(3):
        groups: dict[str, list] = {}
        for i, t in enumerate(terms):
            num, facs = _factors_of(t)
            for k, (b, ex) in facs.items():
                if ex == 1 and _is_fn(b) and _short(b) in ("p", "pbar", "phat") and len(_args(b)) == 2:
                    mom, idx = _args(b)
                    rest = Expression.num(1)
                    for k2, (b2, ex2) in facs.items():
                        if k2 != k:
                            rest = rest * b2 ** ex2
                    key = _short(b) + "|" + idx.to_canonical_string() + "|" + rest.to_canonical_string()
                    groups.setdefault(key, []).append((i, num, mom, rest, idx, _short(b)))
        used = set()
        new_terms = []
        changed = False
        for key, items in groups.items():
            items = [x for x in items if x[0] not in used]
            moms = {x[2].to_canonical_string() for x in items}
            if len(items) < 2 or len(moms) < 2:
                continue
            combo = Expression.num(0)
            for i, num, mom, rest, idx, head in items:
                combo = combo + num * _SYMBOLS["p"](mom)
                used.add(i)
            cf_, inner = _common_factor(list(combo) if combo.get_type() == AtomType.Add else [combo])
            if cf_ is not None:
                combo = sum(inner, Expression.num(0))
            else:
                cf_ = Expression.num(1)
            # the momentum with the smallest number gets a positive coefficient: (p_1 - p_2)
            first = min((x for x in (list(combo) if combo.get_type() == AtomType.Add else [combo])),
                        key=_momentum_number)
            if first.format_plain().startswith("-"):
                combo = (-combo).expand()
                cf_ = -cf_
            head = items[0][5]
            vec = {"p": "vec", "pbar": "vecbar", "phat": "vechat"}[head]
            new_terms.append(cf_ * items[0][3] * _SYMBOLS[vec](combo, items[0][4]))
            changed = True
        if not changed:
            break
        terms = [t for i, t in enumerate(terms) if i not in used] + new_terms
    return sum(terms, Expression.num(0))


# ------------------------------------------------------------------------------------------------ legs and vertices
def _leg_objects(table, legs) -> list[_Obj]:
    return [_Obj("leg", lg.field, lg.labels, momentum=lg.momentum) for lg in legs]


def _leg_letters(table, legs, avoid=()) -> _Letters:
    names = [lg.field for lg in legs]
    L = _Letters(_has_source(table, names), avoid)
    for lg in legs:
        for lab in lg.labels:
            if _space_of_label(lab) != "spinor":
                L.take(lab)
    return L


def leg(name: str, k: int, *, table=None, letters: _Letters | None = None, labels=None,
        outgoing: bool = False, explicit: Iterable[str] | None = None) -> Expression:
    """Leg label ``X(p_k)`` of the package field ``name`` (``legs`` for a whole vertex)."""
    table = table if table is not None else _default_model().fields
    from hepkit_renorm.functional import standard_legs

    if labels is None:
        labels = standard_legs(table, [name] * k)[k - 1].labels if k >= 1 else ()
    L = letters or _leg_letters(table, [type("Lg", (), {"field": name, "labels": labels, "momentum": k})()])
    info = table[name]
    matrix_leg = any(sl.space.kind == "spinor" for sl in info.slots)
    ex = set(explicit or ())
    idx = []
    for sl, lab in zip(info.slots, labels):
        if sl.space.kind == "spinor":
            continue
        if matrix_leg and sl.space.kind == "fund" and f"fund:{sl.space.factor}" not in ex:
            continue
        idx.append(L.take(lab))
    mom = _SYMBOLS["p"](k)
    args = [field(name, idx, table), mom] + ([_SYMBOLS["out"]] if outgoing else [])
    return _SYMBOLS["leg"](*args)


def legs(fields_: Sequence[str], *, table=None, outgoing: Iterable[int] = (),
         explicit: Iterable[str] | None = None) -> list[Expression]:
    """Leg labels of a vertex with the legs ``fields_`` (momenta p_1 ... p_n, the letters of :func:`vertex`).
    Fermion legs (and fermion sources) are matrix legs: they carry group indices only for the spaces listed in
    ``explicit`` (``"fund:G"``, ``"fund:W"``; ``vertex(..., info=...)`` reports them); other legs carry all."""
    table = table if table is not None else _default_model().fields
    from hepkit_renorm.functional import standard_legs

    lgs = standard_legs(table, list(fields_))
    avoid = {table[n].multiplet for n in fields_ if table[n].multiplet}
    L = _leg_letters(table, lgs, avoid)
    out = set(outgoing)
    ex = set(explicit or ())
    return [leg(lg.field, lg.momentum, table=table, letters=L, labels=lg.labels, outgoing=lg.momentum in out,
                explicit=ex) for lg in lgs]


def vertex(obj, legs_: Sequence[str] | None = None, *, model=None, simplifier=None,
           momentum_conservation: bool = False, outgoing: Iterable[int] = (), info: dict | None = None,
           loop: str | None = None) -> Expression:
    """A vertex function in the site notation (module docstring, section 3). If ``info`` (a dict) is given it
    receives ``explicit`` (the fundamental spaces written with explicit indices), ``fields`` and ``legs`` (leg
    labels consistent with the value: fermion legs carry indices only of the explicit spaces)."""
    from hepkit_renorm.functional import LocalFunctional
    from hepkit_renorm.rules import VertexRule, vertex_rule

    if isinstance(obj, VertexRule):
        rule = obj
        table = None
    elif isinstance(obj, LocalFunctional):
        if legs_ is None:
            raise ConversionError("vertex(functional, legs): the legs are needed")
        model = model or _default_model()
        simp = simplifier or model.simplifier(expand_projectors=False)
        rule = vertex_rule(obj, list(legs_), simp, outgoing=tuple(outgoing),
                           eliminate=-1 if momentum_conservation else None)
        table = obj.fields
    else:
        raise ConversionError(f"vertex() needs a VertexRule or a LocalFunctional, not {type(obj).__name__}")
    if table is None:
        table = (model or _default_model()).fields
    value = rule.value.expand()
    names = [lg.field for lg in rule.legs]
    avoid = {table[n].multiplet for n in names if table[n].multiplet}
    letters = _leg_letters(table, rule.legs, avoid)
    leg_objs = _leg_objects(table, rule.legs)
    terms = list(value) if value.get_type() == AtomType.Add else ([] if value.format_plain() == "0" else [value])
    total = Expression.num(0)
    explicit: set[str] = set()
    for term in terms:
        scalar, heads = _split_term(term)
        val, ex = _analyse(table, scalar, heads, [], leg_objs, True, letters)
        total = total + val
        explicit |= ex
    if info is not None:
        info["explicit"] = sorted(explicit)
        info["fields"] = list(names)
        info["legs"] = legs(names, table=table, outgoing=rule.outgoing, explicit=explicit)
    return _with_loop(total, loop)


def _with_loop(total: Expression, loop: str | None) -> Expression:
    """``loop``: ``None`` (as is), ``"pole"`` / ``"finite"`` (the value carries an explicit ``1/π²`` of the
    order store, which becomes ``lfe`` / ``lf`` in front)."""
    if loop is None:
        return _tidy(total)
    if loop not in ("pole", "finite"):
        raise ConversionError(f"loop must be None, 'pole' or 'finite', not {loop!r}")
    pi = E("𝜋")
    scaled = (total * 16 * pi ** 2).expand()
    if scaled.contains(pi):
        raise ConversionError("the value is not proportional to 1/π² (no one-loop factor to take out)")
    return _SYMBOLS["lfe" if loop == "pole" else "lf"] * _tidy(scaled)


def convert(obj, **options) -> Expression:
    """Dispatch: ``VertexRule`` -> :func:`vertex`; ``Density``/``LocalFunctional``/``Monomial`` -> :func:`operator`."""
    from hepkit_renorm.functional import Density, Monomial
    from hepkit_renorm.rules import VertexRule

    if isinstance(obj, VertexRule):
        return vertex(obj, **options)
    if isinstance(obj, (Density, Monomial)):
        return operator(obj, **options)
    raise ConversionError(f"no converter for {type(obj).__name__}")


# ================================================================================================ interactive payload
_SIZES = (5, 10, 25, 50, 100, 250, 500)


def _jdump(x) -> str:
    return json.dumps(x, separators=(",", ":"), ensure_ascii=False)


def answer_key(state: str, msg: dict) -> str | None:
    """The key under which the static adapter looks up the answer to ``msg`` in state ``state`` (the JavaScript
    side builds the same string, see ``assets/viewer/mathview.js``)."""
    a = msg.get("action")
    path = msg.get("value")
    lay = int(msg.get("layout") or 0)
    ctx = msg.get("context") or ""
    if a == "preview":
        popup = bool(msg.get("popup"))
        k = [state, a, path, lay, ctx, 1 if popup else 0]
        if popup:
            k += [int(msg.get("start") or 0), int(msg.get("size") or 0)]
    elif a == "open":
        k = [state, a, path, lay, ctx, int(msg.get("start") or 0), int(msg.get("size") or 0)]
    elif a == "expand":
        k = [state, a, path, lay, ctx, int(msg.get("size") or 0)]
    elif a in ("expand_all", "collapse"):
        k = [state, a, path if "value" in msg else None, lay, ctx]
    elif a == "next":
        k = [state, a, lay]
    elif a in ("previous", "back"):
        k = [state, a]
    else:
        return None
    return _jdump(k)


def _typeset_page(p: dict) -> dict:
    """:func:`_typeset` on the MathML of a pager page and of its layout alternatives."""
    if isinstance(p.get("html"), str):
        p["html"] = _typeset(p["html"])
    if p.get("layouts"):
        p["layouts"] = [dict(alt, html=_typeset(alt["html"])) if isinstance(alt.get("html"), str) else alt
                        for alt in p["layouts"]]
    return p


class _Explorer:
    """Breadth-first exploration of the pager protocol from the initial page (module docstring, section 4)."""

    def __init__(self, expr: Expression, page_size: int | None, budget: int, max_messages: int = 400):
        self.expr = expr
        self.page_size = page_size
        self.budget = budget
        self.max_messages = max_messages
        v = self._fresh()
        p0 = _typeset_page(dict(v.page()))
        self.pages: list[dict] = [p0]
        self.page_index = {self._pkey(p0): 0}
        self.previews: list[dict] = []
        self.preview_index: dict[str, int] = {}
        self.answers: dict[str, list] = {}
        self.states = {"0": {"main": [], "page": p0}}     # state id -> messages that reach it
        self.state_of_page = {self._pkey(p0): "0"}
        self.contexts: dict[str, dict] = {}                # canonical id -> {"state", "chain"}
        self.size = len(_jdump({"v": 1, "pages": [p0], "previews": [], "answers": {}, "complete": False,
                                "bytes": 100000}))
        self.complete = True
        self.count = 0

    def _fresh(self):
        return self.expr.mathml_viewer(self.page_size) if self.page_size else self.expr.mathml_viewer()

    @staticmethod
    def _pkey(p: dict) -> str:
        return json.dumps(p, sort_keys=True, ensure_ascii=False)

    def _replay(self, state: str, chain: Sequence[dict]):
        """A viewer in ``state`` with the popups of ``chain`` open; returns (viewer, canonical -> real context).
        The history of a state is replayed message by message; popups are re-opened to obtain their (new) server
        contexts."""
        v = self._fresh()
        v.page()
        real: dict[str, str] = {}
        for m in list(self.states[state]["main"]) + list(chain):
            mm = {k: x for k, x in m.items() if k != "canonical"}
            if mm.get("context"):
                mm["context"] = real[mm["context"]]
            if mm.get("action") == "preview":
                mm["request"] = "q"
            r = v.handle_message(mm)
            if r is None or "error" in r[1]:
                raise RuntimeError("replay failed")
            if mm.get("action") == "preview":
                real[m["canonical"]] = r[1]["context"]
        return v, real

    def _add(self, key: str, value: list, extra: int) -> bool:
        cost = extra + len(key) + len(_jdump(value)) + 4
        if self.size + cost > self.budget:
            self.complete = False
            return False
        self.size += cost
        self.answers[key] = value
        return True

    def ask(self, state: str, msg: dict, chain: Sequence[dict] = ()) -> tuple | None:
        """Answer one message; returns ``("page", state)`` / ``("preview", canonical ctx, value)`` or ``None``."""
        key = answer_key(state, msg)
        if key is None or key in self.answers:
            return None
        if self.count >= self.max_messages:
            self.complete = False
            return None
        self.count += 1
        try:
            v, real = self._replay(state, chain)
        except Exception:
            return None
        mm = {k: x for k, x in msg.items() if k != "canonical"}
        if mm.get("context"):
            mm["context"] = real.get(mm["context"])
            if mm["context"] is None:
                return None
        if mm.get("action") == "preview":
            mm["request"] = "q"
        r = v.handle_message(mm)
        if r is None:
            return None
        field, value = r
        value = _typeset_page(dict(value))
        if "error" in value:
            return None
        if field == "preview":
            canon = f"c{len(self.contexts)}"
            value["context"] = canon
            value.pop("request", None)
            vk = self._pkey(value)
            extra = 0
            if vk in self.preview_index:
                idx = self.preview_index[vk]
            else:
                idx = len(self.previews)
                extra = len(_jdump(value)) + 1
                if self.size + extra > self.budget:
                    self.complete = False
                    return None
            if not self._add(key, ["v", idx, canon], extra):
                return None
            if vk not in self.preview_index:
                self.preview_index[vk] = idx
                self.previews.append(value)
            self.contexts[canon] = {"state": state, "chain": list(chain) + [dict(msg, canonical=canon)]}
            return ("preview", canon, value)
        if field == "page":
            pk = self._pkey(value)
            extra = 0
            if pk in self.page_index:
                pidx = self.page_index[pk]
            else:
                pidx = len(self.pages)
                extra = len(_jdump(value)) + 1
            new_state = self.state_of_page.get(pk)
            if new_state is None:
                new_state = str(len(self.states))
            if not self._add(key, ["p", pidx, new_state], extra):
                return None
            if pk not in self.page_index:
                self.page_index[pk] = pidx
                self.pages.append(value)
            if new_state not in self.states:
                # the history keeps the popups that were open (a main action can refer to one of them)
                self.states[new_state] = {"main": list(self.states[state]["main"]) + list(chain) + [dict(msg)],
                                          "page": value}
                self.state_of_page[pk] = new_state
            return ("page", new_state, value)
        return None

    # -------------------------------------------------------------------------------------------- candidates
    @staticmethod
    def _layouts(p: dict) -> list[dict]:
        return [p] + [{**p, **alt} for alt in (p.get("layouts") or [])]

    def _marker_msgs(self, page: dict, ctx: str | None, size: int):
        """(priority, message) of the requests a reader can trigger from the markers of ``page``."""
        out = []
        for n, lay in enumerate(self._layouts(page)):
            layout = int(lay.get("layout") or 0)
            alt = 1 if n else 0   # markers of the compact alternatives come after those of the full layout
            for i, (path, _label) in enumerate(lay.get("markers") or []):
                start = (lay.get("marker_starts") or [0] * (i + 1))[i]
                base = {"value": list(path), "layout": layout, "context": ctx}
                if ctx is None:
                    out.append((alt, dict(base, action="preview", popup=False, start=start, size=size)))
                out.append((alt, dict(base, action="preview", popup=True, start=start, size=size)))
                out.append((1 + alt, dict(base, action="open", start=start, size=size)))
        return out

    def run(self) -> dict:
        frontier: list[tuple[int, str, dict, list]] = []
        p0 = self.pages[0]
        size0 = int(p0.get("page_size") or 25)
        for pr, m in self._marker_msgs(p0, None, size0):
            frontier.append((pr, "0", m, []))
        seen_states = {"0"}
        depth = 0
        while frontier and depth < 6:
            frontier.sort(key=lambda x: x[0])
            nxt: list[tuple[int, str, dict, list]] = []
            for pr, state, msg, chain in frontier:
                res = self.ask(state, msg, chain)
                if res is None:
                    continue
                page_size = int(self.states[state]["page"].get("page_size") or 25)
                if res[0] == "preview":
                    _, canon, value = res
                    if not msg.get("popup"):
                        continue
                    newchain = list(chain) + [dict(msg, canonical=canon)]
                    # paging inside the popup
                    if value.get("next") is not None:
                        nxt.append((2, state, dict(msg, start=value["next"]), chain))
                    # nested markers of the popup
                    for pr2, m2 in self._marker_msgs(value, canon, msg.get("size", page_size)):
                        if m2["action"] == "preview" and m2["value"] == msg["value"]:
                            m2 = dict(m2, context=msg.get("context"))
                            nxt.append((2, state, m2, chain))
                            continue
                        nxt.append((3 + pr2, state, m2, newchain))
                    if self.states[state]["page"].get("supports_expansion"):
                        base = {"value": msg["value"], "layout": msg["layout"], "context": msg.get("context")}
                        if value.get("can_expand"):
                            nxt.append((4, state, dict(base, action="expand", size=msg.get("size", page_size)),
                                        chain))
                        if value.get("can_show_all"):
                            nxt.append((4, state, dict(base, action="expand_all"), chain))
                        if value.get("can_collapse"):
                            nxt.append((5, state, dict(base, action="collapse"), chain))
                else:
                    _, new_state, page = res
                    if new_state in seen_states:
                        continue
                    seen_states.add(new_state)
                    ps = int(page.get("page_size") or 25)
                    for pr2, m2 in self._marker_msgs(page, None, ps):
                        nxt.append((2 + pr2, new_state, m2, []))
                    if page.get("back"):
                        nxt.append((1, new_state, {"action": "back", "layout": int(page.get("layout") or 0)}, []))
                        if page.get("next") is not None:
                            nxt.append((2, new_state, {"action": "next", "layout": int(page.get("layout") or 0)},
                                        []))
                    if page.get("previous"):
                        nxt.append((2, new_state, {"action": "previous", "layout": 0}, []))
                    if page.get("can_collapse"):
                        nxt.append((5, new_state, {"action": "collapse", "layout": 0}, []))
            frontier = nxt
            depth += 1
        if frontier:
            self.complete = False
        return {"v": 1, "pages": self.pages, "previews": self.previews, "answers": self.answers,
                "complete": self.complete}


def view_payload(expr: Expression, *, page_size: int | None = None, budget: int = VIEW_BUDGET) -> dict:
    """The ``view`` of an interactive math object (module docstring, section 4): ``pages[0]`` is the initial
    page; ``answers[key]`` is ``["p", page index, new state]`` or ``["v", preview index, context]``."""
    ex = _Explorer(expr, page_size, budget)
    out = ex.run()
    out["bytes"] = len(_jdump(out))
    return out


# ================================================================================================ assets
def pager_javascript() -> str:
    """The pager ES module exactly as Symbolica ships it (``MathMLViewer.javascript()``)."""
    return MathMLViewer.javascript()


def write_assets(site_dir: str | Path) -> list[Path]:
    """Write ``assets/viewer/symbolica_pager.js`` and ``assets/viewer/STIXTwoMath-Regular.woff2`` below
    ``site_dir``; returns the written paths."""
    import symbolica

    d = Path(site_dir) / "assets" / "viewer"
    d.mkdir(parents=True, exist_ok=True)
    js = d / "symbolica_pager.js"
    js.write_text(pager_javascript(), encoding="utf-8")
    font_src = Path(symbolica.__file__).parent / "_data" / "fonts" / "render" / "fonts" / "STIXTwoMath-Regular.woff2"
    font = d / "STIXTwoMath-Regular.woff2"
    if not font.exists() or font.read_bytes() != font_src.read_bytes():
        shutil.copyfile(font_src, font)
    return [js, font]
