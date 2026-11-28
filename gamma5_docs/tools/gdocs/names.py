"""Rich-text forms of the term names of the operator catalogues.

A term keeps a plain-text ``name`` (used for the document title, tooltips and search) and gets a ``name_html``
(rich text R of ``DATA_FORMAT.md``): the same words with the mathematical symbols of the name written as inline
MathML, so that no source spelling such as ``E_A`` or ``R̄_Q¹`` appears on the page.

``name_source(plain) -> str``   prose source (``[[ ... ]]`` spans in the site notation) of a plain name; deterministic,
                                so another module can rebuild the rich form of a name it reads from a data file.
``name_html(plain) -> str``     ``mathview.prose(name_source(plain))``.
``field_chip(head, label)``     site expression of a field symbol with a subscript label, e.g. ω_G, ρ_Z, ζ_γ.
"""

from __future__ import annotations

import html
import re

from gdocs import mathview as mv

_SUP_DIGIT = {"¹": "1", "²": "2"}

_TOKEN = re.compile(
    r"−O_(?P<ok>[ude])(?P<odag>†?)"
    r"|Q⁺_y(?P<qk>[ude])"
    r"|(?P<ev>E_A|Ẽ_F|E_H|Ẽ_H)"
    r"|Δ₁"
    r"|(?P<rb>R̄|R)_(?P<rf>[QUDLE])(?P<rc>[¹²]?)"
    r"|\(with P_(?P<pc>[LR])\)"
)


def _src(e) -> str:
    return mv.M(e)["src"]


def _token(m: re.Match) -> str:
    if m.group("ok"):
        o = _src(mv.symbol("O", sub=m.group("ok")))
        return f"[[-dagger({o})]]" if m.group("odag") else f"[[-{o}]]"
    if m.group("qk"):
        return f"[[{_src(mv.symbol('Q', sub=mv.parse('y_' + m.group('qk')), sup='+'))}]]"
    if m.group("ev"):
        return {"E_A": "[[Eop(A)]]", "Ẽ_F": "[[Etil(F)]]", "E_H": "[[Eop(H)]]", "Ẽ_H": "[[Etil(H)]]"}[m.group("ev")]
    if m.group("rf"):
        head = "Rbar" if m.group("rb") == "R̄" else "R"
        comp = _SUP_DIGIT.get(m.group("rc") or "")
        return f"[[{head}({m.group('rf')}{', up(' + comp + ')' if comp else ''})]]"
    if m.group("pc"):
        return f"(with [[P{m.group('pc')}]])"
    return "[[Delta1]]"


def name_source(plain: str) -> str:
    """Prose source of a plain term name: its symbols as ``[[ ... ]]`` spans, the words HTML-escaped."""
    out, pos = [], 0
    for m in _TOKEN.finditer(plain):
        out.append(html.escape(plain[pos:m.start()], quote=False))
        out.append(_token(m))
        pos = m.end()
    out.append(html.escape(plain[pos:], quote=False))
    return "".join(out)


def name_html(plain: str) -> str:
    """Rich text R of a plain term name."""
    return mv.prose(name_source(plain))


def field_chip(head: str, label) -> "mv.Expression":
    """``head(lo(label))``: a field symbol with a subscript label (ω_G, ω̄_W, b_B, ρ_Z, ζ_γ); ``label`` is a group
    or field letter or a name symbol (``mv.symbol("+")``)."""
    lab = mv.parse(label) if isinstance(label, str) else label
    return mv.parse(head)(mv.parse("lo")(lab))
