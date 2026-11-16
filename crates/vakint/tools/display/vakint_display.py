"""Display forms of vakint's symbols for a Python Symbolica session.

The Python twin of ``crates/vakint/src/display.rs``: both read the same table,
``crates/vakint/src/display_forms.json``, so an expression renders alike in
Rust and Python (MathML for the expression explorer, LaTeX, Typst). Plain,
Mathematica and Sympy printing are never changed: canonical strings stay as
they are.

Usage (before parsing any vakint output in the process):

    import vakint_display
    vakint_display.register()               # once, before E(...) of vakint text
    x = vakint_display.to_spenso_dots(E(text))   # vakint::dot -> spenso::dot
    x._display_()                            # Symbolica's MathML explorer

Symbolica fixes a symbol's print function when the symbol is created, so
``register()`` must run before the first parse that creates these symbols;
``register()`` returns the names it could not decorate (already defined).
Registering up front gives these symbols earlier Symbolica ids, which can change
the factor ORDER of strings printed in this process (Symbolica orders by
creation). Use it in display processes only, never where printed strings are
compared byte by byte; vakint itself (display.rs) creates each decorated
symbol exactly where it created it before, so its output order is unchanged.
The table also carries display forms for symbols vakint emits but does not
own (gammalooprs' rstar_mass_squared, the UFO coupling gs_minus, the beta
projection tags); this module registers them too (display only).
"""

from __future__ import annotations

import json
import pathlib
from functools import reduce
from typing import Any

from symbolica import AtomType, E, Expression, PrintMode, S

TABLE_PATH = pathlib.Path(__file__).resolve().parents[2] / "src" / "display_forms.json"

# Same limits as display.rs (MATHML_MAX_ATOM_BYTES, MATHML_MAX_FRAGMENT_BYTES).
MATHML_MAX_ATOM_BYTES = 1024
MATHML_MAX_FRAGMENT_BYTES = 3500
_MODES = {PrintMode.Latex: "latex", PrintMode.MathML: "mathml", PrintMode.Typst: "typst"}
_PLUS = {"latex": "+", "mathml": "<mo>+</mo>", "typst": "+"}
_MINUS = {"latex": "-", "mathml": "<mo>−</mo>", "typst": "-"}
_COMMA = {"latex": ",", "mathml": "<mo>,</mo>", "typst": ","}
_ATTRIBUTES = {
    "Symmetric": "is_symmetric",
    "Antisymmetric": "is_antisymmetric",
    "Linear": "is_linear",
    "Scalar": "is_scalar",
    "Real": "is_real",
    "Integer": "is_integer",
    "Positive": "is_positive",
}


def load_table(path: str | pathlib.Path | None = None) -> dict[str, Any]:
    table = json.loads(pathlib.Path(path or TABLE_PATH).read_text())
    assert table["schema"] == "vakint.display_forms.v1", table["schema"]
    return table


def _paren(mode: str, inner: str) -> str:
    if mode == "latex":
        return rf"\left({inner}\right)"
    if mode == "mathml":
        return f"<mrow><mo>(</mo>{inner}<mo>)</mo></mrow>"
    return f"({inner})"


def _render(x: Expression, mode: str, pm: PrintMode) -> str:
    out = x.format(mode=pm)
    return f"<mrow>{out}</mrow>" if mode == "mathml" else out


def _negative_number(x: Expression) -> bool:
    return x.get_type() == AtomType.Num and x.to_canonical_string().startswith("-")


def _composite(x: Expression) -> bool:
    t = x.get_type()
    if t == AtomType.Add:
        return True
    if t == AtomType.Num:
        return _negative_number(x)
    if t == AtomType.Mul:
        return any(_negative_number(f) for f in x)
    return False


def _signed(x: Expression, mode: str, pm: PrintMode) -> str:
    if x.get_type() == AtomType.Num:
        negative = _negative_number(x)
        magnitude = -x if negative else x
        sign = _MINUS[mode] if negative else _PLUS[mode]
        if magnitude == E("1"):
            return sign
        return sign + _render(magnitude, mode, pm)
    inner = _render(x, mode, pm)
    if _composite(x) or x.get_type() == AtomType.Mul:
        inner = _paren(mode, inner)
    return _PLUS[mode] + inner


def _parenthesised(x: Expression, mode: str, pm: PrintMode) -> str:
    inner = _render(x, mode, pm)
    return _paren(mode, inner) if _composite(x) else inner


def _short_name(x: Expression) -> str:
    return x.get_name().rsplit("::", 1)[-1]


def _vector(x: Expression, mode: str, pm: PrintMode) -> str:
    if x.get_type() == AtomType.Fn and len(x) == 1:
        base, label = _short_name(x), _render(x[0], mode, pm)
        if mode == "latex":
            return f"{base}_{{{label}}}"
        if mode == "mathml":
            return f"<msub><mi>{base}</mi>{label}</msub>"
        return f"{base}_({label})"
    return _parenthesised(x, mode, pm)


def _join(args: list[Expression], mode: str, pm: PrintMode) -> str:
    return _COMMA[mode].join(_render(a, mode, pm) for a in args)


def _expand(template: str, args: list[Expression], mode: str, pm: PrintMode) -> str:
    out, rest = [], template
    while "⟨" in rest:
        start = rest.index("⟨")
        out.append(rest[:start])
        end = rest.index("⟩", start)
        key, rest = rest[start + 1 : end], rest[end + 1 :]
        if key == "*":
            out.append(_join(args, mode, pm))
        elif key == "init":
            out.append(_join(args[:-1], mode, pm))
        elif key == "last":
            out.append(_render(args[-1], mode, pm))
        elif key.endswith("*"):
            out.append(_join(args[int(key[:-1]) :], mode, pm))
        else:
            index, _, modifier = key.partition(":")
            a = args[int(index)]
            out.append(
                {
                    "": _render,
                    "+": _signed,
                    "p": _parenthesised,
                    "vec": _vector,
                }[modifier](a, mode, pm)
            )
    out.append(rest)
    return "".join(out)


def _accepts(arity: Any, nargs: int | None) -> bool:
    if nargs is None:
        return arity == 0
    return arity == nargs or (arity == "*" and nargs > 0)


def make_printer(form: dict[str, Any]):
    """The print callback of one table entry (None = Symbolica's default)."""

    def printer(x: Expression, mode: PrintMode, **kwargs) -> str | None:
        name = _MODES.get(mode)
        if name is None:
            return None
        # The MathML explorer shows a custom fragment as one atomic block (a hole
        # above 4096 bytes): large arguments keep the native, pageable printing.
        if name == "mathml" and x.get_byte_size() > MATHML_MAX_ATOM_BYTES:
            return None
        out = _bare(x, mode, name)
        if out is not None and name == "mathml" and len(out) > MATHML_MAX_FRAGMENT_BYTES:
            return None
        # A form never carries enclosing parentheses itself: they come only from the
        # printer state (power base, or a sign-led form inside a product).
        if out is not None and (
            (form.get("parenthesise_as_power_base") and kwargs.get("in_exp_base"))
            or (form.get("parenthesise_in_product") and kwargs.get("in_product"))
        ):
            return _paren(name, out)
        return out

    def _bare(x: Expression, mode: PrintMode, name: str) -> str | None:
        try:
            if x.get_type() == AtomType.Var:
                if not _accepts(form["arity"], None):
                    return None
                return _expand(form[name], [], name, mode)
            if x.get_type() == AtomType.Fn:
                args = list(x)
                if not _accepts(form["arity"], len(args)):
                    return None
                templates = form
                if form.get("when_args_equal") and len(args) > 1 and all(a == args[0] for a in args):
                    templates = form["when_args_equal"]
                return _expand(templates[name], args, name, mode)
        except Exception:  # noqa: BLE001  (a display callback never breaks printing)
            return None
        return None

    return printer


def register(include_foreign: bool = True, table: dict[str, Any] | None = None) -> dict[str, list]:
    """Define every table symbol with its display form (call once, before parsing)."""
    table = table or load_table()
    report: dict[str, list] = {"registered": [], "kept_default": []}
    for form in table["forms"]:
        if not form["rust"] and not include_foreign:
            continue
        kwargs = {_ATTRIBUTES[a]: True for a in form.get("attributes", [])}
        try:
            S(form["name"], print=make_printer(form), **kwargs)
            report["registered"].append(form["name"])
        except Exception as error:  # noqa: BLE001  (symbol existed already)
            report["kept_default"].append((form["name"], str(error).splitlines()[0]))
    return report


def to_spenso_dots(x: Expression, dimension: Expression | None = None) -> Expression:
    """vakint::dot(p(i), q(j)) of one-label vectors -> spenso::dot(p(i, mink(D)), q(j, mink(D))).

    The same conversion as display::to_spenso_dot_notation in Rust; spenso's
    printers then show p_i·q_j. D defaults to vakint::d (spenso hides it)."""
    from symbolica.community.tensor import Representation, TensorName, dot

    mink = Representation.mink(dimension if dimension is not None else E("vakint::d"))

    def vector(v: Expression) -> bool:
        return v.get_type() == AtomType.Fn and len(v) == 1 and "spenso::rank1" in v.get_tags()

    def walk(e: Expression) -> Expression:
        t = e.get_type()
        if t == AtomType.Fn:
            args = [walk(a) for a in e]
            if e.get_name() == "vakint::dot" and len(args) == 2 and all(vector(a) for a in args):
                left, right = (TensorName.vector(a.get_name())(a[0], mink) for a in args)
                return dot(left, right).to_expression()
            return e if all(a == b for a, b in zip(args, e)) else S(e.get_name())(*args)
        if t == AtomType.Add:
            return reduce(lambda a, b: a + b, (walk(a) for a in e))
        if t == AtomType.Mul:
            return reduce(lambda a, b: a * b, (walk(a) for a in e))
        if t == AtomType.Pow:
            base, exponent = list(e)
            return walk(base) ** walk(exponent)
        return e

    return walk(x)
