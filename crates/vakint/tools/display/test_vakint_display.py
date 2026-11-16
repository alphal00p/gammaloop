"""Test of the Python twin against the Rust golden file (run with the dashboard's
Python, which has the symbolica wheel): python -I test_vakint_display.py [--bless]

- every table form (vakint-owned and foreign) renders as stored in
  display_forms.python.golden (plain, latex, mathml, typst);
- the vakint-owned forms render exactly like Rust (src/display_forms.golden) in the
  plain and latex modes, and the spenso-dot conversion gives the same canonical string;
  MathML of a whole expression may differ in Symbolica's own wrapper markup between
  the wheel and the Rust pin, so it is compared to the Python golden only.
"""

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import vakint_display as vd  # noqa: E402
from symbolica import E, PrintMode  # noqa: E402

GOLDEN = HERE / "display_forms.python.golden"
RUST_GOLDEN = HERE.parents[1] / "src" / "display_forms.golden"
MODES = [("plain", None), ("latex", PrintMode.Latex), ("mathml", PrintMode.MathML), ("typst", PrintMode.Typst)]
SPENSO_SAMPLE = "dot(p(0),p(0))^2*coef(1/(d-4))+3*dot(p(1),k(2))+x*dot(vec(1),p(1))"


def parse(text):
    return E(text, default_namespace="vakint")


def rendered():
    lines = []
    for form in vd.load_table()["forms"]:
        x = parse(form["sample"])
        for mode, pm in MODES:
            if pm is None:
                # Rust's PrintOptions::file(): namespaced plain text
                text = x.format(mode=PrintMode.Symbolica, show_namespaces=True, max_line_length=None)
            else:
                text = x.format(mode=pm, max_line_length=None)
            lines.append(f"{form['name']}\t{mode}\t{text.replace(chr(10), ' ')}")
    return lines


DOUBLED = ["((", r"\left(\left(", "<mo>(</mo><mo>(</mo>", "<mo>(</mo><mrow><mo>(</mo>"]


def no_enclosing_parentheses():
    """User rule 2026-11-16: no form carries its own parentheses; no doubled parentheses
    around a single symbol as a factor, inside log(...) or as a power base (plain, LaTeX, MathML)."""
    from symbolica import AtomType

    for form in vd.load_table()["forms"]:
        for key in ("latex", "mathml", "typst"):
            assert not form[key].startswith(("(", r"\left(", "<mo>(</mo>", "<mrow><mo>(</mo>")), (form["name"], form[key])
        sample = parse(form["sample"])
        single = sample
        if sample.get_type() in (AtomType.Mul, AtomType.Add, AtomType.Pow):
            for part in [sample] + list(sample) + [q for p in sample for q in (p if p.get_type() == AtomType.Pow else [])]:
                if part.get_type() in (AtomType.Var, AtomType.Fn) and part.get_name() == form["name"]:
                    single = part
                    break
        contexts = {"alone": single, "factor": single * parse("x"), "log": parse("log(x)").replace(parse("x"), single), "power": single**2}
        for label, x in contexts.items():
            for _, pm in MODES:
                text = x.format(mode=pm) if pm is not None else x.format()
                assert not any(d in text for d in DOUBLED), (form["name"], label, text)
                if label in ("alone", "log") and pm is not None and single.get_type() == AtomType.Var:
                    inner = text.split("log", 1)[-1]
                    assert "(" + "(" not in inner, (form["name"], label, text)


def main():
    # The vectors carry spenso tags in vakint; declare them as Rust does before parsing.
    report = vd.register()
    assert not report["kept_default"], report["kept_default"]
    from symbolica.community.tensor import TensorName

    for head in ("vakint::p", "vakint::k"):
        TensorName.vector(head)
    no_enclosing_parentheses()
    actual = rendered()
    if "--bless" in sys.argv:
        GOLDEN.write_text("\n".join(actual) + "\n")
        return
    expected = GOLDEN.read_text().splitlines()
    assert expected == actual, "\n".join(f"- {e}\n+ {a}" for e, a in zip(expected, actual) if e != a)
    rust = {tuple(l.split("\t")[:2]): l for l in RUST_GOLDEN.read_text().splitlines()}
    # Samples that are one atom: the two Symbolica versions may order factors differently.
    from symbolica import AtomType

    single = {f["name"] for f in vd.load_table()["forms"] if parse(f["sample"]).get_type() in (AtomType.Var, AtomType.Fn)}
    compared = 0
    for line in actual:
        key = tuple(line.split("\t")[:2])
        if key in rust and key[1] == "latex" and key[0] in single:
            assert rust[key] == line, f"Rust/Python differ:\n  rust   {rust[key]}\n  python {line}"
            compared += 1
    spenso = vd.to_spenso_dots(parse(SPENSO_SAMPLE)).to_canonical_string()
    rust_spenso = next(l for l in RUST_GOLDEN.read_text().splitlines() if l.startswith("to_spenso_dots\t"))
    assert parse(rust_spenso.split("\t", 1)[1]) == parse(spenso), (rust_spenso, spenso)
    print(f"OK: {len(actual)} renderings, {compared} identical to Rust, spenso conversion identical")


if __name__ == "__main__":
    main()
