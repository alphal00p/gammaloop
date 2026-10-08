"""Canonical energy and momentum notation in expressions and both explorers."""

import json
import re
import sys
import unicodedata
import xml.etree.ElementTree as ET
from html.parser import HTMLParser
from pathlib import Path

from symbolica import E, S
from symbolica.community import hepkit as hep
from symbolica.community import tensor


class MathText(HTMLParser):
    def __init__(self, source):
        super().__init__()
        self.depth = 0
        self.text = ""
        self.feed(source)

    def handle_starttag(self, tag, attrs):
        if tag == "math":
            self.depth += 1

    def handle_endtag(self, tag):
        if tag == "math":
            self.depth -= 1

    def handle_data(self, text):
        if self.depth:
            self.text += text


def check_math(source):
    text = unicodedata.normalize("NFKC", MathText(source).text)
    assert "E" in text and "os" in text, text
    assert "OSE" not in text and "ω" not in text, text
    assert "<merror" not in source


def check_momentum(source, label, index):
    math = ET.fromstring(re.search(r"<math\b.*?</math>", source, re.DOTALL)[0])
    assert label in unicodedata.normalize("NFKC", "".join(math.itertext()))
    # The object ID is a subscript, and the selected component is a superscript.
    assert any(
        node.tag in ("msub", "msubsup") and "".join(node[1].itertext()) == str(index)
        for node in math.iter()
    )
    assert any(
        node.tag in ("msup", "msubsup") and "".join(node[-1].itertext()) == "0"
        for node in math.iter()
    )
    assert not any(node.tag == "mo" and node.text == "(" for node in math.iter())


def main(directory):
    directory.mkdir(parents=True, exist_ok=True)
    ose, q, cind = S("gammalooprs::OSE", "gammalooprs::Q", "spenso::cind")
    for head, label in [
        (q, "q"),
        (S("gammalooprs::K"), "k"),
        (S("gammalooprs::P"), "p"),
    ]:
        component = head(27, cind(0))
        value = tensor.TensorExpression(component)
        assert E(component.format_plain()) == component
        assert value.format_tensor() == f"{label}₂₇^(0)"
        assert value.to_typst() == component.to_typst()
        check_momentum(value._repr_html_(), label, 27)
        assert "lr((" in (component**2).to_typst()
    for edge in (0, 3, 27):
        energy = ose(edge)
        for expression in (energy, energy**2, 1 / (2 * energy) ** 3):
            plain = expression.format_plain()
            assert E(plain) == expression and "OSE" in plain
            assert "OSE" not in expression.to_typst()
            assert r"\mathrm{os}" in expression.to_latex()
            check_math(tensor.to_html(expression))
            check_math(tensor.TensorExpression(expression)._repr_html_())

    diagram = (
        hep.Model.phi3()
        .process(["phi"], ["phi", "phi"])
        .generate_diagrams(loops=1, max_vertices=3, maximum_bridges=0, progress=None)
        .diagrams[0]
    )
    numerator = q(diagram.internal_edges[0].id, cind(0)) ** 2
    for name, result in {
        "scalar": diagram.integrate_energy(method="cff"),
        "quadratic": diagram.integrate_energy(method="cff", numerator=numerator),
        "ltd": diagram.integrate_energy(method="ltd"),
    }.items():
        source = result._repr_html_()
        kind = "ltd" if name == "ltd" else "cff"
        data = json.loads(
            re.search(
                rf'<script type="application/json" data-{kind}>(.*?)</script>', source
            ).group(1)
        )
        for energy in data["energy_html"].values():
            check_math(energy)
        if kind == "cff":
            for surface in data["surfaces"]:
                check_math(surface["expression_html"])
        else:
            for edge, momentum in data["momentum_html"].items():
                check_momentum(momentum, "q", edge)
            for label, field in [
                ("k", "loop_energy_html"),
                ("p", "external_energy_html"),
            ]:
                for index, momentum in enumerate(data[field]):
                    check_momentum(momentum, label, index)
        for energy in result.on_shell_energies.values():
            check_math(energy._repr_html_())
        check_math(result.to_expression()._repr_html_())
        (directory / f"{name}.html").write_text(source)
    print(
        "Canonical OSE and indexed momenta: identity, powers, native/tensor math and CFF/LTD labels passed"
    )


if __name__ == "__main__":
    main(Path(sys.argv[1]))
