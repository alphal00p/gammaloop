"""UFO LaTeX labels survive import and render without changing symbolic identities."""

import json
from pathlib import Path
import unicodedata
import xml.etree.ElementTree as ET

import typst
from symbolica import E, S
from symbolica.community.hep import Model, UfoLoader
from symbolica.community.spenso import Representation, TensorExpression, TensorName


def visible_math(html: str) -> str:
    math = html[html.index("<math") : html.index("</math>") + len("</math>")]
    return unicodedata.normalize("NFKC", "".join(ET.fromstring(math).itertext()))


# Display labels must work even if notebook code created the symbols first.
ee, mass, alpha = S("UFO::ee", "UFO::Me", "UFO::aS")
model = Model.standard_model()
assert model.parameter("ee").texname == "e"
assert model.parameter("Me").texname == r"\text{Me}"
assert model.parameter("aS").texname == r"\alpha _s"
reloaded = Model.from_json(model.to_json())
assert [p.texname for p in reloaded.parameters] == [p.texname for p in model.parameters]

expression = TensorExpression(ee**2 * mass + alpha)
before = expression.to_expression()
source = expression.to_typst()
assert "@preview/mitex:0.2.6" in source
assert r"\text{Me}" in expression.to_latex()
assert "ee" not in expression.to_latex()
assert expression.to_latex(max_line_length=1).startswith(r"$$\begin{gathered}")
visible = visible_math(expression.to_html())
assert "ee" not in visible and "e" in visible and "Me" in visible and "α" in visible
assert "<msub" in expression.to_html()
assert "Me" in visible_math(model.parameter("Me")._repr_html_())
ET.fromstring(expression.to_svg())
ET.fromstring(typst.compile(f"$ {source} $".encode(), format="svg"))
assert expression.to_expression() == before
assert "ee" in str(expression)

vector = TensorName("parameter_label_test::p")(Representation.mink(4))
tensor = vector(E("gammalooprs::hedge(2,1)")) * ee
assert tensor.rank == 1
visible = visible_math(tensor.to_html())
assert "e" in visible and "μ" in visible and "hedge" not in visible

# Two different parameters may legitimately have the same visual label.
custom = json.loads(Model.phi3().to_json())
for parameter in custom["parameters"]:
    parameter["texname"] = r"\alpha_s"
custom_model = Model.from_json(json.dumps(custom))
g, m = S("UFO::g", "UFO::mass")
same_labels = TensorExpression(g + m)
assert same_labels.to_typst().count("mi(") == 2
assert visible_math(same_labels.to_html()).count("α") == 2
assert same_labels.to_expression() == g + m
for parameter in custom["parameters"]:
    parameter.pop("texname")
unlabelled = Model.from_json(json.dumps(custom))
assert all(p.texname is None for p in unlabelled.parameters)
assert "mitex" not in same_labels.to_typst()

root = Path(__file__).resolve().parents[3]
imported = UfoLoader(simplify_model=False).load(root / "assets/models/ufo/sm").model
assert {p.name: p.texname for p in imported.parameters} == {
    p.name: p.texname for p in model.parameters
}
print(
    "Parameter labels: JSON/UFO metadata, MiTeX MathML/SVG/source, tensor indices and exact algebra passed"
)
