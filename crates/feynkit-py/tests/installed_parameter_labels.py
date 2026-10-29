"""Parameter labels survive JSON and UFO import and preserve the exact symbolic algebra."""

import json
import sys
import unicodedata
import xml.etree.ElementTree as ET
from pathlib import Path

from symbolica import E, S
from symbolica.community.hepkit import Model, ModelError, UfoLoader
from symbolica.community.tensor import Representation, TensorExpression, TensorName


def check_optional_svg(expression):
    project = expression.to_typst_project()
    assert {"main.typ", "render.typ", "notation.typ", "tree.cbor"} <= project.keys()
    assert all(isinstance(value, bytes) and value for value in project.values())
    try:
        import typst
    except ImportError:
        try:
            expression.to_svg()
        except ImportError as error:
            assert "optional Python package 'typst'" in str(error)
        else:
            raise AssertionError("SVG export succeeded without its optional compiler")
    else:
        # The optional test environment supplies the real Python Typst compiler.
        assert callable(typst.compile)
        ET.fromstring(expression.to_svg())


def math_element(html: str) -> ET.Element:
    math = html[html.index("<math") : html.index("</math>") + len("</math>")]
    return ET.fromstring(math)


def visible_math(html: str) -> str:
    return unicodedata.normalize("NFKC", "".join(math_element(html).itertext()))


# Display labels must work even if notebook code created the symbols first.
ee, mass, alpha = S("UFO::ee", "UFO::Me", "UFO::aS")
# Electron-mass notation comes from the serialized model, even when the
# original symbol already existed in the kernel.
standard = Model.standard_model()
electron_mass = standard.parameter("Me")
assert electron_mass.mathmlname is not None
expected_subscript = ET.fromstring(electron_mass.mathmlname)
assert expected_subscript.tag == "msub"
assert [node.text for node in expected_subscript] == ["m", "e"]
rendered = math_element(electron_mass._repr_html_())
subscript = rendered.find(".//{*}msub")
assert subscript is not None
assert [node.text for node in subscript] == ["m", "e"]
squared_mass = TensorExpression(mass**2)
squared = math_element(squared_mass.to_html())
power = squared.find(".//{*}msup")
assert power is not None and power[0].find(".//{*}msub") is not None
assert "".join(power[1].itertext()) == "2"
viewer = squared_mass.paged()
try:
    page = ET.fromstring(viewer._render_page(0)["html"])
    assert page.find(".//{*}msup//{*}msub") is not None
finally:
    viewer.close()
assert squared_mass.to_expression() == mass**2

# Choose labels in the fixture so rendering coverage is independent of the
# bundled model's presentation choices.
definition = json.loads(standard.to_json())
labels = {"ee": "e", "Me": r"\text{testmass}", "aS": r"\alpha_s"}
typst_labels = {"ee": "e", "Me": '"testmass"', "aS": "alpha_s"}
mathml_labels = {
    "ee": "<mi>e</mi>",
    "Me": "<mi>displaymass</mi>",
    "aS": "<msub><mi>α</mi><mi>s</mi></msub>",
}
for parameter in definition["parameters"]:
    if parameter["name"] in labels:
        parameter["texname"] = labels[parameter["name"]]
        parameter["typstname"] = typst_labels[parameter["name"]]
        parameter["mathmlname"] = mathml_labels[parameter["name"]]
model = Model.from_json(json.dumps(definition))
reloaded = Model.from_json(model.to_json())
assert [p.texname for p in reloaded.parameters] == [p.texname for p in model.parameters]

assert [p.mathmlname for p in reloaded.parameters] == [
    p.mathmlname for p in model.parameters
]
assert reloaded.parameter("Me").mathmlname == mathml_labels["Me"]

expression = TensorExpression(ee**2 * mass + alpha)
before = expression.to_expression()
source = expression.to_typst()
assert "mitex" not in source
assert "alpha_s" in source
assert labels["Me"] in expression.to_latex()
assert "ee" not in expression.to_latex()
assert expression.to_latex(max_line_length=1).startswith(r"$$\begin{gathered}")
visible = visible_math(expression.to_html())
assert (
    "ee" not in visible
    and "e" in visible
    and "displaymass" in visible
    and "α" in visible
)
assert "<msub" in expression.to_html()
assert "displaymass" in visible_math(model.parameter("Me")._repr_html_())
check_optional_svg(expression)
assert "testmass" in source
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
    parameter["typstname"] = "alpha_s"
    parameter["mathmlname"] = mathml_labels["aS"]
custom_model = Model.from_json(json.dumps(custom))
g, m = S("UFO::g", "UFO::mass")
same_labels = TensorExpression(g + m)
assert same_labels.to_typst().count("alpha_s") == 2
assert visible_math(same_labels.to_html()).count("α") == 2
assert same_labels.to_expression() == g + m
for parameter in custom["parameters"]:
    parameter.pop("texname")
    parameter.pop("typstname")
    parameter.pop("mathmlname")
unlabelled = Model.from_json(json.dumps(custom))
assert all(p.texname is None and p.mathmlname is None for p in unlabelled.parameters)
assert "mass" in visible_math(unlabelled.parameter("mass")._repr_html_())
assert "α" not in visible_math(same_labels.to_html())
custom["parameters"][0]["mathmlname"] = '<mi onclick="run()">x</mi>'
try:
    Model.from_json(json.dumps(custom))
except ModelError as error:
    assert "MathML" in str(error) or "mathmlname" in str(error)
else:
    raise AssertionError("accepted an event handler in a model label")
assert "mitex" not in same_labels.to_typst()

root = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(root / "assets/models/ufo"))
import sm

imported = UfoLoader(simplify_model=False).load(root / "assets/models/ufo/sm").model
source_parameters = {p.name: p for p in sm.all_parameters}
for parameter in imported.parameters:
    source_parameter = source_parameters[parameter.name]
    assert parameter.texname == source_parameter.texname
    assert parameter.mathmlname == getattr(source_parameter, "mathmlname", None)
print(
    "Parameter labels passed JSON/UFO round trips, native MathML, Typst/SVG, tensor indices and exact algebra checks."
)
