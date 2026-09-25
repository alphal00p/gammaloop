"""Model records expose their content in text, notebook HTML and process blobs."""

import html
import json
import xml.etree.ElementTree as ET
from pathlib import Path

import linnet
from IPython.lib.pretty import pretty
from symbolica.community import hep

model = hep.Model.standard_model()
vertex = model.vertex_rule("V_98")
text = repr(vertex)
for value in (*vertex.particles, *vertex.lorentz_structures):
    assert value in text
for row in vertex.couplings:
    for coupling in row:
        if coupling is not None:
            assert coupling in text
rich = vertex._repr_html_()
assert all(
    label in rich
    for label in ("Color structures", "Lorentz structures", "Coupling terms", "Orders")
)
assert pretty(vertex) == repr(vertex)

for members in (
    model.particles,
    model.vertex_rules,
    model.parameters,
    model.couplings,
    model.lorentz_structures,
    model.propagators,
    model.functions,
    model.form_factors,
):
    for member in members:
        rendered = member._repr_html_()
        assert html.escape(member.name) in rendered
        assert "<table" in rendered
        assert member.name in repr(member)
        assert "\x1b" not in repr(member)
        assert pretty(member) == repr(member)

assert "charge=" in repr(model.particle("e-"))
assert "mass=" in repr(model.particle("e-"))
assert "expression=" in repr(model.couplings[0])
assert "structure=" in repr(model.lorentz_structures[0])
assert "numerator=" in repr(model.propagators[0])

process = model.process(["e-", "e+"], ["a", "a"], vertex_allow=[vertex])
svg = process.render()
root = ET.fromstring(svg)
assert root.tag == "{http://www.w3.org/2000/svg}svg"
assert list(root.iter("{http://www.w3.org/2000/svg}path"))
assert "prefers-color-scheme:dark" in svg
assert "vertex_allow=[V_98]" in html.unescape(process._repr_html_())
assert ET.fromstring(process._repr_svg_()).tag == root.tag
assert ET.fromstring(process.render(config=linnet.RenderConfig())).tag == root.tag
for other in (
    model.process(["H"], ["b", "b~"]),
    model.process([], []),
    process.with_final_state_alternatives([["a", "a"], ["mu-", "mu+"]]),
):
    assert ET.fromstring(other.render()).tag == root.tag

# The same asset preparation must preserve existing diagram rendering.
diagram = process.generate_diagrams(progress=None)[0]
assert ET.fromstring(diagram.render()).tag == root.tag

# User model metadata is escaped before being inserted into rich output.
scalar = hep.Model(Path(__file__).parent / "fixtures/scalars_2p_3p.json")
definition = json.loads(scalar.to_json())
definition["vertex_rules"][0]["name"] = 'V_<tag>&"quoted"'
definition["functions"] = [{"name": "square", "arguments": ["z"], "expression": "z^2"}]
definition["form_factors"] = [
    {"name": "FF", "type": "scalar", "value": "1/(1+z^2)"}
]
escaped = hep.Model.from_json(json.dumps(definition))
for member in (escaped.functions[0], escaped.form_factors[0]):
    assert member.name in member._repr_html_()
    assert "z" in repr(member)
    assert pretty(member) == repr(member)

custom = next(v for v in escaped.vertex_rules if "<tag>" in v.name)
assert "V_&lt;tag&gt;&amp;&quot;quoted&quot;" in custom._repr_html_()
assert "<tag>" not in custom._repr_html_()
print(
    "Model displays: all member records, SVG processes, alternatives, vacuum, escaping and diagram rendering passed"
)
