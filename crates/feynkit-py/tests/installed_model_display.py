"""Model records expose their content in text, notebook HTML and process blobs."""

import html
import json
import re
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

# A large process blob must leave four visible legs, with incoming legs left of
# outgoing ones, without inheriting the much wider amplitude centroid target.
width, height = map(float, root.attrib["viewBox"].split()[2:])
assert width / height < 1.5
carriers = [
    path
    for path in root.iter("{http://www.w3.org/2000/svg}path")
    if path.get("fill") == "none"
    and path.get("stroke-linecap") == "round"
    and not path.get("d", "").rstrip().lower().endswith("z")
]
assert len(carriers) == 4
edge_x = [
    float(re.search(r"translate\(([-\d.]+)", p.attrib["transform"])[1])
    for p in carriers
]
assert max(edge_x[:2]) < min(edge_x[2:])

# Radius changes retain visible external legs, and explicit layout settings win.
for config in (
    linnet.RenderConfig(drawing=linnet.DrawOptions(node_radius=5)),
    linnet.RenderConfig(
        drawing=linnet.DrawOptions(node_radius=linnet.AUTO, node_min_radius=5)
    ),
    linnet.RenderConfig(layouts=linnet.LayoutOptions(length_scale=0.3)),
):
    larger = ET.fromstring(process.render(config=config))
    assert float(larger.attrib["viewBox"].split()[2]) > width
    assert (
        sum(
            path.get("fill") == "none"
            and path.get("stroke-linecap") == "round"
            and not path.get("d", "").rstrip().lower().endswith("z")
            and any(
                abs(float(n)) > 1e-6
                for n in re.findall(r"-?\d+(?:\.\d+)?", path.attrib["d"])
            )
            for path in larger.iter("{http://www.w3.org/2000/svg}path")
        )
        == 4
    )

# A single flow has no preferred side: its external legs spread around the blob.
for incoming, outgoing in ((["a"] * 4, []), ([], ["a"] * 4)):
    radial = ET.fromstring(model.process(incoming, outgoing).render())
    radial_width, radial_height = map(float, radial.attrib["viewBox"].split()[2:])
    assert 0.7 < radial_width / radial_height < 1.5

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
definition["form_factors"] = [{"name": "FF", "type": "scalar", "value": "1/(1+z^2)"}]
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
