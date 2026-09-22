"""Check interactive native SVG output for diagrams and physics subgraph views."""

import importlib
import json
import sys
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

from marimo._output.formatting import try_format

fk = importlib.import_module(
    f"symbolica.community.{sys.argv[1] if len(sys.argv) > 1 else 'feynkit'}"
)
model = fk.Model(Path(__file__).parent / "fixtures/scalars_2p_3p.json")
diagram = next(
    candidate
    for candidate in model.generate_diagrams(
        ["scalar_0"],
        ["scalar_0", "scalar_0"],
        loops=1,
        max_vertices=3,
        vertex_allow=["V_3_SCALAR_000"],
        allow_self_loops=True,
    ).diagrams
    if len(
        {
            edge.source if edge.source is not None else edge.target
            for edge in candidate.external_edges
        }
    )
    > 1
)
region = diagram.filter(edge=lambda edge: not edge.data.is_external)
snapshot = diagram.to_json()
namespace = "{http://www.w3.org/2000/svg}"
xlink = "{http://www.w3.org/1999/xlink}"
edges = {edge.id: edge for edge in diagram.edges}
expected = {("node", vertex.id) for vertex in diagram.vertices} | {
    ("edge", edge_id) for edge_id in edges
}
canonical = ET.fromstring(diagram.to_linnet()._repr_html_())
shared_script = canonical.find(namespace + "script").text

for value in (diagram, region):
    svg = value.render()
    root = ET.fromstring(svg)
    assert root.get("data-linnet-interactive") == "true"
    assert "feynkit-diagram-svg" in root.get("class", "")
    assert 'data-theme="dark"' in svg
    targets = root.findall(".//*[@data-linnet-kind]")
    identities = {
        (target.attrib["data-linnet-kind"], int(target.attrib["data-linnet-id"]))
        for target in targets
    }
    assert identities == expected
    assert Counter(
        (target.attrib["data-linnet-kind"], int(target.attrib["data-linnet-id"]))
        for target in targets
        if target.get("tabindex") == "0"
    ) == Counter(dict.fromkeys(expected, 1))
    for target in targets:
        assert target.get("href") is None
        assert target.get(xlink + "href") is None
        assert target.get("role") == "button"
        assert target.get("aria-pressed") == "false"
        if target.get("data-linnet-kind") == "edge":
            edge = edges[int(target.attrib["data-linnet-id"])]
            detail = json.loads(target.attrib["data-linnet-detail"])
            assert detail["particle"] == edge.particle_name
            assert detail["pdg"] == edge.particle_pdg
            assert (detail["source"], detail["sink"]) == (edge.source, edge.target)
            if edge.is_external:
                assert detail["external-state"] == edge.external_state
    scripts = root.findall(namespace + "script")
    assert len(scripts) == 1 and scripts[0].text == shared_script
    output = try_format(value)
    assert output.exception is None, output.traceback
    assert output.mimetype == "text/html"
    assert "<iframe" in output.data
    assert "data-linnet-interactive" in output.data
    assert diagram.to_json() == snapshot
    if value is region:
        assert any(
            element.get("stroke", "").lower() == "#77777773"
            and element.get("stroke-dasharray")
            for element in root.iter()
        ), "inspection must retain the muted, dotted subgraph context"

print("installed FeynmanDiagram and Subgraph SVG interaction checks passed")
