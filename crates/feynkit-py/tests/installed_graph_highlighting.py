"""Exercise native Linnest highlights in an installed FeynKit host with typst-py."""

import importlib
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

fk = importlib.import_module(
    f"symbolica.community.{sys.argv[1] if len(sys.argv) > 1 else 'feynkit'}"
)
model = fk.Model(Path(__file__).parent / "fixtures/scalars_2p_3p.json")
diagrams = (
    fk.Process(model, ["scalar_0"], ["scalar_0", "scalar_0"])
    .generate_diagrams(
        loops=1, max_vertices=3, vertex_allow=["V_3_SCALAR_000"], allow_self_loops=True
    )
    .diagrams
)
diagram = next(
    candidate
    for candidate in diagrams
    if len(
        {
            edge.source if edge.source is not None else edge.target
            for edge in candidate.external_edges
        }
    )
    > 1
)
graph = diagram.to_linnet()
snapshot = diagram.to_json()
source = diagram.to_linnest()
empty = graph.empty_subgraph()
assert diagram.to_linnest(highlight=empty) != source
empty_svg = ET.fromstring(diagram.render(highlight=empty))
assert any(
    element.get("stroke", "").lower() == "#77777773" for element in empty_svg.iter()
)
assert not any(
    element.get("stroke", "").lower() == "#ffd166" for element in empty_svg.iter()
)
internal = graph.filter(edge=lambda edge: not edge.data.is_external)
paired = next(edge for edge in graph.edges() if not edge.data.is_external)
highlights = [
    internal,
    graph.full_subgraph(),
    graph.subgraph(half_edges=[paired.source.index]),
    graph.subgraph(half_edges=[paired.sink.index]),
    graph.filter(edge=lambda edge: edge.data.is_external),
]
for selected in highlights:
    highlighted_source = diagram.to_linnest(highlight=selected)
    assert highlighted_source != source
    assert "fill: none" in highlighted_source
    svg = diagram.render(highlight=selected)
    root = ET.fromstring(svg)
    assert any(
        element.get("stroke", "").lower() == "#ffd166" for element in root.iter()
    ), "selected half-edges must retain gold strokes"
    if selected != graph.full_subgraph():
        assert any(
            element.get("stroke", "").lower() == "#77777773"
            and element.get("stroke-dasharray")
            for element in root.iter()
        ), "unselected half-edges must remain visible, muted and dotted"
    assert "prefers-color-scheme:dark" in svg
    assert 'data-theme="dark"' in svg
    assert diagram.to_json() == snapshot
assert "<figure" in diagram.to_html(highlight=internal)

other = fk.FeynmanDiagram.from_json(model, snapshot)
for invalid in (other.to_linnet().full_subgraph(), object()):
    for render in (diagram.to_linnest, diagram.render, diagram.to_html):
        try:
            render(highlight=invalid)
        except (ValueError, TypeError):
            pass
        else:
            raise AssertionError("foreign highlight selections must be rejected")

graph.reverse_edge(0)
for render in (diagram.to_linnest, diagram.render, diagram.to_html):
    try:
        render(highlight=internal)
    except (ValueError, ReferenceError):
        pass
    else:
        raise AssertionError("stale highlight selections must be rejected")
assert diagram.to_json() == snapshot
assert diagram.to_linnest() == source
cross_section = (
    fk.Process(model, ["scalar_0"], ["scalar_0", "scalar_0"])
    .generate_cross_section(
        loops=1, max_vertices=2, vertex_allow=["V_3_SCALAR_000"], allow_self_loops=True
    )
    .diagrams[0]
)
snapshot = cross_section.to_json()
for selected in (
    cross_section.to_linnet().full_subgraph(),
    cross_section.cuts[0].left.subgraph,
):
    root = ET.fromstring(cross_section.render(highlight=selected))
    assert any(
        element.get("stroke", "").lower() == "#ffd166" for element in root.iter()
    )
assert cross_section.to_json() == snapshot
print("installed graph highlighting checks passed")

# Rendering options and routing share the interactive SVG path, without mutating
# the diagram or a caller-owned configuration.
import json

import linnet as ln

diagram = next(
    candidate
    for candidate in diagrams
    if len(list(candidate.loop_momentum_bases(limit=2))) == 2
)
snapshot = diagram.to_json()
basis = diagram.loop_momentum_basis
config = ln.RenderConfig(
    layouts=ln.LayoutOptions(
        seed=17,
        steps=40,
        epochs=10,
        internal_label_length_scale=0.8,
        external_label_length_scale=0.7,
        external_centroid_bias=1.5,
    ),
    template_options={"show-particle": False, "show-edge-index": True},
)
assert diagram.render(momenta=True) == diagram.render(lmb=basis)
svg = diagram.render(config=config, lmb=basis)
root = ET.fromstring(svg)
assert root.tag == "{http://www.w3.org/2000/svg}svg"
assert 'data-linnet-interactive="true"' in svg
edge_details = [
    json.loads(element.attrib["data-linnet-detail"])
    for element in root.iter()
    if element.get("data-linnet-kind") == "edge"
]
assert edge_details, "configurable SVGs must retain edge hover information"
for signature in basis.edge_signatures.values():
    assert signature.format_momentum() in json.dumps(edge_details)
assert config.template_options == {"show-particle": False, "show-edge-index": True}
assert diagram.to_json() == snapshot

alternative = next(
    candidate
    for candidate in diagram.loop_momentum_bases()
    if candidate.loop_edges != basis.loop_edges
)
assert diagram.to_linnest(lmb=alternative) != diagram.to_linnest(lmb=basis)
assert "<svg" in diagram.render(lmb=alternative)
assert diagram.to_json() == snapshot

foreign = fk.FeynmanDiagram.from_json(model, snapshot).loop_momentum_basis
for render in (diagram.render, diagram.to_html, diagram.to_linnest):
    try:
        render(lmb=foreign)
    except fk.DiagramError as error:
        assert "different diagram" in str(error)
    else:
        raise AssertionError("foreign momentum bases must be rejected")
    try:
        render(config={"steps": 10})
    except TypeError:
        pass
    else:
        raise AssertionError("render config must be a typed RenderConfig")

region = diagram.filter(edge=lambda edge: not edge.data.is_external)
assert "<svg" in region.render(lmb=basis, config=config)
assert "<figure" in diagram.to_html(momenta=True, config=config)
print("installed rendering options and momentum checks passed")

# Index controls must survive the native diagram-to-renderer boundary.
unlabelled = ET.fromstring(
    diagram.render(config=ln.RenderConfig(template_options={"show-particle": False}))
)
indexed = ET.fromstring(
    diagram.render(
        config=ln.RenderConfig(
            template_options={"show-particle": False, "show-node-index": True}
        )
    )
)
assert len(indexed.findall(".//{http://www.w3.org/2000/svg}use")) > len(
    unlabelled.findall(".//{http://www.w3.org/2000/svg}use")
), "node-index labels must be visible when requested"
