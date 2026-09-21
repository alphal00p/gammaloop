"""Exercise native Linnest highlights in an installed FeynKit host with typst-py."""

import importlib
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

fk = importlib.import_module(
    f"symbolica.community.{sys.argv[1] if len(sys.argv) > 1 else 'feynkit'}"
)
model = fk.Model(Path(__file__).parent / "fixtures/scalars_2p_3p.json")
diagrams = model.generate_diagrams(
    ["scalar_0"],
    ["scalar_0", "scalar_0"],
    loops=1,
    max_vertices=3,
    vertex_allow=["V_3_SCALAR_000"],
    allow_self_loops=True,
).diagrams
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
assert diagram.to_linnest(highlight=graph.empty_subgraph()) == source
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
    svg = diagram.to_svg(highlight=selected)
    root = ET.fromstring(svg)
    assert any(
        element.get("stroke", "").lower() == "#ffd166" for element in root.iter()
    ), "Linnest's default underlay is missing"
    assert "prefers-color-scheme:dark" in svg
    assert 'data-theme="dark"' in svg
    assert diagram.to_json() == snapshot
assert "<figure" in diagram.to_html(highlight=internal)

other = fk.FeynmanDiagram.from_json(model, snapshot)
for invalid in (other.to_linnet().full_subgraph(), object()):
    for render in (diagram.to_linnest, diagram.to_svg, diagram.to_html):
        try:
            render(highlight=invalid)
        except (ValueError, TypeError):
            pass
        else:
            raise AssertionError("foreign highlight selections must be rejected")

graph.reverse_edge(0)
for render in (diagram.to_linnest, diagram.to_svg, diagram.to_html):
    try:
        render(highlight=internal)
    except (ValueError, ReferenceError):
        pass
    else:
        raise AssertionError("stale highlight selections must be rejected")
assert diagram.to_json() == snapshot
assert diagram.to_linnest() == source
cross_section = model.generate_diagrams(
    ["scalar_0"],
    ["scalar_0", "scalar_0"],
    kind="cross_section",
    loops=1,
    max_vertices=2,
    vertex_allow=["V_3_SCALAR_000"],
    allow_self_loops=True,
).diagrams[0]
snapshot = cross_section.to_json()
for selected in (
    cross_section.to_linnet().full_subgraph(),
    cross_section.cuts[0].left.subgraph,
):
    root = ET.fromstring(cross_section.to_svg(highlight=selected))
    assert any(
        element.get("stroke", "").lower() == "#ffd166" for element in root.iter()
    )
assert cross_section.to_json() == snapshot
print("installed graph highlighting checks passed")
