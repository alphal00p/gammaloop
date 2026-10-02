"""Native highlights use external metrics; source exports use a fake compiler."""

import xml.etree.ElementTree as ET
import sys
from types import ModuleType
from unittest.mock import patch
from pathlib import Path

from symbolica.community import graph as rendering
from symbolica.community import hepkit as fk


def fixture_linnest(owner, **kwargs):
    # Source-only control: fake compiler returns vector paths with known metrics.
    # Real typography is checked by the explicit optional typst-py harness.
    compiler = ModuleType("typst")
    calls = []

    def compile_labels(main, *, root, format):
        assert format == "svg"
        main, root = Path(main), Path(root)
        assert main.is_absolute() and root.is_absolute()
        assert main == root / "main.typ"
        staged = {
            path.relative_to(root).as_posix(): path.read_bytes()
            for path in root.rglob("*")
            if path.is_file()
        }
        assert staged.keys() == {"main.typ"}
        source = staged["main.typ"].decode()
        assert "linnest-metrics:" in source
        calls.append(source)
        vector_page = (
            '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 12">'
            '<path d="M1 1H23V11H1Z" fill="black"/>'
            '<a href="linnest-metrics:24,9,9,12"><rect width="1" height="1"/></a>'
            "</svg>"
        ).encode()
        return [vector_page] * (source.count("#pagebreak()") + 1)

    compiler.compile = compile_labels
    with patch.dict(sys.modules, {"typst": compiler}):
        result = owner.to_linnest(**kwargs)
    assert calls, "explicit outlined export did not use the fixture compiler"
    return result


def external_metrics_snapshot(view):
    # EXTERNAL METRICS FIXTURE in points, not browser measurements.
    request = view.measurement_request()
    return view.measured(
        measurements={
            "id": request["id"],
            "metrics": [
                {"id": label["id"], "width": 24.0, "height": 12.0, "baseline": 9.0}
                for label in request["labels"]
            ],
        }
    )


model = fk.Model(Path(__file__).parent / "fixtures/scalars_2p_3p.json")
diagrams = (
    model.process(
        ["scalar_0"], ["scalar_0", "scalar_0"], vertex_allow=["V_3_SCALAR_000"]
    )
    .generate_diagrams(loops=1, max_vertices=3, allow_self_loops=True)
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
graph = diagram.to_graph()
snapshot = diagram.to_json()
source = fixture_linnest(diagram)
empty = graph.empty_subgraph()
assert fixture_linnest(diagram, highlight=empty) != source
empty_svg = ET.fromstring(
    external_metrics_snapshot(diagram.render(highlight=empty)).to_svg()
)
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
    highlighted_source = fixture_linnest(diagram, highlight=selected)
    assert highlighted_source != source
    svg = external_metrics_snapshot(diagram.render(highlight=selected)).to_svg()
    root = ET.fromstring(svg)
    assert any(element.get("fill") == "none" for element in root.iter())
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
assert (
    "<figure" in external_metrics_snapshot(diagram.render(highlight=internal)).to_html()
)

other = fk.FeynmanDiagram.from_json(model, snapshot)
for invalid in (other.to_graph().full_subgraph(), object()):
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
assert fixture_linnest(diagram) == source
cross_section = (
    model.process(
        ["scalar_0"], ["scalar_0", "scalar_0"], vertex_allow=["V_3_SCALAR_000"]
    )
    .generate_cross_section(loops=1, max_vertices=2, allow_self_loops=True)
    .diagrams[0]
)
snapshot = cross_section.to_json()
for selected in (
    cross_section.to_graph().full_subgraph(),
    cross_section.cuts[0].left.subgraph,
):
    root = ET.fromstring(
        external_metrics_snapshot(cross_section.render(highlight=selected)).to_svg()
    )
    assert any(
        element.get("stroke", "").lower() == "#ffd166" for element in root.iter()
    )
assert cross_section.to_json() == snapshot
print("installed graph highlighting checks passed")

# Rendering options and routing share the interactive SVG path, without mutating
# the diagram or a caller-owned configuration.
import json


diagram = next(
    candidate
    for candidate in diagrams
    if len(list(candidate.loop_momentum_bases(limit=2))) == 2
)
snapshot = diagram.to_json()
basis = diagram.loop_momentum_basis
config = rendering.RenderSettings(
    layouts=rendering.LayoutSettings(impred_steps=40, impred_pull=1.5),
)
style = fk.DiagramStyle(
    show_particle=False,
    show_edge_index=True,
)
assert (
    external_metrics_snapshot(
        diagram.render(style=fk.DiagramStyle(show_momentum=True, momentum_arrows=True))
    ).to_svg()
    == external_metrics_snapshot(diagram.render(lmb=basis)).to_svg()
)
svg = external_metrics_snapshot(
    diagram.render(config=config, style=style, lmb=basis)
).to_svg()
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
assert style.show_particle is False and style.show_edge_index is True
assert diagram.to_json() == snapshot

alternative = next(
    candidate
    for candidate in diagram.loop_momentum_bases()
    if candidate.loop_edges != basis.loop_edges
)
assert fixture_linnest(diagram, lmb=alternative) != fixture_linnest(diagram, lmb=basis)
assert "<svg" in external_metrics_snapshot(diagram.render(lmb=alternative)).to_svg()
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
        raise AssertionError("render config must be a typed RenderSettings")

region = diagram.filter(edge=lambda edge: not edge.is_external)
assert (
    "<svg"
    in external_metrics_snapshot(
        region.render(lmb=basis, config=config, style=style)
    ).to_svg()
)
assert (
    "<figure"
    in external_metrics_snapshot(
        diagram.render(
            style=fk.DiagramStyle(show_momentum=True, momentum_arrows=True),
            config=config,
        )
    ).to_html()
)
print("installed rendering options and momentum checks passed")

# Index controls must survive the native diagram-to-renderer boundary.
unlabelled = ET.fromstring(
    external_metrics_snapshot(
        diagram.render(style=fk.DiagramStyle(show_particle=False))
    ).to_svg()
)
indexed = ET.fromstring(
    external_metrics_snapshot(
        diagram.render(style=fk.DiagramStyle(show_particle=False, show_node_index=True))
    ).to_svg()
)
assert len(indexed.findall(".//{http://www.w3.org/1998/Math/MathML}math")) > len(
    unlabelled.findall(".//{http://www.w3.org/1998/Math/MathML}math")
), "node-index labels must be visible when requested"
