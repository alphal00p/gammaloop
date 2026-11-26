"""Native highlights use external metrics; source exports use a fake compiler."""

import base64
import hashlib
import re
import xml.etree.ElementTree as ET
import sys
from types import ModuleType
from unittest.mock import patch
from pathlib import Path

from symbolica.community import graph as rendering
from symbolica.community import hepkit as fk


def normalize_svg_fonts(svg):
    """Hash embedded WOFF2 bytes while preserving the SVG and CSS around them."""

    def font_digest(match):
        font = base64.b64decode(match[0], validate=True)
        return "sha256:" + hashlib.sha256(font).hexdigest()

    return re.sub(
        r'(?<=url\("data:font/woff2;base64,)[A-Za-z0-9+/=]+(?="\))', font_digest, svg
    )


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


fixtures = Path(__file__).parent / "fixtures"
has_linnet = globals().get("with_linnet", True)
model = fk.Model(fixtures / "scalars_2p_3p.json")
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
snapshot = diagram.to_json()
source = fixture_linnest(diagram)
empty = diagram.subgraph()
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
internal = diagram.filter(edge=lambda edge: not edge.is_external)
full = diagram.subgraph(half_edges=[half.id for half in diagram.half_edges])
paired = next(edge for edge in diagram.edges if not edge.is_external)
highlights = [
    internal,
    full,
    *(
        diagram.subgraph(half_edges=[half.id])
        for half in diagram.half_edges
        if half.edge.id == paired.id
    ),
    diagram.filter(edge=lambda edge: edge.is_external),
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
    if selected != full:
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
foreign = other.subgraph(half_edges=[half.id for half in other.half_edges])
for invalid in (foreign, object()):
    for render in (diagram.to_linnest, diagram.render, diagram.to_html):
        try:
            render(highlight=invalid)
        except (ValueError, TypeError):
            pass
        else:
            raise AssertionError("foreign highlight selections must be rejected")

if has_linnet:
    graph = diagram.to_graph()
    stale = graph.full_subgraph()
    assert fixture_linnest(diagram, highlight=stale) == fixture_linnest(
        diagram, highlight=full
    )
    graph.reverse_edge(0)
    for render in (diagram.to_linnest, diagram.render, diagram.to_html):
        try:
            render(highlight=stale)
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
    cross_section.subgraph(half_edges=[half.id for half in cross_section.half_edges]),
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

basketball = fk.FeynmanDiagram.from_dot(
    fk.Model.phi4(),
    "digraph basketball {" + 'a -> b [particle="phi"];' * 4 + "}",
)
assert basketball.loop_count == 3
outer = basketball.subgraph(edges=[0, 1, 2, 3])
middle = basketball.subgraph(edges=[0, 1, 2])
inner = basketball.subgraph(edges=[0, 1])
layers = [(outer, "acted-on"), (middle, "middle"), (inner, "treated")]
palette = ("#3b82f6", "#8b5cf6", "#e4504f")


def stroke_stacks(svg):
    stacks = {}
    for element in ET.fromstring(svg).iter():
        paint = element.get("stroke")
        if element.get("d") and paint in palette:
            stacks.setdefault(element.attrib["d"], []).append(
                (paint, float(element.attrib["stroke-width"]))
            )
    return stacks.values()


basketball_svg = basketball.render(
    highlight=layers, style=fk.DiagramStyle(show_particle=False), overlap="casings"
).to_svg()
stacks = [stack for stack in stroke_stacks(basketball_svg) if len(stack) == 3]
assert stacks
for stack in stacks:
    colors, widths = zip(*stack)
    assert colors == palette
    assert widths[0] > widths[1] > widths[2]


def check_palette(svg, prefix, scope):
    theme_rules = {}
    root = ET.fromstring(svg)
    css = root.find("{http://www.w3.org/2000/svg}style").text
    for selector, body in re.findall(r"([^{}]+)\{([^{}]*)\}", css):
        for theme in ("light", "dark"):
            if f'.{scope}[data-theme="{theme}"]' in selector:
                theme_rules[theme] = dict(
                    re.findall(rf"(--{prefix}-[a-z-]+):\s*(#[a-f0-9]+)", body)
                )
    for theme, expected in {
        "light": palette,
        "dark": ("#75a7ff", "#b197fc", "#ff7775"),
    }.items():
        assert [
            theme_rules[theme][f"--{prefix}-{role}"]
            for role in ("acted-on", "middle", "treated")
        ] == list(expected)


check_palette(basketball_svg, "feynkit", "feynkit-diagram-svg")
assert fixture_linnest(basketball, highlight=layers).count("#image(bytes(") == 1

legacy_style = fk.DiagramStyle(show_particle=False)
# SVGs from 055fd0c1 with embedded WOFF2 payloads replaced by their digests.
with patch.dict(sys.modules, {"typst": None}):
    legacy_svg = basketball.render(highlight=inner, style=legacy_style).to_svg()
assert (
    normalize_svg_fonts(legacy_svg).encode()
    == (fixtures / "single-amber-base-feynkit.svg").read_bytes()
), "FeynKit single-amber SVG differs from the base fixture"
assert "#77777773" in legacy_svg

unmuted = ET.fromstring(
    basketball.render(
        highlight=layers[1:],
        mute_outside=False,
    ).to_svg()
)
assert all(element.get("stroke") != "#77777773" for element in unmuted.iter())

owners = [(basketball, [outer, middle, inner])]
if has_linnet:
    assert hasattr(rendering, "Graph"), "expected the Linnet graph bindings"
    generic = basketball.to_graph()
    generic_outer = generic.full_subgraph()
    generic_middle = generic.subgraph(edges=[0, 1, 2])
    generic_inner = generic.subgraph(edges=[0, 1])
    with patch.dict(sys.modules, {"typst": None}):
        legacy_svg = generic_inner.render().to_svg()
    assert (
        normalize_svg_fonts(legacy_svg).encode()
        == (fixtures / "single-amber-base-linnet.svg").read_bytes()
    ), "Linnet single-amber SVG differs from the base fixture"
    owners.append((generic, [generic_outer, generic_middle, generic_inner]))

for owner, selections in owners:

    def svg(highlight):
        return owner.render(highlight=highlight, overlap="casings").to_svg()

    assert svg((selections[-1], "treated")) == svg([(selections[-1], "treated")])
    assert svg([selections[-1]]) == svg(selections[-1])
    for count, roles in (
        (2, ("acted-on", "treated")),
        (3, ("acted-on", "middle", "treated")),
    ):
        assert svg(selections[:count]) == svg(list(zip(selections[:count], roles)))
    a, b, c = selections
    assert svg([a, (b, "middle"), c]) == svg(
        [(a, "acted-on"), (b, "middle"), (c, "treated")]
    )
    assert svg([]) == owner.render().to_svg()
    try:
        svg(selections + selections[:1])
    except ValueError as error:
        assert "explicit highlight styles" in str(error)
    else:
        raise AssertionError("four unstyled layers must require explicit styles")

for owner, selections in owners:
    selected = (selections[-1], "treated")
    styled = owner.render(highlight=selected).to_svg()
    assert styled == owner.render(highlight=selected, overlap="casings").to_svg()
    assert "data-linnet-hull" not in styled
    assert "#77777773" in styled and "#e4504f" in styled
    assert not any(
        node.get("stroke") == "#77777773"
        for node in ET.fromstring(
            owner.render(highlight=selected, mute_outside=False).to_svg()
        ).iter()
    )

if has_linnet:
    generic_svg = generic.to_svg(
        highlight=[generic_outer, generic_middle, generic_inner]
    )
    check_palette(generic_svg, "linnet", "linnet-subgraph")
    disjoint_regions = [generic.subgraph(edges=[0, 1]), generic.subgraph(edges=[2, 3])]
    disjoint = stroke_stacks(
        generic.to_svg(highlight=disjoint_regions, overlap="casings")
    )
    assert disjoint and all(len(stack) == 1 for stack in disjoint)
    partly_empty = stroke_stacks(
        generic.to_svg(
            highlight=[generic.empty_subgraph(), generic_inner], overlap="casings"
        )
    )
    assert partly_empty and all(len(stack) == 1 for stack in partly_empty)
    default_hulls = ET.fromstring(generic.to_svg(highlight=disjoint_regions)).findall(
        ".//*[@data-linnet-hull='true']/{*}path"
    )
    assert {p.get("fill") for p in default_hulls} == {"#3b82f6", "#e4504f"}
    assert len({p.get("d") for p in default_hulls}) == len(default_hulls) == 2
    half_edge = generic.subgraph(half_edges=[generic.edge(0).source.index])
    for selection in (half_edge, (half_edge, "treated"), [half_edge, generic_inner]):
        try:
            generic.render(highlight=selection)
        except ValueError as error:
            assert "requires complete edges" in str(error)
        else:
            raise AssertionError("partial paired edges must be rejected in every layer")
    try:
        generic.render(highlight=inner)
    except TypeError as error:
        assert "symbolica.community.graph.Subgraph" in str(error)
        assert "graph.subgraph(" in str(error)
    else:
        raise AssertionError("a FeynKit region must not be accepted by a Linnet graph")
for highlight in ([(inner, "missing")], [(diagram.subgraph(edges=[0]), "treated")]):
    try:
        basketball.render(highlight=highlight)
    except ValueError:
        pass
    else:
        raise AssertionError("invalid highlight layers must be rejected")
print("installed ordered highlight checks passed")

# Requesting bands explicitly preserves a single explicit amber selection.
with patch.dict(sys.modules, {"typst": None}):
    amber = basketball.render(
        highlight=(inner, "amber"), style=legacy_style, overlap="bands"
    ).to_svg()
assert (
    normalize_svg_fonts(amber).encode()
    == (fixtures / "single-amber-base-feynkit.svg").read_bytes()
)


def svg_paths_at_width(svg, width):
    return [
        element
        for element in ET.fromstring(svg).iter()
        if element.tag.endswith("}path") and element.get("stroke-width") == width
    ]


many = basketball.render(
    highlight=[outer] * 8, style=legacy_style, overlap="bands"
).to_svg()
assert len(svg_paths_at_width(many, "2.4")) == 32
assert len(svg_paths_at_width(many, "2.5")) == 16
assert "--feynkit-slot-8:" in many
assert len({path.get("stroke") for path in svg_paths_at_width(many, "2.4")}) == 8
for overlap in ("bands", "bands-miniature"):
    single = basketball.render(
        highlight=[inner], style=legacy_style, overlap=overlap
    ).to_svg()
    assert len(svg_paths_at_width(single, "2.4")) == 2
    assert len(svg_paths_at_width(single, "2.5")) == 2
    assert any(path.get("stroke-dasharray") for path in ET.fromstring(single).iter())

tadpole = fk.FeynmanDiagram.from_dot(
    fk.Model.phi4(),
    """digraph tadpole {
        ext [style=invis];
        ext -> a [particle="phi"];
        a -> a [particle="phi"];
        a -> ext [particle="phi"];
    }""",
)
region = tadpole.subgraph(edges=[edge.id for edge in tadpole.edges])
svg = tadpole.render(
    overlap="bands", highlight=[region] * 2, style=legacy_style
).to_svg()
assert len(svg_paths_at_width(svg, "2.4")) == 6
assert len(svg_paths_at_width(svg, "2.5")) == 2
assert (
    fixture_linnest(basketball, highlight=layers, overlap="bands").count(
        "#image(bytes("
    )
    == 1
)
if has_linnet:
    svg = generic.to_svg(highlight=[generic_inner] * 2, overlap="bands")
    assert len(svg_paths_at_width(svg, "2.4")) == 4
print("installed parallel band checks passed")

# Reading size chooses H1; a stamp chooses plain recolouring.
hulls_svg = basketball.render(highlight=layers, style=legacy_style).to_svg()
assert len(ET.fromstring(hulls_svg).findall(".//*[@data-linnet-hull='true']")) == 3
stamp = basketball.render(
    highlight=layers, style=fk.DiagramStyle(show_particle=False, miniature=True)
).to_svg()
assert "data-linnet-hull" not in stamp
assert any(
    element.get("stroke") == "#e4504f" for element in ET.fromstring(stamp).iter()
)
assert (
    len(
        ET.fromstring(basketball.render(highlight=[outer] * 5).to_svg()).findall(
            ".//*[@data-linnet-hull='true']"
        )
    )
    == 5
)

label = 'γ₂ & "scheme" <local>'
description = "MS-bar & on-shell"
ct_layers = [
    (outer, "acted-on"),
    (
        inner,
        {
            "contracted": True,
            "label": label,
            "description": description,
        },
    ),
]
ct = basketball.render(highlight=ct_layers, style=legacy_style).to_svg()
root = ET.fromstring(ct)
for key in ("hull", "counterterm"):
    tagged = root.find(f".//*[@data-linnet-{key}='true'][@data-description]")
    if key == "hull":
        tagged = root.find(".//*[@data-linnet-hull='true'][@data-contracted='true']")
    assert tagged.get("data-label") == label
    assert tagged.get("data-description") == description
    assert tagged.find("{*}title").text == f"{label} — {description}"
assert root.find(".//*[@data-linnet-ghost='true']").get("opacity") == "0.35"
ct_stamp = basketball.render(
    highlight=ct_layers, style=fk.DiagramStyle(show_particle=False, miniature=True)
).to_svg()
assert (
    ET.fromstring(ct_stamp).find(".//*[@data-linnet-ghost='true']").get("opacity")
    == "0.62"
)
assert fixture_linnest(basketball, highlight=ct_layers).count("fill-rule") >= 2
print("installed H1 hull and CT2 counterterm checks passed")

mixed = fk.FeynmanDiagram.from_dot(
    fk.Model.qcd(),
    """digraph qse3 {
        i [style=invis]; o [style=invis];
        i -> a [particle="d"]; a -> b [particle="d"]; b -> o [particle="d"];
        a -> c [particle="g"]; c -> d [particle="ghG"]; d -> c [particle="ghG"];
        d -> e [particle="g"]; e -> f [particle="u"]; f -> e [particle="u"];
        f -> b [particle="g"];
    }""",
)
ms_bar = mixed.filter(edge=lambda edge: edge.particle_name == "ghG")
on_shell = mixed.filter(edge=lambda edge: edge.particle_name == "u")
mixed_layers = [
    (ms_bar, {"contracted": True, "scheme": "MS-bar", "label": "γ₁"}),
    (on_shell, {"contracted": True, "scheme": "on-shell", "label": "γ₂"}),
]
for miniature in (False, True):
    source = mixed.render(
        highlight=mixed_layers,
        style=fk.DiagramStyle(show_particle=False, miniature=miniature),
    ).to_svg()
    root = ET.fromstring(source)
    for scheme, fill, cross in (
        ("MS-bar", "#ffffff", "#e4504f"),
        ("on-shell", "#e4504f", "#ffffff"),
    ):
        glyph = root.find(
            f".//*[@data-linnet-counterterm='true'][@data-scheme='{scheme}']"
        )
        assert glyph.find("{*}circle").get("fill") == fill
        assert glyph.find("{*}path").get("stroke") == cross
        assert glyph.find("{*}title").text.endswith(scheme)
        hull = root.find(f".//*[@data-linnet-hull='true'][@data-scheme='{scheme}']")
        assert hull.find("{*}title").text.endswith(scheme)
try:
    mixed.render(highlight=(on_shell, {"contracted": True, "scheme": "wrong"}))
except ValueError as error:
    assert "counterterm scheme" in str(error)
else:
    raise AssertionError("unknown counterterm scheme must fail")
print("installed mixed CT2f scheme checks passed")

with patch.dict(sys.modules, {"typst": None}):
    single_lane = basketball.render(
        highlight=inner, style=legacy_style, overlap="lanes"
    ).to_svg()
assert (
    normalize_svg_fonts(single_lane).encode()
    == (fixtures / "single-amber-base-feynkit.svg").read_bytes()
)

label = 'γ₂ & <subgraph> "name"'
description = "local <integrated> & finite"
metro = basketball.render(
    highlight=[
        (outer, {"style": "acted-on", "label": label, "description": description}),
        (inner, "treated"),
    ],
    style=legacy_style,
    overlap="lanes",
).to_svg()
root = ET.fromstring(metro)
lanes = root.findall(".//*[@data-linnet-lane]")
assert {lane.get("data-linnet-lane") for lane in lanes} == {"0", "1"}
assert (
    root.find(".//*[@data-linnet-lane-crossings]").get("data-linnet-lane-crossings")
    == "0"
)
assert len(root.findall(".//*[@data-linnet-paper]")) == 2
assert "stroke-dasharray" in metro
for lane in lanes:
    path = lane.find("{*}path")
    assert path.get("fill") == "none"
    assert path.get("stroke-width") == "2.1"
    if lane.get("data-linnet-lane") == "0":
        assert lane.get("data-label") == label
        assert lane.get("data-description") == description
        assert label in lane.find("{*}title").text
        assert path.get("stroke") == "#3b82f6"
assert any("C" in lane.find("{*}path").get("d") for lane in lanes)

many_lanes = basketball.render(
    highlight=[outer] * 5, style=legacy_style, overlap="lanes"
).to_svg()
slots = ET.fromstring(many_lanes).findall(".//*[@data-linnet-lane]")
assert {lane.get("data-linnet-lane") for lane in slots} == set("01234")
assert len({lane.find("{*}path").get("stroke") for lane in slots}) == 5

# A loop has one carrier and one turn per slot. The selected external legs
# retain their particle ink and contribute no lanes.
loop_lanes = tadpole.render(
    overlap="lanes", highlight=[region, region], style=legacy_style
).to_svg()
loop_root = ET.fromstring(loop_lanes)
for slot in ("0", "1"):
    assert len(loop_root.findall(f".//*[@data-linnet-lane='{slot}']")) == 2
assert "data-linnet-lane" in fixture_linnest(
    basketball, highlight=layers, overlap="lanes"
)

if has_linnet:
    generic_lanes = generic.to_svg(
        highlight=[generic_outer, generic_inner], overlap="lanes"
    )
    assert ET.fromstring(generic_lanes).find(".//*[@data-linnet-lane]") is not None
    assert "[data-linnet-paper]" in generic_lanes
print("installed transit corner lane checks passed")

for mode in ("hulls", "recolour", "casings", "bands", "bands-miniature", "lanes"):
    drawing = basketball.render(
        highlight=ct_layers, overlap=mode, style=legacy_style
    ).to_svg()
    svg = ET.fromstring(drawing)
    hulls = svg.findall(".//*[@data-linnet-hull='true']")
    assert len(hulls) == (2 if mode == "hulls" else 1), mode
    assert svg.find(".//*[@data-linnet-counterterm='true']") is not None, mode
    if mode == "lanes":
        assert (
            svg.find(".//*[@data-linnet-lane='0']/{*}path[@stroke='#3b82f6']")
            is not None
        )
    elif mode != "hulls":
        assert svg.find(".//{*}path[@stroke='#3b82f6']") is not None, mode


def hull_path(layers):
    svg = basketball.render(
        highlight=layers, overlap="hulls", style=legacy_style
    ).to_svg()
    paths = ET.fromstring(svg).findall(".//*[@data-linnet-hull='true']/{*}path")
    assert len(paths) == 1
    return paths[0].get("d")


assert hull_path([(inner, "treated")]) == hull_path(
    [(inner, "treated"), (basketball.subgraph(), "middle")]
)
crowded = basketball.render(
    highlight=[basketball.subgraph(edges=[0])] * 128,
    overlap="lanes",
    style=legacy_style,
).to_svg()
assert len(ET.fromstring(crowded).findall(".//*[@data-linnet-lane]/{*}path")) == 128

for radius, step, layers_to_draw in (
    (1001, 0, [outer, inner]),
    (999, 2, [outer, inner]),
    (1000, 0, [(inner, {"contracted": True})]),
):
    try:
        basketball.render(
            highlight=layers_to_draw,
            style=fk.DiagramStyle(hull_radius=radius, hull_step=step),
        )
    except ValueError as error:
        assert "(0, 1000] pt" in str(error)
    else:
        raise AssertionError("invalid effective radius must fail before SVG generation")
for radius in (0.01, 0.000001):
    view = basketball.render(
        highlight=(inner, "treated"),
        overlap="hulls",
        style=fk.DiagramStyle(hull_radius=radius, show_particle=False),
    )
    try:
        paths = ET.fromstring(view.to_svg()).findall(".//*[@data-linnet-hull]/{*}path")
    except RuntimeError as error:
        assert "increase hull_radius or select a smaller region" in str(error)
    else:
        assert len(paths) == 1 and paths[0].get("d").count("M") == 1
assert (
    ET.fromstring(
        basketball.render(
            highlight=(inner, {"contracted": True, "scheme": "oN-sHeLl"}),
            style=legacy_style,
        ).to_svg()
    ).find(".//*[@data-linnet-counterterm='true'][@data-scheme='on-shell']")
    is not None
)
for owner, selections in owners:
    try:
        owner.render(
            highlight=(selections[-1], {"contracted": True, "style": "acted-on"})
        )
    except ValueError as error:
        assert "requires style 'treated'" in str(error)
    else:
        raise AssertionError("a contracted layer must use the treated paint")
print("installed highlight review regressions passed")
