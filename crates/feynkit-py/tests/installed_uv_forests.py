"""Installed-wheel checks of UV spinneys, woods, forests and counterterms.

Counts agree with the feynkit-uv Rust tests, which agree with GammaLoop's
hedge-poset regressions for the same diagrams.
"""

import xml.etree.ElementTree as ET
from pathlib import Path

from symbolica.community import hepkit as fk

SVG = "{http://www.w3.org/2000/svg}"


def phi4(dot):
    return fk.FeynmanDiagram.from_dot(fk.Model.phi4(), dot)


def edges_and_degrees(wood):
    return [(tuple(s.edges), s.loop_count, s.degree) for s in wood]


def check_view(view, marker, count):
    svg = view.to_svg()
    root = ET.fromstring(svg)
    assert len(root.findall(f".//*[@{marker}]")) == count, (marker, count)
    assert "prefers-color-scheme:dark" in svg
    assert "<a " not in svg, "stamp graphs carry no inspection links"
    titles = [title.text for title in root.iter(f"{SVG}title")]
    assert len(titles) >= count
    assert "<figure" in view.to_html()
    return svg, titles


# Basketball: 6 one-loop pairs, 4 two-loop triples and the full graph.
basketball = phi4(
    "digraph basketball { edge [particle=phi]; a -> b; a -> b; a -> b; a -> b; }"
)
wood = fk.Wood(basketball)
assert len(wood) == 11 and repr(wood) == "Wood(11 spinneys, max loops 3)"
loops = [spinney.loop_count for spinney in wood]
assert loops == sorted(loops) and [loops.count(n) for n in (1, 2, 3)] == [6, 4, 1]
assert all(s.degree == 4 * s.loop_count - 2 * len(s.edges) for s in wood)
assert all(s.is_divergent and s.approximation == "MUV" for s in wood)
assert wood[-1] == wood[10] and wood[-1].edges == [0, 1, 2, 3]
try:
    wood[11]
except IndexError:
    pass
else:
    raise AssertionError("wood indices are bounded")
covers = wood.covers()
assert len(covers) == 12 + 4
assert all(small.relation(large) == "contained in" for small, large in covers)
forests = wood.forests()
assert len(forests) == 46
assert forests[0].is_empty and forests[0].cover is None and forests[0].sign == 1
assert max(forest.depth for forest in forests) == 3
for forest in forests:
    assert forest.sign == (-1) ** forest.depth
    assert sum(map(len, forest.levels)) == forest.depth == len(forest.spinneys)
    assert all(a.is_compatible(b) for a in forest.spinneys for b in forest.spinneys)
    for child in forest.children:
        assert child.depth == forest.depth + 1 and forest in child.parents
    if not forest.is_empty:
        assert forest.cover == max(forest.spinneys, key=lambda s: s.loop_count)
assert len(forests[0].children) == 11 and len(set(forests)) == 46
check_view(wood.render(), "data-feynkit-spinney", 12)
operation_svg, titles = check_view(wood.operation_graph(), "data-feynkit-forest", 46)
assert len(operation_svg) < 400_000, "stamps are shared through <use>"
assert operation_svg.count("<symbol ") < operation_svg.count("<use ")
assert any("depth 3, sign −" in title for title in titles)
try:
    wood.operation_graph(stamp_size=0)
except fk.UvError as error:
    assert "stamp size" in str(error)
else:
    raise AssertionError("stamp sizes must be positive")
print("installed basketball wood and forest checks passed")

# Sunrise: three overlapping bubbles inside the quadratically divergent graph.
sunrise = phi4(
    """digraph sunrise { edge [particle=phi]; ext [style=invis];
    ext -> a; a -> b; a -> b; a -> b; b -> ext; }"""
)
wood = fk.Wood(sunrise)
assert edges_and_degrees(wood) == [
    ((0, 1), 1, 0),
    ((0, 2), 1, 0),
    ((1, 2), 1, 0),
    ((0, 1, 2), 2, 2),
]
bubble, other, _, full = wood
assert bubble.relation(other) == "overlapping" and not bubble.is_compatible(other)
assert full.relation(bubble) == "contains" and bubble.relation(bubble) == "equal"
assert (full.superficial_degree, full.label) == (2, "γ{0,1,2}")
assert full.diagram.id == sunrise.id and len(full.components) == 1
assert (
    full.subgraph.half_edge_indices()
    == full.subgraph.original.subgraph(edges=[0, 1, 2]).half_edge_indices()
)
assert len(wood.forests()) == 8
assert all(forest.depth <= 2 for forest in wood.forests())
hep_spinney = fk.Spinney(sunrise.subgraph(edges=[0, 1, 2]))
assert hep_spinney == full and hash(hep_spinney) == hash(full)
for edges, reason in (([1], "union of cycles"), ([3], "union of cycles")):
    try:
        fk.Spinney(sunrise.subgraph(edges=edges))
    except fk.UvError as error:
        assert reason in str(error), error
    else:
        raise AssertionError(f"{edges} is not a spinney")
# Externals are dropped: the internal lines of a selection form the spinney.
assert fk.Spinney(sunrise.subgraph(edges=[0, 1, 3])) == bubble

view = full.render()
spinney_svg = ET.fromstring(view.to_svg())
assert any(
    "2 loops, degree 2 (superficial 2), MUV approximation" in (title.text or "")
    for title in spinney_svg.iter(f"{SVG}title")
), "hover carries the spinney details"
nested = next(forest for forest in wood.forests() if forest.depth == 2)
hulls = ET.fromstring(nested.render().to_svg()).findall(
    ".//*[@data-linnet-hull='true']"
)
assert len(hulls) == 2, "a nested forest stacks one hull per spinney"
stamp = nested.render(
    style=fk.DiagramStyle(show_particle=False, miniature=True)
).to_svg()
assert "data-linnet-hull" not in stamp and "#e4504f" in stamp
lanes = nested.render(overlap="lanes").to_svg()
assert "data-linnet-lane" in lanes
print("installed sunrise spinney and forest checks passed")

# Counterterms: ghost in place with one circled cross per outermost spinney,
# open for MS-bar, filled on shell and in a dashed ring while unspecified. The
# default scheme follows the approximation; MUV leaves it unspecified.
for scheme, attribute, fill, dashed in (
    ("MS-bar", "MS-bar", "#ffffff", False),
    ("on-shell", "on-shell", "#e4504f", False),
    (None, None, "#ffffff", True),
    ("unspecified", None, "#ffffff", True),
):
    counterterm = nested.counterterm(scheme)
    assert counterterm.scheme == (attribute or "unspecified")
    assert [s.loop_count for s in counterterm.contracted] == [2, 1]
    root = ET.fromstring(counterterm.render().to_svg())
    glyphs = root.findall(".//*[@data-linnet-counterterm='true']")
    assert len(glyphs) == 1, "only the outermost spinney is crossed"
    (glyph,) = glyphs
    assert glyph.get("data-scheme") == attribute
    circle = glyph.find(f"{SVG}circle")
    assert circle.get("fill") == fill
    assert (circle.get("stroke-dasharray") is not None) == dashed
    title = glyph.find(f"{SVG}title").text
    # The hover lists the chain, inner to outer, with the scheme.
    assert title.startswith(f"{nested.spinneys[0].label} ⊂ {full.label}"), title
    assert "contracted" in title
    assert title.endswith(attribute) if attribute else "MS-bar" not in title
    hulls = root.findall(".//*[@data-linnet-hull='true']")
    assert len(hulls) == 2, "the nested spinney shows as its own hull"
    assert root.find(".//*[@data-linnet-ghost='true']") is not None
local = fk.Counterterm(bubble, scheme="ms-BAR")
assert repr(local) == "Counterterm(γ{0,1}, scheme=MS-bar)"
assert local.with_scheme("on-shell").scheme == "on-shell"
assert fk.Counterterm([full, bubble]).scheme == "unspecified"
for invalid, exception, text in (
    ([bubble, other], fk.UvError, "are overlapping"),
    ([], fk.UvError, "at least one"),
    (wood.forests()[0], fk.UvError, "at least one"),
    (object(), TypeError, "Spinney"),
):
    try:
        fk.Counterterm(invalid)
    except exception as error:
        assert text in str(error), error
    else:
        raise AssertionError(f"invalid counterterm {invalid!r} was accepted")
try:
    bubble.counterterm("minimal")
except ValueError as error:
    assert "MS-bar" in str(error)
else:
    raise AssertionError("unknown schemes must fail")
print("installed counterterm checks passed")

# Nested vertex: the bubble of the quartic eye sits inside the full graph.
eye = phi4(
    """digraph eye { edge [particle=phi]; ext [style=invis]; ext -> a; ext -> a;
    a -> b; a -> c; b -> c; b -> c; b -> ext; c -> ext; }"""
)
forests = fk.Wood(eye).forests()
assert len(forests) == 4
try:
    fk.Spinney(eye.subgraph(edges=[0, 1, 2]))
except fk.UvError as error:
    assert "convergent" in str(error)
else:
    raise AssertionError("the triangle of the eye is convergent")
assert repr(forests[-1]) == "Forest(K_{Γ|a}, k3, depth=2, sign=+1)"

# Disjoint spinneys: two tadpoles form a disconnected spinney and a product forest.
spectacles = fk.FeynmanDiagram.from_dot(
    fk.Model.phi3(), "digraph { edge [particle=phi]; a -> a; a -> b; b -> b; }"
)
wood = fk.Wood(spectacles)
union = next(spinney for spinney in wood if len(spinney.components) == 2)
assert union.loop_count == 2 and union.edges == [0, 2]
left, right = (spinney for spinney in wood if spinney != union)
assert left.relation(right) == "disjoint" and left.is_compatible(right)
forests = wood.forests()
assert len(forests) == 4
product = next(forest for forest in forests if forest.is_disjoint_union)


def labels(spinneys):
    return sorted(spinney.label for spinney in spinneys)


# A product forest lists each factor, never their union.
assert [labels(level) for level in product.levels] == [["γ{0}", "γ{2}"]]
assert labels(product.spinneys) == ["γ{0}", "γ{2}"]
assert (product.depth, product.sign, product.cover) == (2, 1, union)
assert labels(product.counterterm("MS-bar").contracted) == ["γ{0}", "γ{2}"]
assert labels(fk.Counterterm([left, right]).contracted) == ["γ{0}", "γ{2}"]
assert repr(product) == "Forest(K_{a}K_{b}, k3, depth=2, sign=+1)"
glyphs = ET.fromstring(product.counterterm("on-shell").render().to_svg()).findall(
    ".//*[@data-linnet-counterterm='true']"
)
assert sorted(glyph.get("data-label") for glyph in glyphs) == ["γ{0}", "γ{2}"]
operation_svg, titles = check_view(wood.operation_graph(), "data-feynkit-forest", 4)
assert "γ{0,2}" not in operation_svg
# Manually given tadpoles regain their product, as Wood(diagram) has it.
given = fk.Wood.from_spinneys(spectacles, [left, right])
assert len(given) == 3 and len(given.forests()) == 4
assert any(forest.is_disjoint_union for forest in given.forests())

# Two quartic bubbles joined by a line: the full graph is convergent, so the
# only two-element forest is the product of the bubbles.
pair = phi4(
    """digraph pair { edge [particle=phi]; ext [style=invis];
    ext -> a; ext -> a; a -> b; a -> b; b -> c; b -> ext;
    c -> d; c -> d; c -> ext; d -> ext; d -> ext; }"""
)
wood = fk.Wood(pair)
bubbles = [spinney for spinney in wood if len(spinney.components) == 1]
assert labels(bubbles) == ["γ{0,1}", "γ{4,5}"]
assert bubbles[0].relation(bubbles[1]) == "disjoint"
forests = wood.forests()
assert len(forests) == 4
product = next(forest for forest in forests if forest.is_disjoint_union)
assert [labels(level) for level in product.levels] == [["γ{0,1}", "γ{4,5}"]]
assert labels(product.counterterm().contracted) == labels(
    fk.Counterterm(bubbles).contracted
)
hulls = ET.fromstring(product.render().to_svg()).findall(
    ".//*[@data-linnet-hull='true']"
)
assert sorted(hull.get("data-label") for hull in hulls) == ["γ{0,1}", "γ{4,5}"]
assert len({hull.find(f"{SVG}path").get("stroke") for hull in hulls}) == 1
assert len(fk.Wood.from_spinneys(pair, bubbles).forests()) == 4

# Convergent factors: neither the prism's two triangles nor a convergent
# triangle with a divergent tadpole form a spinney, as Wood(diagram) agrees.
prism = phi4(
    """digraph prism { edge [particle=phi]; ext [style=invis];
    a -> b; b -> c; c -> a; d -> e; e -> f; f -> d; a -> d; b -> e; c -> f;
    ext -> a; ext -> b; ext -> c; d -> ext; e -> ext; f -> ext; }"""
)
kite = fk.FeynmanDiagram.from_dot(
    fk.Model.phi3(),
    """digraph kite { edge [particle=phi]; ext [style=invis];
    a -> b; b -> c; c -> a; a -> t; t -> t; ext -> b; c -> ext; }""",
)
assert [spinney.edges for spinney in fk.Wood(kite)] == [[6]]
for diagram, edges in ((prism, [0, 2, 4, 6, 9, 13]), (kite, [0, 2, 3, 6])):
    try:
        fk.Spinney(diagram.subgraph(edges=edges))
    except fk.UvError as error:
        assert "convergent" in str(error), error
    else:
        raise AssertionError(f"{edges} has a convergent component")
assert fk.Spinney(kite.subgraph(edges=[6])).is_divergent

# Half of a tadpole line is not a selection of that line.
tadpole = left.subgraph.half_edge_indices()
try:
    fk.Spinney(spectacles.subgraph(half_edges=tadpole[:1]))
except fk.UvError as error:
    assert "one half of internal edge 0" in str(error), error
else:
    raise AssertionError("half a line is not a spinney")
assert fk.Spinney(spectacles.subgraph(half_edges=tadpole)) == left

# Bubbles sharing a vertex touch: disjoint subgraphs of the forest formula share
# no vertex, so the s-channel chain has no product forest.
chain = phi4(
    """digraph chain { edge [particle=phi]; ext [style=invis]; ext -> a; ext -> a;
    a -> b; a -> b; b -> c; b -> c; c -> ext; c -> ext; }"""
)
first, second, _ = fk.Wood(chain)
assert first.relation(second) == "touching" and not first.is_compatible(second)
try:
    fk.Counterterm([first, second])
except fk.UvError as error:
    assert "are touching" in str(error)
else:
    raise AssertionError("touching spinneys cannot be contracted together")
forests = fk.Wood(chain).forests()
assert len(forests) == 6
assert all(len(level) == 1 for forest in forests for level in forest.levels)
print("installed nested and disjoint forest checks passed")

# QCD: numerators enter the classifying degree once Feynman rules are applied.
rainbow = fk.FeynmanDiagram.from_dot(
    fk.Model.qcd(),
    """digraph rainbow { ext [style=invis];
    ext -> a [particle="d"]; a -> b [particle="d"]; b -> c [particle="d"];
    c -> d [particle="d"]; d -> ext [particle="d"];
    a -> d [particle="g"]; b -> c [particle="g"]; }""",
).apply_feynman_rules()
assert edges_and_degrees(fk.Wood(rainbow)) == [((2, 3), 1, 1), ((0, 1, 2, 3, 4), 2, 1)]
assert len(fk.Wood(rainbow).forests()) == 4

# User-given spinneys: validated against the destination diagram.
pair = fk.Spinney(basketball.subgraph(edges=[0, 1]))
triple = fk.Spinney(basketball.subgraph(edges=[0, 1, 2]))
assert len(fk.Wood.from_spinneys(basketball, [pair, triple, pair]).forests()) == 4


def cubic_bubble(reversed_line):
    second = "b -> a" if reversed_line else "a -> b"
    return fk.FeynmanDiagram.from_dot(
        fk.Model.phi3(),
        "digraph bubble { ext [style=invis]; edge [particle=phi]; "
        f"a -> b [id=0]; {second} [id=1]; ext -> a [id=2]; b -> ext [id=3]; }}",
    )


aligned, reversed_bubble = cubic_bubble(False), cubic_bubble(True)
loop = fk.Spinney(aligned.subgraph(edges=[0, 1]))
assert len(fk.Wood.from_spinneys(aligned, [loop])) == 1
# A structurally identical diagram routes the spinney identically.
assert len(fk.Wood.from_spinneys(cubic_bubble(False), [loop])) == 1
for destination, reason in (
    (reversed_bubble, "its routed sub-basis differs"),
    (basketball, "its width differs"),
):
    try:
        fk.Wood.from_spinneys(destination, [loop])
    except fk.UvError as error:
        assert reason in str(error), error
        assert isinstance(error, fk.FeynkitError)
    else:
        raise AssertionError(f"foreign spinney accepted ({reason})")
try:
    fk.Wood(basketball.subgraph(edges=[0, 1]))
except fk.DiagramError:
    pass
else:
    raise AssertionError("woods need a complete diagram")
try:
    loop.relation(pair)
except ValueError as error:
    assert "different diagrams" in str(error)
else:
    raise AssertionError("relations need one diagram")
print("installed user-given and foreign spinney checks passed")

# The forest-formula notation: Op[m(S|G) · inner] terms with cofactors m(Γ|γ), the whole
# graph named Γ and the other spinneys a, b, …, and stored quantities k0, k1, … in topological order.
vacuum = fk.FeynmanDiagram.from_dot(
    fk.Model.phi3(), "digraph sunrise { edge [particle=phi]; a -> b; a -> b; a -> b; }"
)
wood = fk.Wood(vacuum)
assert [wood.name(spinney) for spinney in wood] == ["a", "b", "c", "Γ"]
terms = sorted(wood.forests(), key=lambda forest: forest.k)
assert [(forest.k, forest.sign, forest.notation()) for forest in terms] == [
    (0, 1, "K[∅]"),
    (1, -1, "K_{Γ}"),
    (2, -1, "K_{a}"),
    (3, -1, "K_{b}"),
    (4, -1, "K_{c}"),
    (5, 1, "K_{Γ|a}"),
    (6, 1, "K_{Γ|b}"),
    (7, 1, "K_{Γ|c}"),
]
assert terms[5].notation("T^{3D}") == "T^{3D}_{Γ|a}"
expansion = ET.fromstring(wood.expansion().to_svg())
assert len(expansion.findall(".//*[@data-feynkit-term]")) == 8
assert len(expansion.findall(f".//{SVG}use")) == 1 + 2 * 4 + 3 * 3
operation = wood.operation_graph(op="δ").to_svg()
assert operation.count("data-feynkit-k=") == 8 and ">δ<" in operation
term = ET.fromstring(terms[5].term().to_svg())
assert len(term.findall(f".//{SVG}use")) == 3
assert "#image(bytes(" in wood.expansion_typst()

dotted = fk.FeynmanDiagram.from_dot(
    fk.Model(Path(__file__).parent / "fixtures" / "scalars_2p_3p.json"),
    """digraph dotted { edge [particle="scalar_0"];
    a -> b; a -> b; a -> c; c -> b; }""",
)
wood = fk.Wood(dotted)
assert wood.k_table() == "\n".join(
    ["k0 = K[∅]", "k1 = K_{Γ}", "k2 = K_{a}", "k3 = K_{Γ|a}[k2]"]
)
heavy, light = (wood.render(stamp_weight=weight).to_svg() for weight in (2.0, 1.0))
assert heavy != light
for weight in (0.0, float("nan")):
    try:
        wood.render(stamp_weight=weight)
    except fk.UvError as error:
        assert "stamp weight" in str(error)
    else:
        raise AssertionError("stamp weights must be positive")

# User-chosen names: by spinney or current name; unknown, invalid or duplicate
# names are rejected and leave the wood unchanged.
named = wood.renamed({wood[0]: "x", "Γ": "G"})
assert named.k_table().splitlines()[3] == "k3 = K_{G|x}[k2]"
assert wood.k_table().splitlines()[3] == "k3 = K_{Γ|a}[k2]"
assert repr(sorted(named.forests(), key=lambda f: f.k)[3]).startswith(
    "Forest(K_{G|x}, k3"
)
assert "K_{G|x} = k3" in named.operation_graph().to_svg()
for names, message in (
    ({"q": "x"}, "no spinney of this wood is named"),
    ({wood[0]: "Γ"}, "reserved for the whole graph"),
    ({wood[0]: "x", "a": "y"}, "more than one new name"),
    ({wood[0]: "∅"}, "must be nonempty"),
    ({wood[0]: "x·y"}, "must be nonempty"),
    ({wood[0]: "k1"}, "k<digits>"),
    ({wood[0]: "x\x01y"}, "control characters"),
    ({wood[0]: "x y"}, "must be nonempty"),
    ({fk.Wood(vacuum)[0]: "x"}, "is not a spinney of this wood"),
):
    try:
        wood.renamed(names)
    except ValueError as error:
        assert message in str(error), str(error)
    else:
        raise AssertionError(f"{names} must be rejected")
print("installed forest notation checks passed")

# Text written into SVG is escaped, and control characters cannot reach it.
assert ET.fromstring(wood.expansion(op="K\x01<&>").to_svg()) is not None

# Review regressions: woods of given spinneys, counterterm normalisation, forest
# order, API consistency, glyph placement and nesting in miniatures.
import math
import re

# A spinney taken from a diagram with Feynman rules applied is convergent in the
# bare diagram, which power counts its own numerators.
triangle = fk.FeynmanDiagram.from_dot(
    fk.Model.qcd(),
    """digraph triangle { i [style=invis]; j [style=invis]; k [style=invis];
    a -> b [particle="u"]; b -> c [particle="u"]; c -> a [particle="u"];
    i -> a [particle="g"]; j -> b [particle="g"]; k -> c [particle="g"]; }""",
)
assert len(fk.Wood(triangle)) == 0
ruled = fk.Wood(triangle.apply_feynman_rules())
assert len(ruled) == 1 and ruled[0].is_divergent
try:
    fk.Wood.from_spinneys(triangle, [ruled[0]])
except fk.UvError as error:
    assert "superficially convergent in this diagram" in str(error), str(error)
else:
    raise AssertionError("a spinney convergent here must be rejected")

# A disconnected given spinney multiplies with a disjoint one: ∅, X, Y and X·Y.
chain3 = phi4(
    """digraph chain3 { edge [particle=phi]; ext [style=invis];
    a -> a [id=0]; b -> b [id=1]; c -> c [id=2]; a -> b [id=3]; b -> c [id=4];
    ext -> a; c -> ext; }"""
)
classified = fk.Wood(chain3)
assert len(classified.forests()) == 8
x = fk.Spinney(chain3.subgraph(edges=[0, 1]))
y = fk.Spinney(chain3.subgraph(edges=[2]))
assert x.relation(y) == "disjoint"
given = fk.Wood.from_spinneys(chain3, [x, y])
assert len(given.forests()) == 4
product = next(forest for forest in given.forests() if forest.is_disjoint_union)
assert (product.sign, labels(product.spinneys)) == (1, ["γ{0,1}", "γ{2}"])
assert given.k_table().splitlines()[-1] == "k3 = k2·k1"
# Giving every classified spinney back reproduces the classified wood.
again = fk.Wood.from_spinneys(chain3, list(classified))
assert again.k_table() == classified.k_table()

# Repeated inputs count once.
assert (
    len(fk.Wood.from_spinneys(spectacles, [left] * 300 + [right] * 300).forests()) == 4
)

# A disconnected spinney contracts as its components, as its product forest
# does; a spinney together with a union holding it is no forest.
by_union = fk.Counterterm([union], "MS-bar")
assert labels(by_union.contracted) == ["γ{0}", "γ{2}"]
assert by_union == fk.Counterterm(
    next(f for f in fk.Wood(spectacles).forests() if f.is_disjoint_union), "MS-bar"
)
assert len({by_union, fk.Counterterm([left, right], scheme="MS-bar")}) == 1
assert by_union != by_union.with_scheme("on-shell")
assert (
    len(
        ET.fromstring(by_union.render().to_svg()).findall(
            ".//*[@data-linnet-counterterm='true']"
        )
    )
    == 2
)
try:
    fk.Counterterm([left, union])
except fk.UvError as error:
    assert "more than once" in str(error)
else:
    raise AssertionError("a tadpole and the union holding it are no forest")

# Forests come in k order.
basket = fk.Wood(basketball)
assert [forest.k for forest in basket.forests()] == list(range(46))


# Each circled cross lies in its own spinney's painted hull, never in a hole,
# at reading size and in miniature.
def polygons(d):
    rings, ring = [], []
    for command, x, y in re.findall(r"([MLZ])\s*(-?[\d.]+)?[ ,]?(-?[\d.]+)?", d):
        if command == "M":
            if ring:
                rings.append(ring)
            ring = [(float(x), float(y))]
        elif command == "L":
            ring.append((float(x), float(y)))
        elif ring:
            rings.append(ring)
            ring = []
    if ring:
        rings.append(ring)
    return rings


def painted(rings, x, y):
    parity = False
    for ring in rings:
        for (x1, y1), (x2, y2) in zip(ring, ring[1:] + ring[:1]):
            if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1) + x1:
                parity = not parity
    return parity


def glyph_placements(svg):
    root = ET.fromstring(svg)
    hulls = {
        group.get("data-label"): polygons(group.find(f"{SVG}path").get("d"))
        for group in root.iter()
        if group.get("data-linnet-hull") == "true"
        and group.get("data-contracted") == "true"
    }
    for glyph in root.iter():
        if glyph.get("data-linnet-counterterm") != "true":
            continue
        match = re.match(
            r"translate\(([-\d.]+) ([-\d.]+)\)(?: scale\(([-\d.]+)\))?",
            glyph.get("transform"),
        )
        cx, cy, scale = float(match[1]), float(match[2]), float(match[3] or 1)
        rings = hulls[glyph.get("data-label")]
        radius = 5.925 * scale
        ring = sum(
            painted(rings, cx + radius * math.cos(t), cy + radius * math.sin(t))
            for t in (2 * math.pi * i / 48 for i in range(48))
        )
        yield painted(rings, cx, cy), ring / 48


placements = []
for diagram in (
    phi4(
        "digraph bubble { edge [particle=phi]; ext [style=invis]; ext -> a; ext -> a; a -> b; a -> b; b -> ext; b -> ext; }"
    ),
    vacuum,
    basketball,
    spectacles,
    phi4(
        """digraph pair { edge [particle=phi]; ext [style=invis];
        ext -> a; ext -> a; a -> b; a -> b; b -> c; b -> ext;
        c -> d; c -> d; c -> ext; d -> ext; d -> ext; }"""
    ),
    phi4(
        "digraph chain { edge [particle=phi]; ext [style=invis]; ext -> a; ext -> a; a -> b; a -> b; b -> c; b -> c; c -> ext; c -> ext; }"
    ),
    chain3,
):
    for forest in fk.Wood(diagram).forests()[1:12]:
        counterterm = forest.counterterm("MS-bar")
        for style in (None, fk.DiagramStyle(miniature=True)):
            placements += [
                (style is None, centre, ring)
                for centre, ring in glyph_placements(
                    counterterm.render(style=style).to_svg()
                )
            ]
assert placements and all(centre for _, centre, _ in placements), placements
# The placement keeps the largest share of the ring inside. At reading size at
# least a quarter of it stays in the hull; a full-size cross on a miniature
# tadpole is larger than the tadpole's whole hull.
reading = [ring for full_size, _, ring in placements if full_size]
assert min(reading) >= 0.25, min(reading)

# Nested contracted levels stay visible in miniatures, one hull per level.
deepest = max(basket.forests(), key=lambda forest: forest.depth)
for style in (None, fk.DiagramStyle(miniature=True)):
    hulls = ET.fromstring(deepest.counterterm().render(style=style).to_svg()).findall(
        ".//*[@data-linnet-hull='true']"
    )
    assert len(hulls) == deepest.depth == 3
print("installed review regressions passed")
