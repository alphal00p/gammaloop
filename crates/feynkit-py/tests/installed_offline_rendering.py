"""Community rendering uses no Python graph or Typst packages, including in WASM."""

import builtins
import importlib.abc
import itertools
import json
import math
import re
import sys
import xml.etree.ElementTree as ET

from symbolica.community import hepkit as hep
from symbolica.community import tensor as spenso
from symbolica.core import S

original_import = builtins.__import__


class NoLinnet(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path, target=None):
        if fullname.split(".")[0] == "linnet":
            raise AssertionError("native graph analysis tried to import linnet")


no_linnet = NoLinnet()
sys.meta_path.insert(0, no_linnet)


def import_without_renderers(name, globals=None, locals=None, fromlist=(), level=0):
    if name.split(".")[0] in {"typst", "linnet", "gammaloop"}:
        raise AssertionError(f"rendering tried to import {name}")
    return original_import(name, globals, locals, fromlist, level)


builtins.__import__ = import_without_renderers
try:
    model = hep.Model.qcd()
    process = model.process(["g"], ["g"])
    assert "<svg" in process._repr_html_()
    inclusive = process.with_final_state_alternatives([["g"], ["u", "u~"]])
    combined = ET.fromstring(inclusive.render())
    figures = combined.findall("{http://www.w3.org/2000/svg}svg")
    assert len(figures) == 2
    for figure in figures:
        assert float(figure.attrib["width"]) == float(
            figure.attrib["viewBox"].split()[2]
        )
    # Process geometry must keep the large interaction blob and visible legs.
    # Check the rendered SVG: seed-only checks miss pins lost during layout.
    ns = "{http://www.w3.org/2000/svg}"
    for incoming, outgoing, config in (
        (["g"], ["g"], {}),
        (["u", "u~"], ["u", "u~"], {}),
        (["g", "g"], ["g", "g", "g"], {}),
        (["u"] * 4, [], {}),
        ([], ["u"] * 4, {}),
        (["u", "u~"], ["u", "u~"], {"drawing": {"node_radius": 5}}),
        (["u", "u~"], ["u", "u~"], {"style": {"node-style": {"radius": 6}}}),
    ):
        drawing = ET.fromstring(model.process(incoming, outgoing).render(config=config))
        blob = drawing.find(ns + "circle")
        cx, cy, radius = (float(blob.attrib[k]) for k in ("cx", "cy", "r"))
        assert radius >= 45, radius
        paths = [p for p in drawing.findall(ns + "path") if p.get("fill") == "none"]
        # A straight particle line can be emitted as several cubic segments.
        legs = []
        for path in paths:
            coords = list(map(float, re.findall(r"-?\d+(?:\.\d+)?", path.attrib["d"])))
            start, end = coords[:2], coords[-2:]
            if legs and math.dist(legs[-1][1], start) < 0.01:
                legs[-1][1] = end
            else:
                legs.append([start, end])
        assert len(legs) == len(incoming) + len(outgoing)
        tips = []
        for i, (start, end) in enumerate(legs):
            tip, contact = (start, end) if i < len(incoming) else (end, start)
            assert math.isclose(math.dist(contact, (cx, cy)), radius, abs_tol=0.01)
            assert math.dist(tip, contact) > 0.5 * radius
            tips.append(tip)
        if incoming and outgoing:
            left, right = tips[: len(incoming)], tips[len(incoming) :]
            assert all(x < cx for x, y in left)
            assert all(x > cx for x, y in right)
            for side in (left, right):
                assert all(a[1] < b[1] for a, b in itertools.pairwise(side))
        else:
            # One-sided processes surround the blob instead of occupying a side.
            xs, ys = zip(*tips)
            assert min(xs) < cx < max(xs) and min(ys) < cy < max(ys)
    ET.fromstring(model.process([], []).render())
    amplitude = process.generate_amplitude(loops=2, progress=None)
    assert len(amplitude.diagrams) == 48
    assert "<svg" in amplitude._repr_html_()
    diagram = amplitude.diagrams[0]
    ET.fromstring(diagram.render())
    ET.fromstring(
        diagram.render(
            config={
                "layouts": {"impred_steps": 2},
                "template_options": {"show-particle": False},
            }
        )
    )
    expression = diagram.numerator_expression()
    assert isinstance(expression, spenso.TensorExpression)
    assert "<math" in expression.to_html()
    ET.fromstring(expression.to_svg())
    # The remaining graph types also render without Typst graph plugins.
    ET.fromstring(process.render(config={"template_options": {"show-particle": False}}))
    network = spenso.TensorNetwork(spenso.TensorExpression(S("direct_svg_test::x") + 2))
    ET.fromstring(network.render(config={"title": "Native network"}))
    assert "#image(bytes(" in network.to_linnest()
    scalar = hep.Model.phi3()
    cross = scalar.process(["phi"], ["phi", "phi"]).generate_cross_section(
        loops=1, max_vertices=2, allow_self_loops=True, progress=None
    )[0]
    original = cross.to_json()
    # Selection and graph algorithms must work without the Python Linnet wheel,
    # including partial half-edges, physics callbacks, and returned result objects.
    full = cross.filter(edge=lambda edge: True)
    assert full.to_json() == original
    assert full.is_connected()
    assert len(full.connected_components()) == 1
    empty = cross.subgraph()
    assert empty.is_connected() and empty.connected_components() == []
    assert full.subgraph(empty).half_edge_indices() == []
    assert (~empty).half_edge_indices() == full.half_edge_indices()
    for half in cross.half_edges:
        assert isinstance(half, hep.DiagramHalfEdge)
        assert half.flow in {"source", "sink"}
        selected = cross.subgraph(half_edges=[half.id])
        assert selected.half_edge_indices() == [half.id]
        assert selected.denominator_expression() == spenso.TensorExpression(1)
        assert cross.filter(
            half_edge=lambda item: item.id == half.id
        ).half_edge_indices() == [half.id]
        assert half.vertex in {half.edge.source, half.edge.target}
    for vertex in cross.vertices:
        assert (
            cross.filter(node=lambda item: item.id == vertex.id).half_edge_indices()
            == cross.subgraph(nodes=[vertex.id]).half_edge_indices()
        )
    internal = cross.filter(edge=lambda edge: not edge.is_external)
    assert not internal.filter(edge=lambda edge: edge.is_external)
    assert isinstance(internal.boundary(), hep.Subgraph)
    assert isinstance(full.bridges(), hep.Subgraph)
    cycles, covered = full.cycle_basis()
    # Sewn initial-state carriers participate in structural cycles but are
    # excluded from the physical momentum-basis loop count.
    assert len(cycles) == len(full.edges) - len(full.vertices) + 1
    assert all(isinstance(cycle, hep.Subgraph) for cycle in cycles)
    assert covered.loop_count == 0 and covered.is_connected()
    assert len(covered.edges) == len(full.vertices) - 1
    forests = full.all_spanning_forests()
    assert forests and all(
        forest.is_connected() and len(forest.edges) == len(full.vertices) - 1
        for forest in forests
    )
    assert full.all_bonds() and full.all_bonds(min_size=99) == []
    partitions = full.all_cuts([0], [1])
    assert partitions and all(isinstance(cut, hep.CutPartition) for cut in partitions)
    for cut in partitions:
        assert cut.source_side.vertices[0].id == 0
        assert cut.target_side.vertices[0].id == 1
        assert {half.edge.id for half in cut.boundary_left.half_edges} == {
            half.edge.id for half in cut.boundary_right.half_edges
        }
        assert not (cut.boundary_left & cut.boundary_right)
    for traverse in (full.depth_first_traverse, full.breadth_first_traverse):
        tree = traverse(0)
        assert isinstance(tree, hep.TraversalTree)
        assert [vertex.id for vertex in tree.nodes] == [0, 1]
        assert tree.parent(0) is None and tree.parent(1).id == 0
        assert [vertex.id for vertex in tree.children(0)] == [1]
        assert [vertex.id for vertex in tree.ancestors(1)] == [0]
        assert len(tree.subgraph.edges) == len(tree.nodes) - 1
        assert tree.covers(full).half_edge_indices() == full.half_edge_indices()
        non_tree = [
            half
            for half in cross.half_edges
            if half.id not in tree.subgraph.half_edge_indices()
        ]
        assert non_tree
        fundamental = tree.fundamental_cycle(non_tree[0].id)
        assert fundamental.is_connected() and len(fundamental.edges) == 2
        assert tree.fundamental_cycle(tree.subgraph.half_edge_indices()[0]) is None
        assert (
            traverse(
                0,
                include=next(half.id for half in cross.half_edges if half.vertex == 0),
            )
            .nodes[0]
            .id
            == 0
        )
    foreign = hep.FeynmanDiagram.from_json(scalar, original).subgraph(nodes=[0])
    for operation, error_type in (
        (lambda: cross.subgraph(object()), TypeError),
        (lambda: cross.subgraph(foreign), ValueError),
        (lambda: tree.covers(foreign), ValueError),
        (lambda: cross.subgraph(nodes=[99]), IndexError),
        (lambda: cross.subgraph(edges=[99]), IndexError),
        (lambda: cross.subgraph(half_edges=[99]), IndexError),
        (lambda: cross.depth_first_traverse(99), IndexError),
        (lambda: cross.depth_first_traverse(0, include=99), IndexError),
        (lambda: empty.depth_first_traverse(0), ValueError),
        (lambda: full.all_bonds(min_size=0), ValueError),
        (lambda: full.all_bonds(min_size=2, max_size=1), ValueError),
        (lambda: full.all_cuts([], [1]), ValueError),
        (lambda: full.all_cuts([0], [0]), ValueError),
        (lambda: full.all_cuts([99], [1]), IndexError),
    ):
        try:
            operation()
        except error_type:
            pass
        else:
            raise AssertionError("invalid native graph input was accepted")
    assert cross.to_json() == original
    # Zero-crown interactions survive selection, component/forest results, and
    # singleton traversals even though they have no half-edge bit to select.
    zero_model = json.loads(scalar.to_json())
    zero_model["vertex_rules"][0]["particles"] = []
    zero_model["lorentz_structures"][0]["spins"] = []
    zero_scalar = hep.Model.from_json(json.dumps(zero_model))
    isolated = hep.FeynmanDiagram.from_dot(
        zero_scalar, "digraph isolated { a [num=2]; b [num=3]; }"
    )
    assert not isolated.is_connected()
    components = isolated.connected_components()
    assert [component.isolated_node_indices() for component in components] == [[0], [1]]
    assert isolated.subgraph(
        nodes=[0]
    ).numerator_expression() == spenso.TensorExpression(2)
    assert isolated.filter(
        node=lambda vertex: vertex.id == 1
    ).isolated_node_indices() == [1]
    assert isolated.all_spanning_forests()[0].isolated_node_indices() == [0, 1]
    for traverse in (isolated.depth_first_traverse, isolated.breadth_first_traverse):
        tree = traverse(0)
        assert [vertex.id for vertex in tree.nodes] == [0]
        assert tree.subgraph.isolated_node_indices() == [0]
        assert tree.covers(isolated.subgraph(nodes=[0, 1])).isolated_node_indices() == [
            0
        ]
        assert tree.parent(0) is None and tree.children(0) == tree.ancestors(0) == []
    tadpole = (
        hep.Model.phi4()
        .process(["phi"], ["phi"])
        .generate_diagrams(
            loops=1, max_vertices=1, allow_self_loops=True, progress=None
        )[0]
    )
    assert len(tadpole.cycle_basis()[0]) == 1
    for traverse in (tadpole.depth_first_traverse, tadpole.breadth_first_traverse):
        tree = traverse(0)
        assert [vertex.id for vertex in tree.nodes] == [0]
        internal_half = next(
            half for half in tadpole.half_edges if not half.edge.is_external
        )
        assert len(tree.fundamental_cycle(internal_half.id).edges) == 1
    sewn_identities = None
    for split in (False, True):
        svg = cross.render(config={"template_options": {"split-initial-state": split}})
        root = ET.fromstring(svg)
        targets = [n for n in root.iter() if "data-linnet-kind" in n.attrib]
        identities = {
            (n.attrib["data-linnet-kind"], n.attrib["data-linnet-id"]) for n in targets
        }
        assert identities
        if sewn_identities is None:
            sewn_identities = identities
        else:
            assert identities == sewn_identities
        assert cross.to_json() == original
        edges = {e.id: e for e in cross.edges}
        for node in targets:
            detail = json.loads(node.attrib["data-linnet-detail"])
            if node.attrib["data-linnet-kind"] != "node":
                edge = edges[detail["edge"]]
                assert (detail["source"], detail["sink"]) == (edge.source, edge.target)
                assert detail["particle"] == edge.particle_name
    custom = json.loads(scalar.to_json())
    custom["particles"][0]["typstname"] = "alpha_1"
    custom["particles"][0]["antitypstname"] = "alpha_1"
    for parameter in custom["parameters"]:
        parameter["typstname"] = "beta_2"
        parameter.pop("texname", None)
    labelled = hep.Model.from_json(json.dumps(custom))
    assert labelled.particle("phi").typstname == "alpha_1"
    assert json.loads(labelled.to_json())["particles"][0]["typstname"] == "alpha_1"
    ET.fromstring(labelled.process(["phi"], ["phi"]).render())
    parameter = labelled.parameter("g").symbol
    assert "beta_2" in parameter.to_typst()
    assert "g" in spenso.TensorExpression(parameter).to_latex()
    assert "mitex" not in parameter.to_typst()
    ET.fromstring(spenso.TensorExpression(parameter).to_svg())
    # Every authored Standard Model label must compile without package access.
    standard = hep.Model.standard_model()
    ET.fromstring(standard.process([p.name for p in standard.particles], []).render())
    all_parameters = sum(p.symbol for p in standard.parameters)
    ET.fromstring(spenso.TensorExpression(all_parameters).to_svg())
    assert not hasattr(sys.modules["symbolica.community"], "linnet")
    assert "linnet" not in sys.modules
    assert "typst" not in sys.modules
    print(
        "Offline QCD process, 48 two-loop diagrams, cross sections, networks, and native Typst labels rendered"
    )
finally:
    builtins.__import__ = original_import
    sys.meta_path.remove(no_linnet)
