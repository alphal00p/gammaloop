"""Community rendering uses no Python graph or Typst packages, including in WASM."""

import builtins
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
