"""Community rendering uses no Python graph or Typst packages, including in WASM."""

import builtins
import json
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
