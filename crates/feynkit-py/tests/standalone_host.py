"""Shared bindings work in a host that installs no symbolica package aliases."""

import builtins
import json
import re
import sys

# Supplied by the Rust host without importing any component module.
core = globals()["core"]
graph = globals()["graph"]
tensor = globals()["tensor"]
hep = globals()["hep"]
linnet_python = globals()["linnet_python"]

assert not any(
    name == "symbolica" or name.startswith("symbolica.") for name in sys.modules
)
original_import = builtins.__import__


def without_symbolica(name, globals=None, locals=None, fromlist=(), level=0):
    if name == "symbolica" or name.startswith("symbolica."):
        raise AssertionError(f"native component tried to import {name}")
    return original_import(name, globals, locals, fromlist, level)


builtins.__import__ = without_symbolica
try:
    model = hep.Model.phi3()
    diagram = hep.FeynmanDiagram.from_dot(
        model,
        """digraph bubble {
            ext [style=invis];
            ext -> a [particle="phi"];
            a -> b [particle="phi", lmb_id=0];
            a -> b [particle="phi"];
            b -> ext [particle="phi"];
        }""",
    )
    numerator = diagram.numerator_expression()
    assert isinstance(numerator, tensor.TensorExpression)
    assert isinstance(numerator.to_expression(), core.Expression)
    basis = diagram.loop_momentum_basis
    rules = basis.momentum_replacements()
    assert len(rules) == len(diagram.edges)
    assert all(isinstance(rule, core.Replacement) for rule in rules)
    q, cind = core.S("gammalooprs::Q", "spenso::cind")
    expression = sum(q(edge.id, cind(0)) ** 2 for edge in diagram.edges)
    assert expression.replace_multiple(rules) == basis.route_expression(expression)

    # A numerator with time components also exercises CFF energy-map formatting.
    cff = diagram.integrate_energy(
        method="cff", numerator=q(diagram.internal_edges[0].id, cind(0)) ** 2
    )
    html = cff._repr_html_()
    assert "data-cff" in html and "<merror" not in html
    data = json.loads(
        re.search(
            r'<script type="application/json" data-cff>(.*?)</script>', html
        ).group(1)
    )
    contributions = [c for o in data["orientations"] for c in o["contributions"]]
    assert contributions
    assert any(c["energy_map"] for c in contributions)
    assert all(
        isinstance(text, str)
        for c in contributions
        for text in c["energy_map"].values()
    )
    assert data["energy_html"]
    assert all("<math" in source for source in data["energy_html"].values())
    assert isinstance(cff.to_expression(), tensor.TensorExpression)

    if linnet_python:
        exported = diagram.to_graph()
        assert isinstance(exported, graph.Graph)
        assert diagram.to_graph() is exported
        selection = diagram.subgraph(edges=[diagram.internal_edges[0].id])
        assert isinstance(selection.to_graph(), graph.Graph)
        assert selection.to_graph() is exported
        assert len(exported.edges()) == len(diagram.edges)
finally:
    builtins.__import__ = original_import

assert not any(
    name == "symbolica" or name.startswith("symbolica.") for name in sys.modules
)
