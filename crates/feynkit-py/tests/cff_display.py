"""Verify rich-display payloads against native CFF objects, including filtered results."""

import json
import os
import re
import xml.etree.ElementTree as ET
from pathlib import Path


def check_display(diagram, result):
    before = result.to_expression()
    html = result._repr_html_()
    assert result.to_expression() == before
    assert "</script><script>throw 'unsafe'</script>" not in html
    data = json.loads(
        re.search(
            r"<script type=\"application/json\" data-cff>(.*?)</script>", html
        ).group(1)
    )
    assert data["name"] == diagram.name
    assert len(data["orientations"]) == len(result.orientations)
    for native, displayed in zip(
        result.orientations, data["orientations"], strict=True
    ):
        assert displayed["id"] == native.id
        assert {int(k): v for k, v in displayed["directions"].items()} == dict(
            native.edge_orientations
        )
        expected = [
            [(s.kind, s.index) for s in term if s.kind != "unit"]
            for term in native.denominator_products()
            if all(s.kind != "infinite" for s in term)
        ]
        actual = [
            [
                (data["surfaces"][key]["kind"], data["surfaces"][key]["index"])
                for key in term
            ]
            for term in displayed["terms"]
        ]
        assert actual == expected
    for native, displayed in zip(result.surfaces, data["surfaces"], strict=True):
        assert displayed["v"] == native.vertices
        assert displayed["e"] == native.positive_energies
        assert displayed["negative"] == native.negative_energies
        assert displayed["q"] == [list(pair) for pair in native.external_shift]
    svg = ET.fromstring(
        re.search(r"<template data-drawing>(.*?)</template>", html, re.DOTALL).group(1)
    )
    edges = {edge.id: edge for edge in diagram.edges}
    carriers = [node for node in svg.iter() if "data-linnet-carrier" in node.attrib]
    assert len(carriers) == len(diagram.half_edges)
    owners = set()
    for carrier in carriers:
        detail = json.loads(carrier.attrib["data-linnet-detail"])
        edge = edges[detail["edge"]]
        assert detail["node"] == (
            edge.source if detail["flow"] == "source" else edge.target
        )
        assert carrier.attrib["tabindex"] == "-1"
        owners.add(detail["half-edge"])
    assert len(owners) == len(diagram.half_edges)
    assert {
        int(n.attrib["data-linnet-node"])
        for n in svg.iter()
        if "data-linnet-node" in n.attrib
    } == {v.id for v in diagram.vertices}
    return html


def check_cff_displays(fk, diagrams):
    assert not hasattr(fk, "CffResult")
    pages = []
    for diagram in diagrams:
        result = diagram.cross_free_family()
        assert type(result) is fk.CrossFreeFamily
        assert type(result).__name__ == "CrossFreeFamily"
        assert type(result).__module__ == "symbolica.community.hepkit"
        assert repr(result).startswith("CrossFreeFamily(")
        assert not hasattr(diagram, "build_cff")
        assert not hasattr(diagram.subgraph(nodes=[0]), "build_cff")
        pages.append(check_display(diagram, result))
        if diagram.loop_count == 3:
            assert len(result.orientations) == 686
            assert len(result) == 3432
            for coefficient in result.pole_coefficients(
                result.raised_surface_groups()[0]
            ):
                pages.append(check_display(diagram, coefficient))
            generator = fk.CffGenerator()
            generator.fix_orientation(diagram.internal_edges[0].id, reversed=True)
            pages.append(check_display(diagram, generator.generate(diagram)))
            generator = fk.CffGenerator()
            generator.contract_edge(diagram.internal_edges[0].id)
            pages.append(check_display(diagram, generator.generate(diagram)))
            pages.append(
                check_display(diagram, diagram.subgraph(nodes=[0]).cross_free_family())
            )
            cycle = fk.CffGenerator()
            for edge, reverse in [(0, True), (1, False), (2, True), (3, False)]:
                cycle.fix_orientation(edge, reversed=reverse)
            zero = cycle.generate(diagram)
            assert not zero.orientations
            pages.append(check_display(diagram, zero))

    # Collections share the same theme/frame helper; keep their preview contract.
    collection = (
        fk.Model.phi3()
        .process(["phi"], ["phi", "phi"])
        .generate_diagrams(loops=0, max_vertices=1, progress=None)
        ._repr_html_()
    )
    assert (
        'class="feynkit-collection" data-feynkit-notebook data-linnet-frame-owner'
        in collection
    )

    if directory := os.environ.get("CFF_DISPLAY_TEST_OUTPUT"):
        target = Path(directory)
        target.mkdir(parents=True, exist_ok=True)
        (target / "collection.html").write_text(collection)
        for index, html in enumerate(pages):
            (target / f"cff-{index}.html").write_text(html)
    print(f"CFF explorer: {len(pages)} native displays passed")
