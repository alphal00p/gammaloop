"""Run against the installed community library; export browser fixtures to argv[1]."""

import json
import re
import sys
from html.parser import HTMLParser
from pathlib import Path

import marimo as mo
from symbolica import E, S
from symbolica.community import hepkit as hep


def payload(source):
    return json.loads(
        re.search(
            r'<script type="application/json" data-ltd>(.*?)</script>', source
        ).group(1)
    )


class FrameParser(HTMLParser):
    def handle_starttag(self, tag, attrs):
        if tag == "iframe":
            self.frame = dict(attrs)


def check(diagram, name, directory):
    ltd = diagram.integrate_energy(method="ltd")
    assert type(ltd) is hep.LtdRepresentation
    assert repr(ltd).startswith("LtdRepresentation(residues=")
    before = ltd.to_expression(expand_surfaces=False)
    source = ltd._repr_html_()
    assert ltd.to_expression(expand_surfaces=False) == before
    data = payload(source)
    assert len(data["residues"]) == len(ltd)
    assert data["edges"] == [e.id for e in diagram.internal_edges]
    assert [
        data["edges"][i] for i in data["reference_chords"]
    ] == diagram.loop_momentum_basis.loop_edges
    parser = FrameParser()
    parser.feed(mo.as_html(ltd).text)
    assert "allow-scripts" in parser.frame.get("sandbox", "allow-scripts")
    assert payload(parser.frame["srcdoc"]) == data

    # Compare the complete displayed scalar products against native Symbolica
    # residues. This catches lost coefficients, multiplicities and ID remapping.
    ose, momentum, cind = S("gammalooprs::OSE", "gammalooprs::Q", "spenso::cind")
    internal = [ose(i) for i in data["edges"]]
    external = [
        momentum(i, cind(0)) for i in diagram.loop_momentum_basis.external_edges
    ]
    surfaces = [
        sum((E(c) * internal[i] for i, c in s["expression"]["internal"]), E("0"))
        + sum(
            (E(c) * external[i] for i, c in s["expression"]["external"]),
            E(s["expression"]["constant"]),
        )
        for s in data["surfaces"]
    ]
    for residue in ltd.residues:
        assert set(residue.pole_signs) == set(residue.cut_edges)
        for edge, sign in residue.pole_signs.items():
            assert residue.energy_map[edge] == sign * ose(edge)
        for edge, pair in residue.surface_pairs.items():
            assert (
                (pair.minus.to_expression() - residue.energy_map[edge] + ose(edge))
                .expand()
                .is_zero()
            )
            assert (
                (pair.plus.to_expression() - residue.energy_map[edge] - ose(edge))
                .expand()
                .is_zero()
            )
    for r, native in zip(
        data["residues"],
        [residue.to_expression() for residue in ltd.residues],
        strict=True,
    ):
        reconstructed = E("0")
        for term in r["terms"]:
            coefficient = E(term["coefficient"])
            for edge in term["energies"]:
                coefficient /= 2 * internal[edge]
            for chain in term["chains"]:
                value = coefficient
                for surface in chain:
                    value /= surfaces[surface]
                reconstructed += value
        assert (reconstructed - native).together().is_zero(), (name, r["id"])
    (directory / f"ltd-{name}.html").write_text(source)
    return ltd, data


def main(directory):
    directory.mkdir(parents=True, exist_ok=True)
    model = hep.Model.phi3()
    triangle = (
        model.process(["phi"], ["phi", "phi"])
        .generate_diagrams(loops=1, max_vertices=3, maximum_bridges=0, progress=None)
        .diagrams[0]
    )
    ltd, data = check(triangle, "triangle", directory)
    assert len(ltd) == 3
    # Dispatch retains the concrete native types, algebra and constructor
    # conventions, rather than wrapping the results in an opaque common object.
    for method, expected in [
        ("cff", triangle.integrate_energy(method="cff")),
        ("ltd", ltd),
    ]:
        result = triangle.integrate_energy(method=method)
        assert type(result) is type(expected)
        assert result.to_expression() == expected.to_expression()
    try:
        triangle.integrate_energy(method="unknown")
    except ValueError as error:
        assert "'cff' or 'ltd'" in str(error)
    else:
        raise AssertionError("an unknown energy representation was accepted")
    cff = triangle.integrate_energy(method="cff").to_expression(expand_surfaces=True)
    momentum, cind = S("gammalooprs::Q", "spenso::cind")
    basis = triangle.loop_momentum_basis
    # CFF keeps all external energies; LTD uses the stored routing, which
    # eliminates the dependent external energy by momentum conservation.
    for edge in triangle.external_edges:
        routed = sum(
            (
                c * momentum(i, cind(0))
                for i, c in zip(
                    basis.external_edges, basis.edge_signatures[edge.id].external
                )
            ),
            E("0"),
        )
        cff = cff.replace(momentum(edge.id, cind(0)), routed)
    # Both conversions include the same contour and on-shell convention.
    assert (ltd.to_expression() - cff).together().is_zero()

    triple = hep.FeynmanDiagram.from_dot(
        model, (Path(__file__).parent / "fixtures/ltd_three_loop.dot").read_text()
    )
    _, data = check(triple, "three-loop", directory)
    assert triple.loop_count == 3 and len(data["residues"]) == 56
    reordered = triple.with_loop_momentum_edges(
        list(reversed(triple.loop_momentum_basis.loop_edges))
    )
    _, other = check(reordered, "reordered", directory)
    assert data["reference_chords"] == list(reversed(other["reference_chords"]))
    assert [r["signs"] for r in data["residues"]] != [
        r["signs"] for r in other["residues"]
    ]

    tree = (
        model.process(["phi", "phi"], ["phi", "phi"])
        .generate_diagrams(loops=0, progress=None)
        .diagrams[0]
    )
    check(tree, "tree", directory)
    contact = (
        model.process(["phi"], ["phi", "phi"])
        .generate_diagrams(loops=0, progress=None)
        .diagrams[0]
    )
    check(contact, "contact", directory)
    tadpole = (
        model.process(["phi"], [])
        .generate_diagrams(
            loops=1,
            max_vertices=1,
            allow_self_loops=True,
            allow_zero_flow_edges=True,
            tadpoles=None,
            zero_snails=None,
            self_energy=None,
            factorized_loop_topologies_count_range=None,
            maximum_bridges=None,
            progress=None,
        )
        .diagrams[0]
    )
    check(tadpole, "tadpole", directory)

    scalar = hep.Model(Path(__file__).parent / "fixtures/scalars_2p_3p.json")
    raised = next(
        d
        for d in scalar.process(
            ["scalar_0"], ["scalar_0"], vertex_allow=["V_2_SCALAR_00", "V_3_SCALAR_000"]
        )
        .generate_diagrams(
            loops=1, max_vertices=3, allow_self_loops=False, progress=None
        )
        .diagrams
        if len(d.internal_edges) == 3
    )
    _, data = check(raised, "raised", directory)
    assert any(
        len(t["energies"]) > len(set(t["energies"]))
        for r in data["residues"]
        for t in r["terms"]
    )

    try:
        triangle.subgraph(edges=[triangle.internal_edges[0].id]).integrate_energy(
            method="ltd"
        )
    except hep.DiagramError:
        pass
    else:
        raise AssertionError("partial subgraph silently used its parent diagram")
    unsafe = hep.FeynmanDiagram.from_dot(
        model,
        (Path(__file__).parent / "fixtures/ltd_three_loop.dot")
        .read_text()
        .replace("triple_box", "unsafe</script><script>throw 1</script>"),
    )
    source = unsafe.integrate_energy(method="ltd")._repr_html_()
    assert "</script><script>throw 1</script>" not in source
    (directory / "ltd-multiple.html").write_text(ltd._repr_html_() + ltd._repr_html_())
    print(
        "LtdRepresentation: native algebra, ordered routing, repeated poles and Marimo dispatch passed"
    )


if __name__ == "__main__":
    main(Path(sys.argv[1]))
