"""Check the public CFF/LTD decomposition, shared energies and contour convention."""

import json
import random
import re
from pathlib import Path

from symbolica import E, Replacement, S
from symbolica.community import hepkit as hep
from symbolica.community.tensor import TensorExpression

FIXTURES = Path(__file__).parent / "fixtures"
OSE, Q, K, P, CIND = S(
    "gammalooprs::OSE",
    "gammalooprs::Q",
    "gammalooprs::K",
    "gammalooprs::P",
    "spenso::cind",
)


def payload(obj, method):
    return json.loads(
        re.search(
            rf'<script type="application/json" data-{method}>(.*?)</script>',
            obj._repr_html_(),
        ).group(1)
    )


def check(diagram):
    cff = diagram.integrate_energy(method="cff")
    ltd = diagram.integrate_energy(method="ltd")
    assert type(cff) is hep.CffRepresentation
    assert type(ltd) is hep.LtdRepresentation
    assert not hasattr(diagram, "cross_free_family")
    assert not hasattr(diagram, "loop_tree_duality")
    assert not hasattr(hep, "CffGenerator")
    assert not hasattr(hep, "LoopTreeDuality")
    assert cff.diagram.to_dot() == ltd.diagram.to_dot() == diagram.to_dot()
    assert list(cff.on_shell_energies) == list(ltd.on_shell_energies)
    basis = diagram.loop_momentum_basis
    for edge, energy in cff.on_shell_energies.items():
        assert energy.symbol == ltd.on_shell_energies[edge].symbol == OSE(edge)
        assert energy.to_expression() == ltd.on_shell_energies[edge].to_expression()
        # Independently reconstruct routed Cartesian components from the basis.
        signature = basis.edge_signatures[edge]
        mass = next(e.particle.mass for e in diagram.internal_edges if e.id == edge)
        squared = mass**2 + sum(
            (
                sum(c * K(e, CIND(axis)) for e, c in enumerate(signature.loops))
                + sum(c * P(e, CIND(axis)) for e, c in enumerate(signature.external))
            )
            ** 2
            for axis in range(1, 4)
        )
        assert (energy.to_expression() ** 2 - squared).expand().is_zero()
    for rep in (cff, ltd):
        assert type(rep.to_expression()) is TensorExpression
        assert rep.to_expression().rank == 0
        for surface in rep.surfaces:
            assert type(surface) is hep.EnergySurface
            reconstructed = (
                surface.constant
                + sum(c * OSE(e) for e, c in surface.energy_coefficients.items())
                + sum(c * Q(e, CIND(0)) for e, c in surface.external_shift.items())
            )
            assert (surface.to_expression() - reconstructed).expand().is_zero()
        replacements = [Replacement(s.symbol, s.to_expression()) for s in rep.surfaces]
        assert (
            rep.to_expression(expand_surfaces=False).replace_multiple(replacements)
            == rep.to_expression()
        )
    assert ltd.report.residues == len(ltd.residues)
    assert ltd.report.trees == len({tuple(r.tree_edges) for r in ltd.residues})
    residues = [r.to_expression() for r in ltd.residues]
    assert sum(residues, TensorExpression(0)) == ltd.to_expression()
    for residue in ltd.residues:
        assert len(residue.cut_edges) == diagram.loop_count
        assert set(residue.tree_edges).isdisjoint(residue.cut_edges)
        assert set(residue.tree_edges) | set(residue.cut_edges) == set(
            ltd.on_shell_energies
        )
        for edge, sign in residue.pole_signs.items():
            assert sign in (-1, 1)
            assert residue.energy_map[edge] == sign * OSE(edge)
        for edge, pair in residue.surface_pairs.items():
            assert edge in residue.tree_edges
            assert (
                (pair.minus.to_expression() - (residue.energy_map[edge] - OSE(edge)))
                .expand()
                .is_zero()
            )
            assert (
                (pair.plus.to_expression() - (residue.energy_map[edge] + OSE(edge)))
                .expand()
                .is_zero()
            )
    prefactor = E("1")
    for energy in cff.on_shell_energies.values():
        prefactor /= -2 * energy.symbol
    # Compare every family's factors without a costly global rational expansion.
    family_expressions = []
    for orientation in cff.orientations:
        assert orientation.edge_signs == {
            edge: {"default": 1, "reversed": -1, "undirected": None}[direction]
            for edge, direction in orientation.edge_orientations
        }
        terms = []
        for family in orientation.families:
            assert type(family) is hep.CrossFreeFamily
            reconstructed = TensorExpression(prefactor)
            for factor in family.factors:
                reconstructed /= factor.to_expression(expand_surfaces=False)
            actual = family.to_expression(expand_surfaces=False)
            assert reconstructed == actual
            terms.append(actual)
        if len(terms) < 10:
            assert (
                (
                    sum(terms, TensorExpression(0))
                    - orientation.to_expression(expand_surfaces=False)
                )
                .together()
                .is_zero()
            )
        family_expressions.extend(terms)
    rng = random.Random(914)
    for sample in range(3):
        replacements = [
            Replacement(OSE(e), E(str(rng.randrange(100, 10000))))
            for e in cff.on_shell_energies
        ]
        for edge in diagram.external_edges:
            value = sum(
                c * (i + sample + 1)
                for i, c in enumerate(basis.edge_signatures[edge.id].external)
            )
            replacements.append(Replacement(Q(edge.id, CIND(0)), E(str(value))))
        left = cff.to_expression().replace_multiple(replacements).together()
        right = ltd.to_expression().replace_multiple(replacements).together()
        assert left == right, (diagram.name, sample)
    # Sum decomposition is certified in local surface coordinates as well.
    surface_values = [
        Replacement(s.symbol, E(str(i + 17))) for i, s in enumerate(cff.surfaces)
    ]
    assert (
        (
            sum(family_expressions, TensorExpression(0))
            - cff.to_expression(expand_surfaces=False)
        )
        .replace_multiple(surface_values)
        .together()
        .is_zero()
    )
    if cff.orientations and cff.orientations[-1].families:
        orientation = cff.orientations[-1]
        family = orientation.families[-1]
        view = payload(orientation, "cff")
        assert view["scope"] == "orientation"
        assert [o["id"] for o in view["orientations"]] == [orientation.id]
        assert len(view["orientations"][0]["terms"]) == len(orientation.families)
        view = payload(family, "cff")
        assert view["scope"] == "family"
        assert len(view["orientations"]) == 1
        assert view["orientations"][0]["family_ids"] == [family.id]
        assert len(view["orientations"][0]["terms"]) == 1
    if ltd.residues:
        view = payload(ltd.residues[-1], "ltd")
        assert view["scope"] == "residue"
        assert [r["id"] for r in view["residues"]] == [ltd.residues[-1].id]
    print(
        diagram.name,
        "shared energies, factors, pole signs, decomposition and CFF=LTD passed",
        flush=True,
    )


def main():
    model = hep.Model.phi3()
    triangle = (
        model.process(["phi"], ["phi", "phi"])
        .generate_diagrams(
            loops=1,
            max_vertices=3,
            maximum_bridges=0,
            progress=None,
        )
        .diagrams[0]
    )
    cases = [
        triangle,
        hep.FeynmanDiagram.from_dot(
            model,
            """digraph bubble {
        edge [particle="phi"]; ext [style=invis];
        ext -> a; a -> b; a -> b [lmb_id=0]; b -> ext;
    }""",
        ),
    ]
    for dot in [
        """digraph tree { edge [particle="phi"]; ext [style=invis];
           ext -> a; ext -> a; a -> b; b -> ext; b -> ext; }""",
        """digraph double_box { edge [particle="phi"]; ext [style=invis];
           ext -> a; ext -> d; c -> ext; f -> ext;
           a -> b; b -> c; d -> e; e -> f;
           a -> d; b -> e [lmb_id=0]; c -> f [lmb_id=1]; }""",
        (FIXTURES / "ltd_three_loop.dot").read_text(),
        """digraph tetrahedron { edge [particle="phi"];
           a -> b; a -> c; a -> d;
           b -> c [lmb_id=0]; b -> d [lmb_id=1]; c -> d [lmb_id=2]; }""",
    ]:
        cases.append(hep.FeynmanDiagram.from_dot(model, dot))
    for diagram in cases:
        check(diagram)
    check(
        cases[-2].with_loop_momentum_edges(
            list(reversed(cases[-2].loop_momentum_basis.loop_edges))
        )
    )
    for options in [
        {"max_orientations": 1},
        {"fixed_orientations": {}},
        {"contracted_edges": []},
        {"initial_state_edges": []},
    ]:
        try:
            triangle.integrate_energy(method="ltd", **options)
        except ValueError as error:
            assert "only" in str(error)
        else:
            raise AssertionError("CFF options silently accepted for LTD")


if __name__ == "__main__":
    main()
