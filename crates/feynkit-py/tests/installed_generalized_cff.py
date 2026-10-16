"""Check numerator-aware CFF against explicit one-loop contour residues."""

import json
import re
import sys
from pathlib import Path

from symbolica import E, Replacement, S
from symbolica.community import hepkit as hep
from symbolica.community.tensor import TensorExpression

Q, OSE, CIND, T = S("gammalooprs::Q", "gammalooprs::OSE", "spenso::cind", "cff_test_t")


def one_loop_oracle(diagram, numerator):
    """Close below: at each simple pole k=E-shift/sign, use -Res."""
    basis = diagram.loop_momentum_basis
    momenta = {}
    for edge in diagram.internal_edges:
        signature = basis.edge_signatures[edge.id]
        sign = signature.loops[0]
        assert sign in (-1, 1)
        shift = sum(
            (
                c * Q(e, CIND(0))
                for e, c in zip(basis.external_edges, signature.external, strict=True)
            ),
            E("0"),
        )
        momenta[edge.id] = (sign, shift)
    numerator = numerator.replace_multiple(
        [
            Replacement(Q(e, CIND(0)), sign * T + shift)
            for e, (sign, shift) in momenta.items()
        ]
    )
    answer = E("0")
    for cut, (sign, shift) in momenta.items():
        root = OSE(cut) - sign * shift
        value = -numerator.replace(T, root) / (2 * OSE(cut))
        for edge, (other_sign, other_shift) in momenta.items():
            if edge != cut:
                value /= (other_sign * root + other_shift) ** 2 - OSE(edge) ** 2
        answer += value
    return answer


def check(diagram, numerator, output=None):
    result = diagram.integrate_energy(method="cff", numerator=numerator)
    assert type(result) is hep.CffRepresentation
    assert result.numerator == numerator
    total = result.to_expression()
    assert type(total) is TensorExpression and total.rank == 0
    orientations = [o.to_expression() for o in result.orientations]
    assert (sum(orientations, TensorExpression(0)) - total).expand().is_zero()
    for orientation in result.orientations:
        terms = []
        for family in orientation.families:
            value = family.coefficient * family.numerator
            for factor in family.factors:
                value /= factor.to_expression()
            assert (family.to_expression() - value).together().is_zero()
            terms.append(value)
            assert set(family.energy_map) == {e.id for e in diagram.internal_edges}
        assert (
            (sum(terms, TensorExpression(0)) - orientation.to_expression())
            .together()
            .is_zero()
        )
    if output is not None:
        source = result._repr_html_()
        data = json.loads(
            re.search(
                r'<script type="application/json" data-cff>(.*?)</script>', source
            ).group(1)
        )
        for surface, displayed in zip(result.surfaces, data["surfaces"], strict=True):
            assert displayed["support"] == [
                edge
                for edge, coefficient in surface.energy_coefficients.items()
                if coefficient != E("0")
            ]
            if surface.origin == "helper" or surface.numerator_only:
                assert not displayed["v"]
        reconstructed = E("0")
        for orientation in data["orientations"]:
            for path, contribution in zip(
                orientation["terms"], orientation["contributions"], strict=True
            ):
                value = E(contribution["coefficient"]) * E(contribution["numerator"])
                for edge in contribution["energies"]:
                    value /= 2 * OSE(edge)
                for surface in path:
                    value /= E(data["surfaces"][surface]["symbol"])
                reconstructed += value
        reconstructed *= E(data["normalization"])
        assert (
            (reconstructed - result.to_expression(expand_surfaces=False))
            .together()
            .is_zero()
        )
        output.write_text(source)
    return result


def main(directory):
    directory.mkdir(parents=True, exist_ok=True)
    model = hep.Model.phi3()
    triangle = (
        model.process(["phi"], ["phi", "phi"])
        .generate_diagrams(loops=1, max_vertices=3, maximum_bridges=0, progress=None)
        .diagrams[0]
    )
    edge = triangle.internal_edges[0].id
    energy = Q(edge, CIND(0))
    for degree in range(5):
        numerator = (energy + E("2/3")) ** degree
        result = check(triangle, numerator, directory / f"cff-degree-{degree}.html")
        assert result.energy_degree_bounds == ({edge: degree} if degree else {})
        oracle = one_loop_oracle(triangle, numerator)
        assert (result.to_expression() - oracle).together().is_zero(), degree
        print(
            f"degree {degree}: exact contour, family sum and display reconstruction passed",
            flush=True,
        )
    double = hep.FeynmanDiagram.from_dot(
        model,
        """digraph double_box {
        edge [particle=phi]; ext [style=invis];
        ext -> a; ext -> d; c -> ext; f -> ext;
        a -> b; b -> c; d -> e; e -> f;
        a -> d; b -> e [lmb_id=0]; c -> f [lmb_id=1];
    }""",
    )
    triple = hep.FeynmanDiagram.from_dot(
        model, (Path(__file__).parent / "fixtures/ltd_three_loop.dot").read_text()
    )
    for diagram in (double, triple):
        selected = diagram.loop_momentum_basis.loop_edges[0]
        numerator = Q(selected, CIND(0)) ** 2
        result = diagram.integrate_energy(method="cff", numerator=numerator)
        ltd = diagram.integrate_energy(method="ltd")
        expected = sum(
            (
                residue.to_expression()
                * numerator.replace_multiple(
                    [
                        Replacement(Q(edge, CIND(0)), value)
                        for edge, value in residue.energy_map.items()
                    ]
                )
                for residue in ltd.residues
            ),
            TensorExpression(0),
        )
        for seed in range(3):
            replacements = [
                Replacement(OSE(edge.id), E(str(17 + 13 * edge.id + seed)))
                for edge in diagram.internal_edges
            ]
            replacements += [
                Replacement(Q(edge.id, CIND(0)), E(str((seed + 1) * (edge.id + 1))))
                for edge in diagram.external_edges
            ]
            assert (
                (result.to_expression() - expected)
                .replace_multiple(replacements)
                .together()
                .is_zero()
            ), (diagram.loop_count, seed)
        print(
            f"{diagram.loop_count} loops: quadratic numerator agrees exactly with weighted LTD at three rational points",
            flush=True,
        )
    for numerator in (
        1 / energy,
        S("cff_test_sin")(energy),
        S("gammalooprs::K")(0, CIND(0)),
    ):
        try:
            triangle.integrate_energy(method="cff", numerator=numerator)
        except ValueError:
            pass
        else:
            raise AssertionError("accepted non-polynomial energy dependence")
    try:
        triangle.integrate_energy(
            method="cff", numerator=energy**2, contracted_edges=[edge]
        )
    except ValueError:
        pass
    else:
        raise AssertionError("silently ignored numerator constraints")
    print("unsupported inputs rejected explicitly", flush=True)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
