"""Validate scalar tensor returns and export scoped native notebook displays."""

import sys
from pathlib import Path

import marimo as mo
from symbolica import Expression
from symbolica.community import hepkit as hep
from symbolica.community.tensor import TensorExpression


def check(directory):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    model = hep.Model.phi3()
    diagram = hep.FeynmanDiagram.from_dot(
        model,
        """digraph double_box { edge [particle="phi"]; ext [style=invis];
        ext -> a; ext -> d; c -> ext; f -> ext;
        a -> b; b -> c; d -> e; e -> f;
        a -> d; b -> e [lmb_id=0]; c -> f [lmb_id=1]; }""",
    )
    cff = diagram.integrate_energy(method="cff")
    orientation = next(o for o in reversed(cff.orientations) if len(o.families) > 1)
    three_loop = hep.FeynmanDiagram.from_dot(
        model, (Path(__file__).parent / "fixtures/ltd_three_loop.dot").read_text()
    )
    family = max(
        three_loop.integrate_energy(method="cff").orientations,
        key=lambda o: len(o.families),
    ).families[-1]
    ltd = three_loop.integrate_energy(method="ltd")
    residue = ltd.residues[-1]
    pair = next(iter(residue.surface_pairs.values()))
    surface = cff.surfaces[0]
    energy = next(iter(cff.on_shell_energies.values()))
    for value in (cff, ltd, orientation, family, residue, surface, energy, pair.minus):
        expression = value.to_expression()
        assert type(expression) is TensorExpression
        assert isinstance(expression, Expression) and expression.rank == 0
        assert type(expression.to_expression()) is Expression
    edge = next(iter(residue.surface_pairs))
    assert (
        (
            pair.minus.to_expression()
            - residue.energy_map[edge]
            + ltd.on_shell_energies[edge].symbol
        )
        .expand()
        .is_zero()
    )

    for name, value in {
        "orientation": orientation,
        "family": family,
        "residue": residue,
    }.items():
        assert "object at" not in repr(value)
        source = value._repr_html_()
        (directory / f"{name}.html").write_text(source)
        assert "iframe" in mo.as_html(value).text

    records = {
        "surface": surface,
        "factor": pair.minus,
        "pair": pair,
        "energy": energy,
        "group": cff.raised_surface_groups()[0],
        "cff-report": cff.report,
        "ltd-report": ltd.report,
    }
    for name, value in records.items():
        assert "object at" not in repr(value)
        source = value._repr_html_()
        assert "<table" in source
        if "report" not in name:
            assert "data-spenso-math" in source, name
        assert "<table" in mo.as_html(value).text
        (directory / f"{name}.html").write_text(source)
    coefficient = cff.pole_coefficients(cff.raised_surface_groups()[0])[0]
    assert "pole coefficient, order 1" in repr(coefficient)
    assert '"pole_order":1' in coefficient._repr_html_()
    print(
        "Scoped native displays, scalar tensor returns, reports and pole-coefficient provenance passed"
    )


if __name__ == "__main__":
    check(sys.argv[1])
