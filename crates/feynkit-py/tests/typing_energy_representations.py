"""Type-check against the installed community stubs with pyright."""

from typing import Literal, assert_type

from symbolica import E, Expression
from symbolica.community import hepkit as hep
from symbolica.community.tensor import TensorExpression


def check_return_types(
    diagram: hep.FeynmanDiagram, method: Literal["cff", "ltd"]
) -> None:
    cff = diagram.integrate_energy(method="cff", max_orientations=1000)
    ltd = diagram.integrate_energy(method="ltd")
    bounded = diagram.integrate_energy(method="cff", numerator=E("1"))
    assert_type(bounded, hep.CffRepresentation)
    assert_type(bounded.energy_degree_bounds, dict[int, int])
    assert_type(bounded.orientations[0].families[0].energy_map, dict[int, Expression])
    assert_type(cff, hep.CffRepresentation)
    assert_type(ltd, hep.LtdRepresentation)
    assert_type(cff.orientations[0].families[0], hep.CrossFreeFamily)
    assert_type(cff.orientations[0].edge_signs, dict[int, int | None])
    assert_type(cff.orientations[0].to_expression(), TensorExpression)
    assert_type(cff.orientations[0].families[0].to_expression(), TensorExpression)
    assert_type(cff.surfaces[0].to_expression(), TensorExpression)
    assert_type(cff.on_shell_energies[0].to_expression(), TensorExpression)
    assert_type(ltd.residues[0].to_expression(), TensorExpression)
    assert_type(
        ltd.residues[0].surface_pairs[0].minus.to_expression(), TensorExpression
    )
    assert_type(ltd.residues[0], hep.LtdResidue)
    assert_type(cff.surfaces[0], hep.EnergySurface)
    assert_type(ltd.surfaces[0], hep.EnergySurface)
    assert_type(cff.on_shell_energies[0], hep.OnShellEnergy)
    assert_type(ltd.residues[0].surface_pairs[0], hep.SurfacePair)
    result = diagram.integrate_energy(method=method)
    assert_type(result, hep.CffRepresentation | hep.LtdRepresentation)
    assert_type(result.to_expression(expand_surfaces=True), TensorExpression)
