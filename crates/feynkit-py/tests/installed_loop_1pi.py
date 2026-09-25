"""Default bridge filtering keeps trees and selects 1PI loop diagrams per graph."""

import inspect
from pathlib import Path
from types import EllipsisType

from symbolica.community import hep

model = hep.Model(Path(__file__).parent / "fixtures/scalars_2p_3p.json")
states = [1000, 1000]


def identities(result: hep.GenerationResult, loops: int) -> set[str]:
    return {diagram.id for diagram in result.diagrams if diagram.loop_count == loops}


for factory_filters in (True, False):

    def generate(
        loops: tuple[int, int] = (0, 1),
        maximum_bridges: int | None | EllipsisType = ...,
    ) -> hep.GenerationResult:
        if factory_filters:
            process = model.process(states, states, vertex_allow=["V_3_SCALAR_000"])
        else:
            process = model.process(states, states).with_filters(
                vertex_allow=["V_3_SCALAR_000"]
            )
        method = process.generate_diagrams
        return method(
            max_vertices=4,
            threads=1,
            loops=loops,
            self_energy=None,
            tadpoles=None,
            zero_snails=None,
            progress=None,
            maximum_bridges=maximum_bridges,
        )

    automatic = generate()
    all_graphs = generate(maximum_bridges=None)
    strict = generate(maximum_bridges=0)
    assert len(identities(automatic, 0)) == 3
    assert identities(automatic, 0) == identities(all_graphs, 0)
    assert not identities(strict, 0)
    assert identities(automatic, 1)
    assert identities(automatic, 1) == identities(strict, 1)
    assert len(identities(all_graphs, 1)) > len(identities(automatic, 1))
    assert identities(generate(loops=(1, 1)), 1) == identities(automatic, 1)
    assert identities(generate(maximum_bridges=...), 0) == identities(automatic, 0)
    assert len(generate(loops=(0, 0), maximum_bridges=1).diagrams) == 3

for method in (
    hep.Process.generate_diagrams,
    hep.Process.generate_amplitude,
    hep.Process.generate_cross_section,
):
    assert inspect.signature(method).parameters["maximum_bridges"].default is Ellipsis

standard_model = hep.Model.standard_model()
diphoton = standard_model.process(["e-", "e+"], ["a", "a"]).generate_diagrams(
    progress=None
)
assert len(diphoton.diagrams) == 2
assert all(diagram.loop_count == 0 for diagram in diphoton.diagrams)
print(
    "Loop-only 1PI default: factory and copied restrictions, per-call loops, mixed ranges, explicit limits and diphoton passed"
)
