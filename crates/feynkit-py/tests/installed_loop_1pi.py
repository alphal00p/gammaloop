"""Explicit bridge policies preserve mixed tree/loop diagram sets."""

import inspect
from pathlib import Path

from symbolica.community import hepkit as hep
from symbolica.community.graph import Bridges

model = hep.Model(Path(__file__).parent / "fixtures/scalars_2p_3p.json")
states = [1000, 1000]


def identities(result: hep.BuiltDiagramSet, loops: int) -> set[str]:
    return {diagram.id for diagram in result if diagram.loop_count == loops}


for merged_restrictions in (True, False):
    restriction = hep.Restrict().only_vertices(["V_3_SCALAR_000"])
    if merged_restrictions:
        restriction = hep.Restrict().merge(restriction)

    def generate(
        loops: tuple[int, int] = (0, 1),
        bridges: Bridges = Bridges.one_pi_in_loops(),
    ) -> hep.BuiltDiagramSet:
        return (
            model.process(states, states)
            .amplitude(
                loops=hep.Loops.total(loops),
                restrict=restriction,
                topology=hep.TopologyRules(max_vertices=4, bridges=bridges),
            )
            .generate(run=hep.Run(threads=1))
            .build(hep.NumeratorOptions(naming=hep.Naming.named("D")))
        )

    loop_1pi = generate()
    all_graphs = generate(bridges=Bridges.any())
    strict = generate(bridges=Bridges.at_most(0))
    assert len(identities(loop_1pi, 0)) == 3
    assert identities(loop_1pi, 0) == identities(all_graphs, 0)
    assert not identities(strict, 0)
    assert identities(loop_1pi, 1)
    assert identities(loop_1pi, 1) == identities(strict, 1)
    assert len(identities(all_graphs, 1)) > len(identities(loop_1pi, 1))
    assert identities(generate(loops=(1, 1)), 1) == identities(loop_1pi, 1)
    assert identities(generate(bridges=Bridges.one_pi_in_loops()), 0) == identities(
        loop_1pi, 0
    )
    assert len(generate(loops=(0, 0), bridges=Bridges.at_most(1))) == 3
    minimal = (
        model.process(states, states)
        .amplitude(
            loops=hep.Loops.total((0, 1)),
            restrict=restriction,
            topology=hep.TopologyRules(max_vertices=4),
        )
        .generate(run=hep.Run(threads=1))
        .build()
    )
    assert {d.id for d in minimal} == {d.id for d in all_graphs}

for method in (hep.Process.amplitude, hep.Process.cross_section):
    assert (
        inspect.signature(method).parameters["loops"].default is inspect.Parameter.empty
    )
    assert "maximum_bridges" not in inspect.signature(method).parameters

standard_model = hep.Model.standard_model()
diphoton = (
    standard_model.process(["e-", "e+"], ["a", "a"])
    .amplitude(loops=hep.Loops.total(0))
    .generate()
    .build()
)
assert len(diphoton) == 2
assert all(diagram.loop_count == 0 for diagram in diphoton)
print(
    "Bridge policies: immutable restrictions, mixed loop ranges, explicit limits and minimal defaults passed"
)
