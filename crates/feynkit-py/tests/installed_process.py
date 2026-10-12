"""Neutral external states create explicit staged physics requests."""

import inspect
import subprocess
import sys
from functools import partial
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

from symbolica import E, S
from symbolica.community import hepkit as hep
from symbolica.community.graph import Bridges

model = hep.Model.standard_model()
process = model.process(["e-", "e+"], ["a", "a"])
assert process.model.name == model.name
assert not hasattr(process, "loop_count")
assert "Process(" in repr(process)
for old in (
    "Generator",
    "GenerationType",
    "GenerationSettings",
    "GenerationOptions",
    "GenerationResult",
):
    assert not hasattr(hep, old)
for old in (
    "generate_diagrams",
    "generate_amplitude",
    "generate_cross_section",
    "with_filters",
):
    assert not hasattr(hep.Process, old)
assert not hasattr(hep.Model, "generate_diagrams")
assert not hasattr(process, "generation_type")

request = process.amplitude(loops=hep.Loops.total(0))
named = hep.NumeratorOptions(naming=hep.Naming.named("D"))
compact = request.generate()
assert len(compact) == request.count().diagrams == 2
result = compact.build(named)
assert len(result) == 2
amplitude = result.into_amplitude()
assert isinstance(amplitude, hep.Amplitude)
assert {d.id for d in amplitude.diagrams} == {d.id for d in result}
assert amplitude.expression() == hep.Amplitude(result.diagrams).expression()
assert amplitude.expression().rank == 4
assert isinstance(amplitude.conjugate(), hep.Amplitude)
assert isinstance(amplitude.squared(), hep.SquaredAmplitude)

dimensional = result.into_amplitude(
    hep.AmplitudeOptions(dimension=S("D"), real=[S("real_coupling")])
)
assert (
    dimensional.expression()
    == hep.Amplitude(
        result.diagrams, dimension=S("D"), real=[S("real_coupling")]
    ).expression()
)

# Make the former cross-section cut and cyclic restrictions explicit.
xs_topology = hep.TopologyRules(
    max_vertices=4,
    self_energies=hep.Veto.all(),
    tadpoles=hep.Veto.all(),
    snails=hep.Veto.massless(),
)
xs_cuts = hep.Cuts(blobs=hep.Range.exact(1), spectators=hep.Range.exact(0))
xs_request = process.cross_section(
    loops=hep.Loops.forward(1), topology=xs_topology, cuts=xs_cuts
)
cross_section = xs_request.generate().build(named)
assert isinstance(cross_section, hep.BuiltDiagramSet)
assert cross_section.report.completed and len(cross_section) > 0
assert all(diagram.loop_count == 1 and diagram.cuts for diagram in cross_section)
assert not hasattr(process, "loop_count")
assert not hasattr(process, "symmetrizes_final")
assert {d.id for d in request.generate().build(named)} == {d.id for d in result}

for method in (process.amplitude, process.cross_section):
    signature = inspect.signature(method)
    assert signature.parameters["loops"].default is inspect.Parameter.empty
    assert (
        not {
            "particle_veto",
            "vertex_allow",
            "vertex_veto",
            "maximum_bridges",
            "progress",
            "filter_zero_color",
        }
        & signature.parameters.keys()
    )
    assert signature.parameters["topology"].default is None
    try:
        method()
    except TypeError:
        pass
    else:
        raise AssertionError("physics request accepted implicit loops")
assert inspect.signature(hep.Run).parameters["progress"].default is None
assert repr(hep.TopologyRules().bridges) == repr(Bridges.any())
assert hep.Cuts().blobs.maximum is None and hep.Cuts().spectators.maximum is None
assert (
    inspect.signature(hep.NumeratorOptions).parameters["filter_zero_color"].default
    is False
)
try:
    hep.Grouping("zeroes")
except ValueError:
    pass
else:
    raise AssertionError("zero filtering is still exposed as a grouping mode")

# Color filtering is independently opt-in at the explicit construction boundary.
for discovery, expected in (
    (request, result.diagrams),
    (request, amplitude.diagrams),
    (xs_request, cross_section.diagrams),
):
    updates = []
    filtered = discovery.generate().build(
        hep.NumeratorOptions(naming=hep.Naming.named("D"), filter_zero_color=True),
        run=hep.Run(progress=updates.append),
    )
    assert [d.to_json() for d in filtered] == [d.to_json() for d in expected]
    assert any(p.stage == hep.Stage.ZeroColor for p in updates)
    assert all(p.zero_numerator_count == 0 for p in updates)
    assert updates[-1].stage == hep.Stage.Done

updates = []
unfiltered_zeroes = compact.build(
    hep.NumeratorOptions(prefactor=E("0"), naming=hep.Naming.named("D"))
)
assert len(unfiltered_zeroes) == len(result)
filtered_zeroes = compact.build(
    hep.NumeratorOptions(
        prefactor=E("0"), naming=hep.Naming.named("D"), filter_zero_color=True
    ),
    run=hep.Run(progress=updates.append),
)
assert len(filtered_zeroes) == 0
assert filtered_zeroes.report.zero_numerator_count == len(unfiltered_zeroes)
assert updates[-1].zero_numerator_count == filtered_zeroes.report.zero_numerator_count
filter_updates = [p for p in updates if p.stage == hep.Stage.ZeroColor]
assert filter_updates and filter_updates[-1].completed == len(unfiltered_zeroes)
assert filter_updates[-1].zero_numerator_count == len(unfiltered_zeroes)
assert all(
    before.zero_numerator_count <= after.zero_numerator_count
    for before, after in zip(filter_updates, filter_updates[1:])
)

# Notebook progress displays the live zero count and closes every indicator.
context = MagicMock()
indicator = context.__enter__.return_value
marimo = SimpleNamespace(
    running_in_notebook=lambda: True,
    status=SimpleNamespace(
        spinner=MagicMock(return_value=context),
        progress_bar=MagicMock(return_value=context),
    ),
    output=SimpleNamespace(_output=SimpleNamespace(flush=MagicMock())),
)
with patch.dict(sys.modules, {"marimo": marimo}):
    compact.build(
        hep.NumeratorOptions(prefactor=E("0"), filter_zero_color=True),
        run=hep.Run(progress="auto"),
    )
filter_subtitles = [
    call.args[-1]
    for call in indicator.update.call_args_list
    if call.args[-2] == "Filtering zero-color diagrams"
]
assert filter_subtitles
assert f"{len(unfiltered_zeroes)} zero graphs filtered" in filter_subtitles[-1]
assert context.__exit__.call_count == context.__enter__.call_count
context.__exit__.assert_called_with(None, None, None)

cancel = hep.CancellationToken()
cancel.cancel()
incomplete = request.generate(run=hep.Run(cancellation=cancel))
assert not incomplete.report.completed
try:
    incomplete.build().into_amplitude()
except hep.GenerationError as error:
    assert "cancelled" in str(error)
else:
    raise AssertionError("cancelled generation produced an amplitude")
try:
    process.amplitude(
        loops=hep.Loops.total(0), topology=hep.TopologyRules(bridges=Bridges.at_most(0))
    ).generate().build().into_amplitude()
except hep.AmplitudeError as error:
    assert "empty" in str(error).lower() or "one diagram" in str(error).lower()
else:
    raise AssertionError("empty generation produced an amplitude")

scalar_model = hep.Model(Path(__file__).parent / "fixtures/scalars_2p_3p.json")
scalar = scalar_model.process([1000], [1000, 1000])
inclusive = scalar.with_final_state_alternatives([[1000], [1000, 1000]])
assert scalar.outgoing_alternatives != inclusive.outgoing_alternatives
try:
    inclusive.amplitude(loops=hep.Loops.total(0))
except hep.GenerationError:
    pass
else:
    raise AssertionError("an amplitude accepted multiple external states")
updates: list[hep.GenerationProgress] = []
scalar_request = scalar.amplitude(
    loops=hep.Loops.total(0), restrict=hep.Restrict().only_vertices(["V_3_SCALAR_000"])
)
generate = partial(
    scalar_request.generate, run=hep.Run(threads=1, progress=updates.append)
)
assert isinstance(generate().build().into_amplitude(), hep.Amplitude)
assert updates[-1].stage == hep.Stage.Done
print(
    "Process: staged diagrams, amplitudes, physical cuts, dimensions, cancellation, defaults and callbacks passed"
)

# Restrictions are immutable request criteria, reused across discovery operations.
qed_vertex = model.vertex_rule("V_98")
restricted_process = model.process([model.particle("e-"), "e+"], ["a", "a"])
restriction = hep.Restrict().only_vertices([qed_vertex])
original_restriction = repr(restriction)
restricted = restricted_process.amplitude(
    loops=hep.Loops.total(0), restrict=restriction
)
assert {d.id for d in restricted.generate().build(named)} == {d.id for d in result}
assert {d.id for d in restricted.generate().build(named).into_amplitude().diagrams} == {
    d.id for d in result
}
assert {
    d.id
    for d in restricted_process.cross_section(
        loops=hep.Loops.forward(1),
        restrict=restriction,
        topology=xs_topology,
        cuts=xs_cuts,
    )
    .generate()
    .build(named)
} == {d.id for d in cross_section}
vetoed = restriction.exclude_particles(["e-"]).exclude_vertices([qed_vertex])
assert repr(restriction) == original_restriction
for blocked in (
    restricted_process.amplitude(loops=hep.Loops.total(1), restrict=vetoed),
    restricted_process.cross_section(
        loops=hep.Loops.forward(1), restrict=vetoed, cuts=xs_cuts
    ),
):
    assert len(blocked.generate()) == 0
try:
    restricted_process.amplitude(
        loops=hep.Loops.total(0), restrict=vetoed
    ).generate().build().into_amplitude()
except hep.AmplitudeError:
    pass
else:
    raise AssertionError("amplitude construction bypassed request restrictions")
assert (
    len(process.amplitude(loops=hep.Loops.total(0), restrict=hep.Restrict()).generate())
    == 2
)
assert (
    len(
        process.amplitude(
            loops=hep.Loops.total(0), restrict=hep.Restrict().only_vertices([])
        ).generate()
    )
    == 0
)
assert {
    d.id
    for d in process.amplitude(
        loops=hep.Loops.total(0), restrict=hep.Restrict().merge(restriction)
    )
    .generate()
    .build()
} == {d.id for d in result}
for name in ("restrict", "topology", "loops"):
    try:
        setattr(restricted, name, None)
    except AttributeError:
        pass
    else:
        raise AssertionError("request criteria are mutable")

# Bad references fail before any generation work begins.
for make in (
    lambda: model.process(["missing_particle"], []),
    lambda: process.amplitude(
        loops=hep.Loops.total(0),
        restrict=hep.Restrict().exclude_particles(["missing_particle"]),
    ),
    lambda: process.amplitude(
        loops=hep.Loops.total(0),
        restrict=hep.Restrict().only_vertices(["missing_vertex"]),
    ),
    lambda: process.amplitude(
        loops=hep.Loops.total(0),
        restrict=hep.Restrict().exclude_vertices(["missing_vertex"]),
    ),
    lambda: process.with_final_state_alternatives([["missing_particle"]]),
):
    try:
        make()
    except hep.GenerationError:
        pass
    else:
        raise AssertionError("request accepted an unknown model member")
try:
    hep.Process(model, [], [])
except TypeError:
    pass
else:
    raise AssertionError("Process must be constructed from Model.process")
print(
    "Model.process: validation, immutable restrictions, unrestricted defaults and shared criteria passed"
)

# Run the shared native/installed contract in its own model session.
subprocess.run(
    [sys.executable, str(Path(__file__).with_name("request_decorator_contract.py"))],
    check=True,
)
