"""A neutral process owns all public generation operations."""

import inspect
import sys
from functools import partial
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

from symbolica import E, S
from symbolica.community import hepkit as hep

model = hep.Model.standard_model()
process = model.process(["e-", "e+"], ["a", "a"])
assert process.model.name == model.name
assert not hasattr(process, "loop_count")
assert "Process(" in repr(process)
assert not hasattr(hep, "Generator")
assert not hasattr(hep, "GenerationType")
assert not hasattr(hep.Model, "generate_diagrams")
assert not hasattr(process, "generation_type")
assert not hasattr(hep.Process, "amplitude")
assert not hasattr(hep.Process, "cross_section")

result = process.generate_diagrams(progress=None)
assert len(result) == 2
amplitude = process.generate_amplitude(progress=None)
assert isinstance(amplitude, hep.Amplitude)
assert {d.id for d in amplitude.diagrams} == {d.id for d in result}
assert amplitude.expression() == hep.Amplitude(result.diagrams).expression()
assert amplitude.expression().rank == 4
assert isinstance(amplitude.conjugate(), hep.Amplitude)
assert isinstance(amplitude.squared(), hep.SquaredAmplitude)

dimensional = process.generate_amplitude(
    dimension=S("D"), real=[S("real_coupling")], progress=None
)
assert (
    dimensional.expression()
    == hep.Amplitude(
        result.diagrams, dimension=S("D"), real=[S("real_coupling")]
    ).expression()
)

cross_section = process.generate_cross_section(loops=1, max_vertices=4, progress=None)
assert isinstance(cross_section, hep.GenerationResult)
assert cross_section.report.completed and len(cross_section) > 0
assert all(diagram.loop_count == 1 and diagram.cuts for diagram in cross_section)
assert not hasattr(process, "loop_count")
assert not hasattr(process, "symmetrizes_final")
assert {d.id for d in process.generate_diagrams(progress=None)} == {
    d.id for d in result
}

for method in (
    process.generate_diagrams,
    process.generate_amplitude,
    process.generate_cross_section,
):
    signature = inspect.signature(method)
    assert signature.parameters["loops"].default == 0
    assert signature.parameters["symmetrize_final"].default is None
    assert (
        not {"particle_veto", "vertex_allow", "vertex_veto"}
        & signature.parameters.keys()
    )
    assert signature.parameters["maximum_bridges"].default is Ellipsis
    assert signature.parameters["progress"].default == "auto"
    assert signature.parameters["filter_zero_color"].default is False

try:
    hep.NumeratorGrouping("zeroes")
except ValueError:
    pass
else:
    raise AssertionError("zero filtering is still exposed as a grouping mode")

# Color filtering is independently opt-in on all generation entry points.
for method, kwargs, expected in (
    (process.generate_diagrams, {}, result.diagrams),
    (process.generate_amplitude, {}, amplitude.diagrams),
    (
        process.generate_cross_section,
        {"loops": 1, "max_vertices": 4},
        cross_section.diagrams,
    ),
):
    updates = []
    filtered = method(filter_zero_color=True, progress=updates.append, **kwargs)
    assert [d.to_json() for d in filtered.diagrams] == [d.to_json() for d in expected]
    assert any(p.stage == "filter_zero_color" for p in updates)
    assert all(p.zero_numerator_count == 0 for p in updates)
    assert updates[-1].stage == "complete"

updates = []
unfiltered_zeroes = process.generate_diagrams(numerator_prefactor=E("0"), progress=None)
assert len(unfiltered_zeroes) == len(result)
filtered_zeroes = process.generate_diagrams(
    numerator_prefactor=E("0"), filter_zero_color=True, progress=updates.append
)
assert len(filtered_zeroes) == 0
assert filtered_zeroes.report.zero_numerator_count == len(unfiltered_zeroes)
assert updates[-1].zero_numerator_count == filtered_zeroes.report.zero_numerator_count
filter_updates = [p for p in updates if p.stage == "filter_zero_color"]
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
    process.generate_diagrams(numerator_prefactor=E("0"), filter_zero_color=True)
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
assert not process.generate_diagrams(
    cancellation_token=cancel, progress=None
).report.completed
try:
    process.generate_amplitude(cancellation_token=cancel, progress=None)
except hep.GenerationError as error:
    assert "cancelled" in str(error)
else:
    raise AssertionError("cancelled generation produced an amplitude")
try:
    process.generate_amplitude(maximum_bridges=0, progress=None)
except hep.AmplitudeError as error:
    assert "empty" in str(error).lower() or "one diagram" in str(error).lower()
else:
    raise AssertionError("empty generation produced an amplitude")

scalar_model = hep.Model(Path(__file__).parent / "fixtures/scalars_2p_3p.json")
scalar = scalar_model.process([1000], [1000, 1000])
inclusive = scalar.with_final_state_alternatives([[1000], [1000, 1000]])
assert scalar.outgoing_alternatives != inclusive.outgoing_alternatives
for method in (inclusive.generate_diagrams, inclusive.generate_amplitude):
    try:
        method(progress=None)
    except hep.GenerationError:
        pass
    else:
        raise AssertionError("an amplitude accepted multiple external states")

updates: list[hep.GenerationProgress] = []
generate = partial(
    scalar.with_filters(vertex_allow=["V_3_SCALAR_000"]).generate_amplitude, threads=1
)
assert isinstance(generate(progress=updates.append), hep.Amplitude)
assert updates[-1].stage == "complete"
print(
    "Process: diagrams, amplitudes, physical cuts, dimensions, cancellation, defaults and callbacks passed"
)

# Restrictions are immutable process state, shared by every generation operation.
qed_vertex = model.vertex_rule("V_98")
restricted = model.process(
    [model.particle("e-"), "e+"], ["a", "a"], vertex_allow=[qed_vertex]
)
assert restricted.vertex_allow is not None
assert [v.name for v in restricted.vertex_allow] == ["V_98"]
assert restricted.particle_veto == [] and restricted.vertex_veto == []
assert process.vertex_allow is None
assert "vertex_allow=[" in repr(restricted)
assert {d.id for d in restricted.generate_diagrams(progress=None)} == {
    d.id for d in result
}
assert {d.id for d in restricted.generate_amplitude(progress=None).diagrams} == {
    d.id for d in result
}
assert {d.id for d in restricted.generate_cross_section(loops=1, progress=None)} == {
    d.id for d in cross_section
}

vetoed = restricted.with_filters(particle_veto=["e-"], vertex_veto=[qed_vertex])
assert vetoed.vertex_allow is not None and vetoed.vertex_allow[0].name == "V_98"
assert vetoed.particle_veto == ["e-"]
assert [v.name for v in vetoed.vertex_veto] == ["V_98"]
assert restricted.particle_veto == [] and restricted.vertex_veto == []
for method in (vetoed.generate_diagrams, vetoed.generate_cross_section):
    assert len(method(loops=1, progress=None)) == 0
try:
    vetoed.generate_amplitude(progress=None)
except hep.AmplitudeError:
    pass
else:
    raise AssertionError("amplitude generation bypassed process restrictions")
cleared = vetoed.with_filters(particle_veto=None, vertex_veto=None, vertex_allow=None)
assert (
    cleared.vertex_allow is None
    and cleared.vertex_veto == []
    and cleared.particle_veto == []
)
assert len(cleared.generate_diagrams(progress=None)) == 2
assert len(process.with_filters(vertex_allow=[]).generate_diagrams(progress=None)) == 0
assert process.with_filters().vertex_allow is None
preserved = restricted.with_filters(vertex_allow=...).vertex_allow
assert preserved is not None and preserved[0].name == "V_98"
for name in ("particle_veto", "vertex_allow", "vertex_veto"):
    try:
        setattr(restricted, name, [])
    except AttributeError:
        pass
    else:
        raise AssertionError("process restrictions are mutable")

# Bad references fail before any generation work begins.
for make in (
    lambda: model.process(["missing_particle"], []),
    lambda: model.process([], [], particle_veto=["missing_particle"]),
    lambda: model.process([], [], vertex_allow=["missing_vertex"]),
    lambda: process.with_filters(vertex_veto=["missing_vertex"]),
    lambda: process.with_final_state_alternatives([["missing_particle"]]),
):
    try:
        make()
    except hep.GenerationError:
        pass
    else:
        raise AssertionError("process accepted an unknown model member")
try:
    getattr(hep, "Process")(model, [], [])
except TypeError:
    pass
else:
    raise AssertionError("Process must be constructed from Model.process")
print(
    "Model.process: validation, immutable restrictions, clearing, and shared generation filters passed"
)
