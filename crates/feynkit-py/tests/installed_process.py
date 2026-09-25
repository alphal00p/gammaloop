"""A neutral process owns all public generation operations."""

import inspect
from functools import partial
from pathlib import Path

from symbolica import S
from symbolica.community import hep

model = hep.Model.standard_model()
process = hep.Process(model, ["e-", "e+"], ["a", "a"])
assert process.model.name == model.name
assert process.loop_count == (0, 0)
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
assert process.loop_count == (0, 0)
assert process.symmetrizes_final is None
assert {d.id for d in process.generate_diagrams(progress=None)} == {
    d.id for d in result
}

for method in (
    process.generate_diagrams,
    process.generate_amplitude,
    process.generate_cross_section,
):
    signature = inspect.signature(method)
    assert signature.parameters["maximum_bridges"].default is Ellipsis
    assert signature.parameters["progress"].default == "auto"

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
scalar = hep.Process(scalar_model, [1000], [1000, 1000])
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
    scalar.generate_amplitude, vertex_allow=["V_3_SCALAR_000"], threads=1
)
assert isinstance(generate(progress=updates.append), hep.Amplitude)
assert updates[-1].stage == "complete"
print(
    "Process: diagrams, amplitudes, physical cuts, dimensions, cancellation, defaults and callbacks passed"
)
