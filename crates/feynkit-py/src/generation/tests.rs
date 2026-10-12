use std::ffi::CString;

use pyo3::types::PyDict;

use super::*;
use linnet_py::PyGraph;

#[test]
fn group_member_wrapper_fits_wasm_allocator_alignment() {
    assert!(std::mem::align_of::<PyGroupMember>() <= 8);
    assert!(
        std::mem::align_of::<<PyGroupMember as pyo3::impl_::pyclass::PyClassImpl>::Layout>() <= 8
    );
}

#[test]
fn automatic_progress_tracks_marimo_lifecycle() {
    Python::initialize();
    Python::attach(|py| {
        let module = PyModule::new(py, "symbolica.community.hepkit").unwrap();
        crate::initialize_feynkit(&module).unwrap();
        let locals = PyDict::new(py);
        locals.set_item("fk", &module).unwrap();
        let graph_module = PyModule::new(py, "symbolica.community.graph").unwrap();
        linnet_py::register(&graph_module).unwrap();
        locals.set_item("graph", &graph_module).unwrap();
        locals
            .set_item(
                "MODEL_JSON",
                include_str!("../../tests/fixtures/scalars_2p_3p.json"),
            )
            .unwrap();
        let code = CString::new(
                r#"
import inspect
import sys
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

model = fk.Model.from_json(MODEL_JSON)
process = model.process([1000], [1000, 1000])
restriction = fk.Restrict().only_vertices(["V_3_SCALAR_000"])
topology = fk.TopologyRules(max_vertices=3, bridges=graph.Bridges.one_pi_in_loops(), self_loops=True,
    self_energies=fk.Veto.all(), tadpoles=fk.Veto.all(), snails=fk.Veto.massless())
context = MagicMock()
indicator = context.__enter__.return_value
marimo = SimpleNamespace(
    running_in_notebook=MagicMock(return_value=True),
    status=SimpleNamespace(
        spinner=MagicMock(return_value=context),
        progress_bar=MagicMock(return_value=context),
    ),
    output=SimpleNamespace(_output=SimpleNamespace(flush=MagicMock())),
)

for via_model in (True, False):
    def generate(progress='auto', cancellation_token=None, **hooks):
        current = model.process([1000], [1000, 1000]) if via_model else process
        request = current.amplitude(loops=fk.Loops.total(1), restrict=restriction, topology=topology, **hooks)
        run = fk.Run(threads=2, progress=progress, cancellation=cancellation_token)
        return request.generate(run=run).build(run=run).group(fk.Grouping("up_to_scalar"), run=run)

    assert inspect.signature(fk.Run).parameters["progress"].default is None
    with patch.dict(sys.modules, {"marimo": None}):
        baseline = generate()
    with patch.dict(sys.modules, {"marimo": marimo}):
        context.reset_mock()
        marimo.status.spinner.reset_mock()
        marimo.status.progress_bar.reset_mock()
        marimo.running_in_notebook.return_value = False
        assert len(generate()) == len(baseline)
        marimo.status.spinner.assert_not_called()
        marimo.status.progress_bar.assert_not_called()
        marimo.running_in_notebook.return_value = True
        assert len(generate(progress=None)) == len(baseline)
        updates = []
        generate(progress=updates.append)
        assert updates[-1].stage == fk.Stage.Done
        marimo.status.spinner.assert_not_called()
        marimo.status.progress_bar.assert_not_called()

        visible = []
        marimo.output._output.flush.side_effect = lambda: visible.append(indicator.update.call_args.args[-2])
        result = generate()
        marimo.status.spinner.assert_called()
        marimo.status.progress_bar.assert_called()
        for call in marimo.status.progress_bar.call_args_list:
            assert call.kwargs["total"] > 0
            assert call.kwargs["remove_on_exit"] is True
            assert call.kwargs["show_eta"] is False
            assert call.kwargs["show_rate"] is False
        context.__enter__.assert_called()
        context.__exit__.assert_called_with(None, None, None)
        assert context.__exit__.call_count == context.__enter__.call_count
        titles = list(dict.fromkeys(call.args[-2] for call in indicator.update.call_args_list))
        assert titles == ["Enumerating skeletons", "Assigning interactions", "Preparing physical diagrams",
            "Diagram generation complete", "Constructing numerators", "Preparing numerators for grouping",
            "Sampling numerators for grouping", "Comparing numerators", "Finalizing diagram groups"]
        # Every update must be flushed, even rapid transitions that Marimo's
        # own 150 ms throttle would suppress until another update arrives.
        assert visible == [call.args[-2] for call in indicator.update.call_args_list]
        assert indicator.update.call_args.args[-1].startswith(f"{len(result)} diagrams retained")
        assert all("elapsed" in call.args[-1] for call in indicator.update.call_args_list)
        assert any(" / " in call.args[-1] for call in indicator.update.call_args_list)

        context.reset_mock()
        empty = generate(progress="auto", prune_skeletons=lambda view: False)
        assert len(empty) == 0
        assert indicator.update.call_args.args[-1].startswith("0 diagrams retained")
        context.__exit__.assert_called_with(None, None, None)
        assert context.__exit__.call_count == context.__enter__.call_count

        context.reset_mock()
        token = fk.CancellationToken()
        token.cancel()
        result = generate(cancellation_token=token)
        assert not result.report.completed
        assert indicator.update.call_args.args[-2] == "Diagram generation cancelled"
        context.__exit__.assert_called_with(None, None, None)
        assert context.__exit__.call_count == context.__enter__.call_count

        for exception in (RuntimeError("filter failed"), KeyboardInterrupt()):
            context.reset_mock()
            def fail(*args):
                raise exception
            try:
                generate(prune_dressings=fail)
            except BaseException as caught:
                assert caught is exception
            else:
                raise AssertionError("generation exception was lost")
            assert context.__exit__.call_count == context.__enter__.call_count
            assert context.__exit__.call_args.args[:2] == (type(exception), exception)
            assert indicator.update.call_args.args[-2] == "Diagram generation stopped"

        context.reset_mock()
        exception = RuntimeError("spinner failed")
        indicator.update.side_effect = exception
        try:
            generate()
        except RuntimeError as caught:
            assert caught is exception
        else:
            raise AssertionError("spinner exception was lost")
        assert context.__exit__.call_count == context.__enter__.call_count
        indicator.update.side_effect = None

        for invalid in ("invalid", True, 1):
            context.reset_mock()
            try:
                generate(progress=invalid)
            except TypeError:
                pass
            else:
                raise AssertionError("invalid progress was accepted")
            context.__enter__.assert_not_called()
"#,
            )
            .unwrap();
        py.run(&code, Some(&locals), Some(&locals)).unwrap();
    });
}

#[test]
fn minimal_output_preserves_physics_and_explicit_id_demands() {
    Python::initialize();
    Python::attach(|py| {
        let module = PyModule::new(py, "symbolica.community.hepkit").unwrap();
        crate::initialize_feynkit(&module).unwrap();
        let locals = PyDict::new(py);
        locals.set_item("fk", &module).unwrap();
        let graph_module = PyModule::new(py, "symbolica.community.graph").unwrap();
        linnet_py::register(&graph_module).unwrap();
        locals.set_item("graph", &graph_module).unwrap();
        locals
            .set_item(
                "MODEL_JSON",
                include_str!("../../tests/fixtures/scalars_2p_3p.json"),
            )
            .unwrap();
        let code = CString::new(r#"
import inspect
import json
model = fk.Model.from_json(MODEL_JSON)
process = model.process([1000, 1000], [1000, 1000])
request = process.amplitude(loops=fk.Loops.total(0), restrict=fk.Restrict().only_vertices(["V_3_SCALAR_000"]))
run = fk.Run(threads=1)
compact = request.generate(run=run)
assert len(compact) == request.count(run=run).diagrams
assert all(isinstance(d, fk.PreparedDiagram) for d in compact)
naming = fk.NumeratorOptions(naming=fk.Naming.named("D"))
named = compact.build(naming, run=run)
minimal = compact.build(run=run)
assert len(named) == len(minimal) > 0
for built in (named, minimal):
    assert len(built.select(fk.Selection(select=[]))) == 0
    assert len(built.select(fk.Selection(veto=[]))) == len(named)
assert all(d.name is not None for d in named)
assert all(d.name is None for d in minimal)
for left, right in zip(named, minimal):
    left_json, right_json = json.loads(left.to_json()), json.loads(right.to_json())
    left_json.pop("name")
    right_json.pop("name")
    assert left_json == right_json
for left, right in zip(named.groups, minimal.groups):
    assert left.master == right.master
    assert len(left.members) == len(right.members)
    for a, b in zip(left.members, right.members):
        assert a.source_id is not None and a.source_name is not None
        assert b.source_id is None and b.source_name is None
        assert (a.source_diagram, a.diagram, a.ratio) == (b.source_diagram, b.diagram, b.ratio)
        assert str(a.overall_factor_expression()) == str(b.overall_factor_expression())
# A short hexadecimal name must remain a name, including keyed momentum bases.
named_choice = named[0]
assert [d.id for d in named.select(fk.Selection(select=[named_choice.name]))] == [named_choice.id]
assert {d.id for d in named.select(fk.Selection(veto=[named_choice.name]))} == {d.id for d in named} - {named_choice.id}
routed = compact.build(fk.NumeratorOptions(naming=fk.Naming.named("D"), loop_momentum_bases=[(named_choice.name, [])]), run=run)
assert [d.id for d in routed] == [d.id for d in named]
try:
    compact.build(fk.NumeratorOptions(naming=fk.Naming.named("D"), loop_momentum_bases=[(named_choice.name, [0])]), run=run)
except fk.GenerationError as error:
    assert "loop-momentum" in str(error)
else:
    raise AssertionError("named momentum-basis override was ignored")
chosen = minimal[0]
assert len(chosen.id) == 32
for selector in (chosen, chosen.id):
    selected = minimal.select(fk.Selection(select=[selector]))
    assert [d.id for d in selected] == [chosen.id] and selected[0].name is None
    assert chosen.id not in [d.id for d in minimal.select(fk.Selection(veto=[selector]))]
    routed = compact.build(fk.NumeratorOptions(loop_momentum_bases=[(selector, [])]), run=run)
    assert len(routed) == len(minimal)
for operation in (
    lambda: minimal.select(fk.Selection(select=[named_choice.name])),
    lambda: minimal.select(fk.Selection(veto=[named_choice.name])),
    lambda: compact.build(fk.NumeratorOptions(loop_momentum_bases=[(named_choice.name, [])])),
):
    try:
        operation()
    except fk.GenerationError as error:
        assert "named" in str(error)
    else:
        raise AssertionError("unnamed build accepted a name requirement")
assert not hasattr(fk.Process, "generate_diagrams")
assert not hasattr(fk.Process, "generate_amplitude")
assert not hasattr(fk.Process, "generate_cross_section")
assert inspect.signature(fk.Process.amplitude).parameters["loops"].default is inspect.Parameter.empty
assert inspect.signature(fk.Process.cross_section).parameters["loops"].default is inspect.Parameter.empty
# Iterators and prepared outputs own their native storage independently of the set wrapper.
iterator = iter(compact)
first = next(iterator)
skeleton = first.skeleton
before = [d.to_json() for d in compact.build()]
first.graph.add_node(graph.Node(data="independent copy"))
assert [d.to_json() for d in compact.build()] == before
assert skeleton.automorphisms > 0
assert len(first.rules) > 0 and all(isinstance(rule, fk.VertexRule) for rule in first.rules)
del compact
assert list(iterator) and first.aut > 0 and skeleton.loops == 0

"#).unwrap();
        py.run(&code, Some(&locals), Some(&locals)).unwrap();
    });
}

#[test]
fn compact_html_reports_discovery_without_demanding_numerators() {
    Python::initialize();
    Python::attach(|py| {
        let module = PyModule::new(py, "symbolica.community.hepkit").unwrap();
        crate::initialize_feynkit(&module).unwrap();
        let locals = PyDict::new(py);
        locals.set_item("fk", &module).unwrap();
        locals
            .set_item(
                "MODEL_JSON",
                include_str!("../../tests/fixtures/scalars_2p_3p.json"),
            )
            .unwrap();
        let code = CString::new(
            r#"
import json
definition = json.loads(MODEL_JSON)
# Discovery needs rules and orders, but this invalid vertex momentum must fail at build.
for lorentz in definition["lorentz_structures"]:
    if lorentz["name"] == "SCALAR_3_LORENTZ_STRUCTURE":
        lorentz["structure"] = "P(1,99)"
model = fk.Model.from_json(json.dumps(definition))
request = model.process([1000, 1000], [1000, 1000]).amplitude(
    loops=fk.Loops.total(0),
    restrict=fk.Restrict().only_vertices(["V_3_SCALAR_000"]),
)
compact = request.generate(run=fk.Run(threads=1))
assert len(compact) == 3 and compact.report.completed
html = compact._repr_html_()
assert "retained diagrams" in html and str(len(compact)) in html and "completed" in html
try:
    compact.build(run=fk.Run(threads=1))
except fk.GenerationError as error:
    assert "references UFO leg 99" in str(error)
else:
    raise AssertionError("the cold-summary control did not require numerator construction")
assert compact._repr_html_() == html and compact.report.completed
# A cancelled empty result must not be presented as a completed zero-diagram search.
token = fk.CancellationToken()
token.cancel()
partial = request.generate(run=fk.Run(threads=1, cancellation=token))
assert len(partial) == 0 and not partial.report.completed
assert "cancelled" in partial._repr_html_()
"#,
        )
        .unwrap();
        py.run(&code, Some(&locals), Some(&locals)).unwrap();
    });
}

#[test]
fn generation_defaults_leave_numerator_grouping_opt_in() {
    Python::initialize();
    Python::attach(|py| {
        let module = PyModule::new(py, "symbolica.community.hepkit").unwrap();
        crate::initialize_feynkit(&module).unwrap();
        let locals = PyDict::new(py);
        locals.set_item("fk", &module).unwrap();
        let graph_module = PyModule::new(py, "symbolica.community.graph").unwrap();
        linnet_py::register(&graph_module).unwrap();
        locals.set_item("graph", &graph_module).unwrap();
        locals
            .set_item(
                "MODEL_JSON",
                include_str!("../../tests/fixtures/scalars_2p_3p.json"),
            )
            .unwrap();
        let code = CString::new(r#"
import inspect
model = fk.Model.from_json(MODEL_JSON)
restriction = fk.Restrict().only_vertices(["V_3_SCALAR_000"])
run = fk.Run(threads=1)
process = model.process([1000], [1000, 1000])
request = process.amplitude(loops=fk.Loops.total(1), restrict=restriction, topology=fk.TopologyRules(max_vertices=3))
stages = []
implicit = request.generate().build(run=fk.Run(progress=lambda p: stages.append(p.stage)))
explicit = process.amplitude(loops=fk.Loops.total(1), restrict=restriction, topology=fk.TopologyRules(max_vertices=3, bridges=graph.Bridges.any(), self_loops=True, zero_flow_edges=False, self_energies=fk.Veto(), tadpoles=fk.Veto(), snails=fk.Veto(), fermion_loops=fk.Range.unrestricted(), factorized_loops=fk.Range.unrestricted())).generate().build()
assert implicit.report.completed and len(implicit) > 0
assert [d.to_json() for d in implicit] == [d.to_json() for d in explicit]
assert all(d.name is None for d in implicit)
assert all(len(group.members) == 1 for group in implicit.groups)
assert fk.Stage.GroupingPreparation not in stages
stages.clear()
grouped = implicit.group(fk.Grouping("up_to_scalar"), run=fk.Run(progress=lambda p: stages.append(p.stage)))
assert grouped.report.completed and len(grouped) > 0
assert fk.Stage.GroupingPreparation in stages
# Bridge restrictions are explicit for both trees and loops.
process = model.process([1000, 1000], [1000, 1000])
base = dict(loops=fk.Loops.total(0), restrict=restriction)
tree = process.amplitude(**base).generate()
assert len(tree) > 0
assert len(process.amplitude(**base, topology=fk.TopologyRules(bridges=graph.Bridges.at_most(0))).generate()) == 0
assert len(process.amplitude(**base, topology=fk.TopologyRules(bridges=graph.Bridges.one_pi_in_loops())).generate()) == len(tree)
# Vacuum requests also require loops and receive the same minimal topology defaults.
vacuum = model.process([], [])
try:
    vacuum.amplitude()
except TypeError:
    pass
else:
    raise AssertionError("vacuum request omitted loops")
a = vacuum.amplitude(loops=fk.Loops.total(2), restrict=restriction).generate().build()
b = vacuum.amplitude(loops=fk.Loops.total(2), restrict=restriction, topology=fk.TopologyRules()).generate().build()
assert [d.to_json() for d in a] == [d.to_json() for d in b]
# Cross-section blob and spectator bounds are explicit; None never installs finite defaults.
forward = process.cross_section(loops=fk.Loops.forward(1), restrict=restriction, topology=fk.TopologyRules(max_vertices=4)).generate().build()
unrestricted = process.cross_section(loops=fk.Loops.forward(1), restrict=restriction, topology=fk.TopologyRules(max_vertices=4), cuts=fk.Cuts(blobs=fk.Range.unrestricted(), spectators=fk.Range.unrestricted())).generate().build()
assert [d.to_json() for d in forward] == [d.to_json() for d in unrestricted]
assert fk.Cuts().blobs.maximum is None and fk.Cuts().spectators.maximum is None
for method in (fk.Process.amplitude, fk.Process.cross_section):
    assert inspect.signature(method).parameters["loops"].default is inspect.Parameter.empty
for make in (lambda: process.amplitude(loops=fk.Loops.forward(0)), lambda: process.cross_section(loops=fk.Loops.total(1))):
    try:
        make()
    except fk.GenerationError:
        pass
    else:
        raise AssertionError("wrong loop domain was accepted")

"#).unwrap();
        py.run(&code, Some(&locals), Some(&locals)).unwrap();
    });
}

#[test]
fn progress_and_partial_filters_use_the_calling_python_thread() {
    Python::initialize();
    Python::attach(|py| {
        let module = PyModule::new(py, "symbolica.community.hepkit").unwrap();
        crate::initialize_feynkit(&module).unwrap();
        let locals = PyDict::new(py);
        locals.set_item("fk", &module).unwrap();
        let graph_module = PyModule::new(py, "symbolica.community.graph").unwrap();
        linnet_py::register(&graph_module).unwrap();
        locals.set_item("graph", &graph_module).unwrap();
        locals.set_item("Graph", py.get_type::<PyGraph>()).unwrap();
        let graph_module = PyModule::new(py, "symbolica.community.graph").unwrap();
        linnet_py::register(&graph_module).unwrap();
        for name in ["Node", "Orientation"] {
            locals
                .set_item(name, graph_module.getattr(name).unwrap())
                .unwrap();
        }
        locals
            .set_item(
                "MODEL_JSON",
                include_str!("../../tests/fixtures/scalars_2p_3p.json"),
            )
            .unwrap();
        locals
            .set_item(
                "SM_MODEL_JSON",
                include_str!("../../../feynkit-model/tests/fixtures/sm.json"),
            )
            .unwrap();
        let code = CString::new(r#"
import threading

model = fk.Model.from_json(MODEL_JSON)
process = model.process([1000], [1000, 1000])
restriction = fk.Restrict().only_vertices(["V_3_SCALAR_000"])
topology = fk.TopologyRules(max_vertices=3, bridges=graph.Bridges.one_pi_in_loops(), self_loops=True,
    self_energies=fk.Veto.all(), tadpoles=fk.Veto.all(), snails=fk.Veto.massless())
caller = threading.get_ident()

for via_model in (True, False):
    def generate(progress=None, cancellation_token=None, **hooks):
        current = model.process([1000], [1000, 1000]) if via_model else process
        request = current.amplitude(loops=fk.Loops.total(1), restrict=restriction, topology=topology, **hooks)
        run = fk.Run(threads=2, progress=progress, cancellation=cancellation_token)
        return request.generate(run=run).build(run=run).group(fk.Grouping("up_to_scalar"), run=run)

    baseline = generate()
    assert baseline.report.completed and len(baseline) > 0
    updates = []
    topologies = []

    def report(p):
        assert threading.get_ident() == caller
        assert isinstance(p, fk.GenerationProgress)
        assert p.completed >= 0
        assert p.total is None or p.completed <= p.total
        updates.append(p)

    snapshots = []
    def keep(view):
        assert threading.get_ident() == caller
        assert type(view).__name__ in ("SkeletonView", "DressingView", "DiagramView")
        assert all(stub is None or stub.data.pdg_code == 1000 for _, _, stub in view.edges())
        topologies.append((view, type(view).__name__))
        if view.is_complete:
            snapshots.append(view.complete_graph())
        return True

    result = generate(progress=report, prune_skeletons=keep, prune_dressings=keep, select=keep)
    assert result.report.completed and len(result) == len(baseline)
    assert {"SkeletonView", "DressingView", "DiagramView"} == {stage for _, stage in topologies}
    assert updates[0].stage == fk.Stage.Skeletons
    assert updates[-1].stage == fk.Stage.Done
    assert updates[-1].completed == len(result)
    stages = [fk.Stage.Skeletons, fk.Stage.Dressing, fk.Stage.Prepare, fk.Stage.Done,
        fk.Stage.Numerators, fk.Stage.GroupingPreparation, fk.Stage.GroupingSamples,
        fk.Stage.GroupingComparison, fk.Stage.Grouping]
    assert list(dict.fromkeys(p.stage for p in updates)) == stages
    for before, after in zip(updates, updates[1:]):
        if before.stage == after.stage:
            assert before.completed <= after.completed
    # Scoped views expire, while explicitly materialized snapshots remain independent.
    for view, _ in topologies:
        try:
            view.n_vertices
        except ReferenceError:
            pass
        else:
            raise AssertionError("callback view escaped its scope")
    assert snapshots
    snapshots[0].add_node(Node(data=42))
    assert len(generate()) == len(baseline)

    updates.clear()
    empty = generate(progress=report, prune_skeletons=lambda view: False)
    assert empty.report.completed and len(empty) == 0
    assert updates[-1].stage == fk.Stage.Done and updates[-1].completed == 0

    for keyword in ("progress", "prune_skeletons", "prune_dressings", "select"):
        expected = RuntimeError("callback failed")
        def fail(*args):
            raise expected
        try:
            generate(**{keyword: fail})
        except RuntimeError as error:
            assert error is expected
        else:
            raise AssertionError("callback exception was lost")
        try:
            generate(**{keyword: 1})
        except TypeError:
            pass
        else:
            raise AssertionError("non-callable callback was accepted")
    for keyword in ("prune_skeletons", "prune_dressings", "select"):
        token = fk.CancellationToken()
        calls = []
        def cancel(view):
            calls.append(type(view).__name__)
            token.cancel()
            return True
        cancelled = generate(cancellation_token=token, **{keyword: cancel})
        assert calls == [{"prune_skeletons": "SkeletonView", "prune_dressings": "DressingView", "select": "DiagramView"}[keyword]]
        assert not cancelled.report.completed
        try:
            generate(**{keyword: lambda view: 1})
        except TypeError:
            pass
        else:
            raise AssertionError("a truthy non-bool verdict was accepted")
    assert len(generate()) == len(baseline)

# Directed snapshot edges carry the electron and point along its flow. External legs are the
# dangling half-edges of their edges: the incoming electron flows into the graph, so its leg is
# the sink, and the incoming positron is an electron flowing out, so its leg is the source.
sm_model = fk.Model.from_json(SM_MODEL_JSON)
legs = set()
def flows(view):
    graph = view.complete_graph()
    if graph is None:
        return True
    for edge in graph.edges():
        if edge.orientation == Orientation.Undirected:
            continue
        assert edge.orientation == Orientation.Default and edge.data.pdg_code == 11
        if edge.source is None:
            legs.add(("sink", edge.sink.data))
        elif edge.sink is None:
            legs.add(("source", edge.source.data))
    return True
annihilation = sm_model.process(["e-", "e+"], ["a", "a"])
assert len(annihilation.amplitude(loops=fk.Loops.total(0), restrict=fk.Restrict().only_particles(["e-",22]), prune_dressings=flows).generate()) == 2
assert legs == {("sink", ("incoming", 0)), ("source", ("incoming", 1))}
"#).unwrap();
        py.run(&code, Some(&locals), Some(&locals)).unwrap();
    });
}

#[cfg(unix)]
#[test]
fn python_signals_interrupt_both_generation_entry_points() {
    Python::initialize();
    Python::attach(|py| {
        let module = PyModule::new(py, "symbolica.community.hepkit").unwrap();
        crate::initialize_feynkit(&module).unwrap();
        let locals = PyDict::new(py);
        locals.set_item("fk", &module).unwrap();
        let graph_module = PyModule::new(py, "symbolica.community.graph").unwrap();
        linnet_py::register(&graph_module).unwrap();
        locals.set_item("graph", &graph_module).unwrap();
        locals
            .set_item(
                "MODEL_JSON",
                include_str!("../../tests/fixtures/scalars_2p_3p.json"),
            )
            .unwrap();
        let code = CString::new(
                r#"
import signal
import threading

model = fk.Model.from_json(MODEL_JSON)
process = model.process([1000], [1000, 1000])
restriction = fk.Restrict().only_vertices(["V_3_SCALAR_000"])
request = process.amplitude(loops=fk.Loops.total(6), restrict=restriction, topology=fk.TopologyRules(max_vertices=13))

def interrupt(signum, frame):
    raise expected_error("generation interrupted")

previous_handler = signal.signal(signal.SIGALRM, interrupt)
try:
    for expected_error in (KeyboardInterrupt, RuntimeError):
        for entry in (request.generate, request.count):
            signal.setitimer(signal.ITIMER_REAL, 0.05)
            try:
                entry(run=fk.Run(threads=2))
            except expected_error as error:
                assert str(error) == "generation interrupted"
            else:
                raise AssertionError("generation ignored its Python signal handler")
            finally:
                signal.setitimer(signal.ITIMER_REAL, 0)
    # The invoking thread must service signals while an opted-in worker releases the GIL.
    released = threading.Event()
    worker_calls = []
    expected = KeyboardInterrupt("slow worker interrupted")
    def release_and_interrupt(signum, frame):
        released.set()
        raise expected
    signal.signal(signal.SIGALRM, release_and_interrupt)
    def slow(view):
        worker_calls.append(1)
        signal.setitimer(signal.ITIMER_REAL, 0.05)
        assert released.wait(2), "calling thread did not service the signal"
        return True
    try:
        process.amplitude(loops=fk.Loops.total(1), restrict=restriction,
            topology=fk.TopologyRules(max_vertices=3), prune_skeletons=slow).generate(
                run=fk.Run(threads=2, thread_safe=True))
    except KeyboardInterrupt as error:
        assert error is expected
    else:
        raise AssertionError("worker callback blocked signal delivery")
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    assert worker_calls == [1]
finally:
    signal.signal(signal.SIGALRM, previous_handler)

# Interruption belongs to one call; subsequent generation must still complete.
result = process.amplitude(loops=fk.Loops.total(0), restrict=restriction).generate(run=fk.Run(threads=2))
assert result.report.completed
assert len(result) == 1
"#,
            )
            .unwrap();
        py.run(&code, Some(&locals), Some(&locals)).unwrap();
    });
}

#[test]
fn process_metadata_and_selectors_are_typed_and_round_trip() {
    Python::initialize();
    Python::attach(|py| {
        let module = PyModule::new(py, "symbolica.community.hepkit").unwrap();
        crate::initialize_feynkit(&module).unwrap();
        let locals = PyDict::new(py);
        locals.set_item("fk", &module).unwrap();
        let graph_module = PyModule::new(py, "symbolica.community.graph").unwrap();
        linnet_py::register(&graph_module).unwrap();
        locals.set_item("graph", &graph_module).unwrap();
        locals
            .set_item(
                "MODEL_JSON",
                include_str!("../../tests/fixtures/scalars_2p_3p.json"),
            )
            .unwrap();
        locals
            .set_item(
                "SM_MODEL_JSON",
                include_str!("../../../feynkit-model/tests/fixtures/sm.json"),
            )
            .unwrap();
        let code = CString::new(
                r#"
by_name = fk.ParticleSelector.by_name("1")
by_pdg = fk.ParticleSelector.by_pdg(1)
model = fk.Model.from_json(MODEL_JSON)
particle = model.particle("scalar_0")

assert by_name.name == "1"
try:
    by_name.pdg
except TypeError:
    pass
else:
    raise AssertionError("named selectors must reject PDG access")
assert by_name.is_name and not by_name.is_pdg
try:
    by_pdg.name
except TypeError:
    pass
else:
    raise AssertionError("PDG selectors must reject name access")
assert by_pdg.pdg == 1
assert by_pdg.is_pdg and not by_pdg.is_name
assert by_name == "1"
assert by_pdg == 1
assert by_name != by_pdg
try:
    by_pdg.pdg = 2
except AttributeError:
    pass
else:
    raise AssertionError("particle selectors must be immutable")

by_name = fk.ParticleSelector.by_name("scalar_1")
by_pdg = fk.ParticleSelector.by_pdg(1001)
mixed_process = model.process([particle, by_name, "scalar_2", 1001], [particle.antiparticle])
assert [
    selector.name if selector.is_name else selector.resolve(model).pdg_code
    for selector in mixed_process.incoming
] == [1000, "scalar_1", "scalar_2", 1001]
assert mixed_process.outgoing_alternatives[0][0].resolve(model).pdg_code == 1000

sm_model = fk.Model.from_json(SM_MODEL_JSON)
gluon = sm_model.particle_by_pdg(21)
gluon_process = sm_model.process([gluon], [gluon.antiparticle])
assert gluon_process.incoming[0].resolve(sm_model).pdg_code == 21
assert gluon_process.outgoing_alternatives[0][0].resolve(sm_model).pdg_code == 21

bottom = sm_model.particle_by_pdg(5)
bottom_process = sm_model.process([bottom], [bottom.antiparticle])
assert bottom_process.incoming[0].resolve(sm_model).pdg_code == 5
assert bottom_process.outgoing_alternatives[0][0].resolve(sm_model).pdg_code == -5

process = model.process([by_name, by_pdg], [fk.ParticleSelector.by_name("scalar_2"), fk.ParticleSelector.by_pdg(1000)])
assert all(isinstance(selector, fk.ParticleSelector) for selector in process.incoming)
assert process.incoming == [by_name, by_pdg]
assert process.outgoing_alternatives[0][0].name == "scalar_2"
assert process.outgoing_alternatives[0][1].pdg == 1000

round_tripped = model.process(process.incoming, process.outgoing_alternatives[0])
assert round_tripped.incoming == process.incoming
assert round_tripped.outgoing_alternatives == process.outgoing_alternatives
assert not hasattr(round_tripped, "particle_veto")
assert not hasattr(process, "loop_count")

selected = sm_model.process(["e-", "e+"], ["a", "a"])
restriction = fk.Restrict().only_particles([sm_model.particle("e-"), 22])
request = selected.amplitude(loops=fk.Loops.total(0), restrict=restriction)
assert len(request.generate()) == 2
assert len(request.generate().build().into_amplitude().diagrams) == 2
assert len(selected.cross_section(loops=fk.Loops.forward(1), restrict=restriction, topology=fk.TopologyRules(max_vertices=4), cuts=fk.Cuts(blobs=fk.Range.exact(1), spectators=fk.Range.exact(0))).generate()) > 0
assert len(selected.amplitude(loops=fk.Loops.total(0), restrict=fk.Restrict().only_particles([])).generate()) == 0
assert len(selected.amplitude(loops=fk.Loops.total(0)).generate()) == 2
assert len(selected.amplitude(loops=fk.Loops.total(0), restrict=restriction.exclude_particles(["e-"])).generate()) == 0
# Resolved aliases share identity; repeated allow clauses intersect after resolution.
assert len(selected.amplitude(loops=fk.Loops.total(0), restrict=restriction.only_particles([11,"a"])).generate()) == 2
assert len(selected.amplitude(loops=fk.Loops.total(0), restrict=restriction.only_particles(["a"])).generate()) == 0
try:
    selected.amplitude(loops=fk.Loops.total(0), restrict=fk.Restrict().only_particles(["missing_particle"]))
except fk.GenerationError:
    pass
else:
    raise AssertionError("particle selection accepted an unknown species")

cross_section = model.process([particle, by_pdg], [by_name])
cross_section = cross_section.with_final_state_alternatives(
    [[particle.antiparticle], [by_name], [by_pdg], ["scalar_2"], [1002]],
)
assert cross_section.incoming[0].resolve(model).pdg_code == 1000
assert cross_section.outgoing_alternatives[0][0].resolve(model).pdg_code == 1000
assert cross_section.outgoing_alternatives[1:] == [
    [by_name],
    [by_pdg],
    [fk.ParticleSelector.by_name("scalar_2")],
    [fk.ParticleSelector.by_pdg(1002)],
]

for particle_veto in ([particle, 1001], ["scalar_0"]):
    assert len(model.process([particle], [particle, particle]).amplitude(loops=fk.Loops.total(0), restrict=fk.Restrict().exclude_particles(particle_veto)).generate()) == 0

foreign_model = fk.Model.from_json(MODEL_JSON.replace("scalar_0", "foreign_scalar_0"))
foreign_particle = foreign_model.particle("foreign_scalar_0")
for make in (
    lambda: model.process([foreign_particle], [foreign_particle.antiparticle]),
    lambda: model.process([particle], [particle,particle]).amplitude(loops=fk.Loops.total(0), restrict=fk.Restrict().only_particles([foreign_particle])),
    lambda: model.process([particle], [particle,particle]).amplitude(loops=fk.Loops.total(0), restrict=fk.Restrict().only_vertices([foreign_model.vertex_rules[0]])),
):
    try:
        make()
    except fk.GenerationError as error:
        assert "model" in str(error)
    else:
        raise AssertionError("foreign model identity was silently discarded")
"#,
            )
            .unwrap();
        py.run(&code, Some(&locals), Some(&locals)).unwrap();
    });
}

/// Typed Python dressing predicates, batched on the calling thread, prune like the same
/// native predicate: the same diagrams with the same names,
/// ids, symmetry and overall factors, in the same order at every thread count.
#[test]
fn python_typed_hooks_prune_like_rust_hooks() {
    use symbolica::atom::AtomCore;

    let model = Arc::new(Model::qcd());
    let ghosts: Vec<_> = model
        .particles()
        .iter()
        .filter(|particle| particle.is_ghost())
        .collect();
    let ghost_ids: Vec<_> = ghosts
        .iter()
        .map(|particle| model.particle_id(&particle.name).unwrap())
        .collect();
    let process = Process::new(["g", "g"], ["g", "g"]);
    let request = Request::amplitude(process.clone(), feynkit_generator::Loops::total(1..=1))
        .with_topology(feynkit_generator::TopologyRules {
            self_loops: false,
            ..Default::default()
        });
    let naming = feynkit_generator::NumeratorOptions::default()
        .with_naming(feynkit_generator::Naming::Named { prefix: "D".into() });
    let summary = |result: &feynkit_generator::BuiltDiagramSet| {
        result
            .diagrams()
            .iter()
            .map(|diagram| {
                (
                    diagram.name().map(str::to_owned),
                    diagram.id().unwrap(),
                    diagram.symmetry_factor(),
                    diagram.overall_factor().to_canonical_string(),
                )
            })
            .collect::<Vec<_>>()
    };
    Python::initialize();
    Python::attach(|py| {
        let locals = PyDict::new(py);
        let pdgs: Vec<_> = ghosts.iter().map(|particle| particle.pdg_code).collect();
        locals.set_item("GHOSTS", pdgs).unwrap();
        let code = CString::new(
                "def keep(view):\n    return all(stub is None or stub.data.pdg_code not in GHOSTS for _, _, stub in view.edges())\n",
            )
            .unwrap();
        py.run(&code, Some(&locals), Some(&locals)).unwrap();
        let keep: Py<PyAny> = locals.get_item("keep").unwrap().unwrap().unbind();
        let mut sequential = None;
        for threads in [1, 8, 32] {
            let ghost_ids = ghost_ids.clone();
            let run = Run::default().threads(threads);
            let rust = request
                .clone()
                .prune_dressings(move |view| {
                    view.edges()
                        .all(|edge| edge.stub.is_none_or(|stub| !ghost_ids.contains(stub.field)))
                })
                .generate(model.clone(), &run)
                .unwrap()
                .build(&naming, &run)
                .unwrap();
            let python_request = PyRequest {
                inner: request.clone(),
                model: model.clone(),
                prune_skeletons: None,
                prune_dressings: Some(keep.clone_ref(py)),
                select: None,
                worker_hooks: Default::default(),
                check_heredity: false,
            };
            let controls = PyRun {
                inner: run.clone(),
                progress: None,
                thread_safe: false,
            };
            let compact = python_request.generate(py, Some(&controls)).unwrap();
            let python = compact.inner.build(&naming, &run).unwrap();
            assert_eq!(summary(&python), summary(&rust), "{threads} threads");
            let sequential = sequential.get_or_insert_with(|| summary(&rust));
            assert_eq!(&summary(&rust), sequential, "{threads} threads");
            assert_eq!(sequential.len(), 277);
        }
    });
}

#[test]
fn worker_permission_keeps_progress_on_caller_and_preserves_failures_and_gc() {
    Python::initialize();
    Python::attach(|py| {
        let module = PyModule::new(py, "symbolica.community.hepkit").unwrap();
        crate::initialize_feynkit(&module).unwrap();
        let locals = PyDict::new(py);
        locals.set_item("fk", &module).unwrap();
        locals
            .set_item(
                "MODEL_JSON",
                include_str!("../../tests/fixtures/scalars_2p_3p.json"),
            )
            .unwrap();
        let code=CString::new(r#"
import gc
import sys
import threading
import weakref
model = fk.Model.from_json(MODEL_JSON)
process = model.process([1000], [1000,1000])
criteria = dict(loops=fk.Loops.total(1), restrict=fk.Restrict().only_vertices(["V_3_SCALAR_000"]), topology=fk.TopologyRules(max_vertices=3))
caller = threading.get_ident()
workers = set()
progress_threads = set()
def keep(view):
    workers.add(threading.get_ident())
    return True
run = fk.Run(threads=2, thread_safe=True, progress=lambda p: progress_threads.add(threading.get_ident()))
result = process.amplitude(**criteria, prune_skeletons=keep, prune_dressings=keep, select=keep).generate(run=run)
assert result.report.completed and len(result) > 0
assert progress_threads == {caller}
assert len(workers) == 1
assert workers == {caller} if sys.platform in ("emscripten","wasi") else caller not in workers
# Reject before allocating an unsendable Graph; retained views expire and collect safely.
retained = []
snapshots = []
def inspect_snapshot(view):
    assert view.n_vertices > 0
    retained.append(view)
    for method in (view.graph, view.complete_graph):
        if sys.platform in ("emscripten", "wasi"):
            snapshot = method()
            if snapshot is not None:
                snapshots.append(snapshot)
        else:
            try:
                method()
            except RuntimeError as error:
                assert "Graph snapshots require the calling thread" in str(error)
            else:
                raise AssertionError("a worker created an unsendable Graph")
    return True
snapshot_result = process.amplitude(**criteria, prune_skeletons=inspect_snapshot,
    prune_dressings=inspect_snapshot, select=inspect_snapshot).generate(run=run)
assert snapshot_result.report.completed and len(snapshot_result) == len(result)
for view in retained:
    try:
        view.graph()
    except ReferenceError:
        pass
    else:
        raise AssertionError("retained worker view did not expire")
for snapshot in snapshots:
    assert snapshot.n_nodes > 0
retained.clear()
del view
gc.collect()
if sys.platform not in ("emscripten", "wasi"):
    nested_errors = []
    def inspect_on_nested_worker(view):
        def inspect():
            try:
                view.graph()
            except RuntimeError as error:
                nested_errors.append(str(error))
        worker = threading.Thread(target=inspect)
        worker.start()
        worker.join()
        assert "Graph snapshots require the calling thread" in nested_errors[-1]
        return False
    assert len(process.amplitude(**criteria, prune_skeletons=inspect_on_nested_worker).generate()) == 0
    assert nested_errors
for hook in ("prune_skeletons","prune_dressings","select"):
    expected = RuntimeError("worker callback failed")
    calls = []
    def fail(view):
        calls.append(1)
        raise expected
    try:
        process.amplitude(**criteria, **{hook: fail}).generate(run=fk.Run(thread_safe=True))
    except RuntimeError as error:
        assert error is expected
    else:
        raise AssertionError("worker exception was lost")
    assert calls == [1]
    token = fk.CancellationToken()
    calls.clear()
    def cancel(view):
        calls.append(1)
        token.cancel()
        return True
    cancelled = process.amplitude(**criteria, **{hook: cancel}).generate(run=fk.Run(thread_safe=True,cancellation=token))
    assert calls == [1] and not cancelled.report.completed
# None continues; explicit False cancels; truthy unrelated values never coerce.
assert process.amplitude(**criteria).generate(run=fk.Run(progress=lambda p: None)).report.completed
assert not process.amplitude(**criteria).generate(run=fk.Run(progress=lambda p: False)).report.completed
# Final progress runs before each native owner seals its completion status.
def stop_at_done(update):
    return update.stage != fk.Stage.Done
request = process.amplitude(**criteria)
assert not request.generate(run=fk.Run(progress=stop_at_done)).report.completed
assert not request.count(run=fk.Run(progress=stop_at_done)).complete
assert not result.build(run=fk.Run(progress=stop_at_done)).report.completed
built = result.build()
assert not built.group(fk.Grouping(), run=fk.Run(progress=stop_at_done)).report.completed
build_updates = []
result.build(run=fk.Run(progress=lambda p: build_updates.append((p.stage, p.completed))))
assert [count for stage, count in build_updates if stage == fk.Stage.Numerators][0] == 0
for invalid in (1,"yes",object()):
    try:
        process.amplitude(**criteria).generate(run=fk.Run(progress=lambda p: invalid))
    except TypeError as error:
        assert "None or bool" in str(error)
    else:
        raise AssertionError("progress accepted a non-bool verdict")
def cycle():
    owner = {}
    def predicate(view):
        assert owner["request"] is not None
        return True
    request = process.amplitude(**criteria, prune_skeletons=predicate)
    owner["request"] = request
    compact = request.generate()
    return weakref.ref(predicate), compact
reference, compact = cycle()
gc.collect()
assert reference() is None and len(compact) > 0
"#).unwrap();
        py.run(&code, Some(&locals), Some(&locals)).unwrap();
    });
}

#[test]
fn staged_python_amplitude_preserves_squared_pairs_after_selection_and_grouping() {
    Python::initialize();
    Python::attach(|py| {
        let module = PyModule::new(py, "symbolica.community.hepkit").unwrap();
        crate::initialize_feynkit(&module).unwrap();
        let locals = PyDict::new(py);
        locals.set_item("fk", &module).unwrap();
        locals
            .set_item(
                "MODEL_JSON",
                include_str!("../../tests/fixtures/scalars_2p_3p.json"),
            )
            .unwrap();
        let code=CString::new(r#"
import copy
import json
definition = json.loads(MODEL_JSON)
low = next(rule for rule in definition["vertex_rules"] if rule["name"] == "V_3_SCALAR_000")
high = copy.deepcopy(low)
high["name"] = "V_3_HIGH"
high["couplings"] = [["HIGH"]]
definition["vertex_rules"] = [low,high]
coupling = copy.deepcopy(definition["couplings"][0])
coupling["name"] = "HIGH"
coupling["orders"] = [["QCD",2]]
definition["couplings"].append(coupling)
model = fk.Model.from_json(json.dumps(definition))
request = model.process(["scalar_0"],["scalar_0","scalar_0"]).amplitude(loops=fk.Loops.total(0), orders=fk.Orders().squared_exact("QCD",3))
compact = request.generate()
assert len(compact) == request.count().diagrams == 2
built = compact.build(fk.NumeratorOptions(naming=fk.Naming.named()))
grouping = fk.Grouping("identical",check_canonical_numerator=True)
grouped = built.group(grouping)
assert len(grouped) == 2
restricted = grouped.into_amplitude().squared().expression().to_expression()
unrestricted = fk.Amplitude(list(grouped)).squared().expression().to_expression()
assert restricted * 2 == unrestricted
for diagram in built:
    selected = built.select(fk.Selection(select=[diagram.name])).group(grouping)
    assert len(selected) == 1
    assert selected.into_amplitude().squared().expression().to_expression() == 0
"#).unwrap();
        py.run(&code, Some(&locals), Some(&locals)).unwrap();
    });
}

#[test]
fn staged_python_selection_uses_run_progress_cancellation_and_errors() {
    Python::initialize();
    Python::attach(|py| {
        let module = PyModule::new(py, "symbolica.community.hepkit").unwrap();
        crate::initialize_feynkit(&module).unwrap();
        let locals = PyDict::new(py);
        locals.set_item("fk", &module).unwrap();
        let code = CString::new(r#"
import threading
model = fk.Model.qed()
built = model.process(["e-","e+"],["a","a"]).amplitude(loops=fk.Loops.total(0)).generate().build()
assert len(built) == 2
calling_thread = threading.get_ident()
updates = []
def record(update):
    assert threading.get_ident() == calling_thread
    updates.append((update.stage,update.completed,update.total))
empty = built.select(fk.Selection(select=[]), run=fk.Run(threads=1,progress=record))
assert len(empty) == 0 and empty.report.completed
assert updates[0] == (fk.Stage.Selection,0,2)
assert updates[-1] == (fk.Stage.Done,0,0)
# Interior updates may coalesce, but the completed scan remains visible.
assert updates[-2] == (fk.Stage.Selection,2,2)
updates.clear()
assert empty.select(fk.Selection(),run=fk.Run(progress=record)).report.completed
assert updates == [(fk.Stage.Selection,0,0),(fk.Stage.Done,0,0)]
def cancel_before_work(update):
    return update.stage != fk.Stage.Selection
cancelled = built.select(fk.Selection(),run=fk.Run(progress=cancel_before_work))
assert len(cancelled) == 0 and not cancelled.report.completed
assert not cancelled.select(fk.Selection()).report.completed
assert len(built) == 2 and built.report.completed
assert not built.select(fk.Selection(),run=fk.Run(progress=lambda p:p.stage != fk.Stage.Done)).report.completed
expected = RuntimeError("selection progress failed")
seen = []
def fail(update):
    seen.append((update.stage,update.completed))
    raise expected
try:
    built.select(fk.Selection(),run=fk.Run(progress=fail))
except RuntimeError as error:
    assert error is expected
else:
    raise AssertionError("selection lost its progress exception")
assert seen == [(fk.Stage.Selection,0)]
try:
    built.select(fk.Selection(),run=fk.Run(progress=lambda _:1))
except TypeError as error:
    assert "None or bool" in str(error)
else:
    raise AssertionError("selection accepted a non-bool progress verdict")
forced = fk.Selection(forced_cuts=[[999]])
try:
    built.select(forced)
except fk.GenerationError:
    pass
else:
    raise AssertionError("complete selection lost its unknown-cut error")
cancelled = built.select(forced,run=fk.Run(progress=cancel_before_work))
assert len(cancelled) == 0 and not cancelled.report.completed
assert not cancelled.select(forced).report.completed
"#).unwrap();
        py.run(&code, Some(&locals), Some(&locals)).unwrap();
    });
}

#[test]
fn request_decorators_preserve_dispatch_gc_and_heredity_contracts() {
    Python::initialize();
    Python::attach(|py| {
        let module = PyModule::new(py, "symbolica.community.hepkit").unwrap();
        crate::initialize_feynkit(&module).unwrap();
        let locals = PyDict::new(py);
        locals
            .set_item("__name__", "request_decorator_contract")
            .unwrap();
        let code = CString::new(include_str!("../../tests/request_decorator_contract.py")).unwrap();
        py.run(&code, Some(&locals), Some(&locals)).unwrap();
        locals
            .get_item("check_request_decorators")
            .unwrap()
            .unwrap()
            .call1((
                &module,
                include_str!("../../tests/fixtures/scalars_2p_3p.json"),
            ))
            .unwrap();
    });
}
