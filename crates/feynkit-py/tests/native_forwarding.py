"""Shared standalone/Community native forwarding contract; injected module objects."""

import builtins
import importlib
import sys

if standalone_host:
    assert not any(
        name == "symbolica" or name.startswith("symbolica.") for name in sys.modules
    ), "standalone forwarding must not fabricate canonical package aliases"

forwarding_model = fk.Model.from_json(MODEL_JSON)
forwarding_process = forwarding_model.process(
    ["scalar_0"],
    ["scalar_0", "scalar_0"],
    vertex_allow=["V_3_SCALAR_000"],
    particle_veto=["scalar_1", "scalar_2"],
)
forwarding_diagram = next(
    item
    for item in forwarding_process.generate_diagrams(
        loops=1,
        max_vertices=3,
        allow_self_loops=True,
        progress=None,
    ).diagrams
    if all(edge.source != edge.target for edge in item.edges)
)
forwarding_cff = forwarding_diagram.integrate_energy(method="cff")
forwarding_ltd = forwarding_diagram.integrate_energy(method="ltd")
forwarding_amplitude = fk.Amplitude([forwarding_diagram])
forwarding_basis = forwarding_diagram.loop_momentum_basis
forwarding_q = core.Expression.parse(
    f"gammalooprs::Q({forwarding_diagram.edges[0].id},spenso::cind(0))"
)
forwarding_source = forwarding_diagram.to_json()

original_import = builtins.__import__
original_import_module = importlib.import_module


def guarded_import(name, *args, **kwargs):
    if name == "symbolica" or name.startswith("symbolica."):
        raise AssertionError(
            f"native forwarding imported canonical Python module {name}"
        )
    return original_import(name, *args, **kwargs)


def guarded_import_module(name, *args, **kwargs):
    if name == "symbolica" or name.startswith("symbolica."):
        raise AssertionError(
            f"native forwarding imported canonical Python module {name}"
        )
    return original_import_module(name, *args, **kwargs)


builtins.__import__ = guarded_import
importlib.import_module = guarded_import_module
try:
    # Prove the guard is active even if Community aliases already exist.
    try:
        __import__("symbolica.core")
    except AssertionError:
        pass
    else:
        raise AssertionError("canonical import guard was not active")

    replacements = forwarding_basis.momentum_replacements()
    assert replacements and all(type(rule) is core.Replacement for rule in replacements)
    assert forwarding_q.replace_multiple(
        replacements
    ) == forwarding_basis.route_expression(forwarding_q)
    for owner in (forwarding_diagram, forwarding_cff, forwarding_ltd):
        expression = (
            owner.numerator_expression()
            if owner is forwarding_diagram
            else owner.to_expression()
        )
        assert type(expression) is tensor.TensorExpression
        assert type(expression.to_expression()) is core.Expression
        before = bytes(expression.to_expression().__getstate__())
        assert "<math" in expression.to_html()
        assert bytes(expression.to_expression().__getstate__()) == before

    for owner in (
        forwarding_diagram,
        forwarding_cff,
        forwarding_ltd,
        forwarding_amplitude,
    ):
        view = owner.render()
        assert type(view) is graph.DiagramView
        request = view.measurement_request()
        assert request["labels"]
        # EXTERNAL METRICS PROTOCOL FIXTURE in points; no browser measurement.
        snapshot = view.measured(
            measurements={
                "id": request["id"],
                "metrics": [
                    {"id": label["id"], "width": 24.0, "height": 12.0, "baseline": 9.0}
                    for label in request["labels"]
                ],
            }
        )
        assert type(snapshot) is graph.DiagramRender
        assert "<math" in snapshot.to_svg()
        assert view.measured().to_svg() == snapshot.to_svg()
        assert view.measured().to_html() == snapshot.to_html()
    assert forwarding_diagram.to_json() == forwarding_source
finally:
    builtins.__import__ = original_import
    importlib.import_module = original_import_module
