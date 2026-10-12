# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo==0.24.0",
#     "symbolica==3.0.1",
#     "typst==0.15.0",
# ]
# ///

import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium", app_title="FeynKit: Models, generation, and graphs")


@app.cell
def _():
    from functools import partial

    import marimo as mo

    table = partial(
        mo.ui.table,
        pagination=False,
        selection=None,
        show_download=False,
    )
    return mo, table


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Models, parameters, and diagram generation

    This tutorial treats a model as physics data: inspect particles and
    parameters, make an immutable parameter update, generate tree and loop
    diagrams, and inspect their graph structure. `Model.phi_3_4()` provides
    one real scalar with independent cubic (`g`) and quartic (`lam`) couplings;
    its mass and both couplings default to one.
    """)
    return


@app.cell
def _():
    from symbolica.community.graph import Bridges
    from symbolica.community.hepkit import (
        FeynmanDiagram,
        Loops,
        Model,
        ModelError,
        Naming,
        NumeratorOptions,
        Run,
        TopologyRules,
        Veto,
    )

    model = Model.phi_3_4()
    return (
        Bridges,
        FeynmanDiagram,
        Loops,
        Model,
        ModelError,
        Naming,
        NumeratorOptions,
        Run,
        TopologyRules,
        Veto,
        model,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Particle content

    Lookups are available by model name or PDG code. Spin follows the UFO
    convention $2s+1$, so the scalar entries have `spin == 1`.
    """)
    return


@app.cell
def _(mo, model, table):
    _particle_rows = [
        {
            "name": particle.name,
            "pdg": particle.pdg_code,
            "spin (2s+1)": particle.spin,
            "color rep": particle.color,
            "mass parameter": particle.mass_parameter,
            "massless": particle.is_massless,
        }
        for particle in model.particles
    ]

    table(_particle_rows)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Parameter cards and immutability

    `ParameterCard` is mutable configuration, while
    `Model.with_parameter_card` returns a new model. Without an evaluator,
    changing an external parameter intentionally invalidates dependent
    internal parameters and couplings rather than leaving stale values.
    """)
    return


@app.cell
def _(ModelError, mo, model, table):
    _card = model.default_parameter_card()
    _card.set("lam", 2.5)
    _updated = model.with_parameter_card(_card)

    try:
        _dependent_coupling = _updated.coupling("SCALAR4_COUPLING").value
    except ModelError:
        _dependent_coupling = "not evaluated after the parameter update"

    table(
        [
            {
                "original lambda": model.parameter("lam").value,
                "updated lambda": _updated.parameter("lam").value,
                "dependent coupling after invalidation": _dependent_coupling,
            }
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Model-owned generation and inclusive loop ranges

    `Model.process(...).amplitude(loops=Loops.total(...))` defines the request.
    External states may use `Particle` objects from that model directly.
    A cross-section request instead uses `.cross_section(loops=Loops.forward(...))`;
    it constructs sewn graphs without integrating phase space. Loop bounds are
    inclusive. `generate()` retains compact prepared assignments; `build()`
    explicitly constructs diagrams and numerator fragments.
    """)
    return


@app.cell
def _(model, table):
    incoming_particles = [model.particle("phi")]
    outgoing_particles = [
        model.particle_by_pdg(9000001),
        incoming_particles[0].antiparticle,
    ]

    table(
        [
            {
                "kind": "amplitude",
                "incoming": ", ".join(particle.name for particle in incoming_particles),
                "outgoing": ", ".join(particle.name for particle in outgoing_particles),
                "loop orders": "0 and 1",
            }
        ]
    )
    return incoming_particles, outgoing_particles


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Reuse immutable `Loops`, `Orders`, and `TopologyRules` specifications.
    Execution settings and progress belong to `Run`. This example explicitly
    retains its loop-only bridge policy and cyclic-subgraph vetoes. Naming is
    requested at build time for displayed diagrams; it is absent by default.
    """)
    return


@app.cell
def _(
    Bridges,
    Loops,
    Naming,
    NumeratorOptions,
    Run,
    TopologyRules,
    Veto,
    incoming_particles,
    mo,
    model,
    outgoing_particles,
    table,
):
    with mo.status.spinner(title="Generating diagrams") as _status:

        def _report(progress):
            _count = f"{progress.completed:,} processed"
            if progress.total is not None:
                _count += f" / {progress.total:,}"
            _status.update(
                title=str(progress.stage).removeprefix("Stage."), subtitle=_count
            )
            return True

        _request = model.process(incoming_particles, outgoing_particles).amplitude(
            loops=Loops.total((0, 1)),
            topology=TopologyRules(
                max_vertices=3,
                self_loops=True,
                bridges=Bridges.one_pi_in_loops(),
                self_energies=Veto.all(),
                tadpoles=Veto.all(),
                snails=Veto.massless(),
            ),
        )
        _run = Run(progress=_report)
        generated = _request.generate(run=_run).build(
            NumeratorOptions(naming=Naming.named("D")),
            run=_run,
        )

    table(
        [
            {
                "retained": generated.report.retained_count,
                "loop orders": ", ".join(
                    str(order)
                    for order in sorted(
                        {diagram.loop_count for diagram in generated.diagrams}
                    )
                ),
            }
        ]
    )
    return (generated,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Graph interchange and validation

    Diagrams round-trip through JSON for lossless storage and through DOT for
    graph-tool interoperability. Deserialization requires the canonical model,
    verifies its fingerprint, and restores a self-contained diagram.
    """)
    return


@app.cell
def _(FeynmanDiagram, generated, mo, model, table):
    _loop_diagram = next(
        diagram
        for diagram in generated.diagrams
        if diagram.loop_count == 1
        and all(edge.source != edge.target for edge in diagram.edges)
    )

    from_json = FeynmanDiagram.from_json(model, _loop_diagram.to_json())
    _from_dot = FeynmanDiagram.from_dot(model, _loop_diagram.to_dot())
    from_json.validate()
    _from_dot.validate()

    table(
        [
            {
                "diagram": mo.as_html(from_json),
                "loops": from_json.loop_count,
                "vertices": len(from_json.vertices),
                "edges": len(from_json.edges),
            }
        ],
        column_widths={"diagram": 440},
    )
    return (from_json,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Select interaction regions

    External states live on dangling edges in an amplitude. Every vertex is
    an interaction. Graph analysis calls Rust Linnet directly, and filters
    receive the diagram's physics objects. `diagram.filter(...)` returns a
    `Subgraph` inheriting `FeynmanDiagram`, so the selected region can compute
    its own numerator, denominator, routing, or CFF expression.
    """)
    return


@app.cell
def _(from_json, mo, table):
    selected_region = from_json.filter(edge=lambda edge: not edge.is_external)
    _basis = selected_region.momentum_basis()
    _numerator = selected_region.numerator_expression()
    mo.vstack(
        [
            selected_region,
            table(
                [
                    {
                        "interaction vertices": len(from_json.vertices),
                        "external momentum carriers": len(from_json.external_edges),
                        "selected components": len(
                            selected_region.connected_components()
                        ),
                        "selected loops": selected_region.loop_count,
                        "selected bridge half-edges": len(selected_region.bridges()),
                    }
                ]
            ),
            mo.md("**Selected numerator in loop and external momenta**"),
            _basis.route_expression(_numerator),
        ]
    )
    return (selected_region,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Loop-momentum bases

    A basis identifies loop edges, tree edges, dependent external momenta, and
    the signed loop/external momentum signature carried by every edge.
    """)
    return


@app.cell
def _(from_json, mo, table):
    _bases = from_json.loop_momentum_bases(limit=8)
    _basis = _bases[0]

    _momentum_rows = [
        {"edge": edge, "momentum": signature.format_momentum()}
        for edge, signature in _basis.edge_signatures.items()
    ]

    mo.vstack(
        [
            table(
                [
                    {
                        "number of bases returned": len(_bases),
                        "loop edges": _basis.loop_edges,
                        "tree edges": _basis.tree_edges,
                        "external edges": _basis.external_edges,
                    }
                ]
            ),
            mo.md("**Momentum carried by each edge**"),
            table(_momentum_rows),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Loading UFO models

    For a raw UFO model directory, use Python 3.11 or newer and the
    Symbolica 3-compatible UFO loader pinned in **04 — Loading UFO models**.
    That notebook installs and exercises the optional boundary explicitly.

    Configure an `UfoLoader`, for example
    `UfoLoader(restriction_name="massless").load(path)`. It returns a
    `LoadedModel` containing the normalized `model`, its `parameters`, and
    detailed loader `diagnostics`. Normalized JSON remains the reproducible,
    dependency-free choice for saved analyses.
    """)
    return


@app.cell(hide_code=True)
def _(from_json, generated, mo):
    mo.Html(
        f'<p data-notebook-ready="01_models_and_diagrams_marimo">Validated a graph round-trip and generated {len(generated)} diagrams across the requested loop orders.</p>'
    )
    return


if __name__ == "__main__":
    app.run()
