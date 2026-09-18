# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo==0.24.0",
#     "symbolica==3.0.0",
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
    diagrams, and inspect their graph structure.
    """)
    return


@app.cell
def _():
    from pathlib import Path

    import symbolica.community.feynkit as fk

    _data_file = (
        Path(__file__).resolve().parents[3]
        / "crates/feynkit-model/tests/fixtures"
        / "scalars_2p_3p.json"
    )
    model = fk.Model(_data_file)
    return fk, model


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
def _(fk, mo, model, table):
    _card = model.default_parameter_card()
    _card.set("lam", 2.5)
    _updated = model.with_parameter_card(_card)

    try:
        _dependent_coupling = _updated.coupling("SCALAR_COUPLING").value
    except fk.ModelError:
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

    `Model.generate_diagrams` is the main entry point. External states may use
    `Particle` objects directly, avoiding a second lookup or a hand-written PDG
    code. Set `kind="cross_section"` to construct cross-section graph
    structures; this does not numerically integrate phase space. Loop bounds
    are inclusive.
    """)
    return


@app.cell
def _(model, table):
    incoming_particles = [model.particle("scalar_0")]
    outgoing_particles = [
        model.particle_by_pdg(1000),
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
    Pass generation settings directly as keyword arguments. Use an integer
    for an exact loop or coupling order, or a pair for an inclusive range.
    Reuse a configuration with a dictionary and `**kwargs`; each call constructs
    fresh settings, so rerunning a cell does not accumulate filters.
    """)
    return


@app.cell
def _(incoming_particles, mo, model, outgoing_particles, table):
    with mo.status.spinner(title="Generating diagrams") as _status:

        def _report(progress):
            _count = f"{progress.completed:,} processed"
            if progress.total is not None:
                _count += f" / {progress.total:,}"
            _status.update(
                title=progress.stage.replace("_", " ").title(), subtitle=_count
            )

        generated = model.generate_diagrams(
            incoming=incoming_particles,
            outgoing=outgoing_particles,
            loops=(0, 1),
            progress=_report,
            max_vertices=3,
            allow_self_loops=True,
            vertex_allow=["V_3_SCALAR_000"],
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
def _(fk, generated, mo, model, table):
    _loop_diagram = next(
        diagram
        for diagram in generated.diagrams
        if diagram.loop_count == 1
        and all(edge.source != edge.target for edge in diagram.edges)
    )

    from_json = fk.FeynmanDiagram.from_json(model, _loop_diagram.to_json())
    _from_dot = fk.FeynmanDiagram.from_dot(model, _loop_diagram.to_dot())
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

    Configure an `fk.UfoLoader`, for example
    `fk.UfoLoader(restriction_name="massless").load(path)`. It returns a
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
