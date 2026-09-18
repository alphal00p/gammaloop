# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo==0.24.0",
#     "symbolica==3.0.0",
#     "typst==0.15.0",
#     "ufo-model-loader @ git+https://github.com/alphal00p/ufo_model_loader.git@70ddee6b416f8c8b340e0d087646d77095c5d24b",
# ]
# ///

import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium", app_title="FeynKit: Loading UFO models")


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
    # Loading a raw UFO model

    This tutorial crosses FeynKit's optional UFO boundary: a conventional Python UFO package is normalized into the same typed `Model` used by the core diagram and expression APIs. We use a tiny scalar theory so every diagnostic count is easy to inspect.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Requirements

    Raw UFO import requires **Python 3.11 or newer** and the Symbolica 3-compatible development revision of `ufo-model-loader`. With the shared Symbolica host installed, use:

    ```bash
    python -m pip install --no-deps 'ufo-model-loader @ git+https://github.com/alphal00p/ufo_model_loader.git@70ddee6b416f8c8b340e0d087646d77095c5d24b'
    ```

    The UFO loader is optional for the other showcases. They share the Marimo and Typst dependencies; normalized JSON loading itself needs no UFO tooling.
    """)
    return


@app.cell
def _():
    import os
    from pathlib import Path

    import symbolica.community.feynkit as fk

    UFO_MODEL = Path(__file__).resolve().parents[3] / "assets/models/ufo/scalars"
    return UFO_MODEL, fk, os


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Normalize the UFO package

    `UfoLoader.load` returns three typed results together: `model` for physics operations, `parameters` for the applied numerical inputs, and `diagnostics` recording exactly how normalization was requested. `restriction_name="default"` applies `restrict_default.dat`.

    The environment variable below belongs to this teaching fixture—not to FeynKit generally. It limits this generated scalar model to two- and three-point interactions, keeping the tutorial fast and deterministic.
    """)
    return


@app.cell
def _(UFO_MODEL, fk, mo, os):
    import importlib.util
    import sys

    _requirements_met = (
        sys.version_info >= (3, 11)
        and importlib.util.find_spec("ufo_model_loader") is not None
    )
    mo.stop(
        not _requirements_met,
        mo.callout(
            mo.md(
                "This optional tutorial needs Python 3.11 or newer and "
                "the pinned `ufo-model-loader` revision above. Install it, then rerun "
                "this cell."
            ),
            kind="warn",
            title="UFO loader unavailable",
        ),
    )

    _interaction_key = "UFO_SCALARS_MODEL_N_POINT_INTERACTIONS"
    _previous_interactions = os.environ.get(_interaction_key)
    os.environ[_interaction_key] = "2,3"
    try:
        loaded = fk.UfoLoader(restriction_name="default").load(UFO_MODEL)
    finally:
        if _previous_interactions is None:
            os.environ.pop(_interaction_key, None)
        else:
            os.environ[_interaction_key] = _previous_interactions

    model = loaded.model
    parameters = loaded.parameters
    diagnostics = loaded.diagnostics
    return diagnostics, model, parameters


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Check normalization diagnostics

    These counts describe the normalized model actually handed to Rust. They are useful provenance when comparing model restrictions or diagnosing an unexpected set of interaction rules.
    """)
    return


@app.cell
def _(diagnostics):
    counts = {
        "orders": diagnostics.order_count,
        "model_parameters": diagnostics.model_parameter_count,
        "particles": diagnostics.particle_count,
        "propagators": diagnostics.propagator_count,
        "lorentz_structures": diagnostics.lorentz_structure_count,
        "couplings": diagnostics.coupling_count,
        "vertices": diagnostics.vertex_rule_count,
        "functions": diagnostics.function_count,
        "form_factors": diagnostics.form_factor_count,
        "parameter_values": diagnostics.parameter_value_count,
    }
    {
        "counts": counts,
        "source": diagnostics.source,
        "restriction": diagnostics.restriction_name,
        "simplified": diagnostics.simplify_model,
        "wrapped Lorentz indices": (diagnostics.wrap_indices_in_lorentz_structures),
    }
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Inspect the typed model and restriction card

    The loader has now left the Python UFO object model behind. Particle, parameter, coupling, and vertex lookups use FeynKit's typed API. Zero widths in the restriction card become internal zero parameters during simplification, while the scalar coupling, masses, and compatibility inputs remain inspectable in the parameter card.
    """)
    return


@app.cell
def _(model):
    particle_rows = [
        {
            "name": particle.name,
            "pdg": particle.pdg_code,
            "mass_parameter": particle.mass_parameter,
            "massless": particle.is_massless,
        }
        for particle in model.particles
    ]
    particle_rows
    return


@app.cell
def _(parameters):
    parameter_values = dict(parameters.items())
    parameter_values
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Reuse the model in diagram generation

    Nothing downstream is UFO-specific. Generate diagrams directly from the
    loaded `model`; selectors, generation options, graph inspection, CFF
    construction, and Symbolica expressions all compose unchanged.
    """)
    return


@app.cell
def _(model):
    generated = model.generate_diagrams(
        incoming=["scalar_0"],
        outgoing=["scalar_0", "scalar_0"],
        loops=0,
        max_vertices=3,
        vertex_allow=["V_3_SCALAR_000"],
    )
    {
        "diagrams": len(generated.diagrams),
        "topologies_considered": generated.report.topology_count,
        "first_diagram": generated.diagrams[0].name,
    }
    return (generated,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Working with your own UFO

    Point `UfoLoader.load` at the UFO package directory and choose a
    `restriction_name` matching `restrict_<name>.dat`; use `None` to accept the
    package's default behavior. Keep `simplify_model=True` for the usual removal
    of zero contributions. Disable `wrap_indices_in_lorentz_structures` only
    when interoperating with code that specifically expects unwrapped UFO
    syntax.

    For reproducible production workflows, normalize once with the UFO bridge and serialize the resulting `Model` to JSON; later sessions can load that JSON without the optional Python dependency.
    """)
    return


@app.cell(hide_code=True)
def _(generated, mo):
    mo.Html(
        f'<p data-notebook-ready="04_ufo_loading_marimo">Loaded the UFO model and generated {len(generated)} tree diagrams.</p>'
    )
    return


if __name__ == "__main__":
    app.run()
