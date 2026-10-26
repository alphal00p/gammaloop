"""Native HEPKit graph and tensor objects in a single Symbolica host."""

# Marimo cell outputs and generated return markers are intentional.
# ruff: noqa: B018, PLR1711
import marimo

__generated_with = "0.25.2"
app = marimo.App(width="full", app_title="HEPKit · Graphs & numerators")


@app.cell
def _():
    from pathlib import Path
    from time import perf_counter

    import marimo as mo
    from shared_gallery import Gallery
    from symbolica import S
    from symbolica.community.hepkit import FeynmanDiagram, Model
    from symbolica.community.tensor import TensorExpression

    return FeynmanDiagram, Gallery, Model, Path, S, TensorExpression, mo, perf_counter


@app.cell(hide_code=True)
def _(Gallery, mo):
    Gallery.header(
        mo,
        "graphs",
        "From a graph to dot products",
        "Follow the native QCD numerator through model dressing, color and Dirac "
        "contractions, and exact scalar products. Every object remains a shared "
        "HEPKit or Symbolica object.",
    )
    return


@app.cell(hide_code=True)
def _(Gallery, mo):
    Gallery.section(
        mo,
        "01 / geometry",
        "Choose the graph",
        "Start small, then inspect the five-loop target.",
    )
    return


@app.cell
def _(mo):
    graph_choice = mo.ui.dropdown(
        {
            "One-loop quark bubble · quick example": "bubble",
            "Five-loop four-rung fixture · geometry": "ladder",
        },
        value="One-loop quark bubble · quick example",
        label="QCD topology",
        full_width=True,
    )
    graph_choice
    return (graph_choice,)


@app.cell
def _(FeynmanDiagram, Gallery, Model, graph_choice):
    model = Model.qcd()
    if graph_choice.value == "bubble":
        topology = """digraph quark_bubble {
            incoming [style=invis]; outgoing [style=invis];
            incoming -> a [particle="g"];
            b -> outgoing [particle="g"];
            a -> b [particle="d", lmb_id=0];
            b -> a [particle="d"];
        }"""
    else:
        topology = (
            Gallery.directory.parents[1] / "tests/fixtures/qcd_ladder/ladder.dot"
        ).read_text()
    diagram = FeynmanDiagram.from_dot(model, topology)
    return diagram, model, topology


@app.cell(hide_code=True)
def _(Gallery, diagram, mo):
    mo.vstack(
        [
            Gallery.metrics(
                mo,
                [
                    (diagram.loop_count, "independent loops"),
                    (len(diagram.internal_edges), "internal particle lines"),
                    (len(diagram.external_edges), "external gluons"),
                ],
            ),
            mo.as_html(diagram.render(momenta=True)),
        ]
    )
    return


@app.cell(hide_code=True)
def _(diagram, mo, topology):
    mo.accordion(
        {
            "Native routing and propagators": mo.ui.table(
                [
                    {
                        "edge": edge.id,
                        "particle": edge.particle_name,
                        "momentum": mo.as_html(edge.momentum_expression(in_lmb=True)),
                    }
                    for edge in diagram.internal_edges
                ],
                selection=None,
                show_download=False,
            ),
            "Compact source DOT": mo.ui.code_editor(
                topology, language="dot", disabled=True
            ),
        }
    )
    return


@app.cell(hide_code=True)
def _(Gallery, mo):
    Gallery.section(
        mo,
        "02 / numerator",
        "Dress, project, contract",
        "The optimized tensor path keeps products factored until contractions need expansion.",
    )
    return


@app.cell
def _(Gallery, graph_choice, mo):
    contract_button = mo.ui.run_button(
        label="Build and contract this numerator", kind="success"
    )
    mo.vstack(
        [
            Gallery.note(
                mo,
                "Explicit computation",
                "The five-loop numerator can take several minutes. This unaveraged metric/color contraction is a tensor demonstration; the next app follows the complete zero-momentum vacuum reduction. The quadratic Taylor fixture is deferred.",
                pending=graph_choice.value == "ladder",
            ),
            contract_button,
        ]
    )
    return (contract_button,)


@app.cell
def _(S, TensorExpression, contract_button, diagram, mo, model, perf_counter):
    mo.stop(not contract_button.value)
    _start = perf_counter()
    dressed = diagram.apply_feynman_rules()
    numerator = dressed.numerator_expression(
        in_lmb=True
    ) * dressed.overall_factor_expression(evaluate=True)
    numerator = (
        model.expand_couplings(numerator)
        .with_lorentz_dimension(S("rustred::d"))
        .collect_factors()
    )
    _slots = numerator.structure.slots()
    if len(_slots) != 4:
        raise ValueError(
            f"Expected two Lorentz and two adjoint slots, found {len(_slots)}"
        )
    projector = TensorExpression.g(_slots[0], _slots[1]) * TensorExpression.g(
        _slots[2], _slots[3]
    )
    projected = projector * numerator
    reduced = projected.simplify_algebra(
        contract="dots", color_substitute_cof_dimension_invariants=True
    )
    contraction_seconds = perf_counter() - _start
    return contraction_seconds, numerator, projected, projector, reduced


@app.cell(hide_code=True)
def _(Gallery, contraction_seconds, mo, numerator, projected, projector, reduced):
    mo.vstack(
        [
            Gallery.metrics(
                mo,
                [
                    (
                        f"{contraction_seconds:.3f} s",
                        "this interactive contraction; not a benchmark sample",
                    )
                ],
            ),
            mo.ui.tabs(
                {
                    "Scalar dot products": Gallery.native_formula(mo, reduced),
                    "Projector": Gallery.native_formula(mo, projector),
                    "Tensor numerator": Gallery.native_formula(mo, numerator),
                    "Projected tensor": Gallery.native_formula(mo, projected),
                }
            ),
            Gallery.note(
                mo,
                "One native tensor workflow",
                "The returned scalar is still a TensorExpression. Model parameters remain symbolic and SU(3) invariants are evaluated. The physical five-loop example specializes the massless quark at zero external momentum and keeps the unaveraged metric/color contraction. Its quadratic Taylor expansion is deferred.",
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(Gallery, mo):
    mo.vstack(
        [
            Gallery.section(mo, "03 / onward", "From tensors to vacuum reduction"),
            mo.md(r"""
        ```python
        from symbolica.community.hepkit.vakint import integral_from_diagram

        family = vacuum_diagram.propagator_family()
        integral = integral_from_diagram(vacuum_diagram, family, scalar_numerator)
        retained = engine.evaluate_integral(integral.to_expression())
        ```

        The vacuum diagram and its family must have matching routing. The complete
        zero-momentum five-loop fixture keeps both current vertices and the
        unaveraged metric/color contraction. Its massless quark numerators are
        specialized at $q=0$; the common mass enters only the vacuum denominators.
        The distinct quadratic Taylor/transverse fixture is deferred and has no
        gallery Run action.
        """),
            Gallery.footer(mo),
        ]
    )
    return


if __name__ == "__main__":
    app.run()
