"""Native HEPKit graph and tensor objects in a single Symbolica host."""

# Marimo cell outputs and generated return markers are intentional.
# ruff: noqa: B018, PLR1711
import marimo

__generated_with = "0.25.2"
app = marimo.App(width="full", app_title="HEPKit · Graphs & numerators")


@app.cell(hide_code=True)
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
        "Start small: the bubble and the three- and four-loop ladders contract "
        "here in under a second; the five-loop ladder is shown as geometry.",
    )
    return


@app.cell(hide_code=True)
def _(mo):
    graph_choice = mo.ui.dropdown(
        {
            "One-loop quark bubble · quick example": 1,
            "Three-loop two-rung ladder (reduced on page 05)": 3,
            "Four-loop three-rung ladder (reduced on page 05)": 4,
            "Five-loop four-rung fixture · geometry (explored on page 03)": 5,
        },
        value="One-loop quark bubble · quick example",
        label="QCD topology",
        full_width=True,
    )
    graph_choice
    return (graph_choice,)


@app.cell(hide_code=True)
def _(FeynmanDiagram, Gallery, Model, graph_choice):
    model = Model.qcd()
    if graph_choice.value == 1:
        topology = """digraph quark_bubble {
            incoming [style=invis]; outgoing [style=invis];
            incoming -> a [particle="g"];
            b -> outgoing [particle="g"];
            a -> b [particle="d", lmb_id=0];
            b -> a [particle="d"];
        }"""
    else:
        topology = Gallery.ladder_dot(graph_choice.value)
    diagram = FeynmanDiagram.from_dot(model, topology)
    return diagram, model, topology


@app.cell(hide_code=True)
def _(Gallery, diagram, mo):
    # The interactive FeynKit drawing runs its own script; inside an iframe it
    # never touches the notebook's DOM, so changing the topology redraws it.
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
            Gallery.diagram_frame(mo, diagram),
        ]
    )
    return


@app.cell(hide_code=True)
def _(Gallery, diagram, mo, topology):
    # Static MathML: an accordion unmounts closed panels, which would dispose
    # a live pager's widget model.
    mo.accordion(
        {
            "Native routing and propagators": Gallery.table(
                mo,
                [
                    {
                        "edge": edge.id,
                        "particle": edge.particle_name,
                        "momentum": Gallery.static_formula(
                            mo, edge.momentum_expression(in_lmb=True)
                        ),
                    }
                    for edge in diagram.internal_edges
                ],
            ),
            # Plain text: marimo's editor has no DOT grammar.
            "Compact source DOT": Gallery.dot_source(mo, topology),
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


@app.cell(hide_code=True)
def _(Gallery, graph_choice, mo):
    _interactive = graph_choice.value <= 4
    contract_button = mo.ui.run_button(
        label="Build and contract this numerator",
        kind="success",
        disabled=not _interactive,
    )
    mo.vstack(
        [
            Gallery.note(
                mo,
                "Explicit computation",
                "This unaveraged metric/color contraction is a tensor demonstration "
                "that runs in this page's kernel. Page 05 reduces the prepared "
                "three- and four-loop ladder inputs end to end in a cancellable "
                "background job."
                if _interactive
                else "The five-loop numerator takes about 16 s of tensor algebra "
                "(longer on a busy server) inside one native call that a page "
                "cannot interrupt, so this page only shows its geometry. Page 03 "
                "explores the five-loop reduction and its recorded run.",
                pending=not _interactive,
            ),
            contract_button,
        ]
    )
    return (contract_button,)


@app.cell(hide_code=True)
def _(
    S,
    TensorExpression,
    contract_button,
    diagram,
    graph_choice,
    mo,
    model,
    perf_counter,
):
    # Never contract the 5L numerator in the kernel, even if re-enabled.
    mo.stop(not contract_button.value or graph_choice.value > 4)
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
    contraction_views = {
        "Scalar dot products": reduced,
        "Projector": projector,
        "Tensor numerator": numerator,
        "Projected tensor": projected,
    }
    view_choice = mo.ui.radio(
        list(contraction_views),
        value="Scalar dot products",
        inline=True,
        label="Show",
    )
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
            view_choice,
        ]
    )
    return contraction_views, view_choice


@app.cell(hide_code=True)
def _(Gallery, contraction_views, mo, view_choice):
    # A fresh live pager per selection, as this cell's own output: prev/next
    # pages, page size, numbered omission markers, + / - / 0 zoom. A single
    # long product does not wrap there (it scrolls sideways): its factors are
    # also listed one per line.
    _value = contraction_views[view_choice.value]
    _factors = Gallery.product_factors(_value)
    mo.vstack(
        [
            Gallery.formula(mo, _value),
            *(
                [
                    mo.accordion(
                        {
                            f"The {len(_factors)} factors of this product, one per line (index names are local to each line)": Gallery.factor_list(
                                mo, _factors
                            )
                        }
                    )
                ]
                if _factors
                else []
            ),
            Gallery.note(
                mo,
                "One native tensor workflow",
                "The returned scalar is still a TensorExpression. Model parameters remain symbolic and SU(3) invariants are evaluated. The physical ladder fixtures specialize the massless quark at zero external momentum (or take its quadratic Taylor coefficient) and keep the unaveraged metric/color contraction.",
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

        The vacuum diagram and its family must have matching routing. The ladder
        fixtures keep both current vertices and the unaveraged metric/color
        contraction. Their massless quark numerators are specialized at $q=0$
        (or Taylor-expanded to second order); the common mass enters only the
        vacuum denominators. Page 05 reduces the prepared 3L and 4L inputs to
        retained masters; page 03 explores the five-loop fixture.
        """),
            Gallery.footer(mo),
        ]
    )
    return


if __name__ == "__main__":
    app.run()
