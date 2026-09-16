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
app = marimo.App(width="medium", app_title="FeynKit: Cross-Free Families")


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
    # Cross-Free Families and Symbolica

    FeynKit turns a loop diagram into a typed Cross-Free Family (CFF): accepted
    acyclic orientations, their denominator products, and the energy or hybrid
    surfaces referenced by those products. The final expression lives in the
    same Symbolica kernel as the rest of the package.
    """)
    return


@app.cell
def _():
    from pathlib import Path

    import symbolica.community.feynkit as fk
    from symbolica import S

    _data_file = (
        Path(__file__).resolve().parents[3]
        / "crates/feynkit-model/tests/fixtures"
        / "scalars_2p_3p.json"
    )
    model = fk.Model(_data_file)
    return S, fk, model


@app.cell
def _(fk, mo, model, table):
    _options = fk.GenerationOptions(max_vertices=3, allow_self_loops=True)
    _options.add_vertex_allow(["V_3_SCALAR_000"])
    _generated = model.generate_diagrams(
        ["scalar_0"],
        [1000, "scalar_0"],
        loops=(0, 1),
        options=_options,
    )

    diagram = next(
        item
        for item in _generated.diagrams
        if item.loop_count == 1
        and all(edge.source != edge.target for edge in item.edges)
    )
    diagram.validate()
    table(
        [
            {
                "diagram": mo.as_html(diagram),
                "loops": diagram.loop_count,
                "vertices": len(diagram.vertices),
                "edges": len(diagram.edges),
            }
        ],
        column_widths={"diagram": 440},
    )
    return (diagram,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Edge IDs are the control surface

    CFF orientation constraints, contractions, and initial-state declarations
    all refer to the stable edge IDs shown below.
    """)
    return


@app.cell
def _(diagram, mo, table):
    _edge_rows = [
        {
            "id": edge.id,
            "from": edge.source,
            "to": edge.target,
            "particle": edge.particle_name,
        }
        for edge in diagram.edges
    ]
    table(_edge_rows)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Generate the CFF

    The diagram owns this conversion. For large graphs, set
    `max_orientations` as an explicit combinatorial guard.
    """)
    return


@app.cell
def _(diagram, mo, table):
    cff = diagram.build_cff(max_orientations=10_000)
    _report = cff.report

    table(
        [
            {
                "candidate orientations": _report.candidate_orientations,
                "acyclic orientations": _report.acyclic_orientations,
                "unfolded terms": _report.unfolded_terms,
                "interned surfaces": _report.interned_surfaces,
            }
        ]
    )
    return (cff,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Inspect surfaces and denominator products

    Energy surfaces contain positive on-shell energies and an external shift.
    Hybrid (H) surfaces may also contain negative energies. Special
    unit/infinite surfaces have no Symbolica symbol.
    """)
    return


@app.cell
def _(cff, mo, table):
    _surface_rows = [
        {
            "symbol": surface.symbol_name,
            "kind": surface.kind,
            "positive energies": surface.positive_energies,
            "negative energies": surface.negative_energies,
            "external shift": surface.external_shift,
            "vertices": surface.vertices,
        }
        for surface in cff.surfaces
    ]
    table(_surface_rows)
    return


@app.cell
def _(cff, mo, table):
    _orientation = cff.orientations[0]
    _product_rows = [
        {
            "product": index,
            "surfaces": " × ".join(
                surface.symbol_name or surface.kind for surface in product
            ),
        }
        for index, product in enumerate(_orientation.denominator_products())
    ]
    mo.vstack(
        [
            mo.md(f"**Orientation {_orientation.id}**"),
            table(
                [
                    {
                        "edge": edge,
                        "orientation": orientation,
                    }
                    for edge, orientation in _orientation.edge_orientations
                ]
            ),
            mo.md("**Denominator products**"),
            table(_product_rows),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Continue symbolically

    `to_expression()` maps CFF surface IDs to names such as `feynkit::E0` and
    `feynkit::H0`. Their typed definitions remain available in `cff.surfaces`,
    so symbolic algebra and physics metadata stay connected.
    """)
    return


@app.cell
def _(S, cff, mo, table):
    _expression = cff.to_expression()
    _weighted = S("g") ** 2 * _expression

    table(
        [
            {
                "expression": mo.as_html(_expression.formatted()),
                "weighted expression": mo.as_html(_weighted.formatted()),
            }
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Advanced controls

    Use `diagram.build_cff(fixed_orientations={edge_id: reversed},
    contracted_edges=[...], initial_state_edges=[...])` for constrained
    constructions. Here `False` keeps an edge's stored direction and `True`
    reverses it. Apply constraints only after inspecting the diagram's edge
    table.
    """)
    return


@app.cell(hide_code=True)
def _(cff, mo):
    mo.Html(
        f'<p data-notebook-ready="02_cff_and_symbolica_marimo">Built {len(cff.orientations)} CFF orientations referencing {len(cff.surfaces)} surfaces.</p>'
    )
    return


if __name__ == "__main__":
    app.run()
