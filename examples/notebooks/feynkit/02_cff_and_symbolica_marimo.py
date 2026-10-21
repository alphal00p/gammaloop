# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo==0.24.0",
#     "symbolica==3.0.1",
#     "typst==0.15.0",
# ]
# ///

# ruff: noqa: B018, PLR1711 -- marimo uses final expressions and explicit cell returns.

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

    FeynKit turns a loop diagram into a typed CFF representation: accepted
    acyclic orientations, their denominator products, and the energy or hybrid
    surfaces referenced by those products. The final expression lives in the
    same Symbolica kernel as the rest of the package.
    """)
    return


@app.cell
def _():
    from symbolica import S
    from symbolica.community.hepkit import Model

    model = Model.phi3()
    return (Model, S, model)


@app.cell
def _(mo, model, table):
    _generated = model.process(["phi"], [9000001, "phi"]).generate_diagrams(
        loops=(0, 1), max_vertices=3, allow_self_loops=True
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
    `max_orientations` as an explicit combinatorial guard. The returned
    `hepkit.CffRepresentation` displays its orientation and family explorer automatically.
    Click an arrowhead to reverse the displayed energy flow, choose a family
    preview, then Shift-click surface factors to compare their graph regions.
    These selections inspect the result; they do not restrict the full sum.
    """)
    return


@app.cell
def _(diagram):
    cff = diagram.integrate_energy(method="cff", max_orientations=10_000)
    cff
    return (cff,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Inspect surfaces and denominator products

    Energy surfaces contain positive on-shell energies and an external shift.
    H surfaces contain mixed energy signs. The signed coefficients refer to
    physical diagram edge IDs.
    """)
    return


@app.cell
def _(cff, mo, table):
    _surface_rows = [
        {
            "symbol": str(surface.symbol),
            "kind": surface.kind,
            "energy coefficients": {
                edge: str(c) for edge, c in surface.energy_coefficients.items()
            },
            "external shift": {
                edge: str(c) for edge, c in surface.external_shift.items()
            },
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
                str(factor.surface.symbol) for factor in family.factors
            ),
        }
        for index, family in enumerate(_orientation.families)
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

    `to_expression()` includes contour and on-shell factors and expands surfaces
    into symbolic `OSE(edge)` and external energies. Use `expand_surfaces=False`
    to retain surface placeholders. Typed definitions stay in `cff.surfaces`,
    and routed square roots are available through `cff.on_shell_energies`.
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

    Use `diagram.integrate_energy(method="cff", fixed_orientations={edge_id: reversed},
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
