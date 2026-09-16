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
app = marimo.App(width="medium", app_title="FeynKit: Vacuum tensor reduction")


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
    # Tensor reduction for vacuum graphs

    Tensor vacuum integrals appear after the Taylor expansion used by the
    $R$-operation. Inside a Lorentz-invariant vacuum integral, covariance
    turns a numerator such as $k^\mu k^\nu$ into metric tensors multiplying
    scalar products. FeynKit's native reducer performs that projection in
    symbolic dimension $D$, while keeping repeated momenta in symmetry
    orbits.

    The input and output use the ecosystem's Spenso notation:

    - rank-one tensors end in `spenso::mink(D, index)`;
    - scalar contractions are `spenso::dot(vector, vector)`;
    - residual free-index tensors use `spenso::g(slot, slot)`.

    All returned values are ordinary Symbolica expressions, so the reduced
    scalar numerators can flow directly into simplification or vacuum-integral
    identification.
    """)
    return


@app.cell
def _():
    import math
    import time
    from pathlib import Path

    import symbolica.community.feynkit as fk
    from symbolica import E, S

    D = S("D")
    momentum = S("FeynKit::Momentum")
    external = S("TensorTutorial::p")
    mink = S("spenso::mink")
    return D, E, Path, S, external, fk, math, mink, momentum, time


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Rank two: the basic projector

    We distinguish the integrated momentum by its complete compact Spenso
    form, `Momentum(10, mink(D))`. This is more precise than selecting every
    momentum with the same head and is useful for multiloop graphs.

    Contracting the free tensor with $p_\mu p_\nu$ gives

    $$
      k^\mu k^\nu p_\mu p_\nu
      \longrightarrow
      \frac{(k\!\cdot\!k)(p\!\cdot\!p)}{D}.
    $$

    Here the arrow denotes equality after the symmetric vacuum integration;
    it is not an identity of the unintegrated numerator.
    """)
    return


@app.cell
def _(D, S, external, fk, mink, momentum):
    _mu = S("tr_mu")
    _nu = S("tr_nu")
    _k = momentum(10, mink(D))

    rank_two_reducer = fk.TensorReducer(D).with_integrated_vector(_k)
    rank_two_input = (
        momentum(10, mink(D, _mu))
        * momentum(10, mink(D, _nu))
        * external(0, mink(D, _mu))
        * external(0, mink(D, _nu))
    )
    rank_two_output = rank_two_reducer.reduce(rank_two_input)
    return rank_two_input, rank_two_output, rank_two_reducer


@app.cell
def _(mo, rank_two_input, rank_two_output, table):
    table(
        [
            {
                "stage": "tensor numerator",
                "expression": mo.as_html(rank_two_input.formatted()),
            },
            {
                "stage": "scalar numerator",
                "expression": mo.as_html(rank_two_output.formatted()),
            },
        ],
        column_widths={"stage": 160, "expression": 620},
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Rank six: exploit internal and projector symmetry

    Now the integrated tensor contains two copies of $k$ and four copies of
    $q$, while all six indices meet the same external momentum $p$. The
    labeled metric basis has $5!!=15$ elements, but only two internal
    contraction orbits survive:

    1. $(k\!\cdot\!k)(q\!\cdot\!q)^2$, representing 3 labeled pairings;
    2. $(k\!\cdot\!q)^2(q\!\cdot\!q)$, representing 12 pairings.

    FeynKit groups these before constructing the answer. The same mechanism is
    what makes repeated loop momenta and symmetric external projectors cheap.
    """)
    return


@app.cell
def _(D, S, external, fk, mink, momentum):
    _indices = [S(f"tr_rank6_mu_{position}") for position in range(6)]
    _k = momentum(20, mink(D))
    _q = momentum(21, mink(D))
    _integrated = [
        momentum(20, mink(D, _indices[0])),
        momentum(20, mink(D, _indices[1])),
        *(momentum(21, mink(D, index)) for index in _indices[2:]),
    ]
    _projector = [external(0, mink(D, index)) for index in _indices]

    _factors = [*_integrated, *_projector]
    rank_six_input = _factors[0]
    for _factor in _factors[1:]:
        rank_six_input *= _factor

    rank_six_reducer = (
        fk.TensorReducer(D).with_integrated_vector(_k).with_integrated_vector(_q)
    )
    rank_six_output = rank_six_reducer.reduce(rank_six_input)
    return rank_six_input, rank_six_output


@app.cell
def _(mo, rank_six_output, table):
    table(
        [
            {
                "rank": 6,
                "labeled pairings": 15,
                "internal orbits": 2,
                "reduced numerator": mo.as_html(rank_six_output.formatted()),
            }
        ],
        column_widths={"reduced numerator": 680},
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Rank twenty: a genuinely high-rank projection

    A naive rank-20 projector has $19!!=654{,}729{,}075$ metric pairings.
    Orthogonal Weingarten coefficients depend only on integer-partition (coset)
    type, so even a completely general coefficient table has only
    $p(10)=42$ entries. Partial symmetry is cheaper still: a fully contracted
    rank-20 projection with multiplicities `[2, 18]` on both sides has two
    contraction orbits per side. FeynKit solves the resulting $2\!\times\!2$
    orbit-Gram system instead of enumerating the square of all labeled
    pairings.

    The displayed example is maximally symmetric: both sides contain one
    repeated vector, so its answer remains a single contraction-orbit term.

    The result remains compact and exact in $D$; no numerical dimension is
    substituted.
    """)
    return


@app.cell
def _(D, S, external, fk, math, mink, momentum, time):
    _indices = [S(f"tr_rank20_mu_{position}") for position in range(20)]
    _compact_loop_momentum = momentum(30, mink(D))
    _rank_twenty_factors = []
    for _index in _indices:
        _rank_twenty_factors.extend(
            [
                momentum(30, mink(D, _index)),
                external(0, mink(D, _index)),
            ]
        )
    _rank_twenty_input = _rank_twenty_factors[0]
    for _factor in _rank_twenty_factors[1:]:
        _rank_twenty_input *= _factor

    _reducer = fk.TensorReducer(D).with_integrated_vector(_compact_loop_momentum)
    _started = time.perf_counter()
    rank_twenty_output = _reducer.reduce(_rank_twenty_input)
    rank_twenty_seconds = time.perf_counter() - _started
    rank_twenty_pairings = math.prod(range(1, 20, 2))
    return rank_twenty_output, rank_twenty_pairings, rank_twenty_seconds


@app.cell
def _(
    mo,
    rank_twenty_output,
    rank_twenty_pairings,
    rank_twenty_seconds,
    table,
):
    table(
        [
            {
                "rank": 20,
                "naive metric terms": f"{rank_twenty_pairings:,}",
                "coefficient classes": 42,
                "output orbits": 1,
                "wall time (s)": round(rank_twenty_seconds, 6),
                "reduced numerator": mo.as_html(rank_twenty_output.formatted()),
            }
        ],
        column_widths={"reduced numerator": 680},
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Applying the reducer to a FeynKit graph

    A finalized vacuum diagram exposes the same operation directly. The cell
    below generates the two-loop pure-gluon theta vacuum graph from FeynKit's
    normalized Standard Model. Its two three-gluon vertices give a rank-two
    momentum numerator. Since the graph has no external legs, every native
    `FeynKit::Momentum` is an integrated vacuum momentum and the convenience
    reducer is unambiguous.

    The reducer dimension must match the graph's `spenso::mink` slots; native
    generated rules currently carry `4`, while a vacuum expression prepared
    directly in symbolic dimension uses `D` as above. Use exact
    `with_integrated_vector(...)` selectors when a numerator also contains
    external objects with the `FeynKit::Momentum` head.

    `reduce_tensor_numerator` reduces the product of the stored numerator and
    external-state projector, returning the compact sum as one Expression. The
    separately stored scalar numerator prefactor is not part of that tensor
    projection.

    `reduce_tensor_graphs` instead returns one native `FeynmanDiagram` per
    nonzero scalar term, with the projector coefficient absorbed into its
    numerator. Every contribution retains the scalar numerator prefactor. The
    external-state projector has been consumed by the reduction and is reset
    to one, so it cannot be applied twice. Graph splitting requires a fully
    contracted projection and raises `fk.TensorReductionError` if Lorentz
    indices remain; use `reduce_tensor_numerator` when tensor-valued output is
    intentional. The original diagram is not mutated. Its topology,
    denominators, and topology ID are preserved in every contribution, while
    deterministic `.tensor[index]` names distinguish them. From there one can
    canonicalize scalar products, identify the vacuum topology, and extract
    the local counterterm required by the $R$-operation.
    """)
    return


@app.cell
def _(E, Path, fk, mo, table):
    _model_path = (
        Path(__file__).resolve().parents[3]
        / "crates/feynkit-model/tests/fixtures"
        / "sm.json"
    )
    _vacuum_model = fk.Model(_model_path)

    _vacuum_options = fk.GenerationOptions(max_vertices=2)
    _vacuum_options.set_coupling_orders({"QCD": (2, 2), "QED": (0, 0)})
    _vacuum_options.add_particle_veto(
        [
            _particle
            for _particle in _vacuum_model.particles
            if abs(_particle.pdg_code) != 21
        ]
    )
    _vacuum_result = _vacuum_model.generate_diagrams(
        incoming=[],
        outgoing=[],
        loops=2,
        options=_vacuum_options,
    )
    if len(_vacuum_result) != 1:
        raise RuntimeError("expected one pure-gluon theta vacuum graph")

    vacuum_diagram = _vacuum_result.diagrams[0]
    vacuum_reducer = fk.TensorReducer.feynkit(E("4"))
    reduced_diagram_numerator = vacuum_diagram.reduce_tensor_numerator(vacuum_reducer)
    scalar_vacuum_graphs = vacuum_diagram.reduce_tensor_graphs(vacuum_reducer)

    table(
        [
            {
                "vacuum graph": mo.as_html(vacuum_diagram),
                "tensor rank": 2,
                "scalar graph terms": len(scalar_vacuum_graphs),
                "Lorentz-scalar numerator": mo.as_html(
                    reduced_diagram_numerator.formatted()
                ),
                "first scalar-graph numerator": mo.as_html(
                    scalar_vacuum_graphs[0].numerator_expression().formatted()
                ),
            }
        ],
        column_widths={
            "vacuum graph": 390,
            "tensor rank": 100,
            "scalar graph terms": 130,
            "Lorentz-scalar numerator": 650,
            "first scalar-graph numerator": 650,
        },
    )
    return (
        reduced_diagram_numerator,
        scalar_vacuum_graphs,
        vacuum_diagram,
        vacuum_reducer,
    )


@app.cell(hide_code=True)
def _(mo, rank_twenty_output, scalar_vacuum_graphs):
    mo.Html(
        f'<p data-notebook-ready="07_tensor_reduction_marimo">Reduced the rank-twenty projector and produced {len(scalar_vacuum_graphs)} scalar vacuum-graph terms.</p>'
    )
    return


if __name__ == "__main__":
    app.run()
