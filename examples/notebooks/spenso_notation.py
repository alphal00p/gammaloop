"""An executable guide to Spenso's symbolic and component notation."""

# ruff: noqa: PLR1711 -- marimo uses explicit cell returns.

import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium", app_title="Reading Spenso tensor notation")


@app.cell
def _():
    import marimo as mo
    from symbolica import Expression, S, Symbol
    from symbolica.community.tensor import (
        AUTO,
        DisplaySettings,
        Representation,
        TensorExpression,
        TensorLibrary,
        TensorName,
        dot,
        to_html,
    )

    def notation(expression: Expression, settings: DisplaySettings | None = None):
        """Use Spenso's renderer for tensors and scalar component expressions."""
        return mo.Html(to_html(expression, settings=settings))

    return (
        AUTO,
        DisplaySettings,
        Representation,
        S,
        Symbol,
        TensorExpression,
        TensorLibrary,
        TensorName,
        dot,
        mo,
        notation,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Reading Spenso tensor notation

    A tensor has **indices in specified vector spaces**, described by Spenso's
    representations. The examples below explain ordinary index notation,
    Schoonschip notation for vector contractions, scalar products and Dirac
    matrices. They use the actual Spenso printers and remain editable expressions.

    Spenso also has its own display conventions for unassigned and contracted
    indices. In its API, a **port** is a tensor axis that can be connected to
    another tensor. The shapes below indicate what has happened to that axis.

    Each mathematical example below is rendered from a Spenso object. The
    captions show how to construct it; the prose explains the resulting display.
    Index letters label tensor axes. Unresolved ports and positions contracted
    with vectors have their own markers, demonstrated below. Contracted positions
    do not count towards the result's rank. A slash denotes a Dirac contraction,
    brackets an ordered matrix product, and a trace sums its matrix indices.

    Display settings change presentation. Calling `contract()` changes the
    symbolic notation while preserving the represented tensor. Calling
    `simplify_algebra()` additionally applies tensor-algebra identities.
    """)
    return


@app.cell
def _(Representation, TensorExpression, TensorName):
    lorentz = Representation.mink(4)
    spinor = Representation.bis(4)
    p = TensorName.vector("notation::p")(lorentz)
    q = TensorName.vector("notation::q")(lorentz)
    Q = TensorName.vector("notation::Q")(spinor)
    A = TensorName("notation::A")(spinor, spinor, lorentz)
    B = TensorName("notation::B")(lorentz, lorentz)
    gamma = TensorExpression.dirac_gamma(4)
    return A, B, Q, gamma, lorentz, p, q, spinor


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 1. Vectors: names, indices and labels

    A vector is a **rank-one tensor**. `TensorName.vector("p")` declares the
    name; applying it to a representation constructs the vector:
    `p = TensorName.vector("p")(Representation.mink(4))`. This is a
    four-dimensional Minkowski vector. Its metric has one positive and three
    negative directions. The same vector constructor works for Euclidean
    vectors, spinors and other index spaces.

    The default head is the name `p`. Its slot carries the index information:
    the unresolved vector has a hollow square; `p("mu")` assigns the abstract
    Lorentz index shown in the second example. This is the whole indexed vector,
    rather than its numerical component at some coordinate.

    A vector family can also have scalar arguments. For
    `P = TensorName.vector("P")`, `P(7, lorentz)` displays a family label
    underneath the head. `P(7, lorentz("mu"))` additionally assigns the index.
    **The 7 identifies a member of the family; the index labels its tensor
    axis.** Component coordinates are introduced later through tensor data.

    A sum of vectors in the same space has one axis too, as shown by
    `(p + q)("mu")`. A dot contracts that axis, while supplying the vector
    to another tensor fills one of its ports. Both operations are shown below.

    Arrows and bold heads are optional print mappings, rather than different
    tensor operations. They keep the same index-space and slot conventions.
    """)
    return


@app.cell
def _(TensorName, lorentz, mo, p, q):
    _P = TensorName.vector("notation::P")
    _family = _P(7, lorentz)
    _family_indexed = _P(7, lorentz("mu"))
    _sum = (p + q)("mu")
    _arrow = TensorName.vector(
        "notation::v_arrow", print={"typst": "arrow(v)", "latex": r"\vec{v}"}
    )(lorentz)
    _bold = TensorName.vector(
        "notation::v_bold", print={"typst": "bold(v)", "latex": r"\mathbf{v}"}
    )(lorentz)
    assert _family_indexed == _family("mu")
    assert _sum.rank == _family.rank == _family_indexed.rank == 1
    mo.vstack(
        [
            mo.hstack(
                [
                    mo.vstack([mo.md("**Unresolved:** `p`"), p]),
                    mo.vstack([mo.md('**Indexed:** `p("mu")`'), p("mu")]),
                    mo.vstack([mo.md("**Family member:** `P(7, lorentz)`"), _family]),
                    mo.vstack(
                        [
                            mo.md('**Family member and index:** `P(7, lorentz("mu"))`'),
                            _family_indexed,
                        ]
                    ),
                ],
                wrap=True,
                align="start",
            ),
            mo.hstack(
                [
                    mo.vstack([mo.md('**Vector sum:** `(p + q)("mu")`'), _sum]),
                    mo.vstack(
                        [
                            mo.md('**Optional arrow:** `print={"typst": "arrow(v)"}`'),
                            _arrow("mu"),
                        ]
                    ),
                    mo.vstack(
                        [
                            mo.md('**Optional bold:** `print={"typst": "bold(v)"}`'),
                            _bold("mu"),
                        ]
                    ),
                ],
                wrap=True,
                align="start",
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 2. Open ports and index labels

    `A` has two bispinor axes and one Lorentz axis. Without assigned labels,
    its three hollow squares identify unresolved ports. Hover a square to see
    its axis and representation. Subscripts on repeated squares distinguish
    the axes; they are **axis numbers, not component coordinates**.

    `A("a", "b", "mu")` assigns labels to all three axes. `AUTO` leaves an axis
    available for a later connection. Both a labelled free index and an
    unresolved port are external axes: each example below has rank three.
    """)
    return


@app.cell
def _(A, AUTO, mo):
    mo.hstack(
        [
            mo.vstack([mo.md("**Unresolved:** `A`"), A]),
            mo.vstack([mo.md('**Labelled:** `A("a", "b", "mu")'), A("a", "b", "mu")]),
            mo.vstack(
                [mo.md('**Partly labelled:** `A("a", AUTO, "mu")'), A("a", AUTO, "mu")]
            ),
        ],
        wrap=True,
        align="start",
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 3. Index alphabets, graph-derived indices and independent copies

    The default alphabet display gives generated Lorentz indices Greek letters
    and bispinor indices Latin letters. These letters label axes; their position
    follows the representation's display convention. In particular, the two
    lower bispinor labels on a gamma matrix mean its row and column axes.

    A Feynman-diagram numerator can use **indices derived from its graph**.
    For example, `hedge(2, 1)` identifies half-edge 2 and local index label 1.
    A half-edge is one endpoint of an edge; the label records where an index
    comes from, rather than a coordinate of a tensor component.

    `index_style="graph"` displays these graph-derived labels in an abbreviated
    form. `index_style="raw"` displays the stored index expressions. The default
    `index_style="alphabet"` gives them readable letters. All three views below
    display the **same tensor**, with the same index connections. Dimensions can
    be included when the vector space or its size needs to be visible.

    `wrap_indices(scope)` gives a copy's indices a separate namespace. Primed
    labels show this separation and prevent accidental contractions between
    independent copies, such as a numerator and its conjugate.
    **A prime does not itself mean conjugation.**
    """)
    return


@app.cell
def _(DisplaySettings, S, gamma, mo, notation, p):
    _hedge = S("gammalooprs::hedge")
    _indexed = gamma(_hedge(2, 1), _hedge(3, 1), _hedge(4, 1))
    _original = p(_hedge(4, 1))
    _copy = _original.wrap_indices(S("notation::copy"))
    mo.vstack(
        [
            mo.hstack(
                [
                    mo.vstack([mo.md("**Alphabet**"), notation(_indexed)]),
                    mo.vstack(
                        [
                            mo.md("**Graph-derived labels**"),
                            notation(_indexed, DisplaySettings(index_style="graph")),
                        ]
                    ),
                    mo.vstack(
                        [
                            mo.md("**Stored index expressions**"),
                            notation(_indexed, DisplaySettings(index_style="raw")),
                        ]
                    ),
                    mo.vstack(
                        [
                            mo.md("**With dimensions**"),
                            notation(_indexed, DisplaySettings(show_dimensions=True)),
                        ]
                    ),
                ],
                wrap=True,
            ),
            mo.md(
                "**Two independent copies:** `original * original.wrap_indices(scope)`"
            ),
            _original * _copy,
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 4. Scalar products in Schoonschip and dot notation

    **Schoonschip notation** writes a vector in the position of an index
    contracted with that vector. This is the convention described in
    [A. Heck, *FORM for Pedestrians*, §1.2.2, pp. 9–10](https://www.nikhef.nl/~form/maindir/documentation/tutorial/book.pdf#page=14).

    The three live displays below start with the indexed Minkowski metric
    multiplied by two vectors. Each repeated index is summed. Next,
    `contract(metrics=False)` moves the vectors into the metric's contracted
    positions: Spenso's default display uses a bra and position markers. Finally,
    `to_dots()` writes a dot product. They represent the **same scalar product**.
    The foldout shows the other display settings for the contracted metric.

    `metrics=False` retains the metric instead of eliminating it.
    Multiplying two unresolved vectors, `p * q`, directly creates this scalar
    product too. Assigning the same explicit index requests its Einstein
    contraction; it does not select a tensor component.

    The last two examples compare a product of scalar products with a tensor
    product of vectors. The small gap separates scalar factors; the dots inside
    them retain their ordinary spacing. Different free labels leave a rank-two
    tensor. They do not form a dot.
    """)
    return


@app.cell
def _(DisplaySettings, TensorExpression, dot, lorentz, mo, notation, p, q):
    _scalar = p * q
    _explicit = (p("mu") * q("mu")).contract().to_dots()
    _metric_source = (
        TensorExpression.g(lorentz("mu"), lorentz("nu")) * p("mu") * q("nu")
    )
    _metric_product = _metric_source.contract(metrics=False)
    _dotted = _metric_product.to_dots()
    _product = (p * q) * (p * p)
    _free = p("mu") * q("nu")
    assert _scalar.rank == _explicit.rank == _metric_product.rank == _product.rank == 0
    assert _free.rank == 2
    assert _scalar == _explicit == _dotted == dot(p, q)
    assert (
        _metric_product.to_tensor().scalar() - _dotted.to_tensor().scalar()
    ).expand() == 0
    mo.vstack(
        [
            mo.hstack(
                [
                    mo.vstack(
                        [mo.md("**Indexed metric and vectors**"), _metric_source]
                    ),
                    mo.vstack(
                        [
                            mo.md("**`contract(metrics=False)`**"),
                            _metric_product,
                        ]
                    ),
                    mo.vstack([mo.md("**Dot: `metric_product.to_dots()`**"), _dotted]),
                ],
                wrap=True,
                align="start",
            ),
            mo.accordion(
                {
                    "Alternative display settings": mo.hstack(
                        [
                            mo.vstack(
                                [
                                    mo.md("**Schoonschip: vector arguments**"),
                                    notation(_metric_product, DisplaySettings.call()),
                                ]
                            ),
                            mo.vstack(
                                [
                                    mo.md("**Schoonschip: index positions**"),
                                    notation(
                                        _metric_product, DisplaySettings.schoonschip()
                                    ),
                                ]
                            ),
                        ],
                        wrap=True,
                        align="start",
                    ),
                    "Spenso-generated Typst: contracted metric": mo.md(
                        "```typst\n" + _metric_product.to_typst() + "\n```"
                    ),
                    "Spenso-generated Typst: dot product": mo.md(
                        "```typst\n" + _dotted.to_typst() + "\n```"
                    ),
                }
            ),
            mo.hstack(
                [
                    mo.vstack([mo.md("**Product of scalar products**"), _product]),
                    mo.vstack(
                        [mo.md('**Two free indices: `p("mu") * q("nu")`**'), _free]
                    ),
                ],
                wrap=True,
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 5. A metric relabels an axis

    `TensorExpression.g(lorentz("mu"), lorentz("nu"))` constructs and indexes
    the metric directly. In the displays below, its contraction with the vector
    replaces the vector's index with the metric's remaining free index.
    Both arguments must describe the same space and dimension; a space paired
    with its dual is allowed.
    """)
    return


@app.cell
def _(TensorExpression, lorentz, mo, p):
    _metric = TensorExpression.g(lorentz("mu"), lorentz("nu"))
    _source = _metric * p("nu")
    _contracted = _source.contract()
    assert _contracted == p("mu")
    mo.hstack(
        [
            mo.vstack([mo.md("**Indexed metric**"), _metric]),
            mo.vstack([mo.md("**Before `contract()`**"), _source]),
            mo.vstack([mo.md("**After `contract()`**"), _contracted]),
        ],
        wrap=True,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 6. Vectors contracted into tensor indices

    The first example contracts a spinor and a Minkowski vector into the
    rank-three tensor `A`. Compare its indexed factors with the result of
    `contract()`: Spenso moves the vectors into a bra/ket and keeps markers at
    their former index positions. Only the second bispinor axis remains free,
    so the result has **rank one**. A bra/ket here records the contraction;
    it does not apply complex conjugation.

    The examples underneath let the renderer show its marker conventions for
    an unresolved port, a self-dual bispinor space, Minkowski space, and a
    user-defined space and its dual. In the latter pair, vectors and covectors
    have opposite orientations. Each marker below is part of the tensor's
    rendered output.
    """)
    return


@app.cell
def _(
    A,
    AUTO,
    B,
    Q,
    Representation,
    Symbol,
    TensorName,
    lorentz,
    mo,
    p,
    spinor,
):
    generic_source = Symbol.I * A(1, AUTO, 1) * Q(1) * p(1)
    generic_compact = generic_source.contract()
    _bra = (p("mu") * B("mu", AUTO)).contract()
    assert generic_compact.rank == _bra.rank == 1
    _spinor_tensor = TensorName("notation::R")(spinor, lorentz)
    _spinor_compact = (_spinor_tensor("a", AUTO) * Q("a")).contract()
    _space = Representation("notation::V", 3, is_self_dual=False)
    _tensor = TensorName("notation::M")(_space.dual(), _space)
    _vector = TensorName.vector("notation::v")(_space)
    _covector = TensorName.vector("notation::w")(_space.dual())
    _ket = (_tensor(AUTO, "i") * _covector("i")).contract()
    _bra_dual = (_vector("i") * _tensor("i", AUTO)).contract()
    assert _spinor_compact.rank == _ket.rank == _bra_dual.rank == 1
    mo.vstack(
        [
            mo.hstack(
                [
                    mo.vstack(
                        [mo.md("**Explicitly indexed factors**"), generic_source]
                    ),
                    mo.vstack([mo.md("**`contract()`**"), generic_compact]),
                ],
                wrap=True,
                align="start",
            ),
            mo.hstack(
                [
                    mo.vstack([mo.md("**Unresolved ports**"), _tensor]),
                    mo.vstack([mo.md("**Bispinor contraction**"), _spinor_compact]),
                    mo.vstack([mo.md("**Minkowski contraction**"), _bra]),
                    mo.vstack([mo.md("**Vector into a dual index**"), _bra_dual]),
                    mo.vstack([mo.md("**Covector into a vector index**"), _ket]),
                ],
                wrap=True,
                align="start",
            ),
        ]
    )
    return generic_compact, generic_source


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 7. Dirac slashes use the same port conventions

    The first display multiplies a Dirac gamma matrix by a vector with a matching
    Lorentz index. `contract()` produces the slash in the next display. It is
    a matrix in bispinor space, with **two** spinor ports and no Lorentz port.
    The last display additionally supplies the spinor `Q` to its row port:
    only the column port remains open. The row marker records that contraction.

    This is the Dirac version of the previous generic rank-three example.
    The slash hides the contracted Lorentz axis; it does not hide either
    of the two spinor ports. Forming a slash does not evaluate gamma identities.
    """)
    return


@app.cell
def _(AUTO, Q, Symbol, gamma, mo, p):
    _slash_source = gamma(AUTO, AUTO, "mu") * p("mu")
    _slash = _slash_source.contract()
    joined = (Symbol.I * gamma(1, AUTO, 1) * Q(1) * p(1)).contract()
    assert _slash.rank == 2
    assert joined.rank == 1
    mo.hstack(
        [
            mo.vstack([mo.md("**Indexed gamma and vector**"), _slash_source]),
            mo.vstack([mo.md("**Slash, two open spinor ports**"), _slash]),
            mo.vstack([mo.md("**After contracting the row spinor Q**"), joined]),
        ],
        wrap=True,
    )
    return (joined,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 8. Ordered chains and traces

    The repeated middle bispinor index connects the two Dirac matrices in the
    first display. `contract()` groups their ordered product in square brackets
    and turns their vector contractions into slashes. Swapping the matrices
    generally changes the tensor. The third display additionally contracts `Q`
    into the row index and retains one open spinor port.

    The second row connects both pairs of bispinor indices of two gamma matrices.
    Spenso displays their sum as a trace, with two free Lorentz indices.
    **A matrix trace is not necessarily a scalar tensor.** `contract()` collects
    the trace; `simplify_algebra(color=False)` evaluates its Dirac identity,
    producing the metric and its coefficient in the last display.
    """)
    return


@app.cell
def _(AUTO, Q, gamma, mo, p, q):
    _open_source = gamma("a", "b", "mu") * p("mu") * gamma("b", "c", "nu") * q("nu")
    _open_compact = _open_source.contract()
    chain_source = (
        Q("a") * gamma("a", "b", "mu") * p("mu") * gamma("b", AUTO, "nu") * q("nu")
    )
    chain_compact = chain_source.contract()
    _closed_source = gamma("a", "b", "mu") * gamma("b", "a", "nu")
    _closed = _closed_source.contract()
    _evaluated = _closed.simplify_algebra(color=False)
    assert chain_compact.rank == 1
    assert _open_source.rank == _open_compact.rank == 2
    assert _closed.rank == _evaluated.rank == 2
    mo.vstack(
        [
            mo.hstack(
                [
                    mo.vstack([mo.md("**Indexed matrix product**"), _open_source]),
                    mo.vstack([mo.md("**`contract()`**"), _open_compact]),
                    mo.vstack([mo.md("**With Q contracted**"), chain_compact]),
                ],
                wrap=True,
            ),
            mo.hstack(
                [
                    mo.vstack([mo.md("**Connected spinor indices**"), _closed_source]),
                    mo.vstack([mo.md("**`contract()`**"), _closed]),
                    mo.vstack(
                        [mo.md("**`simplify_algebra(color=False)`**"), _evaluated]
                    ),
                ],
                wrap=True,
            ),
        ]
    )
    return chain_compact, chain_source


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 9. Displaying contracted vectors

    Schoonschip notation describes the contraction: a vector replaces the
    index contracted with it. Spenso's display settings choose how to print it.

    The default `ports` setting shows contracted vectors in a bra/ket, with
    markers at their former index positions. `schoonschip` places bold vector
    heads in those index positions. `call` writes the tensor with its index
    labels and contracted vectors as function arguments. These are three
    presentations of the **same tensor**, with the same remaining free index.

    Use `to_html(settings=...)` or `to_svg(settings=...)` for alternative layouts.
    `to_typst()` continues to provide pure Typst source for the ports layout.
    """)
    return


@app.cell
def _(DisplaySettings, generic_compact, mo, notation):
    mo.hstack(
        [
            mo.vstack(
                [mo.md("**Indices and bra/ket (`ports`)**"), notation(generic_compact)]
            ),
            mo.vstack(
                [
                    mo.md("**Vectors in index positions (`schoonschip`)**"),
                    notation(generic_compact, DisplaySettings.schoonschip()),
                ]
            ),
            mo.vstack(
                [
                    mo.md("**Vector arguments (`call`)**"),
                    notation(generic_compact, DisplaySettings.call()),
                ]
            ),
        ],
        wrap=True,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 10. Scalar color invariants

    Each column below shows, in order, the fundamental quadratic Casimir, the
    adjoint quadratic Casimir, and the fundamental quadratic Dynkin index.
    These are **scalar invariants**, not tensors with hidden open indices.

    Compare the compact symbols with explicit display, which writes the degree
    and representation. Dimension-aware display also identifies the representation
    size. These settings change notation without replacing the symbolic
    invariant by its numerical value.
    """)
    return


@app.cell
def _(DisplaySettings, Representation, mo, notation):
    _fundamental = Representation.cof(3)
    _adjoint = Representation.coad(8)
    _invariants = [
        _fundamental.casimir(),
        _adjoint.casimir(),
        _fundamental.dynkin_index(),
    ]
    mo.hstack(
        [
            mo.vstack(
                [mo.md("**Compact**"), *[notation(_value) for _value in _invariants]]
            ),
            mo.vstack(
                [
                    mo.md("**Explicit degree**"),
                    *[
                        notation(_value, DisplaySettings(invariant_style="explicit"))
                        for _value in _invariants
                    ],
                ]
            ),
            mo.vstack(
                [
                    mo.md("**With dimensions**"),
                    *[
                        notation(_value, DisplaySettings(show_dimensions=True))
                        for _value in _invariants
                    ],
                ]
            ),
        ],
        wrap=True,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 11. Tensor components are coordinates, not index labels

    The example constructs a tensor of shape 2 by 3, with scalar arguments `x`
    and `7`, then selects `component_tensor[0, 1]`. Both displays show this same
    component: the default uses superscripts, while array mode uses square
    brackets. These numbers are **zero-based basis coordinates**, rather than
    symbolic free indices. The scalar arguments remain in parentheses.

    The concrete tensor explorer starts in Memory grid view; its Matrix toggle
    shows component expressions. Deeper blue means a larger stored component
    payload in bytes, not a larger physical or numerical component. A missing
    sparse entry is an implicit zero.
    """)
    return


@app.cell
def _(DisplaySettings, Representation, S, TensorName, mo, notation):
    component_tensor = TensorName("notation::C")(
        S("notation::x"), 7, Representation.euc(2), Representation.euc(3)
    ).to_tensor()
    _component = component_tensor[0, 1]
    mo.vstack(
        [
            mo.hstack(
                [
                    mo.vstack(
                        [mo.md("**Superscript coordinates**"), notation(_component)]
                    ),
                    mo.vstack(
                        [
                            mo.md("**Array coordinates**"),
                            notation(
                                _component, DisplaySettings(component_style="array")
                            ),
                        ]
                    ),
                ],
                wrap=True,
            ),
            component_tensor,
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 12. Custom tensor names and pure Typst source

    A `TensorName` print mapping customizes the head. The tensor printer still
    places its arguments, ports and component coordinates. Here `macron(J)`
    produces the overline in the live display; the bar is a chosen name, not an operation
    that conjugates the tensor.

    Pure Typst output remains available: `expression.to_typst()` returns the
    math source, without an interactive HTML explorer. The foldout below exposes
    the exact source for the supplied-spinor example.
    """)
    return


@app.cell
def _(TensorName, joined, mo, spinor):
    _Jbar = TensorName(
        "notation::Jbar",
        print={"typst": "macron(J)", "latex": r"\bar{J}"},
    )(spinor)
    mo.vstack(
        [
            _Jbar("a"),
            mo.accordion(
                {
                    "Pure Typst source for the supplied-spinor example": mo.md(
                        "```typst\n" + joined.to_typst() + "\n```"
                    )
                }
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 13. Contraction changes notation, not components or tensor type

    The checks below evaluate the generic tensor and the Dirac matrix product
    before and after `contract()`. Every component agrees. The contracted expressions
    remain `TensorExpression` objects, so further tensor operations can be
    called directly.

    Scalar multiplication also preserves this type in both orders:
    `Symbol.I * gamma` and `gamma * Symbol.I`. Symbolic tensor operations return
    `TensorExpression`; operations involving concrete tensors build a
    `TensorNetwork` for component evaluation.
    """)
    return


@app.cell
def _(
    Symbol,
    TensorExpression,
    TensorLibrary,
    chain_compact,
    chain_source,
    gamma,
    generic_compact,
    generic_source,
    mo,
):
    _library = TensorLibrary.hep_lib_atom()
    for _source, _compact in (
        (generic_source, generic_compact),
        (chain_source, chain_compact),
    ):
        _before = _source.to_tensor(library=_library)
        _after = _compact.to_tensor(library=_library)
        assert _before.shape == _after.shape
        for _a, _b in zip(_before, _after, strict=True):
            assert (_a - _b).expand() == 0
    assert isinstance(Symbol.I * gamma, TensorExpression)
    assert isinstance(gamma * Symbol.I, TensorExpression)
    mo.md(
        "**Verified:** matching components, retained tensor types, and one remaining spinor axis in each supplied example."
    )
    return


if __name__ == "__main__":
    app.run()
