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

    A tensor has **axes in specified representations**, optionally with index
    labels. Its display can keep those labels, expose unresolved ports, or show
    vectors already contracted into its axes. This notebook uses the actual
    Spenso printers throughout; the examples remain editable tensor expressions.

    **Reading key:** a letter is an index label; a hollow **□** is an unresolved
    port; a filled **●**, **■**, **◀︎** or **▶︎** marks a port supplied with a vector.
    Supplied ports are contracted and do not count towards the result's rank.
    A slash denotes a Dirac contraction, brackets an ordered matrix word, and
    `Tr` a closed matrix word.

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
    ## 1. Open ports and index labels

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
    ## 2. Alphabets, graph identities and scoped copies

    The default alphabet display gives generated Lorentz indices Greek letters
    and bispinor indices Latin letters. These letters label axes; their position
    follows the representation's display convention. In particular, the two
    lower bispinor labels on a gamma matrix mean its row and column axes.

    Graph mode exposes graph identities, and raw mode exposes the stored labels.
    All three views below display the **same tensor**. Dimensions can be included
    when the representation or its size needs to be visible.

    Primed labels distinguish indices wrapped in a separate scope. Scoping
    prevents accidental contraction between independent copies, such as a
    numerator and its conjugate; **a prime does not itself mean conjugation**.
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
                            mo.md("**Graph**"),
                            notation(_indexed, DisplaySettings(index_style="graph")),
                        ]
                    ),
                    mo.vstack(
                        [
                            mo.md("**Raw**"),
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
    ## 3. Dots and ordinary multiplication

    A dot is a **complete scalar contraction**. In Minkowski space,
    $p\cdot q=g_{\mu\nu}p^\mu q^\nu$. Multiplying two unresolved vectors has
    an unambiguous pairing and creates this scalar product. Assigning the same
    explicit index requests the corresponding Einstein contraction.

    Two dots next to one another mean a product of scalars,
    $(p\cdot q)(p\cdot p)$. The small gap separates the factors; the dots
    inside them retain their ordinary spacing. Different free labels, as in
    $p^\mu q^\nu$, leave a rank-two tensor. They do not form a dot.
    """)
    return


@app.cell
def _(dot, mo, p, q):
    _scalar = p * q
    _explicit = (p("mu") * q("mu")).contract().to_dots()
    _product = (p * q) * (p * p)
    _free = p("mu") * q("nu")
    assert _scalar.rank == _explicit.rank == _product.rank == 0
    assert _free.rank == 2
    assert _scalar == dot(p, q)
    mo.hstack(
        [
            mo.vstack([mo.md("**`p * q`**"), _scalar]),
            mo.vstack(
                [mo.md('**`(p("mu") * q("mu")).contract().to_dots()`**'), _explicit]
            ),
            mo.vstack([mo.md("**`(p * q) * (p * p)`**"), _product]),
            mo.vstack([mo.md('**`p("mu") * q("nu")`**'), _free]),
        ],
        wrap=True,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 4. A metric relabels an axis

    `TensorExpression.g(lorentz("mu"), lorentz("nu"))` constructs and indexes
    the metric directly. Contracting it with $p^\nu$ gives $p^\mu$ in the
    representation's contraction convention. Both arguments must describe the
    same space and dimension; a space paired with its dual is allowed.
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
    ## 5. Supplied ports: filled shapes and bra/ket notation

    Contracting a vector into an axis can be written as a tensor with that vector
    in its argument, or as a bra/ket around the tensor. For example,
    $\langle p|B$ represents $\sum_\mu p_\mu B^{\mu\nu}$, with only $\nu$ open.
    A bra/ket here records which vectors supply ports; it does not by itself
    apply complex conjugation.

    A filled marker keeps the occupied port's position visible:

    | Marker | Supplied axis |
    | --- | --- |
    | **●** | A self-dual representation, including the bispinor space |
    | **■** | A representation with an inline metric, including Minkowski space |
    | **◀︎** | The ket orientation of a dualizable representation |
    | **▶︎** | The bra orientation of a dualizable representation |
    | **□** | An unresolved external port; no vector has been supplied |

    The rank-three example below represents
    $i\sum_{a,\mu}Q_a A_{ab}{}^\mu p_\mu$. After contraction its **rank is one**:
    only the second bispinor port survives. The filled circle and square
    account for the two supplied axes; the hollow square is the remaining axis.
    """)
    return


@app.cell
def _(A, AUTO, B, Q, Symbol, mo, p):
    generic_source = Symbol.I * A(1, AUTO, 1) * Q(1) * p(1)
    generic_compact = generic_source.contract()
    _bra = (p("mu") * B("mu", AUTO)).contract()
    assert generic_compact.rank == _bra.rank == 1
    mo.hstack(
        [
            mo.vstack([mo.md("**Generic tensor, indexed**"), generic_source]),
            mo.vstack([mo.md("**Generic tensor, supplied ports**"), generic_compact]),
            mo.vstack([mo.md("**One supplied Minkowski vector**"), _bra]),
        ],
        wrap=True,
        align="start",
    )
    return generic_compact, generic_source


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 6. Dirac slashes use the same port conventions

    A slash means $\not p=\gamma^\mu p_\mu$, a matrix in bispinor space.
    The slash has **two** spinor ports and no Lorentz port. Supplying $Q$ to
    its row port gives $\langle Q|\not p$: only the column port remains open.
    Its filled circle means that the row spinor has already been supplied.

    This is the Dirac version of the previous generic rank-three example.
    The slash hides the contracted Lorentz axis; it does not hide either
    of the two spinor ports. Forming a slash does not evaluate gamma identities.
    """)
    return


@app.cell
def _(AUTO, Q, Symbol, gamma, mo, p):
    _slash = (gamma(AUTO, AUTO, "mu") * p("mu")).contract()
    joined = (Symbol.I * gamma(1, AUTO, 1) * Q(1) * p(1)).contract()
    assert _slash.rank == 2
    assert joined.rank == 1
    mo.hstack(
        [
            mo.vstack([mo.md("**Slash, two open spinor ports**"), _slash]),
            mo.vstack([mo.md("**Slash with a supplied row spinor**"), joined]),
        ],
        wrap=True,
    )
    return (joined,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 7. Ordered chains and traces

    $[\not p\,\not q]_{ac}$ is the matrix word
    $\sum_b(\not p)_{ab}(\not q)_{bc}$. Square brackets group the ordered
    factors. Swapping them generally changes the tensor. With a supplied row
    spinor, $\langle Q|[\not p\,\not q]$ retains one open spinor port.

    `Tr` closes the matrix channel: $\operatorname{Tr}(\gamma^\mu\gamma^\nu)$
    sums both spinor indices, while the two Lorentz indices remain free.
    **A matrix trace is not necessarily a scalar tensor.** `contract()` collects
    the trace; `simplify_algebra(color=False)` evaluates its Dirac identity,
    here $\operatorname{Tr}(\gamma^\mu\gamma^\nu)=4g^{\mu\nu}$.
    """)
    return


@app.cell
def _(AUTO, Q, gamma, mo, p, q):
    chain_source = (
        Q("a") * gamma("a", "b", "mu") * p("mu") * gamma("b", AUTO, "nu") * q("nu")
    )
    chain_compact = chain_source.contract()
    _closed = (gamma("a", "b", "mu") * gamma("b", "a", "nu")).contract()
    _evaluated = _closed.simplify_algebra(color=False)
    assert chain_compact.rank == 1
    assert _closed.rank == _evaluated.rank == 2
    mo.hstack(
        [
            mo.vstack([mo.md("**Supplied open chain**"), chain_compact]),
            mo.vstack([mo.md("**Collected trace**"), _closed]),
            mo.vstack([mo.md("**Evaluated trace**"), _evaluated]),
        ],
        wrap=True,
    )
    return chain_compact, chain_source


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 8. Ports, Schoonschip and call layouts

    Ports layout is the default: supplied vectors appear in a bra/ket and
    occupied axes have filled markers. Schoonschip layout writes supplied
    vectors compactly as tensor arguments. Call layout shows the head and its
    arguments explicitly. These are three presentations of the **same
    contracted tensor**, with the same remaining external axis.

    Use `to_html(settings=...)` or `to_svg(settings=...)` for alternative layouts.
    `to_typst()` continues to provide pure Typst source for the ports layout.
    """)
    return


@app.cell
def _(DisplaySettings, generic_compact, mo, notation):
    mo.hstack(
        [
            mo.vstack([mo.md("**Ports**"), notation(generic_compact)]),
            mo.vstack(
                [
                    mo.md("**Schoonschip**"),
                    notation(generic_compact, DisplaySettings.schoonschip()),
                ]
            ),
            mo.vstack(
                [mo.md("**Call**"), notation(generic_compact, DisplaySettings.call())]
            ),
        ],
        wrap=True,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 9. Scalar color invariants

    The compact symbols $C_F$, $C_A$ and $T_R$ are respectively the quadratic
    Casimirs of the fundamental and adjoint representations, and the quadratic
    Dynkin index of the fundamental representation. They are **scalar
    invariants**, not tensors with hidden open indices.

    Explicit display writes the degree and representation: $C_2(F)$, $C_2(A)$
    and $I_2(F)$. Dimension-aware display also identifies the representation
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
    ## 10. Tensor components are coordinates, not index labels

    For a tensor with shape $2\times3$, $C(x,7)^{0,1}$ means its component at
    coordinates $(0,1)$. Array mode writes the same component as $C(x,7)[0,1]$.
    The superscripts here are **zero-based basis coordinates**, rather than
    symbolic free indices. Scalar arguments $x,7$ remain in parentheses.

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
    ## 11. Custom tensor names and pure Typst source

    A `TensorName` print mapping customizes the head. The tensor printer still
    places its arguments, ports and component coordinates. Here `macron(J)`
    produces the overline on $\bar J$; the bar is a chosen name, not an operation
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
    ## 12. Contraction changes notation, not components or tensor type

    The checks below evaluate the generic tensor and the Dirac word before and
    after `contract()`. Every component agrees. The contracted expressions
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
