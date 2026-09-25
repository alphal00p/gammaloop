# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "symbolica==3.0.0",
#     "marimo==0.24.0",
#     "typst==0.15.0",
# ]
# ///

import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium", app_title="Idenso: gamma simplification")


@app.cell
def _():
    from functools import partial
    from math import prod
    from time import perf_counter

    import marimo as mo
    from symbolica import E, S
    from symbolica.community.spenso import (
        AUTO,
        DisplaySettings,
        GammaChainOrdering,
        GammaSimplifySettings,
        Representation,
        TensorExpression,
        TensorName,
        chain,
        trace,
    )

    return (
        AUTO,
        DisplaySettings,
        E,
        GammaChainOrdering,
        GammaSimplifySettings,
        Representation,
        S,
        TensorExpression,
        TensorName,
        chain,
        mo,
        partial,
        perf_counter,
        prod,
        trace,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # How capable is gamma simplification?

    Idenso's `TensorExpression.simplify_gamma()` collects spinor contractions,
    simplifies open Dirac chains and evaluates closed traces. It repeatedly
    applies chain, metric and epsilon rewrites until the expression stops
    changing. Every example below runs the installed implementation; the
    identity cards check the result against an explicit expected expression
    and check that a second pass does nothing.

    | Capability | Scope |
    |:--|:--|
    | Clifford contractions and ordinary traces | Compatible symbolic Lorentz dimensions |
    | Slashed momenta and chain joining | Registered Spenso tensor notation |
    | Chisholm, gamma5, chiral projectors, charge conjugation | Four dimensions |
    | Canonical chain ordering | Opt-in; may produce more terms |
    | Three-gamma epsilon expansion | Opt-in for open chains; automatic in short 4D trace kernels |
    | Short trace kernels | Ordinary and gamma5 traces, 1–14 ordinary gammas |
    | Expanded trace output | Opt-in with `expand_traces=True`; scalar spectators stay factored |

    This is a partial algebraic normal form. It does not impose on-shell
    kinematics, external Dirac equations or a dimensional-regularization
    gamma5 scheme. Run with the **combined Symbolica community host** containing
    Spenso and Idenso, plus Marimo and Typst; the ordinary Symbolica package
    alone does not provide these community bindings.
    """)
    return


@app.cell
def _(AUTO, Representation, S, TensorExpression, chain, partial, trace):
    spin = Representation.bis(4)
    lorentz = Representation.mink(4)
    a, b, c, mu, nu, rho, sigma, D = S(
        *(
            f"gamma_tutorial::{name}"
            for name in ("a", "b", "c", "mu", "nu", "rho", "sigma", "D")
        )
    )
    gamma = TensorExpression.gamma(4)
    gm, gn, gr, gs = (
        gamma(AUTO, AUTO, lorentz(index)) for index in (mu, nu, rho, sigma)
    )
    g5 = TensorExpression.gamma5(4)(AUTO, AUTO)
    plus = TensorExpression.projp(4)(AUTO, AUTO)
    minus = TensorExpression.projm(4)(AUTO, AUTO)
    word = partial(chain, spin(a), spin(b))
    tr = partial(trace, spin)
    identity = spin.g(a, b)
    return (
        D,
        a,
        b,
        c,
        g5,
        gamma,
        gm,
        gn,
        gr,
        gs,
        identity,
        lorentz,
        minus,
        mu,
        nu,
        plus,
        rho,
        sigma,
        spin,
        tr,
        word,
    )


@app.cell(hide_code=True)
def _(mo):
    show_dimensions = mo.ui.switch(value=False, label="Show representation dimensions")
    show_dimensions
    return (show_dimensions,)


@app.cell
def _(DisplaySettings, E, TensorExpression, mo, show_dimensions):
    display_settings = DisplaySettings(show_dimensions=show_dimensions.value)

    def identity_card(title, expression, expected, settings=None):
        result = expression.simplify_gamma(settings)
        expected_atom = (
            expected.to_expression()
            if isinstance(expected, TensorExpression)
            else expected
        )
        # Expand only the small verification difference, not the input numerator.
        assert (result.to_expression() - expected_atom).expand() == E("0"), title
        assert (
            result.simplify_gamma(settings).to_expression() == result.to_expression()
        ), title
        return mo.vstack(
            [
                mo.md(f"**{title}** — identity and fixed point checked"),
                mo.hstack(
                    [
                        mo.vstack(
                            [
                                mo.md("Input"),
                                mo.Html(expression.to_html(settings=display_settings)),
                            ]
                        ),
                        mo.vstack(
                            [
                                mo.md("Result"),
                                mo.Html(result.to_html(settings=display_settings)),
                            ]
                        ),
                    ],
                    widths="equal",
                    align="start",
                ),
            ]
        )

    return display_settings, identity_card


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 1. Open chains and contracted indices

    Gamma's argument order is `(spinor-in, spinor-out, Lorentz)`.
    `AUTO` supplies local ports to the public `chain` and `trace` builders;
    `word` below simply binds the two external spinor indices with `partial`.
    Repeating a Lorentz label contracts it. In four dimensions,
    $\gamma^\mu\gamma_\mu=4\mathbf{1}$ and
    $\gamma^\mu\gamma^\nu\gamma_\mu=-2\gamma^\nu$.

    Explicit products are collected too. For prescribed sewing indices, form
    the Symbolica product explicitly and wrap it in `TensorExpression`.
    """)
    return


@app.cell
def _(
    TensorExpression,
    a,
    b,
    c,
    gamma,
    gm,
    gn,
    gr,
    gs,
    identity,
    identity_card,
    mo,
    mu,
    word,
):
    joined = TensorExpression(
        gamma(a, c, mu).to_expression() * gamma(c, b, mu).to_expression()
    )
    open_checks = mo.vstack(
        [
            identity_card("Adjacent contraction", word(gm, gm), 4 * identity),
            identity_card("Sandwich contraction", word(gm, gn, gm), -2 * word(gn)),
            identity_card(
                "Four-dimensional Chisholm reversal",
                word(gm, gn, gr, gs, gm),
                -2 * word(gs, gr, gn),
            ),
            identity_card("Collect explicitly sewn matrices", joined, 4 * identity),
        ],
        gap=2,
    )
    open_checks
    return (open_checks,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 2. Symbolic dimension and trace normalization

    Change **Lorentz** dimension before simplifying a contracted expression.
    `with_lorentz_dimension(D)` leaves the spinor dimension unchanged, so
    $\gamma^\mu\gamma_\mu=D\mathbf{1}$ while $\mathrm{tr}(\mathbf{1})=4$.
    If you instead explicitly use `bis(D)`, the empty trace is $D$.
    The trace normalization comes from the spin representation; it is not
    inferred as $2^{D/2}$.
    """)
    return


@app.cell
def _(D, Representation, gm, gn, identity, identity_card, mo, mu, nu, tr, trace, word):
    dimension_checks = mo.vstack(
        [
            identity_card(
                "D-dimensional contraction",
                word(gm, gm).with_lorentz_dimension(D),
                D * identity,
            ),
            identity_card(
                "D-dimensional sandwich",
                word(gm, gn, gm).with_lorentz_dimension(D),
                (2 - D) * word(gn).with_lorentz_dimension(D),
            ),
            identity_card(
                "Lorentz D, spinor 4",
                tr(gm, gn).with_lorentz_dimension(D),
                4 * Representation.mink(D).g(mu, nu),
            ),
            identity_card(
                "Explicit bis(D) trace normalization", trace(Representation.bis(D)), D
            ),
        ],
        gap=2,
    )
    dimension_checks
    return (dimension_checks,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 3. Slashed momenta

    Compact momenta inside a gamma factor represent a slash. This example
    constructs that registered Symbolica form directly:
    `gamma(bis(4,a), bis(4,b), p(mink(4)))`.
    No mass-shell relation is assumed: $p^2$ remains a scalar product.
    """)
    return


@app.cell
def _(
    S,
    TensorExpression,
    TensorName,
    a,
    b,
    identity,
    identity_card,
    lorentz,
    mo,
    spin,
    word,
):
    p, q = S("gamma_tutorial::p", "gamma_tutorial::q")
    gamma_head = TensorName.gamma().to_expression()
    metric_head = TensorName.g().to_expression()

    def slash(momentum):
        return TensorExpression(
            gamma_head(spin.to_expression(), spin.to_expression(), momentum)
        )

    p_compact, q_compact = p(lorentz.to_expression()), q(lorentz.to_expression())
    slash_p = slash(p_compact)(a, b)
    slash_q = slash(q_compact)(a, b)
    slash_sandwich = word(slash(p_compact), slash(q_compact), slash(p_compact))
    incoming, outgoing, chain_head = S("spenso::in", "spenso::out", "spenso::chain")
    p_squared = metric_head(p_compact, p_compact)
    p_dot_q = metric_head(p_compact, q_compact)
    slash_check = identity_card(
        "Slash sandwich: 2(p·q) slash(p) − p² slash(q)",
        slash_sandwich,
        2 * p_dot_q * word(slash_p) - p_squared * word(slash_q),
    )
    slash_check
    return chain_head, gamma_head, incoming, outgoing, slash, slash_check


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 4. Ordinary traces

    Odd ordinary traces vanish. Even traces recurse to pairwise metrics;
    four gammas give the familiar three signed pairings. Disable trace
    evaluation when retaining closed spin chains is preferable downstream.

    `GammaSimplifySettings(expand_traces=True)` expands each evaluated trace
    body while keeping surrounding scalar factors factored. The default is
    `False`; `evaluate_traces=False` also disables this expansion. The example
    below expands only the standalone body when constructing its reference.
    """)
    return


@app.cell
def _(
    GammaSimplifySettings,
    gm,
    gn,
    gr,
    gs,
    identity_card,
    lorentz,
    mo,
    mu,
    nu,
    rho,
    sigma,
    tr,
):
    trace_four_expected = 4 * (
        lorentz.g(mu, nu).to_expression() * lorentz.g(rho, sigma).to_expression()
        - lorentz.g(mu, rho).to_expression() * lorentz.g(nu, sigma).to_expression()
        + lorentz.g(mu, sigma).to_expression() * lorentz.g(nu, rho).to_expression()
    )
    trace_checks = mo.vstack(
        [
            identity_card("Empty trace", tr(), 4),
            identity_card("Odd trace", tr(gm, gn, gr), 0),
            identity_card("Four-gamma trace", tr(gm, gn, gr, gs), trace_four_expected),
            identity_card(
                "Keep a trace inert",
                tr(gm, gn),
                tr(gm, gn),
                GammaSimplifySettings(evaluate_traces=False),
            ),
        ],
        gap=2,
    )
    trace_checks
    return (trace_checks,)


@app.cell
def _(
    AUTO,
    GammaSimplifySettings,
    S,
    TensorExpression,
    a,
    b,
    c,
    display_settings,
    gamma,
    lorentz,
    mo,
    mu,
    nu,
    rho,
    sigma,
    tr,
):
    # Expand only the standalone trace body, never the scalar numerator.
    _trace = tr(
        *(gamma(AUTO, AUTO, lorentz(_i)) for _i in (mu, a, b, mu, c, nu, rho, sigma))
    )
    _x, _y = S("expanded_trace_example::x", "expanded_trace_example::y")
    _spectator = (_x + _y) ** 8
    _settings = GammaSimplifySettings(expand_traces=True)
    _expected_body = _trace.simplify_gamma().to_expression().expand()
    _source = TensorExpression(_spectator * _trace.to_expression())
    _expanded = _source.simplify_gamma(_settings)
    assert _expanded.to_expression() == _spectator * _expected_body
    assert (
        _expanded.simplify_gamma(_settings).to_expression() == _expanded.to_expression()
    )
    expanded_trace_check = mo.vstack(
        [
            mo.md(
                "**Expanded trace body; scalar $(x+y)^8$ stays factored.** Exact result and rerun checked."
            ),
            mo.Html(_expanded.to_html(settings=display_settings)),
        ]
    )
    expanded_trace_check
    return (expanded_trace_check,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 5. Conservative versus canonical ordering

    The default `repeated_pairs()` strategy moves gammas toward repeated
    partners and avoids ordering every distinct factor. `canonical()` also
    applies adjacent Clifford swaps. It is useful for cancellations between
    different orders, but can increase intermediate expression size.

    The residual below is the Clifford anticommutator minus its expected
    value. Canonical mode reduces it to zero; the default leaves a residual.
    Neither option promises a globally minimal expression or a complete Fierz basis.
    """)
    return


@app.cell
def _(GammaChainOrdering, mo):
    ordering = mo.ui.radio(
        {
            "Repeated pairs (default)": GammaChainOrdering.RepeatedPairs,
            "Canonical": GammaChainOrdering.Canonical,
        },
        value="Repeated pairs (default)",
        inline=True,
        label="Open-chain ordering",
    )
    ordering
    return (ordering,)


@app.cell
def _(
    E,
    GammaSimplifySettings,
    display_settings,
    gm,
    gn,
    identity,
    lorentz,
    mo,
    mu,
    nu,
    ordering,
    word,
    TensorExpression,
):
    # Use an explicit Symbolica product for the disconnected metric factors.
    anticommutator = TensorExpression(
        word(gm, gn).to_expression()
        + word(gn, gm).to_expression()
        - 2 * lorentz.g(mu, nu).to_expression() * identity.to_expression()
    )
    assert anticommutator.simplify_gamma().to_expression() != E("0")
    assert anticommutator.simplify_gamma(
        GammaSimplifySettings.canonical()
    ).to_expression() == E("0")
    ordering_result = anticommutator.simplify_gamma(
        GammaSimplifySettings(chain_ordering=ordering.value)
    )
    mo.vstack(
        [
            mo.Html(anticommutator.to_html(settings=display_settings)),
            mo.md("**Residual after the selected pass**"),
            mo.Html(ordering_result.to_html(settings=display_settings)),
        ]
    )
    return (ordering_result,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 6. Four-dimensional gamma5 and projectors

    Gamma5 anticommutes through ordinary 4D gammas and squares to one;
    $P_\pm=(1\pm\gamma_5)/2$ obey projector identities. The internal epsilon
    normalization is fixed by
    $\mathrm{tr}(\gamma_5\gamma^\mu\gamma^\nu\gamma^\rho\gamma^\sigma)
    =4\,\epsilon(\mu,\nu,\rho,\sigma)$, **without an extra factor of i**.
    Convert conventions explicitly when comparing to other packages.
    """)
    return


@app.cell
def _(
    S,
    TensorExpression,
    g5,
    gm,
    gn,
    gr,
    gs,
    identity,
    identity_card,
    lorentz,
    minus,
    mo,
    mu,
    nu,
    plus,
    rho,
    sigma,
    tr,
    word,
):
    axial_expected = TensorExpression(
        4
        * S("spenso::epsilon")(
            *(lorentz(index).to_expression() for index in (mu, nu, rho, sigma))
        )
    )
    chiral_checks = mo.vstack(
        [
            identity_card("Gamma5 square", word(g5, g5), identity),
            identity_card("Gamma5 anticommutation", word(g5, gm), -word(gm, g5)),
            identity_card(
                "Axial four-gamma trace", tr(g5, gm, gn, gr, gs), axial_expected
            ),
            identity_card("Orthogonal projectors", word(plus, minus), 0),
            identity_card("Projector trace", tr(plus), 2),
        ],
        gap=2,
    )
    chiral_checks
    return (chiral_checks,)


@app.cell
def _(GammaSimplifySettings, display_settings, gm, gn, gr, mo, word):
    epsilon_settings = GammaSimplifySettings(expand_three_gamma_epsilon=True)
    epsilon_input = word(gm, gn, gr)
    epsilon_result = epsilon_input.simplify_gamma(epsilon_settings)
    assert (
        epsilon_input.simplify_gamma().to_expression() == epsilon_input.to_expression()
    )
    assert epsilon_result.to_expression() != epsilon_input.to_expression()
    assert (
        epsilon_result.simplify_gamma(epsilon_settings).to_expression()
        == epsilon_result.to_expression()
    )
    mo.vstack(
        [
            mo.md(
                "**Optional three-gamma expansion** — disabled by default. This rewrites the word into metric terms and a gamma5–epsilon term; more terms need not mean a better result."
            ),
            mo.Html(epsilon_result.to_html(settings=display_settings)),
        ]
    )
    return (epsilon_result,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 7. Charge conjugation and deliberate stopping points

    In the registered Weyl convention $C=-i\gamma^2\gamma^0$,
    $C^2=-\mathbf{1}$ and $C\gamma^\mu C=(\gamma^\mu)^T$.
    Charge conjugation can produce transposed words; ordinary Clifford rules
    do not treat those words as forward gammas. Mixed dimensions and
    D-dimensional gamma5 also remain explicit. The cards verify that behavior.
    """)
    return


@app.cell
def _(
    D,
    S,
    TensorExpression,
    a,
    b,
    chain_head,
    g5,
    gamma_head,
    gm,
    gn,
    identity,
    identity_card,
    incoming,
    lorentz,
    mo,
    mu,
    outgoing,
    spin,
    word,
):
    _left, _right = spin(a).to_expression(), spin(b).to_expression()
    _c = S("spenso::charge_conjugation")(incoming, outgoing)
    _forward = gamma_head(incoming, outgoing, lorentz(mu).to_expression())
    _transpose = gamma_head(outgoing, incoming, lorentz(mu).to_expression())
    _square = TensorExpression(chain_head(_left, _right, _c, _c))
    _sandwich = TensorExpression(chain_head(_left, _right, _c, _forward, _c))
    _transposed_word = TensorExpression(chain_head(_left, _right, _transpose, _forward))
    _mixed = word(gm, gn.with_lorentz_dimension(D), gm)
    _dimensional_g5 = word(g5, gm).with_lorentz_dimension(D)
    boundary_checks = mo.vstack(
        [
            identity_card("Charge-conjugation square", _square, -identity),
            identity_card(
                "Charge-conjugation sandwich",
                _sandwich,
                TensorExpression(chain_head(_left, _right, _transpose)),
            ),
            identity_card(
                "Transposed gamma word stays explicit",
                _transposed_word,
                _transposed_word,
            ),
            identity_card("Mixed Lorentz dimensions stay explicit", _mixed, _mixed),
            identity_card(
                "No D-dimensional gamma5 prescription", _dimensional_g5, _dimensional_g5
            ),
        ],
        gap=2,
    )
    boundary_checks
    return (boundary_checks,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 8. Trace growth: measure a small example

    Generic-dimensional pairing recursion produces $(2n-1)!!$ terms for
    $2n$ distinct Lorentz indices. Four-dimensional traces now use generated
    kernels through length 14, derived from the three-gamma epsilon identity.
    At lengths 10, 12 and 14 they give **693, 4,383 and 26,931** metric terms,
    matching FORM 5.0.0 `trace4` for these inputs. Repeated indices and momenta
    can simplify earlier through the shared chain contractions.

    Each even arity has a lazily generated ordinary and gamma5 table; odd
    traces return zero. The first call includes initialization if this table
    has not already been used in the session. Warm timing excludes that cost.

    Change the length to time the installed simplifier. The measurement excludes
    rendering and counts expanded terms only for this small diagnostic. It is
    a local observation, not a large-amplitude performance guarantee.
    """)
    return


@app.cell
def _(mo):
    trace_length = mo.ui.slider(
        start=2,
        stop=10,
        step=2,
        value=6,
        label="Distinct gamma factors",
        show_value=True,
    )
    trace_length
    return (trace_length,)


@app.cell
def _(AUTO, S, gamma, lorentz, mo, perf_counter, prod, tr, trace_length):
    _length = trace_length.value
    _input = tr(
        *(
            gamma(AUTO, AUTO, lorentz(S(f"gamma_tutorial::ell_{index}")))
            for index in range(_length)
        )
    )
    _start = perf_counter()
    _result = _input.simplify_gamma()
    trace_seconds = perf_counter() - _start
    trace_terms = len(list(_result.to_expression().expand().terms()))
    assert trace_terms == (1, 3, 15, 105, 693, 4383, 26931)[_length // 2 - 1]
    _samples = []
    for _ in range(3):
        _start = perf_counter()
        _input.simplify_gamma()
        _samples.append(perf_counter() - _start)
    mo.ui.table(
        [
            {
                "gamma factors": _length,
                "expanded metric terms": trace_terms,
                "generic pairing terms": prod(range(1, _length, 2)),
                "first call (ms)": round(1000 * trace_seconds, 3),
                "warm median (ms)": round(1000 * sorted(_samples)[1], 3),
            }
        ],
        selection=None,
        pagination=False,
        show_download=False,
    )
    return trace_seconds, trace_terms


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 9. Schoonschip notation changes the work, not just the display

    Contracting $p_\mu\gamma^\mu$ into $\not p$ exposes repeated momenta
    before the trace generates metric pairings. Compare **the same scalar**
    through three routes:

    | Route | Ordered stages inside the timer |
    |:--|:--|
    | Trace first | Free-index trace → attach momenta → `schoonschip_net()` |
    | Contract first | Indexed input → `schoonschip_net()` → `simplify_gamma()` |
    | Compact input | Constructed slash trace → `simplify_gamma()` |

    The first two start from the same indexed problem; their ratio includes
    the early contraction cost. Compact input starts after that conversion,
    so its time is reported separately. The gamma engine is shared. The first
    difference is whether it sees distinct indices or repeated momenta.
    `DisplaySettings(tensor_layout="schoonschip")` only changes rendering and
    cannot produce this speedup. Here the network pass changes the expression.

    Paired slashes give
    $\mathrm{tr}(\not p\not p\not q^{\,2m-2})=4p^2(q^2)^{m-1}$.
    For alternating slashes, the independent oracle is
    $T_m=2(p\cdot q)T_{m-1}-p^2q^2T_{m-2}$ with
    $T_0=4$ and $T_1=4p\cdot q$.
    These relations do not require on-shell momenta.
    """)
    return


@app.cell
def _(mo):
    slash_length = mo.ui.dropdown([4, 6, 8, 10], value=8, label="Gamma factors")
    slash_pattern = mo.ui.radio(
        {"Adjacent pairs": "paired", "Alternating p, q": "alternating"},
        value="Adjacent pairs",
        inline=True,
        label="Momentum pattern",
    )
    mo.hstack([slash_length, slash_pattern], justify="start")
    return slash_length, slash_pattern


@app.cell
def _(
    AUTO,
    E,
    S,
    TensorExpression,
    TensorName,
    gamma,
    lorentz,
    slash,
    slash_length,
    slash_pattern,
    tr,
):
    from functools import reduce
    from operator import mul

    _metric = TensorName.g().to_expression()
    _p, _q = (
        TensorName.vector(f"gamma_benchmark::{name}").to_expression()
        for name in ("p", "q")
    )
    _n = slash_length.value
    _m = _n // 2
    momentum_names = (
        ["p", "p"] + ["q"] * (_n - 2)
        if slash_pattern.value == "paired"
        else ["p", "q"] * _m
    )
    _momenta = [{"p": _p, "q": _q}[name] for name in momentum_names]
    _indices = [lorentz(S(f"gamma_benchmark::mu{i}")) for i in range(_n)]
    bare_trace = tr(*(gamma(AUTO, AUTO, i) for i in _indices))
    momentum_factors = reduce(
        mul, (p(i.to_expression()) for p, i in zip(_momenta, _indices))
    )
    indexed_trace = TensorExpression(momentum_factors * bare_trace.to_expression())
    compact_trace = tr(*(slash(p(lorentz.to_expression())) for p in _momenta))

    _pp, _qq, _pq = S(
        "gamma_benchmark::pp", "gamma_benchmark::qq", "gamma_benchmark::pq"
    )
    if slash_pattern.value == "paired":
        scalar_oracle = 4 * _pp * _qq ** (_m - 1)
    else:
        _previous, scalar_oracle = E("4"), 4 * _pq
        for _ in range(2, _m + 1):
            _previous, scalar_oracle = (
                scalar_oracle,
                (2 * _pq * scalar_oracle - _pp * _qq * _previous).expand(),
            )
    _p_compact, _q_compact = _p(lorentz.to_expression()), _q(lorentz.to_expression())
    scalar_expected = (
        scalar_oracle.replace(_pp, _metric(_p_compact, _p_compact))
        .replace(_qq, _metric(_q_compact, _q_compact))
        .replace(_pq, _metric(_p_compact, _q_compact))
    )
    return (
        bare_trace,
        compact_trace,
        indexed_trace,
        momentum_factors,
        momentum_names,
        scalar_expected,
        scalar_oracle,
    )


@app.cell
def _(compact_trace, display_settings, indexed_trace, mo):
    contracted_trace = indexed_trace.schoonschip_net()
    assert contracted_trace.to_expression() == compact_trace.to_expression()
    mo.vstack(
        [
            mo.md(
                "**Same input, after contracting the momentum indices** — exact compact expression checked."
            ),
            mo.Html(indexed_trace.to_html(settings=display_settings)),
            mo.Html(contracted_trace.to_html(settings=display_settings)),
        ]
    )
    return (contracted_trace,)


@app.cell
def _(
    E,
    TensorExpression,
    bare_trace,
    compact_trace,
    contracted_trace,
    indexed_trace,
    momentum_factors,
    perf_counter,
    scalar_expected,
):
    from statistics import median

    # Inputs are built outside the timed region. Every timed route computes
    # the complete scalar result; no precomputed trace is reused in a route.
    _routes = {
        "Trace first": lambda: TensorExpression(
            momentum_factors * bare_trace.simplify_gamma().to_expression()
        ).schoonschip_net(),
        "Contract first": lambda: indexed_trace.schoonschip_net().simplify_gamma(),
        "Compact input": lambda: compact_trace.simplify_gamma(),
    }
    _samples = {name: [] for name in _routes}
    benchmark_results = {}
    for _name, _operation in _routes.items():
        _result = _operation()  # one untimed warm-up per route
        assert (_result.to_expression() - scalar_expected).expand() == E("0"), _name
    _names = list(_routes)
    for _round in range(3):
        # Rotate order to reduce systematic first/last-route bias.
        for _name in _names[_round:] + _names[:_round]:
            _start = perf_counter()
            _result = _routes[_name]()
            _samples[_name].append(perf_counter() - _start)
            assert (_result.to_expression() - scalar_expected).expand() == E("0"), _name
            benchmark_results[_name] = _result
    benchmark_seconds = {name: median(samples) for name, samples in _samples.items()}
    benchmark_speedup = (
        benchmark_seconds["Trace first"] / benchmark_seconds["Contract first"]
    )
    return benchmark_results, benchmark_seconds, benchmark_speedup, median


@app.cell
def _(
    bare_trace,
    benchmark_results,
    benchmark_seconds,
    benchmark_speedup,
    display_settings,
    mo,
    prod,
    slash_length,
):
    metric_pairings = len(
        list(bare_trace.simplify_gamma().to_expression().expand().terms())
    )
    assert (
        metric_pairings
        == (1, 3, 15, 105, 693, 4383, 26931)[slash_length.value // 2 - 1]
    )
    mo.vstack(
        [
            mo.ui.table(
                [
                    {
                        "route": name,
                        "median (ms)": round(1000 * seconds, 3),
                        "final scalar terms": len(
                            list(
                                benchmark_results[name].to_expression().expand().terms()
                            )
                        ),
                    }
                    for name, seconds in benchmark_seconds.items()
                ],
                selection=None,
                pagination=False,
                show_download=False,
            ),
            mo.md(
                f"**Measured early-contraction speedup: {benchmark_speedup:.2f}×.** The trace-first route builds **{metric_pairings}** free-index metric terms before seeing the repeated momenta. All three final scalar results agree exactly. Timings are medians of three warm runs, excluding construction, assertions and rendering; they depend on the installed build and hardware."
            ),
            mo.Html(
                benchmark_results["Contract first"].to_html(settings=display_settings)
            ),
        ]
    )
    return (metric_pairings,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 10. Are we at FORM `trace4` parity?

    **The short-trace kernels reach the same output sizes on these inputs.
    FORM remains faster in the amortized native benchmark below.** Term counts
    describe expression size, not correctness: 4D identities allow different
    expressions for the same tensor. Section 11 compares HEP component contractions.

    Closed traces reuse the open-chain adjacent contractions and 4D Chisholm
    identities, including cyclic repeated pairs. A macro dispatches lengths
    1–14 to short kernels. The even kernels cache integer recipes generated
    from the same three-gamma identity as the open-chain epsilon expansion;
    separate tables handle a single gamma5. Unsupported dimensions keep their
    dimension-aware path, and gamma5 remains strictly four-dimensional.

    For ten, twelve and fourteen distinct Lorentz arguments the ordinary
    kernels produce **693, 4,383 and 26,931** terms, respectively. Generic
    pairing gives 945, 10,395 and 135,135. FORM 5.0.0 `trace4` has the same
    smaller counts. Its source uses a general 4D reduction; the generated
    arity-specific cache here is our implementation choice. FORM also has
    cross-trace Chisholm and symmetrization options beyond these kernels.
    See the [FORM Dirac-algebra manual](https://github.com/form-dev/form/blob/master/doc/manual/gamma.tex)
    and [FORM 5.0.0 trace implementation](https://github.com/form-dev/form/blob/v5.0.0/sources/opera.c).

    The optional native comparison checks the selected compact scalar and
    fourteen free indices with FORM `trace4` and `tracen`. It also times the
    installed Idenso fourteen-index trace when explicitly enabled, with its
    first and warm calls reported separately. Rebuilding Python tensor metadata
    for 26,931 terms can be slow, even when the Rust trace kernel is fast. FORM process timings include startup, parsing,
    tracing, sorting and scalar verification; Idenso runs in process. These
    boundaries and potentially different build profiles preclude claiming
    engine parity from one timing ratio.

    For an optimized Rust comparison across all even lengths 2–14, with
    exact scalar checks in both engines and saved FORM programs/logs, run:

    ```sh
    cargo run -p idenso --profile dev-optim --example trace_form -- /path/to/form
    ```
    A recorded native run (2026-09-23, AMD EPYC 9754, Rust 1.98.1,
    Symbolica 3.0.0, `dev-optim`, FORM 5.0.0) measured:

    | Input | Previous Idenso | Current Idenso, warm | FORM `trace4` process | FORM `tracen` process |
    |:--|--:|--:|--:|--:|
    | 10 free indices | 168.05 ms | 1.12 ms | 10.81 ms | 10.63 ms |
    | 12 free indices | Not measured | 8.00 ms | 16.54 ms | 22.57 ms |
    | 14 free indices | Not measured | 63.61 ms | 54.42 ms | 190.45 ms |
    | 14 paired slashes | Not measured | 0.75 ms | 9.12 ms | 8.74 ms |
    | 14 alternating slashes | Not measured | 9.52 ms | 9.58 ms | 9.50 ms |

    **About 150× faster at ten free indices** than the previous Idenso
    implementation in the matching optimized build. The first fourteen-index
    call took 73.12 ms. All native scalar polynomial checks passed. Raw records
    and provenance are in `examples/notebooks/gamma_trace_measurements.json`.
    These are recorded observations on a shared host; the tables below rerun
    the installed software. FORM's roughly 9 ms process overhead hides its
    advantage on small traces, creating an apparent crossover near length 14.

    **Amortizing that overhead changes the comparison.** A follow-up run on
    the same host measured complete FORM processes containing six batches of
    independent free-index traces, including parsing, sorting and disposal
    but excluding scalar verification. Medians of three processes, divided
    by their trace counts, compare with five warm Idenso wall-time samples:

    | Length | Idenso warm | FORM amortized wall | Idenso / FORM |
    |--:|--:|--:|--:|
    | 8 | 0.150 ms | 0.055 ms | 2.75× |
    | 10 | 1.100 ms | 0.386 ms | 2.85× |
    | 12 | 8.094 ms | 2.887 ms | 2.80× |
    | 14 | 61.770 ms | 22.198 ms | 2.78× |

    At length 14, an instrumented copy of the kernel spends about **36.0 ms
    constructing and normalizing products**, **23.5 ms canonicalizing their
    sum**, and **0.15 ms constructing metric atoms**. These are diagnostic
    phase measurements; they do not sum exactly to the uninstrumented time.
    Standalone free traces already bypass the later fixed-point passes, and
    warm timings exclude recipe generation. FORM emits packed metric-index
    records into reusable scratch storage; our kernel constructs general
    Symbolica products and sums. This points to output construction as the
    next optimization target, rather than more cached identities.

    The batch comparison holds more simultaneous expressions in FORM and is
    still a selected workload, not general engine parity. Its raw samples,
    timing boundaries and source hashes are in
    `examples/notebooks/gamma_trace_scaling.json`. Reproduce it with:

    ```sh
    cargo run -p idenso --profile dev-optim --example trace_scaling -- /path/to/form
    ```

    **Where do replacement passes cost us?** The fast ordinary case above is
    only one route. Before the metric-first contractor below, instrumenting
    the complete pipeline found these additional costs in the same optimized
    build (three warm samples unless noted):

    | Input | Full simplification | Main measured cost |
    |:--|--:|:--|
    | 10 ordinary gammas, distinct indices | 1.10 ms | Short kernel; no fixed-point passes |
    | Scalar × that same trace | 150 ms | Schoonschip and normalization traverse the terminal output |
    | 12 gammas with gamma5, distinct indices | 2.56 s | Epsilon stage: 2.36 s (92%) |
    | 14 alternating slashes | 9.66 ms | 11 passes; Schoonschip 3.76 ms, collection 2.02 ms, gamma rewrite 2.50 ms |

    For the axial twelve-factor trace, the epsilon stage visits **44,248 nodes
    per outer pass, twice**. Every nested `term.schoonschip()` call leaves its
    input unchanged, and **zero epsilon identities fire**. About 1.96 s goes
    into those nested calls; actual epsilon-rule checks take about 7.4 ms.
    The first outer pass replaces the trace with a polynomial; the second
    reruns all stages before detecting that the result is unchanged.
    An isolated epsilon-bypass experiment takes about **0.21 s** and produces
    exactly the production result on this input. This is a diagnostic, not a
    valid global instruction to disable epsilon contractions. The first gamma
    rewrite alone takes **2.33 ms** and already equals the final production
    polynomial in this case. A length-14 scalar-prefactor probe takes **8.92 s**,
    compared with **59 ms** for the bare trace (one prefactor sample versus
    three bare-trace samples). Its overhead is again broad traversal.

    This made repeated traversal of terminal outputs and nested Schoonschip
    calls the priority addressed by the cleanup below. Pattern matching itself is a small part of this
    axial cost. For ordinary free traces, output construction remains the
    bottleneck: borrowed factor views improve an isolated length-14 kernel
    from **54.3 to 48.3 ms** (about 11%), while generic product normalization
    and sum merging remain. These separate experiments and their limitations
    are recorded in `examples/notebooks/gamma_trace_profile.json`.

    The trace benchmark CSV distinguishes first-call Idenso time, warm Idenso
    time, FORM `trace4` process time and FORM `tracen` process time. First calls
    in a fresh process expose the lazy table cost; repeated notebook runs may
    already have populated those tables.

    **Earlier metric-first contraction.** These measurements precede the
    shared slot substitution and epsilon cleanup below. At this stage,
    Schoonschip matched
    `g(a_,b_)*remainder___`, then walked the remainder for one compatible
    explicit tensor slot. The replacement joined that slot to the other metric
    endpoint and removed the metric. Compatibility checks representation,
    dimension and duality; scalar parameters and index payloads stay opaque.
    Chain and trace bodies follow the chain-like setting. Substitutions do not
    cross sums or powers, while independent contractions inside them can run.

    The repeated-explicit-index check guards this metric stage. A negative
    result skips metric substitution and chain/vector contraction rules;
    dot normalization and ordinary vector compaction still run. The benchmark reports the raw contractor and
    guarded contractor separately, alongside complete Schoonschip and gamma
    timings. Exact HEP component tests compare evaluations before and after
    rewriting, including antisymmetry and dual representations.

    Five-sample medians on a pinned CPU of the shared host (2026-09-24),
    **isolated contractor only, in microseconds**:

    | Input | Previous tuple rules | Metric first | Speedup |
    |:--|--:|--:|--:|
    | One metric and tensor | 26.75 | 2.31 | 11.6× |
    | Metric and chain body | 37.08 | 4.15 | 8.9× |
    | Metric and cyclic trace | 39.73 | 4.69 | 8.5× |
    | Eight connected metrics | 110.79 | 38.24 | 2.9× |
    | Miss with 48 spectator factors | 58.08 | 43.08 | 1.35× |

    On terminal axial-12 output, the raw failed pass improves from **58.11 to
    23.28 ms**. The repeated-index guard brings the new metric stage to
    **0.382 ms**; the same guard before the old rules takes **0.413 ms**.
    This separates faster contraction from the benefit of avoiding unnecessary
    work. Complete Schoonschip still takes **67.51 ms**, down from **125.44 ms**.

    **Full gamma simplification at that stage versus FORM, in milliseconds:**

    | Axial length | Before metric first | After metric first | FORM `trace4` wall |
    |--:|--:|--:|--:|
    | 6 | 5.79 | 4.78 | 0.00634 |
    | 8 | 45.42 | 36.15 | 0.02023 |
    | 10 | 332.40 | 254.64 | 0.09993 |
    | 12 | 2444.34 | 1791.53 | 0.61944 |

    FORM timings use three measured processes after warmup, six batches per
    process, and include amortized startup, parsing, tracing, sorting and
    cleanup. Rust input construction is outside timing. For axial 12, FORM's
    internal trace-and-sort CPU estimate is **0.593 ms**. At that stage the
    full Idenso pipeline improves **1.36×**; surrounding normalization and
    epsilon processing remain expensive, far from FORM parity.

    `examples/notebooks/metric_contraction.json` records all 25 cases, raw
    samples, build/source identities and FORM programs. At that stage all
    **24 exact HEP tests passed**. The broader suite had 589 passing tests
    and two preexisting failures; complete validation commands and limitations
    are in that record.

    ```sh
    cargo run -p idenso --profile dev-optim --example metric_contraction_benchmark -- /tmp/metric-benchmark
    cargo nextest run -p spenso-hep-lib --test metric_contraction_validation --cargo-profile dev-optim
    ```

    """)
    return


@app.cell
def _(mo):
    time_large_python_trace = mo.ui.checkbox(
        value=False,
        label="Also time fourteen free indices in Python (metadata inference can take minutes)",
    )
    time_large_python_trace
    return (time_large_python_trace,)


@app.cell
def _(momentum_names, scalar_oracle):
    _word = ",".join(momentum_names)
    # The oracle is a short rational polynomial in the declared pp, qq, pq.
    _oracle = str(scalar_oracle)
    compact_form_source = f"""Off Statistics;
Vectors p,q;
Symbols pp,qq,pq;
Local F = g_(1,{_word});
trace4,1;
.sort
id p.p=pp;
id q.q=qq;
id p.q=pq;
.sort
Local Check = F - ({_oracle});
.sort
#$nterms=termsin_(F);
#write "NTERMS=%$",$nterms
#write "CHECK=%E",Check
.end
"""
    _indices = ",".join(f"mu{i}" for i in range(1, 15))
    _projection = "*".join(f"{'p' if i <= 2 else 'q'}(mu{i})" for i in range(1, 15))
    generic_form_source = f"""Off Statistics;
Vectors p,q;
Indices {_indices};
Local F = g_(1,{_indices});
trace4,1;
.sort
#$nterms=termsin_(F);
#write "NTERMS=%$",$nterms
Multiply {_projection};
.sort
Local Check = F - 4*p.p*(q.q)^6;
.sort
#write "CHECK=%E",Check
.end
"""
    return compact_form_source, generic_form_source


@app.cell
def _(
    AUTO,
    S,
    compact_form_source,
    gamma,
    generic_form_source,
    lorentz,
    median,
    perf_counter,
    tr,
    time_large_python_trace,
):
    import os
    import shutil
    import sys

    form_records = []
    form_version = None
    idenso_native_record = None
    # Browser notebooks cannot launch native subprocesses. Native users can
    # select an installed binary explicitly without changing the notebook.
    form_executable = (
        None
        if sys.platform == "emscripten"
        else os.environ.get("FORM_EXECUTABLE") or shutil.which("form")
    )
    if form_executable:
        import re
        import subprocess
        import tempfile
        from pathlib import Path

        form_version = subprocess.run(
            [form_executable, "-v"],
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        ).stdout.strip()
        with tempfile.TemporaryDirectory(prefix="gamma-trace-form-") as _directory:
            _path = Path(_directory) / "trace.frm"
            for _label, _source in (
                ("Compact scalar / trace4", compact_form_source),
                (
                    "Compact scalar / tracen",
                    compact_form_source.replace("trace4,1;", "tracen,1;"),
                ),
                ("Fourteen free indices / trace4", generic_form_source),
                (
                    "Fourteen free indices / tracen",
                    generic_form_source.replace("trace4,1;", "tracen,1;"),
                ),
            ):
                _path.write_text(_source)
                _samples = []
                for _round in range(4):
                    _start = perf_counter()
                    _run = subprocess.run(
                        [form_executable, "-q", str(_path)],
                        cwd=_directory,
                        check=False,
                        capture_output=True,
                        text=True,
                        timeout=30,
                    )
                    _seconds = perf_counter() - _start
                    if _run.returncode:
                        raise RuntimeError(
                            f"FORM failed for {_label}:\n{_run.stdout}\n{_run.stderr}"
                        )
                    assert re.search(
                        r"CHECK=\s*0\s*(?:;|$)", _run.stdout, re.MULTILINE
                    ), _run.stdout
                    _terms = re.search(r"NTERMS=(\d+)", _run.stdout)
                    assert _terms, _run.stdout
                    if _round:
                        _samples.append(_seconds)
                form_records.append(
                    {
                        "case": _label,
                        "terms": int(_terms.group(1)),
                        "process median (ms)": round(1000 * median(_samples), 3),
                        "source": _source,
                        "stdout": _run.stdout,
                    }
                )
        if time_large_python_trace.value:
            _idenso_input = tr(
                *(
                    gamma(AUTO, AUTO, lorentz(S(f"gamma_form::mu{i}")))
                    for i in range(14)
                )
            )
            _start = perf_counter()
            _idenso_result = _idenso_input.simplify_gamma()
            _first = perf_counter() - _start
            _samples = []
            for _ in range(3):
                _start = perf_counter()
                _idenso_result = _idenso_input.simplify_gamma()
                _samples.append(perf_counter() - _start)
            _count = len(list(_idenso_result.to_expression().expand().terms()))
            assert _count == 26931
            idenso_native_record = {
                "case": "Idenso / fourteen free indices",
                "terms": _count,
                "first call (ms)": round(1000 * _first, 3),
                "warm in-process median (ms)": round(1000 * median(_samples), 3),
            }
    return form_executable, form_records, form_version, idenso_native_record


@app.cell(hide_code=True)
def _(
    compact_form_source,
    form_executable,
    form_records,
    form_version,
    idenso_native_record,
    mo,
):
    if form_executable:
        mo.output.append(
            mo.md(
                f"**Native comparison:** `{form_version}`. Every FORM check reduced the scalar difference to zero. Process timings include startup, parsing, tracing, sorting and verification; they are not isolated trace-kernel times."
            )
        )
        mo.output.append(
            mo.ui.table(
                [
                    {
                        key: value
                        for key, value in record.items()
                        if key not in {"source", "stdout"}
                    }
                    for record in form_records
                ],
                selection=None,
                pagination=False,
                show_download=False,
            )
        )
        if idenso_native_record is not None:
            mo.output.append(
                mo.ui.table(
                    [idenso_native_record],
                    selection=None,
                    pagination=False,
                    show_download=False,
                )
            )
        _native_times = {
            record["case"]: record["process median (ms)"] for record in form_records
        }
        _trace4_speedup = (
            _native_times["Fourteen free indices / tracen"]
            / _native_times["Fourteen free indices / trace4"]
        )
        mo.output.append(
            mo.md(
                f"**Within FORM, the fourteen-index pipeline is {_trace4_speedup:.2f}× faster with `trace4` than `tracen`.** This compares the same native executable and projection check. Tiny compact traces may instead be dominated by startup time."
            )
        )
        mo.output.append(
            mo.accordion(
                {
                    record["case"]: mo.md(
                        f"```form\n{record['source']}\n```\n\n```text\n{record['stdout']}\n```"
                    )
                    for record in form_records
                }
            )
        )
    else:
        mo.output.append(
            mo.md(
                "**Native FORM comparison not run.** Browser notebooks cannot launch FORM. Locally, put `form` on `PATH` or set `FORM_EXECUTABLE` to its absolute path before starting Marimo. The program below reproduces the selected scalar check; no FORM timing is invented when it is unavailable."
            )
        )
        mo.output.append(mo.md(f"```form\n{compact_form_source}\n```"))
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 11. Compare tensor values, allowing different simplifications

    The original trace contracts the HEP library's explicit 4×4 gamma
    matrices. The simplified expression contracts its metric and epsilon
    tensors against the **same exact integer four-vectors**. Neither route
    uses the other's symbolic answer as its oracle. Both samples span four
    dimensions, so gamma5 checks cannot pass merely because epsilon vanishes.
    We use the repository's (+−−−) metric and epsilon convention
    $\epsilon_{0123}=-i$, with exact Symbolica components throughout.

    When FORM is available, its ordinary `trace4` and `tracen` outputs are
    imported as metric networks and evaluated by that same HEP library.
    Length 10 is particularly useful: `trace4` has 693 terms and `tracen`
    has 945, yet their four-dimensional component values must agree.
    No comparison requires equal term counts or identical symbolic forms.

    These finite assignments are regression checks, not a proof for all
    tensors. The Rust HEP tests additionally cover every length 1–14,
    gamma5 positions, repeated momenta and cyclic contracted indices.

    ### Recorded comparison: gamma5 and twelve ordinary gammas

    For $\operatorname{Tr}(\gamma_5\gamma^{\mu_0}\cdots\gamma^{\mu_{11}})$,
    both Idenso's **first gamma rewrite** and FORM 5.0.0 `trace4` produce
    **1,029 terms**. The expressions are different: **741 terms coincide**,
    with **288 unique to each**. Their raw symbolic difference has 576 terms.
    The first rewrite equals Idenso's full pipeline output exactly.

    Three separate full-rank integer momentum assignments give these exact
    HEP component contractions of all twelve Lorentz indices:

    | Assignment | Original gamma network | Idenso first rewrite | FORM `trace4` |
    |--:|--:|--:|--:|
    | 1 | −96,536 i | −96,536 i | −96,536 i |
    | 2 | 468,485,024 i | 468,485,024 i | 468,485,024 i |
    | 3 | 114,798,012 i | 114,798,012 i | 114,798,012 i |

    FORM's `d_` and `e_` map to the same HEP metric and epsilon convention.
    The four-gamma axial trace is `4*e_(mu0,mu1,mu2,mu3)` in FORM and
    `4*epsilon(mu0,mu1,mu2,mu3)` in Idenso; no fitted sign is introduced.
    The [FORM metric conventions](https://form-dev.github.io/form-docs/master/manual/#a-few-notes-on-the-use-of-a-metric)
    specify this trace normalization.

    A fresh optimized run measures **2.42 ms** for the first rewrite versus
    **0.654 ms** amortized FORM wall time, about **3.7×** in FORM's favor.
    Full FORM process latency is 11.97 ms; its internal trace-and-sort CPU time
    is 0.624 ms. These are different timing boundaries. The first rewrite is
    already complete on this input, while Idenso's public pipeline still pays
    the previously measured 2.56 s for its extra passes.

    The component assignments, outputs, timings and executable FORM program
    are recorded in `examples/notebooks/gamma_trace_axial_form.json`.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Contraction order and intermediate terms

    This investigation records the behavior before the contraction-order
    update described below.

    Cold Idenso recipe generation and FORM emit the **same raw monomial
    counts** for distinct indices. Warm Idenso uses the collected recipes:

    | Trace | Generated in either engine | Collected recipes / FORM output |
    |:--|:--|:--|
    | Axial 6 / 8 / 10 / 12 | 6 / 33 / 180 / 1053 | 6 / 33 / 180 / 1029 |
    | Ordinary 8 / 10 / 12 | 117 / 801 / 5139 | 105 / 693 / 4383 |

    Axial twelve grows **1 → 1,029** at the first gamma rewrite. Every later
    stage, including the entire second pass, leaves it unchanged. Its slow
    cleanup is not extra trace branching. These are cumulative emissions,
    not peak resident terms or memory.

    **Repeated indices expose a different path.** Let `Tr` denote a 4D trace,
    `a,b` summed Lorentz indices and `p0..p7` distinct slashed vectors:

    ```text
    Tr[a,p0,b,p1,p2,a,p3,p4,b,p5,p6,p7]
      = -2 Tr[a,p0,p4,p3,a,p2,p1,p5,p6,p7]
      =  4 Tr[p3,p4,p0,p2,p1,p5,p6,p7]
    ```

    Idenso selects the shortest cyclic pair, `a`, whose even interior
    branches: ultimately **three length-eight kernels, 315 cached monomials**.
    [FORM prioritizes odd-interior Chisholm contractions](https://github.com/form-dev/form/blob/v5.0.0/sources/opera.c#L828),
    choosing `b` then `a`: **one kernel, 117 raw monomials → 105 collected**.
    Manually taking the first odd reduction changes warm Idenso time from
    **24.73 to 8.33 ms**. Preparing the shorter input is outside timing;
    this preceded the production update below.

    External metrics expose an earlier difference. For
    `g(a,c)*g(b,d)*Tr[a,p0,b,p1,p2,c,p3,p4,d,p5,p6,p7]`, the default
    prepass leaves trace bodies opaque and emits **4,383 cached monomials**.
    Enabling the existing chain-like Schoonschip setting before tracing
    reduces that to **315**, taking **26.64 ms including the prepass** versus
    **656.42 ms** by default. FORM contracts these metrics before tracing,
    then takes the odd-interior path. Timings are five-sample warm medians
    on one pinned CPU of the shared host (2026-09-24).

    Exact HEP network evaluations agree on three nonzero integer assignments
    for the original trace, both odd reductions, the external-metric and
    precontracted inputs, and the repeated-pair and odd-first simplified
    outputs. The large default external-metric output was not HEP-evaluated.
    All 24 instrumented outputs match linked production outputs exactly.
    `examples/notebooks/gamma_trace_intermediates.json` records the counts,
    exact assignments, raw timings, source identities and diagnostic drivers.
    Counts preserve factorization; outer sum arity alone hides nested sums.

    ### Earlier contraction-order update

    The next update made gamma simplification contract metrics into chains and
    traces before **every rewrite pass**, including metrics generated by a neighboring
    trace. The repeated-index scan guards metric substitution and chain/vector
    rules, so terminal expressions without repeated explicit indices skip
    those operations without tensor parsing.

    Cyclic rotation and Chisholm selection use the same priority: adjacent
    pairs, applicable odd interiors, applicable even interiors, then other
    repeated pairs. Chisholm priority requires explicit 4D indices and a
    gamma-only interior. An interval crossing gamma5 cannot hide an applicable
    interval on the other side. Open-chain and generic-dimensional pair
    ordering retain their existing behavior.

    After this update, the ordinary crossing example has **315 → 105 arithmetic
    leaves**, its external-metric version **4,383 → 105**, and the axial crossing
    example **99 → 33**. Counts preserve factorization and do not collect
    terms; they are neither peak-memory measurements nor equality checks.

    These historical measurements precede the cleanup below. The expanded
    native benchmark covers 35 inputs, including 16 complete
    gamma pipelines and 12 matching FORM cases. Five warm Idenso samples
    use the same optimized driver and one pinned CPU before/after; input
    creation is excluded. FORM wall time includes amortized startup, parsing,
    tracing, sorting and cleanup, with six batches in each of three measured
    processes. The host is shared.

    | Input | Before (ms) | After (ms) | FORM wall (ms) |
    |:--|--:|--:|--:|
    | Repeated-index trace, 12 | 19.52 | 6.62 | 0.06043 |
    | External metrics, 12 | 508.55 | 6.74 | 0.06144 |
    | Axial crossing, 12 | 105.96 | 32.22 | 0.02116 |
    | Two-gamma interior, 8 | 0.269 | 0.286 | 0.00491 |
    | Cyclic adjacent pair, 12 | 50.98 | 51.16 | 0.40266 |
    | Free indices, 8 | 0.148 | 0.148 | 0.05719 |
    | Free indices, 12 | 6.24 | 6.32 | 2.91050 |
    | Axial distinct indices, 12 | 1769.61 | 1785.88 | 0.62148 |

    The main ordinary cases improve by **2.95× and 75.4×**, and the axial
    crossing case by **3.29×**. The small two-gamma-interior case costs about
    **17 µs more**. At this stage, axial twelve with distinct indices remains
    dominated by cleanup, and FORM remains substantially faster on the
    contracted cases.

    Ten new unit regressions cover cyclic priority, gamma5 boundaries,
    symbolic dimensions, inert traces and generated-metric idempotence.
    Three new HEP regressions compare exact nonzero component evaluations,
    including the production external-metric output. Raw timings, factorized
    output counts, matching FORM programs, source identities and validation
    results are recorded in `examples/notebooks/gamma_trace_ordering.json`.

    ### Earlier shared slot substitution and epsilon cleanup

    These measurements record the shared-slot update before the local metric
    normalization below.

    Metrics and tagged vectors now use **the same slot substitution**. The
    contractor scans product factors directly, trying metric sources before
    tagged rank-one vectors; finding a source no longer uses the
    metric/remainder pattern. It joins one compatible explicit slot to the
    opposite metric endpoint or compact vector, then removes the source
    factor. Representation, dimension and duality checks, opaque scalar
    parameters, and chain-like and rank-one settings are retained. The
    repeated-index scan still guards this work, preserving sums and powers
    without distributing them.

    Epsilon cleanup stays a separate `simplify_epsilon()` operation. When
    epsilon is present, it retains **one initial whole-expression Schoonschip
    pass**. Candidate selection directly distinguishes powers, products and
    epsilon function heads. It no longer runs Schoonschip at every visited
    node: another whole-expression pass runs only after an epsilon identity
    changed the expression, to consume the generated metrics. Powers and
    epsilon pairs with disjoint explicit indices still receive their identities
    even when the repeated-index predicate is false.

    The fresh comparison runs the same benchmark driver against the saved
    baseline and new libraries: **50 inputs, 316 case/method pairs, all
    before/after output snapshots byte-identical**. Rust timings are medians
    of five warm samples with adaptive repetitions targeting 8 ms per sample,
    pinned to CPU 6 of the shared host; input creation is excluded. Complete
    Schoonschip timings, including normalization and the guard, in microseconds:

    | Input | Before | After | Speedup |
    |:--|--:|--:|--:|
    | Tagged vector and tensor | 44.33 | 19.51 | 2.27× |
    | Tagged vector and trace | 62.62 | 29.76 | 2.10× |
    | Tagged compact metric | 49.34 | 24.46 | 2.02× |

    Isolated epsilon cleanup improves **16.401 → 1.418 ms (11.57×)** for a
    pair with eight distinct explicit indices, and **827.076 → 58.705 ms
    (14.09×)** on the already simplified axial-twelve output. Both include
    the retained initial Schoonschip pass; the trace kernel and its final
    **1,029 terms** are unchanged.

    Full gamma timings and a fresh **FORM 5.0.0** run, in milliseconds:

    | Input | Before | After | FORM wall |
    |:--|--:|--:|--:|
    | Axial distinct indices, 6 | 4.854 | 0.637 | — |
    | Axial distinct indices, 8 | 36.759 | 4.895 | — |
    | Axial distinct indices, 10 | 256.810 | 35.781 | — |
    | Axial distinct indices, 12 | 1805.957 | 259.161 | 0.61493 |
    | Axial crossing, 12 | 32.130 | 3.103 | 0.02139 |
    | External metrics, 12 | 6.846 | 6.105 | 0.06142 |
    | Free indices, 8 | 0.153 | 0.148 | 0.05617 |
    | Free indices, 12 | 6.759 | 6.477 | 2.87345 |

    FORM uses one warmup process and three measured processes, each with
    six batches, pinned to CPU 3. Its wall time includes amortized startup,
    parsing, input normalization, `trace4`, sorting and cleanup. A dash marks
    a case without a fresh matching FORM run.

    Axial twelve improves **6.97×**, and the axial crossing case **10.36×**.
    Ordinary standalone traces change little because they already bypass
    most of this cleanup. FORM remains about **421× faster** for axial twelve
    under these timing boundaries; this does not establish parity.
    `examples/notebooks/gamma_cleanup.json` records all cases, raw samples,
    output comparisons, build/source identities and the fresh FORM programs.

    **Remaining cost on unchanged axial-twelve output:** an instrumented copy
    measures **28.45 ms** for one `normalize_dots()` call, including **17.13 ms**
    in redundant-metric rules and **6.36 ms** in metric-trace rules. The epsilon
    candidate pass takes **1.09 ms**, the repeated-index scan **0.398 ms**, and
    complete Schoonschip **59.80 ms**. These are isolated five-round warm
    medians on CPU 6, with outputs checked against production; they do not add
    up to a decomposition of the full gamma pipeline.

    **HEP validation limitation:** an explicit vector product with dual slots
    evaluates, but both the saved baseline and the new compact-dot output
    discard the dual orientation, so the tensor-network parser rejects that
    output. The diagnostic `tagged_vector_dots_preserve_dual_slot_orientation`
    remains explicitly ignored with this reason; passing component checks do
    not certify this case. The measurement record separates this limitation
    from the unchanged benchmark snapshots.

    ### Applied local metric normalization

    These measurements use **crates.io Symbolica 3.0.0**. The main-branch
    comparison below retains this normalizer and changes the pinned dependency.

    The previous axial-twelve profile spent **17.128 ms** checking redundant
    metrics and **6.361 ms** checking metric traces inside one **28.448 ms**
    `normalize_dots()` call. Every result was unchanged. These stages repeatedly
    walked an already normalized expression and attempted wildcard matches
    at eligible metric heads.

    Intrinsic metric identities now run **locally when Symbolica constructs
    `g`**: compatible repeated abstract slots yield their dimension, and a
    compatible explicit slot with a tagged compact vector yields that vector
    component. Scalar linearity still applies. `normalize_dots()` handles a
    root-level identity directly, then uses a **read-only visitor** to check for
    remaining nested-vector or integral-power identities. If none applies,
    it returns the existing expression without rebuilding it. Otherwise, one
    rewriting walk uses borrowed slot recognition instead of whole-expression
    replacement rules. Scalar parameters and slot payloads stay opaque, and
    numerator sums stay factored.

    This check matters even after removing wildcard rules: assigning an
    unchanged node in Symbolica's rewriting walker to prune its children marks
    that path as changed and rebuilds its ancestors. The visitor can prune the
    same opaque payloads without triggering reconstruction or normalization.

    Recognition requires canonical slot shapes, exact dimensions, compatible
    representations and variance, and a final structural argument for vectors.
    Malformed arguments and dimension mismatches remain explicit. Fixed
    components marked by `cind` or `find` are excluded from abstract-index
    traces and power contractions. For abstract powers, the existing convention
    remains: negative even powers normalize; negative odd and nonintegral
    powers remain explicit.

    Construction also makes scalar composition consistent: **`g(i,i) = D`
    implies `g(i,i)^2 = D^2` and `g(i,i)^3 = D^3`**. The former ordering applied
    the metric-power rule first and incorrectly gave `D` for the square.
    The distinct-index identity `g(i,j)^2 = D` remains a contraction. Exact
    HEP regressions use a separate tensor with metric components and no
    normalization callback, so both sides do not simply collapse during
    construction before the trace-power check.

    The expanded benchmark retains identical source strings and separately
    times parsing, `normalize_dots()` on constructed atoms, and parsing followed
    by normalization. **Parse-only includes constructor callbacks;
    parse-plus-normalize exposes their total cost**, including work moved out
    of the later pass. The run covers **73 inputs and 641 case/method pairs**.
    Each Rust median uses five warm samples, rotating method order and adapting
    repetitions toward 8 ms per sample, pinned to CPU 6 of the shared host.
    Construction-inclusive timings, in microseconds:

    | Input | Parse before | Parse after | Parse + normalize before | Parse + normalize after |
    |:--|--:|--:|--:|--:|
    | Metric self-trace | 3.105 | 3.822 | 12.391 | 3.860 |
    | Metric and tensor | 4.880 | 5.471 | 16.160 | 6.051 |
    | Metric and compact vector | 3.713 | 4.737 | 12.107 | 4.971 |
    | Mismatched dimensions | 3.124 | 3.633 | 13.148 | 4.057 |
    | Compact dot with scalar factors | 5.751 | 6.365 | 11.037 | 6.880 |

    Parsing costs more in these examples because local metric checks now
    happen there. **The combined operation improves with that cost included.**
    Constructed input snapshots can change eagerly; recorded source strings
    remain identical. The trace-power correction and stricter canonical
    recognition also intentionally change some outputs, which the measurement
    record distinguishes from unchanged output comparisons.

    On unchanged axial-twelve output, isolated `normalize_dots()` improves
    **27.794 → 0.858 ms (32.39×)**, and complete Schoonschip improves
    **59.577 → 5.265 ms (11.32×)**. Full gamma timings exclude initial input
    creation but include metrics constructed while evaluating the trace.
    Matching **FORM 5.0.0** wall times include amortized startup, parsing,
    tracing, sorting and cleanup: one warmup process, then three measured
    processes with six batches each on CPU 3. All entries are milliseconds:

    | Input | Before | After | FORM wall |
    |:--|--:|--:|--:|
    | Axial distinct indices, 6 | 0.647 | 0.204 | — |
    | Axial distinct indices, 8 | 4.952 | 1.141 | — |
    | Axial distinct indices, 10 | 36.379 | 7.161 | — |
    | Axial distinct indices, 12 | 264.302 | 49.084 | 0.62948 |
    | Axial crossing, 12 | 3.268 | 1.643 | 0.03008 |
    | External metrics, 12 | 6.291 | 3.524 | 0.06158 |
    | Free indices, 8 | 0.150 | 0.162 | 0.05679 |
    | Free indices, 12 | 6.622 | 6.496 | 2.91656 |

    Axial twelve improves **5.38×**. Free eight becomes about **7.9% slower**;
    its shortcut already avoids cleanup, while constructing already-normal
    metrics can incur the new checks. Free twelve changes little. FORM remains
    about **78× faster** for axial twelve under these timing boundaries.
    `examples/notebooks/dot_normalization.json` preserves all cases, raw
    samples, input/output snapshot hashes, classified changed outputs,
    source identities and matching FORM programs. The benchmark command
    recreates the complete snapshots. A dash marks a case without a fresh
    matching FORM run.

    Eight dot-normalization calls process the large axial-twelve output.
    Multiplying the isolated per-call reduction by eight estimates **215.5 ms
    saved**, close to the **215.2 ms** reduction in the complete pipeline.
    This is a timing-based attribution estimate, not an exact measurement of
    those calls inside the pipeline.

    A separate final profile uses five warm samples on CPU 7 and checks that
    every stage preserves the terminal expression. It measures **0.867 ms**
    for dot normalization, **5.330 ms** for Schoonschip and **6.722 ms** for
    epsilon cleanup. In that measurement, chain collection was the largest stage:
    `chainify` takes **20.957 ms** and `collect_gamma_chains` **24.475 ms**;
    simplifying the already-terminal expression takes **39.275 ms**. These
    stage measurements overlap and cannot be added to reconstruct full gamma
    runtime.

    The selected suites have **632 passing tests**: 336 in Idenso, 261 in
    Spenso with shadowing enabled, and 35 HEP component tests. Two previously
    recorded tests still fail, and 23 are skipped, including the known compact
    dual-dot variance diagnostic above. The record names these limitations
    and retains commands and logs. All three scoped Clippy checks pass with
    warnings denied.

    ### Pinned Symbolica main comparison

    Active Rust workspaces and the Python lock now pin Symbolica main at
    **`06906976bca24fefc5203aee699d90d62ebe08cd`**. Its package version still
    reads 3.0.0; the source revision distinguishes it from the registry release
    above. The unchanged benchmark driver covers **73 inputs and 641
    case/method pairs**. All **787 plain-text source, input and output
    snapshots match** the optimized registry baseline exactly.

    Full gamma medians, in milliseconds. Both executables were measured
    afresh, registry first and pinned main second, with five warm samples on
    CPU 6. Our validation compiler group was paused for both runs and resumed
    afterward. Other shared-host load remains uncontrolled; builds were not
    interleaved.

    | Input | Registry 3.0.0 | Pinned main |
    |:--|--:|--:|
    | Axial distinct indices, 8 | 1.145 | 1.091 |
    | Axial distinct indices, 10 | 7.172 | 6.870 |
    | Axial distinct indices, 12 | 48.962 | 46.879 |
    | Axial crossing, 12 | 1.645 | 1.621 |
    | External metrics, 12 | 3.524 | 3.628 |
    | Free indices, 8 | 0.160 | 0.159 |
    | Free indices, 12 | 6.399 | 6.726 |

    **The differences are small and mixed.** Axial-twelve dot normalization
    takes **0.889 → 0.848 ms**, and its full pipeline is about **4.3% faster**.
    Free twelve is about **5.1% slower**, and the external-metric case about
    **3% slower**. This sequential shared-host comparison does not establish a
    broad upstream speedup.
    `examples/notebooks/symbolica_main_update.json` retains both raw timing
    sets, source and lock identities, build metadata, snapshot hashes and
    validation results. FORM numbers in that record are reused from the
    preceding comparison; FORM was not rerun for this dependency update.

    The binary atom export format changes **0 → 1**, and this revision rejects
    legacy format-0 payloads. The surrounding full state-export format also
    moves **5 → 6**, requiring compatible payload consumers. To migrate saved
    atoms, export expressions as strings with the older version and parse them
    using the new version.
    Matching printed outputs does not establish binary-format compatibility.

    The standalone native notebook host and Tydenso **`wasm32` dependency checks
    pass** at this revision. It includes the upstream fix for decoding a stored
    64-bit length on a 32-bit target. Tydenso also passes
    `preserve_indices = false` at its Dirac-adjoint call site, preserving its
    previous behavior. Separate native and Wasm runtime probes check binary
    round trips, malformed lengths and legacy-format rejection; exact
    cross-platform byte equality covers the zero fixture. Fifteen existing
    Idenso snapshots were updated only for imaginary-unit display; their
    inputs and algebra are unchanged. The selected native suites retain
    **632 passing tests**, the same two known failures and 23 skipped tests.
    All three scoped Clippy checks pass with warnings denied. These native
    results do not certify the separate Tydenso/Tymbolica runtime payload
    exchange.

    That runtime exchange remains **blocked**: the pinned atom-payload
    dependency `18e04916` accepts outer export version 5 and rejects version 6
    before import. The record keeps the reproduced failure and an isolated
    candidate patch; the candidate is not applied to the workspace dependency,
    and the tracked plugin assets retain their previous build. In isolation,
    the candidate passes 22 payload tests against each of registry and
    pinned-main Symbolica, 21 Tydenso runtime checks, and two-way exchange with
    the rebuilt Tymbolica engine.

    ### Further contraction and cleanup improvements

    The next optimization keeps Symbolica pinned at `06906976` and changes how
    Idenso schedules contractions and cleanup. A product finishes its local
    contractions before rebuilding enclosing expressions, while preserving the
    metric-first order. Compact vectors are constructed only after a compatible
    partner is found. Borrowed visitors avoid rebuilding unchanged slot payloads;
    chain, bracket and gamma passes skip work when their required symbols are absent.
    Schoonschip repeats post-contraction cleanup only when a contraction changed
    its input; its outer fixed point still handles initial normalization changes.

    The existing terminal-trace shortcut now also accepts a standalone 4D spin
    trace with exactly one gamma5 and at most fourteen ordinary gammas, whose
    Lorentz indices are distinct explicit `mink(4,index)` slots. Each kernel
    monomial uses every index once and contains at most one epsilon, so it needs
    no further contraction cleanup. Repeated indices, compact slashes and
    spectators retain the full pipeline. The ordinary free-trace shortcut is
    unchanged in scope.

    Full gamma medians below are milliseconds. These are five warm samples on
    CPU 6, using the reversed after-before confirmation run. FORM 5.0.0 was
    remeasured on CPU 3. Its amortized wall timing includes input preparation and
    process overhead; its separate internal timer covers `trace4` plus sorting.
    Rust input preparation is outside timing, and output destruction is included.
    Our builds were quiet, but unrelated shared-host load was uncontrolled.

    | Input | Before | After | FORM wall |
    |:--|--:|--:|--:|
    | Axial distinct indices, 8 | 1.112 | 0.119 | 0.021 |
    | Axial distinct indices, 10 | 6.794 | 0.467 | 0.102 |
    | Axial distinct indices, 12 | 46.343 | 2.512 | 0.646 |
    | External metrics, 12 | 3.332 | 0.930 | 0.062 |
    | Free indices, 8 | 0.162 | 0.166 | 0.060 |
    | Free indices, 12 | 6.542 | 6.776 | 2.939 |

    Axial twelve improves by **18.45 times**. FORM remains **3.89 times faster**
    on that case, so this does not establish performance parity. Ordinary free
    eight and twelve are 2.3% and 3.6% slower in this confirmation run; both were
    slightly faster in the first pair. Those small mixed changes are consistent
    with shared-host variation, rather than an established gain.

    Productive public Schoonschip calls improve too; the following medians are
    microseconds from the same confirmation run.

    | Contraction | Before (µs) | After (µs) |
    |:--|--:|--:|
    | Vector into tensor | 5.115 | 3.356 |
    | Metric into epsilon | 7.758 | 4.935 |
    | Metric into chain | 9.294 | 6.129 |
    | Metric into trace | 9.838 | 6.563 |
    | Eight-metric chain | 26.562 | 23.753 |

    The original before-after run is retained alongside the confirmation. It
    measured axial twelve at **46.658 → 2.383 ms**. Unchanged parsing controls
    exposed a localized timing disturbance in several small cases, which motivated
    the reversed repeat; the confirmation resolves those apparent regressions.
    No confirmation case above one microsecond regresses by more than 10%.

    A separate public-stage profile covers **92 expression shapes and 860
    case/method pairs**, including both original and already-terminal traces.
    On the terminal axial-twelve output, chain collection takes **23.635 → 1.209 ms**,
    Schoonschip **5.280 → 1.567 ms**, epsilon cleanup **6.168 → 2.537 ms**, and full
    gamma simplification **37.017 → 7.038 ms**. These overlapping stage timings
    cannot be added to reconstruct the original trace call. The shortcut explains
    why simplifying the original axial trace is now faster than processing its
    large terminal output through the general pipeline.

    Both full comparisons cover **73 inputs and 641 case/method pairs**. All
    **787 full-driver snapshots and 952 phase snapshots match exactly**. The fresh
    FORM run covers fifteen cases; twelve full polynomials match the historical
    records exactly. The additional axial-six/eight/ten cases have fresh
    counts and timings, without a new FORM component certificate. Term counts
    are never used as an equivalence test.

    The current validation has **383 passing tests**: 344 Idenso tests and 39 HEP
    component tests. The known Idenso parsing snapshot failure remains; 23 tests
    are skipped. Spenso's separate suite was not rerun in this round. Both scoped
    Clippy checks pass with warnings denied. The new HEP checks simplify explicit
    axial traces before attaching spectators, then independently evaluate their
    metric/epsilon tensors against the same vectors as the original gamma-matrix
    network. Lengths four, eight and twelve, three full-rank samples, and several
    gamma5 positions preserve exact Gaussian-integer values.

    The first length-twelve oracle timed out while materializing a tensor with
    twelve open Lorentz slots. Binding each vector inside its raw gamma before
    HEP takes the matrix trace avoids that intermediate without using gamma
    simplification in the oracle. The record retains the timeout and the final
    passing checks.

    The complete record is `examples/notebooks/contraction_performance.json`, including both run orders, raw samples, source identities, FORM programs and validation evidence. Reproduce the current Rust measurements with:

    ```sh
    cargo run --locked -p idenso --profile dev-optim --example metric_contraction_benchmark -- /tmp/contraction-full 5 8
    cargo run --locked -p idenso --profile dev-optim --example contraction_phase_benchmark -- /tmp/contraction-full /tmp/contraction-phases 5 8
    ```
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Generic-dimensional traces and contraction progress

    The 2026-09-24 comparison uses **symbolic Lorentz dimension D** and
    **Tr(1) = 4**. Canonical ordinary traces now reduce summed pairs inside the
    memoized pairing evaluator, including branches with two or more interior
    gammas. Adjacent pairs contribute $D$; one-gamma sandwiches contribute
    $2-D$. Longer sandwiches use the Clifford recurrence, with the existing
    four-dimensional odd-interior Chisholm identity when applicable. The
    branches share subword results and construct metrics only when needed.
    This avoids intermediate trace atoms and repeated general cleanup.

    The direct route requires canonical slots, a uniform dimension, and at most
    two occurrences of each explicit index. Fixed components are not summed.
    Remaining compact slashes reuse adjacent squares and
    $\not p\not q\not p=2(p\cdot q)\not p-p^2\not q$. Scalar spectators retain
    their factorization, and gamma-five remains strictly four-dimensional.

    The latest warm milliseconds use optimized Rust, Symbolica main `06906976`
    and FORM 5.0.0, pinned to CPU 9 on the same shared EPYC host. Idenso uses
    in-process wall time including output destruction; FORM uses its internal
    `tracen`-plus-sort CPU timer, amortized over independent expressions.
    Parsing and startup are excluded. Idenso uses three processes per version
    with five warm batches each, alternating before/after order; FORM uses
    three processes with three retained batches each. Other host activity is
    uncontrolled, so small differences are not established improvements.

    | Free gammas | Idenso, factored | FORM, expanded |
    |--:|--:|--:|
    | 8 | 0.145 | 0.036 |
    | 10 | 0.483 | 0.337 |
    | 12 | 1.997 | 3.900 |
    | 14 | 18.184 | 53.667 |

    **Expanded-output performance is not at parity.** Idenso's factored results
    defer substantial work. A separate five-sample experiment measures trace
    construction, expansion, and destruction of both outputs:

    | Generic-D input | Trace, ms | Expansion, ms | Full lifecycle, ms | FORM, ms |
    |:--|--:|--:|--:|--:|
    | 8 gammas, branching pair | 0.058 | 0.050 | 0.109 | 0.0134 |
    | Pair with four interior gammas | 0.441 | 1.295 | 1.731 | 0.173 |
    | Pair with five interior gammas | 2.230 | 19.728 | 21.894 | 2.550 |
    | 12 gammas, order-sensitive pairs | 1.064 | 4.892 | 5.895 | 0.763 |
    | 14 free gammas | 39.047 | 739.466 | 789.138 | 53.667 |

    The fourteen-gamma polynomial has 135,135 terms and takes about **15×**
    FORM's time including expansion. Each column is a median of its own measured
    quantity; the total includes destruction. This separate experiment has a
    different allocator state from the first-call table. Diagnostic expansion
    applies only to standalone trace polynomials, never to graph numerators or
    spectators. The earlier checkpoint's 776 ms is comparable; the difference
    is not an attributed optimization. FORM's earlier length-ten samples ranged
    from 0.337 to 0.503 ms, another reason to avoid a parity claim there.

    A separate six-case expansion profile isolates Symbolica's expression work.
    At length 14, expansion alone takes **756.9 ms** with tensor slots,
    **635.3 ms** with shallow metrics, and **580.2 ms** with scalar variables.
    Polynomial expansion reduces the variable case to **414.2 ms**, but
    converting it back to metrics takes another **353.6 ms**, before forward
    substitution. All **1,080 inverse-substitution checks** pass. Removing tensor
    syntax therefore does not provide an end-to-end improvement. The length-12
    instruction profile attributes about 19.8% to normalization, 11.7% to
    product/sum buffer extension, and 8.3% each to factor comparison and copying.
    Pattern matching and slot parsing do not appear in this expansion capture.
    These percentages count instructions, not wall time.

    A direct-output prototype emits the final pairings without generic expansion.
    Five samples on CPU 6 measure **266.0 ms** including destruction for length
    14, with exact equality to the production expanded polynomial. Another
    prototype builds an integer polynomial in **19.4 ms** but spends **248.8 ms**
    converting it to an Atom. These are isolated experiments, not production
    improvements, and still leave roughly a **5×** gap to FORM. Caching expanded
    subexpressions regresses most smaller cases and was rejected. The pinned
    standalone program in `examples/reproducers/symbolica-expansion` reproduces
    the expansion bottleneck without Spenso or Idenso.

    Direct branching contraction substantially improves the difficult summed
    words. The baseline already contains the previous adjacent/one-gamma
    reductions; remaining gammas are compact slashes.

    | Generic-D input | Before, ms | Current, ms | FORM, ms |
    |:--|--:|--:|--:|
    | 8 gammas, branching two-gamma interior | 1.142 | 0.053 | 0.0134 |
    | Pair with three interior gammas | 2.085 | 0.096 | 0.018 |
    | Pair with four interior gammas | 18.148 | 0.412 | 0.173 |
    | Pair with five interior gammas | 228.236 | 1.851 | 2.550 |
    | 8 gammas, crossing pairs | 1.876 | 0.042 | 0.0088 |
    | 10 gammas, crossing pairs | 16.863 | 0.192 | 0.070 |
    | 10 gammas, nested branching pairs | 8.794 | 0.112 | 0.0408 |
    | 12 gammas, order-sensitive pairs | 160.308 | 0.964 | 0.763 |

    The branching-eight case improves **21.5×**, with a remaining **4×** gap
    before expansion. Its preceding checkpoint used general rewriting and was
    about 55× slower than that checkpoint's FORM timing. The direct identities
    now cover it and longer interiors. Short compact words still have appreciable
    overhead: paired fourteen slashes take **11.6 µs** versus **0.8 µs**, and
    alternating slashes take **36.0 µs** versus **9.4 µs**. A four-dimensional
    five-interior control improves **0.626 → 0.527 ms**.

    At the shared-scan checkpoint, unchanged generic-D reruns measure:

    | Free gammas | Idenso, ms | FORM, ms |
    |--:|--:|--:|
    | 8 | 0.052 | 0.031 |
    | 10 | 0.483 | 0.277 |
    | 12 | 5.340 | 3.400 |
    | 14 | 72.235 | 46.667 |

    The original fourteen-gamma rerun took about 142 ms. The shared scan first
    reduced it to 81.6 ms; a final paired comparison of that build against the
    leaner observer measures **83.3 → 72.2 ms**. The final check includes cached
    slot metadata, uses three alternating process rounds, and preserves every
    measured output exactly. The earlier 40 ms walker prototype is a separate
    experiment and is not the production result.

    The repeated-explicit-index walk also observes bracket, dot and
    Dirac-operation candidates. It writes head flags only on a hit and avoids
    observing an unwrapped representation head twice. The observer now finishes
    this same syntactic walk after a repeated index, with further index
    bookkeeping disabled. A metric chain can therefore establish that bracket
    and dot normalization have no work. The boolean-only predicate still stops
    at its first hit. Compound opaque payloads leave later passes enabled, and
    any rewrite invalidates the observations. The fourteen-gamma result still
    has 2.68 million tree nodes.

    All **39 case/dimension combinations** pass **480 exact coordinate
    assignments** comparing the previous output, current output and FORM against
    independent Clifford multiplication in four and six dimensions. All
    **21 symbolic-D cases** additionally match FORM's full polynomial exactly.
    Tests cover interiors two through five, nested/crossing pairs, repeated
    compact slashes, free slots and scalar spectators; eliminated summed labels
    are absent from the results. The final scoped optimization preserves every
    output hash from that complete certificate.

    Earlier contraction checkpoints remain useful context: full axial-twelve
    simplification improved **2.356 → 0.621 ms**, ordinary free-twelve
    **6.422 → 2.794 ms**, and public Schoonschip on a 64-metric chain
    **738.4 → 34.7 µs**. FORM's separately measured metric-chain CPU time was
    **3.2 µs**. A single metric hit was slightly slower (**2.90 → 3.11 µs**).
    These are historical measurements, not the branching comparison above.

    A later isolated `SlotMatcher` comparison caches the representation tag and
    variance-wrapper IDs once per matcher. This removes repeated Symbolica
    initialization probes when many distinct tensor heads miss the small cache.
    Public Schoonschip timings, including normalization, are:

    | Input | Before, µs | Cached metadata, µs | FORM, µs |
    |:--|--:|--:|--:|
    | Metric applied to 8 tensor terms | 11.98 | 10.51 | 4.0 |
    | Metric applied to 128 tensor terms | 176.65 | 141.68 | 66.7 |
    | Vector applied to 8 tensor terms | 13.74 | 12.06 | 4.5 |
    | Vector applied to 128 tensor terms | 196.85 | 166.19 | 73.3 |

    These paired Idenso measurements use CPU 7; the separately measured FORM
    reference uses CPU 3 and times contraction plus sorting after input setup.
    Metric chains are essentially unchanged at about 33 µs for 64 metrics at
    that checkpoint. Slot syntax is unchanged.

    The subsequent metric-component change parses endpoints while building the
    incidence map and interns exact representation/dimension pairs. Each metric
    has at most two neighbors, so walking from both sides directly finds the
    path boundaries or closes a loop, without DFS worklists. Private contractor
    times on CPU 7 improve **18.85 → 14.89 µs** for a 64-metric path,
    **16.56 → 12.58 µs** for a loop and **37.53 → 30.62 µs** for disjoint paths.
    All **221 exact-output cases** pass, including shuffled labels, mixed
    dimensions, duality, ambiguous incidences and factored spectators. These
    times exclude public normalization and cannot be compared directly with
    FORM's complete contraction time. Callgrind records **19.3% fewer
    instructions**; endpoint decoding is unchanged and now accounts for 37.3%
    of the private cost. All **36 Schoonschip unit tests** and **46 HEP integration
    tests** pass, with Clippy clean.

    The paired public API comparison, including normalization, improves
    64-metric paths **32.86 → 28.63 µs**, loops **30.13 → 25.34 µs**, and
    disjoint paths **64.17 → 55.64 µs**. All **41 outputs agree exactly** and the
    dependency identities match; sum controls remain effectively unchanged.
    These measurements use 11 alternating process pairs on CPU 7. FORM was not
    rerun for this change: its earlier **3.2 µs** path reference still leaves
    an approximate **9×** gap to the public contractor.

    Completing the shared observer after its first repeated index subsequently
    improves public 64-metric paths **28.09 → 23.44 µs**, loops
    **25.46 → 20.58 µs**, and metrics over 128 tensor terms
    **140.60 → 124.07 µs**. Vectors over sums remain effectively unchanged.
    All **61 public outputs agree exactly** across 11 alternating process pairs
    on CPU 7. Actual fourteen-gamma reruns measure **72.44 → 73.53 ms**, with
    the paired variability interval including no change. This optimization
    avoids unnecessary bracket and dot passes before metric contraction; it
    does not reduce the free-trace scan. Scoped validation passes **12 parser/slot
    tests, 36 Schoonschip tests, 46 HEP integration tests**, and Clippy.

    A later cleanup experiment reuses the contractor's final traversal to detect
    whether its changed output needs bracket or dot normalization. Input
    observations alone are insufficient because tensor callbacks can introduce
    either form. The copied-library trial matches all **249 fixtures under four
    settings**, but is not retained: metrics over 128 terms show a plausible gain
    of about **119 → 107 µs**, while paired runs are noisy and small unchanged
    controls add measurable overhead. Trace controls show no regression. The
    contraction record preserves the phase profile, callback counterexamples and
    full timing ranges.

    Isolated Symbolica patches improve conversion of the length-14 polynomial
    back to an expression **199.0 → 104.1 ms**. Full polynomial expansion
    improves **700.8 → 605.7 ms**, while ordinary expansion remains about
    **670 ms**. A flat-sort addition prototype initially regresses existing-sum
    inputs. Keeping the full merge cursors outside the heap resolves those
    measured regressions: adding the final length-14 monomials takes
    **158.2 → 35.7 ms**, and adding 512 existing sums takes
    **595.3 → 572.5 µs**. All **1,017 raw input/output/rerun cases** and seven
    upstream normalization tests pass. These are isolated primitive measurements;
    the patches are not applied to the production dependency. Their boundaries,
    regressions and paired measurements are recorded in
    `examples/reproducers/symbolica-expansion/performance.typ` and its companion
    `primitive_measurements.json`.

    A separate one-pass dimension/index accessor trial does not survive the
    public comparison. Despite isolated decoding gains, the 64-metric path
    slows **23.24 → 24.28 µs**, with similar **4–5%** regressions in loops and
    gamma reruns. All **61 outputs remain exact**, but the accessor change is
    rejected and the preceding observer implementation is retained.

    The final retained-build checkpoint refreshes three generic-D traces and
    their FORM programs. Times below are **warm milliseconds**; factored tracing
    and the full expanded lifecycle are separate measurements.

    | Input | Idenso, factored | Idenso, fully expanded | FORM, fully expanded |
    |:--|--:|--:|--:|
    | Branching 8 | 0.0506 | 0.1165 | 0.0132 |
    | Order-sensitive 12 | 0.9564 | 5.9803 | 0.7333 |
    | Free 14 | 17.4870 | 798.5513 | 54.0000 |

    The full-output gaps remain approximately **8.8×, 8.2× and 14.8×**. Rust
    uses wall time on CPU 9, including trace construction, expansion and
    destruction; FORM uses its internal tracing-plus-sort CPU timer on CPU 7.
    Parsing and process startup are excluded. For free length fourteen,
    expansion takes **753.20 ms**, about **94%** of the full lifecycle.
    The retained build's
    fourteen-gamma rerun takes **73.58 ms** against FORM's **46.67 ms**.
    Three processes with five samples each pass all **nine comparisons with
    historical factored outputs**, and all reruns are unchanged. Cold
    first calls are retained separately in the record. Different allocator and
    working set states mean the factored-only timing must not be added to an
    independent expansion timing to estimate the full lifecycle.

    A direct expanded constructor avoids converting the factored trace back
    into a polynomial. Starting from an already classified free fourteen-gamma
    word, it takes **256.79 ms** with the frozen production library. A separate
    matched-library experiment improves **248.23 → 157.61 ms** with the
    polynomial-emission patches. These timings include metric construction,
    polynomial assembly, Atom emission and destruction, but exclude public gamma
    recognition and cleanup. They do not replace the complete pipeline timings
    above.

    The direct representation also has mixed downstream effects. With actual
    Spenso metric normalization, the order-sensitive twelve-gamma result shrinks
    **133,818 → 52,401 Atom bytes**, improving its rerun **1.927 → 0.764 ms**.
    The branching-eight result instead grows **1,852 → 3,398 bytes** and its
    rerun slows **29.5 → 52.6 µs**. All **30 downstream comparison processes**
    agree algebraically. These cases rule out expanding every contracted word
    by default; scalar spectators retain their factorization in both routes.

    A subsequent public-library prototype shares the word evaluator between
    factored and sparse output algebras. Its order-sensitive-twelve full expanded
    lifecycle improves **5.995 → 1.472 ms**, against the separate **0.733 ms**
    FORM reference. However, its default first call slows **0.974 → 1.405 ms**,
    while the rerun improves **1.855 → 0.725 ms**.
    The prototype is rejected; the retained timings above remain the production
    checkpoint. These bare traces return directly from the terminal evaluator,
    so an outer cleanup pass cannot explain the first-call regression.

    A frozen-library profile identifies intermediate sparse lists as a material
    cost: the order-sensitive case increases **9.34 → 14.06 million instructions**
    and **3,232 → 11,609 allocations**. Memo-result cloning costs **2.31 million
    instructions**, and sparse finalization costs **6.77 million**. Earlier
    wall-time regressions in unselected controls do not consistently reproduce:
    free lengths twelve and fourteen have identical allocation counts and
    instruction counts within 1%. Duplicate coefficient construction adds only
    **0.03% instructions** in the interior-pair control.

    A scratch follow-up stores shared node handles instead of intermediate
    lists and emits polynomial leaves once, retaining the existing recurrence.
    Its order-sensitive-twelve complete expanded lifecycle measures **1.099 ms**,
    versus **6.276 ms** for the matched factored baseline and **1.570 ms** for
    the sparse-list trial. Instructions fall to **9.56 million**, allocations to
    **3,026**. All **52 trace cases** and **87 boundary/callback files** match
    the validated sparse-list trial exactly. Production was unchanged at that checkpoint; FORM was
    not rerun for this follow-up.

    The retained scalar-power admission fix removes an unnecessary cleanup
    pass. Powers request dot normalization only for a metric-function base;
    rank-one vectors are already detected by the function observer. The scan
    still visits bases and exponents, with conservative handling of opaque
    metadata. A scalar `(x+y)^8` no longer forces Schoonschip over an otherwise
    simplified trace polynomial.

    This matched copied-production checkpoint uses the **default factored API**
    and keeps `S=(x+y)^8` intact. Times are warm milliseconds, including returned
    output destruction: three alternating processes, five samples each.

    | Input | Before | After | Rerun before | Rerun after |
    |:--|--:|--:|--:|--:|
    | S × free trace 8 | 0.3758 | 0.2837 | 0.1452 | 0.0520 |
    | S × free trace 10 | 2.4488 | 1.6274 | 1.2786 | 0.4716 |
    | S × free trace 12 | 23.3361 | 14.6916 | 14.4584 | 5.2234 |

    For twelve gammas, instructions fall **299.24 → 192.20 million** and the
    unnecessary outer Schoonschip walk disappears. Tiny controls remain mixed:
    a scalar variable **0.296 → 0.325 µs**, a scalar sum **0.362 → 0.347 µs**.
    Order-sensitive twelve changes **152.6 → 158.0 ms**, while its rerun changes
    **3.382 → 3.360 ms**. These are workload-specific gains, not FORM parity.

    All **122 power cases × four routes = 488 outputs** match exactly, as do
    **11 production inputs and their reruns**. Compatibility checks also pass
    **52 trace cases, 64 boundary/callback/scope cases, 87 unchanged default
    files and four integration tests** in the experimental expanded API.
    That explicit API was **unintegrated at this checkpoint**: expanding inside a larger
    expression still affected later cleanup costs. The retained power fix is
    independent of it. All **40 targeted repository tests** and Clippy pass. Unit
    fixtures explicitly register their parsed vector namespace and verify a productive power identity; the initial namespace-only
    fixture failure remains recorded. Sources, samples, profiles and validation
    are archived under `retained_scalar_power_admission`.

    The next retained change checks the **whole rebuilt expression** before
    remaining cleanup. If a complete observation proves no contraction,
    normalization or Dirac/epsilon work remains, it returns immediately.
    Otherwise cleanup runs as before; the next iteration reuses the observation
    only if every cleanup result is exactly unchanged. Outer factors and
    callback results remain part of this check.

    After the scalar-power fix, five paired rounds on the final concise
    production source give these default-call milliseconds:

    | Input | Before | After |
    |:--|--:|--:|
    | S × free trace 8 | 0.2844 | 0.2064 |
    | S × free trace 10 | 1.6479 | 0.9773 |
    | S × free trace 12 | 14.3991 | 7.2032 |
    | S × order-sensitive trace 12 | 154.4985 | 154.2121 |

    Free-twelve instructions drop **188.59 → 96.58 million**; its rerun is
    nearly unchanged at **5.379 → 5.179 ms**. An earlier variant added a scan
    on failed admission; carrying unchanged observations removes that duplicate
    work. All **29 final default/disabled/compute fixture modes** and **87
    boundary/callback outputs** match. A separate regression checks a cleanup
    callback that introduces a new trace after the observation. All **29 targeted
    repository Dirac simplification tests** and Clippy pass; this Clippy command
    does not enable warnings-as-errors.

    Small costs remain: disabled free-eight/twelve controls increase **1.8% /
    3.6%**. An untouched Symbolica scalar-expansion control changes **1.7%**,
    with identical instruction counts and an unchanged rerun. This is a large
    free-trace cleanup gain, not a claim of zero overhead or FORM parity.

    The explicit expanded-trace API was **unintegrated at this checkpoint**. With the preceding
    carryover variant on both routes, the same complete expanded free-twelve
    output takes about **50.9 ms** when simplifying first and expanding only
    the trace body, versus **62.7 ms** inside the experimental API. Both preserve
    scalar spectators. The `retained_post_rewrite_completion` record keeps all
    three stages, the untouched-control comparison and the callback regression.

    ### Traces with exact scalar spectators

    The first retained scalar-context dispatch reused the existing terminal evaluator
    for **one ordinary trace times exact function-free scalars**, preserving
    `S=(x+y)^8`. Default-call milliseconds include output destruction:

    | Input | Before | Retained |
    |:--|--:|--:|
    | S × order-sensitive trace 12, D | 159.596 | 2.878 |
    | S × alternating slashes 12, D | 1.636 | 0.0569 |
    | S × paired slashes 12, D | 0.1717 | 0.0163 |
    | S × free trace 12, 4D | 10.995 | 5.493 |
    | S × free trace 12, D | 7.487 | 7.243 |
    | S × free trace 10, D | 0.9855 | 1.0198 |

    The order-sensitive source improves **55.5×**; its rerun changes
    **3.372 → 1.807 ms**. **Expanded-output performance remains below FORM.**
    A separate lifecycle expands only the evaluated trace body and restores S:

    | Trace body | Before, ms | Retained, ms | FORM CPU, ms |
    |:--|--:|--:|--:|
    | Order-sensitive 12, D | 165.193 | 8.104 | 0.7433 |
    | Free 12, D | 55.822 | 53.940 | 3.9333 |

    The Rust matrix uses three alternating process rounds, five samples each,
    on CPU 9 with native **opt-level 2/debug-assertions** libraries. FORM 5.0.0
    reports internal `tracen` plus sorting CPU time after setup; Rust reports
    wall time. The expanded bodies still take about **11× / 14×** the FORM
    reference. Factored default output is a different amount of work.

    Admission requires canonical ordinary words, at most two occurrences per
    explicit index, exact scalar metadata and compact vectors without custom
    normalizers. Axial words, external tensor factors, compound metadata,
    custom vector callbacks and numeric domains outside exact rational/complex-
    rational coefficients keep the established route. Rounded trace units exposed coefficient
    reassociation, and indexed-vector callbacks exposed a different rewrite
    order; the guards preserve their original first and rerun outputs exactly.

    At this checkpoint, the initial shared scan gated one-shot admission. The
    **whole rebuilt product** then had to pass the completion certificate;
    otherwise the same observation entered outer cleanup. The following
    refinement removes that redundant output scan under the same strict
    admission. The explicit expanded-trace API and isolated Symbolica patches
    remained unintegrated.

    Small costs stay explicit: generic-D free-ten is **3.5% slower**, with an
    untouched scalar-expansion control **1.7% slower** and free-ten instruction
    count **0.73% lower**. The wall shift has no established cause. The previous
    9.5% instruction overhead on a 128-metric product is removed; its final
    public time changes **28.36 → 27.77 µs**. Heavy failed admission stays near
    **2.53 ms**.

    Validation passes **20 public cases, 84 boundary cases, 82 trace-body
    polynomial comparisons, 81 default/disabled/compute modes and eight exact
    numeric-dimension status checks**. All **15 raw-versus-cleaned products**
    agree. Independent HEP networks pass **24 case/sample rows**, providing
    **48 comparisons** against original gamma products. All **149 selected
    repository Dirac/Schoonschip tests** and scoped Clippy pass; warnings are
    not treated as errors. Two older assertions now compare exact trace-body
    polynomials while explicitly preserving the spectator; those test
    corrections required no production change. The
    `retained_scalar_context_trace_dispatch` record preserves sources, builds,
    samples, reviews and superseded rounded/admission trials.

    ### Closed terminal products and compact vectors

    The strict scalar-context evaluator now **returns the rebuilt product
    directly**. Each free explicit slot occurs once per trace monomial, and
    the admitted exact scalar spectators cannot add contraction partners.
    Wildcard metadata can retain a diagonal metric, but its eliminated dummy
    has no distinct partner and existing cleanup leaves it inert. Admission,
    callback exclusions and the shared recurrence remain unchanged.

    Matched timings with both changes combined, against the same build before
    either change:

    | Default call | Before, ms | After, ms |
    |:--|--:|--:|
    | S × order-sensitive trace 12, D | 2.907 | 0.964 |
    | S × free trace 10, D | 0.970 | 0.569 |
    | S × free trace 12, D | 7.438 | 2.373 |
    | S × free trace 12, 4D | 5.151 | 2.974 |

    Including expansion of only the trace body gives **8.160 → 6.186 ms** for
    the order-sensitive word and **52.138 → 46.590 ms** for free twelve. FORM
    still takes **0.7433 / 3.9333 ms CPU**; this does not establish parity.
    Rust uses native opt-level 2/debug assertions, three alternating rounds,
    five samples per process on CPU 8, and includes output destruction. The
    isolated closure trial executes **61.2% fewer instructions** on contextual
    free-twelve. Its bare control moved 7.1% slower despite 1.3% fewer
    instructions; the combined run moves about 1% slower. The causes remain
    unestablished. Order-sensitive reruns improve **1.819 → 0.918 ms**;
    generic free-twelve stays near **5.1 ms**. Closure itself does not
    accelerate already simplified reruns.

    The candidate scan also uses the existing strict **SlotMatcher** to recognize
    canonical compact vectors that need no dot normalization. Head observation
    and conservative handling of malformed or opaque shapes are preserved.
    Combined compact-vector products (64 factors) take **22.35 → 17.02 µs**
    through Schoonschip and **33.49 → 18.13 µs** through gamma simplification.
    Malformed fallback grows **0.435 → 0.547 µs** and wrapped compact fallback
    **0.896 → 1.194 µs**. A rejected prefilter recovers those costs but repeats
    decoding and slows successful compact inputs. Paired-slash reruns grow
    **2.465 → 2.858 µs**, disabled free-twelve calls **8.408 → 9.191 µs**,
    and the untouched scalar expansion control is **2.5% slower**. Gains are
    not uniform. Three-vertex networks stay near **9.6 ms** (smallest degree)
    and **9.2 ms** (minimum product terms).

    Closure checks add **392 metadata/symbol cases and 557 wildcard cases**,
    preserving raw, cleaned, contextual and rerun results exactly. Independent
    HEP networks again pass **48 original-gamma comparisons per build**.
    Compact-vector checks preserve **294 mode records**, callback counts and
    **42 public benchmark outputs**. Evidence is archived under
    `retained_terminal_product_closure` and `retained_compact_vector_candidates`.
    The combined build passes **152 selected repository tests**, scoped Clippy,
    **42 compact/gamma/gluon output pairs** and **23 trace lifecycle controls**.
    Its separate measurements are in `retained_closure_compact_combination`.

    ### Expanded trace output: integrated API and measured checkpoints

    The source now includes `GammaSimplifySettings(expand_traces=True)` in
    Rust and Python. It emits the polynomial from the shared trace recurrence
    and existing factored recipe; `S=(x+y)^8` stays outside. Default output
    remains factored. Unsupported sparse shapes use the existing arithmetic
    followed by trace-local expansion, with callback and metadata handling
    preserved. Each independent trace keeps its own expansion boundary.

    The initial matched trial, before strict bare admission and literal-4
    folding, measured these complete calls with S preserved, in milliseconds:

    | Trace body | Simplify then expand body | Prototype | FORM CPU |
    |:--|--:|--:|--:|
    | Order-sensitive 12, D | 6.248 | 0.851 | 0.7533 |
    | Free 10, D | 3.475 | 1.345 | 0.3433 |
    | Free 12, D | 47.868 | 17.153 | 3.9667 |
    | Free 14, D | 871.254 | 418.609 | 53.0 |

    This approaches the order-sensitive reference but **does not establish
    general parity**. Rust reports wall time including output destruction;
    FORM reports trace-plus-sort CPU time with S excluded. Eight fixtures,
    bare and scalar contexts, use three alternating process rounds and five
    samples per mode on CPU 10. Bare free-fourteen improves
    **886.7 → 369.1 ms**. Expanded free-twelve reruns with S improve
    **18.29 → 11.40 ms**, so large reruns remain expensive.

    Exact checks cover defaults, expanded bodies, spectators, rounded boundaries,
    nested callbacks and metadata, including **392 strict leaf and 557 wildcard
    cases**. Nine fresh FORM references retain sorted-word and archived-output
    checks across **27 processes**. Existing HEP certificates remain prior
    independent evidence; they are not new network evaluations.

    The initial trial exposes two concerns: standalone mixed-six grows
    **25.27 → 34.87 µs** because it completes possible callback output, and
    default scalar free-ten has a noisy **24% median increase**. A default-only
    probe finds **5.082 → 5.007 million instructions** but **523 → 556 µs**
    wall medians, with paired ratios **0.978–1.222**. There is no corresponding
    instruction increase, but the wall-time concern remains unexplained.
    A bounded follow-up tries the same strict certificate first for bare traces
    when expansion is requested. Mixed-six improves **34.68 → 24.53 µs**;
    order-sensitive D stays near **0.84 ms**. Exact checks pass across **72
    additional timed processes**. Rejected metadata costs another **0.534 µs**
    and callback/axial controls grow about **4%**. This strict entry path is now
    integrated; the default wall uncertainty remains. Historical trial numbers
    are retained in `expanded_trace_closure_rebase_trial`.

    The integrated sparse backend also folds a **literal exact dimension 4**
    into its integer coefficients. It uses the existing 4D check; rounded,
    wildcard and unsupported dimensions keep their previous routes. The
    matched copied-library comparison below isolates this small change:

    | Expanded trace body | Before, µs | Exact-4 path, µs | FORM CPU, µs |
    |:--|--:|--:|--:|
    | Order-sensitive 12, 4D | 432.27 | 163.82 | 53.33 |
    | Pair with five interior gammas, 4D | 4,086.94 | 1,177.43 | — |
    | Repeated compact 8, 4D | 17.42 | 13.83 | 1.40 |
    | Branching 8, 4D | 23.35 | 21.22 | 2.00 |
    | Alternating 14, 4D | 20.97 | 20.89 | 10.20 |
    | Mixed branching 8, 4D | 24.38 | 24.42 | — |

    With S preserved, order-sensitive twelve takes **443.60 → 172.34 µs**
    and the five-interior case **4.145 → 1.441 ms**. The generic-D control
    stays **827.18 → 827.70 µs** bare and **843.37 → 847.19 µs** with S.
    These are three alternating process pairs, five samples per mode, on
    CPU 10; Rust includes public evaluation, emission and output destruction.
    FORM measures internal `trace4` plus sorting CPU, with setup and S excluded.
    It remains faster. Mixed cases and reruns do not uniformly improve.

    All **84 timed processes**, the reused trace/callback/wildcard matrix,
    **29 additional numeric-domain cases**, and five public Rust tests pass
    against the frozen candidate. The integrated build separately passes
    **158 selected Idenso tests**, **four installed expansion tests**, the
    pipeline check, scoped Clippy and generated-stub check.
    Notebook algebra passes; rich rendering still encounters the Symbolica
    export-version mismatch described above.
    The benchmark preserves exact body equality, reruns and scalar factors;
    it adds no fresh HEP component evaluation.

    **Installed Python follow-up, 2026-09-25:** unchanged sums retain their
    normalized storage; exact no-op transforms reuse validated metadata in a
    fresh object. Multiplicity checks share the existing slot matcher, collect
    all index counts together, and reuse fully explicit built-in leaf interfaces.
    The last two changes have this separate matched comparison, in milliseconds:

    | Expanded body | First call before | First call after | Rerun after | FORM body CPU |
    |:--|--:|--:|--:|--:|
    | Free 2, D | 0.0272 | 0.0262 | 0.00103 | 0.000640 |
    | Free 4, D | 0.1069 | 0.1023 | 0.00183 | 0.001350 |
    | Free 6, D | 0.662 | 0.437 | 0.00738 | 0.005400 |
    | Free 10, D | 70.408 | 30.332 | 0.679 | 0.353333 |
    | Free 12, D | 987.202 | 401.789 | 9.074 | 3.966667 |
    | Free 8, 4D | 5.824 | 2.765 | 0.0591 | 0.048000 |
    | Repeated 8, 4D | 0.163 | 0.133 | 0.00228 | 0.002000 |
    | Repeated 8, D | 1.179 | 0.700 | 0.0126 | 0.014000 |

    Three interleaved process pairs use five samples each, including result
    disposal and retaining S. Setup, references and rendering are excluded.
    These release Python timings are distinct from the native Rust and FORM
    tables. The baseline already has shared slot matching. All eight rerun
    controls remain within **2%**. Changed results retain validation: products
    add explicit multiplicities, sums take their per-index maximum, and the
    interface cache excludes generic functions and unresolved ports.
    Instruction counts fall **35–59%** for free six, ten and twelve. Interface
    inference still takes **39.8%** of the twelve-gamma call, and tensor-power
    lowering takes another **18.3%**. Overall, **97.25%** enters Python result
    construction; those owner percentages overlap. Fresh FORM measurements
    cover the same eight index words across 24 processes with exact polynomial
    checks. FORM times only tracing and sorting the body; it excludes S and
    Python interface construction. The observed ratios therefore compare
    different API boundaries, rather than isolating backend speedups.
    The earlier exact-result
    metadata change reduced twelve-gamma reruns **1,375 → 9.22 ms**; first calls
    still validate changed output. **445,868 differential count queries**,
    **36,864 representation-key checks**, **89 selected Rust tests**, **six
    installed Python checks** and scoped Clippy pass. One wrapped-index Rust
    test fails identically on the baseline; its expectation is unchanged.
    Raw measurements and source identities are in
    `retained_explicit_multiplicity_and_interface_reuse` in the contraction
    parity record.

    **Remaining native gaps:** latest fully expanded native checkpoints are
    0.828 ms versus FORM 0.753 ms for the contracted twelve-gamma D word,
    16.707 versus 3.967 ms for free twelve, and 369.1 versus 53.0 ms for free
    fourteen. Metric chains of length 64 take 23.44 versus 3.20 µs. These are
    earlier matched-output checkpoints, not a new uniform rerun.

    In the fresh three-gluon experiment, the live minimum-product-terms route
    takes 9.221 ms including final scalar emission; a validated but unintegrated
    substitution/dot-fusion candidate takes 8.734 ms. Expand-first remains
    faster on this fixture at about 5.26 ms. FORM's ordered vertex substitution
    and sorting takes 0.325 ms CPU, starting from opaque vertices rather than
    substituted rules. Initial partial parsing is only about 0.064 ms.
    The contraction parity record's `end_to_end_gap_checkpoint` collects all
    available native rows, the fresh Python comparison, and their limitations.

    Late external metrics also contract through factored trace sums. A
    compatible tensor can be carried through a sum when every branch can absorb
    it, including metrics inside the sum with an epsilon outside. Tagged vectors
    follow the same route when vector contraction is enabled; scalar spectators
    retain their factorization. Independent HEP checks cover ordinary, axial and
    symbolic-D traces simplified before the external metric is attached.

    An earlier broad trace/scan checkpoint recorded **363 passing Idenso tests**
    and **46 passing HEP integration tests**. Three Idenso processes initially
    hit the host's open-file limit; the targeted retry passed. The existing tensor-display snapshot
    failure remains, with 23 tests skipped across the two suites. Both scoped
    Clippy checks at that earlier checkpoint passed with warnings denied. The
    Spenso follow-up passes 19 relevant tests; its remaining dual-wrapper filter test fails identically
    with the matcher change removed. All final measured trace outputs match the
    independent certificate.

    Raw samples, generated FORM programs, source hashes, validation and the
    standalone generic-D driver are preserved in
    `examples/notebooks/tensor_contraction_parity.json`.
    """)
    return


@app.cell
def _(E, S, TensorExpression, TensorName, form_executable, g5, lorentz, slash, tr):
    from itertools import permutations

    from symbolica import Expression
    from symbolica.community.spenso import (
        Tensor,
        TensorLibrary,
        TensorNetwork,
    )

    _names = [TensorName.vector(f"gamma_hep::p{i}") for i in range(10)]
    _metric = TensorName.g().to_expression()
    _momenta = [name.to_expression()(lorentz.to_expression()) for name in _names]
    _base = [
        [2, 1, 0, 1],
        [1, 2, 1, -1],
        [3, -1, 2, 1],
        [1, 0, -1, 2],
        [2, -1, -2, 1],
        [-1, 2, 1, 2],
        [3, 1, -1, 0],
        [1, 1, 2, -2],
        [-2, 1, 3, -1],
        [2, 2, -1, 3],
    ]
    _epsilon_name = TensorName("spenso::epsilon", is_antisymmetric=True)
    _epsilon = Tensor.sparse(
        _epsilon_name(lorentz, lorentz, lorentz, lorentz), Expression
    )
    for _perm in permutations(range(4)):
        _inversions = sum(
            _perm[i] > _perm[j] for i in range(4) for j in range(i + 1, 4)
        )
        _epsilon[_perm] = E("-1𝑖" if _inversions % 2 == 0 else "1𝑖")

    _expressions = {}
    for _length in (4, 8, 10):
        for _axial in (False, True):
            _factors = [slash(p) for p in _momenta[:_length]]
            if _axial:
                _factors.insert(0, g5)
            _original = tr(*_factors).to_expression()
            _expressions[_length, _axial] = {
                "Original gamma network": _original,
                "Idenso metric/epsilon network": TensorExpression(_original)
                .simplify_gamma()
                .to_expression(),
            }

    hep_component_form_sources = {}
    if form_executable:
        import re as _re
        import subprocess as _subprocess
        import tempfile as _tempfile
        from pathlib import Path as _Path

        with _tempfile.TemporaryDirectory(prefix="gamma-hep-form-") as _directory:
            _path = _Path(_directory) / "components.frm"
            for _length in (4, 8, 10):
                _word = ",".join(f"p{i}" for i in range(_length))
                for _mode in ("trace4", "tracen"):
                    _source = (
                        f"Off Statistics;\nVectors {_word};\nLocal F=g_(1,{_word});\n"
                        f'{_mode},1;\n.sort\n#write "RESULT=%E",F\n.end\n'
                    )
                    _path.write_text(_source)
                    _run = _subprocess.run(
                        [form_executable, "-q", str(_path)],
                        cwd=_directory,
                        check=True,
                        capture_output=True,
                        text=True,
                        timeout=30,
                    )
                    _polynomial = (
                        _run.stdout.split("RESULT=", 1)[1].strip().removesuffix(";")
                    )
                    # Import each FORM scalar product as a Spenso metric tensor.
                    # Symbolica parses the polynomial; no Python eval is used.
                    _polynomial = _re.sub(
                        r"p(\d+)\.p(\d+)",
                        lambda match: f"gamma_hep::dot{match[1]}x{match[2]}",
                        _polynomial,
                    )
                    _expression = E(_polynomial)
                    for _i in range(_length):
                        for _j in range(_i, _length):
                            _expression = _expression.replace(
                                S(f"gamma_hep::dot{_i}x{_j}"),
                                _metric(_momenta[_i], _momenta[_j]),
                            )
                    _expressions[_length, False][f"FORM {_mode} metric network"] = (
                        _expression
                    )
                    hep_component_form_sources[f"{_length} / {_mode}"] = _source

    hep_component_checks = []
    for _sample in range(2):
        _library = TensorLibrary.hep_lib_atom()
        _library.register(_epsilon)
        for _name, _components in zip(_names, _base, strict=True):
            # Sparse Expression storage forces exact components; automatic
            # dense input conversion can otherwise choose floating-point data.
            _tensor = Tensor.sparse(_name(lorentz), Expression)
            for _axis, _value in enumerate(_components):
                _tensor[_axis] = E(
                    str(_value if _sample == 0 else _value * (1, -1, 2, 1)[_axis])
                )
            _library.register(_tensor)

        def _evaluate(_expression, _library=_library):
            _network = TensorNetwork(_expression, library=_library)
            _network.execute(library=_library)
            return _network.result_scalar()

        assert _evaluate(_epsilon_name.to_expression()(*_momenta[:4])) != E("0")
        for (_length, _axial), _routes in _expressions.items():
            _expected = _evaluate(_routes["Original gamma network"])
            for _route, _expression in _routes.items():
                _value = _evaluate(_expression)
                assert _value == _expected, (_sample, _length, _axial, _route)
                hep_component_checks.append(
                    {
                        "sample": _sample + 1,
                        "gammas": _length,
                        "gamma5": _axial,
                        "route": _route,
                        "exact HEP value": _value.format_plain(),
                        "agrees": True,
                    }
                )
    return hep_component_checks, hep_component_form_sources


@app.cell(hide_code=True)
def _(hep_component_checks, hep_component_form_sources, mo):
    mo.vstack(
        [
            mo.ui.table(
                hep_component_checks,
                selection=None,
                pagination=False,
                show_download=False,
            ),
            mo.md(
                "FORM rows are included only when the executable is available; the HEP checks always run."
            ),
            mo.accordion(
                {
                    f"FORM component check: {name}": mo.md(f"```form\n{source}\n```")
                    for name, source in hep_component_form_sources.items()
                }
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Ordered gluonic-ladder Feynman rules

    This benchmark preserves the supplied eight-vertex routing, the six-term
    three-gluon Lorentz rule, and the substitution order **1 → 8**. It contracts
    the two external Lorentz ports with `k10` and `k20`. No on-shell identities,
    color factors, couplings, or propagator denominators are included.
    The supplied topology has eleven internal edges and eight vertices, hence
    **four loops**, despite the original three-loop label.

    The reference FORM script sorts and contracts after *each* replacement.
    Its intermediate term counts are **6, 34, 192, 1,084, 6,069, 10,947,
    23,937, 9,652**. Changing substitution order changes the intermediate work.

    Two timings answer different questions. **Rule substitution** applies the
    eight replacements with `rhs_cache_size=1000` and a held RHS expansion;
    the vertex product stays factored. **Ordered scalar reduction** additionally
    expands and contracts the accumulated vertices after every replacement.
    Unprocessed vertex factors are held outside that accumulator; they multiply
    every intermediate term and do not change its term count.
    Only the latter has the same 9,652-term endpoint as FORM. This explicitly
    requested polynomial expansion is confined to this standalone benchmark.

    **Recorded validation, 2026-09-24:** an optimized native Spenso build took
    **54.970 s** for the full reduction. FORM 5.0.0 took **0.744 s** warm process
    time (three-run median; **0.72 s** reported CPU time), a wall-time ratio of
    **73.9×**. Rule substitution alone took **0.424 ms** warm (five-run median).
    All eight term counts and every coefficient of the final polynomial agreed.
    These are shared-host measurements, recorded with build information in
    `examples/notebooks/gluon_ladder_timing.json`; the table below measures the
    current session when the full-reduction button is pressed.

    **Refreshed complete route, 2026-09-25:** the current optimized
    Python/TensorExpression build takes **46.474 s** in one complete run.
    Fresh FORM takes **0.742 s** warm process wall (three measured runs), or
    **0.71 s** internal CPU: an observed wall gap of **62.7×**. All eight
    intermediate counts and all final coefficients agree. The last two vertex
    stages take **13.311 s** and **25.147 s**. This uses the same notebook route
    and includes result wrapping; it is separate from the three-vertex
    scoped-substitution experiment below. Full records are retained under
    `full_gluon_ladder_refresh` in the contraction parity archive.

    **Fast constructor inference:** the constructor now shares Spenso's cached
    slot matcher, borrows Atom views, and reads ordinary compact-vector ports
    directly, avoiding temporary dummy-index materialization. All-summand
    validation and custom normalization behavior are retained. Constructor-only
    medians (three calls on each saved, parsed input) are:

    | Input terms | Before | After | Speedup |
    | ---: | ---: | ---: | ---: |
    | 6 | 0.128 ms | 0.095 ms | 1.35× |
    | 36 | 1.230 ms | 0.771 ms | 1.60× |
    | 204 | 8.213 ms | 5.140 ms | 1.60× |
    | 1,152 | 57.788 ms | 35.181 ms | 1.64× |
    | 6,503 | 419.630 ms | 241.470 ms | 1.74× |
    | 50,684 | 3.922 s | 2.242 s | 1.75× |
    | 92,341 | 8.151 s | 4.451 s | 1.83× |
    | 186,516 | 18.917 s | 9.875 s | 1.92× |

    The full ladder's operation sum improves **47.506 → 25.340 s (1.88×)**:
    constructor calls take **36.420 → 16.941 s**, and Schoonschip including
    result wrapping takes **10.291 → 7.648 s**. This harness serializes fixtures
    between operation clocks; its sum is separate from the contiguous run above.
    Fresh FORM takes **0.727 s** process wall (three-run median; **0.70 s** CPU),
    leaving an observed **34.8×** ratio. These are shared-host observations.
    All intermediate Atoms/interfaces and the FORM-certified final polynomial
    agree, as do the **68** adversarial constructor cases and **11** metric/trace
    controls. The `retained_fast_constructor_inference` archive entry retains
    the sources, builds, timings, validation and measurement limits.

    Five-call constructor medians on separate controls improve a 32-metric chain
    **458.7 → 378.4 µs**, compact dot sums **45.2 → 24.2 µs**, and a free-six
    D-dimensional trace **78.5 → 67.7 µs**. At that checkpoint, a native profile
    showed why large
    constructors remain slow: interface/syntax analysis owns **44.7%** of sampled
    cycles, product normalization **20.3%**, multiplicity/placeholder/nesting
    checks **25.3%**, and power lowering **9.4%**. These are sample proportions,
    not exact wall-time phases. Leaf-port extraction is only part of construction.

    **Borrowed product storage and scalar-interface reuse:** the next retained
    change keeps unchanged products/powers in their existing Atom storage and
    reuses inferred scalar dots when their operands are already known to be
    reusable. A fresh before/after comparison improves the full operation sum
    **25.160 → 17.740 s (1.42×)**, with constructors **16.851 → 10.155 s** and
    Schoonschip including wrapping **7.547 → 6.801 s**. The largest standalone
    constructor improves **10.362 → 7.134 s** (three-call medians). Fresh FORM
    takes **0.719 s** process wall, leaving a **24.7×** ratio with the phase-sum
    caveat above. All intermediate Atoms/interfaces and final coefficients agree.

    Complete first trace transformations also improve:

    | Expanded trace | Python before | Python after | FORM body CPU |
    | --- | ---: | ---: | ---: |
    | Free 6, D | 0.398 ms | 0.342 ms | 0.0054 ms |
    | Free 10, D | 29.130 ms | 24.011 ms | 0.3433 ms |
    | Free 12, D | 387.388 ms | 316.376 ms | 4.0667 ms |
    | Free 8, 4D | 2.608 ms | 2.133 ms | 0.0480 ms |
    | Repeated 8, 4D | 0.119 ms | 0.103 ms | 0.0020 ms |
    | Repeated 8, D | 0.651 ms | 0.543 ms | 0.0140 ms |

    Python includes a factored scalar spectator and result construction; FORM
    excludes the spectator. These are different timing boundaries. All **24**
    fresh FORM processes pass exact polynomial checks. Tiny constructor controls
    are mixed: a 32-metric chain improves **378 → 268 µs**, while compact dots
    move **24.2 → 25.8 µs**. Reruns are mostly unchanged. Full records are under
    `retained_unchanged_products_and_scalar_interfaces` in the parity archive.

    **Retaining the interface through expansion:** ordinary changed, nonzero
    typed expansion now validates the smaller factored input and carries its
    ordered interface when the tensor leaves permit this. Callback-sensitive
    leaves, exposed tensor powers, unresolved ports, and other expansion modes
    keep their existing output checks. Unchanged results avoid the extra scan.

    In a baseline/candidate/baseline comparison, the complete typed ladder takes
    **27.740 / 20.358 / 28.684 s**. Its initial expansion phase falls from
    **8.687–9.250 s to 1.849 s (4.70–5.00×)**; total improvement is
    **1.36–1.41×**. All stage Atoms, logical interfaces, and final coefficients
    agree. The **105** expansion controls and **68** constructor cases also
    agree, including rejected inputs and callback behavior. Fresh FORM takes
    **0.751 s** process wall: the typed operation sum remains **27.1×** larger,
    with different timing boundaries. No trace-kernel improvement is claimed.

    This is a typed-route optimization. The raw control moves **21.999 →
    18.253 s** without using this shortcut; that movement is not attributed to
    the change. The earlier **17.740 s** raw result remains a separate checkpoint.
    Multiplication still takes **10.603 s** in the candidate typed ladder. A
    focused profile of its final multiplication assigns **90.3%** of sampled
    cycles to port rewriting, chiefly rebuilding the sum even when indices
    already match. These are cycle samples, not a wall-time partition.

    `retained_factored_expansion_interface` records the measured release,
    controls, profile, and timing limits. The concurrent broader API audit
    modified the wrapper afterward; its newer source is outside this measured
    checkpoint. Those measurements used the then-validated notebook release.

    **Shared symbolic tensor machinery:** Idenso's `SymbolicTensor<S>` now owns
    both explicit network tensors and tensors with a `PartialStructure` that
    retains logical port order and unresolved occurrences. The duplicate
    `StructuredAtom` is removed. Composition, contraction, index substitution,
    fast inference, and checked result reconstruction use this shared owner;
    Spynso handles Python conversion and result wrapping.

    Identity substitutions borrow unchanged branches while still checking for
    the requested index in every sum branch. Genuine relabeling uses bulk
    sum/product construction when exact arithmetic and callback behavior permit.
    Typed zeros keep their interfaces. Rewrites that can invoke a normalizer
    validate the resulting interface: `g(a,b)*T(a)` must not retain a vector
    interface if the callback turns `T(b)` into a scalar.
    Result checks observe the encoded ports without replaying a normalizer on
    synthetic indices. Unchanged results reuse the already-owned result Atom
    and are wrapped directly, preserving metadata without another tensor copy.

    **Partial-network follow-up:** a closed subcase with vertices **1, 2, 8**
    produces **64 scalar terms**. The baseline expand-first route takes about
    **5.56 ms**; network contraction with local sum distribution takes
    **24.7–27.3 ms**, with final polynomial emission separately **0.30–0.48 ms**.
    Initial depth-one parsing is only **0.064 ms**, under **0.3%** of the network
    route. Independently warmed phases are not an exact additive lifecycle.
    Faster factored-only outputs still have uncontracted indices.

    All **72 exact HEP component comparisons** and **nine FORM polynomial
    comparisons** pass. FORM orders **1 → 2 → 8** and **8 → 1 → 2** pass through
    **6 → 34 → 64** and **9 → 45 → 64** terms, taking about **328 µs** and
    **270 µs CPU** respectively. FORM includes rule substitutions and sorting;
    the Rust wall timings start from normalized six-term rules, so these are
    different timing boundaries. Term counts alone do not predict the winner.

    **Retained scalar cleanup improvement:** certified index-free scalar leaves
    now skip recursive network construction, retaining numeric-coefficient
    distribution and the existing handling of parser-owned syntax. At the
    scalar-shortcut-only checkpoint, before residual-query deduplication, a
    matched five-process-pair comparison improves the three-vertex route from
    **23.70 to 13.62 ms**; minimum-product-terms ordering improves from
    **26.63 to 14.08 ms**. The expand-first control stays **5.50–5.53 ms**.
    One-vertex, two-vertex and contracted-sum controls also improve. Parse/merge
    passes drop **487 → 77**, with the same **16 contraction entries**: **four
    consume ports** and **twelve rebuild plain products**. All **3,200** output/error
    comparisons and **39** repository Schoonschip tests pass.

    **Retained residual-query improvement:** deduplicating identical queries
    before running the existing matcher cuts residual-scan instructions by
    about **49%** in both three-vertex ordering profiles. Whole-route instruction
    medians fall **4.9%** with smallest-degree ordering and **2.9%** with
    minimum-product-terms ordering. All **3,712 exact output/error comparisons**
    agree. The wall run is invalid for a speedup claim because unchanged controls
    became strongly bimodal under concurrent host load. The scalar-only timings
    above are a separate checkpoint; no additional wall-time gain is claimed.

    **Retained scalar network-entry shortcut:** the same certificate now avoids
    reconstructing depth-one networks for scalar terms, preserving sum-branch
    scheduling and coefficient cleanup. This matched checkpoint starts after the
    retained scalar-leaf, unique-query and scalar-power changes:

    | Input | Before entry shortcut | After entry shortcut |
    |---|---:|---:|
    | One closed vertex | 0.698 ms | 0.440 ms |
    | Two closed vertices | 2.986 ms | 2.793 ms |
    | Three closed vertices | 13.432 ms | 10.403 ms |
    | Two contracted sums | 0.188 ms | 0.129 ms |

    These are smallest-degree route medians from five paired processes. For three
    vertices, minimum-product-terms ordering improves **13.695 → 9.895 ms**.
    Unchanged expand-first medians stay within **2%**. The clocks include result
    destruction and exclude input parsing, initial normalization and final
    polynomial emission. This checkpoint is not additive to earlier timings or
    the separate initialization-gate trial.

    Parse counts fall **77 → 19** and **117 → 19**. Their original decomposition
    is `77 = 1 + 4*6 + 2*26` and `117 = 1 + 6 + 6 + 6 + 26 + 2*36`; the new route
    keeps `19 = 1 + 3*6`. The same **four port-consuming contractions** and
    **twelve plain-product entries** remain. Instructions fall **19.9%** and
    **23.9%**. Validation matches **16,960 output/error records** plus **192
    caught-panic statuses** from an invalid nonuniform sum boundary at deeper
    parsing, including direct SinglePass and compact-root/callback controls.
    All **41 repository Schoonschip tests** and the scoped Clippy check pass.

    **Retained duplicate-normalization cleanup:** an unchanged scalar term now
    uses its already normalized result directly in the existing coefficient
    cleanup. Parser-bound inputs keep their pre-parse and final normalization.
    A matched five-process-pair checkpoint after the entry shortcut gives:

    | Input and order | Before | After |
    |---|---:|---:|
    | Three vertices, smallest degree | 10.212 ms | 10.022 ms |
    | Three vertices, minimum product terms | 10.033 ms | 9.485 ms |
    | 16 compact dots | 18.60 µs | 16.35 µs |
    | 64 compact dots | 68.62 µs | 59.51 µs |

    The smallest-degree gain is modest. A contracted-sum control is **1.8%
    slower** with overlapping process ranges; its instruction median improves.
    Unchanged expand-first controls stay within **0.7%**. Three-vertex
    instruction medians decrease **4.0%/6.5%**; compact-dot controls, including a
    factored scalar spectator, decrease **10–22%**. All **20,032 output/error
    records** and **192 caught-panic statuses** agree. All **149 selected Dirac and
    Schoonschip tests** and the scoped Clippy check pass. Route clocks exclude final
    polynomial emission and use the frozen dependency family without the separate
    initialization gate; these measurements do not establish FORM parity.

    Direct sum contraction remains a substantial cost. The existing opaque
    tensor boundaries can guide local distribution and directed substitution.
    Sum boundaries must be validated: fast inference reads only the first branch.
    Fixtures, timing scopes and validation are recorded under
    `gluon_network_planning_followup` in `tensor_contraction_parity.json`.
    """)
    return


@app.cell
def _(E, S, TensorName, lorentz, prod):
    from symbolica import T

    ladder_momenta = {
        i: TensorName.vector(f"gluon_ladder::k{i}").to_expression()(
            lorentz.to_expression()
        )
        for i in (0, 1, 2, 3, 4, 10, 20)
    }
    _k = ladder_momenta
    _mu = {i: lorentz(S(f"gluon_ladder::mu{i}")).to_expression() for i in range(1, 12)}
    ladder_metric = TensorName.g().to_expression()
    _vx = S("gluon_ladder::vx")
    # All three momenta enter each vertex. The final three arguments are
    # Lorentz slots, or compact momenta for the two external polarizations.
    _routing = [
        (-_k[0], _k[0] - _k[1], _k[1], _k[10], _mu[1], _mu[8]),
        (-_k[1], _k[2], _k[1] - _k[2], _mu[1], _mu[2], _mu[9]),
        (-_k[2], _k[3], _k[2] - _k[3], _mu[2], _mu[3], _mu[10]),
        (-_k[3], _k[4], _k[3] - _k[4], _mu[3], _mu[4], _mu[11]),
        (-_k[4], _k[0], _k[4] - _k[0], _mu[4], _k[20], _mu[5]),
        (-_k[4] + _k[0], -_k[3] + _k[4], _k[3] - _k[0], _mu[5], _mu[11], _mu[6]),
        (-_k[3] + _k[0], -_k[2] + _k[3], _k[2] - _k[0], _mu[6], _mu[10], _mu[7]),
        (-_k[2] + _k[0], -_k[1] + _k[2], _k[1] - _k[0], _mu[7], _mu[9], _mu[8]),
    ]
    assert all(sum(vertex[:3], E("0")) == E("0") for vertex in _routing)
    ladder_vertices = [_vx(i, *vertex) for i, vertex in enumerate(_routing, 1)]
    ladder_input = prod(ladder_vertices)
    _p1, _p2, _p3, _i1, _i2, _i3 = S(
        *(f"gluon_ladder::{name}_" for name in ("p1", "p2", "p3", "i1", "i2", "i3"))
    )
    _g = ladder_metric
    ladder_rhs = (
        -_g(_i1, _i3) * _g(_p1, _i2)
        + _g(_i1, _i2) * _g(_p1, _i3)
        + _g(_i2, _i3) * _g(_p2, _i1)
        - _g(_i1, _i2) * _g(_p2, _i3)
        - _g(_i2, _i3) * _g(_p3, _i1)
        + _g(_i1, _i3) * _g(_p3, _i2)
    ).hold(T().expand())
    ladder_patterns = [_vx(i, _p1, _p2, _p3, _i1, _i2, _i3) for i in range(1, 9)]
    ladder_expected_counts = (6, 34, 192, 1084, 6069, 10947, 23937, 9652)
    return (
        ladder_expected_counts,
        ladder_input,
        ladder_metric,
        ladder_momenta,
        ladder_patterns,
        ladder_rhs,
        ladder_vertices,
    )


@app.cell
def _(ladder_input, ladder_patterns, ladder_rhs, median, perf_counter):
    ladder_substitution_samples = []
    for _round in range(6):
        _result = ladder_input
        _start = perf_counter()
        for _pattern in ladder_patterns:
            _result = _result.replace(_pattern, ladder_rhs, rhs_cache_size=1000)
        ladder_substitution_samples.append(perf_counter() - _start)
        if _round == 0:
            ladder_factored = _result
        else:
            assert _result == ladder_factored
    ladder_substitution_seconds = median(ladder_substitution_samples[1:])
    return ladder_substitution_samples, ladder_substitution_seconds


@app.cell
def _(mo):
    run_ladder_reduction = mo.ui.run_button(
        label="Time the full ordered ladder reduction (can take minutes)"
    )
    run_ladder_reduction
    return (run_ladder_reduction,)


@app.cell
def _(
    E,
    TensorExpression,
    ladder_expected_counts,
    ladder_patterns,
    ladder_rhs,
    ladder_vertices,
    mo,
    perf_counter,
    run_ladder_reduction,
):
    ladder_stages = []
    ladder_result = None
    if run_ladder_reduction.value:
        _result = E("1")
        _cumulative = 0.0
        with mo.status.progress_bar(
            total=8, title="Applying ladder vertices"
        ) as _progress:
            for _i, (_vertex, _pattern) in enumerate(
                zip(ladder_vertices, ladder_patterns), 1
            ):
                _start = perf_counter()
                _factor = _vertex.replace(_pattern, ladder_rhs, rhs_cache_size=1000)
                # Unprocessed vertices remain spectator factors, just as in FORM.
                # Distribute the current product before contraction, then collect
                # the scalar products and expand the resulting polynomial.
                _result = (
                    TensorExpression((_result * _factor).expand())
                    .schoonschip()
                    .normalize_dots()
                    .to_expression()
                    .expand()
                )
                _seconds = perf_counter() - _start
                _cumulative += _seconds
                _terms = sum(1 for _ in _result.terms())
                assert _terms == ladder_expected_counts[_i - 1], (_i, _terms)
                ladder_stages.append(
                    {
                        "vertex": _i,
                        "Spenso terms": _terms,
                        "Spenso step wall (s)": _seconds,
                        "Spenso cumulative wall (s)": _cumulative,
                    }
                )
                _progress.update()
        ladder_result = _result
    return ladder_result, ladder_stages


@app.cell
def _():
    # Supplied gluon.frm, preserving its rules, routing and sort order.
    ladder_form_source = r"""#-
format nospaces;
Dimension 4;
Auto Index mu;
Auto Vector k;
CF vx, f;
S x;

L F = vx(1,-k0, k0-k1, k1, k10, mu1, mu8)*
    vx(2,-k1, k2, k1-k2, mu1, mu2, mu9)*
    vx(3,-k2, k3, k2-k3, mu2, mu3, mu10)*
    vx(4,-k3, k4, k3-k4, mu3, mu4, mu11)*
    vx(5,-k4, k0, k4-k0, mu4, k20, mu5)*
    vx(6,-k4+k0, -k3+k4, k3-k0, mu5, mu11, mu6)*
    vx(7,-k3+k0, -k2+k3, k2-k0, mu6, mu10, mu7)*
    vx(8,-k2+k0, -k1+k2, k1-k0, mu7, mu9, mu8);

#do i=1,8
    id vx(`i', k1?, k2?, k3?, mu1?, mu2?, mu3?) =
        (- d_(mu1, mu3) * d_(k1, mu2)
                    + d_(mu1, mu2) * d_(k1, mu3)
                    + d_(mu2, mu3) * d_(k2, mu1)
                    - d_(mu1, mu2) * d_(k2, mu3)
                    - d_(mu2, mu3) * d_(k3, mu1)
                    + d_(mu1, mu3) * d_(k3, mu2)
                    );

*    Print +s;
    .sort:gluon-`i';
#enddo

*Print +s;
.end
"""
    return (ladder_form_source,)


@app.cell
def _(
    form_executable, ladder_expected_counts, ladder_form_source, median, perf_counter
):
    import re as _re
    import subprocess as _subprocess
    import tempfile as _tempfile
    from pathlib import Path as _Path

    ladder_form_stages = []
    ladder_form_process_samples = []
    ladder_form_polynomial = None
    ladder_form_stdout = None
    if form_executable:
        with _tempfile.TemporaryDirectory(prefix="gluon-ladder-form-") as _directory:
            _path = _Path(_directory) / "gluon.frm"
            _path.write_text(ladder_form_source)
            for _round in range(4):
                _start = perf_counter()
                _run = _subprocess.run(
                    [form_executable, "-q", str(_path)],
                    cwd=_directory,
                    check=False,
                    capture_output=True,
                    text=True,
                    timeout=120,
                )
                _seconds = perf_counter() - _start
                if _run.returncode:
                    raise RuntimeError(
                        f"FORM ladder failed:\n{_run.stdout}\n{_run.stderr}"
                    )
                _stages = _re.findall(
                    r"Time =\s*([\d.]+) sec\s+Generated terms =\s*(\d+)"
                    r"\s+F\s+Terms in output =\s*(\d+)\s+gluon-(\d+)"
                    r"\s+Bytes used\s*=\s*(\d+)",
                    _run.stdout,
                )
                assert (
                    tuple(int(stage[2]) for stage in _stages) == ladder_expected_counts
                ), _run.stdout
                assert [int(stage[3]) for stage in _stages] == list(range(1, 9))
                if _round:
                    ladder_form_process_samples.append(_seconds)
                    ladder_form_stages.append(_stages)
            ladder_form_stdout = _run.stdout
            ladder_form_stages = [
                {
                    "vertex": _i + 1,
                    "FORM terms": ladder_expected_counts[_i],
                    "FORM generated terms": int(ladder_form_stages[0][_i][1]),
                    "FORM cumulative CPU (s)": median(
                        float(run[_i][0]) for run in ladder_form_stages
                    ),
                }
                for _i in range(8)
            ]
            # A separate, untimed run exports the exact polynomial for validation.
            _path.write_text(
                ladder_form_source.replace(
                    ".end", '#write <polynomial.txt> "%E",F\n.end'
                )
            )
            _check = _subprocess.run(
                [form_executable, "-q", str(_path)],
                cwd=_directory,
                check=False,
                capture_output=True,
                text=True,
                timeout=120,
            )
            if _check.returncode:
                raise RuntimeError(
                    f"FORM polynomial export failed:\n{_check.stdout}\n{_check.stderr}"
                )
            ladder_form_polynomial = (_Path(_directory) / "polynomial.txt").read_text()
    return (
        ladder_form_polynomial,
        ladder_form_process_samples,
        ladder_form_stages,
        ladder_form_stdout,
    )


@app.cell
def _(E, S, ladder_form_polynomial, ladder_metric, ladder_momenta, ladder_result):
    import re as _re

    ladder_exact_match = None
    if ladder_result is not None and ladder_form_polynomial is not None:
        # Compare complete polynomials in independent scalar products. Term
        # counts alone cannot certify signs, coefficients or momentum routing.
        _scalar = ladder_result
        for _i, _p in ladder_momenta.items():
            for _j, _q in ladder_momenta.items():
                if _i <= _j:
                    _scalar = _scalar.replace(
                        ladder_metric(_p, _q), S(f"gluon_ladder::s{_i}x{_j}")
                    )
        _reference = E(
            _re.sub(
                r"k(\d+)\.k(\d+)",
                lambda match: "gluon_ladder::s{}x{}".format(
                    *sorted(map(int, match.groups()))
                ),
                ladder_form_polynomial.strip().rstrip(";"),
            )
        )
        ladder_exact_match = _scalar == _reference
        assert ladder_exact_match, "Spenso and FORM ladder polynomials differ"
    return (ladder_exact_match,)


@app.cell(hide_code=True)
def _(
    form_executable,
    form_version,
    ladder_exact_match,
    ladder_expected_counts,
    ladder_form_process_samples,
    ladder_form_source,
    ladder_form_stages,
    ladder_form_stdout,
    ladder_stages,
    ladder_substitution_samples,
    ladder_substitution_seconds,
    median,
    mo,
):
    _rows = []
    _reference_cpu = (0.00, 0.00, 0.00, 0.00, 0.02, 0.17, 0.36, 0.49)
    for _i, _terms in enumerate(ladder_expected_counts):
        _row = {
            "vertex": _i + 1,
            "reference terms": _terms,
            "supplied FORM cumulative CPU (s)": _reference_cpu[_i],
        }
        if ladder_form_stages:
            _row.update(ladder_form_stages[_i])
        if ladder_stages:
            _row.update(ladder_stages[_i])
        _rows.append(_row)
    _notes = [
        (
            f"**Rule substitution only:** first call {1000 * ladder_substitution_samples[0]:.3f} ms; "
            f"five-run warm median {1000 * ladder_substitution_seconds:.3f} ms. "
            "This result remains factored and is not the 9,652-term scalar polynomial."
        ),
        (
            "Stage times exclude input construction, term counting, validation and display. "
            "Spenso reports a single full run's wall time; FORM statistics report cumulative "
            "CPU time rounded to 0.01 s. The supplied 0.49 s is a reference measurement, "
            "not a measurement of this notebook's machine. Use an optimized community build "
            "for performance comparisons."
        ),
    ]
    if form_executable:
        _notes.append(
            f"**Local FORM:** `{form_version}`. Warm process median "
            f"{median(ladder_form_process_samples):.3f} s over three measured runs, "
            "including startup, parsing and sorting; polynomial export is untimed."
        )
    else:
        _notes.append(
            "Native FORM is unavailable here. The table shows the supplied reference; "
            "set `FORM_EXECUTABLE` in a native notebook session to measure it locally."
        )
    if ladder_stages:
        _seconds = ladder_stages[-1]["Spenso cumulative wall (s)"]
        _notes.append(
            f"**Full Spenso reduction:** {_seconds:.3f} s, with all eight term counts checked."
        )
        if ladder_form_process_samples:
            _notes.append(
                f"Spenso in-process wall / FORM process wall = "
                f"{_seconds / median(ladder_form_process_samples):.2f}×."
            )
    else:
        _notes.append(
            "Run the full reduction above to fill in the Spenso stage timings."
        )
    if ladder_exact_match:
        _notes.append(
            "**Exact check passed:** all coefficients of the 9,652-term scalar polynomial agree with FORM."
        )
    mo.vstack(
        [
            *(mo.md(note) for note in _notes),
            mo.ui.table(_rows, selection=None, pagination=False, show_download=True),
            mo.accordion(
                {
                    "Reference FORM source": mo.md(
                        f"```form\n{ladder_form_source}\n```"
                    ),
                    "Local FORM output": mo.md(
                        f"```text\n{ladder_form_stdout or 'Not run'}\n```"
                    ),
                }
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(
    boundary_checks,
    benchmark_results,
    hep_component_checks,
    metric_pairings,
    form_records,
    chiral_checks,
    dimension_checks,
    epsilon_result,
    expanded_trace_check,
    mo,
    open_checks,
    ordering_result,
    slash_check,
    trace_checks,
    trace_terms,
    ladder_substitution_seconds,
    ladder_form_stages,
    ladder_exact_match,
):
    mo.Html(
        '<p data-notebook-ready="gamma_simplification">All identity, boundary and HEP component checks passed.</p>'
    )
    return


if __name__ == "__main__":
    app.run()
