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
def _(S, TensorExpression, a, b, identity, identity_card, lorentz, mo, spin, word):
    p, q = S("gamma_tutorial::p", "gamma_tutorial::q")
    gamma_head, metric_head = S("spenso::gamma", "spenso::g")
    p_compact, q_compact = p(lorentz.to_expression()), q(lorentz.to_expression())
    slash_p = TensorExpression(
        gamma_head(spin(a).to_expression(), spin(b).to_expression(), p_compact)
    )
    slash_q = TensorExpression(
        gamma_head(spin(a).to_expression(), spin(b).to_expression(), q_compact)
    )
    # These factors carry explicit endpoints, so make the chain-local versions
    # through the registered in/out form rather than inventing sewing labels.
    incoming, outgoing, chain_head = S("spenso::in", "spenso::out", "spenso::chain")
    slash_sandwich = TensorExpression(
        chain_head(
            spin(a).to_expression(),
            spin(b).to_expression(),
            gamma_head(incoming, outgoing, p_compact),
            gamma_head(incoming, outgoing, q_compact),
            gamma_head(incoming, outgoing, p_compact),
        )
    )
    p_squared = metric_head(p_compact, p_compact)
    p_dot_q = metric_head(p_compact, q_compact)
    slash_check = identity_card(
        "Slash sandwich: 2(p·q) slash(p) − p² slash(q)",
        slash_sandwich,
        2 * p_dot_q * word(slash_p) - p_squared * word(slash_q),
    )
    slash_check
    return chain_head, gamma_head, incoming, outgoing, slash_check


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 4. Ordinary traces

    Odd ordinary traces vanish. Even traces recurse to pairwise metrics;
    four gammas give the familiar three signed pairings. Disable trace
    evaluation when retaining closed spin chains is preferable downstream.
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
def _(E, S, TensorExpression, TensorName, slash_length, slash_pattern):
    from functools import reduce
    from operator import mul

    _bis, _mink, _gamma, _trace, _cyclic, _in, _out, _g = S(
        *(
            f"spenso::{name}"
            for name in ("bis", "mink", "gamma", "trace", "cyclic", "in", "out", "g")
        )
    )
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
    _indices = [_mink(4, S(f"gamma_benchmark::mu{i}")) for i in range(_n)]
    bare_trace = TensorExpression(
        _trace(_bis(4), _cyclic(*(_gamma(_in, _out, i) for i in _indices)))
    )
    momentum_factors = reduce(mul, (p(i) for p, i in zip(_momenta, _indices)))
    indexed_trace = TensorExpression(momentum_factors * bare_trace.to_expression())
    compact_trace = TensorExpression(
        _trace(_bis(4), _cyclic(*(_gamma(_in, _out, p(_mink(4))) for p in _momenta)))
    )

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
    scalar_expected = (
        scalar_oracle.replace(_pp, _g(_p(_mink(4)), _p(_mink(4))))
        .replace(_qq, _g(_q(_mink(4)), _q(_mink(4))))
        .replace(_pq, _g(_p(_mink(4)), _q(_mink(4))))
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

    Its CSV distinguishes first-call Idenso time, warm Idenso time, FORM
    `trace4` process time and FORM `tracen` process time. First calls in a
    fresh process expose the lazy table cost; repeated notebook runs may
    already have populated those tables.
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
    """)
    return


@app.cell
def _(E, S, TensorExpression, TensorName, form_executable, lorentz):
    from itertools import permutations

    from symbolica import Expression
    from symbolica.community.spenso import (
        Tensor,
        TensorLibrary,
        TensorNetwork,
    )

    _bis, _mink, _gamma, _g5, _trace, _cyclic, _inside, _outside, _metric = S(
        *(
            f"spenso::{name}"
            for name in (
                "bis",
                "mink",
                "gamma",
                "gamma5",
                "trace",
                "cyclic",
                "in",
                "out",
                "g",
            )
        )
    )
    _names = [TensorName.vector(f"gamma_hep::p{i}") for i in range(10)]
    _momenta = [name.to_expression()(_mink(4)) for name in _names]
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
            _factors = [_gamma(_inside, _outside, p) for p in _momenta[:_length]]
            if _axial:
                _factors.insert(0, _g5(_inside, _outside))
            _original = _trace(_bis(4), _cyclic(*_factors))
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
def _(
    boundary_checks,
    benchmark_results,
    hep_component_checks,
    metric_pairings,
    form_records,
    chiral_checks,
    dimension_checks,
    epsilon_result,
    mo,
    open_checks,
    ordering_result,
    slash_check,
    trace_checks,
    trace_terms,
):
    mo.Html(
        '<p data-notebook-ready="gamma_simplification">All identity, boundary and HEP component checks passed.</p>'
    )
    return


if __name__ == "__main__":
    app.run()
