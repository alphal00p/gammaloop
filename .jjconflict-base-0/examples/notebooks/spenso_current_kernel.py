# Run with a combined Symbolica host exposing symbolica.community.spenso:
#   python -m marimo edit examples/notebooks/spenso_current_kernel.py
# Additional packages: marimo==0.24.0, typst==0.15.0.
# The host build lives in examples/notebooks/symbolica-host; a stock Symbolica
# wheel need not expose this checkout's typed Spenso API.

import marimo

__generated_with = "0.24.0"
app = marimo.App(
    width="medium",
    app_title="Spenso: from a vertex to a current",
)


@app.cell
def _():
    import marimo as mo
    from symbolica import Expression
    from symbolica.community.spenso import AUTO, trace

    return AUTO, Expression, mo, trace


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # From a vertex to a current with Spenso

    Derive the vector current of a fermion–antifermion–vector vertex, with
    the coupling omitted:

    \[
      K^\mu = \sum_{a,b=0}^3 \bar J_a\,(\gamma^\mu)_{ab}\,J_b.
    \]

    This is the example in the **Deriving current kernels with Spenso**
    section of the [pyAmpliCol paper](https://github.com/ValentinHirschi/pyAmpliColPaper/blob/074fa4f78de9a52ef773ee11e452675502474a62/sections/implementation-design.tex#L187).
    We supply all eight symbolic input components, contract the spinor
    indices, and keep the Lorentz index open. Here `Jbar` is an already
    barred current supplied independently of `J`; no conjugation is implicit.

    The current API generates unknown components automatically and provides
    simplification methods on
    `TensorExpression`. Gamma's slots are **(row, column, Lorentz)**:
    `gamma(a, b, mu)` implements the equation above. Reversing `a, b`, as
    in the paper's listing, transposes the matrix in this contraction.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Full paper section — Markdown working draft

    The **Original LaTeX** tab preserves the complete Spenso subsection from
    [paper revision `074fa4f`](https://github.com/ValentinHirschi/pyAmpliColPaper/blob/074fa4f78de9a52ef773ee11e452675502474a62/sections/implementation-design.tex#L187),
    including its original listing and model-support discussion.

    **Read draft** shows the Markdown working revision, a first pass using the executable listing below.
    It clarifies automatic symbolic components, corrects gamma-slot ordering
    and component extraction, and describes Idenso without adding a redundant
    simplification to the single-gamma example. The original model-support
    claims and cross-references are retained for our next prose pass.

    Use **Edit Markdown** to revise the text; **Read draft** previews the result.
    Download the Markdown to keep editor changes. Equations use math notation,
    citation keys are retained as `[@key]`, and references to other paper
    sections use descriptive names. We can convert back to LaTeX later.
    """)
    return


@app.cell(hide_code=True)
def _():
    from inspect import cleandoc as _cleandoc

    paper_section_original = (
        _cleandoc(r"""
    \subsection{Deriving current kernels with \spenso{}}
    \label{sec:spenso}
    \lucien{Improve prose}
    UFO Lorentz rules are indexed expressions rather than lists of numerical
    component formulae.  \pyac{} chooses the current joins using the rules of
    section~\ref{sec:general-contact-decomposition}; \spenso{}~\cite{SpensoSoftware}
    supplies the index spaces and tensor contractions which derive their component
    expressions and hence their numerical kernels.
    Automatic derivation of these component rules is also used in
    \aloha{}~\cite{deAquino:2011ub} and the UFO extension of
    \comix{}~\cite{Hoeche:2014kca}.  In \pyac's checked four-point
    decompositions, the uncontracted tensor indices determine the components of
    the auxiliary current.

    The translation is explicit.  Lorentz indices belong to a four-dimensional
    Minkowski-vector representation and fermion indices to a four-component
    bispinor representation; fundamental and adjoint $SU(3)$ representations are
    used for colour.  \spenso{}'s high-energy-physics tensor library provides the
    metric, Dirac matrices, $\gamma^5$ and chiral projectors.  Momentum vectors and
    the components of the currents entering a join are supplied as component
    tensors.  A tensor network then contracts every repeated index while leaving
    the index of the chosen result field open.  The \idenso{} simplification routines
    distributed with \symbolica{}, \texttt{simplify\_metrics},
    \texttt{simplify\_gamma} and \texttt{simplify\_color}, reduce the corresponding
    metric, Dirac and colour identities~\cite{SymbolicaSoftware}.

    For example, the Standard Model UFO expression \texttt{Gamma(3,1,2)} denotes
    $(\gamma^{\mu_3})_{s_1s_2}$.  With an antifermion current $\bar J_1$ and a
    fermion current $J_2$, choosing the vector as the result field and omitting
    the coupling gives
    \begin{equation}
      \bigl[K_3(J_1,J_2)\bigr]^{\mu_3}
      =\sum_{s_1,s_2}(\bar J_1)_{s_1}
       (\gamma^{\mu_3})_{s_1s_2}(J_2)_{s_2} .
      \label{eq:spenso-gamma-example}
    \end{equation}
    \Needspace{25\baselineskip}
    The essential construction is illustrated below; the arrays holding the two
    input currents and uninteresting naming details are elided.
    \spenso{}'s gamma primitive takes its spinor indices in input--output order,
    which explains their order in the code.

    \begin{lstlisting}[language=Python,caption={A vector current constructed with \spenso{} and simplified with \idenso{}.},label={lst:spenso-current}]
    from symbolica.community.idenso import simplify_gamma, simplify_metrics
    from symbolica.community.spenso import (
        LibraryTensor, Representation, TensorLibrary, TensorName, TensorNetwork,
    )
    lib = TensorLibrary.hep_lib_atom()
    S = Representation.bis(4)
    s1, s2 = S("s1"), S("s2")
    Jbar, J = TensorName("Jbar"), TensorName("J")
    lib.register(LibraryTensor.dense(Jbar(S), Jbar_components))
    lib.register(LibraryTensor.dense(J(S), J_components))
    gamma = lib[TensorName.gamma()]
    expr = (Jbar(s1).to_expression() *
            gamma("s2", "s1", "mu") * J(s2).to_expression())
    expr = simplify_metrics(simplify_gamma(expr))
    net = TensorNetwork(expr, lib)
    net.execute(library=lib)
    K = net.result_tensor(lib)
    K.to_dense()  # four components
    \end{lstlisting}

    The call to \texttt{result\_tensor} retains the unmatched index $\mu_3$;
    \texttt{to\_dense} returns four ordered component expressions, one for each of
    its values.  For a momentum-independent four-point rule handled by the checked
    split, two indices may instead remain open.  This produces the 16
    components of $X_{\rho\sigma}$ in the examples of
    section~\ref{sec:general-contact-decomposition}.  \spenso{} performs the indexed
    contraction; \pyac{} then uses \symbolica{} to identify exact zeros and equality
    or sign relations.  For the pairings shown there, this reduces the stored
    components to one representative for $g_{\mu\nu}g_{\rho\sigma}$ and to six
    for the antisymmetric combination in
    eq.~\eqref{eq:antisymmetric-contact}.  More generally, \pyac{} chooses how a
    rule is split into currents, \spenso{} carries out its indexed contractions, and
    \symbolica{} simplifies and compares the resulting scalar expressions.  The
    same machinery also handles the supported higher-point colour-singlet
    rules: when no reduced tensor factorisation is used, \spenso{} contracts the
    original interaction with its input currents to derive the final join.
    The component expressions are then turned into numerical kernels for both
    recurrence and compiled execution.  Thus the tensor algebra serves both the
    analysis of a supported decomposition and its numerical implementation.
    These tools are used during model processing and kernel construction, not
    ordinary double-precision event evaluation.

    This division also makes the limits of model support explicit.  Unresolved
    UFO tensor functions are rejected.  For momentum-independent four-point
    rules, the open-component ordering, zeros and exact sign relations are
    recorded for every inequivalent result field.  Supported coloured contact
    rules additionally have to match one of the explicitly reconstructed colour
    classes described in section~\ref{sec:general-contact-decomposition}.  New
    Lorentz structures can be added by defining their components and contraction
    rules.  New colour representations, such as sextets and their Clebsch--Gordan
    tensors, would also require the
    associated colour flows, projections and contractions.  Other restrictions, including
    Majorana fermion flow and unverified coloured high-point contacts, are listed
    in table~\ref{tab:ufo-coverage}.

    Within these model boundaries, the same current construction supports
    multiple open quark lines and LC, NLC and full-colour calculations.
    """)
        + "\n"
    )
    paper_section_draft = (
        _cleandoc(r"""
    ## Deriving current kernels with Spenso
    To turn a UFO Lorentz rule into a current kernel, pyAmpliCol selects an
    output leg and supplies currents for the remaining legs. Spenso contracts
    their indices, leaving the output leg’s indices open.

    For each component of the output current, Spenso produces a scalar symbolic
    expression in the components of the input currents and any momenta or model
    parameters. For a vector current, these are the expressions for
    $K^0$, $K^1$, $K^2$ and $K^3$ in a chosen basis. The contracted indices
    have been summed over; each expression corresponds to a fixed value of the
    remaining index. Symbolica simplifies these expressions before they are
    turned into numerical kernels.

    This construction is similar in purpose to ALOHA [@deAquino:2011ub]
    and the UFO extension of Comix [@Hoeche:2014kca]. Spenso [@SpensoSoftware] provides a
    generic tensor-contraction framework in which index spaces and tensor
    component data are specified independently of the contraction algorithm.
    Different explicit tensor representations, such as a different basis for
    the Dirac matrices, can therefore be used by supplying the corresponding
    tensor and current components, without rewriting the contraction code.

    Spenso assigns a representation to each tensor index: a four-dimensional
    Minkowski-vector space for Lorentz indices, a four-component bispinor space
    for fermion indices, and fundamental or adjoint $SU(3)$ spaces for colour.
    Its high-energy-physics tensor library supplies the metric, Dirac matrices,
    $\gamma^5$ and chiral projectors. Named tensors absent from the library
    acquire symbolic components when the network is parsed; explicit component
    data can instead be registered when required. Repeated compatible indices
    are contracted, while unmatched indices form the result's tensor interface.
    The Idenso operations exposed through `TensorExpression`, including
    `contract`, `simplify_gamma` and
    `simplify_color`, can reduce metric, Dirac and colour identities
    before component evaluation [@SymbolicaSoftware].

    For example, the Standard Model UFO expression `Gamma(3,1,2)` denotes
    $(\gamma^{\mu_3})_{s_1s_2}$.  With an antifermion current $\bar J_1$ and a
    fermion current $J_2$, choosing the vector as the output leg and omitting
    the coupling gives

    $$
      \bigl[K_3(J_1,J_2)\bigr]^{\mu_3}
      =\sum_{s_1,s_2}(\bar J_1)_{s_1}
       (\gamma^{\mu_3})_{s_1s_2}(J_2)_{s_2} .
    $$

    The following complete example derives this current with symbolic inputs.
    Here $\bar J_1$ is supplied as an independent barred current; the code does
    not construct it by conjugating $J_2$. The built-in HEP library uses Weyl
    matrices and the metric $g=\mathrm{diag}(1,-1,-1,-1)$.
    Spenso's gamma slots have the order `(row, column, Lorentz)`, so the
    first spinor slot is contracted with $\bar J_1$ and the second with $J_2$.
    A single gamma matrix requires no Dirac-algebra simplification.

    ```python
    from symbolica.community.spenso import (
        Representation,
        TensorExpression,
        TensorLibrary,
        TensorName,
    )

    spinor = Representation.bis(4)
    Jbar = TensorName("Jbar")(spinor)
    J = TensorName("J")(spinor)
    gamma = TensorExpression.gamma(4)

    # Gamma slots: row, column, Lorentz. Only mu remains open.
    current = Jbar("a") * gamma("a", "b", "mu") * J("b")
    library = TensorLibrary.hep_lib_atom()
    network = current.to_network(library=library)
    network.execute(library=library)
    kernel = network.result_tensor(library=library)
    bar_components, components = Jbar.components(), J.components()
    print(kernel)  # Components in the order mu = 0, 1, 2, 3.
    ```

    The result `kernel` retains the unmatched Lorentz index $\mu_3$.
    Its slice `kernel[:]` returns the four component expressions in
    logical order, $\mu_3=0,1,2,3$. The calls to `components()`
    retrieve each input current's generated scalar expressions, with the same
    identities used in the contracted result. These expressions can be passed
    directly as numerical-evaluator parameters; no separate array of input
    symbols is needed. By contrast, `to_dense()` changes a tensor's
    storage in place and does not return its components.

    For a momentum-independent four-point rule handled by the checked
    split, two indices may instead remain open.  This produces the 16
    components of $X_{\rho\sigma}$ in the examples of
    the section on general contact decomposition.  Spenso performs the indexed
    contraction; pyAmpliCol then uses Symbolica to identify exact zeros and equality
    or sign relations.  For the pairings shown there, this reduces the stored
    components to one representative for $g_{\mu\nu}g_{\rho\sigma}$ and to six
    for the antisymmetric combination in
    the antisymmetric-contact equation.  More generally, pyAmpliCol chooses how a
    rule is split into currents, Spenso carries out its indexed contractions, and
    Symbolica simplifies and compares the resulting scalar expressions.  The
    same machinery also handles the supported higher-point colour-singlet
    rules: when no reduced tensor factorisation is used, Spenso contracts the
    original interaction with its input currents to derive the final join.
    The component expressions are then turned into numerical kernels for both
    recurrence and compiled execution.  Thus the tensor algebra serves both the
    analysis of a supported decomposition and its numerical implementation.
    These tools are used during model processing and kernel construction, not
    ordinary double-precision event evaluation.

    This division also makes the limits of model support explicit.  Unresolved
    UFO tensor functions are rejected.  For momentum-independent four-point
    rules, the open-component ordering, zeros and exact sign relations are
    recorded for every inequivalent output leg.  Supported coloured contact
    rules additionally have to match one of the explicitly reconstructed colour
    classes described in the section on general contact decomposition.  New
    Lorentz structures can be added by defining their components and contraction
    rules.  New colour representations, such as sextets and their Clebsch--Gordan
    tensors, would also require the
    associated colour flows, projections and contractions.  Other restrictions, including
    Majorana fermion flow and unverified coloured high-point contacts, are listed
    in the UFO-coverage table.

    Within these model boundaries, the same current construction supports
    multiple open quark lines and LC, NLC and full-colour calculations.
    """)
        + "\n"
    )
    return paper_section_draft, paper_section_original


@app.cell(hide_code=True)
def _(mo, paper_section_draft):
    section_editor = mo.ui.code_editor(
        paper_section_draft,
        language="markdown",
        min_height=500,
        max_height=900,
        label="Working draft (download to keep editor changes)",
    )
    return (section_editor,)


@app.cell(hide_code=True)
def _(mo, paper_section_original, section_editor):
    mo.vstack(
        [
            mo.ui.tabs(
                {
                    "Read draft": mo.md(section_editor.value),
                    "Edit Markdown": section_editor,
                    "Original LaTeX": mo.ui.code_editor(
                        paper_section_original,
                        language="latex",
                        disabled=True,
                        min_height=500,
                        max_height=900,
                        label="Full original subsection",
                    ),
                }
            ),
            mo.download(
                section_editor.value.encode("utf-8"),
                filename="spenso-section.md",
                mimetype="text/markdown",
                label="Download Markdown draft",
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Suggested paper listing — copy the next cell

    The next cell is self-contained Python. Copy its body into the paper;
    no other notebook cell or Marimo helper is needed. Unknown currents acquire
    symbolic components automatically. The exact HEP library supplies the Weyl
    matrices; `components()` exposes the same generated input symbols afterward.

    **Suggested caption:** Deriving the four symbolic components of
    $K^\mu=\bar J\gamma^\mu J$ with Spenso. Repeated spinor indices are
    contracted, while the Lorentz index remains open.
    """)
    return


@app.cell
def _():
    from symbolica.community.spenso import (
        Representation,
        TensorExpression,
        TensorLibrary,
        TensorName,
    )
    from symbolica import PrintMode

    spinor = Representation.bis(4)
    Jbar = TensorName("Jbarrr", print={"typst": "macron(J)"})(spinor)
    J = TensorName("J")(spinor)
    gamma = TensorExpression.gamma(4)

    # Gamma slots: row, column, Lorentz. Only mu remains open.
    current = Jbar("a") * gamma("a", "b", "mu") * J("b")
    library = TensorLibrary.hep_lib_atom()
    network = current.to_network(library=library)
    network.execute(library=library)
    kernel = network.result_tensor(library=library)
    kernel  # Components in the order mu = 0, 1, 2, 3.
    return Representation, TensorExpression, kernel, network


@app.cell(hide_code=True)
def _(kernel, mo, network):
    mo.vstack(
        [
            mo.md("## The indexed rule and its four component expressions"),
            mo.Html(network.expression().to_html()),
            mo.Html(kernel.to_html()),
            mo.md(
                "`kernel[:]` reads the components in logical order, μ = 0, 1, 2, 3. "
                "`to_dense()` only changes storage in place; it does not return a list."
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Check every component exactly

    The built-in HEP library uses Weyl matrices with metric
    $g=\mathrm{diag}(1,-1,-1,-1)$:
    $\gamma^0=\left(\begin{smallmatrix}0&I\\I&0\end{smallmatrix}\right)$ and
    $\gamma^i=\left(\begin{smallmatrix}0&\sigma^i\\-\sigma^i&0\end{smallmatrix}\right)$.
    An independent matrix calculation checks all four symbolic polynomials,
    including their signs and the imaginary entries of $\gamma^2$.
    """)
    return


@app.cell
def _(Expression, Representation, bar_components, components, kernel):
    assert kernel.structure.slots == (Representation.mink(4)("mu"),)
    assert len(kernel) == 4
    _i = Expression.I
    _weyl = (
        ((0, 0, 1, 0), (0, 0, 0, 1), (1, 0, 0, 0), (0, 1, 0, 0)),
        ((0, 0, 0, 1), (0, 0, 1, 0), (0, -1, 0, 0), (-1, 0, 0, 0)),
        ((0, 0, 0, -_i), (0, 0, _i, 0), (0, _i, 0, 0), (-_i, 0, 0, 0)),
        ((0, 0, 1, 0), (0, 0, 0, -1), (-1, 0, 0, 0), (0, 1, 0, 0)),
    )
    for _mu, _matrix in enumerate(_weyl):
        _expected = sum(
            bar_components[_a] * _matrix[_a][_b] * components[_b]
            for _a in range(4)
            for _b in range(4)
        )
        assert (kernel[_mu] - _expected).expand() == 0
    symbolic_check = "All four components agree exactly with the Weyl matrices."
    return (symbolic_check,)


@app.cell(hide_code=True)
def _(mo):
    _intro = mo.md("""
    ## Reuse the numerical kernel

    Build the evaluator once, then vary the input without repeating the tensor
    contraction. For `Jbar = [1, 2, 3, 4]` and `J = [5, 6, 7, 8]`, the result
    is `[62, -16, 4j, 0]`. The slider rescales `J`.
    """)
    scale = mo.ui.slider(-2, 2, step=0.25, value=1, label="Scale of J", show_value=True)
    mo.vstack([_intro, scale])
    return (scale,)


@app.cell
def _(bar_components, components, kernel):
    evaluator = kernel.evaluator(
        constants={},
        funs={},
        params=[*bar_components, *components],
        n_cores=1,
    )
    return (evaluator,)


@app.cell
def _(evaluator, mo, scale, symbolic_check):
    _inputs = [1, 2, 3, 4, *[scale.value * _x for _x in (5, 6, 7, 8)]]
    numeric_current = evaluator.evaluate_complex([_inputs])[0]
    assert numeric_current[:] == [scale.value * _x for _x in (62, -16, 4j, 0)]
    mo.vstack(
        [
            mo.Html(numeric_current.to_html()),
            mo.callout(
                mo.md(symbolic_check + " Numerical evaluation also checked."),
                kind="success",
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## An Idenso identity, before component evaluation

    A single gamma matrix needs no Dirac simplification. To show a useful
    symbolic rewrite, close a two-gamma chain instead:
    $\mathrm{tr}(\gamma^\mu\gamma^\nu)=4g^{\mu\nu}$.
    `AUTO` supplies the spinor endpoints that `trace` closes.
    """)
    return


@app.cell
def _(AUTO, Representation, TensorExpression, mo, trace):
    lorentz = Representation.mink(4)
    _gamma = TensorExpression.gamma(4)
    gamma_trace = trace(
        Representation.bis(4),
        _gamma(AUTO, AUTO, lorentz("mu")),
        _gamma(AUTO, AUTO, lorentz("nu")),
    )
    reduced_trace = gamma_trace.simplify_gamma().to_expression()
    assert (
        reduced_trace.to_expression() - 4 * lorentz.g("mu", "nu").to_expression()
    ).expand() == 0
    mo.vstack([mo.Html(gamma_trace.to_html()), mo.Html(reduced_trace.to_html())])
    return


if __name__ == "__main__":
    app.run()
