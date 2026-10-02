# ruff: noqa: B018, PLR1711 -- marimo uses final expressions and explicit cell returns.

import marimo

__generated_with = "0.24.0"
app = marimo.App(
    width="medium",
    app_title="Deriving current kernels with Spenso",
)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Deriving current kernels with Spenso
    To turn a UFO Lorentz rule into a current kernel, pyAmpliCol selects an
    output leg and supplies currents for the remaining legs. Spenso
    contracts their indices, leaving the output leg's indices open.

    For each component of the output current, Spenso produces a scalar
    symbolic expression in the components of the input currents and any
    momenta or model parameters. For the spinor current $Q(2,4)_c$ constructed
    below, these are four expressions, one for each value $c=0,1,2,3$ in a
    chosen basis. The contracted indices have been summed over, while $c$
    labels the remaining spinor index. Symbolica can then optimise these
    expressions before they are turned into numerical kernels.

    This construction is similar in purpose to
    ALOHA [@deAquino:2011ub] and the UFO extension of
    Comix [@Hoeche:2014kca]. Spenso [@SpensoSoftware]
    provides a generic tensor-contraction framework in which index spaces
    and tensor component data are specified independently of the
    contraction algorithm. Different explicit tensor representations,
    such as a different basis for the Dirac matrices, can therefore be
    used by supplying the corresponding tensor and current components,
    without rewriting the contraction code.

    Spenso assigns a representation to each tensor index: a
    four-dimensional Minkowski-vector space for Lorentz indices, a
    four-component bispinor space for fermion indices, and fundamental
    or adjoint $SU(3)$ spaces for colour. Its high-energy-physics tensor
    library supplies the metric, Dirac matrices, $\gamma^5$ and chiral
    projectors. For tensors absent from the library, Spenso generates symbolic
    components automatically. Users can instead supply their own symbolic or
    numerical components through the tensor library. Repeated compatible
    indices are contracted, while unmatched indices remain open. The Idenso
    operations exposed through `TensorExpression`, including `simplify_algebra`,
    can apply metric, Dirac and colour identities before component evaluation.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    To illustrate the quark–gluon join
    $Q(2,4)=P_d V_{dgd}[Q(2),G(4)]$ in equation (3.6), in the section on
    off-shell-current recurrence, we define its tensors directly. For this
    standalone example we use a four-component row spinor and choose the vertex
    normalization $V^\mu=i\gamma^\mu$. Let $q_{24}=q_2+q_4$ and
    $s_{24}=q_{24}^2$. With a massless quark propagator,

    $$
    Q(2,4)_c = \frac{i}{s_{24}}\sum_b
    \left[i\sum_{a,\mu}Q(2)_a(\gamma^\mu)_{ab}G(4)_\mu\right]
    (\not q_{24})_{bc}.
    $$

    The two factors of $i$ come from the chosen vertex and propagator. The colour
    matrix and strong coupling are omitted.
    The propagator acts on the right because $Q$ is a row spinor. Components of
    vectors are supplied contravariantly; contraction in the Minkowski space
    supplies the signs of $g=\mathrm{diag}(1,-1,-1,-1)$.

    The following listings assemble this expression using the tensor API and
    the library's Weyl gamma matrices. We build the amputated join, attach the
    propagator, and evaluate the components. For a current with several vertex
    contributions, these would be summed before applying the propagator once.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    We define the input currents $Q(2)$ and $G(4)$, and the momentum $q_{24}$.
    `TN.vector` declares a one-rank tensor; the representation specifies
    whether that axis belongs to the bispinor or Lorentz-vector space.
    The integer arguments before the representation identify the associated
    external legs: `2` for $Q(2)$, `4` for $G(4)$, and `2, 4` for $q_{24}$.
    Tensor indices are assigned when these objects enter a contraction.
    For example, `Q2(1)` assigns the spinor index label `1` to $Q(2)$.
    """)
    return


@app.cell
def _():
    from symbolica import Symbol
    from symbolica.community.tensor import Representation as Rep, TensorExpression as T, TensorName as TN

    spinor, vector = Rep.bis(4), Rep.mink(4)
    Q2 = TN.vector("Q")(2,spinor)
    G4 = TN.vector("G")(4,vector)
    q24 = TN.vector("q")(2,4,vector)

    Q2, G4, q24
    return G4, Q2, Symbol, T, q24


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Contract the input currents with the vertex $i\gamma^\mu$ to form the
    amputated join. The gamma slots are ordered `(row, column, Lorentz)`.
    The label `1` pairs the bispinor slot of `Q2` with the first bispinor slot of the gamma tensor, and label `2` does the same for the Lorentz slot of `G4`. `AUTO`, written `_`, leaves
    the remaining spinor port available for the next multiplication.
    """)
    return


@app.cell
def _(G4, Q2, Symbol, T):
    from symbolica.community.tensor import AUTO as _
    gamma = T.dirac_gamma(4)
    i = Symbol.I
    join = i * gamma(1, _, 2) * Q2(1) * G4(2)
    join
    return gamma, i, join


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    The propagator is $i\not q_{24}/q_{24}^2$. Multiplying the unlabelled vector
    `q24` by itself forms its Minkowski scalar product. The first spinor port
    of `gamma(_, 3, 3)` will connect to the join; its second spinor port is
    also indexed `3` and becomes the output index. Note that since the bispinor and lorentz occupy different spaces, we can reuse the same index (3) without fear of overlap.
    """)
    return


@app.cell
def _(gamma, i, q24):
    from symbolica.community.tensor import AUTO as _

    s24 = q24 * q24
    propagator = i * gamma(_, 3, 3) * q24(3) / s24
    propagator
    return (propagator,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Attach the propagator on the right. Multiplication connects the two
    unlabelled spinor ports, leaving the labelled output port open.
    """)
    return


@app.cell
def _(join, propagator):
    current = join * propagator
    current
    return (current,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Evaluate the contractions to obtain the four component expressions.
    The result `Q24` retains one open bispinor index and includes the
    propagator denominator.
    """)
    return


@app.cell
def _(current):
    Q24 = current.to_tensor()
    Q24
    return (Q24,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Build a numerical evaluator from these expressions. `components()` provides
    the input symbols in component order, so they can be used directly as the
    evaluator's parameters. The denominator $s_{24}$ is computed from the
    momentum components and needs no separate input. Once constructed, the
    evaluator can be reused for different currents and momenta.
    """)
    return


@app.cell
def _(G4, Q2, Q24, q24):
    parameters = [*Q2.components(), *G4.components(), *q24.components()]
    evaluator = Q24.evaluator(constants={}, funs={}, params=parameters, n_cores=1)

    value = evaluator.evaluate_complex([[1, 0, 0, 0, 0, 0, 1, 0, 2, 1, 0, 1]])[0]
    value
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    The numerical example uses $q_2=(1,0,0,1)$, $q_4=(1,1,0,0)$,
    $Q(2)=(1,0,0,0)$ and
    $G(4)=(0,0,1,0)$. These inputs satisfy $q_2^2=q_4^2=0$,
    $Q(2)\not q_2=0$ and $q_4\cdot G(4)=0$; the spinor's overall normalisation
    is arbitrary. Here $q_{24}=(2,1,0,1)$ and $s_{24}=2$, giving
    $Q(2,4)=(-i/2,-i/2,0,0)$. This is a local current evaluation, not a complete
    five-leg phase-space point.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    When pyAmpliCol decomposes a momentum-independent four-point interaction
    into successive current joins, as described in
    the section on general contact decomposition, an intermediate
    current can carry two open Lorentz indices.  Spenso performs
    the indexed contraction; this gives the 16
    components of $X_{\rho\sigma}$, the intermediate current defined there by contracting the four-point tensor with its first two input currents, $\rho$ and $\sigma$ being its open Lorentz indices. pyAmpliCol then uses Symbolica to identify
    exact zeros and equality or sign relations. For the pairings shown
    there, this reduces the stored components to one representative for
    $g_{\mu\nu}g_{\rho\sigma}$ and to six for the antisymmetric combination
    in the antisymmetric-contact equation. More generally, pyAmpliCol
    chooses how a rule is split into currents, Spenso carries out its
    indexed contractions, and Symbolica simplifies and compares the
    resulting scalar expressions. The same machinery also handles the
    supported higher-point colour-singlet rules: when no reduced tensor
    factorisation is used, Spenso contracts the original interaction
    with its input currents to derive the final join. The component
    expressions are then turned into numerical kernels for both
    recurrence and compiled execution. Thus the tensor algebra serves
    both the analysis of a supported decomposition and its numerical
    implementation. These tools are used during model processing and
    kernel construction, not during ordinary double-precision event evaluation.

    This division also makes the limits of model support explicit.
    Unresolved UFO tensor functions are rejected. For
    momentum-independent four-point rules, the open-component ordering,
    zeros and exact sign relations are recorded for every inequivalent
    output leg. Supported coloured contact rules additionally have to
    match one of the explicitly reconstructed colour classes described
    in the section on general contact decomposition. New Lorentz
    structures can be added by defining their components and
    contraction rules. New colour representations, such as sextets and
    their Clebsch–Gordan tensors, would also require the associated
    colour flows, projections and contractions. Other restrictions,
    including Majorana fermion flow and unverified coloured high-point
    contacts, are listed in the UFO-coverage table.

    Within these model boundaries, the same current construction supports
    multiple open quark lines and LC, NLC and full-colour calculations.
    """)
    return


@app.cell(hide_code=True)
def _():
    import marimo as mo

    return (mo,)


if __name__ == "__main__":
    app.run()
