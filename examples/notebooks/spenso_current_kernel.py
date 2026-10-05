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

    The numerical kernels described in the preceding subsections require
    scalar expressions for the components of the off-shell currents. In
    the recurrence introduced above, these expressions are implicit in
    the vertex rules and propagators.
    Spenso [@SpensoSoftware] provides the tensor contractions needed to
    obtain them. During model processing, pyAmpliCol selects an output leg
    of a UFO Lorentz rule and supplies symbolic currents for the remaining
    legs. Spenso sums over their indices in a chosen basis, leaving one
    scalar expression for each component of the output current. These are
    the expressions passed to Symbolica for optimisation and compilation.

    To illustrate this step, return to the $d\bar d\to Zgg$ example and
    consider the last current in equation (3.6),

    $$Q(2,4)=P_d V_{dgd}[Q(2),G(4)].$$

    Here $Q(2)=\bar u_{h_2}(q_2)$ and $G(4)$ are the one-leg currents
    introduced there. We use Spenso's tensor API to derive the components
    of their quark–gluon join and then apply the propagator $P_d$, keeping
    the input components symbolic until numerical evaluation.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    For this example, take a massless quark and choose the colour-stripped
    vertex normalisation $V^\mu=i\gamma^\mu$ in four-component notation,
    with the strong coupling omitted as in equation (3.6). We use the Weyl gamma
    matrices supplied by Spenso's high-energy-physics tensor library and
    the metric $g=\mathrm{diag}(1,-1,-1,-1)$. With
    $q_{24}=q_2+q_4$ and $s_{24}=q_{24}^2$, the current is then

    $$
    Q(2,4)_c = \frac{i}{s_{24}}\sum_b
    \left[i\sum_{a,\mu}Q(2)_a(\gamma^\mu)_{ab}G(4)_\mu\right]
    (\not q_{24})_{bc}.
    $$

    Here $a,b,c$ are spinor indices and $\mu$ is a Lorentz index, each
    taking four values, with $\not q_{24}=\gamma^\nu q_{24,\nu}$.
    The bracket is $V_{dgd}[Q(2),G(4)]$; since $Q(2)$ is a row spinor,
    the propagator $i\not q_{24}/s_{24}$ acts on its right. Summing the
    repeated indices therefore leaves four expressions labelled by $c$.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    To express this contraction in code, we first specify the spaces to
    which the indices belong: `Rep.bis(4)` for the four spinor components
    and `Rep.mink(4)` for Minkowski vectors. The latter includes the metric
    used in contractions, so vector components can be supplied
    contravariantly. With these representations, `TN.vector` declares
    a tensor carrying one index. The arguments `2`, `4` and `2, 4` retain
    the leg labels of $Q(2)$, $G(4)$ and $q_{24}$ in their names; tensor
    indices are assigned when the contractions are formed below.
    """)
    return


@app.cell
def _():
    from symbolica import Symbol
    from symbolica.community.tensor import Representation as Rep
    from symbolica.community.tensor import TensorExpression as T
    from symbolica.community.tensor import TensorName as TN

    spinor, vector = Rep.bis(4), Rep.mink(4)
    Q2 = TN.vector("Q")(2, spinor)
    G4 = TN.vector("G")(4, vector)
    q24 = TN.vector("q")(2, 4, vector)

    Q2, G4, q24
    return G4, Q2, Symbol, T, q24


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    The first contraction forms the bracket in the component expression.
    Spenso orders the indices of a gamma matrix as `(row, column, Lorentz)`.
    Assigning the same label to compatible indices contracts them, so
    `Q2(1)` joins the row index and `G4(2)` the Lorentz index. We mark the
    remaining spinor index with `AUTO`, imported as `_`, to connect it to
    the propagator by multiplication in a later step.
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
    The same gamma matrix supplies the propagator numerator when its
    Lorentz index is contracted with $q_{24}$, while `q24 * q24` gives
    the denominator $s_{24}$. Here the label `3` is used for both the
    contracted Lorentz index and the open output spinor index. Since they
    belong to different spaces, these indices remain distinct; the first
    spinor index is again marked with `_` to receive the vertex contraction.
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
    Multiplication now connects the two spinor indices marked with `_`,
    completing the expression for $Q(2,4)$ with only its output index open.
    This current has a single vertex contribution. For the three-leg
    currents following equation (3.6), the two contributions would first
    be added and the propagator applied to their sum.
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
    Up to this point, `current` records the tensors and their contractions.
    Parsing it with `to_network()` turns this expression into an executable
    tensor network, using the built-in tensor library for the gamma matrices
    and symbolic components for the named inputs. The network retains the
    contractions and the scalar operations needed for the propagator
    denominator; these are evaluated in the next step.
    """)
    return


@app.cell
def _(current):
    network = current.to_network()
    network
    return (network,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Executing the network carries out the component sums and scalar
    operations, leaving the input components symbolic. Since `execute()`
    changes the network in place, we execute a copy so that the original
    network remains available for inspection.
    """)
    return


@app.cell
def _(network):
    network.execute()
    network
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    The completed network now contains the current's component expressions.
    `result_tensor()` extracts them as the tensor `Q24`, whose open spinor
    index labels the four expressions for $Q(2,4)_c$. Each includes the
    momentum-dependent propagator denominator and is ready for numerical
    evaluation.
    """)
    return


@app.cell
def _(network):
    Q24 = network.result_tensor()
    Q24
    return (Q24,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    These expressions can now be evaluated through the Symbolica interface
    introduced above. The tensor's `evaluator` method accepts the same
    options as `Expression.evaluator` and optimises all four components
    together. We use `components()` to list its inputs in the order
    $Q(2)$, $G(4)$, $q_{24}$; each numerical input row then produces a
    four-component tensor. As in the earlier scalar example, SymJIT
    compiles the evaluator on its first numerical use.
    """)
    return


@app.cell
def _(G4, Q2, Q24, q24):
    parameters = [*Q2.components(), *G4.components(), *q24.components()]
    evaluator = Q24.evaluator(parameters)

    value = evaluator.evaluate_complex([[1, 0, 0, 0, 0, 0, 1, 0, 2, 1, 0, 1]])[0]
    value
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    The numerical inputs above correspond to lightlike momenta
    $q_2=(1,0,0,1)$ and $q_4=(1,1,0,0)$, with
    $Q(2)=(1,0,0,0)$ and $G(4)=(0,0,1,0)$. They satisfy
    $Q(2)\not q_2=0$ and $q_4\cdot G(4)=0$, with an arbitrary overall
    spinor normalisation. Since $q_{24}=(2,1,0,1)$ and $s_{24}=2$, the
    result is $Q(2,4)=(-i/2,-i/2,0,0)$. Only the inputs to this local
    join are needed for the evaluation.

    Returning to the recurrence, $Q(2,4)$ supplies both
    $V_{dZd}[Q(2,4),Z(3)]$ in $Q(2,3,4)$ and
    $V_{dgd}[Q(2,4),G(5)]$ in $Q(2,4,5)$. The shared-current construction
    ensures that its value is computed once for these two uses. Each later
    join has its own component rule derived in the same way, until the
    final sum is contracted with the closing-leg spinor $v_{h_1}(q_1)$
    without a propagator, as in the amplitude expression following (3.6).
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Deriving numerical routines from interaction tensors is also the purpose
    of ALOHA [@deAquino:2011ub] and the UFO extension of Comix
    [@Hoeche:2014kca]. Spenso provides a generic contraction algorithm in
    which the index spaces and tensor entries are supplied separately.
    Thus the example can be expressed in a different Dirac basis by
    supplying the corresponding matrices and current components, while
    retaining the same contractions. The tensor library also provides
    $\gamma^5$, chiral projectors and fundamental or adjoint $SU(3)$
    representations; additional tensors can be defined through their
    symbolic or numerical entries. Before evaluating components, the
    symbolic tensor-algebra library Idenso can apply metric, Dirac and
    colour identities through `TensorExpression.simplify_algebra`.

    The open indices need not describe a single spinor or vector. In the
    momentum-independent four-point construction of the section on general
    contact decomposition, contracting two inputs leaves the intermediate
    tensor $X_{\rho\sigma}$ with two Lorentz indices. Spenso obtains its
    16 component expressions in the same way as the four expressions
    above. Symbolica then identifies exact zeros and components related by
    equality or a sign, giving the one stored representative for the
    metric-pair structure and the six for the antisymmetric four-gluon
    structure discussed there. In the four-gluon case, only the six entries
    with $\rho<\sigma$ need to be stored: the others follow from
    $X_{\rho\rho}=0$ and $X_{\sigma\rho}=-X_{\rho\sigma}$. pyAmpliCol
    records these relations so that the second join can contract with the
    remaining input current using the six stored entries, with exactly
    the same result as using all 16.

    For the supported higher-point colour-singlet interactions where no
    compact intermediate current is constructed, pyAmpliCol retains the
    input currents and momenta until the final join. Spenso derives the
    component expressions for that join by contracting the original UFO
    interaction with all its inputs. Thus pyAmpliCol specifies which
    intermediate currents to build, while Spenso supplies the expressions
    needed to evaluate the chosen joins. The section on general contact
    decomposition describes these constructions and their validity checks;
    the UFO-coverage table lists the interactions currently supported.
    """)
    return


@app.cell(hide_code=True)
def _():
    import marimo as mo

    return (mo,)


if __name__ == "__main__":
    app.run()
