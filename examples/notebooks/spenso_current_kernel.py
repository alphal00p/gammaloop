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

    A UFO Lorentz rule describes an interaction through its tensor indices.
    To turn it into a current kernel, pyAmpliCol selects an output leg and
    supplies currents for the others. Spenso contracts their indices and
    leaves those of the output leg open.

    Each component of the resulting current is a scalar symbolic expression
    in the input components, momenta and model parameters. The spinor current
    $Q(2,4)_c$ below has four such expressions, one for each value
    $c=0,1,2,3$. All contracted indices have been summed over. Symbolica can
    then optimise these expressions and turn them into numerical kernels.

    The purpose is similar to that of ALOHA [@deAquino:2011ub] and the UFO
    extension of Comix [@Hoeche:2014kca]. Spenso [@SpensoSoftware] separates
    the index spaces and tensor components from the contraction algorithm.
    Changing the basis of the Dirac matrices, for example, requires supplying
    the corresponding tensor and current components; the contraction code
    remains the same.

    Each index belongs to a specified representation: Minkowski vectors for
    Lorentz indices, four-component bispinors for Dirac indices, and
    fundamental or adjoint $SU(3)$ representations for colour. The HEP library
    provides the metric, Dirac matrices, $\gamma^5$ and chiral projectors.
    For tensors absent from the library, Spenso generates symbolic components;
    users may instead supply their own symbolic or numerical entries.
    Before component evaluation, Idenso's `simplify_algebra` method on
    `TensorExpression` can apply metric, Dirac and colour identities.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Consider the quark–gluon join $Q(2,4)=P_d V_{dgd}[Q(2),G(4)]$ from
    equation (3.6) in the section on off-shell-current recurrence. We will
    build it directly with the tensor API. For this standalone example,
    we use a four-component row spinor, a massless quark propagator and the
    vertex normalization $V^\mu=i\gamma^\mu$, omitting the colour matrix and
    strong coupling. With $q_{24}=q_2+q_4$ and $s_{24}=q_{24}^2$,

    $$
    Q(2,4)_c = \frac{i}{s_{24}}\sum_b
    \left[i\sum_{a,\mu}Q(2)_a(\gamma^\mu)_{ab}G(4)_\mu\right]
    (\not q_{24})_{bc}.
    $$

    The bracket contains the vertex contraction. The propagator acts on its
    right because the current is a row spinor; it supplies the second factor
    of $i$. We use the library's Weyl gamma matrices and the metric
    $g=\mathrm{diag}(1,-1,-1,-1)$. Vector components are supplied
    contravariantly, and Minkowski contractions provide the metric signs.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    We begin with the two input currents and their combined momentum.
    `TN.vector` declares a tensor with one index, whose space is set by the
    representation. The integers preceding that representation identify the
    external legs: `2` for $Q(2)$, `4` for $G(4)$, and `2, 4` for $q_{24}$.
    These labels are part of the tensor's name. A later call such as `Q2(1)`
    assigns its tensor index, here the spinor label `1`.
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
    The vertex joins the spinor and gluon currents. Its gamma slots follow
    the order `(row, column, Lorentz)`: index `1` contracts with `Q2`, and
    index `2` with `G4`. The remaining spinor slot, marked by `AUTO` as `_`,
    stays unlabelled so that it can connect to the propagator.
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
    Next comes the propagator, $i\not q_{24}/q_{24}^2$. The product
    `q24 * q24` forms the Minkowski scalar product. In the numerator,
    the Lorentz index `3` contracts the gamma matrix with the momentum,
    while the spinor index `3` remains open. The same label can serve both
    roles because the indices belong to different spaces. The unlabelled
    first spinor slot will receive the join.
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
    Multiplying the join by the propagator connects their unlabelled spinor
    slots and leaves the output index open. If several vertex contributions
    fed this current, we would add them before applying the propagator once.
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
    The expression is now complete. Evaluating its contractions gives four
    scalar component expressions, including the propagator denominator.
    They form the tensor `Q24`, with its one open bispinor index.
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
    These expressions become the numerical kernel. `components()` returns
    the input symbols in component order, ready to serve as evaluator
    parameters. The denominator is already expressed in the momentum
    components. Once built, the evaluator can be reused at new input values
    without repeating the tensor contraction.
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
    For the evaluation above, take $q_2=(1,0,0,1)$, $q_4=(1,1,0,0)$,
    $Q(2)=(1,0,0,0)$ and $G(4)=(0,0,1,0)$. The input momenta are lightlike,
    the spinor satisfies $Q(2)\not q_2=0$, and the gluon polarization obeys
    $q_4\cdot G(4)=0$. The spinor's overall normalization is arbitrary.
    Their combined momentum is $q_{24}=(2,1,0,1)$, with $s_{24}=2$, and the
    current evaluates to $Q(2,4)=(-i/2,-i/2,0,0)$. This specifies the local
    join; the remaining legs of the scattering process are not needed here.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    The same construction extends to currents with more open indices.
    When pyAmpliCol decomposes a momentum-independent four-point interaction
    into successive joins, as described in the section on general contact
    decomposition, the intermediate current may carry two Lorentz indices.
    Contracting the first two inputs gives the tensor $X_{\rho\sigma}$,
    with 16 components before reduction. Symbolica identifies exact zeros
    and equality or sign relations among them. For the pairings discussed
    there, one independent component remains for
    $g_{\mu\nu}g_{\rho\sigma}$, and six for the combination in the
    antisymmetric-contact equation.

    The division of work remains the same: pyAmpliCol chooses the joins,
    Spenso contracts the tensors, and Symbolica simplifies and compares
    the component expressions. Supported higher-point colour-singlet rules
    follow this route too. Where no reduced tensor factorisation is used,
    the original interaction is contracted directly with its input currents.
    The resulting expressions supply kernels for recurrence and compiled
    execution. The symbolic work takes place during model processing and
    kernel construction, before ordinary double-precision event evaluation.

    These tools do not by themselves establish model support. Unresolved
    UFO tensor functions are rejected. For momentum-independent four-point
    rules, pyAmpliCol records the component ordering, zeros and exact sign
    relations for each inequivalent output leg. Coloured contact rules must
    also match an explicitly reconstructed colour class from the section on
    general contact decomposition. Adding a Lorentz structure requires its
    components and contraction rules; adding a colour representation, such
    as a sextet, also requires the associated Clebsch–Gordan tensors, colour
    flows, projections and contractions. Further restrictions, including
    Majorana fermion flow and unverified coloured high-point contacts,
    appear in the UFO-coverage table.

    Within these boundaries, the construction supports multiple open quark
    lines and LC, NLC and full-colour calculations.
    """)
    return


@app.cell(hide_code=True)
def _():
    import marimo as mo

    return (mo,)


if __name__ == "__main__":
    app.run()
