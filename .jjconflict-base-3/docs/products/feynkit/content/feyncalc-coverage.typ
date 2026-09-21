#let feyncalc-coverage = [
= FeynCalc example coverage

This audit follows the #link("https://feyncalc.github.io/examples")[FeynCalc example gallery]
as retrieved on 21 September 2026. The machine-readable inventory is
`docs/products/feynkit/feyncalc-gallery.json`: 86 distinct example links, of which
85 were available. `EW/Tree/AnelEl-QubarQd` returned HTTP 404. The inventory
records function calls and requirements for every example; it deliberately does
not treat the presence of a primitive as an end-to-end validation.

== Shared ownership

#table(
  columns: (1fr, 1fr, 2fr),
  [Capability], [Owner], [Coverage and remaining work],
  [Models and diagrams], [`feynkit-model`, `feynkit-ufo`, `feynkit-generator`],
  [Existing validated models, UFO import and generation. A generated QED
   electron-positron to muon-pair benchmark checks the sewn graph, physical cut,
   graph factors and uncut propagators. Other gallery processes remain pending.],
  [Dirac, color and Lorentz algebra], [Idenso and Spenso],
  [Existing gamma, color, metric, epsilon and adjoint operations. Keep their
   Symbolica expression interface; do not implement a second algebra in FeynKit.],
  [External spin sums], [`feynkit-generator::SpinSum`],
  [Scalar, Dirac and vector completeness tensors and external wavefunction-pair
   replacements. GammaLoop delegates both formulas and replacement construction
   here. Higher-spin states and polarized projectors
   are still missing.],
  [Symbolic kinematics], [`feynkit-kinematics::Kinematics`],
  [Scoped scalar products in integer or symbolic dimensions and four-dimensional
   unequal-mass Mandelstam substitutions. Scalar products expand bilinearly in
   declared momenta using Symbolica. General momentum elimination still needs
   coverage.],
  [Covariant tensor reduction], [`feynkit-tensor`],
  [Symmetry-aware vacuum projection and external-basis reduction. Symbolica
   inverts the external Gram matrix; the existing projector averages the
   transverse components. Integral-family reduction remains a separate step.],
  [UV expansion], [`feynkit-graph` and Vakint],
  [Graph and subgraph expansion and vacuum-integral infrastructure exist.
   The gallery's complete renormalization constants are not yet validated.],
  [Integral families and mappings], [`feynkit-graph::IntegralFamily`],
  [Generated diagrams expose families through their shared denominator builder
   and momentum routing. Affine propagator rank, partial fractions, scalar-product completion and
   numerator mappings use Symbolica linear algebra. Verified affine momentum
   shifts and subtopology embeddings are available, including automatic search
   from quadratic propagators. Symanzik incidence graphs support parameter
   permutations beyond fixed-external-momentum shifts. Family discovery and
   subtopology minimization remain outstanding.],
  [Feynman parameters], [`feynkit-graph::IntegralFamily`],
  [Symanzik polynomials use Symbolica determinant and inverse operations.
   One- and two-loop reference polynomials are validated. Complete parameter
   measures and tensor numerators remain outstanding. Polynomial equivalence
   alone does not establish contour or prescription equivalence.],
  [Scaleless sectors], [`feynkit-graph::IntegralFamily`],
  [Positive-power sector selection and exact parametric scaling certificates
   reuse the Symanzik polynomials and Symbolica linear solves. Automatic sector
   discovery, degenerate quadratic forms, UV/IR splitting and generator filter
   integration remain outstanding.],
  [IBP reduction], [External reducer adapters needed],
  [22 available examples call Kira and four call FIRE. These are external
   dependencies used through FeynHelpers; they are not FeynCalc's own solvers.
   Export, invocation and result import still need implementation.],
  [Analytic loop evaluation], [Shared integration components needed],
  [General non-vacuum one-loop scalar masters, epsilon expansions and UV/IR
   separation remain outstanding. Reuse existing integration packages where
   their conventions and capabilities match.],
  [Cross sections and decay rates], [FeynKit kinematics and process APIs],
  [Shared symbolic/numerical initial-state flux and four-dimensional two-body
   phase space. The generated QED benchmark includes its angular distribution
   and total unpolarized cross section. General phase space, identical-particle
   bookkeeping and the remaining gallery observables need further coverage.],
)

Among the available pages, 49 use fermion or vector spin sums, 46 use
Mandelstam or scalar-product expansion helpers, 55 use color simplification,
76 use Dirac simplification, 35 use non-vacuum tensor reduction, and 30 use
integral-family mappings. These counts are requirements, not passed tests.
The inventory recognizes bracketed, prefix and postfix Mathematica calls.

== Use the shared external-state tensors

`Particle.spin_sum(momentum, left, right)` returns an ordinary Symbolica
expression. Supply bare index symbols and an unindexed momentum name. The
fermion result uses Spenso bispinor slots; the vector result uses Minkowski
slots. `average=True` divides by the physical four-dimensional spin count.
The antiparticle record selects the negative mass term in the Dirac
completeness relation.

Massive vector sums use the Proca projector. A massless vector without a
reference uses the covariant projector, appropriate for a gauge-invariant
amplitude. Supplying `reference=n` selects the axial projector, including
the term proportional to the reference's squared norm. Its scalar product
with the external momentum must be nonzero. `covariant=True` explicitly
selects the Feynman-gauge numerator for a massive vector as well.

For generated sewn graphs, `Particle.sum_spins(expression, momentum, edge=id)`
applies the completeness relation directly to the matching wavefunction pair.
Use `diagram.projector_expression()` as the starting expression and apply the
sum once for each external edge. `average=True` supplies initial-state spin
averaging. Other edge labels and unpaired wavefunctions remain unchanged;
this operation does not conjugate an amplitude or sum color. Both FeynKit and
GammaLoop use the same replacement construction.

The generator retains provenance annotations for graph signs and multiplicities.
`diagram.overall_factor_expression(evaluate=True)` evaluates these annotations
through the shared graph implementation. It preserves arbitrary symbolic
factors. Include this weight once when assembling a squared matrix element;
`diagram.symmetry_factor` alone omits fermion signs and grouping weights.

Idenso's `to_dots` handles explicit vector slots even when the momentum head
was created as an ordinary Symbolica symbol without tensor tags. GammaLoop
uses the same contraction operation.

== Keep kinematics local to a calculation

`Kinematics.mandelstam([p1, p2, p3, p4], masses_squared, [s, t, u])` uses
two incoming and two outgoing momenta. It sets the four on-shell products
and six pair products without modifying Symbolica's global state. Apply
the context to the output of Idenso contractions with `kinematics.apply(expr)`.
The same context recognizes Spenso dot products and metric shorthand.

The mass arguments are squared masses. The invariants satisfy
`s + t + u = sum(masses_squared)`. Eliminate a chosen invariant using a
normal Symbolica substitution. Separate contexts can describe separate
processes without clearing a global scalar-product table.

`kinematics.scalar_product(p1 + p2, p1 + p2)` expands to `s`. Setting an
assumption declares the two momentum names. For combinations without assumptions,
use `Kinematics(momenta=[p, q])`; other symbols in the combination are scalar
coefficients. For example, `scalar_product(a*p + q, p)` produces
`a*(p.p) + q.p`. Nonlinear momentum expressions are rejected instead of being
interpreted as vectors. Labeled names such as `Q(1)` work in the same way.

For dimensional regularization, construct `Kinematics(D)` with a Symbolica
dimension symbol and set individual products with `with_scalar_product`.
These assumptions do not match four-dimensional slots. Substitute
`D = 4 - 2 eps` after contraction, since Spenso dimensions are symbols or
nonnegative integers rather than compound expressions.

The installed-host test `crates/feynkit-py/tests/installed_feyncalc_tree.py`
combines particle spin sums, Idenso Dirac traces, and Mandelstam substitutions
to reproduce the massive and massless unpolarized and chiral-projected
squared-current results for
#link("https://feyncalc.github.io/FeynCalcExamples/QED/Tree/ElAel-MuAmu")[electron-positron annihilation to muons].
The separate `installed_feyncalc_generated_qed.py` regression starts from the
Standard Model fixture, restricts generation to the electron-photon and
muon-photon vertices, and reproduces the same massive and massless squared
matrix elements. The sewn graph has one loop, while each amplitude side has
zero loops. Incoming wavefunctions are summed through `Particle.sum_spins`;
cut fermion numerators already contain the final-state completeness factors.
The test excludes cut denominators from the squared amplitude and retains the
uncut photon denominators and complete graph weight. Its massless result is
integrated with the two-body phase-space measure and Symbolica polynomial
integration, reproducing the unpolarized total cross section. Polarized
production from the generated sewn graph remains a separate check before marking
the full example complete: inserting the same right projector at every generated
vertex gives the opposite angular correlation to the physical current benchmark.
The external-state and cut-flow conventions still need to be reconciled.

Idenso now evaluates four-dimensional traces containing Spenso's `projp` and
`projm` by reducing them through its existing gamma-five trace identities. Use
`TensorExpression.projp(4)` or `TensorExpression.projm(4)` with `chain` to compose
projected currents. Symbolic-dimensional traces remain inert; this does not
choose a dimensional-regularization gamma-five scheme. Chirality specifies
helicity only in the massless limit.

== Normalize and integrate two-body observables

`Kinematics.flux(p1, p2)` returns the denominator
$4 sqrt((p_1 dot p_2)^2 - m_1^2 m_2^2)$ of an invariant squared amplitude.
`kin.two_body_phase_space(k1, k2)` returns the four-dimensional
$d Phi_2 / d Omega$ in the final pair's rest frame. The measure includes the
$(2 pi)^4$ multiplying the momentum-conservation delta function. Their ratio
multiplies the squared matrix element to give a differential cross section in
natural units. These conventions agree with the
#link("https://pdg.lbl.gov/2025/reviews/rpp2025-rev-kinematics.pdf")[PDG kinematics review].

Supply future-directed on-shell momenta in the physical region. Declare positive
invariant symbols when known, so Symbolica can resolve square-root branches.
A symbolic measure is not a numerical phase-space generator and does not insert
threshold step functions. Spin/color averages and final-state identical-particle
factors remain explicit; do not apply a factorial already included in a generated
graph weight a second time.

```python
from symbolica import E, Expression, S
from symbolica.community import feynkit as fk

p1, p2, k1, k2, t, u, c, alpha = S(
    "p1", "p2", "k1", "k2", "t", "u", "cos_theta", "alpha"
)
s = S("s", is_positive=True)
pi = Expression.PI
kin = fk.Kinematics.mandelstam(
    [p1, p2, k1, k2], [E("0")] * 4, [s, t, u]
)
squared = 2 * (4*pi*alpha)**2 * (t**2 + u**2) / s**2
angular = squared.replace(t, -s*(1-c)/2).replace(u, -s*(1+c)/2)
differential = (angular * kin.two_body_phase_space(k1, k2) / kin.flux(p1, p2)).expand()
primitive = differential.to_polynomial().integrate(c).to_expression()
total = 2*pi*(primitive.replace(c, E("1")) - primitive.replace(c, E("-1")))
assert (total - 4*pi*alpha**2/(3*s)).together() == E("0")
```

Polynomial integration is available directly in Symbolica and suffices for this
angular dependence; it does not require an optional general symbolic integrator.
The installed generated regression obtains `squared` from the graph instead of
entering it by hand.

For a decay, `kin.flux(parent)` returns the rest-frame denominator $2 M$.
A constant squared amplitude into two distinct massless particles then gives
$Gamma = abs(cal(M))^2/(16 pi M)$. Numerical `FourMomentum.flux()` instead uses
the supplied frame's energy, returning $2 E$; `p.flux(q)` is Lorentz invariant.
Both numerical and symbolic APIs, and GammaLoop's integrand, use the same
`feynkit-kinematics::InitialStateFlux` implementation. GammaLoop retains its
runtime model masses and barn conversion outside that shared boundary.

== Reduce tensors with external momentum dependence

`TensorReducer(D).with_integrated_vector(k(mink(D)))` selects a loop vector.
Add independent denominator directions with `with_external_vector(p(mink(D)))`.
The reducer retains longitudinal components and applies its existing vacuum
projection only in the transverse space. This supports free indices, multiple
loop vectors and odd total ranks without a second projector implementation.
See the #link("guides/tensor-reduction/")[tensor-reduction guide] for a complete
example and the treatment of null external directions.

The installed-host test `crates/feynkit-py/tests/installed_feyncalc_tensor.py`
checks covariant moments through rank four, scalar prefactor preservation and
an auxiliary null basis. Rust tests also cover two-direction Gram inversion,
mixed loop momenta, and a basis spanning the full Lorentz space. These tests
validate the projection identities, not the gallery's completed loop integrals.

== Validation standard

The #link("guides/integral-families/")[integral-family guide] describes the shared
rank and completion APIs. Rust and installed-host tests cover a massive bubble,
a two-loop incomplete family, eikonal forms, dependent propagators and invalid
momentum declarations. These are algebraic checks before external IBP reduction.
`installed_diagram_integral_families.py` generates a massless scalar bubble,
extracts families in each routing, compares the Symanzik polynomials and checks
the on-shell scaleless limit. A Rust regression covers a massive bubble.
`installed_partial_fractions.py` additionally reconstructs the original rational
functions and verifies independent propagator support after decomposition,
including homogeneous relations and repeated powers.
Mapping regressions include mixed loop bases, eikonal shifts, subtopology power
maps, and a published three-loop example. These do not establish full topology
minimization or the gallery's complete two-loop IBP workflows.
Parameter-space regressions cover a published family pair that requires an
external-momentum transformation for an affine map, all permutations of a
massive vacuum family, and rejection of changed masses. These compare both
Symanzik polynomials and exercise the same shared power-reordering implementation.

An example counts as reproduced only after exercising the relevant shared
components and comparing its final observable or symbolic identity with the
published result. Remaining work includes tree-level QED/QCD/EW observables,
polarized amplitudes, anomaly conventions, loop-renormalization examples,
external reduction jobs, and the two-loop topology-minimization example.
Full FeynCalc gallery parity is not yet achieved.
]
