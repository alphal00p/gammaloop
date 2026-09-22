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
   electron-positron to muon-pair benchmark and a different-flavor QCD quark
   annihilation benchmark check the sewn graph, physical cut, graph factors
   and uncut propagators. Other gallery processes remain pending.],
  [Dirac, color and Lorentz algebra], [Idenso and Spenso],
  [Existing gamma, color, metric, epsilon and adjoint operations. Keep their
   Symbolica expression interface; do not implement a second algebra in FeynKit.],
  [External spin sums], [`feynkit-generator::SpinSum`],
  [Scalar, Dirac and vector completeness tensors and external wavefunction-pair
   replacements. GammaLoop delegates both formulas and replacement construction
   here. Higher-spin states and polarized projectors
   are still missing.],
  [External color sums], [`feynkit-generator::ColorSum`],
  [Singlet, fundamental, sextet and adjoint completeness tensors reuse the
   generator's representation mapping and Spenso metrics. The Python particle
   API exposes optional initial-state averaging. Different-flavor QCD
   annihilation validates the generated massive SU(N) result.],
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
  [Symanzik polynomials use Symbolica determinant and cofactor operations.
   One- and two-loop reference polynomials and singular quadratic forms are
   validated. Singular output is algebraic data, not an integration formula.
   Complete parameter
   measures and tensor numerators remain outstanding. Polynomial equivalence
   alone does not establish contour or prescription equivalence.],
  [Scaleless sectors], [`feynkit-graph::IntegralFamily`],
  [Positive-power sector selection and exact parametric scaling certificates
   reuse the Symanzik polynomials and Symbolica linear solves. Automatic sector
   discovery, degenerate quadratic forms, UV/IR splitting and generator filter
   integration remain outstanding.],
  [IBP reduction], [Integration deferred],
  [22 available examples call Kira and four call FIRE. These are external
   dependencies used through FeynHelpers; they are not FeynCalc's own solvers.
   Integration is on hold at the maintainer's request pending upcoming IBP
   software; no Kira or FIRE adapter is being added.],
  [Analytic loop evaluation], [Shared OneLOop integration in the HEP host],
  [The HEP namespace exposes A0, B0, dB0, C0 and D0, including evaluable
   Symbolica Laurent coefficients. A generated massive photon self-energy
   validates the transverse form factor, UV pole and finite part through A0/B0.
   General reduction to these masters, higher epsilon orders and separate UV/IR
   bookkeeping remain outstanding.],
  [Cross sections and decay rates], [FeynKit kinematics and process APIs],
  [Shared symbolic/numerical initial-state flux and four-dimensional two-body
   phase space. The generated QED and QCD annihilation benchmarks include their
   angular distributions and total unpolarized cross sections. General phase space, identical-particle
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
reference uses the covariant projector. It reproduces physical polarizations
when unphysical states decouple; summing several external gluons covariantly
can require ghost subtraction, as demonstrated below. Supplying `reference=n`
selects the axial projector, including
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
integration, reproducing the unpolarized total cross section.

The live `symbolica-community/examples/hep_showcase.py` notebook also constructs
chiral-projected production directly from the generated vertices and propagators.
Its massive and massless results agree with the reference. Physical final-state
charges come from `cut.particles`: the finalized left side is the conjugate
amplitude, so the stored species of a cut edge need not be the outgoing species.
`cut.orientations` supplies the sign relating a selected loop coordinate to the
positive-energy outgoing momentum. Apply that sign when assigning physical
Mandelstam labels. These conventions leave native graph momentum routing intact;
changing which cut edge carries the loop coordinate must not interchange the
physical angular invariants.

Idenso now evaluates four-dimensional traces containing Spenso's `projp` and
`projm` by reducing them through its existing gamma-five trace identities. Use
`TensorExpression.projp(4)` or `TensorExpression.projm(4)` with `chain` to compose
projected currents. Symbolic-dimensional traces remain inert; this does not
choose a dimensional-regularization gamma-five scheme. Chirality specifies
helicity only in the massless limit.

== Generated Compton amplitudes

The #link("https://feyncalc.github.io/FeynCalcExamples/QED/Tree/ElGa-ElGa")[massive Compton reference]
is reproduced by generating the two ordinary tree amplitudes, summing them,
and contracting the operator with its `TensorExpression.dirac_adjoint` and
physical `Particle.spin_sum` tensors. The installed-host regression is
`crates/feynkit-py/tests/installed_feyncalc_compton.py`. It retains interference
between the two diagrams and averages only the incoming spin states. Electron
and positron amplitudes both give the same massive and massless result, with
both covariant photon sums and axial sums using the incoming fermion momentum
as a reference. The nonzero reference norm is retained. At $e = m_e = 1$,
$s = 3$ and $u = 0$, all four calculations give $3$.

External ports are aligned by matching generated wavefunctions with Symbolica;
the calculation does not assume internal half-edge numbers. Idenso retains
conjugation of scalar quantities explicitly, so the example declares its
physical momenta, charge, mass and invariants real before contracting the
adjoint. `wrap_indices` and `CookSettings.indices()` keep the adjoint's summed
indices separate. Dirac adjunction exchanges the input/output matrix roles;
fermion completeness tensors connect ket and Dirac-adjoint indices accordingly.
The physical completeness relation remains the one shared with GammaLoop.

== Massive Compton scattering: unresolved sewn-state convention

The ordinary-amplitude result does not validate conversion from a sewn forward
graph into a squared amplitude. Generation supplies four sewn
contributions for `e gamma -> e gamma`: both diagonal terms and both
interferences. Select the electron cut edge with
`with_loop_momentum_edges`, obtain its physical charge and momentum sign from
the cut, exclude cut denominators, and contract the initial-state spin sums.
This calculation reproduces the massless result
$-2 e^4 (s/u + u/s)$ for electrons and positrons.

The massive result currently fails when the physical `Particle.sum_spins`
tensor is applied directly to the native sewn projector. At
$e = m_e = 1$, $s = 3$, $t = -1$, $u = 0$, the reference squared matrix element
is $3$. The generated result is $-17$ with covariant photon sums and $1$ with
an axial incoming-photon sum whose reference is the incoming electron momentum.
The reference is timelike and has nonzero contraction with the photon, so these
are two admissible polarization sums for the same physical process.

Reversing only the mass term in the incoming fermion density matrix is a
diagnostic: it restores the complete massive reference expression and agreement
between both photon sums, for electrons and positrons. This is not a change to
the physical completeness relation. The missing boundary is between a physical
spin density and the sewn initial-state momentum convention. The unsquared
Compton amplitude routes its fermion propagator with the expected physical
momentum. Resolve the sewing boundary in shared code before treating this
sewn calculation as validated; do not compensate by changing a reference formula
or inserting a process-specific mass replacement into an example.

GammaLoop already owns the relevant workflow:
`CrossSectionGraph.apply_spin_sum` calls `ParticleTrait.polarization_sum`, which
now delegates to the shared `SpinSum`. Its original `GeneratePolarizations`
implementation selects wavefunctions from each half-edge's flow after sewing.
The physical completeness relation must remain shared and unchanged.

A comparison with GammaLoop's original `ParseGraph` sewing callback found that
FeynKit finalization reverses the incoming attachment, while the original
callback preserves it. A diagnostic restoring that callback, its vertex-slot
assignments, and the corresponding physical cut side restored covariant/axial
agreement for massive electron and positron Compton scattering. Both then gave
the negative of the complete reference expression, while annihilation retained
the correct sign. A common cut-phase multiplier cannot correct both results.

The external ordering calculation counts the closed fermion cycles created by
sewing. For $n$ open chains joined into $c$ cycles, the relative permutation
parity is $(-1)^(n-c)$, rather than the virtual-loop factor $(-1)^c$. An isolated
correction retaining the incoming attachment and including this chain-count
parity reproduces six exact full-mass comparisons: electron and positron Compton
scattering with both covariant and timelike axial incoming-photon sums,
electron-positron annihilation into muons, and electron-muon scattering. The
physical completeness tensors and the existing antifermion factor are unchanged.

This correction remains isolated from the production implementation pending
review of two native fixtures that explicitly encode the reversed-carrier
convention. With the correction, 142 of 144 graph/generator tests pass; the two
convention-dependent fixtures fail. Clippy passes with warnings denied. The
public Python regression also passes all six comparisons in a separate host
extension, including the massless Compton limit. The single-chain ordering
regression passes, as does GammaLoop's native-to-runtime parity regression for
both Compton charges. Two older Python checks still stop at hard-coded routing
assertions on both the current and corrected implementations; their previously
proposed routing-independent drafts pass the unpolarized angular distribution,
total cross section, and both polarized momentum choices without changing the
reference formulas. No separate sewing-specific spin-sum formula has been
introduced.

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
checks covariant moments through rank four in both explicit-index and compact-dot
notation, scalar prefactor preservation and an auxiliary null basis. Compact dots
between loop vectors and spectators are projected in the shared reducer without
a separate index-expansion step; loop invariants and declared external-basis
products stay scalar weights. Rust tests also cover two-direction Gram inversion,
mixed loop momenta, and a basis spanning the full Lorentz space. These tests
validate the projection identities, not the gallery's completed loop integrals.

== Generated photon self-energy

`crates/feynkit-py/tests/installed_feyncalc_photon_self_energy.py` runs in the
installed `symbolica.community.hep` host, where FeynKit and OneLOop share the
Symbolica runtime. It generates the electron loop, promotes Lorentz slots to
symbolic dimension before the Dirac trace, projects the tensor numerator and
rewrites scalar products with the diagram's integral family. Existing family
maps and vacuum projection check the shifted tadpole moments. This equal-mass
bubble needs only A0 and B0; no IBP adapter or second scalar evaluator is used.

With $s=p^2$, $D=4-2 epsilon$ and $tr(1)=4$, the result is
$(s g^(mu nu)-p^mu p^nu) F(D)$, where
$ F(D) = frac(2 e^2, s (D-1))
  (2 (D-2) A_0 - ((D-2)s+4m^2) B_0). $
The expressions omit the common loop factor $i/(16 pi^2)$.
The pole of $F$ is $-4 e^2/(3 epsilon)$, agreeing with the photon component of
#link("https://feyncalc.github.io/FeynCalcExamples/QED/OneLoop/Renormalization")[the gallery's QED renormalization example].
Keeping the dimension symbolic preserves the finite rational term from the
product of dimension-dependent coefficients and scalar poles.

Five finite-part checks use independently computed Feynman-parameter integrals,
including spacelike momenta, points on both sides of the pair threshold, and
unequal positive mass and renormalization scales. The above-threshold imaginary
part is checked separately against the two-particle cut. Native OneLOop
coefficient functions evaluate directly as Symbolica expressions. The benchmark
assumes positive squared mass and scale and nonzero $s$; it does not establish
the degenerate Gram limit or reproduce the lepton self-energy, vertex and
counterterm parts of the full renormalization example.

== Generated QCD annihilation

`Particle.color_sum(left, right, average=False)` constructs the identity in
that particle's color representation. Supply bare indices: the left slot carries
its representation and the right slot its dual. Antiquarks and antisextets
reverse the dual orientation. Singlets contribute one; averaging divides by
1, 3, 6 or 8 according to the UFO representation. Closing existing color slots
requires the dual slots. This operation uses the same representation mapping
as generated vertices and returns a standard Symbolica expression for Spenso
and Idenso to simplify. GammaLoop runtime color ports also use this table,
through `ColorRepresentation::from_ufo(...).representation()`. The shared
conversion returns a typed Spenso representation; source/sink flow selects
whether to dualize it. A runtime regression compares both flows with the
shared completeness tensors and checks their traces for every supported color
code, including singlets and conjugate sextets.

`installed_feyncalc_qcd_annihilation.py` generates
$b bar(b) -> t bar(t)$ through gluon exchange with both masses retained.
The cut supplies final-state completeness; `Particle.sum_spins` and
`Particle.color_sum` supply initial-state averages. No color algebra is
reimplemented in Python or GammaLoop. The full massive result agrees with the
#link("https://feyncalc.github.io/FeynCalcExamples/QCD/Tree/QiQibar-QjQjbar")[different-flavor FeynCalc example],
with $T_R=1/2$ and four-dimensional external spin states.

For symbolic SU(N), Spenso takes named fundamental and adjoint dimensions.
The relation $d_A=N_c^2-1$ is imposed after contraction. The massless result is
$ abs(cal(M))^2 = frac((N_c^2-1) g_s^4, 2 N_c^2 s^2) (t^2+u^2). $
At $N_c=3$, the existing flux and two-body phase-space APIs give
$ frac(d sigma, d Omega) = frac(alpha_s^2, 18s) (1+cos^2 theta), quad
  sigma = frac(8 pi alpha_s^2, 27s). $
The separate `hep/qcd_annihilation.py` notebook runs these exact comparisons
in the same Marimo instance as the other examples. Identical-flavor channels,
quark-gluon scattering and general QCD observables remain unvalidated.

== Generated quark annihilation into gluons

`installed_feyncalc_qcd_gluons.py` generates the three ordinary tree amplitudes
for $b bar(b) -> g g$, aligns their external ports, and forms the Dirac adjoint
of their sum. All interference terms are retained. Initial spin and color
averages use `Particle.spin_sum` and `Particle.color_sum`; final colors and
physical polarizations are summed with the same shared APIs.

The full massive SU(N) result is compared with
#link("https://feyncalc.github.io/FeynCalcExamples/QCD/Tree/QQbar-GlGl")[the FeynCalc two-gluon example].
Choosing the opposite final gluon or the incoming massive quark as each
polarization reference gives the same result, including the nonzero
reference-norm terms for the latter choice. The massless SU(3) limit is
$ abs(cal(M))^2 = frac(32 g_s^4, 27) frac(t^2+u^2, t u)
  - frac(8 g_s^4, 3) frac(t^2+u^2, s^2). $
The squared amplitude is symmetric under exchanging the labeled final gluons.
Integrating over both labels requires the usual $1/2!$ identical-particle factor;
this example does not supply a finite total massless cross section.

The calculation exposed two shared Idenso gaps: the existing color-conjugation
helper was disconnected from `spenso_conjugate` and `dirac_adjoint`, and color
metrics were not contracted inside collected traces before applying terminal
trace identities. Both paths now reuse their existing shared implementations.
Contracted symmetric traces through degree four reuse Spenso's normalized
projector expansion and Idenso's Casimir rules. Scalar representation labels
remain unchanged under conjugation. Higher-degree symmetric invariants remain
symbolic.

Compact fundamental generator chains and ordered color traces now conjugate
consistently with explicit networks, reversing the generator sequence and the
open endpoints. `installed_color_conjugation.py` checks SU(2), SU(3) and SU(5),
including the summed norm $(N_c^2-1)^2/(4 N_c)$ and conjugation twice. The live
`hep/color_algebra.py` notebook demonstrates the existing tensor constructors,
chain collection, conjugation and color contraction. Symmetric, antisymmetric
and cyclic generator groups, including nested groups and complex coefficients,
conjugate without expanding their permutations. Tests compare against Spenso's
explicit projector expansion and verify reversal parity for three-generator
traces. No separate FeynKit or GammaLoop conjugation implementation is introduced.

`hep/qcd_gluons.py` presents the generated diagrams and both gauge-reference
checks in a separate notebook on the same Marimo server.

== Fierz contractions with closed color traces

The existing shared Idenso Fierz implementation contracts generators between
open fundamental chains, between a chain and a trace, and between two traces.
A trace is cut at the contracted generator, preserving the cyclic order of its
remaining factors. This reduction runs before terminal trace decomposition so
that short traces do not hide a reducible contraction in symmetric invariants.
Only matching fundamental representations use this identity.

The #link("https://feyncalc.github.io/FeynCalcBook/SUNSimplify.html")[FeynCalc mixed trace example]
is covered by exact Rust and public Python regressions. They check both trace
evaluation settings, repeated simplification, and the option to keep separate
color lines. `installed_color_fierz.py` checks SU(2), SU(3), and SU(5), including
$ sum_(a,b,c) abs(op("Tr")(T^a T^b T^c))^2
  = frac((N^2-1)(N^2-2), 8N). $
An independent contraction with Spenso's numerical SU(3) matrices gives $7/3$.
`hep/color_fierz.py` presents the mixed identity and closed norm on the existing
Marimo server. No separate FeynKit or GammaLoop Fierz formulas are introduced.

== Generated symmetric color vertices

The shared generator accepts the UFO color tensor `d(1,2,3)`. It uses the same
adjoint-slot validation as `f` and lowers to four times Idenso's normalized
symmetric fundamental generator trace. This implements the
#link("https://feyncalc.github.io/FeynCalcBook/SUND.html")[standard symmetric SU(N) tensor]
with $T_R=1/2$, without adding another color algebra implementation.
The #link("https://link.springer.com/article/10.1140/epjc/s10052-023-11780-9")[UFO format]
assigns positive indices to vertex legs and negative indices to summed slots;
the existing shared index localization applies to `d` as well.

The generator regression verifies permutation symmetry, reality, the SU(3)
contractions $d^(a b c) d^(a b c)=40/3$ and
$d^(a b c) d^(a b e)=5/3 delta^(c e)$, a vanishing repeated-index contraction,
and rejection of an incompatible fundamental external slot.
`installed_ufo_symmetric_color.py` generates a cubic adjoint-scalar test vertex
through the public API and checks the three-port interface and exact norm.
`hep/symmetric_color_vertex.py` presents it in a separate live notebook.
The scalar records are reused from the stored model to define a toy adjoint
interaction; this is a color-factor validation.

Color epsilon and sextet interaction tensors remain unsupported by generation.
Their color-space completeness tensors do not establish interaction-tensor
coverage. The older GammaLoop reindexer also leaves those tensor heads
unlowered and rejects them.

== Covariant gluon sums and ghost subtraction

`installed_feyncalc_qcd_ghosts.py` validates
#link("https://feyncalc.github.io/FeynCalcExamples/QCD/Tree/QQbar-GlGl-2")[the ghost-subtraction variant]
using the same shared generation, color, spin and tensor APIs. Both ghost
orderings are generated independently, each with one diagram. Their external
color ports are aligned through Linnet half-edge metadata, since ghost legs
have no polarization wavefunctions.

With initial spin and color averages, each ghost contribution is
$ cal(G) = frac((N_c^2-1) g_s^4 (u-m^2)(t-m^2), 4 N_c s^2). $
The physical result is the covariant two-gluon square minus the two ghost
contributions. The exact massive SU(N) comparison, massless SU(3) limit and
Bose exchange all pass. The ghost correction is explicitly checked to be
nonzero: covariant gluon sums alone would give the wrong result for this
process. No extra ghost spin multiplicity or final-state symmetry factor is
inserted.

`hep/qcd_ghosts.py` displays the generated gluon and ghost diagrams, both ghost
contributions and the subtraction in its own notebook on the same server.

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
Singular-form regressions compare cofactor results with a nonsingular regulator
limit. Scaling and mapping APIs reject singular forms rather than interpreting
vanishing polynomials as a proof.

An example counts as reproduced only after exercising the relevant shared
components and comparing its final observable or symbolic identity with the
published result. Remaining work includes tree-level QED/QCD/EW observables,
polarized amplitudes, anomaly conventions, loop-renormalization examples,
external reduction jobs, and the two-loop topology-minimization example.
Full FeynCalc gallery parity is not yet achieved.
]
