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
   and uncut propagators. Ordinary generated amplitudes also validate Compton,
   diphoton, Bhabha, Møller, selected QCD channels and on-shell Higgs and Z
   decays. The inventory records the validated scope of each process.],
  [Dirac, color and Lorentz algebra], [Idenso and Spenso],
  [Existing gamma, color, metric, epsilon and adjoint operations. Keep their
   Symbolica expression interface; do not implement a second algebra in FeynKit.],
  [External spin sums], [`feynkit-generator::SpinSum`],
  [Scalar, Dirac and vector completeness tensors and external wavefunction-pair
   replacements. GammaLoop delegates both formulas and replacement construction
   here. Massive Dirac density matrices select a physical spin vector through
   the same owner. Higher-spin states and massless helicity projectors remain
   outstanding.],
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
   The one-loop QED and QCD renormalization workflows combine generated
   self-energies, fermion vertices and the ghost-gluon vertex with symbolic gauge dependence and
   counterterm linear solves. Generated two-loop massless electron and scalar
   self-energies validate bare UV poles with supplied counterterm sums. The
   two-loop photon example retains symbolic gauge dependence and calculates
   its four one-loop counterterm insertions using signed propagator powers.
   Analytic vacuum values remain explicit inputs.
   Automatic forest generation stays separate.],
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
  [IBP reduction], [RustRed's native HEP bridge],
  [`rustred-feynkit` is now integrated into the Symbolica Community host.
   `hep.IBPFamily` consumes the existing FeynKit family and exposes symbolic
   identities, bounded Laporta elimination and parametric recurrences. Residual
   integrals at a finite search depth are not certified masters. Two-loop
   scalar and Feynman-gauge massless electron self-energies, symbolic-gauge photon
   renormalization, the electron Pauli form factor and unequal-mass bubble
   regressions pass in the installed host.
   Separate notebooks run on the existing Marimo instance. Of the gallery pages,
   22 call Kira and four call FIRE through FeynHelpers; those external
   interfaces are not FeynCalc-owned solvers.],
  [Analytic loop evaluation], [Shared OneLOop integration in the HEP host],
  [The HEP namespace exposes A0, B0, dB0, C0 and D0, including evaluable
   Symbolica Laurent coefficients. A generated massive photon self-energy
   validates the transverse form factor, UV pole and finite part through A0/B0.
   The generated electron vertex also reproduces the Pauli form factor after
   native IBP reduction. General reduction to these masters, higher epsilon orders
   and separate UV/IR bookkeeping remain outstanding.],
  [Cross sections and decay rates], [FeynKit kinematics and process APIs],
  [Shared symbolic/numerical initial-state flux and four-dimensional two-body
   phase space. The generated QED and QCD annihilation benchmarks include their
   angular distributions and total unpolarized cross sections. Chiral Z decays
   validate massive two-body widths for all four fermion classes; Higgs decays
   cover charged leptons, quarks, WW and ZZ. General phase space, identical-particle
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

For a massive Dirac fermion, `spin_vector=s` selects one physical spin state
in either `spin_sum` or `sum_spins`:
$ rho_u = frac(1,2) (slash(p)+m)(1+gamma^5 slash(s)), quad
  rho_v = frac(1,2) (slash(p)-m)(1+gamma^5 slash(s)). $
Here $s$ is the dimensionless boosted rest-frame spin direction, with
$p dot s=0$ and $s^2=-1$; the caller supplies these on-shell constraints to
`Kinematics`. Particle and antiparticle use the same physical spin direction,
as in the #link("https://sites.ualberta.ca/~gingrich/courses/phys512/node61.html")[covariant spin-projector convention].
A selected state cannot also be spin averaged. The model must declare a
nonzero mass: a zero-valued mass parameter is treated as massless even when
its symbolic parameter name exists. This does not implement massless helicity
projectors or dimension-generic gamma-five schemes.

`installed_polarized_spin_density.py` checks both charges, opposite-spin
completeness, traces, rank-one purity, left and right Dirac equations, and
edge-specific wavefunction replacement. A GammaLoop regression compares all
128 complex density-matrix entries across two momentum directions, both
helicities and both charges with its existing numerical spinors. It uses the
shared HEP tensor library and adds no second spinor implementation.
`hep/polarized_spin.py` displays the density matrix and constructs the spin
vector with the existing shared `Boost` and `FourMomentum` APIs.

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

The older `installed_feyncalc_polarized_qed.py` fixture still imports the former
`community.feynkit` namespace and assumes a positive cut orientation for either
loop coordinate. It does not pass unchanged in the current host. A scratch
correction using `community.hep` and the physical momentum sign passes both
massive and massless reference assertions for both coordinates; the repository
fixture update remains subject to maintainer approval.

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

== Generated Bhabha and Møller amplitudes

The #link("https://feyncalc.github.io/FeynCalcExamples/QED/Tree/ElAel-ElAel")[Bhabha]
and #link("https://feyncalc.github.io/FeynCalcExamples/QED/Tree/ElEl-ElEl")[Møller]
workflows generate both ordinary tree amplitudes, retain their graph signs and
propagators, average incoming spins and sum outgoing spins. The massive square
includes interference between the two channels, with $s+t+u=4m_e^2$.
`installed_feyncalc_identical_leptons.py` compares the massive reference
expressions, their massless limits, separate interference terms, and Møller's
final-electron exchange symmetry.

`TensorExpression.dirac_adjoint(preserve_indices=True)` keeps labels attached
to the same physical external legs. Unlike a matrix-adjoint endpoint convention,
it distributes over sums whose diagrams pair the external fermions differently.
The complete amplitude can therefore be conjugated in one call before applying
shared particle spin sums. Conjugation, gamma-zero boundary factors and index
handling remain in Idenso; the example needs no channel-specific endpoint swaps.
The default `preserve_indices=False` retains the existing matrix convention and
factored output. Shared Rust regressions verify conjugated complex coefficients,
additivity, involution and the default behavior. The installed-host regression
passes both complete massive references, the interference checks, massless
limits and angular-cut cross sections through this shared adjoint path.

For equal external masses, the shared two-body measure divided by the flux is
$1/(64 pi^2 s)$. The massless event densities, in units of $alpha^2/s$ and with
$x=cos theta$, are $(3+x^2)^2/(4(1-x)^2)$ for Bhabha scattering and
$(3+x^2)^2/(2(1-x^2)^2)$ for Møller scattering. The latter includes $1/2!$ on
the full sphere; using the labeled density on one hemisphere gives the same
event count. This factor does not remove either exchange amplitude.
The regression integrates a symmetric angular cut with Symbolica and compares
the two Møller counting conventions. `hep/identical_leptons.py` displays both
diagrams, the massive square, interference and angular densities, with reaction,
CM-speed and scattering-angle controls. It retains the ordinary labeled
amplitudes; external-fermion representative symmetrization transfers permutation
sign handling to the caller and is not used in this calculation. Five live
control configurations pass for both reactions, including equal Møller results
at opposite scattering cosines with a nonzero electron mass.

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

== Generated electron–positron annihilation into photons

`installed_feyncalc_diphoton.py` generates both tree diagrams for
$e^- e^+ -> gamma gamma$ and retains their interference. The full-mass
spin-averaged result agrees with the
#link("https://feyncalc.github.io/FeynCalcExamples/QED/Tree/ElAel-GaGa")[FeynCalc diphoton example].
It is unchanged by choosing covariant photon sums, the other photon as a null
reference, or the incoming electron as a timelike reference. The spin-summed
square vanishes when either photon polarization tensor is replaced by its
longitudinal momentum product. The massless limit and exchange of the two
photons are checked separately.

The ordinary amplitude labels both photons. Its spin-averaged square contains
no final-state factorial; neither does `Kinematics.two_body_phase_space`.
Integrating over both photon labels therefore requires $1/2!$. Equivalently,
select the forward photon and integrate its angle over a single hemisphere.
For $abs(cos(theta)) < C < 1$, the massless event cross section is
$ sigma = frac(2 pi alpha^2, s) (ln frac(1+C, 1-C) - C). $
This also agrees with the
#link("https://arxiv.org/pdf/hep-ex/0409058")[DELPHI Born cross section, Eq. (2)].
Symbolica integrates the rational angular distribution directly. The example
also integrates the massive distribution at incoming speed
$ beta = sqrt(1 - 4 m_e^2/s) $, checks the threshold normalization and recovers
the massless limit. It selects the positive above-threshold flux branch
explicitly. `hep/diphoton.py` offers polarization-reference, angular-cut and
incoming-speed controls on the existing Marimo server.

The massive sewn-forward conversion remains unresolved. At
$m_e=e=1$, $s=10$, $t=-1$ and $u=-7$, its Bose-completed result is $59/8$ with
covariant photon sums and $265/16$ with an incoming-electron reference; the
ordinary-amplitude reference is $83/8$. Both sewn massless limits agree with
the reference. This gauge dependence cannot be repaired by an
identical-particle normalization factor and supplies another regression target
for the pending shared sewing correction. Existing graph weights must not be
multiplied by a second inverse automorphism factor.

== Generated chiral Z decays

`installed_feyncalc_z_decay.py` reproduces all four classes in the
#link("https://feyncalc.github.io/FeynCalcExamples/EW/Tree/Z-FFbar")[FeynCalc Z-decay example]:
neutrinos, charged leptons, up-type quarks and down-type quarks. It generates
ordinary amplitudes, retains the massive chiral interference, averages the
three initial Z polarizations with the full Proca projector, and sums final
spins and colors through the shared particle APIs.

The calculation exposed missing special-matrix conjugation in Idenso.
The existing conjugation operation now handles four-dimensional gamma-five
and chiral projectors with its gamma-zero machinery. The physical adjoints
are $overline(gamma^5)=-gamma^5$ and $overline(P_L)=P_R$.
Exact Rust identities and 64 numerical matrix components check the adjoints.
Dimension-generic gamma-five conventions remain unspecified.

The shared two-body measure and rest-frame decay flux give
$ Gamma = frac(N_c G_F M^3 beta, 6 pi sqrt(2))
  (c_V^2 (1+2r) + c_A^2 (1-4r)), $
where $r=m_f^2/M^2$, $beta=sqrt(1-4r)$,
$c_V=T_3-2Q_f sin^2(theta_W)$ and $c_A=T_3$.
The fermion and antifermion are distinct. Tests check the positive
above-threshold phase-space branch, massless limits, and vector versus axial
threshold powers. `hep/z_decay.py` offers all four channels in a separate
notebook on the same server.

== Generated massive Higgs decays

`installed_feyncalc_higgs_decay.py` reproduces the gallery's
#link("https://feyncalc.github.io/FeynCalcExamples/EW/Tree/H-FFbar")[fermion],
#link("https://feyncalc.github.io/FeynCalcExamples/EW/Tree/H-WW")[WW] and
#link("https://feyncalc.github.io/FeynCalcExamples/EW/Tree/H-ZZ")[ZZ] decays.
All five representative channels use generated amplitudes, model parameter
expressions, Idenso adjoints, shared particle spin/color sums, and the shared
two-body phase space and decay flux. No extra amplitude evaluator is introduced.

The fermion calculation explicitly identifies each UFO Yukawa mass with its
pole mass for this tree-level comparison. The vector channels retain all
three physical Proca polarizations. Only the identical ZZ final state receives
an explicit phase-space factor of $1/2!$; it is not part of the labeled
amplitude or spin sum. The tests compare the fully massive squared amplitudes
and widths. `hep/higgs_decay.py` exposes all five channels in a separate
notebook on the existing server.

These are above-threshold, on-shell two-body decays. The WW and ZZ formulas
do not cover off-shell vector decays at the physical 125 GeV Higgs mass.

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
the degenerate Gram limit. The one-loop renormalization workflow below adds
the lepton self-energy, vertex UV pole and counterterm matching. Complete finite
lepton self-energy and renormalized vertex form factors remain outside that
workflow.

== Generated one-loop QED renormalization

`hep/qed_renormalization.py` and
`installed_feyncalc_qed_renormalization.py` assemble the local UV structures of
#link("https://feyncalc.github.io/FeynCalcExamples/QED/OneLoop/Renormalization")[the massive one-loop QED renormalization example].
They generate the electron and photon self-energies, electron-photon vertex,
and its tree normalization from the model. The internal photon propagator
retains the symbolic covariant-gauge parameter $xi$; the electron mass stays
symbolic throughout. Generated wavefunctions identify the open tensor ports.
Only the external Wick-order sign is removed from the amputated electron
self-energy; internal fermion-loop signs and other graph weights are retained.
The vertex is normalized against the generated tree vertex.

`diagram.uv_expansion` retains every term through the graph's UV degree.
Idenso traces and the shared vacuum `TensorReducer` project the local
structures, and `IntegralFamily` rewrites them in the common vacuum denominator
$k^2-M$, where $M=m_"UV"^2$. Native `IBPFamily.reduce_laporta` reduces raised
powers to the tadpole. Its analytic pole $A_0(M)=M/epsilon+O(1)$ is a supplied
input, checked with OneLOop; it is not inferred by the IBP solver. Restoring the
loop measure $i (4 pi)^(epsilon-2)$ and using $a_4=e^2/(16 pi^2)$ fixes the
normalization. Lorentz contractions keep $D$ symbolic and use $tr(1)=4$;
substitute $D=4-2 epsilon$ after reduction.

The reference structures compared by the calculation are
$ Sigma_"UV" = frac(i a_4,epsilon)
  (xi slash(p)-(xi+3)m), quad
  Pi_"UV"^(mu nu) = -frac(4 i a_4 N_f,3 epsilon)
  (p^2 g^(mu nu)-p^mu p^nu), $
and the vertex pole is $a_4 xi/epsilon$ times its tree value. The photon loop
receives a symbolic flavor multiplicity $N_f$. Physical-mass and auxiliary-mass
correction terms in the UV expansion cancel local $m^2 g^(mu nu)$ and
$M g^(mu nu)$ poles; the auxiliary mass does not survive in the final poles.

A six-by-six Symbolica linear solve matches these generated poles to the
counterterm operators, with $Z_j=1+a_4 delta Z_j$. The structures follow from
$psi_0=sqrt(Z_psi) psi$, $m_0=Z_m m$, $A_0=sqrt(Z_A) A$,
$xi_0=Z_xi xi$ and $e_0=Z_e e$. The electron counterterm is proportional to
$delta Z_psi slash(p)-m(delta Z_psi+delta Z_m)$, while the vertex coefficient is
$delta Z_psi+delta Z_e+delta Z_A/2$. The photon longitudinal equation is
$xi B+(xi-1)delta Z_A+delta Z_xi=0$, where $B$ multiplies $p^mu p^nu$ in the
loop pole divided by $i a_4$. Solving with symbolic $xi$ gives a regular
continuation to Landau gauge, $xi=0$.

The comparison values are
$ delta Z_psi=-frac(xi,epsilon), quad delta Z_m=-frac(3,epsilon), quad
  delta Z_A=delta Z_xi=-frac(4 N_f,3 epsilon), quad
  delta Z_e=frac(2 N_f,3 epsilon), quad delta Z_(A m)=0. $
They obey the Ward relation $delta Z_e+delta Z_A/2=0$, or
$delta Z_1=delta Z_psi$ for the vertex renormalization constant
$Z_1=Z_psi Z_e sqrt(Z_A)$. Here $delta Z_(A m)$ multiplies the additive
auxiliary operator $a_4 M A_mu A^mu/2$; its zero follows from retaining the
mass-correction terms above. These counterterm structures are supplied from
the Lagrangian and their coefficients are solved, rather than generated from
a counterterm model. Automatic counterterm insertions and subtraction forests,
finite off-shell form factors and higher-loop renormalization remain separate
work.

The installed regression passes the four UV-structure comparisons and all six
counterterm references with symbolic $xi$ in the rebuilt HEP host. It also
checks the exact linear-system residual, Ward relations, cancellation of
physical and auxiliary mass dependence, and absence of double poles. An
independent electron self-energy calculation verifies that the longitudinal
photon denominator cancels after trace and vacuum projection, leaving no
hidden loop-momentum dependence in the scalar-family coefficients. All 78
generator tests and ten installed physics regressions pass. The notebook
passes strict Marimo checks and headless export; the existing live instance
executes its default and nine combinations of gauge parameter and flavor count
without cell errors. The native extension is installed in that same host.

== Generated one-loop QCD renormalization

`hep/qcd_renormalization.py` and
`installed_feyncalc_qcd_renormalization.py` assemble the quark, gluon and ghost
self-energies and the quark-gluon vertex for the
#link("https://feyncalc.github.io/FeynCalcExamples/QCD/OneLoop/Renormalization")[massive one-loop QCD renormalization reference].
The gluon two-point function includes gluon, ghost, quark and four-gluon
tadpole diagrams. The vertex includes both its Abelian and non-Abelian
contributions, normalized against the generated tree vertex. The internal
gluon propagator keeps symbolic $xi$, and the quark mass remains symbolic.
Shared particle color tensors close external color slots; Idenso's color
simplification and representation-aware Casimir conversion retain $C_F$ and
$C_A$. Quark-loop multiplicity is $N_f$, with $T_F=1/2$.

The calculation reuses the graph UV expansion, vacuum tensor reduction and
native IBP path described above. Four vacuum-integral powers reduce to the
single tadpole with supplied pole $A_0(M)=M/epsilon+O(1)$, checked with OneLOop.
The total gluon pole is transverse. In units of $a_4=g_s^2/(16 pi^2)$ relative
to the tree vertex, the two vertex poles are
$ frac(xi(C_F-C_A/2),epsilon) quad "and" quad
  frac(3C_A(xi+1),4epsilon). $
Their sum fixes the coupling counterterm together with the quark and gluon
field counterterms.

An eight-by-eight Symbolica solve matches the generated local structures to
Lagrangian counterterm operators, using $Z_j=1+a_4 delta Z_j$. The six
physical renormalization constants are compared with
$ delta Z_q=-frac(C_F xi,epsilon), quad
  delta Z_m=-frac(3C_F,epsilon), $
$ delta Z_A=delta Z_xi=frac(C_A(13-3xi)-4N_f,6epsilon), quad
  delta Z_c=frac(C_A(3-xi),4epsilon), $
$ delta Z_g=-frac(11C_A-2N_f,6epsilon). $
Here $Z_c$ renormalizes the ghost field and $Z_g$ the strong coupling. The
quark-gluon counterterm coefficient is $delta Z_q+delta Z_g+delta Z_A/2$;
the mass and coupling results are independent of $xi$.

The two auxiliary mass coefficients are $delta Z_(A m)=delta Z_(c m)=0$ in
this complete UV expansion. Terms correcting the auxiliary mass are retained
through logarithmic order; even the expanded massless tadpole has cancelling
UV poles. FeynCalc's nonzero auxiliary gluon mass counterterm belongs to its
selective infrared rearrangement and is a different prescription. The physical
renormalization constants agree between the two prescriptions. Analytic master
values and Lagrangian counterterm structures remain explicit inputs; this
workflow does not generate counterterm diagrams or subtraction forests.

The installed regression passes the symbolic-gauge pole comparisons, all eight
counterterm equations and the exact linear-system residual. It checks that
scalar-family coefficients contain no hidden loop momentum, the projected UV
coefficients contain neither physical nor auxiliary mass dependence, and no
double poles remain. Ruff, strict Marimo validation and a headless HTML export
pass. The notebook also passes in the existing live instance for its default
state and 27 combinations of gauge parameter, quark-flavor count and SU(N)
color group. Its interactive counterterm table checks the vertex cancellation
and the gauge-independent one-loop coefficient $beta_0=(11C_A-2N_f)/3$.

== Generated ghost-gluon vertex and crossing

`hep/qcd_ghost_vertex.py` and `installed_feyncalc_qcd_ghost_vertex.py` generate
both one-loop vertex topologies and the ghost self-energy for incoming ghosts
and antighosts. They reproduce the UV result of
#link("https://feyncalc.github.io/FeynCalcExamples/QCD/OneLoop/GhGl-Gh")[the separate ghost-gluon vertex example]
with symbolic $xi$. Every external color and Lorentz index remains open;
the full pole tensor is compared with its generated tree tensor before taking
a ratio. The two diagrams contribute, in units of $a_4$ relative to the tree,
$ frac(C_A xi,8epsilon) quad "and" quad frac(3C_A xi,8epsilon). $
Their sum is $C_A xi/(2epsilon)$ and vanishes in Landau gauge.

The shared SM model supplies the canonical UFO `UUV1` rule
`P(3,2) + P(3,3)`, equal to $-p_1$ when all vertex momenta are incoming.
Both FeynKit and GammaLoop consume that rule; the notebook only specializes
the gluon propagator's gauge parameter. Shared graph finalization remaps
indices inside momentum arguments before transporting their momentum carriers.
Each signed carrier is transported once, including when a sign change causes
a product to collapse. This preserves linearity of sums for incoming neutral
vectors as well as the charged-vector and ghost crossings.

The model correction was also validated in all three crossings of the twelve
electroweak and one QCD SM interactions using `UUV1`: 36 electroweak and three
QCD cases. Those checks compared the generated tensor, particle order and
couplings with the original UFO convention, retained the QCD
structure-constant sign, and passed 15 additional Lorentz-linearity probes.
They used the portable SM fixture without an external MadGraph installation.

The vertex workflow reuses graph UV expansion, Spenso/Idenso color algebra,
vacuum tensor reduction and native IBP. The two integral targets $I(2)$ and
$I(3)$ reduce to $I(1)$; the tadpole pole $I(1)=M/epsilon+O(1)$ is an explicit
analytic input with $M=m_("UV")^2$ and loop measure $i/(16pi^2)$.
The regression checks the full crossed tensors, the absence of residual loop
momentum and double poles, and cancellation of auxiliary mass dependence.
The computed self-energy gives $delta Z_c=C_A(3-xi)/(4epsilon)$ with no
auxiliary ghost mass pole.

The gluon field counterterm $delta Z_A$ is supplied from the separate QCD
renormalization calculation above. Combining it with the calculated ghost
and vertex poles yields
$ delta Z_g=-frac(C_A xi,2epsilon)-delta Z_c-frac(delta Z_A,2)
  =-frac(11C_A-2N_f,6epsilon). $
The notebook displays both crossing choices, the IBP reductions and the
counterterm cancellation as the gauge parameter, color count and flavor
count vary. The live notebook passes its default state and 18 combinations
of crossing, gauge, color and flavor. This checks the coupling result
independently of the quark-gluon
vertex. Finite vertex form factors, other QCD vertices and automatic
counterterm-diagram or subtraction-forest generation remain separate coverage.

== Native IBP reduction through RustRed

The Symbolica Community host now registers `rustred-feynkit` alongside FeynKit
and OneLOop in `symbolica.community.hep`. `hep.IBPFamily(family)` reads the
existing `IntegralFamily` boundary: denominator order and signs, masses,
dimension and external Gram products are retained. The bridge passes native
Symbolica expressions into RustRed in the same extension and kernel. Family
construction, momentum mappings and tensor projection stay with their existing
FeynKit owners; elimination and recurrence discovery stay with RustRed.

Use `ibp_identities()` to inspect the ordinary integration-by-parts equations,
`reduce_laporta(targets, max_depth=...)` for a bounded target search, and
`solve_parametric(sector, ...)` for recurrence rules. The returned `residuals`
are integrals unresolved at that search depth, not a proof of a minimal master
basis. Parametric rules retain nonzero conditions and exceptional index loci.
Neither reduction method supplies the analytic values of the residual integrals.

Pass `integral=I` to the existing solution or rule method to obtain a native
Symbolica sum of coefficients times `I(*powers)`. Omitting it retains the
`(powers, coefficient)` list. For the two-denominator bubble family below:

// docs-example: compile
```python
from symbolica import S

I = S("I")
laporta = bubble_ibp.reduce_laporta([[2, 1]], max_depth=2)
reduced_expression = laporta.reduce([2, 1], integral=I)
recurrence = bubble_ibp.solve_parametric([True, True], fixed=[None, 1], max_depth=1)
one_step = recurrence.rules[0].apply([2, 1], integral=I)
```

This option changes the return form without bypassing sector, power or
exception checks. Conditions still symbolic in kinematic parameters must remain
nonzero when evaluated. Laporta rules already include back-substitution through
solved targets; a parametric solution or rule performs one recurrence step per
call. Returning an expression does not recursively reduce the remaining
integrals or certify them as masters.

`hep/ibp_phi4.py` and `installed_feyncalc_ibp_phi4.py` construct the two scalar
integrands from the #link("https://feyncalc.github.io/FeynCalcExamples/Phi4/TwoLoops/Renormalization-SS")[two-loop scalar self-energy reference].
A Taylor expansion through the external momentum squared uses the shared vacuum
tensor reducer. The example computes Laporta reductions of the resulting
massive vacuum family, identifies equivalent residuals through six verified
loop-momentum mappings, and obtains a doubled-tadpole recurrence parametrically.
Analytic tadpole and equal-mass vacuum Laurent coefficients are supplied as
reference inputs. Together with the stated loop measure and counterterms,
the notebook checks the UV poles and renormalization constants; the IBP solver
does not evaluate those analytic integrals itself.

`hep/ibp_bubble.py` constructs a bubble with unequal nonzero masses and nonzero
external momentum. The targets $I_(2 1)$, $I_(1 2)$ and $I_(2 2)$ reduce to the
two tadpoles and the scalar bubble $B_(1 1)$. Expanding the dimension-dependent
coefficients before removing dimensional regularization retains finite terms
from the individual residual integrals' UV poles. OneLOop supplies their scalar
values; independent Feynman-parameter quadrature below threshold checks the
finite raised-power results. Generic kinematics avoid exceptional coefficient
denominators; this example does not establish threshold continuation or
exceptional-mass and Gram limits.

`installed_feyncalc_ibp_phi4.py` and `installed_feyncalc_ibp_bubble.py` pass in
the rebuilt, installed native host. All 36 bridge tests and 76 generator tests
also pass. Both IBP notebooks execute without cell errors in the existing
Marimo instance; the bubble checks all twelve combinations of its three raised
integrals and four kinematic presets. Headless exports pass for both IBP
notebooks and the Bhabha/Møller notebook. These checks validate the stated
reductions and observables, not a general master-basis certification or the
remaining gallery's IBP coverage.

== Generated electron anomalous magnetic moment

`hep/gminus2.py` and `installed_feyncalc_gminus2.py` generate the tree and
one-loop electron-photon vertices from the same model. They reproduce the
#link("https://feyncalc.github.io/FeynCalcExamples/QED/OneLoop/El-GaEl")[gallery's one-loop anomalous magnetic moment]
and retain a generic spacelike photon momentum. External electrons satisfy
$p^2=p'^2=m^2$ and $t=(p-p')^2<0$. Shared particle spin sums close the two
projector traces, and Symbolica solves their two-by-two Gram system for the
coefficients of $gamma^mu$ and $(p+p')^mu/(2m)$. Normalizing against the generated
tree vertex retains every graph weight and fixes the external-fermion phase.

The diagram's integral family rewrites the projected numerator. Native Laporta
reduction resolves nine targets to two shifted tadpoles and an equal-mass
bubble. The tadpole identification uses their shift equivalence; residuals at
this bounded search depth are not a general master-basis certification.
With the common $e^2/(16 pi^2)$ removed, the coefficient of the second basis
vector is
$ b(D) = frac(2(D-5)(-(D-2)A_0+2m^2(D-3)B_0), (D-3)(t-4m^2)). $
Expanding at $D=4-2 epsilon$ cancels its UV pole and retains the finite term
$ b = frac(4(A_0^"fin"+m^2-m^2 B_0^"fin"), t-4m^2). $
The Gordon decomposition gives $F_2=-alpha b/(4 pi)$, so taking the limit only
after integration yields $F_2(0)=alpha/(2 pi)$ exactly, independently of the
mass and renormalization scale. Taking $t=0$ before solving the Gram system
would make those projectors degenerate.

Five spacelike points, including unequal mass and scale, agree with independent
96-node Feynman-parameter quadrature and OneLOop's bubble derivative to better
than $2 times 10^(-12)$. The live notebook varies $-t/m^2$ over five decades;
its default and five slider states execute without cell errors. This covers
the Pauli form factor, not the complete renormalized Dirac form factor or
analytic continuation through timelike thresholds.

== Generated two-loop massless electron self-energy

`hep/electron_two_loop.py` and `installed_feyncalc_electron_two_loop.py` generate
all three one-particle-irreducible two-loop QED electron diagrams. They validate
the Feynman-gauge specialization $xi=1$ of the
#link("https://feyncalc.github.io/FeynCalcExamples/QED/TwoLoops/Renormalization-LeAle-Massless")[massless electron self-energy reference].
The graph's shared auxiliary-mass UV expansion, Idenso traces and vacuum tensor
reduction produce 22 distinct powers in a common massive vacuum family.
Verified momentum maps identify denominator orderings; native Laporta reduction
leaves the equal-mass sunset and three equivalent products of tadpoles.

The amputated two-point kernel removes only the generator's named external
Wick-order sign. Internal fermion-loop signs and every other graph weight remain
intact. The installed regression independently checks this kernel convention
with a one-loop projector. The single closed fermion loop receives a symbolic
flavor multiplicity $N_f$.

Writing $a_4=e^2/(16 pi^2)$ and $M=m_"UV"^2$, the coefficient of
$i a_4^2 slash(p)$ in the bare UV poles is
$ frac(1,2 epsilon^2)
  + frac(log(4 pi)-log(M)-17/12-7N_f/3,epsilon). $
Analytic vacuum Laurent coefficients and the summed one-loop counterterm
insertions are explicitly supplied reference inputs, including the auxiliary
photon-mass counterterm. They yield
$ Z_psi = 1-frac(a_4,epsilon)
  + a_4^2 (frac(1,2 epsilon^2)+frac(4N_f+3,4 epsilon)), $
with the auxiliary mass and loop-measure logarithms cancelling exactly.
The notebook and installed-host regression pass; automatic subtraction forests,
arbitrary gauge parameter and the finite off-shell self-energy remain separate
work. The displayed reference counterterms are not inferred by the IBP solver.

== Generated two-loop photon renormalization

`hep/photon_two_loop.py` and `installed_feyncalc_photon_two_loop.py` reproduce
#link("https://feyncalc.github.io/FeynCalcExamples/QED/TwoLoops/Renormalization-GaGa")[the massless two-loop photon example]
with symbolic gauge parameter $xi$. All three native diagrams retain their
closed-fermion-loop signs and both open Lorentz indices. Native IBP reduces
112 vacuum targets to the sunset and equivalent tadpole products.

The shared `uv_expansion` and `uv_counterterm` accept signed `edge_powers`,
using the denominator API's edge IDs and selection semantics. The Feynman
photon term uses power one; the longitudinal term uses a polynomial numerator
and power two. Freezing the massive denominators before retaining the explicit
auxiliary-mass zeroth coefficient implements the reference's selective infrared
rearrangement. The default UV expansion still retains its mass compensation terms.

Four local insertions on a generated one-loop bubble calculate the counterterm
sum: two vertex factors and two fermion kinetic insertions with squared
propagators. The supplied one-loop input is $delta Z_psi=delta Z_1=-xi/epsilon$.
Five one-loop targets reduce to the tadpole. Its finite term and the analytic
two-loop master poles are explicit inputs; coefficients are checked to be
regular at $D=4$ before using these truncated series. Tensor basis elements
remain fixed during the Laurent expansion.

With $Z=1+a_4 delta Z_(1)+a_4^2 delta Z_(2)$, the calculated bubble and insertions
give
$ delta Z_(A,1)=-frac(4N_f,3epsilon), quad delta Z_(A,2)=-frac(2N_f,epsilon), $
$ delta Z_("Am",1)=-frac(2N_f,epsilon), quad
  delta Z_("Am",2)=frac(N_f xi,2epsilon)-frac(N_f(2N_f+xi),epsilon^2). $
The auxiliary operator is $i M(Z_("Am")^2-1)g^(mu nu)$, including the square of the
one-loop coefficient at second order. Both tensor coefficients cancel exactly;
physical field renormalization is gauge independent and all logarithms cancel.
These explicit insertions do not provide arbitrary counterterm models or forest
enumeration. Finite two-loop amplitudes and general analytic master evaluation
remain separate work. The notebook passes its default state and 18 live
gauge/flavor/auxiliary-mass combinations; the installed regression and
headless export pass as well.

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
open endpoints, as in the
#link("https://feyncalc.github.io/FeynCalcBook/ComplexConjugate.html")[FeynCalc color-trace conjugation example]. `installed_color_conjugation.py` checks SU(2), SU(3) and SU(5),
including the summed norm $(N_c^2-1)^2/(4 N_c)$ and conjugation twice. The live
`hep/color_algebra.py` notebook demonstrates the existing tensor constructors,
chain collection, conjugation and color contraction. Symmetric, antisymmetric
and cyclic generator groups, including nested groups and complex coefficients,
conjugate without expanding their permutations. Tests compare against Spenso's
explicit projector expansion and verify reversal parity for three-generator
traces. Scalar powers, inverse factors and scalar functions inside compact words
use Symbolica's existing conjugation semantics, including symbolic exponents.
Unknown matrix factors remain explicitly conjugated. Public Python regressions
construct weighted words with the existing `chain` and `trace` APIs, check
involution and Dirac adjoints, and independently evaluate a weighted SU(3) norm
with explicit matrices. No separate FeynKit or GammaLoop conjugation
implementation is introduced.

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
Numerical validation also checks every free-index component of the mixed
trace/chain identity. The identity uses the existing typed metric factory as
`TensorExpression.g(fund, fund.dual())("i", "j")`; the optional second
representation preserves its supplied logical port order and validates exact
dimension equality. The shared Spenso syntax classifier recognizes dual
representation slots inside metrics, so a fundamental identity retains its
oriented tensor ports when added to generator products. Regression tests cover
all parser filters, expanded and opaque parsing, and scalar precontraction;
a scalar cannot be added to an open identity tensor.
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

== Generated antisymmetric color vertices

The shared UFO lowering accepts `Epsilon(i,j,k)` for three triplet slots and
`EpsilonBar(i,j,k)` for three antitriplet slots. Both use Idenso's existing
antisymmetric epsilon symbol and determinant expansion. Color conjugation
exchanges the two representations without reversing their index order; the
Lorentz epsilon convention is unchanged. Pair reduction requires dual spaces
with matching dimensions. Two same-orientation color epsilons, mismatched
spaces and rank-three tensors in a different fundamental dimension do not
acquire a determinant identity.

Summed fundamental indices may join dual slots while external slots retain the
orientation fixed by their particle records. Color identity chains transmit
this orientation before lowering; this also permits contracted generator
products such as `T(a,-1,-2)*T(b,-2,-1)`.

`installed_ufo_color_epsilon.py` generates a conjugate pair of cubic interactions
between three distinct complex scalar species. They are a toy color model, not
a Standard Model interaction. The regression checks the three open color ports,
conjugation twice and the positive norm $epsilon_(i j k) epsilon^(i j k)=6$.
Shared Rust tests also verify $epsilon_(i j k) epsilon^(i j l)=2 delta_k^l$,
permutation signs, incompatible representations and dummy-index contractions.
`hep/color_epsilon.py` shows the generated diagrams, both tensors and the
incoming-triplet color average $6/3=2$ in the same Marimo instance.
This establishes symbolic epsilon interaction support; numerical tensor-library
components and sextet interaction tensors remain outstanding.

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
momentum declarations. These validate the algebraic family boundary separately
from the RustRed reduction workflows described above.
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

Charge conjugation now has a shared algebra and tensor-data implementation.
Idenso's `GammaLibrary.charge_conjugation` registers
`spenso::charge_conjugation` with the ALOHA Weyl convention
$C = -i γ^2 γ^0$. It is real and antisymmetric, with $C^2 = -1$.
Explicit-index products and compact chains use the same Idenso simplifier;
charge-conjugation sandwiches transpose each supported matrix with its correct
sign, preserving factor order. Four-dimensional gamma and slash arguments are
supported, including `P(label,mink(4))`; symbolic-D gamma arguments remain
opaque. Self-dual common-end chain joins transpose one ordered word without
conjugating scalar coefficients.

`spenso-hep-lib` owns the single sparse component definition used by the numeric
and symbolic HEP libraries and generator grouping. Both the shared generator
and GammaLoop's UFO reindexer lower `C(i,j)` to this registered tensor. The
focused symbolic and common-end component tests pass, as do the shared Cargo
check and Clippy checks. `installed_charge_conjugation.py` adds explicit-index,
compact-slash and independent Weyl-component comparisons. Its initial symbolic
identities pass in the rebuilt public host; full installed-host validation
awaits fixture setup corrections. This primitive support does not close the
Majorana-generation gap: general Majorana and fermion-number-violating diagrams
still require changes to external fermion-flow normalization, which currently
accepts particle/antiparticle pairs.

An example counts as reproduced only after exercising the relevant shared
components and comparing its final observable or symbolic identity with the
published result. Remaining work includes tree-level QED/QCD/EW observables,
polarized amplitudes, anomaly conventions, loop-renormalization examples,
broader IBP reductions and analytic master evaluation, and the two-loop
topology-minimization example.
Full FeynCalc gallery parity is not yet achieved.
]
