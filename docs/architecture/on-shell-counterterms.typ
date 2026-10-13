= On-shell counterterms for massive fermions
<on-shell-counterterms>

#let slash(body) = math.cancel(body, stroke: 0.55pt)

This maintained design specifies the massive-fermion component operation and its
independent validation oracles. The target scheme has on-shell masses and
MS-bar fields and couplings. A support claim requires the separate algebraic,
local-limit and integrated checks below; a successful numerical comparison does
not establish any of the other claims.

The prescription is based on
#link("https://arxiv.org/pdf/2203.11038v1")[Capatti, Hirschi and Ruijl,
arXiv:2203.11038v1], sections 4.3–4.4 and 5.1.3–5.2.3. Equations 4.22–4.24
specify the local operator, equations 5.10–5.11 identify its mass correction,
and figures 10–11 motivate the nested tests. The calculations below give an
independent one-loop oracle in the repository conventions.

The comparison requirements below are maintained test contracts. Historical
pass records and measurements describe the original OS stack before its
rebase onto FeynKit. They do not certify the rebased implementation. Current
validation records identify the tested revision, selection and run; exact
reconstruction, local limits and integrated renormalization remain separate
claims.

== Ordered component operation

Use the Minkowski metric $(+, -, -, -)$, $d=4-2 epsilon$ and
$slash(q)=gamma^mu q_mu$. A fermion self-energy is an ordered linear map from
its incoming spinor slot to its outgoing slot. Define

$ p_("OS") = (m,0,0,0), quad P_sigma = (1 + sigma gamma^0)/2, $
$ O[Sigma] = sum_(sigma=plus.minus 1) P_sigma Sigma(sigma p_("OS")), $
$ T_("OS")[Sigma] = U_1[Sigma] + O[Sigma] - U_1[O[Sigma]]. $

$U_1$ is the existing ordinary degree-one UV operation, including its auxiliary
vacuum mass prescription. First complete both on-shell branches, including the
projectors. Then apply $U_1$ to those complete expressions. Matrix products are
ordered: the left projector contracts the outgoing open spinor index with a
fresh internal spinor index. A commuting scalar product of gamma symbols is
not an implementation of this operation.

For a completed child expression $X$, the recursive operation creates
$C_gamma[X]=-T_("OS"),gamma[X]$. This is the existing single subtraction sign.
The relative minus sign in $-U_1 O$ is internal to the operator and does not
supply the forest sign. Thus one nested family is
$C_Gamma[C_gamma[X]]=+T_Gamma[T_gamma[X]]$, with both complete operations.
No sign is reconstructed from forest depth or multiplied at the terminal node.

The substitution belongs only to the active component and replaces its entire
incoming boundary momentum $P_gamma mapsto sigma p_("OS")$. If the ambient
routing writes $P_gamma=r+p$, both pieces belong to that boundary value: the
operation does not replace only $p$ and retain $r$. An internal component edge
carrying $k+P_gamma$ becomes $k+sigma p_("OS")$, with the component loop $k$
fixed. The cograph's separate occurrence of the boundary momentum remains the
live $r+p$, and other components retain their own momenta. A global replacement
of an ambient symbol would change the surrounding graph and is invalid. The
component boundary index, directed fermion flow and physical mass identify the
operation; a generic two-point scalar graph is insufficient.

== Scaling and ownership

The rest-frame point is the expression $m(1,0,0,0)$, never an unrelated frozen
four-vector. In a subsequent UV Taylor jet, physical masses and external
momenta have weight one, so $m mapsto t m$ implies
$p_("OS") mapsto t p_("OS")$. The completed $O$ branches consequently contribute
their first-order mass and energy-shift coefficients to $U_1 O$. In the
implementation's inverse-hard convention this is expressed through its existing
hard-loop rescaling and physical-mass bookkeeping. The vacuum mass keeps its
existing hard weight; it must not be identified with the physical mass.

The implementation distinguishes the expansion bookkeeping mass `mUVexp` from
the terminal auxiliary vacuum mass `mUV` during recursion and outer Taylor
operations. The completed physical integrand identifies `mUVexp = mUV` and sets
`epsilon = 0` at its final normalization boundary. Local UV and soft-limit
certificates use that same completed physical copy; raw operator and source
reconstruction identities retain the independent bookkeeping symbols. This
identification never sets the auxiliary mass equal to the physical fermion mass.

A surrounding soft operation scales its own boundary momenta, leaving the
physical mass and an already constructed child rest-frame point fixed. A
surrounding OS operation evaluates only its own boundary momentum; inner
rest-frame points remain related to their own physical masses. Distinct
components using the same mass symbol retain distinct owners.

Completed child loops are fixed generators in a compatible sub-LMB. An outer
operation acts on the complete child expression with that basis; it does not
shift a child loop coordinate by an external momentum to fit another chart.
The existing hard, physical-source and soft provenance rules remain in force.
Integrated coefficients retain their source owners and consumed-loop weights;
a containing operation may act on their mass dependence. Only normalized
localization kernels are frozen outside subsequent operations. In particular,
freezing an entire integrated coefficient would discard needed outer terms.

== Independent one-loop local oracle

Strip the common QCD color, coupling and Feynman-rule phase. Route the gluon
momentum as $k$ and the fermion momentum as $k+p$. In Feynman gauge the open-chain
identity $gamma_mu gamma_nu gamma^mu=(2-d)gamma_nu$ gives

$ F(p,k) = frac((2-d)(slash(k)+slash(p))+d m,
  k^2 ((k+p)^2-m^2)), quad D=k^2-M^2. $

Every denominator carries the same causal $+i 0$ prescription; it is suppressed
in these algebraic formulas. The source numerator is one graph-local factor.
The finite expansion below is an oracle, not permission to distribute an
arbitrary production numerator across a forest sum.

Taylor expansion through degree one gives

$ U_1 F = frac((2-d)(slash(k)+slash(p))+d m,D^2)
  - frac(2(2-d)slash(k)(k dot.op p),D^3), $
$ O F = sum_(sigma=plus.minus 1)
  frac(P_sigma ((2-d)(slash(k)+sigma m gamma^0)+d m),
    k^2(k^2+2 sigma m k^0)), $
$ U_1 O F = frac((2-d)slash(k)+2m,D^2)
  -frac(2(2-d)m gamma^0 slash(k) k^0,D^3). $

In particular, $sum_sigma P_sigma=1$ and
$sum_sigma sigma P_sigma=gamma^0$. These identities are applied after the
projector has multiplied the numerator on the left. Define the vacuum remainder

$ V = U_1 F-U_1 O F
  = (2-d) (frac(slash(p)-m,D^2)
  -frac(2(slash(k)(k dot.op p)-m gamma^0 slash(k) k^0),D^3)). $

Then $T_("OS")F=O F+V$. Direct differentiation of the source gives

$ V=(p-gamma^0 p_("OS"))^mu
  [partial_(p^mu) F]_(p=0,m=0,"vacuum mass"=M). $

This proves the composed formula agrees with equation 4.24 by exact ordered
algebra. Moving $gamma^0$ through $slash(k)$ changes the spatial gamma terms
and fails this identity. At $m=0$, both branches evaluate the same $F(0,k)$,
the projectors sum to the identity, and the operator becomes
$U_1+S_0-U_1 S_0$, the existing degree-one soft-refined operation. This is an
algebraic specialization; a massive on-shell integral cannot be made massless
by naively setting $m=0$ inside its logarithms.

The checked-in validation graph chooses the quark momentum as $k$ and the
massless-gluon momentum as $k-p$. For that original homogeneous basis, write
$F'=N(k,m)/((k^2-m^2)(k-p)^2)$, with $N_0=N(k,0)$. The independent compact
remainder used by the production test is instead

$ V'=frac(2(N_0(k dot.op p)-m gamma^0 N_0 k^0),D^3). $

The two OS branches have denominators
$(k^2-m^2)(k^2-2 sigma m k^0+m^2)$. This is a separately derived expression
in the fixture basis, not a post-subtraction loop shift of the previous oracle.
The tests compare the source QCD numerator as a finite open-spinor tensor,
retain identical denominator products separately and check that attaching the
projector to the opposite spinor port changes the answer.

For the UV check, compare the exact source and completed subtraction under
$k mapsto Lambda k$. $O-U_1 O$ has no terms at the retained divergent orders,
so $T_("OS")$ has the same UV jet as $U_1$. Check this statement before energy
integration and after summing CFF orientations. Infrared claims are separate:
individual on-shell branches possess their physical soft denominators, and
need not be pointwise soft finite on their own.

== Integrated policy and normalization

Let $[X]$ mean integration with the repository Minkowski measure. For one loop
its MS-bar continuation is

$ [X]= (mu^2 exp(gamma_("E"))/(4pi))^epsilon
  integral frac(d^d k,(2pi)^d) X. $

Thus
$ I_2(M)=[D^(-2)]=i/(16pi^2) exp(gamma_("E")epsilon)
  Gamma(epsilon)(mu^2/M^2)^epsilon $.
This is Vakint's `MSbar` normalization with `additional_normalization = "1"`.
The CFF energy convention is $d k^0/(2pi i)$, leaving $i/(2pi)^3$ per
three-dimensional loop. Forest signs, the Feynman-rule phase and the contour
orientation are separate factors; no extra minus sign per integrated loop is
introduced.

Analytically integrate only the vacuum remainder $V$, not the full $O+V$.
Writing $"Pole"$ for negative Laurent powers and $"Fin"=1-"Pole"$ for all
retained nonnegative powers, the signed counterterm and signed addback are

$ C_("local")=-(O+V), quad
 C_("scheme")=+N "Fin"[V], quad [N]=1. $

Equivalently, the generic signed-series machinery receives $-V$, so that its
negative nonnegative projection gives $+"Fin"[V]$. Its integrated series is a
scheme payload, not the integral of the complete local counterterm. The result
at the integrated boundary is

$ [C_("local")+C_("scheme")]=-[O]-"Pole"[V]. $

The standalone `renormalize` output uses the terminal analytic projection,
matching the existing MUV command semantics. For OS it therefore exports the
finite vacuum scheme correction $+"Fin"[V]$; it is not an export of the full
OS mass counterterm. The physical mass term remains in the local $-O$ branch
and enters the assembled integrand through that branch.

The retained $-[O]$ supplies the OS mass subtraction. Passing the complete
$-(O+V)$ to the generic finite-addback policy would cancel its finite mass
correction and implement a different mass scheme. This distinction must be
explicit in recursive state and integration dispatch. The localized coefficient
has the same owner and epsilon-depth rules as ordinary integrated components.

Equations 5.10–5.11 determine this policy without an ambiguous reuse of the
paper's superscript $T^"os"$: applying an idempotent $U_1$ to the complete
$T_("OS")$ yields $U_1$, so $U_1-U_1 T_("OS")$ would incorrectly vanish.
The vacuum payload is $U_1-U_1 O$.

Tensor integration provides an independent simplification:

$ [k^mu k^nu D^(-3)] = g^(mu nu)/d [k^2 D^(-3)], quad
  [k^2 D^(-3)] = d/4 I_2(M), $
$ [V] = (2-d)/2 I_2(M)(slash(p)-m). $

Keep $d$ symbolic through both gamma algebra and tensor reduction. The external
spinor trace normalization remains four. The identity $gamma^0 gamma^0=1$
then gives the scalar mass term. Setting $d=4$ before multiplying the pole loses
a finite rational contribution.

For comparison with conventional self-energy coefficients, remove the common
$i$ from the loop measure and include $g_s^2 C_F$. Set
$a=g_s^2/(16pi^2)=alpha_s/(4pi)$ and $L_M=log(mu^2/M^2)$. The vacuum answer is

$ [V]_(Sigma) = -a C_F (slash(p)-m)
  (1/epsilon+L_M-1
  +epsilon (L_M^2/2-L_M+pi^2/12) \
  +epsilon^2 (L_M^3/6-L_M^2/2+(L_M-1)pi^2/12-zeta(3)/3)
  +O(epsilon^3)). $

The finite $-1$ inside this bracket comes from the dimensional numerator. One
enclosing $1/epsilon$ pole needs the child's $epsilon$ coefficient; two enclosing
poles need its $epsilon^2$ coefficient. At $mu=M$, the signed finite scheme
coefficient, after removing its common loop normalization, coupling and
$(slash(p)-m)$ matrix, is

$ 1-pi^2/12 epsilon+(pi^2/12+zeta(3)/3)epsilon^2+O(epsilon^3). $

Its nonzero quadratic coefficient contributes to the finite part of a
three-loop product and is lost if the child is truncated at $epsilon^1$.
The existing whole-forest loop bound therefore requests four Vakint terms for
a three-loop problem: a one-loop child retains powers $-1,0,1,2$, while a
three-loop master retains powers $-3,-2,-1,0$. This is the required depth, not
an additional expansion margin. Do not truncate every completed child to
epsilon zero. Backend insufficiency must produce an error. The higher orders
do not change the scheme: the analytic addback still contains only
$+"Fin"[V]$, and the physical OS mass branch $-O$ remains local.

== Independent integrated mass and field references

A direct Feynman-parameter calculation avoids using the production UV
integration path as its own oracle. Combine the two source denominators with
$x in [0,1]$ and shift $ell=k+x p$. At $p^2=m^2$ the combined mass is $x^2m^2$.
The projected numerator is $m(2+(d-2)x)$, so

$ [O]_(Sigma)=a C_F m exp(gamma_("E")epsilon)Gamma(epsilon)
 (mu^2/m^2)^epsilon integral_0^1 (2+(2-2epsilon)x)x^(-2epsilon) d x $
$ =a C_F m exp(gamma_("E")epsilon)Gamma(epsilon)
 (mu^2/m^2)^epsilon frac(3-2epsilon,1-2epsilon). $

With $L_m=log(mu^2/m^2)$ this is

$ delta m_("loop")=a C_F m (3/epsilon+3L_m+4
 +epsilon(3L_m^2/2+4L_m+8+pi^2/4) \
 +epsilon^2(L_m^3/2+2L_m^2+8L_m+16+L_m pi^2/4+pi^2/3-zeta(3))
 +O(epsilon^3)). $

The coefficient check integrates the parameter weight independently, using

$ integral_0^1 x^r log(x)^n d x = frac((-1)^n n!, (r+1)^(n+1)). $

The moments at $r=0,1$ give $3+4epsilon+8epsilon^2+16epsilon^3$.
The cubic weight is necessary because it multiplies the scalar master's
$1/epsilon$ pole. Combining it with
$log(exp(gamma_("E")epsilon)Gamma(1+epsilon))
=pi^2 epsilon^2/12-zeta(3)epsilon^3/3+O(epsilon^4)$
gives the displayed quadratic mass coefficient without expanding a production
Taylor result.

The extended analytic tests request exactly three or four terms, compare the
full bare-fermion normalization in both tensor-integration modes at two
auxiliary scales, and select AlphaLoop and MATAD separately with symbolic
$zeta(3)$. The one- and two-enclosing-pole products include negative controls
that discard the required positive-epsilon child term. Both backends pass
through $epsilon^2$. MATAD's symbolic mode retains $Gamma(1)$ and
polygamma constants; the comparison boundary uses the exact identities
$Gamma(1)=1$, $psi(1)=-gamma_("E")$, $psi^((1))(1)=pi^2/6$ and
$psi^((2))(1)=-2zeta(3)$, without numerical substitution. All three extended
reference tests pass after this exact change of constant basis.

Define the propagator inverse as $slash(p)-m-Sigma(p)$ and the bare parameter
as $m_0=Z_m m$. Then the conventional one-loop constants are

$ delta Z_m^"OS"=-a C_F(3/epsilon+3L_m+4), quad
 delta Z_2^"MS-bar"=-a C_F/epsilon. $

The integrated unsigned subtraction is
$K=delta m_("loop")+(slash(p)-m)delta Z_2^"MS-bar"$.
It reproduces the OS pole position with MS-bar field renormalization, and the
signed insertion is $-K$. This explicitly distinguishes the positive loop
mass shift from the negative bare-mass counterterm. A comparison of the scalar
coefficient alone must include the $-m delta Z_2$ piece. Field renormalization
here is not the OS residue condition; no derivative of the full on-shell
self-energy is subtracted.

The auxiliary scale $M$ drops out after the local remainder and its scheme
addback are integrated and combined. It need not cancel pointwise in an
individual local term or CFF orientation. The physical scale $m$ remains in
the OS mass relation.

== Exact CFF reconstruction and residue oracle

Retain the raw completed four-dimensional sectors with their source occurrence
owners. Reconstruct only a separate projection view, using the production
assignment plan. The following are independent exact certificates:

+ Substitute the source factorized numerator and reconstructed numerator into
  common formal loop four-momenta, retaining the fixed component shifts, and
  prove their difference zero by factor-preserving rewrites.
+ Match denominator momenta up to the permitted overall sign, masses, powers,
  full causal prescriptions and component domains. The physical occurrence
  witness still determines incidence and provenance.

The two poles of the OS branch involve $k^0+sigma m$. They are fixed affine
energy shifts, while the fermion remains a physical-mass propagator. A matching
spatial momentum alone does not certify a denominator. Neither certificate
may be replaced by evaluation at sample coordinates.

The first milestone validates the point bridge for an isolated component.
It lifts the fixed point into a temporary source-crown coordinate and requires
the inverse substitution to reproduce the exact source. This is injective only
when that same coordinate does not also represent a live surrounding momentum.
A nested numerator can contain both the child's fixed point and a parent's
live crown momentum with the same edge owner; merging them would change the
graph. Two disjoint point owners likewise require separate assignments.
The exact inverse-reconstruction guard remains mandatory.

The nested and disjoint implementation therefore uses a per-component source
assignment with distinct neutral coordinates for the live crown and the fixed
point. Other components' points must remain affine spectators until their own
contour step. Exact reconstruction and CFF projection now pass for all three
physical nested fixtures and the disjoint fixture below. Those checks certify
the source assignment and projection boundary. Independent complete-family
Taylor identities now also pass. The separate small-graph limit checks cover
all four hard subspaces below and the relevant constant/linear IR soft jets;
the OS/OS, OS-inside-IR, IR-inside-OS and disjoint radial CFF controls also
pass in the domain below. Counting points per component
alone does not solve the live-crown collision.

As a separate contour oracle, let $omega=abs(bold(k))$,
$E=sqrt(bold(k)^2+m^2)$, and
$N_sigma(z)=(2-d)(gamma^0(z+sigma m)-bold(gamma) dot.op bold(k))+d m$.
Closing the real energy contour below gives

$ integral frac(d z,2pi i) P_sigma frac(N_sigma(z),
 (z^2-omega^2+i 0)((z+sigma m)^2-E^2+i 0)) = $
$ -P_sigma (frac(N_sigma(omega),2omega((omega+sigma m)^2-E^2))
  +frac(N_sigma(E-sigma m),2E((E-sigma m)^2-omega^2))). $

Apply the residue derivative formula to the squared and cubed vacuum poles in
$V$ and combine them with both branches. Avoid coincident poles when using
pointwise samples; exact symbolic pole cancellation is a separate check. The
local signed counterterm introduces its one forest minus sign after this
linear energy integration. Compare complete orientation sums first, then
orientation-local values when the representation exposes compatible labels.

== Validation graphs and support boundaries

#import "../assets/typst/portal-graphs/layout.typ": layout as draw-graphs

// Draw the actual fixture incidence. Splitting the invisible external node is
// a presentation-only transformation; internal vertices and edge IDs survive.
#let show-fixture(path) = {
  let input = read(path)
    .replace("ext [style=invis]", "ext_in [style=invis]; ext_out [style=invis]")
    .replace("ext ->", "ext_in ->")
    .replace("-> ext [", "-> ext_out [")
  context draw-graphs(input, columns: 1, unit: 1,
    amplitude-mode: true,
    edge-style-options: (momentum-arrows: false, show-edge-index: true,
      show-momentum: false, show-particle: true),
    style-options: (node-label: none),
    additional-data: (steps: 350, seed: 42))
}

#figure(show-fixture("../../tests/resources/graphs/local_os_top_self_energy.dot"),
  caption: [One-loop validation graph from `local_os_top_self_energy.dot`.
  The open fermion line has physical mass $m$; the loop contains one quark and
  one massless gluon. The independent oracle uses the equivalent routing
  $k$ for the gluon and $k+p$ for the quark.])

#figure(show-fixture("../../tests/resources/graphs/local_nested_soft_top_self_energy.dot"),
  caption: [Existing nested validation topology from
  `local_nested_soft_top_self_energy.dot`. Internal edge 2 and child edge 3
  provide the compatible homogeneous loop basis. This is the topology needed
  for the figure-11 OS-inside-OS extension; its pre-existing soft tests do not
  establish OS support.])

#figure(show-fixture("../../tests/resources/graphs/paper_appendix_b1_nested_gluon_self_energy.dot"),
  caption: [Figure-10 topology from
  `paper_appendix_b1_nested_gluon_self_energy.dot`: OS quark child on edges
  3–4 inside an IR gluon parent on edges 2–6. The additional logarithmic
  fermion-box component on edges 2, 3, 5, 6 must also be accounted for.])

#figure(show-fixture("../../tests/resources/graphs/local_soft_gluon_in_top_self_energy.dot"),
  caption: [New physical fixture `local_soft_gluon_in_top_self_energy.dot`:
  the IR quark-loop gluon self-energy on edges 3–4 lies inside an OS quark
  self-energy. Gluon carrier edge 2 is the outer loop coordinate, keeping
  its soft point independent of the parent's external momentum.])

#figure(show-fixture("../../tests/resources/graphs/local_disjoint_os_top_self_energies.dot"),
  caption: [New physical fixture `local_disjoint_os_top_self_energies.dot`:
  two disjoint OS components on edges 2–3 and 5–6 share a physical quark
  line. Bridge edge 4 retains its live momentum while the two subtraction
  points remain local to their respective components.])

#figure(show-fixture("../../tests/resources/graphs/local_os_three_loop_top_self_energy.dot"),
  caption: [Larger all-OS fixture `local_os_three_loop_top_self_energy.dot`.
  The primitive child on edges 4–5 lies inside the proper two-loop component
  on edges 3–7, itself inside the three-loop quark self-energy. Carriers
  $p=Q_2$, $k=Q_3$, $l=Q_4$ are homogeneous. The proper component's OS point
  belongs to edge 2; its spectator cograph on edges 2, 8, 9 stays live.])

The larger example first tests the complete four-family proper two-loop
component at fixed, formal $p$: independent bounded $U_1$, $O$, $U_1 O$
certificates, the owned mass points on edges 2 and 3, and the four hard
subspaces of $k,l$. Exact cograph attachment then carries these component
identities into the three-loop graph without distributing its numerator.
This proper-component check passes (5.78 seconds). Separate tests assemble
$(1-T_G)(1-T_Gamma)(1-T_gamma)X$ as eight signed families. They certify each
outer $U_1$, $O$, $U_1 O$ branch after its inner prefix, then compose the
identities by linearity; this eight-family certificate passes (73.89 seconds).
The separate hard matrix contains the fourteen nontrivial
partitions of the carriers $p,k,l$ and a fixed-momentum block: seven ordinary
hard subsets and seven correlated assignments, including fixed $k-p$ and
$l-k$. Each independent hard block contributes four measure powers. Bare
nonzero controls accompany the superficially divergent strata. All fourteen
strata pass: the ordinary subsets take 580.22 seconds and the correlated
assignments 26.43 seconds. Their exact records are retained in
`tests/artifacts/os_three_loop_repeated_pairs_validation.json`. A separate
open-spinor logarithmic coefficient certificate passes in 241.53 seconds
(`tests/artifacts/os_logarithmic_coefficient_validation.json`). These large
checks belong to the full main-push selection, preserving the ordinary
four-minute bound. The open-spinor oracle uses repeated-pair gamma reduction
and exact finite tensor contraction; production gamma canonicalization is
unchanged. It preserves the graph numerator's factorization, including in
diagnostic copies. These local UV certificates are independent of the
eight-family operator identity and of integrated renormalization conditions.
The example uses OS for every fermion self-energy. A three-loop IR parent
with a two-loop OS child would also contain forests that skip the middle component; the existing one-loop
insertion Ward certificate does not justify those mixed forests, whose
public guard remains in force.

The first milestone is validated on the one-loop graph: the local ordered OS
operator, its vacuum-only analytic correction, exact CFF reconstruction,
independent contour checks and the integrated mass/field references above.
The local identities are checked in both fermion-loop and gluon-loop input
bases. The dimensional analytic checks retain the positive epsilon coefficients
through quadratic order and exercise both integration modes at two auxiliary
scales with both analytic backends.
The physical Dirac
CFF check uses exact rational control points in addition to the general source
reconstruction certificates. Production generation and evaluation checks also
pass for local and integrated OS configurations.

The nested and disjoint projection checks, the certified nested analytic
payload checks and the disjoint owned-product/bridge checks now pass. The
figure-10 topology also passes production generation and evaluation. Independent
complete-family Taylor comparisons now pass for all three nested fixtures.
Their small-graph UV matrix covers child-only, parent-only, simultaneous and
correlated hard momenta; the relevant constant/linear IR soft jets also pass.
Small exact child scalar-residue controls also pass at the rest-frame OS
points and at moving physical on-shell points. The complete OS-inside-OS
CFF radial matrix passes for the physical child-gluon, rest-frame child-gluon,
outer-gluon and simultaneous-gluon rays (299.55 seconds). The OS-inside-IR
matrix also passes its physical and rest-frame child-gluon rays (722.76
seconds). The retained record is `tests/artifacts/os_full_radial_validation.json`.
The IR-inside-OS radial check also passes on the faster explicit-direct route
(57.74 seconds, `tests/artifacts/os_proper_guard_validation.json`). Configuration
acceptance does not establish those claims. The public guard
currently admits OS/OS nesting, an OS parent with an IR child, and disjoint OS
components. An IR parent with an OS child additionally requires the actual
minimal-QCD rules and propagators to pass
`Graph::has_qcd_self_energy_insertion_ward_identity`. Other mixed prescriptions
remain unsupported. Both direct 3D modes support a complete isolated one-loop
OS component: the massive fixture agrees at exact control points, the massless
specialization agrees with $H_1$, and the integrated-backend matrix passes.
Physical proper one-loop children, completed two-loop parents in the full
graph domain, their finite analytic prefixes and the two-component disjoint
product pass comparisons in both direct orientation modes. The complete
OS/OS and IR-inside-OS sums and the figure-10 families also agree with projected
4D at the control points. Public direct support additionally admits a proper
two-loop OS component with exactly one classified connected proper descendant,
which must be a degree-one, one-loop OS component. The bare source, complete
child and finite-prefix comparisons pass in that domain with a live cograph.
Other proper multiloop histories remain guarded. The complete public
three-loop generation and orientation-summed route comparison passes in
1,773.62 seconds with the validated orientation grouping change. Its retained
record is `tests/artifacts/os_three_loop_grouped_routes_validation.json`.
The exhaustive numerical UV matrix also passes, separately from the exact
operator and reconstruction certificates.

#table(columns: (1.25fr, 1.15fr, 2fr), inset: 5pt,
  table.header([Component or forest], [Prescription], [Validation and guard]),
  [One-loop massive quark self-energy], [OS mass; MS-bar field],
  [Validated through local 4D, analytic correction and projected CFF.],
  [Massless quark specialization], [OS at $m=0$],
  [Validated exact agreement with degree-one soft refinement.],
  [Quark child inside gluon parent], [OS child; IR parent],
  [Ward certificate, nested analytic input, complete-family Taylor identities,
  exact CFF projection, all four hard subspaces, the external soft jet and
  production generation/evaluation pass. Complete CFF physical/rest-frame
  child-gluon radial controls also pass; arbitrary numerators excluded.],
  [Quark child inside quark parent], [OS child and parent],
  [Independent complete-family Taylor identities, exact nested CFF projection
  and all four hard subspaces pass. The complete four-stratum radial CFF soft
  matrix also passes in the stated below-threshold domain.],
  [Gluon child inside quark parent], [IR child; OS parent],
  [Independent complete-family Taylor identities, exact nested CFF projection
  and all four hard subspaces pass. The gluon-child constant/linear soft jet
  and the complete radial CFF soft control pass.],
  [Disjoint quark self-energies], [OS on each component],
  [Independent component identities, owned product, live bridge, exact CFF
  projection and left/right/both/correlated UV limits pass. Full CFF left,
  right and joint radial soft checks also pass.],
  [Proper two-loop quark component in a three-loop graph], [OS at both inner levels],
  [Exact component, mass-point, cograph and four-stratum UV certificates
  pass. Bare, completed-child and finite-prefix component-domain checks agree
  across all three routes (72.70, 316.35 and 374.70 seconds). The complete
  eight-family certificate and all fourteen ordinary/correlated hard strata
  pass. Public direct dispatch admits exactly one primitive OS descendant.
  Full public generation and orientation-summed agreement pass across all
  three routes. The exhaustive numerical UV matrix also passes for all
  26 LMBs and 182 nonempty hard-coordinate subsets.],
  [Ordinary UV or pole parent above OS], [MUV or PolePart],
  [Public dispatch guarded: a general mixed integrated policy has not been established.],
  [Massive scalar/vector or analytic chiral subgraph], [OS],
  [Unsupported fermionic projector or dimensional numerator prescription.],
  [Direct 3D, primitive and small forests], [Supported OS/IR combinations above],
  [Both modes pass primitive, proper one-loop child, complete two-loop parent,
  finite-prefix and two-component disjoint-product comparisons. A proper
  two-loop OS component with exactly one primitive OS descendant is admitted;
  other proper multiloop histories, disjoint products containing a multiloop
  component, and unions of three or more OS components remain guarded.],
)

Write $H_Gamma=U_2+S_1-U_2 S_1$ for the soft-refined gluon operation.
For a quark child $gamma$ in a gluon parent $Gamma$, test each of
$X$, $-T_gamma X$, $-H_Gamma X$, and $+H_Gamma T_gamma X$ before their sum.
For a quark child inside a quark parent, replace $H_Gamma$ by $T_("OS"),Gamma$.
The soft-refined gluon-inside-quark case has the opposite operator order.
Disjoint components commute only on their owned factors; their signed completed
counterterms and integrated prefixes use the established disconnected replay.
No global OS substitution is permitted on the shared direct-CFF root.

Every OS parent passes only its signed vacuum remainder
$-(U_1-U_1 O)[X_("completed")]$ to the analytic scheme machinery, including
when its completed child used an IR prescription. Its local operation still
contains the complete $-(U_1+O-U_1 O)[X_("completed")]$. The specialized IR
parent policy below has a separate certificate; neither policy follows solely
from the child's prescription.

== Derived Ward policy for the figure-10 insertion

This derivation applies to the completed one-loop quark counterterm inserted
into an equal-mass quark loop with two external gluons and canonical vector QCD
vertices. It is not a general permission to discard a gluon soft jet. Keep the
Lorentz dimension $d=4-2epsilon$ symbolic and the spinor trace normalization
four throughout.

After integrating the child, spatial parity removes its spatial-vector and
$gamma^0 gamma^i$ terms. Reversing the energy and exchanging the two OS branches
removes the remaining $gamma^0$ term. The child therefore contributes the two
insertion structures

$ [O_gamma]=delta m_("loop") 1, quad
  [V_gamma]=c(slash(r)-m), quad c=(2-d)/2 I_2(M). $

The on-shell mass integral remains an unknown scalar in this proof. Its value
is not needed. The fundamental color factor is proportional to the identity
and commutes with the two outer QCD generators.

Let $B^(mu nu)(q)$ denote the stripped outer one-loop gluon bubble. At zero
external momentum, define $D=r^2-m^2$ and $J_n=[D^(-n)]$. Direct Dirac algebra
and angular integration give

$ B^(mu nu)(0)=4g^(mu nu)
  (2/d (J_1+m^2 J_2)-J_1)=0. $

The last equality follows from the vacuum integration-by-parts identity
$2m^2J_2=(d-2)J_1$. At $q=0$ the mass insertion on either of the two fermion
segments is equal by cyclicity and loop reversal. Each is half the mass
derivative of the bubble:

$ 1/2 partial_m B^(mu nu)(0)
 =4m g^(mu nu)((4/d-1)J_2+4m^2/d J_3)=0, $

using $4m^2 J_3=(d-4)J_2$. The field insertion cancels one adjacent inverse
fermion propagator and reduces to $c B^(mu nu)(q)$, with the original common
Feynman-rule phases retained. The linear external-momentum coefficient of
either integrated insertion vanishes as a Lorentz tensor of rank three with
no preferred vector. In particular, the child's rest-frame direction has
already disappeared through the parity argument.

These are exact meromorphic identities in $d$, so they hold at every required
Laurent order. Thus the parent's $S_1$ contribution to the completed OS-child
insertion integrates to zero without evaluating a non-vacuum child mass
integral. For this certified nested family the analytic payload can use the
signed vacuum remainder $-(U_2-U_2 S_1)[X_("completed")]$. The complete
local parent operation remains $-(U_2+S_1-U_2 S_1)[X_("completed")]$; the Ward
identity changes only the analytic scheme input. A normalized frozen child
kernel integrates to one; its owned finite coefficient retains the
positive epsilon orders required by subsequent poles.

This argument does *not* set $[S_1 X_("bare")]=0$ for the bare two-loop
self-energy-insertion graph. Its inner bare self-energy has momentum-dependent
form factors, and transversality may require other diagrams. The existing IR
policy for the bare parent family therefore retains the generic full
$-H_Gamma[X_("bare")]$ analytic input, without the Ward omission. The topology,
canonical numerator, color flow and equal masses are prerequisites for the
specialized nested policy. Tests now validate the rule/propagator certificate,
the independent all-$d$ bubble and mass-insertion identities, the signed nested
analytic input and the positive-epsilon finite child prefix. Complete-family
Taylor identities, the four hard-subspace checks and the relevant external
soft jet now pass. The complete CFF physical/rest-frame child-gluon radial
controls also pass. Local forest identities alone do not establish the
integrated policy.

The logarithmic box component on edges 2, 3, 5, 6 of the figure-10 fixture has
a cograph on edge 4 independent of the outer external momentum. Its child-only
and IR-parent-composed forest families cancel because the parent soft constant
coefficient reproduces that complete insertion. The independent fixture test
now certifies this cancellation, including equality of the parent's $U$ and
$U S$ branches, before reducing the full forest to the four displayed OS/IR families;
the correlated hard region $k=p+r$ with $p$ hard and $r$ fixed keeps the
child gluon momentum finite and is covered independently in the UV checks. The
other independent strata are hard child momentum at fixed parent momentum,
hard parent momentum at fixed child momentum, and simultaneous generic hard
scaling. Gluon soft limits require their own checks; exact agreement of the
four local families is not a UV or IR certificate.

The small scalar contour checks distinguish the physical soft gluon
$k-p -> 0$ from an OS branch's rest-frame soft gluon $k-sigma p_("OS") -> 0$.
They first match the actual child pole polynomials and its affine energy
numerator. With the clockwise $d k^0/(2pi i)$ contour and unit spatial ray,
the OS-point scalar coefficients satisfy
$rho^2 R[1] -> 1/(4m)$ and $rho^2 R[k^0] -> sigma/4$.
Moving physical controls at $m=12$, $p=(15,9,0,0)$ also retain nonzero exact
limits on parallel and transverse rays. These certify the child scalar
kernels. The complete OS/OS CFF sum separately passes the four radial strata
below after all parent residues. The OS-inside-IR sum passes its two
child-gluon strata; the IR-inside-OS outer-gluon check also passes.
The scalar checks alone do not claim full IR finiteness.

For the nested quark forest, the complete inner remainder has a stronger
identity at the #emph[rest-frame] parent subtraction point. Write
$n=(1,0,0,0)$ and $p=sigma m n$. The ordered projectors and the graded relation
$U_1 O=O U_1$ imply

$ P_sigma (1-T_("OS"),gamma) Sigma(sigma m n)=0. $

In particular, the vacuum remainder vanishes under this projection. Its finite
analytic addback is proportional to $slash(p)-m$, so it obeys the same identity,
also when multiplied by its frozen scalar localization kernel. At the outer
soft point, the two adjacent quark propagators would otherwise supply a
double near-on-shell pole. The projected inner remainder removes its leading
coefficient and hence one power of the soft denominator: an isolated insertion
may behave as $rho^(-3)$ after energy integration, while the completed forest
is bounded by $rho^(-2)$ on this single-soft ray, with the other loop fixed
away from its singularities.

The numerical rays use $q=(5,0,0,0)$ and $m=12$, below the bare graph's massive
threshold, and keep the nonsoft loop generic. This argument does not assert a
pointwise boosted OS condition: the local operator is defined at rest. The
moving physical child-gluon limit belongs to the complete CFF orientation and
forest sum; its boundedness must not be inferred from one isolated residue.

Mixed-prescription guards can be relaxed only for those cases actually derived
and tested. Required coverage includes all relevant hard subspaces and the
corresponding gluon soft limits, rather than only uniform scaling of every loop.
Exact small forest terms precede the large graph tests.

The direct representations implement the same isolated-component operation
and agree with the projected-four-dimensional route after summing orientations.
The small proper-child, completed-prefix and full mixed-forest checks have
separate ownership and equivalence certificates. The proper two-loop component
with one primitive OS child also passes these three-route certificates in the
larger graph. Route measurements choose the supported representation for the
large scaling matrix; small equivalence checks remain for the other routes.

The completed scalar evaluator groups conditional summands by their exact
orientation condition at each additive node. It preserves products, powers
and lazy branch boundaries; each changed pass removes at least one conditional,
so the fixed-point iteration terminates. Four focused tests validate exact
symbolic selector truth tables, same-key sums, factorization and idempotence.
The three selector regressions and the primitive per-host comparison across
all three routes also pass. These checks certify the grouping transformation
independently of the completed large public generation and route comparison.

== Performance measurements and route selection

The isolated `ci-optim` measurements from 2026-10-13 use one fresh process per
fixture and route, binary SHA-256 prefix `c46b91da6603`, with orientation branch
grouping and generation profiling enabled. All six cases pass. The primitive
fixture includes local and integrated OS terms; the figure-10 fixture measures
its local forest. The host is shared, so process isolation does not imply an
idle machine. The complete record, profile hashes and raw reports are retained
in `tests/artifacts/os_isolated_performance_grouped_20261013/summary.json`.

The table separates counterterm orchestration from CFF generation and the two
evaluator construction stages. All times are milliseconds. The CT column
reports inclusive orchestration followed by its remainder after subtracting
profiled CFF work inside that boundary. This remainder includes forest
traversal, Taylor algebra, analytic integration, surface conversion and
assembly; it is not a pure Taylor timer. CFF totals include raw requests,
exact-source preparation, numerator templates, certificate construction and
assignment selection, including rejected native proposals exactly once.
Native CFF time is a subset of that total. These boundaries must not be added
to the inclusive CT timer. Compiled evaluator construction is zero in this run.

#table(
  columns: (1.45fr, 1.25fr, 1.1fr, 1fr, 1fr),
  table.header([Fixture / route], [CT inclusive / excluding CFF],
    [CFF total / native], [Spenso build], [Symbolica build]),
  [Primitive / projected], [154.89 / 132.77], [22.36 / 1.06], [24.52], [2.80],
  [Primitive / explicit], [146.44 / 146.32], [0.35 / 0.23], [33.22], [2.59],
  [Primitive / orientation], [150.19 / 150.05], [0.39 / 0.25], [32.68], [9.77],
  [Figure 10 / projected], [1,695.23 / 859.90], [836.58 / 53.39], [1,782.31], [283.84],
  [Figure 10 / explicit], [3,067.14 / 3,057.90], [10.51 / 1.69], [93,115.01], [2,007.27],
  [Figure 10 / orientation], [2,922.10 / 2,907.60], [15.76 / 1.91], [91,188.26], [2,588.99],
)

Expression sizes count the stored pre-evaluator roots and retained function
definitions without unfolding them. The Linux `VmHWM` field covers process
start through report capture, including model import. The separate watchdog
samples process-tree RSS every 250 ms through the full test and may miss brief
peaks; it is not interchangeable with the process high-water mark. Memory
columns are MiB; expression columns are exact byte counts and definition counts.

#table(
  columns: (1.45fr, 1fr, 1fr, 1.15fr, 1.25fr),
  table.header([Fixture / route], [VmHWM], [Sampled tree RSS],
    [Root bytes], [Definition bytes / count]),
  [Primitive / projected], [47.19], [47.61], [1,647], [7,976 / 6],
  [Primitive / explicit], [49.25], [49.67], [4,486], [6,621 / 10],
  [Primitive / orientation], [50.08], [50.47], [4,600], [6,621 / 10],
  [Figure 10 / projected], [269.49], [278.94], [305,815], [104,169 / 28],
  [Figure 10 / explicit], [2,663.73], [2,679.01], [1,031,437], [96,223 / 75],
  [Figure 10 / orientation], [2,695.26], [2,706.63], [1,139,938], [96,223 / 75],
)

The existing `bench` command uses `SingleParametric` evaluation in
`minimal-integrand` mode, which disables point caching, stability rotations
and event processing. Each of two fixed nonsingular momenta has five repeats
of twenty batches. The two-second target calibrates a sample count from ten
warmups; it does not guarantee two seconds of measured work. The last column
reports actual timed wall time summed over all ten repeats. The evaluator and
total columns give the range of the two point means in µs per sample, with
sample-weighted batches before averaging repeats.

#table(
  columns: (1.45fr, 1.2fr, 1.2fr, 1fr),
  table.header([Fixture / route], [Evaluator (µs)], [Total (µs)], [Timed wall (s)]),
  [Primitive / projected], [1.699–1.703], [26.70–26.80], [1.26],
  [Primitive / explicit], [1.328–1.330], [24.48–24.49], [1.24],
  [Primitive / orientation], [2.683–2.686], [27.15–27.23], [1.33],
  [Figure 10 / projected], [68.02–69.09], [102.95–104.37], [4.16],
  [Figure 10 / explicit], [95.67–96.92], [132.70–134.91], [4.88],
  [Figure 10 / orientation], [900.60–912.96], [990.16–1,002.34], [13.68],
)

Explicit direct has the fastest observed primitive kernel and total means.
Projected 4D is preferred for the extensive figure-10 matrix: it has lower
construction cost, memory, stored size and evaluation cost at both points.
Its evaluator repeat standard deviations are 2.29–3.69 µs, versus 4.18–5.63 µs
for explicit direct. The orientation-local nested evaluator has standard
deviations of 33.17 and 6.48 µs, with observed repeat means spanning
882–958 µs. Repeat variation within one process is not independent-run
uncertainty, and these diagnostics do not replace the existing paired
performance gate. The retained earlier run is a separate baseline; changes
between these single observations do not isolate the grouping change's cost.

The three other local-only forest matrices also pass with binary
`c46b91da6603`: OS inside OS in 41.96 seconds, IR gluon inside OS in
43.77 seconds, and the disjoint OS pair in 18.16 seconds. Their pass record is
`tests/artifacts/os_remaining_forests_grouped_validation.json`. They use the
same two-point, five-repeat, twenty-batch protocol. The table gives the range
of the two point means in µs per evaluator call.

#table(
  columns: (1.5fr, 1fr, 1fr, 1fr),
  table.header([Local forest], [Projected 4D], [Explicit direct], [Orientation-local]),
  [OS inside OS], [53.10–54.01], [53.87–54.28], [780.65–784.81],
  [IR gluon inside OS], [281.46–282.18], [37.89–38.37], [795.22–799.92],
  [Disjoint OS pair], [7.86–8.22], [15.51–16.07], [25.30–25.63],
)

The complete generation profiles include the CT and CFF events. Using the
same stage boundaries as above, all entries below are milliseconds; the
last column separates Spenso from Symbolica evaluator construction. Compiled
evaluator construction is zero for all nine routes.

#table(
  columns: (1.4fr, 1.3fr, 1.15fr, 1.35fr),
  table.header([Forest / route], [CT inclusive / excluding CFF],
    [CFF total / native], [Spenso / Symbolica build]),
  [OS/OS / projected], [2,169.88 / 815.50], [1,355.62 / 14.39], [967.05 / 307.62],
  [OS/OS / explicit], [1,211.11 / 1,202.21], [10.09 / 1.62], [5,212.16 / 808.89],
  [OS/OS / orientation], [1,181.03 / 1,172.15], [10.08 / 1.64], [6,224.14 / 1,340.64],
  [IR/OS / projected], [1,940.80 / 980.62], [961.44 / 18.06], [1,012.88 / 2,772.19],
  [IR/OS / explicit], [3,011.78 / 3,011.78], [1.20 / 1.16], [2,020.58 / 184.21],
  [IR/OS / orientation], [3,055.08 / 3,055.08], [1.20 / 1.16], [2,647.69 / 339.02],
  [Disjoint / projected], [703.51 / 317.09], [386.89 / 2.48], [217.95 / 39.47],
  [Disjoint / explicit], [521.15 / 460.17], [61.34 / 1.94], [3,893.88 / 165.73],
  [Disjoint / orientation], [500.81 / 439.97], [61.20 / 1.95], [3,944.57 / 262.52],
)

Each fixture runs its three routes sequentially in one process. Generation RSS
and lifetime high-water marks can retain prior-route allocations; these are
not isolated memory comparisons. The complete per-route stage breakdowns,
actual timed evaluation durations and repeat variation, sampled RSS and
`VmHWM`, and factorized root/definition byte counts are retained in each
`summary.json` below
`tests/artifacts/os_remaining_forests_grouped_20261013/`: the fixture directories
are `os_quark_nested_route_matrix`, `os_soft_gluon_route_matrix`, and
`os_disjoint_route_matrix`. None of these size measurements expands stored
numerator functions.

Explicit direct remains the clear evaluation choice for IR inside OS;
projected 4D is the clear evaluator choice for the disjoint pair. For
OS inside OS, the evaluator means are close: projected 4D is 0.5–1.4 percent
lower, within the measured variation. Its generation takes 3.49 seconds
versus 7.28 for explicit direct, and its total sample means are lower.
Projected 4D therefore remains the provisional choice for that extensive
matrix. These choices retain small equivalence checks for the other routes
and make no isolated-memory claim from the shared-process matrices.

The full three-loop local forest also passes public generation and
orientation-summed pointwise agreement across all three routes with binary
`c46b91da6603`. Its complete test takes 1,773.62 seconds. The profile and
route summaries are retained in
`tests/artifacts/os_three_loop_route_matrix/summary.json`. Construction times
below are seconds; the evaluator column gives the range of the two point means
in milliseconds per sample, using the same repeat protocol. The CT and CFF
columns retain the nonadditive boundaries defined above. Compiled evaluator
construction is zero.

#table(
  columns: (1fr, 1.2fr, 1.15fr, 1.35fr, 1fr),
  table.header([Route], [CT inclusive / excluding CFF], [CFF total / native],
    [Spenso / Symbolica build], [Evaluator (ms)]),
  [Projected 4D], [65.88 / 25.04], [40.857 / 0.257], [53.58 / 31.31], [2.898–2.900],
  [Explicit direct], [73.50 / 73.44], [0.073 / 0.026], [559.58 / 49.70], [2.180–2.184],
  [Orientation direct], [71.98 / 71.92], [0.079 / 0.027], [599.55 / 119.08], [346.47–354.47],
)

Total generation takes 150.97, 683.10 and 790.75 seconds for projected,
explicit and orientation-direct routes respectively. The stored factorized
roots plus definitions occupy 11,109,441, 29,292,094 and 32,500,232 bytes,
with 88, 152 and 152 definitions. The artifact retains their separate root
and definition counts without unfolding them. All routes share one process,
in that order: recorded lifetime `VmHWM` values of 3.30, 11.96 and 17.20 GiB
include allocations from earlier routes. They are observed run peaks, not
isolated per-route memory requirements. Actual timed evaluation totals are
18.18, 17.35 and 71.07 seconds; the expensive orientation evaluator reaches
the minimum twenty samples per repeat.

The scalar orientation grouping step merges 36,030 branches in two passes,
reducing its input from 124,346,232 to 123,750,806 bytes in 3.688 seconds.
The previous pairwise matcher spent more than an hour in repeated merging;
that diagnostic was stopped after 5,193.13 seconds and is recorded as
superseded in `tests/artifacts/os_three_loop_routes_validation.json`.
It supplies no completed baseline timing or whole-run speedup ratio. The new
grouping allows the full route comparison to complete, while the roughly
350 ms orientation-local evaluator remains a substantial cost.

The exhaustive numerical UV matrix passes on projected 4D: 26 physical LMBs
times seven nonempty hard subsets give 182 tested limits. Each uses 25 scales
from $10^8$ to $10^12$, for 4,550 samples plus a separate matched bare control.
The successful record is `tests/artifacts/os_three_loop_uv_scaling_validation.json`.
This establishes numerical local UV behavior across those coordinate limits;
the fourteen exact four-dimensional ordinary/correlated hard-stratum
certificates establish their algebraic statements independently.

Explicit direct has the faster evaluator at
the generic benchmark points, but projected generation saves approximately
532 seconds. Multiplying the observed total-sample difference by 4,550 gives
only about 3.4 seconds, so projected 4D is the provisional choice for the
fastest complete check. This estimate is not a measurement of UV-ray
throughput: hard rays and arbitrary-precision retries can have different
costs. The first passing run takes 211.58 seconds with an unrestricted worker
pool on the 384-CPU host. One observation records 41.0 GiB of resident memory
(42,993,496 KiB); it is a sampled lower bound on the maximum, not a measured
peak. The maintained test scopes the profiler to one worker without changing
any limit, scale or precision retry. This bounded-resource run passes in
542.56 seconds with binary `63a324ef017d`; both the subtracted and matched
bare profiles cover all 182 limits. The 30 GB watchdog completes successfully
in 543.50 seconds and records a sampled process-tree maximum of 3,632,726,016
bytes (3.383 GiB). Its sampling can miss brief peaks; this is a different
measurement scope from the earlier single-process RSS observation.
The pass record and memory trace are
`tests/artifacts/os_three_loop_uv_single_worker_validation.json` and
`tests/artifacts/os_three_loop_single_worker.memory.jsonl`. The separate
`subtracted.json` and `bare.json` reports are retained beneath
`tests/artifacts/os_three_loop_single_worker/os_uv_scaling_local_os_three_loop_top_self_energy/`.

Ordinary local and NixCI cases retain the four-minute failure bound. Larger
physical and nested integration workloads belong to the main-only `ci_slow`
selection: add their exact names both to the fast profile's long-case exclusions
and to that slow selection, while keeping the tests nonignored. Main runs the
full core selection and four long integration partitions. This scheduling does
not alter the paired performance gate: relevant, non-draft, same-repository
`final-review` PRs and relevant merge-group commits run the existing comparison;
main pushes skip that comparison, and manual dispatch remains available. The
trigger and workload policy is maintained in
#link("ci.typ")[CI maintenance and measurements].

The physical end-to-end acceptance covers $gamma^* arrow.r t overline(t)$
through NLO with the self-energy and vertex master graphs, all $t overline(t)$
and $t overline(t)g$ cuts, and all orientations. The independent analytic LO
reference is 1.387581507745258 and the NLO correction is 0.4328300867769855 at
$E_("cm")=400$, pole mass $m=173$, $mu_r=173$, $alpha_s=0.118$ and
$alpha_("EW")^(-1)=132.507$. The LO control gives
$1.40582 plus.minus 0.02767$ from 10,000 samples.

All three routes pass the integrated reference and auxiliary-scale checks at
$(M,L)=(173,173),(346,173),(173,346)$. Each case uses 320,000 samples with
independent seeds 1337–1345, relative-spherical coordinates, map-density channel
weights, and explicit sums over graphs, orientations and sampling channels.
The Vegas continuous-dimension learning rate is 0.25; batches start at 20,000
samples and increase by 20,000, with 20,000 samples required before a grid
update. Ten workers evaluate each case. The complete acceptance takes 486.70
seconds and belongs to the full main-push lane, outside the ordinary
four-minute test selection.

#table(
  columns: (1fr, 1fr, 1fr, 1fr),
  table.header([Route], [$(M,L)=(173,173)$], [$(346,173)$], [$(173,346)$]),
  [Projected 4D], [$0.43331 plus.minus 0.01225$],
    [$0.42651 plus.minus 0.01250$], [$0.41210 plus.minus 0.01771$],
  [Explicit direct], [$0.43121 plus.minus 0.00869$],
    [$0.43260 plus.minus 0.00659$], [$0.43580 plus.minus 0.00795$],
  [Orientation direct], [$0.42854 plus.minus 0.01081$],
    [$0.42292 plus.minus 0.01094$], [$0.43058 plus.minus 0.01057$],
)

The largest uncertainty is 4.10% of the analytic NLO reference, below the
unchanged 5% requirement. The largest discrepancies from the reference,
within a route's scale scan, and between routes at the same scales are below
1.18, 0.99 and 1.23 standard deviations respectively; all satisfy the
three-standard-deviation requirements. Final measurements, settings and the
successful binary identity are retained under
`tests/artifacts/gamma_star_ttx_os_sampling_diagnostics/relative_spherical_learning025_320k/`.
Successful-test cleanup removes the per-case integration checkpoints, so these
results make no retained per-sample outlier diagnostic claim.

The routes also agree in arbitrary precision at momentum-space controls for
each master graph and at normalized-coordinate probes chosen from the large
weights. The latter separately compare actual mapped momenta, Jacobians and
channel partitions, and require the two graph-group contributions to sum to
the complete value. These are numerical controls of the sampling and evaluation
boundaries; the exact algebraic reconstruction certificates above establish
the symbolic identities independently.

The integrated observable tests the primitive OS mass prescription,
inclusive real/virtual cancellation and auxiliary-scale independence only
at the summed and integrated boundary. It cannot certify nested forests or
independently separate field and coupling schemes; their local and analytic
certificates above remain necessary.
