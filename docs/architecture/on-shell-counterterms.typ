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

The substitution $p mapsto sigma p_("OS")$ belongs only to the active component.
An edge crossing its boundary may carry $r+p$ in the compatible sub-LMB; only
that component's occurrence is evaluated at $r+sigma p_("OS")$. External
cograph factors and other components keep their own momenta. A global
replacement of the symbol $p$ would change the surrounding graph and is invalid.
The component boundary index, directed fermion flow and physical mass identify
the operation; a generic two-point scalar graph is insufficient.

== Scaling and ownership

The rest-frame point is the expression $m(1,0,0,0)$, never an unrelated frozen
four-vector. In a subsequent UV Taylor jet, physical masses and external
momenta have weight one, so $m mapsto t m$ implies
$p_("OS") mapsto t p_("OS")$. The completed $O$ branches consequently contribute
their first-order mass and energy-shift coefficients to $U_1 O$. In the
implementation's inverse-hard convention this is expressed through its existing
hard-loop rescaling and physical-mass bookkeeping. The vacuum mass keeps its
existing hard weight; it must not be identified with the physical mass.

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
  +epsilon (L_M^2/2-L_M+pi^2/12)+O(epsilon^2)). $

The finite $-1$ inside this bracket comes from the dimensional numerator. An
outer $1/epsilon$ pole requires the displayed positive epsilon coefficient;
the existing forest loop-order bound determines the actual required depth.
Do not truncate every completed child to epsilon zero. Backend insufficiency
must produce an error.

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
 +epsilon(3L_m^2/2+4L_m+8+pi^2/4)+O(epsilon^2)). $

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

The initial acceptance target is the first graph: local OS subtraction,
explicit vacuum-only scheme correction, exact CFF reconstruction, an independent
residue comparison and the integrated mass/field references above. The support
boundary must remain narrower than the algebraic proposal until these checks
run through production. In particular, configuration acceptance or a
scalarized toy identity alone does not establish open-spinor support.

#table(columns: (1.15fr, 1.55fr, 2fr), inset: 5pt,
  table.header([Component], [Prescription], [Required boundary]),
  [Massive quark self-energy], [OS mass; MS-bar field],
  [Degree one, one identified open fermion chain, component-local point.],
  [Massless quark self-energy], [Massless OS specialization],
  [Must agree exactly with the existing soft-refined degree-one operator.],
  [Gluon/ghost self-energy and supported vertices], [Existing soft-refined or ordinary UV],
  [Retain the existing component-specific prescription.],
  [Massive scalar or vector self-energy], [OS],
  [Unsupported: no fermionic projector prescription applies.],
  [Chiral or gamma-five analytic subgraph], [Integrated OS],
  [Unsupported until a dimensional numerator prescription is supplied.],
  [Mixed or nested components], [OS with soft refinement or ordinary UV],
  [Enable only combinations with completed-child proofs and limit coverage.],
)

For a quark child $gamma$ in a gluon parent $Gamma$, test each of
$X$, $-T_gamma X$, $-H_Gamma X$, and $+H_Gamma T_gamma X$ before their sum.
For a quark child inside a quark parent, replace $H_Gamma$ by $T_("OS"),Gamma$.
The soft-refined gluon-inside-quark case has the opposite operator order.
Disjoint components commute only on their owned factors; their signed completed
counterterms and integrated prefixes use the established disconnected replay.
No global OS substitution is permitted on the shared direct-CFF root.

Mixed-prescription guards can be relaxed only for those cases actually derived
and tested. Required coverage includes all relevant hard subspaces and the
corresponding gluon soft limits, rather than only uniform scaling of every loop.
Exact small forest terms precede the large graph tests.

The direct representations must implement the same component operation and
agree with the validated projected-four-dimensional route after summing
orientations. Only then choose the fastest route for the large scaling matrix;
retain small equivalence checks for the other routes. Ordinary local and NixCI
cases remain below four minutes. Larger physical and nested workloads belong
to the full main-push suite, following the existing performance-gate trigger.
Measure counterterm construction, CFF generation, evaluator construction,
evaluation, peak memory and expression sizes separately for every route. The
physical end-to-end target begins with the required contributions to
$gamma^* arrow.r t overline(t)$; auxiliary-scale independence is checked at its
summed and integrated boundary.
