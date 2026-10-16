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
pass records and measurements belong to their stated source revisions,
including the original OS stack ending at `ce58d97de847` and the earlier
FeynKit rebases recorded below. They do not certify the retained-tensor port.
Current validation records must identify the tested revision, test selection
and run; exact reconstruction, local limits and integrated renormalization
remain separate claims.

At retained-tensor source `bb4ee859bb121796015db6ba3e4a9faa40089cff`, the
finite-density case passed in run `a9ef128e-92bc-4870-b449-b0eb9aac73a0`
(25.472 seconds): eighteen grouped-master comparisons retain both charge-mate
source weights, their independent Born sum and both vertex spectator charts.
The density fixture explicitly selects its authored input bases before grouping:
GL0 uses edges $(1,2)$ and both raw self-energy sources GL1/GL2 use $(1,4)$.
The retained grouping, signed-routing and exact basis assertions remain in force.

The expanded route-map sampling, published absolute MS-bar acceptance and its
independent reference anchor passed in run `925d775b-4c03-4248-938f-bbcaad9995e2`
(42.737 seconds). The sampling matrix contains two configurations, three points
and three routes with the complete native role/channel inventory. All four
cases use the unchanged ordinary 240-second per-case limit. Their source,
binary and exact comparison records are retained under
`target/os-feynkit/retained-tensor-os-rebase/08-finite-reinsertion/all-source-chart/`.
The nine-case inclusive OS workload remains deferred to the existing full-main
gate; these ordinary checks and the earlier source/run identities below do
not qualify that workload on this port.

At retained-tensor source `b4c11e911594acd1d3a91a992efc0a156bbc9055`, the
three open-component controls passed in run `36c16895-4b5e-4d12-8fb3-52d0e6806127`
(44.174 seconds). Exact old/new fixture imports also retain every component
tensor and logical interface, local factor, denominator and signed basis;
only native external fermion attachment closes the final scalar interface.
The soft-gluon scalar route case passed in run
`a198f3ef-a685-4544-8454-a97a663e9a56` (70.441 seconds), comparing the three
local-only routes at two points under the existing full-main 7200-second cap.
The other five affected scalar-route and complete-UV cases remain deferred.
Current source, projector identities and route reports are bound under
`target/os-feynkit/retained-tensor-os-rebase/09-scalar-fixtures/`.

The rebased boosted implementation at `5dc4d91fa9ea` passed 34 focused cases
in nextest run `1ccd928f-3bdf-4c13-a511-c9d377138076`. These cover the independent
one-loop operator and CFF certificates, local UV coefficients and complete
rainbow soft limits, and the integrated pole condition for all nine pairs of
subtraction and external velocities. The three-route comparisons additionally
cover a proper child with a nonaxial point and unrelated live external momenta;
these pointwise comparisons supplement the exact certificates. The slowest
case, the nonaxial complete rainbow, took 69.7 seconds under the ordinary
240-second limit. This record does not certify the physical inclusive-rate or
large performance and scaling workloads on the rebased stack. The later
public-generation and performance records below identify their own revision
and measurement scopes.

The historical temporal-trace stack `6a5a3c3faf7bb662998a8b402549aea7f0573d19`
passed the one-loop mechanism, maintained nested algebra and local limits,
all six CLI three-route comparisons, and the nine-case physical acceptance.
The detailed historical record below keeps these claims separate. At that
checkpoint, exhaustive numerical UV scans and final isolated performance
measurements remained pending the consolidated FeynKit repair rebase; those
records did not certify that later revision or the current retained-tensor port.

== Ordered component operation

Use the Minkowski metric $(+, -, -, -)$, $d=4-2 epsilon$ and
$slash(q)=gamma^mu q_mu$. A fermion self-energy is an ordered linear map from
its incoming spinor slot to its outgoing slot. Define

$ p_("OS") = (m,0,0,0), quad P_sigma = (1 + sigma gamma^0)/2, $
$ O[Sigma] = sum_(sigma=plus.minus 1) P_sigma Sigma(sigma p_("OS")), $
$ T_("OS")[Sigma] = U_1[Sigma] + O[Sigma] - U_1[O[Sigma]]. $

The rest frame is the default of a fixed-frame prescription. The generation
setting `uv.on_shell_subtraction_velocity` accepts four exact rational strings
in the order $(v^0,v^1,v^2,v^3)$, with $v^0>0$ and $v^2=1$. For example,
`["5/4", "3/4", "0", "0"]` selects an axial boost and
`["5/4", "1/4", "1/2", "1/2"]` a non-axial boost. The point and projector
are changed together:

$ p_*=m v, quad P_sigma(v)=(1+sigma slash(v))/2, quad
  O_v[Sigma]=sum_(sigma=plus.minus 1)P_sigma(v)Sigma(sigma m v), $
$ T_v=U_1+O_v-U_1 O_v. $

Here $slash(v)=v^0 gamma^0-sum_(i=1)^3 v^i gamma^i$. Its four fixed
components are embedded in the symbolic $d$-dimensional Lorentz algebra;
the spinor slots remain four-dimensional. The velocity is fixed
generation data, not a function of a live external or loop momentum. Every
component retains its own point owner and physical mass even when components
share the velocity. The default strings are `["1", "0", "0", "0"]`.
Nonconstant, nonrational, past-directed and nonunit input is rejected. This
exact-input restriction avoids normalizing an approximately on-shell vector.

Rotational stability probes must also respect this fixed frame. Rotating the
live momenta while holding a nonzero $bold(v)$ at a different orientation
changes the local integrand; agreement is restored by jointly rotating the
reference or by using rotations that preserve it. The current representation
keeps the reference fixed. Generated graph payloads retain its exact spatial
components, and runtime validation accepts only supported probes that preserve
every retained reference. For the non-axial example, use
`rotation_axis = [{ type = "axis_angle", axis = [1.0, 2.0, 2.0], cos_theta = 0.0 }]`
in the runtime stability settings. This is a nontrivial quarter-turn about
$bold(v)$; the axis is normalized at the active precision. Its parallelism
with the generated reference is certified exactly, including the binary64
values supplied in the rotation settings. Identity and aligned coordinate-axis
rotations are also supported. Euler rotations are admitted when each nonzero
elementary rotation preserves the reference. The rest-frame default has no
spatial restriction. Precision escalation and accuracy thresholds are unchanged.

The reference metadata belongs to generated graph terms, so graph selection,
cloning and loading preserve the same constraint. Only an OS component with a
massive boundary fermion creates a fixed reference. A massive internal bubble
inside a massless external-quark self-energy does not create this constraint: its
subtraction point is still zero and its two projectors sum to the identity.
Combined state manifest version 22 records this layout on top of the fully
indexed version-21 interfaces. Both the older indexed-only version 21 and the
earlier partial-interface OS version 21 require regeneration; no positional
compatibility decoder is provided.

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

The fixed point is the expression $m v$, never an unrelated frozen four-vector;
the default has $v=(1,0,0,0)$. In a subsequent UV Taylor jet, physical masses and external
momenta have weight one, so $m mapsto t m$ implies
$p_("OS") mapsto t p_("OS")$. The completed $O$ branches consequently contribute
their first-order mass and fixed-momentum shift coefficients to $U_1 O$. In the
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
physical mass and an already constructed child fixed point unchanged. A
surrounding OS operation evaluates only its own boundary momentum; inner
points remain related to their own physical masses. Distinct
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

For a fixed velocity $v$, the same independent differentiation replaces
$k^0$ by $k dot.op v$ and the ordered $gamma^0$ by $slash(v)$:

$ V'_v=frac(2(N_0(k dot.op p)-m slash(v)N_0(k dot.op v)),D^3). $

The on-shell branch denominators become
$(k^2-m^2)((k-sigma m v)^2)$. Their spatial shifts belong to the same
component-owned point as their energy shift. CFF reconstruction must retain
both shifts; restoring an OS crown coordinate to a zero spatial vector is
valid only for the default rest velocity.

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
spinor trace normalization remains four. The identity $slash(v)^2=1$
gives the scalar mass term for every fixed unit velocity. Setting $d=4$ before multiplying the pole loses
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
that discard the required positive-epsilon child term. Both backends must agree
through $epsilon^2$. MATAD's symbolic mode retains $Gamma(1)$ and
polygamma constants; the comparison boundary uses the exact identities
$Gamma(1)=1$, $psi(1)=-gamma_("E")$, $psi^((1))(1)=pi^2/6$ and
$psi^((2))(1)=-2zeta(3)$, without numerical substitution. The extended reference
tests require exact agreement after this change of constant basis.

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

=== Independent finite addback in the photon forward graph

At NLO there is no enclosing UV pole in the selected physical cuts, so this
reference uses the finite coefficient at $epsilon^0$. Write
$a_C=a C_F=C_F alpha_s/(4 pi)$ and $L_M=ln(M^2/mu^2)$. The native quark
propagator is $S=i(slash(p)+m)/(p^2-m^2+i 0)$. The field insertion has coefficient
$-i a_C(1+L_M)(slash(p)-m)$. Its numerator identity yields the causal pinch
$S(slash(p)-m)S=i S$, so each self-energy contributes
$a_C(1+L_M)$ times its Born cograph. Its two adjacent virtual cuts are one
raised-residue group and do not supply an extra multiplicity.
This equality is taken in the causal limit; the exact numerator identity
contains $p^2-m^2$ and does not cancel a finite nonzero $i 0$.

The independent vertex average is
$ gamma_mu slash(q) gamma^nu slash(q) gamma^mu
  arrow.r (d-2)^2/d q^2 gamma^nu. $
The scalar vacuum integral with numerator $q^2$ and three equal denominators
has expansion $i/(16 pi^2)(1/epsilon-L_M-1/2)$, while
$(d-2)^2/d=1-3epsilon/2+O(epsilon^2)$. With the actual propagator phase
$i^2(-i)=i$ and quark-gluon coupling $i g_s$, each vertex therefore contributes
$-a_C(2+L_M)$ times its Born cograph. The two vertex components and two
self-energy components have the same external averaging and closed-loop
factor as the Born graph. Their complete finite addback integrates to

$ [2a_C(1+L_M)-2a_C(2+L_M)] "Born"=-2a_C "Born". $

For $m=mu=173$ and $E_("cm")=400$, this is $-0.0347455650616886$ in the
acceptance test's normalization. The OS mass contribution remains in the local
on-shell branch and is not included a second time in this finite reference.
Before integration, the density has an extra child-loop factor
$t_*^3 rho_L(t_* q)$ relative to the Born LU measure. The vertex graph has
two distinct Born momentum charts; they must both be included. The integral
of this factor is one, independently of $L$.

The grouped-density fixture checks both retained self-energy source weights
before comparison. Its independent reference sums the Born densities at $-p$
and $p$ with a nonzero parity control, rather than multiplying by an observed grouped
weight. Both vertex spectator charts remain separate. The test makes eighteen
grouped-master density comparisons at three points and three scale pairs,
retaining both self-energy source weights and both vertex charts. The historical
twenty-seven-comparison record below used three ungrouped masters. These are
the maintained test's coverage requirements; their retained-tensor
qualification is recorded separately.

== Fixed boosts and the rainbow soft condition

The one-loop integrated condition is independent of the chosen fixed velocity.
In the fixture's quark-loop routing, the Feynman-parameter moments are
$J_0=1/(1-2epsilon)$ and
$J_1=1/((1-2epsilon)(2-2epsilon))$. They multiply the constant numerator and
the difference between its value at $k=p$ and $k=0$, respectively. Applying
these moments directly to the production on-shell branches must give the
same scalar mass correction for every $v$. After the signed field pole is
included, $Sigma_R(p)P_+(p/m)=0$ for every future on-shell $p$, even when
$p != m v$. The whole matrix need not vanish: a term proportional to
$slash(p)-m$ is consistent with MS-bar field renormalization.

The integrated reference checks all nine combinations of rest, axial and
non-axial subtraction and external velocities. The local projected-shell,
UV and massless-specialization checks use each of the three subtraction
velocities. The production vacuum reference also repeats the two analytic
backends, two integration modes, two auxiliary masses and both retained
positive-epsilon depths at each velocity. These checks are independent
of the three-route comparison; agreement at sampled momenta does not replace
the analytic pole condition.

Local soft cancellation is a stronger, different statement. In the two-loop
rainbow fixture, write the outer gluon momentum as $k=r-p=rho K$ and leave
the child loop $ell$ symbolic. Two adjacent quark propagators carry $r$.
For a generic ray with $K^2 != 0$ and $p dot.op K != 0$, their product
with the gluon denominator starts at $rho^4 K^2 (2p dot.op K)^2$.
Its leading numerator contains

$ gamma_mu (slash(p)+m) R_("child")(p;ell)
  (slash(p)+m) gamma^mu. $

Here $R_("child")$ includes its local subtraction and its actual finite
analytic coefficient with the frozen localization kernel. The complete
forest also contains the overall parent subtraction and the parent operation
on the completed child. The overall analytic addback has a regular
localization kernel on this stratum; it cannot supply a negative Laurent
power. A matched fixed point imposes the local projected child condition.
A different fixed point imposes the integrated mass condition but does not
generically annihilate this numerator at each $ell$.

The exact tests compare the complete $rho^(-4)$ coefficient against this
independently constructed cograph insertion. The matched cases retain a
nonzero $rho^(-3)$ coefficient as a control. They cover rest, axial and non-axial
external velocities with rest and matching subtraction points. These are
four-dimensional Laurent certificates. They do not by themselves assert
absolute integrability at intersecting soft strata or establish the result
after CFF energy integration. The boosted one-loop reconstruction and route
checks certify that component's representation boundary; boosted two-loop
CFF soft behavior requires its own limit check.

The complete forest tests require the following exact coefficients at $m=12$ with
$p=(12,0,0,0)$, $(15,9,0,0)$ and $(15,3,6,6)$, retaining symbolic $K$
and $ell$ away from the other denominator zeros:

#table(columns: (1.5fr, 1fr, 1fr), inset: 5pt,
  table.header([External point / subtraction], [$rho^(-4)$], [$rho^(-3)$]),
  [Rest / rest], [Exactly zero], [Nonzero],
  [Axial boost / rest], [Nonzero], [Not needed for the leading limit],
  [Axial boost / matched], [Exactly zero], [Nonzero],
  [Non-axial boost / rest], [Nonzero], [Not needed for the leading limit],
  [Non-axial boost / matched], [Exactly zero], [Nonzero],
)

These tests include the actual finite child correction, its parent operation,
frozen localization ownership and the regular overall analytic addback.
Direct finite-component certificates preserve the same assertions without
expanding the graph-level diagnostic numerator. The fixed-frame public
generation test requires agreement across all three routes and preservation
of the stability constraint after state loading.

The original OS run recorded 7.1, 17.5 and 59.9 seconds for the rest, axial and
non-axial cases in `tests/artifacts/on_shell_boosted_validation.json`. These
are historical timings; the rebased stack requires its own ordinary-lane run.

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

An OS branch carries the fixed four-vector shift $sigma m v$: its energy
shift is $sigma m v^0$, and its spatial shift is $sigma m bold(v)$, with the
routing sign retained. At rest this reduces to the $k^0+sigma m$ shift in the
residue example below. The fermion remains a physical-mass propagator. A matching
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
contour step. Exact reconstruction and CFF projection must be certified for
all three physical nested fixtures and the disjoint fixture below. These
checks establish the source assignment and projection boundary. Independent
complete-family Taylor identities establish operator composition. Separate
small-graph limit checks cover all four hard subspaces below and the relevant
constant/linear IR soft jets; the OS/OS, OS-inside-IR, IR-inside-OS and disjoint
radial CFF controls test the stated domains. Counting points per component
alone does not solve the live-crown collision.

As a separate rest-frame contour oracle, let $omega=abs(bold(k))$,
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

The boosted Dirac oracle instead uses the fixture's quark routing. Define
$omega_sigma=abs(bold(k)-sigma m bold(v))$; the lower-half-plane poles are
$k^0=E$ and $k^0=sigma m v^0+omega_sigma$. Apply the ordered $P_sigma(v)$
to each residue and use the second derivative at the cubed vacuum pole for
$V'_v$. Exact control points use $m=8$, $M=20$, $E=17$ and vacuum energy
$25$. At rest, $bold(k)=(9,plus.minus 12,0)$ gives $omega_sigma=15$.
For the axial boost, $bold(k)=plus.minus(15,0,0)$; for the non-axial boost,
$bold(k)=plus.minus(5,10,10)$. The two branch energies are $(9,21)$ or
$(21,9)$, respectively. The six matrix-valued residue comparisons are required
in addition to the general exact numerator and denominator certificates.

== Algebra and operation regression boundaries

The OS stack builds on the copied FeynKit branch rebased onto main. Generic
algebra fixes and their independent regressions are separate changes below
the OS operations. Run
`direnv exec .. nix develop -c just test-fixed-frame-algebra`
from a workspace containing those changes. The gate covers Spenso, Idenso, GammaLoop's analytic
bridge and Vakint without requiring the OS component operation. It retains
the ordinary four-minute per-test limit.

#table(columns: (1fr, 2.7fr), inset: 5pt,
  table.header([Algebra boundary], [Exact regression]),
  [Fixed basis metric and Clifford algebra],
  [All sixteen basis products have signature $(+, -, -, -)$ in symbolic $d$.
  Each physical spatial slash anticommutes with $gamma^0$, squares to $-1$,
  and changes sign under temporal conjugation. General vectors, unsupported
  components and formal spin dimensions retain their guards.],
  [Ordered chain interfaces],
  [Dummy-label changes preserve ordered temporal joins. General chain
  transposition uses its explicit native operation; temporal-matrix symmetry
  does not authorize transposition of arbitrary gamma words.],
  [Native vacuum tensor ports],
  [Actual gamma and fixed-vector leaves retain connected Lorentz ports through
  angular projection. Restoring their native interface precedes symbolic-$d$
  algebra; dimensional unification precedes allocation of fresh dummy indices.],
  [Dimensional numerator],
  [Open and closed chains retain the contributions from $(d-4)$ and $(d-4)^2$
  multiplying single and double poles. Both integration policies preserve
  the enclosing operation's required positive epsilon orders.],
  [Compact coefficients],
  [Weighted dots, powered scalar spectators and native tensor coefficients
  retain their factorization in both tensor-reduction modes. Shallow exact
  numeric-content normalization includes implicit unit coefficients. A bounded
  collection in literal scalar factors preserves their complete bodies and
  rejects hidden or unsupported dependence.],
  [Residual contractions],
  [Unsupported internal Lorentz contractions remain errors per analytic
  branch. Scalar dimensional substitution leaves Lorentz-slot metadata intact;
  free tensor interfaces survive exact zero and equality checks.],
)

OS-specific tests sit above that algebra boundary. They check the production
projector on its outgoing component port against independent Weyl matrices,
including a noncommuting matrix unit; both projector signs, idempotence,
orthogonality and completeness; physical-mass scaling of every fixed-point
component; point lifting and restoration with an untouched live cograph;
and the composed local operator, integrated pole condition and complete
rainbow soft coefficients. These operation tests cannot be replaced by a
standalone Clifford identity.

The analytic bridge uses the native typed algebra and contraction owners.
It does not restore the superseded GammaLoop scalar-product alias machinery
or append a separate legacy gamma-chain pass after projection. Finite
coefficient certificates use the same exact literal-factor collection owner;
only complete scalar sums with a certified scalar interface may become
literal keys. Diagnostic expressions preserve graph numerator factors,
function bodies, inverse scopes and fractional-power branches. Numerical
agreement and coefficient sampling do not supply an exact zero certificate.

The finite route-comparison oracle retains the completed integrand's numerator
catalog. For compiled comparisons, one native tensor probe closes the graph's
free interface before shared-family preparation. Its independent component
parameters select each output at evaluation time, with the native dual-slot
order and metric signs. Only this small probe is materialized; graph numerator
bodies retain their factorization. Exact native component identities and an
independent numeric tensor contraction check this extraction separately from
the three physics routes. Compiling an evaluator and contracting a tensor at a
fixed point have different costs, so their construction and evaluation times
are recorded separately.

== Validation graphs and support boundaries

#import "../assets/typst/portal-graphs/layout.typ": layout as draw-graphs

// Draw the actual fixture incidence. Splitting the invisible external node is
// a presentation-only transformation; internal vertices and edge IDs survive.
#let show-fixture(path, native-sm: false) = {
  let input = read(path)
    .replace("ext [style=invis]", "ext_in [style=invis]; ext_out [style=invis]")
    .replace("ext ->", "ext_in ->")
    .replace("-> ext [", "-> ext_out [")
  // Native runtime DOT stores model-local IDs rather than display names.
  // These two fixtures use sm-default: 0=photon, 4=gluon, 25=top.
  // Add presentation metadata only; incidence, arrows and source slots survive.
  if native-sm {
    for (id, particle) in (("0", "a"), ("4", "g"), ("25", "t")) {
      input = input.replace(
        "particle_id=\"" + id + "\"",
        "particle_id=\"" + id + "\" particle=\"" + particle + "\"",
      )
    }
  }
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

#figure(show-fixture("../../tests/resources/graphs/local_os_native_gl1_self_energy.dot",
    native-sm: true),
  caption: [Actual native GL1 self-energy contribution to the photon forward
  graph. The OS component comprises edges 1–2; its matrix ports are hedges
  5 and 8, and incoming momentum belongs to edge 4. The runtime DOT retains
  the generated Feynman factors and local UFO port ordering.])

#figure(show-fixture("../../tests/resources/graphs/local_os_native_gl2_self_energy.dot",
    native-sm: true),
  caption: [Actual native GL2 charge-conjugate contribution. Its OS component
  again comprises edges 1–2; its matrix ports are hedges 7 and 6, with incoming
  momentum on edge 3. Each component uses its own directed fermion chart.])

#figure(show-fixture("../../tests/resources/graphs/local_nested_soft_top_self_energy.dot"),
  caption: [Existing nested validation topology from
  `local_nested_soft_top_self_energy.dot`. Internal edge 2 and child edge 3
  provide the compatible homogeneous loop basis. The OS-inside-OS forest
  tests and the fixed-boost rainbow soft tests use this figure-11 topology.])

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

The larger example tests the complete four-family proper two-loop component
at fixed formal $p$: independent bounded $U_1$, $O$, $U_1 O$ certificates,
the owned mass points on edges 2 and 3, and the four hard subspaces of $k,l$.
Exact cograph attachment carries those identities into the three-loop graph
without distributing its numerator. Separate tests assemble
$(1-T_G)(1-T_Gamma)(1-T_gamma)X$ as eight signed families, certify each outer
operator branch after its complete inner prefix, and compose the identities
by linearity.

The exact hard matrix contains fourteen nontrivial partitions of the carriers
$p,k,l$ and a fixed-momentum block: seven ordinary hard subsets and seven
correlated assignments, including fixed $k-p$ and $l-k$. Each independent hard
block contributes four measure powers. Bare nonzero controls accompany the
superficially divergent strata. A separate open-spinor logarithmic coefficient
certificate uses repeated-pair gamma reduction and exact finite contraction,
without changing production gamma canonicalization or expanding diagnostic
numerators. These large tests belong to the full main-push selection. The
operator identity, local UV limits and integrated renormalization conditions
remain independent validation obligations.

The first milestone consists of the one-loop local ordered OS operation, its
vacuum-only analytic correction, exact CFF reconstruction, independent contour
checks and integrated mass/field references. Local identities use both
fermion-loop and gluon-loop input bases. Analytic checks retain epsilon terms
through quadratic order, use both integration modes, two auxiliary scales and
two analytic backends. The physical Dirac CFF test adds exact rational residue
controls after the general numerator and denominator certificates. Production
generation and evaluation tests exercise local and integrated configurations.

For the three nested fixtures, the independent complete-family jets and
projection certificates are followed by child-only, parent-only, simultaneous
and correlated hard limits. The relevant constant/linear IR soft jets have
separate checks. Child scalar-residue controls distinguish rest-frame OS
points from moving physical on-shell points. Complete radial CFF tests cover
physical child-gluon, rest-frame child-gluon, outer-gluon and simultaneous-gluon
rays for OS inside OS; both child-gluon rays for OS inside IR; and the outer
gluon for IR inside OS. Configuration acceptance alone proves none of these
claims.

The public guard admits OS/OS nesting, an OS parent with an IR child, and
disjoint OS components. An IR parent with an OS child additionally requires
`Graph::has_qcd_self_energy_insertion_ward_identity` to certify the actual
native minimal-QCD vertex and propagator tensors, topology, color flow and
equal masses. An autogenerated source label is insufficient. The negative
fixture changes a rule through the immutable model and must be rejected.
Other mixed prescriptions remain unsupported.

Both direct 3D modes implement the isolated operation and share the exact
source certificate infrastructure for proper components. Their tests compare
physical proper one-loop children, completed two-loop parents in the full
graph domain, finite analytic prefixes and a two-component disjoint product
with the projected route. Public direct dispatch also admits a proper two-loop
OS component with exactly one classified connected proper descendant, which
must be a degree-one one-loop OS component. Bare, completed-child and
finite-prefix tests retain the live cograph in that domain. Other proper
multiloop histories remain guarded.

#table(columns: (1.25fr, 1.15fr, 2fr), inset: 5pt,
  table.header([Component or forest], [Prescription], [Required certificates and guard]),
  [One-loop massive quark self-energy], [OS mass; MS-bar field],
  [Local operator, vacuum scheme correction, exact projected CFF and
  independent local and integrated references.],
  [Massless quark specialization], [OS at $m=0$],
  [Exact agreement with degree-one soft refinement.],
  [Quark child inside gluon parent], [OS child; IR parent],
  [Actual-rule Ward certificate, nested analytic input, complete-family jets,
  exact CFF projection, four hard subspaces, external soft jet and complete
  physical/rest-frame child-gluon controls; arbitrary numerators excluded.],
  [Quark child inside quark parent], [OS child and parent],
  [Complete-family jets, exact nested CFF, four hard subspaces and the complete
  four-stratum radial soft matrix in the stated below-threshold domain.],
  [Gluon child inside quark parent], [IR child; OS parent],
  [Complete-family jets, exact nested CFF, four hard subspaces, gluon-child
  constant/linear soft jet and complete radial CFF soft control.],
  [Disjoint quark self-energies], [OS on each component],
  [Independent component identities, owned product, live bridge, exact CFF,
  left/right/both/correlated hard limits and left/right/joint soft controls.],
  [Proper two-loop quark component in a three-loop graph], [OS at both inner levels],
  [Component/mass-point/cograph certificates, all eight forest families and
  fourteen ordinary/correlated hard strata. Direct dispatch admits exactly
  one primitive OS descendant. Public route comparisons precede the full
  numerical matrix of 26 LMBs and 182 nonempty hard-coordinate subsets.],
  [Ordinary UV or pole parent above OS], [MUV or PolePart],
  [Guarded: a general mixed integrated policy has not been established.],
  [Massive scalar/vector or analytic chiral subgraph], [OS],
  [Unsupported fermionic projector or dimensional numerator prescription.],
  [Direct 3D, primitive and small forests], [Supported OS/IR combinations above],
  [Primitive, proper-child, complete-parent, finite-prefix and two-component
  comparisons in both modes. Other proper multiloop histories, disjoint
  products containing a multiloop component, and unions of at least three
  OS components remain guarded.],
)

The larger graph uses OS for every fermion self-energy. A three-loop IR parent
with a two-loop OS child also has forests that skip the middle component;
the one-loop insertion Ward argument does not establish that mixed policy,
so its public guard remains in force.

=== Original OS validation records

The following timings belong to the original OS stack, not the FeynKit rebase.
They preserve the earlier run evidence while the maintained tests above define
the current obligations:

- `tests/artifacts/os_three_loop_repeated_pairs_validation.json`: proper
  two-loop component 5.78 seconds; eight-family identity 73.89 seconds;
  seven ordinary hard strata 580.22 seconds and seven correlated strata
  26.43 seconds.
- `tests/artifacts/os_logarithmic_coefficient_validation.json`: exact
  open-spinor logarithmic coefficient 241.53 seconds.
- `tests/artifacts/os_full_radial_validation.json`: complete OS/OS radial
  matrix 299.55 seconds and OS-inside-IR child-gluon matrix 722.76 seconds.
- `tests/artifacts/os_proper_guard_validation.json`: IR-inside-OS radial
  control 57.74 seconds. The proper two-loop bare, completed-child and
  finite-prefix route comparisons recorded 72.70, 316.35 and 374.70 seconds.
- `tests/artifacts/os_three_loop_grouped_routes_validation.json`: full
  three-loop generation and summed route comparison 1,773.62 seconds.

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
actual native vertex and propagator tensors, color flow and equal masses are
prerequisites for the specialized nested policy. Admission compares those
tensors exactly through the typed algebra owner; autogenerated metadata alone
cannot establish the Ward identity. The independent all-$d$ bubble and
mass-insertion identities, signed nested analytic input and positive-epsilon
finite child prefix require their own tests. Complete-family Taylor identities,
four hard-subspace checks, the relevant external soft jet and complete CFF
physical/rest-frame child-gluon controls establish separate claims. Local
forest identities alone do not establish the integrated policy.

The logarithmic box component on edges 2, 3, 5, 6 of the figure-10 fixture has
a cograph on edge 4 independent of the outer external momentum. Its child-only
and IR-parent-composed forest families cancel because the parent soft constant
coefficient reproduces that complete insertion. The independent fixture test
requires this cancellation, including equality of the parent's $U$ and
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
kernels. The complete OS/OS CFF sum separately tests the four radial strata
below after all parent residues. OS-inside-IR tests its two child-gluon
strata; IR-inside-OS tests the outer-gluon ray.
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
threshold, and keep the nonsoft loop generic. These radial controls use the default rest-frame
subtraction and do not assert a pointwise condition at a different boosted
point. The matching fixed-boost four-dimensional test is described separately
above. The moving physical child-gluon limit belongs to the complete CFF orientation and
forest sum; its boundedness must not be inferred from one isolated residue.

Mixed-prescription guards can be relaxed only for those cases actually derived
and tested. Required coverage includes all relevant hard subspaces and the
corresponding gluon soft limits, rather than only uniform scaling of every loop.
Exact small forest terms precede the large graph tests.

The direct representations implement the same isolated-component operation;
their complete orientation sums must agree with the projected-four-dimensional
route. Small proper-child, completed-prefix and full mixed-forest checks retain
separate ownership and equivalence certificates. The proper two-loop component
with one primitive OS child requires the same three-route certificates in the
larger graph. After equivalence is established, measurements choose the fastest
supported route for the large scaling matrix; small equivalence checks remain
for the other routes.

The completed scalar evaluator groups conditional summands by their exact
orientation condition at each additive node. It preserves products, powers
and lazy branch boundaries; each changed pass removes at least one conditional,
so the fixed-point iteration terminates. Focused tests require exact symbolic
selector truth tables, same-key sums, factorization and idempotence. The
primitive comparison checks both direct constructions per production host,
then all three routes after summing hosts. These checks establish the grouping
transformation independently of large public generation and route comparisons.

== Rebased FeynKit validation and measurements

=== Historical frozen temporal-trace stack

This section records the original frozen source's qualification. Its runs and
measurements are historical; the current retained-tensor checks above have
their own source identities, selections and evidence.

The frozen reference stack is
`6a5a3c3faf7bb662998a8b402549aea7f0573d19`, named
`on-shell-validated-temporal-stack`. Its shared temporal-trace patch is exactly
`fd1fef9ff2d4d176e5658ad2ef8f74a63c552551`, copied above the other shared
repairs as `7b9380de7b668c2fbec8a884d9a2bd4fb0796d81`. The physical
four-dimensional trace recurrence preserves formal-dimension and disabled-trace
boundaries. The associated evaluator change releases reduction caches before
component graph construction. The shared patch passes all 152 Dirac tests,
11 evaluator controls and strict Idenso/GammaLoop Clippy.

The OS/core gates pass 79 bounded cases and 9 large cases, respectively:
`c09e9887-e3a8-4ecf-9c78-ca94a873f6dc` and
`89ecea37-c8aa-4d49-8906-8fa1004161d9`. The bounded total is 226.480 seconds;
the longest individual case is 62.151 seconds under the unchanged 240-second
limit. Large checks remain in the existing full-main lane. They include the
independent complete three-loop component jets, proper two-loop cographs,
generic and correlated hard limits, open finite coefficients and the maintained
OS/OS and OS/IR radial soft strata. All 485 guarded sources and the core binary
remain unchanged across these runs. The complete record is
`target/os-feynkit/core-validation-after-temporal-trace/`.

All six scalar CLI fixtures pass the three-route comparison at both nonzero
points: the integrated primitive, quark in gluon, nested quark, soft gluon in
quark, disjoint components and the three-loop chain. The three-loop run is
`a46e7cd4-fe58-41a2-a867-be921fcd5642`; the primitive run is
`c319bc7c-d0ba-4b82-9cf4-f7892efd2fdd`; the four two-loop cases are
`a531d84b-ba4e-46e4-9f7b-484d14c62d88`. Numerator and denominator reconstruction
certificates remain mandatory during generation. The final numerical comparison
uses the existing stability-aware complex-norm criterion; it is not a symbolic
identity proof. These shared-process, concurrent correctness runs do not supply
isolated performance or per-route peak-memory measurements.

The small physical checks pass in
`a89026cf-9e23-4b83-bfce-cbf919f90925` (Born-density addbacks, route/sampling
probes and the published reference) and
`d5edb5bc-bb3a-40b0-8cbc-98ad25294ca1` (localization-kernel jets). The complete
nine-case physical acceptance then passes in
`50e78246-d7d6-494c-b55a-0b8257cc95bd`, taking 2952.964 seconds under the
unchanged 7200-second limit. It retains the same 320,000 samples per case,
seeds 1337 through 1345, ten workers, spherical coordinates and all-LMB proposal.
All reference, imaginary-part, auxiliary-scale and cross-route assertions pass.
The maximum reference distance is 1.616 reported standard errors; the maintained
auxiliary-scale and cross-route distances are at most 0.919 and 1.585 combined
standard errors. These are the test's comparisons, without a claim that adjacent
seed offsets produce independent random streams.

Every case has zero unstable samples. Quad precision is used for 10,354 samples;
no arbitrary-precision escalation is used. All 1306 source guards and physical
test binary `46f6c74a116d4c5a05d6098806737f79a865906043e9514a5983621f9ae3fea2`
remain unchanged. This fresh run covers the Quad paths that byte-identical
compiled double-precision objects alone could not certify. The LO and all nine
NLO checkpoints, 31 compiled objects, results and source/binary manifests remain
at their original paths in
`target/os-feynkit/physical-acceptance-after-temporal-trace/`. In particular,
the earlier final-checkpoint gap is closed.

Shared repairs retain their own tests and immutable named commits. Numeric
normalization `6480e4cd` distributes only the current numeric coefficient at
each product/sum boundary. Parser batching `eab3d857` combines exact scalar
summands after tensor admission while preserving the ordered fold in other
coefficient domains. Printer registration `6eda5d7a` serializes first callback
creation across subset/index types, with unchanged formatter output and error
propagation after unlocking. The latter passes 261 Linnet unit tests and 256
concurrent registrations in the unchanged reproducer. Strict production-library
Clippy passes; 39 inherited all-target `energy-timing` configuration warnings
remain separately recorded. No new cyclic-word Taylor repair is claimed here.

Encoded Atom bytes are a factorized size diagnostic, not an algebraic invariant:
pre-evaluator family tags come from a process-wide scope counter, and Symbolica
uses packed symbol identifiers. An isolated quark-in-gluon measurement matches
all ten numerical values and sizes of the older isolated measurement exactly,
while the two sequential-route records share another size pattern. The audit
is in `target/os-feynkit/performance-after-registration-fix/qg-history-review/`.
The measurement tooling preserves both sizes and literal-equality diagnostics,
checks graph/count inventory, and uses the strict machine-epsilon floor of the
unchanged numerical criterion. It does not infer algebraic equivalence from
matching sizes or a numerical sample. Graph numerators remain factorized.

The following earlier records retain their original revision and scope. Five
isolated measurements at `91ab7745` are historical observations; they are not
relabelled as temporal-trace or consolidated-rebase results. At the time of
this historical record, the exhaustive UV inventory of 259 limits per control
and 39 isolated replays was prepared with its existing reference, seed,
precision, timeout and performance-gate policies. Final measurements awaited
the consolidated FeynKit repairs, including component factor collection and
the restored antighost derivative. No full three-loop radial soft matrix was
supplied by those two-loop checks.

=== Earlier rebased records

At source revision `1121da7ae5f2f68312bff6e4341175de084e1f45`, the fresh native
$gamma^* arrow.r t overline(t)$ capture has the same complete topology inventory
as before the paired-wavefunction sampling fix, identically in all three
routes. GL0 is the vertex with four physical cuts; GL1 and GL2 are self-energy
masters with three cuts each. The conservative numerator comparison retains
both self-energy representatives. The acceptance classifies these as two
physical roles, retaining every native master, group member and physical cut.
MS-bar reference rows select all masters of their role through the existing
normalized runtime subset owner; no manual multiplicity is introduced. OS
acceptance retains the fully inclusive graph sum. The capture and comparison
are in `target/os-feynkit/physical-topology-after-sampling/topology-summary.json`.

Eight small public checks passed in 119.452 seconds, run
`17bb38eb-c91a-48f3-b056-2a44382b4ea6`: the independent published-reference
anchors, primitive projected dispatch, rest and boosted integrated three-route
generation, boosted state reload, primitive orientation UV controls, and both
appendix-B1 dispatch checks. The all-master physical x-space probe passed in
17.674 seconds, run `ee74544c-c29a-4de7-a06a-499e224ff933`. It compares the
three routes at three fixed points and checks that all native master
contributions sum to the complete weighted value. These numerical boundary
checks supplement the exact certificates; they do not establish the inclusive
integrated OS reference or auxiliary-scale independence.

The native component/CFF and physical localization jet regressions are
preserved in change `1bcacabc9743c5c2bc07b1899b2d2ef723ba82da`. They passed
together in run `dc8847e1-cb83-489f-87ac-2a04e8fb4a8a` (two cases, 3.734
seconds), against the same production code. Their source hashes begin
`e42c2c3c4141` and `2d796c3ad592`. The component test
certifies the symbolic-dimensional numerator and both rest/nonaxial signed
OS and vacuum expressions, then the analytic field pole and finite addback.
It uses the exact boundary signatures to distinguish a point's edge owner
from its momentum coordinate. The adjacent native propagator numerators
satisfy $N(p) F(p) N(p)=i(p^2-m^2)N(p)$ with their actual spin/color interfaces;
the additional integrated coefficient supplies the second factor of $i$.
This numerator identity does not cancel a denominator's nonzero causal
prescription as an ordinary rational function.

Both the raised native cograph and its independent pinched source pass general
numerator and denominator reconstruction certificates in a compatible parent
chart. Pinching removes one original basis carrier, so the replacement chart
is certified against every active signed source signature before either
quotient is constructed. Each native graph produces eight raised source maps
and two pinched maps. Their complete CFF orientation sums agree at two exact
rational component points, giving four exact specializations. These supplement
the general source certificates; they are not a proof of a general symbolic
identity between the final CFF sums. Core Clippy passes with warnings denied.

The separate runtime test executes the existing parameter binding, dual loop
rescaling and LU measure on a small raised-cut scalar cograph. At two values
of $t$ and four independent UV/localization scale pairs it checks the zeroth
and first jets of $t^3 rho_L(t q)$, including the kernel derivative and the
measure derivative. The UV and localization scales remain constant under this
physical rescaling. This small runtime check does not certify the complete
photon graph's integrated rate.

The six isolated performance tests passed in 221.243 seconds, run
`19d2c848-1924-4bd9-beba-79e9ce68c1c5`, under the ordinary 240-second per-test
limit. Each route has a distinct process, one generation measurement and ten
evaluation repeats: two nonsingular momentum points, five repeats each, twenty
batches per repeat. The primitive scope includes local and integrated terms;
the quark-in-gluon scope is local only. The table reports generation wall time
without stage logging and evaluator construction in seconds, followed by the
mean and sample standard deviation of the ten repeat means in µs per sample.
The total includes the benchmark wrapper; the evaluator column measures its
kernel category. Compiled evaluator construction is zero in these runs.

#table(
  columns: (1.45fr, 0.8fr, 0.9fr, 1.2fr, 1.2fr),
  table.header([Fixture / route], [Generate (s)], [Evaluator build (s)],
    [Evaluator (µs)], [Total (µs)]),
  [Primitive / projected], [0.243], [0.0222], [1.500 ± 0.005], [28.661 ± 0.259],
  [Primitive / explicit], [0.236], [0.0255], [1.338 ± 0.010], [28.485 ± 0.617],
  [Primitive / orientation], [0.243], [0.0310], [2.511 ± 0.103], [34.439 ± 2.814],
  [Quark in gluon / projected], [5.490], [1.838], [42.277 ± 0.446], [83.242 ± 5.819],
  [Quark in gluon / explicit], [98.389], [93.481], [91.279 ± 0.919], [130.361 ± 1.368],
  [Quark in gluon / orientation], [96.151], [92.512], [222.738 ± 2.489], [276.861 ± 8.350],
)

Generation RSS is sampled process memory, not incremental allocations or the
whole test's lifetime peak. Separate Linux `VmHWM` readings, which include
model import, remain in the artifacts. Stored sizes count the factorized
pre-evaluator root and retained definitions without unfolding them.

#table(
  columns: (1.45fr, 0.85fr, 1fr, 1.4fr),
  table.header([Fixture / route], [Sampled RSS (MiB)], [Root bytes],
    [Definition bytes / count]),
  [Primitive / projected], [55.32], [1,554], [6,733 / 7],
  [Primitive / explicit], [56.34], [4,025], [6,493 / 11],
  [Primitive / orientation], [56.48], [4,165], [6,493 / 11],
  [Quark in gluon / projected], [283.00], [180,809], [76,473 / 28],
  [Quark in gluon / explicit], [2,681.79], [1,005,151], [94,145 / 75],
  [Quark in gluon / orientation], [2,693.15], [1,094,096], [94,145 / 75],
)

The baseline has no stage logfile. Separate generation-only replays enable
module-scoped profiling and reproduce the same factorized expression sizes.
Their times below are not substituted into the baseline wall measurements.
The forest timer includes traversal, completed local operations, analytic
integration when enabled, projection and assembly. Its remainder subtracts
only the profiled CFF work inside that boundary; it still includes surface
conversion and unprofiled CFF-related work, so it is not pure Taylor time.
CFF totals include routing, source preparation and assignment selection,
including rejected proposals once. Native engine time is a subset. The
production root CFF lies outside the forest timer. These columns are
overlapping measurements and must not be added to the inclusive forest time.

#table(
  columns: (1.45fr, 1.6fr, 1.4fr),
  table.header([Fixture / route], [Forest inclusive / remainder (ms)],
    [CFF total / native (ms)]),
  [Primitive / projected], [193.078 / 175.333], [18.045 / 0.902],
  [Primitive / explicit], [185.022 / 184.832], [0.534 / 0.292],
  [Primitive / orientation], [179.146 / 178.973], [0.462 / 0.254],
  [Quark in gluon / projected], [3,618.258 / 2,939.564], [680.060 / 47.461],
  [Quark in gluon / explicit], [4,884.653 / 4,873.615], [12.418 / 1.727],
  [Quark in gluon / orientation], [3,676.512 / 3,665.832], [12.034 / 1.700],
)

The extraction requires one matched forest interval, matching outer-routing
boundaries, and one numerator-template record per admitted exact assignment.
Nested preparation and raw proposals are not charged again within their
enclosing records; negative remainders are rejected. Component-operation and
inclusive projection timers remain separate diagnostics in the artifacts.

The requested two seconds per repeat calibrates a fixed workload from ten
warmup samples; it is not a minimum measured duration. The six baseline cases
execute 1,833–4,667 samples per repeat, with actual batch-wall totals of about
0.12–0.78 seconds per repeat. Repeat spread is not independent-run uncertainty,
and a single generation has no measured generation-time spread. Source hashes,
per-point results, exact counts and stage records are retained in
`target/os-feynkit/performance-review/six-route-summary.json` and
`complete-phase-summary.json` in the same directory. The accounting recipe is
`complete-phase-methodology.typ` there.

A separate three-process primitive replay disables integrated terms to match
the UV scaling workload. Its evaluator and total columns retain the same
repeat-mean convention; construction times are milliseconds. The final column
gives sampled generation RSS and stored root-plus-definition bytes. Its
standalone-process RSS baseline must not be compared directly with the earlier
integration-test executable's absolute RSS.

#table(
  columns: (1.15fr, 1fr, 1.15fr, 1.15fr, 1.2fr),
  table.header([Local-only route], [Generate / build (ms)], [Evaluator (µs)],
    [Total (µs)], [RSS (MiB) / bytes]),
  [Projected 4D], [150.54 / 20.19], [1.263 ± 0.008], [29.336 ± 2.264], [75.61 / 7,663],
  [Explicit direct], [125.82 / 24.42], [1.060 ± 0.033], [31.395 ± 2.728], [81.79 / 9,895],
  [Orientation direct], [127.78 / 28.58], [2.377 ± 0.007], [32.285 ± 2.575], [80.45 / 10,043],
)

This matched replay supports retaining explicit direct for the bounded primitive
scan: its observed generation is smallest and its evaluator is faster at both
points. Total wrapper means overlap in repeat spread and favour projected 4D
in their point estimates. The quark-in-gluon measurements strongly favour
projected 4D for construction, memory and evaluation. Neither comparison
establishes a universal fastest route for the other forests. The local-only
record is `target/os-feynkit/performance-review/primitive-local-only-summary.json`.

The inclusive run `be16f780-3c24-488c-a161-520ed7b59c75` completed all nine
320,000-sample integrations in 597.153 seconds, then failed the unchanged
reference assertion for orientation-direct at doubled UV mass:
$0.42349060 plus.minus 0.00294377$ versus $0.43283009$ (3.17 standard errors).
An offline check also finds the orientation-direct localization case 3.05
standard errors below the reference. These are failed validation results,
not evidence of an algebraic mismatch. The final assertion loop stopped at
the first miss; unexecuted scale and route assertions are not counted as passes.
The existing integrator initializes iteration $j$ with seed $s+j$, so the
consecutive requested seeds also reuse RNG streams across runs. Their claimed
statistical independence is therefore not established.

The retained-extrema control evaluates eight original inputs in all
three routes, using production precision and forced arbitrary precision.
All 48 cold versus full-warm values agree exactly; the largest relative
between-route difference in forced arbitrary precision is $2.02 times 10^(-296)$.
The largest Double-versus-Arb difference is $3.05 times 10^(-11)$. These are
unit-weight API values including their physical parameterization Jacobian;
the original Monte Carlo extrema also include a Vegas weight and are not
equated to these values. The control excludes a percent-level route or warm-up
discrepancy at the tested inputs, but neither proves a general identity nor
establishes the cause of the failed integrations.

The original Vegas-weight boundary was checked separately against all eight
retained extrema. Each stored extremum is reproduced bitwise using the native
physical value, its accumulation precision and its original sample weight.
All 96 weighted replay comparisons satisfy the requested Double tolerance.
No weight from the final adapted grid replaces an original sample's weight.

Fresh 1,280,000-sample diagnostics retain the two failing orientation-direct
seeds and settings, increasing only the diagnostic sample budget. The doubled
UV-mass result is $0.42520132 plus.minus 0.00167351$ (4.56 standard errors below
the reference); the doubled-localization result is
$0.42484646 plus.minus 0.00240693$ (3.32 standard errors below). Their first
300,000 samples share the original trajectories, so these are not independent
replications. The committed acceptance bounds, seeds and budgets are unchanged.

A separate probe uses 20,000 predetermined common samples from the original
baseline's frozen grid, on a distinct native RNG stream. Subtracting settings
per sample gives a doubled-UV-mass difference of
$+0.00095338 plus.minus 0.00375058$, and a doubled-localization difference of
$-0.00653574 plus.minus 0.00092152$. Seventeen of twenty fixed blocks have
negative localization differences. Its
largest contribution supplies 34.7% of the centered second moment and its ten
largest supply 64.6%. Thus the reported localization discrepancy is a
statistical observation, not an exact normalization defect. Every sample is
retained; none is clipped or replaced. The complete records are under
`target/os-feynkit/physical-mc-triage/`, in `statistics-fresh/`,
`vegas-weight-replay/` and `paired-scale-corrected-run/`.

Subsequent checks identify an infinite-variance tail in that sampling setup.
For a gluon/spectator channel $(k,p)$, relative spherical coordinates sample
$(k+p)/2$ and $k-p$ independently. At fixed spectator momentum and large
$abs(k)=R$, the linear radial maps give proposal density $O(R^(-8))$, whereas
the localized field density falls as $R^(-4)$. The native self-energy graphs
each have one such channel. The vertex graph's alternative channel has
density $O(R^(-4))$ in the same physical region; its map-density partition
suppresses the poorly matched contribution and cannot cancel the
self-energy estimator's leading $R^4$ growth. Its nonzero coefficient follows
from the field/Born identity above. The resulting second moment contains
$integral^infinity R^2 dif R$ and diverges. A finite Vegas histogram changes
the prefactor, not that asymptotic power.

Seven complete x-space probes confirm the noncancellation and the approach to
the predicted fourth power, with all graph/channel contributions retained.
Consequently, the quoted standard-error multiples from this proposal are not
calibrated Gaussian significances. Increasing its sample count does not repair
that variance defect. The derivation and native maps are retained in
`sampling-tail-derivation.typ` and `localization-complete-tail/` under the
triage directory. The recorded baseline and doubled-localization complete
values also approach the fourth power on this ray; other individual settings
are not certified by that observation.

The replacement uses ordinary spherical coordinates with the existing
`auto:lmb` catalogue: eight native channels for GL0 and seven for each
self-energy. In 122 controlled evaluations, both fermion and gluon soft
directions have sufficient proposal density, and the hard-child body and
proposal fall as $q^(-4)$. The complete baseline, doubled-localization and
contrast weights approach constants on the UV sequence. These checks cover
the identified regions, not every possible singular region. The retained
record is `spherical-atlas-generic-rays/` under the same triage directory.

The finite physics coefficient was checked independently of that sampler.
At three momenta and three UV/localization scale pairs, all 27 differences
between integrated-enabled and integrated-disabled generation agree with the
field and vertex Born-density references above. The largest relative residual
is $1.32 times 10^(-298)$ using the stored model coupling; the separate
$alpha_s$ reference differs by $3.79 times 10^(-16)$ from that coupling's stored
rounding. These are pointwise runtime controls, alongside the exact native
field identity, not a numerical proof of an integrated equality. Their record
is `finite-addback-reference/` in the same triage directory.

The maintained finite-density test, the three-route sampler comparison and
the published-reference anchors passed in run
`862d2f07-5a4c-451e-80f0-3dbd6d2e91c2` at source
`9768730e262343a55e74c0ddc66db99f558c8687`. Total runtime was 46.118 seconds;
the two generation cases took 23.317 and 22.786 seconds. The sampler comparison
retains both old and replacement configurations, every native master and the
complete graph sum. It also requires a full-rank native channel exposing each
OS child edge with an independent cograph carrier. Recorded Arb diagnostic
maps use exact binary64 input promotion and are distinguished from the
evaluator's fixed proposal precision. Integration-test Clippy and the
documentation registry check pass.

The maintained inclusive acceptance then passed all nine route/scale cases in
run `4ff2a464-957b-4529-b6e6-05026b2b389f`, with the replacement proposal and
unchanged references, tolerances, seeds, learning rates and 320,000 samples per
case. The launch revision record identifies `ee8e27131e80`; the acceptance
source SHA-256 is
`f6430735a93065d5d5811a9c2eef1d47a45df6193829cbdf41b07bb3e35cb485`.
The retained binary and full run manifest identify the executable independently
of later documentation edits. At the analytic reference $0.432830086776986$:

#table(
  columns: (1.1fr, 1fr, 1fr, 1fr),
  table.header([Route], [$(M,L)=(173,173)$], [$(346,173)$], [$(173,346)$]),
  [Projected 4D], [$0.432276 plus.minus 0.002624$],
    [$0.433514 plus.minus 0.002703$], [$0.434514 plus.minus 0.002675$],
  [Explicit direct], [$0.436554 plus.minus 0.002623$],
    [$0.433403 plus.minus 0.002702$], [$0.433117 plus.minus 0.002666$],
  [Orientation direct], [$0.430494 plus.minus 0.002619$],
    [$0.430878 plus.minus 0.002695$], [$0.428544 plus.minus 0.002652$],
)

Every reference, imaginary-part, auxiliary-scale and cross-route assertion
passed. The largest reference distance is 1.616 reported standard errors;
the largest maintained combined-error comparison is 0.919 for auxiliary scales
and 1.585 against the projected route. These reproduce the test's existing
comparisons; the adjacent seed offsets still do not establish independent
RNG streams. The LO check also passed. The complete test took 3029.939 seconds
under the existing full-main deadline; concurrent validation workloads mean this is
validation timing, not a new performance benchmark.

The exact results, JUnit, tested source and binary, and hash manifests are in
`target/os-feynkit/spherical-physical-acceptance-run/`. The normal successful
test cleanup removed its original state directory. Review copies retain the LO
and first eight NLO checkpoints, all 31 compiled shared objects, and the on-disk
model state; the final case's central result survives but its checkpoint was
removed before copying. The generated process collection was not serialized,
so these copies are not a self-contained resumable state. Replay requires the
maintained generation inputs and original paths, with fingerprint guards intact.

This pass establishes the inclusive reference and auxiliary-scale checks for
this rest-frame physical calculation. It does not turn the numerical proposal
ray checks into a proof of global finite variance, or replace exact algebraic
reconstruction and local-limit certificates. The earlier failed proposal and
all diagnostic samples remain preserved rather than being relabelled as passes.

A separate nine-case large-core run,
`3a4a2f9d-44ee-4374-abf3-b368806c4f9f`, passed eight cases in 1323.896 seconds:
both outstanding OS-in-OS and OS-in-IR radial soft controls, all three proper
two-loop live-cograph comparisons, the generic and correlated three-loop hard
limits, and the open-spinor finite-coefficient certificate. The independent
three-loop component-jet comparison initially failed after 193.880 seconds,
at class 16 of the final nested outer operation. Logs and unchanged-source
identities are in `target/os-feynkit/large-algebra-soft-validation/`.

Each of the twelve captured finite entries contains the following zero factor:
$ (a+b)(c+d)+(a+b)(c-d)+(a-b)(c+d)+(a-b)(c-d)-4 a c = 0 $
Pairwise factor collection proves this identity while the graph spectator
remains untouched. The test certificate had selected the related signed
binomials as independent opaque keys, hiding their relation. Shared commit
`3d166ca1105a7e4d32980ec1b3ade67710da72d2` changes only that test-certificate
selection: distinct summed keys with equal coefficient-free term support
leave their atomic constituents available to the existing bounded exact
collector. Positive-power bases remain protected, as do function bodies and
inverse denominators. No production operator or normalizer changes.

The native momentum-component regression fails with the old certificate and
passes with the correction, including nonzero and retained-factor controls.
All fifteen certificate tests pass in canonical run
`c2c26211-903c-438a-9095-dd8f55b08a0e` (0.945 seconds), and all twelve captured
entries independently replay to zero. The complete original three-loop
assertion then passes unchanged in run
`7324cfbb-3d2f-4961-bbdf-debb764cc1b1` (194.950 seconds), with its input hashes
unchanged during the run. These are exact certificates, separate from the
hard-limit and physical integration results. The full rerun is retained in
`target/os-feynkit/three-loop-identity-fixed/`.
All twenty-six bounded consumers then pass in run
`3f780ba2-8fb7-42d0-b317-d2a13607f767`: one-loop and native-component identities,
small nested forests, integrated boosted poles, complete rainbow soft limits,
primitive direct-route controls and overlap-proof guards. The serial total is
225.193 seconds; the longest case takes 61.710 seconds under the unchanged
240-second per-case limit. Strict Clippy and workspace formatting also pass.

The scalar CLI fixtures for a soft gluon inside a quark self-energy, disjoint
quark self-energies and the three-loop chain attach native external fermion
states at final integrand assembly. Their component numerators remain open
through every counterterm operation. An explicit unit projector in the old
fixtures suppressed that final contraction and correctly failed the scalar
evaluator's interface check. Removing only that override leaves the complete
factorized component numerators, every edge and vertex factor, denominator
atoms, particles, routing and signed loop bases exactly unchanged, as checked
in `target/os-feynkit/projector-fixture-ownership/`.

The disjoint fixture keeps $p=(250,0,0,50)$ so its live bridge remains off its
physical pole. Its external spinors define a scalar diagnostic functional of
this off-shell Green function, with the supplied timelike momentum determining
their components; this is not a physical on-shell top amplitude. Open-matrix
component certificates and the integrated on-shell conditions remain separate
checks. No momentum, helicity, scale, seed or numerical tolerance changes.

The remaining CLI route comparisons and their extensive numerical UV matrices
are still open. No full-three-loop radial soft matrix is supplied by the
two-loop controls. The following historical
results do not discharge those obligations for this revision. Performance-gate
triggers and the ordinary four-minute versus full-main test policy remain
unchanged.

== Original OS baseline measurements and route selection

All numerical results, timings, memory figures and route choices in this
section are historical observations from the original OS stack ending at
`ce58d97de847`; the binary and artifact identities below identify its individual
runs. They are retained as a baseline, not as validation or performance data
for the rebased FeynKit stack. New measurements must identify the new revision
and repeat the same boundaries and protocols before updating a route choice.

The isolated `ci-optim` measurements from 2026-10-13 use one fresh process per
fixture and route, binary SHA-256 prefix `c46b91da6603`, with orientation branch
grouping and generation profiling enabled. All six cases passed. The primitive
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

The three other local-only forest matrices also passed with binary
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

These measurements favored explicit direct for IR inside OS and projected
4D for the disjoint pair. For
OS inside OS, the evaluator means are close: projected 4D is 0.5–1.4 percent
lower, within the measured variation. Its generation takes 3.49 seconds
versus 7.28 for explicit direct, and its total sample means are lower.
Projected 4D was therefore the provisional choice for that extensive
matrix. These choices retain small equivalence checks for the other routes
and make no isolated-memory claim from the shared-process matrices.

The full three-loop local forest also passed public generation and
orientation-summed pointwise agreement across all three routes with binary
`c46b91da6603`. Its complete test took 1,773.62 seconds. The profile and
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

The exhaustive numerical UV matrix passed on projected 4D: 26 physical LMBs
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
any limit, scale or precision retry. This bounded-resource run passed in
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

All three routes passed the integrated reference and auxiliary-scale checks at
$(M,L)=(173,173),(346,173),(173,346)$. Each case uses 320,000 samples with
independent seeds 1337–1345, relative-spherical coordinates, map-density channel
weights, and explicit sums over graphs, orientations and sampling channels.
The Vegas continuous-dimension learning rate is 0.25; batches start at 20,000
samples and increase by 20,000, with 20,000 samples required before a grid
update. Ten workers evaluate each case. The complete acceptance took 486.70
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

Those runs also agreed in arbitrary precision at momentum-space controls for
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
