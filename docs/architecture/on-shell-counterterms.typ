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

Current qualification as of 20 October 2026 is summarized below and detailed
in #link(<latest-feynkit-replay-qualification>)[Latest FeynKit replay and
current-source qualification]. The combined source has its own closed
forty-case qualification, followed by the vector-reference and independent
signed-center qualifications below. Earlier source-79 results and the
separately bound numerical reference retain their original source and run
identities.

#table(columns: (1.25fr, 0.9fr, 2.85fr),
  table.header([Source or claim], [Closed result], [Boundary and remaining work]),
  [Earlier OS `79b315b0`], [69 selected cases],
    [Planned 48, five route matrices, six mandatory soft and ten full-main
     core cases pass; all retain this exact frozen source. Exact, sampled-limit
     and integrated
     comparisons retain their separate evidence.],
  [Endpoint proof `8588b9af`], [6 cases],
    [Twelve exact ordered left-projector endpoint zeros, including both signs
     and rest/nonaxial points; raw bare source is a nonzero control.],
  [Shared memo `d2375ca5`], [13 cases],
    [Exact catalogue and hyperdual controls; bounded unit requests 64,
     cache hits 62, unique inputs 2. No OS speedup or peak-RAM claim.],
  [Combined `c8eade36`], [40 fresh selected cases],
    [Includes all 259 matched UV limits and the 48-group, 240-point
     hierarchical soft scan with 96 tail inequalities.],
  [Direct vector `c9532785`], [5 cases + one integral],
    [Exact source/chart admission and a separate numerical bare-vector
     reference pass; no simultaneous-CFF integration claim.],
  [Vector union `b762b717`], [7 fresh selected cases],
    [Combines the qualified memo, native endpoints and direct vector source;
     exact numerical input equality binds the preceding integral.],
  [Maintained reference `edefbb2c`], [8 fresh selected cases],
    [Seven ordinary cases and the existing main-only numerical reference:
     nine mass plus one vector integral, with 47 mass and four vector checks.],
  [Signed-center union `5b6a297a`], [9 fresh selected cases],
    [Eight ordinary cases include the exact coordinate certificates; the
     full-main scan covers 104 groups, 520 points and 208 tail inequalities.
     Exact export equality binds the maintained numerical reference.],
  [Performance qualification], [Six ordinary cases pass],
    [Both nested direct routes meet the unchanged 240-second cap on
     `c8eade36`. Quiet paired measurements remain separate.],
)

== Historical retained-tensor continuation records

The following early checkpoints retain their original source and run
identities. Their pending-work statements describe those checkpoints;
the current status is the table above and the explicitly linked latest section.

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
At that checkpoint the nine-case inclusive OS workload remained deferred to
the full-main gate; those ordinary checks did not qualify that workload.

At retained-tensor source `b4c11e911594acd1d3a91a992efc0a156bbc9055`, the
three open-component controls passed in run `36c16895-4b5e-4d12-8fb3-52d0e6806127`
(44.174 seconds). Exact old/new fixture imports also retain every component
tensor and logical interface, local factor, denominator and signed basis;
only native external fermion attachment closes the final scalar interface.
The soft-gluon scalar route case passed in run
`a198f3ef-a685-4544-8454-a97a663e9a56` (70.441 seconds), comparing the three
local-only routes at two points under the existing full-main 7200-second cap.
At that checkpoint the other five affected scalar-route and complete-UV
cases remained deferred. Its source, projector identities and reports are bound under
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

== Independent full-dimensional two-loop mass reference

The independent reference below uses the original native export at `289611be`.
Its numerical results are historical; current retained-tensor source and runtime
qualification are recorded separately in this note. It does not certify the production two-loop OS
forest or the emitted four-dimensional local-plus-scheme result. Those remain
separate comparisons with the proper and overall counterterms defined below.

The input has one massive top flavor, five explicitly massless light flavors,
zero widths and zero chemical potentials. Ungrouped, one-particle-irreducible
$t arrow.r t$ generation at QCD order four gives twelve sources: ghost and
gluon bubbles (D00–D01), a three-gluon diagram (D02), a four-gluon tadpole
(D03), five light-quark bubbles (D04–D08), the rainbow (D09), the massive-quark
bubble (D10), and crossed exchange (D11). The five
light sources are retained individually. Their complete scalar, denominator,
phase and weight inputs agree exactly, allowing one integral to be reused
five times. D11 has the exact color factor $C_F(C_F-C_A/2)$.
D03 has a nonzero numerator and a separate scaleless-integral certificate:
its isolated massless loop is absent from the other propagators and numerator,
and the native parametric scaling satisfies the exact Euler identity.
It must still be retained when constructing local UV and soft forests.

=== Normalization and independent insertions

Write $a=alpha_s/(4pi)$, $d=4-2epsilon$, $m_0=Z_m m$ and
$psi_0=sqrt(Z_psi)psi$, with the full-theory MS-bar coupling and
$n_f=n_l+n_h$. Use the same loop measure as above and the mass projection

$ Phi[X]=frac("tr"_("spin,color")[(slash(p)+m)X],4m^2 N_c),
  quad p^2=m^2. $

The spin trace has $"tr"(1)=4$ while Lorentz contractions retain formal $d$.
The generated ordered amplitude contains a named external-fermion Wick sign;
remove only that factor to obtain the amputated kernel $Gamma_2=-i Sigma$.
Retain closed-loop signs, symmetry factors, numerator prefactors and the exact
model couplings, including $"GC"_(11)=i g_s$. In an $ell$-loop PySecDec input the
resulting scalar prefactor is $i^(ell+1)$ times the kernel weight and mass
projection, divided by $g_s^(2ell)$. The backend retains its own propagator
sign and Gamma normalization; the MS-bar factor $exp(ell gamma_("E")epsilon)$
is inserted once at $mu=m=1$. Native numerator/denominator reconstruction and
color decomposition are exact; finite tensor contraction does not expand a
whole forest numerator.

Let $F=C_F f$ be the one-loop mass projection divided by $a$, and let
$H=C_F g$ be $m$ times its derivative with respect to the internal physical
mass, holding the external momentum and projector fixed. Give that internal mass a separate symbol
$m_("int")$, keep external $p$ and the projector mass $m$ fixed, and evaluate
the derivative at $m_("int")=m$. The fixed external relation $p^2=m^2$ may
already be imposed; it must not vary with $m_("int")$. At $mu=m$ the independent
anchors are

$ f=3/epsilon+4+epsilon(8+pi^2/4)+O(epsilon^2), quad
  g=6/epsilon+2+epsilon(8+pi^2/2)+O(epsilon^2). $

Exactly, $g=2(1-epsilon)f$. Differentiating the already on-shell-composed answer
would instead give the wrong coefficient. Direct PySecDec anchors and the
native one-loop export agree through $epsilon^1$, including the coupling and
phase normalization. Both positive-epsilon terms are needed when multiplied
by an enclosing simple pole.

Define $Z_m=1+a z_1+a^2 z_2$, $Z_psi=1+a w_1+a^2 w_2$,
$z_1=-F$, $w_1=-C_F/epsilon$ in Feynman gauge, and
$beta_0=11C_A/3-4T_F(n_l+n_h)/3$. If $B=Phi[Sigma_("bare",2)]/a^2$,
then the bare-parameter relation is

$ z_2=-B-z_1 H+beta_0 F/epsilon. $

The reference is independently obtained from the two-loop OS constant in
#link("https://arxiv.org/pdf/hep-ph/0005131")[Melnikov and van Ritbergen,
hep-ph/0005131, equations 39–41], with the coupling converted to $a$ above.
The finite mass relation in
#link("https://arxiv.org/pdf/hep-ph/9912391")[Melnikov and van Ritbergen,
hep-ph/9912391, equations 10–11] supplies a separate cross-check.
At $mu=m$, writing
$B=C_F(C_F B_F+C_A B_A+T_F n_l B_l+T_F n_h B_h)$ gives:

#table(columns: (0.7fr, 0.7fr, 0.7fr, 3fr), inset: 5pt,
  table.header([Channel], [$epsilon^(-2)$], [$epsilon^(-1)$], [Finite coefficient]),
  [$B_F$], [$27/2$], [$75/4$],
    [$441/8+29pi^2/4+12zeta(3)-8pi^2 log(2)$],
  [$B_A$], [$11/2$], [$91/4$],
    [$605/8-5pi^2/12-6zeta(3)+4pi^2 log(2)$],
  [$B_l$], [$-2$], [$-7$], [$-45/2-5pi^2/3$],
  [$B_h$], [$-2$], [$-7$], [$-69/2+7pi^2/3$],
)

=== Historical numerical result and proof scope

Seven independent nonzero two-loop scalar kernels suffice after the exact
light-flavor, D11 color and D03 zero certificates. PySecDec 1.6.6 uses the
on-shell relation before sector decomposition, seed 1337, requested relative
accuracy $10^(-7)$, `minn=10000` and `maxeval=1000000`. Generation uses two
processes with two FORM threads each; C++ compilation and numerical evaluation
each use two workers. No numerical budget, seed or comparison tolerance was
tuned after seeing a result. All Laurent
coefficients through the finite term satisfy the fixed comparison bound
$10^(-6)(1+abs("reference"))$.

#table(columns: (0.65fr, 1.65fr, 1.4fr), inset: 5pt,
  table.header([Channel], [Numerical finite coefficient], [Independent target]),
  [$B_F$], [$86.37560907 plus.minus 0.00000185$], [$86.37560703$],
  [$B_A$], [$91.66467772 plus.minus 0.00000198$], [$91.66467727$],
  [$B_l$], [$-38.94934099 plus.minus 0.00000061$], [$-38.94934067$],
  [$B_h$], [$-11.47092302 plus.minus 0.00000014$], [$-11.47092306$],
)

The displayed errors are linear sums of absolute backend error estimates for
the independent scalar kernels; they are not Gaussian confidence intervals.
No imaginary coefficient survives. At $N_c=3$, $n_l=5$, $n_h=1$, the bare
reference is
$B=38/epsilon^2+289/(3epsilon)+382.7369262019$.
The numerical result is
$37.99999999/epsilon^2+96.33333716/epsilon+382.73693057$.
Adding the independently known one-loop mass and coupling insertions reconstructs
$z_2=22/epsilon^2-17/(3epsilon)-90.1810702379$;
the corresponding numerical finite value is $-90.18107461$ with summed backend
error $0.00001248$. This reconstruction is the same bare comparison expressed
as a renormalization constant, not a second independent two-loop integration.

The exact native source certificates and numerical results are retained under
`target/os-feynkit/physics-validation/scalar-power-refresh/export/` and
`target/os-feynkit/physics-validation/pysecdec/` in the retained
`feynkit-algebra-regressions` task workspace.
`full-two-loop-reference-check.json` records the twelve-source incidence,
input identities, coefficient-by-coefficient comparisons and the explicit
false flag for a production-forest integrated-condition claim. Numerical
agreement does not establish an algebraic identity or promote a later OS
implementation whose own checks are incomplete.

=== Production forest target and the four-dimensional boundary

An actual forest comparison must keep the proper insertions separate from the
overall OS counterterm. The Lagrangian mass vertex contains $Z_psi Z_m-1$,
so its two-loop on-shell mass coefficient is $K_2=z_2+w_1z_1$.
The complete proper one-loop insertions have the independent target

$ P_2=z_1 H-beta_0 F/epsilon+w_1 F, quad B+P_2+K_2=0. $

Thus extract the parameter constant as $z_2=K_2-w_1z_1$.
At $mu=m$, the mixed term is
$w_1z_1=C_F^2(3/epsilon^2+4/epsilon+8+pi^2/4)$.
Omitting it fails already in the double pole. Compare the actual proper
prefixes with $P_2$ before applying the overall OS operator; do not define an
expected counterterm as minus the production bare-plus-forest output.
The QCD logarithmic vertex children remain MUV. OS-over-MUV admission requires
the concrete canonical-vertex, completed-child and analytic-policy certificates;
the full bare reference does not authorize relaxing the general mixed guard.

The classified twelve-source inventory confirms two overlapping logarithmic
vertex children in each of D02 and D11, giving six forest operations per source.
The other ten sources have four operations each. The two overlapping children
never occur together in a proper forest. D02 and D11 are the only rejected
prescription combinations in that inventory; the scaleless bare D03 source
remains present before all forest operations.

The actual four isolated MUV children independently pass their open-tree
Laurent comparison through $epsilon^1$. With the native tree phase $i g_s$ and
$S_M=exp(gamma_("E")epsilon)Gamma(1+epsilon)(mu^2/M^2)^epsilon$, the unsigned
Taylor-integral coefficients relative to the tree are

$ V_A=(C_F-C_A/2) S_M (1-epsilon)^2/epsilon, quad
  V_N=C_A S_M (3-2epsilon)/(2epsilon). $

The first is the two-quark/one-gluon triangle in D11; the second is the
one-quark/two-gluon triangle in D02. These follow from the dimensional vacuum
tensor integral of the actual hard-loop numerator. Angular averaging is not a
pointwise substitution in the unintegrated cograph. At $mu=M$, their series are

$ V_A=(C_F-C_A/2)[1/epsilon-2+epsilon(1+pi^2/12)]+O(epsilon^2), quad
  V_N=C_A[3/(2epsilon)-1+epsilon pi^2/8]+O(epsilon^2). $

The tests compare the signed integral of the actual stored $-U_0$ with
$-g_s^2/(16pi^2)$ times the native tree and these coefficients at $N_c=3$.
The source vertices, their rule IDs and field order, the crown attachment, and
the exact child-times-cograph numerator are certified before integration.
The vacuum integration retains three Laurent terms; both the pole and finite
addback owners enter the comparison. Both historical test cases pass in run
`e0c37b1d-69dd-42c4-ad32-12ff251ecff4` (2.330 seconds), covering two distinct
attachments of each triangle. Their records are retained under
`target/os-feynkit/physics-validation/native-vertex-and-scoped-oracle/isolated-children/`
in `on-shell-counterterms`. This certifies the standalone children; the complete
child under its OS parent, including the consumed measure and frozen kernel,
requires separate qualification on the retained-tensor source.

The pre-finalizer sectors retain full-dimensional numerator algebra and the
nonnegative-epsilon analytic child coefficients. Frozen localization kernels
have unit normalization and must not be integrated again as dimensional loops.
The final four-dimensional emitted integrand is a different boundary. Already
at one loop, setting only the spin-numerator coefficients to four dimensions,
with the same formal-$d$ denominators, loop measure and angular reduction, gives
$f_4=3/epsilon+5+epsilon(9+pi^2/4)+O(epsilon^2)$, whereas the full-dimensional
$f_D$ has finite coefficient 4. Their exact difference is
$ f_D-f_4=-frac(exp(gamma_("E")epsilon)Gamma(1+epsilon)
  (mu^2/m^2)^epsilon,1-epsilon). $
The finite difference $-1$ is an expected evanescent control, not a defect to
remove. At two loops the complete local subtraction and analytic finite
correction must reproduce the corresponding dimensional renormalized quantity,
including evanescent coefficients multiplied by enclosing poles. The matched emitted mass condition and the stronger off-shell dimensional/R2
comparison require separate checks on the complete result. The bare four-color
success alone establishes neither of them nor the integrated condition of the
production OS forest.

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
checks certify that component's representation boundary. The additional
two-loop CFF checks below test the separate energy-integrated limit on explicit
nonpinching rays.

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

The historical repaired-production run
`837a9f7d-2b4c-4ef2-847a-d80f830f573f` passed these three additional tests.
This record predates the retained-tensor port and does not qualify its newly
ported test bodies. The exact scalar residue check takes
0.136 seconds, the full-forest radial matrix 67.055 seconds and the two-point
three-route comparison 33.233 seconds. All remain within the ordinary
240-second per-test limit. These additions change tests, not the production
operator.

On retained-tensor source `6f558ab5979b6e7d591d343559b9d9cf1604458e`, all three
cases passed in run `1e969f55-be49-41ce-a2f4-d4476513a33a` (199.116 seconds total),
with the same exact residue coefficients, nonzero ray controls and ordinary
240-second per-case bounds. Only authored scalar admission and retained
source/definition access were adapted; the production operators and all
reference mathematics remain unchanged. This current run supplements the
historical timings above without certifying intersecting soft strata or a
full integrated two-loop mass coefficient.

The scalar contour oracle retains a doubled massive propagator and a massless
gluon propagator, with numerator $a_0+a_1(z-q^0)$ times an opaque powered
spectator. It certifies the source numerator, denominator routing, masses,
powers and incidence before CFF generation, then compares the result exactly
with the sum of the independently differentiated fermion double-pole residue
and the gluon residue. With the common $i/(2pi)^3$ factor removed, unit spatial
ray $bold(n)$ and $q^0=sigma E_0$, the leading coefficients are

$ lim_(rho arrow.r 0) rho^3 R[1]
  =-1/(8(E_0+bold(q) dot.op bold(n))^2), quad
  lim_(rho arrow.r 0) rho^2 R[z-q^0]
  =sigma/(8(E_0+bold(q) dot.op bold(n))^2). $

Both energy signs, both boosted velocities and four signed unit rays give
32 exact leading-coefficient checks. Scalar contour algebra does not unfold
the spectator or a graph numerator.

The full rainbow uses the projected route with the actual finite child
coefficient and the parent acting on that completed child. The overall
analytic addback has a regular spatial localization kernel and cannot change
the negative powers under test. Keep the child spatial momentum at $(0,9,0)$
and take $bold(r)=bold(p)+rho bold(n)$ for the two unit rays along the first
and second axes. Five radii, $rho=2^(-3),dots,2^(-7)$, at rest and at both
boosts give 50 open-matrix evaluations across the five point/prescription pairs
in the table above. In that historical run, on the last two halvings the $rho^2$-weighted tensor
norm had ratios 0.9956–1.0049 for matching points and 1.9987–2.0246 for boosted
external points with rest subtraction. Thus the sampled matching limits are
consistent with $rho^(-2)$, while the mismatched controls retain $rho^(-3)$.
The assertions retain nonzero controls and exclude an extra integer pole.

The separate complete-forest comparison uses a fixed nonaxial subtraction
point and both matching and unrelated on-shell external momenta. All three
routes agree component by component, including the finite child correction.
This checks that fixing a component point does not replace the surrounding
graph's live momentum. The radial checks concern the open spinor matrix;
they do not establish an external-spinor-projected identity, general soft
integrability at intersecting strata or the full two-loop integrated mass
coefficient. Logs, source guards and the exact residue derivation are retained
under `target/os-feynkit/physics-validation/boosted-cff/` and
`target/os-feynkit/physics-validation/boosted-rainbow-cff-leading-coefficient.typ`.

The original OS run recorded 7.1, 17.5 and 59.9 seconds for the rest, axial and
non-axial cases in `tests/artifacts/on_shell_boosted_validation.json`. These
are historical timings; the rebased stack requires its own ordinary-lane run.

=== Iterated analytic child near the physical shell

The native J0/J1 source above also fixes the off-shell vector coefficient
after the child has been integrated analytically. Measure momenta in m, set
$mu=m$, and remove the common native amplitude phase and the one-loop factor
$a C_F m$, where $a=g_s^2/(16pi^2)$. For the remaining dimensionless matrix
write $A_R(s)slash(q)+B_R(s)I$, with $q=p/m$ and $s=q^2$.
This convention is calibrated by the native mass reference
$3/epsilon+4+epsilon(8+pi^2/4)$ and field endpoint below; it introduces no
new forest sign. These are analytic reference units, not a substitution for
the surrounding graph's live momentum.

Set $delta=1-s$, $d=4-2epsilon$ and
$ H=exp(gamma_("E")epsilon)Gamma(epsilon), quad
  Delta_delta=x[x+(1-x)delta], $
$ J_0(delta)=integral_0^1 Delta_delta^(-epsilon) dif x, quad
  J_1(delta)=integral_0^1(1-x)Delta_delta^(-epsilon) dif x. $
The actual affine numerator $d M+(2-d)K$ and its shift by $(1-x)q$ give
$ A_R=(2-d)H J_1+1/epsilon, $
$ B_R=H[d J_0-(3-2epsilon)J_0(0)]-1/epsilon. $
Here $J_0(0)=1/(1-2epsilon)$ and
$J_1(0)=1/((1-2epsilon)(2-2epsilon))$. Thus
$A_R(1)=-2-epsilon(4+pi^2/12)+O(epsilon^2)$ and $B_R(1)=-A_R(1)$,
as independently certified by the native prefix owner. At epsilon zero,
$ A_R=1+2integral_0^1(1-x)ln Delta_delta dif x, quad
  B_R=-6-4integral_0^1 ln Delta_delta dif x. $

The coefficient needed between the two neighbouring quark propagators is
$C(s)=(s+1)A_R(s)+2B_R(s)$, rather than $A_R+B_R$. For $X=slash(q)$,
$X^2=s$ gives the exact ordered identity
$ (X+I)(A_R X+B_R I)(X+I)
   =C(s)(X+I)+(s-1)(A_R+B_R)I. $
No unintegrated right-projector zero is assumed in this identity.
With
$ F_j=integral_0^1(1-x)^j ln(1+delta(1-x)/x)dif x, quad j=0,1, $
elementary primitives give
$ F_0=-delta ln delta/(1-delta), quad
  F_1=-delta/(2(1-delta))
      -delta(2-delta)ln delta/(2(1-delta)^2). $
Consequently
$ C(1-delta)=2delta+2(2-delta)F_1-8F_0
  =delta(4-4delta-delta^2)ln delta/(1-delta)^2
    -delta^2/(1-delta)
  =4delta ln delta+O(delta^2(1+abs(ln delta))). $
Meanwhile $A_R+B_R=2delta ln delta-delta+O(delta^2(1+abs(ln delta)))$.
The causal continuation uses $delta=1-s-i 0$ and
$"Log"(-eta-i 0)=ln eta-i pi$; the following radial estimate needs only
the positive-delta side.

The rate extends to each fixed epsilon coefficient. Put
$a_s=(2-d)H J_1(1-s)$ and $b_s=d H J_0(1-s)$. Constant counterterms cancel in
differences, so the exact endpoint condition yields
$ C_epsilon(s)=-delta A_R(1,epsilon)
    +(2-delta)(a_s-a_1)+2(b_s-b_1). $
For $0<delta<=1/2$, split every logarithmic moment difference at $x=delta$.
On the inner interval its absolute integral is bounded by a constant times
$delta(1+abs(ln delta)^k)$, since powers of $ln x$ are integrable.
On the outer interval
$0<=ln(Delta_delta/x^2)<=delta/x$; the difference-of-powers inequality
and integration of $abs(ln x)^(k-1)/x$ give the same bound.
The simple pole in H means epsilon degree n requires only $k<=n+1$.
Thus each fixed coefficient of C, and likewise of $A_R+B_R$, is
$O(delta(1+abs(ln delta)^(n+1)))$. This includes epsilon one for one
enclosing pole and epsilon two when two enclosing poles require it.
It is not a bound uniform in epsilon.

For the specified future outer-gluon residue, let $k^0=rho$ and
$abs(bold(k))=rho$, with matched future timelike $p^2=m^2$. Then
$delta=2p dot.op k/m^2$ is bounded above and below by positive constants
times rho, uniformly in this gluon's direction. The already-integrated
child sandwich divided by the two adjacent quark denominators is at worst
$O((1+abs(ln rho)^(n+1))/rho)$. The gluon residue supplies $1/(2rho)$,
so its three-dimensional radial measure has the integrable majorant
$O(1+abs(ln rho)^(n+1))dif rho$ in this iterated channel.

This does not justify integrating the original child first by Fubini.
It supplies no bound for other residue/cut groupings, no majorant uniform
in the third physical/fixed soft normal or the remaining child variables,
and no absolute simultaneous-CFF convergence or UV/IR limit exchange.
The full unregulated physical two-loop R2 claim remains separate.

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

=== Retained private oracle and isolated native references

At retained-tensor source `137afa91448d3c0479ff327e51facbe23fe0853d`, eight
private forest-oracle controls pass in run `01af0491-9e04-44a2-a0d9-1b56021d3cd6`.
The oracle retains typed RHS carriers and reserves existing dummy names from
prepared roots, RHS logical interfaces and the graph projector before
allocating one shared execution-output frame. Local metric bridges preserve
independent additive scopes. Native selectors open valid fixed spinor
components at their graph-owned ports. A ranked matrix zero retains its
ordered interface, while a closed scalar zero stays scalar. Three existing
Idenso reservation controls and one Spenso reservation control also pass.
These private-oracle rules do not modify production counterterms.

The same source passes the native thirteen-source mass-projection control and
both isolated qg-child controls in run `988f0f4b-c1a9-47de-b93b-a1ace0ab8c54`.
Actual typed diagram and child/cograph carriers retain their complete poles,
source factors, phases and formal-dimensional interfaces. The isolated
abelian and nonabelian child results agree with the independent open-tree
Laurent references through epsilon one. The complete boosted finite-child
three-route consumer passes in run `1aa340cb-2a30-4350-8c59-ed183990c98e`
(39.618 seconds). All sixteen selected controls retain the ordinary
240-second per-case cap. Exact source, interface and assertion evidence is in
`target/os-feynkit/retained-tensor-os-rebase/oracle-prerequisite/` and
`target/os-feynkit/retained-tensor-os-rebase/19-native3/`.

The later vertex-prefix checkpoint `01dd0484` retains each certified native
diagram, graph, child, local and integrated payload and generation setting for
its next proof. The new ordinary density control uses those same four actual
attachments (two in D02 and two in D11), with their exact quotient poles,
source/cograph factors and single remaining outer measure. It closes the
formal-dimensional mass trace before comparison with the independent SU(3)
one-loop density. The one-loop integral through epsilon one then fixes the
finite term after multiplication by the vertex pole. The three source and
isolated-child controls plus this proper-density control pass together in
`baf4c286` (2.101 seconds), with all source and binary guards unchanged. This
is an exact density transfer, not a new two-loop numerical integration or a
completed parent-operator proof. Evidence is retained under
`target/os-feynkit/retained-tensor-os-rebase/23-24-vertex-prefix/`.

The private parent-input checkpoint `83fce52c` carries the same four native
children through ordered OS parent operations at rest and at the fixed
nonaxial velocity $(5/4,1/4,1/2,1/2)$. Two scalar-certificate controls preserve
factored cancellation, ranked zeros and nonzero perturbations (`73ec61c0`).
The rest and nonaxial parent runs pass in `3f3e7a3a` and `135b9a28`
(24.629 and 90.296 seconds), all within the ordinary 240-second cap.

The proof retains its 144-component physical external-block certificate,
formal-dimensional internal source, ordered logical ports, complete signed
parent composition and separate vacuum-remainder input. Actual output sectors
retain the full frozen child LMB, in addition to the input and sector-count
checks. Exact common-chart comparisons group complete denominator polynomials
with their powers and prescriptions; they do not establish a new general
per-frozen-domain or CFF/localization commutation theorem. All source and binary
guards remain unchanged. Evidence is retained under
`target/os-feynkit/retained-tensor-os-rebase/parent-input/`. Parent finite analytic
integration, the full26 condition and emitted four-dimensional/R2 equivalence
are not evaluated by these input tests, and public mixed-prescription rejection
remains an explicit control.

These source and isolated-child certificates neither admit a mixed-prescription
parent nor prove the full integrated two-loop OS condition. The mandatory
nonignored numerical bare reference and its pinned backend/Nix qualification
are recorded separately in #link(<mandatory-bare-mass-numerical-reference>)[Mandatory bare-mass numerical reference]. No
three-loop geometry or extensive soft sweep is included in this qualification.

The later retained-tensor checkpoint `f02380d4` passes the exact
four-dimensional finite-field-vertex/live-bridge pinch certificate before
energy projection (`ada21d70`). Its pointwise test then compares all sixteen
components of the full, local and finite contributions across three routes,
including the independently projected pinched sibling (`6ebdf4ff`). The actual
carriers retain their source factors, full poles and frozen child domains;
only independent finite-component proof copies export their expressions.
The pinched scalar roots use the same final mass and epsilon substitutions
as the production projected builder.

Three existing disjoint finite-prefix/product controls, two shared output-frame
controls and the complete raised-LU derivative regression pass in `5ef9fa21`.
The physical GL1 self-energy master also passes its complete-cut, two-seed soft
control in `2e099092`. All nine cases pass the unchanged ordinary 240-second
cap; the largest case takes 202.012 seconds. Their source and binary guards
remain unchanged. The detailed evidence is retained under
`target/os-feynkit/retained-tensor-os-rebase/21-22-durable/`. These exact
reconstruction and finite pointwise controls remain separate from an inclusive
cross-section calculation or the full production OS renormalization condition.


== Validation graphs and support boundaries

The retained-tensor complete three-loop uniform-soft case is compiled and
listed at source `c0aaace381705b33a21fb3a8507c726080716243`. It requires the
production Wood's eight local forest histories and the six certified physical
and component-owned massless planes. Its declared matrix has twenty affine
flats of ranks one through three, two rays per flat and five radii, plus three
rest-center controls: 230 evaluations and 92 tail bounds with radial measure
$rho^(3r-1)$. The fixed-rest, off-shell test uses external points $(5,0,0,0)$
and $(5,1,2,3)$ at mass $m=12$. It does not include finite integrated prefixes,
hierarchical approaches or an all-angle integrability claim.

This current case has *not been executed*. Its native inventory confirms the
ordinary-lane exclusion and existing 7200-second full-main selection. Execution
is deferred until the final batched/hierarchy source and fresh route-cost
review. Compilation and source-owner checks are recorded under
`target/os-feynkit/retained-tensor-os-rebase/17-three-loop-geometry/`; they do
not transfer a historical uniform-matrix pass to this source.


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

#figure(
  grid(columns: (1fr, 1fr), gutter: 8pt,
    [#show-fixture("../assets/typst/portal-graphs/graphs/os-mass-two-loop-D00.dot")
      #align(center)[D00: ghost bubble]],
    [#show-fixture("../assets/typst/portal-graphs/graphs/os-mass-two-loop-D01.dot")
      #align(center)[D01: gluon bubble]],
    [#show-fixture("../assets/typst/portal-graphs/graphs/os-mass-two-loop-D02.dot")
      #align(center)[D02: three-gluon vertex]],
    [#show-fixture("../assets/typst/portal-graphs/graphs/os-mass-two-loop-D03.dot")
      #align(center)[D03: four-gluon tadpole]],
    [#show-fixture("../assets/typst/portal-graphs/graphs/os-mass-two-loop-D04.dot")
      #align(center)[D04: light bubble (five flavors)]],
    [#show-fixture("../assets/typst/portal-graphs/graphs/os-mass-two-loop-D09.dot")
      #align(center)[D09: rainbow]],
    [#show-fixture("../assets/typst/portal-graphs/graphs/os-mass-two-loop-D10.dot")
      #align(center)[D10: heavy bubble]],
    [#show-fixture("../assets/typst/portal-graphs/graphs/os-mass-two-loop-D11.dot")
      #align(center)[D11: crossed exchange]],
  ),
  caption: [Actual native two-loop QCD self-energy incidence. D04 represents
  the five separately retained massless flavors D04–D08. The other seven
  diagrams are distinct sources. Edge IDs, endpoints, particle assignments,
  directions and local source slots agree with the original ungrouped exports.
  These small DOT files are presentation views; numerical algebra uses the
  complete native graphs. In particular, D03 is retained despite its scaleless
  bare integral.],
)
The three-loop fixture tests the complete four-family proper two-loop component
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

=== Historical rebase onto the repaired FeynKit review

The earlier OS functionality stack was replayed onto the exact review snapshot
`e3adefb1fde053256f7c99c1c446e5835069333c`. All eleven separately named
shared algebra repairs were included in that base. The replay retained its
explicit CT-history arguments, native tensor reductions before Taylor
operations, and medium/vacuum CFF options. Each implementation prefix compiled.
Its later runtime validation is recorded in the historical subsection below;
the temporal-trace results retain their original revision only.

The historical combined state manifest version 19 joined that repaired
FeynKit layout with generated fixed OS frame metadata; both earlier
version-18 layouts required regeneration. The current Indexed OS port
uses version 22 as documented above, with the same strict regeneration policy.
The earlier rebase used owned copies when changing borrowed CFF generation
options and propagated errors from the now-fallible parameter builder in the
kernel-jet reference test. Its boosted soft-control test retained zero
finite-field and numerical samples and disabled exact fallback when using the
then-current zero-test API, preserving the previous exact quick-check semantics.

Massless fixtures apply the explicit light-quark massless restriction before
graph import and assert an exact zero quark mass. Raw SM defaults contain a
nonzero down-quark mass, so selecting that graph alone does not define the
massless specialization.

The repaired generator combines the two charge-mate self-energy sources into
one master. The density contract above retains both unit grouping ratios,
source signs, the summed factor, physical cuts and component momentum charts.
Its independent self-energy oracle sums Born densities at $-bold(p)$ and
$+bold(p)$; in the rest-current kinematics, $t(p)$ and $t^3 rho_L(t q)$ are even
under $(bold(q), bold(p)) arrow.r (-bold(q), -bold(p))$. Both vertex spectator
charts remain. The eighteen grouped-master comparisons and historical
twenty-seven ungrouped comparisons are distinguished above; the grouped test
does not separately evaluate an actual GL2 master. Declared grouping metadata
and numerical density agreement do not establish a symbolic identity of the
grouped graphs or supply the observed weight as an expected-value multiplier.

=== Historical validation on the repaired review

The following qualification belongs only to the frozen `faad517b` source and
its original shared repairs. It does not certify the later state-21 retained-tensor
port, the state-22 Indexed OS port, or their deferred runtime inventories.

The separately named shared repair `feynkit-nested-vector-admission`
(`e2d55fbfaef7db3d6f699c10a360bb280a4bfb96`) precedes all thirteen replayed
OS changes. Nested-vector Schoonschip such as
`δ(cind(0),Q(metadata,mink(d)))` is a valid contraction. The typed entry point
now uses its existing dot normalizer before power and interface inference.
This narrow pass preserves explicit notation, opaque scalar payloads and
callback boundaries; full validation still follows normalization. Independent
regressions cover symbolic and concrete dimensions, all four metric signs,
powers, open spectator ports and factorization. The reproducing case failed
before the fix; both new tests and all 992 ordinary Idenso tests pass afterward.
The existing 24 exclusions remain unchanged. Strict Idenso Clippy, formatting,
the documentation checker and the architecture-note render pass.

The frozen runtime candidate is `faad517b089022402dbdc1e01c43dfd39fb73274`.
It passes all 80 bounded OS/core checks in
`52cb8e0e-7307-47ff-9009-47a170ca07e4` (241.969 seconds total; 62.241 seconds
for the longest case under the unchanged 240-second per-test cap). The earlier
mixed soft-gluon/OS-quark component-jet admission failure is included and passes.
Five core boundary controls, six Idenso controls and four state controls also
pass with their exact inventories. These tests distinguish exact algebraic
certificates, local limits and integrated scheme conditions.

The primitive three-route CLI test passes in
`71122de6-9aa6-4cfa-8d35-839960e0b89c`; the boosted integrated comparison
and fixed-velocity state reload pass in
`d988082e-9c40-4b14-9153-f31f1fcc3328`. The three small physical checks and
localization-kernel jet check pass. All 18 grouped finite-density comparisons
retain their original tolerance; the largest relative arbitrary-precision
difference is $1.315 times 10^(-298)$. Runtime evidence, private test reports,
and the 2117-file source guard are retained under
`target/os-feynkit/rebased-review-e3adefb1/nested-vector-refresh/`.
All nine large core checks pass in
`f1b971ab-3e3a-4028-ac2d-a23bb8b6b021` (1340.726 seconds total; longest
case 374.878 seconds under the unchanged full-main limit). These include complete
nested radial soft controls, live cographs, generic and correlated three-loop
hard limits, independent component jets and exact finite coefficients.
The four two-loop CLI cases pass in
`e818ea89-13c8-4290-9770-8ca0197822b6`, retaining the three-route comparison
for quark-in-gluon, nested-quark, soft-gluon-in-quark and disjoint components.
The three-loop CLI comparison also passes in
`643c3576-817b-45eb-906a-51746d99661f` (2516.170 seconds), completing all
eight fresh CLI checks. All three routes agree at both maintained nonzero
control points. All guarded source, binary, configuration and metadata
identities remain unchanged.

The full physical acceptance run passes in
`7fb5f4f7-9639-40c6-8ae3-63bb2d592e95` (2144.550 seconds). It retains
nine NLO cases of 320,000 samples, the 10,000-sample LO control, seeds 1337
through 1345, ten workers and the 7200-second cap. The maximum reference,
auxiliary-scale and cross-route deviations are respectively 1.399, 1.272 and
1.453 reported standard errors, using combined errors for the latter two.
All imaginary-part and uncertainty assertions pass. There are 8323 Quad
samples, no arbitrary-precision samples and no unstable samples. All ten
original checkpoints are retained, including 110 files and 22 compiled objects;
their hashes, source guard and binary identity pass final verification.
Concurrent correctness timings do not supply isolated performance measurements.

=== Historical repaired-review route measurements

All measurements and UV verdicts in this subsection belong to the historical
`faad517b089022402dbdc1e01c43dfd39fb73274` source. They do not qualify the current
retained-tensor port or change its explicit runtime deferrals. That stack completed
39 fresh processes: 18 logging-off
baselines, 18 separate construction profiles and three primitive local-only
baselines. Routes run separately and serially within this task on a shared
host. The main primitive includes integrated terms; the other fixtures are
local-only. These are observations of that historical source, not paired speedups over an older
revision. One- and two-loop measurements are configured with sixteen Rayon
workers; the three-loop measurement uses one. The UV profiles below all use one.
All 180 baseline point/repeat comparisons pass the unchanged strict bound,
with the largest discrepancy using 3.48 percent of that bound.

For compactness, slash-separated cells below give *projected / explicit /
orientation* in that order. Generation and evaluator-build times are seconds;
evaluation means are µs per sample. There are five repeats at each of two points,
with twenty batch means per repeat. Repeat spread and per-point values remain
in the artifacts; they are not independent-run uncertainty. Evaluator
compilation time is zero in every row.

#table(
  columns: 5,
  table.header([Fixture], [Generate (s)], [Build (s)], [Kernel (µs)], [Total (µs)]),
  [Primitive], [0.27 / 0.254 / 0.255], [0.0386 / 0.0433 / 0.047], [1.66 / 1.6 / 2.83], [10.8 / 10.4 / 11.7],
  [Quark in gluon], [8.92 / 227 / 226], [2.27 / 223 / 222], [33.5 / 93.9 / 220], [46.8 / 108 / 237],
  [OS in OS], [5 / 12.3 / 13.2], [1.18 / 10.3 / 11.6], [22.6 / 52.1 / 533], [36 / 67 / 551],
  [Soft gluon in OS], [10.8 / 8.33 / 8.68], [2.83 / 4.38 / 4.88], [67.5 / 30.3 / 766], [82.5 / 44.4 / 785],
  [Disjoint OS], [2.42 / 5.41 / 5.14], [0.974 / 4.62 / 4.41], [24.5 / 16.3 / 39.3], [37.8 / 28 / 53.7],
)

The three-loop evaluation columns use *milliseconds*, with the sample standard
deviation of the ten repeat means shown explicitly.

#table(
  columns: 5,
  table.header([Three-loop route], [Generate (s)], [Build (s)], [Kernel (ms)], [Total (ms)]),
  [projected], [159.22], [46.36], [0.473 ± 0.032], [0.494 ± 0.033],
  [explicit], [1240.32], [1116.30], [4.036 ± 0.732], [4.090 ± 0.756],
  [orientation], [1072.24], [975.42], [284.951 ± 38.069], [285.280 ± 38.112],
)

The following profiles are separate from the baselines. CT orchestration is
inclusive. CFF envelopes include the production root and rejected proposals
once; native CFF is a subset. These columns overlap and must not be added.
The retained CT remainder excludes only profiled CFF work inside orchestration,
so it is not pure Taylor time. Route order remains projected / explicit /
orientation.

#table(
  columns: 4,
  table.header([Fixture], [CT inclusive (s)], [CFF envelopes (s)], [Native CFF (ms)]),
  [Primitive], [0.193 / 0.18 / 0.178], [0.018 / 0.000473 / 0.000448], [0.767 / 0.279 / 0.26],
  [Quark in gluon], [6.66 / 4.32 / 3.44], [0.601 / 0.0118 / 0.0119], [39.3 / 1.35 / 1.39],
  [OS in OS], [3.76 / 1.93 / 1.6], [1.01 / 0.0117 / 0.0116], [11.2 / 1.35 / 1.33],
  [Soft gluon in OS], [7.74 / 3.87 / 3.71], [0.885 / 0.00109 / 0.00111], [13.5 / 1 / 1.01],
  [Disjoint OS], [1.44 / 0.749 / 0.685], [0.312 / 0.0891 / 0.0891], [2.26 / 1.79 / 1.96],
  [Three-loop OS], [111 / 124 / 92.9], [28.3 / 0.0699 / 0.0707], [206 / 20.9 / 21],
)

Peak memory below means process-tree RSS sampled every two seconds through
import, generation and evaluation; brief peaks may be missed. Generation RSS
and Linux VmHWM have separate scopes. Stored sizes count factorized roots plus
retained definitions without unfolding; their root/definition decomposition
remains in JSON. Encoded-byte equality is diagnostic, not an algebraic identity.

#table(
  columns: 4,
  table.header([Fixture], [Sampled RSS (MiB)], [Stored size (KiB)], [Definitions]),
  [Primitive], [89.1 / 93.8 / 93.8], [8.1 / 10.3 / 10.4], [7 / 11 / 11],
  [Quark in gluon], [300.9 / 2908.5 / 2934.9], [282.6 / 1073.5 / 1160.3], [28 / 75 / 75],
  [OS in OS], [99.5 / 214.6 / 280.2], [210.2 / 409.4 / 426.6], [24 / 38 / 38],
  [Soft gluon in OS], [131.2 / 133.0 / 139.4], [359.4 / 918.4 / 1313.2], [26 / 46 / 46],
  [Disjoint OS], [92.6 / 166.8 / 167.3], [119.9 / 111.0 / 136.5], [28 / 62 / 62],
  [Three-loop OS], [1702.6 / 16865.1 / 19310.6], [6506.2 / 28057.6 / 31059.5], [88 / 152 / 152],
)

The matched primitive local-only supplement disables integrated terms. Its
twenty direct-versus-projected point/repeat checks pass the existing strict
bound, with maximum relative difference $5.029 times 10^(-14)$. This numerical
check supplements the separate exact and integrated certificates.

#table(
  columns: 4,
  table.header([Local-only route], [Generate / build (s)], [Kernel (µs)], [Total (µs)]),
  [projected], [0.1678 / 0.0344], [1.339], [10.032],
  [explicit], [0.1586 / 0.0409], [1.363], [10.129],
  [orientation], [0.1475 / 0.0435], [2.835], [11.706],
)

Projected 4D has lower observed generation, evaluation and memory costs for
quark-in-gluon, OS-in-OS and the three-loop chain. Explicit direct evaluates
soft-gluon-in-OS faster at both generic points. Disjoint and primitive local-only
cases retain construction/evaluation tradeoffs. Hard UV rays and precision
retries require their own route-cost controls.

The three-loop orientation benchmark uses the Double-only minimal settings.
Its twenty-sample minimum per repeat gives 200 timed calls and 100 warmups:
57.056 and 28.601 seconds respectively, accounting for most of the 88.349
seconds outside generation. The two-second target calibrates workload; it is
not a measured-duration promise. The JSON has no precision histogram; the
Double-only scope follows the frozen benchmark settings.

Full per-point spread, Spenso/Symbolica build stages, CT/CFF accounting,
component/projection counts, actual workload and memory/size records are in
`target/os-feynkit/uv-perf-review-faad517b/performance/summary.json`.
These diagnostics do not replace the existing paired CI performance gate.

=== Historical complete numerical UV matrix

The six maintained tests pass on the same source in run
`6bace624-e180-48fe-9be2-acd91e130b78`, taking 275.278 seconds in total.
Each test constructs the local-only OS forest and a matched unsubtracted
control. All physical loop-momentum bases and all nonempty hard subsets are
included, with 25 scales from $10^8$ to $10^12$ and seed 1337. The native
degree-of-divergence threshold, fit criteria and precision retries are unchanged.
One test and one Rayon worker run at a time, using the existing full-main
7200-second per-test limit; ordinary tests retain their 240-second limit.

#table(
  columns: (1.4fr, 1fr, 0.7fr, 0.9fr, 0.9fr),
  table.header([Fixture], [Route], [LMBs], [Limits per control], [Whole test (s)]),
  [Primitive], [explicit], [2], [2], [0.607],
  [Quark in gluon], [projected], [7], [21], [17.979],
  [OS in OS], [projected], [7], [21], [6.513],
  [Soft gluon in OS], [explicit], [7], [21], [11.159],
  [Disjoint OS], [projected], [4], [12], [4.459],
  [Three-loop OS], [projected], [26], [182], [234.559],
)

The retained reports contain exactly 259 limits per control, with identical
scales, LMBs and subsets between the OS and bare profiles. All OS limits pass.
The bare profiles have 122 resolved degree-of-divergence failures, four
unstable-fit verdicts and 133 passes; each fixture has a resolved divergence as
required. There are no missing fits.
The nominal 12,950 initial profile evaluations are separate from 147 retained
full-ray arbitrary-precision retries: 21 in OS profiles and 126 in bare profiles.
Each retry repeats 25 points. These counts do not count the internal precision
evaluations made by each profile call.

Route choices use the construction and evaluation measurements together.
Six representative one-worker profiles compare projected and explicit routes
for the soft-gluon, disjoint and primitive fixtures; all twelve native limits
pass with matching paired inventories. The soft-gluon construction-plus-profile
costs are 11.302 versus 8.493 seconds, and the disjoint costs are 2.701 versus
5.771 seconds. These representative profiles do not measure complete bare-control
costs. A separate byte-exact copy of the maintained primitive helper measures
both controls and all LMBs, including the orientation checks where supported:
projected, explicit and orientation routes pass in 0.514, 0.522 and 0.686 seconds.
The first two differ by only 8.1 ms in one shared-host observation, so the
existing explicit choice is retained without claiming a reliable ranking.

The complete three-loop projected test, including both controls, takes less
time than either direct route's measured subtracted construction alone
(1072 and 1240 seconds with one worker). This supports retaining projected 4D
for the extensive matrix. These are observed workload costs, not a proof of a
globally fastest route. Small three-route equivalence tests remain in place.
The primitive orientation-local control also passes separately in
`666270b5-dc8a-4d2a-bbb4-301a27d90a0e`, covering two summed limits and four
orientation entries in each control.

All 2117 guarded runtime sources and the native binary, metadata and test
configuration remain unchanged throughout this run. The twelve original
profiles, native verdicts and inventory audit are retained under
`target/os-feynkit/uv-perf-review-faad517b/uv/`. This numerical UV coverage
supplements the separate exact reconstruction and integrated OS certificates;
it does not establish a complete three-loop radial soft matrix.

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
relabelled as temporal-trace or consolidated-rebase results. The completed
historical repaired-review UV matrix and isolated replays are recorded above;
this supersedes only that earlier source's pending measurement record. No
full three-loop radial soft matrix follows from those numerical UV controls
or their two-loop soft checks.

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
the separate `Performance main` workflow measures every main push against the
previous main SHA, independently of those gates. Manual dispatch remains
available. The trigger and workload policy is maintained in
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
deterministic base seeds 1337–1345, relative-spherical coordinates, map-density channel
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

== Mandatory bare-mass numerical reference
<mandatory-bare-mass-numerical-reference>

The native QCD source certificates feed a separate numerical reference through
an already closed finite scalar codec. The nonignored
`native_two_loop_qcd_bare_mass_matches_primary_reference` test is selected only
by `ci_os_mass_reference`, a dedicated x86_64-linux full-main check with the
existing 7200-second limit. Ordinary and general full-suite profiles exclude
that one external-runtime consumer; their existing selections and time limits
otherwise remain unchanged.

The pinned Python 3.13 runtime retains PySecDec 1.6.6, NumPy 2.5.4, SymPy 1.14.0,
mpmath 1.3.0 and nest-asyncio 1.6.0 with wheel hashes. Two generation, FORM,
compilation and QMC workers, seed 1337, the original evaluation bounds and
coefficient tolerance remain fixed. The two one-loop anchors and seven distinct
two-loop kernels consume all thirteen certified native sources, retaining the
scaleless D03 certificate and exact flavor/color reuse. The reconstructed mass
renormalization constant uses independent analytic insertions; this bare-reference
check does not evaluate the production OS forest's integrated condition.

The embedded Python source belongs to the core compile-time test input; the
external runtime is supplied only to the dedicated main consumer. Missing
runtime or an empty test selection is an error. Failed callbacks retain their
private inputs and backend log, while successful callbacks clean generated
workspaces. Launcher identity and all source, policy and output hashes must be
recorded for a fresh numerical qualification on the migrated stack. Adding this
wiring alone makes no new numerical pass claim.

The retained-tensor preparation at `0ffc3845` compiled the native consumer,
verified its exact profile selection, and passed all three native source and
isolated-child controls (`0e6c0aa0`). Canonical CI regeneration and graph checks
preserved the generated configuration. The resolved runtime derivation is
identical to the previously qualified two-anchor runtime, so its seven cached
checks remain evidence about that immutable backend and input. The subsequent
fresh run on `da40a0b4` (`4ac0dfbd`) completed all nine kernels and all 47
coefficient comparisons in 515.412 seconds. Its native scalar input hash is
`e7ce7a27c02d203386e2a69a809fc033b88c384239fa609b5ed010e55ffab5bd`.
The largest component difference was 0.095561 of the unchanged absolute
tolerance. All 3522 tracked source guards and native/runtime artifact guards
remained unchanged. A focused `uv, integrated, summary` logfile retained all
nine backend reports and their input and library hashes before successful
callback cleanup. The production OS forest condition remains a separate,
unexecuted condition on this retained-tensor checkpoint.

== Indexed OS continuation
<indexed-os-continuation>

The 26 retained OS functionality changes were copied individually above the
qualified shared Indexed checkpoint `ccbba18a`. Each Rust prefix compiled before
the next change; documentation-only prefixes retained exact code bridges. The
preserved replay tip is `2524f931`. Native GammaLoop carriers use explicit
Indexed ports, while generic FeynKit diagrams and independent coefficient-row
work retain their existing Partial owners. Owned checked conversions preserve
payloads, logical order and ranked zeros; they do not assign missing indices.
Combined full-state version 22 first enters with the boosted-frame payloads and
rejects both earlier version-21 layouts. The upstream rational-shell fixture,
shared algebra contracts and public prescription guards are unchanged.

The additive proper-self-energy successor `89e6043c` preserves every existing
native-reference statement and adds only private proof consumers. Its fresh
producer passed 135 distinct baseline tests before all four new proper-prefix
tests. This includes ordered one-loop local and analytic references, exact
completed 4D/CFF reconstruction, direct-route and finite-history comparisons,
boosted rainbow controls, native parent inputs at rest and a nonaxial point,
ranked-zero and logical-port checks, both nonempty fixed-frame persistence
paths, and strict state-version rejection. All 3533 tracked source files and
all five producer binaries remained unchanged across the runs. Every selected
case retained the ordinary 240-second cap; the longest observed case was the
finite-disjoint route comparison at 199.330 seconds. These are correctness
measurements on a shared host, not route-performance rankings.

The physical GL1 self-energy master complete-cut soft test retained both seed
profiles and the existing native summary. The reported measure-weighted scalings are
1.999999999457756 and 1.9999999996998061, each with fit coefficient $R^2 = 1$;
the real-cut negative controls remain near zero. The three routes agree at the
two prescribed points. These are the native test's finite-ray fit and assertion
results, not an independent refit or an all-angle soft proof.

The four additional tests certify the longitudinal odd remainder, the D03
local cancellation before its scaleless integral, the D09 mass/field proper
insertion, and the six fermion-bubble field-pole densities. They retain formal
internal dimension, exact source polynomials, consumed-loop ownership and the
required finite-child epsilon depth. Their success does not establish the
pure-gauge/P2 sum, actual integrated overall parent condition, or emitted
four-dimensional/R2 bridge. The nine-kernel numerical result above remains
historical on its stated source: no new scalar-input bridge or backend run is
claimed here. The held main-only geometry and performance inventories also
remain unexecuted on this Indexed checkpoint.

=== Complete proper-prefix sum and first actual overall input

The subsequent private proof source `544f1e0b` passed six affected and new cases
before its actual D09 parent test. The pure-gauge child proof combines the
native ghost, gluon and tadpole sources before requiring a transverse pole;
its inverse crown map retains the complete tensor value, including the ordered
ports of zero. Together with the quark, vertex and six fermion-bubble records,
the complete fourteen-occurrence proper-prefix sum agrees exactly with the
independent mass reference through the finite coefficient. The unchanged Ward
and D03 controls reuse their exact `89e6043c` passes; their bodies and shared
finite-row owner are byte-identical. No new external numerical backend run is
part of this step.

The D09 test then evaluates the actual parent operator and analytic integral on
both its bare source and its completed OS child, using the same retained graph,
child, settings and frozen sectors. Both finite mass projections vanish after
formal-dimensional contraction and Laurent truncation, without substituting
auxiliary scales. The native log records both input-completion certificates;
the test passed in 26.146 seconds under the unchanged 240-second limit. This is
a two-input qualification. The full twenty-six-input mass condition, public
mixed-scheme admission, and the emitted four-dimensional/R2 comparison remain
separate gates.

=== Actual native twenty-six-input mass condition

The private proof source `c6d1d6c5` passes the complete native inventory and
all twenty-six overall inputs at the matched rest point. The inventory binds
thirteen native sources, fourteen proper occurrences and fifty-two classified
nodes. Independently composed source views can expose the same named ports in
a different logical order: the private reference is transported by the native
axis-permutation owner, with exact canonical and logical interfaces, unchanged
expression, inverse-value and ranked-zero checks. The actual runtime carrier
and the final exact value comparison remain unchanged.

The inventory's 23477-byte scalar payload has the same SHA-256
`e7ce7a27c02d203386e2a69a809fc033b88c384239fa609b5ed010e55ffab5bd`
as the qualified nine-kernel numerical input. The embedded consumer and pinned
runtime launcher also match. This certifies reuse of that numerical result for
the identical input; no new backend integration was performed.

The full-main run `81cba6eb` completed all twelve bare and fourteen proper
inputs, retaining both D03 occurrences. All twenty-six actual integrated finite
mass projections vanish without auxiliary-scale substitution. The exact
operator and source identities then transfer the independently qualified bare
reference and the actual proper densities to the independent source-level
mass condition. For SU(3), five massless and one massive flavor, with the
reference scales and physical mass set to one, the result through the finite
coefficient is
$ K_2 = 82/(3 epsilon^2) + 13/(9 epsilon) - 1285/18
  + 8 pi^2/9 + 8 zeta(3)/3 - 16 pi^2 log(2)/9. $
This equals the independently specified $z_2+w_1 z_1$; the expected value is not
defined by subtracting the produced forest input.

The run took 416.070 seconds under the existing 7200-second full-main cap.
The ordinary inventory retains its 240-second cap. Public complete-forest
admission, emitted four-dimensional/R2 equivalence and the two-loop field pole
are separate outstanding claims. This result does not qualify the deferred
three-loop soft or performance inventories.

=== Physical finite-addback and inclusive NLO continuation

The retained `89e6043c` NLO executable was used with an exact source-closure
bridge to `c6d1d6c5`: production code, API, NLO tests, fixtures and locked inputs
are unchanged; the intervening Rust additions are private core-test modules.
No freshly compiled NLO binary is claimed. The finite-addback test passed all
eighteen master-graph, point and auxiliary-scale comparisons, retaining both
charge-conjugate Born contributions and the original model-coupling rounding
control.

The subsequent full-main run `4b7abc1c` passed LO and all nine NLO measurements
at center-of-mass energy 400, pole mass 173 and renormalization scale 173.
Each NLO measurement retains 320000 evaluations and ten workers; LO retains
10000 evaluations. The analytic NLO reference is 0.43283008677698553.
The table records the actual measurements, with $(M,L)$ denoting the UV and
localization scales.

#table(
  columns: (1fr, 1fr, 1fr, 1fr),
  table.header([Route], [$(M,L)=(173,173)$], [$(346,173)$], [$(173,346)$]),
  [Projected 4D], [$0.430876 plus.minus 0.002979$],
    [$0.431886 plus.minus 0.003062$], [$0.433900 plus.minus 0.003027$],
  [Explicit direct], [$0.436993 plus.minus 0.002976$],
    [$0.433092 plus.minus 0.003061$], [$0.431606 plus.minus 0.003015$],
  [Orientation direct], [$0.429829 plus.minus 0.002976$],
    [$0.430986 plus.minus 0.003057$], [$0.428941 plus.minus 0.003005$],
)

The largest discrepancy from the reference is 1.400 standard deviations.
The largest uncertainty is 0.708% of that reference, below the unchanged 5%
bound. Prescribed scale and route comparisons are below 1.272 and 1.453 standard
deviations, respectively. The original seeds, including their documented stream
overlap, and all tolerances are unchanged. All ten result files and the exact
source/binary mapping are retained. The run took 2438.704 seconds under the
7200-second cap on a shared host; it establishes no route-performance ranking
or all-angle symbolic soft proof.

=== Primitive emitted-dimensional bridge qualification

The primitive regression combines the exact local spin-four source certificate
with the same production dimensional scheme addback. Its continuation retains
$d$-dimensional denominators and angular measure while setting only the already
contracted spin coefficients to four. After stripping $a C_F$, the finite
differences obey

$ Delta B=slash(p)-2m, quad Delta O=-m, quad Delta V=slash(p)-m, $
$ Delta B-Delta O-Delta V=0. $

The actual addback is $+"Fin"[V_d]$. A spin-four addback loses the nonzero open
field matrix $slash(p)-m$; adding it again as an independent primitive R2 term
would double count it. The separate mass difference $-m$ is nonzero. The
regression retains rest, axial and nonaxial subtraction and external velocities,
two auxiliary masses, a one-sided wrong-angular-average control, and the
positive-epsilon difference multiplied by a scalar pole. That last control is
mechanism-only, not a nested-graph certificate.

This extension was qualified on historical source `0b4c6a9c`. The Indexed
source `ea715048` now passes the three extended integrated references and four
primitive supports: the compact local operator, independent Dirac CFF residues,
fixed-point denominator contractions, and live-crown/prescription reconstruction.
All seven retain the ordinary 240-second cap, their exact formulas and negative
controls. The same production integral supplies both pole and finite parts;
no extra backend call or counterterm is introduced. Exact numerator and complete
denominator certificates precede CFF aggregation in the existing primitive
projection consumer. Scalar dispatch roots and retained tensor bodies keep their separate
normalization owners. The actual D09 finite-child dimensional contribution has
its own certificate; neither result establishes the complete two-loop emitted
condition, general off-shell/R2 matching, or public mixed-prescription support.

=== Owning Indexed FeynKit rebase qualification

The replay descends from FeynKit `8959f063dd3a`, preserving the original OS
input `35d011195068` and its earlier qualified checkpoint `170773df7e1c`.
The two named bounded-coordinate proof changes precede the thirty-eight OS
functionality changes. All twenty-seven Rust-changing OS prefixes compile;
the eleven remaining prefixes retain the preceding checked Rust sources.
Transient Indexed products move their operands, while retained or shared
operands are cloned explicitly. The generic PartialStructure reference API
keeps its borrowed operation. The final API-state formatting uses that crate's
declared Rust edition 2021. An independent source review preserves the physics
assertions, ordering, factorization and ownership rules.

Fresh producer source `8fb7c81b45ef` passes 235 distinct ordinary cases and
two full-main cases. These cover the local and analytic one-loop references,
primitive dimensional bridge, completed CFF reconstruction, fixed boosted
points and ordered projectors, complete rainbow soft controls, nested and
disjoint route comparisons, native proper-prefix sums, physical complete-cut
soft behavior, and state-format 22 persistence. Exact source, binary, selector,
JUnit and profile bindings are retained. No case is skipped or retried. The
ordinary cap remains 240 seconds, with a longest observed case of 205.896
seconds; both larger cases retain the 7200-second full-main cap.

The standalone denominator-source regression also passes its previously
compiled fixtures. Five vacuum occurrences agree after their explicit
auxiliary-mass identification, with the differing literal comparisons recorded
separately. Two retained physical-crown occurrences agree literally as well.
Raw-prescription drift and six ownership or metadata corruptions remain
rejected. This fixture does not certify the separate opt-in Wood-hook consumer.

The full nested CLI case passes in 427.341 seconds. Its two points, three
routes and five repeats yield thirty emitted values; repeats are identical
within each route. An independent review verifies every direct/projected pair
against the stricter binary64-epsilon floor of the original accuracy bound.
This numerical check supplements the separate exact source certificates.

The twenty-six-input mass case passes in 408.745 seconds. Independent review
checks all twelve bare and fourteen proper inputs in order, all 104 literal
projector attachments, and all twenty-six exact zero finite mass projections
without auxiliary-scale substitution. The final scalar coefficients agree
exactly with the independently specified $z_2+w_1 z_1$ reference above.
The matched-rest, formal-dimensional SU(3) reference scope remains unchanged.

The fresh 23477-byte scalar payload, numerical consumer and pinned runtime
launcher match the qualified nine-kernel input exactly. All forty-seven
retained numerical reference checks are independently rechecked; the backend
integrals are reused without a new numerical integration. The fresh physical
finite-addback case passes all eighteen density cells, with an independent
review of their emitted distances and coupling-rounding control. Individual
Born and separate on/off inputs are not emitted by this artifact, so that
review does not independently reconstruct the Born-density oracle.

Strict Clippy, formatting, documentation checks and canonical CI regeneration
pass. Generated CI configuration is unchanged and implementation CI remains
disabled. Timings and sampled memory are shared-host diagnostics, not controlled
route rankings. The inclusive Monte Carlo result remains historical; its
separate full-main case is retained but was not rerun for this rebase. Public
complete-forest admission, complete two-loop emitted/R2 matching, the two-loop
field pole, and the larger scaling and performance matrices remain outstanding.



=== Physics continuation on the owning stack

The native D11 and D02 Wood regressions retain six histories each: the empty
history, either overlapping child, the parent, and either child followed by its
parent. They exclude the incompatible simultaneous children. The four local
and completed-finite-prefix cases pass in 58.985 seconds. Their 108 actual
source preparations and 336 denominator occurrences receive exact certificates;
the completed children retain nonzero pole, finite and positive-epsilon parts.
Parent operations consume each complete child in the same fresh context, with
its original frozen frame and coefficient owner. Public mixed-prescription
admission remains unchanged. These exact construction tests do not replace UV,
soft, emitted or integrated renormalization checks.

The D09 rational-finalizer proof passes for the bare source and completed OS
child. It matches full active/frozen domains before comparing the source and
signed local result, normalizes each term separately, and verifies the actual
integrated finite carrier through recursion before the mass projection.
Independent unit-mass and nonzero kernel controls reject wrong local and finite
mass insertions. Complete native representation slots own cooked index
admission. The ordinary thirteen-source regression independently checks both
native open and closing interfaces, logical ordering, inverse transport and
ranked zeros, rejecting contextless index conversion through the same guard.
It takes 0.367 seconds without a backend call.

Fresh source `442c9b49` also passes the actual D09 parent-CFF construction
prefix in 44.325 seconds. Bare and completed-child projection use fresh contexts
with four and twenty-eight full denominator certificates, respectively. Each
count equals the actual source-preparation count. The existing canonical
numerator and immutable selected-assignment certificates run before CFF
aggregation. The grouped full domains retain their exact frozen kernels.
The actual finite localizer's ordered call arguments reconstruct its typed
carrier against that branch's authoritative energy map before finalization.
The final builder consumes the actual local and integrated operands while
retaining numerator bodies separately from dispatch roots. These are exact
construction and transport claims. Emitted source/local mass cancellation and
complete dimensional/R2 matching are still separate gates.

The inclusive physical acceptance was rerun after the rebase. Its guarded
retained `8fb7c81b` executable has exactly the current production and runtime
inputs; only private cfg(test) code/declarations and documentation changed.
No new NLO executable is claimed.
The new run passes LO and all nine route/auxiliary-scale measurements in
2458.497 seconds under the unchanged 7200-second full-main cap. An independent
420-digit review binds all ten raw results and recomputes the analytic
reference and all twenty-one comparisons. The largest reference, scale and
route differences are 1.400, 1.272 and 1.453 standard deviations. The largest
reported uncertainty is 0.708% of the NLO reference, below the unchanged 5%
bound. Seeds, sample counts and tolerances are unchanged. This is integrated
physical acceptance with its prescribed statistical comparison; it establishes
neither an exact algebraic identity nor an isolated performance ranking.

=== Generic-momentum two-loop field reference

The independent MS field convention uses $a=alpha_s/(4 pi)$,
$psi_B=sqrt(Z_psi)psi$, and the renormalized Feynman gauge parameter
$xi_R=1$, with $xi_B=Z_3 xi_R$. Appendix A.18 of
#link("https://arxiv.org/pdf/1909.00697")[Falcioni--Gardi--Milloy] gives
$w_2=C_F((C_A+C_F/2)/epsilon^2
+(-17 C_A/4+3 C_F/4+T_F n_f)/epsilon)$.
For SU(3) and six active flavors this is
$w_2=44/(9 epsilon^2)-35/(3 epsilon)$.
The gauge convention is independently fixed by
#link("https://arxiv.org/pdf/1709.08541")[Chetyrkin et al., section 2].
Holding the bare gauge at one instead would shift the double pole by
$4/(3 epsilon^2)$; the reference does not admit that alternative.

The inverse-propagator Lagrangian gives
$i a^2(w_2 slash(p)-(w_2+K_2)m)$ with $K_2=z_2+w_1 z_1$.
The separately local OS mass branch owns $-i a^2 K_2 m$.
Here the signed vacuum term is $-(U_1-U_1O)$. Its integrated pole must be
$i a^2 w_2(slash(p)-m)$, independently of the produced mass projection. The finite addback remains
$+"Fin" integral (U_1-U_1O)$ in the convention defined above; its sign is not
the sign of this pole term.
The native one-loop source fixes the phase and loop measure with
$w_1=-C_F/epsilon$. Equation 5.11 of
#link("https://arxiv.org/pdf/2203.11038v1")[Capatti--Hirschi--Ruijl]
provides the same local-mass/vacuum-field separation.

The actual-native check reuses each of the twenty-six existing
integrations, preserves the full generic external momentum, and proves an
exact affine decomposition before any on-shell specialization. Complete native
spin/color slots, inverse interface maps and ranked zeros align the actual
open matrices. Scalar color invariants are evaluated by the existing Idenso
owner in the declared SU(3) convention. The captured pre-normalization scalar
regression proves its cancellation for an arbitrary auxiliary logarithm.
The original expected coefficients and all sensitivity controls are unchanged.
The full-main test at source `6a895deb` passes in 430.658 seconds. It checks
the complete generic-p kinetic matrix first, then the independent constant
$-i a^2 w_2 m I$ and their full sum. Physical mass remains symbolic and no
auxiliary-scale substitutions enter these comparisons. It retains the
independent OS mass reference through its finite coefficient. The wrong bare
gauge convention, simple-pole shift, constant mass shift, matched-point-only
matrix and quadratic dependence remain rejected; the first four produce actual
nonzero finite-scalar witnesses rather than scope or budget refusals.
An independent reader verifies all twenty-six input extractions, their native
interfaces, 497 full polynomial denominator classes and 104 literal projector
attachments, then checks the final field and mass references and control
records. This qualifies the formal-dimensional integrated pole matrix in the
specified SU(3), six-flavor, fixed-rest-subtraction convention. Finite off-shell,
emitted dimensional/R2 and public mixed-prescription claims remain separate.

=== Actual D09 source branches and complete finite products

The D09 source-construction extension at `a8ec1ab7` passes in a 46.296-second
ordinary run. It builds the bare source through the actual direct root and
final builder. The completed-child source retains its active local sectors
and actual finite completion separately, with exact recomposition, consumed
one-loop measure and frozen-frame checks. The finite completion passes through
the real localizer and the same-key outside-numerator multiplication.
No analytic backend call is repeated.

Only a cloned projection graph has its external projector replaced by one,
leaving the original graph and actual operator results unchanged. The parent
is built in this same open external-matrix view from the already certified
projected coefficients and actual Integrated object. Before materialization,
every nonzero-coefficient finite-source monomial must contain its actual
retained child family and cograph family, with no family hidden in a scalar
coefficient. Exact residue keys and complete ordered arguments bind both typed
bodies. Their native tensor product must reconstruct the complete open
external matrix of that product, with its original spin/color interface,
under that exact source map.
The run certifies eight complete products under eight distinct full residue
keys, although only six literal call-pair monomials are distinct. Identical
calls therefore do not determine their complete source map. The earlier
first-monomial prefix remains a separate historical qualification.

At this construction checkpoint these are transport certificates, not an
emitted mass condition. A subsequent mass projector must act once on each complete product;
projecting its child and cograph independently would change the operation.
Full cancellation still requires parent key/domain association and exact
component-aware transport of compact energies, including the implicit momentum
dependence of physical on-shell-energy functions. No map is inferred from a
low-rank numerator image. The source construction retains compact energy and
surface functions and separate numerator stores. It does not clear denominators
across a complete sum or resolve the graph numerator catalog.

This construction checkpoint retains its stated source and scope. Later
qualifications below separately record the formal-dimensional full26 recheck,
open-prefix/angular mechanisms and the actual emitted D09 mass condition.

=== Qualified open-prefix and angular mechanisms

Source `1c3d8c31` freshly passes five ordinary support cases in 1.366 seconds
of JUnit suite time: the open D09 proper prefix, captured finite-Laurent row,
fixed-component Wick projection, joint physical-four angular moments and exact
Laurent-grid admission. They use the same retained executable and 2700-file
source closure as the emitted-mass case below, with the unchanged 240-second
failure bound.

The actual dimensional child is attached to its doubled-quark cograph before
projection. At the unit mass and scale point, independent F/H moments fix the
signed proper insertion. Its two open coefficients obey
$ P_d-P_4=C_F^2[(-5/epsilon-17)slash(p)+(10/epsilon+24)m I]. $
Their mass trace is $C_F^2(5/epsilon+7)$ and does not define those coefficients.
The mass-invisible $slash(p)-I$ control is independently nonzero: closing with
$I-slash(p)$ gives $-24$, or $-2$ after division by twelve. Dropping the actual
finite child's epsilon-one term loses $-C_F^2 pi^2/4=-4pi^2/9$ through its
cograph pole. Both spin companions share that loss.

Here $P_4$ changes the spin coefficient within the same formal kk/pk invariant
family. It does not lift physical-four loop components into the dimensional
angular integral. Complete native words and representation slots remain
outside scalar Laurent truncation; the captured finite scalar coefficients
retain their factorization and exact nonzero controls.

The separate angular cases use the existing dimensional tensor projector to
check fixed-component Wick signs, the $4/d$ moment and the joint rank-four
moments of $k_4^2$, $l_4^2$ and $k_4 dot l_4$. Their early-four controls retain
both epsilon-one/simple-pole and epsilon-two/double-pole finite losses, with a
factorized spectator. These qualify the angular mechanism; they do not bind
that lift to the actual complete D09 vacuum inputs.

=== Fresh full26 field and mass closure

After the finite-Laurent proof-owner extension, source `65c8048b` reruns the
complete twenty-six-input field and mass case in 411.467 seconds of case time,
under the unchanged 7200-second full-main bound. All twelve bare and fourteen
proper inputs pass, with 497 full denominator/frame certificates and 104
literal projector attachments. Every actual integrated finite parent mass
projection is exactly zero. The generic-p kinetic and constant coefficients
still give the independently specified $i a^2 w_2(slash(p)-m I)$ pole above,
with physical mass symbolic and no auxiliary-scale substitution. The wrong
gauge, simple-pole, constant-mass and matched-point controls retain actual
nonzero witnesses; the quadratic control reaches the non-affine guard.
Independent exact-rational review also reproduces the proper sum and $K_2$
through its finite term. The 23477-byte bare reference input remains identical;
no numerical backend integration is repeated.

The subsequent changes through `1c3d8c31` affect only the private D09 test
consumer and its grid regression. All other tracked files, the full26 driver
and its shared parent/frame methods retain the qualified source closure.
The full26 result belongs to its retained source-65 executable; it is not a
fresh full26 run on the later child.

=== Qualified D09 emitted matched mass condition

Source `1c3d8c31` passes the actual bare and completed-child emitted D09 mass
case in 118.283 seconds of JUnit suite time, run
`70d1d7d8-f6b8-4eb2-a2c4-ec23c99d40dd`. Independent review binds its exact
selector, retained executable, all 2700 source records and runtime artifacts.
The qualified point has unit physical mass and rest velocity, with
$p=m v$. This is a private D09 condition, not public mixed-prescription or
complete emitted all-source admission.

Before Phi, the actual raw frozen child grown by its parent cograph must equal
the constructed physical-four source: one global numerator, identical full
prescribed denominators, their loop/tree partition and the original frozen
frame. The F/G call witness is a separate transport check. The actual final
builder retains complete key/domain inventories and ordered arguments before
cancellation; its captured chunks recompose the factorized final payload.
The mass projector acts once on each complete child/cograph product.

Every actual row is covered once. Bare coverage is 49 rows: 33 exact zeros and
16 members in five signed coefficient classes. Completed coverage is 287 rows:
195 exact zeros, 52 rows in eleven direct and one reversible-invariant group,
and 40 rows in two exact linear triples. The five bare and seventeen completed
parent chunks retain their full key incidences and ordered row ranges. Actual
scalar-coefficient cross-relations select one disjoint two-triple cover; no
key or body independence is inferred from a numerator image.

Reversible whole-term views retain the complete invariants and products.
The physical energy pair satisfies $C^2-B^2=1$, with $C>=1$, so the known unit
rescaling and $t=B+C>0$, $B=(t-t^(-1))/2$, $C=(t+t^(-1))/2$ define an exact
rational chart. Recursive syntax bounds give degrees at most two in each of
the first two variables, four in the coupling, and t exponents from minus
four to four. Conceptually multiplying by $t^4$ gives separate degrees
$(2,2,4,8)$. Each full $3 times 3 times 5 times 9=405$ distinct-node exact QQ(i)
grid is unisolvent. The two original factorized inputs are therefore exactly
zero, with 810 zero checks and 810 separate original-input-plus-one equals-one
checks. These controls concern the unit-rescaled chart; the independent native
mass insertion supplies the physical nonvacuity check.

Original group-sum-plus-one controls are also retained before compaction.
Where whole-expression collection obscures the unit, the existing exact
one-coupling owner collects only the actual registered coupling G. Its constant
row is literally one and its complete G-free coefficient of $G^4$ is exactly
zero; all other graph factors remain opaque. Exact row recomposition binds
this proof to the original input. No graph numerator or Phi polynomial is
materialized, and no common CFF denominator is cleared.

The separate observation child `28bdf6aa` passes the same mass case in
121.351 seconds. Before any physical diagnostic normalization it retains
exactly three original signed analytic inputs: the bare two-loop sector, the
proper-local two-loop sector, and a one-loop proper-finite sector with the
child frame frozen. Their complete active/frozen frames and provenance remain
beside their factorized numerators. All three retain formal $d$; only the finite
sector contains explicit epsilon, the renormalization scale and the
consumed-loop marker. All three analytic
payloads differ from their physical local expressions. This records the actual
integration policy without introducing another integration path.

Full R2 still requires the physical-four angular lift on the actual bare,
proper-local and frozen-finite vacuum inputs, including their dimensional
angular projectors, physical kernel measures and required child epsilon depth.
A private D09 W observation can retain those inputs and their provenance; it
is not that identity. An independent IR-safe rearrangement or uniform UV/IR
exchange condition is also required. Neither the emitted mass zero nor the
formal spin companion establishes general off-shell HV/CDR equivalence.
Public complete-forest admission and the larger scaling/performance matrices
remain separate. Shared-host case times are correctness records, not rankings.

=== Actual vacuum ranks and sixth-order angular checks

The retained source-28 inputs also pass a factorization-preserving rank
inventory through the original rational-shell parser and their complete frame
maps. It retains eight numerator terms and 42 denominator occurrences. The
largest simultaneous active-degree bounds are four for the bare input, six
for the proper-local input, and four for the proper-finite input. The one
degree-six bound has separate loop bounds two and four; a conservative bound
does not prove a nonzero leading coefficient. The finite child's frozen
coordinate remains outside the active numerator rank.

The bare denominators include the coupled channel $(l-k)^2-m_"UV"^2$. The
proper-local denominators instead depend separately on k and l, with effective
powers $(3,3)$, $(3,4)$, $(2,3)$ and $(2,4)$. The finite input has one active
channel with powers three and four. An independent proper-local reference can
therefore use sequential Wick averages and one-loop tadpoles; the raw total
degree bound alone does not require an irreducible joint sixth moment.
For the dimensional comparison, physical-four projection must act on the
separated numerator. The negative
denominator shells retain their original dimensional momenta and prescriptions;
a physical-four numerator norm cannot be pinched as a full-dimensional norm.

Source `3f0be4d6` additionally passes the extended angular mechanism and the
fixed-component Wick case in 0.895 seconds of JUnit suite time, with strict
Clippy passing. The same native coupled-vacuum tensor projector now checks
$K L^2$ and $M^2 L$, where $K=k_4^2$, $L=l_4^2$ and $M=k_4 dot l_4$.
For a rank-r projector, the independent moments
$E_2=r(r+2)/(d(d+2))$ and $E_3=r(r+2)(r+4)/(d(d+2)(d+4))$ give
$B=(r E_2-E_3)/(d-1)$ and $C=(E_2-E_3)/(d-1)$ through its trace and
idempotence. Their r=4 contractions agree exactly with the native reduction.
The single-pole and double-pole controls retain the independent epsilon-one
and epsilon-two coefficients, with the original factorized spectator intact.
These are additional angular checks. They establish neither the complete
actual D09 angular lift nor its open R2 residual or UV/IR limit exchange.

=== Actual proper-local D09 vacuum scalar reference

Source `d52a04b9` passes the original proper-local D09 vacuum comparison in
48.380 seconds of JUnit suite time, run
`d0d1f020-b783-4fb5-a2f1-eca80f03c756`. Six supporting cases pass in 1.851
seconds: the homogeneous external frame, native child mass/field insertion,
factored Laurent grid, captured finite coefficient and both angular cases.
The frame test admits its retained formal-d or physical-four representation
and rejects mixed or unsupported dimensions. The emitted bare/completed-child mass proof also passes freshly on this source
in 112.074 seconds, run `a29efe15-6d66-4af5-a714-93460f12fdbf`.
The earlier full26 result remains
bound to its source-65 executable; it is not a fresh run of this extended frame
owner.

The input is the ORIGINAL signed analytic W with two active components and no
frozen child. Exact native comparisons establish parent-U absorption of the
child O-minus-UO terms, UO=OU, and direct equality of W with the absorbed-child
construction. These are separate certificates on complete numerator words and
full denominator polynomials. The four source-edge power inventories remain
$(1,2,1,1,1)$, $(1,2,1,1,2)$, $(1,1,1,1,1)$ and $(1,1,1,1,2)$ on edges two
through six. Independent zero-momentum source algebra checks
$i G^4 C_F^2 d^2$ before the new scalar integrals. Neither the graph phase nor
its source weight is fitted to their outputs.

The source-derived identity and independently authored gamma-zero dual share
the actual named spin/color ports. Their complete two-by-two Gram matrix is
exactly the identity, including both zero cross traces. One shared endpoint
owner preserves the inverse map, original interface and ranked zeros in both
this test and the earlier generic-p child test. The existing finite component
oracle closes the finite Gram and the physical-four numerator; negative
propagator shells never enter it. No scalar d or epsilon coefficient is frozen.

Write $T_2=exp(gamma_E epsilon)op("Gamma")(epsilon)$ at unit scales.
Independent radial recurrences give
$I[k^2/(k^2-1)^3]=d T_2/4$ and
$I[k^4/(k^2-1)^4]=d(d+2)T_2/24$.
Keeping the child sign and ordered Clifford product gives
$c_d=(d-2)^2(d+2)/24$ for the proper-local W. For the physical-four numerator,
the rank-four Wick factor $6/(d(d+2))$ instead gives $c_4=1$; the angular
integration dimension remains d. This is not obtained by setting d to four in
an already integrated scalar coefficient.

With $a=G^2/(16 pi^2)$ and the independently retained source weight $omega$,
put $Omega=-i a^2 C_F^2 omega$. At the unit rest point, the independent reference matrices are
$Omega c_d T_2^2(gamma^0-I)$ and $Omega T_2^2(gamma^0-I)$. This gate compares
only their TWO scalar duals: the identity dual gives minus each scalar factor
and the temporal dual gives plus it.
The absolute references, including finite terms, are
$c_d T_2^2=epsilon^(-2)-7/(3epsilon)+5/3+pi^2/6+O(epsilon)$ and
$T_2^2=epsilon^(-2)+pi^2/6+O(epsilon)$.
All four exact Laurent differences vanish, and four additional checks retain
the independently nonzero double-pole coefficients. Native SU(3) invariants
are reduced by their existing owner before comparison to numerical
$C_F=4/3$. The finite scalar certificate has no color or Clifford fallback.

Only separated numerator representations are projected. All negative
propagator polynomials, masses, powers, prescriptions and active frames stay
dimensional. The physical mass and its linked subtraction point are fixed
only after the actual parent operators. The existing Integrated owner receives
signed W once; its pole minus finite-counterterm payload reconstructs I(W).
Original local and analytic inputs remain unchanged. This establishes the
proper-local scalar channels, not completeness of an open two-word basis.

The bare connected vacuum input and the frozen finite-child input still need
their own native angular/reference admissions. In particular, the frozen
child's epsilon-one term can multiply the outer pole, and its consumed-loop
marker and vacuum mass must undergo the parent's hard scaling before any unit
specialization. A complete open reconstruction and the three shifted soft
surfaces remain separate gates. No full emitted R2 or general HV/CDR claim
follows from these scalar comparisons or from the matched mass identity.

=== Independent references and qualified OS mechanism

The final two changes separate complete finite-matrix/pole certificates
(`e7036dfc`) from the source-bound analytic child reference (`3b04cc25`).
They sit above OS base `3dc81ff7`, shared dual-port fix `0a24ecb9` and
FeynKit `d39d9bf2`. That shared repair recognizes explicit and compact
represented dual ports in the existing direct-argument tensor filter;
foreign function scopes remain opaque. Original named references are retained.
State format 22, all 74 protected scopes and the existing CI/performance
policy pass the pinned audit. The complete final tree equals `1abbe1ba`.

The first change passes its actual core-only library-test build in 389.255
seconds and strict Clippy in 19.755 seconds. Four ordinary regressions pass
in 0.349 seconds, run `ce33faac-c49e-46b3-bf02-03e468ac4430`: both new
certificate cases, the dimensional-radial finite basis and captured Laurent
recomposition. The final source passes Clippy in 22.497 seconds and its full
producer in 48.171 seconds. Original pole source, momentum, mass,
prescription and power stay literal. Matrix weights and consumed-child scale
powers have independent wrong-input controls; graph factors remain intact.

The source-bound one-loop reference uses the actual quark/gluon poles,
$Delta=x m^2+(1-x)lambda^2-x(1-x)p^2$, and the affine numerator
$d M+(2-d)K(k)$, with $F=K(p)-M$ on its original ports and phase.
At $lambda=0$, $p^2=m^2$, and unit scales,
$J_0=1/(1-2epsilon)$, $J_1=1/((1-2epsilon)(2-2epsilon))$ and
$G_e=e^(gamma_E epsilon) op("Gamma")(epsilon)$ independently give
$G_e(3-2epsilon)J_0=3/epsilon+4+epsilon(8+pi^2/4)+O(epsilon^2)$ and
$G_e(2-d)J_1+1/epsilon=-2-epsilon(4+pi^2/12)+O(epsilon^2)$.
The second value is the MS-subtracted endpoint coefficient, not an OS field
renormalization condition. These anchors retain formal-dimensional numerator
algebra, one Minkowski measure and one forest sign. Depth one serves one
enclosing simple pole; the separate depth-two gate covers two poles.
Positive lambda certifies a different source denominator, not these
lambda-zero Beta values.

The child gate on `1abbe1ba` passes in 2.963 seconds, run
`1a303bee-759a-4b70-8826-0d4e8641efc9`. Eight ordered products cover both
signs, rest and nonaxial boost, and both fermion endpoints. Only after the
complete ordered product has no unconstrained Lorentz direction does the
existing finite matrix owner admit components zero through three. Scalar D
is retained. Adding a native mass carrier must fail on the same endpoint.
The scale/pole certificate passes in 0.109 seconds, run
`d7411cb5-eeac-4b76-b2fd-8b0a47305202`, including changed source, momentum,
mass, prescription, power and consumed-child scale controls. Fixed on-shell
points remain local to their component.

At unit mass, unit scale and rest, divide the independent original-W rows
below by $Omega Y$, where $Omega=-i a^2 C_F^2 omega$ retains the actual source
weight and $Y=gamma^0-I$. Every $W=U O-U$ already owns its forest minus.
Both integrations remain D dimensional; subscript four changes the active
physical numerator, not the denominators or angular measure.
#table(columns: (2.4fr, 2fr),
  table.header([Vacuum-policy difference], [Coefficient through finite order]),
  [Bare $I_D(W_4-W_D)$], [$-7/(6epsilon)+23/36$],
  [Proper-local $I_D(W_4-W_D)$], [$7/(3epsilon)-5/3$],
  [Frozen-finite active $I_D(W_4[C_e]-W_D[C_e])$], [$4/3$],
  [Frozen child $I_D(W_4[C_0]-W_4[C_e])$], [$pi^2/12$],
  [Complete UV-vacuum sum], [$7/(6epsilon)+11/36+pi^2/12$],
)
Here $C_e=1-pi^2 epsilon/12+O(epsilon^2)$ and $C_0=1$. The unknown absolute
bare finite master cancels in its difference; no absolute bare finite
reference is inferred. Omitting the bare angular epsilon-two term shifts
its finite row by $5/36$; omitting the child's epsilon-one term loses the
distinct $pi^2/12$.

Actual bare, proper-local and proper-finite W comparisons on `1abbe1ba`
pass in runs `4983332c-e5a3-453e-b648-a3989bbab6ec` (17.095 seconds),
`21e386b3-1451-4a10-a9e8-66b5b7b1b0b7` (135.625 seconds) and
`8710cfe0-5860-4679-aedf-11770b73e443` (50.704 seconds).
Source, pre-master open-span and fourteen scalar-reference checks bind the
independent coefficients to all three actual inputs and both calibrated
duals. The aggregate is checked from their actual rows; the connected finite
master cancels before summing sources. This is a UV-vacuum policy difference,
not the physical bare two-loop R2 or a general-scale independence claim.
With $r=11/36+pi^2/12$ and $A=-"Fin" I(W)$, the finite addback difference is
$A_4-A_D=-Omega r Y$. The existing analytic addback already uses D, so no
additional counterterm is inserted. Its minus is distinct from the forest sign.

The earlier 39-case ledger is complete on physics source `3b04cc25`, producer
41. Five native mechanism cases retain their original `1abbe1ba` source,
UUIDs and times: all 2702 guarded inputs and the actual core executable
are byte-identical. The public case retains source `3dc81ff7` through its
identical UV executable and unchanged guarded runtime inputs. The remaining
33 cases execute freshly on `3b04cc25`; exact case identities count once.
The closed fresh batches and their suite times are:
#table(columns: (1.7fr, 2.7fr, 0.7fr),
  table.header([Batch and case count], [Run UUID], [Seconds]),
  [Ordinary references, 17], [`9659a679-be58-4c21-865d-90b20d15bba8`], [35.620],
  [Full26 field and mass, 1], [`69a43e37-759c-4c04-bb71-56ce104ef371`], [443.606],
  [Emitted mass and grid, 2], [`fee073ec-8fb8-4072-a061-b8f9e599aca6`], [132.205],
  [Soft hierarchy, 1], [`01a01b4b-f8fc-4765-a773-1240d9e4646c`], [124.517],
  [Complete UV matrices, 6], [`de88b39d-bdb1-41e4-97bf-c7360ae75f0c`], [528.037],
  [Graph admission, 6], [`c19bab24-75f2-4d5a-93da-9ff8a793259e`], [0.088],
)
The ordinary batch includes the eight rest/boost references, exact soft
coordinate geometry, finite matrix/Laurent regressions and two regulator/frame
checks. These are actual batch times, not inferred subgroup or performance
measurements. Ordinary cases retain 240 seconds; existing full-main cases
retain 7200 seconds. The independent consumers and exact case ledger preserve
every source, selection, runtime artifact and retained executable identity.

The fresh full26 case contains twelve bare and fourteen proper inputs,
497 denominator/frame certificates and 104 literal projector attachments.
All 26 actual finite matrices are nonzero while their mass projections vanish,
without auxiliary-scale specialization. The generic MS field polynomial is
$w_2=44/(9epsilon^2)-35/(3epsilon)$, with its independent wrong-reference
controls. The dimensional source/loop normalization is checked independently;
this gate does not repeat the historical bare numerical integration.

The complete emitted matched-rest certificate proves $Phi[R]=0$ before
integration, with original factorized numerators and full prescribed poles.
All 336 rows are covered once. Bare coverage is 49 rows: 33 native-certified
zeros and 16 signed-class members; only 17 are raw literal zeros. Completed
coverage is 287 rows: 195 native-certified zeros, 52 direct/reversible members
and 40 members in two exact linear triples; only 63 are raw literal zeros.
Two 405-node QQ(i) grids certify both original triples with matching
original-plus-one controls. Source/parent incidences remain 8+41 and 54+233,
in five and seventeen parent chunks respectively.

A well-defined linear integration of this complete zero remains zero; the
matched mass condition needs neither removal of lambda nor exchange of
D-to-four with integration. It assigns no values to divergent isolated terms.
The independent dimensional full26/mass references fix the physical QCD
normalization. The complementary field direction $Y=gamma^0-I$ is invisible to Phi.
The calibrated native diagnostic $Psi=(gamma^0-I)/(8N_c)$ has $Psi[Y]=1$;
MS-bar fields do not require it to vanish. The source-bound gate
`native_d09_emitted_vacuum_policy_field_matches_actual_source_rows`
on `79b315b05cc0e786198426e8d44774b35bc6c4a9` passes in run
`af45fbc2-c8e6-4066-9f40-f7dc2deb462f` (216.405 seconds). At unit rest, with
unit physical, UV and renormalization scales in the coefficients, it binds
both actual production $A_D$ matrices to
their dimensional W rows and reconstructs $A_4$ independently from the actual
physical-four/$C_0$ rows. With $r=11/36+pi^2/12$, the complete source sum is

$ Delta I=Omega(7/(6epsilon)+r)Y, quad A_D-A_4=Omega r Y. $

Omega uses the amputated Grassmann-kernel source weight. The finalizer retains
the independently source-derived external field-ordering factor $q_e=-1$;
this is separate from the forest sign. Within each history, both policies
reuse the same projected local payload.
A separate exact supplement compares the bare and completed-child histories:
their frozen parent localizer K, tree factor t and source ordering are
literally equal. The observed pointwise difference is
$Delta F=q_e K t (A_D-A_4)$, with the original K and common t retained. Thus
$Phi[Delta F]=0$ and $Psi[Delta F]=-K t Omega r$. The full matched-rest
mass identity $Phi[R]=0$ remains the separately qualified complete-zero claim.

Exact registered-binding bijections, their inverses, typed RHS/port checks and
complete dispatch-root/catalogue equality cancel the shared local emissions.
The remaining full-key carrier sum retains $q_e K t$ once. A doubled actual
nonzero mapped local carrier, a zeroed field coefficient and replacing Psi
by Phi must fail; the source epsilon-depth controls remain active. Neither
carrier is defined by the final reference, and no extra production R2 term
is inserted. This qualifies the emitted unit-rest/unit-scale UV-vacuum policy
correction, not the full physical unregulated two-loop R2, a general-scale
identity or a UV/IR integration exchange. Psi remains a native test closure.

The fresh local-only UV matrix resolves all 259 spanning-tree-LMB/nonempty
hard-subset limits at 25 scales from $10^8$ to $10^12$. Every subtracted fit
passes; matched bare domains show 122 UV failures, four unstable fits and
133 finite limits, with a resolved divergence in each fixture. Eighteen
subtracted and 126 bare ArbPrec retries remain recorded. Complete orientation
sums are tested; representation-local checks retain their small separate gates.
The hierarchy checks 48 groups, 240 finite nonzero ArbPrec values and 96 tail
inequalities for fixed signed rest/nonaxial endpoints. Worst ratios 0.514505
and 0.252150 meet bounds 0.75 and 0.375. These finite grids do not prove a
uniform dimensional or angular limit or a UV/IR integration exchange.

The retained public run on `3dc81ff7` takes 101.066 seconds, UUID
`bbd7513f-9153-408e-8a79-de702fcabd4a`, for three routes, two loop points and
five repeats. All values are finite and nonzero; maximum recorded relative
route difference is `1.826e-15`. Its completed production clocks are:
#table(columns: (1.8fr, 1fr, 1fr, 1fr),
  table.header([Duration (s)], [Projected 4D], [Explicit direct], [Orientation direct]),
  [4D CT and scheme], [5.21742], [5.12271], [4.93814],
  [Direct Taylor / finite prefix], [No calls], [1.11669 / 0.19016], [1.07866 / 0.16486],
  [Raw CFF requests], [0.0165560], [0.00335726], [0.00304752],
  [Exact CFF preparation / requests], [0.0886278 / 1.98562], [0.000699892 / 0.00933384], [0.000748433 / 0.00978161],
  [CFF surface conversion], [0.00217806], [0.000191281], [0.000199432],
  [Evaluator construction], [1.30314], [10.47935], [11.86201],
)
Scopes are inclusive and can overlap. Medians of ten evaluator summary means
are 30.1215, 60.3942 and 635.110 µs. Factorized roots plus definitions occupy
232637, 452594 and 471417 bytes; sampled RSS is 94236672, 251195392 and
310034432 bytes including prior route allocations. These are correctness-run
diagnostics, not isolated peaks or a paired performance gate. Source and
run identity remain `3dc81ff7` even though its exact inputs qualify reuse.


=== Latest FeynKit replay and current-source qualification
<latest-feynkit-replay-qualification>

The earlier replay onto FeynKit `58814885` preserved ninety ordered copies
and its eighty-two-change OS prefix. Physics source
`79b315b05cc0e786198426e8d44774b35bc6c4a9`, producer 45, first closed its
48-case ledger, then the full 69-case selected qualification set below.
Every admitted run binds all 2702
immutable source inputs and its retained executable; no earlier-source runtime
is counted as a fresh result. The earlier 39-case ledger and its timings above
retain their original source and run identities. The new ledger includes the
actual emitted complementary-field certificate and the physical NLO gates.

The inclusive $gamma^* arrow.r t overline(t)$ NLO gate passes in run
`07de4b07-bf9c-4026-a746-e24208c3477a` (2431.954 seconds). At centre-of-mass
energy 400 GeV, pole mass and renormalization scale 173 GeV, its three routes
are projected 4D, explicit direct 3D and orientation direct 3D. Each route
samples the UV/localization scales $(M,L)=(173,173),(346,173),(173,346)$ GeV,
giving nine actual NLO cells. Ten raw result files retain one 10000-evaluation
LO result and nine 320000-evaluation NLO results: 2.88 million NLO evaluations,
with ten workers. A separate 420-digit consumer recomputes the pole-mass reference
and the original five-percent error and three-sigma comparisons from the
persisted estimates. Maximum reference, auxiliary-scale and route distances
are respectively 1.3991591, 1.2716299 and 1.4529333 sigma; all pass. This is
statistical agreement of the complete summed and integrated NLO observable,
not an algebraic identity. Adjacent deterministic seeds and iteration offsets
reuse some random streams; independent streams are not asserted.

Five further full-main route-matrix cases pass freshly on the same source
in run `09cd113f-163e-4e2f-ad2c-200e5fe07ccb` (3463.259 seconds): a quark OS
child in a gluon self-energy, a nested quark self-energy, a soft-refined gluon
child, disjoint OS components and the three-loop nested quark fixture. The
three routes at two momentum points and five repeats retain 150 finite
nonzero values. All 100 cross-route comparisons pass the stricter bound
$|a-b|<=1000 epsilon_(64) max(|a|,|b|)$; the worst relative difference is
$1.007387 times 10^(-13)$, or 0.453687 of this bound. Repeats agree literally.
These local-only comparisons sum orientations; the disjoint fixture is an
off-shell Green-function scalar diagnostic. Numerical agreement does not
replace the separate exact reconstruction certificates.

Six mandatory soft-counterterm cases additionally pass freshly on the same
source in run `8b75f683-694b-4ec8-8d68-bb0dafa35910` (7792.200 suite seconds).
Every case retains its own 7200-second cap; the slowest is 6894.169 seconds.
The independent consumer admits all 45 retained states, 36 UV profiles with
all spanning-tree LMBs/nonempty hard subsets and fifteen matched control
groups, plus nineteen soft reports. The full and child-only top-bubble
comparisons give ordinary-U slopes about -2.000198799 and H slopes about
0.000390166: all six summed route comparisons improve by more than 2.0005889
powers on the identical routed rays. The two orientation-local inventories
each retain thirty orientations and exactly eight resolved divergent-U
controls, all removed by H. Original slope, fit-quality and scale policies
are unchanged; the retained positive fit subsets remain distinct from the
twenty-five requested UV scales. These are finite-profile and source-bound
forest checks, not a uniform soft or integration-exchange theorem.

Ten further full-main native core cases pass in
`5f427b53-03f7-43aa-9413-24bd90fd61c9` (8059.841 suite seconds).
Each case retains its original 7200-second cap; the slowest is 3616.801 seconds.
The closed consumer checks 58 finite nonzero soft-ray groups, 290 ArbPrec
evaluations and 116 weighted tail inequalities; the largest ratio is
1.0123844 below 1.5. Independent rational geometry recovers all twenty
generic affine flats (six of rank one, ten of rank two and four of rank three)
and 200 off-plane sample coordinates.
The native exact hard tests cover all fourteen nontrivial leading-carrier
partitions, with fifteen completed inverse-hard coefficient certificates,
the eight ordered-forest component jets and the open logarithmic coefficient.
Three proper-parent/finite-prefix cases additionally preserve the live
cograph across three routes and two external points. Numerical component
values not serialized by this suite remain source/JUnit-bound; sampled soft
rays and exact local UV certificates do not imply a simultaneous CFF majorant.

The frozen-source ledger covers 69 distinct fresh selected qualification cases
(48 + 5 + 6 + 10), not the entire repository main suite. On that source,
the explicit nested telemetry case timed out at its unchanged 240-second cap;
the orientation case was not run after that failure. These remain historical
results on `79b315b0`. All six ordinary route telemetry cases pass freshly
on qualified combined source `c8eade36`; its nested explicit and
orientation cases take 131.398 and 126.365 suite seconds. Correctness-run
clocks under shared-host load establish no quiet paired performance gate.

Two extensions are qualified on separate immutable sources. The new native
pointwise child endpoint test on `8588b9af` passes in
`dc99fda6-6b52-4e51-a780-2eb8305b9f4d` (6.384 suite seconds); five existing
references pass in `7ae8d197-9c9b-42fd-84a1-6f8030da2c23` (6.458 seconds).
The test proves twelve exact LEFT $P_sigma$ zeros for the actual bare plus
once-signed local child, frozen analytic child and complete child at both
signs and rest/nonaxial subtraction velocities. Physical mass and point stay
tied, the actual sub-LMB and consumed $S^4$ remain owned, and the bare O branch
is an exact nonzero control. This four-dimensional endpoint gate supplies
neither an unintegrated right-zero identity nor a simultaneous causal-surface
majorant or integration-limit exchange.

The named shared component-collection memoization fix `d2375ca5` passes
thirteen native cases in `d26d8e89-a246-410b-b3ce-2bcc61d8995d` (0.918
seconds), together with strict Clippy, formatting and the existing docs
builder checks. Its exact catalogue/inverse-binding and signed hyperdual
controls preserve six logical axes, formal and alias roles and foreign scalar
heads. The bounded unit case records 64 requests, 62 cache hits and two unique
inputs, versus 64 requests and no hits for eager collection. These are
functional cache observations, not OS performance measurements.

Qualified combined source `c8eade36` retains eighty-six ordered OS copies
above that shared prefix. Its five-file delta from `79b315b0` contains
the byte-identical qualified cache and endpoint proofs plus their documentation
registration. Forty selected cases execute freshly from independently retained
binaries: nineteen cache/endpoint cases, six ordinary telemetry cases, the
public integrated three-route rainbow, three ordinary rest/boost route checks,
four physical-small cases, six exhaustive UV cases and one hierarchical soft
case. This is a selected qualification set, not the entire main suite.

The UV run `2dcf2fe7-33ad-41a8-bb3b-6c4c53b60518` (555.724 seconds)
covers all 259 spanning-tree-LMB/nonempty-hard-subset domains, with matched
bare controls and 25 scales per domain. Primitive and soft-gluon fixtures
use explicit direct; the other four use projected 4D. The hierarchy run
`12c1736c-1a98-4a80-9478-2b5b7dd03a10` (127.362 seconds) checks
48 groups, 240 finite nonzero ArbPrec values and 96 tail inequalities.
Its largest ratios remain 0.514505 and 0.252150 below 0.75 and 0.375.
These finite scaling tests establish neither uniform angular/epsilon bounds
nor an integration-limit exchange. Earlier `79b315b0`, `8588b9af` and
`d2375ca5` results retain their separate identities.


=== Six-route construction, evaluation and memory diagnostics

The six ordinary telemetry cases on `c8eade36` all pass the unchanged
240-second cap. Primitive rows include local and integrated counterterms;
nested rows use the Appendix-B1 OS-child graph with local counterterms only.
P4 denotes projected 4D, ED explicit direct and OD orientation direct.
The CT column records raw local-4D/direct-3D construction calls. Raw CFF,
exact preparation/lifting, assignment and surface-conversion scopes are
recorded separately. Every phase is inclusive and can overlap another:
the columns must not be summed. A dash means no completed observed call.

#table(columns: (1.3fr, 0.85fr, 0.85fr, 0.85fr, 0.85fr, 0.85fr, 0.9fr),
  table.header([Fixture / route], [Local CT (s)], [Raw CFF (ms)],
    [Exact prep. (ms)], [Assignment (ms)], [Surfaces (ms)], [Evaluator (s)]),
  [Primitive P4], [0.0230175], [0.767602], [2.24379], [26.9557], [0.098431], [0.0305905],
  [Primitive ED], [0.0497085], [0.491882], [—], [—], [0.04552], [0.0457731],
  [Primitive OD], [0.0413434], [0.387422], [—], [—], [0.03617], [0.0463836],
  [Nested P4], [1.98762], [10.9265], [32.2076], [876.628], [2.04225], [2.09325],
  [Nested ED], [2.09194], [3.33874], [1.33365], [16.2789], [0.33613], [115.813],
  [Nested OD], [2.16424], [1.3839], [0.690452], [15.4038], [0.090183], [111.305],
)

Evaluation columns below are arithmetic means of five repeat means at each
of two fixed momentum points. Evaluator time is separate from the Total
category of the minimal-integrand benchmark. Expression sizes retain
factorized roots and function definitions; the parenthesized count is the
number of definitions. Memory columns are
actual test-process observations: sampled generation VmRSS and separately
reported lifetime VmHWM, in MiB. The sample monitor can miss short peaks and
the VmHWM capture includes earlier process activity. Two rows have sampled
VmRSS above the independently reported VmHWM; both original values are kept.
GNU-time launcher-chain RSS is not used as a guaranteed test-process or
simultaneous process-tree peak.

#table(columns: (1.25fr, 1.05fr, 1.05fr, 1.65fr, 0.9fr, 0.9fr),
  table.header([Fixture / route], [Evaluator p0 / p1 (µs)],
    [Benchmark Total p0 / p1 (µs)], [Root / definition bytes (count)],
    [Sampled VmRSS (MiB)], [Reported VmHWM (MiB)]),
  [Primitive P4], [1.58425 / 1.59032], [10.5958 / 10.6143], [1554 / 5739 (7)], [51.176], [53.785],
  [Primitive ED], [1.67203 / 1.67525], [10.4646 / 10.4764], [4025 / 5233 (11)], [53.570], [55.996],
  [Primitive OD], [3.039 / 3.01143], [12.1257 / 11.9781], [4165 / 5233 (11)], [54.164], [56.809],
  [Nested P4], [33.2359 / 32.0745], [47.4969 / 45.7302], [180809 / 96366 (28)], [274.133], [273.434],
  [Nested ED], [91.3237 / 91.5843], [106.071 / 106.408], [1005163 / 78909 (75)], [2654.793], [2677.043],
  [Nested OD], [234.867 / 233.262], [254.187 / 252.612], [1094064 / 78909 (75)], [2719.715], [2713.641],
)

The six cases retain sixty comparisons with archived source-79 point values.
Actual nested direct collection records 125013 requests, 108296 hits and
16717 unique inputs for each direct representation. These observations do
not establish a quiet paired speedup, a route ranking or an isolated peak-RAM
result. The generation, evaluator, evaluation and memory scopes remain
separate; the integrated public rainbow has its own independent case.


=== Direct dimensional D09 vector reference

The separate native vector reference `c9532785` passes five ordinary cases
in run `8e78ccac-adc3-4365-a796-0428a7173eb9` (3.981 suite seconds).
The direct $slash(p)/(4N_c p^2)$ trace agrees with an independently authored
ordered Clifford contraction before any mass-reference comparison. At the
unit external shell, with k the outer quark momentum and l the inner one,
the bare scalar numerator is

$ n_A=(2-d)^2[2(k dot l)(p dot k)+(m^2-k^2)(p dot l)]
      +2d(2-d)m^2(p dot k). $

The raw trace keeps internal MT symbolic; the certified finite scalar export
then sets physical masses to one. The exact source phase is
$+i G^4 C_F^2 n_A$, and the common $i^(L+1)$ bridge retains
$Gamma_2=-i Sigma$. The generated external field-ordering sign remains a
separate factor. Vector/identity and vector/kinetic Gram checks give zero and
one; early d-to-four and replacing the vector by the mass trace are detected.
An algebraic kinetic-shift witness uses the actual source's declared logical
frame with its checked inverse; it is not a second full source execution.

All five original propagator occurrences retain their momenta, masses and
unit powers. The current inner/outer chart $(k_0,k_1)$ and archived
$(k_1,k_2)$ chart have an exact simultaneous inverse, including scalar-product
metadata. The scalar family certifies zero widths and absent custom
denominators; the common Feynman $+i 0$ prescription is explicitly imposed at
the independent backend boundary. It is not extracted from a
prescription-free scalar family.

The unchanged mass branch exports all thirteen sources byte for byte:
23477 bytes, SHA256 `e7ce7a27c02d203386e2a69a809fc033b88c384239fa609b5ed010e55ffab5bd`.
The independently known proper-insertion vector reference retains epsilon
one before multiplication and gives
$5/epsilon^2+12/epsilon+32+11pi^2/12$
in units of $a^2 C_F^2$, with $a=G^2/(16pi^2)$, at unit physical,
auxiliary and renormalization scales. The actual once-signed dimensional W
kinetic rows independently give the parent counterterm
$a^2 C_F^2[1/(2epsilon^2)-1/(4epsilon)]$ after the same amputated
phase conversion.

One separate bare-vector backend integral now closes on `c9532785`.
The guarded backend run takes 47.319 seconds. Adding the exact proper
reference and the once-signed UV poles gives residuals
$2.4885 times 10^(-10)$ and $3.79630087 times 10^(-7)$ at epsilon powers -2
and -1; each meets the unchanged componentwise absolute bound $10^(-6)$.
The prepared finite kinetic coefficient at $m=mu=1$, in units of $a^2 C_F^2$,
is
$ A_("prepared")^(0)=-19.0743405536087186
   plus.minus 6.3862627 times 10^(-7). $
The uncertainty is the backend's reported bare-integral error; proper terms
are analytic. The raw epsilon-minus-three output is reported as zero, which
is not promoted to an exact algebraic identity. The fixed absolute pole
tolerance is separate from the reported integration error. These numerical
reference data do not establish ordinary four-dimensional CFF absolute
integrability, Fubini or a full unregulated physical R2 condition.
Combined source `b762b717` passes seven ordinary cases in run
`167a1ddc-8ad9-4a3a-98fa-43ade5a77575` (10.100 suite seconds), with fresh
Clippy, formatting and producer checks. Its exact backend inputs match the
qualified standalone vector source, so the numerical reference retains the
preceding integral's source identity.

The maintained optional `vector_source` path is qualified separately at
`edefbb2c`, named `on-shell-maintained-d09-vector-reference-20261020`.
Seven ordinary cases pass in `a4088c53-8544-4a0a-886e-a68d4c03096d`
(10.012 seconds), together with strict Clippy, Rust formatting and Ruff.
The existing main-only reference passes in
`ad69d70c-fe9d-4e99-96b3-81c2046d5321` (505.995 seconds), under its
unchanged 7200-second limit. All thirteen mass sources, nine original backend
inputs and 47 comparison centres and bounds remain unchanged. The optional
vector contributes one further integral and four coefficient checks. Malformed
source, projector, dimensional numerator, occurrence, ordering and prescription
controls reject before integration; separately perturbed output coefficients
fail their comparisons.

The fresh prepared finite value is $-19.07434047342197$ with reported bare
error $7.5664505 times 10^(-7)$. The preceding reference error and their linear
sum $1.3952713 times 10^(-6)$ are recorded separately, without a statistical
independence claim or a change to the centre-based acceptance tolerance.
The fresh signed pole residuals at epsilon powers -2 and -1 are
$2.4884983 times 10^(-10)$ and $-5.3733297 times 10^(-7)$; both satisfy
the unchanged $10^(-6)$ zero-reference bound. The epsilon-minus-three backend
output remains a numerical zero. The source, retained binary, all 3538 tracked
input hashes and embedded reference owner are bound to this run. The finite
reference uses the same pinned backend and policy as the retained datum; its
independence is from production counterterm construction. Integration of the
complete production OS forest and the physical four-dimensional R2 condition
remain separate validation boundaries.

=== Independent signed centers in the boosted rainbow

The original 48 hierarchical groups put both loop spatial momenta near the
signed live external point. Distinct component branches also require
independent centers. Keep $p=s m v$ and $v$ fixed, and define

$ bold(r)=bold(k)-sigma m bold(v), quad
  bold(w)=bold(l)-tau m bold(v), quad
  bold(q)=(tau-sigma)m bold(v)+bold(w)-bold(r). $

Here $s,sigma,tau in {+1,-1}$, and q is the spatial normal of the physical
inner-gluon soft plane.
For the nonaxial velocity and mass twelve, the mixed-center offset has squared
spatial norm 324. Thus q remains hard when both r and w tend to zero at mixed
centers. At rest the spatial centers coincide. The coordinate certificates
retain the original edge routing, parent and child owners, physical mass and
fixed velocity; only these small coordinate polynomials are expanded.

The extra finite-grid contract covers $sigma=tau=-s$ and both mixed choices
$(sigma,tau)=(s,-s),(-s,s)$ for each signed live external point. Coincident
centers retain all three hierarchies; mixed centers use only the outer and fixed-child
hierarchies. With exponents $a in {2,3}$, the radial measure is
$rho^(3a+2) dif rho$. The coincident-center candidate has margin $a-1$.
For mixed branches, the adjacent ordered projectors
$P_sigma P_tau=0$ remove the leading doubled-outer-quark coefficient of the
child O branch. The regular-in-w part additionally requires the complete
child identity $P_sigma R_c(sigma m v,l)=0$, including its frozen finite
contribution. A residue derivative of the outer numerator retains only the
weaker outer pole; differentiating the denominator retains this endpoint
zero. Together these give the candidate $|r|^(-2)|w|^(-2)$ and margin a.
The projectors remain on their declared child row throughout. These are
predeclared finite-ray hypotheses for the complete forest sum. A uniform
angular CFF remainder and an epsilon-uniform integration bound require
separate proofs.

The expanded numerical grid passes on `92d15824` in run
`8171d04c-c21e-4a6a-a592-2a190c37cb3d` (248.339 suite seconds). Its
unchanged original 48 groups and 56 additional groups give 520 finite nonzero
ArbPrec evaluations and 208 independently recomputed tail inequalities, using
two cached fixed-velocity evaluators. Exact rational arithmetic independently
recovers the centers and physical-gluon offsets at all 280 additional sample
coordinates. These numerical results retain this source and run identity;
qualification of the separate symbolic-coordinate assertion is recorded
independently.

Corrected source `5b6a297a` qualifies that symbolic certificate using exact
expansion of the finite coordinate polynomials. Eight ordinary cases pass in
`4310c67b-dfef-4384-b53e-bc484d35f322` (9.976 seconds), together with
strict Clippy and formatting. Its fresh full-main scan passes in
`267dca30-d643-493d-9aed-26d56d4486fd` (245.550 seconds), retaining the
same 104 groups, 520 finite nonzero evaluations, 208 tail comparisons and
280 independent rational-coordinate checks. All 3538 tracked inputs and
the retained executable stay unchanged throughout each run. The original
48-group source and all predeclared bounds remain intact. The larger scan
uses the existing full-main limit; ordinary cases retain their 240-second cap.

Fresh mass and vector exports are byte-identical to the qualified `edefbb2c`
exports, and both Rust reference owners and the embedded Python owner are
unchanged. The nine mass and one vector integrations therefore keep their
original `edefbb2c` run identity. The final union adds local signed-center and
exact geometry coverage, with no additional production counterterm.
