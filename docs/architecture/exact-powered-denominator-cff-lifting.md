# Exact powered-denominator lifting and factorized EMR assignment

## Status and scope

This document describes the production design used when a completed local 4D
term is projected to a causal-flow representation (CFF) after UV Taylor
operations have raised one or more propagators. It expands two tightly coupled
requirements:

1. carry each original `source_edge` through the UV Taylor operator so the
   source minors provide the exact term's skeleton and provenance, then lift
   every rewritten denominator power into an owner-independent rational
   occurrence graph; and
2. rewrite the still-factorized numerator over those exact occurrences while
   minimizing the largest energy-degree bound supplied to generalized CFF.

This exact-source reconstruction and minimax-dispatch machinery is exclusive to
the projected local-4D route. Both direct local-3D representations instead
complete the loop-energy integration first and apply every UV Taylor operator
to the complete/global CFF expression. Their bodies are identical; the
explicit-sum form only omits the localized form's orientation selectors.

The two requirements cannot be solved independently. A denominator power first
creates the occurrence-local energy namespace in which the generalized residue
problem is posed. The numerator then has to use exactly that namespace, and the
same factor-to-occurrence assignment must control both CFF generation and
numerical numerator evaluation.

This is GammaLoop production machinery. It is not a requirement that the
diagnostic `3Drep` CLI prepare its inputs or expressions in the same way. It is
also deliberately independent of LTD; no LMB variable is ever used as an
energy-identity or energy-capacity index here.

## Equality is defined for the complete fixed-cut functional

The direct local-3D and projected local-4D constructions need not agree on an
individual generalized residue-map entry or an individual raised-order slot.
Their authoritative equality is

\[
  \mathcal R_C^{\mathrm{direct}}=\mathcal R_C^{\mathrm{projected}}
\]

for each fixed physical Cutkosky cut \(C\), after both sides have summed all
generalized residue maps and combined all derivative pieces from every raised
order required by that cut. The number, labeling, and individual values of the
intermediate entries may differ. They are coordinates of a residue
decomposition, not observables.

For example, a double-pole result can be represented schematically by

\[
  \mathcal R_C=\left.\partial_\eta f_2(\eta)\right|_{\eta=0}+f_1(0).
\]

For any regular \(h\), the change

\[
  f_2\mapsto f_2+\eta h,
  \qquad
  f_1\mapsto f_1-h
\]

(with the second \(h\) evaluated at the expansion point) changes the two
stored pieces but leaves \(\mathcal R_C\) exactly invariant. Higher raised
powers have the corresponding triangular freedom between adjacent derivative
slots. A test which pairs `lu_cut_order` values across two construction routes
therefore imposes a nonphysical canonicalization unless that particular
fixture independently proves the stronger identity. The safe comparison sums
first, with the Cutkosky cut held fixed.

This does not weaken orientation-local UV subtraction in the direct-3D route.
Each direct residue-key selector must still give its required local UV
behavior. It only means that projected local-4D, which is defined after an
explicit residue sum, has no canonical one-to-one key or raised-order pairing
with that selector-local decomposition.

### Factorized staged map composition

Projection integrates disconnected or nested Taylor components in stages. A
state is retained as a CFF carrier times a still-factorized numerator, rather
than as one expanded polynomial. At component \(a\), its exact source map acts
only on the numerator factors owned by \(a\); the resulting carrier/numerator
pair is passed unchanged to component \(a+1\). After all UV components have
been integrated, the outer CFF map acts once on the still-unmapped soft and
cograph factors. Symbolically, the implementation realizes products of mapped
factors and sums of CFF carriers without expanding the full numerator.

Genuinely independent unsummed component residue states form the mathematical
Cartesian product of their sums. A child coefficient which has already summed
its source-local residues may be replicated under several compatible
production hosts solely so those hosts can supply the later outer map; those
replicas are not independent states and must be consumed once per outer host,
not cross-multiplied again. This distinction preserves both the complete
fixed-cut functional and the factorized numerator representation.

## Executive summary

For a term

\[
  I(q)=\frac{N(q)}{\prod_j D_j(q)^{r_j}},
  \qquad
  D_j(q)=P_j(q)^2-m_j^2,
\]

the production projection performs the following operations:

1. A negative integer power `den(...)^-r` becomes `r` denominator
   occurrences. The numerator is not expanded.
2. Every occurrence retains its `source_edge`. Contracted cograph and UV source
   minors use those IDs to provide the already known skeleton, topology domain,
   and external attachment. Matching momentum signatures never solve an
   inverse incidence or Kirchhoff problem.
3. The `r` occurrences of a rewritten denominator are represented by a minimal
   occurrence graph. A raised source wrapper is realized by serial subdivision,
   with `r-1` auxiliary two-valent vertices; a powered self-loop becomes an
   `r`-edge cycle. The final rational incidence and loop rank are canonical in
   the denominator algebra, rather than in the arbitrary physical owner labels.
4. Exact signatures and masses are canonicalized into algebraic channels. This
   is the safe `D(Q)=D(-Q)` equivalence used by repeated-channel algebra. Source
   minors seed the construction, but relabeling algebraically identical
   occurrences with another compatible physical owner cannot change the CFF
   residue or loop rank.
5. Physical owner IDs survive separately for occurrence-energy provenance,
   physical-surface projection, diagnostics, and mapping back to the parent
   graph. Cut support is the union of all physical owners which instantiate an
   algebraic channel, with raised-line representatives added without discarding
   that union.
6. The numerator is analyzed in physical EMR variables without expansion.
   Addition takes a maximum degree, multiplication and multilinear slots add
   degrees, and a nonnegative integer power repeats a factorized base.
7. Exact denominator occurrences which carry the same algebraic on-shell
   energy as one physical EMR edge form a certified candidate set.
8. Each numerator factor is assigned to one candidate so as to minimize the
   largest per-occurrence energy degree. Ties are resolved by canonical
   occurrence order.
9. One immutable assignment plan owns both the bounds passed to generalized CFF
   and the substitutions used later to evaluate the numerator.
10. If an occurrence has momentum `-Q` rather than `Q`, that literal algebraic
    sign is retained by the numerator mapper. It is composed with any later
    coherent reversal of the complete rational-routing component by comparing
    the original exact signature with the final parsed occurrence signature.
    The resulting sign converts the occurrence energy back to the physical
    numerator convention; denominator evenness never causes the numerator to
    be treated as even.
11. Once incidence, power subdivisions, and non-vacuum source-crown hedges are
    complete, Graphica canonically relabels the skeleton with internal
    propagators treated as undirected and external crown incidences left
    directed. Connected rational-routing components then receive one
    deterministic coherent direction. This supplies deterministic equality and
    cache keys without inferring topology.
12. All additive Taylor terms are planned before CFF generation. Terms with the
    same canonical topology share one channel-wise maximal capacity envelope,
    while each term retains its own immutable minimax numerator assignment.

The resulting graph is therefore a source-backed occurrence representation of
the rewritten rational function. Provenance is not discarded: it supplies the
source-minor scaffold and all maps back to the physical graph. It is also not an
extra argument of the rational function. Once algebraically identical
occurrences have been lifted, changing only their compatible owner labels must
not change their residue, loop rank, or owner-free CFF cache identity.

## Terminology and identity layers

Several IDs coexist in this path. Conflating signature equivalence with physical
incidence, or occurrence-local capacity with physical EMR ownership, produces
incorrect topology or numerator sampling.

| Concept | Meaning | May determine rational topology? | May own a numerator bound? |
| --- | --- | --- | --- |
| Physical graph edge | An edge of the original amplitude or forward-scattering graph | It supplies a known source-minor attachment; its label is not part of the owner-free rational identity | Yes, during the first physical-EMR analysis |
| Provenance or source owner | The original edge recorded in a `den(edge, momentum, mass, value)` wrapper | It selects the source-minor scaffold/domain and physical projections, but cannot distinguish algebraically identical rational functions | No, by itself |
| Denominator occurrence | One copy produced by a negative denominator power | Its algebraic channel and multiplicity determine the canonical rational occurrence graph on the source scaffold | Yes, after certified physical-to-exact lifting |
| Exact CFF energy ID | The occurrence-local OSE/EMR label used by the temporary exact source | It labels an already source-constructed occurrence | Yes |
| Physical EMR energy | The temporal component of `Q(edge, ...)` in the parent graph | No; this is numerator vocabulary | Yes |
| LMB coordinate | A coordinate used to route loop momenta | No; it only expresses exact momentum signatures and active coordinates | Never |

An occurrence-local exact energy is not a new physical momentum. If a double
propagator is represented by exact energies `E_a` and `E_b`, both can be
certified representations of the same physical energy `E_5`. They remain
distinct variables only because a generalized residue calculation needs one
energy slot per denominator occurrence.

## Why UV processing creates this problem

The unexpanded graph normally has one denominator wrapper per propagator. UV
Taylor operations differentiate and expand rational factors. Even if the input
topology contains only simple propagators, a completed local 4D term can contain

```text
den(owner, P, m2, D(P))^-2
```

or a product such as

```text
den(owner_a, P,  m2, D(P))^-1
den(owner_b, -P, m2, D(P))^-2.
```

At that point, using the original graph unchanged is insufficient: the energy
integral has a double or triple pole, whereas an ordinary graph edge encodes
only one denominator occurrence. It is wrong to discard the source edge and
try to reconstruct a graph from the rewritten signatures: the UV/cograph
minors and their attachment to the unshrunk graph are already known from the
original graph. It is equally wrong to promote that retained label into
algebra. Two occurrences with the same domain, mass, and denominator
`D(Q)=D(-Q)` represent the same rational channel even when they came from
different physical edges. A compatible owner relabeling may change physical
provenance, but cannot change the contour residue or the number of active loops.

The numerator poses a second problem. UV differentiation can simultaneously
produce higher powers of temporal momenta. For example, the rational term may
contain

\[
  \frac{(Q^0+c)^2}{D(Q)^2}
  \quad\text{or}\quad
  \frac{(Q^0+c_1)(Q^0+c_2)(Q^0+c_3)}{D(Q)^2}.
\]

After creating two occurrence energies, the numerator can be evaluated through
either occurrence. Choosing just the first is algebraically possible on the
physical diagonal, but it supplies unnecessarily high generalized-CFF bounds
and can make CFF generation much more expensive. Worse, advertising balanced
bounds while evaluating all factors through the first occurrence would be
incorrect. This is why structural lifting and numerator assignment are one
coupled operation.

## Production data flow

```text
completed local 4D UV term
  |
  | factorized term projection
  v
numerator Atom + Vec<FourDDenominator>
  |
  | source-minor lookup plus algebraic-channel validation
  v
known UV/cograph skeleton and occurrence channels keyed by
  (topology domain, signature up to sign, mass^2)
  |
  | source-backed attachment + owner-independent occurrence lift
  v
source-backed occurrence graph
  |
  | explicit source-crown completion for non-vacuum boundaries
  | post-construction Graphica canonical relabeling
  v
canonical temporary exact ParsedGraph
  |                         physical factorized numerator
  |                                      |
  |                                      | EMR degree analysis
  |                                      v
  |                         physical edge -> physical degree
  |                                      |
  | exact literal +/-Q candidates        |
  +--------------------+-----------------+
                       |
                       | deterministic minimax planning
                       v
            exact occurrence -> energy bound
            factor -> exact occurrence assignment
                       |
             +---------+---------+
             |                   |
             v                   v
      generalized CFF       exact numerator mapper
      generation            used by every orientation,
                            residue, and contact sector
```

The assignment plan is created once, before generalized CFF generation. It is
then retained with the exact-source numerator mapper. This is the certification
boundary: the expression which is eventually sampled is the expression whose
bounds were supplied to CFF.

Projection is a two-pass batch. The first pass constructs every canonical exact
source and joins the required degrees of equal topologies. The join is performed
in canonical occurrence coordinates and by algebraically repeated energy
channel: it takes the maximum total degree requested on that channel and
redistributes that total by the same quotient/remainder minimax rule. The second
pass generates one CFF expression per topology with that envelope and maps each
term with its original plan. Thus batching can increase available capacity but
can neither change a term's factor assignment nor understate its rank.

## Non-negotiable common-LMB commutative reconstruction invariant

The projected local-4D route is not required to preserve an individual CFF
orientation. It is therefore free to choose any convenient EMR spelling of a
completed Taylor numerator on the reconstructed factorized UV graph. That
freedom does **not** make the rewrite heuristic. Every chosen spelling must
pass the following exact commutative check before its generalized-CFF image is
trusted.

Let the Taylor operator be defined in the compatible hard sub-LMB
\(\ell=(\ell_1,\ldots,\ell_L)\). Provenance attached before differentiation
retains, for every original edge \(e\), both its immutable owner and its hard
momentum

\[
  H_e(\ell)=\sum_i A_{ei}\ell_i.
\]

After taking the requested Taylor coefficient, let its still-factorized
numerator be \(N_T(\{H_e\})\). Build the UV skeleton from the retained original
owners; never infer its incidence from a signature matrix. A differentiated
denominator may add serial copies of its own line, but no other operation
changes the source skeleton. Choose a deterministic LMB
\(\kappa=(\kappa_1,\ldots,\kappa_L)\) for this reconstructed graph and write
each of its EMRs as

\[
  \widehat Q_a(\kappa,p_{\mathrm{fixed}})
    =\sum_i B_{ai}\kappa_i+C_a p_{\mathrm{fixed}}.
\]

For an MSbar vacuum counterterm the fixed-external part is absent,
\(C_a=0\). A future on-shell counterterm may retain a fixed insertion such as
\((m,0,0,0)\); it is simply included on both sides of the same check.

The mapper proposes a factorized UV-graph numerator
\(\widehat N(\{\widehat Q_a\})\). Correctness means

\[
  \boxed{
  N_T\bigl(\{H_e(\kappa,p_{\mathrm{fixed}})\}\bigr)
  =
  \widehat N\bigl(\{\widehat Q_a(\kappa,p_{\mathrm{fixed}})\}\bigr)
  }
\]

as an exact symbolic identity in one common set of formal loop four-vectors,
masses, and fixed external data. In a test or diagnostic assertion, provenance
tags may be erased and the difference expanded or put over a common
denominator to prove that it is exactly zero. That expansion is an oracle only:
the production numerator on both sides stays factorized.

This check has four important consequences:

1. The original numerator is never minimax-dispatched. Its retained owner and
   hard-momentum payload determine an exact UV-EMR representative directly.
2. Only new numerator energy factors produced by differentiating a denominator
   may be assigned among the serial copies created for that same physical line.
3. Such an assignment is correct before it is optimal. A deterministic first
   valid copy is sufficient; minimax redistribution is allowed only after the
   boxed identity passes and only among those degenerate copies.
4. Canonicalizing \(D(Q)\) as \(D(-Q)\) cannot discard an odd-numerator sign.
   The complete signed coordinate map from the tagged hard momentum to the
   selected UV occurrence participates in the boxed equality. Denominator
   evenness changes the denominator channel, never the numerator transformation.

The LMB has a deliberately limited role here. It is a common coordinate chart
used to prove that two EMR expressions denote the same function. It is not an
energy identity, an occurrence owner, a CFF bound label, or a fallback when
provenance is missing. Once the exact identity is certified, only the candidate
UV-graph EMR expression and its physical/occurrence metadata are passed to
generalized CFF.

The denominator side has a parallel certificate: after substituting the same
common coordinates, the reconstructed occurrence multiset must reproduce the
Taylor term's denominator momenta, masses, powers, and topology domains exactly
(up to the safe even relation \(D(Q)=D(-Q)\)). Numerator equality and
denominator equality are checked separately so denominator evenness cannot
hide a numerator sign error.

This commutative identity is the first diagnostic boundary for every projected
local-4D mismatch. If it fails, the defect is in provenance preservation,
owner-backed skeleton construction, occurrence assignment, or signed EMR
rewriting. If it passes but the residue-summed result differs from direct
local-3D UV, the defect is downstream—in CFF input normalization, prefactor or
component composition, or residue aggregation. Rank-envelope optimization must
never be used to repair a failed identity.

### Worked diagnostic: GL04 temporal-square `T0` sector

The minimal GL04 acceptance probe localizes the original numerator on edge 5,

\[
  N_0=i g^5\,(Q_5^0)^2.
\]

With integrated and threshold counterterms disabled, the first differing
forest contains the Taylor-vacuum components labelled `1zs` and `8Ok`.
Certifying its derivative-free outer `8Ok/T0` sector first cleanly separates
identity reconstruction from the derivative-bearing `1zs/T2` sector. The
typed local-4D `T0` sector entering projection reports exactly

```text
component 1zs:
  owner 5, momentum  Q5
  owner 6, momentum -Q5

component 8Ok:
  owner 2, momentum  Q2
  owner 1, momentum -Q2
  owner 7, momentum -Q2

numerator:
  i*g^5*Q(provenance(owner=5, role=fixed, hard=Q5), 0)^2
```

Every occurrence has the Taylor-vacuum mass and no external shift. Introduce
common formal four-vector carriers

\[
  q\equiv Q_5,\qquad r\equiv Q_2,
  \qquad D(k)=k^2-M_{\mathrm{UV}}^2.
\]

Erasing only the test copy of the provenance metadata, the authoritative
post-T rational function becomes

\[
  I_T(q,r)=
  \frac{i g^5(q^0)^2}
       {D(q)D(-q)D(r)D(-r)D(-r)}
  =\frac{i g^5(q^0)^2}{D(q)^2D(r)^3}.
\]

The owner-built reconstructed UV graph has the same five occurrences and the
same deterministic component LMBs. Because `T0` generated no new numerator
factor, the only admissible candidate is the unchanged fixed-owner spelling

\[
  \widehat N=i g^5(\widehat Q_5^0)^2,
  \qquad \widehat Q_5=q.
\]

Substituting the reconstructed graph LMB therefore gives

\[
  \widehat I(q,r)=
  \frac{i g^5(q^0)^2}
       {D(q)D(-q)D(r)D(-r)D(-r)},
  \qquad I_T-\widehat I=0
\]

exactly. There is no minimax choice in this comparison: the original fixed
factor stays on owner 5, and all five denominator powers are one. At the inner
CFF residue, the actual source map sends

\[
  Q_5^0\mapsto \pm E_5,
  \qquad E_5=\sqrt{M_{\mathrm{UV}}^2+\boldsymbol q_5^2},
\]

so both ordinary pole branches map the numerator to \(E_5^2\), while the
zero-sampling contact maps it to zero. The outer component leaves that already
sampled factor untouched.

An independent diagnostic using the production generalized-CFF engine then
compared the disconnected five-denominator source with the product of its
separately embedded three- and two-denominator components. The combined source
had
\(6\times5=30\) branches—six scalar three-denominator orientations times the
five quadratic two-denominator key families `contact:+`, `contact:-`,
`contact:0`, `remainder:+`, and `remainder:-`. Its total agreed with the
separately embedded component product at the tested points to approximately
\(10^{-16}\).

This worked sector therefore passes the pre-CFF common-LMB certificate, the
generalized-CFF component-product diagnostic, and a static trace of the
GammaLoop host/prefactor/outer-factor coefficient. It is not the first unequal
boundary in the enclosing forest. The next mandatory comparison is therefore
the derivative-bearing `1zs/T2` sector. Changing this fixed numerator
assignment would violate the exact certificate and cannot be a valid repair.

### Worked live reproducer: GL04 temporal-square `1zs/T2`

The full-graph generation LMB used by the minimal reproducer contains

\[
\begin{aligned}
 Q_3 &= -K_2-P_0+K_0, & Q_5 &= K_1,\\
 Q_6 &= -K_1-K_2-P_0+K_0=-Q_5+Q_3.
\end{aligned}
\]

For the `1zs` Taylor component, set

\[
 q=Q_5,\qquad p=Q_3,\qquad
 A=D_5(q),\qquad B=D_6(-q),\qquad
 x=p\mathbin{\cdot}(-q),\qquad U=M_{\mathrm{UV,exp}}^2.
\]

The completed order-two Taylor trace has the exact denominator occurrence
multiset

```text
owner 5: +q, +q       (production occurrence ids 9 and 12)
owner 6: -q, -q, -q   (production occurrence ids 10, 11, and 13)
```

and the still-factorized numerator is algebraically

\[
 N_T=-g^2(q^0)^2
 \left[
   -A B^2 +(U+p^2)A B-4x^2A+U B^2+2xA B
 \right].
\]

This is the authoritative left-hand side.  It is obtained before any CFF call,
and no numerator/denominator cancellation is performed or needed.  The
original factor \((q^0)^2\) retains `TaylorFixed(owner=5, hard=+q)` provenance.
The occurrences beyond the first copy of each owner and every factor generated
by differentiating a denominator retain `DenominatorDerived` provenance.

To perform the actual reconstruction certificate, introduce *neutral proof
coordinates* rather than treating occurrence ids as independent physical
energies.  Let every parsed occurrence have a formal momentum $P_j$.  The
source records two signs for each selected occurrence:

\[
 H=hR,\qquad P_j=rR,\qquad H^0=hrP_j^0,
\]

where $H$ is the immutable hard momentum, $R$ is the raw rewritten
denominator momentum, and $P_j$ is the canonical parsed momentum.  In this
reproducer the selected owner-5 occurrences have $h=r=+1$, while the
owner-6 occurrences have $H=R=-q$, $P_j=-R=+q$, hence $h=+1,r=-1$.
Consequently the production mapper must give

\[
\begin{aligned}
 (Q_5^0)^2_{\rm fixed} &\longmapsto (P_9^0)^2,\\
 D_5(q) &\longmapsto (P_{12}^0)^2-\boldsymbol P_{12}^{2}-U,\\
 D_6(-q) &\longmapsto (P_j^0)^2-\boldsymbol P_j^{2}-U,
     &&j\in\{10,11,13\},\\
 p\mathbin{\cdot}(-q) &\longmapsto p\mathbin{\cdot}(-P_j)
     =-p^0P_j^0+\boldsymbol p\mathbin{\cdot}\boldsymbol P_j.
\end{aligned}
\]

Only now impose the common reconstructed-graph LMB chart,

\[
 P_9=P_{10}=P_{11}=P_{12}=P_{13}=q.
\]

Every mapped $D_5$ becomes $A$, every mapped $D_6$ becomes $B$, every
mapped dot becomes $x$, and the candidate numerator becomes

\[
 \widehat N\big|_{P_j=q}=-g^2(q^0)^2
 \left[-A B^2 +(U+p^2)A B-4x^2A+U B^2+2xA B\right]
 =N_T.
\]

That last zero-difference identity is the required numerator certificate.  The
parallel denominator certificate is the multiset $A^2B^3$, with the common
Taylor-vacuum mass $U$, the same component domain, and momentum equality
using only $D(q)=D(-q)$.  The LMB is used only to establish these identities;
it does not define occurrence ownership or CFF rank capacity.

The actual immutable production plan has now been enumerated term by term.  Set

\[
 d_j=(P_j^0)^2-\boldsymbol P_j^{2}-U,
 \qquad s_j=p\mathbin{\cdot}(-P_j).
\]

It makes the following deterministic choices:

```text
original fixed Q5^0 * Q5^0:                   9, 9
every positive owner-5 D5 wrapper:            12
the D6*D6 factors in the cubic term:           10, 11
the single D6 in the (U+p^2) term:             10
both p·(-Q5) factors in the quadratic-dot term: 10, 10
the p·(-Q5), D5, D6 factors in the linear term: 10, 12, 11
```

Thus the exact expression that GammaLoop asks the mapper to sample is

\[
 \widehat N=
 -g^2(P_9^0)^2\left[
   -d_{12}d_{10}d_{11}
   +(U+p^2)d_{12}d_{10}
   -4s_{10}^2d_{12}
   +U d_{10}d_{11}
   +2s_{10}d_{12}d_{11}
 \right].
\]

After applying the common-LMB diagonal $P_9=P_{10}=P_{11}=P_{12}=q$,

\[
 d_9=d_{10}=d_{11}=d_{12}=d,
 \qquad s_{10}=x,
\]

and therefore

\[
 \widehat N-N_T=0
\]

exactly.  This is the real post-plan certificate: it includes every selected
occurrence and every signed hard/raw/parsed conversion, rather than comparing
an earlier pre-plan expression.  The parallel denominator certificate is also
unchanged and exact.  Values of the formal (P_j) away from the common-LMB
diagonal are not an alternative physical routing and are not used to judge
the reconstruction identity.

The executable regression
`gl04_t2_planned_lift_matches_post_t_numerator_in_common_loop_coordinates`
freezes this production trace and runs both proofs independently.  It builds
the actual owner-5/owner-6 candidate sets, asks the production
`EnergyPowerAnalyzer` for the immutable plan, maps the factorized numerator
with `ExactSourceEnergyMapper`, and only then erases tags and expands a
test-only copy.  It first keeps $P_9^0,P_{10}^0,P_{11}^0,P_{12}^0,P_{13}^0$
distinct and proves that the mapper produced the documented term-by-term
assignment.  It then imposes the common chart and proves the numerator identity
above.  Its final exact-zero assertion compares the source multiset
`D(+q)^2 D(-q)^3` with the five reconstructed occurrences.  The focused test
passed on 2026-09-03.  This is a frozen executable certificate of the live
trace, not a replacement for the enclosing three-route acceptance: the latter
is still needed to detect downstream CFF composition or evaluator defects.

### Reconstruction correctness does not determine the CFF consumer boundary

The two certificates above end before GammaLoop converts a generated CFF into
its historical evaluator convention. This separation matters: a perfectly
reconstructed numerator can still acquire an overall wrong sign at that later
boundary.

Generalized variants retain occurrence-local positive half-edge factors

\[
  \prod_{j\in c}\frac{1}{2E_j},
\]

while the ordinary GammaLoop CFF path removes those factors from each variant
and appends the source-global convention

\[
  \prod_{j\in c}\frac{1}{-2E_j}.
\]

For one independently generated exact rational component `c`, the fresh-source
conversion is therefore

\[
  (-1)^{N_c}\,B_{\mathrm{den},c}\,B_{\mathrm{core},c},
\]

where `N_c` is the number of source denominator occurrences, `B_den` is the
scalar-denominator frame recorded by the generator, and `B_core` cancels the
generated causal core's uniform convention. This conversion is
component-local but **generation-context invariant**. Marking a source as an
`EmbeddedCffFactor` may select one equivalent terminal residue; it cannot
change the value or sign of the same rational energy integral.

This is not the contract for an already persisted production root. The stored
root already represents the complete production residue functional; changing
between global and variant-local energy-factor ownership only changes where
those factors are stored. Materializing that root consumes `B_core` once and
must not apply `(-1)^N B_den` again. Thus the two adapters are deliberately
separate:

\[
  B_{\rm fresh\ exact}=\prod_c(-1)^{N_c}B_{{\rm den},c}B_{{\rm core},c},
  \qquad
  B_{\rm stored\ root}=B_{\rm core}.
\]

The minimal clean-tree oracle is an uncancelled, still-factorized quotient

\[
  \frac{D(Q)(Q^0+c)}{D(Q)^3}
  =\frac{Q^0+c}{D(Q)^2}.
\]

At `|q|=0`, `E=1`, and `c=2`, the historical core-only fresh-source adapter at
`c8e763173` returned exact Arb values `+2.015720902074968...e-3 i` and
`-2.015720902074968...e-3 i`, respectively. The raw generalized-CFF identity
was already correct; only the production values were opposite. The typed
component conversion above makes this oracle and its scalar, reversed-routing,
quintic-to-quartic, and raised-residue extensions agree. This is decisive
evidence that the gross sign is not caused by the EMR reconstruction and does
not justify modifying generalized-CFF residue recursion.

The complementary stored-root oracle uses a scalar self-energy `T0` child and
the factorized outer numerator `Q3^0 Q4^0`. The child reconstructed from 4D is
exactly `-1/[D(q)D(-q)]`, its generated CFF agrees with the analytic
lower-contour integral of `-1/D(q)^2`, and the outer LU energy map is coherent.
The root metadata are `N=5`, `B_den=+1`, `B_core=-1`. Reusing the fresh-source
formula at the stored-root boundary changes `-1` to `+1` and flips the complete
direct result. Consuming only `B_core` makes direct post-CFF Taylor,
projected-child-times-outer, and whole exact-source calculations agree for the
complete selected-cut functional. In this deliberately canonical `T0` oracle,
the separately stored LU-order pieces also happen to agree bit-for-bit at
1000-bit Arb precision, including the first Taylor derivative. That stronger
fixture-specific result must not be generalized into a per-order route
contract. Since the generated residues, numerator map, and projected assembly
are held fixed in this one-line A/B, the gross sign is a GammaLoop
double-conversion defect, not a `generalized_3drep` contact sign.

The plan supplies bounds `(9,2),(10,2),(11,2),(12,2)`.  Keeping each positive
`GS.den` wrapper wholly on one same-owner occurrence is a conservative lifting
convention, but the stable residue-summed A/B remains numerically unchanged by
it.  A component-only exact audit subsequently proved that all five occurrences
form one repeated channel with one mass key and that its eight `+/- n`,
`n=1,...,4`, generalized-CFF families reproduce the analytic `d^-2`, `d^-3`,
and `d^-4` lower-sector moments term by term.  Both standalone and
`EmbeddedCffFactor` generation give the same `3/8` total at the rational audit
point.  The first unequal boundary is therefore later still: attachment and
mapping of this certified child CFF into the outer factor, final GammaLoop
aggregation/evaluation, or the corresponding direct-side Taylor term.  A
former diagnostic which compared numerators before applying the plan was
removed: it could pass without testing the object actually sampled by CFF, and
its unconditional symbolic expansion violated the factorized production
boundary.

A trace of that next boundary has excluded three further false leads.  The
three observed `216 -> 60` projection waves are the three distinct physical LU
Cutkosky-cut invocations, not three copies of one Taylor sector.  With
integrated counterterms disabled, the first-level `1zs` sector has no frozen
LMB and its frozen localizer is exactly one.  Finally, on a concrete selector-0
branch the outer source map replaces only unwrapped cograph/soft energies.  It
leaves all 24 owner-5 and all 32 owner-6 provenance-wrapped child atoms
unchanged, so the already sampled child numerator is not being sampled a
second time.  The 216 child entries are production-host selector ids carrying
the same already-residue-summed child coefficient and no child source map; they
are not 216 independent residue families.  Matching them to the 60 outer source
maps broadcasts exactly one complete child coefficient into every surviving
outer host, so the `216 -> 60` boundary neither drops nor multiplies a child
residue.

The smallest unresolved comparison is consequently one physical cut: sum its
projected `outer CFF * mapped child T2` branches and compare that exact
expression with the sum of its direct post-CFF `T2` branches.  Raw branch counts
cannot be paired.  Every projected wave materializes 60 Cartesian host records,
of which 40 are exactly zero for physical cut `[1,2]` and 48 are exactly zero
for each of `[3,4]` and `[1,7,8]`; the direct lane instead prunes and coalesces
its generalized residue keys to 56, 48, and 32 records.  For cut `[1,2]`, the 20
nonzero projected records split further into 14 ordinary `Q3^0 = +/- E3`
source-map branches and six pure `Q3^0 = 0` contact branches.  Here `lu_cut_1`
denotes a first-order Cutkosky residue, not physical cut id 1.  The only valid
next oracle is therefore equality of that complete 20-term projected outer
functional and the complete 56-term direct functional after exact summation,
not equality of individual keys or branch counts.

## Stage 1: retain denominator multiplicity without expanding the numerator

The local-4D term projector traverses the Symbolica expression structurally.
For a negative integer power of a denominator wrapper, it repeats the
`FourDDenominator` record `abs(power)` times. Non-denominator powers remain in
the numerator. The first argument of every typed `GS.den` wrapper is the
original `EdgeIndex`; the UV Taylor operator differentiates and rewrites the
other arguments while carrying that owner unchanged. Consequently, the
completed Taylor coefficient still has a direct lookup into the original graph
and never has to recover an edge from the rewritten momentum.

Taylor collection can also expose a positive typed `GS.den` factor inside a
numerator multiplying a negative occurrence of the same wrapper. Projection
does not cancel those factors or distribute that nested numerator sum. It
splits only an outer additive Taylor result whose addends need separate
denominator topologies. For example,

```text
((A+B)*den_0 + den_1) / (den_0*den_1)
```

is projected with numerator `(A+B)*den_0 + den_1` and denominator occurrences
`[den_0,den_1]`. The numerator remains one factorized atom, and generalized CFF
performs any resulting pinches internally. By contrast, a completed outer sum
`T_0/den_0 + T_1/den_1` is parsed as two terms because its addends genuinely
carry different denominator topologies. The focused regression
`term_projection_keeps_typed_numerator_factors_uncancelled` protects that
distinction.

For example,

```text
first^-2 * second^-1 * N
```

becomes

```text
numerator:    N
denominators: [first, first, second]
```

The focused test
`term_projection_preserves_dots_and_distinct_same_edge_expressions` uses two
wrappers with the same physical edge and momentum but different retained
`full_expr` provenance. It verifies that the first wrapper occurs twice and the
second once. The companion test
`term_projection_keeps_factorized_numerator_atoms` verifies that an expression
such as

\[
  \bigl(A(k_0)+A(k_1)\bigr)^2(k_0+k_1)
\]

survives term projection exactly as a factorized atom.

This step supports arbitrary negative integer powers. It does not contain a
special case for a double pole or for a self-energy graph.

## Stage 2: recover the source skeleton and classify algebraic channels

Two keys are retained here for different purposes. The physical `source_edge`
locates an occurrence in the already known cograph or UV source minor. It tells
GammaLoop which original component, attachment, and topology domain produced
the factor, so no graph has to be inferred from a matrix of rewritten momentum
signatures. Separately, GammaLoop computes an exact momentum signature in the
temporary source coordinates and canonicalizes it up to an overall sign. The
algebraic-channel key is

```text
(
    topology_domain,
    canonical_momentum_signature_up_to_sign,
    canonical_mass_squared,
)
```

The retained `full_expr` payload, input order, and physical owner label do not
change this channel. The owner remains the authoritative link back to the
source skeleton and physical graph, but it is not a variable of the rewritten
rational function.

This gives three separate invariants:

- every exact occurrence can be attached to the appropriate original source
  minor without an algebraic graph-reconstruction search;
- repeated copies of one source wrapper form a powered line only when they
  represent the same rewritten denominator up to `D(Q)=D(-Q)`; and
- after that source-backed construction, the rational residue and loop rank of
  an algebraic channel are invariant under a compatible relabeling of its
  physical owners.

A compatible relabeling changes only the `source_edge` provenance among
occurrences which instantiate the same normalized denominator in the same
topology domain and for which the source-backed lift remains valid. It does not
move a factor between UV and cograph domains, alter a mass or affine external
shift, or change a non-vacuum boundary. Those changes describe a different
rational problem rather than a provenance relabeling.

If one source edge appears with incompatible normalized signatures or masses in
one additive 4D term, generation fails instead of choosing an incidence by
heuristic. UV-vacuum and cograph denominators are constructed in disjoint source
minor domains even when their algebraic data happen to look alike.

### Provenance is retained, not discarded

Owner provenance remains on every `FourDDenominator`. Across construction and
the later projection back from the owner-free `ParsedGraph`, it is used to:

- select the source-minor component and attachment which seed the exact
  skeleton;
- recognize repeated copies produced by raising one physical wrapper;
- project occurrence-local energy IDs back to physical edge IDs;
- retain distinct occurrence-energy provenance;
- build Cutkosky and LU cut support as the union of physical owners in each
  algebraic channel;
- recognize which exact surfaces are original physical surfaces;
- report and diagnose the source of a UV-expanded factor;
- recover parent-graph spatial momenta and affine external shifts.

The rule is therefore:

> Use source provenance to recover the known skeleton and all physical maps;
> use normalized denominator algebra to define the owner-independent rational
> occurrence graph on that skeleton.

## Stage 3: instantiate the source-minor scaffold and lift powers

The original graph is never reconstructed. Instead, each occurrence first
receives a deterministic attachment in the appropriate source minor.

For cograph denominators, a disjoint-set union over the original graph nodes
contracts source edges absent from the additive exact term. Each surviving
owner identifies the endpoints of its original `EdgeIndex` in that quotient.
UV denominators use a separate disjoint-set union for the UV source minor:
absent UV edges are contracted there, and every active UV owner inherits its
endpoints from that minor. The two node domains remain distinct.

There is no Kirchhoff solver, signature-matroid reconstruction, or
signature-to-endpoint fallback. A rewritten momentum such as `k1+k2` still has
a `source_edge`; its source-minor attachment identifies where that factor came
from. Exact signatures are subsequently checked against the scaffold. If the
source-backed construction and supplied momenta cannot realize a valid exact
occurrence graph, generation fails explicitly.

For a raised source wrapper, once its base attachment is known, multiplicity
`r` is subdivided as

```text
tail --occ_0-- x_1 --occ_1-- ... --occ_(r-2)-- x_(r-1) --occ_(r-1)-- head
```

This introduces precisely `r-1` auxiliary vertices. If `tail == head`, the
result is an `r`-edge cycle. Every occurrence keeps its own exact energy ID. Its
`ParsedGraph` signature is normalized up to sign, and every segment is directed
coherently.

That source-backed subdivision is a construction device, not an
owner-sensitive definition of the integral. Equal algebraic occurrence
multisets are subsequently canonicalized without owner labels. Thus two copies
of `D(Q)` must yield the same one-loop, double-pole CFF whether their provenance
is `(edge 0, edge 1)` or `(edge 0, edge 0)`, provided both owner assignments are
compatible with the same source skeleton and topology domain. The source IDs
avoid a combinatorial endpoint search; they do not grant two spellings of the
same rational function different residues.

Here, “owner-independent” does **not** mean that GammaLoop throws away the
source graph and reconstructs another graph from canonical signatures. The only
topological operations are contraction of known source-minor edges, local
subdivision for denominator multiplicity, and canonical relabeling of the
result. Owner invariance is a contract on that deterministic lift. There is no
enumeration of candidate incidences, Kirchhoff-matrix inversion, or
signature-matroid search anywhere in this path.

The UV Taylor operator can remove a soft/cograph shift from a denominator. For
example, an original source edge may carry `Q_hard+Q_soft`, while its rewritten
UV denominator carries only `Q_hard`. Its scaffold attachment still comes
directly from the original edge. To determine only whether that attachment is
used forward or backward, GammaLoop tests both signs modulo the span of the
known opposite topology domain:

```text
exact_row - sign * source_row in span(opposite_domain_rows), sign in {+1,-1}.
```

Exactly one sign must pass. This rank calculation is a routing-coordinate
validation after both graph minors are already fixed; it never chooses nodes,
constructs a graph, or repairs momentum balance. The focused hard-plus-soft
triangle regression checks that the correct sign leaves only a soft residual
and that the opposite sign raises the quotient-space rank.

Serial subdivision is the graph-theoretic representation of a dotted line: the
new two-valent vertices do not introduce independent momentum constraints. Raw
occurrence signs and physical owners are retained separately by the exact
numerator/provenance maps rather than encoded as alternating segment directions.

### Non-vacuum boundaries use explicit source-crown hedges

After internal incidence is fixed, a non-vacuum exact signature can leave an
external-momentum imbalance at a source-backed node. GammaLoop represents that
boundary explicitly with source-crown hedges carrying the corresponding
external coefficients. This completes the already known source topology; it
does not use external or loop signatures to infer internal endpoints. Contracted
components containing no exact denominator or cut carrier are multiplicative
identities and do not manufacture causal surfaces.

The implemented non-vacuum path currently wires pure-external source crowns:
the relevant parent crown hedge has no parent-loop coordinate and carries its
known external coefficients into the exact source. This is already sufficient
for ordinary non-vacuum amplitude boundaries. A future on-shell UV scheme may
instead leave a two-point exact topology with a fixed insertion such as
`(m,0,0,0)`. That extension should add the fixed four-vector as an explicit
boundary payload and wire it to the same source-derived endpoints. It does not
change the internal disjoint-set construction and does not require an incidence
or Kirchhoff reconstruction.

### Determinism

Occurrence input is first put in a stable order so raised-wrapper subdivisions
are reproducible. After all source attachments, auxiliary power vertices, cut
carriers, and crown hedges exist, Graphica canonically relabels the edge-coloured
graph. Internal propagator incidences are undirected during this relabeling
because `D(Q)=D(-Q)`; external crown incidences remain directed. In the
canonical node namespace, exact occurrences are grouped only for routing by
connected component and the key

```text
(UV/cograph domain, canonical signature up to sign, mass, power marker).
```

For each such component, GammaLoop compares the sorted multiset of directed
`(tail,head,signature,...)` edge keys with the multiset obtained by reversing
every incidence and negating every complete momentum signature. It retains the
lexicographically smaller whole-component spelling. Reversing the component as
one unit is essential: a same-owner dotted chain, or equal rational
denominators carried by successive distinct source edges, must not acquire
alternating segment directions. Disconnected equal denominators remain
separate routing components, and owner IDs remain provenance rather than an
algebraic routing key.

Exact edges are then sorted in the canonical node namespace and the
occurrence-to-original-factor map is reordered with them. Source-node and
power-node names survive as aliases of the canonical nodes; cache keys
deliberately clear those names and edge labels. Reversing the order of Symbolica
factors therefore does not change the parsed exact graph or its CFF cache entry.
Crucially, this pass only relabels and orients an already source-constructed
graph: it never derives an endpoint from a momentum signature.

### Rational incidence and physical support cross the boundary separately

The owner-free exact CFF and the maps back to the physical graph deliberately
cross separate interfaces:

- **Energy provenance is occurrence-local.** Each canonical exact energy maps
  back to the physical owner carried by its original denominator occurrence.
  Distinct owners therefore remain distinguishable to numerator evaluation and
  diagnostics even when their denominators are algebraically identical.
- **Cut support is channel-wide.** For each algebraic denominator channel,
  GammaLoop collects the sorted, deduplicated union of every physical
  `source_edge` which instantiated it. Any exact occurrence used as residue
  support is projected to that full union. This prevents an arbitrary canonical
  occurrence representative from erasing a valid physical Cutkosky carrier.
- **Raised-line compatibility is additive.** Downstream LU selection may name a
  canonical representative of a physical raised-edge group. GammaLoop adds
  those representative aliases to the already unioned physical support and
  merges equal support-sign entries. It does not replace the original owner
  set. Thus existing raised-cut matching continues to work without losing
  provenance.

The channel used for this neutral cut-support union includes the topology
domain, mass, and exact signature canonicalized up to sign. Consequently
`D(Q)` and `D(-Q)` share support, but equal-looking UV and cograph factors, or
factors of different mass, do not. None of these support operations modifies
the rational incidence, exact energy map, or numerator routing sign.

## Stage 4: analyze the factorized numerator in physical EMR variables

The numerator is analyzed before any occurrence-local substitution. The
analysis works solely in physical `Q(edge, index)` variables.

For each physical edge `e`, it computes a conservative exact polynomial degree
`d_e` using the following composition rules:

| Expression | Degree rule |
| --- | --- |
| Constant, parameter, or spatial-only value | zero |
| `Q(e,0)` or the temporal component of `Q(e,mu)` | one in edge `e` |
| Sum | componentwise maximum |
| Product | componentwise sum |
| Nonnegative integer power | multiply the base degree by the exponent |
| Dot product | add the degrees of the two vector slots |
| Declared multilinear function | add the degrees of its argument slots |
| `den(..., full_expr)` retained positively in a numerator | analyze `full_expr`, not provenance metadata |

These rules avoid algebraic expansion. In particular,

\[
  (Q^0+c)^r

\]

is traversed as `r` factor slots containing the same base. It is not expanded
into `r+1` monomials. Additive branches reuse capacity because the degree of a
sum is a maximum, not a sum.

Opaque nonlinear functions of EMR energies and negative energy-dependent
powers are rejected. Production does not guess a bound for an expression it
cannot certify.

## Stage 5: certify equivalent exact occurrence energies

For every physical EMR edge which has nonzero numerator degree, the exact
source mapper finds denominator occurrences whose rewritten momentum is
literally

\[
  +Q_e \quad\text{or}\quad -Q_e.
\]

The exact on-shell energy of every candidate is then computed after source
replacement. Candidates are partitioned into classes of equal algebraic
on-shell energy. The largest class is selected; equal-size ties are resolved by
the lowest canonical occurrence ID.

Candidate sets for different physical edges must be disjoint. A missing,
empty, duplicated, or overlapping candidate set is an error. In particular,
GammaLoop does not treat an arbitrary shifted momentum `Q_e+p` as an alias for
`Q_e`, and it does not use an LMB coordinate as a fallback.

The literal `+/-Q_e` restriction is important for sign safety. It supplies both
an algebraic energy identity and an exact routing sign. Equality of positive
on-shell energies alone would not be sufficient to map an odd numerator.

## Stage 6: compute the deterministic minimax assignment

Suppose physical edge `e` has total numerator degree `d` and `n` certified
equivalent exact occurrences. Any valid assignment has nonnegative loads

\[
  \ell_0,\ldots,\ell_{n-1},
  \qquad
  \sum_i \ell_i=d.
\]

The smallest possible maximum load is

\[
  \min \max_i \ell_i = \left\lceil\frac{d}{n}\right\rceil.
\]

The implementation obtains this optimum by quotient/remainder distribution:

\[
  q=\left\lfloor\frac{d}{n}\right\rfloor,
  \qquad r=d\bmod n,
\]

and assigns load `q+1` to the first `r` canonical candidates and `q` to the
rest. This is both minimax-optimal and deterministic.

Examples are:

| Physical degree | Equivalent occurrences | Exact loads |
| ---: | ---: | --- |
| 2 | 2 | `(1,1)` |
| 3 | 2 | `(2,1)` |
| 4 | 2 | `(2,2)` |
| 5 | 3 | `(2,2,1)` |
| 7 | 3 | `(3,2,2)` |

The planner assigns individual factor slots, not merely total degrees. For

\[
  (Q^0+c_1)(Q^0+c_2)(Q^0+c_3)

\]

over two equivalent occurrences `a,b`, the canonical plan is conceptually

\[
  (E_a+c_1)(E_a+c_2)(E_b+c_3),

\]

with bounds `(2,1)`. Another balanced assignment would be algebraically valid,
but canonical ordering makes generation reproducible.

For

\[
  (Q^0+c)^2,

\]

the internal planned representation repeats the unexpanded base and maps it as

\[
  (E_a+c)(E_b+c).

\]

The original expanded polynomial is never constructed.

## Stage 7: use one immutable plan for generation and evaluation

The assignment plan contains both:

- `exact occurrence -> certified degree bound`; and
- `factor-local physical edge -> chosen exact occurrence`.

The first map is passed to the shared generalized-CFF generator. The second map
is retained by the exact-source numerator evaluator. Each orientation, LU
residue, threshold residue, and lower-sector contact term maps the numerator
using the same assignments.

This prevents a subtle but serious inconsistency. It would be wrong to balance
the advertised bound as `(1,1)` but later replace every physical `Q^0` by
`E_a`. The CFF expression would then have been generated for a lower rank than
the expression it actually samples. Because one immutable plan owns both
operations, that state is not representable.

Sampled directions also survive lower-sector pinching. If a contact sector
assigns a zero sample to the occurrence which owns a factor, that exact zero
map remains authoritative even if an equivalent dotted occurrence survives.
The numerator is not silently reassigned to the surviving occurrence after the
pinch.

## Stage 8: batch canonical topologies without weakening numerator bounds

Completed Taylor coefficients are prepared in two passes. The first pass builds
each source-backed exact graph, computes that term's immutable minimax plan, and
registers its requested capacity under a key consisting of the complete
canonical `ParsedGraph` and all non-bound generation options. The second pass
generates CFF once for each key and evaluates every term with its own original
plan.

Capacity is joined in the algebraic coordinates understood by generalized CFF.
For a repeated on-shell-energy channel, the relevant invariant is the **total**
degree across all equivalent occurrences. The batch therefore takes the
largest channel total requested by any term and redistributes only that total
over the canonical occurrences with quotient/remainder minimax balancing. It
does not take a pointwise maximum which could inflate `(1,0)` and `(0,1)` into
the unnecessary channel rank `(1,1)`. Bounds on non-repeated occurrences use
the ordinary componentwise maximum. The shared envelope is generator capacity
only: it never replaces a term's factor-to-occurrence assignment.

The production-shaped regression
`dod_one_triangle_keeps_one_uncancelled_exact_source_and_matches_lower_sectors`
applies the real degree-one UV Taylor operator to the three-edge UV triangle in
the double-triangle skeleton. The collected coefficient remains one
common-denominator exact source with owner multiplicities

```text
(2,1,2).
```

The positive typed denominator factors which would algebraically expose the
base and singly dotted lower sectors remain in the additive factorized
numerator. No parser-side distribution or cancellation occurs. The source uses
one canonical generator-cache entry, and the regression verifies that its
exact CFF orientation sum equals a test-only oracle in which those lower sectors
are reduced explicitly. The fixed cograph component is part of the canonical
key, so this reuse does not conflate different surrounding graphs or future
non-vacuum boundaries.

## Why `D(Q)=D(-Q)` cannot introduce numerator sign mistakes

This deserves a separate proof because the denominator and numerator have
different parity properties.

Let the literal exact denominator carry momentum

\[
  P_{\mathrm{lit}}=s_{\mathrm{lit}}Q,
  \qquad s_{\mathrm{lit}}\in\{+1,-1\}.
\]

For a scalar propagator,

\[
  D(P_{\mathrm{lit}})=P_{\mathrm{lit}}^2-m^2
  =(s_{\mathrm{lit}}Q)^2-m^2=D(Q).
\]

It is therefore correct to give `+Q` and `-Q` the same normalized denominator
signature and algebraic channel. If they are repeated occurrences of one source
wrapper, they validate as copies of that powered line. If they belong to
distinct source edges, those owners remain in physical provenance and unioned
cut support, while the rational occurrence graph remains owner-invariant.
Neither case reconstructs topology from the sign.

It does **not** replace `Q` by `-Q` in the numerator. After owner-free graph
canonicalization, a complete rational-routing component may also have been
reversed coherently. Write its final parsed occurrence momentum as

\[
  P_{\mathrm{parsed}}=s_{\mathrm{route}}P_{\mathrm{lit}},
  \qquad s_{\mathrm{route}}\in\{+1,-1\}.
\]

The mapper determines `s_route` by comparing the final parsed signature with
the original exact signature; it must not reuse only the earlier
`canonical_up_to_sign` result because the later component canonicalization may
have reversed that direction again. The occurrence-local CFF energy map is the
temporal component `p_parsed^0`, so a physical numerator written in terms of
`Q` is reconstructed as

\[
  q^0=s_{\mathrm{lit}}s_{\mathrm{route}}p_{\mathrm{parsed}}^0,
\]

because both signs are their own inverses. The certified occurrence sign stored
by the mapper is therefore the product
`s_literal * s_route`; an explicit negation is performed exactly when that
product is `-1`.

Consequently,

\[
  Q^0+c
  \longmapsto
  \begin{cases}
    p_{\mathrm{parsed}}^0+c,
      & s_{\mathrm{lit}}s_{\mathrm{route}}=+1,\\
    -p_{\mathrm{parsed}}^0+c,
      & s_{\mathrm{lit}}s_{\mathrm{route}}=-1.
  \end{cases}
\]

An odd numerator therefore retains its sign. A vector numerator is handled in
the same way: spatial components remain in the physical parent-graph EMR
basis, while the corrected temporal energy is inserted with the temporal unit
vector. The mapper never performs a blanket substitution `Q -> -Q` on the
physical numerator.

There are thus two separate routing layers:

1. The parsed graph uses one coherent canonical direction for each completed
   rational-routing component. Its initial attachment came from the source
   minor, but owner labels are absent from the owner-free rational identity.
2. The numerator mapper retains the literal sign, derives the final component
   routing sign from the parsed occurrence, and composes both with the
   occurrence-local temporal map to recover the physical `Q` convention.

The first layer is even denominator algebra on a source-backed skeleton; the
second is signed physical-numerator mapping. Denominator evenness affects only
signature equivalence and never silently turns an odd numerator into an even
one. Likewise, unioning `+Q` and `-Q` cut support is a neutral provenance
operation and never rewrites a numerator factor.

### A focused sign oracle

The test `exact_energy_bounds_balance_across_routing_signs_of_one_energy_channel`
uses three equivalent occurrences with sign patterns `(+,-,-)` and `(+,+,-)`.
For the first pattern, a degree-two physical numerator is assigned to the first
two occurrences, one positive and one negative. The CFF maps supplied to the
mapper are correspondingly `(+E,-E)`. The explicit routing correction produces

\[
  (+E)\,[-(-E)] = E^2.
\]

Without the numerator-side sign restoration, the same test would produce
`-E^2`; the even final power does not hide the error because only one assigned
factor uses the reversed occurrence.

The end-to-end exact-CFF identity test also retains an explicitly odd factor
`Q^0+c` while cancelling one denominator. It verifies at two rational phase-
space points that

\[
  \frac{D(Q)(Q^0+c)}{D(S)D(Q)^2}
  =
  \frac{Q^0+c}{D(S)D(Q)}

\]

matches both the lower exact source and the ordinary CFF source. This protects
the odd numerator mapping independently of a purely even `Q^2` probe.

## Complete example 1: a same-owner cubic powered line

The unit fixture `exact_source_serializes_dotted_same_edge_occurrences` starts
from the physical graph

```text
        edge 0
    a ----------> b -> c -> a
                         |
                         +----> d
```

and constructs three exact denominator records for physical edge 0. Two records
carry one retained `full_expr` label and the third another, but all have the
same rewritten momentum and mass.

The exact source contains three internal occurrence edges and two auxiliary
power-chain nodes. Because the base incidence closes onto itself after source
contraction, the three segments form a three-edge cycle. The test checks:

- all three occurrences survive;
- every cycle node has one incoming and one outgoing segment, independent of
  canonical edge-vector order;
- exactly two power-chain nodes exist;
- vertex momentum balance has no violation;
- occurrence-local energy IDs remain distinct; and
- all three occurrence IDs project back to physical edge 0.

This is the simplest complete realization of

\[
  D(Q)^{-3}

\]

as a graph accepted by the shared CFF recursion.

## Complete example 2: owner relabeling preserves the rational CFF

The regression
`exact_source_owner_relabeling_preserves_residue_rank_and_cut_provenance`
uses a physical graph with two parallel edges:

```text
        edge 0
    a ----------> b
        edge 1
    a ----------> b
```

It constructs two spellings of the same rewritten denominator multiset:

```text
distinct provenance: (owner 0, owner 1) -> D(Q)^-2
relabeled provenance: (owner 0, owner 0) -> D(Q)^-2
```

The source graph and retained IDs make both constructions immediate: no
algorithm searches momentum signatures for a two-node bubble or solves a
Kirchhoff incidence problem. Their normalized signatures and masses identify
the same algebraic double channel. The owner-free exact graphs are both valid,
have one topological loop, and give the same CFF residue at exact rational
on-shell points after the neutral physical-energy aliases are identified.

What changes is only the physical projection. The distinct spelling retains
energy provenance `{0,1}` and cut support `{0,1}`. The relabeled spelling
retains energy provenance `{0}` and cut support `{0}`. In the first case, each
exact occurrence projects to the union `{0,1}` for cut eligibility, even though
its occurrence-local energy still maps to the one owner which produced it.

This fixture protects the precise boundary: source IDs recover the known
skeleton and survive as provenance, but owner identity cannot control the exact
rational incidence, loop rank, or contour residue. It also demonstrates why
cut support must be unioned after canonicalization: choosing one canonical
occurrence must not hide the other physical owner.

## Complete example 3: opposite routing within one powered line

The fixture `exact_source_normalizes_opposite_spelling_inside_one_power_chain`
uses two occurrences of one rational denominator with momenta `+Q(0)` and
`-Q(0)`, plus a balancing denominator of a different mass.

The `+Q` and `-Q` occurrences enter the same owner-local power group because
their canonical signatures and masses agree. Their parsed signatures are both
normalized, and both serial segments are attached to the source-edge-0 scaffold
before receiving one coherent canonical direction. The raw opposite signs
remain in the exact numerator mapper. The separate physical edge completes the
known source scaffold, so the whole exact graph is momentum-balanced without a
signature-derived sign bridge or an owner-sensitive rational residue.

A separate end-to-end test,
`exact_cff_handles_opposite_repeated_routing_without_a_sign_bridge`, uses two
distinct physical source edges whose exact signatures normalize to one
double-pole channel. Their physical provenance remains separate while the
owner-free occurrence graph represents the common rational channel. The test
finds two explicit orientations, keeps the original physical edges undirected,
and verifies the complete orientation sum

\[
  \frac{i}{32\pi^3 E^3}

\]

in GammaLoop's production normalization. No signature-derived incidence bridge
is needed: source routing and numerator sign restoration are handled at their
proper layers.

## Complete example 4: the production-shaped cubic UV rewrite

The most representative scalar fixture is
[`scalar_cubic_exact_uv_rewrite.dot`](../../tests/resources/graphs/uv_tests/scalar_cubic_exact_uv_rewrite.dot).
Its full graph is

```text
incoming -> v5

v0 -> v2   [edge 1, outer LMB carrier]
v3 -> v0   [edge 2]
v0 -> v5   [edge 3]

v2 -> v1   [edge 4, numerator Q(4)^0 + m]
v1 -> v3   [edge 5, inner LMB carrier, numerator Q(5)^0 - 2m]
v2 -> v3   [edge 7]

v1 -> v4 -> outgoing
```

Edges 4, 5, and 7 are the only UV-active cycle. In its leading UV routing they
become

\[
  Q,\quad Q,\quad -Q.
\]

Their owner edges have genuinely different physical parent-graph momentum
signatures. After the local 4D UV expansion, all three UV denominators are the
same massive UV propagator up to routing sign. Each owner nevertheless retains
its edge as physical provenance, and the UV source minor supplies the known
triangle scaffold without reconstruction. The shared CFF layer recognizes the
three normalized signatures as one cubic repeated channel and canonicalizes its
rational occurrence graph without using those owners as algebraic labels.

The integration fixture's numerator probe is deliberately factorized:

\[
  \bigl(Q_4^0+m\bigr)\bigl(Q_5^0-2m\bigr).

\]

The physical analyzer retains these as two separate linear edge dependencies;
it does not first expand their product or replace them by an LMB polynomial.

The analytic companion test then uses the sharper single-channel probe

\[
  (Q_5^0)^2.
\]

All three exact cubic occurrences are certified aliases of physical edge 5's
on-shell energy. The minimax assignment places the two factor slots on two
different occurrences and supplies bounds `(1,1)`. No expanded quadratic
polynomial is constructed.

The unit test `exact_cff_cubic_uv_rewrite_matches_production_convention` uses the
source-compatible routing `(4:+Q,5:+Q,7:-Q)` and permutes the input factor order.
For every order it verifies:

- a valid momentum-balanced exact graph;
- one active topological loop;
- no spurious causal surface from the zero-denominator cograph component;
- canonical ParsedGraph and bare CFF orientation sums independent of input
  factor order;
- the analytic cubic-pole production contour
  \(3i/(128\pi^3E^5)\);
- factorized quadratic bounds `(1,1)`;
- the correct lower-sector parity for a `q_0^2` numerator,
  \(-i/(128\pi^3E^3)\); and
- the temporal/spatial cancellation of a full `q^2` numerator.

The mirrored fixture
[`scalar_mirrored_cubic_exact_uv_rewrite.dot`](../../tests/resources/graphs/uv_tests/scalar_mirrored_cubic_exact_uv_rewrite.dot)
reverses edge 5 and changes the inner carrier so that the leading routing is

\[
  Q,\quad -Q,\quad -Q.

\]

Both fixtures are exercised through three production routes:

1. orientation-local direct local 3D UV;
2. direct local 3D UV with an explicit orientation sum; and
3. local 4D UV followed by exact CFF projection and an explicit orientation
   sum.

The first two routes are parity comparators and do not use the lifting described
in this document: they apply Taylor operators after complete/global CFF energy
integration and differ only by selector omission. Only route 3 reconstructs
exact source occurrences and performs minimax EMR dispatch.

They are compared at several small and large momentum points. The strongest
comparison evaluates the explicit-local-3D and projected-local-4D routes
directly as 1000-bit Arb values and compares the Arb numbers without conversion
through `f64`. This makes the fixture sensitive both to structural errors and
to precision loss in the precision-escalation path. Because the two real graph
numerator factors are odd affine functions of `Q(4)^0` and `Q(5)^0`, while the
original and mirrored UV cycles contain `Q/Q/-Q` and `Q/-Q/-Q` respectively,
their passing three-route comparison is also an end-to-end sign oracle—not only
a denominator-evenness test.

## Complete example 5: UV and cograph topology remain separate

The fixture `exact_cff_separates_uv_topology_from_the_cograph` uses two parallel
bubbles joined by a bridge:

```text
    edges 0,1              edge 4              edges 2,3
 a ========== b -------------------------- c ========== d
```

Edges 0 and 1 are represented with the UV expansion mass; edges 2 and 3 use
their physical masses. Even if canonical momentum signatures coincide after a
contraction, the UV and cograph occurrences live in distinct topology domains.

The exact source therefore preserves both factors rather than accidentally
joining them into one owner- or signature-based power chain. The related
`exact_cff_lu_residue_factorizes_from_quadratic_cubic_spectator` fixture selects
an LU residue on the cograph bubble and combines it with a cubic UV spectator
whose routings are `(+,+,-)` and whose numerator is quadratic. It verifies that
the combined LU residue equals the cograph LU factor times the independent UV
factor, including the production prefactor bridge.

This is important for nested and disconnected UV structures: source-minor
attachment keeps the domains structurally separate. Any later signature
equivalence is repeated-channel algebra on the completed graph and cannot merge
UV with cograph topology merely because their momenta look alike.

## Complete example 6: numerator/denominator cancellation, including LU

The exact-CFF test `exact_cff_uncancelled_powered_denominator_matches_lower_source`
constructs

\[
  \frac{D(Q)\,(Q^0+c)}{D(S)D(Q)^2}.

\]

The positive `D(Q)` remains a provenance-wrapped factor in the numerator. Its
polynomial value is analyzed, while its owner/momentum/mass metadata are not
double-counted as numerator powers. The exact powered source is compared with

\[
  \frac{Q^0+c}{D(S)D(Q)}

\]

and with the ordinary unpowered production CFF. The equality is checked at two
exact rational points. This simultaneously tests:

- a raised denominator;
- a numerator factor which algebraically cancels one copy;
- an odd temporal factor `Q^0+c`;
- exact occurrence mapping; and
- normalization against the ordinary source.

The companion test
`exact_cff_uncancelled_powered_denominator_matches_lower_lu_residues` performs
the same comparison after selecting a physical LU raised residue. It checks the
powered exact source, the lower exact source, and the ordinary CFF after
assembling the complete fixed-cut functional: all generalized residue-map
entries and all raised-order derivative pieces are summed before comparison.
Any stronger per-order agreement observed in this fixture is diagnostic only.

At acceptance level,
[`raised_cut_numerator_cancellation.dot`](../../tests/resources/graphs/raised_cut_numerator_cancellation.dot)
contains two LU cross-section graphs:

- `powered_cancel`, whose edge numerator is
  \(Q_1^\mu Q_{1\mu}-m^2=D(Q_1)\); and
- `lower_bubble`, the corresponding graph with the cancelled propagator
  removed and the compensating graph sign.

The test
`raised_cut_numerator_cancels_one_propagator_in_both_orientation_modes`
evaluates two nontrivial momentum points in all three local-UV routes listed
above. It requires the powered and lower graphs to cancel to relative
`1e-12`, verifies the combined process sum, explicitly sums orientations, and
uses OSE-weighted, explicitly summed LMB channels. Integrated UV terms are
disabled so the oracle isolates the local raised-cut and CFF machinery.

This is the direct LU analogue of the self-energy mechanism in which a
factorized `p_slash` from the integrated UV result combines with a neighboring
propagator numerator to form `p^2` and cancel a denominator.

## Complete example 7: the `epem_a_ddx` GL0 self-energy source

The scalar fixtures isolate the mechanism, but the motivating production case
is the NLO `epem_a_ddx` self-energy graph GL0.

Structured generation traces show two different stages:

| Source | Physical source-edge EMR bounds |
| --- | --- |
| GL0 ordinary full graph | `[(2,1),(3,1),(5,1),(7,1)]` |
| GL2 ordinary full graph | `[(2,1),(3,1),(5,1),(6,1)]` |
| GL0 contracted integrated-UV source | `[(2,1),(3,2),(5,1)]` |
| GL2 first contracted integrated-UV source | `[(2,1),(3,1)]` |
| GL2 second contracted integrated-UV source | `[(5,1),(6,1)]` |

Thus the only genuinely raised physical EMR dependency is the quadratic bound
on edge 3 in the contracted GL0 source. The two propagator edges adjacent to
the self-energy insertion are physically related, but after the vakint result
and its external EMR shift are remapped, edge 3 deterministically owns the
remaining quadratic dependency; edge 5 remains linear. This is a physical-EMR
provenance decision, not an LMB choice.

When a completed local 4D term contains multiple exact occurrences capable of
carrying edge 3's on-shell energy, the minimax plan then performs the second
lift—from the physical degree-two bound to occurrence-local exact bounds. For
two certified equivalent occurrences, it supplies `(1,1)` rather than placing
degree two on one occurrence.

GL2 supplies a useful contrast. Both of its contracted integrated-UV sources
remain linear, as expected for the triangle correction which factorizes the
Born graph. The GL0/GL2 comparison therefore tests that the generalized energy
machinery is activated only where the self-energy numerator actually requires
it.

The finalized `ProcessIntegrand` currently retains the ordinary graph and its
evaluators, not these transient contracted-source bound records. Consequently,
the table above is visible in structured generation traces rather than through
the ordinary post-generation graph-bound API. That inspectability limitation
does not alter the production assignment, but it matters when deciding where a
future acceptance assertion should read the contracted-source metadata.

## Complete example 8: multi-loop and rank-deficient exact sources

The approach is not restricted to one-loop bubbles.

The `sunrise_pow4.dot` shared-core fixture exercises a multi-loop repeated
topology with degree-five numerator support on one occurrence. It verifies that
generalized CFF constructs nonempty known-factor completion sectors and retains
only energy surfaces in denominator trees. Other fixtures exercise a
single-quadratic lower sector and high-contact completions on multi-loop penta
topologies.

The GammaLoop fixture
`exact_rank_deficient_source_keeps_the_complete_active_direction` begins with a
two-loop graph formed by two triangles sharing a vertex. Its exact term retains
only one denominator with momentum equal to the sum of the two parent loop
carriers. The exact source correctly identifies one inactive loop direction and
one complete active direction. It does not invent a physical owner by selecting
an arbitrary LMB coordinate.

The companion initial-cut fixture verifies that a cut carrier remains a literal
external alias even when a stored noncanonical signature contains a component
outside the active exact-source span.

These tests show why LMB data and EMR ownership have distinct responsibilities:
LMB coordinates express and reduce the exact signature space, but physical EMR
provenance determines numerator identity and bounds.

## Generalized-CFF behavior after the lift

Once the exact parsed graph and occurrence bounds have been produced, the
shared `three-dimensional-reps` crate handles generalized residues. Its focused
fixtures cover:

- a repeated box with branching denominator trees;
- one-loop quadratic, cubic, and quartic completions;
- `box_pow3` with repeated quadratic recursive contact and remainder sectors;
- repeated high-power channels with a common auxiliary sampling scale `M`;
- multi-loop sunrise powers and penta contact sectors;
- several choices of which `Q/Q/-Q` occurrence owns an equivalent bound; and
- retention of distinct affine maps even when two terms have the same coarse
  orientation.

In particular, the shared-core fixture
`cff_repeated_quadratic_channel_is_invariant_under_bound_ownership` generates
the same encoded expression for all of

```text
(3,2), (4,2), (5,2),
(3,1)+(4,1), (3,1)+(5,1), (4,1)+(5,1)
```

on one repeated `Q/Q/-Q` channel. Its direct `ParsedGraph` input deliberately
retains signed signatures and verifies that the reversed occurrence receives the
negated sample map. This is a shared repeated-channel and bound-assignment
oracle, not a physical-owner or incidence oracle. GammaLoop's exact source has
already normalized the parsed signatures and restores the raw sign in its
separate numerator mapper.

When a bound exceeds the direct interpolation range, generalized CFF may sample
an EMR energy as a signed integer multiple `a*M` of the common auxiliary scale.
This remains an occurrence-local EMR substitution. It is never an LMB rewrite,
and the final physical expression is required to be independent of the nonzero
value chosen for `M`.

## Why the construction is generic

The confidence claim is deliberately scoped: the construction is generic for
the factorized polynomial EMR numerators and rational denominator powers
produced by GammaLoop's local UV machinery. It is not a claim that every
arbitrary symbolic function can be projected without certification.

Within that production class, genericity follows from the following properties.

### 1. Denominator multiplicity is arbitrary

Any representable negative integer power becomes the corresponding number of
occurrences of its rewritten denominator. The serial-chain constructor accepts
an arbitrary raised-wrapper multiplicity and adds exactly one fewer auxiliary
vertices. Algebraically coincident occurrences are handled by the same channel
machinery regardless of how many distinct owners supplied them. No power-two or
GL0 branch exists in the implementation.

### 2. Skeleton recovery is source-backed and rational identity is owner-free

The source edge and its cograph/UV minor provide the original skeleton,
component domain, and attachment for every process and UV-spinney shape. This
sidesteps the generic inverse problem of deriving graph incidence from momentum
signatures. Within one raised source wrapper, domain, exact momentum signature
up to sign, and mass must agree. Across owners, the same normalized signature
and mass identify the same rational channel. Relabeling compatible owners can
change physical provenance, but not the channel's residue or loop rank.

### 3. Canonicalization is topology-preserving

Cograph and UV graph minors are constructed directly from the physical source,
including contraction of absent edges. Explicit crown hedges retain non-vacuum
external boundaries. Only after this scaffold exists does the owner-blind
occurrence lift and Graphica relabeling establish deterministic rational
equality and caching. No signature matroid, Kirchhoff system, or incidence
search is part of exact-source topology construction.

The only rank solve near this boundary is the two-candidate quotient-space test
for the routing sign of an already attached owner: `exact_row -/+ source_row`
must lie in the span of the opposite source domain. It returns one unique
`+1` or `-1`, or fails. It cannot propose endpoints, connect components, or
synthesize an edge to repair vertex balance.

### 4. Numerator degree composition is structural

The analyzer implements polynomial composition rules over sums, products,
powers, dot products, and multilinear functions. It does not recognize a
particular numerator string or expand it into process-specific monomials.

### 5. The assignment is the exact minimax solution

For disjoint equivalent candidate sets, quotient/remainder balancing is the
mathematical optimum, not a heuristic. Candidate order supplies a deterministic
tie-break without altering the optimum.

### 6. Bound and evaluation cannot diverge

One assignment object drives both generalized-CFF bounds and later numerator
substitutions. Generic future terms cannot accidentally take a different
mapping path after their bounds have been certified.

Batching does not weaken this statement. A shared topology is generated with a
channel-total envelope which dominates every registered term, while each term
still evaluates through its own factor-local plan. Canonical topology identity
includes external edges and affine boundary data, so reuse remains valid for
disconnected and non-vacuum exact sources as well as vacuum UV factors.

### 7. Routing signs are first-class data

Canonicalization up to sign normalizes the denominator signature stored in the
parsed graph. A rational-routing component then receives one coherent canonical
direction, which can reverse that intermediate spelling. The numerator mapper
composes the separate literal occurrence sign with the final parsed-routing
sign and restores the physical numerator convention. This works for odd, even,
scalar, and tensor numerator factors.

### 8. Provenance projection is channel-complete

Owner relabeling does not erase physical information. Occurrence-energy maps
retain the particular source owner, while neutral cut support is the set union
of all distinct owners in the algebraic channel. The union automatically
deduplicates repeated copies of one owner and extends to any channel
multiplicity. Adding raised-group representatives preserves compatibility with
physical cut selection without replacing that owner set. None of these rules
depends on a bubble, triangle, GL0, or a particular number of raised edges.

### 9. Unsupported cases fail closed

The implementation returns contextual errors for:

- a non-integer or unsupported energy-dependent denominator power;
- an opaque nonlinear energy function;
- a physical EMR degree with no certified literal exact occurrence;
- duplicate or overlapping candidate ownership;
- one source edge instantiating incompatible normalized denominators in an
  additive term;
- a source-backed exact topology which violates loop-momentum conservation; and
- a rank which has no valid source-carrier basis.

Failing instead of selecting a convenient owner or LMB coordinate is part of
the genericity guarantee: adding a more complicated nested self-energy cannot
silently fall through to a test-specific approximation.

## Minimality and runtime consequences

The design minimizes work at three levels.

First, the graph lift is minimal: a power `r` requires exactly `r` occurrences,
and raising one source wrapper requires `r-1` auxiliary subdivision vertices.
It does not clone the surrounding graph or algebraically cancel
numerator/denominator pairs upstream.

Second, the numerator lift minimizes the maximum exact energy rank. Generalized
CFF cost grows with the required reconstruction/contact sectors, so replacing a
degree-four bound `(4,0)` by `(2,2)` can materially reduce generation and
evaluation work. The factorized planner obtains that reduction without
expanding the numerator.

Third, retaining a collected common denominator avoids repeating the CFF
recursion for algebraic lower sectors that generalized CFF can generate by
pinching. Canonical batching still joins genuinely outer additive terms without
the rank inflation of a pointwise maximum. The real degree-one triangle
coefficient therefore performs one CFF generation for its uncancelled source;
its base and singly dotted lower sectors are internal CFF sectors rather than
separate upstream graph topologies.

Upstream algebraic cancellation such as `D(Q)/D(Q)^2 -> 1/D(Q)` is deliberately
not performed. Supplying the proper occurrence bounds lets generalized CFF
perform the required pinches and lower-sector reconstruction internally. This
keeps UV orchestration simple and avoids making the production graph topology
depend on process-specific numerator simplification.

## Maintained invariants

Future changes to this path should preserve all of the following:

1. A denominator occurrence is never dropped merely because another occurrence
   has the same owner.
2. The original `EdgeIndex` survives the Taylor operator as `source_edge`.
3. Source owners and the original graph determine the contracted cograph/UV
   scaffold and physical attachments; signatures never reconstruct that
   skeleton or choose original endpoints.
4. The owner-free rational occurrence graph, loop rank, and CFF residue are
   invariant under compatible owner relabeling of denominators with the same
   topology domain, mass, and `D(Q)=D(-Q)` channel.
5. Raising one source wrapper produces the requested occurrence multiplicity
   through minimal serial subdivision; the owner's label is not subsequently
   part of rational CFF identity.
6. UV and cograph denominator node domains never merge.
7. Repeated occurrences of one owner must have the same normalized signature
   and mass, or generation fails.
8. `+Q` and `-Q` share a normalized denominator signature, while the literal
   sign and any later coherent component reversal are composed separately for
   physical-numerator mapping.
9. The quotient-space rank calculation can choose only the unique routing sign
   on fixed source endpoints; it cannot reconstruct or repair incidence.
10. Occurrence-energy projection retains the particular physical owner which
    supplied each occurrence, independently of owner-free rational incidence.
11. Cut support for an exact algebraic channel is the sorted, deduplicated union
    of all distinct physical owners which instantiate that channel.
12. Raised-edge representative aliases are added to unioned cut support for LU
    compatibility; they never replace or erase the original physical owners.
13. A physical numerator `Q` is reconstructed with the physical sign even when
   assigned to a `-Q` occurrence.
14. Numerator energy degrees are computed solely in physical EMR variables.
15. LMB coordinates never own an energy bound or serve as an identity fallback.
16. The numerator remains factorized through analysis and mapping.
17. The same immutable per-term plan owns bounds and substitutions.
18. A shared CFF batch envelope uses the maximum total degree of each repeated
    algebraic energy channel and componentwise maxima only for non-repeated
    occurrences; it never replaces a term's plan.
19. Additive branches reuse capacity; multiplicative and multilinear slots
    consume capacity.
20. Lower-sector pinching does not reassign a factor away from its certified
    occurrence.
21. Input factor order does not alter the canonically relabeled exact source;
    changing only a compatible physical owner changes provenance, not the
    rational graph or residue.
22. Non-vacuum external boundaries remain explicit source-crown hedges.
23. No internal edge or external-balance edge is synthesized from momentum
    signatures; an inconsistent source-backed graph fails validation.
24. Unsupported mappings fail explicitly rather than choosing a convenient
    edge.

## Code map

The principal implementation sites are:

- `crates/gammalooprs/src/uv/approx/local_4d.rs`
  - factorized term projection;
  - negative-power extraction into repeated `FourDDenominator` records;
  - outer Taylor-sum splitting without canceling positive typed `GS.den`
    numerator factors against negative denominator occurrences.
- `crates/gammalooprs/src/graph/three_d_source.rs`
  - exact source coordinates and occurrence ordering;
  - disjoint-set cograph/UV source-minor contraction from original
    `source_edge` attachments;
  - minimal raised-wrapper subdivision and owner-independent rational
    canonicalization;
  - normalized denominator validation and repeated-channel signatures;
  - explicit source-crown boundary completion;
  - unique `+/-` routing validation on already fixed endpoints;
  - post-construction Graphica node/edge canonicalization;
  - literal `+/-Q` candidate certification;
  - occurrence-to-physical energy mapping and routing-sign restoration;
  - channel-wide physical cut-support projection.
- `crates/gammalooprs/src/numerator/energy_degree.rs`
  - factorized physical-EMR degree analysis;
  - equivalent-candidate validation;
  - deterministic minimax assignment;
  - immutable factor-local mapping plan.
- `crates/gammalooprs/src/cff/generation.rs`
  - physical degree extraction;
  - physical-to-exact planning;
  - canonical topology registration;
  - channel-total batch-envelope joining and exact occurrence bounds passed to
    the shared CFF generator.
- `crates/gammalooprs/src/uv/approx/local_3d.rs`
  - in the projected-local4D branch only, two-pass registration then generation
    of completed local-4D Taylor terms;
  - in that branch, retention of each term's original exact numerator plan.
- `crates/gammalooprs/src/cff/mod.rs`
  - exact-source CFF assembly;
  - retention and use of the planned exact-source numerator;
  - physical cut and surface projection;
  - unioned physical-owner support with additive raised-edge representatives.
- `crates/three-dimensional-reps/src/generation.rs`
  - generalized residue, finite-pole, and lower-sector CFF generation from the
    supplied exact occurrence bounds.

## Production acceptance checkpoint

The retained 2026-08-31 10x campaign passes all four physical DD/TT
acceptances (`4/4`) after exercising orientation-local 3D, explicit-sum 3D, and
projected local 4D where the acceptance compares routes. Pulls are signed
differences from the published target in units of the Monte Carlo error; ratio
pulls include the LO uncertainty.

| Acceptance | 10x LO result | 10x NLO result | Graph and ratio evidence |
| --- | --- | --- | --- |
| direct `gamma* -> d d~` | `0.5068703962 +/- 0.0025987972` (`+1.449 sigma`) | `0.01966009810 +/- 0.00053595339` (`+1.424 sigma`) | `GL0=-0.03132123586 +/- 0.00023922726` (`+0.729 sigma`), `GL2=+0.05112479005 +/- 0.00046213299` (`+1.584 sigma`); `alpha_s/pi` pull `+1.141 sigma` |
| converted `e+e- -> gamma* -> d d~` | `0.1950499744 +/- 0.0010326753 pb` (`+1.479 sigma`) | `0.007824189766 +/- 0.000340601513 pb` (`+1.630 sigma`) | signed MC components `GL0=-0.01996339254 +/- 0.00015479207`, `GL2=+0.02786745983 +/- 0.00028425324`; no separate published component targets; `alpha_s/pi` pull `+1.453 sigma` |
| direct `gamma* -> t t~` | `2.901968994 +/- 0.015639978` (`+1.641 sigma`) | `0.2079169992 +/- 0.0042953541` (`+1.489 sigma`) | `GL0=-0.1443600613 +/- 0.0035809931` (`+0.669 sigma`), `GL2=+0.3522770605 +/- 0.0023720361` (`+1.687 sigma`); paper-ratio pull `+1.037 sigma` |
| converted `e+e- -> gamma* -> t t~` | `0.3307052414 +/- 0.0018004843 pb` (`+1.603 sigma`) | `0.02356890542 +/- 0.00056839205 pb` (`+1.058 sigma`) | summed-graph integration has no persisted GL0/GL2 rows; paper-ratio pull `+0.685 sigma` |

The converted DD graph components are closure diagnostics, not separately
published observables, so this design does not assign them artificial targets.
All LO integrations used 100,000 samples; the direct DD NLO used 400,000, the
direct TT graph rows used 400,000 each, and each converted NLO central slot used
200,000.

The scalar LU 15-case matrix covers the same three production routes with local
UV, integrated UV, and threshold counterterms all enabled. As retained
pre-reversal evidence, its 2026-08-31 `dev-optim` / `test_gammaloop` rerun passed
`15/15`, with `235` skipped, in `70.438 s`. The restored post-CFF direct route
requires a current matrix rerun before merge readiness. Four near-zero cases use
the authorized f64-input `1e-14` unit-scale
fallback; the Arb-to-Arb comparison remains run and reports non-scaling. This
is test-oracle handling only and required no production change.

## Validation map

The most relevant focused tests are:

| Layer | Test or fixture | Protected property |
| --- | --- | --- |
| Local-4D term parsing | `term_projection_keeps_factorized_numerator_atoms` | Numerator remains factorized |
| Local-4D term parsing | `term_projection_preserves_dots_and_distinct_same_edge_expressions` | Arbitrary denominator multiplicity survives |
| Local-4D term parsing | `term_projection_keeps_typed_numerator_factors_uncancelled` | Only the outer Taylor sum is split; nested numerator sums and positive typed denominators remain factorized |
| Real Taylor projection | `dod_one_triangle_keeps_one_uncancelled_exact_source_and_matches_lower_sectors` | One common-denominator DOD1 source retains `(2,1,2)` owner multiplicities and factorized typed numerator factors, uses one cache topology, and matches the explicitly reduced lower-sector oracle |
| Exact graph | `exact_source_serializes_dotted_same_edge_occurrences` | Same-owner cubic serial chain |
| Exact graph | `exact_source_owner_relabeling_preserves_residue_rank_and_cut_provenance` | Owner-invariant exact residue and loop rank, distinct energy provenance, and unioned physical cut support |
| Exact graph | `exact_source_keeps_source_instantiated_domains_and_masses_separate` | Source-component and UV/cograph separation |
| Exact graph | `exact_source_normalizes_opposite_spelling_inside_one_power_chain` | Normalized `Q/-Q` signatures on one source-routed chain |
| Exact graph | `exact_uv_component_inherits_source_minor_and_rejects_wrong_provenance` | Multi-loop UV source-minor incidence and validation |
| Exact graph | `exact_source_routes_a_hard_uv_row_modulo_the_soft_cograph_span` | Quotient-space solve chooses only the unique sign on source-fixed endpoints |
| Exact graph | `exact_shifted_factor_reversal_preserves_projected_affine_maps` | Post-construction canonical labels are input-order independent |
| Exact graph | `exact_uv_triangle_reuses_only_three_and_four_edge_cff_topologies` | Canonical topology cache reuse across locally dotted owners |
| Exact graph | `exact_uv_source_retains_a_non_vacuum_two_point_shift` | Pure-external crowns preserve a non-vacuum two-point boundary |
| Numerator plan | `minimax_assignment_balances_two_linear_factors` | `(2 over 2) -> (1,1)` |
| Numerator plan | `minimax_assignment_reuses_capacity_across_additive_branches` | Sum uses maximum degree |
| Numerator plan | `minimax_assignment_repeats_power_bases_without_expansion` | Factorized powers stay unexpanded |
| Numerator plan | `minimax_assignment_uses_canonical_lexicographic_tie_break` | Deterministic `(3 over 2) -> (2,1)` |
| Sign mapping | `exact_energy_bounds_balance_across_routing_signs_of_one_energy_channel` | Negative occurrence restores physical numerator sign |
| Pinching | `exact_energy_mapper_keeps_canonical_zero_sample_when_duplicate_survives` | Assigned zero sample survives a lower sector |
| Exact CFF | `exact_cff_cubic_uv_rewrite_matches_production_convention` | Analytic cubic pole, quadratic rank, and input-order invariance |
| Exact CFF | `exact_cff_uncancelled_powered_denominator_matches_lower_source` | Dotted cancellation with odd numerator |
| Exact LU | `exact_cff_uncancelled_powered_denominator_matches_lower_lu_residues` | Same identity for the complete fixed-cut residue functional |
| Shared CFF | `cff_repeated_quadratic_channel_is_invariant_under_bound_ownership` | Common `Q/Q/-Q` channel independent of bound owner |
| Amplitude acceptance | `scalar_amplitudes_match_across_local_uv_routes` | Local 3D/local 4D equivalence, including native Arb comparison |
| LU acceptance | `raised_cut_numerator_cancels_one_propagator_in_both_orientation_modes` | `q^2-m^2` cancellation in all three routes |

Together these tests cover parsing, source-backed graph construction, canonical
relabeling, sign handling, factorized rank planning, generalized-CFF contact
sectors, exact analytic
residues, local-3D/local-4D equivalence, and LU raised cuts. No single numerical
process test is being used as a substitute for the structural invariants.

## Design conclusion

The essential distinction is:

> Original source edges recover the known UV/cograph skeleton and retain
> physical provenance; normalized denominator algebra defines the
> owner-independent rational occurrence graph, while numerator factors are
> assigned according to certified physical EMR provenance.

Canonical `D(Q)=D(-Q)` normalization is safe because owner labels remain in
energy provenance and cut-support metadata rather than contaminating the
rational residue, while the raw routing sign remains separate data used to
reconstruct the physical numerator. Raised-wrapper subdivision supplies the
correct number of residue variables. The physical cut projection unions all
owners of an algebraic channel and then adds raised-line representatives, so LU
compatibility does not cost provenance. Explicit crown hedges preserve
non-vacuum boundaries, and post-construction Graphica relabeling supplies
deterministic cache identity without any topology inference from Kirchhoff or
incidence matrices. The minimax plan keeps the numerator factorized and supplies
the smallest possible maximal occurrence rank. Because the same plan controls
both CFF generation and evaluation, the construction extends to higher powers,
nested self-energies, disconnected UV components, LU residues, and multi-loop
sources without a process-specific patch.
