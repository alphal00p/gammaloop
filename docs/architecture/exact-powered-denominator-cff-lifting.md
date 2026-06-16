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
2. rewrite the still-factorized numerator over those exact occurrences, then
   select the fewest native generated source-map rows among at most three
   certified proposals, with rank supplying proposal order and count tie-breaks.

This exact-source reconstruction and bounded-dispatch machinery is exclusive to
the projected local-4D route. Both direct local-3D representations instead
complete the loop-energy integration first and apply every UV Taylor operator
to the complete/global CFF expression. Their bodies are identical; the
explicit-sum form only omits the localized form's orientation selectors.

The two requirements cannot be solved independently. A denominator power first
creates the occurrence-local energy namespace in which the generalized residue
problem is posed. The numerator then has to use exactly that namespace, and the
same factor-to-occurrence assignment must control both CFF generation and
numerical numerator evaluation.

This is GammaLoop production machinery, independently of LTD; no LMB variable
is ever used as an energy-identity or energy-capacity index here.

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
its source-local residues is stored once, together with its frozen localizing
factor. Different Taylor topologies may have child maps of different lengths;
those maps are consumed inside their own sums and are never paired with outer
maps. Each actual outer CFF branch maps the completed coefficient times the
remaining numerator once, then multiplies the frozen factor and contributes to
the checked cut-indexed sum. No artificial child production hosts or selector
copies are constructed. This distinction preserves both the complete fixed-cut
functional and the factorized numerator representation.

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
8. Original factors retain their owner occurrence. Rank analysis proposes a
   baseline assignment and at most two deterministic alternatives for admissible
   new factors. The proposal producing fewer actual native CFF map rows wins;
   the rank/canonical proposal order breaks count ties.
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
12. Each additive Taylor term prepares its proposals once. Counts and generated
    payloads use the same canonical topology/options/per-occurrence-capacity key.
    Only winning payloads enter the full cache; a losing contender leaves a
    count-only memo. Each term retains its selected immutable numerator plan.

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
                       | rank-ordered proposals, K <= 3
                       v
            exact occurrence -> energy bound
            factor -> exact occurrence assignment
                       |
                       | generate/lookup actual native map counts
                       v
             chosen payload + matching immutable plan
                       |
             +---------+---------+
             |                   |
             v                   v
      selected CFF          exact numerator mapper
      payload               used by every orientation,
                            residue, and contact sector
```

Each candidate plan is immutable before its generalized CFF request. Selection
retains the winning generated payload with that same plan and exact-source
numerator mapper. This is the certification boundary: the expression eventually
sampled is the expression whose bounds were supplied to CFF.

Projection constructs each canonical exact source and its bounded proposal set
once. It reuses an expression only when canonical topology, occurrence-local
capacity and generation options match, and maps each term with its selected
plan. Equal physical energies do not allow the cache to redistribute their
individual bounds. Independent requests never combine into a larger Cartesian
capacity, so no preliminary registration pass is needed.

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
\(C_a=0\), in a centered hard chart. A compatible chart built from retained
physical edge momenta may instead represent the same vacuum integral by a
balanced affine loop translation. Their signed incidence sum must vanish at
every UV vertex; they do not represent momentum entering the vacuum graph.
The exact source certificate retains them on both sides rather than rejecting each
shifted denominator separately. A future on-shell counterterm may retain a
fixed insertion such as \((m,0,0,0)\); it is simply included on both sides of
the same check.

For nested Taylor operators, the enclosing UV node fixes the hard/soft split.
If a compatible carrier \(Q\) has soft part \(S\) in that node's reference
LMB, its rescaling in the projected route's inverse-hard scale is
\(Q\mapsto (Q-S)/t+S\). The same affine transport
applies to a retained child's hard-momentum payload; its owner and provenance
role remain fixed. Literal external-coordinate carriers stay soft even when
their out-of-domain LMB rows vanish. At the vacuum-integration boundary these
external coordinates must also stay fixed when solving for canonical loop
momenta. Finite MUV integration retains one Laurent term beyond the maximal
pole order; the legacy forest owner requests `wood.max_loops + 1`, matching
its finite-term extraction contract.

Vacuum masses within the component under Taylor expansion retain the hard
rescaling `mUV -> mUV/t`. Direct local-3D Taylor expansion applies the equivalent
forward scaling `mUV -> lambda*mUV` to an integrated coefficient only when the
current subgraph contains its owner. Disjoint owners remain fixed; a later
enclosing operation can scale them. The direct route tags each connected
coefficient's mass as transient `mUV(owner)` before multiplying coefficients,
using the existing symbol and subgraph encoding. Output copies restore physical
`mUV`, while stored sectors retain ownership for later Taylor operations. The
separate normalized localization kernel stays frozen throughout.

Four-dimensional and `IntegratedCts` runtime behavior is unchanged.
`IntegratedCts` projects the coefficient as `C(mUV) * s^(4*r)`, where
`s = uvIntegratedLoopScale` and `r` counts consumed loops. The four-dimensional
`s -> s*t` transformation restores only the consumed-measure power; it does not
cancel the coefficient's mass rescaling in an enclosing operation. Integration
and physical projection set `s = 1`. The active denominator rewrite introduces
the enclosing vacuum mass independently. See
[Signs and Integrated Projections](uv-renormalization.md#signs-and-integrated-projections)
for the integrated-owner and frozen-kernel distinction.

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
masses, and fixed external data. Tests may erase provenance tags for this
comparison, but must preserve the graph numerator's factors, including in
diagnostic copies. Use exact factor cancellation and normalize scalar CFF
denominators before inserting the numerator. The separate analytic UV-subgraph
expansion boundary is described in [UV renormalization](uv-renormalization.md).

This check has four important consequences:

1. The original numerator is never redistributed by the dispatch search. Its retained owner and
   hard-momentum payload determine an exact UV-EMR representative directly.
2. Only new hard numerator energy factors produced by differentiating a
   denominator may be assigned among the serial copies of that same line.
3. Such an assignment is correct before it is optimal. A deterministic first
   valid copy is sufficient; cost-driven redistribution is allowed only after
   the boxed identity passes and only among those degenerate copies.
4. Canonicalizing \(D(Q)\) as \(D(-Q)\) cannot discard an odd-numerator sign.
   The complete signed coordinate map from the tagged hard momentum to the
   selected UV occurrence participates in the boxed equality. Denominator
   evenness changes the denominator channel, never the numerator transformation.

New soft Taylor factors are outside this child-contour assignment. Their
separate, exactly certified cograph routing occurs at the outer CFF boundary
described in Stage 6; it never moves an original or child hard factor.

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
boundary in the enclosing forest. The isolated `1zs/T2` comparison below
locates the remaining discrepancy. Changing this fixed numerator assignment
would violate the exact certificate and cannot be a valid repair.

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

The frozen reconstruction fixture represents the complete truncation
`T0+T1+T2`, rather than the isolated `T2` coefficient. Its exact denominator
occurrence multiset is

```text
owner 5: +q, +q       (production occurrence ids 9 and 12)
owner 6: -q, -q, -q   (production occurrence ids 10, 11, and 13)
```

and the still-factorized numerator is algebraically

\[
 N_{\leq 2}=-g^2(q^0)^2
 \left[
   -A B^2 +(U+p^2)A B-4x^2A+U B^2+2xA B
 \right].
\]

This is the authoritative full-truncation left-hand side. It is obtained
before any CFF call, and no numerator/denominator cancellation is performed
or needed.  The
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
 =N_{\leq 2}.
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
 \widehat N-N_{\leq 2}=0
\]

exactly.  This is the real post-plan certificate: it includes every selected
occurrence and every signed hard/raw/parsed conversion, rather than comparing
an earlier pre-plan expression.  The parallel denominator certificate is also
unchanged and exact.  Values of the formal (P_j) away from the common-LMB
diagonal are not an alternative physical routing and are not used to judge
the reconstruction identity.

The executable regression
`gl00_gl04_planned_lifts_match_post_t_numerators_in_common_loop_coordinates`
runs both reconstruction proofs on the production source. It builds the actual
owner-5/owner-6 denominator records, asks the production `EnergyPowerAnalyzer`
for an immutable plan, and applies it with `ExactSourceEnergyMapper`. Distinct
occurrence witnesses first check that original factors stay on their base
owner and derivative-created factors stay within that owner's copies; the
test does not prescribe the particular allocation shown in this recorded
trace. It then imposes the common chart and proves the complete factorized
numerator identity. The denominator certificate compares the retained source
multiset `D(+q)^2 D(-q)^3` with the reconstructed occurrence product, including
mass and UV-domain ownership. The GL04 case certifies the full `T0+T1+T2`
truncation; isolated coefficient extraction and the complete three-route
comparison require the additional checks below.

### Isolated coefficient and complete-source certificate

Writing $D=q^2-M^2$ and stripping the common numerator coupling, the actual
second Taylor coefficient of the child is

\[
 I_2(q,p)=(q^0)^2\left[
   -\frac{p^2+2M^2}{D^3}
   +\frac{4(q\cdot p)^2}{D^4}
 \right].
\]

This equality is a diagnostic expression for the coefficient, not a production
cancellation or expansion. The production source retains its factorized
numerator and owner-backed occurrences. The recorded 2026-09-05 diagnostic
extraction separately certified the child and the complete child-plus-outer source after applying
all immutable occurrence assignments and the signed $H=hR$, $P=rR$ maps.
Both numerator differences vanished exactly in independent neutral
four-vectors, and both denominator multisets matched. The complete source has
physical occurrence counts

```text
{1:1, 2:1, 3:1, 4:1, 5:2, 6:3, 7:1, 8:1}.
```

The proof must resolve every EMR component, retain independent crown
momenta, and never install a replacement for `emr_vec_index(edge,cind(0))`:
that spatial temporal component is already literal zero. Replacing it would
silently turn every zero in the diagnostic into an energy variable.

An independent child contour uses $E=\sqrt{|\boldsymbol q|^2+M^2}$ and
$a=\boldsymbol q\cdot\boldsymbol p$. Its positive-pole residue sums are

\[
 R_0=\frac{1}{4E},\qquad
 R_1=\frac{a}{8E^3},\qquad
 R_2=\frac{2M^2-|\boldsymbol p|^2-(p^0)^2}{16E^3}
       +\frac{a^2}{8E^5}.
\]

Each real-axis energy integral contributes $-iR_k$ with the usual Below
closure. Expanding the exact energy-integrated two-propagator child gives the
same coefficients. For the complete cut `[1,2]`, the additional moving child
pole in the outer contour begins at order $t^8$ before the $t^{-3}$ spatial
measure, so it cannot change these three coefficients.

At the rational sample

\[
 P_0=(1,\boldsymbol 0),\quad
 \boldsymbol K_0=(3,4,0)/10,\quad
 \boldsymbol K_1=(1,2,3)/10,\quad
 \boldsymbol K_2=(-2,1,1)/10,
\]

with $M=0$, physical masses one on edges 0 and 8 and zero on every other
edge, and all graph numerator factors one except $(Q_5^0)^2$, the complete isolated `T2` contour evaluates in GammaLoop's
normalization to

\[
 -1.8354223786936939047533107951846543304391732719501496
 \times10^{-6}\,i.
\]

Before the scalar-base correction, the complete one-call exact CFF agreed
with projected local-4D at approximately 300 digits, giving
$-2.79328040709378476003397681050696066768826159163\times10^{-6}i$.
The independent analytic contour agreed with direct local-3D. Thus the
reconstruction and factorized child/outer composition passed before the
shared generalized scalar contact was changed. In that recorded post-correction diagnostic, direct,
projected, and one-call results agree with the analytic value: one-call
minus direct is numerically zero at 1024-bit precision, and projected minus
direct is about $-1.25\times10^{-306}i$. The unchanged production GL04
`q5_temporal_square` acceptance also passed with the complete Taylor
truncation, integrated UV counterterms, and thresholds enabled.

### GL00 odd-numerator child certificate

The same regression includes the complete post-Taylor child source for the
physical GL00 edge-5 temporal probe used by
`scalar_3l_cross_section_gl00_q5_temporal_component_inspects_match`. Its
one-loop UV child has physical owners 4 and 5, with $Q_4=q$, $Q_3=p$ and
$Q_5=-q-p$. The original $Q_5^0$ stays on owner 5; $p$ is soft boundary data.
Both child vertices contribute the retained factor $g^2$, where
`g = UFO::SCALAR_COUPLING`.

Write $D_v(q)=q^2-M_v^2$ with the complete Minkowski square. For a soft
Taylor parameter $\lambda$, the source numerator is $-q^0-\lambda p^0$,
while its denominators are $D_v+O(\lambda^2)$ and
$D_v+2\lambda(p\cdot q)+O(\lambda^2)$. The independent complete `T0+T1`
oracle is therefore

\[
 g^2\,\frac{(-q^0-p^0)D_v(q)+2q^0(p\cdot q)}{D_v(q)^3},
 \qquad p\cdot q=p^0q^0-\sum_{i=1}^{3}p^iq^i.
\]

Before normalization, production stores distinct mass fields: the full
propagator expression uses $M_v^2$ (`mUV`), while exact-occurrence metadata
uses $M_e^2$ (`mUVexp`). The regression preserves this distinction. Its
numerator proof compares the complete factorized tagged payload
$g^2[d_5(H_{5,\mathrm{fixed}}^0-p^0)+
2(p\cdot H_{5,\mathrm{derived}})H_{5,\mathrm{fixed}}^0]$, with $H_5=-q$,
to the numerator above. Separately, its denominator proof checks the full
product $D_e(+q)D_e(-q)^2$, $D_e(q)=q^2-M_e^2$, including physical owner
multiplicities `{4:1, 5:2}`, masses and the child UV domain. The original
factor stays on the base occurrence; the denominator derivative uses only
owner 5's new occurrence. The odd terms retain their sign independently of
the denominator's evenness.

`local_4d.rs` introduces the two mass fields, and final integrand
simplification identifies `mUVexp = mUV` after exact mapping. Applying that
specified substitution yields the normalized rational oracle above. The
GL04 fixture also preserves these fields separately. Neither numerator is
distributed into a polynomial for the comparison.

This is a complete child-sector reconstruction certificate. It does not
certify the full cut sum or every GL00 forest. The GL00 numerical regression
complements it with complete physical-route evaluations including local and
integrated UV and threshold counterterms.

### One source convention across recursive scalar bases

Reconstruction correctness ends before CFF normalization. Let $G_c$ denote a
raw independently generated component, including positive local $1/(2E_j)$
factors. Its target source-frame coefficient is

\[
 \beta_c=(-1)^{N_c}
 \begin{cases}
 B_{\mathrm{core},c},&\text{ordinary component},\\
 B_{\mathrm{den},c},&\text{generalized component}.
 \end{cases}
\]

For generalized and ordinary public components, the signed contour is
$J_c=\beta_c G_c$, using $dq^0/(2\pi i)$ for each loop. Internal pure-CFF
catalogues retain their native orientation normalization, distinguished below.
Here $N_c$ counts denominator occurrences and $B_{\mathrm{den},c}$ is the
numerator-bound-independent scalar denominator frame. The generalized core
sign is already encoded in its Laurent functional; multiplying it again
would let an unused degree allowance change a scalar integral. In the
ordinary GammaLoop path, surface conversion removes the positive half-edge
factors and restores $1/\prod_j(-2E_j)$, which already supplies $(-1)^{N_c}$.
A generalized expression retains its positive factors on each variant, so its
typed adapter supplies that parity. More precisely, the actual GammaLoop
adapter is `product_c [(-1)^(rho*N_c) * frame_c]`: `rho` is one when the whole
expression retains variant-local positive factors, while each component chooses
its own core or denominator frame. Thus mixed expressions also apply this parity
to their ordinary components. Fresh exact sources and stored production roots
obey the same rule; they do not require separate fresh-versus-stored formulas.

The same convention must survive denominator deletion *inside* a generalized
functional. `BoundedCffBuilder` and `KnownFactorCffBuilder` carry the original
source prefactor through every remainder, contact, and recursive reentry.
`LowerSectorCffBuilder` converts a scalar base by
$\beta_{\mathrm{native}}/\beta_{\mathrm{source}}$ once, at its embedding
boundary. It derives the native factor from the base that was actually built:

- A connected nonterminal pure component with $n$ denominator occurrences,
  rank $r$, and repeated-channel excess $d$ has native factor
  $(-1)^{n+r-1+d}$.
- A free one-line pure component contains both scalar orientations and evaluates
  to $1/E$. Its signed Below contour is $-1/(2E)$, so its native factor is
  $-1/2$. Both generalized scalar bases and public ordinary sources apply
  this conversion once. The public standalone source retains both maps with
  half weight; the embedded source retains one complete Below residue.
- A selected residue basis already includes its contour closure. Its additional
  stored inherited coordinate-Jacobian sign must be removed on conversion.
  For $w=p-q$, reversal of the real integration limits cancels $dq=-dw$;
  the residue generator already accounts for the reversed closure.

Component products use their actual native factors and their assembly sign,
not metadata inferred from an unrelated fresh child source. This requires a
rational prefactor, not merely a sign. Energy convergence of the starting
source is part of this contour contract: a surviving free one-denominator
component can act only on a scalar quotient in that energy variable.
Public ordinary sources also require their denominator rows to span every
loop energy; initial-state cut aliases cannot supply missing denominator
rank. Formal lower-sector contacts keep their internal rank-reduction path.

The GL04 outer graph exposes the missing conversion particularly clearly.
Edges 3 and 8 carry the same momentum with different masses, so they remain
distinct denominator channels. The complete physical cut obeys both

\[
 (p^0)^2=D_3+E_3^2,\qquad
 \frac{(p^0)^2}{D_3D_8}
 =\frac{E_3^2}{E_3^2-E_8^2}\frac1{D_3}
  -\frac{E_8^2}{E_3^2-E_8^2}\frac1{D_8}.
\]

Deleting edge 3 changes the denominator-count parity while leaving the
nonterminal core metadata unchanged. Before the fix, ordinary scalar and
both lower-topology contours passed independently, but the generalized
quadratic differed by minus twice the physical edge-3 contact. The coefficient
of $(p^0)^2$ in $R_2$ is $-1/(16E^3)$, giving exactly the full live `T2`
discrepancy without fitting a sign to the UV result.

The permanent regressions
`gl04_outer_unequal_mass_quadratic_obeys_complete_cut_identity`,
`mixed_carrier_contact_matches_independent_convergent_contour`, and
`attached_tadpole_repeated_bubble_matches_independent_convergent_contour`
compare against literal pole residues in both generation contexts. Their
scalar, lowered, and terminal controls distinguish a source-frame error from
an incorrect contour oracle. Complete physical-cut sums are the comparison
boundary; branch counts and individual generalized keys are not route
invariants. The recorded initial-cut repeated-spectator contour check also
passed unchanged, checking the inherited-closure sign needed by the original
root contact.

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
`term_projection_preserves_complete_factorized_values` protects that
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

The regression `term_projection_preserves_complete_factorized_values` rebuilds
the complete expression from every projected term. Its cases include wrappers
with the same physical edge and momentum but distinct `full_expr` provenance,
repeated powers, positive numerator denominators, and distinct outer topologies.
It compares exact factorized values instead of prescribing storage counts. A
factorized expression such as

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

> Use the original source graph and retained owners to recover the known
> skeleton and all physical maps. Add only derivative-created serial copies
> of those same lines. Normalized denominator algebra classifies channels on
> that certified scaffold; it never supplies incidence or reassigns original
> numerator owners.

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

The source constructor certifies each owner's same-momentum, same-mass serial
occurrences. The mapper exposes the retained base as the fixed-factor candidate
and derivative-created copies as the new-hard-factor candidates. This uses
source incidence and provenance, not a search over equal on-shell expressions
or a largest equality class. Candidate order is canonical and deterministic.

Candidate sets for different physical edges must be disjoint. A missing,
empty, duplicated, or overlapping candidate set is an error. In particular,
GammaLoop does not treat an arbitrary shifted momentum `Q_e+p` as an alias for
`Q_e`, and it does not use an LMB coordinate as a fallback.

The literal `+/-Q_e` restriction is important for sign safety. It supplies both
an algebraic energy identity and an exact routing sign. Equality of positive
on-shell energies alone would not be sufficient to map an odd numerator.

## Stage 6: propose certified assignments and select by actual map count

Original factors have a single retained owner occurrence. Newly introduced hard
factors may choose only among the certified serial copies of their own line.
The factorized expression supplies a Pareto frontier of admissible loads: sums
take componentwise maxima, while products and multilinear slots add loads.
The rank baseline compares the complete descending envelope lexicographically,
then uses deterministic assignment tie-breaks. It prefers `(4,2,2)` to `(4,3,1)`
even though their maximal degree is equal. An original quartic factor stays
fixed. This rank ordering proposes candidates; actual map count selects the
production assignment.

### Rank-ordered hard proposals

For the special case of `d` freely assignable unit factors and `n` equivalent
occurrences with no fixed load, any valid assignment has nonnegative loads

\[
  \ell_0,\ldots,\ell_{n-1},
  \qquad
  \sum_i \ell_i=d.
\]

The smallest possible maximum load is

\[
  \min \max_i \ell_i = \left\lceil\frac{d}{n}\right\rceil.
\]

For this rank objective, the exact frontier optimum is the
quotient/remainder distribution:

\[
  q=\left\lfloor\frac{d}{n}\right\rfloor,
  \qquad r=d\bmod n,
\]

and assigns load `q+1` to the first `r` canonical candidates and `q` to the
rest. This is whole-envelope optimal and deterministic; it does not predict
native map count. With fixed loads or shared factors inside sums, the frontier
retains the actual factor structure instead of applying this special-case
formula to a summed degree.

Rank-baseline examples are:

| Physical degree | Equivalent occurrences | Exact loads |
| ---: | ---: | --- |
| 2 | 2 | `(1,1)` |
| 3 | 2 | `(2,1)` |
| 4 | 2 | `(2,2)` |
| 5 | 3 | `(2,2,1)` |
| 7 | 3 | `(3,2,2)` |

The planner assigns individual admissible factor slots, not merely total
degrees. For three newly denominator-derived unit factors,

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

The original expanded polynomial is never constructed. These examples do not
permit redistributing a pre-existing numerator factor across other owners.

### Native map-row selection

The hard search proposes the rank baseline and at most two alternatives, each
changing a single owner's assignment. Each owner's first two distinct
alternative occurrence bounds come from the existing frontier; the complete
descending envelope and deterministic assignment order select up to two global
challengers. Equal sorted envelopes with different occurrence positions remain
distinct proposals.

For each proposal, the ordinary exact source generator receives the same source
and options with that plan's occurrence bounds. The score is the actual stored
row count `generated.expression.orientations.len()`. This is the generator's
operational output after its normal complete-map emission/coalescing; retained
zero rows count too. It is neither the sum of key lengths nor a degree-product
formula, and it has no extra physical-orientation multiplier. Different Taylor
terms keep their own source maps even when key lengths differ.

The smaller count wins; an equal count retains the earlier rank-ordered proposal.
The same generated payload, assignment and bound report move forward together.
Existing rank/Pareto pruning only defines the proposal class: no monotonicity
of map count under rank reduction or independence between owners is assumed.
K<=3 bounds the number of candidates, not each generation's time or memory, and
does not promise a global minimum. A generation error remains an error rather
than silently discarding an unsupported contender. Stage 8 describes how the
hard cache avoids rebuilding or retaining losing trial payloads.

### Soft Taylor energies at the outer CFF boundary

Denominator differentiation also creates soft factors `S^0`. Their
`DenominatorDerivedSoft` provenance records the actual crown carrier and literal
payload. A later enclosing Taylor operator consumes a newly hard part on that
carrier, while any remaining soft part keeps its provenance. This is separate
from assigning new hard factors among serial copies of a UV denominator.

At final projected assembly, the completed child coefficient is multiplied by
the untouched outer numerator before choosing soft routings. For each tagged
soft carrier, the active cograph momentum rows provide exact rational basis
alternatives. Every candidate must reproduce all loop coefficients, with the
fixed external shift explicitly restored. No on-shell energy equality or
change of original numerator ownership supplies this certificate. Concrete
spatial components retain their original routing and contribute zero energy
rank.

The same factorized sum/product rules combine candidate loads into a Pareto
frontier. Up to three proposals with distinct ordered active-capacity vectors,
ranked by descending-envelope/deterministic order, enter the raw outer-source
generator. Selection uses the same native row count as the hard path, before
surface conversion, cut selection or hosting;
proposal order breaks ties. Final selected `(host, source map)` branches can
have a different count and are not the objective of this raw-row search.

Only the incumbent and one challenger are retained at a time. The loser is
dropped before any persistent graph surface, bound-report or cut/hosting effect.
The winning raw payload is converted once and reused with its matching
factorized Atom for numerator mapping. This soft boundary has no persistent
count memo. When an empty contraction reuses an existing production root
expression, the rank baseline follows that existing path without regenerating
or changing the root's established capacity.

Child contour maps are already consumed and play no role in this routing. A
loop-dependent carrier outside the active span is an error, not a constant in
the rank analysis. The search is best-of-three within the existing certified
basis/Pareto proposal class, not a global count optimum. Basis enumeration and
frontier construction can still grow combinatorially; the three-candidate
generation budget does not bound that preparation work. Neither preparation
nor scoring expands the numerator or searches for cancellation by expansion. The direct-3D
and LTD paths retain their existing behavior.

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

## Stage 8: cache exact topology and per-request capacity

Each completed Taylor coefficient builds its source-backed exact graph and
bounded immutable proposal set once. Counts and generated expressions use the
same complete canonical `ParsedGraph`, energy-edge map and generation-options
key, including the unchanged per-occurrence bounds. Every term is evaluated
with its selected plan.

`ExactCffGenerationCache` keeps a small count-only memo for each scored key and
admits a full generated payload only when that key wins a request. Known
contenders are compared by count without cloning their payloads. A key which
previously lost can later win a different comparison; if its payload was not
retained, it is generated only after selection. Previously winning cache entries
remain reusable. Cold comparisons retain only an incumbent and one challenger,
dropping a losing trial immediately. A request makes at most three core calls,
including any recovery of a count-only winner, and never regenerates a freshly
generated winner solely to convert or map it.

No maximum is taken across independent requests: such a join can manufacture
combinations of high ranks on different channels which no individual numerator
needs. Nor does the cache sum or redistribute bounds between occurrences of
equal physical energies. Only requests with identical occurrence capacities
share a generated expression. This keeps cache reuse independent of the
factor-to-occurrence assignment and avoids repeating source construction and
rank planning in a separate registration pass.

Provisioning one outer CFF for independently evaluated residue branches is a
different boundary. That single source requires the per-edge maximum over the
factorized branch products, including their common remaining numerator. The
existing physical-EMR analyzer consumes these products one at a time and merges
only degree maps, with initial-cut and tree-edge energies excluded. Equal bodies
can share an analysis, but opposite bodies cannot cancel across branches. No
tagged symbolic sum or expanded numerator is constructed, and the original
branch expressions, residue keys and numerator maps remain unchanged.

The production-shaped regression
`dod_one_triangle_keeps_separate_denominator_topologies`
applies the real degree-one UV Taylor operator to the three-edge UV triangle in
the double-triangle skeleton. The coefficient retains its natural denominator
topologies with owner multiplicities

```text
(1,1,1), (2,1,1), (1,1,2).
```

The base and singly dotted terms are separate exact sources. No common
denominator is imposed across them, so denominator clearing does not manufacture
positive typed factors in their numerators. The parser preserves each original
factorized numerator and can group terms only when their typed denominators
already agree. Genuine positive typed denominator factors remain supported;
their CFF pinches are distinct from artificial factors introduced by regrouping.
The fixed cograph component is part of each canonical key, so cache reuse does
not conflate different surrounding graphs or future non-vacuum boundaries.

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

The test
`exact_energy_bounds_keep_original_on_canonical_occurrence_across_routing_signs`
uses three equivalent occurrences with sign patterns `(+,-,-)` and `(+,+,-)`.
An original degree-two numerator stays on its retained base occurrence and
maps to `E^2`; changing the other routing signs cannot authorize dispatch to a
serial copy.

The complementary test
`denominator_derived_dispatch_composes_each_occurrence_routing_sign`
dispatches two derivative-created temporal factors over derived copies. It
varies the immutable hard sign `h`, the literal occurrence signs, and the parsed
routing sign. After the signed conversion, the complete mapped numerator is

\[
  (hE+a)(hE+b).
\]

The independent shifts retain odd powers of `E`, so a lost numerator-side sign
cannot be hidden by an even final power. The test also forces an abstract vector
factor onto each eligible derived occurrence and checks the same physical
spatial and temporal hard-momentum convention.

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

The unit fixture `exact_source_preserves_the_complete_cubic_propagator_contour` starts
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

The source represents the repeated propagator by serial occurrences on its
contracted physical incidence. The regression checks valid momentum balance,
retained physical ownership, and the complete signed contour against the
independent cubic-pole residue $-3/(16E^5)$. That oracle detects lost or duplicated
powers without fixing auxiliary-node counts or canonical edge ordering.

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
`exact_cff_keeps_opposite_source_routing_without_a_sign_bridge`, uses two
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
triangle scaffold without reconstruction. The shared CFF layer keeps all three
occurrences and their numerator energies independent, even though their physical
denominators coincide. Canonical graph labels do not merge their rank bounds or
numerator samples.

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

All three exact cubic occurrences share the relevant algebraic energy, but
this quadratic probe is an original numerator factor. It retains its canonical
source occurrence and supplies one degree-two bound, as the companion test
explicitly requires. Only newly denominator-derived factors may enter the
bounded dispatch search. No expanded quadratic polynomial is constructed.

The unit test `exact_cff_cubic_uv_rewrite_matches_production_convention` uses the
source-compatible routing `(4:+Q,5:+Q,7:-Q)` and permutes the input factor order.
For every order it verifies:

- a valid momentum-balanced exact graph;
- one active topological loop;
- no spurious causal surface from the zero-denominator cograph component;
- complete signed CFF values independent of input factor order;
- the analytic cubic-pole production contour
  \(3i/(128\pi^3E^5)\);
- complete contours for original and derivative-created quadratic factors,
  retaining their distinct ownership domains;
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
exact source occurrences and performs bounded EMR dispatch.

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

The physical degree-two bound is not itself permission to redistribute its
factors. A completed local-4D term retains original factors on their source
occurrences. If new denominator-derived factors have two eligible serial
copies, the rank baseline may propose `(1,1)`; the native map-row comparison
then selects among the bounded proposals. The physical-source trace above
neither records that selected occurrence plan nor proves its map count.

GL2 supplies a useful contrast. Both of its contracted integrated-UV sources
remain linear, as expected for the triangle correction which factorizes the
Born graph. The GL0/GL2 comparison therefore tests that the generalized energy
machinery is activated only where the self-energy numerator actually requires
it.

`State::generation_summary(...)` exposes contracted-source diagnostics in each
graph report's `stats.cff_energy_degree_bound_reports` field.
The DDx acceptance reads those reports and distinguishes the physical parent
bounds from exact-source assigned bounds. They are generation diagnostics,
separate from the ordinary graph retained by `ProcessIntegrand`; acceptance
checks should use this summary rather than infer transient bounds from the
final evaluator.

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
`cff_repeated_occurrence_bounds_preserve_equivalent_physical_numerators`
compares the complete signed scalar and quadratic functions requested by

```text
(3,2), (4,2), (5,2),
(3,1)+(4,1), (3,1)+(5,1), (4,1)+(5,1)
```

on a repeated `Q/Q/-Q` denominator at two deterministic points. The fixture
reverses both endpoints and the full momentum signature, verifies the graph,
and includes the reversed numerator sign. Each representation acts on the
factor-local polynomial class specified by its bounds; equivalent physical
numerators must yield the same complete function. Separate analytic quadratic
and quartic contours retain independent mathematical oracles.
GammaLoop's exact source normalizes denominator signatures and restores the raw
sign in its separate numerator mapper.

Higher-degree reconstruction interpolates one occurrence at a time and divides
each interpolation basis polynomial only by that occurrence's denominator.
Contracting that edge passes an exact, factorized quotient to the existing
lower-sector CFF completion. It does not replace several equal denominators by a
single powered channel or impose a common sampling node on their numerators.

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
vertices. Algebraically coincident occurrences retain separate bounds and
samples regardless of how many distinct owners supplied them. No power-two or
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

The rank solve at source-topology validation is the two-candidate quotient-space test
for the routing sign of an already attached owner: `exact_row -/+ source_row`
must lie in the span of the opposite source domain. It returns one unique
`+1` or `-1`, or fails. It cannot propose endpoints, connect components, or
synthesize an edge to repair vertex balance.

### 4. Numerator degree composition is structural

The analyzer implements polynomial composition rules over sums, products,
powers, dot products, and multilinear functions. It does not recognize a
particular numerator string or expand it into process-specific monomials.

### 5. Actual native map rows select among certified rank proposals

Rank/Pareto planning retains fixed owner loads and the factorized sum/product
rules. Comparing descending envelopes proposes balanced lower degrees even when
another occurrence fixes the maximum; quotient/remainder balancing is its
special case for freely assignable unit factors. The production choice uses
actual native source-map rows among K<=3 certified proposals, with rank and
deterministic order breaking count ties. This restricted proposal class does
not imply a globally minimal map count.

### 6. Bound and evaluation cannot diverge

One assignment object drives both generalized-CFF bounds and later numerator
substitutions. Generic future terms cannot accidentally take a different
mapping path after their bounds have been certified.

Batching does not weaken this statement. Reuse requires equal canonical
topology and identical per-occurrence capacity, while each term still evaluates
through its own factor-local plan. Canonical topology identity includes
external edges and affine boundary data, so reuse remains valid for disconnected
and non-vacuum exact sources as well as vacuum UV factors.

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

Second, the numerator lift compares actual native map rows for a small set of
rank-ordered proposals. Assigning four eligible new hard unit factors as
`(2,2)` instead of `(4,0)` can reduce reconstruction/contact work, but topology
and complete-map coalescing determine the actual count. The bounded search may
prefer a less balanced proposal when it produces fewer rows. Original factors
remain fixed, and soft routing uses its separate exact off-shell certificate.
The additional contender has a generation cost of its own; a net runtime or
memory benefit requires measurements rather than an envelope-only argument.
All proposal and selection work preserves the factorized numerator.

Third, retaining natural Taylor denominator topologies avoids artificially
raised propagator powers and compensating numerator factors. Canonical batching
reuses expressions only for equal topology and identical per-occurrence capacity,
avoiding rank inflation across independent terms. The real degree-one triangle
therefore presents its base and singly dotted terms as separate source
topologies, each with its own factorized numerator and capacity.

For genuine numerator factors, upstream algebraic cancellation such as
`D(Q)/D(Q)^2 -> 1/D(Q)` is deliberately not performed. Supplying the proper
occurrence bounds lets generalized CFF perform the required pinches and
lower-sector reconstruction internally. Preserving natural Taylor topologies
does not require such cancellation: it avoids creating the compensating factors
in the first place, without making source incidence depend on process-specific
numerator simplification.

## Maintained invariants

Future changes to this path should preserve all of the following:

1. A denominator occurrence is never dropped merely because another occurrence
   has the same owner.
2. The original `EdgeIndex` survives the Taylor operator as `source_edge`.
3. Source owners and the original graph determine the contracted cograph/UV
   scaffold and physical attachments; signatures never reconstruct that
   skeleton or choose original endpoints.
4. Compatible owner-relabeling comparisons require valid source-backed lifts
   on both sides, with the same topology domain, mass and `D(Q)=D(-Q)` channel.
   Their equal rational residue is not permission to infer incidence from
   denominators or to move an original numerator factor off its retained owner.
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
18. CFF cache reuse requires equal topology and identical per-occurrence bounds.
    Bounds of degenerate occurrences are never summed or redistributed;
    independent terms retain their own plans.
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
25. Public initial-state cut slots carry the canonical `Default` direction;
    their signed fixed energies remain in the exact energy maps. This label
    normalization occurs after fusion and preserves unfiltered residue IDs and
    their order; it cannot supply missing denominator rank or re-sign a
    numerator factor.

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
  - deterministic rank-ordered hard and certified soft proposals (K<=3);
  - immutable factor-local mapping plan.
- `crates/gammalooprs/src/cff/generation.rs`
  - physical degree extraction;
  - physical-to-exact planning and raw native map-row selection;
  - count-only contender memoization and winning-payload reuse, with keys
    retaining exact topology, options and occurrence bounds;
  - raw ordinary-source generation before persistent surface conversion.
- `crates/gammalooprs/src/uv/approx/projected_4d.rs`
  - single-pass proposal preparation and cached generation of completed
    local-4D Taylor terms;
  - retention of each term's selected exact numerator plan and matching payload.
- `crates/gammalooprs/src/uv/approx/local_3d/residue_localizer.rs`
  - projected soft-proposal selection before surface conversion;
  - shared cut/hosting conversion for the winning payload and ordinary paths.
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
pulls include the LO uncertainty. These are historical measurements, not a
fresh validation of the current tree. The converted `epem_a_ttx` acceptance
compares complex magnitudes, so its pass alone does not certify the signed
result or overall phase. The direct `gamma* -> t t~` test separately retains
fixed signed LO and graph-NLO component checks.

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
requires fresh validation after subsequent source changes. Four cases in that
historical run used an f64-input `1e-14` unit-scale fallback. Current nonzero
comparisons require finite values and precision-scaled relative agreement;
only an exact source-zero certificate permits a separate absolute bound.

## Validation map

The core source certificates and graph fixtures enter with local UV
reconstruction. The acceptance tests
`scalar_amplitudes_match_across_local_uv_routes` and
`raised_cut_numerator_cancels_one_propagator_in_both_orientation_modes` are
introduced after that core boundary with the command execution workflows.

The most relevant focused tests are:

| Layer | Test or fixture | Protected property |
| --- | --- | --- |
| Local-4D term parsing | `term_projection_preserves_complete_factorized_values` | Complete reconstruction with repeated powers, distinct wrappers, outer Taylor sums and positive numerator denominators; no storage-count requirement |
| Real Taylor projection | `dod_one_triangle_keeps_separate_denominator_topologies` | DOD1 retains natural `(1,1,1)`, `(2,1,1)` and `(1,1,2)` denominator topologies with factorized numerators and no manufactured clearing factors |
| Exact graph | `exact_source_preserves_the_complete_cubic_propagator_contour` | Complete signed cubic-pole contour and physical ownership |
| Exact graph | `exact_source_owner_relabeling_preserves_residue_rank_and_cut_provenance` | Owner-invariant exact residue and loop rank, distinct energy provenance, and unioned physical cut support |
| Exact graph | `exact_source_keeps_source_instantiated_domains_and_masses_separate` | Source-component and UV/cograph separation |
| Exact graph | `exact_source_normalizes_opposite_spelling_inside_one_power_chain` | Normalized `Q/-Q` signatures on one source-routed chain |
| Exact graph | `exact_uv_component_inherits_source_minor_and_rejects_wrong_provenance` | Multi-loop UV source-minor incidence and validation |
| Exact graph | `exact_source_routes_a_hard_uv_row_modulo_the_soft_cograph_span` | Quotient-space solve chooses only the unique sign on source-fixed endpoints |
| Exact graph | `exact_shifted_factor_reversal_preserves_projected_affine_maps` | Reversing shifted-denominator order preserves the complete projected residue |
| Exact graph | `exact_uv_triangle_cached_and_uncached_residues_agree` | Cached and uncached complete residues agree for the base and every locally dotted owner |
| Exact graph | `exact_uv_source_retains_a_non_vacuum_two_point_shift` | Pure-external crowns preserve a non-vacuum two-point boundary |
| Numerator plan | `assignment_plans_preserve_owner_domains_and_complete_factorized_values` | Every proposal retains fixed/derived domains, sufficient energy capacity and the complete factorized numerator across sums, powers and multilinear functions |
| Sign mapping | `exact_energy_bounds_keep_original_on_canonical_occurrence_across_routing_signs` | Original factors retain their base owner across routing signs |
| Sign mapping | `denominator_derived_dispatch_composes_each_occurrence_routing_sign` | Dispatched derivative-created factors preserve signed physical hard momentum |
| Fixed-owner sampling | `exact_energy_mapper_uses_selected_fixed_owner_sample_with_uv_mass` | The selected occurrence sample is consumed before its physical UV-mass on-shell replacement |
| Dispatch proposals | `assignment_plans_preserve_owner_domains_and_complete_factorized_values` | Every bounded proposal retains fixed owners and exact factors |
| Dispatch and cache parity | `bounded_exact_cff_dispatch_preserves_owner_contour` | Complete selected-plan and cached/uncached functions against an independent contour |
| Exact CFF | `exact_cff_cubic_uv_rewrite_matches_production_convention` | Analytic cubic pole, quadratic rank, and input-order invariance |
| Exact CFF | `exact_cff_uncancelled_powered_denominator_matches_lower_source` | Dotted cancellation with odd numerator |
| Exact LU | `exact_cff_uncancelled_powered_denominator_matches_lower_lu_residues` | Same identity for the complete fixed-cut residue functional |
| Shared CFF | `cff_repeated_occurrence_bounds_preserve_equivalent_physical_numerators` | Complete equivalent physical numerators across every `Q/Q/-Q` bound distribution |
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
incidence matrices. The bounded search keeps the numerator factorized and
chooses fewer native map rows among its certified proposals. Because the same
selected plan controls both CFF generation and evaluation, the construction
extends to higher powers,
nested self-energies, disconnected UV components, LU residues, and multi-loop
sources without a process-specific patch.
