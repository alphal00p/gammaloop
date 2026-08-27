# Exact powered-denominator lifting and factorized EMR assignment

## Status and scope

This document describes the production design used when a completed local 4D
term is projected to a causal-flow representation (CFF) after UV Taylor
operations have raised one or more propagators. It expands two tightly coupled
requirements:

1. lift every denominator power onto the incidence of its physical `source_edge`,
   with one occurrence-local edge for every denominator copy; and
2. rewrite the still-factorized numerator over those exact occurrences while
   minimizing the largest energy-degree bound supplied to generalized CFF.

The two requirements cannot be solved independently. A denominator power first
creates the occurrence-local energy namespace in which the generalized residue
problem is posed. The numerator then has to use exactly that namespace, and the
same factor-to-occurrence assignment must control both CFF generation and
numerical numerator evaluation.

This is GammaLoop production machinery. It is not a requirement that the
diagnostic `3Drep` CLI prepare its inputs or expressions in the same way. It is
also deliberately independent of LTD; no LMB variable is ever used as an
energy-identity or energy-capacity index here.

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
2. Every occurrence inherits the incidence of its `source_edge` in either the
   contracted cograph source minor or the separately contracted UV source minor.
   Matching momentum signatures never choose endpoints or merge owners.
3. The `r` occurrences of one source edge are represented by a minimal serial
   subdivision of that edge, with `r-1` auxiliary two-valent vertices. A powered
   self-loop becomes an `r`-edge cycle. Equal denominators on distinct source
   edges retain their distinct incidences.
4. The exact signature and mass of repeated occurrences of one owner are
   canonicalized and validated. This is the safe `D(Q)=D(-Q)` equivalence used
   by repeated-channel algebra, not a graph-reconstruction rule.
5. Physical owner IDs survive separately for cut support, physical-surface
   projection, diagnostics, and mapping back to the parent graph.
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
10. If an occurrence has momentum `-Q` rather than `Q`, that raw algebraic sign
    is retained by the numerator mapper and explicitly inverted when the
    occurrence is used to evaluate a physical numerator `Q`. Parsed-graph
    routing comes from the source edge instead; denominator evenness never
    causes the numerator to be treated as even.
11. Once incidence, power subdivisions, and non-vacuum source-crown hedges are
    complete, Graphica canonically relabels nodes and exact edges. This supplies
    deterministic equality and cache keys without inferring topology.
12. All additive Taylor terms are planned before CFF generation. Terms with the
    same canonical topology share one channel-wise maximal capacity envelope,
    while each term retains its own immutable minimax numerator assignment.

The resulting graph is therefore a source-backed occurrence representation of
the rewritten rational function. Provenance is not merely descriptive metadata:
it selects the source incidence before the owner-free graph is passed to the
shared CFF implementation.

## Terminology and identity layers

Several IDs coexist in this path. Conflating signature equivalence with physical
incidence, or occurrence-local capacity with physical EMR ownership, produces
incorrect topology or numerator sampling.

| Concept | Meaning | May determine rational topology? | May own a numerator bound? |
| --- | --- | --- | --- |
| Physical graph edge | An edge of the original amplitude or forward-scattering graph | Yes, through its source-minor graph attachment | Yes, during the first physical-EMR analysis |
| Provenance or source owner | The original edge recorded in a `den(edge, momentum, mass, value)` wrapper | Yes; it selects the inherited incidence and the local power group | No, by itself |
| Denominator occurrence | One copy produced by a negative denominator power | It subdivides its owner's already selected incidence | Yes, after certified physical-to-exact lifting |
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
only one denominator occurrence. The required operation is local: repeated
wrappers with the same source edge subdivide that edge's incidence. It is also
wrong to discard the source edge and try to reconstruct incidence from the
rewritten signatures. Two different physical edges can carry algebraically
equal denominators while still occupying different places in the source graph.

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
  | source-owner grouping and algebra validation
  v
occurrence groups keyed by source_edge
  validated by (UV/cograph domain, signature up to sign, mass^2)
  |
  | source-minor incidence + owner-local serial subdivision
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

## Stage 1: retain denominator multiplicity without expanding the numerator

The local-4D term projector traverses the Symbolica expression structurally.
For a negative integer power of a denominator wrapper, it repeats the
`FourDDenominator` record `abs(power)` times. Non-denominator powers remain in
the numerator. The first argument of every typed `GS.den` wrapper is the
original `EdgeIndex`; the UV Taylor operator differentiates and rewrites the
other arguments while carrying that owner unchanged. Consequently, the
completed Taylor coefficient still has a direct lookup into the original graph
and never has to recover an edge from the rewritten momentum.

Taylor collection can also expose a positive typed `GS.den` factor multiplying
a negative occurrence of the same wrapper. Projection distributes only an
addition which contains such typed denominator provenance, then cancels the
matching positive and negative occurrences. It does not distribute an ordinary
numerator sum. For example,

```text
((A+B)*den_0 + den_1) / (den_0*den_1)
```

is projected as the two factorized terms `(A+B)/den_1` and `1/den_0`; it is not
expanded into `A/den_1 + B/den_1 + 1/den_0`. The focused regression
`term_projection_cancels_provenance_without_expanding_ordinary_sums` protects
that distinction.

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

## Stage 2: form owner-local power groups and validate their algebra

The structural group key is the physical `source_edge`. For each occurrence,
GammaLoop also computes an exact momentum signature in the temporary source
coordinates and canonicalizes it up to an overall sign. Repeated occurrences of
one source edge are accepted only if the following validation tuple agrees:

```text
(
    uses_uv_loop_basis,
    canonical_momentum_signature_up_to_sign,
    canonical_mass_squared,
)
```

The retained `full_expr` payload and the input order of denominator factors do
not affect this validation. The source owner deliberately does affect the
structural group: it is the only authoritative link back to graph incidence.

This gives two separate invariants:

- repeated copies of one source edge form one powered line only when they
  represent the same rewritten denominator up to `D(Q)=D(-Q)`; and
- equal rewritten denominators on different source edges remain different graph
  edges, even though the shared CFF layer may later classify them as one
  repeated algebraic channel.

If one source edge appears with incompatible normalized signatures or masses in
one additive 4D term, generation fails instead of choosing an incidence by
heuristic. UV-vacuum and cograph denominators are constructed in disjoint source
minor domains even when their algebraic data happen to look alike.

### Provenance is retained, not discarded

Owner provenance remains on every `FourDDenominator`. Before the owner-free
`ParsedGraph` is produced, it is used to:

- select the inherited source-minor incidence;
- decide which repeated occurrences subdivide one physical line;
- project occurrence-local energy IDs back to physical edge IDs;
- retain Cutkosky and LU cut support;
- recognize which exact surfaces are original physical surfaces;
- report and diagnose the source of a UV-expanded factor;
- recover parent-graph spatial momenta and affine external shifts.

The rule is therefore:

> Use source provenance to construct the occurrence topology, and use normalized
> signatures to validate and process the algebra on that fixed topology.

## Stage 3: inherit source-minor incidence and subdivide powers locally

Each owner-local group receives the incidence `(tail, head)` of its physical
source edge.

For cograph denominators, a disjoint-set union over the original graph nodes
contracts source edges absent from the additive exact term. Each surviving
owner then uses the endpoints of its original `EdgeIndex` in that quotient.
UV denominators use a separate disjoint-set union for the UV source minor:
absent UV edges are contracted there, and every active UV owner inherits its
endpoints from that minor. The two node domains remain distinct.

There is no Kirchhoff solver, signature-matroid reconstruction, or
signature-to-endpoint fallback. A rewritten momentum such as `k1+k2` still has
a `source_edge`; its source-minor attachment determines where it lives. Exact
signatures are subsequently checked against that choice. If source-backed
incidence and the supplied momenta violate loop-momentum conservation,
generation fails explicitly.

Once the base incidence is known, a group of multiplicity `r` is subdivided as

```text
tail --occ_0-- x_1 --occ_1-- ... --occ_(r-2)-- x_(r-1) --occ_(r-1)-- head
```

This introduces precisely `r-1` auxiliary vertices. If `tail == head`, the
result is an `r`-edge cycle. Every occurrence keeps its own exact energy ID. Its
`ParsedGraph` signature is normalized up to sign, and every segment is directed
consistently with the physical source-edge routing.

The UV Taylor operator can remove a soft/cograph shift from a denominator. For
example, an original source edge may carry `Q_hard+Q_soft`, while its rewritten
UV denominator carries only `Q_hard`. Its endpoints still come directly from
the original edge. To determine only whether that inherited incidence is used
forward or backward, GammaLoop tests both signs modulo the span of the known
opposite topology domain:

```text
exact_row - sign * source_row in span(opposite_domain_rows), sign in {+1,-1}.
```

Exactly one sign must pass. This rank calculation is a routing-coordinate
validation after both graph minors are already fixed; it never chooses nodes or
repairs momentum balance. The focused hard-plus-soft triangle regression checks
that the correct sign leaves only a soft residual and that the opposite sign
raises the quotient-space rank.

Serial subdivision is the graph-theoretic representation of a dotted line: the
new two-valent vertices do not introduce independent momentum constraints. Raw
occurrence signs are retained separately by the exact numerator mapper rather
than encoded as alternating segment directions.

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

Occurrence input is first put in a stable order so owner-local chains are
reproducible. After all source incidences, auxiliary power vertices, cut
carriers, and crown hedges exist, Graphica canonically relabels the directed,
edge-coloured graph. Exact edges are then sorted in that canonical node
namespace and the occurrence-to-original-factor map is reordered with them.
Source-node and power-node names survive as aliases of the canonical nodes;
cache keys deliberately clear those names and edge labels. Reversing the order
of Symbolica factors therefore does not change the parsed exact graph or its CFF
cache entry. Crucially, this pass only relabels an already constructed graph: it
never derives an endpoint from a momentum signature.

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
`dod_one_triangle_taylor_terms_keep_source_owners_and_reuse_two_cffs` applies
the real degree-one UV Taylor operator to the three-edge UV triangle in the
double-triangle skeleton. Its four additive terms retain owner multiplicities

```text
(1,1,1), (1,1,1), (1,1,2), (2,1,1).
```

The two undotted terms share the base three-denominator UV vacuum component;
the two owner-local derivatives share its four-denominator dotted subdivision.
Although their sparse numerator plans place capacity on different canonical
occurrences, the batch generates exactly those two CFF topologies. The fixed
cograph component is part of the canonical key, so this reuse does not conflate
different surrounding graphs or future non-vacuum boundaries.

## Why `D(Q)=D(-Q)` cannot introduce numerator sign mistakes

This deserves a separate proof because the denominator and numerator have
different parity properties.

Let an exact occurrence carry momentum

\[
  P_r=s_r Q,
  \qquad s_r\in\{+1,-1\}.
\]

For a scalar propagator,

\[
  D(P_r)=P_r^2-m^2=(s_r Q)^2-m^2=D(Q).
\]

It is therefore correct to give `+Q` and `-Q` the same normalized denominator
signature. If they are repeated occurrences of one source edge, they validate
as copies of that owner's powered line. If they belong to distinct source edges,
they keep those distinct incidences while remaining algebraically equivalent to
the repeated-channel machinery. Neither case derives topology from the sign.

It does **not** replace `Q` by `-Q` in the numerator. The relative sign `s_r`
is stored next to the occurrence when literal exact candidates are certified.
The occurrence-local CFF energy map represents the temporal component of
`P_r`,

\[
  p_r^0=s_r q^0.
\]

When that occurrence is assigned to a physical numerator factor written in
terms of `Q`, the mapper reconstructs

\[
  q^0=s_r p_r^0,
\]

because `s_r^{-1}=s_r`. In the implementation this is the explicit negation
performed when the certified occurrence sign is `-1`.

Consequently,

\[
  Q^0+c
  \longmapsto
  \begin{cases}
    p_r^0+c, & P_r=+Q,\\
    -p_r^0+c, & P_r=-Q.
  \end{cases}
\]

An odd numerator therefore retains its sign. A vector numerator is handled in
the same way: spatial components remain in the physical parent-graph EMR
basis, while the corrected temporal energy is inserted with the temporal unit
vector. The mapper never performs a blanket substitution `Q -> -Q` on the
physical numerator.

There are thus two separate routing layers:

1. The parsed graph uses one direction derived from the physical source-edge
   routing. All segments of an owner-local power chain use that direction after
   their signatures have been normalized.
2. The numerator mapper retains each occurrence's raw sign and composes it with
   the occurrence-local temporal map to recover the physical `Q` convention.

The first layer is graph provenance; the second is algebraic numerator mapping.
Denominator evenness affects only signature equivalence and never silently turns
an odd numerator into an even one.

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

## Complete example 2: coincident denominators retain distinct-owner incidence

The fixture
`exact_source_keeps_coincident_distinct_owner_incidences`
uses a physical graph with two parallel edges:

```text
        edge 0
    a ----------> b
        edge 1
    a ----------> b
```

It creates one denominator record from each owner. In the production basis their
signatures normalize to the same momentum and their masses agree, so the shared
repeated-channel detector correctly recognizes one algebraic double channel.

The source-backed exact graph nevertheless contains the two original parallel
incidences and no auxiliary power node: neither owner is individually raised.
The repeated-group detector returns `[0,1]` only after those incidences are
fixed. The physical projection separately records which canonical exact
occurrence came from owner 0 and which came from owner 1.

This fixture separates two concepts which must not be conflated: source
provenance determines graph incidence, while normalized signatures determine
denominator equivalence for generalized repeated-channel algebra.

## Complete example 3: opposite routing within one powered line

The fixture `exact_source_normalizes_opposite_spelling_inside_one_power_chain`
uses two occurrences of one rational denominator with momenta `+Q(0)` and
`-Q(0)`, plus a balancing denominator of a different mass.

The `+Q` and `-Q` occurrences enter the same owner-local power group because
their canonical signatures and masses agree. Their parsed signatures are both
normalized, and both serial segments follow the routing of source edge 0. The
raw opposite signs remain in the exact numerator mapper. The separate physical
edge supplies the reverse source incidence needed to close the graph, so the
whole exact graph is momentum-balanced without a signature-derived sign bridge.

A separate end-to-end test,
`exact_cff_handles_opposite_repeated_routing_without_a_sign_bridge`, uses two
distinct physical source edges whose exact signatures normalize to one
double-pole channel. Their source incidences remain separate. The test finds two
explicit orientations, keeps the original physical edges undirected, and
verifies the complete orientation sum

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
its edge in the UV source-minor triangle. The shared CFF layer recognizes the
three normalized signatures as one cubic repeated channel on that fixed
triangle; it is not reconstructed as a source-independent serial line.

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
incidence keeps the domains structurally separate. Any later signature
equivalence is repeated-channel algebra on the completed graph and cannot merge
their nodes or power chains.

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
powered exact source, the lower exact source, and the ordinary CFF residue
per-residue rather than only after a global sum.

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
occurrences of that source edge. The serial-chain constructor accepts an
arbitrary owner-local group length and adds exactly one fewer auxiliary
vertices. No power-two or GL0 branch exists in the implementation.

### 2. Incidence is source-backed and algebra is validated separately

The source edge and its cograph/UV minor determine endpoints for every process
and UV-spinney shape. Within one owner-local power group, domain, exact momentum
signature up to sign, and mass must agree. Across owners, the same normalized
signature and mass may identify a repeated algebraic channel without changing
either incidence.

### 3. Canonicalization is topology-preserving

Cograph and UV graph minors are constructed directly from the physical source,
including contraction of absent edges. Explicit crown hedges retain non-vacuum
external boundaries. Only after this graph exists does Graphica relabel it for
deterministic equality and caching. No signature matroid, Kirchhoff system, or
incidence search is part of exact-source topology construction.

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
parsed graph. Chain incidence follows the source-edge routing. The separate raw
occurrence sign remains in the numerator mapper and restores the physical
numerator convention. This works for odd, even, scalar, and tensor numerator
factors.

### 8. Unsupported cases fail closed

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

First, the graph lift is minimal: a power `r` of one source edge requires exactly
`r` occurrences and `r-1` auxiliary vertices on that incidence. It does not
clone the surrounding graph or algebraically cancel numerator/denominator pairs
upstream.

Second, the numerator lift minimizes the maximum exact energy rank. Generalized
CFF cost grows with the required reconstruction/contact sectors, so replacing a
degree-four bound `(4,0)` by `(2,2)` can materially reduce generation and
evaluation work. The factorized planner obtains that reduction without
expanding the numerator.

Third, canonical topology batching avoids repeating the CFF recursion for
Taylor terms which differ only in which equivalent owner was dotted or where a
factor-local bound was assigned. Channel-total joining avoids the rank inflation
of a pointwise maximum. The real degree-one triangle coefficient therefore
performs two CFF generations for four additive terms: one base three-denominator
source and one four-denominator dotted source.

Upstream algebraic cancellation such as `D(Q)/D(Q)^2 -> 1/D(Q)` can still be an
optional optimization, but correctness does not depend on it. Supplying the
proper occurrence bounds lets generalized CFF perform the required pinches and
lower-sector reconstruction internally. This keeps UV orchestration simple and
avoids making the production graph topology depend on process-specific
numerator simplification.

## Maintained invariants

Future changes to this path should preserve all of the following:

1. A denominator occurrence is never dropped merely because another occurrence
   has the same owner.
2. The original `EdgeIndex` survives the Taylor operator as `source_edge`.
3. Every occurrence inherits the cograph- or UV-minor incidence of that
   `source_edge`; signatures never choose endpoints.
4. Only repeated occurrences of one owner subdivide one incidence into a power
   chain. Distinct owners never form a cross-owner chain.
5. UV and cograph denominator node domains never merge.
6. Repeated occurrences of one owner must have the same normalized signature
   and mass, or generation fails.
7. `+Q` and `-Q` share a normalized denominator signature, while the raw sign is
   retained separately for physical-numerator mapping.
8. The quotient-space rank calculation can choose only the unique routing sign
   on fixed source endpoints; it cannot reconstruct or repair incidence.
9. A physical numerator `Q` is reconstructed with the physical sign even when
   assigned to a `-Q` occurrence.
10. Numerator energy degrees are computed solely in physical EMR variables.
11. LMB coordinates never own an energy bound or serve as an identity fallback.
12. The numerator remains factorized through analysis and mapping.
13. The same immutable per-term plan owns bounds and substitutions.
14. A shared CFF batch envelope uses the maximum total degree of each repeated
    algebraic energy channel and componentwise maxima only for non-repeated
    occurrences; it never replaces a term's plan.
15. Additive branches reuse capacity; multiplicative and multilinear slots
    consume capacity.
16. Lower-sector pinching does not reassign a factor away from its certified
    occurrence.
17. Input factor order does not alter the canonically relabeled exact source.
    Reassigning a factor to a different physical owner may intentionally change
    topology because the owner supplies incidence.
18. Non-vacuum external boundaries remain explicit source-crown hedges.
19. No internal edge or external-balance edge is synthesized from momentum
    signatures; an inconsistent source-backed graph fails validation.
20. Unsupported mappings fail explicitly rather than choosing a convenient
    edge.

## Code map

The principal implementation sites are:

- `crates/gammalooprs/src/uv/approx/local_4d.rs`
  - factorized term projection;
  - negative-power extraction into repeated `FourDDenominator` records;
  - typed positive/negative `GS.den` cancellation without expanding ordinary
    numerator sums.
- `crates/gammalooprs/src/graph/three_d_source.rs`
  - exact source coordinates and occurrence ordering;
  - disjoint-set cograph/UV source-minor contraction from original
    `source_edge` incidence;
  - owner-local serial power-chain construction;
  - normalized denominator validation and repeated-channel signatures;
  - explicit source-crown boundary completion;
  - unique `+/-` routing validation on already fixed endpoints;
  - post-construction Graphica node/edge canonicalization;
  - literal `+/-Q` candidate certification;
  - occurrence-to-physical energy mapping and routing-sign restoration.
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
  - two-pass registration then generation of completed local-4D Taylor terms;
  - retention of each term's original exact numerator plan.
- `crates/gammalooprs/src/cff/mod.rs`
  - exact-source CFF assembly;
  - retention and use of the planned exact-source numerator;
  - physical cut and surface projection.
- `crates/three-dimensional-reps/src/generation.rs`
  - generalized residue, finite-pole, and lower-sector CFF generation from the
    supplied exact occurrence bounds.

## Validation map

The most relevant focused tests are:

| Layer | Test or fixture | Protected property |
| --- | --- | --- |
| Local-4D term parsing | `term_projection_keeps_factorized_numerator_atoms` | Numerator remains factorized |
| Local-4D term parsing | `term_projection_preserves_dots_and_distinct_same_edge_expressions` | Arbitrary denominator multiplicity survives |
| Local-4D term parsing | `term_projection_cancels_provenance_without_expanding_ordinary_sums` | Typed denominator cancellation does not expand an ordinary numerator sum |
| Real Taylor projection | `dod_one_triangle_taylor_terms_keep_source_owners_and_reuse_two_cffs` | Four DOD1 terms retain source owners but generate only base and dotted CFF topologies |
| Exact graph | `exact_source_serializes_dotted_same_edge_occurrences` | Same-owner cubic serial chain |
| Exact graph | `exact_source_keeps_coincident_distinct_owner_incidences` | Distinct source incidence plus repeated-signature equivalence |
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
| Exact LU | `exact_cff_uncancelled_powered_denominator_matches_lower_lu_residues` | Same identity per LU residue |
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

> Denominator powers subdivide the incidence supplied by their physical source
> edge; normalized denominator algebra validates that lift, while numerator
> factors are assigned according to certified physical EMR provenance.

Canonical `D(Q)=D(-Q)` normalization is safe because it cannot change
source-backed incidence, while the raw sign remains separate data used to
reconstruct the physical numerator. Owner-local serial subdivision supplies the
correct number of residue variables. Explicit crown hedges preserve non-vacuum
boundaries, and post-construction Graphica relabeling supplies deterministic
cache identity without topology inference. The minimax plan keeps the numerator
factorized and supplies the smallest possible maximal occurrence rank. Because
the same plan controls both CFF generation and evaluation, the construction
extends to higher powers, nested self-energies, disconnected UV components, LU
residues, and multi-loop sources without a process-specific patch.
