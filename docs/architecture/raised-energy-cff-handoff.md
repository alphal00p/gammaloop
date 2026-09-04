# Raised-energy CFF rebase handoff

All discrepancy investigations in this work must follow the mandatory
[bird's-eye discrepancy triage](../../CONTRIBUTING.md#discrepancy-triage) before
entering CFF, UV, or evaluator internals. In particular, first identify the last
shared pipeline stage, compare standalone or forest artifacts at the first
differing boundary, and reduce totals by graph, forest, cut, residue order, LMB,
orientation, and term. A selector-local versus explicit-sum mismatch is a
selector partition-of-unity problem until a direct branch comparison proves
otherwise.

For direct local-3D UV counterterms, `explicit_orientation_sum_only` is an
orientation-representation choice, not a UV-construction strategy. Both values
first complete the loop-energy integration and then apply every UV Taylor
operator to the complete/global CFF expression. They build the same resulting
local Taylor/CFF bodies: ordinary mode multiplies them by orientation selectors,
where each selector identifies one complete generalized residue-map key, while
explicit-sum mode retains each keyed body once without those selectors. A
physical edge-orientation signature is only metadata inside that key and must
never be used to merge or re-host distinct generalized maps.
Expanded local-4D projection, including exact-source reconstruction and minimax
EMR dispatch, remains a separate construction selected only by
`local_uv_cts_from_expanded_4d_integrands`. Any future direct
physical-shell algorithm must have its own setting and must support both
orientation representations identically.

## Authoritative route-comparison boundary

**Equality between the direct local-3D and projected local-4D routes is an
equality of the complete residue functional for one fixed Cutkosky cut.** It is
checked only after summing every generalized residue-map entry and assembling
every derivative piece belonging to all raised-order slots of that cut. A
single generalized key, contact/remainder family, or `lu_cut_order` entry is an
internal decomposition term, not a separately invariant physical quantity.
The direct and projected routes may legitimately use different numbers of such
terms and redistribute finite contributions among them.

For a double-pole residue, write the assembled functional schematically as

\[
  \mathcal R_C=\left.\partial_\eta f_2(\eta)\right|_{\eta=0}+f_1(0).
\]

The redefinition

\[
  f_2(\eta)\longmapsto f_2(\eta)+\eta h(\eta),
  \qquad
  f_1(\eta)\longmapsto f_1(\eta)-h(0)
\]

leaves \(\mathcal R_C\) unchanged while changing both stored raised-order
pieces. Higher powers admit the analogous triangular redistributions between
adjacent derivative slots. Therefore a disagreement in one `lu_cut_order`
followed by compensation in another is not a route mismatch. Conversely,
agreement of one order cannot certify the full cut. Keep the physical
Cutkosky-cut identity fixed, assemble its required jets and raised orders, sum
all residue maps, and compare only that result. Per-key UV finiteness remains a
separate required property of the orientation-local direct-3D construction;
it does not create a key-by-key comparison with projected local-4D.

Projected component composition must preserve numerator factorization while
forming that sum. Each independently integrated component carries a pair
`(C_i, N_i)`: its CFF carrier `C_i` and its still-factorized numerator `N_i`.
The component's exact source map acts only on the factors it owns. The mapped
pair is then passed to the next component, and the outer production map acts
once on the still-unmapped soft/cograph factors. Products remain products;
the implementation must not expand a full numerator to align residue keys.
Cartesian multiplication is required only for genuinely independent,
unsummed component states. Repeated production hosts of an already summed
child coefficient are mapping metadata and must not multiply that coefficient
again.

The projected local-4D route has one non-negotiable reconstruction invariant.
For every completed Taylor term, the numerator retained in the original hard
sub-LMB and the candidate numerator written in the EMRs of the reconstructed
UV graph **must become exactly the same symbolic function after both are
expressed in one common set of formal loop momenta**. The original edge owner
and hard-LMB payload carried through the Taylor operator provide the left-hand
side; the owner-built UV skeleton, including only derivative-created serial
copies, provides the right-hand side. This exact identity is the correctness
gate before generalized CFF generation. Any EMR spelling which passes it is
valid; minimizing its energy-rank envelope is a later performance optimization,
not part of the definition. The comparison may expand expressions in tests or
debug assertions, but production numerator construction and evaluation must
remain factorized. LMB coordinates are used here only to prove equality of two
representations; they never become EMR identities or CFF capacity labels. See
the dedicated “common-LMB commutative reconstruction invariant” section in
[`exact-powered-denominator-cff-lifting.md`](exact-powered-denominator-cff-lifting.md).

This is a **stop rule**, not an optional diagnostic.  For every projected-4D
discrepancy, do not change CFF signs, contact terms, prefactors, or test targets
until the *actual production candidate* has been checked in this order:

1. take the completed post-`T` numerator in its stored compatible hard sub-LMB;
2. build the UV skeleton solely from retained source-edge owners, adding only
   derivative-created serial copies of those same lines;
3. apply the immutable production energy-assignment plan, including its exact
   selected occurrences and the complete `H = h R`, `P = r R`, hence
   `H^0 = h r P^0`, sign conversion;
4. substitute both the post-`T` numerator and that final mapped candidate into
   one neutral set of formal `K`/`P` four-momenta and prove their exact
   difference is zero;
5. independently certify denominator momentum (modulo only the even identity
   `D(Q) = D(-Q)`), mass, multiplicity, and component domain;
6. only after both certificates pass, compare generalized-CFF output and final
   residue aggregation with residue-summed direct local-3D.

A comparison made before the assignment plan or before occurrence-routing
signs is not this certificate.  Nor is equality of only the numerical final
answer.  Algebraic expansion is permitted for a disposable test/debug copy at
step 4; the production numerator must remain factorized throughout.

Last updated: 2026-09-03

This is the machine-independent continuation note for the raised-energy CFF
work. The detailed chronological record, evidence, and decision history remain
in [`kysvnqlq-rebase-review.md`](kysvnqlq-rebase-review.md). Read that ledger
before changing architecture or tests; this file is the shorter operational
view.

The ledger's **Current readiness** table is authoritative. Its later sections
are chronological evidence and deliberately retain some superseded “pending”
or old-version wording; do not reopen resolved rebase or top-level-state work
from a historical checkpoint alone.

## Checkout and canonical revisions

The published remote bookmark `raised_energy_cff_wip` now points at commit
`5c3e7f`. This is the current recoverable remote checkpoint. A working tree may
still contain newer validation or test-contract cleanup, so compare it with
that checkpoint before assuming the local tree is fully published. Use the
bookmark rather than reviving one of the old conflicted source revisions:

```bash
jj git fetch --remote origin
jj bookmark track raised_energy_cff_wip@origin
jj new raised_energy_cff_wip
```

If the bookmark is already tracked but points elsewhere, first inspect the
divergence and then set it explicitly to `raised_energy_cff_wip@origin`. Do not
merge a hidden, stale remote tip by assumption.

The original published starting point for this continuation was:

- change `lktlsypu`, commit `2bbc6f59460c`, bookmark
  `raised_energy_cff_wip`;
- rebased base `main`, change `llvxkzxn`, commit `395610143576`;
- deferred proper-LTD source reserve `proper_ltd_core`, change `rpspmzvw`,
  commit `377dd390366f`.

`proper_ltd_core` is a pushed sibling, not an ancestor or descendant of
`raised_energy_cff_wip`. The old `ltd_in_gammaloop_symbolica_update` and
`kysvnqlq` revisions are historical source material, not continuation points.
The unrelated conflicted bookmark `push-wwvzyzvysvro` is outside this ancestry
and must not be resolved as part of this work.

The non-LTD continuation extended that original starting point with exact EMR
capacity mapping (`ozntmkwv`), typed stored production roots (`toxllnzs`),
integration-workspace versioning (`onotmrxw`), sampling-scale validation
(`qyprnryk`), interface maintenance (`qykwxoxl`), typed amplitude archive reload
coverage (`trrvmwrw`), stored-root selector preservation (`zroqlyxz`), mapped
numerator scalarization regression coverage (`uukzotwz`), the initial production
sign bridge (`yyksullr`), exact source-local EMR carriers through integrated
finite-UV localization (`opmnpruw`), and the full typed CFF global-prefactor
bridge (`kypmvxsq`). Those early focused checks have since been superseded by
direct-`gamma*` DDx and TTX acceptances, three-route scalar coverage, and
end-to-end UV-profile tests described below. The current local workspace also
contains validation newer than the `5c3e7f` remote checkpoint; publish it only
after the remaining test-contract and final-suite gates are resolved.

Before editing, read `AGENTS.md`, `CONTRIBUTING.md`, this handoff, and the full
decision ledger. Use Jujutsu and start a new described change for each approved
unit of work.

## Status at a glance

| Area | Status | Continuation meaning |
| --- | --- | --- |
| Rebase onto current `main` | **Complete** | Fourteen retained changes are replayed directly onto `39561014`; the selected ancestry is conflict-free. The obsolete `kysvnqlq` pin was abandoned because it only downgraded main's Symbolica stack. |
| Latest-main evaluator/archive overlap | **Resolved and validated** | Main's Symbolica, native-code, O2, evaluator construction, and command ordering are retained. Deferred explicit CFF orientation sums are layered onto that architecture. |
| Direct-3D raised-energy CFF | **Post-CFF Taylor route restored; matrix rerun pending** | Both direct orientation representations now apply every Taylor operator after complete/global CFF energy integration and differ only by selector materialization. Physical EMR/source-edge bounds and the repeated-channel CFF construction are retained. The protected GL11 `Q1^2`/`Q7^2` snapshots remain unchanged; A73 records why their denominator-only cancellation residues are not generalized-CFF value oracles. |
| 4D-local UV then CFF with raised energy | **Source-backed dotted lifting and minimax EMR dispatch implemented** | A79 carries original edge provenance through the Taylor operator, constructs the UV skeleton from those owners, serially subdivides raised powers, and minimizes factorized numerator-energy rank before canonical CFF batching. Sign-normalized denominator signatures validate channels but never reconstruct incidence or rewrite the sign-sensitive numerator. The empty forest remains the ordinary factorized production root in both routes. |
| Integration-workspace persistence | **Implemented and focused-validation complete** | Integration workspaces have an independent version-1 manifest checked before both resume decoders. Incompatible workspaces return contextual `integrate --restart` guidance instead of decoding positional bincode or panicking (A60). |
| Other GammaLoop completeness | **Physical acceptances retained; restored-route scalar matrix pending** | The retained 10x campaign passes direct and converted DD/TT normalization (`4/4`) against the published absolute or inclusive-ratio targets. DD also checks per-graph `M_UV` and localization-scale independence plus opposite GL0/GL2 `mu_r` slopes; its converted graph rows remain targetless closure components. The pre-reversal scalar LU 15-case all-enabled rerun passed `15/15`, with `235` skipped, in `70.438 s` under `dev-optim` / `test_gammaloop`; that is historical evidence, not a current gate result. The restored post-CFF direct route must complete the matrix rerun before merge readiness. The default divergent-only profiler covers amplitude and LU domains. Remaining gates are that rerun, stale shared-CFF test-contract cleanup, final curated tests, and latest-main synchronization. |
| Computed UV-forest export | **Deferred; not exact yet** | Persisted graph terms retain only coarse orientations. Production generation/evaluation already uses exact maps, but export must regenerate the production catalogue or persist exact map IDs before it can claim exact affine-map identity. |
| Proper LTD | **Explicitly deferred** | Keep the `Cff`/`Ltd` selection; `Ltd` must immediately return the typed unsupported error. Port proper LTD later from the sibling reserve onto the finalized CFF plus 4D architecture. |
| Auxiliary `3Drep` CLI | **Diagnostic only; not a production contract** | Its input normalization, expression preparation, and display may differ from GammaLoop. CLI-only regressions or QoL gaps do not block GammaLoop readiness unless they expose shared CFF mathematics used by production. |
| Other split material | **Preserved outside this WIP** | The NNLO example-card delta and proper-LTD comparisons are still source material. Review them against the final architecture instead of copying their historical revisions wholesale. |

The phrase “rebase complete” refers to the conflict-free replay. It does not
mean the raised-energy feature branch is ready to merge.

## Locked architecture decisions

These decisions were discussed and approved. Do not reopen or silently bypass
them while resolving the remaining work.

1. There is one public CFF implementation, owned by the
   `three-dimensional-reps` crate. GammaLoop owns production graph/source
   construction, UV orchestration, exact mapping, and evaluator preparation.
   The `3Drep` command/evaluator pipeline is diagnostic and need not prepare the
   same inputs or expressions.

2. Every CFF power and capacity question is answered solely in physical
   EMR/source-edge energies; numerator support is represented as
   `Vec<(edge_id, degree)>`. LMB coordinates are routing data and are never
   consulted as an energy identity, capacity index, or fallback. Exact/local
   sources must translate physical EMR capacity explicitly. Production analysis
   accepts `Q(edge, index)` and rejects `K(loop, index)` until an upstream
   provenance-aware producer normalizes it to a physical edge.

3. Preserve main's factorized, incremental numerator ownership. Scaling
   analysis may inspect the factorized numerator first, but 3D UV recursion
   must continue to grow it from spinney-owned numerator fragments. Preserve
   parametric energy signs and shifts; do not expand/concretize tensors as a
   shortcut. Higher-power interpolation may substitute the relevant EMR energy
   by signed integer `a*M`; `M` is one common auxiliary scale, is required to be
   nonzero only when a finalized evaluator uses it, and the result is invariant
   under changing its nonzero value. This substitution never uses LMB variables.

4. There are two UV routes. Both retain the complete 4D integrated
counterterms. The local subtraction then either:

   - completes the loop-energy integration first and applies every UV Taylor
     operator to the complete/global CFF expression. Its localized and
     explicit-sum forms use the same transformed bodies and differ only because
     the latter omits orientation selectors; or
   - performs 4D-local UV first and projects each finished proper UV term to CFF
     after its required series and limits. Exact-source reconstruction and
     minimax EMR dispatch are exclusive to this projected route.

   The empty forest is the ordinary factorized production root in both routes;
   the expanded-4D setting governs only proper, nonempty UV approximations.

5. The projected 4D-to-CFF route permits only an explicit source-local sum over
   all orientations. It must reject individual orientation selection, patterns,
   and orientation sampling. Production `OrientationID`s are not expected to
   match exact-source IDs. Keep the sum deferred through evaluator-style
   function maps to limit expression growth.
   Per-orientation UV profiling is correspondingly meaningful only for the
   orientation-parametric localized direct-3D form. Selector-free explicit-sum
   direct 3D and projected local 4D must reject that profiling mode.

6. Exact projection keeps its identity layers separate. Occurrence-local exact
   energy IDs serve source algebra, while physical EMR edge IDs serve initial
   numerator-capacity analysis and physical provenance; neither can be
   substituted for the other. The original `EdgeIndex` survives the Taylor
   operator as `source_edge`. Disjoint-set contractions of absent original
   edges construct the cograph and UV source minors, so that `source_edge`
   directly determines incidence. Only repeated occurrences of that
   owner subdivide its incidence into a serial power chain; distinct owners are
   never joined by matching signatures. Exact signatures are canonicalized up
   to sign only to validate a source-local denominator and expose algebraically
   repeated CFF channels. A quotient-space rank solve selects only the unique
   `+/-` routing sign for already fixed endpoints; it is not an endpoint solver.
   The raw sign is composed back in when mapping an occurrence energy to the
   physical numerator convention. After incidence and source-crown hedges are
   complete, Graphica canonically relabels the graph for order-independent cache
   reuse. No incidence/Kirchhoff reconstruction or balance-edge synthesis is
   present. Owner IDs
   remain available for cuts, physical surfaces, and parent-graph mapping. One
   immutable minimax plan owns both occurrence-local CFF bounds and factorized
   numerator substitutions. Projection plans every Taylor term before
   generation. Equal canonical topologies share a capacity envelope formed by
   the maximum total degree of each repeated algebraic energy channel, followed
   by minimax redistribution; non-repeated occurrence bounds use componentwise
   maxima. Each genuinely outer term nevertheless evaluates with its own plan.
   The real DOD1 triangle remains one common-denominator exact source with
   UV-owner multiplicities `(2,1,2)` and one canonical generator entry. Its
   positive typed `GS.den` numerator factors remain uncancelled, and generalized
   CFF supplies the base and singly dotted lower sectors internally. Pure-external
   non-vacuum crowns are implemented; a future OS two-point insertion such as
   `(m,0,0,0)` needs only an explicit fixed-boundary payload on the same
   source-built topology. The complete design is in
   [`exact-powered-denominator-cff-lifting.md`](exact-powered-denominator-cff-lifting.md).
   In direct-3D UV, a production `OrientationID` is the opaque selector for one
   complete generalized residue-map key, while
   an optional source-local `edge_energy_map` remains authoritative for every
   incrementally attached numerator factor. A contracted contact map need not
   be an affine extension of any production map. The projected-4D route remains
   the separate explicit source-local sum described above. The shared CFF
   global-prefactor metadata stays attached to each generated source and
   GammaLoop consumes it exactly once for root, reduced, and exact production
   CFF sources.

7. Structural non-cut tree/bridge edges stay outside the CFF graph. GammaLoop
   projects them separately to 3D in caller-owned `tree_denoms`, as main does.
   The 3D route still has shrunken-edge semantics; do not erase that concept.

8. Use main's Symbolica revision and main's parameter/command ordering. Keep
   useful QoL changes, but do not retain dead alternate implementations.

9. LTD is deferred, not discarded. Do not restore the old hybrid residue,
   H-surface, or second “confluent CFF” implementation under a hidden code path.

10. A64's future numerator-call deduplication should key production identity by
    canonical authoritative `edge_q0` and retain one deterministic `loop_q0`
    diagnostic. It was explicitly deferred and is not a rebase blocker.

## Generalized-CFF inherited-contour sign

The generalized-CFF recursion has one sign rule which belongs entirely to the
`three-dimensional-reps` crate. GammaLoop must neither reconstruct nor
compensate for it.

When a powered-pole contact removes a contour coordinate, the recursive lower
sector can factor into `C` rational denominator components. The lower-sector
builder evaluates those components independently and multiplies their CFFs.
Relating that product of independently closed contours to the parent contour
requires the relative component-product convention

\[
  (-1)^{C-1}
\]

exactly when the lower sector starts, or inherits, that component-product
frame. At the top level it always starts a new frame. For a recursive contact,
the deleted contour rows are stored from outermost to innermost. The inherited
frame survives precisely when the **outermost** deleted row still has a carrier
among the surviving denominator rows. If that row has no surviving carrier,
the parent contact has completely consumed that contour direction; an inner
row may fix a component-local closure, but it cannot reopen the already
consumed product frame.

This predicate deliberately ignores the sign of the carrier. Reversing a
routing changes the Above/Below closure and the residue Jacobian together;
those are component-local oriented-residue data already applied separately.
It does not change whether the denominator still carries the deleted
direction, and therefore cannot change the component-product bridge.

The rule corrects an assumption introduced in `d54de2049`: a deleted row lying
wholly in one powered component was assumed not to need the bridge because the
component's public CFF supposedly owned the relevant duplicate parity. That
conflated two different signs. Duplicate-denominator parity belongs to the
component's own rational CFF, whereas `(-1)^(C-1)` relates the product of
component contours to its parent contour. A powered component can therefore
need both. The correction is internal to generalized CFF and must not be
mirrored by a GammaLoop "production bridge".

The protected direct-core oracles include:

- `repeated_cubic_pole_with_squared_denominator_numerator_pinches_to_simple_pole`,
  which checks `D(Q)^2 / D(Q)^3 = 1 / D(Q)` without a GammaLoop convention;
- the repeated-quintic lowering tests, including multiple pinches, numerator
  sampling scales, and routing reversals. The mixed-cograph oracle checks
  `D(Q) / D(Q)^5 = 1 / D(Q)^4` for both disconnected and incidence-connected
  two-component products, with both `Q` and `-Q` routing;
- the lower-sector component-product and nested inherited-contour tests, which
  separately cover an outer surviving carrier, its routing reversal, and an
  outer row with no surviving carrier.

An isolated A/B check at the pre-correction `c8e763173` source boundary is a
particularly useful diagnostic. With only the incidence-connected,
oppositely-routed mixed-cograph fixture added, the old predicate returned
`-1/9` times the independently generated quartic lower channel at all three
kinematic points. Applying only the inherited-component bridge correction made
the same test agree locally. This establishes that the correction acts inside
generalized-CFF recursion rather than through GammaLoop prefactor handling.

The historical CFF-versus-LTD suite remains the strongest future independent
regression gate when proper LTD is ported. It validated the predecessor
higher-power implementation, but it predates `inherited_contour_rows` and
`requires_inherited_component_bridge`; it must therefore not be cited as proof
of this post-`d54de2049` correction. When LTD support is absorbed, rerun at
least `ltd_multiloop_high_power_boundaries_match_cff` and its independent
`3Drep test-cff-ltd` counterpart before changing this rule again.

The theory reference used for the audit is
[`generalised_ltd.tex` at revision `65ab03e4`](https://github.com/alphal00p/generalised_ltd/blob/65ab03e4fb7d442ff362392df2c2a59ef323d28c/docs/generalised_ltd.tex).
In particular, CFF retains denominator powers and applies repeated-channel
normal form; the numerator remains a black box with one argument and degree cap
per original EMR edge.

## GammaLoop CFF normalization boundaries

The inherited-contour rule above is internal to generalized CFF. A separate,
typed boundary converts each **freshly generated source-local rational
component** to GammaLoop's source-denominator convention. For a component `c`
with `N_c` denominator occurrences, GammaLoop applies

\[
B_c =
\begin{cases}
B_{\mathrm{core},c},
  & \text{global source product},\\
(-1)^{N_c} B_{\mathrm{den},c} B_{\mathrm{core},c},
  & \text{variant-local source}.
\end{cases}
\]

The complete fresh-source conversion is `product_c B_c`. The three factors have distinct
origins: `B_core` is the generalized-CFF shared contour convention,
`B_den` is the scalar denominator-source convention, and `(-1)^N` converts the
variant-local positive `1/(2E)` factors to GammaLoop's `1/(-2E)` convention.
`EmbeddedCffFactor` retains one deterministic representative when a terminal
contour would otherwise be duplicated; that changes the catalogue only.  The
retained representative therefore uses the same component bridge as the
corresponding standalone residue.  This is a hard invariant: generation
context may not enter the production-prefactor bridge.  In particular, a
repeated source such as `D(Q)^-5` is nonterminal, so standalone and embedded
generation both call the same bounded-CFF builder and return byte-for-byte the
same raw rational functional.  Assigning those identical expressions opposite
signs merely because one is later multiplied by an outer graph has no contour,
measure, or BPHZ justification.
For a one-loop repeated channel, this fresh-source bridge is invariant under algebraic
lowering.  At `N=2` it is `(+1)(-1)(+1)=-1`; at `N=3` it is
`(-1)(+1)(+1)=-1`; and the same alternation continues for quartic and quintic
sources.  The occurrence conversion and denominator-frame parity cancel each
other's change when one copy of `D(Q)` moves between numerator and denominator.
Consequently neither the aggregate core sign nor the aggregate denominator
sign is a valid **fresh-source** bridge on its own.

This rule must not be applied a second time when consuming the persisted
production root. That object already is the graph's complete production
residue functional. Its `GlobalSourceProduct`/`VariantLocal` field records
where the energy factors live; it does not request another denominator-frame
conversion. The stored-root boundary therefore consumes exactly

\[
  B_{\mathrm{stored\ root}} = B_{\mathrm{core}}.
\]

The distinction is by object lifetime and ownership, not by a physics mode:

- a newly reconstructed exact/local-4D source crosses the componentwise
  source-convention bridge once;
- a stored production root has already crossed its construction boundary and
  reuses only its recorded core sign when materialized for direct post-CFF
  Taylor subtraction.

Do not use the fresh-source helper in `cff_from_production_expression`. Doing
so double-converts the stored residue frame. Conversely, do not remove the
componentwise bridge from exact 4D sources; that reintroduces the powered-source
sign error described below.

An exact clean-tree A/B at `c8e763173` freezes why the fresh-source bridge belongs
to the GammaLoop adapter rather than to generalized-CFF recursion. With a
core-only **fresh exact-source** adapter, the already-existing uncancelled powered-source
oracle evaluated

```text
D(Q) * (Q0 + c) / D(Q)^3  ->  +2.015720902074968...e-3 i
(Q0 + c) / D(Q)^2         ->  -2.015720902074968...e-3 i
```

at the exact point `|q|=0`, `E=1`, `c=2`: the values were exact opposites at
Arb precision.  The raw generalized-CFF lowering identities were already
green.  Applying the component bridge above makes this production comparison
green, as well as its scalar, globally reversed-routing, quintic-to-quartic,
and raised-LU-residue variants.  Thus the gross powered-source sign was a
latent ownership-conversion defect at the GammaLoop boundary: ordinary CFF
clears local half-edge factors and appends `prod 1/(-2E)`, whereas generalized
variants retain `prod 1/(2E)`.  It is not evidence for changing the residue
map or its contact recursion.

The fresh-source fold must remain component-local. A disconnected product can mix a
variant-local powered bubble with a global ordinary component, and aggregate
ownership cannot reconstruct that mixture. `GeneratedThreeDExpression`
therefore persists `energy_factor_components`; its source-bound report remains
transient analysis metadata. This changes the positional bincode layout, so
states written before this metadata was persisted must be regenerated rather
than decoded as the new format. Neither GammaLoop boundary rule may be copied
back into generalized CFF, and neither may be conflated with the internal
`(-1)^(C-1)` inherited-component rule.

The minimal stored-root regression is the scalar self-energy `T0` source with
the untouched outer numerator `Q3^0 Q4^0`. Its production root has five
variant-local denominator occurrences with `B_den=+1` and `B_core=-1`.
Applying the fresh-source formula would give
`(-1)^5 B_den B_core=+1`, while the stored-root contract gives `-1`.
Changing only that consumer boundary to the fresh-source value makes the
assembled direct result the exact negative of the independently reconstructed
child-times-outer result. In this deliberately canonical `T0` fixture, keeping
the stored core sign also happens to make the separately stored LU-order pieces
agree bit-for-bit at 1000-bit Arb precision, for both the values and their
required first Taylor derivatives. That fixture-specific stronger equality is
not the route contract: in general only the sum of all raised-order derivative
pieces and residue-map entries for the fixed cut is invariant. The raw
generalized residue family and projected child remain unchanged in this A/B,
which locates the defect in GammaLoop's stored-root adapter rather than in
generalized-CFF recursion.

## Current NLO normalization acceptance

A79 supersedes A61's unique-literal restriction. Original edge IDs survive the
UV Taylor operator and determine the contracted UV skeleton directly. Raised
denominators become serial dotted edges belonging to that owner; canonical
`+Q/-Q` signatures only validate algebraic channel equivalence. The still
factorized numerator is then dispatched over occurrence-local EMR energies with
a deterministic minimax plan. No LMB coordinate is used as an energy identity,
and no Kirchhoff, Laplacian, or signature-to-incidence graph reconstruction is
performed.

The retained 2026-08-31 10x campaign completes all four physical DD/TT
acceptances (`4/4`). Direct-current tests use the summed Feynman-gauge projector
`-g^{mu nu}` and no unit conversion; the `e+e-` tests use the published
direct-current conversion and inclusive ratio. Pulls are signed differences
from the published target in units of the Monte Carlo error, with LO uncertainty
included in ratio pulls.

| Acceptance | 10x LO result | 10x NLO result | Graph and ratio evidence |
| --- | --- | --- | --- |
| direct `gamma* -> d d~` | `0.5068703962 +/- 0.0025987972` (`+1.449 sigma`) | `0.01966009810 +/- 0.00053595339` (`+1.424 sigma`) | `GL0=-0.03132123586 +/- 0.00023922726` (`+0.729 sigma`), `GL2=+0.05112479005 +/- 0.00046213299` (`+1.584 sigma`); `alpha_s/pi` pull `+1.141 sigma` |
| converted `e+e- -> gamma* -> d d~` | `0.1950499744 +/- 0.0010326753 pb` (`+1.479 sigma`) | `0.007824189766 +/- 0.000340601513 pb` (`+1.630 sigma`) | signed MC components `GL0=-0.01996339254 +/- 0.00015479207`, `GL2=+0.02786745983 +/- 0.00028425324`; no separate published component targets; `alpha_s/pi` pull `+1.453 sigma` |
| direct `gamma* -> t t~` | `2.901968994 +/- 0.015639978` (`+1.641 sigma`) | `0.2079169992 +/- 0.0042953541` (`+1.489 sigma`) | `GL0=-0.1443600613 +/- 0.0035809931` (`+0.669 sigma`), `GL2=+0.3522770605 +/- 0.0023720361` (`+1.687 sigma`); paper-ratio pull `+1.037 sigma` |
| converted `e+e- -> gamma* -> t t~` | `0.3307052414 +/- 0.0018004843 pb` (`+1.603 sigma`) | `0.02356890542 +/- 0.00056839205 pb` (`+1.058 sigma`) | summed-graph integration has no persisted GL0/GL2 rows; paper-ratio pull `+0.685 sigma` |

All LO integrations used 100,000 samples. The direct DD NLO used 400,000,
the direct TT graph rows used 400,000 each, and each converted NLO central slot
used 200,000. The converted DD acceptance also verifies pointwise equality of
all three routes, the GL0 quadratic and GL2 linear EMR bounds, per-graph `M_UV`
and localization-scale independence, and opposite GL0/GL2 `mu_r` logarithms.
Its graph rows are not separately Ward-transverse and therefore are not assigned
invented published targets.

The auxiliary `3Drep` input/evaluation pipeline remains diagnostic only and is
not a production architecture to imitate.

The protected GL11 `Q1^2` and `Q7^2` fixtures are a separate test-intent issue.
Their stored values are respectively
`-1.17703539069413922e-13 - 2.37175489451059046e-30 i` and
`-2.22459688841192310e-11 - 4.48261675062501648e-28 i`. The validated baseline
cleared numerator-energy bounds before constructing those roots, so the values
are denominator-only cancellation residues rather than generalized raised-EMR
oracles. Do not edit either fixture or snapshot without explicit approval.

## Supporting GammaLoop work

### Integration resume boundary (A60)

GammaLoop writes `<workspace>/manifest.json` and
`<workspace>/state/integration_state.bin`. `IntegrationState` contains
`StatisticsCounter` directly and in `slot_stats`; adding
`numerical_stability` changed the nested positional bincode layout. The
independent integration-workspace boundary is now implemented: new manifests
write version 1, and both the normal and summary resume paths validate that
version before decoding state. Missing, older, and future versions return
contextual errors with `integrate --restart` guidance. There is deliberately no
positional migration for unversioned workspaces.

### Remaining product checks

- A66: add fixed serialized-row and fixed-kinematics numerical CFF goldens from
  the external reference. Do not alter existing tests or fixtures without
  approval.
- A67: the mechanical schema, live-example, explicit read-only-output, Python
  API migration, and BNL-documentation repairs are complete. The four
  integration-helper dead-code suppressions remain necessary while their only
  caller is inside the legacy commented method. Default output placement, the
  late-parsing benchmark fixture, and the wider benchmark suite remain
  unchanged because they require test or tooling-policy decisions.
- A65's conditional zero-scale validation and nonzero-scale invariance are
  complete. The existing degree-one `All` activation remains deliberately
  unchanged because altering that test-backed policy is not required for
  correctness.
- A72 retains the no-op `serde=[]` feature as documented downstream
  compatibility vocabulary; serde dependencies and derives remain
  unconditional.
- Archive coverage: the amplitude JSON archive now has a typed reload and real
  eager-evaluation test. The cross-section export/loader smoke test still lacks
  the equivalent typed reload/evaluation coverage.
- Computed UV-forest export: choose between regenerating the production exact
  CFF catalogue and persisting exact map IDs. Do not infer affine shifts from
  coarse stored orientations.
- A69 records the lost pre-rebase auxiliary `3Drep graph-from-signatures`
  command. Restore it only in a separate diagnostic-tooling task if desired; it
  is not GammaLoop readiness work.
- A68/A71 and other auxiliary preparation/display differences are nonblocking
  unless they expose shared CFF mathematics or affect a GammaLoop value.

## Validation already completed

The following checks passed at their documented rebased checkpoints.
Symbolica-backed tests must be run serially: concurrent identical binaries have
hit the restricted license/global-initialization boundary and aborted, while
the same tests pass individually. Rerun the final suite after the newer exact
carrier and production-prefactor changes before declaring readiness.

```bash
nix develop -c cargo fmt --all -- --check
nix develop -c cargo check --locked -p gammalooprs --lib
nix develop -c cargo metadata --locked --offline --no-deps --format-version 1
nix develop -c cargo hakari verify
```

Focused passing coverage includes:

- `deferred_explicit_sum_lowering_is_local_and_matches_materialized_value`;
- `deferred_explicit_sum_preprocesses_tensor_function_bodies`;
- `integrands::process::tests::explicit_orientation_sum_rejects_runtime_filters_and_ltd`;
- `test_integrated_uv_cts::orientation_monte_carlo_rejects_non_parametric_evaluators`;
- `exact_cff_keeps_dotted_same_edge_occurrences`;
- `exact_cff_handles_opposite_repeated_routing_without_a_sign_bridge`;
- `term_projection_keeps_typed_numerator_factors_uncancelled`;
- `exact_source_routes_a_hard_uv_row_modulo_the_soft_cograph_span`;
- `dod_one_triangle_keeps_one_uncancelled_exact_source_and_matches_lower_sectors`;
- the GammaLoop amplitude standalone export, typed JSON reload, and real eager
  evaluation test;
- the cross-section archive/loader export smoke test;
- integration-workspace version rejection before both resume decoders;
- conditional zero-`M` rejection and higher-power invariance under two nonzero
  sampling scales;
- an independent generated amplitude-loader type check.

The latest focused production evidence after the lower-sector CFF repair is:

- the scalar LU 15-case matrix exercises all three local-UV routes at f64 and
  native Arb with local UV, integrated UV, and threshold counterterms enabled.
  Its unit-numerator graphs are left as generated and companion numerator lanes
  attach only Feynman-rule-local edge factors. There is no process-name or
  graph-specific production branch. As retained pre-reversal evidence, the
  2026-08-31 `dev-optim` / `test_gammaloop` rerun passed `15/15`, with `235`
  skipped, in `70.438 s`. The restored post-CFF direct route still requires a
  current matrix rerun. Four near-zero cases
  use the authorized f64-input `1e-14` unit-scale fallback; Arb-to-Arb remains
  run and reports non-scaling. This is test-oracle handling only and required no
  production change. The curated suite also selects the focused DOD0/1/2
  orientation-local bubble regression. All six base scalar-matrix graph tests
  (GL00, GL02, GL04, GL08, GL09, and GL24) invoke localized-direct3D
  per-orientation profiling; their current post-restoration run is pending;
- the direct and converted DD/TT physical acceptances pass together (`4/4`) in
  the retained 10x campaign, with every absolute and ratio pull below `1.7
  sigma`; the exact results are recorded in the table above;
- the end-to-end UV-profile CLI tests pass (`2/2`), and the lower-level
  amplitude/LU selection and fitting tests pass (`15/15`), observed at 1.213
  and 0.660 seconds respectively. The default
  `only-divergent` mode profiles every expected cycle union with DOD >= 0,
  using the generation LMB when suitable and otherwise the first suitable LMB
  in the deterministically sorted complete list. `all` remains the exhaustive
  opt-in mode; graph and Cutkosky-cut selectors are supported for LU input.

The reported GammaLoop generation command has reached completion at documented
checkpoints:

```bash
./gammaloop --release --clean-state \
  ./examples/cli/epem_a_ddx/NLO/epem_a_ddx_NLO.toml \
  run generate_diagrams generate_integrands -c "quit -n"
```

It generates GL0 and GL2 with five evaluators each. This generation-only card
is not by itself the inclusive NLO proof; the direct-current and converted-LU
acceptances above provide that numerical evidence. The current card has
`generate_integrated = true`, so it exercises direct-3D local subtraction plus
the integrated finite UV addback and its exact source-local EMR carrier. The
default `local_uv_cts_from_expanded_4d_integrands = false` means this does
**not** validate the alternate 4D-local-to-CFF route.

The current full shared-CFF crate inventory is `98/100`. Its two remaining
failures are stale test-contract assertions rather than a demonstrated
GammaLoop numerical regression:

- `eval::tests::lower_sector_powered_pole_contact_reconstructs_numerator_derivatives`;
- `generation::causal_generation_tests::inherited_mixed_theta_basis_uses_correlated_parent_residues`.

The powered-pole fixture needs a physical/unimodular value-oracle redesign; the
mixed-theta assertion still names a superseded private selector origin while
its retained value comparison agrees. Change those tests only after confirming
their intended contract, as required by repository policy.

## Safe continuation sequence

1. Treat remote commit `5c3e7f` as the published checkpoint, confirm
   `jj log -r 'conflicts() & ancestors(@)'` is empty and inspect the complete
   local continuation above that checkpoint and its `main` base.
2. Preserve the four passing DD/TT three-route normalization acceptances, their
   pointwise local-3D/local-4D equivalence checks, and the DD scale invariants;
   leave phase and overall sign outside the criterion for now.
3. Keep every higher-power capacity, finite-pole sample, and integer-multiple
   `a*M` substitution inside EMR. Preserve the factorized numerator and do not
   reproduce the diagnostic `3Drep` input/evaluation pipeline in GammaLoop.
4. Do not change the protected GL11 tests or snapshots while their
   fixture/oracle intent is
   unresolved. If approval is obtained, run its battery serially:

   ```bash
   nix develop -c cargo test --locked -p gammaloop-integration-tests \
     --test test_runs \
     scalar_3l_cross_section_gl11_quadratic_energy_inspects_match::q7_squared \
     -- --test-threads=1
   ```

5. Rerun the scalar LU 15-case all-enabled matrix as a final gate for the
   restored post-CFF direct route. Its pre-reversal 2026-08-31 run passed
   `15/15`, with `235` skipped, in `70.438 s` under `dev-optim` /
   `test_gammaloop`; retain that result as historical evidence. If the current
   matrix regresses, first deconstruct the
   generated expressions and use the generalized `3Drep` CLI only as the
   independent 3D diagnostic lane; do not use a handwritten higher-loop contour
   calculation as an oracle.
   Resolve the two shared-CFF stale test contracts, then rerun the physical
   acceptances plus `just test_gammaloop`.
6. Record every architectural finding, approval, test, and remaining failure in
   `kysvnqlq-rebase-review.md`, then refresh this operational summary if its
   status changes.
7. Fetch/rebase the final stack on the latest `main`, repeat the conflict and
   artifact audit, and push only the intended bookmark. Do not include the
   unrelated conflicted bookmark.

## Things not to infer

- “Conflict-free” does not mean “ready to merge.”
- A passing direct-3D GammaLoop command does not prove projected 4D raised
  support.
- `source_edge`, an occurrence ID, an orientation ID, a loop carrier, and a
  physical EMR energy ID are not interchangeable.
- The factorized EMR rank envelope provisions generalized-CFF numerator
  reconstruction; it is not an overall loop-energy convergence test. GammaLoop
  separately rejects a source graph whenever the existing all-cycle UV DOD
  analysis gives `DOD_E = DOD_4D - 3 L >= 0`. This changes only the loop
  measure from four powers to one per energy integration and does not introduce
  a bespoke exact-source contour certificate.
- The auxiliary `3Drep` CLI is diagnostic only. Its input and expression
  preparation need not match GammaLoop, and CLI-only regressions do not block
  production unless shared CFF mathematics or a GammaLoop value is affected.
- Proper LTD is preserved for later work; its deferral is not permission to
  delete the setting, source reserve, tests, or architectural intent.

## Stale source warnings

Do not resolve, merge, or branch from the conflicted local
`ltd_in_gammaloop_symbolica_update` line (`swnpzvnu` / `f46f9f56`), its older
remote source (`cc8c5979`), or the remote historical monolith
`ltd_in_gammaloop_symbolica_update_monolith` (`uksvkwzw` / `1b73243b`). The
detached old ledger change `zzxptpwq` is historical as well. Relevant intent is
already captured in the live ledger or the pushed `proper_ltd_core` reserve.
