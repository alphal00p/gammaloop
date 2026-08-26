# Raised-energy CFF rebase handoff

Last updated: 2026-08-26

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

The last published remote bookmark `raised_energy_cff_wip` is the starting point
for the newer local continuation, not its current implementation tip. Until the
final stack is cleaned, validated, and pushed, fetching it on another machine
recovers only that published starting point. Use it rather than reviving one of
the old conflicted source revisions:

```bash
jj git fetch --remote origin
jj bookmark track raised_energy_cff_wip@origin
jj new raised_energy_cff_wip
```

If the bookmark is already tracked but points elsewhere, first inspect the
divergence and then set it explicitly to `raised_energy_cff_wip@origin`. Do not
merge a hidden, stale remote tip by assumption.

The published starting point for this continuation was:

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

The local non-LTD continuation now extends that starting point with exact EMR
capacity mapping (`ozntmkwv`), typed stored production roots (`toxllnzs`),
integration-workspace versioning (`onotmrxw`), sampling-scale validation
(`qyprnryk`), interface maintenance (`qykwxoxl`), typed amplitude archive reload
coverage (`trrvmwrw`), stored-root selector preservation (`zroqlyxz`), mapped
numerator scalarization regression coverage (`uukzotwz`), the initial production
sign bridge (`yyksullr`), exact source-local EMR carriers through integrated
finite-UV localization (`opmnpruw`), and the full typed CFF global-prefactor
bridge (`kypmvxsq`). Focused coverage and successful generation at these
boundaries do not constitute the standalone inclusive NLO acceptance introduced
under `skvrlunw`; that acceptance remains open.

Before editing, read `AGENTS.md`, `CONTRIBUTING.md`, this handoff, and the full
decision ledger. Use Jujutsu and start a new described change for each approved
unit of work.

## Status at a glance

| Area | Status | Continuation meaning |
| --- | --- | --- |
| Rebase onto current `main` | **Complete** | Fourteen retained changes are replayed directly onto `39561014`; the selected ancestry is conflict-free. The obsolete `kysvnqlq` pin was abandoned because it only downgraded main's Symbolica stack. |
| Latest-main evaluator/archive overlap | **Resolved and validated** | Main's Symbolica, native-code, O2, evaluator construction, and command ordering are retained. Deferred explicit CFF orientation sums are layered onto that architecture. |
| Direct-3D raised-energy CFF | **Implemented; protected historical oracle conflicts with generalized semantics** | Physical EMR/source-edge bounds and the repeated-channel CFF construction are restored. The protected GL11 `Q1^2`/`Q7^2` snapshots remain unchanged; A73 records why their denominator-only cancellation residues are not generalized-CFF value oracles. |
| 4D-local UV then CFF with raised energy | **Unique-literal EMR mapping restored for proper UV nodes** | A61 maps physical EMR bounds only through one literal exact occurrence and rejects missing, ambiguous, shifted, or affine cases. The empty forest remains the ordinary factorized production root in both routes; only proper UV nodes take the expanded-4D projection route. |
| Integration-workspace persistence | **Implemented and focused-validation complete** | Integration workspaces have an independent version-1 manifest checked before both resume decoders. Incompatible workspaces return contextual `integrate --restart` guidance instead of decoding positional bincode or panicking (A60). |
| Other GammaLoop completeness | **NLO normalization acceptance still open** | Conditional zero-`M` validation, nonzero sampling-scale invariance, exact source-local EMR carrier preservation, route-independent empty-root identity, schemas, live example pins, explicit read-only guards, package maintenance, and typed amplitude archive reload/evaluation are implemented. The production-prefactor bridge has focused coverage, but no current inclusive acceptance pass has been recorded. Fixed external goldens, the wider benchmark-suite policy, and typed cross-section reload/evaluation remain separate gaps. |
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

4. There are two UV routes. Both build the complete 4D local and integrated
   counterterms. The local subtraction then either:

   - runs the direct-3D UV recursion on CFF expressions, retaining localized
     orientation-parametric subtraction; or
   - performs 4D-local UV first and projects each finished proper UV term to CFF
     after its required series and limits.

   The empty forest is the ordinary factorized production root in both routes;
   the expanded-4D setting governs only proper, nonempty UV approximations.

5. The projected 4D-to-CFF route permits only an explicit source-local sum over
   all orientations. It must reject individual orientation selection, patterns,
   and orientation sampling. Production `OrientationID`s are not expected to
   match exact-source IDs. Keep the sum deferred through evaluator-style
   function maps to limit expression growth.

6. Exact projection keeps its identity layers separate. Occurrence-local exact
   energy IDs serve source algebra, while physical EMR owner IDs serve
   capacity/ownership/topology; neither can be substituted for the other. In
   direct-3D UV, a production `OrientationID` partitions theta sectors while an
   optional source-local `edge_energy_map` remains authoritative for every
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

The theory reference used for the audit is
[`generalised_ltd.tex` at revision `65ab03e4`](https://github.com/alphal00p/generalised_ltd/blob/65ab03e4fb7d442ff362392df2c2a59ef323d28c/docs/generalised_ltd.tex).
In particular, CFF retains denominator powers and applies repeated-channel
normal form; the numerator remains a black box with one argument and degree cap
per original EMR edge.

## Immediate continuation: NLO normalization acceptance

A61 is implemented. `ExactSourceEnergyMapper` translates a physical bound only
through one literal `Q(edge)` occurrence; missing, repeated, shifted, and
general-affine cases still hard-error. This preserves the factorized numerator
and does not use LMB coordinates or ownership as energy identity.

The remaining production acceptance is the standalone test based on
`examples/cli/epem_a_ddx/NLO`. Its inclusive GL0+GL2 magnitude must reproduce
`(alpha_s/pi) * LO`; phase and overall sign are deliberately not acceptance
criteria yet. The live card enables integrated UV generation while retaining
the direct-3D local-subtraction route. It therefore exercises exact source-local
EMR carrier preservation and the one-time GammaLoop production-prefactor
bridge, but not the alternate 4D-local-to-CFF projection. Those implementation
boundaries have focused coverage; the inclusive acceptance has no recorded
passing result and remains open. The auxiliary `3Drep` input/evaluation pipeline
is diagnostic only and is not a production architecture to imitate.

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
- the GammaLoop amplitude standalone export, typed JSON reload, and real eager
  evaluation test;
- the cross-section archive/loader export smoke test;
- integration-workspace version rejection before both resume decoders;
- conditional zero-`M` rejection and higher-power invariance under two nonzero
  sampling scales;
- an independent generated amplitude-loader type check.

The reported GammaLoop generation command has reached completion at documented
checkpoints:

```bash
./gammaloop --release --clean-state \
  ./examples/cli/epem_a_ddx/NLO/epem_a_ddx_NLO.toml \
  run generate_diagrams generate_integrands -c "quit -n"
```

It generates GL0 and GL2 with five evaluators each. This proves generation, not
the still-open inclusive NLO normalization acceptance. The current card has
`generate_integrated = true`, so it exercises direct-3D local subtraction plus
the integrated finite UV addback and its exact source-local EMR carrier. The
default `local_uv_cts_from_expanded_4d_integrands = false` means this does
**not** validate the alternate 4D-local-to-CFF route.

The feature-gated shared-core diagnostic inventory is 53/59. Its six failures
are classified under A59 in the detailed ledger. That count is not a GammaLoop
product gate; change an assertion only after establishing regression and test
intent, and ask first as required by repository policy.

## Safe continuation sequence

1. Until the bookmark is moved to the final tip, confirm
   `jj log -r 'conflicts() & ancestors(@)'` is empty and inspect the complete
   local continuation above its `main` base.
2. Finish the standalone `epem_a_ddx` GL0+GL2 acceptance and require its
   inclusive magnitude to agree with `(alpha_s/pi) * LO`; leave phase and
   overall sign outside the criterion for now.
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

5. Rerun the exact `epem_a_ddx` generation command, its standalone acceptance,
   and the focused validation above serially.
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
- A spatial UV profile after generation cannot prove the loop-energy contour
  condition that permits CFF generation.
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
