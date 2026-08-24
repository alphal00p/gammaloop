# Raised-energy CFF rebase handoff

Last updated: 2026-08-24

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

The working branch is the remote bookmark `raised_energy_cff_wip`. On another
machine, fetch it rather than reviving one of the old conflicted source
revisions:

```bash
jj git fetch --remote origin
jj bookmark track raised_energy_cff_wip@origin
jj new raised_energy_cff_wip
```

If the bookmark is already tracked but points elsewhere, first inspect the
divergence and then set it explicitly to `raised_energy_cff_wip@origin`. Do not
merge a hidden, stale remote tip by assumption.

The code tip immediately before this documentation-only handoff is:

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

Before editing, read `AGENTS.md`, `CONTRIBUTING.md`, this handoff, and the full
decision ledger. Use Jujutsu and start a new described change for each approved
unit of work.

## Status at a glance

| Area | Status | Continuation meaning |
| --- | --- | --- |
| Rebase onto current `main` | **Complete** | Fourteen retained changes are replayed directly onto `39561014`; the selected ancestry is conflict-free. The obsolete `kysvnqlq` pin was abandoned because it only downgraded main's Symbolica stack. |
| Latest-main evaluator/archive overlap | **Resolved and validated** | Main's Symbolica, native-code, O2, evaluator construction, and command ordering are retained. Deferred explicit CFF orientation sums are layered onto that architecture. |
| Direct-3D raised-energy CFF | **Implemented; final historical oracle still pending** | Physical EMR/source-edge bounds and the repeated-channel CFF construction have been restored. The reported `epem_a_ddx` GammaLoop generation command passes. Rerun the historical GL00 value before changing any snapshot (A73). |
| 4D-local UV then CFF with raised energy | **Known pre-rebase regression; next priority** | GL00 `Q7^2` stops at the conservative A52 exact-source guard. The narrow old capability is understood, but no repair has been approved or implemented (A61). |
| Integration-workspace persistence | **Product issue, decision pending** | The numerical-stability field changed positional bincode data without versioning the independent integration workspace. This predates the replay and is not a rebase conflict, but it affects GammaLoop resume behavior (A60). |
| Other GammaLoop completeness | **Open** | Uniform-scale validation, fixed numerical goldens, committed schemas, read-only behavior, benchmark adaptation, and typed archive round-trip coverage remain. |
| Computed UV-forest export | **Deferred; not exact yet** | Persisted graph terms retain only coarse orientations. Production generation/evaluation already uses exact maps, but export must regenerate the production catalogue or persist exact map IDs before it can claim exact affine-map identity. |
| Proper LTD | **Explicitly deferred** | Keep the `Cff`/`Ltd` selection; `Ltd` must immediately return the typed unsupported error. Port proper LTD later from the sibling reserve onto the finalized CFF plus 4D architecture. |
| Auxiliary `3Drep` CLI | **Not the product gate** | Only demonstrated regressions against the selected pre-rebase baseline matter here. `graph-from-signatures` is such a regression (A69), but it is below GammaLoop product work. |
| Other split material | **Preserved outside this WIP** | The NNLO example-card delta and proper-LTD comparisons are still source material. Review them against the final architecture instead of copying their historical revisions wholesale. |

The phrase “rebase complete” refers to the conflict-free replay. It does not
mean the raised-energy feature branch is ready to merge.

## Locked architecture decisions

These decisions were discussed and approved. Do not reopen or silently bypass
them while resolving the remaining work.

1. There is one public CFF implementation, owned by the
   `three-dimensional-reps` crate. GammaLoop owns graph/source construction,
   UV orchestration, surface mapping, and evaluator integration around it.

2. Numerator energy support is indexed by physical EMR/source-edge energy IDs:
   `Vec<(edge_id, degree)>`. It is not support in loop-momentum-basis variables.
   Exact/local sources must translate that physical capacity explicitly.

3. Preserve main's factorized, incremental numerator ownership. Scaling
   analysis may inspect the factorized numerator first, but 3D UV recursion
   must continue to grow it from spinney-owned numerator fragments. Preserve
   parametric energy signs and shifts; do not expand/concretize tensors as a
   shortcut.

4. There are two UV routes. Both build the complete 4D local and integrated
   counterterms. The local subtraction then either:

   - runs the direct-3D UV recursion on CFF expressions, retaining localized
     orientation-parametric subtraction; or
   - performs 4D-local UV first and projects each finished 4D term/root to CFF
     after its required series and limits.

5. The projected 4D-to-CFF route permits only an explicit source-local sum over
   all orientations. It must reject individual orientation selection, patterns,
   and orientation sampling. Production `OrientationID`s are not expected to
   match exact-source IDs. Keep the sum deferred through evaluator-style
   function maps to limit expression growth.

6. Exact projection uses dual identities: occurrence-local exact energy IDs
   for algebra and physical owner IDs for ownership/topology. Neither ID can be
   substituted for the other.

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

## Immediate continuation: A61

The next product blocker is the raised GL00 `Q7^2` comparison between direct
3D and 4D-local-then-CFF. The present projected route fails before CFF
generation because physical bound `(7, 2)` has not been translated into the
exact occurrence namespace.

This is a lost pre-rebase GammaLoop capability:

- `b8bec53446615f64` introduced
  `expanded_4d.rs::energy_degree_bounds_for_expanded_source`;
- the adapter remained unchanged through integrated source `cc8c5979` and the
  historical monolith `1b73243b4`;
- the pre-rebase comparison passed at
  `-2.22459688841192310e-11 + 4.48261675062501648e-28 i`;
- current GL00 has nine physical graph edges and literal exact occurrences
  `Q1:9` through `Q8:16`, so physical `Q7` maps uniquely to occurrence 15;
- `ExactSourceEnergyMapper` already evaluates the factorized Q7 numerator from
  occurrence 15. Only its capacity-bound translation is missing.

The old adapter is evidence, not code to copy wholesale: it copied a physical
degree to every occurrence sharing `source_edge`, conflating ownership with
energy identity and overcounting repeated owners.

The immediate guard is in `crates/gammalooprs/src/cff/generation.rs`; the
existing occurrence/affine mapping owner is `ExactSourceEnergyMapper` in
`crates/gammalooprs/src/graph/three_d_source.rs`.

The narrow proposed repair, which still requires maintainer approval, is:

1. Map a physical EMR bound to exactly one exact occurrence only when that
   occurrence momentum is literally the same `Q(edge)`.
2. Preserve the factorized numerator and its existing mapper.
3. Hard-error on missing, ambiguous, shifted/general-affine, or otherwise
   unproved mappings.
4. Do not invert `physical_energy_edge_index_map`, use loop carriers as energy
   identity, or copy bounds using `source_edge` ownership.
5. Require the existing GL00 battery and historical numerical oracle to pass;
   do not update its expectation to fit a new result.

A possible rule for multiple literal copies of one physical repeated channel
is deliberately separate and unapproved. A full numerator-row API is still
needed for physical bounded energies with no literal exact occurrence,
shifted/UV/general-affine sources, and shrunken non-bridge cases. Keep A52's
typed failures for those cases until that contract is designed.

After A61, run the live direct GL00 historical oracle under A73. The stale bad
snapshot is not evidence and must not be refreshed.

## Other open GammaLoop work

Prioritize these after the A61/parity work unless a failure exposes a more
fundamental dependency.

### Integration resume boundary (A60)

GammaLoop writes `<workspace>/manifest.json` and
`<workspace>/state/integration_state.bin`. `IntegrationState` contains
`StatisticsCounter` directly and in `slot_stats`; adding
`numerical_stability` changed the nested positional bincode layout. Main and
the selected pre-rebase WIP both lack an independent workspace version, so
this is not a replay regression, but one normal resume path can panic while
decoding.

The minimal proposed fix is a required version-1 field in
`IntegrationWorkspaceManifest`, validated before either decode, plus contextual
`Result` errors for serialization/deserialization. Reject missing, old, and
future versions and recommend `integrate --restart`; do not attempt positional
migration. The remaining non-mechanical decision is whether all unversioned
workspaces may be rejected. Repository policy says old on-disk compatibility
is a non-goal, but obtain explicit approval before implementing this boundary.

### Remaining product checks

- A65: reject `numerator_sampling_scale = 0` when generated CFF actually uses
  uniform sampling. Separately decide whether `All` should activate only for
  repeated-channel degree greater than one; an existing test currently blesses
  degree-one activation.
- A66: add fixed serialized-row and fixed-kinematics numerical CFF goldens from
  the external reference. Do not alter existing tests or fixtures without
  approval.
- A67: regenerate committed GammaLoop schemas, update the live Rust example to
  main's Symbolica pin, restore read-only guards, and adapt the useful CFF
  benchmark suite to GammaLoop's standalone top-level `bench` command. Default
  output placement and test-affecting benchmark choices still need approval.
- A72: decide whether to remove the public but no-op `serde` Cargo feature or
  retain it as downstream compatibility vocabulary.
- Archive coverage: amplitude and cross-section export tests pass, but neither
  performs a typed archive reload and evaluation. Add that coverage before
  claiming the combined archive-version change is end-to-end validated.
- Computed UV-forest export: choose between regenerating the production exact
  CFF catalogue and persisting exact map IDs. Do not infer affine shifts from
  coarse stored orientations.
- A69: later restore the validated pre-rebase auxiliary
  `3Drep graph-from-signatures` command independently of LTD. This is a real
  regression, but it is lower priority than GammaLoop product work.
- A68/A71 and other auxiliary display issues are nonblocking unless a
  pre-rebase regression or GammaLoop effect is demonstrated.

## Validation already completed

The following checks passed on the current rebased code. Symbolica-backed tests
must be run serially: concurrent identical binaries have hit the restricted
license/global-initialization boundary and aborted, while the same tests pass
individually.

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
- the GammaLoop amplitude standalone export and cross-section archive/loader
  export tests;
- an independent generated amplitude-loader type check.

The exact reported GammaLoop command passes twice on the current direct-3D
configuration:

```bash
./gammaloop --release --clean-state \
  ./examples/cli/epem_a_ddx/NLO/epem_a_ddx_NLO.toml \
  run generate_diagrams generate_integrands -c "quit -n"
```

It generates GL0 and GL2 with five evaluators each. The card has
`generate_integrated = false` and selects direct-3D local subtraction, so this
does **not** validate the raised 4D-local-to-CFF route blocked by A61.

The feature-gated shared-core diagnostic inventory is 53/59. Its six failures
are classified under A59 in the detailed ledger. That count is not a GammaLoop
product gate; change an assertion only after establishing regression and test
intent, and ask first as required by repository policy.

## Safe continuation sequence

1. Confirm `jj log -r 'conflicts() & ancestors(raised_energy_cff_wip)'` is empty
   and that `raised_energy_cff_wip` is based on `main@395610143576`.
2. Read A52, A61, A73, and the broad theory audit in the detailed ledger.
3. Ask for approval of the narrow unique-literal A61 mapping before editing
   production code. Do not silently expand it to repeated-owner or affine cases.
4. Make that repair in its own described Jujutsu change and preserve the
   factorized numerator/dual-ID boundary.
5. Run the unchanged GL00 projected battery serially:

   ```bash
   nix develop -c cargo test --locked -p gammaloop-integration-tests \
     --test test_runs \
     scalar_3l_cross_section_gl00_quadratic_energy_inspects_match::q7_squared \
     -- --test-threads=1
   ```

   Also obtain the direct-only historical numerical oracle. Do not refresh
   snapshots to make a mismatch pass.
6. Rerun the exact `epem_a_ddx` command and the focused validation above.
7. Record every architectural finding, approval, test, and remaining failure in
   `kysvnqlq-rebase-review.md`, then refresh this operational summary if its
   status changes.
8. Push only the intended bookmark after inspection. Do not include the
   unrelated conflicted bookmark.

## Things not to infer

- “Conflict-free” does not mean “ready to merge.”
- A passing direct-3D GammaLoop command does not prove projected 4D raised
  support.
- `source_edge`, an occurrence ID, an orientation ID, a loop carrier, and a
  physical EMR energy ID are not interchangeable.
- A spatial UV profile after generation cannot prove the loop-energy contour
  condition that permits CFF generation.
- The auxiliary `3Drep` CLI is not the main product, but a demonstrated loss of
  its validated pre-rebase behavior should still be tracked.
- Proper LTD is preserved for later work; its deferral is not permission to
  delete the setting, source reserve, tests, or architectural intent.

## Stale source warnings

Do not resolve, merge, or branch from the conflicted local
`ltd_in_gammaloop_symbolica_update` line (`swnpzvnu` / `f46f9f56`), its older
remote source (`cc8c5979`), or the remote historical monolith
`ltd_in_gammaloop_symbolica_update_monolith` (`uksvkwzw` / `1b73243b`). The
detached old ledger change `zzxptpwq` is historical as well. Relevant intent is
already captured in the live ledger or the pushed `proper_ltd_core` reserve.
