# `kysvnqlq` rebase architecture and review log

For a concise cross-machine continuation guide, see
[`raised-energy-cff-handoff.md`](raised-energy-cff-handoff.md). This file remains
the authoritative chronological decision and evidence ledger.

## Purpose

This file records the rebase of the change DAG rooted at `kysvnqlq` onto
`main`, including the architectural intent preserved from each change,
conflicts encountered, how mechanical conflicts were resolved, questions that
require a semantic decision, and validation of the rebased result.

The log was initially maintained in a separate Jujutsu change (`zzxptpwq`) so
conflict resolution could not silently alter the intent or ownership of feature
changes. At the maintainer's request it is now included directly in the
`raised_energy_cff_wip` change (`lktlsypu`) and serves as that change's live
readiness record. The detached revision is only the old source copy and is not
part of the intended final DAG.

## Current readiness

Last updated 2026-09-01. The original published feature-code parent was change
`lktlsypu`; the remote `raised_energy_cff_wip` checkpoint is now commit
`5c3e7f`, and the current local non-LTD continuation includes the
A60/A61/A65/A67/A70/A72 maintenance, typed
amplitude archive reload coverage, stored-root selector preservation, mapped
numerator scalarization regression coverage, exact source-local EMR carriers
through integrated finite-UV localization, and the typed production CFF global
prefactor bridge, source-backed local-4D projection repairs, and the current
lower-sector CFF construction. The direct-local3D architecture is now restored:
both direct orientation representations apply Taylor operators only after the
complete/global CFF energy integration and differ only by selector omission.
Exact-source reconstruction and minimax dispatch remain exclusive to projected
local4D. All four direct/converted DD/TT physical acceptances pass in the
retained 10x campaign. The 2026-08-31 scalar LU 15-case all-enabled rerun passed
`15/15`, with `235` skipped, in `70.438 s` under `dev-optim` /
`test_gammaloop`; that remains historical evidence. A current matrix rerun after
the direct-route restoration is still required for merge readiness.
This section is the live status summary; the later sections retain the
chronological decision and validation evidence.

Across this summary, every CFF power and capacity question is posed solely in
physical EMR/source-edge energies; LMB coordinates are never consulted as an
energy-power identity or fallback. The production numerator stays factorized
and grows incrementally by ownership. Higher-power interpolation may substitute
an EMR energy by signed integer `a*M`, where the common scale `M` must be nonzero
only when a finalized evaluator uses it, and results are invariant under
changing its nonzero value.

In the live table below, the recorded 2026-08-31 `15/15` scalar-matrix result is
retained historical evidence. It does not certify the just-restored post-CFF
direct-local3D route; a current rerun remains a merge-readiness gate.

| Area | State | Evidence and remaining work |
| --- | --- | --- |
| Main and file-level rebase | **Complete and conflict-free on latest main** | The selected `raised_energy_cff_wip` ancestry has been replayed from its old `dea00567` base onto local/remote `main` at `39561014`. The obsolete selected `kysvnqlq` dependency-pin change was abandoned because its only remaining effect was to replace main's newer Symbolica stack, contrary to the approved main-owned dependency policy. The evaluator/archive overlap with main #95 is resolved under A76 and the validated resolution is squashed into `lktlsypu`; there are no conflicts in the selected ancestry. The unrelated conflicted bookmark `push-wwvzyzvysvro` remains outside this ancestry. |
| Raised-energy CFF core | **Repeated-channel, finite-pole EMR, factorized rank, source-backed dotted lifting, normalization, and lower-sector frame repairs implemented** | Physical edge-indexed bounds, source/local remapping, factorized analysis, and explicit `None` versus `Some([])` are retained. A55 routes the full occurrence graph through ordinary CFF or the restored channel-normal-form bounded builder. A56 restores independent finite-pole EMR sampling without rewriting diagnostic loop coordinates. A78 removes the bespoke contour certificate, retains the required physical-EMR rank envelope, and adds a source-integral hard gate from the existing all-cycle UV DOD with `DOD_E = DOD_4D - 3 L`. A53 excludes structural non-cut bridges from the CFF denominator domain and leaves their existing main-style 3D projection/mapping ownership intact. The corrected A79 path carries each original `EdgeIndex` through Taylor expansion so cograph/UV source minors can seed and validate the known skeleton, retain cut and physical-map provenance, and bookkeep repeated copies of one wrapper. The rewritten rational denominator channels and their multiplicities then fix the owner-independent exact residue topology: a compatible owner relabeling cannot change its incidence, loop rank, or residue. One sign-safe minimax plan supplies exact bounds and factorized substitutions. A quotient-space solve chooses only a routing sign on the source-backed scaffold. Graphica canonically relabels this already-built topology for caching; no Kirchhoff/incidence reconstruction or balance synthesis remains. Two-pass topology batching joins repeated-channel **totals** while retaining every genuinely outer term's plan. The real DOD1 triangle instead remains one uncancelled common-denominator source with owner multiplicities `(2,1,2)`; generalized CFF supplies its lower sectors internally. A63's complete contour normalization remains valid, while its earlier same-owner/distinct-owner incidence distinction is superseded by this owner-invariant residue contract. A58 restores the validated full-rank top-level pure-CFF boundary while retaining contextual rank-deficient component/particular lifts only inside lower-sector construction. The current full shared crate is `98/100`: the two failures are stale private-selector/fixture contracts, not a demonstrated GammaLoop value regression. The four physical NLO acceptances pass; the final scalar LU 15-case all-enabled rerun passed `15/15`, with `235` skipped, in `70.438 s` under `dev-optim` / `test_gammaloop`. The feature-gated diagnostic evaluator remains comparison evidence only; CLI-only failures are not a product gate unless they expose shared CFF mathematics used by GammaLoop. A64's numerator-call deduplication idea is retained as a deliberate post-rebase follow-up and is not a readiness blocker. |
| Production CFF carrier and prefactor bridges | **Implemented with focused and inclusive coverage** | Integrated finite UV terms retain their exact source-local EMR maps while each production selector identifies one complete generalized residue-map entry. Physical edge directions are metadata within that entry and never define the selector partition. The shared core carries its connected-loop and pure duplicate-denominator global sign as typed metadata. GammaLoop consumes the bridge exactly once for root, reduced, and exact production sources, cancelling that core-local uniform convention and restoring latest-main's complete-integrand convention. The DDx and TTX three-route inclusive acceptances now exercise these bridges end to end. |
| `epem_a_ddx` NLO generation | **Exact reported command reaches integrated direct-3D generation** | `./gammaloop --release --clean-state ./examples/cli/epem_a_ddx/NLO/epem_a_ddx_NLO.toml run generate_diagrams generate_integrands -c "quit -n"` generates and compiles GL0 and GL2 with five evaluators each. The live card has `generate_integrated = true`, so it exercises direct-3D local subtraction and the integrated finite UV addback with exact source-local EMR carriers. The default `local_uv_cts_from_expanded_4d_integrands = false` means it does not validate the alternate 4D-local-to-CFF route. |
| `epem_a_ddx` NLO normalization | **10x physical acceptance passes** | The converted DD run gives `LO=0.1950499744 +/- 0.0010326753 pb` against `0.1935223910 pb` (`+1.479 sigma`) and `NLO=0.007824189766 +/- 0.000340601513 pb` against `0.007268847158 pb` (`+1.630 sigma`). Its `alpha_s/pi` pull, including LO uncertainty, is `+1.453 sigma`. The signed graph-MC components are `GL0=-0.01996339254 +/- 0.00015479207` and `GL2=+0.02786745983 +/- 0.00028425324`; they are closure components with no separate published targets. The route, scale-law, and EMR-bound checks complete rather than stopping on the former stale occurrence-shape assertion. |
| Direct `gamma* -> d d~` MSbar NLO | **10x published absolute and graph values pass** | With the summed Feynman-gauge projector `-g^{mu nu}` and no picobarn conversion, `LO=0.5068703962 +/- 0.0025987972` has pull `+1.449 sigma`, while `NLO=0.01966009810 +/- 0.00053595339` has pull `+1.424 sigma`. `GL0=-0.03132123586 +/- 0.00023922726` and `GL2=+0.05112479005 +/- 0.00046213299` have pulls `+0.729` and `+1.584 sigma`; the `alpha_s/pi` pull is `+1.141 sigma`. The acceptance exercises all three local-UV routes, including precise graphwise comparisons. |
| `epem_a_ttx` fully MSbar NLO normalization | **10x converted acceptance passes** | At `600 GeV`, the run gives `LO=0.3307052414 +/- 0.0018004843 pb` against `0.3278184406 pb` (`+1.603 sigma`) and `NLO=0.02356890542 +/- 0.00056839205 pb` against `0.02296767591 pb` (`+1.058 sigma`). The paper-ratio pull, including LO uncertainty, is `+0.685 sigma`. This summed-graph integration persists no separate GL0/GL2 rows. |
| Direct `gamma* -> t t~` fully MSbar NLO | **10x published absolute and graph values pass** | The direct-current run gives `LO=2.901968994 +/- 0.015639978` (`+1.641 sigma`) and graph-summed `NLO=0.2079169992 +/- 0.0042953541` (`+1.489 sigma`). `GL0=-0.1443600613 +/- 0.0035809931` and `GL2=+0.3522770605 +/- 0.0023720361` have pulls `+0.669` and `+1.687 sigma`; the paper-ratio pull is `+1.037 sigma`. This is the published direct-current setup with the summed `-g^{mu nu}` projector and no conversion. |
| UV profiling | **Amplitude/LU divergent-limit defaults and restored-route coverage implemented; rerun pending** | `--selected-limits only-divergent` profiles every expected cycle union with DOD >= 0. It uses the generation LMB when suitable, otherwise the first suitable basis in the deterministically sorted complete LMB list; `all` remains opt-in. Graph and Cutkosky-cut selectors work for LU input, and failures are summarized in a colored table. Per-orientation profiling is valid only for orientation-parametric localized direct local3D; explicit-sum direct local3D and projected local4D are summed representations and reject it. The selected focused regression covers DOD0/1/2 bubbles, and every selected base scalar-matrix graph (GL00, GL02, GL04, GL08, GL09, and GL24) now invokes the localized-direct3D per-orientation profile. Their post-restoration rerun is pending. Current end-to-end CLI coverage is `2/2` and lower-level API/unit coverage is `15/15`. |
| Direct-3D versus 4D-local-then-CFF parity | **Architecture restored; current parity rerun pending** | Both direct forms apply every Taylor operator after complete/global CFF energy integration; their local bodies are identical and explicit-sum only omits selectors. A79's exact reconstruction and derivative-created minimax dispatch apply exclusively to projected local4D. A77 keeps the empty forest as the ordinary factorized production root in both routes. All four physical DD/TT results and the 2026-08-31 scalar `15/15` run remain historical validation; the current restored-route scalar matrix has not yet been rerun. |
| Direct-3D historical parity | **Protected GL11 oracle reclassified; snapshots unchanged** | Earlier audit text mislabeled the protected graph as GL00. The stored GL11 `Q7^2` value came from a validated path that explicitly cleared numerator-energy bounds before root generation, so its cancellation residue is not a generalized raised-EMR value oracle. A73 records the evidence and the approval boundary. |
| Scalar snapshot audit | **Historical 15-case pass retained; restored-route rerun pending** | Historical proper-LTD-vs-CFF runs validate the live GL00 `Q7^2` magnitude `1.317240193281177e-9`; the negative-imaginary `.snap.new` remains stale because a fresh rebuild gives the positive phase. GL02 without a numerator is Arb-stable and equals the historical LTD value up to the current global `-i` convention at cancellation magnitude `4.2e-16`. The earlier GL24 investigation used a nonlocal global numerator and a handwritten contour calculation, so its mismatch claim was not authoritative. The replacement scalar LU 15-case matrix starts from generated scalar graphs, retains unmodified unit-numerator lanes, and uses companion Feynman-rule-local edge numerators. It compares all three routes at f64/native Arb with local UV, integrated UV, and threshold counterterms enabled. Its 2026-08-31 pre-reversal rerun passed `15/15`, with `235` skipped, in `70.438 s` under `dev-optim` / `test_gammaloop`; a current rerun is required after restoring post-CFF direct local3D. Four near-zero cases use the authorized f64-input `1e-14` unit-scale fallback, while Arb-to-Arb remains run and reports non-scaling; this is test-oracle handling only and required no production change. No graph-specific production repair or graph-name branch is present. The protected GL11 expectations remain a separate test-intent issue. |
| GammaLoop bridge ownership and pre-rebase `3Drep` residual compatibility | **Implemented and focused-validation complete** | GammaLoop continues to contract structural non-cut bridges and retain their factors in caller-owned `tree_denoms`. Separately, the auxiliary `3Drep build` command preserves its validated pre-rebase behavior: it selects paired, non-dummy `tree_edges - initial_state_cut`; the shared core contracts those requested zero-loop edges from CFF, lifts the active result back to physical IDs, and returns their affine maps plus typed `ResidualDenominator`s. Pure trees return one identity orientation plus residuals. This is retained as a regression oracle, not as the GammaLoop production boundary. |
| Saved-state compatibility | **Top-level and integration-workspace boundaries implemented** | Top-level saved state writes manifest version 4 and rejects incompatible versions before decode or overwrite. Version 4 fences the positional generated-CFF layout that persists the typed core global-prefactor sign. Integration workspaces independently write version 1 and reject missing, old, or future versions before either resume decoder, with contextual `integrate --restart` guidance and no positional-bincode migration. |
| GammaLoop QoL and committed interfaces | **Mechanical maintenance complete; larger test/tooling choices remain deferred** | Committed schemas match the live CLI, the Rust example uses main's Symbolica source, explicit read-only output guards are restored, the Python API again compiles against typed stored 3D expressions, and BNL documentation uses the top-level `bench` vocabulary. The four integration-helper dead-code suppressions remain necessary because their only caller is still inside a legacy commented method. Default-output policy, the late-parsing benchmark fixture, and restoration of the broader benchmark suite remain unchanged because they require test or tooling-policy decisions. Auxiliary `3Drep` tooling remains outside GammaLoop readiness. |
| Computed UV-forest export | **Deferred; not ready** | Persisted graph terms contain only coarse orientations and cannot reconstruct exact affine map identity. Production generation/evaluation uses exact maps and is unaffected, but computed-forest export must either regenerate the production catalogue or persist exact map IDs before it can be declared exact. |
| Deferred proper LTD | **Preserved separately** | Live `RepresentationMode::Ltd` selection returns a typed `NotImplemented` error. `proper_ltd_core` / `377dd390` is pushed as source preservation, is not an ancestor of the WIP chain, and must later be ported onto the finalized CFF + 4D-UV head rather than replayed wholesale. |
| Remaining split material | **Preserved, not yet landed** | The NNLO example-card delta and the separately audited QoL candidates remain source material. Proper-LTD comparison/evaluation tests remain assigned to the later LTD delta. The generic CFF residual-denominator contract and its CFF-only tests are now part of this WIP under A49; no LTD runtime was restored with them. None of the remaining material is silently rejected or part of the current WIP tree. |
| Documentation, cleanup, and push | **Remote checkpoint refreshed; final local publication pending** | This ledger remains the authoritative chronological record and `raised-energy-cff-handoff.md` is its concise companion. The published `raised_energy_cff_wip` bookmark now resolves to commit `5c3e7f`; newer local lower-sector validation and test-contract cleanup must be compared against that checkpoint and pushed only after final cleanup, fetch/rebase, validation, and artifact audit. `proper_ltd_core` remains a pushed sibling reserve. |

## Rebase contract

- Source root before rebase: `kysvnqlq` / `5c2d0c19` (`update workspace dependency pins`).
- Original parent: `lwolvxrp` / `4fdbf430` (`Validated Muv (#88)`).
- Destination at start: `main` / `21cefaee` (`Disconnected UV (#93)`).
- Current selected base: `main` / `39561014` (`Improve Typst graph drawing and physics edge styles (#73)`).
- Approved scope: the source and the split/reviewable descendants. The alternate
  monolithic change and bookmark are retired rather than rebased.
- GammaLoop is the product boundary and supplies the terminology used for
  readiness. The shared `three-dimensional-reps` library remains in scope where
  GammaLoop calls it.
- The separate `3Drep` CLI and its feature-gated evaluator are diagnostic only,
  not production contracts. Their input normalization, expression preparation,
  and display need not match GammaLoop. CLI-only regressions and enhancements do
  not affect readiness unless they expose shared CFF mathematics used by
  GammaLoop.
- Mechanical resolutions may be applied directly and must be documented below.
- Any resolution that changes behavior, public or internal architecture, dependency
  selection, test intent, or the relative meaning of both sides is non-mechanical
  and requires a maintainer decision before it is applied.
- Existing comments move with their code and are not deleted during conflict
  resolution.

## Original topology and intent

The DAG has two heads. `ltd_in_gammaloop_symbolica_update` is the integrated
split implementation; `ltd_in_gammaloop_symbolica_update_monolith` is an
alternate monolithic implementation branching after the shared split
infrastructure.

The two heads have byte-identical trees before rebase, but different histories.
The maintainer does not need the monolithic alternative, so only the
split/reviewable history is retained. The final split merge `swnpzvnu` is not
empty despite its description: it contains a small
`gammaloop-api/src/state.rs` reconciliation. That reconciliation now belongs to
`tmlszmts`; the old structural merge is superseded by the fresh narrow deltas.

The required feature order after porting is architectural, not merely a testing
priority:

1. raised-energy support for CFF;
2. four-dimensional UV support on top of that CFF foundation;
3. proper LTD support only afterward.

Conflict resolutions must not make the initial CFF/4D-UV layers depend on the
later proper-LTD route.

```text
kysvnqlq  update workspace dependency pins
└─ xvxsluoq  support split infrastructure
   ├─ uksvkwzw  ltd in gammaloop symbolica update [monolith head]
   └─ yqmtxtlv  update tensor and symbolic support crates
      └─ mwpprnsu  add CLI state output helpers
         ├─ sttlspyn → wrzsxkku → tmlszmts → qontqzwm  [QoL line]
         └─ yqnvuymp → kqqlnrnx → nrzomklu           [3D-representation line]
            ╲                         ╱
             lrznwzuo  merge QoL and 3Drep branches
             ├─ tormmyww  update example run cards
             ├─ wskyzxur  add repeated-mass CFF regression tests
             └─ ltnoyvll  add LU scalar cross-section regression battery
                ╲              │              ╱
                 swnpzvnu  integrate split feature branches [integrated head]
```

The per-change review below will be expanded from commit messages and diffs
before conflicts are resolved.

## Conflict ledger

| Change | Files / subsystem | Classification | Resolution or pending question | Review evidence |
| --- | --- | --- | --- | --- |
| _Preflight_ | Rebase scope | Non-mechanical, decided | Retire the monolithic change/bookmark and rebase only the split/reviewable history. | `kysvnqlq::` contained the split head `swnpzvnu` and the byte-identical alternate 419-file monolithic head `uksvkwzw`. The maintainer explicitly does not need the monolithic history. |
| `kysvnqlq` | Symbolica dependency policy | Non-mechanical, decided and refreshed on latest main | Use main's Symbolica generation instead of restoring any historical branch pin. | The original target used Symbolica 2.1 / `9135b409`; the first rebase resolution used main's then-current `46d690719`. Main #95 subsequently advanced the resolved Symbolica 2.2 stack to `0441bd7a`, enabled `native_code_generation`, and added O2 compilation forwarding/defaults. The latest replay keeps all of that main-owned behavior. |
| `kysvnqlq` | Dependency files | Superseded by latest main; selected change abandoned | Drop the selected rebased `kysvnqlq` instance because its only live diff was the now-stale `46d690719` pin. Keep main's `Cargo.toml`, `Cargo.lock`, `crate-hashes.json`, and Hakari manifest rather than manufacturing generated-file merges for a downgrade. | Abandoning that selected divergent instance removed all three dependency/generated-file conflicts and replayed its fourteen descendants directly onto `39561014`. The historical source instance remains available outside the selected live stack. Offline locked metadata resolves, and the Hakari manifest/lockfile were mechanically regenerated for the retained shared crate and pass `cargo hakari verify`. |
| `xvxsluoq` | Split infrastructure scope | Non-mechanical, decided | Prune to feature-essential infrastructure: drop unrelated README, global timeout/stack, test-recipe, and macOS linker rewrites unless later evidence requires them; retain feature output ignore rules; move 3D-crate registration to the crate-introduction change. | The cleanly applied commit changes nextest timeouts, ignore rules, README badges, 15 build scripts, Nix/test stack size, and justfile behavior. Registering `three-dimensional-reps` before that crate exists makes intermediate clippy history invalid. |
| `yqmtxtlv` | `EvaluationDomain for Complex<T>` ownership | Non-mechanical, decided | Discard the duplicate implementation added to `complex.rs` and the already-present `map_coeff` bound, then remove only the unnecessary `+ FloatLike` bound from main's owning implementation in `complex/symbolica_traits.rs`. | Main already moved the impl and already provides the branch's default `resolve_function` behavior. Applying the patch verbatim creates overlapping impls and imports Symbolica `FloatLike` explicitly, which repository guidance forbids. |
| `yqnvuymp` | Core-crate feature boundaries | Non-mechanical, resolved | The initial core stage is CFF-only and locked-buildable. Its inert one-variant `RepresentationMode::Cff` scaffold was removed there; A28 later introduces a real `Cff`/`Ltd` setting while keeping LTD generation explicitly unavailable. | An isolated archive passed formatting, Hakari, locked all-feature check, Clippy with warnings denied, and all 42 original tests. The integrated CFF stage also passes display/eval package checks. No LTD residue, evaluator, or comparison code is present; A42 retains one public CFF owner and A55 restores its validated repeated-channel organization. |
| `xvxsluoq` | Approved infrastructure pruning | Non-mechanical, resolved | Reduced the change to five ignore rules for nested integration workspaces and local benchmark outputs. | The final tree relative to its parent changes only `.gitignore`; all broad build/test/linker/README policy churn and premature crate registration were removed. |
| `yqmtxtlv` | Approved main-owned trait port | Non-mechanical, resolved | Reduced the change to removal of the unnecessary `+ FloatLike` bound in main's existing `complex/symbolica_traits.rs` implementation. | The duplicate implementation and redundant `map_coeff` edit were removed. `cargo check --offline -p spenso` passes. |
| `wrzsxkku` | Symbolica-state import ordering | Non-mechanical, decided | Retain main's `State::new` → model/input-card → Symbolica-state order and discard the branch's extra tracing override. Add fresh-process reload coverage for the formerly false warning and preserved UFO custom printing. | The branch reorder originated in `4afcaa3` as a Symbolica 1.5 workaround for a false duplicate-UFO warning and explicitly accepted plain formatting after reload. Symbolica `24fe7dad` fixed that exact warning before both the old 2.1 pin and the approved 2.2 pin. Current import warns only when it must create a genuinely new non-exportable symbol. Model-first registers the UFO symbol's Spenso/Typst callback and import remaps the saved ID; import-first creates a plain symbol that later model loading cannot upgrade. Main's `State::new` also deliberately installs GammaLoop tracing before legitimate import diagnostics. The maintainer approved main's order. |
| `wrzsxkku` | Mixed change ownership | Non-mechanical, decided | Keep only parameterized-command-block behavior here; move `GraphImportOptions`, inline-DOT/process-spec parsing and coverage to `tmlszmts`; move 3D command completion coverage to the later CFF/LTD CLI stages; preserve unrelated main behavior. | The original split commit removes main's Bench, changes preprocess arguments, introduces an import API before its callers, contains tests for later `--inline-dot`/`--process-spec`, and asserts nonexistent `3Drep`/`test-cff-ltd` commands. The maintainer approved repairing topic ownership and intermediate buildability. |
| `wrzsxkku` | Exact topic partition | Non-mechanical, resolved | Kept all `display.rs`, `run.rs`, and `session.rs` changes; kept only placeholder/run-`-D` behavior in `repl.rs`; kept only template serialization/discovery/filtering and run constructors in `state.rs`; kept the six run-variable CLI tests and late-parsing template test. Moved every triple-quoted inline-DOT parser/API/test hunk to `tmlszmts`; deferred 3D/LTD completion assertions; removed state-access and exhaustive-settings coverage from this topic. | The final cumulative patch contains exactly the six approved files. Main's model-first load order, top-level GammaLoop `bench` command, and generation preprocess arguments are unchanged. Formatting, API check, all three completion tests, the late-parsing test, and the stateful CLI workflow pass. |
| `tmlszmts` | Inline-DOT and process-definition import boundary | Non-mechanical, resolved | Own triple-quoted command parsing, `GraphImportOptions`, process-definition validation, inline-DOT/process-spec CLI flow, every API caller migration, inline-import completion, and focused cut-filtering coverage here. Pulled the Python and second state-test caller migrations backward so this intermediate commit builds independently. Preserved main's path/reference disambiguation and comments. | The final cumulative patch excludes the original split's unrelated `renormalize.rs`, `save.rs`, macOS build, read-only, benchmark-removal, and 3D/LTD changes. A final parser review caught and fixed a missing index advance after an escaped character in a double-quoted command block, with focused regression coverage. |
| `tmlszmts` | Physical-cut filtering test ownership | Mechanical history repair, resolved | Keep the original inline process-spec fixture, full 25-to-9 cut-filter regression, and count metadata with the feature they test. Extend main's existing `compute_process_valid_cuts` helper rather than restoring the branch's duplicate selection helper. Drop their later duplicate additions from `kqqlnrnx` and `ltnoyvll`. | The aggregate source introduced the fixture, import path, test, and count API together; the split accidentally scattered them. They use shared cross-section/CFF types and do not require 3D, UV, or LTD generation. |
| `qontqzwm` | Inspect bench versus detailed output | Non-mechanical, implemented | Keep main's `show_detailed_output` and generated-event policy for normal single evaluation. Move the branch's richer fixed-point benchmark engine to the top-level GammaLoop `bench` command and remove `inspect --bench`. | The GammaLoop bench engine is on the clean main-based line; all six focused benchmark tests pass. |
| `qontqzwm` | Inspect JSON under read-only state | Non-mechanical, decided | Keep JSON output, validate at the `Commands::Inspect` dispatch through existing `ensure_write_target_outside_active_state`, and add focused regression coverage. | Read-only session policy is owned by `CLISettings` at command dispatch. This guards both normal and bench JSON paths before integrand lookup without changing lower-level Rust inspection APIs or adding a helper. |
| `qontqzwm` | Stale main regressions and deferred LTD coverage | Non-mechanical, resolved under A5/A7 | Preserve main's Python CFF fields and `generate_cff*` APIs, approved Symbolica pin, profile/event-order/mass-approach tests, and typed-UV behavior. Split branch-local LTD profiles/routes/tests into the later LTD-owned change rather than placing them in the QoL parent. | The clean QoL line retains main's APIs/tests and intentionally contains no LTD implementation. Every LTD source hunk remains tracked for the deferred LTD delta; none is rejected merely because the current CFF/4D stages do not implement it. |
| `qontqzwm` | Bench command ownership | Non-mechanical, decided | Keep one top-level GammaLoop `bench` command, replace its legacy engine with the branch's fixed-point benchmark engine, and remove `inspect --bench` rather than retaining two benchmark boundaries. | The branch path is a benchmark subsystem: warm-up, duration-based sizing, batches and uncertainty, stage timing, progress, fixed-point selectors, minimal-integrand settings, JSON, and propagated errors. Keeping it inside `Inspect` adds a second action and makes `Inspect::run` return an arbitrary warm-up value. Main's top-level GammaLoop command is the clearer CLI boundary, but its implementation must be replaced: it ignores `n_cores` and evaluation errors, permits zero samples, and has no warm-up, statistics, selectors, or structured output. A future random-throughput mode can live under `bench` if needed. The maintainer approved the GammaLoop `bench` boundary. |
| `xywpoyxv` | QoL resolution scope | Non-mechanical, implemented | Keep all audited QoL behavior: main's detailed inspect/event behavior, one fixed-point GammaLoop `bench` command, inspect JSON with the read-only guard, stability reporting, nonzero CLI status for runtime failures with subprocess coverage, and the four Python read-only-state/access-folder getters. Relocate the A13 fresh-process state-load regression to its owning compatibility change rather than dropping it. | The QoL changes are folded into clean `qontqzwm` on current main. The A13 regression and load-order comment now belong to `wrzsxkku`; its fresh-process test passes. |
| `kqqlnrnx` | Main-owned model/state fixtures | Mechanical under A3/A13, resolved after final rebase | Preserve main's scalar-model parameter order, retain required `goldstone` schema fields, and regenerate `state_map.bin` against the final approved Symbolica/model tree. | A clean-process export produced the 22,753-byte SHA-256 `dca35a6485e22e528fc15b627ae2af496fa7f7a410a81e9fa21f7b1f2ac521be`. The pre-existing generator test then fails after export because it imports from its write-only `File::create` handle (`EBADF`); no test was changed. |
| Remaining QoL audit | Additional branch behavior | Non-mechanical, pending maintainer confirmation before application | Review as separate changes: constant/plateau fitting, numeric UFO `complexconjugate` evaluators, and the read-only diagram-generation guard. Keep main's explicit error for missing threshold kinematics. Retain the independent exact-zero/factorial/factor-two/comment/test precision corrections mechanically where they still apply. | These behaviors are outside the raised-CFF stage and do not leak into `smlxpnry`. The Symbolica-dependent derivative-order and rescaling-symbol changes are not QoL and must be re-proved against main's approved Symbolica revision before reuse. |
| `kqqlnrnx` | Raised-energy CFF integration | Non-mechanical, implemented under A14–A24 and revised by A42–A68 | Keep main's `cff::orientations::GraphOrientation`; implement that owning trait for external 3D-representation orientation data; retain all physical cuts and `CutGroup` terminology; port richer E-surface-to-cut mappings without restoring `RaisedCut*`, the left-stride bug, or the old UV bypass. | Main's representative/internal-only orientation selectors and UV orchestration remain authoritative. The port keeps exact affine numerator maps, spinney ownership, source-local edge-indexed support, and map-indexed UV localization through one public CFF owner. A55 replaces the intermediate mass-jet organization; A56–A66 record the remaining audited semantic and validation defects. |
| `kqqlnrnx` | Raised-CFF numerator ownership in 3D UV | Non-mechanical, maintainer-confirmed | Preserve main's spinney-local numerator growth exactly. The complete prepared numerator may be inspected before generation to determine the energy-degree envelope, but a UV transition attaches only `graph.numerator(&current.reduced_subgraph(given), given.subgraph())`; the outside numerator and global factors remain delayed until final assembly. Never restore the branch's full-production-numerator root override. | Main's `local_3d::start`/`t_tilde`, `Numerator::fill_in_reduced`, and `FinalIntegrandBuilder::build_3d` already partition numerator factors by spinney ownership. Historical commit `35b33f1f3`, retained by `626cd990e` and `b8bec5344`, inserted an expanded full numerator only in the root final output while recursive children still used denominator-only CFF, so it cannot implement raised 3D UV consistently. |
| `kqqlnrnx` | Factorized numerator and parametric energy maps | Non-mechanical, maintainer-confirmed | Keep numerator factors multiplicative and paired with the exact `loop_energy_map`/`edge_energy_map` through residue selection, UV recursion, disconnected joins, and terminal materialization. Energy signs remain parametric and the complete affine numerator-energy map remains explicit; do not reduce a raised term to only an `EdgeVec<Orientation>` or call the full-numerator `.expand()` path. UV rescaling and localization remain transformations of these mapped factors rather than part of numerator-map identity. | `OrientationExpression` and `LinearEnergyExpr` already retain numerator-map identity, internal OSE coefficients, external-energy shifts, uniform-`M` shifts, and constants. `Graph::cff` currently discards these maps when it creates `CFFTerm`, which is sound only for denominator-only CFF. Multiple shifted numerator maps can share one coarse orientation and must remain distinct while later spinney factors are grown with the same map. |
| `kqqlnrnx` | Exact-map UV localization carrier | Non-mechanical, maintainer-approved; corrected after integrated-NLO diagnosis | Retain one UV forest traversal with each additive local-3D branch identified by both a complete production residue-map-key selector (`OrientationID`) and, where the CFF source owns the algebra, its exact source-local `edge_energy_map`. The production ID selects the full generalized map entry; it is not a coarse physical-direction selector and does not replace or have to reproduce a contracted source map. Carry that source map through finite-CT localization, incremental `start`/`t_tilde` numerator growth, active-sector replay, and final complement/global-numerator assembly, then collapse to ordinary integrands. Stored production-owned branches use their exact typed ID/map; only a reduced branch without its own source map uses strict affine matching, never a coarse fallback. | The integrated `epem_a_ddx` GL0 `Fyy` source contains a legitimate quadratic contact map with `Q3 -> 0`, while every direction-compatible root production map has `Q3 -> +/-E3`; no full affine extension exists. Requiring equality therefore discarded real source capacity. Reusing `LinearEnergyExpr` in a private `(selector ID, optional source map)` carrier preserves that capacity without a second numerator abstraction or another forest traversal. Distinct source maps sharing one selector remain separate additive branches, while identical pairs coalesce; no Cartesian product is formed. |
| `smlxpnry` | Exact-map implementation review | Non-mechanical under A17/A22, corrected for source-local integrated maps | Keep typed production root maps and contracted source maps through `Graph::cff`, defer full orientation-pattern filtering until exact projection, and retain local-3D branches by the dual selector/source-map identity. Union distinct source maps even when their selector IDs agree; zip-add only an identical pair. Map each newly owned numerator fragment and integrated finite CT exactly once with that branch's source map, leave root selector extensions unnormalized, and apply `1/N` when one collapsed local or integrated source residue is partitioned over `N` compatible production selectors. Disable the final coarse sign fallback in exact mode. Restore an omitted contracted-parent map only from exact zero-inner source coordinates. | Independent review first corrected root `1/N` scaling, premature reduced-pattern filtering/no-match conflation, unmapped finite CTs, and zero-filled reconstructible parent edges. The later integrated-NLO diagnosis showed that affine equality to a production map is not a valid requirement for a source-owned contact map. Production IDs now supply only compatible theta sectors for such a branch, while its exact source EMR map remains authoritative for all factorized numerator owners. |
| `smlxpnry` | Contracted raised-energy capacity | A18 representation approved; coordinate semantics repaired under A50 and audited under A51 | Keep factorized `Atom` analysis, but express its conservative degree envelope in original EMR/source-edge IDs and remap those IDs through every contraction/component embedding. A declared loop symbol is accepted only as an alias for its explicit LMB carrier edge, so the stored capacity still belongs to that EMR edge. | The historical edit has restored the edge-indexed adapter and strict source-edge mapping. The theory audit confirms the ownership model, while final validation remains blocked by the lower-sector/repeated-channel issues in A51. |
| `smlxpnry` | Sparse raised-CFF reconstruction | Historical implementation checkpoint; carrier resampling invalidated by A56 | Keep correlated support global and use sparse interpolation/contact weights rather than Cartesian boxes. | Sparse Vandermonde weights, recursive remainder/contact completion, and disconnected denominator products were implemented. The later choice to alter an LMB carrier and recompute every dependent EMR map is not part of the approved edge-owned semantics: finite-pole division must replace only the selected `edge_q0`. Historical origin assertions are catalogued in A59. |
| `smlxpnry` | Uniform-`M` raised CFF | Historical implementation checkpoint, superseded by A55/A65 | Sample each nonlinear repeated channel at signed `q_e=a_e M`, retain sparse support correlations, and leave residual affine variables exact. | Sparse `{0,2,5}` and correlated `{K0 K1,K0^2 K1}` reconstruction established the required uniform/nonuniform equivalence and affine-shift behavior. The one-CFF channel-normal-form builder now owns this mechanism. A65 records the remaining degree-one activation policy and missing GammaLoop-side nonzero validation. |
| `smlxpnry` | Last-denominator polynomial contact | Non-mechanical boundary approved and implemented; valid-topology focused coverage pending | A degree-two last-carrier quotient is represented as a formal Unit-denominator constant contact. Degree three or higher leaves a positive power of a free loop energy, which the current 3D expression does not represent. | Sparse nonuniform/uniform, legacy quadratic, known-factor, and repeated-channel paths retain scalar contacts and return `UnlocalizedEnergyContact` before partially emitting a positive free-energy quotient. The terminal quartic Unit/contact production regression passes, but the earlier dedicated None/All degree-two-through-four test used a one-edge self-loop fixture and is no longer a trustworthy live topology oracle. A balanced valid-topology replacement requires explicit approval. `HybridSurfaceID::Infinite` remains rejected as a fake substitute. |
| `smlxpnry` | Rank-deficient lower sectors | Historical A24 implementation checkpoint; boundary retained by A58/A62 | The rank-projected implementation proved the nonintegral-denominator, null-support, initial-cut, rational-coordinate, and correlated-numerator boundaries before its raw builder was removed. | The old `LowerSectorRankProjection` and repeated-residue construction are not retained as a second production path. A58 restores full active-denominator rank at the public top-level boundary while the current lower-sector component builder retains its contextual particular lift; A62 separately certifies scaling subspaces inside that active span. The replacement cut tests exercise both sides explicitly. |
| `smlxpnry` | Core/GammaLoop migration cleanup | Implemented after zero-caller audit and completed by A44 | Remove genuinely dead duplicate migration APIs while preserving live causal conventions, source coordinates, comments, and all deferred LTD source in its separate snapshot. | GammaLoop cleanup removes an amputated duplicate generation island; A44 removes the remaining raw repeated-residue structs and unused provenance/sign machinery. Proper LTD source is preserved outside the live CFF tree. |
| `smlxpnry` | Generated Nix workspace graph | Mechanical, resolved | Regenerate `nix/ci-workspace-graph.json` after adding the 3D crate dependency rather than accepting the stale generated file. | `nix build .#ci-workspace-graph-json` regenerated the file and `nix build .#checks.x86_64-linux.gammaloop-guppy-workspace-graph --no-link -L` passes. |
| `smlxpnry` | Computed UV-forest export | Non-mechanical, deferred | Production amplitude/cross-section generation uses exact maps, but `uv/export.rs` still receives coarse orientations from persisted process-integrand graph terms. Before declaring computed-forest export exact, decide whether export regenerates the production exact CFF catalogue or the persisted term format stores the exact map identity. | The current persisted graph-term schema has no full `OrientationExpression`/`OrientationID`, so the export-only `OrientationProjection::new` path cannot reconstruct affine shifts from an `EdgeVec<Orientation>`. This does not justify a coarse fallback in production UV. |
| `rpspmzvw` | Deferred 4D/LTD staging | Non-mechanical repair required under A7/A8 | Preserve this sibling full-tree reserve as LTD source material, but do not replay the full snapshot into the current CFF/4D ancestry. Later extract proper LTD as a delta on the post-4D descendant. Review its generic preserved-denominator contract independently as shared CFF material; retain shared comments/tests and leave the unrelated macOS linker workaround outside the feature. | The reserve is a sibling of `yqnvuymp`, re-adds the entire crate, lacks its lockfile entry, combines `ResidualDenominator`/`preserve_internal_edges_as_four_d_denominators` with LTD, and deletes existing tree comments and seven unit tests relative to the CFF core. The residual contract does not belong to proper LTD, is not required by A48's global-catalogue fix, and remains a separate pending port for generic shared-library CFF output. |
| `kqqlnrnx` | E-surface cache ownership | Non-mechanical, maintainer-approved and implemented under A34 | Keep exact topology-inventory IDs separate from expression-local raised-residue IDs. Persist dense `TopologicalThresholdId`s and exact topology surfaces; translate only eligible counterterms to the frozen expression's `RaisedEsurfaceId` by exact-first, raised-edge-normalized-second matching. | Topology classification, persistence, diagnostics, and API edge reporting no longer depend on the expression cache. Numerical counterterms remain expression-owned. This replaces the earlier blanket cache-lookup proposal and intentionally permits dense public IDs to differ from historical expression-cache offsets. |
| `kqqlnrnx` | Separate direct-3D and 4D-local UV routes | Non-mechanical, implemented; focused route tests pass, scalar parity blocked at A45 | Preserve main's direct `Local3DApproximation` recursion as one production route. Add a separately selectable 4D-local route using main's typed `Local4dCts`/`Full4dCts`; replace only the branch's competing `Expanded4DApprox`/`ApproxOp` implementation, not the route or its parity-test intent. | Change `ypnxsozo` restores main's local-3D state, dependency-frontier replay, numerator-growth kernel, symbol builders, and final-integrand localization as the default. Setting `local_uv_cts_from_expanded_4d_integrands=true` instead composes typed 4D local terms and projects them to CFF only afterward. Direct and typed union/spectacles regressions coexist and pass, but they did not cover GL11's distinct-physical/coincident-denominator topology. |
| Post-CFF 4D bridge | Residual-denominator ownership, exact normalization, and orientation ownership | Non-mechanical under A19/A27/A38; focused tests pass, scalar parity blocked at A45 | Keep `Local4dCts`/`Full4dCts` as the sole owner of complete 4D local/root atoms on the projected route. Per additive term, preserve exact external-energy-only `GS.den(edge,momentum,mass,full_expr)` occurrences as residual factors, send loop-dependent occurrences through raised CFF, and recombine each residual exactly once. Multiplicity comes from the actual term, never graph metadata. Powered-pole variants own retained half-edge normalization; an ordinary global-source result receives the global inverse-energy product. Each post-4D CFF explicitly sums its source-local orientations and never matches them to production `OrientationID`s. | The projected carrier stores compact tagged calls and source-local bodies separately. Every body maps the complete factorized numerator under one exact source map before evaluator-local lowering. Focused exact-CFF, projected-local, disconnected, tensor, and deferred-carrier suites pass, but GL11 proves that source-wide parallel topology selection is too coarse. This term-local projected-4D ownership does not provide generic shared-library CFF preservation of caller-unowned tree propagators. Proper LTD remains a separate preserved delta. |
| Post-4D representation setting | Deferred LTD selection | Non-mechanical, implemented under A28 | Expose serializable `Cff`/`Ltd` representation vocabulary and `3Drep build --representation`/`--family`, defaulting to CFF. Reject LTD immediately and explicitly; do not introduce LTD generation, residue metadata, evaluation, or comparison support. | Production CFF option constructors set `RepresentationMode::Cff` explicitly. CLI LTD selection returns the core typed `NotImplemented` error before graph selection or artifact writes, and each public core generation entry point enforces the same boundary. |
| `nrzomklu` | CLI feature boundary | Non-mechanical, implemented under A26/A28 | Create a fresh CFF-only delta with `3Drep validate` and `3Drep build`; use exact graph support/generation, JSON artifacts, and CFF display. An explicit per-build numerator-sampling mode overrides the global setting for that invocation without mutating it. A later narrow change exposes the representation selector but errors on LTD. Defer Evaluate, PureLTD/LTD comparison, manual degree bounds, reconstruction, and mixed demo cards. | The source is a roughly 13k-line mixed change. The fresh 610-line delta uses one private borrowed `GraphCatalog` because no existing selector spans imported/generated amplitude and cross-section graphs. It always regenerates from live graph/numerator state, retains artifacts and `latest`, and has no unsafe historical cache or `--clean` path. Focused end-to-end coverage passes. |
| `tormmyww` | NNLO example card | Non-mechanical, deferred until reached | Use main's card as the base; review branch-only benchmark blocks/assets against the staged architecture instead of taking either 600+ line file wholesale. | The two cards differ in graph selection, evaluator, threshold/UV modes, samples, stability, and integration targets. |
| `wskyzxur` | Repeated-mass CFF regressions | Non-mechanical boundary implemented and validated | Port the two threshold-disabled split/equal mass-limit tests and all four genuinely new graph fixtures as a fresh CFF-only delta. Remove the obsolete `global.3d_representation=CFF` setting, retain the original sign and strict convergence assertions, and ask before changing either if they fail. | Fresh change `xxwzulrk` reuses existing CLI/workspace, `inspect_xspace_process`, and `complex_distance` behavior plus the two approved private test helpers. Both additional-parameter and model-parameter tests pass unchanged. The third unit-localizer/threshold/event/Arb test remains deferred to the later 4D/threshold stage. |
| `ltnoyvll` | Raised-power scalar route-comparison battery | CFF delta implemented; semantic validation blocked; LTD leg preserved and deferred | Preserve its CFF direct-3D versus CFF 4D-local comparison, adapting the latter from `Expanded4DApprox` to main's typed 4D route. Split the LTD comparison and its setup into the later LTD-owned test change; do not delete it. Retain the UV-disabled fallback cases as CFF/threshold coverage while distinguishing them from raised-power UV support. | Fresh change `lktlsypu` ports the 198 live CFF cases. A later source-history audit found 195 stored snapshots byte-identical, the approved GL02 `Q1^4` exact-zero change, an unapproved cancellation-scale GL02 no-numerator phase change, and a material unapproved GL11 `Q7^2` change. No additional snapshot update is approved. |
| Standard saved states | Positional bincode schema version | Non-mechanical, maintainer decision pending | The projected carrier and explicit-sum metadata add fields to saved process/evaluator layouts, while `state_manifest.toml` still reports version 1. Decide between a narrow version-2 boundary with an early regeneration error for v0/v1, or a substantially larger frozen-v1 migration layer. | Standalone amplitude/cross-section archives were independently version-bumped and compile against main's Symbolica. Ordinary v1 state folders currently pass manifest validation before failing opaquely in bincode decode; leaving version 1 would contradict the documented schema-marker contract. |
| `ltnoyvll`, `nrzomklu` | Test-module ownership | Mechanical history repair, implemented for the current CFF ancestry | Declare each later module only in the fresh leaf that owns its file; keep CFF-only CLI/repeated-mass coverage before 4D and LTD comparisons after proper LTD. | Fresh leaves `tqsrsqrk`, `xxwzulrk`, and `lktlsypu` each own and declare their CFF CLI, repeated-mass, and scalar-battery modules. Only the proper-LTD comparison module remains deferred. |
| `kqqlnrnx` and descendants | CFF orientation ownership | Non-mechanical, decided | Port capabilities into main's owning `cff::orientations` abstraction rather than restoring the branch's parallel `cff::expression::GraphOrientation` implementation. | Main deliberately moved orientation ownership; the target expands the superseded location substantially. Repository guidance prefers one semantic boundary per concept. |
| `kqqlnrnx` and descendants | Feature dependency order | Non-mechanical, decided | Land raised-energy CFF support first, add 4D UV on top, and introduce proper LTD only later. | The maintainer requires CFF and 4D UV to stand independently of the later LTD route. |
| `kqqlnrnx` and descendants | UV, subtraction, and cut grouping | Non-mechanical, decided | Preserve main's disconnected-UV, E-surface classification, physical-cut exclusion, overlap validity, LMB normalization, and `CutGroup` invariants while porting the new routes. | Main's #93 rewrites these same paths after the branch point; taking either side wholesale would discard behavior. |
| `swnpzvnu` | Final split merge | Existing discrepancy, resolved by relocation | Keep the `GraphImportOptions` caller migration in its owning `tmlszmts` change; do not replay the old structural merge into the fresh linear ancestry. | The source merge's state reconciliation is already covered by `tmlszmts`. Its other topology reconciliation is superseded by the narrow CFF, QoL, 4D-route, and battery deltas; the stale source revision remains available until all deferred material is extracted and the obsolete history is retired. |
| _Rebase operation_ | Retained split DAG | Mechanical replay and conflict resolution complete; feature sign-off remains open | Retired the local monolithic change and forgot its local/remote-tracking bookmark without scheduling a remote deletion; rebased the retained source and replaced mixed descendants with fresh narrow changes on current `main`. | The initial operation reported conflicts in 15 of 16 changes. After abandoning the obsolete dependency-pin root, the selected fourteen-change ancestry through `lktlsypu` is conflict-free; obsolete source descendants remain separate and may still carry conflicts until final retirement. |

## Architectural decisions

| ID | Status | Decision | Rationale and consequences |
| --- | --- | --- | --- |
| A1 | Decided | Retain only the split/reviewable DAG; retire the alternate monolithic change and bookmark. | Both heads had identical trees. The maintainer does not need the monolithic review alternative, and retaining it would duplicate semantic conflict resolution. |
| A2 | Superseded | Do not preserve the original two-head topology. | The approved output has the split head only. |
| A3 | Superseded by A47 | The ledger was initially kept in its own change so early conflict resolution could not rewrite the review record. | Once the selected WIP ancestry became conflict-free, the maintainer explicitly requested that the live ledger be included in `raised_energy_cff_wip`. |
| A4 | Decided | Pin main's current Symbolica generation. | This preserves reproducibility and main compatibility without restoring `9135b409` and downgrading APIs/runtime behavior used by newer main commits. |
| A5 | Decided | Preserve main's owning CFF-orientation abstraction and disconnected-UV invariants. | The branch capabilities are ported onto the newer mainline model; branch-era parallel ownership is not restored. This is the primary semantic review area. |
| A6 | Decided and implemented | Relocate the non-empty final structural merge's `GraphImportOptions` reconciliation into its owning `tmlszmts` change and supersede the old merge. | This makes the history match its documented review structure; the remaining mixed merge content is replaced by the fresh narrow CFF, QoL, 4D-route, and battery deltas. |
| A7 | Decided | Enforce the dependency order CFF raised energy → 4D UV → proper LTD. | Earlier layers must remain usable and reviewable without the later LTD implementation. |
| A8 | Decided, revised by A28 | Keep only genuinely shared internal mathematics in the CFF layer. The public representation vocabulary may expose `Ltd` early only as an explicit unsupported setting; every LTD residue/builder route, metadata schema, evaluator, comparison tool, document, and numerical test remains in the later proper-LTD layer. | This keeps CFF and 4D UV independent of LTD without duplicating shared mathematics, while allowing stable configuration and CLI syntax to report a precise forward-compatible error. |
| A9 | Decided, clarified by A29 | The separately selectable 4D-local CFF route supports hedge-poset disconnected unions immediately through main's typed UV architecture. | Replace the branch's parallel local-4D forest implementation with main's implementation while retaining the route itself. Full component composition is required rather than an explicit-rejection first stage, and main's independent direct-3D route remains available. |
| A10 | Decided | Keep inspect JSON output and enforce main's read-only-state write boundary with the existing CLI helper and a focused regression test. | This retains the feature without introducing a second output-policy abstraction. |
| A11 | Decided | Repair QoL commit ownership so each intermediate change builds and tests independently. | Move graph-import plumbing/tests from `wrzsxkku` to `tmlszmts` and move 3D/LTD completion assertions to their actual feature stages. Symbolica load ordering remains a separate evidence-driven decision. |
| A12 | Decided | Make the top-level GammaLoop `bench` command the sole benchmark action, replace its legacy implementation with the branch's fixed-point benchmark engine, and keep `inspect` a single-evaluation action. | The richer engine has benchmark-specific lifetime, reporting, and failure semantics. This GammaLoop command preserves the established boundary, avoids an arbitrary benchmark return value from `Inspect::run`, and leaves room for explicit future modes without parallel command boundaries. Moving the selector layer requires checking for an existing reusable abstraction before introducing one. |
| A13 | Decided and implemented | Keep model-first state reload with main's tracing initialization; do not carry the obsolete Symbolica-first workaround. | The load-order comment and fresh-process regression live in `wrzsxkku`. The child reload preserves quoted one-character UFO printing, emits none of the obsolete duplicate-function warning, and passes against the final main-based tree. |
| A14 | Decided | Raised CFF must preserve main's spinney-owned numerator-growth model exactly. | The complete prepared numerator may be used only for the up-front scaling analysis. Recursive UV state starts without an eagerly assembled full numerator; at each `given` to `current` transition it grows by exactly the newly owned graph numerator fragment, then attaches the outside/global factors once at final assembly. The historical root-only injection of `production_numerator_atom_for_full_3d_expression()` is rejected because it bypasses this ownership model and gives recursive children different numerator semantics. |
| A15 | Decided | Raised-CFF numerators remain factorized and are evaluated with parametric energy signs and complete affine energy shifts. | Every factor in one product uses the same exact numerator map; exact maps remain distinct even when they share an `EdgeVec` orientation. Internal OSE coefficients/signs, external-energy shifts, uniform sampling-scale shifts, and constants survive until a UV series or terminal evaluator mathematically requires materialization. UV rescaling/localization state remains a separate transformation layer. Expansion of the full numerator at the CFF/UV ownership boundary is forbidden. |
| A16 | Decided | Generate against the degree envelope of the eventual complete numerator while allowing the attached numerator to grow incrementally. | Capacity is independent of the factors currently owned by a spinney. The up-front analyzer consumes the factorized full product compositionally rather than expanding it. Growth may validate that the accumulated strict degree stays within the declared envelope; this validation is principally a debugging invariant, not a production regeneration mechanism. |
| A17 | Decided; source-local identity corrected after integrated-NLO diagnosis | Add one private exact-map projection operation and retain map-indexed local-3D terms through a single UV forest traversal. | `Graph::numerator` remains the only physical-factor ownership mechanism. Each additive carrier branch holds a production `OrientationID` selecting one complete generalized residue-map key and, when generated by a contracted source, the source-local `edge_energy_map` that owns numerator substitution. Physical edge directions remain metadata and cannot merge distinct production IDs. Distinct `(selector, source map)` pairs are unioned rather than cross-multiplied; the same pair may be zip-added. Root production-owned branches retain their selector partition and strict exact-map identity without normalization. A collapsed local or integrated source residue is averaged over its selected direction-compatible production entries by `1/N`, but its affine/contact map need not equal any production map. This preserves maps through localization, active-sector joins, and final assembly without recomputing the forest once per map. |
| A18 | Partially approved; semantic implementation reopened by A50 | Keep `Atom` as the sole expression tree and derive a flat canonical `EnergyPowerSupport` without `.expand()`. Sparse correlated support, factorization, parametric signs/shifts, and spinney-owned numerator growth were knowingly approved. | A18 incorrectly bundled representation shape with variable identity. The LMB rewrite and removal of edge-bound remapping were presented as an exact/lossless refinement, without disclosing that numerator energy powers would cease to be indexed by EMR/source-denominator energies. That semantic replacement is not treated as informed approval. |
| A19 | Required by A7/A8; staged route and CFF battery implemented, final validation blocked at A45 | Build later stages as deltas in the order CFF core → separate typed 4D-local UV route → proper LTD; never land the sibling `rpspmzvw` tree wholesale. | Main's `local_4d.rs` replaces the implementation of the alternate `expanded_4d.rs` route without replacing main's direct-3D route. The 4D bridge splits exact term-level denominator occurrences: zero-active-loop factors remain typed-4D-owned residuals, loop-dependent factors stay active raised CFF, and exact multiplicity/form prevents double counting. The CFF battery is present in `lktlsypu`; proper LTD follows as a preserved separate semantic delta. |
| A20 | Superseded by A78 | Keep the nullity-aware certificate `N-2|A| < -d`, checked arithmetic, and typed arrangement-size errors, and fail conservatively when the energy-contour condition is not proved. A post-generation spatial UV profile is corroborating evidence, not a replacement for the energy-contour condition. | This was the historical bespoke hard-gate design, subsequently refined by A52/A62. A78 separates the required EMR rank envelope from source convergence and replaces that machinery with the existing all-cycle UV-DOD calculation. |
| A21 | Numerical boundary retained; mechanism superseded by A55 | Uniform-`M` sampling remains sparse and channel-local, preserves correlated support, and leaves residual affine variables exact. | The live one-CFF implementation realizes this through repeated-channel interpolation and normal form. A56 distinguishes that legitimate simultaneous channel sample from finite-pole division, where unrelated EMR arguments must stay fixed; A65 records the remaining activation and `M != 0` boundary issues. |
| A22 | Decided and implemented | Restore a contracted parent-edge energy map only when its exact source coordinates have zero inner-loop component; construct it from the orientation's outer loop map plus the parent external signature. Never invent an outer map for an inner-dependent contracted edge. | Generic index remapping uses zero both as a valid expression and as a missing-slot fill. `GraphThreeDSource` now exposes only the reconstructible outer suffix; the lift uses exact rational scaling and actual parent external-edge IDs. Production tree contraction has no inner source loops, while nested inner-dependent holes remain unavailable. A focused source test checks both sides of this boundary. |
| A23 | Decided and implemented; valid-topology focused coverage pending | Represent only the scalar quotient of a final quadratic energy denominator as a Unit lower sector. If consuming the last carrier leaves a positive free-loop-energy power, return `UnlocalizedEnergyContact` before emitting any partial terms. | This is an algebraic local-contact boundary, not a prescription for assigning an integrated value to `∫dq0 1`. The implementation and terminal quartic Unit/contact path remain live, but the historical dedicated None/All degree-two-through-four regression used a one-edge self-loop. It was not restored as a topology oracle; a balanced replacement fixture remains to be approved. |
| A24 | Superseded by A44/A55; top-level/lower-sector rank boundary retained by A58 | The earlier rank-projected lower-sector implementation established the required nonintegral-projection, null-direction, initial-cut, and correlated-support boundaries. Its private raw repeated-residue machinery is no longer the production architecture. | `LowerSectorRankProjection`, `RepeatedResidueBuilder`, and their raw `q±E`/H-surface construction remain intentionally removed. A58 restores full rank in the ordinary top-level builder and keeps deterministic particular lifts only in the live lower-sector component builder, without restoring a second residue engine. |
| A25 | Deferred proper-LTD schema boundary, updated by A44 | Current CFF threshold selection owns only canonical positive-energy residue conventions. Proper LTD must later introduce its own metadata owner for Cutkosky combinatorial coordinates, local-series surface signs, per-cut denominator alternatives, and prefactor bridges. | Threshold and true LTD/LU bridges are not interchangeable. Live `denominator_edge_support_signs` remains neutral CFF merge/cut-support provenance, not a partial LTD implementation. No unused LTD builder is retained; the pushed `proper_ltd_core` snapshot preserves deferred source material until it is ported as a reviewed delta. |
| A26 | Decided and implemented | The first `3Drep` CLI surface is CFF-only: `validate` and `build` consume current GammaLoop graph state and exact numerator-energy support. A direct sampling-normalization option has invocation-local precedence over the global setting. Every build regenerates; artifact directories and the `latest` pointer are outputs, not a semantic cache. | One private borrowed graph-catalog adapter spans generated/imported amplitude and cross-section owners without adding a public graph abstraction. The historical cache did not key the actual graph or numerator and is removed, together with `--clean`. Evaluate, reconstruction, manual bounds, and LTD implementation remain absent. A28 later adds representation selection solely to return an explicit LTD error. |
| A27 | Energy-factor ownership retained; exact repeated normalization corrected and validated by A63 | Generated bounded variants that retain half-edges own their complete local normalization; only an ordinary global-source CFF result receives the caller-owned `prod(-2 E)^-1`. Extension averaging remains part of direct-3D production-orientation localization and never applies to post-4D exact sources. | The former same-routing/opposite-routing coefficient claim came from the now-removed mass-jet organization. A63 derives the replacement ordinary-CFF oracle from the contour: two `+i/(64 pi^3 omega^3)` orientations give the explicit sum `+i/(32 pi^3 omega^3)`. The corrected tests now validate that complete sum rather than assigning the total to each raw orientation. |
| A28 | Decided and implemented | Raised-power CFF remains the only implemented representation. Add serializable `RepresentationMode::{Cff, Ltd}` and the matching CLI build setting now, default to CFF, and return a typed `NotImplemented` error for LTD before doing representation work. Do not add proper-LTD metadata or numerical paths yet. | This preserves a stable setting without implying partial LTD support. GammaLoop's production CFF constructors state `Cff` explicitly; the CLI rejects LTD before graph selection and output writes; the core rejects it at every public generator entry point. Focused core and CLI regressions pass. |
| A29 | Maintainer correction, implemented; final route parity blocked at A45 | Preserve two CFF UV construction routes: main's direct local-3D recursion and an independently selectable 4D-local forest implemented with main's typed `Local4dCts`, followed by CFF projection. Replacing `Expanded4DApprox` replaces only the old implementation. | The default false setting retains main's recursion and spinney-owned numerator growth. The true setting projects only after typed 4D local composition. Focused construction tests pass, but the GL11 scalar battery currently exposes the distinct-physical/coincident-denominator projected-route defect recorded under A45. |
| A30 | Maintainer-confirmed terminology and scope | Name the second route “4D-local UV, then CFF”; do not conflate it with a final evaluable pure-4D parametric integrand. | Source `b8bec534` computes an expanded 4D local forest but explicitly calls `expanded_4d_terms_to_3d_parametric_integrands`; main's `FinalIntegrandDimension::FourD` remains rejected by `UVOrchestrator::parametric_integrands`. Final pure-4D evaluation is outside this stage. |
| A31 | Boundary retained; terminology and mechanism superseded by A55 | Raised-energy requirements are a property of the active CFF source and its exact source-owned numerator-energy support, never of whether an LU selector is present. Production and direct-3D local/LU CFF therefore use the same causal source family and main's canonical positive-energy `select_esurface_residue(lu_cut)`. | The earlier “confluent versus nonconfluent” split described an intermediate implementation, not two supported CFF representations. A55 keeps one causal CFF implementation and handles raised repeated sectors through channel interpolation/normal form, while preserving spinney-owned numerator growth and adding no LTD metadata. |
| A32 | Abstract-index normalization regression, mechanically restored to main | Preserve main's distinction between an abstract `mink(4,n)` dummy-index label and an explicit `cind(n)` component selector. Only `Q3(e,cind(n))` may normalize to a concrete scalar component; `Q3(e,mink(4,n))` remains an indexed spatial vector. Reject the earlier terminal-dot-materialization proposal and keep exact energy mapping at tensor level. | `kqqlnrnx` had broadened the `Q3` normalizer so the numeric dummy label in `mink(4,1)` was mistaken for component 1. Exact mapping then malformed `Q3(e,μ)+A_e δ(cind(0),μ)` into `Q(e,cind(1))+A_e δ(cind(0),mink(4,1))`, a scalar-plus-vector sum. Main, every preserved source revision, and the current energy-support analysis use `cind`-only component semantics. Restoring that invariant lets the GL24 quartic direct-3D integrand complete evaluator construction; the run proceeds to the separately tracked 4D exact-map-extension failure. A focused normalization regression parses the exact `spenso::mink(4,1)` syntax and structurally verifies that `Q3` and its abstract slot survive; it passes 1/1. No early dot expansion is justified. |
| A33 | Physical cut-group membership, mechanically implemented under the retained-all-cuts contract | Build each raised E-surface group's `cuts` from the explicit intersection of all physical `cut_esurface_map` entries with all E-surface IDs in the group; do not infer physical membership from only `esurface_ids[0]`. | The raised generator can place a nonphysical generated E-surface first in a group that also contains a physical Cutkosky surface. The old first-element heuristic then omitted that physical cut, while its one-to-one reverse `HashMap` also collapsed multiple physical cuts mapped to one E-surface. `CutGroupData` now filters the ordered physical cut map against the complete raised group and skips only an empty intersection. A focused regression covers both failure modes and passes; this removed the GL02 missing-cut failure before A42 replaced the intermediate raised-carrier mechanism. |
| A34 | Maintainer-approved and implemented | Use two ID domains: `TopologicalThresholdId` indexes a persisted exact topology inventory, while `RaisedEsurfaceId` indexes expression-owned residue groups. Build the topology inventory independently of settings; translate eligible threshold counterterms locally against the frozen expression cache by exact match first and repeated-edge-normalized match second, rejecting ambiguity after conversion to raised-group identity. | Main's classic CFF happened to enumerate every threshold topology and hid the ownership distinction. Raised production does not. The dual boundary keeps exact surfaces and dense IDs for classification/API while counterterm construction uses only expression-owned groups. Stored process/integrand binaries intentionally change format under repository policy; JSON/Python field names remain stable. All-target checks and focused API association/dense-ID/cut-exclusion regressions pass. |
| A35 | Superseded by A77; earlier ordinary-root rejection retained as chronology | Both UV paths build the complete typed 4D forest because integrated counterterms require it. The earlier decision required the projected path to map both local counterterms and the root through exact 4D CFF. | Main and the validated pre-rebase path instead assemble the empty forest from the ordinary production root expression. A77 records the later product-semantics correction: expanded-4D projection changes proper UV nodes, not the original integrand identity. |
| A36 | Source-coordinate repair implemented; production-ID objective superseded by A38 | Exact 4D sources use the production coordinate frame, deterministic occurrence-local IDs, exact masses/OSEs, and unique loop-carrier labels. Preserve full affine signs and shifts, but do not infer production-map identity from that shared coordinate frame. | The coordinate and dual-ID repair makes shifted, duplicate, and UV-basis exact sources deterministic without collapsing occurrence semantics. GL24 then demonstrates the stronger boundary: a transformed source can have a valid source-local residue decomposition with no meaningful production `OrientationID`. A38 keeps the source repair and removes the invalid identity requirement. |
| A37 | Maintainer-approved and implemented under A38 | Analyze each completed 4D term with the factorized product of its term-owned numerator and the outside/current-spinney plus global numerator. Grow the outside factors only after the 4D forest has completed, immediately before term-local CFF projection; then map that complete factorized product with one source-local affine energy map. Never expand it and never remap the outside factors later through a production ID. | Main's spinney-owned numerator growth remains unchanged during 4D construction. At projection, one owned exact mapper rewrites the whole product, including abstract temporal tensor components, while leaving concrete spatial components spatial. Focused coverage verifies multiplicative factorization and every affine sign, OSE, external shift, uniform sampling shift, and constant. |
| A38 | Maintainer-approved, implemented and focused tests passing | `4D-local UV, then CFF` is strictly an explicit-orientation-sum route. Every transformed exact CFF owns its complete source-local orientation sum; there is no source-to-production `OrientationID` match, extension averaging, or uniform outer embedding. Require `explicit_orientation_sum_only`, reject a generation orientation pattern and runtime individual-orientation sampling through the existing validation boundary, and leave direct 3D fully production-orientation-parametric. | Projected branches remain separate as serialized deferred bodies through symbolic rewrites and are lowered only at evaluator construction. Production orientation patterns do not filter the source-local sum. Amplitude, cross-section, threshold, and LU consumers use the same boundary; per-production-orientation profiling is rejected. The direct-3D route retains its production-map IDs and ordinary orientation selection. |
| A39 | Maintainer-approved and implemented | Exact post-4D energy mapping uses owned source coordinates and occurrence-local energy IDs, not production identities. Validate numerator support in the complete rewritten source basis before discarding inactive coordinates. Initial-state cuts are literal signed aliases of their physical external edges even when a stored loop row has inactive components. | Complete-factor rewriting permits inactive placeholders to cancel without expansion; genuine inactive dependence still errors. Exact OSE occurrences, rational parent coordinates, external signs/shifts, factorized products, and tensor index structure survive. A zero-active-loop pure tree has mathematical loop-signature rank zero and is supported without invoking a zero-column Symbolica matrix. |
| A40 | Maintainer-approved and implemented | Represent finalized projected sums as a compact expression plus serialized source-local bodies, then register evaluator-local tagged `GS.projected_cff_sum(tag)` functions. Run every deferred body through the same color, Dirac, metric, dot, tensor-network, aliasing, and contraction preprocessing as an ordinary evaluator atom before function registration. | This avoids materializing a potentially large orientation sum while preserving full spin/color/tensor semantics. The one private `EvaluatorStack::preprocess_atom` method is the former inline production path, reused by both ordinary atoms and deferred bodies. A closed bispinor metric trace proves deferred and ordinary evaluators both reduce to 4. |
| A41 | Implemented under A28/A38 and refreshed by A76 | Preserve `Ltd` configuration as an explicit future setting but reject `runtime.general.use_ltd` immediately while the backend is CFF. Explicit-sum stacks evaluate exactly once regardless of the requested summed evaluator, reject runtime orientation filters/sampling, and archive this mode in GammaLoop-generated evaluator outputs. | Proper LTD is neither deleted nor partially emulated. Main #95 and this WIP independently changed the same positional archive layouts under versions 5/8; the combined amplitude and cross-section formats are therefore version 6 and 9. Generated loaders use main's current Symbolica stack and reject per-orientation requests clearly. |
| A42 | Public boundary retained; repeated-pole mechanism superseded by A55 | Keep exactly one public CFF implementation in `three-dimensional-reps`; `gammalooprs::cff` constructs sources and adapts surfaces, cuts, and evaluators but does not own a second generator. Numerators remain factorized with parametric affine energy signs and shifts. | The mass-derivative/numerator-jet organization recorded here was an intermediate repair and is now rejected for CFF repeated sectors. A55 restores the validated CFF organization: retain occurrence powers and apply channel interpolation/normal form before ordinary lower-CFF generation. The raw repeated-residue/H-surface engine remains removed and proper LTD remains deferred. |
| A43 | Ownership model retained; post-A55 assertions under review | Make inverse-energy normalization ownership explicit in `GeneratedThreeDExpression`: `GlobalSourceProduct` means the GammaLoop adapter appends the source-wide product, while `VariantLocal` means every generated variant already owns its half-edge factors. | This prevents applying the inverse-energy product twice. After A55, ordinary non-raised repeated sectors correctly report `GlobalSourceProduct`, while bounded channel-normal-form sectors remain variant-local. Two new tests still assert the removed occurrence-jet ownership and require explicit review; A63 separately reopens the exact repeated sign/magnitude oracle. |
| A44 | Maintainer-approved and implemented | Remove dead raw repeated-residue code rather than retaining a second hidden implementation. Preserve future LTD through the explicit unsupported setting and the separate pushed source snapshot, not through unused production structs or adapters. | `RepeatedResidueBuilder`, `LowerSectorRankProjection`, raw `q±E`/H-surface construction, unused `rising`, zero-caller adapters/types, and the dead public duplicate-sign setting were removed. Live internal routing/support-sign mathematics remains with its causal owner. Proper LTD remains available in `proper_ltd_core` and must later be ported onto the post-4D head. |
| A45 | Historical topology model; owner-sensitive incidence distinction superseded by A79 | Parser-expanded denominator powers are canonical serial chains. Exact energy identity remains occurrence-local, with an explicit many-to-one projection to physical source/cut ownership. | The earlier source-wide “all exact occurrences are parallel” rule was too coarse because powers of one source wrapper require serial subdivision on the known source scaffold. The later claim that coincident denominators on distinct owners must retain owner-distinct exact incidence is withdrawn. A79 instead retains those owners as provenance while normalized rewritten denominator channels and multiplicity determine an owner-invariant exact residue topology. Cograph/UV minors, same-owner power-chain bookkeeping, crown hedges, and post-construction canonical relabeling implement that boundary without generic graph reconstruction. |
| A46 | Boundary retained; exact owner-invariance assertion refined by A57/A79 | Test only the live causal generator boundary, preserve topology meaning independent of list rotation, treat isolated zero-edge components as multiplicative identity, and make parser-powered repeated fixtures actual serial chains. | The serial assertion compares topology without depending on edge-vector order, and CFF recursion filters isolated vertices instead of manufacturing a zero-valued causal surface. Repeated copies of one source wrapper use its known scaffold incidence for deterministic subdivision, but compatible relabeling among owners of the same rewritten rational channel cannot change the final exact incidence, loop rank, or residue. |
| A47 | Maintainer-requested and implemented | Include this review/readiness ledger directly in `raised_energy_cff_wip` instead of rebasing it later as a separate final documentation change. | The file now travels with the WIP it evaluates. Detached `zzxptpwq` is retained only until cleanup so no historical text is lost; it is not intended as a second final ledger revision. |
| A48 | Mechanical rebase repair, implemented and smoke-tested | For the global cross-section CFF catalogue, contract `tree_edges - initial_state_cut` before generalized CFF generation, exactly as current main and the retained pre-rebase caller do. Keep those propagators in the caller-owned `tree_denoms`; do not weaken the generator's zero-loop lower-sector guard. | The raised-CFF integration had replaced the contraction list with `&[]`. In `epem_a_ddx` GL0 this admitted external-only tree edges 4 and 6, whose zero loop rows reached the generic lower-sector guard. Restoring the caller boundary makes the command complete without double-counting tree denominators, but does not validate the LMB-indexed numerator support reported by the current implementation. |
| A49 | Maintainer-approved residual-denominator restoration implemented and focused-validation complete | Keep GammaLoop's production bridge ownership unchanged: graph-theoretic bridges are caller-owned `tree_denoms`, not causal-graph denominators. Separately preserve the validated pre-rebase auxiliary `3Drep build` behavior by returning typed `ResidualDenominator`s and projected tree factors rather than rejecting those graphs. | The Graph-backed auxiliary caller, rather than arbitrary `ParsedGraph` core input, explicitly selects paired, non-dummy `tree_edges - initial_state_cut` by physical ID. `generate_3d_expression` remaps those IDs locally, rejects cut or loop-carrying requests, contracts them from the causal graph, projects bounds to the active graph, then lifts maps/surfaces and typed residual factors back to physical IDs. Pure trees emit one identity orientation plus residuals. Display, serialization, remapping, and numeric evaluation of `prod_e (q0_e^2-E_e^2)^(-p_e)` are live. GammaLoop production rejects residuals crossing its separate caller-owned `tree_denoms` boundary. The exact pre-rebase external-tree and pure-tree CFF regressions, a new numerical residual-factor regression, formatting, checks, and warnings-denied shared-core Clippy pass; no existing assertion changed. Historical symbolic `to_atom` methods remain CFF-only and do not absorb the typed residual list. |
| A50 | Maintainer-approved historical repair implemented; raised-channel validation remains | Undo the unapproved semantic part of `58a1634e`: numerator energy-power bounds belong to original EMR/source-edge IDs and must be remapped through `EnergyEdgeIndexMap` plus every component/contraction embedding. Keep factorized `Atom` analysis without expansion, but do not use LMB monomials as the capacity contract. | Optional edge-indexed bounds, exact source-edge remapping, partial supplied numerators, and strict rejection of genuinely unrepresentable shrunken-edge dependence are restored. `None` means the legacy affine class; `Some([])` means explicitly energy-constant. A declared `K(loop,0)` is normalized only through its explicit LMB carrier while stored ownership remains that physical EMR edge. Production/test code no longer references the rejected LMB support API. |
| A51 | EMR-identity boundary retained; repeated implementation superseded by A55 | Treat each original EMR `edge_q0` entry as the authoritative black-box numerator argument and identity. Repeated-denominator algebra never collapses numerator arguments; generalized samples may lack a global `loop_q0` lift. Keep deleted-edge finite-pole samples separate from lifted affine values. | The pinned [theory note at revision `65ab03e4`](https://github.com/alphal00p/generalised_ltd/blob/65ab03e4fb7d442ff362392df2c2a59ef323d28c/docs/generalised_ltd.tex) and dispatcher confirm that CFF retains denominator powers and performs channel/contact normal form; it does not use a separate repeated-pole mass derivative. A55 implements that correction and A56 restores edge-owned finite-pole sampling. A64 records a remaining pre-existing numerator-call deduplication mismatch that the maintainer deliberately deferred beyond this rebase. |
| A52 | Maintainer-approved and implemented, deliberately conservative | Fail generation when the contour-at-infinity certificate is not proved, and fail exact/shrunken occurrence generation when nonconstant physical-EMR numerator capacity lacks an affine numerator-only channel. Keep both failures explicit and revisit them when stronger proofs/maps exist. | The certificate is edge-bound and nullity aware. Exact 4D terms analyze their actual factorized owned numerator and allow explicit energy-constant `Some([])` terms. The exact `epem_a_ddx_NLO.toml` command passed after a fresh release build, while its configuration audit records that it does not exercise integrated 4D-to-CFF projection. |
| A53 | Maintainer-corrected production bridge boundary approved and implemented; auxiliary-CLI boundary resolved by A49 | Keep main's structural ownership boundary: `tree_edges - initial_state_cut` are graph-theoretic bridges, are contracted before CFF, and are projected separately to 3D with their denominators retained in caller-owned `tree_denoms`. Exclude their physical-EMR degree bounds only from the reduced CFF capacity domain; keep the complete factorized numerator and external-affine physical `edge_q0` map for later projection/evaluation. Preserve `Some([])` if this removes the last bound. Do not create a CFF numerator-only/shrunken exception for bridges. | Direct and exact-4D generation classify bridges structurally from non-cut `tree_edges`, verify that a contracted bridge has exactly zero reconstructed outer-loop coordinates, and remove only cut/bridge IDs from the reduced CFF capacity set. The denominatorless exact-term branch now applies the same exclusion before using its empty-source affine mapper; residual exact tree occurrences remain outside CFF. The complete numerator and existing physical map restoration/projection remain untouched. Bounded contracted **non-tree** edges still reject, including coincidentally zero-coordinate rows. The finite-CT map, external-tree, pure-tree, and strict nonbridge tests pass without assertion changes. GammaLoop retains `tree_denoms`; auxiliary `3Drep build` generation uses A49's explicit typed residuals. |
| A54 | Superseded by A55; rank boundary retained by A58 | The temporary occurrence-collapse repair used `LowerSectorCffBuilder` when collapse lowered rank, removed obsolete LMB/confluent tests, and retained EMR-edge coverage. | A55 removed occurrence collapse entirely, so the binary ordinary-versus-collapsed repair is no longer the production architecture. A58 now requires the ordinary top-level active denominator rows to span the declared loop space and confines particular lifts to genuine lower-sector construction. Historical assertions changed during the intermediate phase are catalogued in A59 and have not been silently accepted. |
| A55 | Maintainer-approved and implemented; integration validation incomplete | For **all** repeated CFF sectors, keep the full occurrence graph. Non-raised sectors use the ordinary causal CFF builder; raised sectors use `BoundedCffBuilder`, whose channel interpolation/normal form acts before deleting copies and then generates ordinary lower CFF trees. Do not restore `RepeatedResidueBuilder`, H-surface/LTD machinery, or a second public “confluent CFF” entry point. Remove the CFF occurrence-collapse/mass-jet path rather than retaining two repeated-pole organizations. | The pinned note's `channel-interpolation` and `channel-normal-form` equations explicitly require this ordering. Its `sunrise_pow4`, `edge 2:5` example gives `y^5/D^4 = omega^4 y/D^4 + 2 omega^2 y/D^3 + y/D^2` and lists the expected rows. The pre-`58a1634e` dispatcher and validated pre-rebase shared CFF core use this organization. The current implementation restores the old logical-channel/known-factor normal form, removes the mass-jet/occurrence-projection dead route, compiles, passes warnings-denied Clippy, and generates the former degree-five sunrise blocker. A56 now also restores the validated finite-pole sampling semantics. The refreshed evaluator-enabled inventory remains 53/59 after A49/A56/A58/A62; A59, the strict A52 fixture boundary, and A63 normalization remain open. |
| A56 | Maintainer-approved and implemented; focused validation complete | Finite-pole interpolation varies only the selected physical EMR argument. All other `edge_q0` entries remain fixed, even when the selected edge is the named exact LMB carrier; `loop_q0` remains diagnostic and need not reconstruct generalized contact samples. | The theory's quadratic division explicitly holds every other EMR edge energy fixed, and its JSON semantics make `edge_q0` authoritative because generalized samples need not have a global loop-energy lift. The validated pre-rebase core (`b9c4e2ff`, also integrated unchanged in `cc8c5979`) cloned the existing edge map and replaced only the active edge in quadratic remainder/contact branches. Sibling snapshot `0c8d166b` introduced the incorrect carrier rewrite independently of `58a1634e`. The repair removes that helper, selector, error, and sampled-variable plumbing; remainder/contact branches retain their diagnostic loop map, copy their authoritative edge map, and overwrite only the selected sample. The approved A23 zero-edge scalar contact no longer requires an exact carrier, and terminal known-factor replacements again remain edge-only. Repeated-channel interpolation is deliberately unchanged because it replaces all signed aliases of one physical channel together. With approval, the carrier-recomputation test was replaced by a production-path regression that compares the full spectator edge-map vector and loop map against the degree-reduced source. That regression plus focused quadratic remainder/contact, quartic Unit/contact, and mixed known-factor tests pass. The `eval` feature's separate loop-vector fallback is not part of this production repair; see A64. |
| A57 | Historical owner-aware chain key; provenance bookkeeping retained and owner-sensitive incidence superseded by A79 | Keep the validated shared CFF core's topology contract of local occurrence IDs plus explicit serial incidence, with no physical-owner field added to `ParsedGraph`. | The source boundary, not the shared `ParsedGraph`, owns provenance. GammaLoop uses `source_edge` to seed and validate the cograph/UV-minor scaffold, retain cut support and occurrence-to-physical maps, and collect repeated copies of one wrapper for deterministic serial subdivision. Exact signature, mass, domain, and multiplicity define the rewritten rational channel presented to the owner-free shared core. Changing only compatible owner labels may change provenance maps, but cannot change that channel's exact incidence, loop rank, or contour residue. |
| A58 | Maintainer-approved top-level/lower-sector rank boundary implemented and focused-validation complete | Ordinary top-level pure CFF must choose a full residue basis from active denominator rows `D`. Initial-state-cut entries are external aliases: they neither add denominator rank nor supply equations that complete a deficient residue basis. Rank-deficient particular lifts remain valid only inside lower-sector builders whose surrounding residue construction already owns the removed directions. | Top-level pure CFF now greedily selects the same first active matroid basis as the pre-rebase ordered-combination search, requires `basis.len() == loop_names.len()`, restores the full-loop sign, and calls `solve_loop_energy_from_target_edge_exprs`; a missing direction returns `SingularBasis` before orientation enumeration. `LowerSectorCffBuilder` still projects each component to its own rank and uses `solve_loop_energy_particular_from_target_edge_exprs` only while lifting those proven components back to parent coordinates. The invalid momentum-imbalanced top-level cut fixture was replaced, with approval, by a balanced non-self-loop source: public generation rejects even though cut rows span the missing coordinate, while direct lower-sector construction succeeds, zeros only its diagnostic free coordinate, and keeps both cut energies as external aliases. Both tests, shared-core checks/Clippy, the GammaLoop library check, the 53/59 diagnostic inventory, and the exact NLO command pass at this boundary. |
| A59 | Diagnostic inventory recorded; regression classification pending; not a GammaLoop product gate | Use the evaluator-enabled serial inventory to compare the shared core with the selected pre-rebase baseline. Only newly introduced failures or failures reaching a GammaLoop call path block the rebase; resolve production defects first and ask before changing historical or newly introduced assertions. | The refreshed `cargo nextest run -p three-dimensional-reps --features eval --no-fail-fast --test-threads 1` inventory runs 59 tests in isolated serial processes: 53 pass and 6 fail. The raw count alone does not imply GammaLoop unreadiness. The approved A58 replacement removed the invalid top-level cut/ownership failure and added two passing boundary tests. Of the six remaining failures, two ownership assertions expect the removed occurrence-jet ownership where ordinary CFF reports `GlobalSourceProduct`; two pre-A55 high-power tests generate successfully but still assert the post-`58a1634e` recursive origin instead of their original known-factor origin; and two evaluator tests use cap vectors rejected by the approved strict A52 contour gate before their analytic internal-algebra checks. The approved A56 carrier test has now been replaced and passes; it was not one of those six failures. No remaining failing expected value, ownership assertion, origin assertion, or graph fixture has been edited. Four exhaustive test-option literals received only the empty restored A49 option field. The A49/A56/A58/A62 regressions pass. |
| A60 | Maintainer-authorized integration-workspace version boundary implemented and validated | Version integration-workspace resume data independently of `state_manifest.toml`. Reject old, missing, or future workspace schema versions before decoding `state/integration_state.bin`, return contextual errors instead of panicking, and do not add a positional-bincode migration. | `IntegrationWorkspaceManifest` now has a required version-1 boundary checked before both normal resume and summary decoding. Missing, older, and future formats report contextual errors with `integrate --restart` guidance; state I/O and bincode failures retain their paths instead of being treated as absent state or panicking. Focused manifest, normal-resume, and summary-resume rejection tests pass. |
| A61 | Maintainer-authorized unique-literal EMR mapping implemented and focused-validation complete | Translate a physical EMR bound to exactly one exact occurrence only when that occurrence momentum is literally the same `Q(edge)`. Keep missing, ambiguous, shifted, and general-affine mappings as typed errors; do not use LMB coordinates or ownership as energy identity. | `ExactSourceEnergyMapper` now owns the physical-to-occurrence capacity translation, and exact 4D-to-CFF generation passes its mapped occurrence bounds to the shared CFF core while preserving the complete factorized numerator. Focused coverage accepts one literal occurrence and rejects absent and repeated literals. The projected GL11 path now reaches evaluation; its disputed historical numerical oracle remains a separate A73 test-intent issue, and no GL11 snapshot or fixture changed under A61. |
| A62 | Historical implementation superseded by A78 | Separate the actual denominator arrangement from numerator-owned EMR capacity. Let `D` contain only active CFF denominator occurrences, excluding initial-state cuts and caller-owned bridges while retaining occurrence multiplicity; derive active rank/nullity from `D`. Let `N` contain physical EMR/source rows with degree bounds after cut-alias, bridge separation, and shrunken-edge semantics. For every denominator flat `F`, use `s = rank(D) - rank(F)` and certify `sum(active N bounds) + s < 2 * count(active D occurrences)`. | This records the former bespoke contour-certificate implementation and its tests. A78 removes that machinery, retains the independently required factorized EMR energy-bound analysis, and gates source integrals through the simpler existing UV-DOD calculation. |
| A63 | Contour normalization retained; owner-sensitive incidence distinction withdrawn | Retain the ordinary pure-CFF duplicate-excess sign for algebraically repeated denominators and validate the **explicit orientation sum**, not each raw orientation against the total. | The contour result `+i/(32 pi^3 omega^3)` for the complete double-pole sum remains valid. Repeated copies of one source wrapper still use same-owner serial-chain bookkeeping on the source scaffold. Equal rewritten rational denominator multisets must, however, have the same exact incidence, loop rank, and explicit-sum residue after neutral owner aliases are identified; distinct owners survive only in cut support and physical projection metadata. |
| A64 | Maintainer-deferred beyond the rebase; shared-core design retained | Future repair: key numerator-call labels and compatible-variant grouping only by canonical authoritative `edge_q0`; retain the first deterministic `loop_q0` as diagnostic data when equal edge maps merge. Keep different edge maps distinct. Do not add an LTD path or a GammaLoop-only duplicate coalescer. | The mismatch predates the rebase and exists in the validated pre-rebase shared core: `assign_numerator_map_labels` and three live variant-push paths use both loop and edge maps, although GammaLoop evaluates numerator energies from physical edge maps. No wrong GammaLoop value is currently demonstrated; the known cost is duplicate numerator calls and expression swell. The retained shared-core proposal would leave `eval.rs` unchanged, but an unusual auxiliary feature-gated evaluator `loops[]` expression without a named carrier could observe only the retained first diagnostic lift after rows merge. This pre-existing auxiliary-evaluator concern is not a rebase blocker. The maintainer chose to defer that behavior decision while keeping the proposal in this ledger. No A64 production code, evaluator code, or test has changed. LTD still fails immediately through the approved unsupported-mode boundary. |
| A65 | Runtime uniform-sampling validation implemented and focused-validation complete | Reject `M = 0` exactly when a finalized generated evaluator uses the auxiliary numerator sampling scale. Preserve the existing `All` activation behavior rather than changing a test-backed, algebraically harmless degree-one policy during this repair. | Finalized evaluators record whether their emitted expressions use `M`, including through archive serialization. Warm-up accepts any nonzero scale, rejects zero only for evaluators that use `M`, and leaves independent expressions unaffected. Shared-core and GammaLoop tests verify higher-power value invariance under two nonzero scales, archive preservation of the usage bit, and the conditional zero-scale error. |
| A66 | Validation gap and test-intent decisions pending | Port fixed serialized-row and fixed-kinematics numerical CFF goldens from the external generalized-LTD reference, including unique-edge-map/one-call, quadratic and multiloop cases, repeated cubic/quartic/sunrise/iterated-bubble cases, UV-convergent cap coverage, uniform-scale invariance, and the exact double-quintic UV-negative case after A62. Keep this CFF-only; no proper-LTD runtime dependency is required. | Relative to the CFF-only Rust baseline, no structural CFF test was lost; the renamed repeated-bubble case has a causal replacement, and A49 has restored the external-tree/pure-tree regressions. Relative to the full validated source, the remaining non-LTD losses are the graph-signature tests in A69, the double-quintic rejection, and numerical comparisons that can be replaced by fixed CFF goldens. The external validation evidence was never made a repository CI gate, so strategy assertions can miss value changes. Two analytic contact tests now stop at A52 before their intended algebra; making their fixtures UV-convergent with spectator denominators is preferred to a public unchecked route, but every test edit requires explicit approval. A59 lists the immediate assertion changes. |
| A67 | Mechanical interface maintenance complete; test-affecting benchmark expansion deferred | Keep committed schemas, live example pins, read-only output guards, Python lint suppressions, and user-facing command documentation synchronized with the implemented API. Do not silently alter default-output or benchmark-test intent. | The three committed schemas were regenerated, the live Rust example uses main's Symbolica source, and read-only guards cover explicit renormalize and approach outputs. The four Python integration-helper suppressions were audited and retained: their only caller remains inside the legacy commented method, and warnings-denied Python compilation confirms they are not stale. The obsolete Python subgraph `generate_cff*` endpoints were removed with their retired GammaLoop-side generator, and orientation inspection now traverses the typed stored expression. The BNL note documents the top-level `bench` syntax and explicitly records that its old fixed point and nonexistent command block cannot reproduce the archived timings. The late-parsing fixture and broader 23-file benchmark-suite restoration remain unchanged because they affect tests or introduce a larger tooling policy choice. |
| A68 | Diagnostic/production preparation difference accepted; nonblocking | `3Drep` is not required to reproduce GammaLoop's production numerator, initial-cut treatment, global-factor preparation, or display. | The command prints `Graph::full_numerator_atom()`, while GammaLoop production uses its separately prepared factorized numerator. That difference is intentional product-boundary separation, not an alignment task. No action is required unless shared CFF mathematics or a GammaLoop value is affected. |
| A69 | Diagnostic-only tooling regression; deferred and nonblocking | `3Drep graph-from-signatures` may be restored in a separate tooling task independently of LTD if desired. | Its former module, CLI surface, and tests are useful diagnostic source material, but the command is not required for GammaLoop generation, evaluation, or readiness. Any restoration should honor `minimize_externals`; it is not part of the production merge gate. |
| A70 | Shared-package maintenance implemented; auxiliary evaluator remains a manual regression oracle | Keep the package-local macOS `EXTRA_MACOS_LIBS_FOR_GNU_GCC` linker hook aligned with the repository workaround. Describe `eval` as an opt-in diagnostic oracle, document its serial command, and do not present it as GammaLoop production or normal workspace-test coverage. | `three-dimensional-reps/build.rs` now emits the standardized `gcc_s`, `gcc`, and `gfortran` links on macOS when requested. Its README documents `cargo nextest run -p three-dimensional-reps --features eval --no-fail-fast --test-threads 1` and correctly separates that manual inventory from GammaLoop evaluator construction. |
| A71 | Auxiliary display issues recorded; deferred unless a pre-rebase regression is demonstrated | Retain `loop_q0` as diagnostic/fallback and `edge_q0` as the authoritative numerator map. Possible auxiliary follow-ups are correcting orientation coloring and preserving deliberate generalized `0`/`x` labels across physical remapping. | The terminology comes from the pinned note. The current display can gray internal characters using an external-edge count and can overwrite a generalized base label, although `|N#` still keeps selection IDs distinct. No GammaLoop value issue or pre-rebase auxiliary-CLI regression has been demonstrated, so these output changes do not block the rebase and no test has been changed. |
| A72 | No-op public feature retained as compatibility vocabulary | Keep the default `serde=[]` feature name so downstream manifests do not break, while documenting that serialization dependencies and derives are unconditional. | No `cfg(feature = "serde")` exists, so toggling the feature changes neither build nor API behavior. The README now states this explicitly instead of implying that the no-op flag controls serialization support. |
| A73 | Protected GL11 denominator-only oracle isolated; snapshots unchanged | Earlier audit text mislabeled this graph as GL00. Preserve the protected GL11 `Q1^2` and `Q7^2` snapshots until their test intent is explicitly reclassified, but do not force generalized CFF to reproduce a denominator-only cancellation residue. | The stored values are `-1.17703539069413922e-13 - 2.37175489451059046e-30 i` and `-2.22459688841192310e-11 - 4.48261675062501648e-28 i`. The validated `cc8c5979` root path explicitly cleared `energy_degree_bounds` before generation and attached the numerator after residues. Current production deliberately retains physical EMR bounds and generalized contact maps, as required for `a*M` higher powers. Reproducing the old scalar by clearing bounds or arranging floating-point cancellation would undo that contract; no snapshot or fixture has been changed. Actual GL00 `Q7^2` is the distinct purely imaginary `-4.81280681117012732e-9` snapshot. |
| A74 | Maintainer-approved active-variable analyticity repair; EMR-only production boundary tightened | CFF capacity analysis needs two distinct facts: a nonnegative polynomial degree bound for UV/contact reconstruction, and validation that every energy-dependent factor assigned to CFF has no additional numerator poles. Negative powers of cut/bridge energies may remain in their separately projected external sector, but a negative power of an active CFF EMR affine form must fail rather than be assigned degree zero. | Production Graph analysis now uses `EnergyPowerAnalyzer::for_physical_emr_edges(active_edges)`. In one factorized traversal it rejects negative active powers, energy-dependent exponents, opaque nonlinear active-energy functions, and every `K(loop,index)` atom; a provenance-aware producer must normalize `K` to physical `Q(edge,index)` first. The legacy K-to-carrier constructors remain only for non-production analyzer diagnostics. This removes every LMB lookup from production CFF capacity analysis while leaving cut/bridge `Q` factors caller-owned. |
| A75 | Maintainer-approved product and diagnostic boundary | Use GammaLoop terminology and judge completion by GammaLoop generation, projection, evaluation, persistence, and user-facing commands. Audit the shared CFF library wherever GammaLoop depends on it. Treat the separate `3Drep` CLI and feature-gated evaluator as diagnostic tools, not production contracts; their inputs and prepared expressions need not match GammaLoop. | CLI-only defects and improvements, including A69, are optional tooling follow-ups and do not block this rebase. Shared-library mathematical defects remain in scope when they can affect GammaLoop, regardless of where focused tests live. |
| A76 | Latest-main evaluator/archive replay resolved, validated, and squashed | Preserve main #95's unconditional evaluator construction, removal of `override_if`/override input, `0441bd7a` Symbolica stack with native-code support, and O2 compilation behavior. Layer the approved deferred explicit CFF orientation sum through evaluator-local tagged function maps. Explicit-sum stacks evaluate once from representative input, reject every individual-orientation request, and do not resurrect main-deleted conditional tests. | The only semantic source conflict in `kqqlnrnx` was resolved by retaining main's ordinary orientation guard and adding the explicit-sum dispatch. The top replay retained only the two new deferred-CFF evaluator tests, updated `collect_orientation_if` to main's one-argument API, removed all stale override references, and combined archive layouts as amplitude version 6 / cross-section version 9. Symbolica references use main's `0441bd7a`; O2 forwarding/defaults remain. Formatting, locked GammaLoop library compilation, both deferred-sum tests run serially, the explicit-sum runtime/LTD boundary, main's ordinary orientation-Monte-Carlo guard, offline metadata, Hakari verification, and an independent generated-amplitude-loader compile pass. GammaLoop's amplitude `save standalone --json` and cross-section archive/loader export tests pass. A follow-up now deserializes the amplitude JSON into its typed archive, loads it, and performs a real finite eager evaluation; equivalent typed cross-section reload/evaluation coverage remains open. The validated resolution was squashed into `lktlsypu`; no conflict remains in the selected ancestry. |
| A77 | Empty-forest production identity and single root-capacity analysis implemented | Treat the empty UV forest as the original factorized production integrand in both local-UV routes. Apply `local_uv_cts_from_expanded_4d_integrands` only to proper, nonempty UV nodes. Provision the root's physical EMR capacity once in `production_cff_3d_expression_options`; do not pass the full numerator again merely to overwrite the same bounds. | Main's ordinary root and the validated pre-rebase `root_expanded_4d` both assemble the final empty forest from the production root expression. The approximation, final-integrand, and hedge-poset guards implement that identity consistently, while proper UV nodes retain exact finished-4D projection. A focused regression verifies that toggling the route leaves both the root local term and finalized integrand identical. The numerator remains factorized and caller-owned; only its EMR degree envelope enters root generation. This supersedes A35's earlier projected-root requirement. |
| A78 | Maintainer-corrected convergence boundary implemented | Keep the factorized EMR rank-envelope analysis required to tell generalized CFF the maximal numerator power in each energy channel. Remove the separate bespoke contour machinery, but reject source integrals whose existing all-cycle UV DOD analysis finds a nonconvergent loop-energy integration. | These remain different questions: the rank envelope provisions polynomial reconstruction and finite-pole sampling, while source convergence reuses the existing LMB/UV analysis. For each nonempty potentially failing cycle union, GammaLoop computes `DOD_E = DOD_4D - 3 L`, which replaces the four-dimensional `4 L` measure contribution by one power per energy integration, and requires `DOD_E < 0`. The check runs once on the original amplitude or cross-section source before CFF generation; it does not infer convergence from EMR caps or recreate a per-Taylor-term exact-source contour analyzer. |
| A79 | Maintainer-corrected source-backed dotted lifting, minimax numerator assignment, and canonical topology batching implemented | The typed `GS.den` wrapper carries the original `EdgeIndex` unchanged through the UV Taylor operator. After derivatives raise propagators, repeat every negative denominator power as occurrence-local exact edges. Disjoint-set contractions of absent original edges build the known cograph and UV source-minor scaffold without inferring a graph from signatures. The retained `source_edge` seeds and validates an occurrence's attachment, records cut and physical-map provenance, and groups repeated copies of one wrapper for deterministic serial subdivision. Canonical momentum signature up to sign, mass, topology domain, and multiplicity define the rewritten rational denominator channels and therefore the owner-independent exact residue topology on that scaffold. A compatible owner relabeling alone cannot change incidence, loop rank, or residue. The only rank solve checks which unique `+/-` source routing agrees modulo the opposite domain after the source scaffold is fixed; it never chooses endpoints. Retain the raw exact sign for physical-numerator mapping and retain owner IDs for cut-support unions, surfaces, parent mapping, and source-minor validation. Term projection splits only the completed expression's genuine outer Taylor sum when its addends remain distinct rational topologies. Nested numerator sums and positive typed `GS.den` factors remain factorized and uncancelled; generalized CFF owns any resulting pinches. Complete pure-external non-vacuum boundaries with explicit source-crown hedges, and reserve an explicit fixed-boundary payload for future OS two-point insertions such as `(m,0,0,0)`. Canonically relabel the finished owner-free rational topology with Graphica for deterministic cache keys. Analyze the still-factorized numerator in physical EMR variables and use one immutable quotient/remainder minimax plan for both exact CFF bounds and factor-local substitutions. Plan all genuinely outer Taylor terms first; equal canonical topologies then share a generator envelope obtained by maximizing each repeated algebraic channel's total degree and minimax-redistributing it, while non-repeated bounds use componentwise maxima. | This supersedes A45/A57/A63's owner-sensitive incidence distinction while preserving their source-backed scaffold and explicit-sum normalization, and removes every need for Kirchhoff, matroid, signature-to-graph reconstruction, or synthesized balance edges. Per-term numerator plans remain authoritative even when generation capacity is shared. `term_projection_keeps_typed_numerator_factors_uncancelled`, `exact_source_routes_a_hard_uv_row_modulo_the_soft_cograph_span`, `exact_source_owner_relabeling_preserves_residue_rank_and_cut_provenance`, and the real-Taylor `dod_one_triangle_keeps_one_uncancelled_exact_source_and_matches_lower_sectors` cover the new boundaries; the latter retains one common-denominator DOD1 source with owner multiplicities `(2,1,2)`, factorized positive typed denominator factors, and one canonical cache topology, and matches an explicitly reduced lower-sector oracle. Same-owner cubic, compatible distinct-owner relabeling, UV/cograph product, rank-deficient, powered-cancellation, canonical-order, non-vacuum crown, and analytic cubic-pole tests cover the remaining structural layers. The three-route scalar matrix compares localized local 3D, explicit local 3D, and projected local 4D over dotted, double-dotted, spectacles, and mirrored `Q/Q/-Q` cubic fixtures; its native 1000-bit Arb comparison passes in roughly 3.5 seconds. The LU `q^2-m^2` cancellation passes all three routes in roughly 2.4 seconds. Full design and sign proof: [`exact-powered-denominator-cff-lifting.md`](exact-powered-denominator-cff-lifting.md). |
| A80 | Maintainer-directed direct-local3D restoration implemented; final matrix rerun pending | Both direct orientation representations first write the complete/global CFF as `sum_k delta_k C_k`, where `k` is one complete generalized residue-map key, and then apply every UV Taylor operator as `sum_k delta_k T(C_k)`. The selector is opaque to `T`; `explicit_orientation_sum_only=true` performs the identical construction with every `delta_k` replaced by one. Physical edge directions are key metadata, not the partition identity. Exact-source reconstruction and derivative-created minimax EMR dispatch are exclusive to projected local4D. Per-key UV finiteness is required and profileable only for the selector-parametric localized direct local3D; selector-free explicit-sum direct local3D and projected local4D are summed representations. | This supersedes every earlier implication that either direct route is reverse-engineered from completed local-4D Taylor structures. The focused DOD0/1/2 regression and all six base scalar-matrix graph tests (GL00, GL02, GL04, GL08, GL09, and GL24), including their localized-direct3D per-key profiles, are selected by `test_gammaloop`. The pre-reversal 2026-08-31 `15/15` scalar matrix remains historical evidence; the current six-graph matrix must be rerun before final readiness. |

## Per-change review

| Change | Intended responsibility | Post-rebase review |
| --- | --- | --- |
| `kysvnqlq` | Pin the Symbolica stack and related workspace dependency resolution before feature changes. | Superseded on latest main. Its selected rebased instance was abandoned because the maintainer-approved policy is to keep main's Symbolica and its only remaining diff would have downgraded main #95's resolved `0441bd7a` stack to `46d690719`. The remaining fourteen selected changes now descend directly from `39561014`. |
| `xvxsluoq` | Introduce infrastructure shared by the split feature lines. | Resolved to feature-output ignore rules only. |
| `uksvkwzw` | Carry the alternate monolithic LTD implementation. | Excluded by maintainer decision; change and bookmark will be retired. |
| `yqmtxtlv` | Update tensor and symbolic support crates for the split implementation. | Resolved to one main-owned bound relaxation; targeted compile check passes. |
| `mwpprnsu` | Add shared CLI state output helpers. | Rebased conflict-free as narrow change `9bf4495b`; its later CLI consumers remain in their owning descendants. |
| `sttlspyn` | Add command-template placeholder support. | Rebased conflict-free as narrow change `e23f1790`; parameterized command execution follows in `wrzsxkku`. |
| `wrzsxkku` | Parameterize command blocks. | Resolved as the parameterized-command/state-compatibility change on current main. It now owns the approved A13 load-order comment and fresh-process regression; the regression passes. |
| `tmlszmts` | Update CLI generation and state plumbing. | Resolved and validated as process-aware inline graph import plumbing in change `tmlszmts`. The fixture and physical cut-count API were moved back from later split changes to the feature they test, using main's existing cut-selection helper. |
| `qontqzwm` | Add stability and inspect quality-of-life reporting. | Folded and rebased conflict-free onto current main after raised CFF. Bench/inspect/read-only JSON, stability, nonzero CLI failure status, and Python state-access behavior are retained; A13 was moved to `wrzsxkku`. |
| `yqnvuymp` | Add the three-dimensional-representations core crate. | Resolved as a CFF-only, locked-buildable core. The inert one-variant public mode was removed at this stage; `xxouqtmp` later adds the real `Cff`/`Ltd` selection boundary while keeping LTD explicitly unimplemented. |
| `kqqlnrnx` | Route raised-energy CFF through the 3D-representation implementation while leaving LTD explicitly deferred. | A50 restored source-EMR bounds and edge remapping, and A55 restored repeated-channel CFF normal form. Main graph/orientation/cut/UV and spinney-owned numerator growth remain required. The broader audit found the open EMR-sampling, certificate-domain, active-rank, occurrence-provenance, and map-identity defects in A56–A66; the current change is therefore not ready for final replay/push. |
| `nrzomklu` | Integrate 3D representations into the CLI. | Replaced by fresh conflict-free CFF-only change `tqsrsqrk`: Validate+Build, exact support, invocation-local sampling override, always-regenerated artifacts, CFF display, and one end-to-end regression. |
| `lrznwzuo` | Merge the QoL and 3D-representation feature lines. | Its empty structural merge is superseded by the selected linear stack, with the audited QoL change placed after raised CFF. No replay is needed. |
| `tormmyww` | Update example run cards for the merged behavior. | Preserved as separate source material and not yet landed; its NNLO card must be reviewed against the staged CFF/4D architecture rather than copied wholesale. |
| `wskyzxur` | Add repeated-mass CFF regression coverage. | Replaced by fresh conflict-free change `xxwzulrk`: two threshold-disabled CFF mass-limit tests, four new fixtures, original signs and strict convergence assertions, and no deferred threshold/localizer machinery. |
| `ltnoyvll` | Add LU scalar cross-section regression coverage. | Its two CFF routes and 198 live snapshot cases are ported in `lktlsypu`; semantic sign-off is blocked by GL11 and the proper-LTD third leg remains preserved and deferred. |
| `swnpzvnu` | Integrate the example and regression branches into the split head. | Superseded by fresh narrow linear changes. Its `GraphImportOptions` reconciliation belongs to `tmlszmts`; remaining example-card and proper-LTD work is tracked independently. |

## Validation record

This is a chronological record. A passing historical checkpoint does not
override the live readiness table or a later regression finding.

- Pre-rebase working copy: clean.
- Pre-rebase source topology: captured with `jj log`.
- Pre-rebase overlap preview: 67 files touched by both sides; approximately 34
  textual conflict files plus the generated workspace-hack manifest are expected.
- Rebase result: 16 retained changes were initially rebased onto `main`; 15
  contained conflicts. The selected fresh ancestry through `lktlsypu` is now
  conflict-free. Conflicted obsolete source revisions remain only as material
  for the not-yet-extracted example-card, QoL, localizer-test, and proper-LTD
  deltas.
- `wrzsxkku`: `cargo fmt --all -- --check`, `cargo check --offline -p
  gammaloop-api`, the three `run -D` completion tests, the late-parsing run-card
  test, and `cli_stateful_workflow_behaviors` pass. The API check reports only
  pre-existing dead-code warnings for state-output helpers awaiting their later
  consumers. The first test build exhausted the workspace quota; clearing only
  Cargo's reproducible test-profile artifacts freed 34 GiB, after which all
  focused tests passed.
- Approved model-first reload regression: implemented in the compatibility-owning
  `wrzsxkku` change. Its fresh-process check probes one-character `UFO::G` Typst
  output and the exact obsolete duplicate-symbol warning and passes.
- `tmlszmts`: `cargo fmt --all -- --check`, `cargo check --offline -p
  gammaloop-api`, the same check with `python_api`, triple-quoted parser tests,
  inline-import completion, and `cli_stateful_workflow_behaviors` pass. The
  stateful test verifies the original physical process-spec behavior, including
  reduction from 25 multi-edge candidates to 9 compatible selected cuts. The
  final escaped-quote/semicolon parser regression also passes. Checks report
  only pre-existing dead-code warnings for state-output helpers awaiting later
  consumers.
- Raised-CFF numerator-ownership audit: main's `Local3DCts::root` is
  denominator-only; `local_3d::start` and `t_tilde` grow by the exact
  `current - given` numerator fragment; final assembly adds the graph complement
  and global factors once. The historical raised-root override instead used a
  fully assembled `.expand().collect_factors()` numerator only at the root and
  is rejected by A14/A15.
- Raised-CFF map audit: `OrientationExpression` retains exact loop/edge energy
  maps and `numerator_map_index`, and `LinearEnergyExpr` retains internal,
  external, uniform-scale, and constant pieces. The pre-port `CFFTerm`
  conversion dropped that information; `smlxpnry` now keeps structured map
  identity and factor growth through the forest boundary for production paths.
- Exact-map implementation: `cargo check --offline -p gammalooprs --lib`,
  `cargo fmt --all -- --check`, `git diff --check`, the affine two-factor map
  test, the distinct-map-generation test, and focused exact-projection/root
  weighting tests pass. `scalars_integrated_cts_compare_legacy_and_hedge_poset`
  also passes. Clippy completes with known dead-code/argument-width warnings;
  two new `filter_map(bool::then)` warnings were removed.
- Exact-map review found and fixed three P1 issues: root branches were divided by
  their extension count, full patterns were applied to reduced undirected
  edges, and integrated finite CTs were left coarse. Exact mode now maps every
  owned factor once and does not run final coarse sign replacement.
- The zero-inner contracted-map lift is implemented and independently reviewed.
  Source coordinates solve `B^T c=r_e` in `[inner, outer]` order; reconstructible
  maps are the exact rational outer combination plus the original external
  signature. A focused source regression accepts a contracted tree edge and
  rejects an inner-dependent carrier. The approved finite-CT regression now
  uses a parallel loop plus a real external-flow bridge, proves the bridge is in
  `contract_edges`, and verifies both exact localization maps; it passes.
- A18 implementation: `EnergyPowerAnalyzer` remains the sole structural `Atom`
  walker; `normal_emr_replacement`, compatible inner LMB construction, exact
  rational row solving, and the existing affine maps remain the owning
  primitives. `EnergyPowerSupport` is the one shared coefficient-blind support
  type; `KnownPolynomial` and Symbolica polynomial materialization remain
  rejected as the boundary. Source-local factorized analysis, shifted affine
  carriers, and global correlated disconnected support have focused coverage.
- Sparse raised generation implements noncontiguous support interpolation,
  finite-pole remainder/contact recursion with independent authoritative EMR
  samples, and complete uniform-`M` normal forms. Numerical eval coverage compares
  sparse `{0,2,5}` and correlated shifted multiloop support against nonuniform
  construction at two `M` values and verifies `M` independence.
- The approved last-denominator boundary is implemented. A fully contracted
  scalar contact is a Unit lower sector; positive free-`q` powers fail with
  `UnlocalizedEnergyContact` before any partial contact is emitted. The live
  terminal quartic Unit/contact regression passes, but the historical dedicated
  None/All degree-two-through-four fixture was a one-edge self-loop and is not
  retained as a topology oracle. A balanced focused replacement remains open.
- The former A24 rank-projection conclusion is invalidated by A51. Its analytic
  oracle used `loops[]` on a fixture with no named carriers, so it tested the
  diagnostic `loop_q0` fallback rather than the authoritative bounded EMR map.
  Rewriting the same polynomial as `edges[]` and removing the false common-lift
  rejection makes generation complete, but the residue path returns
  `0.06064661065525115` instead of the analytic `0.05657291696421757`.
  The constant-contact parity case still passes. No oracle was updated.
- Historical convergence checkpoint: the report's nullity arithmetic and its
  ten focused module tests, including `d=2`, exponent overflow, and
  arrangement-mask overflow, passed before A52 made it a hard generation gate.
  A62 now shows that its single signature set wrongly serves as both denominator
  arrangement and numerator EMR domain; the arithmetic result is not meaningful
  until cuts, occurrence multiplicity, active rank, and numerator-only rows are
  separated. The previously cited penta bounds do not establish a false
  rejection.
- Compatibility cleanup removed the last `max_degrees()` collapse, all
  automatic edge-degree bridge APIs, and the unreachable affine/bounds and old
  finite-pole island: about 1.3k lines and all 15 core migration warnings.
  Comments on Cartesian correlation and contracted-contact signs were moved to
  their live owners. Formatting, affected package checks, sparse tests, the
  canonical clippy command, and an independent before/after capability audit
  pass.
- The inert one-value `RepresentationMode::Cff` scaffold is removed from the
  public options, display/eval, tests, and GammaLoop adapter. CFF display text is
  unchanged, and repository search finds no remaining mode field or enum.
- At the pre-A42 core checkpoint with `display,eval`, all 65 then-current tests
  passed. The two stale
  known-factor-origin assertions now verify sparse remainder/contact structure,
  and the then-named explicit-confluent assertion verified the fourth-order cut
  instead of a correctly empty same-routing sign map. The new A24 analytic, parity,
  rational-lattice, and initial-cut regressions pass. Warnings-denied all-target
  core Clippy passes.
- `cargo check --offline -p gammalooprs --all-targets --locked` passes with zero
  warnings after removal of the seven audited dead GammaLoop migration groups. The
  exact outer-source support regression, its zero-inner reconstruction checks,
  and `contracted_raised_generation_uses_exact_source_owned_support` pass. The
  full GammaLoop library run still reaches pre-existing red fixtures and aborts on
  the unrelated cross-section stack-overflow test; no such tests were edited.
- Core residue cleanup removed roughly 300 net lines of zero-caller cut-aware
  selection/sign plumbing while preserving canonical/generated-basis residue
  behavior, comments, and neutral powered-pole support provenance. Formatting and
  `git diff --check` pass on the combined tree.
- The generated Nix workspace graph now includes `three-dimensional-reps`;
  `gammaloop-guppy-workspace-graph` validates the regenerated JSON.
- Core-stage audit: an isolated `yqnvuymp` archive passes formatting, Hakari,
  locked all-feature check, clippy with warnings denied, test compilation, and
  all 42 unit/doc tests. The later `rpspmzvw` reserve is structurally invalid as
  a stage and must be rebuilt as described by A19.
- QoL-stage audit: the resolution is folded into `qontqzwm` and stacked after
  raised CFF on current main. Formatting, diff check, `cargo check --offline -p
  gammaloop-api`, all six benchmark unit tests, and the relocated A13
  fresh-process reload regression pass. The constant/plateau fitter, numeric UFO conjugation, and read-only
  generation guard are tracked separately pending maintainer confirmation;
  exact arithmetic/style corrections are queued as mechanical follow-ups.
- Main rebase checkpoint: `smlxpnry` was folded into clean `kqqlnrnx`; the
  reviewed QoL resolution was folded into `qontqzwm`; the ten selected changes
  were moved from `xpqplvru` onto current main `pzptorkt`, and the QoL line was
  stacked after raised CFF. The selected ancestry contains no conflicts. The
  old mixed `nrzomklu`/merge/test descendants remain source material and will
  be replaced by fresh narrow deltas rather than replayed.
- Main's `local_4d.rs` is the selected implementation for the 4D-local route;
  the alternate `expanded_4d.rs` and `Expanded4DApprox` implementation are
  absent. Main's direct-3D route remains the default and the typed 4D-local
  route is separately selectable in `ypnxsozo`. `rpspmzvw` remains preserved
  source material for the proper-LTD delta. `wskyzxur` supplies CFF-only
  coverage, while `lktlsypu` owns the two CFF scalar-battery routes and leaves
  their LTD comparison intact for the deferred LTD test change.
- Workspace quota recovery: after integration-test linking exhausted the
  filesystem, `cargo clean --profile dev-optim` removed 33.2 GiB of generated,
  reproducible artifacts and restored 11 GiB free without touching source or
  Jujutsu changes.
- The pre-battery raised-CFF scope audit found no scope blocker at that
  checkpoint. This statement is superseded for final readiness by the later
  GL11 findings. `expanded_4d.rs` is replaced by
  main's typed `local_4d.rs`, while its route intent and comments remain inputs
  to the separate 4D-local stage. The two root diagnostic Markdown files and
  all LTD-specific material remain assigned to deferred feature changes rather
  than being deleted. No isolated comment was removed from retained logic, and
  pending QoL does not leak into this stage. After the move to current main,
  `scalars.json` was restored to main's parameter order while retaining the
  required `goldstone` fields.
- `state_map.bin` was exported in a clean process against the current tree as a
  22,753-byte file with SHA-256
  `dca35a6485e22e528fc15b627ae2af496fa7f7a410a81e9fa21f7b1f2ac521be`.
  The existing `process::tests::failing::test_proc_definition_encode` fails
  after the successful export because it calls `State::import` on the same
  write-only `File::create` handle (`EBADF`). Per repository policy, the test
  was not changed.
- Conflict-free DAG check: `conflicts() & ancestors(lktlsypu)` is empty on the
  selected current-main path. Obsolete source descendants remain conflicted by
  design until their semantic replacements are complete and they are retired.
- Formatting: `cargo fmt --all -- --check` passes on the final staged tree.
- `cargo check`: `cargo check --offline -p gammalooprs --all-targets --locked`
  and `cargo check --offline -p gammaloop-api --all-targets` pass. The CFF CLI
  consumes the state-output guard; the remaining retained QoL output helper is
  explicitly reserved for later output-producing commands.
- Clippy: warnings-denied all-target core Clippy passes.
- Targeted tests: all 65 core tests pass; the final exact-localization finite-CT
  regression passes (1/1), as does the production-valid out-of-span initial-cut
  regression (1/1). The fresh-process A13 reload test passes (1/1), and the
  top-level GammaLoop benchmark tests pass (6/6). `git diff --check` passes.
- CFF-only CLI implementation: fresh change `tqsrsqrk` wires only the display
  feature, `3Drep validate`, `3Drep build`, exact support/generation, artifact
  output, and one new end-to-end box regression using main's existing
  `scalar_box.dot`. The test verifies graph validation, constant support, 14
  CFF orientations, and invocation-local `M_for_all` precedence over an
  unchanged global `none`. Formatting, Hakari, API all-target check, API-only
  warnings-denied Clippy, generated Nix graph equality, focused test, and diff
  checks pass. Full API Clippy currently stops in unchanged parent 3D-UV code
  on two eight-argument functions. Evaluate/LTD/reconstruction and mixed demo
  cards remain deferred.
- Repeated-mass implementation: fresh change `xxwzulrk` ports four new fixtures
  and the two threshold-disabled mass-limit tests with their original signs,
  three delta samples, monotonic absolute/relative convergence, and final
  `1e-3` bound. Both tests pass. Integration-test check, test-target
  warnings-denied Clippy, formatting, scope search, and diff checks pass. The
  source's unit-localizer/event/Arb test remains deferred.
- 4D bridge audit: residual powers must be extracted from exact term-level
  denominator occurrences, not graph metadata. Zero-active-loop occurrences
  remain typed-4D-owned residuals, loop-dependent occurrences remain active
  raised CFF, and factors sharing an edge ID but different full expressions are
  never grouped. No new public helper is required. This audit covers projected
  4D terms; it does not supply the generic shared-library CFF residual-denominator
  contract tracked in A49.
- Historical exact-normalization checkpoint, superseded by A55/A63/A79: the
  mass-jet implementation distinguished variant-local half-edge normalization
  from the ordinary global inverse-energy product and made four then-focused
  `exact_cff_` tests pass serially. Once A55 removed that CFF organization, the
  same-owner dotted and distinct-owner opposite-routing expected values ceased
  to be validated oracles. The tests are new to this rebase battery, and the
  broader audit led to the source-backed topology and derived explicit-sum
  normalization recorded in A79 before either production code or expectations
  changed. The
  earlier parallel run's global UFO initialization race remains environmental
  rather than semantic evidence.
- Typed-4D kernel validation: a clean scalar-bubble generation prints main's local
  `-i/(64 pi^3 E^3)` coefficient for both production extensions, with the graph
  symmetry factor outside it. The unchanged `dod0_bubble_uv` integration test
  passes: 1 run, 1 passed, 1703 skipped.
- Two-route 4D-stage checks: locked offline `gammalooprs` all-target check passes
  without warnings. The combined direct/typed local-3D suite passes 12/12,
  local-4D passes 2/2, direct-versus-typed union and spectacles coverage passes
  4/4, and restored dependency-frontier/triple-union coverage passes 2/2.
  Formatting and `git diff --check` pass. Strict warnings-denied all-target
  Clippy now reaches the two eight-argument `t_tilde`/`start` functions restored
  verbatim with main's direct-3D engine and reports main's existing
  `too_many_arguments` lint; no suppression or refactor was added in this
  feature change. The broad hedge module also retains main's stale inline
  `T`-versus-`Approx` snapshot in its explicitly named `failing` fixture, and
  its heterogeneous terminal test needs a larger test-thread stack; the latter
  passes unchanged with `RUST_MIN_STACK=67108864`. Neither test was edited.
- Proper-LTD source audit: `rpspmzvw` is a preserved 26k-line sibling source
  snapshot and remains unsuitable for wholesale replay. The interim live
  `RepeatedResidueBuilder` and rank-projection machinery were subsequently
  removed under A42/A44 rather than retained as a hidden second CFF/LTD engine.
  The remaining proper-LTD work is a new reviewed delta behind the reserved
  representation setting: LTD provenance/schema, numerical dispatch,
  comparison/evaluation surfaces, and LTD-specific coverage on top of the typed
  4D bridge.
- Deferred-LTD setting: fresh change `xxouqtmp` exposes serializable
  `RepresentationMode::{Cff, Ltd}` and `3Drep build --representation` (alias
  `--family`). CFF is the serde and CLI default; LTD returns the typed
  `GenerationError::NotImplemented` before graph selection or artifact writes.
  Production CFF options state their mode explicitly. The core error regression
  and the existing CFF CLI end-to-end test with its added LTD-error assertion
  both pass 1/1. Locked all-target/all-feature core check, API all-target check,
  warnings-denied core/API Clippy, formatting, and diff checks pass.
- Diagnostic cleanup: the undocumented `GAMMALOOP_DISABLE_CONFLUENT_CFF`
  escape hatch was removed from all three repeated-channel decisions in the 4D
  bridge. The serial exact-CFF suite passed 4/4 under the then-current mass-jet
  organization. A55 later removed that organization and restored the validated
  repeated-channel normal form within the single CFF implementation; A63
  reopens the new exact-source normalization assertions.
- Scalar-battery audit and implementation: source `b8bec534` added
  `local_uv_cts_from_expanded_4d_integrands` specifically to select between
  main's direct-3D UV recursion and a 4D-local forest before CFF projection.
  `ypnxsozo` now preserves both routes with the typed 4D implementation, and
  `lktlsypu` ports their CFF rich-output comparison. Its first raised-power run
  reaches A31; its LTD third leg and three-way helper remain intact in the
  deferred LTD test delta.
  A final evaluable 4D parametric integrand is a distinct feature and remains
  explicitly unsupported. The 223 source snapshots include 25 files no longer
  referenced by the current matrix; only the live 198-case inventory recorded
  above is ported. The maintainer approved changing the GL02 `Q1^4`, UV-on,
  threshold-on result from its historical floating residue to exact complex
  zero. A later source-history audit found 195 of 198 live snapshots
  byte-identical, not 197: GL02 without a numerator also has an unapproved
  cancellation-scale phase change, and GL11 `Q7^2` has an unapproved material
  change.
- Scalar-battery implementation compiles cleanly against the rebased APIs. The
  port excludes the 25 stale threshold/UV-flag variants and duplicates no
  logical snapshot key. Its initial default run exposed three independent
  raised-CFF defects in the intermediate implementation: inconsistent raised
  carrier selection, incomplete physical
  cut-group membership, and the `mink(4,n)`/`cind(n)` normalization regression.
  A31--A34 corrected those boundaries. GL24 then established A38: its
  transformed exact source has a valid source-local pole decomposition but no
  meaningful production-map identity. The final projected route therefore
  sums those source-local maps explicitly after the complete proper 4D forest,
  rather than weakening exact matching. Its earlier projected-root conclusion
  is superseded by A77's empty-forest production identity.
- A38 exact-source validation: the initial-cut physical-alias regression, the
  inactive-source support boundary, and the factorized affine/tensor mapper
  regression pass. The latter checks cancellation before inactive placeholders
  are zeroed, exact OSE replacement, external signs/shifts, uniform sampling
  scale, constants, abstract Minkowski components, and concrete spatial
  components without expanding the product.
- A38 projected-route validation: all six `typed_4d_projection` tests pass,
  including source-local pattern independence, per-body finite-CT mapping,
  spinney-owned numerator growth, external trees, and zero-loop pure trees.
  The exact-CFF suite passes 4/4; projected disconnected union/spectacles pass
  2/2; local-4D tensor preservation passes 2/2; deferred-carrier tests pass 2/2.
  The eight maintainer-approved existing tests were migrated without changing
  their physics/normalization intent. The pure-tree run exposed and fixed the
  mathematically mechanical rank of a zero-column loop-signature matrix. These
  are projected-4D residual tests, not substitutes for the absent generic
  shared-library CFF residual tests from the sibling source. They also did not
  contain the later GL11 distinct-physical/coincident-denominator group and
  therefore do not close A45.
- Deferred evaluator validation: both the scalar lowering/value regression and
  the closed bispinor metric-trace regression pass. The latter would previously
  leave an undefined tensor function inside the `FunctionMap`; it now follows
  the ordinary preprocessing path and evaluates to 4. Runtime/profile guards,
  evaluator-local tag/key alignment, and single-evaluation behavior in
  GammaLoop-generated standalone evaluators
  received an independent source audit with no remaining P0/P1 numerical or
  orientation-normalization finding. Generated amplitude and cross-section
  GammaLoop-generated standalone loaders compile against main's exact Symbolica revision.
- Saved-state review remains open: same-build serialization includes deferred
  bodies and explicit-sum metadata, but ordinary process/evaluator bincode
  layouts changed while `state_manifest.toml` remains version 1. A maintainer
  decision is pending between a narrow v2 rejection boundary and a full frozen
  v1 migration layer.
- Dual threshold-ID validation: locked offline all-target checks pass for core
  and API. The general cross-section info test verifies dense topology IDs and
  association resolution; the physical-cut exclusion test passes through
  generation and evaluation. Fourteen focused cross-section unit tests pass;
  the remaining large Symbolica helper test hits its known default-thread stack
  overflow before reaching the dual-ID code, and was not changed.
- Historical A42 validation, superseded by A55: the temporary mass-derivative
  organization had focused passing regressions for uniform/nonuniform samples,
  factorized numerators, orientation inheritance, and affine shifts. The theory
  and pre-rebase audit later showed that CFF must retain repeated-channel powers
  and perform channel normal form before lower-CFF generation. A55 implements
  that correction and keeps `RepeatedResidueBuilder`, `LowerSectorRankProjection`,
  raw `q±E`/H-surface construction, and their unused adapters absent.
- A43 normalization-ownership validation: generated expressions state either
  `GlobalSourceProduct` or `VariantLocal`; the GammaLoop adapter branches on this
  type rather than inferring ownership from empty half-edge lists. Focused
  mixed disconnected, pure-axis, sampled-numerator-orientation, and
  no-external-shift derivative tests pass.
- A45 pure-case topology validation: parser-expanded serial powers remove
  redundant edges and contract intermediate endpoints; exact parallel
  occurrence copies remove duplicate records without endpoint contraction.
  Both focused tests pass, including an order-independent undirected endpoint
  assertion, and disconnected recursion receives the topology setting. The
  subsequent GL11 interpretation that distinct physical owners create a third
  owner-sensitive exact-incidence category is withdrawn. Source owners still
  recover and validate the known minor, support same-wrapper power-chain
  bookkeeping, and survive in physical maps, but coincident rewritten rational
  channels have one owner-invariant exact residue topology. A neutral owner
  relabeling cannot change its incidence, loop rank, or residue.
- A46 live-path validation: isolated zero-edge vertices are removed before CFF
  recursion and therefore contribute the multiplicative identity rather than a
  `1/0` causal surface. The repeated bubble fixture is a genuine serial chain.
  The formerly private lower-sector cases now run through public generation and
  retain their nonintegral projection, initial-cut, energy-dependent numerator,
  and raised-LU assertions. Focused tests pass.
- A48 NLO regression diagnosis and validation: the reported GL0 support is the
  current implementation's LMB-coordinate sector `K0^0..3 x K1^0..1`; its
  cubic term was not the immediate cause of the reported failure. The global
  cross-section caller passed no contracted edges after the
  raised-CFF integration, unlike both current main and retained pre-rebase
  source. External-only tree edges 4 and 6 therefore entered the active CFF
  graph with zero loop signatures and reached the generic lower-sector guard.
  Restoring `tree_edges - initial_state_cut` contraction preserves the existing
  `tree_denoms` ownership boundary. `cargo fmt --check` and `cargo check -p
  gammalooprs --lib` pass. A release `gammaloop-api` binary was built, and the
  maintainer's exact `epem_a_ddx_NLO.toml` command then completed with GL0 and
  GL2 compiled, five evaluators each, saved state, and exit status 0.
- A50 semantic audit invalidates the raised-energy conclusion previously drawn
  from that successful command. Before `58a1634e`, numerator energy bounds were
  indexed by EMR/source-edge energy IDs and remapped into every local source.
  That commit instead passed sparse support indexed by source-LMB loop variables
  unchanged. The historical `kqqlnrnx` edit now restores edge-indexed bounds,
  remapping, factorized max/add/power analysis, partial supplied numerators, and
  strict failure rather than bound filtering or carrier reassignment for an
  unrepresentable shrunken edge. A declared loop-energy symbol is normalized to
  the matching explicit LMB carrier edge and still produces an EMR-edge bound.
  The adapter compiles and its focused tests pass individually; this is not yet
  an end-to-end raised-CFF validation.
- A51 theory audit: the exact source reviewed is
  [`generalised_ltd.tex` at repository revision `65ab03e4`](https://github.com/alphal00p/generalised_ltd/blob/65ab03e4fb7d442ff362392df2c2a59ef323d28c/docs/generalised_ltd.tex)
  (TeX blob `39bb90c37e295ebffd6f38780dd4a470097daaba`). The collapsed
  repeated-channel equation is a denominator identity only: the black-box
  numerator retains one argument and one degree cap for every original EMR
  edge. This independently confirms A50 and the restored edge-remapping
  boundary. The factorized GammaLoop envelope is compatible with that black-box cap
  contract, but the paper does not define GammaLoop's 4D UV forest or spinney-owned
  incremental numerator growth and therefore cannot replace A14's ownership
  requirements.
- A51 repeated/contact finding: physical repeated-copy sign samples, including
  mixed signs, need not integrate to one global loop-energy assignment.
  `edge_q0` is authoritative; `loop_q0` is a diagnostic particular solution and
  auxiliary data used by the feature-gated diagnostic evaluator. Consequently `LowerSectorRankProjection` must not
  recompute the edge map from one loop lift and reject a valid sample as
  `SingularBasis`. Numerator-call/orientation identity is the unique EMR edge
  map; variants may share that call even when auxiliary loop diagnostics differ.
  An auxiliary `3Drep` graph with no named carrier currently evaluates the diagnostic
  evaluator's `loops[]` syntax through `loop_energy_map`, but that legacy
  behavior must not be used to infer an EMR bound or split GammaLoop numerator
  identities. The corresponding shared-core call-deduplication proposal is
  retained under A64 but deliberately not implemented during this rebase.
- A51 contact-map finding: a contracted lower CFF minor supplies a lifted
  affine value `qhat_e` for known polynomial factors, whereas the deleted
  edge's black-box sample remains an independent finite-pole coordinate such as
  `+omega_e`, `0`, or `-omega_e`. Future pre-shrunken bounded-edge support must
  carry both concepts plus original source identity; an affine lift or a
  surviving-denominator carrier alone is insufficient. The current strict
  rejection is therefore an honest incomplete boundary.
- A51 bounded-state and branch finding: absence of bounds means the legacy
  affine black-box class, while omitted edges in an explicitly bounded vector
  normalize to degree zero. That mode bit must survive lower-sector and
  component projection even when every remaining cap is zero. The historical
  edit now represents this exactly as `Option<Vec<(edge_id, degree)>>`; the
  component regression distinguishes `None` from `Some([])` and passes. This is
  serde-visible: an old serialized `[]` now means an explicitly constant class,
  while an omitted/null field means legacy affine support.
- A51 dispatcher correction after comparing prose with executable reference:
  multiaffinity is per original EMR argument, not per collapsed denominator
  channel. The current Python pure-CFF dispatcher keeps `[1,1,1]`, `[2,1,0]`,
  and `[2,2,0]` repeated-copy caps in ordinary/quadratic pure CFF; it selects
  channel/known-factor normal form for an individual cap above two or the
  requested uniform mode. Inside that builder the interpolation degree is the
  summed remaining channel degree `d_C = sum(d_e)`. The earlier ledger claim
  that a summed degree above two must itself pre-empt per-edge quadratic
  recursion was incorrect and is withdrawn.
- A51 architecture and test finding: pure CFF keeps repeated denominator powers
  as E-surface factors and performs channel/contact normal form; the raw
  repeated-pole residue/finite-difference construction belongs to hybrid/LTD
  organization. This reinforces the approved A42/A44 removal instead of
  establishing a second “confluent CFF.” GammaLoop currently selects that second
  engine whenever a bounded production edge belongs to a repeated group. The
  approved test rewrite now uses `edges[]`; after removing the invalid
  global-lift rejection, its energy-dependent case evaluates to
  `0.06064661065525115` instead of the analytic `0.05657291696421757`. The
  constant-contact case passes, as does the isolated named-carrier versus
  no-carrier fallback test. No expected value has been changed: the numerical
  failure is evidence to restore pure-CFF dispatch, not to repair the temporary
  residue engine.
- A51 historical UV-policy discrepancy, superseded by A78: the reference note makes the loop-energy UV check
  a hard precondition and refuses representation construction on failure. For
  a scaling subspace with active edges `A`, it tests the cap envelope
  `delta = sum_{e in A}(d_e) - 2|A|`; one-dimensional contour closure requires
  `delta < -1`. The Rust nullity correction is independently sound for absolute
  radial scaling in an `r`-dimensional subspace: require `delta < -r`. A62
  implemented that correction before A78 removed this bespoke certificate and
  replaced it with the source-level `DOD_E = DOD_4D - 3 L` gate.
- A51 historical UV implementation audit, superseded by A78, corrected A20's wording: generation did not
  currently consult or serialize `energy_divergence_report`, so “advisory” means
  unchecked, not “build with an explicit uncertified warning.” A failed cap
  envelope is inconclusive for GammaLoop's particular factorized numerator because
  coefficients, affine correlations, and cross-term cancellations are lost,
  but it is a sound refusal for the shared CFF builder's advertised universal
  black-box class: that class contains cap-saturating numerators with a residue
  at infinity. The correct boundary may therefore be a hard shared-library default
  plus a typed external proof/unchecked-sum contract for GammaLoop, not a silent
  global bypass. The ledger's earlier “formerly blocked penta contact” evidence
  is also withdrawn: the live penta bounds are certified by both the
  one-dimensional and nullity-aware criteria and do not reproduce that claim.
  A78 resolves the formerly proposed hard-versus-advisory boundary by retaining
  the distinct EMR rank-envelope contract and making the existing source-level
  energy DOD a hard precondition.
- A51 post-generation-profile audit, retained as a distinction rather than a
  production contour requirement: GammaLoop's existing `UVProfileable::profile`
  requires a built `ProcessIntegrand` and samples simultaneous large *spatial*
  loop-momentum scalings of that emitted 3D representation. Its amplitude-only
  analytic mode likewise starts from `all_mighty_integrand` after CFF
  construction; analytic cross-section profiling is not implemented. It can
  therefore catch bad 3D UV scaling and bad local-counterterm cancellation, but
  it cannot recover or detect a complex loop-energy residue at infinity that
  CFF generation has already discarded. A78 replaces the separate contour
  certificate with the existing source-level UV-DOD gate rather than asking
  this spatial profile to replace it. Automatically running the latter is also much
  heavier than the cap check (evaluator construction, kinematics, multiple
  scale points, and nonempty loop-momentum subsets), and it is unavailable to
  the shared representation crate without GammaLoop runtime data.
- The approved narrow regression reuses
  `tests/resources/run_cards/uv/epem_a_ddx_xs_nlo.toml` and adds no helper,
  fixture, or snapshot. It unsets that older card's stale hard-coded production
  orientation pattern only in the new test, leaving the shared card unchanged.
  The focused run passes with GL0 generated, five evaluators, compilation
  disabled, and one test passed.
- Generic residual-denominator boundary: later source history contains
  external-tree and pure-tree regressions plus a typed `ResidualDenominator`
  contract. That is reusable CFF source material, not proper-LTD code. It is
  unnecessary for A48's catalogue caller, but generic auxiliary `3Drep` evaluated CFF
  output cannot contract such tree factors until explicit residual ownership
  is restored. No part of that non-mechanical port has been applied.
- A52 provisional exact-source boundary: the maintainer approved failing rather
  than guessing whenever nonconstant physical-EMR numerator bounds cannot be
  represented in an exact/shrunken occurrence source, while warning that this
  may be more conservative than necessary. Exact 4D generation now analyzes
  the factorized numerator actually owned by the term. An energy-constant term
  proceeds as the explicit class `Some([])`; a nonempty bound set fails with a
  message that the numerator-only affine EMR-to-occurrence contract is missing.
  Denominator occurrence IDs, loop carriers, and `source_edge` provenance are
  never substituted for that absent contract. This representability failure is
  independent of A51's contour-at-infinity condition: passing the latter does
  not construct the missing affine map. Revisit the strict boundary when that
  map is implemented rather than treating the current failure as a theorem that
  the 4D projection is impossible.
- A52 historical acceptance-route audit: at that checkpoint the maintainer's
  exact `epem_a_ddx_NLO.toml` command selected `final_integrand = ThreeD` and
  `local_uv_cts_from_expanded_4d_integrands = false`. It constructed full 4D
  local counterterms, but local subtraction used the direct 3D CFF path; it did
  not call `project_4d`, `cff_from_4d_denominators`, or the provisional exact
  source guard. The command did exercise the restored EMR-edge bounds, strict
  contour certificate, contracted-tree-edge handling, and incremental
  spinney-owned numerator growth. After a fresh fat-LTO release build, the exact
  command completed with exit status zero on 2026-08-21 and was rerun
  successfully on 2026-08-24 after A49/A74: GL0 and GL2 each built five
  evaluators, both compiled, the semi-inclusive integrand was duplicated, and
  state was saved. Configuration inspection and the earlier debug run were not
  counted as either pass.
- A52 historical direct reduced-CFF capacity audit: expression ownership was
  main-compatible. Root CFF has no eagerly assembled numerator; each UV
  transition maps only its newly owned numerator fragment; final assembly
  multiplies the outside/current-spinney numerator and global factor exactly
  once. At that checkpoint, a nonzero integrated counterterm generated its
  reduced CFF with the complete production-numerator bounds twice: once through
  stored production options and again through full-numerator reanalysis. The
  actual mapped product can instead contain the integrated finite counterterm
  together with spinney-owned and outside/global factors. The then-current
  capacity could therefore over-reject a bound on a contracted edge or
  underbound new outer-energy dependence introduced by the integrated
  counterterm. The
  then-live acceptance card had `generate_integrated = false` and returned
  before this path, so that checkpoint could not validate or invalidate the
  pending non-mechanical fix. `opmnpruw` later set the live card to
  `generate_integrated = true`, retained exact source-local EMR maps through
  finite-UV localization, and unblocked integrated generation; the A52 result
  is historical route evidence, not current numerical acceptance.
  A53 later resolved the graph-theoretic bridge subset by retaining main's
  separately projected `tree_denoms` ownership, without relaxing contracted
  non-tree numerator capacity.
- A55 compile/lint checkpoint before A58/A62: `cargo check -p
  three-dimensional-reps --tests`, warnings-denied core Clippy, formatting, and
  `git diff --check` pass. The restored channel-normal-form builder generates
  the former degree-five sunrise blocker. The evaluator-enabled serial nextest
  inventory runs in isolated processes and gives 46/53; A59 classifies all
  seven failures, and no assertion, expected value, or test fixture has been
  changed. A single-process libtest run still encounters the global Symbolica
  initialization race and is not used as semantic evidence.
- A52 projected-route boundary, superseded by A61: the GL11 `Q7^2` comparison
  stopped at the approved A52 exact-source guard before CFF generation because
  its nonconstant physical-EMR numerator capacity had no numerator-only
  occurrence map. A later forensic audit found that the still-earlier exact
  zero arose after healthy generation when massless exact surfaces compared
  normalized `0` against raw `UFO::ZERO^2` and were all marked nonphysical.
  Current normalized mass handling fixes that separate issue; A61 records the
  physical-to-occurrence capacity boundary.
- A61's pre-rebase comparison proved that GL11's A52 stop was a GammaLoop
  regression, not a prospective auxiliary feature. The old expanded-4D adapter
  copied physical EMR degree caps to every exact denominator with the same
  owner. That was sufficient for GL11, whose `Q7` has the single literal exact
  occurrence 15, but it was not a general affine contract and overcounted
  repeated owner copies. The runtime mapper already evaluates the factorized Q7
  numerator from occurrence 15. Restoring only the unambiguous
  capacity translation is therefore separable from the larger numerator-row
  API; ownership-only, shifted/UV, and genuinely affine cases must continue to
  fail until their energy identity is explicit.
- Independent direct-route evidence, later reclassified by A73: GL11 `Q7^2`
  has the same pre-rebase value in the shared-core and integrated GammaLoop
  routes already on Symbolica 2.2. Forensics then showed that both validated
  roots deliberately cleared numerator-energy bounds before residues. Their
  agreement therefore protects the old denominator-only representation, not a
  generalized raised-EMR value. The unchanged snapshots still require an
  explicit test-intent decision.
- Cancellation-noise audit: GL02 without a numerator agrees between both live
  routes but is the exact negative of the stored value at magnitude at most
  about `4.2e-16`. Selected GL16 no-numerator, `Q1^2`, and `Q1^4` cases also
  agree between routes while drifting around zero, at most about `2.31e-15` in
  the inspected values. A component clamp at `1e-14` would affect 188 of 198
  live snapshots, so no global clamp or snapshot rewrite is justified.
- Current A49/A53/A56/A58/A62/A74 checkpoint: structural non-cut bridges are
  removed from the reduced production CFF denominator/capacity domain and
  retain main's caller-owned `tree_denoms` projection. Standalone `3Drep build`
  now supplies the same Graph-owned bridge set explicitly and receives lifted
  affine maps plus typed residual denominators instead. The contour report
  receives denominator occurrences and physical numerator EMR bounds as
  distinct `D` and `N` inputs. That certificate is valid inside `span(D)`;
  A58 now independently requires a top-level public CFF source to span every
  declared LMB direction, while lower-sector component lifts remain
  rank-deficient only in their proven context. Formatting, core checks,
  `gammalooprs` library check, warnings-denied shared-core Clippy, two
  focused certificate regressions, the finite-counterterm bridge-map
  regression, strict contracted-nonbridge rejection, the two restored
  auxiliary `3Drep` tree regressions, and numerical residual-factor reconstruction
  pass. A56's production-path regression proves that finite-pole sampling
  preserves the diagnostic loop map and complete spectator EMR map; focused
  quadratic remainder/contact, quartic Unit/contact, and mixed known-factor
  tests also pass. All six focused typed-4D projection tests pass after A74 made
  the active-domain atom traversal strict; the formerly failing assertion was
  not changed. Warnings-denied `gammalooprs` Clippy remains blocked only by the
  separately existing `too_many_arguments` warnings in
  `local_3d::{t_tilde,start}`. The refreshed evaluator-enabled serial inventory
  remains 53/59 after A56: the approved A58 replacement contributes two passing
  tests and removes the invalid cut/ownership failure; the six remaining
  failures are classified by A59. After a fresh release CLI rebuild, the
  maintainer's exact `epem_a_ddx` command completes with A56 and builds ten
  evaluators across GL0 and GL2.
- Current DAG/push checkpoint: the ledger is included in `lktlsypu`, and the
  concise cross-machine handoff is its documentation-only child on the
  published `raised_energy_cff_wip` bookmark. The fourteen retained feature
  changes have been replayed directly onto local/remote `main` at `39561014`;
  the obsolete selected dependency-pin root was abandoned. The
  evaluator/archive conflicts with main #95 are resolved, validated, and
  squashed into `lktlsypu`. `proper_ltd_core` remains a fully pushed sibling
  source reserve, not a descendant. The unrelated conflicted bookmark
  `push-wwvzyzvysvro` remains outside the selected ancestry.

## Broad pre-rebase and theory audit after A55

The comparison baselines are the validated pre-rebase shared CFF core `b9c4e2ff` /
`8a9dcea`, the later pre-`58a1634e` snapshot `0c8d166b`, the integrated source
snapshot `cc8c5979`, and
[`generalised_ltd.tex` at revision `65ab03e4`](https://github.com/alphal00p/generalised_ltd/blob/65ab03e4fb7d442ff362392df2c2a59ef323d28c/docs/generalised_ltd.tex).
The audit deliberately separates a text-conflict mistake from a latent source
defect that becomes observable only after the rebase, and from a new test whose
expected value has no historical oracle.

- A56 was a rebase-chain regression already present by `0c8d166b`, independent
  of `58a1634e`. Finite-pole division now again replaces only the selected
  physical EMR value: the carrier-rewrite helper and its passing wrong test are
  gone, and the production regression verifies the complete spectator map and
  diagnostic loop map against the degree-reduced source. Repeated-channel
  uniform sampling remains different and correctly changes all signed aliases
  of the one channel variable together.
- A62's certificate rank and A58's source-rank boundary are now repaired as two
  separate invariants. Initial cuts are inactive external-energy aliases, not
  constraints that own or fix a free loop-energy direction. The certificate
  receives denominator occurrences and physical numerator variables as
  separate `D` and `N` inputs and checks scaling subspaces inside the active
  denominator span. Independently, public top-level pure CFF now requires
  `rank(D) == loop_names.len()` so no declared integration coordinate is
  silently set to zero; only the lower-sector component lift may choose a
  deterministic particular solution after projecting to its proven rank. The
  generic report defect predated the rebase; the approved A52 hard gate made it
  a production acceptance bug. A53 further classifies non-cut `tree_edges`
  structurally as bridges projected outside CFF: they belong to neither CFF
  denominator suppression nor active loop-energy numerator degree, although
  their complete external-affine physical map and numerator factor remain in
  the production caller-owned projection or A49's auxiliary-CLI typed-residual
  result.
- A57/A63 established same-wrapper serial-power bookkeeping and the
  contour-correct explicit double-pole sum `+i/(32 pi^3 omega^3)`. A79 sharpens
  the provenance boundary: source owners seed and validate the known minor,
  retain cut support and physical maps, and identify copies of one wrapper for
  deterministic subdivision, but they are not variables of the rewritten
  rational function. Canonical signature, mass, domain, and multiplicity fix
  the owner-invariant exact residue topology, so a compatible owner relabeling
  cannot change incidence, loop rank, or residue. A79's sign-safe minimax
  numerator plan and concrete fixtures are documented in
  `exact-powered-denominator-cff-lifting.md`.
- A64 is a pre-existing note/code discrepancy, not a demonstrated rebase value
  regression. The note defines the authoritative production numerator identity
  by `edge_q0`, while labels and grouping also include diagnostic `loop_q0`.
  GammaLoop production already evaluates numerator energies from physical edge
  maps. The retained future design groups by canonical edge maps and keeps one
  deterministic loop diagnostic, but the maintainer deferred it because the
  shared-core shape change could affect the eager-f64 evaluator's orphan
  `loops[]` fallback. No GammaLoop-only duplicate layer is planned. LTD remains a
  selectable name that immediately reports the approved unsupported-mode error.
- A65 finds two uniform-scale boundaries: GammaLoop accepts `M=0` even though the
  generated expression can divide by powers of `M`, and the `All` policy
  activates for degree one although the note only needs repeated-channel
  `d_C > 1`. The first is a correctness defect; the second is unnecessary
  expression growth but not a nonzero-`M` value error.
- A66 records the missing validation evidence rather than claiming deleted Rust
  goldens. The external generalized-LTD reference validated unique edge maps, numerical
  quadratic/multiloop/repeated cases, 210 UV-convergent box caps, and uniform
  scale invariance; the Rust suite is mostly structural. Two present analytic
  tests also became unreachable behind A52. No test has been changed during
  this audit.
- A73 corrects the earlier graph label to GL11 and reclassifies its protected
  `Q7^2` scalar. The old value is stable because the validated root path cleared
  numerator-energy bounds before residues; it is not a value oracle for the
  generalized EMR/contact representation. The protected `Q1^2` and `Q7^2`
  snapshots remain unchanged pending an explicit test-intent decision.
- The GammaLoop ownership review found no second numerator-growth regression.
  Direct 3D still attaches only each spinney-owned fragment. The projected route
  builds the complete typed 4D term and maps `term numerator × outside numerator
  × global factor` once; final assembly does not regrow it. Exact occurrence IDs
  own the owner-invariant rational residue, while the separate physical-owner
  maps retain cut support, parent projection, and source-minor provenance.
  A63's explicit-sum normalization oracle is unchanged by compatible owner
  relabeling. Canonical node labels are assigned only after the source-backed
  scaffold, power subdivisions, and source-crown hedges are complete.
- A60/A67/A68 cover non-CFF integration maintenance and product-boundary
  findings: unversioned integration workspace bincode after the stability field
  addition; stale committed
  schemas; one obsolete Symbolica example pin; missing read-only write guards
  and dropped default-output QoL; unported live benchmark tooling/docs; and a
  `3Drep` diagnostic that prints a different numerator from production. The
  schema regeneration, exact dependency pin, and guards are mechanical. The
  diagnostic/production preparation difference is intentional and requires no
  alignment. Default-output, benchmark API, and every affected test remain
  approval-gated.
- 2026-08-26 production-bridge audit: integrated finite UV terms retain their
  exact source-local EMR maps; production selector IDs identify complete
  generalized residue-map entries without replacing source algebra. Physical
  direction data remains metadata inside those entries. The shared connected-loop and pure
  duplicate-denominator global sign is typed metadata consumed exactly once for
  root, reduced, and exact GammaLoop CFF sources. At that checkpoint, focused
  coverage and integrated generation established those ownership boundaries,
  while inclusive `epem_a_ddx` NLO normalization had not yet been rerun; the
  following validation entries supersede that historical gap.
- 2026-08-27 inclusive validation supersedes that open status. DDx passes the
  `alpha_s/pi` normalization, three-route pointwise and integrated comparisons,
  per-graph `M_UV`/localization-scale invariance, and opposite GL0/GL2 `mu_r`
  slopes. Fully MSbar TTX passes all three routes after applying the published
  direct-photon-to-`e+e-` conversion. The current default UV profiler selects
  one physical-LMB representation per expected divergent spinney on the same
  amplitude/LU domains as generation, retries poor full-range fits in Arb, and
  requires a long high-scale suffix with `R^2 >= 0.99`.
- The same day's tentative GL24 blocker is withdrawn. Its fixture used a
  nonlocal global numerator, and the handwritten one-dimensional contour
  calculation used to challenge the surviving cut is not an authoritative
  higher-loop oracle. GL24 must instead use a Feynman-rule-local numerator and
  be validated across the three production routes, with the generalized
  `3Drep` implementation and CLI as the independent 3D diagnostic lane.
- 2026-08-31 validation supersedes the remaining tentative GL24 wording. The
  replacement scalar LU matrix generates its graphs from the scalar model;
  unit-numerator processes are not mutated, and companion numerators are
  attached through local edge ownership. The 15-case all-enabled matrix covers
  orientation-local 3D, explicit-sum 3D, and projected local 4D at f64/native
  Arb with local UV, integrated UV, and threshold counterterms enabled. The
  final `dev-optim` / `test_gammaloop` rerun passed `15/15`, with `235` skipped,
  in `70.438 s`. Four near-zero cases use the authorized f64-input `1e-14`
  unit-scale fallback, while Arb-to-Arb remains run and reports non-scaling;
  this is test-oracle handling only and required no production change. The
  retained 10x physical campaign is independently `4/4` green: direct and
  converted DD/TT LO and NLO
  results reproduce their published absolute targets within `1.7 sigma`, and
  the `alpha_s/pi` or paper-ratio pulls are respectively `+1.141`, `+1.453`,
  `+1.037`, and `+0.685 sigma`. The exact values, errors, graph rows, and pulls
  are recorded in the readiness table; converted DD graph rows have no
  separately published targets. UV-profile CLI coverage is `2/2` (observed
  1.213 seconds) and API/unit coverage is `15/15` (observed 0.660 seconds);
  divergent-only selection uses the
  generation LMB when suitable and otherwise the first suitable
  deterministically sorted basis. The full shared CFF crate is `98/100`; its
  remaining failures are the stale contracts
  `lower_sector_powered_pole_contact_reconstructs_numerator_derivatives` and
  `inherited_mixed_theta_basis_uses_correlated_parent_residues`. The approved
  rank-projected origin assertion has been updated and now passes. The published
  remote checkpoint for this status is `raised_energy_cff_wip` at `5c3e7f`;
  any newer local repair remains subject to final-suite validation before push.
- 2026-09-01 architecture restoration supersedes only the live readiness claim,
  not the historical measurements above. Both direct-local3D forms now apply UV
  Taylor operators after complete/global CFF energy integration and differ only
  by orientation-selector omission; exact-source reconstruction and minimax EMR
  dispatch remain projected-local4D-only. The focused DOD0/1/2 bubble regression
  and all six base scalar-matrix graph tests, including their localized-direct3D
  per-orientation profiles, are selected by `test_gammaloop`. Their current
  post-restoration rerun remains pending.
