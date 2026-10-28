= Combined generation integration plan
<generation-integration-plan>
#quote(block: true)[
#strong[Status:] accepted for execution on 2026-10-27. Sol implements each
bounded step; parent Astra assigns independent Astra review after each step.
WP0 evidence preservation and source availability was approved by Astra round 2
(P2 closed, no remaining findings). WP1 corrected certificate `f473709c` and WP2 correction `5226894e` were
approved by Astra round 2. B2 is closed for the original GL028/GL034 pair,
not all grouping physics. WP3 implementation is authorized. Routine development and validation within accepted
packages need no repeated approval. Publication and destructive historical cleanup retain their
explicit gates. The earlier stopped v5 rebuild remains historical; no rebuild
is performed in WP0.

#strong[Scope:] turning the combined Feynman-generation line (the reviewed
v4 stack) into a publishable FeynKit integration on current upstream, then
completing the agreed API stages (names, `route()`, Classes gate), closing the
correctness blockers and the performance and memory regressions, and
recording milestones with gammabench. The drawing and licence lane is
tracked here only where it can block publication.

#strong[Evidence labels:] #emph[reported] means a historical claim from the
handoff or the 2026-10-22 review; #emph[observed] means read in this
repository on the date stated; #emph[future] means work this plan proposes.
WP0 checks availability and archive integrity only. Historical figures remain
reported, not reproduced. See #link("generation-integration-wp0.typ")[WP0 execution record]
and its evidence inventory for surviving receipts and inaccessible raw evidence.
]

== Sources and precedence
<sources>

All handoff paths are under
`/common/dev/gammaloop/handoffs/linnet-generation-combined/`.

+ The 2026-10-27 execution authorization and workspace/base instructions take
  precedence over older draft and source-repository restrictions. User decisions
  in `HANDOFF.md` section 9 remain binding where not superseded.
+ The ground rules and decisions in `HANDOFF.md` sections 1 and 6, and the
  plan in section 4 where section 9 does not override it.
+ `decisions/names-routed-classes.md` for the design of names, `route()` and
  the Classes gate.
+ `review/plan-review-2026-10-22.md`: findings B1–B2, M1–M11, m1–m12 and
  questions U-1…U-11. Its recommendations are inputs, not decisions; where
  section 9 answered a question, the answer wins.
+ Receipts in `prepare/`, especially `root-integration-b-quotient.json`.

For factual implementation status, inspected repository state takes precedence
over an older report; it does not override a user decision. Before each work
package, re-read the refs it depends on instead of trusting the figures below.

== Known baselines
<baselines>

#table(
  columns: 3,
  [*Item*], [*Value*], [*Basis*],
  [Upstream FeynKit], [`3f0df2a8`, including merged \#209 and \#210], [current accepted base, observed 2026-10-27],
  [\#210], [merged colour-only coefficient completion; supersedes \#177], [observed commit; HANDOFF §9],
  [Combined v4], [`kpslpwqz e0a77980`, 46 changes on `58814885`, unpushed, Codex repo], [Git objects imported read-only and pinned in WP0; historical tests not rerun],
  [Stage-1 merge / transplant / v3], [`f25961e2` / `4e008517` / `dfe63876`], [reported],
  [Tensor/MathML line], [merged upstream through \#209 at `3f0df2a8`; admission \#193 comes through that merge], [observed 2026-10-27; no separate dependency choice],
  [SM three-loop filtered snapshot], [3150 diagrams], [reported correct (U-2 classification)],
  [`four_loop_numerator`], [4664 / 306], [reported, at 1/8/32 threads on v4],
  [Oracle fast tier], [364/364 at 1, 8, 32 threads], [reported at v4],
  [Full-mode Ir, v4 ÷ stage 1], [gggg_L2 1.748, uugg_L2 2.023, ee4f_L2 1.367, eemm_L2 1.913], [reported (callgrind, 3 runs)],
  [Full-mode peak (DHAT), stage 1 → v4], [gggg 163 → 365 MB, eemm 687 → 1696 MB, ee4f 116 → 249 MB], [reported],
  [Wall, qcd_gggg_L2 full, v4 ÷ stage 1], [2.17 [1.50, 2.45] at 1 thread; 1.56 [1.12, 2.21] at 8], [reported, 11 rounds],
  [Generation mode, v4 ÷ stage 1], [0.980–0.996 QCD/QED; 1.146 sm_eemm_L2 (more diagrams)], [reported],
)

The review's "about 40% fixed" is the share of the instruction gap between the
transplant and stage 1 closed on four workloads (35–52%). It is not a
completion measure, and it compares outputs whose diagram sets differ on SM
(69697 vs 81586 on eemm_L2). Earlier admission-lever estimates (option A on
`217c0695`) are superseded by \#193 and are not acceptance evidence.

== Settled decisions
<settled>

These are carried forward and not reopened.

- User decisions 1–11 and the binding list in HANDOFF section 6: ghost model A
  with the momentum-conservation reduction, names by canonical form with
  `name@id` pins and errors on unknown names, a separate `route()` stage with
  construction after enumeration, skeleton design 1, Classes::Auto default
  with a gate, mandatory loop counts, author Lucien Huber `<im@lcnbr.ch>`,
  publish only when happy, the Symbolica rev pin while crates.io lacks the
  needed commits, no compatibility shims, minimal filtering defaults, no
  vendored CeTZ, no GPL ec-layout or NewCM.
- Ground rules from HANDOFF section 1: work only in own workspaces, never
  rewrite others' commits, `nix-ci.nix` top-level `enable = false` until the
  CI-readiness sequence, shared cargo cache untouched, callgrind Ir and
  allocation counts as decision metrics (wall time only as ratios with spread,
  for workloads of at least 100 µs), gammabench for milestones and the
  ablation.
- U-1: the user personally readied and merged \#175 and \#176, redesigned
  \#177 and opened \#202. These are authorized.
- U-3: \#193 answers the admission cost and replaces levers A and Bx; C2 and G
  remain possible ablations on top of it.
- U-4: the separate lazy rule-set lane is dropped. Hosting \#193's template
  cache at model scope may be a small follow-up only if gammabench shows it is
  worth its code.
- U-2 is settled upstream by \#210: drop the first \#177 copy (`nxmztorq`),
  adapt the zero-colour path of the determinism lane (`yrpmmsun`) to \#210,
  and re-check 3150 and 4664/306.
- The integration must not duplicate \#160/\#193: take the tensor/MathML line
  from merged upstream \#209 at `3f0df2a8`; do not carry separate copies.

- U-6 (2026-10-27): preserve symmetrized gg→ggg total 57/8, including
  the 2!·3! symmetry division. WP4 uses this settled normalization.
- U-8 (2026-10-27): measure the current integration's full-mode memory first,
  then propose a justified limit. No ≤ 1.2× ceiling is pre-authorized.
- U-5 (2026-10-27): retain Clarabel and ship all required notices.
  Replacement-solver work stays separate from this integration.

== Open choices
<open>

These are the remaining decisions and review checkpoints. Technical findings
should first produce a concrete recommendation; only unresolved policy or
scope choices need a user answer.

+ #strong[U-7] U3 snapshot rows: which rows to re-certify under the standard in
  WP9 and which to mark failing with reasons. Blocks WP9 completion.
+ #strong[B2 remedy], only if WP1 finds a defect: the fix design follows from
  the diagnosis and needs review before it lands.
+ #strong[U-9, U-10] out-of-repo and housekeeping items (symbolica-fonts STIX,
  symbolica-community NOTICE, spynso3 plugin sources, cooked symbols behind 7
  feynkit-py scripts, U1 PR, extra \#191 test, gammabench `qxvqktlr` PR, d22
  issue). None blocks the integration.
+ #strong[Publication shape]: one squashed PR or a feature-reorganised series,
  decided at WP14.

== Work packages
<work-packages>

Each package lands as its own described change (or bounded evidence record)
with a durable bookmark, and finishes with Sol implement → Astra independent review →
fix. Exit criteria are necessary, not sufficient: reviewers may add checks.

=== WP0 Evidence preservation and source availability (no physics code)
<wp0>
- Preserve bounded integration evidence in the task-owned durable directory
  `docs/architecture/evidence/generation-integration/`; leave the historical
  handoff read-only. Index inaccessible raw evidence and its regeneration work.
- Pin the frozen v4, stage-1, transplant, upstream, v2/v3 and M1 source commits
  under new `generation-integration-source-*` bookmarks.
- #strong[2026-10-27 result:] approved by independent Astra round 2 (P2 closed);
  see #link("generation-integration-wp0.typ")[execution and availability record].
- #strong[Exit:] no decision relies exclusively on scratch: surviving relevant
  receipts are archived; inaccessible evidence remains reported and requires
  replacement checks in its owning package. Archiving is not certification.

=== WP1 B2 diagnosis on frozen v4 (independent, early)
<wp1>
Run before any rebuild churn, on `e0a77980` and upstream `58814885`, in its
own workspace. Bounded: one regression, one bisect, one written verdict.
- Treat sample-looking coefficients as a symptom, not proof. Decide whether
  any sample dependence survives exact cancellation.
- Exact comparisons set grouping `samples = 0`. The quotient receipt records a
  default of 5 sampled fallbacks after an inconclusive exact `compare_atoms`;
  leaving it on mixes the two paths.
- Primary check: take the original numerator pair of Group(28) and prove or
  refute `N_a = r · N_b` exactly, preserving numerator factorization, and
  compare `r` with the stored ratio. Evaluating the stored ratio alone is
  insufficient. Fresh-point evaluation of the pair is supplementary.
- Separately reproduce the production sampled/fallback path with its recorded
  budget and seed, then vary the sampling points. Check each accepted stored
  ratio against the original pair. A `samples = 0` test alone cannot exercise
  a defect confined to sampled fallback.
- Check MT dependence against the grouping's permitted scalar-ratio contract
  (read it from the code and docs): sample coordinates are never permitted;
  whether model parameters are permitted must be stated, not assumed.
- If the defect is reproduced, bisect `58814885..e0a77980` across the whole chain, including the
  transplant. The "nondeterministic Group(28)" entry was logged at v2 as
  identical at `75b6c7a8`, before the sampled-quotient changes entered v3, so
  `mpyquuut`/`msurxwzp` are suspects, not presumed causes.
- #strong[Outputs:] a minimal regression test (failing if the defect is real),
  the first bad commit when established, and a verdict: defect, permitted,
  test artefact, or unresolved. An unresolved correctness question blocks
  publication; it is not certified by matching samples.
- #strong[2026-10-27 result:] corrected certificate `f473709c` approved by Astra
  round 2. B2 is closed for the original GL028/GL034 pair under 4D SU(3),
  $N_28=(7/2)N_34$; stored sample and MT dependence cancels. No production fix
  or bisect is justified. Port the original fixtures, signed provenance,
  factor-preserving certificate and zero-sample/sampled controls in WP3.
- #strong[Exit:] verdict reviewed; ignored original-pair regression explicitly
  rerun on rebuilt (R) and final (F) tips.

=== WP2 Base and dependency
<wp2>
- Use the accepted fetched base `3f0df2a8` containing merged \#209 and \#210.
  Re-read its dependency contents; no unresolved \#209 choice remains.
- Fix the base: upstream tip containing \#210, including the merged tensor/MathML
  dependency, as one layer. Never carry separate
  copies of \#160 or \#193.
- List which v4 changes are now upstream (at least \#175 and \#176), which
  are superseded (the first \#177 by \#210), and which open PRs
  (\#166, \#172–\#174, \#178, \#191, \#195–\#197,
  \#204) are still needed, with their current heads.
- #strong[2026-10-27 result:] fixed target and all 46 source dispositions in
  #link("generation-integration-wp2.typ")[WP2 base manifest]; approved by Astra round 2 including `5226894e`.
  WP1 corrected certificate `f473709c` is approved by Astra round 2. The user-selected merged base
  is settled; no renewed base approval is required.
- #strong[Exit:] exact base/dependency manifest independently reviewed.

=== WP3 Rebuild v5
<wp3>
- #strong[Prerequisite slice status 2026-10-27:] Astra coordinator approved
  source03–07 at `9e2c4c23`, closing the sole P3 macro-import finding. The
  independent semantic review passed 418 tests with no semantic findings; the
  coordinator verified the exact correction and retained check hashes without
  inventing a delegated review task. Initial raw-admission protections remain
  approved and intact.
  The cumulative approved numeric08/09 patch (including production `19e0e9c2`
  and precision correction `c7272601`, evidence tip `25442024`) and scoped scalar14
  patch `872e84b2` are integrated on this base. Source13 now asserts the actual
  rendered `ExecutionMode.Single` attribute docstring. Sources10/11 are accounted
  for as upstream drops and12 as superseded by \#210, per approved WP2 evidence.
  Combined checks passed: 572 owner/control tests and two Spenso exporter tests,
  relevant all-target Clippy, formatting, kernel isolation, Typst and integrity
  checks. Astra approved the combined prerequisite slice at `438b6150`, with
  no new must-fix findings. Seven inherited registry digests remain explicitly
  pending source31/final documentation reconciliation.
  Fixed-key hashing and lazy fresh-index allocation remain deferred; only zero
  versus one discovery calls were established, with no end-to-end speedup claim.
  Approved cumulative drawing15/16 plus common notices are integrated at
  `8430e701`; the atomic resolved source17/18 transplant with isolated37
  formatting is reconciled against the admitted template owners. The
  original-pair B2 certificate passed on the staged tree with exact original
  fixtures and two explicit test-only presentation-tag bindings. This bounded candidate
  received aggregate REQUEST_CHANGES from all five independent Astra round1
  reviews at frozen `817b5886`. E4 at `edabf052` and all cumulative code
  corrections were subsequently approved and integrated, including D `b1061d7`
  and graph/Python `4b60161d`. Both independent combined integration reviews
  APPROVE frozen `d48f9d1a` / code `d109e898` within their assigned scopes;
  historical code and validation identities remain unchanged. The coordinator
  released source19–21 as the next bounded slice. It requires NEW independent
  Astra review before22; no whole-WP3 or performance/distribution approval is claimed.
  Source38 awaits the reviewed transplant and prerequisites through port36.
  Approved input25/36 and cumulative sampling27–29 remain separate for later
  integration; no benchmark milestone is claimed. Full rebuilt-tip validation
  and M-rebuild remain pending. See
  #link("generation-integration-wp3-staged.typ")[the staged integration record]
  and `evidence/generation-integration/wp3/prerequisites08-14.json`.

- `jj duplicate` the v4 chain (`kursouxx::kpslpwqz`) onto the WP2 base; do not
  rebase in place, so the 19 live workspaces on older lines stay valid.
- Drop the empty \#175/\#176 copies and the first \#177 (`nxmztorq`); adapt
  `yrpmmsun`'s zero-colour path to \#210 (it relied on
  `PreparedNumerator::color_products_vanish`).
- Squash the conflicted transplant `sptyyttx` with its resolution `vwptmtxs`
  and the rustfmt fixup `roxklnrz`; rewrite its description (what and why,
  4664/306, no local paths). Fix the stale descriptions and missing
  Co-Authored-By lines (m1).
- Reconcile \#191 with \#193; use one Symbolica rev in every lock; port \#193's
  generator templates into the staged API.
- Port the WP1 fix if one is agreed.
- #strong[Exit:] `jj resolve --list` empty on every commit; each commit
  formatted; CI clippy with `--deny warnings`; fast lane green; oracle fast
  tier at 1, 8, 32 threads; SM three-loop filtered 3150 at 1 and 8 threads;
  `four_loop_numerator` 4664/306 and the perturbation regressions at 1, 8, 32;
  `test_generate_sm_a_ddx` in 5 fresh processes; BRST quartet, epem NLO,
  CoolQM gate. Then milestone #strong[M-rebuild] (WP10). Any changed count
  needs the WP9 standard, not a re-record.

=== WP4 Names (decision 2)
<wp4>
Follow `decisions/names-routed-classes.md`, option (a1)+(d).
- Assign GL names in `Request::generate` after a parallel sort by linnet
  canonical form; accept `name@<DiagramId prefix>` pins in selections and
  named-LMB keys; unknown names, ids, keys and pins are errors that say where a
  pinned diagram went.
- Move `classes` out of the persisted `Request`, drop the serde fallback and
  its test (m3, P4-4), bump the state version once, fix the changelog versions
  (m4).
- P0 compares Auto and None names; fix the gg_hhh cards (P1-1), B1's scalar-3L
  tests and A5; re-record sm_jets with the settled 57/8 normalization including 2!·3!.
- #strong[Exit:] names identical across threads, Auto/None and grouping on the
  oracle; card selections resolve to their intended diagrams; pin errors
  tested; cost recorded (reported design cost +4.1% ungrouped, ≈0 grouped).

=== WP5 Routed stage (decision 3)
<wp5>
- `DiagramSet::route() → RoutedSet`, `RoutedSet::numerators() →
  BuiltDiagramSet`, `build() = route + numerators`; promote the test-only
  signed routing to the routed element; selection by id or name before
  numerators; Python and stub types.
- #strong[Exit:] `build()` output unchanged on the oracle; routed-only peak
  and Ir recorded; gammabench signed-routed adapter enabled.

=== WP6 Classes::Auto gate (decision 6)
<wp6>
- Control hash variance first: fix the hasher seeds for measurement, or
  measure the spread over several seeds and require the gate to exceed it
  (reported eager-region noise 0.4% against a 1.005 gate).
- Skip user dressing hooks on states with class edges unless they opt in,
  with a heredity test.
- #strong[Exit:] Auto/None ≤ 1.005 Ir in generation and full mode on the eight
  reference workloads plus qcd_tad_L3, qcd_gse_L4, qed_ase_L4, qcd_vac_L4;
  DHAT bytes ≤ +5%; heavy tier passes in both modes. If a metric fails,
  investigate and fix it or present the measured tradeoff for a user decision.
  Gate failure does not itself authorize changing the agreed Auto default or
  weakening the threshold.

=== WP7 Full-mode time
<wp7>
- Re-measure on M-rebuild against stage 1 and v4 with identical output
  contracts. Where sets differ, report the difference; per-diagram figures
  are supplementary and do not control for diagram complexity.
- Attribute the remaining gap by stage (admission, numerator, finalization,
  staging overhead).
- Remove the serial numerical-sampling loop in grouping and the per-call
  thread pool (`Run::in_pool`), or record them as dropped with numbers.
- Re-record `generation_ir`'s baseline at the tip and keep the Feynman-rules
  stage visible (`#[inline(never)]` on `construct_numerators`) (m9).
- Measure model-scope hosting of \#193's template cache; keep only under the
  ≳2% rule.
- #strong[Exit:] a stage table at the tip; each kept optimisation justified in
  WP10's ablation. No target ratio is pre-agreed; the user judges the result.

=== WP8 Memory and prepared-set retention
<wp8>
Expanded from the plan's retention item, which quoted 1.22–1.31× against a
different base.
- Diagnose with DHAT on gggg_L2, eemm_L2, ee4f_L2 in full mode: the retained
  prepared set (`prepared: Some(Arc::new(self.clone()))`), the `group(None)`
  clone, per-diagram SymbolicTensor fragment clones versus one structure per
  rule binding, and `HedgeGraphBuilder::build_with_map` growth (112 MB vs
  37 MB reported).
- Fix what has a clear owner; record the rest with numbers.
- #strong[Exit:] peak ratios versus stage 1 recorded at the tip, with the
  holders named; then propose a justified limit from those measurements; no pre-authorized
  ≤ 1.2× ceiling applies.

=== WP9 Slow lane, snapshots and heavy tier
<wp9>
- Classify every slow-lane failure as inherited (also fails upstream at the
  same base), new, or a justified expectation change. Only the last changes
  expectations.
- Re-certifying a snapshot needs group membership and weights, or full
  expression equivalence, plus an independent correctness argument. Equal
  group counts and summed weights are not enough; re-recorded hashes are not
  proof. Commit only expectations that pass (m2).
- Inputs: 25 U3 inline-snapshot mismatches and 21 uncertified rows,
  sm_jets (U-6), `finite_part_photon_3loop_no_ghost` d22 MD² residual
  (inherited; file upstream on the user's go), `photonic_amplitudes` near its
  cap (split it or give it budget), determinism residuals (sampled
  `UpToScalar` printing, open-network chain spelling) as named known-failing
  tests.
- Heavy oracle tier: record qcd_gse_L5 at 32 threads and the four 32-only rows
  at 1, 8, 32 outside nextest; run the tier once at the final tip.
- #strong[Exit:] slow lane green or every remaining failure listed with class,
  reason and owner; heavy tier recorded at the tip.

=== WP10 gammabench milestones and ablation
<wp10>
- Update the staged-API adapters to the rebuilt API (signed-routed after WP5).
- Milestones under both lazy and eager contracts: M-rebuild, one per decision
  stage (names, route, Classes), M-final at the commit to be pushed; recheck
  qcd_budget_2pt_l3 and qcd_2pt_l4 where M1 lost to FeynKit.
- Ablation with `gammabench variants` (balanced alternation, ≥ 21 rounds,
  counters on): one variant commit per optimisation off; keep if it gains
  ≳2% Ir or time on some workload class for its code, else drop; one row each
  with paired ratio and spread, Ir ratio, non-test lines, new caches or
  concurrency. Candidates: N1, F1–F4, P2 masks, P3 budgets, P4 classes, memo
  maps, routing cache, C2, G, template hosting.
- Run `.github/scripts/run_performance.py` locally as the CI gate.
- #strong[Exit:] registered milestones and the keep/drop table; dropped
  variants become removal changes.

=== WP11 Drawing and licence
<wp11>
Separate lane; it must not block generation delivery.
- Retain Clarabel and ship the required Clarabel and NetworkX notices in
  linnet-py, symbolica-host and gammaloop-api's embedded files (M7).
- Integrate only approved drawing layers: `kkumutpy`, then
  `particle(momentum)` labels. Replacement-solver work stays separate from
  this integration.
- #strong[Exit:] required notices present in every shipped artefact, or nothing ships.

=== WP12 Minor cleanup change
<wp12>
One small change for api-F7, api-F8, api-F9, A3's count basis against 7549
(m11), and the m10 leftovers that are in this repository (tydenso wasm rebuild,
two invalid DOT fixtures, batch-driver allocation, maturin sdist on
`vendor/clarabel`). Out-of-repo items stay with U-9.

=== WP13 Final validation and adversarial review
<wp13>
Full matrix below at the final tip, an independent adversarial review of the
whole stack, `just docs-check`, and a crates.io check for Symbolica/Numerica
(drop the `[patch]` only if a release has the needed commits).

=== WP14 Publication and cleanup (user's go only)
<wp14>
See #link(<gates>)[Gates for publication and cleanup].

=== Order and first slice
<order>

WP0 → (WP1 ∥ WP2) → WP3 → WP4 → WP5 → WP6 → (WP7 ∥ WP8) → WP9 → WP10
milestones run after WP3–WP6 and the ablation after WP8 → WP12 → WP13 →
WP14. WP11 runs alongside throughout.

#strong[First execution slice,] accepted 2026-10-27: WP0, WP1 and WP2.
They do not rebuild the integration stack. WP1 may require targeted historical
builds for diagnosis; the slice resolves B2 and the dependency base before
integration work starts.

== Validation matrix
<validation>

R = rebuild (WP3), N/Ro/C = names, route, Classes stages, F = final (WP13).

#table(
  columns: 3,
  [*Check*], [*When*], [*Pass condition*],
  [fmt, CI clippy `--deny warnings`, no jj conflicts], [every commit], [clean],
  [Fast lane], [R N Ro C F], [green; failures classified per WP9],
  [Oracle fast tier, 1/8/32 threads], [R N Ro C F], [all rows; headroom to 240 s noted (m8)],
  [SM three-loop filtered], [R F], [3150 at 1 and 8 threads],
  [`four_loop_numerator` + perturbation regressions], [R F], [4664/306 at 1, 8, 32],
  [`test_generate_sm_a_ddx`], [R F], [5 of 5 fresh processes],
  [BRST quartet, epem NLO, CoolQM 3σ], [R F], [as recorded in v4],
  [B2 exact comparison and production sampled/fallback regression], [R F], [WP1 correctness question resolved; accepted ratios valid for the original pair],
  [Name stability and pin errors], [N F], [identical across threads, Auto/None, grouping],
  [Routed vs built equivalence], [Ro F], [identical diagrams and weights],
  [Slow lane, heavy tier], [F (heavy also before M-final)], [WP9 exit],
  [Callgrind Ir, generation and full], [R N Ro C F], [recorded vs stage 1 and previous milestone; same output contract],
  [DHAT peak and blocks], [R C F], [recorded; U-8 threshold if set],
  [Classes Auto/None], [C F], [≤ 1.005 Ir, ≤ +5% DHAT bytes],
  [gammabench paired wall], [milestones, ablation], [≥ 21 balanced rounds; ratios with spread],
  [Docs check, wasm and SVG reproduction], [R F], [clean; tydenso rebuilt],
)

== Parallelism and owners
<owners>

- #strong[User:] unresolved policy and scope choices; authorization of
  execution; the explicit publication and historical-cleanup gates. Approved
  work packages include their routine implementation and validation steps.
- #strong[Coordinator:] plan upkeep, briefs, receipts, HANDOFF entries,
  integration of reviewed changes, durable refs.
- #strong[Implementer] (a delegated coding agent, per the user's current
  model preference): one work package at a time in its own workspace.
- #strong[Independent reviewers:] adversarial review per package; at least one
  not involved in the implementation.
- #strong[Drawing coordinator:] WP11 only.

Useful parallel lanes, each in its own workspace: WP1 beside WP2; heavy-tier
recording (hours) in the background after WP3; gammabench adapter updates
beside WP5; WP8's DHAT diagnosis beside WP7; WP11 throughout. Do not run two
lanes that edit the same files of the generator at the same time.

== Gates for publication and cleanup
<gates>

- No push, fetch, PR or label without the user's explicit go; then the
  CI-readiness sequence of `CONTRIBUTING.typ` (enable, validate, upload, push,
  `final-review`).
- Before any push: every commit conflict-free, described in plain words without
  local paths, authored by Lucien Huber; blockers B1/B2 closed; WP13 done;
  licence notices present if anything ships.
- Evidence lives under durable named bookmarks and the task-owned evidence directory,
  never only in `/tmp`.
- Cleanup is selective and needs the user's go: only our own stale bookmarks,
  scratch and empty-`@` workspaces after a per-item check; keep the M1 source
  (`generation-m1-source-9540f5bb`) and any baseline still cited; leave other
  contributors' and the drawing coordinator's refs alone; never touch
  `/common/build-cache`.
- Each package's report lists what it retained and why.

== Coverage of review findings and agreed items
<coverage>

#table(
  columns: 2,
  [*Item*], [*Where*],
  [B1 v4 not publishable (conflicts, base, first \#177)], [WP2, WP3],
  [B2 sample-dependent vacuum ratio], [WP1, WP3, matrix],
  [M1 full-mode time], [WP7, WP10],
  [M2 peak memory], [WP8, U-8],
  [M3 names, `route()`, Classes gate], [WP4, WP5, WP6],
  [M4 slow lane red], [WP9, U-6, U-7],
  [M5 superseded \#177], [settled by \#210; WP3],
  [M6 publication trail], [settled (U-1); record in HANDOFF],
  [M7 licence notices], [WP11, U-5],
  [M8 admission option A], [settled by \#193 (U-3); C2/G in WP10],
  [M9 lazy rule set], [dropped (U-4); template hosting in WP7],
  [M10 gammabench milestones and ablation], [WP10],
  [M11 heavy oracle tier], [WP9],
  [m1 history hygiene; m2 failing expectations], [WP3; WP9],
  [m3 classes serde fallback; m4 changelog], [WP4],
  [m5, m10, m11 minors], [WP12, U-9],
  [m6 retention, serial sampling, per-call pool], [WP8, WP7],
  [m7 determinism and U3 residuals], [WP9],
  [m8 oracle headroom; m9 `generation_ir`], [matrix; WP7],
  [m12 HANDOFF out of date], [coordinator, each package],
  [Agreed: kkumutpy, `particle(momentum)`], [WP11],
  [Agreed: A3 count basis], [WP12],
  [Info: `Plan::count` BigRational cost], [WP7 measurement; keep unless the user objects],
)
