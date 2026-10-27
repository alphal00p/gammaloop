= Generation integration base and dependency manifest
<generation-integration-wp2>

== Fixed target and review state

WP2 is ready for independent Astra review, not approved. Sol produced this
manifest on 2026-10-27 without replaying code. WP0 at `12582cbd` plus
`dc6f0b6b79f77c482e8747af0c7f12b817ffdaca` was independently approved by Astra
in round 2 (P2 closed, no remaining findings), as supplied by the coordinator.
WP1 is running separately in `generation-b2-diagnosis`; this record makes no
claim about its verdict. WP3 waits for both reviewed packages.

The agreed, already rebased integration target is
`3f0df2a870806076820786994d4883decec8e9e3`, tree
`3c91ba47525c627fbc02d79f5afa1d58aa7c5b66`.
The live GitHub `feynkit` ref matched that target when inspected; the timestamp
and response are retained. Do not fetch, chase a moving tip or rebase again.
Any later movement belongs to a coordinator decision.

The durable build input is
`evidence/generation-integration/wp2/base-manifest.json`. It records the exact
commit, tree, parents and JJ change identity of all 46 changes in
`588148851e3e989ecbb7478279db7d46727e23a8..e0a779808174c64057fecfd062c948b8de8e5c1f`,
with a disposition, reason, upstream match where applicable, planned unit and
ordered dependencies. Its source baselines retain all nine WP0 identities and
closures. Original source refs and the plan bookmark remain frozen. Historical
receipts, including inaccessible review counters, remain reported evidence.

== Dependency decisions

- Drop the v4 copies of merged \#175 and \#176. The \#176 patch is identical.
  \#175 has the same exact `flake.nix` patch, but the merged commit additionally
  changes nextest policy; its whole-commit patch ID differs. Preserve both
  upstream changes rather than describing the whole commits as identical.
- Drop the #emph[first] \#177 copy `nxmztorq`. GitHub reports redesigned \#177
  merged, but that does not make this first patch equivalent. The final owner
  is \#210 at `f02faa90d55c1d0a875b5345c010f03746d29549`, already in the target.
  Adapt `yrpmmsun` to its completed colour algebra beside indexed spectators;
  do not restore `PreparedNumerator::color_products_vanish` as a second proof.
- \#209 is merged at the fixed target. \#193 itself is closed #emph[unmerged];
  its authoritative tensor and admission-template work arrives via \#209,
  together with the tensor/MathML line. Neither \#160 nor \#193 is a separate
  layer to replay.
- The inspected PRs \#166, \#172–\#174, \#178, \#191, \#195–\#197 and \#204
  remain open. Their exact heads, bases, update dates and GitHub merge flags
  are in `wp2/pull-requests.jsonl` and the manifest. They are background
  dependencies, not PRs being modified by this task. Keep their required
  behavior through the mapped source units; do not stack their full heads
  over the same v4 copies. \#195 depends on \#178; \#197's published base is
  \#196. Generator registration follows planner/scoping; label placement
  follows plugin removal. The manifest preserves that order.
- \#191's scalar inference shortcut is not present on the target. Port its
  syntax eligibility and shortcut/full-admission equivalence tests against
  current metadata and proofs. Do not paste its obsolete proof construction
  over \#209's observed interface/history. This is reconciliation, not an
  upstream-drop claim.

== Assembly contract for WP3

Duplicate the frozen source range into a new task-owned integration line at
this exact target. Never rebase the frozen sources in place, alter the plan
bookmark, switch the WP1 workspace or rewrite another workspace's references.
Use the manifest order: prerequisites, one resolved staged-transplant unit,
then deterministic grouping, quotient/fallback handling, routing and API
fixes, docs/fixtures and colour-order completion. There are 41 planned units:
46 source changes less three dropped/superseded changes and two folds.
`retain` preserves intent and tests; it is not certification of automatic
patch application. `reconcile` explicitly requires adaptation to the target.

The transplant unit combines `sptyyttx 8cf6ceca`, its child resolution
`vwptmtxs 2b1ef4f0` and `roxklnrz cfe9e8d1`'s assertion formatting. Start from
the net resolved change across the first two, not Git's serialized JJ conflict
container. Fold only the final formatting delta from the later fixup; do not
pull intervening commits into this unit. Resolve against the target in the
#emph[new] copy and preserve source identities. Give each resulting commit a
plain description, correct Lucien Huber author identity and applicable
Co-Authored-By trailers. Check every resulting commit for conflicts, not just
the tip. No conflict-free rebuild is claimed by WP2.

Keep \#209's admitted interaction/propagator templates and metadata-preserving
`NumeratorInterface::translate` in the staged construction owner. Admission
happens outside cache write guards; publication retains one shared template.
Preserve actual-owner dynamic-rule admission, declared port order, hidden zero
history, injective label localization, reversed propagator interfaces and
concurrent publication tests. Carry these tests into the staged module rather
than losing them when replacing monolithic generation. Do not add a separate
lazy-rule-set lane or a parallel compatibility implementation.

The accepted target locks Symbolica, Numerica and Graphica to
`1a1829b1f3d0a88b53b44085b3b9f34ef77945b3` in root, Tydenso and notebook-host
locks; the manifest records each exact source string. Preserve that locked
alignment and host-subset rule with `--locked` builds; do not replay v4's old
kernel pin or wasm blob. Cargo and uv source declarations remain on `main`;
that branch spelling is not permission to update the resolved kernel. Any
required pin change must update all owners together and rerun feature
isolation, Clippy and GammaLoop tests. Rebuild generated assets only at their
owning implementation step, retaining Clarabel and required notices.

The WP1 reviewed verdict controls any quotient correction in `mpyquuut` and
`msurxwzp`. Matching sample values is insufficient: the original numerator
pair must satisfy every accepted stored scalar ratio. WP2 introduces no
physics correction and does not decide B2.

== WP3 validation and corrections needed

`wp2/validation-matrix.json` gives commands, owners, thread settings and pass
conditions. Every command runs from the integration workspace through
`nix develop -c` (or explicit `direnv exec .. nix develop -c` if a shared envrc
is later configured and allowed). Preserve an exported `SYMBOLICA_LICENSE`,
otherwise load it silently from the current user's XDG config. Never print it.
Use locked dependencies, `INSTA_UPDATE=no`, no forced-pass variables and no
oracle record mode. Preserve test assertions and current fast-lane limits.

Run fmt, CI Clippy with denied warnings and a conflict check for every resulting
commit. At the rebuilt tip run the fast lane, oracle at 1/8/32 threads, SM
filtered three-loop 3150 at 1/8, four-loop 4664/306 and registration/colour
perturbation tests at 1/8/32, five fresh SM a→ddx processes, BRST quartet,
epem NLO, fixed-seed CoolQM at its accepted three-sigma gate, and both WP1
exact and production-fallback regressions. Port and run current target
admission/template regressions as well as source staged tests. A passing
command that selects zero tests is a failed validation attempt.

Two harness discrepancies require bounded implementation work in WP3:
+ `four_loop_numerator` currently asserts only that retained plus zero equals
  4970. Its logged split must additionally be checked as 4664 retained and
  306 zero with colour filtering enabled; do not accept 4970 alone.
+ The source SM filtered test explicitly sets `n_cores.feyngen=1` inside its
  body. An external thread variable cannot certify eight threads. Parameterize
  its execution setting or add an equivalent eight-thread invocation retaining
  the entire physical selection, `3150 | -9861/16` assertion and the two
  Lorentz-zero graph checks. Do not edit the expected values.

The source oracle recipe mentioned in its module comment is absent from the
accepted target recipes. Use the explicit cargo-nextest command in the matrix
rather than claiming a nonexistent recipe passed. Ported oracle baselines and
historically refreshed ghost/sunset and aa3 rows still require WP9's exact
justification. Do not regenerate snapshots to make the rebuild green. Keep
symmetrized gg→ggg 57/8 including 2!·3!, and report differing output contracts
before comparing measurements. M-rebuild records Callgrind generation/full
and DHAT counters without a preauthorized 1.2× memory ceiling; the coordinator
runs the milestone through WP10. Slow/heavy certification remains WP9.

== Verification and retention

The WP2 focused receipt records JSON/TOML structure, exact 46-change coverage,
commit/tree/parent identities, ordered dependencies, baseline closures, source
and plan ref integrity, registry/title contracts and Typst compilation.
No code, snapshot or kernel lock changes, physics experiments, full builds,
source rewrites, PR mutations, publication or historical cleanup were performed.
Top-level NixCI enable remains false. Retain this workspace, the named
implementation ref, manifest, PR/ref observations, matrix and verification
receipt for Astra and WP3. One-off scripts and rendered PDFs are removed;
shared caches remain untouched and no background process is retained.
