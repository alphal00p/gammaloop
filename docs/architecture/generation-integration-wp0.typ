= Combined generation integration WP0 execution record
<generation-integration-wp0>

== Accepted scope and ownership

On 2026-10-27 the user authorized execution of the accepted combined generation
integration plan: Sol implements each bounded step and parent Astra assigns a
separate Astra review after each step. This record covers WP0 only. Astra
independently approved round 2 after
`dc6f0b6b79f77c482e8747af0c7f12b817ffdaca` closed P2; no findings remain.
This record does not begin WP1, WP2 or the v5 rebuild.
Routine steps within accepted packages need no renewed approval; publication
and destructive historical cleanup remain gated.

The task workspace is
`/shared/lcnbr/dev/gammaloop/generation-integration-plan`. The preserved plan
change is `kskntrlv 4374f720`; implementation is its child under bookmark
`generation-integration`. The accepted base is `3f0df2a8`, with \#209 merged;
\#193 admission arrives through that merge. There is no unresolved \#209
dependency choice. U-1 authorization is settled, \#210 supersedes \#177, and
the separate lazy-rule-set lane remains dropped. The following are settled as of 2026-10-27: preserve symmetrized gg→ggg
total 57/8 including the 2!·3! symmetry division; measure current integration
full-mode memory first and then propose a justified limit (no pre-authorized
≤ 1.2× ceiling); retain Clarabel and ship required notices. Replacement-solver
work stays separate from this integration. Do not reopen these choices.
Snapshot certification remains a later WP9 decision.

== Durable evidence

The task-owned archive is `docs/architecture/evidence/generation-integration/`.
`evidence-inventory.json` records exact original paths, byte lengths, SHA-256
hashes and archive-relative paths. Raw JSON receipts are preserved byte for
byte. They are historical reports, including superseded findings and embedded
scratch paths; they are not new test results or authorization. The canonical
plan and this record determine current scope. The original handoff, review and
decision prose remain read-only at their durable external locations, with
hashes in the inventory. No historical handoff was overwritten.

The archive covers the ghost decision, stage-1 and rebase provenance, planner
and determinism preparation, v2/v3/linearization and quotient receipts, bounded
U3 provenance, M0/M1 and counter-adapter provenance, and preparation for names,
routing and Classes. Large binaries, caches, raw Callgrind/DHAT output and
unrelated drawing-coordinator work are excluded. No cleanup of external
scratch or historical refs was performed.

The 2026-10-22 review's original `review-v4` scratch and the original naming
prototype are inaccessible with permission denied. The inventory gives their
exact paths and regeneration packages. Consequently the review's 364/364,
3150, 4664/306, instruction ratios, memory peaks and B2 symptom remain
#emph[reported], not recertified by WP0. WP1 must recover the original pair
from the frozen source and perform the exact and production-fallback checks;
WP7/WP8 must record new counters and WP9 must run the slow/heavy matrix.
No decision treats inaccessible evidence as successful validation. Archiving
receipts removes the exclusive dependence on scratch without closing these
correctness or performance gates.

== Historical source availability

`sources.json` records full commit, tree, parent and local JJ identities and
new task-specific bookmarks for v4 `e0a77980`, stage 1 `f25961e2`, transplant
`4e008517`, upstream `58814885`, v2 `d9ff71a6`, pre-v2 `75b6c7a8`, v3
`dfe63876`, original v3 `6d207cfc` and M1 source `9540f5bb`. The two v3
commits have the same tree `1fc62584e241d44d0e2352aa848201f99830a0ee`;
their histories must not be conflated. Historical v4 change id `kpslpwqz`
is a source-repo identity; local imported identities are recorded separately.

Both source repositories reject JJ's secure-config read at
`.jj/repo/config-id` with permission denied, even with
`--ignore-working-copy`. The Codex checkout's guidance also rejects direct
file reads; its Git HEAD guidance was read through command-scoped
`safe.directory` object access. That same genuinely read-only access resolves
all required GammaLoop commits in `/shared/codex/dev/gammaloop/main`.
The requested commits were not found in the Git store at
`/common/dev/gammaloop/feynkit-rich-display`, whose guidance was read first.
No source metadata, permissions, working copies or refs were changed.

`object-import.json` records the local pack transfers. Objects were copied
into our repo, without remote fetch or push. New
`generation-integration-source-*` refs retain the baselines; existing refs
remain untouched. Reachable-object checks verify source closure availability,
not buildability. JJ import recovers the original change identities, including `kpslpwqz`,
and recognizes the historical conflict in `sptyyttx 8cf6ceca`. WP3 must
handle that conflict as planned rather than assume the historical stack is clean.

== Blockers and independent work

No broad missing GammaLoop source blocker remains for WP1 or WP3 once the
local closure checks pass. Builds on historical pins can still require their
recorded dependencies and license; WP0 does not test those environments.
Raw review evidence is inaccessible, so future gates require fresh evidence.
The gammabench M1 source and snapshot availability are recorded separately;
M1 remains local historical evidence and is never relabelled as a fresh run.

WP1 can proceed from the pinned v4 and upstream sources after parent review.
WP2 can inventory the accepted merged base independently. WP3 remains a later
package after WP1/WP2, including its B2 verdict and dependency reconciliation.
The settled normalization, memory-measurement policy and Clarabel decision
do not prevent those independent checks.

== Verification and retention

Focused Nix-backed verification checks archive hashes and JSON syntax, source
object closure, documentation registry/title contracts and Typst syntax.
`verification.json` records commands and outcomes, including any limitation.
No physics tests, B2 experiments, v5 rebuild, snapshot edits or heavy suite run
belongs to this step. Top-level NixCI enable remains false.

Retained for review and continuation: this workspace, the implementation and
source bookmarks, imported source objects, concise evidence archive and check
receipts. Shared caches and external scratch remain untouched. Task-owned
one-off inventory scripts and rendered check outputs are removed after use;
no background process is intentionally retained.
