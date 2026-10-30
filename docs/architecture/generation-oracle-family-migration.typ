= Four certified numerator oracle migrations

This is a #emph[separate migration preparation candidate], not MAIN integration.
The mathematical certificate at `87fa13a93b1c0b3e5792fd56937486ad962b2b55`
is approved by root aggregation of independent round-three A and B reviews.
The new migration and its runtime evidence require a NEW independent migration/
test review. No review is self-started here.

== Exact change and approval boundary

The approved certificate is preserved at `generation-oracle-family-approved`
and `generation-oracle-family-round3`. All earlier refs, reports, receipts,
checkers, captures and reviewer evidence remain unchanged, including the large
round-three reviewer pair/control trees. Historical source remains pristine
over `2b1ef4f07a44bc061bcc97c378bbdd9c41cb0e94`, with diagnostic `5e61ef89`
preserved. Both approval JSON files are copied byte-exact into this evidence
directory and bound by full SHA-256 in the receipt.

Exactly four numerator literals change in
`crates/feynkit-generator/tests/generation_oracle/main.rs`; no other code line
or Row field changes. All three `xs_qed_*` literals remain as they were in THIS
preparation base, deliberately not imported from MAIN. The independent FNV
recomputation uses the actual certified Classes None payloads and the immutable
source's exact JSON/field ordering, little-endian integers and multiset encoding.
It reproduces all historical Rows and every non-numerator field of the current
Rows: counts, rational weights, inverse symmetry, histograms, content, named and
ID digests. Hash equality is inventory/digest evidence, not mathematical proof.
Physical weight handling admits only these observed exact AutG and fermion-sign
forms; it never substitutes arbitrary physical parameters at one.

#let receipt = json("evidence/generation-oracle-family-migration/receipt.json")
#table(
  columns: (2fr, 2fr, 2fr),
  table.header([Independent case], [Old literal], [Current literal]),
  ..receipt.rows.map(r => (raw(r.case), raw(r.old_literal), raw(r.new_literal))).flatten(),
)

== Producer qualification and actual runtime matrix

Here “old” means the #emph[current-route producer with old expectations], not
the historical c78 kernel. The old retained ELF is SHA-256 `133899560973c3b824a1fde41061c0cf32995438434436567f1bfaa035d9d7ba`.
Its source/compiler/kernel provenance is corroborated by the approved reviews,
and its actual 341-unit fingerprint closure was retained before the cache
artifact was replaced by normal candidate compilation. Original producer final
build-log gaps remain gaps; no original transcript is invented.

The source-qualified new ELF is SHA-256 `49ed27bd4a26c706d1785a7444e2de4001669a22e77d5d545b8987fe594a9df2`.
The qualified build has actual source hashes before and after, compiler/Cargo
outputs, emitted artifact metadata, full fingerprints and retained ELF. A
byte-preserving mtime touch forced a genuine compiler invocation; no cache was
deleted or transplanted and no CARGO_TARGET_DIR was overridden. Old/new root
compiler, test profile, features, flags/config and all 341 dependency units
compare exactly. “Root default/native” is not used as a substitute for dependency
feature qualification. Both assertion ELFs use locked kernel
`1a1829b1f3d0a88b53b44085b3b9f34ef77945b3` and rustc 1.98.1 (48a229cea).

#table(
  columns: 5,
  table.header([Actual route], [Selections], [Passed], [Failed], [Ignored]),
  ..("old", "new").map(name => {
    let m = receipt.matrix.at(name)
    (name, str(m.selected), str(m.passed), str(m.failed), str(m.ignored))
  }).flatten(),
)

Each route executes exactly these four real test names independently at actual
pools 1, 8 and 32. Every ordinary run selects one nonignored top-level test;
no nested-child, zero-selection or record-mode result receives pass credit.
The old failures are confined to the numerator literal: complete observed
left/right assertion Summaries are compared with independently recomputed Rows.
All twelve new runs pass the original full-Row assertion and class/Named/Minimal
checks. The final matrix is `runtime-final/`; the earlier repeated matrix is
supporting development evidence, not silently substituted for final runs.

Each runtime uses a private tmp/state/cache and no-network read-only sandbox,
with only its output root writable. Record/diagnostic/update/force-pass flags
are cleared; only the actual thread selector and explicitly authorized dump
flag are set. Every Nix evaluation/develop/build has BOTH license variables
unset. Only the direct runtime exec receives the real private configuration
license; no value, secret-derived hash or secret-bearing derivation is recorded.
Invalid or absent licenses can yield correct restricted-mode output: license
validity is not an equality predicate or a correctness theorem. Policy failures
and infrastructure failures never count as inequalities.

== Captured owner, interface and mode correspondence

The bounded observer reuses the approved actual model/request/options, Run,
Census and typed-tensor owners in a private approved-source copy. Its sole
observational overlay adds a pool argument and actual pool folder names; it
changes no generation option or numerator algorithm. The emitted ELF and its
actual source/fingerprint/build log are retained. It is the same observer for
old/new routes because assertion literals are not compiled into it. Its build
root profile is distinct from the assertion test profile; complete records,
not generic profile equivalence, qualify that boundary.

Twenty-four fresh observer executions supply 48 mode bindings. Complete actual
None/Auto assertion dumps and source companions are joined to their observer
payloads, not merely by FNV. Each model, parameter/request/options/grouping,
graph key/ID/multiplicity, physical field, source factor, owner structure,
expression, typed flag, interface #emph[presence and order], and ordered port
list is compared type-strictly with the corresponding approved current capture.
Run debug text in capture.json is explicitly inspected for the actual
Some(1/8/32), classes, Grouping None, prepared and built counts. It is not falsely
claimed to be consumed by the mathematical validator. All actual same-mode
payloads bind the approved data at every pool; their exact observed byte relation
is retained in validation.json, alongside full-field comparisons.

CT None/Auto are not globally byte-identical: complete canonical-key/ID joins
retain the certified output-order permutation, source indices, groups and
factors. None positions 109–116 map to 113,109,114,110,115,111,116,112; other
positions are fixed. Names and group/member identities follow that mapping on
both routes and at each pool. No new normalization is invented. Mathematical
approval remains four separate cases × None/Auto, Grouping None, #emph[pool 1].
Pools 8/32 are additional current-route runtime observations, not historical
mathematical certification.

Observer typed metadata is captured within the observer producer, then linked
to each assertion producer through the complete public payload. It is not
in-memory instrumentation of the assertion ELF. Port typing against expression
free indices, factor locality and downstream positional consumer correctness
remain source-qualified/capture-bound; the theorem is named multiset transport.

== Independent replay and negative controls

The approved portable entrypoints are executed unchanged, with fresh outputs
and no dependence on migrated literals: 8 baselines, 531 positive/negative
controls and two separately classified infrastructure checks reproduce their
intended exact stages. Both original exact-checker ELFs and the original code/
raw inputs remain bound. Baseline anchors are self-comparisons, not extra
algebra. B's additional independent held-outs are retained and cited, not
presented as permanently exercised by the committed 531 suite.

Reapplying the original sealed root-manifest guard yields an #emph[expected
failure] at the changed registry only, with zero validation credit. Its original
recorded registry hash matches the immutable approved donor. No guard is
disabled, no frozen output is overwritten, and the whole outer closure is not
called green. Portable semantic replay and separate immutable preservation
inventories provide the appropriate scoped evidence.

The actual new ELF with a nonexistent exact selector exits zero but selects
zero tests; that is rejected as pass evidence. The old ELF in deliberate record
mode exits zero with one test by skipping its Row assertion; RECORD presence
disqualifies it. Ordinary old-ELF failures are meaningful wrong-literal controls.
All such negatives retain argv, flags, selections, exits and raw outputs.

== Review qualifications and preservation chronology

#table(
  columns: (1fr, 5fr),
  table.header([Qualification], [Disposition / retained limit]),
  ..receipt.qualifications.pairs().map(pair => (raw(pair.at(0)), pair.at(1))).flatten(),
)

The strictly pre-run round-three development-before-v1 inventory and its
per-run guards are distinguished from frozen-before, whose hashing overlapped
initial runs. Both development inventories remain superseded, not rewritten.
The old mode-order artifact and its real pre-final-common.py producing command
are attributed explicitly; migration rederives the correspondence independently.
New global preservation snapshots are not retroactively described as pre-run:
actual per-run source/input/ELF guards and immutable approved indexes provide
that strict consumed-input boundary. Unrecorded old edits/renames or missing
logs are never reconstructed as recorded commands.

Failed and preliminary preparation remains retained: three stalled Nix clients
were interrupted before compilation/runtime; the first qualified inventory
encountered a tracked gitlink; the first observer fingerprint lookup assumed
the emitted path held fingerprints, whereas the normal cache forwards the ELF.
Retries record the actual normal-cache binding. None of these failures supplies
pass credit. Helper versions and pre-format copies retain source chronology.

== Static source ownership, docs and handoff

Root's read-only MAIN30f preflight is SHA-256
`3dc3a7f7b3bff69249feec35e5c7f54439b0b03ae7db711777dbb8f6820a22e9`.
It is context, not import authority: at that immutable head the four complete
old Rows match the approved base, are disjoint from the three migrated QED rows,
and fall outside existing 113 scope ranges. Future reconciliation must refresh
that basis and apply #emph[only the four-line patch], never transplant main.rs.
No physical MAIN file/index is changed here, and no full MAIN suite is run.

Archive paths/bytes remain unchanged. Being outside Cargo auto-target discovery
and scoped generator/Idenso/core compile roots does not exclude archived Rust
from broad commonCargoSources or broad build/fmt/test/WASM/docs source sets.
Docs runtime includes the entire docs tree. These are static ownership findings,
not evaluated derivation/NAR identity theorems; actual future combined sources
require qualification. No rename, filter repair or policy change is made.

The new registry owner is inserted without changing any of the three certificate
classifications/scopes. Canonical titles, pinned PDF/semantic HTML, edition-2024
formatting, scoped Clippy and the exact code diff are checked. Full docs is
#emph[not green]: the inherited eight HTML paths (plus preserved reviewer-owned
commands.md in the live tree) are the first gate, followed by masked unknown
sunset owners and source46 classification. No unrelated docs repair is imported.

All other 74 of the remaining 78 MAIN rows stay unchanged and uncertified.
No kernel/lock/manifest/NixCI, production algorithm/options, physics, grouping
or other expectation changes occur. No whole-graph expansion/distribution,
trace, numerical or parameter-at-one substitute is used. Work stops frozen for
root's NEW independent migration/test review, then only later refreshed MAIN
reconciliation and a separate combined review. No push/fetch/upload/PR/CI or
publication is performed.
