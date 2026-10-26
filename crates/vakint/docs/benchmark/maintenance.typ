#set document(title: "Benchmark report maintenance")
#set page(paper: "a4", margin: 2cm)
#set text(font: "Libertinus Serif", size: 10pt)
#show raw: set text(font: "DejaVu Sans Mono", size: 8pt)

= Benchmark report maintenance

The four-loop collector and refresh command remain in
#link("../../AGENTS.md")[Vakint's contribution instructions]. Preserve every
historical root field, raw sample, source identity and selected four-loop run.
The five-loop adapter #raw("merge_five_loop_results.py") extends those records;
it neither runs reductions nor authenticates owner proofs. Its file-contract
tests use synthetic temporary files, never measured or reference physics data.

== Read-only native cache suite measurements

#raw("measure_cache --cache PATH") accepts the existing compressed single file,
a closed native bundle directory, or a closed standalone native pack directory.
Held-out suites forbid #raw("--running-cache"): ordinary #raw("load_cache") may fill
in-process memory but never captures records or modifies the supplied cache.
Preserve all prefixes separately and benchmark an immutable completed snapshot.

Single-file records retain their existing #raw("cache_sha256") field unchanged.
Directories instead set #raw("cache_native_index_identity"), with SHA-256 and byte
lengths for the top bundle index and every child native index, or the standalone
pack index. Its aggregate SHA-256 hashes the compact JSON schema/layout/index
record before adding the aggregate itself; relative names are sorted and cache
location is excluded. The existing before/after identity check rejects changed
indexes. This commits to native bindings, generations, ordered runs and referenced
chunk hashes without hashing a multi-gigabyte payload for every worker. Native
loading and lookup still validate bindings, framing and used chunks; this is not
a claim that unused chunk payloads have been independently checked.

Identity checks remain outside the reduction timer. The existing
#raw("cache_loading_seconds") measures cache attachment separately, and native
phase timings retain their documented combined scope: pure IBP and terminal
substitution are not independently measured. Full public result comparison and
the canonical sorted, deduplicated guard export remain unchanged; that guard
export is not the native ordered guard sequence. Final performance acceptance
still requires the separately recorded serialized window and matched repetitions.

Public five-loop data loading validates and activates the staged saved application
even when prepared-cache loading is disabled. Report that setup separately from
cache import and attachment: import validates and installs a snapshot without
attaching it. Individual formulas may still be prepared lazily during reduction.

== Prepare the physical input once

The complete quadratic numerator preparation is a separate measured stage. Build
#raw("five_loop_qcd_ladder") with the normal release profile, then save its complete
scalar input once:

```sh
five_loop_qcd_ladder --work-dir /path/to/workspace/quadratic-input \
  --fresh-prepare --prepare-only
five_loop_qcd_ladder --work-dir /path/to/workspace/uncached-run \
  --prepared-input /path/to/workspace/quadratic-input --artifacts /path/to/external-data
five_loop_qcd_ladder --work-dir /path/to/workspace/cached-run \
  --prepared-input /path/to/workspace/quadratic-input --artifacts /path/to/external-data \
  --cache /path/to/selected-cache.rrcache.gz
```

Each reduction remains a fresh process. All 21 quadratic contributions, mass
factors and projector conditions come from the same verified input; this is input
reuse, separate from the optional persistent reduction cache. The manifest
#raw("prepared-input.json") authenticates exact scalar/condition bytes, projection,
scientific definitions and producer source identities. The collector automatically
reuses a valid bundle in its work directory when no explicit input is selected.
An invalid bundle, including payloads with a missing manifest, errors before any
native regeneration. #raw("--fresh-prepare") requires a new work
directory and explicitly validates native generation again; it never overwrites a
completed bundle. #raw("--prepare-only") can also validate an existing bundle.

Prepared-input reports set #raw("input_mode: prepared") and
#raw("fresh_native_preparation: false"); current #raw("preparation") and
#raw("preparation_seconds") are null. The measured #raw("input_loading_seconds")
covers verification and parsing. The original producer's preparation timings and
source hashes remain under #raw("prepared_input.manifest") as historical evidence,
excluded from the current end-to-end duration. The manifest path and SHA-256 bind
that provenance. Current consumer fixture hashes are recorded separately. Never
combine fresh and prepared-input samples into one comparison, or describe a reused
input's end-to-end time as fresh graph-to-result time. The reduction timer itself
keeps its existing boundary after external-data/cache setup and before result
classification/serialization.

The ignored slow acceptance test follows the same reuse policy in
#raw("VAKINT_QCD_LADDER_EVIDENCE"). Optionally select an external bundle with
#raw("VAKINT_QCD_LADDER_PREPARED_INPUT"). Force independent regeneration using
#raw("VAKINT_QCD_LADDER_FRESH_PREPARE=1") with a new evidence directory and no explicit
prepared-input environment variable. Census, complete exact reduction, warm parity,
and full nonzero conditions remain checked in both modes.

== Append evidence without changing history

```sh
python3 crates/vakint/docs/benchmark/merge_five_loop_results.py \
  --base crates/vakint/docs/benchmark/data/results.json \
  --evidence /path/to/workspace/evidence-index.json \
  --output /path/to/workspace/new-results.json --select-current
```

The output must not already exist. Omit #raw("--select-current") to append while
retaining the old selection. Completed strict cache preparation can be included
with #raw("--build build.json --cache cache.gz"); completed suites additionally
use #raw("--held-out suite.json") and/or #raw("--tuning suite.json"). Suites require
their exact cache build and proof/cohort bindings. Duplicate evidence runs fail.
No root value or historical run is rewritten; only #raw("five_loop_runs") and an
explicit #raw("selected_five_loop_run") are added. Review the new JSON and rendered
PDF before replacing the canonical trio together. Do not update canonical data
while a measurement or source-freeze window is active.

The same #raw("--build RECEIPT --cache PATH") accepts an authenticated
#raw("vakint.cumulative-native-cache-report.v1") receipt with an immutable native
cache directory. That receipt binds the collector's native index identity,
semantic program binding, data-manifest hash, compressed-byte/chunk inventory,
raw stored-reference count, original source receipts and hashed native validation.
Stored references include shadowed entries and are never presented as unique keys.
A retained incomplete source remains incomplete, even after successful opaque
preservation or union into a new snapshot.

Cumulative report rows distinguish #raw("preserved-committed-prefix") from
#raw("validated-committed-snapshot"). They list individually validated completed
physical examples; #raw("physical_examples_complete") requires both the
zero-momentum diagnostic and the complete quadratic fixture. The separate
#raw("all_source_runs_completed") also accounts for retained incomplete prefixes.
Snapshot preservation alone has #raw("performance_eligible: false") and never
supplies missing timings or pretends to be a completed finite cache workload.

Final cumulative-cache suites additionally require an independently accepted
#raw("vakint.cumulative-cache-consumer-requalification.v1") receipt. It binds the
original source starting-key inventories and route/alias exclusion against the
unchanged frozen owner sets, explicitly classifying exposed starting inputs.
The existing exact result/guard, binding/cohort, twelve held-out sets, three fresh
matched samples and serialized-window requirements remain in force. Useful
incidental descendant and root hits stay visible. There is no aggregate
50/95 MiB rejection, and none of the original four-loop history is rewritten.

== File bindings and functional records

An evidence index has schema #raw("vakint.five-loop-report-evidence.v1") and a
#raw("records") array. Each entry contains #raw("case_id"), #raw("record"),
#raw("launcher"), #raw("artifact"), and an optional #raw("window"). Every reference
is a #raw("{\"path\": \"...\", \"sha256\": \"...\"}") object. Relative paths resolve
against the index directory. All JSON references are rehashed on every invocation;
missing evidence is an error. Retain raw source receipts alongside the published
report; paths alone are not portable archives.

Primary case IDs are #raw("quadratic-four-rung-uv") and
#raw("zero-momentum-diagnostic"). Input reports come from
#raw("five_loop_qcd_ladder") with the corresponding projection. Hashes bind the
actual executable, exact input, result, complete nonzero-condition file (including
preparation conditions), and full unknown-master inventory. Retained success and
catalogue closure are distinct. Unknowns are preserved explicitly, and prevent a
primary performance claim. Closed retained masters still need not have numerical
coefficient data.

Historical #raw("vakint.q0-full-weighted-native-diagnostic.v2") records are also
accepted as exploratory functional evidence, with their original compiled cohort,
source/data identities and launcher bindings. Their b72/opt2 provenance is never
relabeled as a final release or converted into a median. Ordinary concurrent
functional launchers also remain functional-only, regardless of elapsed time,
release optimization, affinity or a zero exit code.

== Final performance eligibility

The adapter consumes an *existing, independently accepted completion receipt*;
it does not issue eligibility decisions or launch/schedule workers. An index's
window reference must use schema #raw("vakint.serialized-performance-window.v1")
and explicitly record all of:

- #raw("accepted"), #raw("final_performance_eligible"), #raw("serialized"),
  #raw("all_processes_drained") and #raw("source_unchanged") as true;
- actual numeric #raw("started_unix_seconds") and #raw("finished_unix_seconds"),
  nonempty #raw("cpu_affinity"), and #raw("thread_environment") fixing
  #raw("RAYON_NUM_THREADS"), #raw("OMP_NUM_THREADS") and #raw("MKL_NUM_THREADS") to
  the string #raw("1");
- a #raw("samples") array with each actual report, launcher and artifact SHA-256,
  case ID, variant, integer sample index, start/end timestamps, fresh-process
  observation and successful exit code. Intervals must lie inside the window
  and never overlap. The report/launcher/artifact references must match exactly.

Each final launcher must itself explicitly declare
#raw("final_performance_eligible: true"), #raw("fresh_process: true"), the actual
#raw("binary_sha256"), case/variant/sample identity, matching affinity and thread
environment, and the manual #raw("cache_sha256") for cached primary samples.
Do not add these declarations to a concurrent functional launch after the fact.
The accepted receipt must describe the observed serialized execution.

For suite samples, each launch/window entry must exactly match the measurement's
actual #raw("supervision") fields: #raw("child_pid"), #raw("supervisor_pid"),
#raw("started_unix_seconds"), #raw("finished_unix_seconds") and
#raw("interval_scope"). That scope is “supervisor pre-spawn setup through
supervisor-observed worker completion”. These are observation boundaries;
the existing wait polling does not identify an exact operating-system exit time.
Missing observations are not reconstructed from deadlines, elapsed values or
filesystem timestamps. Primary direct launchers retain their own observed process
facts. The serialized-window acceptance remains external to both collectors.

Compiler artifacts may use the existing single-primary or three-driver artifact
shape. They must identify the actual normal optimization-3, fat-LTO executable,
release profile, unchanged source, Main/Native/kernel revisions and both frozen
source snapshots. Available immutable ELF files are rehashed; an absent archived
ELF does not erase its mandatory identity. Runtime Git HEAD alone is insufficient.

A primary comparison needs at least three matched fresh samples *per mode*, one
artifact/cohort, representation, input, fixture/data identities, result/guard set,
cache file and accepted window. All automatic 1L–5L backings must be absent;
only the cached mode manually attaches the selected 5L file. An otherwise valid
functional record without a window is retained, with no aggregate timings.

Suite qualification binds consumer revisions/source snapshots to its own
owner proof and requires one adopted window. The cache producer retains its
original provenance. Final held-out
results need at least twelve distinct nonempty protected owner sets. Tuning and
held-out rows stay separate. Missing closure fields in an older collector are an
error; a successful old retained result is not silently reclassified. Updating a
collector/cohort or admitting separately validated classification is a distinct
review decision, never an adapter fallback. Primary performance also requires
its completed strict cache build, with the identical cache and data-manifest
hashes; functional records need no cache-build receipt.

== Requalifying an existing cache for a new consumer

Identical producer/consumer workload, proof and revision follow the existing
strict path. A changed consumer additionally requires the index field
#raw("cache_requalification: {path, sha256}"), referencing an independently
accepted #raw("vakint.cache-consumer-requalification.v1") receipt. Its
#raw("accepted") flag and these #raw("checks") must all be true:
#raw("canonical_consumer_proof"), #raw("unchanged_selection"),
#raw("twelve_distinct_held_out_sets"),
#raw("training_raw_and_nonzero_owner_disjoint"), and
#raw("cache_import_validated"). No acceptance is issued by this adapter.

The receipt binds hashed producer #raw("cache"), #raw("build_report"),
#raw("workload"), #raw("owner_proof") and #raw("compiled_cohort") references;
consumer references are #raw("collector_artifacts"), #raw("compiled_cohort"),
#raw("qualification_workload"), #raw("final_workload"), #raw("owner_proof"),
#raw("protected_union"), #raw("training_route_qualification") and
#raw("alias_exclusion"). Paths resolve relative to that receipt. Cache bytes are
hashed as an opaque binary file; the other receipts are JSON.

#raw("immutable_selection") records #raw("source_census_sha256") plus ordered
#raw("training"), #raw("tuning") and #raw("held_out") arrays of
#raw("{id, powers, public_input_source_sha256}"). Input hashes use exact UTF-8
source strings; only an originally absent input has a null hash. The adapter
compares this selection against both original and new workloads. Only qualified
lowered/protected vectors and proof metadata can change. Cache SHA, census,
D40 data-manifest SHA and exact #raw("semantic_binding") stay equal. Both original
and new compiled cohorts/proofs remain referenced without relabeling revisions.

The acceptance is produced only after actual new public qualification, canonical
consumer proof validation, and fresh training route/alias exclusion succeed.
Hash and selection checks here do not repeat or replace native authentication.
A later published manifest or cache-inventory update is a separate validation;
it must never be substituted for the D40 identity used in these timings.

== Explicit projection of unused qualification training metadata

Without #raw("qualification_projection"), an accepted compatibility bridge keeps
the exact three-way #raw("immutable_selection") equality between producer,
final consumer workload and qualification workload. The narrow optional relation
#raw("unexecuted-training-prefix-pruning") instead records
#raw("checks.unchanged_selection: false"): the qualification metadata is different,
while the original producer and final two-root selection remain exactly equal.

The hashed relation has schema #raw("vakint.qualification-workload-projection.v1")
and references #raw("evidence"), #raw("independent_review"),
#raw("validator_scope") and #raw("incomplete_attempt"). It preserves the original
complete 42fa qualification and its thirteen successful public receipts, and
names f022 only as the compatibility target. The target's two training roots
must be the original producer's unchanged ordered list and the first two of the
qualification's three unused declarations. Public case IDs, roles, order,
powers, exact input strings, data, census and actual consumer cohort stay fixed.
The native worker never executes those training declarations during uncached
qualification; its reviewed source identity remains bound to the artifact.

The adapter rehashes all original complete result and guard payloads, checks
proof-derived raw/owner sets, twelve distinct nonempty held-out owner sets,
the full 51-key protected union and every selected nonzero direct owner endpoint.
It also requires the actual target OwnerProof validator result and strict cache
import/equality evidence. The validator checks supplied target metadata against
the original proof; its thirteen checks are not thirteen fresh reductions.
The separately failed target attempt remains ten successes, one deadline failure
and two unstarted cases, with no completed failed-case worker receipt or new
mapper proof. Matching partial output bytes never complete that attempt.

This relation is native compatibility evidence only. It transfers no timings,
selects no cache and grants no execution or performance acceptance. All source
worker identities remain attached to their original qualification hash. A final
performance sample must still carry its own actual target workload, successful
supervision, exact output/guards and accepted serialized window. Default
compatibility, strict cache completion and final-selection checks are unchanged.

== Selecting a completed pilot cache

An original #raw("phase: pilot") build remains a pilot build after a later
selection. To retain that completed pack, the existing evidence index must also
pin #raw("final_selection: {path, sha256}"), with schema
#raw("vakint.cache-final-selection.v1"). The receipt records
#raw("accepted: true"), #raw("decision: retain_complete_candidate"),
#raw("producer_phase: pilot"), an actual UTC #raw("decided_utc"), and a nonempty
#raw("reason"). It binds hashed #raw("cache"), #raw("build_report"),
#raw("producer_workload"), #raw("consumer_workload"),
#raw("candidate_acceptance") and #raw("cache_requalification") references.
The latter must be the same accepted bridge pinned by the evidence index.

The original candidate acceptance binds the same complete build, cache and
selected training roots. Every strict-completion, successful-supervision,
size, hash and failure check still applies. The accepted bridge checks the
original and current workloads, proofs, own cohorts, D40 data, census and Native
semantic binding, even if their revisions match. Qualification of an additional
route may be a conservative exclusion check; it never increases the selected
training count. Ordinary completed #raw("phase: final") builds retain their
existing admission path.

Root selection, consumer compatibility, execution release and final performance
acceptance are separate decisions. The adapter issues none of them. It retains
#raw("raw_report") unchanged, including the pilot phase, original source/proof
and preparation time, and stores the later selection separately. Selecting an
existing cache creates no new preparation or performance measurement.

== Connected control and warm reuse

The separate #raw("warm_reuse") supplement uses the existing evidence index
with case ID #raw("connected-five-loop-warm-reuse-control"). Each record is one
#raw("pair.json") from #raw("five_loop_saved_application --timed-pair-only"),
not two independently launched samples. The launcher pins hashed
#raw("input_manifest"), #raw("scientific_fixture") and #raw("source_snapshot")
references. The scientific fixture binds the accepted connected ladder-pinch
DOT, its original native-order powers and historical science receipt; it does
not transfer old timings or catalogue classification to the new run.

Final display requires exactly three fresh invocations with indices 0–2 and
variant #raw("cold-warm-pair") in one accepted serialized window. Actual
#raw("child_pid"), #raw("supervisor_pid"), #raw("child_start_ticks"),
start/end timestamps and observation scope must agree between launcher and
window. The first and warm durations remain nested within that process; no
inner timestamp is reconstructed. Rehash input, data/application, source and
all complete value/ordered-guard/family payloads. The original phase records
and all exact parity checks remain intact.

Report the first call and same-engine repeat separately. Their whole public API
timers include deferred owner preparation; loading/input preparation,
serialization, cache snapshots and parity checks are outside. The observed
attachment policy is prepared 1L–4L caches present and persistent 5L absent at
all three snapshots. This is neither a physical-rung cache speedup nor a
held-out cold median. Pure IBP, owner-only, terminal-only and catalogue closure
remain null. Unqualified completed pairs stay functional-only. Existing
primary/suite aggregates and their acceptance requirements are unchanged.

== Timing and rendering

=== Cumulative physical capture and fresh seen-input reuse

The final cumulative cache includes strict completed roots and intermediate
descendants produced while running the physical examples. The finite
#raw("build_cache") workload remains separately identified evidence; it is not
the sole production source. Preserve every durable prefix and immutable base.
There is no 50 MiB publication target or 95 MiB aggregate limit for native packs.
Record actual sizes and operational allowances; an older single-file format's
bound must not become a campaign-wide truncation policy. See
#link("../../cache_preparation.typ")[cache preparation and cumulative capture].

The physical capture producer and its fresh seen-input reuse have a separate
append-only report path:

```sh
python3 crates/vakint/docs/benchmark/merge_physical_capture_results.py \
  --base crates/vakint/docs/benchmark/data/results.json \
  --campaign /path/to/workspace/physical-capture-campaign \
  --output /path/to/workspace/physical-capture-results.json --select-current
```

The campaign must already contain the accepted producer and reuse completions,
their immutable configurations, full output/guard/unknown-master files, terminal
cache receipts and #raw("seen-reuse-comparison.json"). The adapter performs no
reduction, cache write, restart or eligibility decision. Missing completion or
failed flush/close is an error, not a partially populated timing row. While these
receipts are absent, the report remains explicitly pending; live counters and
source code are not substitute measurements.

The adapter rehashes the bound evidence and verifies identical executable,
input, data and settings, exact complete retained result and ordered guards,
family inventory and unknown-master inventory. It checks fresh-process
observations, successful completion, no newly captured reuse records, and
unchanged native snapshot generation and digest. Compiler/codegen receipts
retain the actual source snapshot and ELF identity; configured profile labels
alone do not establish optimization. Producer elapsed time includes capture
overhead. Pure IBP time remains unmeasured.

Append only under #raw("physical_capture_runs") and optionally update
#raw("selected_physical_capture_run"). Preserve every historical four-loop sample
and protocol, primary five-loop run, finite cache selection and held-out result.
This supplement reports *seen physical input*, separately from frozen held-out
sets. One matched fresh pair remains functional evidence, even with an optimized
executable: report raw observations and actual cache hits, without a median,
performance colour, final speedup claim or transfer into another benchmark.
The original comparison receipt is retained verbatim, including any diagnostic
ratio and its non-isolated scope. It does not relax the accepted serialized-window
and repeated-sample requirements above.

```sh
python3 -m unittest discover -s crates/vakint/docs/benchmark \
  -p 'test_merge_physical_capture_results.py'
```

=== Numerical summaries

Medians and full observed ranges are recomputed from actual per-case samples.
The coloured ratio is always *cached ÷ uncached* (four-loop cells remain
*RustRed ÷ FMFT*), so green means faster. Do not average different physical inputs
into an unlabeled “suite median”. Raw counters distinguish requested-root and
intermediate cache hits; neither extends campaign coverage or proves master closure.

Loading, cache import, fresh graph preparation, lowering, deferred owner
preparation, saved-application intervals, complete reduction and end-to-end time
retain their original boundaries. Events may overlap. Owner preparation is
already included in inclusive saved-application/reduction intervals and must not
be added again. #raw("ibp_excluding_owner_seconds") retains its historical field
name but means *saved application minus deferred-owner preparation*, including
transport, exact combination and terminal substitution. Pure IBP, terminal-only,
and primary owner-subtracted timings remain null. Missing events are null, not
zero. Primary E2E ends before final report/exit; suite E2E spans worker spawn
through exit. Representation reload is inside primary reduction; initial loading
uses Sparse. Preserve that setup caveat.

Render a workspace copy of #raw("benchmark.typ") beside the candidate
#raw("data/results.json") and inspect every page. Current pending 5L data must
still render as Pending. For layout tests, synthetic populated JSON/PDF must be
marked prominently as test fixtures and must never replace measured data.

```sh
python3 -m unittest discover -s crates/vakint/docs/benchmark \
  -p 'test_merge_five_loop_results.py'
typst compile /path/to/workspace/preview/benchmark.typ \
  /path/to/workspace/preview/benchmark.pdf
```
