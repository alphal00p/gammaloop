#set document(title: "Benchmark report maintenance")
#set page(paper: "a4", margin: 2cm)
#set text(font: "Libertinus Serif", size: 10pt)
#show raw: set text(font: "DejaVu Sans Mono", size: 8pt)

= Benchmark report maintenance

== Ordered-engine physical rung comparisons

Completed functional observations that cannot establish a qualified cohort are
kept under #raw("ordered_rung_diagnostics"), selected separately through
#raw("selected_ordered_rung_diagnostics"). The narrow
#raw("append_diagnostic") helper verifies report/result/completion hashes,
successful exits, MSbar conventions, catalogue closure and actual timing scope.
It does not run an equality check. Only the explicitly named passing functional
CI receipt can mark its saved pair and numerical checks passed; unpaired
catalogue-closed observations remain exact-comparison pending. These entries
always set performance eligibility false and both efficiency ratios null.
Their amber table keeps math, public-call and guard durations distinct, with
source cohorts and precision warnings preserved. No incomplete no-memo run is
converted into a timing or a completed three-arm benchmark.

The generic collector #raw("examples/ordered_reduction_rungs.rs") compares
ordinary DFS, DFS without completed-result memoization/bookkeeping, and Ordered
on the same authenticated physical scalar input. Prepared/saved caches and
previous-run result state remain unloaded in every arm. Only ordinary DFS enables
primary-engine completed-result memoization. Ordered may use a fresh private memo
inside opaque checked factorization; the no-memo DFS arm bypasses dependent
result memos too. *All new rung records use MSbar normalization.*
The adapter rejects new pySecDec records; earlier pySecDec fixtures, goldens
and measurements remain historical and are never relabeled. See #link("../ordered_reduction.typ")[the engine
guide] for the scientific fixture, strategy selection and progress semantics.

Through four loops, use #raw("--engine all") for exact and tracked-numerical
functional parity. Use separate #raw("--engine dfs"),
#raw("--engine dfs-no-memo") and #raw("--engine ordered") invocations,
with #raw("--no-monitor"), the same prepared input and normal release executable, for fresh-process
measurements. Do not promote same-process functional timings or historical
five-loop baselines into matched fresh-process speedup claims.
Five-loop physical runs use ordered reduction only, first zero momentum and
then quadratic, with the saved historical zero-momentum result as reference.
Every reduction invocation declares its actual #raw("--build-profile").

After separate workers finish, use the same collector with
#raw("--work-dir OUT --compare-reports DFS_DIRECTORY DFS_NO_MEMO_DIRECTORY ORDERED_DIRECTORY") to
authenticate the retained payloads and compare exact expressions and serialized
numerical centers/error radii. This cheap comparison performs no IBPs or
numerical reevaluation and writes the bound #raw("parity.json") consumed below.

The adapter #raw("merge_ordered_rung_results.py") consumes bound parity receipts:

```sh
python3 crates/vakint/docs/benchmark/merge_ordered_rung_results.py \
  --base crates/vakint/docs/benchmark/data/results.json \
  --id ordered-rungs-RECORD-ID \
  --parity /absolute/workspace/L1-q0-parity.json \
  --parity /absolute/workspace/L1-quadratic-parity.json \
  --output /absolute/workspace/benchmark-with-ordered.json \
  --select-current
```

It writes a *new* file, appending under #raw("ordered_rung_runs"). Selecting
the run sets #raw("selected_ordered_rung_run") without changing the historical
four-loop, five-loop, ingress or cache selections. Review the output before
replacing the canonical measured JSON and rendering the PDF. All 48 original
four-loop measurements and their original provenance remain untouched.
Synthetic adapter tests never supply published timing data.
To display balanced and asymmetric records from the same measured source
together, set #raw("selected_ordered_rung_runs") to their existing run IDs.
Their raw protocol and statistical counts remain separate; the renderer checks
matching source/Cargo identities and preserves each row's actual per-arm counts.
The singular selector remains the default when that list is absent. Measurements
from later normalization/materialization changes must not be inferred from an
earlier source cohort.

New parity receipts bind all three executed arms:

```json
{
  "schema": "vakint.ordered-rung-parity.v2",
  "input_sha256": "<prepared scalar input digest>",
  "reports": {
    "dfs": {"path": "<report.json>", "sha256": "<digest>"},
    "dfs_no_memo": {"path": "<no-memo report.json>", "sha256": "<digest>"},
    "ordered": {"path": "<report.json>", "sha256": "<digest>"}
  },
  "retained": {"dfs": "<result digest>", "dfs_no_memo": "<result digest>", "ordered": "<result digest>"},
  "intrinsic_memoization": {"dfs": true, "dfs_no_memo": false, "ordered": false},
  "exact_retained_equality": true,
  "numerical_enclosures_agree": null,
  "equality_scope": "<actual exact and numerical checks performed>"
}
```

Paths resolve relative to the parity file. Numerical equality is #raw("true")
when every report contains checked numerical references and #raw("null") when
none does. The adapter checks file hashes, payload sizes, input identities,
source/data/normalization agreement, completed catalogue closure, exact
dimension/evaluated-input bindings and cache exclusion. Numerical comparisons
share requested precision, effective finite source cap, evaluation point and
Laurent orders; equivalent exact expressions may expose different source-symbol
inventories, which are preserved separately. It does not perform
the symbolic comparison itself or certify the origin of an arbitrary receipt.
New raw reports retain #raw("vakint.ordered-rung-comparison.v1") and add an
explicit #raw("intrinsic_memoization") boolean, matching
#raw("cache_policy.depth_first_intrinsic_memoization"). The expected strategies
are #raw("DepthFirst"), #raw("DepthFirstWithoutMemoization") and #raw("Ordered").
The no-memo arm requires the actual dedicated walker, not a zero cache budget,
zero cache counters or a subtraction estimate. Receipt key #raw("dfs_no_memo")
maps to CLI/directory #raw("dfs-no-memo").
The boolean describes the primary engine's completed-result memo, not Ordered's
isolated auxiliary factorization implementation.
Existing two-arm version-one records remain as recorded and are not relabeled
or requalified as no-memo results.
New appended rows use #raw("vakint.ordered-rung-report.v2") and retain separate
Ordered/ordinary-DFS and Ordered/DFS-without-memoization ratios.

Without #raw("--protocol PATH"), measured times are diagnostic observations:
the report publishes no ordered/DFS speed ratio. A protocol supplied by the
measurement supervisor qualifies the matched fresh-process scope:

```json
{
  "schema": "vakint.ordered-rung-protocol.v2",
  "fresh_processes": true,
  "serialized_execution": true,
  "build_profile": "release",
  "profile_overrides": false,
  "cpu_affinity": "<fixed worker affinity>",
  "thread_environment": {"RAYON_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"},
  "executable": {"path": "<release ELF>", "sha256": "<digest>"},
  "build_receipt": {"path": "<closed build receipt>", "sha256": "<digest>"},
  "sources": {
    "main_revision": "<40 hexadecimal digits>",
    "main_source_tree_sha256": "<64 hexadecimal digits>",
    "native_revision": "<40 hexadecimal digits>",
    "native_source_tree_sha256": "<64 hexadecimal digits>",
    "kernel_revision": "<40 hexadecimal digits>",
    "cargo_lock_sha256": "<64 hexadecimal digits>"
  },
  "samples": [{
    "dfs": {
      "report_sha256": "<digest>", "pid": 123, "start_ticks": 456,
      "exit_code": 0, "executable_sha256": "<digest>",
      "launch": {"path": "<guard launch.json>", "sha256": "<digest>"},
      "completion": {"path": "<guard completion.json>", "sha256": "<digest>"}
    },
    "dfs_no_memo": {
      "report_sha256": "<digest>", "pid": 125, "start_ticks": 458,
      "exit_code": 0, "executable_sha256": "<digest>",
      "launch": {"path": "<guard launch.json>", "sha256": "<digest>"},
      "completion": {"path": "<guard completion.json>", "sha256": "<digest>"}
    },
    "ordered": {
      "report_sha256": "<digest>", "pid": 124, "start_ticks": 457,
      "exit_code": 0, "executable_sha256": "<digest>",
      "launch": {"path": "<guard launch.json>", "sha256": "<digest>"},
      "completion": {"path": "<guard completion.json>", "sha256": "<digest>"}
    }
  }]
}
```

The numeric process identifiers above illustrate a file schema, not real
measurements. Every supplied sample must bind a distinct successful original
worker, the same release executable and the matching producer identity.
The adapter checks bound #raw("run_guarded.py") launch/completion receipts for
the original worker identity, successful unsignaled exit, matching engine/case/
profile command, affinity, threads and executable, and non-overlapping UTC
windows. The supervisor owns the actual build and source-snapshot verification.
The closed build receipt has the following contract:

```json
{
  "schema": "vakint.ordered-rung-build.v1",
  "success": true,
  "build_profile": "release",
  "profile_overrides": false,
  "command": ["cargo", "build", "--locked", "--release", "<other actual arguments>"],
  "build_log": {"path": "<actual build log>", "sha256": "<digest>"},
  "executable": {"path": "<same ELF>", "sha256": "<same digest>"},
  "sources": {"<same source fields as the protocol>": "<same identities>"}
}
```

The executable and source objects must match the protocol exactly. Record every
actual argument separately in #raw("command"), including any local configuration
used by a diagnostic build; qualified final builds use coherent published pins.
The adapter rejects invalid supplied protocols instead of silently downgrading
them to diagnostic records. DEV reports with explicit profile and nullable
self-hash remain admissible without a performance protocol.
Current qualification uses *one* actual ordinary-DFS run, *one* actual
no-memo DFS run, and three Ordered runs. Earlier balanced version-two cohorts
retain their three original samples per arm unchanged. The asymmetric contract
changes only the protocol/report version; exact comparison still uses the
existing three-directory comparator and version-two parity receipts:

```json
{
  "schema": "vakint.ordered-rung-protocol.v3",
  "sample_counts": {"dfs": 1, "dfs_no_memo": 1, "ordered": 3},
  "samples": [
    {"dfs": "<first actual worker>", "dfs_no_memo": "<first actual worker>", "ordered": "<first actual worker>"},
    {"ordered": "<second actual worker>"},
    {"ordered": "<third actual worker>"}
  ]
}
```

The worker objects use the same guarded file contract above; all other build,
source, cache and timing requirements remain unchanged. This is one case per
asymmetric protocol. Run the cheap comparator three times, pointing to the
*same immutable* DFS/no-memo report directories and each distinct Ordered
directory. Repeated baseline references are equality evidence, not additional
measurements. The adapter deduplicates measured samples by report hash, requires
exactly one/one/three unique workers, and writes
#raw("vakint.ordered-rung-report.v3") with per-arm #raw("sample_counts").
It never repeats baseline numbers to create fake medians or ranges.

The table labels both baselines #raw("n = 1 · single measurement"); Ordered
retains its measured median/range. Ratios use each single baseline observation
and carry no estimate of baseline run-to-run variability. All five workers run
serially. Each arm's boolean policy stays fixed. No new physical five-loop DFS
run is included. Existing version-one/two cohorts remain as recorded, including
their original sample counts and qualification.
Release status is not inferred from #raw("debug_assertions=false").

The headline is #raw("reduction_math_seconds"), using fixed scope
#raw("vakint.rustred.reduction-math.v1"): non-overlapping IBP/routing and exact
weighted terminal-coefficient assembly/grouping spans, including required
terminal substitutions, projections, lower-factor reductions and request-local
reducer/dispatch state construction in both engines. Immutable artifact loading,
deserialization/program import are excluded, along with graph/data preparation,
matching, initial scalar lowering, MATAD normalization,
master numerical evaluation and serialization. This is not a pure Native-kernel
IBP timer. Gross math time minus measured lazy owner-preparation wall time gives
the checked net. The collector fixes owner preparation to one worker; the
subtracted native wall-time counters must not be interpreted as CPU time.
The gross is a semantic sum of disjoint blocks, not a contiguous interval ending
at normalization. Later exact combination and natural request-state cleanup
inside the call remain counted; normalization work and engine-lifetime teardown
after the call remain excluded. No forced cache clear is added.

The report contract preserves unavailable measurements explicitly:

```json
{
  "reduction_math_seconds": null,
  "reduction_math_timing": {
    "scope_id": "vakint.rustred.reduction-math.v1",
    "complete": false,
    "gross_seconds": null,
    "measured_owner_preparation_seconds": null,
    "interval_count": 0,
    "incompleteness_reason": "math capture feature disabled"
  }
}
```

With complete telemetry, gross, preparation and net are measured numbers; net
must equal the checked difference, with valid non-overlapping interval scope.
An invalid preparation delta produces incomplete telemetry and a null net, not
a physics failure. Optional capture requires #raw("experimental-rustred").
The adapter keeps gross/deduction/count metadata, rejects scope or deduction
inconsistency, and requires complete math telemetry for every qualified sample.
Missing measurements stay null: no subtraction of unrelated public stages,
no partial-sample median, and no fallback to the whole public timer.
Qualified headline comparisons require #raw("monitoring_enabled=false") in
both reports and #raw("--no-monitor") in each bound launch command. The command
must not enable terminal or snapshot observation. Monitoring and logging inside
math spans contribute real wall time; no observer-overhead subtraction is
invented. Monitored functional/dashboard runs, including five-loop examples,
remain diagnostic and cannot supply a qualified ratio. Independent OS/RAM
guarding is unaffected.

Record #raw("worker_policy.owner_preparation_workers=1") and the actual
#raw("cache_policy"): prepared caches not loaded, previous-run caches not
restored, explicit intrinsic DFS memoization policy, and ordered result cache
unsupported. This is measurement metadata, not a new cache toggle. Keep the
runtime policy unchanged unless the user explicitly chooses otherwise.

Retain #raw("retained_reduction_seconds") separately as *public retained
evaluation*, including frontend lowering and finalization. Loading is separate;
#raw("completed_seconds") starts after scalar preparation and
includes subsequent numerical reference generation when enabled. Whole-process
peak RSS includes earlier stages, including previous engines in a same-process
functional pair. Unavailable pure-IBP, isolated routing/lowering, terminal
substitution and ordered-progress observations remain #raw("null"), never zero.
Phases refer to one native request and are not a wall-time percentage.

The benchmark keeps the existing large colored cells. Its ordered multiplier is
*ordered net reduction math ÷ DFS net reduction math*; original RustRed/FMFT
multipliers retain their original meaning.
Three-arm Ordered cells print both ratios. Their color compares Ordered with
the faster DFS variant, so green means it beats both baselines; a win over only
the no-memo arm cannot hide a regression against ordinary DFS.
An independent historical 40-hour five-loop baseline remains labeled with its
original source, cache and input cohort. It cannot qualify through this matched
new-report protocol or silently become a new uncached sample.
Its reduction-math time remains unknown; the historical 40-hour public-call
duration must never populate that field.

For a completed new ordered zero-momentum run, use
#raw("--historical-parity PATH") to append a separate historical comparison.
This option may be supplied without #raw("--parity"); it never qualifies for
a fresh-process protocol or a speed ratio. Its owner-produced receipt binds:

```json
{
  "schema": "vakint.ordered-rung-historical-parity.v1",
  "exact_retained_equality": true,
  "input_sha256": "<identical authenticated scalar input>",
  "historical": {
    "acceptance": {"path": "<old parent acceptance>", "sha256": "<digest>"},
    "report": {"path": "<old benchmark.json>", "sha256": "<digest>"},
    "retained": {"path": "<old result.symbolica>", "sha256": "<digest>"},
    "definition": {"path": "<legacy-definition-identity.json>", "sha256": "<digest>"}
  },
  "ordered_report": {"path": "<completed new report.json>", "sha256": "<digest>"},
  "normalization_conversion": {
    "schema": "vakint.ordered-rung-normalization-conversion.v1",
    "source_normalization": "pySecDec", "target_normalization": "MSbar",
    "api": "VakintSettings::get_integral_normalization_factor_atom",
    "ratio_identity": "N_target/N_source", "loop_count": 5,
    "epsilon_symbol": "<configured epsilon>", "mu_r_sq_symbol": "<configured scale>",
    "dimension_binding": {"<same fields as the new report>": "<same values>"},
    "factor": {"canonical_expression": "<exact API-derived Atom>", "sha256": "<UTF8 digest>"},
    "source_retained_sha256": "<old raw digest>",
    "target_report_sha256": "<new report digest>",
    "converted_retained": {"path": "<separate converted payload>", "bytes": 123, "sha256": "<digest>"},
    "exact_conversion_verified": true,
    "applied_before_laurent_truncation": true
  },
  "retained": {
    "historical": "<old raw digest>", "historical_converted": "<converted digest>",
    "ordered": "<new MSbar result digest>"
  },
  "equality_scope": "<actual exact comparison performed>"
}
```

The adapter checks the old report/result against its accepted hashes, byte
lengths, catalogue closure and recorded public-call duration, then checks the
new ordered report's input, projection and MSbar normalization. It additionally
binds the historical pySecDec definition, API-derived exact target/source factor,
dimension binding and separate converted historical payload. Apply the full
factor before any Laurent truncation; never compare differently normalized raw
results or alter the original old payload. The existing normalization API owns
all factors, phases and scale dependence; the adapter performs no handwritten
normalization algebra. The exact
comparator still owns mathematical equality; merely constructing this receipt
does not perform it. No incomplete ordered result can be published through this
path. The accepted 55-master public-pipeline reference remains distinct from
the earlier 63-master Native-only result. The historical baseline used running
cache capture and concurrent engineering jobs; it is retained as one historical
observation, not a new matched no-cache sample. The adapter never opens native
cache payloads or the large historical guard export.

Run the lightweight adapter contract checks from the benchmark directory with
#raw("python3 -m unittest test_merge_ordered_rung_results.py"). These checks
use synthetic temporary files and do not run reductions.

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

Both #raw("measure_cache") and #raw("five_loop_qcd_ladder") accept explicit
single-file grants through #raw("--cache-max-compressed-bytes N"),
#raw("--cache-max-decoded-bytes N"), and #raw("--cache-max-index-bytes N").
Omitting them retains the public API's compact defaults. Set grants from the
authenticated cache inventory and available resources; a refusal never raises
them automatically or retries a reduction. Decoded bytes include the native
header, while the index grant covers retained native indices rather than total
process RAM. These flags require #raw("--cache") and do not apply to directory
packs. Record-count diagnostics are not a cache sizing control.

For a large independently hosted file, also pass
#raw("--cache-decoded-directory PATH") with an existing workspace directory.
The native reader decompresses into private owned disk storage and shares one
lazy symbolic-state import across its reader forks. The collector resolves the
directory to an absolute path before loading or spawning workers. Resolved byte
grants and directory are recorded under #raw("cache_file_options"); suite workers
receive the same policy, and resume rejects a different one. Uncached comparison
samples retain this configured policy in their identity while attaching no cache.
Cache loading remains outside the reduction timer.

Python exposes the same authenticated single-file attachment through the frozen
#raw("CacheFileIdentity") record. Obtain the identity document from an
independently authenticated manifest or a trusted native export report; parsing
the document alone does not authenticate its origin. The four fields use the
native schema: compressed and decoded byte lengths plus their 32-byte BLAKE3
digests. #raw("to_document()") preserves its JSON byte arrays, while the read-only
digest properties are Python #raw("bytes"). For example, with the authenticated
document and an existing decoded directory:

```python
from symbolica.community.hepkit.vakint import CacheFileIdentity, CacheFileOptions

identity = CacheFileIdentity.from_document(authenticated_identity_document)
options = CacheFileOptions(
    max_compressed_bytes=identity.compressed_bytes,
    max_decoded_bytes=identity.decoded_bytes,
    decoded_directory=decoded_directory,
)
engine.load_cache(5, cache_path, file_options=options, file_identity=identity)
```

The native loader validates both content digests and lengths before attachment.
The identity does not increase resource grants: set any larger index/storage
grants explicitly from the inventory and available resources. Omission retains
the existing trusted-generated-file API; official external-data loading retains
its own authenticated identity handling. Single-file identities are rejected for
lower-loop packs before their paths are read.

Python #raw("evaluate"), #raw("evaluate_rustred_symbolic"), and
#raw("prepare_cache") expose #raw("cache_hit_diagnostic=False"). Setting it
to true selects the existing lossless native observer; DEBUG verbosity alone
does not enable that observer. It requires #raw("experimental-rustred") and
cannot be combined with terminal presentation or progress snapshots. The host
initializes the existing Rust-to-Python logger once and preserves a logger
already owned by an embedding application. Configure Python handlers and DEBUG
levels before the first native message because the bridge caches enablement:

```python
import logging
from symbolica import core

core.use_custom_logger()
logging.getLogger("vakint").setLevel(logging.DEBUG)
# Attach a caller-owned logging.Handler before evaluating.
value, guards, families = engine.evaluate_rustred_symbolic(
    integral, cache_hit_diagnostic=True,
)
```

The native #raw("vakint::cache_workload") and
#raw("vakint::benchmark_cache_hits") targets become Python
#raw("vakint.cache_workload") and #raw("vakint.benchmark_cache_hits"). Their
messages retain the existing JSON documents. Retain and reconcile every event;
a latest monitor snapshot or aggregate hit counter cannot establish per-seed
Route coverage. The flag performs no numerical reevaluation or cache capture.

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
  --physical-scope /path/to/workspace/physical-scope.json \
  --output /path/to/workspace/new-results.json --select-current
```

Both existing report adapters accept the optional #raw("--physical-scope PATH")
artifact. Its schema is #raw("vakint.benchmark-physical-scope.v1"); the latest user
override requires only the zero-momentum diagnostic and defers the quadratic
fixture:

```json
{
  "schema": "vakint.benchmark-physical-scope.v1",
  "physical_scope": {
    "required": ["zero-momentum-diagnostic"],
    "deferred": ["quadratic-four-rung-uv"]
  },
  "reason": "User requested zero-momentum-only physical completion; quadratic remains deferred."
}
```

The shared validator checks the partition and nonempty reason. Each appended run
records #raw("physical_scope") plus #raw("physical_scope_authorization"), retaining
the artifact's path, SHA-256 and raw receipt. Reuse the same artifact for both
adapters. Without it, the original requirement for both physical examples remains.
Deferred primary rows remain visible as *Deferred*, with null timings, closure and
parity; their original #raw("functional_records") stay intact. The override changes
which physical completions are required, not the evidence needed to measure a
required input. Preserve all 48 historical four-loop measurements and the frozen
twelve-input held-out suite.

The output must not already exist. Omit #raw("--select-current") to append while
retaining the old selection. Completed strict cache preparation can be included
with #raw("--build build.json --cache cache.gz"); completed suites additionally
use #raw("--held-out suite.json") and/or #raw("--tuning suite.json"). Suites require
their exact cache build and proof/cohort bindings. Duplicate evidence runs fail.
No historical run is rewritten. Evidence is appended under #raw("five_loop_runs");
the explicit selection fields choose which runs are displayed. Review the new JSON and rendered
PDF before replacing the canonical trio together. Do not update canonical data
while a measurement or source-freeze window is active.

=== Single-observation engineering comparisons

The physical adapter also accepts #raw("--engineering INDEX.json") instead of
#raw("--campaign"). The #raw("vakint.q0-engineering-evidence-index.v1") index contains
exactly #raw("before"), #raw("after_off") and #raw("after_on") observations. Each binds
hashed benchmark, execution and verification JSON; the two new observations also
bind their complete coefficient-proof JSON. The adapter checks successful actual
waits, source/input/result/guard/cache identities, all scalar product proofs and
raw stage samples without reading cache or symbolic payloads or running reductions.
OFF/ON means cache-hit diagnostics disabled/enabled; both may use the same cache.

This appends a separate #raw("engineering_optimization_runs") record and optional
#raw("selected_engineering_optimization_run"). It never changes physical/final suite
selections, acceptance scope, historical samples or performance medians. Retain
raw receipts and source hashes with each observation. Ratios use matching public
or prepared-input-through-serialization fields; supervised startup/exit wall
time is distinct. Missing phases and pure IBP remain unmeasured, never zero.
Different cohorts and background work qualify these as engineering observations,
not final matched cache comparisons. For example:

```sh
python3 crates/vakint/docs/benchmark/merge_physical_capture_results.py \
  --base crates/vakint/docs/benchmark/data/results.json \
  --engineering crates/vakint/docs/benchmark/provenance/q0_engineering_evidence_index.json \
  --output /path/to/workspace/new-results.json --select-current
```

The additional #raw("--file-engineering INDEX.json") mode imports completed
#raw("vakint.cache-file-engineering-evidence-index.v1") metadata: a full cold FILE
handoff, LIGHT utility report, and accepted Route/Apply export-budget receipts.
It resolves only bounded hashed JSON records, including native count/reimport
proofs and the original tuning baseline; it never hashes or opens referenced
cache/result payloads. The separate #raw("selected_cache_file_engineering_run")
displays this supplement without replacing the prior engineering or final-suite
selection. Raw receipts, stage arrays, closed captures and source cohorts remain
under #raw("engineering_optimization_runs").

Keep FILE loading/authentication, collector whole-file SHA scans and controller
startup/close boundaries explicit. A shorter public call can accompany slower
cold end-to-end behavior, and lower rule counts need not produce lower runtime.
Concurrent engineering observations do not establish a causal format advantage,
a three-run median or a LIGHT utility winner. The 95,000,000-byte LIGHT cap and
complete native export/reimport checks remain required; neither export counts nor
typed exporter identities establish typed application attachment or performance.

The #raw("--delivery-checkpoint INDEX.json") mode appends hashed user decisions,
the completed monitor acceptance, and the closed HEAVY source inventory without
altering historical measurements or their selections. Its
#raw("vakint.q0-delivery-evidence-index.v1") index binds bounded JSON records for
#raw("light_freeze"), #raw("authorization"), #raw("monitor"), and
#raw("heavy_sources"), plus the existing historical physical and cache-file
engineering run IDs. It validates the accepted LIGHT bytes against that export,
keeps the historical uncached cohort separate, and leaves final HEAVY export and
published-download fields null until actual receipts exist. A user waiver of
new repetitions never creates a matched median or a cache-only speedup claim.
Optional paired #raw("publication") and #raw("release_assets") receipts record
the actual source push and immutable-before-alias publication. They must bind
the same source revision and successful deterministic pack preflight. Publication
alone leaves fresh-download and offline-reuse validation pending.
The optional #raw("published_consumers") receipt requires those publication
bindings, the frozen LIGHT identity, all five cumulative loops, fresh immutable
and latest downloads, both fresh-process offline reuses, cache opt-out and all
640 uncomputed master getters. Zero TCP connections apply only to each official
loader window, rather than the whole licensed process.
Optional #raw("heldout_correctness") records one completed strict control and
fresh replay separately from timing suites. Validate full-value bytes, public
canonical guard sets and used-master/family inventories, retain their limited
comparison scope, and require every new nonempty capture in the sealed final
HEAVY source manifest. This evidence never changes frozen LIGHT selection.
Final HEAVY evidence is an all-or-none set of #raw("heavy_census"),
#raw("heavy_export"), #raw("heavy_reimport"),
#raw("heavy_reimport_controller") and #raw("heavy_reimport_terminals"). The
adapter uses the existing Native phase/controller/PIDFD schemas, matching the
closed source manifest, all reference visits, typed file identity and exhaustive
coefficient/original-ordered-guard equality. It does not parse cache payloads or
reconstruct Native records. Cache-finalization durations retain their original
phase scope and never become reduction timings.

The subsequent single cached final-file q0 observation supplies another complete
set: #raw("final_q0_benchmark"), #raw("final_q0_execution"),
#raw("final_q0_verification"), #raw("final_q0_coefficient_proof"),
#raw("final_q0_monitor"), #raw("final_q0_terminals") and
#raw("final_q0_capture"). It requires accepted final HEAVY evidence, actual
original process terminals, complete exact scalar-product/public-guard proofs,
verified raw stage arrays and a clean saved capture. Unexpected new captured
records require preservation and a new final-HEAVY disposition before delivery.
Suite collectors require an explicit null #raw("reduction_error"); the physical
collector may omit that field only with its separately verified successful
execution receipt. An explicitly reported error is never ignored.
Display these actual observations as #raw("n = 1"), retaining the historical cold
reference, separate source/cache cohorts, null pure-IBP/hash timings and the
explicit waiver. Do not leave completed science labeled as pending matched
benchmarks; matched medians and suite timings that were not obtained remain
explicitly unclaimed. Append a new checkpoint/index; preserve its earlier
pending predecessor and all historical four-loop samples.

```sh
python3 crates/vakint/docs/benchmark/merge_physical_capture_results.py \
  --base crates/vakint/docs/benchmark/data/results.json \
  --delivery-checkpoint crates/vakint/docs/benchmark/provenance/q0_delivery_evidence_index.json \
  --output /path/to/workspace/new-results.json --select-current
```

=== Display both cache tiers

Use #raw("--select-tier heavy") for an authenticated cumulative native cache, or
#raw("--select-tier light") for the actual selected gzip and its strict export
proof. The adapter records #raw("cache_tier") on each run and updates only that
tier's entry in #raw("selected_five_loop_runs: {heavy: runID, light: runID}"). It
checks the actual cache representation; the flag does not relabel a finite
training cache or create a missing measurement.

```sh
python3 crates/vakint/docs/benchmark/merge_five_loop_results.py \
  --base /path/to/workspace/heavy-results.json \
  --build /path/to/workspace/strict-light-export-validation.json \
  --cache /path/to/immutable/light.rrcache.gz \
  --tuning /path/to/workspace/light-tuning.json \
  --held-out /path/to/workspace/light-held-out.json \
  --evidence /path/to/workspace/light-evidence-index.json \
  --physical-scope /path/to/workspace/physical-scope.json \
  --select-tier light --output /path/to/workspace/both-tier-results.json
```

The input above already selects a genuine HEAVY run. Adding LIGHT preserves that
selection, the legacy singular selector and all 48 four-loop observations.
#raw("--select-current") continues to control #raw("selected_five_loop_run"). When
the plural selector is absent, the report uses that existing single-run view.
When present, each tier gets its own uncached/cached comparison and full timing
ledger. Each ratio uses that tier's matched baseline. The same actual direct q0
uncached references may serve both comparisons only when the existing validation
accepts their identical source, data, input, limits and execution window; never
copy or restamp a baseline. Suite baselines keep their tier-specific identities.

LIGHT uses the genuine #raw("vakint.strict-light-cache-export-validation.v1")
receipt as #raw("cache.raw_report"). Its compressed bytes and SHA-256 come from
#raw("output.cache"), and its selected unique-record count from #raw("selection").
The original validated cumulative cache stays under #raw("cache.source_producer");
its source data and index identity are provenance, not the LIGHT consumer dataset.
The exporter keeps its own compiled cohort. Export, native-validation,
selection-admission and metadata hashes remain separately visible. Reporting the
authenticated LIGHT file does not require reopening the full HEAVY payload.

The LIGHT ledger has no invented finite-workload root counts. Its optional raw
Native #raw("phase_seconds") map records finite, nonnegative monotonic durations
for completed named intervals. Each named phase runs at most once per invocation. The adapter
retains that map under #raw("cache.export_validation.phase_seconds") and projects
the following #raw("cache.production_seconds") fields:

- #raw("selection"): null; selection happens outside this helper.
- #raw("preparation"): #raw("load_saved_application"), when measured.
- #raw("export"): #raw("encode_selected_completed_cache") plus
  #raw("write_selected_cache"), only when *both* phases are present.
- #raw("reimport"): #raw("reimport_selected_cache"), when measured.
- #raw("comparison"): #raw("compare_selected_native_visitor"), when measured.

An absent phase or map leaves its projected cost null; render null as “—”. The
helper does not invent #raw("selection") or #raw("output_publication") phases.
Every duration belongs to the raw native-validation receipt's own producer cohort,
even when the export wrapper or consumer has a different cohort. Never infer cost
from file size, record count, a validation flag or a heavy producer's duration. Actual
data/cache loading and reduction costs belong to fresh consumer samples,
separately. Adding this optional timing contract does not establish that an
export has run or that any phase has been measured.

=== Cumulative cache provenance

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
physical examples; #raw("physical_examples_complete") requires every input in the
recorded #raw("physical_scope.required"). An absent scope retains the original
requirement for both the zero-momentum diagnostic and complete quadratic fixture.
Deferred examples are labelled *Deferred*, never measured or silently omitted.
The separate #raw("all_source_runs_completed") also accounts for retained incomplete prefixes.
Snapshot preservation alone has #raw("performance_eligible: false") and never
supplies missing timings or pretends to be a completed finite cache workload.

A generic completed tuning capture uses #raw("kind: running-capture") and
#raw("scope: completed-capture"). Its producer remains the original measurement;
separate hashed completion evidence joins the accepted functional result to the
exact embedded capture, successful flush/close, zero pending work, Native
snapshot and frozen index. The snapshot barrier binds that completion and records
#raw("producer_completed: true"). The accepted legacy tuning receipt binds its
actual worker and terminal; the CLI monitor chain additionally binds its
postflight, result and retained-capture evidence. Original provisional and
starting-key-authentication flags remain unchanged. This scope contributes only
to #raw("all_source_runs_completed"), never to completed physical examples or
scientific/performance qualification. A failed or unfinished source remains a
#raw("committed-prefix").

#raw("compose_cumulative_cache_report.py --union-validation ACTUAL_RECEIPT")
composes an existing ordered Native opaque-union receipt with repeated
#raw("--source") descriptors. It requires actual per-source Native censuses and
frozen snapshots, exact source order, bindings, copied/reference arithmetic and
the producer's own compiled artifact and resource admission. Completion receipts
alone cannot supply those Native preservation checks. All original-source
validators run before the new report directory is published; rejected evidence
leaves no accepted report behind. The single-source strict-preparation command
retains its existing behavior.
An existing #raw("native-index-prefix") source joins its immutable binary index
to the Native source snapshot and the historical proof's authenticated per-role
census, rather than the earlier union's aggregate counts. Its failed strict
producer remains incomplete and its scope remains #raw("committed-prefix").
Every historical captured role must still be represented exactly once. Omitting
an earlier q0 prefix requires separate actual same-writer containment authority;
a later generation, larger inventory or matching path does not establish it.

The union command also requires #raw("--data-manifest PATH") and
#raw("--proof-root PATH"). The manifest's authenticated bytes supply report data
context, not a claim that the standalone Native producer loaded that dataset or
linked Vakint. They do not alter its compiled cohort or rebase original source
data roots. The proof root resolves the actual producer's artifact references.
Current consumer, application-binding, import, exposure and performance joins
remain separate requirements. A preservation receipt is not final HEAVY/LIGHT
delivery, an application evaluation, a distinct-key census or source eligibility.

Final cumulative-cache suites additionally require an independently accepted
#raw("vakint.cumulative-cache-consumer-requalification.v1") receipt. It binds the
original source starting-key inventories and route/alias exclusion against the
unchanged frozen owner sets, explicitly classifying exposed starting inputs.
The #raw("captured-tuning-workload") inventory kind binds an original provisional
measurement to its hashed workload, exact tuning case and input, complete raw
lowered-key list and current Native owner proof. Its original cohort and
#raw("native_keys_authenticated: false") flag remain recorded; postcombine seeds cannot replace raw
source roots or remove an overlap with protected held-out inputs.
An optional hashed #raw("source_owner_proof") authenticates an exposed historical
tuning case independently of the current eligible tuning selection. It must use
the full current compiled cohort and artifacts, bind its own successful uncached
worker and qualification workload, and preserve the captured source's exact input
and raw roots. This source authority does not add a known-seen case to the final
tuning or held-out sets; current route and protected-owner exclusion checks still apply.
Compatible multi-source unions preserve each source's original data-manifest
identity rather than rebasing it to the union report's data manifest.
The existing #raw("physical-census") kind requires an actual
#raw("vakint.qcd-ladder.campaign-census.v1") source. Pathless input metadata resolves
through the configuration's hashed prepared manifest and its established
#raw("vakint-input.symbolica") payload, with matching projection, hash and byte
length. An older selected-witness census with the same input/family does not
establish the active public matcher's dispatch. Its historical schema remains
unsupported here; actual producer/witness equivalence is still required.
Its #raw("cache_import") reference must bind the successful public-import
observation for the same descriptor, native indexes, data and kernel; import
validation remains separate from attachment and reduction measurements.
The existing exact result/guard, binding/cohort, twelve held-out sets, three fresh
matched samples and serialized-window requirements remain in force. Useful
incidental descendant and root hits stay visible. The full heavy collection has
no aggregate size cap. A separate five-loop light cache must fit within
95,000,000 compressed bytes and is the only tier published in external data.
Neither tier rewrites the original four-loop history.

Publishing the selected cache changes release-manifest hashes. A consumer may
therefore supply #raw("data_metadata_requalification: {path, sha256}") pointing to
a #raw("vakint.external-data-metadata-requalification.v1") relation. It binds the
old and new manifests and the frozen workload, owner proof and collector artifacts.
The relation selects exactly one representation: the existing
#raw("runtime_cache_snapshot") descriptor addition, or
#raw("prepared_cache_replacement: {loop_count: 5, export_validation: {path, sha256}}")
for the strict light file. A light replacement requires no selected runtime
snapshot that could override that prepared file.

Only the existing 5L prepared-cache file and #raw("cache.json") may be replaced;
the inline native identity and decoded checksum must match the independently
validated export. Binding, format and paths remain exact. The compressed cache
must contain 1–95,000,000 bytes. All other loop-pack payloads, scientific fields,
1L–4L caches and version-bound root provenance remain unchanged, apart from the
existing release-version/archive metadata changes.

The light acceptance binds the frozen heavy source, selected keys, eligible usage
profile, actual exporter/reader artifacts, strict reimport counts and complete
record/ordered-guard equality. Native validation owns record interpretation;
the external source validator owns selection eligibility. Source proofs live
under repository-root #raw("provenance/") and are not added to 5L runtime packs.
Their portable references in #raw("cache.json") resolve from that root; individual
source files remain below 100 MiB, with compressed key/profile summaries permitted.
The raw native validation does not reference #raw("cache.json"), avoiding a hash
cycle. The top acceptance separately binds both metadata and native/source proofs.

Source captures and import evidence retain their original identities. Fresh
consumer proofs must belong to the new dataset and their own authenticated
artifacts. This relation never retargets old receipts, qualifies light timings from
heavy evidence, or permits mixing datasets within a fresh measurement comparison.

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
Keep that full heavy collection without an aggregate size limit and deliver its
compressed local path for independent hosting. Build a separate light cache of
at most 95,000,000 compressed bytes for the external-data repository and release.
Use measured completed-record hits, input diversity and demonstrated benefit to
select records; aggregate monitor snapshots do not establish per-record frequency.
Keep final held-out results out of selection and report loading, reduction and
end-to-end measurements for both immutable tiers. Record actual sizes and
operational allowances; trimming never deletes the heavy cache or its prefixes. See
#link("../../cache_preparation.typ")[cache preparation and cumulative capture].

The physical capture producer and its fresh seen-input reuse have a separate
append-only report path:

```sh
python3 crates/vakint/docs/benchmark/merge_physical_capture_results.py \
  --base crates/vakint/docs/benchmark/data/results.json \
  --campaign /path/to/workspace/physical-capture-campaign \
  --physical-scope /path/to/workspace/physical-scope.json \
  --output /path/to/workspace/physical-capture-results.json --select-current
```

The scope artifact uses the same validation and hashed provenance described above.
The zero-momentum-only override leaves the quadratic example deferred and retains
its historical evidence; it does not qualify a missing physical completion.

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

The #raw("measure_cache") collector has no per-worker deadline by default.
An explicit positive #raw("--max-case-seconds N") enables a diagnostic budget
over the whole fresh worker, including data/cache loading and result/guard
serialization. Only that opt-in arms its watchdog and permits deadline-based
termination. It is independent of the fast closure tests' 300-second ceiling.
Without a completed measurement receipt, an interrupted unbounded attempt stays
pending: a later supervisor refuses to start a duplicate or infer completion
from elapsed time. Preserve it for independent terminal review. A bounded
interrupted attempt retains the existing expired-attempt failure without retry.
Final comparisons keep the complete frozen suite and a declared serialized
window; a slow case is never silently dropped or retried.

Medians and full observed ranges are recomputed from actual per-case samples.
The coloured ratio is always *cached ÷ uncached* (four-loop cells remain
*RustRed ÷ FMFT*), so green means faster. Do not average different physical inputs
into an unlabeled “suite median”. Preserve raw cache counters and the collector's
#raw("backing_reuse") label as provenance; derive the displayed scope separately.
In the audited Native #raw("4190f2c6f315") implementation,
#raw("requested_root_hits") recognizes requested *Route* records only.
#raw("intermediate_hits") also includes requested-owner *Apply* hits, so a zero
root counter does not prove intermediate-only or proper-descendant-only reuse.
A positive root counter proves some exact requested-root reuse; other-hit scope
remains unclassified. Report this conservatively unless separate authenticated
initial-frontier evidence establishes the finer scope. A revision string, claimed
label or unvalidated additional field does not establish that evidence. Neither
counter extends campaign coverage or proves master closure. Incidental cached
descendants remain distinct from explicitly selected starting training integrals
and their direct owner aliases; this display qualification does not change the
frozen held-out exclusions.

The optional #raw("measure_cache --classify-cache-hits --repetitions 1")
diagnostic requires #raw("experimental-rustred") and a read-only #raw("--cache").
The flag explicitly sets #raw("RustRedMonitorOptions.cache_hit_diagnostic");
its default is false. Registering a log sink alone never enables this observer.
It runs one cached worker per frozen case, without a presenter or running
capture. Its interval-zero observer records exact surviving postcombine roots
and validated backing-hit identities, with complete per-dispatch reconciliation.
The existing #raw("authenticate_owner_roots") mapper binds that successful
receipt to an accepted same-cohort uncached receipt, follows only authenticated
Route-to-Route/Apply transitions, and joins the saved hits offline. Apply-rule
children are excluded. No second reduction is needed to create the certificate.
Missing events, an incomplete routing frontier or an unknown transport leave
classification unavailable. The versioned recursive certificate supplements
the original direct-owner proof; claimed flags never upgrade historical proofs.
Diagnostic timings include observation overhead and supply no performance
medians, final comparison or training-selection authority. The full collector
guard files retain their canonical sorted/deduplicated policy; they are not
evidence of the Native internal condition-vector order.

The physical #raw("five_loop_qcd_ladder --classify-cache-hits") diagnostic
uses that same lossless logger and reconciliation with one retained evaluation
of an authenticated prepared input. It requires #raw("--cache") and
#raw("experimental-rustred"), forbids fresh preparation, warm repetitions,
running capture and presentation, and supplies no performance claim. Preserve
its full trace, source/helper hashes and immutable cache identities before/after.
Its reference must be a separately supervised classifier-off cached replay of
the same final executable, prepared input, scientific definition, data and cache.
Join both actual launch/completion records and exact retained result, family and
full canonical guard-file hashes. The public guard file includes preparation
conditions and is sorted/deduplicated; Native cache record guards retain their
own original order. Historical physical results remain independent references,
not substitutes for that same-cohort pair. The external
#raw("vakint.physical-cache-usage-receipt.v1") adapter authenticates these
source-specific records before reusing the existing per-dispatch classification
contract; it does not create a synthetic Workload or widen campaign coverage.

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

A previously completed finite preparation enters a cumulative cache as
#raw("kind: strict-preparation, scope: completed-preparation"). Keep its original
producer data identity and cohort, accepted output, preserved native indexes,
and a separately authenticated full Native/GZIP source-parity run. That run
must validate every source record and its ordered guards under the current
owner; a selected subset cannot qualify the entire source. Its actual successful
guarded supervision and own compiler artifact remain required evidence.
Finite preparation completion never completes a physical example or qualifies
final held-out exclusions, cache selection, publication, or performance.

#raw("compose_cumulative_cache_report.py") composes this existing evidence into
the ordinary cumulative report and validation receipt. It checks the actual
read-only source indexes but neither copies chunks nor imports coefficients.
The source JSON names the seven hashed documents #raw("producer"),
#raw("acceptance"), #raw("preservation"), #raw("native_validation"),
#raw("native_supervision"), #raw("native_artifacts"), and #raw("inventory"),
beside its unique ID, kind, and scope. Use a new output directory; preserve
failed attempts and original receipts. An unchanged source index/inventory
allows an already validated selected export to retain its original identity.

The #link("provenance/producer_sources.typ")[producer source supplement] preserves
the exact private exporter, build recorders, admission helpers and controller
examples used for the qualified optimized pilot. Its checksummed source archive
and index live beside that note. These are development source records; relocation
requires new paths, source hashes, compiler identity and guard registration.
The archive contains no cache payload or credentials and makes no fresh-host
rebuild or final-release claim.

```sh
python3 crates/vakint/docs/benchmark/compose_cumulative_cache_report.py \
  --source /path/to/workspace/provenance/source.json \
  --cache /path/to/workspace/preserved-native-bundle \
  --output-directory /path/to/workspace/provenance/new-cumulative-report
```

```sh
python3 -m unittest discover -s crates/vakint/docs/benchmark \
  -p 'test_merge_five_loop_results.py'
typst compile /path/to/workspace/preview/benchmark.typ \
  /path/to/workspace/preview/benchmark.pdf
```
