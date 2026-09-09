# Measuring CI changes

See the [recorded experiment results](ci-measurement-results.md) for current
measurements and outstanding acceptance gates.

The [shared test-feature experiment](ci-test-feature-sharing.md) has separate
local validation; its effect is not included in the earlier NixCI comparisons.

Run `just ci-report MANIFEST OUTPUT_DIR` to collect existing NixCI jobs and
GitHub Actions attempts. The underlying command is
`node .github/scripts/ci-report.mjs MANIFEST OUTPUT_DIR`. It uses Node built-ins
and the authenticated `gh` CLI, and never launches, retries, cancels, or changes
CI. The just recipe selects the repository's pinned Node runtime.

## Manifest

Each entry identifies one exact application layout, CI implementation, scenario,
and NixCI suite. Use full commit SHAs and the commit-level suite URL, without a
trailing job UUID. For example, this reads an existing historical suite:

```json
{
  "repository": "alphal00p/gammaloop",
  "runs": [
    {
      "layout": "historical-raised-energy",
      "variant": "baseline",
      "scenario": "warm",
      "sha": "91142139e1fde44cc4709a35ad090631d27ca451",
      "suiteUrl": "https://nix-ci.com/gh:alphal00p:gammaloop/raised_energy_cff_wip/91142139e1fde44cc4709a35ad090631d27ca451"
    }
  ]
}
```

The labels describe the experiment; the collector does not infer them from
branch names. Replace the example with the actual measured layout and scenarios.
`variant` is `baseline` or `candidate`. Pairs share `layout`, `scenario`, and
an optional `pair` label (default `"1"`). Give a repeat a different pair label.
Duplicate suite URLs and ambiguous pairs are rejected.

By default the collector retrieves all Actions runs for the SHA and all their
attempts. Optional `actionRunIds: [34063777721]` limits the selected workflow
runs, retaining every attempt of each selected run. There is no workflow
consolidation. GitHub check runs use `filter=all` and are joined to NixCI jobs by
exact details URL; another job attempt on the same SHA cannot supply its clock.

Required-check latency defaults to the latest attempt of each NixCI test plus
`gammaloop-clippy`, `gammaloop-fmt`, and `gammaloop-guppy-workspace-graph`.
Packaging, documentation, and WASM checks are excluded from this primary set. For a
controlled comparison of different scheduling layouts, provide the same
explicit `requiredAttributes` array in both entries. Include every intended
test group, Clippy, doctest, formatting, and graph/configuration check. Missing
required attributes leave required latency unknown and make the suite incomplete
and ineligible for paired percentages. Requested, observed, and missing
attributes are recorded separately. Required `nix-ci-check-*` groups must have
their test phase; the package build alone cannot satisfy them. The final success
job is measured separately.

If NixCI recovers an abandoned worker under the same job URL, its API and logs
may replace the earlier attempt. Preserve the incident when first observed and
add an explicit annotation to that suite's manifest entry:

```json
{
  "observedInterruptions": [
    {
      "jobUrl": "https://nix-ci.com/gh:owner:repo/branch/FULL_SHA/JOB_UUID",
      "observedAt": "2026-09-08T13:16:29Z",
      "evidenceFile": "incidents/abandoned-worker.json"
    }
  ]
}
```

The evidence JSON must contain the same `jobUrl` and `status: "abandoned"`.
Keep its earlier observed metrics, such as `workerSeconds` and
`downloadReportedBytes`, and an optional `file` pointing to the preserved raw
log. Unknown measurements remain null. Log paths resolve relative to that
evidence file; evidence paths resolve relative to the manifest.

The report copies the incident JSON and linked raw log into
`raw/N/interruptions/`. Earlier measurements remain in
`observedInterruptions`, separately from current-log resource totals, to avoid
double counting streams that might overlap. The annotation makes the suite
incomplete even after a successful retry. Final-only collection cannot detect
interruption history that the service has already overwritten.

A log stream can also be replaced while its status stays `started` and its job
URL and check ID stay the same. Record these separately as `observedLogReplacements`, using the
same `jobUrl`, `observedAt`, and `evidenceFile` annotation fields. The evidence
JSON must repeat that exact URL and observation time, with `file` pointing to
the prior NDJSON and `replacement.file` to the later NDJSON. Both paths resolve
relative to the evidence JSON. No abandoned status is required or inferred.

The collector parses both saved streams and compares their records. Identical
or append-only streams do not establish replacement; changed worker starts or
loss/replacement of the prior record prefix do. It copies both raw logs and
recomputed separate metrics into `raw/N/log-replacements/` and emits a
`log-replacement` incident in JSON, CSV, and Markdown. A shortened stream proves
lost visible history, not worker termination. The cause remains unknown, current
resource totals become lower bounds, and the suite stays ineligible for paired
percentages even if all required checks later succeed. Earlier stream totals
are never added to current totals because their overlap is unresolved. Missing
or invalid requested evidence also makes the suite incomplete.

Optional submission clocks come from the local push call, never the commit's
author timestamp. Set `submittedAt` immediately before calling `git push` and,
when available, `submissionCompletedAt` immediately after it returns. The report
keeps suite/configuration latency and separately measures
`submissionToRequiredSeconds` through the last required check. The two recorded
push bounds also produce `submissionDurationSeconds` and a range for elapsed
time after server receipt; the exact receipt time remains unknown.

To inspect earlier polling evidence, add `observationFiles` containing saved
snapshot paths to the suite's manifest entry. This reuses the existing snapshot
shape, including files saved as `raw/N/snapshot.json` by live collection:

```json
{
  "at": "2026-09-08T15:59:32Z",
  "suites": [
    {
      "spec": { "sha": "FULL_SHA", "suiteUrl": "EXACT_SUITE_URL" },
      "suite": { "commit": "FULL_SHA", "runs": [] }
    }
  ]
}
```

Retain the actual suite API response in `suite`; the empty runs array above is
only a shape example. For queue observations, also set
`dependencySnapshot: { "sha": "FULL_SHA", "file": "generated-config.json" }`.
The file must contain the generated `dependencies` mapping for that exact
revision. The report records a queued job only when every declared prerequisite
has exactly one successful/cached build row in the same snapshot. Missing
dependency evidence, a running prerequisite, duplicate prerequisite rows, or a pre-worker
gap alone cannot establish this incident. Hidden dependencies and the scheduling
cause remain unknown; separate snapshots do not prove continuous queue time.

Snapshots also retain visible restarts under one job URL and distinct URLs for
the same type and attribute. A normal build/test pair is not a repeated attempt.
The trigger for a repeated attempt remains unknown. Replaced or absent attempt history makes resource totals lower bounds
and prevents accepted paired percentages. This collector never retries jobs or
runs an unbounded watcher.

## Outputs and replay

The output directory contains:

- `manifest.json`: the exact input manifest, including explicit check expectations.
- `report.json`: all metrics, per-job evidence, Actions jobs/steps, paired
  comparisons, and collection errors.
- `report.md`: readable suite, Actions-attempt, and paired tables.
- `suites.csv`, `jobs.csv`, `actions.csv`: aggregate and individual clocks.
- `transfers.csv`: completed transfers with reported sizes, durations, and
  timestamps.
- `compilations.csv`: compiled crate descriptions with job/context identity.
- `repeated-derivations.csv`: exact derivation paths built in multiple jobs.
- `incidents.json`, `incidents.csv`: timestamped symptoms, job/log links, saved
  evidence links, and known versus unknown causes; also included in report JSON
  and its Markdown incident table.
- `pairs.csv`: baseline/candidate percentage changes. Negative is an
  improvement; a zero baseline gives an undefined percentage.
- `raw/N/`: suite, check, and Actions JSON plus fetched NDJSON logs and a job
  index, in manifest order.

New raw evidence files are created with owner-only file permissions. Review raw logs
before publishing them; the summary does not include arbitrary log messages.
Existing output files are replaced when collecting again. Use a separate
directory to retain a previous snapshot.

To replay a snapshot, add this object to its manifest entry:

```json
{
  "offline": {
    "suite": "saved/raw/1/suite.json",
    "checks": "saved/raw/1/checks.json",
    "jobs": "saved/raw/1/jobs.json",
    "actions": "saved/raw/1/actions.json"
  }
}
```

Offline mode makes no network calls. These paths are relative to the manifest.
When replaying interruption annotations, point each `evidenceFile` at its copied
`raw/N/interruptions/N.json`; its linked log path is relative and self-contained.
For log replacement replay, point each `evidenceFile` at its copied
`raw/N/log-replacements/N.json`; both linked raw logs are self-contained.
For observation replay, point `observationFiles` at the copied
`raw/N/observations/N.json` snapshots and `dependencySnapshot.file` at the copied
`raw/N/dependency-snapshot.json`. These retain the original snapshot schema.
The job index is an array of objects containing `url` and `file`; absolute
log paths are accepted, and relative log paths resolve against the job index.
This also accepts the earlier scratch collector's `manifest.json` or
`jobs.json` index directly. Check evidence may be an array of check objects,
one GitHub response object, or an array of paginated response objects.
Actions evidence is an array of REST workflow-attempt objects with a `jobs`
array attached; use `[]` when no matching workflow ran.

A missing, empty, or malformed worker log produces an explicit error and leaves
that job's measurements unknown. Cached jobs need no worker log. Skipped or
unstarted jobs are not counted as missing workers. The command preserves its
partial report and exits **2** for incomplete evidence; invalid invocation or
manifest exits **1**. Resource totals in an incomplete report are observed
lower bounds. Abandoned workers retain all visible measurements, but their logs
may omit trailing work: `interruptedJobs` counts them and
`observedResourceLowerBound` remains true even if a later retry succeeds. Such a
recovered suite cannot supply paired performance percentages. A queued, running,
or cancelled suite is also incomplete, even when every currently available log
has been collected. A pair with incomplete evidence, unmet required checks, or
an unsuccessful suite has no performance percentage comparison.

Run the offline regression checks with:

```sh
node --test .github/scripts/ci-report.test.mjs
```

## Metric definitions

**Suite latency** runs from the configuration check's start to the last
completed actual NixCI job. Required latency ends at the last intended check.
Final aggregate latency ends at the completed deploy/success check; final tail
is the gap after required checks. A queued deploy placeholder is not a
completion timestamp. Failed suites retain diagnostic elapsed time but are not
successful timing samples. Cancelled/skipped suites have no completed duration;
abandoned job timestamps cannot extend completed-job latency.

Each job also records completion measured from the suite's configuration start.
Paired group comparisons use that completion latency, preserving their own
check durations separately so delayed scheduling remains visible.

**Worker minutes** sum each observed log's first-to-last worker interval,
preferring its monotonic relative nanoseconds. Check duration, pre-worker gap,
and post-worker gap are separate fields. The pre-worker gap can include
dependency waiting, scheduling, or setup; the available logs do not isolate
queue time. The observed span can omit work before the first or after the last
record. Summing parallel workers measures resource occupancy, not suite latency
or billed CPU. No CHF estimate is inferred.

**Transfers** count only completed `Downloaded cached … (size) in duration`
and `Uploaded … in duration` records. Announcements do not establish completed
transfers. Each completion timestamp and duration defines an interval. Union
overlapping intervals within each worker, then sum across workers. Download
and upload times can overlap, so a combined transfer union is also reported.
Upload sizes remain unknown if the logs omit them. Reported artifact sizes are
not established compressed wire bytes. Timers may include network,
decompression, store import, and locking; they do not distinguish those phases.

`transferTimeouts` records low-speed upload/download warnings separately, with
reported thresholds and retry attempt numbers. Request URLs and arbitrary error
text are omitted. `transferTimeoutCount` also appears in job/suite totals; it
does not invent durations or byte counts for unfinished attempts. Use these
events when assessing whether a pair had comparable service conditions.

`transferErrors` additionally records HTTP/2 stream/framing errors and other
failed cache transfers, retaining direction, HTTP status and numeric retry
fields. An HTTP status of 200 does not override a logged framing failure. Neither
these failures nor timeout warnings count as completed transfers. Incident text
uses fixed descriptions and omits request URLs and arbitrary error messages;
the underlying backend/network cause is unknown.

For a job currently reported as failed or abandoned, the incident timestamp
prefers its recorded check completion. If absent, it uses the earliest saved
snapshot reporting that exact job URL and status; otherwise it is null. `timeSource` records
`check-completion`, `status-snapshot`, or null. A snapshot timestamp records when
the status was observed, not precise worker termination. The last worker log
record is never used to infer failure time.

The incident table also flags service-reported job failures/interruption, exact
derivations built under multiple job URLs, multiple check attempts sharing one
URL, missing/reversed clocks, and a successful final aggregate completing more
than 60 seconds after required checks. Clock defects retain available evidence
but invalidate timing comparisons. A final-success tail is a measured gap; it
does not identify why finalization took that time.

`intermediateDownloadReportedBytes` totals content restored by jobs classified
as artifact producers from their attribute names: per-crate dependency/test
artifacts, Cargo-artifact roots, prebuild, `ci-test-inputs`, and test-binary archives. This is an
operational proxy for intermediate traffic, not an identification of every
downloaded path's role. Inspect `transfers.csv` when producer boundaries change;
moving the same transfers to a differently named job is not a saving. Total
reported restores and worker minutes are retained alongside this subset.

**Compilation** reports Cargo `Compiling` and `Checking` messages separately by
job/context and the
`Finished … in …` command wall times. These are not per-crate CPU profiles.
A repeated crate name can represent different features, target kinds, profiles,
or compiler modes. Only the identical full `.drv` path appearing in build-start
messages across different job URLs counts as repeated identical Nix work.
Repeated messages within one job do not count as multiple workers. Job URLs
identify worker invocations, not physical machines.

**Test evidence** is `executed` when test result/summary lines are visible,
`reused` when the service marks the test cached or the actual result derivation
is visibly substituted, and `unknown` otherwise. Substituting the small
`nix-ci-check-*` wrapper does not establish successful-result reuse. A reused
successful result is valid coverage when its inputs and expected inventory
match. Summary counts and visible test IDs are retained, but are not a complete
inventory: skipped tests and successful cached runs may omit names. Therefore
`coverageEquivalent` stays null in paired output. Validate test inventories,
filters, features, Python behavior, licenses, and doctest/Clippy coverage
separately before accepting performance. Nextest summary `failed` includes
ordinary failures and timed-out tests; `timedOut` records the timeout subset.
Do not add these two fields together. Individual timeout rows retain `TIMEOUT`.

License mode is a benchmark observation, not an additional Python CI gate.
Use comparable license settings in both variants; existing tests retain their
license handling. If explicit confirmation is needed, the optional diagnostic
below uses the existing Python environment and packaged extension. Replace the
example module path with the actual built extension path. This diagnostic is
separate from ordinary test execution and its elapsed time must be recorded.

```bash
# Report the packaged runtime's license mode without exposing library logs.
if python3 - "/absolute/path/to/gammaloop/_gammaloop.so" >/dev/null 2>&1 <<'PYTHON'
import ctypes
import sys

library = ctypes.CDLL(sys.argv[1])
check = library.is_licensed
check.argtypes = []
check.restype = ctypes.c_bool
# Distinguish an unlicensed result from interpreter or library failures.
raise SystemExit(0 if check() else 2)
PYTHON
then
  echo '{"symbolica_runtime_licensed": true}'
else
  case $? in
    2) echo '{"symbolica_runtime_licensed": false}' ;;
    *)
      echo '{"symbolica_runtime_licensed": null}'
      echo "Symbolica license status probe failed" >&2
      ;;
  esac
  exit 1
fi
```

True means the pinned runtime accepts unrestricted mode; false cannot distinguish
expiry, an unknown key, or server failure. A probe failure reports null.
Timestamped keys may validate locally before background registration, so true
does not establish a fresh synchronous server check or the state of every later
test process. No credential values or library logs are printed.

**Actions** retain all attempts separately. Creation-to-update duration matches
the historical comparison, but GitHub keeps the original `created_at` on a
rerun. Use attempt-start-to-last-job for rerun latency. Worker minutes use
observable job start/completion, including completed cancelled workers as
resource consumption. Cancelled/skipped attempts have no completed workflow
timing. Job runner labels and step durations remain available in JSON; elapsed
time alone does not establish equal hardware or cache state.

## Controlled experiment and acceptance

Keep every existing Actions workflow automatic. Compare CI implementations
within each application layout; main versus FeynKit is not a causal measurement
of crate splitting. Use these twelve NixCI suites:

| Layout | Scenarios, each baseline and candidate | Suites |
|---|---|---:|
| Main | Initial; unchanged follow-up commit | 4 |
| FeynKit integration | Initial; unchanged follow-up; GammaLoop-only source edit; FeynKit CFF source edit | 8 |

Start each source-edit scenario independently from its clean application
revision. Preserve code, toolchain/lockfiles, features, license availability, and
resource limits within each pair. Record the exact revisions and suite URLs
before comparing. Do not call the first remote run cold: shared caches may
already contain its inputs. Allow one ambiguous pair to repeat, for **14 suites
maximum**; the collector itself launches none.

In addition, run project-cold builds locally in disposable stores with identical
limits for both variants of both layouts. Inspect derivations and use targeted
builds to verify GammaLoop-only, Python-binding, CFF, model, embedded-fixture,
and manifest/lockfile invalidation. Preserve unaffected per-crate artifacts and
successful result caching. Whole-workspace Clippy/doctest costs remain separate
from package test compilation.

Acceptance requires equivalent coverage first. Warm unchanged runs require
zero real Cargo recompilation, targeting 80% less intermediate restored content
and 50% fewer worker minutes. Changed-source runs target 30% fewer worker
minutes without repeated identical test-artifact builds or invalidation of
unaffected packages. Project-cold completion must be no more than 10% slower.
Inspect individual-group latency so grouping does not delay independent cached
checks, and require final success within 60 seconds of the last intended check.
Performance targets are measurements to establish, not promised savings.

The initial candidate trial grouped final test archives and the Python module
behind one producer, retaining internal per-crate compilation. Its acceptance
remains unproven: ready groups waited behind that shared prerequisite, and NixCI
still uploaded intermediate outputs. Cache transfer timeouts affected baseline
and candidate jobs, so their stalled timings cannot isolate a grouping effect.
The proposed configuration therefore schedules archive groups independently;
the grouping trial remains in the recorded benchmark history. Its raw compiled
artifacts are unchanged by removing the scheduling bundle.

The compact-archive implementation was also checked with a tiny local Cargo
round trip: both dependency and application remained fresh after restoration,
with epoch-1 timestamps. Its smaller archive validates the representation, not
production speed or transfer savings. See the
[artifact reuse audit](architecture/nix-crane-cache-reuse.md) for that boundary.
