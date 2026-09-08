# CI measurements, 8 September 2026

Faster completion of tests, Clippy and doctests remains unproven. Per-crate cache
boundaries and test inventories are preserved, and local merged artifact output
sizes fell about 75%, but size alone does not meet the acceptance criterion.
The first four remote suites had cache-upload timeouts and interrupted workers.
An unchanged baseline reported failure while work remained queued, then continued
executing checks for another hour. Repeated identical builds and replaced logs
make that run unsuitable as a speed baseline. Validation has resumed: the final
ungrouped main candidate is the sixth submitted suite, within the twelve-suite
plan and fourteen-suite maximum. No performance change has been accepted for merge.

[Machine-readable observations](ci-measurement-results.json) record revisions,
derivations, clocks, limitations and remote log links. The
[measurement procedure](ci-measurement.md) defines reproduction and acceptance.

## Local compilation

| Layout | Baseline | Candidate | Merged artifact NAR sum | Selected tests |
|---|---:|---:|---:|---:|
| Main | 42m34s | 38m38s | 87.44 → 21.92 GiB | 1,633 |
| FeynKit | 45m57s | 36m54s | 126.00 → 31.59 GiB | 1,834 |

Separate project-cold stores were seeded only with external outputs. Both
variants used cores=8, max-jobs=4 and sandboxing on the same shared host.
These advisory settings are not hard CPU isolation; baseline compression used
unrestricted threads. Elapsed times are observations, not guaranteed speedups.
NAR sums count distinct store outputs, not network bytes, worker minutes or CHF.
These timings use the frozen candidate revisions recorded in JSON, before
grouping removal and later correctness/reporting follow-ups. The final scheduler
has identity/graph validation. The compiled main roots retain their frozen
identities; a remote run of the final configuration is now underway.

The 11 main / 12 FeynKit roots cover archive groups, Python packaging, Clippy,
formatting and graph validation. Licensed execution and doctests are
excluded from those cold timings and are measured separately below. All 17/26 archive inventories match names, filters and ignored flags.
All 58/91 artifact contexts have identical observed compilation messages
(738/827 each), with zero compilation at the final archive stage. Both subsequent
revision-only runs performed zero Cargo compilation.

Seven independent FeynKit edits passed targeted builds and exact repeats:
GammaLoop, Python, CFF, model, embedded fixture, manifest and lockfile. Every
repeat started zero builders. Independent source copies were built sequentially
in the populated FeynKit candidate store; later cases can reuse earlier outputs.
Rust/manifest/lockfile probes append comments, and the fixture adds whitespace.
These compile without running tests and check invalidation/reuse, rather than
providing independent cold or downstream-consumer timings. All five unchanged source/fixture sentinels
retained their hashes and references; manifest and lockfile edits invalidate
broadly. The full derivation matrix preserves all nine extracted production
crates for a GammaLoop edit and all eight Python siblings for a Python edit.

Python retains an inherited inefficiency: its source belongs to its own
dependency bundle, so a leaf edit repeats ancestor compilation inside that
bundle. Its 13 dependency compilation messages match the frozen baseline and
candidate. Separating those context-specific dependencies is a future experiment;
this change preserves existing features and anchor commands.

One root-only encoding sample illustrates the limit of the NAR metric. Using
identical xz-6 settings, the expanded baseline encoded to 557.62 MiB while the
candidate encoded to 775.16 MiB, despite raw NAR size falling from 3.411 to
0.767 GiB. Encoding took 170.858s versus 48.512s. The candidate already contains
a zstd archive. Its registered closure also shrank from 19.782 to 1.758 GiB,
retaining one merged bundle instead of nine. Neither this encoding proxy nor the
closure size establishes provider network traffic or faster test completion.

## Scheduling and remote evidence

The shared final-archive/Python producer experiment was removed. Ready groups
waited for remaining builds and uploads, and intermediate uploads continued.
The proposal selects 55 main / 75 FeynKit build outputs with synchronous
discovery and explicit edges, plus runtime checks and service stages. Removing the bundle preserves all compiled artifact identities.

Five-minute low-speed upload retries occurred in the
[main candidate](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency/bf1910fdc8e185e312a1451d7ce21730571627f5/db2344e9-2865-4487-b852-d41fca6ed636),
[FeynKit baseline](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-benchmark-feynkit-baseline/e4f0a1838f922df270f05d5d858d3950cc05a3be/5cac03d6-4989-4f7e-93b6-ed0dc5c7c3e5),
and [FeynKit candidate](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency-feynkit/8a07fa8ad6be4733ec537d6ac4de8eae8d154e52/5d55885a-8255-4c8a-9858-ffca53cf554e).
The underlying cause is unknown. Main's producer eventually succeeded after
1h57m07s, with 1h45m04s of completed upload intervals overlapping other work.
FeynKit's producer failed at two hours. Upload success messages do not establish
visibility to later workers. Three abandoned workers also recovered under the
same URLs with overwritten logs; saved incident annotations keep those pairs
ineligible for improvement percentages.

The completed main grouped suite took 3h15m05s to finish required checks,
compared with 1h09m55s for its baseline. Its final success followed 9m50s later; its final worker ran for about 19s. FeynKit baseline Python was
observed queued after its displayed prerequisites had passed, with gaps of
48m05s and 10m54s. The later FeynKit candidate attempt reported 15 HTTP/2 cache
download errors; its completed restores covered about 6.00 GiB in 49.923s, with
zero compilation messages. These are symptoms with unknown causes, recorded
separately from compiler work. Repeated job phases are distinguished from retries.

The unchanged [main baseline](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-benchmark-main-baseline/d44bf2ed781d73840175b349d2b7f15acf1a7589)
failed after its GammaLoop dependency job became
[hopeless](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-benchmark-main-baseline/d44bf2ed781d73840175b349d2b7f15acf1a7589/0a535a40-ea11-4343-acb7-dd547f530050).
At 18:30:17 UTC it had 31 successful, 19 cached, one hopeless and 28 queued jobs,
with no active workers. Doctest passed 43 cases, Clinnet five and Vakint 74; five
other runtime groups had not run at that observation. The failed job's last
exposed log contains only
the initial build command, with no diagnostic explaining the failure.

Later saved worker and GitHub completion clocks establish real execution after
that failed-suite snapshot: the core group passed 814 tests at 19:07 UTC and the
integration group passed 104 at 19:29 UTC. Their workers compiled two and 54
prerequisite units respectively before executing archived tests. By 20:13 UTC
the suite had 42 successful, 27 cached, one hopeless and nine queued jobs. The
18:30 metrics above remain historical lower bounds, not a final resource total.

The first actual application worker log began 55m58s after submission,
54m08s after evaluation. Sixty-nine sampled states preserve readiness and waiting
observations. Another [phase-fix suite](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fraised_energy_cff_wip_optimized_phase_fix/32f27225460fbf5a2c9ad507355a8b8423fd237a)
had active workers during this period and finished successfully at 18:15:48.
This overlap may confound timing; the queue cause and worker limits are unknown.

Eighteen verified log replacements across 11 job URLs, plus an observed WASM
abandonment, leave resource history incomplete. Current streams contain at least
52.67 worker minutes, 16.75 GiB of reported downloaded content and 164 compilation
messages. Earlier streams are preserved separately to avoid double counting.
These are lower bounds, not billed cost or network bytes. No transfer-timeout
warning was observed in this fifth run's saved streams.

The [initial doctest](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-benchmark-main-baseline/7198fb48ea2bb5fe4baea4f904944b35bf857c89/13557053-744d-40be-939a-04c880e036d9)
and [unchanged doctest](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-benchmark-main-baseline/d44bf2ed781d73840175b349d2b7f15acf1a7589/d7963413-0a33-4c65-96a4-f80a3627bc9c)
produced the exact same full derivation and output paths. Both passed 43 cases
and compiled the same 73 units. The earlier cache-copy hook completed successfully,
yet the later run rebuilt the output. Cargo preparation took 6m54s then 13m07s;
total worker spans were 8m19s then 14m42s. Those are distinct clocks, and the
slower compilation has no established cause. Exact result identity rules out a
source, feature or license-input identity change for this pair without inspecting
credentials. Three raw package producers—Clinnet, Kurvst and tracing-filter
macros—also rebuilt exact derivations that earlier logs reported as uploaded.
The JSON includes paired job links, identities and completed-upload timestamps.

Our shared doctest Cargo archive does omit additional compilation contexts, so
it compiles when the result check actually executes. This inherited repository
limitation is separate from failure to reuse the identical completed Nix result.
Adding another large producer is not justified by this evidence. Cache visibility,
substitution behavior and replaced worker histories need diagnosis first.
The final ungrouped main candidate is now
[running on NixCI](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency/9f7b4392be6d8694a845315631505714fcc26b1b).
Independent FeynKit source scenarios are undergoing local full-consumer trials
before submission. No Actions workflow was manually dispatched.

Both frozen Nix Actions workflows passed:
[main, 2h10m52s](https://github.com/alphal00p/gammaloop/actions/runs/34227037456)
and [FeynKit, 2h04m32s](https://github.com/alphal00p/gammaloop/actions/runs/34227238762).
They include packaging and other work beyond the primary NixCI selection, so
these are not equivalent-coverage speed comparisons. Their doctest jobs passed
43/62 cases and logged 73/81 compilation messages separately.

## Resumed local validation

For the same integration-test input context, exporting its complete registered
closure to a local binary cache reduced raw content from 19.782 to 1.758 GiB
and encoded content from 3.874 to 1.149 GiB. Both variants used xz-6 and the same
hard eight-CPU affinity. Export took 1,185.048s versus 130.102s; two fresh-store
imports took 203.693s / 202.824s versus 18.228s / 15.424s. This includes the
older merged bundles retained by the baseline, which the root-only sample above
excluded. The local restore improvement is real for this context; it does not
measure WAN transfers or establish faster completion of the whole test suite.

Both variants pass all selected Rust/Python tests, doctests and Clippy locally,
using the supplied license. The full comparison includes an initial phase and
a recovery phase:

| Layout | Final selected tests / doctests, each variant | Integration recovery, baseline → candidate | Exact full cached repeat |
|---|---:|---:|---:|
| Main | 1,633 / 43 | 176.864s → 175.002s | 0.253s → 0.263s |
| FeynKit | 1,834 / 62 | 154.516s → 159.827s | 0.276s → 0.306s |

Fresh required-check results used cached build inputs in separate stores, with
identical eight-CPU affinity and four concurrent builders within each layout.
Both implementations hit the same two existing integration-test timeouts while
Clippy and doctests compiled. Rerunning the identical required roots after those
compiler checks finished rebuilt only the failed integration group; all tests
then passed, with zero Cargo compilation. This supports contention as a possible
cause, without establishing it uniquely. No test, timeout or filter was changed.

All 11 main / 12 FeynKit required outputs are registered as successful. Exact
full repeats started zero builders and executed no tests, demonstrating local
result reuse. Final coverage uses the latest successful result for each group
and does not add duplicate executions. Initial failed timings remain excluded
from speedup claims. Nextest runtime phases compiled nothing; Clippy/doctest
compilation-message multisets match between variants: 45/73 on main, 44/81 on
FeynKit. Individual completion observations use Nix registration, not output
directory appearance. Local result reuse does not prove remote cache visibility.

## Remaining acceptance

The boolean-only license probe remains an optional benchmark diagnostic; its
extra Python CI gate was removed before the final-configuration comparison.
Its unlicensed negative control and six historical error/redaction cases pass.
Historical green jobs alone do not establish unrestricted license status;
compare equivalent license settings without adding a new ordinary test gate.

The current public Linux CLI builds and loads its native libraries. Strict
help/version assertions fail: both commands exit 1 after printing expected text
through unchanged Clap error handling. Baseline parser/dependency sources are
identical, but its binary was not retested. No application change or relaxed
assertion is included. Main native Actions also exposes a preexisting Clippy
warning under Rust 1.98; pinned Nix Clippy passes.

Before merging, complete equivalent runtime validation and the remaining remote
unchanged/source comparisons for the final configuration. Faster completion of
tests, Clippy and doctests is the primary requirement. Verify worker minutes,
restores, repeated builds, individual-group latency and the final-success tail.
The 80% restore, 50% warm-worker and 30% changed-worker targets remain goals.
