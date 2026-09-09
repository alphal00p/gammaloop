# CI measurements, 9 September 2026

The completed unchanged-main pair finished required checks in **8m27s instead of
37m57s**, with **84% fewer observed worker minutes** and **93% less reported
intermediate restored content**. Both variants passed all 11 required checks.
This is one pair on shared infrastructure, not a guaranteed improvement or a
CHF estimate. The zero-recompilation goal remains unmet.

The proposal keeps independent per-crate compilation, narrows automatic NixCI
work to the required checks and their producers, and publishes self-contained
compressed Cargo merges. The shared final producer experiment was removed.
A separate fix preserves sibling crate fingerprints; its full FeynKit validation
passes. Eleven of at most fourteen experimental NixCI suites have been submitted.
The fixed FeynKit remote run finished with a retry; source comparisons remain pending;
no merge acceptance has been declared.

[JSON observations](ci-measurement-results.json) retain exact revisions, clocks,
limitations and job links, including the earlier unsuccessful experiments.
[The measurement procedure](ci-measurement.md) defines reproduction and acceptance.

## Unchanged main on NixCI

| Measurement | Baseline | Candidate |
|---|---:|---:|
| Required checks | 37m57s | 8m27s |
| Whole suite | 38m31s | 8m54s |
| Observed worker minutes | 145.85 | 23.63 |
| Reported intermediate restored content | 70.93 GiB | 4.69 GiB |
| Cargo compilation messages | 404 | 103 |
| Derivations repeated across workers in the suite | 26 | 0 |
| Final success after last required check | 34s | 27s |

The [baseline](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-benchmark-main-baseline/64cdbfff2dffcae199d9ea41d1e5a3cd1debee4a)
and [candidate](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency/7e6dc57327fbd6ab3e6181151d66963e8e93675c)
have complete observed histories, without replaced logs or transfer errors.
Clinnet passed with insufficient log detail to classify execution versus reuse;
its status must remain unknown. Coverage parity is independently checked locally.
The baseline's precise submission clock is unavailable, so the table measures
from the first check. Candidate submission to required completion was 8m36s.

The candidate's 103 compilation messages comprise 73 in doctests, 16 in the
Linnet runtime job's preparation, and 14 across individual artifact producers.
An identical completed doctest result and previously published Kurvst artifacts
were rebuilt across commits. Those are separate from within-suite duplication.
These observations meet the restore, worker-time and final-tail goals for this
pair, but fail the zero-Cargo goal.

## Local cold and source-change builds

| Layout | Project-cold baseline | Candidate | Change | Selected tests |
|---|---:|---:|---:|---:|
| Main | 63m00s | 54m23s | −13.66% | 1,633 |
| FeynKit | 56m05s | 52m48s | −5.88% | 1,834 |

Each variant used a fresh store seeded only with external outputs, hard affinity
to eight CPUs, cores=8, max-jobs=4, sandboxing, no remote project cache and no
upload hook. Both pairs pass exact archive-inventory and compiler-message parity;
exact repeats start zero builders. Memory, storage and other host workloads were
shared. These are single local observations. They cover archive groups, Python,
Clippy, formatting and graph validation, excluding runtime tests and doctests.
Both stay within the no-more-than-10%-slower cold goal for this scope.

The cold comparisons use the ungrouped configurations before the fingerprint fix.
They supersede the earlier advisory-limit timings (main 42m34s → 38m38s;
FeynKit 45m57s → 36m54s), retained in JSON. All 17/26 archive inventories match;
58/91 artifact contexts have matching compiler-message multisets (738/827), with
zero compilation at the final archive stage. Clippy is counted separately.

Earlier independent FeynKit source comparisons, also before the fingerprint fix:

| Independent source edit | Baseline preparation | Candidate preparation | Change |
|---|---:|---:|---:|
| GammaLoop | 26m07s | 21m38s | −17.14% |
| CFF | 23m11s | 27m29s | +18.48% |

Each measures all eight archive groups, Python and Clippy in its own seeded store;
runtime and doctests are excluded. All 26 inventories and compiler-message
multisets match within each pair, and exact repeats start no builders. Host load
differed, so these mixed results do not establish a general source-change gain.
The fixed GammaLoop-edit candidate finishes preparation in 21m28s and repeats
with zero builders. It removes four compiler messages from three contexts,
including all repeated CFF compilation; the inherited `feynkit-ufo` compilation
in Python dependencies remains. Its five unaffected archive groups are reused.
The fixed CFF-edit candidate finishes in 24m02s with zero builders on repeat.
Both fixed scenarios match all 26 baseline archive inventories and 1,834 selected
tests. These are candidate-only validations, not new paired timing results.
Each edit starts independently from its clean application revision.

## Smaller cache closures and preserved boundaries

A complete integration-test input closure restored from a local binary cache in
15–18s instead of about 203s. With identical xz-6 encoding and eight-CPU affinity,
raw closure content fell from 19.782 to 1.758 GiB, encoded content from 3.874 to
1.149 GiB, and export time from 1,185.048s to 130.102s. The candidate no longer
retains eight older merged bundles through references. This measures local cache
encoding and import, not WAN throughput or NixCI billing.

Root-only encoding gave a different result: 557.62 MiB → 775.16 MiB, despite raw
NAR size falling 3.411 → 0.767 GiB. That sample excludes the baseline's referenced
older bundles. Both observations are retained; root size alone is insufficient.
Across the cold project outputs, distinct merged NAR sums fell approximately
75% (main 87.44 → 21.92 GiB; FeynKit 126.00 → 31.59 GiB).

Derivation inspection of the fixed FeynKit code preserves all nine extracted
production libraries after a GammaLoop edit and all eight siblings after a
Python-binding edit. CFF and model edits affect their dependency consumers;
embedded fixture changes affect test artifacts while production stays stable.
Manifest and lockfile edits still invalidate broadly. Seven earlier targeted
builds and exact repeats passed; those sequential probes test reuse boundaries,
not independent source timings. The fixed seven-edit derivation matrix preserves
the same boundaries.

The fingerprint fix corrects a repository bug: stripping `feynkit` previously
also removed `feynkit-cff` fingerprints. Matching the complete hash suffix stops
that deletion. A real GammaLoop consumer no longer compiles CFF (135.615s →
112.740s in one diagnostic pair). Full fixed-FeynKit validation removes five
compilation messages across four contexts, with no added compilation messages.
This does not turn the separate workspace-wide Clippy/doctest inputs into
per-crate checks. Python's inherited dependency-bundle invalidation also remains
broader than an ideal isolated leaf build.

## Correctness

Both earlier variants pass 1,633 main / 1,834 FeynKit selected tests, 43/62
doctests, Python behavior and pinned Clippy. They initially hit the same two UV
integration timeouts while compiler checks ran; identical-root recovery passed
once those compilers were cached. No test, timeout or filter was changed, and
failed timings are excluded from speed claims. Exact full repeats start no
builders or tests.

The fixed FeynKit revision `19a6460c658f411e7362060603c87d1d1e990af4` passes all
1,834 tests and 62 doctests on its first runtime attempt, with all 12 required
outputs registered and all 26 inventories unchanged. Its runtime phase takes
553.265s including 81 doctest compilation messages; an exact repeat takes 0.262s
and starts zero builders. This is unpaired correctness validation, not a cold
speed comparison. The supplied local license is used privately; no expired or
restricted-license diagnostic was observed. All 23 reporter regressions pass.

The Linux CLI builds and loads native libraries. Strict help/version exit-code
assertions still fail through unchanged Clap handling; no relaxed assertion or
application fix is included. Native Actions also exposes a preexisting Clippy
warning under Rust 1.98; pinned Nix Clippy passes.

## NixCI anomalies and excluded experiments

The [latest FeynKit baseline](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-benchmark-feynkit-baseline/2cf939655e8925a016770614219b7bf1b4ff0758)
passed all 12 required checks in 51m37s; final success followed 167s later.
An API test-dependency worker's log was replaced after a completed 1.4 GiB
restore. Both histories are preserved. Resource totals are lower bounds, and
this run is ineligible for paired performance acceptance. Its separate CLI
packaging artifact occupied a worker for 52m04s, including 40m38s of publication
intervals; it is absent from every primary check's dependency closure and is
excluded by the candidate. Necessary Python and test producers remain selected.

The [fixed FeynKit candidate](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency-feynkit/19a6460c658f411e7362060603c87d1d1e990af4)
is the first remote run with the fingerprint correction. Its build identities
change once, so earlier candidate runs do not establish a warm cache for it. Its
integration group failed one test before passing all 106 tests on the provider's
second attempt (125s then 99s of test execution; zero compilation in either
worker). The failure was a missing `feyn_gen_generation_test/model.json` while
saving initial test state. Eighteen existing tests use and clean that same
directory, consistent with a repository test-isolation race. These files are
unchanged between baseline and candidate. The failed worker and retry remain in
resource totals; no retry was dispatched by this task. A separate path-isolation
patch is prepared for review, with no assertions or test selection changed.

[Failed integration attempt](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency-feynkit/19a6460c658f411e7362060603c87d1d1e990af4/a222820f-9f79-492e-8abe-a478c273a4e7)
and [successful second attempt](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency-feynkit/19a6460c658f411e7362060603c87d1d1e990af4/69f787dc-0497-4809-8a91-082fd653e4bf)
retain the observed outcomes. Concurrent deletion is inferred from the code and
failure, not directly traced. The Python group subsequently passed and the whole suite finished in 79m29s.
All twelve latest service outcomes passed, but required-check timing remains
unknown: GitHub kept the first failed integration job URL while updating its
check to the successful retry's result and clocks. The reporter now preserves
explicit attempt numbers and marks this mismatch as incomplete evidence.
Offline replay preserves the original logs and main-pair metrics. This initial
FeynKit run is not an accepted speedup comparison.

The [GammaLoop dependency producer](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency-feynkit/19a6460c658f411e7362060603c87d1d1e990af4/0ede29a8-8d18-4460-87b7-265ca3dfc144)
occupied a worker for 46m36s, with 42m03s of completed publication intervals.
Of 159 publication operations, 64 Cargo bundles and merges account for 39m05s;
94 source/generated inputs account for 2m57s. The grouping uses logged names,
and these intervals do not measure pure network time or exclusive ownership.
The new fingerprint identities make this an initial-run observation; it does
not establish how much publication the independent source edits will avoid.

Repeated exact builds occur even when expected ordering edges exist and consumers
start after producers report successful uploads. In the unchanged main
[Linnet job](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency/7e6dc57327fbd6ab3e6181151d66963e8e93675c/ae96ed97-3d57-48cd-81eb-999dfde8b283),
Kurvst recompiles before the runner executable is built. It belongs to outer
`nix run` preparation; nested test-script settings cannot explain that case.
The repeated completed doctest result is a separate nested-build observation.
Cache visibility and substitution causes remain unknown; successful upload
messages alone do not establish availability to later workers. The fixed FeynKit
run contains a stronger example: the exact workspace dependency derivation
`h3fwjfsra58risirfsn5f08b589zi9is` and all 219 compilation messages recur in the
[tensor dependency worker](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency-feynkit/19a6460c658f411e7362060603c87d1d1e990af4/f1d07b96-965e-4897-a0d9-fafd25bf500e),
786s after the [producer published it](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency-feynkit/19a6460c658f411e7362060603c87d1d1e990af4/c6858f76-c99b-4c5d-b48b-4e31fd546c3e).
Its direct scheduling dependency is present. The cause of failed reuse remains
unknown. Ninety-three derivations repeat across the suite; many are cheap
source/metadata derivations, so that total must not be described as 93 repeated
compilations.

The earlier grouped-main producer took 1h57m07s, and required completion worsened
from 1h09m55s to 3h15m05s. Grouping was removed. Earlier suites also had low-speed
upload retries, HTTP/2 download errors and overwritten worker logs. The fifth
main baseline continued executing checks after a failed-suite snapshot, so that
snapshot was not terminal. Full details, exact identities and paired job links
remain in JSON. Queued jobs with observed ready prerequisites are waiting
observations, not proof of a provider stall. Upload intervals can overlap other
work and are not pure network time.

Both older Nix Actions runs passed:
[main, 2h10m52s](https://github.com/alphal00p/gammaloop/actions/runs/34227037456)
and [FeynKit, 2h04m32s](https://github.com/alphal00p/gammaloop/actions/runs/34227238762).
They include packaging and other outputs beyond the NixCI primary selection;
these are not equivalent-coverage speed comparisons. No Actions run was manually
dispatched during these experiments.

## Maintainer changes and remaining acceptance

At the recorded review, [PR 104](https://github.com/alphal00p/gammaloop/pull/104)
(`cf733769…`) enables synchronous discovery and retains manual ordering only for
impure runners and formatting. Our candidate already uses synchronous discovery,
with narrower job selection and additional generated explicit edges. Its graph
simplification is separate from the frozen experiment.
[PR 105](https://github.com/alphal00p/gammaloop/pull/105) (`13d6826d…`) simplifies
the local cache hook. Both reviewed hook versions use temporary metadata to
avoid unsigned local cache metadata shadowing signed server metadata. Saved
repeat-build logs contain no explicit signature diagnostic, so this mechanism
has not been established as their cause. Neither PR was merged by this task.

Complete the four source suites (the GammaLoop baseline is now running),
tracking required-check and individual-group latency, repeated compilation,
worker minutes and final-success delay. The restore and worker targets are
measured goals; correctness is mandatory. Report unmeasurable comparisons as
such, and stay within fourteen experimental suites. Draft suppression on NixCI
remains unsupported according to the maintainer; Actions draft/ready/manual and
cancellation changes remain available independently.


The eleventh suite, the [GammaLoop-source baseline](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-benchmark-feynkit-baseline-gammaloop/58909af3eaf0e42b8ecbb6c39526f1ebc42b6858),
is waiting after its [first artifact worker was abandoned](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-benchmark-feynkit-baseline-gammaloop/58909af3eaf0e42b8ecbb6c39526f1ebc42b6858/9d21561c-635c-4977-808d-15c7957fe2fb)
without a usable log. Its worker execution time remains unknown. A single-check
GitHub retry request at 07:54 UTC returned HTTP 404; no second request or new
push followed. Three experimental suites remain unsubmitted. The bounded log
collector stays active; no merge or performance acceptance follows from this
incomplete baseline. Four other branch heads had older completed suites and do
not establish a wider outage today.
