# CI measurements, 8 September 2026

Merged artifact output sizes fell about 75% in both layouts, with compilation
inputs and test inventories preserved. Remote worker savings remain unproven:
the first four suites had cache-upload timeouts and interrupted workers, which
prevent an accepted paired comparison. The final ungrouped comparison is now
running sequentially; five of twelve planned suites have been submitted.

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
has identity/graph validation, without a new cold or remote performance pair.

The 11 main / 12 FeynKit roots cover archive groups, Python packaging, Clippy,
formatting and graph validation. Local licensed execution and doctests are
excluded. All 17/26 archive inventories match names, filters and ignored flags.
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

The new [main baseline](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-benchmark-main-baseline/d44bf2ed781d73840175b349d2b7f15acf1a7589)
was submitted at 17:13:55–17:13:57 UTC on 8 September. It and the final main
candidate compare unchanged application code; the candidate CI configuration has
changed since its initial run. Actual cache warmth must be established from logs.
By 17:47:03 UTC, it had no application worker running, more than 31 minutes
after evaluation. Thirteen nodes appeared ready in 26 saved snapshots. Another
[phase-fix suite](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fraised_energy_cff_wip_optimized_phase_fix/32f27225460fbf5a2c9ad507355a8b8423fd237a)
had three running jobs at 17:41:24. This overlap may confound elapsed timings;
worker limits and the queue cause are unknown. Further benchmark pushes are held
for this baseline to finish. The remaining scenarios use independently prepared
FeynKit source edits, with no new Actions dispatch.

Both frozen Nix Actions workflows passed:
[main, 2h10m52s](https://github.com/alphal00p/gammaloop/actions/runs/34227037456)
and [FeynKit, 2h04m32s](https://github.com/alphal00p/gammaloop/actions/runs/34227238762).
They include packaging and other work beyond the primary NixCI selection, so
these are not equivalent-coverage speed comparisons. Their doctest jobs passed
43/62 cases and logged 73/81 compilation messages separately.

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
