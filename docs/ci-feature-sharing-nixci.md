# CI measurements, 9 September 2026

## Pushed feature-sharing changes

The main candidate is [983fc214f](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency/983fc214fc27a78b08c79f31cccd14380910d3a7).
The FeynKit candidate is [1a2e7a199](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency-feynkit/1a2e7a1998d87e5dcdddaf95c54875293379bf5c).
Main also includes the intervening artifact-fingerprint correction, so its timing
cannot isolate feature sharing. FeynKit compares the same application revision.

FeynKit passed all 12 required checks. Required-check latency was 9,545s
(2h39m05s), whole-suite latency 9,558s, and final-success tail 13s. About 70 minutes
elapsed before build workers started. Observed worker time was 173.68 minutes;
reported restored content was 29,693,769,912 bytes, including 26,256,069,637 bytes
for intermediate producers. The logs contain 376 Cargo `Compiling` progress
messages. These are content and activity proxies, not network bytes or CHF costs.
One abandoned worker had no log, so worker accounting remains a lower bound.
This sample does not demonstrate an additional feature-sharing speedup.

Main is still running at the time of this interim note. Its completed report will
be collected separately; the result above is not a main/FeynKit crate-split
comparison.

## Targeted follow-up

The final reserved NixCI sample is
[07d0efc91](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency-feynkit/07d0efc910cf4c3175142c808e476fe7e393f815),
pushed at 13:20:05 UTC. It removes build-time Python extension inputs and stops
scheduling the production GammaLoop artifact separately from the Python module.
The local tests, dependency proof and rationale are in
[ci-python-runtime-inputs.md](ci-python-runtime-inputs.md). This uses suite 14,
the campaign's last reserved suite; no further experiments were pushed.

Both FeynKit revisions passed all 12 required checks, with matching reported
runtime counts in all nine executing test/doctest groups. The changed archives
also have matching local inventories. The application code is identical.

| Observed metric | Feature sharing, 1a2e7a199 | Python dependency fix, 07d0efc91 |
|---|---:|---:|
| Push start to required checks | 2h40m14.5s | 32m10.5s |
| NixCI suite start to required checks | 2h39m05s | 24m46s |
| Whole suite, from NixCI suite start | 2h39m18s | 24m51s |
| Push start to first selected worker | 70m28.5s | 12m06.3s |
| Observed worker minutes, lower bounds | 173.68 | 54.45 |
| Reported intermediate restored content | 24.45 GiB | 23.49 GiB |
| Publication intervals summed across workers | 90m44.9s | 8m02.8s |
| Cargo `Compiling` progress messages | 376 | 148 |
| Final success after last required check | 13s | 5s |

This is a promising observed improvement, with two distinct contributors: less
waiting before workers started, and less preparation after they started. The
[Python module job](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency-feynkit/07d0efc910cf4c3175142c808e476fe7e393f815/55746db1-b23e-4a5d-9fa3-4cf0faef3b9f)
reused its final output in 17.83s with no Cargo compilation, without scheduling
the preceding 55m44s production job. That dependency behavior matches the local
probe. The final critical check was doctests, at 12m01s of worker execution.

The full worker reduction cannot be attributed to the code change alone: remote
cache state and queues differ, one baseline worker was abandoned without a log,
and a candidate worker's log was replaced. The reporter marks both samples as
incomplete evidence and leaves formal percentage comparisons unset. All visible
jobs finished successfully; no missing work or network bytes are estimated.
Restored intermediate content remains large, and zero-Cargo warm runs have not
been achieved. Our per-commit scheduling wrappers intentionally realize stable
artifacts again to work around previously memoized producer results whose outputs
were no longer available to consumers. The symlinks themselves are tiny, but
realizing their targets can still restore large closures. Removing that remaining
preparation safely needs reliable cache availability between workers, or scheduling
that skips intermediate preparation when the final test archive is available. This pair does not establish the cold-build or changed-source
acceptance targets.

The combined JSON, CSV and readable comparison is in
`/tmp/gammaloop-ci-validation/python-runtime-only/remote/reports/watched-suites-2026-09-09T13-53-17.122Z/`.
`paired-runtime-counts.json` records the matching execution summaries separately
from timing, and `first-worker-timings.json` records the initial waiting boundary.

## Provider observations

A controlled local file-cache round trip of the exact repeatedly rebuilt
workspace-hack artifact (`6y3l3h9ipzzmmj345mz74ynhsqq3g0z2`) succeeded. Its
191,959,144-byte NAR restored into a fresh store in 3.87s with `--max-jobs 0`,
no builders and an identical NAR hash. Export to the local Zstd file cache took
6.54s. This used local Nix 2.34.8, no network or uploads to NixCI, and a controlled
unsigned file cache. It does not validate NixCI authentication, signing, retention
or cross-worker visibility. Evidence is in
`/tmp/gammaloop-ci-validation/cache-roundtrip-main/`.

The Python archive already uses nextest's binary filter. In the locally validated
FeynKit candidate it is 8,909,815 compressed file bytes; its Nix closure is one
8,910,152-byte store object. The ordinary integration archive is 290,687,419 file
bytes, also with no referenced store paths. A larger hidden closure does not
explain publication of these particular final archives. Raw evidence is in
`/tmp/gammaloop-ci-validation/python-runtime-only/feynkit-archive-closures.json`.

Main's [Python archive worker](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency/983fc214fc27a78b08c79f31cccd14380910d3a7/8e8b3bdc-24f7-4cbd-8813-115c2c6051c7)
was still rebuilding and publishing vendored `cargo-package-*` source outputs at
14:08 UTC, before its archive build. Its publication intervals therefore must
not be described as time spent uploading the final Python archive.


- Main's [dependency-discovery job](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency/983fc214fc27a78b08c79f31cccd14380910d3a7/9fcea729-c325-4080-bee7-9c9d2f8969e6)
  reported six cache HTTP 502 retries, was observed abandoned, and later exposed a
  replacement worker log at the same URL. Both streams were retained; the retry
  trigger is unknown.
- FeynKit's [kinematics test-binary job](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency-feynkit/1a2e7a1998d87e5dcdddaf95c54875293379bf5c/249f68a1-e877-43d1-a73a-3067adc48343)
  was observed abandoned with an empty log. Its missing work cannot be inferred.
- Main [uploaded its workspace-hack dependency result](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency/983fc214fc27a78b08c79f31cccd14380910d3a7/3f83b350-d788-4bbf-8e85-5e33f2daaa58)
  by 11:55:56 UTC, but the [GammaLoop producer](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency/983fc214fc27a78b08c79f31cccd14380910d3a7/50a9ea92-0711-4044-9a3f-3f26c0b3b93e)
  rebuilt the identical `6y3l3h9ipzzmmj345mz74ynhsqq3g0z2` derivation at 12:14:17.
  It also repeated the workspace dependency derivation already
  [built and uploaded by cargoArtifacts](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency/983fc214fc27a78b08c79f31cccd14380910d3a7/b73fcfea-3bf2-4acb-a0d5-bb3fe7f1c844).
  These are exact derivation matches, not just identical crate names. The five
  inspected repeated dependency derivations do not disable substitution, force
  local builds or declare impure builds. Nix
  [documents `allowSubstitutes = false`](https://nix.dev/manual/nix/2.34/language/advanced-attributes.html#adv-attr-allowSubstitutes)
  as the switch that would intentionally skip cache substitution; it is absent
  here. Cache access and publication still need provider-side investigation.
- The [previous FeynKit Spenso build](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency-feynkit/19a6460c658f411e7362060603c87d1d1e990af4/524e55ef-cad3-4652-85d0-da7610afbbbf)
  uploaded `nzrd4y6a7v26930ai5r5q30h9n0dh0w3`, but the
  [new Spenso job](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-efficiency-feynkit/1a2e7a1998d87e5dcdddaf95c54875293379bf5c/b7f76a76-2eed-482f-9fc6-aa38c169d4ca)
  rebuilt it. Cache retention, visibility and substitution behavior need
  investigation; these logs do not establish the cause.

The completed/interim JSON, CSV and readable reports and retained worker logs are
under `/tmp/gammaloop-ci-validation/test-feature-unification/remote/`. The paired
follow-up manifest and reports are under
`/tmp/gammaloop-ci-validation/python-runtime-only/remote/`.
