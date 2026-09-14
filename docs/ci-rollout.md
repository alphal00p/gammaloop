# Consolidated CI changes

`codex/ci-consolidated` collects the validated main implementation from
`dd6ba768692302ca4855727110c4d8422b0979f9`, based on main
`395610143576507503fd2c785db3ba62340f4277`. The 58 development and measurement
commits are organized into five core review commits: mechanical flake extraction,
reporting, CI/artifact behavior, Actions triggers, and evidence/documentation.
Contributor instructions and the Justfile organization are follow-ups. The Nix
build definitions and CI selection are unchanged from the successful NixCI revision.

## What is included

- Rust build machinery in `nix/rust-workspace.nix`, with an explicit repository
  root; CI selection and groups in `nix/ci.nix`; generated `nix-ci.nix`.
- Automatic NixCI work limited to tests, Clippy, doctests, formatting, graph
  validation, necessary producers and final success. Packaging, docs, WASM and
  development outputs remain available through the flake.
- Per-crate artifacts, corrected sibling fingerprints, compatible compiler-state
  reuse, compressed handoffs, fewer redundant artifact merges, stable publication
  wrappers, shared check dependencies and Python preparation from test artifacts.
- Lightweight runners that look up the successful result before needing compiler
  state or test binaries. Test groups remain independently reported and scheduled.
- Synchronous dependency discovery, as discussed with Syd in PR #104. We retain
  our explicit selection and ordering; this is not a verbatim merge of his full
  dependency-graph simplification.
- Scoped local uploads adapted from Syd's PR #105. `just ci-checks-and-upload`
  checks first, then publishes selected outputs as the calling user. Entering a
  development shell does not enable a global upload hook.
- The existing Actions workflows, with draft guards, ready-for-review and manual
  triggers, and cancellation of superseded PR runs. Their branch filters and
  merge-queue coverage are preserved. NixCI itself still has no draft suppression.
- `just ci-report MANIFEST OUTPUT_DIR`, with JSON, CSV and readable comparisons,
  including incomplete/replaced log histories and queue observations.

The shared final-producer experiment and rejected artifact-sharing/storage modes
are excluded. The earlier standalone Clippy-cache trial is superseded by the
shared `cargoCheckArtifacts` implementation. Historical reports are retained;
their dated statements describe earlier stages, not the current rollout status.

## Final warm NixCI comparison

[Baseline](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-clan-cache-fixed-validation/f1ed70f585d293ecb851cf39f28821a69b33ac05)
versus [validated implementation](https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-lazy-check-runners/dd6ba768692302ca4855727110c4d8422b0979f9).

| Measurement | Before | After |
|---|---:|---:|
| Required checks | 8m26s | 3m46s |
| Whole suite | 8m37s | 4m05s |
| Observed worker minutes | 6.98 | 3.96 |
| Test-worker seconds, summed | 339.52 | 160.83 |
| Logged artifact restores | 1.31 GiB | 8,160 bytes |
| Final success after required checks | 11s | 19s |
| Build jobs cached | 53/53 | 53/53 |

All eight test jobs passed, with no compilation or executed test summaries.
Four show zero builds/one copied result, two show a copied 0.1 KiB output, and two
omit a final copy counter. Full copied paths are truncated; the generic reporter
conservatively labels their execution state unknown. Preflight verified all 61
selected outputs before the push.

This is one warm comparison. Artifact sizes are content proxies, not total
network traffic: one inner Nix command also displayed 13.4 MiB of download
progress during source lookup. Worker minutes are not billing or CHF. Different
waiting conditions contributed to the whole-suite improvement. No cold or
changed-source timing is inferred from this run.

Two observations remain for NixCI: 46 unchanged cache roots returned HTTP 404
before republishing, and the doctest job appeared started with empty logs before
a 78-second gap in observed suite worker activity ended. Their causes are
unknown. No abandoned workers, repeated builds, replaced logs or transfer errors
were observed. [Compact final measurements](ci-final-warm-comparison.json).

## Existing local evidence

The final artifact implementation passed 1,633 nextest tests and 43 doctests on
main, and 1,834 tests and 62 doctests in the temporary FeynKit integration, with
matching inventories and the existing skips. Its matched source-edit runtime
comparison improved main from 10m12s to 4m53s and FeynKit from 9m58s to 5m00s.
See [artifact preparation and correctness](ci-artifact-preparation.md).

Cold runtime testing still trails direct Cargo: 21m46s with Cargo versus 28m48s
through Nix for the same 1,633 selected main tests, including Python preparation.
That comparison excludes Clippy and doctests. See
[the controlled cold comparison](ci-cold-cargo-comparison.md).

The final runner change separately passed all eight cached checks with builds,
downloads, compiler artifacts and test binaries absent. A genuine cache miss ran
five Clinnet tests; a missing license still caused failure. Consolidation retains
those same check and artifact definitions; it does not change test coverage.

## Consolidation validation

All 63 compared derivation/output identities match the validated revision. The
11 intended checks (seven runtime groups, doctests, Clippy, formatting and graph
validation) pass in the preserved local store with builders and downloads
disabled. This reuses successful results; it does not rerun the test bodies.
All 24 reporter regression tests pass, and Just parses all CI commands. Runtime
files are byte-for-byte identical to the successful NixCI revision. The generated
configuration check passes, so no branch metadata regeneration is needed for main.

## Using and rolling out the change

Run `just ci-checks` for the intended local CI scope. Use
`just ci-checks-and-upload` to publish matching results before a push. Licensed
checks need `SYMBOLICA_LICENSE`; cache authentication comes from `NIXCI_NETRC` or
`~/.netrc`. Credentials stay outside Git. The already-activated Clan credential
configuration is separate from this repository; user upload access does not
configure the shared Nix daemon to trust/download from the cache.

After merging this branch into main, active branches should rebase onto main,
resolve their application-specific test groups and source filters, then run:

```sh
just ci-update
just ci-checks
```

For FeynKit, preserve its eighth group, nine extracted crates, embedded fixtures,
optional Python dependencies, boundary checks and stub checks when resolving the
rebase. Do not overwrite its group definition with main's seven groups. Regenerate
both the workspace graph and CI scheduling instead of retaining stale edges.
The shared `feynkit` branch was not changed by the experiments or consolidation.

The cache seed remains pinned to
`5181661ec340ebfb181a0045dac79fcec4f35525`. Keep its published history reachable;
consolidating commit history is not a reason to refresh this pin or invalidate
its cache. Refresh deliberately after a green compatible revision with
`just ci-cache-base REVISION`. Old published experiment references can be retired
once no active branch, cache pin or measurement depends on them.

[Branch provenance](ci-branch-provenance.json) identifies included, superseded and
layout-specific work. Original branches/worktrees and local archive refs are
preserved, including unfinished experiment patches. The main checkout and shared
FeynKit branch are untouched. Consolidation does not itself claim a new hosted
performance measurement; the linked run identifies the exact validated code.
