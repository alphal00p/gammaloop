= CI maintenance and measurements

Use `just ci-checks` for the selected local suite and `just ci-checks-and-upload`
before pushing. The latter checks first, then publishes matching outputs to
NixCI. See #link("../CONTRIBUTING.typ#nixci-cache")[CONTRIBUTING.typ] for credentials,
licenses and retries. `just check` runs Cargo checking only.

== Configuration and cache boundaries

#table(
  columns: (1fr, 1fr),
  table.header([File], [Responsibility]),
  [`flake.nix`], [Compose public packages, checks, apps and development shells],
  [`nix/rust-workspace.nix`], [Rust sources, features, build profiles and artifact reuse],
  [`nix/ci.nix`], [Test groups, required producers and NixCI scheduling],
  [`nix/ci-workspace-graph.json`, `nix-ci.nix`], [Generated workspace graph and self-contained NixCI configuration],
  [`just/ci.just`], [Local checks, uploads, regeneration and reporting],
  [`.github/scripts/ci-report.mjs`], [Collect and compare existing CI runs],
)

Automatic NixCI work covers tests, Clippy, doctests, formatting, graph validation,
necessary producers and final success. Packaging, documentation and WASM remain
available through the flake; `nix flake check --impure` includes those extra checks.

Compilation stays per crate. Source filtering preserves unaffected packages;
manifest and lockfile changes still invalidate broadly. Test dependencies avoid
multi-crate cycles: Linnet enables its own optional test features, and Spenso owns
the macro integration test. Regenerate Hakari after dependency/feature changes so
the shared cache does not retain unused build dependencies. Compatible test and
Python builds share dependencies, while Python retains its required ABI, interpreter
and features. Clippy and doctests share check dependencies but retain workspace-wide
source inputs, so their rebuild isolation differs from package test artifacts.

Merged artifacts are self-contained compressed archives. Recursive inheritance,
writable extraction, Cargo fingerprints and epoch-1 timestamps must survive
changes to archive handling. Stable publication outputs depend on the artifacts,
not the commit. Lightweight test runners look up the successful check result
before fetching compiler state or test binaries. Necessary artifact producers
remain explicitly selected for publication; an empty successful test output alone
does not retain them. Test groups are scheduled and reported independently.

Synchronous dependency discovery incorporates Syd's PR \#104 suggestion while
retaining the explicit graph. Local uploads adapt \#105 through the Just command;
entering a shell installs no global upload hook. The existing Actions workflows
skip draft PRs, react to ready-for-review/manual events and cancel superseded PR
runs, retaining branch filters and merge-queue coverage. NixCI draft suppression
remains a separate integration issue.

After rebasing an active branch onto the CI changes, preserve its test groups and
source filters, run `just ci-update`, and commit both generated files. FeynKit
needs its eighth group, nine extracted crates, embedded fixtures, optional Python
dependencies, boundary checks and Python stub checks. Its extracted libraries
must remain independent of GammaLoop and its private workspace-hack dependency.
The shared FeynKit branch was not changed by the experiments.

The `ci-cache-base` input pins a compatible compiler-state seed. Keep that commit
reachable in published history. Refresh deliberately after a green revision with
`just ci-cache-base FULL_COMMIT_SHA`; rebasing or cleaning history alone is not a
reason to invalidate it. See the #link("architecture/nix-crane-cache-reuse.typ")[cache-reuse audit]
for earlier investigation and implementation details.

== Measured results

The Cargo dependency cleanup was compared locally with `246c03693` on Rust
1.98.1, using eight jobs, `ci-optim`, no incremental compilation, and the same 631
selected tests from the Linnet and Spenso groups (658 listed, identical filters).

#table(
  columns: (1fr, 1fr, 1fr),
  table.header([Local Cargo workflow], [Before], [After]),
  [Empty-target test build and execution], [6m52s], [6m13s],
  [Empty-target check followed by tests], [10m31s], [7m12s],
)

These are single pairs on a shared host with vendored dependencies available,
not full-suite or remote measurements. Direct test CPU time fell 4.4%; most of
the larger combined gain comes from checking. Both unchanged retries compiled
nothing. Nix archive inventories matched all 631 tests; the selected Nix suite
passed 1,633 runtime tests (including five Python cases), 43 doctests and static
checks. FeynKit graph evaluation preserved its groups and extracted-crate cache
boundaries; its runtime suite was not rerun for this cleanup.

The final warm comparison used the same main application layout:
#link("https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-clan-cache-fixed-validation/f1ed70f585d293ecb851cf39f28821a69b33ac05")[baseline]
and #link("https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-lazy-check-runners/dd6ba768692302ca4855727110c4d8422b0979f9")[candidate].

#table(
  columns: (1fr, 1fr, 1fr),
  table.header([Warm NixCI measurement], [Before], [After]),
  [Required checks], [8m26s], [3m46s],
  [Whole suite], [8m37s], [4m05s],
  [Observed worker minutes], [6.98], [3.96],
  [Logged artifact restores], [1.31 GiB], [8,160 bytes],
  [Final success after required checks], [11s], [19s],
  [Cached build jobs], [53/53], [53/53],
)

All eight test jobs passed without compilation or executed test summaries. Their
truncated output paths prevent the generic reporter from proving result reuse
for every job. This is one warm pair; waiting conditions also differed. Reported
artifact sizes exclude other traffic, including 13.4 MiB shown by one inner Nix
command. Worker minutes are observed log spans, not billed time or CHF.

Two unresolved observations accompany this pair: 46 previously readable unchanged
cache roots returned HTTP 404 before republishing, and the doctest appeared started
with empty logs before a 78-second gap in observed worker activity ended. Their
causes are unknown. No abandoned workers, repeated builds, replaced logs or
transfer errors were observed in the completed pair.

Earlier local comparisons measured different scopes:

#table(
  columns: (1fr, 1fr, 1fr),
  table.header([Controlled local runtime comparison], [Before / Cargo], [After / Nix]),
  [Main source edit, including evaluation/setup], [10m12s], [4m53s],
  [FeynKit source edit, including evaluation/setup], [9m58s], [5m00s],
  [Main project-cold build and tests: Cargo versus Nix], [21m46s], [28m48s],
)

These include Python preparation and all runtime groups, but exclude Clippy and
doctests. Each comparison matches application code within its own layout; it does
not measure the crate split itself. The cold pair used identical toolchains and
eight CPUs, an empty Cargo target and a disposable Nix store without project
outputs. External toolchains, native libraries and crate sources were available;
remote builders, downloads and uploads were disabled. Nix still took 32.3% longer.

The complete local correctness runs passed 1,633 tests/43 doctests on main and
1,834 tests/62 doctests on FeynKit, with matching inventories, filters and existing
skips, including five Python cases per layout. FeynKit's boundary checks passed;
its Python-stub workflow was preserved but stub generation was not rerun. A
pre-existing probe expecting a `feynkit-py` edit to affect only that crate failed:
its API consumer also depends on it. That assertion was not changed. An unchanged
retry passed two intermittent main state-file failures; that failed attempt is
excluded from successful timings. The final runner fix separately passed eight
cached checks without compiler artifacts/downloads, a real five-test Clinnet cache
miss, and a missing-license failure. Consolidation preserved all 63 compared Nix
output identities and all 11 required local checks; those checks reused results.

Detailed, immutable evidence remains in the pre-cleanup commit:
#link("https://github.com/alphal00p/gammaloop/blob/0cc8cd2bb42a2e5d1bebe468df347d1d5420d897/docs/ci-final-warm-comparison.json")[final warm data], #link("https://github.com/alphal00p/gammaloop/blob/0cc8cd2bb42a2e5d1bebe468df347d1d5420d897/docs/ci-artifact-preparation-results.json")[local artifact and coverage results],
#link("https://github.com/alphal00p/gammaloop/blob/0cc8cd2bb42a2e5d1bebe468df347d1d5420d897/docs/ci-cold-cargo-comparison.json")[cold commands and results], and #link("https://github.com/alphal00p/gammaloop/tree/0cc8cd2bb42a2e5d1bebe468df347d1d5420d897/docs")[experiment history]. This
keeps generated reports and experiment diaries out of the maintained source tree.

== Collecting another comparison

Run `just ci-report MANIFEST OUTPUT_DIR`. It selects pinned Node and uses an
authenticated `gh` for GitHub. The reporter only reads existing runs; it never
launches, retries or cancels jobs. Put generated output outside the source tree.
A minimal manifest for reading an existing suite is:

```json
{
  "repository": "alphal00p/gammaloop",
  "runs": [{
    "layout": "main",
    "variant": "candidate",
    "scenario": "warm",
    "sha": "dd6ba768692302ca4855727110c4d8422b0979f9",
    "suiteUrl": "https://nix-ci.com/gh:alphal00p:gammaloop/codex%2Fci-lazy-check-runners/dd6ba768692302ca4855727110c4d8422b0979f9"
  }]
}
```

Add a `baseline` entry with matching `layout`, `scenario` and optional `pair`
(default `"1"`) for paired percentages. Use a new pair label for repeats. Each
entry requires a full SHA and exact commit-level suite URL. Duplicate suites and
ambiguous pairs are rejected. Optional `actionRunIds` limits workflows; otherwise
all Actions runs for the SHA and every attempt are retained.

For controlled comparisons, supply the same explicit `requiredAttributes` in both
entries: every test group, doctest, Clippy, formatting and graph check. Otherwise
the default is all observed test jobs plus the three static checks. A runner's
build alone does not satisfy its required test phase. Missing groups, incomplete
logs, failed suites or ambiguous clocks prevent accepted paired percentages.
Coverage equivalence is never inferred: compare inventories, filters, features,
Python behavior, license settings and static-check results separately.

The output contains `report.md`, `report.json`, CSV tables for suites, jobs,
Actions, transfers, compilations, repeated derivations, incidents and pairs, plus
`raw/N/` evidence and the input manifest. Use a fresh directory for each snapshot;
existing files are overwritten. New raw files are owner-readable/writable; review
them before sharing. Exit status is 0 for complete evidence, 2 for a partial report,
and 1 for invalid invocation or manifest.

Metric interpretation:

- Required latency ends at the last intended check; whole-suite latency ends at
  the last completed job. Final-success tail is reported separately. Recorded
  `submittedAt`/`submissionCompletedAt` push bounds can additionally measure elapsed
  time from submission; commit author dates are not submission clocks.
- Worker minutes sum observed worker spans. Check duration and pre/post-worker
  gaps stay separate; a pre-worker gap does not prove queue time.
- Only completed transfer messages supply sizes/times. Intervals are unioned per
  worker before summing, because downloads and uploads can overlap. Content sizes
  are not wire bytes; timers may include decompression, import and locking.
- `Compiling`, `Checking` and Cargo command durations remain separate by context.
  Repeated crate names can mean different features or profiles. Only the same
  full derivation path built under different job URLs establishes repeated Nix work.
- Test execution is `executed`, `reused` or `unknown`. Fetching a runner wrapper
  does not prove reuse of its test result. Nextest `failed` includes `timedOut`;
  do not add those counters together.

To replay without network calls, add `offline` to an entry, with `suite`, `checks`,
`jobs` and `actions` paths pointing at a saved `raw/N/` snapshot's corresponding
JSON files. Paths resolve relative to the manifest; log paths in `jobs.json`
resolve relative to that index. The reporter accepts the earlier collector's job
index too. Keep raw evidence alongside the report so replays remain possible.

To preserve unusual service behavior, retain snapshots when observed; final logs
may overwrite earlier attempts. The optional manifest fields are:

#table(
  columns: (1fr, 1fr),
  table.header([Field], [Evidence to save]),
  [`observationFiles`], [Snapshot JSON with `at` and `suites: [{spec: {sha, suiteUrl}, suite: API_RESPONSE}]`],
  [`dependencySnapshot`], [`{sha, file}` for that exact revision's generated config containing `dependencies`],
  [`observedInterruptions`], [Entries `{jobUrl, observedAt, evidenceFile}`; evidence repeats the URL with `status: "abandoned"` and optional prior log `file`],
  [`observedLogReplacements`], [Same entry fields; evidence repeats URL/time with earlier `file` and later `replacement.file`],
)

Evidence paths resolve relative to the manifest, linked logs relative to their
evidence file. The reporter copies them into `raw/N/` for offline replay; point
future manifests at those copies. Replaced logs and abandoned work make current
resource totals lower bounds, even after a successful retry. Prior totals are not
added because streams may overlap. Queue observations require every declared
prerequisite to be successful in the same snapshot; they do not prove continuous
queue time or its cause. Transfer errors, invalid clocks, mismatched outcomes and
repeated attempts retain evidence without attributing a backend cause.

Run the reporter's offline regression checks with
`node --test .github/scripts/ci-report.test.mjs`.
