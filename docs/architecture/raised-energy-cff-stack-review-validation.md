# Reviewed CFF stack: current validation

All seven boundaries and all required final runtime selections have complete passing collector receipts. C1/C2 reuse exact unchanged revision certificates; C3–C7 and final runtime are freshly executed after the accepted GL00/GL04 certificate amendment. Document rendering and final history closure are separate evidence, described below.

Boundary snapshot: `2026-09-11T19:15:44.234810+00:00`. Completed passing boundaries: 1, 2, 3, 4, 5, 6, 7.

Boundary executions, repeated feature configurations and the final runtime selections remain separate counts. A skipped identity reported by nextest is outside that selected run; it is not a passing execution.

Raw commands, lists, logs, guard records and source audits belong in the [fresh closeout evidence archive](raised-energy-cff-closeout-evidence.tar.gz) and its [JSON receipt](raised-energy-cff-closeout-evidence.json). Independent verification confirms all 2,484 archive members, their hashes and metadata, and their recorded source links. Final document and history closure have separate receipts. Paths below are relative to the closeout capture directory and must be resolved through their manifest. No ignored target-directory log is used as a report link.

## Execution identities

| Boundary | Source revision | Tree | Status |
| --- | --- | --- | --- |
| 1 — Symbolic and tensor foundations | `97eb06cbd3b09cd62c242f27da668accf68173c8` | `5c270b8f8e491416413667f00e257a7e4429b62d` | Passed — unchanged revision certificate reused |
| 2 — Exact generalized CFF | `2ba79605128090281bf4ebffdfa0f1c52fb4ba7f` | `e28fd22e70c1737f02f8951ad8c3872afb623fa1` | Passed — unchanged revision certificate reused |
| 3 — Raised-energy CFF and local UV | `1f0d4fd6141b1e572fe780824dc2579142f2343f` | `d4171b93510ce534115058b86a4db29a3ee73763` | Passed — fresh execution |
| 4 — Command and evaluation workflows | `48104fa453f8a8a35077ec8a5e1d9ad94e7a9c16` | `c2e8c3349caccfb6e1479a2b6c5f48baa1b3a32f` | Passed — fresh execution |
| 5 — Physical phases and model sewing | `32a57c5dd5beedcdf7e42c17091baf7a6e9a6fec` | `b91576e77e880b1aa277a409a9be81420f0a1d7b` | Passed — fresh execution |
| 6 — D-dimensional integrated UV | `ee798ce279f280000537e4d8c05ecc83fc6ec2cf` | `98fa4f223071d411bdd942eaaf66cb6addf5f011` | Passed — fresh execution |
| 7 — Final stack and documentation | `9e260684581ab2893eff4032c5ca8ea1fc2c18cb` | `9cebc793fe813699b9d689f1ab245715902ac04e` | Passed — fresh execution |

Final runtime execution pin: `9e260684581ab2893eff4032c5ca8ea1fc2c18cb`; tree `9cebc793fe813699b9d689f1ab245715902ac04e`. Snapshot: `2026-09-11T19:15:44.870994+00:00`.

Completed boundary receipts require pre/post tracked-source audits, matching revision/tree records, unchanged contents during source-mtime refresh, and final head equality. C1 refreshes all tracked regular-file mtimes because prior cache provenance is unknown; subsequent boundaries refresh changed blobs/modes against the preceding audited checkout. The first amended C3 checkout records its explicit pre-v2 source reference; it does not infer that reference from the retained runtime identifier. This guards against stale build reuse across checkouts. It does not infer identity from a requested revision alone.

## Per-commit gates and feature checks

Each completed boundary has formatting before locked workspace checking, a locked all-target `dev-optim` build, and Clippy. Guard time includes the complete command; it is distinct from nextest execution time. Peak memory is process-tree RSS in decimal GB.

| Boundary | Step | Result | Guard time | Peak RSS | Limit | Command |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | fmt | Pass | 3.796 s | 0.121 GB | 30 GB | [recorded](#command-c1-fmt) |
| 1 | check | Pass | 17.450 s | 2.560 GB | 30 GB | [recorded](#command-c1-check) |
| 1 | build | Pass | 288.729 s | 7.496 GB | 30 GB | [recorded](#command-c1-build) |
| 1 | clippy | Pass | 28.756 s | 3.394 GB | 30 GB | [recorded](#command-c1-clippy) |
| 2 | fmt | Pass | 3.757 s | 0.121 GB | 30 GB | [recorded](#command-c2-fmt) |
| 2 | check | Pass | 1.941 s | 0.516 GB | 30 GB | [recorded](#command-c2-check) |
| 2 | build | Pass | 9.006 s | 0.890 GB | 30 GB | [recorded](#command-c2-build) |
| 2 | clippy | Pass | 26.547 s | 2.805 GB | 30 GB | [recorded](#command-c2-clippy) |
| 2 | shared-all | Pass | 2.905 s | 0.516 GB | 30 GB | [recorded](#command-c2-shared-all) |
| 2 | shared-no-default | Pass | 1.517 s | 0.235 GB | 30 GB | [recorded](#command-c2-shared-no-default) |
| 3 | fmt | Pass | 4.609 s | 0.140 GB | 30 GB | [recorded](#command-c3-fmt) |
| 3 | check | Pass | 13.867 s | 2.338 GB | 30 GB | [recorded](#command-c3-check) |
| 3 | build | Pass | 454.618 s | 7.991 GB | 30 GB | [recorded](#command-c3-build) |
| 3 | clippy | Pass | 24.170 s | 2.591 GB | 30 GB | [recorded](#command-c3-clippy) |
| 3 | shared-all | Pass | 0.621 s | 0.074 GB | 30 GB | [recorded](#command-c3-shared-all) |
| 3 | shared-no-default | Pass | 0.623 s | 0.074 GB | 30 GB | [recorded](#command-c3-shared-no-default) |
| 4 | fmt | Pass | 4.517 s | 0.151 GB | 30 GB | [recorded](#command-c4-fmt) |
| 4 | check | Pass | 28.088 s | 3.871 GB | 30 GB | [recorded](#command-c4-check) |
| 4 | build | Pass | 1198.346 s | 12.766 GB | 30 GB | [recorded](#command-c4-build) |
| 4 | clippy | Pass | 38.850 s | 4.349 GB | 30 GB | [recorded](#command-c4-clippy) |
| 4 | shared-all | Pass | 0.610 s | 0.074 GB | 30 GB | [recorded](#command-c4-shared-all) |
| 4 | shared-no-default | Pass | 0.612 s | 0.074 GB | 30 GB | [recorded](#command-c4-shared-no-default) |
| 4 | python-api | Pass | 11.333 s | 1.321 GB | 30 GB | [recorded](#command-c4-python-api) |
| 5 | fmt | Pass | 4.446 s | 0.143 GB | 30 GB | [recorded](#command-c5-fmt) |
| 5 | check | Pass | 25.002 s | 3.397 GB | 30 GB | [recorded](#command-c5-check) |
| 5 | build | Pass | 1131.240 s | 12.361 GB | 30 GB | [recorded](#command-c5-build) |
| 5 | clippy | Pass | 35.865 s | 3.988 GB | 30 GB | [recorded](#command-c5-clippy) |
| 5 | shared-all | Pass | 0.599 s | 0.074 GB | 30 GB | [recorded](#command-c5-shared-all) |
| 5 | shared-no-default | Pass | 0.600 s | 0.074 GB | 30 GB | [recorded](#command-c5-shared-no-default) |
| 6 | fmt | Pass | 4.851 s | 0.144 GB | 30 GB | [recorded](#command-c6-fmt) |
| 6 | check | Pass | 27.887 s | 3.686 GB | 30 GB | [recorded](#command-c6-check) |
| 6 | build | Pass | 1304.772 s | 11.517 GB | 30 GB | [recorded](#command-c6-build) |
| 6 | clippy | Pass | 40.040 s | 4.096 GB | 30 GB | [recorded](#command-c6-clippy) |
| 6 | shared-all | Pass | 0.621 s | 0.074 GB | 30 GB | [recorded](#command-c6-shared-all) |
| 6 | shared-no-default | Pass | 0.616 s | 0.074 GB | 30 GB | [recorded](#command-c6-shared-no-default) |
| 6 | vakint-community | Pass | 12.323 s | 2.417 GB | 30 GB | [recorded](#command-c6-vakint-community) |
| 7 | fmt | Pass | 4.479 s | 0.143 GB | 30 GB | [recorded](#command-c7-fmt) |
| 7 | check | Pass | 0.614 s | 0.074 GB | 30 GB | [recorded](#command-c7-check) |
| 7 | build | Pass | 0.616 s | 0.074 GB | 30 GB | [recorded](#command-c7-build) |
| 7 | clippy | Pass | 0.617 s | 0.095 GB | 30 GB | [recorded](#command-c7-clippy) |
| 7 | shared-all | Pass | 0.614 s | 0.074 GB | 30 GB | [recorded](#command-c7-shared-all) |
| 7 | shared-no-default | Pass | 0.612 s | 0.074 GB | 30 GB | [recorded](#command-c7-shared-no-default) |

Shared CFF is absent at C1. From C2 onward, all-feature and no-default-feature all-target checks are required. The Python API check belongs to C4 and the Vakint community-module check to C6; all of these required checks have passed.

## Per-commit behavioral selections

| Boundary | Selection | Passed/executed | Nonpass | Skipped | Test time | List/run match | Command |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | commit1-foundations | 693/693 | 0 | 28 | 9.532 s | Verified | [recorded](#command-c1-commit1-foundations-run) |
| 2 | shared-all | 108/108 | 0 | 0 | 21.176 s | Verified | [recorded](#command-c2-shared-all-run) |
| 2 | shared-default | 37/37 | 0 | 0 | 5.310 s | Verified | [recorded](#command-c2-shared-default-run) |
| 2 | shared-no-default | 37/37 | 0 | 0 | 5.167 s | Verified | [recorded](#command-c2-shared-no-default-run) |
| 3 | commit3-cff-uv | 219/219 | 0 | 486 | 140.906 s | Verified | [recorded](#command-c3-commit3-cff-uv-run) |
| 3 | shared-all | 108/108 | 0 | 0 | 21.871 s | Verified | [recorded](#command-c3-shared-all-run) |
| 3 | shared-default | 37/37 | 0 | 0 | 5.256 s | Verified | [recorded](#command-c3-shared-default-run) |
| 3 | shared-no-default | 37/37 | 0 | 0 | 5.182 s | Verified | [recorded](#command-c3-shared-no-default-run) |
| 4 | commit4-scalar-smoke | 19/19 | 0 | 242 | 133.831 s | Verified | [recorded](#command-c4-commit4-scalar-smoke-run) |
| 4 | commit4-workflows | 599/599 | 0 | 829 | 187.727 s | Verified | [recorded](#command-c4-commit4-workflows-run) |
| 4 | shared-all | 108/108 | 0 | 0 | 21.444 s | Verified | [recorded](#command-c4-shared-all-run) |
| 4 | shared-default | 37/37 | 0 | 0 | 5.354 s | Verified | [recorded](#command-c4-shared-default-run) |
| 4 | shared-no-default | 37/37 | 0 | 0 | 5.158 s | Verified | [recorded](#command-c4-shared-no-default-run) |
| 4 | signed-acceptances | 4/4 | 0 | 0 | 120.172 s | Verified | [recorded](#command-c4-signed-acceptances-run) |
| 4 | source-certificates | 1/1 | 0 | 691 | 0.125 s | Verified | [recorded](#command-c4-source-certificates-run) |
| 5 | commit5-phases | 211/211 | 0 | 1199 | 31.746 s | Verified | [recorded](#command-c5-commit5-phases-run) |
| 5 | shared-all | 108/108 | 0 | 0 | 21.465 s | Verified | [recorded](#command-c5-shared-all-run) |
| 5 | shared-default | 37/37 | 0 | 0 | 5.224 s | Verified | [recorded](#command-c5-shared-default-run) |
| 5 | shared-no-default | 37/37 | 0 | 0 | 5.173 s | Verified | [recorded](#command-c5-shared-no-default-run) |
| 5 | signed-acceptances | 4/4 | 0 | 0 | 119.769 s | Verified | [recorded](#command-c5-signed-acceptances-run) |
| 5 | source-certificates | 1/1 | 0 | 718 | 0.087 s | Verified | [recorded](#command-c5-source-certificates-run) |
| 5 | vertex-rules | 1/1 | 0 | 718 | 0.302 s | Verified | [recorded](#command-c5-vertex-rules-run) |
| 6 | commit6-vakint-uv | 241/241 | 0 | 619 | 73.850 s | Verified | [recorded](#command-c6-commit6-vakint-uv-run) |
| 6 | shared-all | 108/108 | 0 | 0 | 21.475 s | Verified | [recorded](#command-c6-shared-all-run) |
| 6 | shared-default | 37/37 | 0 | 0 | 5.217 s | Verified | [recorded](#command-c6-shared-default-run) |
| 6 | shared-no-default | 37/37 | 0 | 0 | 5.113 s | Verified | [recorded](#command-c6-shared-no-default-run) |
| 6 | source-certificates | 1/1 | 0 | 725 | 0.089 s | Verified | [recorded](#command-c6-source-certificates-run) |
| 7 | commit5-phases | 212/212 | 0 | 1205 | 30.817 s | Verified | [recorded](#command-c7-commit5-phases-run) |
| 7 | commit6-vakint-uv | 241/241 | 0 | 619 | 73.977 s | Verified | [recorded](#command-c7-commit6-vakint-uv-run) |
| 7 | curated | 2065/2065 | 0 | 297 | 465.870 s | Verified | [recorded](#command-c7-curated-run) |
| 7 | scalar-all | 166/166 | 0 | 107 | 922.523 s | Verified | [recorded](#command-c7-scalar-all-run) |
| 7 | shared-all | 108/108 | 0 | 0 | 21.324 s | Verified | [recorded](#command-c7-shared-all-run) |
| 7 | shared-default | 37/37 | 0 | 0 | 5.191 s | Verified | [recorded](#command-c7-shared-default-run) |
| 7 | shared-no-default | 37/37 | 0 | 0 | 5.116 s | Verified | [recorded](#command-c7-shared-no-default-run) |
| 7 | signed-acceptances | 4/4 | 0 | 0 | 120.076 s | Verified | [recorded](#command-c7-signed-acceptances-run) |
| 7 | ufo-model-parity | 1/1 | 0 | 378 | 0.620 s | Verified | [recorded](#command-c7-ufo-model-parity-run) |
| 7 | vertex-rules | 1/1 | 0 | 725 | 0.356 s | Verified | [recorded](#command-c7-vertex-rules-run) |

The accepted C3 amendment extends the existing exact source certificate to two fixtures, GL00 and GL04. It compares the complete signed numerator and powered-denominator assignment in a common loop chart, beginning from the complete post-Taylor child source. The retained coupling, physical owner multiplicities and distinct raw vacuum/expansion masses are part of that contract. One test identity exercises the two fixtures; fixture count is not a second test execution. The C3 owning selection and final curated selection must contain `graph::three_d_source::tests::gl00_gl04_planned_lifts_match_post_t_numerators_in_common_loop_coordinates`. Numerical GL00 routes remain complementary end-to-end coverage.

The collector reconciles every executed binary/test identity against its recorded nextest JSON selection. Repeated final-status lines are deduplicated only when identity, outcome, duration and progress index agree. Separate boundaries and feature runs are not added as unique test coverage.

## Final runtime selections

| Selection | Passed/executed | Nonpass | Skipped | Test time | Status | Command |
| --- | --- | --- | --- | --- | --- | --- |
| curated | 2065/2065 | 0 | 297 | 465.870 s | completed | [recorded](#command-final-curated-run) |
| commit5-phases | 212/212 | 0 | 1205 | 30.817 s | completed | [recorded](#command-final-commit5-phases-run) |
| signed-acceptances | 4/4 | 0 | 0 | 120.076 s | completed | [recorded](#command-final-signed-acceptances-run) |
| commit6-vakint-uv | 241/241 | 0 | 619 | 73.977 s | completed | [recorded](#command-final-commit6-vakint-uv-run) |
| scalar-all | 166/166 | 0 | 107 | 922.523 s | completed | [recorded](#command-final-scalar-all-run) |
| shared-default | 37/37 | 0 | 0 | 5.191 s | completed | [recorded](#command-final-shared-default-run) |
| shared-all | 108/108 | 0 | 0 | 21.324 s | completed | [recorded](#command-final-shared-all-run) |
| shared-no-default | 37/37 | 0 | 0 | 5.116 s | completed | [recorded](#command-final-shared-no-default-run) |
| ufo-model-parity | 1/1 | 0 | 378 | 0.620 s | completed | [recorded](#command-final-ufo-model-parity-run) |
| vertex-rules | 1/1 | 0 | 725 | 0.356 s | completed | [recorded](#command-final-vertex-rules-run) |

Complete final-runtime aggregation: **2872 executions across 2207 distinct binary/test identities**. Status counts: `{"PASS": 2872}`. Identities with any nonpass: 0. Distinct execution configurations: 8.

Missing or incomplete required selections: none.

The complete scalar selection executes 26 standard and 140 slow cases (166 total). Its recorded command selects the two scalar namespaces directly in release mode, with ignored tests enabled and four workers; it does not use the wrapper's forced single-worker setting.

## Execution controls and warnings

Completed test contexts record four nextest workers, zero retries, disabled snapshot updates (`INSTA_UPDATE=no`) and a 30 GB process-tree RSS limit. The boundary driver uses four Cargo build jobs. Exact argv, guards and per-job context hashes remain in the evidence. Tests use `--no-fail-fast`; this reports remaining selected outcomes without retrying a failure. Optional feature changes and the release scalar selection use their own recorded Cargo configurations.

Passing Clippy or Cargo exits do not mean warning-free compilation. The recorded runtime driver explicitly exports `VALIDATION_ALLOW_WARNINGS=1`: the curated and scalar selections therefore omit their optional `-Dwarnings` injection. This is an explicit warning-policy exception. The recorded workspace Clippy command also has no `-D warnings` argument. Compiler and dependency warnings remain visible in the raw logs. The following unique warning lines are taken only from hash-matching completed invocation logs:

- `warning: gammalooprs (lib test) generated 1 warning (1 duplicate)`
- `warning: gammalooprs (lib) generated 1 warning`
- `warning: the following packages contain code that will be rejected by a future version of Rust: proc-macro-error2 v2.0.1`
- `warning: very complex type used. Consider factoring parts into type definitions`

## Interrupted or failed prerequisite attempts

This collector snapshot contains the completed boundary runs. The pre-gate startup diagnostic is recorded separately below and contributes no compilation or test outcome.

A separate pre-gate driver startup attempt is retained in `gl00-certificate-v2/diagnostics/initial-functional-driver/diagnostic-and-fix.json`. It concerns the prefix-workspace start precondition and is excluded from production test outcomes. Its raw log and disposition remain in the archive.

## Separate drafting diagnostics and approval

Closeout drafting diagnostics retained for transparency; excluded from reconstructed-boundary and final-runtime pass/nonpass totals. No new commands were executed when this receipt was assembled.

**json-serde-deserialize-bound.** The initial pair-list adapter check failed with E0277 because serde(with) suppressed the inferred Deserialize bound for T. The field-level bound restores that deserialization constraint without restricting GenericAdditionalWeightInfo<T> itself. Added bound(deserialize = "T: Deserialize<'de>") on the serialized field. The subsequent check exited 0; that diagnostic check does not replace fresh reconstructed-boundary validation.

Initial check: exit 101; corrected check: exit 0. The initial failure occurred before runtime testing.

```sh
cargo check -p gammaloop-integration-tests --test test_differential --test test_runs --profile dev-optim --locked
```

Captured diagnostic sources: `/tmp/raised-stack/closeout/json-fix/check.log`, `/tmp/raised-stack/closeout/json-fix/check.watchdog.jsonl`, `/tmp/raised-stack/closeout/json-fix/with-bound/check.log`, `/tmp/raised-stack/closeout/json-fix/with-bound/check.watchdog.jsonl`. Their hashes are pinned in `validation-diagnostic-attempts.json`; the final archive manifest maps each captured path to its actual member without reusing obsolete namespace labels.

**new-opaque-scalar-test-expectations.** Two newly drafted assertions required the parser to reject ambiguous tensor products; the public parser instead preserves their complete opaque scalar function values. Production behavior was unchanged. With explicit approval, the new tests now check the complete preserved values through public parsing and execution. The corrected draft selection passed 44/44; those executions remain separate from the reconstructed stack totals.

The draft run executed 44 tests: 42 passed and 2 failed. The approval receipt records the answer “Yes, check the preserved value”. These were incorrect expectations in two new draft tests; they were not production regressions or failures of the reconstructed stack.

Captured diagnostic sources: `/tmp/stack-closeout-compact/review.json`, `/tmp/stack-closeout-compact/test-correction-approval.json`, `/tmp/stack-closeout-compact/tests.log`, `/tmp/stack-closeout-compact/tests.watchdog.jsonl`, `/tmp/stack-closeout-compact/failing-draft-parsing-test.rs`, `/tmp/stack-closeout-compact/final-tests.log`, `/tmp/stack-closeout-compact/final-tests.watchdog.jsonl`. Their hashes are pinned in `validation-diagnostic-attempts.json`; the final archive manifest maps each captured path to its actual member without reusing obsolete namespace labels.

**draft-watchdog-option.** Watchdog invoked with unsupported --limit-bytes; no cargo check ran. Corrected invocation uses --limit-gb 15 and has separate baseline-check-corrected receipts. This unsupported watchdog option prevented Cargo from starting; it is not a compilation or test verdict.

Captured diagnostic sources: `/tmp/stack-closeout-compact/baseline-check.log`, `/tmp/stack-closeout-compact/review.json`. Their hashes are pinned in `validation-diagnostic-attempts.json`; the final archive manifest maps each captured path to its actual member without reusing obsolete namespace labels.


## Document closeout diagnostics

The archive helper exited 1 only while printing its final size summary: it called the integer `stat().st_size` as a function. Publication, its built-in member verification, and the ledger write had already completed. Separate root and independent checks then verified all published members and source links without overwriting the archive. The invocation error and verification are retained in external `diagnostics/archive-publication-summary-error.json` and `independent-published-archive-review.json` receipts; they are excluded from test totals.

## Scope and limits

- Manual physical PySecDec integrations remain excluded; passing Rust adapters does not certify those integrations.
- Ordinary selections exclude legacy failing-class tests and `aa_aa::important::aa_aa_local_inspect_backend_consistency`. The exact reviewed vertex-rule case is an explicit separate selection.
- Three scalar triangle, double-triangle and TBT reference points remain outside automatic coverage; their decimal targets lack an independently recorded derivation.
- Optional installed-Python subprocess tests are not certified by these runtime selectors. Model parity and the Python API compile check have distinct scopes.
- Sampled numerator grouping and Monte Carlo acceptances are bounded numerical evidence. Local electroweak/Ward checks do not establish complete electroweak LU acceptance; complex-mass cutting rules and arbitrary non-Hermitian interactions remain unsupported.
- General serialized-graph validation and a changed singular-surface contract for the diagnostic f64 evaluator are outside this review.
- This generated page certifies only completed records in its collector inputs. Final PDF compilation/visual inspection, source-to-bookmark closure and final artifact integrity belong to separate external final-closure receipts referencing the archive and final report hashes; those receipts are not embedded in the tracked archive.

## Recorded command appendix

These are the exact command-file contents captured for completed invocations; they are not proposed commands or reconstructed argv. Source file paths identify entries in the fresh closeout capture, whose archive manifest retains their hashes.

<a id="command-c1-fmt"></a>
<details>
<summary>C1 fmt</summary>

Raw command record: `validation-prep/boundary-1/fmt.command`.

```sh
cargo fmt --all -- --check 
```

</details>

<a id="command-c1-check"></a>
<details>
<summary>C1 check</summary>

Raw command record: `validation-prep/boundary-1/check.command`.

```sh
cargo check --workspace --all-targets --locked 
```

</details>

<a id="command-c1-build"></a>
<details>
<summary>C1 build</summary>

Raw command record: `validation-prep/boundary-1/build.command`.

```sh
cargo build --workspace --all-targets --locked --profile dev-optim 
```

</details>

<a id="command-c1-clippy"></a>
<details>
<summary>C1 clippy</summary>

Raw command record: `validation-prep/boundary-1/clippy.command`.

```sh
cargo clippy --workspace --all-targets --locked 
```

</details>

<a id="command-c2-fmt"></a>
<details>
<summary>C2 fmt</summary>

Raw command record: `validation-prep/boundary-2/fmt.command`.

```sh
cargo fmt --all -- --check 
```

</details>

<a id="command-c2-check"></a>
<details>
<summary>C2 check</summary>

Raw command record: `validation-prep/boundary-2/check.command`.

```sh
cargo check --workspace --all-targets --locked 
```

</details>

<a id="command-c2-build"></a>
<details>
<summary>C2 build</summary>

Raw command record: `validation-prep/boundary-2/build.command`.

```sh
cargo build --workspace --all-targets --locked --profile dev-optim 
```

</details>

<a id="command-c2-clippy"></a>
<details>
<summary>C2 clippy</summary>

Raw command record: `validation-prep/boundary-2/clippy.command`.

```sh
cargo clippy --workspace --all-targets --locked 
```

</details>

<a id="command-c2-shared-all"></a>
<details>
<summary>C2 shared-all</summary>

Raw command record: `validation-prep/boundary-2/features/shared-all.command`.

```sh
cargo check -p three-dimensional-reps --all-targets --all-features --locked 
```

</details>

<a id="command-c2-shared-no-default"></a>
<details>
<summary>C2 shared-no-default</summary>

Raw command record: `validation-prep/boundary-2/features/shared-no-default.command`.

```sh
cargo check -p three-dimensional-reps --all-targets --no-default-features --locked 
```

</details>

<a id="command-c3-fmt"></a>
<details>
<summary>C3 fmt</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-3/fmt.command`.

```sh
cargo fmt --all -- --check 
```

</details>

<a id="command-c3-check"></a>
<details>
<summary>C3 check</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-3/check.command`.

```sh
cargo check --workspace --all-targets --locked 
```

</details>

<a id="command-c3-build"></a>
<details>
<summary>C3 build</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-3/build.command`.

```sh
cargo build --workspace --all-targets --locked --profile dev-optim 
```

</details>

<a id="command-c3-clippy"></a>
<details>
<summary>C3 clippy</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-3/clippy.command`.

```sh
cargo clippy --workspace --all-targets --locked 
```

</details>

<a id="command-c3-shared-all"></a>
<details>
<summary>C3 shared-all</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-3/features/shared-all.command`.

```sh
cargo check -p three-dimensional-reps --all-targets --all-features --locked 
```

</details>

<a id="command-c3-shared-no-default"></a>
<details>
<summary>C3 shared-no-default</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-3/features/shared-no-default.command`.

```sh
cargo check -p three-dimensional-reps --all-targets --no-default-features --locked 
```

</details>

<a id="command-c4-fmt"></a>
<details>
<summary>C4 fmt</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-4/fmt.command`.

```sh
cargo fmt --all -- --check 
```

</details>

<a id="command-c4-check"></a>
<details>
<summary>C4 check</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-4/check.command`.

```sh
cargo check --workspace --all-targets --locked 
```

</details>

<a id="command-c4-build"></a>
<details>
<summary>C4 build</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-4/build.command`.

```sh
cargo build --workspace --all-targets --locked --profile dev-optim 
```

</details>

<a id="command-c4-clippy"></a>
<details>
<summary>C4 clippy</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-4/clippy.command`.

```sh
cargo clippy --workspace --all-targets --locked 
```

</details>

<a id="command-c4-shared-all"></a>
<details>
<summary>C4 shared-all</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-4/features/shared-all.command`.

```sh
cargo check -p three-dimensional-reps --all-targets --all-features --locked 
```

</details>

<a id="command-c4-shared-no-default"></a>
<details>
<summary>C4 shared-no-default</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-4/features/shared-no-default.command`.

```sh
cargo check -p three-dimensional-reps --all-targets --no-default-features --locked 
```

</details>

<a id="command-c4-python-api"></a>
<details>
<summary>C4 python-api</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-4/features/python-api.command`.

```sh
cargo check -p gammaloop-api --features python_api --all-targets --locked 
```

</details>

<a id="command-c5-fmt"></a>
<details>
<summary>C5 fmt</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-5/fmt.command`.

```sh
cargo fmt --all -- --check 
```

</details>

<a id="command-c5-check"></a>
<details>
<summary>C5 check</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-5/check.command`.

```sh
cargo check --workspace --all-targets --locked 
```

</details>

<a id="command-c5-build"></a>
<details>
<summary>C5 build</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-5/build.command`.

```sh
cargo build --workspace --all-targets --locked --profile dev-optim 
```

</details>

<a id="command-c5-clippy"></a>
<details>
<summary>C5 clippy</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-5/clippy.command`.

```sh
cargo clippy --workspace --all-targets --locked 
```

</details>

<a id="command-c5-shared-all"></a>
<details>
<summary>C5 shared-all</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-5/features/shared-all.command`.

```sh
cargo check -p three-dimensional-reps --all-targets --all-features --locked 
```

</details>

<a id="command-c5-shared-no-default"></a>
<details>
<summary>C5 shared-no-default</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-5/features/shared-no-default.command`.

```sh
cargo check -p three-dimensional-reps --all-targets --no-default-features --locked 
```

</details>

<a id="command-c6-fmt"></a>
<details>
<summary>C6 fmt</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-6/fmt.command`.

```sh
cargo fmt --all -- --check 
```

</details>

<a id="command-c6-check"></a>
<details>
<summary>C6 check</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-6/check.command`.

```sh
cargo check --workspace --all-targets --locked 
```

</details>

<a id="command-c6-build"></a>
<details>
<summary>C6 build</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-6/build.command`.

```sh
cargo build --workspace --all-targets --locked --profile dev-optim 
```

</details>

<a id="command-c6-clippy"></a>
<details>
<summary>C6 clippy</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-6/clippy.command`.

```sh
cargo clippy --workspace --all-targets --locked 
```

</details>

<a id="command-c6-shared-all"></a>
<details>
<summary>C6 shared-all</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-6/features/shared-all.command`.

```sh
cargo check -p three-dimensional-reps --all-targets --all-features --locked 
```

</details>

<a id="command-c6-shared-no-default"></a>
<details>
<summary>C6 shared-no-default</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-6/features/shared-no-default.command`.

```sh
cargo check -p three-dimensional-reps --all-targets --no-default-features --locked 
```

</details>

<a id="command-c6-vakint-community"></a>
<details>
<summary>C6 vakint-community</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-6/features/vakint-community.command`.

```sh
cargo check --workspace --all-targets --features vakint/symbolica_community_module --locked 
```

</details>

<a id="command-c7-fmt"></a>
<details>
<summary>C7 fmt</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/fmt.command`.

```sh
cargo fmt --all -- --check 
```

</details>

<a id="command-c7-check"></a>
<details>
<summary>C7 check</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/check.command`.

```sh
cargo check --workspace --all-targets --locked 
```

</details>

<a id="command-c7-build"></a>
<details>
<summary>C7 build</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/build.command`.

```sh
cargo build --workspace --all-targets --locked --profile dev-optim 
```

</details>

<a id="command-c7-clippy"></a>
<details>
<summary>C7 clippy</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/clippy.command`.

```sh
cargo clippy --workspace --all-targets --locked 
```

</details>

<a id="command-c7-shared-all"></a>
<details>
<summary>C7 shared-all</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/features/shared-all.command`.

```sh
cargo check -p three-dimensional-reps --all-targets --all-features --locked 
```

</details>

<a id="command-c7-shared-no-default"></a>
<details>
<summary>C7 shared-no-default</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/features/shared-no-default.command`.

```sh
cargo check -p three-dimensional-reps --all-targets --no-default-features --locked 
```

</details>

<a id="command-c1-commit1-foundations-list"></a>
<details>
<summary>C1 commit1-foundations list</summary>

Raw command record: `validation-prep/boundary-1/tests/commit1-foundations.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(package\(=spenso\)\ or\ package\(=spenso-hep-lib\)\ or\ package\(=idenso\)\ or\ package\(=vakint\)\ or\ package\(=linnet\)\ or\ package\(=linnest\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c1-commit1-foundations-run"></a>
<details>
<summary>C1 commit1-foundations run</summary>

Raw command record: `validation-prep/boundary-1/tests/commit1-foundations.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(package\(=spenso\)\ or\ package\(=spenso-hep-lib\)\ or\ package\(=idenso\)\ or\ package\(=vakint\)\ or\ package\(=linnet\)\ or\ package\(=linnest\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c2-shared-all-list"></a>
<details>
<summary>C2 shared-all list</summary>

Raw command record: `validation-prep/boundary-2/tests/shared-all.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --lib --all-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c2-shared-all-run"></a>
<details>
<summary>C2 shared-all run</summary>

Raw command record: `validation-prep/boundary-2/tests/shared-all.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --lib --all-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c2-shared-default-list"></a>
<details>
<summary>C2 shared-default list</summary>

Raw command record: `validation-prep/boundary-2/tests/shared-default.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --lib -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c2-shared-default-run"></a>
<details>
<summary>C2 shared-default run</summary>

Raw command record: `validation-prep/boundary-2/tests/shared-default.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --lib -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c2-shared-no-default-list"></a>
<details>
<summary>C2 shared-no-default list</summary>

Raw command record: `validation-prep/boundary-2/tests/shared-no-default.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --lib --no-default-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c2-shared-no-default-run"></a>
<details>
<summary>C2 shared-no-default run</summary>

Raw command record: `validation-prep/boundary-2/tests/shared-no-default.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --lib --no-default-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c3-commit3-cff-uv-list"></a>
<details>
<summary>C3 commit3-cff-uv list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-3/tests/commit3-cff-uv.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(\(package\(=gammalooprs\)\ and\ \(\(binary\(=gammalooprs\)\ and\ test\(/\^\(cff::\|uv::\|graph::three_d_source::\)/\)\)\ or\ binary\(=test_renormalization\)\)\)\ or\ \(package\(=gammaloop-integration-tests\)\ and\ binary\(=uv\)\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c3-commit3-cff-uv-run"></a>
<details>
<summary>C3 commit3-cff-uv run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-3/tests/commit3-cff-uv.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(\(package\(=gammalooprs\)\ and\ \(\(binary\(=gammalooprs\)\ and\ test\(/\^\(cff::\|uv::\|graph::three_d_source::\)/\)\)\ or\ binary\(=test_renormalization\)\)\)\ or\ \(package\(=gammaloop-integration-tests\)\ and\ binary\(=uv\)\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c3-shared-all-list"></a>
<details>
<summary>C3 shared-all list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-3/tests/shared-all.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --lib --all-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c3-shared-all-run"></a>
<details>
<summary>C3 shared-all run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-3/tests/shared-all.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --lib --all-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c3-shared-default-list"></a>
<details>
<summary>C3 shared-default list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-3/tests/shared-default.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --lib -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c3-shared-default-run"></a>
<details>
<summary>C3 shared-default run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-3/tests/shared-default.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --lib -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c3-shared-no-default-list"></a>
<details>
<summary>C3 shared-no-default list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-3/tests/shared-no-default.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --lib --no-default-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c3-shared-no-default-run"></a>
<details>
<summary>C3 shared-no-default run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-3/tests/shared-no-default.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --lib --no-default-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c4-commit4-scalar-smoke-list"></a>
<details>
<summary>C4 commit4-scalar-smoke list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-4/tests/commit4-scalar-smoke.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(binary\(=test_runs\)\ and\ \(test\(/\^scalar_3l_cross_section_inspects::default_scalar_3l_cross_section_inspects::scalar_3l_cross_section_gl\(00\|04\)_/\)\ or\ test\(=scalar_3l_cross_section_inspects::scalar_3l_gl00_higher_power_cff_is_invariant_under_nonzero_sampling_scale\)\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c4-commit4-scalar-smoke-run"></a>
<details>
<summary>C4 commit4-scalar-smoke run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-4/tests/commit4-scalar-smoke.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(binary\(=test_runs\)\ and\ \(test\(/\^scalar_3l_cross_section_inspects::default_scalar_3l_cross_section_inspects::scalar_3l_cross_section_gl\(00\|04\)_/\)\ or\ test\(=scalar_3l_cross_section_inspects::scalar_3l_gl00_higher_power_cff_is_invariant_under_nonzero_sampling_scale\)\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c4-commit4-workflows-list"></a>
<details>
<summary>C4 commit4-workflows list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-4/tests/commit4-workflows.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(package\(=gammaloop-api\)\ or\ \(package\(=gammalooprs\)\ and\ binary\(=gammalooprs\)\ and\ test\(/\^\(integrands::\(evaluation\|process\)::\|integrate::\|observables::\|uv::profile::\)/\)\)\ or\ \(package\(=gammaloop-integration-tests\)\ and\ \(binary\(/\^\(test_cli\|test_evaluation_api\|test_feyngen\|uv\)\$/\)\ or\ \(binary\(=test_differential\)\ and\ test\(=additional_event_weights_roundtrip_preserves_all_keys_and_values\)\)\ or\ \(binary\(=test_runs\)\ and\ test\(/\^\(inspect\|events\|approach\|integrations\|multi_integrand\|profile_bulk\|runtime_graph_subset\|repeated_masses\|test_3d_reps\|smoke\)::/\)\)\)\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c4-commit4-workflows-run"></a>
<details>
<summary>C4 commit4-workflows run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-4/tests/commit4-workflows.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(package\(=gammaloop-api\)\ or\ \(package\(=gammalooprs\)\ and\ binary\(=gammalooprs\)\ and\ test\(/\^\(integrands::\(evaluation\|process\)::\|integrate::\|observables::\|uv::profile::\)/\)\)\ or\ \(package\(=gammaloop-integration-tests\)\ and\ \(binary\(/\^\(test_cli\|test_evaluation_api\|test_feyngen\|uv\)\$/\)\ or\ \(binary\(=test_differential\)\ and\ test\(=additional_event_weights_roundtrip_preserves_all_keys_and_values\)\)\ or\ \(binary\(=test_runs\)\ and\ test\(/\^\(inspect\|events\|approach\|integrations\|multi_integrand\|profile_bulk\|runtime_graph_subset\|repeated_masses\|test_3d_reps\|smoke\)::/\)\)\)\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c4-shared-all-list"></a>
<details>
<summary>C4 shared-all list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-4/tests/shared-all.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --lib --all-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c4-shared-all-run"></a>
<details>
<summary>C4 shared-all run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-4/tests/shared-all.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --lib --all-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c4-shared-default-list"></a>
<details>
<summary>C4 shared-default list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-4/tests/shared-default.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --lib -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c4-shared-default-run"></a>
<details>
<summary>C4 shared-default run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-4/tests/shared-default.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --lib -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c4-shared-no-default-list"></a>
<details>
<summary>C4 shared-no-default list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-4/tests/shared-no-default.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --lib --no-default-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c4-shared-no-default-run"></a>
<details>
<summary>C4 shared-no-default run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-4/tests/shared-no-default.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --lib --no-default-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c4-signed-acceptances-list"></a>
<details>
<summary>C4 signed-acceptances list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-4/tests/signed-acceptances.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(binary\(/\^\(test_epem_a_ddx_nlo_acceptance\|test_epem_a_ttx_nlo_acceptance\|test_gamma_star_ddx_nlo_acceptance\|test_gamma_star_ttx_nlo_acceptance\)\$/\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c4-signed-acceptances-run"></a>
<details>
<summary>C4 signed-acceptances run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-4/tests/signed-acceptances.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(binary\(/\^\(test_epem_a_ddx_nlo_acceptance\|test_epem_a_ttx_nlo_acceptance\|test_gamma_star_ddx_nlo_acceptance\|test_gamma_star_ttx_nlo_acceptance\)\$/\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c4-source-certificates-list"></a>
<details>
<summary>C4 source-certificates list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-4/tests/source-certificates.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(package\(=gammalooprs\)\ and\ binary\(=gammalooprs\)\ and\ test\(=graph::three_d_source::tests::gl00_gl04_planned_lifts_match_post_t_numerators_in_common_loop_coordinates\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c4-source-certificates-run"></a>
<details>
<summary>C4 source-certificates run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-4/tests/source-certificates.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(package\(=gammalooprs\)\ and\ binary\(=gammalooprs\)\ and\ test\(=graph::three_d_source::tests::gl00_gl04_planned_lifts_match_post_t_numerators_in_common_loop_coordinates\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c5-commit5-phases-list"></a>
<details>
<summary>C5 commit5-phases list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-5/tests/commit5-phases.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(\(package\(=gammaloop-api\)\ and\ test\(/\^\(commands::generate::\|commands::import::model::\)/\)\)\ or\ \(package\(=gammalooprs\)\ and\ \(\(binary\(=gammalooprs\)\ and\ test\(/\^\(cff::\|feyngen::\|model::\|processes::cross_section::\|graph::feynman_graph::\|graph::lmb::\|uv::export::\|uv::settings::\)/\)\)\ or\ binary\(=test_renormalization\)\)\)\ or\ \(package\(=gammaloop-integration-tests\)\ and\ \(binary\(=test_feyngen\)\ or\ \(binary\(=test_runs\)\ and\ test\(/\^\(amplitude_phase_conventions\|scalar_phase_conventions\|scalar_virtual_phase_conventions\)::/\)\)\)\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c5-commit5-phases-run"></a>
<details>
<summary>C5 commit5-phases run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-5/tests/commit5-phases.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(\(package\(=gammaloop-api\)\ and\ test\(/\^\(commands::generate::\|commands::import::model::\)/\)\)\ or\ \(package\(=gammalooprs\)\ and\ \(\(binary\(=gammalooprs\)\ and\ test\(/\^\(cff::\|feyngen::\|model::\|processes::cross_section::\|graph::feynman_graph::\|graph::lmb::\|uv::export::\|uv::settings::\)/\)\)\ or\ binary\(=test_renormalization\)\)\)\ or\ \(package\(=gammaloop-integration-tests\)\ and\ \(binary\(=test_feyngen\)\ or\ \(binary\(=test_runs\)\ and\ test\(/\^\(amplitude_phase_conventions\|scalar_phase_conventions\|scalar_virtual_phase_conventions\)::/\)\)\)\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c5-shared-all-list"></a>
<details>
<summary>C5 shared-all list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-5/tests/shared-all.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --lib --all-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c5-shared-all-run"></a>
<details>
<summary>C5 shared-all run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-5/tests/shared-all.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --lib --all-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c5-shared-default-list"></a>
<details>
<summary>C5 shared-default list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-5/tests/shared-default.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --lib -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c5-shared-default-run"></a>
<details>
<summary>C5 shared-default run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-5/tests/shared-default.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --lib -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c5-shared-no-default-list"></a>
<details>
<summary>C5 shared-no-default list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-5/tests/shared-no-default.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --lib --no-default-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c5-shared-no-default-run"></a>
<details>
<summary>C5 shared-no-default run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-5/tests/shared-no-default.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --lib --no-default-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c5-signed-acceptances-list"></a>
<details>
<summary>C5 signed-acceptances list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-5/tests/signed-acceptances.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(binary\(/\^\(test_epem_a_ddx_nlo_acceptance\|test_epem_a_ttx_nlo_acceptance\|test_gamma_star_ddx_nlo_acceptance\|test_gamma_star_ttx_nlo_acceptance\)\$/\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c5-signed-acceptances-run"></a>
<details>
<summary>C5 signed-acceptances run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-5/tests/signed-acceptances.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(binary\(/\^\(test_epem_a_ddx_nlo_acceptance\|test_epem_a_ttx_nlo_acceptance\|test_gamma_star_ddx_nlo_acceptance\|test_gamma_star_ttx_nlo_acceptance\)\$/\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c5-source-certificates-list"></a>
<details>
<summary>C5 source-certificates list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-5/tests/source-certificates.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(package\(=gammalooprs\)\ and\ binary\(=gammalooprs\)\ and\ test\(=graph::three_d_source::tests::gl00_gl04_planned_lifts_match_post_t_numerators_in_common_loop_coordinates\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c5-source-certificates-run"></a>
<details>
<summary>C5 source-certificates run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-5/tests/source-certificates.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(package\(=gammalooprs\)\ and\ binary\(=gammalooprs\)\ and\ test\(=graph::three_d_source::tests::gl00_gl04_planned_lifts_match_post_t_numerators_in_common_loop_coordinates\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c5-vertex-rules-list"></a>
<details>
<summary>C5 vertex-rules list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-5/tests/vertex-rules.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --ignore-default-filter --run-ignored all --workspace -E \(package\(=gammalooprs\)\ and\ binary\(=gammalooprs\)\ and\ test\(=graph::parse::tests::failing::vertex_rules\)\)\ and\ \(all\(\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c5-vertex-rules-run"></a>
<details>
<summary>C5 vertex-rules run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-5/tests/vertex-rules.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --ignore-default-filter --run-ignored all --workspace -E \(package\(=gammalooprs\)\ and\ binary\(=gammalooprs\)\ and\ test\(=graph::parse::tests::failing::vertex_rules\)\)\ and\ \(all\(\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c6-commit6-vakint-uv-list"></a>
<details>
<summary>C6 commit6-vakint-uv list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-6/tests/commit6-vakint-uv.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(package\(=vakint\)\ or\ \(package\(=gammalooprs\)\ and\ \(\(binary\(=gammalooprs\)\ and\ test\(/\^uv::/\)\)\ or\ binary\(=test_renormalization\)\)\)\ or\ \(package\(=gammaloop-integration-tests\)\ and\ binary\(=uv\)\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c6-commit6-vakint-uv-run"></a>
<details>
<summary>C6 commit6-vakint-uv run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-6/tests/commit6-vakint-uv.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(package\(=vakint\)\ or\ \(package\(=gammalooprs\)\ and\ \(\(binary\(=gammalooprs\)\ and\ test\(/\^uv::/\)\)\ or\ binary\(=test_renormalization\)\)\)\ or\ \(package\(=gammaloop-integration-tests\)\ and\ binary\(=uv\)\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c6-shared-all-list"></a>
<details>
<summary>C6 shared-all list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-6/tests/shared-all.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --lib --all-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c6-shared-all-run"></a>
<details>
<summary>C6 shared-all run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-6/tests/shared-all.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --lib --all-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c6-shared-default-list"></a>
<details>
<summary>C6 shared-default list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-6/tests/shared-default.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --lib -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c6-shared-default-run"></a>
<details>
<summary>C6 shared-default run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-6/tests/shared-default.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --lib -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c6-shared-no-default-list"></a>
<details>
<summary>C6 shared-no-default list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-6/tests/shared-no-default.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --lib --no-default-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c6-shared-no-default-run"></a>
<details>
<summary>C6 shared-no-default run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-6/tests/shared-no-default.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --lib --no-default-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c6-source-certificates-list"></a>
<details>
<summary>C6 source-certificates list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-6/tests/source-certificates.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(package\(=gammalooprs\)\ and\ binary\(=gammalooprs\)\ and\ test\(=graph::three_d_source::tests::gl00_gl04_planned_lifts_match_post_t_numerators_in_common_loop_coordinates\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c6-source-certificates-run"></a>
<details>
<summary>C6 source-certificates run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-6/tests/source-certificates.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(package\(=gammalooprs\)\ and\ binary\(=gammalooprs\)\ and\ test\(=graph::three_d_source::tests::gl00_gl04_planned_lifts_match_post_t_numerators_in_common_loop_coordinates\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c7-commit5-phases-list"></a>
<details>
<summary>C7 commit5-phases list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/commit5-phases.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(\(package\(=gammaloop-api\)\ and\ test\(/\^\(commands::generate::\|commands::import::model::\)/\)\)\ or\ \(package\(=gammalooprs\)\ and\ \(\(binary\(=gammalooprs\)\ and\ test\(/\^\(cff::\|feyngen::\|model::\|processes::cross_section::\|graph::feynman_graph::\|graph::lmb::\|uv::export::\|uv::settings::\)/\)\)\ or\ binary\(=test_renormalization\)\)\)\ or\ \(package\(=gammaloop-integration-tests\)\ and\ \(binary\(=test_feyngen\)\ or\ \(binary\(=test_runs\)\ and\ test\(/\^\(amplitude_phase_conventions\|scalar_phase_conventions\|scalar_virtual_phase_conventions\)::/\)\)\)\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c7-commit5-phases-run"></a>
<details>
<summary>C7 commit5-phases run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/commit5-phases.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(\(package\(=gammaloop-api\)\ and\ test\(/\^\(commands::generate::\|commands::import::model::\)/\)\)\ or\ \(package\(=gammalooprs\)\ and\ \(\(binary\(=gammalooprs\)\ and\ test\(/\^\(cff::\|feyngen::\|model::\|processes::cross_section::\|graph::feynman_graph::\|graph::lmb::\|uv::export::\|uv::settings::\)/\)\)\ or\ binary\(=test_renormalization\)\)\)\ or\ \(package\(=gammaloop-integration-tests\)\ and\ \(binary\(=test_feyngen\)\ or\ \(binary\(=test_runs\)\ and\ test\(/\^\(amplitude_phase_conventions\|scalar_phase_conventions\|scalar_virtual_phase_conventions\)::/\)\)\)\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c7-commit6-vakint-uv-list"></a>
<details>
<summary>C7 commit6-vakint-uv list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/commit6-vakint-uv.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(package\(=vakint\)\ or\ \(package\(=gammalooprs\)\ and\ \(\(binary\(=gammalooprs\)\ and\ test\(/\^uv::/\)\)\ or\ binary\(=test_renormalization\)\)\)\ or\ \(package\(=gammaloop-integration-tests\)\ and\ binary\(=uv\)\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c7-commit6-vakint-uv-run"></a>
<details>
<summary>C7 commit6-vakint-uv run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/commit6-vakint-uv.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(package\(=vakint\)\ or\ \(package\(=gammalooprs\)\ and\ \(\(binary\(=gammalooprs\)\ and\ test\(/\^uv::/\)\)\ or\ binary\(=test_renormalization\)\)\)\ or\ \(package\(=gammaloop-integration-tests\)\ and\ binary\(=uv\)\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c7-curated-list"></a>
<details>
<summary>C7 curated list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/curated.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked -p gammaloop-api -p gammalooprs -p idenso -p linnest -p linnet -p spenso -p spenso-hep-lib -p spenso-macros -p three-dimensional-reps -p vakint -p gammaloop-integration-tests -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c7-curated-run"></a>
<details>
<summary>C7 curated run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/curated.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked -p gammaloop-api -p gammalooprs -p idenso -p linnest -p linnet -p spenso -p spenso-hep-lib -p spenso-macros -p three-dimensional-reps -p vakint -p gammaloop-integration-tests -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c7-scalar-all-list"></a>
<details>
<summary>C7 scalar-all list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/scalar-all.list.command`.

```sh
cargo nextest list --release --profile test_gammaloop --locked --test test_runs --ignore-default-filter --run-ignored all --workspace -E \(test\(/\^scalar_3l_cross_section_inspects::default_scalar_3l_cross_section_inspects::/\)\ or\ test\(/\^scalar_3l_cross_section_inspects::slow::/\)\)\ and\ \(not\ test\(/\(\^\|::\)failing::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c7-scalar-all-run"></a>
<details>
<summary>C7 scalar-all run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/scalar-all.run.command`.

```sh
cargo nextest run --release --profile test_gammaloop --locked --test test_runs --ignore-default-filter --run-ignored all --workspace -E \(test\(/\^scalar_3l_cross_section_inspects::default_scalar_3l_cross_section_inspects::/\)\ or\ test\(/\^scalar_3l_cross_section_inspects::slow::/\)\)\ and\ \(not\ test\(/\(\^\|::\)failing::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c7-shared-all-list"></a>
<details>
<summary>C7 shared-all list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/shared-all.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --lib --all-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c7-shared-all-run"></a>
<details>
<summary>C7 shared-all run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/shared-all.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --lib --all-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c7-shared-default-list"></a>
<details>
<summary>C7 shared-default list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/shared-default.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --lib -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c7-shared-default-run"></a>
<details>
<summary>C7 shared-default run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/shared-default.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --lib -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c7-shared-no-default-list"></a>
<details>
<summary>C7 shared-no-default list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/shared-no-default.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --lib --no-default-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c7-shared-no-default-run"></a>
<details>
<summary>C7 shared-no-default run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/shared-no-default.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --lib --no-default-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c7-signed-acceptances-list"></a>
<details>
<summary>C7 signed-acceptances list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/signed-acceptances.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(binary\(/\^\(test_epem_a_ddx_nlo_acceptance\|test_epem_a_ttx_nlo_acceptance\|test_gamma_star_ddx_nlo_acceptance\|test_gamma_star_ttx_nlo_acceptance\)\$/\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c7-signed-acceptances-run"></a>
<details>
<summary>C7 signed-acceptances run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/signed-acceptances.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(binary\(/\^\(test_epem_a_ddx_nlo_acceptance\|test_epem_a_ttx_nlo_acceptance\|test_gamma_star_ddx_nlo_acceptance\|test_gamma_star_ttx_nlo_acceptance\)\$/\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c7-ufo-model-parity-list"></a>
<details>
<summary>C7 ufo-model-parity list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/ufo-model-parity.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --features gammaloop-api/ufo_support --workspace -E \(package\(=gammaloop-api\)\ and\ binary\(=gammaloop_api\)\ and\ test\(=commands::import::model::tests::fresh_sm_ufo_preserves_virtual_gauge_and_covariant_cut_states\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c7-ufo-model-parity-run"></a>
<details>
<summary>C7 ufo-model-parity run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/ufo-model-parity.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --features gammaloop-api/ufo_support --workspace -E \(package\(=gammaloop-api\)\ and\ binary\(=gammaloop_api\)\ and\ test\(=commands::import::model::tests::fresh_sm_ufo_preserves_virtual_gauge_and_covariant_cut_states\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-c7-vertex-rules-list"></a>
<details>
<summary>C7 vertex-rules list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/vertex-rules.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --ignore-default-filter --run-ignored all --workspace -E \(package\(=gammalooprs\)\ and\ binary\(=gammalooprs\)\ and\ test\(=graph::parse::tests::failing::vertex_rules\)\)\ and\ \(all\(\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-c7-vertex-rules-run"></a>
<details>
<summary>C7 vertex-rules run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/vertex-rules.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --ignore-default-filter --run-ignored all --workspace -E \(package\(=gammalooprs\)\ and\ binary\(=gammalooprs\)\ and\ test\(=graph::parse::tests::failing::vertex_rules\)\)\ and\ \(all\(\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-final-curated-list"></a>
<details>
<summary>Final curated list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/curated.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked -p gammaloop-api -p gammalooprs -p idenso -p linnest -p linnet -p spenso -p spenso-hep-lib -p spenso-macros -p three-dimensional-reps -p vakint -p gammaloop-integration-tests -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-final-curated-run"></a>
<details>
<summary>Final curated run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/curated.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked -p gammaloop-api -p gammalooprs -p idenso -p linnest -p linnet -p spenso -p spenso-hep-lib -p spenso-macros -p three-dimensional-reps -p vakint -p gammaloop-integration-tests -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-final-commit5-phases-list"></a>
<details>
<summary>Final commit5-phases list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/commit5-phases.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(\(package\(=gammaloop-api\)\ and\ test\(/\^\(commands::generate::\|commands::import::model::\)/\)\)\ or\ \(package\(=gammalooprs\)\ and\ \(\(binary\(=gammalooprs\)\ and\ test\(/\^\(cff::\|feyngen::\|model::\|processes::cross_section::\|graph::feynman_graph::\|graph::lmb::\|uv::export::\|uv::settings::\)/\)\)\ or\ binary\(=test_renormalization\)\)\)\ or\ \(package\(=gammaloop-integration-tests\)\ and\ \(binary\(=test_feyngen\)\ or\ \(binary\(=test_runs\)\ and\ test\(/\^\(amplitude_phase_conventions\|scalar_phase_conventions\|scalar_virtual_phase_conventions\)::/\)\)\)\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-final-commit5-phases-run"></a>
<details>
<summary>Final commit5-phases run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/commit5-phases.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(\(package\(=gammaloop-api\)\ and\ test\(/\^\(commands::generate::\|commands::import::model::\)/\)\)\ or\ \(package\(=gammalooprs\)\ and\ \(\(binary\(=gammalooprs\)\ and\ test\(/\^\(cff::\|feyngen::\|model::\|processes::cross_section::\|graph::feynman_graph::\|graph::lmb::\|uv::export::\|uv::settings::\)/\)\)\ or\ binary\(=test_renormalization\)\)\)\ or\ \(package\(=gammaloop-integration-tests\)\ and\ \(binary\(=test_feyngen\)\ or\ \(binary\(=test_runs\)\ and\ test\(/\^\(amplitude_phase_conventions\|scalar_phase_conventions\|scalar_virtual_phase_conventions\)::/\)\)\)\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-final-signed-acceptances-list"></a>
<details>
<summary>Final signed-acceptances list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/signed-acceptances.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(binary\(/\^\(test_epem_a_ddx_nlo_acceptance\|test_epem_a_ttx_nlo_acceptance\|test_gamma_star_ddx_nlo_acceptance\|test_gamma_star_ttx_nlo_acceptance\)\$/\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-final-signed-acceptances-run"></a>
<details>
<summary>Final signed-acceptances run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/signed-acceptances.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(binary\(/\^\(test_epem_a_ddx_nlo_acceptance\|test_epem_a_ttx_nlo_acceptance\|test_gamma_star_ddx_nlo_acceptance\|test_gamma_star_ttx_nlo_acceptance\)\$/\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-final-commit6-vakint-uv-list"></a>
<details>
<summary>Final commit6-vakint-uv list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/commit6-vakint-uv.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(package\(=vakint\)\ or\ \(package\(=gammalooprs\)\ and\ \(\(binary\(=gammalooprs\)\ and\ test\(/\^uv::/\)\)\ or\ binary\(=test_renormalization\)\)\)\ or\ \(package\(=gammaloop-integration-tests\)\ and\ binary\(=uv\)\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-final-commit6-vakint-uv-run"></a>
<details>
<summary>Final commit6-vakint-uv run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/commit6-vakint-uv.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E \(package\(=vakint\)\ or\ \(package\(=gammalooprs\)\ and\ \(\(binary\(=gammalooprs\)\ and\ test\(/\^uv::/\)\)\ or\ binary\(=test_renormalization\)\)\)\ or\ \(package\(=gammaloop-integration-tests\)\ and\ binary\(=uv\)\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-final-scalar-all-list"></a>
<details>
<summary>Final scalar-all list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/scalar-all.list.command`.

```sh
cargo nextest list --release --profile test_gammaloop --locked --test test_runs --ignore-default-filter --run-ignored all --workspace -E \(test\(/\^scalar_3l_cross_section_inspects::default_scalar_3l_cross_section_inspects::/\)\ or\ test\(/\^scalar_3l_cross_section_inspects::slow::/\)\)\ and\ \(not\ test\(/\(\^\|::\)failing::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-final-scalar-all-run"></a>
<details>
<summary>Final scalar-all run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/scalar-all.run.command`.

```sh
cargo nextest run --release --profile test_gammaloop --locked --test test_runs --ignore-default-filter --run-ignored all --workspace -E \(test\(/\^scalar_3l_cross_section_inspects::default_scalar_3l_cross_section_inspects::/\)\ or\ test\(/\^scalar_3l_cross_section_inspects::slow::/\)\)\ and\ \(not\ test\(/\(\^\|::\)failing::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-final-shared-default-list"></a>
<details>
<summary>Final shared-default list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/shared-default.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --lib -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-final-shared-default-run"></a>
<details>
<summary>Final shared-default run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/shared-default.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --lib -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-final-shared-all-list"></a>
<details>
<summary>Final shared-all list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/shared-all.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --lib --all-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-final-shared-all-run"></a>
<details>
<summary>Final shared-all run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/shared-all.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --lib --all-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-final-shared-no-default-list"></a>
<details>
<summary>Final shared-no-default list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/shared-no-default.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --lib --no-default-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-final-shared-no-default-run"></a>
<details>
<summary>Final shared-no-default run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/shared-no-default.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --lib --no-default-features -p three-dimensional-reps -E \(all\(\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-final-ufo-model-parity-list"></a>
<details>
<summary>Final ufo-model-parity list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/ufo-model-parity.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --features gammaloop-api/ufo_support --workspace -E \(package\(=gammaloop-api\)\ and\ binary\(=gammaloop_api\)\ and\ test\(=commands::import::model::tests::fresh_sm_ufo_preserves_virtual_gauge_and_covariant_cut_states\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-final-ufo-model-parity-run"></a>
<details>
<summary>Final ufo-model-parity run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/ufo-model-parity.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --features gammaloop-api/ufo_support --workspace -E \(package\(=gammaloop-api\)\ and\ binary\(=gammaloop_api\)\ and\ test\(=commands::import::model::tests::fresh_sm_ufo_preserves_virtual_gauge_and_covariant_cut_states\)\)\ and\ \(not\ test\(/\(\^\|::\)\(slow\|failing\)::/\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

<a id="command-final-vertex-rules-list"></a>
<details>
<summary>Final vertex-rules list</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/vertex-rules.list.command`.

```sh
cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --ignore-default-filter --run-ignored all --workspace -E \(package\(=gammalooprs\)\ and\ binary\(=gammalooprs\)\ and\ test\(=graph::parse::tests::failing::vertex_rules\)\)\ and\ \(all\(\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --message-format json 
```

</details>

<a id="command-final-vertex-rules-run"></a>
<details>
<summary>Final vertex-rules run</summary>

Raw command record: `gl00-certificate-v2/validation-prep/boundary-7/tests/vertex-rules.run.command`.

```sh
cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --ignore-default-filter --run-ignored all --workspace -E \(package\(=gammalooprs\)\ and\ binary\(=gammalooprs\)\ and\ test\(=graph::parse::tests::failing::vertex_rules\)\)\ and\ \(all\(\)\)\ and\ \(not\ test\(/\^aa_aa::important::aa_aa_local_inspect_backend_consistency\$/\)\)\ and\ \(not\ \(binary\(/pysecdec/\)\ or\ test\(=test_integrate_1l_decorated_indices_pysecdec\)\)\) --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never 
```

</details>

## Collector inputs and durable evidence

The [fresh archive](raised-energy-cff-closeout-evidence.tar.gz) and [receipt](raised-energy-cff-closeout-evidence.json) are the durable evidence targets. The receipt records payload paths and hashes. C1/C2 use the retained certificates only because their full revision and tree are unchanged. Initial a474 C3–C7 and runtime evidence remains preserved separately and does not certify the amended source.

| Input | SHA256 |
| --- | --- |
| `gl00-certificate-v2/validation-prep/boundary-evidence-complete.json` | `1747059f8b30b7518a0992beaab28fc086aeab7333143a48eb24ec0b88ffc9a4` |
| `gl00-certificate-v2/reconstructed-candidate.json` | `7b18a59891d084875c7c0843a3bfe5714ce239d3ed51c10f785039c054cbccf7` |
| `gl00-certificate-v2/amendment.json` | `1609514078a49453593d452dc381b6eac0a2c50a0400c24020f48dcbd3a7c897` |
| `gl00-certificate-v2/validation-inputs.json` | `b9740c1a68baaeef17cfe860b927ea7766d65cfa93bf9044adab6dd79c5a76a1` |
| `validation-diagnostic-attempts.json` | `a5397cec275ea8f9b28e51d8240fbf3453a41b4b82083dea48b824db0600404f` |
| `gl00-certificate-v2/validation-prep/final-runtime-job` | `439b183a02984ebc11d29c4a9ef33a3eff85bd58c463597588746aaa2e7a6de9` |
| `gl00-certificate-v2/validation-prep/remaining-validation` | `0c88375d7067be00458b248ab25bbcd79420442e2522bba71bc2a14a1bd8078c` |
| `gl00-certificate-v2/validation-prep/runtime-evidence-complete.json` | `a1c47b8f2077c4693addc75638f37dd1d1c99978a9b29260dd808bbcfe1f6de0` |
| `gl00-certificate-v2/runtime-candidate.json` | `6b7b4556c202c305e43241e682eb459ce533fc03c3cc48efb26624d6ffc56e86` |
| `gl00-certificate-v2/diagnostics/initial-functional-driver/diagnostic-and-fix.json` | `9c1babb009b59020c0ed7b6e24a7f26f82b15886558196fcb248aa53509bb4b9` |

Regenerate this page with the bundled `refresh-validation-report.py` using a current boundary snapshot and, when available, the final runtime snapshot. Collection and rendering execute no tests.
