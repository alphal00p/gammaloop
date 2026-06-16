# Reconstructed-stack validation

Current generated status: **COMPLETE**. Completed commands must match the exact planned argument vector and source revision. Reuse from a tested C8 candidate requires an explicit verified document-equivalence certificate; executed revisions and counts are never relabelled. The latest matching attempt determines its status; an earlier passing attempt cannot hide a later failure. Unfinished, unaudited and stale-source receipts remain explicit.

| Commit | Format | Check | All-target build | Clippy | Other planned stages |
| --- | --- | --- | --- | --- | --- |
| C1 | PASS | PASS | PASS | PASS | 8 pass |
| C2 | PASS | PASS | PASS | PASS | 8 pass |
| C3 | PASS | PASS | PASS | PASS | 2 pass |
| C4 | PASS | PASS | PASS | PASS | 4 pass |
| C5 | PASS | PASS | PASS | PASS | 8 pass |
| C6 | PASS | PASS | PASS | PASS | 2 pass |
| C7 | PASS | PASS | PASS | PASS | 4 pass |
| C8 | PASS | PASS | PASS | PASS | 25 pass |

The final runtime currently has 10/10 audited passing selections: 2706 executions and 2279 distinct binary/test identities. These are partial collector values until every final selection completes. They exclude the immutable incoming baseline and owning-prefix runs; no imported or earlier-prefix suite total is added.

## Final runtime selections

| Selection | Status | Selected / passed | Excluded | Receipt |
| --- | --- | --- | --- | --- |
| curated | PASS | 2135 / 2135 | 300 | c8-results/curated-v5/result.json |
| scalar-all | PASS | 167 / 167 | 106 | c8-results/scalar-all-v2/result.json |
| vertex-rules | PASS | 1 / 1 | 756 | c8-results/vertex-rules-v2/result.json |
| ufo-model-parity | PASS | 1 / 1 | 378 | c8-results/ufo-model-parity-v2/result.json |
| physical-local-uv-routes | PASS | 2 / 2 | 32 | c8-results/physical-local-uv-routes/result.json |
| spenso-shadowing-tests | PASS | 168 / 168 | 0 | c8-results/spenso-shadowing-tests/result.json |
| spenso-no-default-tests | PASS | 50 / 50 | 0 | c8-results/spenso-no-default-tests/result.json |
| three-dimensional-reps-default-tests | PASS | 37 / 37 | 0 | c8-results/three-dimensional-reps-default-tests/result.json |
| three-dimensional-reps-all-tests | PASS | 108 / 108 | 0 | c8-results/three-dimensional-reps-all-tests/result.json |
| three-dimensional-reps-no-default-tests | PASS | 37 / 37 | 0 | c8-results/three-dimensional-reps-no-default-tests/result.json |

## Current boundary receipts

| Boundary / stage | Status | Executed revision / final document target | Seconds | Receipt |
| --- | --- | --- | --- | --- |
| C1/fmt | PASS | `665658b16898` | 4.146 | c1-results/fmt-v4/result.json |
| C1/check | PASS | `665658b16898` | 2.788 | c1-results/check-v3/result.json |
| C1/build | PASS | `665658b16898` | 43.730 | c1-results/build-v2/result.json |
| C1/clippy | PASS | `665658b16898` | 10.927 | c1-results/clippy-v2/result.json |
| C1/tensor-foundations-v2-list | PASS | `665658b16898` | 1.379 | c1-results/tensor-foundations-v2-list/result.json |
| C1/tensor-foundations-v2 | PASS | `665658b16898` | 12.267 | c1-results/tensor-foundations-v2/result.json |
| C1/spenso-shadowing-check | PASS | `665658b16898` | 10.716 | c1-results/spenso-shadowing-check/result.json |
| C1/spenso-shadowing-tests-list | PASS | `665658b16898` | 75.409 | c1-results/spenso-shadowing-tests-list/result.json |
| C1/spenso-shadowing-tests | PASS | `665658b16898` | 3.615 | c1-results/spenso-shadowing-tests/result.json |
| C1/spenso-no-default-check | PASS | `665658b16898` | 9.325 | c1-results/spenso-no-default-check/result.json |
| C1/spenso-no-default-tests-list | PASS | `665658b16898` | 27.221 | c1-results/spenso-no-default-tests-list/result.json |
| C1/spenso-no-default-tests | PASS | `665658b16898` | 1.391 | c1-results/spenso-no-default-tests/result.json |
| C2/fmt | PASS | `a9397ee0f86b` | 4.081 | c2-results/fmt/result.json |
| C2/check | PASS | `a9397ee0f86b` | 3.629 | c2-results/check/result.json |
| C2/build | PASS | `a9397ee0f86b` | 9.269 | c2-results/build/result.json |
| C2/clippy | PASS | `a9397ee0f86b` | 31.675 | c2-results/clippy/result.json |
| C2/shared-cff-default-list | PASS | `a9397ee0f86b` | 8.163 | c2-results/shared-cff-default-list/result.json |
| C2/shared-cff-default | PASS | `a9397ee0f86b` | 6.776 | c2-results/shared-cff-default/result.json |
| C2/three-dimensional-reps-all-check | PASS | `a9397ee0f86b` | 4.101 | c2-results/three-dimensional-reps-all-check/result.json |
| C2/three-dimensional-reps-all-tests-list | PASS | `a9397ee0f86b` | 10.295 | c2-results/three-dimensional-reps-all-tests-list/result.json |
| C2/three-dimensional-reps-all-tests | PASS | `a9397ee0f86b` | 25.249 | c2-results/three-dimensional-reps-all-tests/result.json |
| C2/three-dimensional-reps-no-default-check | PASS | `a9397ee0f86b` | 3.151 | c2-results/three-dimensional-reps-no-default-check/result.json |
| C2/three-dimensional-reps-no-default-tests-list | PASS | `a9397ee0f86b` | 7.736 | c2-results/three-dimensional-reps-no-default-tests-list/result.json |
| C2/three-dimensional-reps-no-default-tests | PASS | `a9397ee0f86b` | 6.803 | c2-results/three-dimensional-reps-no-default-tests/result.json |
| C3/fmt | PASS | `3ea313789a1a` | 5.457 | c3-results/fmt-v2/result.json |
| C3/check | PASS | `3ea313789a1a` | 10.601 | c3-results/check-v2/result.json |
| C3/build | PASS | `3ea313789a1a` | 526.097 | c3-results/build-v2/result.json |
| C3/clippy | PASS | `3ea313789a1a` | 42.104 | c3-results/clippy-v2/result.json |
| C3/commit3-cff-uv-list | PASS | `3ea313789a1a` | 1.435 | c3-results/commit3-cff-uv-list/result.json |
| C3/commit3-cff-uv | PASS | `3ea313789a1a` | 144.107 | c3-results/commit3-cff-uv/result.json |
| C4/fmt | PASS | `04637884f24f` | 5.462 | c4-results/fmt/result.json |
| C4/check | PASS | `04637884f24f` | 32.840 | c4-results/check/result.json |
| C4/build | PASS | `04637884f24f` | 1361.123 | c4-results/build/result.json |
| C4/clippy | PASS | `04637884f24f` | 49.413 | c4-results/clippy/result.json |
| C4/commit4-workflows-list | PASS | `04637884f24f` | 1507.575 | c4-results/commit4-workflows-v2-list/result.json |
| C4/commit4-workflows | PASS | `04637884f24f` | 188.640 | c4-results/commit4-workflows-v2/result.json |
| C4/commit4-scalar-smoke-list | PASS | `04637884f24f` | 1.982 | c4-results/commit4-scalar-smoke-list/result.json |
| C4/commit4-scalar-smoke | PASS | `04637884f24f` | 114.178 | c4-results/commit4-scalar-smoke/result.json |
| C5/fmt | PASS | `1f2cf6d8236d` | 4.957 | c5-results/fmt/result.json |
| C5/check | PASS | `1f2cf6d8236d` | 34.345 | c5-results/check/result.json |
| C5/build | PASS | `1f2cf6d8236d` | 1278.204 | c5-results/build/result.json |
| C5/clippy | PASS | `1f2cf6d8236d` | 39.038 | c5-results/clippy/result.json |
| C5/commit5-phases-list | PASS | `1f2cf6d8236d` | 2.065 | c5-results/commit5-phases-list/result.json |
| C5/commit5-phases | PASS | `1f2cf6d8236d` | 34.035 | c5-results/commit5-phases/result.json |
| C5/signed-acceptances-list | PASS | `1f2cf6d8236d` | 1.465 | c5-results/signed-acceptances-list/result.json |
| C5/signed-acceptances | PASS | `1f2cf6d8236d` | 125.464 | c5-results/signed-acceptances/result.json |
| C5/vertex-rules-list | PASS | `1f2cf6d8236d` | 1.418 | c5-results/vertex-rules-list/result.json |
| C5/vertex-rules | PASS | `1f2cf6d8236d` | 1.824 | c5-results/vertex-rules/result.json |
| C5/ufo-model-parity-list | PASS | `1f2cf6d8236d` | 465.673 | c5-results/ufo-model-parity-list/result.json |
| C5/ufo-model-parity | PASS | `1f2cf6d8236d` | 1.832 | c5-results/ufo-model-parity/result.json |
| C6/fmt | PASS | `8e4e4664f3f9` | 5.253 | c6-results/fmt/result.json |
| C6/check | PASS | `8e4e4664f3f9` | 28.346 | c6-results/check/result.json |
| C6/build | PASS | `8e4e4664f3f9` | 1391.270 | c6-results/build/result.json |
| C6/clippy | PASS | `8e4e4664f3f9` | 45.579 | c6-results/clippy/result.json |
| C6/commit6-vakint-uv-list | PASS | `8e4e4664f3f9` | 1.420 | c6-results/commit6-vakint-uv-list/result.json |
| C6/commit6-vakint-uv | PASS | `8e4e4664f3f9` | 79.734 | c6-results/commit6-vakint-uv/result.json |
| C7/fmt | PASS | `6d7a9220129c` | 5.916 | c7-results/fmt-v3/result.json |
| C7/check | PASS | `6d7a9220129c` | 12.614 | c7-results/check-v3/result.json |
| C7/build | PASS | `6d7a9220129c` | 1192.457 | c7-results/build-v2/result.json |
| C7/clippy | PASS | `6d7a9220129c` | 26.856 | c7-results/clippy-v2/result.json |
| C7/uv-performance-outward-list | PASS | `6d7a9220129c` | 1.376 | c7-results/uv-performance-outward-v2-list/result.json |
| C7/uv-performance-outward | PASS | `6d7a9220129c` | 19.312 | c7-results/uv-performance-outward-v2/result.json |
| C7/physical-local-uv-routes-list | PASS | `6d7a9220129c` | 155.360 | c7-results/physical-local-uv-routes-v2-list/result.json |
| C7/physical-local-uv-routes | PASS | `6d7a9220129c` | 57.991 | c7-results/physical-local-uv-routes-v2/result.json |
| C8/fmt | PASS | `449a59f6f173` | 5.468 | c8-results/fmt-v2/result.json |
| C8/check | PASS | `449a59f6f173` | 0.928 | c8-results/check-v2/result.json |
| C8/build | PASS | `449a59f6f173` | 124.544 | c8-results/build-v2/result.json |
| C8/clippy | PASS | `449a59f6f173` | 0.928 | c8-results/clippy-v2/result.json |
| C8/curated-list | PASS | `449a59f6f173` | 27.086 | c8-results/curated-v5-list/result.json |
| C8/curated | PASS | `449a59f6f173` | 431.177 | c8-results/curated-v5/result.json |
| C8/scalar-all-list | PASS | `449a59f6f173` | 1240.018 | c8-results/scalar-all-v2-list/result.json |
| C8/scalar-all | PASS | `449a59f6f173` | 674.722 | c8-results/scalar-all-v2/result.json |
| C8/vertex-rules-list | PASS | `449a59f6f173` | 157.080 | c8-results/vertex-rules-v2-list/result.json |
| C8/vertex-rules | PASS | `449a59f6f173` | 1.820 | c8-results/vertex-rules-v2/result.json |
| C8/ufo-model-parity-list | PASS | `449a59f6f173` | 413.069 | c8-results/ufo-model-parity-v2-list/result.json |
| C8/ufo-model-parity | PASS | `449a59f6f173` | 1.786 | c8-results/ufo-model-parity-v2/result.json |
| C8/physical-local-uv-routes-list | PASS | `449a59f6f173` | 67.358 | c8-results/physical-local-uv-routes-list/result.json |
| C8/physical-local-uv-routes | PASS | `449a59f6f173` | 57.875 | c8-results/physical-local-uv-routes/result.json |
| C8/spenso-shadowing-check | PASS | `449a59f6f173` | 0.898 | c8-results/spenso-shadowing-check/result.json |
| C8/spenso-shadowing-tests-list | PASS | `449a59f6f173` | 1.332 | c8-results/spenso-shadowing-tests-list/result.json |
| C8/spenso-shadowing-tests | PASS | `449a59f6f173` | 3.102 | c8-results/spenso-shadowing-tests/result.json |
| C8/spenso-no-default-check | PASS | `449a59f6f173` | 0.931 | c8-results/spenso-no-default-check/result.json |
| C8/spenso-no-default-tests-list | PASS | `449a59f6f173` | 1.387 | c8-results/spenso-no-default-tests-list/result.json |
| C8/spenso-no-default-tests | PASS | `449a59f6f173` | 1.343 | c8-results/spenso-no-default-tests/result.json |
| C8/three-dimensional-reps-default-check | PASS | `449a59f6f173` | 2.643 | c8-results/three-dimensional-reps-default-check/result.json |
| C8/three-dimensional-reps-default-tests-list | PASS | `449a59f6f173` | 1.333 | c8-results/three-dimensional-reps-default-tests-list/result.json |
| C8/three-dimensional-reps-default-tests | PASS | `449a59f6f173` | 6.539 | c8-results/three-dimensional-reps-default-tests/result.json |
| C8/three-dimensional-reps-all-check | PASS | `449a59f6f173` | 0.899 | c8-results/three-dimensional-reps-all-check/result.json |
| C8/three-dimensional-reps-all-tests-list | PASS | `449a59f6f173` | 1.333 | c8-results/three-dimensional-reps-all-tests-list/result.json |
| C8/three-dimensional-reps-all-tests | PASS | `449a59f6f173` | 24.265 | c8-results/three-dimensional-reps-all-tests/result.json |
| C8/three-dimensional-reps-no-default-check | PASS | `449a59f6f173` | 0.890 | c8-results/three-dimensional-reps-no-default-check/result.json |
| C8/three-dimensional-reps-no-default-tests-list | PASS | `449a59f6f173` | 1.337 | c8-results/three-dimensional-reps-no-default-tests-list/result.json |
| C8/three-dimensional-reps-no-default-tests | PASS | `449a59f6f173` | 6.520 | c8-results/three-dimensional-reps-no-default-tests/result.json |

## Other completed or interrupted attempts

These receipts remain evidence of actual work. Their source or command differs from the final gate, or they belong to the immutable incoming baseline. They contribute no current-boundary passing count. A failed listing runs no test selection; compiler failures are not failing behavioral assertions.

| Attempt | Revision | Exit / status | Seconds | Scope or disposition |
| --- | --- | --- | --- | --- |
| c1/build | `2b713449b9ba` | 0 / PASS | 1108.954 | Different revision/arguments or superseded attempt; retained separately |
| c1/check | `9d3c3a1cfe58` | 101 / FAIL | 29.617 | E0624 in new sum fixtures: use existing public graph builder; production helper stays private and expected values stay unchanged. |
| c1/check-v2 | `2b713449b9ba` | 0 / PASS | 8.358 | Different revision/arguments or superseded attempt; retained separately |
| c1/clippy | `2b713449b9ba` | 101 / FAIL | 52.542 | Different revision/arguments or superseded attempt; retained separately |
| c1/fmt | `4dc330d044dd` | 1 / FAIL | 4.702 | Initial Linnet formatting mismatch; corrected at the owning C1 boundary. |
| c1/fmt-v2 | `9d3c3a1cfe58` | 0 / PASS | 4.047 | Different revision/arguments or superseded attempt; retained separately |
| c1/fmt-v3 | `2b713449b9ba` | 0 / PASS | 4.115 | Different revision/arguments or superseded attempt; retained separately |
| c1/tensor-foundations | `2b713449b9ba` | 0 / PASS | 11.059 | Different revision/arguments or superseded attempt; retained separately |
| c1/tensor-foundations-list | `2b713449b9ba` | 0 / PASS | 1.380 | Different revision/arguments or superseded attempt; retained separately |
| c3/build | `3ff52658fb63` | 0 / PASS | 1057.945 | Different revision/arguments or superseded attempt; retained separately |
| c3/check | `3ff52658fb63` | 0 / PASS | 32.484 | Different revision/arguments or superseded attempt; retained separately |
| c3/clippy | `3ff52658fb63` | 101 / FAIL | 37.003 | Existing SoftMomentumAssignment alias moved into C3 and nested test conditional collapsed; Clippy diagnostics, no changed expected values or runtime assertion failure. |
| c3/fmt | `3ff52658fb63` | 0 / PASS | 5.102 | Different revision/arguments or superseded attempt; retained separately |
| c4/commit4-workflows | `04637884f24f` | 100 / FAIL | 179.049 | Different revision/arguments or superseded attempt; retained separately |
| c4/commit4-workflows-list | `04637884f24f` | 0 / PASS | 2.637 | Different revision/arguments or superseded attempt; retained separately |
| c7/build | `bbb8ff1db499` | 0 / PASS | 1502.497 | Different revision/arguments or superseded attempt; retained separately |
| c7/check | `a3af927afc63` | 101 / FAIL | 30.024 | Different revision/arguments or superseded attempt; retained separately |
| c7/check-v2 | `bbb8ff1db499` | 0 / PASS | 28.112 | Different revision/arguments or superseded attempt; retained separately |
| c7/clippy | `bbb8ff1db499` | 0 / PASS | 40.699 | Different revision/arguments or superseded attempt; retained separately |
| c7/fmt | `a3af927afc63` | 0 / PASS | 5.002 | Different revision/arguments or superseded attempt; retained separately |
| c7/fmt-v2 | `bbb8ff1db499` | 0 / PASS | 6.255 | Different revision/arguments or superseded attempt; retained separately |
| c7/mre-abstract-control | `a3af927afc63` | 0 / PASS | 0.501 | Standalone MRE build/control evidence; excluded from workspace totals, see raised-energy-cff-mre-validation.md |
| c7/mre-build | `a3af927afc63` | 0 / PASS | 315.154 | Standalone MRE build/control evidence; excluded from workspace totals, see raised-energy-cff-mre-validation.md |
| c7/mre-check | `a3af927afc63` | 0 / PASS | 42.099 | Standalone MRE build/control evidence; excluded from workspace totals, see raised-energy-cff-mre-validation.md |
| c7/mre-clippy | `a3af927afc63` | 0 / PASS | 0.978 | Standalone MRE build/control evidence; excluded from workspace totals, see raised-energy-cff-mre-validation.md |
| c7/mre-control | `a3af927afc63` | 0 / PASS | 0.498 | Standalone MRE build/control evidence; excluded from workspace totals, see raised-energy-cff-mre-validation.md |
| c7/mre-fmt-current | `bbb8ff1db499` | 0 / PASS | 0.479 | Standalone MRE build/control evidence; excluded from workspace totals, see raised-energy-cff-mre-validation.md |
| c7/mre-import-read | `a3af927afc63` | 101 / FAIL | 0.493 | Standalone import diagnostic; expected assertion exit 101, separately source-scoped in raised-energy-cff-mre-validation.md |
| c7/mre-import-write | `a3af927afc63` | 0 / PASS | 0.495 | Standalone MRE build/control evidence; excluded from workspace totals, see raised-energy-cff-mre-validation.md |
| c7/physical-local-uv-routes | `bbb8ff1db499` | 0 / PASS | 56.586 | Different revision/arguments or superseded attempt; retained separately |
| c7/physical-local-uv-routes-list | `bbb8ff1db499` | 0 / PASS | 523.209 | Different revision/arguments or superseded attempt; retained separately |
| c7/uv-performance-outward | `bbb8ff1db499` | 100 / FAIL | 22.178 | Different revision/arguments or superseded attempt; retained separately |
| c7/uv-performance-outward-list | `bbb8ff1db499` | 0 / PASS | 1.848 | Different revision/arguments or superseded attempt; retained separately |
| c8/build | `1bc2316109c4` | 0 / PASS | 154.019 | Different revision/arguments or superseded attempt; retained separately |
| c8/check | `1bc2316109c4` | 0 / PASS | 1.431 | Different revision/arguments or superseded attempt; retained separately |
| c8/clippy | `1bc2316109c4` | 0 / PASS | 1.374 | Different revision/arguments or superseded attempt; retained separately |
| c8/curated-list | `449a59f6f173` | 125 / FAIL | 155.618 | Different revision/arguments or superseded attempt; retained separately |
| c8/curated-v2-list | `449a59f6f173` | 125 / FAIL | 90.879 | Different revision/arguments or superseded attempt; retained separately |
| c8/curated-v3-list | `449a59f6f173` | 101 / FAIL | 337.987 | Different revision/arguments or superseded attempt; retained separately |
| c8/curated-v4-list | `449a59f6f173` | 101 / FAIL | 1263.955 | Different revision/arguments or superseded attempt; retained separately |
| c8/fmt | `1bc2316109c4` | 0 / PASS | 5.378 | Different revision/arguments or superseded attempt; retained separately |
| c8/scalar-all | `1bc2316109c4` | 0 / PASS | 710.838 | Different revision/arguments or superseded attempt; retained separately |
| c8/scalar-all-list | `1bc2316109c4` | 0 / PASS | 1274.329 | Different revision/arguments or superseded attempt; retained separately |
| c8/ufo-model-parity | `1bc2316109c4` | 100 / FAIL | 1.874 | Different revision/arguments or superseded attempt; retained separately |
| c8/ufo-model-parity-list | `1bc2316109c4` | 0 / PASS | 1679.748 | Different revision/arguments or superseded attempt; retained separately |
| c8/vertex-rules | `1bc2316109c4` | 0 / PASS | 1.898 | Different revision/arguments or superseded attempt; retained separately |
| c8/vertex-rules-list | `1bc2316109c4` | 0 / PASS | 1.472 | Different revision/arguments or superseded attempt; retained separately |
| latest/affected-tests-list | `a1140c90c334` | 104 / FAIL | 22.150 | All-target nextest listing required unavailable gungraun-runner; corrected list/run target selection, retaining all targets for build/check/Clippy. |
| latest/affected-tests-v2 | `a1140c90c334` | 0 / PASS | 20.187 | Immutable incoming baseline |
| latest/affected-tests-v2-list | `a1140c90c334` | 0 / PASS | 2.004 | Immutable incoming baseline |
| latest/build | `a1140c90c334` | 0 / PASS | 1648.562 | Immutable incoming baseline |
| latest/check | `a1140c90c334` | 0 / PASS | 362.186 | Immutable incoming baseline |
| latest/fmt | `a1140c90c334` | 0 / PASS | 5.553 | Immutable incoming baseline |

## Execution controls and limits

The recorded commands use the licensed environment by invocation only. The report generator never reads it or executes commands from receipts. Recorded settings include a 30 GB process-tree guard, four Cargo/nextest/Rayon workers, disabled retries and snapshot updates, and `GL_NO_HARD_WARNINGS=1`; individual argument vectors and overrides remain authoritative. Clippy uses its separately recorded warning policy. Timing runs with different controls retain their own scope.

New receipts record complete pinned tracked-file Git blobs, symlink target bytes and executable modes in tracked_before/after, plus HEAD and the Spenso source fingerprint. Earlier receipts without those fields check only HEAD and crates/spenso/src/network/mod.rs; they retain that narrower scope. The report does not upgrade old guards or claim to inspect untracked files. Selected binary hashes are separate evidence.

Nextest counts require an audit matching the result fingerprint, revision, selected identities and actual outcomes. Selected binary hashes belong to their original audit time; later builds do not refresh that evidence. Changes to selection normalization are retained in the audit correction receipt. Listing and test commands are distinct executed operations.

Manual PySecDec integrations remain outside automatic selections. GL262 has a source-pinned Symbolica evaluator-construction panic and incomplete paired performance evidence. The planned physical-route selection names GL00 and GL01; it does not certify GL262. A standalone nonsymmetric function import MRE is a separate defect. No incomplete diagnostic, skipped test or guard termination is counted as passing.

Source-pinned physical performance and older benchmark/CI totals are separate evidence. The C8 source identity is an observed snapshot; copying reports and rendering PDFs requires a final documentation/closure receipt. This generator does not certify PDF compilation or visual inspection.

## Executed command appendix

Each block below is an actual receipt argument vector, rendered with shell quoting for review; it is not executed by this generator. Receipt JSON retains exact argument and environment arrays.

Tool versions: rustc: rustc 1.97.0 (2d8144b78 2026-07-07); cargo: cargo 1.97.0 (c980f4866 2026-06-30); nextest: cargo-nextest 0.9.140; typst: typst 0.15.0 (3ae52774); poppler: pdftoppm version 26.06.0.

Motivation PDF: 15 pages, inspection covers pages [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]; findings []. PDF SHA-256 `8660b1ab2481b3a9da90ff103b7a5da12bd1aaada6c20cc66138144ea186f665`, Typst SHA-256 `b7ef97198ce6b95185106a1e0cb40e5b31bbd90850dee4a6539246e6fe365dd9`. Evidence: `motivation-render/receipt.json` and `motivation-render/visual-inspection.json`. This is the parent reviewer’s recorded inspection, not a rerender by this updater.

### c1/build

Revision `2b713449b9baa63b972f1dc13cc2a973b901be6e`; PASS; receipt `c1-results/build/result.json`; SHA-256 `f2a62a230199103d0a4281cee7e01278a430fdf1063d8a89a1ff9b066108d40a`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c1-results/build/memory.jsonl --limit-gb 30 -- cargo build --workspace --all-targets --locked --profile dev-optim
```

### c1/build-v2

Revision `665658b168981c5c508e6cd51d3b57bbb55a27fb`; PASS; receipt `c1-results/build-v2/result.json`; SHA-256 `671afb7e764245f4772b854760c471a7e35b9e870423bd8b6807a8594f7ec24a`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c1-results/build-v2/memory.jsonl --limit-gb 30 -- cargo build --workspace --all-targets --locked --profile dev-optim
```

### c1/check

Revision `9d3c3a1cfe587f884f060d251e61167334c87fd0`; FAIL; receipt `c1-results/check/result.json`; SHA-256 `0a898ff638a9954cb23cd14d4e7c4db2073554cc78b45a3719ba88141887b818`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c1-results/check/memory.jsonl --limit-gb 30 -- cargo check --workspace --all-targets --locked --profile dev-optim
```

### c1/check-v2

Revision `2b713449b9baa63b972f1dc13cc2a973b901be6e`; PASS; receipt `c1-results/check-v2/result.json`; SHA-256 `f543e56b8f90fdda32a95e959a62fcea2bc7ab32a4e8692c4f8cf8e06b7fc1e7`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c1-results/check-v2/memory.jsonl --limit-gb 30 -- cargo check --workspace --all-targets --locked --profile dev-optim
```

### c1/check-v3

Revision `665658b168981c5c508e6cd51d3b57bbb55a27fb`; PASS; receipt `c1-results/check-v3/result.json`; SHA-256 `0f29ef0a36a3c5a7c838a03222023791fc08e555f87c8688749eb4d58cea7bee`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c1-results/check-v3/memory.jsonl --limit-gb 30 -- cargo check --workspace --all-targets --locked --profile dev-optim
```

### c1/clippy

Revision `2b713449b9baa63b972f1dc13cc2a973b901be6e`; FAIL; receipt `c1-results/clippy/result.json`; SHA-256 `87740eeb6e2578c051681fc677ee28ac40534dce1d91a1517b2bd969f1c74ef6`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c1-results/clippy/memory.jsonl --limit-gb 30 -- cargo clippy --workspace --all-targets --locked --profile dev-optim -- -D warnings
```

### c1/clippy-v2

Revision `665658b168981c5c508e6cd51d3b57bbb55a27fb`; PASS; receipt `c1-results/clippy-v2/result.json`; SHA-256 `392525531d99cbf77ba67de40eabdfc232eb51b0e43786943c8eaa7f706101a7`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c1-results/clippy-v2/memory.jsonl --limit-gb 30 -- cargo clippy --workspace --all-targets --locked --profile dev-optim -- -D warnings
```

### c1/fmt

Revision `4dc330d044dde11658405ae40ef3c142e3f55150`; FAIL; receipt `c1-results/fmt/result.json`; SHA-256 `900608e95be176111b8953139950d921012e5bd9364e164d11cdb661f1b91016`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c1-results/fmt/memory.jsonl --limit-gb 30 -- cargo fmt --all -- --check
```

### c1/fmt-v2

Revision `9d3c3a1cfe587f884f060d251e61167334c87fd0`; PASS; receipt `c1-results/fmt-v2/result.json`; SHA-256 `12ce391c0b6d3393622de05a58b69bc5570348f4157ebb1022fee8d4c36a7774`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c1-results/fmt-v2/memory.jsonl --limit-gb 30 -- cargo fmt --all -- --check
```

### c1/fmt-v3

Revision `2b713449b9baa63b972f1dc13cc2a973b901be6e`; PASS; receipt `c1-results/fmt-v3/result.json`; SHA-256 `4aa22194a291f0981f742523a9e8b0537edde1192009447b89abe20be5a09121`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c1-results/fmt-v3/memory.jsonl --limit-gb 30 -- cargo fmt --all -- --check
```

### c1/fmt-v4

Revision `665658b168981c5c508e6cd51d3b57bbb55a27fb`; PASS; receipt `c1-results/fmt-v4/result.json`; SHA-256 `bbe86a86ac3cb6b95450bb73a402711601d60237fa3fc76f5906d69d5e9dc573`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c1-results/fmt-v4/memory.jsonl --limit-gb 30 -- cargo fmt --all -- --check
```

### c1/spenso-no-default-check

Revision `665658b168981c5c508e6cd51d3b57bbb55a27fb`; PASS; receipt `c1-results/spenso-no-default-check/result.json`; SHA-256 `c720492f6a0195880121425564d5dc8cf62e8ec177db7008fbfa6045c1cb49a3`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c1-results/spenso-no-default-check/memory.jsonl --limit-gb 30 -- cargo check --locked --profile dev-optim -p spenso --all-targets --no-default-features
```

### c1/spenso-no-default-tests

Revision `665658b168981c5c508e6cd51d3b57bbb55a27fb`; PASS; receipt `c1-results/spenso-no-default-tests/result.json`; SHA-256 `76d808217bc1b949bcec6b156e8e22f82824b407075bb04484445c83f0fe7a1a`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c1-results/spenso-no-default-tests/memory.jsonl --limit-gb 30 -- cargo nextest run --locked --cargo-profile dev-optim --profile test_gammaloop -p spenso --lib --no-default-features -E 'all() and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c1/spenso-no-default-tests-list

Revision `665658b168981c5c508e6cd51d3b57bbb55a27fb`; PASS; receipt `c1-results/spenso-no-default-tests-list/result.json`; SHA-256 `eee587c0fe06256671058a32b49bebcc060b375b1db4e10ee4ccf14519296976`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c1-results/spenso-no-default-tests-list/memory.jsonl --limit-gb 30 -- cargo nextest list --locked --cargo-profile dev-optim --profile test_gammaloop -p spenso --lib --no-default-features -E 'all() and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --message-format json
```

### c1/spenso-shadowing-check

Revision `665658b168981c5c508e6cd51d3b57bbb55a27fb`; PASS; receipt `c1-results/spenso-shadowing-check/result.json`; SHA-256 `e07917596f54a9a870a87af7e333ef1cf3029f20c1fc9fb3eba2712303e1f489`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c1-results/spenso-shadowing-check/memory.jsonl --limit-gb 30 -- cargo check --locked --profile dev-optim -p spenso --all-targets --features shadowing
```

### c1/spenso-shadowing-tests

Revision `665658b168981c5c508e6cd51d3b57bbb55a27fb`; PASS; receipt `c1-results/spenso-shadowing-tests/result.json`; SHA-256 `e06ce1d0db1b0f8113e6e8b433507baf3154e06a88452cdbb0ba57ac5e36ac09`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c1-results/spenso-shadowing-tests/memory.jsonl --limit-gb 30 -- cargo nextest run --locked --cargo-profile dev-optim --profile test_gammaloop -p spenso --lib --features shadowing -E 'all() and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c1/spenso-shadowing-tests-list

Revision `665658b168981c5c508e6cd51d3b57bbb55a27fb`; PASS; receipt `c1-results/spenso-shadowing-tests-list/result.json`; SHA-256 `f17becb1aa4d261c137ef0d3f9412e34f7a9039871cd0748aaf9cdd3373de9e1`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c1-results/spenso-shadowing-tests-list/memory.jsonl --limit-gb 30 -- cargo nextest list --locked --cargo-profile dev-optim --profile test_gammaloop -p spenso --lib --features shadowing -E 'all() and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --message-format json
```

### c1/tensor-foundations

Revision `2b713449b9baa63b972f1dc13cc2a973b901be6e`; PASS; receipt `c1-results/tensor-foundations/result.json`; SHA-256 `d1c58c59d770122b2d3abdf8a980b724613c7fd37f0ea06298ad6f04c972e3c8`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c1-results/tensor-foundations/memory.jsonl --limit-gb 30 -- cargo nextest run --locked --cargo-profile dev-optim --profile test_gammaloop --workspace -E '(package(=spenso) or package(=idenso) or package(=linnet) or package(=spenso-hep-lib) or package(=spenso-macros)) and not test(/(^|::)(slow|failing)::/) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c1/tensor-foundations-list

Revision `2b713449b9baa63b972f1dc13cc2a973b901be6e`; PASS; receipt `c1-results/tensor-foundations-list/result.json`; SHA-256 `4e9cebc3d9c797edfb46120d627b341fdc425a3074ce1cf729c5dbada02aa3c2`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c1-results/tensor-foundations-list/memory.jsonl --limit-gb 30 -- cargo nextest list --locked --cargo-profile dev-optim --profile test_gammaloop --workspace -E '(package(=spenso) or package(=idenso) or package(=linnet) or package(=spenso-hep-lib) or package(=spenso-macros)) and not test(/(^|::)(slow|failing)::/) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --message-format json
```

### c1/tensor-foundations-v2

Revision `665658b168981c5c508e6cd51d3b57bbb55a27fb`; PASS; receipt `c1-results/tensor-foundations-v2/result.json`; SHA-256 `2f6a2e616e73436b65b091d3eb45cea7f369b2ecd49c36faa376eca9ec7265c0`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c1-results/tensor-foundations-v2/memory.jsonl --limit-gb 30 -- cargo nextest run --locked --cargo-profile dev-optim --profile test_gammaloop --workspace -E '(package(=spenso) or package(=idenso) or package(=linnet) or package(=spenso-hep-lib) or package(=spenso-macros)) and not test(/(^|::)(slow|failing)::/) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c1/tensor-foundations-v2-list

Revision `665658b168981c5c508e6cd51d3b57bbb55a27fb`; PASS; receipt `c1-results/tensor-foundations-v2-list/result.json`; SHA-256 `0d5cd06fb9d15a7bad6d5ee133abe2f41915885474b626b8301725ec52bc7d9c`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c1-results/tensor-foundations-v2-list/memory.jsonl --limit-gb 30 -- cargo nextest list --locked --cargo-profile dev-optim --profile test_gammaloop --workspace -E '(package(=spenso) or package(=idenso) or package(=linnet) or package(=spenso-hep-lib) or package(=spenso-macros)) and not test(/(^|::)(slow|failing)::/) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --message-format json
```

### c2/build

Revision `a9397ee0f86b115fdfab50e81c48dde2137692e7`; PASS; receipt `c2-results/build/result.json`; SHA-256 `dbe26cd922a386250ffaad7ab7fd974c05edfdbefacc59b84eeaa6c0bcbedc4c`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c2-results/build/memory.jsonl --limit-gb 30 -- cargo build --workspace --all-targets --locked --profile dev-optim
```

### c2/check

Revision `a9397ee0f86b115fdfab50e81c48dde2137692e7`; PASS; receipt `c2-results/check/result.json`; SHA-256 `6a7c38f1862160e52cc51d770575732700bd7bb2aca07dc01c4bc920fb5f31ed`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c2-results/check/memory.jsonl --limit-gb 30 -- cargo check --workspace --all-targets --locked --profile dev-optim
```

### c2/clippy

Revision `a9397ee0f86b115fdfab50e81c48dde2137692e7`; PASS; receipt `c2-results/clippy/result.json`; SHA-256 `2cd65b5a88622638a1b611ea056ea24a4cde69ae9b0c29e17c9b49143a80d5be`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c2-results/clippy/memory.jsonl --limit-gb 30 -- cargo clippy --workspace --all-targets --locked --profile dev-optim -- -D warnings
```

### c2/fmt

Revision `a9397ee0f86b115fdfab50e81c48dde2137692e7`; PASS; receipt `c2-results/fmt/result.json`; SHA-256 `936ea314dc8067c5626fdeb99ff02fb00d5feea01917fdc773b06477dac4db61`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c2-results/fmt/memory.jsonl --limit-gb 30 -- cargo fmt --all -- --check
```

### c2/shared-cff-default

Revision `a9397ee0f86b115fdfab50e81c48dde2137692e7`; PASS; receipt `c2-results/shared-cff-default/result.json`; SHA-256 `f821eeb6b57aa6b9a0a9325f65876a927cf9288b0b15ad4ad5907a20ff4c4678`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c2-results/shared-cff-default/memory.jsonl --limit-gb 30 -- cargo nextest run --locked --cargo-profile dev-optim --profile test_gammaloop -p three-dimensional-reps --lib -E 'package(=three-dimensional-reps) and not test(/(^|::)(slow|failing)::/) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c2/shared-cff-default-list

Revision `a9397ee0f86b115fdfab50e81c48dde2137692e7`; PASS; receipt `c2-results/shared-cff-default-list/result.json`; SHA-256 `d4b975487d17d9b23992a56cace939aff0b9c2ee818841d4b4be4808f3b20c6d`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c2-results/shared-cff-default-list/memory.jsonl --limit-gb 30 -- cargo nextest list --locked --cargo-profile dev-optim --profile test_gammaloop -p three-dimensional-reps --lib -E 'package(=three-dimensional-reps) and not test(/(^|::)(slow|failing)::/) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --message-format json
```

### c2/three-dimensional-reps-all-check

Revision `a9397ee0f86b115fdfab50e81c48dde2137692e7`; PASS; receipt `c2-results/three-dimensional-reps-all-check/result.json`; SHA-256 `ff533104458988a008e1cabbb766cedcd5c65eeddf1dc835b1b3fbd7d6c70cd3`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c2-results/three-dimensional-reps-all-check/memory.jsonl --limit-gb 30 -- cargo check --locked --profile dev-optim -p three-dimensional-reps --all-targets --all-features
```

### c2/three-dimensional-reps-all-tests

Revision `a9397ee0f86b115fdfab50e81c48dde2137692e7`; PASS; receipt `c2-results/three-dimensional-reps-all-tests/result.json`; SHA-256 `5c0fbf2bf302f5dd27229f24f0dad0ef6c22791145b23667cac908db06669c04`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c2-results/three-dimensional-reps-all-tests/memory.jsonl --limit-gb 30 -- cargo nextest run --locked --cargo-profile dev-optim --profile test_gammaloop -p three-dimensional-reps --lib --all-features -E 'all() and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c2/three-dimensional-reps-all-tests-list

Revision `a9397ee0f86b115fdfab50e81c48dde2137692e7`; PASS; receipt `c2-results/three-dimensional-reps-all-tests-list/result.json`; SHA-256 `22784f932748b05eba001a62a0b81589ce5174e413aab93acfc72f363893b5e3`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c2-results/three-dimensional-reps-all-tests-list/memory.jsonl --limit-gb 30 -- cargo nextest list --locked --cargo-profile dev-optim --profile test_gammaloop -p three-dimensional-reps --lib --all-features -E 'all() and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --message-format json
```

### c2/three-dimensional-reps-no-default-check

Revision `a9397ee0f86b115fdfab50e81c48dde2137692e7`; PASS; receipt `c2-results/three-dimensional-reps-no-default-check/result.json`; SHA-256 `7c784320455a65dbc8fc5db7bbc04a49ca0ba2840bcd753c12aed14ee2eae122`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c2-results/three-dimensional-reps-no-default-check/memory.jsonl --limit-gb 30 -- cargo check --locked --profile dev-optim -p three-dimensional-reps --all-targets --no-default-features
```

### c2/three-dimensional-reps-no-default-tests

Revision `a9397ee0f86b115fdfab50e81c48dde2137692e7`; PASS; receipt `c2-results/three-dimensional-reps-no-default-tests/result.json`; SHA-256 `e3c093ec7c55d48f3734447bfbba69b47f88caca7c4f7018271154ca52b94bdd`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c2-results/three-dimensional-reps-no-default-tests/memory.jsonl --limit-gb 30 -- cargo nextest run --locked --cargo-profile dev-optim --profile test_gammaloop -p three-dimensional-reps --lib --no-default-features -E 'all() and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c2/three-dimensional-reps-no-default-tests-list

Revision `a9397ee0f86b115fdfab50e81c48dde2137692e7`; PASS; receipt `c2-results/three-dimensional-reps-no-default-tests-list/result.json`; SHA-256 `5cf12fcf2fa08760d56ed5534eafa3332206dacfd155b4dd7fe9791b038a3392`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c2-results/three-dimensional-reps-no-default-tests-list/memory.jsonl --limit-gb 30 -- cargo nextest list --locked --cargo-profile dev-optim --profile test_gammaloop -p three-dimensional-reps --lib --no-default-features -E 'all() and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --message-format json
```

### c3/build

Revision `3ff52658fb6308ef545734b65cebb926adf61dea`; PASS; receipt `c3-results/build/result.json`; SHA-256 `ad8cc8ffdb22358e7a483d74fd6a893f55c51fcac15913540a064383e4f5eae4`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c3-results/build/memory.jsonl --limit-gb 30 -- cargo build --workspace --all-targets --locked --profile dev-optim
```

### c3/build-v2

Revision `3ea313789a1a5290093579febd7d1c601b92594e`; PASS; receipt `c3-results/build-v2/result.json`; SHA-256 `b46fda40a98e1d6325fa9dbe55951ae17a905dd37a657fe64556510607930c1e`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c3-results/build-v2/memory.jsonl --limit-gb 30 -- cargo build --workspace --all-targets --locked --profile dev-optim
```

### c3/check

Revision `3ff52658fb6308ef545734b65cebb926adf61dea`; PASS; receipt `c3-results/check/result.json`; SHA-256 `5d92c80ef9660804edee05a2a12b7c44ed12214b8c7de0758bf8423529c2a4a5`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c3-results/check/memory.jsonl --limit-gb 30 -- cargo check --workspace --all-targets --locked --profile dev-optim
```

### c3/check-v2

Revision `3ea313789a1a5290093579febd7d1c601b92594e`; PASS; receipt `c3-results/check-v2/result.json`; SHA-256 `62b82ac7c9d19b28ab924b78d587d91a374a5b2ec6c623f9325fde38b75bb192`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c3-results/check-v2/memory.jsonl --limit-gb 30 -- cargo check --workspace --all-targets --locked --profile dev-optim
```

### c3/clippy

Revision `3ff52658fb6308ef545734b65cebb926adf61dea`; FAIL; receipt `c3-results/clippy/result.json`; SHA-256 `5f978aa3b12facfedd2764cd6d30d78c887674f5516695ff9671457221d28572`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c3-results/clippy/memory.jsonl --limit-gb 30 -- cargo clippy --workspace --all-targets --locked --profile dev-optim -- -D warnings
```

### c3/clippy-v2

Revision `3ea313789a1a5290093579febd7d1c601b92594e`; PASS; receipt `c3-results/clippy-v2/result.json`; SHA-256 `ee59fad8addc68d911264ee881f50ea2a88df36e136c6f8da0daafa8e7360c66`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c3-results/clippy-v2/memory.jsonl --limit-gb 30 -- cargo clippy --workspace --all-targets --locked --profile dev-optim -- -D warnings
```

### c3/commit3-cff-uv

Revision `3ea313789a1a5290093579febd7d1c601b92594e`; PASS; receipt `c3-results/commit3-cff-uv/result.json`; SHA-256 `5b6576da7e67e03abdea98a027d8ba608d217dd1539d94314d203145e0e6f8ad`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c3-results/commit3-cff-uv/memory.jsonl --limit-gb 30 -- cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E '((package(=gammalooprs) and ((binary(=gammalooprs) and test(/^(cff::|uv::|graph::three_d_source::|numerator::energy_degree::|integrands::process::evaluators::|utils::symbols::)/)) or binary(=test_renormalization))) or (package(=gammaloop-integration-tests) and binary(=uv))) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c3/commit3-cff-uv-list

Revision `3ea313789a1a5290093579febd7d1c601b92594e`; PASS; receipt `c3-results/commit3-cff-uv-list/result.json`; SHA-256 `f9f01267cf6bced7b0d65260ec7d3ce1a778a1e0dec1a4960f6b2dee276f0b69`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c3-results/commit3-cff-uv-list/memory.jsonl --limit-gb 30 -- cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E '((package(=gammalooprs) and ((binary(=gammalooprs) and test(/^(cff::|uv::|graph::three_d_source::|numerator::energy_degree::|integrands::process::evaluators::|utils::symbols::)/)) or binary(=test_renormalization))) or (package(=gammaloop-integration-tests) and binary(=uv))) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --message-format json
```

### c3/fmt

Revision `3ff52658fb6308ef545734b65cebb926adf61dea`; PASS; receipt `c3-results/fmt/result.json`; SHA-256 `c3fd14cf4b6dff308df80e55617d865215122e17e5613e75219f4bc8c54cab31`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c3-results/fmt/memory.jsonl --limit-gb 30 -- cargo fmt --all -- --check
```

### c3/fmt-v2

Revision `3ea313789a1a5290093579febd7d1c601b92594e`; PASS; receipt `c3-results/fmt-v2/result.json`; SHA-256 `be9a6a8ecf718d4b1e0bb49ab25b61abb5bf7e2a9b7602bdf2e34cadc6892bb9`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c3-results/fmt-v2/memory.jsonl --limit-gb 30 -- cargo fmt --all -- --check
```

### c4/build

Revision `04637884f24fb28a247574f04e03308cf4c04afb`; PASS; receipt `c4-results/build/result.json`; SHA-256 `fb80395154edcce58d1b29708a6dabc6c57068832c1c36a37d5ef8b487e98407`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c4-results/build/memory.jsonl --limit-gb 30 -- cargo build --workspace --all-targets --locked --profile dev-optim
```

### c4/check

Revision `04637884f24fb28a247574f04e03308cf4c04afb`; PASS; receipt `c4-results/check/result.json`; SHA-256 `bfca021cfb09136db2229a16a3ce6087b84d99746812daee5357a7a8dbf765f1`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c4-results/check/memory.jsonl --limit-gb 30 -- cargo check --workspace --all-targets --locked --profile dev-optim
```

### c4/clippy

Revision `04637884f24fb28a247574f04e03308cf4c04afb`; PASS; receipt `c4-results/clippy/result.json`; SHA-256 `0aa56c3f4f0d8bb807d4f21f1545919243e6ad79f7be961f3a0f45710699344e`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c4-results/clippy/memory.jsonl --limit-gb 30 -- cargo clippy --workspace --all-targets --locked --profile dev-optim -- -D warnings
```

### c4/commit4-scalar-smoke

Revision `04637884f24fb28a247574f04e03308cf4c04afb`; PASS; receipt `c4-results/commit4-scalar-smoke/result.json`; SHA-256 `8bca9f4ef957dc1ca373211b1712980771d1cd6d41041475da572abb582f15fc`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c4-results/commit4-scalar-smoke/memory.jsonl --limit-gb 30 -- cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E '(binary(=test_runs) and (test(/^scalar_3l_cross_section_inspects::default_scalar_3l_cross_section_inspects::scalar_3l_cross_section_gl(00|04)_/) or test(=scalar_3l_cross_section_inspects::scalar_3l_gl00_higher_power_cff_is_invariant_under_nonzero_sampling_scale))) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c4/commit4-scalar-smoke-list

Revision `04637884f24fb28a247574f04e03308cf4c04afb`; PASS; receipt `c4-results/commit4-scalar-smoke-list/result.json`; SHA-256 `b13f807c73c10fee23ede62af8924e004fc6ad61c641f6a071bfd63f444036c0`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c4-results/commit4-scalar-smoke-list/memory.jsonl --limit-gb 30 -- cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E '(binary(=test_runs) and (test(/^scalar_3l_cross_section_inspects::default_scalar_3l_cross_section_inspects::scalar_3l_cross_section_gl(00|04)_/) or test(=scalar_3l_cross_section_inspects::scalar_3l_gl00_higher_power_cff_is_invariant_under_nonzero_sampling_scale))) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --message-format json
```

### c4/commit4-workflows

Revision `04637884f24fb28a247574f04e03308cf4c04afb`; FAIL; receipt `c4-results/commit4-workflows/result.json`; SHA-256 `4b5078dcac290ca436a7d8d1264eb44f55c311e123776873365acb9553370214`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c4-results/commit4-workflows/memory.jsonl --limit-gb 30 -- cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E '(package(=gammaloop-api) or (package(=gammalooprs) and binary(=gammalooprs) and test(/^(integrands::(evaluation|process)::|integrate::|observables::|uv::profile::)/)) or (package(=gammaloop-integration-tests) and (binary(/^(test_cli|test_evaluation_api|test_feyngen|uv)$/) or (binary(=test_differential) and test(=additional_event_weights_roundtrip_preserves_all_keys_and_values)) or (binary(=test_runs) and test(/^(inspect|events|approach|integrations|multi_integrand|profile_bulk|runtime_graph_subset|repeated_masses|test_3d_reps|smoke)::/))))) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c4/commit4-workflows-list

Revision `04637884f24fb28a247574f04e03308cf4c04afb`; PASS; receipt `c4-results/commit4-workflows-list/result.json`; SHA-256 `445a4a54e38f84ad3dd39699f061ad2035c65c5abd5de870ae1d74b76bb46af8`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c4-results/commit4-workflows-list/memory.jsonl --limit-gb 30 -- cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E '(package(=gammaloop-api) or (package(=gammalooprs) and binary(=gammalooprs) and test(/^(integrands::(evaluation|process)::|integrate::|observables::|uv::profile::)/)) or (package(=gammaloop-integration-tests) and (binary(/^(test_cli|test_evaluation_api|test_feyngen|uv)$/) or (binary(=test_differential) and test(=additional_event_weights_roundtrip_preserves_all_keys_and_values)) or (binary(=test_runs) and test(/^(inspect|events|approach|integrations|multi_integrand|profile_bulk|runtime_graph_subset|repeated_masses|test_3d_reps|smoke)::/))))) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --message-format json
```

### c4/commit4-workflows-v2

Revision `04637884f24fb28a247574f04e03308cf4c04afb`; PASS; receipt `c4-results/commit4-workflows-v2/result.json`; SHA-256 `fb3826be6a100ff70ad377ccbec978885652a95638ceb22164251c364adc08b7`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c4-results/commit4-workflows-v2/memory.jsonl --limit-gb 30 -- cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E '(package(=gammaloop-api) or (package(=gammalooprs) and binary(=gammalooprs) and test(/^(integrands::(evaluation|process)::|integrate::|observables::|uv::profile::)/)) or (package(=gammaloop-integration-tests) and (binary(/^(test_cli|test_evaluation_api|test_feyngen|uv)$/) or (binary(=test_differential) and test(=additional_event_weights_roundtrip_preserves_all_keys_and_values)) or (binary(=test_runs) and test(/^(inspect|events|approach|integrations|multi_integrand|profile_bulk|runtime_graph_subset|repeated_masses|test_3d_reps|smoke)::/))))) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c4/commit4-workflows-v2-list

Revision `04637884f24fb28a247574f04e03308cf4c04afb`; PASS; receipt `c4-results/commit4-workflows-v2-list/result.json`; SHA-256 `1602dd16f016447af50bab62da5369dfdba8eb4a894ccfbe48d88688108f9c46`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c4-results/commit4-workflows-v2-list/memory.jsonl --limit-gb 30 -- cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E '(package(=gammaloop-api) or (package(=gammalooprs) and binary(=gammalooprs) and test(/^(integrands::(evaluation|process)::|integrate::|observables::|uv::profile::)/)) or (package(=gammaloop-integration-tests) and (binary(/^(test_cli|test_evaluation_api|test_feyngen|uv)$/) or (binary(=test_differential) and test(=additional_event_weights_roundtrip_preserves_all_keys_and_values)) or (binary(=test_runs) and test(/^(inspect|events|approach|integrations|multi_integrand|profile_bulk|runtime_graph_subset|repeated_masses|test_3d_reps|smoke)::/))))) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --message-format json
```

### c4/fmt

Revision `04637884f24fb28a247574f04e03308cf4c04afb`; PASS; receipt `c4-results/fmt/result.json`; SHA-256 `c09430cc7dc6b0593eb97b1268f63e171dca6bd974efd252a3d0bedcd495e9cb`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c4-results/fmt/memory.jsonl --limit-gb 30 -- cargo fmt --all -- --check
```

### c5/build

Revision `1f2cf6d8236de91138951af8380bf1fd1eb4e783`; PASS; receipt `c5-results/build/result.json`; SHA-256 `fa6179274d28a0bd667113e9940fe9ba52e3fafefd5c8369adcae6db248893c3`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c5-results/build/memory.jsonl --limit-gb 30 -- cargo build --workspace --all-targets --locked --profile dev-optim
```

### c5/check

Revision `1f2cf6d8236de91138951af8380bf1fd1eb4e783`; PASS; receipt `c5-results/check/result.json`; SHA-256 `3062df2fe3d1f3450d86cf224194e3e0818208ba4a0dc20e3717ca5a6268ef42`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c5-results/check/memory.jsonl --limit-gb 30 -- cargo check --workspace --all-targets --locked --profile dev-optim
```

### c5/clippy

Revision `1f2cf6d8236de91138951af8380bf1fd1eb4e783`; PASS; receipt `c5-results/clippy/result.json`; SHA-256 `c04562a0cf91613e22144893cd7a58d34ea230076d617b4f30fb7585e97085c2`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c5-results/clippy/memory.jsonl --limit-gb 30 -- cargo clippy --workspace --all-targets --locked --profile dev-optim -- -D warnings
```

### c5/commit5-phases

Revision `1f2cf6d8236de91138951af8380bf1fd1eb4e783`; PASS; receipt `c5-results/commit5-phases/result.json`; SHA-256 `c1b4bb27b93805ab769c6ec1fc772470f290eace0e04af831070de8b9353d3b7`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c5-results/commit5-phases/memory.jsonl --limit-gb 30 -- cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E '((package(=gammaloop-api) and test(/^(commands::generate::|commands::import::model::)/)) or (package(=gammalooprs) and ((binary(=gammalooprs) and test(/^(cff::|feyngen::|model::|processes::cross_section::|graph::feynman_graph::|graph::lmb::|uv::export::|uv::settings::)/)) or binary(=test_renormalization))) or (package(=gammaloop-integration-tests) and (binary(=test_feyngen) or (binary(=test_runs) and test(/^(amplitude_phase_conventions|scalar_phase_conventions|scalar_virtual_phase_conventions)::/))))) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c5/commit5-phases-list

Revision `1f2cf6d8236de91138951af8380bf1fd1eb4e783`; PASS; receipt `c5-results/commit5-phases-list/result.json`; SHA-256 `8dd2902c2de22bbd85790e9f0b6a14a20132174fa16bab127cc0a174388721e9`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c5-results/commit5-phases-list/memory.jsonl --limit-gb 30 -- cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E '((package(=gammaloop-api) and test(/^(commands::generate::|commands::import::model::)/)) or (package(=gammalooprs) and ((binary(=gammalooprs) and test(/^(cff::|feyngen::|model::|processes::cross_section::|graph::feynman_graph::|graph::lmb::|uv::export::|uv::settings::)/)) or binary(=test_renormalization))) or (package(=gammaloop-integration-tests) and (binary(=test_feyngen) or (binary(=test_runs) and test(/^(amplitude_phase_conventions|scalar_phase_conventions|scalar_virtual_phase_conventions)::/))))) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --message-format json
```

### c5/fmt

Revision `1f2cf6d8236de91138951af8380bf1fd1eb4e783`; PASS; receipt `c5-results/fmt/result.json`; SHA-256 `8899de122affaecad2cffde4330ea756e3e1e0fd9c508b7242912d4df884368f`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c5-results/fmt/memory.jsonl --limit-gb 30 -- cargo fmt --all -- --check
```

### c5/signed-acceptances

Revision `1f2cf6d8236de91138951af8380bf1fd1eb4e783`; PASS; receipt `c5-results/signed-acceptances/result.json`; SHA-256 `3fa9dde3c1b7594cee3cc480a99694bcde2d78f940f76cc40e802ed863653873`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c5-results/signed-acceptances/memory.jsonl --limit-gb 30 -- cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E '(binary(/^(test_epem_a_ddx_nlo_acceptance|test_epem_a_ttx_nlo_acceptance|test_gamma_star_ddx_nlo_acceptance|test_gamma_star_ttx_nlo_acceptance)$/)) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c5/signed-acceptances-list

Revision `1f2cf6d8236de91138951af8380bf1fd1eb4e783`; PASS; receipt `c5-results/signed-acceptances-list/result.json`; SHA-256 `769afbeb6c3d5f37676d89d8b63b3c0bcdafbb17c8fffe66f81baaa598da1673`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c5-results/signed-acceptances-list/memory.jsonl --limit-gb 30 -- cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E '(binary(/^(test_epem_a_ddx_nlo_acceptance|test_epem_a_ttx_nlo_acceptance|test_gamma_star_ddx_nlo_acceptance|test_gamma_star_ttx_nlo_acceptance)$/)) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --message-format json
```

### c5/ufo-model-parity

Revision `1f2cf6d8236de91138951af8380bf1fd1eb4e783`; PASS; receipt `c5-results/ufo-model-parity/result.json`; SHA-256 `24927ea8b2da37ed3dc85d8cb733a5d0b730fb9975bb6a49939283568ba1fc29`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c5-results/ufo-model-parity/memory.jsonl --limit-gb 30 -- cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --features ufo_support -p gammaloop-api --lib -E '(package(=gammaloop-api) and binary(=gammaloop_api) and test(=commands::import::model::tests::fresh_sm_ufo_preserves_virtual_gauge_and_covariant_cut_states)) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c5/ufo-model-parity-list

Revision `1f2cf6d8236de91138951af8380bf1fd1eb4e783`; PASS; receipt `c5-results/ufo-model-parity-list/result.json`; SHA-256 `3c88a509e0c75e97e2cb82b884cb6afa14395e8c3d260b0ec336eafc5b6831c1`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c5-results/ufo-model-parity-list/memory.jsonl --limit-gb 30 -- cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --features ufo_support -p gammaloop-api --lib -E '(package(=gammaloop-api) and binary(=gammaloop_api) and test(=commands::import::model::tests::fresh_sm_ufo_preserves_virtual_gauge_and_covariant_cut_states)) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --message-format json
```

### c5/vertex-rules

Revision `1f2cf6d8236de91138951af8380bf1fd1eb4e783`; PASS; receipt `c5-results/vertex-rules/result.json`; SHA-256 `99ebb00b8d43d536595af91444da0bbd9b1565db8f1b7eae41340ffb75a2f62d`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c5-results/vertex-rules/memory.jsonl --limit-gb 30 -- cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --ignore-default-filter --run-ignored all --workspace -E '(package(=gammalooprs) and binary(=gammalooprs) and test(=graph::parse::tests::failing::vertex_rules)) and (all()) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c5/vertex-rules-list

Revision `1f2cf6d8236de91138951af8380bf1fd1eb4e783`; PASS; receipt `c5-results/vertex-rules-list/result.json`; SHA-256 `521f676d70d17e0f02dd4ebd258bb50d30e4a5ae83d861800233c2812a2c6977`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c5-results/vertex-rules-list/memory.jsonl --limit-gb 30 -- cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --ignore-default-filter --run-ignored all --workspace -E '(package(=gammalooprs) and binary(=gammalooprs) and test(=graph::parse::tests::failing::vertex_rules)) and (all()) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --message-format json
```

### c6/build

Revision `8e4e4664f3f9303760eb18d7972b91ff867b12ea`; PASS; receipt `c6-results/build/result.json`; SHA-256 `466de6371ff73c80176b90eb0753b12dd221670b0d7785966ce4a56f1d82c008`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c6-results/build/memory.jsonl --limit-gb 30 -- cargo build --workspace --all-targets --locked --profile dev-optim
```

### c6/check

Revision `8e4e4664f3f9303760eb18d7972b91ff867b12ea`; PASS; receipt `c6-results/check/result.json`; SHA-256 `4b4f2fcd0e73e9087a2bfb15dc9e59d7190a92e5a6e8ffae937f0d08f948027a`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c6-results/check/memory.jsonl --limit-gb 30 -- cargo check --workspace --all-targets --locked --profile dev-optim
```

### c6/clippy

Revision `8e4e4664f3f9303760eb18d7972b91ff867b12ea`; PASS; receipt `c6-results/clippy/result.json`; SHA-256 `3851226eba3c431812c4a41b7c64b9989bb4306c6a4c0d231047ea497878aa33`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c6-results/clippy/memory.jsonl --limit-gb 30 -- cargo clippy --workspace --all-targets --locked --profile dev-optim -- -D warnings
```

### c6/commit6-vakint-uv

Revision `8e4e4664f3f9303760eb18d7972b91ff867b12ea`; PASS; receipt `c6-results/commit6-vakint-uv/result.json`; SHA-256 `a459f2efaa8c62fd66e8f29da580aa966dda606c97c0f9501d1798d85613230d`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c6-results/commit6-vakint-uv/memory.jsonl --limit-gb 30 -- cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E '(package(=vakint) or (package(=gammalooprs) and ((binary(=gammalooprs) and test(/^uv::/)) or binary(=test_renormalization))) or (package(=gammaloop-integration-tests) and binary(=uv))) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c6/commit6-vakint-uv-list

Revision `8e4e4664f3f9303760eb18d7972b91ff867b12ea`; PASS; receipt `c6-results/commit6-vakint-uv-list/result.json`; SHA-256 `95283477ed1a54689a247d61d98cdef7e4fa5f63934e9e83c70cc9ad9a9de318`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c6-results/commit6-vakint-uv-list/memory.jsonl --limit-gb 30 -- cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --workspace -E '(package(=vakint) or (package(=gammalooprs) and ((binary(=gammalooprs) and test(/^uv::/)) or binary(=test_renormalization))) or (package(=gammaloop-integration-tests) and binary(=uv))) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --message-format json
```

### c6/fmt

Revision `8e4e4664f3f9303760eb18d7972b91ff867b12ea`; PASS; receipt `c6-results/fmt/result.json`; SHA-256 `1ea184dfcb125caa7482a9692e19db628ec06a2ba8469dd6699e033cda8ad603`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c6-results/fmt/memory.jsonl --limit-gb 30 -- cargo fmt --all -- --check
```

### c7/build

Revision `bbb8ff1db499a78b3b25c41f774c9a034acb5b4d`; PASS; receipt `c7-results/build/result.json`; SHA-256 `4233b9f24a495e385edfffc33da6a0c544f54c12f707d96deef752512d458c8c`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/build/memory.jsonl --limit-gb 30 -- cargo build --workspace --all-targets --locked --profile dev-optim
```

### c7/build-v2

Revision `6d7a9220129c0fac37bdc0e04ebefccb1688c7be`; PASS; receipt `c7-results/build-v2/result.json`; SHA-256 `428976e46f524fe2366fcda327fe614c751dc5568951e85f8cb6e22ecde49659`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/build-v2/memory.jsonl --limit-gb 30 -- cargo build --workspace --all-targets --locked --profile dev-optim
```

### c7/check

Revision `a3af927afc633c6e6c2474d3f03fa02e0242eaa4`; FAIL; receipt `c7-results/check/result.json`; SHA-256 `2159319209d267d883d5b98345e4fbcefe92559aba7df632650d3e9aeaeca319`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/check/memory.jsonl --limit-gb 30 -- cargo check --workspace --all-targets --locked --profile dev-optim
```

### c7/check-v2

Revision `bbb8ff1db499a78b3b25c41f774c9a034acb5b4d`; PASS; receipt `c7-results/check-v2/result.json`; SHA-256 `db1cd69751dc9ecfd179dd55f62ca978d4f0c509fc462d800c59f449fc352e96`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/check-v2/memory.jsonl --limit-gb 30 -- cargo check --workspace --all-targets --locked --profile dev-optim
```

### c7/check-v3

Revision `6d7a9220129c0fac37bdc0e04ebefccb1688c7be`; PASS; receipt `c7-results/check-v3/result.json`; SHA-256 `6591ba498e0fd975d9764901fe95f81ad5a029213c38793c729c59e377e44fac`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/check-v3/memory.jsonl --limit-gb 30 -- cargo check --workspace --all-targets --locked --profile dev-optim
```

### c7/clippy

Revision `bbb8ff1db499a78b3b25c41f774c9a034acb5b4d`; PASS; receipt `c7-results/clippy/result.json`; SHA-256 `bfc144b31926c3185d78f2e9e3676562f84755970cd571aef4b529051ee0fc30`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/clippy/memory.jsonl --limit-gb 30 -- cargo clippy --workspace --all-targets --locked --profile dev-optim -- -D warnings
```

### c7/clippy-v2

Revision `6d7a9220129c0fac37bdc0e04ebefccb1688c7be`; PASS; receipt `c7-results/clippy-v2/result.json`; SHA-256 `2efdde4ed94a3444639f4637f693a5cc05dd225888fd20b8146049ec6f7ea875`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/clippy-v2/memory.jsonl --limit-gb 30 -- cargo clippy --workspace --all-targets --locked --profile dev-optim -- -D warnings
```

### c7/fmt

Revision `a3af927afc633c6e6c2474d3f03fa02e0242eaa4`; PASS; receipt `c7-results/fmt/result.json`; SHA-256 `fb3c3dcbbd14f61233ff2599474a3db1660302af5f950ad488bea4ddc3a9458f`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/fmt/memory.jsonl --limit-gb 30 -- cargo fmt --all -- --check
```

### c7/fmt-v2

Revision `bbb8ff1db499a78b3b25c41f774c9a034acb5b4d`; PASS; receipt `c7-results/fmt-v2/result.json`; SHA-256 `bbb968f2ec7f0c134f4542bf1144e0186821c18a33306e2da245cee1da3f5195`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/fmt-v2/memory.jsonl --limit-gb 30 -- cargo fmt --all -- --check
```

### c7/fmt-v3

Revision `6d7a9220129c0fac37bdc0e04ebefccb1688c7be`; PASS; receipt `c7-results/fmt-v3/result.json`; SHA-256 `082ec4cdf2c02afd2625d832133f6db35f1dac4eb3a39de6ba5851952e3f88cb`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/fmt-v3/memory.jsonl --limit-gb 30 -- cargo fmt --all -- --check
```

### c7/mre-abstract-control

Revision `a3af927afc633c6e6c2474d3f03fa02e0242eaa4`; PASS; receipt `c7-results/mre-abstract-control/result.json`; SHA-256 `88e13a422634631cc1ef5abd56adbe56889fd078197d6174e9cff8b4c3ccf59c`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/mre-abstract-control/memory.jsonl --limit-gb 30 -- /common/dev/gammaloop/higher-power-energies-review/tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/target/dev-optim/symbolica-evaluator-mre --abstract-parameters --control
```

### c7/mre-build

Revision `a3af927afc633c6e6c2474d3f03fa02e0242eaa4`; PASS; receipt `c7-results/mre-build/result.json`; SHA-256 `a6333b267e0a0529ab76f1d8a72bd972bb3568072074f5f1cadc661d54415112`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/mre-build/memory.jsonl --limit-gb 30 -- cargo build --manifest-path tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/Cargo.toml --target-dir /common/dev/gammaloop/higher-power-energies-review/tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/target --locked --all-targets --profile dev-optim
```

### c7/mre-check

Revision `a3af927afc633c6e6c2474d3f03fa02e0242eaa4`; PASS; receipt `c7-results/mre-check/result.json`; SHA-256 `5cfcb30e54db391e205d0df341fdf993927246f7f1ea5e075037ad33d95c2c59`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/mre-check/memory.jsonl --limit-gb 30 -- cargo check --manifest-path tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/Cargo.toml --target-dir /common/dev/gammaloop/higher-power-energies-review/tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/target --locked --all-targets --profile dev-optim
```

### c7/mre-clippy

Revision `a3af927afc633c6e6c2474d3f03fa02e0242eaa4`; PASS; receipt `c7-results/mre-clippy/result.json`; SHA-256 `25334ae9d988694f1568dfb4640bba08b944d72c2f83ff96a5362d2fd07d68d6`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/mre-clippy/memory.jsonl --limit-gb 30 -- cargo clippy --manifest-path tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/Cargo.toml --target-dir /common/dev/gammaloop/higher-power-energies-review/tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/target --locked --all-targets --profile dev-optim -- -D warnings
```

### c7/mre-control

Revision `a3af927afc633c6e6c2474d3f03fa02e0242eaa4`; PASS; receipt `c7-results/mre-control/result.json`; SHA-256 `f48b907210faa2f50366643f318e0379f095dd660b14415d8316257e75548e2d`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/mre-control/memory.jsonl --limit-gb 30 -- /common/dev/gammaloop/higher-power-energies-review/tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/target/dev-optim/symbolica-evaluator-mre --control
```

### c7/mre-fmt-current

Revision `bbb8ff1db499a78b3b25c41f774c9a034acb5b4d`; PASS; receipt `c7-results/mre-fmt-current/result.json`; SHA-256 `792986cbc3717c7d4760b3d272e149fc2a606abd4dd6932e366fb7e71e5b4665`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/mre-fmt-current/memory.jsonl --limit-gb 30 -- cargo fmt --manifest-path tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/Cargo.toml --all -- --check
```

### c7/mre-import-read

Revision `a3af927afc633c6e6c2474d3f03fa02e0242eaa4`; FAIL; receipt `c7-results/mre-import-read/result.json`; SHA-256 `1587b2056dbd303dd98fc2eb3d6e7ee82a368dfc0df8577a8972ef7f1f207a3b`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/mre-import-read/memory.jsonl --limit-gb 30 -- /common/dev/gammaloop/higher-power-energies-review/tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/target/dev-optim/examples/import_remap read /tmp/raised-stack-latest-review-20260914/mre-tiny.symbolica
```

### c7/mre-import-write

Revision `a3af927afc633c6e6c2474d3f03fa02e0242eaa4`; PASS; receipt `c7-results/mre-import-write/result.json`; SHA-256 `d4d821be1cb21c213bfd99f89f491790ebda9390d103af0df0c5287fd532cf1e`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/mre-import-write/memory.jsonl --limit-gb 30 -- /common/dev/gammaloop/higher-power-energies-review/tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/target/dev-optim/examples/import_remap write /tmp/raised-stack-latest-review-20260914/mre-tiny.symbolica
```

### c7/physical-local-uv-routes

Revision `bbb8ff1db499a78b3b25c41f774c9a034acb5b4d`; PASS; receipt `c7-results/physical-local-uv-routes/result.json`; SHA-256 `e126be4f6119bc21d0466c25e495fff8463428642ec08f5da53582ceaf71f27c`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/physical-local-uv-routes/memory.jsonl --limit-gb 30 -- cargo nextest run --locked --cargo-profile dev-optim --profile test_gammaloop -p gammaloop-integration-tests --test uv --ignore-default-filter --run-ignored all -E 'package(=gammaloop-integration-tests) and binary(=uv) and (test(=slow::aa_aa_2l_gl00_matches_across_local_uv_routes) or test(=slow::aa_aa_2l_gl01_matches_across_local_uv_routes)) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c7/physical-local-uv-routes-list

Revision `bbb8ff1db499a78b3b25c41f774c9a034acb5b4d`; PASS; receipt `c7-results/physical-local-uv-routes-list/result.json`; SHA-256 `fe8d7b71190f69ab97a5b9656408f1c450c8b3b8eda811ef9179e6e38f4acc60`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/physical-local-uv-routes-list/memory.jsonl --limit-gb 30 -- cargo nextest list --locked --cargo-profile dev-optim --profile test_gammaloop -p gammaloop-integration-tests --test uv --ignore-default-filter --run-ignored all -E 'package(=gammaloop-integration-tests) and binary(=uv) and (test(=slow::aa_aa_2l_gl00_matches_across_local_uv_routes) or test(=slow::aa_aa_2l_gl01_matches_across_local_uv_routes)) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --message-format json
```

### c7/physical-local-uv-routes-v2

Revision `6d7a9220129c0fac37bdc0e04ebefccb1688c7be`; PASS; receipt `c7-results/physical-local-uv-routes-v2/result.json`; SHA-256 `0c954aa86237b344663d974d299c9f57f6f97a26332bdca25568459a544f0b21`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/physical-local-uv-routes-v2/memory.jsonl --limit-gb 30 -- cargo nextest run --locked --cargo-profile dev-optim --profile test_gammaloop -p gammaloop-integration-tests --test uv --ignore-default-filter --run-ignored all -E 'package(=gammaloop-integration-tests) and binary(=uv) and (test(=slow::aa_aa_2l_gl00_matches_across_local_uv_routes) or test(=slow::aa_aa_2l_gl01_matches_across_local_uv_routes)) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c7/physical-local-uv-routes-v2-list

Revision `6d7a9220129c0fac37bdc0e04ebefccb1688c7be`; PASS; receipt `c7-results/physical-local-uv-routes-v2-list/result.json`; SHA-256 `c4d5e4b2b2c8ec57cb395312a4dd0aba4a1a9249a1b5dcd327a9eb43ac33b244`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/physical-local-uv-routes-v2-list/memory.jsonl --limit-gb 30 -- cargo nextest list --locked --cargo-profile dev-optim --profile test_gammaloop -p gammaloop-integration-tests --test uv --ignore-default-filter --run-ignored all -E 'package(=gammaloop-integration-tests) and binary(=uv) and (test(=slow::aa_aa_2l_gl00_matches_across_local_uv_routes) or test(=slow::aa_aa_2l_gl01_matches_across_local_uv_routes)) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --message-format json
```

### c7/uv-performance-outward

Revision `bbb8ff1db499a78b3b25c41f774c9a034acb5b4d`; FAIL; receipt `c7-results/uv-performance-outward/result.json`; SHA-256 `6937c71c0fea22b9873948276f16368140cae643f44ac10080fe47ac62374739`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/uv-performance-outward/memory.jsonl --limit-gb 30 -- cargo nextest run --locked --cargo-profile dev-optim --profile test_gammaloop --workspace -E '(package(=gammalooprs) and binary(=gammalooprs) and test(/^(cff::|graph::three_d_source::|numerator::energy_degree::|uv::)/)) and not test(/(^|::)(slow|failing)::/) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c7/uv-performance-outward-list

Revision `bbb8ff1db499a78b3b25c41f774c9a034acb5b4d`; PASS; receipt `c7-results/uv-performance-outward-list/result.json`; SHA-256 `6ec9ccd96b9be6f909e802a4890271b50e08c553dc774f42876d0d857aa085d7`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/uv-performance-outward-list/memory.jsonl --limit-gb 30 -- cargo nextest list --locked --cargo-profile dev-optim --profile test_gammaloop --workspace -E '(package(=gammalooprs) and binary(=gammalooprs) and test(/^(cff::|graph::three_d_source::|numerator::energy_degree::|uv::)/)) and not test(/(^|::)(slow|failing)::/) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --message-format json
```

### c7/uv-performance-outward-v2

Revision `6d7a9220129c0fac37bdc0e04ebefccb1688c7be`; PASS; receipt `c7-results/uv-performance-outward-v2/result.json`; SHA-256 `d3299ba448d5f79e10f539027dbe2479d713ebefcee0c43137897bf3d92cd29b`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/uv-performance-outward-v2/memory.jsonl --limit-gb 30 -- cargo nextest run --locked --cargo-profile dev-optim --profile test_gammaloop --workspace -E '(package(=gammalooprs) and binary(=gammalooprs) and test(/^(cff::|graph::three_d_source::|numerator::energy_degree::|uv::)/)) and not test(/(^|::)(slow|failing)::/) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c7/uv-performance-outward-v2-list

Revision `6d7a9220129c0fac37bdc0e04ebefccb1688c7be`; PASS; receipt `c7-results/uv-performance-outward-v2-list/result.json`; SHA-256 `9cc48c81031dc544fc903e39789fc69774666b78dc6f2c8b1b73f7f9631aa0bb`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c7-results/uv-performance-outward-v2-list/memory.jsonl --limit-gb 30 -- cargo nextest list --locked --cargo-profile dev-optim --profile test_gammaloop --workspace -E '(package(=gammalooprs) and binary(=gammalooprs) and test(/^(cff::|graph::three_d_source::|numerator::energy_degree::|uv::)/)) and not test(/(^|::)(slow|failing)::/) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --message-format json
```

### c8/build

Revision `1bc2316109c415e9b08020fb836c614847f5d926`; PASS; receipt `c8-results/build/result.json`; SHA-256 `4ef996b20687a39d75f0aa6c51e8ccdfae61d84328bea915a3820e0460325db0`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c8-results/build/memory.jsonl --limit-gb 30 -- cargo build --workspace --all-targets --locked --profile dev-optim
```

### c8/build-v2

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/build-v2/result.json`; SHA-256 `7e5071f4f94aace139de656e45327cc4193986858610af8909e5e7dff5e6b96c`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c8-results/build-v2/memory.jsonl --limit-gb 30 -- cargo build --workspace --all-targets --locked --profile dev-optim
```

### c8/check

Revision `1bc2316109c415e9b08020fb836c614847f5d926`; PASS; receipt `c8-results/check/result.json`; SHA-256 `e5fe5a4f97e004a6eb63d09e73c6ed6ed77a4075de248ec35c4d24ac630c5d7a`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c8-results/check/memory.jsonl --limit-gb 30 -- cargo check --workspace --all-targets --locked --profile dev-optim
```

### c8/check-v2

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/check-v2/result.json`; SHA-256 `15275041a28d4a081dbbb63f4d20f2eac249fa193aabf48a278d84c9bd3c94d3`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c8-results/check-v2/memory.jsonl --limit-gb 30 -- cargo check --workspace --all-targets --locked --profile dev-optim
```

### c8/clippy

Revision `1bc2316109c415e9b08020fb836c614847f5d926`; PASS; receipt `c8-results/clippy/result.json`; SHA-256 `519fd1b8191285c0a24c87b16ea9c842209436cf478bd02879049c7bb4accd2f`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c8-results/clippy/memory.jsonl --limit-gb 30 -- cargo clippy --workspace --all-targets --locked --profile dev-optim -- -D warnings
```

### c8/clippy-v2

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/clippy-v2/result.json`; SHA-256 `34c93a374e08c39b311c1b28804ba68e85d1d98cb39958cd00b20021a4a3f6be`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c8-results/clippy-v2/memory.jsonl --limit-gb 30 -- cargo clippy --workspace --all-targets --locked --profile dev-optim -- -D warnings
```

### c8/curated-list

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; FAIL; receipt `c8-results/curated-list/result.json`; SHA-256 `fd2c79dba982ca97da13b3654b1fc382ea040951ec9a6845d032dc4e3b9a86bb`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c8-results/curated-list/memory.jsonl --limit-gb 30 -- cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked -p gammaloop-api -p gammalooprs -p idenso -p linnest -p linnet -p spenso -p spenso-hep-lib -p spenso-macros -p three-dimensional-reps -p vakint -p gammaloop-integration-tests -E '(all()) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --message-format json
```

### c8/curated-v2-list

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; FAIL; receipt `c8-results/curated-v2-list/result.json`; SHA-256 `f78145018e49c1120eeb5e858e8dbd4e250672f1ca2d219fd561272d485a38b0`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c8-results/curated-v2-list/memory.jsonl --limit-gb 30 -- cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked -p gammaloop-api -p gammalooprs -p idenso -p linnest -p linnet -p spenso -p spenso-hep-lib -p spenso-macros -p three-dimensional-reps -p vakint -p gammaloop-integration-tests -E '(all()) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --message-format json
```

### c8/curated-v3-list

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; FAIL; receipt `c8-results/curated-v3-list/result.json`; SHA-256 `263d9b970cad663fdaa92fcb0323c911e7698305d6dcae1c2b6a93820b658a99`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/curated-v3-list/memory.jsonl --limit-gb 30 -- cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked -p gammaloop-api -p gammalooprs -p idenso -p linnest -p linnet -p spenso -p spenso-hep-lib -p spenso-macros -p three-dimensional-reps -p vakint -p gammaloop-integration-tests -E '(all()) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --message-format json
```

### c8/curated-v4-list

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; FAIL; receipt `c8-results/curated-v4-list/result.json`; SHA-256 `32370ff3d68dcc29c4f1feddfc2a50d8f65d395ad4ae6cc3c47a39a4c15a44e1`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/curated-v4-list/memory.jsonl --limit-gb 30 -- cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked -p gammaloop-api -p gammalooprs -p idenso -p linnest -p linnet -p spenso -p spenso-hep-lib -p spenso-macros -p three-dimensional-reps -p vakint -p gammaloop-integration-tests -E '(all()) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --message-format json
```

### c8/curated-v5

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/curated-v5/result.json`; SHA-256 `6c2bf68a7fd060cdc2e5f0027c4ada727b6c6c642ce03120b72ee4ddad5d28c4`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/curated-v5/memory.jsonl --limit-gb 30 -- cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked -p gammaloop-api -p gammalooprs -p idenso -p linnest -p linnet -p spenso -p spenso-hep-lib -p spenso-macros -p three-dimensional-reps -p vakint -p gammaloop-integration-tests -E '(all()) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c8/curated-v5-list

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/curated-v5-list/result.json`; SHA-256 `15c19005f7b36322b949e6efc92e01608a240a39676fe1ca6f48adf504dad7f7`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/curated-v5-list/memory.jsonl --limit-gb 30 -- cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked -p gammaloop-api -p gammalooprs -p idenso -p linnest -p linnet -p spenso -p spenso-hep-lib -p spenso-macros -p three-dimensional-reps -p vakint -p gammaloop-integration-tests -E '(all()) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --message-format json
```

### c8/fmt

Revision `1bc2316109c415e9b08020fb836c614847f5d926`; PASS; receipt `c8-results/fmt/result.json`; SHA-256 `62750d326eca7f7318426fb594ce84ff1f085035bf6a11a333befffec121795a`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c8-results/fmt/memory.jsonl --limit-gb 30 -- cargo fmt --all -- --check
```

### c8/fmt-v2

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/fmt-v2/result.json`; SHA-256 `eac617d0326b00bb78e85218ea57babf15c8738a216b33f30fb90cc915bbd582`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c8-results/fmt-v2/memory.jsonl --limit-gb 30 -- cargo fmt --all -- --check
```

### c8/physical-local-uv-routes

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/physical-local-uv-routes/result.json`; SHA-256 `f0e2c937a6409619a7dace25ac88fa6f5699dc3de4558fa170829f3016810f9c`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/physical-local-uv-routes/memory.jsonl --limit-gb 30 -- cargo nextest run --locked --cargo-profile dev-optim --profile test_gammaloop -p gammaloop-integration-tests --test uv --ignore-default-filter --run-ignored all -E 'package(=gammaloop-integration-tests) and binary(=uv) and (test(=slow::aa_aa_2l_gl00_matches_across_local_uv_routes) or test(=slow::aa_aa_2l_gl01_matches_across_local_uv_routes)) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c8/physical-local-uv-routes-list

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/physical-local-uv-routes-list/result.json`; SHA-256 `82b4b27e77e07f7bedd7b3c20c28110e9416985924ecaf6f40ced083f31f7d11`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/physical-local-uv-routes-list/memory.jsonl --limit-gb 30 -- cargo nextest list --locked --cargo-profile dev-optim --profile test_gammaloop -p gammaloop-integration-tests --test uv --ignore-default-filter --run-ignored all -E 'package(=gammaloop-integration-tests) and binary(=uv) and (test(=slow::aa_aa_2l_gl00_matches_across_local_uv_routes) or test(=slow::aa_aa_2l_gl01_matches_across_local_uv_routes)) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --message-format json
```

### c8/scalar-all

Revision `1bc2316109c415e9b08020fb836c614847f5d926`; PASS; receipt `c8-results/scalar-all/result.json`; SHA-256 `74ae99f4dac4b6049f16d18bac5e948be480e32a8b4c76af07ca4ab2e1960e17`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c8-results/scalar-all/memory.jsonl --limit-gb 30 -- cargo nextest run --release --profile test_gammaloop --locked --test test_runs --ignore-default-filter --run-ignored all --workspace -E 'binary(=test_runs) and test(/^scalar_3l_cross_section_inspects::/) and not test(/(^|::)failing::/) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c8/scalar-all-list

Revision `1bc2316109c415e9b08020fb836c614847f5d926`; PASS; receipt `c8-results/scalar-all-list/result.json`; SHA-256 `6e7ce4fa3b25ebebd96de89014bcde71b1d6b08c6e562c7ae4af4e223c38c359`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c8-results/scalar-all-list/memory.jsonl --limit-gb 30 -- cargo nextest list --release --profile test_gammaloop --locked --test test_runs --ignore-default-filter --run-ignored all --workspace -E 'binary(=test_runs) and test(/^scalar_3l_cross_section_inspects::/) and not test(/(^|::)failing::/) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --message-format json
```

### c8/scalar-all-v2

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/scalar-all-v2/result.json`; SHA-256 `96ca18c328c5753848652014fc923d31092396b6c3b72d55eb2a54982d6a4595`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/scalar-all-v2/memory.jsonl --limit-gb 30 -- cargo nextest run --release --profile test_gammaloop --locked --test test_runs --ignore-default-filter --run-ignored all --workspace -E 'binary(=test_runs) and test(/^scalar_3l_cross_section_inspects::/) and not test(/(^|::)failing::/) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c8/scalar-all-v2-list

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/scalar-all-v2-list/result.json`; SHA-256 `9693d9078d45ac0595783ecde95b660a728f8c580730a5992a177af8b187d83f`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/scalar-all-v2-list/memory.jsonl --limit-gb 30 -- cargo nextest list --release --profile test_gammaloop --locked --test test_runs --ignore-default-filter --run-ignored all --workspace -E 'binary(=test_runs) and test(/^scalar_3l_cross_section_inspects::/) and not test(/(^|::)failing::/) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --message-format json
```

### c8/spenso-no-default-check

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/spenso-no-default-check/result.json`; SHA-256 `45cf85717878483df8ff230967bbadc5f820eb1b1cbaba0583835af1fa5aaf85`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/spenso-no-default-check/memory.jsonl --limit-gb 30 -- cargo check --locked --profile dev-optim -p spenso --all-targets --no-default-features
```

### c8/spenso-no-default-tests

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/spenso-no-default-tests/result.json`; SHA-256 `16503506dc7aa5d4308de9cc3e6433301b7ca83d4ed6df122384bd7a2251ad9a`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/spenso-no-default-tests/memory.jsonl --limit-gb 30 -- cargo nextest run --locked --cargo-profile dev-optim --profile test_gammaloop -p spenso --lib --no-default-features -E 'all() and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c8/spenso-no-default-tests-list

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/spenso-no-default-tests-list/result.json`; SHA-256 `a2d7f663bf576b6b23c79c1a3de15e4a4c8c8a6cc0fd24a83dbc95ba5ab73613`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/spenso-no-default-tests-list/memory.jsonl --limit-gb 30 -- cargo nextest list --locked --cargo-profile dev-optim --profile test_gammaloop -p spenso --lib --no-default-features -E 'all() and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --message-format json
```

### c8/spenso-shadowing-check

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/spenso-shadowing-check/result.json`; SHA-256 `9eb3d53767cbe372102831d6c4059a2086c3bc07d76ac1ebc7cf130062796515`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/spenso-shadowing-check/memory.jsonl --limit-gb 30 -- cargo check --locked --profile dev-optim -p spenso --all-targets --features shadowing
```

### c8/spenso-shadowing-tests

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/spenso-shadowing-tests/result.json`; SHA-256 `f43618e38d511cf700cf90d8b4956c29fa5354135e36c49231f61bdb6e465309`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/spenso-shadowing-tests/memory.jsonl --limit-gb 30 -- cargo nextest run --locked --cargo-profile dev-optim --profile test_gammaloop -p spenso --lib --features shadowing -E 'all() and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c8/spenso-shadowing-tests-list

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/spenso-shadowing-tests-list/result.json`; SHA-256 `5c0b1576821049846883001f92319c19da0ed13554e43f41b6a3b922829e24eb`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/spenso-shadowing-tests-list/memory.jsonl --limit-gb 30 -- cargo nextest list --locked --cargo-profile dev-optim --profile test_gammaloop -p spenso --lib --features shadowing -E 'all() and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --message-format json
```

### c8/three-dimensional-reps-all-check

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/three-dimensional-reps-all-check/result.json`; SHA-256 `207bfe31e1c947923ccff668f29ed4f7411b0268b6938dcf5260496e0e3629b9`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/three-dimensional-reps-all-check/memory.jsonl --limit-gb 30 -- cargo check --locked --profile dev-optim -p three-dimensional-reps --all-targets --all-features
```

### c8/three-dimensional-reps-all-tests

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/three-dimensional-reps-all-tests/result.json`; SHA-256 `9fde9cf82ee98681507aa5db057cd0e121797518e3ff3b8e5da9f49667cb2bac`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/three-dimensional-reps-all-tests/memory.jsonl --limit-gb 30 -- cargo nextest run --locked --cargo-profile dev-optim --profile test_gammaloop -p three-dimensional-reps --lib --all-features -E 'all() and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c8/three-dimensional-reps-all-tests-list

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/three-dimensional-reps-all-tests-list/result.json`; SHA-256 `ed6e9307f7adca02794e5da0bdac7c34e7a56bc943e0444fbfc1a7adb83d0bc9`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/three-dimensional-reps-all-tests-list/memory.jsonl --limit-gb 30 -- cargo nextest list --locked --cargo-profile dev-optim --profile test_gammaloop -p three-dimensional-reps --lib --all-features -E 'all() and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --message-format json
```

### c8/three-dimensional-reps-default-check

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/three-dimensional-reps-default-check/result.json`; SHA-256 `a27faede0a38f413b3ca32e11b199adb65652623dbbe217db19ca0afbed1d45d`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/three-dimensional-reps-default-check/memory.jsonl --limit-gb 30 -- cargo check --locked --profile dev-optim -p three-dimensional-reps --all-targets
```

### c8/three-dimensional-reps-default-tests

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/three-dimensional-reps-default-tests/result.json`; SHA-256 `4b38e98888c4b03e6ffc1779294e1702a8aec0ba1f0f2c353b819c3fbaa15378`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/three-dimensional-reps-default-tests/memory.jsonl --limit-gb 30 -- cargo nextest run --locked --cargo-profile dev-optim --profile test_gammaloop -p three-dimensional-reps --lib -E 'all() and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c8/three-dimensional-reps-default-tests-list

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/three-dimensional-reps-default-tests-list/result.json`; SHA-256 `786b07d5319886149f985fc9026f7b3fa5161816c9f4cc1079f0ae398634fe4a`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/three-dimensional-reps-default-tests-list/memory.jsonl --limit-gb 30 -- cargo nextest list --locked --cargo-profile dev-optim --profile test_gammaloop -p three-dimensional-reps --lib -E 'all() and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --message-format json
```

### c8/three-dimensional-reps-no-default-check

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/three-dimensional-reps-no-default-check/result.json`; SHA-256 `27f88eb43dbead20228fce0fa6332645fe5f9c616ab69db8be78591b572bc9b3`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/three-dimensional-reps-no-default-check/memory.jsonl --limit-gb 30 -- cargo check --locked --profile dev-optim -p three-dimensional-reps --all-targets --no-default-features
```

### c8/three-dimensional-reps-no-default-tests

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/three-dimensional-reps-no-default-tests/result.json`; SHA-256 `9b9b90cdfa80d0e7013ab49e306790d3194e233f5204d047b00eaa220f85a6e0`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/three-dimensional-reps-no-default-tests/memory.jsonl --limit-gb 30 -- cargo nextest run --locked --cargo-profile dev-optim --profile test_gammaloop -p three-dimensional-reps --lib --no-default-features -E 'all() and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c8/three-dimensional-reps-no-default-tests-list

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/three-dimensional-reps-no-default-tests-list/result.json`; SHA-256 `049adca9c7d7dbf934aad689124e75f428627d6289937b2bbed55431e2f39a9c`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/three-dimensional-reps-no-default-tests-list/memory.jsonl --limit-gb 30 -- cargo nextest list --locked --cargo-profile dev-optim --profile test_gammaloop -p three-dimensional-reps --lib --no-default-features -E 'all() and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))' --message-format json
```

### c8/ufo-model-parity

Revision `1bc2316109c415e9b08020fb836c614847f5d926`; FAIL; receipt `c8-results/ufo-model-parity/result.json`; SHA-256 `4642c404b270b3e58c140f0e620954352b8cf3a5c624d699bd49a82e74f48421`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c8-results/ufo-model-parity/memory.jsonl --limit-gb 30 -- cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --features gammaloop-api/ufo_support --workspace -E '(package(=gammaloop-api) and binary(=gammaloop_api) and test(=commands::import::model::tests::fresh_sm_ufo_preserves_virtual_gauge_and_covariant_cut_states)) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c8/ufo-model-parity-list

Revision `1bc2316109c415e9b08020fb836c614847f5d926`; PASS; receipt `c8-results/ufo-model-parity-list/result.json`; SHA-256 `85b9bcff7f4e4c345d3f32bb88e2d7b4461ba2925c338c80edc6e48d9f9323f3`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c8-results/ufo-model-parity-list/memory.jsonl --limit-gb 30 -- cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --features gammaloop-api/ufo_support --workspace -E '(package(=gammaloop-api) and binary(=gammaloop_api) and test(=commands::import::model::tests::fresh_sm_ufo_preserves_virtual_gauge_and_covariant_cut_states)) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --message-format json
```

### c8/ufo-model-parity-v2

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/ufo-model-parity-v2/result.json`; SHA-256 `fac3c3ac335b2aa77afa31080818b60e843e8d84c8a1529e3bf46af4dc5ea365`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/ufo-model-parity-v2/memory.jsonl --limit-gb 30 -- cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --features ufo_support -p gammaloop-api --lib -E '(package(=gammaloop-api) and binary(=gammaloop_api) and test(=commands::import::model::tests::fresh_sm_ufo_preserves_virtual_gauge_and_covariant_cut_states)) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c8/ufo-model-parity-v2-list

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/ufo-model-parity-v2-list/result.json`; SHA-256 `36e9683e6d1ca72ab6a6933516f1be55e23dfdb69ec3c2ee20c566c2d8ab43d0`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/ufo-model-parity-v2-list/memory.jsonl --limit-gb 30 -- cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --features ufo_support -p gammaloop-api --lib -E '(package(=gammaloop-api) and binary(=gammaloop_api) and test(=commands::import::model::tests::fresh_sm_ufo_preserves_virtual_gauge_and_covariant_cut_states)) and (not test(/(^|::)(slow|failing)::/)) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --message-format json
```

### c8/vertex-rules

Revision `1bc2316109c415e9b08020fb836c614847f5d926`; PASS; receipt `c8-results/vertex-rules/result.json`; SHA-256 `a8e9425299616d97a62f91b8837f8435966c6ed24488f6cc29387a4f4a5d0dc7`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c8-results/vertex-rules/memory.jsonl --limit-gb 30 -- cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --ignore-default-filter --run-ignored all --workspace -E '(package(=gammalooprs) and binary(=gammalooprs) and test(=graph::parse::tests::failing::vertex_rules)) and (all()) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c8/vertex-rules-list

Revision `1bc2316109c415e9b08020fb836c614847f5d926`; PASS; receipt `c8-results/vertex-rules-list/result.json`; SHA-256 `355bba119eda69574daa36aff37239bbae1c711d3d14903a58d45ec8f406b558`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/c8-results/vertex-rules-list/memory.jsonl --limit-gb 30 -- cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --ignore-default-filter --run-ignored all --workspace -E '(package(=gammalooprs) and binary(=gammalooprs) and test(=graph::parse::tests::failing::vertex_rules)) and (all()) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --message-format json
```

### c8/vertex-rules-v2

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/vertex-rules-v2/result.json`; SHA-256 `2bcf88e994909ccf1bb8627671e9cb4309a972c5e0a963fef42199788ce590fd`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/vertex-rules-v2/memory.jsonl --limit-gb 30 -- cargo nextest run --cargo-profile dev-optim --profile test_gammaloop --locked --ignore-default-filter --run-ignored all --workspace -E '(package(=gammalooprs) and binary(=gammalooprs) and test(=graph::parse::tests::failing::vertex_rules)) and (all()) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### c8/vertex-rules-v2-list

Revision `449a59f6f173d100c2da858c541faa1e0dccca6e`; PASS; receipt `c8-results/vertex-rules-v2-list/result.json`; SHA-256 `fbe47cc3b4e5a4b6e6e6679097615f21c7b5a878e071da4eafb2217b0cac2839`.

```sh
/tmp/raised-stack/run env PYTHONPATH=/tmp/raised-stack/python-deps /tmp/raised-stack/python /tmp/raised-stack-latest-review-20260914/ram-watchdog-ps10.py --log /tmp/raised-stack-latest-review-20260914/c8-results/vertex-rules-v2-list/memory.jsonl --limit-gb 30 -- cargo nextest list --cargo-profile dev-optim --profile test_gammaloop --locked --ignore-default-filter --run-ignored all --workspace -E '(package(=gammalooprs) and binary(=gammalooprs) and test(=graph::parse::tests::failing::vertex_rules)) and (all()) and (not test(/^aa_aa::important::aa_aa_local_inspect_backend_consistency$/)) and (not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec)))' --message-format json
```

### latest/affected-tests-list

Revision `a1140c90c334ff58a3b040ff3eae06d3645bf0bb`; FAIL; receipt `latest-results/affected-tests-list/result.json`; SHA-256 `25a32743c646a9fdd3237e57e5fd17e36c56b6556477b1cc854e084abaa60b2f`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/latest-results/affected-tests-list/memory.jsonl --limit-gb 30 -- cargo nextest list --workspace --all-targets --locked --cargo-profile dev-optim --profile test_gammaloop -E '((package(=spenso) or package(=idenso) or package(=linnet) or package(=spenso-hep-lib) or package(=spenso-macros)) and not test(/(^|::)(slow|failing)::/) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))) or ((package(=gammalooprs) and binary(=gammalooprs) and test(/^(cff::|graph::three_d_source::|numerator::energy_degree::|uv::)/)) and not test(/(^|::)(slow|failing)::/) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))) or (package(=gammalooprs) and binary(=gammalooprs) and test(/^integrands::process::evaluators::/))' --message-format json
```

### latest/affected-tests-v2

Revision `a1140c90c334ff58a3b040ff3eae06d3645bf0bb`; PASS; receipt `latest-results/affected-tests-v2/result.json`; SHA-256 `c3254bfcf9636cb1cfb193388072c9cc6d9c90bf3c64d51907ce7362affa8b5e`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/latest-results/affected-tests-v2/memory.jsonl --limit-gb 30 -- cargo nextest run --workspace --locked --cargo-profile dev-optim --profile test_gammaloop -E '((package(=spenso) or package(=idenso) or package(=linnet) or package(=spenso-hep-lib) or package(=spenso-macros)) and not test(/(^|::)(slow|failing)::/) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))) or ((package(=gammalooprs) and binary(=gammalooprs) and test(/^(cff::|graph::three_d_source::|numerator::energy_degree::|uv::)/)) and not test(/(^|::)(slow|failing)::/) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))) or (package(=gammalooprs) and binary(=gammalooprs) and test(/^integrands::process::evaluators::/))' --test-threads 4 --retries 0 --no-fail-fast --status-level pass --final-status-level all --failure-output immediate-final --success-output never --no-tests fail
```

### latest/affected-tests-v2-list

Revision `a1140c90c334ff58a3b040ff3eae06d3645bf0bb`; PASS; receipt `latest-results/affected-tests-v2-list/result.json`; SHA-256 `4bc8d0e29c07ca9d47df11a9a36afb0152c4fbb795970adc5dfe18b0d0c350f1`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/latest-results/affected-tests-v2-list/memory.jsonl --limit-gb 30 -- cargo nextest list --workspace --locked --cargo-profile dev-optim --profile test_gammaloop -E '((package(=spenso) or package(=idenso) or package(=linnet) or package(=spenso-hep-lib) or package(=spenso-macros)) and not test(/(^|::)(slow|failing)::/) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))) or ((package(=gammalooprs) and binary(=gammalooprs) and test(/^(cff::|graph::three_d_source::|numerator::energy_degree::|uv::)/)) and not test(/(^|::)(slow|failing)::/) and not (binary(/pysecdec/) or test(=test_integrate_1l_decorated_indices_pysecdec))) or (package(=gammalooprs) and binary(=gammalooprs) and test(/^integrands::process::evaluators::/))' --message-format json
```

### latest/build

Revision `a1140c90c334ff58a3b040ff3eae06d3645bf0bb`; PASS; receipt `latest-results/build/result.json`; SHA-256 `b162d334fd46624e8e7ab738dc78ae8c2cabac9168d9f4a4f08586c9bbe3d950`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/latest-results/build/memory.jsonl --limit-gb 30 -- cargo build --workspace --all-targets --locked --profile dev-optim
```

### latest/check

Revision `a1140c90c334ff58a3b040ff3eae06d3645bf0bb`; PASS; receipt `latest-results/check/result.json`; SHA-256 `1cefbed33a375c559fb8f67599b63302349d5f299bf5e19e7ec97eca0b5af089`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/latest-results/check/memory.jsonl --limit-gb 30 -- cargo check --workspace --all-targets --locked --profile dev-optim
```

### latest/fmt

Revision `a1140c90c334ff58a3b040ff3eae06d3645bf0bb`; PASS; receipt `latest-results/fmt/result.json`; SHA-256 `cc9a27cc401b746eb1951b48dc037ac3d943726200c9de44b3f9515c38d22a1b`.

```sh
/tmp/raised-stack/run /tmp/raised-stack/python /common/dev/gammaloop/higher-power-energies-review/bin/ram_watchdog.py --log /tmp/raised-stack-latest-review-20260914/latest-results/fmt/memory.jsonl --limit-gb 30 -- cargo fmt --all -- --check
```
