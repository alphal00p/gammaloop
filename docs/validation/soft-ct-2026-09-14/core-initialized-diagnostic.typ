= Initialized-process diagnostic for the isolated logarithmic test

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.
<initialized-process-diagnostic-for-the-isolated-logarithmic-test>
Exactly two unchanged tests from the previously built GammaLoop library test binary were selected with two full-name `--exact` filters and run together using `--test-threads 1`. Libtest help documents multiple filters as an OR and default alphabetical ordering; a preliminary `--list` returned exactly these two tests in the observed execution order. No rebuild, source/test edits, license extraction/export, or alternate activation was used. The existing campaign environment was sourced unchanged.

+ `uv::approx::local_4d::tests::local_soft_provenance_records_the_selected_route_and_branch_sizes` --- passed previously in the isolated core campaign (0.062 s); calls normal `test_initialise().unwrap()` first at #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:4008")[local\_4d.rs:4008];.
+ `uv::approx::local_4d::tests::logarithmic_tilde_is_zero_so_hat_is_t` --- its existing zero assertion is unchanged at #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3998")[local\_4d.rs:3998];. This test does not call the initialization helper itself.

#strong[Observed result: 2 passed, 0 failed.] The runner measured 0.061044 seconds wall time, 0.048621 seconds user CPU, 0.016678 seconds system CPU, and 27,820 KiB largest-child peak RSS (about 27.17 MiB). This is combined process cost, not an independently measured duration for the second test. Libtest reported its own total as 0.05 s.

The helper uses the application\'s existing `initialise()` and `init_vakint()` calls: #link("../../../crates/gammalooprs/src/initialisation.rs:41")[initialisation.rs:41];. The result establishes that the exact original logarithmic assertion passes in the normal initialized process. Combined with its fresh-process startup abort and absent per-test initialization, this supports a missing bootstrap boundary in isolated execution, rather than a reproduced functional failure of that assertion. It does not change the original Nextest result and does not address the five other core assertion failures or the seven separate support-binary startup failures.

Evidence: unchanged-test execution log;, runner command and resources;. The diagnostic is a separate two-test execution with one deliberately reused initializer; it contributes no new unique test to the original campaign denominator.
