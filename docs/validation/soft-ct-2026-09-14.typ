= Soft-counterterm rebase: validation and performance, 14 September 2026

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.
<soft-counterterm-rebase-validation-and-performance-14-september-2026>
Follow-up: the #link("soft-ct-expression-growth-2026-09-14.typ")[matched expression-growth investigation] reproduces large scalar materialization with ordinary UV alone and documents the zero central-child contribution in the original selected orientation. The 314 GiB failure must not be interpreted as an established intrinsic soft-counterterm cost.

The jj rebase is mechanically complete: the current stack has no conflicted ancestors. Validation is #strong[not fully green];. The supported soft-counterterm routes have substantial passing symbolic, structural and numerical coverage, including the measured two-power top-bubble soft improvement. Supplementary Appendix B1 scans expose a projected-route numerical UV failure, in addition to assertion failures and Symbolica startup issues. The large full-orientation Figure B1 generator panicked after 105 minutes with a 313.962 GiB peak resident-memory cost. A complete acceptance sign-off is not justified.

This report describes the current revision and observed results. A failing Rust test remains a failure even when an independent symbolic certificate explains its residual.

== Execution outcome
<execution-outcome>
#figure(
  align(center)[#table(
    columns: 7,
    align: (auto,right,right,right,right,right,right,),
    table.header([Group], [Selected], [Pass], [Test failure], [Timeout], [Startup abort], [JUnit wall sum],),
    table.hline(),
    [Core units], [122], [116], [5], [0], [1], [22.570 s],
    [Supporting crates], [7], [0], [0], [0], [7], [7.055 s],
    [Regular integration], [17], [17], [0], [0], [0], [130.695 s],
    [Mandatory acceptance], [6], [4], [2], [0], [0], [6,733.246 s],
  )]
  , kind: table
  )

#strong[152/152 unique original test executions completed: 137 passed, 7 failed (6 assertion failures and 1 generation panic), 0 timed out, and 8 aborted at startup.] Retries and supplementary diagnostics are reported separately, without double-counting cases. JUnit wall sums add serial testcase walls within a group; they exclude build/setup cost and are not total task elapsed time.

== Revision and changes
<revision-and-changes>
Production input: `efaf2da158c5b01bddbc5435592e7cfea7a45a91`, jj change `wxtyqpvsvkxssrlwmymsuquowsyoyxtq`, rebased onto `faster_uv_from_4d@origin` at `c08b0d9f1605f7e3ec7445508c69edc07a4b599f`. The validated source snapshot after two user-approved API repairs is `9cfaac7c4a380ad9fa1cb904b8930257950609cc`, in new jj change `ptpvxxoxsyztwnypstosvlnlnwlunrkq`. Adding this report changes the working-copy commit ID but not the tested Rust source.

The only source edits made for this validation are:

- #link("../../crates/gammalooprs/src/uv/approx/direct_3d/soft_tests.rs:173")[soft\_tests.rs:173];: pass a one-element array to the rebased `zip_add` API.
- #link("../../tests/tests/uv.rs:2228")[uv.rs:2228];: call the existing `evaluate_momentum_sample` with `orientation=None` to request the summed evaluation.

Both edits were explicitly approved. Assertions and production code are unchanged. Initial compilation failures from these obsolete calls are retained in the command evidence; all listed test executions use the repaired source.

== Method and scope
<method-and-scope>
The campaign selected 122 core unit tests, seven supporting tensor/backend regressions, 17 regular integration tests and the six mandatory soft-CT acceptance tests defined in #link("../../.config/nextest.toml:121")[.config/nextest.toml];. The regular selection includes the 14 dedicated soft-CT cases plus scalar route agreement, integrated spectacles route agreement and MUV/PolePart sunrise controls. It is a focused soft-CT/reconstruction campaign, not the entire workspace suite.

Every nextest run used `--locked --offline --cargo-profile dev-optim --retries 0 --no-fail-fast --no-capture`. `--no-capture` forces serial execution; nextest consequently warns that the additional `--test-threads 1` flag is redundant. Regular tests used `ci_gammaloop`; mandatory tests used `soft_ct_acceptance`. The latter has a 60-minute slow period and termination after two periods. The `local_test` profile was deliberately avoided because it can mark timeouts as passes. No timeout was counted as success.

Three local routes are compared: #strong[localized 3D] (`explicit=false, local4D=false`), #strong[explicit 3D] (`true,false`) and #strong[projected 4D → 3D] (`true,true`). Projection and explicit orientation summation are separate settings. Local asymptotic acceptance disables integrated UV addbacks and threshold subtraction. Separate integrated tests enable addbacks and compare complete arbitrary-precision point evaluations. “Integrated” here does not mean a Monte Carlo integral or a physical cross section was evaluated.

Environment: NixOS 26.11, Linux 6.18.45, AMD EPYC 9754, 384 available logical CPUs, approximately 1.1 TiB host RAM; Rust/Cargo 1.97.0, nextest 0.9.140, Python 3.14.6 and FORM 5. Build jobs were capped at 16, Rayon at eight, and the default Rust thread stack at 128 MiB (`RUST_MIN_STACK`); some tests explicitly use 64 MiB. The cards report one generation core. `dev-optim` uses optimization level 2 with debug information and dependency-specific overrides. Cached Nix tools supplied the build environment after the offline development-shell attempt could not realize missing dependencies.

This was a shared host with high, changing load. There is no matched pre-rebase benchmark. Timings describe this run and cannot establish a rebase speedup or regression. Compilation, test execution, graph/evaluator construction and repeated sample evaluation are reported separately.

Each run used a fresh directory under `/tmp/soft-ct-validation-2026-09-14`, with successful states retained. Existing generation/profile logging recorded summaries without distributing graph numerators. Test-created `global_settings.toml` files often precede in-memory route overrides; replay must restore those overrides explicitly. Command metadata records exact argv, timestamps, CPU time, command wall time and largest-child peak RSS. JUnit records the actual per-test wall time.

== What the test cases establish
<what-the-test-cases-establish>
#figure(
  align(center)[#table(
    columns: 3,
    align: (auto,right,auto,),
    table.header([Group], [Cases], [Tested behavior],),
    table.hline(),
    [Direct 3D soft operator], [16], [Exact `H=U+S−US` at degrees 0/1/2; factorized `S+U(X−S)`; forest signs; one signed completed atom; boundary scaling; nested masses and fixed cograph poles.],
    [Local 4D operator], [56], [Taylor/soft jets and paper formulas; physical versus UV mass grading; nested compatible charts; completed-child replay; factorized momentum provenance; integrated massless, massive and mixed-mass cases.],
    [Analytic integration], [17], [Dirac reduction before epsilon truncation, retained child Laurent coefficients, squared masses, loop incidence, tensor closure and backend routing.],
    [Scheme orchestration], [7], [Integrated IR acceptance and contextual rejection of unsupported OS/mixed-scheme woods.],
    [Forest construction], [22], [Dependency order, disconnected products/unions, independently divergent components and ordinary graph-topology baselines.],
    [Projected 4D reconstruction], [4], [Exact powered-component sign, retention/cache invariants, batched versus sequential projection and reuse across residue states.],
    [Spenso/Idenso], [4], [Closed lazy scalar sums, tensor powers, trace-projector slots and the GammaLoop parser fixture; execution blocked before assertions.],
    [Vakint], [3], [Partly massless sunset matching and independent gamma-integral oracles, plus first-use dot attributes; execution blocked before assertions.],
  )]
  , kind: table
  )

The #link("soft-ct-2026-09-14/unit-case-purposes.typ")[unit case inventory] explains individual cases. The #link("soft-ct-2026-09-14/integration-case-purposes.typ")[integration inventory] explains the fixtures, expected forest contents and assertions. The #link("soft-ct-2026-09-14/test-cases.typ")[execution appendix] lists every selected case, result and elapsed time.

The integration matrix deliberately contains both positive and negative controls. Bare profiles must fail; subtracted profiles must pass on matching rays. Disconnected and consecutive scalar fixtures require nonzero component and union atoms. Persistence starts fresh subprocesses to save/reload/regenerate generated soft CFF state. Vacuum tests compare six variants to a nonzero analytic reference, `Im = −9.760529078735244e−4`, with absolute tolerance `1e−12`. Thirteen integrated-IR variants cover quark, gluon, spectacles and nested top fixtures across the local routes, including disabled integrated tensor projection for the nested explicit route. Complete arbitrary-precision point sums agree within 1,000 times reported numerical accuracy, with the source\'s serialization floor; spectacles and nested results must be nonzero.

== UV and soft scaling results
<uv-and-soft-scaling-results>
UV slopes are degrees of divergence fitted after multiplying the sampled integrand by `λ^(3 n_hard)`, where `n_hard` is the number of scaled loop momenta (#link("../../crates/gammalooprs/src/uv/profile.rs:2854")[implementation];). A reported one-hard slope of +2 therefore corresponds to a raw integrand exponent near −1. Standard UV scans use 25 points from `10^8` to `10^12`, seed 1337; DGSE uses 33 points through `10^16`. Passing fits require slope at most −0.9 and finite `R² ≥ 0.99`. Resolved fits, certified vanishing rays and missing fits are counted separately in the numerical evidence. An expected failing bare control is not a test failure.

All #strong[1,217 regular subtracted/ordinary-MUV UV records pass];: 228 summed and 989 orientation records. All have fitted slopes; none relies on a vanishing-fit exemption. In particular:

- DGSE resolves all 27 summed rays in each route, and all 810 localized orientation/ray combinations across 30 orientations. The least-negative summed slopes are −0.999811725, −0.999811792 and −0.998761370 for localized, explicit and projected routes. The matching bare control fails 17/27 summed and 494/810 orientation records.
- The generation-filtered DGSE orientation resolves and passes 27/27 individual rays, with value agreement to the same orientation in the unfiltered integrand.
- Massless quark slopes range from −1.027403362 to −0.935352940, with minimum `R²=0.992963692`; gluon slopes are approximately −1.00001. All three routes pass.
- Consecutive and disconnected scalar components each pass all 12 summed rays per route and 48 localized orientation/ray combinations. Summed slopes are −2 for single-hard and −4 for double-hard limits; the corresponding bare controls grow as +2/+4.

The regular DGSE soft test also passes its fixture-specific bounds: same routed ray, finite `R²≥0.9`, subtracted scaling greater than 3, bare scaling at most 3 and improvement at least 1. Its exact soft fit values were not persisted because the test requests no soft output file. These DGSE bounds use the same measure-adjusted soft-scaling definition as the top-bubble values below, but impose different acceptance requirements. Full regular measurements and control failure categories are in #link("soft-ct-2026-09-14/regular-numerical.typ")[regular numerical results];.

For infrared profiles, reported `scaling = raw fitted exponent + 3 n_soft` (#link("../../crates/gammalooprs/src/integrands/process/ir.rs:700")[implementation];). For the top bubble, the same `S(e6)` ray uses 25 points from `10^-2` to `10^-5`, seed 1337, canonical LMB `[e6,e8]` and fingerprint `94ac2449745d853c`. Here `n_soft=1`. All three routes give the same printed summed scaling values:

#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,right,right,right,),
    table.header([Fixture], [Ordinary U scaling], [Soft H scaling], [H−U improvement],),
    table.hline(),
    [Child-only top bubble], [−2.001104448869], [−0.000573378786], [2.000531070083],
    [Complete top-bubble forest], [−2.001104448869], [−0.000572806050], [2.000531642819],
  )]
  , kind: table
  )

All summed fits have `R²>0.99999996`. The full bare slope is −2.001114489716. Localized mode checks all 30 soft orientations: exactly eight ordinary-U orientations have slope below −1, every H orientation has slope above −0.5, and the eight improvements range from 1.99979001 to 2.00018732 across the two fixtures. For those eight divergent orientations, the minimum `R²` exceeds 0.9999998186. The other 22 are not required to have high-quality soft power fits; their minimum H-fit `R²` is approximately 0.556 (child) and 0.862 (full). Thus the measured improvement is approximately two powers, stronger than the test\'s required minimum of 1.5. The corresponding raw fitted exponents are approximately −5.0011 (U) and −3.00057 (H); adding the same measure power leaves their two-power difference unchanged. This is a physical scaling improvement, not a runtime speedup.

The generic CLI IR verdict requires strictly positive scaling and reports `passed=false` / `all_passed=false` for these negative, near-zero summed H fits. The Rust tests instead require H scaling above −0.5 and improvement of at least 1.5 powers, which they satisfy. The result establishes that targeted improvement near the logarithmic boundary; it does not establish strict IR convergence, pointwise boundedness or the generic CLI\'s positive-scaling criterion.

Across the child/full fixtures, all 19 saved soft-profile reports have `all_passed=false`: 15 summed reports and four localized per-orientation reports. Each orientation report has 22 strictly passing and eight nonpassing orientations; the same eight remain nonpositive after H despite their measured improvement. These generic verdicts are distinct from the passing Rust improvement assertions.

Mandatory-case per-route UV counts, vanishing counts, bare controls and soft fits are retained in #link("soft-ct-2026-09-14/acceptance-numerical.typ")[acceptance numerical results];. Full top-bubble H and U each pass 27 summed rays per route and 810 localized orientation/ray combinations; nested top passes 21 summed and 294 localized combinations; the double-triangle control passes 24 summed and 432 localized combinations.

== Projected Appendix B1 numerical failure
<projected-appendix-b1-numerical-failure>
The original Appendix B1 Rust acceptance stopped at its root-provenance assertion before the full numerical matrix. Six supplementary fresh CLI states completed that matrix without changing the test: full and bare integrands in all three routes, using the original card, matching saved settings, runtime configuration, model and profile identities. Every route uses the complete orientation sum. The full profiles use the standard 25-point, `10^8–10^12`, seed-1337 UV scan.

#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,right,right,auto,),
    table.header([Route], [Full summed UV fits], [Full orientation fits], [Matching bare control],),
    table.hline(),
    [Localized 3D], [21/21 pass], [294/294 pass], [17/21 summed fail (13 slope, 4 unstable); 178/294 orientations fail (100 slope, 78 unstable), as expected],
    [Explicit 3D], [21/21 pass], [Not separately requested], [17/21 summed fail (13 slope, 4 unstable), as expected],
    [Projected 4D → 3D], [#strong[8/21 pass; 13 fail];], [Not separately requested], [17/21 summed fail (13 slope, 4 unstable), as expected],
  )]
  , kind: table
  )

Six projected failures are gamma1 child-hard rays with positive slopes #strong[+1.0254 to +1.0662] (`R²=0.9911–0.9985`), against the required slope ≤−0.9. Seven all-hard fits are unstable. All 21 rays resolve, and all thirteen failures persist after arbitrary-precision retries (`ArbPrec`, 1,000 bits). This is an additional numerical acceptance failure; it cannot be dismissed as the root-provenance assertion or counted as an expected bare failure. #link("soft-ct-2026-09-14/appendix-b1-supplementary.typ")[Full matrix and fit evidence] record failure categories and route/configuration checks.

Two existing controls narrow the discrepancy. With only the ordinary U1 child retained, explicit and projected routes both pass all six target child-hard rays with slope approximately −1 and minimum `R² > 0.99999999999999`. The remaining fifteen limits are intentionally unsubtracted in this control and fail as expected. #link("soft-ct-2026-09-14/appendix-child-diagnostic.typ")[Child-only comparison] records the exact prescription changes and all six matched fits.

Restoring the second ordinary child while keeping the outer H2 disabled gives a shared divergence: both routes have slope approximately +2 on all six target rays, with `R²=1` and arbitrary-precision retries. All 21 recorded fits, including their sampled point data, match across routes in this partial forest. The two children overlap, so its structural forest contains the root and two child branches, with no product branch. This is a diagnostic partial forest, not a positive six-ray acceptance control. The #link("soft-ct-2026-09-14/appendix-two-child-diagnostic.typ")[two-child comparison] shows that the complete forest is where the route-dependent acceptance failure is observed. Matching partial-forest observations do not prove equal symbolic coefficients or exclude a smaller difference hidden beneath shared leading growth; they do not locate a defect solely inside H2.

Separately, ten matched canonical points were evaluated in the full and bare explicit/projected states at 1,000-bit precision. All ten bare pairs agree after f64 output serialization, and full values agree at the ordinary point. Hard-scaled full values diverge: scaling the second loop by `10^12` gives imaginary parts `−5.44049e−54` (explicit) and `+1.98819e−48` (projected). These are additional canonical rays, not the seeded candidate-LMB rays; raw point inspection does not report an accuracy bound. #link("soft-ct-2026-09-14/appendix-point-comparisons.typ")[Point comparison evidence] retains all values and limitations.

The full routes share the input graph, physics settings, requested subtraction and final 3D evaluator machinery; direct local-3D construction versus local-4D construction/projection is the first different pipeline. Bare agreement and the passing isolated U1 control locate the reproduced failure in the complete counterterm-dependent result beyond that isolated child. They do not prove which added term or representation step causes it. The next comparison must certify exact factorized numerator and denominator reconstruction for the first unequal contribution, retaining edge ownership and provenance. No such certificate has yet been established for this full-forest discrepancy, so no attribution to CFF recursion, contour signs or residue aggregation is made.

== Failures and incomplete validation
<failures-and-incomplete-validation>
No failing assertion was edited. The findings below locate the first differing boundary before proposing deeper algebra changes, following repository discrepancy guidance.

#figure(
  align(center)[#table(
    columns: 3,
    align: (auto,auto,auto,),
    table.header([Failing case], [Observed boundary], [Implication],),
    table.hline(),
    [`paper_appendix_b_6_b_7_production_outer_h2_matches_the_routed_rhs`], [Production reduces color in `grow_factor`; the paper oracle keeps raw color generators. The first S1 physical comparison retains a raw-chain versus Casimir/trace residual.], [Later UV-quadratic/full-H checks were not reached. Matching color normalization is the next comparison.],
    [`paper_appendix_b_10_b_12_massless_child_and_local_b_16_nested_branch`], [Same color boundary, already at massless B.11 physical S0.], [Nested B.16 assertions were not reached; this is not evidence of a B.16 operator error.],
    [`production_outer_u_grades_a_free_fermion_mass_but_keeps_muv_fixed`], [Raw numerator/denominator oracle disagrees with color-reduced production grow before any U rescaling.], [The current failure precedes mass grading.],
    [`top_bubble_actual_fermion_child_h2_removes_both_soft_jet_coefficients`], [Structural zero check leaves a factorized expression.], [Both captured constant/linear residuals cancel exactly under factor-preserving regrouping, but the Rust assertion still fails.],
    [`collective_three_component_region_keeps_only_its_divergent_join`], [Production excludes the convergent disconnected component; the manually constructed excluded-region oracle is a factorized zero that structural collection does not recognize.], [Four exact scalar coefficients are `[F,−F,−F,F]` with one immutable numerator. The certificate proves this captured cancellation, not a passing Rust test.],
    [`slow::paper_appendix_b1_ir_specialization_generates_and_profiles`], [Test requires nonempty `projection_paths` for every term; production\'s identity `Root` intentionally returns an empty list.], [The six-node forest and isolated-child checks passed; full UV checks and later routes were not reached by this test. This contract contradiction is present in earlier source history as well.],
    [`slow::paper_figure_b1_massless_bubble_uses_a_complete_soft_wood`], [Explicit-3D generation panicked in `bytes 1.11.1`: buffer length 1,576,152,287, requested advance 5,871,119,583.], [Failed after 6,326.706 s, with peak RSS 313.962 GiB. Localized checks passed; explicit UV checks and the projected route were not reached. No timeout occurred.],
  )]
  , kind: table
  )

For the three color cases, both routes share the same graph, routing and local operation; color representation first differs immediately after numerator construction. For the two exact-zero cases, diagnostics retained factorized tensor/graph numerators and distinct original denominator owners. No CFF recursion, contour-sign, integration-backend or evaluator change is justified by these local comparison failures.

Detailed source references, the full route comparisons, smallest reproducers and exact certificates are in #link("soft-ct-2026-09-14/local4d-failures.typ")[local 4D triage];, #link("soft-ct-2026-09-14/forest-failure.typ")[forest triage];, forest certificate;, top-bubble certificate and #link("soft-ct-2026-09-14/appendix-b1-failure.typ")[Appendix B1 triage];. The standalone proof scripts and original expressions remain under `/tmp`; they do not import or modify the Rust implementation.

#strong[Eight original Nextest executions abort at Symbolica startup.] Seven are the standalone Spenso/Idenso/Vakint cases; a retry using the repository CI fallback was rejected as “Unknown license.” A valid standalone license environment is still needed, and the retry is not seven additional unique cases. The eighth is `logarithmic_tilde_is_zero_so_hat_is_t`, which omits normal per-test initialization. A separate two-test libtest run initialized GammaLoop through an existing passing test, then executed that exact unchanged logarithmic assertion successfully (2/2, combined 0.061 s). This establishes a missing bootstrap boundary in isolated execution; its original Nextest abort remains recorded. No license extraction, source edit or rebuild was used. See #link("soft-ct-2026-09-14/core-initialized-diagnostic.typ")[initialization diagnostic];. No license key is retained in artifacts.

== Limits and remaining work
<limits-and-remaining-work>
A full green sign-off requires resolving the projected Appendix B1 numerical failure, the six comparison/provenance assertions, the large explicit-3D generation panic, and the isolated-test bootstrap issue, with maintainer agreement where tests change. The seven standalone support cases still require a valid Symbolica environment. All affected complete checks must then be rerun. Exact diagnostics narrow the work; they do not substitute for those executions.

The curated suite excludes the two pre-existing `uv::hedge_poset::tests::failing` diagnostics, unrelated slow amplitude cases, pySecDec backends and the broader workspace. No Python API end-to-end or Monte Carlo cross-section acceptance is claimed. Local OS remains deferred; integrated OS, IR/PolePart in the same wood, and a positive-degree MUV parent above a nontrivial positive-degree IR child remain deliberately rejected. Passing their rejection tests confirms the policy, not support for those modes. Unsupported gamma5/chiral contractions and unsupported backend mass patterns remain outside this result.

== Selected integration outcomes
<selected-integration-outcomes>
Updated 2026-09-14T11:08:47+00:00. Exactly 23 selected cases: 21 passed, 2 failed and 0 timed out.

Durations are exact Nextest JUnit `testcase@time` seconds; they include generation, profiling and teardown and are not throughput measurements. Three routes means localized direct 3D, explicit direct 3D and projected local 4D. All mandatory soft cases use local counterterms with integrated generation and threshold subtraction disabled; integrated addbacks have separate tests below.

#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,auto,right,auto,),
    table.header([Selected test], [Result], [Wall seconds], [Assertion / purpose],),
    table.hline(),
    [#link("../../tests/tests/uv.rs:2557")[local\_ir\_integrated\_generation\_succeeds\_across\_local\_uv\_routes];], [PASS], [35.984], [Generate integrated addbacks for quark/gluon self-energies, disconnected spectacles and nested top self-energy in 13 variants, requiring finite complete sums, nonzero spectacles/nested results and accuracy-aware route agreement.],
    [#link("../../tests/tests/uv.rs:2733")[genuine\_vacuum\_tadpole\_matches\_all\_local\_uv\_routes\_and\_vakint\_inputs];], [PASS], [4.058], [Check all three routes with tensor-integral projection on/off against a nonzero vacuum tadpole reference, Im = −9.760529078735244×10⁻⁴ within 10⁻¹².],
    [#link("../../tests/tests/uv.rs:2904")[local\_os\_integrated\_generation\_is\_rejected\_by\_run\_generate];], [PASS], [0.730], [Require the intended integrated-OS policy error and no partially generated integrand.],
    [#link("../../tests/tests/uv.rs:2941")[local\_ir\_and\_pole\_part\_wood\_is\_rejected\_by\_run\_generate];], [PASS], [0.722], [Require contextual rejection of a mixed IR/PolePart wood and no partially generated integrand.],
    [#link("../../tests/tests/uv.rs:2994")[local\_ir\_child\_under\_positive\_degree\_muv\_parent\_is\_rejected\_by\_run\_generate];], [PASS], [0.672], [Require rejection of a positive-degree MUV parent over a positive-degree IR child, with the locality-policy explanation and no partially generated integrand.],
    [#link("../../tests/tests/uv.rs:4235")[dgse\_local\_ir\_profiles\_are\_non\_vacuous];], [PASS], [35.451], [Check DGSE forest identities, all 27 summed and 810 localized orientation UV fits, matched failing bare UV controls and resolved soft-scaling improvement across all three routes.],
    [#link("../../tests/tests/uv.rs:4259")[dgse\_local\_ir\_generation\_orientation\_filter\_profiles\_one\_orientation];], [PASS], [3.761], [Verify that filtering one DGSE production orientation preserves its value and yields exactly 27 resolved, passing orientation UV fits.],
    [#link("../../tests/tests/uv.rs:4411")[massless\_quark\_self\_energy\_local\_ir\_profile\_is\_non\_vacuous];], [PASS], [3.927], [Identify the intended quark self-energy soft component and require passing subtracted UV profiles with matched failing bare controls in all three routes.],
    [#link("../../tests/tests/uv.rs:4486")[massless\_gluon\_self\_energy\_local\_ir\_profile\_is\_non\_vacuous];], [PASS], [4.216], [Identify the degree-two gluon self-energy soft component and require passing subtracted UV profiles with matched failing bare controls in all three routes.],
    [#link("../../tests/tests/uv.rs:4685")[soft\_cff\_state\_survives\_a\_fresh\_process\_reload];], [PASS], [3.565], [Save and reload generated soft CFF states in fresh processes for the explicit 3D and projected 4D routes, checking regeneration, forest metadata and finite evaluation.],
    [#link("../../tests/tests/uv.rs:4732")[massive\_top\_self\_energy\_local\_os\_dispatch\_is\_deferred];], [PASS], [0.636], [Catch the exact intentional deferred-OS panic for massive top self-energy; passing confirms rejection, not an implemented OS subtraction.],
    [#link("../../tests/tests/uv.rs:4757")[paper\_appendix\_b1\_os\_child\_reaches\_deferred\_dispatch];], [PASS], [0.632], [Catch the exact intentional deferred-OS panic for the Appendix B.1 child topology.],
    [#link("../../tests/tests/uv.rs:4788")[local\_ir\_disconnected\_completed\_components\_form\_nonzero\_union];], [PASS], [11.977], [Check two completed H₂ components and their nonzero neutral union, with H and MUV controls passing matched UV limits and bare controls failing across three routes.],
    [#link("../../tests/tests/uv.rs:4966")[consecutive\_soft\_components\_have\_exact\_forest\_and\_cli\_profile];], [PASS], [8.409], [Require the exact four-node forest for two consecutive H₂ bubbles and their nonzero union, with passing subtracted and matched failing bare UV profiles across three routes.],
    [#link("../../tests/tests/uv.rs:5123")[sunrise\_pole\_part\_matches\_muv\_inspect];], [PASS], [4.703], [Generate scalar sunrise with MUV and PolePart prescriptions under the Compare orchestrator and require summed inspection agreement at one deterministic momentum point to relative tolerance 10⁻¹⁰.],
    [#link("../../tests/tests/uv.rs:1249")[scalar\_amplitudes\_match\_across\_local\_uv\_routes];], [PASS], [5.295], [Compare finite, nonzero ordinary local-UV amplitudes at ten points across six scalar graphs in three routes, plus an explicit-3D/projected-4D Arb comparison retaining about 125 matching bits.],
    [#link("../../tests/tests/uv.rs:1698")[scalar\_spectacles\_integrated\_matches\_across\_local\_uv\_routes];], [PASS], [5.957], [Compare nonzero ordinary integrated-addback evaluations for spectacles and its overall-logarithmic variant in three routes at two fixed points with f64 relative tolerance 10⁻¹⁴.],
    [#link("../../tests/tests/uv.rs:5556")[slow::paper\_figure\_b1\_double\_triangle\_is\_a\_non\_vacuous\_uv\_fixture];], [PASS], [58.997], [Use the ordinary Figure B.1 double triangle as a control, requiring a nonempty UV forest and non-vacuous passing subtracted/failing bare UV profiles in all three routes.],
    [#link("../../tests/tests/uv.rs:5621")[slow::local\_nested\_soft\_fixture\_generates];], [PASS], [19.888], [Require the exact H₁ child/parent chain and four-node nested forest, with passing subtracted and matched failing bare UV profiles in all three routes.],
    [#link("../../tests/tests/uv.rs:6165")[slow::gamma\_star\_ddbar\_top\_bubble\_child\_only\_has\_two\_power\_soft\_improvement];], [PASS], [57.020], [Isolate the degree-two massive-top child and require matched U→H soft improvement of at least 1.5 powers in three routes, including all 30 localized orientations and exactly eight U orientations with reported scaling below −1.],
    [#link("../../tests/tests/uv.rs:5747")[slow::paper\_appendix\_b1\_ir\_specialization\_generates\_and\_profiles];], [FAIL], [72.271], [Check two ordinary UV children under an H₂ parent, six isolated child-hard rays, six-node computed provenance and full-forest UV profiles; execution stopped at a localized provenance assertion before the full-forest profiles.],
    [#link("../../tests/tests/uv.rs:6300")[slow::gamma\_star\_ddbar\_top\_bubble\_has\_two\_power\_soft\_improvement];], [PASS], [198.364], [Require ordinary-U and completed-H full forests to pass matched UV controls and retain the top-child soft improvement after outer UV subtraction across all three routes.],
    [#link("../../tests/tests/uv.rs:5527")[slow::paper\_figure\_b1\_massless\_bubble\_uses\_a\_complete\_soft\_wood];], [FAIL], [6326.706], [Check the central massless H₂ bubble and containing wood with matched UV controls; localized profiles completed, then explicit generation panicked in bytes 1.11.1 before its UV profile and before the projected route.],
  )]
  , kind: table
  )

Profile exponents include the measure shift: UV reports fit lambda^(3#emph[n\_hard) times the evaluated magnitude, and IR reports export s=b+3];n\_soft. Both top-bubble tests pass their targeted two-power-improvement criteria, while their six summed H reports and eight H orientations per localized fixture fail the separate strict CLI infrared criterion s\>0. DGSE uses the same convention with the stronger bound sH\>3.

The selected scalar integrated-spectacles test exercises the f64 path only; its separate ignored 1000-bit Arb diagnostic was not selected and documents a roughly 4×10⁻¹⁵ localized/explicit discrepancy. Integrated-addback point comparisons do not establish a Monte Carlo integrated cross section.

Evidence files: `/tmp/soft-ct-validation-2026-09-14/{regular,acceptance}-<test-name>.junit.xml` (the `slow::` prefix is omitted from filenames).

== Build and generation performance
<build-and-generation-performance>
All relevant repaired Cargo checks and builds completed. `cargo fmt --all -- --check` passed after the API repairs. Targeted Clippy passed for `gammalooprs` and the UV integration target, with one existing `too_many_arguments` warning at #link("../../crates/gammalooprs/src/uv/profile.rs:3096")[profile.rs:3096];. Cargo also reported the existing future-incompatibility warning for `proc-macro-error2 2.0.1`. This is not an all-workspace, warnings-denied lint claim.

#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,right,right,auto,),
    table.header([Command stage], [Wall time], [Largest-child peak RSS], [Result],),
    table.hline(),
    [Repaired core check], [31.453 s], [2.278 GiB], [Pass],
    [Repaired UV integration check], [1.756 s], [0.316 GiB], [Pass],
    [Spenso/Idenso/Vakint check], [7.978 s], [0.579 GiB], [Pass],
    [CLI check], [33.487 s], [1.588 GiB], [Pass],
    [Targeted Clippy], [57.463 s], [1.604 GiB], [Pass, warning above],
    [Core/support test compilation], [564.536 s], [10.321 GiB], [Pass],
    [UV integration test compilation], [328.263 s], [9.528 GiB], [Pass],
    [CLI compilation], [274.851 s], [8.477 GiB], [Pass],
  )]
  , kind: table
  )

These are observed command costs with an evolving build cache. The initial failed check spent 277.232 seconds compiling/checking dependencies before finding the obsolete unit-test call; it is separate from the repaired checks. Some setup commands overlapped and waited for Cargo locks, so summing their walls does not give campaign elapsed time. Peak RSS is the largest observed child-process high-water mark, not total concurrent memory. Optimized dependency/test compilation dominated setup cost.

The successful integrated-IR test recorded 13 generation summaries; the vacuum test recorded six. Approximate sums of the displayed per-graph expression, Spenso-construction, Symbolica-evaluator-construction and compilation stages are:

#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,right,right,right,),
    table.header([Integrated fixture], [Localized 3D], [Explicit 3D], [Projected 4D → 3D],),
    table.hline(),
    [Quark self-energy], [0.204 s], [0.179 s], [0.177 s],
    [Gluon self-energy], [0.215 s], [0.212 s], [0.317 s],
    [Scalar spectacles], [0.832 s], [0.799 s], [0.695 s],
    [Nested top self-energy], [6.830 s], [6.110 s], [4.939 s],
  )]
  , kind: table
  )

Nested explicit 3D with integrated tensor projection disabled takes approximately 7.472 seconds of displayed stages. The projected nested route spends 3.85 seconds constructing expressions and 1.089 seconds constructing Spenso/Symbolica evaluators; explicit 3D spends 3.61 and 2.500 seconds respectively. In this run, reduced evaluator-construction cost explains the lower projected total. The gluon fixture shows that projected generation is not universally faster.

Vacuum expression construction spans 59--120 ms across its six variants; both evaluator-construction stages are sub-millisecond. Rounded `0ms` cells can mean sub-millisecond work. The original “Symbolica eval” column measures evaluator #strong[construction];, not repeated evaluation. Backend code compilation is disabled for these acceptance cards.

Generation memory is sampled every 100 ms in the process. The nested localized/explicit/projected entries are 202.68/213.34/129.59 MiB, but routes share a process and retain allocator/cache state. These values cannot establish intrinsic per-route memory savings. Full stage tables, model/route associations, precision settings and measurement caveats are in #link("soft-ct-2026-09-14/integrated-generation.typ")[integrated generation performance];.

The complete-soft-wood Figure B1 acceptance #strong[failed during explicit-3D generation after 6,326.706 seconds (105.445 minutes)];. The wrapper wall was 6,327.581 seconds, with 6,416.209 user CPU seconds and 18,080.822 system CPU seconds. Largest-child peak RSS was #strong[329,213,524 KiB = 313.962 GiB];. The captured panic was `advance out of bounds: the len is 1576152287 but advancing by 5871119583`, in `bytes 1.11.1` at `src/lib.rs:170`; the test propagated the worker panic at #link("../../tests/tests/uv.rs:5549")[uv.rs:5549];. This was a generation failure, not a timeout or a reported out-of-memory kill. No backtrace was captured, and the exact upstream caller is not proven by this captured log;. The two lengths differ by exactly `2^32`; the pinned Symbolica source contains unchecked function/product length casts to u32. This supports a 32-bit encoded-length truncation hypothesis, but does not prove the originating atom or call site without the missing backtrace. Panic evidence records the exact arithmetic and source references.

This is a one-graph, three-loop/eight-internal-edge workload, compared with two loops/five internal edges in the ordinary double-triangle control. Its localized route intentionally filters to one audited nonzero orientation, whereas explicit 3D and projected 4D request the full source-orientation sum. Localized versus explicit costs therefore represent unequal work. The localized 196 summed and 196 single-orientation UV fits passed. Explicit generation produced no completed generation summary or UV profile, and projected generation was not reached. The sampled resource history and bounded CPU profiles describe the failed explicit preparation, not a successful runtime evaluation.

A bounded CPU profile during this long generation reached final single-parametric evaluator preparation. In a 15-second, 49 Hz user-CPU stack sample (818 samples, zero lost), about 74.94% was under `parametrize_residue_map_selectors`/`AtomView::contains_symbol`, with 62.96% self time in the containment search. About 24.45% was in sysinfo/Rayon process-and-task enumeration used by the RAM monitor. These are sampled user-CPU fractions, not whole-run wall-time attribution; inclusive stack rows overlap. An earlier 352-sample window had 45.74% in memory copying, whose caller was not captured.

Source inspection finds two concrete investigation points. Explicit-sum construction strips orientation selectors before shared evaluator construction, which nevertheless scans for those selectors again; absence after tensor preprocessing still needs certification before skipping the scan. Separately, the pinned sysinfo implementation traverses host processes and their tasks before applying the monitor\'s PID filter. Neither observation justifies a claimed speedup without a matched measurement. #link("soft-ct-2026-09-14/large-fixture-performance.typ")[Detailed performance diagnosis] records the exact call path, existing ownership fast paths, remaining uncertainty and the next comparison. Profiler notes document sampling and discarded profiler attempts. The test itself was never paused or modified by these samples.

The six benchmark states were also generated in separate fresh CLI processes, with complete sums in every route. Preparation ran sequentially alongside the long Figure B1 acceptance on the shared host; timed evaluation ran after the acceptance campaign and all supplementary generation jobs finished. Their full-precision graph-generation summaries give:

#figure(
  align(center)[#table(
    columns: 6,
    align: (auto,right,right,right,right,right,),
    table.header([Fixture / route], [Reported graph generation], [Expression], [Spenso build], [Symbolica build], [Sampled peak RAM],),
    table.hline(),
    [nested-integrated-explicit], [6.370 s], [3.806 s], [1.841 s], [0.724 s], [144.2 MiB],
    [nested-integrated-localized], [7.331 s], [4.101 s], [1.884 s], [1.346 s], [197.4 MiB],
    [nested-integrated-projected], [4.886 s], [3.819 s], [0.848 s], [0.219 s], [102.5 MiB],
    [top-child-local-explicit], [9.259 s], [1.865 s], [2.018 s], [5.376 s], [464.9 MiB],
    [top-child-local-localized], [9.540 s], [1.815 s], [1.564 s], [6.161 s], [476.7 MiB],
    [top-child-local-projected], [6.805 s], [2.897 s], [3.084 s], [0.824 s], [260.3 MiB],
  )]
  , kind: table
  )

These totals aggregate the recorded per-graph preprocessing and integrand-construction stages. They exclude CLI startup, state serialization and unreported surrounding process work; they are not full command wall times. “Expression” is the production residual `total − Spenso build − Symbolica build − compile`. Backend compilation is zero for all six. Unlike the earlier multi-route test memory table, these fresh CLI processes do not retain earlier route allocations; the 100 ms sampling limitation and shared-host timing limitations still apply. Exact durations and original summaries are in benchmark generation data;.

== Fixed-point evaluation throughput
<fixed-point-evaluation-throughput>
Two fixtures were generated in fresh CLI processes: integrated nested top self-energy and local child-only top bubble. Each ran all three local routes. A read-only graph audit verifies identical graph structures, runtime settings and ordered loop coordinates within each trio: nested `(k0,k1)=(e2,e3)`, top child `(e6,e8)`. See coordinate audit;. The benchmark repeatedly evaluates the complete summed momentum-space integrand at `(0.11, -0.29, 0.37, 0.59, -0.43, 0.71)`. Each measurement uses ten warmup samples followed by ten batches. Initial nominal five-second targets produced only 0.167--6.111 seconds of timed samples because the ten-sample calibration amortizes per-request setup over fewer samples than the larger timed batches. Setup occurs inside both requests. A same-process primer pilot confirmed that setup repeats. The nine windows below three seconds were rerun with proportionally increased nominal targets (capped at 90 seconds), aiming for roughly five seconds of measured samples; the three configured top-child windows were retained. Original results are retained under `benchmark-initial`, and repeat metadata records the adjustment. The table reports actual timed windows, not requested duration.

“Raw” is the existing `--minimal-integrand` mode: default double precision, caches/events/selectors/observables and stability rotations/escalation disabled. “Configured” retains the card runtime settings. The reported uncertainty is the standard error across equally weighted batch means (sample standard deviation divided by the square root of ten), not a confidence interval over independent machines or input points. The ten warmup evaluations are excluded from timed sample totals; batch sizes differ by at most one. Throughput is the reciprocal of mean sample wall time. The measured Total includes evaluation call overhead; nested category timers are not added twice. The primary-evaluator column records only the first precision level and its primary rotation. Additional rotations and higher-precision retries are excluded by the timer gate and can appear under “other/overhead”; that category is not all framework overhead. Total includes the complete sample cost.

#figure(
  align(center)[#table(
    columns: 7,
    align: (auto,auto,right,right,right,right,right,),
    table.header([Fixture / route], [Mode], [Total µs/sample ± SE], [Primary evaluator µs], [Samples/s], [Timed samples], [Timed window s],),
    table.hline(),
    [nested-integrated-explicit], [raw], [40.107 ± 0.492], [31.235], [24,933], [119,194], [4.781],
    [nested-integrated-explicit], [configured], [74.130 ± 0.665], [30.053], [13,490], [55,576], [4.120],
    [nested-integrated-localized], [raw], [140.311 ± 2.416], [128.665], [7,127], [149,980], [21.044],
    [nested-integrated-localized], [configured], [281.570 ± 3.504], [129.651], [3,552], [17,800], [5.012],
    [nested-integrated-projected], [raw], [18.797 ± 0.057], [12.438], [53,200], [182,812], [3.436],
    [nested-integrated-projected], [configured], [35.149 ± 0.260], [12.301], [28,450], [107,762], [3.788],
    [top-child-local-explicit], [raw], [122.012 ± 0.582], [111.672], [8,196], [52,398], [6.393],
    [top-child-local-explicit], [configured], [33493.111 ± 1052.403], [284.072], [30], [150], [5.024],
    [top-child-local-localized], [raw], [364.303 ± 2.620], [346.750], [2,745], [13,373], [4.872],
    [top-child-local-localized], [configured], [92275.330 ± 5872.305], [769.907], [11], [67], [6.111],
    [top-child-local-projected], [raw], [42.675 ± 0.335], [33.684], [23,433], [137,736], [5.878],
    [top-child-local-projected], [configured], [8469.332 ± 243.074], [75.482], [118], [601], [5.091],
  )]
  , kind: table
  )

At the stable nested point, the projected route has the lowest observed configured cost: 35.149 µs/sample, versus 74.130 for explicit 3D and 281.570 for localized 3D. The observed reciprocal-time ratios are 2.11× and 8.01× respectively. These compare current routes at one point; they are not a measured rebase speedup.

All six first-warmup returned values for nested-integrated are finite and nonzero; maximum relative difference across routes and modes is 7.76603e-16. All six first-warmup returned values for top-child-local are finite and nonzero; maximum relative difference across routes and modes is 6.15084e-09. The JSON returned value is the first warmup evaluation; the benchmark does not persist or check every timed result. These point comparisons are a non-vacuity/agreement check, not an arbitrary-precision certificate. The short preliminary command smoke is excluded from these timings.

These are fixed-point microbenchmarks on the shared host, run after the acceptance campaign. Fixed route order and cache evolution remain possible confounders. They do not measure integration convergence, random-point distributions, throughput at every momentum point or a pre-/post-rebase change. See benchmark metadata and individual JSON files for batch times, returned values and category breakdowns.

== Configured benchmark stability
<configured-benchmark-stability>
The configured top-child timings include repeated precision escalation and should be read as the cost of evaluating a point that exhausts the configured stability checks. Six supplementary read-only, single-point inspections used the same saved states, full summed selectors, canonical momenta and runtime settings as the benchmarks, with process debug logging enabled. Each inspection returned exactly the same f64-serialized complex value as its configured benchmark\'s first warmup result. All three nested-integrated routes were accepted at f64. All three top-child routes attempted f64, f128 and 1,000-bit Arb and remained unstable at the final level.

For the top-child fixture, the Arb identity and Euler-rotated results have norms `2.0929513612552366e-8` and `2.025300031403487e-8`. The implemented check compares each norm to their mean: `max(abs(norm - mean_norm) / mean_norm) = 0.01642719771122922`, or #strong[1.64271977%];, against the configured #strong[0.001%] bound. The same mismatch persists in Quad and Arb and in all three routes. Unrotated agreement and finite returned values therefore do not establish stability under the configured rotation check.

The evaluator timer records only the primary rotation at the first attempted precision level (#link("../../crates/gammalooprs/src/integrands/process/mod.rs:2949")[timing gate];). Further rotations, higher-precision work and API setup largely enter `other/overhead`, computed as the residual of total wall minus recorded categories (#link("../../crates/gammaloop-api/src/commands/bench.rs:573")[benchmark accounting];). The large difference between configured Total and evaluator time is consequently not all framework overhead. These matching inspections establish the exercised stability path for the diagnostic point; the benchmark does not persist per-sample stability histories, so attributing every timed repetition to that path remains an inference from the unchanged deterministic input/configuration.

The top-child logs also warn that its external fixed-helicity photon is off shell (`p²=90000`, expected `0`) and that such states may be ill-defined at this point. This is a relevant fixture limitation, not a demonstrated cause of the rotation mismatch. The result is shared by all three routes and does not diagnose the separate projected Appendix failure or alter the Rust test counts. Stability evidence retains exact commands, values, precision log lines and norm calculations, with the six complete diagnostic logs and JSON outputs beside it. These logged single-point inspections are not timing samples.

== Reproducing and auditing this run
<reproducing-and-auditing-this-run>
Selection metadata contains the exact 122-case core filter, seven supporting case names and all 23 integration cases. Command metadata records the actual commands, recorded process timings and exit codes, including unsuccessful compile/parser attempts and expected-failing bare/partial controls. Those command failures are distinct from the unique Rust testcase results. Environment metadata records the tools and resource settings; license secrets are excluded.

For a selected integration case, the essential invocation is:

```bash
cargo nextest run --locked --offline --cargo-profile dev-optim \
  --profile soft_ct_acceptance --retries 0 --no-fail-fast --no-capture \
  -p gammaloop-integration-tests --test uv \
  -E 'test(=slow::gamma_star_ddbar_top_bubble_has_two_power_soft_improvement)'
```

Use `ci_gammaloop` for the listed regular integration cases and `test_gammaloop` for the curated unit/support filters. The campaign supplied a copied #link("soft-ct-2026-09-14/nextest.toml")[nextest configuration];, with only JUnit output paths redirected to the evidence directory. Each invocation received its own generated-state and log directory. Set `GAMMALOOP_TESTS_NO_CLEAN_STATE=1` to retain those states and use fresh `TESTS_GAMMALOOP_STATE_PATH` / `GL_TEST_LOG_DIR` paths when repeating a case. Restore the captured Nix/native-library environment and available Symbolica setup before comparing results.

The standalone CLI diagnostics and microbenchmarks use saved states from fresh generation, not test-state files whose settings may predate in-memory overrides. Their commands retain the exact global route/subtraction settings. Inline momentum coordinates must be whitespace-separated. When replaying generation, use `save state --path /absolute/intended/state` explicitly: restored settings can carry `state.folder="."`. Bare save commands in this campaign left additional workspace copies; those generated copies were verified absent from the input revision and preserved under `/tmp/soft-ct-validation-2026-09-14/cwd-state-spill`, outside the final source/documentation change. The intended benchmark states were separately verified and successfully loaded read-only. `commands.json` links to the original state/output paths under `/tmp/soft-ct-validation-2026-09-14`; the human-readable fit summaries and diagnoses are retained beside this report; raw certificates, JUnit results and benchmark JSON are not distributed. Large generated states and full debug traces remain under `/tmp` and are not included in the repository diff.

#link("soft-ct-2026-09-14/test-cases.typ")[Every selected testcase and its result];, #link("soft-ct-2026-09-14/unit-case-purposes.typ")[unit purposes] and #link("soft-ct-2026-09-14/integration-case-purposes.typ")[integration purposes] provide the detailed audit trail. The result concerns the recorded source snapshot and sampled fixtures; it does not establish correctness of untested modes or a historical performance change.
