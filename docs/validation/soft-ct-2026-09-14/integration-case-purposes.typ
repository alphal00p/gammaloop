= Soft-counterterm acceptance inventory (read-only inspection, 2026-09-14)

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.
<soft-counterterm-acceptance-inventory-read-only-inspection-2026-09-14>
This document describes existing test assertions, not test outcomes. Sources are the current `tests/tests/uv.rs`, `tests/resources/run_cards/*.toml`, `.config/nextest.toml`, and the production profile/report schemas. No repository files were edited and no builds were started by the inventory agent.

== Route matrix and common numerical contract
<route-matrix-and-common-numerical-contract>
All six mandatory acceptance tests run three configurations in order:

#figure(
  align(center)[#table(
    columns: 3,
    align: (auto,right,right,),
    table.header([Label], [explicit\_orientation\_sum\_only], [local\_uv\_cts\_from\_expanded\_4d\_integrands],),
    table.hline(),
    [Localized direct 3D], [false], [false],
    [Explicit direct 3D], [true], [false],
    [Projected local 4D], [true], [true],
  )]
  , kind: table
  )

Final integrands are 3D in every case. These local fixtures explicitly require `generate_integrated=false`, with threshold subtraction disabled. Therefore 18 test/route combinations do not constitute validation of threshold subtraction or full integrated physical cross sections. The integrated-generation and vacuum tests below provide separate addback coverage.

The common CLI UV helper (uv.rs:3541) reads each card\'s `profile_uv`, adds `--selected-limits all`, and checks graph/LMB/subset/orientation identity against the matched bare/control run. Standard windows are 25 points from 10^8 to 10^12, seed 1337. DGSE uses 33 points from 10^8 to 10^16. Localized direct 3D profiles the sum and each production orientation. Explicit modes remove `--per-orientation` because they do not expose independently selectable production channels.

UV acceptance calls `pass_fail(-0.9)`: a resolved fit must have finite R² \>= 0.99 and slope \<= -0.9 (production uv/profile.rs:61,2218). Missing fits can be accepted only for certified vanishing rays when the profile explicitly permits it. Most cases additionally require at least one resolved summed and, where applicable, orientation fit; they do not universally require every fit to resolve. DGSE does require complete resolution. Bare controls must exhibit at least one resolved excessive-slope failure; this prevents vacuous success. Counts of selected, resolved, failed, and certified vanishing limits should be kept separate in the report.

Common pointwise route comparisons use the reported evaluation accuracy with f64 serialization floor and factor 1000: delta \<= 1000#emph[accuracy];scale, with 1000#emph[accuracy \< 1. Counterterm-visibility checks invert this inequality, requiring a difference \>= 1000];accuracy\*scale. These are accuracy-aware checks, not a fixed arbitrary relative tolerance.

== Six mandatory slow acceptances
<six-mandatory-slow-acceptances>
Configuration: `.config/nextest.toml:121`. Exact allow-list below, one nextest test at a time, zero retries via inherited CI profile, fail-fast disabled. Slow warning after 60 minutes, termination after two periods. The configuration explains that some symbolic-plus-numerical cases can exceed an hour. Test time includes repeated generation, symbolic forest exports, numerical profiling, and teardown.

#figure(
  align(center)[#table(
    columns: 3,
    align: (auto,auto,auto,),
    table.header([Existing test (`slow::` prefix)], [Fixture and purpose], [Existing expectations],),
    table.hline(),
    [`paper_figure_b1_massless_bubble_uses_a_complete_soft_wood` (uv.rs:5526)], [Three-loop double triangle containing central massless quark-loop gluon self-energy; complete containing soft wood], [Exactly one matching H2 central component plus a positive-degree containing IR component; 9 integration dimensions; independent parent-CFF orientation inventory; localized mode keeps deterministic first nonzero orientation, while stored source CFF retains all source multiplicities. Structural forest has \>=4 families. Finite matched bare/subtracted point differs above numerical accuracy. Subtracted UV profiles pass, matched bare UV profiles fail.],
    [`paper_figure_b1_double_triangle_is_a_non_vacuous_uv_fixture` (5564)], [Base Figure-B.1 double triangle as ordinary UV control], [\>=6 integration dimensions; \>=1 connected UV component; non-vacuous passing subtracted UV profile and failing matched bare control in all three routes.],
    [`local_nested_soft_fixture_generates` (5620)], [Two-loop nested top self-energy; nested H/H replay], [Exactly H1 child edges \[3,4\] and H1 parent \[2,3,4,5,6\], plus inert empty root. Exactly four forest nodes (root, child, parent, chain). \>=6 dimensions. Subtracted UV profiles pass and matched bare controls fail.],
    [`paper_appendix_b1_ir_specialization_generates_and_profiles` (5745)], [Nested gluon self-energy; ordinary UV children beneath H2 outer component], [Exactly U1 child \[3,4\], U0 child \[2,3,5,6\], H2 outer \[2,3,4,5,6\]. Six structural/computed forest nodes. Isolated U1 control tests exactly six child-hard rays: free e3/e4 × fixed e2/e5/e6, each slope \< -0.9. Computed direct-3D provenance checks unique node identity, parent existence, child-to-parent order, canonical/actual routes, and U(+1), S(+1), US(-1) branch coefficients. Full subtraction passes all UV limits and matches control identities; bare control fails.],
    [`gamma_star_ddbar_top_bubble_child_only_has_two_power_soft_improvement` (6164)], [Two-loop gamma\*→ddbar vertex containing massive-top gluon self-energy; isolates child without overall CT], [U-only and H variants each contain exactly degree-2 child edges \[7,8\], external PDGs \[-21,21\], internal \[6\], with no overall vertex CT. Identical routed soft rays. Summed U scaling \< -1, H scaling \> -0.5, H−U \>=1.5; both R²\>=0.98. In localized route all30 orientations are checked: every H slope \>-0.5; exactly8 U orientations have slope\<-1, each improved \>=1.5 with R²\>=0.98.],
    [`gamma_star_ddbar_top_bubble_has_two_power_soft_improvement` (6299)], [Full massive-top bubble vertex forest; proves outer ordinary UV subtraction preserves soft improvement], [Q=300 GeV, MT=173 GeV, WT=0, explicit color closure delta\_i^j/Nc, canonical LMB \[e6,e8\]. Compare bare, U2 child+U0 overall, H2 child+U0 overall. Soft forest exactly root/child/parent/chain. Both U-only and H variants pass UV, bare fails on matched limits. Summed soft U slope\<-1, H\>-0.5, improvement\>=1.5, all three fit R²\>=0.98 and identical ray fingerprints. Localized route checks same30/eight orientation behavior as isolated child.],
  )]
  , kind: table
  )

Top-bubble soft profiles use CLI `profile bulk`, graph `top_bubble_vertex`, limit `S(e6)`, 25 points 10^-2 to 10^-5, seed1337. Fingerprint requires LMB \[6,8\] and a16-character digest. These assertions accept an improvement of at least1.5 powers; the report should call approximately two powers rather than claim exact2.000 from the assertion alone. Top-bubble cards use `SingleParametric`, m\_uv=mu\_r=20, iterative orientation optimization disabled.

== Fourteen directly related regular integration tests
<fourteen-directly-related-regular-integration-tests>
The exact14 are listed in `.config/nextest.toml:158` with a serial nextest group and 30-minute slow period, two periods before termination in `ci_gammaloop`. The three-route matrix applies unless noted otherwise.

#figure(
  align(center)[#table(
    columns: 2,
    align: (auto,auto,),
    table.header([Test], [Coverage/expectations],),
    table.hline(),
    [`dgse_local_ir_profiles_are_non_vacuous` (4234)], [Matching massless quark soft component, computed forest identities. Exactly30 production orientations and27 hard limits, all summed fits resolved; localized route exactly810 orientation fits, all resolved and passing. Soft limit `se S(e0)` uses20 points10^-2→10^-3, seed1337; bare/subtracted R²\>=0.9 and identical routed ray; subtracted scaling\>3, bare\<=3, improvement\>=1. DGSE and top-bubble use the same reported scaling s=b+3\*n\_soft; both selected limits have n\_soft=1, so DGSE scaling\>3 means raw exponent b\>0, while the top-bubble test targets an approximately two-power improvement near s=0. Bare must fail UV.],
    [`dgse_local_ir_generation_orientation_filter_profiles_one_orientation` (4258)], [Localized direct route only. Physical orientation --00+++0+ at unfiltered index3 becomes visible slot0 after generation filtering. Values agree within1000\*reported accuracy; exactly27 selected/resolved orientation fits and no UV failure.],
    [`massless_quark_self_energy_local_ir_profile_is_non_vacuous` (4410)], [Correct IR identifier \[-1,1\]/\[1,21\]; passing local UV profiles with failing bare controls on matched routes.],
    [`massless_gluon_self_energy_local_ir_profile_is_non_vacuous` (4485)], [Correct IR identifier \[-21,21\]/\[1\]; degree-two massless-gluon local subtraction passes UV, bare fails, matched routes.],
    [`soft_cff_state_survives_a_fresh_process_reload` (4684)], [Two routes (explicit direct3D, projected4D), each spawned fresh save/load processes. Save generated soft CFF, reload, regenerate, validate computed forest identities, finite normal evaluation and metadata. Child test is ignored for standalone invocation and launched by parent.],
    [`massive_top_self_energy_local_os_dispatch_is_deferred` (4731)], [Pass means expected intentional deferred-OS panic was caught with exact message; does not mean OS works.],
    [`paper_appendix_b1_os_child_reaches_deferred_dispatch` (4756)], [Same expected intentional deferred-OS rejection for Appendix topology.],
    [`local_ir_disconnected_completed_components_form_nonzero_union` (4787)], [Dotted scalar spectacles: exactly two H2 components; neutral MUV union tag; root, two completed components, union are all nonzero. H and ordinary-MUV controls pass UV on same routes; bare fails.],
    [`consecutive_soft_components_have_exact_forest_and_cli_profile` (4965)], [Connected carrier has H2 bubbles \[0,1\] and\[3,4\], one neutral union \[0,1,3,4\], exactly four forest nodes in one wood, each nonzero. Passing UV vs failing matched bare.],
    [`local_ir_integrated_generation_succeeds_across_local_uv_routes` (2556)], [Integrated addback enabled for massless quark, massless gluon, disconnected scalar spectacles, nested top self-energy. Standard three routes plus nested explicit direct without tensor-integral projection (13 fixture/route variants). Correct component counts, finite evaluations; spectacles/nested values must be nonzero; arbitrary-precision evaluations of complete sums agree within1000\*reported accuracy. Quark/gluon total allowed zero; quark uncontracted atom has separate unit oracle. This tests generation/one-point agreement, not a Monte Carlo integrated cross section.],
    [`genuine_vacuum_tadpole_matches_all_local_uv_routes_and_vakint_inputs` (2732)], [Massive vacuum tadpole, all3 local routes × tensor projection on/off (6 variants). Physical degree2 vacuum component distinct from inert empty root. All values finite/nonzero and agree; reference Re≈0, Im=-9.760529078735244e-4 within1e-12.],
    [`local_os_integrated_generation_is_rejected_by_run_generate` (2903)], [Integrated OS policy rejection with corrective setting; no partial generated integrand.],
    [`local_ir_and_pole_part_wood_is_rejected_by_run_generate` (2940)], [Mixed IR/PolePart wood rejected with graph/cut/component context; no partial integrand.],
    [`local_ir_child_under_positive_degree_muv_parent_is_rejected_by_run_generate` (2993)], [Positive-degree MUV parent over positive-degree IR child rejected with explicit locality/scheme policy explanation; no partial integrand.],
  )]
  , kind: table
  )

== Performance and durable evidence
<performance-and-durable-evidence>
Set `GAMMALOOP_TESTS_NO_CLEAN_STATE=1` before execution or success deletes the numerical/profile artifacts (tests/src/lib.rs:287). Set `TESTS_GAMMALOOP_STATE_PATH` to a fresh report-specific directory to avoid stale artifacts; note no-clean also disables harness pre-run cleanup. State names encode `explicit_{bool}_project_local_4d_{bool}`. Repeated cases reuse some fixture names, so distinct test-run artifact roots are safest.

Nextest JUnit `testcase@time` gives test wall seconds, incorporating all route generations, profiles, symbolic exports, and other work. `.config/nextest.toml` writes mandatory results to `target/nextest/soft_ct_acceptance/junit.xml`. Treat nextest timeout/crash as incomplete/failing execution, and never interpret previous terminated JUnit as fresh evidence. For wall/CPU/RSS of the whole command, use external time/resource measurements if available. No throughput benchmark is inherent in this acceptance suite.

Generated numerical files:

- `<state>/local_ct_uv_profile/uv_profile.json`: `.scales`, `.graphs[].lmbs[].subsets[]`; each subset has `.free`, `.fixed`, `.initial_dod`, `.analysis.inspect_level.result.{slope,r_squared}`, optional `.per_orientation_inspect_entries[].analysis.result.{slope,r_squared}`, and fit-status details inside analyses.
- `<state>/<integrand>_soft_profile.json` and `<integrand>_soft_per_orientation_profile.json`: `.settings`, `.graphs[].single_limit_reports[]` with `limit_name`, `orientation_label`, `ray_fingerprint`, `scaling`, `r_squared`, `passed`.
- DGSE IR helper returns its fit in memory and does not request a soft JSON output file, unlike top-bubble helpers.

Generation stats schema (`processes/generation_report.rs`; API state.rs:90) can persist to `<state>/processes/amplitudes/<process>/<integrand>/generation_summary.json` when a generated state is saved: `peak_ram_bytes`, `reports[].{graph_name,integrand_name,stats}`. Stats include `evaluator_count`, `total_time`, `evaluator_spenso_time`, `evaluator_symbolica_time`, `evaluator_compile_time`. Rust `Duration` JSON uses seconds/nanoseconds. Expression construction time is max(0,total−spenso−symbolica−compile). The generation table logs the same quantities. Caveat: these tests generally save an empty state in `get_test_cli`, then generate without saving again. Generation summaries may therefore be absent except in persistence tests; do not claim stage timings from missing JSON. Most cards set both logging directives off, so enabling existing logging filters may be required to retain generation summaries.

Performance conclusions must distinguish current absolute test costs, current per-route generation costs (if recorded), and a true pre/post-rebase benchmark. There is no freshly run pre-rebase baseline supplied by the inventory. Route comparisons also change orientation representation; localized vs explicit differences are not attributable solely to the 4D implementation. State reuse, warmed compilation caches, machine contention, and repeated-run variance should be identified before claiming a speedup.
