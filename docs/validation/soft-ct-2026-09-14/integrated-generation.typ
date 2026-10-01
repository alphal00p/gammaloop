= Integrated soft-counterterm generation performance

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.
<integrated-soft-counterterm-generation-performance>
Both completed integration tests passed: #strong[35.984 s] for the 13 integrated-IR variants and #strong[4.058 s] for the 6 vacuum variants (JUnit testcase wall time). Overall command walls were 36.945 s and 5.256 s. The generation summaries below account for graph/evaluator construction, not the full test wall time. No tests were rerun for this extraction.

All 19 summaries report one generation core. Localized 3D uses explicit\_orientation\_sum\_only=false/local4D=false; explicit 3D uses true/false; local4D uses true/true. Every variant enables integrated UV counterterms and produces a final 3D integrand. “Integrated projection” refers independently to project\_integrated\_uv\_cts\_onto\_tensor\_integrals. Route association follows the source loops, verified against exact summary counts, matching graph names, the expected state directories, and local\_projection\_total events appearing only in the assigned local4D variants.

The saved global\_settings.toml files are pre-override snapshots and #strong[do not record these in-memory route settings];. They must not be treated as a faithful replay of the test without restoring the source overrides.

== Integrated IR routes
<integrated-ir-routes>
The loop order in #link("../../../tests/tests/uv.rs:2557")[uv.rs] is quark, gluon, spectacles, nested top; each runs localized 3D, explicit 3D, then local4D. Nested top additionally runs explicit 3D with integrated tensor projection disabled. Quark/gluon/nested top use sm-default; scalar spectacles uses scalars-default.json. Quark uses the compare orchestrator; the other fixtures use hedge\_poset.

#figure(
  align(center)[#table(
    columns: 10,
    align: (right,auto,auto,auto,right,right,right,right,right,right,),
    table.header([\#], [Fixture], [Local route], [Integrated projection], [Evaluators], [Expression], [Spenso build], [Symbolica build], [Compile], [Sampled RAM],),
    table.hline(),
    [1], [Quark self-energy], [Localized 3D], [on], [2], [192ms], [7ms], [5ms], [0ms], [40.08 MiB],
    [2], [Quark self-energy], [Explicit 3D], [on], [1], [173ms], [5ms], [1ms], [0ms], [42.31 MiB],
    [3], [Quark self-energy], [Local4D], [on], [1], [169ms], [7ms], [1ms], [0ms], [43.62 MiB],
    [4], [Gluon self-energy], [Localized 3D], [on], [1], [186ms], [24ms], [5ms], [0ms], [44.75 MiB],
    [5], [Gluon self-energy], [Explicit 3D], [on], [1], [186ms], [23ms], [3ms], [0ms], [45.25 MiB],
    [6], [Gluon self-energy], [Local4D], [on], [1], [235ms], [80ms], [2ms], [0ms], [46.83 MiB],
    [7], [Scalar spectacles], [Localized 3D], [on], [1], [748ms], [54ms], [30ms], [0ms], [48.34 MiB],
    [8], [Scalar spectacles], [Explicit 3D], [on], [1], [731ms], [53ms], [15ms], [0ms], [43.74 MiB],
    [9], [Scalar spectacles], [Local4D], [on], [1], [668ms], [22ms], [5ms], [0ms], [47.08 MiB],
    [10], [Nested top self-energy], [Localized 3D], [on], [1], [3.70s], [1.82s], [1.31s], [0ms], [202.68 MiB],
    [11], [Nested top self-energy], [Explicit 3D], [on], [1], [3.61s], [1.80s], [700ms], [0ms], [213.34 MiB],
    [12], [Nested top self-energy], [Local4D], [on], [1], [3.85s], [868ms], [221ms], [0ms], [129.59 MiB],
    [13], [Nested top self-energy], [Explicit 3D], [off], [1], [4.72s], [2.01s], [742ms], [0ms], [187.19 MiB],
  )]
  , kind: table
  )

The nested-top displayed stage sums are approximately #strong[6.830 s localized 3D, 6.110 s explicit 3D, 4.939 s local4D, and 7.472 s explicit 3D with integrated projection off];. For the projected-integrated comparison, local4D has a slightly larger expression stage than explicit 3D (3.85 s vs 3.61 s) but substantially smaller Spenso plus Symbolica evaluator construction (1.089 s vs 2.500 s). Its displayed total is about 19% lower than explicit 3D and 28% lower than localized 3D in this run. These are within-run observations, not a pre-/post-rebase performance comparison. The local4D route is not universally fastest: the gluon example displays 317 ms of stages versus 212 ms explicit 3D.

Each route also passed finite arbitrary-precision point evaluation and agreement of complete physical sums with the localized reference, using 1000 times the reported numerical relative accuracy (with an f64 serialization floor). Spectacles and nested top must be nonzero. The quark and gluon point values may vanish. Expected IR-component identities/counts are asserted. These are deterministic point checks, not Monte Carlo integration or evaluation-throughput benchmarks.

== Genuine vacuum tadpole routes
<genuine-vacuum-tadpole-routes>
The #link("../../../tests/tests/uv.rs:2732")[vacuum test source] runs the three local routes first with integrated tensor projection off and then with it on. This is a one-loop scalar self-loop with numerator 1 and mass 2, no external momenta/helicities, MUV subtraction, and m\_uv=mu\_r=localization scale=1000.

#figure(
  align(center)[#table(
    columns: 9,
    align: (right,auto,auto,right,right,right,right,right,right,),
    table.header([\#], [Local route], [Integrated projection], [Evaluators], [Expression], [Spenso build], [Symbolica build], [Compile], [Sampled RAM],),
    table.hline(),
    [1], [Localized 3D], [off], [1], [120ms], [0ms], [0ms], [0ms], [40.38 MiB],
    [2], [Explicit 3D], [off], [1], [108ms], [0ms], [0ms], [0ms], [42.02 MiB],
    [3], [Local4D], [off], [1], [97ms], [0ms], [0ms], [0ms], [44.29 MiB],
    [4], [Localized 3D], [on], [1], [59ms], [0ms], [0ms], [0ms], [44.93 MiB],
    [5], [Explicit 3D], [on], [1], [62ms], [0ms], [0ms], [0ms], [44.98 MiB],
    [6], [Local4D], [on], [1], [79ms], [0ms], [0ms], [0ms], [45.02 MiB],
  )]
  , kind: table
  )

All six routes passed the physical-component/inert-root classification, nonzero finite point evaluation, mutual numerical agreement, and the reference result: real part within 1e-12 of zero and imaginary part -9.760529078735244e-4 within 1e-12. The table shows generation expression costs from 59 to 120 ms. Every Spenso and Symbolica build is submillisecond: “0ms” is truncation, not zero work (the original summaries report nonzero percentages). The second projection block benefits from later position in the same process, so its shorter times alone do not establish a causal projection speedup.

== Measurement definitions and limits
<measurement-definitions-and-limits>
- #link("../../../crates/gammalooprs/src/processes/generation_report.rs:39")[GraphGenerationStats] defines expression build as total generation time minus Spenso build, Symbolica evaluator build, and backend compile, using saturating subtraction. “Symbolica eval” in the original table is evaluator construction. Evaluator count is object count, not samples or events.
- #link("../../../crates/gammaloop-api/src/commands/generate.rs:400")[Duration formatting] truncates durations below one second to integer milliseconds and rounds seconds to two decimals. Sums of these displayed columns are approximate; underlying nanoseconds cannot be recovered from these logs. Backend compilation is disabled for these tests.
- #link("../../../crates/gammaloop-api/src/state.rs:614")[GenerationMonitor] samples only the current process memory every 100 ms. Its sample window restarts per generation, but the process, libraries, allocator state, and caches persist across routes. A late route can retain earlier allocations or release them. The 202.68/213.34/129.59/187.19 MiB nested-top entries therefore do not prove intrinsic per-route memory requirements. Sub-100ms spikes may be missed, child processes are excluded, and values are rounded to 0.01 MiB.
- Test-command max RSS is a separate runner measurement: 228192 KiB (222.844 MiB) for integrated IR and 104624 KiB (102.172 MiB) for vacuum. It includes command/process lifetime scope and must not be substituted for an individual generation RAM row.
- Host load was high: integrated-IR start/end load1 211.62/176.42, vacuum 176.42/185.92. The fixed route order, evolving caches, lack of repeated trials, and other host activity limit performance claims. No historical baseline was collected for these two cases.
- Overall test walls include fixture imports, generation setup, arbitrary-precision evaluation, assertions, and cleanup policy. Gaps between successive generation-summary timestamps are not standalone route-generation timings. No runtime samples/second are available from these logs.

== Artifacts
<artifacts>
- local\_ir\_integrated\_generation\_succeeds\_across\_local\_uv\_routes: JUnit;, runner metadata;, generation JSONL;.
- genuine\_vacuum\_tadpole\_matches\_all\_local\_uv\_routes\_and\_vakint\_inputs: JUnit;, runner metadata;, generation JSONL;.

The structured extraction retains all original stage cells including percentages, source JSONL line numbers/timestamps, full route flags, graph/model association, approximate displayed-stage sums, memory and command metadata.
