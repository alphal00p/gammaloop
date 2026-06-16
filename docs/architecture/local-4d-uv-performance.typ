= Local 4D UV validation and performance
<local-4d-uv-validation-and-performance>
The implementation canonicalizes completed hard UV denominators before source reconstruction, merges compatible forest terms and reuses dynamically prepared numerator mappings. Its signed source certificates and CFF allocation are described in #link("architecture-current.typ")[the architecture];. The original plan and dated implementation history are preserved in #link("../../SPEED_UP_UV_CTS_FROM_4D.typ")[SPEED\_UP\_UV\_CTS\_FROM\_4D.md];.

The measurements below are the source-pinned performance evidence imported with `a1140c90c334ff58a3b040ff3eae06d3645bf0bb`. They belong to the exact numerical sources and binaries identified in their receipts; they are not fresh validation of the reconstructed stack or its review corrections. Current build and test outcomes belong to the #link("raised-energy-cff-stack-review-validation.typ")[stack validation receipt];. GL262 has an evaluator-construction panic and incomplete paired performance evidence; it is not counted as passing here.

At the measured source, the GL00/GL01 slowdown is reduced to the reported gates: generation premiums over erased 3D are 4.45% and 0.82%, and sampling is faster. This does not resolve every broader scalar case. The final scalar matrix retains some generation regressions with integrated and threshold counterterms enabled; earlier-source GL21 diagnostics identify a numerator-dependent cost, described separately below.

=== Measurement boundary
<measurement-boundary>
All physical generations, saved-state timings and six separate profiles in the main measurement cohort use numerical source `911f1824f994611d723dff0f50416cb301718617`, including the shared scalar-product and odd tensor-power fixes. All three routes were rebuilt and remeasured after those shared changes. The immutable CLI SHA-256 is `ef6811aa40114590a0fd72318f99f91850960239ff212abcaf272d7182e48691`; the matching evaluator-counter SHA-256 is `3e8912fdae918061719b0e8e6e8499ece14ced57a3482d5af7353677e9becaf8`. The `f938f9c` results remain a separately pinned cohort. Embedded receipts retain their capture-time status and source identity; proposed next actions in those artifacts are evidence context, not current execution instructions. Review corrections to zero projection, value-based tests and graph bookkeeping are not remeasured by these frozen binaries.

Physical measurements use the same model, physical numerator, external kinematics, LMB and eager, uncompiled evaluator settings in all three routes. Horner iterations are 1, common-pair-elimination rounds are 5, and worker limits are 1. Integrated and threshold counterterms are disabled for these isolated local-UV measurements; the scalar correctness matrix retains its enabled subtraction settings. The `dev-optim` build retains debug assertions. Generation and runtime run serially under an exclusive measurement lock, with no concurrent build or profiling. Each process has a 500 GB guard; no guard intervention occurred.

The #link("local-4d-uv-benchmarks.json")[benchmark receipt] contains all six representative parsed generation cards, saved global settings, complete runtime defaults and per-integrand overrides, the full model-parameter map, literal DOT inputs and card/state hashes. The two runtime-setting layers are both necessary because saved overrides omit default-valued fields. The constant external momenta are specified independently of saved `e_cm=300`: the incoming energy is 173.0, not 300. Both graphs use loop edges `[4, 9]` with their supplied signed routings. The DOT's `num=1` is an additional graph factor; the physical Standard Model Feynman numerator and overall factors −4/−8 are retained. Across routes, only explicit orientation summation and the local-4D route flag change.

Each route has three fresh-process generations. At each fixed momentum-space point, every saved state is primed before three accepted passes of twenty nonempty batches. The original point is `[0.11, -0.07, 0.19, -0.13, 0.05, 0.29]`; the second is its 100× counterpart. The existing minimal-integrand benchmark setting is unchanged. Ratios use route medians. Observed spread and twice the reported within-pass SEM supply conservative repeat triggers, not statistical confidence intervals.

The initial three-second and subsequent uniform fifteen-second rounds passed all median gates but retained inconclusive uncertainty envelopes from isolated long batches. Both rounds remain reported separately, without outlier removal. A predeclared matched confirmation repeated GL00/base, GL01/base and GL01/scaled across all three routes and all three saved generations, requesting 300 seconds to obtain at least 30 actual measured seconds per pass. GL00/scaled retains its complete fifteen-second cohort, initially requesting 100 seconds. Thus each final graph/point uses one complete matched cohort; no individual favorable passes are selected. The same sample accounting, formulas, workers and 1.15 threshold apply throughout. Undersized attempts are preserved and excluded by the existing rule.

A separate four-pass CPU diagnostic did not reproduce the earlier multi-second bursts: process CPU advanced for 112.10 of 112.92 observed wall seconds. It does not establish the cause of the earlier delays and supplies no acceptance samples. Profiling runs are likewise excluded from final timing gates.

The native generation summary distinguishes expression construction, Spenso tensor preprocessing, evaluator orchestration and optional compilation. `evaluator_symbolica_time` includes function-map/expression preparation, actual Symbolica builds and conversion to numerical programs. It is broader than the literal `.build()` call. Subtracting it retains tensor preprocessing and other graph work; subtracting both Spenso and orchestration gives the earlier expression boundary. Neither difference isolates only the UV forest. Separate profiling records literal build-call intervals and RAM at pre/post-Spenso and build entry.

=== Physical results at the final measured revision
<physical-results-at-the-final-measured-revision>
GL00 and GL01 are two-loop four-photon amplitudes with a top-quark loop and a top/gluon self-energy insertion. Each input has six interaction vertices and seven internal lines (six top-quark lines and one gluon), plus four external photons. They differ in external-photon ordering and momentum routing. The production forest exports each contain four computation nodes: the bare root, the full-graph DOD0 subtraction, the self-energy DOD1 subtraction and the nested self-energy/full-graph subtraction. Thus the count includes the bare state; it is not four distinct UV-divergent regions. Inputs are #link("../../examples/cli/aa_aa/2L/graphs/GL00.dot")[`GL00.dot`] and #link("../../examples/cli/aa_aa/2L/graphs/GL01.dot")[`GL01.dot`];.

GL00 and GL01 pass all ten generation/runtime gates against orientation-erased 3D. Each graph is assessed separately; the maximum accepted median ratio is 1.15. Exact observations, uncertainty envelopes and operation counts are retained in the #link("local-4d-uv-benchmarks.json")[benchmark receipt];.

#figure(
  align(center)[#table(
    columns: (20%, 26.67%, 26.67%, 26.67%),
    align: (auto,right,right,right,),
    table.header([Graph], [Generation 4D / erased], [Evaluator base / scaled], [Total sample base / scaled],),
    table.hline(),
    [GL00], [1.044523], [0.764529 / 0.811538], [0.791531 / 0.834682],
    [GL01], [1.008161], [0.783770 / 0.779277], [0.810855 / 0.800211],
  )]
  , kind: table
  )

All observed range and within-pass repeat-trigger envelopes are below 1.15; the largest runtime envelope upper endpoint is 1.052729. The phase split explains the remaining generation premium:

#figure(
  align(center)[#table(
    columns: (15.79%, 21.05%, 21.05%, 21.05%, 21.05%),
    align: (auto,right,right,right,right,),
    table.header([Graph / stage], [Localized 3D (s)], [Erased 3D (s)], [Direct 4D (s)], [4D / erased],),
    table.hline(),
    [GL00 / Expression construction], [3.064179], [3.081590], [1.933256], [0.627357],
    [GL00 / Spenso preprocessing], [1.880577], [1.883711], [3.381613], [1.795187],
    [GL00 / Native evaluator orchestration], [7.193893], [1.358823], [1.084717], [0.798277],
    [GL00 / Outside native orchestration], [4.962263], [4.965301], [5.520979], [1.111912],
    [GL00 / Complete graph generation], [12.126443], [6.324124], [6.605696], [1.044523],
    [GL01 / Expression construction], [3.102451], [3.111336], [1.948065], [0.626118],
    [GL01 / Spenso preprocessing], [1.904321], [1.912664], [3.384367], [1.769452],
    [GL01 / Native evaluator orchestration], [7.122215], [1.363317], [1.090730], [0.800056],
    [GL01 / Outside native orchestration], [5.007115], [5.026130], [5.343508], [1.063146],
    [GL01 / Complete graph generation], [12.130374], [6.392516], [6.444683], [1.008161],
  )]
  , kind: table
  )

Phase medians are independent and need not sum to the pipeline median. Expression construction is 37.26% faster for GL00 and 37.39% faster for GL01. Including tensor preprocessing but excluding native evaluator orchestration leaves 11.19% and 6.31% premiums respectively. These measurements do not show every preparation stage getting faster in 4D. All eighteen compilation timers are zero: generated evaluator code compilation is disabled.

Native graph time excludes CLI startup and state serialization: the localized/erased/direct child-wall medians are 13.008/8.005/8.005 seconds for GL00 and 13.007/8.005/8.005 seconds for GL01, with one-second supervisor polling. The 17.137-second native GL00 localized generation remains in the three-run sample; no generation outlier was removed.

#figure(
  align(center)[#table(
    columns: (15.79%, 21.05%, 21.05%, 21.05%, 21.05%),
    align: (auto,right,right,right,right,),
    table.header([Graph / route], [Evaluator μs/sample, base / scaled], [Total μs/sample, base / scaled], [Generation VmHWM (MB)], [Original instructions],),
    table.hline(),
    [GL00 / 3d\_local], [78.026 / 79.728], [87.531 / 89.437], [538.878--564.609], [87,768],
    [GL00 / 3d\_erased], [78.225 / 76.906], [87.596 / 86.554], [324.346--360.272], [25,840],
    [GL00 / 4d\_direct], [59.805 / 62.413], [69.335 / 72.245], [384.369--455.496], [20,210],
    [GL01 / 3d\_local], [75.248 / 75.363], [84.740 / 84.907], [513.069--562.774], [89,639],
    [GL01 / 3d\_erased], [75.325 / 75.626], [84.830 / 85.143], [342.286--363.442], [24,789],
    [GL01 / 4d\_direct], [59.038 / 58.933], [68.785 / 68.132], [356.819--385.733], [20,029],
  )]
  , kind: table
  )

Memory uses decimal MB and observed process high-water marks; sampled RSS and the graph-reported peak are separate fields in the receipt. The localized route also retains a summed-function-map program (26,041 instructions for GL00; 24,987 for GL01). Retained programs are not additive per-sample operation counts. All operation counts agree across the three independently generated states and with the earlier numerical source.

The final matched cohorts contain 108 accepted passes, 2,160 batches and 54,219,459 samples over 4,364.548388 actual timed seconds. GL00 contributes 22,513,882 samples / 1,839.368333 seconds with pass lengths 15.112271--55.182439 seconds; GL01 contributes 31,705,577 / 2,525.180054 seconds with pass lengths 32.741366--56.092744 seconds. Five undersized attempts associated with the retained GL00/scaled cohort remain excluded; superseded rounds retain their own rejected attempts separately. All 226 child receipts across the original and confirmation rounds succeeded with nonoverlapping recorded intervals. The final states retain their original hashes. All 36 pointwise comparisons pass the unchanged 1e−9 tolerance; worst relative complex-norm differences are 4.3294e−15 for GL00 and 4.4759e−15 for GL01.

=== Dispatch and mapping profile
<dispatch-and-mapping-profile>
Six separate fresh final-source runs cover every physical graph/route with the same cards and H1/CPE5 settings, minimal structured profiling and the exclusive measurement lock. Their generation times are excluded from unprofiled gates. The #link("local-4d-uv-dispatch-profile.json")[dispatch receipt] certifies complete accounting, successful immutable receipts and cold projection caches.

#figure(
  align(center)[#table(
    columns: (9.68%, 12.9%, 12.9%, 12.9%, 12.9%, 12.9%, 12.9%, 12.9%),
    align: (auto,right,right,right,right,right,right,right,),
    table.header([Graph], [Complete nonroot UV (ms)], [Allocation, selection and cache (ms)], [Fraction], [Winning hard CFF (ms)], [Winning template (ms)], [Losing template (ms)], [Row mapping (ms)],),
    table.hline(),
    [GL00], [1646.428], [65.951], [4.0057%], [20.973], [62.042], [29.474], [22.539],
    [GL01], [1650.696], [68.635], [4.1579%], [21.009], [64.956], [31.486], [22.802],
  )]
  , kind: table
  )

Both fractions pass the strict 10% limit. The denominator sums the three nonroot nodes of each four-node forest: Taylor construction, component projection and outer-CT assembly. Root production CFF, final forest summation and later tensor/evaluator work are excluded. Dispatch includes degree/allocation, certificates, all losing native generation and template construction, raw physical-degree reports, cache work and outer routing selection. Outer selection subtracts only the winner's mandatory native generation, source reconstruction and postprocessing. Nested timers are charged once; losing preparation is never free. Other normalization, source, composition and cache intervals are preserved individually in the receipt.

Each direct run has eight source requests, ten admitted candidates and ten template builds, producing 196 selected hard residue rows (at most 62 per request). Six unique native keys are certified by equal independent lower/upper bounds. Final cache snapshots record 12 CFF hits / 6 misses and 1,652 subtree hits / 602 misses, with no LRU evictions and component-local row results cleared. These counts describe the recorded implementation and workload. Request/build counts do not establish unique source bindings or full sample tuples, and are not correctness expectations for the current tests. The allocator's independent sunset regression compares signed complete contours across distributions; cache and projection regressions compare full coefficients and frozen-domain values.

=== Evaluator build intervals and RAM
<evaluator-build-intervals-and-ram>
These are single fresh diagnostic runs per route, not additional acceptance medians. The literal build column sums every actual Symbolica `.build()` call; localized 3D builds two programs. Native orchestration also contains work outside those calls, including function-map preparation. Numeric conversion and expression preparation intervals are separate in the profile receipt.

#figure(
  align(center)[#table(
    columns: (20%, 26.67%, 26.67%, 26.67%),
    align: (auto,right,right,right,),
    table.header([Graph / route], [Literal builds (s)], [Native orchestration (s)], [Graph time minus literal builds (s)],),
    table.hline(),
    [GL00 / 3d\_local], [3.718931], [7.167359], [8.765949],
    [GL00 / 3d\_erased], [1.112125], [1.354999], [5.285839],
    [GL00 / 4d\_direct], [0.862837], [1.100293], [5.579780],
    [GL01 / 3d\_local], [3.781439], [7.143527], [8.468376],
    [GL01 / 3d\_erased], [1.142256], [1.385334], [5.375652],
    [GL01 / 4d\_direct], [0.875562], [1.110370], [5.786078],
  )]
  , kind: table
  )

Excluding only literal builds in these profiles leaves direct/erased ratios 1.05561 for GL00 and 1.07635 for GL01. These broader preparation costs remain slightly higher in 4D, despite faster UV expression construction.

The following cells give #strong[RSS / cumulative VmHWM in MB (observation lag in ms)];. Both localized build entries are listed in order. RSS is sampled after observing the structured event, not synchronously at the instruction boundary; VmHWM is a process lifetime peak, not this phase's allocation. Every observation and its exact timestamps are retained in the linked profile JSON.

#figure(
  align(center)[#table(
    columns: (20%, 26.67%, 26.67%, 26.67%),
    align: (auto,right,right,right,),
    table.header([Graph / route], [Before Spenso], [After Spenso], [At literal build entry],),
    table.hline(),
    [GL00 / 3d\_local], [47.837 / 47.837 (32.615)], [163.430 / 186.655 (39.549)], [265.400 / 405.090 (16.152); 352.887 / 563.630 (34.450)],
    [GL00 / 3d\_erased], [49.635 / 49.635 (39.057)], [167.993 / 177.394 (51.266)], [167.993 / 177.394 (17.508)],
    [GL00 / 4d\_direct], [51.696 / 51.696 (37.974)], [195.777 / 247.493 (37.485)], [195.777 / 247.493 (9.246)],
    [GL01 / 3d\_local], [46.801 / 46.801 (31.984)], [184.750 / 184.750 (36.600)], [262.967 / 424.108 (15.486); 297.046 / 526.062 (38.398)],
    [GL01 / 3d\_erased], [48.407 / 48.407 (32.997)], [166.367 / 173.548 (7.064)], [166.367 / 173.548 (24.593)],
    [GL01 / 4d\_direct], [52.023 / 52.023 (205.917)], [187.339 / 229.552 (49.998)], [187.339 / 229.552 (22.751)],
  )]
  , kind: table
  )

For comparison, whole-process peaks in these six profiles are 563.630/350.990/388.420 MB for GL00 and 526.062/357.020/398.221 MB for GL01 (localized/erased/direct). They are distinct from the unprofiled generation high-water marks and the pre-evaluator values above. The monitor also retains forest completion, expression-preparation and numeric-conversion events.

=== Source-pinned correctness evidence
<source-pinned-correctness-evidence>
The receipt for numerical source `911f1824` records a passing dedicated 183-test integration selection using frozen binaries with assertions enabled, one worker, no retries and a 30 GB process-tree guard. The unchanged selection contains 167 scalar checks (2732.365 s), two physical GL00/GL01 three-route comparisons (105.018 s), four UV-composition checks (22.294 s), seven cut/threshold checks (18.212 s), two API checks (3.987 s) and one analytic check (0.378 s). Exact selections, executable identities and descriptive route timings are retained in the #link("local-4d-uv-correctness.json")[correctness receipt];. The separately pinned `f938f9c` source passed 630 focused unit tests and the same 183-test selection. Neither cohort certifies the reconstructed stack.

The final scalar run supplies 166 single-run setup/generation observations, with a median direct/erased ratio of 0.833248: 122 are below 1.00 and twelve exceed 1.15. The largest are GL21/base (1.431860), GL21/quadratic `q7` (1.376004) and GL17/base (1.368080). These include evaluator construction and enabled integrated/threshold subtraction. They are descriptive observations, not repeated physical acceptance benchmarks or expression-only timings.

On the earlier source, all 167 scalar tests passed. Its 166 single-run setup/generation observations have a median direct/erased ratio of 0.814800; 131 are below 1.00 and eight exceed 1.15. The largest are GL21 (1.417311) and GL17 (1.368180). These observations include evaluator construction and enabled integrated/threshold subtraction; they are not measurements of expression construction alone or repeated physical acceptance benchmarks. A separate frozen GL21 base diagnostic passed and locates its representative cost. The no-numerator UV/forest orchestration improves from 972.364 ms to 405.324 ms; the graph-owned product numerator instead takes 5,944.078 ms in 4D versus 3,377.569 ms in erased 3D. New hard projection accounts for only 55.666 ms and outer routing selection for 82.152 ms. The largest remaining gaps follow integrated-addback localization and precede the next component or final assembly. The shared final-integrand simplification and forest-aggregation boundary is therefore the next generic optimization target; these logs do not identify one dominant internal function. Numerator tensor preprocessing is also higher (1,005.779 versus 644.301 ms), while its Symbolica stack construction is slightly lower (172.676 versus 176.267 ms). This is a single diagnostic with integrated/threshold terms enabled, not a new physical gate. Exact intervals and limitations accompany the #link("local-4d-uv-dispatch-profile.json")[dispatch profile];.

The scalar selection contains 166 route-comparison cases across 49 graph labels, including quadratic/quartic numerator variants, plus one sampling-scale check. It is not 167 distinct topologies. The physical comparisons use the Compare orchestrator to exercise both forest implementations.

The current physical route tests require finite, nonzero complete values and use a fixed `1e-9` relative complex-norm tolerance at three points. Reported numerical accuracy is diagnostic data; it cannot enlarge the acceptance bound. The scalar suite retains its own existing precision and physics contracts.

The imported change has source base `78395e3ab3ddd8d8f62b2f674d7488484eace197`. Its receipt identities preserve that comparison. Current formatting, locked workspace checks, all-target builds, Clippy and runtime selections require their own reconstructed-stack evidence; no PR or remote CI status is inferred from this document.

The imported CI evidence records three failures among 2,137 tests: `finite_part_ghost_2loop`, `se1l_uv` and `epem_a_bbx_amp_uv`. The first constructs analytically integrated 4D renormalization terms and does not enter the canonical 3D projection, CFF dispatch or evaluator stages. A matched original-base comparison identified the first differing expression before Vakint: analytic spin expansion changed `(p.k)^2` into `p^2 k^2` by reusing a contraction index. The analytic owner now protects certified scalar products through existing aliases while expanding open tensor contractions. A related shared executor defect in odd tensor powers is also corrected.

The corresponding rebuilt focused selection records twenty passing analytic/physical checks, including all three CI failures and the new powered-product/compact-metric regression, with unchanged expectations and tolerances. The preceding broader run passes all Spenso/Idenso checks, including the new five-leaf odd-power coverage. The source-pinned physical timings and dispatch profiles meet their recorded gates. The imported CI receipt records 2,140 Linux passes, 2,139 macOS passes with its platform exclusion, and 173 GitHub Nix integration passes. The separate external Nix CI failed while uploading its cache, independently of the completed test checks. These totals describe that measured source only.

=== GL262 and automatic coverage limits
<gl262-and-automatic-coverage-limits>
In the matched CLI18 diagnostic, GL262 forest construction took 900.274 s in direct 4D and 2,254.471 s in erased 3D, a ratio of 0.3993. Direct tensor preprocessing completed, but the erased run hit its time cap during preprocessing. There is no completed paired measurement of the full pipeline. The later direct run failed inside Symbolica evaluator construction with compilation disabled; no complete GL262 evaluator/runtime result is claimed.

The retained CLI21 direct diagnostic observed 46.050947 GB RSS and 68.617658 GB cumulative VmHWM after Spenso (42.844 ms event lag). Actual Symbolica build entry, following 175.972 seconds of expression preparation, had 46.050968 GB RSS and the same VmHWM (8.255 ms lag). The later observed peak was 120.250585 GB RSS / 120.480551 GB VmHWM. Before Spenso it observed 1.487663 GB RSS / 3.001602 GB VmHWM (17.912 ms lag). The erased diagnostic has no completed post-Spenso boundary. These historical failure observations remain separate from final GL00/GL01 evidence in the #link("local-4d-uv-dispatch-profile.json")[profile receipt];.

The separate nonsymmetric `f(x,y)` import regression is reproduced in the tracked #link("../../tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/README.md")[MRE package];. Its 2,480-byte input imports as `f(y,x)` when the reader registers the argument symbols in reverse order. That confirmed import defect is distinct from the large evaluator panic, which has no confirmed standalone reproducer. The large local artifacts are preserved separately and are not required to read this report.

Manual PySecDec integrations remain outside the automatic selections. Neither the GL00/GL01 local-UV measurements, which disable integrated and threshold terms, nor the interrupted GL262 diagnostics establish a complete large-graph physical acceptance with every subtraction enabled. Current regression tests also cover typed-zero projection after same-frame cancellation and complete values under cache eviction; no benchmark gain is claimed for those correctness changes.
