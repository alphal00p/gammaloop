= Large Figure B1 fixture: bounded performance diagnosis

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.
<large-figure-b1-fixture-bounded-performance-diagnosis>
The massless complete-soft-wood acceptance #strong[failed with a byte-buffer bounds panic during explicit 3D production generation];, after 6326.706 seconds reported by Nextest (about 105.45 minutes). It did not complete the explicit UV profile or reach the projected 4D iteration. This failure is not a timeout or an OOM termination. A prior 15-second CPU sample identified a large-expression selector search as dominant sampled production work, but no backtrace was captured for the final panic. No source/test changes, reruns, generation or new profiling were performed for this final diagnosis.

== Observations and their limits
<observations-and-their-limits>
The target is `slow::paper_figure_b1_massless_bubble_uses_a_complete_soft_wood`. It runs localized 3D, explicit 3D, then projected 4D sequentially. Its localized route retains one independently audited nonzero orientation. Explicit routes clear that filter and retain the complete parent source-orientation inventory, including repeated direction signatures with distinct source maps. Comparing their elapsed times cannot yield an equal-work route speedup. See #link("../../../tests/tests/uv.rs:5527")[fixture ordering];, #link("../../../tests/tests/uv.rs:5216")[explicit filter reset];, and #link("../../../tests/tests/uv.rs:5280")[inventory assertion];.

The localized subtracted/bare profiles completed before the test entered explicit 3D. The final log contains the explicit-mode announcement and the bounds panic, with no completed explicit generation summary or UV profile. Source order places `local_ct_cli(...)` before structural spinney assertions and numerical profiling. The earlier CPU sample narrows its sampled production phase to preparation of the final single-parametric evaluator; it does not establish the final panic's call chain. See final test log and #link("../../../tests/tests/uv.rs:5311")[generation-before-assertions boundary];.

The successful stack sample used `cpu-clock:u`, 49 Hz, 15 seconds, DWARF stacks with an 8192-byte capture, and reports 818 samples with zero lost. Its principal entries are:

#figure(
  align(center)[#table(
    columns: 3,
    align: (auto,right,right,),
    table.header([Sampled call path], [Inclusive share], [Self share],),
    table.hline(),
    [`GenericEvaluator::new_from_builder` → lazy iterator consumption → `parametrize_residue_map_selectors` → `AtomView::contains_symbol`], [74.94%], [62.96% in `contains_symbol`],
    [Rayon work executing sysinfo process/task enumeration], [about 24.45%], [spread across directory walking, path allocation, and process/task collection],
  )]
  , kind: table
  )

These rows describe this sampled user-CPU window across sampled threads. Parent/child inclusive rows overlap; do not add them or reinterpret them as wall-time, whole-run, or kernel-CPU fractions. The sample excludes kernel execution. The earlier independent 15-second self-only sample had 352 samples, zero lost, and 45.74% in `__memmove_avx512_unaligned_erms`; without call stacks it cannot attribute that copying to a particular expression operation. Profiling/post-processing itself adds some overhead. See stack report;, self-only report;, and profiler metadata;.

Saved `/proc` observations show 107,867,784 KiB current RSS/high-water at 09:51:13 UTC (102.871 GiB). At 09:55:09 UTC the high-water had grown to 131,995,848 KiB (125.881 GiB), while current RSS was 70,272,824 KiB (67.017 GiB). Report current RSS and high-water separately: a high-water remains high after memory has been released. These are process observations during a shared multi-route test, not isolated expression sizes, final resource totals, or proof of a memory leak. A later saved observation at 10:28:38 UTC records a 248,299,296 KiB high-water (236.797 GiB) and 148,412,352 KiB current RSS (141.537 GiB) while the test remained active. See resource observations;.

== First established production boundary
<first-established-production-boundary>
+ The amplitude builder supplies no runtime orientation catalog when `explicit_orientation_sum_only` is true: #link("../../../crates/gammalooprs/src/integrands/process/amplitude/mod.rs:238")[amplitude/mod.rs:238];.
+ That dispatch selects `new_explicit_sum_with_timings`, which disables the three optional orientation evaluator variants, replaces `OrientationID(...)` selectors with one, then enters shared `new_with_timings` with empty orientation arrays. The explicit-only flag is set only #strong[after] the shared constructor returns: #link("../../../crates/gammalooprs/src/integrands/process/evaluators.rs:555")[dispatch:555];, #link("../../../crates/gammalooprs/src/integrands/process/evaluators.rs:528")[selector stripping and shared call:528];, #link("../../../crates/gammalooprs/src/integrands/process/evaluators.rs:298")[replacement:298];.
+ Shared construction preprocesses the atoms, then always transfers the resulting scalars into `new_single_parametric`: #link("../../../crates/gammalooprs/src/integrands/process/evaluators.rs:881")[preprocessing:881];, #link("../../../crates/gammalooprs/src/integrands/process/evaluators.rs:975")[single evaluator:975];.
+ The single evaluator lazily maps each atom through `parametrize_residue_map_selectors` and `collect_orientation_if`. The first helper searches the expression for `OrientationID` even on this explicit route: #link("../../../crates/gammalooprs/src/integrands/process/evaluators.rs:333")[single evaluator map:333];, #link("../../../crates/gammalooprs/src/integrands/process/evaluators.rs:284")[containment guard:284];.
+ `GenericEvaluator::new_from_raw_params` consumes that lazy iterator while assembling `exprs`. Consequently the sampled `new_from_builder` frame identifies #strong[expression preparation];, before that expression\'s later Symbolica `.evaluator(...).build()` call. It is not evidence of slow compiled evaluator execution: #link("../../../crates/gammalooprs/src/integrands/process/evaluators.rs:1497")[preparation:1497];, #link("../../../crates/gammalooprs/src/integrands/process/evaluators.rs:1537")[later Symbolica build:1537];.

The compiled Symbolica dependency is version 2.2.0 from patched git revision `4d0a833eb8e059d1f95bdae5abed2559830b235f`, verified against checkout HEAD, #link("../../../Cargo.lock:4682")[Cargo.lock];, and the historical dev-optim dependency record `target/dev-optim/deps/symbolica-f7a31ebb3c676ab1.d`;. Its `contains_symbol` implementation walks the expression with an explicit stack, entering function arguments, powers, products and sums. An absent symbol requires a full traversal; this implementation has no cross-call presence cache. The captured power-base/exponent and packed-rational-reading frames agree with that traversal. See #link("/home/lcnbr/.cargo/git/checkouts/symbolica-7e92af2066df6214/4d0a833/src/id.rs:1231")[Symbolica source];.

== Bounded interpretation: avoidable work candidates
<bounded-interpretation-avoidable-work-candidates>
#strong[Selector scans.] The explicit path strips the selector and then re-enters a path that searches for it. That is a concrete candidate for unnecessary full-expression work. The sample proves that the search is costly in the captured window; it does not record the search\'s return value or independently certify that later tensor preprocessing cannot reintroduce a selector. Certifying selector absence after preprocessing is the missing boundary before declaring the scan redundant. The source invokes this guard once per mapped atom; neither the sample nor source inspection proves an accidental retry loop, superlinear traversal, or how many atoms/calls have completed during the long run.

A bounded static pass over `preprocess_atom` and its directly referenced machinery found no explicit `OrientationID` or `residue_map_id` constructor in preprocessing, Idenso, or Spenso. Color normalization rewrites color invariants; scalar aliases restore values from the input scalar store; the configured function library provides conjugation. This makes selector reintroduction unlikely and strengthens the optimization candidate, but it is not a certificate for this actual large post-preprocessing atom or every dynamically populated tensor-library entry. See #link("../../../crates/idenso/src/color/casimir.rs:130")[color rewrite dispatch];, #link("../../../crates/spenso/src/network/mod.rs:1755")[input scalar restoration];, and #link("../../../crates/gammalooprs/src/utils/mod.rs:4943")[configured tensor/function libraries];.

After this search, `collect_orientation_if` can perform two more containment searches for `theta` and `IF` before returning an unchanged atom. If both symbols are absent, these are two further full walks; whether they occurred or mattered in this window is not established. If an `IF` or orientation condition is present, its transformations can be required for correct inactive-branch handling. Do not remove them on the strength of a time sample alone. See #link("../../../crates/gammalooprs/src/utils/symbols.rs:275")[conditional collector];.

#strong[Copying and retained large inputs.] The explicit selector-removal helper starts from a borrowed atom and produces a new owned atom; that stripped vector remains the shared constructor\'s borrowed input while preprocessing constructs scalar output atoms. This provides plausible large simultaneous owners to measure, but does not quantify their memory or attribute the earlier `memmove` sample. Existing code already transfers `parsed_atoms` into the final evaluator and skips empty function-map replacement passes. An owned no-selector atom returned through `AtomOrView::into_owned` is moved, not cloned. Therefore the sampled containment guard should not itself be labelled a full-Atom clone. See #link("../../../crates/gammalooprs/src/integrands/process/evaluators.rs:975")[existing ownership transfer];, #link("../../../crates/gammalooprs/src/integrands/process/evaluators.rs:1503")[empty-substitution guard];, and #link("/home/lcnbr/.cargo/git/checkouts/symbolica-7e92af2066df6214/4d0a833/src/atom.rs:2677")[Symbolica owned/view distinction];.

#strong[Memory-monitor overhead.] Production starts a named `generation-ram-monitor` thread and requests memory information for just the current PID, then sleeps 100 ms between refreshes: #link("../../../crates/gammaloop-api/src/state.rs:614")[state.rs:614];. However, pinned sysinfo 0.33.1 first walks `/proc`, recursively enumerates each entry\'s tasks, and only #strong[afterward] filters the collected entries to the requested PID. Its multithread feature uses Rayon for this directory walk. This explains why a seemingly single-process RSS query produces substantial host process/task enumeration in the sample. The work is monitoring rather than soft-CT algebra. Its whole-run cost and its effect on production elapsed time remain unmeasured. See #link("/home/lcnbr/.cargo/registry/src/index.crates.io-1949cf8c6b5b557f/sysinfo-0.33.1/src/unix/linux/process.rs:772")[enumeration-before-filter implementation] and #link("/home/lcnbr/.cargo/registry/src/index.crates.io-1949cf8c6b5b557f/sysinfo-0.33.1/src/unix/linux/process.rs:658")[recursive task collection];.

== Single next comparison
<single-next-comparison>
On the same explicit post-preprocessing atom, record its byte size, the `OrientationID` guard\'s result and elapsed time, and the elapsed time/result identity of the subsequent orientation-conditional collector, preserving all factorization. This resolves the smallest uncertainty: whether the observed long traversal is proving absence on an already selector-free scalar and how much additional unchanged work follows it. Existing semantic generation timing tags provide the surrounding phase boundaries; no graph-numerator expansion is needed. A subsequent optimization comparison must use the same atom and complete source-orientation inventory and preserve exact output identity. Separately measuring or reducing monitor refresh cost would require a matched measurement; the captured 24.45% cannot be advertised as an achievable speedup.

The acceptance outcome is a production-generation failure. The sampled hotspot and final bounds panic are distinct observations; neither establishes a historical regression or speedup. This fixture supplies no completed explicit/projected numerical UV verdict or runtime throughput measurement.

== Recorded workload size and missing counters
<recorded-workload-size-and-missing-counters>
The graph inventory makes the scale change concrete: the ordinary Figure B1 double triangle has one graph, two loops and five internal edges; the massless self-energy insertion adds a loop and three internal edges, producing one graph with three loops and eight internal edges. The generated integration dimensions are correspondingly 6 and 9. See #link("../../../tests/resources/graphs/paper_figure_b1_double_triangle.dot")[ordinary graph];, #link("../../../tests/resources/graphs/paper_figure_b1_double_triangle_soft_ir.dot")[dressed graph];, #link("../../../tests/tests/uv.rs:5572")[ordinary dimension assertion];, and #link("../../../tests/tests/uv.rs:5357")[dressed dimension assertion];.

#figure(
  align(center)[#table(
    columns: 6,
    align: (auto,right,right,right,right,right,),
    table.header([Existing positive profile], [LMBs], [Hard subsets per LMB], [Summed fits], [Retained localized profile labels], [Orientation fit records],),
    table.hline(),
    [Ordinary double triangle], [8], [3], [24], [18], [432],
    [Dressed massless, localized route], [28], [7], [196], [1], [196],
  )]
  , kind: table
  )

These are measured profile inventories, not estimates of expression size or source-map multiplicity. The ordinary explicit/projected profiles each also contain 24 summed fits and intentionally no per-orientation records. The dressed explicit 3D profile was not reached before the panic, and projected 4D was not reached. Sources: ordinary localized profile;, dressed localized profile;.

The dressed localized run completed one evaluator with rounded stage costs of 0.600 s expression construction, 1.85 s Spenso construction and 3.20 s Symbolica construction (about 5.65 s combined), with a 522.41 MiB reported generation peak. Its bare control reports 0.128/0.014/0.003 s and one evaluator. The ordinary subtracted localized/explicit/projected runs likewise each report one evaluator, with combined displayed stage times of approximately 10.12/7.61/14.06 s. These stage sums are rounded and routes share their test process; they do not establish an equal-work generation ratio. Sources: dressed generation log, lines 17 and 19 and ordinary generation log;.

The final dressed log does #strong[not] record a numerical full parent source-map count, multiplicities of equal direction signatures, exact forest-node count, or evaluator input/output atom-byte sizes. Source assertions preserve full source-map multiplicities while the localized evaluator keeps only one audited orientation, and require a central soft child plus a positive-degree containing soft component; they do not print the complete inventories. The ordinary projected log contains component-level row/cache counters, but these are neither the full evaluator size nor directly comparable source-map/forest totals. Process RSS is not a substitute for any of these missing counters. No heavy generation, binary deserialization, new profiling, or expression dumps were performed to fill these gaps.

== Fixed-point benchmark limitation to verify
<fixed-point-benchmark-limitation-to-verify>
The campaign benchmark script uses one handpicked six-coordinate momentum-space point `(0.11, -0.29, 0.37, 0.59, -0.43, 0.71)` for each of the two benchmark fixtures, repeating it for a warmup-calibrated nominal 5-second target in ten batches. It supplies no independent point-selection criterion establishing representative difficulty. Its JSON `displayed_result` comes from the first warmup sample; the benchmark command itself does not compare values across routes or prove that a near-zero result is materially resolved. The read-only saved-state audit below verifies matching ordered physical loop coordinates across routes. Before using route throughput ratios, also verify that returned values are finite and materially nonzero, and configured results agree across routes to their reported accuracy. Check raw-versus-configured value agreement and stability/cache behavior separately. See benchmark plan;, #link("../../../crates/gammaloop-api/src/commands/bench.rs:325")[warmup and first-result output];.

Minimal mode disables the integrand cache, events, selectors and several stability mechanisms, using a single double-precision stability level; it is a distinct runtime policy, not a numerical reference. The repeated point also keeps its evaluator/precision path highly predictable. Report its timing as throughput at that fixed input and stored state, rather than a representative phase-space average, integrated Monte Carlo throughput, or numerical validation of the large fixture that failed generation. See #link("../../../crates/gammaloop-api/src/commands/bench.rs:481")[minimal settings];.

The six prepared benchmark states have now been checked with read-only integrand displays and small graph exports with UV and full numerators disabled. Their generation LMB sets, complete LMB catalogs, default runtime settings and integrand runtime settings match across routes within each fixture. The exported generated graphs confirm the #strong[ordered] physical assignments: nested top self-energy uses `k0 → e2`, `k1 → e3`; child-only top bubble uses `k0 → e6`, `k1 → e8`. The graph structures, including edge directions, are identical across each fixture\'s three routes. This resolves the coordinate-mapping concern for the chosen benchmark point. Evidence: saved-state LMB audit and exact exported edge lines;. Exports and display logs are referenced individually there. These queries did not regenerate or evaluate the integrands.

The extractor\'s mean/SE labels agree with the implementation. Each timed batch contributes its wall time around `evaluate_samples` divided by its actual number of samples. The reported mean gives these batch means equal weight; batch sizes differ by at most one. The SE is their unbiased sample standard deviation divided by the square root of the batch count. `total_samples` excludes the ten warmup samples, and samples/s is the reciprocal of the reported mean, rather than a separately fitted quantity. These arithmetic relationships were independently checked against the parent\'s two-batch smoke JSON; smoke timings are excluded from performance claims. The smoke also returns identical finite, nonzero nested-localized raw/configured values, with imaginary part about `-2.130253309e-8`; it does not establish cross-route value agreement. See #link("../../../crates/gammaloop-api/src/commands/bench.rs:570")[batch timing];, #link("../../../crates/gammaloop-api/src/commands/bench.rs:805")[SE implementation];, raw smoke;, and configured smoke;.

== Final failure and resources
<final-failure-and-resources>
#figure(
  align(center)[#table(
    columns: 2,
    align: (auto,right,),
    table.header([Measurement], [Final value],),
    table.hline(),
    [Nextest test duration], [6326.706 s],
    [Runner wall time], [6327.580893 s],
    [Child user CPU], [6416.208603 s],
    [Child system CPU], [18080.821573 s],
    [Largest-child peak RSS], [329,213,524 KiB (313.962 GiB)],
    [Runner exit code], [100],
  )]
  , kind: table
  )

The runner started at 09:18:38.538 UTC and completed at approximately 11:04:06 UTC. CPU totals and RSS are runner child-resource observations for the complete command, not attribution to a specific symbolic operation or sums of all simultaneous process memory. The final failure precedes the configured two-hour termination limit. See final runner metadata and Nextest JUnit;.

The recorded panic is:

```text
advance out of bounds: the len is 1576152287 but advancing by 5871119583
```

It originates from #link("/home/lcnbr/.cargo/registry/src/index.crates.io-1949cf8c6b5b557f/bytes-1.11.1/src/lib.rs:169")[bytes 1.11.1's panic formatter];. The slice `Buf::advance` guard reports exactly these available/requested values when the requested advance exceeds the slice: #link("/home/lcnbr/.cargo/registry/src/index.crates.io-1949cf8c6b5b557f/bytes-1.11.1/src/buf/buf_impl.rs:2901")[buf\_impl.rs:2901];. The log explicitly says that the backtrace was omitted.

The size arithmetic gives a specific width-boundary clue: `5,871,119,583 - 1,576,152,287 = 4,294,967,296 = 2^32`, and the available length is exactly the low 32 bits of the requested length. Pinned Symbolica stores function/product byte lengths in `u32` fields and writes larger `usize` lengths using unchecked `as u32` casts: #link("/home/lcnbr/.cargo/git/checkouts/symbolica-7e92af2066df6214/4d0a833/src/atom/representation.rs:894")[function body extension];, #link("/home/lcnbr/.cargo/git/checkouts/symbolica-7e92af2066df6214/4d0a833/src/atom/representation.rs:1148")[product extension];, #link("/home/lcnbr/.cargo/git/checkouts/symbolica-7e92af2066df6214/4d0a833/src/atom/representation.rs:1171")[product replacement];. Additions use a packed `u64` byte length instead: #link("/home/lcnbr/.cargo/git/checkouts/symbolica-7e92af2066df6214/4d0a833/src/atom/representation.rs:1308")[addition header];.

The atom list iterator consumes these mixed-width encodings: it advances by a decoded `u32` for functions/products and by a decoded `u64` for additions: #link("/home/lcnbr/.cargo/git/checkouts/symbolica-7e92af2066df6214/4d0a833/src/atom/representation.rs:2109")[list traversal];. The final requested advance is larger than `u32::MAX`, so it cannot be the direct unchanged `u32` argument at that function/product reader. One compatible mechanism is a large child carrying a 64-bit size inside a parent whose 32-bit byte length was truncated. Other 64-bit advance paths also exist, including coefficient decoding: #link("/home/lcnbr/.cargo/git/checkouts/symbolica-7e92af2066df6214/4d0a833/src/atom/coefficient.rs:685")[coefficient traversal];.

This makes overflow/truncation at a 32-bit encoded-length boundary a #strong[strong hypothesis];, not a captured root cause. No failing atom, originating header write, or final stack was retained, so the exact container/caller and how it became inconsistent remain unproven. It is also not established that the earlier sampled selector scan was the eventual panic site. The next useful certificate is the exact failing call stack plus the relevant atom/container's declared and actual lengths, followed by checking the associated narrowing write; no numerator expansion is required. The arithmetic and source locations are preserved in bounds-panic evidence;.
