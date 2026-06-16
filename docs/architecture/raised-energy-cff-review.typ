// Self-contained typesetting of raised-energy-cff-review.md.
// Code links target the reviewed commit; no external Typst packages are needed.
#set document(date: datetime(year: 2026, month: 9, day: 8), title: "Raised-energy CFF review", author: "Codex",
  description: "Functionality, test contracts, Rust idioms, and KISS; completed 2026-09-08.")
#set page(paper: "a4", margin: (x: 18mm, top: 19mm, bottom: 18mm),
  header: align(right)[#text(font: "DejaVu Sans", size: 8pt, fill: rgb("667085"))[GammaLoop / Engineering review]],
  footer: context [#line(length: 100%, stroke: 0.4pt + rgb("d0d5dd"))
    #v(2mm)
    #text(font: "DejaVu Sans", size: 8pt, fill: rgb("667085"))[Completed 8 September 2026 #h(1fr) #counter(page).display("1 / 1", both: true)]])
#set text(font: "Libertinus Serif", size: 10.5pt, fill: rgb("202939"))
#set par(leading: 0.58em, spacing: 0.78em, justify: false)
#set heading(numbering: none)
#show heading.where(level: 1): set text(font: "DejaVu Sans", size: 16pt, fill: rgb("163e59"))
#show heading.where(level: 2): set text(font: "DejaVu Sans", size: 11pt, weight: "bold")
#show link: set text(fill: rgb("146c94"))
#set list(indent: 1em, body-indent: 0.5em, spacing: 0.7em)
#set table(stroke: 0.4pt + rgb("d0d5dd"), inset: 6pt,
  fill: (x, y) => if y == 0 { rgb("e8eff5") } else if calc.odd(y) { rgb("f7f9fb") } else { none })

#text(font: "DejaVu Sans", size: 9pt, weight: "bold", fill: rgb("146c94"))[ENGINEERING REVIEW]
#v(4mm)
#text(font: "DejaVu Sans", size: 27pt, weight: "bold", fill: rgb("163e59"))[Raised-energy CFF]
#v(2mm)
#text(font: "DejaVu Sans", size: 12pt)[Functionality · Test contracts · Rust idioms · KISS]
#v(3mm)
#text(font: "DejaVu Sans Mono", size: 9pt)[main::raised_energy_cff_wip]
#v(2mm)
#text(size: 9pt, fill: rgb("667085"))[Review completed: 8 September 2026 · Full code diff]
#v(4mm)
#block(fill: rgb("f7f9fb"), inset: 8pt)[
  #strong[Historical record.] The counts, APIs and findings below apply to the
  recorded revisions. See the #link("raised-energy-cff-stack-review.md")[stack review]
  and #link("raised-energy-cff-stack-review-validation.md")[fresh validation record]
  for the reconstructed stack. The original body is preserved.
]

Reviewed base #text(font: "DejaVu Sans Mono", size: 0.82em, "395610143") through branch tip #text(font: "DejaVu Sans Mono", size: 0.82em, "91142139e"): 392 changed files, 76,394 added lines and 7,976 deleted lines. This continuation completes the code-diff review, including tests, removed code, CLI/API, shared CFF generation, UV reconstruction, runtime integration, Spenso/Vakint, and build tooling. Every changed file is accounted for: 186 complete changed-hunk reviews, 205 generated-artifact audits and one binary inventory; generated artifacts and the binary state map have separately stated audit limits. Implementation and tests were left unchanged.

The branch has substantial, useful mathematical safeguards, and no incorrect CFF value was reproduced for a valid physics input. Nevertheless, the tests do #strong[not] meet the requested outward-functionality-only standard. Several numerical assertions can miss invalid or relatively large errors. The broader review also reproduced CLI scope/state-output defects and diagnostic-library input problems. The twelve findings below distinguish those concrete issues from design comments. I would address the P2 findings before using this suite as a merge gate.

The Rust is often idiomatic locally. It is only partly KISS overall: typed ownership and completion stages justify complexity, while duplicated builders, divergent command preparation and correlated optional state add avoidable work. Clean Clippy output does not settle either simplicity or test quality.

= Actionable findings

== 1. \[P2\] Nested inline runs expand variables before entering the inner scope.

#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/gammaloop-api/src/commands/run.rs#L268")[#text("run.rs:268")] sends any placeholder-containing inline command through template expansion before recognizing a nested #text(font: "DejaVu Sans Mono", size: 0.82em, "run"). A direct invocation with its own #text(font: "DejaVu Sans Mono", size: 0.82em, "-D level=warn") and a body using #text(font: "DejaVu Sans Mono", size: 0.82em, "$(level)") succeeds; wrapping that same invocation inside #text(font: "DejaVu Sans Mono", size: 0.82em, "run -c") fails with “Missing command-block variable 'level'”. Both cases were run through the real CLI. Source tracing also shows an outer definition can capture the inner placeholder before the inner override applies. Prepare the inner invocation in its own environment, as the ordinary block path already does at line 477. One environment-aware preparation boundary would remove this divergence.

== 2. \[P2\] Loading nested command blocks demands variables intended for invocation.

#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/gammaloop-api/src/session.rs#L238")[#text("session.rs:238")] skips static validation only when a command's own raw text contains a placeholder. For #text(font: "DejaVu Sans Mono", size: 0.82em, "outer = [\"run inner\"]") and a parameterized #text(font: "DejaVu Sans Mono", size: 0.82em, "inner"), it recursively prepares #text(font: "DejaVu Sans Mono", size: 0.82em, "inner") with an empty environment before the top-level #text(font: "DejaVu Sans Mono", size: 0.82em, "run outer -D level=warn") can execute. A real boot-card reproduction exits 1; the equivalent single-block control exits 0. Reusable nested block libraries are affected even before execution. Keep static name/cycle validation, but defer invocation-variable resolution until the inherited environment exists.

== 3. \[P2\] Default 3Drep output violates the explicit state-write policy.

#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/gammaloop-api/src/commands/threedreps/mod.rs#L476")[#text("threedreps/mod.rs:476")] selects #text(font: "DejaVu Sans Mono", size: 0.82em, "<active-state>/threed_workspace") in ordinary mode, and lines 395–410 save JSON and an output pointer by default. I imported a scalar box, ran #text(font: "DejaVu Sans Mono", size: 0.82em, "3Drep build --no-pretty"), and ended with #text(font: "DejaVu Sans Mono", size: 0.82em, "quit -n"); both files were created inside the state directory despite no save command. #link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/CONTRIBUTING.md#L369")[#text("CONTRIBUTING.md:369")] reserves state-directory writes for explicit save/quit-output operations and enabled logs. Use the existing cwd-output policy consistently, or obtain a deliberate repository-policy change. This is a filesystem behavior issue, not a CFF-value discrepancy.

#pagebreak(weak: true)
== 4. \[P2\] The public diagnostic numerator parser silently changes arithmetic.

#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/three-dimensional-reps/src/eval.rs#L756")[#text("eval.rs:756")] bypasses its right binding power for exponentiation; #link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/three-dimensional-reps/src/eval.rs#L777")[#text("eval.rs:777")] binds unary minus more tightly than powers; and #link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/three-dimensional-reps/src/eval.rs#L836")[#text("eval.rs:836")] narrows an exponent to #text(font: "DejaVu Sans Mono", size: 0.82em, "i32") without checking its range. A compiled public #text(font: "DejaVu Sans Mono", size: 0.82em, "evaluate_expression") reproduction on a generated unit graph gives #text(font: "DejaVu Sans Mono", size: 0.82em, "-2**2 = 4"), #text(font: "DejaVu Sans Mono", size: 0.82em, "-(2**2) = -4"), #text(font: "DejaVu Sans Mono", size: 0.82em, "2**3**2 = 64"), and #text(font: "DejaVu Sans Mono", size: 0.82em, "2**4294967296 = 1"). The first and chained-power cases contradict conventional power precedence; the last silently wraps the exponent to zero. Use checked conversion and define or reject chained powers explicitly. Tests should compare evaluated arithmetic at this public boundary. This parser belongs to the shared crate's diagnostic evaluation feature; the production compiled GammaLoop integrand does not use it.

== 5. \[P2\] Route-comparison tolerances can accept order-one relative errors.

#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/tests/tests/test_runs/scalar_3L_cross_section_inspects.rs#L969")[#text("scalar_3L_cross_section_inspects.rs:969")] applies an absolute #text(font: "DejaVu Sans Mono", size: 0.82em, "1e-14") floor whenever both Arb results have norm below one. The alternative precision-scaled relative assertion consequently does not constrain small nonzero results: #text(font: "DejaVu Sans Mono", size: 0.82em, "actual = 1e-18") and #text(font: "DejaVu Sans Mono", size: 0.82em, "reference = 2e-18") pass despite a 50% normalized difference. The comment describes cancellation to zero, but the exception is applied to every case. The common f64 event and total comparisons also use #text(font: "DejaVu Sans Mono", size: 0.82em, "max(norm, 1)") with #text(font: "DejaVu Sans Mono", size: 0.82em, "1e-10"), so they do not restore sensitivity at small scales (#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/tests/tests/test_runs/utils.rs#L457")[#text("utils.rs:457")], #link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/tests/tests/test_runs/utils.rs#L641")[#text("utils.rs:641")]). The same exception occurs in #link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/tests/tests/test_gamma_star_ttx_nlo_acceptance.rs#L315")[#text("test_gamma_star_ttx_nlo_acceptance.rs:315")]. Use a relative criterion for nonzero results at identical inputs. Handle certified zeros separately, with a justified absolute or propagated error bound. Supplying coordinates as f64 does not by itself justify a universal absolute error in the observable.

== 6. \[P2\] Numerical failure predicates can accept NaN and infinity.

#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/tests/tests/test_runs/utils.rs#L635")[#text("utils.rs:635")] only panics when #text(font: "DejaVu Sans Mono", size: 0.82em, "distance > tolerance"). With a NaN total, the distance is NaN and the comparison is false. With an infinite total, both distance and tolerance can be infinite, so the comparison is false again. This can pass when the remaining event data agree, including empty event groups. The explicit finite checks in the scalar test are inside the optional #text(font: "DejaVu Sans Mono", size: 0.82em, "localized_3d_results") branch; they do not protect every caller path. Require finite components before an affirmative #text(font: "DejaVu Sans Mono", size: 0.82em, "distance <= tolerance") assertion. This also affects failure-list aggregation in the independent shared-crate tests (#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/three-dimensional-reps/src/eval.rs#L2066")[#text("eval.rs:2066")], #link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/three-dimensional-reps/src/generation.rs#L6818")[#text("generation.rs:6818")]) and Arb raised-LU comparisons (#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/gammalooprs/src/cff/mod.rs#L4418")[#text("cff/mod.rs:4418")], repeated at 4510, 4626, 4692 and 4791). A compiled reproduction confirms that the actual #text(font: "DejaVu Sans Mono", size: 0.82em, "F<ArbPrec>") NaN compares neither greater than nor less than or equal to a finite tolerance. An empty failure list therefore does not prove valid agreement. No current production NaN at those sample points is claimed. Using #text(font: "DejaVu Sans Mono", size: 0.82em, "hypot") would also avoid avoidable overflow in the norm calculation.

#pagebreak(weak: true)
== 7. \[P2\] Multiple tests freeze private construction choices.

The following would fail under valid implementation simplifications:

#block[
#set text(size: 9.4pt)
#set par(leading: 0.48em)
#table(
  columns: (0.8fr, 1.3fr, 1.45fr),
  table.header(
    [#strong[Location]],
    [#strong[Incidental commitment]],
    [#strong[Outward replacement]],
  ),
  table.cell(breakable: false)[#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/three-dimensional-reps/src/generation.rs#L10639")[#text("generation.rs:10639")]],
  table.cell(breakable: false)[A private component-product builder must return #text(font: "DejaVu Sans Mono", size: 0.82em, "None").],
  table.cell(breakable: false)[Generate through the public entry point and compare the complete residue against an independent reference.],
  table.cell(breakable: false)[#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/three-dimensional-reps/src/generation.rs#L10912")[#text("generation.rs:10912")]],
  table.cell(breakable: false)[Exactly one variant per orientation and a branching denominator tree.],
  table.cell(breakable: false)[Compare the represented denominator/residue function after summing its representation.],
  table.cell(breakable: false)[#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/three-dimensional-reps/src/generation.rs#L11503")[#text("generation.rs:11503")]],
  table.cell(breakable: false)[The exact origin string #text(font: "DejaVu Sans Mono", size: 0.82em, "bounded_degree_known_factor_cff").],
  table.cell(breakable: false)[Check the high-power reconstruction identity; an internal strategy name is not a physics result.],
  table.cell(breakable: false)[#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/gammalooprs/src/uv/approx/local_3d/tests.rs#L625")[#text("local_3d/tests.rs:625")]],
  table.cell(breakable: false)[A completed independent source sum is hosted at exactly #text(font: "DejaVu Sans Mono", size: 0.82em, "OrientationID(0)").],
  table.cell(breakable: false)[Verify compatible mapping, exactly-once selection, and unchanged complete sum under a different valid host.],
  table.cell(breakable: false)[#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/gammalooprs/src/graph/three_d_source.rs#L6058")[#text("three_d_source.rs:6058")]],
  table.cell(breakable: false)[Exact edge-count sequence and cache population of two.],
  table.cell(breakable: false)[Compare cached/uncached values and sufficient capacity. Put an explicitly required cache resource bound in a separate performance check.],
  table.cell(breakable: false)[#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/gammalooprs/src/uv/approx/mod.rs#L2066")[#text("approx/mod.rs:2066")]],
  table.cell(breakable: false)[Exactly one sector/frame and particular loop-carrier IDs.],
  table.cell(breakable: false)[Check coordinate compatibility and complete reconstructed value across valid charts.],
  table.cell(breakable: false)[#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/gammalooprs/src/graph/three_d_source.rs#L3393")[#text("three_d_source.rs:3393")]],
  table.cell(breakable: false)[Synthetic node-name prefixes; nearby tests fix complete parsed graphs and canonical occurrence maps.],
  table.cell(breakable: false)[Check physical ownership, valid references and complete mapped values under relabeling.],
  table.cell(breakable: false)[#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/gammalooprs/src/cff/mod.rs#L3122")[#text("cff/mod.rs:3122")]],
  table.cell(breakable: false)[Zip two orientation lists and require equal individual maps and carriers.],
  table.cell(breakable: false)[Compare complete owner-invariant sums; separately test any promised public selector semantics.],
  table.cell(breakable: false)[#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/gammaloop-api/src/commands/bench.rs#L1041")[#text("bench.rs:1041")]],
  table.cell(breakable: false)[Field-by-field settings changes, including a private sentinel.],
  table.cell(breakable: false)[Execute benchmarking and subsequent evaluation; verify values and restoration on success and error.],
  table.cell(breakable: false)[#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/gammaloop-api/src/state.rs#L3618")[#text("state.rs:3618")], #link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/gammalooprs/src/processes/process.rs#L1608")[#text("process.rs:1608")]],
  table.cell(breakable: false)[Bincode equality of entire regenerated CFF, graph and cut-group structures.],
  table.cell(breakable: false)[Compare semantic exports, ownership and evaluated save/load results. Byte equality unnecessarily fixes representation.],
  table.cell(breakable: false)[#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/spenso/src/network/tests.rs#L311")[#text("network/tests.rs:311")], #link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/spenso/src/tensors/parametric.rs#L3347")[#text("parametric.rs:3347")]],
  table.cell(breakable: false)[Exact internal tensor storage and allocation choices.],
  table.cell(breakable: false)[Check contraction/evaluation and explicit resource budgets if those are required.],
)
]

This does not mean every structural assertion is wrong. Original numerator ownership, denominator momentum/mass/multiplicity, and preserved factorization are explicit repository requirements. Test those semantic invariants without prescribing an arbitrary tree, cache layout, diagnostic label, or host. A deterministic result need not retain one particular lexicographic tie-break, and public threshold IDs need not be dense unless that is a documented API guarantee. Private-builder #text(font: "DejaVu Sans Mono", size: 0.82em, "None"), diagnostic origin names, host zero and storage layouts should not become indirect correctness requirements. Unit-test placement is also fine: a private module can still be tested by its observable result. Do not replace independent correctness tests with comparisons that call the same implementation twice.

== 8. \[P2\] The GL04 denominator certificate does not inspect reconstructed denominators.

#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/gammalooprs/src/graph/three_d_source.rs#L5208")[#text("three_d_source.rs:5208")] checks two owner-list lengths in a manually constructed mapper, then compares handwritten #text(font: "DejaVu Sans Mono", size: 0.82em, "D(q)^2 D(-q)^3") with handwritten #text(font: "DejaVu Sans Mono", size: 0.82em, "D(q)^5"). It never obtains denominators from the exact source builder. The assertion proves the evenness of that chosen denominator, but cannot detect a generated wrong mass, multiplicity or routing. The preceding numerator identity is useful within its narrower mapper scope. Derive one side from actual reconstructed source records and compare its complete denominator with the retained Taylor target. The earlier report overstated this test's denominator coverage; that claim is corrected here.

== 9. \[P2\] The 198 added scalar snapshots are not exercised.

The new scalar snapshot files under #link("https://github.com/alphal00p/gammaloop/tree/91142139e1fde44cc4709a35ad090631d27ca451/tests/tests/test_runs/snapshots")[#text("snapshots")] have no snapshot assertion or file reader in #link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/tests/tests/test_runs/scalar_3L_cross_section_inspects.rs")[#text("scalar_3L_cross_section_inspects.rs")]. The route loop ends with comparisons between generated routes at line 992. Editing those snapshot values cannot fail these tests. A common error shared by the compared routes can therefore pass even if it changes the recorded values. Restore a small, independently justified outward baseline where appropriate, or explicitly retain these files as historical artifacts rather than calling them regression coverage. These historical values were not treated as current numerical evidence in this review.

== 10. \[P3\] The preserved-tree fast path skips energy-bound validation.

#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/three-dimensional-reps/src/generation.rs#L668")[#text("generation.rs:668")] returns before checking the bound IDs at line 689. A public call with a one-edge, zero-loop #text(font: "DejaVu Sans Mono", size: 0.82em, "ParsedGraph"), preserved edge #text(font: "DejaVu Sans Mono", size: 0.82em, "[0]"), and bounds #text(font: "DejaVu Sans Mono", size: 0.82em, "[(999, 2)]") succeeds and retains those invalid source bounds. I reproduced this against the compiled library. #text(font: "DejaVu Sans Mono", size: 0.82em, "ParsedGraph") uses the default absent edge-index map, so the earlier optional remapping validation does not catch it. Validate all input bound IDs before the early return. This is an input validation inconsistency; no valid-input CFF algebra error is claimed.

== 11. \[P3\] The public graph validator panics on malformed signature dimensions.

#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/three-dimensional-reps/src/validator.rs#L39")[#text("validator.rs:39")] allocates vertex balances from declared loop/external counts, then indexes them using unchecked signature lengths at lines 52–53. A compiled public-API call with an edge signature #text(font: "DejaVu Sans Mono", size: 0.82em, "[1]") and no declared loops panics with a zero-length out-of-bounds index instead of returning an invalid #text(font: "DejaVu Sans Mono", size: 0.82em, "GraphValidation"). Public #text(font: "DejaVu Sans Mono", size: 0.82em, "ParsedGraph") fields allow this construction. Validate dimensions before balance accumulation and return a structured rejection. This concerns malformed external Rust input; no valid internally generated graph failure was demonstrated.

== 12. \[P3\] Arb stability overflow renders an impossible accuracy bound.

#link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/gammalooprs/src/integrands/evaluation.rs#L443")[#text("evaluation.rs:443")] computes #text(font: "DejaVu Sans Mono", size: 0.82em, "10.0_f64.powf(-1000.0)") for the Arb histogram's overflow bound. That underflows to zero. An Arb stability result with estimated relative accuracy zero is placed in this overflow bin, and the displayed median becomes #text(font: "DejaVu Sans Mono", size: 0.82em, "<0.0e0") via lines 500–513. Keep the logarithmic exponent through formatting, or otherwise represent the bound without first converting it to an unrepresentable f64 number. This affects reporting, not the integrand.

#pagebreak(weak: true)
= Functional architecture and its evidence

The important pipeline comparison is between complete results for the same physical cut:

#block[
#set text(size: 9.4pt)
#set par(leading: 0.48em)
#table(
  columns: (0.8fr, 1.3fr, 1.45fr),
  table.header(
    [#strong[Route]],
    [#strong[Ordered stages]],
    [#strong[Correct comparison boundary]],
  ),
  table.cell(breakable: false)[Direct local 3D],
  table.cell(breakable: false)[Parsed physical graph and numerator → bounded CFF → Taylor operations on the complete CFF → keyed localization → evaluator.],
  table.cell(breakable: false)[Each direct key must have its required UV behavior; the complete cut includes all maps and raised-order derivative pieces.],
  table.cell(breakable: false)[Explicit direct 3D],
  table.cell(breakable: false)[Same direct construction → retain each keyed contribution once without runtime selectors → evaluator.],
  table.cell(breakable: false)[Selector-weighted sum must equal the explicit sum.],
  table.cell(breakable: false)[Projected local 4D],
  table.cell(breakable: false)[Retained source owners → 4D Taylor sectors → owner-preserving exact graph and immutable assignment → component CFFs → outer composition → evaluator.],
  table.cell(breakable: false)[First certify numerator and denominator reconstruction; then compare the complete assembled cut with direct 3D.],
)
]

The physical graph, intended numerator and completed cut are shared inputs; local Taylor construction first differs at #text(font: "DejaVu Sans Mono", size: 0.82em, "Direct3dApproximation::run") versus #text(font: "DejaVu Sans Mono", size: 0.82em, "Local4dCts") reconstruction and #text(font: "DejaVu Sans Mono", size: 0.82em, "Projected4dCts") composition. The shared CFF engine and integrated add-back are downstream boundaries, so they should not be blamed for a route difference until their actual inputs have been compared. The review checked ownership/denominator contracts before assessing low-level algebra. It did not establish a new valid-input route mismatch.

Different per-key decompositions or individual #text(font: "DejaVu Sans Mono", size: 0.82em, "lu_cut_order") pieces are not, by themselves, evidence of a physical mismatch. The handoff correctly explains redistribution between derivative slots. Conversely, agreement of just one such piece is insufficient. Tests should reflect this boundary.

The exact-source mapper retains sign-sensitive numerator provenance, and the generation metadata retains component-local prefactor conventions rather than reconstructing them from final algebra. The explicit-sum evaluator rejects individual orientation selection; the runtime group checks distinguish exact residue-map catalogs from physical sign patterns. UV profiling clones its evaluator before fallible work, so errors do not leave the process without its production integrand. It retains graph/cut identity, transforms candidate LMB samples into generation coordinates, and evaluates the summed limit even when individual orientations pass. Versioned state/workspace loading rejects incompatible positional data before decoding. These are useful functional safeguards, not gratuitous abstractions.

Strong existing outward oracles should remain:

- #link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/three-dimensional-reps/src/eval.rs#L3145")[#text("eval.rs:3145")] derives simple/double-pole residues independently and compares full CFF values.
- #link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/three-dimensional-reps/src/eval.rs#L3256")[#text("eval.rs:3256")] compares independent bubble contour moments across seeds and numerator forms. Keep the value assertions while relaxing incidental branch-layout assertions.
- #link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/three-dimensional-reps/src/eval.rs#L1255")[#text("eval.rs:1255")] checks a denominator-cancellation identity under a selected cut.
- #link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/three-dimensional-reps/src/eval.rs#L1730")[#text("eval.rs:1730")] checks invariance under a nonzero auxiliary sampling scale.
- #link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/crates/gammalooprs/src/graph/three_d_source.rs#L4986")[#text("three_d_source.rs:4986")] checks a common-loop numerator identity for a manually constructed mapper. Its denominator section does not certify production reconstruction; see finding 8.
- #link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/tests/tests/test_epem_a_ddx_nlo_acceptance.rs#L30")[#text("test_epem_a_ddx_nlo_acceptance.rs:30")] uses physical normalization and uncertainty requirements, supplying an oracle beyond agreement between implementations.
- #link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/tests/tests/test_cli.rs#L253")[#text("test_cli.rs:253")] exercises command substitution, defaults, nested overrides, validation before partial execution, and persisted command history through outward behavior.

#pagebreak(weak: true)
= Rust patterns and effectiveness

#block[
#set text(size: 9.4pt)
#set par(leading: 0.48em)
#table(
  columns: (0.83fr, 1.17fr),
  table.header(
    [#strong[Pattern and example]],
    [#strong[Assessment]],
  ),
  table.cell(breakable: false)[Sum types: #text(font: "DejaVu Sans Mono", size: 0.82em, "Local3DCts::{Direct, Projected4d}"), #text(font: "DejaVu Sans Mono", size: 0.82em, "Direct3dCts::{Root, Sectors}").],
  table.cell(breakable: false)[Effective. Their variants correspond to real completion boundaries and constrain which operations make sense. Retain these distinctions.],
  table.cell(breakable: false)[Newtypes: #text(font: "DejaVu Sans Mono", size: 0.82em, "Local4dCts"), #text(font: "DejaVu Sans Mono", size: 0.82em, "Full4dCts"), #text(font: "DejaVu Sans Mono", size: 0.82em, "FinalIntegrands"), #text(font: "DejaVu Sans Mono", size: 0.82em, "OrientationID"), #text(font: "DejaVu Sans Mono", size: 0.82em, "CffGlobalPrefactorSign").],
  table.cell(breakable: false)[Mostly effective. Completion stages, index namespaces, and sign parity are meaningful invariants. The sign type's small #text(font: "DejaVu Sans Mono", size: 0.82em, "product")/#text(font: "DejaVu Sans Mono", size: 0.82em, "factor") API is particularly concise.],
  table.cell(breakable: false)[Borrowed graph adapters: #text(font: "DejaVu Sans Mono", size: 0.82em, "ThreeDGraphSource") and #text(font: "DejaVu Sans Mono", size: 0.82em, "GraphThreeDSource<'a>").],
  table.cell(breakable: false)[Effective separation: the shared crate owns CFF generation while GammaLoop owns graph parsing and physical provenance. Borrowing avoids transferring ownership of the full graph.],
  table.cell(breakable: false)[Immutable plan with owned mapping state: #text(font: "DejaVu Sans Mono", size: 0.82em, "EnergyPowerAssignmentPlan"), #text(font: "DejaVu Sans Mono", size: 0.82em, "ExactSourceEnergyMapper"), shared through #text(font: "DejaVu Sans Mono", size: 0.82em, "Arc") in #text(font: "DejaVu Sans Mono", size: 0.82em, "CFFTerm").],
  table.cell(breakable: false)[Effective. One certified assignment controls generation and numerator mapping, and it can survive temporary reconstructed graphs without copying it into every orientation.],
  table.cell(breakable: false)[#text(font: "DejaVu Sans Mono", size: 0.82em, "Result"), #text(font: "DejaVu Sans Mono", size: 0.82em, "thiserror"), iterator #text(font: "DejaVu Sans Mono", size: 0.82em, "collect::<Result<_>>()"), #text(font: "DejaVu Sans Mono", size: 0.82em, "Option::transpose()").],
  table.cell(breakable: false)[Idiomatic and useful at fallible graph/parameter boundaries. Optional integrated localization in #text(font: "DejaVu Sans Mono", size: 0.82em, "direct_3d/forest.rs:280") is a concise example.],
  table.cell(breakable: false)[Deterministic #text(font: "DejaVu Sans Mono", size: 0.82em, "BTreeMap")/#text(font: "DejaVu Sans Mono", size: 0.82em, "BTreeSet") keys and #text(font: "DejaVu Sans Mono", size: 0.82em, "Entry") aggregation.],
  table.cell(breakable: false)[Appropriate for reproducible capacities, ownership, and branch identity. Tests should assert determinism where promised without prescribing the particular deterministic ordering.],
  table.cell(breakable: false)[Explicit encoding of persistent metadata and omission of transient bounds.],
  table.cell(breakable: false)[Justified by the persisted source-frame contract. Custom encode/decode is more maintenance than derive, so semantic save/load evaluation is more valuable than exact encoded-byte equality.],
  table.cell(breakable: false)[Strategy options and small enums.],
  table.cell(breakable: false)[Generally clear, but booleans still permit invalid configuration pairs. Existing validation correctly enforces projected-4D requirements. Keep mode resolution at a single boundary as options grow.],
  table.cell(breakable: false)[Correlated #text(font: "DejaVu Sans Mono", size: 0.82em, "Option") fields in #text(font: "DejaVu Sans Mono", size: 0.82em, "GraphThreeDSource").],
  table.cell(breakable: false)[Less effective: separate frame and coordinate-LMB options require a runtime both-or-neither check at #text(font: "DejaVu Sans Mono", size: 0.82em, "three_d_source.rs:869"). #text(font: "DejaVu Sans Mono", size: 0.82em, "Option<(LoopMomentumBasis, ExactUvSubLmbFrame)>") expresses that particular invariant directly without adding another helper type.],
  table.cell(breakable: false)[Generic candidate equivalence machinery in #text(font: "DejaVu Sans Mono", size: 0.82em, "energy_degree.rs:96").],
  table.cell(breakable: false)[Unnecessary generality. The only #text(font: "DejaVu Sans Mono", size: 0.82em, "try_new<K>") caller supplies #text(font: "DejaVu Sans Mono", size: 0.82em, "K = ()"), so lines 141–164 always create one equivalence class. Direct validation of already-certified occurrence IDs would preserve behavior with less code.],
  table.cell(breakable: false)[Duplicated orchestration in #text(font: "DejaVu Sans Mono", size: 0.82em, "direct_3d/forest.rs").],
  table.cell(breakable: false)[KISS weakness: #text(font: "DejaVu Sans Mono", size: 0.82em, "run") at line 287 and #text(font: "DejaVu Sans Mono", size: 0.82em, "run_local") at line 395 repeat root/sector conversion, coordinate-frame extension, and signed Taylor mapping. Consolidate through the existing owner and preserve its comments/invariants.],
  table.cell(breakable: false)[Free helpers and forwarding layers.],
  table.cell(breakable: false)[Mixed. Small mathematical functions are readable, but #text(font: "DejaVu Sans Mono", size: 0.82em, "contains_placeholder") constructs a whole set merely to test existence, and placeholder scanning is repeated in #text(font: "DejaVu Sans Mono", size: 0.82em, "placeholder_specs")/#text(font: "DejaVu Sans Mono", size: 0.82em, "expand"). Simplify within the existing module when changing this area; avoid adding another abstraction layer.],
  table.cell(breakable: false)[Public compatibility and unused code.],
  table.cell(breakable: false)[The no-op #text(font: "DejaVu Sans Mono", size: 0.82em, "serde = []") feature deliberately preserves downstream Cargo feature vocabulary; deleting it breaks callers selecting that feature. Its maintenance cost is small. The unused #text(font: "DejaVu Sans Mono", size: 0.82em, "default_active_state_output_path"), suppressed with #text(font: "DejaVu Sans Mono", size: 0.82em, "allow(dead_code)"), is a clearer deletion candidate.],
  table.cell(breakable: false)[Term aggregation: #text(font: "DejaVu Sans Mono", size: 0.82em, "CFFTerm") holds expression, orientation and mapper together.],
  table.cell(breakable: false)[Effective improvement over parallel expression/orientation vectors and a truncating #text(font: "DejaVu Sans Mono", size: 0.82em, "zip"); related data now travel as one value.],
  table.cell(breakable: false)[Typed indices with #text(font: "DejaVu Sans Mono", size: 0.82em, "TiVec"), including #text(font: "DejaVu Sans Mono", size: 0.82em, "TopologicalThresholdId") and #text(font: "DejaVu Sans Mono", size: 0.82em, "OrientationID").],
  table.cell(breakable: false)[Effective: topology IDs, expression indices and residue-map keys have different meanings. Non-dense archive tests exercise that distinction through selected values.],
  table.cell(breakable: false)[Deferred functions and shared #text(font: "DejaVu Sans Mono", size: 0.82em, "preprocess_atom") lowering.],
  table.cell(breakable: false)[Effective: avoids expanding large symbolic expressions and reuses tensor lowering. Materialized-value and hyperdual tests are the right evidence.],
  table.cell(breakable: false)[#text(font: "DejaVu Sans Mono", size: 0.82em, "GraphImportSource"), borrowed #text(font: "DejaVu Sans Mono", size: 0.82em, "GraphCatalog"), and #text(font: "DejaVu Sans Mono", size: 0.82em, "GraphImportOptions").],
  table.cell(breakable: false)[Idiomatic sum type, adapter and labeled argument object. #text(font: "DejaVu Sans Mono", size: 0.82em, "graph_name_by_id") can delegate to existing graph lookup instead of duplicating the ownership match.],
  table.cell(breakable: false)[#text(font: "DejaVu Sans Mono", size: 0.82em, "PreparedRun") with inherited ordered environments.],
  table.cell(breakable: false)[Good all-or-nothing preparation model, undermined by separate block, raw-template and inline preparation routes. The two scope bugs show the cost of duplicating the semantic boundary.],
  table.cell(breakable: false)[Hidden #text(font: "DejaVu Sans Mono", size: 0.82em, "CommandTemplate") plus optional raw text.],
  table.cell(breakable: false)[Weak: a sentinel command and separate optional payload allow a template without its text. A private parsed-or-template enum carrying its payload would express the state directly.],
  table.cell(breakable: false)[Repeated shared CFF builder operations.],
  table.cell(breakable: false)[#text(font: "DejaVu Sans Mono", size: 0.82em, "BoundedCffBuilder"), #text(font: "DejaVu Sans Mono", size: 0.82em, "KnownFactorCffBuilder") and #text(font: "DejaVu Sans Mono", size: 0.82em, "LowerSectorCffBuilder") duplicate surface copying, interning and variant insertion at generation.rs:3165, 4637 and 5751. Reuse an existing owner/interner for those operations instead of adding a fourth builder framework.],
  table.cell(breakable: false)[Repeated runtime evaluator strategy dispatch.],
  table.cell(breakable: false)[Amplitude/cross-section terms and their counterterms repeat deferred/explicit-sum/parametric mode selection. Put compatibility and construction in the existing evaluator owner so four sites cannot drift.],
  table.cell(breakable: false)[Rational coefficients stored as unrestricted #text(font: "DejaVu Sans Mono", size: 0.82em, "Atom") or canonical strings.],
  table.cell(breakable: false)[Less effective: surface.rs:69 repeatedly asserts rationality; expression.rs:576 and 593 stringify and reparse rational fusion keys. Keeping the exact value typed would remove invalid states and conversion plumbing.],
  table.cell(breakable: false)[Boxed lazy Cartesian-product iterator in lower-sector generation.],
  table.cell(breakable: false)[Reasonable generation-time tradeoff: limits intermediate storage. Dynamic dispatch alone is not a reason to replace this with a more elaborate generic design.],
  table.cell(breakable: false)[Fallible functions without an error path.],
  table.cell(breakable: false)[#text(font: "DejaVu Sans Mono", size: 0.82em, "KnownLinearExpr::mul_rational"), #text(font: "DejaVu Sans Mono", size: 0.82em, "rational_to_coefficient") and #text(font: "DejaVu Sans Mono", size: 0.82em, "scale_linear_energy_expr_rational") return #text(font: "DejaVu Sans Mono", size: 0.82em, "Result") unnecessarily. Direct values make the real fallible boundaries easier to see.],
  table.cell(breakable: false)[Exact rational basis search.],
  table.cell(breakable: false)[Clear and correct for small loop counts, but candidate row/column combinations can grow combinatorially. No performance conclusion is claimed without measurement.],
)
]

The branch is idiomatic in many local expressions but only partly KISS at the architectural level. The useful complexity is ownership, stage distinction, and exact algebra. The avoidable complexity is duplicated orchestration, impossible optional states, unused generic machinery, and tests requiring incidental layouts. Splitting a large module can help navigation, but it does not itself remove this complexity. Prefer the concrete reductions above before introducing more traits or builders.

The existing #text(font: "DejaVu Sans Mono", size: 0.82em, "series(...).unwrap()") calls in the direct Taylor kernel remain robustness debt: symbolic expansion should ideally propagate contextual errors. Those calls were present in the old local-3D path, so they are not reported as new regressions. No unsafe-code defect was identified in the inspected changes.

#pagebreak(weak: true)
= Suggested test contract

Cover simple/repeated poles and scalar/quadratic/high-power numerators using independent residues; compare complete cuts across all three UV routes; exercise edge relabeling, routing reversal, alternative valid LMBs, equivalent numerator factorizations, and nonzero sampling scales. Verify exactly-once selector coverage and evaluated save/load parity. Require finite results and meaningful relative error bounds, isolating certified-zero cases. Keep resource budgets in explicit performance checks rather than exact internal cache counts.

Use complete emitted source records for reconstruction certificates. Keep non-expansion, physical ownership and factorization tests because those are explicit repository contracts. Relax exact generated names, origin tags, graph layouts, variant counts and private cache populations unless each has a public reason to remain fixed. Tests may live inside a module and still exercise an outward result; moving every test into another crate would not itself improve its oracle.

Exercise nested inline scopes and boot-card loading with invocation variables through the CLI. The existing late-template fixture at #link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/tests/tests/test_cli.rs#L1106")[#text("test_cli.rs:1106")] merely retains a string containing obsolete #text(font: "DejaVu Sans Mono", size: 0.82em, "bench --samples ... -c 1") syntax. Use a valid current command and execute it after substitution. For benchmarks, check subsequent user-visible evaluation after both successful and failed runs, rather than mirroring every temporary settings assignment.

The large CFF test module repeats full-orientation summation, evaluator setup and error calculations, including a roughly 700-line spectator test. Reuse existing comparison/evaluation support and split independent outward cases into named tests. Avoid a new general test framework or helpers that reproduce the production decomposition. No tests were changed during this review.

= Documentation and generated-artifact consistency

The architecture documents preserve useful mathematical contracts, but several “current” statements have drifted. #link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/docs/architecture/architecture-current.md#L514")[#text("architecture-current.md:514")] still names state schema version 4 while the API uses 5; #link("https://github.com/alphal00p/gammaloop/blob/91142139e1fde44cc4709a35ad090631d27ca451/docs/architecture/kysvnqlq-rebase-review.md#L76")[#text("kysvnqlq-rebase-review.md:76")] has readiness/A79 statements about a common-denominator DOD1 topology that conflict with the newer separate natural-topology implementation. The run-history schema lacks the new UV profile selection fields and enum. Reconcile current contracts and regenerate schemas; retain historical reports with a clear historical scope. These are consistency comments, not evidence of a numerical defect.

All 200 changed snapshots were audited as artifacts. The 198 scalar snapshots have the coverage problem in finding 9. The two other snapshots are referenced by assertions. JSON syntax and local schema references were checked; fixture text and source references were inspected. The changed #text(font: "DejaVu Sans Mono", size: 0.82em, "state_map.bin") was inventoried by size/hash, not semantically decoded. Historical acceptance values were not independently rederived from their publications.

#pagebreak(weak: true)
= Review coverage and executed validation

The earlier scoped report did not establish whole-diff coverage. This continuation read all changed code hunks, including removals, and all new code files in full. The very large new shared generation/evaluation files and source mapper were reviewed in full, including their tests. The coverage appendix records every one of the 392 paths, reviewer coverage and the limits of generated/binary audits. A file inventory alone is not treated as source inspection.

Actual validation performed on the reviewed branch:

#block[
#set text(size: 9.4pt)
#set par(leading: 0.48em)
#table(
  columns: (0.83fr, 1.17fr),
  table.header(
    [#strong[Check]],
    [#strong[Observed result]],
  ),
  table.cell(breakable: false)[Formatting; shared all-feature check and Clippy across targets; API check/Clippy across targets, locked dependencies.],
  table.cell(breakable: false)[Passed. These are build/lint checks, not proof of mathematical correctness.],
  table.cell(breakable: false)[Shared representation library, all features.],
  table.cell(breakable: false)[122 tests passed in the earlier licensed run.],
  table.cell(breakable: false)[Earlier focused GammaLoop run plus follow-up.],
  table.cell(breakable: false)[68 tests passed, then 6 additional targeted tests passed. These overlap the broader run below and must not be added to it as distinct tests.],
  table.cell(breakable: false)[Broader CFF, source mapper, LMB, UV approximation and runtime filters.],
  table.cell(breakable: false)[210 tests executed: 206 initially passed; four aborted on test-thread stack overflow. All four passed when rerun with #text(font: "DejaVu Sans Mono", size: 0.82em, "RUST_MIN_STACK=67108864"). No test edits.],
  table.cell(breakable: false)[CLI build and command parsing/template/bench/profile/completion tests.],
  table.cell(breakable: false)[CLI built successfully; 32 targeted API tests passed.],
  table.cell(breakable: false)[Targeted outward integration tests.],
  table.cell(breakable: false)[7 tests passed: 3Drep CLI output/read-only behavior, repeated masses, raised-cut cancellation, benchmark restoration and stability histogram output.],
  table.cell(breakable: false)[Real CLI reproductions and controls.],
  table.cell(breakable: false)[Two scope failures reproduced with successful direct controls; unsaved 3Drep output reproduced.],
  table.cell(breakable: false)[Compiled public Rust reproductions.],
  table.cell(breakable: false)[Invalid preserved-tree bound accepted; malformed validator input panicked; diagnostic parser arithmetic changed; Arb NaN comparisons were unordered.],
  table.cell(breakable: false)[Python dependency lock consistency.],
  table.cell(breakable: false)[#text(font: "DejaVu Sans Mono", size: 0.82em, "uv lock --check --offline") passed (122 resolved packages).],
  table.cell(breakable: false)[Watchdog CLI scenarios.],
  table.cell(breakable: false)[Successful child, propagated child failure, memory-limit termination and invalid-limit rejection produced expected exits 0, 7, 137 and 2.],
)
]

The broad core filter skips 486 of the 696 discovered core tests, including profile exclusions. The four stack failures are recorded as initial failures, not hidden behind the successful rerun; the repository already uses larger stacks for symbolic workloads. Parallel workers were used with the supplied license. The secret was supplied only to child process environments and is not in these documents or repository files.

The full workspace suite, complete scalar release matrix, NLO numerical campaigns, and Python/standalone end-to-end campaigns were not rerun. Source inspection of those tests is not equivalent to executing them. Exact commands, reproduction inputs, logs and overlapping-test-count caveats are in the validation appendix.

Companion files: #link("raised-energy-cff-review-validation.md")[#text("validation evidence")] and #link("raised-energy-cff-review-coverage.md")[#text("file coverage ledger")].
