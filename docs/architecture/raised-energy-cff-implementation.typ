// Generated from raised-energy-cff-implementation.md by typeset-raised-implementation.py.
// Self-contained: no external Typst packages, imports, or content files are required.
// Links preserve the Markdown targets; this is a companion to the historical review.
#set document(title: "Raised-energy CFF review implementation", author: "Codex",
  description: "Review fixes, Rust patterns, outward test contracts, and post-change verification.")
#set page(paper: "a4", margin: (x: 18mm, top: 19mm, bottom: 18mm),
  header: align(right)[#text(font: "DejaVu Sans", size: 8pt, fill: rgb("667085"))[GammaLoop / Implementation companion]],
  footer: context [#line(length: 100%, stroke: 0.4pt + rgb("d0d5dd"))
    #v(2mm)
    #text(font: "DejaVu Sans", size: 8pt, fill: rgb("667085"))[Raised-energy CFF implementation #h(1fr) #counter(page).display("1 / 1", both: true)]])
#set text(font: "Libertinus Serif", size: 10.5pt, fill: rgb("202939"))
#set par(leading: 0.58em, spacing: 0.78em, justify: false)
#set heading(numbering: none)
#show heading.where(level: 1): set text(font: "DejaVu Sans", size: 16pt, fill: rgb("163e59"))
#show heading.where(level: 2): set text(font: "DejaVu Sans", size: 11pt, weight: "bold")
#show link: set text(fill: rgb("146c94"))
#set list(indent: 1em, body-indent: 0.5em, spacing: 0.7em)
#set table(stroke: 0.4pt + rgb("d0d5dd"), inset: 6pt,
  fill: (x, y) => if y == 0 { rgb("e8eff5") } else if calc.odd(y) { rgb("f7f9fb") } else { none })

#text(font: "DejaVu Sans", size: 9pt, weight: "bold", fill: rgb("146c94"))[IMPLEMENTATION COMPANION]
#v(4mm)
#text(font: "DejaVu Sans", size: 27pt, weight: "bold", fill: rgb("163e59"))[Raised-energy CFF]
#v(2mm)
#text(font: "DejaVu Sans", size: 12pt)[Review fixes · Rust patterns · Test contracts]
#v(3mm)
#text(size: 9pt, fill: rgb("667085"))[Companion to the historical review · Verification status in the final section]
#v(4mm)
#block(fill: rgb("f7f9fb"), inset: 8pt)[
  #strong[Historical record.] The counts, APIs and findings below apply to the
  recorded revisions. See the #link("raised-energy-cff-stack-review.md")[stack review]
  and #link("raised-energy-cff-stack-review-validation.md")[fresh validation record]
  for the reconstructed stack. The original body is preserved.
]

This companion records the implementation following the #link("raised-energy-cff-review.md")[#text("historical review")] of #text(font: "DejaVu Sans Mono", size: 0.82em, "main...raised_​energy_​cff_​wip"). The original findings and their reproduced failures remain historical evidence; this document describes the resulting changes. The #link("raised-energy-cff-review-validation.md")[#text("validation appendix")] records the earlier review commands and reproductions. The #link("raised-energy-cff-implementation-validation.md")[#text("implementation validation record")] gives the actual post-change commands, counts and verification scope.

= Changes corresponding to the twelve findings

#block[
#set text(size: 9.4pt)
#set par(leading: 0.48em)
#table(
  columns: (0.66fr, 1.9fr),
  table.header(
    [#strong[Review finding]],
    [#strong[Resulting behavior and verification boundary]],
  ),
  table.cell(breakable: false)[1. Nested inline variable scope],
  table.cell(breakable: false)[Inline commands pass through the same structured preparation boundary as named blocks. An inner #text(font: "DejaVu Sans Mono", size: 0.82em, "run -D level=error") establishes its environment before expanding its body; absent overrides inherit the outer environment. CLI tests execute both cases and inspect the resulting setting.],
  table.cell(breakable: false)[2. Reusable nested blocks at boot],
  table.cell(breakable: false)[Loading reusable blocks checks literal callee names, recursion and depth without demanding invocation-time variables. Actual invocation still prepares the complete command sequence before executing it. Tests load parameterized nested blocks, invoke them with different values, and reject missing or recursive literal callees before merging their definitions.],
  table.cell(breakable: false)[3. Default 3Drep output location],
  table.cell(breakable: false)[The default workspace uses the existing current-directory output policy in both session modes. A CLI regression builds a scalar box, exits without saving state, and checks the exported expression, its latest-expression pointer and absence of a workspace inside the unsaved state. Explicit state-path protections continue to apply.],
  table.cell(breakable: false)[4. Diagnostic numerator arithmetic],
  table.cell(breakable: false)[The existing parser gives exponentiation precedence over unary signs and associates chained powers to the right. Signed parenthesized integer exponents are supported; arithmetic and narrowing overflow return parse errors. Public evaluation cases include #text(font: "DejaVu Sans Mono", size: 0.82em, "-2**2 = -4"), #text(font: "DejaVu Sans Mono", size: 0.82em, "2**3**2 = 512"), integer limits and oversized exponent chains.],
  table.cell(breakable: false)[5. Comparisons accepting large relative errors],
  table.cell(breakable: false)[Nonzero route values use relative agreement without a universal #text(font: "DejaVu Sans Mono", size: 0.82em, "max(scale,​ 1)") or #text(font: "DejaVu Sans Mono", size: 0.82em, "1e-14") escape. Only the existing exact source-zero certificate permits its separately checked absolute bound. The shared f64 comparison scales components before subtraction and uses #text(font: "DejaVu Sans Mono", size: 0.82em, "hypot"), retaining sensitivity for tiny values without overflow for large finite values.],
  table.cell(breakable: false)[6. NaN/infinity passing comparisons],
  table.cell(breakable: false)[Numerical comparisons require finite components and affirmative closeness. Failure accumulation uses the negation of that complete acceptance condition. Regression pairs reject NaN, infinity, a factor-two discrepancy at #text(font: "DejaVu Sans Mono", size: 0.82em, "1e-18"), and large finite mismatches; equal zero, tiny and near-limit finite values remain accepted.],
  table.cell(breakable: false)[7. Tests fixing private construction choices],
  table.cell(breakable: false)[Tests compare complete expressions, mapped residues, evaluations, exported results and selector behavior. Cache sizes, synthetic names, positional allocation and private builder decisions cease to define success. Persistence checks compare semantic functions and computed exports rather than bincode byte layouts. The benchmark is checked by evaluated values and settings restoration after both success and an output error.],
  table.cell(breakable: false)[8. GL04 denominator certificate],
  table.cell(breakable: false)[The certificate constructs the actual production source from the retained owner-5/owner-6 #text(font: "DejaVu Sans Mono", size: 0.82em, "A²B³") denominator multiset. It derives routing, masses and occurrence metadata from the reconstructed source, checks ownership and UV domain membership, and compares the reconstructed denominator product and mapped numerator against the retained Taylor target.],
  table.cell(breakable: false)[9. Unused scalar snapshots],
  table.cell(breakable: false)[The 198 unreferenced scalar snapshots are preserved byte-for-byte under #text(font: "DejaVu Sans Mono", size: 0.82em, "tests/​resources/​historical/​scalar_​3l_​inspects/​"), with a README identifying them as historical artifacts. Route comparisons and analytic/exact oracles remain the active checks; these files are not presented as exercised golden references.],
  table.cell(breakable: false)[10. Energy-bound validation bypass],
  table.cell(breakable: false)[Bounds are normalized before preserved-denominator and rational-component early returns. Invalid edge IDs therefore fail consistently even when no recursive CFF construction is needed. Diagnostics describe the exclusive #text(font: "DejaVu Sans Mono", size: 0.82em, "0..n") range.],
  table.cell(breakable: false)[11. Malformed graph signatures],
  table.cell(breakable: false)[Validation checks internal loop/external signature lengths and external coefficient lengths before accumulating vertex balances. Malformed rows produce edge-ID diagnostics and #text(font: "DejaVu Sans Mono", size: 0.82em, "ok = false") instead of indexing outside a balance vector.],
  table.cell(breakable: false)[12. Arb median display underflow],
  table.cell(breakable: false)[Median values retain their base-10 logarithm through formatting. Scientific notation is assembled from a bounded mantissa and a separate exponent, so an overflow-bin bound below f64 range displays as #text(font: "DejaVu Sans Mono", size: 0.82em, "<1.0e-1000"). A subnormal in-range case also retains a positive mantissa.],
)
]

= Additional defects exposed by the stronger tests

The unchanged depth-three scalar banana regression exposed a sign mismatch when an ordinary terminal CFF is embedded in a larger UV construction. The shared generator used its component frame sign while advertising the core loop-count convention. The embedded ordinary terminal now supplies the source prefactor consistent with that advertised convention, as the standalone route already does. The generalized normalization remains separate.

An independent two-loop contour regression covers joined and disconnected components in both standalone and embedded contexts. The previously failing depth-two/depth-three core test passes without changing its assertions or tolerance.

The strict GL04 comparisons exposed a second defect in mixed integrated/local UV subtraction. The direct 3D Taylor kernel omitted hard vacuum-mass scaling in a finite integrated coefficient when the Taylor subgraph encloses its owner. Broader GL00 and GL04 comparisons also establish why applying that scaling globally is incorrect: an integrated coefficient from a disjoint component must stay fixed while the other component is expanded.

The direct route now tags each connected coefficient's vacuum mass as transient #text(font: "DejaVu Sans Mono", size: 0.82em, "mUV(owner)") before multiplying coefficients. A Taylor operation scales that mass only when its subgraph contains the owner; it leaves disjoint owners fixed and rejects partial overlaps. A later enclosing operation can scale those retained owners. #text(font: "DejaVu Sans Mono", size: 0.82em, "DirectSector:​:​combine") restores physical #text(font: "DejaVu Sans Mono", size: 0.82em, "mUV") on output copies, preserving the stored owners for further forest steps. Normalized localization kernels stay frozen throughout, and the auxiliary OSE deformation mass retains its separate role.

Identical captured inputs, independent 90–330 digit replays and exact source algebra isolate the initial GL04 discrepancy to one coefficient of a one-loop integrated localizer. The default localization scale of 1000 suppresses its visible size; it is not an Arb arithmetic limit. The independent two-loop renormalization oracle passes on the untouched pre-implementation baseline and rejects changing the established four-dimensional mass scaling. The correction therefore stays in the direct route; the four-dimensional Taylor and integrated-coefficient code retain their established mass scaling.

= Rust patterns and KISS

#strong[Reuse the owner of an invariant.] #text(font: "DejaVu Sans Mono", size: 0.82em, "PreparedCommand:​:​prepare_​with_​environment") now owns inline and named command preparation. Clap uses the existing #text(font: "DejaVu Sans Mono", size: 0.82em, "RunVariable:​:​from_​str") implementation directly. Graph-name lookup delegates to graph lookup. These changes remove parallel paths that had drifted in variable expansion and storage selection.

The existing surface interner becomes an owned #text(font: "DejaVu Sans Mono", size: 0.82em, "ExpressionAssembler"). The three CFF builders reuse its surface copying, interning, variant insertion and numerator-label finalization. Moving these operations together removes repeated assembly code without adding a second builder framework. #text(font: "DejaVu Sans Mono", size: 0.82em, "Direct3dApproximation:​:​run") similarly reuses its existing local and integrated execution methods, keeping Taylor transitions and frozen-factor handling in their established owners.

#text(font: "DejaVu Sans Mono", size: 0.82em, "DirectSector") keeps the active CFF and coefficient separate from normalized localization kernels. The direct correction also carries each integrated factor with its existing subgraph owner, using a tuple at the existing execution boundary. A transient one-argument #text(font: "DejaVu Sans Mono", size: 0.82em, "mUV(owner)") application preserves that ownership after symbolic multiplication; existing subgraph encoding and containment checks determine which masses scale. This adds a small internal symbolic convention, with validation at Taylor expansion and removal on output copies. It resolves the disjoint-versus-enclosing distinction without adding a helper, Rust type, mode or symbol, or changing four-dimensional integration.

#strong[Use types for real alternatives.] #text(font: "DejaVu Sans Mono", size: 0.82em, "CommandTemplate(String)") carries the text required to expand a template; a missing optional history string can no longer turn into a synthetic command. It is skipped by Clap. Template equality is ordinary command equality, and the existing history serialization path reads the payload. The late-parsing test uses a real boolean #text(font: "DejaVu Sans Mono", size: 0.82em, "quit") option and verifies its executed exit behavior after save/load.

Private exact-frame and sub-LMB inputs travel as one #text(font: "DejaVu Sans Mono", size: 0.82em, "Option") tuple, and exact projection options belong to the existing #text(font: "DejaVu Sans Mono", size: 0.82em, "Exact") variant. These representations prevent mismatched optional inputs. The generic energy-candidate classification was only used with #text(font: "DejaVu Sans Mono", size: 0.82em, "()"); its validation now lives directly at the existing source-occurrence boundary.

#strong[Keep shared work in one implementation.] The existing placeholder parser now provides one lazy iterator with source ranges. Discovery, detection and expansion share it; existence checks stop at the first placeholder. Expansion preserves unrecognized text through the untouched ranges instead of maintaining a second scanner.

#text(font: "DejaVu Sans Mono", size: 0.82em, "EvaluatorStack:​:​from_​integrand_​with_​timings") selects the existing orientation-local, explicit-sum or deferred-body constructor. Four callers supply their data instead of repeating three branches. Borrowed optional bodies avoid copies, and a supplied orientation catalog remains distinguishable from explicit mode. The invalid combination of deferred bodies and a catalog returns an error. Existing catalog-length checks still apply. This consolidation preserves the counterterm callers' existing #text(font: "DejaVu Sans Mono", size: 0.82em, "unwrap") behavior; it does not create a universal recoverable-error API for evaluator construction.

#strong[Express rational invariants in the types.] The four #text(font: "DejaVu Sans Mono", size: 0.82em, "LinearEnergyExpr") coefficient fields and #text(font: "DejaVu Sans Mono", size: 0.82em, "CFFVariant:​:​prefactor") now use #text(font: "DejaVu Sans Mono", size: 0.82em, "Rational"), including public constructors, arithmetic, fusion, cut handling and adapters. Symbolic energies remain indexed terms; constants are rational numbers. Generated coefficients already obey this restriction. Mapping tests retain their arbitrary symbolic witnesses through existing energy substitutions.

Native rational arithmetic replaces runtime rationality assertions, repeated #text(font: "DejaVu Sans Mono", size: 0.82em, "Atom") conversion and coefficient string parsing. Integer-only surface conversion still rejects fractions and out-of-range integers. #text(font: "DejaVu Sans Mono", size: 0.82em, "Atom:​:​num") marks the symbolic output boundary; numerical evaluation calls #text(font: "DejaVu Sans Mono", size: 0.82em, "Rational:​:​to_​f64") at its existing f64 boundary. Existing scalar constructors and scaling methods share their rational implementations. This removes helper layers and invalid coefficient states without introducing a new abstraction.

Native serde and bincode support replace the custom coefficient serializers. The coefficient storage format changes; no old-format compatibility layer is retained. State-free coefficient, surface and orientation types no longer require a Symbolica state map to decode. The existing archive test now compares complete expressions and energy maps after binary and JSON round-trips, using negative, fractional and larger-than-i64 coefficients instead of checking re-encoded bytes.

These changes reuse existing functions and methods or local test closures. They do not add a helper framework. Explanatory comments are moved or updated with the code they describe.

= Public interfaces affected

#text(font: "DejaVu Sans Mono", size: 0.82em, "GraphValidation") adds two #text(font: "DejaVu Sans Mono", size: 0.82em, "Vec<usize>") fields:

- #text(font: "DejaVu Sans Mono", size: 0.82em, "internal_​signature_​dimension_​violations") identifies internal edges whose loop or external signature lengths disagree with the graph's declared dimensions.
- #text(font: "DejaVu Sans Mono", size: 0.82em, "external_​signature_​dimension_​violations") identifies external edges with a mismatched external coefficient length.

Both contribute to #text(font: "DejaVu Sans Mono", size: 0.82em, "ok"). They also appear in serialized diagnostic output. Rust callers constructing this public struct must provide the fields; no compatibility fallback for older serialized validation objects is added.

#text(font: "DejaVu Sans Mono", size: 0.82em, "NumericalStabilityMedian") retains its three variants, with their payloads now expressed in base-10 logarithms: #text(font: "DejaVu Sans Mono", size: 0.82em, "Underflow { log10_​lower_​bound }"), #text(font: "DejaVu Sans Mono", size: 0.82em, "InRange { log10_​relative_​accuracy }") and #text(font: "DejaVu Sans Mono", size: 0.82em, "Overflow { log10_​upper_​bound }"). Rust callers matching or constructing these fields must use the new names and units. #text(font: "DejaVu Sans Mono", size: 0.82em, "formatted_​relative_​accuracy()") remains the presentation boundary and does not require materializing an unrepresentably small f64 accuracy.

The hidden Rust #text(font: "DejaVu Sans Mono", size: 0.82em, "Commands:​:​CommandTemplate") variant now requires a #text(font: "DejaVu Sans Mono", size: 0.82em, "String") payload and is excluded from command-line parsing. Command-history text remains the persistence boundary; there is no new template serialization format or public sentinel command.

= Test contracts

Symbolic tests expand copies only to compare exact identities; this normalization does not enter the production path that preserves factorized numerators. Persistence tests in the core consume the stored source with the existing production CFF conversion, including energy factors, component conventions and selectors. API persistence tests additionally re-import computed UV node DOT exports and compare their summed numerator function after loading. They do not fix how many internal nodes are used to represent that function.

Cached versus uncached construction and owner relabeling are checked by complete mapped residues. Host selection is checked by the full selector truth table: every allowed host must produce the same complete contribution exactly once. The combined depth-two/depth-three banana test retains complete Arb evaluations and numerator-restoration checks. Completion-shape, private-cut and builder-layout smoke tests are consolidated into complete function comparisons across graph/bound cases, sampling modes and deterministic seeds; independent contour and partial-fraction oracles remain separate.

Identifiers are checked for uniqueness and valid references, and physical cuts retain their compatibility and selected-result checks. Tests no longer require dense threshold IDs, particular extra-batch positions, intermediate cut-candidate counts or generated name prefixes. Public diagnostic structure and physical provenance remain valid contracts where consumers need them.

The benchmark regression compares a fixed prepared x-space point before, during and after minimal benchmarking. It also directs JSON output to an existing directory, producing a write error after evaluation, then checks restored settings and the same evaluated value. This replaces the private settings mirror with observable behavior.

= Post-change verification

All 805 distinct selected tests are verified passing for the #text(font: "DejaVu Sans Mono", size: 0.82em, "Atom") → #text(font: "DejaVu Sans Mono", size: 0.82em, "Rational") migration. The total combines the concurrent production-code runs with the approved fixture's focused rerun; repeated executions are counted once. The suites used four workers and the supplied license, retaining the repository's per-test resource limits. No assertion or numerical tolerance was weakened. Exact commands, results and scope are in the #link("raised-energy-cff-implementation-validation.md")[#text("validation record")].

#block[
#set text(size: 9.4pt)
#set par(leading: 0.48em)
#table(
  columns: (0.66fr, 1.9fr),
  table.header(
    [#strong[Executed check]],
    [#strong[Observed result]],
  ),
  table.cell(breakable: false)[Shared crate, all features],
  table.cell(breakable: false)[All 106 tests passed in 10.579 seconds, including the extended binary/JSON coefficient round-trip.],
  table.cell(breakable: false)[Core CFF, source mapping, UV, evaluator and numerator selection],
  table.cell(breakable: false)[All 276 selected tests are verified: 275 passed the original run; the corrected affine-circulation fixture passed its focused rerun in 0.634 seconds. Coverage includes exact mapping certificates, cache comparisons, depth-three and independent two-loop pole oracles.],
  table.cell(breakable: false)[API library],
  table.cell(breakable: false)[All 375 tests passed in 13.231 seconds.],
  table.cell(breakable: false)[CLI and selected integrations],
  table.cell(breakable: false)[All 48 tests passed in 801.824 seconds. Coverage includes GL00/GL04, repeated masses, 3Drep export/load, benchmark/settings restoration, scalar cross-section routes, UV pointwise comparisons and full NLO acceptance.],
  table.cell(breakable: false)[Formatting, checks and Clippy],
  table.cell(breakable: false)[All passed for the shared crate and affected core/API/integration packages. The shared crate also passed all-target checks with default features disabled. That configuration retains three pre-existing unused-import/helper warnings; the dependency #text(font: "DejaVu Sans Mono", size: 0.82em, "proc-macro-error2") future-compatibility notice remains.],
)
]

The full NLO acceptance test is included in the 48 integration checks; it requires strict Arb agreement among three routes and the published absolute LO/NLO cross-section results. No full workspace test pass is claimed. Earlier validation recorded 807 distinct tests, including two Spenso tests unaffected by this migration; those earlier executions are kept separate in the validation record.

Earlier checks verified schema export from an unsaved state, 205 internal schema references and all 198 byte-preserving snapshot moves. The untouched baseline also passed the independent #text(font: "DejaVu Sans Mono", size: 0.82em, "finite_​part_​ghost_​2loop") oracle in 69.121 seconds, helping reject the discarded four-dimensional rescaling change.

Cache and label tests retain their nonzero/distinct-map witnesses. Reversed-channel tests keep consistent endpoints and signatures, while triangle residues may legitimately cancel. All original affine-fixture assertions remain. Both mappers now receive the same symbolic witness replacement while retaining their own physical replacements. The approved correction passed its focused test, formatting, core all-target check and Clippy. Persistence tests compare recovered expressions and energy maps.
