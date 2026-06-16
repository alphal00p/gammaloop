// Self-contained typesetting of raised-energy-cff-motivation-audit.md.
// Current implementation, validation results, and coverage limits.
#set document(date: datetime(year: 2026, month: 9, day: 11), title: "Raised-energy CFF motivation audit", author: "Codex",
  description: "Current change motivations, controlled physical and representative benchmarks, reproduced defects, Rust patterns, KISS and coverage limits.")
#set page(paper: "a4", margin: (x: 18mm, top: 19mm, bottom: 18mm),
  header: align(right)[#text(font: "DejaVu Sans", size: 8pt, fill: rgb("667085"))[GammaLoop / Change motivation and benchmarks]],
  footer: context [#line(length: 100%, stroke: 0.4pt + rgb("d0d5dd"))
    #v(2mm)
    #text(font: "DejaVu Sans", size: 8pt, fill: rgb("667085"))[Current motivation audit #h(1fr) #counter(page).display("1 / 1", both: true)]])
#set text(font: "Libertinus Serif", size: 10.5pt, fill: rgb("202939"))
#set par(leading: 0.58em, spacing: 0.78em, justify: false)
#set heading(numbering: none)
#show heading.where(level: 1): set text(font: "DejaVu Sans", size: 16pt, fill: rgb("163e59"))
#show heading.where(level: 2): set text(font: "DejaVu Sans", size: 11pt, weight: "bold")
#show link: set text(fill: rgb("146c94"))
#set table(stroke: 0.4pt + rgb("d0d5dd"), inset: 6pt,
  fill: (x, y) => if y == 0 { rgb("e8eff5") } else if calc.odd(y) { rgb("f7f9fb") } else { none })

#text(font: "DejaVu Sans", size: 9pt, weight: "bold", fill: rgb("146c94"))[ENGINEERING REVIEW]
#v(4mm)
#text(font: "DejaVu Sans", size: 26pt, weight: "bold", fill: rgb("163e59"))[CFF changes: evidence and cost]
#v(2mm)
#text(font: "DejaVu Sans", size: 12pt)[Physical benchmarks · Motivating bugs · Rust patterns · KISS]
#v(3mm)
#text(size: 9pt, fill: rgb("667085"))[11 September 2026 · 73 logical changes across seven commits]
#v(4mm)

= #text("Assessment")

#text("The changes do not all have the same justification. Several repair concrete wrong values, index capture, invalid-input handling or output contracts. Others add requested capabilities or simplify ownership. The measurements support grouped dense/mixed contraction, no-chain collection, direct alias registration and exact base caching on their reached workloads. They also expose cost without demonstrated benefit on the tested inputs: CFF compaction and unproductive assignment contenders deserve reconsideration. Powered scalar-dot freshening has a reproduced correctness requirement at the whole-numerator FORM boundary. No universally best contraction ordering is established.")

#text("The inventory covers ")#strong("73 logical categories, all 579 commit/path records and 505 distinct paths")#text(" in the seven reconstructed commits. A category is a semantic group, not an assertion that each changed hunk has a separate causal reproducer. The complete inventory includes source revisions, concrete inputs, current conclusions and a path-to-owner map. A successful current test alone is not evidence that the old implementation failed.")

#text("Benchmark measurements and the 579-record path map are pinned to ")#text(font: "DejaVu Sans Mono", size: 0.86em, "d55a0f3e9d25e5c64bdb9ef9fb57749ecf680390")#text(" (functional parent ")#text(font: "DejaVu Sans Mono", size: 0.86em, "1a747039")#text("); the original tensor baseline is ")#text(font: "DejaVu Sans Mono", size: 0.86em, "39561014")#text(". Most A/B experiments restore one old method or disable one optimization on those dependencies. They are not whole historical-release comparisons. Current closeout dispositions are identified separately: inspection JSON uses a pair-list adapter, the three compact-parser tests exercise public execution, and powered-dot freshening remains required. Final reconstructed-stack validation is supplied in the ")#link("raised-energy-cff-stack-review.md")[#text("stack review")]#text("; these benchmark counts do not certify it.")

= #text("Inputs and measurement contract")

#text("The requested physical input is an isolated two-loop gluon two-point double triangle: a three-gluon vertex at 0, gluons 0–1 and 0–2, and a closed top-quark loop 1→2→3→1. Vertices use the SM ")#text(font: "DejaVu Sans Mono", size: 0.86em, "V_36")#text(" and three ")#text(font: "DejaVu Sans Mono", size: 0.86em, "V_137")#text(" rules, with the fermion-loop sign and complete propagator/vertex numerators retained. An external Lorentz/color trace closes the tensor. This is one physical diagram, not a complete gauge-invariant amplitude or integrated cross section. The 140,005-byte captured contraction input includes its complete CFF expression; UV and thresholds are disabled for that contraction benchmark.")

#text("The ttH capture is the physical GL38 NLO numerator (the source's GL20 label), with UV and threshold terms. BNL supplies four scalar aliases; nine legacy delta-head names are migrated to the current owner without changing arguments, coefficients or factorization. The factorized rank-four stress input supplies large scalar sums specifically to reach automatic deferred broadcasting. Controlled 4D/8D DD/DS/SD/SS tensors isolate storage effects. The CFF-only double-triangle-derived cases use owner-compatible energy-polynomial witnesses, not the full spin/color numerator or certified Taylor sectors. Quark-bubble Vakint inputs retain actual tensor and color algebra.")

#text("Root contraction, foundation, shared CFF and core variants have one excluded warm-up followed by six balanced forward/reverse rounds. Vakint and Feyngen body timing instead use six interleaved rounds with unequal pair-leading order. Feyngen end-to-end uses two exact forward/reverse rounds after a warm-up, a resource-based design revision declared before its measurements. Those limitations constrain small-difference claims. Rust builds use ")#text(font: "DejaVu Sans Mono", size: 0.86em, "dev-optim")#text(" (GammaLoop optimization level 2, dependency package overrides at 3). Executables, inputs and patches are hash-pinned. Root, foundation, core and Vakint timings use one worker on CPU 8; shared CFF uses CPU 9. Feyngen body timing uses CPU 8. Feyngen end-to-end uses the same CPU affinity 8–11 with either one or four workers; a one-worker process may migrate within that mask. A shared lock serializes the timing families and heavy comparisons; builds were coordinated to finish before measurements. Other users run on the host, an AMD EPYC 9754; there is no exclusive-machine or fixed-frequency guarantee. Small differences with overlapping ranges are descriptive observations.")

#text("Stage times exclude parsing/output when named. Whole-child RSS includes initialization, parsing, contraction and output, and can include a launcher high-water floor. It is not isolated kernel allocation. Binary hashing streams bytes; an invalid preliminary whole-binary-hash RSS experiment is excluded. Process-tree watchdog caps are 8GB for the main workflow families and 15GB for foundation probes; builds have separate caps within the coordinated 30GB aggregate budget. Failed observations remain visible and are not retried. Snapshot updates are disabled. Resource-guard termination is distinct from a returned wrong value.")

= #text("Tensor scheduling")

#text("Cells show median network-execution milliseconds / median whole-child RSS in decimal MB. All six scheduling choices use the same complete input, alias threshold and execution mode. SmallestDegree is the existing strategy selected by the diagnostic switch, even though its underlying CLI argument names the default preset.")

#block[
#set text(size: 8.6pt)
#table(columns: (0.31fr, 0.23fr, 0.23fr, 0.23fr),
table.header([#text("Strategy")], [#text("Double triangle")], [#text("ttH")], [#text("BNL")]),
[#text("intermediate-cost")],
[#text("94.564 / 34.7")],
[#text("393.707 / 123.8")],
[#text("7.810 / 149.2")],
[#text("sparse-atom-aware")],
[#text("80.164 / 37.8")],
[#text("289.702 / 251.3")],
[#text("7.871 / 489.0")],
[#text("atom-aware")],
[#text("80.223 / 37.8")],
[#text("300.307 / 251.0")],
[#text("7.545 / 487.5")],
[#text("result-rank-only")],
[#text("89.591 / 50.4")],
[#text("3032.963 / 2699.8")],
[#text("42.700 / 3830.7")],
[#text("entry-aware")],
[#text("83.464 / 37.8")],
[#text("303.174 / 251.3")],
[#text("7.440 / 490.3")],
[#text("smallest-degree")],
[#text("170.297 / 50.4")],
[#text("8GB guard; no result")],
[#text("8GB guard; no result")],
)
]

#text("Whole diagnostic-process medians in seconds include parsing, complete alias resolution and output dumping. They are not full GammaLoop generation or integration times.")

#block[
#set text(size: 8.6pt)
#table(columns: (0.31fr, 0.23fr, 0.23fr, 0.23fr),
table.header([#text("Strategy")], [#text("Double triangle s")], [#text("ttH s")], [#text("BNL s")]),
[#text("intermediate-cost")],
[#text("0.200")],
[#text("0.675")],
[#text("0.793")],
[#text("sparse-atom-aware")],
[#text("0.191")],
[#text("0.804")],
[#text("1.993")],
[#text("atom-aware")],
[#text("0.207")],
[#text("0.836")],
[#text("2.036")],
[#text("result-rank-only")],
[#text("0.212")],
[#text("7.962")],
[#text("19.526")],
[#text("entry-aware")],
[#text("0.204")],
[#text("0.808")],
[#text("1.977")],
[#text("smallest-degree")],
[#text("0.298")],
[#text("8GB guard; no result")],
[#text("8GB guard; no result")],
)
]

#text("Intermediate-cost has the lowest whole-process median on ttH and BNL despite slower network execution than some alternatives. In BNL, alias resolution alone takes 0.281s with intermediate-cost versus approximately 1.44–1.46s with the atom/entry alternatives and 18.81s with result-rank-only. Choosing from the network timer alone would miss this measured downstream cost.")

#text("Intermediate-cost uses less whole-child memory on ttH and BNL, with some network-execution cost relative to faster presets. These measurements do not separate live intermediate allocation from parsing and output. Result-rank-only is especially costly on ttH and BNL. SmallestDegree crosses the 8GB guard during warm-up on both, at 8.33GB and 8.45GB process-tree RSS; those runs have no valid completion time or full result. Its twelve later observations are deliberately unattempted. This compares the provided choices, not every possible contraction tree, cost-weight configuration or threading strategy.")

#text("The scheduling/kernel matrix records 177 observations: 150 measured successes, 25 successful warm-ups and two guard-terminated warm-ups. All ten distinct successful output hashes map to passed complete-value comparisons at three shared exact algebraic points. One output prints 78 instances of the exactly integral literal 2.00000000000000; an explicit 2.0→2 conversion ledger covers that representation boundary while retaining the unchanged oracle and original dump. Representation differences are not correctness failures. Raw distributions retain minima, maxima and sample standard deviations.")

= #text("Multi-index contraction and scalar aliases")

#text("Dense–dense and mixed multi-index contraction already worked through the generic API. The change extends grouped symbolic accumulation to those storage combinations. It does not add previously missing contraction semantics. The physical g–g trace reaches dense contractions over two, three and four common axes; both versions also pass the independent storage-component oracle.")

#block[
#set text(size: 8.6pt)
#table(columns: (0.31fr, 0.23fr, 0.23fr, 0.23fr),
table.header([#text("Input")], [#text("Old network ms")], [#text("Current network ms")], [#text("Old / current")]),
[#text("gg-double-triangle")],
[#text("140.115")],
[#text("94.564")],
[#text("1.48×")],
[#text("ttH-GL38")],
[#text("488.379")],
[#text("393.707")],
[#text("1.24×")],
[#text("BNL-alias-pressure")],
[#text("8.136")],
[#text("7.810")],
[#text("1.04×")],
)
]

#text("The strongest physical kernel gains are on gg and ttH. BNL has a smaller approximately 4% median difference; raw paired samples and spread remain in the evidence. These are network-stage comparisons, not equivalent total-workflow speedups.")

#block[
#set text(size: 8.6pt)
#table(columns: (0.31fr, 0.23fr, 0.23fr, 0.23fr),
table.header([#text("Controlled input")], [#text("Old ms")], [#text("Current ms")], [#text("Old / current")]),
[#text("d4-density1-DD")],
[#text("0.4494")],
[#text("0.3138")],
[#text("1.43×")],
[#text("d4-density1-DS")],
[#text("0.3954")],
[#text("0.2860")],
[#text("1.38×")],
[#text("d4-density1-SD")],
[#text("0.3803")],
[#text("0.2903")],
[#text("1.31×")],
[#text("d4-density1-SS")],
[#text("0.2843")],
[#text("0.2739")],
[#text("1.04×")],
[#text("d4-density16-DD")],
[#text("0.2328")],
[#text("0.0100")],
[#text("23.35×")],
[#text("d4-density16-DS")],
[#text("0.0252")],
[#text("0.0130")],
[#text("1.94×")],
[#text("d4-density16-SD")],
[#text("0.0265")],
[#text("0.0133")],
[#text("2.00×")],
[#text("d4-density16-SS")],
[#text("0.0135")],
[#text("0.0134")],
[#text("1.01×")],
[#text("d8-density1-DD")],
[#text("8.1510")],
[#text("4.1105")],
[#text("1.98×")],
[#text("d8-density1-DS")],
[#text("8.1740")],
[#text("4.2139")],
[#text("1.94×")],
[#text("d8-density1-SD")],
[#text("8.2313")],
[#text("4.1515")],
[#text("1.98×")],
[#text("d8-density1-SS")],
[#text("4.1439")],
[#text("4.1941")],
[#text("0.99×")],
[#text("d8-density16-DD")],
[#text("3.9245")],
[#text("0.0870")],
[#text("45.13×")],
[#text("d8-density16-DS")],
[#text("0.3519")],
[#text("0.1266")],
[#text("2.78×")],
[#text("d8-density16-SD")],
[#text("0.3997")],
[#text("0.1526")],
[#text("2.62×")],
[#text("d8-density16-SS")],
[#text("0.1959")],
[#text("0.1372")],
[#text("1.43×")],
)
]

#text("Each controlled tensor has two common Euclidean axes and one free axis per operand. Density16 means one populated coordinate in sixteen, including deliberately zero-filled dense storage. Every output component is compared with an independent coordinate sum using expanded subtraction on finite components. This does not distribute a graph numerator. The 32 old/current storage smokes and all 224 timed/warm-up storage observations agree. Sparse–sparse is mostly unchanged; neither the old sparse fiber nor the new loop is a nonzero-only traversal. These dimensions/densities do not establish asymptotic performance for enormous sparse domains.")

#text("On BNL, direct registration changes the alias-registration stage from ")#strong("4.669 ms to 0.079 ms")#text(" (58.8×), with four aliases and equivalent complete results. That stage is a small part of total runtime; the total-wall distributions do not establish an end-to-end speedup. The double triangle and ttH create no aliases, so their old-registration flag timings are inactive controls.")

= #text("Chain collection, factorization and powered dots")

#block[
#set text(size: 8.6pt)
#table(columns: (0.31fr, 0.23fr, 0.23fr, 0.23fr),
table.header([#text("Operation/input")], [#text("Old µs")], [#text("Current µs")], [#text("Old / current")]),
[#text("chain_physical")],
[#text("14483.235")],
[#text("83.504")],
[#text("173.44×")],
[#text("chain_scalar")],
[#text("31.408")],
[#text("0.239")],
[#text("131.19×")],
[#text("prune_nozero")],
[#text("7430.241")],
[#text("360.509")],
[#text("20.61×")],
[#text("gl06")],
[#text("195.927")],
[#text("229.500")],
[#text("0.85×")],
[#text("nested")],
[#text("54.603")],
[#text("115.079")],
[#text("0.47×")],
[#text("square")],
[#text("39.461")],
[#text("55.606")],
[#text("0.71×")],
)
]

#text("The no-chain early return removes work while preserving the complete physical numerator. The no-zero result times full canonicalization under a scoped multi-file ablation; it cannot apportion the improvement solely to pruning. Separately, the old zero-pruning path distributes an unrelated spectator, violating the explicit factorization contract even though its expanded algebra is equal.")

#text("Powered scalar-dot freshening is required by the whole-numerator FORM route. For ")#text(font: "DejaVu Sans Mono", size: 0.86em, "(k·p + k·q)^2")#text(", the correct vacuum average is ")#text(font: "DejaVu Sans Mono", size: 0.86em, "k^2*(p^2 + 2*p·q + q^2)/D")#text(", with ")#text(font: "DejaVu Sans Mono", size: 0.86em, "D=4-2*eps")#text(". Disabling only the final freshening pass gives ")#text(font: "DejaVu Sans Mono", size: 0.86em, "k^2*(p^2+q^2) + 2*k^2*p·q/D")#text(": both diagonal terms lose ")#text(font: "DejaVu Sans Mono", size: 0.86em, "1/D")#text(". A four-case public-API comparison gives correct complete values for both current modes and unfreshened projected mode; unfreshened whole-numerator mode is wrong. The ")#link("raised-energy-cff-powered-dot-freshening.md")[#text("powered-dot follow-up")]#text(" retains the actual backend inputs and independent exact residual checks.")

#text("The square, nested and GL06 Rust round trips remain valid controls. They do not establish backend equivalence: the Rust Lorentz parser contracts a scalar base before its exponent, whereas FORM receives indexed copies. Freshening costs more on the measured conversions, but this counterexample justifies retaining it for correctness. No manual PySecDec integration is claimed.")

= #text("Deferred scalar-weight broadcasting")

#block[
#set text(size: 8.6pt)
#table(columns: (0.31fr, 0.23fr, 0.23fr, 0.23fr),
table.header([#text("Variant")], [#text("Network ms")], [#text("Whole wall s")], [#text("Child RSS MB")]),
[#text("deferred")],
[#text("2.943")],
[#text("3.337")],
[#text("18.9")],
[#text("old-eager")],
[#text("465.080")],
[#text("4.409")],
[#text("550.4")],
[#text("aliases")],
[#text("1.560")],
[#text("3.758")],
[#text("189.3")],
[#text("old-alias-registration")],
[#text("1.566")],
[#text("3.676")],
[#text("189.8")],
)
]

#text("The controlled factorized input is ")#text(font: "DejaVu Sans Mono", size: 0.86em, "(c0 q0⊗q0⊗q0⊗q0 + c1 q1⊗q1⊗q1⊗q1) : q2⊗q2⊗q2⊗q2")#text(", where each coefficient sums 12,000 simple rational terms. The independent value is ")#text(font: "DejaVu Sans Mono", size: 0.86em, "c0(q0·q2)^4 + c1(q1·q2)^4")#text(". Python Fraction arithmetic checks three complete values, including full alias definitions. It is a representative scalar-broadcast stress input, not a physical UV-locality certificate or an identity proof for all scalar parameters. Current versus old-eager with aliases disabled isolates automatic deferral; enabling aliases is a separate intervention. None of the ordinary physical captures reaches automatic deferral.")

= #text("Shared CFF generation")

#text("Generation medians in milliseconds; complete output serialization is outside the timer. The complete distributions and all five variants are retained in the evidence bundle.")

#block[
#set text(size: 8.6pt)
#table(columns: (0.26fr, 0.148fr, 0.148fr, 0.148fr, 0.148fr, 0.148fr),
table.header([#text("Input")], [#text("Current")], [#text("Eager layers")], [#text("Eager fanout")], [#text("No base cache")], [#text("No compaction")]),
[#text("gg affine")],
[#text("0.993")],
[#text("1.001")],
[#text("1.002")],
[#text("1.003")],
[#text("1.002")],
[#text("gg dotted")],
[#text("6.831")],
[#text("6.926")],
[#text("6.897")],
[#text("7.024")],
[#text("6.996")],
[#text("cubic box")],
[#text("53.243")],
[#text("52.968")],
[#text("54.906")],
[#text("109.092")],
[#text("60.488")],
[#text("quartic sunrise")],
[#text("20.668")],
[#text("19.777")],
[#text("19.304")],
[#text("34.098")],
[#text("18.129")],
[#text("two bubbles")],
[#text("1.141")],
[#text("1.118")],
[#text("1.183")],
[#text("1.119")],
[#text("1.085")],
[#text("gg rank-6 stress")],
[#text("1433.976")],
[#text("1298.905")],
[#text("1297.506")],
[#text("2376.280")],
[#text("1156.223")],
)
]

#text("All 210 observations complete (180 measured, 30 warm-ups), and every dumped expression matches its fully evaluated trace for the same variant. Exact base caching helps the reached box, sunrise and larger gg topology. On the larger gg diagnostic, disabling compaction is faster in every paired round; mean generation falls from 1.475s to 1.154s, with the same 22.4MB expression and approximately 44MB peak child RSS. This is a concrete simplification candidate for this path, not proof that compaction is unnecessary for products of multiple components.")

#text("Streaming reaches only one-component products in that larger stress case; the two-component bubble is tiny. Neither eager/streaming comparison demonstrates a resident-memory saving here. Source evidence records a GL16 30GB failure, but these isolated experiments do not reproduce or apportion that resource failure. No isolated Rational-versus-Atom or merged-versus-independent-occurrence benchmark is claimed.")

#text("The complete CFF checks include 285 numerical observations and 15 independent denominator-pinch oracles, with maximum pinch error 1.507e-15. The 28 fresh public contracts pass. Exact old-parser-method replay distinguishes ")#text(font: "DejaVu Sans Mono", size: 0.86em, "-2**2")#text(" (4→−4), ")#text(font: "DejaVu Sans Mono", size: 0.86em, "2**3**2")#text(" (64→512), and an oversized exponent (wrapped value→error). These are small public-API counterexamples, not physical benchmark inputs.")

= #text("Core UV caches and assignment scoring")

#text("The core experiments keep physical numerator owners, request-local capacities, exact reconstruction certificates and selected payload/assignment association intact. Total-cache disabling removes both payload and count reuse; count-only disabling retains completed winners. The one-proposal variant scores the first certified rank-ordered candidate but still constructs the proposal list. It is not a complete recreation of the old planner or union-envelope cache.")

#block[
#set text(size: 8.6pt)
#table(columns: (0.31fr, 0.23fr, 0.23fr, 0.23fr),
table.header([#text("Input / variant")], [#text("Wall mean ± SD s")], [#text("Wall median s")], [#text("Peak child MB")]),
[#text("GL04 / current")],
[#text("15.601 ± 0.291")],
[#text("15.622")],
[#text("55.0")],
[#text("GL04 / no-count-memo")],
[#text("15.801 ± 0.392")],
[#text("15.643")],
[#text("54.3")],
[#text("GL04 / no-exact-cache")],
[#text("15.729 ± 0.576")],
[#text("15.552")],
[#text("54.5")],
[#text("GL04 / one-proposal")],
[#text("15.438 ± 0.394")],
[#text("15.416")],
[#text("54.4")],
[#text("GL04-temporal-square / current")],
[#text("5.746 ± 0.285")],
[#text("5.802")],
[#text("44.3")],
[#text("GL04-temporal-square / no-count-memo")],
[#text("5.766 ± 0.300")],
[#text("5.793")],
[#text("44.3")],
[#text("GL04-temporal-square / no-exact-cache")],
[#text("5.778 ± 0.313")],
[#text("5.829")],
[#text("44.3")],
[#text("GL04-temporal-square / one-proposal")],
[#text("5.615 ± 0.515")],
[#text("5.746")],
[#text("44.4")],
)
]

#text("All 56 observations pass (48 measured, 8 warm-ups). The measured unit includes fresh import/generation, three inspections and computed forest export, including work repeated by that exporter. Mean differences of 0.3–2.3% are small relative to the 0.28–0.58s sample spread; these runs do not establish a convincing total-workflow speed or memory benefit.")

#text("Both GL04 inputs have exact equality of all three complete inspection JSON values and all nine computed forest/node DOT exports across variants. Their small real values are compared directly, without a unit-sized absolute tolerance floor. All alternative proposals lose: eleven GL04 groups have map counts ")#text(font: "DejaVu Sans Mono", size: 0.86em, "[14,20,20]")#text("; seven temporal-square groups have ")#text(font: "DejaVu Sans Mono", size: 0.86em, "[32,56,44]")#text(". Count-only disabling retains the same hit counts as current code, so these inputs do not demonstrate losing-count reuse. They exercise scoring overhead but do not justify its benefit. A partial GL16 trace does reach a winning third proposal: ")#text(font: "DejaVu Sans Mono", size: 0.86em, "[1200,1224,1104]")#text(" reduces the rank-first map count by 8%. All four GL16 variants time out at 180s before inspection or computed exports; none supplies full-value parity or a total-runtime comparison. The traced peak process-tree RSS is 545MB, so this is a time limit, not a reproduction of the source's 30GB failure. Matching diagnostic prefixes show far fewer generated contenders with memoization, and 52 paired parsed-source/bounds records match exactly. Cache hits omit the source and some canonical-key fields, so the traces are not an independent proof of every hit's key identity.")

= #text("Vakint and numerator matching")

#text("Both Vakint input modes and both forest owners produce the same complete nonzero Laurent expression for the physical gluon–top-quark bubble. Exact dummy-name normalization and rational cancellation check pole and finite terms; mode flags are verified at the actual backend boundary. The scalar self-energy is a control. All 56 measured/warm-up workflows complete and their full computed exports match the validated results. Scalar medians are 0.783–0.809s; quark-bubble medians are 1.382–1.405s, with 35MB/69MB peak child RSS. Paired monolithic/projected ratios are about 0.990 and 1.003 for the two owners, with overlapping ranges: no consistent speed advantage is established. Projected mode leads three of six paired rounds for quark/hedge and four of six for the other mode pairs; the schedule is interleaved, not exactly balanced.")

#text("For the full two-loop double triangle, projected processing exceeds 180 seconds for both owners, while whole-numerator processing reaches the 8GB guard for both. These are explicit coverage limits, not failed numerical comparisons. Both modes still require d-dimensional closure and adequate Laurent depth; choosing an input mode does not replace those correctness requirements.")

#text("The numerator-matching input is a two-loop three-gluon vertex family with top-quark loops: 36 candidates in 24 topologies. Twelve nonsingleton buckets contain two candidates each; all twelve compared pairs are rejected. The inputs carry five numerical and five polynomial samples, with no canonical numerator. Twenty repeated calls use each actual physical pair. This reaches early rejection but neither long bucket scans nor the canonical-numerator path. The smaller physical control additionally exercises a matching pair with complete result ")#text(font: "DejaVu Sans Mono", size: 0.86em, "Some(1)")#text(".")

#text("The body experiment completes all 28 observations (24 measured, four warm-ups). Every pair key is unique; all complete comparison results agree between old and current methods, and every 36-graph DOT map agrees with the validated physical reference. Times below are medians per comparison call, averaged across the twelve physical pairs within each process. The timer includes black-box argument/result handling and a push into a preallocated result vector; setup, allocation and equality assertions are outside it. These microbenchmarks repeat the comparator, not the full generator.")

#block[
#set text(size: 8.6pt)
#table(columns: (0.31fr, 0.23fr, 0.23fr, 0.23fr),
table.header([#text("Comparison")], [#text("Old ms/call")], [#text("Current ms/call")], [#text("Paired old/current ratio [range]")]),
[#text("sign")],
[#text("2.510")],
[#text("0.355")],
[#text("7.074 [6.875, 7.349]")],
[#text("scalar")],
[#text("172.269")],
[#text("171.988")],
[#text("0.999 [0.984, 1.009]")],
)
]

#text("Ratios are computed within each paired round before taking their median; they need not equal the ratio of the two marginal time medians.")

#text("Sign-only rejection is faster on this reached workload; scalar-factor comparison shows no benefit beyond run-to-run variation. The six body rounds are interleaved with current first in four sign rounds and five scalar rounds, so small differences must not be presented as precise causal effects. Repetition-inclusive CLI medians are about 67s for sign and 110s for scalar: most of that work is numerator preparation, and twenty artificial comparator repetitions must not be confused with default generation cost.")

#text("The global-lock intervention retains current bucket storage, clone removal and zero-mask construction. It restores serialization around lookup/search/insertion rather than rebuilding the complete old architecture. The larger four-gluon case reaches 567 candidates but times out during preparation/grouping at 180 seconds; it supplies no completed graph parity or timing ratio. A 5,000-repetition calibration also times out, which is a limit of the artificially amplified body experiment, not a demonstrated timeout of default generation.")

#text("The end-to-end experiment accounts for 48 executed observations out of 48 planned: 48 successful, 32 measured, 0 failed and 0 omitted after a failed variant. Successful observations retain all 36 physical graph exports and agree directly with the validated default. Each variant has one excluded warm-up and two exact forward/reverse measured rounds; the repetition hook is disabled.")

#text("The cells below contain the range of the two complete CLI wall times in seconds. Preparation, lookup and grouping remain included; these are descriptive observations on a shared host.")

#block[
#set text(size: 8.6pt)
#table(columns: (0.32fr, 0.17fr, 0.17fr, 0.17fr, 0.17fr),
table.header([#text("Mode / workers")], [#text("Current s")], [#text("Old comparator s")], [#text("Global lock s")], [#text("Both s")]),
[#text("sign / 1 workers")],
[#text("66.002–66.507")],
[#text("66.206–67.265")],
[#text("66.208–66.614")],
[#text("66.516–66.776")],
[#text("sign / 4 workers")],
[#text("19.048–20.435")],
[#text("18.634–18.810")],
[#text("18.146–18.170")],
[#text("17.737–18.867")],
[#text("scalar / 1 workers")],
[#text("68.633–68.913")],
[#text("68.789–68.800")],
[#text("68.526–69.190")],
[#text("69.236–69.276")],
[#text("scalar / 4 workers")],
[#text("19.641–19.725")],
[#text("18.622–18.761")],
[#text("18.560–18.595")],
[#text("19.005–19.526")],
)
]

#block[
#set text(size: 8.6pt)
#table(columns: (0.31fr, 0.23fr, 0.23fr, 0.23fr),
table.header([#text("Mode / workers")], [#text("Old/current ratio")], [#text("Global/current ratio")], [#text("Both/current ratio")]),
[#text("sign / 1 workers")],
[#text("0.995–1.019")],
[#text("1.002–1.003")],
[#text("1.004–1.008")],
[#text("sign / 4 workers")],
[#text("0.920–0.978")],
[#text("0.888–0.954")],
[#text("0.868–0.991")],
[#text("scalar / 1 workers")],
[#text("0.998–1.002")],
[#text("0.994–1.008")],
[#text("1.005–1.009")],
[#text("scalar / 4 workers")],
[#text("0.948–0.951")],
[#text("0.943–0.945")],
[#text("0.968–0.990")],
)
]

#text("No end-to-end gain is established. In the four-worker observations, current takes longer than each ablation in both rounds; these results must not be hidden by the faster isolated sign-comparison body. The small shared-host sample does not establish the cause or a general regression.")

#text("The numerator-preparation/grouping milestones, whole-child RSS and full samples are retained in the evidence. The milestone is an inclusive preparation/grouping duration, not isolated lock time. Two observations per variant cannot establish a precise speedup; worker scaling mostly reflects the complete parallel preparation workload. With at most two candidates per bucket, this experiment cannot certify the motivating large-bucket/global-lock claim.")

= #text("Reproduced defects and current dispositions")

#text("On the pinned benchmark source, the focused foundation suite has 22/22 passes. Reversing eight production files on current dependencies gives 12 passes and 10 failures: seven index/normalization/namespace/factorization failures and three newly supported compact-syntax probes. These are not ten independent physics bugs. The baseline regression build additionally enables Idenso serialization/reference-case features; the named arithmetic paths are unchanged, but the full build feature sets are not identical. The nested-sum Symbolica panic is a pinned-source counterexample; the fixed dependency remains in this local ablation. No failing product test is edited.")

#strong("Inspection JSON preserves structured additional-weight identifiers.")#text(" ")#text(font: "DejaVu Sans Mono", size: 0.86em, "GenericAdditionalWeightInfo.weights")#text(" uses ")#text(font: "DejaVu Sans Mono", size: 0.86em, "#[serde(with = \"vectorize\", bound(deserialize = \"T: Deserialize<'de>\"))]")#text(", reusing the existing adapter to serialize ")#text(font: "DejaVu Sans Mono", size: 0.86em, "[key, value]")#text(" pairs, including ")#text(font: "DejaVu Sans Mono", size: 0.86em, "[]")#text(" for no weights. The explicit bound restores Serde's original deserialization requirement without strengthening the generic type. The motivating saved-GL04 probe isolates JSON serialization after successful numerical evaluation: the benchmark source rejects a structured threshold-counterterm map key with ")#text(font: "DejaVu Sans Mono", size: 0.86em, "key must be a string")#text(". Current outward tests cover all four key variants and empty weights through JSON/bincode round trips, plus physical inspection with weight retention both enabled and disabled. They compare retained weights, event grouping and evaluation values; the CLI comparison does not claim every event metadata field. Final execution evidence is supplied in the ")#link("raised-energy-cff-stack-review.md")[#text("stack review")]#text(".")

#text("A diagnostic boundary also needs care: reparsing Python Symbolica's canonical string for ")#text(font: "DejaVu Sans Mono", size: 0.86em, "(2+3i)/7*x")#text(" changes coefficient grouping. The plain printer round-trips in the reproducer. This is not evidence of a GammaLoop evaluation defect; the comparison harness retains live expressions and parses original Rust plain-string dumps, never intermediate canonical strings. The exact reproducer is retained.")

#text("The three compact-parser regressions in ")#link("../../crates/spenso/src/network/parsing/test.rs")[#text("parsing/test.rs")]#text(" use public parsing, execution and complete scalar results. The weighted-vector test compares against an independent four-component Minkowski sum with signature ")#text(font: "DejaVu Sans Mono", size: 0.86em, "(+,-,-,-)")#text(" for weights ")#text(font: "DejaVu Sans Mono", size: 0.86em, "-2")#text(", ")#text(font: "DejaVu Sans Mono", size: 0.86em, "a")#text(" and ")#text(font: "DejaVu Sans Mono", size: 0.86em, "a+b")#text(". The two ambiguous arguments—two compact vectors, or a compact vector times an explicit tensor—retain their complete opaque function values. Exact expanded zero residuals check these values without fixing intermediate products, allocated slots or traversal order. Separate external-slot contracts remain. Exact factorization assertions remain appropriate where preservation of factorization is itself the intended public behavior.")

= #text("Rust patterns and KISS")

#block[
#set text(size: 8.6pt)
#table(columns: (0.35fr, 0.65fr),
table.header([#text("Pattern")], [#text("Effectiveness supported by the audit")]),
[#text("Shared ")#text(font: "DejaVu Sans Mono", size: 0.86em, "Rc<Cell<_>>")#text(" allocation plus explicit reservations")],
[#text("Shared allocation is established behavior; reservations fix a reproduced capture. Keep ownership local to one parse and do not claim this introduces sharing.")],
[#text("Borrowed input views and grouped iterator accumulation")],
[#text("Avoid repeated symbolic normalization/copying; dense/mixed measurements support the specialization. General sparse-coordinate traversal remains a limit.")],
[#text("Owned scalar aliases and direct registration")],
[#text("The alias already exists, so registering it avoids a redundant root replacement search. The four-alias BNL stage supports this simplification.")],
[#text("Explicit deferred-result alternatives")],
[#text("Preserve a factorized sum until terminal materialization when eager broadcasting is expensive. The added branch must earn its complexity on an activated workload.")],
[#text("Concrete ")#text(font: "DejaVu Sans Mono", size: 0.86em, "Rational")#text(" coefficients and boundary conversions")],
[#text("Excludes unsupported symbolic coefficient states and makes exact arithmetic intent clear. Previous rational Atoms were already exact; no rounding or speed claim follows from the type alone.")],
[#text("One expression assembler and fallible public conversions")],
[#text("Shared invariant ownership reduces duplicate construction logic. Small malformed-input counterexamples justify Result-based rejection before indexing.")],
[#text("Immutable assignment plans and typed occurrence/owner IDs")],
[#text("Keep capacity, routing, physical owner and selected expression together. Reconstruction certificates justify these distinctions; LMB coordinates do not replace physical ownership.")],
[#text("Exact-key caches and bounded heuristic search")],
[#text("Base caching pays on reached high-rank inputs. No benefit from losing-count memoization or extra contenders is demonstrated on the completed GL04 inputs; a smaller map is only a proxy for total cost.")],
[#text("Mode enums and builder-owned settings")],
[#text("Clarify supported alternatives and keep migration at an owner. Both Vakint modes are requested functionality; no measured universal advantage justifies presenting one as inherently superior.")],
[#text("Short bucket locks and staged exact/sample comparisons")],
[#text("The reached sign-rejection body is faster; scalar matching has no demonstrated gain. The physical family has only two-candidate buckets, so the large-bucket motivation remains unproved. Canonical-numerator matching is not reached, and sampled matches remain probabilistic.")],
[#text("Serde field adapter with an explicit deserialization bound")],
[#text("Pair-list serialization preserves structured identifiers without string parsing. Reusing vectorize keeps the fix at the owning field; its explicit bound retains the original generic deserialization requirement.")],
)
]

#text("Keep the demonstrated correctness fixes and required feature boundaries. Treat combined root/product compaction on the reached single-component path, unproductive proposal scoring and unproved large-domain sparse behavior as review findings with explicit evidence limits. Do not remove mathematically necessary ownership distinctions merely to reduce line count, or retain a performance mechanism only because current tests pass.")

= #text("Per-change motivation inventory")

#text("Each row states the current justification or measurement limit. The machine-readable inventory retains pinned benchmark evidence and the original path ownership, with separate ")#text(font: "DejaVu Sans Mono", size: 0.86em, "closeout_disposition")#text(" fields for the completed JSON/test changes and retained powered-dot requirement. Recorded counterexamples identify a concrete failure in retained source evidence; they are not relabelled as fresh executions. Capability, policy and consumer migrations need not invent a motivating bug.")

= #text("Commit 1: motivation rows")

#block[
#set text(size: 8.6pt)
#table(columns: (0.35fr, 0.65fr),
table.header([#text("Change")], [#text("Evidence and conclusion")]),
[#text("C1-01 — Reserve explicit compact-parser dummy names")],
[#text("Reproduced index capture: a compact contraction collides with an explicit Dummy(1000000). The source also records the physical ttH overcontraction. Reservation is the new behavior; shared Rc allocation already existed.")],
[#text("C1-02 — Reserve external and dual indices during canonical dummy naming")],
[#text("Reproduced external/dual index capture during canonical relabeling. Preserve external index identity; no timing claim is needed.")],
[#text("C1-03 — Ignore canonical names from canceled representation groups")],
[#text("Reproduced canonical-identity failure after one representation cancels. A canceled group must not consume names used by the surviving contraction.")],
[#text("C1-04 — Preserve completed contractions across nested tensor sums in pinned Symbolica")],
[#text("Pinned source records a nested-sum contraction-counter panic and a reduced expression. Current probes pass; the old Symbolica dependency itself was not rebuilt in this audit.")],
[#text("C1-05 — Prune antisymmetric zeros locally without distributing spectators")],
[#text("Reproduced unwanted distribution of a spectator factor. The no-zero fixture also runs faster under current full canonicalization; that timing does not isolate pruning from all naming changes.")],
[#text("C1-06 — Skip chain collection on expressions with no chains")],
[#text("Measured benefit: the captured physical numerator has no remaining chains. Complete output is unchanged, and the early return removes almost all collection work.")],
[#text("C1-07 — Allow one implicit compact axis beside explicit spectator slots")],
[#text("Supported syntax extension: one implicit axis can coexist with explicit spectator slots. Public probes check external indices; this is not evidence that dense/mixed contraction was broken.")],
[#text("C1-08 — Allow scalar-weighted compact vectors while retaining ambiguity checks")],
[#text("Supported syntax extension for scalar-weighted vectors, useful for momentum differences in a 3g vertex. Public parse/execute regressions now compare an independent weighted Minkowski sum and complete opaque values for ambiguous arguments.")],
[#text("C1-09 — Freshen independent copies of powered scalar dots")],
[#text("Reproduced correctness requirement: disabling only freshening breaks the whole-numerator FORM angular average of (k·p+k·q)^2, dropping 1/D from both diagonal terms. Both current modes and unfreshened projected mode match the complete independent oracle. Rust-only round trips miss this backend failure; the measured conversion cost does not justify removal.")],
[#text("C1-10 — Use numeric imaginary coefficients in Vakint normalization")],
[#text("Reproduced wrong symbolic normalization when a user registers the imaginary-name symbol. A numeric complex coefficient preserves the intended value independently of the symbol registry.")],
[#text("C1-11 — Preserve canonical user namespaces across FORM")],
[#text("Reproduced FORM parsing failure with canonical user namespaces. Complete tensor-reduction and adapter probes pass with the fix.")],
[#text("C1-12 — Generalize grouped Atom contraction across tensor storage kinds")],
[#text("Measured optimization. Old dense/mixed contractions already return correct values. Generalized grouping speeds the physical captures and most controlled storage cases; sparse–sparse gains are limited.")],
[#text("C1-13 — Score contraction candidates by intermediate symbolic cost")],
[#text("Measured tradeoff across all five exposed presets plus SmallestDegree. Intermediate-cost has the lowest whole diagnostic time and memory on ttH/BNL, although other presets have faster network steps. It is not a universal optimum. SmallestDegree crosses 8GB on two captures.")],
[#text("C1-14 — Register already-created scalar aliases without recompressing the root")],
[#text("Measured simplification on the four-alias BNL capture: direct registration avoids an unnecessary root search. The physical double triangle and ttH create no aliases and are inactive controls.")],
[#text("C1-15 — Defer large coefficient broadcasting while preserving ordinary terminal results")],
[#text("Measured on a factorized rank-four scalar-weighted tensor that reaches automatic deferral. The three ordinary captures do not reach it. Complete Fraction component oracles check both routes.")],
[#text("C1-16 — Select clang for macOS Rust unwinding")],
[#text("Pinned source records the macOS linker/panic-unwinding problem. Linux execution cannot reproduce that platform failure; no cross-platform speed claim is made.")],
[#text("C1-17 — Guard heavy process trees and keep local runtime outputs out of version control")],
[#text("Resource policy and tooling support. The audit exercises the process-tree guard and retains its terminations. Ignore patterns have no runtime performance justification.")],
[#text("C1-18 — Keep manual PySecDec integrations outside automated selections and require supported nextest")],
[#text("Explicit test-selection policy and supported nextest requirement. Manual PySecDec integrations remain outside this automatic audit; this is not a numerical algorithm change.")],
[#text("C1-19 — Pin compatible Python Symbolica and UFO loader")],
[#text("Pinned source records Python Symbolica import/fixed-point trouble and a UFO loader requirement. Dependency compatibility is the motivation; no independent physics speedup is claimed.")],
[#text("C1-20 — Regenerate shared dependency feature metadata")],
[#text("Required generated dependency metadata. It keeps feature resolution/builds consistent with the owning API changes; it is not an independent optimization.")],
[#text("C1-21 — Remove excess trait bounds and migrate sum-capable consumers")],
[#text("API simplification and necessary consumers of sum-capable terminal results. Removing unused bounds narrows conceptual requirements; no separate wrong-value bug is established.")],
[#text("C1-22 — Preserve physical numerator factors in diagnostic consumers")],
[#text("Required factorization-preserving diagnostics and outward assertions. Tests must compare complete values or an explicit factorization contract, not formatting or intermediate storage.")],
[#text("C1-23 — Regenerate current optional defaults and schema formatting")],
[#text("Generated schema/default alignment. Acceptance is agreement with serialized public settings; no contraction benchmark applies.")],
[#text("C1-24 — Clarify the scope of older large-expression observations")],
[#text("Documentation accuracy. Earlier memory observations do not certify a current speedup; this audit replaces broad claims with workload-specific evidence.")],
)
]

= #text("Commit 2: motivation rows")

#block[
#set text(size: 8.6pt)
#table(columns: (0.35fr, 0.65fr),
table.header([#text("Change")], [#text("Evidence and conclusion")]),
[#text("C2-01 — Shared generalized CFF engine and selectable representation")],
[#text("Requested raised-power capability. Shared CFF contracts and complete signed/pinched values provide acceptance; an unsupported parent input is not an old wrong-result bug.")],
[#text("C2-02 — Retain independent repeated propagator occurrences")],
[#text("Independent occurrence ownership is covered by current complete-value probes. No isolated merged-versus-independent performance comparison or old wrong-value counterexample is established here.")],
[#text("C2-03 — Stream component products and shorten intermediate lifetimes")],
[#text("Source records a GL16 30GB resource failure. Current streaming/eager comparisons show no speed or RSS advantage on the reached workloads; the larger stress products have only one component.")],
[#text("C2-04 — Reuse exact bases and compact completed products")],
[#text("Exact base caching helps the reached high-rank cases. Compaction adds cost without reducing the larger double-triangle-derived output; necessity for other multi-component inputs remains open.")],
[#text("C2-05 — Correct shared terminal/embedded contour normalization")],
[#text("Source records an embedded-versus-standalone signed contour mismatch exposed by the unchanged depth-three banana. Correct common normalization is a mathematical requirement, independent of speed.")],
[#text("C2-06 — Correct diagnostic numerator arithmetic")],
[#text("Fresh exact old-method replay: ")#text(font: "DejaVu Sans Mono", size: 0.86em, "-2**2")#text(" gives 4 instead of -4, ")#text(font: "DejaVu Sans Mono", size: 0.86em, "2**3**2")#text(" gives 64 instead of 512, and an oversized exponent wraps. Current public contracts reject/correct these.")],
[#text("C2-07 — Validate bounds before fast paths and reject malformed graph signatures")],
[#text("Source records accepted unknown bound 999 and a malformed-signature panic. Fresh current public probes reject both; this audit does not rebuild every old validation branch.")],
[#text("C2-08 — Additional input arithmetic/shape validation")],
[#text("Fresh current arithmetic, shape, scale and graph-balance contracts pass. Some guards address source-level boundary counterexamples rather than an observed physical failure.")],
[#text("C2-09 — Rational coefficient storage across shared engine and consumers")],
[#text("Rational storage expresses the coefficient domain and removes unsupported symbolic states. Rational-valued Atoms were already exact. No rounding bug or isolated storage speedup is claimed.")],
[#text("C2-10 — Shared expression assembly, surface copying and builder consolidation")],
[#text("KISS consolidation: one assembler owns interning, copying and finalization. Complete-value coverage supports behavior; there is no standalone measured benefit or old bug for each moved line.")],
[#text("C2-11 — Preserve zero values for filtered and infinite denominator trees")],
[#text("Concrete numeric/symbolic zero-contract mismatch: absent/Infinite denominator trees could panic or evaluate near 1e80. Fresh public zero/fusion checks pass.")],
)
]

= #text("Commit 3: motivation rows")

#block[
#set text(size: 8.6pt)
#table(columns: (0.35fr, 0.65fr),
table.header([#text("Change")], [#text("Evidence and conclusion")]),
[#text("C3-01 — Integrate shared raised-energy CFF into GammaLoop and evaluator persistence")],
[#text("Necessary GammaLoop integration of the new shared engine, with generated evaluators and persistence. Physical captures and exact scalar/source certificates exercise the boundary.")],
[#text("C3-02 — Exact owner-preserving local-4D UV source reconstruction")],
[#text("Source contains complete GL04 numerator and denominator reconstruction certificates isolating owner/routing errors. Equal denominator squares alone cannot certify an odd numerator sign.")],
[#text("C3-03 — Typed direct/local/projected UV routes and explicit selector-free sums")],
[#text("Requested direct/projected routes plus removal of stale state. Explicit selector truth tables are the behavioral contract; API cleanup has no separate performance claim.")],
[#text("C3-04 — Owner-specific integrated UV mass scaling")],
[#text("Source records the GL04 finite-localizer coefficient discrepancy and enclosing/disjoint controls. Owner-specific vacuum-mass scaling is a correctness fix, not a global mass rescaling.")],
[#text("C3-05 — Keep natural Taylor denominator topologies and factorized numerator ownership")],
[#text("Source records expensive common-denominator reconstruction and the natural Taylor-topology requirement. This audit does not isolate old-versus-current topology collection performance.")],
[#text("C3-06 — Factorized degree analysis and independent outer rank envelopes")],
[#text("Source records slow factorized degree/rank-envelope analysis. No isolated timing of old degree-analysis methods is established here; combined generation timing cannot apportion its benefit.")],
[#text("C3-07 — Bound exact-CFF cache keys to each request capacity")],
[#text("Request-local capacity prevents union envelopes from enlarging unrelated requests. Core cache-removal measurements keep capacities fixed, so they do not benchmark the old union-envelope policy.")],
[#text("C3-08 — Choose up to three certified derivative-energy assignments by generated map size")],
[#text("Core experiment compares complete scoring against the first certified proposal, including downstream cost. It retains proposal planning and certificates; fewer generated rows alone are not a speedup certificate.")],
[#text("C3-09 — Route child-soft carriers into parent hard chart")],
[#text("Source identifies GL08 retaining an unintegrated child momentum. Compatible parent-chart routing is the recorded fix; no fresh reduced old/current GL08 replay is claimed here.")],
[#text("C3-10 — Resolve symbolic constants before SymJIT compilation")],
[#text("Pinned source records unresolved pi^-3 producing infinity in SymJIT while the mapped evaluator is finite. Resolve constants at the compiler boundary; no speed motivation is needed.")],
[#text("C3-11 — Validate source integer arithmetic and raised residue selections")],
[#text("Concrete invalid-input and overflow boundary cases motivate fallible selection/conversion. No typical physical graph reaching the integer limits is claimed.")],
[#text("C3-12 — Uniform nonzero M evaluator parameter and finalized expression ownership")],
[#text("Uniform nonzero M and finalized expression ownership simplify the evaluator contract. This is an explicit API choice, not a demonstrated rounding fix or isolated optimization.")],
)
]

= #text("Commit 4: motivation rows")

#block[
#set text(size: 8.6pt)
#table(columns: (0.35fr, 0.65fr),
table.header([#text("Change")], [#text("Evidence and conclusion")]),
[#text("C4-01 — Command templates and invocation scope")],
[#text("Source contains real CLI old-fail controls for nested inline/boot template scope. Lexical inheritance and local overrides are observable behavior; graph timing is irrelevant.")],
[#text("C4-02 — Explicit output, state detection/versioning and history persistence")],
[#text("Source reproduces unsaved scalar-box 3D output writing into the state directory. Explicit write/version boundaries protect persistence. A Serde pair-list adapter now preserves structured additional-weight keys in inspect JSON; round-trip and physical inspection tests cover the public output contract.")],
[#text("C4-03 — 3D representation build/validate/export CLI")],
[#text("Requested 3D inspection/export workflow. Complete exported-function evaluations are the acceptance criterion; output paths and saved/reloaded functions are public contracts.")],
[#text("C4-04 — Benchmark command, timing reports and temporary minimal evaluation settings")],
[#text("Requested benchmark observability and temporary-settings restoration. Existing public success/error restoration checks support it; no old settings-loss reproduction or large-batch memory bound is claimed.")],
[#text("C4-05 — UV profiling and stable finite-range numerical reports")],
[#text("Concrete extreme-value counterexamples motivate nonzero arbitrary-precision bounds and finite UV sample validation. Small physical integrands cannot exercise every reporting exponent boundary.")],
[#text("C4-06 — Numerical-stability histograms, dashboard and integration checkpoint outputs")],
[#text("New statistics/checkpoint workflow and required consumers. Public batch, persistence and restoration contracts motivate it; no separate physical bug or reporting-overhead improvement is established.")],
[#text("C4-07 — Resource-aware test execution, Nix inputs and tracing")],
[#text("Recorded resource/build issues and explicit tooling policy. Platform-specific memory accounting is not inferred from Linux performance measurements.")],
[#text("C4-08 — Behavioral numerical tests and archive unconsumed scalar snapshots")],
[#text("Test-detection defects are concrete: scale-dependent tolerances accept a factor-two error, nonfinite predicates can pass, and unused snapshots assert nothing. Independent value oracles are the remedy.")],
[#text("C4-09 — Remove obsolete wrapped model and migrate examples/metadata")],
[#text("Cleanup of an unreferenced wrapped model plus necessary examples/schemas. Reference and live-model audits justify removal; no independent speedup is claimed.")],
[#text("C4-10 — Preserve precise stability digits before f64 narrowing and report inclusive bounds")],
[#text("Concrete precise-output errors: narrowing 1e-400 to f64 zero loses information; an inclusive overflow bin must not be labelled with a strict bound. Current boundary tests cover the contracts.")],
[#text("C4-11 — Process import/reference contracts and named-integrand provenance")],
[#text("Public import/provenance contracts. The audit itself encountered inferred-all-cuts versus explicit-process inputs, confirming why the process specification matters; this is not a core CFF mismatch.")],
)
]

= #text("Commit 5: motivation rows")

#block[
#set text(size: 8.6pt)
#table(columns: (0.35fr, 0.65fr),
table.header([#text("Change")], [#text("Evidence and conclusion")]),
[#text("C5-01 — Physical amplitude/LU and threshold phase normalization")],
[#text("Independent signed scalar/left-right contour counterexamples establish the phase requirement. A common wrong phase can survive strategy-equality checks; no speedup is needed.")],
[#text("C5-02 — Inverse-process sewing, complex couplings and explicit CP opt-in")],
[#text("Complex charged-current and scalar spin-matrix oracles reject the proposed extra-adjoint route; this is not a fresh reproduction of an original-main defect. Real-coupling gluon graphs cannot certify this contract. CP optimization remains an explicit assumption.")],
[#text("C5-03 — Consistent W/Z virtual propagators, covariant cut multiplets and physical event labels")],
[#text("Full covariant versus physical polarization sums and Ward identities motivate gauge sewing/model consistency. QCD-only benchmarks do not establish the electroweak contract.")],
[#text("C5-04 — Reject unsupported complex masses/custom pole rules and use ghost anticommutation")],
[#text("Unsupported real-energy/complex-mass and custom-pole boundaries require validation. Ghost loop/exchange signs share the anticommutation predicate; this is a supported-domain contract.")],
[#text("C5-05 — Feyngen accepted-numerator lookup optimization")],
[#text("Physical repeated-comparator measurements show a 7.07× median sign-rejection benefit and no scalar-matching benefit on twelve rejecting pairs. All 48 end-to-end observations preserve complete exports, but no workflow gain is established; current four-worker runs are slower than all ablations in both measured rounds. With only two rounds and two-candidate buckets, neither a precise speed effect nor long-bucket scaling is established. Canonical-numerator matching remains unexercised. The global mutex is a serialization ablation, not the complete old map/clone architecture.")],
[#text("C5-06 — Retain fixed-momentum internal bridges during initial-state subtraction")],
[#text("Recorded reproducible GL3 bridge bug: an internal fixed-momentum edge touches neither initial-cut endpoint, yet the old helper removes a vertex and destroys cut-side loop rank.")],
)
]

= #text("Commit 6: motivation rows")

#block[
#set text(size: 8.6pt)
#table(columns: (0.35fr, 0.65fr),
table.header([#text("Change")], [#text("Evidence and conclusion")]),
[#text("C6-01 — Projected tensor kernels/scalar coefficient preservation and whole-numerator alternative")],
[#text("Both Vakint input modes and forest owners agree on complete nonzero quark-bubble Laurent output. Timings show no consistent mode advantage there; full double-triangle UV exceeds declared limits.")],
[#text("C6-02 — Complete d-dimensional Dirac/metric closure before Laurent truncation")],
[#text("Independent d-dimensional Dirac/metric and scalar-master oracles expose finite evanescent terms lost by premature 4D specialization. Both input modes need the same closure boundary.")],
[#text("C6-03 — Account for coefficient/normalization epsilon poles and reject insufficient depth")],
[#text("Laurent order budgeting must add numerator-coefficient and normalization pole orders. Exact single/double-pole counterexamples justify the change without a performance claim.")],
[#text("C6-04 — Validate Lorentz domains, opaque tensor slots and scoped gamma5 errors")],
[#text("Source records odd-open/renamed-dummy boundary failures; gamma5 rejection is a scoped unsupported-domain policy. Cancellation cannot hide invalid tensor slots before validation.")],
[#text("C6-05 — Remove exactly-zero factorized scalar coefficient subtrees")],
[#text("Source records GL29/35/38/40/46 generation failures and a concrete exactly-zero scalar coefficient. Removing that zero while retaining nonzero spectator factors is the value/factorization contract.")],
[#text("C6-06 — Canonical backend serialization and post-restoration notation")],
[#text("Source gives corrupt FORM input (1+2i)*(a+b)*tensor losing parentheses, canonical-name corruption and invalid decimal-complex tokens. Adapter correctness, not speed, motivates the fixes.")],
[#text("C6-07 — Preserve signed numerical momenta, complex parameters and combined estimator")],
[#text("Signed real momenta, invalid complex inputs and transactional settings have concrete adapter contracts. Manual PySecDec integration evidence is not fresh automatic coverage in this audit.")],
[#text("C6-08 — Normalize completed analytic return before reinsertion and migrate forest owners")],
[#text("Source gives a sunrise analytic-return normalization discrepancy with equal exact Laurent algebra. Both forest owners and modes require the same reinsertion normalization; not a new integral value.")],
)
]

= #text("Commit 7: motivation rows")

#block[
#set text(size: 8.6pt)
#table(columns: (0.35fr, 0.65fr),
table.header([#text("Change")], [#text("Evidence and conclusion")]),
[#text("C7-01 — Current reports, Typst/PDF and durable source/validation provenance")],
[#text("Requested documentation and provenance. The report records actual commands, complete-result coverage, measured costs, current dispositions and unsupported motivation claims. No runtime speedup applies.")],
)
]

= #text("Pinned benchmark validation evidence")

#text("Counts describe the accepted command families in this audit. Warm-ups are retained separately from measured observations; repeated calls within one process are not independent test executions. The report does not combine these heterogeneous checks into a headline test total.")

#block[
#set text(size: 8.6pt)
#table(columns: (0.42fr, 0.29fr, 0.29fr),
table.header([#text("Family")], [#text("Fresh result")], [#text("Complete-result evidence")]),
[#text("Foundation regressions")],
[#text("Current 22/22; scoped baseline 12 passes and 10 observed distinguishing failures")],
[#text("Seven defect/factorization cases and three supported-syntax probes; unchanged product tests")],
[#text("Foundation timing and storage")],
[#text("308 successful observations: 264 measured, 44 warm-ups; 32 storage smokes pass")],
[#text("Independent component/value oracles outside the timed operation")],
[#text("Physical contraction matrix")],
[#text("175 successful observations: 150 measured, 25 warm-ups; two guard terminations")],
[#text("All ten successful output hashes checked at three exact shared algebraic points")],
[#text("Deferred scalar-weight stress")],
[#text("28 successful observations: 24 measured, four warm-ups")],
[#text("Three distinct full output hashes checked with the independent Fraction oracle")],
[#text("Shared CFF")],
[#text("210/210 observations: 180 measured, 30 warm-ups; 28/28 public contracts")],
[#text("Full-expression trace linkage, 285 numerical observations and 15 independent pinch oracles")],
[#text("Core GL04 reconstruction")],
[#text("56/56 observations: 48 measured, eight warm-ups; four GL16 timeouts")],
[#text("Three full inspections and nine complete computed forest/node exports")],
[#text("Vakint modes and owners")],
[#text("56/56 observations: 48 measured, eight warm-ups; four full two-loop limits")],
[#text("Full nonzero Laurent oracle and direct linkage of every timed export")],
)
]

#text("Feyngen's body and end-to-end counts are reported with their timing tables and actual schedule. The reproduction bundle records formatting, locked checks and builds of the diagnostic executables, exact compiler/dependency hashes, commands and failures. These counts belong to the pinned benchmark sources and isolated interventions. Final reconstructed-stack commands and results are supplied in the ")#link("raised-energy-cff-stack-review.md")[#text("stack review")]#text(", without reusing benchmark or earlier workspace-suite totals as fresh validation.")

= #text("Evidence and reproduction")

#text("The ")#link("raised-energy-cff-motivation-inventory.json")[#text("complete inventory")]#text(" and ")#link("raised-energy-cff-motivation-evidence.tar.gz")[#text("evidence bundle")]#text(" accompany this report. The bundle contains source pins, exact intervention patches, public-API probes, physical graph/model captures, input hashes, full timing manifests, per-output comparison receipts, resource failures and independent methodology review. It excludes executable binaries and license material. The reproducibility guide distinguishes valid final observations from discarded setup measurements and states how to provide a licensed local environment.")

#text("The validation counts in this audit come from its own pinned commands and receipts. No manual PySecDec integration is run, no full gauge-invariant two-loop gg acceptance is claimed, and no result demonstrates every possible contraction tree or arbitrary sparse domain. Unsupported performance necessity is a finding, not a silently completed benchmark.")
