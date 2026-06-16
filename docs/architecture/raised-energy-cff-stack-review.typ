// Self-contained typesetting of raised-energy-cff-stack-review.md.
// Current validation status is recorded with explicit counts.
#set document(date: datetime(year: 2026, month: 9, day: 14), title: "Raised-energy CFF stack review", author: "Lucien Huber",
  description: "Eight-commit stack: current functionality, Rust idioms, KISS, outward test contracts and source-pinned validation.")
#set page(paper: "a4", margin: (x: 18mm, top: 19mm, bottom: 18mm),
  header: align(right)[#text(font: "DejaVu Sans", size: 8pt, fill: rgb("667085"))[GammaLoop / Reconstructed stack review]],
  footer: context [#line(length: 100%, stroke: 0.4pt + rgb("d0d5dd"))
    #v(2mm)
    #text(font: "DejaVu Sans", size: 8pt, fill: rgb("667085"))[Current implementation review #h(1fr) #counter(page).display("1 / 1", both: true)]])
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
#text(font: "DejaVu Sans", size: 26pt, weight: "bold", fill: rgb("163e59"))[Raised-energy CFF stack]
#v(2mm)
#text(font: "DejaVu Sans", size: 12pt)[Functionality · Test contracts · Rust idioms · KISS]
#v(3mm)
#text(size: 9pt, fill: rgb("667085"))[2026-09-14 · Seven functional commits and documentation]
#v(4mm)

// BEGIN GENERATED MARKDOWN BODY

= #text("Assessment")

#text("The stack contains seven functional commits and one documentation commit. It combines exact generalized CFF generation, local and integrated UV reconstruction, physical phase conventions and command workflows. The additional optimization commit owns canonical UV algebra and reuse of immutable graph preparation. All eight pinned commits use ")#text(font: "DejaVu Sans Mono", size: 0.86em, "Lucien Huber <im@lcnbr.ch>")#text(".")

#strong("Fresh validation: COMPLETE — 93 PASS of 93 planned stages.")#text(" This report distinguishes static review, accepted test executions and measurements on explicitly recorded sources. Passing strategy comparisons support correctness; independent signed contours and complete reconstruction certificates provide mathematical oracles. Detailed commands, attempted executions and exclusions are in the ")#link("raised-energy-cff-stack-review-validation.md")[#text("validation ledger")]#text(".")

= #text("Functional ownership")

#block[
#set text(size: 8.6pt)
#table(columns: (0.10fr, 0.23fr, 0.67fr),
table.header([#text("Commit")], [#text("Pinned identity")], [#text("Responsibility")]),
[#text("C1")],
[#text(font: "DejaVu Sans Mono", size: 0.86em, "665658b16898")],
[#text("Symbolic/tensor foundations and graph bookkeeping")],
[#text("C2")],
[#text(font: "DejaVu Sans Mono", size: 0.86em, "a9397ee0f86b")],
[#text("Exact shared generalized CFF generation")],
[#text("C3")],
[#text(font: "DejaVu Sans Mono", size: 0.86em, "3ea313789a1a")],
[#text("GammaLoop adapters, local UV reconstruction and evaluator preparation")],
[#text("C4")],
[#text(font: "DejaVu Sans Mono", size: 0.86em, "04637884f24f")],
[#text("Command, persistence, evaluation and benchmark workflows")],
[#text("C5")],
[#text(font: "DejaVu Sans Mono", size: 0.86em, "1f2cf6d8236d")],
[#text("Physical phases and model sewing")],
[#text("C6")],
[#text(font: "DejaVu Sans Mono", size: 0.86em, "8e4e4664f3f9")],
[#text("D-dimensional integrated UV algebra and analytic scalar-product protection")],
[#text("C7")],
[#text(font: "DejaVu Sans Mono", size: 0.86em, "6d7a9220129c")],
[#text("Canonical UV algebra, bounded occurrence allocation and graph-owned reuse")],
[#text("C8")],
[#text(font: "DejaVu Sans Mono", size: 0.86em, "oxsmszsroxrt")],
[#text("Current reports, source attribution and validation evidence")],
)
]

#text("Commit 7, ")#text(font: "DejaVu Sans Mono", size: 0.86em, "6d7a9220129c0fac37bdc0e04ebefccb1688c7be")#text(", pins the executable implementation. Commit 8 is identified by its stable jj change, avoiding a reference to its own final hash. The ")#link("raised-energy-cff-stack-source-mapping.md")[#text("source mapping")]#text(" records full revisions and trees. The final bookmark and documentation delta are verified separately during the final report update and closure.")

= #text("1. Symbolic and tensor foundations")

#text("Tensor parsing preserves scalar factors, open slots and independent dummy indices. Supported compact vectors accept scalar weights; ambiguous vector products inside scalar functions remain complete opaque arguments. Scalar aliases resolve through the existing result boundary. Odd powers repeatedly contract one fixed base square with the remaining tensor, preserving the requested exponent across all five leaf forms.")

#text("Graph rewiring retains slot order, edge flow and complete node/forest membership. Completed execution waves release unused tensor payloads while retaining live aliases and handles. Tests compare full contractions with independent component sums across dense, sparse and mixed storage, metric signatures and execution strategies. Factorization tests preserve scalar spectators. No universal contraction ordering or isolated speedup follows from these contracts.")

= #text("2. Exact generalized CFF generation")

#text("The shared engine stores coefficients as ")#text(font: "DejaVu Sans Mono", size: 0.86em, "Rational")#text(" and treats repeated propagator occurrences independently, including their numerator capacities. Parsing handles arithmetic precedence, signs and checked exponents; source validation checks IDs, endpoints, dimensions and cut aliases before indexing. Zero and infinite-denominator inputs have explicit value contracts.")

#text("Complete signed residues cover simple/repeated poles, disconnected components, reversed routing and execution strategies. Reuse keys retain the relevant topology, capacities and contour state. Rational storage expresses the coefficient domain; it is not itself a measured performance claim.")

= #text("3. Source adapters and evaluator preparation")

#text("GammaLoop owns physical incidence, denominator occurrences, energy assignments and residue-map construction. Immutable assignments bind capacity generation to the actual numerator mapping. Denominator momentum, mass, multiplicity and domain are certified separately from the complete signed numerator. A loop-momentum basis supplies coordinates, not ownership.")

#text("GL00/GL04 child certificates compare complete factorized numerators with independent Taylor oracles in a neutral chart; nested/disjoint tests cover surrounding forest composition. Finite tensor-sum boundaries may close compatible leaves before the normal contraction schedule, preserving separate sums and scalar spectators. Evaluator preparation stays at its existing numerical execution boundary and does not distribute graph numerators.")

= #text("4. Commands and evaluation workflows")

#text("One placeholder parser defines command templates and lexical scope. Validation catches references and cycles; temporary settings and active blocks are restored after returned errors. Manifest-based saved-state handling keeps persistence behind explicit save operations. Unsupported versions fail before payload decoding or replacement, and CLI errors produce failing exits.")

#text("A Serde pair-list adapter preserves structured additional-weight keys through inspection output. Tests execute commands, reload exported functions and saved state, check complete weights and history, and verify restoration after success/failure. Timing fields distinguish expression construction, tensor preprocessing, orchestration and literal Symbolica evaluator builds; inclusive intervals are not interchangeable.")

= #text("5. Physical phases and model sewing")

#text("Normalization retains complete propagator/vertex factors and applies the shared conversion at the cut-group boundary. Connected right-hand components include cut hairs; conjugate couplings and spin matrices enter once. Fixed-momentum bridge vertices remain unless their endpoints touch the initial cut.")

#text("Model multiplets validate spin, statistics, charge, color, mass, propagators and conjugation. Independent polarization sums, Ward identities, left/right mirrors and signed LO/NLO acceptances constrain complementary parts of the convention. Exact numerator matching precedes sampled fallback. Sampling remains probabilistic; arbitrary non-Hermitian interactions and complex-mass cutting rules are outside the supported contract.")

= #text("6. D-dimensional integrated UV algebra")

#text("Projected-tensor and complete-numerator Vakint modes share Lorentz validation, coefficient restoration and Laurent-depth accounting. Coefficient/normalization poles determine the necessary master-integral depth. D-dimensional algebra closes before scalar dimension substitution, preserving evanescent finite contributions.")

#text("Collision-safe scalar aliases protect completed dot products during analytic spin expansion while exposing open tensor slots. Powered-dot freshening also remains necessary at the whole-numerator FORM boundary: the independent average of ")#text(font: "DejaVu Sans Mono", size: 0.86em, "(k·p+k·q)^2")#text(" is ")#text(font: "DejaVu Sans Mono", size: 0.86em, "k^2*(p^2+2*p·q+q^2)/D")#text(", with ")#text(font: "DejaVu Sans Mono", size: 0.86em, "D=4-2*eps")#text(". Forest composition retains physical mass owners: contained owners participate in Taylor expansion, disjoint owners remain fixed and partial overlaps are rejected.")

= #text("7. Canonical local UV projection and reuse")

#text("Raw sectors retain physical owners and provenance for subsequent Taylor operations. A separate canonical projection view groups exact signed routing, mass, denominator polynomial/prescription and component domain. Physical incidence remains an independent witness. Only certified serial paths or pure cycles permit equivalent channels to merge into powered quotients; nonadjacent channels retain their topology.")

#text("Completed hard factors may use an unambiguous certified class. Physical-source-fixed/soft factors keep their restrictions; hard factors without a surviving class remain affine carriers. Odd numerator signs are independent of denominator evenness. Complete numerator and denominator reconstruction certificates precede native CFF generation.")

#text("The compressed allocator admits a baseline, an eligible packed alternative and one placement challenger. Native source-map rows score the bounded proposals; rank and stable order break ties. This is not a global optimum. One nonserialized graph-owned context reuses immutable preparations, assignment templates and mapped subtrees under complete keys and bounded retention. Component waves and bulk cut-key addition combine complete coefficient values while preserving factorization.")

#text("Same-frame ")#text(font: "DejaVu Sans Mono", size: 0.86em, "+N/-N")#text(" sectors may cancel during grouping. The grouped result returns the existing typed zero before requiring production maps, preserving every supported final cut order. The full projection/build regression tests that value rather than sector layout.")

= #text("Behavioral contracts and findings")

#text("Small scalar identities use exact subtraction and zero checks. Factorization-specific tests compare the complete retained factorized expression. Tensor oracles permit dummy-index renaming while preserving physical labels. Cache/build counts, storage layout, traversal/proposal order and private row-count goldens do not define numerical correctness; public retention bounds and supported value/error behavior remain valid contracts.")

#block[
#set text(size: 8.6pt)
#table(columns: (0.30fr, 0.25fr, 0.45fr),
table.header([#text("Finding")], [#text("Current contract")], [#text("Fresh execution evidence")]),
[#text("Projection cancellation")],
[#text("Complete typed zero at every supported final cut order")],
[#text("PASS: C7/uv-performance-outward, C8/curated")],
[#text("Private layout/cache assertions")],
[#text("Full values, supported errors, factorization and public bounds")],
[#text("Exact reviewed correction fingerprints in the source mapping; owning selections below")],
[#text("Physical route tolerance")],
[#text("Fixed 1e-9 complex-relative acceptance for GL00/GL01")],
[#text("PASS: C7/physical-local-uv-routes, C8/physical-local-uv-routes")],
[#text("Sunset allocation")],
[#text("Independent complete degree-five/seven signed contour values")],
[#text("PASS: C7/uv-performance-outward, C8/curated")],
)
]

#text("The signed sunset oracle covers all 126 distributions of degree five over five repeated occurrences and admitted degree-five/seven plans. Independent successive clockwise contours at energies ")#text(font: "DejaVu Sans Mono", size: 0.86em, "(2,3,7)")#text(" give ")#text(font: "DejaVu Sans Mono", size: 0.86em, "1213/16387080192")#text(" and ")#text(font: "DejaVu Sans Mono", size: 0.86em, "-365/47775744")#text(". The retained mathematical audit is separate from fresh application-test execution. GL00/GL01 route comparisons require finite, nonzero complex results at three fixed points within ")#text(font: "DejaVu Sans Mono", size: 0.86em, "1e-9")#text(" relative complex-norm tolerance. Reported stability estimates cannot relax that threshold; integrated and threshold counterterms are disabled to isolate this comparison.")

= #text("Rust patterns and KISS")

#block[
#set text(size: 8.6pt)
#table(columns: (0.32fr, 0.68fr),
table.header([#text("Pattern")], [#text("Effectiveness and limit")]),
[#text("Borrowed ")#text(font: "DejaVu Sans Mono", size: 0.86em, "AtomView")#text(", owned ")#text(font: "DejaVu Sans Mono", size: 0.86em, "Atom")#text(", concrete ")#text(font: "DejaVu Sans Mono", size: 0.86em, "Rational")],
[#text("Avoid copying during inspection, retain immutable witnesses and enforce the coefficient domain. Domain clarity does not prove speed.")],
[#text("Typed references and domain enums")],
[#text("Separate owners, canonical classes, occurrences, residue maps and valid input modes; raw/canonical UV types protect the recursion/projection boundary.")],
[#text("Immutable ")#text(font: "DejaVu Sans Mono", size: 0.86em, "Arc")#text(" payloads in a graph-owned context")],
[#text("Keep assignments and generated payloads together; scoped lifetime permits reuse without serializing caches into state.")],
[#text("Existing builders and fallible conversions")],
[#text("Keep validation and numerical preparation at established owners. Validate complete inputs before mutation or dispatch.")],
[#text(font: "DejaVu Sans Mono", size: 0.86em, "BTreeMap::entry")#text(" and consuming iterators")],
[#text("Accumulate exact values and release consumed state. Deterministic order aids diagnostics without becoming an oracle.")],
[#text("Fixed points and bounded proposals")],
[#text("Structural convergence or decreasing state counts support termination; proposal limits do not bound total runtime.")],
[#text("Result accessors and ")#text(font: "DejaVu Sans Mono", size: 0.86em, "Cow")],
[#text("Unify literal/lazy scalar interpretation and borrow when possible; unsupported interpretation returns an explicit error.")],
[#text("Temporary scalar aliases")],
[#text("Preserve completed scalar contractions during tensor expansion; collision checks prevent capture.")],
)
]

#text("Shared preparation, immutable bindings, one scalar-result boundary and one bulk accumulation path remove duplication. Canonical source reconstruction and tensor-boundary rewiring remain substantial algorithms with explicit invariants. Further cache layers and diagnostic state need representative evidence to justify their complexity. Strategy comparisons do not establish universal optimality.")

= #text("Fresh validation")

#block[
#set text(size: 8.6pt)
#table(columns: (0.12fr, 0.46fr, 0.42fr),
table.header([#text("Boundary")], [#text("Format / check / build / Clippy")], [#text("Owning behavioral selections")]),
[#text("C1")],
[#text("Four gates PASS")],
[#text("3/3 selections PASS; 811 accepted executions")],
[#text("C2")],
[#text("Four gates PASS")],
[#text("3/3 selections PASS; 182 accepted executions")],
[#text("C3")],
[#text("Four gates PASS")],
[#text("1/1 selections PASS; 275 accepted executions")],
[#text("C4")],
[#text("Four gates PASS")],
[#text("2/2 selections PASS; 623 accepted executions")],
[#text("C5")],
[#text("Four gates PASS")],
[#text("4/4 selections PASS; 217 accepted executions")],
[#text("C6")],
[#text("Four gates PASS")],
[#text("1/1 selections PASS; 242 accepted executions")],
[#text("C7")],
[#text("Four gates PASS")],
[#text("2/2 selections PASS; 291 accepted executions")],
[#text("C8")],
[#text("Four gates PASS")],
[#text("10/10 selections PASS; 2706 accepted executions")],
)
]

#block[
#set text(size: 8.6pt)
#table(columns: (0.30fr, 0.25fr, 0.45fr),
table.header([#text("Final selection")], [#text("Accepted passed / executed")], [#text("Status")]),
[#text("Curated GammaLoop (includes phase/Vakint subsets)")],
[#text("2135 / 2135")],
[#text("PASS")],
[#text("Full scalar matrix, including slow cases")],
[#text("167 / 167")],
[#text("PASS")],
[#text("Vertex rules")],
[#text("1 / 1")],
[#text("PASS")],
[#text("UFO model parity (API library)")],
[#text("1 / 1")],
[#text("PASS")],
[#text("Physical local UV routes: GL00/GL01")],
[#text("2 / 2")],
[#text("PASS")],
[#text("Spenso: shadowing")],
[#text("168 / 168")],
[#text("PASS")],
[#text("Spenso: no default features")],
[#text("50 / 50")],
[#text("PASS")],
[#text("Shared CFF: default features")],
[#text("37 / 37")],
[#text("PASS")],
[#text("Shared CFF: all features")],
[#text("108 / 108")],
[#text("PASS")],
[#text("Shared CFF: no default features")],
[#text("37 / 37")],
[#text("PASS")],
)
]

#text("Accepted final runtime: ")#strong("2706 executions over 2279 distinct binary/test identities")#text(", from 10/10 passing selections. Accepted configurations report 1572 excluded entries, which may overlap. Baseline and owning-prefix executions are excluded from these final totals. Phase/Vakint cases inside the curated selection are not added again as independent suites. Configuration overlap and owning-prefix executions remain distinct in the ledger. Pending, running, failed and unaudited stages prevent a complete result; zero accepted failures does not establish overall success.")

#text("The recorded configuration uses ")#text(font: "DejaVu Sans Mono", size: 0.86em, "--locked")#text(" and the ")#text(font: "DejaVu Sans Mono", size: 0.86em, "test_gammaloop")#text(" nextest profile. Workspace check/build/Clippy gates and most suites use the ")#text(font: "DejaVu Sans Mono", size: 0.86em, "dev-optim")#text(" Cargo profile; the full scalar matrix uses ")#text(font: "DejaVu Sans Mono", size: 0.86em, "--release")#text(". Gates include formatting, all-target workspace check/build and Clippy; separate shared-CFF/Spenso feature configurations have their own receipts. UFO parity selects the API library with ")#text(font: "DejaVu Sans Mono", size: 0.86em, "ufo_support")#text("; this is package-scoped optional-feature coverage. Cargo, nextest and Rayon use four workers under the recorded 30 GB process-tree guard. Final runtime selections allow up to ten seconds per process query; the guard still stops the pipeline on monitoring failure. Snapshot updates and automatic test retries are disabled. The full scalar matrix selects its underlying namespaces directly, including slow cases. Exact arguments, warning policy, overrides and integrity guards remain authoritative in the ledger.")

#text("Dependency: Symbolica ")#text(font: "DejaVu Sans Mono", size: 0.86em, "2.2.0")#text(", with ")#text(font: "DejaVu Sans Mono", size: 0.86em, "graphica")#text(" and ")#text(font: "DejaVu Sans Mono", size: 0.86em, "numerica")#text(", from ")#text(font: "DejaVu Sans Mono", size: 0.86em, "alphal00p/symbolica")#text(" at ")#text(font: "DejaVu Sans Mono", size: 0.86em, "4d0a833eb8e059d1f95bdae5abed2559830b235f")#text(". This is a specific locked snapshot. Tool versions: rustc 1.97.0 (2d8144b78 2026-07-07); cargo 1.97.0 (c980f4866 2026-06-30); cargo-nextest 0.9.140. License material is excluded from reports and archives.")

= #text("Source preservation and review scope")

#text("The mapping retains 37 remote and five review source attributions. Incoming ")#text(font: "DejaVu Sans Mono", size: 0.86em, "a1140c90c334ff58a3b040ff3eae06d3645bf0bb")#text(" is reviewed against ")#text(font: "DejaVu Sans Mono", size: 0.86em, "78395e3ab3ddd8d8f62b2f674d7488484eace197")#text(": 68 paths (65 text, 3 binary), 792 text hunks, with explicit reviewer attribution. The three serialized Symbolica blobs received byte/hash and source-provenance review; their contents were not semantically decoded. At the tested candidate, exact correction fingerprints cover 21 changed files and the whole-stack inventory contains 647 commit/path records over 549 paths; inventory does not establish an unperformed fresh rereview of unchanged hunks.")

#text("Immutable corrected source ")#text(font: "DejaVu Sans Mono", size: 0.86em, "392c98c3d2f70896ef94267addfcf610e2df2c88")#text(" and reconciled C8 ")#text(font: "DejaVu Sans Mono", size: 0.86em, "449a59f6f173d100c2da858c541faa1e0dccca6e")#text(" share tree ")#text(font: "DejaVu Sans Mono", size: 0.86em, "c87cc95044663a68413010aceb81f9475be342c1")#text(" exactly. Later final-report changes require separate exact, reviewed documentation deltas from both the reconciled C8 and the executed candidate. Production/build/test/config content must remain identical. Candidate receipt hashes and execution counts are preserved; documentation closure does not count as another test run.")

= #text("Limits and evidence")

#text("GL00/GL01 automatic route selections do not certify GL262. The source-scoped GL262 diagnostic reaches a Symbolica evaluator-construction panic with compilation disabled and supplies no complete matched runtime/performance result. No confirmed standalone reproducer of that large panic is established. The nonsymmetric ")#text(font: "DejaVu Sans Mono", size: 0.86em, "f(x,y)")#text(" import MRE demonstrates a separate argument-order defect and does not identify the panic's cause.")

#text("Fresh standalone evidence at executed revision ")#text(font: "DejaVu Sans Mono", size: 0.86em, "a3af927afc633c6e6c2474d3f03fa02e0242eaa4")#text(" has locked all-target check/build/Clippy passes, ")#strong("two positive controls (39 and 4)")#text(", and ")#strong("one reproduced import defect")#text(": the writer exits 0 and the reader exits the expected 101 for ")#text(font: "DejaVu Sans Mono", size: 0.86em, "f(y,x)")#text(" versus ")#text(font: "DejaVu Sans Mono", size: 0.86em, "f(x,y)")#text(". These are not workspace test counts or a GL262 replay. Source-equivalent standalone package; no additional execution. The ")#link("raised-energy-cff-mre-validation.md")[#text("standalone receipt")]#text(" records exact package/manifest/lock/config identity and original binary/result hashes; carrying that evidence adds no execution. Separately recorded formatting passes at ")#text(font: "DejaVu Sans Mono", size: 0.86em, "bbb8ff1db499a78b3b25c41f774c9a034acb5b4d")#text("; its tracked package/build inputs match observed C7.")

#text("Manual physical PySecDec integrations remain outside automatic coverage. Installed-Python subprocess tests and other platforms have no separate acceptance unless an explicit receipt records them. Finite coordinate probes and Monte Carlo acceptances have numerical/statistical limits; local Ward checks do not establish complete electroweak LU acceptance.")

#text("The ")#link("raised-energy-cff-motivation-audit.md")[#text("motivation audit")]#text(" and ")#link("local-4d-uv-performance.md")[#text("local UV measurements")]#text(" retain workload/source/binary identities. On source ")#text(font: "DejaVu Sans Mono", size: 0.86em, "911f1824f994611d723dff0f50416cb301718617")#text(", the local-only GL00/GL01 generation premiums are 4.45%/0.82%. A separate single-run scalar cohort includes integrated and threshold subtraction: GL21/base setup/generation takes 1.431860 times erased 3D. Other measured comparisons also show regressions; no general speedup is established. Those source-scoped observations do not certify this reconstructed stack. Stage timings, complete generation and saved-state sampling remain distinct.")

#text("Markdown and Typst are generated from the same report body, with non-whitespace visible-text parity checked. PDF compilation, every-page visual inspection, final bookmark and final document hashes belong to the external closure evidence produced after the final report update. This report embeds no hash of its own final content.")
