// Self-contained typesetting of raised-energy-cff-stack-review.md.
// Current implementation, validation results, and coverage limits.
#set document(date: datetime(year: 2026, month: 9, day: 11), title: "Raised-energy CFF stack review", author: "Codex",
  description: "Seven-commit reconstruction: correctness, Rust idioms, KISS and outward test contracts; current implementation and completed validation.")
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
#text(size: 9pt, fill: rgb("667085"))[11 September 2026 · Seven functional commits]
#v(4mm)

= #text("Assessment")

#text("The stack is reconstructed as seven functional commits, with the documentation in jj change ")#text(font: "DejaVu Sans Mono", size: 0.86em, "vmrowquy")#text(" and ")#text(font: "DejaVu Sans Mono", size: 0.86em, "codex/raised-energy-cff-reviewed")#text(" as the final local bookmark. Inspection JSON preserves structured additional-weight identifiers, and the three compact-parser regressions check public execution results. Powered-dot freshening remains required by an independent FORM-boundary counterexample. ")#strong("All seven reconstructed boundaries have passing certificates, and all ten final runtime selections pass at the amended source.")#text(" Complete path accounting and the execution evidence below refer to these reconstructed source revisions; final document and history review is recorded separately in the external closure receipt.")

#text("The implementation combines exact rational CFF coefficients, independent repeated occurrences, owner-specific integrated UV mass handling, shared contour normalization, and physical phase conventions. Complete signed contour oracles, separate numerator and denominator reconstruction certificates, and independent dimensional and physical references provide the strongest correctness evidence. Agreement between execution strategies is supporting evidence; a shared normalization can give both strategies the same phase error.")

#text("The ")#link("raised-energy-cff-motivation-audit.md")[#text("change-motivation audit")]#text(" (")#link("raised-energy-cff-motivation-audit.pdf")[#text("PDF")]#text(") covers 73 logical categories at its pinned benchmark source. It supplies controlled measurements on the physical gluon double triangle with a three-gluon vertex and top-quark loop, additional physical captures, and targeted inputs for otherwise inactive mechanisms. Its observations distinguish measured improvements, concrete defects, requested capabilities and unproved optimization benefits. Those measurements and source identities remain intact; they do not supply fresh validation or complete path coverage for the reconstructed stack.")

#text("Rust idiomaticity is good at the domain boundaries. Borrowed expressions, owned results, concrete rationals, typed IDs, explicit alternatives, and fallible conversions express useful invariants. The design keeps expression assembly, command preparation, source assignment, evaluator construction, and Laurent restoration in their existing owners. Source reconstruction and tensor algebra remain substantial because they represent distinct mathematical requirements.")

= #text("Functional commits")

#text("Each commit groups its signature changes, callers, settings, schemas and required tests under one functional owner, with no temporary compatibility layer. Formatting, locked workspace checks, all-target builds, Clippy and owning behavioral selections pass at every commit with its ancestors. C1 and C2 retain their exact unchanged-revision certificates; C3–C7 have fresh checks after the GL00 certificate amendment. All seven commits use ")#text(font: "DejaVu Sans Mono", size: 0.86em, "Lucien Huber <im@lcnbr.ch>")#text(" as author and committer.")

#block[
#set text(size: 8.6pt)
#table(columns: (0.08fr, 0.22fr, 0.70fr),
table.header([#text("#")], [#text("Commit")], [#text("Scope")]),
[#text("1")],
[#text(font: "DejaVu Sans Mono", size: 0.86em, "97eb06cb")],
[#text("Symbolic and tensor foundations")],
[#text("2")],
[#text(font: "DejaVu Sans Mono", size: 0.86em, "2ba79605")],
[#text("Exact generalized CFF generation")],
[#text("3")],
[#text(font: "DejaVu Sans Mono", size: 0.86em, "1f0d4fd6")],
[#text("Raised-energy CFF and local UV reconstruction")],
[#text("4")],
[#text(font: "DejaVu Sans Mono", size: 0.86em, "48104fa4")],
[#text("Command execution and evaluation workflows")],
[#text("5")],
[#text(font: "DejaVu Sans Mono", size: 0.86em, "32a57c5d")],
[#text("Physical phases and model sewing")],
[#text("6")],
[#text(font: "DejaVu Sans Mono", size: 0.86em, "ee798ce2")],
[#text("D-dimensional integrated UV algebra")],
[#text("7")],
[#text(font: "DejaVu Sans Mono", size: 0.86em, "vmrowquy")],
[#text("Review, architecture, and validation documents")],
)
]

#text("The ")#link("raised-energy-cff-stack-source-mapping.md")[#text("source mapping")]#text(" identifies the full functional commit hashes and source ownership. The first six commits pin the executable implementation. Documentation uses stable jj change ")#text(font: "DejaVu Sans Mono", size: 0.86em, "vmrowquy")#text("; the final local bookmark is ")#text(font: "DejaVu Sans Mono", size: 0.86em, "codex/raised-energy-cff-reviewed")#text(".")

= #text("1. Symbolic and tensor foundations")

#text("Symbolic numerators retain their factorization through dummy allocation, tensor materialization, and contraction. Parser clones share fresh-index allocation. Explicit-name reservations and collision skipping prevent compact dummy names from capturing existing indices; shared allocation already existed at the parent boundary. Canonicalization preserves external and dual slots. Compact-vector materialization accepts supported, unambiguous axes. Scalar aliases and deferred sums are completed at the existing result boundary in both execution strategies.")

#text(font: "DejaVu Sans Mono", size: 0.86em, "MinIntermediateCost")#text(" estimates symbolic copying and normalization with saturating arithmetic. It is a scheduling heuristic and does not change coefficients. Independent coordinate enumeration checks complete contractions across strategies, sparse and dense inputs, metric signs, and both operand orders. Dense–dense and mixed multi-index contractions already worked through the generic contraction API; the generalized grouped specialization changes their cost. The complete-contraction assertions permit different intermediate scores, storage variants and allocation indices. Three public parse/execute regressions compare a weighted Minkowski component sum for weights ")#text(font: "DejaVu Sans Mono", size: 0.86em, "-2")#text(", ")#text(font: "DejaVu Sans Mono", size: 0.86em, "a")#text(" and ")#text(font: "DejaVu Sans Mono", size: 0.86em, "a+b")#text(", and preserve complete opaque function values for two ambiguous tensor arguments. They do not constrain the private materializer or its allocated slots.")

#text("The symbolic dependency and Vakint prerequisites support exact imaginary constants and namespace-preserving FORM round trips. Powered scalar-dot copies require independent dummy indices at the whole-numerator FORM boundary: otherwise ")#text(font: "DejaVu Sans Mono", size: 0.86em, "(k·p+k·q)^2")#text(" loses ")#text(font: "DejaVu Sans Mono", size: 0.86em, "1/D")#text(" from both diagonal angular-average terms. Commit 1 includes a public tensor-reduction regression against ")#text(font: "DejaVu Sans Mono", size: 0.86em, "k^2*(p^2+2*p·q+q^2)/D")#text(", with ")#text(font: "DejaVu Sans Mono", size: 0.86em, "D=4-2*eps")#text(", using the whole-numerator API available there. Rust scalar round trips alone do not expose this error. The ")#link("raised-energy-cff-powered-dot-freshening.md")[#text("powered-dot follow-up")]#text(" retains the independent counterexample; the measured conversion overhead serves a correctness requirement. Both old sparse fibers and the generalized kernel can visit absent coordinates. The motivation audit measures all six exposed scheduling choices, storage kernels, aliases and deferred sums on physical captures and controlled stress inputs.")

= #text("2. Exact generalized CFF generation")

#text("The shared engine keeps repeated occurrences independent, including when physical energies coincide. Numerator capacity belongs to an occurrence; cache reuse cannot redistribute it among equal-energy copies. The four ")#text(font: "DejaVu Sans Mono", size: 0.86em, "LinearEnergyExpr")#text(" coefficient domains and ")#text(font: "DejaVu Sans Mono", size: 0.86em, "CFFVariant::prefactor")#text(" use ")#text(font: "DejaVu Sans Mono", size: 0.86em, "Rational")#text(". Conversion to symbolic expressions occurs through ")#text(font: "DejaVu Sans Mono", size: 0.86em, "Atom::num")#text(" at output boundaries.")

#text(font: "DejaVu Sans Mono", size: 0.86em, "ExpressionAssembler")#text(" owns surface interning, source copying, variant insertion, and label finalization. Cache keys include source topology and relevant contour and capacity state. Complete residue maps determine which contributions can be combined. The terminal builder enforces its rank requirement directly.")

#text("Public APIs validate compact edge IDs, endpoints, signature dimensions, and cut aliases before indexing. Balance accumulation uses i64 for sums of valid i32 routes. Numerical evaluation rejects incompatible shapes, zero uniform scale, and unsupported integer powers. Arithmetic parsing preserves precedence, right association, unary signs, and checked exponents. Empty trees and trees with infinite denominators obey their symbolic zero contract; fusion preserves zero contributions.")

#text("The tests cover signed simple and repeated poles, connected and disconnected components, both occurrence owners, reversed routing, lower sectors, scale invariance, and execution strategies. Tree transformations compare complete denominator algebra; selectors use an explicit truth table; cut checks compare complete weighted sums. Fresh default, all-feature and no-default-feature test selections pass; all-feature and no-default-feature all-target checks also pass from commit 2 onward.")

= #text("3. Raised-energy CFF and local UV reconstruction")

#text("GammaLoop owns physical source construction, denominator multiplicity, energy assignments, UV orchestration, and evaluator preparation. The shared engine owns exact CFF algebra. Physical edges, exact occurrences, production residue maps, and topological thresholds have distinct IDs. Physical source energies determine capacity; loop-momentum-basis coordinates route momenta without defining energy identity.")

#text("One immutable assignment plan controls generation capacity and numerator mapping. Original factors keep their source occurrences, and derived factors use only eligible copies. Exact reconstruction retains numerator signs under routing reversal. GL00 and GL04 child certificates independently check complete factorized common-chart numerators and denominator momenta, multiplicities, masses, physical owners and UV domains. The GL00 edge-5 temporal case uses an independent signed T0+T1 oracle; the GL04 case covers T0+T1+T2. Both retain the vacuum mass in full numerator expressions and the expansion mass in denominator metadata until the documented final mass identification. These checks start from the completed post-Taylor child source and do not certify every forest or the Taylor-producing stage. Separate GL04 source evidence also includes the cograph domain. Equality of denominators alone cannot establish the sign of an odd numerator.")

#text("Direct local 3D Taylor operations act on the complete generated CFF. Projected local 4D reconstruction uses completed raised-denominator terms and their complete source sums. The routes share evaluator assembly, production orientation catalogs, runtime validation, and exact source conversion. Grouped orientation sampling rejects incompatible map catalogs.")

#text("Integrated finite coefficients carry their physical owner through symbolic multiplication. A Taylor operation scales the vacuum mass of a contained owner, leaves a disjoint owner fixed, and rejects partial overlaps. The normalized localization kernel stays fixed. Output copies restore the physical mass symbol while stored sectors retain ownership for subsequent forest operations.")

#text("Behavioral tests check complete nested and disjoint UV functions, both local bubble contours plus one analytic finite addback, and exported forest values. DOT exports are re-imported and aggregated by physical forest and residue identity. Renormalization checks use signed Laurent expressions and a scale derivative. The nested vacuum fixture declares unit loop normalization explicitly, and the quark oracle aligns the external color basis using an independent color identity.")

#text("Saved-state manifest version 7 describes rational CFF storage. Incompatible manifests are rejected before decoding or overwriting a payload. Named generation records shared settings and requires consistent provenance for generated siblings. Persistence tests exercise rejected-operation preservation, public exports, and save/load behavior. Arbitrary-precision baseline expectations follow decimal promotion exactly.")

= #text("4. Command execution and evaluation workflows")

#text("Command preparation uses one placeholder parser and one environment-aware boundary. Nested blocks inherit lexical variables; local definitions override them. Template text is an enum payload, so reusable blocks can contain variables supplied at invocation. Preparation validates names and cycles and restores active-block tracking on returned errors.")

#text("Graph import options collect source, process, and naming policy at one boundary. Inline DOT is a literal argument. Explicit process specifications must agree with the graph kind and existing process definition. Rust schemas describe the serialized string-or-command contract, and generated schema references resolve.")

#text("Read-only and explicit state-write policies cover evaluation, inspection, benchmarks, profiling, and 3D outputs. Diagnostic 3D output defaults to a workspace in the current directory. CLI failures produce a failing process exit status. Resume distinguishes absent state from corrupt or unreadable state and checks workspace versions. Benchmark settings are restored after success and returned errors; panic recovery is outside the contract.")

#text("Retained additional weights serialize as ")#text(font: "DejaVu Sans Mono", size: 0.86em, "[key, value]")#text(" pairs, including ")#text(font: "DejaVu Sans Mono", size: 0.86em, "[]")#text(" for no weights. ")#text(font: "DejaVu Sans Mono", size: 0.86em, "GenericAdditionalWeightInfo.weights")#text(" reuses the existing ")#text(font: "DejaVu Sans Mono", size: 0.86em, "vectorize")#text(" Serde adapter; an explicit ")#text(font: "DejaVu Sans Mono", size: 0.86em, "T: Deserialize<'de>")#text(" bound retains the original deserialization requirement without strengthening the generic type. Structured threshold-counterterm identifiers remain values, so ")#text(font: "DejaVu Sans Mono", size: 0.86em, "inspect --json-output")#text(" can preserve them. JSON/bincode round-trip tests cover all four key variants and empty weights. A physical inspection regression compares numerical values, retained weights and grouping with retention disabled and enabled; it does not claim every event metadata field. The separate fresh GL04 replay checks complete retained event kinematics, cut metadata, weights and grouping under this toggle, exercises structured threshold-counterterm keys in JSON, and verifies the saved state stays unchanged. This replay supplies serialization and retention evidence; the independent signed contours and physical oracles described elsewhere supply physics evidence.")

#text("Stability output retains batch statistics until assembly. Median formatting preserves a decimal exponent without underflowing an arbitrary-precision value through f64. Timing rows distinguish evaluator time from inclusive integrand time. Long benchmark requests retain a full batch and can consume substantial memory.")

#text("Tests execute nested commands and inspect saved histories, public state, exported and reloaded functions, output paths, restoration, and errors. The 3D CLI regression evaluates its written artifact. UV fail-fast checks its reported failure against exhaustive results without choosing traversal order. Public completion counts remain valid contracts.")

= #text("5. Physical phases and model sewing")

#text("Physical normalization retains complete UFO propagator and vertex factors and applies common conversion at the cut-group boundary for ordinary and counterterm terms. Connected right-hand components include cut hairs, preserving contact interactions. Marked inverse-process vertices supply Hermitian partner couplings and spin matrices exactly once. Left/right threshold signs are checked against fixed pole-distribution derivatives and signed contour references.")

#text("Initial-state subtraction removes an endpoint only when it actually touches the initial-state cut. Fixed-momentum internal bridges retain their vertices and loop structure. The endpoint search returns ")#text(font: "DejaVu Sans Mono", size: 0.86em, "Option")#text(", expressing the valid absence of such an endpoint. The worker-count regression exercises this production behavior.")

#text("Gauge sewing uses model-declared covariant quartets with validated spin, statistics, charge, color, mass, propagator, and conjugation contracts. It closes complete final-state multisets before filtering graphs and reports physical event labels. Independent polarization sums in rest and boosted frames check the complete vector, Goldstone, and ghost result. Matching numerical masses alone does not identify gauge partners.")

#text("One anticommuting-particle predicate covers fermions and ghosts. Optional CP symmetrization is an explicit user assumption in generated data. Mass checks use current model values and the used propagator denominators. Complex-mass cutting rules and arbitrary non-Hermitian interactions are outside the supported contract.")

#text("Numerator grouping uses short global lookups and atomic search/insertion within topology buckets. Exact symbolic matching precedes sampled fallback. Worker-count tests compare complete values within each grouping strategy. Different strategies may select different momentum routing representatives, so pointwise equality between strategies is not a general contract. Sampled fallback remains probabilistic.")

#text("Signed amplitudes, phase space, flux, bubble discontinuity, virtual mirrors, spin matrices, Ward identities, and LO/NLO acceptances provide complementary oracles. Acceptances constrain real and imaginary parts. Graphwise scale checks require the analytic sign explicitly. Seven vertex-rule fixtures compare complete signed tensors after dummy-index alignment.")

= #text("6. D-dimensional integrated UV algebra")

#text("Projected-tensor and complete-numerator Vakint inputs share Lorentz validation, coefficient restoration, and Laurent truncation. Epsilon-dependent and complex coefficients retain poles, finite terms, and evanescent effects. Kernel depth accounts for coefficient and normalization poles before truncation; a factor such as ")#text(font: "DejaVu Sans Mono", size: 0.86em, "1/(4-2ε)")#text(" must participate while the required scalar-master orders are available.")

#text("Commit 6 extends the complete powered-scalar angular-average regression from commit 1 to both ")#text(font: "DejaVu Sans Mono", size: 0.86em, "project_onto_tensor_integrals")#text(" settings. The same independent ")#text(font: "DejaVu Sans Mono", size: 0.86em, "1/D")#text(" oracle checks both routes when the two-mode API becomes available, with no compatibility helper in the earlier commit.")

#text("The GammaLoop boundary restores projected slots, completes supported d-dimensional gamma algebra, substitutes scalar dimensions, and requests Laurent coefficients. Open representation slots retain symbolic dimensions. Unsupported tensor operators and traces return errors. Input mode and finite depth propagate consistently through cut forests, hedge-poset forests, and MUV renormalization.")

#text("Numerical settings preserve signed real momentum components and reject nonreal momentum inputs. General complex scalar coefficients have a distinct contract from pole masses. Fallible parameter collection completes before settings mutation. The numerical route combines sectors into one kernel so its uncertainty describes the complete integral.")

#text("Tests check full Laurent values, signed Minkowski components, alpha-renamed contractions, malformed-input errors, and mode parity. Reconstructed expressions and public Laurent accessors permit omitted or stored exact zeros while rejecting unexpected nonzero powers. Independent Dirac identities and scalar-master oracles cover nonzero evanescent finite terms. Public adapters exercise both production modes; injected analytic kernels provide additional bounded algebra checks.")

= #text("Behavioral test contracts")

#text("The reviewed tests assert complete values, errors, physical ownership, exported functions, and state restoration. For small symbolic expressions, equality uses the complete expanded difference and requires zero. Factorization-specific assertions use direct equality. Tensor comparisons align contracted dummy indices while preserving physical labels. Independent oracle values, signs, and tolerances remain explicit.")

#text("Private cache counts, allocation indices, planner ordering, storage layouts, and incidental expression formatting are outside these contracts. Deliberate public display output, physical coordinate identities, and required factorization remain observable behavior. Fixed fixture counts ensure zipped comparisons cover every expected item. Some unchanged repository tests still inspect internals; the review does not certify every test in the repository as implementation-independent.")

= #text("Rust patterns and KISS")

#block[
#set text(size: 8.6pt)
#table(columns: (0.32fr, 0.68fr),
table.header([#text("Pattern")], [#text("Effectiveness and tradeoff")]),
[#text("Borrowed ")#text(font: "DejaVu Sans Mono", size: 0.86em, "AtomView")#text(", owned ")#text(font: "DejaVu Sans Mono", size: 0.86em, "Atom")],
[#text("Inspection avoids copying expression trees; completed values and cache keys own their data. Clones are justified where source provenance must survive mutation.")],
[#text("Concrete ")#text(font: "DejaVu Sans Mono", size: 0.86em, "Rational")#text(" and checked conversion")],
[#text("The type enforces the coefficient domain. Explicit symbolic and precision boundaries, checked dimensions, exponents, and Laurent depth prevent silent narrowing.")],
[#text("Typed IDs and domain enums")],
[#text("Owners, occurrences, residue maps, thresholds, parser values, and UV routes stay distinct. These types do not validate arbitrary deserialized indices by themselves.")],
[#text("Paired optional data and enum payloads")],
[#text("An exact frame travels with its basis, and template text belongs to its template variant. Correlated state is represented together.")],
[#text("Existing builders and immutable plans")],
[#text("Expression assembly, source assignment, function construction, and evaluator preparation each have an owner responsible for their invariants.")],
[#text("Scoped shared state")],
[#text(font: "DejaVu Sans Mono", size: 0.86em, "Rc<Cell<_>>")#text(" and ")#text(font: "DejaVu Sans Mono", size: 0.86em, "Rc<RefCell<_>>")#text(" fit the recursive single-threaded parser; ")#text(font: "DejaVu Sans Mono", size: 0.86em, "Arc<Mutex<_>>")#text(" fits parallel topology grouping. Parser borrows are scoped, and parallel lookups release the global map lock before bucket comparisons.")],
[#text(font: "DejaVu Sans Mono", size: 0.86em, "Result")#text(", ")#text(font: "DejaVu Sans Mono", size: 0.86em, "?")#text(", ")#text(font: "DejaVu Sans Mono", size: 0.86em, "transpose")#text(", fallible collection")],
[#text("Absence remains distinct from failure. Exposed malformed input gets contextual errors, and collect-before-mutate protects settings. Internal invariant panics are not a universal recovery boundary.")],
[#text("Iterators and explicit loops")],
[#text("Lazy component products limit retained Cartesian state. Explicit loops remain clear for stateful contractions and polynomial accumulation. Local boxed dispatch handles changing component depth.")],
[#text("Semantic cache keys and deterministic maps")],
[#text("Reuse depends on graph, frame, and capacity state. Ordered maps provide reproducible inventories; behavioral tests allow storage strategies to change.")],
[#text("Const generics and saturating heuristics")],
[#text("Contraction policy reuses established rules, and overflow cannot make expensive candidates look cheap. The heuristic provides no global optimality or physical memory-bound guarantee.")],
[#text("Serde field adapter and explicit bound")],
[#text("Existing pair-list serialization preserves structured keys without a string grammar or custom helper. An explicit deserialization bound restores inference suppressed by the field adapter and keeps the public generic type unchanged.")],
)
]

#text("The KISS strengths are shared preparation and assembly paths, existing builders, direct invariant checks, and narrowly scoped ownership. No generic testing framework or compatibility layer is needed. Specialized source/UV algebra and independent certificates remain readability costs. Their owner, sign, denominator, and dimensional distinctions are mathematically necessary. The motivation audit does not demonstrate a benefit from compaction on its reached single-component product or from extra losing assignment proposals on its completed GL04 inputs. These remain bounded simplification candidates.")

= #text("Validation")

#strong("All seven reconstructed boundaries and all ten final runtime selections pass.")#text(" The executable implementation is pinned by commits 1–6, ending at ")#text(font: "DejaVu Sans Mono", size: 0.86em, "ee798ce279f280000537e4d8c05ecc83fc6ec2cf")#text(". Final runtime ran at ")#text(font: "DejaVu Sans Mono", size: 0.86em, "9e260684581ab2893eff4032c5ca8ea1fc2c18cb")#text(" (tree ")#text(font: "DejaVu Sans Mono", size: 0.86em, "9cebc793fe813699b9d689f1ab245715902ac04e")#text("). Documentation remains in stable jj change ")#text(font: "DejaVu Sans Mono", size: 0.86em, "vmrowquy")#text("; the external closure receipt records its final content/render review and exact source-to-bookmark checks.")

#text("The ")#link("raised-energy-cff-stack-review-validation.md")[#text("validation reference")]#text(" records exact fresh commands, selected identities, counts, configurations, diagnostics and warnings. The ")#link("raised-energy-cff-closeout-evidence.tar.gz")[#text("closeout evidence archive")]#text(" and ")#link("raised-energy-cff-closeout-evidence.json")[#text("archive receipt")]#text(" retain the underlying records. The retained initial-candidate runtime and motivation-audit measurements are not counted as amended-source runtime executions. C1/C2 reuse exactly matching revision certificates; C3–C7 are revalidated after the GL00 amendment.")

#block[
#set text(size: 8.6pt)
#table(columns: (0.56fr, 0.24fr, 0.20fr),
table.header([#text("Boundary checks")], [#text("Passed / recorded")], [#text("Guard time")]),
[#text("Formatting, locked workspace check, all-target build and Clippy")],
[#text("28 / 28")],
[#text("4727.474 s")],
[#text("Required shared-CFF and optional API feature checks")],
[#text("14 / 14")],
[#text("34.207 s")],
)
]

#text("Guard times include complete commands and build work. The following times are nextest test-execution times and are not end-to-end benchmark timings.")

#block[
#set text(size: 8.6pt)
#table(columns: (0.56fr, 0.24fr, 0.20fr),
table.header([#text("Owning boundary")], [#text("Passed / selected executions")], [#text("Test time")]),
[#text("Commit 1")],
[#text("693 / 693")],
[#text("9.532 s")],
[#text("Commit 2")],
[#text("182 / 182")],
[#text("31.653 s")],
[#text("Commit 3")],
[#text("401 / 401")],
[#text("173.215 s")],
[#text("Commit 4")],
[#text("805 / 805")],
[#text("473.811 s")],
[#text("Commit 5")],
[#text("399 / 399")],
[#text("183.766 s")],
[#text("Commit 6")],
[#text("424 / 424")],
[#text("105.744 s")],
)
]

#text("C1–C6 owning selections account for 2904 passing executions, zero failing executions and 6255 reported skipped entries. C1/C2 figures describe their retained unchanged-revision certificates, while C3–C6 figures describe fresh amended-boundary executions. Repeated boundaries and feature configurations remain separate executions. C7’s behavioral records are the same ten final runtime selections below and are counted once.")

#block[
#set text(size: 8.6pt)
#table(columns: (0.56fr, 0.24fr, 0.20fr),
table.header([#text("Final runtime selection")], [#text("Passed / selected")], [#text("Test time")]),
[#text("Curated GammaLoop")],
[#text("2065 / 2065")],
[#text("465.870 s")],
[#text("Phase regressions")],
[#text("212 / 212")],
[#text("30.817 s")],
[#text("Signed acceptances")],
[#text("4 / 4")],
[#text("120.076 s")],
[#text("Vakint / UV")],
[#text("241 / 241")],
[#text("73.977 s")],
[#text("Full scalar matrix, including slow cases")],
[#text("166 / 166")],
[#text("922.523 s")],
[#text("Shared CFF: default")],
[#text("37 / 37")],
[#text("5.191 s")],
[#text("Shared CFF: all features")],
[#text("108 / 108")],
[#text("21.324 s")],
[#text("Shared CFF: no default features")],
[#text("37 / 37")],
[#text("5.116 s")],
[#text("UFO model parity")],
[#text("1 / 1")],
[#text("0.620 s")],
[#text("Vertex rules")],
[#text("1 / 1")],
[#text("0.356 s")],
)
]

#text("Final runtime: ")#strong("2872 passing executions across 2207 distinct binary/test identities")#text(", zero failing executions, 3331 reported skipped entries, and 8 recorded Cargo/environment configurations. Distinct identities are binary/test pairs across all runs; configurations and repeated executions are reported separately. Skipped entries are outside the selected run and are not passing executions. Overlapping selections and feature configurations are counted separately as executions.")

#text("The full scalar matrix includes 26 standard and 140 slow cases (166 total), selected directly in release mode with ignored tests enabled. Licensed tests run with four workers, zero retries, disabled snapshot updates and a 30 GB process-tree RSS limit; Cargo uses four build jobs. The scalar command does not use the wrapper’s forced single-test setting. Clippy and the explicitly warning-tolerant curated/scalar commands pass with compiler/dependency warnings retained in the validation reference; passing exits do not mean warning-free compilation.")

#text("The fresh GL04 replay completes four inspection calls: plain and JSON output with additional weights disabled and enabled. Its frozen CLI SHA-256 ")#text(font: "DejaVu Sans Mono", size: 0.86em, "53f80aa6ca611a8566335c488ffb028d73d3d1ed3057e741a750edefa9aff0a4")#text(" exactly matches the binary copied from the successful C7 build. Disabled output matches the retained successful evaluation; the toggle preserves complete non-additional-weight content, including every event’s kinematics, cut metadata and complex weight, with group membership and multiplicity retained. Every emitted additional-weight entry is a well-formed key/value pair, and the output contains a structured ")#text(font: "DejaVu Sans Mono", size: 0.86em, "ThresholdCounterterm")#text("; the replay does not independently determine the expected physical key set. Only documented timing metadata and the empty-object-to-pair-list weight encoding are normalized. The saved state and binary hashes are unchanged. This is serialization/retention evidence, independent of the signed-contour and analytic physics oracles; it is not a performance measurement.")

#text("Coverage accounts for ")#strong("592 commit/path records and 518 distinct paths")#text(" at inventory revision ")#text(font: "DejaVu Sans Mono", size: 0.86em, "197b29f67ae065d2cf4d5d88b7658c29a0712af6")#text(". C1–C6 coverage consists of 558 records over 485 paths. Relative to the retained initial review, 554 exact before/after transitions reuse their existing review, 2 owning architecture patches replay exactly on the amended parent, and 2 C3 source/document records have fresh GL00 review. Earlier fresh-hunk and original-patch review methods remain linked in that retained lineage. C7 contains 34 documentation/evidence paths. Final C7 content and render review is recorded separately in the external closure receipt; inventory alone is not semantic review. These totals come from complete immutable-tree inventories, including the new evidence files. The motivation audit does not certify new paths, and coverage does not mean every unchanged repository line was reread.")

= #text("Scope and remaining limitations")

#text("Manual physical PySecDec integrations are outside automatic coverage. Other failing-class diagnostics and ")#text(font: "DejaVu Sans Mono", size: 0.86em, "aa_aa::important::aa_aa_local_inspect_backend_consistency")#text(" are excluded; the vertex-rule case is selected explicitly. The three documented triangle, double-triangle, and three-loop TBT decimal fixture points are outside the automatic suite, and their derivation is unverified.")

#text("The retained UFO/JSON model-asset audit covers 22 Lorentz structures, 108 couplings, 153 vertices, 72 parameter declarations, 43 particle contracts, three quartets and three vector propagators, with five deterministic substitutions for sampled expression checks. All 25 audited asset entries retain their modes and blobs. The fresh ")#text(font: "DejaVu Sans Mono", size: 0.86em, "fresh_sm_ufo_preserves_virtual_gauge_and_covariant_cut_states")#text(" regression checks covariant multiplets and twelve named W/Z, Goldstone and ghost propagator/state contracts. Local electroweak and Ward certificates do not establish a complete electroweak LU acceptance. Finite routing probes, sampled grouping, and Monte Carlo acceptances have mathematical or statistical limits.")

#text("The diagnostic f64 evaluator retains its documented singular-surface behavior. General serialized-graph validation, complex-mass cutting rules, and arbitrary non-Hermitian interactions are outside the supported guarantees. The motivation audit bounds its contraction results to the measured dimensions, densities and physical captures; extremely sparse large domains remain unmeasured. Large benchmark batches still require workload-specific resource bounds. Inspection JSON, the three compact-parser execution regressions and the powered-dot angular-average regression are included in their fresh owning selections. The GL04 replay certifies serialization and retention at its recorded point; it provides no performance measurement or independent physics oracle.")
