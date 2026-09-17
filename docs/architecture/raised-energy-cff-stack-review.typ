= Raised-energy CFF stack review
<raised-energy-cff-stack-review>
Current functionality, Rust patterns, KISS and behavioral contracts · 2026-09-14 · Lucien Huber

=== Assessment
<assessment>
The stack contains seven functional commits and one documentation commit. It combines exact generalized CFF generation, local and integrated UV reconstruction, physical phase conventions and command workflows. The additional optimization commit owns canonical UV algebra and reuse of immutable graph preparation. All eight pinned commits use `Lucien Huber <im@lcnbr.ch>`.

#strong[Fresh validation: COMPLETE --- 93 PASS of 93 planned stages.] This report distinguishes static review, accepted test executions and measurements on explicitly recorded sources. Passing strategy comparisons support correctness; independent signed contours and complete reconstruction certificates provide mathematical oracles. Detailed commands, attempted executions and exclusions are in the #link("raised-energy-cff-stack-review-validation.typ")[validation ledger];.

=== Functional ownership
<functional-ownership>
#figure(
  align(center)[#table(
    columns: (33.33%, 33.33%, 33.33%),
    align: (auto,auto,auto,),
    table.header([Commit], [Pinned identity], [Responsibility],),
    table.hline(),
    [C1], [`665658b16898`], [Symbolic/tensor foundations and graph bookkeeping],
    [C2], [`a9397ee0f86b`], [Exact shared generalized CFF generation],
    [C3], [`3ea313789a1a`], [GammaLoop adapters, local UV reconstruction and evaluator preparation],
    [C4], [`04637884f24f`], [Command, persistence, evaluation and benchmark workflows],
    [C5], [`1f2cf6d8236d`], [Physical phases and model sewing],
    [C6], [`8e4e4664f3f9`], [D-dimensional integrated UV algebra and analytic scalar-product protection],
    [C7], [`6d7a9220129c`], [Canonical UV algebra, bounded occurrence allocation and graph-owned reuse],
    [C8], [`oxsmszsroxrt`], [Current reports, source attribution and validation evidence],
  )]
  , kind: table
  )

Commit 7, `6d7a9220129c0fac37bdc0e04ebefccb1688c7be`, pins the executable implementation. Commit 8 is identified by its stable jj change, avoiding a reference to its own final hash. The #link("raised-energy-cff-stack-source-mapping.typ")[source mapping] records full revisions and trees. The final bookmark and documentation delta are verified separately during the final report update and closure.

=== 1. Symbolic and tensor foundations
<symbolic-and-tensor-foundations>
Tensor parsing preserves scalar factors, open slots and independent dummy indices. Supported compact vectors accept scalar weights; ambiguous vector products inside scalar functions remain complete opaque arguments. Scalar aliases resolve through the existing result boundary. Odd powers repeatedly contract one fixed base square with the remaining tensor, preserving the requested exponent across all five leaf forms.

Graph rewiring retains slot order, edge flow and complete node/forest membership. Completed execution waves release unused tensor payloads while retaining live aliases and handles. Tests compare full contractions with independent component sums across dense, sparse and mixed storage, metric signatures and execution strategies. Factorization tests preserve scalar spectators. No universal contraction ordering or isolated speedup follows from these contracts.

=== 2. Exact generalized CFF generation
<exact-generalized-cff-generation>
The shared engine stores coefficients as `Rational` and treats repeated propagator occurrences independently, including their numerator capacities. Parsing handles arithmetic precedence, signs and checked exponents; source validation checks IDs, endpoints, dimensions and cut aliases before indexing. Zero and infinite-denominator inputs have explicit value contracts.

Complete signed residues cover simple/repeated poles, disconnected components, reversed routing and execution strategies. Reuse keys retain the relevant topology, capacities and contour state. Rational storage expresses the coefficient domain; it is not itself a measured performance claim.

=== 3. Source adapters and evaluator preparation
<source-adapters-and-evaluator-preparation>
GammaLoop owns physical incidence, denominator occurrences, energy assignments and residue-map construction. Immutable assignments bind capacity generation to the actual numerator mapping. Denominator momentum, mass, multiplicity and domain are certified separately from the complete signed numerator. A loop-momentum basis supplies coordinates, not ownership.

GL00/GL04 child certificates compare complete factorized numerators with independent Taylor oracles in a neutral chart; nested/disjoint tests cover surrounding forest composition. Finite tensor-sum boundaries may close compatible leaves before the normal contraction schedule, preserving separate sums and scalar spectators. Evaluator preparation stays at its existing numerical execution boundary and does not distribute graph numerators.

=== 4. Commands and evaluation workflows
<commands-and-evaluation-workflows>
One placeholder parser defines command templates and lexical scope. Validation catches references and cycles; temporary settings and active blocks are restored after returned errors. Manifest-based saved-state handling keeps persistence behind explicit save operations. Unsupported versions fail before payload decoding or replacement, and CLI errors produce failing exits.

A Serde pair-list adapter preserves structured additional-weight keys through inspection output. Tests execute commands, reload exported functions and saved state, check complete weights and history, and verify restoration after success/failure. Timing fields distinguish expression construction, tensor preprocessing, orchestration and literal Symbolica evaluator builds; inclusive intervals are not interchangeable.

=== 5. Physical phases and model sewing
<physical-phases-and-model-sewing>
Normalization retains complete propagator/vertex factors and applies the shared conversion at the cut-group boundary. Connected right-hand components include cut hairs; conjugate couplings and spin matrices enter once. Fixed-momentum bridge vertices remain unless their endpoints touch the initial cut.

Model multiplets validate spin, statistics, charge, color, mass, propagators and conjugation. Independent polarization sums, Ward identities, left/right mirrors and signed LO/NLO acceptances constrain complementary parts of the convention. Exact numerator matching precedes sampled fallback. Sampling remains probabilistic; arbitrary non-Hermitian interactions and complex-mass cutting rules are outside the supported contract.

=== 6. D-dimensional integrated UV algebra
<d-dimensional-integrated-uv-algebra>
Projected-tensor and complete-numerator Vakint modes share Lorentz validation, coefficient restoration and Laurent-depth accounting. Coefficient/normalization poles determine the necessary master-integral depth. D-dimensional algebra closes before scalar dimension substitution, preserving evanescent finite contributions.

Collision-safe scalar aliases protect completed dot products during analytic spin expansion while exposing open tensor slots. Powered-dot freshening also remains necessary at the whole-numerator FORM boundary: the independent average of `(k·p+k·q)^2` is `k^2*(p^2+2*p·q+q^2)/D`, with `D=4-2*eps`. Forest composition retains physical mass owners: contained owners participate in Taylor expansion, disjoint owners remain fixed and partial overlaps are rejected.

=== 7. Canonical local UV projection and reuse
<canonical-local-uv-projection-and-reuse>
Raw sectors retain physical owners and provenance for subsequent Taylor operations. A separate canonical projection view groups exact signed routing, mass, denominator polynomial/prescription and component domain. Physical incidence remains an independent witness. Only certified serial paths or pure cycles permit equivalent channels to merge into powered quotients; nonadjacent channels retain their topology.

Completed hard factors may use an unambiguous certified class. Physical-source-fixed/soft factors keep their restrictions; hard factors without a surviving class remain affine carriers. Odd numerator signs are independent of denominator evenness. Complete numerator and denominator reconstruction certificates precede native CFF generation.

The compressed allocator admits a baseline, an eligible packed alternative and one placement challenger. Native source-map rows score the bounded proposals; rank and stable order break ties. This is not a global optimum. One nonserialized graph-owned context reuses immutable preparations, assignment templates and mapped subtrees under complete keys and bounded retention. Component waves and bulk cut-key addition combine complete coefficient values while preserving factorization.

Same-frame `+N/-N` sectors may cancel during grouping. The grouped result returns the existing typed zero before requiring production maps, preserving every supported final cut order. The full projection/build regression tests that value rather than sector layout.

=== Behavioral contracts and findings
<behavioral-contracts-and-findings>
Small scalar identities use exact subtraction and zero checks. Factorization-specific tests compare the complete retained factorized expression. Tensor oracles permit dummy-index renaming while preserving physical labels. Cache/build counts, storage layout, traversal/proposal order and private row-count goldens do not define numerical correctness; public retention bounds and supported value/error behavior remain valid contracts.

#figure(
  align(center)[#table(
    columns: (33.33%, 33.33%, 33.33%),
    align: (auto,auto,auto,),
    table.header([Finding], [Current contract], [Fresh execution evidence],),
    table.hline(),
    [Projection cancellation], [Complete typed zero at every supported final cut order], [PASS: C7/uv-performance-outward, C8/curated],
    [Private layout/cache assertions], [Full values, supported errors, factorization and public bounds], [Exact reviewed correction fingerprints in the source mapping; owning selections below],
    [Physical route tolerance], [Fixed 1e-9 complex-relative acceptance for GL00/GL01], [PASS: C7/physical-local-uv-routes, C8/physical-local-uv-routes],
    [Sunset allocation], [Independent complete degree-five/seven signed contour values], [PASS: C7/uv-performance-outward, C8/curated],
  )]
  , kind: table
  )

The signed sunset oracle covers all 126 distributions of degree five over five repeated occurrences and admitted degree-five/seven plans. Independent successive clockwise contours at energies `(2,3,7)` give `1213/16387080192` and `-365/47775744`. The retained mathematical audit is separate from fresh application-test execution. GL00/GL01 route comparisons require finite, nonzero complex results at three fixed points within `1e-9` relative complex-norm tolerance. Reported stability estimates cannot relax that threshold; integrated and threshold counterterms are disabled to isolate this comparison.

=== Rust patterns and KISS
<rust-patterns-and-kiss>
#figure(
  align(center)[#table(
    columns: (50%, 50%),
    align: (auto,auto,),
    table.header([Pattern], [Effectiveness and limit],),
    table.hline(),
    [Borrowed `AtomView`, owned `Atom`, concrete `Rational`], [Avoid copying during inspection, retain immutable witnesses and enforce the coefficient domain. Domain clarity does not prove speed.],
    [Typed references and domain enums], [Separate owners, canonical classes, occurrences, residue maps and valid input modes; raw/canonical UV types protect the recursion/projection boundary.],
    [Immutable `Arc` payloads in a graph-owned context], [Keep assignments and generated payloads together; scoped lifetime permits reuse without serializing caches into state.],
    [Existing builders and fallible conversions], [Keep validation and numerical preparation at established owners. Validate complete inputs before mutation or dispatch.],
    [`BTreeMap::entry` and consuming iterators], [Accumulate exact values and release consumed state. Deterministic order aids diagnostics without becoming an oracle.],
    [Fixed points and bounded proposals], [Structural convergence or decreasing state counts support termination; proposal limits do not bound total runtime.],
    [Result accessors and `Cow`], [Unify literal/lazy scalar interpretation and borrow when possible; unsupported interpretation returns an explicit error.],
    [Temporary scalar aliases], [Preserve completed scalar contractions during tensor expansion; collision checks prevent capture.],
  )]
  , kind: table
  )

Shared preparation, immutable bindings, one scalar-result boundary and one bulk accumulation path remove duplication. Canonical source reconstruction and tensor-boundary rewiring remain substantial algorithms with explicit invariants. Further cache layers and diagnostic state need representative evidence to justify their complexity. Strategy comparisons do not establish universal optimality.

=== Fresh validation
<fresh-validation>
#figure(
  align(center)[#table(
    columns: (33.33%, 33.33%, 33.33%),
    align: (auto,auto,auto,),
    table.header([Boundary], [Format / check / build / Clippy], [Owning behavioral selections],),
    table.hline(),
    [C1], [Four gates PASS], [3/3 selections PASS; 811 accepted executions],
    [C2], [Four gates PASS], [3/3 selections PASS; 182 accepted executions],
    [C3], [Four gates PASS], [1/1 selections PASS; 275 accepted executions],
    [C4], [Four gates PASS], [2/2 selections PASS; 623 accepted executions],
    [C5], [Four gates PASS], [4/4 selections PASS; 217 accepted executions],
    [C6], [Four gates PASS], [1/1 selections PASS; 242 accepted executions],
    [C7], [Four gates PASS], [2/2 selections PASS; 291 accepted executions],
    [C8], [Four gates PASS], [10/10 selections PASS; 2706 accepted executions],
  )]
  , kind: table
  )

#figure(
  align(center)[#table(
    columns: (33.33%, 33.33%, 33.33%),
    align: (auto,auto,auto,),
    table.header([Final selection], [Accepted passed / executed], [Status],),
    table.hline(),
    [Curated GammaLoop (includes phase/Vakint subsets)], [2135 / 2135], [PASS],
    [Full scalar matrix, including slow cases], [167 / 167], [PASS],
    [Vertex rules], [1 / 1], [PASS],
    [UFO model parity (API library)], [1 / 1], [PASS],
    [Physical local UV routes: GL00/GL01], [2 / 2], [PASS],
    [Spenso: shadowing], [168 / 168], [PASS],
    [Spenso: no default features], [50 / 50], [PASS],
    [Shared CFF: default features], [37 / 37], [PASS],
    [Shared CFF: all features], [108 / 108], [PASS],
    [Shared CFF: no default features], [37 / 37], [PASS],
  )]
  , kind: table
  )

Accepted final runtime: #strong[2706 executions over 2279 distinct binary/test identities];, from 10/10 passing selections. Accepted configurations report 1572 excluded entries, which may overlap. Baseline and owning-prefix executions are excluded from these final totals. Phase/Vakint cases inside the curated selection are not added again as independent suites. Configuration overlap and owning-prefix executions remain distinct in the ledger. Pending, running, failed and unaudited stages prevent a complete result; zero accepted failures does not establish overall success.

The recorded configuration uses `--locked` and the `test_gammaloop` nextest profile. Workspace check/build/Clippy gates and most suites use the `dev-optim` Cargo profile; the full scalar matrix uses `--release`. Gates include formatting, all-target workspace check/build and Clippy; separate shared-CFF/Spenso feature configurations have their own receipts. UFO parity selects the API library with `ufo_support`; this is package-scoped optional-feature coverage. Cargo, nextest and Rayon use four workers under the recorded 30 GB process-tree guard. Final runtime selections allow up to ten seconds per process query; the guard still stops the pipeline on monitoring failure. Snapshot updates and automatic test retries are disabled. The full scalar matrix selects its underlying namespaces directly, including slow cases. Exact arguments, warning policy, overrides and integrity guards remain authoritative in the ledger.

Dependency: Symbolica `2.2.0`, with `graphica` and `numerica`, from `alphal00p/symbolica` at `4d0a833eb8e059d1f95bdae5abed2559830b235f`. This is a specific locked snapshot. Tool versions: rustc 1.97.0 (2d8144b78 2026-07-07); cargo 1.97.0 (c980f4866 2026-06-30); cargo-nextest 0.9.140. License material is excluded from reports and archives.

=== Source preservation and review scope
<source-preservation-and-review-scope>
The mapping retains 37 remote and five review source attributions. Incoming `a1140c90c334ff58a3b040ff3eae06d3645bf0bb` is reviewed against `78395e3ab3ddd8d8f62b2f674d7488484eace197`: 68 paths (65 text, 3 binary), 792 text hunks, with explicit reviewer attribution. The three serialized Symbolica blobs received byte/hash and source-provenance review; their contents were not semantically decoded. At the tested candidate, exact correction fingerprints cover 21 changed files and the whole-stack inventory contains 647 commit/path records over 549 paths; inventory does not establish an unperformed fresh rereview of unchanged hunks.

Immutable corrected source `392c98c3d2f70896ef94267addfcf610e2df2c88` and reconciled C8 `449a59f6f173d100c2da858c541faa1e0dccca6e` share tree `c87cc95044663a68413010aceb81f9475be342c1` exactly. Later final-report changes require separate exact, reviewed documentation deltas from both the reconciled C8 and the executed candidate. Production/build/test/config content must remain identical. Candidate receipt hashes and execution counts are preserved; documentation closure does not count as another test run.

=== Limits and evidence
<limits-and-evidence>
GL00/GL01 automatic route selections do not certify GL262. The source-scoped GL262 diagnostic reaches a Symbolica evaluator-construction panic with compilation disabled and supplies no complete matched runtime/performance result. No confirmed standalone reproducer of that large panic is established. The nonsymmetric `f(x,y)` import MRE demonstrates a separate argument-order defect and does not identify the panic's cause.

Fresh standalone evidence at executed revision `a3af927afc633c6e6c2474d3f03fa02e0242eaa4` has locked all-target check/build/Clippy passes, #strong[two positive controls (39 and 4)];, and #strong[one reproduced import defect];: the writer exits 0 and the reader exits the expected 101 for `f(y,x)` versus `f(x,y)`. These are not workspace test counts or a GL262 replay. Source-equivalent standalone package; no additional execution. The #link("raised-energy-cff-mre-validation.typ")[standalone receipt] records exact package/manifest/lock/config identity and original binary/result hashes; carrying that evidence adds no execution. Separately recorded formatting passes at `bbb8ff1db499a78b3b25c41f774c9a034acb5b4d`; its tracked package/build inputs match observed C7.

Manual physical PySecDec integrations remain outside automatic coverage. Installed-Python subprocess tests and other platforms have no separate acceptance unless an explicit receipt records them. Finite coordinate probes and Monte Carlo acceptances have numerical/statistical limits; local Ward checks do not establish complete electroweak LU acceptance.

The #link("raised-energy-cff-motivation-audit.typ")[motivation audit] and #link("local-4d-uv-performance.typ")[local UV measurements] retain workload/source/binary identities. On source `911f1824f994611d723dff0f50416cb301718617`, the local-only GL00/GL01 generation premiums are 4.45%/0.82%. A separate single-run scalar cohort includes integrated and threshold subtraction: GL21/base setup/generation takes 1.431860 times erased 3D. Other measured comparisons also show regressions; no general speedup is established. Those source-scoped observations do not certify this reconstructed stack. Stage timings, complete generation and saved-state sampling remain distinct.

Markdown and Typst are generated from the same report body, with non-whitespace visible-text parity checked. PDF compilation, every-page visual inspection, final bookmark and final document hashes belong to the external closure evidence produced after the final report update. This report embeds no hash of its own final content.
