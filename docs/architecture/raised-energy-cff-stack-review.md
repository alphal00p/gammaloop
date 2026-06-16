# Raised-energy CFF stack review

Functionality, Rust idioms, KISS, and behavioral test contracts · 11 September 2026

## Assessment

The stack is reconstructed as seven functional commits, with the documentation in jj change `vmrowquy` and `codex/raised-energy-cff-reviewed` as the final local bookmark. Inspection JSON preserves structured additional-weight identifiers, and the three compact-parser regressions check public execution results. Powered-dot freshening remains required by an independent FORM-boundary counterexample. **All seven reconstructed boundaries have passing certificates, and all ten final runtime selections pass at the amended source.** Complete path accounting and the execution evidence below refer to these reconstructed source revisions; final document and history review is recorded separately in the external closure receipt.

The implementation combines exact rational CFF coefficients, independent repeated occurrences, owner-specific integrated UV mass handling, shared contour normalization, and physical phase conventions. Complete signed contour oracles, separate numerator and denominator reconstruction certificates, and independent dimensional and physical references provide the strongest correctness evidence. Agreement between execution strategies is supporting evidence; a shared normalization can give both strategies the same phase error.

The [change-motivation audit](raised-energy-cff-motivation-audit.md) ([PDF](raised-energy-cff-motivation-audit.pdf)) covers 73 logical categories at its pinned benchmark source. It supplies controlled measurements on the physical gluon double triangle with a three-gluon vertex and top-quark loop, additional physical captures, and targeted inputs for otherwise inactive mechanisms. Its observations distinguish measured improvements, concrete defects, requested capabilities and unproved optimization benefits. Those measurements and source identities remain intact; they do not supply fresh validation or complete path coverage for the reconstructed stack.

Rust idiomaticity is good at the domain boundaries. Borrowed expressions, owned results, concrete rationals, typed IDs, explicit alternatives, and fallible conversions express useful invariants. The design keeps expression assembly, command preparation, source assignment, evaluator construction, and Laurent restoration in their existing owners. Source reconstruction and tensor algebra remain substantial because they represent distinct mathematical requirements.

## Functional commits

Each commit groups its signature changes, callers, settings, schemas and required tests under one functional owner, with no temporary compatibility layer. Formatting, locked workspace checks, all-target builds, Clippy and owning behavioral selections pass at every commit with its ancestors. C1 and C2 retain their exact unchanged-revision certificates; C3–C7 have fresh checks after the GL00 certificate amendment. All seven commits use `Lucien Huber <im@lcnbr.ch>` as author and committer.

| # | Commit | Scope |
| --- | --- | --- |
| 1 | `97eb06cb` | Symbolic and tensor foundations |
| 2 | `2ba79605` | Exact generalized CFF generation |
| 3 | `1f0d4fd6` | Raised-energy CFF and local UV reconstruction |
| 4 | `48104fa4` | Command execution and evaluation workflows |
| 5 | `32a57c5d` | Physical phases and model sewing |
| 6 | `ee798ce2` | D-dimensional integrated UV algebra |
| 7 | `vmrowquy` | Review, architecture, and validation documents |

The [source mapping](raised-energy-cff-stack-source-mapping.md) identifies the full functional commit hashes and source ownership. The first six commits pin the executable implementation. Documentation uses stable jj change `vmrowquy`; the final local bookmark is `codex/raised-energy-cff-reviewed`.

## 1. Symbolic and tensor foundations

Symbolic numerators retain their factorization through dummy allocation, tensor materialization, and contraction. Parser clones share fresh-index allocation. Explicit-name reservations and collision skipping prevent compact dummy names from capturing existing indices; shared allocation already existed at the parent boundary. Canonicalization preserves external and dual slots. Compact-vector materialization accepts supported, unambiguous axes. Scalar aliases and deferred sums are completed at the existing result boundary in both execution strategies.

`MinIntermediateCost` estimates symbolic copying and normalization with saturating arithmetic. It is a scheduling heuristic and does not change coefficients. Independent coordinate enumeration checks complete contractions across strategies, sparse and dense inputs, metric signs, and both operand orders. Dense–dense and mixed multi-index contractions already worked through the generic contraction API; the generalized grouped specialization changes their cost. The complete-contraction assertions permit different intermediate scores, storage variants and allocation indices. Three public parse/execute regressions compare a weighted Minkowski component sum for weights `-2`, `a` and `a+b`, and preserve complete opaque function values for two ambiguous tensor arguments. They do not constrain the private materializer or its allocated slots.

The symbolic dependency and Vakint prerequisites support exact imaginary constants and namespace-preserving FORM round trips. Powered scalar-dot copies require independent dummy indices at the whole-numerator FORM boundary: otherwise `(k·p+k·q)^2` loses `1/D` from both diagonal angular-average terms. Commit 1 includes a public tensor-reduction regression against `k^2*(p^2+2*p·q+q^2)/D`, with `D=4-2*eps`, using the whole-numerator API available there. Rust scalar round trips alone do not expose this error. The [powered-dot follow-up](raised-energy-cff-powered-dot-freshening.md) retains the independent counterexample; the measured conversion overhead serves a correctness requirement. Both old sparse fibers and the generalized kernel can visit absent coordinates. The motivation audit measures all six exposed scheduling choices, storage kernels, aliases and deferred sums on physical captures and controlled stress inputs.

## 2. Exact generalized CFF generation

The shared engine keeps repeated occurrences independent, including when physical energies coincide. Numerator capacity belongs to an occurrence; cache reuse cannot redistribute it among equal-energy copies. The four `LinearEnergyExpr` coefficient domains and `CFFVariant::prefactor` use `Rational`. Conversion to symbolic expressions occurs through `Atom::num` at output boundaries.

`ExpressionAssembler` owns surface interning, source copying, variant insertion, and label finalization. Cache keys include source topology and relevant contour and capacity state. Complete residue maps determine which contributions can be combined. The terminal builder enforces its rank requirement directly.

Public APIs validate compact edge IDs, endpoints, signature dimensions, and cut aliases before indexing. Balance accumulation uses i64 for sums of valid i32 routes. Numerical evaluation rejects incompatible shapes, zero uniform scale, and unsupported integer powers. Arithmetic parsing preserves precedence, right association, unary signs, and checked exponents. Empty trees and trees with infinite denominators obey their symbolic zero contract; fusion preserves zero contributions.

The tests cover signed simple and repeated poles, connected and disconnected components, both occurrence owners, reversed routing, lower sectors, scale invariance, and execution strategies. Tree transformations compare complete denominator algebra; selectors use an explicit truth table; cut checks compare complete weighted sums. Fresh default, all-feature and no-default-feature test selections pass; all-feature and no-default-feature all-target checks also pass from commit 2 onward.

## 3. Raised-energy CFF and local UV reconstruction

GammaLoop owns physical source construction, denominator multiplicity, energy assignments, UV orchestration, and evaluator preparation. The shared engine owns exact CFF algebra. Physical edges, exact occurrences, production residue maps, and topological thresholds have distinct IDs. Physical source energies determine capacity; loop-momentum-basis coordinates route momenta without defining energy identity.

One immutable assignment plan controls generation capacity and numerator mapping. Original factors keep their source occurrences, and derived factors use only eligible copies. Exact reconstruction retains numerator signs under routing reversal. GL00 and GL04 child certificates independently check complete factorized common-chart numerators and denominator momenta, multiplicities, masses, physical owners and UV domains. The GL00 edge-5 temporal case uses an independent signed T0+T1 oracle; the GL04 case covers T0+T1+T2. Both retain the vacuum mass in full numerator expressions and the expansion mass in denominator metadata until the documented final mass identification. These checks start from the completed post-Taylor child source and do not certify every forest or the Taylor-producing stage. Separate GL04 source evidence also includes the cograph domain. Equality of denominators alone cannot establish the sign of an odd numerator.

Direct local 3D Taylor operations act on the complete generated CFF. Projected local 4D reconstruction uses completed raised-denominator terms and their complete source sums. The routes share evaluator assembly, production orientation catalogs, runtime validation, and exact source conversion. Grouped orientation sampling rejects incompatible map catalogs.

Integrated finite coefficients carry their physical owner through symbolic multiplication. A Taylor operation scales the vacuum mass of a contained owner, leaves a disjoint owner fixed, and rejects partial overlaps. The normalized localization kernel stays fixed. Output copies restore the physical mass symbol while stored sectors retain ownership for subsequent forest operations.

Behavioral tests check complete nested and disjoint UV functions, both local bubble contours plus one analytic finite addback, and exported forest values. DOT exports are re-imported and aggregated by physical forest and residue identity. Renormalization checks use signed Laurent expressions and a scale derivative. The nested vacuum fixture declares unit loop normalization explicitly, and the quark oracle aligns the external color basis using an independent color identity.

Saved-state manifest version 7 describes rational CFF storage. Incompatible manifests are rejected before decoding or overwriting a payload. Named generation records shared settings and requires consistent provenance for generated siblings. Persistence tests exercise rejected-operation preservation, public exports, and save/load behavior. Arbitrary-precision baseline expectations follow decimal promotion exactly.

## 4. Command execution and evaluation workflows

Command preparation uses one placeholder parser and one environment-aware boundary. Nested blocks inherit lexical variables; local definitions override them. Template text is an enum payload, so reusable blocks can contain variables supplied at invocation. Preparation validates names and cycles and restores active-block tracking on returned errors.

Graph import options collect source, process, and naming policy at one boundary. Inline DOT is a literal argument. Explicit process specifications must agree with the graph kind and existing process definition. Rust schemas describe the serialized string-or-command contract, and generated schema references resolve.

Read-only and explicit state-write policies cover evaluation, inspection, benchmarks, profiling, and 3D outputs. Diagnostic 3D output defaults to a workspace in the current directory. CLI failures produce a failing process exit status. Resume distinguishes absent state from corrupt or unreadable state and checks workspace versions. Benchmark settings are restored after success and returned errors; panic recovery is outside the contract.

Retained additional weights serialize as `[key, value]` pairs, including `[]` for no weights. `GenericAdditionalWeightInfo.weights` reuses the existing `vectorize` Serde adapter; an explicit `T: Deserialize<'de>` bound retains the original deserialization requirement without strengthening the generic type. Structured threshold-counterterm identifiers remain values, so `inspect --json-output` can preserve them. JSON/bincode round-trip tests cover all four key variants and empty weights. A physical inspection regression compares numerical values, retained weights and grouping with retention disabled and enabled; it does not claim every event metadata field. The separate fresh GL04 replay checks complete retained event kinematics, cut metadata, weights and grouping under this toggle, exercises structured threshold-counterterm keys in JSON, and verifies the saved state stays unchanged. This replay supplies serialization and retention evidence; the independent signed contours and physical oracles described elsewhere supply physics evidence.

Stability output retains batch statistics until assembly. Median formatting preserves a decimal exponent without underflowing an arbitrary-precision value through f64. Timing rows distinguish evaluator time from inclusive integrand time. Long benchmark requests retain a full batch and can consume substantial memory.

Tests execute nested commands and inspect saved histories, public state, exported and reloaded functions, output paths, restoration, and errors. The 3D CLI regression evaluates its written artifact. UV fail-fast checks its reported failure against exhaustive results without choosing traversal order. Public completion counts remain valid contracts.

## 5. Physical phases and model sewing

Physical normalization retains complete UFO propagator and vertex factors and applies common conversion at the cut-group boundary for ordinary and counterterm terms. Connected right-hand components include cut hairs, preserving contact interactions. Marked inverse-process vertices supply Hermitian partner couplings and spin matrices exactly once. Left/right threshold signs are checked against fixed pole-distribution derivatives and signed contour references.

Initial-state subtraction removes an endpoint only when it actually touches the initial-state cut. Fixed-momentum internal bridges retain their vertices and loop structure. The endpoint search returns `Option`, expressing the valid absence of such an endpoint. The worker-count regression exercises this production behavior.

Gauge sewing uses model-declared covariant quartets with validated spin, statistics, charge, color, mass, propagator, and conjugation contracts. It closes complete final-state multisets before filtering graphs and reports physical event labels. Independent polarization sums in rest and boosted frames check the complete vector, Goldstone, and ghost result. Matching numerical masses alone does not identify gauge partners.

One anticommuting-particle predicate covers fermions and ghosts. Optional CP symmetrization is an explicit user assumption in generated data. Mass checks use current model values and the used propagator denominators. Complex-mass cutting rules and arbitrary non-Hermitian interactions are outside the supported contract.

Numerator grouping uses short global lookups and atomic search/insertion within topology buckets. Exact symbolic matching precedes sampled fallback. Worker-count tests compare complete values within each grouping strategy. Different strategies may select different momentum routing representatives, so pointwise equality between strategies is not a general contract. Sampled fallback remains probabilistic.

Signed amplitudes, phase space, flux, bubble discontinuity, virtual mirrors, spin matrices, Ward identities, and LO/NLO acceptances provide complementary oracles. Acceptances constrain real and imaginary parts. Graphwise scale checks require the analytic sign explicitly. Seven vertex-rule fixtures compare complete signed tensors after dummy-index alignment.

## 6. D-dimensional integrated UV algebra

Projected-tensor and complete-numerator Vakint inputs share Lorentz validation, coefficient restoration, and Laurent truncation. Epsilon-dependent and complex coefficients retain poles, finite terms, and evanescent effects. Kernel depth accounts for coefficient and normalization poles before truncation; a factor such as `1/(4-2ε)` must participate while the required scalar-master orders are available.

Commit 6 extends the complete powered-scalar angular-average regression from commit 1 to both `project_onto_tensor_integrals` settings. The same independent `1/D` oracle checks both routes when the two-mode API becomes available, with no compatibility helper in the earlier commit.

The GammaLoop boundary restores projected slots, completes supported d-dimensional gamma algebra, substitutes scalar dimensions, and requests Laurent coefficients. Open representation slots retain symbolic dimensions. Unsupported tensor operators and traces return errors. Input mode and finite depth propagate consistently through cut forests, hedge-poset forests, and MUV renormalization.

Numerical settings preserve signed real momentum components and reject nonreal momentum inputs. General complex scalar coefficients have a distinct contract from pole masses. Fallible parameter collection completes before settings mutation. The numerical route combines sectors into one kernel so its uncertainty describes the complete integral.

Tests check full Laurent values, signed Minkowski components, alpha-renamed contractions, malformed-input errors, and mode parity. Reconstructed expressions and public Laurent accessors permit omitted or stored exact zeros while rejecting unexpected nonzero powers. Independent Dirac identities and scalar-master oracles cover nonzero evanescent finite terms. Public adapters exercise both production modes; injected analytic kernels provide additional bounded algebra checks.

## Behavioral test contracts

The reviewed tests assert complete values, errors, physical ownership, exported functions, and state restoration. For small symbolic expressions, equality uses the complete expanded difference and requires zero. Factorization-specific assertions use direct equality. Tensor comparisons align contracted dummy indices while preserving physical labels. Independent oracle values, signs, and tolerances remain explicit.

Private cache counts, allocation indices, planner ordering, storage layouts, and incidental expression formatting are outside these contracts. Deliberate public display output, physical coordinate identities, and required factorization remain observable behavior. Fixed fixture counts ensure zipped comparisons cover every expected item. Some unchanged repository tests still inspect internals; the review does not certify every test in the repository as implementation-independent.

## Rust patterns and KISS

| Pattern | Effectiveness and tradeoff |
| --- | --- |
| Borrowed `AtomView`, owned `Atom` | Inspection avoids copying expression trees; completed values and cache keys own their data. Clones are justified where source provenance must survive mutation. |
| Concrete `Rational` and checked conversion | The type enforces the coefficient domain. Explicit symbolic and precision boundaries, checked dimensions, exponents, and Laurent depth prevent silent narrowing. |
| Typed IDs and domain enums | Owners, occurrences, residue maps, thresholds, parser values, and UV routes stay distinct. These types do not validate arbitrary deserialized indices by themselves. |
| Paired optional data and enum payloads | An exact frame travels with its basis, and template text belongs to its template variant. Correlated state is represented together. |
| Existing builders and immutable plans | Expression assembly, source assignment, function construction, and evaluator preparation each have an owner responsible for their invariants. |
| Scoped shared state | `Rc<Cell<_>>` and `Rc<RefCell<_>>` fit the recursive single-threaded parser; `Arc<Mutex<_>>` fits parallel topology grouping. Parser borrows are scoped, and parallel lookups release the global map lock before bucket comparisons. |
| `Result`, `?`, `transpose`, fallible collection | Absence remains distinct from failure. Exposed malformed input gets contextual errors, and collect-before-mutate protects settings. Internal invariant panics are not a universal recovery boundary. |
| Iterators and explicit loops | Lazy component products limit retained Cartesian state. Explicit loops remain clear for stateful contractions and polynomial accumulation. Local boxed dispatch handles changing component depth. |
| Semantic cache keys and deterministic maps | Reuse depends on graph, frame, and capacity state. Ordered maps provide reproducible inventories; behavioral tests allow storage strategies to change. |
| Const generics and saturating heuristics | Contraction policy reuses established rules, and overflow cannot make expensive candidates look cheap. The heuristic provides no global optimality or physical memory-bound guarantee. |
| Serde field adapter and explicit bound | Existing pair-list serialization preserves structured keys without a string grammar or custom helper. An explicit deserialization bound restores inference suppressed by the field adapter and keeps the public generic type unchanged. |

The KISS strengths are shared preparation and assembly paths, existing builders, direct invariant checks, and narrowly scoped ownership. No generic testing framework or compatibility layer is needed. Specialized source/UV algebra and independent certificates remain readability costs. Their owner, sign, denominator, and dimensional distinctions are mathematically necessary. The motivation audit does not demonstrate a benefit from compaction on its reached single-component product or from extra losing assignment proposals on its completed GL04 inputs. These remain bounded simplification candidates.

## Validation

**All seven reconstructed boundaries and all ten final runtime selections pass.** The executable implementation is pinned by commits 1–6, ending at `ee798ce279f280000537e4d8c05ecc83fc6ec2cf`. Final runtime ran at `9e260684581ab2893eff4032c5ca8ea1fc2c18cb` (tree `9cebc793fe813699b9d689f1ab245715902ac04e`). Documentation remains in stable jj change `vmrowquy`; the external closure receipt records its final content/render review and exact source-to-bookmark checks.

The [validation reference](raised-energy-cff-stack-review-validation.md) records exact fresh commands, selected identities, counts, configurations, diagnostics and warnings. The [closeout evidence archive](raised-energy-cff-closeout-evidence.tar.gz) and [archive receipt](raised-energy-cff-closeout-evidence.json) retain the underlying records. The retained initial-candidate runtime and motivation-audit measurements are not counted as amended-source runtime executions. C1/C2 reuse exactly matching revision certificates; C3–C7 are revalidated after the GL00 amendment.

| Boundary checks | Passed / recorded | Guard time |
| --- | --- | --- |
| Formatting, locked workspace check, all-target build and Clippy | 28 / 28 | 4727.474 s |
| Required shared-CFF and optional API feature checks | 14 / 14 | 34.207 s |

Guard times include complete commands and build work. The following times are nextest test-execution times and are not end-to-end benchmark timings.

| Owning boundary | Passed / selected executions | Test time |
| --- | --- | --- |
| Commit 1 | 693 / 693 | 9.532 s |
| Commit 2 | 182 / 182 | 31.653 s |
| Commit 3 | 401 / 401 | 173.215 s |
| Commit 4 | 805 / 805 | 473.811 s |
| Commit 5 | 399 / 399 | 183.766 s |
| Commit 6 | 424 / 424 | 105.744 s |

C1–C6 owning selections account for 2904 passing executions, zero failing executions and 6255 reported skipped entries. C1/C2 figures describe their retained unchanged-revision certificates, while C3–C6 figures describe fresh amended-boundary executions. Repeated boundaries and feature configurations remain separate executions. C7’s behavioral records are the same ten final runtime selections below and are counted once.

| Final runtime selection | Passed / selected | Test time |
| --- | --- | --- |
| Curated GammaLoop | 2065 / 2065 | 465.870 s |
| Phase regressions | 212 / 212 | 30.817 s |
| Signed acceptances | 4 / 4 | 120.076 s |
| Vakint / UV | 241 / 241 | 73.977 s |
| Full scalar matrix, including slow cases | 166 / 166 | 922.523 s |
| Shared CFF: default | 37 / 37 | 5.191 s |
| Shared CFF: all features | 108 / 108 | 21.324 s |
| Shared CFF: no default features | 37 / 37 | 5.116 s |
| UFO model parity | 1 / 1 | 0.620 s |
| Vertex rules | 1 / 1 | 0.356 s |

Final runtime: **2872 passing executions across 2207 distinct binary/test identities**, zero failing executions, 3331 reported skipped entries, and 8 recorded Cargo/environment configurations. Distinct identities are binary/test pairs across all runs; configurations and repeated executions are reported separately. Skipped entries are outside the selected run and are not passing executions. Overlapping selections and feature configurations are counted separately as executions.

The full scalar matrix includes 26 standard and 140 slow cases (166 total), selected directly in release mode with ignored tests enabled. Licensed tests run with four workers, zero retries, disabled snapshot updates and a 30 GB process-tree RSS limit; Cargo uses four build jobs. The scalar command does not use the wrapper’s forced single-test setting. Clippy and the explicitly warning-tolerant curated/scalar commands pass with compiler/dependency warnings retained in the validation reference; passing exits do not mean warning-free compilation.

The fresh GL04 replay completes four inspection calls: plain and JSON output with additional weights disabled and enabled. Its frozen CLI SHA-256 `53f80aa6ca611a8566335c488ffb028d73d3d1ed3057e741a750edefa9aff0a4` exactly matches the binary copied from the successful C7 build. Disabled output matches the retained successful evaluation; the toggle preserves complete non-additional-weight content, including every event’s kinematics, cut metadata and complex weight, with group membership and multiplicity retained. Every emitted additional-weight entry is a well-formed key/value pair, and the output contains a structured `ThresholdCounterterm`; the replay does not independently determine the expected physical key set. Only documented timing metadata and the empty-object-to-pair-list weight encoding are normalized. The saved state and binary hashes are unchanged. This is serialization/retention evidence, independent of the signed-contour and analytic physics oracles; it is not a performance measurement.

Coverage accounts for **592 commit/path records and 518 distinct paths** at inventory revision `197b29f67ae065d2cf4d5d88b7658c29a0712af6`. C1–C6 coverage consists of 558 records over 485 paths. Relative to the retained initial review, 554 exact before/after transitions reuse their existing review, 2 owning architecture patches replay exactly on the amended parent, and 2 C3 source/document records have fresh GL00 review. Earlier fresh-hunk and original-patch review methods remain linked in that retained lineage. C7 contains 34 documentation/evidence paths. Final C7 content and render review is recorded separately in the external closure receipt; inventory alone is not semantic review. These totals come from complete immutable-tree inventories, including the new evidence files. The motivation audit does not certify new paths, and coverage does not mean every unchanged repository line was reread.

## Scope and remaining limitations

Manual physical PySecDec integrations are outside automatic coverage. Other failing-class diagnostics and `aa_aa::important::aa_aa_local_inspect_backend_consistency` are excluded; the vertex-rule case is selected explicitly. The three documented triangle, double-triangle, and three-loop TBT decimal fixture points are outside the automatic suite, and their derivation is unverified.

The retained UFO/JSON model-asset audit covers 22 Lorentz structures, 108 couplings, 153 vertices, 72 parameter declarations, 43 particle contracts, three quartets and three vector propagators, with five deterministic substitutions for sampled expression checks. All 25 audited asset entries retain their modes and blobs. The fresh `fresh_sm_ufo_preserves_virtual_gauge_and_covariant_cut_states` regression checks covariant multiplets and twelve named W/Z, Goldstone and ghost propagator/state contracts. Local electroweak and Ward certificates do not establish a complete electroweak LU acceptance. Finite routing probes, sampled grouping, and Monte Carlo acceptances have mathematical or statistical limits.

The diagnostic f64 evaluator retains its documented singular-surface behavior. General serialized-graph validation, complex-mass cutting rules, and arbitrary non-Hermitian interactions are outside the supported guarantees. The motivation audit bounds its contraction results to the measured dimensions, densities and physical captures; extremely sparse large domains remain unmeasured. Large benchmark batches still require workload-specific resource bounds. Inspection JSON, the three compact-parser execution regressions and the powered-dot angular-average regression are included in their fresh owning selections. The GL04 replay certifies serialization and retention at its recorded point; it provides no performance measurement or independent physics oracle.