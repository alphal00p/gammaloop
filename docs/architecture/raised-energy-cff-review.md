> Historical record for the revisions and dates stated below. Its test counts,
> API descriptions and continuation instructions do not certify or govern the
> reconstructed stack. See the [stack review](raised-energy-cff-stack-review.md)
> and [fresh validation record](raised-energy-cff-stack-review-validation.md),
> together with [current architecture](architecture-current.md) and
> [CONTRIBUTING.md](../../CONTRIBUTING.md). The original body is preserved.

**Review of `main::raised_energy_cff_wip`, completed 2026-09-08**

Reviewed base `395610143` through branch tip `91142139e`: 392 changed files,
76,394 added lines and 7,976 deleted lines. This continuation completes the
code-diff review, including tests, removed code, CLI/API, shared CFF generation,
UV reconstruction, runtime integration, Spenso/Vakint, and build tooling.
Every changed file is accounted for: 186 complete changed-hunk reviews,
205 generated-artifact audits and one binary inventory;
generated artifacts and the binary state map have separately stated audit limits.
Implementation and tests were left unchanged.

The branch has substantial, useful mathematical safeguards, and no incorrect
CFF value was reproduced for a valid physics input. Nevertheless, the tests do
**not** meet the requested outward-functionality-only standard. Several numerical
assertions can miss invalid or relatively large errors. The broader review also
reproduced CLI scope/state-output defects and diagnostic-library input problems.
The twelve findings below distinguish those concrete issues from design comments.
I would address the P2 findings before using this suite as a merge gate.

The Rust is often idiomatic locally. It is only partly KISS overall: typed
ownership and completion stages justify complexity, while duplicated builders,
divergent command preparation and correlated optional state add avoidable work.
Clean Clippy output does not settle either simplicity or test quality.

**Actionable findings**

1. **[P2] Nested inline runs expand variables before entering the inner scope.**
   [run.rs:268](/common/dev/gammaloop/higher-power-energies-review/crates/gammaloop-api/src/commands/run.rs:268)
   sends any placeholder-containing inline command through template expansion
   before recognizing a nested `run`. A direct invocation with its own
   `-D level=warn` and a body using `$(level)` succeeds; wrapping that same
   invocation inside `run -c` fails with “Missing command-block variable 'level'”.
   Both cases were run through the real CLI. Source tracing also shows an outer
   definition can capture the inner placeholder before the inner override applies.
   Prepare the inner invocation in its own environment, as the ordinary block
   path already does at line 477. One environment-aware preparation boundary
   would remove this divergence.

2. **[P2] Loading nested command blocks demands variables intended for invocation.**
   [session.rs:238](/common/dev/gammaloop/higher-power-energies-review/crates/gammaloop-api/src/session.rs:238)
   skips static validation only when a command's own raw text contains a
   placeholder. For `outer = ["run inner"]` and a parameterized `inner`, it
   recursively prepares `inner` with an empty environment before the top-level
   `run outer -D level=warn` can execute. A real boot-card reproduction exits 1;
   the equivalent single-block control exits 0. Reusable nested block libraries
   are affected even before execution. Keep static name/cycle validation, but
   defer invocation-variable resolution until the inherited environment exists.

3. **[P2] Default 3Drep output violates the explicit state-write policy.**
   [threedreps/mod.rs:476](/common/dev/gammaloop/higher-power-energies-review/crates/gammaloop-api/src/commands/threedreps/mod.rs:476)
   selects `<active-state>/threed_workspace` in ordinary mode, and lines 395–410
   save JSON and an output pointer by default. I imported a scalar box, ran
   `3Drep build --no-pretty`, and ended with `quit -n`; both files were created
   inside the state directory despite no save command. [CONTRIBUTING.md:369](/common/dev/gammaloop/higher-power-energies-review/CONTRIBUTING.md:369)
   reserves state-directory writes for explicit save/quit-output operations
   and enabled logs. Use the existing cwd-output policy consistently, or obtain
   a deliberate repository-policy change. This is a filesystem behavior issue,
   not a CFF-value discrepancy.

4. **[P2] The public diagnostic numerator parser silently changes arithmetic.**
   [eval.rs:756](/common/dev/gammaloop/higher-power-energies-review/crates/three-dimensional-reps/src/eval.rs:756)
   bypasses its right binding power for exponentiation;
   [eval.rs:777](/common/dev/gammaloop/higher-power-energies-review/crates/three-dimensional-reps/src/eval.rs:777)
   binds unary minus more tightly than powers; and
   [eval.rs:836](/common/dev/gammaloop/higher-power-energies-review/crates/three-dimensional-reps/src/eval.rs:836)
   narrows an exponent to `i32` without checking its range. A compiled public
   `evaluate_expression` reproduction on a generated unit graph gives
   `-2**2 = 4`, `-(2**2) = -4`, `2**3**2 = 64`, and
   `2**4294967296 = 1`. The first and chained-power cases contradict conventional
   power precedence; the last silently wraps the exponent to zero. Use checked
   conversion and define or reject chained powers explicitly. Tests should
   compare evaluated arithmetic at this public boundary. This parser belongs
   to the shared crate's diagnostic evaluation feature; the production compiled
   GammaLoop integrand does not use it.

5. **[P2] Route-comparison tolerances can accept order-one relative errors.**
   [scalar_3L_cross_section_inspects.rs:969](/common/dev/gammaloop/higher-power-energies-review/tests/tests/test_runs/scalar_3L_cross_section_inspects.rs:969)
   applies an absolute `1e-14` floor whenever both Arb results have norm below
   one. The alternative precision-scaled relative assertion consequently does
   not constrain small nonzero results: `actual = 1e-18` and `reference = 2e-18`
   pass despite a 50% normalized difference. The comment describes cancellation
   to zero, but the exception is applied to every case. The common f64 event
   and total comparisons also use `max(norm, 1)` with `1e-10`, so they do not
   restore sensitivity at small scales
   ([utils.rs:457](/common/dev/gammaloop/higher-power-energies-review/tests/tests/test_runs/utils.rs:457),
   [utils.rs:641](/common/dev/gammaloop/higher-power-energies-review/tests/tests/test_runs/utils.rs:641)).
   The same exception occurs in
   [test_gamma_star_ttx_nlo_acceptance.rs:315](/common/dev/gammaloop/higher-power-energies-review/tests/tests/test_gamma_star_ttx_nlo_acceptance.rs:315).
   Use a relative criterion for nonzero results at identical inputs. Handle
   certified zeros separately, with a justified absolute or propagated error
   bound. Supplying coordinates as f64 does not by itself justify a universal
   absolute error in the observable.

6. **[P2] Numerical failure predicates can accept NaN and infinity.**
   [utils.rs:635](/common/dev/gammaloop/higher-power-energies-review/tests/tests/test_runs/utils.rs:635)
   only panics when `distance > tolerance`. With a NaN total, the distance is
   NaN and the comparison is false. With an infinite total, both distance and
   tolerance can be infinite, so the comparison is false again. This can pass
   when the remaining event data agree, including empty event groups. The
   explicit finite checks in the scalar test are inside the optional
   `localized_3d_results` branch; they do not protect every caller path. Require
   finite components before an affirmative `distance <= tolerance` assertion.
   This also affects failure-list aggregation in the independent shared-crate
   tests ([eval.rs:2066](/common/dev/gammaloop/higher-power-energies-review/crates/three-dimensional-reps/src/eval.rs:2066),
   [generation.rs:6818](/common/dev/gammaloop/higher-power-energies-review/crates/three-dimensional-reps/src/generation.rs:6818))
   and Arb raised-LU comparisons
   ([cff/mod.rs:4418](/common/dev/gammaloop/higher-power-energies-review/crates/gammalooprs/src/cff/mod.rs:4418),
   repeated at 4510, 4626, 4692 and 4791). A compiled reproduction confirms that
   the actual `F<ArbPrec>` NaN compares neither greater than nor less than or
   equal to a finite tolerance. An empty failure list therefore does not prove
   valid agreement. No current production NaN at those sample points is claimed.
   Using `hypot` would also avoid avoidable overflow in the norm calculation.

7. **[P2] Multiple tests freeze private construction choices.**
   The following would fail under valid implementation simplifications:

   | Location | Incidental commitment | Outward replacement |
   | --- | --- | --- |
   | [generation.rs:10639](/common/dev/gammaloop/higher-power-energies-review/crates/three-dimensional-reps/src/generation.rs:10639) | A private component-product builder must return `None`. | Generate through the public entry point and compare the complete residue against an independent reference. |
   | [generation.rs:10912](/common/dev/gammaloop/higher-power-energies-review/crates/three-dimensional-reps/src/generation.rs:10912) | Exactly one variant per orientation and a branching denominator tree. | Compare the represented denominator/residue function after summing its representation. |
   | [generation.rs:11503](/common/dev/gammaloop/higher-power-energies-review/crates/three-dimensional-reps/src/generation.rs:11503) | The exact origin string `bounded_degree_known_factor_cff`. | Check the high-power reconstruction identity; an internal strategy name is not a physics result. |
   | [local_3d/tests.rs:625](/common/dev/gammaloop/higher-power-energies-review/crates/gammalooprs/src/uv/approx/local_3d/tests.rs:625) | A completed independent source sum is hosted at exactly `OrientationID(0)`. | Verify compatible mapping, exactly-once selection, and unchanged complete sum under a different valid host. |
   | [three_d_source.rs:6058](/common/dev/gammaloop/higher-power-energies-review/crates/gammalooprs/src/graph/three_d_source.rs:6058) | Exact edge-count sequence and cache population of two. | Compare cached/uncached values and sufficient capacity. Put an explicitly required cache resource bound in a separate performance check. |
   | [approx/mod.rs:2066](/common/dev/gammaloop/higher-power-energies-review/crates/gammalooprs/src/uv/approx/mod.rs:2066) | Exactly one sector/frame and particular loop-carrier IDs. | Check coordinate compatibility and complete reconstructed value across valid charts. |
   | [three_d_source.rs:3393](/common/dev/gammaloop/higher-power-energies-review/crates/gammalooprs/src/graph/three_d_source.rs:3393) | Synthetic node-name prefixes; nearby tests fix complete parsed graphs and canonical occurrence maps. | Check physical ownership, valid references and complete mapped values under relabeling. |
   | [cff/mod.rs:3122](/common/dev/gammaloop/higher-power-energies-review/crates/gammalooprs/src/cff/mod.rs:3122) | Zip two orientation lists and require equal individual maps and carriers. | Compare complete owner-invariant sums; separately test any promised public selector semantics. |
   | [bench.rs:1041](/common/dev/gammaloop/higher-power-energies-review/crates/gammaloop-api/src/commands/bench.rs:1041) | Field-by-field settings changes, including a private sentinel. | Execute benchmarking and subsequent evaluation; verify values and restoration on success and error. |
   | [state.rs:3618](/common/dev/gammaloop/higher-power-energies-review/crates/gammaloop-api/src/state.rs:3618), [process.rs:1608](/common/dev/gammaloop/higher-power-energies-review/crates/gammalooprs/src/processes/process.rs:1608) | Bincode equality of entire regenerated CFF, graph and cut-group structures. | Compare semantic exports, ownership and evaluated save/load results. Byte equality unnecessarily fixes representation. |
   | [network/tests.rs:311](/common/dev/gammaloop/higher-power-energies-review/crates/spenso/src/network/tests.rs:311), [parametric.rs:3347](/common/dev/gammaloop/higher-power-energies-review/crates/spenso/src/tensors/parametric.rs:3347) | Exact internal tensor storage and allocation choices. | Check contraction/evaluation and explicit resource budgets if those are required. |

   This does not mean every structural assertion is wrong. Original numerator
   ownership, denominator momentum/mass/multiplicity, and preserved factorization
   are explicit repository requirements. Test those semantic invariants without
   prescribing an arbitrary tree, cache layout, diagnostic label, or host.
   A deterministic result need not retain one particular lexicographic tie-break,
   and public threshold IDs need not be dense unless that is a documented API
   guarantee. Private-builder `None`, diagnostic origin names, host zero and
   storage layouts should not become indirect correctness requirements.
   Unit-test placement is also fine: a private module can still be tested by
   its observable result. Do not replace independent correctness tests with
   comparisons that call the same implementation twice.

8. **[P2] The GL04 denominator certificate does not inspect reconstructed denominators.**
   [three_d_source.rs:5208](/common/dev/gammaloop/higher-power-energies-review/crates/gammalooprs/src/graph/three_d_source.rs:5208)
   checks two owner-list lengths in a manually constructed mapper, then compares
   handwritten `D(q)^2 D(-q)^3` with handwritten `D(q)^5`. It never obtains
   denominators from the exact source builder. The assertion proves the
   evenness of that chosen denominator, but cannot detect a generated wrong
   mass, multiplicity or routing. The preceding numerator identity is useful
   within its narrower mapper scope. Derive one side from actual reconstructed
   source records and compare its complete denominator with the retained Taylor
   target. The earlier report overstated this test's denominator coverage;
   that claim is corrected here.

9. **[P2] The 198 added scalar snapshots are not exercised.**
   The new scalar snapshot files under
   [snapshots](/common/dev/gammaloop/higher-power-energies-review/tests/tests/test_runs/snapshots)
   have no snapshot assertion or file reader in
   [scalar_3L_cross_section_inspects.rs](/common/dev/gammaloop/higher-power-energies-review/tests/tests/test_runs/scalar_3L_cross_section_inspects.rs).
   The route loop ends with comparisons between generated routes at line 992.
   Editing those snapshot values cannot fail these tests. A common error shared
   by the compared routes can therefore pass even if it changes the recorded
   values. Restore a small, independently justified outward baseline where
   appropriate, or explicitly retain these files as historical artifacts rather
   than calling them regression coverage. These historical values were not
   treated as current numerical evidence in this review.

10. **[P3] The preserved-tree fast path skips energy-bound validation.**
   [generation.rs:668](/common/dev/gammaloop/higher-power-energies-review/crates/three-dimensional-reps/src/generation.rs:668)
   returns before checking the bound IDs at line 689. A public call with a
   one-edge, zero-loop `ParsedGraph`, preserved edge `[0]`, and bounds
   `[(999, 2)]` succeeds and retains those invalid source bounds. I reproduced
   this against the compiled library. `ParsedGraph` uses the default absent
   edge-index map, so the earlier optional remapping validation does not catch
   it. Validate all input bound IDs before the early return. This is an input
   validation inconsistency; no valid-input CFF algebra error is claimed.

11. **[P3] The public graph validator panics on malformed signature dimensions.**
   [validator.rs:39](/common/dev/gammaloop/higher-power-energies-review/crates/three-dimensional-reps/src/validator.rs:39)
   allocates vertex balances from declared loop/external counts, then indexes
   them using unchecked signature lengths at lines 52–53. A compiled public-API
   call with an edge signature `[1]` and no declared loops panics with a
   zero-length out-of-bounds index instead of returning an invalid
   `GraphValidation`. Public `ParsedGraph` fields allow this construction.
   Validate dimensions before balance accumulation and return a structured
   rejection. This concerns malformed external Rust input; no valid internally
   generated graph failure was demonstrated.

12. **[P3] Arb stability overflow renders an impossible accuracy bound.**
   [evaluation.rs:443](/common/dev/gammaloop/higher-power-energies-review/crates/gammalooprs/src/integrands/evaluation.rs:443)
   computes `10.0_f64.powf(-1000.0)` for the Arb histogram's overflow bound.
   That underflows to zero. An Arb stability result with estimated relative
   accuracy zero is placed in this overflow bin, and the displayed median
   becomes `<0.0e0` via lines 500–513. Keep the logarithmic exponent through
   formatting, or otherwise represent the bound without first converting it
   to an unrepresentable f64 number. This affects reporting, not the integrand.

**Functional architecture and its evidence**

The important pipeline comparison is between complete results for the same
physical cut:

| Route | Ordered stages | Correct comparison boundary |
| --- | --- | --- |
| Direct local 3D | Parsed physical graph and numerator → bounded CFF → Taylor operations on the complete CFF → keyed localization → evaluator. | Each direct key must have its required UV behavior; the complete cut includes all maps and raised-order derivative pieces. |
| Explicit direct 3D | Same direct construction → retain each keyed contribution once without runtime selectors → evaluator. | Selector-weighted sum must equal the explicit sum. |
| Projected local 4D | Retained source owners → 4D Taylor sectors → owner-preserving exact graph and immutable assignment → component CFFs → outer composition → evaluator. | First certify numerator and denominator reconstruction; then compare the complete assembled cut with direct 3D. |

The physical graph, intended numerator and completed cut are shared inputs;
local Taylor construction first differs at `Direct3dApproximation::run` versus
`Local4dCts` reconstruction and `Projected4dCts` composition. The shared CFF
engine and integrated add-back are downstream boundaries, so they should not
be blamed for a route difference until their actual inputs have been compared.
The review checked ownership/denominator contracts before assessing low-level
algebra. It did not establish a new valid-input route mismatch.

Different per-key decompositions or individual `lu_cut_order` pieces are not,
by themselves, evidence of a physical mismatch. The handoff correctly explains
redistribution between derivative slots. Conversely, agreement of just one
such piece is insufficient. Tests should reflect this boundary.

The exact-source mapper retains sign-sensitive numerator provenance, and the
generation metadata retains component-local prefactor conventions rather than
reconstructing them from final algebra. The explicit-sum evaluator rejects
individual orientation selection; the runtime group checks distinguish exact
residue-map catalogs from physical sign patterns. UV profiling clones its evaluator before fallible work, so errors do not
leave the process without its production integrand. It retains graph/cut identity,
transforms candidate LMB samples into generation coordinates, and evaluates the
summed limit even when individual orientations pass. Versioned state/workspace
loading rejects incompatible positional data before decoding. These are useful
functional safeguards, not gratuitous abstractions.

Strong existing outward oracles should remain:

- [eval.rs:3145](/common/dev/gammaloop/higher-power-energies-review/crates/three-dimensional-reps/src/eval.rs:3145)
  derives simple/double-pole residues independently and compares full CFF values.
- [eval.rs:3256](/common/dev/gammaloop/higher-power-energies-review/crates/three-dimensional-reps/src/eval.rs:3256)
  compares independent bubble contour moments across seeds and numerator forms.
  Keep the value assertions while relaxing incidental branch-layout assertions.
- [eval.rs:1255](/common/dev/gammaloop/higher-power-energies-review/crates/three-dimensional-reps/src/eval.rs:1255)
  checks a denominator-cancellation identity under a selected cut.
- [eval.rs:1730](/common/dev/gammaloop/higher-power-energies-review/crates/three-dimensional-reps/src/eval.rs:1730)
  checks invariance under a nonzero auxiliary sampling scale.
- [three_d_source.rs:4986](/common/dev/gammaloop/higher-power-energies-review/crates/gammalooprs/src/graph/three_d_source.rs:4986)
  checks a common-loop numerator identity for a manually constructed mapper.
  Its denominator section does not certify production reconstruction; see finding 8.
- [test_epem_a_ddx_nlo_acceptance.rs:30](/common/dev/gammaloop/higher-power-energies-review/tests/tests/test_epem_a_ddx_nlo_acceptance.rs:30)
  uses physical normalization and uncertainty requirements, supplying an oracle
  beyond agreement between implementations.
- [test_cli.rs:253](/common/dev/gammaloop/higher-power-energies-review/tests/tests/test_cli.rs:253)
  exercises command substitution, defaults, nested overrides, validation before
  partial execution, and persisted command history through outward behavior.

**Rust patterns and effectiveness**

| Pattern and example | Assessment |
| --- | --- |
| Sum types: `Local3DCts::{Direct, Projected4d}`, `Direct3dCts::{Root, Sectors}`. | Effective. Their variants correspond to real completion boundaries and constrain which operations make sense. Retain these distinctions. |
| Newtypes: `Local4dCts`, `Full4dCts`, `FinalIntegrands`, `OrientationID`, `CffGlobalPrefactorSign`. | Mostly effective. Completion stages, index namespaces, and sign parity are meaningful invariants. The sign type's small `product`/`factor` API is particularly concise. |
| Borrowed graph adapters: `ThreeDGraphSource` and `GraphThreeDSource<'a>`. | Effective separation: the shared crate owns CFF generation while GammaLoop owns graph parsing and physical provenance. Borrowing avoids transferring ownership of the full graph. |
| Immutable plan with owned mapping state: `EnergyPowerAssignmentPlan`, `ExactSourceEnergyMapper`, shared through `Arc` in `CFFTerm`. | Effective. One certified assignment controls generation and numerator mapping, and it can survive temporary reconstructed graphs without copying it into every orientation. |
| `Result`, `thiserror`, iterator `collect::<Result<_>>()`, `Option::transpose()`. | Idiomatic and useful at fallible graph/parameter boundaries. Optional integrated localization in `direct_3d/forest.rs:280` is a concise example. |
| Deterministic `BTreeMap`/`BTreeSet` keys and `Entry` aggregation. | Appropriate for reproducible capacities, ownership, and branch identity. Tests should assert determinism where promised without prescribing the particular deterministic ordering. |
| Explicit encoding of persistent metadata and omission of transient bounds. | Justified by the persisted source-frame contract. Custom encode/decode is more maintenance than derive, so semantic save/load evaluation is more valuable than exact encoded-byte equality. |
| Strategy options and small enums. | Generally clear, but booleans still permit invalid configuration pairs. Existing validation correctly enforces projected-4D requirements. Keep mode resolution at a single boundary as options grow. |
| Correlated `Option` fields in `GraphThreeDSource`. | Less effective: separate frame and coordinate-LMB options require a runtime both-or-neither check at `three_d_source.rs:869`. `Option<(LoopMomentumBasis, ExactUvSubLmbFrame)>` expresses that particular invariant directly without adding another helper type. |
| Generic candidate equivalence machinery in `energy_degree.rs:96`. | Unnecessary generality. The only `try_new<K>` caller supplies `K = ()`, so lines 141–164 always create one equivalence class. Direct validation of already-certified occurrence IDs would preserve behavior with less code. |
| Duplicated orchestration in `direct_3d/forest.rs`. | KISS weakness: `run` at line 287 and `run_local` at line 395 repeat root/sector conversion, coordinate-frame extension, and signed Taylor mapping. Consolidate through the existing owner and preserve its comments/invariants. |
| Free helpers and forwarding layers. | Mixed. Small mathematical functions are readable, but `contains_placeholder` constructs a whole set merely to test existence, and placeholder scanning is repeated in `placeholder_specs`/`expand`. Simplify within the existing module when changing this area; avoid adding another abstraction layer. |
| Public compatibility and unused code. | The no-op `serde = []` feature deliberately preserves downstream Cargo feature vocabulary; deleting it breaks callers selecting that feature. Its maintenance cost is small. The unused `default_active_state_output_path`, suppressed with `allow(dead_code)`, is a clearer deletion candidate. |
| Term aggregation: `CFFTerm` holds expression, orientation and mapper together. | Effective improvement over parallel expression/orientation vectors and a truncating `zip`; related data now travel as one value. |
| Typed indices with `TiVec`, including `TopologicalThresholdId` and `OrientationID`. | Effective: topology IDs, expression indices and residue-map keys have different meanings. Non-dense archive tests exercise that distinction through selected values. |
| Deferred functions and shared `preprocess_atom` lowering. | Effective: avoids expanding large symbolic expressions and reuses tensor lowering. Materialized-value and hyperdual tests are the right evidence. |
| `GraphImportSource`, borrowed `GraphCatalog`, and `GraphImportOptions`. | Idiomatic sum type, adapter and labeled argument object. `graph_name_by_id` can delegate to existing graph lookup instead of duplicating the ownership match. |
| `PreparedRun` with inherited ordered environments. | Good all-or-nothing preparation model, undermined by separate block, raw-template and inline preparation routes. The two scope bugs show the cost of duplicating the semantic boundary. |
| Hidden `CommandTemplate` plus optional raw text. | Weak: a sentinel command and separate optional payload allow a template without its text. A private parsed-or-template enum carrying its payload would express the state directly. |
| Repeated shared CFF builder operations. | `BoundedCffBuilder`, `KnownFactorCffBuilder` and `LowerSectorCffBuilder` duplicate surface copying, interning and variant insertion at generation.rs:3165, 4637 and 5751. Reuse an existing owner/interner for those operations instead of adding a fourth builder framework. |
| Repeated runtime evaluator strategy dispatch. | Amplitude/cross-section terms and their counterterms repeat deferred/explicit-sum/parametric mode selection. Put compatibility and construction in the existing evaluator owner so four sites cannot drift. |
| Rational coefficients stored as unrestricted `Atom` or canonical strings. | Less effective: surface.rs:69 repeatedly asserts rationality; expression.rs:576 and 593 stringify and reparse rational fusion keys. Keeping the exact value typed would remove invalid states and conversion plumbing. |
| Boxed lazy Cartesian-product iterator in lower-sector generation. | Reasonable generation-time tradeoff: limits intermediate storage. Dynamic dispatch alone is not a reason to replace this with a more elaborate generic design. |
| Fallible functions without an error path. | `KnownLinearExpr::mul_rational`, `rational_to_coefficient` and `scale_linear_energy_expr_rational` return `Result` unnecessarily. Direct values make the real fallible boundaries easier to see. |
| Exact rational basis search. | Clear and correct for small loop counts, but candidate row/column combinations can grow combinatorially. No performance conclusion is claimed without measurement. |

The branch is idiomatic in many local expressions but only partly KISS at the
architectural level. The useful complexity is ownership, stage distinction,
and exact algebra. The avoidable complexity is duplicated orchestration,
impossible optional states, unused generic machinery, and tests requiring
incidental layouts. Splitting a large module can help navigation, but it does
not itself remove this complexity. Prefer the concrete reductions above before
introducing more traits or builders.

The existing `series(...).unwrap()` calls in the direct Taylor kernel remain
robustness debt: symbolic expansion should ideally propagate contextual errors.
Those calls were present in the old local-3D path, so they are not reported as
new regressions. No unsafe-code defect was identified in the inspected changes.

**Suggested test contract**

Cover simple/repeated poles and scalar/quadratic/high-power numerators using
independent residues; compare complete cuts across all three UV routes; exercise
edge relabeling, routing reversal, alternative valid LMBs, equivalent numerator
factorizations, and nonzero sampling scales. Verify exactly-once selector
coverage and evaluated save/load parity. Require finite results and meaningful
relative error bounds, isolating certified-zero cases. Keep resource budgets in
explicit performance checks rather than exact internal cache counts.

Use complete emitted source records for reconstruction certificates. Keep
non-expansion, physical ownership and factorization tests because those are
explicit repository contracts. Relax exact generated names, origin tags, graph
layouts, variant counts and private cache populations unless each has a public
reason to remain fixed. Tests may live inside a module and still exercise an
outward result; moving every test into another crate would not itself improve
its oracle.

Exercise nested inline scopes and boot-card loading with invocation variables
through the CLI. The existing late-template fixture at
[test_cli.rs:1106](/common/dev/gammaloop/higher-power-energies-review/tests/tests/test_cli.rs:1106) merely retains a string
containing obsolete `bench --samples ... -c 1` syntax. Use a valid current
command and execute it after substitution. For benchmarks, check subsequent
user-visible evaluation after both successful and failed runs, rather than
mirroring every temporary settings assignment.

The large CFF test module repeats full-orientation summation, evaluator setup
and error calculations, including a roughly 700-line spectator test. Reuse
existing comparison/evaluation support and split independent outward cases into
named tests. Avoid a new general test framework or helpers that reproduce the
production decomposition. No tests were changed during this review.

**Documentation and generated-artifact consistency**

The architecture documents preserve useful mathematical contracts, but several
“current” statements have drifted.
[architecture-current.md:514](/common/dev/gammaloop/higher-power-energies-review/docs/architecture/architecture-current.md:514) still names state
schema version 4 while the API uses 5;
[kysvnqlq-rebase-review.md:76](/common/dev/gammaloop/higher-power-energies-review/docs/architecture/kysvnqlq-rebase-review.md:76) has readiness/A79 statements about a common-denominator DOD1 topology that conflict with the newer
separate natural-topology implementation. The run-history schema lacks the new
UV profile selection fields and enum. Reconcile current contracts and regenerate
schemas; retain historical reports with a clear historical scope. These are
consistency comments, not evidence of a numerical defect.

All 200 changed snapshots were audited as artifacts. The 198 scalar snapshots
have the coverage problem in finding 9. The two other snapshots are referenced
by assertions. JSON syntax and local schema references were checked; fixture
text and source references were inspected. The changed `state_map.bin` was
inventoried by size/hash, not semantically decoded. Historical acceptance values
were not independently rederived from their publications.

**Review coverage and executed validation**

The earlier scoped report did not establish whole-diff coverage. This continuation
read all changed code hunks, including removals, and all new code files in full.
The very large new shared generation/evaluation files and source mapper were
reviewed in full, including their tests. The coverage appendix records every
one of the 392 paths, reviewer coverage and the limits of generated/binary audits.
A file inventory alone is not treated as source inspection.

Actual validation performed on the reviewed branch:

| Check | Observed result |
| --- | --- |
| Formatting; shared all-feature check and Clippy across targets; API check/Clippy across targets, locked dependencies. | Passed. These are build/lint checks, not proof of mathematical correctness. |
| Shared representation library, all features. | 122 tests passed in the earlier licensed run. |
| Earlier focused GammaLoop run plus follow-up. | 68 tests passed, then 6 additional targeted tests passed. These overlap the broader run below and must not be added to it as distinct tests. |
| Broader CFF, source mapper, LMB, UV approximation and runtime filters. | 210 tests executed: 206 initially passed; four aborted on test-thread stack overflow. All four passed when rerun with `RUST_MIN_STACK=67108864`. No test edits. |
| CLI build and command parsing/template/bench/profile/completion tests. | CLI built successfully; 32 targeted API tests passed. |
| Targeted outward integration tests. | 7 tests passed: 3Drep CLI output/read-only behavior, repeated masses, raised-cut cancellation, benchmark restoration and stability histogram output. |
| Real CLI reproductions and controls. | Two scope failures reproduced with successful direct controls; unsaved 3Drep output reproduced. |
| Compiled public Rust reproductions. | Invalid preserved-tree bound accepted; malformed validator input panicked; diagnostic parser arithmetic changed; Arb NaN comparisons were unordered. |
| Python dependency lock consistency. | `uv lock --check --offline` passed (122 resolved packages). |
| Watchdog CLI scenarios. | Successful child, propagated child failure, memory-limit termination and invalid-limit rejection produced expected exits 0, 7, 137 and 2. |

The broad core filter skips 486 of the 696 discovered core tests, including
profile exclusions. The four stack failures are recorded as initial failures,
not hidden behind the successful rerun; the repository already uses larger
stacks for symbolic workloads. Parallel workers were used with the supplied
license. The secret was supplied only to child process environments and is not
in these documents or repository files.

The full workspace suite, complete scalar release matrix, NLO numerical campaigns,
and Python/standalone end-to-end campaigns were not rerun. Source inspection of
those tests is not equivalent to executing them. Exact commands, reproduction
inputs, logs and overlapping-test-count caveats are in the validation appendix.

Companion files: [validation evidence](/common/dev/gammaloop/higher-power-energies-review/docs/architecture/raised-energy-cff-review-validation.md)
and [file coverage ledger](/common/dev/gammaloop/higher-power-energies-review/docs/architecture/raised-energy-cff-review-coverage.md).
