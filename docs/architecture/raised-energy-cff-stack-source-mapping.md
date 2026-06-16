# Source-to-commit mapping

Seven functional commits and one documentation commit retain the original 42 source attributions and integrate the fully reviewed incoming change. Counts below are generated from the pinned Git objects; execution outcomes are recorded separately in the [validation report](raised-energy-cff-stack-review-validation.md).

| Commit | Stable jj change | Observed revision | Responsibility |
| --- | --- | --- | --- |
| C1 | `pmwlrxpkzyzsoulmpvnqyxsxvnkomwor` | `665658b168981c5c508e6cd51d3b57bbb55a27fb` | Symbolic/tensor foundations and graph bookkeeping |
| C2 | `rslqmrwvpzykwysxopvzrosppxtwrzno` | `a9397ee0f86b115fdfab50e81c48dde2137692e7` | Exact shared generalized CFF generation |
| C3 | `nysmrrzntwywmrppouuyssvyqourzwnx` | `3ea313789a1a5290093579febd7d1c601b92594e` | GammaLoop adapters, local UV reconstruction and evaluator preparation |
| C4 | `qpmwptwkyxnrtyoqtuvozvnyvmlxlwvq` | `04637884f24fb28a247574f04e03308cf4c04afb` | Command, persistence, evaluation and benchmark workflows |
| C5 | `lkqlqpmsqyokutpqrvvmzxossmozsxkv` | `1f2cf6d8236de91138951af8380bf1fd1eb4e783` | Physical phases and model sewing |
| C6 | `otnumvznxsrzplxvoqorwlyoolyurzwr` | `8e4e4664f3f9303760eb18d7972b91ff867b12ea` | D-dimensional integrated UV algebra and analytic scalar-product protection |
| C7 | `vylkrvskwvkuqxozltmsnztvkqpwpmlp` | `6d7a9220129c0fac37bdc0e04ebefccb1688c7be` | Canonical UV algebra, bounded occurrence allocation and graph-owned reuse |
| C8 | `oxsmszsroxrtxtzpowrowplwlqnowzyn` | Stable change; executed candidate is recorded in validation | Current reports, source attribution and validation evidence |

C8 is identified by its stable jj change. Executed candidate revisions remain in validation receipts. A final documentation revision may reuse those results only through an explicit certificate of identical production/build/test/config content and the exact allowed document delta. Final bookmark and rendered-document closure remain external.

The immutable pre-document source `392c98c3d2f70896ef94267addfcf610e2df2c88` and reconciled C8 `449a59f6f173d100c2da858c541faa1e0dccca6e` share tree `c87cc95044663a68413010aceb81f9475be342c1`: **verified**. This is the exact split-preservation proof before final reports. Later documentation changes require a separate reviewed delta from this reconciled C8 and from the executed test candidate; the final full tree need not equal this pre-document tree. Untouched input duplication/squashing is recorded separately in `squash-tree-preservation.json`.

The observed eight boundaries contain **647 commit/path records and 549 distinct paths**. Exact before/after Git blobs appear in `raised-energy-cff-stack-current-mapping.json`. These counts are an inventory, not a fresh review claim for unchanged earlier hunks.

The preserved reconciliation chain retains each immutable source/tree pair and its explicitly reviewed source followup.

| Checkpoint | Corrected source | Reconciled C8 | Equal tree |
| --- | --- | --- | --- |
| 1 | `d0242b4714a5a6aec54bdb38fcf46b82f26e325a` | `a1cf5b3f0ee68b065225e6a18a264a333a0006ab` | `c74a276735c3142a6db8cafe1bbd08127a6a7681` |
| 2 | `3a917074a0278f8b275f493b0728ca50dc288ea7` | `b15d068f52a853f88a1dac2818e30af44fbcbc9c` | `0ff2af98e5b10287282f8552e405afb6c2659270` |
| 3 | `8b2b355a76f563cd68dc0a3424bacba10e1eb865` | `1bc2316109c415e9b08020fb836c614847f5d926` | `792c75a898c7b60f92ce5446256153c3c4ea7e88` |
| 4 | `392c98c3d2f70896ef94267addfcf610e2df2c88` | `449a59f6f173d100c2da858c541faa1e0dccca6e` | `c87cc95044663a68413010aceb81f9475be342c1` |

## Original source attribution

All 37 remote and five review Git objects retain their recorded trees. Their 42 source descriptions and attributions are preserved below, with the former documentation owner 7 mapped to owner 8. The empty remote marker contributes no payload; no empty review intermediary is included.

| # | Source | Kind | Source description | Owners | Functional disposition |
| --- | --- | --- | --- | --- | --- |
| 1 | `42ce9317` | remote | ignore split feature outputs | 1 | Generated split outputs are ignored by foundational tooling. |
| 2 | `e4dbd169` | remote | relax complex evaluation-domain bound | 1 | Complex EvaluationDomain bound relaxation stays with tensor foundations. |
| 3 | `0ee54199` | remote | add CLI state output helpers | 4 | CLI state/output helpers belong to command execution. |
| 4 | `8a3512ed` | remote | add three-dimensional-reps core crate | 2 | Introduces the complete shared engine and its dependency metadata. |
| 5 | `939af81f` | remote | stage raised-energy CFF integration | 1, 2, 3 | Workspace metadata is foundational; shared-generator changes belong to 2; core adapters, source lifting, evaluator consumers and outward regressions belong to 3. Later workflow/phase edits to the same paths do not move these original hunks. |
| 6 | `bb60b02f` | remote | add command-template placeholder support | 4 | Command-template parsing and late argument expansion. |
| 7 | `49999ca0` | remote | parameterize command blocks | 4 | Command-block scoping, completion, persistence and CLI tests. |
| 8 | `b58d1ef0` | remote | update CLI generation and state plumbing | 3, 4 | Core generation consumers required by 3 stay together; general CLI/import/state behavior belongs to 4. |
| 9 | `a4ea80a9` | remote | add stability and standalone benchmark QoL | 4 | Benchmarks, stability summaries, dashboard and evaluation output. |
| 10 | `cf36d272` | remote | add CFF-only 3Drep CLI | 4 | Public 3Drep commands, API dependency and command-level tests. |
| 11 | `28b16c61` | remote | add repeated-mass CFF regressions | 3, 4 | Four graph fixtures are introduced with core sources in 3; command-driven repeated-mass behavioral regressions arrive with workflows in 4. |
| 12 | `181f77af` | remote | add selectable 4D-local UV before CFF | 3 | Selectable projected local 4D reconstruction and its core adapters. |
| 13 | `481a4301` | remote | expose deferred 3D representation selection | 2, 3, 4 | Shared representation selection belongs to 2; core adapter use to 3; CLI mode selection and exported behavior to 4. Superseded deferred machinery is not restored. |
| 14 | `2bbc6f59` | remote | add CFF direct-3D versus 4D-local scalar battery | 1, 2, 3, 4, 8 | Foundation lock/feature metadata 1; shared engine 2; core source/evaluator consumers 3; command-driven scalar battery, archived numerical fixtures and CLI 3D-representation/export workflows 4; historical rebase account 8. |
| 15 | `bfedaacc` | remote | document raised-energy CFF handoff | 8 | Historical handoff and rebase evidence. |
| 16 | `bd513be5` | remote | checkpoint local 4d uv diagnostics | 1, 2, 3, 4, 5, 8 | Spenso materialization 1; shared evaluator/generator 2; UV/source diagnostics 3; commands/profile/stability/schema and early command-driven acceptances 4; final signed acceptance/phase expectations 5; historical handoffs 8. |
| 17 | `5c3e7fca` | remote | fix(uv): rebuild exact CFF sources from edge provenance | 2, 3, 4 | Shared bound/source support 2; ownership-certified reconstruction and current architecture 3; generation/state/profile entrypoints 4. Required consumers/fixtures remain with their API owner. |
| 18 | `d54de204` | remote | checkpoint raised-energy CFF and typed local UV routes | 1, 2, 3, 4, 5, 8 | Tooling and dependency prerequisites 1; shared expression engine 2; typed UV routes, source certificates and evaluator migration 3; command-driven local-UV battery, workflow/profile/tracing 4; later signed acceptance refinement 5; historical records 8. Generated state_map.bin delta is explicitly excluded. |
| 19 | `c8e76317` | remote | refactor UV approximation lanes | 3 | Splits UV responsibilities into direct/local/projected owners; deleted intermediate diagnostic modules remain superseded. |
| 20 | `89a52f17` | remote | checkpoint: raised-energy UV/CFF handoff (not green) | 2, 3, 4, 8 | Shared evaluator/generation 2; source/UV correctness 3; profiling/state 4; durable handoff history 8. Temporary HANDOFF_REMOVE_FOR_PR files are net-absent in source tip. |
| 21 | `e30d799e` | remote | fix local UV counterterms and prepare validation candidate | 1, 2, 3, 4, 8 | Tensor and early Vakint compile prerequisites 1; shared engine fixes 2; local UV/source/mass correctness 3; tracing/workflows/resource execution 4; historical documentation 8. |
| 22 | `0ad34794` | remote | bound exact CFF cache capacity and correct macOS memory guarding | 1, 3, 4, 8 | Foundational parser/resource tooling 1; bounded core cache retention 3; CI/profile policy 4; historical notes 8. Current architecture accompanies each relevant functional owner. |
| 23 | `568679e1` | remote | satisfy new Rust 1.98 clippy lints | 1, 3, 4 | Clippy changes follow foundational utility consumers 1, graph/core 3 and evaluation/profile 4. |
| 24 | `1b760f63` | remote | apply CI resource overrides to the selected nextest profile | 4 | Resource overrides affect selected nextest execution profile; prerequisite floor/manual-exclusion metadata is separately staged 1. |
| 25 | `c1a839f5` | remote | retain numerical CFF oracles in isolated Nix tests | 2 | Shared numerical oracles and Nix eval feature metadata stay together. |
| 26 | `d3dd17d9` | remote | stream lower-sector CFF component products | 2, 3, 8 | Streamed shared component products 2; current CFF architecture 3; historical accounts 8. |
| 27 | `7b0b09ed` | remote | resolve JIT constants and shorten CFF intermediate lifetimes | 2, 3 | Shared generator lifetime shortening 2; core JIT constant resolution 3. |
| 28 | `e4f5f778` | remote | reuse exact CFF bases and compact completed products | 2 | Exact shared basis reuse and completed-product compaction. |
| 29 | `8a197b8c` | remote | preserve separate UV Taylor denominator topologies | 3, 8 | Owner/denominator-preserving Taylor composition and current guide 3; historical handoff 8. |
| 30 | `5a5a1bfc` | remote | stream independent outer CFF rank bounds | 3, 5, 8 | Independent outer ranks and complete UV source consumers 3; final signed acceptance 5; historical handoff 8. |
| 31 | `241cb46e` | remote | refresh verified GL20 snapshot and current architecture | 4, 5 | Saved event payload accompanies evaluation workflow 4 and final signed phase/event-label correction 5; current architecture is split by actual behavior. |
| 32 | `91142139` | remote | ALL GREEN (ALL_GREEN) | none (empty) | Empty remote marker: duplicated for provenance, contributes no content after squash. |
| 33 | `ac99d32d` | remote | ALL_GREEN: preserve independent CFF occurrences and simplify UV generation | 1, 2, 3, 4, 8 | Contraction foundations 1; independent shared occurrences 2; optimized source/UV generation and consumers 3; profiling/CLI 4; historical integration narrative 8. |
| 34 | `dd134b71` | remote | fix amplitude and LU phase conventions | 3, 4, 5, 8 | Physical phases, model conventions/sewing and signed oracles 5; shared core scaffolding 3 and workflow-only example syntax 4 are separated by hunk; historical phase record 8. |
| 35 | `5c8b46f3` | remote | fix factorized UV conversion and gauge-state sewing | 1, 3, 5, 6, 8 | Symbolica/Idenso and tooling 1; factor-preserving source/UV fixes 3; gauge-state sewing/models and signed fixtures 5; full Vakint API 6; historical record 8. |
| 36 | `26b57b78` | remote | fix factorized tensor bookkeeping and numerator lookup | 1, 3, 5, 6, 8 | Tensor/bookkeeping and dependency prerequisites 1; denominator-only incidence correction 3; numerator lookup and physical consumers 5; full Vakint input/coefficient API 6; historical record 8. |
| 37 | `32f27225` | remote | complete d-dimensional UV algebra in both Vakint input modes | 1, 4, 6, 8 | Spenso materialization and test-runner prerequisites 1; execution recipe/resource policy 4; complete d-dimensional API, settings, both owners and tests 6; historical phase evidence 8. Current Vakint architecture belongs to 6. |
| 38 | `3906e661` | review | review raised energy cff functionality tests and Rust patterns | 8 | Historical original review. |
| 39 | `6b9c51b9` | review | typeset raised energy CFF review | 8 | Historical review Typst and PDF. |
| 40 | `6f9e41de` | review | complete full-diff raised-energy CFF review | 8 | Historical full-diff coverage/validation plus Typst/PDF updates. |
| 41 | `9c8aa0f9` | review | implement raised energy CFF review fixes | 1, 2, 3, 4, 5, 8 | Tensor tests 1; shared validation/generator/evaluator tests 2; exact reconstruction, owner-specific integrated mass and outward core tests 3; CLI/persistence/evaluation tests 4; final signed acceptance 5; historical implementation reports 8. The 198 numeric snapshots are archived with the command-driven scalar matrix 4. |
| 42 | `ce6654cc` | review | store CFF coefficients as rationals | 2, 3, 8 | Rational storage across shared expression/surface/generation/evaluation 2 and core adapters/source consumers 3; historical implementation report 8. Current rational architecture is owned by 3. |

## Incoming review coverage

Incoming `a1140c90c334ff58a3b040ff3eae06d3645bf0bb` is reviewed against `78395e3ab3ddd8d8f62b2f674d7488484eace197`: **68 paths**, including 65 text and 3 binary artifacts. The retained ledger records 792 text hunks, 30321 added lines and 1718 deleted lines. The three serialized Symbolica blobs have byte/hash and source-provenance review; their contents were not semantically decoded. Bounded reviewer attribution is explicit; filename overlap alone does not assign functionality.

| Path | Reviewed owners | Reviewer(s) | Review correction |
| --- | --- | --- | --- |
| `.github/workflows/continuous-integration.yml` | C4 | /root | no |
| `.github/workflows/nix.yml` | C4 | /root | no |
| `CONTRIBUTING.md` | C7 | /root | no |
| `SPEED_UP_UV_CTS_FROM_4D.md` | C7 | /root/latest_cff_review | no |
| `bin/README.md` | C4 | /root | no |
| `bin/ram_watchdog.py` | C4 | /root | no |
| `crates/gammalooprs/src/cff/generation.rs` | C7 | /root/latest_cff_review | yes |
| `crates/gammalooprs/src/cff/mod.rs` | C7 | /root/latest_cff_review | yes |
| `crates/gammalooprs/src/graph/lmb.rs` | C7 | /root/latest_cff_review | no |
| `crates/gammalooprs/src/graph/three_d_source.rs` | C7 | /root/latest_cff_review | yes |
| `crates/gammalooprs/src/integrands/process/evaluators.rs` | C3 | /root | yes |
| `crates/gammalooprs/src/numerator/energy_degree.rs` | C7 | /root/latest_cff_review | yes |
| `crates/gammalooprs/src/utils/symbols.rs` | C3, C7 | /root/latest_cff_review | no |
| `crates/gammalooprs/src/uv/approx/direct_3d/branches.rs` | C7 | /root/latest_uv_review | no |
| `crates/gammalooprs/src/uv/approx/direct_3d/kernel.rs` | C7 | /root/latest_uv_review | no |
| `crates/gammalooprs/src/uv/approx/final_integrand.rs` | C7 | /root/latest_uv_review | no |
| `crates/gammalooprs/src/uv/approx/integrated.rs` | C6 | /root/latest_uv_review | no |
| `crates/gammalooprs/src/uv/approx/local_3d/residue_localizer.rs` | C7 | /root/latest_uv_review | no |
| `crates/gammalooprs/src/uv/approx/local_3d/tests.rs` | C7 | /root/latest_uv_review | no |
| `crates/gammalooprs/src/uv/approx/local_4d.rs` | C7 | /root/latest_uv_review | yes |
| `crates/gammalooprs/src/uv/approx/mod.rs` | C7 | /root/latest_uv_review | no |
| `crates/gammalooprs/src/uv/approx/projected_4d.rs` | C7 | /root/latest_uv_review | yes |
| `crates/gammalooprs/src/uv/forest.rs` | C7 | /root/latest_uv_review | no |
| `crates/gammalooprs/src/uv/hedge_poset.rs` | C7 | /root/latest_uv_review | no |
| `crates/gammalooprs/src/uv/mod.rs` | C7 | /root/latest_uv_review | no |
| `crates/gammalooprs/src/uv/orchestrator.rs` | C7 | /root/latest_uv_review | no |
| `crates/gammalooprs/src/uv/tests.rs` | C7 | /root/latest_uv_review | no |
| `crates/idenso/src/color/simplify.rs` | C1 | /root/latest_shared_review | no |
| `crates/idenso/src/color/test/mod.rs` | C1 | /root/latest_shared_review | yes |
| `crates/idenso/src/shorthands/chain.rs` | C1 | /root/latest_shared_review | no |
| `crates/idenso/src/shorthands/schoonschip/normalize_dots.rs` | C1 | /root/latest_shared_review | no |
| `crates/linnet/src/half_edge.rs` | C1 | /root/latest_shared_review | no |
| `crates/linnet/src/half_edge/hedgevec.rs` | C1 | /root/latest_shared_review | yes |
| `crates/linnet/src/half_edge/nodestore/test.rs` | C1 | /root/latest_shared_review | no |
| `crates/linnet/src/tree/child_pointer.rs` | C1 | /root/latest_shared_review | yes |
| `crates/spenso/src/network/contract.rs` | C1 | /root/latest_shared_review | no |
| `crates/spenso/src/network/graph.rs` | C1 | /root/latest_shared_review | yes |
| `crates/spenso/src/network/mod.rs` | C1 | /root/latest_shared_review | no |
| `crates/spenso/src/network/parsing/tensor_from_expression.rs` | C1 | /root/latest_shared_review | no |
| `crates/spenso/src/network/store.rs` | C1 | /root/latest_shared_review | yes |
| `crates/spenso/src/network/tests.rs` | C1 | /root/latest_shared_review | yes |
| `crates/spenso/src/shadowing/collect.rs` | C1 | /root/latest_shared_review | no |
| `crates/spenso/src/shadowing/tests.rs` | C1 | /root/latest_shared_review | yes |
| `docs/architecture/architecture-current.md` | C7 | /root | yes |
| `docs/architecture/exact-powered-denominator-cff-lifting.md` | C7 | /root | yes |
| `docs/architecture/local-4d-uv-benchmarks.json` | C7 | /root/latest_cff_review | no |
| `docs/architecture/local-4d-uv-correctness.json` | C7 | /root/latest_cff_review | no |
| `docs/architecture/local-4d-uv-dispatch-profile.json` | C7 | /root/latest_cff_review | no |
| `docs/architecture/local-4d-uv-performance.md` | C7 | /root | yes |
| `flake.nix` | C4 | /root | no |
| `tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/Cargo.lock` | C7 | /root/latest_cff_review | no |
| `tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/Cargo.toml` | C7 | /root/latest_cff_review | no |
| `tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/README.md` | C7 | /root | no |
| `tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/evidence/original_cli21_failure.json` | C7 | /root/latest_cff_review | no |
| `tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/evidence/original_cli21_failure.md` | C7 | /root/latest_cff_review | no |
| `tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/evidence/original_cli21_panic.txt` | C7 | /root/latest_cff_review | no |
| `tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/examples/import_remap.rs` | C7 | /root | no |
| `tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/import_remap_r1/receipt.json` | C7 | /root/latest_cff_review | no |
| `tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/import_remap_r1/tiny.symbolica` | C7 | /root/latest_cff_review | no |
| `tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/inputs/function_map_entries.symbolica` | C7 | /root/latest_cff_review | no |
| `tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/inputs/params.symbolica` | C7 | /root/latest_cff_review | no |
| `tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/inputs/provenance.json` | C7 | /root/latest_cff_review | no |
| `tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/src/main.rs` | C7 | /root | no |
| `tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/validation_summary.json` | C7 | /root/latest_cff_review | no |
| `tests/tests/test_runs/amplitude_phase_conventions.rs` | C5 | /root/latest_uv_review | no |
| `tests/tests/test_runs/scalar_3L_cross_section_inspects.rs` | C7 | /root/latest_uv_review | no |
| `tests/tests/test_runs/scalar_phase_conventions.rs` | C5 | /root/latest_uv_review | no |
| `tests/tests/uv.rs` | C7 | /root/latest_uv_review | yes |

## Deliberate review corrections

The current incoming-to-corrected delta contains **21 files**. Independent review covers the 15 production/test correction files, followed by complete public node-set checks, C1 formatting and public sum-fixture construction, plus current documentation sources and any separately verified rendered artifact. Initial formatting and compiler failures remain in the validation attempt ledger. Test oracles are preserved or strengthened; complete values replace incidental cache/class/order assertions.

| Path | Semantic owner | Review receipt |
| --- | --- | --- |
| `crates/gammalooprs/src/cff/generation.rs` | C7 | independent-corrections-review.json:coverage; /tmp/raised-stack-latest-review-20260914/c7-imports-fix-review.json |
| `crates/gammalooprs/src/cff/mod.rs` | C7 | independent-corrections-review.json:coverage |
| `crates/gammalooprs/src/graph/three_d_source.rs` | C7 | independent-corrections-review.json:coverage |
| `crates/gammalooprs/src/integrands/process/evaluators.rs` | C3 | independent-corrections-review.json:coverage; /tmp/raised-stack-latest-review-20260914/c3-clippy-fix-review.json |
| `crates/gammalooprs/src/numerator/energy_degree.rs` | C7 | independent-corrections-review.json:coverage |
| `crates/gammalooprs/src/uv/approx/local_4d.rs` | C7 | independent-corrections-review.json:coverage; /tmp/raised-stack-latest-review-20260914/c7-numeric-sign-fix-review.json |
| `crates/gammalooprs/src/uv/approx/projected_4d.rs` | C7 | independent-corrections-review.json:coverage |
| `crates/idenso/src/color/test/mod.rs` | C1 | independent-corrections-review.json:coverage |
| `crates/linnet/src/half_edge/hedgevec.rs` | C1 | independent-corrections-review.json:coverage |
| `crates/linnet/src/tree/child_pointer.rs` | C1 | independent-corrections-review.json:coverage; independent-corrections-review.json:followup; c1-format-fix.patch; fresh fmt receipts |
| `crates/spenso/src/network/graph.rs` | C1 | independent-corrections-review.json:coverage; independent-corrections-review.json:followup |
| `crates/spenso/src/network/store.rs` | C1 | independent-corrections-review.json:coverage; independent-corrections-review.json:followup |
| `crates/spenso/src/network/tests.rs` | C1 | independent-corrections-review.json:coverage; bulk-sum-public-builder-fix-review.json; /tmp/raised-stack-latest-review-20260914/c1-clippy-readability-fix-review.json |
| `crates/spenso/src/shadowing/tests.rs` | C1 | independent-corrections-review.json:coverage |
| `docs/architecture/architecture-current.md` | C7 | docs-edits-receipt.json |
| `docs/architecture/exact-powered-denominator-cff-lifting.md` | C7 | docs-edits-receipt.json |
| `docs/architecture/local-4d-uv-performance.md` | C7 | docs-edits-receipt.json |
| `docs/architecture/raised-energy-cff-motivation-audit.md` | C8 | docs-edits-receipt.json |
| `docs/architecture/raised-energy-cff-motivation-audit.pdf` | C8 | motivation-render/visual-inspection.json |
| `docs/architecture/raised-energy-cff-motivation-audit.typ` | C8 | docs-edits-receipt.json |
| `tests/tests/uv.rs` | C7 | independent-corrections-review.json:coverage |

C1 owns shared tensor/graph prerequisites; C2 the shared CFF engine; C3 adapters/evaluators; C4 workflows; C5 phases; C6 integrated algebra; C7 canonical UV optimization; C8 reports. Required signatures, consumers and tests remain with the compiling owner. Source-pinned benchmark artifacts retain their own scope; manual PySecDec and unresolved GL262 evaluator diagnostics do not become current acceptance through attribution.
