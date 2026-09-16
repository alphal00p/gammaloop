# RustRed acceptance coverage through three loops

This is a source-level inventory and audited frozen-process test report.
The legacy input inventory has passing native peers; the two supplemental
all-class failures below are resolved by the signed MATAD routing correction.
It covers every
existing end-to-end scalar comparison/reference input with one common mass and
at most three loops. Tensor-bearing native peers use the existing FeynKit
prepass, followed by RustRed scalar reduction and master substitution. The
shared harness gives both native stages an invalid FORM path; FORM is available
only to the explicitly separate AlphaLoop/MATAD oracle lanes.

## Typed four-loop routing preparation (2026-09-16)

The existing K6 signed parent-slot and simultaneous momentum-transport checks
now live in one crate-private validator reusable by future parent-artifact
descriptors. K1/K3 behavior, defaults, dependencies and shipped assets are
unchanged. There is no new graph match, topology-name dispatch or IBP applier.

Native typed tests inspect all 19 registered four-loop classes, 123 surviving
physical slots and six nonidentity witnesses. They also reject missing/short
witnesses, altered coordinates, incorrect surviving momentum, changed
pinched-out parent slots and sign-flipped descriptors. Every four-loop class
still rejects RustRed admission: routing evidence alone is not a closed artifact.

The fresh scoped gate passes **39 tests**: 23 native library tests, with one
unchanged offline MATAD-catalog test ignored, plus all 16 K6 pipeline tests.
The native process and strict native pipeline settings use invalid FORM paths;
separate oracle lanes use FORM5. Formatting, `cargo check --locked -p vakint
--lib --tests` and scoped Clippy with `-D warnings` pass. The exact commands,
frozen hashes and outputs are retained at
`/tmp/vakint-four-loop-routing.dkHXAz/`.

This is not a rerun of the complete 83-test selection below and does not
establish four-loop scalar reduction or numerical parity. It replaces the
earlier diagnostic-only matcher inspection with durable typed regression tests.

## Symbolica 3 / current-artifact migration rerun (2026-09-16)

The same complete **83-test selection passes again: 83 passed, 0 failed,
0 ignored**, after the atomic dependency/artifact migration:

- RustRed `09cef8e3cf7487dded7803edc26df5d51bb9f500`;
- Symbolica, Graphica and Numerica
  `953e26e2754e9a4b918404fbdea590725d8e863d` (version 3.0.0);
- regenerated K6 bytes with source-port parent-plan tag `0x704`, the unrestricted
  root domain, and the unchanged set of 38 exact terminal keys.

The current shipped K6 artifact is 8,916,759 bytes, with 623 rules and 5,639
cells; its SHA256 is
`bcf45f876695c99516a7cefec2ae3cd3c5886a89d3da156391a6ea75237ff2ff`.

Every dependency is pinned to a published Git revision; no development-only
local path is required. K1 and K3 were also regenerated with the pinned
producer's existing CLI presets, cold-inspected and applied to canaries; their
bytes are identical to the previous shipped copies. K6's metadata label now
uses underscores, matching the generic RustRed CLI/Python input. Its private
expected fingerprint and embedded bytes change together; the mathematics and
terminal catalog do not change, and there is no obsolete-schema reader.

Spenso forwards Symbolica's native `SolutionSet`, including global coverage and
guards. Vakint accepts a routing witness only after complete, guard-free
coverage and a unique unconditional point are established. Native matching
depth setters replace the removed builder API without changing matching bounds.
Existing numerical assertions, tolerances, defaults and FORM oracle lanes are
unchanged. Native FeynKit/RustRed stages still use invalid FORM paths.

The table below remains the exact test inventory. New evidence is in
`/tmp/vakint-symbolica3-main.bMELju/`: published dependency resolution,
`cargo check --locked -p vakint --no-default-features --lib --tests`, frozen
binary hashes, the unchanged `selected.tsv`, and `acceptance/results.tsv` with
all nine passing batches. The debug correctness batches total 146.64 seconds
on one pinned core; this is **not** a release performance comparison. The
1m56s check and 6m24s test compilation are excluded from that runtime total.

The full actual workspace passes `cargo check --locked --workspace --lib` and
`cargo check --locked --workspace --tests --benches`. Wider migration adapters
use native Symbolica arithmetic/solvers and preserve zero-dimensional
continuous samples through a narrowly reviewed point-domain adapter. Symbolica
3 uses runtime signed-license/environment activation; the removed compiled OEM
scheme is not emulated (see the root contributor guide).

The fresh focused runtime gate also passes: **28 GammaLoop tests**, **4 Spenso
native-solver forwarding tests**, and **52 FeynKit tensor tests**. FeynKit's
pre-existing ignored explicit rank-ten, 893,025-pair cross-check remains
ignored. The GammaLoop selection covers six point-domain tests, runtime-license
initialization, exact high-precision conversion, complete/free/guarded momentum
routing (including cancellation of an input denominator), LMB mappings,
implicit-function derivatives and grid monitoring. All focused processes use
an invalid FORM path. Logs and frozen executable hashes are in the same
evidence directory under `focused-*`; the 12m25s compilation is not a runtime
benchmark. Scoped `cargo clippy --no-deps --locked -p vakint
--no-default-features --lib --test rustred_k6_pipeline_tests
--test integral_alphaloop_vs_matad_tests -- -D warnings` also passes; this is
not a whole-workspace lint claim.
These selected checks do not certify every GammaLoop runtime workflow, and do
not provide four-loop RustRed closure.

## Previous complete selected-matrix rerun (2026-09-16)

The complete recorded **83-test selection now passes: 83 passed, 0 failed,
0 ignored** on revision `433e42d3`. This rerun follows the signed MATAD routing
correction described below. It uses freshly built, hash-frozen binaries and
the exact union of the earlier native/catalog, three-loop-peer,
legacy-complement and offline-catalog manifests; no failing entry was dropped.

| Binary | Selected / passed |
| --- | ---: |
| Vakint library (native catalog/defaults and offline MATAD catalog) | 16 |
| Tensor reduction modes | 1 |
| Freeform evaluations | 2 |
| RustRed scalar evaluations | 18 |
| Analytic references | 19 |
| RustRed acceptance variants | 6 |
| Input matching | 5 |
| AlphaLoop/MATAD comparisons | 13 |
| Supplemental K6 all-class pipelines | 3 |
| **Total** | **83** |

Each binary is invoked once with all of its exact selected test names,
`--include-ignored --test-threads=1 --nocapture`; the runner checks the selected
name count, successful result count and executable hashes. Numerical assertions
and tolerances are unchanged. Native FeynKit/RustRed settings still contain an
invalid FORM path; only the explicitly separate oracle evaluations use FORM.
This includes the offline exact MATAD comparison for all 38 shipped terminals.
It is a complete rerun of the recorded selection, not a claim that every test
in every binary belongs to this matrix. The separate 34-test targeted gate
below includes additional routing diagnostics and overlaps this selection.

Build manifests and hashes are in `/tmp/vakint-orientation-full-gate.Ny652I/`;
the successful runner, all nine result rows and per-binary logs/timings are in
`/tmp/vakint-batched-83.baAvzc/`. An earlier one-process-per-test attempt hit
the Symbolica activation rate limit and was stopped; its partial logs are not
counted as mathematical failures or as this successful run. Batching after the
requested cooldown completed without activation errors.

Clippy now also passes with `-D warnings` for the Vakint library,
`rustred_k6_pipeline_tests` and `integral_alphaloop_vs_matad_tests`, using the
available Nix Clippy executable. These are correctness gates, not release
performance measurements. Four-loop RustRed artifacts and four-loop numerical
acceptance remain unfinished.

## Routing correction and targeted rerun (2026-09-16)

The contracted basketball basis contains reversed parent edges. MATAD's
defining identities require Vakint parent slots to map as
`[p4,p5,p6,-p1,p2,-p3]`. The old numerator adapter retained only the unsigned
permutation. It therefore gave the wrong sign to mixed products involving a
reversed loop coordinate, while squared momenta and scalar inputs were
unaffected. Independent scalar probes and two source-level mathematical
audits established this without treating a majority of reducers as proof.

The correction includes edge orientation in the existing numerator
Wick-rotation sign. Denominator numbering, one/two-loop maps, defaults, APIs,
RustRed artifacts, numerical tolerances and existing assertions are unchanged.

| Fresh selection | Passed | Failed | Ignored | Debug process seconds |
| --- | ---: | ---: | ---: | ---: |
| New `matad::routing_tests` | 3 | 0 | 0 | 0.49 |
| Entire `rustred_k6_pipeline_tests` | 16 | 0 | 0 | 56.13 |
| Entire `integral_alphaloop_vs_matad_tests` | 15 | 0 | 0 | 39.50 |

Both original all-class comparisons now reach all five three-loop classes.
The pipeline includes unit/nonunit mass, individual mixed scalar products,
and both terms of the original tensor numerator. The shared harness retains
an invalid FORM path for FeynKit and RustRed. Legacy oracle lanes alone use
FORM. An independent auditor reran the routing and full pipeline selections.
These are correctness timings, not release performance benchmarks.

Reproduce with `cargo test --locked -p vakint --no-default-features --lib
matad::routing_tests`, then the complete integration targets
`rustred_k6_pipeline_tests` and `integral_alphaloop_vs_matad_tests`; configure
the licensed Symbolica runtime and the existing FORM oracle executable.
That historical checkpoint used the coherent RustRed `ce92d3a7` / Symbolica
2.2 stack. Its check and build passed; Clippy was attempted but unavailable in
that development shell. The current Symbolica 3 pins and passing Clippy gate
are recorded at the top of this document. Raw logs, command evidence and
binary/input hashes for the older run are locally in
`/tmp/vakint-matad-orientation.yRdp6a/`.

That targeted checkpoint resolved every previously identified failure, but did
not itself rerun all 83 obligations. The subsequent complete selected-matrix
run above now supplies that missing evidence. Four-loop RustRed artifacts and
end-to-end four-loop acceptance remain unfinished.

## Historical inventory and audit

The five legacy files contain **40 test entries / 46 input executions**. These
are not 46 distinct integrals: power loops expand two entries, and several
backend-specific tests share workloads. Originally, 21 entries / 27 executions
had explicit enabled native peers and 11 three-loop entries had pending peers.
The remaining eight entries comprise four workload aliases and four genuinely
unrepresented settings/input variants. The latter now have isolated peers in
`rustred_acceptance_variants_tests.rs`; the original tests remain unchanged.

Before K6 activation, 25 legacy entries / 31 executions mapped to enabled peers
and 15 entries / 15 executions mapped to pending K6 peers. Genuine K6 bytes are
now shipped and the existing loader/routing adapter and all peer entrypoints
are enabled. All 40 entries / 46 executions have active native mappings.
Enabled means configured to run. The original report at revision `98550cc0`
counted **76 PASS records**, but its three-loop selection manifest contained
24 tests and its result table stopped after 17. The next test had a failing
stdout/stderr log without a corresponding result row. Those 76 records were
therefore not evidence that the complete selected matrix passed.

An independent rerun of the seven missing result entries used the same frozen
executables and FORM oracle, unchanged assertions, one CPU and one thread per
compute pool. Five passed and two failed. Joining those results with the
original completed rows gives **83 distinct selected tests: 81 passed,
2 failed, 0 ignored**:

| Group | Passed | Failed |
| --- | ---: | ---: |
| Native catalog/default | 15 | 0 |
| Three-loop peers and supplemental class checks | 22 | 2 |
| Legacy-complement checks | 43 | 0 |
| Offline 38-terminal MATAD catalog oracle | 1 | 0 |

The three previously unrecorded legacy-setting variants and both basketball
finite-part companions now have fresh passing evidence. The two failures are
`feynkit_three_loop_class_scalar_oracles_agree` and
`feynkit_rustred_three_loop_class_numerical_peers`. Both stop at
`I3L_pinch_1_6` on an AlphaLoop-versus-MATAD difference, before the harness
compares the RustRed numerical result or reaches the fifth class. At that
checkpoint, complete supplemental all-class parity was not established; the
fresh rerun above resolves these failures.

The original evidence is in `target/vakint-short-pattern.OAgEGT/`; the seven
independent reruns and their exact runner are in
`target/vakint-three-loop-audit.JA4fsa/`. Frozen executable and FORM hashes
match. The original source-input hash check differs only at the currently
edited workspace `Cargo.toml`; no binary was rebuilt. Native FeynKit and
RustRed lanes retain their invalid FORM paths. FORM remains available only
to separate comparison-oracle lanes and the offline catalog oracle. These
counts complement, rather than replace, the 40-entry / 46-input inventory.

## Complete mapping

File abbreviations below are relative to this directory:

- **AM**: `integral_alphaloop_vs_matad_tests.rs`
- **AR**: `integral_evaluation_analytic_tests.rs`
- **FF**: `integral_evaluation_freeform_tests.rs`
- **PC**: `integral_comparison_vs_pysecdec_tests.rs`
- **PR**: `integral_evaluation_pysecdec_tests.rs`
- **RS**: `rustred_scalar_evaluation_tests.rs`
- **RV**: `rustred_acceptance_variants_tests.rs`

“Same entry” already invokes the common multi-lane harness. `RS::pysecdec peers`
means `rustred_covers_applicable_pysecdec_peer_inputs_through_two_loops`, which
uses native oracles independently of PySecDec availability. Pending peers
were explicitly ignored, not silently skipped successful tests. Those K6
attributes are now removed; all numerical comparison bodies remain unchanged.

| Legacy entry | Executions | Native mapping | State |
|---|---:|---|---|
| AM::test_integrate_1l_no_numerator | 6 | Same entry; powers 1–6 | Enabled |
| AM::test_integrate_1l_no_numerator_squared_mass | 2 | Same entry; powers 1–2 | Enabled |
| AM::test_integrate_2l_no_numerator | 1 | Same entry | Enabled |
| AM::test_integrate_3l_basketball_a | 1 | AM::rustred_numerical_parity_3l_basketball_a | Passed |
| AM::test_integrate_3l_basketball_b | 1 | AM::rustred_numerical_parity_3l_basketball_b | Passed |
| AM::test_integrate_3l_no_numerator | 1 | AM::rustred_numerical_parity_3l_no_numerator | Passed |
| AM::test_integrate_3l_rank_4 | 1 | AM::rustred_numerical_parity_3l_rank_4 | Passed |
| AM::test_integrate_3l_rank_4_different_scales | 1 | AM::rustred_numerical_parity_3l_rank_4_different_scales | Passed |
| AR::test_integrate_1l_a | 1 | Same entry | Enabled |
| AR::test_integrate_1l_simple | 1 | Same entry | Enabled |
| AR::test_integrate_1l_simple_squared_mass | 1 | Same entry | Enabled |
| AR::test_integrate_1l_cross_product | 1 | Same entry | Enabled |
| AR::test_integrate_1l_cross_product_with_additional_symbols_numerator | 1 | Same entry | Enabled |
| AR::test_integrate_1l_dot_product_external | 1 | Same entry | Enabled |
| AR::test_integrate_2l | 1 | Same entry | Enabled |
| AR::test_integrate_3l | 1 | AR::rustred_numerical_parity_3l | Passed |
| AR::test_integrate_3l_rank_4 | 1 | AR::rustred_numerical_parity_3l_rank_4 | Passed |
| AR::test_integrate_3l_rank_4_additional_symbols_numerator | 1 | AR::rustred_numerical_parity_3l_rank_4_additional_symbols_numerator | Passed |
| AR::test_integrate_3l_rank_4_matad | 1 | AR::rustred_numerical_parity_3l_rank_4_matad | Passed |
| AR::test_integrate_3l_rank_4_matad_additional_symbols_numerator | 1 | AR::rustred_numerical_parity_3l_rank_4_matad_additional_symbols_numerator | Passed |
| AR::test_integrate_3l_matad | 1 | AR::rustred_numerical_parity_3l_matad | Passed |
| FF::test_integrate_1l_decorated_indices_alphaloop | 1 | Same entry | Enabled |
| FF::test_integrate_1l_decorated_indices_matad | 1 | Same entry | Enabled |
| FF::test_integrate_1l_decorated_indices_pysecdec | 1 | RV::rustred_decorated_one_loop_five_terms | Passed |
| PC::test_integrate_1l_pysecdec | 1 | RS::pysecdec peers | Enabled |
| PC::test_integrate_1l_pysecdec_non_unit_mass | 1 | RS::pysecdec peers | Enabled |
| PC::test_integrate_1l_pysecdec_non_unit_scale | 1 | RS::pysecdec peers | Enabled |
| PC::test_integrate_1l_pysecdec_num_rank_two | 1 | RS::pysecdec peers | Enabled |
| PC::test_integrate_1l_pysecdec_dot_product_external | 1 | RS::pysecdec peers | Enabled |
| PC::test_integrate_2l_pysecdec | 1 | RS::pysecdec peers | Enabled |
| PC::test_integrate_2l_pysecdec_pinched | 1 | RS::pysecdec peers | Enabled |
| PC::test_integrate_2l_pysecdec_pinched_other_lmb | 1 | RS::pysecdec peers | Enabled |
| PC::test_integrate_2l_pysecdec_rank_four_num | 1 | RS::pysecdec peers | Enabled |
| PC::test_integrate_3l_pysecdec | 1 | RV::rustred_three_loop_scalar_nonunit_scale | Passed |
| PC::test_integrate_3l_rank_4 | 1 | RV::rustred_three_loop_rank_four_large_external_vectors | Passed |
| PC::test_integrate_3l_rank_4_matad | 1 | RV::rustred_three_loop_rank_four_five_terms | Passed |
| PR::test_integrate_1l_simple | 1 | Alias of AR::test_integrate_1l_simple | Enabled alias |
| PR::test_integrate_1l_cross_product | 1 | Alias of AR::test_integrate_1l_cross_product | Enabled alias |
| PR::test_integrate_1l_cross_product_with_additional_symbols_numerator | 1 | Alias of AR::test_integrate_1l_cross_product_with_additional_symbols_numerator | Enabled alias |
| PR::test_integrate_3l_o_eps | 1 | Alias of AR::test_integrate_3l_matad and its peer | Passed alias |

The aliases preserve the same integral, mass/scale substitutions, external
vectors, normalization and epsilon depth. Their original optional numerical
references use lower precision; native peers use the corresponding existing
32-digit analytic references, not the rounded PySecDec numbers. Alias entries
are counted as legacy obligations, not additional independently executed tests.

## Four previously missing variants

The new RV peers reuse existing comparison helpers; they introduce no backend,
tensor algorithm, reusable helper API or fallback. The one-loop case retains
decorated indices/symbol attributes, the custom renormalization-scale symbol,
MSbar normalization, power two and five epsilon terms. The three-loop cases
retain the optional comparisons' nonunit renormalization scale, large versus
small external vectors, and four versus five epsilon terms. Native parity uses
32 digits and relative tolerance `1e-25`; loose Monte Carlo tolerances are not
transferred to RustRed acceptance. Five-term cases use the existing analytic/
MATAD lanes, which can select an eligible method rather than forcing AlphaLoop
beyond its supported epsilon depth. PySecDec itself remains supplemental.

## Finite parts, topology coverage and exclusions

Both original basketball comparisons request three epsilon terms at three
loops: they stop at epsilon^-1. They are preserved verbatim. Separate
RV companions `rustred_basketball_a_finite_part` and
`rustred_basketball_b_finite_part` request four terms, including epsilon^0.
These two added executions are not included in the 40-entry legacy census.
They pass in the independent seven-test frozen-binary rerun recorded above.

The five scalar matcher-class fixtures in RS and the five tensor-bearing
class inputs in `rustred_k6_pipeline_tests.rs` supplement this matrix; they do
not replace the original eleven three-loop acceptance bodies. A numerical-only
K6 terminal basis is permitted, but invalid-FORM-path checks remain mandatory.
In the initial pre-correction audit, the scalar matcher-class peers passed and
the tensor prepasses matched FORM for all five classes, but the two subsequent
all-class scalar comparisons failed at the fourth class before reaching the
fifth. At that fourth class,
the epsilon^-3 real coefficient is approximately `0.30116343610153126`
for AlphaLoop and `0.30102068072626876` for MATAD, a relative difference
of about `4.741e-4` against the unchanged `1e-25` tolerance. The original logs
alone did not identify the responsible reducer. Subsequent signed-routing
diagnostics isolated and corrected the MATAD numerator orientation, and the
complete 83-test reruns above now pass without weakening these assertions.
The old explicit unsupported-K6 inventory assertion now requires a successful
nonzero parent reduction, and its existing numerical peer inventory is live.
Before acceptance, execute every peer: merely removing ignore attributes or
passing a source-level inventory is not evidence of numerical success.

Excluded from the 40 entries: four-loop cases, master-evaluator-only checks,
matching-only/tensor-only tests, and offline MATAD catalog utilities.
`FF::test_integrate_1l_decorated_indices_fmft` is actually a four-tadpole,
four-loop input despite its name. `PR::test_integrate_2l_different_masses` has
three independent symbolic masses; equal numerical substitutions do not turn
that input into a common-mass family without changing it. Defaults and all
historical FORM-backed methods remain unchanged.

## Historical native dependency and K6 routing preparation

The following preparation checkpoints predate the completed migration and
acceptance gates at the top of this document. Their recorded intermediate
failures and older artifact identities are retained as historical evidence,
not as current build or acceptance status.

The approved pinned-API migration now passes a combined Vakint/Spenso library
and test-target check. Fresh frozen libtest processes passed three Spenso
forwarding tests, five native K6 candidate-data tests, six existing RustRed
artifact/catalog tests, and four backend/default tests. Native processes used
an invalid FORM path. These are focused unit gates, not the end-to-end matrix
above. Evidence is in `target/vakint-approved-routing.C5fu9J`.

The six routing tests did execute but did not reach their mathematical
assertions: the first calls topology construction before the existing Vakint
symbol initializer, defining `uedge` without its symmetric attribute and
poisoning the shared symbol initializer for the other five. A separate normal
full-lib run of the unchanged frozen binary had
22 passes, one failure and one ignored test: ambient initialization allowed
five routing tests to pass, but the five-class test exposed a missing native
`allow_new_wildcards_on_rhs(true)` setting while constructing component
templates. The user explicitly approved both setup corrections: existing symbol
initialization first in all six tests and the native flag only on intermediate
template construction. Those seven lines are now applied, with all assertions
unchanged, and all six now pass in isolated fresh processes in the activation
gate. The historical failed logs remain preserved; ambient initialization was
not used as a substitute for that isolated rerun.
The separate all-38 offline MATAD oracle also passes at the configured
four-term Laurent order, with actual
FORM 5.0.1 used only on the oracle side and an invalid FORM path on candidate
master finalization. This checks candidate data, not K6 scalar integration.

The pinned Symbolica revision no longer exposes the former atom-level
`solve_linear_system`. The required Vakint basis-forcing call now validates
linearity natively and accepts only one unconditional, fully determined point
solution, in requested-variable order. FeynKit's two existing projector solvers
instead use the native rational-polynomial matrix interface, which is the same
field solve used previously; symbolic dimension poles remain exact rational
functions. No tensor algorithm or public signature is changed. Seven focused
unit tests cover negative orientation, nonlinear/dependent/conditional routing
rejection, generic dimension poles and singular/inconsistent field solves.
The routing subset now passes after the approved setup repairs; the three
FeynKit projector tests also pass in separate fresh processes.

An additional topology-module regression derives all five three-loop classes
from the existing matcher census, rather than another hand-maintained topology
list. It checks all six parent slots, exact unit-Jacobian basis transport,
simultaneous scalar-numerator substitution, nonunit shared mass, and unchanged
shifted/dotted/negative powers. Its component assertions now pass after the
approved native template-construction correction. Contracted classes call
`force_an_lmb`, so their numerator coordinates cannot simply be assumed to be
the parent artifact coordinates even when propagator IDs agree.

The same test now also exercises simultaneous component-wise transport and
conversion back to atomic scalar products, against an independent fixed-index
component construction. This matters because substituting a sum directly into
`dot(k(i), k(j))` does not perform bilinear expansion: the existing scalar
lowerer expects atomic vector arguments. External-vector and shared-mass
spectators are retained. No second graph match is introduced. Production
retention of the basis witness is now implemented in private `Integral`
metadata. The contraction constructor retains the defining parent expression
and the selected LMB momenta `k_new(i) = q_old(edge_i)` before the existing
solve rewrites the propagators in the opposite direction. The public
`force_an_lmb` signature is unchanged. The all-five-class test compares these
retained coordinates against its independent solved witness, then uses the
retained coordinates for the simultaneous scalar/component checks. A separate
test checks negative orientation, identity routing and constructor retention.
The historical full-lib run passed five routing tests but exposed the missing
component-template wildcard permission; the earlier isolated run also exposed
missing symbol initialization. Both approved setup fixes preserve the original
assertions. All six routing tests now pass in separate fresh processes, including
the all-class component check. Production K6 scalar steering is wired through
the retained witness; its numerical peers pass in the corrected frozen matrix.
The short-form matching fixtures now use only actually bound power wildcards,
and their adapter assertions pass.

The existing generic loader, typed terminal catalog, memoized reducer and
homogeneity materialization now consume a genuine V5 K6 artifact and its exact
terminal manifest. The numerator uses the matcher's retained simultaneous
parent routing before the existing RustRed scalar service. No new loader,
recursive applier or graph match is introduced.
Exact offline MATAD-basis terminal projections are preferred over unnecessary
20,000-digit literal tables.

### Earlier K6 asset and exact terminal manifest, catalog runtime gate passed

The earlier canonical RustRed producer generated byte-identical one/six-worker artifacts:
8,911,462 bytes, 623 rules, 5,639 cells, 38 typed corner terminals and 26 zero
masks. The existing CLI independently cold-loaded both files and reduced five
targets, including dotted and negative-power inputs. The exact 38-key catalog
join and disjoint 64-mask partition passed independent audit. Its then-shipped SHA256
was `53bb589f98beaa735332cffbd080b174cfff5fc7a664e3499d5f667ad8fb8434`;
the actual algorithm ID is `rustred.source-port-original-domain.v1`. Evidence:
RustRed `target/spired-k6-producer.007TpB/RESULTS.md`. These are artifact/CLI
gates were joined by passing legacy-input native peers and, at that checkpoint,
two failing supplemental class comparisons. Both failures are now resolved;
the current regenerated artifact hash and complete passing matrix are above.

`src/rustred_evaluation/terminal/k6.rs` now records the actual 38 canonical
corner keys checked by RustRed's 623-rule unit-mass program. Six authenticated
symmetry orbits need five distinct expressions: two orbits give the existing
one-loop tadpole cubed; the remaining product orbit gives tadpole times sunset;
the three irreducible anchors use the existing four-, five- and six-line MATAD
records. The full five-line expression retains its outer minus sign and its
`miD5`, `miBN`, `miT111` and `Gam` terms. No new numerical table is introduced.

These are unit-mass Minkowski values in Vakint's internal `ep` convention.
They contain neither the terminal's later classical mass power nor the
per-loop `(m²)^(-3*ep)` normalization. The keys already use parent slots, so the
earlier source-reference slot exchange must not be applied again.

The module was deliberately `cfg(test)` while genuine K6 bytes were unavailable.
The same unchanged records are now consumed by the existing authenticated
artifact/catalog loader, bound to the observed schema/algorithm/fingerprint and
complete actual key set. Four tests check the complete 38-key set and six
symmetry orbits, native exact comparison to the three recorded MATAD anchors,
the existing tadpole/sunset product signs, and absence of mass variables or
numerical-tail symbols. All four pass in both the historical and current frozen
GammaLoop libtest binaries. No fake `ClosedArtifact` or alternate loader is added;
`TerminalCatalog::compile` authenticates the shipped artifact and requires exact
equality with its complete typed terminal set before use. The current 15-test
native catalog/default gate passes, including master finalization and the
existing artifact contract checks. The audited union also contains 22 passing
three-loop selections and 43 passing legacy-complement checks, with invalid
FORM paths for RustRed scalar evaluation. Two supplemental all-class
tensor/scalar selections fail as described above; they must not be hidden by
counting only existing PASS rows.

Two additional checks are prepared. The enabled native-only test now passes: it
materializes the five different terminal expressions through the finite
three-loop term at `m²=13/10`, with an invalid FORM path and the existing MATAD
master finalizer. The separate ignored offline test now also passes when
explicitly selected: it compares all 38 candidate values, including their
symbolic mass restoration, against actual MATAD reductions at the configured
four-term Laurent order. The exact symbolic differences vanish after that
finalization; this is not an all-dimension rational-function identity check.
The latter requires `VAKINT_K6_ORACLE_FORM_PATH` or
`FORM_PATH`; it is a catalog-data oracle, not RustRed runtime acceptance or
proof of artifact closure. Neither test constructs a replacement catalog or
another master-evaluation implementation.

At that preparation checkpoint, the Spenso forwarding wrapper used the pinned
native solve builder and returned complete solution branches, retaining
conditions and free variables. Its existing matrix conversion supplied the
linearity gate; no local solver was added. Three focused tests covered unique
exact values, conditional/free branches and nonlinear rejection. The K1/K3
asset contract declared schema V5 and rejected V4; their bytes were regenerated,
not patched or read through a compatibility shim, and the unchanged exact
producer-versus-embedded-byte test passed. The then-current combined library
and test-target check was not a full GammaLoop workspace build. A census found
three obsolete GammaLoopRS calls outside that dependency slice, which were
left untouched at that time. The current migration adapts those calls, adds
the cancellation-guard regression, and passes the whole-workspace compilation
and focused runtime gates recorded above.

Artifact decoding and recursive IBP application remain exclusively in the
existing RustRed library. Vakint contributes matcher/routing steering and its
existing terminal/master evaluation; these changes introduce neither another
artifact reader nor another reduction engine.
