# RustRed acceptance coverage through three loops

This is a source-level inventory and corrected frozen-process passing-test
report. It covers every
existing end-to-end scalar comparison/reference input with one common mass and
at most three loops. Tensor-bearing native peers use the existing FeynKit
prepass, followed by RustRed scalar reduction and master substitution. The
shared harness gives both native stages an invalid FORM path; FORM is available
only to the explicitly separate AlphaLoop/MATAD oracle lanes.

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
Enabled means configured to run. After correcting the three native-pattern
conditions and the contracted-class fixtures without changing numerical
assertions, the fresh frozen matrix completed **76 selected invocations: 76
passed, 0 failed, 0 ignored**: 15 catalog/default, 17 three-loop peers, 43
legacy-complement, and one offline 38-terminal MATAD oracle. The offline test
uses FORM only to validate the terminal catalog; all scalar RustRed peers use
an invalid FORM path. Evidence is in
`target/vakint-short-pattern.OAgEGT/`; the invocation count complements,
rather than replaces, the 40-entry / 46-input inventory below.

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
They are enabled and pass in the corrected frozen matrix.

The five scalar matcher-class fixtures in RS and the five tensor-bearing
class inputs in `rustred_k6_pipeline_tests.rs` supplement this matrix; they do
not replace the original eleven three-loop acceptance bodies. A numerical-only
K6 terminal basis is permitted, but invalid-FORM-path checks remain mandatory.
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

## Native dependency and K6 routing preparation

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

### Genuine K6 asset and exact terminal manifest, runtime gate passed

The canonical RustRed producer generated byte-identical one/six-worker artifacts:
8,911,462 bytes, 623 rules, 5,639 cells, 38 typed corner terminals and 26 zero
masks. The existing CLI independently cold-loaded both files and reduced five
targets, including dotted and negative-power inputs. The exact 38-key catalog
join and disjoint 64-mask partition passed independent audit. The shipped SHA256
is `53bb589f98beaa735332cffbd080b174cfff5fc7a664e3499d5f667ad8fb8434`;
the actual algorithm ID is `rustred.source-port-original-domain.v1`. Evidence:
RustRed `target/spired-k6-producer.007TpB/RESULTS.md`. These are artifact/CLI
gates joined by the complete Vakint three-loop scalar parity matrix below.

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
existing artifact contract checks. The 17 three-loop peer tests and 43
legacy-complement tests also pass from fresh frozen processes with invalid FORM
paths for RustRed scalar evaluation.

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

The approved Spenso forwarding-wrapper migration now uses the pinned native
solve builder and returns its complete solution branches, retaining conditions
and free variables. Its existing matrix conversion supplies the linearity
gate; no local solver is added. Three focused tests cover unique exact values,
conditional/free branches and nonlinear rejection; all three now pass.
The existing K1/K3 asset contract now declares schema V5 and explicitly rejects
V4. Their bytes were regenerated by the current RustRed producers, not patched
or read through a compatibility shim; the unchanged exact
producer-versus-embedded-byte test now passes. The current combined library and
test-target check passes, but this is not a full GammaLoop workspace build.
A source-text census also found three obsolete calls in
GammaLoopRS; those are outside this Vakint/FeynKit test dependency slice and
remain untouched. All of these observations are distinct from successful
compilation or end-to-end RustRed acceptance.

Artifact decoding and recursive IBP application remain exclusively in the
existing RustRed library. Vakint contributes matcher/routing steering and its
existing terminal/master evaluation; these changes introduce neither another
artifact reader nor another reduction engine.
