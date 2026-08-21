= API documentation debt register
<api-documentation-debt-register>

#quote(block: true)[
  *Status:* Open register

  *Review date:* 2026-08-21

  *Reviewed baseline:* `4e77d8804aa3394992572575dab5d3f9dedf65a9`

  *Scope:* Native Rustdoc and generated Python reference surfaces for GammaLoop,
  Linnet, Spenso, Idenso, and Vakint, including adjacent public Rust crates whose
  support status still needs classification. CLI and settings entries are included
  only where they contradict an API contract.
]

== Purpose

This is the canonical backlog for known API-reference deficiencies. The broader
#link("documentation-improvement-plan.typ")[documentation improvement plan]
sets the audience and product-content strategy; this register records the concrete
reference surfaces that still lack accurate contracts, useful examples, or an
explicit support boundary.

The companion #link("manual-content-debt-register.typ")[manual content debt
register] owns findings about audience, scientific explanation, complete workflows,
terminology, citations, and result interpretation. API items link those workflows
without duplicating their prose requirements here.

A generated signature is inventory, not complete documentation. A supported API
item is complete only when a reader can determine:

- its domain purpose and the conditions under which it is appropriate;
- parameter meaning, conventions, units, defaults, and related preconditions;
- return-value structure, ownership or mutation behavior, and interpretation;
- errors, panics, side effects, feature or license requirements, and stability;
- at least one minimal, independently runnable example with a checked result; and
- the conceptual guide to read when correct use requires more than local prose.

Rust references should additionally state important invariants, complexity or
allocation behavior, feature gates, and safety or panic conditions. Python examples
must include imports and setup, must not depend on unnamed variables or private
state, and must distinguish illustrative fragments from verified journeys.

== Maintenance contract

Identifiers in this register are stable. Use `Open`, `In progress`, `Blocked`, or
`Resolved` as the status. When API work exposes a new deficiency, add or refine an
entry in the same change. Mark an entry resolved only when its completion condition
is backed by the source documentation and the relevant example-verification tier.
Keep resolved entries through the next documentation review with the resolving
change ID, then remove them from the active register.

The counts below are an audit snapshot, not quality scores or permanent targets.
They measure whether generated module-index rows have summaries and whether Python
members have nonempty prose; they do not measure scientific correctness. Future
automation may regenerate these counts, but the semantic and example-quality gaps
remain a maintained review responsibility.

== Audit snapshot

=== Native Rustdoc

#table(
  columns: (1.2fr, 0.8fr, 2.4fr),
  table.header([*Component*], [*Described index rows*], [*Immediate interpretation*]),
  [`gammalooprs`], [66 / 884 (7%)], [No crate orientation and very sparse contracts across the scientific core.],
  [`gammaloop-api`], [11 / 160 (6%)], [The supported state/session facade has no end-to-end Rust workflow.],
  [Linnet], [94 / 216 (43%)], [Best baseline, but several prominent claims are stale or incorrect.],
  [Spenso], [113 / 480 (23%)], [Central structure, contraction, and network contracts remain fragmented.],
  [`spenso-hep-lib`], [4 / 25 (16%)], [Scientifically sensitive conventions and identities are mostly absent.],
  [`spenso-macros`], [0 / 1], [Macro availability exists without an orienting contract or example.],
  [Idenso], [43 / 109 (39%)], [Item prose exists unevenly, while crate/module and rewrite-stage orientation is missing.],
  [Vakint], [0 / 84], [The public Rust engine and result types are effectively an undocumented inventory.],
)

The curated Rust support catalog attaches an example to each supported item, but
many examples only import a type. Those examples prove name resolution, not behavior
or semantics, and do not satisfy the completion standard above.

=== Python

#table(
  columns: (1.2fr, 1fr, 2.2fr),
  table.header([*Package*], [*Surface snapshot*], [*Immediate interpretation*]),
  [`gammaloop-python`], [40 exports; 262 members], [65 member descriptions are blank, chiefly around public integration-result records; core examples depend on an unspecified state and point.],
  [`linnet-py`], [21 exports; 178 members], [54 member descriptions are blank and the generated reference contains no runnable examples.],
  [`spynso3`], [14 exports; 119 members], [34 member descriptions are blank; numerous examples are fragments or violate the documented registration workflow.],
  [`idenso-community`], [16 functions], [All exports have prose, but only nine have examples and none is runtime-verified.],
  [`vakint-community`], [4 classes; 20 members], [Member prose exists, but examples are non-standalone and exception contracts are absent.],
)

The Python stub checks currently verify exports and signatures. Docstring examples are
embedded prose rather than structured examples and are not executed, so example
presence must not be reported as example validity.

== Integrity defects

=== APIDOC-001 · GammaLoop sample coordinates and precision

*Priority:* P0 · *Status:* Resolved · *Resolved in:* `umvyopkm` · *Surface:* GammaLoop CLI, Rust, and Python

`inspect` and `approach` describe momentum input as `(p0 px ...)`, while
#link("../../crates/gammaloop-api/src/commands/evaluate_samples.rs")[the shared input builder]
requires flattened `(px, py, pz)` triplets. The same CLI options describe
`use_arb_prec` as f128 while Python calls it arbitrary precision.

The resolved contract now distinguishes integration coordinates from flattened spatial loop
momenta and rejects incomplete triplets consistently. It also records the legacy precision
behavior without changing it: despite its name, `use_arb_prec` selects configured f128 and
falls back to Arb only when no f128 level exists; ordinary CLI, Python, and Rust results remain
f64. Source-level regressions cover the help text, malformed shapes, and selected stability level.

*Complete when:* every interface states the same coordinate layouts and precision
levels, rejects malformed input consistently, and links to one canonical evaluation
contract with a checked example.

=== APIDOC-002 · Incorrect Linnet Rustdoc claims

*Priority:* P0 · *Status:* Resolved · *Resolved in:* `umvyopkm`, completed by `rtmnuuuy` ·
*Surface:* Linnet Rustdoc

#link("../../crates/linnet/src/parser/mod.rs")[The parser module] names a nonexistent
`HedgeGraphSet`, uses the nonexistent `linnet::dot_parser` path, and hides its example
behind `ignore`. The half-edge module calls the structure a DCEL and claims face
traversal although no embedding or rotation-system API supplies that operation.

The corrected Rustdoc now describes the implemented incidence/involution model, explicitly
excludes embeddings and face traversal, names `DotGraphSet` and `linnet::parser`, and runs a
parse → inspect → serialize → reparse example. Linnet check, clippy, and all doctests pass.
The full all-product Rustdoc build subsequently exposed the exported `dot!` macro link as
module-relative; follow-up `rtmnuuuy` qualifies it through the crate root, where Rustdoc can
resolve it under warnings-as-errors.

*Complete when:* the public description matches the implemented half-edge model, all
named paths and types resolve, and the conceptual parser example runs.

=== APIDOC-003 · Contradictory Spynso network examples

*Priority:* P0 · *Status:* Resolved · *Resolved in:* `umvyopkm` · *Surface:* `spynso3` Python reference

The generated `TensorNetwork` example executes named tensor expressions without
registering their data or supplying a library. The maintained Python guide correctly
requires registration before execution.

The source docstring, generated reference, and authored guide now share the same complete
Euclidean identity-matrix journey: create the structure and library, register dense data,
construct and execute the network, retrieve the tensor, and assert its four components.

*Complete when:* generated and authored examples share one independently runnable
registration-to-execution journey with an asserted result.

=== APIDOC-004 · Idenso identity claims require scientific review

*Priority:* P0 · *Status:* In progress · *Surface:* `idenso-community` Python reference

The `simplify_color` prose contains index-inconsistent contractions and a Fierz formula
that is not consistent with the stated normalization. Gamma-trace and gamma-five
claims also lack the dimensional-scheme qualifications needed to interpret them.

*Progress (`umvyopkm`):* the public formulas now use consistent free and dummy indices, keep
$T_R$, $C_F$, and $C_A$ distinct, state the four-dimensional gamma-five boundary, and explain
that unsupported open structures form a partial normal form. The public convention page is
separated from its developer source map. Exact Python-runtime verification of every advertised
example remains open. Eight exact Rust reference cases covering the advertised ordinary/gamma-five
traces, generator normalization, Fierz, Casimir, and structure-constant contraction pass when
run as separate processes under the restricted Symbolica license.

*Complete when:* the identities are checked against the implemented conventions,
free and dummy indices are consistent, dimensional assumptions are explicit, and
each advertised identity has a runtime-verified example.

=== APIDOC-005 · Vakint reference examples and public terminology

*Priority:* P0 · *Status:* In progress · *Surface:* Vakint Rust and Python references

Python snippets assume undefined expressions, masses, momenta, engines, or Symbolica
objects. Public prose also contains misspellings such as `susbstitute_masters`,
`direct_numerical_substition`, “espilon,” and “vaking integral.” The crate example
named `twoloop_matching` configures an evaluation backend and probes FORM despite the
matching-only tutorial's stated dependency boundary.

*Progress (`umvyopkm`):* matching/canonicalization now use an empty evaluation order, source
docstrings contain complete imports and setup, evaluation examples name their FORM or external
backend requirements, terminology is corrected, and legacy misspelled keyword parameters are
identified as compatibility constraints. Runtime checks for the matching and Laurent-series
journeys and a complete exception contract remain open.

*Complete when:* matching and evaluation examples are separately named, standalone,
and checked; public terminology is corrected or explicitly preserved as a compatibility
constraint; and each example declares its external-tool requirements.

=== APIDOC-006 · Spenso parallel state and error-label integrity

*Priority:* P1 · *Status:* Open · *Surface:* Spenso Rust API

The Rustdoc review found two API-integrity boundaries that prose alone cannot close.
`Network::execute` refreshes the cached `Network::state`, whereas `execute_parallel` rewrites
the graph and store without refreshing that field. In addition, the display messages attached
to `InvalidDotFunction` and `NonSelfDualTensorPower` are exchanged. The current Rustdoc states
the parallel-state caveat so callers are not misled, but it remains surprising runtime behavior.

*Complete when:* sequential and parallel execution expose one tested state contract, the two
error variants render their own conditions, and the result/state/error Rustdoc and regressions
agree.

=== APIDOC-007 · GammaLoop accepted-but-unsupported calculation options

*Priority:* P0 · *Status:* Open · *Surface:* GammaLoop generation and runtime settings

Several public inputs currently cross the parsing or deserialization boundary without a safe
implemented meaning. Powered coupling constraints such as `QCD^2==n` retain the coupling name
and order but discard the parsed power. `integrated_phase = both` is accepted by the settings
type although both integration branches are unimplemented. Complex particle masses reach an
explicit panic, and the default flux path implements only one- and two-particle incoming states.
These are capability limits as well as reference hazards: documenting the general syntax without
the runtime boundary would invite silently different selections or process termination.

*Complete when:* each unsupported value is rejected with a local actionable error or gains a
tested implementation; the generated settings/CLI reference states the supported subset; and
maintained examples cover both one accepted and one rejected value for each boundary.

=== APIDOC-008 · Spenso example isolation under restricted Symbolica

*Priority:* P1 · *Status:* Resolved · *Resolved in:* `rtmnuuuy` · *Surface:* Spenso Rust
examples and verification harness

The new `Network` contraction example compiles and produces `[104, 236]` when run in a fresh
process. In the combined documentation-example binary, however, the existing direct-contraction
tutorial starts the one permitted restricted Symbolica instance and the later `Network`
construction aborts while attempting to start another. Serial test execution does not change
the outcome. This prevents the aggregate harness from treating both examples as runnable in one
process and may also expose an initialization boundary relevant to unlicensed users of the
`shadowing` feature.

*Resolution and build evidence (2026-08-21):* the generated harness now runs every executable
manual example and the Spenso `Network` catalog examples in fresh child processes. Its shared
libtest path serializes those child launches; every runner lets process-global Symbolica state end
with each isolated example. With `SYMBOLICA_LICENSE` explicitly absent, both Cargo and Nextest pass
all 87 examples, including the direct-contraction tutorial and both network cases. The licensed
suite passes the same assertions. Pages therefore remains a pure, cacheable derivation and the
credential is not part of its derivation metadata. Restricted single-instance initialization
remains a runtime constraint, not a documentation-harness defect.

*Complete when:* maintainers determine whether the two API paths should share one idempotent
Symbolica initialization; otherwise the verification harness isolates the examples in separate
processes or explicitly records the network example as compile-only in the aggregate tier. The
Rustdoc example must remain independently compiled and executed with its asserted result.

== Native Rustdoc gaps

=== APIDOC-101 · GammaLoop scientific-core contract

*Priority:* P1 · *Status:* Open · *Surface:* `gammalooprs`

Add crate-level scope and non-goals, a supported-versus-internal boundary, and module
orientation. Document the sampling, weighting, precision, merging, and observable
invariants around `HasIntegrand` and `evaluate_samples_raw`. Explain global/runtime
settings lifecycle, precedence, units, defaults, errors, and cancellation-context
semantics.

*Complete when:* a Rust user can evaluate and interpret one maintained integrand
without learning the execution contract from implementation code.

=== APIDOC-102 · GammaLoop Rust facade workflow

*Priority:* P1 · *Status:* Open · *Surface:* `gammaloop-api`

At the reviewed baseline, the crate had no crate/module orientation or Rustdoc
`Examples` sections. Central state-load options, loaded-state types, and `cli_session`
needed contracts for filesystem effects, read-only behavior, ownership, errors, and
persistence.

*Source-audit findings (2026-08-21):* loading initializes process-global GammaLoop/Symbolica and
tracing state and applies boot commands before returning; it is not pure deserialization. Dropping
a loaded state does not save it. Read-only mode prevents GammaLoop-managed writes to the active
state tree but permits explicit exports elsewhere and cannot constrain external processes such as
the `!` shell command. `model_file` affects only an existing manifested state, and an unmanifested
directory starts a blank state rather than restoring its contents. A borrowed CLI session prevents
direct access to the loaded bundle until it is dropped. Most commands expose effects through
state, settings, logs, or files and return `CommandOutput::None`; only evaluation and integration
currently return a structured Rust value.

*Progress (`rtmnuuuy`):* crate and command-session orientation now document the load → borrowed
session → command → typed result lifecycle, ownership, persistence, filesystem effects, read-only
scope, process-global initialization, history recording, control flow, and error boundaries. A
compile-checked `no_run` example exercises a blank read-only state and typed command result without
pretending that licensed global initialization is safe inside the aggregate doctest process. A
maintained runnable evaluation/integration result and explicit settings-precedence table remain
open.

*Complete when:* one runnable workflow demonstrates load or create → command or
generation → evaluation or integration → interpreted output, with license and
side-effect boundaries stated.

=== APIDOC-103 · Vakint native engine and result model

*Priority:* P1 · *Status:* Open · *Surface:* `vakint`

The crate root, exported modules, momentum representations, error type, integral,
backend options, evaluation order, settings, engine, expression, and numerical result
types lack native Rustdoc contracts.

*Complete when:* the reference defines mathematical notation and normalization,
epsilon-order semantics, backend selection and dependencies, stage transitions,
temporary-file/error behavior, and a checked Laurent-series calculation.

=== APIDOC-104 · Spenso HEP conventions

*Priority:* P1 · *Status:* Open · *Surface:* `spenso-hep-lib`

Dirac and Weyl gamma matrices, gamma five, projectors, sigma tensors, and library
constructors lack basis, metric, index-order, normalization, and required-structure
contracts. An implementation comment notes that a required structure shape is not
checked.

*Complete when:* the public preconditions and failure behavior are explicit and
checked Clifford, projector, and SU(3) identities demonstrate the chosen conventions.

=== APIDOC-105 · Spenso structure and network invariants

*Priority:* P1 · *Status:* In progress · *Surface:* `spenso`

`TensorStructure`, direct contraction, `Network`, execution modes/results/errors, and
strategy types lack one cohesive contract. Coordinate order, free/dummy matching,
dense/sparse guarantees, planning versus execution ownership, feature gates, and
complexity are distributed or absent.

*Progress (`umvyopkm`):* central Rustdoc now defines free/dummy matching, dual compatibility,
row-major coordinates, dense/sparse result behavior, strategy ownership, the in-place/no-stored-
plan lifecycle, result ownership and kinds, bounded execution, and the error taxonomy. A Rustdoc
three-tensor example contracts $A_(i j) B_(j k) c_k$ to the surviving $i$ slot and asserts
components `[104, 236]`; it passes in both licensed and restricted aggregate runs, with the
process-isolation boundary recorded by APIDOC-008. Sparse/dummy/dual examples, lower-level
graph/store/library contracts, feature orientation, and broader structure-method coverage remain.

*Complete when:* those contracts are local to the relevant items and one checked
three-tensor example demonstrates contraction order and result structure.

=== APIDOC-106 · Idenso rewrite-stage contracts

*Priority:* P1 · *Status:* Open · *Surface:* `idenso`

Add crate and exported-module orientation, then document `IndexTooling`, cooking,
selective expansion, representation setup, dummy-index freshness, termination or
idempotence expectations, and errors. Distinguish currently implemented FORM-derived
rules from target specifications.

*Complete when:* one checked Dirac/color workflow makes each stage, invariant, and
implemented rule boundary visible from Rustdoc.

=== APIDOC-107 · Linnet behavioral Rust examples

*Priority:* P1 · *Status:* Open · *Surface:* `linnet`

After APIDOC-002, document graph, builder, storage, and subgraph invariants; explain
compound algorithm returns and mutation invalidation; and replace ignored or import-only
examples with parse → algorithm → serialize and mutation workflows.

*Complete when:* the examples run under the documented features and assert graph
properties rather than merely resolving names.

=== APIDOC-108 · Linnest Rust support boundary

*Priority:* P1 · *Status:* Open · *Surface:* `linnest`

Linnest is excluded from published Rust components but publicly re-exports byte/archive
operations and types without crate orientation or examples.

*Complete when:* the surface is explicitly classified as an internal Typst/WASM bridge
or a supported Rust API. A supported surface must document wire compatibility, schemas,
determinism, errors/panics, and a checked byte round trip.

== Python reference gaps

=== APIDOC-201 · GammaLoop integration-result boundary

*Priority:* P1 · *Status:* Open · *Surface:* `gammaloop-python`

Integration output, estimates, statistics, and table records are public, unsupported,
and largely undocumented, while the authored guide says Python has no supported
structured integration producer.

*Complete when:* either a supported producer and full result contract exist, including
units and uncertainties, or these records are removed from the curated public surface
and clearly marked provisional at their source boundary.

=== APIDOC-202 · GammaLoop evaluation and event journeys

*Priority:* P1 · *Status:* Open · *Surface:* `gammaloop-python`

Current examples rely on an unspecified `./state`, undefined points, and no expected
values. Events, correlated event groups, histogram accumulators, snapshots, merging,
and rebinning have no complete examples.

*Statistical-contract finding (2026-08-21):* each call to `fill_continuous_sample` or
`fill_discrete_sample` increments the sample count once, but multiple entries landing in the same
bin are currently sent to the bin accumulator separately. Their squared weights therefore add as
$w_1^2 + w_2^2$ rather than $(w_1 + w_2)^2$. Native observable processing first groups all
same-bin contributions from the complete correlated event-group list. The Python helper cannot yet
faithfully reproduce that path from raw correlated entries.

*Progress (`qtvpkmto`):* `evaluate_sample`, `evaluate_samples`, and `get_integrand_info` now have
self-contained repository-fixture examples with structural assertions, and the maintained Python
event script checks graph groups, event groups, histogram shapes, and sample counts. A separate
distinct-bin accumulator example verifies pending-state merge, raw sums, squared sums, and rebinning
without overstating the unsafe correlated same-bin path. A discoverable integration-space
dimension, a scheduled known-value evaluation, and the statistical-contract repair remain open.

*Complete when:* a shipped state supports two verified journeys: open → inspect
dimension → evaluate known point → interpret value/Jacobian/weight/stability; and batch
evaluate → preserve event groups → fill/commit/merge/rebin → assert a bin result. The manual
accumulator first matches native same-bin correlation semantics with a regression covering two
correlated entries in one bin.

=== APIDOC-203 · Linnet Python algorithms and proxies

*Priority:* P1 · *Status:* Open · *Surface:* `linnet-py`

The reference has no examples. It does not explain the spanning-forest result from
`cycle_basis`, the left/cut/right parts of `all_cuts`, traversal flags, callback
signatures, mutation failures, or proxy lifetime and write-through behavior.

*Complete when:* DOT parse/build → indexing → components/cycle basis → checked result
is runnable, and compound returns, callbacks, proxies, and mutation behavior have
local semantic contracts.

=== APIDOC-204 · Spynso standalone examples and member semantics

*Priority:* P1 · *Status:* Open · *Surface:* `spynso3`

Many example headings exist, but snippets depend on variables such as `evaluator`,
`large_input_batch`, `some_expression`, or `expr1`. Important dual, structure,
naming, bracket, and broadcast members do not state mutation/copy or compatibility
semantics.

*Complete when:* examples contain imports, setup, operation, and assertion; placeholder
variables are eliminated; and constructors/evaluators document shape, lookup,
compilation, mutation, and result-kind errors.

=== APIDOC-205 · Idenso Python transformation coverage

*Priority:* P1 · *Status:* In progress · *Surface:* `idenso-community`

`dirac_adjoint`, bispinor/color/metric/Minkowski expansion, initialization, and color
simplification lack checked examples. Existing authored examples are syntax-compiled
only.

*Progress (`umvyopkm`):* `dirac_adjoint`, initialization, and all five selective expansion
families now have self-contained imports, registered representations, operations, and assertions.
Their prose distinguishes distribution from component substitution and from later metric, Dirac,
or color simplification. Runtime execution and the combined color → Dirac HEP journey remain.

*Complete when:* a convention-explicit HEP expression passes through color and Dirac
rewrite stages with asserted free indices and expected identities, and each remaining
transformation has at least one runtime-verified example.

=== APIDOC-206 · Vakint Python journeys and failures

*Priority:* P1 · *Status:* Open · *Surface:* `vakint-community`

No declaration provides a useful raises contract, and examples are fragments. Parsing,
unknown topology, missing external tool, backend failure, missing numerical parameters,
and comparison failure are not distinguished.

*Complete when:* one matching-only journey asserts an exact topology, one configured
evaluation journey checks Laurent coefficients, and public operations document their
exception conditions.

== Cross-reference infrastructure gaps

=== APIDOC-301 · Curated support boundaries

*Priority:* P1 · *Status:* Open

GammaLoop and Spynso use curated support lists, while Linnet, Idenso, and Vakint
currently classify every exported Python item as supported. Some Rust implementation
surfaces are published without an equivalent stability decision.

*Complete when:* every generated reference item is classified as supported,
provisional, or internal from implementation-owned metadata, and authored interface
guides use the same boundary.

=== APIDOC-302 · Structured, verified examples

*Priority:* P1 · *Status:* Open

Python docstring examples render as prose and populate no structured example field.
Rust catalog examples often prove only imports. Signature checks therefore cannot
distinguish a runnable journey from a fragment.

*Complete when:* supported items can attach structured examples with a verification
tier, dependencies, expected result, and source location; publication reports example
coverage separately from signature coverage; and selected no-license examples execute
in the pull-request tier.

=== APIDOC-303 · Semantic completeness gate

*Priority:* P2 · *Status:* Open

Current generation rejects some missing descriptions but does not require conventions,
returns, failures, side effects, or useful examples. Enabling strict missing-docs checks
across every exposed implementation item immediately would create undifferentiated
noise.

*Complete when:* the curated supported surface is checked first against the completion
standard in this register, with actionable reporting by component and gap category.
Coverage of provisional and internal surfaces may then be expanded deliberately.
