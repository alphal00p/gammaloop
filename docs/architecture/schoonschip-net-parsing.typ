= Schoonschip network architecture
<schoonschip-network-parsing>

The original note was audited on 2026-09-14. This revision records source review
on 2026-10-27 against `4a9ff115`; no tests were rerun for this note update.

Spenso owns parsing, logical tensor interfaces and the operation graph. Idenso
assigns metric/vector and selected tensor identities to that graph through its
shared typed reduction engine.

== Entry Points
<entry-points>

#link("../../crates/idenso/src/tensor/contract.rs")[`SymbolicTensor<PartialStructure>::contract(ContractSettings)`]
performs structural contractions and optional chain/trace collection.
`simplify_algebra(AlgebraSettings)` selects gamma, color and epsilon identities
and additional structural work. Both return `Result<Self, TensorInferenceError>`
and preserve the tensor's logical interface. A zero can retain rank and port
order independently of its Atom payload.

== Network Algorithm
<network-algorithm>

Both APIs lower to the same outside-in domain planner. Borrowed shallow Spenso
graphs keep tensorial leaves opaque; Idenso's factorized slot contractor operates
on selected connections without distributing unrelated scalar factors. Changed
regions are reconsidered, while observations belong to the same immutable
payload and interface.

=== Driver outline
<driver-outline>

```text
checked symbolic tensor + logical interface
  -> contract(settings) or simplify_algebra(settings)
  -> shared capability-driven domain planner
  -> selected structural connections and enabled identity kernels
  -> typed tensor with exact payload and policy-relative reduction status
```

This is the maintained architecture, rather than a whole-expression loop over
removed raw-Atom network methods.

== Settings
<settings>

`ContractSettings` selects an optional factor order, representation filters,
metric/rank-one contraction, chain/trace collection and local expansion.
`expand = false` retains independent sum alternatives while allowing permitted
port substitutions. Trace collection forms notation and does not evaluate it.

`AlgebraSettings` independently selects gamma, color and epsilon capabilities.
`AlgebraContraction::{Fully, Dots, Minimal, Selected, None}` controls extra
structural work and output form. An enabled identity still permits the
contractions required by that identity, including under `Minimal` or `None`.
Both settings types have optional `max_passes`; the default is unlimited.

== Contraction Step
<contraction-step>

The structural contractor uses the admitted graph's connections and occurrence
interfaces. Representation filters preserve other ports. Selected contraction
work retains factorization outside the affected arithmetic scopes; local
distribution required by an enabled identity is distinct from expanding the
complete numerator. The logical interface remains attached to the result.

== Pattern Path
<pattern-path>

`shorthands/schoonschip/api.rs` and its raw `Schoonschip` trait are compiled only
as unit-test oracles. Removed `schoonschip_net`, `schoonschip_with_net`,
`NetworkSchoonschip` and `Schoonschipify` are not live public APIs.
`SchoonschipSettings` is crate-private configuration for a local ordered
metric/vector fallback: metrics, representation selection, chain-like functions
and rank-one handling. It has no public traversal, depth or contraction-order modes.

== Preserved structure and ownership
<why-the-network-path-should-be-better>

=== Structure before syntax
<structure-before-syntax>

A tensor retains its logical interface even when its expression is zero.
Checked constructors and changed reconstructions validate the payload against
that interface; printed slots alone cannot recover a ranked zero.

=== Controlled order
<controlled-order>

Factor order and representation selection are explicit contraction policy.
Neither a particular order nor retaining a shallow graph is a performance claim;
measure the actual workload and report the settings used.

=== Generic indices
<generic-indices>

The ordered specialization supports native index-aware network execution.
`PartialStructure` additionally retains logical port order and occurrence-local
unresolved axes for typed composition. Operations must transport that interface,
rather than infer anonymous axis order from a reconstructed Atom.

=== Separate operations
<less-double-meaning>

`to_dots` normalizes compact scalar-product notation. `undo_dots`, `undo_chain`
and `undo_trace` expose the corresponding notation with fresh compatible dummy
indices. These operations do not implicitly request full algebra reduction.
The parser no longer offers `Fast`/`Expanded` opaque-inference modes; shared
construction and observed interfaces own leaf admission.

== Examples
<examples>

Use `contract` to select structural work and chain/trace collection. Use
`simplify_algebra` with explicit family settings when gamma, color or epsilon
identities are required. Inspect the returned tensor before converting its
payload to a plain Atom, particularly for zero or unresolved logical axes.

== Expected Outcome
<expected-outcome>

`reduction_status()` reports `Complete`, `Deferred` or `Capped` relative to the
requested capabilities and policy. Deferred or capped results retain the exact
remaining expression. Successful return does not mean every unsupported tensor
head or identity has been eliminated. The maintained entry points are the typed
tensor's `contract` and `simplify_algebra`, not a second raw-Atom network pipeline.
