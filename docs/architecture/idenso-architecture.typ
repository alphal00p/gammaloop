= Idenso implementation architecture

#quote(block: true)[
#strong[Status:] Current implementation architecture; source reviewed on 2026-10-27
against `4a9ff115`. The original source audit was dated 2026-09-14. This update did
not rerun the validation commands described below.

This note covers the `idenso` Rust crate: representation-aware Symbolica transformations built
on Spenso's tensor syntax and network parser. It does not describe concrete tensor-component
evaluation, which belongs to Spenso and its tensor libraries.
]

== Boundaries

Idenso is an in-process symbolic identity layer. Its public modules divide responsibility as
follows:

- `representations` and `rep_symbols` define spin, bispinor, and color representations plus the
  namespaced Symbolica symbols used by patterns;
- `tensor` adapts a Symbolica expression plus an ordered Spenso structure into `SymbolicTensor`
  and `SymbolicNet`;
- `shorthands` owns chain, trace, dot-product, metric, and Schoonschip-style normalization;
- `dirac`, `color`, and `epsilon` own domain-specific identities and settings;
- `IndexTooling` and `cook` provide expression preparation, index hygiene,
  compact symbol encodings, and canonicalization;
- `spynso3` owns the unified `symbolica.community.tensor` Python facade; its
  `TensorExpression` wraps the shared native symbolic tensor.

Idenso always depends on Spenso with `shadowing`, as well as Linnet and Symbolica. It reuses
Spenso's representation tags, slot parser, tensor-network graph, and contraction scheduling;
it does not fork those abstractions. It has no FORM process, numerical backend, filesystem
workspace, network service, or database.

== Symbol and representation boundary

`initialize()` forces registration of Spenso's standard Lorentz-family representations,
Idenso's `SpinFundamental`, `Bispinor`, `ColorFundamental`, `ColorSextet`, and `ColorAdjoint`
types, their duals, and the metric, gamma, epsilon, and color symbols used by rewrite rules. The
crate-level Symbolica initializer performs the same registration when the community module is
loaded.

These representation types implement Spenso's `RepName` contract through
`SimpleRepresentation`: their Symbolica spelling, dual orientation, ordering, and tags become
the structural vocabulary recognized by Spenso parsing. Rewrite patterns rely on those tags and
namespaced symbols. A visually similar untagged function is therefore not automatically the
same tensor object.

Symbolica owns the atom data and its process-wide symbol registry. Idenso's lazily initialized
symbol sets borrow that registry; callers should initialize before parsing or applying rules and
must not treat printed names alone as a complete serialized registry.

== Core symbolic representation

`SymbolicTensor<S, E = Atom>` owns a structure `S` and a symbolic payload `E`.
The default structure is `OrderedStructure<LibraryRep, AbstractIndex>`; explicit
ordered tensors support network execution, while `SymbolicTensor<PartialStructure>`
retains logical port order and occurrence-local unresolved indices for composition.
An owned `Atom` is the usual payload; borrowed atom views use the same carrier.
The `is_metric` and `is_composite` operation flags are part of value identity.
Private proof and observation caches belong to that immutable payload/interface;
changed or untrusted reconstructions establish or invalidate the applicable facts.

The interface remains meaningful when the expression is zero. Checked constructors
validate nonzero payloads against their declared ports, and typed operations carry
the logical interface through composition, rewriting and conjugation. Callers inspect
`expression()` and `structure()`; a plain Atom alone cannot encode the rank of zero.

`SymbolicNet<Aind, E = Atom>` reuses Spenso's graph and store. The ordered specialization
implements direct structural contraction by multiplying the symbolic payloads and
merging their interfaces. Shared typed `contract` and `simplify_algebra` operations
return a tensor, retaining its interface and exact unprocessed expression.

== Rewrite families and execution flow

Public algebra reduction is selected through `AlgebraSettings` on
`SymbolicTensor<PartialStructure>::simplify_algebra`. Gamma, color and epsilon
capabilities are explicit; each enabled identity permits its structural prerequisites.
`AlgebraContraction` selects additional structural work and desired output notation.
The same domain planner also serves `contract(ContractSettings)`. It reconsiders
changed regions and retains unrelated scalar factorization rather than imposing a
fixed sequence of whole-expression simplifiers.

`IndexTooling` still supplies raw-Atom canonicalization, conjugation and index hygiene.
Canonicalization admits the full source before antisymmetric-zero pruning, preserves
representation duality wrappers, and reconstructs the validated graph without executing
unrequested tensor algebra. Closed scalar powers are scoped independently before
index labeling. `wrap_indices` uses `spenso::index_scope(header, index)` so wrapped
indices remain valid native abstract indices.

Cooking is separate from algebra selection. `CookSettings` controls function, index and
dimension payload encodings, source filters and output tags. Notation conversions
`to_dots`, `undo_dots`, `undo_chain` and `undo_trace` are explicit tensor operations;
collecting a trace does not evaluate it. Dimension-gated gamma and epsilon identities
continue to require the appropriate representation dimension.

== Schoonschip contraction and the shared planner

`contract(ContractSettings)` and `simplify_algebra(AlgebraSettings)` use the same
outside-in domain planner. Structural work reads a borrowed, shallow Spenso network
with opaque tensorial leaves and uses Idenso's factorized slot contractor. A contraction
order or representation filter selects eligible work; unaffected factors stay intact.
`expand = false` retains independent sum alternatives while allowing permitted index
substitutions. Algebra identities may still perform distribution required by that identity.

`max_passes` is an optional explicit budget; both public settings types default to no
limit. `ReductionStatus::{Complete, Deferred, Capped}` reports completion relative
to the requested policy without dropping the exact remaining expression.

The old raw `Schoonschip` trait in `shorthands/schoonschip/api.rs` is now compiled
only for unit-test oracle use. `SchoonschipSettings` is crate-private configuration
for a local metric/vector fallback; it no longer owns public traversal, parser-depth
or contraction-order modes. The detailed current boundary is described in the
#link("schoonschip-net-parsing.typ")[Schoonschip contraction note].

== Parsing, shorthand, and index invariants

Spenso representation syntax is the shared contract. Slots carry a representation, dimension,
abstract index, and dual orientation. Products contract matching dual slots; additions must
expose compatible external structure. `chain`, `trace`, `dot`, metrics, and bracket-like syntax
are parser-owned shorthands rather than arbitrary opaque functions.

The explicit `undo_dots`, `undo_chain` and `undo_trace` operations select the corresponding
Spenso expansion mode and reconstruct that notation with fresh parse-local dummies. The shared
contraction engine instead keeps tensorial leaves opaque and validates their observed interfaces;
there is no `Fast`/`Expanded` opaque-inference selector. Construction, notation expansion and
algebra reduction are separate policies.

Index wrapping is an ownership operation. `list_dangling` discovers external slots from a parsed
network; `wrap_dummies` changes only non-external index payloads; `wrap_indices` changes every
recognized slot. Independently created expressions should be wrapped before multiplication when
same-spelled dummy names must not contract.

The shared concrete syntax and which crate owns each rewrite are specified in the
#link("spenso-symbolica-syntax-and-rewrites.typ")[Spenso/Idenso Symbolica syntax note].

== Features and serialization

Idenso defaults to `native`, forwarding GMP/MPFR support; `wasm` selects the Wasm
backend with default features disabled. The Rust algebra layer always uses Spenso's
`shadowing` support. `reference-cases` exposes curated identity fixtures. Python and
stub generation belong to `spynso3`'s unified tensor facade; Idenso has no separate
`python` or `python_stubgen` Cargo feature or community module.

With `bincode`, `SymbolicTensor<PartialStructure>` encodes its Atom payload, logical
slot order, and metric/composite flags. Decode uses Symbolica's state map and validates
the payload/interface before returning a tensor. Proof, observation and planner caches
are never serialized; this is not a versioned workspace or cross-version checkpoint promise.

Reversible function/index cooking normally stores an Atom in Symbolica symbol user data
and requires the corresponding registry for restoration. Compound dimension cooking is
different: reversible mode uses the owner's portable textual dimension encoding, accepts
only supported exact rational arithmetic, and rejects function/callback-dependent or
metadata-bearing dimensions; the registered scalar color parameter `Nc` remains
  explicitly supported. Ordinary atomic dimensions stay unchanged.

== Ownership and error boundaries

The shared symbolic tensor owns its logical interface and payload; settings are borrowed
for a transformation. Immutable tensors may retain validated interface facts, regional
observations and a shallow graph. A changed result must invalidate or re-establish those
facts rather than treating cached observations as independent authority.

Typed construction, composition, contraction and algebra operations return
`TensorInferenceError` on incompatible interfaces or unsupported admission. Raw network
parsing/execution retains `NetworkToolingError`; adjunction, cooking, canonicalization
and index tooling preserve their respective typed errors. Unmatched identity syntax may
remain exactly represented. Inspect `reduction_status()` and the remaining expression;
a successful call alone is not a promise that every possible algebraic identity was applied.

== Maintained invariants

- A tensor's logical interface is authoritative for ranked zero and unresolved ports;
  nonzero payloads must agree with the checked interface.
- Representation, dimension, dual orientation and logical port order matter independently
  of printed function names. Index transformations transport the interface and payload together.
- Dummy allocation is scoped; independent expressions need explicit index hygiene.
- Four-dimensional identities require an established four-dimensional representation.
- Planner completion is relative to the requested capabilities and output policy. Deferred
  or budget-capped work preserves the exact remaining expression.
- Function/index user-data cooking requires the matching symbol registry. Portable cooked
  dimensions follow their separate validated textual encoding contract.

== Verification and related documentation

Tests live beside tensor parsing and canonicalization, cooking, index tooling, typed
composition and interface validation, contraction planning, notation conversion, and the
Dirac, color and epsilon kernels. Curated FORM and FeynCalc examples are reference fixtures,
not runtime calls to those systems. The benchmark sources remain separate from the
source-review evidence recorded here.

The native boundary is exercised with `cargo test -p idenso`; optional tensor/representation
encoding uses `cargo test -p idenso --features bincode`. Python and stub coverage belong to the
shared `spynso3` host. `reference-cases` makes curated fixtures available to non-test consumers.
These commands describe coverage entry points; this source review did not execute them.

For supported workflows, start with the
#link("../../../products/idenso/latest/tutorial/")[controlled identity tutorial], continue with
the #link("../../../products/idenso/latest/guides/algebra/")[algebra and index-hygiene guide],
and consult the
#link("../../../products/idenso/latest/reference/form-color-dirac/")[source-backed rule
specification]. Exact public signatures are in the
#link("../../../products/idenso/latest/reference/rust/idenso/")[native Idenso Rustdoc]
and the
#link("../../../products/idenso/latest/reference/python/spynso3/")[shared tensor Python
module reference].
