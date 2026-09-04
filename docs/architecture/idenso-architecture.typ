= Idenso implementation architecture

#quote(block: true)[
#strong[Status:] Current implementation architecture, audited against the Idenso source on
2026-08-18.

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
- `selective_expand`, `IndexTooling`, and `cook` provide expression preparation, index hygiene,
  compact symbol encodings, and canonicalization;
- the optional `python` module registers a Symbolica community-module facade.

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

`SymbolicTensor<Aind>` carries four fields:

- an `OrderedStructure<LibraryRep, Aind>` containing the external slots;
- the owned Symbolica `Atom` expression;
- `is_metric`, used by contraction-specialized identities;
- `is_composite`, distinguishing a direct tensor function from an expression-backed leaf.

It implements both `TensorStructure` and `HasStructure`. Structural methods delegate to the
ordered slots, while permutations rewrite slot atoms and introduce identity tensors when
needed. Its generic `Contract` implementation multiplies the two atoms and merges their
structures; domain-specific cleanup is a later rewrite stage rather than component-wise
evaluation.

`SymbolicNet<Aind>` is a Spenso `Network` whose local tensors are `SymbolicTensor`, whose scalars
are Symbolica atoms, and whose function keys are Symbolica symbols. `SymbolicNetParse` forces
Spenso's `ContainsReps` tensor filter, so representation-bearing functions can become tensor
leaves even without a separate tensor tag. `SymbolicNetExt::simple_execute` uses sequential
Spenso execution and reconstructs the final atom from zero, one, or the result tensor.

This representation is deliberately ephemeral. It is a topology-aware vehicle for parsing,
canonicalization, and contractions; the authoritative public result of Idenso transformations is
normally an `Atom`.

== Rewrite families and execution flow

Idenso does not expose one universal simplification pipeline. Each family selects and owns its
normal form:

- `IndexTooling` parses enough structure to canonicalize tensor indices, list dangling slots,
  wrap all indices or only dummies, form adjoints, and alias tensor subexpressions;
- `SelectiveExpand` distributes only factors carrying selected representation families, keeping
  unrelated sectors factorized;
- `MetricSimplifier` and `Chain` normalize metrics, dot products, open chains, and closed traces;
- `GammaSimplifier` collects bispinor chains and applies dimension-gated Clifford, trace,
  projector, gamma5, and optional four-dimensional epsilon rules;
- `ColorSimplifier` collects fundamental color lines, closes traces, applies generator,
  structure-constant, Fierz, and Casimir rules, prunes antisymmetric zero terms, and iterates to a
  fixed point;
- `EpsilonSimplifier` owns epsilon contractions and reductions;
- `Cookable` replaces selected functions or representation-index payloads with compact symbols,
  either as readable flattened names or reversible Symbolica `UserData::Atom` encodings.

Settings objects are part of the semantics. Gamma ordering and trace evaluation, color Fierz and
invariant substitutions, Schoonschip depth and traversal, and cooking source/tag filters can all
change the result form. Reproducible callers should record the exact settings and the order in
which independent rewrite families ran.

Most pattern engines use a local fixed-point loop: transform the current atom, compare it with
the previous atom, and stop when unchanged. That guarantees termination only for the implemented
oriented rule set; custom rules composed by a caller remain the caller's responsibility.

== Schoonschip network path

The Schoonschip-style simplifier is the main bridge from domain rewrites to Spenso network
execution:

```text
Atom or one term of an Add
  -> normalize dot-product syntax
  -> parse a SymbolicNet with opaque shorthand leaves
  -> limit parser depth and composite-scalar recursion from settings
  -> execute sequentially with a Schoonschip contraction strategy
  -> choose tensor pairs by degree or expression-size policy
  -> simplify the contracted symbolic boundary
  -> extract an Atom and normalize dots again
  -> repeat until unchanged, unless SinglePass was selected
```

Additions are processed term by term. `SchoonschipSettings` chooses single-pass or recursive
depth-/breadth-first traversal, parse depth, contracted-sum expansion, chain-like function and
rank-one handling, and one of several topology- or expression-size contraction orders. The
custom `Schoonschipify` strategies plug into Spenso's `ContractionStrategy`; they do not replace
Spenso's graph or store.

The direct pattern path exposed by `schoonschip_with_settings` applies the local metric/function
replacement families without building the full network. The network path is preferable when
contraction topology and ordering matter. Their exact behavior and settings are detailed in the
#link("schoonschip-net-parsing.typ")[Schoonschip network parsing note].

== Parsing, shorthand, and index invariants

Spenso representation syntax is the shared contract. Slots carry a representation, dimension,
abstract index, and dual orientation. Products contract matching dual slots; additions must
expose compatible external structure. `chain`, `trace`, `dot`, metrics, and bracket-like syntax
are parser-owned shorthands rather than arbitrary opaque functions.

`UndoShorthands` selects Spenso `ShorthandParsing::Expand` modes, parses a symbolic network, and
executes it to reconstruct explicit tensor products with fresh parse-local dummies. The
Schoonschip path instead requests opaque shorthand with fast structure inference so it can
simplify contraction boundaries selectively. These two modes are intentionally different:
expansion exposes topology but can grow expressions, while opaque inference preserves compact
syntax and trusts the inferred external slots.

Index wrapping is an ownership operation. `list_dangling` discovers external slots from a parsed
network; `wrap_dummies` changes only non-external index payloads; `wrap_indices` changes every
recognized slot. Independently created expressions should be wrapped before multiplication when
same-spelled dummy names must not contract.

The shared concrete syntax and which crate owns each rewrite are specified in the
#link("spenso-symbolica-syntax-and-rewrites.typ")[Spenso/Idenso Symbolica syntax note].

== Features and serialization

Idenso has no default Cargo features. The core Rust rewrite layer still always includes
Symbolica and Spenso's `shadowing` support. `python` enables the community-module functions and
representation initializer; `python_stubgen` adds stub metadata and enables Symbolica's stub
surface. `reference-cases` exposes the otherwise test-only curated identity cases.

The optional `bincode` feature derives binary encoding only for Idenso's zero-sized
representation marker types. `SymbolicTensor`, rewrite settings, symbol registries, networks,
and cooked expressions do not form an Idenso checkpoint format. Reversible cooking stores the
source atom in Symbolica symbol user data and derives a stable-looking hash name, but recovery
still depends on the matching Symbolica registry and cooking settings. Printed or binary names
alone are not a portable physics-result archive.

== Ownership and error boundaries

Most high-level identity traits take an `Atom` or `AtomView` and return a newly owned `Atom`;
Symbolica owns internal expression memory and symbol metadata. Settings are borrowed for one
transformation. Parsed `SymbolicNet` values, dummy libraries, and contraction plans are local to
the call and are discarded after atom extraction.

Fallible structural entry points remain visible: `SymbolicNetParse` returns
`TensorNetworkError`, structure inference returns `StructureError`, `dirac_adjoint` returns an
`eyre::Result`, and `CookSettings::try_cook`/`try_cook_indices` report `CookingError`. By
contrast, convenience APIs such as `schoonschip`, shorthand expansion, dangling-index queries,
and `simple_execute` unwrap parsing or network-execution errors internally. They assume an
initialized, Spenso-compatible expression and may panic on malformed or inconsistent tensor
syntax; use the fallible parsing boundary first when inputs are not trusted.

Domain simplifiers usually leave unmatched syntax unchanged rather than diagnosing it as an
error. A successful return therefore means the configured rewrite reached its fixed point, not
that every physics object was recognized or eliminated. Verification should inspect residual
representation-bearing factors when completeness matters.

== Maintained invariants

- Tensor structure is inferred from tagged representation arguments, including dimension and
  dual orientation; function spelling alone is insufficient.
- `SymbolicTensor.structure` and the slot atoms inside `expression` must describe the same
  external indices. Permutation operations rewrite both sides together.
- Dummy indices introduced by parsing or shorthand expansion are fresh only within that parse
  state. Cross-expression hygiene remains explicit through wrapping or canonicalization.
- Dimension-specific identities fire only when the required representation dimension can be
  established; four-dimensional gamma/epsilon rules are not generic-dimensional rules.
- Fixed-point simplifiers compare whole Symbolica atoms, so canonical ordering and normalization
  are part of their termination contract.
- A reversible cooked symbol is meaningful only with its Symbolica user data and matching tag
  policy; a flattened cooked name intentionally cannot reconstruct its source.
- Schoonschip contraction order is semantic/performance policy over one parsed network. Changing
  it must not be confused with changing the parser's shorthand or depth policy.

== Verification and related documentation

Tests live beside tensor parsing and canonicalization, cooking, index tooling, shorthand
expansion, Schoonschip traversal and contraction orders, metric/chain normalization, and the
Dirac, color, and epsilon rules. Snapshot tests pin canonical Symbolica strings. Curated FORM
and FeynCalc examples are reference fixtures checked by tests; they are not runtime calls to
those external systems. Benchmarks separately cover Schoonschip modes and vertex-algebra paths.

The default boundary is exercised with `cargo test -p idenso`. Optional representation encoding
uses `cargo test -p idenso --features bincode`; community-module and stub coverage require their
respective Python features. The `reference-cases` feature makes the curated cases available to
non-test consumers but does not add a second simplifier.

For supported workflows, start with the
#link("../../../products/idenso/latest/tutorial/")[controlled identity tutorial], continue with
the #link("../../../products/idenso/latest/guides/algebra/")[algebra and index-hygiene guide],
and consult the
#link("../../../products/idenso/latest/reference/form-color-dirac/")[source-backed rule
specification]. Exact public signatures are in the
#link("../../../products/idenso/latest/reference/rust/idenso/")[native Idenso Rustdoc]
and the
#link("../../../products/idenso/latest/reference/python/idenso-community/")[Python community
module reference].
