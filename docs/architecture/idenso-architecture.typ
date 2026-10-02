= Idenso implementation architecture

#quote(block: true)[
#strong[Status:] Current implementation architecture, audited against the Idenso source on
2026-09-14.

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
- typed collection, `IndexTooling`, and `cook` provide expression preparation, index hygiene,
  compact symbol encodings, and canonicalization;
- `spynso3` owns the unified Python facade, exposing Idenso algebra on `TensorExpression`.

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

`contract(ContractSettings)` and `simplify_algebra(AlgebraSettings)` lower into one
domain scheduler. Their settings are independent: structural permissions and
connection filters belong to contraction, while identity families, conventions
and algebraic output forms belong to algebra reduction. An enabled identity
authorizes its required contractions without enabling unrelated identities.
Each family retains its internal kernel:

- `IndexTooling` parses enough structure to canonicalize tensor indices, list dangling slots,
  wrap all indices or only dummies, and form adjoints;
- Typed `collect` and `coefficient_list` use the shared graph/tape owner to retain selected
  representation sectors while keeping unrelated coefficients factored;
- the shared component contractor and internal chain/dot normalization own metric
  contraction, compact scalar products, open chains and closed traces;
- the Dirac pass collects bispinor chains and applies dimension-gated Clifford, trace,
  projector, gamma5, and optional four-dimensional epsilon rules;
- `ColorSimplifier` collects fundamental color lines, closes traces, applies generator,
  structure-constant, Fierz, and Casimir rules, prunes antisymmetric zero terms, and iterates to a
  fixed point;
- the epsilon kernel owns epsilon contractions and reductions;
- `Cookable` replaces selected functions or representation-index payloads with compact symbols,
  either as readable flattened names or reversible Symbolica `UserData::Atom` encodings.

Settings objects are part of the semantics. Gamma ordering and trace evaluation, color Fierz and
invariant substitutions, contraction budgets, and cooking source/tag filters can all
change the result form. Reproducible callers record the exact settings; the
planner chooses eligible work using candidate counts, prerequisites and estimated
growth. The expression receives one initial candidate scan. Regional replacements
update occurrence counts; unchanged regions are not rescanned. Per-family settled
regions and explicit change flags drive local fixed points.
Completed identical scope observations share immutable storage even during the
initial walk. Each occurrence is still inspected and contributes its counts;
sharing does not skip admission or change the planner's budget. Parent records
retain the shared completed record instead of duplicate observation subtrees.
The same observations retain explicit metric and vector source ports by
representation. Their absence can certify a completed structural fallback even
when compact vectors or repeated gamma ports remain. Head occurrence alone
does not establish eligible work. Validated external-port counts also exclude
false connections between same-variance or dimension-incompatible ports.
The observation records the maximum explicit-index count across sum alternatives
and adds counts across product factors. When that maximum equals the established
external boundary, no branch contains an internal pair. This fact survives
regional updates without rescanning an unchanged sum; powered scopes retain
their separate multiplicity checks.

Algebra prerequisites use these same observations to discharge free-port
regions before contraction setup or budget accounting. Unselected scalar
spectators and unrelated tensor factors remain opaque to the collector.

Positive powers of a dot with a native fixed-component selector use the same
intrinsic normalization proof when their scope is copied. An enclosing scalar
spectator therefore does not turn that supported coefficient collection into
a refusal. Callback-sensitive descendants retain the existing refusal path;
this changes neither scalar zero-proof policy nor copy-local index ownership.

The same domain retains its depth-one graph and exact opaque-leaf interfaces.
When a trusted intrinsic rewrite changes topology, unchanged leaves reuse their
logical ports through Spenso's admitted parser entry. That entry skips redundant
admission and dummy-reservation walks; only newly selected boundaries need
interface discovery. Callback-sensitive rewrites discard invalidated interface
facts and retain checked result boundaries.

Identity kernels may publish a terminal polynomial without enumerating every
new scalar coefficient. A generated power such as $M^2$ therefore need not
have a regional observation yet. When a later contraction encounters that
power, the existing scalar-spectator certificate checks its intrinsic metadata,
exact coefficients and empty interface, then retains it as one opaque leaf.
It does not re-admit the whole tensor or invalidate the domain cache. Indexed,
callback-sensitive and inexact powers retain their existing scoped or deferred
paths, and transformation budgets remain unchanged.

== One symbolic contraction path

Prerequisite regions carry the enclosing carrier's conservative index history
through their local contraction. They return their payload through
`with_identity_result`, retaining its release-time checks for callback-sensitive
publication; region construction does not re-admit the enclosing tensor.

The tensor contractor reads the shared structural graph once, keeping opaque
leaves and their logical ports. It selects an incident factor, visits its
alternatives and applies port substitutions through the existing component
reducer. Equal remaining states merge before result emission. Metrics alone
and metrics with rank-one tensors use this same owner. Powers whose bases
contain local dummy pairs finish that scalar scope before the power is applied.

An explicit metric port may be renamed through native unary conjugation while
preserving the conjugation wrapper, coefficients and slot variance. This does
not permit vector or scalar substitution through conjugation, and arbitrary
broadcast functions remain opaque. Refused vector substitutions do not invoke
their normalizers. Conjugation and contraction remain separate operations.

Callback-sensitive syntax uses the existing ordered substitution schedule and
checked result construction. It does not acquire a trusted interface merely
because a substitution was algebraically a contraction. The example
`g(a,b)*T(a)` whose normalizer turns the free `T(b)` into a scalar violates
the normalization contract and remains a rejection on the checked contractor path.

Ordered substitutions publish the checked surviving interface, including when
the ordered fallback exposes more work for the component reducer. Discovery
probes and refused linear branches do not publish speculative observations.
Exact imaginary coefficients remain scalar leaves of the same graph.

The retained repeated-index observations decide which factor boundaries can
contract. Disjoint open factors stay opaque; their sums are not visited merely because they
contain metrics. Sum alternatives and powers retain the scanner's established
scope rules. Public `ReductionStatus` distinguishes `Complete`, `Deferred` and
`Capped` relative to the requested work. Deferred work is checked again only
when its inputs change; an exhausted frontier is retained exactly. Neither
state is certified complete merely because its algebraic output is valid.

Epsilon simplification uses the same selected collector. It exposes epsilon
bodies and their incident metric/vector factors. Unrelated factors and
callback-sensitive templates remain opaque.

The four-dimensional Minkowski epsilon convention is fixed by
`Tr[gamma5 gamma^mu gamma^nu gamma^rho gamma^sigma] = 4 epsilon^{mu nu rho sigma}`.
With the mostly-minus metric and the existing Weyl gamma5 this means
`epsilon^{0123} = -i`. `contract_epsilon4` reuses the determinant kernel:
its rows are covector components, and four contravariant vector inputs require
lowering their spatial entries first. UFO's real epsilon convention is converted
at generator import by `-i` times Idenso epsilon. Concrete HEP gamma matrices and
symbolic gamma-zero identities keep their existing owners and conventions.
A vector-source sum can also contract through a partner sum when distribution
is permitted and every term exposes exactly one compatible slot. Unsupported
or ambiguous frontiers remain deferred; unrelated factors remain opaque.

`contract_ports` is the separate logical binary composition operation.
`to_dots`, `undo_dots`, `undo_chain` and `undo_trace` change symbolic notation.
`to_dots` converts surviving index-free compact scalar products and performs no
index contraction; undo operations introduce fresh compatible dummy scopes and
never evaluate traces. The intrinsic `g` normalizer exploits dot symmetry at
construction. A finite `expand_dots` request crosses Spenso's component
execution boundary and leaves unrelated factors untouched.

The former network Schoonschip engine, its traversal/order strategy family and
duplicate metric collector have been removed. The shared tensor scheduler runs
one fixed point over affected expression regions, carrying observations
and completion proofs when their invariants allow it. It retains the exact
completed settings for an unchanged expression. An actual expression
change invalidates those completion facts; callback boundaries retain
their existing invalidation rules. Status describes the current request, and
capped or deferred requests never enter the completed-settings cache. Budgets
retain exact unfinished work and do not certify it complete.

The graph, substitution and callback boundaries are detailed in the
#link("schoonschip-net-parsing.typ")[shared contraction architecture note].

== Parsing, shorthand, and index invariants

Spenso representation syntax is the shared contract. Slots carry a representation, dimension,
abstract index, and dual orientation. Products contract matching dual slots; additions must
expose compatible external structure. `chain`, `trace`, `dot`, metrics, and bracket-like syntax
are parser-owned shorthands rather than arbitrary opaque functions.

Index canonicalization temporarily gives composite dimensions such as `Nc^2-1`
reversible symbols so the tensor parser can recognize their slots. Each symbol
stores its defining expression inside an `idenso::representation_dimension`
wrapper in Symbolica user data. Canonicalization restores only those dimension
payloads afterward; separately cooked spectator functions remain cooked. This
preserves formal dimensions without expanding the graph numerator.

`UndoShorthands` selects Spenso `ShorthandParsing::Expand` modes, parses a symbolic network, and
executes it to reconstruct explicit tensor products with fresh parse-local dummies. The
shared contractor requests opaque shorthand with fast structure inference so it can
simplify contraction boundaries selectively. These two modes are intentionally different:
expansion exposes topology but can grow expressions, while opaque inference preserves compact
syntax and validates its declared boundary.

Compact inner products may retain explicit spectator ports. The shared slot
matcher identifies the contracted implicit axis in each operand; both ordered
and partial inference retain the other ports. Gamma chain assembly exposes only
inner products containing a selected gamma or chain endpoint, through Spenso's
existing shorthand materializer. Scalar dot coefficients and unrelated sectors
remain opaque.

Index wrapping is an ownership operation. `list_dangling` discovers external slots through the
shared parser's structure-only construction; `wrap_indices(scope, dummies_only=True)` scopes non-external explicit indices; its default scopes all explicit indices while leaving unresolved
open-port markers intact. Independently created expressions should be wrapped before multiplication when
same-spelled dummy names must not contract.

The shared concrete syntax and which crate owns each rewrite are specified in the
#link("spenso-symbolica-syntax-and-rewrites.typ")[Spenso/Idenso Symbolica syntax note].

== Features and serialization

Idenso defaults to `native`, forwarding GMP/MPFR support; `wasm` selects the Wasm backend with
default features disabled. The core Rust rewrite layer always includes Symbolica and Spenso's
`shadowing` support. The `spynso3` crate provides the Python methods and automatic representation
initialization; its `python_stubgen` feature adds metadata for the unified Spenso module.
`reference-cases` exposes the otherwise test-only curated identity cases.

Binary encoding is always derived for Idenso's zero-sized representation marker types.
`SymbolicTensor`, rewrite settings, symbol registries, networks,
and cooked expressions do not form an Idenso checkpoint format. Reversible cooking stores the
source atom in Symbolica symbol user data and derives a stable-looking hash name, but recovery
still depends on the matching Symbolica registry and cooking settings. Printed or binary names
alone are not a portable physics-result archive.

== Ownership and error boundaries

The public tensor operations retain the logical interface, typed zeros and
metadata in `SymbolicTensor`. Spynso converts arguments, dispatches to this
owner, translates errors and wraps the result. Symbolica owns atom memory,
normalization; Spenso owns graph topology,
logical layouts and component execution.

Structural entry points are fallible. Inference and certified tensor rewrites
report `TensorInferenceError`; network parsing and execution keep their typed
network errors. `dirac_adjoint` reports `AdjointError`, cooking reports
`CookingError`, canonicalization reports `CanonicalizationError`, and raw index
queries report `IndexToolingError`. A raw-expression escape hatch does not
retain a separately declared logical layout automatically. Callers restoring
that layout must supply the original structure.

Composition, inference and cooking preserve underlying admission, slot and
index errors as typed sources. Python bindings translate those errors at the
presentation boundary. The Bincode custom-error interface accepts a diagnostic
string, so formatting is confined to that serialization boundary.

Fallible structural entry points use typed errors: `SymbolicNetParse` returns
`TensorNetworkError`, structure inference returns `StructureError`, `dirac_adjoint` reports
`AdjointError`, and `CookSettings::try_cook`/`try_cook_indices` report `CookingError`.
Canonicalization returns `CanonicalizationError`; dangling-index queries and dummy wrapping
return `IndexToolingError`. Shorthand expansion, network Schoonschip passes, and
`simple_execute` return `Result<Atom, NetworkToolingError>` rather than unwrapping parser or
execution failures. The direct pattern-based `schoonschip` path still returns an `Atom`.

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
- Contraction completeness is an explicit certificate. A valid retained frontier, metric-only
  result, or unrecognized residual tensor is not automatically fully contracted.

== Verification and related documentation

Tests live beside tensor parsing and canonicalization, cooking, index tooling, shorthand
expansion, shared factor-graph contraction, metric/chain normalization, and the
Dirac, color, and epsilon rules. Snapshot tests pin canonical Symbolica strings. Curated FORM
and FeynCalc examples are reference fixtures checked by tests; they are not runtime calls to
those external systems. Benchmarks separately cover Schoonschip modes and vertex-algebra paths.

The default boundary is exercised with `cargo test -p idenso`. Optional representation encoding
uses `cargo test -p idenso --features bincode`; community-module and stub coverage live in
`spynso3` with its `python_stubgen` feature. The `reference-cases` feature makes the curated cases available to
non-test consumers but does not add a second simplifier.

For supported workflows, start with the
#link("../../../products/idenso/latest/tutorial/")[controlled identity tutorial], continue with
the #link("../../../products/idenso/latest/guides/algebra/")[algebra and index-hygiene guide],
and consult the
#link("../../../products/idenso/latest/reference/form-color-dirac/")[source-backed rule
specification]. Exact public signatures are in the
#link("../../../products/idenso/latest/reference/rust/idenso/")[native Idenso Rustdoc]
and the
#link("../../../products/idenso/latest/reference/python/spynso3/")[Python community
module reference].

== Factorized canonicalization and symbolic targets

Index canonicalization protects independent scalar factors and reserves explicit external indices
before allocating contraction dummies. Canceled representation groups do not consume canonical dummy
names. The fallible `Concretize` implementation preserves symbolic dimensions and uses the supplied
canonical layout to restore logical slot order in the symbolic expression.


== Historical M1 consolidation measurements

The historical M1 comparison used three alternating fresh-interpreter runs
against the immutable M0 release on CPU 8. The harness is
`examples/notebooks/fermion_ladder.py`. Each run has
one unmeasured warmup. The primary clock includes rule application, contraction
and explicit scalar materialization; FORM uses batched body CPU time, excluding
process startup. These measurements are pinned among cooperating jobs, not on an
exclusive host.

#table(columns: 4,
  table.header([Case], [M0 CPU ms], [M1 CPU ms], [FORM CPU ms]),
  [Historical ladder, early order], [451.95], [150.05], [265.00],
  [Historical ladder, original order], [1112.81], [753.63], [848.00],
  [Physical four-loop gluon, 4D], [763.43], [327.09], [759.00],
  [Physical four-loop gluon, D], [1200.43], [753.08], [1148.00],
  [Physical three-loop fermion, 4D], [0.685], [0.686], [0.360],
  [Physical three-loop fermion, D], [2.144], [2.588], [1.640],
  [Physical four-loop fermion, 4D], [5.592], [5.481], [6.000],
  [Physical four-loop fermion, D], [43.87], [43.94], [55.40],
  [Free length-14 trace, D], [416.37], [456.33], [61.33],
)

The plan's historical reverse/early order is `5,4,6,3,7,2,8,1`. Separate
diagnostic calls give 39.84 ms for contraction and 108.14 ms for explicit
materialization, with 180 literal alias definitions occupying 91,194 bytes. All
three M1 thresholds pass. Diagnostic phase medians are independent samples, so
their sum is not the primary end-to-end measurement. All ladder gains held in
every paired run. Traces retain their existing route at M1; the long free
D-dimensional trace remains a substantial gap and did not improve.

The cohort retains 102 timing records, 34 correctness/diagnostic records and
FORM scripts/results. Exact comparisons cover compatible polynomial bases;
HEP component evaluations and D-to-4 specialization cover the differing
four-dimensional trace bases. Component samples do not prove a general symbolic
identity. Source and installed-core identities remained unchanged throughout.
