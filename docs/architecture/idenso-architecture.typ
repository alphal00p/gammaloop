= Idenso implementation architecture

#quote(block: true)[
#strong[Status:] Current implementation architecture, audited against the Idenso source on
2026-09-29.

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

`SymbolicTensor<S, E = Atom>` owns the tensor interface and symbolic payload:

- the structure `S`, defaulting to `OrderedStructure<LibraryRep, AbstractIndex>`;
- the symbolic payload `E`, normally an owned Symbolica `Atom`;
- `is_metric`, used by contraction-specialized identities;
- `is_composite`, distinguishing a direct tensor function from an expression-backed leaf;
- private interface, observation and completion proofs, invalidated when the operation cannot
  justify carrying them to its result.

The explicit ordered specialization implements `TensorStructure` and network contraction.
Structural methods delegate to the ordered slots, while permutations rewrite slot atoms and
introduce identity tensors when needed. Its `Contract` implementation multiplies the atoms and
merges their structures; domain-specific cleanup follows separately. `HasStructure` supports
structure types implementing `TensorStructure`, and mapping the structure retains the
expression and classification flags.

The same symbolic tensor with `PartialStructure` retains canonical-to-logical layout and
occurrence-local unresolved ports. It owns positional composition, index substitution,
chain/trace assembly, and the fast interface inference used by Python. Unresolved ports remain
distinct under ordered multiplication, and a zero retains its declared shape. This replaces
the former Spynso `StructuredAtom`; Python bindings retain argument conversion, dispatch,
error translation, and presentation metadata. `PartialStructure` does not pretend to implement
the canonical-storage `TensorStructure` contract.

Canonical singleton `cind(n)` and `find(n)` markers select a concrete component,
so a rank-one leaf ending in one of these markers has a scalar interface. The
shared slot matcher validates a nonnegative integer component and its final
argument position; malformed markers, extra structural arguments and unrelated
functions with the same short name are rejected. This agrees with Spenso's
existing concrete-component interpretation without admitting arbitrary
zero-port vector syntax. The built-in basis vector `δ(cind(k), port)` is
an explicit exception for metadata: its first marker selects the basis component,
while its final argument is the sole tensor port. This does not admit a second
component marker on arbitrary rank-one functions.

Aliased results use the same tensor with Symbolica's `AliasedAtom` payload and
`AliasInterfaces` structure. Symbolica owns the single registry of literal
definitions; Idenso retains the root layout and each handle/body's declared
logical layout. Re-inferring those layouts from normalized syntax would lose
port order. Registration checks encoded interfaces, branch consistency, index
multiplicity, conflicting definitions and cycles. A fresh literal handle is
registered for each port labelling; a relabelled use is not a parametric lookup.
Resolution is explicit and checks the resulting interface, including callback
rank changes. Expansion is a separate explicit operation; the evaluator consumes
the alias registry directly. The Python `AliasedTensorExpression` only converts
typed arguments, invokes callbacks, retains presentation descriptors, and wraps
these shared operations. Scalar outputs use the same class.

`TensorRule` owns a Symbolica pattern, RHS, conditions and matching settings.
Wildcard closure and eligible literal RHS interfaces are checked once. A rule
whose wildcard bindings determine the tensor shape still checks each distinct
instantiated RHS and each actual target; equal bindings can match alternatives
with different interfaces. Matcher state and callback caches belong to one
application, so a reusable rule does not retain source views or suppress
callbacks across calls. A zero-sized RHS cache preserves per-match callbacks.

Factorized metric/vector contraction is owned by
`SymbolicTensor<PartialStructure>::contract`. It uses the existing component
reducer: only one selected factor contributes alternatives at a time, while
remaining factors carry occurrence-local port overrides. Equal remaining states
share a coefficient; identical coefficient bodies share a literal alias. Scalar
spectators, disconnected components and foreign sectors remain factorized.
The current state key retains original factor identity and port bindings, so it
can under-merge different bindings that normalize to equal factors. This is a
performance limitation, not an algebraic equivalence assumption.

Contraction, rule application and polynomial materialization are separate
operations. `AliasedTensorExpression.expand()` evaluates the existing definition
DAG forward with a variable table local to each definition. An additive root
merges its emitted summands in one bulk construction, avoiding a global dense
variable table. Alias uses inside opaque functions and callback-
sensitive bodies retain the ordinary checked resolution path before the
explicit expansion. An evaluator consumes the DAG directly.

Frontier growth is bounded before processing the next factor. The current guard
counts generated alternatives and estimates state/definition bytes; it is not a
strict heap cap on coefficient buffers. If the guard fires, the result retains
exact pending factors, and its internal completion flag remains false. It must
not be certified as fully contracted. Callback-sensitive inputs continue through
the existing checked contractor, including the rank-loss rejection for a
normalizer that turns `T(b)` in `g(a,b)*T(a)` into a scalar.

Checked reconstruction preserves established interfaces where the operation justifies it.
Public fields and low-level constructors are not validity certificates: inference can retain
raw repeated ports until `checked_parts` merges their explicit contractions. This intermediate
boundary also lets bindings discard stored-data identity when a declaration becomes a
contraction, even if its input and final result have the same rank.
Callback-sensitive rewrites validate their output because normalization can remove tensor
ports. Result validation observes the existing ports without temporary index materialization
or callback replay; constructor inference retains its materialization semantics within the
same inference owner. Index substitution distinguishes graph storage identities from logical unresolved
ports, checks every sum branch, and retains necessary multiplicity checks. Unchanged branches
are borrowed; changed exact sums and products use bulk construction when callback ordering
permits it.

Closed gamma traces emit factored local identities into the existing literal alias
registry. Four-dimensional short-word and symbolic-dimensional recurrence rules
share that result owner. The interface proof checks dimensions, surviving ports,
unresolved identities and normalization behavior. Callback-sensitive output is
validated rather than inheriting the original interface blindly. Neither trace
simplification nor its surrounding scheduler expands the full result.

Uniformly reversed closed gamma and gamma5 words use the trace transpose
identity before entering those same kernels: reverse the factor order and
restore forward endpoints. The normalization remains local to the trace, so
surrounding callbacks observe the evaluated result in the same pass. Mixed
orientations and open chains remain opaque, and gamma5 retains its existing
four-dimensional restrictions.

Lorentz-dimension changes map both the Atom and its declared interface here. Newly coincident
explicit indices contract; excess occurrences fail. Callback validation compares the encoded
interfaces before and after the change, separately from the retained logical layout. A no-op
normalizer must accept mixed representations stored in a different order, while a normalizer
that removes a nonzero tensor's ports must fail. Python only converts the requested dimension,
translates errors, and updates presentation metadata if the rank changes.

`SymbolicNet<Aind>` is a Spenso `Network` whose local tensors are `SymbolicTensor`, whose scalars
are Symbolica atoms, and whose function keys are Symbolica symbols. `SymbolicNetParse` forces
Spenso's `ContainsReps` tensor filter, so representation-bearing functions can become tensor
leaves even without a separate tensor tag. `SymbolicNetExt::simple_execute` uses sequential
Spenso execution and reconstructs the final atom from zero, one, or the result tensor.

The owned network remains the boundary for concrete component execution and
structural consumers which need owned tensors. Symbolic contraction instead
reads Spenso's borrowed operation graph. Public tensor transformations return
the shared typed tensor or its aliased specialization; raw `Atom` utilities are
explicit lower-level boundaries.

== Rewrite families and execution flow

The shared tensor scheduler runs the selected identity families. Each family owns its
normal form:

- `IndexTooling` parses enough structure to canonicalize tensor indices, list dangling slots,
  wrap all indices or only dummies, form adjoints, and alias tensor subexpressions;
- Typed `collect` and `coefficient_list` use the shared graph/tape owner to retain selected
  representation sectors while keeping unrelated coefficients factored;
- the shared component contractor and internal chain/dot normalization own metric
  contraction, compact scalar products, open chains and closed traces;
- the Dirac pass collects bispinor chains and applies dimension-gated Clifford, trace,
  projector, gamma5, and optional four-dimensional epsilon rules;
- `ColorSimplifier` collects fundamental color lines, closes traces, applies generator,
  structure-constant, Fierz, and Casimir rules, prunes antisymmetric zero terms, and iterates to a
  fixed point;
- `EpsilonSimplifier` owns epsilon contractions and reductions;
- `Cookable` replaces selected functions or representation-index payloads with compact symbols,
  either as readable flattened names or reversible Symbolica `UserData::Atom` encodings.

Settings objects are part of the semantics. Gamma ordering and trace evaluation, color Fierz and
invariant substitutions, contraction budgets, and cooking source/tag filters can all
change the result form. Reproducible callers should record the exact settings and the order in
which independent rewrite families ran.

Most pattern engines use a local fixed-point loop: transform the current atom, compare it with
the previous atom, and stop when unchanged. That guarantees termination only for the implemented
oriented rule set; custom rules composed by a caller remain the caller's responsibility.

== One symbolic contraction path

The tensor contractor reads the shared structural graph once, keeping opaque
leaves and their logical ports. It selects an incident factor, visits its
alternatives and applies port substitutions through the existing component
reducer. Equal remaining states merge before result emission. Metrics alone
and metrics with rank-one tensors use this same owner. Powers whose bases
contain local dummy pairs finish that scalar scope before the power is applied.

Callback-sensitive syntax uses the existing ordered substitution schedule and
checked result construction. It does not acquire a trusted interface merely
because a substitution was algebraically a contraction. The example
`g(a,b)*T(a)` whose normalizer turns `T(b)` into a scalar remains a rejection.

Ordered substitutions carry committed alias-port changes into the existing
literal-use registry. Discovery probes and refused linear branches discard
those observations; a product that vanishes publishes none. The alias owner
registers the rewritten definition with its checked surviving interface before
publishing the result, including when the ordered fallback exposes more work
for the component reducer. Literal relabellings use the pass's original
uncontracted templates. Completed bodies are staged separately and published
together, so traversal order cannot compare definitions at different rewrite
stages. Exact imaginary coefficients remain scalar leaves of the same graph.

After the original alias templates finish, the existing repeated-index scan
decides whether a domain can contract across literal alias boundaries. Disjoint
open definitions stay opaque; their sums are not visited merely because they
contain metrics. Sum alternatives and powers retain the scanner's established
scope rules. The private contraction status distinguishes deferred opaque work
from an exhausted budget. Deferred work can expose its registered definitions
and is checked again by the existing simplification fixed point; an exhausted
frontier is retained exactly and never restarted by this exposure. Neither
state is certified complete merely because its algebraic output is valid.

Epsilon simplification uses the same selected collector with the literal
registry. This exposes epsilon bodies and their incident metric/vector factors,
so antisymmetry also acts across alias boundaries. Unrelated definitions and
callback-sensitive templates remain opaque.

`contract_ports` is the separate logical binary composition operation.
`to_dots` and `undo_dots` change symbolic notation; `to_dots` does not silently
run a full contraction. The intrinsic `g` normalizer exploits dot symmetry at
construction. A finite `expand_dots` request crosses Spenso's component
execution boundary and leaves unrelated factors untouched.

The former network Schoonschip engine, its traversal/order strategy family and
duplicate metric collector have been removed. The shared tensor scheduler runs
one fixed point over the root and reachable definitions, carrying observations
and completion proofs when their invariants allow it. A default rerun of a
certified complete aliased result reuses its sealed allocation. Budgets retain
exact unfinished work and do not certify it complete.

The graph, substitution and callback boundaries are detailed in the
#link("schoonschip-net-parsing.typ")[shared contraction architecture note].

== Parsing, shorthand, and index invariants

Spenso representation syntax is the shared contract. Slots carry a representation, dimension,
abstract index, and dual orientation. Products contract matching dual slots; additions must
expose compatible external structure. `chain`, `trace`, `dot`, metrics, and bracket-like syntax
are parser-owned shorthands rather than arbitrary opaque functions.

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
remain opaque. Alias definitions are handled by the typed alias owner rather
than guessed by the raw chain walker.

Index wrapping is an ownership operation. `list_dangling` discovers external slots through the
shared parser's structure-only construction; `wrap_dummies` changes only non-external index payloads; `wrap_indices` changes explicit index payloads while leaving unresolved
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

The optional `bincode` feature derives binary encoding only for Idenso's zero-sized
representation marker types. `SymbolicTensor`, rewrite settings, symbol registries, networks,
and cooked expressions do not form an Idenso checkpoint format. Reversible cooking stores the
source atom in Symbolica symbol user data and derives a stable-looking hash name, but recovery
still depends on the matching Symbolica registry and cooking settings. Printed or binary names
alone are not a portable physics-result archive.

== Ownership and error boundaries

The public tensor operations retain the logical interface, typed zeros and alias
metadata in `SymbolicTensor`. Spynso converts arguments, dispatches to this
owner, translates errors and wraps the result. Symbolica owns atom memory,
normalization and literal alias definitions; Spenso owns graph topology,
logical layouts and component execution.

Structural entry points are fallible. Inference and certified tensor rewrites
report `TensorInferenceError`; network parsing and execution keep their typed
network errors. `dirac_adjoint` reports `AdjointError`, cooking reports
`CookingError`, canonicalization reports `CanonicalizationError`, and raw index
queries report `IndexToolingError`. A raw-expression escape hatch does not
retain a separately declared logical layout automatically. Callers restoring
that layout must supply the original structure.

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

The saved `examples/notebooks/tensor_benchmark_m1.json` compares three alternating
fresh-interpreter runs against the immutable M0 release on CPU 8. Each run has
one unmeasured warmup. The primary clock includes rule application, contraction
and explicit scalar materialization; FORM uses batched body CPU time, excluding
process startup. These measurements are pinned among cooperating jobs, not on an
exclusive host.

#table(columns: 4,
  [Case], [M0 CPU ms], [M1 CPU ms], [FORM CPU ms],
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
