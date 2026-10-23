= Idenso implementation architecture

#quote(block: true)[
#strong[Status:] Current implementation architecture, audited against the Idenso source on
2026-10-18. This date records source review; qualification is recorded separately below.

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

Tensor metadata is authoritative after initial admission. Raw Atom input is checked against
the tensor grammar, declared or inferred interface, placeholder scopes and index multiplicities
at that boundary. Ordinary typed algebra derives its result from admitted operands and carries
their structure and scoped index usage; it does not reconstruct or revalidate the resulting
Atom to recover facts already established by the operation. Symbolica normalization may erase
visible factors without erasing their construction history. Callbacks used by ordinary algebra
must preserve the admitted tensor meaning; arbitrary interface-changing transformations belong
at a checked raw-expression boundary.

`SymbolicTensor<S, E = Atom>` owns the tensor interface and symbolic payload:

- the structure `S`, defaulting to `OrderedStructure<LibraryRep, AbstractIndex>`;
- the symbolic payload `E`, normally an owned Symbolica `Atom`;
- `is_metric` and `is_composite`, retaining operation classification independently of the
  current Atom's outer syntax; these flags participate in equality and serialization,
  and are not interface certificates;
- authoritative scoped explicit-index usage, with entries retaining their original `Slot`
  orientation while usage keys match representation base, dimension and index. An aligned
  dummy mask records multiplicity: one use has an unset bit, two uses have a set bit, and a
  third use is rejected. Two uses commonly form a contracted pair, but can also be two free
  ports with the same orientation; the logical free interface determines their actual role;
- compact per-variance maxima alongside that total-use mask. These retain both orientations
  when a dimension mapping separates a hidden use from a public port. Alternatives combine
  maxima; products add uses. Their sum can exceed the alternative's maximum total;
- private admission, observation and completion facts. Construction metadata travels with
  certified operations; payload-specific graphs, observations and completion facts are retained
  only when the operation justifies carrying them to its result.

The usage metadata is conservative within its multiplicative scope. A product that normalizes
to `1` or `0` can still retain a contracted pair, so a later factor using the same index cannot
silently become a third occurrence. Raw admission of the printed `1` has no such history and is
a different construction. Addition exports its shared external interface while each summand
binds its internal pairs; a closed scalar power likewise owns its internal scope. Logical
structure separately retains orientation, axis order and the declared shape of zero.

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

FeynKit retains this carrier through instantiated graph fragments, diagram
numerators, projectors, amplitudes and tensor-reduction results. GammaLoop uses
`SymbolicTensor<IndexedStructure>` for runtime graph fragments, numerator states,
UV sectors and shared numerator definitions. `IndexedStructure` is
`Canonicalized<OrderedStructure<LibraryRep, AbstractIndex>>`: it retains logical
axis order while requiring explicit public indices. Bare `OrderedStructure`
alone would discard the outer layout needed by ranked zeros and component axes.
Model recipes and genuinely scalar coefficients remain
Atoms. Borrowing `expression()` for
a scalar algorithm leaves the source interface available; rebuilding an unrelated
tensor from the returned Atom would lose that interface, conservative index history
and reduction facts.
Mathematical equality alone cannot decide whether an old shared carrier may be
reused: `same_reduction_state` also compares its validation and completion state.

`try_into_indexed` checks an unvalidated carrier before converting its interface,
rejects unresolved public ports (including scoped open indices), and transfers
its payload, operation flags and proof handles without allocating indices or
re-sorting canonical slots. Closed scalar shorthand remains valid. New source
values use `infer_indexed`; declared indexed interfaces, including ranked zeros,
use `from_indexed_parts`. Both constructors retain the shared validation rules.
Strict admission also follows reversible cooked index payloads through the
existing cooking owner: an encoded Open or scoped Open cannot become an explicit
port merely by being stored as a Symbol. Only public index metadata is inspected;
opaque named labels and closed compact contractions retain their meaning.
Use receiver syntax for shared methods such as `tensor.multiply(...)`.
Unqualified associated-call syntax can be ambiguous between the indexed and
partial specializations; static operations require the explicit structure type.

Indexed algebra methods reuse the existing Idenso engine. Consuming conversion
moves the payload, layout and proofs. A borrowed conversion to an owned partial
carrier copies the Atom's packed bytes as well as rank-sized metadata; Symbolica
Atoms are not shared handles. Checked-result publication instead borrows the
original Atom and proof cells through a private `BorrowedTensorParts` view, so
new validation facts update the source and the supplied result moves into its
carrier. Dummy reservation likewise borrows the payload and preserves logical
index order, including the order observable by named-index callbacks.
The generic payload parameter alone does not provide this borrowed engine API:
the existing carrier still owns its proof cells. The private view borrows both,
so observing a source does not update a detached copy of its validation cache.

The reduction planner accepts borrowed or owned input through one `Cow` core.
Validation and wildcard-label handling precede ownership acquisition. An indexed
adapter's owned projection moves into that core; a borrowed partial input is
copied only when a completed result or mutable plan needs ownership. This removes
the extra planner-entry payload copy without changing public method signatures
or identity rules. Completion-policy metadata is retained separately.
Indexed `multiply(self, right: Self)` consumes both operands and moves them into
the existing partial multiplication owner. This avoids copying an already owned
graph-fold accumulator and preserves the same product planning. Compatibility and
multiplicity checks use the operands' retained metadata before Atom normalization;
the product receives the derived structure, sorted usage and dummy mask without
result-payload admission. Callers clone only inputs whose owning graph, shared definition or
Cartesian product needs them again. The partial-interface multiplication method
retains its borrowed API. New products receive new result facts; input reduction
proofs are not certificates for their product. Collection, dimension/index
transforms, notation changes and adjunction still use an owned projection and
therefore copy their inputs. No conversion exports,
decodes or re-infers an admitted payload. Actual workload costs require separate
measurement; explicit indexing is an interface invariant, not a speed guarantee.
A result with unresolved public ports is rejected instead of materialized.
Native FeynKit export uses the standard indexed-to-partial `From` conversion to
transfer the same tensor and proofs into its more general interface. Coefficient
listing remains a partial-engine API because selected subexpressions can have
unresolved ports even when the complete input is scalar; independent diagnostic
oracles may use it without changing GammaLoop's runtime storage contract.

External declared interfaces enter through `from_declared_parts`, which checks
the payload against the logical slots once, rejecting a nonzero payload with an
unrelated declared interface. Bincode instead transports an admitted tensor:
it stores the payload, logical slots, operation flags, sorted explicit-index
usage, dummy mask and optional product-readiness fact. Decode remaps symbols and
checks metadata consistency, then restores construction history, including pairs
absent from the normalized payload. Future products reuse established readiness
instead of classifying the payload again. Planner graphs, observations and
reduction-completion caches are not serialized. Ranked zero and unresolved
logical axis order survive both boundaries.

Opening dot or chain notation and canonicalizing dummy names establishes new index identities.
Those materialization owners record the new uses and retain conservative source history; they
do not re-infer the free interface or reject a valid typed product because its normalized Atom
has coalesced into a power. Ordinary scalar algebra carries the established history directly.

The shared `feynkit_graph::NumeratorInterface` codec owns source admission and
export. Diagram DTOs and DOT use ordinary expression strings with optional
`num_interface`, `numerator_prefactor_interface` and `projector_interface`
metadata. The declaration preserves axes that cannot be recovered from source
text alone, especially ranked zero and unresolved logical order. Admission
validates the declaration; it is not a way to certify an arbitrary payload.
Ordinary source text remains unchanged when no metadata is needed. Whole-payload
decoding belongs at an explicit source or backend boundary. Graph transport
instead decodes individual encoded index names, computes their simultaneous
rename and uses the native checked tensor operation.

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

Reduction and collection return the same symbolic tensor with a plain Atom
payload and ordered partial interface. Generated sums remain factorized inside
that expression; no alias registry or separate result wrapper is created.
Python returns `TensorExpression` for both scalar and open-tensor results.
`to_expression()` exposes the ordinary Symbolica Atom and `expand()` explicitly
distributes it. The optional aliased implementation and its measurements are
preserved separately on the `tensor-aliases-preserved` jj bookmark.

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
share a coefficient. Scalar
spectators, disconnected components and foreign sectors remain factorized.
The current state key retains original factor identity and port bindings, so it
can under-merge different bindings that normalize to equal factors. This is a
performance limitation, not an algebraic equivalence assumption.

Contraction, rule application and polynomial materialization are separate
operations. `TensorExpression.expand()` explicitly materializes arithmetic
sums and products while preserving the checked external tensor interface.
User callbacks retain their checked boundaries. Kernel observations avoid
re-inferring known interfaces solely because a trusted transformation returned.

Frontier growth is bounded before processing the next factor. The current guard
counts generated alternatives and estimates state bytes; it is not a
strict heap cap on coefficient buffers. If the guard fires, the result retains
exact pending factors, and its internal completion flag remains false. It must
not be certified as fully contracted. Callback-sensitive inputs continue through
the existing checked contractor, including the rank-loss rejection for a
normalizer that turns `T(b)` in `g(a,b)*T(a)` into a scalar.

Source-attached algebra preserves established interfaces and construction metadata when its
operation owns the invariant. Attaching unrelated raw parts does not establish that invariant:
inference can retain raw repeated ports until `checked_parts` merges their explicit
contractions. This intermediate boundary also lets bindings discard stored-data identity when
a declaration becomes a contraction, even if its input and final result have the same rank.
Arbitrary expression replacement and callback-sensitive rewrites which cannot certify their
interface check the actual output, because normalization can remove tensor ports. That boundary
observes existing ports without temporary index materialization or callback replay;
constructor inference retains its materialization semantics within the same inference owner.
Index substitution distinguishes graph storage identities from logical unresolved
ports, checks every sum branch, and retains necessary multiplicity checks. Unchanged branches
are borrowed; changed exact sums and products use bulk construction when callback ordering
permits it.

Fully explicit external ports can be compared by their complete slot identities
when normalized factor order changes. This retains the source logical layout
without treating syntax order as axis order. Representation, dimension,
variance, index identity and multiplicity must all agree. Anonymous axes retain
their stricter ordered correspondence: equal rank or equal representations alone
do not justify swapping them. The checked-expression fallback collects tensor
witnesses while leaving unstructured scalar weights opaque, just as ordinary
product inference does. Its final scope and multiplicity validation still runs.
The special affine rank-one rewrite requires a uniquely identified changed
factor with the remaining tensor witnesses fixed.

Before power and interface admission, the existing dot normalizer lowers
nested rank-one notation such as `δ(cind(0),p(mink(d)))` to
`g(δ(cind(0),mink(d)),p(mink(d)))`, retaining the symbolic dimension.
The same contraction accepts scalar-weighted vector sums when the existing
materializer certifies one compatible compact axis and a scalar coefficient.
It retains the summed operand and its scalar spectators without distributing them.
This admission pass leaves explicit metrics, dots, indexed contractions and
scalar payloads unchanged. Callback-sensitive ancestors and metadata remain
on the existing checked path. The native fixed-component delta shares the
intrinsic-normalizer classification used by metric, dot and dual notation:
its normalizer only compares literal components. Descendants are still checked
for callbacks; this classification is not an interface certificate. Full admission
validates the normalized expression before publishing it.

Closed gamma traces emit factored local identities into the tensor expression.
Four-dimensional short-word and symbolic-dimensional recurrence rules
share that result owner. The interface proof checks dimensions, surviving ports,
unresolved identities and normalization behavior. Callback-sensitive output is
validated rather than inheriting the original interface blindly. Neither trace
simplification nor its surrounding scheduler expands the full result.

Exact scalar weights in Dirac words are extracted without distributing their
products. A weighted sum of ordinary slashed vectors enters the same trace
recurrence as one vector argument; it does not generate a Cartesian product of
matrix words. Explicit-index matrix sums remain outside that rule. Weighted
compact metrics and dots share complete-interface admission before bilinear
normalization. Uncontracted mixed physical `gamma0` traces require explicit
four-dimensional spin and Lorentz representations; formal dimensions and
disabled trace evaluation retain their words. A repeated abstract Lorentz
index around `gamma0` has the separate dimension-generic Clifford contraction
`gamma_mu gamma0 gamma^mu = (2-d) gamma0` in physical four-component spinor
space. Open words and enabled traces share this pair-contraction owner. Fixed
components, compact slashes, distinct or incompatible ports, and formal spin
dimensions do not enter that summed-index rule.

Uniformly reversed closed gamma and gamma5 words use the trace transpose
identity before entering those same kernels: reverse the factor order and
restore forward endpoints. The normalization remains local to the trace, so
surrounding callbacks observe the evaluated result in the same pass. Mixed
orientations and open chains remain opaque to that trace-transpose identity,
and gamma5 retains its existing four-dimensional restrictions.
An adjacent pair of ordinary gamma matrices with identical vector arguments
may also have both endpoint orders reversed: its product is the transpose of
the ordinary equal-argument product and has the same metric coefficient.
This local contraction preserves the original spin representation, including
formal spin dimensions, and never treats a mixed-orientation pair as equal.
The same local owner handles equal abstract-index gammas on either side of
physical `gamma0` in a four-dimensional spin representation. Since `gamma0`
is symmetric, transposing the contracted sandwich preserves `(2-D)*gamma0`;
an extra neighboring `gamma0` gives `(2-D)` times the spin metric. The summed
Lorentz dimension remains formal. Fixed/vector arguments, mixed orientations,
incompatible representations and generic spin dimensions keep their existing
rules, and the global parser still leaves reversed open words opaque.

Fixed components of the physical Minkowski basis have a separate local rule.
Their metric contractions use the physical signature while retaining the
symbolic Lorentz dimension. In a four-dimensional spinor representation, a
slash of the temporal basis vector becomes `gamma0`, and `gamma0` anticommutes
with each fixed physical spatial slash. These identities preserve ordered
spinor multiplication and do not set the Lorentz dimension to four. Formal
spin dimensions and unsupported component representations retain their word.

Analytic UV admissibility is separate from general index-pair observation.
A gamma matrix projected onto a validated fixed physical basis vector does not
require an unresolved D-dimensional tensor reduction. The analytic boundary
excludes only recognized physical Minkowski selector leaves (coordinates zero
through three), using borrowed component and representation bounds; it leaves
scalar D and epsilon dependence untouched. Its per-leaf predicate distinguishes
valid selectors from unchecked uses of the same delta head. The query first
requires the original unfiltered admitted index metadata; excluding a leaf cannot
hide a malformed source. Filtered discovery uses a separate fresh cache. General pair
observation still counts these connections when no exclusion is requested, and
unknown/out-of-range selectors and genuine summed Lorentz contractions remain
visible. The borrowed bounds query neither invokes normalization callbacks nor
materializes indices; the compact-vector recognizer keeps its existing grammar.

`map_representation_dimensions` maps representation dimensions in the payload,
logical interface and retained scoped index usage, including ports carried only
by a ranked zero and indices absent after cancellation. Its callback receives
the representation kind and the decoded dimension; portable dimension payloads
are decoded and encoded individually. It does not substitute unrelated scalar
coefficients. `with_lorentz_dimension` delegates to this owner for its
physical-four-dimensional Lorentz replacement. Newly coincident compatible
explicit ports contract; invalid dimensions, incompatible internal pairs and
excess occurrences fail. The operation derives its result metadata from the
mapped source metadata without admitting the rebuilt Atom again. Normalizers
triggered while rebuilding must preserve that mapped tensor meaning. This
operation does not provide component tables for a new formal dimension.

`rename_symbols` applies an injective simultaneous map to the native payload,
logical interface and retained scoped index usage, including encoded dimensions
and labels absent from the simplified payload. Swaps are allowed; collisions
with existing targets that are not themselves renamed are rejected. Mapped
metadata establishes the result's interface and index roles; normalizers
triggered while rebuilding must preserve that mapped tensor meaning.

The native `dirac_adjoint(preserve_indices)` shares preparation and cleanup with
the raw expression operation while retaining the tensor carrier. It dualizes
the external slots and preserves their logical order. Anonymous axes first
become explicit fresh indices through `materialized()`, including for ranked
zero, so their identity survives adjunction. The option retains the historical
boundary-endpoint transpose convention; it does not relax interface validation.
Spinor words receive the required gamma-zero construction and cleanup. Pure
scalar and color expressions retain the existing conjugation semantics. The
typed operation uses `AbstractIndex` directly without a whole-payload
decode-and-reinfer round trip.

`SymbolicNet<Aind>` is a Spenso `Network` whose local tensors are `SymbolicTensor`, whose scalars
are Symbolica atoms, and whose function keys are Symbolica symbols. `SymbolicNetParse` forces
Spenso's `ContainsReps` tensor filter, so representation-bearing functions can become tensor
leaves even without a separate tensor tag. `SymbolicNetExt::simple_execute` uses sequential
Spenso execution and reconstructs the final atom from zero, one, or the result tensor.
Native symbolic consumers can select `AbstractIndex` directly. A legacy index
type is not needed merely to materialize compact notation for a structural
check; decoding and recooking can otherwise introduce collisions between
generated dummies and encoded source names.

The owned network remains the boundary for concrete component execution and
structural consumers which need owned tensors. Symbolic contraction instead
reads Spenso's borrowed operation graph. Public tensor transformations return
the shared typed tensor; raw `Atom` utilities are
explicit lower-level boundaries.

The consuming `SymbolicTensor<PartialStructure>::into_network` specialization
hands off to the existing Spenso component parser. It validates the retained
interface, materializes unresolved ports, and releases its own planning-cache
references before parsing components. Other shared owners can still retain the
same caches. The optional configuration must select reversible index capture;
its `dimensions=false` source filter controls cooking, not inverse restoration.
`uncook` also restores portable dimension payloads, including dimensions that
have become scalar coefficients. Actual representation dimensions must be
specialized through the checked dimension owner before finite component
execution. Restoring a symbolic dimension does not provide a component table.
This handoff still parses and materializes component tensors. The result includes
`CanonicalLayout`, which downstream code must use
for component coordinates and shadow axes; matching only the port multiset does
not preserve logical axis order. A ranked zero becomes an empty sparse tensor
with its original ports. Concrete components and contraction execution remain
Spenso responsibilities.

The indexed `into_network` entry point consumes the same retained carrier and
uses this handoff with an already explicit public interface. Its logical layout,
ranked-zero handling and index-decoding checks are unchanged.

The Atom-like API distinguishes raw replacement from source-attached operations.
Symbolica seals `AtomCore`, preventing a downstream implementation on this
carrier. `with_checked_expression` checks an arbitrary replacement against the
retained interface; `with_transformed_expression` infers an arbitrary transformed
interface. Both are admission boundaries, including for scalar outputs.
`with_rewritten_expression` retains its checked path for callback-sensitive
identity results. Component realization also checks its retained interface at
the network boundary.
`replace` accepts a `TensorRule` whose actual matched interfaces and index bounds
are checked locally. Addition, multiplication and powers use fallible tensor
operations: they check operand metadata and derive result metadata before payload
normalization. `with_algebra_result` and certified identity finishers instead
require the producer's source-attached algebraic invariant. They carry the admitted
structure and conservative index history when that invariant applies; they do not
certify unrelated replacement Atoms. Producer-certified simultaneous index
relabeling maps retained usage as well as public ports, so injectivity also covers
identities absent from the simplified payload.
There is no typed Taylor-series or general derivative API yet, so GammaLoop
restores scalar-engine results against the retained interface. `expand_dots`
returns the tensor carrier, preserving untouched open ports while checking
transformations that could change them.

An identity dimension mapping retains completed reductions and cached
observations through `finish_parts`, after running the mapper and establishing
that the native payload, logical layout and oriented index-usage metadata are
unchanged. Admitted carriers reuse their facts; raw carriers still undergo
admission. An unchanged zero expression or free interface alone cannot establish
that hidden index history is unchanged.
The color-construction macros have a separate syntax trap: a bare `object.field`
is representation-pattern shorthand. Literal chain placeholder symbols must use
`Atom::var`, rather than that field shorthand.

The typed coefficient collector is not a drop-in replacement for every scalar
parser. GammaLoop's rational UV shell must keep positive-power numerator blocks
such as `(den + x)^5` factorized while extracting inverse denominator factors.
Its existing scalar parser still uses temporary Atoms and intermediate
inference; `FourDTerm::from_tensor` checks each completed numerator against the
original carrier to restore its logical interface. That residual conversion is
an API limitation, not a claim that all tensor work now avoids parsing.

`CookSettings::dimension_payload` exposes the existing portable dimension
decoder to consumers of retained tensors. It returns `None` for ordinary
symbols and index encodings. A reducer can compare the decoded scalar dimension
while retaining the original native slots; a scalar-ratio validator can expose
`Nc^2 - 1` after all tensor factors cancel. Neither use requires decoding the
complete tensor. Mapping represented dimensions and interpreting an encoded
dimension that has become a scalar coefficient are distinct boundaries.

Exact compact vector selectors have their own authored-source boundary. FeynKit's
`with_integrated_vector` and `with_external_vector` admit each small selector once
through the shared source codec and return a fallible result. The retained
numerator goes directly to `reduce`; matching must not compare a raw compound
dimension against its encoded native spelling. Graph reduction promotes source
Lorentz dimensions through `map_representation_dimensions`, retaining ordered
ports and ranked zero, while projector coefficients use the supplied scalar
dimension. This selector admission is distinct from numerator re-inference and
from a finite component-execution boundary.

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

Formal SU(N) color reduction recognizes complete invariant fundamental networks
with zero or one external adjoint. Singlet/adjoint projection and the local
cubic symmetric-trace contractions retain other open lines and scalar factors.
Known dimension relations are decoded from portable representation payloads
only for identity selection; the original port payloads remain intact. An
adjoint generator contracted with its incident chain endpoint vanishes by
structure-constant antisymmetry before network materialization. Trace evaluation
and cross-line Fierz expansion are independent settings.
Adjacent repeated raw adjoint matrices also reduce in their retained word:
`f(i,k,a) f(k,j,a) = -C_A delta(i,j)`. The existing orientation reader supplies
the matrix phases; explicit structure-constant products independently check the
sign. This closes internal color labels before scalar coefficients reach an
analytic backend.

Two invariant symmetric rank-four color tensors sharing two distinct adjoint
ports can close their remaining ports through two structure constants or a raw
two-matrix adjoint word. Invariance gives `C_A/3` times their scalar Gram
contraction; the raw word has the opposite sign. Representation phases,
structure-constant orientation, matching dimensions and complete closed wiring
remain explicit guards. The bounded matcher retains spectator factors rather
than expanding a complete graph numerator. Scalar metadata carrying
representation arguments is excluded from both matrix-chainification paths;
genuine tensor children still participate through their existing owner.

A retained adjoint chain can also connect to a separate fundamental trace or
color line. The product owner exposes that one chain through Spenso's reserved
dummy allocator only when an existing cross-factor identity can consume the
connection. A notation-only opening is discarded, preventing an open/recollect
cycle. This lets the same structure-constant identities act on explicit products
and compact chains without distributing spectator factors or losing the logical
ports of a zero result.

A two-matrix adjoint word can also close a loop through its opposite endpoint:
`chain(A,B; F(B),F(C))` represents `f(A,u,B) f(u,B,C) = C_A g(A,C)`.
The chain owner exposes only that local word and reuses the existing two-f
contraction, including its matrix phases. Shared structure-constant positions
represent summed edges only when their indices are abstract; equal fixed
components do not authorize a sum. Fixed uncontracted spectators are preserved
by the local identity. Public symbolic tensor admission currently rejects the
raw fixed-color chain and structure-constant forms used in these controls, so
their fixed-index guard is tested directly at the identity owner.
Identical materialized factors can normalize to a square. The ordinary and
squared forms share this contraction owner, retaining scalar weights and the
same summed-index guards without distributing sums.

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
Unbound chain endpoints do not yet have named indices, so `list_dangling`
rejects them. Scalar metadata stays opaque to that check.

Index canonicalization reconstructs the admitted operation graph without
executing tensor algebra. Closed scalar inverse powers keep their own dummy
scope while their interior indices are canonicalized. These preparation and
conjugation operations remain orthogonal to the reduction planner.
Dual-index wrappers remain in the canonical graph. The base index links both
occurrences, while `dind(...)` retains the variance of each occurrence. Erasing
variance and restoring it by function name can turn oppositely oriented metrics
into an invalid tensor square. A closed fundamental metric cycle stays a scalar
with trace equal to its dimension, independent of the dummy labels.

The shared concrete syntax and which crate owns each rewrite are specified in the
#link("spenso-symbolica-syntax-and-rewrites.typ")[Spenso/Idenso Symbolica syntax note].

== Notebook presentation

Spenso's existing `network::tags::mathml` printer emits bounded semantic leaves
directly as MathML, reusing the tensor, slot, representation and chain
recognizers. Gamma and slash notation, ordered chains and traces, bra/ket ports,
color invariants and display settings use native MathML. Idenso retains ownership
of color-invariant semantics. Presentation preserves the source Atom bytes and
the typed logical interface; it does not rewrite tensor algebra.

Spynso's viewer adapter retains typed metadata, display-local index aliases and
a bounded cache of leaf fragments. Alias resolution is scoped to the viewer and
supplied to the printer for each bounded leaf; it requires no complete-source
walk. Symbolica's official `MathMLViewer::expression` owns arithmetic, cursor
caching, responsive layouts, paging and navigation. The pager has no CBOR render
tree, Typst compilation or notation worker.

Explicit expression `to_html()` emits complete direct MathML; it does not page,
truncate or compile Typst. `formatted()`, `paged()` and `to_html()` reject
`notation_source` with `ValueError`. `TensorName` head mappings accept `mathml`
fragments, and callable printers handle `PrintMode.MathML` for both rich display
and complete HTML exports. Typst source is never parsed into MathML.

Bare `to_typst()` remains source-only.
`to_typst_project(..., notation_source=...)` returns a `dict[str, bytes]`
containing `main.typ`, `render.typ`, `notation.typ` and `tree.cbor`, without
invoking a compiler.
Expression `to_svg()` compiles that project through optional Python typst-py
(distribution `typst`); the native extension contains no Rust Typst compiler.
Custom Typst notation belongs to these project and vector exports.

Tensor-network `render()` returns a shared Linnet `DiagramView`. Its live
AnyWidget measures direct MathML labels in the notebook browser before the
existing native `Scene::render` draws the geometry. `view.measured()` retrieves
the completed immutable `DiagramRender`; SVG and HTML retain the measured boxes
and embedded font. Scripts explicitly request `measured(browser="chromium")`
with optional Playwright and Chromium, or return a receipt from the same view's
`measurement_request()` to `measured(measurements=...)`. Pending measurement
raises an error on export; there are no estimated bounds or automatic headless
browser. Default graph SVG uses XHTML `foreignObject` MathML and targets browsers.
For PDF, PNG or Typst SVG import, `view.outlined()` uses optional typst-py to
produce vector labels.

== Features and serialization

Idenso defaults to `native`, forwarding GMP/MPFR support; `wasm` selects the Wasm backend with
default features disabled. The core Rust rewrite layer always includes Symbolica and Spenso's
`shadowing` support. The `spynso3` crate provides the Python methods and automatic representation
initialization; its `python_stubgen` feature adds metadata for the unified Spenso module.
`reference-cases` exposes the otherwise test-only curated identity cases.

Binary encoding is always derived for Idenso's zero-sized representation marker types.
The optional `bincode` implementations for partial and indexed symbolic tensors
encode admitted tensor metadata together with the payload, operation flags and
optional product-readiness fact. Restored readiness lets future products reuse
that fact without classifying the decoded payload again. Decode uses the caller's
state map, restores sorted index usage after symbol remapping, and checks archive
version, unique index keys, mask alignment and
agreement between retained index roles and the logical interface. Indexed
archives additionally reject unresolved public ports. This is a transfer of
authoritative construction metadata, not raw-expression admission: cancelled
pairs survive cloning and archive round trips even when the payload is `0` or
`1`. Disposable planner, observation and reduction-completion caches are not
serialized. The tensor archive carries its own format version; GammaLoop's
already lowered standalone evaluator formats are unchanged.
These wire implementations, rewrite settings, symbol
registries, networks and cooked expressions do not together define an Idenso
checkpoint format. Reversible cooking stores the
source atom in Symbolica symbol user data and derives a stable-looking hash name, but recovery
still depends on the matching Symbolica registry and cooking settings. Printed or binary names
alone are not a portable physics-result archive.

Explicit representation-dimension cooking is a narrower exception: reversible
mode encodes supported exact rational arithmetic in ordinary symbols into a
portable name. Atomic dimensions remain unchanged, and decoding checks the
payload before normalization. Function calls, callback-bearing expressions and
unsupported symbol metadata are rejected. Index/function cooking keeps its
existing registry-dependent behavior; construction and cooking remain separate
from reduction.

Index-only cooking preserves validated native named, Open and scoped index forms
in both flattened and reversible modes. It interns opaque function-valued labels
with expression-local reuse, retaining their shared identity across representation
kinds. Tagged third-party heads that do not parse as native indices are interned
as opaque labels. Malformed Spenso-owned Open and scope forms remain errors.
The separate all-function cooking mode retains its own source/tag policy.

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

Admission establishes the logical interface and scoped index metadata once.
Typed algebra uses those facts as the tensor's meaning, even when the resulting
Atom no longer displays the same index occurrences. Raw admission still rejects
malformed signatures, incompatible sum interfaces, invalid placeholder scopes,
excess index multiplicities and mismatched declared ports; cancellation in a
generated payload is not permission to weaken those checks. A new arbitrary Atom
replacement must establish its own admitted facts.

Gamma identities preserve the source-attached invariant instead of inferring and
validating their generated results again. Generated homogeneous trace sums read the
interface of one representative term through Spenso's existing syntactic
structure reader. This relies on the trace recipe's algebraic invariant;
arbitrary user rewrites cannot claim it. Callback-sensitive results retain
their checked boundary, and unresolved positional ports retain their identity
checks.

Unchanged regions carry the same certificates through planner operations.
Admission does not imply contraction completion or absence of user normalizers: those
are separate facts. Raw mutable payload or interface access invalidates construction
metadata and derived facts; certified algebra carries construction metadata while
discarding caches that describe an earlier payload.
Graph scratch interfaces containing encoded open-port identities are converted
to public logical ports before acquiring a validation certificate.

The historical gamma and alias-domain measurement, including its exact-output
qualification protocol, is preserved with the optional alias implementation on
`tensor-aliases-preserved`.

An earlier full Idenso qualification passed 1,028 tests across five binaries,
with 24 existing skips. It predates the native adjoint, public dimension and
rename operations, anonymous scalar-witness correction and frontend retention
changes described here. It therefore does not qualify those additions. Current
source identities, commands and subsequent results belong to the
#link("feynkit-merge-readiness.typ")[FeynKit merge-readiness record]; the source
review date above does not imply a completed test run.

Domain simplifiers usually leave unmatched syntax unchanged rather than diagnosing it as an
error. A successful return therefore means the configured rewrite reached its fixed point, not
that every physics object was recognized or eliminated. Verification should inspect residual
representation-bearing factors when completeness matters.

== Maintained invariants

- Tensor structure is inferred from tagged representation arguments, including dimension and
  dual orientation; function spelling alone is insufficient.
- Initial admission establishes agreement between the payload and its declared or inferred
  interface. Thereafter structure, sorted scoped index usage and the dummy mask are authoritative;
  ordinary typed operations derive them from operands without re-admitting the result Atom.
  Cancellation can remove visible index occurrences while their collision history remains.
  Ranked zero retains its declared interface. Permutations preserve logical axis correspondence.
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
