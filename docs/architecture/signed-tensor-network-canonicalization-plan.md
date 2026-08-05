# Signed Tensor-Network Canonicalization Plan

## Objective

Replace the current two-stage canonicalization path

```text
signed zero-detection graph -> bare Symbolica tensor canonicalization
```

with one network-owned canonicalization path:

```text
normalized Atom
  -> canonical-policy validation and parse
  -> CanonicalPolicyNet { network, normalized source Atom }
  -> [unified expression-and-incidence graph
      -> one whole-network Graphica canonical form, vertex map, and generators
      -> group-aware reconstruction with signed leaf values and exact-scope Neg
      -> ordinary symbolic network execution
      -> canonical Atom and normalized network reparse]
     skip another Graphica call only when reconstruction certifies stability
     or mandatory execution returns the exact iteration input;
     otherwise repeat on the complete reparsed network
  -> stable CanonicalPolicyNet
  -> canonical SymbolicNet and Atom output projections
```

The unified graph contains both the network operation tree and the flat tensor
incidence data. Products, sums, powers, functions, tensor ports, and index lines
therefore remain visible to the same global labeling problem. The graph is not
contracted into opaque intermediate nodes and Graphica is not invoked
recursively on subgraphs. If reconstruction cannot certify that zero
substitution, signed placement, canonical transport, and normalized execution
preserved the visible problem, the complete reparsed network is projected and
globally canonicalized again; there is still no termwise or local Graphica
pass.

Graphica canonicalizes the unsigned structural graph. Its selected
input-to-canonical `vertex_map` does not by itself determine a signed tensor
representative: equally valid maps can differ by an odd stabilizer. Mapping the
canonical graph back to a network is therefore a group-aware signed
reconstruction. It fills canonical slots into tensor layouts, carries induced
parity only while interpreting the returned automorphism generators, and then
resolves that parity into ordinary network values. A leaf-owned phase negates a
fresh scalar `Atom` or `SymbolicTensor.expression`; a nested partial-group phase
edits its exact Atom path. A whole-value sign in a Product may be sunk to one
deterministic fresh leaf in that same multiplicative scope, preferring an
existing scalar Atom and otherwise a tensor expression. For a phase owned by a
larger operation result, create or reuse the unary `NetworkOp::Neg` variant at
that exact scope when no such sign-transparent leaf exists, because it cannot in
general be pushed through a Sum, nonlinear Function, or Power. Bottom-up
execution evaluates that scope and then negates the scalar or tensor leaf it
produces. There is no operation-node sign/coefficient payload and no new signed
network value type.

Existing scalar expressions remain ordinary visible network leaves. After the
whole graph has been canonized and rebuilt, execute a clone through the ordinary
symbolic network execution path. Symbolic contraction multiplies
`SymbolicTensor.expression` Atoms, and the existing Sum, Function, Neg, and
Power operations likewise produce the normalized Symbolica Atom. Execution may
collapse subgraphs into opaque leaves at this point because global labeling is
already complete. This must produce both a deterministic canonical network and
zeros caused by negative stabilizers.

The graph and reconstruction must encode function symmetries with the same
head-local semantics as Symbolica: an attribute acts on every immediate
argument of the head that owns it. Spenso adds a stricter grammar for
intrinsically symmetric tensor heads so that those immediate arguments are
exactly the tensor slots. Bare `symbolica::AtomCore::canonize_tensors` must no
longer be part of the production path.

## Design Invariants

- A Symbolica symmetry attribute acts on all immediate arguments of its owning
  head. It never skips an argument based on whether Spenso recognizes that
  argument as a slot.
- A tensor-tagged head carrying `Symmetric`, `Antisymmetric`, or
  `Cyclesymmetric` has one direct slot per immediate argument. Reject scalar
  parameters, nested expressions, and `aind(...)` bundles on that head instead
  of giving it a different slot-only interpretation from Symbolica.
- An ordered tensor may mix parameters and slots. Partial slot symmetry is
  written explicitly with a nested `sym(...)`, `antisym(...)`, or
  `cyclic(...)` head, for example `F(p, antisym(a, b))`. When one of these heads
  is used as a structural partial-slot group, each of its immediate children is
  one direct slot; the existing projector use over chain/trace factors remains
  unchanged.
- Ordered argument positions constrain graph isomorphisms. The immediate
  children of symmetric, antisymmetric, and cyclesymmetric heads have their
  positions hidden from graph comparison but retained for reconstruction and
  sign calculation. Nesting never projects the inner symmetry onto the outer
  ordered tensor.
- A cyclesymmetric head has a directed cycle through all of its immediate
  arguments, so rotations are permitted and reflections are not.
- `Linear` is opt-in. The canonicalizer and tensor constructors must not add it
  implicitly. A user-declared linear outer tensor may lift the sign produced by
  a nested `antisym(...)`; a nonlinear outer tensor keeps that sign inside the
  argument. `spenso::bracket(...)` keeps a composite argument such as `p + q`
  atomic with respect to that outer linearity. Rust users opt in through the
  existing `tensor_symbol!("F"; Linear)` attribute forwarding, and Python users
  through `TensorName(..., is_linear=True)`.
- External indices retain their complete `LibrarySlot`; contracted indices lose
  only their dummy name and retain endpoint representation and contraction-group
  information.
- Structured index identity comes from the parsed network and reversible index
  representation. Flattened symbol names must never be used as proof identity.
- Every fixed-point iteration constructs one complete
  expression-and-incidence graph and calls Graphica once on that whole graph.
  Operator nodes remain explicit; no execution or local canonicalization may
  replace a subgraph with an opaque canonical key before the global call.
  The initial input must be the opaque result of the canonical parser, paired
  with the exact normalized Atom from which it was parsed. A raw
  `SymbolicNet`, including one manually assembled or parsed with different
  boundaries, is not a canonicalization input. There is no pre-Graphica
  execute/reparse exception and no API for asserting that an arbitrary network
  is normalized.
- Graph colors contain semantic identity only. Tensor/store indices, hedges,
  original dummy names, and layout member positions are hidden reconstruction
  data and must not influence graph equality or canonical ordering. Visible
  semantic colors must not derive ordering from process-local Symbol interner
  IDs.
- Existing scalar leaves and Product structure remain visible
  semantic graph data. Do not copy Symbolica's term-local optimization that
  strips every index-free factor before graph construction: the unified graph
  needs those values to distinguish global Sum and Function contexts. Signed
  group analysis carries parity transiently after Graphica, but reconstruction
  resolves it into ordinary leaf Atom values, nested Atom paths, or
  `NetworkOp::Neg` at the exact scope before returning a network. `Neg` is a
  supported reconstruction/execution operation, but the canonical parser may
  already absorb an input minus into a scalar or tensor Atom leaf; the plan does
  not require input `Neg` syntax to survive that parse.
- `NetworkLeaf` variants are storage and laziness choices, not mathematical
  distinctions. For a value representable as `LibraryKey`, `LocalTensor`,
  `TensorSum`, `ScaledTensor`, or `ScaledTensorSum`, eager unary execution and
  self-loop tracing must produce the same result through one shared fallible
  materialization boundary. `Power(0)` and identity `Power(1)` are handled
  before that boundary so they do not unnecessarily destroy the original lazy
  or library representation.
- Do not apply that eager boundary to `Neg`, Sum, or Product. Their existing
  lazy tensor sums and scaled terms are intentional, but results must remain
  independent of the chosen leaf representation and of commutative child order.
  Appending one network store to another must likewise rebase every scalar and
  tensor reference exactly once through one shared path.
- The input-to-canonical phase is placed during reconstruction at the symmetry
  owner that produced it. An intrinsic antisymmetric tensor contributes a
  whole-tensor sign; a nested `antisym(...)` contributes one inside its stored
  argument path. Input transport compares the original member frame with a
  canonical local port frame; canonical automorphisms compare canonical source
  and target frames. Hidden original member positions must not define the
  canonical generator parity. Signs combine by parity only at the same semantic
  sink. After that ownership is established, a whole-value sign may move only
  within its sign-transparent Product scope to a deterministic scalar or tensor leaf.
  Reconstruction must not collapse signs across other scopes, mutate an aliased
  store entry, or treat them as general scalar coefficients.
- A negative automorphism zeros the smallest semantic expression scope on which
  its induced action is `-identity`. A zero summand is removed, while a zero
  function argument does not erase a nonlinear enclosing function. An
  identically zero base or complete materialized magnitude scope of a negative
  Power remains a zero denominator and does not make the reciprocal zero;
  ordinary symbolic execution owns the resulting singular value. An
  automorphism may span an outer product and a permutation of sum branches, so
  this scope is derived from the global action rather than from a preselected
  flat subgraph.
- Canonicalization preserves factorization. It must not implicitly distribute
  products of sums merely to discover additional zero monomials. Ordinary
  Symbolica normalization is nevertheless authoritative for forms it rewrites
  without distributive expansion, including an integer Power of a Product; if
  such a rewrite changes the visible network form, reparse it before declaring
  the result stable.
- Materialized Power occurrences have identical visible `PowerCopy` roots but
  retain complete-copy membership. The exact signed exponent is visible on the
  outer Power-result scope, while a typed inner `PowerMagnitude` scope owns the
  unsigned-magnitude copies. This makes a zero denominator scope distinct from
  its reciprocal result, and a negative exponent represents division without
  making `B^n` isomorphic to `B^-n`. Self-dual boundary ports are paired by
  whole copy with two-ended lines according to the existing even/odd magnitude
  interface; no multiway line or per-slot partner choice is allowed. Reject an
  expansion that exceeds the preflight whole-graph resource budget rather than
  hiding it.
- Every reconstructed iteration has one complete-network execution through the
  ordinary symbolic executor, and its normalized Atom is reparsed with the
  canonical parse policy. Exact group proofs may additionally execute temporary
  decorated scope fragments through that same fallible path. Zero removal,
  signed term cancellation, term combination, and any other normalization that
  removes or coarsens a visible
  distinction can expose automorphisms absent from the previous graph. A
  visible execution/reparse rewrite may also establish a different required
  network normal form. Skip another whole-network Graphica iteration only when
  reconstruction produces a complete one-sided stability certificate or the
  executed normalized Atom is exactly the iteration input. Exact equality after
  the mandatory execution is itself a complete extensional fixed-point
  certificate because parsing and projection are deterministic. The one-sided
  certificate proves that canonical transport preserved visible topology,
  semantic colors, scope and Power structure, and the equality partition of
  Symbolica-visible leaves and operation children while changing only
  association, commutative ordering, declared unsigned symmetry order, or a
  bijective dummy assignment. A sign, zero, residual phase, payload edit,
  equality-class merge or split, or any ambiguity declines to certify stability
  and therefore retries after the unconditional execution/reparse. Idempotence
  is a reconstruction-and-execution invariant and a required test.

## Target Ownership and API Flow

The implementation remains network-owned, but canonicalization is not a method
on the raw `SymbolicNet` type alias. After checking the existing parser and
result types for an equivalent proof-carrying type, add the smallest opaque
wrapper required by the migrated callers:

```rust
struct CanonicalPolicyNet<Aind> {
    network: SymbolicNet<Aind>,
    normalized_atom: Atom,
}

CanonicalPolicyNet::parse(normalized_atom, /* parser context */)
    -> Result<CanonicalPolicyNet<Aind>, CanonicalizationError>

canonical_net.canonize(new_dummy)
    -> Result<CanonicalPolicyNet<Aind>, CanonicalizationError>
```

Both fields are private. Construction is possible only by validating a
normalized Atom and parsing it with the canonical policy; there is no
`From<SymbolicNet<_>>`, unchecked constructor, mutable network projection, or
`assume_normalized` escape hatch. Give the wrapper only the visibility required
by migrated callers. A read-only network view and consuming output projections
may be exposed where required, but a projected raw network cannot be fed back
into canonicalization. This makes parser provenance a type invariant instead
of a runtime guess about settings that `SymbolicNet` does not record.

Use one driver that retains this normalized network/Atom pair through the
fixed point. The canonical network and Atom returned to existing facades are
two projections of that one result, not two canonicalization modes, and neither
projection triggers another execution. The first iteration's input-normal-form
precondition is certified by the opaque canonical-parser result; this does not
prejudge its post-reconstruction stability certificate. Every later input is a fresh
`CanonicalPolicyNet` built from the normalized Atom produced by the preceding
post-reconstruction execution. Before adding this wrapper or a traversal
helper, confirm that the existing parser, expression-tree traversal, and result
types cannot supply the same invariant; the current codebase search must be
recorded with the implementation.

Implementation search record (2026-08-05): the existing abstractions were the
raw `Network`/`NetworkStore`, parser settings and transient parser state, and
execution result types. None retained the exact normalized source `Atom`
together with the network parsed from that Atom under one fixed policy. The
private `CanonicalPolicyNet` pair is therefore the minimum wrapper that makes
that provenance an invariant.

The wrapper owns one deterministic parse policy: fully recursive
`ParseSettings::default()` with every Sum term enabled, together with the same
fixed tensor-symmetry validator used at entry. Preserve every parser-supported
signed integer Power. Under the existing grammar, any scalar-valued base may
have a negative exponent, and a self-dual tensor may have a negative even
exponent because its complete-copy contraction is scalar before inversion;
negative odd non-scalar and non-self-dual tensor Powers remain parse errors.
Derive a materialized exponent's magnitude with an unsigned conversion such as
`i8::unsigned_abs()`, never signed `abs()`, so `i8::MIN` means 128 copies and is
supported subject to the ordinary graph-size budget. The parser currently falls
back to an opaque scalar leaf when an exponent is fractional, symbolic, or
outside `i8`; the canonical validator may accept that fallback only when its
base contains no recognized tensor topology. Otherwise return a typed Power
grammar error before Graphica rather than hiding tensor ports. A parser
configured with depth limits, opaque shorthand boundaries, or other
noncanonical settings returns a raw `SymbolicNet` and therefore cannot construct
this wrapper. Ordinary algebraic normalization already represented by the
source Atom is accepted. There is no pre-Graphica whole-network execution or
arbitrary-network round-trip case because arbitrary networks never enter the
API.

Wrap `new_dummy` in a request-local memoizing allocator: each numeric position
calls the caller's `FnMut` at most once, and every fixed-point iteration reuses
the cached `Aind`. Recompute the representation/ordinal-to-position assignment
from the surviving canonical line order on every iteration while retaining the
position-indexed callback cache. If an earlier line disappears after zero
substitution, a surviving line moves into that earlier position instead of
retaining a stale representation-specific number. Product factors and
materialized Power copies use disjoint simultaneous ranges. Mutually exclusive
Sum branches reserve the shared interface prefix and may reuse the same
branch-local suffix, sized to the maximum surviving branch. This keeps a
stateful factory deterministic across the iterations that are actually
required.

`IndexTooling::canonize` remains the Atom-facing facade, but its implementation
becomes:

1. reversibly encode structured indices;
2. validate the cooked Atom's head-local tensor symmetry grammar;
3. parse it once into a `CanonicalPolicyNet`, retaining that exact cooked Atom;
4. run the network canonicalization driver through its conditional
   Graphica/reconstruct/execute/reparse iterations;
5. restore the reversibly encoded structured indices in the stable Atom
   retained from the final execution;
6. return that uncooked Atom.

Validation must cover every canonical Atom ingestion route, including dynamically
created tagged symbols and the Python-facing constructors. An intrinsically
symmetric tensor occurrence whose immediate arguments are not all direct slots
returns a typed error. Parsing must preserve the attributes the user declared
on each head and must never manufacture `Linear`.

The internal graph projection and reconstruction stay crate-private. Do not add
a raw-network entry point, second public canonicalization mode, or old-path
fallback. Migrate all in-tree callers to obtain their canonicalization input
from the canonical Atom parser. A caller that possesses only a manually built
network must first obtain a normalized source Atom through its ordinary owning
workflow and let the canonical parser construct a new wrapper from that Atom.
It must discard the old network rather than pairing it with the Atom;
canonicalization itself does not execute that network to manufacture
provenance.

Within each canonicalization iteration, ordinary `Network::execute` runs only
after the complete network has been projected, globally canonized, and rebuilt.
Execute a clone with the same `Sequential`/`SmallestDegree<()>`, dummy library,
and `Wrap` semantics as `simple_execute()` so the unexecuted reconstruction
remains available until its normalized reparse is accepted. Product contraction
and ready-operation replacement may make the executed clone opaque, but they
cannot affect a labeling decision that has already been made. For symbolic
tensors this execution is precisely Atom multiplication, addition, negation,
Function construction, and Power evaluation. The current `simple_execute()`
convenience method unwraps errors; the fallible driver should call the same
existing `Network::execute` plus result-extraction path directly and map failures
to `CanonicalizationError::Execution`, not introduce a second renderer.

## Algorithmic Boundary

The design follows the global target of Butler--Portugal and graph-based tensor
canonicalizers: every tensor and contraction relevant to one canonical object
must be visible when its canonical labeling is chosen. Butler--Portugal/xPerm
make that global statement for one monomial and normally handle sums termwise.
SeQuant globally rebuilds each complete Product region but makes Sum and
Adjoint results opaque to a parent Product. GammaLoop instead keeps the
expression nodes visible, so one automorphism can include an index swap, an
outer antisymmetric tensor, and a permutation of two Sum branches.

The relevant baselines are the global signed double-coset formulation of
[Butler--Portugal/xPerm](https://arxiv.org/abs/0803.0862), the refinement and
symmetry-propagation improvements of
[Niehoff](https://arxiv.org/abs/1702.08114), and canonical graph labeling in the
[McKay--Piperno](https://arxiv.org/abs/1301.1493) family used by Graphica. In
the inspected SeQuant implementation, a Product recollects every tensor factor
below that Product but stops at Sum and Adjoint boundaries before rebuilding a
complete flat tensor network
([factor collection](https://github.com/ValeevGroup/SeQuant/blob/c44e41abfe7d2a18f30916dc31db4cead970ea10/SeQuant/core/eval/eval_expr.cpp#L359-L380),
[Product rebuild](https://github.com/ValeevGroup/SeQuant/blob/c44e41abfe7d2a18f30916dc31db4cead970ea10/SeQuant/core/eval/eval_expr.cpp#L538-L588)).
SeQuant therefore respects globality within each Product region but represents
a Sum as an opaque nonsymmetric intermediate to its parent
([Sum IR](https://github.com/ValeevGroup/SeQuant/blob/c44e41abfe7d2a18f30916dc31db4cead970ea10/SeQuant/core/eval/eval_expr.cpp#L432-L479)).
In that revision, its Bliss call requests no automorphism callback and records
only the phase of the selected canonical slot map
([canonical form](https://github.com/ValeevGroup/SeQuant/blob/c44e41abfe7d2a18f30916dc31db4cead970ea10/SeQuant/core/tensor_network/v3.cpp#L714-L721),
[phase](https://github.com/ValeevGroup/SeQuant/blob/c44e41abfe7d2a18f30916dc31db4cead970ea10/SeQuant/core/tensor_network/v3.cpp#L843-L890)).
It is therefore not the source of GammaLoop's negative-stabilizer algorithm;
GammaLoop keeps the operation context and consumes Graphica's generators.

Symbolica's own tensor canonicalizer strips top-level factors containing none
of the declared tensor indices before canonizing one Product, multiplies those
constants back afterward, and handles a top-level Sum term-by-term
([Product constants](https://github.com/symbolica-dev/symbolica/blob/9135b4091505276baa11ef262fe1165f2b95ea33/src/tensors.rs#L180-L215),
[top-level Sum](https://github.com/symbolica-dev/symbolica/blob/9135b4091505276baa11ef262fe1165f2b95ea33/src/tensors.rs#L123-L154)).
GammaLoop does not generalize that term-local optimization into a hidden scalar
payload: scalar leaves stay visible because one global automorphism may cross a
Sum or Function boundary. Only the mathematical `±1` phase induced after
unsigned graph canonization is carried transiently during reconstruction; it is
resolved into ordinary leaf values or a created/reused Neg scope before
execution.

Graphica supplies:

- the canonical unsigned graph;
- the input-to-canonical vertex map;
- generators for the graph automorphism group.

The network canonicalizer supplies the algebra that Graphica intentionally does
not know:

- the phase of every induced antisymmetric slot permutation;
- the location at which that phase is reconstructed;
- the signed action of graph automorphisms on products, sums, functions, and
  powers;
- the semantic scope at which a negative stabilizer implies zero.

Network execution must not build the Graphica input: it would walk and collapse
the operation tree before global labeling. Use one direct projection instead.
Ordinary symbolic execution is restricted to two disposable contexts:
temporary post-Graphica scope clones used for exact group comparisons, and one
clone of each canonical rebuilt network used to materialize its complete
normalized Atom. The initial Graphica input is projected directly from the
proof-carrying canonical-parser result, not from an executed network. No
executed leaf is projected directly. Opaque fragment compression and
boundary-coset caching are out of scope; no executed leaf or cache entry
substitutes for visible descendants in the unified graph.

## Unified Expression-and-Incidence Graph Projection

Use Symbolica's tensor graph implementation as the behavioral specification,
not as a requirement to copy its Atom syntax tree. The parser-produced network
inside `CanonicalPolicyNet` and the induced Graphica graph stay at Spenso's
existing granularity:

- explicit graph nodes represent products, sums, negations, powers, enclosing
  functions, and scalar leaves;
- one `LocalTensor` leaf owns a flat `OrderedStructure`, its slot hedges, and the
  exact `SymbolicTensor.expression` Atom;
- one projected tensor occurrence has one tensor vertex and one vertex per flat
  structural port;
- one typed line vertex represents each index component derived from network
  connectivity.

Product and Sum children are unordered but remain distinguished by their
operator context. `NetworkOp::Function` is unary and has one visible Argument
edge; multiargument Symbolica head symmetry belongs only to a `LocalTensor`'s
slot-erased layout. Negation is unary. A materialized Power has a signed result
node, a typed `PowerMagnitude` child, and typed base/copy roles; copy incidence
depends only on the unsigned magnitude. These operation nodes
prevent tensor occurrences in alternative Sum
branches or nonlinear arguments from being mistaken for simultaneous Product
factors, while keeping cross-branch automorphisms visible.

A visible root marker and directed, typed parent-to-child edges make the
expression orientation recoverable. The projector takes a non-mutating
flattened view of associative Product and Sum nodes so that parser or builder
association does not color the graph. Normalize an empty Product to one, an
empty Sum to zero, unary Product/Sum nodes to their child, `Power(0)` to one, and
`Power(1)` to its base; anchor scalar-only roots with the same root marker.

Line vertices come from the network's hedge connectivity, never from equal
printed index names. An external line retains its complete `LibrarySlot`; an
internal line retains the normalized contraction group and the endpoint
representations on its ports. Uses of one interface in alternative branches may
share the corresponding line vertex when the network connects them through the
operator boundary, but the branch nodes preserve their additive meaning.

Global visibility means every `NetworkOp` and every recognized structural slot
is present. Scalar leaves and opaque tensor-layout fields remain atomic only
after validation proves that they own no network-contracted slot. A tagged
ordered argument or multiargument scalar function that would hide recognized
tensor topology is rejected with a typed pre-Graphica grammar error unless an
explicit bracket/projector opacity boundary contains that topology. A lexical
symbol that merely resembles a dummy remains ordinary opaque data. Their
visible colors use a stable semantic Atom/symbol key based on qualified names,
attributes, and recursive payload structure, never raw Symbol interner IDs.
Apply the same stable ordering rule to function heads and representation colors
so graph canonicalization is independent of interner history.

Use one shared stable semantic-key contract for graph colors,
stability-certificate syntax classes, and signed decoration classes after first
confirming that no existing helper provides it. Give Atom variants a fixed explicit precedence;
encode exact canonical numbers, qualified symbol names and semantic attributes,
ordered Function arguments, and Power base/exponent recursively; sort Sum and
Product child keys by this same total order rather than Symbolica's internal
term order. Parameterize structured slots by canonical line handles for scope
comparisons. A hash may index these keys but never substitutes for full key
equality or ordering.

Build line components as an explicit equivalence closure over network slot
hedges:

1. union the two halves of every paired `NetworkEdge::Slot`;
2. at each operation node, union only crown slot hedges that represent the same
   declared operation-interface `LibrarySlot`;
3. never union tensor ports or branch-local components merely because their
   printed index names compare equal elsewhere in the graph.

A component containing a dangling/self-paired network slot is external and must
carry one consistent complete `LibrarySlot`. Every internal component must have
one compatible normalized representation group, while its ports retain their
endpoint orientations. Conflicting external identities, representations, or
operation-interface declarations return a typed projection error. This closure
is connectivity plumbing only; its representative IDs remain hidden.

The arguments inside a `LocalTensor` do not become additional expression nodes.
Instead, the tensor-vertex color carries a slot-erased argument layout, and
hidden origin data retains the information required to fill its holes again.
In particular, `p`, `cyclic(...)`, and `spenso::bracket(...)` remain in that
layout; they are not separate graph vertices.

### Tensor layout and tensor-local incidence fragment

For example, parse

```text
F(p, cyclic(a, b))
```

as one `LocalTensor` with structure `[a, b]` and exact expression
`F(p, cyclic(a, b))`. The tensor-local portion of its Graphica projection is
equivalent to:

```text
tensor vertex color: F(p, cyclic(□, □))
tensor -> port(a)  [group = argument 1 / cyclic]
tensor -> port(b)  [group = argument 1 / cyclic]
port(a) -> port(b) -> port(a)  [cyclic marker]
port(a) -- line(a)
port(b) -- line(b)
```

There is no intermediate `cyclic` vertex. The slot-erased layout distinguishes
this tensor from `F(q, cyclic(□, □))`, `F(cyclic(□, □), p)`, and
`F(p, sym(□, □))` without putting dummy-index identities into the tensor color.
A bracket boundary is likewise retained in the layout so that reconstruction
cannot flatten it.

A mixed tensor makes the division between outer ordering and group symmetry
explicit. For

```text
F(p, a, sym(b, c), d, antisym(e, f), cyclic(g, h, i))
```

the tensor color contains the ordered, slot-erased layout

```text
arg[0] = Opaque(p)
arg[1] = DirectSlot(□)
arg[2] = SymmetricGroup(□, □)
arg[3] = DirectSlot(□)
arg[4] = AntisymmetricGroup(□, □)
arg[5] = CyclicGroup(□, □, □)
```

and its tensor-local incidence fragment contains

```text
tensor -> port(a) [visible = Ordered(arg[1])]
tensor -> port(b) [visible = SymmetricGroup(arg[2]), hidden member = 0]
tensor -> port(c) [visible = SymmetricGroup(arg[2]), hidden member = 1]
tensor -> port(d) [visible = Ordered(arg[3])]
tensor -> port(e) [visible = AntisymmetricGroup(arg[4]), hidden member = 0]
tensor -> port(f) [visible = AntisymmetricGroup(arg[4]), hidden member = 1]
tensor -> port(g) [visible = CyclicGroup(arg[5]), hidden member = 0]
tensor -> port(h) [visible = CyclicGroup(arg[5]), hidden member = 1]
tensor -> port(i) [visible = CyclicGroup(arg[5]), hidden member = 2]
port(g) -> port(h) -> port(i) -> port(g)  [cyclic markers]
each port(x) -- line(x)
```

There is no ordering edge from `a` to `b` or `c`. Outer order is absolute:
the layout places `DirectSlot(arg[1])` before `SymmetricGroup(arg[2])`, and
those distinct visible incidence roles prevent any automorphism from exchanging
`a` with `b` or `c`. The two group members share the same visible `arg[2]`
role, so they may exchange. Hidden member positions are used only for
the input member frame and exact reconstruction path, not as the canonical
coordinate system for automorphism parity. Only the cyclic group adds
port-to-port edges.

Each tensor-to-port incidence edge has a typed role:

- an ordered slot has a visible argument path or position;
- members of one symmetric, antisymmetric, or cyclic group share a visible,
  stable group key derived from the layout path;
- the original member position is hidden input-frame and reconstruction data;
- distinct partial groups have distinct group keys, so their slots cannot mix;
- cyclic members additionally have directed cycle-marker edges;
- antisymmetric groups additionally contribute their induced permutation to
  the sign calculation.

Group keys must come from stable structural paths, not occurrence IDs. Symmetry
and parity belong to the group recorded on the tensor-local incidence edges;
they are not inferred from all ports incident to the outer tensor vertex.

An external line vertex retains its complete `LibrarySlot`. An internal line
vertex is colored by normalized representation group and never by its dummy
name. Port vertices retain endpoint representation and orientation. A line may
touch ports in mutually exclusive Sum branches; it represents network index
identity and interface flow, not a simultaneous multiway contraction.

### Mapping Spenso tensors into the graph

For each `LocalTensor` occurrence:

1. classify the tensor's own head and validate the direct-slot-only rule when it
   carries intrinsic symmetry;
2. scan `SymbolicTensor.expression` to build a slot-erased layout and record each
   supported slot occurrence's exact argument path, ordered/group role, group
   kind, and original member position;
3. extend generic structure inference so direct children of a recognized
   structural `sym`/`antisym`/`cyclic` group contribute to the flat
   `OrderedStructure`, without creating a nested network node;
4. match those occurrences positionally to `external_structure_iter()` and to
   the occurrence's slot hedges sorted by `network.graph.slot_order`; do not use
   Atom equality as occurrence identity because repeated slots are valid;
5. create one tensor vertex plus its flat port vertices, add the typed incidence
   and cyclic-marker edges, and connect every port to the single line vertex
   derived from the existing Spenso hedge and operation-boundary connectivity;
6. retain the tensor/store origin, slot hedge, expression path, layout hole,
   symmetry-group scope, and original member position needed for signed
   reconstruction.

The existing trace/projector interpretation of `sym`/`antisym`/`cyclic` remains
unchanged. Only their recognized use as a structural partial-slot argument is
flattened into a `LocalTensor` structure. A slot hidden under any other nested
head or under `spenso::bracket(...)` is not an exposed structural slot.

The current `TensorColor { head, arguments }` shape is not enough for this: it
would retain `cyclic(a, b)`, including dummy identities, as opaque color data.
Replace that role with the reversible slot-erased layout after first confirming
that existing Spenso structure metadata cannot already express it.

Sums and enclosing functions remain explicit in this same graph; they are not
canonicalization barriers. For every parser-supported `Power(n)` with `n != 0`,
let `m = usize::from(n.unsigned_abs())` and project `m` copies of its complete
base fragment. Retain the exact signed `n` on the outer Power-result node, add
one typed `PowerMagnitude` child representing
`D = product_(k=1..m) B_k`, and attach the interchangeable copy roots beneath
that magnitude scope. Thus the incidence symmetries of multiplication and
reciprocal multiplication are shared without identifying their semantics or
their zero scopes. Give each occurrence an explicit `PowerCopy` root connected
to exactly one complete cloned base root.
Expression membership must force every automorphism to map the complete vertex
set below one copy root to the complete vertex set below one target copy root,
and reconstruction validates this property before folding. Normalize
`Power(0)` to one before materialization.

With canonical `ParseSettings::default()`, a purely scalar division such as
`x^-2` is normally precontracted into one visible scalar Atom leaf; preserve its
exact normalized Atom payload rather than manufacturing a structural Power.
The copy construction applies to materialized Powers of tensor-derived scalar
networks and supported self-dual tensors.

Clone base-internal lines independently per occurrence. Materialize the current
Spenso self-dual Power interface exactly, using one complete-copy pairing across
every exposed base slot:

- for even `m`, pair copies `(0,1)`, `(2,3)`, and so on, connecting every pair's
  corresponding boundary ports through fresh two-ended internal line vertices;
- for odd `m`, pair all but one copy identically and connect the remaining
  copy's boundary ports to the original external interface;
- scalar bases have no boundary-pairing fragment;
- non-self-dual tensor Powers remain unsupported, and the parser admits a
  negative self-dual Power only when `m` is even and the paired result is
  scalar before inversion.

The initial pair/copy ordinals are hidden because equivalent pairings of
identical copies are isomorphic, while the resulting two-ended line incidence
is visible. Never identify all corresponding ports into one multiway line or
choose partner copies independently per slot. This exposes the valid pair,
pair-of-pairs, internal-base, and odd external-copy actions without distributing
any Sum inside the base.

Retain hidden Power/copy origins. For every canonical copy, derive a descriptor
containing its canonical copy-root vertex, complete vertex membership,
canonical boundary pairing, normalized base `ClassKey`, and residual whole-copy
phase. Reconstruction may fold the copies back only after verifying that every
copy-root image is complete, all operation/tensor topology and internal/boundary
line pairings agree, and their reduced base decorations represent one common
semantic base up to a whole-base scalar sign. Choose that common base by the
same exact signed-orbit rule used for map independence. A non-scalar, localized,
partial-copy, or interface mismatch is a typed Power reconstruction error. If
copy `k` differs from the common base by `q_k`, rebuild one network Power node
with the original signed exponent `n` and, when `XOR_k q_k` is odd, wrap that
Power in a created/reused exact-scope Neg.
This residual must remain outside the base: for even `m`, one odd copy still
gives `-B^n`, whereas `(-B)^n` would lose the minus for either sign of `n`.
Inversion does not alter a boolean phase, so the XOR rule is identical for
positive and negative exponents. A sign already common to every copy may remain
inside the chosen base, with the outer residual adjusted by the parity of `m`.

Estimate the complete expanded vertex and edge counts before allocating the
Graphica graph. Apply a deterministic whole-canonical-graph resource budget,
because cost depends on `m * base_size`, not on the sign of the exponent.
Compute every count with checked unsigned arithmetic. Exceeding
the budget returns a typed `CanonicalizationError::GraphSizeLimit` with Power
context, signed exponent, unsigned magnitude, requested vertex/edge counts, and
configured limits; it must never silently leave the Power opaque. Allow smaller
injectable test limits and select the production limits from Phase 6
graph-size/search benchmarks. A future compact algebraic repeated-copy action
is an optimization only; the initial implementation uses literal complete-copy
expansion.

Execute the reconstructed Power through the ordinary network executor; do not
add a canonicalization-specific Power path. Correct generic execution is
therefore a prerequisite: make `Power(0)` return a valid one for every leaf kind
before `NetworkState::Tensor` or non-self-dual rejection and without reading
`scalar(0)` from a possibly empty tensor-only store. `Network::pow(0)` and
`NetworkState::pow(0)` must likewise return the existing pure-scalar one path
before constructing a Power graph or matching on the base state. Fix the
duplicated odd-magnitude paths whose current repeated-squaring loop maps
magnitude three to five. Compute the positive magnitude with
`i8::unsigned_abs()` and, for a negative exponent, invert the resulting scalar
exactly once through the existing `ref_one() / value` semantics. Preserve the
parser's scalar-result restriction instead of attempting tensor inversion.
Lock signed exponents `-5..=5` for scalar leaves, the supported nonnegative and
negative-even cases for `LibraryKey`, local-tensor, tensor-sum, scaled-tensor,
and scaled-tensor-sum leaves, and supported scalar-result `i8::MIN` cases before
canonicalization relies on them. The normalized executed form is reparsed;
`Power(0)` becomes one, `Power(1)` becomes its base, and negative Powers become
ordinary Symbolica division rather than an unsupported subtree or panic.

After building the complete graph, call `Graph::canonize()` exactly once for
that whole-network iteration. The canonical order of operation nodes, tensor
occurrences, ports, and line vertices is therefore chosen together. Signed
reconstruction returns a fresh canonical network. Resolve transient parity into
fresh leaf values, nested Atom paths, or exact-scope Neg nodes, assign
deterministic dummies, and execute a clone through the ordinary symbolic network
executor.

Reparse the still reversibly encoded executed Atom on every iteration. Do not
attempt to decide stability by hashing or comparing two independently labeled
visible graphs: ignoring association, commutative ordering, dummy spelling, and
hidden origins while retaining complete global wiring would itself require a
graph-isomorphism decision. Instead, reconstruction builds a one-sided
stability certificate from Graphica's selected `vertex_map`. It may certify
stability only when all of the following hold:

1. the iteration input is a `CanonicalPolicyNet` paired with the exact
   normalized Atom accepted by the fixed canonical parser; on retries that Atom
   is the result of the preceding ordinary execution;
2. canonical transport preserved visible operation topology, semantic colors,
   scope boundaries, exact signed Power-result structure, and every
   `PowerMagnitude`/complete-copy relation;
3. reconstruction introduced no sign, zero substitution, residual Power phase,
   exact-scope Neg, or scalar/tensor/nested leaf-payload edit; and
4. the selected map transports bijectively the equality partitions of
   contracted-index occurrences, Symbolica-visible scalar and tensor leaves,
   and bottom-up Product, Sum, Function, and Power child syntax classes.

The fourth condition catches an unsigned slot or dummy transport that makes two
formerly distinct Sum children syntactically equal and therefore combinable by
Symbolica. This check is a verifier under the already selected map, not a
second graph canonicalizer or an Atom renderer. A merge, split, disappearance,
value change, or ambiguous class correspondence declines to certify stability;
failure to certify is not an error and conservatively starts another complete
iteration after execution/reparse unless execution returned the exact input
Atom. Association, commutative ordering, declared
unsigned symmetry transport, and a bijective dummy rename are harmless when the
certificate proves that every visible class is otherwise preserved. Reuse the
request-local dummy allocator and canonical parse policy. Exact executed-Atom
equality proves a fixed point, while a terminal zero or one returns immediately.

Canonical vertex transport and dummy renaming do not by themselves cause a
second Graphica call. A generated sign may: direct negation of a leaf or
execution of an exact-scope Neg can reparse as a different visible scalar or
Product form, such as `2 * (-A) -> -2 * A`. In that case rerun the whole graph;
do not promise that every sign-only input is a one-call case. An unsigned input
requiring only relabeling may still use one Graphica call. Do not run an
unconditional second Graphica verification pass.

Retain the exact canonical-problem identity after every Graphica call; a hash
may index candidates, but full canonical graph/key equality confirms a match.
If a conservative retry immediately produces the same identity as the preceding
iteration, that Graphica call has supplied the proof the one-sided certificate
could not: return the normalized network and Atom retained from the preceding
execution without reconstructing or executing again. If an identity matches an
older non-consecutive state, determinism of canonical frames, orbit minima,
dummy allocation, execution, and reparsing makes it a genuine cycle; return a
typed `CanonicalizationError::ConvergenceCycle` with the first/repeated
iterations, cycle length, and recorded retry reasons.

Do not claim a monotone fixed-point measure across sign recoloring, zero/Sum
reduction, Function linearity, and Power normalization without proving one.
Enforce a deterministic iteration budget as an additional fail-safe and return
`CanonicalizationError::IterationLimit` with the limit and last retry reason if
distinct states exhaust it. Allow the crate-private/test driver to inject a
smaller budget; select and document the production default from the Phase 6
iteration benchmarks before migrating callers. Record diagnostic retry reasons
such as sign edit, zero substitution, residual Power phase, visible-class
merge/split, execution-topology change, and incomplete stability certificate.

This repetition is necessary because a branch later proven zero can initially
break an automorphism of its siblings; that automorphism is absent from the
original `Aut(G)` and becomes discoverable only after rebuilding the reduced
whole graph.

## Group-Aware Signed Reconstruction

Signed reconstruction has two distinct inputs and responsibilities:

1. the `vertex_map` determines how this input is transported to Graphica's
   selected canonical labeling;
2. the automorphism generators prove that this selected signed decoration is
   independent of the chosen map, or identify a negative stabilizer and its
   zero scope.

Do not conflate the second responsibility with checking whether one returned
generator is globally odd. That is sufficient only when the relevant action is
already a one-dimensional sign character, as for one flat multiplicative
scope. A generator may permute local sign scopes, and an odd element fixing one
scope can be a word in the generators.

### Transporting the selected canonical map

For every intrinsic or partial antisymmetric group, record a stable sign site
and its exact semantic sink. After Graphica assigns canonical vertex numbers,
give every canonical symmetry site a fixed local port frame. For symmetric and
antisymmetric sites, sort that site's port vertices by canonical vertex ID. For
a cyclic site, start at its smallest canonical port vertex and follow the
visible directed cycle. Identify the site by its canonical tensor-header vertex
and stable group-layout path; repeated slots remain distinct because their port
vertices are distinct even when they meet the same line.

For selected-map transport, the input frame is the group's original member
order. Map those ports through `vertex_map`, express them as positions in the
target canonical frame, and take that permutation's parity. For a canonical
automorphism, express the image of the canonical source frame in the canonical
target frame. Thus input occurrence order contributes exactly once during
input-to-canonical transport and never becomes the coordinate system of the
canonical generated-group action. Reconstruction values carry the selected-map
boolean only transiently and XOR parities at the same sink before changing the
fresh rebuilt network. Reconstruction then:

- negates the fresh scalar `Atom` or `SymbolicTensor.expression` for a phase
  owned by that leaf occurrence, never a shared store entry;
- edits a partial-group sign at its exact stored Atom path, for example
  `F(p, -antisym(a, b))`;
- after preserving nested/nonlinear ownership, flattens associative Product-only
  descendants and considers only leaf occurrences reached entirely through
  Product nodes; order eligible occurrences by `(scalar before tensor,
  canonical vertex ID)`, sink the sign to the first, and clone that store value;
- represents a phase owned by an operation result by creating or reusing
  `NetworkOp::Neg` at that exact scope when no eligible Product-local leaf
  exists or sinking would cross Sum, nonlinear Function, or Power; ordinary
  execution later lowers it onto the scalar or tensor result leaf;
- combines signs already represented by an explicit Neg and removes even
  pairs;
- leaves a sign inside a nonlinear function argument;
- allows ordinary Symbolica semantics to lift it only through a head declared
  `Linear`;
- folds equivalent materialized Power copies to one common base and puts the
  XOR of their residual whole-copy phases in a Neg outside the Power, where it
  cannot be lost for an even unsigned magnitude, including a negative even
  exponent;
- rebuilds each Sum child with its own leaf-local or exact-branch sign, while
  retaining every pre-existing scalar child as an ordinary visible Atom that
  may itself receive a phase when it is the selected sign owner.

This placement makes ordinary operator scope structural rather than a separate
phase-bubbling or coefficient algorithm. Never move a sign merely to reach an
arbitrary descendant leaf: Product is sign-transparent, but `-(x+y)`, `-f(x)`
for nonlinear `f`, and `-(x^2)` demonstrate why the other scope boundaries must
be preserved. The canonicalizer owns sign-site transport, exact
negative-stabilizer proofs, and zero placement. Ordinary Product multiplication,
Sum coefficient arithmetic and cancellation, declared linearity, and Power
normalization belong to symbolic network execution.

First reconstruct with canonical line handles rather than concrete dummy names.
Preserve external indices exactly. Fill ordered layout holes by their visible
role, symmetric and antisymmetric group holes in canonical-frame order, and
cyclic group holes in their canonical frame from the smallest port along the
directed cycle.
Rebuild every partial-group wrapper at its stored argument path. Assign concrete
dummies only after innermost negative-stabilizer zero substitution: shared
interface lines keep one enclosing-scope name, while branch-local internal
lines use a canonical per-term namespace whose counters reset for each mutually
exclusive addend. This makes alpha-equivalent executed terms syntactically
identical before Symbolica combines them.

### Interpreting the automorphism group

Retain every Graphica generator as a paired action consisting of its complete
permutation of canonical graph vertices and its induced signed permutation of
sign sites:

```text
g = (pi_g, sigma_g, epsilon_g)
pi_g: graph vertex -> graph vertex
site s -> (sigma_g(s), epsilon_g(s))
```

If `F_s` is the canonical frame of site `s`, the induced member permutation
`q_(g,s)` is defined by expressing `pi_g(F_s)` in the target frame
`F_(sigma_g(s))`; `epsilon_g(s)` is its parity. Never compute this canonical
generator parity from hidden original member positions.

Compose both parts exactly. The vertex part is ordinary permutation
composition, `pi_(h o g) = pi_h o pi_g`. If `g` maps `s` to
`(sigma_g(s), epsilon_g(s))` and `h` maps a site to
`(sigma_h(s), epsilon_h(s))`, then `h o g` maps `s` to
`(sigma_h(sigma_g(s)), epsilon_g(s) XOR epsilon_h(sigma_g(s)))`. Derive
inverses from the same paired law. This is the affine signed-site action on
mathematical boolean decorations attached to the full graph action; it does not
define a stored network payload. Generator order must not affect its generated
subgroup, and no composition, Schreier generator, or stabilizer operation may
discard the full vertex permutation.

The fixed frames make the cocycle law structural:
`q_(h o g,s) = q_(h,sigma_g(s)) o q_(g,s)`. Directly computing the induced
permutation of a composed element must therefore give the same parity as XORing
the two stored phases. Apply the corresponding inverse-frame identity when
deriving inverse elements.

`vertex_map` maps input vertices to canonical vertices; Graphica's
`orbit_generators` are already expressed in canonical numbering. Materialize
all disjoint cycles returned for one generator as one complete permutation of
the canonical vertex domain, rather than as independent generators, and attach
the translated site descriptors once. This is the contract of the Graphica
revision currently pinned in `Cargo.lock`
(`9135b4091505276baa11ef262fe1165f2b95ea33`) and must be protected by a
dependency-contract test.

The expression graph maps each site to its owning tensor, partial-group path,
term, and surrounding operator ancestry. For a semantic scope, its exposed
boundary is the ordered tuple of line vertices incident both inside and outside
that scope, together with external lines. Local reasoning uses the subgroup
that fixes the scope root and this boundary pointwise. A setwise boundary
permutation belongs to an enclosing action unless an explicit interface
symmetry declares it local.

The pointwise test is performed on `pi_g`, not on the signed-site projection.
Two elements may induce the same odd self-action on one antisymmetry site while
one fixes its exposed lines and the other exchanges them. For example, swapping
the two free lines of `A_antisym(i,j)` expresses antisymmetry and does not make
that tensor zero; an odd action whose complete scope boundary is fixed can prove
zero. Keep the full action until the relevant scope subgroup has been selected,
and only then project to the child/site action used for scalar `+identity` or
`-identity`.

Process scopes bottom-up and analyze the generated action, not only the raw
generators. Here "scalar action" means `+identity` or `-identity` on the complete
semantic value; it never means extracting an existing scalar graph leaf:

- Product retains the induced permutation of equal child values. Once that
  permutation is quotiented by ordinary commutative Product semantics, child
  signs XOR to `+identity` or `-identity` on the complete Product.
- Sum retains an exact signed permutation of live child scopes relative to its
  boundary solely for group-theoretic `-identity`, map-independence, and outward
  propagation. It does not strip coefficients or aggregate terms; direct
  symbolic network execution owns that algebra.
- Neg transports the child's action unchanged relative to the negated value.
- A nonlinear Function is a boundary for outward sign and zero propagation.
  Analyze its argument relative to the pointwise-fixed function interface, but
  reconstruct a proven zero argument as `f(0)` rather than zeroing the function.
- A declared `Linear` Function is transparent for sign and zero propagation.
- For either sign of `Power(n)`, analyze the magnitude scope
  `D = product_(k=1..m) B_k` on `m = |n|` materialized copies. A stabilizer may
  prove `D = 0` without proving an individual base copy zero, so do not
  strengthen a magnitude-scope result into a base result merely to fold it.
- For `n > 0`, replace a proven-zero `D` by zero for every positive exponent,
  including even `n`.
- For `n < 0`, the value is `D^-1`; a proven-zero `D` is a denominator
  singularity. Reconstruct a canonical scalar Atom zero leaf, not the empty-Sum
  `Network::zero()` form, under an ordinary reciprocal Power, decline the
  stability certificate, and let whole-network symbolic execution produce the
  backend-defined normalized singular value.
  It must never become a zero result merely because the same magnitude scope
  would zero a positive Power. A negative Power is therefore a barrier to
  outward zero propagation, not to signed group analysis inside its
  denominator.
- When the magnitude scope is not proven zero, fold copies that differ only by
  scalar phase to one exact common base and retain the XOR of copy-relative
  phases as an exact Neg outside the signed Power. Ordinary execution owns
  inversion and the parity of any sign retained uniformly inside the common
  base; that parity is `m mod 2`.

Compute the required pointwise stabilizer generators exactly. First check for
existing permutation-group support; if it is insufficient, add the minimum
composition, inverse, orbit/transversal, strong-generator, and
stabilizer-chain/Schreier operations after maintainer confirmation. Never fall
back to checking only Graphica's raw generators. Every transversal, inverse,
strong generator, and Schreier word retains the paired full-vertex and signed
site action. Form the subgroup fixing the scope root and each exposed boundary
vertex pointwise using the full vertex permutation; only within that subgroup
project to the induced child/sink classes and sign character. Use the exact
generated action to obtain the subgroup stabilizing the normalized child/sink
classes; checking each strong generator independently is not a substitute for
this subgroup calculation. An odd generator that merely maps one sink, scope,
or boundary vertex to another proves neither zero. A local zero requires an
element fixing the scope root and boundary whose induced action on the entire
reduced semantic value of that scope is scalar `-identity`. When several
ancestors qualify, replace the innermost qualifying scopes first and then
execute the whole-network iteration.

For map independence, let `f` be the selected input-to-canonical map and `H` the
relevant pointwise stabilizer. Canonicalize the exact orbit of normalized
sink-decoration classes with two distinct stable keys:

```text
ClassKey = (
  canonical scope kind and root,
  ordered canonical boundary-line descriptors,
  stable semantic key of the executed normalized scope value,
)

RealizationKey = sorted [(canonical SinkKey, phase), ...]
```

`SinkKey` identifies a leaf by canonical leaf/header vertex plus its exact
stable layout or Atom path, or an operation result by canonical operation
vertex plus scope role. Before forming `RealizationKey`, XOR phases at one sink,
apply deterministic Product sign sinking, combine exact-scope Neg pairs, and
preserve Sum, nonlinear Function, and Power boundaries. To form `ClassKey`,
reconstruct that normalized decoration as an ordinary scope fragment on
canonical line handles, execute a clone fallibly with deterministic canonical
line placeholders, and encode the resulting Atom with the shared stable
semantic key. Map execution/result-extraction failures to
`CanonicalizationError::Execution` with the scope identity; do not introduce a
separate scope renderer.

Choose the exact lexicographic minimum `(ClassKey, RealizationKey)` over
`{h . d | h in H}`. The primary key selects the normalized mathematical value;
the secondary key chooses one concrete sign placement when several decorations
execute to the same value. Orbit enumeration, generator order, generating-set
choice, hash iteration, original dummy names, hidden origins, concrete
`new_dummy` results, and raw Symbol interner IDs must not affect this minimum.
A sign site inside an already proven-zero group, expression, or Power magnitude
is unobservable: erase those coordinates before enumeration and require the
active-site mask to be invariant under every retained site transport. Enumerate
the resulting observable quotient orbit, not all affine assignments to erased
coordinates, so independent signs below zero scopes cannot cause exponential
work. Preserve the full-length decoration with erased coordinates fixed to
false for stable site addressing during reconstruction.
A stabilizer-chain or branch-and-bound implementation may avoid explicit orbit
enumeration only if it is proven to return the same minimum as exhaustive
ordering.

Compute the exact class stabilizer

```text
K = {h in H | ClassKey(normalize(h . d*)) == ClassKey(normalize(d*))}
```

using primary semantic equality, not the secondary concrete placement. An
element that moves a Product-wide minus between equivalent leaves may change
`RealizationKey` while remaining in `K`. Only on `K`, where the remaining action
is scalar, is the parity map a one-dimensional sign character; checking that
character on a strong generating set for `K` proves it for every word. An odd
character value is the `-identity` case at that scope. Do not use the first
orbit element encountered, the full representative key as the class
stabilizer, or base-versus-generator comparisons before forming this exact
semantic class stabilizer. Return `AmbiguousSignScope` only when exact semantic
class representatives remain inequivalent, not merely because a sign site or
concrete sign sink moves.

## Sums, Functions, Execution, and Zero Propagation

A Sum and all of its descendants participate in the global Graphica call.
During bottom-up signed reconstruction, analyze its exact signed child action
relative to the subgroup fixing the Sum root and exposed boundary pointwise.
An action that permutes that boundary is a symmetry of the whole Sum value to be
propagated to its parent; it must not be used as a local fixed-boundary action.
Generator compositions remain essential for proving `-identity`, even though
ordinary addition is delegated to Symbolica.

Do not construct coefficient-free body keys, strip scalar factors, transport
coefficients with parity orbits, or run a local graph canonicalizer. Existing
scalar children remain visible in the whole graph. Reconstruct every surviving
addend with its signed leaf values or exact-branch Neg and canonical line
handles, then materialize branch-local dummies from a counter namespace reset
for each addend. Ordinary execution adds the resulting scalar/tensor leaves;
Symbolica normalization therefore combines identical terms, performs its
supported coefficient arithmetic, cancels opposites, and removes zero addends.
Do not promise coefficient collection beyond Symbolica's actual normal form.

If executed addition can delete or merge visible branches, reconstruction
cannot certify stability and the normalized reparse starts another whole-network
iteration because the smaller graph may have new automorphisms. If canonical
transport only orders the Sum and preserves its bottom-up child equality
partition, the certificate may skip another Graphica call. This must work at
any depth, giving

```text
f(x + zero) -> f(x)
```

without requiring `f` to be linear.

This additive simplification is distinct from zero in a function argument:

```text
f(x + zero)                  -> f(x)
F(p, antisym(a, a))          -> F(p, zero)   # F is nonlinear
F(p, antisym(a, a))          -> zero         # F is user-declared Linear
F(p, antisym(b, a))          -> F(p, -antisym(a, b))
F(p, antisym(b, a))          -> -F(p, antisym(a, b))  # F is Linear
```

The final two examples describe the same nested sign with different outer-head
attributes. `F(spenso::bracket(p + q), antisym(a, b))` keeps `p + q` atomic even
when `F` is linear; without the bracket, the user has explicitly requested
Symbolica's ordinary linear distribution.

An odd boundary-fixing stabilizer of one addend zeros that addend, not the
complete Sum. A boundary-fixed signed exchange of addends can instead make two
reconstructed terms cancel. An exchange that also moves the Sum interface is
retained as an action of the complete Sum and propagated outward. Such a global
automorphism can span the Sum and its surrounding context; for example,
swapping `i` and `j` and the two branches proves

```text
A_antisym(i, j) * (B(i) * C(j) + B(j) * C(i)) -> zero
```

because the outer Product receives one uniform minus sign. Never XOR sign
contributions across additive children without first deriving this scoped
action.

Do not globally expand products containing sums during canonicalization. The
example above is visible without distribution because the Sum node, both
branches, their shared interface lines, and the outer tensor are in one graph.
Zeros requiring an algebraic distributive rewrite beyond graph automorphism
remain the responsibility of an explicit expansion pass.

## Applying the Canonical Result to the Network

The canonical graph, `vertex_map`, paired full-vertex/signed-site action, and
hidden reconstruction origins should produce a fresh canonical `SymbolicNet`
rather than partially mutating the original network. Boolean sign decorations
exist only in the group/reconstruction calculation. Resolve them while building
fresh store
entries: negate a scalar `Atom` or `SymbolicTensor.expression` owned by a leaf,
sink a whole-value Product sign to its deterministic existing scalar/tensor
leaf, edit a nested partial group at its exact Atom path, or retain an
operation-value phase by creating or reusing Neg at that exact scope. Do not
mutate a shared tensor/scalar store slot in place, invent an operation payload,
or synthesize a scalar `-1` Product factor or special signed value type.
Traverse explicit canonical operation vertices and rebuild, in one coherent
operation:

- tensor expressions by filling the reversible slot-erased layouts;
- intrinsic versus ordered head classification and all user-declared
  attributes;
- explicit partial-symmetry heads and bracket boundaries;
- flat tensor structures, their ordered slots, and `slot_order`;
- operation-aware line connectivity, recreating tensor-port and
  operation-boundary hedges before emitting any dummy labels;
- existing network operations for products, sums, Neg, functions, and powers,
  with every local sign resolved at its owning leaf, nested path, or exact Neg
  scope;
- external slots and representation orientation.

Use existing Spenso network constructors and builders where they already own
these invariants. Before adding any new constructor or relabel helper, confirm
that the existing network construction API cannot express the required atomic
operation.

One canonical line may touch outer context and ports in mutually exclusive Sum
branches. Reconstruct that scoped interface through the existing network
builders; do not turn all incident ports into one simultaneous contraction.
Reserve names for shared boundary lines in enclosing canonical-line order, then
reset representation-group counters for purely branch-local lines in each
canonical addend. Alpha-equivalent executed terms consequently receive
identical local dummy names without merging their graph components. The
allocator maps each
canonical `(scope, representation group, ordinal)` to a deterministic numeric
position and reuses the request-local cached factory result across retry
iterations.

Sort commutative Product and Sum children by canonical target vertex, and sort
ports and line vertices by explicit canonical IDs and visible roles. Never rely
on adjacency iteration order. The projector's associative normalization makes
parser or builder association irrelevant.

After every signed reconstruction, clone the rebuilt network and execute that
clone with the ordinary symbolic `Network::execute` path and its existing
contraction strategy. Product contraction multiplies Atom expressions, Sum adds
them, Function uses `Wrap`, Neg clones and negates the resulting scalar or
tensor leaf, and corrected generic Power execution evaluates the exponent.
Ordinary Symbolica normalization consequently combines coefficients and terms,
cancels opposites, and applies sign lifting authorized by a declared `Linear`
attribute. Respect operation dependencies when a negative Power has a
proven-zero denominator: do not short-circuit an enclosing Product to zero
before the reciprocal is executed, because Symbolica must decide expressions
such as `0 * (0^-1)`. No canonicalizer-side lifting through nonlinear heads is
allowed.

Extract the executed result Atom fallibly and reparse it into a fresh
`CanonicalPolicyNet` with the canonical policy. The fresh rebuilt network
contains only ordinary negated leaf values, exact nested Atom edits, and
exact-scope Neg nodes, not a signed payload type. If reconstruction produced the
complete stability certificate, return that normalized network/Atom pair.
Otherwise reproject the complete reparsed pair for the next whole-network
iteration, except that terminal zero or one returns immediately.
Before another reconstruction, compare the new exact canonical problem identity
with the seen states: consecutive equality verifies the preceding result,
non-consecutive equality is a convergence cycle, and the iteration budget
guards a long sequence of distinct states.
The clone consumed by execution is disposable, so its opaque contracted leaves
never replace the stable proof-carrying result or any Graphica input. No subsequent
bare tensor canonicalization, second full-network execution, bespoke rendering
pass, or post-execution hash comparison is allowed. Exact group proofs may
execute temporary scope fragments through the same ordinary path for quotient
comparison.

## Implementation Phases

### 1. Lock reference behavior

- Pin the canonical parser contract before graph work: each supported normalized
  Atom produces a deterministic network with fixed recursive settings and the
  tensor-symmetry validator, and executing/reparsing supported fixtures preserves
  the exact semantic Atom. Confirm by codebase search that no existing type
  already couples that source Atom to its parsed network.
- Keep bare Symbolica `canonize_tensors` as a test-only oracle with deterministic
  dummies. Require exact equality only for the tensor-local overlap; compound
  graphs use semantic equivalence, idempotence, and isomorphic-input equality.
- Cover intrinsic and partial symmetry, stable opaque payloads, line
  representations, nonlinear/`Linear`/bracket behavior, and typed grammar
  rejection.
- Pin ordinary symbolic `Network::execute` behavior for Product, Sum, Neg,
  Function, and Power. Fix `Power(0)` for every leaf kind and every duplicated
  `LibraryKey`/tensor-like odd-Power path first. Handle zero before parser
  non-self-dual rejection, `NetworkState::pow`, or Power-graph construction so
  every base kind produces `Network::one()`. Replace signed `i8::abs()` by an
  unsigned magnitude, make negative scalar-result Powers invert exactly once,
  and add execution regressions for signed exponents `-5..=5` plus low-level
  scalar `i8::MIN`; canonicalization depends on this shared executor. Correct
  the currently swapped `InvalidDotFunction` and
  `NonSelfDualTensorPower` display messages before asserting the typed parser
  failures.
- The codebase audit found no general tensor-leaf materializer or
  contraction-Power helper, and maintainer confirmation for the narrow helpers
  was recorded during plan review. Add one context-aware eager materialization
  boundary for `Function`, nontrivial Power, and self-loop tracing, reusing the
  existing tensor-sum/scaled-tensor primitives and preserving library lookup
  errors. Handle `Power(0)` and `Power(1)` before it. Feed the resulting concrete
  tensor to one shared signed-integer Power evaluator so the five tensor-like
  leaf variants cannot diverge.
- Consolidate the duplicate tensor-sum/scaled-tensor materializers in
  `network/mod.rs` and `network/contract.rs` around one core. Preserve reduction
  order as an explicit caller policy while tests pin the current balanced versus
  linear behavior; deduplication must not silently change floating-point
  associativity or the lazy-sum threshold.
- Classify every Sum child by semantic scalar/tensor shape before choosing an
  accumulator, so mixed `Scalar`, scalar-valued tensor, and scaled-tensor sums
  are independent of child order. Keep Sum's lazy/eager policy. Extract the
  duplicated Product scalar scan used by normal and in-place contraction into
  one analysis result, while leaving their distinct graph rewrites separate.
- Add the tensor-index counterpart of `NetworkLeaf::map_scalar_refs` and one
  shared append/rebase-store path. Use them for tensor shifting, parallel-result
  rebasing, binary and n-ary Sum/Product assembly, and scalar convenience
  operations. Fix the existing `Add<i8>` path that shifts tensor references
  twice and never rebases the appended scalar reference.
- Add the global regressions before refactoring: shared-boundary sums, signed
  branch exchange, associative tree variants, and identical fragments with
  different global wiring.

### 2. Build the unified unsigned projection

- Add the opaque `CanonicalPolicyNet` parser result with private fields and no
  raw-network conversion or mutable projection. Make it the only input accepted
  by the canonicalization driver; do not execute it before projection.
- Keep the reversible tensor-layout and typed-incidence work, but make it the
  tensor-local part of one graph containing the root, every operation, scalar
  and tensor leaves, ports, and line components.
- Keep existing scalar children visible. The unsigned projection has no general
  coefficient or sign payload; mathematical sign decorations do not exist until
  post-Graphica reconstruction and are resolved into ordinary network data.
- Add total associative/unary/root normalization, explicit line-equivalence
  validation, stable semantic color keys, and the specified complete-rooted
  Power-copy encoding. Retain the signed exponent on the scope while deriving
  copy count, parity, boundary pairing, and checked resource estimates from its
  unsigned magnitude. Return `GraphSizeLimit` when the injectable resource
  budget is exceeded.
- Merge `ContextGraph` and `TensorGraph`; remove `has_barrier`,
  `canonize_flat`, canonical boundary IDs, and all per-scope Graphica calls.
- Call Graphica once per whole-network iteration and retain the canonical graph,
  input-to-canonical `vertex_map`, and canonical-numbered generators.

### 3. Rebuild the unsigned network

- Reconstruct operation topology, tensor layouts, port/line connectivity, and
  external orientation directly into a fresh `SymbolicNet`.
- Derive and validate the rebuilt `NetworkState` from graph connectivity before
  applying a folded Power, reusing the existing n-ary Product/state path where
  possible. Do not inherit the stale `SelfDualTensor` state that binary network
  multiplication can retain for a closed scalar contraction.
- Restore operation-boundary hedges before assigning shared and branch-local
  dummy namespaces; sort only by explicit canonical IDs and visible roles.
- Preserve exact opaque payloads, brackets, declared attributes, and signed
  Power structure through reconstruction. A proven-zero negative-Power
  magnitude scope is the explicit exception: rebuild a zero denominator under
  an ordinary reciprocal so whole-network execution, not the canonicalizer,
  determines the singular normalized value.
- Execute a clone of each rebuilt network through the existing fallible
  `Network::execute`/result-extraction path and reparse its Atom. Build the
  one-sided stability certificate under the selected vertex map, including the
  bottom-up visible equality-class check; do not add a structural renderer,
  general graph matcher, or canonicalization-specific Power evaluator.

### 4. Add signed reconstruction and conditional fixed-point reduction

- Map every intrinsic/partial antisymmetric site to its semantic sink, compute
  canonical local port frames and selected-map transport, and retain each
  complete canonical-numbered Graphica vertex permutation together with its
  frame-relative signed site action. Original member positions define only the
  input frame and reconstruction path.
- After checking existing helpers, implement the exact composition, inverse,
  orbit/transversal, strong-generator, and pointwise-stabilizer operations
  required by a stabilizer-chain/Schreier algorithm over that paired action.
  Restrict scope roots and boundaries using the full vertex permutation before
  projecting to child/site actions. Raw-generator-only or signed-site-only zero
  detection is not an allowed intermediate production behavior.
- Process scopes bottom-up using pointwise boundaries, exact signed Product and
  Sum child actions, Function barriers, sign-sensitive Power zero semantics,
  and magnitude-copy-relative Power phases. Positive zero magnitudes yield
  zero; negative zero denominators remain under reciprocal execution. Carry
  boolean signs only in the group calculation; resolve them into fresh
  scalar/tensor Atom payloads, exact nested paths, or created/reused exact-scope
  Neg operations.
- Validate that canonical Power automorphisms and reconstruction descriptors map
  whole rooted copies, preserve the exact unsigned-magnitude even/odd
  boundary-pairing plan and signed Power polarity, and differ only by a common
  base plus scalar copy phases before folding.
- Canonicalize the exact normalized decoration orbit by the lexicographic
  `(ClassKey, RealizationKey)` contract on canonical line handles. Compute the
  class stabilizer using only primary executed-semantic `ClassKey` equality,
  then evaluate the resulting sign character on strong generators. Return
  `AmbiguousSignScope` only for a genuine semantic-class inequivalence.
  Symbolica, not the group layer, owns coefficient arithmetic and ordinary
  Sum/Product normalization.
- Execute and reparse every reconstructed iteration. Skip reprojection only
  when the complete one-sided stability certificate succeeds or execution
  returns the exact normalized iteration input; otherwise retry the reparsed
  whole network, with terminal zero/one as the remaining immediate return.
  Treat immediate equality of exact canonical problem identities as successful
  conditional verification, an older non-consecutive equality as
  `ConvergenceCycle`, and exhaustion of the injectable deterministic iteration
  budget as `IterationLimit`, without an unconditional Graphica verification
  pass.

### 5. Migrate integrations

- Route `IndexTooling`, metric, and color canonicalization through the
  network-owned path by constructing `CanonicalPolicyNet` only from each
  normalized source Atom, preserving color fixed-point re-entry. Remove the raw
  `SymbolicNetExt::canonize` entry point rather than emulating provenance with a
  preflight execution.
- Apply the same validation to Rust, dynamic-symbol, and Python ingestion; do
  not manufacture `Linear` or restrict ordinary untagged Symbolica heads.
- Remove `remove_antisymmetric_zero_terms`, duplicate parsing/expansion, all
  production `canonize_tensors` calls, and caller-level duplicates of the one
  reversible cooking/uncooking boundary retained inside `IndexTooling`.

### 6. Validate cost and downstream physics

- Benchmark graph size, fixed-point iteration count, and Graphica search on
  large color products, symmetric sums, repeated branches, and powers without
  implicit distribution.
- Use those iteration distributions to select and document the production
  convergence budget before migrating callers; benchmark thresholds are an
  acceptance decision rather than an unspecified telemetry result.
- Benchmark predicted versus actual expanded Power vertices/edges and Graphica
  search time, then select and document the production whole-graph resource
  limits before migrating callers.
- Run the three-loop RQFT signed six-`f` comparison and the disconnected scalar
  spectacles pointwise and Monte Carlo factorization checks.

The production budgets were measured on 2026-08-05 in jj change `rxtulqsp` on
an Intel Xeon W-2135 (restricted single-core Symbolica, Rust `test` profile).
The ignored `phase6_canonicalization_budget_report` uses the same exact
preflight estimator as production, asserts predicted and allocated graph sizes
agree, times Graphica directly, and then measures the complete driver. Timings
are diagnostic rather than test thresholds:

| case | vertices | edges | Graphica ms | iterations | driver ms |
|---|---:|---:|---:|---:|---:|
| signed K3,3 | 35 | 43 | 6.042 | 1 | 27.855 |
| three disconnected signed K3,3 components | 101 | 127 | 198.362 | 1 | 538.713 |
| repeated symmetric Sum branches | 42 | 49 | 2.765 | 2 | 57.269 |
| `Power(4)` | 23 | 26 | 1.489 | 1 | 20.761 |
| `Power(-4)` | 23 | 26 | 1.470 | 1 | 50.165 |
| `Power(24)` | 123 | 146 | 647.347 | 1 | 1,352.517 |

The observed maximum is two whole-graph iterations. The production limit of
eight leaves fourfold iteration headroom, including future cascading-zero
fixtures, while still giving a deterministic failure bound. The production
whole-graph limits are 128 vertices and 160 edges. The disconnected K3,3
product and `Power(24)` deliberately exercise those limits closely with two
different high-symmetry families; the exact preflight rejects a larger graph
before allocating any Graphica nodes. A fast non-ignored guard keeps the corpus
below all three production limits; the ignored report can be rerun with:

```text
nix develop -c cargo nextest run -p idenso --profile test_gammaloop \
  --run-ignored all \
  -E 'test(/^tensor::canonicalize::projection::tests::phase6_canonicalization_budget_report$/)' \
  --no-capture
```

The Power scaling probe was deliberately incremental. `Power(32)` already
needed 163 vertices, 194 edges, 1,583.129 ms in Graphica, and 3,368.133 ms in
the complete driver. A minimal symmetric two-port `Power(i8::MIN)` has an exact
643-vertex/770-edge projection; its direct sufficient-injected-budget Graphica
run was stopped without completion after 209 seconds. The signed `i8` grammar
therefore remains supported, but this pathological literal expansion returns
`GraphSizeLimit` under the production budget. A non-ignored regression proves
the unsigned magnitude is 128, constructs the exact projection under an
unbounded test budget, and checks rejection on both budget axes before Graphica
allocation. The ignored
`minimum_signed_power_projects_with_a_sufficient_injected_budget` test retains
the full direct probe for future Graphica improvements. The 128/160 limits are
thus a measured conservative acceptance envelope, not a claim that graph size
alone predicts canonical-search time.

The larger downstream checks remain separate from the timing corpus: the K3,3
checks pass, and Idenso's
`rqft_ghost_three_loop_d2_six_f_terms_vanish_independently` regression contains
the exact two structured-index six-`f` contractions extracted from d2 and
asserts that each vanishes independently. The
`scalar_spectacles_integrated_uv_factorizes_over_bridge` passes both its
pointwise and independent Monte Carlo comparisons. The enclosing historical
three-loop test currently proceeds through d10 before an unrelated d11 RQFT
snapshot mismatch. On the current debug test runner that ignored three-loop
case requires `RUST_MIN_STACK=67108864`; the default test-thread stack is
exhausted before completion, while the enlarged-stack run reaches the same d11
mismatch in 431 seconds.

## Test Plan

### Differential tests without antisymmetry

Use bare Symbolica as a test-only semantic reference with the same dummy
factory. Require exact Atom equality only for isolated tensor-local encodings
whose canonical-order contract is intentionally shared:

- ordered functions with scalar arguments before, between, and after slots;
- intrinsically symmetric and cyclesymmetric tensor heads containing direct
  slots only;
- untagged Symbolica symmetric and cyclesymmetric normalization controls over
  complete mixed immediate-argument lists, retained as exact non-network
  payloads rather than projected argument graphs;
- ordered tensors with parameters and nested structural `sym(...)` or
  `cyclic(...)` slot groups;
- `F(p, cyclic(a, b))` parsing as one `LocalTensor` with flat structure
  `[a, b]`, one tensor header, two port vertices, their line vertices, and no
  argument/group graph vertices;
- several partial groups in one tensor remaining distinct through stable layout
  paths;
- user-declared linear tensors with both bracketed and unbracketed composite
  arguments;
- free indices, multiple dummy groups, and dual representations;
- exact invariance under internal dummy relabeling with external slots fixed;
- equivariance, rather than exact equality, when the caller correspondingly
  renames free external indices.

For products, nested functions, sums, and parser-supported signed integer
powers, require that all isomorphic inputs produce the same Spenso canonical
output and that a test-only Symbolica canonicalization proves semantic
equivalence. Do not require two different canonical-labeling algorithms to
choose the same factor or dummy order unless that ordering is separately
specified as part of the API.

### Ordinary executor and leaf-representation tests

- represent the same tensor value as `LibraryKey`, `LocalTensor`, `TensorSum`,
  `ScaledTensor`, and `ScaledTensorSum`; require exact equal results from a
  Function, every supported signed Power case, and self-loop tracing through the
  shared eager materializer;
- instrument `Power(0)` and `Power(1)` to prove they return one or the original
  leaf before eager materialization, including library and lazy-sum leaves;
- permute mixed scalar Sum children across `Scalar`, scalar-valued local tensor,
  tensor sum, scaled tensor, and scaled-tensor sum representations; every
  permutation must produce the same result and error behavior;
- pin the shared tensor-sum materializer's explicit balanced and linear
  reduction policies on symbolic and floating-point fixtures, together with the
  existing lazy-sum threshold, so consolidation causes no accidental
  associativity change;
- feed identical Product operands through normal and in-place contraction and
  require the shared scalar analysis to identify the same accumulated scalar,
  scalar positions, and tensor positions before their distinct rewrites;
- append networks with already nonempty scalar and tensor stores and exercise
  every leaf variant, binary and n-ary Sum/Product, parallel replacement, and
  scalar convenience operations; all scalar and tensor references are rebased
  exactly once and still address the original values;
- specifically regress `network + i8` with a pre-existing scalar store entry,
  proving the appended scalar reference is shifted while tensor references are
  not shifted twice;
- verify that Neg preserves lazy sums/scales, Sum retains its lazy/eager choice,
  and Product retains scaled-term contraction rather than passing any of them
  through the eager unary materializer.

### Unified graph and globality tests

- instrument the canonicalizer to assert one root Graphica call per
  whole-network iteration and no per-Sum, per-Function, or per-Product calls;
- `2*A(i)` projects as an explicit Product with visible scalar `2` and tensor
  child `A(i)`; no input coefficient is stored in group sign metadata;
- an unsigned canonical-parser-produced input requiring only relabeling and
  ordering uses exactly one Graphica call because its complete stability
  certificate succeeds;
- `2*A_antisym(b,a)` first derives the whole-tensor sign, then sinks it within
  the Product to the existing scalar leaf as `2 -> -2`; ordinary execution
  produces `-2*A_antisym(a,b)` without hiding that scalar before Graphica, and
  the visible sign rewrite gives this fixture exactly two Graphica calls before
  stability;
- execution and normalized reparse run after every reconstruction, while only
  a complete one-sided stability certificate or exact executed-Atom equality
  may skip another Graphica call;
- association, commutative ordering, declared even/unsigned slot transport, and
  a bijective dummy rename preserve every visible equality class and certify
  stability under the selected vertex map;
- canonical transport that merges or splits contracted-index, leaf, Product,
  or Sum child syntax classes declines to certify stability even when it
  introduces no sign;
- two globally different line wirings with the same coarse node/color/degree
  statistics cannot be declared stable by a hash or multiset fingerprint;
- an incomplete or ambiguous stability witness conservatively retries instead
  of returning an early result;
- instrumentation proves there is no full-network execution before the first
  Graphica call and distinguishes the one full-network execution per rebuilt
  iteration from temporary scope executions used by exact group proofs;
- when direct leaf negation and execution move the sign into an existing scalar,
  for example `2 * (-A) -> -2 * A`, the changed semantic color retries the
  complete graph; a merge of previously distinct scalar/term colors is also a
  genuine coarsening invalidation;
- a locally zero addend whose removal exposes a new global odd automorphism
  invalidates the stability certificate and requires another whole-network
  iteration;
- terminal zero/one results return without a redundant verification call;
- a reconstruction for which the stability certificate conservatively declined,
  but whose next Graphica result has the same exact canonical problem identity,
  terminates successfully as a conditional verification without another
  reconstruction/execution;
- synthetic two- and three-state cycles return `ConvergenceCycle` with the
  correct first/repeated iterations and cycle length;
- a low injected budget over distinct synthetic states returns
  `IterationLimit`, while a legitimate cascading zero-removal fixture completes
  below the production budget;
- a synthetic problem-identity hash collision is rejected by full canonical
  graph/key equality;
- convergence diagnostics retain the retry reason that prevented stability
  certification;
- associative Product and Sum tree shapes produce identical graphs and output;
- commutative factor and addend order does not affect the result;
- identical local tensor fragments with different global line wiring remain
  distinguishable;
- dummy names reused independently in Sum branches do not merge line vertices;
- shared operator-boundary lines remain shared and receive one canonical dummy;
- a retry that removes an earlier representation group rebuilds the surviving
  line-to-position assignment while reusing the request-local callback value at
  that position;
- affine sign coordinates below proven-zero scopes are quotiented before orbit
  enumeration, while a mask that is not invariant under site transport is
  rejected rather than silently changing the group action;
- a root operator and an otherwise identical nested operator cannot exchange;
- an automorphism spanning an outer antisymmetric tensor and a Sum branch
  permutation is visible without distributing the Product;
- nonlinear Function boundaries retain local signs even when the global graph
  permutes equivalent arguments or surrounding factors;
- materialized Power copies expose cross-copy symmetries, validate as one base,
  and execute through the corrected ordinary Network Power operation;
- `Power(m)` and `Power(-m)` have the same complete-copy incidence but distinct
  signed Power-scope colors and exact canonical-problem identities;
- every materialized Power has one typed `PowerMagnitude` root with exactly
  `m` interchangeable `PowerCopy` children; a denominator-zero proof targets
  that root rather than conflating it with the reciprocal result node;
- every materialized occurrence has one identically colored `PowerCopy` root,
  and each canonical automorphism maps the complete descendant set of one copy
  root to exactly one complete target copy;
- copied Power bases clone internal lines and reproduce the existing Power
  interface pairing without creating a multiway line;
- supported self-dual exponents `-4`, `-2`, and two through five use one
  consistent complete-copy partner across every boundary slot: even magnitudes
  pair all copies, while positive odd magnitudes leave exactly one complete copy
  on the original external interface;
- pair-member swaps and complete pair-of-pairs swaps are visible, while a
  synthetic per-slot partner mismatch is rejected;
- equivalent hidden initial copy/pair ordinals produce the same canonical graph;
- the expanded whole-graph vertex/edge estimate is checked before allocation,
  agrees with constructed small fixtures, and returns `GraphSizeLimit` under a
  small injected budget rather than projecting an opaque Power. The same check
  uses unsigned magnitude for negative Powers, including 128 copies for
  `i8::MIN`;
- executing the reconstructed clone may contract Products into opaque leaves,
  but the retained canonical network still contains its separate tensor leaves
  and no executed leaf enters a Graphica projection;
- `A*(B+C)` remains factored, keeps one unified graph, and is never distributed
  merely to expose termwise tensor zeros.

### Signed unit tests

- odd and even antisymmetric cycles;
- repeated antisymmetric arguments;
- intrinsically antisymmetric tensor heads containing direct slots only;
- untagged Symbolica antisymmetric and cyclic normalization controls with mixed
  immediate arguments, retained exactly and never reinterpreted as slot-only
  network symmetry;
- ordered tensors with nested structural `antisym(...)` groups;
- negative stabilizers spanning several antisymmetric tensor factors;
- even stabilizers that must survive;
- two odd tensor contributions forming one even Product action;
- a returned generator moving a sign site rather than fixing it;
- an odd local stabilizer obtainable only by composing generators;
- two generated elements with the same signed-site projection but different
  full boundary permutations are distinguished by the pointwise stabilizer;
- the odd port exchange of `A_antisym(i,j)` with free `i` and `j` expresses
  antisymmetry without zeroing the tensor, while the corresponding repeated or
  pointwise-fixed-boundary case is zero;
- `A_antisym(i,j) * (B(i)*C(j) + B(j)*C(i))` is not zeroed at the inner tensor
  scope because its exposed lines move, but is zeroed at the complete outer
  Product scope;
- a composed Schreier generator retains both its complete canonical-vertex
  permutation and the correctly composed signed-site parity;
- exchanging two canonical antisymmetry sites whose input occurrences used
  opposite member orders has the same canonical generator parity in either
  input; the occurrence-order phase appears only in selected-map transport;
- direct frame-relative parity of a composed three-site map equals the signed
  cocycle composition, and inverse construction gives the inverse member
  permutation with the same parity;
- changing hidden occurrence origins or original member labels without changing
  the input frames' mathematical expressions does not change the canonical
  signed action;
- repeated group slots connected to one line remain distinct canonical-frame
  members;
- cyclic reconstruction starts at the smallest canonical port and follows the
  directed cycle, independently of original member rotation;
- signed Sum-branch exchange and cancellation;
- a three-term negative signed cycle and a signed term relation that appears
  only after generator composition;
- Sum combination and cancellation are produced by ordinary network execution
  plus Symbolica normalization, with no coefficient-free body-key reducer;
- `2*A + 3*A`, `A - A`, and symbolic-coefficient variants match Symbolica's
  actual normal form; a visible merge or cancellation in these cases triggers
  retry under the general visible-invalidation rule;
- `2*A_antisym(b,a) + 2*A_antisym(a,b)` executes to zero and covers visible
  scalars, scalar-leaf sign sinking, cancellation, conditional retry, and
  terminal return in one fixture;
- independence from generator order and generating-set choice in direct tests
  of the signed-action layer;
- breadth-first versus depth-first orbit traversal and reordered hash insertion
  choose the same lexicographic `(ClassKey, RealizationKey)` minimum rather than
  the first representative encountered;
- Product decorations `(-A)*B` and `A*(-B)` have one executed-semantic
  `ClassKey`, while deterministic Product sign sinking gives one concrete
  `RealizationKey`; the class stabilizer retains movements between those
  equivalent placements;
- nonlinear `f(-x)` versus `-f(x)` and `-B^2` versus `(-B)^2` have different
  semantic class keys;
- alpha-equivalent scope fragments with different concrete dummy names, Symbol
  interner histories, hidden origins, and caller dummy values have identical
  class keys, while a synthetic hash collision between unequal full keys does
  not compare equal;
- for small generated groups, exhaustive orbit enumeration and the optimized
  stabilizer-chain minimum select the same semantic class and realization;
- a nested antisymmetric sign remaining inside a nonlinear outer tensor;
- the same sign lifting through an explicitly `Linear` outer tensor, while the
  same undeclared head remains nonlinear by default;
- `spenso::bracket(p + q)` remaining atomic under a linear outer tensor while
  an unbracketed `p + q` follows ordinary linear distribution;
- a zero nested antisymmetric group remaining an argument of a nonlinear outer
  tensor and zeroing an explicitly linear outer tensor;
- zero summands inside ordinary nonlinear functions:
  `f(x + zero) == f(x)`;
- external indices and incompatible representation groups blocking invalid
  automorphisms;
- closed, open, and nested signed-Power cases, including executed results for
  `(-A)^2`, `(-A)^3`, scalar `(-s)^-3`, self-dual `(-A)^-2`, and positive and
  negative proven-zero magnitude scopes;
- exact-result generic execution of scalar, `LibraryKey`, and every tensor-like
  leaf kind for its supported cases across `-5..=5`, guarding the
  empty-scalar-store zero-Power case, every duplicated odd-magnitude regression,
  and inversion exactly once after a scalar magnitude result;
- `Power(0)` of scalar, self-dual, and non-self-dual bases takes the parser,
  `Network::pow`, and `NetworkState::pow` early-one path without constructing a
  Power graph, inspecting dangling tensor state, or requiring an existing
  scalar store entry;
- a low-level explicitly built scalar `Power(i8::MIN)` executes with magnitude
  128 without overflow. Canonical parsing of pure scalar `x^-128` instead keeps
  one exact opaque scalar Atom leaf and creates no Power copies. A materialized
  tensor-derived scalar or self-dual `Power(i8::MIN)` exercises the even
  128-copy contraction before inversion when it fits the graph budget;
- a negative odd open self-dual tensor remains a typed parser error, while
  negative odd closed scalar tensor networks and all scalar leaves remain valid;
- reconstructed closed Products recompute `NetworkState::Scalar` from their
  graph before a negative odd Power; a stale builder-side
  `NetworkState::SelfDualTensor` cannot change admissibility or output state;
- a synthetic even Power whose canonical materialized copies differ by one
  residual scalar phase executes as `-B^2`, not `(-B)^2`, while two residual
  phases cancel;
- the reciprocal counterpart executes as `-B^-2`, not `(-B)^-2`; uniform base
  signs satisfy `(-B)^-2 = B^-2` and scalar `(-B)^-3 = -B^-3`;
- a stabilizer-proven zero denominator reconstructs as a scalar zero under an
  ordinary reciprocal, declines stability, and lets Symbolica normalize the
  singular value (`0^-1` becomes Symbolica's complex infinity); neither it nor
  an enclosing nonlinear Function is replaced by zero;
- whole-network execution, without an outer-Product zero shortcut, owns cases
  such as `0 * (zero denominator)^-1`, reparses Symbolica's indeterminate result,
  and exposes that normalized form to the next complete graph iteration;
- Power folding rejects partial-copy, boundary-pairing, localized-sign, and
  non-scalar base mismatches instead of choosing an arbitrary regrouping;
- `(A*B)^2` follows the executor's Symbolica-normalized topology, retries when
  that topology prevents a stability certificate, and is then idempotent;
- `(A*B)^-2` likewise follows Symbolica's integer-Power-of-Product
  normalization, preserves division polarity, retries on any visible topology
  rewrite, and is then idempotent;
- an output containing a direct leaf sign or exact-scope Neg remains idempotent
  after execution and reparsing;
- direct scalar-leaf, whole-tensor-leaf, and nested partial-group signs are
  copy-on-write and do not mutate another occurrence sharing the original store
  entry;
- a Product with multiple eligible scalar/tensor leaves sinks its sign by
  `(scalar before tensor, canonical vertex ID)` across Product-only descendants;
  a Product containing only Sum, Function, or Power children instead receives
  an exact-scope Neg;
- whole-Sum, nonlinear Function-result, and residual even-Power signs remain at
  their exact Neg scopes until execution; tests distinguish `-(x+y)` from one
  negated addend, `-f(x)` from `f(-x)`, and `-B^2` from `(-B)^2`;
- induced signs never synthesize a scalar `-1` Product factor. Pin how the
  canonical parser represents an explicit source-level `-1` factor and require
  that exact parser-produced scalar representation, rather than unnormalized
  input syntax, to remain visible and be preserved;
- two selected canonical maps differing by a stabilizer execute to the same
  signed representative.

For small ranks and networks, add a brute-force reference that enumerates
intrinsic and partial tensor symmetries, identical-factor permutations,
dummy-line relabelings, and signed Sum actions in addition to declared Function
semantics. Graphica's complete canonical graph fixes the unsigned vertex
numbering first; the reference exhaustively checks the signed transformation
closure and chooses the minimal scoped decoration within that fixed unsigned
class, rather than introducing a competing raw-Atom factor/dummy ordering. It
detects whether both signs reach one representative. Explicitly limit each
generated case to the transformations covered by that oracle.

### Structured-index and failure tests

- reject scalar, nested-function, `aind(...)`, and partial-group immediate
  arguments on intrinsically symmetric tensor heads;
- accept the parameterized counterpart as an ordered tensor containing an
  explicit partial-slot group;
- require every child of a structural partial-slot group to be a direct slot,
  without changing the established projector use of the same heads;
- preserve unrestricted mixed arguments on ordinary untagged Symbolica symmetry
  heads without projecting them as tensor topology;
- confirm that neither tensor construction nor canonicalization adds `Linear`;
- define slots hidden inside `spenso::bracket(...)` as outside the enclosing
  tensor's direct structure and reject an occurrence that would otherwise rely
  on them as structural slots;
- structurally distinct indices that previously shared a flattened name;
- nested functions, namespaces, underscores, and zero-argument functions;
- index payloads containing sums, products, powers, and supported numeric
  forms;
- invalid or unsupported inputs returning the documented result rather than
  panicking;
- negative scalar-result Powers are accepted as division. A pure scalar
  `x^-128` remains one opaque scalar Atom leaf; a materialized tensor-derived
  scalar or self-dual `Power(i8::MIN)` either creates exactly 128 copies under a
  sufficient injected budget or returns `GraphSizeLimit` before allocation
  under a smaller one. A separate low-level scalar-Power executor test covers
  magnitude 128; neither path calls signed `abs()` or overflows;
- negative odd open self-dual Powers and every non-self-dual tensor Power return
  their existing typed parser errors, while sign alone is never an
  `UnsupportedPower` reason;
- a fractional, symbolic, or out-of-`i8` exponent may remain an opaque scalar
  Atom only when validation proves that its base contains no recognized tensor
  topology; otherwise canonical parsing rejects it before Graphica rather than
  hiding tensor ports in a scalar leaf;
- the canonical parser returns `CanonicalPolicyNet` with the exact normalized
  source Atom, while a parser using depth limits, shorthand/opaque boundaries,
  or other noncanonical settings returns only a raw `SymbolicNet` and has no
  canonicalization entry point;
- the proof wrapper has no raw-network conversion, unchecked constructor,
  mutable network projection, or other API by which the source Atom and network
  can become inconsistent;
- instrument the first iteration to prove that Graphica receives the network
  directly from `CanonicalPolicyNet` without a pre-Graphica ordinary execution;
- injected failures from post-reconstruction full-network execution, temporary
  scope execution, and result extraction return
  `CanonicalizationError::Execution` with the relevant whole-network or scope
  context rather than panicking;
- a non-consecutive exact canonical-state repetition returns
  `CanonicalizationError::ConvergenceCycle`, while exhaustion of the
  deterministic injectable iteration budget returns
  `CanonicalizationError::IterationLimit`;
- a Power whose predicted literal expansion exceeds an injected whole-graph
  budget returns `CanonicalizationError::GraphSizeLimit` before allocation,
  without falling back to an opaque Power projection; prediction uses the
  unsigned magnitude and reports both `n` and `|n|`;
- repeated calls remaining independent of Symbolica's global symbol-interner
  history;
- scalar and opaque colors containing index-like lexical symbols either remain
  validated atomic data or fail if they own a recognized contracted slot;
- a dependency-contract test pins Graphica's input-to-canonical `vertex_map`
  direction and canonical-numbered complete-generator cycle representation,
  including fixed points when materialized over the full vertex domain.

### Invariants and integration tests

- canonicalizing the final proof-carrying wrapper again is idempotent in both
  its network and Atom projections, and the Atom facade is independently
  idempotent;
- the network and Atom projections of the final `CanonicalPolicyNet` are the
  same result: executing the network projection produces exactly the retained
  Atom used by `IndexTooling::canonize`, modulo the facade's reversible
  cook/uncook boundary;
- dummy relabeling and commutative factor/addend reordering do not change the
  result;
- tensor-header colors contain slot-erased layouts and never retain contracted
  dummy identities;
- hidden origins, hedge IDs, layout member positions, and union-find
  representatives do not change visible graph identity;
- stable semantic colors, not raw Symbol interner IDs, determine graph order;
- equivalent inputs receive the same canonical sign;
- valid intrinsic symmetry, generic mixed-argument symmetry, and nested partial
  symmetry each retain their distinct head-local semantics;
- color K3,3 and the two exact RQFT d2 six-`f` contractions remain zero;
- normal-profile color, metric, Dirac, and feyngen canonicalization tests pass;
- `scalar_spectacles_integrated_uv_factorizes_over_bridge` validates
  pointwise `bubble_left * bubble_right / (2s)` and the Monte Carlo spectacles
  integral against the independently sampled `bubble^2 / (2s)` target within
  the configured statistical criterion.

## Acceptance Criteria

- Production idenso code no longer calls bare Symbolica `canonize_tensors`.
- Every iteration projects the complete normalized operation tree, recognized
  tensor ports, and index lines into one graph and calls Graphica once. The
  `ContextGraph`/flat-scope split and barriers are gone; neither local Graphica
  calls, opaque pre-labeling contractions, nor a custom graph-building executor
  remain. Opaque contraction is allowed only inside temporary post-Graphica
  group-proof scope clones and the post-reconstruction full-network execution
  clone; every such result is reparsed or reduced to an Atom comparison before
  further use.
- Every reconstruction has one complete-network execution through ordinary
  symbolic network execution and is reparsed. Fallible temporary scope
  executions are allowed only for exact group quotient proofs. Skip another
  Graphica iteration only when reconstruction certifies that the normalized
  input changed solely by association, commutative ordering, declared unsigned
  symmetry transport, and bijective dummy renaming while preserving visible
  topology, colors, scopes, Power structure, and all bottom-up Symbolica-visible
  equality classes, or when mandatory execution returns the exact normalized
  iteration input and thereby proves an extensional fixed point. Any sign, zero,
  payload edit, residual phase, class merge/split, or ambiguity otherwise
  retries the complete reparsed problem; terminal zero/one returns immediately.
  Consecutive equality of exact canonical problem
  identities is successful conditional verification; an older repetition is a
  typed convergence cycle, and a benchmarked deterministic iteration budget
  bounds distinct states. Hash equality alone proves neither case. Ordinary
  canonical transport has no unconditional final Graphica verification pass.
  Network and Atom idempotence remain required.
- Initial Atom parsing and every post-execution reparse use the same fixed,
  fully recursive canonical parse policy and tensor-symmetry validator. The
  opaque parser result pairs the exact normalized source Atom with its network,
  keeps both fields consistent, and is the canonicalizer's only input type. A
  raw `SymbolicNet`, whether manually built or parsed with noncanonical
  settings, has no canonicalization entry point or unchecked conversion into
  that wrapper. A non-`i8`, fractional, or symbolic exponent may be retained as
  an opaque scalar Atom only when its base has no recognized tensor topology or
  an explicit bracket/projector opacity boundary owns the complete payload; it
  cannot otherwise hide tensor ports from the unified graph.
- Root and operation roles, unary `NetworkOp::Function`, signed Power-result
  nodes, `PowerMagnitude`/`PowerCopy` nodes, line connectivity, and semantic
  colors are explicit and stable across parser association, dummy names, hidden
  origins, and Symbol interner history. A Power's exact signed exponent is part
  of its visible color and canonical problem identity; magnitude topology,
  parity, and pairing use only its unsigned magnitude.
- Existing scalar expressions remain visible leaves and Product children. The
  signed group layer carries only mathematical boolean decorations. It stores no
  scalar/sign payload in network operations: reconstruction resolves each sign
  into a deterministic fresh scalar/tensor Atom value within the same
  sign-transparent Product, an exact nested Atom path, or a created/reused
  exact-scope Neg before execution, never a synthesized scalar `-1` Product
  factor.
- Each `LocalTensor` uses one reversible slot-erased header plus typed ports and
  cyclic markers. Intrinsic direct-slot symmetry and explicit nested partial
  groups match the declared Symbolica overlap; invalid intrinsic arguments fail
  cleanly and untagged opaque heads are not reinterpreted.
- Network connectivity, not printed index equality, defines lines.
  Reconstruction restores operation boundaries before applying shared-interface
  and branch-local dummy namespaces, preserving external `LibrarySlot`s and
  structured indices.
- Selected-map phases are placed at their exact intrinsic or partial-group
  sinks. Local zero requires scalar `-identity` on the complete reduced scope
  under its pointwise boundary stabilizer, including generator compositions.
  Every canonical symmetry site has a local port frame derived only from
  canonical vertex IDs and, for cyclic groups, the visible directed cycle.
  Selected-map parity compares the original member frame with that canonical
  frame; generator parity compares canonical source and target frames.
  Every original and derived group element retains its complete permutation of
  canonical graph vertices together with its signed-site action; scope roots
  and exposed boundaries are fixed pointwise using the vertex permutation
  before projecting to the semantic child/site action.
- Exact Product and Sum child actions are used only for signed group proofs,
  map independence, and outward propagation. There is no custom coefficient
  extraction, coefficient-free Sum body key, or coefficient aggregation;
  ordinary symbolic network execution owns multiplication, addition,
  cancellation, and declared linearity. Nonlinear Function scopes retain local
  signs and zeros.
- Alternative maps produce the same exact normalized decoration-orbit
  representative on canonical line handles after innermost zero substitution.
  The representative is the explicit lexicographic minimum of a primary stable
  executed-semantic `ClassKey` and secondary canonical-sink `RealizationKey`,
  independent of generator traversal, dummy names, hidden origins, hashes, and
  Symbol interner IDs. The class stabilizer uses primary semantic equality,
  allowing movement among equivalent concrete sign placements, before its
  one-dimensional sign character is tested on strong generators. Power
  magnitude scopes validate their materialized copy-line semantics and retain
  the XOR of copy-relative scalar phases in a Neg outside the signed Power.
  Proven-zero positive magnitudes produce zero; proven-zero negative
  denominators remain under an ordinary reciprocal, decline stability, and are
  normalized by whole-network Symbolica execution rather than being replaced
  by zero. Each materialized occurrence has an identical visible copy root with
  enforced whole-copy membership; self-dual boundary lines follow the existing
  even/odd unsigned-magnitude complete-copy pairing, and folding validates
  topology, pairings, semantic base classes, and original exponent polarity.
  Generic Power execution uses unsigned magnitude, supports division and
  `i8::MIN`, and never attempts non-scalar tensor inversion. Literal expansion
  is preflighted with checked arithmetic against benchmarked whole-graph
  vertex/edge limits and fails with `GraphSizeLimit` rather than becoming
  opaque.
- `Linear` remains opt-in, brackets remain atomic, and canonicalization never
  distributes a Product over a Sum merely to find tensor zeros. It does accept
  ordinary non-distributive Symbolica normal forms, such as rewriting an
  integer Power of a Product, and reparses them before declaring stability.
- Equivalent `NetworkLeaf` storage forms have identical ordinary-executor
  semantics. `Function`, nontrivial Power, and self-loop tracing share one
  fallible eager tensor-leaf materializer, while `Power(0)` and `Power(1)`
  bypass it and Neg, Sum, and Product retain their lazy behavior. Mixed scalar
  Sum results are independent of commutative child order.
- Tensor and scalar reference traversal and network-store append/rebase use
  shared paths that rebase each reference exactly once. The duplicate tensor-sum
  materializers share one core with an explicit reduction-order policy, and
  normal and in-place Product contraction share operand classification without
  conflating their distinct graph rewrites.
- Only temporary post-Graphica scope clones used by exact group proofs and one
  post-reconstruction full-network clone per iteration are passed to ordinary
  `Network::execute`. The initial parser-produced network is projected without
  execution. Every executed value is reparsed or reduced to an Atom comparison;
  no opaque executed leaf enters Graphica. The normalized full-network reparse
  and its exact source Atom form the stable `CanonicalPolicyNet`; existing
  facades expose network or Atom projections. Canonicalization contains no
  bespoke structural renderer or Power evaluator.
- Targeted unit and integration checks, the K3,3 checks, Idenso's independent
  assertions over the two exact d2 signed six-`f` RQFT contractions, and the
  disconnected Monte Carlo spectacles validation pass. The enclosing
  historical RQFT test still has the unrelated d11 snapshot mismatch recorded
  in Phase 6; full downstream validation is not yet complete.
