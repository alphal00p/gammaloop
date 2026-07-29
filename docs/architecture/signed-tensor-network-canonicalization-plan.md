# Signed Tensor-Network Canonicalization Plan

## Objective

Replace the current two-stage canonicalization path

```text
signed zero-detection graph -> bare Symbolica tensor canonicalization
```

with one network-owned canonicalization path:

```text
Atom
  -> SymbolicNet
  -> flat slot-symmetry graph
  -> Graphica canonical form and stabilizer generators
  -> signed canonical SymbolicNet
  -> network execution
  -> canonical Atom
```

The graph must encode function symmetries with the same head-local semantics as
Symbolica: an attribute acts on every immediate argument of the head that owns
it. Spenso adds a stricter grammar for intrinsically symmetric tensor heads so
that those immediate arguments are exactly the tensor slots. Canonicalization
must produce both the canonical network labeling and its sign, including zeros
caused by negative stabilizers. Bare
`symbolica::AtomCore::canonize_tensors` must no longer be part of the production
path.

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
- Signed reasoning is local to the subtree being canonicalized. A negative
  automorphism zeros that subtree: a zero summand is removed, while a zero
  function argument does not erase a nonlinear enclosing function.
- Canonicalization preserves factorization. It must not implicitly distribute
  products of sums merely to discover additional zero monomials.
- Applying canonicalization twice produces the same network and Atom.

## Target Ownership and API Flow

The canonicalizer belongs to `SymbolicNet`, exposed through the existing
`SymbolicNetExt` boundary in `crates/idenso/src/tensor/mod.rs`. Its owning method
should consume a network and return a canonical network, using the caller's
dummy-index factory:

```rust
network.canonize(new_dummy) -> Result<SymbolicNet<Aind>, CanonicalizationError>
```

`IndexTooling::canonize` remains the Atom-facing facade, but its implementation
becomes:

1. parse the Atom into one `SymbolicNet`;
2. validate the head-local tensor symmetry grammar;
3. canonize that network;
4. call `simple_execute()` on the canonical network;
5. restore any reversibly encoded structured indices.

Validation must cover every Atom/network ingestion route, including dynamically
created tagged symbols and the Python-facing constructors. An intrinsically
symmetric tensor occurrence whose immediate arguments are not all direct slots
returns a typed error. Parsing must preserve the attributes the user declared
on each head and must never manufacture `Linear`.

The internal graph projection and reconstruction stay crate-private. Do not add
a second public canonicalization mode or retain an old-path fallback. Migrate
all in-tree callers to the network-owned behavior.

## Flat Graph Projection

Use Symbolica's tensor graph implementation as the behavioral specification,
not as a requirement to copy its Atom syntax tree. `SymbolicNet` and the induced
Graphica graph stay at Spenso's existing granularity:

- the network operation tree owns products, sums, powers, and enclosing
  functions;
- one `LocalTensor` leaf owns a flat `OrderedStructure`, its slot hedges, and the
  exact `SymbolicTensor.expression` Atom;
- one projected tensor occurrence has one tensor vertex and one vertex per flat
  structural slot.

The arguments inside a `LocalTensor` do not become network or Graphica nodes.
Instead, the tensor-vertex color carries a slot-erased argument layout, and
hidden origin data retains the information required to fill its holes again.
In particular, `p`, `cyclic(...)`, and `spenso::bracket(...)` remain in that
layout; they are not separate graph vertices.

### Tensor layout and flat slot graph

For example, parse

```text
F(p, cyclic(a, b))
```

as one `LocalTensor` with structure `[a, b]` and exact expression
`F(p, cyclic(a, b))`. Its flat Graphica projection is equivalent to:

```text
tensor vertex color: F(p, cyclic(□, □))
tensor -> a  [group = argument 1 / cyclic]
tensor -> b  [group = argument 1 / cyclic]
a -> b -> a  [cyclic marker]
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

and the flat incidence graph contains

```text
tensor -> a     [visible = Ordered(arg[1])]
tensor -> b     [visible = SymmetricGroup(arg[2]), hidden member = 0]
tensor -> c     [visible = SymmetricGroup(arg[2]), hidden member = 1]
tensor -> d     [visible = Ordered(arg[3])]
tensor -> e     [visible = AntisymmetricGroup(arg[4]), hidden member = 0]
tensor -> f     [visible = AntisymmetricGroup(arg[4]), hidden member = 1]
tensor -> g     [visible = CyclicGroup(arg[5]), hidden member = 0]
tensor -> h     [visible = CyclicGroup(arg[5]), hidden member = 1]
tensor -> i     [visible = CyclicGroup(arg[5]), hidden member = 2]
g -> h -> i -> g  [cyclic markers]
```

There is no ordering edge from `a` to `b` or `c`. Outer order is absolute:
the layout places `DirectSlot(arg[1])` before `SymmetricGroup(arg[2])`, and
those distinct visible incidence roles prevent any automorphism from exchanging
`a` with `b` or `c`. The two group members share the same visible `arg[2]`
role, so they may exchange. Hidden member positions are used only for
reconstruction and antisymmetric parity. Only the cyclic group adds
slot-to-slot edges.

Each tensor-to-slot incidence edge has a typed role:

- an ordered slot has a visible argument path or position;
- members of one symmetric, antisymmetric, or cyclic group share a visible,
  stable group key derived from the layout path;
- the original member position is hidden reconstruction data;
- distinct partial groups have distinct group keys, so their slots cannot mix;
- cyclic members additionally have directed cycle-marker edges;
- antisymmetric groups additionally contribute their induced permutation to
  the sign calculation.

Group keys must come from stable structural paths, not occurrence IDs. Symmetry
and parity belong to the group recorded on the flat incidence edges; they are
not inferred from all slots incident to the outer tensor vertex.

External slot vertices retain their complete `LibrarySlot`. Internal slot
vertices retain endpoint representation, while undirected contraction edges
are colored by normalized representation group.

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
5. create one tensor vertex plus its flat slot vertices, add the typed incidence
   and cyclic-marker edges, and connect contractions from the existing Spenso
   hedge pairs;
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

Sums and enclosing functions are canonicalized recursively at the network
operation level, outside the per-scope Graphica incidence graph. Positive powers
may materialize repeated flat tensor-and-slot copies. Unsupported powers return
a typed error or retain the original subtree according to the Atom-facing API
contract; they must not panic.

## Signed Canonicalization

For each additive term or recursively canonicalized multiplicative scope:

1. project the flat unsigned tensor-and-slot graph;
2. call `Graph::canonize()` once;
3. evaluate the induced permutation separately for every antisymmetric tensor
   or partial-group key on each returned automorphism generator;
4. if a generator is negative at a valid local sign scope, replace that local
   tensor/group subtree with zero, lifting to the whole term only through
   products or user-declared linear heads;
5. otherwise calculate the scoped signs of the input-to-canonical `vertex_map`;
6. rebuild the canonical network by filling the canonical slots into each
   tensor layout and attaching every sign at its recorded group path.

The signed output must be independent of which map to the canonical graph is
chosen. Two maps differ by a stabilizer; after local negative stabilizers have
been replaced by zero, the remaining stabilizers must induce the same sign at
every nonlinear scope. Signs may be combined only after crossing a product or a
user-declared linear boundary. If the retained metadata cannot establish that
scope, canonicalization must return a typed error instead of choosing an
arbitrary global sign.

Assign fresh dummy indices in canonical contraction-edge order, separated by
representation group. Preserve external indices exactly. Fill ordered layout
holes by their visible role, symmetric and antisymmetric group holes by
canonical member order, and cyclic group holes by following the directed cycle
from its canonical start. Rebuild the original wrapper at its stored argument
path; do not construct an intermediate graph subtree.

Attach the sign at the scope where it was produced. A sign arising inside a
function argument must remain inside that argument; it cannot be moved to the
global coefficient of a nonlinear enclosing function. If the user declared the
enclosing tensor head `Linear`, normal Symbolica semantics may lift that local
coefficient during reconstruction or execution. The canonicalizer must neither
infer linearity nor hoist a coefficient through a nonlinear head itself.

## Sums, Functions, and Zero Propagation

Canonicalize summands recursively and independently, remove zero summands, and
then rebuild the enclosing sum. This must work at any depth, giving

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

Do not interpret an odd automorphism confined to one addend as a negative
automorphism of the complete sum. Do not globally expand products containing
sums during canonicalization. Detecting zeros that exist only after
distribution remains the responsibility of an explicit expansion pass.

## Applying the Canonical Result to the Network

The canonical graph should produce a fresh canonical `SymbolicNet` rather than
round-tripping through an Atom or partially mutating the original network.
Rebuild or replace, in one coherent operation:

- tensor expressions by filling the reversible slot-erased layouts;
- intrinsic versus ordered head classification and all user-declared
  attributes;
- explicit partial-symmetry heads and bracket boundaries;
- flat tensor structures, their ordered slots, and `slot_order`;
- contraction edges and canonical dummy labels;
- existing network operations for products, sums, functions, and powers, plus
  local signs;
- external slots and representation orientation.

Use existing Spenso network constructors and builders where they already own
these invariants. Before adding any new constructor or relabel helper, confirm
that the existing network construction API cannot express the required atomic
operation.

Once the signed canonical network is complete, use `simple_execute()` to
materialize the final Atom. Symbolic contraction then combines already
canonical tensor expressions and preserves their local signs. Standard
Symbolica normalization may perform sign lifting authorized by a declared
`Linear` attribute; no canonicalizer-side lifting through nonlinear heads is
allowed. No subsequent bare tensor canonicalization is allowed.

## Implementation Phases

### 1. Capture the unsigned reference behavior

- Add a test-only reference helper that calls bare Symbolica
  `canonize_tensors` with deterministic dummy names and representation groups.
- Split the reference corpus between valid intrinsic tensor symmetry over
  direct slots and ordered tensors containing explicit partial-symmetry groups.
- Keep ordinary untagged Symbolica heads with arbitrary immediate arguments as
  normalization controls: their exact Atoms must remain unchanged by the
  network projection, but they are not expanded into tensor/argument graph
  nodes under the default tagged-tensor parser.
- Record representative ordered, symmetric, cyclesymmetric, nested-function,
  sum, power, external-index, and dual-representation cases, plus nonlinear,
  user-declared `Linear`, and bracketed normalization cases.
- Record typed rejection of mixed intrinsic tensor heads separately from the
  Symbolica overlap corpus.
- Confirm that relabeled or factor-permuted versions produce identical
  reference Atoms.

### 2. Replace the projection with the flat graph

- Replace `TensorColor { head, arguments }` with a slot-erased tensor layout that
  keeps non-slot data, ordered argument boundaries, partial-group wrappers, and
  brackets while excluding dummy identities.
- Keep the projection flat: one tensor header, one vertex per structural slot,
  typed incidence roles, and cyclic-marker edges. Do not add graph nodes for
  arguments or partial-group wrappers.
- Extend structure inference to expose direct children of structural
  `sym`/`antisym`/`cyclic` groups as flat slots, while retaining their layout
  paths for reconstruction and sign scope.
- Validate intrinsic tensor heads and reject slots hidden under unsupported
  nested heads or brackets.
- Use Spenso hedges for contractions and remove flattened cooking from graph
  identity.
- Keep this phase unsigned and compare its canonical ordering against the
  Symbolica reference suite.

### 3. Rebuild and execute a canonical network

- Convert the canonical graph directly into a fresh `SymbolicNet`.
- Fill each stored layout from the canonical flat slot order, then rebuild the
  tensor expression, `OrderedStructure`, `slot_order`, dummy assignment,
  external indices, and representation orientation together.
- Preserve intrinsic tensor heads, exact untagged Atom payloads, explicit
  partial groups, brackets, and declared attributes without adding `Linear`.
- Execute the resulting network with `simple_execute()`.
- Require exact Atom equality with bare Symbolica for every supported reference
  case that contains no antisymmetric function, covering valid intrinsic and
  nested partial-group forms. Check the untagged generic controls separately so
  the flat projection does not reinterpret them.

### 4. Add signed transport and stabilizer zeros

- Calculate parity from the flat incidence edges belonging to each
  antisymmetric tensor or structural partial-group key. Keep the recorded layout
  path as that parity's sign scope.
- Detect negative stabilizer generators and produce zero at the correct local
  scope.
- Calculate the scoped input-to-canonical transport signs when the remaining
  stabilizers are even at those scopes.
- Apply signs and zeros at their owning subtree and canonicalize sums
  recursively.
- Verify that nested signs stay local under nonlinear outer tensors, lift only
  under user-declared `Linear` heads, and retain bracketed composite arguments.

### 5. Migrate integrations and remove the old path

- Replace the prepass in `IndexTooling::canonize` with network
  canonicalization followed by execution.
- Enforce the same intrinsic-symmetry validation when ingesting expressions
  produced by Rust builders/macros, dynamically created tagged symbols, and
  Python-facing constructors. Do not restrict ordinary untagged Symbolica heads
  or change the default `Linear` setting.
- Route metric and color canonicalization through the same network-owned
  implementation.
- Preserve color fixed-point re-entry when canonicalization removes a color
  term and exposes more color rewrites.
- Remove `remove_antisymmetric_zero_terms`, its duplicate parsing/expansion,
  and all production calls to bare `canonize_tensors`.
- Centralize reversible structured-index handling at the Atom/network boundary;
  remove caller-level pre-cooking from the RQFT diagnostic.

### 6. Validate performance and downstream physics

- Confirm that canonicalization does not increase term count solely because a
  product contains sums.
- Benchmark representative large color networks and positive powers against a
  recorded pre-change baseline.
- Run the three-loop RQFT color comparison with its signed six-`f` zeros.
- Run the disconnected scalar spectacles Monte Carlo test, including both the
  pointwise factorization check and the independent integral comparison against
  `bubble^2 / (2s)`.

## Test Plan

### Differential tests without antisymmetry

Compare the new network path exactly against bare Symbolica using the same
dummy factory:

- ordered functions with scalar arguments before, between, and after slots;
- intrinsically symmetric and cyclesymmetric tensor heads containing direct
  slots only;
- untagged Symbolica symmetric and cyclesymmetric normalization controls over
  complete mixed immediate-argument lists, retained as exact non-network
  payloads rather than projected argument graphs;
- ordered tensors with parameters and nested structural `sym(...)` or
  `cyclic(...)` slot groups;
- `F(p, cyclic(a, b))` parsing as one `LocalTensor` with flat structure
  `[a, b]`, one tensor header, two slot vertices, and no argument/group graph
  vertices;
- several partial groups in one tensor remaining distinct through stable layout
  paths;
- user-declared linear tensors with both bracketed and unbracketed composite
  arguments;
- nested functions and opaque index-free subexpressions;
- products, nested sums, and positive powers;
- free indices, multiple dummy groups, and dual representations;
- deterministic relabelings and product-factor permutations.

The comparison is direct Atom equality, not string similarity or numerical
evaluation.

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
- closed, open, and nested positive-power cases.

For small ranks and networks, add a brute-force reference that enumerates the
declared function symmetries, chooses the minimal representative, and detects
whether both signs reach that representative.

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
- repeated calls remaining independent of Symbolica's global symbol-interner
  history.

### Invariants and integration tests

- canonicalization is idempotent at both network and Atom levels;
- dummy relabeling and commutative factor reordering do not change the result;
- tensor-header colors contain slot-erased layouts and never retain contracted
  dummy identities;
- equivalent inputs receive the same canonical sign;
- valid intrinsic symmetry, generic mixed-argument symmetry, and nested partial
  symmetry each retain their distinct head-local semantics;
- color K3,3 and RQFT six-`f` zeros remain zero;
- normal-profile color, metric, Dirac, and feyngen canonicalization tests pass;
- `scalar_spectacles_integrated_uv_factorizes_over_bridge` validates
  pointwise `bubble_left * bubble_right / (2s)` and the Monte Carlo spectacles
  integral against the independently sampled `bubble^2 / (2s)` target within
  the configured statistical criterion.

## Acceptance Criteria

- Production idenso code no longer calls bare Symbolica `canonize_tensors`.
- The network retains its existing operation boundaries, while each
  `LocalTensor` projects to one flat tensor header plus structural slot vertices.
- A reversible slot-erased tensor layout preserves non-slot arguments,
  direct-slot positions, partial-group paths, and bracket boundaries without
  introducing argument/group Graphica vertices or retaining dummy identities in
  tensor colors.
- Typed incidence roles and cyclic-marker edges reproduce the allowed ordered,
  symmetric, and cyclic slot permutations on the unsigned overlap suite.
- Valid intrinsic tensor heads and ordered tensors with explicit
  partial-symmetry groups produce exactly the same unsigned Atom as Symbolica;
  ordinary untagged generic heads remain exact normalized payloads and are not
  reinterpreted as network tensor topology.
- Intrinsically symmetric tensor heads with any non-direct-slot argument fail
  cleanly, including scalar parameters, nested groups, and `aind(...)` bundles.
- Partial tensor-slot symmetry is represented only by an explicit nested
  `sym(...)`, `antisym(...)`, or `cyclic(...)` group under an ordered tensor.
- Negative antisymmetric stabilizers produce zero; nonzero representatives have
  a deterministic canonical sign.
- A nested antisymmetric sign or zero remains local under a nonlinear outer
  tensor and lifts through the outer head only when the user declared it
  `Linear`; `spenso::bracket(...)` remains the atomic-argument escape hatch.
- No construction or canonicalization path adds `Linear` implicitly, and
  ordinary Symbolica heads retain unrestricted mixed-immediate-argument
  behavior.
- Structured indices cannot collide through a flattened encoding and supported
  structured payloads do not panic.
- Canonicalization preserves factorization unless the caller explicitly
  expands the expression.
- The signed canonical `SymbolicNet` is the object executed to obtain the final
  Atom.
- Unit, integration, RQFT, and disconnected Monte Carlo spectacles validation
  all pass.
