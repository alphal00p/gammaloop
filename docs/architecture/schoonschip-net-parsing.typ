= Schoonschip network architecture
<schoonschip-network-parsing>

The shared symbolic tensor owns index contraction. Spenso owns parsing,
logical port layouts, the operation graph, and concrete network execution;
Idenso assigns metric/vector algebra to that graph. Python bindings convert
arguments and wrap the shared result.

== One symbolic carrier

`SymbolicTensor` carries an expression together with its logical interface.
Reduction and collection retain a plain Symbolica Atom in that same abstraction;
there is no separate Python composition tensor or contraction-specific tensor.

General composition, index substitution, contraction, and interface checks live
in #link("../../crates/idenso/src/tensor/mod.rs")[the shared tensor module].
Composition computes one product plan and batches operand relabeling. Unchanged
branches are reused, and arithmetic reconstruction uses bulk sum/product
construction where normalization semantics permit it.

== Parsing and contraction

#link("../../crates/idenso/src/shorthands/schoonschip/slot_contraction/components/factorized.rs")[The factorized contraction implementation]
reads Spenso's borrowed symbolic network. Arithmetic scopes are parsed once;
tensorial shorthand remains opaque. Graph leaves retain borrowed atom views,
and the graph's ordered layout supplies occurrence-local port identities.
Unselected arithmetic subtrees remain opaque terminals of the existing
component reducer.

The same admitted depth-one Spenso parser supplies this shallow graph. Exact
leaf-interface observations survive trusted replacements, so unchanged large
leaves do not repeat interface discovery. Selected scope opening and callback
invalidation follow the ordinary parser boundaries rather than a second
contraction path. `contract` and `simplify_algebra` schedule this implementation;
the latter permits prerequisite connections only for its selected identity.
Representation filters select compatible connections and preserve other ports.
Chain and trace collection share notation primitives, and collecting a trace
never evaluates it. `reduction_status` distinguishes complete, deferred, and
budget-capped work without discarding the exact remaining expression.

The contraction frontier chooses an incident factor, visits its alternatives,
and carries index overrides to adjacent leaves. Equivalent remaining states
merge before emitting an expression. A distinct relabeled tensor is constructed
once; generated scalar sums remain factorized. Neither rule application
nor contraction materializes the complete distributed numerator.

A powered scalar subexpression owns its internal dummy pairs. Contract that
scope before applying its power; copies must not accidentally identify their
dummies. Open powers retain the existing cross-copy contraction semantics.

`contract(...)` returns a `TensorExpression` in Python.
Order and rank-one handling are
policies of this operation, not separate contraction engines. Metric-only
callers retain explicit vector contractions until their chosen later stage.
Algebraic validity and completion are separate: a frontier budget can leave an
exact pending contraction, so an incomplete result must never acquire the
completed-result shortcut. The read-only `contraction_complete` property means
that the default metric/vector contractor is certified complete. `false` also
covers an unproved result; metric-only work does not establish default
completion. It says nothing about arbitrary tensor component contraction.

== Callbacks and interface proofs

Intrinsic normalization permits the component reducer to use graph incidence
without rebuilding tensor syntax after every substitution. User normalizers,
custom metric behavior, and opaque payloads require the existing ordered
substitution schedule and checked result construction.

For example, contracting `g(a,b)*T(a)` can invoke a normalizer that turns
`T(b)` into a scalar. The original interface cannot be retained blindly.
Logical port identities and their order matter independently of a multiset of
representations; unresolved ports, metadata, and typed zero interfaces remain
part of the shared invariant.

A completed certified result can reuse its observations on a default rerun.
Changed regions and callback-sensitive results invalidate the relevant proofs.
Unselected factors retain their callbacks without being evaluated as a side
effect of unrelated contraction.

== Dots and concrete components

The intrinsic `g` hook and the shared dot normalizer own compact vector and
metric normalization. `to_dots` and `undo_dots` are symbolic representation
operations. A caller that also needs index contraction requests `contract`
explicitly before the representation change.

The typed `expand_dots()` returns an `Atom` at the explicit finite component
Spenso execution boundary, rather than sealing concrete component markers as
symbolic ports. It realizes only selected dots; unrelated numerator factors
remain untouched and symbolic dimensions remain unsupported.
Symbolically opening a dot is insufficient when a vector normalizer creates
an internal contraction. Batch execution may share construction only where
normalization and execution metadata establish that doing so preserves
callbacks; otherwise preserve the original parse/execute ordering. Each input
keeps its dummy scope and its own unsupported-input outcome.

== Explicit materialization and validation

`to_expression()` exposes the ordinary Symbolica expression in Python.
`expand()` explicitly materializes its polynomial. Neither conversion runs
implicitly during tensor inference or contraction.

Validation covers logical ordering, unresolved ports, callback-induced rank
changes, copied dummy scopes, and identity/genuine index substitutions. FORM
and exact HEP components provide independent algebra checks. Captured graph
numerators remain factorized even in diagnostics; component comparisons record
their momentum assignments and dimension scope.

The retired `NetworkSchoonschip` route, its ordering strategies, and the old
pattern-rule names are not alternative public APIs. The ordinary
`TensorNetwork` remains the concrete execution and rendering object.
