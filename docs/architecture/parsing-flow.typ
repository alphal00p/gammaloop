= Spenso parsing flow
<spenso-parsing-flow>
This document describes how symbolic atoms are parsed into tensor
networks in `crates/spenso/src/network/parsing`. It focuses on the
control flow, shorthand expansion, opaque leaves, structure inference,
and edge cases that affect Schoonschip-style notation.

#strong[Audit status:] source and documentation reviewed 2026-10-27 against
`4a9ff115`. This attestation does not claim builds or tests were run.
Lifecycle: current implementation architecture.

== Entry Points
<entry-points>
Checked external parsing starts from `NetworkParse` methods such as
`parse_to_atom_net`. `parse_root` rejects nested chain/trace placeholder
consumers, creates a fresh `ParseState`, and reserves written slot-index
spellings before shorthand allocates dummies. Clones share recognition
state, the allocator, and reservations.

`try_from_view_impl` appends to one operation-local `Construction` and
returns a subtree handle. `finish` emits the existing Linnet graph and
store. `try_external_slots` shares that parser with structure-only
construction instead of emitting graph topology.

`try_from_admitted_view` requires already validated syntax with
`ShorthandParsing::Opaque` and `depth_limit = Some(1)`. Its caller owns
syntax, sum-interface, index-scope and chain-nesting admission. An optional
`LeafInterfaceCache` reuses exact expressions' classification and logical
slots only while those facts and classification settings remain valid.
Callbacks and materialized rewrites leave the trusted path.

#link("../../crates/spenso/src/network/parsing/mod.rs#L539")[Checked root parsing],
#link("../../crates/spenso/src/network/parsing/mod.rs#L454")[admitted intake], and
#link("../../crates/spenso/src/network/parsing/construction.rs#L40")[shared construction].

#table(
    columns: 2,
    align: (auto,auto,),
    table.header([*Concept*], [*Role*]),
    [`ParseState`], [Shares arithmetic depth, recognition state, the
    dummy allocator and reserved index names through one parse.],
    [Tensor library], [Maps parsed structures to library tensors when a
    matching key exists.],
    [Function library], [Used by the opaque tensor-expression boundary
    for target tensor types that need it.],
    [`ParseSettings`], [Controls scalar precontraction, depth leaves,
    shorthand mode, and sum parsing behavior.],
)

== Settings
<settings>
#table(
    columns: 2,
    align: (auto,auto,),
    table.header([*Setting*], [*Behavior*]),
    [`precontract_scalars`], [When parsing products or sums, combine
    factors/summands that parse as pure scalars instead of keeping them
    as separate graph nodes.],
    [`take_first_term_from_sum`], [Parse only the first summand. This is
    mainly for structure discovery when all summands are expected to
    expose the same slots.],
    [`depth_limit`], [Products, sums, and powers each consume one nesting
    level. At the limit, composite expressions become leaves. Bracket and
    broadcast roots honor the current depth before recursive handling.],
    [`shorthand_parsing`], [`Expand { schoonschip, trace, chain }`
    independently controls the shorthand families lowered into explicit
    graph structure. `Opaque` keeps recognized shorthand roots as tensor
    leaves and infers their exposed structure.],
    [`parse_composite_scalars_as_tensors`], [Passed into structure
    parsers that want to keep composite scalar expressions inspectable
    instead of immediately storing them as pure scalars.],
    [`strict_tensor_filter`], [`Tagged` (the default) accepts ordinary
    heads tagged as tensors. `TaggedChecked` additionally requires
    representation syntax except for a zero-argument tensor head, while `ContainsReps` also accepts untagged
    heads containing representation syntax. Parser-owned shorthand,
    metric, representation, and broadcast syntax keeps its fixed
    meaning.],
)

== Option Examples
<option-examples>
#table(
    columns: 2,
    align: (auto,auto,),
    table.header([*Option*], [*Behavior*]),
    [`ShorthandParsing::expand_all()`], [`F(p(mink(D)), mink(D,i))` is
    rewritten before generic leaf parsing as
    `F(mink(D,d), mink(D,i)) * p(mink(D,d))`, where `d` is a fresh
    dummy.],
    [`ShorthandParsing::Opaque`], [`chain(bis(D,i), bis(D,j), F(in,out))`
    becomes one tensor leaf. Syntactic inference observes its exposed
    endpoints, then `TensorFromExpression` realizes the leaf. There is no
    selectable `Fast` or `Expanded` inference variant.],
    [`pure_scalar(...)`], [A scalar escape hatch requiring one argument.
    Expansion consumes the wrapper; opaque intake retains the literal
    for domain rewriting. Neither path treats it as a tensor leaf.],
    [`dot(a,b,c)`], [Dot requires two compatible rank-one operands. Wrong
    arity or non-vector operands are rejected before scalar fallback,
    opaque inference, or materialization. Partial contractions with a
    higher-rank tensor's exposed port use metric syntax.],
    [`precontract_scalars = true`], [In `(a + b) * F(mink(D,i))`,
    `a + b` can be kept as one scalar factor next to the tensor network
    instead of a separate graph node.],
    [`take_first_term_from_sum = true`], [`F(mink(D,i)) + G(mink(D,i))`
    is represented by the first summand\'s structure. Use this only when
    the caller already knows the summands expose the same slots.],
    [`depth_limit = Some(1)`], [A nested product can be handed to the
    opaque tensor-expression boundary as a leaf instead of recursively
    expanding every child. This is the same boundary used by opaque
    parsing.],
)

== Main Dispatch
<main-dispatch>
`try_from_view_impl` dispatches by atom kind:

```
Mul  -> try_from_mul
Fun  -> try_from_fun
Add  -> try_from_add
Pow  -> try_from_pow
else -> scalar fallback
```

Products, sums, and powers can hit the depth leaf boundary before
recursing. Function parsing validates dot operands, handles depth-aware
bracket and broadcast roots, classifies tensor syntax, and uses the
opaque tensor-expression boundary for recognized shorthand roots.

Dispatch appends to the shared `Construction`; the following is an
abbreviated control-flow sketch from
#link("../../crates/spenso/src/network/parsing/mod.rs#L588")[the common dispatcher]:

```
match value {
    AtomView::Mul(m) => Self::try_from_mul(construction, m, ...),
    AtomView::Fun(f) => Self::try_from_fun(construction, f, ...),
    AtomView::Add(a) => Self::try_from_add(construction, a, ...),
    AtomView::Pow(p) => Self::try_from_pow(construction, p, ...),
    a => Ok(construction.scalar(Self::scalar_from_expression(retain(a))?)),
}
```

#context if target() == "paged" {
[
== Flow Diagrams
<flow-diagrams>
=== Top Level
<top-level>
#strong[Input];`parse_to_*_net(atom, settings)`
↓
#strong[Validation];Reject nested chain/trace placeholder consumers.
↓
#strong[State];Create shared construction and `ParseState`; reserve written indices.
↓
#strong[Decision];Which `AtomView` variant is this?
Mul
==== `AtomView::Mul`
<atomviewmul>
+ Check depth leaf.
+ Parse factors recursively.
+ Optionally precontract pure scalars.
+ Multiply factor networks.

Add
==== `AtomView::Add`
<atomviewadd>
+ Check depth leaf.
+ Parse summands recursively.
+ Require compatible states.
+ Optionally precontract scalar summands.

Pow
==== `AtomView::Pow`
<atomviewpow>
+ Check depth leaf.
+ Parse the base.
+ Apply integer-power network rules.
+ Use scalar fallback for non-integer powers.

Fun
==== `AtomView::Fun`
<atomviewfun>
+ Check dot arity and compatible rank-one operands.
+ Handle bracket and broadcast roots, honoring the current depth first.
+ Apply the strict tensor-syntax filter.
+ Keep recognized shorthand roots opaque when requested.
+ Materialize generic shorthand before leaf parsing.

Other
==== Other atoms
<other-atoms>
Parse directly as scalars.

=== Function Parsing
<function-parsing>
#strong[Function];`try_from_fun(fun)`
↓
#strong[Decision];`dot` without two compatible rank-one operands?
Yes
==== Syntax error
<syntax-error>
Return `InvalidDotFunction`. This is a parser syntax check, not a tensor
fallback.

No
==== Continue
<continue>
Bracket and broadcast roots are considered next.

↓
#strong[Decision];Is this a bracket or broadcast root?
No
==== Continue
<continue-1>
Apply the strict tensor-syntax filter.

Yes
==== Bracket or broadcast boundary
<broadcast-wrapper>
Check the current depth first; at the limit, return an opaque leaf.
Otherwise a bracket parses its arguments as a product. A broadcast
requires one argument and parses it recursively. A scalar-tagged head
returns that inner result; another broadcast head wraps it. These
branches return before ordinary tensor classification.

↓
#strong[Decision];Does `is_tensorial(strict_tensor_filter)` accept the
root?
Yes
==== Tensor syntax
<tensor-syntax>
Continue to shorthand dispatch.

No
==== Scalar path
<scalar-path>
Parse the function as a scalar. `pure_scalar(x)` requires one argument;
expansion unwraps it, while opaque intake retains the literal.

↓
#strong[Decision];Is `shorthand_parsing` opaque and is this a recognized
shorthand root?
Yes
==== Opaque tensor leaf
<opaque-tensor-leaf>
+ Infer exposed structure from the original atom.
+ Observe syntactic slots, or reuse valid admitted observations.
+ Call `TensorFromExpression`.

No
==== Expand-mode parsing
<expand-mode-parsing>
Enter symbol dispatch and materialization.

↓
#strong[Decision];Is the symbol special expand syntax?
Yes
==== Special syntax
<special-syntax>
- Materialize enabled `chain` and `trace` roots.
- Materialize selected compact Schoonschip families.
- Bracket and broadcast roots returned at the earlier boundary.

No
==== Ordinary function candidate
<ordinary-function-candidate>
Run shorthand materialization before ordinary leaf parsing.

↓
#strong[Boundary];`materialize_shorthand` returns a shared-construction subtree handle.
Topology
==== Materialized network
<materialized-network>
+ `chain` and `trace` return parsed topology.
+ Compact dot/Schoonschip rewrites recurse through normal parsing.
+ This prevents compact arguments from becoming metadata.

Fixed point
==== Ordinary tensor leaf
<ordinary-tensor-leaf>
+ If the Schoonschip atom rewriter makes no change, do not recurse.
+ Run `S::structure_from_parser` with existing recognition and admission state.
+ Use a library tensor if a key exists.
+ Otherwise retain a rank-zero result as a scalar, or call `tensor_from_leaf`.

Empty structure
==== Scalar fallback
<scalar-fallback>
If structure parsing reports `StructureError::EmptyStructure`, parse the
function as a scalar expression. Other structure errors propagate.

=== Shorthand Materialization
<shorthand-materialization>
#strong[Root];`materialize_shorthand(construction, fun, ...) -> NodeIndex`
↓
#strong[Decision];Root is an enabled `chain` or `trace`?
Yes
==== Topology materialization
<topology-materialization>
+ Materialize compact chain endpoints when configured.
+ Allocate chain or cyclic trace links.
+ Replace `in`/`out`.
+ Materialize compact Schoonschip syntax inside factors.
+ Append the parsed factors and return their subtree handle.

No
==== Compact syntax check
<compact-syntax-check>
Use `SchoonschipMaterializer` to rewrite compact dot or compact vector
arguments.

Dot
==== Compact scalar product
<compact-scalar-product>
`dot(p(rep), q(rep))` or `g(p(rep), q(rep))`

+ Allocate one fresh dummy slot `d`.
+ Build `p(d) * q(d)`.
+ Parse that expression into the same construction.

Arg
==== Compact vector argument
<compact-vector-argument>
`F(..., p(rep), ...)`

+ Allocate fresh dummy slot `d`.
+ Rebuild `F(..., d, ...)`.
+ Parse `F(..., d, ...) * p(d)`.

No change
==== Ordinary leaf candidate
<ordinary-leaf-candidate>
Parse as an ordinary tensor leaf. This is the recursion stop.

=== Chain And Trace
<chain-and-trace>
chain
==== `materialize_chain_shorthand(chain(start, end, factors...))`
<materialize_chain_shorthandchainstart-end-factors...>
+ Parse concrete endpoints, or materialize compact endpoints when
  Schoonschip expansion permits it.
+ Branch: no factors parses `g(start, end)`.
+ Branch: factors allocate links, replace `in`/`out`, materialize
  shorthand, parse, and multiply.

trace
==== `materialize_trace_shorthand(trace(rep, factors...))`
<materialize_trace_shorthandtracerep-factors...>
+ Parse `rep` as a representation.
+ Branch: no factors parses `dim(rep)`.
+ Branch: non-empty traces allocate one cyclic link per factor.
+ Replace placeholders, optionally materialize factor shorthand, parse,
  and multiply.
]
}

== Function Flow
<function-flow>
Function parsing has the most important ordering rules.

+ Reject malformed `dot` syntax first: require two compatible rank-one
  operands.
+ Handle bracket and broadcast roots before ordinary classification.
  Honor the current depth leaf first. A bracket parses its operands as a
  product. A broadcast requires one argument and preserves its function
  head only when the head is not tagged scalar.
+ Classify tensor syntax under `strict_tensor_filter`, reusing an
  admitted observation only when valid. Rejected roots take the scalar
  path. Expansion unwraps `pure_scalar(x)`; opaque intake retains it.
+ If `shorthand_parsing` is `Opaque` and the root is recognized
  shorthand, infer its canonicalized structure and realize one leaf.
+ Otherwise enter expanded parsing. The `Expand` fields independently
  select chain, trace and compact Schoonschip materialization. A disabled
  shorthand remains an inferred leaf; an unchanged ordinary head uses
  ordinary leaf parsing.
+ Ordinary leaves infer with `S::structure_from_parser`, use a library
  tensor and canonical layout when available, retain rank-zero results
  as scalars, or call the target's fallible `tensor_from_leaf` boundary.
  `EmptyStructure` is a scalar fallback; other errors propagate.

The current checked boundaries are abbreviated below; omitted arguments
carry the same construction, parser state, libraries, settings and source
retention policy. See
#link("../../crates/spenso/src/network/parsing/mod.rs#L906")[function dispatch] and
#link("../../crates/spenso/src/network/parsing/mod.rs#L1027")[scalar function intake].

```
if symbol == dot
    && matcher.compact_inner_product_parts::<Aind>(value).is_none()
{
    return Err(TensorNetworkError::InvalidDotFunction(...));
}

if symbol == bracket || symbol.has_tag(&broadcast) {
    // Honor the current depth leaf before recursive bracket/broadcast parsing.
    return Self::parse_expanded_function(construction, value, ...);
}

// TensorialSyntax and an optional valid boundary cache classify this root.
if !tensorial {
    return Self::parse_scalar_function(
        construction, value, settings.shorthand_parsing, retain,
    );
}
if settings.shorthand_parsing == ShorthandParsing::Opaque
    && Self::is_shorthand_function(value, &state)
{
    return Self::as_inferred_leaf(construction, value.as_view(), ...);
}
Self::parse_expanded_function(construction, value, ...)
```

Important: materialization must run before the generic leaf parser.
Otherwise `F(p(rep), slot)` can be accepted as a lower-rank tensor with
`p(rep)` stored as metadata, losing the required `p(dummy)` factor.

== Expand vs Opaque
<expand-vs-opaque>
=== Expand
<expand>
Selected shorthand families are rewritten into explicit parser input.
The `schoonschip`, `trace`, and `chain` fields can be configured
independently. Expanded chain and trace roots allocate links; compact
vector arguments allocate fresh slots. All append to the same
operation-local graph construction.

=== Opaque
<opaque>
A recognized shorthand root that survives strict tensor filtering and
broadcast handling becomes one leaf. Structure is inferred first, then
`TensorFromExpression` constructs the target tensor from the original
expression and inferred structure. Ordinary tensor heads still use
regular leaf parsing.

Opaque mode is checked before chain or trace expansion, but after strict
tensor classification and broadcast handling. Non-tensor functions stay
scalar in opaque mode; depth-limited composite expressions use the
separate leaf boundary described below.

== Opaque Leaf Boundary
<opaque-leaf-boundary>
Opaque shorthand parsing separates interface observation from realization.
See #link("../../crates/spenso/src/network/parsing/mod.rs#L726")[the inferred-leaf boundary]
and #link("../../crates/spenso/src/network/parsing/tensor_from_expression.rs#L36")[target realization].

+ `StructureFromAtom::structure_from_parser` observes syntactic slots
  using parser-local recognition and admission state. Admitted shallow
  parsing may reuse exact boundary observations; there is no selectable
  expanded-inference mode.
+ `TensorFromExpression::tensor_from_expression` receives the source as
  `AtomOrView` and a canonicalized structure, and returns the realized
  tensor together with its layout. Symbolic targets can retain the
  expression and source interface. Eager targets expand and execute a
  nested network, then obtain its tensor result with scales intact.

Abbreviated parser boundary:

```
let structure = S::structure_from_parser(state)?;
let realized = T::tensor_from_expression(
    expression, structure, library, function_library, settings,
)?;
let (tensor, layout) = realized.into_parts();
construction.tensor(tensor, layout)
```

Eager realization uses `Sequential` with `MinIntermediateCost` and
`result_tensor`, including completed lazy scaled tensors and sums. Its
returned layout describes the realized coordinate order. Internal ports
can disappear during eager contraction; the parser must not reattach the
larger source layout afterward. Symbolic targets preserve their source
layout because their contractions still belong to the outer network.

== Structure Inference
<structure-inference>
Fast structure inference is syntactic and intentionally non-semantic.

#table(
    columns: 2,
    align: (auto,auto,),
    table.header([*Syntax*], [*Inferred structure*]),
    [Function], [Direct slot arguments expose slots. `aind(...)` bundles
    are flattened; malformed bundles return a slot error. Other arguments are metadata.],
    [Product], [Merge the structures of all factors. Scalar factors
    contribute an empty structure.],
    [Sum], [Use the first summand\'s structure. The parser separately
    checks compatibility when it parses full sums.],
    [Power], [Scalars stay scalar. Even powers of fully self-dual
    tensors become scalar. Odd integer powers keep the base structure.
    Other tensor powers are rejected.],
    [`chain`], [Open endpoints are external slots. Factor slots are
    scanned recursively. `in` and `out` placeholders are wiring labels,
    not visible slots.],
    [`trace`], [The representation argument is not an external slot.
    Factor slots are scanned recursively.],
)

== Materialization
<materialization>
#link("../../crates/spenso/src/network/parsing/materialization.rs#L895")[`materialize_shorthand`]
appends to the shared `Construction` and returns
a subtree handle. The root parser later finalizes the network. Enabled
chain and trace roots materialize their topology; compact syntax rewrites
re-enter dispatch only after an actual change. Disabled shorthand stays
an inferred leaf and ordinary fixed points use ordinary leaf parsing.

`SchoonschipMaterializer` is the narrower atom rewriter. It keeps the
current atom and additional factors, and the completed expression returns
to the shared parser. Materialized rewrites clear trusted admission/cache
state before recursively parsing generated syntax.

Abbreviated control flow:

```
fn materialize_shorthand(construction: &mut Construction<...>, ...)
    -> Result<NodeIndex, TensorNetworkError<K, Symbol>>
{
    // Handle enabled chain/trace topology and configured compact rewrites.
    if no_actual_rewrite {
        if root_chain_disabled || root_trace_disabled || has_shorthand {
            return Self::as_inferred_leaf(construction, ...);
        }
        return Self::parse_regular_function_leaf(construction, ...);
    }
    Self::try_from_view_impl(
        construction, materialized.as_view(), state.materialized(), ...,
    )
}
```

#link("../../crates/spenso/src/network/parsing/materialization.rs#L135")[Atom-level materializer]

```
pub(super) fn materialize_shorthand(&self, value: AtomView<'_>) -> Atom {
    self.materialize_shorthand_root(value)
        .map(SchoonschipMaterialization::into_expression)
        .unwrap_or_else(|| value.to_owned())
}
```

#link("../../crates/spenso/src/network/parsing/materialization.rs#L280")[Compact scalar-product materialization]

```
let (lhs, rhs, lhs_rep, rhs_rep) = self
    .state
    .matcher
    .borrow_mut()
    .compact_inner_product_parts::<Aind>(value)?;
let lhs_slot = self.state.slot(&lhs_rep);
let rhs_slot = rhs_rep.slot::<Aind, _>(lhs_slot.aind());
```

Chain and trace materialization first chooses concrete `in`/`out`
replacements. The same compact-syntax helper then expands configured
shorthand inside the factor; rewritten syntax re-enters checked parsing.

=== Compact Vector Argument
<compact-vector-argument-1>
```
F(..., p(rep), ...)
-> F(..., rep(dummy), ...) * p(rep(dummy))
```

The dummy slot is allocated from the active `ParseState`, so it shares
the namespace of the surrounding network parse.

=== Compact Scalar Product
<compact-scalar-product-1>
```
g(p(rep), q(rep.dual()))
-> p(rep(dummy)) * q(rep.dual()(dummy))
```

The compact representations must match under duality, including dimension.
For self-dual representations, `rep.dual()` is `rep`.
`dot(p(rep), q(rep.dual()))` follows the same rule.

The stripped representation in `p(rep)` means that the vector replaces
an omitted tensor slot. It is valid only where such a slot exists:
inside another tensor argument, or inside a recognized compact scalar
product head such as `g` or `dot`. For example, `p!(q!(mink!(4)))` is
valid shorthand and materializes like `p(mink(4, d)) * q(mink(4, d))`. A
bare product such as `p!(mink!(4)) * q!(mink!(4))` is malformed input,
because the stripped vectors are not replacing slots.

=== Compact Vector Detection
<compact-vector-detection>
A compact vector is a function that:

- has a tensor-tagged head;
- is not `metric` or `dot`;
- is not itself a representation;
- has no explicit slot argument;
- has exactly one direct argument matching the representation wildcard
  convention.

Sums are accepted only when every summand is a compact vector with the
same representation. Unary broadcasts, projectors, and products preserve
compact-vector syntax when they contain exactly one compatible vector and
otherwise scalar factors. Powers are not compact vector syntax.

== Chain Expansion
<chain-expansion>
When chain expansion is enabled, `chain(start, end, factors...)` first
tries to parse each endpoint as a concrete slot. If Schoonschip
expansion is enabled, a compact vector endpoint can instead materialize
into a fresh endpoint slot plus an additional rank-one factor. The
network-level materializer then chooses internal links, rewrites each
factor\'s `in`/`out` placeholders, and optionally expands compact
arguments inside the factor.

+ Materialize compact endpoints when configured and retain their
  additional factors.
+ With no chain factors, add `g(start, end)`.
+ Each factor gets a left and right link slot.
+ `in` and `out` placeholders are replaced syntactically.
+ The factor is passed through shorthand materialization.
+ Each materialized factor is parsed into one or more factor networks.
+ The factor networks are multiplied in chain order.

If a materialized chain factor is a product with scalar factors and
exactly one tensor factor, scalar factors are split out so they can be
handled separately. Otherwise the whole factor product is parsed as one
expression.

```
chain(s, e, F(in, out, p(rep)))
-> replace placeholders: F(s, e, p(rep))
-> materialize factor:   F(s, e, d) * p(d)
-> parse product into the chain network
```

#link("../../crates/spenso/src/network/parsing/materialization.rs#L1013")[Chain materialization]

The materializer obtains endpoint slots through the operation-local
`SchoonschipMaterializer`, replaces factor placeholders using the shared
`SlotMatcher`, and appends parsed factors to the shared construction.

== Trace Expansion
<trace-expansion>
When trace expansion is enabled, `trace(rep, factors...)` starts from a
representation, not endpoint slots. The network-level shorthand
materializer returns the representation dimension for an empty trace and
otherwise assigns cyclic links uniformly across all factors. Depending
on the Schoonschip expansion settings, it can materialize compact
arguments inside each factor after placeholder replacement and before
parsing.

+ No factors: parse the representation dimension.
+ For a non-empty trace, allocate one cyclic link per factor.
+ Replace each factor\'s placeholders with its link and the dual of the
  next cyclic link.
+ Materialize shorthand inside each factor when configured.
+ Multiply the factor networks.

```
trace(rep, F(in, out, p(rep)))
-> allocate one cyclic link
-> replace placeholders: F(link, link*, p(rep))
-> materialize factor:   F(link, link*, d) * p(d)
```

#link("../../crates/spenso/src/network/parsing/materialization.rs#L1149")[Trace materialization]

Each factor receives its cyclic link and the dual of the next link.
`ChainExpansion::replace_placeholders` uses the shared matcher while
rewriting those slots. Parsed factors append to the same construction;
a trace with no factors parses the representation dimension.

== Products, Sums, Powers
<products-sums-powers>
#table(
    columns: 2,
    align: (auto,auto,),
    table.header([*Syntax*], [*Parsing behavior*]),
    [Product], [Parse every factor recursively. If scalar precontraction
    is enabled, pure scalar factors are combined. The old special case
    `p(rep) * q(rep)` is intentionally not a shorthand materialization
    path.],
    [Sum], [Parse summands recursively and require compatible network
    states. Pure scalar summands may be precontracted into one scalar
    sum.],
    [Power], [Integer powers recurse into the base. Tensor powers are
    allowed only when the network state supports the exponent.
    Non-integer powers fall back to scalar parsing.],
)

Scalar precontraction retains an entirely scalar source sum unchanged.
For mixed sums, child admission precedes coefficient assembly. Exact
`Natural`/`Large` leading coefficients permit one `Atom::add_many` call;
other coefficient domains keep the source's existing rest-then-first
fold. Function and power bodies remain literal, and inexact arithmetic
is not reassociated.

Positive powers that lower shorthand with fresh internal dummies parse
additional copies from the original base with the shared allocator.
A self-dual tensor square also uses a product so the selected contraction
strategy controls that contraction. Other supported cases use an ordinary
power node. Negative odd tensor powers are rejected unless the base is
scalar. Scalar extraction under repeated powers preserves contraction
boundaries through `scalar_power_base`. See
#link("../../crates/spenso/src/network/parsing/mod.rs#L1223")[power parsing].

== Edge Cases
<edge-cases>
#table(
    columns: 2,
    align: (auto,auto,),
    table.header([*Input*], [*Behavior*]),
    [`F(p(rep), slot)`], [Must materialize before ordinary function
    parsing. Expected expansion is `F(dummy, slot) * p(dummy)`, leaving
    `slot` external and the dummy internal.],
    [`p(rep)`], [A standalone compact vector root is not materialized.
    It can only become a vector when used in a materializable argument
    or scalar-product position.],
    [`p(rep) * q(rep)`], [Not valid shorthand in the parser. The product
    path does not materialize this into a dot product.],
    [`dot(a,b,c)`], [Invalid syntax. Dot requires two compatible
    rank-one operands, checked before any scalar or opaque fallback.],
    [`g(slot, p(rep))`], [Materializes to `g(slot, dummy) * p(dummy)`.
    The contraction is represented by the network, not by a pre-parser
    algebra rewrite.],
    [`g(p(rep), q(rep))`], [For matching self-dual reps, materializes to
    `p(dummy) * q(dummy)`.],
    [Dualizable compact scalar products], [Matching dual representations
    share one fresh abstract index while retaining their respective
    orientations: `g(p(rep), q(rep.dual()))` becomes
    `p(rep(dummy)) * q(rep.dual()(dummy))`.],
    [Compact vector endpoints in chains], [When Schoonschip expansion is
    enabled, compact endpoints materialize into fresh slots and
    additional rank-one factors. With that expansion disabled, endpoints
    must be concrete slots.],
    [Metadata false positives], [Compact vector heads must be tensor-tagged,
    with exactly one direct representation argument and no slots.
    Untagged metadata wrappers do not qualify.],
    [Nested chains and traces], [Rejected before parsing or structure
    inference because nested containers would share ambiguous `in`/`out`
    placeholders.],
    [Opaque interface discovery], [Syntactic inference observes exposed
    slots without building a complete arithmetic network. Admitted
    observations reuse unchanged interfaces. Component execution is a
    separate correctness oracle, not an inference-mode setting.],
    [Opaque scalar result], [A recognized shorthand can infer tensorial
    structure and still finalize to a scalar tensor. Opaque mode sends
    that shorthand through `TensorFromExpression`; ordinary function
    heads continue through expanded parsing.],
    [Depth leaves], [Depth leaves use the opaque tensor-expression
    boundary. Non-tensorial expressions still fall back to pure scalars
    unless composite scalar tensor leaves are explicitly enabled.],
)

== Short Algorithm
<short-algorithm>
This sketch omits generic types and error propagation. Every returned
handle belongs to the same operation-local construction.

```
parse_checked(atom):
  validate structural chain/trace scope
  reserve written index spellings
  root = dispatch(construction, atom, state)
  return construction.finish(root)

parse_fun(fun):
  if dot lacks two compatible rank-one operands:
    return InvalidDotFunction
  if bracket or broadcast:
    honor current depth leaf
    return parsed product or unary broadcast boundary
  if tensor classification rejects the root:
    consume pure_scalar only during expansion
    return scalar handle
  if opaque mode and recognized shorthand:
    infer structure using parser state
    realize original expression and returned canonical layout
    return tensor handle
  return materialize_shorthand(construction, fun)

materialize_shorthand(construction, fun):
  if an enabled chain or trace:
    materialize endpoints or representation dimension as appropriate
    allocate chain links or cyclic trace links
    replace placeholders through shared matcher
    parse configured factor rewrites into this construction
    return completed subtree handle
  select effective compact expansion mode
  rewrite enabled compact syntax
  if no actual change:
    return inferred shorthand leaf or ordinary tensor leaf
  clear trusted admission/cache path for generated syntax
  return dispatch(construction, rewritten expression, state)
```
