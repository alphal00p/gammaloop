= Spenso parsing flow
<spenso-parsing-flow>
This document describes how symbolic atoms are parsed into tensor
networks in `crates/spenso/src/network/parsing`. It focuses on the
control flow, shorthand expansion, opaque leaves, structure inference,
and edge cases that affect Schoonschip-style notation.

#strong[Audit status:] reviewed 2026-08-17 against `c9f4e32acd2c`.
Lifecycle: current implementation architecture.

== Entry Points
<entry-points>
Parsing starts from `NetworkParse` methods such as `parse_to_atom_net`.
The parser creates a fresh `ParseState` and calls `try_from_view_impl`
on the input `AtomView`.

#table(
    columns: 2,
    align: (auto,auto,),
    [`ParseState`], [Tracks recursion depth and owns the dummy-index
    allocator used by shorthand materialization.],
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
    [`precontract_scalars`], [When parsing products or sums, combine
    factors/summands that parse as pure scalars instead of keeping them
    as separate graph nodes.],
    [`take_first_term_from_sum`], [Parse only the first summand. This is
    mainly for structure discovery when all summands are expected to
    expose the same slots.],
    [`depth_limit`], [Turns composite expressions into leaves once the
    parser reaches the configured depth. The current checks live on add,
    mul, and pow parsing.],
    [`depth_is_product_depth`], [When true, only products increment
    depth. When false, sums and powers increment depth too.],
    [`shorthand_parsing`], [`Expand { schoonschip, trace, chain }`
    independently controls the shorthand families lowered into explicit
    graph structure. `Opaque` keeps recognized shorthand roots as tensor
    leaves and infers their exposed structure.],
    [`parse_composite_scalars_as_tensors`], [Passed into structure
    parsers that want to keep composite scalar expressions inspectable
    instead of immediately storing them as pure scalars.],
    [`strict_tensor_filter`], [`Tagged` (the default) accepts ordinary
    heads tagged as tensors. `TaggedChecked` additionally requires
    representation syntax, while `ContainsReps` also accepts untagged
    heads containing representation syntax. Parser-owned shorthand,
    metric, representation, and broadcast syntax keeps its fixed
    meaning.],
)

== Option Examples
<option-examples>
#table(
    columns: 2,
    align: (auto,auto,),
    [`ShorthandParsing::expand_all()`], [`F(p(mink(D)), mink(D,i))` is
    rewritten before generic leaf parsing as
    `F(mink(D,d), mink(D,i)) * p(mink(D,d))`, where `d` is a fresh
    dummy.],
    [`ShorthandParsing::Opaque { inference: Fast }`], [`chain(bis(D,i), bis(D,j), F(in,out))`
    becomes one tensor leaf. Its visible structure is inferred
    syntactically as the two endpoints `bis(D,i)` and `bis(D,j)`, then
    `TensorFromExpression` realizes the leaf.],
    [`ShorthandParsing::Opaque { inference: Expanded }`], [The same
    expression is expanded in a temporary network only to read dangling
    slots. The final parse still stores one opaque leaf.],
    [`pure_scalar(...)`], [This wrapper is a hard scalar escape hatch.
    Tensorial-syntax classification rejects the wrapper, then the scalar
    parser unwraps its single argument, so opaque mode never sees it as
    a tensor leaf.],
    [`dot(a,b,c)`], [Dot syntax is accepted only with two arguments. A
    three-argument dot is rejected before scalar fallback, opaque
    inference, or expand-mode materialization.],
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
    [`depth_is_product_depth = false`], [Sums and powers count toward
    the depth limit too. With the default `true`, only product nesting
    increases parse depth.],
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

Products, sums, and powers may hit the depth leaf boundary before they
recurse. Function parsing validates dot arity, applies
`is_tensorial(strict_tensor_filter)`, handles broadcast wrappers, and
uses the opaque tensor-expression boundary only for recognized shorthand
roots.

#link("../../crates/spenso/src/network/parsing/mod.rs#L362")[crates/spenso/src/network/parsing/mod.rs:362]

```
match value {
    AtomView::Mul(m) => Self::try_from_mul(...),
    AtomView::Fun(f) => Self::try_from_fun(...),
    AtomView::Add(a) => Self::try_from_add(...),
    AtomView::Pow(p) => Self::try_from_pow(...),
    a => Ok(Network::from_scalar(a.try_into()?)),
}
```

== Flow Diagrams
<flow-diagrams>
=== Top Level
<top-level>
#strong[Input];`parse_to_*_net(atom, settings)`
↓
#strong[State];Create `ParseState` with depth and dummy allocator.
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
+ Run dot-arity and strict tensor-syntax checks.
+ Handle broadcast wrappers.
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
#strong[Decision];`dot` with anything other than two args?
Yes
==== Syntax error
<syntax-error>
Return `InvalidDotFunction`. This is a parser syntax check, not a tensor
fallback.

No
==== Continue
<continue>
Scalar escape hatches are considered next.

↓
#strong[Decision];Does `is_tensorial(strict_tensor_filter)` accept the
root?
Yes
==== Tensor syntax
<tensor-syntax>
Continue to broadcast and shorthand dispatch.

No
==== Scalar path
<scalar-path>
Parse the whole function as a scalar, except that `pure_scalar(x)`
unwraps its one argument.

↓
#strong[Decision];Is the head tagged as `broadcast`?
No
==== Continue
<continue-1>
Consider opaque shorthand handling.

Yes
==== Broadcast wrapper
<broadcast-wrapper>
Require one argument, parse it recursively, then apply the broadcast
head to the result.

↓
#strong[Decision];Is `shorthand_parsing` opaque and is this a recognized
shorthand root?
Yes
==== Opaque tensor leaf
<opaque-tensor-leaf>
+ Infer exposed structure from the original atom.
+ `Fast`: syntactic structure inference.
+ `Expanded`: expanded net dangling slots.
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
- `bracket`: parse args and multiply.
- `chain`, `trace`, and selected Schoonschip families according to their
  expansion settings.

No
==== Ordinary function candidate
<ordinary-function-candidate>
Run shorthand materialization before ordinary leaf parsing.

↓
#strong[Boundary];`materialize_shorthand(fun)` returns a network.
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
+ Run `S::parse` for fast structure inference.
+ Use a library tensor if a key exists.
+ Otherwise concretize the parsed tensor shell.

Empty structure
==== Scalar fallback
<scalar-fallback>
If structure parsing reports `StructureError::EmptyStructure`, parse the
function as a scalar expression. Other structure errors propagate.

=== Shorthand Materialization
<shorthand-materialization>
#strong[Root];`materialize_shorthand(fun) -> Network`
↓
#strong[Decision];Root is an enabled `chain` or `trace`?
Yes
==== Topology materialization
<topology-materialization>
+ Materialize compact chain endpoints when configured.
+ Allocate chain or cyclic trace links.
+ Replace `in`/`out`.
+ Materialize compact Schoonschip syntax inside factors.
+ Parse and return the resulting network.

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
+ Parse that expression into a network.

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

== Function Flow
<function-flow>
Function parsing has the most important ordering rules.

+ Reject malformed `dot` syntax first. Dot has exactly two arguments in
  parser syntax.
+ Apply `is_tensorial(settings.strict_tensor_filter)`. If it rejects the
  root, parse it as a scalar; `pure_scalar(x)` is unwrapped in that
  scalar path.
+ If the head has the `broadcast` tag, require one argument, parse that
  argument recursively, and apply the head to the result.
+ If `shorthand_parsing` is `Opaque` and the root is recognized
  shorthand (`chain`, `trace`, `dot`, or compact Schoonschip syntax),
  parse the whole function as one inferred tensor leaf. Ordinary tensor
  heads do not take this branch.
+ Otherwise enter expanded-function parsing. `bracket` parses its
  arguments and multiplies them; other roots enter the shorthand
  materializer.
+ The `Expand` fields independently decide whether `chain`, `trace`, and
  compact Schoonschip syntax materialize. A disabled shorthand root is
  retained as a fast-inferred leaf. An unchanged ordinary head uses
  regular leaf parsing.
+ Try `S::parse`. If it succeeds, use a library tensor when possible,
  otherwise concretize a tensor shell.
+ If structure parsing reports `StructureError::EmptyStructure`, parse
  the function as a scalar expression. Propagate other structure errors.

#link("../../crates/spenso/src/network/parsing/mod.rs#L606")[crates/spenso/src/network/parsing/mod.rs:606]

```
if symbol == SPENSO_TAG.dot && value.get_nargs() != 2 {
    return Err(TensorNetworkError::InvalidDotFunction(...));
}

if !value.as_view().is_tensorial(settings.strict_tensor_filter) {
    return Self::parse_scalar_function(value);
}
```

#link("../../crates/spenso/src/network/parsing/structure_inference.rs#L50")[crates/spenso/src/network/parsing/structure\_inference.rs:50]

```
pub trait AtomStructureExt {
    fn infer_structure<S: StructureFromAtom>(...);
    fn is_tensorial(&self, filter: StrictTensorFilter) -> bool;
}
```

#link("../../crates/spenso/src/network/parsing/mod.rs#L626")[crates/spenso/src/network/parsing/mod.rs:626]

```
if let Some(inference) = settings.shorthand_parsing.opaque_inference()
    && Self::is_shorthand_function(value)
{
    return Self::as_inferred_leaf(...);
}
```

#link("../../crates/spenso/src/network/parsing/mod.rs#L729")[crates/spenso/src/network/parsing/mod.rs:729]

```
Self::materialize_shorthand(
    value,
    state,
    library,
    function_library,
    settings,
)
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
vector arguments allocate fresh slots.

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
Opaque shorthand parsing is split into two steps.

+ `StructureFromAtom::structure_from_atom` determines the exposed slots.
  `Fast` uses a syntactic walk. `Expanded` builds an expanded shorthand
  network and reads dangling slots.
+ `TensorFromExpression::tensor_from_expression` turns the original atom
  plus inferred structure into a tensor. Symbolic tensors keep the
  original expression on a tensor with the inferred structure. Concrete
  and parametric tensor classes parse the same atom again with
  `ShorthandParsing::expand_all()`, execute that nested network, and
  return the resulting tensor or tensor scalar.

#link("../../crates/spenso/src/network/parsing/structure_inference.rs#L162")[crates/spenso/src/network/parsing/structure\_inference.rs:162]

```
match mode {
    StructureInferenceMode::Fast => Self::leaf_structure_from_atom(value),
    StructureInferenceMode::Expanded =>
        Self::expanded_shorthand_structure_from_atom(value),
}
```

#link("../../crates/spenso/src/network/parsing/tensor_from_expression.rs#L25")[crates/spenso/src/network/parsing/tensor\_from\_expression.rs:25]

```
/// Builds a tensor leaf from an opaque expression once its structure is known.
pub trait TensorFromExpression<S, Sc, K, FK, Aind, Lib, FunLib> {
    fn tensor_from_expression(...);
}
```

#link("../../crates/spenso/src/network/parsing/tensor_from_expression.rs#L120")[crates/spenso/src/network/parsing/tensor\_from\_expression.rs:120]

```
let mut expanded_settings = settings.clone();
expanded_settings.shorthand_parsing = ShorthandParsing::expand_all();

let mut network = Network::try_from_view_with_function_library(
    expression,
    tensor_library,
    function_library,
    &expanded_settings,
)?;
network.execute::<Sequential, SmallestDegree, ...>(...)?;
```

#link("../../crates/idenso/src/tensor/mod.rs#L118")[crates/idenso/src/tensor/mod.rs:118]

```
let mut tensor = structure.structure;
...
tensor.expression = expression.to_owned();
Ok(tensor)
```

== Structure Inference
<structure-inference>
Fast structure inference is syntactic and intentionally non-semantic.

#table(
    columns: 2,
    align: (auto,auto,),
    [Function], [Direct slot arguments expose slots. `aind(...)` bundles
    are flattened. Other arguments are metadata.],
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
The parser\'s `materialize_shorthand` method is the configurable
shorthand boundary. It always returns a network. Enabled `chain` and
`trace` roots materialize topology directly. Enabled compact
dot/Schoonschip syntax rewrites to ordinary parser input and recurses. A
disabled shorthand root becomes a fast-inferred leaf; an unchanged root
containing disabled compact shorthand also remains a leaf. An unchanged
ordinary function uses regular leaf parsing.

`SchoonschipMaterializer` is the narrower atom rewriter used by that
network boundary. It stores a `current` atom plus `additional_factors`.
Additional factors are multiplied beside the current atom and are not
inspected by the Schoonschip materializer during that same pass; they
are parsed when the parser recurses on the complete expression.

#link("../../crates/spenso/src/network/parsing/materialization.rs#L457")[crates/spenso/src/network/parsing/materialization.rs:457]

```
fn materialize_shorthand(...)
    -> Result<Self, TensorNetworkError<K, Symbol>>
{
    if symbol == SPENSO_TAG.chain && !root_chain_disabled {
        return Self::materialize_chain_shorthand(...);
    }
    if symbol == SPENSO_TAG.trace && !root_trace_disabled {
        return Self::materialize_trace_shorthand(...);
    }

    let materialized = SchoonschipMaterializer::with_mode(&state, mode)
        .materialize_shorthand(value.as_view());
    if materialized == value.as_view().to_owned() {
        if root_chain_disabled || root_trace_disabled || has_schoonschip_shorthand {
            return Self::as_inferred_leaf(..., StructureInferenceMode::Fast, ...);
        }
        return Self::parse_regular_function_leaf(...);
    }
    Self::try_from_view_impl(materialized.as_view(), ...)
}
```

#link("../../crates/spenso/src/network/parsing/materialization.rs#L114")[crates/spenso/src/network/parsing/materialization.rs:114]

```
pub(super) fn materialize_shorthand(&self, value: AtomView<'_>) -> Atom {
    self.materialize_shorthand_root(value)
        .map(SchoonschipMaterialization::into_expression)
        .unwrap_or_else(|| value.to_owned())
}
```

#link("../../crates/spenso/src/network/parsing/materialization.rs#L171")[crates/spenso/src/network/parsing/materialization.rs:171]

```
let rep = Self::compact_vector_rep(*lhs)?;
if rep != Self::compact_vector_rep(*rhs)? || !rep.rep.is_self_dual() {
    return None;
}

let slot = self.state.slot(&rep).to_atom();
```

#link("../../crates/spenso/src/network/parsing/materialization.rs#L20")[crates/spenso/src/network/parsing/materialization.rs:20]

```
// Chain and trace materialization chooses the symbolic in/out
// replacements, then lets this Schoonschip helper expand compact
// arguments inside each factor.
```

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
g(p(rep), q(rep))
-> p(rep(dummy)) * q(rep(dummy))
```

This currently requires both compact vectors to have the same self-dual
representation. `dot(p(rep), q(rep))` follows the same rule.

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

- is not `metric` or `dot`;
- is not itself a representation;
- has no explicit slot argument;
- has exactly one direct argument matching the representation wildcard
  convention.

Sums are accepted only when every summand is a compact vector with the
same representation. Products and powers are not compact vector syntax.

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

#link("../../crates/spenso/src/network/parsing/materialization.rs#L540")[crates/spenso/src/network/parsing/materialization.rs:540]

```
let start = Self::materialize_chain_endpoint(args[0], "start", ...)?;
let end = Self::materialize_chain_endpoint(args[1], "end", ...)?;
...
let factor = ChainExpansion::replace_placeholders(...);
let factor = SchoonschipMaterializer::with_mode(&state, factor_mode)
    .materialize_shorthand(factor.as_view());
```

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

#link("../../crates/spenso/src/network/parsing/materialization.rs#L701")[crates/spenso/src/network/parsing/materialization.rs:701]

```
let rep = Representation::try_from(*rep_view)?;
...
let links = (0..factors.len())
    .map(|_| state.slot(&rep))
    .collect::<Vec<_>>();
let left = links[position].to_atom();
let right = links[(position + 1) % factors.len()].dual().to_atom();
let factor = ChainExpansion::replace_placeholders(factor, &left, &right);
```

== Products, Sums, Powers
<products-sums-powers>
#table(
    columns: 2,
    align: (auto,auto,),
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

== Edge Cases
<edge-cases>
#table(
    columns: 2,
    align: (auto,auto,),
    [`F(p(rep), slot)`], [Must materialize before ordinary function
    parsing. Expected expansion is `F(dummy, slot) * p(dummy)`, leaving
    `slot` external and the dummy internal.],
    [`p(rep)`], [A standalone compact vector root is not materialized.
    It can only become a vector when used in a materializable argument
    or scalar-product position.],
    [`p(rep) * q(rep)`], [Not valid shorthand in the parser. The product
    path does not materialize this into a dot product.],
    [`dot(a,b,c)`], [Invalid parser syntax. Dot has two arguments; the
    parser rejects this before any scalar or opaque fallback.],
    [`g(slot, p(rep))`], [Materializes to `g(slot, dummy) * p(dummy)`.
    The contraction is represented by the network, not by a pre-parser
    algebra rewrite.],
    [`g(p(rep), q(rep))`], [For matching self-dual reps, materializes to
    `p(dummy) * q(dummy)`.],
    [Dualizable compact scalar products], [Current gap: compact
    scalar-product materialization requires self-duality. Dualizable
    cases produced by algebraic Schoonschip rules are not handled by
    this specific materializer rule.],
    [Compact vector endpoints in chains], [When Schoonschip expansion is
    enabled, compact endpoints materialize into fresh slots and
    additional rank-one factors. With that expansion disabled, endpoints
    must be concrete slots.],
    [Metadata false positives], [The compact-vector convention is
    syntactic. Any function argument with exactly one direct
    representation argument and no slots can be treated as a compact
    vector in expand mode.],
    [Opaque fast vs expanded inference], [Fast inference is cheap and
    syntactic. Expanded inference builds an expanded network and reads
    dangling slots, so it is a validation oracle but more expensive.],
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
```
parse(atom):
  dispatch by atom kind

parse_fun(fun):
  if fun is dot with arg count != 2:
    return InvalidDotFunction

  if !fun.is_tensorial(strict_tensor_filter):
    return scalar(fun)

  if fun is a broadcast wrapper:
    parse its single argument with broadcast handling

  if opaque mode and fun is recognized shorthand:
    infer structure
    return TensorFromExpression(original expression, inferred structure)

  return materialize_shorthand(fun)

materialize_shorthand(fun):
  if chain expansion is enabled and fun is chain:
    materialize compact endpoints when configured
    allocate links
    replace in/out placeholders
    materialize configured shorthand inside each factor
    parse factor expressions
    return chain network

  if trace expansion is enabled and fun is trace:
    return the dimension for zero factors
    allocate one cyclic link per factor
    replace in/out placeholders
    materialize configured shorthand inside each factor
    parse factor expressions
    return trace network

  select the effective Schoonschip expansion mode
  rewrite the function if that mode is enabled
  if an enabled shorthand remains unchanged:
    return an inferred tensor leaf
  if an ordinary function remains unchanged:
    return a regular tensor leaf
  parse the rebuilt expression recursively
```
