= Schoonschip network architecture
<schoonschip-network-parsing>
The network-backed Schoonschip path always parses tensorial shorthand
expressions as opaque tensor leaves. Parser shorthand expansion is
deliberately not part of this path.

#strong[Audit status:] reviewed 2026-08-17 against `c9f4e32acd2c`.
Lifecycle: current implementation architecture.

== Entry Points
<entry-points>
There are two Schoonschip implementations in
#link("../../crates/idenso/src/shorthands/schoonschip/api.rs")[crates/idenso/src/shorthands/schoonschip/api.rs];.
`schoonschip_with_settings` calls the older pattern-based implementation
in
#link("../../crates/idenso/src/shorthands/schoonschip/with_settings.rs")[with\_settings.rs];.
`schoonschip_with_net` calls the symbolic-network implementation.
`schoonschip_net` is a convenience wrapper for the default network pass
only. `schoonschip_with_net_full` is a separate convenience entry point
that enables contracted-sum expansion; neither wrapper subsequently
calls `schoonschip()` on the full result.

```
fn schoonschip_net<Aind: AbsInd + DummyAind + ParseableAind + 'static>(&self) -> Atom {
    self.schoonschip_with_net::<false, Aind>(&SchoonschipSettings::default_network())
}
```

== Network Algorithm
<network-algorithm>
#strong[Fixed Point] `NetworkSchoonschip::run` repeats until unchanged
unless the mode is `SinglePass`.
#strong[Term Split] `Add` nodes are processed term-by-term, then rebuilt
as a sum.
#strong[Opaque Network] Each term is dot-normalized and parsed into a
symbolic tensor network with opaque tensorial shorthand.
#strong[Contract] The selected contraction order chooses graph-adjacent
tensors and contracts them through `Schoonschipify`.
#strong[Normalize] The result expression is extracted, dot-normalized,
and fed back into the fixed point loop.
The parser settings are built in
#link("../../crates/idenso/src/shorthands/schoonschip/api.rs")[NetworkSchoonschip::run\_once];.
The important convention is that tensorial shorthand is never expanded
by the parser in this path:

```
shorthand_parsing: ShorthandParsing::Opaque {
    inference: StructureInferenceMode::Fast,
},
```

=== Driver outline
<driver-outline>
```
current = input
loop:
    next = apply(current)
    if settings.mode == SinglePass || next == current:
        return next
    current = next

apply(expr):
    if expr is a sum:
        return sum(apply(term) for term in expr)

    normalized = normalize_dots(expr)
    net = parse_to_symbolic_net(normalized, opaque_shorthand, depth_limit)
    execute net with selected Schoonschip contraction strategy
    return net.result_expression().normalize_dots()
```

Pure scalar atoms still use the parser\'s scalar override. Tensorial
shorthand such as dot products, traces, chains, and Schoonschip vector
slots becomes a symbolic tensor leaf with inferred structure. That keeps
graph construction structural; algebraic meaning is assigned later by
the contraction code in
#link("../../crates/idenso/src/shorthands/schoonschip/contraction.rs")[contraction.rs];.

== Settings
<settings>
#table(
    columns: 2,
    align: (auto,auto,),
    table.header([*Setting*], [*Behavior*]),
    [`mode`], [Selects `SinglePass` or a recursive fixed point.
    Recursive mode is either depth-first or breadth-first; the const
    parameters `RECURSE` and `DEPTH_FIRST` are chosen from this.],
    [`depth_limit`], [Passed to parsing. It limits how deeply composite
    expressions are turned into network structure.],
    [`contraction_order`], [Chooses the graph edge selection strategy
    used during execution.],
    [`expand_contracted_sums`], [Enables sum-by-sum contraction
    expansion inside `Schoonschipify`. It does not change parser
    shorthand mode.],
    [`simplify_chain_like_functions`], [Controls the optional chain-like
    metric simplifier in the pattern implementation. The network
    implementation does not read this field.],
    [`schoonschip_rank1_tensors`], [Enables rank-one rewrite rules in
    the pattern implementation. The network implementation does not read
    this field.],
)

== Contraction Step
<contraction-step>
#table(
    columns: 2,
    align: (auto,auto,),
    table.header([*Step*], [*Behavior*]),
    [Scalar pre-pass], [The selected contraction strategy builds a
    `ProductContraction`, simplifies composite scalar tensors, and calls
    `product.contract_scalars` before, between, and after tensor-pair
    contractions. Recursive modes pass composite scalar expressions
    through the same network algorithm.],
    [Edge choice], [The strategy is selected by
    `SchoonschipContractionOrder`: smallest degree, largest degree, or
    expression-size heuristics such as minimum largest operand bytes,
    product terms, or product bytes.],
    [Structure merge], [The two tensor structures are merged first. The
    merge reports which slots were contracted on each side, so the
    algebra uses graph edge information rather than guessing from
    function argument shapes.],
    [Special contractions], [Metrics rewrite the target slot with the
    metric\'s free slot. Rank-one tensors strip their consumed
    representation slot and are inserted into the contracted target
    slot.],
    [General product], [If no shortcut applies, the expressions are
    multiplied and wrapped as a composite tensor. For genuine contracted
    sum-by-sum products, the expanded mode first tries a direct
    contraction of the smaller expanded side with slot-aware residual
    cleanup. It falls back to distribution followed by a network
    Schoonschip pass. This is contraction behavior, not parser shorthand
    expansion.],
    [Finish], [`finish_contract` can recurse on the result in
    breadth-first mode, then dot-normalizes the expression before
    returning the new symbolic tensor.],
)

== Pattern Path
<pattern-path>
The pattern implementation in
#link("../../crates/idenso/src/shorthands/schoonschip/with_settings.rs")[with\_settings.rs]
is a fixed point over global Symbolica replacements. One pass applies:
dot normalization, repeated metric-function contractions, optional
rank-one vector and vector-dot rules, optional chain/trace rules, and a
final dot normalization.

```
let metric_simplified = view
    .normalize_dots()
    .replace_multiple_repeat(&*METRIC_FUNCTION_CONTRACTIONS);

let simplified = if settings.schoonschip_rank1_tensors {
    metric_simplified
        .replace_multiple(&*VECTOR_DOT_PRODUCTS)
        .replace_multiple_repeat(&*SCHOONSCHIP_VECTOR)
} else {
    metric_simplified
};
```

That path is still useful for local cleanup and for simple expressions.
It is not a good main engine for contraction-heavy tensor expressions,
because it repeatedly scans the whole atom and must encode structural
facts as wildcard patterns and conditions.

== Why The Network Path Should Be Better
<why-the-network-path-should-be-better>
=== Structure before syntax
<structure-before-syntax>
The network path knows which tensor slots are connected. A contraction
such as `g(mu, nu) * F(..., nu, ...)` is driven by the merged edge
positions, not by matching a function with a slot-looking second
argument.

=== Controlled order
<controlled-order>
Pattern replacement has no global contraction-order model. The network
path can choose the next edge by degree or expression-size estimates,
which is expected to reduce intermediate expression growth.

=== Generic indices
<generic-indices>
Slots carry `Aind` through the graph. This is less brittle than relying
on concrete Symbolica patterns when GammaLoop-owned abstract index types
are valid indices.

=== Less double meaning
<less-double-meaning>
Opaque parsing prevents parser materialization from creating dummy
indices or expanded subgraphs before Schoonschip chooses contractions.
Shorthand is represented once, then simplified by the network rules.

== Examples
<examples>
#table(
    columns: 2,
    align: (auto,auto,),
    table.header([*Example*], [*Result*]),
    [Metric into tensor], [`g(mu, nu) * F(a, nu, b)` contracts through
    the graph edge on `nu`. The metric\'s remaining slot `mu` is
    substituted into the tensor expression, yielding `F(a, mu, b)`. The
    pattern path needs explicit wildcard forms for the same idea.],
    [Rank-one into tensor], [`p(rep(d, i)) * F(a, rep(d, i), b)` becomes
    `F(a, p(rep(d)), b)`. The consumed vector slot is known from the
    structure merge, so the rule does not have to inspect unrelated
    arguments.],
    [Dot product], [In `p(mu) * q(mu)`, each vector parses as a rank-one
    leaf and the shared slot creates the contraction edge. The rank-one
    shortcut inserts one stripped vector into the other; the finishing
    dot normalization then rewrites that nested form to a compact scalar
    such as `g(p(rep), q(rep))`.],
    [Chain or trace shorthand], [`chain(...)` and `trace(...)` are
    tensorial shorthand. In the network path they parse as opaque leaves
    with inferred external structure. Chain-like pattern rules are
    available only through the pattern implementation; the corresponding
    settings field is ignored by the network path.],
    [Composite scalar recursion], [A tensor product can contract to a
    scalar while still carrying a composite expression. Recursive modes
    pass that scalar expression back through `schoonschip_with_net`,
    using the same opaque parser, until the configured traversal reaches
    a fixed point.],
    [Contracted sums], [With contracted-sum expansion enabled, genuine
    sum-by-sum contractions such as `(A(mu)+B(mu)) * (C(mu)+D(mu))` are
    expanded inside the contraction step. One-sided sums such as
    `g(mu, nu) * (A(nu)+B(nu))` are left to the metric shortcut instead
    of being blindly distributed.],
)

== Expected Outcome
<expected-outcome>
The network path should be the default for expressions where contraction
topology matters: many tensors, mixed metrics, rank-one substitutions,
composite scalar products, and sums where expansion order can dominate
the cost. The pattern path remains valuable as a simple local rewrite
pass and can be selected explicitly with `schoonschip_with_settings`;
`schoonschip_net` does not append it as a whole-expression final pass.
