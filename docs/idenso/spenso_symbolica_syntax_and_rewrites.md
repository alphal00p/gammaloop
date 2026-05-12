# Spenso Symbolica Syntax And Rewrite Idioms

This note describes the Symbolica syntax Spenso uses for tensor expressions and
the Rust idioms for implementing algebraic rewrites over that syntax. It is
meant as a general reference for tensor, metric, Schoonschip, Dirac, and color
simplifiers; the FORM-derived identities in
`docs/idenso/form_symbolica_color_and_dirac.typ` are examples of this style, not
the definition of the syntax.

## Tensor Syntax

Spenso tensor expressions are Symbolica atoms. The names in the surface syntax
mirror the Rust model:

- a representation kind implements `RepName`;
- a concrete representation is a `Representation<LibraryRep>`;
- an indexed tensor slot is a `Slot<LibraryRep, Aind>`;
- a symbolic tensor leaf carries a structure made from those slots.

In Symbolica syntax, a tensor slot is a function whose head is the
representation name. The first argument is always the dimension, and the second
argument is the abstract index:

```text
mink(D, mu)
```

Here `mink` is the representation head, `D` is the dimension, and `mu` is the
abstract index. The same representation without its index is the stripped
representation:

```text
mink(D)
```

Tensor leaves are functions whose head is the tensor symbol and whose structural
arguments are representation slots.

```text
F(mink(D, mu), cof(NC, i), dind(cof(NC, j)))
```

The standard slot spellings are:

```text
mink(D)          stripped Minkowski representation
mink(D, mu)      indexed Minkowski slot
bis(D, i)        indexed self-dual bispinor slot
cof(NC, i)       indexed color-fundamental slot
dind(cof(NC, j)) dual color-fundamental slot
```

`dind(...)` is the dual wrapper from `AIND_SYMBOLS`. Self-dual representations
such as `mink`, `bis`, and `coad` use the same representation on both sides of a
contraction. Dualizable representations such as `cof` use the unwrapped slot for
one orientation and `dind(...)` for the dual orientation.

Direct slot arguments define the exposed tensor structure. Non-slot arguments
are metadata on the symbolic tensor leaf:

```text
P(1, mink(D, mu))     rank-one tensor with label 1 and one Lorentz slot
T(a, b, cof(NC, i))   symbolic tensor with scalar labels a, b and one slot
```

The same convention is used in patterns: structural arguments are represented by
slot or representation patterns, and scalar labels are ordinary Symbolica atoms.

## Supported Surface Forms

The parser and simplifiers recognize the following Spenso surface forms. The
same row shows how to build the form, which tag or symbol makes it recognizable,
and what to use when matching it in a replacement rule.

| Surface form | Meaning | Concrete builder | Tag or symbol | Pattern form |
| --- | --- | --- | --- | --- |
| `f(args...)` with no representation-tagged syntax | pure scalar Symbolica expression | `function!(f, args...)` | no Spenso tensor tag | ordinary Symbolica pattern |
| `pure_scalar(expr)` | force `expr` to parse as a scalar | `function!(T.pure_scalar, expr)` | `SPENSO_TAG.pure_scalar` symbol | `function!(T.pure_scalar, W_.x_)` |
| `rep(dim)` | stripped representation, for compact notation or traces | `rep.to_symbolic([])` or `IntoAtom` on `Representation` | representation head tagged `representation` | `T.rep_::<N, _>([W_.d_])` |
| `rep(dim, i)` | indexed representation slot | `slot!(rep, i)` | representation head tagged `representation`; maybe `self_dual` or `dualizable` | `T.rep_::<N, _>([W_.d_, W_.i_])` |
| `dind(rep(dim, i))` | dual slot | `slot.dual().to_atom()` or explicit `AIND_SYMBOLS.dind` | `dind` dual wrapper plus representation tags | `T.dualizable_dual_::<N, _>([W_.d_, W_.i_])` |
| `aind(slot...)` | bundle of structural slots inside one argument | `function!(AIND_SYMBOLS.aind, slots...)` | `AIND_SYMBOLS.aind` symbol | `function!(AIND_SYMBOLS.aind, W_.x___)` |
| `F(..., slot, ...)` | ordinary tensor leaf; direct slots define structure | `tensor!(F, args...)` | head tagged `tensor` | `T.tensor_::<N, _>([W_.a___])` or `function!(W_.f_, args...)` |
| `p(..., rep(dim))` | compact rank-one tensor shorthand | `vector!(p, args...)`, `p!(args...)`, `q!(args...)` | head tagged `tensor` and `rank1` | `T.rank1_::<N, _>([W_.a___, T.rep_::<M, _>([W_.d_])])` |
| `g(slot_a, slot_b)` | metric tensor syntax | `function!(ETS.metric, a, b)` | `ETS.metric` symbol | `function!(ETS.metric, a, b)` |
| `g(p(rep), q(rep))` | compact scalar product shorthand | `function!(ETS.metric, p, q)` | `ETS.metric` plus rank-one compact arguments | `function!(ETS.metric, T.rank1_::<0, _>(...), T.rank1_::<1, _>(...))` |
| `dot(p(rep), q(rep))` | two-argument compact dot shorthand | `function!(T.dot, p, q)` | `SPENSO_TAG.dot` symbol | `function!(T.dot, a, b)` |
| `chain(start, end, factors...)` | ordered open chain | `chain!(start, end, factors...)` | `SPENSO_TAG.chain` symbol; ordered arguments | `chain!(start, end, factors...)` |
| `trace(rep, factors...)` | ordered closed trace | `trace!(rep, factors...)` | `SPENSO_TAG.trace` symbol; ordered arguments | `trace!(rep, factors...)` |
| `bracket(expr...)` | product-like parser grouping of network factors | `function!(T.bracket, expr...)` | `SPENSO_TAG.bracket` symbol | `function!(T.bracket, W_.x___)` |
| `broadcast(expr)` | apply a scalar/broadcast function to the parsed inner tensor | `broadcast_symbol!(f)` then `function!(f, expr)` | head tagged `broadcast` | `function!(W_.f_, W_.x_)` with a broadcast-head condition |
| `sum`, `product`, `integer power` | ordinary Symbolica arithmetic, parsed recursively | `+`, `*`, `.pow(...)` | Symbolica arithmetic nodes | ordinary Symbolica pattern, with conditions for exponent cases |

`dot` is a two-argument shorthand. A three-argument spelling such as
`dot(rep, p, q)` is not parser syntax.

`chain` and `trace` are ordered and non-symmetric. Their factor order is part of
the algebra. A factor inside a chain uses the symbols `in` and `out` as local
placeholder slots:

```text
chain(bis(D, i), bis(D, j),
  gamma(in, out, mink(D, mu)),
  gamma(in, out, mink(D, nu)))

trace(bis(D),
  gamma(in, out, mink(D, mu)),
  gamma(in, out, mink(D, nu)))
```

In expanded parsing, chain materialization replaces `in` and `out` with actual
slots and creates intermediate dummy slots between adjacent factors. In opaque
parsing, the chain or trace remains a leaf and its exposed structure is inferred
from the endpoints and visible external slots.

Compact Schoonschip-style rank-one notation is local shorthand:

```text
F(..., p(mink(D)), ...) -> F(..., mink(D, d), ...) * p(mink(D, d))
g(p(mink(D)), q(mink(D))) -> p(mink(D, d)) * q(mink(D, d))
dot(p(mink(D)), q(mink(D))) -> p(mink(D, d)) * q(mink(D, d))
```

The dummy `d` is fresh and local to the materialized expansion.

## Tagged Builders And Patterns

Concrete expression builders, parser detection, and replacement patterns are the
three uses of the same `SPENSO_TAG` mechanism shown in the table above. Spenso
does not identify tensor syntax by string names alone. A symbol for `mink`
carries the `representation` tag; a self-dual representation also carries
`self_dual`; a dualizable representation also carries `dualizable`. Tensor heads
such as `F`, `gamma`, and `t` carry `tensor`, while vector heads such as `p`,
`q`, `u`, and `v` carry both `tensor` and `rank1`. Broadcast heads carry
`broadcast`, index symbols may carry `index`, and `upper`/`lower` are available
for oriented index syntax.

The expression macros are the concrete builders for the rows in the surface-form
table. `slot!(rep, mu)` builds a concrete `Slot<LibraryRep, Aind>` and prints as
`rep(dim, mu)`. `tensor!(F, args...)` creates a tensor-tagged head.
`vector!(p, args...)`, and the short `p!(...)` and `q!(...)` forms, create
rank-one tensor heads. `chain!` and `trace!` build the ordered containers,
including iterator forms `chain!(start, end; factors)` and
`trace!(rep; factors)`. `sym!` and `antisym!` build inert ordered projectors.
`function!(head, args...)` remains the fallback when there is no Spenso-specific
builder, and `s!(mu)` is the lightweight plain-symbol helper for labels and
indices.

When adding new heads, use the symbol-constructor macro for the semantic class
you need: `tensor_symbol!(F)`, `vector_symbol!(p)`,
`representation_symbol!(rep)`, `self_dual_symbol!(mink)`,
`dualizable_symbol!(cof)`, `index_symbol!(i)`, or
`broadcast_symbol!(sqrt)`. These macros are the concrete side of the same tag
system used by the pattern constructors.

For concrete expressions in tests and examples, prefer the expression macros:

```rust
let mu = slot!(mink4, mu);
let expr = chain!(
    slot!(bis4, i),
    slot!(bis4, j),
    gamma!(mu),
    gamma!(slot!(mink4, nu)),
);
let dot = function!(T.dot, p!(&mink4), q!(&mink4));
```

For replacement rules, use wildcard heads with the same tags instead of naming
one concrete symbol. `T.rep_::<N, _>(args)` matches any representation syntax;
`T.self_dual_::<N, _>(args)` restricts that to self-dual representations;
`T.dualizable_::<N, _>(args)` and `T.dualizable_dual_::<N, _>(args)` match the
two orientations of dualizable representations; `T.tensor_::<N, _>(args)`
matches any tensor-tagged head; and `T.rank1_::<N, _>(args)` matches any
rank-one tensor head. Numbered constructors create independent wildcard heads,
so `rank1_::<0, _>` and `rank1_::<1, _>` can match two different vector symbols
in the same rule.

```rust
let slot = T.self_dual_::<0, _>([W_.d_, W_.i_]);
let stripped = T.self_dual_::<0, _>([W_.d_]);
let lhs = T.rank1_::<0, _>([W_.a___, &slot])
    * T.rank1_::<1, _>([W_.b___, &slot]);
let rhs = function!(
    T.dot,
    T.rank1_::<0, _>([W_.a___, &stripped]),
    T.rank1_::<1, _>([W_.b___, &stripped]),
);
```

Use `function!(W_.f_, args...)` when the head itself is arbitrary Symbolica
syntax rather than a tagged tensor class. Use `function!(ETS.metric, a, b)` for
metric rules and `function!(T.dot, a, b)` for compact dot rules. `chain!` and
`trace!` are also the preferred pattern builders for ordered containers.

`W_` contains the plain Symbolica wildcards used inside these constructors:
`x_` matches one atom, `x__` matches one or more atoms, and `x___` matches zero
or more atoms. Repeating a wildcard name enforces equality; different names mean
independent matches. Add conditions when tags and wildcard shape are still too
broad, for example `filter_single` for even/odd powers, `filter_match` for
structural predicates, or local predicates such as `not_slot(...)` to prevent
variadic tails from swallowing structural slots.

idenso adds domain-specific expression macros for common Dirac objects:
`gamma!`, `gamma0!`, `gamma5!`, `u!`, and `v!`. Use these in tests and examples
when they make the expression read like the algebra. In generic rewrite rules,
still match the structural class with Spenso tag patterns when the rule is not
specific to one concrete head.

## Rust Rewrite Idioms

Prefer a small pass type over a pile of free helpers:

```rust
struct DotNormalizer;

impl DotNormalizer {
    fn run(view: AtomView<'_>) -> Atom {
        view.to_owned()
            .replace_multiple_repeat(&VECTOR_RULES[..])
            .replace_multiple_repeat(&METRIC_RULES[..])
    }

    fn even_power(view: AtomView<'_>) -> bool {
        // predicate body
    }
}
```

For settings-driven simplifiers, keep the settings on the pass object and make
the pass order explicit:

```rust
struct AlgebraSimplifier<'a> {
    settings: &'a AlgebraSettings,
}

impl AlgebraSimplifier<'_> {
    fn run(&self, view: AtomView<'_>) -> Atom {
        let mut current = view.to_owned();
        loop {
            let next = self.apply_once(current.as_view());
            if next == current {
                return next;
            }
            current = next;
        }
    }

    fn apply_once(&self, view: AtomView<'_>) -> Atom {
        // ordered passes
        view.to_owned()
    }
}
```

Use private static replacement tables for local, static identities:

```rust
static VECTOR_DOT_PRODUCTS: LazyLock<[Replacement; 1]> = LazyLock::new(|| {
    let slot = T.self_dual_::<0, _>([W_.d_, W_.i_]);
    let stripped = T.self_dual_::<0, _>([W_.d_]);

    [Replacement::new(
        (T.rank1_::<0, _>([W_.a___, &slot])
            * T.rank1_::<1, _>([W_.b___, &slot]))
        .to_pattern(),
        function!(
            T.dot,
            T.rank1_::<0, _>([W_.a___, &stripped]),
            T.rank1_::<1, _>([W_.b___, &stripped]),
        ),
    )]
});
```

This style is appropriate when the right-hand side is determined entirely by the
matched wildcard substitution.

Use `replace_map` or `ReplaceBuilder::with_map` when the right-hand side must be
computed from the match:

- reversing or rotating ordered chain factors;
- flipping `in` and `out` in reversed chain factors;
- Chisholm rules depending on the parity or length of the matched middle;
- recursive trace expansion;
- cyclic color trace normalization;
- loop-detection algorithms such as FORM `ReplaceLoop`;
- fresh dummy allocation.

Disable RHS caching when a map creates match-local fresh symbols:

```rust
let settings = MatchSettings {
    rhs_cache_size: 0,
    ..Default::default()
};
```

Implement FORM-like simplifiers as ordered passes, not as one unordered
`replace_multiple` table. Put cheap terminal rules first, then structural
normalization, then broader recursive searches. If a pass can expose another
rule in the same family, run to a fixed point with an equality guard.

Keep semantic boundaries clear:

- `chain` and `trace` stay opaque ordered containers except inside chain/trace
  materialization or their dedicated algebraic simplifiers.
- Four-dimensional identities must match a four-dimensional representation in
  the pattern or check it in the map.
- Fresh dummies belong at parser materialization or explicit dynamic rewrite
  boundaries, not in static replacement tables.
- Static rule tables should be private; expose only the small module entry point
  needed by callers.
- Unit tests for narrow rule families should live next to the implementation.
  Larger FORM/FeynCalc/reference cases should stay in the broader test modules.

The goal is source-backed named rule families plus a small ordered pass driver.
Avoid a universal rewrite engine unless the algebraic strategy itself is shared.
