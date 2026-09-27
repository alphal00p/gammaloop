= Symbolica alias and streaming requests
<symbolica-alias-and-streaming-requests>
#quote(block: true)[
#strong[Lifecycle:] standalone reproducers for three Symbolica requests raised by
#link("../../../docs/architecture/tensor-engine-consolidation-plan.typ")[the tensor engine consolidation plan]
(decision D4, register item R18). Each script prints the limitation it
demonstrates and the request it motivates. Recorded output is from
2026-09-27 against Symbolica `06906976` (Rust) and the community build 3.0.0
(Python), on a shared host; timings are indicative.
]

== Running

Each Rust file is a `rust-script` with its dependency block in the header and
is also a binary of the package in this directory, which shares one build:

```bash
rust-script parametric_alias.rs
cargo run --release --bin parametric_alias
cargo run --release --bin streaming_expansion -- 8 6
cargo run --release --bin aliased_evaluator -- 18
python aliased_evaluator_python.py 18
```

The Python script needs an interpreter with the Symbolica community build
and NumPy installed; numerical evaluator output uses NumPy arrays.

M0 revalidation on 2026-09-27 ran all four examples successfully against
the preserved `cfe4e655` host and source-matching release Rust executables.
#link("qualification.json")[The qualification record] retains executable and
source hashes, commands, exact assertions and raw output. Its clocks are
diagnostics, not engine performance comparisons. NumPy was installed in a
separate reproducer environment, leaving the benchmark environments unchanged.

== Request 1: parametric aliases
<parametric-aliases>

`parametric_alias.rs`. An `AliasedAtom` handle that carries slot arguments,
`T(1, mink(4,mu), mink(4,nu))`, stands for a tensor-valued definition with two
free ports. Contraction rewrites the ports of such a handle without opening it,
so the same alias is used as `T(1, mink(4,rho), mink(4,nu))` after
`g(mink(4,rho), mink(4,mu))` acted on it. `AliasedAtom::into_inner` resolves
handles by literal atom equality, so the relabeled use stays unresolved. A
one-line pattern replacement (`T(1, mink(4,a_), mink(4,b_))` with a body in
`a_`, `b_`) resolves every use, which is what the alias table should be able
to store. The available workaround, one literal alias per labelling, repeats
the definition per use.

Request: `AliasedAtom::register_pattern_alias(pattern, body)`, honoured by
`into_inner`, `get_aliases` and `EvaluatorBuilder::add_aliases`, so a handle
with arguments resolves like a function of its ports.

Recorded output:

```text
root:                 T(1,mink(4,mu),mink(4,nu))+T(1,mink(4,rho),mink(4,nu))
alias:                T(1,mink(4,mu),mink(4,nu)) = m*g(mink(4,mu),mink(4,nu))+p(mink(4,mu))*q(mink(4,nu))
into_inner():         m*g(mink(4,mu),mink(4,nu))+p(mink(4,mu))*q(mink(4,nu))+T(1,mink(4,rho),mink(4,nu))
relabeled handle left unresolved by literal lookup: true
pattern resolution:   m*g(mink(4,mu),mink(4,nu))+m*g(mink(4,rho),mink(4,nu))+p(mink(4,mu))*q(mink(4,nu))+p(mink(4,rho))*q(mink(4,nu))
handles left after pattern resolution: false
bytes: one pattern alias 182 | 8 literal aliases 966 (5x the definition)
per-labelling into_inner() resolves the relabeled use: true

ask: AliasedAtom::register_pattern_alias(pattern, body) used by into_inner,
     get_aliases and EvaluatorBuilder::add_aliases, so a handle with slot
     arguments resolves like a function of its ports.
```

== Request 2: streaming expansion and replacement
<streaming-expansion>

`streaming_expansion.rs`. `expand()` is the only way to see the terms a
product of sums generates, and it materializes all of them into one normalized
atom. A contraction engine needs to visit each generated term as borrowed
factor views, contract it and merge it into a small table; the expanded atom is
never wanted. The script times `expand()` on a product of eight six-term sums
with distinct variables (1,679,616 terms, none merging) against a walk over the
factor views that visits every generated term and builds no atom. Idenso's
collector re-implements that walk (its tape and distributor) because Symbolica
has no streaming primitive.

Request: an `expand`, and a `replace`, that hand each generated term to a
consumer as borrowed factor views with Symbolica's coefficient and power
folding applied, without assembling the sum.

Recorded output:

```text
product of 8 sums with 6 terms: 1679616 generated terms
expand():      6.179 s, 1679616 terms, 52.3 MB
walk views:    0.034 s, 1679616 terms visited, coefficient checksum 1, 0 atoms built
materialization costs 180x the traversal

ask: an `expand`/`replace` that hands each generated term to a consumer as
     borrowed factor views (with Symbolica's coefficient and power folding),
     so an engine can contract and merge terms without the expanded atom.
```

== Request 3: evaluators from aliased atoms in Python
<aliased-evaluators>

`aliased_evaluator.rs` and `aliased_evaluator_python.py`. A contraction result
is a DAG: every node is an alias whose definition refers to earlier aliases.
In Rust, `AliasedAtom::evaluator_multiple` builds an evaluator straight from
the alias table. The only alternative is to resolve the aliases into one tree
first, which duplicates every shared node; the scripts use a chain whose
definitions use their predecessor twice, so the tree doubles per level while
the alias table grows by one line. Python exposes only the tree route:
`Expression.evaluator` has no alias input and there is no aliased expression
type, which is what `gammalooprs` would have to use to consume a contraction
result without expanding it.

Request: an `AliasedExpression` (root plus `(alias, body)` pairs) with
`evaluator(...)` in the Python bindings, or an `aliases=` input on
`Expression.evaluator`, mirroring the Rust entry point.

Recorded output (Rust):

```text
alias table: 19 definitions, 883 bytes
evaluator from aliases:  built in 0.002 s, operations OperationCount { additions: 18, multiplications: 54, inversions: 0, function_calls: 0 }
into_inner():            1.054 s, 7.1 MB tree (46 bytes per definition otherwise)
evaluator from the tree: built in 0.354 s, operations OperationCount { additions: 18, multiplications: 54, inversions: 0, function_calls: 0 }
same value: true (3.3225074082222934e-6)

ask: expose the alias route in Python — an `AliasedExpression` (root + alias
     pairs) with `evaluator(...)`, or an `aliases=` input on `Expression.evaluator`,
     so results are consumed as DAGs without resolving them to trees.
```

Recorded output (Python):

```text
symbolica 3.0.0
Expression.evaluator parameters: self, params, functions, iterations, cpe_iterations, n_cores, verbose, jit_compile, direct_translation, jit_direct_translation, jit_optimization_level, jit_options, max_horner_scheme_variables, max_common_pair_cache_entries, max_common_pair_distance
alias input available: False
aliased expression type available: False
resolved tree for depth 18: 0.048 s, 7.1 MB
evaluator from the tree: built in 0.133 s
value at x = 0.5: [[3.32250741e-06]]

ask: an AliasedExpression (root + (alias, body) pairs) with .evaluator(...) in the
     Python bindings, mirroring AliasedAtom::evaluator_multiple in Rust.
```
