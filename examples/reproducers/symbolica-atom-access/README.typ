= Symbolica atom-access performance reproducers

These standalone programs depend only on Symbolica 3.0.0 and the Rust standard
library. They construct their inputs at startup, assert that both operations
return the same result, and exclude parsing and validation from the timings.
The optimized loops consume their results with `black_box` and alternate method
order. `Cargo.lock` pins the dependency tree.

== Run

From this directory:

```sh
cargo build --release
./target/release/iterator > iterator.jsonl
./target/release/header > header.jsonl
./target/release/last_argument > last_argument.csv
```

On Linux, optionally prefix each executable with `taskset -c N`, selecting an
available CPU. Repeat the executable in fresh processes. The release profile uses
optimization level 2 with LTO disabled; change `opt-level` to 3 to check sensitivity
to compilation. Results are timings, not performance assertions with fixed
thresholds. The Symbolica integer and float backends use GMP and MPFR.

== Second argument through nth

In `iterator.rs`, the essential comparison is:

```rust
arguments.nth(1).unwrap()
```

against:

```rust
arguments.next();
arguments.next().unwrap()
```

Both consume identical copies of an opaque `ListIterator` from `f(4,mu)`.
Symbolic, rational, function, and sum first arguments exercise additional atom
shapes. Both inline-eligible readers and equally outlined readers are measured.
This exposes a code-generation difference in equivalent public operations. The
standalone assembly does not contain the out-of-line `try_fold` call observed in
the earlier full scanner, so that earlier explanation does not transfer directly.

== Function header access

In `header.rs`, both operations return exactly the same head ID and arity:

```rust
(fun.get_symbol_id(), fun.iter().len())
(fun.get_symbol_id(), fun.get_nargs())
```

The cases include a short header, mixed heads and arities, a 300-argument header,
and a wider head ID. The direct accessors can share their header decode after
inlining; the iterator path can retain a separate `FunView::iter` call. This
reproduces header-access overhead. It does not show that a consumer needing the
arguments can remove its iterator, or measure a hypothetical fused accessor.

== Final-argument boundary

In `last_argument.rs`, both readers check for exactly two arguments and decode
the first argument. The ordinary reader then calls `next()` again. The comparison
borrows the exact remaining function bytes as the final argument, using only
public `get_data` and `AtomView::from` methods. It assumes a valid exact `FunView`;
it is not a parser or validator for arbitrary serialized bytes. No packed-format
constants, private readers, or cached argument offsets are used.

For `f(4,mu)` this changes only a small constant cost. A nested power as the final
argument exposes scaling: the ordinary iterator walks its encoding to locate an
end that the function boundary already determines. Input construction is outside
timing. Nested powers are a deliberate boundary case; axial-trace slots normally
contain simple symbol indices, so the large speed ratio is not a prediction for
the complete axial scanner.

== Measurements and limits

`measurements.json` records the independent package build and timing samples.
These are three separate reproducers, not additive portions of the gamma
simplifier runtime. All operations borrow immutable atoms and return identical
views; no algebraic rewrite or tensor-network construction occurs in their loops.
