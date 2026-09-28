= Symbolica atom-access performance reproducers

The three Cargo targets depend only on Symbolica 3.0.0 and the Rust standard
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

== Global initializer completion

`initialization_lifecycle.rs` is a separate, untimed lifecycle regression
executable. It requires the two patches in `patches/` and matching Symbolica and
Spenso libraries. It is intentionally not a Cargo target here: this package's
unchanged pinned dependency does not expose the proposed completion query.

The existing Spenso lazy-symbol handle probes Symbolica initialization before
taking a bundle's lazy lock. That cold ordering prevents recursive initialization
deadlocks. Repeating the probe after initialization also takes a state read lock
and hashes a dummy builtin name. The proposed Symbolica query reads the existing
initializer's completion state; the Spenso patch uses it to skip completed probes.
It adds no separate readiness cache. A reentrant initializer callback can return
from a probe before other callbacks finish, so successful probing alone cannot
certify global completion.

Apply the upstream query patch independently to an isolated Symbolica checkout
at `06906976bca24fefc5203aee699d90d62ebe08cd`:

```sh
git -C "$SYMBOLICA_CHECKOUT" apply --check \
  "$REPRODUCER_DIR/patches/symbolica-state-completion.patch"
git -C "$SYMBOLICA_CHECKOUT" apply \
  "$REPRODUCER_DIR/patches/symbolica-state-completion.patch"
```

`REPRODUCER_DIR` is the absolute path to this directory. The query returns false
throughout registered callbacks, including reentrant calls, and after an
initializer panic. It does not initiate initialization. Unsafe symbol-state reset
does not rerun the global callbacks or clear their completion status.

For the Spenso comparison, use two isolated copies of the same FeynKit source;
apply `patches/spenso-initialization-gate.patch` only to the candidate. Build both
against that same query-patched Symbolica. Rebuild transitive consumers such as
Symbolica-utils and Linnet against it too: mixing two Symbolica artifacts can
create distinct global states and invalidates the test. Preserve the same
features, optimization settings and remaining dependencies on both sides.

Link the lifecycle source with the resulting matching libraries. For a direct
`rustc` build, provide the complete dependency search directories and native link
configuration from that family:

```sh
CARGO_CRATE_NAME=symbolica_init_lifecycle rustc --edition=2024 \
  --crate-name symbolica_init_lifecycle \
  -L "dependency=$MATCHED_DEPS" \
  --extern "symbolica=$MATCHED_SYMBOLICA_RLIB" \
  --extern "spenso=$MATCHED_SPENSO_RLIB" \
  initialization_lifecycle.rs -o initialization_lifecycle

for mode in bundle-first state-first concurrent initializer-panic \
            bundle-panic marker-unwind; do
  timeout 20s ./initialization_lifecycle "$mode"
done
```

Repeat for the unchanged Spenso control. Every mode runs in a fresh process;
`State::reset` cannot reset the initializer's one-time state. The callback
deliberately accesses an already warmed Spenso bundle without Spenso's local
initializer marker. The concurrent case holds that callback until another
thread attempts bundle access, then checks that the reader waits for completion.
Caught-panic diagnostics are expected in the panic and unwind modes. Also retain
the existing Spenso `symbolica_lazy_deadlock_mwe --safe` check.

The isolated checkpoint in
`../../notebooks/tensor_contraction_parity.json`, under
`symbolica_initialization_gate_checkpoint`, records 14 passing lifecycle/MWE
runs, 3,712 exact network output/error comparisons, 61 metric controls, 11 gamma
fixtures with rerun checks on both sides, and four structure/compute controls.
It includes matching build identities, raw paired timings and instruction
profiles. The copied Idenso source predates the refined scalar-power admission
change; these are gate-only measurements and cannot be added to that separate
optimization's savings.

Across three paired process rounds, the three-vertex gluon network improved from
11.96 to 10.74 ms for smallest-degree order and from 12.55 to 11.22 ms for
minimum-product-terms order. The corresponding instruction counts fell about
14.5%; the removed builtin-probe work accounts for nearly all that reduction.
The untouched scalar-expansion control used exactly the same instruction count
and had a paired timing ratio of 0.989–1.035. Short metric chains gained little,
and some small controls and terminal gamma timings remained noisy. These results
establish the initialization owner as useful work to remove; they do not show
FORM parity or combined current-pipeline performance. The dependency pins and
production initialization code remain unchanged.
