# Standalone Symbolica evaluator-construction diagnostic

This package has one direct dependency, Symbolica. It imports a supplied scalar
expression and its exact ordered parameters/function definitions, then calls
`expression.evaluator(&params).function_map(map).optimization_settings(settings).build()`.
There are no GammaLoop or Spenso dependencies, tensor processing, numerator
expansion, C++ export, generated-code compilation, or JIT calls. Cargo builds
the Rust diagnostic executable; it does not compile the generated evaluator.
The enabled `native_code_generation` library feature matches the original
environment, but this program never calls that capability.

Status: the full-input evaluator panic has **not yet been reproduced by this
standalone program**. A separate, 2480-byte Symbolica import regression is
reproduced below. The original supplied-context path encounters that import
problem on the 19 MB physical control; it is not validated for physical replay.
Successful subsets are controls, not reproducers. Preserve the actual failing
input and logs before calling a case an evaluator-construction MRE.

Investigation was deferred by the user on 2026-09-13 so the UV feature could
proceed to merge. All 3,050 physical scalar components were recovered. A
5.36 GB abstract input built successfully in 232.995 s (43.697 GB VmHWM).
A 13.887 GB combined dump was saved, but its subsequent import was interrupted
by the watchdog's two-second `ps` timeout before evaluator construction.
That infrastructure failure is not the original Symbolica panic. The large
dumps and detailed run receipts remain local, outside version control.

## Build and small control

Symbolica and its Numerica/Graphica dependencies are pinned to
`4d0a833eb8e059d1f95bdae5abed2559830b235f`. The local Cargo.lock is seeded from
the original workspace lock and subsequently resolved for this independent
package. The `dev-optim` profile keeps debug assertions enabled. In the shared
workspace the build uses `/tmp/gammaloop-aa-aa-uv-build.lock` and this package's
own target directory; do not overwrite the workspace target or dependency
checkout.

From this directory, after configuring `SYMBOLICA_LICENSE` and the compiler:

```sh
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export RUST_MIN_STACK=67108864 RUST_BACKTRACE=full
export CARGO_TARGET_DIR="$PWD/target"
cargo check --locked --profile dev-optim
cargo build --locked --profile dev-optim
target/dev-optim/symbolica-evaluator-mre --control
```

The tiny control exercises a tagged function, a constant alias, ordered
parameters and the same tuple decoder used for captured data. It evaluates
`f(7,x)+c`, with `f(7,t)=(t+1)^2+(t+1)^3` and `c=3`, at `x=2`; the result must
be 39. It is a correctness check, not a performance measurement.

## Supplied physical expression

```sh
target/dev-optim/symbolica-evaluator-mre \
  /path/to/expression.symbolica \
  inputs/params.symbolica inputs/function_map_entries.symbolica
```

The two small files in `inputs` are byte-for-byte copies from the successful
CLI18 physical-component replay. `inputs/provenance.json` identifies their
original paths/hashes and the physical scalar controls. They contain 220
ordered parameters and 63 definitions, before the two imaginary-unit aliases.
The program accepts other supplied files using the same generic format.

Each file uses the public `Atom::export` / `Atom::import` format. The parameter
container is `evaluator_probe::params(p0,...)`. Definitions are an ordered
`evaluator_probe::definitions(...)` list of
`evaluator_probe::definition(lhs,rhs,evaluator_probe::tags(...),evaluator_probe::args(...))`.
Variable keys become aliases; function keys use the supplied tags and arguments
in their original order. Insertion order matters because FunctionMap retains
the first definition of a duplicate key. The program appends
`vakint::𝑖 → Atom::i()` and `symbolica::𝑖 → Atom::i()` exactly as the measured
evaluator owner does. It adds no loader-specific aliases or replacements.

Context imports occur before creating application symbol names. Atom imports
merge the exported symbol state and remap identifiers; they do not preserve
custom application callbacks or guarantee the original registration order.
The ordinary evaluator builder does not request derivatives. Matching a
successful captured physical control is still required to check this boundary
empirically. Equal expression byte counts alone are not equality certificates.

Portable context SHA256 identities:

- `params.symbolica`: `75b09ec28d4c9fb1f4d92ce487705a4a90782216ed4328c0642e4977926f91b6`
- `function_map_entries.symbolica`: `7a54415d0414f5b949d4ac6e43d6fcddb909faed70d26f7d1537e25fc7304380`

The original failure is preserved in
[`evidence/original_cli21_failure.typ`](evidence/original_cli21_failure.typ), with
its machine-readable metadata and panic excerpt. This is historical application
evidence; new independent standalone outcomes are recorded separately.

The existing 3,389,407,083-byte Atom is
`../captured_tensor_replay_v5_candidate/largest_cli18_export/scalar_4.symbolica`.
It is one component, not the full 3050-component expression. The existing
GammaLoop constituent replay succeeded with H1/CPE5 and counted 80,026
additions, 67,213 multiplications, 21 inversions and 11 function calls. Those
counts are a useful control for the standalone importer; they do not establish
that the full expression will succeed.

The smaller physical control is
`../captured_tensor_replay_v5_candidate/full_cli18_export_mre_r1/scalar_0.symbolica`
(19,479,381 Atom bytes). Its source is the newly recovered CLI18 scalar capture.
Control executions that overlap that capture are labelled correctness-only.
Large executions must wait for the existing measurement lock and a supervised
wall/RSS budget; they must not overlap final physical benchmark timing.

## Combining an actual supplied subset or complete capture

```sh
target/dev-optim/symbolica-evaluator-mre --combine combined.symbolica \
  /path/to/scalar_0.symbolica /path/to/scalar_4.symbolica
target/dev-optim/symbolica-evaluator-mre combined.symbolica \
  inputs/params.symbolica inputs/function_map_entries.symbolica
```

Alternatively, the single input after OUTPUT may be a directory; the program
selects `scalar_INDEX.symbolica` files and sorts them by numeric index. The
launch manifest must verify the expected capture completeness and exact hashes
before treating this as a full-graph case. The program contains no fixed graph
name, component count, or expression-size threshold.

Combination uses Symbolica's `Atom::add_many`, preserving factorized
subexpressions without expansion. It exports the resulting Atom separately
before any evaluator build, writing a `.partial` file and flushing/syncing it
before the final rename. Existing output or partial files are not overwritten.
The combined dump is useful independently of the
capture program. The two commands isolate combination/import failures from
evaluator-construction failures and allow the builder to run in a fresh process.

A full input must match the expression supplied to the original evaluator
builder, including any prior selector postprocessing. This program deliberately
performs no GammaLoop postprocessing. A physical partial sum that reproduces
the same panic is also useful; identify it as a subset, retaining its component
list and hashes. The original complete GL262 direct run failed during
Symbolica construction with `advance out of bounds: the len is 1721315294 but
advancing by 10311249876`, with evaluator compilation disabled. A standalone
reproduction, its exact stage, and backtrace remain required.

The settings are H1, CPE5, one core, direct translation enabled, no hot start,
abort level zero, maximum 500 Horner variables, 1,000,000 common-pair cache
entries, maximum common-pair distance 1000 and verbosity off. The abort
callback always returns false, matching uninterrupted construction. Stages
report elapsed time, RSS and lifetime VmHWM on Linux (null on other systems).
The `symbolica_build_start` sample is the separate pre-evaluator memory boundary.
`count_operations()` is printed only after a successful
build. No numerical Monte Carlo sampling occurs.

## Explicit abstract-parameter reduction

```sh
target/dev-optim/symbolica-evaluator-mre --abstract-parameters --control
target/dev-optim/symbolica-evaluator-mre --abstract-parameters combined.symbolica
```

This mode interprets every distinct variable and whole function subtree of the
**imported** expression as an independent parameter. It calls the existing
`get_all_indeterminates(false)`, sorts its AtomViews by Symbolica's exact `Ord`,
and uses an empty FunctionMap. The input expression is never rewritten or
expanded. Symbolica's direct evaluator checks exact Var/Fun parameter matches
before builtin or function-map dispatch, so the result is a well-defined
abstract scalar. The settings and `.build()` call are shared with the context
path above.

That lookup happens after the builder's own Horner pass. A function argument
that Horner changes could stop matching its original parameter key; successful
simple controls do not prove the builder supports every possible abstraction.
Do not silently rewrite the input or alter the optimization settings to hide
such a failure.

This is an intentional mathematical reduction, not a repair of the import bug
and not a claim of physical equivalence. Import may already have changed symbol
identities. A matching evaluator panic on this valid abstract expression would
still be independently useful to the Symbolica author; retain the exact dump,
mode, counts, settings, stage and backtrace. The tiny abstraction control uses
`sin(x)+f(x)^2/cos(y)`, treats its three function subtrees as independent inputs,
and verifies the result is 4 when all three parameters equal 2. It checks that
builtins are also consumed through the exact parameter lookup.

## Separate small import regression

```sh
cargo build --locked --profile dev-optim --example import_remap
target/dev-optim/examples/import_remap write import_test.symbolica
target/dev-optim/examples/import_remap read import_test.symbolica
```

The writer registers `x,y,f` and exports `f(x,y)`; the separate reader registers
`y,x,f`. It imports `f(y,x)` instead of `f(x,y)` and fails the equality assertion.
`f` has no symmetry attributes: its ordered arguments must retain their identities.
The function keeps its original symbol ID, while its argument IDs need remapping.
Pinned `AtomView::rename_no_norm` only visits a function's arguments when that
function's own ID appears in the remapping table. No application dependencies
or physics are involved. Original evidence is in `import_remap_r1/receipt.json`,
`write.log`, `read.log`, `source.rs` and the 2480-byte `tiny.symbolica` input
(SHA256 `630e1b3c5a67d8536bb3e604f39af259f6e267429bf6a823e175b79f6acd19d3`).

The separate context inspection command is
`--inspect-context PARAMS DEFINITIONS`. `context_inspection.log` shows that the
captured first parameters import as unrelated symbols, before the scalar is
read. `physical_scalar0_control_r1` preserves the resulting `UndefinedVariable`
error for `gammalooprs::mUV`; this is not the original large `.build()` panic.

## Validation on 2026-09-13

`cargo check`, `cargo build` and `cargo clippy --all-targets --no-deps -- -D warnings`
passed under the build lock, using the package target and debug assertions.
`build_receipt.json` binds the source, dependency lock and executable identities.
The executable SHA256 is
`2a1a680f267fa56578d6a99aa983ef4ab4567661ce657c07181faa86eddaca79`.

Three checks passed with this executable: `tiny_control_r2` (39),
`tiny_abstract_control_r1` (4), and `physical_scalar0_abstract_r1` (evaluator
construction succeeds). These overlapped ongoing scalar capture and provide
correctness evidence only. The abstract scalar0 case imports 20,073,529 Atom
bytes, derives 34 parameters, and counts 3132 additions, 3383 multiplications,
one inversion and five function calls. Build time was 0.539004993 seconds.
Immediately before `.build()`, RSS was 80,584,704 bytes and lifetime VmHWM was
130,838,528 bytes; after the build lifetime VmHWM was 244,830,208 bytes. Full
logs and guarded command receipts are retained in the respective directories.

The 3.389 GB physical context control and original large evaluator-construction
panic remain unvalidated in this standalone package. The 5.36 GB abstract
combined-input control passed as recorded above. The supplied-context
import issue is unresolved; no physical equivalence is claimed for the
abstract-parameter mode.
