# Generalised LTD Implementation Logs

## Context

This file records the migration of the Python `generalised_ltd` prototype into the
GammaLoop Rust workspace.

The implementation is being developed in staged steps so that GammaLoop remains
buildable and testable throughout the migration. The main regression command for
GammaLoop compatibility is:

```bash
just test_gammaloop
```

At the start of Step II, the selected GammaLoop test baseline was known to pass:

```text
980 tests run: 980 passed, 128 skipped
```

## Completed Before This Log

### Step I

Step I generalized the existing GammaLoop CFF expression structure in place.

Main changes:

- Introduced a generalized orientation expression model with orientation data,
  variants, numerator energy maps, rational variant prefactors, generalized
  half-edge prefactors, and generalized surfaces.
- Added support in the expression structure for the sampling scale `M` where it
  can appear structurally.
- Preserved compatibility with the existing CFF generation path.
- Verified the change with:
  - `cargo fmt`
  - `cargo check`
  - `just clippy`
  - `just test_gammaloop`

The Step I work was committed locally as:

```text
5d78c02f generalize cff expression for step i
```

## Step II Agreed Plan

Create a new workspace crate:

```text
crates/generalise-ltd
```

The Cargo package name is `generalise-ltd`; the Rust library crate name is
`generalise_ltd`.

### Public Representation Modes

The public generation modes are:

- `ltd`
- `cff`

The former Python/prototype name `hybrid` is removed. The derivative-free
repeated-propagator representation is called `ltd`.

### Internal and Diagnostic Modes

The crate may also expose the following modes for diagnostics, testing, and
standalone CLI exploration:

- `pure_ltd`: available in the `bin` CLI and tests. Repeated propagators are
  allowed, but the generation emits a warning that numerator derivatives would
  be required and that the result is diagnostic only.
- `old_cff`: available behind an `old_cff` feature. This is the migrated legacy
  GammaLoop CFF generator, used for parity testing. It must cleanly error when
  higher-power energy numerator dependence is requested.

`old_cff` is not exposed to GammaLoop integration code.

### Expression Ownership

Move the generalized `CFFExpression` model and related expression/surface/tree
objects out of `gammalooprs` and into `generalise_ltd`.

The neutral canonical type is:

```rust
ThreeDExpression<O>
```

with:

```rust
type CFFExpression<O> = ThreeDExpression<O>;
```

kept as a compatibility alias.

### Graph Handling

All graph manipulation must use `linnet`. No custom graph topology
implementation should be introduced.

For repeated propagator grouping, the key is:

- momentum signature canonicalized up to an overall sign,
- resolved mass expression identity, or exact numerical value when the resolved
  mass expression is numerical only.

The overall sign is retained and propagated into the energy-map construction.

Numerical mass identity is exact, not tolerance based.

### Symbolic Handling

All symbolic manipulation must use `symbolica`.

If an efficient targeted operation is missing from Symbolica or linnet, implement
the smallest local helper in `crates/generalise-ltd/src/utils.rs`, using
Symbolica/linnet primitives, and record the missing upstream functionality here.

### Feature-Gated Components

The GammaLoop library path should compile only the expression and generation
logic it needs.

Optional components:

- `bin`: standalone CLI.
- `display`: colored/tabled pretty printing.
- `eval`: standalone evaluator support.
- `symbolica-eval`: Symbolica-based evaluator construction.
- `compiled-eval`: compiled complex-valued evaluator support.
- `old_cff`: migrated legacy CFF generator for parity tests and CLI exploration.
- `test-support`: helpers for `pure_ltd` and multi-way comparison tests.

### Validation Commands

The implementation should be validated incrementally with:

```bash
cargo fmt
cargo check
cargo test -p generalise-ltd --no-default-features
cargo test -p generalise-ltd --all-features
just clippy
just test_gammaloop
```

## Live Implementation Log

### 2026-04-29

- Confirmed no further clarification is required before starting implementation.
- Confirmed the working tree was clean before beginning Step II changes.
- Started by creating this implementation log so the plan and implementation
  history remain visible in the top-level repository.
- Confirmed all non-ignored Python `generalised_ltd` tests are intended to be
  ported to Rust and included in the `just test_gammaloop` gate.
- Confirmed `just test_gammaloop` currently enumerates packages explicitly, so
  `generalise-ltd` must be added to that package list when its tests are added.
- Added the initial `crates/generalise-ltd` crate skeleton with feature-gated
  CLI support.
- Added core expression, surface, tree, symbols, generation option, and utility
  modules. The expression model is independent of `gammalooprs` by making the
  E-surface/H-surface payload types generic while keeping shared IDs and
  expression structure in `generalise_ltd`.
- Rewired `gammalooprs` to depend on `generalise-ltd`.
- Added `generalise-ltd` to the explicit `just test_gammaloop` package list.
- Replaced GammaLoop's local CFF expression/surface/tree definitions with
  re-exports and aliases to the new crate. GammaLoop keeps the concrete
  `Esurface` and `Hsurface` payload implementations locally.
- Added GammaLoop implementations of the new orientation and raised-esurface
  view traits so existing filtering and residue logic can continue to work with
  the moved expression structure.
- Verified `cargo check -p generalise-ltd -p gammalooprs` after the structural
  move; the first pass succeeded after removing migration warnings.
- Verified `RUSTFLAGS=-Dwarnings cargo check -p generalise-ltd -p gammalooprs`.
- Verified `RUSTFLAGS=-Dwarnings cargo check -p gammalooprs --tests`.
- Verified `cargo test -p generalise-ltd --no-default-features`.
- Verified `cargo test -p generalise-ltd --all-features`.
- Added a feature-gated `old_cff` mode selector and high-energy-power rejection
  path in the new crate. The legacy generator body still needs to be migrated
  from GammaLoop into the `old_cff` feature module.
- Verified `just clippy`.
- The first `just test_gammaloop` attempt after the structural move failed
  during compilation because the local `target/` directory had exhausted the
  available disk space. `target/` was about 38 GiB, so `cargo clean` was run and
  removed 47.4 GiB of build artifacts.
- After rebuilding from a clean target directory, `just test_gammaloop` reached
  the test phase but failed 60 tests. The common root cause is that expression
  atoms emitted from the moved `generalise_ltd` code now use namespaced
  Symbolica symbols such as `generalise_ltd::theta(generalise_ltd::sigma(6))`,
  while the GammaLoop evaluators and parameter builders still expect the
  existing unqualified GammaLoop symbols such as `theta(sigma(6))`.
- Next action: restore GammaLoop-compatible Symbolica symbol identities for the
  moved expression code without reintroducing a `gammalooprs` dependency into
  `generalise_ltd`.
- Added GammaLoop-side extension traits that render moved orientation
  expressions, CFF variants, and surface substitutions through GammaLoop's
  existing `GS` symbol set. This keeps `generalise_ltd` standalone symbols
  intact while preserving GammaLoop evaluator compatibility during the
  migration.
- Updated GammaLoop CFF, evaluator, settings, and UV profile/localization call
  sites to use the GammaLoop-symbol rendering methods.
- Verified `RUSTFLAGS=-Dwarnings cargo check -p gammalooprs --tests` after the
  symbol-compatibility fix.
- Re-verified `cargo test -p generalise-ltd --no-default-features`.
- Re-verified `cargo test -p generalise-ltd --all-features`.
- Re-verified `just clippy` after the GammaLoop-symbol compatibility update.
- Re-verified the full curated gate with `just test_gammaloop` after adding
  `generalise-ltd` to its package list:

```text
974 tests run: 974 passed, 128 skipped
```
- Resumed Step II implementation after the compatibility checkpoint. No
  clarification blocker exists; the earlier stop was only a checkpoint after
  restoring `just test_gammaloop`.
- Copied the Python prototype's example DOT graphs and example scripts into
  `crates/generalise-ltd/examples/`.
- Started porting the Python test surface into Rust with the DOT graph layer:
  added `graph_signatures`, `graph_io`, and `validator` modules using
  `linnet::parser::DotGraph` for DOT parsing and graph topology access.
- Added Rust tests corresponding to the Python validation and graph-shape tests:
  all valid examples except `noisy_example.dot`, noisy rejection, one-loop
  10/15-external shapes, box graph shape, and repeated-group detection.
- Added Rust coverage for the copied example script executability checks,
  five-loop no-repeat graph metadata, and the diagnostic warning row helper for
  repeated-propagator `ltd` profiling rows.
- Verified the current graph/test-surface port with:

```text
cargo test -p generalise-ltd --no-default-features --lib -- --test-threads=1 --nocapture
12 passed
```
- Ported the `graph_from_signatures` layer using Symbolica for expression
  parsing. The Rust code now extracts `prop(q,m)` calls from the Symbolica atom
  tree, reconstructs integer momentum signatures, and emits DOT in both the
  default and `vakint` modes.
- Added corresponding Rust tests for the two Python `graph_from_signatures`
  cases and verified:

```text
cargo test -p generalise-ltd --no-default-features --lib -- --test-threads=1 --nocapture
15 passed
cargo test -p generalise-ltd --features bin --lib -- --test-threads=1 --nocapture
15 passed
cargo test -p generalise-ltd --features bin --bins -- --test-threads=1 --nocapture
0 passed
```
- Smoke-tested the feature-gated CLI command:

```text
cargo run -p generalise-ltd --features bin --bin generalise-ltd -- graph-from-signatures --signatures 'prop(k1+p1,mA)*prop(k1+p1-q1,mB)*prop(k1-p2+q2,mC)*prop(k1,mD)'
```
- Ported the first cut-structure/generation layer:
  - exact-sign LTD residue enumeration from loop-line signatures,
  - non-repeated `ltd` generation from parsed GammaLoop DOT graphs into
    `ThreeDExpression<OrientationID>`,
  - linear surface interning,
  - loop/edge energy-map solving for unimodular bases.
- Repeated-propagator `ltd` generation currently reaches an explicit
  `LtdRepeatedPropagators` error until the confluent repeated-propagator kernel
  is ported.
- Ported the pure-CFF graph recursion:
  - acyclic edge-orientation enumeration,
  - source/sink boundary surface construction,
  - exhaustive CFF vertex-contraction branches,
  - repeated-topology support at the old pure-CFF level.
- CFF contraction branches are stored as multiple `CFFVariant` entries under
  the same orientation/numerator energy map.
- Wired the `old_cff` feature mode to the migrated pure-CFF generator for parsed
  DOT inputs and kept the explicit higher-energy-power rejection path.
- Replaced the placeholder `build` CLI with a feature-gated implementation that:
  - accepts GammaLoop DOT input,
  - accepts the prototype-style `--family`, `--energy-degree-bounds`,
    `--numerator-expr`, `--uniform-numerator-sampling-scale`, `--pretty`, and
    `--show-details-for-orientation` flags,
  - emits simplified Rust JSON containing graph metadata, validation, and the
    serialized `ThreeDExpression`.
- Ported the energy-degree UV convergence diagnostic, including the
  non-coordinate direction scan from the Python prototype. Generation now
  rejects non-convergent energy-degree bounds before building a representation.
- Ported the automatic numerator-expression helper for sparse energy-degree
  bounds and wired `--numerator-expr auto` into the feature-gated `build` CLI.
- Current major remaining implementation block: the confluent
  repeated-propagator `ltd` kernel and the bounded-degree generalized `cff`
  kernel. The public `ltd` path still explicitly errors on repeated
  propagators until that kernel is ported; the old pure-CFF recursion is
  available and feature-gated `old_cff` parity support is wired for parsed DOT
  inputs.
- Verified the expanded no-default-feature test set:

```text
cargo test -p generalise-ltd --no-default-features --lib -- --test-threads=1 --nocapture
23 passed
cargo test -p generalise-ltd --features bin --lib --bins -- --test-threads=1 --nocapture
21 passed, 0 bin tests
cargo test -p generalise-ltd --all-features --lib --bins -- --test-threads=1 --nocapture
26 passed, 0 bin tests
RUSTFLAGS=-Dwarnings cargo check -p generalise-ltd --all-features
passed
RUSTFLAGS=-Dwarnings cargo check -p gammalooprs --tests
passed
just clippy
passed
just test_gammaloop
996 tests run: 996 passed, 128 skipped
```

- Smoke-tested the showcase-style CLI command:

```text
cargo run -q -p generalise-ltd --features bin --bin generalise-ltd -- build --dot crates/generalise-ltd/examples/graphs/box_pow3.dot --family cff --energy-degree-bounds 0:1,1:1,2:1,3:4 --numerator-expr auto --uniform-numerator-sampling-scale none --pretty --show-details-for-orientation '0++xxx|N2'
Cff 6 62
cargo run -q -p generalise-ltd --features bin --bin generalise-ltd -- build --dot crates/generalise-ltd/examples/graphs/box_pow3.dot --family cff --energy-degree-bounds 0:1,1:1,2:0,3:4 --numerator-expr auto
dot(edges[0], ext[0]) * dot(edges[1], ext[1]) * dot(edges[3], ext[0]) * dot(edges[3], ext[1]) * dot(edges[3], ext[2]) * dot(edges[3], ext[0])
```

## Missing Upstream Functionality To Report

This section will be kept current as implementation proceeds.

### Symbolica

- A small exact-rational linear algebra layer remains local in
  `crates/generalise-ltd/src/utils.rs` for rank checks and rational solves over
  the integer loop-energy signature matrices. A compact Symbolica/numerica API
  for exact matrix rank and exact linear solve over rationals would remove this
  local helper.
- Polynomial/interpolation combinatorics for bounded numerator sampling are
  currently local Rust helpers over exact rationals. If Symbolica grows a
  targeted API for small multivariate interpolation/contact-term coefficient
  extraction, the local implementation can be simplified.

### Linnet

- Contact-sector and lower-sector reconstruction still need local adapters over
  `ParsedGraph` metadata. A linnet helper to contract selected DOT/parsed graph
  edges while preserving deterministic edge-id maps and momentum/mass metadata
  would make this path cleaner.

### Local Interim Limitations

- The Rust loop-energy solver currently accepts the unimodular integer-basis
  cases needed by the shipped graph tests. If a future topology requires
  rational coefficients in energy maps, `LinearEnergyExpr` must be generalized
  from integer coefficients to exact rational coefficients.

## Resume: Repeated-Propagator LTD Kernel

- No additional clarification blocker was found. The next implementation block
  is the repeated-propagator `ltd` kernel, corresponding to the Python
  `_build_confluent_hybrid_bundle` routine, with the public Rust-facing name
  kept as `ltd`.
- Planned edits for this block:
  - add the numerator sampling-scale mode to `Generate3DExpressionOptions`,
  - port the logical-channel collapse and confluent denominator derivative
    expansion,
  - port the affine physical numerator-map reconstruction and the
    bounded-degree finite-difference numerator samples, including the optional
    uniform scale `M`,
  - group denominator variants by common numerator energy map in
    `ThreeDExpression`,
  - organize stateful Rust generation logic behind builder structs and methods
    rather than mirroring the Python implementation as free-floating
    functions,
  - replace the current repeated-LTD rejection test with structural tests for
    repeated `ltd`, multiloop repeated `ltd`, and M-dependent finite-difference
    maps.

### Progress 2026-04-29 10:02 CEST

- Implemented the repeated-propagator public `ltd` path with a Rust builder
  organization:
  - `RepeatedLtdBuilder` owns the repeated-channel construction,
  - `EnergySolver` owns loop/edge energy substitutions and derivative
    matrices,
  - `LinearSurfaceInterner` owns surface-cache insertion,
  - `DenominatorDerivativeExpander` owns confluent denominator derivatives,
  - `AffineWeightSystem` owns the physical numerator-map interpolation.
- Added `NumeratorSamplingScaleMode` to `Generate3DExpressionOptions` and wired
  the CLI `--uniform-numerator-sampling-scale` values `none`,
  `beyond-quadratic`, and `all` into generation.
- Ported the repeated-channel confluent LTD construction:
  - logical repeated channels,
  - cut and dual-surface derivative factors,
  - affine physical numerator sample orbit for the no-bounds case,
  - bounded-degree finite-difference numerator samples,
  - optional uniform sampling scale `M` in energy maps and inverse powers in
    variant prefactors,
  - grouped variants by common numerator map and common prefactor/half-edge
    structure.
- Replaced the repeated-propagator rejection test with structural repeated-LTD
  tests and added checks for:
  - one-loop repeated `box_pow3.dot`,
  - multiloop repeated `sunrise_pow4.dot`,
  - bounded repeated LTD maps using the uniform scale `M`.
- Local CLI shape check:

```text
cargo run -q -p generalise-ltd --features bin --bin generalise-ltd -- build \
  --dot crates/generalise-ltd/examples/graphs/box_pow3.dot \
  --family ltd --pretty --json-out /tmp/rs_box_pow3_ltd.json
17 orientations, 73 grouped variants, first label +xxxxx|N1

cargo run -q -p generalise-ltd --features bin --bin generalise-ltd -- build \
  --dot crates/generalise-ltd/examples/graphs/box_pow3.dot \
  --family ltd --energy-degree-bounds 3:4 \
  --uniform-numerator-sampling-scale beyond-quadratic \
  --json-out /tmp/rs_box_pow3_ltd_M.json
M appears in edge energy maps and variant uniform_scale_power entries
```

- Current focused test result:

```text
cargo test -p generalise-ltd --all-features --lib --bins -- --test-threads=1 --nocapture
28 passed, 0 failed
```

- Follow-up gates after formatting and Clippy fixes:

```text
RUSTFLAGS=-Dwarnings cargo check -p generalise-ltd --all-features
passed
cargo test -p generalise-ltd --all-features --lib --bins -- --test-threads=1 --nocapture
28 passed, 0 failed
cargo test -p generalise-ltd --no-default-features --lib -- --test-threads=1 --nocapture
26 passed, 0 failed
RUSTFLAGS=-Dwarnings cargo check -p gammalooprs --tests
passed
just clippy
passed
just test_gammaloop
998 tests run: 998 passed, 128 skipped
Nextest run ID 2484ea8f-6443-4067-9456-a9235dee55cd
```

- Clarification received after this checkpoint: the migrated legacy
  pure-CFF/`old_cff` implementation should remain as close as possible to the
  original GammaLoop implementation, because it is the comparison anchor for
  the affine/linear numerator regime. I will not refactor that anchor into
  linnet or otherwise structurally rewrite it. New generalized `ltd`/`cff`
  functionality remains subject to the modern Symbolica/linnet requirement.

### Progress 2026-04-29 10:13 CEST

- Changed `cff` generation with energy-degree bounds above the legacy affine
  regime to fail cleanly with
  `CffHigherEnergyPowerNotImplemented` instead of silently returning the old
  pure-CFF structure. This keeps the current Rust port honest until the
  bounded-degree generalized CFF kernel is implemented.
- The earlier CFF CLI smoke log with `--energy-degree-bounds ...3:4` is now
  obsolete by design: that command should fail until the new bounded-degree
  CFF kernel lands.

### Progress 2026-04-29 10:24 CEST

- Reconfirmed the latest implementation constraint: the migrated legacy
  pure-CFF/`old_cff` anchor must not be modernized or structurally refactored.
  It remains the affine/linear comparison anchor. New generalized `ltd`/`cff`
  paths remain the place for the modern Symbolica/linnet-backed implementation.
- Removed an incomplete in-progress bounded-CFF hook that referenced an
  undefined builder. The active higher-power `cff` path is back to the clean
  `CffHigherEnergyPowerNotImplemented` error until bounded CFF lands as a
  complete implementation slice.
- Verified the repository is back to a buildable state for the affected crate:

```text
cargo fmt
cargo check -p generalise-ltd --all-features
passed
```
- Cleaned the copied example scripts so they exercise the current Rust
  feature-gated CLI (`cargo run -p generalise-ltd --features bin --bin
  generalise-ltd -- ...`) instead of stale Python prototype commands. The
  scripts no longer use the retired `hybrid` mode; repeated-propagator
  derivative-free examples use `ltd`, and the diagnostic reference uses
  `pure-ltd` where applicable.
- Updated the repeated-propagator profiling warning helper to warn on
  `pure_ltd` diagnostic rows rather than public `ltd` rows, matching the new
  naming policy.
- Verified the script/diagnostic cleanup with:

```text
cargo fmt
cargo test -p generalise-ltd --all-features --lib --bins -- --test-threads=1
29 passed, 0 failed
```

### In Progress 2026-04-29 10:35 CEST

- Starting the feature-gated human display block for the Rust CLI:
  - add a `display` module behind the existing `display` feature,
  - render summary, surface, orientation/variant, and optional orientation
    detail tables from `ThreeDExpression`,
  - wire `generalise-ltd build --pretty` to print the display while preserving
    JSON output through `--json-out`,
  - keep the renderer independent from `gammalooprs`.

### Progress 2026-04-29 10:48 CEST

- Added the feature-gated `display` module and re-exported
  `DisplayOptions`/`render_expression_summary` behind the `display` feature.
- The renderer now emits:
  - a summary table with graph counts, representation, repeated groups,
    energy-degree bounds, numerator metadata, and sampling-scale mode,
  - a linear-surface table with kind/origin/numerator-only flags,
  - an orientation/variant table with labels, origins, rational prefactors,
    half-edge prefactors, inverse `M` powers, numerator surfaces, and compact
    denominator-tree shapes,
  - optional detailed loop/edge energy-map and denominator-tree tables for a
    selected orientation.
- Wired `generalise-ltd build --pretty` to print this display. When
  `--json-out` is supplied, the JSON is still written to that path and the
  display is printed to stdout.
- Added `--no-color` for deterministic plain terminal output.
- Fixed feature wiring so `cargo check -p generalise-ltd --features cli`
  compiles: `cli` now enables `serde_json` and `display`, while `pure-ltd`
  remains available only when `bin` or `test-support` is enabled.
- Verified this block with:

```text
cargo fmt
cargo check -p generalise-ltd --features cli
cargo check -p generalise-ltd --no-default-features
cargo run -q -p generalise-ltd --features bin --bin generalise-ltd -- build \
  --representation ltd --dot crates/generalise-ltd/examples/graphs/box_pow3.dot \
  --pretty --no-color --show-details-for-orientation '+xxxxx|N1' \
  --json-out /tmp/rs_box_pow3_display.json
cargo test -p generalise-ltd --all-features --lib --bins -- --test-threads=1
30 passed, 0 failed
```

### In Progress 2026-04-29 11:05 CEST

- Continuing Step II after the display checkpoint.
- Current focus is the bounded-degree generalized `cff` kernel. The legacy
  migrated `old_cff`/pure-CFF anchor remains untouched except for integration
  boundaries needed by the new crate.
- Before coding the bounded-CFF path, I am re-reading the Python
  `orientation_bundle.py` implementation and the relevant `linnet`/Symbolica
  APIs to avoid reimplementing graph or symbolic operations that already exist.
- Expected next implementation slice:
  - introduce a dedicated bounded-CFF builder owned by a struct/method API,
  - port only the necessary symbolic helpers using Symbolica primitives in
    `utils.rs`,
  - use linnet-parsed graph topology and document any missing graph ergonomics
    that force a local adapter,
  - keep the public `cff` higher-power path failing cleanly until the Rust
    bounded-CFF builder has enough parity tests to replace that error.

### Progress 2026-04-29 11:28 CEST

- Added the first bounded-degree generalized `cff` implementation slice:
  one-loop, non-repeated, quadratic-or-lower energy caps using the E-surface
  contact-completion construction.
- The implementation is owned by `BoundedCffBuilder`, keeping the new
  generalized code method-oriented rather than a flat transduction of Python
  helper functions.
- Public `cff` now dispatches higher-than-affine bounds to this builder. Cases
  outside this slice, including repeated channels and cubic-or-higher caps, keep
  failing cleanly with `CffHigherEnergyPowerNotImplemented`.
- Added exact rational interpolation/contact polynomial helpers for this slice.
  These are local rational-polynomial utilities, not a replacement symbolic
  engine; Symbolica remains the only symbolic-expression layer.
- Added `LinearEnergyExpr` helpers for zero/one checks and internal-edge
  remapping.
- Added a bounded-CFF Rust unit test for `box.dot` with `--energy-degree-bounds
  0:2`, checking that the structure grows beyond ordinary CFF, remains
  E-surface-only, and marks bounded-CFF contact-completion variants.
- CLI smoke check for the new slice:

```text
cargo run -q -p generalise-ltd --features bin --bin generalise-ltd -- build \
  --representation cff --dot crates/generalise-ltd/examples/graphs/box.dot \
  --energy-degree-bounds 0:2 --json-out /tmp/rs_box_cff_b2.json

shape: 32 orientation maps, 38 variants, 14 linear surfaces
all generated linear surfaces are E-surfaces
```

- Focused verification:

```text
cargo fmt
cargo test -p generalise-ltd --all-features --lib --bins -- --test-threads=1
31 passed, 0 failed
RUSTFLAGS=-Dwarnings cargo check -p generalise-ltd --all-features
passed
cargo check -p generalise-ltd --no-default-features
passed
```

- Missing linnet ergonomic noted for later upstreaming: this slice currently
  uses a small metadata-level contraction adapter on `ParsedGraph` for the
  contact sectors. The adapter is built from the linnet-parsed graph metadata,
  but a direct "contract selected parsed/DOT edges and preserve edge-id mapping"
  helper in linnet would remove this local code and make the intent clearer.

### In Progress 2026-04-29 11:38 CEST

- Inventoried the Python `tests/test_cli.py` suite: 73 non-ignored tests.
- Rust coverage currently includes the validation/graph-shape/script
  executability tests, graph-from-signatures tests, energy-UV diagnostics,
  display smoke coverage, pure `ltd`/`cff` construction tests, repeated `ltd`
  confluent construction, uniform-M repeated LTD maps, and the first
  bounded-CFF one-loop quadratic structural test.
- Remaining large blocks are:
  - Symbolica/evaluator compile/evaluate tests,
  - numerical CFF/LTD parity tests,
  - recursive bounded-CFF tests for repeated channels, multiloop lower sectors,
    and cubic-or-higher caps,
  - CLI `test`/profile/stability commands,
  - runtime uniform-scale validation and complex-only evaluator policy.

### Progress 2026-04-29 12:02 CEST

- Started porting the recursive quadratic bounded-CFF algorithm needed for
  one-loop repeated channels:
  - recursive remainder/contact split,
  - deletion sectors,
  - helper numerator surfaces,
  - one-loop lower-sector sign experiment.
- A temporary numerical adapter from the Rust JSON expression to the Python
  prototype evaluator confirmed:
  - the non-repeated `box.dot`, bound `0:2`, agrees with the Python LTD value at
    high precision (`~5e-79` absolute difference for `edges[0][0]**2`),
  - the repeated `box_pow3.dot`, bounds `0:2,1:1,2:1,3:2`, does **not** yet
    agree with the Python bounded CFF/hybrid value.
- Because of that mismatch, recursive repeated-channel bounded CFF dispatch is
  intentionally disabled again. The scaffold remains in the source as
  in-progress code, but public `cff` still returns the clean unsupported error
  for repeated/multiloop/cubic bounded-CFF cases until the lower-sector
  component-product construction is ported correctly.
- Re-verified the crate after disabling the unsafe dispatch:

```text
cargo fmt
RUSTFLAGS=-Dwarnings cargo check -p generalise-ltd --all-features
cargo test -p generalise-ltd --all-features --lib --bins -- --test-threads=1
31 passed, 0 failed
RUSTFLAGS=-Dwarnings cargo check -p gammalooprs --tests
passed
cargo test -p generalise-ltd --no-default-features --lib -- --test-threads=1
28 passed, 0 failed
just test_gammaloop
1000 tests run: 1000 passed, 128 skipped
Nextest run ID 6064bd7e-7dc0-4262-902c-d99e1c45a36e
```

### In Progress 2026-04-29 12:18 CEST

- Continuing with the bounded-CFF recursive lower-sector component-product
  port. This is the current blocker for repeated-channel and multiloop
  bounded-CFF numerical parity tests.
- Immediate goal: replace the disabled one-loop recursive bounded-CFF dispatch
  scaffold with a correct lower-sector base construction, starting from the
  Python `build_lower_sector_cff_bundle` logic.
- The legacy `old_cff` anchor remains untouched.

### Progress 2026-04-29 12:43 CEST

- Ported the lower-sector component-product base used by recursive quadratic
  bounded CFF into a dedicated `LowerSectorCffBuilder`.
- Re-enabled one-loop quadratic recursive bounded-CFF dispatch, including
  repeated channels with quadratic-or-lower caps.
- Added helpers for:
  - vector-matroid component decomposition of loop-energy rows,
  - component projection onto a local loop basis,
  - lower-sector component product assembly,
  - rank-deficient particular loop-energy solves with free loop energies set to
    zero.
- Temporary high-precision numerical adapter checks against the Python
  prototype:

```text
box.dot, cff, bounds 0:2, numerator edges[0][0]**2:
  |rust bounded CFF - python LTD| ~= 5e-79

box_pow3.dot, cff, bounds 0:2,1:1,2:1,3:2,
numerator edges[0][0]**2 * edges[1][0] * edges[2][0] * edges[3][0]**2:
  |rust bounded CFF - python bounded CFF| ~= 8e-84
  |rust bounded CFF - python LTD/hybrid| ~= 2e-76
```

- Added a Rust structural regression for one-loop repeated quadratic CFF
  completion.
- Focused verification:

```text
cargo fmt
RUSTFLAGS=-Dwarnings cargo check -p generalise-ltd --all-features
cargo test -p generalise-ltd --all-features --lib --bins -- --test-threads=1
32 passed, 0 failed
```

- Probed broadening the same recursive path to multiloop quadratic caps. It
  builds, but the `sunrise_pow4.dot` probe produced H-surfaces, while the
  Python bounded-CFF multiloop quadratic path is expected to remain
  E-surface-only. Dispatch is therefore intentionally kept to the validated
  one-loop quadratic class until the auxiliary graph reconstruction path for
  multiloop lower sectors is ported.
- Probed the direct finite-pole E-surface lift for one-loop non-repeated
  higher-power caps. A single isolated high-power edge matches the Python
  bounded-CFF/LTD value at high precision (`box.dot`, bounds `0:3`, numerator
  `edges[0][0]**3`). Mixed high-power plus spectator caps do not match with the
  direct lift alone, so dispatch is restricted to:
  - one-loop non-repeated all caps <= 2,
  - one-loop non-repeated exactly one nonzero cap and that cap > 2,
  - one-loop repeated quadratic-or-lower caps through the recursive
    lower-sector component-product path.
- Added a Rust structural regression for the single isolated high-power
  one-loop CFF completion.
- Re-verified after the dispatch restriction:

```text
cargo fmt
RUSTFLAGS=-Dwarnings cargo check -p generalise-ltd --all-features
cargo test -p generalise-ltd --no-default-features --lib -- --test-threads=1
30 passed, 0 failed
cargo test -p generalise-ltd --all-features --lib --bins -- --test-threads=1
33 passed, 0 failed
RUSTFLAGS=-Dwarnings cargo check -p gammalooprs --tests
passed
just test_gammaloop
1002 tests run: 1002 passed, 128 skipped
Nextest run ID 6e040115-fa43-4502-925e-fca110624966
```

### In Progress 2026-04-29 13:08 CEST

- Continuing with bounded-CFF mixed high-power/affine-spectator support.
- Immediate target is the Python `_known_recursive_terms` /
  `build_known_factor_bounded_cff_bundle` logic, but organized as Rust methods
  on dedicated structs.
- Dispatch will stay conservative: new cases are enabled only after structural
  and temporary numerical checks against the Python prototype pass.

### In Progress 2026-04-29 12:34 CEST

- Resumed implementation of the known-factor bounded-CFF path.
- The next edit will add a Rust `KnownLinearExpr` and a dedicated
  `KnownFactorCffBuilder` so the port keeps recursive state behind methods
  instead of accumulating free-floating helper functions.
- Initial dispatch will remain limited to one-loop non-repeated mixed
  high-power cases until numerical checks show the same parity as the Python
  prototype.

### Progress 2026-04-29 12:44 CEST

- Added `KnownLinearExpr` and a method-oriented `KnownFactorCffBuilder`.
- Enabled one-loop non-repeated mixed high-power bounded-CFF dispatch through
  the known-factor recursion.
- Added Rust structural regressions for:
  - normal box bounds `0:1,1:1,2:1,3:3`,
  - normal box bounds `0:4,1:2`.
- Temporary adapter checks against the Python prototype evaluator:

```text
box.dot, cff, bounds 0:1,1:1,2:1,3:3,
numerator edges[0][0] * edges[1][0] * edges[2][0] * edges[3][0]**3:
  |rust bounded CFF - python bounded CFF| = 0
  |rust bounded CFF - python LTD| ~= 2.4e-80

box.dot, cff, bounds 0:4,1:2,
numerator edges[0][0]**4 * edges[1][0]**2:
  |rust bounded CFF - python bounded CFF| ~= 1.6e-81
  |rust bounded CFF - python LTD| ~= 4.7e-79
```

- Focused verification:

```text
cargo fmt
cargo test -p generalise-ltd --all-features --lib --bins -- --test-threads=1
35 passed, 0 failed
```

- Next target is repeated-channel high-power bounded CFF. This requires porting
  the Python `_channel_recursive_terms` logic while preserving the current
  `old_cff` anchor untouched.

### Progress 2026-04-29 12:58 CEST

- Ported the one-loop repeated-channel bounded-CFF recursion:
  - interpolation nodes,
  - channel normal-form decomposition,
  - channel deletion sectors,
  - relative signature signs for repeated members,
  - optional uniform-scale channel substitutions.
- Broadened the one-loop known-factor fallback so high-power non-repeated
  edges with repeated spectator channels are handled by the same method-oriented
  recursion.
- Added Rust structural regressions for:
  - `box_pow3.dot`, bounds `0:1,1:1,3:4`,
  - `box_pow3.dot`, bounds `0:3,1:1,2:1,3:1`.
- Updated the unsupported-case regression to use the still-unported multiloop
  high-power `sunrise_pow4.dot` sector.
- Temporary adapter checks against the Python prototype evaluator:

```text
box_pow3.dot, cff, bounds 0:1,1:1,3:4,
numerator edges[0][0] * edges[1][0] * edges[3][0]**4:
  |rust bounded CFF - python bounded CFF| ~= 2.2e-81

box_pow3.dot, cff, bounds 0:3,1:1,2:1,3:1,
numerator edges[0][0]**3 * edges[1][0] * edges[2][0] * edges[3][0]:
  |rust bounded CFF - python bounded CFF| ~= 2.1e-82
```

- Focused verification:

```text
cargo fmt
cargo test -p generalise-ltd --all-features --lib --bins -- --test-threads=1
37 passed, 0 failed
```

### Progress 2026-04-29 13:04 CEST

- Added structural coverage for uniform-scale `M` in repeated high-power CFF
  channel recursion.
- Temporary Python adapter checks for `box_pow3.dot`, bounds `0:1,1:1,3:4`,
  `--uniform-numerator-sampling-scale all`, and
  numerator `edges[0][0] * edges[1][0] * edges[3][0]**4`:

```text
M = 1.7:
  |rust bounded CFF - python bounded CFF| ~= 1.5e-80

M = 5.3:
  |rust bounded CFF - python bounded CFF| ~= 8.9e-79
```

- Focused verification:

```text
cargo fmt
cargo test -p generalise-ltd --all-features --lib --bins -- --test-threads=1
38 passed, 0 failed
```

### Progress 2026-04-29 13:17 CEST

- Ported the graphic auxiliary lower-sector reconstruction used when a
  vector-matroid component has more denominator edges than rank.
- The implementation now uses the crate's Symbolica/linnet-backed
  `reconstruct_dot` path instead of the synthetic-cycle fallback for those
  components, then parses the reconstructed DOT back through the standard
  GammaLoop-compatible DOT parser.
- Re-enabled bounded-CFF recursive quadratic support beyond one loop only for
  the currently validated conservative class: all caps <= 2 and at most one
  genuinely quadratic edge outside one-loop topologies.
- Added a structural regression for `sunrise_pow4.dot` with bounds `0:2` that
  checks denominator surfaces remain E-surfaces.
- Temporary adapter check against the Python prototype evaluator:

```text
sunrise_pow4.dot, cff, bounds 0:2, numerator edges[0][0]**2:
  |rust bounded CFF - python bounded CFF| ~= 6.6e-83
```

- Probed broader multiloop quadratic combinations
  (`proper_iterated_sandwiched_bubble.dot` and `four_loop_stress.dot`). They
  build with E-only denominator surfaces, but the temporary numerical adapter
  does not yet match the Python bounded-CFF values for all cases. Dispatch is
  therefore kept conservative until the remaining lower-sector/product assembly
  mismatch is found.
- Focused verification:

```text
cargo fmt
cargo test -p generalise-ltd --all-features --lib --bins -- --test-threads=1
39 passed, 0 failed
cargo test -p generalise-ltd --no-default-features --lib -- --test-threads=1
36 passed, 0 failed
RUSTFLAGS=-Dwarnings cargo check -p gammalooprs --tests
passed
cargo clean
removed 42.8 GiB after the full gate initially failed with "No space left on device"
just test_gammaloop
1008 tests run: 1008 passed, 128 skipped
Nextest run ID 7ac98bfc-1a86-4f1c-887d-e5eaab38e4fc
just clippy
passed
```

### In Progress 2026-04-29 14:26 CEST

- Starting the standalone evaluator/CLI slice.
- The first evaluator implementation is an eager f64 evaluator, kept behind the
  optional `eval`/`bin` path. It is intended for CLI validation, parity tests,
  and development diagnostics, not for GammaLoop's production evaluator path.
- Supported numerator syntax will initially match the prototype tests that use
  edge energies and dot products: constants, `edges[i][j]`, `loops[i][j]`,
  `ext[i][j]`, `dot(a,b)`, unary signs, `+`, `-`, `*`, `/`, and integer powers
  with `**`.

### Progress 2026-04-29 15:03 CEST

- Added the feature-gated `eval` module and wired it into the `bin` feature.
- Implemented an eager f64 expression evaluator for standalone validation. It
  evaluates the generated 3D expression over deterministic or user-provided
  external momenta, loop three-momenta, resolved mass values, and optional
  uniform sampling scale `M`.
- The evaluator follows the generated orientation structure directly:
  half-edge prefactors, uniform-scale powers, energy substitutions, linear
  surfaces, numerator-surface factors, denominator trees, variant prefactors,
  and numerator-call caching are all handled in the same pass.
- Added CLI commands:
  - `evaluate` evaluates a serialized generated expression JSON against a DOT
    graph and a numerator expression.
  - `test-cff-ltd` builds CFF and LTD expressions from one DOT input and
    compares their local f64 values for the same generated sample.
- Added lightweight numerator parsing for CLI diagnostics:
  constants, `edges[i][j]`, `loops[i][j]`, `ext[i][j]`, vector references,
  `dot(a,b)`, unary signs, `+`, `-`, `*`, `/`, and integer powers with `**`.
- Smoke-tested the CLI:

```text
generalise-ltd build --representation cff --dot box.dot --json-out /tmp/gl_eval_box_cff.json
generalise-ltd evaluate --orientation-json /tmp/gl_eval_box_cff.json \
  --dot box.dot \
  --numerator-expr 'dot(edges[0], ext[0]) + edges[1][0]'

value = -0.028314663517129336
numerator_calls = 14

generalise-ltd test-cff-ltd --dot box.dot --numerator-expr '1'

cff = 0.08668240169697802
ltd = 0.0866824016920873
|cff-ltd| = 4.890726712503124e-12

generalise-ltd test-cff-ltd --dot box.dot \
  --energy-degree-bounds 0:2 \
  --numerator-expr 'edges[0][0]**2'

cff = -0.005037506621746662
ltd = -0.005037506634835154
|cff-ltd| = 1.3088492095691961e-11
```

- Focused verification:

```text
cargo fmt
RUSTFLAGS=-Dwarnings cargo check -p generalise-ltd --features bin
cargo test -p generalise-ltd --features bin --lib --bins -- --test-threads=1
38 passed, 0 failed
cargo test -p generalise-ltd --no-default-features --lib -- --test-threads=1
36 passed, 0 failed
RUSTFLAGS=-Dwarnings cargo check -p gammalooprs --tests
passed
just clippy
passed
```

### Progress 2026-04-29 15:23 CEST

- Added a `validate` CLI command that emits graph shape and validation JSON and
  can optionally fail with `--expect-valid`.
- Added `test` as an alias for the current `test-cff-ltd` diagnostic command.
- Corrected the diagnostic comparison path so common `--energy-degree-bounds`
  are passed to both CFF and derivative-free LTD generation.
- Added CLI integration coverage behind the `bin` feature:
  - `validate` reports the expected box shape and marks the noisy example
    invalid.
  - `build --numerator-expr auto` records numerator metadata, and
    `evaluate --numerator-expr auto` consumes those JSON energy bounds.
  - the `test` alias compares a bounded non-repeated CFF/LTD box case.
- Focused verification:

```text
cargo fmt
cargo test -p generalise-ltd --features bin --test cli -- --test-threads=1
3 passed, 0 failed
RUSTFLAGS=-Dwarnings cargo check -p generalise-ltd --features bin
passed
cargo test -p generalise-ltd --features bin --lib --bins --tests -- --test-threads=1
41 passed, 0 failed
```

- Follow-up noted: the current standalone f64 evaluator still shows a local
  mismatch for the repeated `box_pow3.dot` bounded CFF versus derivative-free
  LTD diagnostic comparison:

```text
box_pow3.dot, bounds 0:1,1:1,2:0,3:4,
numerator edges[0][0] * edges[1][0] * edges[3][0]**4:
  cff = -0.0040451133544909465
  ltd = -0.0044228533349723875
  |cff-ltd| ~= 3.8e-4
```

  This is not asserted as a passing test yet. The bounded CFF side was already
  checked against the Python prototype through the temporary adapter, so the
  next investigation should focus on the Rust eager evaluator's repeated LTD
  variant/energy-map handling or on the exact diagnostic mode being compared.

### Progress 2026-04-29 15:34 CEST

- Updated the shipped example scripts to exercise the Rust `generalise-ltd`
  binary rather than stale Python-prototype commands.
- The one-loop runtime comparison scripts now build LTD/CFF JSONs, evaluate
  them with the feature-gated eager f64 evaluator, and run a lightweight
  `test-cff-ltd` scalar comparison.
- The five-loop script now builds the no-repeat and repeated structure JSONs
  and runs a lightweight eager f64 diagnostic on the no-repeat graph. It still
  explicitly notes that the compiled Symbolica runtime path is not ported.
- The `box_pow3` showcase now builds/evaluates bounded CFF and evaluates the
  uniform-scale LTD structure with an explicit `M`.
- Focused verification:

```text
GENERALISED_LTD_DEMO_DIR=/tmp/gl_box_showcase \
  crates/generalise-ltd/examples/scripts/box_pow3_cli_showcase.sh
passed
cargo test -p generalise-ltd --features bin \
  validator::tests::copied_runtime_scripts_are_executable -- --exact
passed
GENERALISED_LTD_10EXT_DIR=/tmp/gl_10ext \
  crates/generalise-ltd/examples/scripts/one_loop_10_external_runtime_compare.sh
passed, with |cff-ltd| ~= 3.9e-13 for the scalar diagnostic
cargo test -p generalise-ltd --all-features --lib --bins --tests -- --test-threads=1
43 passed, 0 failed
cargo test -p generalise-ltd --no-default-features --lib -- --test-threads=1
36 passed, 0 failed
```

### Progress 2026-04-29 15:45 CEST

- Added `crates/generalise-ltd/README.md` documenting the crate ownership,
  GammaLoop-facing library API, feature flags, standalone CLI commands, and the
  current bounded-CFF support limits.
- Added `compare` as a CLI alias for the current CFF/LTD diagnostic comparison,
  alongside the `test` alias.
- Added eager-evaluator coverage that expressions using the runtime uniform
  sampling scale `M` reject missing `M`, and that CLI parsing rejects `M=0`.
- Added public DOT convenience entry points:
  `generate_3d_expression_from_dot_graph(...)` and
  `generate_3d_expression_from_dot_path(...)`, both using the same linnet DOT
  parser path as the CLI before dispatching to
  `generate_3d_expression_from_parsed(...)`.
- Focused verification:

```text
cargo fmt
RUSTFLAGS=-Dwarnings cargo check -p generalise-ltd --features bin
passed
cargo test -p generalise-ltd --features bin \
  generation::ltd_tests::dot_path_api_builds_cff_box_expression -- --exact
passed
cargo test -p generalise-ltd --features bin \
  eval::tests::evaluator_requires_nonzero_uniform_scale_when_expression_uses_m -- --exact
passed
cargo test -p generalise-ltd --all-features --lib --bins --tests -- --test-threads=1
45 passed, 0 failed
cargo test -p generalise-ltd --no-default-features --lib -- --test-threads=1
37 passed, 0 failed
RUSTFLAGS=-Dwarnings cargo check -p gammalooprs --tests
passed
```

### Verification 2026-04-29 16:02 CEST

- Re-ran the full repo gates after the CLI/evaluator/example-script updates and
  after adding the DOT convenience API.
- Cleaned up stale generation error text so unsupported sectors report the
  current conservative support boundary instead of claiming repeated LTD is
  unimplemented.

```text
git diff --check
passed
just clippy
passed
just test_gammaloop
1009 tests run: 1009 passed, 128 skipped
Nextest run ID 82fe9a47-9999-4231-88a7-46c81e62bedb
RUSTFLAGS=-Dwarnings cargo check -p generalise-ltd --features bin
passed
cargo test -p generalise-ltd --features bin --lib --tests -- --test-threads=1
43 passed, 0 failed
cargo clean
removed 35.6 GiB after final verification to recover workspace disk headroom
```

### Refactor Plan 2026-04-29 Symbolica-native `three-dimensional-reps`

- The generalized representation crate will be renamed from
  `generalise-ltd` / `generalise_ltd` to `three-dimensional-reps` /
  `three_dimensional_reps`.
- The standalone binary will be split into a separate workspace package named
  `three-dimensional-reps-cli` / `three_dimensional_reps_cli`. This package may
  depend on both `three-dimensional-reps` and `gammalooprs`, avoiding a cyclic
  dependency while still allowing the CLI to parse GammaLoop DOT graphs through
  GammaLoop's rich graph importer.
- The core library API must be graph-first. `DotGraph` and the current
  `ParsedGraph` path are no longer the common entry point for
  `generate_3d_expression`; they are CLI/test conveniences only. The
  GammaLoop-facing generator should consume the already-parsed graph data via
  traits implemented in `gammalooprs`.
- All symbolic/algebraic manipulation in the new implementation path must use
  Symbolica-native objects: `Atom`, `AtomView`, Symbolica rational
  coefficients, Symbolica/numerica matrices, `AtomCore::solve_linear_system`,
  `Matrix::{det, rank, row_reduce, solve}`, and Symbolica GCD primitives.
- The migrated `old_cff` path remains a comparison anchor and should not be
  modernized beyond the minimum adapters required by shared data types.
- `CFFExpression` / `ThreeDExpression` must keep bincode support using the
  Symbolica state map as decode context, while serde/JSON should render any
  `Atom` payloads as canonical strings.
- All symbols used by `three-dimensional-reps` must be registered in a lazy
  static symbol collection. Symbols intended to interoperate with existing
  crates must use their exact existing namespaces, for example spenso concrete
  indices such as `spenso::cind`; newly introduced internal symbols use the
  `three_dimensional_reps` namespace.

### Progress 2026-04-29 Rename and CLI split

- Renamed the core package directory to `crates/three-dimensional-reps` and the
  Rust library crate to `three_dimensional_reps`.
- Updated `gammalooprs` imports and dependency wiring to use
  `three-dimensional-reps`.
- Added the separate `crates/three-dimensional-reps-cli` package with the
  `three-dimensional-reps` binary, leaving the core package as a library-only
  crate.
- Updated the example scripts and `just test_gammaloop` package list to target
  `three-dimensional-reps-cli`.
- Verified the structural split with:

```text
cargo fmt
RUSTFLAGS=-Dwarnings cargo check -p three-dimensional-reps -p three-dimensional-reps-cli
passed
```

### Progress 2026-04-29 CLI Split Verification

- Verified the renamed standalone CLI package after the
  `three-dimensional-reps` / `three-dimensional-reps-cli` split:

```text
cargo test -p three-dimensional-reps-cli --features old_cff -- --test-threads=1
3 passed
```

### Verification 2026-04-29 Minimal Core Feature Check

- Checked the core crate with default features disabled:

```text
cargo test -p three-dimensional-reps --no-default-features --lib -- --test-threads=1
37 passed
```

### Verification 2026-04-29 Full GammaLoop Test Selection

- Ran the repository sanity gate with the renamed core and CLI crates included
  in the `just test_gammaloop` package list:

```text
just test_gammaloop
1015 tests run: 1015 passed, 128 skipped
```

### Progress 2026-04-29 Additional Graph Signature Test Port

- Ported more isolated Python graph-signature tests:
  - two-loop Symbolica-style propagator extraction,
  - unsatisfiable signature reconstruction error handling,
  - direct parsed-graph reconstruction preserving internal-edge mass keys.
- I also checked the Python `vakint` total-valence cap rejection test. The Rust
  reconstruction currently applies `max_degree` to the internal placement search
  and does not yet reproduce the Python test's total-valence semantics for
  vakint output, so I did not add that as a passing test. This remains an
  explicit graph-signature parity gap.
- Verified the graph-signature test subset:

```text
cargo test -p three-dimensional-reps graph_signatures::tests -- --test-threads=1
8 passed
```

### Verification 2026-04-29 Core Tests After Graph Signature Ports

```text
cargo fmt
RUSTFLAGS=-Dwarnings cargo check -p three-dimensional-reps -p three-dimensional-reps-cli -p gammalooprs
passed

cargo test -p three-dimensional-reps --features diagnostics,old_cff,test-support --lib --tests -- --test-threads=1
42 passed
```

### Progress 2026-04-29 Affine Coefficient Refactor Start

- The remaining symbolic-native gap in the current expression layer is
  `LinearEnergyExpr`: scalar prefactors and matrix solves now use Symbolica
  rationals, but affine energy maps still store edge and `M` coefficients as
  native `i64`.
- I am changing those coefficients to the existing Symbolica-backed
  `RationalCoefficient` wrapper while keeping edge identifiers as Rust indices.
  This preserves the current compact structured representation needed by the
  evaluator and display code, but removes native integer arithmetic from the
  affine coefficient layer and allows rational weights for OSE, external-energy,
  and numerator-scale terms.
- This is intentionally not applied to the legacy `old_cff` generation logic
  beyond the adapters needed to keep it compiling as a comparison anchor.

### Progress 2026-04-29 Affine Coefficient Refactor

- Changed `LinearEnergyExpr` internal-edge, external-energy, and numerator-scale
  coefficients from native `i64` to `RationalCoefficient`, the local wrapper
  around Symbolica's rational coefficient type.
- Added rational-coefficient constructors for OSE, external energy, and
  numerator-scale terms so generalized surfaces can carry non-integral weights
  without conversion back through fixed-width integers.
- Updated display, Atom conversion, evaluator evaluation, coefficient-map
  canonicalization, and GammaLoop surface substitution adapters to consume the
  Symbolica-backed coefficients directly.
- Left the legacy `old_cff` generation algorithm structurally unchanged; it now
  only converts its existing integral coefficients through the new
  `LinearEnergyExpr` adapters.
- Verified after formatting with:

```text
cargo fmt
RUSTFLAGS=-Dwarnings cargo check -p three-dimensional-reps -p three-dimensional-reps-cli -p gammalooprs
passed
```

### Verification 2026-04-29 Projected Graph Cleanup

- Re-ran the focused tests after removing the DOT payload and projected
  subproblem DOT roundtrip:

```text
cargo test -p three-dimensional-reps --features diagnostics,old_cff,test-support --lib --tests -- --test-threads=1
39 passed

cargo test -p three-dimensional-reps-cli --features old_cff -- --test-threads=1
3 passed
```

### Progress 2026-04-29 Projected Graph Reconstruction Cleanup

- Replaced the projected CFF lower-sector DOT roundtrip with a direct
  `ParsedGraph` reconstruction helper.
- The graph-from-signatures CLI can still render DOT for standalone exploration,
  but the generator no longer builds a DOT string and parses it back for the
  projected subproblem.
- Verified after formatting with:

```text
cargo fmt
RUSTFLAGS=-Dwarnings cargo check -p three-dimensional-reps -p three-dimensional-reps-cli -p gammalooprs
passed
```

### Verification 2026-04-29 Symbolica-Rational Affine Layers

- Re-ran the focused core and CLI tests after changing the affine coefficient
  storage:

```text
cargo test -p three-dimensional-reps --features diagnostics,old_cff,test-support --lib --tests -- --test-threads=1
39 passed

cargo test -p three-dimensional-reps-cli --features old_cff -- --test-threads=1
3 passed
```

### Progress 2026-04-29 Parsed Graph DOT Payload Removal Start

- `ParsedGraph` still stores an optional `DotGraph` from the first port, which
  makes DOT look like the shared generator representation even for GammaLoop's
  native `Graph` adapter.
- I am removing that payload so `DotGraph` remains only a CLI/test parser input.
  The generator will carry the compact parsed edge/node/signature data only.

### Progress 2026-04-29 Parsed Graph DOT Payload Removal

- Removed the optional `DotGraph` payload from `ParsedGraph`.
- Updated all contracted/projected subgraph constructors and the GammaLoop
  `Graph` adapter to build `ParsedGraph` from parsed edge/node/signature data
  only.
- DOT parsing still exists as a CLI/test input path through `parse_dot_graph`
  and `ThreeDGraphSource for DotGraph`, but the generator no longer propagates a
  DOT graph through its common representation.
- Verified after formatting with:

```text
cargo fmt
RUSTFLAGS=-Dwarnings cargo check -p three-dimensional-reps -p three-dimensional-reps-cli -p gammalooprs
passed
```

### Progress 2026-04-29 Known Affine Expression Refactor

- Converted the internal generalized-generator `KnownLinearExpr` normal form
  from native integer coefficients to Symbolica rationals as well.
- Removed the integral-coefficient restriction from rational scaling of known
  factors; sampled-variable, OSE, external-energy, constant, and numerator-scale
  pieces now stay in Symbolica's rational domain throughout this layer.
- Added a rational scaling method on `LinearEnergyExpr` so known-factor
  substitution no longer needs to round or reject non-integral affine
  coefficients when turning sampled-variable expressions into surface
  expressions.
- Verified after formatting with:

```text
cargo fmt
RUSTFLAGS=-Dwarnings cargo check -p three-dimensional-reps -p three-dimensional-reps-cli -p gammalooprs
passed
```

### Progress 2026-04-29 GammaLoop Graph Source Adapter

- Added a `ThreeDGraphSource` implementation for `gammalooprs::graph::Graph`.
- The adapter extracts node connectivity, edge momentum signatures from
  `loop_momentum_basis`, edge labels, and mass keys directly from GammaLoop's
  parsed graph structures. It does not serialize the graph back to DOT.
- The extracted `ParsedGraph` deliberately carries `graph: None`; DOT is now an
  optional diagnostic/CLI payload rather than a required library input.
- Split edges are currently reported as a clean source-extraction error for the
  full-graph adapter; explicit subgraph/cut adapters should be added separately
  when the generalized generator needs that topology.
- Verified with:

```text
cargo fmt
RUSTFLAGS=-Dwarnings cargo check -p three-dimensional-reps -p gammalooprs
passed
```

### Verification 2026-04-29 Core Generalized Tests

- Ran the current Rust port tests for the core crate, including diagnostics,
  `old_cff`, and test-support paths:

```text
cargo test -p three-dimensional-reps --features diagnostics,old_cff,test-support --lib --tests -- --test-threads=1
39 passed
```

### Progress 2026-04-29 Graph-Source Generator Entry

- Introduced `ThreeDGraphSource` in `three-dimensional-reps::graph_io`.
- Changed `generate_3d_expression` to consume a `ThreeDGraphSource` instead of
  a bare generic `HedgeGraph` placeholder.
- Implemented the trait for `ParsedGraph` and for `DotGraph` as a CLI/test
  convenience. This leaves a direct path for `gammalooprs::graph::Graph` to
  implement the trait without converting back to DOT.
- Made `ParsedGraph` carry `Option<DotGraph>` so future GammaLoop-native
  sources can construct parsed input without requiring a DOT graph payload.
- Verified after formatting with:

```text
cargo fmt
RUSTFLAGS=-Dwarnings cargo check -p three-dimensional-reps -p three-dimensional-reps-cli
passed
```

### Progress 2026-04-29 Expression Solve Refactor

- Replaced the hand-written `solve_expr_system_unimodular` Gaussian elimination
  with calls to the Symbolica-backed `solve_rational_system` helper.
- Added a narrow `scale_linear_energy_expr_rational` adapter that applies
  Symbolica rational solution coefficients to the current linear energy
  expression storage while still enforcing integral coefficients for edge and
  uniform-scale terms.
- Verified after formatting with:

```text
cargo fmt
RUSTFLAGS=-Dwarnings cargo check -p three-dimensional-reps -p three-dimensional-reps-cli
passed
```

### Progress 2026-04-29 Symbolica Rational Surface Coefficients

- Replaced `RationalCoefficient { numerator, denominator }` with a thin wrapper
  around Symbolica's rational type.
- Removed the local `gcd_i64` normalization path; all coefficient reduction now
  comes from Symbolica's rational field.
- Updated expression/evaluator plumbing to borrow/clone rational coefficients
  explicitly.
- `rational_to_coefficient` now preserves arbitrary-size Symbolica rational
  values instead of requiring conversion back to `i64`.
- Verified this change with:

```text
RUSTFLAGS=-Dwarnings cargo check -p three-dimensional-reps -p three-dimensional-reps-cli
passed
```

### Progress 2026-04-29 Symbolica Algebra Refactor Start

- Beginning the first mechanical refactor pass that removes the custom exact
  rational implementation and the hand-written determinant/rank/linear-solve
  routines from `three-dimensional-reps`.
- The first pass keeps the current module boundaries and helper names where
  practical, but changes their internals to Symbolica/numerica primitives. This
  limits blast radius while making later Atom-backed surface and affine-map
  refactors easier to review.
- Targeted replacements for this pass:
  - `utils::Rational` becomes Symbolica's
    `symbolica::domains::rational::Rational`.
  - `rank_i64`, `rank_rational`, and `solve_rational_system` delegate to
    `symbolica::tensors::Matrix` over the rational field.
  - determinant sign/non-zero tests delegate to Symbolica matrices rather than
    recursive integer expansion.
  - combinatorial helpers return Symbolica rational coefficients instead of
    native integer products.

### Progress 2026-04-29 Symbolica Rational/Matrix Adapter Pass

- Replaced the custom `utils::Rational` struct with Symbolica's
  `symbolica::domains::rational::Rational`.
- Reimplemented rank and solve helpers through `symbolica::tensors::Matrix`
  over Symbolica's rational field.
- Replaced recursive determinant expansion in `cut_structure` and the remaining
  generation basis checks with Symbolica matrix determinant queries.
- Updated derivative/interpolation coefficient plumbing to clone Symbolica
  rational atoms explicitly instead of relying on `Copy` semantics.
- Verified this adapter pass with:

```text
RUSTFLAGS=-Dwarnings cargo check -p three-dimensional-reps -p three-dimensional-reps-cli
passed
```

### Progress 2026-04-29 Symbol registry groundwork

- Replaced the ad hoc symbol constructors in `three-dimensional-reps` with a
  lazy `ThreeDimensionalRepSymbols` registry and `SYMBOL_REGISTRY`.
- Existing GammaLoop-facing symbols now use explicit `gammalooprs::...`
  namespaces; the concrete index helper now uses `spenso::cind`; internal
  placeholders use the `three_dimensional_reps` namespace.
- Added `StringSerializedAtom`, a bincode-compatible `Atom` wrapper decoded
  with a `symbolica::state::HasStateMap` context and serialized to JSON as a
  canonical string.
- Verified the registry/helper changes with:

```text
cargo fmt
RUSTFLAGS=-Dwarnings cargo check -p three-dimensional-reps -p three-dimensional-reps-cli
passed
```

### Current Checkpoint 2026-04-29

- The crate rename requested by the user is now in place:
  `three-dimensional-reps` for the core library crate and
  `three-dimensional-reps-cli` for the standalone CLI package.
- The core expression layer uses Symbolica rationals for scalar prefactors,
  affine `LinearEnergyExpr` coefficients, and internal known-factor affine
  coefficients. Fixed-width integer coefficients remain only for graph/signature
  incidence data and edge indices.
- `ParsedGraph` no longer stores a `DotGraph`; DOT is only an input/rendering
  path for CLI and tests. The projected CFF lower-sector path now reconstructs a
  `ParsedGraph` directly instead of writing and reparsing DOT.
- The GammaLoop `Graph` adapter implements `ThreeDGraphSource` and extracts
  loop/external signatures, masses, nodes, and edge directions directly from the
  parsed GammaLoop graph.
- Additional graph-signature tests have been ported. A remaining parity gap is
  the Python `vakint` total-valence cap behavior: Rust currently applies
  `max_degree` to the internal placement search and does not yet reproduce that
  exact vakint total-valence rejection.
- Python-test-suite port status: incomplete. The Rust crates currently have 42
  core tests and 3 CLI tests passing, plus the full example graph/script set
  copied into `crates/three-dimensional-reps/examples` and covered by validation
  tests. The Python prototype has 73 pytest tests in `tests/test_cli.py` and 12
  unittest-style graph-signature tests in `src/graph_signatures.py`; many of the
  large numerical comparison, compiled-symbolica evaluator, all-example
  parametrized, split-mass, and slow/stress tests still need Rust equivalents.
- Latest verification:

```text
cargo fmt
RUSTFLAGS=-Dwarnings cargo check -p three-dimensional-reps -p three-dimensional-reps-cli -p gammalooprs
passed

cargo test -p three-dimensional-reps --features diagnostics,old_cff,test-support --lib --tests -- --test-threads=1
42 passed

cargo test -p three-dimensional-reps-cli --features old_cff -- --test-threads=1
3 passed

cargo test -p three-dimensional-reps --no-default-features --lib -- --test-threads=1
37 passed

just test_gammaloop
1018 tests run: 1018 passed, 128 skipped
```

### Revised Plan 2026-04-29 GammaLoop 3Drep command

- The standalone `three-dimensional-reps-cli` package is no longer the right
  shape. The next refactor is to remove that package and move its useful command
  surface into `crates/gammaloop-api/src/commands/threedreps` as a new
  GammaLoop command named `3Drep`.
- `three-dimensional-reps` remains the core library crate. It should stay free
  of GammaLoop API/evaluator dependencies and expose graph-source traits,
  expression construction, expression data structures, display helpers, and
  optional diagnostic/old_cff support.
- The new `3Drep` command should obtain graphs through GammaLoop's existing
  state/model/process/integrand machinery, using the same `-p`, `-i`, and graph
  selection conventions as existing commands. Direct DOT parsing in the
  three-dimensional-reps layer is no longer acceptable for GammaLoop-facing CLI
  flows.
- Standalone evaluator support should be implemented in the GammaLoop API
  command layer, using existing GammaLoop `Evaluator` and `ParamBuilder`
  infrastructure. It is diagnostic/testing-only and must be able to evaluate an
  oriented JSON 3D expression together with the selected graph/model context.
- A robust multi-way comparison harness is the main validation target. It should
  drive mass-shift pure-LTD diagnostics, LTD, CFF, and multiple numerator-scale
  `M` choices from the same selected graph and numerator setup, then report
  total and focused discrepancies for the `3Drep test-cff-ltd` command and Rust
  integration tests.
- Pure LTD remains exploratory/testing-only; no derivative evaluator is planned.
- The old CFF implementation remains an anchor for parity tests and should not
  be structurally refactored beyond adapters needed by the surrounding crate
  moves.

### Revised Decisions 2026-04-29 Fully integrated 3Drep

- `3Drep` is now fully integrated into the GammaLoop CLI ecosystem. It should
  only operate on the active state and on already-loaded/generated integrands;
  there is no standalone model/DOT ingestion mode anymore.
- The command may expose aliases such as `3drep` and `threedreps`, but the
  canonical command name remains `3Drep`.
- Graph selection is always for one individual graph at a time, by graph index
  or graph name. Graph-group selection is explicitly out of scope for 3D
  representation commands.
- `graph-from-signatures` remains as an orthogonal subcommand for now. It takes
  a signature/numerator string, emits a GammaLoop-loadable DOT graph string, and
  should be tested by reloading the emitted graph and checking that edge
  signatures and numerator information round-trip.
- The numerator interpolation scale `M` should be read from
  `default-runtime-settings.general`, at the same level as `m_uv` and `mu_r`.
  The comparison harness should test the active value, twice that value, and
  half that value.
- Heavy and integration-like 3D representation tests belong in
  `tests/tests/test_threedreps.rs` with fixtures under
  `tests/resources/graphs/threedreps/`, not in `gammaloop-api` unit tests.
- `three-dimensional-reps` should be purged of user harness/parsing concerns
  except for minimal library support needed by `generate_3d_expression` and the
  retained graph-from-signatures generator.

### Revised Decisions 2026-04-29 Scale and 3Drep workspace

- `default-runtime-settings.general.numerator_interpolation_scale` should be an
  optional value. When it is `None`, GammaLoop should choose reconstruction
  strategies that do not use that scale, and `3Drep test-cff-ltd` should skip
  comparisons over scaled `M` variants.
- `graph-from-signatures` should only reconstruct graph signatures. It should
  not accept or propagate a numerator; emitted GammaLoop DOT should set the
  global numerator to `"1"`.
- Evaluation-oriented 3Drep subcommands must use a separate artifact workspace,
  defaulting to `threed_workspace` under the active state when writable and to a
  cwd path with state-name suffix under read-only-state mode, following the
  integration workspace convention. This workspace is for generated evaluator
  files, oriented JSON output, Symbolica expressions, parameters, and other
  inspectable artifacts; it must not be part of the GammaLoop state and must not
  mutate the state.
- Reusing GammaLoop evaluator machinery means reusing the builder/evaluator
  structs and backend infrastructure to create diagnostic 3Drep evaluators, not
  reusing already-generated state integrand evaluators.

### Revised Decisions 2026-04-29 Numerator energy degree caps

- `3Drep` generation should accept an optional user override map for numerator
  energy-power caps per edge. Missing entries should be filled from automatic
  analysis of the selected graph numerator.
- The automatic degree analysis belongs in `gammalooprs`, not in the
  `three-dimensional-reps` crate, because it is useful beyond the 3Drep command
  and should operate directly on GammaLoop's parsed numerator and graph data.
- The analyzer must use Symbolica and avoid expanding the numerator in all
  variables. It should compute weighted energy degrees by traversing the
  Symbolica atom structure, using max over sums and additive degrees over
  products/powers.
- The analyzer should support both EMR numerator dependence through temporal
  `Q(edge, 0)` components and LMB numerator dependence through temporal
  `K(lmb, 0)` components, mapping LMB loop energies back to the corresponding
  LMB-carrier edges when building the per-edge cap map.
- `3Drep` should expose CLI overrides for any set of edges; overrides replace
  the automatically deduced value for those edges only.
- `build` should support options to suppress JSON saving and/or suppress pretty
  terminal output, while preserving the default build-to-workspace chain used by
  `evaluate`.

### Clarification 2026-04-29 Energy index counting

- Explicitly spatial vectors such as `Q3(edge, ...)` do not contribute to
  energy-power scaling.
- Four-vector momentum atoms such as `Q(edge, index)` and `K(lmb, index)` do
  contribute one energy power when `index` is either a Lorentz placeholder of
  the form `mink(4, symbol)` or a concrete temporal index `cind(0)`.
- Concrete non-temporal indices such as `cind(n)` with `n > 0` do not contribute
  to energy-power scaling.

### Progress 2026-04-29 GammaLoop-side energy cap analysis

- Added `general.numerator_interpolation_scale: Option<f64>` to runtime general
  settings with default `None`.
- Added `gammalooprs::numerator::energy_degree`, including a Symbolica
  structural energy-degree analyzer and `Graph::numerator_energy_power_caps()`.
- The analyzer walks Symbolica atoms without global expansion: sums take
  edge-wise maxima, products add edge-wise degrees, and positive integer powers
  scale degrees. Non-polynomial energy powers are reported as errors.
- The implementation recognizes temporal EMR and LMB four-vector entries
  according to the clarified `Q/K` + `mink(...)`/`cind(0)` rules and maps LMB
  loop-energy dependence to the corresponding loop-basis carrier edge.
- Verified this step with:

```text
RUSTFLAGS=-Dwarnings cargo check -p gammalooprs
passed

cargo test -p gammalooprs numerator::energy_degree -- --test-threads=1
6 passed
```

### Progress 2026-04-29 Integrated 3Drep command skeleton

- Removed the obsolete standalone `three-dimensional-reps-cli` package from the
  workspace.
- Added `crates/gammaloop-api/src/commands/threedreps` and registered the
  canonical `3Drep` command with aliases `3drep` and `threedreps`.
- The new command resolves graphs only from the active GammaLoop state using
  `-p`, `-i`, and `-g`, where `-g` accepts graph ids, graph names, and inspect
  display labels of the form `#id : name`.
- Implemented `3Drep validate`, `3Drep build`, and `3Drep graph-from-signatures`
  in the integrated command module. `build` uses automatic numerator
  energy-degree caps plus edge-specific overrides, writes oriented JSON to the
  3Drep workspace by default, and supports disabling JSON saving or pretty
  output.
- Added the default 3Drep workspace path convention: writable states use
  `<state>/threed_workspace`; read-only states use a cwd `threed_workspace`
  path with state-name suffix when available.
- `3Drep evaluate` and `3Drep test-cff-ltd` are registered but still return
  explicit errors pending the evaluator-backed adapter and multi-way comparison
  harness.
- Verified this command migration checkpoint with:

```text
RUSTFLAGS=-Dwarnings cargo check -p gammaloop-api
passed
```

### Progress 2026-04-29 Integrated command cleanup and GammaLoop DOT output

- Split graph selection options from JSON/output options in the integrated
  `3Drep` command, so `-p/-i/-g` are now represented by a dedicated selector
  and artifact paths live on the subcommands that actually write artifacts.
- Changed `graph-from-signatures` to request the new GammaLoop-style DOT
  renderer by default instead of the old standalone visualization format.
- Added `ReconstructDotFormat::Gammaloop` in `three-dimensional-reps`. The
  renderer emits graph/node/edge numerator defaults, mass-only scalar edges,
  LMB-style `K(i,a___)`/`P(i,a___)` momentum labels, and `lmb_id` annotations
  for loop variables whose input signature has an explicit unit carrier edge.
- Kept the previous hybrid and vakint DOT renderers available for internal
  roundtrip tests and diagnostics.
- Known limitation: for a signature set that does not contain a pure unit
  carrier edge for every loop variable, the GammaLoop DOT output can preserve
  the human-readable `lmb_rep` label but cannot yet force GammaLoop's parser to
  reconstruct exactly the same loop-basis signatures. Supporting that robustly
  likely requires teaching GammaLoop's DOT parser to consume `lmb_rep` as an
  actual loop-momentum-basis constraint instead of only a serialization label.
- Verified this checkpoint with:

```text
RUSTFLAGS=-Dwarnings cargo check -p three-dimensional-reps
passed

RUSTFLAGS=-Dwarnings cargo check -p gammaloop-api
passed

cargo test -p three-dimensional-reps graph_from_signatures -- --test-threads=1
3 passed
```

### Progress 2026-04-29 Broader current Rust test checkpoint

- Re-ran the three-dimensional-reps crate test suite with the currently relevant
  optional implementation/test features enabled after the integrated command
  cleanup.
- Verified this checkpoint with:

```text
cargo test -p three-dimensional-reps --features diagnostics,old_cff,test-support --lib --tests -- --test-threads=1
43 passed
```

### Next implementation checkpoint 2026-04-29

- The next focus is the GammaLoop graph-source boundary. The current core
  generator works with compact local internal-edge ids, while GammaLoop graphs,
  numerators, parameter builders, and user-facing graph selectors use the
  original linnet `EdgeIndex` values. This can silently mismatch expressions
  whenever paired internal edges are not numbered as `0..n_internal`.
- Planned fix: keep the compact ids inside the pure three-dimensional-reps
  algorithms, but let rich graph sources provide an edge-index remapping. The
  public `generate_3d_expression(graph, options)` entry point will translate
  energy-degree bounds to compact ids before generation and then remap the
  resulting oriented expression, surfaces, half-edge factors, and numerator
  energy maps back to GammaLoop's original internal/external energy ids.

### Progress 2026-04-29 GammaLoop edge-id remapping boundary

- Added `EnergyEdgeIndexMap` and a default hook on `ThreeDGraphSource` so rich
  graph sources can declare how compact internal/external 3D-generation ids map
  back to the source graph's energy ids.
- Updated `generate_3d_expression(graph, options)` to translate
  energy-degree bounds into compact local ids before generation, then remap the
  generated orientations, edge-energy maps, loop-energy maps, half-edge factors,
  and linear surface cache back to the source ids.
- Implemented the hook for GammaLoop `Graph`: internal terms map to original
  linnet `EdgeIndex` values, external energy terms map through the graph's loop
  momentum basis external edges, and orientation vectors are expanded to the
  source graph edge count.
- Added a regression test covering a deliberately nontrivial compact-to-source
  map, including bound translation and external energy remapping.
- Verified this checkpoint with:

```text
cargo fmt --check
passed

RUSTFLAGS=-Dwarnings cargo check -p three-dimensional-reps
passed

cargo test -p three-dimensional-reps rich_graph_source_remaps -- --test-threads=1
1 passed

RUSTFLAGS=-Dwarnings cargo check -p gammaloop-api
passed
```

### Next validation checkpoint 2026-04-29 Simple box numerator-one

- The immediate focus has been narrowed to the simplest common reference case:
  the one-loop box with no repeated propagators and numerator `1`.
- Python reference commands were run in `IGNORE/generalised_ltd` for `cff`,
  old `hybrid` (now `ltd`), and the numerical multi-way test. The expected
  values are:
  - CFF/LTD numerical value:
    `0.5132881569466608530785378140887218534301002662092683450316784064319557537623902`
    for the deterministic Python sample.
  - CFF structure: 14 orientations, 12 surfaces, one variant per orientation,
    half-edges `[0, 1, 2, 3]`, and orientation labels
    `-+++`, `+-++`, `--++`, `++-+`, `-+-+`, `+--+`, `---+`, `+++-`,
    `-++-`, `+-+-`, `--+-`, `++--`, `-+--`, `+---`.
  - LTD structure: 4 orientations `+xxx`, `x+xx`, `xx+x`, `xxx+`, one
    half-edge per orientation, and prefactor `-1` for each orientation.
- Plan for this checkpoint: first restore a clean Rust compile state for the
  in-progress `3Drep` harness, then add a focused integration test in
  `tests/tests/test_runs/test_3d_reps.rs` that exercises the current Rust
  generation/evaluation path for this box reference and compares it against the
  Python deterministic result.

### Progress 2026-04-29 Simple box numerator-one validation

- Added the first focused Rust parity test at
  `tests/tests/test_runs/test_3d_reps.rs`.
- The test builds the simple box CFF and LTD expressions from the same dot graph
  used by the Python reference, checks the expected orientation labels,
  half-edge factors, CFF/LTD e-surface denominator collections, and evaluates
  both representations at the deterministic Python sample.
- Fixed the migrated pure-CFF generator shape: multiple CFF surface chains for a
  single orientation are now represented as one variant with a branching
  denominator tree, matching the old GammaLoop/Python oriented structure instead
  of incorrectly materializing one variant per chain.
- Verified the focused test with:

```text
cargo test -p gammaloop-integration-tests --test test_runs simple_box_numerator_one_matches_python_cff_and_ltd_reference -- --nocapture
1 passed
```
- Updated the `test_gammaloop` just recipe package list after the
  `three-dimensional-reps-cli` removal so that the selected test profile can run
  the integrated `3Drep` code path.
- Verified the broader checkpoint with:

```text
cargo fmt --check
passed

RUSTFLAGS=-Dwarnings cargo check -p gammaloop-api
passed

RUSTFLAGS=-Dwarnings cargo check -p gammaloop-integration-tests --test test_runs
passed

cargo test -p three-dimensional-reps --features diagnostics,old_cff,test-support,eval --lib --tests -- --test-threads=1
46 passed

just test_gammaloop
1024 tests run: 1024 passed, 128 skipped
```

### Next validation checkpoint 2026-04-29 Repeated box numerator-one

- The next focused target is `box_pow3.dot`: the one-loop box where the last
  channel is represented by three repeated propagator edges.
- Python reference commands were run in `IGNORE/generalised_ltd` for `cff`,
  old `hybrid` (now public `ltd`), pure `ltd`, and the numerical test with
  numerator `1` and zero energy-degree bounds.
- Expected numerical targets from the Python deterministic sample:
  - CFF and derivative-free LTD/hybrid:
    `0.33810887306602043713023885730518228244342876238882636557147668349464819846425127`.
  - Pure-LTD distinct-mass split proxy sequence for epsilons
    `0.1, 0.05, 0.025, 0.0125`:
    `0.34211618742478167`, `0.33910608053087872`,
    `0.33835788637720207`, `0.33817110836350058`.
- Expected high-level structures:
  - CFF: 6 internal edges, 27 e-surfaces, 62 orientations, one variant per
    orientation with half-edges `[0, 1, 2, 3, 4, 5]`.
  - LTD/hybrid: 24 linear surfaces, 4 orientations, labels `+xxxxx`, `x+xxxx`,
    `xx+xxx`, `xxx+xx`; the repeated channel orientation has three variants
    with prefactors `-1`, `-3`, `-6` and increasing powers of the repeated
    representative half-edge.
  - Pure LTD diagnostic/mass-shift structure: 57 linear surfaces and 6 cut
    orientations. It is only used here through distinct-mass shifts as a
    diagnostic proxy, not as the repeated-propagator representation.

### Progress 2026-04-29 Repeated box numerator-one validation

- Extended `tests/tests/test_runs/test_3d_reps.rs` with a `box_pow3.dot`
  parity case using numerator `1` and zero energy-degree bounds.
- The test checks:
  - CFF structure counts, orientation labels, prefactors, and half-edge factors.
  - Repeated-propagator LTD structure counts and the three variants of the
    repeated-channel pseudo-orientation, including prefactors `-1`, `-3`, `-6`.
  - Pure-LTD diagnostic structure counts for the distinct-mass split proxy.
  - CFF/LTD numerical agreement with the Python reference at the deterministic
    sample.
  - The pure-LTD mass-shift sequence for epsilons
    `0.1, 0.05, 0.025, 0.0125`, including monotonic convergence to the CFF/LTD
    value.
- The repeated-channel LTD f64 evaluator differs from the high-precision Python
  reference by about `7.4e-12`; the test therefore uses a `1e-10` tolerance for
  the repeated-box LTD and pure-LTD split numerical comparisons while keeping
  structural checks exact.
- While running the full `dev-optim` nextest profile, the repeated-channel
  variants exposed a nondeterministic order from a `HashMap` used in denominator
  derivative accumulation. This was fixed by using a `BTreeMap` there and by
  deriving ordering for linear/hybrid surface ids, so the oriented structure is
  stable across profiles.
- Verified the focused test with:

```text
cargo test -p gammaloop-integration-tests --test test_runs repeated_box_numerator_one_matches_python_cff_ltd_and_mass_shift_references -- --nocapture
1 passed
```
- Verified the combined and broader checkpoints with:

```text
cargo test -p gammaloop-integration-tests --test test_runs test_3d_reps -- --nocapture
2 passed

cargo fmt --check
passed

RUSTFLAGS=-Dwarnings cargo check -p gammaloop-api
passed

RUSTFLAGS=-Dwarnings cargo check -p gammaloop-integration-tests --test test_runs
passed

cargo test -p three-dimensional-reps --features diagnostics,old_cff,test-support,eval --lib --tests -- --test-threads=1
46 passed

just test_gammaloop
1025 tests run: 1025 passed, 128 skipped
```

### Progress 2026-04-29 Old Python Test Matrix Through GammaLoop CLI

- Added GammaLoop-format equivalents of the old Python
  `IGNORE/generalised_ltd/examples/graphs` fixtures under
  `tests/resources/graphs/threedreps/`. These fixtures are intentionally not
  copied verbatim from the Python repository: they use GammaLoop graph syntax,
  explicit `node [num="1"]` defaults, scalar-model-compatible masses, and are
  only loaded through the normal `import graphs` command in tests.
- Added import-shape coverage for the old examples:
  - `box.dot`
  - `box_pow3.dot`
  - `sunrise_pow4.dot`
  - `kite_double_nested_repeats.dot`
  - `kite_nested_repeats.dot`
  - `kite_sandwich_repeats.dot`
  - `proper_iterated_sandwiched_bubble.dot`
  - `mercedes_multi_repeats.dot`
  - `four_loop_stress.dot`
  - `five_loop_no_repeats.dot`
  - `one_loop_10_external.dot`
  - `one_loop_15_external.dot`
  - the five `five_loop_ultimate_basis*.dot` slow fixtures.
- Added a CLI-only numerator energy-power reconstruction test. The test injects
  actual GammaLoop numerator expressions into imported DOT fixtures and checks
  that `3Drep test-cff-ltd` records the same automatic edge-energy caps as the
  intended old Python bounds, including both EMR-style `Q(edge,cind(0))` powers
  and LMB-style `K(loop,cind(0))` powers.
- `3Drep test-cff-ltd` now records generation failures per representation case
  in `test_cff_ltd_manifest.json` instead of aborting the whole diagnostic
  command. This lets the harness inspect graph metadata and automatic numerator
  bounds even when one representation sector is not implemented yet.
- Added `graph-from-signatures` closure tests through the GammaLoop CLI. The
  generated DOT is immediately loaded back with `import graphs` and the parsed
  signatures are compared to the requested signatures. The non-trivial closure
  case is a three-loop four-point Mercedes-like graph.
- Added an old-Python-inspired case matrix in
  `tests/tests/test_runs/test_3d_reps.rs`. The matrix imports all active
  GammaLoop-format fixtures once through `import graphs`, then probes:
  - the named old Python CFF/LTD/LTD-with-repeated-propagator scenarios,
  - uniform numerator sampling scale modes already represented by the current
    Rust structures,
  - all edge-linear numerator probes for every repeated old fixture,
  - all 210 one-loop box energy-degree caps with total degree at most 6, matching
    the old exhaustive one-loop CFF/LTD convergence sweep.
- Added the old slow five-loop ultimate-basis matrix as an ignored Rust test,
  matching the old `HYBRID3D_RUN_SLOW` behavior. It imports all five basis
  fixtures through the CLI and records CFF/LTD comparison status when explicitly
  run.
- The active expanded matrix currently records the following status:

```text
old Python 3Drep case matrix: 280 ok, 29 failed, 309 total
```

- Current known failing buckets from that matrix:
  - 8 cases fail with
    `this generalized CFF higher energy-numerator sector is not supported by the current Rust port`.
    These are high-power multiloop/repeated-channel sectors still left for the
    later generalized CFF completion.
  - 7 `cff_pureltd` repeated-propagator diagnostic comparisons fail, as expected:
    pure LTD without numerator derivatives is only a diagnostic branch and is
    not the correct repeated-propagator reference.
  - 14 cases build but numerically mismatch. These are currently:
    - one explicit quartic single-edge box probe,
    - the four-loop stress linear probe,
    - the pure one-loop box exhaustive bound entries with a single edge degree
      of 4, 5, or 6.
- Verified the expanded matrix with:

```text
cargo test -p gammaloop-integration-tests --test test_runs cli_old_python_case_matrix_records_current_three_drep_status -- --nocapture --test-threads=1
1 passed
```

- Verified the full 3D-reps integration-test selection with:

```text
cargo test -p gammaloop-integration-tests --test test_runs test_3d_reps -- --nocapture --test-threads=1
7 passed; 1 ignored
```

- Verified formatting and warning-as-error checks with:

```text
cargo fmt --check
passed

RUSTFLAGS=-Dwarnings cargo check -p gammaloop-api
passed

RUSTFLAGS=-Dwarnings cargo check -p gammaloop-integration-tests --test test_runs
passed
```

- Verified the selected GammaLoop gate with the expanded test selection:

```text
just test_gammaloop
1024 tests run: 1024 passed, 128 skipped
```

### Progress 2026-04-29 Removed obsolete 3Drep DOT parser path

- Removed the old Python-style DOT graph parser from `three-dimensional-reps`.
  The crate no longer exposes or uses `load_dot_graph`, `parse_dot_graph`,
  `DotGraph` as a `ThreeDGraphSource`, or `generate_3d_expression_from_dot_*`.
  3D expression generation now accepts `ParsedGraph` or a rich graph source
  supplied by GammaLoop.
- Removed the stale standalone `three-dimensional-reps-cli` references and the
  obsolete example DOT/script directory from the `three-dimensional-reps` crate.
  Diagnostic commands are documented as GammaLoop CLI commands under `3Drep`.
- Kept signature-expression parsing only in the `graph-from-signatures` path.
  Removed the old DOT-label momentum parser and old Hybrid/Vakint DOT round-trip
  tests. Remaining crate tests either use direct `ParsedGraph` fixtures for
  local unit coverage or go through the GammaLoop CLI.
- Updated `graph-from-signatures` DOT emission to produce GammaLoop-format edges
  with explicit scalar particles, global numerator `1`, and external edge names
  that preserve the physical external symbols. This allows CLI closure tests to
  project duplicate external edge instances back to the original signature
  symbols after GammaLoop import.
- Added a CLI closure test for a nontrivial 3-loop four-point Mercedes-style
  signature set:
  - runs `3Drep graph-from-signatures` with `--signatures-file` and
    `--num-vertices 4`,
  - imports the emitted DOT through the normal `import graphs` command,
  - extracts the parsed graph through GammaLoop's rich graph representation,
  - compares imported internal signatures against the original signature input,
    after grouping external edge instances by their physical names.
- Verified:

```text
cargo test -p gammaloop-integration-tests --test test_runs test_3d_reps -- --nocapture --test-threads=1
3 passed

cargo test -p three-dimensional-reps --features diagnostics,old_cff,test-support,eval,display --lib --tests -- --test-threads=1
40 passed

cargo fmt --check
passed

RUSTFLAGS=-Dwarnings cargo check -p gammaloop-api
passed

RUSTFLAGS=-Dwarnings cargo check -p gammaloop-integration-tests --test test_runs
passed

just test_gammaloop
1019 tests run: 1019 passed, 128 skipped
```

### Progress 2026-04-29 CLI imported-graph 3Drep tests

- Added GammaLoop-format DOT fixtures for the simple box and repeated-edge
  `box_pow3` case under `tests/resources/graphs/threedreps/`.
- Added CLI-session integration tests in
  `tests/tests/test_runs/test_3d_reps.rs` that:
  - start a normal test CLI state,
  - import `scalars-default.json`,
  - import the graph fixture through the regular `import graphs` command,
  - run `3Drep test-cff-ltd` on the imported graph,
  - read the generated 3Drep manifest and assert the same high-level CFF/LTD
    structures as the direct core tests.
- Initial targeted run showed that `3Drep` was still requiring generated
  integrands and therefore rejected graph lists imported through `import
  graphs`. I refactored `3Drep` graph selection so it can catalog either a
  generated `ProcessIntegrand` or the imported amplitude/cross-section graph
  lists already stored in the active state.
- The repeated-edge imported fixture also exposed the expected parser/model
  requirement for explicit vertex numerators on artificial two-point vertices;
  the new fixtures now use the same `node [num="1"]` convention as existing
  graph fixtures.
- The first CLI run also showed that the default iterative evaluator builder is
  not valid for LTD pseudo-orientations with undirected `x` entries: selecting
  one pseudo-orientation leaves theta factors from other pseudo-orientations
  unresolved. The `3Drep` diagnostic/evaluate commands now build with
  `iterative_orientation_optimization = false`, which preserves the orientation
  signs as evaluator parameters and supports the generalized oriented
  structures.
- The repeated-box CLI run then showed that numerator `1` was producing an
  empty automatic degree-bound list. That is insufficient for repeated
  propagators because the generator needs explicit zero caps for each loop edge
  to group repeated signatures with the intended edge powers. The command layer
  now initializes all internal paired propagator edges with degree `0` and
  overlays the numerator-analysis result on top.
- With the GammaLoop-format `box_pow3` fixture, the repeated channel is the
  signature class represented by imported edge `6` rather than the compact
  Python example's last edge. The CLI test now asserts the imported graph's
  deterministic structure directly: one LTD pseudo-orientation has the
  `-1/-3/-6` repeated-propagator variants on edge `6`, while the other
  pseudo-orientations remain single-variant.
- Verified the imported-graph command path and the full selected profile with:

```text
cargo test -p gammaloop-integration-tests --test test_runs cli_imported_box -- --nocapture --test-threads=1
2 passed

cargo test -p gammaloop-integration-tests --test test_runs test_3d_reps -- --nocapture --test-threads=1
4 passed

cargo fmt --check
passed

RUSTFLAGS=-Dwarnings cargo check -p gammaloop-api
passed

RUSTFLAGS=-Dwarnings cargo check -p gammaloop-integration-tests --test test_runs
passed

just test_gammaloop
1027 tests run: 1027 passed, 128 skipped
```

### Progress 2026-04-30 3Drep runtime settings and comparison verdicts

- Reworked `3Drep evaluate` and `3Drep test-cff-ltd` to use the active
  GammaLoop generation/runtime settings when selecting diagnostic evaluator
  behavior. The commands now report the requested runtime precision, selected
  evaluator method/backend, active numerator-sample normalization strategy,
  seed, input scale, and the definite value of the numerator interpolation
  scale `M`.
- Added the new user-facing generation option spelling
  `--numerator-samples-normalization` with values `never_M`, `M_for_all`, and
  `M_for_beyond_quadratic_only`. The old CLI spelling remains accepted as a
  compatibility alias for now, but the demo and completions use the new name.
- Fixed the generalized CFF construction path so `M_for_all` also routes the
  quadratic numerator-sampling nodes through the uniform `M` normalization,
  instead of reusing the special quadratic OSE-only path.
- Added precision-aware diagnostic evaluation for `Double`, `Quad`, and
  `ArbPrec`. `Double` can use the compiled backend selected by the global
  generation settings; `Quad` and `ArbPrec` force the eager evaluator. The
  command output now prints full relevant digits for the chosen precision.
- Extended diagnostic input generation with deterministic `--seed` and
  `--scale` controls, and made the printed parameter table include all evaluator
  inputs: external energies, external spatial momenta, loop spatial components,
  additional parameters, `m_uv`, `mu_r^2`, and `M`.
- Split timing output into `evaluator build time` and `sample evaluation time`,
  with dynamically scaled units. Reused cached evaluators leave the build-time
  entry empty.
- Changed `test-cff-ltd` so it only evaluates the numerator-normalization mode
  active in the global generation settings, reports CFF/LTD/PureLTD values with
  relative `Δ` against CFF, and ends with a colored `SUCCESS`/`FAIL` verdict.
  The CLI integration tests now assert the verdict result, which intentionally
  exposes the current representation gaps instead of merely checking that rows
  were produced.
- Added `--mass-shift` to control the starting pure-LTD split-mass diagnostic
  scale, defaulting to `1.0e-4 * --scale`, and display the per-shift values with
  three significant digits.
- Updated the demo run card to use the new normalization spelling and to show
  deterministic precision/seed/scale options for `evaluate` and `test-cff-ltd`.

Focused verification completed:

```text
cargo fmt
passed

cargo check -p gammaloop-api
passed

cargo check -p gammaloop-integration-tests --tests
passed

just clippy -- -D warnings
passed

cargo test -p three-dimensional-reps all_sampling_scale_mode_affects_quadratic_reconstruction --features diagnostics,old_cff,test-support,eval,display -- --nocapture
passed

cargo test -p gammaloop-api completion_offers_3drep_nested_selectors_and_enum_values -- --nocapture
passed

just build-cli
passed

./gammaloop --dev-optim ./examples/cli/scalar_topologies/three_d_expr_demo.toml run demo -c "quit -n"
passed; the comparison summary now correctly reports FAIL for the still-open
CFF/LTD numerical gaps
```

## 2026-04-30: make 3Drep comparisons verdict-bearing

- Added an explicit verdict object to the `3Drep test-cff-ltd` manifest. The
  command now compares direct finite CFF/LTD/pure-LTD evaluations against the
  CFF reference using a half-double-style threshold, rejects exactly identical
  formatted values across different direct representations, and applies a
  deliberately looser couple-digit threshold to the split-mass pure-LTD
  diagnostics.
- Added per-check records with absolute/relative differences, tolerances, and
  human-readable messages. The terminal summary now ends with a colored
  `SUCCESS` or `FAIL` status and includes a table of the individual checks.
- Existing cached comparison manifests without a verdict are no longer reused,
  so the new comparison outcome is always available after the next command run.
- Updated the 3Drep CLI integration-test harness so `test-cff-ltd` invocations
  assert the verdict outcome instead of only checking that evaluation rows exist.
  This intentionally exposes the current implementation gaps before we fix
  them.

Verification completed:

```text
cargo fmt
passed

cargo check -p gammaloop-api
passed

cargo check -p gammaloop-integration-tests --tests
passed
```

## 2026-04-30: 3Drep diagnostic evaluator NaN cleanup

- Fixed the `3Drep test-cff-ltd` diagnostic evaluator setup so it no longer
  samples the trivial zero kinematic point. The `ParamBuilder` now exposes
  grouped evaluator-input metadata and a deterministic nonzero diagnostic
  kinematic initializer used by `3Drep evaluate` and `test-cff-ltd`.
- Updated the printed and JSON evaluation records to keep reporting the input
  parameters with their source classification, including the automatically
  chosen diagnostic external energies, spatial momenta, loop momenta, and
  additional parameters.
- Fixed repeated-propagator `pure_ltd` diagnostics. The comparison harness now
  follows the old Python reference route: split repeated masses by several
  epsilon values, generate an ordinary LTD expression on each split-mass graph,
  and evaluate those split expressions instead of directly evaluating the
  degenerate repeated-mass pure-LTD expression.
- Added `mass_shift_values` to each evaluation record. Repeated-propagator
  pure-LTD rows now include the individual base and shifted masses for every
  split edge in JSON and in the table's mass-shift label.
- Added regression coverage that rejects `NaN` values in 3Drep comparison
  manifests and checks that repeated-propagator pure-LTD diagnostics report
  explicit split-mass values.

Verification completed:

```text
cargo fmt
passed

cargo check -p gammaloop-api
passed

just clippy -- -D warnings
passed

cargo test -p gammaloop-integration-tests --test test_runs 3drep -- --nocapture --test-threads=1
2 passed; 70 filtered out

cargo test -p gammaloop-integration-tests --test test_runs test_3d_reps:: -- --nocapture --test-threads=1
8 passed; 1 ignored; 63 filtered out
```

## 2026-04-30: richer `3Drep evaluate` and `test-cff-ltd` evaluation output

Follow-up implementation after the previous display/cache pass:

- `3Drep evaluate` now actually evaluates the cached oriented expression instead
  of only building the Symbolica evaluator artifacts. The command writes the
  resulting complex value, wall-clock evaluation time, evaluator-only timing,
  and the evaluator input parameter values into `evaluate_manifest.json` and
  displays the same information in colored tables.
- `3Drep test-cff-ltd` now evaluates every successfully built comparison case
  and records a per-case evaluation table with representation, sampling-scale
  mode, active `M` value when applicable, mass-shift label, status, complex
  value, wall-clock timing, and evaluator-only timing. When no interpolation
  scale is configured, the comparison still evaluates each method once and
  displays `M` as absent.
- Added small public `ParamBuilder`/`EvaluationMetaData` accessors so the CLI can
  reuse GammaLoop's existing `EvaluatorStack` machinery without duplicating the
  private evaluator input layout.
- Old cached `test-cff-ltd` manifests that do not contain evaluation rows are no
  longer reused; they are regenerated so the new summary table is populated.
- Extended the CLI integration tests to assert that `evaluate_manifest.json`
  contains an evaluated value, timings, and parameter records, and that every
  successful `test-cff-ltd` case contains at least one evaluation record with
  value and timing fields.

Verification run so far:

```text
cargo fmt
passed

cargo check -p gammaloop-api
passed

cargo test -p gammaloop-integration-tests --test test_runs 3drep -- --nocapture --test-threads=1
2 passed

cargo test -p gammaloop-integration-tests --test test_runs cli_validate_build_and_evaluate_use_gammaloop_graph_state -- --nocapture --test-threads=1
1 passed

just clippy -- -D warnings
passed

just build-cli
passed

./gammaloop --dev-optim ./examples/cli/scalar_topologies/three_d_expr_demo.toml run demo -c "quit -n"
passed; verified the new evaluate value/timing/parameter table and the
test-cff-ltd per-method/per-M timing table in the normal CLI path

just test_gammaloop
1030 tests run: 1030 passed, 128 skipped
```

## 2026-04-30: second display/selector pass for `3Drep`

- Fixed denominator-tree rendering in the overview and detail views. The display
  now reports node IDs and child node IDs with the form
  `<node_id>: <surface-or-1> -> [child_node_ids]`, while surfaces remain
  attached to the node data. The table truncation helper is now ANSI-aware, so
  colored `den tree` cells no longer leave unterminated color spans when
  shortened.
- Reworked orientation display metadata. Orientation labels are sorted by their
  base tag and then by numerator/energy-map index, display as `tag|N0`,
  `tag|N1`, ..., and detail selectors now support all of:
  `tag`, `tag|N0`, and `tag|N0|variant_id`. A final global label pass is
  applied after remapping and variant fusion so CLI-generated expressions do
  not inherit duplicated local `N0` indices from lower-sector components.
- External half-edges in orientation tags are rendered in gray independently of
  their `x`, `+`, `-`, or `0` assignment; internal entries keep their semantic
  colors.
- Changed variant surface lists to report all distinct multiplicative powers
  encountered on root-to-leaf denominator branches. For example, if a surface
  appears once on one branch and twice on another it displays as `5,5^2`;
  numerator-only occurrences remain parenthesized.
- Changed top-level energy-degree bounds so edge IDs are green, kept high
  powers red, and changed `3Drep build` numerator display to use the existing
  GammaLoop/Symbolica pretty printer instead of canonical strings.
- Corrected `NumeratorSamplingScaleMode::All` so it activates for every
  positive numerator energy degree. Verified that remaining factors such as
  `4^2` in `--numerator-sampling-scale-mode all` come from physical/denominator
  half-edge factors in the known-factor and lower-sector decomposition, not from
  numerator sampling-node normalization.
- Updated `examples/cli/scalar_topologies/three_d_expr_demo.toml` to showcase
  the new `tag|N0|variant` detail selector.

Verification completed:

```text
cargo fmt
passed

cargo test -p three-dimensional-reps --features display,diagnostics,test-support,eval,old_cff display::tests -- --nocapture
3 passed

cargo check -p gammaloop-api --features default
passed

just clippy -- -D warnings
passed

just build-cli
passed

./gammaloop --dev-optim -s /tmp/gammaloop-3drep-display-state ./examples/cli/scalar_topologies/three_d_expr_demo.toml run demo -c "quit -n"
passed; verified cached detail selector `xxxx-++xxx|N0|3`, node-id denominator trees, pretty numerator output, green energy-bound edge IDs, and cache reuse output

just test_gammaloop
1030 tests run: 1030 passed, 128 skipped
```

### Progress 2026-04-30 graph-from-signatures closure hardening

- Added a CLI-driven closure harness for `3Drep graph-from-signatures` in
  `tests/tests/test_runs/test_3d_reps.rs`. Each case specifies only a list of
  `(momentum_signature, mass)` factors; the helper derives the `prop(q,m)`
  input expression, runs `3Drep graph-from-signatures`, asserts the generated
  DOT contains no explicit `lmb_rep`, reloads that DOT via
  `import graphs --inline-dot ...` with the scalars model, and compares the imported
  internal signatures and masses against the input multiset up to momentum
  reversal.
- The new closure coverage includes:
  - a one-loop repeated-signature box-like chain,
  - mixed `p`/`q` external prefixes,
  - a two-loop sunrise,
  - a three-loop Mercedes-style case,
  - four-loop and five-loop banana cases.
- `graph-from-signatures` no longer emits explicit `lmb_rep` attributes in
  GammaLoop DOT output. It emits only topology, masses, external leg names, and
  `lmb_id` on loop-basis carrier edges, so normal GammaLoop import must recover
  the momentum signatures from the graph itself.
- The reconstruction search now tracks external-momentum imbalance while
  placing internal edges and, for CLI output, minimizes the number of external
  half-edges before rendering. This avoids accepting the first loop-balanced
  topology when it creates unnecessary external attachments.
- Added an explicit connected-graph consistency check:
  `L = E - V + 1` for the explicit-edge representation. The command
  `3Drep graph-from-signatures --signatures
  'prop(k1,m1)*prop(k1,m11)*prop(k1+p1,m2)*prop(k1+p1+p2,m3)*prop(k1+p1+p2+p3,m4)'
  --num-vertices 4` now errors cleanly, because five explicit propagator
  signatures with one loop momentum require five internal vertices. With
  `--num-vertices 5`, the command emits the expected one-loop explicit repeated
  edge chain with no `lmb_rep`.

Verification completed:

```text
cargo fmt
passed

just clippy -- -D warnings
passed

cargo test -p three-dimensional-reps --features diagnostics,old_cff,test-support,eval,display graph_signatures -- --nocapture --test-threads=1
5 passed; 35 filtered out

cargo test -p gammaloop-integration-tests --test test_runs graph_from_signatures -- --nocapture --test-threads=1
2 passed; 70 filtered out

cargo check -p gammaloop-api
passed

cargo check -p gammaloop-integration-tests --test test_runs
passed

just build-cli
passed

just test_gammaloop
1027 tests run: 1027 passed, 128 skipped
```

### Progress 2026-04-29 CLI-only box and box_pow3 coverage

- Removed the older direct core tests for the simple box and repeated
  `box_pow3` examples from `tests/tests/test_runs/test_3d_reps.rs`.
- The two remaining tests now exercise only the intended GammaLoop command path:
  start a CLI state, `import model`, import the DOT fixture through the regular
  `import graphs` command, run `3Drep test-cff-ltd`, read the generated manifest,
  and load the saved oriented expressions from the 3Drep workspace.
- Added a test-side rich run object carrying the CLI state, parsed graph,
  generated manifest, source edge ids, and the CFF/LTD/pure-LTD expressions in
  both source-edge and compact-edge numbering. The tests use this to keep the
  same structural and numerical assertions that the direct tests had before.
- Updated the GammaLoop-format `box.dot` and `box_pow3.dot` fixtures so their
  imported topology reproduces the Python reference signature chain
  `k, k+p1, k+p1+p2, k+p1+p2+p3`. The tests now assert these parsed internal
  signatures explicitly, because imported GammaLoop graphs derive their loop
  momentum basis from topology and `lmb_id`; the printed `lmb_rep` attribute is
  not the import-time source of truth.
- The CLI tests cover:
  - manifest metadata and automatic zero energy-power caps,
  - CFF/LTD/pure-LTD orientation and surface counts,
  - CFF orientation labels, prefactors, half-edge prefactors, and denominator
    surface ids for the simple box,
  - repeated-propagator LTD variants, prefactors, repeated half-edge factors,
    and denominator surfaces for `box_pow3`,
  - pure-LTD mass-shift split values and monotonic convergence for `box_pow3`,
  - numerical agreement with the deterministic Python references for CFF and
    LTD.
- Verified the focused CLI-only suite with:

```text
cargo test -p gammaloop-integration-tests --test test_runs test_3d_reps -- --nocapture --test-threads=1
2 passed
```

- Verified the formatting, warning-as-error checks, and selected full profile
  after the CLI-only migration with:

```text
cargo fmt --check
passed

RUSTFLAGS=-Dwarnings cargo check -p gammaloop-api
passed

RUSTFLAGS=-Dwarnings cargo check -p gammaloop-integration-tests --test test_runs
passed

just test_gammaloop
1025 tests run: 1025 passed, 128 skipped
```

### Progress 2026-04-30 Remaining old-Python case matrix fixes

- Resumed the old-Python parity matrix after the first CLI-only coverage pass.
  The starting diagnostic matrix had 29 non-ok rows:
  - 8 unsupported generalized CFF high-power sectors,
  - 14 one-loop high-contact numerical mismatches,
  - 7 pure-LTD repeated-propagator diagnostic mismatches.
- Re-read the generalized LTD notes around high-contact one-loop CFF sectors,
  especially the terminal-contact completion and lower-sector normal-form
  recursion. The failing degree-4/5/6 single-edge box sectors matched the
  documented case that should be handled by the lower-sector recursive builder
  rather than by the quadratic E-surface-only shortcut.
- Updated `BoundedCffBuilder` so isolated high-power one-loop sectors no longer
  enter the quadratic E-surface-only path. That path is now restricted to
  non-repeated one-loop sectors with all edge bounds at most quadratic; higher
  one-loop sectors fall through to `KnownFactorCffBuilder`.
- Verified the targeted diagnostic matrix after the routing fix:

```text
cargo fmt && cargo test -p gammaloop-integration-tests --test test_runs cli_old_python_case_matrix_records_current_three_drep_status -- --nocapture --test-threads=1
old Python 3Drep case matrix: 293 ok, 16 failed, 309 total
```

- The remaining non-ok rows are now:
  - one real CFF/LTD numerical sign mismatch in `four_loop_stress_linear`,
  - 8 unsupported multiloop/repeated generalized CFF high-power sectors,
  - 7 expected pure-LTD repeated-propagator diagnostic mismatches, where the
    derivative-free pure-LTD diagnostic is not expected to agree for repeated
    propagators.

### Progress 2026-04-30 Old-Python case matrix closure

- Generalized the high-power CFF dispatch further. `BoundedCffBuilder` now sends
  every sector with energy degree greater than one, or repeated logical-channel
  total degree greater than two, through the lower-sector known-factor recursive
  builder before reporting `CffHigherEnergyPowerNotImplemented`. The quadratic
  E-surface-only builder remains restricted to non-repeated one-loop sectors with
  all bounds at most quadratic.
- This removed the unsupported generalized CFF rows from the old-Python case
  matrix. The interim diagnostic matrix was:

```text
cargo test -p gammaloop-integration-tests --test test_runs cli_old_python_case_matrix_records_current_three_drep_status -- --nocapture --test-threads=1
old Python 3Drep case matrix: 300 ok, 9 failed, 309 total
```

- Traced the remaining `four_loop_stress_linear` and
  `four_loop_stress_quadratic` sign/magnitude mismatches to the GammaLoop-format
  fixture, not to expression generation. The external `p3` leg was attached to
  `K0`, which made the imported graph expose only one repeated group of size 3.
  The old Python topology has repeated groups of sizes 2 and 3. Moving `p3` to
  the central `R` vertex restores the intended signatures and repeated-channel
  structure.
- Updated the import-closure assertion for `four_loop_stress.dot` to expect
  repeated groups `[2, 3]`, and verified the fixture through the standard
  GammaLoop `import graphs` path:

```text
cargo test -p gammaloop-integration-tests --test test_runs cli_imports_old_python_equivalent_threedrep_fixtures -- --nocapture --test-threads=1
passed
```

- The only remaining matrix failures were the seven repeated-propagator
  pure-LTD rows. These are not valid numerical equalities in the derivative-free
  implementation: pure LTD is kept as a diagnostic/mass-shift comparison mode,
  and without numerator derivatives it is not expected to agree for repeated
  propagators. The matrix now treats those rows as build-only diagnostics while
  retaining CFF/LTD numerical comparisons for the derivative-free repeated
  representation.
- Verified the complete old-Python-inspired case matrix after the routing,
  fixture, and diagnostic-harness updates:

```text
cargo test -p gammaloop-integration-tests --test test_runs cli_old_python_case_matrix_records_current_three_drep_status -- --nocapture --test-threads=1
old Python 3Drep case matrix: 309 ok, 0 failed, 309 total
```

- Updated the stale core unit test that still expected
  `CffHigherEnergyPowerNotImplemented` for the `sunrise_pow4` high-power CFF
  sector. That sector is now supported through the known-factor recursive
  builder, so the test now asserts the generated expression shape and the
  presence of numerator-surface known-factor completion terms.
- Re-ran the focused 3Drep integration suite and the core crate tests:

```text
cargo test -p gammaloop-integration-tests --test test_runs test_3d_reps -- --nocapture --test-threads=1
7 passed; 1 ignored

cargo test -p three-dimensional-reps --features diagnostics,old_cff,test-support,eval,display --lib --tests -- --test-threads=1
40 passed
```

- Re-ran formatting, warning-as-error checks, and the full selected GammaLoop
  gate:

```text
cargo fmt --check
passed

RUSTFLAGS=-Dwarnings cargo check -p gammaloop-api
passed

RUSTFLAGS=-Dwarnings cargo check -p gammaloop-integration-tests --test test_runs
passed

just test_gammaloop
1024 tests run: 1024 passed, 128 skipped
```

### Progress 2026-04-30 CLI coverage audit after matrix closure

- Audited the old Python pytest surface against the Rust `test_3d_reps` coverage
  after the old-Python case matrix reached `309 ok, 0 failed`. The generation
  and eager-evaluation parity lanes are now covered by:
  - GammaLoop CLI DOT import and graph-shape closure for all old Python example
    topologies ported to GammaLoop DOT syntax,
  - graph-from-signatures closure through the real `import graphs` command,
  - automatic EMR/LMB numerator energy-power cap extraction,
  - the old-Python-inspired CFF/LTD case matrix, including high-contact
    one-loop CFF, multiloop/repeated high-power CFF, uniform numerator sampling
    scale modes, and repeated-propagator LTD diagnostics,
  - explicit box and `box_pow3` oriented-structure checks through
    `3Drep test-cff-ltd`,
  - the slow five-loop ultimate-basis parity probe.
- Ran the slow probe explicitly even though it stays ignored in the default
  suite, matching the old Python `HYBRID3D_RUN_SLOW` convention:

```text
cargo test -p gammaloop-integration-tests --test test_runs cli_old_python_slow_five_loop_ultimate_basis_case_matrix -- --ignored --nocapture --test-threads=1
1 passed
```

- Added a CLI integration test for the diagnostic path not covered by
  `test-cff-ltd`: `3Drep validate`, `3Drep build`, and `3Drep evaluate` now run
  on an imported box graph through the GammaLoop state and assert the generated
  validation JSON, oriented-expression artifact, Symbolica expression artifact,
  and parameter-builder artifact.
- Updated `crates/three-dimensional-reps/README.md` to remove the stale
  high-power CFF limit statement. The README now describes the validated
  old-Python generation/eager-evaluation surface and keeps the important
  distinction that compiled Symbolica runtime evaluation is CLI-side evaluator
  construction, not part of the library crate's eager evaluator.
- Re-ran focused and full verification:

```text
cargo test -p gammaloop-integration-tests --test test_runs test_3d_reps -- --nocapture --test-threads=1
8 passed; 1 ignored

cargo test -p three-dimensional-reps --features diagnostics,old_cff,test-support,eval,display --lib --tests -- --test-threads=1
40 passed

cargo fmt --check
passed

RUSTFLAGS=-Dwarnings cargo check -p gammaloop-api
passed

RUSTFLAGS=-Dwarnings cargo check -p gammaloop-integration-tests --test test_runs
passed

just test_gammaloop
1025 tests run: 1025 passed, 128 skipped
```

### Progress 2026-04-30 inline graph import and 3D expression demo card

- Added an inline DOT import mode to the existing GammaLoop graph importer:
  `import graphs --inline-dot '<dot graph content>'`. The old path-based
  `import graphs <path>` mode is still the default, and inline imports are
  routed through the same `Graph::from_string(...)` GammaLoop parsing path used
  by normal graph ingestion.
- Added `3d_expr` and `3d-expr` aliases for the existing `3Drep` command so
  examples can use the terminology of the oriented 3D-expression workflow
  without changing the underlying command implementation.
- Updated the interactive REPL completion machinery and tests so tab completion
  covers:
  - `import graphs --inline-dot`,
  - the `3d_expr` command alias,
  - nested `3d_expr` subcommands,
  - process/integrand selectors and representation/scale-mode enum values.
- Added `examples/cli/scalar_topologies/three_d_expr_demo.toml`, a compact
  scalar-topology run card that imports a `box_pow3`-style graph from an inline
  DOT string with an energy-dependent numerator, then demonstrates
  `3d_expr validate`, `build`, `evaluate`, `test-cff-ltd`, and
  `graph-from-signatures`.
- Updated the existing `cli_validate_build_and_evaluate_use_gammaloop_graph_state`
  integration test to exercise the new inline import path before running the
  3D-expression validate/build/evaluate workflow.

Focused verification completed so far:

```text
cargo fmt --check
passed

RUSTFLAGS=-Dwarnings cargo check -p gammaloop-api
passed

cargo test -p gammaloop-api completion_offers -- --nocapture
52 passed; 244 filtered out

cargo run -p gammaloop-api --bin gammaloop -- --clean-state examples/cli/scalar_topologies/three_d_expr_demo.toml
passed

cargo test -p gammaloop-integration-tests --test test_runs test_3d_reps -- --nocapture --test-threads=1
8 passed; 1 ignored

RUSTFLAGS=-Dwarnings cargo check -p gammaloop-integration-tests --test test_runs
passed

just test_gammaloop
1027 tests run: 1027 passed (1 leaky), 128 skipped
```

### Progress 2026-04-30 3Drep display, caching, and demo polish

- Removed the temporary `3d_expr`, `3d-expr`, `3drep`, and `threedreps`
  aliases from the top-level command registration. The canonical user command
  is now only `3Drep`; completion tests were updated accordingly while keeping
  `import graphs --inline-dot` completion coverage.
- Updated `examples/cli/scalar_topologies/three_d_expr_demo.toml` so the inline
  DOT graph no longer specifies `lmb_rep`. The demo controls the loop-momentum
  basis only through `lmb_id`, uses only `3Drep` commands, demonstrates cache
  reuse, and shows a concrete orientation/variant selector using
  `xxxx-++xxx|3`.
- Reworked the `three-dimensional-reps` display layer:
  - surface IDs now render as `S0`, `S1`, ... for generalized linear surfaces;
  - surface rows merge class/origin/numerator-only into a single colored `type`
    column with helper surfaces carrying the `_sp` suffix;
  - linear surface and energy-map expressions are sorted, explicitly signed,
    aligned between OSE and shift groups, and color-coded;
  - variant rows now show colored origin, rational prefactor, combined
    M/half-edge factors such as `[M^2,4,7^3]`, sorted surface IDs with powers,
    per-variant node counts, and colored denominator trees using `→`;
  - the global summary includes total node count;
  - detail mode now suppresses the full structure and supports
    `orientation|variant_id`, filtering out matching orientations that do not
    contain the requested variant.
- Added variant fusion in the 3D expression representation for variants with
  identical rational prefactor, half-edge factors, uniform-scale power, and
  numerator-surface set. Fused variants are represented by a new root node whose
  children are the original denominator trees, and mixed origins are labelled
  `mixed`.
- Removed the custom `RationalCoefficient` type from the implementation. Variant
  prefactors and linear-expression coefficients now use Symbolica `Atom`
  coefficients directly, with serde support through canonical strings in JSON
  and the existing Symbolica-state-aware bincode derivations.
- Added per-process/integrand/graph 3Drep workspace subfolders. `build`,
  `evaluate`, and `test-cff-ltd` now reuse matching cached artifacts by default,
  emit the exact relative path being reused, and honor `--clean` to force
  recomputation. Metadata mismatches fall back to recomputation and overwrite
  the JSON when saving is enabled.
- Added colored tabled summaries for `3Drep evaluate` and
  `3Drep test-cff-ltd`, while still saving the JSON manifests for inspection.
- During verification, the first dev-optim rebuild hit a local
  `No space left on device` error while writing generated incremental compiler
  artifacts. Only generated `target/*/incremental` directories were removed,
  after which the dev-optim build completed.

Verification completed:

```text
cargo fmt
passed

just clippy -- -D warnings
passed

cargo check -p three-dimensional-reps --features diagnostics,display,eval,test-support,old_cff
passed

cargo check -p gammaloop-api
passed

cargo check -p gammaloop-integration-tests --test test_runs
passed

cargo test -p three-dimensional-reps --features diagnostics,old_cff,test-support,eval,display --lib --tests -- --test-threads=1
40 passed

cargo test -p gammaloop-api completion_offers -- --nocapture
52 passed; 244 filtered out

cargo test -p gammaloop-integration-tests --test test_runs 3drep -- --nocapture --test-threads=1
2 passed; 70 filtered out

cargo test -p gammaloop-integration-tests --test test_runs graph_from_signatures -- --nocapture --test-threads=1
2 passed; 70 filtered out

just build-cli
passed

./gammaloop --dev-optim ./examples/cli/scalar_topologies/three_d_expr_demo.toml run demo -c "quit -n"
passed; verified cache reuse, detail-only selector output, evaluate summary, and comparison summary

just test_gammaloop
1027 tests run: 1027 passed, 128 skipped
```

### Progress 2026-04-30 post-clippy cleanup for runtime-output pass

- Addressed the clippy cleanup from the runtime-output pass by grouping
  evaluator-build inputs, grouping comparison-distance data, using direct vector
  initialization, and simplifying the status-color branch.
- Removed the generated demo state/workspace artifacts after verification. The
  only remaining untracked demo-adjacent file is the pre-existing local
  `three_d_expr_demo_vh_do_not_push.toml`, which remains untouched.
- The comparison harness now asserts the `test-cff-ltd` verdict, so the broader
  3Drep integration tests are expected to fail until the known representation
  discrepancies are fixed in the next implementation pass.

Verification completed:

```text
cargo fmt
passed

cargo check -p gammaloop-api
passed

cargo check -p gammaloop-integration-tests --tests
passed

just clippy -- -D warnings
passed

cargo test -p three-dimensional-reps all_sampling_scale_mode_affects_quadratic_reconstruction --features diagnostics,old_cff,test-support,eval,display -- --nocapture
1 passed

cargo test -p gammaloop-api completion_offers_3drep_nested_selectors_and_enum_values -- --nocapture
1 passed
```

### Progress 2026-04-30 multi-way comparison discrepancy pass

- Investigated the first real `test-cff-ltd` discrepancy on the one-loop box
  with unit numerator and no repeated propagators. The LTD and PureLTD values
  were exactly four times the CFF value because the diagnostic CLI was building
  the full orientation-theta-gated parametric atom and then evaluating it once
  per partial LTD orientation. In the existing theta convention, undirected
  entries evaluate as zero and therefore pass the `IF(theta + 1)` gate, so each
  partial orientation activated the full LTD sum.
- Added a diagnostic 3D-expression atom path that sums orientation
  contributions without orientation theta gates and uses the orientation energy
  substitution maps directly. The normal GammaLoop parametric atom path remains
  theta-gated for the production orientation machinery.
- Updated `3Drep evaluate` and `3Drep test-cff-ltd` to evaluate that diagnostic
  theta-free atom once. The simple box now compares CFF, LTD and PureLTD at the
  expected numerical precision instead of showing the factor-four overcount.
- Reworked the PureLTD split-mass verdict for repeated propagators. The
  diagnostic now reports all split-mass trials but accepts the best finite mass
  shift within the loose mass-shift threshold, matching the intended
  couple-digit role of the mass-shift comparison rather than requiring
  monotonic convergence of the default split sequence.
- Ran the old Python equivalent case matrix through the Rust harness and
  obtained `309 ok, 0 failed, 309 total`; the deliberately slow five-loop
  ultimate-basis parity probe remains implemented and ignored, matching the old
  Python slow-test treatment.
- Checked the automatic energy-power cap test with a high-power LMB numerator.
  The old Python reference also shows a CFF/LTD mismatch for
  `loops[0][0]**3`, while the corresponding EMR case `edges[0][0]**3` agrees.
  This was initially treated as a build-only cap-extraction parity gap, but the
  later carrier-edge substitution investigation below identified and fixed the
  Rust-side cause.

Verification completed:

```text
cargo fmt
passed

cargo test -p gammaloop-integration-tests --test test_runs cli_imported_box_3drep_test_uses_gammaloop_graph_path -- --nocapture --test-threads=1
passed

cargo test -p gammaloop-integration-tests --test test_runs cli_imported_box_pow3_3drep_test_uses_gammaloop_graph_path -- --nocapture --test-threads=1
passed

cargo test -p gammaloop-integration-tests --test test_runs cli_reconstructs_energy_power_caps_from_imported_graph_numerators -- --nocapture --test-threads=1
passed

cargo test -p gammaloop-integration-tests --test test_runs test_3d_reps:: -- --nocapture --test-threads=1
8 passed, 0 failed, 1 ignored; old Python 3Drep case matrix reported 309 ok, 0 failed, 309 total

cargo check -p gammaloop-api
passed

cargo check -p gammaloop-integration-tests --tests
passed

just clippy -- -D warnings
passed

just test_gammaloop
1031 tests run: 1031 passed, 128 skipped

python3 IGNORE/generalised_ltd/generalised_ltd.py test --dot IGNORE/generalised_ltd/examples/graphs/box.dot --energy-degree-bounds 0:3,1:0,2:0,3:0 --numerator-expr 'loops[0][0]**3' --uniform-numerator-sampling-scale none
confirmed old Python CFF/LTD mismatch for the analogous LMB high-power numerator

python3 IGNORE/generalised_ltd/generalised_ltd.py test --dot IGNORE/generalised_ltd/examples/graphs/box.dot --energy-degree-bounds 0:3,1:0,2:0,3:0 --numerator-expr 'edges[0][0]**3' --uniform-numerator-sampling-scale none
confirmed old Python CFF/LTD agreement for the EMR high-power numerator
```

### Progress 2026-04-30 LMB carrier-energy substitution fix

- Investigated the remaining conceptual discrepancy between a numerator written
  as `K(0,spenso::cind(0))^3` and the same numerator written as
  `Q(e_lmb,spenso::cind(0))^3` for the carrier edge of the first loop momentum
  basis vector.
- Confirmed that the automatic energy-degree analyzer already assigns the
  `K(i,0)` power to the corresponding LMB carrier edge. The failure was in the
  diagnostic numerator substitution: `K(i,0)` was using the generated
  `loop_energy_map`, while generalized CFF numerator sampling is fundamentally
  expressed by the `edge_energy_map`. Some finite-difference/contact sample maps
  are not globally integrable into a single loop-energy assignment, so
  `loop_energy_map[i]` can legitimately differ from `edge_energy_map[e_lmb]`.
- Fixed the GammaLoop numerator substitution so `K(i,...)` is lowered through
  the carrier edge: `K(i,0)` now uses exactly the same energy map as
  `Q(loop_edges[i],0)`, and `K(i,mink(...))` uses the carrier edge spatial vector
  plus that same mapped energy component.
- Updated the standalone `three-dimensional-reps` evaluator similarly for
  `loops[...]` numerator expressions: if the parsed graph identifies a loop
  carrier edge by name, loop four-vectors are assembled from that carrier edge's
  energy map and spatial momentum.
- Re-enabled the LMB high-power numerator parity check in the CLI integration
  cap-extraction test and added an explicit generated matrix subcase for
  `loops[0][0]**3`.

Verification completed:

```text
cargo fmt
passed

cargo check -p gammaloop-api
passed

cargo test -p gammaloop-integration-tests --test test_runs cli_reconstructs_energy_power_caps_from_imported_graph_numerators -- --nocapture --test-threads=1
passed; the K(0,0)^3 box cap/parity case now reports SUCCESS

cargo test -p gammaloop-integration-tests --test test_runs cli_old_python_case_matrix_records_current_three_drep_status -- --nocapture --test-threads=1
1 passed; generated matrix reported 310 ok, 0 failed, 310 total

cargo test -p gammaloop-integration-tests --test test_runs test_3d_reps:: -- --nocapture --test-threads=1
8 passed, 0 failed, 1 ignored; generated matrix reported 310 ok, 0 failed, 310 total

cargo check -p gammaloop-api
passed

cargo check -p gammaloop-integration-tests --tests
passed

just clippy -- -D warnings
passed

just test_gammaloop
1031 tests run: 1031 passed, 128 skipped
```

### Progress 2026-04-30 `3Drep test-cff-ltd` display and mass-shift audit

- Started a follow-up pass on the `test-cff-ltd` report formatting and pure-LTD
  split-mass diagnostics after the demo run showed the mass-shift row failing
  while CFF and LTD agreed.
- Planned changes for this pass:
  - display final-check distances with three significant digits,
  - assign stable IDs to each evaluation row and reference those IDs from the
    final check table,
  - change the default split-mass ladder from `1.0e-4 * scale` with half-steps
    to `1.0e-2 * scale` followed by three decade-smaller shifts,
  - rerun the scalar-topology demo with cache invalidation and inspect whether
    the apparent best-at-largest-epsilon behavior is numerical instability or a
    construction/evaluation bug.

Implemented and investigated:

- Added stable sequential IDs to `ThreeDrepEvaluationRecord`, rejected cached
  comparison manifests that do not contain ordered IDs, added the `id` column to
  the evaluation table, and included the same ID in final-check labels. Final
  direct and mass-shift checks now explicitly refer to the compared evaluation
  and the CFF reference, e.g. `#1 ... vs #0`.
- Shortened final-check `abs diff`, `rel diff`, and `tolerance` fields to three
  significant digits.
- Changed the default pure-LTD split-mass starting shift to `1.0e-2 * scale`
  and the ladder to `[eps, eps/10, eps/100, eps/1000]`.
- Removed the explicit too-small `--mass-shift 2.0e-4` from the demo run card
  and added `--clean` to the demo comparison command so the example cannot
  silently reuse stale comparison manifests.
- Reproduced the old demo failure with the explicit small starting shift
  `2.0e-4` in Double precision: the best row had relative difference
  `1.64e-2`, just outside the loose `1.0e-2` mass-shift threshold, and smaller
  shifts rapidly lost precision. With the new default `2.0e-2` starting shift,
  the Double diagnostic succeeds with best relative difference `2.54e-5`.
- Repeated the same split-mass sequence in Quad precision. The relative
  difference improves monotonically as the shift decreases down to `2.0e-5`,
  reaching `2.71e-11`. This points to cancellation/conditioning of the
  split-mass diagnostic in Double precision, not to a wrong split-mass
  construction.

Verification completed:

```text
cargo fmt
passed

cargo check -p gammaloop-api
passed

cargo check -p gammaloop-integration-tests --tests
passed

just build-cli
passed

./gammaloop --dev-optim ./examples/cli/scalar_topologies/three_d_expr_demo.toml run demo -c "quit -n"
passed after cleaning the previously generated demo state/workspace; comparison status SUCCESS

./gammaloop --dev-optim -n ./examples/cli/scalar_topologies/three_d_expr_demo.toml 3Drep test-cff-ltd ... --precision Double --mass-shift 2.0e-4 --clean
reproduced the old small-shift Double failure: best mass-shift rel diff 1.64e-2

./gammaloop --dev-optim -n ./examples/cli/scalar_topologies/three_d_expr_demo.toml 3Drep test-cff-ltd ... --precision Quad --mass-shift 2.0e-2 --clean
confirmed smaller shifts converge in Quad: best mass-shift rel diff 2.71e-11

cargo test -p gammaloop-integration-tests --test test_runs test_3d_reps::cli_imported_box_pow3_3drep_test_uses_gammaloop_graph_path -- --nocapture --test-threads=1
passed

cargo test -p gammaloop-integration-tests --test test_runs test_3d_reps:: -- --nocapture --test-threads=1
8 passed, 0 failed, 1 ignored; generated matrix reported 310 ok, 0 failed, 310 total

just clippy -- -D warnings
passed

post-label cleanup:
cargo fmt
cargo check -p gammaloop-api
just clippy -- -D warnings
passed

just test_gammaloop
1031 tests run: 1031 passed, 128 skipped
```

## 2026-05-01: 3Drep display polish and stronger CLI comparison coverage

Follow-up request:

- shorten the `den tree` preview in `3Drep build` to 30 visible characters and
  keep the truncation ellipsis in normal terminal color,
- shorten verbose origin labels with UTF-8 beta/gamma and compact LTD
  abbreviations,
- spell out the `orientation` table header,
- rename the final comparison table's `check` column to `evaluation id` and
  show only `#<id>`,
- change the default pure-LTD split-mass start to `1.0 * --scale`,
- add CLI tests for the most demanding repeated-propagator/high-energy-power
  case across all numerator-sample normalization modes and across eager,
  symjit, and assembly double evaluators,
- extend the scalar-topology demo card with Double, Quad, ArbPrec, symjit, and
  assembly examples for `evaluate` and `test-cff-ltd`.

Implementation notes:

- Updated the colored 3D expression display code so the denominator tree is
  truncated at 30 visible characters and resets ANSI styling before appending
  `...`.
- Added compact origin rendering (`ltd_confluent` to `ltd_cflt`,
  `bounded_degree` to `bd`, and `beta`/`gamma` to `β`/`γ`) while keeping the
  raw origin data unchanged in the serialized expression.
- Added test helpers that run `3Drep test-cff-ltd` through a real CLI session,
  read the generated manifest, assert `SUCCESS`, check stable evaluation IDs,
  parse direct CFF/LTD values, and compare complex values numerically.
- Added one normalization-strategy test for `box_pow3` with a high-power
  repeated-propagator numerator. It runs `never_M`,
  `M_for_beyond_quadratic_only`, and `M_for_all`, confirms the global setting
  reaches the command, verifies the default mass-shift start is exactly
  `--scale`, and compares direct CFF/LTD values as well as CFF values across
  strategies.
- Added one compiled-backend test for the same case. It compares eager, symjit,
  and assembly Double evaluations and asserts the compiled backends recorded in
  the manifest match the requested backend.
- Expanded `examples/cli/scalar_topologies/three_d_expr_demo.toml` to run the
  same oriented-expression evaluation/comparison in Double, Quad, ArbPrec,
  symjit Double, and assembly Double modes. The demo now uses `--scale 0.25`
  so the new default `mass_shift_start = scale` remains a clean split-mass
  diagnostic without producing a negative shifted mass in the first row.

Verification completed:

```text
cargo fmt
passed

cargo check -p gammaloop-api
passed

cargo check -p gammaloop-integration-tests --tests
passed

cargo test -p gammaloop-integration-tests --test test_runs \
  test_3d_reps::cli_box_pow3_high_power_agrees_for_all_numerator_sample_normalizations \
  -- --nocapture --test-threads=1
passed

cargo test -p gammaloop-integration-tests --test test_runs \
  test_3d_reps::cli_box_pow3_high_power_compiled_double_backends_agree_with_eager \
  -- --nocapture --test-threads=1
passed

cargo test -p gammaloop-integration-tests --test test_runs test_3d_reps:: \
  -- --nocapture --test-threads=1
10 passed, 0 failed, 1 ignored; generated matrix reported 310 ok, 0 failed, 310 total

just build-cli
passed

./gammaloop --dev-optim ./examples/cli/scalar_topologies/three_d_expr_demo.toml run demo -c "quit -n"
passed

just clippy -- -D warnings
passed

just test_gammaloop
1033 tests run: 1033 passed, 128 skipped
```

## 2026-05-01: Dense numerator probes and high-precision deltas

Follow-up request:

- confirm that multi-way tests use numerators matching the intended energy
  powers,
- make those numerators dense four-dimensional contractions rather than pure
  temporal-component monomials,
- compute `3Drep test-cff-ltd` deltas/differences without downcasting the
  compared high-precision values to `f64`,
- commit and push the accumulated changes.

Implementation notes:

- Reworked the imported-graph CLI test numerators used for automatic energy
  power extraction. They now build explicit four-component Minkowski
  contractions such as `Q(e,0) Q(p,0) - Q(e,1) Q(p,1) - ...`, so the temporal
  component still drives the requested cap while the numerator depends on all
  four components.
- Added the analogous dense LMB numerator construction with
  `K(loop,component)` contracted against an external-edge `Q`; this keeps the
  LMB cap extraction path covered while avoiding evaluator-internal symbolic
  indices that are not parameterized in the standalone diagnostic evaluator.
- Updated the local old-Python parity harness so power-bound probe cases marked
  by their bounds are evaluated with dense `dot(edges[i], ext[j])` or
  `dot(loops[0], ext[j])` numerators. This keeps the intended energy caps but
  avoids temporal-only local probe numerators for the non-vacuum cases.
- Fixed automatic numerator energy cap extraction in `gammalooprs` so EMR
  `Q(edge, energy-index)` factors only contribute caps for paired/internal
  edges when the graph is available. External-edge `Q` factors used as dense
  contraction partners no longer create invalid generation bounds.
- Added a direct unit test for the new external-edge filtering behavior in the
  energy-degree analyzer.
- Replaced the `test-cff-ltd` comparison/delta path that parsed
  `value_re/value_im` into `f64`. It now parses the formatted decimal strings,
  performs real/imaginary subtraction on decimal digit vectors first, and only
  reduces the final norm/ratio to a compact scientific display value. This
  avoids displaying exact `+0.00e0` for Quad/Arb comparisons solely because the
  difference was below double precision.

Verification completed:

```text
cargo fmt
passed

cargo check -p gammaloop-api
passed

cargo check -p gammaloop-integration-tests --tests
passed

cargo nextest run -p gammaloop-integration-tests --cargo-profile dev-optim \
  -P local_test -E 'test(test_3d_reps::cli_reconstructs_energy_power_caps_from_imported_graph_numerators) | test(test_3d_reps::cli_box_pow3_high_power_agrees_for_all_numerator_sample_normalizations) | test(test_3d_reps::cli_box_pow3_high_power_compiled_double_backends_agree_with_eager)'
3 tests run: 3 passed

cargo nextest run -p gammaloop-integration-tests --cargo-profile dev-optim \
  -P local_test -E 'test(test_3d_reps::)'
10 tests run: 10 passed

cargo nextest run -p gammalooprs --cargo-profile dev-optim \
  -P local_test energy_degree
7 tests run: 7 passed

just clippy -- -D warnings
passed
```

Transient note:

- The first attempt to run the focused `gammalooprs` dev-optim unit tests
  failed with `No space left on device` while writing target artifacts. I
  removed only the transient `target/debug` build cache, then reran the focused
  unit test successfully.

## 2026-05-01: ArbPrec visible deltas in `3Drep test-cff-ltd`

Follow-up request:

- Running `3Drep test-cff-ltd ... --precision ArbPrec` still displayed
  `+0.00e0` in the final comparison table for non-identical CFF/LTD values.

Implementation notes:

- Tightened the per-evaluation `Delta` display to use three significant digits,
  matching the final comparison table.
- Added a regression test using the concrete ArbPrec CFF/LTD decimal strings
  from the box-pow3 diagnostic run. The test verifies that the decimal
  comparison path keeps the sub-double difference visible as
  `+1.04e-299` absolute and `+3.99e-295` relative, instead of collapsing to
  zero through an `f64` boundary.

Verification completed:

```text
cargo fmt
passed

cargo test -p gammaloop-api decimal_complex_distance_keeps_sub_double_precision_differences
passed

cargo check -p gammaloop-api
passed

just build-cli
passed

./gammaloop --dev-optim ./examples/cli/scalar_topologies/three_d_expr_demo.toml run demo \
  -c "3Drep test-cff-ltd -p three_d_expr_box_pow3 -i default -g 0 --precision ArbPrec \
      --workspace-path ./examples/cli/scalar_topologies/three_d_expr_demo_workspace \
      --clean; quit -n"
passed; final table displayed abs diff +1.04e-299 and rel diff +3.99e-295
```

## 2026-05-01: `3Drep evaluate` profiling and precision-dependent comparison thresholds

Follow-up request:

- Add `3Drep evaluate --profile <duration>` to repeatedly call the evaluator on
  the same diagnostic input and report the timing per sample instead of the
  evaluation value.
- Add `3Drep evaluate --eager` to force the eager backend for Double precision.
- Make `test-cff-ltd` direct CFF/LTD relative tolerance track half of the
  requested precision, tighten pure-LTD mass-shift tolerances for Quad/ArbPrec,
  add one more default mass-shift epsilon, and expose `--n-epsilon-steps`.

Implementation notes:

- `evaluate --profile` performs a 10-call warmup, derives the number of profile
  calls from the requested target duration, evaluates the same input repeatedly,
  and records target, warmup, total, call count, and per-sample timing in the
  evaluate manifest.
- `evaluate --eager` disables compiled backend selection in the 3Drep evaluator
  settings and reports the forced eager setting in the summary.
- `test-cff-ltd` now uses direct relative tolerances `1e-7`, `1e-16`, and
  `1e-150` for Double, Quad, and ArbPrec respectively. Pure-LTD mass-shift
  relative tolerances are now `1e-2`, `1e-3`, and `1e-4` for Double, Quad, and
  ArbPrec respectively.
- The default pure-LTD mass-shift sequence now has five epsilon values
  `eps, eps/10, ..., eps/10000`; `--n-epsilon-steps` overrides that count and is
  included in cache matching and the rendered settings table.

Verification completed:

```text
cargo fmt
cargo fmt --check
cargo check -p gammaloop-api
cargo test -p gammaloop-api commands::threedreps::tests
just clippy -- -D warnings
just build-cli

./gammaloop --dev-optim -n -s ./examples/cli/scalar_topologies/three_d_expr_demo_state \
  3Drep evaluate --workspace-path ./examples/cli/scalar_topologies/three_d_expr_demo_workspace \
  --precision Double --seed 11 --scale 0.25 --profile 10ms --eager
passed; summary reported eager backend, profile calls, total time and per-sample timing

./gammaloop --dev-optim -n -s ./examples/cli/scalar_topologies/three_d_expr_demo_state \
  3Drep test-cff-ltd -p three_d_expr_box_pow3 -i default -g 0 \
  --workspace-path ./examples/cli/scalar_topologies/three_d_expr_demo_workspace \
  --precision Double --seed 11 --scale 0.25 --n-epsilon-steps 2 --clean
passed; rendered epsilon steps and two pure-LTD mass-shift evaluations

./gammaloop --dev-optim -n -s ./examples/cli/scalar_topologies/three_d_expr_demo_state \
  3Drep test-cff-ltd -p three_d_expr_box_pow3 -i default -g 0 \
  --workspace-path ./examples/cli/scalar_topologies/three_d_expr_demo_workspace \
  --precision Quad --seed 11 --scale 0.25 --n-epsilon-steps 1 --clean
passed; direct tolerance rendered as +1.00e-16 and mass-shift tolerance as +1.00e-3

./gammaloop --dev-optim -n -s ./examples/cli/scalar_topologies/three_d_expr_demo_state \
  3Drep evaluate --workspace-path ./examples/cli/scalar_topologies/three_d_expr_demo_workspace \
  --precision Double --seed 11 --scale 0.25 --profile 1ms --eager
passed after the final dev-optim rebuild; summary still reported eager backend and profile timing
```

## 2026-05-01: Relocated 3Drep demo cards

Follow-up request:

- Move the scalar box-pow3 3Drep demo card out of
  `examples/cli/scalar_topologies` into a dedicated
  `examples/cli/three_d_expr_demo` folder.
- Add a similar demo card for the `aa_aa_triple_box.dot` graph without running
  that new card yet.

Implementation notes:

- Moved the box demo to
  `examples/cli/three_d_expr_demo/box_one_edge_pow3.toml`.
- Renamed the box demo state, process, graph, workspace, and generated artifact
  paths to use `box_one_edge_pow3`.
- Added `examples/cli/three_d_expr_demo/aa_aa_triple_box.toml`, importing
  `examples/cli/three_d_expr_demo/aa_aa_triple_box.dot` through the normal
  `import graphs` file path and applying the same validate/build/evaluate/
  test-cff-ltd demo sequence with isolated state/workspace paths.

Verification completed:

```text
python3 - <<'PY'
import tomllib
from pathlib import Path
for path in [
    Path('examples/cli/three_d_expr_demo/box_one_edge_pow3.toml'),
    Path('examples/cli/three_d_expr_demo/aa_aa_triple_box.toml'),
]:
    with path.open('rb') as f:
        tomllib.load(f)
    print(f'parsed {path}')
PY
passed for both TOML cards
```

## 2026-05-01: Truncated long 3Drep evaluator diagnostics

Follow-up request:

- `3Drep evaluate` on the `aa_aa_triple_box` demo can surface very large
  Symbolica evaluator-build diagnostics, including full expressions and
  function maps. These should be middle-truncated before being propagated to
  the CLI output.

Implementation notes:

- Added a shared diagnostic formatter for `3Drep` errors that keeps the first
  200 and last 200 lines and replaces the middle with a truncation marker.
- Applied the formatter to the `3Drep evaluate` evaluator-build error boundary,
  so build errors no longer escape through `?` with their full unbounded debug
  body.
- Reused the same formatter for 3Drep generation/evaluator/evaluation error
  strings written to manifests or shown in summaries.
- Added a unit test confirming the 200-line head/tail truncation behavior.

Verification completed:

```text
cargo fmt
cargo fmt --check
cargo test -p gammaloop-api commands::threedreps::tests::diagnostic_error_truncation_keeps_first_and_last_200_lines
cargo check -p gammaloop-api
just clippy -- -D warnings
passed
```

## 2026-05-01: 3Drep physical-graph evaluator preparation

Follow-up request:

- The `3Drep evaluate` path for the `aa_aa_triple_box` demo failed while
  building/evaluating the diagnostic expression for physical graphs because
  color algebra had not been simplified before evaluator construction and the
  standalone 3Drep `ParamBuilder` setup did not populate dynamic external
  polarization parameters.

Implementation notes:

- Added a public `ParamBuilder::set_external_kinematics_and_polarizations(...)`
  helper in `gammalooprs` that reuses the same runtime external momenta and
  polarization construction machinery as the standard amplitude/cross-section
  integrand build path.
- Made `ParamBuilder::update_model_values(...)` public so 3Drep diagnostics can
  complete model-parameter population without duplicating internals in
  `gammaloop-api`.
- Updated `3Drep evaluate` and `3Drep test-cff-ltd` evaluator preparation to:
  - color-simplify the generated Symbolica atom before evaluator building;
  - resolve the graph's model from the active state;
  - populate runtime parameters, model parameters, runtime external kinematics,
    and derived polarization entries before creating evaluator inputs;
  - fall back to randomized external diagnostic inputs only for graphs without
    polarization slots when the active runtime settings do not provide a usable
    external-kinematics construction.
- The saved `param_builder.txt` for `3Drep evaluate` now reflects the prepared
  diagnostic builder used by the evaluator rather than the raw graph builder.

Verification in progress:

```text
cargo fmt
cargo check -p gammalooprs -p gammaloop-api
passed
cargo test -p gammaloop-api commands::threedreps::tests::diagnostic_error_truncation_keeps_first_and_last_200_lines
passed
just build-cli
passed
./gammaloop --dev-optim ./examples/cli/three_d_expr_demo/aa_aa_triple_box.toml run demo -c "quit -o"
passed with the currently commented one-evaluate demo card; the 3Drep evaluate
summary now shows runtime kinematics and polarization entries in the parameter
table and no longer trips over the previous color/polarization evaluator setup
failure.
just clippy -- -D warnings
passed
just build-cli
passed after the final parameter-source display tweak
./gammaloop --dev-optim ./examples/cli/three_d_expr_demo/aa_aa_triple_box.toml run demo -c "quit -o"
passed after the final parameter-source display tweak
```

## 2026-05-01: 3Drep numerator-only evaluation mode

Follow-up request:

- Add `3Drep evaluate --numerator-only` to build/evaluate only a single numerator
  call, with all internal edge energy components substituted by their positive
  on-shell energies, without introducing a fake new representation in the 3Drep
  oriented-expression layer.

Implementation notes:

- Added `numerator_with_positive_internal_ose_gs(graph)` in the gammaloop CFF
  expression bridge. It applies the same Symbolica replacement conventions used
  by 3D orientation energy maps, but only for internal edges: `Q(e,0) -> OSE(e)`,
  `Q(e,mink(...)) -> Q3(e,...) + OSE(e) delta_0(...)`, and `K(l,0)` follows the
  corresponding LMB edge. External momenta remain runtime/input parameters.
- Added `--numerator-only` to `3Drep evaluate`. In this mode the command loads
  the same artifact only to resolve the process/integrand/graph, then builds the
  evaluator from the on-the-fly numerator atom and a one-entry all-positive
  internal-edge orientation vector.
- Kept the numerator-only artifacts separate from the standard evaluate outputs:
  `numerator_only_symbolica_expression.txt`, `numerator_only_param_builder.txt`,
  and `numerator_only_evaluate_manifest.json`.
- Added an `evaluation mode` row to the evaluate summary so the output clearly
  distinguishes full 3D-expression evaluation from numerator-only evaluation.

Verification in progress:

```text
cargo fmt
cargo check -p gammalooprs -p gammaloop-api
passed
cargo test -p gammaloop-api commands::threedreps::tests::diagnostic_error_truncation_keeps_first_and_last_200_lines
passed
just build-cli
passed
./gammaloop --dev-optim -n -s ./examples/cli/three_d_expr_demo/aa_aa_triple_box_state \
  3Drep evaluate --workspace-path ./examples/cli/three_d_expr_demo/aa_aa_triple_box_workspace \
  --precision Double --seed 11 --scale 0.25 --numerator-only --clean
passed; summary reported evaluation mode `numerator-only` and saved the
separate numerator-only Symbolica expression, ParamBuilder, and manifest files.
cargo fmt --check
cargo check -p gammalooprs -p gammaloop-api
just clippy -- -D warnings
passed after the final needless-borrow cleanup
just build-cli
passed after the final needless-borrow cleanup
./gammaloop --dev-optim -n -s ./examples/cli/three_d_expr_demo/aa_aa_triple_box_state \
  3Drep evaluate --workspace-path ./examples/cli/three_d_expr_demo/aa_aa_triple_box_workspace \
  --precision Double --seed 11 --scale 0.25 --numerator-only --clean
passed after the final needless-borrow cleanup
```

## 2026-05-01: explicit 3Drep evaluate artifact selection

Follow-up request:

- Make `3Drep evaluate` explicit about which oriented expression it evaluates,
  because a workspace can contain CFF, LTD, and PureLTD artifacts for the same
  graph and the old latest-pointer behavior was too implicit.
- Update existing tests and demo cards to use the explicit form, but do not run
  the test suite in this pass.

Implementation plan:

- Add optional `-p/-i/-g` graph-selection flags to `3Drep evaluate`, plus:
  - `--json-in <path>` for an explicit oriented-expression JSON;
  - `--representation cff|ltd|pure-ltd` for cached artifact lookup;
  - `--numerator-samples-normalization ...` to select the cached CFF/LTD
    sampling-normalization artifact when `--representation` is used.
- Resolve the artifact path in this order:
  1. `--json-in`, which takes precedence over all cached lookup options;
  2. `--representation` plus `-p/-i/-g`, resolving the standard
     process/integrand/graph workspace build path;
  3. the existing `oriented_expression.json` / latest-pointer fallback for
     backwards-compatible ad-hoc use.
- Keep `--numerator-only` using the selected artifact only to identify the
  graph/model context.

Implementation notes:

- Added an optional graph selector specifically for `evaluate`, so `evaluate`
  remains usable with `--json-in` or the latest-pointer fallback while still
  supporting explicit state-based artifact lookup.
- Added clear errors when a user supplies graph-selection flags or
  `--numerator-samples-normalization` without either `--representation` or
  `--json-in`.
- Updated the CLI integration test and the `box_one_edge_pow3` /
  `aa_aa_triple_box` demo cards to specify the representation being evaluated.
- Extended the completion regression coverage so `3Drep evaluate` offers the
  same process selector and representation / normalization enum completions as
  `3Drep build`.

Verification status:

- Per the user request, the test suite has not been rerun for this pass yet.
- `cargo fmt` passed.
- `cargo check -p gammalooprs -p gammaloop-api` passed.
- `just clippy -- -D warnings` passed as a static check; no tests were executed.

## 2026-05-01: reorganized 3Drep example cards

Follow-up request:

- Reorganize `examples/cli/three_d_expr_demo` so it keeps only:
  `aa_aa.toml`, `box_one_edge_pow3.toml`, and `aa_aa.dot`.
- Keep top-level run-card `commands = []`, requiring users to explicitly run
  named command blocks.
- Add one `load` block per card and split the remaining examples into focused,
  descriptive command blocks.
- For the `aa_aa` card, duplicate the showcase blocks for 1L, 2L, 3L, and 4L,
  targeting graph ids 0, 1, 2, and 3 respectively, and add `1L_demo` through
  `4L_demo` grouping blocks.
- Verify the `1L_demo` block goes through, without running the test suite.

Implementation notes:

- Replaced the old single-block `aa_aa` and `box_one_edge_pow3` cards with
  load/build/evaluate/profile/compare-oriented command blocks.
- Kept the box-pow3 graph inline in `box_one_edge_pow3.toml`, so no extra DOT
  file is needed in the demo directory.
- Used `aa_aa.dot` as the single shared physical-graph DOT file containing the
  box, double-box, triple-box, and quadruple-box topologies.
- Fixed the `aa_aa.dot` graph metadata so each graph has its own `group_id`;
  otherwise importing all four graphs at once failed with multiple group
  masters for group 0.
- Removed generated state/workspace directories and obsolete scratch DOT/card
  files from the demo folder.
- Staged the old scalar-topologies demo-card deletion and the new
  `three_d_expr_demo` files so the move is represented in the pending commit
  set.

Verification completed:

```text
python3 tomllib parse of aa_aa.toml and box_one_edge_pow3.toml
passed; both cards have commands = [] and a load block

just build-cli
passed; needed so the local ./gammaloop binary knows the new explicit
3Drep evaluate selector flags

./gammaloop --dev-optim --clean-state ./examples/cli/three_d_expr_demo/aa_aa.toml run 1L_demo -c "quit -n"
passed after assigning unique group ids in aa_aa.dot
```

The full test suite was not run for this pass.

## 2026-05-03: remove GammaLoop-side legacy CFF mode

Follow-up correction:

- `global.generation.uv.use_legacy` is no longer a supported setting. UV
  renormalization always uses the hedge-poset forest implementation, and the old
  branch-selection setting was removed from defaults, tests, and example cards.
- The old CFF generator is no longer exposed through `gammalooprs`. Production
  CFF construction now goes through `generate_3d_expr`; denominator-only helper
  options remain available only for CFF denominator/orientation construction and
  do not form a legacy mode.
- The `gammalooprs` crate no longer has an `old_cff` feature. The optional
  comparison anchor remains only in `crates/three-dimensional-reps` behind the
  `old_cff` feature.
- Removed the Python API helper methods that imported GammaLoop-side
  subgraph-CFF generation directly.
- Dropped the now-unused old forest pole-part helper and renamed stale comments
  so the UV stack no longer carries a legacy-mode concept.

Verification completed:

```text
cargo fmt --all
cargo check -p three-dimensional-reps -p gammalooprs -p gammaloop-api --tests --locked
cargo check -p three-dimensional-reps --tests --features old_cff --locked

just test_gammaloop
Summary [169.177s] 1054 tests run: 1054 passed, 125 skipped
```

## 2026-05-02: external tree preservation in generalized CFF generation

Step III.a external-tree follow-up:

- `generate_3d_expr` now supports preserving selected internal edges as
  residual four-dimensional denominators instead of feeding them into CFF
  partial fractioning.
- The generated expression carries these residual denominators explicitly, with
  stable edge ids and powers, and remaps them alongside orientation edge maps.
- Pure-tree inputs are handled as a degenerate CFF expression: one undirected
  orientation, no CFF E-surfaces, and all internal tree propagators listed as
  residual four-dimensional denominators.
- GammaLoop's uncontracted full-amplitude CFF path now calls the generalized
  `generate_3d_expr` path. External tree edges are detected from the graph
  tree edges whose loop-momentum signature has no loop dependence; these are
  excluded from the contracted CFF source graph and reintroduced through the
  residual-denominator metadata.
- The residual denominator factor is materialized centrally in `Graph::cff`,
  including the integrated/local UV CFF paths that call through it, so the old
  manual four-dimensional tree-denominator multiplications are no longer
  duplicated in UV approximations.
- The legacy pure-CFF generator has not yet been eliminated from `gammalooprs`.
  It remains a bounded fallback for rank-deficient reduced sources,
  initial-state cut graphs, and contracted-subgraph CFFs. A strict removal pass
  was tried during this update, but it regressed the mass-approach inspect
  checks, so retiring this fallback still requires generalized support for those
  reduced/initial-cut cases.

Coverage added:

- 3D-reps unit coverage for preserving an external tree branch attached to a
  loop and for the pure-tree degenerate expression.
- Top-level end-to-end inspect comparison tests for imported scalar graphs:
  a triangle with one external tree, a box with two external trees, and a pure
  tree. Each compares normal CFF generation against
  `explicit_orientation_sum_only=true` using local inspect only.

Verification completed for this pass so far:

```text
cargo test -p three-dimensional-reps cff_generation_ --locked
cargo check -p three-dimensional-reps -p gammalooprs -p gammaloop-api --tests --locked
cargo test -p gammaloop-integration-tests imported_external_tree_inspects_match_explicit_orientation_sum --locked -- --nocapture
cargo fmt --all --check
just test_gammaloop
```

`just test_gammaloop` passed with warnings enforced:

```text
1057 tests run: 1057 passed, 128 skipped
```

The recipe now supplies `RUST_MIN_STACK=33554432` by default, while preserving a
caller-provided override, because the curated GammaLoop test suite otherwise
stack-overflows in deep symbolic generation tests on the default Rust test
thread stack.

## 2026-05-02: production numerator energy bounds for `generate_3d_expression`

Follow-up implementation after the initial Step III.c pass:

- Moved automatic numerator energy-degree bound construction from the `3Drep`
  command implementation into `gammalooprs`.
- Added `Graph::automatic_numerator_energy_degree_bounds()`, which returns a
  bound entry for every paired internal edge and overlays the detected numerator
  energy powers from `Graph::numerator_energy_power_caps()`.
- Updated `3Drep build` and `3Drep test-cff-ltd` to call the new `Graph` method
  instead of carrying a local copy of the bound-building logic in the API crate.
- Updated production `Graph::generate_3d_expression_for_integrand(...)` so the
  caller supplies `Generate3DExpressionOptions` instead of using hard-coded empty
  energy bounds.
- Added production CFF option construction from `GenerationSettings`, mapping
  `global.generation.uniform_numerator_sampling_scale` to the generalized
  `NumeratorSamplingScaleMode`.
- Updated amplitude and cross-section production CFF generation to pass automatic
  numerator energy bounds and the selected numerator-sampling scale mode into
  the generalized 3D expression generator.
- Added a guard at the legacy `Graph::generate_cff(...)` entry point: if the
  graph numerator has any detected EMR/loop-energy power above linear, legacy
  CFF generation now errors with a clear message instead of silently proceeding
  through an unsupported path.
- Added a focused regression test that constructs a nonlinear energy numerator
  and verifies the legacy CFF guard reports the unsupported edge degree.

Verification completed before the full selected suite:

```text
cargo fmt
cargo check -p gammalooprs -p gammaloop-api --tests --locked
cargo test -p gammalooprs \
  cff::generation::tests_cff::legacy_cff_rejects_nonlinear_energy_numerator \
  --locked -- --nocapture
```

## 2026-05-02: commit summary for `ltd_in_gammaloop`

This working-set commit collects the Step III production-integration work and
the supporting CLI/test/documentation changes needed to keep the branch
reviewable.

Included changes:

- Added the production Step III settings surface:
  - `global.3d_representation`, with parsed values `CFF` and `LTD`;
  - `global.generation.explicit_orientation_sum_only`;
  - `global.generation.uv.local_uv_cts_from_expanded_4d_integrands`.
- Kept unsupported Step III paths explicit:
  - production `LTD` generation still errors as not implemented;
  - `local_uv_cts_from_expanded_4d_integrands` still errors outside explicit
    orientation-sum mode.
- Added the generalized CFF bridge
  `Graph::generate_3d_expression_for_integrand(...)`, including:
  - GammaLoop-to-`ThreeDGraphSource` contracted graph extraction;
  - GammaLoop edge-id preservation through `EnergyEdgeIndexMap`;
  - conversion from generalized linear E/H surfaces back into GammaLoop concrete
    surface caches;
  - stripping generalized half-edge inverse-energy factors to preserve the
    existing production CFF convention;
  - legacy CFF fallback for initial-state-cut and contracted-subgraph cases
    until those parity boundaries are closed.
- Wired production amplitude, cross-section, and lower-level CFF call sites to
  the generalized CFF bridge while preserving the default CFF behavior.
- Implemented `explicit_orientation_sum_only = true` for production CFF:
  - graph-term atoms are explicitly summed over all selected orientations at
    generation time;
  - evaluator construction uses a direct summed evaluator path with no runtime
    orientation branching;
  - runtime settings must use `general.evaluator_method = Summed`;
  - individual-orientation Monte Carlo sampling and explicit orientation
    selection are rejected cleanly;
  - event/observable orientation ids collapse to a single effective orientation
    id `0`, displayed as all `x` markers.
- Propagated explicit orientation-sum semantics through local UV, integrated UV,
  threshold counterterm, amplitude-counterterm, and LU-counterterm evaluator
  construction.
- Added an exhaustive saved-state default-settings regression test ensuring
  `global_settings.toml` and `default_runtime_settings.toml` retain all
  schema-visible default options when `SHOWDEFAULTS` is active.
- Updated `AGENTS.md` to document `just test_gammaloop` as the default selected
  GammaLoop suite and to note that it runs with warnings as errors.
- Extended `3Drep evaluate` with:
  - `--manifest-name` for selecting the evaluate manifest file name;
  - `--no-show-parameters` for suppressing the parameter table in command
    output.
- Added the `summarize_3drep_workspace.py` demo helper and ignored generated
  nested `*_workspace/` directories in `examples/`.
- Added/updated tests for:
  - Step III setting parsing and guard behavior;
  - exhaustive default settings serialization;
  - generalized CFF bridge parity on a scalar bubble fixture;
  - inspect-level equality between ordinary CFF and explicit orientation-sum CFF;
  - explicit runtime rejection of non-`Summed`, orientation Monte Carlo, and
    explicit orientation selection;
  - stable event-order expectations in graph-grouped event tests.

Final verification before commit:

```text
just test_gammaloop
passed: 1055 tests run, 1055 passed, 128 skipped
```

## 2026-05-02: Step III kickoff - production GammaLoop integration

Step III targets the production GammaLoop integrand-generation path. The branch
for this work is:

```text
ltd_in_gammaloop
```

The intended merge target is `ltd_expr`, not `main`.

### Step III subgoals

a. Replace the production GammaLoop `generate_cff` call by the generalized
`three-dimensional-reps::generate_3d_expression` entry point, while preserving
the existing production behavior for the default CFF mode.

b. Add a global top-level `3d_representation` setting with variants:

- `CFF` (default)
- `LTD`

For the initial guard pass, `LTD` is parsed but raises an explicit
not-implemented error when production generation is reached.

c. Add a generation-level `explicit_orientation_sum_only` setting. The intended
semantics are to remove Monte Carlo sampling over individual orientations and
force an explicit sum over all orientations per sample. This is required for
LTD, while CFF should eventually support both the old sampled-orientation mode
and explicit-sum mode.

For the initial guard pass, this setting defaults to `false` and raises an
explicit not-implemented error when evaluator/generation construction is
reached.

d. Add `global.generation.uv.local_uv_cts_from_expanded_4d_integrands`. The
intended semantics are to build local UV counterterms from the expanded
four-dimensional integrand and then localize them, instead of expanding the
already generated three-dimensional expression of the original integrand.

This setting defaults to `false`, is only valid together with
`global.generation.explicit_orientation_sum_only = true`, and raises an
explicit not-implemented error when reached.

### Initial guard implementation

- Added `global.3d_representation` as a top-level serialized setting, backed by
  a Rust `ThreeDRepresentation` enum.
- Added `global.generation.explicit_orientation_sum_only`.
- Added
  `global.generation.uv.local_uv_cts_from_expanded_4d_integrands`.
- Added guard checks at process preprocessing and graph-term construction
  boundaries so unsupported Step III modes fail explicitly before silently using
  the old CFF implementation.
- Added a direct UV forest guard so direct UV calls with
  `local_uv_cts_from_expanded_4d_integrands = true` do not silently use the
  old local-3D UVCT construction.
- Added settings tests covering parsing and guard errors for the new Step III
  settings.

Verification completed:

```text
cargo fmt --all -- --check
cargo test -p gammalooprs step_iii_3d_settings_parse_and_guard_pending_modes
cargo test -p gammalooprs local_uv_4d_source_requires_explicit_orientation_sum_mode
cargo check -p gammalooprs --tests --locked
just check
just clippy -- -D warnings
just test_gammaloop
passed: 1051 tests run, 1051 passed, 128 skipped
```

### Plan for subgoal III.a only

1. Add a GammaLoop graph adapter audit for `ThreeDGraphSource`.
   The `3Drep` CLI already calls `generate_3d_expression` from GammaLoop graph
   state. The production path must reuse that adapter rather than serializing
   DOT or rebuilding graph topology. The first check is that the generated edge
   index map, mass expressions, external-energy shifts, initial-state cut
   treatment, and contracted-tree-edge handling match the current
   `Graph::generate_cff` production call.

2. Introduce a production helper next to `Graph::generate_cff`, provisionally
   `Graph::generate_3d_expression_for_integrand(...)`.
   It should accept the same contraction/canonization context as the current
   CFF helper, translate `GenerationSettings` into
   `Generate3DExpressionOptions`, call `generate_3d_expression`, and return the
   existing GammaLoop concrete alias `CFFExpression<OrientationID>` /
   `ThreeDExpression<OrientationID, Esurface, Hsurface>`.

3. Keep default behavior CFF-only for this subgoal.
   `global.3d_representation = CFF` should be the only enabled mode during
   III.a. The old `generate_cff` implementation should remain available as a
   parity fallback until the new CFF path passes existing production tests.

4. Port surface-cache compatibility deliberately.
   Production threshold subtraction and UV code currently depend on
   `graph.surface_cache.esurface_cache` and concrete GammaLoop `Esurface` /
   `Hsurface` payloads. The generalized crate returns linear-surface cache
   entries. The helper must either construct equivalent GammaLoop concrete
   surface payloads or add a narrow conversion layer so existing
   `determine_raised_esurfaces_from_expression`, residue selection, and
   threshold-counterterm builders keep seeing the same physical E-surface
   identities.

5. Rebuild numerator substitution around orientation energy maps.
   Current production CFF substitution is mostly orientation-sign based. The
   generalized expression invariant is stronger: one orientation corresponds to
   one explicit EMR edge-energy numerator map. Numerator replacement must read
   `OrientationExpression.edge_energy_map` directly and replace every EMR
   temporal component with that map. This must cover ordinary `+/- OSE`, zero
   samples, finite integer samples, external-energy shifts, and future
   uniform-scale `M` samples without special-casing the representation family.

6. Rebuild threshold residues from denominator E-surfaces generically.
   Current code selects residues by matching raised GammaLoop E-surfaces in CFF
   denominator trees. The new path must collect denominator E-surface
   occurrences directly from each `CFFVariant.denominator` tree, ignoring
   numerator-only surfaces, and map them to the same raised-cut data. This keeps
   threshold subtraction tied to actual denominator singularities rather than to
   CFF-specific orientation assumptions.

7. Generalize orientation marker handling before enabling generated expressions
   in production.
   The generalized representation uses:

   - `+` and `-` for `+/- OSE` numerator samples;
   - `0` for a genuine zero-energy numerator sample;
   - `X`/unspecified for no orientation marker, corresponding to the old
     tree-level/undirected placeholder that was previously encoded as `0`.

   Production `OrientationPattern`, selector rendering, runtime filtering, and
   integrated-UV orientation localization currently conflate "undirected" with
   `0`. III.a must separate these semantics. Internally this likely still maps
   to `Orientation::Undirected`, but user-facing patterns and labels must stop
   treating `0` as "no orientation".

8. Preserve integrated UV orientation localization.
   The existing localization picks a representative full orientation compatible
   with reduced CFF terms. With generalized expressions, compatibility must use
   denominator/support orientation constraints and tolerate `X`/unspecified
   entries, while still rejecting incompatible representatives. This should be
   covered by the existing integrated-UV tests plus targeted tests for zero
   numerator samples.

9. Add parity tests before deleting or bypassing the old CFF path.
   The minimum gate for III.a is:

   - CFF production generation for existing amplitude and cross-section tests;
   - threshold-subtraction tests using raised E-surfaces;
   - UV local/integrated tests, including orientation localization;
   - a direct comparison between old `generate_cff` and new
     `generate_3d_expression(... CFF ...)` on representative graph fixtures;
   - `just test_gammaloop`.

### Main difficulties for III.a

- The generalized crate has a neutral linear-surface model, while production
  GammaLoop still uses concrete `Esurface`/`Hsurface` payloads and a mutable
  graph-level surface cache. Maintaining stable E-surface ids is the highest
  risk for threshold subtraction.
- Numerator evaluation must stop deriving energy substitutions from compact
  orientation labels. Labels are diagnostic; the source of truth is
  `edge_energy_map`.
- `0` now has physical meaning as a zero-energy numerator sample. Existing
  GammaLoop pattern parsing and orientation-delta display historically used
  `0` for an unspecified/undirected marker, so this is an easy place to create
  silent selector bugs.
- Existing evaluator code assumes an orientation-parametric integrand evaluated
  with orientation variables. The explicit-orientation-sum mode is deliberately
  not part of III.a, but III.a must not make it harder to add later.
- Local and integrated UVCT code builds reduced local 3D terms from CFF-like
  expressions. The generalized CFF default must keep those terms equivalent
  before LTD or 4D-expanded UVCT construction is attempted.

### Subgoal III.a implementation progress

Implemented the first production bridge from GammaLoop CFF generation to the
generalized `three-dimensional-reps::generate_3d_expression` entry point while
keeping the enabled mode CFF-only.

- Added `Graph::generate_3d_expression_for_integrand(...)` next to the legacy
  `Graph::generate_cff(...)`.
- Added a contracted GammaLoop graph source adapter for the generalized
  `ThreeDGraphSource` trait. It contracts the same tree edges as the legacy
  production path, preserves GammaLoop edge ids through
  `EnergyEdgeIndexMap`, and keeps the original external-energy basis.
- Replaced the production CFF generation calls in:
  - amplitude graph generation;
  - cross-section graph generation;
  - the lower-level `Graph::cff(...)` helper.
- Added a linear-surface conversion layer from generalized
  `LinearSurfaceKind::{Esurface,Hsurface}` entries back into GammaLoop
  concrete `Esurface` / `Hsurface` cache entries. This keeps downstream raised
  E-surface discovery and threshold-residue code operating on the existing
  concrete surface ids.
- Lowered paired initial-state cut terms in generated linear surfaces into
  GammaLoop external shifts with the cross-section cut convention sign. This is
  needed for the generated CFF surface cache to contain the same s-channel cut
  E-surfaces used by `CrossSectionCut::new_from_cut_left(...)`.
- Preserved the current GammaLoop convention that inverse on-shell-energy
  factors live outside the CFF denominator tree. The generalized crate stores
  these factors in `CFFVariant::half_edges`, so the bridge strips them after
  surface remapping to avoid double-counting in the existing production
  pipeline.
- Kept numerator-energy substitution tied to the generalized
  `OrientationExpression::edge_energy_map`; this had already been introduced
  in the GammaLoop expression layer and now becomes the production source of
  truth for CFF-generated orientations.
- Added a parity regression comparing legacy CFF generation with the new
  generalized CFF bridge for a scalar bubble, both without contractions and
  with the same amplitude-production tree contractions.

Current deliberately unresolved Step III pieces:

- `global.3d_representation = LTD` remains guarded as not implemented in
  production generation.
- `global.generation.explicit_orientation_sum_only = true` remains guarded.
- `global.generation.uv.local_uv_cts_from_expanded_4d_integrands = true`
  remains guarded.
- Genuine zero-energy numerator sample markers are still only exercised inside
  the generalized crate. The production CFF bridge does not enable those modes
  yet, so the full user-facing `0` versus `X` selector split remains part of
  the later explicit-orientation/numerator-sampling work.

Verification completed for this pass:

```text
cargo fmt --all -- --check
cargo check -p gammalooprs --tests --locked
cargo clippy -p gammalooprs --tests --locked -- -D warnings
cargo test -p gammalooprs \
  cff::generation::tests_cff::generated_cff_matches_legacy_cff_for_scalar_bubble \
  --locked -- --nocapture
cargo test -p gammalooprs \
  processes::amplitude::test::generation_orientation_pattern_filters_evaluator_orientations \
  --locked -- --nocapture
cargo test -p gammaloop-integration-tests \
  cli_stateful_workflow_behaviors --test test_cli --locked -- --nocapture
cargo test -p gammaloop-integration-tests \
  cross_section_vector_spin_sum_matches_explicit_incoming_helicity_average \
  --test test_runs --locked -- --nocapture
cargo test -p gammaloop-integration-tests \
  save_state_writes_exhaustive_default_settings_files \
  --test test_cli --locked -- --nocapture
```

`cargo test -p gammalooprs --locked` still aborts in the existing
`cff::esurface::tests::failing::test_from_cut_left_box` test before reaching the
new bridge-specific tests when run as a whole.

## 2026-05-02: Step III.c explicit orientation sum mode

Implemented `global.generation.explicit_orientation_sum_only = true` for the
current production CFF mode.

Semantics implemented:

- Explicit-sum generation now builds graph-term evaluators that contain the
  direct sum over all valid CFF orientations in the generated atom, so runtime
  evaluation does not sample or branch over orientation choices.
- Runtime settings are validated for this mode:
  - `general.evaluator_method` must be `Summed`;
  - individual-orientation Monte Carlo sampling is rejected;
  - explicit runtime orientation selection is rejected.
- The effective event/observable orientation space is collapsed to one
  orientation, id `0`, rendered as an all-`x` signature.
- Generation evaluator options that only select between iterative/summed
  orientation implementations are ignored in this mode, except that an
  explicitly requested iterative generation path still remains unsupported.
- `global.generation.uv.local_uv_cts_from_expanded_4d_integrands` is accepted
  when explicit orientation summation is enabled, but still guarded as
  not-implemented without explicit orientation summation.
- Local UV and threshold counterterm evaluator construction uses the same
  explicit-sum evaluator path when this mode is active.
- Superseded on 2026-05-03: the generalized CFF production helper no longer
  keeps any GammaLoop-side legacy CFF fallback. Initial-state-cut and
  contracted-subgraph support must be handled by `generate_3d_expr` directly;
  the old comparison generator remains only in the `three-dimensional-reps`
  `old_cff` feature.

Test coverage added:

- inspect-level comparison between ordinary CFF generation and
  `explicit_orientation_sum_only = true` for scalar amplitude topologies;
- explicit runtime rejection tests for non-`Summed` evaluator mode,
  individual-orientation Monte Carlo sampling, and explicit orientation
  selection;
- settings coverage for the revised
  `local_uv_cts_from_expanded_4d_integrands` guard semantics.

Verification completed:

```text
just test_gammaloop
passed: 1055 tests run, 1055 passed, 128 skipped
```

## 2026-05-01: aa_aa profiling command blocks

Follow-up implementation:

- Extended `examples/cli/three_d_expr_demo/aa_aa.toml` so the `load` block now
  imports both:
  - the physical numerator demo graphs from `aa_aa.dot` into process `aa_aa`;
  - the scalar no-numerator variants from `aa_aa_no_numerator.dot` into process
    `aa_aa_no_numerator`.
- Added `aa_aa_no_numerator.dot`, mirroring the 1L/2L/3L/4L `aa_aa` box
  topologies with graph, edge, node, and projector numerators set to `1`.
- Added long-form profiling command blocks:
  - `1L_profile`, `2L_profile`, `3L_profile`, `4L_profile`;
  - `1L_profile_numerator_only`, `2L_profile_numerator_only`,
    `3L_profile_numerator_only`, `4L_profile_numerator_only`;
  - `1L_profile_no_numerator`, `2L_profile_no_numerator`,
    `3L_profile_no_numerator`, `4L_profile_no_numerator`.
- Each profile block builds/reuses LTD and CFF artifacts in that order and runs
  `3Drep evaluate --profile 10s` for:
  - Double eager;
  - Double compiled with assembly;
  - Double compiled with symjit;
  - Quad eager;
  - ArbPrec eager.
- The numerator-only blocks exercise the direct `--numerator-only` evaluator
  path without requiring prebuilt oriented-expression JSON.

Verification completed:

```text
cargo fmt
cargo check -p gammaloop-api
just build-cli
cargo clippy -p gammaloop-api -- -D warnings

python3 - <<'PY'
from pathlib import Path
import tomllib
with Path("examples/cli/three_d_expr_demo/aa_aa.toml").open("rb") as f:
    data = tomllib.load(f)
blocks = [b["name"] for b in data["command_blocks"]]
assert all(name in blocks for name in [
    "load",
    "1L_profile", "1L_profile_numerator_only", "1L_profile_no_numerator",
    "2L_profile", "2L_profile_numerator_only", "2L_profile_no_numerator",
    "3L_profile", "3L_profile_numerator_only", "3L_profile_no_numerator",
    "4L_profile", "4L_profile_numerator_only", "4L_profile_no_numerator",
])
PY

./gammaloop --dev-optim --clean-state \
  ./examples/cli/three_d_expr_demo/aa_aa.toml run load -c "quit -n"

./gammaloop --dev-optim --clean-state \
  ./examples/cli/three_d_expr_demo/aa_aa.toml \
  run load \
  1L_profile 1L_profile_numerator_only 1L_profile_no_numerator \
  2L_profile 2L_profile_numerator_only 2L_profile_no_numerator \
  -c "quit -n"
```

The x=1 and x=2 profile, numerator-only profile, and no-numerator profile
blocks all completed successfully. The full `just test_gammaloop` suite was not
run for this pass.

## 2026-05-01: 3Drep workspace summary helper

Follow-up implementation:

- Added the top-level helper script `summarize_3drep_workspace.py`.
- The script scans a supplied 3Drep workspace for all `evaluate_manifest.json`
  files and renders a colored `prettytable` summary with:
  - graph name;
  - representation;
  - numerator category (`no_numerator`, `only_numerator`, or `full`);
  - precision;
  - evaluator mode (`eager`, `compiled assembly`, `compiled symjit`, etc.);
  - evaluator build time in seconds;
  - profile/evaluation time per sample in microseconds.
- Timings are printed with three significant digits. When a Double/eager
  baseline exists for the same graph, representation, and numerator category,
  timing ratios are shown in parentheses.
- If multiple manifests share the same displayed category, the script adds a
  final `details` column containing the path and relevant settings needed to
  distinguish the entries.

Verification completed:

```text
python3 -m py_compile summarize_3drep_workspace.py

./summarize_3drep_workspace.py \
  ./examples/cli/three_d_expr_demo/aa_aa_workspace --no-color
```

The helper found and summarized the 13 evaluate manifests currently present in
the local `aa_aa_workspace`.

## 2026-05-02: named 3Drep evaluate manifests

Follow-up implementation:

- Added `3Drep evaluate --manifest-name <name>`.
- The option selects the evaluate summary JSON filename inside the selected
  3Drep artifact directory. The default remains `evaluate_manifest.json` for
  existing workflows. Names without a `.json` suffix receive it automatically;
  path separators are rejected so the option cannot escape the artifact folder.
- Added `3Drep evaluate --no-show-parameters`.
- The option keeps the full `parameters` array in the JSON manifest, but skips
  the large input-parameter table in the terminal output.
- Updated `examples/cli/three_d_expr_demo/aa_aa.toml` so every `3Drep evaluate`
  command writes a unique manifest name and passes `--no-show-parameters`.
- Updated `summarize_3drep_workspace.py` to scan all JSON files under a
  workspace and retain only those matching the evaluate-manifest schema. This
  lets it pick up arbitrary manifest names supplied with `--manifest-name`.

Verification completed:

```text
python3 - <<'PY'
from pathlib import Path
import tomllib
with Path("examples/cli/three_d_expr_demo/aa_aa.toml").open("rb") as f:
    data = tomllib.load(f)
names = []
for block in data["command_blocks"]:
    for cmd in block["commands"]:
        if cmd.startswith("3Drep evaluate"):
            assert "--manifest-name " in cmd
            assert "--no-show-parameters" in cmd
            parts = cmd.split()
            names.append(parts[parts.index("--manifest-name") + 1])
assert len(names) == len(set(names))
PY

./gammaloop --dev-optim --clean-state -s /tmp/gammaloop_state_named_manifest \
  ./examples/cli/three_d_expr_demo/aa_aa.toml run load \
  -c "3Drep evaluate -p aa_aa -i default -g 0 --representation cff \
      --workspace-path /tmp/gammaloop_3drep_named_manifest --precision Double \
      --seed 11 --scale 0.25 --profile 1ms --eager --clean \
      --numerator-only --manifest-name smoke_eval.json --no-show-parameters; quit -n"

./summarize_3drep_workspace.py /tmp/gammaloop_3drep_named_manifest --no-color
```

The smoke test wrote the named `smoke_eval.json` manifest, suppressed the
parameter table in the terminal output, and the summary helper discovered the
non-default manifest name.

## 2026-05-02: wrapped 3Drep summary details

Follow-up implementation:

- Updated `examples/cli/three_d_expr_demo/summarize_3drep_workspace.py` so the
  optional `details` column wraps at 40 visible characters.
- This keeps long manifest paths readable without allowing duplicate-category
  rows to make the summary table excessively wide.
- The summary now deduplicates manifests with identical visible category and
  identical detail properties.
- By default, entries that differ only in detail properties are collapsed to the
  first representative row; `--show-duplicates` expands them and highlights
  differing detail fields in red.
- Added a final `manifest` column with the manifest basename.
- Timing ratios are now measured per graph, representation, and numerator kind
  against the first available baseline in this priority order within that
  subset:
  1. Double compiled assembly;
  2. Double eager.
- The selected reference row is sorted first inside each
  graph/representation/kind subset and its manifest basename is colored green.
- Added horizontal separators between graph groups.

Verification completed:

```text
python3 -m py_compile examples/cli/three_d_expr_demo/summarize_3drep_workspace.py

./examples/cli/three_d_expr_demo/summarize_3drep_workspace.py \
  ./examples/cli/three_d_expr_demo/aa_aa_workspace --no-color

./examples/cli/three_d_expr_demo/summarize_3drep_workspace.py \
  ./examples/cli/three_d_expr_demo/aa_aa_workspace --no-color --show-duplicates
```

## 2026-05-01: aa_aa 3Drep integration test coverage

Follow-up implementation completed:

- Copied the physical `aa_aa.dot` fixture into
  `tests/resources/graphs/threedreps/aa_aa.dot` so the 3Drep tests no longer
  depend on demo files.
- Added a graph-id aware wrapper around the existing CLI
  `3Drep test-cff-ltd` manifest harness.
- Added SM-model test setup mirroring the `aa_aa` demo kinematics and
  helicities, using normal `import graphs` loading from the test resource.
- Added multi-way comparison coverage for the aa -> aa one-loop box and
  two-loop double box with:
  - Double precision + assembly compiled evaluator;
  - Double precision + symjit compiled evaluator;
  - Quad precision eager evaluator.
- Added an `evaluate` cache-reuse test for both graphs that:
  - builds and evaluates the CFF expression once with `--clean`;
  - runs `3Drep evaluate` again without `--clean` and asserts the evaluator
    build time is absent, proving the cached evaluator was reused;
  - exercises `--numerator-only` directly from the graph, with an explicit
    `--numerator-q0 4:1.25e-1` override;
  - asserts the numerator-only evaluator cache is reused on the second call.

Verification completed:

```text
cargo fmt
cargo check -p gammalooprs -p gammaloop-api -p three-dimensional-reps
cargo clippy -p gammalooprs -p gammaloop-api -p three-dimensional-reps -- -D warnings

just test_gammaloop -- cli_aa_aa
passed: 2 tests run, 2 passed, 1166 skipped

just test_gammaloop
passed: 1040 tests run, 1040 passed, 128 skipped
```

## 2026-05-01: PR review follow-up fixes

Follow-up implementation completed:

- Replaced the reserved positional inline graph import selector
  `import graphs string <DOT>` with an explicit `import graphs --inline-dot
  <DOT>` option, so a file literally named `string` can still be imported as a
  normal path.
- Updated the scalar 3Drep demo card, CLI integration tests, and completion
  tests to use/cover `--inline-dot`.
- Documented that `ParamBuilder::ensure_additional_input_parameters` must be
  called before constructing any evaluator or function map, since existing
  evaluator parameter indexing cannot be updated by mutating the builder.
- Changed repeated-mass splitting for PureLTD diagnostics to propagate mass
  evaluation errors instead of silently using a zero mass.
- Changed numerator energy-power analysis to report malformed `Q(edge, index)`
  and `K(loop, index)` arguments as structured errors instead of using sentinel
  ids.

Verification completed:

```text
cargo fmt
cargo check -p gammalooprs -p gammaloop-api -p three-dimensional-reps
cargo test -p gammaloop-api completion_offers_inline_graph_import_and_3drep_command
cargo test -p gammalooprs numerator::energy_degree
just test_gammaloop -- cli_validate_build_and_evaluate_use_gammaloop_graph_state
just test_gammaloop -- cli_graph_from_signatures_closes_without_explicit_lmb_rep
cargo clippy -p gammalooprs -p gammaloop-api -p three-dimensional-reps -- -D warnings
just test_gammaloop
passed: 1042 tests run, 1042 passed, 128 skipped
```

## 2026-05-01: numerator-only q0 inputs and evaluator-cache work

Follow-up request:

- In `3Drep evaluate --numerator-only`, do not replace internal edge energies
  with `+OSE(edge)`. Keep them as explicit evaluator parameters using the
  natural Symbolica form `Q(edge, spenso::cind(0))`.
- Add a way for the user to override those numerator-only energy inputs; random
  values must otherwise be assigned after the usual parameters so existing
  seeded inputs remain stable.
- Cache 3Drep evaluators themselves in the workspace when `--clean` is not
  supplied, reusing the normal `EvaluatorStack` serialization machinery. For
  `symjit`, reload the serialized evaluator and call `jit_compile()` again,
  since serialized JIT payloads are not available yet.
- Add a TODO for top-level graph display: introduce a `simplified_numerator`
  field next to `numerator`, and use the expression obtained after
  `simplify_metrics()` and `simplify_color()` as the evaluator input.

Implementation notes:

- Added GammaLoop CFF helpers that replace EMR/LMB temporal components by
  explicit internal-energy parameter atoms instead of positive OSEs:
  `numerator_with_internal_energy_parameters_gs`,
  `internal_energy_parameter_atom_gs`, and
  `internal_energy_parameter_replacements_gs`.
- Added `ParamBuilder::ensure_additional_input_parameters`, so diagnostic-only
  input parameters can be registered without hand-editing ranges from the CLI
  layer.
- Added `3Drep evaluate --numerator-q0 EDGE:VALUE` (also accepting comma lists
  and `EDGE=VALUE`) for numerator-only energy overrides.
- `--numerator-only` now appends all internal `Q(edge, cind(0))` parameters to
  the ParamBuilder, assigns seeded diagnostic values to them by default, and
  applies user overrides last. The parameter table marks overridden entries as
  `user override`.
- Added `LoopExtSignature::is_loop_dependent()` and use it to ensure
  `--numerator-q0` overrides apply only to loop-dependent internal edges, not
  external energy parameters.
- Evaluator construction now receives the exact ParamBuilder used for the
  evaluation, which is necessary for numerator-only q0 parameters.
- Added binary evaluator caches under the local 3Drep workspace
  `evaluators/` subdirectory, with a small JSON manifest keyed by evaluator
  name, precision/backend settings, parameter list, input atoms, and orientation
  data. Cached loads leave the evaluator build-time field empty and emit the
  exact cache paths used.
- Added `EvaluatorStack::activate_cached_f64_backend`, which reactivates eager,
  symjit, or external compiled backends after deserialization.
- Updated the box-pow3 demo numerator-only command to showcase explicit
  `--numerator-q0` overrides.

Display TODO recorded:

- Add `simplified_numerator` to the top-level graph/build summary table and
  ensure it prints the `simplify_metrics().simplify_color()` expression used for
  evaluator construction, while keeping the original `numerator` field for
  comparison. This was implemented in the subsequent display pass.

Verification completed:

```text
cargo fmt
cargo check -p gammalooprs -p gammaloop-api
cargo clippy -p gammalooprs -p gammaloop-api -- -D warnings
just build-cli

./gammaloop --dev-optim --clean-state -s /tmp/gammaloop_state_q0_test \
  ./examples/cli/three_d_expr_demo/aa_aa.toml run load \
  -c "3Drep evaluate -p aa_aa -i default -g 0 \
      --workspace-path /tmp/gammaloop_3drep_q0_test --precision Double \
      --seed 11 --scale 0.25 --profile 1ms --eager --clean \
      --numerator-only --numerator-q0 4:0.125; quit -n"
passed; the Symbolica expression contains q₄[0] style parameters rather than
OSE substitutions and the parameter table reports q4[0] as a user override.

./gammaloop --dev-optim --clean-state -s /tmp/gammaloop_state_q0_test2 \
  ./examples/cli/three_d_expr_demo/aa_aa.toml run load \
  -c "3Drep evaluate -p aa_aa -i default -g 0 \
      --workspace-path /tmp/gammaloop_3drep_q0_test --precision Double \
      --seed 11 --scale 0.25 --profile 1ms --eager \
      --numerator-only --numerator-q0 4:0.125; quit -n"
passed; the second run reused `numerator_only.evaluator.bin` and its evaluator
manifest from the 3Drep workspace.
```

The full test suite was not run for this pass.

## 2026-05-01: numerator-only evaluate no longer requires cached 3Drep JSON

Follow-up request:

- `3Drep evaluate --numerator-only -p ... -i ... -g ...` should work even when
  no CFF/LTD oriented-expression JSON has been built yet. The selected graph is
  enough context for the numerator-only evaluator.
- Keep the numerator-only Symbolica expression artifact available for user
  inspection, like the standard evaluate artifacts.

Implementation notes:

- Split `3Drep evaluate` input resolution into:
  - direct graph selection for `--numerator-only` when `-p/-i/-g` is supplied;
  - oriented-expression JSON loading for standard evaluation, or for
    numerator-only only when the user intentionally relies on `--json-in` /
    latest-artifact context.
- In direct numerator-only mode, `--representation`, `--json-in`, and cached
  oriented-expression lookup are ignored because they are not needed to build
  the trivial positive-OSE numerator call.
- Store direct numerator-only artifacts under the deterministic graph workspace
  subfolder:
  `process_.../graph_.../numerator_only/`.
- Use the standard artifact filenames in that folder:
  `symbolica_expression.txt`, `param_builder.txt`, and
  `evaluate_manifest.json`.
- The rendered evaluate summary now shows `oriented expression = not used for
  numerator-only` when no JSON was loaded.

Verification completed:

```text
cargo fmt
cargo check -p gammalooprs -p gammaloop-api
just build-cli

./gammaloop --dev-optim --clean-state ./examples/cli/three_d_expr_demo/aa_aa.toml run load \
  -c "3Drep evaluate -p aa_aa -i default -g 3 --representation cff \
      --numerator-samples-normalization M_for_beyond_quadratic_only \
      --workspace-path ./examples/cli/three_d_expr_demo/aa_aa_workspace/1L \
      --precision Double --seed 11 --scale 0.25 --profile 1ms --eager \
      --clean --numerator-only; quit -n"
passed; no cached oriented-expression JSON was required, and the summary pointed
to `.../graph_0003_QuadrupleBox/numerator_only/symbolica_expression.txt`.

just clippy -- -D warnings
passed
```

The full test suite was not run for this pass.

## 2026-05-01: simplified numerator display

Follow-up implementation:

- `3Drep build` now passes both the original graph numerator and a
  `simplified_numerator` string to the pretty renderer.
- The simplified numerator is computed with the same simplification order used
  for evaluator preparation in this path: `simplify_metrics()` followed by
  `simplify_color()`.
- The top-level pretty summary now has both `numerator` and
  `simplified_numerator` rows, so users can inspect what expression is being
  fed into evaluator construction after metric/color simplification.

Verification completed:

```text
cargo fmt
cargo check -p three-dimensional-reps -p gammaloop-api
cargo clippy -p three-dimensional-reps -p gammaloop-api -- -D warnings
just build-cli

./gammaloop --dev-optim --clean-state -s /tmp/gammaloop_state_simpl_num \
  ./examples/cli/three_d_expr_demo/box_one_edge_pow3.toml run load \
  -c "3Drep build -p box_one_edge_pow3 -i default -g 0 \
      --representation cff --workspace-path /tmp/gammaloop_3drep_simpl_num \
      --no-color --clean; quit -n"
passed; the build summary shows both `numerator` and `simplified_numerator`.
```

The full test suite was not run for this pass.

## 2026-05-03: final Step III cleanup pass before staging

Global review outcome:

- Confirmed the retired `generate_cff` entry points are not present in active
  GammaLoop source anymore. The only `old_cff` gate left is in
  `crates/three-dimensional-reps`, and it is disabled by default.
- Removed the GammaLoop crate-level `old_cff` feature so production GammaLoop
  cannot accidentally enable or route through the old generator.
- Renamed current GammaLoop CFF production wrappers from `generate_cff...` to
  `build_cff_expression...` so the remaining production path is clearly a
  `generate_3d_expr` bridge, not the retired generator.
- Removed the Python API's direct subgraph-CFF helpers. The Python API no
  longer imports the removed GammaLoop-side subgraph generator.
- Moved `SurfaceCache` imports to `cff::surface`; `cff::generation` no longer
  re-exports surface storage.
- Removed stale `use_legacy` UV settings, example comments, and the old UV
  branch. UV renormalization now has a single hedge-poset path.
- Dropped the unused old forest pole-part helper that became dead code after
  removing the UV branch.
- Removed the last dead GammaLoop-side CFF graph-contraction implementation
  (`CFFGenerationGraph`, its edge/vertex types, and its self-contained tests).
  The module now contains only the compact `VertexSet` bitset still used by
  surface bookkeeping.
- Renamed that surviving surface bitset module from `cff_graph` to
  `vertex_set` and re-exported `VertexSet` through `cff`, avoiding imports that
  suggest a remaining GammaLoop-side CFF graph builder.
- Removed the unregistered `benches_old` Criterion files that still called the
  retired `Graph::generate_cff()` API. Keeping them would have left broken
  source files in the tree and an implicit second CFF construction surface.
- Added a public GammaLoop-side `Graph::three_d_expression_options(...)` helper
  for the standalone `3Drep` CLI so auto-detected numerator energy bounds,
  initial-state-cut exclusions, and preserved 4D tree denominators are selected
  from the same place as production CFF generation.
- Centralized preserved tree-edge handling in the GammaLoop
  `generate_3d_expression` bridge. Production callers no longer pass external
  tree edges as contracted edges; the bridge selects preserved 4D denominators,
  filters them out of virtual graph contraction, and the denominator-only CFF
  wrapper explicitly excludes them from GammaLoop's surrounding inverse-OSE
  normalization.
- Changed the loop-signature rank helper used while extracting reduced
  generalized-3D graph sources from a floating-point row reduction to exact
  rational row reduction.
- Replaced derived `from_f64` constants in CFF/evaluator/subtraction hot paths
  with active-precision constants: exact integer signs now use representative
  `.from_i64(...)` values, derivative factorial normalization uses
  `.from_usize(...)`, and orientation/override zeros use `T::new_zero()`.
- Unified GammaLoop graph extraction for generalized 3D generation behind a
  single `GraphThreeDSource` adapter. The direct `Graph` API and production CFF
  bridge now share dummy-edge filtering, initial-state-cut externalization,
  contracted-edge handling, and active loop-energy column reduction.
- Tightened automatic numerator energy-bound detection and generalized-3D graph
  extraction so dummy/non-denominator edges do not participate in numerator
  energy analysis or parsed 3D graph sources.
- Normalized logbook feature spelling from `old_cff`'s earlier dashed spelling
  to the current Cargo feature name.
- Simplified the preserved-edge graph contraction helper to avoid repeatedly
  cloning the union-find parent map while computing roots.

Design observations:

- Current CFF production in `gammalooprs` always calls
  `three_dimensional_reps::generate_3d_expression`; there is no fallback to the
  retired CFF generator.
- GammaLoop no longer carries a hidden old-CFF graph builder behind
  `allow(dead_code)`, so the remaining CFF-related code is limited to surface
  representation/evaluation and the generalized 3D-expression bridge.
- A repository-wide source scan now finds no active `generate_cff` references
  outside historical logs/artifacts and the explicitly gated old-CFF mode in
  `three-dimensional-reps`.
- `3Drep build` and `3Drep test-cff-ltd` no longer hand-construct incomplete
  generalized-3D options; they reuse GammaLoop graph defaults and then apply
  only user-supplied energy-bound overrides.
- Standalone 3Drep CLI graph validation/build paths now use the same GammaLoop
  parsed-source extraction rules as production CFF generation.
- Denominator-only CFF options still exist for UV/threshold denominator
  construction, but they are explicit zero-degree generalized-3D options, not a
  legacy mode.
- The only old-CFF comparison anchor is `RepresentationMode::OldCff` in the
  `three-dimensional-reps` crate, compiled only with `--features old_cff`.

Verification completed:

```text
cargo fmt --all
cargo check -p three-dimensional-reps -p gammalooprs -p gammaloop-api --tests --locked
cargo clippy -p three-dimensional-reps -p gammalooprs -p gammaloop-api --tests --locked -- -D warnings
cargo check -p three-dimensional-reps --tests --features old_cff --locked
git diff --check

just test_gammaloop
Summary [152.605s] 1039 tests run: 1039 passed, 125 skipped
```

## 2026-05-03: final bridge cleanup, detailed test target, and state-artifact hardening

Follow-up cleanup after the previous global pass:

- Added `just test_gammaloop_detailed`, an alias for the curated GammaLoop
  suite with nextest reporter flags:
  `--show-progress=bar --status-level=slow --final-status-level=fail
  --max-progress-running=8`.
- Moved the GammaLoop-to-3D graph adapter into
  `crates/gammalooprs/src/graph/three_d_source.rs`, keeping `graph/mod.rs`
  focused on graph data and general graph operations.
- Kept graph contraction virtual inside `GraphThreeDSource`. I checked the
  linnet contraction helpers during this pass; they mutate/reindex graph
  structure, while this adapter must preserve original `EdgeIndex` values for
  energy caches, residual denominators, and orientation vectors.
- Made `GraphThreeDSource` own the initial-state-cut edge set and refuse to
  virtually contract those edges even if a future caller accidentally includes
  them in a contracted-edge request.
- Normalized preserved 4D denominator selection in one helper on the GammaLoop
  generation bridge. Initial-state-cut edges and preserved 4D denominator edges
  are now treated uniformly as non-contractible before the reduced source graph
  reaches `three_dimensional_reps::generate_3d_expression`.
- Adjusted denominator-only CFF construction so preserved 4D denominator edges
  are also excluded from the separate inverse on-shell-energy product assembled
  by GammaLoop around the generated expression.
- Replaced several exact zero/one/two/three/four constants near the new path
  with active-precision constructions (`T::new_zero()`, `.one()`,
  `.from_i64(...)`, `.from_usize(...)`) so internally derived constants do not
  pass through `f64` in generic-precision evaluation paths.
- Hardened the mass-approach run card with an explicit
  `tests/artifacts/mass_approach_scalar_self_energy` state folder and removed
  stale top-level state artifacts that had been left from an earlier ad hoc run.

Verification completed:

```text
cargo fmt --all
cargo check -p gammalooprs -p gammaloop-api --tests --locked
cargo clippy -p gammalooprs -p gammaloop-api --tests --locked -- -D warnings
git diff --check

just test_gammaloop_detailed
Summary [153.382s] 1039 tests run: 1039 passed, 125 skipped
```
