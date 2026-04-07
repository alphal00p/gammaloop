# Localizing The Integrated UV CT Over Existing Graph Orientations

## Goal

When GammaLoop generates an integrated UV counterterm, the integrated term should not be duplicated over all valid orientations of the subgraph that has been integrated out. Instead, for each valid orientation of the reduced graph, the integrated UV term should be dispatched to exactly one representative full-graph orientation that is already part of the canonical global acyclic orientation basis.

This implementation deliberately does **not** create new orientations. It reuses the orientation list that is already generated for the full graph and already consumed by evaluator generation, sampling, export, and load.

## Why The Orientation List Is Not Modified

The canonical full-graph orientation list is fixed early in amplitude preprocessing:

1. `AmplitudeGraph::generate_cff()` builds the full-graph acyclic orientation basis.
2. `AmplitudeGraph::build_integrands()` constructs UV forests and parametric integrands from that fixed basis.
3. The process evaluator copies the same orientation list into the compiled evaluator stack.

That same list is then assumed by:

- orientation parameter allocation
- evaluator construction
- amplitude / cross-section export
- amplitude / cross-section load
- discrete orientation sampling

Because of that, adding extra partially-undirected orientations later would be a much larger architectural change. The current implementation therefore localizes the integrated UV symbolically while keeping the global orientation basis unchanged.

## Matching Rule

The integrated UV term is built from a reduced graph where the edges of the integrated subgraph have been contracted. In the reduced-graph CFF generation, those contracted edges are marked as `Undirected`.

For a reduced valid orientation, the compatible full-graph candidates are therefore the existing valid global orientations that agree with it on every edge whose reduced orientation is not `Undirected`.

Equivalently:

- edges external to the integrated subgraph must match
- contracted edges are ignored when matching representatives

## Deterministic Representative Selection

Among all compatible full-graph candidates, we choose one deterministic representative by maximizing, in order:

1. the number of `Undirected` entries on the internal edges of the integrated subgraph
2. then the number of `Default` entries on those same edges
3. then the lexicographically maximal `is_default` pattern when the internal edge ids are ordered increasingly

The first criterion is currently future-proof only. In the present code, the canonical full-graph acyclic orientation basis is built from `Default` / `Reversed` assignments on paired internal edges, so `Undirected` entries are not expected there unless the global orientation-generation step is changed in the future.

## Avoiding Duplicate Selectors

The reduced-graph CFF factor already carries orientation selectors for the edges that remain explicit in the reduced graph. Since contracted edges are marked `Undirected`, they do not contribute selectors there.

The integrated-UV localization therefore adds selectors **only** for the internal edges of the integrated subgraph. It never adds a full `orientation_thetas()` factor for the representative full orientation, because that would duplicate selectors on edges that are already selected by the reduced-graph CFF term.

## Insertion Point

The localization is inserted in the integrated-UV assembly path, not in the global evaluator orientation machinery.

Concretely, the reduced-graph CFF is kept orientation-resolved for the integrated branch long enough to attach the representative selector before summing over reduced orientations.

The ordinary local CFF path remains unchanged.

## Current Implementation

The implementation lives in:

- `crates/gammalooprs/src/uv/approx/orientation_localization.rs`
- `crates/gammalooprs/src/uv/approx/mod.rs`
- `crates/gammalooprs/src/uv/settings.rs`

The helper module provides three pieces:

1. `select_representative_orientation(...)`
2. `internal_orientation_selector(...)`
3. `localize_reduced_orientation_term(...)`

The UV assembly path now builds the integrated reduced factor through a dedicated helper in
`uv/approx/mod.rs`:

- `internal_paired_edges_of_subgraph(...)`
- `localized_integrated_reduced_factor(...)`

That helper:

1. calls `graph.cff(...)` on the reduced graph
2. keeps the individual reduced orientations instead of collapsing them immediately
3. localizes each reduced orientation term onto one representative full orientation
4. only then sums the localized terms back into one atom per reduced CFF term

The generated forest call chain was updated so that the fixed global orientation basis is threaded
down explicitly into the UV approximation code:

- `AmplitudeGraph::build_integrands(...)`
- `AmplitudeGraph::build_threshold_counterterm_parametric_integrand(...)`
- `AmplitudeGraph::renormalization_part(...)`
- `CrossSectionProcess` UV-forest call sites
- `CutForests::compute(...)`
- `Forest::compute(...)`
- `Approximation::root(...)`
- `Approximation::compute(...)`
- `Approximation::final_integrand(...)`

This keeps the change local to UV generation while reusing the canonical global orientation basis
that the evaluator stack already knows about.

## Integrated-Only Debug Generation Mode

The earlier ad hoc environment-variable hook has been replaced by a real UV generation setting:

```toml
[cli_settings.global.generation.uv]
generate_only_integrated_uv_with_n_loops = 1
```

or interactively:

```text
set global kv global.generation.uv.generate_only_integrated_uv_with_n_loops=1
```

This mode is only active when integrated UV terms are being generated (`generate_integrated=true`).
When enabled, UV-forest generation suppresses every non-integrated contribution, including the
root/original graph contribution, and keeps only forest terms whose chain consists exclusively of
integrated UV counterterms with the requested loop count.

Implementation-wise:

- the root is treated as a neutral ancestor so that first-level integrated UV terms can still be kept
- the root final integrand is forced to zero
- for matching descendants, the local UV part is forced to zero and only the integrated branch is kept
- descendants whose integrated-CT chain contains a different loop count are forced to zero

This makes the feature usable as a stable debugging mode without relying on environment variables or
temporarily replacing the raw integrated counterterm by `1`.

## Tests

The standalone helper module should cover:

- the toy dispatch example discussed during debugging
- the `Undirected` / `Default` / leading-`Default` tie-breaking rule
- the absence of representatives for nonexistent external-orientation classes
- the fact that the added selector only depends on integrated-subgraph internal edges

The standalone helper module now contains those unit tests.

The standalone integration-test coverage now lives in:

- `tests/tests/test_runs/test_integrated_uv_cts.rs`

It contains two fully standalone scalar-bubble tests:

- `scalar_bubble_inspect`
- `scalar_bubble_integrated`

Those tests do not depend on `examples/cli/scalar_topologies/bubble.toml`. They build the bubble
processes directly inside the test workspace, set the same kinematics/runtime settings, and compare
against fixed inspect/integration benchmarks.

The remaining end-to-end validation target is:

```bash
./gammaloop --clean-state ./examples/cli/scalar_topologies/bubble.toml run generate integrate_bubble_below_threshold -c "quit -n"
```

Before this change, the integrated UV term for the bubble is duplicated over both valid bubble orientations, and the current result is:

```text
bubble@scalar_bubble_below_thres = -9.272(86)e-3
```

After the localization, that result should change, and the integrated UV contribution in the standalone/exported expression should appear on only one representative orientation per external-orientation class.

## Current Verification

Using a temporary copy of the bubble run card with the same process definition, the same seed, and
the same 1M-sample first iteration, but with the state folder redirected to `/tmp` to avoid
clashing with an older long-running workspace, the localized implementation gives:

```text
bubble_no_integrated_UV@scalar_bubble_below_thres = +1.4737(23)e-2
bubble@scalar_bubble_below_thres                  = +2.764(65)e-3
```

This replaces the previous wrong-sign / duplicated-orientation behavior:

```text
bubble@scalar_bubble_below_thres = -9.272(86)e-3
```

The standalone export confirms the structural change as well: the integrated UV branch now sits
only under the `σ(2)=σ(3)=+` selector for the bubble, while the local UV structure still contains
both bubble orientations.

## How To Run The Relevant Tests

Helper-unit tests for the localization logic:

```bash
cargo test -p gammalooprs orientation_localization
```

Standalone inspect regression:

```bash
cargo test -p gammaloop-integration-tests --test test_runs scalar_bubble_inspect -- --nocapture
```

Standalone integrated regression:

```bash
cargo test -p gammaloop-integration-tests --test test_runs scalar_bubble_integrated -- --nocapture
```

The integrated regression intentionally runs a full 1M-sample benchmark to compare against the same
targets used in `bubble.toml`, so it is substantially slower than the inspect regression.
