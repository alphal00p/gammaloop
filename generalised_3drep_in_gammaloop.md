# Generalised 3D Representations in GammaLoop

## Scope

GammaLoop builds production three-dimensional integrands through
`three_dimensional_reps::generate_3d_expression`.  The supported production
representations are:

- `CFF`: the default causal representation.
- `LTD`: the loop-tree-duality representation used only with explicitly summed
  orientations.

The old GammaLoop-local `generate_cff` implementation is not part of the
production path.  The feature-gated `old_cff` generator lives only inside the
`three-dimensional-reps` crate for crate-local parity checks.

The implementation is intentionally representation-driven: graph conversion,
energy-bound detection, external-tree preservation, initial-state cut handling,
UV counterterm construction, threshold residues, and runtime evaluator assembly
all consume the same oriented expression model.

## Main Settings

### Global Representation

The top-level setting is:

```toml
3d_representation = "CFF" # or "LTD"
```

`CFF` is the default.  `LTD` is valid only with
`global.generation.explicit_orientation_sum_only = true`.

### Generation Settings

The relevant generation fields are under `[generation]`:

```toml
[generation]
explicit_orientation_sum_only = false
uniform_numerator_sampling_scale = "none"
orientation_pattern = []
```

- `explicit_orientation_sum_only`
  - `false` by default.
  - Required for `3d_representation = "LTD"`.
  - Represents the full orientation sum as a single generated contribution.
  - Rejects a non-default `orientation_pattern`.
- `orientation_pattern`
  - Supported for ordinary CFF orientation-localized generation.
  - Unset in explicit-orientation-sum mode.
- `uniform_numerator_sampling_scale`
  - Serialized values: `none`, `beyond_quadratic`, `all`.
  - `none`: no uniform numerator sampling scale is introduced.
  - `beyond_quadratic`: use the runtime scale only for energy-degree sectors
    above quadratic.
  - `all`: use the runtime scale for every non-constant energy-degree sector.

The runtime scale value is:

```toml
[runtime.general]
numerator_sampling_scale = 1.0
```

This value is an explicit runtime `f64` setting boundary.  Precision-generic
evaluation code must not create derived constants by downcasting through `f64`.

### UV Settings

The relevant UV-generation fields are under `[generation.uv]`:

```toml
[generation.uv]
subtract_uv = true
generate_integrated = true
local_uv_cts_from_expanded_4d_integrands = false
```

- `subtract_uv` controls local UV subtraction.
- `generate_integrated` controls integrated UV counterterm generation.
- `local_uv_cts_from_expanded_4d_integrands`
  - `false`: build local UV counterterms from 3D expansions where supported.
  - `true`: build local UV counterterms by first expanding the original 4D
    integrand and then applying the selected 3D representation to each local
    term.

For cross-section supergraphs:

- `CFF + local_uv_cts_from_expanded_4d_integrands = false` is supported.
- `CFF + local_uv_cts_from_expanded_4d_integrands = true` is supported.
- `LTD + local_uv_cts_from_expanded_4d_integrands = true` is supported.
- `LTD + local_uv_cts_from_expanded_4d_integrands = false` is not supported.

`LTD + local UV from 3D expansions` is intentionally unsupported.  Its local
counterterms would be tied to LTD dual/H-surface structure rather than to the
representation-neutral local 4D UV forest, with no practical benefit over the
expanded-4D construction.

### Threshold Settings

The threshold settings are under `[generation.threshold_subtraction]`:

```toml
[generation.threshold_subtraction]
enable_thresholds = true
check_esurface_at_generation = false
skip_thresholds_that_are_cuts = true
disable_integrated_ct = false
```

Physical threshold counterterms are selected from canonical CFF-like
`E`-surface groups.  LTD dual `H`-surfaces are not threshold surfaces; they are
spurious for the threshold-subtraction problem after the explicit orientation
sum.

### Runtime Evaluator Mode

Explicit-orientation-sum generation requires:

```toml
[runtime.general]
evaluator_method = "Summed"
```

When `explicit_orientation_sum_only = true`:

- `SingleParametric`, `Iterative`, and `SummedFunctionMap` runtime modes are
  rejected for production evaluation.
- Individual-orientation Monte Carlo sampling is rejected.
- Explicit orientation selection in the evaluation input is rejected.
- Orientation-specific generation settings such as `summed` and
  `summed_function_map` are ignored by construction.

## Graph Conversion

GammaLoop converts its `Graph` objects into the `ThreeDGraphSource` interface
implemented in `crates/gammalooprs/src/graph/three_d_source.rs`.  The
three-dimensional crate receives:

- internal denominator edges with momentum signatures and masses,
- external half-edge information,
- initial-state cut-edge annotations,
- preserved 4D denominator edge annotations,
- numerator energy-degree bounds,
- energy-edge index remapping back to GammaLoop edge ids.

Graph reconstruction from momentum signatures is a CLI/testing convenience in
the `3Drep` command.  It is not used by GammaLoop production generation or by
the CFF lower-sector builder.

### External Trees and Pure Trees

External tree structures attached to loop graphs are preserved as residual 4D
denominators.  They are not partial-fractioned by CFF and are not cut by LTD.
Their energy components are fixed by the external momentum decomposition.

Pure tree graphs are handled as the zero-loop limit of the same mechanism: the
generated 3D expression contains the residual 4D denominator product and no
ordinary loop-energy residues.

### Initial-State Cut Edges

Forward-scattering cross-section graphs can contain internal edges that encode
initial-state cuts.  GammaLoop passes these as initial-state cut edges rather
than ordinary internal denominators.  The 3D-representation crate externalizes
them as correlated external energy shifts, preserving the sign of the
corresponding half-edge assignment.

This makes a forward-scattering graph behave like a graph with external legs on
both sides of the initial-state cut, while still using the original GammaLoop
graph as the source of momentum routing.

## Energy-Degree Bounds

GammaLoop automatically detects numerator energy-degree bounds from the EMR
numerator before calling `generate_3d_expression`.

For production CFF:

- Bounds are assigned to active denominator edges that actually participate in
  the 3D representation.
- Initial-state cut edges and preserved external-tree/4D denominator edges are
  excluded.
- The generator first attempts loop-momentum-basis bounds when they are
  convergent and otherwise falls back to source-edge bounds.
- Empty bounds are treated as trivially convergent.

For production LTD:

- The energy-degree bound vector is still recorded and passed through the same
  interface.
- LTD generation itself does not require CFF contact completion for ordinary
  non-repeated dual residues.
- Repeated-propagator and high-energy sectors use the generalized
  derivative-free construction in `three-dimensional-reps`.

The energy-bound convergence check rejects bound choices with non-vanishing
residues at infinity.  This is a genuine mathematical boundary, not a fallback
to the removed GammaLoop-local CFF implementation.

## Oriented Expression Model

The central representation type is:

```rust
ThreeDExpression<OrientationID>
```

It contains:

- a graph metadata block,
- a surface cache,
- a list of orientation expressions,
- residual 4D denominators,
- numerator energy maps,
- denominator variants grouped by shared numerator map,
- half-edge prefactors,
- numerator-only helper surfaces,
- factorized denominator surface trees.

An orientation corresponds to one EMR numerator call.  Its variants are summed
inside that numerator call.  In explicit-orientation-sum mode the outer
orientation selector is omitted and all orientations contribute directly to one
summed evaluator.

Orientation markers use:

- `+` and `-` for ordinary directed orientation choices,
- `0` for the explicit-sum effective orientation where applicable,
- `X` for no orientation marker, for example external/preserved 4D tree
  denominators.

## CFF Generation

The production CFF path calls `generate_3d_expression` with
`RepresentationMode::Cff`.

Supported features:

- ordinary CFF causal denominator trees,
- explicitly summed orientations,
- orientation-localized CFF generation when explicit summing is disabled,
- higher EMR energy powers through generalized bounded CFF contact sectors,
- repeated-signature denominators,
- preserved external tree denominators,
- initial-state cut-edge externalization,
- automatic numerator-bound detection,
- local UV from 3D expansions,
- local UV from expanded 4D integrands,
- integrated UV counterterms,
- threshold and LU cut residue selection through canonical `E`-surface groups.

The generalized CFF higher-energy algorithm keeps denominator surfaces
`E`-surface-only.  Non-repeated high-power contacts are reduced through
lower-sector CFF recursion.  Repeated-channel high-power contacts use the
channel/contact normal form.  The implementation can become expensive for large
bound vectors, but no production GammaLoop fallback uses the removed local
`generate_cff`.

The feature-gated `old_cff` generator rejects any numerator bound outside the
legacy affine/linear energy regime.

## LTD Generation

The production LTD path calls `generate_3d_expression` with
`RepresentationMode::Ltd` and requires explicit orientation summing.

Supported features:

- LTD residues for ordinary internal denominator edges,
- derivative-free repeated-propagator handling through the generalized
  repeated-channel construction,
- high EMR energy powers in accepted UV-convergent sectors,
- preserved external tree denominators,
- initial-state cut-edge externalization,
- threshold/LU cut alignment through the CFF `E`-surface projection catalogue,
- local UV from expanded 4D integrands,
- integrated UV counterterms folded into the same spatial measure as CFF.

Unsupported production modes:

- LTD with individual-orientation Monte Carlo generation.
- LTD with a non-trivial generation `orientation_pattern`.
- LTD local UV counterterms from 3D expansions.

`PureLTD` is a diagnostic mode in the 3D-representation crate and CLI.  It is
not the production formula for unsplit repeated propagators with numerator
dependence because ordinary unsplit LTD would require numerator derivatives.

## Threshold and LU Cut Uniformisation

GammaLoop uses canonical `E`-surface groups as the representation-neutral
description of physical threshold and LU cut surfaces.

For CFF:

- the generated expression exposes causal `E`-surfaces directly;
- threshold and LU residue selectors operate on these surfaces;
- cut-edge sets restrict residues to the appropriate Cutkosky cut structure.

For LTD:

- the expression may contain additional dual `H`-surfaces;
- `H`-surfaces are ignored for threshold subtraction;
- a companion CFF projection expression is built from the same graph data;
- raised LU cut groups and left/right threshold residues are selected from the
  CFF projection catalogue;
- the selected residue data are mapped back to the requested LTD integrand
  before local numerator maps and event weights are assembled.

Surface canonicalization normalizes signs such as
`-E_1 - E_2 - ... + p0_shift` into
`E_1 + E_2 + ... - p0_shift`.  The sign introduced by this canonicalization is
carried by the residue/surface map.  It is not patched per graph.

Two global convention corrections remain:

- LTD LU residues carry the loop-measure parity `(-1)^(L-1)`, with `L` the
  forward-scattering loop count after initial-state cut edges are removed.
- Expanded-4D CFF local UV residues are corrected by the difference between the
  full-graph and reduced-source generalized-CFF sign exponents.  This accounts
  for duplicate-signature excess in reduced UV sources.

Both corrections are formulaic and representation-convention based.  The code
contains no graph-label sign patches.

## UV Counterterms

### Local UV from 3D Expansions

This path is available for CFF.  The local UV forest is expanded directly in
the 3D representation.  It preserves the existing CFF local cancellation
structure and remains the default CFF local-UV mode.

When explicit-orientation-sum CFF hits a duplicate-leading-denominator boundary
inside a UV-rescaled subgraph, the corresponding forest term is routed through
the expanded-4D bridge.  This is a structural boundary of that local forest
term, not a graph-specific fallback.

### Local UV from Expanded 4D Integrands

This path is available for CFF and LTD and is required for LTD cross-section
supergraphs.

The construction:

- expands the relevant UV forest term in four dimensions,
- builds a parsed 4D source for each local term,
- keeps the local numerator in strict four dimensions,
- does not run `gamma_simplify` on the local numerator, because that would
  destroy local cancellations,
- applies the selected 3D representation to the local source,
- selects the same threshold/LU residues as the original integrand through the
  representation-neutral cut catalogue,
- reconstructs the local counterterm in the original event-weight structure.

The source can be disconnected.  Its parsed representation carries enough
information for the 3D generator to build each component consistently while
keeping numerator contractions intact at the GammaLoop expression level.

### Integrated UV Counterterms

Integrated UV counterterms remain representation-independent functions of the
spatial loop momenta.  The denominator-like normalization used to fold them into
the original integral measure is not a CFF or LTD representation and is not
modified by the 3D-representation choice.

## Cross-Section Integrands

Cross-section supergraphs are handled as forward-scattering graphs with
initial-state cut edges externalized before 3D generation.

The local inspect output is required to match between:

- `CFF + local UV from 3D expansions`,
- `CFF + local UV from expanded 4D integrands`,
- `LTD + local UV from expanded 4D integrands`.

The rich comparison checks:

- total inspect weight,
- event-group count and grouping,
- event weights,
- additional weights such as original and threshold-counterterm components.

Known cross-section generation boundaries:

- UV forest unions with disconnected union nodes are not supported.
- Raised threshold counterterms for cross-section supergraphs are not
  supported.

The scalar 3L all-graph test harness exposes both boundaries through per-graph
fallback flags so future work can turn them into hard failures when the missing
features are implemented.

## Amplitude Integrands

Amplitude integrands use the same representation setting and the same
underlying 3D-expression generator.

The inspect tests compare CFF and LTD for amplitude topologies in
explicit-orientation-sum mode, including imported scalar examples with repeated
propagators and higher-power EMR numerators.  CFF behavior remains unchanged
when the representation setting is left at the default.

## 3Drep CLI

The `3Drep` CLI command exposes the 3D-representation crate for diagnostics and
standalone evaluator experiments.

Relevant modes include:

- `build`: build and cache oriented JSON expressions.
- `display`: print graph and expression structure.
- `evaluate`: build/evaluate a cached or explicit oriented expression.
- `test-cff-ltd`: compare CFF, LTD, and diagnostic pure-LTD/proxy modes where
  appropriate.
- `graph-from-signatures`: construct a diagnostic graph from momentum
  signatures; this is not used by GammaLoop production logic.
- `validate`: validate generated representation artifacts.

Relevant `evaluate` artifacts:

- `symbolica_expression_pretty.txt`: formatted `expr.log_print(Some(80))`.
- `symbolica_expression.txt`: canonical Symbolica expression string.
- `symbolica_expression_raw.json`: the complete input archive for Symbolica
  evaluator construction after Spenso numerator treatment.
- `symbolica_expression_raw.rs`: executable standalone Rust script that reads
  the JSON archive, builds the evaluator call named `3d_rep`, and performs one
  representative evaluation.

`3Drep evaluate --standalone-rust-only` writes the expression/debug artifacts
and standalone script, then skips evaluator construction, evaluation, and
manifest writing.

## Testing Surface

The selected default GammaLoop suite is:

```bash
just test_gammaloop
```

It runs with warnings as errors and excludes explicit `slow` and `failing`
classes by default.

The detailed nextest display alias is:

```bash
just test_gammaloop_detailed
```

Slow scalar cross-section 3L coverage is available with:

```bash
just test_gammaloop slow
```

Current coverage includes:

- unit and crate tests in `three-dimensional-reps` for CFF/LTD generation,
  initial-state cut externalization, external trees, repeated channels, lower
  sectors, high-power numerator sectors, and sampling-scale behavior;
- CLI integration tests for `3Drep` build/evaluate/test-cff-ltd workflows;
- amplitude inspect parity tests comparing CFF and LTD;
- cross-section inspect parity tests comparing the three production local-UV
  modes;
- a scalar 3L all-graph slow harness over `GL00` through `GL48`, with baseline,
  quadratic, and targeted quartic numerator probes;
- six scalar cross-section rich-inspect anchors in the default test profile.

As of the current branch state, there is no retained scalar cross-section rich
local-inspect parity failure among the three production modes.  The exhaustive
slow scalar sweep documents per-graph UV/threshold fallbacks in
`scalar_3L_cross_section_all_graphs_test.md`.

## Current Unsupported Cases

- `3d_representation = "LTD"` without
  `generation.explicit_orientation_sum_only = true`.
- Explicit orientation selection or individual-orientation Monte Carlo runtime
  sampling in explicit-orientation-sum mode.
- Non-default `generation.orientation_pattern` in explicit-orientation-sum
  mode.
- LTD local UV counterterms from 3D expansions.
- Cross-section raised threshold counterterms.
- UV forest unions with disconnected union nodes.
- Diagnostic `PureLTD` as a production formula for unsplit repeated
  propagators with numerator dependence.
- Very high CFF energy-degree bounds can be impractical even when supported by
  the algebraic construction.

## Review Notes

The branch-wide review found no production graph-label conditionals or
per-topology sign patches in the 3D-representation production paths.  The
remaining sign corrections are:

- the LTD LU loop-measure parity,
- the expanded-4D CFF source/full-graph sign-exponent correction,
- ordinary surface-canonicalization signs carried by residue selection.

Unchecked threshold-residue selection in the cut CFF helper now reports a
structured invariant error if a threshold residue does not produce exactly one
expression.
