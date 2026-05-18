# Multi-Way Comparison Debug Status

This note tracks the scalar cross-section rich local inspect parity effort
across:

1. CFF with local UV counterterms from 3D expansions,
2. CFF with local UV counterterms from expanded 4D integrands,
3. LTD with local UV counterterms from expanded 4D integrands.

The invariant being protected is event-level parity of `Original`, threshold
counterterms, full multiplicative factors, cut ids, event grouping, and final
weights. Overall common phases are less important than relative parity between
the three modes, but graph-specific sign modifiers and loop-number fudge
factors are not allowed.

## Current Checkpoint

The slow scalar cross-section local-inspect matrix is the current production
anchor. At this revision it passes as a true no-force snapshot run after
accepting deterministic final-weight snapshot updates that were exposed by the
principled sign/coordinate fixes. The rich parity assertions run before the
snapshot assertion and still compare the three independent GammaLoop production
modes, not a shared 3Drep diagnostic expression.

The current implementation changes since the previous pushed checkpoint are:

- Forward-scattering initial-state cut carriers are handled as external
  half-edge coordinates. If the generated selected LU `E`-surface and the
  omitted positive-energy Cutkosky edge have identical derivatives with respect
  to several initial-state cut carriers, these candidates are equivalent local
  charts for the same pole. The local-coordinate builder now selects a
  canonical LMB external-coordinate representative and lets the resolved
  coordinate determinant supply the physical Jacobian.
- Direct LTD original/root LU residues no longer include the auxiliary CFF
  full-graph projection bridge. Those terms are extracted from the requested
  LTD production expression; the CFF projection is only metadata for identifying
  and canonicalizing LU surfaces. The global LTD loop-measure bridge is still
  applied by the embedding code.
- The ordinary simple direct-original bridge is now expressed directly as the
  LU surface-family orientation times the local-series denominator-coordinate
  signs and the full-source same-routing repeated-channel source bridge.
  Repeated selected residues still use their own repeated-pole derivative and
  edge-support bridges.
- Scalar inspect snapshots that drifted only in the final benchmark value were
  accepted after the rich multi-way checks passed. Some old references were
  subnormal placeholders from previous behavior; the new references are the
  deterministic values produced by the current production construction.

The current working tree has no graph-name branches, no loop-count sign
exceptions, no production use of `graph_from_signature`, and no temporary
`GAMMALOOP_TRACE_LTD_LU_SIGNS` diagnostics in Rust code.

## Validation Matrix

Already run on this revision:

```text
env -u RUSTFLAGS EXTRA_MACOS_LIBS_FOR_GNU_GCC=T RUST_BACKTRACE=0 \
  RUST_MIN_STACK=33554432 \
  cargo nextest run -p gammaloop-integration-tests \
  --test test_evaluation_api --cargo-profile dev-optim \
  --run-ignored all --ignore-default-filter --no-capture --retries 0

12 tests run: 12 passed
```

```text
env -u RUSTFLAGS EXTRA_MACOS_LIBS_FOR_GNU_GCC=T RUST_BACKTRACE=0 \
  RUST_MIN_STACK=33554432 \
  cargo nextest run -p gammaloop-integration-tests --test test_runs \
  --cargo-profile dev-optim --run-ignored all --ignore-default-filter \
  -E 'test(/differential::lu/) | test(/spin_sums/) | test(/ltd_generated_forward_cross_section/)' \
  --no-capture --retries 0

9 tests run: 9 passed
```

```text
env -u RUSTFLAGS EXTRA_MACOS_LIBS_FOR_GNU_GCC=T RUST_BACKTRACE=0 \
  RUST_MIN_STACK=33554432 \
  cargo nextest run -p gammaloop-integration-tests --test test_runs \
  --cargo-profile dev-optim --run-ignored all --ignore-default-filter \
  -E 'test(/scalar_3l_cross_section_inspects::default_scalar_3l_cross_section_inspects/)' \
  --no-capture --retries 0

6 tests run: 6 passed
```

```text
env -u RUSTFLAGS EXTRA_MACOS_LIBS_FOR_GNU_GCC=T RUST_BACKTRACE=0 \
  RUST_MIN_STACK=33554432 INSTA_FORCE_PASS=1 \
  cargo nextest run -p gammaloop-integration-tests --test test_runs \
  --cargo-profile dev-optim --run-ignored all --ignore-default-filter \
  -E 'test(/scalar_3l_cross_section_inspects::slow/)' \
  --no-capture --retries 0

143 tests run: 143 passed
```

The forced run generated scalar inspect `.snap.new` files. They were inspected:
the diffs were limited to `assertion_line` metadata and the final benchmark
value (`re`, plus tiny imaginary-rounding drift in a few cases). The rich
multi-way parity assertions had passed for every generated snapshot. Those
snapshot updates were accepted, then the real no-force anchor was run:

```text
env -u RUSTFLAGS EXTRA_MACOS_LIBS_FOR_GNU_GCC=T RUST_BACKTRACE=0 \
  RUST_MIN_STACK=33554432 \
  cargo nextest run -p gammaloop-integration-tests --test test_runs \
  --cargo-profile dev-optim --run-ignored all --ignore-default-filter \
  -E 'test(/scalar_3l_cross_section_inspects::slow/)' \
  --no-capture --retries 0

143 tests run: 143 passed (2 slow), 124 skipped
```

## Sensitive Construction Points

- `Graph::active_three_d_repeated_channels()` is the single repeated-channel
  descriptor used for global 3D sign excess, residue bridges, source duplicate
  bridges, and raised-edge support normalization.
- `CrossSectionGraph::residue_selector_for_raised_cut_group` is the central
  cross-section entry point for LU residue planning. All simple, raised,
  repeated, threshold, and expanded-4D UV residue selectors must flow through
  this graph/LMB-derived data instead of separate per-case sign algorithms.
- `CrossSectionGraph::lu_cut_edge_sets_with_cutkosky_signs` maps raised LU
  groups to positive-energy Cutkosky edge support and signs.
- `CrossSectionGraph::simple_ltd_lu_cut_local_coordinate_signs` computes the
  generated selected `E`-surface sign, local-series denominator signs, and
  LU surface-family bridge from the resolved coordinate determinant. It now
  treats equivalent initial-state cut half-edge coordinates as the same local
  chart.
- `ResidueSelector` stores selected signs, local-series signs, ordinary residue
  prefactors, direct-original prefactors, and threshold prefactors. Direct
  original/root residues intentionally omit the auxiliary CFF projection bridge
  in LTD production mode.
- `cff/expression.rs::ltd_lu_local_series_coefficients_from_parametric_atom`
  consumes the LTD LU local-series coordinates and signs. The expressions
  passed to it for LTD production are LTD expressions.
- `uv/approx/mod.rs` and `uv/approx/expanded_4d.rs` dispatch threshold residue
  selection by representation and must remain LTD-expression-facing for LTD
  output.
- `three-dimensional-reps::generate_cff_ltd_comparison_expression` is the only
  diagnostic CFF/LTD comparison entry point. It is intentionally not the
  GammaLoop production CFF local-UV generator.

## Remaining Plan

1. Commit and push this green scalar checkpoint before touching the remaining
   full-suite failures.
2. Re-run the selected full GammaLoop suite with the macOS build-link
   environment:

   ```text
   env -u RUSTFLAGS EXTRA_MACOS_LIBS_FOR_GNU_GCC=T just test_gammaloop
   ```

   Compare the resulting failure buckets against the previous `1101/1131`
   status. The coordinate-blocked API, differential LU, spin-sum, and generated
   forward inspect tests are expected to be resolved by the current changes.
3. If failures remain in the divergent/amplitude-like scalar bubble local or
   integrated UV bucket, handle them separately from the scalar cross-section
   LU bridge:
   - first confirm the original bubble graph evaluates consistently in CFF and
     LTD at the same sample;
   - then compare CFF 3D local UV, CFF expanded-4D local UV, and LTD
     expanded-4D local UV before event assembly;
   - classify the mismatch as source-basis projection, repeated/contact
     collapse, local/integrated UV normalization, or numerator localization;
   - fix only the shared mathematical source, then rerun the full slow scalar
     sweep immediately.
4. Revisit any 3Drep diagnostic failures only after production GammaLoop
   failures are classified. Fix them through the centralized diagnostic entry
   point or shared graph/LMB coordinate construction, not by redirecting
   production CFF local-UV generation.
5. Final acceptance after all buckets are green:
   - `cargo fmt`
   - `env -u RUSTFLAGS EXTRA_MACOS_LIBS_FOR_GNU_GCC=T cargo check`
   - `env -u RUSTFLAGS EXTRA_MACOS_LIBS_FOR_GNU_GCC=T just test_gammaloop`
   - full slow scalar cross-section sweep without `INSTA_FORCE_PASS`
   - broad 3Drep integration sweep
   - no generated `.snap.new` files unless deliberately accepted
   - commit, push, and update the PR body with the implemented changes and
     validation matrix.
