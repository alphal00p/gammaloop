# gammaLoop Architecture Ideas

This document proposes architecture improvements beyond the current implementation.

## Proposed Improvements

### P1: Introduce a dedicated application service layer
- Problem: clap command types currently execute domain logic directly and mutate `State` in-place.
- Improvement: add a service boundary (for example `AppService`) that owns use-cases (`generate`, `integrate`, `inspect`, `save`). Keep commands as parsing/DTO only.
- Impact: reduces coupling, improves testability, and makes Python/CLI behavior easier to keep equivalent.

### P1: Replace manual cache-ID protocol with typed cache sessions
- Problem: `ProcessIntegrandImpl` cache semantics rely on manual counter updates and conventions spread across many methods.
- Improvement: introduce a `CacheSession` / `ExternalConfigEpoch` abstraction that centralizes increment/revert semantics and can be validated at compile-time boundaries.
- Impact: fewer cache-invalidations bugs and cleaner evaluation code.

### P1: Unify amplitude/cross-section integrand pipeline shape
- Problem: both paths implement similar lifecycle steps (`preprocess`, `build`, `warm_up`, `compile`, `save/load`) with duplicate orchestration patterns.
- Improvement: extract shared pipeline traits/utilities for process and integrand lifecycle while retaining domain-specific internals.
- Impact: lower maintenance cost and more predictable feature rollout across both generation types.

### P2: Harden unsupported and partially-implemented runtime paths
- Problem: several user-facing paths still rely on `todo!()` or non-implemented branches (e.g. batch command, some event-processing and cross-section standalone behavior).
- Improvement: convert unsupported paths to explicit typed errors and feature flags; track with test coverage for error contracts.
- Impact: avoids panic-style failures in production workflows.

### P2: Remove panic-oriented error paths in sampling/inspect
- Problem: some paths still panic when settings are invalid instead of returning structured errors.
- Improvement: propagate `Result` consistently through inspect/integrand parameterization and command handlers.
- Impact: safer automation and clearer diagnostics in REPL/run cards.

### P3: Process-level orchestration parallelism
- Problem: outer loops over processes/integrands are mostly sequential even when independent.
- Improvement: parallelize process-level preprocessing/build/compile where deterministic ordering is not required.
- Impact: better throughput for multi-process workloads.

### P3: Improve observability for long-running pipelines
- Problem: tracing is already strong, but instrumentation is uneven across some high-cost transitions.
- Improvement: define a small set of canonical spans and metrics for `load -> preprocess -> build -> warm_up -> evaluate/integrate` and emit stable fields.
- Impact: easier profiling and regression triage.

## Suggested Implementation Order
1. Done: Fix persistence format mismatch and introduce versioned manifest.
2. Extract `AppService` use-case boundary from command handlers.
3. Refactor cache handling into typed session semantics.
4. Consolidate amplitude/cross-section lifecycle scaffolding.
5. Replace remaining `todo!()`/panic paths in user-facing flows with typed errors.
