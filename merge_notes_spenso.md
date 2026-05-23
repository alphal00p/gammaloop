# Merge notes: `improve-parsing-and-execution` into `ltd_in_gammaloop_symbolica_update`

This file records non-trivial decisions made while attempting the merge on the
study branch `codex/spenso-merge-study`.

The merge was started with:

```bash
git merge --no-commit --no-ff origin/improve-parsing-and-execution
```

The intent is to keep the current branch's CFF/LTD/3Drep parity and recent CLI
state/runtime work as authoritative, while bringing in the Spenso/Idenso parsing
and execution refactor from `improve-parsing-and-execution`.

## Initial direct conflicts

- `AGENTS.md`
- `Cargo.lock`
- `crates/gammaloop-api/src/python.rs`
- `crates/gammalooprs/Cargo.toml`
- `crates/gammalooprs/src/cff/generation.rs`
- `crates/gammalooprs/src/processes/amplitude.rs`
- `crates/gammalooprs/src/utils/symbols.rs`
- `crates/gammalooprs/src/uv/approx/integrated.rs`
- `crates/gammalooprs/src/uv/approx/mod.rs`
- `crates/gammalooprs/src/uv/forest.rs`
- `crates/gammalooprs/state_map.bin`
- `crates/gammalooprs/tests/test_renormalization.rs`
- `crates/idenso/src/color.rs`
- `crates/idenso/src/gamma.rs`
- `crates/spenso/src/network/library/symbolic.rs`
- `crates/spenso/src/network/parsing.rs`
- `crates/spenso/src/structure/representation.rs`
- `crates/spenso/src/tensors/parametric.rs`

## Running decision log

### Overall merge strategy

Decision: preserve the current branch's GammaLoop UV/CFF/local-inspect
architecture and port the incoming Spenso/Idenso canonicalization and execution
changes into it, instead of taking the incoming older UV control flow.

Reason: the current branch contains the newer keyed `CutCFFIndex` construction,
expanded-4D local UV parity work, 3Drep CLI support, inspect benchmarking,
placeholder syntax, process-aware graph import, and numerical stability
aggregation. Replacing those flows with incoming older versions would likely
regress already validated scalar cross-section parity.

Open risk: several incoming simplification calls may be semantically required
for the new Spenso execution model, but the correct insertion point in the
current UV pipeline is not obvious in every case.

### `crates/gammalooprs/src/cff/generation.rs`

Decision: keep the current branch version for the first compile pass.

Reason: this file contains the current branch's 3Drep/CFF surface-term
normalization work, including the initial-state cut-edge treatment used by the
local-inspect parity path. The incoming conflict mostly reintroduced an older
test/helper region around orientation generation and did not obviously contain
new production logic that should override the current branch.

Discussion point: after compilation, compare the incoming branch's CFF tests and
helpers to decide whether any of them should be reintroduced as tests against the
current implementation.

### `crates/gammalooprs/src/uv/approx/mod.rs`

Decision: keep the current branch version for the first compile pass.

Reason: the current branch has the keyed `CutCFFIndex` / expanded-4D / local-3D
parity implementation. The incoming branch's conflicted section used the older
order-zipped integrated-UV flow. Reintroducing it would risk undoing the scalar
cross-section parity work.

Discussion point: the incoming Spenso/Schoonschip simplification may still need
to be inserted into current UV helper functions. This must be done deliberately,
not by accepting the older incoming UV control flow.

### `crates/gammalooprs/src/uv/approx/integrated.rs`

Decision: keep the current helper-based `integrated_uv_start(...)` flow and do
not inline the incoming older construction.

Reason: the incoming conflict block referenced `graph` and `reduced` variables
from its older local context and divided by the graph denominator inline. The
current branch factors that logic into helpers that are shared with the local UV
parity implementation.

Discussion point: the incoming chain
`collect_metrics().simplify_metrics().simplify_gamma().schoonschip_net::<Aind>().to_dots().normalize_dots()`
is probably important for the new parser/execution branch. It should be audited
and inserted into the current helper flow only if it preserves the current UV
local-limit semantics.

### Idenso and Spenso module restructuring

Decision: accept deletion of the old flat files:

- `crates/idenso/src/color.rs`
- `crates/idenso/src/gamma.rs`
- `crates/spenso/src/network/parsing.rs`

Reason: the incoming branch replaces these with structured module directories.
Keeping the old files would duplicate modules and conflict with the incoming
public API.

Discussion point: current-branch fixes to Symbolica print closure signatures and
GammaLoop-specific parser assumptions must be ported into the new module tree
where needed.

### `crates/gammalooprs/src/utils/symbols.rs`

Decision: combine current branch's Symbolica 3-argument print closures and
concrete Lorentz-index normalization with incoming Spenso tensor/scalar tags.

Reason: the current branch is compatible with the newer Symbolica print callback
signature and handles both `cind(i)` and `Minkowski(4,i)` concrete indices. The
incoming branch adds tags that the new Spenso parser uses to identify rank-1
tensors, scalar functions, and tensor symbols.

Discussion point: this is one of the most important semantic resolutions. The
merged symbol set should be checked against the actual parser expectations,
especially for `Q`, `Q3`, `OSE`, external/loop momenta, spinors, and polarization
vectors.

### Spenso symbol definitions

Decision: for `crates/spenso/src/network/library/symbolic.rs` and
`crates/spenso/src/structure/representation.rs`, keep current branch's
three-argument Symbolica print callbacks and combine them with incoming tags such
as `Linear` and `SPENSO_TAG.upper`.

Reason: these are compatibility and metadata changes that should both survive.

Discussion point: whether every representation should still receive the `upper`
tag is worth checking against the incoming parser assumptions; the first merge
keeps the current branch's behavior.

### Symbolica API compatibility in incoming Spenso/Idenso files

Decision: update incoming two-argument `print = |a, opt|` callbacks to the
current Symbolica three-argument `print = |a, opt, _state|` form in the new
Spenso/Idenso module tree, and update `collect_symbol::<i16>(symbol, None, None)`
to the current one-argument API.

Reason: these are API-compatibility fixes required for the merged tree to
compile against the Symbolica revision already used by the current branch.

Discussion point: this is mechanically necessary, but it confirms that the
incoming branch had not been rebased over the current Symbolica API changes.

Follow-up compile fixes in the same category:

- Changed incoming diagnostic tests from `Atom::parse(wrap_input!(...), settings)`
  to `Atom::parse_with_default_namespace(wrap_input!(...), settings)`.
- Updated old `idenso::metric::MetricSimplifier` imports to
  `idenso::shorthands::metric::MetricSimplifier`.

## Compile and test status

After resolving the conflicts, Symbolica/Idenso/Spenso API drift, and the
current branch's local-4D UV projection refactor, the merged study branch
compiles with:

```bash
cargo check --workspace --all-targets --locked --profile dev-optim
```

Targeted scalar local-inspect parity remains green:

```bash
env EXTRA_MACOS_LIBS_FOR_GNU_GCC=T RUST_BACKTRACE=0 RUST_MIN_STACK=33554432 \
  RUSTFLAGS=-L/opt/local/lib/libgcc INSTA_FORCE_PASS=1 \
  cargo nextest run -p gammaloop-integration-tests --test test_runs \
  --cargo-profile dev-optim --run-ignored all --ignore-default-filter \
  -E 'test(/scalar_3l_cross_section_inspects::slow/)' --no-capture --retries 0
```

Result: 143 tests run, 143 passed. The accepted snapshot drift was restricted to
last-displayed-digit changes in scalar real/imaginary values.

The remaining selected-suite failures after the merge were:

- example-card loading discovered historical benchmark output TOMLs under
  `examples/cli/bench/OLD_outputs/...`;
- the AA AA symjit local-inspect backend hit a SymJIT direct-translator abort on
  `_func_power_`;
- the `utils::test_utils::normalization` inline snapshot needed to reflect that
  `Q` now carries Spenso rank-one/tensor tags;
- ignored `large_spenso_actual` diagnostics were pulled into
  `just test_gammaloop` by `--run-ignored all` and initially left scalar
  expressions as non-scalar tensors when parsing tag-stripped inputs.

The first three have targeted passing coverage. The large Spenso diagnostic
group now also passes after using representation-slot based tensor recognition
for those tag-stripped diagnostic expressions.

## High-risk auto-merged areas still needing review

These files did not all have direct conflict markers, but the merge result should
be reviewed before considering this branch semantically correct:

- `crates/gammalooprs/src/integrands/process/evaluators.rs`
- `crates/gammalooprs/src/model/mod.rs`
- `crates/gammalooprs/src/numerator/*`
- `crates/gammaloop-api/src/state.rs`
- `crates/gammaloop-api/src/tracing.rs`
- `crates/gammaloop-tracing-filter/*`
- `crates/spenso/src/network/*`
- `crates/spenso/src/shadowing/*`
- `crates/idenso/src/color/*`
- `crates/idenso/src/dirac/*`
- `flake.nix`

The compile result only proves that the merged API surface is coherent. It does
not prove that the new Spenso execution path preserves the current branch's
CFF/LTD/3Drep local-inspect parity or existing UV-renormalization behavior.

### `crates/gammalooprs/src/uv/forest.rs`

Decision: combine both sides by keeping the current color simplification result
and preserving incoming debug output before simplification.

Discussion point: the incoming branch comments suggest dimensional replacement
questions around color simplification. I left this as diagnostic-only behavior
for the compile pass.

### `crates/gammalooprs/tests/test_renormalization.rs`

Decision: keep the current branch's cached/uncached exact-expression comparison
for the compile pass.

Reason: the incoming branch changed this test to compare after a normalization
pipeline, which is likely reasonable for the new simplifier but changes the test
oracle. I did not want to silently weaken or redefine the invariant while doing
a merge-only compile pass.

Discussion point: this deserves an explicit decision. If the new Spenso
canonicalization intentionally changes expression presentation while preserving
mathematical equivalence, the test should probably compare a canonicalized form.
If the current branch's exact cache-equivalence invariant is still important, the
implementation must keep exact equality between cached and uncached paths before
presentation normalization.

### `crates/gammalooprs/state_map.bin`

Decision: temporarily take the incoming binary for the compile pass.

Reason: the incoming branch introduces a larger Spenso/Idenso symbol/tag set and
the binary is required to resolve the merge. This is not a principled final
resolution.

Discussion point: the correct final resolution is to regenerate `state_map.bin`
from the final merged symbol set, because both branches changed symbol/state
assumptions.

### `Cargo.lock`

Decision: regenerate `Cargo.lock` with `cargo generate-lockfile` after combining
the merged manifests, rather than choosing either side of the lockfile conflict.

Reason: both branches add workspace dependencies. The resulting lockfile should
reflect the merged `Cargo.toml` graph.

Discussion point: this updates the lockfile mechanically but should still be
reviewed for unintended dependency version shifts, because `cargo
generate-lockfile` resolves the full workspace graph from the merged manifests.

### Local 4D UV projection boundary

Decision: refactor the expanded-4D local UV path so the forest recursion remains
in pure 4D `Atom` space and the 3D representation map happens at the projection
boundary.

Reason: this follows the target described in `TARGET_UV_STRUCTURE.md`. The new
`Local4DApproximation` implements the same `ApproximationKernel<UVCtx<'_>>`
shape as the integrated UV flow: `start` builds the integrated-style source,
`t` applies the local 4D rescaling and projection, and the caller maps the final
local atom into top-level `CutCFFIndex` keyed 3D representation terms once. This
prevents factorized subgraph residue details from leaking across the UV
orchestrator boundary.

Discussion point: this is the intended long-term structure, but the remaining
semantic validation is the full selected GammaLoop suite and another unforced
slow scalar cross-section sweep after all merge fallout is fixed.

Follow-up: restored the incoming branch's tensor simplification chain in the
shared UV source construction used by both integrated UV and local expanded-4D
UV:

```text
collect_metrics()
simplify_metrics()
simplify_gamma()
schoonschip_net::<Aind>()
to_dots()
normalize_dots()
```

This lives before the integrated and local-4D paths diverge, so the local-4D
projection boundary remains intact while both paths see the same simplified 4D
source numerator.

### SymJIT direct translator

Decision: add local `symjit.toml` files with `options.direct = false`.

Reason: the merged selected suite exposed an upstream SymJIT direct-translator
regression: the direct translator emits a call to the `power` builtin but does
not install the `_func_power_` label in the generated function table. The
indirect translator registers that builtin and the AA AA symjit local-inspect
precision/backend test passes with `direct = false`. This is a backend
configuration change, not a GammaLoop sign or normalization workaround.

Discussion point: this should be reported upstream if not already known. A
confirmed SymJIT translator bug should not be hidden behind expression-level
rewrites in GammaLoop.

### Tag-stripped Spenso diagnostic parsing

Decision: set `StrictTensorFilter::ContainsReps` in the `large_spenso_actual`
diagnostic parser settings.

Reason: those diagnostics deliberately strip serialized `::{...}::` Symbolica
tag annotations from large saved expressions before parsing. With the merged
Spenso parser's default `Tagged` filter, untagged heads such as `Q(...)`,
`ϵ(...)`, and `ϵbar(...)` are treated as scalar even though they expose
representation slots. Metrics and gamma tensors still parse tensorially, so
mathematically scalar contractions are left as rank-two or higher tensor
results. `ContainsReps` is the principled mode for tag-stripped historical
GammaLoop expressions: any ordinary head with direct representation syntax is a
tensor leaf.

Discussion point: production GammaLoop numerator parsing currently mostly relies
on tagged symbols or explicit parser-owned Spenso heads. It may still be worth
centralizing a GammaLoop numerator parse setting that uses the same
representation-slot filter for user-provided or tag-stripped numerator strings,
but I kept this change scoped to the failing diagnostics for the first merge
stabilization pass.

### Example-card discovery

Decision: skip both `outputs` and `OLD_outputs` under `examples/cli/bench` when
loading example run cards.

Reason: these folders contain generated benchmark outputs and saved state TOMLs,
not top-level example cards. Treating them as cards makes the loadability test
fail for data that is not intended to be run-card input.

### AA AA 3Drep CFF/LTD assembly diagnostic

Decision: make `3Drep test-cff-ltd` use the existing iterative/componentized
build strategy when the selected evaluator backend is `assembly` and the user did
not already request another strategy.

Reason: the failing
`cli_aa_aa_box_and_double_box_multibackend_threedrep_comparisons` case is the
AA->AA Box CFF/LTD diagnostic in `Double + assembly`. Fresh eager and SymJIT
diagnostics agree with LTD/PureLTD at roughly machine precision, and
`Double + assembly --iterative` does too. The mismatch is therefore not a CFF/LTD
construction error: it is specific to lowering one very large cancellation-heavy
CFF orientation sum into a single inline-assembly atom. The iterative strategy
uses the same generated representation terms and still exercises the assembly
backend, but leaves the explicit numerator/orientation components visible at the
evaluator boundary before summing them.

Discussion point: this is intentionally scoped to the `test-cff-ltd`
representation-equivalence diagnostic. It should not be generalized into a
per-graph sign/normalization change or into production CFF/LTD generation.
