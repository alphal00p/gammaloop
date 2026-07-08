# Renormalization forest audit

## Reproduction

Baseline: `main` at `4fdbf430`. Current implementation: `69b3078b` (the
parent of the empty investigation change).

The baseline was checked out with `jj new main`; the current change was restored
with `jj edit xtopnouw`. Both revisions used the same serial command:

```bash
GL_DISPLAY_FILTER=off \
GL_LOGFILE_FILTER='off,gammalooprs::uv::forest[{generation,uv,graph,term,!summary}]=debug,gammalooprs::uv::approx::integrated[{uv,integrated,vakint,trace,input}]=debug,gammalooprs::uv::approx::integrated[{uv,integrated,vakint,trace,raw}]=debug' \
GL_TEST_LOG_DIR=<log-dir> \
cargo nextest run -p gammalooprs --test test_renormalization \
  --cargo-profile dev-optim --test-threads 1 --no-fail-fast
```

`XDG_STATE_HOME` was also placed under `<log-dir>` because the sandbox cannot
write nextest's default run store. Raw captures are in:

- `/tmp/gammaloop-renormalization-main-4fdbf430`
- `/tmp/gammaloop-renormalization-current-69b3078b`
- `/tmp/gammaloop-renormalization-pole-sign`
- `/tmp/gammaloop-rqft-integrated-only`

The records below omit timestamps, retries, backtraces, timings, ANSI styling,
dummy-index differences, and finite terms that cannot affect the requested pole
part.

## Test results

| Test | `main` | Before correction | Pole-sign correction | Integrated-only RQFT recursion |
| --- | --- | --- | --- | --- |
| `scalar_pole_part` | pass (no assertion) | pass (no assertion) | pass (no assertion) | pass (no assertion) |
| `finite_part_quark_lo` | pass | fail | pass | pass |
| `finite_part_ghost_2loop` | pass | fail on `d1` | fail on `d2` | pass with native GammaLoop snapshot |
| `finit_part_ghlo` | pass | fail | pass | pass |
| `failing::finite_part_ghost_3loop` | fail on `d1` | fail on `d1` | fail on `d1` | same failure as `main` |

The two-loop `d1` snapshot is

```text
main:              -3i/16 eps^-2 + 5i/32 eps^-1
before correction: -9i/16 eps^-2 - 5i/32 eps^-1
```

The three-loop snapshot expects

```text
-14 eps^-3 + 65/3 eps^-2 + (-91/4 - 7*pi^2/3) eps^-1
```

but `main` already produces

```text
-14 eps^-3 + 23/3 eps^-2 + (-203/36 - 7*pi^2/3) eps^-1
```

and the implementation before the pole-sign correction produces

```text
-14 eps^-3 + 94 eps^-2 + (1169/12 - 7*pi^2/3) eps^-1.
```

The three-loop baseline failure is therefore independent of the current sign
change. The pole-sign correction returns it to the `main` value exactly.

## Snapshot normalization

Let `T1` and `T2` denote the RQFT targets written in the comments above graphs
`d1` and `d2`. The snapshot and helper changes have different effects:

| Graph | `main` helper | `main` snapshot | Current helper | Current snapshot | Expected raw result |
| --- | --- | --- | --- | --- | --- |
| `d1` | `-atom` | `T1` | `atom` | `T1` | changes from `-T1` to `T1` |
| `d2` | `-atom` | `-T2` | `atom` | `T2` | remains `T2` |

The `d2` snapshot was changed relative to `main`, from
`+1/16 eps^-2 - 1/32 eps^-1` to
`-1/16 eps^-2 + 1/32 eps^-1`. That is exactly the global sign needed to
compensate for removing `-atom` from `align_to_rqft`; it does not change the
underlying expected renormalization expression. By contrast, leaving the `d1`
snapshot unchanged while removing the helper minus intentionally flips its raw
expectation.

## Two-loop forest

For graph `d1`, factor out
`H = i*GC_10^4*(q0.q4)*CA^2`. The terminal pole contributions are:

| Forest leaf | `main / H` | Before correction `/ H` |
| --- | --- | --- |
| `140 -> 0` | `-3/16 eps^-2 - 5/32 eps^-1` | same |
| `140 -> 132 -> 0` | `+3/16 eps^-2` | `-3/16 eps^-2` |
| `140 -> DM -> 0` | `+3/16 eps^-2` | `-3/16 eps^-2` |

The primitive leaf is unchanged. Both leaves containing one nested forest
operation differ in sign from `main`. Their double poles consequently add to the
root instead of cancelling part of it in the current snapshot.

The one-loop traces show the same primitive expressions on both revisions. The
one-loop failures expose the convention change from removing the overall minus
from `align_to_rqft`.

## Three-loop forest

For graph `d1`, factor out
`G = GC_10^2*GC_12^2*(q0.q2)*CA^3`. All eight reached terminal leaves are:

| Forest leaf | `main / G` | Before correction `/ G` |
| --- | --- | --- |
| `4GC -> 0` | `-91/6 eps^-2 - 1295/36 eps^-1` | same |
| `4GC -> 4C8 -> 0` | `+35/2 eps^-2 + 91/3 eps^-1` | negative of `main` |
| `4GC -> 3Oa -> 0` | `+28/3 eps^-2 + 52/9 eps^-1` | negative of `main` |
| `4GC -> zw -> 0` | `-14 eps^-3 + 44/3 eps^-2 + (52/9 - 7*pi^2/3) eps^-1` | negative of `main` |
| `4GC -> 3zk -> 0` | `+28/3 eps^-2 + 4 eps^-1` | negative of `main` |
| `4GC -> 4C8 -> 3Oa -> 0` | `-28/3 eps^-2 - 52/9 eps^-1` | same |
| `4GC -> 4C8 -> zw -> 0` | `-28/3 eps^-2 - 52/9 eps^-1` | same |
| `4GC -> 4C8 -> 3zk -> 0` | `-28/3 eps^-2 - 4 eps^-1` | same |

The pattern is exact: depth-one terms equal `main`, depth-two terms are negated,
and depth-three terms equal `main` again. This is a sign-composition issue, not a
Vakint coefficient error.

## Audit

The intended convention distinguishes the actual integrated counterterm from a
pole-only output projection:

1. Every operation is signed: `local_4d::uv_limit` constructs the next local
   counterterm as `-T(...)` in
   `crates/gammalooprs/src/uv/approx/local_4d.rs`. Nested operations must obtain
   their alternating signs from this recursion.
2. The MS-bar integrated counterterm is negative relative to that signed local
   term: `-finite(integral(local))` removes the finite part so the local
   counterterm retains the poles.
3. `pole_part=true` asks for the pole output itself. That branch must preserve
   the sign already supplied by the `T` operations rather than adding the
   integrated-counterterm minus.

No extra terminal parity should be applied. Doing so would move the operation
sign out of `T` and double-sign correctly composed terms. `SimpleApprox.sign` is
currently unused metadata, and `OperationNode.key.op_count()` should not be
applied by either terminal aggregator.

The correction is consequently limited to `Integrated::run`:

```text
pole_part=true:   pole(integral(local))
pole_part=false: -finite(integral(local))
```

This preserves operation parity in pole mode without changing the actual
MS-bar addback. It advances the suite from one passing test to three; graph `d1`
also reaches and passes its legacy/hedge equality and forest-size checks.

With the pole-sign correction but before selecting the RQFT recursion input,
graph `d2` gives

```text
expected: -1i/16 eps^-2 + 1i/32 eps^-1
actual:   -1i/16 eps^-2 - 1i/32 eps^-1
```

Its three terminal leaves, after factoring out
`i*GC_10^4*(q0.q4)*CA^2`, are:

| Forest leaf | Pole contribution |
| --- | --- |
| `140 -> 0` | `-1/16 eps^-2 - 1/32 eps^-1` |
| `140 -> 132 -> 0` | `+1/16 eps^-2` |
| `140 -> oW -> 0` | `-1/16 eps^-2` |

The raw RQFT signs in that snapshot are not the native GammaLoop convention for
this graph. GammaLoop's `+i` ghost propagator and RQFT's `-i` ghost propagator
give a full-graph factor `GammaLoop/RQFT = -1` because `d2` has three internal
ghost lines; the vertex convention changes cancel in the full product. The
native GammaLoop reference is therefore direct
`-1/16 eps^-2 - 1/32 eps^-1` and two nested terms of `+1/16 eps^-2`. The direct
and `132` leaves already match it; only `oW` is wrong. The explicit RQFT
per-term convention map now lives with the `d2` snapshot in
[`test_renormalization.rs`](../crates/gammalooprs/tests/test_renormalization.rs).

A controlled parent decomposition identifies the wrong recursion operand.
Feeding only the dependent integrated pole into the outer operation gives
`+1/16 eps^-2` for both nested leaves, while feeding only the dependent local
term gives zero for both. In RQFT pole mode the local term must be discarded
after integration, so both orchestrators now use only the integrated pole for a
non-root dependency. Ordinary MSbar generation continues to recurse through
`local + integrated`, where `integrated = -finite(integral(local))`.

After this change, the three `d2` leaves are the native GammaLoop reference:

| Forest leaf | Pole contribution |
| --- | --- |
| `140 -> 0` | `-1/16 eps^-2 - 1/32 eps^-1` |
| `140 -> 132 -> 0` | `+1/16 eps^-2` |
| `140 -> oW -> 0` | `+1/16 eps^-2` |

Their sum is `+1/16 eps^-2 - 1/32 eps^-1`. The inline `d2` forest breakdown now
records both this native GammaLoop value and the graph's
`GammaLoop / RQFT = -1` convention factor. With that expectation,
`finite_part_ghost_2loop` passes completely: all six graphs, every
legacy/hedge-poset equality, and every forest-size assertion pass.

There is also a separate hedge-poset composition risk. Its symbolic node label
multiplies each commuting leaf operation with that leaf's own dependency
frontier. `compute_4d_for_node`, however, loads only the first leaf's frontier
and sequentially mutates that value for every later leaf. A legacy/hedge
mismatch confined to multi-leaf Foata levels after the pole-sign correction
would point to this order-dependent computation.

## Coverage gaps

- `finite_part_ghost_2loop` asserts each legacy snapshot before calling
  `assert_new_paths_match_legacy`. Graph `d1` now checks both paths, but the `d2`
  snapshot failure prevents its hedge-poset comparison and all later graphs.
- `scalar_pole_part` contains no assertion.
- The three-loop settings use the default legacy orchestrator, so that test does
  not exercise the hedge-poset path.
- The three-loop test aborts at graph `d1`, so later graph checks are not reached.
  For indices 6 through 77, 72 name assertions also read
  `amp.graphs[4].graph.name` instead of the graph just evaluated. The first
  guaranteed mismatch compares that `d5` name with `""`; the accompanying
  repeated `gs^4` expression snapshots are placeholders in a `gs^6` corpus.

The implementation changes separate the pole output sign from the negative
finite-part MSbar addback and discard local dependencies when recursing in RQFT
pole mode. The `d2` inline snapshot records the resulting native GammaLoop
convention.
