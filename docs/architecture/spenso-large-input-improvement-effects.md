# Spenso Large Input Improvement Effects

This note summarizes the effects of the local improvements stacked on top of
`origin/improve-parsing-and-execution` for `spenso_eval_input_0.txt`. The large
input is still too slow for a broad full-expression matrix, so the measurements
below use term 18 of the 46-way sum, as a representative heavy branch.

Unless stated otherwise, commands used:

```bash
SPENSO_SUM46_START=18 SPENSO_SUM46_LIMIT=1 \
cargo test --profile dev-optim --package gammalooprs \
  --test large_spenso_actual \
  spenso_eval_input_0_sum46_individual_term_execute \
  -- --ignored --nocapture
```

Timings are `dev-optim` diagnostics, not release benchmarks. They are useful for
relative behavior and failure mode identification, not final throughput claims.

## Term 18 Baseline Shape

The selected branch has:

- expression bytes after replacement: `5,220,365`
- selected term bytes: `5,228,463`
- `δ` occurrences: `319`
- `Q3` occurrences: `453`
- parsed graph: `2608` nodes, `4071` graph edges
- tensors: `772` rank-1 param tensors
- scalars: `649`, total scalar payload about `5.18 MB`
- library leaves: `95`

The largest product shape is:

```text
scalar * sum * sum * g * g * g * gamma^6
```

That is important: the remaining problem is not raw parsing. It is late tensor
sums interacting with small tensor products and scalar atom growth.

## Individual Improvements

### Sparse Shadow Tensors

Change: `lxoxtyzm`, controlled by `SPENSO_SPARSE_SHADOW_TENSORS=1`.

Effect on term 18 in `NETWORK_ONLY` mode:

| Mode | Parse | Param entries | Dense tensors | Sparse tensors | Nonzero entries |
| --- | ---: | ---: | ---: | ---: | ---: |
| dense shadows | `352.6 ms` | `3088` | `772` | `0` | `1678` |
| sparse shadows | `359.2 ms` | `1678` | `0` | `772` | `1678` |

This is a clean storage win: stored param entries drop by about `46%`, while the
graph shape does not change. It also exposes the real sparse structure of the
`δ` and `Q3` tensors: `δ` has one nonzero entry, and `Q3` has three.

Separate execution effect: not enough by itself. With sparse shadows but without
lazy tensor sums, term 18 currently fails quickly:

```text
parallel min-result-rank actual execution failed:
fused tensor sum contract requires at least one shared slot
```

That failure is from a fast path declining incorrectly by returning an error for
a no-shared-slot pair.

### One-Hot Selector Contraction

Change: `mqwtplvp`.

This recognizes rank-1 one-nonzero tensors and directly selects the matching
component from the target tensor. Term 18 contains `319` `δ` tensors, so the
applicability is high.

Expected effect:

- `δ(i)` contractions become slice/select operations instead of generic sparse
  tensor contractions.
- The improvement is local and cheap; it does not change parse shape.
- It synergizes strongly with sparse shadows, because sparse storage makes the
  one-hot structure explicit.

Measured separately: no standalone timing is available from the current knobs,
because this path is always active once the code is present. It also does not
solve late `sum * tensor` materialization by itself.

### Fast Sparse Tensor Sum Materialization

Change: `pnznryzy`.

This adds `FastTensorSum` for sparse symbolic `ParamTensor`s. Instead of adding
whole tensors entry-by-entry in a dense loop, it groups sparse entries by output
flat index and only calls `Atom::add_many` for entries that collide.

Expected effect:

- Helps when a lazy sum must be materialized.
- Avoids adding zeros and avoids repeated intermediate atom additions.
- It is most useful after sparse shadows expose the real nonzero support.

Measured separately: parse/network-only shape is unchanged. In full term 18 this
is hit at late materialization points, but term 18 still does not complete within
the current timeout because later pair selection and result materialization
dominate.

### Lazy Scaled Tensor Terms

Change: `wxytorxt`.

This adds `TensorTerm` and `TensorTermSum`, representing `scalar * tensor`
without multiplying the scalar into every tensor entry immediately.

Expected effect:

- Keeps common scalar factors outside tensor entries.
- Prevents scalar payload replication across many sparse tensor entries.
- Enables lazy tensor sums to carry scaled terms forward.

Observed effect with sparse shadows and lazy sums enabled:

```bash
SPENSO_SPARSE_SHADOW_TENSORS=1 \
SPENSO_NETWORK_LAZY_TENSOR_SUMS=1 \
SPENSO_NETWORK_PROFILE=1 \
SPENSO_SUM46_START=18 SPENSO_SUM46_LIMIT=1 ...
```

The executor reaches much later stages and produces many lazy tensor sums, but
term 18 still exceeded `60s`. The late profile showed expensive result-rank pair
selection:

```text
product.result_rank_slow_pair elapsed_ms=28145
product.result_rank_slow_pair elapsed_ms=20104
```

So this moves the bottleneck forward, but does not finish term 18 yet.

### Fused Numeric Tensor-Sum Contraction

Change: `wxytorxt` plus `ptmztwpo`.

This attempts to contract `(T1 + ... + Tn) * N` without first materializing the
sum, when `N` is a sparse numeric tensor. For larger sums it can return scaled
lazy terms instead of a materialized tensor.

Positive effect:

- This is the right shape for many `Q3`/`δ`/metric/gamma interactions.
- It avoids constructing a huge intermediate tensor before a contraction that
  may shrink it.

Current failure mode:

With sparse shadows alone, and also with a very low lazy-distribution cap, term
18 fails quickly with:

```text
fused tensor sum contract requires at least one shared slot
```

This should be a fast-path miss, not a hard error. The fast path needs to return
`None` when the tensor pair has no shared contraction slot.

### Grouped Sparse Atom Contraction

Change: `ukzyqptw`.

This groups all products contributing to the same output sparse index and then
calls `Atom::add_many` once per output index.

Expected effect:

- Reduces incremental atom growth during sparse-sparse contraction.
- Works best with sparse shadows, because dense shadow tensors hide sparsity.
- Particularly relevant for small metric/`δ`/`Q3` products whose entries collide
  into the same output slot.

Measured separately: no clean standalone timing yet, because the fused numeric
fast path currently fails before a useful isolated sparse-only execution can
complete.

### Lazy Distribution Cap

Change: `ulrtqvwk`, controlled by
`SPENSO_NETWORK_MAX_LAZY_DISTRIBUTED_TERMS`.

Default cap is `96`. With sparse shadows and lazy sums, the profile reaches late
products such as:

```text
product.lazy_tensor_sum_materialize left_terms=2 right_terms=583 distributed_terms=1166 max_distributed_terms=96
```

The default cap prevents unconstrained distribution, but the forced
materialization then leads into slow result-rank pair selection and term 18 does
not finish quickly.

Cap experiments:

| Cap | Result |
| ---: | --- |
| `1` | Fails quickly with no-shared-slot fused contraction error |
| `96` | Runs past `60s`; late `result_rank_slow_pair` calls around `20-28s` |
| `100000` | Runs past `60s`, then fails with `result is an unmaterialized scaled tensor sum` |

Conclusion: the cap is necessary, but the current strategy needs two fixes:

- fused contraction must decline unsupported pairs instead of erroring
- final/result extraction must materialize `TensorTermSum` when it remains at the
  result boundary

### Horner Contract Entries

Change: `topxkqvv`, controlled by
`SPENSO_NETWORK_HORNER_CONTRACT_ENTRIES`.

This applies `collect_horner` to atom entries produced by tensor contractions.
It does not affect parsing or network shape, only scalar payloads after
contraction.

Current status: not benchmarked to completion on term 18 because the execution
path fails or exceeds the timeout before a stable final comparison can be made.
It is still a plausible scalar-payload cleanup pass, but it is downstream of the
current orchestration problems.

### Metric Prepass Diagnostic

Diagnostic knob: `SPENSO_SUM46_PREPASS=metrics`.

Effect on term 18 in `NETWORK_ONLY` mode:

| Mode | Prepass time | Graph nodes | Library leaves | Max product children |
| --- | ---: | ---: | ---: | ---: |
| none | n/a | `2608` | `95` | `12` |
| metrics | `60-70 ms` | `2605` | `92` | `9` |

This removes three top-level metric library leaves and shortens the largest
product, but it barely changes the overall branch. It is useful cleanup, not the
main performance lever.

## Synergies

### Sparse Shadows + One-Hot Selectors

This is a good pairing. Sparse shadows expose that `δ` is one-hot, and one-hot
selector contraction makes those `319` `δ` contractions cheap local selections.
This should stay.

### Sparse Shadows + Fast Sparse Sums

This is also a good pairing. Sparse shadows reduce stored entries from `3088` to
`1678` for term 18 and make sparse sum materialization operate on the actual
support. Without sparse shadows, the sum code still sees dense rank-1 vectors and
does unnecessary work.

### Sparse Shadows + Lazy Tensor Terms + Fused Numeric Contraction

This is the main intended execution strategy. It reaches deeper than the default
dense path, but the current implementation is incomplete:

- unsupported no-shared-slot pairs currently error
- high lazy caps can leave an unmaterialized scaled tensor sum at the result
- default cap eventually materializes a large `2 * 583` lazy product and hits
  expensive result-rank selection

The direction is right, but term 18 still does not execute in reasonable time.

### Metric Prepass + Sparse Shadows

This combines cleanly but is only a small structural cleanup. It reduces library
leaves from `95` to `92` and max product width from `12` to `9`; sparse storage
effects are unchanged.

## Current Conclusion

The improvements do help, but not yet enough for the heavy large-input branch.
The most concrete wins are:

- sparse shadows: `3088 -> 1678` stored param entries for term 18
- metric prepass: small product/library cleanup
- lazy scaled tensor terms: avoids immediate scalar replication and reaches much
  later execution stages

The current blockers are:

- fused numeric tensor-sum contraction treats some fast-path misses as errors
- `TensorTermSum` can survive to the result boundary without materialization
- default lazy distribution still causes late large materializations
- result-rank pair selection remains very expensive once large lazy sums survive
  deep into execution

The next useful fix is to make fused tensor-sum contraction a true optional fast
path: unsupported pairs should return `None`, not `Err`. After that, result
boundary materialization for scaled tensor sums should be fixed so high lazy caps
can be tested honestly.
