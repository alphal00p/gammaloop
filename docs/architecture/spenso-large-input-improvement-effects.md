# Spenso Large Input Improvement Effects

This note summarizes the effects of the local improvements stacked on top of
`origin/improve-parsing-and-execution` for `spenso_eval_input_0.txt`. The large
input is still too slow for a broad full-expression matrix, so the measurements
below use term 18 of the 46-way sum, as a representative heavy branch.

Unless stated otherwise, commands used a five minute cap:

```bash
timeout 300s env SPENSO_SUM46_START=18 SPENSO_SUM46_LIMIT=1 \
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

## Five Minute Timing Matrix

The previous measurements were inconclusive because several configurations were
only watched for about one minute. With a `300s` timeout and the fast-path bugs
fixed, term 18 now gives a useful matrix:

| Mode | Parse | Execute | Total | Result |
| --- | ---: | ---: | ---: | --- |
| dense shadows, no lazy sums | `360.8 ms` | `>300s` | timeout | did not finish |
| sparse shadows | `361.7 ms` | `99.405s` | `113.159s` | scalar, `9.25 GB` atom payload |
| sparse shadows + lazy sums, cap `1` | `365.2 ms` | `108.211s` | `119.691s` | scalar, `9.46 GB` atom payload |
| sparse shadows + lazy sums, cap `96` | `363.8 ms` | `101.374s` | `112.171s` | scalar, `8.99 GB` atom payload |
| sparse shadows + lazy sums, cap `100000` | `363.2 ms` | `57.842s` | `60.156s` | scalar from rank-0 tensor, `2966` terms, `252 MB` atom payload |
| sparse shadows + metric prepass | `363.5 ms` | `116.344s` | `128.772s` | scalar, `10.00 GB` atom payload |
| sparse shadows + Horner contract entries | `363.6 ms` | `151.539s` | `160.749s` | scalar, `9.25 GB` atom payload |

Two bugs had to be fixed to get this matrix:

- unsupported fused tensor-sum/numeric contractions with no shared slot now
  decline the fast path instead of erroring
- result extraction now materializes `TensorSum`, `TensorTerm`, and
  `TensorTermSum` leaves when lazy terms survive to the result boundary

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

Separate execution effect: this is the main practical win so far. The dense
shadow path did not finish within `300s`; sparse shadows complete term 18 in
`99.405s` execution time. The stored-entry reduction also makes the one-hot and
grouped sparse contraction paths visible to the executor.

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

Measured separately: parse/network-only shape is unchanged and there is no clean
standalone runtime knob for only this improvement. In the current stack, sparse
shadows plus the always-on sparse sum materialization complete term 18 in
`99.405s`; with high-cap lazy distribution the same branch completes in
`57.842s`, plus `1.510s` to force the final atom.

### Lazy Scaled Tensor Terms

Change: `wxytorxt`.

This adds `TensorTerm` and `TensorTermSum`, representing `scalar * tensor`
without multiplying the scalar into every tensor entry immediately.

Expected effect:

- Keeps common scalar factors outside tensor entries.
- Prevents scalar payload replication across many sparse tensor entries.
- Enables lazy tensor sums to carry scaled terms forward.

Observed effect with sparse shadows and lazy sums enabled: default cap `96`
completes, but is roughly neutral on execution time (`101.374s` versus
`99.405s`). It does reduce final scalar payload size from about `9.25 GB` to
about `8.99 GB`, so it is helping representation size without yet improving this
branch's runtime.

### Fused Numeric Tensor-Sum Contraction

Change: `wxytorxt` plus `ptmztwpo`.

This attempts to contract `(T1 + ... + Tn) * N` without first materializing the
sum, when `N` is a sparse numeric tensor. For larger sums it can return scaled
lazy terms instead of a materialized tensor.

Positive effect:

- This is the right shape for many `Q3`/`δ`/metric/gamma interactions.
- It avoids constructing a huge intermediate tensor before a contraction that
  may shrink it.

Fixed failure mode: this used to error on tensor pairs with no shared slot. That
case is now treated as a fast-path miss, so sparse-only and low-cap lazy
execution both complete. The low cap `1` run is slower than sparse-only
(`108.211s`), which suggests forced early materialization is not useful for term
18.

### Grouped Sparse Atom Contraction

Change: `ukzyqptw`.

This groups all products contributing to the same output sparse index and then
calls `Atom::add_many` once per output index.

Expected effect:

- Reduces incremental atom growth during sparse-sparse contraction.
- Works best with sparse shadows, because dense shadow tensors hide sparsity.
- Particularly relevant for small metric/`δ`/`Q3` products whose entries collide
  into the same output slot.

Measured separately: no clean standalone timing is available from a runtime knob,
because this path is active once the code is present. In the current stack, the
closest combined measurement is sparse shadows without lazy sums: `99.405s`
execution time.

### Lazy Distribution Cap

Change: `ulrtqvwk`, controlled by
`SPENSO_NETWORK_MAX_LAZY_DISTRIBUTED_TERMS`.

Default cap is `96`. With sparse shadows and lazy sums, the profile reaches late
products such as:

```text
product.lazy_tensor_sum_materialize left_terms=2 right_terms=583 distributed_terms=1166 max_distributed_terms=96
```

The default cap prevents unconstrained distribution, but the forced
materialization still leads into slow scalar growth. It now completes, but is not
faster than sparse-only.

Cap experiments:

| Cap | Result |
| ---: | --- |
| `1` | completes in `108.211s`; scalar result payload about `9.46 GB` |
| `96` | completes in `101.374s`; scalar result payload about `8.99 GB` |
| `100000` | completes in `57.842s`; final rank-0 tensor extraction to an atom takes `1.510s`, producing `2966` terms and about `252 MB` |

Conclusion: the high cap is the best timing here, even after forcing a final
atom. Avoiding forced intermediate materialization matters much more than the
`1.510s` final rank-0 tensor-to-atom extraction.

### Horner Contract Entries

Earlier experiment: `topxkqvv`, originally controlled by the now-removed
`SPENSO_NETWORK_HORNER_CONTRACT_ENTRIES` environment variable. The current
version selects the same idea explicitly through the contraction strategy type,
for example `MinResultRank<HornerAtomComponents<4096>>`.

This applies `collect_horner` to atom entries produced by tensor contractions.
It does not affect parsing or network shape, only scalar payloads after
contraction.

Measured with sparse shadows enabled, this makes term 18 slower:

```text
execute=151.539s total=160.749s
```

It leaves final payload size roughly unchanged for this branch. This is not a
low-hanging win in its current placement.

### Metric Prepass Diagnostic

Diagnostic knob: `SPENSO_SUM46_PREPASS=metrics`.

Effect on term 18 in `NETWORK_ONLY` mode:

| Mode | Prepass time | Graph nodes | Library leaves | Max product children |
| --- | ---: | ---: | ---: | ---: |
| none | n/a | `2608` | `95` | `12` |
| metrics | `60-70 ms` | `2605` | `92` | `9` |

Full execution with sparse shadows enabled:

```text
prepass=60.859ms execute=116.344s total=128.772s
```

This removes three top-level metric library leaves and shortens the largest
product, but it makes term 18 slower. It is useful cleanup, not the main
performance lever.

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

This is the main intended execution strategy. With the fast-path and result
boundary bugs fixed, it completes. The cap matters:

- cap `1`: `108.211s`
- cap `96`: `101.374s`
- cap `100000`: `57.842s` execution plus `1.510s` final atom extraction

The high-cap run is the first result here that is plausibly close to reasonable
for one heavy branch. It also shows that premature materialization is currently
more harmful than allowing a large lazy sum to survive until the result boundary,
where it can be collapsed to an atom cheaply relative to the contraction.

### Metric Prepass + Sparse Shadows

This combines cleanly but is only a small structural cleanup. It reduces library
leaves from `95` to `92` and max product width from `12` to `9`; sparse storage
effects are unchanged. Full execution is slower for term 18 (`116.344s` versus
`99.405s` sparse-only).

## Current Conclusion

The improvements now make term 18 executable under a five minute cap. The most
concrete wins are:

- sparse shadows: `3088 -> 1678` stored param entries for term 18
- sparse shadows alone: dense timeout `>300s` becomes `99.405s`
- high-cap lazy distribution: `99.405s` becomes `57.842s`, with another
  `1.510s` to force the final atom

The non-wins are also clear:

- default lazy cap `96` is neutral/slightly slower than sparse-only
- cap `1` is slower, so forced early materialization is bad here
- metric prepass and Horner contract entries are both slower on term 18

The next useful direction is not more scalar cleanup at the end. The strongest
signal is that orchestration should preserve laziness through contraction
boundaries for longer, and only materialize when the contraction has genuinely
reduced the tensor/sum shape.
