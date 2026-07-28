# Spenso Large Expression Pathology

This note tracks the current status of the large root-level Spenso diagnostic
inputs:

- `spenso_eval_input_0.txt`
- `symbolica_expression.txt`
- `examples/cli/BNL/profiling/bnl_integrated_evaluator_atom_unfiltered_pre_network.sym`

The focus is the actual GammaLoop tensor-network path. The current reruns use
the sparse-aware contraction order through `MinResultRank` and scalar aliasing
before execution where execution is attempted.

## Current Status

The situation has improved compared to the original investigation.

- `symbolica_expression.txt` now completes end-to-end in the aliased
  sparse-aware path in about 4.6 seconds.
- `spenso_eval_input_0.txt` can be parsed, aliased, and analyzed for graph
  boundary candidates in a few seconds. Full execution currently fails
  structurally with a disconnected product, both with and without
  component-level Hornering.
- The BNL four-term MWE completes in milliseconds when scalar aliases are kept
  unresolved, and cleanly reproduces the benefit of not carrying large parsed
  scalar payloads through every tensor entry.
- Component-level Hornering is the Spenso-owned optimization: generated tensor
  components can be Horner-collected after concrete/parametric component
  contraction, using the contraction strategy marker
  `HornerAtomComponents<N>`.
- The late tensor-sum MWE is no longer a 10s-of-ms direct execution problem
  under the current order; direct execution is sub-millisecond to
  low-millisecond in these diagnostics. Component-level Hornering remains a
  useful comparison.
- A naive `result-rank-only` order is still bad enough to get killed on
  `symbolica_expression.txt`, so contraction order does matter.

What remains unresolved is scalar materialization during execution. Pre-execution
aliasing only aliases scalar refs that already exist after parsing. Large scalar
atoms can still be generated later by contractions, so an aliased result can
still have a very large root.

## Rerun Context

The current diagnostics were run from the `large_spenso_actual` test target with
bounded settings:

```bash
RAYON_NUM_THREADS=1
CARGO_BUILD_JOBS=2
SPENSO_ALIAS_SCALAR_THRESHOLD=4096
```

The contraction order is the current default `MinResultRank`, which is the
sparse-aware pair score. `spenso_eval_input_0.txt` was attempted both with and
without component-level Hornering; both runs failed before producing a scalar
result because the product executor reached disconnected tensor operands.

## `spenso_eval_input_0.txt`

Parsing is not the current bottleneck:

```text
input bytes=53331469
symbolica_parse: 1.441 s, terms=1, bytes=39268617
actual_network_parse: 1.129 s
graph_nodes=23711
graph_edges=38052
tensors=6559
scalars=5751
logical_entries=26236
param_entries=26236
```

Applying pre-execution scalar aliasing changes the scalar-pressure accounting
substantially:

```text
scalar_aliases threshold_bytes=4096
  created=369
  terms=3138
  bytes=37515435
  max_bytes=1856190
post-alias leaf scalar_bytes=1392399
```

Plain and component-Horner execution both fail after parsing and aliasing with
the same residual product shape:

```text
error: product contraction did not collapse to one leaf; 3 operands remain
operands=[
  #0:Scalar:no_tensor_structure,
  #1:LocalTensor:order=1,dense_size=4,
  #2:LocalTensor:order=2,dense_size=16
]
matching_pairs=[]
max_rss=646053888 bytes
```

So for this input, the failure should not be attributed to component-level
Hornering. The product has already reached a disconnected tensor-valued state
with no matching index pairs.

The graph-derived product-boundary diagnostic still finds the same small set of
strong top-level candidates:

```text
boundary_candidates total=15
top candidates: six edges crossing from a huge sum-heavy subtree
  left_nodes=23700
  left_sums=3393
  left_max_sum_children=46
  left_tensor_entries=49852
  left_scalar_bytes=1392389
to single tensor leaves with 16 or 64 tensor entries
```

This is useful, but it should not be read as proof that staged execution is
currently solved. It only says the graph still exposes clear boundary
candidates once large parsed scalars are aliased.

## `symbolica_expression.txt`

This input now executes end-to-end with pre-execution aliasing and the current
sparse-aware order:

```text
input bytes=6529054
stripped bytes=4720678
symbolica_parse: 109 ms, terms=1, bytes=3016471
actual_network_parse: 710 ms
graph_nodes=69975
graph_edges=123483
tensors=13720
scalars=11663
logical_entries=54880
```

The execution finishes:

```text
scalar_aliases threshold_bytes=4096
  created=22
  terms=22
  bytes=100397
  max_bytes=5733
sparse-aware execute: 4.567 s
result construction: 3.345 ms
result aliases=22
root_terms=1
root_bytes=204416860
```

The important part is the last line. Execution finishes quickly enough, but the
unresolved aliased root is still about 204 MB. That large root is not explained
by scalar refs present immediately after parsing: the alias pass only found 22
scalars above the 4096-byte threshold. The large scalar payload is generated
during execution.

Lowering the alias threshold to 64 does not materially change the result:

```text
threshold_bytes=64
aliases=686
aliased bytes=1136803
execute=4.113 s
root_bytes=203392391
```

That confirms the issue is not just that the threshold was too high. More
pre-execution aliases help only marginally because most of the final root is
created after the alias pass.

## BNL Four-Term MWE

The BNL standalone evaluator input used while debugging aliasing is:

```text
examples/cli/BNL/profiling/bnl_integrated_evaluator_atom_unfiltered_pre_network.sym
```

With sparse-aware order and aliases at threshold 4096, it runs quickly and keeps
the unresolved aliased root small:

```text
input bytes=2038141
symbolica_parse: about 45 ms, terms=1, bytes=1410951
network_parse: about 150 to 310 ms
graph_nodes=98
graph_edges=179
tensors=24
scalar_aliases threshold_bytes=4096
  created=4
  terms=4
  bytes=1405119
  max_bytes=953984
network_execute: about 4 to 7 ms
tensor_entry_scalars after_execute:
  tensors=103
  entries=1669
  total_terms=2460
  total_bytes=1493281
  max_bytes=545943
result_aliased_root:
  aliases=4
  terms=1
  bytes=545960
```

Without scalar aliases, the tensor topology is essentially the same, but the
scalar payload carried in tensor entries is much larger:

```text
network_execute: about 77 ms
tensor_entry_scalars after_execute:
  tensors=103
  entries=1669
  total_terms=2460
  total_bytes=128217701
  max_bytes=29162247
result_atom:
  terms=1
  bytes=29167519
```

This MWE is therefore a good alias-pressure reproducer. It is not a good
contraction-order reproducer: with aliases enabled, `atom-aware` and
`sparse-aware` produce the same post-execution tensor-entry profile and the same
unresolved root size.

## Contraction Order

Ordering matters, but the current sparse-aware order does not improve this
specific `symbolica_expression.txt` case relative to `atom-aware`:

```text
result-rank-only:
  killed by SIGKILL before producing an execution summary

atom-aware:
  execute=3.926 s
  root_bytes=204416860

sparse-aware:
  execute=4.072 s
  root_bytes=204416860
```

`atom-aware` and `sparse-aware` also finish with the same post-execution tensor
profile:

```text
graph_nodes=1
tensors=70658
scalars=17839
logical_entries=1416590
param_entries=823886
orders={0: 686, 1: 30184, 2: 26068, 3: 13720}
```

So the improvement is not that sparse-aware finds a smaller final scalar for
this input. The improvement is that the current execution path finishes for the
safe orders, while naive result-rank-only ordering remains pathological. The
remaining large-root problem needs a materialization/aliasing fix inside or
after execution, not only a lower pre-execution alias threshold.

## Component-Level Hornering

Spenso's Hornering hook is now component-level. A contraction strategy such as
`MinResultRank<HornerAtomComponents<4096>>` keeps the same parsed graph and
only Horner-collects generated `Atom` entries inside parametric tensor
components after contraction. Lowering the threshold, for example
`HornerAtomComponents<1>`, forces the behavior on small MWEs.

On the distilled late tensor-sum MWE, component-level Hornering is compared
against the expanded baseline:

```text
expanded:
  aliases=0
  execute=2.000 ms
  root_terms=4
  root_bytes=10374
component-Hornered:
  aliases=0
  execute=0.050 ms
  root_terms=1
  root_bytes=2389
```

For tensor-valued products like

```text
(x1*Q1(mu) + x2*Q3(mu)) * (Q2(rho) + Q4(rho))
```

the graph is valid: it is an outer product with free `mu` and `rho`
components. Component Hornering should act on the generated entries, for
example

```text
(x1*Q1_i + x2*Q3_i) * (Q2_j + Q4_j)
```

not on the full raw expression before the graph is built.

For `symbolica_expression.txt`, the component-Horner diagnostic uses the same
aliased sparse-aware network:

```text
actual_network_parse: 668 ms
graph_nodes=69975
graph_edges=123483
tensors=13720
scalars=11663
scalar_aliases threshold_bytes=4096
  created=22
  terms=22
  bytes=100397
  max_bytes=5733
execute=5.387 s
result=3.632 ms
root_terms=1
root_bytes=200970362
max_rss=3338059776 bytes
```

Compared to the plain aliased sparse-aware run, component Hornering is slower
here (`5.387 s` vs `4.567 s`) and only slightly reduces the unresolved root
payload (`200.97 MB` vs `204.42 MB`). The large generated root remains the same
issue.

For `spenso_eval_input_0.txt`, the normal parse with scalar aliasing was rerun
through parse, aliasing, plain execution attempt, component-Horner execution
attempt, and boundary-candidate diagnostics:

```text
symbolica_parse: 1.441 s
actual_network_parse: 1.129 s
graph_nodes=23711
graph_edges=38052
tensors=6559
scalars=5751
logical_entries=26236
scalar_aliases threshold_bytes=4096
  created=369
  terms=3138
  bytes=37515435
  max_bytes=1856190
post-alias leaf scalar_bytes=1392399
boundary_candidates total=15
top candidates: six edges crossing from a huge sum-heavy subtree
  left_nodes=23700
  left_sums=3393
  left_max_sum_children=46
  left_tensor_entries=49852
  left_scalar_bytes=1392389
to single tensor leaves with 16 or 64 tensor entries
```

The strong boundary candidates remain large-subtree-to-small-leaf edges.
Component-level Hornering does not change this graph shape. The plain and
component-Horner execution attempts both fail before reaching a single scalar
result with the same disconnected residual product.

## Staged Disconnection

The old staged-disconnection experiment is stale under the current order.

The direct late tensor-sum MWE executes:

```text
direct parse+execute:
  graph_nodes=51
  execution about 0.5 ms
  result terms=4, bytes=10374
```

But the disconnected tensor stage now fails:

```text
disconnected first stage:
  graph_nodes=51
  error: product contraction did not collapse to one leaf; 2 operands remain
```

The graph-selected version still chooses the expected boundary:

```text
late_tensor_sum_mwe_graph_selected boundary_candidates total=13
chosen slot=mink4|nu
```

but it fails at the same disconnected execution stage. This means boundary
selection is still useful diagnostic information, but staged execution is not
currently a working mitigation unless the executor can represent partially
contracted tensor-valued products.

## Lazy Tensor Sums

`NetworkLeaf::TensorSum` is still available behind:

```bash
SPENSO_NETWORK_LAZY_TENSOR_SUMS=1
```

This remains diagnostic-only. The earlier conclusion still applies: lazy tensor
sums help when later contractions shrink each summand enough, but can be harmful
when they postpone dense addition into many expensive products.

## Current Diagnosis

The current pathological ingredients are:

1. Dense-ish tensor support with many simple entries, so scalar simplification
   does not see a small number of large top-level sums to compress.
2. Contraction sequences that generate large scalar atoms during execution.
3. Pre-execution scalar aliasing that successfully removes parsed scalar
   pressure, but cannot alias scalar atoms that do not exist yet.
4. Some contraction orders that are still catastrophically bad
   (`result-rank-only`), even though the current safe orders finish.
5. Staged-disconnection ideas that need a richer executor representation before
   they can be used again.

The current improvement is real: the `symbolica_expression.txt` diagnostic now
finishes, and `spenso_eval_input_0.txt` can be parsed and aliased without the
diagnostic itself being dominated by parsed scalar bytes.

The remaining work is to control scalar materialization during or immediately
after execution. Good next directions are:

- re-alias generated scalar entries during execution or before result
  extraction;
- make `AliasedAtom` construction avoid inlining newly generated scalar payloads;
- add profiling for where the largest post-execution scalar entries are
  produced;
- revisit staged execution only after products that do not collapse to one leaf
  can be represented or resumed safely.
