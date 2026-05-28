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
  sparse-aware path in about 4 seconds.
- `spenso_eval_input_0.txt` can be parsed, aliased, and analyzed for graph
  boundary candidates in a few seconds without exercising the old runaway
  execution path.
- The BNL four-term MWE completes in milliseconds when scalar aliases are kept
  unresolved, and cleanly reproduces the benefit of not carrying large parsed
  scalar payloads through every tensor entry.
- Boundary factoring was rerun on the large inputs with scalar aliasing. It
  remains useful on the distilled MWE, but it does not currently solve the
  large inputs: `symbolica_expression.txt` fails during boundary-factored
  execution, and `spenso_eval_input_0.txt` remains only a bounded
  parse/candidate diagnostic.
- The late tensor-sum MWE is no longer a 10s-of-ms direct execution problem
  under the current order; direct execution is around 0.5 ms, and boundary
  factoring/Hornering still improves it further.
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
sparse-aware pair score. Full execution of `spenso_eval_input_0.txt` was not
rerun in this pass; only parse, aliasing, and boundary-candidate diagnostics
were rerun for that larger input.

## `spenso_eval_input_0.txt`

Parsing is not the current bottleneck:

```text
input bytes=53331469
symbolica_parse: 1.348 s, terms=1, bytes=39268617
actual_network_parse: 1.061 s
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
symbolica_parse: 100 ms, terms=1, bytes=3016471
actual_network_parse: 623 ms
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
sparse-aware execute: 4.110 s
result construction: 3.365 ms
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

## Boundary Factoring

`ParseSettings::factor_add_boundaries` runs Horner collection on large Add
nodes before parsing them as tensor-network sums. It is off by default.

On the distilled late tensor-sum MWE, boundary factoring and Hornering remain
effective. This was rerun with sparse-aware order and the scalar alias pass at
the 4096-byte threshold; aliasing creates no aliases for this MWE, so these
numbers are effectively the current sparse-aware execution costs:

```text
parse.add_boundary_factor terms=12 after_terms=1
expanded:
  aliases=0
  execute=0.983 ms
  root_terms=4
  root_bytes=10374
boundary-factored:
  aliases=0
  execute=0.054 ms
  root_terms=1
  root_bytes=2389
factored:
  aliases=0
  execute=0.041 ms
  root_terms=1
  root_bytes=2389
Hornered:
  aliases=0
  execute=0.036 ms
  root_terms=1
  root_bytes=2389
```

On `spenso_eval_input_0.txt`, the pass still fires:

```text
parse.add_boundary_factor terms=46 before_bytes=39256877 after_terms=4 after_bytes=39254919
```

Historically, this did not make the full parsed graph much smaller:

```text
graph_nodes roughly 23750
graph_edges roughly 38169
tensors 6559
scalars roughly 5764
```

So boundary factoring is a good local improvement and a useful MWE fix, but it
does not by itself solve the large root-level inputs.

For `symbolica_expression.txt`, boundary factoring with scalar aliasing parses
to a smaller graph than the unfactored form:

```text
actual_network_parse: 3.208 s
graph_nodes=61757
graph_edges=111171
tensors=10980
scalars=11664
logical_entries=43920
scalar_aliases threshold_bytes=4096
  created=22
  terms=22
  bytes=100045
  max_bytes=5717
```

but execution currently fails:

```text
error: product contraction did not collapse to one leaf; 9 operands remain
```

So this is not an improvement over the unfactored aliased sparse-aware run,
which executes in about 4 seconds.

For `spenso_eval_input_0.txt`, boundary factoring with scalar aliasing was
rerun only through parse and boundary-candidate diagnostics:

```text
symbolica_parse: 1.416 s
actual_network_parse: 1.293 s
graph_nodes=23750
graph_edges=38169
tensors=6559
scalars=5764
logical_entries=26236
scalar_aliases threshold_bytes=4096
  created=369
  terms=3138
  bytes=37513085
  max_bytes=1856190
post-alias leaf scalar_bytes=1392651
boundary_candidates total=15
top candidates: six edges crossing from a huge sum-heavy subtree
  left_nodes=23739
  left_sums=3406
  left_max_sum_children=23
  left_tensor_entries=49852
  left_scalar_bytes=1392641
to single tensor leaves with 16 or 64 tensor entries
```

The candidate shape remains essentially the same as the unfactored aliased
diagnostic. The top-level maximum sum arity drops from 46 to 23, but the graph
does not become smaller and the strong boundary candidates remain the same kind
of large-subtree-to-small-leaf edges.

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
finishes, and `spenso_eval_input_0.txt` can be parsed and analyzed with aliases
without the diagnostic itself being dominated by parsed scalar bytes.

The remaining work is to control scalar materialization during or immediately
after execution. Good next directions are:

- re-alias generated scalar entries during execution or before result
  extraction;
- make `AliasedAtom` construction avoid inlining newly generated scalar payloads;
- add profiling for where the largest post-execution scalar entries are
  produced;
- revisit staged execution only after products that do not collapse to one leaf
  can be represented or resumed safely.
