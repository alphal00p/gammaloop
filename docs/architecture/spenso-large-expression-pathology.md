# Spenso Large Expression Pathology

This note summarizes the current diagnostics for the large root-level files
`spenso_eval_input_0.txt` and `symbolica_expression.txt`. The focus is the
actual GammaLoop tensor path, not the symbolic parse-only network tests.

## Baseline Shape

For `spenso_eval_input_0.txt`, parsing is not the main bottleneck:

- file read plus tag stripping is well below one second;
- Symbolica parsing is about 3 seconds;
- actual tensor-network parsing is about 3 to 5 seconds, depending on
  diagnostics enabled;
- execution is the part that still fails to complete quickly.

The late expensive tensor sum has this profile in the default path:

```text
leaves=46
logical_entries=188416
atom_entries=188416
zero_entries=78976
nonzero_entries=109440
nonzero_density=0.580842
top_level_sum_entries=0
total_terms=188416
max_terms=1
```

The important point is that the entries are not already top-level Symbolica
sums. They are many simple atom entries spread over a fairly dense tensor
support. The cost is therefore not "simplify a few big scalar sums"; it is
entry-wise materialization over a large support, followed by contractions that
can keep the expression expanded.

## Lazy Tensor Sums

`NetworkLeaf::TensorSum` is available behind:

```bash
SPENSO_NETWORK_LAZY_TENSOR_SUMS=1
```

It stores a tensor-valued sum as a list of tensor leaves and distributes through
one-sided contractions. This proves useful as a diagnostic but is not a default
fix:

- the distilled late-sum MWE remains equivalent;
- one-sided lazy distribution avoids the immediate dense add;
- on `spenso_eval_input_0.txt`, earlier lazy sub-sums accumulate and the run
  reaches a late `103 x 1` distributed product that does not finish quickly.

The conclusion is that lazy sums only help when the later contraction shrinks
each summand enough. They are harmful when they postpone a dense addition into
many expensive products.

## Boundary Factoring

`ParseSettings::factor_add_boundaries` runs Horner collection on large Add
nodes before parsing them as tensor-network sums. It is off by default.

On the distilled MWE it is exactly the right transformation:

```text
parse.add_boundary_factor terms=12 after_terms=1
direct execution: about 28 ms
boundary-factored execution: below 1 ms
```

On `spenso_eval_input_0.txt` the pass does fire:

```text
parse.add_boundary_factor terms=46 before_bytes=39256877 after_terms=4 after_bytes=39254919
```

However, the parsed graph is not smaller in practice:

```text
graph_nodes roughly 23750
graph_edges roughly 38169
tensors 6559
scalars roughly 5764
```

The collected form has fewer top-level terms, but it moves the complexity into
nested products and sums. Boundary factoring is therefore a good local tool and
a good MWE fix, but not sufficient for the full input.

## Staged Index Disconnection

The manual orchestration experiment disconnects a boundary index before
parsing, executes the tensor-valued intermediate, then reconnects with the
corresponding metric and executes the final scalar.

For the late tensor-sum MWE:

```text
direct parse+execute:
  graph_nodes=51
  execution about 28 ms
  result terms=4

disconnected first stage:
  graph_nodes=51
  execution about 1 ms
  tensor result order=2, logical_entries=16

reconnect metric stage:
  execution about 0.15 ms
  expanded scalar difference from direct result is zero
```

This is the clearest evidence so far that orchestration can matter as much as
the local tensor kernels. By delaying one contraction boundary, the executor
contracts the repeated scalar/tensor structure before forming the dense scalar
sum.

For the full `spenso_eval_input_0.txt` case, the same idea needs a reliable way
to choose cut indices from the generated graph or network. Blind text-level
renaming is not enough, because the raw expression contains many occurrences of
the same hedge and edge labels and only some of them are on a useful boundary.

## Graph-Derived Boundary Selection

The current diagnostic path finds staged-disconnection candidates from the
parsed network itself:

1. Walk Product expression nodes.
2. For each immediate Product child, collect its expression subtree.
3. Scan paired slot half-edges.
4. Keep the slot edges whose endpoints land in different Product child
   subtrees.
5. Score the boundary by large-side pressure and small-side cost. The current
   pressure heuristic combines maximum sum fanout, tensor logical entries, and
   scalar byte size.

This selects the same useful MWE boundary without text rewriting:

```text
late_tensor_sum_mwe_graph_selected boundary_candidates total=13
chosen slot=mink4|nu
direct execution: about 28 ms
disconnected stage: about 1 ms
reconnect metric stage: about 0.12 ms
expanded scalar difference from direct result is zero
```

The graph-side transformation has two parts. Splitting the half-edge pair is not
enough by itself, because the tensor store still carries the original slot
labels. The selected cheap-side expression subtree must also have all local
tensor structures reindexed from the original slot to a fresh dummy slot. After
that, execution produces a tensor-valued intermediate, and the final scalar is
obtained by multiplying it by the reconnecting metric
`g(original_slot, fresh_slot)`.

On `spenso_eval_input_0.txt`, the parse-only graph diagnostic now finds a small
number of strong top-level candidates in seconds:

```text
symbolica_parse: about 3.2 s
actual_network_parse: about 2.3 s
boundary_candidates total=15
top candidates: six edges crossing from a huge sum-heavy subtree
  left_nodes=23700
  left_sums=3393
  left_max_sum_children=46
  left_tensor_entries=49852
  left_scalar_bytes=38892668
to single tensor leaves with 16 or 64 tensor entries
```

That is the key distinction from blind renaming. The candidate is a specific
paired slot edge at a specific Product boundary, not a global textual label.
This should let staged execution choose only the boundary that separates the
large repeated subtree from a small external tensor, while leaving other equal
labels inside the large subtree untouched.

## Current Diagnosis

The pathological ingredients are:

1. A large, late tensor-valued sum with dense-ish support.
2. Leaves whose entries are individually simple, so scalar simplification does
   not see an obvious top-level sum to compress.
3. Contraction order that can materialize the sum before the contraction that
   would shrink it.
4. Global Horner collection that reduces some top-level term counts but leaves
   nested tensor-network complexity almost unchanged.

The next structural direction is execution orchestration: identify contraction
boundaries where delaying an index identification keeps an intermediate tensor
small, then reconnect with explicit metrics once the dense repeated structure
has been reduced.
