# Network-Based Symbolic Simplification Status

## Scope

This note summarizes the current state of the `idenso` network-based
Schoonschip-style symbolic simplification path, with emphasis on vertex
algebra benchmarks and the knobs available for further investigation.

The relevant implementation is mainly in:

- `crates/idenso/src/schoonschip/api.rs`
- `crates/idenso/src/schoonschip/contraction.rs`
- `crates/idenso/src/schoonschip/settings.rs`
- `crates/idenso/src/schoonschip/utils.rs`
- `crates/idenso/benches/vertex_algebra_once.rs`

## Current Model

The network simplifier parses a Symbolica expression into a tensor network of
`SymbolicTensor`s, executes contractions in a configurable order, and uses
Schoonschip-style pattern simplification at contraction boundaries.

The important distinction is between:

- **bare symbolic substitution**: pattern matching over already-expanded
  expressions, used as the reference behavior in the vertex benchmarks.
- **network-informed contraction**: uses network slots to identify dummy pairs
  and attempts direct replacement of the matching index/vector/metric before
  falling back to explicit distribution.

The network path should be faster for shapes such as:

```text
p(mu) * sum(... mu ...)
g(mu, nu) * sum(... mu ...)
```

because the network already knows the contracted slot and can rewrite the
target side directly. These are not optional shortcuts; they are the mechanism
that lets the network path avoid expanding one-sided sums.

## Current Direct Sum Contraction Path

For sum-by-sum contractions, the current path is:

1. Only enter the expanded-sum logic when both contracted operands are sums.
2. Pick the smaller expanded side term-by-term.
3. For each term, parse tensor factors and match them against contracted slot
   pairs.
4. Build direct replacements:
   - rank-1 tensor factor: replace target slot by the stripped vector
   - metric with one matched slot: replace target slot by the other metric slot
   - metric with both matched slots: identify the two target slots
5. Apply replacements transitively to the target expression.
6. Run compact cleanup with `schoonschip()`.
7. If removed dummy slots remain, try bounded local residual-boundary cleanup.
8. If residual slots still remain, fall back to
   `distribute_smallest_expanded_sum_side(...).schoonschip_with_net::<false, false, ...>(...)`.

The fallback is the expensive path. It explicitly distributes the smaller sum
side, then runs network simplification again on the expanded boundary.

## Current Findings

The direct replacement path works for genuinely local boundaries. A small MWE
where compact cleanup leaves `mu2` but full expansion removes it is now handled
by bounded local target-boundary expansion.

The 7-vertex full-algebra case is still not solved by this local recursion. The
trace reaches a residual boundary roughly of this shape:

```text
target_terms=36
target_bytes=5812309
residual slot=mink(4, mu4)
```

The local recursion is capped and immediately skips this target:

```text
local_boundary skip_target_size target_bytes=5812309 max_bytes=500000
```

Earlier uncapped experiments showed that the relevant local residual product
can expand to 1296 terms and still fail on at least one residual-bearing term.
So the current conclusion is:

- local residual-boundary recursion is useful for small cases;
- it is not the structural fix for the 7-vertex slowdown;
- the remaining slowdown is the large fallback after compact direct cleanup
  leaves a residual dummy.

## Known Reproduction Commands

Run the one-shot benchmark harness:

```bash
cargo bench --package idenso --bench vertex_algebra_once --profile dev-optim -- network_full_algebra_7 smallest_degree
```

Available full-algebra filters:

```text
network_full_algebra_5
network_full_algebra_6
network_full_algebra_7
network_full_algebra_8
```

Reference bare substitutions:

```text
bare_vertex_substitution_5
bare_vertex_substitution_6
bare_vertex_substitution_8
```

The harness asserts that no internal network indices remain unless `no_assert`
is passed as the third argument:

```bash
cargo bench --package idenso --bench vertex_algebra_once --profile dev-optim -- network_full_algebra_7 smallest_degree no_assert
```

## Public Settings Knobs

`SchoonschipSettings` exposes these main controls:

- `depth_limit`: parse depth for symbolic network parsing.
- `mode`:
  - `single_pass`
  - recursive depth-first
  - recursive breadth-first
- `parse_inner_products`: internal control; recursion disables this for scalar
  sub-simplification.
- `expand_contracted_sums`: forces contracted sums through the expanded
  contraction path.
- `simplify_chain_like_functions`: enables chain-like metric product
  simplification.
- `contraction_order`: selects the network contraction ordering strategy.

Convenience constructors:

```text
SchoonschipSettings::single_pass(...)
SchoonschipSettings::depth_first(...)
SchoonschipSettings::breadth_first(...)
SchoonschipSettings::partial()
SchoonschipSettings::full()
SchoonschipSettings::with_expanded_contracted_sums()
SchoonschipSettings::with_chain_like_functions()
SchoonschipSettings::with_contraction_order(...)
```

## Contraction Ordering Knobs

The current ordering enum is:

```text
SmallestDegree
LargestDegree
MinLargestOperandBytes
MinProductTerms
MinProductBytes
SmallestDegreeMinLargestOperandBytes
SmallestDegreeMinProductTerms
SmallestDegreeMinProductBytes
```

The one-shot benchmark accepts the same names:

```text
smallest_degree
largest_degree
min_largest_operand_bytes
min_product_terms
min_product_bytes
smallest_degree_min_largest_operand_bytes
smallest_degree_min_product_terms
smallest_degree_min_product_bytes
```

The expression-order metrics score candidate contractions using combinations of:

- whether the edge is internal or degree-zero;
- contraction degree;
- largest operand byte size;
- sum of operand byte sizes;
- product of term counts;
- product of byte sizes.

Current investigation status:

- `smallest_degree` is the baseline, matching normal network execution logic.
- pure product metrics can pick contraction orders that look small locally but
  may fail to simplify fully.
- degree-first plus product/byte tie-breakers are safer candidates, but still
  need residual-index checks on 6, 7, and 8 vertices.

## Direct Sum and Boundary Knobs

Direct sum contraction can be disabled at runtime:

```bash
IDENSO_DISABLE_DIRECT_SUM_CONTRACTIONS=1 ...
```

This is useful for comparing direct network-informed replacement against the
fallback explicit-distribution path.

Current hard-coded local residual-boundary caps:

```text
LOCAL_BOUNDARY_MAX_EXPANDED_TERMS = 512
LOCAL_BOUNDARY_MAX_PRODUCT_BYTES = 100000
LOCAL_BOUNDARY_MAX_TARGET_BYTES = 500000
```

These caps deliberately keep the local recursive cleanup from becoming another
large expansion path. Increasing them may make small repros pass, but the
uncapped 7-vertex experiment showed that this does not address the main
slowdown.

## Trace and Profiling Knobs

Useful environment variables:

```text
IDENSO_TRACE_SUM_CONTRACTIONS
IDENSO_TRACE_CONTRACTION_ORDERING
IDENSO_TRACE_DIRECT_SUM_TERMS
IDENSO_TRACE_DIRECT_SUM_TERM_EXPRESSIONS
IDENSO_TRACE_FINISH_CONTRACTS
IDENSO_TRACE_SCHOONSCHIP_PATTERNS
IDENSO_TRACE_SCHOONSCHIP_PATTERN_MISSES
```

Common trace command:

```bash
IDENSO_TRACE_CONTRACTION_ORDERING=1 \
IDENSO_TRACE_SUM_CONTRACTIONS=1 \
IDENSO_TRACE_DIRECT_SUM_TERMS=1 \
cargo bench --package idenso --bench vertex_algebra_once --profile dev-optim -- network_full_algebra_7 smallest_degree
```

Instruction-count profiling can be done with the Linux-only Valgrind/Callgrind
flake additions from the current investigation branch. Keep those dependencies
Linux-gated because macOS devices using the flake do not have Valgrind.

## Current Bottleneck Shape

The key bad shape is not plain dot normalization. The traces point to compact
direct cleanup producing a large target expression with a residual dummy slot.
Once that happens, the algorithm must use the fallback distribution path to
recover full simplification.

The difficult case is a contraction that requires looking inside a product of
sums after direct slot replacement. Compact cleanup can preserve the product of
sums as one term, while full simplification needs enough boundary expansion to
pair all residual occurrences correctly.

For small products this is now handled locally. For the 7-vertex case the
boundary is already too large, so "recursively analyze then locally expand" is
not sufficient.

## Open Directions

The most promising remaining directions are:

1. **Change contraction order to avoid the bad boundary.**
   The 7-vertex slowdown may be avoidable if the network never creates the
   36-term, multi-megabyte residual target before the dummy is eliminated.

2. **Represent residual sum contraction without expanding pairings blindly.**
   Direct slot replacement works because the network knows the dummy. The next
   step would need to preserve that information across a product of sums and
   contract the residual slot structurally rather than by fallback expansion.

3. **Improve failure diagnostics for ordering methods.**
   Every ordering should report whether it fully removed internal indices, not
   only elapsed time and output term count.

4. **Keep local boundary cleanup bounded.**
   It fixes the small MWE and should remain a guarded cleanup pass, not a broad
   expansion mechanism.

5. **Compare against bare symbolic output by term count and residual indices.**
   Timing alone is not enough: some network results are faster because they are
   under-simplified.

## Practical Checklist for Next Experiments

For each candidate ordering or boundary rewrite:

1. Run 5 and 6 vertices first and compare against bare symbolic term counts.
2. Run 7 vertices in release/dev-optim with residual-index assertions enabled.
3. If it fails, rerun with `no_assert` and trace residual slots.
4. If it succeeds, compare output terms and bytes against the bare symbolic
   result.
5. Only then run 8 vertices.
