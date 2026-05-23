# Target UV Structure

This note records the desired structure for the local expanded-4D UV path.

## Core principle

The UV forest recursion should remain in pure 4D atom space.

The CFF/LTD/3Drep projection must happen only after the relevant local 4D forest
term has been fully assembled. In particular, `map_to_3drep` must not be called
inside every recursive `ApproximationKernel::kernel` invocation.

## Desired kernel shape

Introduce a first-class local 4D approximation kernel:

```rust
pub struct Local4DApproximation;

impl ApproximationKernel<UVCtx<'_>> for Local4DApproximation {
    fn kernel<S: ForestNodeLike>(
        &self,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        integrand: &Atom,
    ) -> Result<Atom> {
        match current.renormalization_scheme() {
            ApproximationType::MUV => self.t(
                ctx,
                current,
                given,
                &self.start(ctx, current, given, integrand)?,
            ),
            ApproximationType::IR => Err(eyre!("Not yet implemented IR")),
            ApproximationType::VaccuumLimit => Err(eyre!("Not yet implemented VaccuumLimit")),
            ApproximationType::OS => Err(eyre!("Not yet implemented OS")),
            ApproximationType::Unsubtracted => {
                panic!("should have been kept out of the wood");
            }
        }
    }
}
```

This should mirror `Integrated<'_>` as closely as possible, except that it stops
before `integrate_and_truncate(...)`.

Schematically:

```text
Integrated:
    start -> t -> integrate_and_truncate via Vakint

Local4D:
    start -> t
```

where:

```text
start = integrated_uv_start(...)
t     = integrated_uv_rescale_series(...) followed by projection to strict 4D
```

The strict-4D projection should set `dim_epsilon -> 0`, `dim -> 4`, and map
Minkowski dimensions to 4.

## Iterated T operators

For nested UV spinneys, for example:

```text
gamma_1 subset gamma_2 subset G
```

the local 4D forest term should be assembled recursively in 4D atom space:

```text
T_gamma2^4D( T_gamma1^4D(I_gamma1) * I_{gamma2/gamma1} ) * I_{G/gamma2}
```

Therefore the recursive approximation kernel must return an `Atom`, not a
projected CFF/LTD residue map.

Calling `map_to_3drep` inside every kernel invocation would project too early and
would break this recursive structure.

## Projection boundary

The 3Drep projection happens after the local 4D forest atom has been assembled
and after the remaining factorized graph contribution has been inserted
multiplicatively:

```text
map_to_3drep( local_4d_forest_atom * factorized_cograph_atom )
```

This projection expands the product into CFF/LTD orientation and residue
components:

```text
map_to_3drep(A / product_e D_e)
  =
  sum_orientations sum_residue_components
      weight(orientation, residue) * A_on_that_component
```

The output should be keyed in the top-level residue selector convention:

```rust
BTreeMap<CutCFFIndex, Atom>
```

`CutCFFIndex` must refer only to top-level selected residue axes:

```text
left_threshold_order
right_threshold_order
lu_cut_order
```

It must not encode factorized source topology, generated local E-surface IDs, or
any reduced-source graph identity.

## Boundary invariant

Inside `map_to_3drep`, it is allowed to build reduced/factorized source graphs,
generate local source CFF/LTD expressions, and use source-local E-surface IDs.

Outside `map_to_3drep`, none of that source-local information should leak.

The only exposed result should be:

```rust
BTreeMap<CutCFFIndex, Atom>
```

with keys generated from the caller/top-level `CutSet` convention.

## Top-level graph mutability

The local 4D projection does not fundamentally require mutating the top-level
graph. Current code may require `&mut Graph` because existing 3Drep/CFF
conversion helpers populate or consult graph-owned surface caches.

This is an API artifact, not a mathematical requirement. It is acceptable for an
initial refactor to keep `&mut Graph`, but the long-term target is:

```rust
fn map_to_3drep(graph: &Graph, ...) -> Result<BTreeMap<CutCFFIndex, Atom>>
```

with any mutable state confined to local source-projection builders.

## Combination rule

The factorized graph contribution should be inserted multiplicatively before
projection:

```text
local_4d_forest_atom * factorized_cograph_atom
```

Then this product is projected once:

```text
map_to_3drep(local_4d_forest_atom * factorized_cograph_atom)
```

The preferred construction is not:

```text
map_to_3drep(local_4d_forest_atom) * map_to_3drep(factorized_cograph_atom)
```

because independently projecting the two factors can put them in different local
surface bases and can hide or introduce residue-selection mismatches.

For finite integrated subtraction terms, project each already assembled product
into the same top-level `CutCFFIndex` key space and combine additively key by
key.

## Practical implementation target

1. Add `Local4DApproximation` with a kernel structurally matching
   `Integrated<'_>` but without Vakint integration.
2. Store the recursively assembled local 4D atom in the UV orchestrator.
3. Move expanded-4D CFF/LTD/3Drep projection into a finalization function,
   conceptually `map_to_3drep`.
4. Ensure `map_to_3drep` returns a top-level `BTreeMap<CutCFFIndex, Atom>`.
5. Keep source-local graph, source-local E-surface, and generated 3Drep details
   private to `map_to_3drep`.
