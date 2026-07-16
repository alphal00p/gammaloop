# CFF Surface-Cache Ownership Proposal

## Status

This document proposes an ownership change. It describes a target architecture,
not the current implementation.

## Problem

CFF generation discovers energy and hybrid surfaces while recursively building
orientation trees. The tree stores compact `EsurfaceID` and `HsurfaceID` values,
so generation needs a mutable interning table to deduplicate each discovered
surface and assign its ID.

The current `Graph::generate_cff` implementation uses `Graph::surface_cache` as
that table. It then clones the complete cache into the returned
`CFFExpression::surfaces`. A later CFF generation on the same graph, including
one requested while constructing local UV counterterms, can append more
surfaces to the graph cache.

Appending does not renumber existing IDs. The resulting problem is instead an
ownership mismatch: the graph cache no longer denotes exactly the surfaces of
an earlier CFF expression. Code that combines an earlier expression with the
current graph cache can consequently include surfaces that never occur in the
expression.

For example, raised-surface analysis currently receives a `CFFExpression` but
enumerates `Graph::surface_cache`:

```text
global CFF:       tree IDs {0, 1}, expression cache [eta0, eta1]
later UV CFF:     graph cache becomes [eta0, eta1, eta2, eta3]
raised analysis:  groups IDs {0, 1, 2, 3} against a tree using {0, 1}
```

The old IDs remain valid as indices into the append-only graph cache, but the
extra IDs do not belong to the global expression. Correctness then depends on
when the analysis runs, which is a fragile pipeline-ordering constraint.

## Ownership Invariant

A surface ID should be meaningful only within the `SurfaceCache` owned by the
same `CFFExpression`:

```text
CFFExpression = orientation trees + their surface-ID arena
```

The expression and arena should be created together and remain immutable after
generation. A consumer of an `EsurfaceID` or `HsurfaceID` from an expression
must resolve it through that expression's cache, never through a cache on the
source graph or another expression.

## Proposed Design

Use one generation-local `SurfaceCache` for each CFF expression. All
orientations produced by that generation share the local cache, preserving
deduplication and coherent IDs within the expression. When generation finishes,
move the cache into `CFFExpression` rather than retaining and cloning a mutable
copy on `Graph`.

Conceptually, the API should have this shape:

```rust
impl Graph {
    fn generate_cff(&self, /* generation options */) -> Result<CFFExpression<OrientationID>> {
        let mut surfaces = SurfaceCache::new();
        // Generate every orientation using `surfaces` as the interning arena.
        // Return the trees together with `surfaces`.
    }
}
```

The mutation remains necessary during construction, but it is confined to the
new result. Calling CFF generation no longer mutates unrelated graph state, and
previously generated expressions cannot be affected by later generations.

`Graph::surface_cache` should be removed if an audit finds no independent graph
invariant that requires it. Cross-expression surface deduplication alone is not
a sufficient reason to retain it: expression-local IDs are small, generation
already needs to visit the surfaces, and shared interning introduces an
ownership protocol that every consumer must obey.

## Immediate Correction

Before the full ownership refactor, consumers should consistently use the
snapshot already carried by `CFFExpression`.

In particular:

- `determine_raised_esurfaces_from_expression` should enumerate
  `expr.surfaces.esurface_cache`, not `self.surface_cache.esurface_cache`.
- CFF surface substitution should use the cache belonging to the expression
  being substituted.
- Runtime amplitude data should retain that expression cache and use it for
  threshold counterterms, overlap checks, diagnostics, and sampling.
- Grouped derived data should resolve raised IDs through each graph term's CFF
  expression cache.

This removes dependence on whether raised-surface analysis happens before or
after UV integrand construction. Moving the analysis earlier can avoid the
current symptom, but it does not establish the ownership invariant.

## Migration

1. Change the graph-level CFF entry point to create a fresh `SurfaceCache` and
   share it only across orientations in that call.
2. Move the completed cache into the returned `CFFExpression`; avoid cloning it
   solely to preserve graph-global state.
3. Replace every lookup through `Graph::surface_cache` with a lookup through
   the `CFFExpression` or runtime graph term that supplied the ID.
4. Adapt residue selection and CFF-to-atom substitution so transformed
   expressions retain and use their original cache.
5. Remove `Graph::surface_cache` after its remaining uses have been audited and
   migrated.

The subgraph and UV CFF entry points should follow the same rule. Some already
construct a local cache; those paths can provide the reference behavior for the
graph-level API.

## Validation

The refactor should include tests for these invariants:

- Generating a second, different CFF on a graph does not change the surface set
  or raised-surface data of the first expression.
- Every surface ID appearing in an expression is in bounds for that
  expression's cache.
- Raised-surface groups contain only surfaces represented by the expression.
- CFF substitution gives the same atom before and after another CFF is
  generated from the same graph.
- Threshold construction and evaluation use the same expression-owned cache.

## Rejected Alternative

A graph-global arena can be made correct if IDs explicitly carry an arena or
generation identity and every expression records its referenced subset. That
design preserves cross-expression interning but adds typed-session and subset
bookkeeping throughout the CFF and threshold APIs. It should only be preferred
if profiling demonstrates that cross-expression deduplication is materially
valuable. No such requirement is currently established.
