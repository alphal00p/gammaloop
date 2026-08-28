# CFF Surface-Cache Ownership

## Status

Implemented. `feynkit-cff` owns the only structural CFF expression, surface,
identifier, tree, and cache types used by FeynKit and GammaLoop.

## Invariant

A surface identifier is meaningful only with the `SurfaceCache` supplied with
the expression that contains it. Consumers must discover an expression's
surface set from its trees; they must not treat every entry in a shared arena as
part of that expression.

```text
CffResult
  expression -> trees containing SurfaceId values
  surfaces   -> arena that resolves those IDs
  report     -> generation statistics
```

The arena may contain additional entries when related expressions deliberately
share it. IDs are append-only and stable, while `CffExpression` exposes the
exact referenced subset. Raised-surface analysis, residue selection, symbolic
lowering, threshold construction, and diagnostics therefore start from the
expression's IDs and resolve only those IDs through its arena.

## Generation Modes

`CffGenerator::generate` creates a fresh arena and returns it in `CffResult`.
This is the normal standalone and Python API.

`CffGenerator::generate_into` accepts a caller-owned `SurfaceCache` and returns
an expression whose IDs use that arena. GammaLoop uses this mode for related
root, contracted, subgraph, and UV expressions that must share stable IDs. The
runtime graph stores the canonical FeynKit cache; it does not define a second
surface cache or translate through GammaLoop-specific surface IDs.

When independently generated results must be combined, `SurfaceCache::merge`
returns a `SurfaceIdMap`, and `CffExpression::remap_surfaces` applies that map to
the expression. This makes arena changes explicit instead of relying on index
coincidence.

## Downstream Extensions

GammaLoop adds numerical and Symbolica-specific behavior with extension traits
implemented directly for FeynKit CFF types. These extensions may evaluate a
surface, lower it to a runtime expression, classify thresholds, or construct
subtraction data, but they do not wrap or copy the structural CFF IR.

The topology conversion is implemented once by `HedgeGraphCffExt`, directly on
Linnet graphs. GammaLoop classifies runtime edges as standard, omitted dummy, or
sewn initial-state edges through that extension trait; it does not build an
adapter object or reconstruct a `CffGraph` itself. It then consumes the returned
FeynKit expression and shared cache directly. There is no separate UV CFF
engine.

## Required Checks

Tests and the ownership-boundary guard enforce these properties:

- every surface ID in an expression resolves in the associated arena;
- generating another expression into a shared arena never changes existing
  IDs or the referenced surface subset of an earlier expression;
- raised-surface groups contain only IDs present in the source expression;
- merging arenas and remapping preserves symbolic lowering;
- root, contracted, subgraph, and UV call paths all use `feynkit-cff`;
- GammaLoop does not define compatibility aliases or duplicate CFF structures.
