# UV Renormalization Architecture

## Scope

This document records the implemented UV-renormalization invariants shared by
the legacy DAG-forest and hedge-poset orchestrators. Historical revision
comparisons, temporary log locations, and external-reference convention tables
belong in tests or investigation records rather than in the current
architecture.

## Orchestrators

`UVgenerationSettings::orchestrator` selects one of three execution modes:

- `legacy_dag_forest` computes the established `Approximation` DAG and remains
  the reference implementation for connected graphs.
- `hedge_poset` unfolds commuting UV operations into trace levels and stores
  their results in a compute store. It supports disconnected unions.
- `compare` computes both connected-graph backends, compares normalized
  expressions, and returns the legacy result. It is a validation mode, not a
  request to execute the returned integrand with the hedge-poset backend.

The backends share analytic operations and result types. Scheduling,
dependency lookup, caching, and disconnected-component composition belong to
the orchestrator because their traversal models differ.

## Scheme and Computation Ownership

The renormalization prescription belongs to the `Spinney`. A compute node does
not copy that policy; it stores computed values:

- the signed local four-dimensional counterterm;
- the integrated counterterm projections;
- cut-dependent local three-dimensional and final integrands.

This separation lets the same integrated result be consumed differently
without storing scheme-specific copies.

## Signs and Integrated Projections

Each UV operation supplies its own subtraction sign. In particular, the local
four-dimensional limit constructs the next counterterm as a signed `-T(...)`
operation. Nested signs therefore arise from recursive operation composition,
not from terminal depth or parity.

For a connected counterterm, integration stores a Laurent series from which two
physical projections are obtained:

```text
pole projection:               pole(integral(local))
finite-counterterm projection: -finite(integral(local))
```

No extra terminal parity is applied. Adding a sign derived from operation count
would duplicate signs already supplied by the local operations.

The projections are used as follows:

- a non-root `PolePart` dependency recurses through the completed integrated
  pole only;
- a root `PolePart` result combines its local term with the integrated pole;
- an `MUV` dependency combines its local term with the signed finite
  counterterm;
- terminal renormalization output selects the pole or finite projection from
  the source Spinney's prescription;
- final three-dimensional integrand assembly always localizes the signed finite
  counterterm, regardless of the terminal prescription.

The last distinction is essential: changing terminal output policy must not
remove the finite addback required when assembling an integrable 3D expression.

## Disconnected Composition

Disconnected composition depends on which object is being combined.

### Four-dimensional local terms

Each connected component first obtains its complete recursive counterterm,
including the projection selected by that component's Spinney. The union's
local counterterm is the product of those full component counterterms.

### Integrated terms

Integrated counterterms are scalars and factorize over connected components.
The stored disconnected value provides the product of component pole
projections and the product of component finite-counterterm projections.

The aggregate is sufficient when every component uses the same terminal
prescription. At terminal use, the hedge-poset backend nevertheless returns to
the component nodes, selects each component's own pole or finite projection,
and only then multiplies them. This also handles mixed prescriptions without
storing scheme-specific variants of the aggregate.

### Three-dimensional local terms

Cut CFF structures need not factorize over disconnected components. Multiplying
complete per-component local results would also multiply the common root and
can duplicate or incorrectly separate intertwined CFF structure.

The hedge-poset backend therefore replays component-local operation paths from
their common root. At a union it:

1. obtains the replay states for every component;
2. forms compatible combinations of their integrated prefixes;
3. joins those prefixes in the trace/Foata representation;
4. applies each component's local operations to the shared incoming state.

The construction operates over an arbitrary number of components and is not
special-cased for a two-component spectacles graph.

## Marker Representation

UV markers are diagnostic symbolic structure. Approximation, integration,
series, and truncation are distinct canonical function heads applied to the
history of current/given subgraphs. Spenso printing supplies the compact
`K[...]`, angled-bracket, `Σ(...)`, and truncation notation while structured DOT
uses the same underlying atoms. Typst presentation combines Symbolica's Typst
arithmetic mode with those Spenso shorthands; truncation is emitted as
`op("Tr")(…)`, which is valid Typst math, while compact terminal output remains
`Tr(…)`.

DOT `full_num` fields contain Typst math fragments because the drawing template
evaluates them in math mode. Per-graph renormalization `.typ` files wrap the same
fragment in a display-math document, while the accompanying `.txt` files remain
the round-trippable Symbolica representation.

Subgraph labels must be created through the subgraph symbol factory so every
consumer receives the same symbol metadata and print behavior.

## Current Boundaries

- Parametric integrand generation currently supports final 3D output only.
  `FourD` is used by integrated renormalization internally but is rejected by
  the parametric orchestrator.
- The legacy DAG executor constructs union-shaped nodes but does not execute
  multi-parent unions. Disconnected generation therefore requires the
  hedge-poset backend.
- `compare` cannot validate disconnected unions until the legacy backend has a
  corresponding execution path.
- External systems such as RQFT may use different propagator or vertex sign
  conventions. Those graph-specific conversion factors are recorded next to
  the reference expressions in renormalization tests; they are not part of the
  internal recursive sign convention.
