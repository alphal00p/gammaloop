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
- `hedge_poset` is the default. It unfolds commuting UV operations into trace
  levels and stores their results in a compute store. It supports disconnected
  unions.
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

### Direct three-dimensional local terms

Both direct local-3D representations perform the complete loop-energy
integration first and apply every UV Taylor operator to the resulting
complete/global CFF expression. They therefore share exactly the same local
Taylor-transformed CFF bodies. If the generalized residue map is
`{ k -> C_k }`, the orientation-parametric form is represented as
`sum_k sigma(k) C_k`, where `sigma(k)` selects the complete residue-map key
(loop map, edge map, numerator sampling map, and physical direction metadata).
The Taylor operators treat `sigma(k)` as opaque. The sparse implementation
therefore applies them branchwise without expanding the factorized numerator.
`explicit_orientation_sum_only=true` merely replaces every `sigma(k)` by one
and explicitly sums the same bodies. Neither direct form
uses the exact-source reconstruction or minimax dispatch of completed local-4D
Taylor terms.

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

The projected local-4D route is separate: it completes the Taylor expansion in
four dimensions, reconstructs the exact source occurrence graph, performs the
factorized minimax EMR dispatch needed for derivative-created occurrences, and
only then projects that completed term to CFF. That reconstruction machinery is
exclusive to projected local4D and is not a replacement for the direct replay
above.

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

## UV profiling from the CLI

`profile ultra-violet` probes the generated integrand along rays on which one
or more loop momenta are scaled to infinity. It works for both amplitudes and
Local-Unitarity cross sections. The ordinary, fast invocation is:

```text
profile ultra-violet -p <process> -i <integrand>
```

By default, `--selected-limits only-divergent` is active. For each graph, this
mode uses the graph's expected bare UV degrees to retain only cycle unions with
degree of divergence greater than or equal to zero. It profiles each distinct
cycle union once in one generated physical loop-momentum basis, preferring the
generation LMB when that basis represents the cycle. Thus the default avoids
the old Cartesian scan over every loop subset in every LMB while still testing
every expected divergent cycle. The target cycles are enumerated on the same
physical domain as UV generation: dummy edges are excluded for amplitudes; LU
first removes the incoming side and then enumerates separately in the
complement of every compatible exact physical cut. If the stored physical LMBs
cannot represent an expected divergent cycle exactly, profiling reports an
error instead of silently claiming complete coverage.

The exhaustive behavior remains available explicitly:

```text
profile ultra-violet -p <process> -i <integrand> \
  --selected-limits all
```

`all` visits every nonempty loop-coordinate subset of every generated physical
LMB. The same graph-theoretic cycle may consequently appear in several rows;
that repetition is intentional because this mode audits all LMB
representations. Cutkosky cuts do not create duplicate rows: in LU mode they
are compatibility alternatives for a limit, not an additional profiling
dimension.

One graph can be selected by its name, numeric ID, or `#`-prefixed numeric ID:

```text
profile ultra-violet -p <process> -i <integrand> --graph GL0
profile ultra-violet -p <process> -i <integrand> --graph '#2'
```

For a cross section, a particular physical Cutkosky cut can additionally be
selected by its edge IDs. A graph selection is required:

```text
profile ultra-violet -p <process> -i <integrand> \
  --graph GL0 --cutkosky-cut 2,5
```

`--cut-edges` is a visible alias for `--cutkosky-cut`. The edge order is
immaterial. GammaLoop rejects a missing, ambiguous, or amplitude-mode cut
selection and reports the available physical cuts for a missing match.

### Physical cuts and LU residue groups

Two identities must remain distinct when propagators are raised. Several
physical Cutkosky cuts can correspond to the same normalized E-surface
residue. Production LU evaluates that residue group once and emits an event
tagged with the group's deterministic first cut. Therefore selecting any
non-representative physical cut in that group correctly projects the weight of
the representative event; evaluating another independent residue would
double-count the group.

The selected physical cut nevertheless remains authoritative for deciding
which UV limits exist. A cycle is compatible with a cut precisely when it does
not contain an edge of that cut. Consequently:

- with `--cutkosky-cut` (or `--cut-edges`), compatibility is tested against
  exactly the requested physical edge set, even if another raised-edge alias
  belongs to the same residue group;
- without a cut selection, an LU limit is retained when it is compatible with
  at least one exact physical cut of the graph.

The second rule applies in both selected-limit modes. In exhaustive `all` mode
it means that every nonempty LMB subset which is physically compatible with at
least one cut is present, while the subset is still evaluated only once as
part of the graph's summed LU integrand. In `only-divergent` mode the same
compatibility filter is applied before retaining the unique expected-divergent
cycle representatives.

### Result and failure summaries

The command first prints one detailed fitted-slope table per selected graph.
The final `UV limit tests` line reports the passed and total counts. If any
limit fails, it is followed by one colored rounded table containing only the
failures, with columns for graph, selected cut edges, LMB, fixed and scaled
loop edges, orientation, and failure reason. This final table is intended to
make failures actionable without searching all successful rows.
`--per-orientation` is available only for the orientation-parametric, localized
direct-local3D representation; it includes orientation-local fits and failures
in addition to the summed fit. Selector-free explicit-sum direct local3D and
projected local4D are summed representations and reject this option rather than
claiming per-production-orientation finiteness.

Each finite-precision ray is first fitted over its complete scale range. A
missing or non-linear double-precision fit triggers a complete Arb reevaluation
of that ray. The reported power is then fitted on the longest high-scale suffix
with at least half of the samples (and at least five) whose
\(R^2\) is at least `0.99`; this permits a short pre-asymptotic transient without
discarding isolated points. A fit which still has no such asymptotic suffix is
reported as `unstable_fit`, even when its raw slope happens to lie below the UV
threshold. At least five scale points are consequently required.
The minimum and maximum scale exponents must also be finite and strictly
increasing, so the fitted suffix is always the large-scale end of the ray.

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
