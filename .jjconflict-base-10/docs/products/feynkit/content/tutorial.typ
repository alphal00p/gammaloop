#import "../../shared.typ": callout, source-link

#let tutorial = [
= From a model to finalized diagrams

Begin with a normalized JSON model, as in the quickstarts. Model validation resolves particle,
parameter, interaction, propagator, and Lorentz-structure references to stable IDs. A diagram
retains the same model so particle and interaction IDs remain meaningful after generation.

== Import a UFO model

Raw UFO import is an optional adapter around Python's `ufo_model_loader`. Install that package
in the interpreter used by the Symbolica host, then point the loader at your UFO directory:

// docs-example: compile
```python
from symbolica.community.feynkit import UfoLoader

loaded = UfoLoader(restriction_name="massless").load("models/sm")
model = loaded.model
```

The directory and restriction are model-specific inputs. Rust callers enable `feynkit/ufo`
and call `UfoLoader` with their own attached interpreter; the adapter does not initialize Python
or replace its logging and standard streams.

== Inspect and configure before generating

`Model` exposes particles, parameters, couplings, vertex rules, Lorentz structures, propagators,
functions, and form factors through collections and lookups. `particle.antiparticle` resolves the
paired record in the same model. `ParameterCard` changes external inputs; derived expression
evaluation is an explicit boundary. Python callers can pass an evaluator to
`with_parameter_card()` or use `recompute_with()`. The evaluator receives `EvaluationRequest`
and returns `EvaluatedValues`; incomplete results raise `ModelError` without publishing a
partially updated model.

Python accepts generation settings directly as keyword arguments. Exact coupling orders use
integers; pairs specify inclusive bounds, with `None` for an unbounded coupling maximum.
Filter objects specify which subgraphs to reject. Python's defaults are ported from the
GammaLoop CLI: self-loops are permitted, zero-flow edges are rejected, and non-vacuum
processes filter self-energies, tadpoles and zero-momentum snails. Vacuum processes leave
those filters off and require one factorized loop topology. Cross sections default to one
cut blob, no spectators, and symmetrized final states.

`maximum_bridges=0` is the Python default for every process. Use `maximum_bridges=None`
to include unrestricted bridge topologies, such as tree-level exchange channels.
Omitting `numerator_grouping` groups up to scalar rescaling; explicit `numerator_grouping=None`
disables numerator comparison. Diagrams still contain their vertex, propagator and aggregate
numerators. As in the current GammaLoop `--only-diagrams` path, generation does not construct
an evaluable integrand. Numerator validation first checks identical factored expressions;
it expands only when needed to compare differing representations.

For process-dependent filters and grouping, `...` in the signature means automatic defaults.
Pass an options object to customize a filter, or `None` to disable it. The same distinction
applies to factorized-loop, cut-blob and spectator ranges.

// docs-example: compile
```python
import symbolica.community.feynkit as fk

result = model.generate_diagrams(
    incoming=["g"],
    outgoing=["g"],
    loops=1,
    coupling_orders={"QCD": 2, "QED": 0},
    particle_veto=["c", "t", "s", "u", "d"],
    zero_snails=fk.SnailFilterOptions(veto_attached_to_massless=True),
    threads=4,
)
```

Reuse keyword settings with a Python dictionary and `**kwargs`; each call constructs a fresh
Rust configuration. A process can also be configured independently of generation:

// docs-example: compile
```python
import symbolica.community.feynkit as fk

process = fk.Process.amplitude(["e-", "e+"], ["mu-", "mu+"])
process = process.with_loop_count(1, 1)
result = fk.Generator(model).generate(process, max_vertices=6, maximum_bridges=None)
```

Here `model` must provide those particles. Loop range, vertex bounds, allowed interactions,
particle vetoes, coupling orders, and numerator grouping are explicit generation choices. An
empty result can therefore mean the chosen constraints admit no graph. This is a completed
`GenerationResult` with zero retained diagrams, not a configuration error. Inspect the generation
report and relax a specific constraint before enlarging the search indiscriminately.

Both generation entry points check Python signals while they run. Ctrl-C or a notebook
interrupt stops generation and raises `KeyboardInterrupt`; exceptions from custom Python
signal handlers also propagate. Cancelling an explicit `CancellationToken` instead returns
an incomplete result, preserving that API's existing partial-result behavior.

== Observe progress and prune partial topologies

Both generation methods default to `progress="auto"`. In Marimo, they check
`marimo.running_in_notebook()` and maintain one spinner showing the current stage,
processed counts (including a total when known), and elapsed time. The spinner closes
on completion, cancellation, or error. Outside a running notebook, automatic progress
stays quiet and does not import Marimo. Pass `progress=None` to disable it.

For custom reporting, pass `progress=callback`. The callback receives an immutable
`GenerationProgress` with `stage`, `completed`, and `total`. Counts restart at each stage;
`total is None` means the amount of work is unknown, so show a spinner instead of a
percentage. Topology enumeration counts discovered unique graphs, not every explored branch.
Later stages count processed inputs, which may be rejected. `complete` reports the final
retained count; `cancelled` marks a partial result. The stages are `topologies`,
`topology_filters`, `interactions`, `interaction_filters`, `numerators`, `selection`,
`grouping`, and the terminal stage.

Native workers publish snapshots; Python callbacks run on the thread that called generation,
so Marimo keeps its cell context. Updates within a stage are coalesced to about 150 ms;
stage transitions and completion are delivered promptly. Browser kernels invoke callbacks on
their calling thread. Callback delivery itself does not force a browser paint or UI refresh.

// docs-example: compile
```python
result = model.generate_diagrams(
    ["e-", "e+"], ["mu-", "mu+"], loops=1, maximum_bridges=None
)
```

Use `filter=callback` for live pruning during enumeration. It receives a Symbolica `Graph`
snapshot and `completed_vertices`: the first N vertices have all their incident edges fixed,
while later vertices may still gain edges. Returning `False` rejects that search branch,
without cancelling the run. Reject only conditions that cannot become acceptable as the
partial graph grows; a numerator or UV-degree test requires a finished `FeynmanDiagram`.

// docs-example: compile
```python
def no_self_edges(topology, completed_vertices):
    return all(source != target for source, target, directed, pdg in topology.edges())

result = model.generate_diagrams(
    ["e-", "e+"], ["mu-", "mu+"], loops=1,
    allow_self_loops=True, maximum_bridges=None, filter=no_self_edges,
)
```

This example illustrates the callback contract; for this particular condition, prefer the
built-in `allow_self_loops=False` option to avoid Python callback overhead.

Node data is a Symbolica integer: 0 for an interaction vertex, `-(index+1)` for an incoming
external leg, and `+(index+1)` for an outgoing leg, using the generator's external-leg index.
Edge data is the base particle's PDG code; directed edges distinguish particle flow from
antiparticle flow. These snapshots can be retained or modified without changing generation.
Keep the filter cheap: it runs synchronously for enumeration candidates, and the current
Symbolica API requires building each Python snapshot through `Graph.add_node` and
`Graph.add_edge`. No FeynKit topology wrapper or independent graph algorithms are involved.
Exceptions from either callback cancel generation and propagate unchanged.

== Analyze interaction regions with Linnet

Finalized diagrams use GammaLoop's half-edge conventions. Amplitude external states are
dangling half-edges, incoming at a sink and outgoing at a source. Cross-section initial
states are sewn paired edges; their `is_external` metadata still identifies external
momentum carriers. Every item in `diagram.vertices` is an interaction. A missing edge
endpoint is `None`, and external names, indices, and states belong to the edge.

Install the matching `linnet-py` extension to use the graph analysis interface.
`to_linnet()` returns that module's canonical `Graph`, with `DiagramVertex` and
`DiagramEdge` payloads. Its selection algebra and graph algorithms therefore work directly:

// docs-example: compile
```python
graph = diagram.to_linnet()
selected = graph.filter(edge=lambda edge: not edge.data.is_external)
numerator = diagram.numerator_expression(subgraph=selected)
components = diagram.connected_components(selected)
boundary = diagram.boundary(selected)
basis = diagram.momentum_basis(subgraph=selected)
routed = basis.route_expression(numerator)
```

Optional `subgraph` filters also apply to denominators, momentum-basis enumeration,
parent-compatible and contracted routing, superficial divergence, CFF construction, and
tensor reduction. Numerator `without` excludes edges and their incident vertex factors,
using the shared GammaLoop traversal. Full-graph defaults retain whole-diagram behavior.
Partial tensor reduction requires an explicit `projector`; a boundary momentum is never
implicitly integrated. Overall factors, numerator prefactors, and projectors remain separate.

Selections belong to one exported graph revision. A structural edit to that graph leaves
the physics diagram unchanged and invalidates its selections. The next `to_linnet()` call
creates a fresh analysis graph. Foreign or stale selections raise an error. Half-edge
views carry their native diagram half-edge ID in `data`; the bridge preserves this mapping
even when the exported builder assigns different half-edge numbers.

`diagram.cuts` contains generated physical cuts with reusable `left.subgraph` and
`right.subgraph`, coupling orders, side loop counts, crossing edges, and momentum signatures.
`diagram.topology_threshold_candidates` records the distinct process-independent partitions.
Linnet's `all_bonds` and `all_cuts` describe structural cuts, without assigning physical
final-state validity.

== Inspect the propagator denominator

// docs-example: compile
```python
denominator = diagram.denominator_expression()
integrand = diagram.numerator_expression() / denominator
```

The denominator is a scalar Spenso `TensorExpression` using GammaLoop's
`denom(edge, momentum, mass_squared, quadratic)` annotation. The quadratic is
$q_e^2 - m_e^2$, with the shared symbolic dimension by default; pass `dimension=4`
for a fixed dimension. `edge_powers` maps edge IDs to signed powers, including zero.
Momentum labels match the numerator, masses remain symbolic model parameters, and the
UFO `ZERO` mass becomes zero. The default region excludes external carriers and dummy
edges; an explicit selection includes every paired edge, matching GammaLoop. Widths and
an imaginary prescription are excluded. Custom UFO denominator formulas are not instantiated.
A selection without internal propagators returns one.
The ratio above still requires the diagram's separate overall factor and any unapplied
projector or numerator prefactor when assembling a complete integrand.

== Inspect superficial UV power counting

// docs-example: compile
```python
degree = diagram.superficial_degree_of_divergence()
degree_in_six_dimensions = diagram.superficial_degree_of_divergence(dimension=6)
```

The local superficial degree adds the loop integration measures, momentum powers in the
stored vertex and internal-edge numerators, and minus two for each quadratic propagator
denominator. External legs, projectors, and global prefactors do not enter this count.
Zero denotes a logarithmic superficial divergence, a positive value a power divergence,
and a negative value superficial convergence. This local bound scales all momenta at each
vertex together; it does not account for tensor cancellations or rule out UV subdivergences.

== Use finalized output

Generation performs topology expansion, interaction assignment, filters, canonicalization,
fermion flow, numerator and projector construction, routing, zero detection, and grouping.
The returned diagrams retain vertex/edge numerator fragments, the aggregate numerator,
external projector, and separate scalar prefactors. Keep group multiplicities and prefactors
when forming a complete expression; a diagram's numerator alone is not the whole amplitude.

Named coefficients such as `UFO::GC_11` stay in generated numerators. For analytic work, expand
them using the model that owns the diagram:

// docs-example: compile
```python
analytic_numerator = model.expand_couplings(diagram.numerator_expression())
```

This returns a new expression and preserves the stored named coefficients. Serialize with
`to_json()` or `to_dot()` and restore with `FeynmanDiagram.from_json(model, text)` or
`from_dot(model, text)`. The model fingerprint and structural validation protect the
interpretation of stored IDs.

#callout("CFF IDs belong to their arena", [
  `diagram.build_cff()` returns orientations, surfaces, generation statistics, and symbolic
  denominators. Rust `generate()` creates a fresh arena; `generate_into()` accepts a shared one.
  Resolve only the surface IDs referenced by an expression through its associated arena.
])

The #link("reference/interfaces/")[API map] identifies the owning crate for each stage.
For executable source examples, see
#source-link("crates/feynkit/examples/generate.rs", label: "the Rust generation program") and
#source-link("crates/feynkit-py/examples/ufo_generation.py", label: "the Python UFO program").
]
