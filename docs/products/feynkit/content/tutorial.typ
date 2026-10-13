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
from symbolica.community.hepkit import UfoLoader

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

Create an immutable request from `model.process(...)`. Every physics request requires
explicit finite loops: `fk.Loops.total(0)` or `fk.Loops.total((0, 1))` for amplitudes,
including vacuum amplitudes; `fk.Loops.forward(...)` counts sewn forward-graph loops.
Use `Orders().exact("QCD", 2)` or inclusive `Range` values for coupling constraints.
Repeated bounds intersect; conflicting bounds describe an empty request.

`Restrict` owns particle and vertex clauses. `TopologyRules` owns bridge and cyclic-subgraph
constraints; `Cuts` owns physical cut content. These defaults are shared and unrestricted:
there is no automatic non-vacuum or vacuum preset. Self-loops are allowed, zero-flow edges
are rejected, and all external symmetrizations are disabled. Optional count ranges default
to unrestricted, including factorized loops, cut blobs, and spectators.

The shared `graph.Bridges` constructors distinguish `any()`, `at_most(n)`, `one_pi()`, and
`one_pi_in_loops()`. The last retains tree exchange channels and requires one-particle
irreducibility only for graphs with loops, even in a mixed loop interval. External legs
do not count as bridges. A bridge limit on a cross-section request acts on the forward
graph and can remove physical s-channel contributions, so choose it explicitly.

This example preserves a filtered one-loop gluon two-point search by specifying the
cyclic vetoes and loop-only bridge policy that it needs:

// docs-example: compile
```python
from symbolica.community import graph
import symbolica.community.hepkit as fk

request = model.process(["g"], ["g"]).amplitude(
    loops=fk.Loops.total(1),
    restrict=fk.Restrict().exclude_particles(["c", "t", "s", "u", "d"]),
    orders=fk.Orders().exact("QCD", 2).exact("QED", 0),
    topology=fk.TopologyRules(
        bridges=graph.Bridges.one_pi_in_loops(),
        self_energies=fk.Veto.all(),
        tadpoles=fk.Veto.all(),
        snails=fk.Veto.massless(),
    ),
)
compact = request.generate(run=fk.Run(threads=4))
result = compact.build()
```

Reuse the specification objects directly. A request can explain or count its discovery
without constructing Feynman diagrams or numerators:

// docs-example: compile
```python
process = model.process(["e-", "e+"], ["mu-", "mu+"])
request = process.amplitude(
    loops=fk.Loops.total(1),
    topology=fk.TopologyRules(
        max_vertices=6,
        self_energies=fk.Veto.all(),
        tadpoles=fk.Veto.all(),
        snails=fk.Veto.massless(),
    ),
)
print(request.explain())
counts = request.count()
assert counts.complete
compact = request.generate()
assert len(compact) == counts.diagrams
result = compact.build()
```

Counts distinguish skeletons, complete dressings, and physically retained diagrams.
`sum_inv_aut` is an exact `fractions.Fraction`; `physical_weight` also includes signed
process factors. `per_loop` holds the exact totals at each loop order. Compact results
own their model, prepared assignments, and shared skeletons independently of the request.
Accessing their metadata or iterating does not construct numerators, names, or diagram IDs.

The same process can produce a symbolic amplitude or sewn cross-section diagrams:

// docs-example: compile
```python
process = model.process(["e-", "e+"], ["a", "a"])
amplitude = process.amplitude(loops=fk.Loops.total(0)).generate().build().into_amplitude()
squared = amplitude.squared().sum_spins(average_initial=True)
forward = process.cross_section(
    loops=fk.Loops.forward(1),
    topology=fk.TopologyRules(
        self_energies=fk.Veto.all(),
        tadpoles=fk.Veto.all(),
        snails=fk.Veto.massless(),
    ),
    cuts=fk.Cuts(blobs=fk.Range.exact(1), spectators=fk.Range.exact(0)),
).generate().build()
```

The cross-section result contains sewn forward diagrams and cut metadata; phase-space
integration is a later operation. A two-particle tree cut requires one forward loop.
The example explicitly preserves a one-blob, zero-spectator cut policy; these restrictions
are not defaults. `into_amplitude()` rejects incomplete or empty output rather than
presenting an incomplete coherent sum. It also rejects already-sewn forward diagrams.

Identical final-state particles of a cross section are cut edges, already counted by graph
automorphisms, so cross-section requests reject `Symmetrize.FINAL`.
`Symmetrize.INITIAL` retains one graph per exchange of identical beams and multiplies its
overall factor by the represented orbit size. `Cuts(symmetrize_left_right=True)` identifies
forward graphs related by a global exchange of the amplitude sides and sums their
multiplicities. Cut enumeration is independent of the mirror frame, including when
combined with `Symmetrize.INITIAL`.
Evaluating only retained graphs is exact after integration when observables and selectors
respect beam exchange and beam helicities are summed or equal; exchange-fixed graphs are
exact in every event. Fewer graphs need not lower the integration error: for
$Z Z -> e^+ e^-$ it was 19 times larger at equal sample counts. For amplitudes,
`Symmetrize.INITIAL` and `Symmetrize.FINAL` make the corresponding identical legs unlabelled,
and the caller sums their permutations.

External states are validated in `Model.process`. Request restrictions resolve against the
same model: allow clauses intersect, exclusions unite, and `.only_vertices([])` permits no
interactions. A completed empty discovery can mean no graph satisfies the constraints;
inspect `request.explain()` and the result report before enlarging the search. A cancelled
count has `complete=False`; cancelled discovered or built output has `report.completed=False`.

Names are absent by default. Build with `fk.NumeratorOptions(naming=fk.Naming.named("D"))`
when names are required; the no-argument named constructor uses prefix `GL`. Name selection
and named momentum-basis overrides require that build demand. Stable IDs are computed only
when requested. Selector strings of exactly 32 hexadecimal digits denote IDs; other strings
denote names. `built.select(fk.Selection(select=[]))` keeps no diagrams. Optional numerator
comparison is a later `built.group(fk.Grouping("identical"))` operation, not a discovery
setting. Grouping checks equal factored expressions first and expands only when comparison
requires it. Graph generation and construction do not create an evaluable integrand.

Pass `Run` to each execution stage, including `built.select(selection, run=run)`, for
threads, cancellation, progress, and hereditary checks. Selection reports the number of
input diagrams scanned, including vetoed diagrams. Cancellation returns an incomplete
subset; unmatched forced cuts are diagnosed only after a complete scan. Python signals are checked during execution: Ctrl-C raises `KeyboardInterrupt`, and
custom signal-handler exceptions propagate. An explicit `CancellationToken` instead returns
an incomplete result. Request criteria and stored compact discoveries can be reused with
independent runs.

== Observe progress and prune partial topologies

Progress is disabled by default. Pass `fk.Run(progress="auto")` for a Marimo status display
or `fk.Run(progress=callback)` for custom reporting. Automatic progress checks
`marimo.running_in_notebook()` and stays quiet outside a running notebook. Known totals use
a progress bar; unknown or empty totals use a spinner. The display closes on completion,
cancellation, or error.

A custom callback receives an immutable `GenerationProgress` with typed `Stage`, exact
`completed`, optional `total`, and `zero_numerator_count` values. Return `True` to continue
or `False` to cancel. Counters restart at each stage: `Skeletons` counts complete untyped
topologies, `Dressing` counts complete interaction assignments, and `Prepare` counts inputs
processed by physical preparation. Later stages are `Numerators`, `Selection`, `ZeroColor`,
`GroupingPreparation`, `GroupingSamples`, `GroupingComparison`, `Grouping`, and `Done`.
The grouping preparation, sample, and comparison phases count processed input diagrams;
`Done` reports retained output. Check the result's completion flag to distinguish cancellation.

Native workers publish snapshots; Python callbacks run on the thread invoking the operation,
including browser kernels. Marimo updates retain the calling cell context. Output flushing
alone does not yield to the browser for painting.

// docs-example: compile
```python
events = []
run = fk.Run(progress=lambda update: events.append(update) or True)
compact = request.generate(run=run)
assert compact.report.completed
```

Both `process.amplitude(..., classes=fk.Classes.auto())` and
`process.cross_section(..., classes=fk.Classes.auto())` default to safe flavour-class
refinement. Import `graph` from `symbolica.community` for the shared view types:
`fk.Classes is graph.Classes`, and `request.classes` uses this same type.
`graph.Classes.none()` assigns particles directly while retaining the same completed
particle/rule assignments and automorphism weights. In dressing callbacks,
`view.field(edge_index)` returns an owned `graph.Dressed` tagged `"class"` or `"particle"`;
`members` contains the original model `Particle` records for a class, while `particle` holds
a concrete record. `view.flow(edge_index)` can already be `graph.Flow.Source` or `Sink` when
`edges()[edge_index][2]` is still None. A pruning predicate must remain valid for every class
member and every later refinement. A retained `Dressed` outlives the callback, but its
expired view cannot be queried again.

Use `prune_skeletons` for hereditary pruning of untyped topologies and
`prune_dressings` for hereditary pruning of field assignments. Returning `False`
rejects that state and every refinement, without cancelling the run. `select`
receives only complete assignments and may apply an arbitrary predicate.

// docs-example: compile
```python
def no_self_edges(view):
    return all(view.self_loops(vertex) == 0 for vertex in range(view.n_vertices))

request = model.process(["e-", "e+"], ["mu-", "mu+"]).amplitude(
    loops=fk.Loops.total(1),
    topology=fk.TopologyRules(
        self_energies=fk.Veto.all(),
        tadpoles=fk.Veto.all(),
        snails=fk.Veto.massless(),
    ),
    prune_skeletons=no_self_edges,
)
result = request.generate().build()
```

This example illustrates the callback contract; for this particular condition, prefer the
built-in `TopologyRules(self_loops=False)` option to avoid Python callback overhead.

Hooks receive read-only `graph.SkeletonView`, `graph.DressingView` or `graph.DiagramView` values. They
provide degrees, open stubs, loops, edges, external legs, rules, field domains and
coupling-order bounds without building a graph. Assigned fields are `Particle`
objects; external labels are `(state, index)` pairs. Use `view.graph()` or
`view.complete_graph()` explicitly when a graph snapshot is needed. These snapshots
can be retained or modified independently; the typed view itself expires when its
callback returns. Numerator and UV-degree tests require the later built diagrams.

The search batches independent pending predicates and calls Python on the calling
thread after the native workers return. Every predicate must return a `bool`.
Exceptions propagate unchanged, and cancellation stops the remaining callbacks.

Register or replace a hook using `@request.prune_skeletons`,
`@request.prune_dressings` or `@request.select`; each decorator returns the original
callable. `@request.prune_dressings(thread_safe=True)` permits only that hook on the
existing predicate worker. Unmarked hooks and progress stay on the invoking thread;
`Run(thread_safe=True)` explicitly permits all hooks. WebAssembly runs them inline.
Registering the same hook again replaces both its callable and worker permission.

Set `request.check_heredity = True` to debug hereditary pruning. Either this request
preference or `Run(check_heredity=True)` enables the native check; neither can disable
an explicitly enabled check on the other. The request preference is not persisted.
`print(request)` shows physical defaults, active hooks, worker permissions and this
preference.

== Assemble and square amplitudes

`Amplitude` owns a coherent collection of finalized diagrams, aligning external tensor
ports by physical leg identity. It retains the source diagrams and individual weighted
terms. The operator includes each diagram's overall factor and numerator prefactor
exactly once; the diagnostic automorphism order is not a second symmetry factor.

// docs-example: compile
```python
amplitude = result.into_amplitude(fk.AmplitudeOptions(dimension=4))
operator = amplitude.expression()
adjoint = amplitude.conjugate().expression()
squared = amplitude.squared()
unpolarized = squared.sum_spins(average_initial=True).sum_colors(average_initial=True)
tensor = unpolarized.expression().simplify_algebra(gamma=True, color=True, epsilon=True)
tensor = tensor.contract(collect_chains=False, collect_traces=False)
```

These expressions are amputated and unintegrated. External wavefunctions remain open
spin/color ports until explicitly summed. `legs` exposes each physical particle,
incoming/outgoing state, momentum, and bare tensor index. `terms` exposes individual
`TensorExpression` operators. Conjugation uses Idenso's physical-leg Dirac adjoint and
color conjugation. Model-declared real parameters and physical momenta are real;
additional assumptions can be supplied with `AmplitudeOptions(real=[...])`. Other scalar symbols,
including complex couplings, retain their conjugation. Named model couplings are expanded
into their analytic parameter expressions before conjugation.

Without squared-order constraints, squaring includes every interference. With
`Orders().squared_exact(...)` or other squared bounds, an amplitude contribution needs one
common partner satisfying every bound jointly. The resolved class-pair relation is retained
through selection, zero removal, and grouping. `built.into_amplitude()` transfers it so that
squaring includes only permitted pairs. Constructing `Amplitude(list(built))` explicitly
starts an unrestricted amplitude and loses that request policy.

Ket and bra indices remain separate, as do their
loop integration variables. `sum_spins([0, 1], average_initial=True)` closes selected
legs; omitting the labels selects all remaining legs. An explicitly repeated sum raises
an error, preventing repeated averaging. Vector `references={leg: momentum}` and Dirac
`spin_vectors={leg: spin}` use the existing completeness API. Omitting a massless vector
reference chooses the covariant sum, which assumes a gauge-invariant amplitude.
Color sums close the actual dual representations of the open ports.

The Lorentz dimension defaults to four and can be set with `AmplitudeOptions(dimension=D)`; bispinor
spaces retain dimension four. Denominators follow the diagram's quadratic convention,
without widths or an imaginary prescription. Flux, identical-final-state factors, loop
measures, and phase-space integration are not supplied by `SquaredAmplitude`.
Already-sewn forward diagrams are rejected: their inverse-process right side must not
receive a second amplitude conjugation. GammaLoop's CFF/subtraction/sampling pipeline
continues to own those diagrams separately.

== Analyze interaction regions with Linnet

Finalized diagrams use GammaLoop's half-edge conventions. Amplitude external states are
dangling half-edges, incoming at a sink and outgoing at a source. Cross-section initial
states are sewn paired edges; their `is_external` metadata still identifies external
momentum carriers. Every item in `diagram.vertices` is an interaction. A missing edge
endpoint is `None`, and external names, indices, and states belong to the edge.

Graph analysis calls Rust Linnet directly and needs no Python `linnet` installation.
Filters receive FeynKit `DiagramVertex`, `DiagramEdge`, and `DiagramHalfEdge` objects.
Select native diagram IDs with `diagram.subgraph(nodes=..., edges=..., half_edges=...)`,
or use `diagram.filter(...)` to obtain a FeynKit `Subgraph`:

// docs-example: compile
```python
region = diagram.filter(edge=lambda edge: not edge.is_external)
numerator = region.numerator_expression()
components = region.connected_components()
boundary = region.boundary()
basis = region.momentum_basis()
routed = basis.route_expression(numerator)
```

`Subgraph` inherits `FeynmanDiagram` and retains shared ownership of its original diagram.
Its inherited operations use the selected region: denominators, momentum-basis enumeration,
parent-compatible and contracted routing, superficial divergence, CFF construction, UV
expansion, and tensor reduction. These methods have no `subgraph` argument.
Numerator `without` excludes edges and their incident vertex factors, using the shared
GammaLoop traversal. `region.original` accesses the full original diagram.
Partial tensor reduction requires an explicit `projector`; a boundary momentum is never
implicitly integrated. Overall factors, numerator prefactors, and projectors remain separate.

Nested `region.subgraph(...)` and `region.filter(...)` intersect with the current region.
Set operations `&`, `|`, `^`, `-`, and `~` preserve the original diagram owner.
Cycles are FeynKit `Subgraph` objects. `CutPartition` exposes `source_side`,
`target_side`, `boundary_left`, and `boundary_right` regions. `TraversalTree` returns
physics vertices in discovery order and supports parent, child, ancestor, and
fundamental-cycle queries. All result regions retain the immutable physics owner.
Structural cycles include sewn initial-state carriers, while the physical
`loop_count` excludes those carriers when constructing a momentum basis.

For optional interoperability, install the matching Python `linnet` extension.
`diagram.to_graph()` returns its canonical `Graph` with physics payloads, and
`diagram.subgraph(selection)` imports graph-bound Linnet selections. Exported half-edge
indices may differ from native IDs; the exported half-edge's `data` contains its native ID.
`region.to_graph()` returns the complete analysis graph; `region.linnet_selection` supplies
its corresponding canonical selection for direct Linnet algorithms. Importing a canonical
selection checks the graph owner and topology revision. A structural edit to the analysis
graph invalidates those canonical selections, while existing physics views retain their
immutable topology and remain usable.

Call `region.excise()` to obtain an independent `FeynmanDiagram`. Excision retains selected
vertices and their complete interaction slots, preserving selected paired edges and opening
other incident halves into external boundary legs. It remaps local edge, vertex, half-edge,
tensor, and momentum identities consistently. A proper extraction starts with unit global
factors, numerator prefactor, projector, and symmetry factor, and does not inherit whole-graph
physical cuts or add external wavefunctions. A full selection gives an independent exact
copy, including the original factors and routing. Serialization and operations requiring a
complete topology, such as integral-family construction, require explicit excision first.

`diagram.cuts` contains generated physical cuts with reusable `left.subgraph` and
`right.subgraph` physics views, coupling orders, side loop counts, crossing edges, and momentum signatures.
`diagram.topology_threshold_candidates` records the distinct process-independent partitions.
Linnet's `all_bonds` and `all_cuts` describe structural cuts, without assigning physical
final-state validity.

== Inspect the propagator denominator

Each edge exposes its model particle and its instantiated momentum and denominator:

// docs-example: compile
```python
edge = diagram.internal_edges[0]
mass = edge.particle.mass
width_parameter = edge.particle.width_parameter
momentum = edge.momentum_expression(dimension=4, in_lmb=True)
signature = edge.momentum_signature()
denominator = edge.denominator_expression(dimension=4, in_lmb=True)
```

The momentum is a rank-one Spenso vector; the signature exposes integer `loops`
and `external` coefficients. Both methods accept `lmb=basis` for an alternative
routing from the same diagram. `edge.propagator` returns the model's propagator
template for an internal physical line, or `None` when unavailable. Its expressions
contain model placeholders; the edge expression methods supply diagram factors.
External and dummy lines have no propagator denominator and raise `DiagramError`.
Use `power=2`, for example, for a squared edge denominator.

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

Use `in_lmb=True` on both expressions to replace edge momenta by the diagram's
stored loop and external momentum coordinates. Python already exposes the routing as
`LoopMomentumBasis`; obtain the stored basis with `diagram.loop_momentum_basis`, or
select another basis from `diagram.loop_momentum_bases()`.

// docs-example: compile
```python
numerator = diagram.numerator_expression(in_lmb=True)
denominator = diagram.denominator_expression(in_lmb=True)
basis = diagram.loop_momentum_bases()[0]
numerator = diagram.numerator_expression(lmb=basis)
denominator = diagram.denominator_expression(lmb=basis)
```

Supplying `lmb` enables routing and takes precedence over `in_lmb`. The basis must
belong to the same diagram instance; bases for a subgraph or contracted region are
accepted. Region selection, ignored numerator regions, denominator dimensions, and
propagator powers are applied before routing. Both methods return `TensorExpression`.

`basis.route_expression(expression)` also preserves `TensorExpression` inputs and
their ordered interfaces, including when routing cancels the tensor to zero.
Plain Symbolica expressions and scalar constants return a plain `Expression`.

Pass `loop_momenta` and `external_momenta` to name the independent basis vectors
using Spenso `TensorName.vector` objects. Loop names follow `basis.loop_edges`;
external names follow `basis.external_edges` with `basis.dependent_externals`
omitted. Each supplied list must cover its independent coordinates. Omitted lists
retain the canonical indexed `K` or `P` calls. The same naming applies to vectors
already expressed in the basis, preserving indexed and compact vector syntax.

// docs-example: compile
```python
from symbolica.community.tensor import TensorName

# For a one-loop, two-point graph:
K = TensorName.vector("example::K")
P = TensorName.vector("example::P")
routed = basis.route_expression(numerator, loop_momenta=[K], external_momenta=[P])
```

`kinematics.apply(expression)` follows the same convention: scalar-product
substitution preserves a tensor's ordered ports, even when an on-shell relation
makes it zero. Tensor methods can therefore be chained directly on the result.
Contraction and momentum-conservation rewrites remain explicit operations;
`apply` only substitutes the scalar products declared in the kinematics context.

For the local massive Taylor expansion, `diagram.uv_expansion(mUV, edge_powers=powers)`
and `region.uv_counterterm(mUV, edge_powers=powers)` accept the same signed power map.
Omitted internal edges have power one; entries outside the selected internal edges are ignored.
Powers change the integrand, while the selected region and its loop integration measure stay fixed.
When a prepared numerator contains an extra inverse propagator, move that factor into
`edge_powers` before UV expansion. For example, multiply a longitudinal massless-vector
numerator by $q_e^2$ and assign power two to that edge; this preserves the original integrand
and gives both denominator factors the shared auxiliary-mass expansion.

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

== Cancel reversed fermion loops before integration

The host notebook
#link("https://github.com/symbolica-dev/symbolica-community/blob/main/examples/hep/odd_photons.py")[hep/odd_photons.py]
uses typed particles and interaction rules to generate the one-, three- and
five-photon electron loops. Keep the one-point topology by relaxing the bridge
and tadpole constraints. The staged defaults impose neither; permit zero-flow edges explicitly:

// docs-example: compile
```python
electron, photon = (model.particle(name) for name in ("e-", "a"))
vertices = [
    vertex for vertex in model.vertex_rules
    if sorted(vertex.particles) == sorted([electron.name, electron.antiname, photon.name])
]
photon_count = 3
result = model.process([photon], [photon] * (photon_count - 1)).amplitude(
    loops=fk.Loops.total(1),
    restrict=fk.Restrict().only_vertices(vertices),
    topology=fk.TopologyRules(max_vertices=photon_count, zero_flow_edges=True),
).generate().build()
```

Contract the external Lorentz slots with independent symbolic probe vectors,
promote their Lorentz dimension to $D$, and use Idenso to trace the compact
graph-edge momenta. Apply `diagram.momentum_basis().route_expression` afterwards
to expand scalar products into the loop basis. Keep
`overall_factor_expression(evaluate=True)` and `numerator_prefactor_expression()`
in every numerator. Add the probes to the diagram family's external-vector
list, retaining its generated denominators, so an `IntegralMapping` transports
both loop-momentum and loop-probe scalar products.

With generic off-shell momenta, `find_mapping` identifies the unique reversed
partner of each triangle or pentagon. Verify its denominator permutation and
unit propagator powers, then add the mapped weighted numerators. Their exact
zero proves the full tensor identity without IBP or master integrals. The
one-photon numerator instead remains nonzero and odd in the loop momentum;
the shared vacuum tensor reducer integrates it to zero.

The notebook pairs the families before tracing and evaluates a selected pair.
The #source-link("crates/feynkit-py/tests/installed_feyncalc_odd_photons.py", label: "companion regression")
checks both triangles and all 24 pentagons, retaining symbolic dimension and
electron mass throughout. These vector-current examples do not fix general
gamma-five or Majorana conventions. The corrected shared Graphica automorphism
count reserves the endpoint-exchange factor for undirected self-loops; directed
fermion loops retain their orientation. The Symbolica Community host selects
that owner, while standalone adoption awaits its registry publication. The
#source-link("crates/feynkit-py/tests/installed_tadpole_normalization.py", label: "nonzero Higgs-fermion tadpole regression")
checks normalization, which a zero Furry amplitude cannot establish.

For a nonzero example, the
#link("https://github.com/symbolica-dev/symbolica-community/blob/main/examples/hep/tadpole_mass_insertions.py")[tadpole mass-insertion notebook]
differentiates the generated massive top-quark Higgs tadpole integrand with
respect to its propagator mass. Hold the Yukawa coupling and scale fixed, use
`IntegralFamily.rewrite_numerator` to identify the raised propagator powers,
and reduce powers one through three with native `IBPFamily.reduce_laporta`.
Keep the dimension symbolic until combining the $A_0$ pole with its
coefficients, then evaluate the finite coefficient using shared OneLOop.
The notebook also compares against numerical differentiation of the full
finite tadpole. It does not impose a tadpole renormalization condition.
HTML export and the live default plus twelve mass/scale states pass, including
zeros of the finite tadpole and its derivative; the scaled finite-difference
error stays below $2 times 10^(-9)$.

== Carry a generated triangle through its finite amplitude

The host notebook
#link("https://github.com/symbolica-dev/symbolica-community/blob/main/examples/hep/higgs_gluons.py")[hep/higgs_gluons.py]
starts from the model's top-quark interactions for $H -> g g$. Select typed
`VertexRule` objects by their particle content instead of relying on model-specific
vertex numbers:

// docs-example: compile
```python
top, higgs, gluon = (model.particle(name) for name in ("t", "H", "g"))
allowed_particles = [
    sorted([top.antiname, top.name, higgs.name]),
    sorted([top.antiname, top.name, gluon.name]),
]
vertices = [v for v in model.vertex_rules if sorted(v.particles) in allowed_particles]
result = model.process([higgs], [gluon, gluon]).amplitude(
    loops=fk.Loops.total(1),
    restrict=fk.Restrict().only_vertices(vertices),
    topology=fk.TopologyRules(
        max_vertices=3,
        bridges=graph.Bridges.at_most(0),
        self_energies=fk.Veto.all(),
        tadpoles=fk.Veto.all(),
        snails=fk.Veto.massless(),
    ),
).generate().build()
```

Keep both orientations and each diagram's complete weight. Promote Lorentz slots
to a dimension symbol before tracing, then give `TensorReducer` the loop vector
and both independent external directions. `IntegralFamily.rewrite_numerator`
produces scalar integral targets for native `IBPFamily.reduce_laporta`; OneLOop
supplies their Laurent coefficients. Reduction and master evaluation are distinct
steps. Substituting $D=4$ before the Laurent expansion would lose the finite
rational $2$ from the dimension-dependent bubble coefficient.

Use the generated tree Yukawa vertex to fix the phase. The
#link("https://feyncalc.github.io/FeynCalcExamples/QCD/OneLoop/H-GlGl")[FeynCalc example]
chooses `PreFactor -> -1`, which reverses the native $i cal(M)$ convention.
Shared physical gluon spin and color sums, two-body phase space and decay flux
then form the width, with one explicit $1/2!$ for identical gluons. Above the
quark-pair threshold, square the complex amplitude with its conjugate. The
#source-link("crates/feynkit-py/tests/installed_feyncalc_higgs_gluons.py", label: "companion regression")
checks this single-flavor leading-order calculation; additional flavors must be
summed coherently before forming a rate.

== Generate local counterterm insertions

The
#link("https://github.com/symbolica-dev/symbolica-community/blob/main/examples/hep/qed_renormalization.py")[QED renormalization notebook]
derives one-loop counterterm couplings by expanding the bare factors
$Z_psi$, $Z_psi Z_m$, $Z_A$, $Z_A/Z_xi$ and $Z_psi Z_e sqrt(Z_A)$ with
$Z_j=1+a_4 delta Z_j$. The local operator basis remains a model input.
Add these normalized rules with a bookkeeping coupling order `CT`, then use
the existing typed generation interface:

// docs-example: compile
```python
counterterms = ct_model.process([ct_electron], [ct_electron]).amplitude(
    loops=fk.Loops.total(0),
    restrict=fk.Restrict().only_vertices(ct_vertices),
    orders=fk.Orders().exact("QED", 2).exact("CT", 1),
    topology=fk.TopologyRules(max_vertices=1),
).generate().build()
```

The loop count describes topology; the `CT` order describes perturbative
bookkeeping. A two-point insertion adds no topological loop, so choose an
explicit vertex bound. Trace the generated numerators with the same
projectors and native factors as the loop diagrams, extract coefficients
linear in the unknown counterterms, and let Symbolica solve the matching
system. The auxiliary mass rule derives from $M(Z_(A m)^2-1)$, which
supplies $2 delta Z_(A m)$ at first order.

The notebook compares the massive full UV expansion with the massless
reference's direct massification prescription and verifies the MS/MSbar
measure conversion against shared OneLOop masters. No separate generator,
Dirac algebra, IBP solver or master evaluator is introduced for counterterms.
Automatic import of NLO-UFO counterterm metadata remains separate.

The
#link("https://github.com/symbolica-dev/symbolica-community/blob/main/examples/hep/qcd_renormalization.py")[QCD renormalization notebook]
extends this construction to eight quark, gluon and ghost insertions. Use
`Particle.color_sum` to close two-point color slots and preserve the generated
external ordering factors separately for fermions and ghosts. Copy the
original quark-gluon rule's full color/Lorentz coupling matrix before dressing
it with $Z_q Z_g sqrt(Z_A)-1$. The distinct ghost fields give an auxiliary-mass
coefficient $delta Z_(c m)$, whereas the identical gluons give
$2 delta Z_(A m)$. The generated coefficients determine the linear system
and verify both MS and MSbar subtraction. The notebook exposes the gauge,
flavor count, color count and subtraction scheme as live controls.
Its direct massive infrared-rearrangement path adds an auxiliary mass only to
massless denominators and keeps each longitudinal gluon's extra propagator
power explicit. After the external Taylor expansion and vacuum tensor
reduction, use Symbolica partial fractions to separate tadpoles at the physical
and auxiliary masses. Reconstructing the scalar integrand before native IBP
guards against dropping either contribution. The physical pole constants agree
between the two prescriptions, while the auxiliary gluon mass constant changes.
The separate
#link("https://github.com/symbolica-dev/symbolica-community/blob/main/examples/hep/qcd_massless_renormalization.py")[massless QCD notebook]
sets the quark mass to zero before this massification step. It uses Marimo's
`app.embed()` to access the shared calculation and exposes the seven remaining
counterterms in its own controls. This composition keeps the symbolic workflow
in one notebook while giving each physical case a focused presentation.

== Keep finite terms in an off-shell self-energy

The host notebook
#link("https://github.com/symbolica-dev/symbolica-community/blob/main/examples/hep/electron_self_energy.py")[hep/electron_self_energy.py]
generates the massive electron self-energy with a symbolic covariant-gauge
photon numerator. Two trace projectors extract the coefficients of
$slash(p)$ and $m$ without imposing an on-shell equation. Rewrite these
coefficients in the generated family's denominator coordinates before native
IBP reduction. The longitudinal photon contributes an additional inverse
photon denominator, so the target powers must retain that factor.

Keep each diagram's native weight. Its generated tree vertex establishes the
external ordering convention, which is then explicitly converted to the
amputated inverse-propagator insertion. The two residual masters are the
massive tadpole and massive-massless bubble; use their OneLOop Laurent
coefficients and keep $D$ symbolic until the series expansion. This retains
the finite rational terms from $D$ multiplying UV poles.

The notebook compares the result with independent parameter quadrature as
momentum, mass, scale and gauge change. It splits at the physical branch cut,
checks Landau gauge and the on-shell gauge-independent mass combination,
and evaluates the symbolic zero-momentum limit separately from the generic
$1/p^2$ projector. The
#source-link("crates/feynkit-py/tests/installed_feyncalc_electron_self_energy.py", label: "companion regression")
checks exact-dimensional identities, limits and sixty numerical points.
The finite answer is unrenormalized. The QED notebook above demonstrates
generated counterterm matching separately from scalar-integral reduction.

== Use finalized output

Discovery performs topology expansion, interaction assignment, pruning, and physical
preparation, including cuts, signs, and process weights. Count-only execution includes those
same decisions. Building then constructs numerator fragments, projectors, and routing;
zero-color removal is requested through `NumeratorOptions`. Selection and grouping are
explicit subsequent operations. Canonical diagram IDs and names are computed only when
an operation requests them. Compact traversal is deterministic, while grouping retains its
canonical comparison order. Requested names can change between revisions; select using
names from the current built result.

The returned diagrams retain vertex/edge numerator fragments, an external projector and
separate scalar prefactors. A generated aggregate numerator is built and cached on first use;
numerator grouping, expression access and serialization may request it. Reading `diagram.id`
or serializing the diagram similarly computes and caches its ID. JSON and annotated DOT
preserve an absent name separately from an empty name. Keep group multiplicities and prefactors
when forming a complete expression; a diagram's numerator alone is not the whole amplitude.

Named coefficients such as `UFO::GC_11` stay in generated numerators. For analytic work, expand
them using the model that owns the diagram:

// docs-example: compile
```python
analytic_numerator = model.expand_couplings(diagram.numerator_expression())
```

This returns a new expression and preserves the stored named coefficients. Serialize with
`to_json()` or `to_dot()` and restore with `FeynmanDiagram.from_json(model, text)` or
`from_dot(model, text)`. For a document containing multiple digraphs,
`FeynmanDiagram.from_dot_set(model, text)` returns a list in document order.
It accepts compact physics DOT and annotated exports; exports preserve each
diagram's momentum basis and physical cuts. An invalid graph fails the whole import.
The model fingerprint and structural validation protect the
interpretation of stored IDs. Model serialization keeps `UFO` implicit and
preserves all other Symbolica namespaces, so custom gauge symbols and tensor
functions retain their identities through `Model.to_json()` and
`Model.from_json()`. Distinct namespaces also remain distinct in the model
fingerprint.

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
