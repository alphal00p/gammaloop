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
Filter objects specify which subgraphs to reject. Omitting a filter leaves it disabled.

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
result = fk.Generator(model).generate(process, max_vertices=6)
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

== Inspect the propagator denominator

// docs-example: compile
```python
denominator = diagram.denominator_expression()
integrand = diagram.numerator_expression() / denominator
```

The denominator is a scalar Spenso `TensorExpression`: the product of
$q_e^2 - m_e^2$ over internal edges, with four-dimensional Minkowski scalar products.
Its momentum labels match the numerator and its masses remain symbolic model parameters;
the UFO `ZERO` mass becomes zero. External legs, widths, and an imaginary prescription
are excluded. This follows FeynKit's quadratic-propagator convention; custom UFO
denominator formulas are not instantiated. A diagram without internal edges returns one.
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
