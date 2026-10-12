#import "../../shared.typ": callout, source-link

#let quickstart-python = [
= Using FeynKit from Python

FeynKit is a Symbolica community module. Use a combined Symbolica distribution built with
`feynkit-py`; a standalone FeynKit wheel is not produced by this repository. The
#link("guides/community-host/")[host integration guide] gives the registration and packaging
steps. Do not assume that an arbitrary published Symbolica wheel already includes this module.

// docs-example: syntax
```sh
python -c "import symbolica.community.hepkit as fk; print(fk.__name__)"
```

With that host installed, run the following from the GammaLoop repository root. It uses the same
small normalized model as the Rust quickstart and needs neither UFO import nor an external
integral-evaluation backend.

// docs-example: compile feynkit-community-quickstart
```python
import symbolica.community.hepkit as fk

model = fk.Model("crates/feynkit-model/tests/fixtures/scalars_2p_3p.json")
request = model.process(["scalar_0"], ["scalar_0", "scalar_0"]).amplitude(
    loops=fk.Loops.total(0),
    topology=fk.TopologyRules(max_vertices=3),
)
counts = request.count()
assert counts.complete and counts.diagrams > 0
compact = request.generate()
assert len(compact) == counts.diagrams
result = compact.build()

diagram = result[0]
diagram.validate()
assert diagram.loop_count == 0
restored = fk.FeynmanDiagram.from_json(model, diagram.to_json())
assert restored.name == diagram.name
print(diagram.to_dot())

cff = diagram.build_cff()
print(cff.to_expression())
```

`model.process(incoming, outgoing)` validates external states against the model. Particle
names, PDG codes, and particles obtained from that model can select those states.
Restrictions belong to the request: pass `restrict=fk.Restrict().only_vertices([...])`
or `.exclude_particles([...])` to `process.amplitude(...)` or `process.cross_section(...)`.
Allow clauses intersect and exclusions unite; `.only_vertices([])` permits no interactions.

Every physics request requires a finite loop interval: `fk.Loops.total(0)` for a tree
amplitude, or `fk.Loops.forward(1)` for a forward graph with a two-particle tree cut.
`request.explain()` describes planning constraints; `request.count()` returns exact discovery
totals without constructing numerators. `request.generate()` returns a compact `DiagramSet`,
whose indexing and iteration expose prepared assignments. `compact.build(fk.NumeratorOptions(...))`
constructs a `BuiltDiagramSet`; use its `select`, `group`, and `into_amplitude` methods for
those later operations. Threading, cancellation, and progress belong to `fk.Run(...)` passed
to the operation. None of these calls changes the process or request.

Shared defaults apply no particle, coupling-order, bridge, or cyclic-subgraph restriction.
Self-loops are allowed and zero-flow edges are rejected. Names are not computed by default:
pass `fk.NumeratorOptions(naming=fk.Naming.named("D"))` to `build` when names are needed.
The example above leaves `diagram.name` as `None`, including after serialization.

In notebooks, displaying a `Process` draws the external states around a central blob.
`process.render()` returns a displayable `DiagramRender` and accepts typed settings,
for example `process.render(config=graph.RenderSettings(drawing=graph.DrawOptions(node_radius=5)))`
after importing `graph` from `symbolica.community`. Use `help(graph.RenderSettings)`
to discover the options.
Use the result's `to_svg()` to export the schematic as text. Alternative final
states are drawn separately.

Model members display their defining data and expressions, rather than only their names:

Particle masses and physical electric charges are symbolic expressions. The
particle charge includes its sign and is inferred from photon interaction vertices:

```python
from symbolica.community import hepkit as hep

qed_model = hep.Model.qed()
electron = qed_model.particle("e-")
mass, charge = electron.mass, electron.electric_charge
e = -electron.electric_charge
assert charge == -e
assert electron.charge == -1  # dimensionless charge in units of e
```

The charge accessor uses the normalization of minimal photon–fermion or
photon–scalar vertices, without assuming parameter names. A neutral particle
has zero electric charge even in a model without photon interactions.

```python
model.vertex_rules[0]  # ordered legs, color/Lorentz structures, and coupling terms
model.particle("scalar_0")  # quantum numbers, mass, width, and antiparticle
print(model.vertex_rules[0])  # a complete plain-text summary
```

The same rich and plain-text printers cover `model.parameters`, `model.couplings`,
`model.lorentz_structures`, `model.propagators`, `model.functions`, and `model.form_factors`.
For a collection in marimo, use `mo.vstack(model.vertex_rules)`; in a terminal, iterate and
`print` the members or use `pprint` on the list.

#callout("Keep the same Symbolica kernel", [
  Numerators are Spenso `TensorExpression` values, subclasses of Symbolica's `Expression`
  with an ordered tensor interface; CFF expressions are ordinary `Expression` values. Algebraic
  methods such as `expand_num()` and `collect_factors()` preserve the tensor type.
  Inherited Symbolica methods remain available directly. Use `to_expression()` to explicitly
  discard the tensor interface and select ordinary expression behavior. The FeynKit module
  must share `symbolica.core` with Spenso, Idenso, and other community modules; independently
  linking a second Symbolica extension breaks that ownership contract.
])

== Exact particle quantum numbers

`Particle.charge` is a Symbolica expression backed by Rust's `Rational`.
JSON stores exact integers and fractions as strings, such as `"-1"` and `"2/3"`.
Integer numeric literals are accepted too; fractional floating-point literals
are rejected because they do not specify the intended exact charge.

// docs-example: compile
```python
from symbolica import E, S
from symbolica.community import hepkit as fk

charm = fk.Model.standard_model().particle("c")
assert charm.charge == E("2/3")
assert charm.y_charge == E("1/3")
assert charm.y_charge_right == E("4/3")
assert charm.weak_isospin == E("1/2")
assert charm.weak_isospin_right == E("0")
assert charm.mass == S("UFO::MC")
```

Hypercharge follows $Q = T_3 + Y/2$, as in the
#link("https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf")[PDG electroweak review].
For fermions, `y_charge` and `weak_isospin` refer to the left-handed component;
`y_charge_right` and `weak_isospin_right` refer to the right-handed component.
Charge conjugation exchanges chirality and reverses the charges: the positron
has left/right hypercharges $2, 1$ and weak-isospin components $0, 1/2$.
`None` means missing metadata, an absent chiral field, or a state without a
definite hypercharge. In particular, the SM has no right-handed neutrino field,
and the neutral real Higgs and Goldstone fields are not hypercharge eigenstates.
Weak isospin is derived from the model's charges, without a PDG-code lookup.

The UFO adapter preserves `Y`, `YRight`, and the optional `charge_exact`
attribute. Use fraction strings in those attributes to avoid rounding.
Floating-point UFO attributes retain their decimal values as rationals; the
adapter does not guess a small-denominator fraction from a rounded number.

`mass` returns exact zero for the UFO `ZERO` parameter and a symbol
for every other mass parameter, including parameters whose current value is
zero. The shipped SM values match `restrict_default.json`: electron, muon,
and charm masses and their Yukawa inputs are zero, and the CKM matrix is
identity. These are parameter-card choices; the complete particle and vertex
content is retained. The disk JSON, fixture, and embedded model share the
same Feynman-gauge propagators, Goldstone flags, and Lorentz structures.

Try the #link("guides/showcases/first-diagram/")[interactive first-diagram showcase] or choose
a component from the #link("guides/showcases/")[notebook gallery].

Continue with #link("tutorial/")[model inspection and generation], the
#link("guides/tensor-reduction/")[tensor selector example], or
#link("guides/notebooks/")[notebook rendering]. The complete signature and docstring inventory
is in the #link("reference/python/feynkit-community/")[Python API reference].
]
