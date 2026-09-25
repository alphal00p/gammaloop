#import "../../shared.typ": callout, source-link

#let quickstart-python = [
= Using FeynKit from Python

FeynKit is a Symbolica community module. Use a combined Symbolica distribution built with
`feynkit-py`; a standalone FeynKit wheel is not produced by this repository. The
#link("guides/community-host/")[host integration guide] gives the registration and packaging
steps. Do not assume that an arbitrary published Symbolica wheel already includes this module.

// docs-example: syntax
```sh
python -c "import symbolica.community.feynkit as fk; print(fk.__name__)"
```

With that host installed, run the following from the GammaLoop repository root. It uses the same
small normalized model as the Rust quickstart and needs neither UFO import nor an external
integral-evaluation backend.

// docs-example: compile feynkit-community-quickstart
```python
import symbolica.community.feynkit as fk

model = fk.Model("crates/feynkit-model/tests/fixtures/scalars_2p_3p.json")
result = model.process(["scalar_0"], ["scalar_0", "scalar_0"]).generate_diagrams(loops=0, max_vertices=3)
assert result.diagrams

diagram = result.diagrams[0]
diagram.validate()
assert diagram.loop_count == 0
restored = fk.FeynmanDiagram.from_json(model, diagram.to_json())
assert restored.name == diagram.name
print(diagram.to_dot())

cff = diagram.build_cff()
print(cff.to_expression())
```

`model.process(incoming, outgoing)` binds external states and particle/vertex restrictions
to the model. Use `particle_veto`, `vertex_allow`, and `vertex_veto` at construction, or
`process.with_filters(...)` to create a restricted copy. Omitted filters are preserved;
`None` clears a filter, while `vertex_allow=[]` allows no interactions.
Use `process.generate_diagrams()` for a `GenerationResult`, `process.generate_amplitude()`
for an `Amplitude`, and `process.generate_cross_section()` for sewn diagrams and their cuts.
All three accept loop order (`loops=0` by default), coupling orders, symmetrization,
topology filters, and execution settings as keywords. These settings belong to each call. Particle names, PDG codes, and particles
obtained from that model can select external states.

In notebooks, displaying a `Process` draws the external states around a central blob.
`process.render()` exports the same schematic and accepts a `linnet.RenderConfig`.
Alternative final states are drawn separately.

Model members display their defining data and expressions, rather than only their names:

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
  Numerators are Spenso `TensorExpression` values, which extend Symbolica `Expression`; CFF
  expressions are ordinary `Expression` values. The FeynKit module
  must share `symbolica.core` with Spenso, Idenso, and other community modules; independently
  linking a second Symbolica extension breaks that ownership contract.
])

== Exact particle quantum numbers

`Particle.charge` is a Symbolica expression backed by Rust's `Rational`.
JSON stores exact integers and fractions as strings, such as `"-1"` and `"2/3"`.
Integer numeric literals are accepted too; fractional floating-point literals
are rejected because they do not specify the intended exact charge.

// docs-example: syntax
```python
from symbolica import E, S
from symbolica.community import feynkit as fk

charm = fk.Model.standard_model().particle("c")
assert charm.charge == E("2/3")
assert charm.y_charge == E("1/3")
assert charm.y_charge_right == E("4/3")
assert charm.weak_isospin == E("1/2")
assert charm.weak_isospin_right == E("0")
assert charm.mass_expression == S("UFO::MC")
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

`mass_expression` returns exact zero for the UFO `ZERO` parameter and a symbol
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
