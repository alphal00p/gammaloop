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
options = fk.GenerationOptions(max_vertices=3)
result = model.generate_diagrams(
    ["scalar_0"], ["scalar_0", "scalar_0"], loops=0, options=options
)
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

`Model`, `Generator`, and `FeynmanDiagram` own their respective operations. Use
`model.generate_diagrams(...)` for a short workflow, or `Generator(model).generate(process,
options)` when reusing a configured `Process`. Particle names, PDG codes, and particles obtained
from that model can select external states.

#callout("Keep the same Symbolica kernel", [
  Numerators and CFF expressions are native Symbolica `Expression` values. The FeynKit module
  must share `symbolica.core` with Spenso, Idenso, and other community modules; independently
  linking a second Symbolica extension breaks that ownership contract.
])

Continue with #link("tutorial/")[model inspection and generation], the
#link("guides/tensor-reduction/")[tensor selector example], or
#link("guides/notebooks/")[notebook rendering]. The complete signature and docstring inventory
is in the #link("reference/python/feynkit-community/")[Python API reference].
]
