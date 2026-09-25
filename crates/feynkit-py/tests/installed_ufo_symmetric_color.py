"""Generate a toy adjoint-scalar vertex with the standard UFO d^{abc} tensor.

The existing SM fixture supplies scalar mass/coupling records only; the modified
model is a toy color interaction, not the Standard Model Higgs interaction.
Reference normalization: https://feyncalc.github.io/FeynCalcBook/SUND.html
"""

import json
from pathlib import Path

from symbolica import E, S
from symbolica.community import hep as fk
from symbolica.community.spenso import TensorExpression

source = json.loads(
    (Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json").read_text()
)
source["name"] = "adjoint_scalar_example"
for particle in source["particles"]:
    if particle["pdg_code"] == 25:
        particle["color"] = 8
source["vertex_rules"] = [v for v in source["vertex_rules"] if v["name"] == "V_9"]
source["vertex_rules"][0]["color_structures"] = ["d(1,2,3)"]
model = fk.Model.from_json(json.dumps(source))
generated = model.process([25, 25], [25]).generate_diagrams(
    loops=0,
    max_vertices=1,
    maximum_bridges=None,
    numerator_grouping=None,
    progress=None,
)
assert len(generated.diagrams) == 1
diagram = generated.diagrams[0]
color = TensorExpression(
    diagram.numerator_expression().to_expression().replace(S("UFO::GC_69"), E("1"))
)
assert len(color.structure.slots) == 3
assert color.spenso_conjugate().to_expression() == color.to_expression()
# The Symbolica product contracts all three matching explicit adjoint indices.
norm = (
    TensorExpression(color.to_expression() * color.spenso_conjugate().to_expression())
    .simplify_color()
    .to_cof_dimension_invariants()
)
assert norm.is_scalar
assert norm.to_expression() == E("40/3")
# Initial-state color averaging is independently supplied by the particle API.
particle = model.particle_by_pdg(25)
assert TensorExpression(
    particle.color_sum(S("x"), S("x"), average=True)
).simplify_metrics().to_expression() == E("1")
print("Generated UFO d tensor: real adjoint interface and SU(3) norm 40/3 passed")
