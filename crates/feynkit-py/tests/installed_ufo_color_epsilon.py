"""Generate conjugate three-triplet scalar vertices through the shared UFO path.

Three distinct complex scalars avoid an identically vanishing antisymmetric
interaction of identical bosons. The SM fixture supplies scalar parameters and
Lorentz structures only; these color interactions define a toy model.
"""

import json
from pathlib import Path

from symbolica import E, S
from symbolica.community import hep as fk
from symbolica.community.spenso import TensorExpression

source = json.loads(
    (Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json").read_text()
)
source["name"] = "triplet_scalar_epsilon_example"
scalar = next(p for p in source["particles"] if p["pdg_code"] == 251)
for flavor in range(1, 4):
    for sign in (1, -1):
        name = f"X{flavor}" + ("~" if sign < 0 else "")
        antiname = f"X{flavor}" + ("~" if sign > 0 else "")
        source["particles"].append(
            dict(
                scalar,
                name=name,
                antiname=antiname,
                texname=name,
                antitexname=antiname,
                pdg_code=sign * (9000000 + flavor),
                color=sign * 3,
                charge=0.0,
            )
        )
vertex = next(v for v in source["vertex_rules"] if v["name"] == "V_9")
source["vertex_rules"] = [
    dict(
        vertex,
        name="epsilon" + ("_bar" if conjugate else ""),
        particles=[f"X{flavor}" + ("~" if conjugate else "") for flavor in range(1, 4)],
        color_structures=[("EpsilonBar" if conjugate else "Epsilon") + "(1,2,3)"],
    )
    for conjugate in (False, True)
]
model = fk.Model.from_json(json.dumps(source))
for sign in (1, -1):
    generated = fk.Generator(model).generate(
        fk.Process.amplitude(
            [sign * 9000001], [-sign * 9000002, -sign * 9000003]
        ).with_loop_count(0, 0),
        max_vertices=1,
        maximum_bridges=None,
        numerator_grouping=None,
        progress=None,
    )
    assert len(generated.diagrams) == 1
    color = TensorExpression(
        generated.diagrams[0]
        .numerator_expression()
        .to_expression()
        .replace(S("UFO::GC_69"), E("1"))
    )
    assert len(color.interface) == 3
    conjugate = color.spenso_conjugate()
    assert conjugate.spenso_conjugate().to_expression() == color.to_expression()
    norm = TensorExpression(
        color.to_expression() * conjugate.to_expression()
    ).simplify_epsilon()
    assert norm.is_scalar
    assert norm.to_expression() == E("6")
    # One incoming triplet is averaged; the two distinct outgoing species are summed.
    average = TensorExpression(
        model.particle_by_pdg(sign * 9000001).color_sum(S("i"), S("i"), average=True)
    ).simplify_metrics()
    assert average.to_expression() == E("1")
    print(
        f"Generated {'Epsilon' if sign > 0 else 'EpsilonBar'}: dual conjugate, norm 6 passed"
    )
