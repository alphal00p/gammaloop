"""Check SM ghost vertices against the original UFO in all three crossings.

Upstream records: https://github.com/mg5amcnlo/mg5amcnlo/tree/3.x/models/sm
The particle order and couplings below come from vertices.py/couplings.py.
lorentz.py defines UUV1 = P(3,2) + P(3,3), which equals -p1 with all
three momenta incoming. No external MadGraph installation is required.
"""

import copy
import json
from pathlib import Path

from symbolica import E, S
from symbolica.community import hep
from symbolica.community.spenso import TensorExpression

model = hep.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
specification = json.loads(model.to_json())
vertices = {
    "V_12": (["ghA", "ghWm~", "W-"], "GC_3"),
    "V_13": (["ghA", "ghWp~", "W+"], "GC_4"),
    "V_15": (["ghWm", "ghA~", "W+"], "GC_3"),
    "V_18": (["ghWm", "ghWm~", "a"], "GC_4"),
    "V_19": (["ghWm", "ghWm~", "Z"], "GC_53"),
    "V_21": (["ghWm", "ghZ~", "W+"], "GC_52"),
    "V_23": (["ghWp", "ghA~", "W-"], "GC_4"),
    "V_26": (["ghWp", "ghWp~", "a"], "GC_3"),
    "V_27": (["ghWp", "ghWp~", "Z"], "GC_52"),
    "V_29": (["ghWp", "ghZ~", "W-"], "GC_53"),
    "V_31": (["ghZ", "ghWm~", "W-"], "GC_52"),
    "V_33": (["ghZ", "ghWp~", "W+"], "GC_53"),
    "V_35": (["ghG", "ghG~", "g"], "GC_10"),
}
couplings = {
    "GC_3": E("-1i*UFO::ee"),
    "GC_4": E("1i*UFO::ee"),
    "GC_52": E("-1i*UFO::cw*UFO::ee/UFO::sw"),
    "GC_53": E("1i*UFO::cw*UFO::ee/UFO::sw"),
    "GC_10": E("-UFO::G"),
}
rules = {
    vertex["name"]: vertex
    for vertex in specification["vertex_rules"]
    if vertex["lorentz_structures"] == ["UUV1"]
}
assert set(rules) == set(vertices)
p2 = "UFO::P(UFO::idx(1,3),UFO::idx(1,2))"
p3 = "UFO::P(UFO::idx(1,3),UFO::idx(1,3))"
uuv = next(
    item for item in specification["lorentz_structures"] if item["name"] == "UUV1"
)
assert E(uuv["structure"]) == E(p2 + "+" + p3)
P, mink, hedge, wave, index, mu = S(
    "gammalooprs::P",
    "spenso::mink",
    "gammalooprs::hedge",
    "ghost_crossing::wave_",
    "ghost_crossing::index_",
    "ghost_crossing::mu",
)
zero, one = E("0"), E("1")
for name, (particles, coupling_name) in vertices.items():
    rule = rules[name]
    assert rule["particles"] == particles
    assert rule["couplings"] == [[coupling_name]]
    if name != "V_35":
        assert rule["color_structures"] == ["1"]
    coupling = couplings[coupling_name]
    assert (
        model.expand_couplings(S("UFO::" + coupling_name)) - coupling
    ).together() == zero
    ghost, antighost, vector = [model.particle(particle) for particle in particles]
    # P0 = P1 + P2 for the physical process. The expected tensor is -p1
    # of the ordered UFO ghost leg, independent of the stored Lorentz expression.
    # slots records the external indices of (ghost, antighost, vector).
    for incoming, outgoing, slots, momentum in (
        (
            [ghost],
            [model.particle(vector.antiname), model.particle(antighost.antiname)],
            (0, 2, 1),
            -P(0, mink(4, mu)),
        ),
        (
            [vector],
            [model.particle(ghost.antiname), model.particle(antighost.antiname)],
            (1, 2, 0),
            P(1, mink(4, mu)),
        ),
        (
            [antighost],
            [model.particle(vector.antiname), model.particle(ghost.antiname)],
            (2, 0, 1),
            P(0, mink(4, mu)) - P(1, mink(4, mu)),
        ),
    ):
        diagrams = (
            hep.Process(model, incoming, outgoing)
            .generate_diagrams(
                loops=0,
                max_vertices=1,
                maximum_bridges=0,
                vertex_allow=[name],
                numerator_grouping=None,
                progress=None,
            )
            .diagrams
        )
        assert len(diagrams) == 1, (name, slots)
        diagram = diagrams[0]
        edges = {edge.external_index: edge for edge in diagram.external_edges}
        assert edges[0].external_state == "incoming"
        assert edges[0].source is None and edges[0].target == 0
        assert all(edges[i].external_state == "outgoing" for i in (1, 2))
        assert all(edges[i].source == 0 and edges[i].target is None for i in (1, 2))
        port = dict(
            next(
                diagram.projector_expression().match(
                    wave(edges[slots[2]].id, mink(4, index)), max_level=0
                )
            )
        )[index]
        actual = model.expand_couplings(
            diagram.numerator_expression(in_lmb=True).to_expression()
        ).replace(mink(4, port), mink(4, mu))
        color = one
        if name == "V_35":
            # Keep the full upstream f(a_ghost, a_antighost, a_gluon), including
            # its crossing sign. Only obtain index labels from graph metadata.
            half_edges = json.loads(diagram.to_json())["half_edge_order"]
            assert len(half_edges) == 3
            indices = {
                half["edge"]: hedge(number, 1) for number, half in enumerate(half_edges)
            }
            color = TensorExpression.f(8)(
                *(indices[edges[external].id] for external in slots)
            ).to_expression()
        assert diagram.overall_factor_expression(evaluate=True) == one
        assert diagram.numerator_prefactor_expression() == one
        residual = (actual - coupling * color * momentum).expand()
        assert residual == zero, (name, slots, residual)

# Generation must also be linear before applying momentum conservation. Signed
# remapping of Q inside a sum once collapsed -(-Q) and visited the new Q again;
# a coefficient of two prevented that collapse and exposed the inconsistency.
linearity = {}
for label, tensor in (
    ("P2", p2),
    ("P3", p3),
    ("P2+P3", p2 + "+" + p3),
    ("P2-P3", p2 + "-" + p3),
    ("P2+2P3", p2 + "+2*" + p3),
):
    spec = copy.deepcopy(specification)
    next(item for item in spec["lorentz_structures"] if item["name"] == "UUV1")[
        "structure"
    ] = tensor
    probe_model = hep.Model.from_json(json.dumps(spec))
    for incoming, outgoing, name in (
        (["a"], ["ghWm~", "ghWm"], "V_18"),
        (["Z"], ["ghWm~", "ghWm"], "V_19"),
        (["g"], ["ghG~", "ghG"], "V_35"),
    ):
        diagrams = (
            hep.Process(probe_model, incoming, outgoing)
            .generate_diagrams(
                loops=0,
                max_vertices=1,
                maximum_bridges=0,
                vertex_allow=[name],
                numerator_grouping=None,
                progress=None,
            )
            .diagrams
        )
        assert len(diagrams) == 1, (name, label)
        linearity[name, label] = diagrams[0].numerator_expression().to_expression()
for name in ("V_18", "V_19", "V_35"):
    for label, coefficient in (("P2+P3", 1), ("P2-P3", -1), ("P2+2P3", 2)):
        residual = (
            linearity[name, label]
            - linearity[name, "P2"]
            - coefficient * linearity[name, "P3"]
        ).expand()
        assert residual == zero, (name, label, residual)

print("Passed 36 EW and 3 QCD ghost crossings, plus 15 Lorentz-linearity probes.")
