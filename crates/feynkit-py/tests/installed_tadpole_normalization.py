"""Generated fermion and neutral scalar tadpoles agree with Wick contractions."""

from math import factorial, pi
from pathlib import Path

from symbolica import E, S, Symbol
from symbolica.community import hep
from symbolica.community.hep import oneloop
from symbolica.community.spenso import ColorSimplifySettings, TensorExpression

model = hep.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
top, higgs = (model.particle_by_pdg(pdg) for pdg in (6, 25))
yukawas = [
    vertex
    for vertex in model.vertex_rules
    if sorted(vertex.particles) == sorted([top.antiname, top.name, higgs.name])
]
cubics = [
    vertex for vertex in model.vertex_rules if vertex.particles == [higgs.name] * 3
]
quartics = [
    vertex for vertex in model.vertex_rules if vertex.particles == [higgs.name] * 4
]
assert len(yukawas) == len(cubics) == len(quartics) == 1

y, m, Nc, D, index = S("tadpole::y", "UFO::MT", "tadpole::Nc", "tadpole::D", "index_")
cof, mink, metric, K = S("spenso::cof", "spenso::mink", "spenso::g", "gammalooprs::K")
lam, vev, mh = S("UFO::lam", "UFO::vev", "UFO::MH")
left, right = S("left_", "right_")

# Independent Wick counts: psi and psibar can contract in only one way,
# with the closed-fermion-loop minus sign. Three identical scalars allow
# 3 external attachments at a 1/3! vertex; a quartic selfenergy permits
# 4*3 assignments of its two labelled external legs at a 1/4! vertex.
wick_weights = [E("-1"), E("3") / factorial(3), E("12") / factorial(4)]

# Establish the UFO vertex phases separately using generated trees.
for label, vertices, incoming, outgoing, coupling in [
    ("top", yukawas, [higgs], [top, top.antiparticle], y),
    ("scalar_cubic", cubics, [higgs], [higgs, higgs], 6 * vev * lam),
    ("scalar_quartic", quartics, [higgs, higgs], [higgs, higgs], 6 * lam),
]:
    trees = model.generate_diagrams(
        incoming,
        outgoing,
        max_vertices=1,
        vertex_allow=vertices,
        numerator_grouping=None,
        progress=None,
    ).diagrams
    assert len(trees) == 1, label
    tree = trees[0]
    kernel = (
        model.expand_couplings(tree.numerator_expression().to_expression())
        .replace(S("UFO::yt") * E("1/2").sqrt(), y)
        .replace(metric(left, right), E("1"))
    )
    assert (kernel + Symbol.I * coupling).expand() == E("0"), (label, kernel)
    assert tree.overall_factor_expression(evaluate=True) == E("1"), label
    assert tree.numerator_prefactor_expression() == E("1"), label

coefficients = {}
for label, vertices, outgoing, mass, expected_trace, wick_weight in [
    ("top", yukawas, [], m, 4 * m * y * Nc, wick_weights[0]),
    ("scalar_cubic", cubics, [], mh, 6 * vev * lam, wick_weights[1]),
    ("scalar_quartic", quartics, [higgs], mh, 6 * lam, wick_weights[2]),
]:
    generated = model.generate_diagrams(
        [higgs],
        outgoing,
        loops=1,
        max_vertices=1,
        # Momentum conservation makes the single external momentum zero.
        allow_zero_flow_edges=True,
        maximum_bridges=None,
        vertex_allow=vertices,
        tadpoles=None,
        zero_snails=None,
        self_energy=None,
        numerator_grouping=None,
        progress=None,
    )
    assert len(generated.diagrams) == 1, label
    diagram = generated.diagrams[0]
    assert diagram.overall_factor_expression(evaluate=True) == wick_weight, label
    assert diagram.numerator_prefactor_expression() == E("1"), label
    numerator = (
        model.expand_couplings(
            diagram.numerator_expression(in_lmb=True).to_expression()
        )
        .replace(S("UFO::yt") * E("1/2").sqrt(), y)
        .replace(cof(3, index), cof(Nc, index))
        .replace(mink(4, index), mink(D, index))
    )
    trace = (
        TensorExpression(numerator.expand())
        .simplify_gamma()
        .expand()
        .simplify_color(ColorSimplifySettings(substitute_cof_dimension_invariants=True))
        .simplify_metrics()
        .to_dots()
        .to_expression()
    )
    assert (trace - expected_trace).expand() == E("0"), (label, trace)
    if label == "top":
        # The Lorentz dimension is D while the Dirac identity has trace four.
        assert (trace / (m * y * Nc)).cancel() == E("4")
    kin = hep.Kinematics(D, momenta=[K(0)])
    family = diagram.integral_family(kinematics=kin)
    if label != "scalar_quartic":
        # The one-point families need no external scalar products or extra ISPs.
        assert family.rewrite_numerator(trace, [S("tadpole::d0")]) == trace
    assert len(family.denominators) == 1, label
    # The generated one-line family is precisely the ordinary massive A0.
    assert (
        family.denominators[0] - (kin.scalar_product(K(0), K(0)) - mass**2)
    ).expand() == E("0"), label
    coefficient = (
        trace
        * diagram.overall_factor_expression(evaluate=True)
        * diagram.numerator_prefactor_expression()
    ).expand()
    assert (coefficient / (expected_trace * wick_weight)).cancel() == E("1"), label
    coefficients[label] = coefficient

# Coefficients multiply i*A0/(16*pi^2). OneLOop uses
# A0 = m^2[1/eps + 1 - log(m^2/mu^2)] + O(eps).
# Check the shared master at mu^2=m^2=4 independently of graph weights.
master = oneloop.master_coefficients(S("oneloopmaster::A0")(E("4"), E("4")))
assert abs(complex(master[2].evaluate({}))) < 1.0e-14
assert abs(complex(master[1].evaluate({})) - 4) < 1.0e-14
assert abs(complex(master[0].evaluate({})) - 4) < 1.0e-14

# Nc=3, y=1, m=2 gives coefficient -24 and A0 pole/finite part four:
# both the pole and finite amplitude coefficient are -6i/pi^2 at mu=m.
top_coefficient = complex(coefficients["top"].evaluate({Nc: 3.0, y: 1.0, m: 2.0}))
for term in (master[1], master[0]):
    amplitude = 1j * top_coefficient * complex(term.evaluate({})) / (16 * pi**2)
    assert abs(amplitude + 6j / pi**2) < 1.0e-14

print(
    "Generated tadpole normalization: fermion Wick factor, scalar controls, and A0 passed"
)
