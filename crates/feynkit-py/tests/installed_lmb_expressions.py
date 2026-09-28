"""Exercise optional expression routing in an installed FeynKit community host.

Pass the community module name as the first argument (``hep`` for that host).
"""

import importlib
import json
import sys
from pathlib import Path

from symbolica import E, Expression, S
from symbolica.community.spenso import Representation, TensorExpression, TensorName, dot

fk = importlib.import_module(
    f"symbolica.community.{sys.argv[1] if len(sys.argv) > 1 else 'feynkit'}"
)
model = fk.Model(Path(__file__).parent / "fixtures/scalars_2p_3p.json")
diagrams = (
    model.process(
        ["scalar_0"], ["scalar_0", "scalar_0"], vertex_allow=["V_3_SCALAR_000"]
    )
    .generate_diagrams(loops=1, max_vertices=3, allow_self_loops=True)
    .diagrams
)
diagram = next(
    candidate
    for candidate in diagrams
    if len(
        {
            edge.source if edge.source is not None else edge.target
            for edge in candidate.external_edges
        }
    )
    == 3
)
loop_edge = next(
    edge
    for edge in diagram.internal_edges
    if edge.id == diagram.loop_momentum_basis.loop_edges[0]
)

# Give one vertex an indexed momentum on a real edge and keep every local
# numerator consistent with the aggregate serialized numerator.
payload = json.loads(diagram.to_json())
payload["numerator"] = f"gammalooprs::Q({loop_edge.id},spenso::mink(4,lmb_test::mu))"
for vertex in payload["vertices"]:
    vertex["numerator"] = "1"
for _, edge in payload["edges"]:
    edge["numerator"] = "1"
payload["vertices"][loop_edge.source]["numerator"] = payload["numerator"]
diagram = fk.FeynmanDiagram.from_json(model, json.dumps(payload))
basis = diagram.loop_momentum_basis
assert isinstance(basis, fk.LoopMomentumBasis)
alternate = next(
    candidate
    for candidate in diagram.loop_momentum_bases()
    if candidate.loop_edges != basis.loop_edges
)
raw_numerator = diagram.numerator_expression()
assert raw_numerator.rank == 1
assert raw_numerator == diagram.numerator_expression(in_lmb=False, lmb=None)
assert basis.route_expression(raw_numerator) != raw_numerator
assert alternate.route_expression(raw_numerator) != basis.route_expression(
    raw_numerator
)

# Routing preserves the caller's expression type, including open ports and
# ordered interfaces which cannot be recovered from a zero atom alone.
lorentz = Representation.mink(4)
mu, nu = lorentz("routing_type_mu"), lorentz("routing_type_nu")
Q = TensorName.vector("gammalooprs::Q", tags=S("gammalooprs::Q").get_tags())
for selected_basis in (basis, alternate):
    indexed = Q(loop_edge.id, mu)
    # Expose numerical signs before routing, which only substitutes momenta.
    cancelling = indexed - selected_basis.route_expression(indexed.to_expression())
    cancelling = TensorExpression(
        cancelling.to_expression().expand_num(), structure=cancelling.structure
    )
    for tensor in (
        indexed,
        Q(loop_edge.id, lorentz),
        Q(loop_edge.id, nu) * indexed,
        TensorExpression(3),
        0 * indexed,
        cancelling,
    ):
        routed = selected_basis.route_expression(tensor)
        assert isinstance(routed, TensorExpression)
        assert routed.structure.slots == tensor.structure.slots
        plain = selected_basis.route_expression(tensor.to_expression())
        assert isinstance(plain, Expression) and not isinstance(plain, TensorExpression)
        assert routed.to_expression() == plain
    cancelled = selected_basis.route_expression(cancelling)
    assert not cancelled
    assert cancelled.rank == 1
    for scalar in (E("3"), 3, 2.5):
        routed = selected_basis.route_expression(scalar)
        assert isinstance(routed, Expression) and not isinstance(
            routed, TensorExpression
        )
        assert routed == scalar

for expression in (diagram.numerator_expression, diagram.denominator_expression):
    raw = expression()
    assert raw == expression(in_lmb=False, lmb=None)
    for options, selected_basis in (
        ({"in_lmb": True}, basis),
        ({"lmb": basis}, basis),
        ({"lmb": alternate}, alternate),
        ({"in_lmb": False, "lmb": alternate}, alternate),
        ({"in_lmb": True, "lmb": alternate}, alternate),
    ):
        routed = expression(**options)
        assert isinstance(routed, TensorExpression)
        assert routed.structure.slots == raw.structure.slots
        assert routed == selected_basis.route_expression(raw)
        assert routed != raw

graph = diagram.to_linnet()
full = diagram.subgraph(graph.full_subgraph())
empty = diagram.subgraph()
region = diagram.filter(edge=lambda edge: edge.data.id == loop_edge.id)
region_basis = region.momentum_basis()
for name in ("numerator_expression", "denominator_expression"):
    assert getattr(full, name)(in_lmb=True) == getattr(diagram, name)(in_lmb=True)
    for selection, selected_basis in ((region, region_basis), (empty, basis)):
        expression = getattr(selection, name)
        raw = expression()
        assert expression(in_lmb=True) == basis.route_expression(raw)
        routed = expression(lmb=selected_basis)
        assert isinstance(routed, TensorExpression)
        assert routed == selected_basis.route_expression(raw)
    assert getattr(empty, name)(in_lmb=True) == TensorExpression(1)
    assert getattr(empty, name)(lmb=alternate) == TensorExpression(1)

without_region = diagram.numerator_expression(without=region)
assert without_region == TensorExpression(1)
assert diagram.numerator_expression(without=region, in_lmb=True) == (
    basis.route_expression(without_region)
)
assert diagram.numerator_expression(without=region, lmb=alternate) == (
    alternate.route_expression(without_region)
)

powers = {
    edge.id: power
    for edge, power in zip(diagram.internal_edges, (2, -1, 0), strict=True)
}
for dimension in (4, S("lmb_test::D")):
    options = {"edge_powers": powers, "dimension": dimension}
    raw = diagram.denominator_expression(**options)
    for selected_basis in (basis, alternate):
        routed = diagram.denominator_expression(**options, lmb=selected_basis)
        assert isinstance(routed, TensorExpression)
        assert routed.rank == 0
        assert routed == selected_basis.route_expression(raw)
    assert diagram.denominator_expression(**options, in_lmb=True) == (
        basis.route_expression(raw)
    )
    partial = region.denominator_expression(**options)
    assert region.denominator_expression(
        **options, lmb=region_basis
    ) == region_basis.route_expression(partial)

restored = fk.FeynmanDiagram.from_json(model, diagram.to_json())
foreign_basis = restored.loop_momentum_basis
for expression in (diagram.numerator_expression, diagram.denominator_expression):
    for in_lmb in (False, True):
        try:
            expression(in_lmb=in_lmb, lmb=foreign_basis)
        except fk.DiagramError as error:
            assert "momentum basis belongs to a different diagram" in str(error)
        else:
            raise AssertionError("a supplied basis must belong to the diagram instance")


# Edge access retains the original diagram's identity and uses its shared
# denominator builder, even when the edge was obtained through a subgraph.
for dimension in (4, S("lmb_test::D")):
    for selected_basis in (None, basis, alternate):
        product = TensorExpression(1)
        for edge in diagram.internal_edges:
            options = {"dimension": dimension, "lmb": selected_basis}
            denominator = edge.denominator_expression(power=powers[edge.id], **options)
            product *= denominator
            selected = diagram.subgraph(edges=[edge.id])
            assert denominator == selected.denominator_expression(
                edge_powers=powers, **options
            )
            assert selected.edges[0].denominator_expression(**options) == (
                edge.denominator_expression(**options)
            )
        assert product == diagram.denominator_expression(
            edge_powers=powers, dimension=dimension, lmb=selected_basis
        )

for edge in diagram.edges:
    assert edge.particle.name == edge.particle_name
    assert edge.particle.mass_expression == 0
    assert edge.particle.width_parameter == "ZERO"
    if edge.is_external:
        assert edge.propagator is None
        try:
            edge.denominator_expression()
        except fk.DiagramError as error:
            assert "not an internal propagator" in str(error)
        else:
            raise AssertionError("external carriers do not supply propagators")
    else:
        assert isinstance(edge.propagator, fk.Propagator)
        assert edge.propagator.particle == edge.particle_name
        assert edge.denominator_expression(in_lmb=True) == basis.route_expression(
            edge.denominator_expression()
        )
    raw = edge.momentum_expression(dimension=4)
    assert isinstance(raw, TensorExpression) and raw.rank == 1
    assert raw == Q(edge.id, lorentz)
    assert edge.momentum_expression(dimension=4, in_lmb=True) == (
        basis.route_expression(raw)
    )
    for selected_basis in (basis, alternate):
        routed = edge.momentum_expression(dimension=4, lmb=selected_basis)
        assert routed == selected_basis.route_expression(raw)
        assert routed.rank == 1
        assert edge.momentum_signature(lmb=selected_basis).integer_coefficients() == (
            selected_basis.edge_signatures[edge.id].integer_coefficients()
        )
    assert edge.momentum_signature().integer_coefficients() == (
        basis.edge_signatures[edge.id].integer_coefficients()
    )

edge = diagram.internal_edges[0]
for expression in (
    edge.denominator_expression,
    edge.momentum_expression,
    edge.momentum_signature,
):
    try:
        expression(lmb=foreign_basis)
    except fk.DiagramError as error:
        assert "momentum basis belongs to a different diagram" in str(error)
    else:
        raise AssertionError("edge accepted another diagram's basis")

# A massive propagator must keep its symbolic mass rather than its numerical
# parameter value, and must agree with the diagram-level denominator.
massive = next(
    iter(
        model.process(
            ["scalar_1", "scalar_1"],
            ["scalar_1", "scalar_1"],
            vertex_allow=["V_3_SCALAR_111"],
        ).generate_diagrams(loops=0, max_vertices=2, progress=None)
    )
)
edge = massive.internal_edges[0]
assert edge.particle.mass_expression == model.particle("scalar_1").mass_expression
assert edge.particle.mass_expression != 0
assert edge.particle.width_parameter == "width_scalar_1"
assert edge.denominator_expression() == massive.denominator_expression()
assert (
    edge.particle.mass_expression
    in edge.denominator_expression().to_expression().get_all_symbols()
)
a, b, c, quadratic = S(
    "edge_test::a_", "edge_test::b_", "edge_test::c_", "edge_test::q_"
)
explicit = (
    edge.denominator_expression(dimension=4)
    .to_expression()
    .replace(S("gammalooprs::denom")(a, b, c, quadratic), quadratic)
)
momentum = edge.momentum_expression(dimension=4)
assert (
    explicit
    == (dot(momentum, momentum) - edge.particle.mass_expression**2).to_expression()
)
denominator = edge.denominator_expression(in_lmb=True)
del massive
assert edge.denominator_expression(in_lmb=True) == denominator

# Momentum conservation sets a one-point function's external momentum to zero.
# It must still be a vector, so it can participate in subsequent contractions.
tadpole = next(
    iter(
        model.process(
            ["scalar_0"],
            [],
            vertex_allow=["V_3_SCALAR_000"],
        ).generate_diagrams(
            loops=1,
            max_vertices=1,
            allow_self_loops=True,
            allow_zero_flow_edges=True,
            tadpoles=None,
            zero_snails=None,
            self_energy=None,
            maximum_bridges=None,
            progress=None,
        )
    )
)
zero_momentum = tadpole.external_edges[0].momentum_expression(in_lmb=True)
assert not zero_momentum and zero_momentum.rank == 1


# Custom names label independent coordinates without changing the routing.
# Exercise multiple loops, alternate bases, indexed vectors, and compact dots.
two_loop = next(
    iter(
        model.process(
            ["scalar_0"], ["scalar_0"], vertex_allow=["V_3_SCALAR_000"]
        ).generate_diagrams(
            loops=2, max_vertices=4, allow_self_loops=False, progress=None
        )
    )
)
for selected_basis in (basis, alternate, two_loop.loop_momentum_basis):
    loop_names = [
        TensorName.vector(f"named_routing::k{i}")
        for i in range(len(selected_basis.loop_edges))
    ]
    independent_external = [
        i
        for i, edge in enumerate(selected_basis.external_edges)
        if edge not in selected_basis.dependent_externals
    ]
    external_names = [
        TensorName.vector(f"named_routing::p{i}") for i in independent_external
    ]
    options = {"loop_momenta": loop_names, "external_momenta": external_names}
    zero = 0 * Q(selected_basis.loop_edges[0], mu)
    routed_zero = selected_basis.route_expression(zero, **options)
    assert isinstance(routed_zero, TensorExpression)
    assert not routed_zero and routed_zero.structure.slots == zero.structure.slots
    for edge, signature in selected_basis.edge_signatures.items():
        for port in (mu, lorentz):
            raw = Q(edge, port)
            expected = sum(
                (
                    coefficient * name(port).to_expression()
                    for coefficient, name in [
                        *zip(signature.loops, loop_names),
                        *(
                            (signature.external[i], name)
                            for i, name in zip(independent_external, external_names)
                        ),
                    ]
                ),
                E("0"),
            )
            named = selected_basis.route_expression(raw, **options)
            assert isinstance(named, TensorExpression)
            assert named.structure.slots == raw.structure.slots
            assert named.to_expression() == expected
            assert selected_basis.route_expression(named, **options) == named
            assert (
                selected_basis.route_expression(
                    selected_basis.route_expression(raw), **options
                )
                == named
            )
            plain = selected_basis.route_expression(raw.to_expression(), **options)
            assert isinstance(plain, Expression) and not isinstance(
                plain, TensorExpression
            )
            assert plain == expected
        compact = dot(Q(edge, lorentz), Q(edge, lorentz))
        named_dot = selected_basis.route_expression(compact, **options)
        assert isinstance(named_dot, TensorExpression) and named_dot.is_scalar
        assert not {
            S("gammalooprs::Q"),
            S("gammalooprs::K"),
            S("gammalooprs::P"),
        } & set(named_dot.to_expression().get_all_symbols())

    for invalid, message in (
        ({"loop_momenta": []}, "loop_momenta requires"),
        (
            {"external_momenta": [*external_names, loop_names[0]]},
            "external_momenta requires",
        ),
        (
            {
                "loop_momenta": [TensorName("named_routing::not_a_vector")]
                * len(loop_names)
            },
            "TensorName.vector",
        ),
    ):
        try:
            selected_basis.route_expression(
                Q(selected_basis.loop_edges[0], mu), **invalid
            )
        except ValueError as error:
            assert message in str(error)
        else:
            raise AssertionError(f"invalid momentum names accepted: {invalid}")
