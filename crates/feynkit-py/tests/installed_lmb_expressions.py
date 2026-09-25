"""Exercise optional expression routing in an installed FeynKit community host.

Pass the community module name as the first argument (``hep`` for that host).
"""

import importlib
import json
import sys
from pathlib import Path

from symbolica import E, S, Expression
from symbolica.community.spenso import Representation, TensorExpression, TensorName, dot

fk = importlib.import_module(
    f"symbolica.community.{sys.argv[1] if len(sys.argv) > 1 else 'feynkit'}"
)
model = fk.Model(Path(__file__).parent / "fixtures/scalars_2p_3p.json")
diagrams = (
    fk.Process(model, ["scalar_0"], ["scalar_0", "scalar_0"])
    .generate_diagrams(
        loops=1, max_vertices=3, vertex_allow=["V_3_SCALAR_000"], allow_self_loops=True
    )
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
    cancelling = (
        indexed - selected_basis.route_expression(indexed.to_expression())
    ).expand_num()
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
    assert cancelled == 0
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
    assert getattr(empty, name)(in_lmb=True) == 1
    assert getattr(empty, name)(lmb=alternate) == 1

without_region = diagram.numerator_expression(without=region)
assert without_region == 1
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


# Custom names label independent coordinates without changing the routing.
# Exercise multiple loops, alternate bases, indexed vectors, and compact dots.
two_loop = next(
    iter(
        fk.Process(model, ["scalar_0"], ["scalar_0"]).generate_diagrams(
            loops=2,
            max_vertices=4,
            vertex_allow=["V_3_SCALAR_000"],
            allow_self_loops=False,
            progress=None,
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
    options = dict(loop_momenta=loop_names, external_momenta=external_names)
    zero = 0 * Q(selected_basis.loop_edges[0], mu)
    routed_zero = selected_basis.route_expression(zero, **options)
    assert isinstance(routed_zero, TensorExpression)
    assert routed_zero == 0 and routed_zero.structure.slots == zero.structure.slots
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
        } & set(named_dot.get_all_symbols())

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
