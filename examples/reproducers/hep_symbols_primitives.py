"""Audit HEP symbol construction against the existing typed tensor API.

Run with a Symbolica community host exposing community.tensor.
Conversions below are equality-oracle boundaries, not construction requirements.
Canonical raw heads below serve as independent expression oracles.
"""

from symbolica import S
from symbolica.community.tensor import (
    AUTO,
    BroadcastFunction,
    Representation,
    TensorExpression,
    TensorLibrary,
    TensorName,
    chain,
    dot,
    trace,
)


def main():
    checks = 0

    def same(value, expected, *, rank=None):
        nonlocal checks
        if rank is not None:
            assert isinstance(value, TensorExpression)
            assert value.rank == rank
            value = value.to_expression()
        assert value == expected, (value, expected)
        checks += 1

    mink = Representation.mink(4)
    spin = Representation.bis(4)
    fund = Representation.cof(3)
    adj = Representation.coad(8)
    index = S("hep_symbols_audit::i")
    for representation, head, dimension in (
        (mink, S("spenso::mink"), 4),
        (spin, S("spenso::bis"), 4),
        (fund, S("spenso::cof"), 3),
        (adj, S("spenso::coad"), 8),
    ):
        same(representation.to_expression(), head(dimension))
        same(representation(index).to_expression(), head(dimension, index))

    mu, nu = (mink(index).to_expression() for index in ("mu", "nu"))
    a, b = (spin(index).to_expression() for index in ("a", "b"))
    same(TensorExpression.g(mink)("mu", "nu"), S("spenso::g")(mu, nu), rank=2)
    same(
        TensorExpression.dirac_gamma(4)("a", "b", "mu"),
        S("spenso::gamma")(a, b, mu),
        rank=3,
    )
    same(
        TensorExpression.color_f(8)("a", "b", "c"),
        S("spenso::f")(*(adj(index).to_expression() for index in ("a", "b", "c"))),
        rank=3,
    )
    same(fund.casimir(), S("spenso::cas")(2, fund.to_expression()))
    same(fund.dynkin_index(), S("spenso::idx")(2, fund.to_expression()))

    p = TensorName.vector("hep_symbols_audit::p")(mink)
    q = TensorName.vector("hep_symbols_audit::q")(mink)
    same(dot(p, q), S("spenso::dot")(p.to_expression(), q.to_expression()), rank=0)
    same(BroadcastFunction.conj()(p), S("spenso::conj")(p.to_expression()), rank=1)

    factors = (
        TensorExpression.dirac_gamma(4)(AUTO, AUTO, "mu"),
        TensorExpression.dirac_gamma(4)(AUTO, AUTO, "nu"),
    )
    raw_factors = (
        S("spenso::gamma")(S("spenso::in"), S("spenso::out"), mu),
        S("spenso::gamma")(S("spenso::in"), S("spenso::out"), nu),
    )
    same(
        chain(spin("a"), spin("b"), *factors),
        S("spenso::chain")(a, b, *raw_factors),
        rank=4,
    )
    same(
        trace(spin, *factors),
        S("spenso::trace")(spin.to_expression(), S("spenso::cyclic")(*raw_factors)),
        rank=2,
    )

    library = TensorLibrary.hep_lib_atom()
    for key, head in (
        ("spenso::gamma0", S("spenso::gamma0")),
        ("spenso::charge_conjugation", S("spenso::charge_conjugation")),
    ):
        reference = library[key].expression()
        same(reference("a", "b"), head(a, b), rank=2)

    # An explicitly indexed epsilon works through the generic tensor primitive.
    # The named factory additionally supplies distinct unresolved ports.
    slots = [mink(index) for index in ("mu", "nu", "rho", "sigma")]
    epsilon = TensorName("spenso::epsilon")(*slots)
    same(
        epsilon, S("spenso::epsilon")(*(slot.to_expression() for slot in slots)), rank=4
    )

    # Introspection and raw patterns have a different boundary from construction.
    for name, head in (
        (TensorName.g(), S("spenso::g")),
        (TensorName.dirac_gamma(), S("spenso::gamma")),
        (TensorName.color_f(), S("spenso::f")),
    ):
        same(name.to_expression(), head)
    same(BroadcastFunction.conj().to_expression(), S("spenso::conj"))

    print(f"{checks} typed construction, interface, and head comparisons passed")


if __name__ == "__main__":
    main()
