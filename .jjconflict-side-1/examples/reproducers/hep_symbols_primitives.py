"""Audit HEP symbol construction against the existing typed tensor API.

Run with a combined Symbolica community host exposing community.tensor and hep.
Conversions below are equality-oracle boundaries, not construction requirements.
The accompanying architecture note records unsupported cases separately.
"""

from symbolica import S
from symbolica.community.hep import Symbols
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
        (mink, Symbols.lorentz, 4),
        (spin, Symbols.spinor, 4),
        (fund, Symbols.color_fundamental, 3),
        (adj, Symbols.color_adjoint, 8),
    ):
        same(representation.to_expression(), head(dimension))
        same(representation(index).to_expression(), head(dimension, index))

    mu, nu = (mink(index).to_expression() for index in ("mu", "nu"))
    a, b = (spin(index).to_expression() for index in ("a", "b"))
    same(TensorExpression.g(mink)("mu", "nu"), Symbols.metric(mu, nu), rank=2)
    same(TensorExpression.gamma(4)("a", "b", "mu"), Symbols.gamma(a, b, mu), rank=3)
    same(
        TensorExpression.f(8)("a", "b", "c"),
        Symbols.color_f(*(adj(index).to_expression() for index in ("a", "b", "c"))),
        rank=3,
    )
    same(fund.casimir(), Symbols.casimir(2, fund.to_expression()))
    same(fund.dynkin_index(), Symbols.color_index(2, fund.to_expression()))

    p = TensorName.vector("hep_symbols_audit::p")(mink)
    q = TensorName.vector("hep_symbols_audit::q")(mink)
    same(dot(p, q), Symbols.dot(p.to_expression(), q.to_expression()), rank=0)
    same(BroadcastFunction.conj()(p), Symbols.conjugate(p.to_expression()), rank=1)

    factors = (
        TensorExpression.gamma(4)(AUTO, AUTO, "mu"),
        TensorExpression.gamma(4)(AUTO, AUTO, "nu"),
    )
    raw_factors = (
        Symbols.gamma(Symbols.chain_in, Symbols.chain_out, mu),
        Symbols.gamma(Symbols.chain_in, Symbols.chain_out, nu),
    )
    same(
        chain(spin("a"), spin("b"), *factors),
        Symbols.gamma_chain(a, b, *raw_factors),
        rank=4,
    )
    same(
        trace(spin, *factors),
        Symbols.trace(spin.to_expression(), Symbols.cyclic(*raw_factors)),
        rank=2,
    )

    library = TensorLibrary.hep_lib_atom()
    for key, head in (
        ("spenso::gamma0", Symbols.gamma_zero),
        ("spenso::charge_conjugation", Symbols.charge_conjugation),
    ):
        reference = library[key].expression()
        same(reference("a", "b"), head(a, b), rank=2)

    # An explicitly indexed epsilon works through the generic tensor primitive.
    # A named epsilon factory and safe unresolved construction are still missing.
    slots = [mink(index) for index in ("mu", "nu", "rho", "sigma")]
    epsilon = TensorName("spenso::epsilon")(*slots)
    same(
        epsilon, Symbols.levi_civita(*(slot.to_expression() for slot in slots)), rank=4
    )

    # Introspection and raw patterns have a different boundary from construction.
    for name, head in (
        (TensorName.g(), Symbols.metric),
        (TensorName.gamma(), Symbols.gamma),
        (TensorName.f(), Symbols.color_f),
    ):
        same(name.to_expression(), head)
    same(BroadcastFunction.conj().to_expression(), Symbols.conjugate)

    print(f"{checks} typed construction, interface, and head comparisons passed")


if __name__ == "__main__":
    main()
