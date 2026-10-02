"""Verify tensor-owned replacements, also usable in a tensor-only host."""

from symbolica import E, S
from symbolica.community.tensor import (
    AUTO,
    BroadcastFunction,
    FactorProjector,
    Nc,
    PortPattern,
    Representation,
    TensorExpression,
    TensorLibrary,
    TensorName,
    TensorPattern,
    chain,
    dot,
    trace,
)


def main():
    mink, spin = Representation.mink(4), Representation.bis(4)
    dimension, degree, index, x, y, factors = S(
        "removal_audit::D_",
        "removal_audit::k_",
        "removal_audit::i_",
        "removal_audit::x_",
        "removal_audit::y_",
        "removal_audit::factors___",
    )
    arbitrary_rep = S("removal_audit::rep_")
    fundamental = Representation.cof(3)
    fundamental_pattern = Representation.cof(dimension)
    replacements = (
        (
            mink("mu").to_expression(),
            PortPattern.exact(Representation.mink(dimension), index),
        ),
        (fundamental.casimir(), fundamental_pattern.casimir(degree)),
        (fundamental.dynkin_index(), fundamental_pattern.dynkin_index(degree)),
        (fundamental.casimir(), TensorPattern.casimir(degree, arbitrary_rep)),
        (fundamental.dynkin_index(), TensorPattern.dynkin_index(degree, arbitrary_rep)),
        (BroadcastFunction.conj()(S("removal_audit::z")), BroadcastFunction.conj()(x)),
        (
            TensorExpression.g(mink)("mu", "nu").to_expression(),
            TensorPattern(TensorName.g(), ports=[factors]),
        ),
    )
    for target, pattern in replacements:
        assert target.replace(pattern, E("7")) == E("7"), (target, pattern)

    assert mink.to_expression().get_head() == S("spenso::mink")
    assert fundamental.casimir().get_head() == S("spenso::cas")
    assert fundamental.dynkin_index().get_head() == S("spenso::idx")
    assert Nc == S("spenso::Nc")
    assert Nc.evaluate({}) == 3

    p = TensorName.vector("removal_audit::p")(mink)
    slash = (
        (TensorExpression.dirac_gamma(4)(AUTO, AUTO, "mu") * p("mu"))
        .contract()
        .to_expression()
    )
    word = chain(spin("i"), spin("j"), slash)
    expected = S("spenso::chain")(
        spin("i").to_expression(),
        spin("j").to_expression(),
        S("spenso::gamma")(S("spenso::in"), S("spenso::out"), p.to_expression()),
    )
    assert word.to_expression() == expected

    gamma = TensorExpression.dirac_gamma(4)
    gamma_factors = [gamma(AUTO, AUTO, label) for label in ("mu", "nu")]
    projected_trace = trace(spin, FactorProjector.symmetric(*gamma_factors))
    compact = (
        (dot(p, p), TensorPattern.dot(x, y)),
        (word, TensorPattern.chain(x, y, factors)),
        (trace(spin, *gamma_factors), TensorPattern.trace(arbitrary_rep, factors)),
        (
            projected_trace,
            TensorPattern.trace(arbitrary_rep, TensorPattern.symmetric(factors)),
        ),
        (TensorExpression.gamma0(4)("i", "j"), TensorPattern.gamma0(dimension, x, y)),
        (
            TensorExpression.charge_conjugation(4)("i", "j"),
            TensorPattern.charge_conjugation(dimension, x, y),
        ),
        (
            TensorExpression.levi_civita(mink)("i", "j", "k", "l"),
            TensorPattern(TensorName.levi_civita(), ports=[factors]),
        ),
    )
    for target, pattern in compact:
        assert target.to_expression().replace(pattern, E("7")) == E("7"), (
            target,
            pattern,
        )
    assert TensorName.charge_conjugation()(spin, spin).to_expression() == E("0")
    assert TensorName.levi_civita()(mink, mink, mink, mink).to_expression() == E("0")
    assert TensorExpression.charge_conjugation(4)(
        "i", "j"
    ) == TensorLibrary.hep_lib_atom()["spenso::charge_conjugation"].expression()(
        "i", "j"
    )
    assert projected_trace.expand().to_expression() == projected_trace.to_expression()
    assert (
        projected_trace.expand_projectors().expand_projectors()
        == projected_trace.expand_projectors()
    )
    print(
        "Tensor-owned construction, wildcard patterns, metadata and projector checks passed"
    )


if __name__ == "__main__":
    main()
