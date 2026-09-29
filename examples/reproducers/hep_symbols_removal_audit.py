"""Report tensor API removal gaps without asserting that bugs should persist.

Run in the same combined community host as hep_symbols_primitives.py. Positive
checks exercise existing replacements; observations describe missing surfaces.
"""

import json

from symbolica import E, S
from symbolica.community.hep import Symbols
from symbolica.community.tensor import (
    AUTO,
    BroadcastFunction,
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
    fundamental = Representation.cof(3)
    fundamental_pattern = Representation.cof(dimension)
    replacements = (
        (
            mink("mu").to_expression(),
            PortPattern.exact(Representation.mink(dimension), index),
        ),
        (fundamental.casimir(), fundamental_pattern.casimir(degree)),
        (fundamental.dynkin_index(), fundamental_pattern.dynkin_index(degree)),
        (BroadcastFunction.conj()(S("removal_audit::z")), BroadcastFunction.conj()(x)),
        (
            TensorExpression.g(mink)("mu", "nu").to_expression(),
            TensorPattern(TensorName.g(), ports=[factors]),
        ),
    )
    for target, pattern in replacements:
        assert target.replace(pattern, E("7")) == E("7"), (target, pattern)

    assert mink.to_expression().get_head() == Symbols.lorentz
    assert fundamental.casimir().get_head() == Symbols.casimir
    assert fundamental.dynkin_index().get_head() == Symbols.color_index
    assert Symbols.color_number.evaluate({}) == 3

    p = TensorName.vector("removal_audit::p")(mink)
    slash = (
        (TensorExpression.dirac_gamma(4)(AUTO, AUTO, "mu") * p("mu"))
        .contract()
        .to_expression()
    )
    word = chain(spin("i"), spin("j"), slash)
    expected = Symbols.gamma_chain(
        spin("i").to_expression(),
        spin("j").to_expression(),
        Symbols.dirac_gamma(Symbols.chain_in, Symbols.chain_out, p.to_expression()),
    )
    assert word.to_expression() == expected

    observations = {}

    def observe(label, operation):
        try:
            value = operation()
            result = {"type": type(value).__name__, "value": str(value)}
            if isinstance(value, TensorExpression):
                result.update(rank=value.rank, is_zero=value.to_expression() == E("0"))
            observations[label] = result
        except (TypeError, ValueError, AttributeError) as error:
            observations[label] = {"error": type(error).__name__, "message": str(error)}

    observe("dot_operand_wildcards", lambda: dot(x, y))
    observe("chain_sequence_wildcard", lambda: chain(spin("i"), spin("j"), factors))
    observe("trace_sequence_wildcard", lambda: trace(spin, factors))
    observe(
        "pattern_auto_ports",
        lambda: TensorPattern.dirac_gamma(dimension, AUTO, AUTO, index),
    )
    observe(
        "generic_unresolved_charge_conjugation",
        lambda: TensorName("spenso::charge_conjugation")(spin, spin),
    )
    observe(
        "generic_unresolved_epsilon",
        lambda: TensorName("spenso::epsilon")(mink, mink, mink, mink),
    )
    observe(
        "library_unresolved_charge_conjugation",
        lambda: TensorLibrary.hep_lib_atom()["spenso::charge_conjugation"].expression(),
    )
    raw_factors = [
        Symbols.dirac_gamma(
            Symbols.chain_in, Symbols.chain_out, mink(label).to_expression()
        )
        for label in ("mu", "nu")
    ]
    projector = Symbols.symmetric(*raw_factors)
    observe("standalone_projector", lambda: TensorExpression(projector))
    observe("projector_in_typed_chain", lambda: chain(spin("i"), spin("j"), projector))
    projected_trace = TensorExpression(Symbols.trace(spin.to_expression(), projector))
    observe("complete_projected_trace", lambda: projected_trace)
    observe(
        "ordinary_expand_preserves_projector",
        lambda: (
            projected_trace.expand().to_expression() == projected_trace.to_expression()
        ),
    )
    print(
        "10 existing replacement, wildcard, metadata, and compact-chain checks passed"
    )
    print(json.dumps(observations, indent=2))


if __name__ == "__main__":
    main()
