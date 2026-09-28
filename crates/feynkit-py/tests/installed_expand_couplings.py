"""Coupling expansion preserves tensor ports and its precise Python return type."""

import json
from typing_extensions import assert_type

from symbolica import E, Expression, S
from symbolica.community import hep
from symbolica.community.spenso import Representation, TensorExpression, TensorName

model = hep.Model.standard_model()
coupling = S("UFO::GC_11")
analytic = model.coupling("GC_11").expression
assert_type(model.expand_couplings(coupling), Expression)
plain = model.expand_couplings(coupling)
assert type(plain) is Expression
assert plain == analytic

# Deliberately use noncanonical logical port order, including typed zero tensors.
lorentz = Representation.mink(4)
vector = TensorName.vector("coupling_test::v")
ordered = vector(lorentz("nu")) * vector(lorentz("mu"))
for source in (
    coupling * ordered,
    (coupling + 1) * ordered,
    (coupling - analytic) * ordered,
    0 * ordered,
    TensorExpression(coupling),
    ordered,
):
    assert isinstance(source, TensorExpression)
    before = source.to_expression()
    result = model.expand_couplings(source)
    assert_type(result, TensorExpression)
    assert type(result) is TensorExpression
    assert result.structure.slots == source.structure.slots
    assert result.to_expression() == model.expand_couplings(before)
    assert source.to_expression() == before
    assert model.expand_couplings(result) == result

cancelled = model.expand_couplings((coupling - analytic) * ordered)
assert isinstance(cancelled, TensorExpression)
assert not cancelled and cancelled.rank == 2
assert model.expand_couplings(coupling * ordered) == analytic * ordered

# A user-defined vanishing coupling must retain the rank of an open tensor too.
definition = json.loads(model.to_json())
next(c for c in definition["couplings"] if c["name"] == "GC_11")["expression"] = "0"
zero_model = hep.Model.from_json(json.dumps(definition))
zero = zero_model.expand_couplings(coupling * ordered)
assert isinstance(zero, TensorExpression)
assert not zero and zero.structure.slots == ordered.structure.slots

# Exercise a factored generated numerator without expanding its local factors.
diagram = model.process(["e-", "e+"], ["a", "a"]).generate_diagrams(progress=None)[0]
numerator = diagram.numerator_expression(in_lmb=True)
expanded = model.expand_couplings(numerator)
assert expanded.structure.slots == numerator.structure.slots
assert expanded.to_expression() == model.expand_couplings(numerator.to_expression())
assert not expanded.to_expression().matches(S("UFO::GC_3"))
assert model.expand_couplings(E("0")) == 0
print(
    "Coupling expansion: overloads, coefficients, ordered ports, zeros and diagram passed"
)
