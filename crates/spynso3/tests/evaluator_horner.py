"""Exercise the shared native wrappers without installing a Python package."""

import inspect
import itertools

# Supplied by the Rust host without importing a second native module.
core = globals()["core"]
tensor = globals()["tensor"]

assert inspect.signature(tensor.Tensor.evaluator) == inspect.signature(core.Expression.evaluator)
x, y = core.S("horner_tensor_test::x", "horner_tensor_test::y")
name = tensor.TensorName("horner_tensor_test::values")(tensor.Representation.euc(2))
values = tensor.Tensor.dense(name, [(x + y) ** 3 + x * y, x**2 - y**2])
points = [[2.0, 3.0], [-1.0, 2.0]]
expected = [[131.0, -5.0], [-1.0, -3.0]]

for direct_translation, horner_use_cpe in itertools.product((False, True), repeat=2):
    options = {
        "iterations": 4,
        "n_cores": 1,
        "jit_compile": False,
        "direct_translation": direct_translation,
        "horner_use_cpe": horner_use_cpe,
    }
    evaluator = values.evaluator([x, y], **options)
    scalar = core.Expression.evaluator_multiple(values[:], [x, y], **options)
    assert isinstance(evaluator.scalar_evaluator, core.Evaluator)
    scalar_rows = scalar.evaluate(points)
    if hasattr(scalar_rows, "tolist"):
        scalar_rows = scalar_rows.tolist()
    assert scalar_rows == expected
    for result, row in zip(evaluator.evaluate(points), expected, strict=True):
        assert result[:] == row
        assert result.structure == values.structure
        assert result.expression().to_expression() == values.expression().to_expression()

# Use the upstream candidate-scoring regression polynomial to make forwarding
# observable: accepting the keyword but always passing False must fail here.
params = list(core.S("horner_forward::a", "horner_forward::b", "horner_forward::c"))
different_modes = False
for a, b, c in itertools.permutations(params):
    polynomial = (
        b * c**3 - 3 * b**2 - a * c**2 + 2 * a**2 * c
        - 2 * a**2 * b**2 - 3 * a**3 * c + 2 * a**3 * c**2
        + 3 * a**3 * b**2 * c - 2 * a**3 * b**2 * c**2
    )
    source = tensor.Tensor.dense(name, [polynomial, polynomial + 1])
    instructions = []
    for horner_use_cpe in (False, True):
        options = {
            "iterations": 6,
            "n_cores": 1,
            "jit_compile": False,
            "horner_use_cpe": horner_use_cpe,
        }
        wrapped = source.evaluator(params, **options)
        scalar = core.Expression.evaluator_multiple(source[:], params, **options)
        expected_instructions = scalar.get_instructions().instructions
        assert wrapped.scalar_evaluator.get_instructions().instructions == expected_instructions
        instructions.append(expected_instructions)
    different_modes |= instructions[0] != instructions[1]
assert different_modes
