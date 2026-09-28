"""Shared four-dimensional charge conjugation through the installed tensor API.

The ALOHA Weyl convention is C = -i gamma(2) gamma(0).
"""

from symbolica import E, S
from symbolica.community.spenso import (
    AUTO,
    Representation,
    TensorExpression,
    TensorName,
    chain,
)

spin = Representation.bis(4)
mink = Representation.mink(4)
i, j, k, mu, nu, dimension = S(
    "charge_python::i",
    "charge_python::j",
    "charge_python::k",
    "charge_python::mu",
    "charge_python::nu",
    "charge_python::D",
)
bis_i, bis_j, bis_k = (spin(index).to_expression() for index in (i, j, k))
mink_mu, mink_nu = (mink(index).to_expression() for index in (mu, nu))
# Construct C through the existing generic expression boundary. Its registered
# symbol owns real/antisymmetric normalization and the shared numerical data.
c_head, g_head, chain_head, incoming, outgoing = S(
    "spenso::charge_conjugation",
    "spenso::gamma",
    "spenso::chain",
    "spenso::in",
    "spenso::out",
)
charge = TensorExpression(c_head(bis_i, bis_j))
assert charge.spenso_conjugate().to_expression() == charge.to_expression()
assert c_head(bis_j, bis_i) == -charge.to_expression()
identity = TensorExpression.g(spin)(i, j)
square = TensorExpression(c_head(bis_i, bis_k) * c_head(bis_k, bis_j))
assert (
    square.simplify_gamma().to_expression().to_expression() == -identity.to_expression()
)

c = c_head(incoming, outgoing)
gamma_mu = g_head(incoming, outgoing, mink_mu)
transposed_mu = g_head(outgoing, incoming, mink_mu)
transposed_nu = g_head(outgoing, incoming, mink_nu)
sandwich = TensorExpression(chain_head(bis_i, bis_j, c, gamma_mu, c))
transposed_sandwich = TensorExpression(chain_head(bis_i, bis_j, c, transposed_mu, c))
assert sandwich.simplify_gamma().to_expression().to_expression() == chain_head(
    bis_i, bis_j, transposed_mu
)
assert (
    transposed_sandwich.simplify_gamma().to_expression().to_expression()
    == chain_head(bis_i, bis_j, gamma_mu)
)

# The public chain builder selects each matrix channel; no sewing indices are
# manufactured here. A product of transposes retains the original factor order.
scalar = E("2+3𝑖")
open_charge = TensorExpression(
    chain_head(spin.to_expression(), spin.to_expression(), c)
)
word = scalar * chain(
    spin(i),
    spin(j),
    open_charge,
    TensorExpression.gamma(4)(AUTO, AUTO, mu),
    TensorExpression.gamma(4)(AUTO, AUTO, nu),
    open_charge,
)
expected = -scalar * chain_head(bis_i, bis_j, transposed_mu, transposed_nu)
reduced = word.simplify_gamma().to_expression()
assert reduced.to_expression() == expected
assert (
    reduced.simplify_gamma().to_expression().to_expression() == reduced.to_expression()
)

# Slash dimension inference also accepts a momentum label followed by mink(D).
momentum = TensorName.vector("charge_python::P").to_expression()
mink_head = S("spenso::mink")
for dim in (E("4"), dimension):
    slash = g_head(incoming, outgoing, momentum(7, mink_head(dim)))
    compact = TensorExpression(chain_head(bis_i, bis_j, c, slash, c))
    if dim == E("4"):
        expected_slash = chain_head(
            bis_i, bis_j, g_head(outgoing, incoming, momentum(7, mink_head(dim)))
        )
        assert (
            compact.simplify_gamma().to_expression().to_expression() == expected_slash
        )
    else:
        assert (
            compact.simplify_gamma().to_expression().to_expression()
            == compact.to_expression()
        )
    ordinary = g_head(incoming, outgoing, mink_head(dim, mu))
    ordinary_word = TensorExpression(chain_head(bis_i, bis_j, c, ordinary, c))
    if dim == dimension:
        assert (
            ordinary_word.simplify_gamma().to_expression().to_expression()
            == ordinary_word.to_expression()
        )

# Evaluate original and rewritten words through the shared Weyl data. Explicit
# densification includes implicit sparse zeros. Assert axis labels before using
# row-major components, rather than assuming canonical storage order.
values = {}
for name, expression, ports in (
    ("C", charge, (bis_i, bis_j)),
    ("gamma", TensorExpression.gamma(4)(i, j, mu), (bis_i, bis_j, mink_mu)),
    ("square", square, (bis_i, bis_j)),
    ("sandwich", sandwich, (bis_i, bis_j, mink_mu)),
    ("transposed_sandwich", transposed_sandwich, (bis_i, bis_j, mink_mu)),
    ("word", word, (bis_i, bis_j, mink_mu, mink_nu)),
    ("reduced", reduced, (bis_i, bis_j, mink_mu, mink_nu)),
):
    network = expression.to_network()
    network.execute()
    tensor = network.result_tensor()
    assert tuple(slot.to_expression() for slot in tensor.structure.slots) == ports
    tensor.to_dense()
    values[name] = [complex(value) for value in tensor[:]]
    assert len(values[name]) == 4 ** len(ports)

gamma = values["gamma"]
for row in range(4):
    for column in range(4):
        # Derive C independently from existing gamma components, without copying
        # the newly registered charge-conjugation matrix entries.
        expected_c = -1j * sum(
            gamma[16 * row + 4 * inner + 2] * gamma[16 * inner + 4 * column]
            for inner in range(4)
        )
        assert values["C"][4 * row + column] == expected_c
        assert values["square"][4 * row + column] == -int(row == column)
        for lorentz in range(4):
            position = 16 * row + 4 * column + lorentz
            assert (
                values["sandwich"][position] == gamma[16 * column + 4 * row + lorentz]
            )
            assert values["transposed_sandwich"][position] == gamma[position]
            for second_lorentz in range(4):
                expected_word = -(2 + 3j) * sum(
                    gamma[16 * inner + 4 * row + lorentz]
                    * gamma[16 * column + 4 * inner + second_lorentz]
                    for inner in range(4)
                )
                position = 64 * row + 16 * column + 4 * lorentz + second_lorentz
                assert values["word"][position] == expected_word
                assert values["reduced"][position] == expected_word

print(
    "Shared C: explicit square, gamma/multiword/slash identities, dimension guards and Weyl components passed"
)
