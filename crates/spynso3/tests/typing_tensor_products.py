"""Tensor multiplication typing against an installed Symbolica community package."""

from typing import assert_type

from symbolica import Expression
from symbolica.community.tensor import Tensor, TensorExpression, TensorNetwork


def check_tensor_products(
    tensor: TensorExpression,
    scalar: Expression,
    components: Tensor,
    network: TensorNetwork,
) -> None:
    assert_type(scalar * tensor, TensorExpression)
    assert_type(tensor * scalar, TensorExpression)
    assert_type(tensor.__mul__(scalar), TensorExpression)
    assert_type(tensor * tensor, TensorExpression)
    assert_type(tensor * 2, TensorExpression)
    assert_type(tensor * 2.0, TensorExpression)
    assert_type(tensor * 1j, TensorExpression)
    assert_type(tensor * components, TensorNetwork)
    assert_type(tensor * network, TensorNetwork)
    assert_type(components * tensor, TensorNetwork)
    assert_type(network * tensor, TensorNetwork)
