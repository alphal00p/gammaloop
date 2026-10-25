"""Exercise complete and selected numerator interfaces in the combined native host."""

import json


def check_ordered_numerator_ports(
    *,
    core,
    feynkit,
    spenso,
    model,
    ordered_left,
    ordered_right,
    ordered_diagram,
    inconsistent_diagram,
    repaired_diagram,
    replaced_diagram,
):
    expected = ordered_left * ordered_right
    expected_axes = ordered_left.axes + ordered_right.axes
    complete = ordered_diagram.numerator_expression()
    assert complete == expected
    assert complete.axes == expected_axes
    assert core.E(ordered_diagram.numerator) == complete.to_expression()

    source = ordered_diagram.to_json()
    stored = json.loads(source)
    assert len(stored["num_interface"]) == len(expected_axes)
    assert all(
        axis.to_expression() == core.E(port)
        for axis, port in zip(expected_axes, stored["num_interface"])
    )
    restored = feynkit.FeynmanDiagram.from_json(model, source)
    assert restored.numerator_expression() == complete
    assert restored.numerator_expression().axes == expected_axes

    full = ordered_diagram.subgraph(edges=[edge.id for edge in ordered_diagram.edges])
    left = ordered_diagram.subgraph(nodes=[0])
    right = ordered_diagram.subgraph(nodes=[1])
    empty = ordered_diagram.subgraph()
    assert full.numerator_expression() == complete
    assert left.numerator_expression() == ordered_left
    assert right.numerator_expression() == ordered_right
    assert ordered_diagram.numerator_expression(without=right) == ordered_left
    assert full.numerator_expression(without=left) == ordered_right
    assert left.numerator_expression(without=left) == spenso.TensorExpression(1)
    assert empty.numerator_expression() == spenso.TensorExpression(1)

    # Explicit exclusion requests retain the existing local-fragment traversal,
    # including the empty exclusion. This traverses right before left in this graph.
    reconstructed = ordered_diagram.numerator_expression(without=empty)
    assert reconstructed.axes == ordered_right.axes + ordered_left.axes
    assert reconstructed.to_expression() == complete.to_expression()

    assert replaced_diagram.numerator_expression() == ordered_left
    assert repaired_diagram.numerator_expression() == 2 * ordered_right
    try:
        inconsistent_diagram.numerator_expression()
    except feynkit.DiagramError:
        pass
    else:
        raise AssertionError(
            "complete getter accepted a stale aggregate after map_data"
        )
