"""Run the canonical four-loop notebook with the installed community host."""

import importlib.util
from collections import Counter
from pathlib import Path
from types import SimpleNamespace


def main():
    path = (
        Path(__file__).resolve().parents[3]
        / "examples/notebooks/feynkit/05_four_loop_numerator_marimo.py"
    )
    spec = importlib.util.spec_from_file_location("four_loop_showcase", path)
    notebook = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(notebook)
    _, initial = notebook.app.run()
    result = initial["generation_result"]
    assert result.report.completed and len(result) > 0
    assert result.report.retained_count == len(result.diagrams)
    print(result.report, flush=True)
    assert all(diagram.loop_count == 4 for diagram in result.diagrams)
    counts = [
        Counter(edge.particle_name for edge in diagram.internal_edges)
        for diagram in result.diagrams
    ]
    assert all(set(content) <= {"g", "ghG", "d", "d~"} for content in counts)
    # Select by particle content rather than historical dropdown positions.
    cases = {
        "ghost": lambda c: "ghG" in c and not ({"d", "d~"} & c.keys()),
        "gluon": lambda c: set(c) == {"g"},
        "quark": lambda c: bool({"d", "d~"} & c.keys()) and "ghG" not in c,
        "mixed": lambda c: bool({"d", "d~"} & c.keys()) and "ghG" in c,
    }
    selected = {
        label: next(i for i, content in enumerate(counts) if predicate(content))
        for label, predicate in cases.items()
    }
    selected["contract-regression"] = next(
        i
        for i, diagram in enumerate(result.diagrams)
        if diagram.name == "D2830" and counts[i] == {"g": 10}
    )
    for label, index in selected.items():
        _, values = notebook.app.run(
            defs={
                "generation_result": result,
                "model": initial["model"],
                "process": initial["process"],
                "topology_choice": SimpleNamespace(value=index),
            }
        )
        assert values["numerator"].structure.rank == 4
        assert values["projected_numerator"].structure.is_scalar
        assert values["reduced"].structure.is_scalar
        assert len(values["reduced"].expand().to_expression()) > 0
        print(
            f"Four-loop {label}: {result[index].name}, {dict(counts[index])}",
            flush=True,
        )


if __name__ == "__main__":
    main()
