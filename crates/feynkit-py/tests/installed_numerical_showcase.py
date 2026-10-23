"""Execute the numerical showcase with the installed native community host."""

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from symbolica import Float


def main():
    path = (
        Path(__file__).resolve().parents[3]
        / "examples/notebooks/feynkit/08_numerical_cff_ltd_marimo.py"
    )
    spec = importlib.util.spec_from_file_location("numerical_showcase", path)
    notebook = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(notebook)
    settings = {
        "batch_size": 1000,
        "batches": 4,
        "scale": 1.0,
        "alpha": 2.0,
        "seed": 2026,
    }

    for regime, mass, integrable in [
        ("On-shell external legs", 1.0, True),
        ("Spacelike external legs", 1.0, True),
        ("Spacelike external legs", 0.0, True),
        ("On-shell external legs", 0.25, False),
    ]:
        _, definitions = notebook.app.run(
            defs={
                "regime": SimpleNamespace(value=regime),
                "mass": SimpleNamespace(value=mass),
                "integration_settings": SimpleNamespace(value=settings),
            }
        )
        assert definitions["integrable"] == integrable
        # The physical kinematics give an independent oracle for the E/H
        # classification, including factors that originally had a minus sign.
        spacelike = regime == "Spacelike external legs"
        grid = np.linspace(-7 if spacelike else -2, 7 if spacelike else 2, 81)
        x, z = np.meshgrid(grid, grid)
        points = np.column_stack([x.ravel(), np.zeros(x.size), z.ravel()])
        surfaces = definitions["surface_evaluator"].evaluate(
            definitions["kinematic_inputs"](points)
        )
        crossing = {
            kind: any(
                k == kind and column.min() < 0 < column.max()
                for k, column in zip(definitions["surface_kinds"], surfaces.T)
            )
            for kind in ["E", "H"]
        }
        assert crossing == {"E": not spacelike and mass < 0.5, "H": spacelike}
        if not integrable:
            assert "summaries" not in definitions
            # At s=1, the equal-mass two-particle threshold is a circle of
            # radius sqrt(s/4-m²) about the common propagator center.
            points = definitions["channel_centers"] + [np.sqrt(0.25 - mass**2), 0, 0]
            surfaces = definitions["surface_evaluator"].evaluate(
                definitions["kinematic_inputs"](points)
            )
            e_columns = np.array(definitions["surface_kinds"]) == "E"
            assert np.min(np.abs(surfaces[:, e_columns])) < 1e-12
            continue

        # Compare the native estimator with the independently retained weights,
        # including the complete discrete/continuous density for channel runs.
        for row in definitions["summaries"]:
            weights = definitions["distributions"][row["method"]]
            assert row["evaluations"] == len(weights) == 4000
            assert np.isclose(row["estimate"], weights.mean(), rtol=1e-12)
            assert np.isclose(
                row["standard error"],
                weights.std(ddof=1) / np.sqrt(len(weights)),
                rtol=1e-11,
            )
            assert abs(row["deviation / error"]) < 8, (regime, row)

        inputs = definitions["kinematic_inputs"]([[1.0, 2.0, 3.0]])[0]
        evaluator = definitions["cff_evaluator"]
        precise = evaluator.evaluate_with_prec(
            [Float(value, decimal_digits=50) for value in inputs], 50
        )[0]
        reference = evaluator.evaluate_with_prec(
            [Float(value, decimal_digits=70) for value in inputs], 70
        )[0]
        assert precise.precision > 150
        assert abs(float((precise - reference) / reference)) < 1e-45
        values = definitions["comparison_evaluator"].evaluate([inputs])[0]
        assert np.isclose(values[0], values[1:].sum(), rtol=1e-10)

    # A finite-difference determinant verifies both coordinate-map Jacobians.
    for name, evaluator in definitions["coordinate_maps"].items():
        point = np.array([0.31, 0.62, 0.27, 1.2])
        jacobian = evaluator.evaluate([point])[0, 3]
        columns = []
        for axis in range(3):
            delta = np.zeros(4)
            delta[axis] = 1e-6
            columns.append(
                (
                    evaluator.evaluate([point + delta])[0, :3]
                    - evaluator.evaluate([point - delta])[0, :3]
                )
                / 2e-6
            )
        assert np.isclose(abs(np.linalg.det(columns)), jacobian, rtol=1e-7), name
    print(
        "Showcase: native estimates, OneLOop, precision, Jacobians, E/H contours and threshold gate passed"
    )


if __name__ == "__main__":
    main()
