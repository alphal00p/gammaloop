import os
from pathlib import Path

import numpy as np
from gammaloop import GammaLoopAPI


def main() -> None:
    example_dir = Path(__file__).resolve().parent
    run_card = example_dir / "run.toml"
    state_dir = example_dir / "state"

    # The typed Python API is the preferred entry point for stateful workflows.
    # `run(...)` remains available as the generic fallback for CLI-only commands.
    api = GammaLoopAPI(
        state_folder=state_dir,
        boot_commands_path=run_card,
        clean_state=True,
    )
    # The example card includes a reusable block that shows the new named
    # process-setting display commands through the Python `run(...)` bridge.
    api.run("run display_named_settings_examples")

    info = api.get_integrand_info()
    assert info.kind == "cross section"
    assert info.graph_count == 2
    assert info.graph_group_count == len(info.graph_groups)
    assert all(
        sum(graph.is_master for graph in group.graphs) == 1
        for group in info.graph_groups
    )

    point = np.array([0.17, 0.31, 0.53, 0.23, 0.41, 0.67], dtype=float)
    result = api.evaluate_sample(point.tolist(), return_events=True)
    assert result.parameterization_jacobian is not None
    assert result.stability_results
    assert result.event_groups

    batch_points = np.array(
        [
            [0.17, 0.31, 0.53, 0.23, 0.41, 0.67],
            [0.11, 0.29, 0.47, 0.19, 0.37, 0.59],
        ],
        dtype=float,
    )
    batch_result = api.evaluate_samples(batch_points, return_events=True)
    assert len(batch_result.samples) == 2
    assert all(sample.event_groups for sample in batch_result.samples)
    leading_jet_pt = batch_result.observables["leading_jet_pt_hist"]
    jet_count = batch_result.observables["jet_count_hist"]
    assert len(leading_jet_pt.bins) == 8
    assert len(jet_count.bins) == 6
    assert leading_jet_pt.sample_count == 2
    assert jet_count.sample_count == 2

    # Momentum-space input has one flattened (px, py, pz) triplet per loop;
    # neither loop energies nor external momenta belong in this array.
    loop_momenta = [
        0.11,
        -0.07,
        0.19,  # k1 = (px, py, pz)
        -0.13,
        0.05,
        0.29,  # k2 = (px, py, pz)
    ]
    momentum_result = api.evaluate_sample(loop_momenta, momentum_space=True)

    print("== x-space evaluate_sample ==\n")
    print(result)
    print(f"integrand_result = {result.integrand_result}")
    print(f"observables = {result.observables}")

    print("\n== x-space evaluate_samples ==\n")
    print(batch_result)
    for idx, sample in enumerate(batch_result.samples):
        print(f"sample[{idx}] integrand_result = {sample.integrand_result}")
    print(f"batch observables = {batch_result.observables}")

    print("\n== momentum-space evaluate_sample ==\n")
    print(momentum_result)


if __name__ == "__main__":
    main()
