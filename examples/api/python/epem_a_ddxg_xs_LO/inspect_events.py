import os
from pathlib import Path

import numpy as np
from gammaloop import GammaLoopAPI


def main() -> None:
    example_dir = Path(__file__).resolve().parent
    run_card = example_dir / "run.toml"
    state_dir = example_dir / "state"

    # Temporary until the proper LU numerator / theta / tree-denominator path
    # is fixed.
    os.environ.setdefault("GL_LU_E2E_HACK", "1")

    # The typed Python API is the preferred entry point for stateful workflows.
    # `run(...)` remains available as the generic fallback for CLI-only commands.
    api = GammaLoopAPI(
        state_folder=state_dir,
        boot_commands_path=run_card,
        fresh_state=True,
    )
    # The example card includes a reusable block that shows the new named
    # process-setting display commands through the Python `run(...)` bridge.
    api.run("run display_named_settings_examples")

    point = np.array([0.17, 0.31, 0.53, 0.23, 0.41, 0.67], dtype=float)
    result = api.evaluate_sample(point.tolist())

    batch_points = np.array(
        [
            [0.17, 0.31, 0.53, 0.23, 0.41, 0.67],
            [0.11, 0.29, 0.47, 0.19, 0.37, 0.59],
        ],
        dtype=float,
    )
    batch_result = api.evaluate_samples(batch_points)

    momentum_result = api.evaluate_sample(
        [0.11, -0.07, 0.19, -0.13, 0.05, 0.29],
        momentum_space=True,
    )

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
