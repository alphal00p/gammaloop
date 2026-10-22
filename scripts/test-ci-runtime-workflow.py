#!/usr/bin/env python3
"""Exercise Actions runtime selection and verdicts with a cached fixture runner."""

import argparse
import json
import os
import subprocess
import tempfile
from pathlib import Path

import yaml


def execute(script, env):
    return subprocess.run(
        ["bash", "-e", "-o", "pipefail", "-c", script],
        env=env,
        capture_output=True,
        check=False,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("launcher", type=Path)
    parser.add_argument("--report-dir", type=Path, required=True)
    args = parser.parse_args()
    args.report_dir.mkdir(parents=True, exist_ok=True)
    with open(".github/workflows/nix.yml") as handle:
        jobs = yaml.safe_load(handle)["jobs"]
    assert "runtime" in jobs["flake-checks"]["needs"]
    select = next(
        step["run"]
        for step in jobs["prepare"]["steps"]
        if step.get("id") == "runtime-matrix"
    )
    aggregate = jobs["flake-checks"]["steps"][0]["run"]
    body = next(
        step["run"] for step in jobs["runtime"]["steps"] if step.get("id") == "execute"
    )
    with tempfile.TemporaryDirectory(prefix="gammaloop-workflow-contract-") as tmp:
        root = Path(tmp)
        (root / "bin").mkdir()
        matrix = root / "matrix.json"
        matrix.write_text(
            json.dumps(
                {
                    "include": [
                        {"name": "fixture", "package": "nix-ci-check-runtime-fixture"}
                    ]
                }
            )
        )
        fake_nix = root / "bin/nix"
        fake_nix.write_text('#!/usr/bin/env bash\nprintf "%s\\n" "$FIXTURE_MATRIX"\n')
        fake_nix.chmod(0o755)
        output = root / "output"
        env = dict(
            os.environ,
            PATH=str(root / "bin") + ":" + os.environ["PATH"],
            GITHUB_OUTPUT=str(output),
            LICENSE="public-workflow-fixture",
            FIXTURE_MATRIX=str(matrix),
        )
        result = execute(select.replace("${{ github.ref }}", "refs/heads/main"), env)
        assert result.returncode == 0, "Runtime selection failed"
        assert "any-jobs=true" in output.read_text(), (
            "Cached pure matrix suppressed runtime rows"
        )
        selected = json.loads(output.read_text().split("matrix=", 1)[1].splitlines()[0])
        assert len(selected["include"]) == 1
        verdict_env = dict(
            env,
            PREPARE_RESULT="success",
            BUILD_RESULT="skipped",
            ANY_JOBS="false",
            ANY_RUNTIME_JOBS="true",
        )
        assert (
            execute(aggregate, dict(verdict_env, RUNTIME_RESULT="success")).returncode
            == 0
        )
        assert (
            execute(aggregate, dict(verdict_env, RUNTIME_RESULT="failure")).returncode
            != 0
        )
        # The same already-built executable is invoked by the actual Actions
        # execution step twice; no Nix build or Symbolica process is mocked here.
        runner = root / "runtime-launcher/bin"
        runner.mkdir(parents=True)
        (runner / "nix-ci-check-runtime-fixture").symlink_to(args.launcher.resolve())
        env.update(
            RUNNER_TEMP=str(root),
            RUNNER_PACKAGE="nix-ci-check-runtime-fixture",
            CHECK_NAME="fixture",
            SYMBOLICA_LICENSE_SIGNED="public-workflow-fixture",
            CI_RUNTIME_EXPECTED="public-workflow-fixture",
            CI_RUNTIME_COUNTER=str(root / "counter"),
            HOME=str(root / "home"),
            XDG_CONFIG_HOME=str(root / "config"),
            GITHUB_SHA="fixture-commit",
            GITHUB_SERVER_URL="https://example.invalid",
            GITHUB_REPOSITORY="fixture/repo",
            GITHUB_RUN_ID="1",
            GITHUB_RUN_ATTEMPT="1",
        )
        env.pop("SYMBOLICA_LICENSE", None)
        for _ in range(2):
            result = execute('cd "$RUNNER_TEMP"\n' + body, env)
            assert result.returncode == 0, (
                "Actions did not execute cached runtime fixture"
            )
        assert (root / "counter").read_text().splitlines() == ["executed", "executed"]
        env["CI_RUNTIME_EXIT"] = "37"
        result = execute('cd "$RUNNER_TEMP"\n' + body, env)
        assert result.returncode == 37, (
            "Actions replaced runtime failure with report status"
        )
    (args.report_dir / "receipt.json").write_text(
        json.dumps(
            {
                "empty_pure_matrix_keeps_runtime": True,
                "runtime_failure_fails_aggregate": True,
                "cached_launcher_executions": 2,
                "runtime_exit_37_preserved": True,
            },
            indent=2,
        )
        + "\n"
    )
    print(
        "Workflow contract passed: cached runtime selection, two executions, failure aggregate and exit 37."
    )


if __name__ == "__main__":
    main()
