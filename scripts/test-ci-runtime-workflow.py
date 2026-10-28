#!/usr/bin/env python3
"""Exercise Actions runtime selection and verdicts with a cached fixture runner."""

import argparse
import json
import os
import shutil
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
        workflow = yaml.load(handle, Loader=yaml.BaseLoader)
    assert workflow["on"] == {"push": {"branches": ["main"]}, "workflow_dispatch": ""}
    with open(".github/workflows/docs-pages.yml") as handle:
        docs_workflow = yaml.load(handle, Loader=yaml.BaseLoader)
    docs_body = next(
        step["run"]
        for step in docs_workflow["jobs"]["build"]["steps"]
        if step.get("id") == "docs"
    )
    jobs = workflow["jobs"]
    assert "runtime" in jobs["flake-checks"]["needs"]
    assert "runtime-guards" in jobs["flake-checks"]["needs"]
    assert "SYMBOLICA_LICENSE" not in json.dumps(jobs["runtime-guards"])
    select = next(
        step["run"]
        for step in jobs["prepare"]["steps"]
        if step.get("id") == "runtime-matrix"
    )
    guard = next(
        step["run"]
        for step in jobs["runtime-guards"]["steps"]
        if "run" in step and "ci-runtime-license-boundary" in step["run"]
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
        bash = shutil.which("bash")
        assert bash is not None, "Workflow fixture requires bash"
        fake_nix.write_text(
            f'#!{bash}\nprintf "%s\\n" "$@" > "$FIXTURE_ARGUMENTS"\n'
            'printf "%s\\n" "$FIXTURE_MATRIX"\n'
        )
        fake_nix.chmod(0o755)
        output = root / "output"
        env = dict(
            os.environ,
            PATH=str(root / "bin") + ":" + os.environ["PATH"],
            GITHUB_OUTPUT=str(output),
            FIXTURE_MATRIX=str(matrix),
            FIXTURE_ARGUMENTS=str(root / "nix-arguments"),
        )
        for key in ("SYMBOLICA_LICENSE", "SYMBOLICA_LICENSE_SIGNED"):
            env.pop(key, None)
        for ref, selection in (
            ("refs/heads/main", "ci-runtime-jobs-main"),
            ("refs/heads/manual", "ci-runtime-jobs"),
        ):
            output.write_text("")
            result = execute(select.replace("${{ github.ref }}", ref), env)
            assert result.returncode == 0, (
                f"Runtime selection failed: {result.stderr.decode()}"
            )
            assert (root / "nix-arguments").read_text().splitlines()[
                -1
            ] == ".#" + selection
            assert "any-jobs=true" in output.read_text(), (
                "Cached pure matrix suppressed runtime rows"
            )
            selected = json.loads(
                output.read_text().split("matrix=", 1)[1].splitlines()[0]
            )
            assert len(selected["include"]) == 1
        # Exercise the actual guard step with public source/export fixtures.
        # Reusing --flake-source must reach the probe and a failed probe must fail
        # the job, even when its preceding self-test derivation is cached.
        fake_python = root / "bin/python3"
        fake_python.write_text(
            f"#!{bash}\n"
            'test "$*" = "scripts/tracked_source.py --print-source" || exit 91\n'
            'printf "exported\\n" >> "$FIXTURE_EXPORTS"\n'
            'printf "%s\\n" "path:/nix/store/workflow-fixture-source"\n'
        )
        fake_python.chmod(0o755)
        fake_nix.write_text(
            f"#!{bash}\n"
            'printf "%s\\n" "$*" >> "$FIXTURE_GUARDS"\n'
            'if [ "$1" = run ]; then\n'
            '  test "$4" = --flake-source || exit 92\n'
            '  test "$5" = path:/nix/store/workflow-fixture-source || exit 93\n'
            '  exit "${FIXTURE_GUARD_STATUS:-0}"\n'
            "fi\n"
            'test "$1" = build || exit 94\n'
        )
        env.update(
            FIXTURE_EXPORTS=str(root / "exports"), FIXTURE_GUARDS=str(root / "guards")
        )
        assert execute(guard, env).returncode == 0
        assert (root / "exports").read_text().splitlines() == ["exported"]
        assert (root / "guards").read_text().splitlines() == [
            "build path:/nix/store/workflow-fixture-source#checks.x86_64-linux.gammaloop-ci-runtime --print-build-logs",
            "run path:/nix/store/workflow-fixture-source#ci-runtime-license-boundary -- --flake-source path:/nix/store/workflow-fixture-source --report-dir target/ci-results/license-boundary",
        ]
        assert execute(guard, dict(env, FIXTURE_GUARD_STATUS="37")).returncode == 37
        fake_python.unlink()
        verdict_env = dict(
            env,
            PREPARE_RESULT="success",
            GUARDS_RESULT="success",
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
        assert (
            execute(
                aggregate,
                dict(verdict_env, RUNTIME_RESULT="success", GUARDS_RESULT="failure"),
            ).returncode
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
        # Execute the current Pages branch with pure Nix fixtures. The real
        # step must strip both keys for eval/build and retain them for runtime.
        fake_git = root / "bin/git"
        fake_git.write_text(f"#!{bash}\nprintf '1\\n'\n")
        fake_git.chmod(0o755)
        fake_nix.write_text(
            f"#!{bash}\n"
            'test -z "${SYMBOLICA_LICENSE+x}" || exit 95\n'
            'test -z "${SYMBOLICA_LICENSE_SIGNED+x}" || exit 96\n'
            'printf "%s\\n" "$1" >> "$FIXTURE_DOCS_COMMANDS"\n'
            'case "$1" in\n'
            "  eval) printf 'yes\\n' ;;\n"
            '  build) test -x "$RUNNER_TEMP/alphal00p-docs-launcher/bin/alphal00p-docs-build" ;;\n'
            "  *) exit 97 ;;\n"
            "esac\n"
        )
        docs_launcher = root / "alphal00p-docs-launcher/bin/alphal00p-docs-build"
        docs_launcher.parent.mkdir(parents=True)
        docs_launcher.write_text(
            f"#!{bash}\n"
            'test "$SYMBOLICA_LICENSE" = public-docs-explicit || exit 98\n'
            'test "$SYMBOLICA_LICENSE_SIGNED" = public-docs-signed || exit 99\n'
            'mkdir -p "$RUNNER_TEMP/alphal00p-docs"\n'
        )
        docs_launcher.chmod(0o755)
        docs_env = dict(
            env,
            SYMBOLICA_LICENSE="public-docs-explicit",
            SYMBOLICA_LICENSE_SIGNED="public-docs-signed",
            CHANNEL="latest",
            SNAPSHOT_TAG="",
            FIXTURE_DOCS_COMMANDS=str(root / "docs-commands"),
        )
        result = execute(docs_body, docs_env)
        assert result.returncode == 0, result.stderr.decode()
        assert (root / "docs-commands").read_text().splitlines() == ["eval", "build"]
        assert (root / "alphal00p-docs").is_dir()
        assert "public-docs-" not in (result.stdout + result.stderr).decode()
    (args.report_dir / "receipt.json").write_text(
        json.dumps(
            {
                "empty_pure_matrix_keeps_runtime": True,
                "runtime_failure_fails_aggregate": True,
                "guard_failure_fails_aggregate": True,
                "guard_reuses_source": True,
                "guard_probe_failure_preserved": True,
                "main_and_manual_runtime_selections": True,
                "cached_launcher_executions": 2,
                "runtime_exit_37_preserved": True,
                "pages_credentials_only_at_runtime": True,
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
