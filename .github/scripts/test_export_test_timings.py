"""Check the metadata-only boundary of JUnit timing exports."""

import csv
import json
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


class ExportTestTimingsTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.report = self.directory / "junit.xml"
        self.output = self.directory / "metrics.csv"

    def export(self, xml):
        self.report.write_text(xml)
        return subprocess.run(
            [
                sys.executable,
                str(Path(__file__).with_name("export-test-timings.py")),
                "--check",
                "public-check",
                "--output",
                str(self.output),
                str(self.report),
            ],
            capture_output=True,
            text=True,
            check=False,
        )

    def test_exports_metadata_and_statuses_without_private_output(self):
        private = "PRIVATE_CAPTURED_PAYLOAD"
        self.output.write_text("stale CSV")
        result = self.export(f"""<testsuite>
          <testcase classname="public-binary" name="pass" time="1.230">
            <system-out>{private}</system-out></testcase>
          <testcase classname="public-binary" name="fail" time="2.5">
            <failure message="{private}">{private}</failure></testcase>
          <testcase classname="public-binary" name="error" time="3">
            <error>{private}</error><system-err>{private}</system-err></testcase>
          <testcase classname="public-binary" name="skip">
            <skipped message="{private}"/></testcase>
        </testsuite>""")
        self.assertEqual(result.returncode, 0, result.stderr)
        with self.output.open() as output:
            self.assertEqual(
                list(csv.reader(output)),
                [
                    ["check", "test_binary", "test_name", "status", "seconds"],
                    ["public-check", "public-binary", "pass", "passed", "1.230"],
                    ["public-check", "public-binary", "fail", "failed", "2.5"],
                    ["public-check", "public-binary", "error", "error", "3"],
                    ["public-check", "public-binary", "skip", "skipped", ""],
                ],
            )
        self.assertNotIn(private, self.output.read_text())
        self.assertNotIn(private, result.stdout + result.stderr)

    def test_invalid_reports_remove_stale_csv_without_partial_rows(self):
        reports = ["<testsuite>"] + [
            f'<testsuite><testcase name="valid" time="1"/>'
            f'<testcase name="invalid" time="{seconds}"/></testsuite>'
            for seconds in ["not-a-number", "NaN", "Infinity", "-1"]
        ]
        for report in reports:
            with self.subTest(report=report):
                self.output.write_text("stale CSV")
                result = self.export(report)
                self.assertEqual(result.returncode, 2)
                self.assertFalse(self.output.exists())


class NixWorkflowTests(unittest.TestCase):
    """Exercise the GitHub timing consumer with harmless Nix and runner stubs."""

    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory()))
        self.bin = self.root / "bin"
        self.bin.mkdir()
        self.prepared = self.root / "prepared"
        (self.prepared / "bin").mkdir(parents=True)
        scripts = self.root / ".github/scripts"
        scripts.mkdir(parents=True)
        for name in ("run-nix-ci-check.sh", "export-test-timings.py"):
            shutil.copyfile(ROOT / ".github/scripts" / name, scripts / name)
        for name in ("env", "mkdir", "rm", "tee", "find", "jq"):
            (self.bin / name).symlink_to(shutil.which(name))
        (self.bin / "python3").symlink_to(sys.executable)
        self.env = {
            "PATH": str(self.bin),
            "RUNNER_TEMP": str(self.root),
            "GITHUB_OUTPUT": str(self.root / "github-output"),
            "GITHUB_SHA": "test-commit",
            "GITHUB_SERVER_URL": "https://example.invalid",
            "GITHUB_REPOSITORY": "test/repository",
            "GITHUB_RUN_ID": "1",
            "GITHUB_RUN_ATTEMPT": "1",
            "TEST_PREPARED": str(self.prepared),
            "TEST_ROOT": str(self.root),
            "SYMBOLICA_LICENSE": "runtime-sentinel",
            "SYMBOLICA_LICENSE_PATH": "unused-sentinel",
            "SYMBOLICA_LICENSE_SIGNED": "signed-sentinel",
            "SYMBOLICA_LICENSE_FILE": "unused-sentinel",
        }
        self.stub(
            self.bin / "nix",
            "assert not any(name in os.environ for name in "
            "('LICENSE', 'SYMBOLICA_LICENSE', 'SYMBOLICA_LICENSE_PATH', "
            "'SYMBOLICA_LICENSE_SIGNED', 'SYMBOLICA_LICENSE_FILE'))\n"
            "assert '--impure' in sys.argv\n"
            "(root / 'nix-call.json').write_text(json.dumps(sys.argv[1:]))\n"
            "if sys.argv[1] == 'eval':\n"
            "    print(os.environ['TEST_CHECK_NAMES'])\n"
            "    sys.exit(0)\n"
            "assert sys.argv[1] == 'build'\n"
            "status = int(os.environ.get('TEST_BUILD_EXIT', '0'))\n"
            "if status: sys.exit(status)\n"
            "link = Path(sys.argv[sys.argv.index('--out-link') + 1])\n"
            "link.unlink(missing_ok=True)\n"
            "link.symlink_to(os.environ['TEST_PREPARED'])\n"
            "print(os.environ['TEST_PREPARED'])\n",
        )

    def stub(self, path, body):
        path.write_text(
            f"#!{sys.executable}\nimport json, os, sys\nfrom pathlib import Path\n"
            "root = Path(os.environ['TEST_ROOT'])\n" + body
        )
        path.chmod(0o755)

    def script(self, step, workflow="nix.yml"):
        workflow = (ROOT / ".github/workflows" / workflow).read_text()
        block = workflow.split(f"      - id: {step}\n", 1)[1].split("\n      - ", 1)[0]
        return textwrap.dedent(block.split("        run: |\n", 1)[1])

    def invoke(self, script, **environment):
        return subprocess.run(
            [shutil.which("bash"), "--noprofile", "--norc", "-c", script],
            cwd=self.root,
            env=self.env | environment,
            capture_output=True,
            text=True,
            check=False,
            timeout=10,
        )

    def runner(self, name):
        self.stub(
            self.prepared / "bin" / f"nix-ci-check-{name}",
            "assert os.environ['SYMBOLICA_LICENSE'] == 'runtime-sentinel'\n"
            "assert len(sys.argv) == 2\n"
            "reports = Path(sys.argv[1])\n"
            "assert reports == root / 'target/ci-results/reports'\n"
            "with (root / 'runs.jsonl').open('a') as log:\n"
            "    log.write(json.dumps(sys.argv) + '\\n')\n"
            "if not os.environ.get('TEST_NO_REPORTS'):\n"
            "    for name in ('package.xml', 'browser/context with spaces/client.xml'):\n"
            "        report = reports / name\n"
            "        report.parent.mkdir(parents=True, exist_ok=True)\n"
            "        report.write_text(os.environ.get('TEST_XML', "
            '\'<testsuite><testcase name="test" time="1"/></testsuite>\'))\n'
            "sys.exit(int(os.environ.get('TEST_RUN_EXIT', '0')))\n",
        )

    def consume(self, attr, **environment):
        return self.invoke(
            ". .github/scripts/run-nix-ci-check.sh", CHECK_ATTR=attr, **environment
        )

    def test_cached_preparations_execute_exact_runner_and_export_nested_reports(self):
        rows = [
            f"{collection}.x86_64-linux.{name}"
            for collection in ("hydraJobs", "hydraJobsMain")
            for name in (
                "alphal00p-docs",
                "gammaloop-doctest",
                "gammaloop-nextest-core",
                "gammaloop-nextest-integration",
            )
        ] + [
            f"hydraJobsMain.x86_64-linux.gammaloop-nextest-long-{index}"
            for index in range(1, 5)
        ]
        for attr in rows:
            name = attr.rsplit(".", 1)[1]
            self.runner(name)
            with self.subTest(attr=attr):
                for _ in range(2):
                    result = self.consume(attr)
                    self.assertEqual(result.returncode, 0, result.stderr)
                runs = (self.root / "runs.jsonl").read_text().splitlines()
                self.assertEqual(len(runs), 2)
                self.assertTrue(
                    all(
                        json.loads(run)[0].endswith(f"nix-ci-check-{name}")
                        for run in runs
                    )
                )
                (self.root / "runs.jsonl").unlink()
                with (self.root / "target/ci-results/metrics.csv").open() as output:
                    records = list(csv.DictReader(output))
                self.assertEqual(len(records), 2)
                self.assertEqual({row["check"] for row in records}, {attr})
                call = json.loads((self.root / "nix-call.json").read_text())
                self.assertEqual(call[-1], f".#{attr}")

    def test_build_only_rows_do_not_execute_available_binaries(self):
        for name in (
            "gammaloop-doc",
            "alphal00p-docs-watch",
            "gammaloop-nextest-binaries",
            "gammaloop-nextest-binaries-core",
        ):
            with self.subTest(name=name):
                self.runner(name)
                result = self.consume(f"hydraJobs.x86_64-linux.{name}")
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertFalse((self.root / "runs.jsonl").exists())
                self.assertFalse((self.root / "target/ci-results/metrics.csv").exists())

    def test_runtime_failure_preserves_exit_and_exports_available_reports(self):
        self.runner("gammaloop-doctest")
        for run_exit, xml, expected in (
            (
                23,
                '<testsuite><testcase name="failed"><failure/></testcase></testsuite>',
                23,
            ),
            (23, "<malformed>", 23),
            (0, "<malformed>", 2),
        ):
            with self.subTest(run_exit=run_exit, xml=xml):
                result = self.consume(
                    "hydraJobs.x86_64-linux.gammaloop-doctest",
                    TEST_RUN_EXIT=str(run_exit),
                    TEST_XML=xml,
                )
                self.assertEqual(result.returncode, expected, result.stderr)
                provenance = (
                    self.root / "target/ci-results/provenance.txt"
                ).read_text()
                self.assertIn(f"run-exit-code={run_exit}", provenance)
                metrics = self.root / "target/ci-results/metrics.csv"
                self.assertEqual(metrics.exists(), xml != "<malformed>")
                if metrics.exists():
                    self.assertIn(",failed,", metrics.read_text())

    def test_failed_preparation_and_missing_runner_cannot_pass(self):
        attr = "hydraJobs.x86_64-linux.gammaloop-doctest"
        result = self.consume(attr)
        self.assertEqual(result.returncode, 127, result.stderr)
        self.runner("gammaloop-doctest")
        result = self.consume(attr, TEST_BUILD_EXIT="17")
        self.assertEqual(result.returncode, 17, result.stderr)
        self.assertFalse((self.root / "runs.jsonl").exists())

    def test_previous_runtime_reports_are_not_reused(self):
        self.runner("alphal00p-docs")
        attr = "hydraJobs.x86_64-linux.alphal00p-docs"
        self.assertEqual(self.consume(attr).returncode, 0)
        result = self.consume(attr, TEST_NO_REPORTS="1")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse((self.root / "target/ci-results/metrics.csv").exists())

    def test_matrix_restores_cached_runtime_rows_once_and_filters_all_licensed_rows(
        self,
    ):
        runtime = ["alphal00p-docs", "gammaloop-doctest", "gammaloop-nextest-core"]
        for collection in ("hydraJobs", "hydraJobsMain"):
            names = runtime + (
                [f"gammaloop-nextest-long-{index}" for index in range(1, 5)]
                if collection == "hydraJobsMain"
                else []
            )
            collection += ".x86_64-linux"
            static = [
                {"attr": f"{collection}.{name}", "name": name}
                for name in ("gammaloop-doc", "gammaloop-nextest-binaries-core")
            ]
            for licensed in (False, True):
                for cached in (False, True):
                    with self.subTest(
                        collection=collection, licensed=licensed, cached=cached
                    ):
                        existing = static.copy() if not cached else []
                        if not cached:
                            existing += [
                                {"attr": f"{collection}.{name}"} for name in names
                            ]
                        result = self.invoke(
                            self.script("filter-matrix"),
                            CHECK_COLLECTION=collection,
                            MATRIX=json.dumps({"include": existing}),
                            LICENSE="license-sentinel" if licensed else "",
                            TEST_CHECK_NAMES=json.dumps(
                                names
                                + ["gammaloop-doc", "gammaloop-nextest-binaries-core"]
                            ),
                        )
                        self.assertEqual(result.returncode, 0, result.stderr)
                        outputs = dict(
                            line.split("=", 1)
                            for line in Path(self.env["GITHUB_OUTPUT"])
                            .read_text()
                            .splitlines()
                        )
                        matrix = json.loads(outputs["matrix"])["include"]
                        expected = (
                            [f"{collection}.{name}" for name in names]
                            if licensed
                            else []
                        )
                        if not cached:
                            expected = [row["attr"] for row in static] + expected
                        self.assertEqual([row["attr"] for row in matrix], expected)
                        self.assertEqual(
                            outputs["any-jobs"], str(bool(expected)).lower()
                        )
                        for row in matrix[len(static) if not cached else 0 :]:
                            self.assertEqual(row["system"], "x86_64-linux")
                            self.assertEqual(row["os"], ["ubuntu-24.04"])

    def test_matrix_uses_selected_collection_and_secret_free_impure_evaluator(self):
        self.stub(
            self.bin / "hestia",
            "assert not any(name.startswith('SYMBOLICA_LICENSE') for name in os.environ)\n"
            "(root / 'hestia-call.json').write_text(json.dumps(sys.argv[1:]))\n",
        )
        for ref, collection in (
            ("refs/heads/main", "hydraJobsMain.x86_64-linux"),
            ("refs/heads/feature", "hydraJobs.x86_64-linux"),
        ):
            with self.subTest(ref=ref):
                script = self.script("matrix").replace("${{ github.ref }}", ref)
                result = self.invoke(
                    script,
                    HESTIA_BIN=str(self.bin / "hestia"),
                    HESTIA_SOCKET="unused-socket",
                    HESTIA_DRAIN_TIMEOUT="1",
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                call = json.loads((self.root / "hestia-call.json").read_text())
                self.assertEqual(call[call.index("--flake") + 1], f".#{collection}")
                self.assertEqual(call[call.index("--attr-prefix") + 1], collection)
                self.assertIn("-- --impure", call[call.index("--nix-eval-jobs") + 1])
                self.assertIn(
                    f"check-collection={collection}",
                    Path(self.env["GITHUB_OUTPUT"]).read_text().splitlines(),
                )

    def prepare_docs(self):
        (self.bin / "grep").symlink_to(shutil.which("grep"))
        self.stub(
            self.bin / "git",
            "print('test-commit' if 'rev-parse' in sys.argv else '1')\n",
        )
        self.stub(
            self.prepared / "bin/alphal00p-docs-pages",
            "assert os.environ['SYMBOLICA_LICENSE'] == 'runtime-sentinel'\n"
            "assert len(sys.argv) == 2\n"
            "destination = Path(sys.argv[1])\n"
            "assert destination == Path(os.environ['RUNNER_TEMP']) / 'alphal00p-docs'\n"
            "destination.mkdir()\n"
            "(destination / 'index.html').write_text('rendered fixture')\n"
            "sys.exit(int(os.environ.get('TEST_RUN_EXIT', '0')))\n",
        )

    def test_docs_pages_executes_prepared_latest_and_snapshot_renderers(self):
        self.prepare_docs()
        for channel, tag in (("latest", ""), ("snapshot", "v0.3.4")):
            with self.subTest(channel=channel):
                temporary = self.root / channel
                temporary.mkdir()
                result = self.invoke(
                    self.script("docs", "docs-pages.yml"),
                    CHANNEL=channel,
                    SNAPSHOT_TAG=tag,
                    GITHUB_EVENT_NAME="workflow_dispatch",
                    SNAPSHOT_RUNNER="true",
                    RUNNER_TEMP=str(temporary),
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(
                    (temporary / "alphal00p-docs/index.html").read_text(),
                    "rendered fixture",
                )
                self.assertIn(
                    f"path={temporary}/alphal00p-docs",
                    Path(self.env["GITHUB_OUTPUT"]).read_text().splitlines(),
                )

    def test_docs_pages_preserves_renderer_failure(self):
        self.prepare_docs()
        result = self.invoke(
            self.script("docs", "docs-pages.yml"),
            CHANNEL="latest",
            SNAPSHOT_TAG="",
            GITHUB_EVENT_NAME="push",
            TEST_RUN_EXIT="27",
        )
        self.assertEqual(result.returncode, 27, result.stderr)
        self.assertFalse(Path(self.env["GITHUB_OUTPUT"]).exists())


if __name__ == "__main__":
    unittest.main()
