"""Read-only publication coverage contracts with a mocked local Nix store."""

import contextlib
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import cache_publication as publication


class PublicationTests(unittest.TestCase):
    def setUp(self):
        self.source = Path(self.enterContext(tempfile.TemporaryDirectory()))
        (self.source / "flake.nix").touch()
        self.enterContext(
            mock.patch.dict(os.environ, {"GAMMALOOP_FLAKE_SOURCE": str(self.source)})
        )

    def inspect(self, missing=(), source_dependent=False):
        metadata = {
            "sourceRevisionDependent": ["docs"] if source_dependent else [],
            "local": ["producer", "launcher"] + (["docs"] if source_dependent else []),
            "github": ["launcher"],
            "producers": ["producer"],
            "githubMaxNarBytes": 8589934592,
        }
        outputs = {
            "producer": [{"output": "out", "path": "/nix/store/producer"}],
            "launcher": [{"output": "out", "path": "/nix/store/launcher"}],
        }
        if source_dependent:
            outputs["docs"] = [{"output": "out", "path": "/nix/store/docs"}]
        calls = []

        def run(command, **kwargs):
            calls.append(command)
            self.assertNotIn("SYMBOLICA_LICENSE", kwargs["env"])
            self.assertNotIn("SYMBOLICA_LICENSE_SIGNED", kwargs["env"])
            self.assertNotIn("NIX_CONFIG", kwargs["env"])
            self.assertNotIn("NIXCI_NETRC", kwargs["env"])
            self.assertNotIn("NIX_USER_CONF_FILES", kwargs["env"])
            self.assertNotIn("GAMMALOOP_NIX_TEST_ARGS", kwargs["env"])
            self.assertNotIn("GAMMALOOP_DEVELOP_SPEC", kwargs["env"])
            self.assertNotIn("NIX_REMOTE", kwargs["env"])
            self.assertFalse(set(command) & {"build", "copy", "--realise", "run"})
            if "eval" in command:
                self.assertIn("allow-import-from-derivation", command)
                self.assertIn("false", command)
                self.assertIn("--offline", command)
                index = command.index("--extra-experimental-features")
                self.assertEqual(command[index + 1], "nix-command flakes")
                result = json.dumps({"metadata": metadata, "outputs": outputs})
            elif "--check-validity" in command:
                return subprocess.CompletedProcess(
                    command, int(command[-1] in missing), "", ""
                )
            else:
                self.assertIn("path-info", command)
                self.assertIn("--offline", command)
                # The same requisite occurs under both roots, once in JSON.
                result = json.dumps(
                    {
                        path: {"narSize": size}
                        for path, size in [
                            ("/nix/store/launcher", 10),
                            ("/nix/store/shared", 20),
                        ]
                    }
                )
            return subprocess.CompletedProcess(command, 0, result, "")

        with (
            mock.patch.dict(
                os.environ,
                {
                    "SYMBOLICA_LICENSE": "sentinel",
                    "SYMBOLICA_LICENSE_SIGNED": "sentinel",
                    "NIX_CONFIG": "secret",
                    "NIXCI_NETRC": "secret",
                    "NIX_USER_CONF_FILES": "secret",
                    "NIX_REMOTE": "ssh://never",
                    "GAMMALOOP_NIX_TEST_ARGS": "stale",
                    "GAMMALOOP_DEVELOP_SPEC": "stale",
                },
            ),
            mock.patch.object(publication.subprocess, "run", side_effect=run),
            contextlib.redirect_stdout(io.StringIO()) as output,
        ):
            publication.preview(Path.cwd(), ["both"])
        return output.getvalue(), calls

    def test_complete_deduplicated_closure_and_dropped_intermediate(self):
        output, _ = self.inspect()
        self.assertIn("2 deduplicated closure paths; 30 raw NAR bytes", output)
        self.assertIn("confirmed absent from full set", output)
        self.assertIn("producer.out /nix/store/producer (realized)", output)
        self.assertNotIn("sentinel", output)

    def test_missing_roots_are_partial_and_full_absence_unknown(self):
        output, calls = self.inspect(["/nix/store/launcher"])
        self.assertIn("PARTIAL coverage", output)
        self.assertIn("full set unknown", output)
        self.assertIn("launcher.out /nix/store/launcher missing", output)
        for call in calls:
            if "path-info" in call:
                self.assertNotIn("/nix/store/launcher", call)

    def test_all_missing_performs_no_closure_query(self):
        output, calls = self.inspect(["/nix/store/launcher", "/nix/store/producer"])
        self.assertIn("0 deduplicated closure paths; 0 raw NAR bytes", output)
        self.assertFalse(any("path-info" in call for call in calls))

    def test_source_dependent_root_excluded_from_coverage_and_size(self):
        for missing in ((), ("/nix/store/docs",)):
            with self.subTest(missing=missing):
                output, calls = self.inspect(missing, source_dependent=True)
                self.assertIn("local: 2/2 comparable root outputs realized", output)
                self.assertNotIn("PARTIAL", output)
                self.assertIn("complete comparable coverage", output)
                self.assertIn(
                    "docs.out /nix/store/docs not comparable across source identities",
                    output,
                )
                self.assertIn("Entire publication size unknown", output)
                self.assertIn(
                    "Comparable realized-root subtotal: 2 deduplicated closure paths; 30 raw NAR bytes",
                    output,
                )
                self.assertTrue(all("/nix/store/docs" not in call for call in calls))
                self.assertIn("confirmed absent from full set", output)

    def test_malformed_scope_never_starts_queries(self):
        with mock.patch.object(publication.subprocess, "run") as run:
            with (
                contextlib.redirect_stderr(io.StringIO()),
                self.assertRaises(SystemExit),
            ):
                publication.preview(Path.cwd(), ["garbage"])
            run.assert_not_called()

    def test_missing_or_invalid_source_never_queries(self):
        for source in ("", str(self.source / "missing")):
            with (
                mock.patch.dict(os.environ, {"GAMMALOOP_FLAKE_SOURCE": source}),
                mock.patch.object(publication.subprocess, "run") as run,
            ):
                with self.assertRaisesRegex(ValueError, "GAMMALOOP_FLAKE_SOURCE"):
                    publication.preview(Path.cwd(), [])
                run.assert_not_called()

    def test_default_and_explicit_scopes(self):
        with mock.patch.object(publication.PublicationPreview, "inspect") as inspect:
            for arguments, scope in (
                ([], "both"),
                (["local"], "local"),
                (["github"], "github"),
            ):
                publication.preview(Path.cwd(), arguments)
                inspect.assert_called_with(scope)

    def test_evaluation_failure_does_not_fall_back_to_build(self):
        with (
            mock.patch.object(
                publication.subprocess,
                "run",
                return_value=subprocess.CompletedProcess([], 1, "", "IFD disabled"),
            ) as run,
        ):
            with self.assertRaisesRegex(ValueError, "IFD disabled"):
                publication.preview(Path.cwd(), [])
            self.assertEqual(run.call_count, 1)


class UploadSelectionTests(unittest.TestCase):
    def test_exact_flake_schedule_drives_selection_and_builds(self):
        repository = Path(__file__).resolve().parents[2]
        recipe = (
            (repository / "just/ci.just")
            .read_text()
            .split("ci-checks-and-upload:\n", 1)[1]
            .split("\n# Check the workspace graph", 1)[0]
        )
        script = "\n".join(line.removeprefix("    ") for line in recipe.splitlines())
        for caller_jobs, selected_jobs in [
            (["old"], ["old", "new"]),
            (["old", "unselected"], ["old"]),
        ]:
            with self.subTest(caller=caller_jobs, selected=selected_jobs):
                with tempfile.TemporaryDirectory() as temporary:
                    root = Path(temporary)
                    caller, selected = root / "caller", root / "selected"
                    caller.mkdir()
                    selected.mkdir()
                    jobs = lambda names: [f"packages.x86_64-linux.{n}" for n in names]
                    for directory, names in [
                        (caller, caller_jobs),
                        (selected, selected_jobs),
                    ]:
                        values = " ".join(json.dumps(name) for name in jobs(names))
                        (directory / "flake.nix").write_text(
                            "{ outputs = { self }: { lib.x86_64-linux.cachePublication.local = [ "
                            + values
                            + " ]; }; }\n"
                        )
                    program = root / "upload.sh"
                    program.write_text(script)
                    binary = root / "bin"
                    binary.mkdir()
                    log = root / "calls.jsonl"
                    stub = binary / "nix"
                    stub.write_text(
                        f"#!{sys.executable}\n"
                        "import json, os, sys\n"
                        "from pathlib import Path\n"
                        "a = sys.argv[1:]\n"
                        "record = {'args': a}\n"
                        "if a[0] == 'copy': record['stdin'] = sys.stdin.read()\n"
                        "with open(os.environ['UPLOAD_TEST_LOG'], 'a') as f: f.write(json.dumps(record) + '\\n')\n"
                        "if a[0] == 'eval':\n"
                        "    if '--expr' in a and a[a.index('--expr')+1] == 'builtins.currentSystem': print('x86_64-linux')\n"
                        "    elif os.environ['UPLOAD_TEST_REF'] + '#lib.x86_64-linux.cachePublication.local' in a: print('\\n'.join(json.loads(os.environ['UPLOAD_TEST_SELECTED'])))\n"
                        "    elif '--expr' in a: print('\\n'.join(json.loads(os.environ['UPLOAD_TEST_CALLER'])))\n"
                        "    else: raise AssertionError(a)\n"
                        "elif a[0] == 'build':\n"
                        "    for i, target in enumerate(x for x in a if '#' in x): print('/nix/store/upload-fixture-' + str(i))\n"
                        "elif a[0] != 'copy': raise AssertionError(a)\n"
                    )
                    stub.chmod(0o755)
                    netrc = root / "netrc"
                    netrc.touch()
                    reference = f"path:{selected}"
                    environment = os.environ | {
                        "PATH": f"{binary}:{os.environ['PATH']}",
                        "NIXCI_NETRC": str(netrc),
                        "GAMMALOOP_CI_FLAKE_REF": reference,
                        "UPLOAD_TEST_REF": reference,
                        "UPLOAD_TEST_LOG": str(log),
                        "UPLOAD_TEST_CALLER": json.dumps(jobs(caller_jobs)),
                        "UPLOAD_TEST_SELECTED": json.dumps(jobs(selected_jobs)),
                    }
                    result = subprocess.run(
                        [shutil.which("bash"), str(program)],
                        cwd=caller,
                        env=environment,
                        text=True,
                        capture_output=True,
                        check=False,
                    )
                    self.assertEqual(result.returncode, 0, result.stderr)
                    calls = [json.loads(line) for line in log.read_text().splitlines()]
                    build = next(
                        row["args"] for row in calls if row["args"][0] == "build"
                    )
                    self.assertEqual(
                        [arg for arg in build if "#" in arg],
                        [f"{reference}#{name}" for name in jobs(selected_jobs)],
                    )
                    copy = next(row for row in calls if row["args"][0] == "copy")
                    self.assertEqual(
                        len(copy["stdin"].splitlines()), len(selected_jobs)
                    )
                    self.assertIn("CI upload completed", result.stdout)


if __name__ == "__main__":
    unittest.main()
