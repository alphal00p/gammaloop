"""Promotion must verify and publish before changing the reviewed seed input."""

import io
import json
import os
import stat
import subprocess
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest.mock import Mock, patch

from test_dev import SOURCE, TemporaryTestCase, dev


class PromotionTests(TemporaryTestCase):
    def setUp(self):
        super().setUp()
        self.root = self.temporary / "workspace"
        self.root.mkdir()
        self.cache = self.temporary / "cache"
        self.netrc = self.temporary / "netrc"
        self.netrc.write_text("fixture credentials; never used for a real upload")
        self.manifest = self.root / "nix" / "development-seeds.json"
        self.environment = patch.dict(
            os.environ,
            {
                "GAMMALOOP_DEV_CACHE": str(self.cache),
                "NIXCI_NETRC": str(self.netrc),
            },
        )
        self.environment.start()
        self.addCleanup(self.environment.stop)
        self.key = "test-binaries-example-" + "a" * 64
        self.descriptor = {
            "owner": "test-binaries-example",
            "compatibility": {"opaque": [1, {"native": "context"}]},
            "source": SOURCE,
            "drvPath": self.store_path("actual.drv"),
            "unseededDrvPath": self.store_path("unseeded.drv"),
            "out": self.store_path("compiled"),
            "incremental": self.store_path("queries"),
            "seed": {
                "drvPath": self.store_path("seed.drv"),
                "incremental": self.store_path("seed-queries"),
            },
        }
        self.paths = [
            SOURCE,
            self.descriptor["drvPath"],
            self.descriptor["unseededDrvPath"],
            self.descriptor["out"],
            self.descriptor["incremental"],
            self.descriptor["seed"]["drvPath"],
            self.descriptor["seed"]["incremental"],
        ]
        self.api = Mock()
        self.api.prepare_ci.side_effect = lambda state: {
            "state": state,
            "paths": self.paths,
        }

    @staticmethod
    def store_path(name):
        return "/nix/store/" + "0" * 32 + "-" + name

    def seed(self, *, successful=True, contexts=None):
        registry = dev.CacheRegistry(self.cache)
        wrapper = self.make_result(
            "known-result", contexts=contexts or {self.key: self.descriptor}
        )
        with registry.generation() as generation:
            link = generation / "result-0"
            link.symlink_to(wrapper, target_is_directory=True)
            registry.publish(generation, SOURCE, [link], successful=successful)

    def invoke(self, *arguments, copy_error=None, source_changed=False):
        stdout, stderr = io.StringIO(), io.StringIO()
        sources = (
            [SOURCE, self.store_path("changed-source")]
            if source_changed
            else [SOURCE, SOURCE]
        )
        with (
            patch.object(dev.SourceSnapshot, "store", side_effect=sources),
            patch.object(dev, "NixDevelopment", return_value=self.api),
            patch.object(dev.subprocess, "run", side_effect=copy_error) as copy,
            redirect_stdout(stdout),
            redirect_stderr(stderr),
        ):
            status = dev.prepare_ci(
                dev.parse_args(["prepare-ci", *arguments]), self.root
            )
        return status, stdout.getvalue(), stderr.getvalue(), copy

    def test_full_actual_and_seed_roots_publish_before_atomic_manifest(self):
        self.seed()

        def copying(command, **kwargs):
            self.assertFalse(self.manifest.exists())
            self.assertEqual(command[:2], ["nix", "copy"])
            self.assertEqual(
                command[command.index("--netrc-file") + 1], str(self.netrc)
            )
            self.assertEqual(kwargs["input"].splitlines(), self.paths)
            self.assertTrue(kwargs["check"])
            self.assertNotIn("shell", kwargs)

        status, stdout, _, copy = self.invoke(copy_error=copying)
        self.assertEqual(status, 0)
        self.assertEqual(stdout.strip(), str(self.manifest))
        self.assertEqual(
            json.loads(self.manifest.read_text()),
            {"version": 1, "contexts": {self.key: self.descriptor}},
        )
        self.assertEqual(stat.S_IMODE(self.manifest.stat().st_mode), 0o644)
        self.assertEqual(
            self.api.prepare_ci.call_args.args[0]["contexts"][self.key], self.descriptor
        )
        copy.assert_called_once()

    def test_dry_run_needs_no_credentials_and_does_not_publish_or_write(self):
        self.seed()
        self.netrc.unlink()
        before = dev.CacheRegistry(self.cache).read()
        status, stdout, _, copy = self.invoke("--dry-run")
        self.assertEqual(status, 0)
        self.assertEqual(json.loads(stdout)["paths"], self.paths)
        self.assertFalse(self.manifest.exists())
        copy.assert_not_called()
        self.assertEqual(dev.CacheRegistry(self.cache).read(), before)

    def test_missing_and_failed_only_state_are_not_promotable(self):
        with self.assertRaisesRegex(dev.DevError, "no successful rooted"):
            self.invoke()
        self.seed(successful=False)
        with self.assertRaisesRegex(dev.DevError, "no successful rooted"):
            self.invoke()
        self.api.prepare_ci.assert_not_called()
        self.assertFalse(self.manifest.exists())

    def test_unknown_requested_context_fails_before_nix(self):
        self.seed()
        with self.assertRaisesRegex(
            dev.DevError, "no successful rooted state for contexts"
        ):
            self.invoke("--context", "unknown-context")
        self.api.prepare_ci.assert_not_called()

    def test_static_owners_are_reported_and_omitted(self):
        static = "check-narrow-" + "b" * 64
        self.seed(
            contexts={self.key: self.descriptor, static: {"owner": "check-narrow"}}
        )
        self.api.prepare_ci.side_effect = lambda state: {
            "state": {
                "version": 1,
                "contexts": {self.key: state["contexts"][self.key]},
            },
            "paths": self.paths,
        }
        _, _, stderr, _ = self.invoke()
        self.assertIn(f"Omitted {static}", stderr)
        self.assertEqual(
            set(json.loads(self.manifest.read_text())["contexts"]), {self.key}
        )

    def test_no_matching_ci_consumer_never_writes_a_manifest(self):
        self.seed(contexts={"build-only": {"owner": "build-only"}})
        self.api.prepare_ci.return_value = {
            "state": {"version": 1, "contexts": {}},
            "paths": [],
        }
        self.api.prepare_ci.side_effect = None
        with self.assertRaisesRegex(
            dev.DevError, "no identical selected NixCI consumers"
        ):
            self.invoke()
        self.assertFalse(self.manifest.exists())

    def test_upload_failure_keeps_existing_manifest_exactly(self):
        self.seed()
        self.manifest.parent.mkdir()
        self.manifest.write_text("existing reviewed seed data\n")
        with self.assertRaises(subprocess.CalledProcessError):
            self.invoke(copy_error=subprocess.CalledProcessError(1, ["nix", "copy"]))
        self.assertEqual(self.manifest.read_text(), "existing reviewed seed data\n")
        self.assertEqual(list(self.manifest.parent.iterdir()), [self.manifest])

    def test_source_edit_during_upload_prevents_manifest_write(self):
        self.seed()
        with self.assertRaisesRegex(
            dev.DevError, "workspace changed during publication"
        ):
            self.invoke(source_changed=True)
        self.assertFalse(self.manifest.exists())

    def test_nix_identity_rejection_never_publishes(self):
        self.seed()
        self.api.prepare_ci.side_effect = subprocess.CalledProcessError(
            1, ["nix", "eval"]
        )
        with self.assertRaises(subprocess.CalledProcessError):
            self.invoke()
        self.assertFalse(self.manifest.exists())

    def test_endpoint_cannot_replace_or_add_opaque_contexts(self):
        self.seed()
        for contexts in (
            {self.key: self.descriptor | {"source": "wrong"}},
            {"unexpected": self.descriptor},
        ):
            with self.subTest(contexts=contexts):
                self.api.prepare_ci.side_effect = None
                self.api.prepare_ci.return_value = {
                    "state": {"version": 1, "contexts": contexts},
                    "paths": self.paths,
                }
                with self.assertRaisesRegex(dev.DevError, "changed an opaque context"):
                    self.invoke()
        self.assertFalse(self.manifest.exists())

    def test_invalid_store_roots_cannot_be_injected_into_stdin(self):
        self.seed()
        for path in (
            "/tmp/not-a-store-path",
            self.paths[0] + "\nextra",
            "--arbitrary-option",
        ):
            with self.subTest(path=path):
                self.api.prepare_ci.side_effect = lambda state, path=path: {
                    "state": state,
                    "paths": [path],
                }
                with self.assertRaisesRegex(dev.DevError, "invalid closure roots"):
                    self.invoke()
        self.assertFalse(self.manifest.exists())

    def test_license_field_is_never_published(self):
        self.seed(
            contexts={
                self.key: self.descriptor | {"SYMBOLICA_LICENSE": "fixture-secret"}
            }
        )
        with self.assertRaisesRegex(dev.DevError, "license fields are forbidden"):
            self.invoke()
        self.assertFalse(self.manifest.exists())

    def test_prepare_transport_uses_fixed_apply_and_private_json_file(self):
        state = {"version": 1, "contexts": {self.key: self.descriptor}}
        api = dev.NixDevelopment(self.root, SOURCE, "x86_64-linux")
        files = []

        def evaluate(command, **kwargs):
            request = Path(kwargs["env"]["GAMMALOOP_DEV_REQUEST_FILE"])
            files.append(request)
            self.assertEqual(json.loads(request.read_text()), state)
            self.assertEqual(stat.S_IMODE(request.stat().st_mode), 0o600)
            self.assertEqual(
                command[command.index("--apply") + 1], dev.PREPARE_CI_APPLY
            )
            self.assertNotIn("shell", kwargs)
            return json.dumps({"state": state, "paths": self.paths})

        with patch.object(dev.subprocess, "check_output", side_effect=evaluate):
            self.assertEqual(api.prepare_ci(state)["state"], state)
        self.assertTrue(files)
        self.assertTrue(all(not path.exists() for path in files))
