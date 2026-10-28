"""Provisioned shared policy has authority over per-command retention flags."""

import grp
import json
import os
from unittest.mock import patch

from test_dev import TemporaryTestCase, dev


class SharedPolicyTests(TemporaryTestCase):
    def setUp(self):
        super().setUp()
        self.shared = self.temporary / "shared"
        self.shared.mkdir()
        self.shared.chmod(0o2770)
        self.policy_file = self.temporary / "policy.json"
        self.policy = {
            "version": 1,
            "directory": str(self.shared),
            "group": grp.getgrgid(self.shared.stat().st_gid).gr_name,
            "keepSuccessfulGenerations": 6,
            "maxGenerations": 73,
            "maxAgeDays": 11,
            "preserveActiveLeases": True,
        }
        self.write_policy()
        environment = patch.dict(
            os.environ,
            {
                "GAMMALOOP_DEV_CACHE": str(self.shared),
                "GAMMALOOP_DEV_CACHE_POLICY": str(self.policy_file),
            },
            clear=True,
        )
        environment.start()
        self.addCleanup(environment.stop)

    def write_policy(self):
        self.policy_file.write_text(json.dumps(self.policy))
        self.policy_file.chmod(0o644)

    def registry(self, *arguments):
        return dev.CacheRegistry.from_environment(dev.parse_args(["test", *arguments]))

    def test_policy_limits_replace_private_defaults_on_matching_provisioned_root(self):
        registry = self.registry()
        self.assertEqual(registry.directory, self.shared)
        self.assertTrue(registry.shared)
        self.assertEqual(registry.keep, 6)
        self.assertEqual(registry.maximum, 73)
        self.assertEqual(registry.max_age, 11 * 86400)
        self.assertEqual(self.registry("--keep-generations", "6").keep, 6)

    def test_each_conflicting_shared_retention_override_is_rejected(self):
        for flag in ("--keep-generations", "--max-generations", "--max-age-days"):
            with (
                self.subTest(flag=flag),
                self.assertRaisesRegex(dev.DevError, "per-command override rejected"),
            ):
                self.registry(flag, "2")

    def test_shared_target_must_match_policy_and_group(self):
        with (
            patch.dict(
                os.environ, {"GAMMALOOP_DEV_CACHE": str(self.temporary / "elsewhere")}
            ),
            self.assertRaisesRegex(dev.DevError, "differs from shared policy"),
        ):
            self.registry()
        self.policy["group"] = "gammaloop-nonexistent-fixture-group"
        self.write_policy()
        with self.assertRaisesRegex(dev.DevError, "group does not exist"):
            self.registry()

    def test_shared_root_must_already_have_exact_provisioned_mode(self):
        for mode in (0o770, 0o2775):
            self.shared.chmod(mode)
            with (
                self.subTest(mode=mode),
                self.assertRaisesRegex(dev.DevError, "mode 2770"),
            ):
                self.registry()
            self.assertEqual(self.shared.stat().st_mode & 0o7777, mode)
        self.assertFalse((self.shared / "runs").exists())

    def test_policy_root_and_version_fields_are_strict(self):
        original = self.policy.copy()
        invalid = (
            {"version": 2},
            {"version": True},
            {"preserveActiveLeases": False},
            {"keepSuccessfulGenerations": 0},
            {"keepSuccessfulGenerations": True},
            {"maxGenerations": -1},
            {"maxAgeDays": float("nan")},
            {"directory": "relative"},
            {"group": ""},
            {"unrecognized": True},
        )
        for changes in invalid:
            with self.subTest(changes=changes):
                self.policy = original | changes
                self.write_policy()
                with self.assertRaises(dev.DevError):
                    self.registry()
        self.policy = original.copy()
        del self.policy["maxGenerations"]
        self.write_policy()
        with self.assertRaises(dev.DevError):
            self.registry()

    def test_policy_cannot_be_group_or_world_writable(self):
        self.policy_file.chmod(0o666)
        with self.assertRaisesRegex(dev.DevError, "not writable by group/others"):
            self.registry()
        self.assertFalse((self.shared / "runs").exists())

    def test_policy_directory_cannot_redirect_through_a_symlink(self):
        alias = self.temporary / "alias"
        alias.symlink_to(self.shared, target_is_directory=True)
        self.policy["directory"] = str(alias)
        self.write_policy()
        with self.assertRaisesRegex(dev.DevError, "real provisioned directory"):
            self.registry()

    def test_private_override_is_explicit_and_does_not_read_shared_policy(self):
        self.policy_file.write_text(
            "invalid policy deliberately bypassed for private storage"
        )
        private = self.temporary / "private"
        registry = self.registry(
            "--private-cache",
            str(private),
            "--keep-generations",
            "2",
            "--max-generations",
            "3",
        )
        self.assertFalse(registry.shared)
        self.assertEqual(registry.directory, private)
        self.assertEqual((registry.keep, registry.maximum), (2, 3))
        self.assertEqual(private.stat().st_mode & 0o7777, 0o700)
        self.assertFalse((self.shared / "runs").exists())

    def test_private_override_cannot_use_shared_or_relaxed_storage(self):
        with self.assertRaisesRegex(dev.DevError, "private storage must be private"):
            self.registry("--private-cache", str(self.shared))
        private = self.temporary / "relaxed-private"
        private.mkdir(mode=0o755)
        with self.assertRaisesRegex(dev.DevError, "mode 0700"):
            self.registry("--private-cache", str(private))
        self.assertEqual(private.stat().st_mode & 0o7777, 0o755)

    def test_shared_env_directory_requires_policy(self):
        with (
            patch.dict(
                os.environ, {"GAMMALOOP_DEV_CACHE": str(self.shared)}, clear=True
            ),
            self.assertRaisesRegex(
                dev.DevError, "requires a matching GAMMALOOP_DEV_CACHE_POLICY"
            ),
        ):
            self.registry()

    def test_unconfigured_xdg_default_keeps_final_private_suffix(self):
        xdg = self.temporary / "xdg"
        with patch.dict(os.environ, {"XDG_CACHE_HOME": str(xdg)}, clear=True):
            registry = self.registry()
        self.assertEqual(registry.directory, xdg / "gammaloop" / "development")
        self.assertEqual(
            (registry.keep, registry.maximum, registry.max_age), (4, 64, 30 * 86400)
        )
        self.assertFalse(registry.shared)
