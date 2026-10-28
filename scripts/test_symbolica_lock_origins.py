"""Small declaration/lock controls; no Cargo, Nix, native kernel or network calls."""

import copy
import importlib.util
import json
import shlex
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).with_name("check-symbolica-lock-origins.py")
SPEC = importlib.util.spec_from_file_location("lock_origins", SCRIPT)
origins = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(origins)
REV = "1" * 40
OTHER = "2" * 40
URL = "https://github.com/symbolica-dev/symbolica"


def documents(revision=REV, url=URL):
    manifest = {
        "workspace": {
            "dependencies": {"symbolica": {"version": "3.0.0"}, "numerica": "^3.0.1"}
        },
        "patch": {
            "crates-io": {
                name: {"git": url, "rev": revision} for name in origins.FAMILY
            }
        },
    }
    lock = {
        "package": [
            {
                "name": name,
                "version": "3.0.1",
                "source": f"git+{url}?rev={revision}#{revision}",
            }
            for name in origins.FAMILY
        ]
    }
    return manifest, lock


class OriginTests(unittest.TestCase):
    def test_each_workspace_uses_its_own_coherent_revision(self):
        root = origins.validate(*documents())
        tydenso = origins.validate(
            *documents(OTHER, "https://github.com/ValentinHirschi/symbolica")
        )
        self.assertNotEqual(root[1], tydenso[1])

    def test_registry_origin_and_cargo_caret_versions(self):
        manifest, lock = documents()
        manifest.pop("patch")
        for entry in lock["package"]:
            entry["source"] = origins.CRATES_IO
        origins.validate(manifest, lock)
        for requirement, version, accepted in [
            ("3.0.0", "3.0.1", True),
            ("^3.0.0", "4.0.0", False),
            ("=3.0.0", "3.0.1", False),
            ("0.3.0", "0.4.0", False),
            ("0.0.3", "0.0.4", False),
        ]:
            with self.subTest(requirement=requirement):
                self.assertEqual(origins.permits(requirement, version), accepted)
        with self.assertRaises(origins.InvalidFamily):
            origins.permits(">=3.0.0", "3.0.1")
        for malformed in ("03.0.0", "3.00.0", "٣.0.0"):
            with (
                self.subTest(malformed=malformed),
                self.assertRaises(origins.InvalidFamily),
            ):
                origins.permits(malformed, "3.0.1")

    def test_duplicate_or_missing_family_member_rejected(self):
        for mutate in (
            lambda entries: entries.append(copy.deepcopy(entries[0])),
            lambda entries: entries.pop(),
        ):
            manifest, lock = documents()
            mutate(lock["package"])
            with self.assertRaises(origins.InvalidFamily):
                origins.validate(manifest, lock)

    def test_wrong_source_selector_or_resolved_commit_rejected(self):
        wrong = [
            f"git+{URL}?rev={OTHER}#{OTHER}",
            f"git+{URL}?rev={REV}#{OTHER}",
            f"git+{URL}?rev={REV}",
            f"git+{URL}?branch=main#{REV}",
            f"git+{URL}?tag=v3.0.1#{REV}",
            f"git+{URL}?rev={REV}&rev={REV}#{REV}",
            f"git+https://github.com/other/symbolica?rev={REV}#{REV}",
            origins.CRATES_IO,
        ]
        for source in wrong:
            with self.subTest(source=source):
                manifest, lock = documents()
                lock["package"][1]["source"] = source
                with self.assertRaises(origins.InvalidFamily):
                    origins.validate(manifest, lock)

    def test_mixed_version_patch_or_direct_requirement_rejected(self):
        manifest, lock = documents()
        lock["package"][1]["version"] = "3.0.2"
        with self.assertRaises(origins.InvalidFamily):
            origins.validate(manifest, lock)
        for change in [
            "partial",
            "mixed",
            "wrong-table",
            "branch",
            "direct-git",
            "incompatible-numerica",
        ]:
            with self.subTest(change=change):
                manifest, lock = documents()
                if change == "partial":
                    del manifest["patch"]["crates-io"]["graphica"]
                elif change == "mixed":
                    manifest["patch"]["crates-io"]["graphica"]["rev"] = OTHER
                elif change == "wrong-table":
                    manifest["patch"][URL] = manifest["patch"].pop("crates-io")
                elif change == "branch":
                    spec = manifest["patch"]["crates-io"]["symbolica"]
                    spec.pop("rev")
                    spec["branch"] = "main"
                elif change == "direct-git":
                    manifest["workspace"]["dependencies"]["symbolica"]["git"] = URL
                else:
                    manifest["workspace"]["dependencies"]["numerica"] = "=3.0.0"
                with self.assertRaises(origins.InvalidFamily):
                    origins.validate(manifest, lock)


class ExistingFeatureGuardTests(unittest.TestCase):
    def test_actual_recipe_keeps_native_and_wasm_feature_guards(self):
        root = SCRIPT.parent.parent
        text = (root / "justfile").read_text()
        body = text.split("check-symbolica-feature-isolation:\n", 1)[1].split(
            "\n\n", 1
        )[0]
        body = "\n".join(line.removeprefix("    ") for line in body.splitlines())
        for wrong_mix in ("portable", "missing-native"):
            with (
                self.subTest(wrong_mix=wrong_mix),
                tempfile.TemporaryDirectory() as temporary,
            ):
                fixture = Path(temporary)
                hack = fixture / "crates/gammaloop-workspace-hack"
                hack.mkdir(parents=True)
                features = ["integer-gmp", "float-mpfr", "native_code_generation"]
                if wrong_mix == "portable":
                    features.append("wasm")
                else:
                    features.remove("native_code_generation")
                (hack / "Cargo.toml").write_text("features = " + json.dumps(features))
                command = body.replace(
                    'root="{{ justfile_directory() }}"',
                    "root=" + shlex.quote(str(fixture)),
                )
                result = subprocess.run(
                    ["bash", "-c", command], capture_output=True, text=True, check=False
                )
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("native GammaLoop workspace hack", result.stderr)
                self.assertNotIn("check-symbolica-lock-origins.py", result.stderr)


if __name__ == "__main__":
    unittest.main()
