from __future__ import annotations

import hashlib
import importlib.util
import io
import json
import subprocess
import sys
import tempfile
import unittest
from argparse import Namespace
from contextlib import redirect_stdout
from pathlib import Path

MODULE_PATH = Path(__file__).with_name("zenodo.py")
SPEC = importlib.util.spec_from_file_location("gammaloop_zenodo", MODULE_PATH)
assert SPEC and SPEC.loader
zenodo = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = zenodo
SPEC.loader.exec_module(zenodo)


def record(
    identifier: int,
    concept: int,
    version: str,
    *,
    marker: str | None = None,
    published: bool = True,
    files: dict[str, str] | None = None,
):
    relations = []
    if marker:
        relations.append({"identifier": marker, "relation": "isAlternateIdentifier"})
    return {
        "id": identifier,
        "conceptrecid": concept,
        "submitted": published,
        "state": "done" if published else "unsubmitted",
        "metadata": {"version": version, "related_identifiers": relations},
        "files": [
            {"filename": name, "checksum": f"md5:{checksum}"}
            for name, checksum in (files or {}).items()
        ],
    }


class VersionTests(unittest.TestCase):
    def test_normalizes_supported_forms(self):
        for value in ("0.17.0", "v0.17.0", "linnet-v0.17.0"):
            with self.subTest(value=value):
                self.assertEqual(zenodo.normalize_version(value), "0.17.0")

    def test_preserves_semver_build_metadata(self):
        self.assertEqual(
            zenodo.normalize_version("spenso-v0.6.0+sandbox.123.2"),
            "0.6.0+sandbox.123.2",
        )

    def test_rejects_ambiguous_versions(self):
        for value in ("release-0.17.0", "vv0.17.0", "0.17", "01.2.3"):
            with self.subTest(value=value), self.assertRaises(zenodo.ReleaseError):
                zenodo.normalize_version(value)


class LineageTests(unittest.TestCase):
    def test_marker_selects_one_concept_lineage(self):
        marker = "https://example.test/repo#family"
        records = [
            record(1, 10, "v1.0.0", marker=marker),
            record(2, 10, "v1.1.0", marker=marker),
            record(3, 20, "v9.0.0"),
        ]
        self.assertEqual(
            [item["id"] for item in zenodo.select_lineage(records, marker)], [1, 2]
        )

    def test_explicit_concept_adopts_unmarked_linnet(self):
        records = [record(18429583, 15494393, "linnet-v0.17.0")]
        selected = zenodo.select_lineage(
            records, "https://example.test/linnet", "15494393"
        )
        self.assertEqual([item["id"] for item in selected], [18429583])
        self.assertEqual(zenodo._record_version(selected[0]), "0.17.0")

    def test_marker_collision_fails(self):
        marker = "https://example.test/repo#family"
        records = [
            record(1, 10, "v1.0.0", marker=marker),
            record(2, 20, "v2.0.0", marker=marker),
        ]
        with self.assertRaises(zenodo.CollisionError):
            zenodo.select_lineage(records, marker)

    def test_marker_and_explicit_concept_must_agree(self):
        marker = "https://example.test/repo#linnet"
        with self.assertRaisesRegex(zenodo.CollisionError, "not configured concept"):
            zenodo.select_lineage(
                [record(1, 10, "linnet-v0.17.0", marker=marker)], marker, "20"
            )


class ArchiveTests(unittest.TestCase):
    def test_archive_is_byte_deterministic(self):
        entries = [
            zenodo.ArchiveEntry("pkg/b.txt", b"b"),
            zenodo.ArchiveEntry("pkg/a.txt", b"a", 0o755),
        ]
        with tempfile.TemporaryDirectory() as temporary:
            first = Path(temporary) / "first.tar.gz"
            second = Path(temporary) / "second.tar.gz"
            first_hash = zenodo.write_archive(first, entries, 1_700_000_000)
            second_hash = zenodo.write_archive(second, reversed(entries), 1_700_000_000)
            self.assertEqual(first_hash, second_hash)
            self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_git_tree_materializes_safe_symlink_and_rejects_escape(self):
        with tempfile.TemporaryDirectory() as temporary:
            repository = Path(temporary)
            subprocess.run(["git", "init", "-q"], cwd=repository, check=True)
            subprocess.run(
                ["git", "config", "user.email", "test@example.test"],
                cwd=repository,
                check=True,
            )
            subprocess.run(
                ["git", "config", "user.name", "Test"], cwd=repository, check=True
            )
            (repository / "assets/models").mkdir(parents=True)
            (repository / "assets/models/data.json").write_text("{}\n")
            (repository / "crates/api/assets").mkdir(parents=True)
            (repository / "crates/api/assets/models").symlink_to(
                "../../../assets/models"
            )
            subprocess.run(["git", "add", "."], cwd=repository, check=True)
            subprocess.run(
                ["git", "commit", "-qm", "fixture"], cwd=repository, check=True
            )
            commit = zenodo._commit(repository, "HEAD")
            tree = zenodo._git_tree(repository, commit)
            self.assertEqual(tree["crates/api/assets/models/data.json"][0], b"{}\n")

            (repository / "crates/api/assets/models").unlink()
            (repository / "crates/api/assets/models").symlink_to("../../../../outside")
            subprocess.run(["git", "add", "."], cwd=repository, check=True)
            subprocess.run(
                ["git", "commit", "-qm", "unsafe"], cwd=repository, check=True
            )
            with self.assertRaisesRegex(zenodo.ReleaseError, "escapes repository"):
                zenodo._git_tree(repository, zenodo._commit(repository, "HEAD"))

    def test_crate_assets_use_independent_component_versions(self):
        tree = {
            "crates/spenso/Cargo.toml": (
                b'[package]\nname="spenso"\nversion="0.6.0"\n',
                0o644,
                None,
            ),
            "crates/macros/Cargo.toml": (
                b'[package]\nname="spenso-macros"\nversion="0.4.1"\n',
                0o644,
                None,
            ),
            "crates/clinnet/Cargo.toml": (
                b'[package]\nname="clinnet"\nversion="0.8.2"\n',
                0o644,
                None,
            ),
        }
        versions = zenodo._crate_versions(tree, ["spenso", "spenso-macros", "clinnet"])
        self.assertEqual(
            versions,
            {"spenso": "0.6.0", "spenso-macros": "0.4.1", "clinnet": "0.8.2"},
        )

    def test_prepare_builds_monorepo_and_scoped_archives_without_wheels(self):
        with tempfile.TemporaryDirectory() as temporary:
            repository = Path(temporary)
            subprocess.run(["git", "init", "-q"], cwd=repository, check=True)
            subprocess.run(
                ["git", "config", "user.email", "test@example.test"],
                cwd=repository,
                check=True,
            )
            subprocess.run(
                ["git", "config", "user.name", "Test"], cwd=repository, check=True
            )
            (repository / ".zenodo/metadata").mkdir(parents=True)
            (repository / "crates/api").mkdir(parents=True)
            (repository / "crates/core").mkdir(parents=True)
            (repository / "crates/api/Cargo.toml").write_text(
                '[package]\nname="api"\nversion="1.2.3"\n'
            )
            (repository / "crates/core/Cargo.toml").write_text(
                '[package]\nname="core"\nversion="1.2.3"\n'
            )
            (repository / "LICENSE.md").write_text("MIT fixture\n")
            citation = "cff-version: 1.2.0\ntitle: fixture\ntype: software\n"
            (repository / "CITATION.cff").write_text(citation)
            (repository / ".zenodo/monorepo.cff").write_text(citation)
            metadata = {
                "title": "fixture",
                "upload_type": "software",
                "access_right": "open",
                "license": "mit",
                "creators": [{"name": "Tester, Test"}],
                "grants": [{"id": "grant"}],
                "communities": [{"identifier": "alphaloop"}],
                "related_identifiers": [],
            }
            for name in ("monorepo", "gamma"):
                (repository / f".zenodo/metadata/{name}.json").write_text(
                    json.dumps(metadata)
                )
            family_common = {
                "title": "fixture",
                "tag": "v{version}",
                "marker": "https://example.test/repo#family",
                "concept_recid": None,
                "full_repository": False,
                "archive_paths": [],
                "crates": [],
                "include_sdist": False,
                "production_baseline_version": "1.2.3",
            }
            registry = {
                "schema": 1,
                "repository": "https://example.test/repo",
                "community": "alphaloop",
                "grant": "grant",
                "families": {
                    "monorepo": {
                        **family_common,
                        "version_source": "release",
                        "metadata_file": ".zenodo/metadata/monorepo.json",
                        "citation_file": ".zenodo/monorepo.cff",
                        "full_repository": True,
                        "marker": "https://example.test/repo#monorepo",
                    },
                    "gamma": {
                        **family_common,
                        "version_source": "crates/api/Cargo.toml",
                        "metadata_file": ".zenodo/metadata/gamma.json",
                        "citation_file": "CITATION.cff",
                        "archive_paths": ["crates/api", "crates/core", "LICENSE.md"],
                        "crates": ["api", "core"],
                        "include_sdist": True,
                        "marker": "https://example.test/repo#gamma",
                    },
                },
            }
            (repository / ".zenodo/records.json").write_text(json.dumps(registry))
            subprocess.run(["git", "add", "."], cwd=repository, check=True)
            subprocess.run(
                ["git", "commit", "-qm", "fixture"], cwd=repository, check=True
            )
            subprocess.run(
                ["git", "tag", "-a", "v1.2.3", "-m", "fixture 1.2.3"],
                cwd=repository,
                check=True,
            )
            commit = zenodo._commit(repository, "HEAD")
            assets = repository / "release-dist"
            assets.mkdir()
            payloads = {
                "api-1.2.3.crate": b"api",
                "core-1.2.3.crate": b"core",
                "gammaloop-1.2.3.tar.gz": b"sdist",
                "gammaloop-1.2.3.whl": b"wheel",
            }
            for name, body in payloads.items():
                (assets / name).write_bytes(body)
            (assets / "SHA256SUMS").write_text(
                "".join(
                    f"{hashlib.sha256(body).hexdigest()}  {name}\n"
                    for name, body in sorted(payloads.items())
                )
            )
            output = repository / "zenodo-payload"
            with redirect_stdout(io.StringIO()):
                zenodo._prepare(
                    Namespace(
                        repository_root=str(repository),
                        registry=".zenodo/records.json",
                        release_version="1.2.3",
                        release_tag="v1.2.3",
                        release_commit=commit,
                        assets_dir=str(assets),
                        output_dir=str(output),
                        plan=str(output / "plan.json"),
                        version_suffix=None,
                    )
                )
            plan = json.loads((output / "plan.json").read_text())
            self.assertEqual(
                plan["families"]["gamma"]["component_versions"],
                {"api": "1.2.3", "core": "1.2.3"},
            )
            self.assertEqual(
                {item["source"] for item in plan["families"]["gamma"]["assets"]},
                {"api-1.2.3.crate", "core-1.2.3.crate", "gammaloop-1.2.3.tar.gz"},
            )
            for family in ("monorepo", "gamma"):
                archive = output / plan["families"][family]["base_archive"]
                names = [entry.name for entry in zenodo.read_archive(archive)]
                self.assertEqual(len(names), len(set(names)))
                self.assertTrue(any(name.endswith("/LICENSE.md") for name in names))
            loaded_registry = zenodo.load_registry(repository / ".zenodo/records.json")
            final = zenodo._finalize_payload(
                output / "plan.json",
                plan["families"]["gamma"],
                loaded_registry.families["gamma"],
                "10.5072/zenodo.2",
                "10.5072/zenodo.1",
                citation,
            )
            self.assertEqual(
                set(final),
                {
                    "gamma-source-1.2.3.tar.gz",
                    "api-1.2.3.crate",
                    "core-1.2.3.crate",
                    "gammaloop-1.2.3.tar.gz",
                    "CITATION.cff",
                    "LICENSE.md",
                    "PROVENANCE.json",
                    "SHA256SUMS",
                },
            )
            self.assertNotIn(".whl", "".join(final))
            self.assertIn("10.5072/zenodo.2", final["CITATION.cff"].read_text())
            checksums = final["SHA256SUMS"].read_text()
            self.assertNotIn("SHA256SUMS", checksums)
            self.assertEqual(len(checksums.splitlines()), len(final) - 1)


class MetadataTests(unittest.TestCase):
    def test_relations_use_component_versions_and_umbrella_concept(self):
        concepts = {
            "gammaloop-monorepo": "10.1/concept-root",
            "spenso": "10.1/concept-s",
        }
        versions = {"gammaloop-monorepo": "10.1/root-v", "spenso": "10.1/spenso-v"}
        umbrella = zenodo.merge_relations(
            [], "https://marker/root", "gammaloop-monorepo", concepts, versions
        )
        self.assertIn(zenodo._relation("10.1/spenso-v", "hasPart"), umbrella)
        self.assertNotIn(zenodo._relation("10.1/concept-s", "hasPart"), umbrella)
        component = zenodo.merge_relations(
            [], "https://marker/s", "spenso", concepts, versions
        )
        self.assertIn(zenodo._relation("10.1/concept-root", "isPartOf"), component)

    def test_relation_merge_does_not_duplicate_legacy_relation(self):
        legacy = zenodo._relation("10.5281/zenodo.15913113", "isDerivedFrom")
        result = zenodo.merge_relations([legacy], "https://marker/s", "spenso", {}, {})
        self.assertEqual(result.count(legacy), 1)

    def test_citation_merges_existing_concept_identifier(self):
        template = 'cff-version: 1.2.0\ntitle: linnet\nidentifiers:\n  - type: doi\n    value: "10.1/concept"\n'
        rendered = zenodo.render_citation(
            template, "0.17.0", "2026-08-14", "10.1/version", "10.1/concept"
        )
        self.assertEqual(rendered.count("10.1/concept"), 1)
        self.assertIn("10.1/version", rendered)
        self.assertIn('version: "0.17.0"', rendered)

    def test_stable_metadata_requires_community_and_relations(self):
        expected = {
            "version": "spenso-v0.6.0",
            "publication_date": "2026-08-14",
            "title": "spenso",
            "creators": [{"name": "Tester, Test"}],
            "license": "mit",
            "grants": [{"id": "grant"}],
            "communities": [{"identifier": "alphaloop"}],
            "related_identifiers": [
                {"identifier": "https://marker/s", "relation": "isAlternateIdentifier"}
            ],
        }
        deposited = {
            "metadata": json.loads(json.dumps(expected)),
            "submitted": True,
            "state": "done",
        }
        zenodo.assert_stable_metadata(deposited, expected)
        deposited["metadata"]["communities"] = []
        with self.assertRaisesRegex(zenodo.CollisionError, "communities"):
            zenodo.assert_stable_metadata(deposited, expected)

    def test_sandbox_metadata_never_emits_production_doi_relations(self):
        metadata = {
            "related_identifiers": [
                {"identifier": "10.5281/zenodo.1", "relation": "isDerivedFrom"},
                {"identifier": "10.5072/zenodo.2", "relation": "isPartOf"},
                {
                    "identifier": "https://example.test/marker",
                    "relation": "isAlternateIdentifier",
                },
            ]
        }
        sandbox = zenodo.metadata_for_environment(metadata, False)
        self.assertFalse(
            any(
                relation["identifier"].startswith("10.5281")
                for relation in sandbox["related_identifiers"]
            )
        )
        self.assertEqual(
            len(zenodo.metadata_for_environment(metadata, True)["related_identifiers"]),
            3,
        )


class CollisionTests(unittest.TestCase):
    def test_same_version_same_files_is_idempotent(self):
        existing = record(1, 10, "linnet-v0.17.0", files={"source.tar.gz": "abcd"})
        zenodo.assert_same_version_files(existing, "0.17.0", {"source.tar.gz": "abcd"})

    def test_same_version_different_files_fails(self):
        existing = record(1, 10, "v1.0.0", files={"source.tar.gz": "old"})
        with self.assertRaises(zenodo.CollisionError):
            zenodo.assert_same_version_files(
                existing, "1.0.0", {"source.tar.gz": "new"}
            )

    def test_published_match_with_leftover_draft_fails_closed(self):
        lineage = [
            record(1, 10, "v1.0.0"),
            record(2, 10, "v1.0.0", published=False),
        ]
        with self.assertRaisesRegex(zenodo.CollisionError, "leftover draft"):
            zenodo._matching_published(lineage, "1.0.0")

    def test_production_preflight_requires_the_exact_migration_baseline(self):
        family = zenodo.Family(
            "fixture",
            "fixture",
            "release",
            "v{version}",
            "https://example.test/repo#fixture",
            None,
            "metadata.json",
            "CITATION.cff",
            True,
            (),
            (),
            False,
            production_baseline_version="1.0.0",
        )
        registry = zenodo.Registry(
            "https://example.test/repo", "test", "grant", {"fixture": family}
        )
        with self.assertRaisesRegex(zenodo.ReleaseError, "fixture 1.0.0"):
            zenodo._production_preflight(
                registry, {"fixture": [record(2, 10, "v2.0.0")]}
            )
        zenodo._production_preflight(
            registry,
            {"fixture": [record(1, 10, "v1.0.0"), record(2, 10, "v2.0.0")]},
        )

    def test_explicit_linnet_baseline_is_the_only_legacy_skip(self):
        family = zenodo.Family(
            "linnet",
            "linnet",
            "crates/linnet/Cargo.toml",
            "linnet-v{version}",
            "https://example.test/repo#linnet",
            "15494393",
            "metadata.json",
            "CITATION.cff",
            False,
            ("crates/linnet",),
            ("linnet",),
            False,
        )
        legacy = record(1, 15494393, "linnet-v0.17.0")
        self.assertTrue(zenodo.is_legacy_linnet_match(family, legacy, "0.17.0", True))
        self.assertFalse(zenodo.is_legacy_linnet_match(family, legacy, "0.17.0", False))


class FakeTransport:
    def __init__(self, replies):
        self.replies = list(replies)
        self.calls = []

    def request(self, method, url, headers, body):
        self.calls.append((method, url, headers, body))
        expected_method, expected_url, response = self.replies.pop(0)
        if method != expected_method or url != expected_url:
            raise AssertionError(
                f"expected {expected_method} {expected_url}, got {method} {url}"
            )
        return response


def response(status, document=None):
    return zenodo.Response(
        status, {}, json.dumps(document).encode() if document is not None else b""
    )


class ClientTests(unittest.TestCase):
    def test_plan_environment_must_match_api(self):
        production = zenodo.ZenodoClient("token", zenodo.PRODUCTION_API)
        sandbox = zenodo.ZenodoClient("token", zenodo.SANDBOX_API)
        production.validate_plan_environment({"release": {"sandbox_suffix": None}})
        sandbox.validate_plan_environment({"release": {"sandbox_suffix": "123.1"}})
        with self.assertRaisesRegex(zenodo.ReleaseError, "sandbox plan"):
            production.validate_plan_environment(
                {"release": {"sandbox_suffix": "123.1"}}
            )
        with self.assertRaisesRegex(zenodo.ReleaseError, "unsuffixed plan"):
            sandbox.validate_plan_environment({"release": {"sandbox_suffix": None}})

    def test_newversion_follows_latest_draft(self):
        base = zenodo.SANDBOX_API
        latest = {"id": 4}
        transport = FakeTransport(
            [
                (
                    "POST",
                    f"{base}/deposit/depositions/4/actions/newversion",
                    response(
                        201,
                        {"links": {"latest_draft": f"{base}/deposit/depositions/5"}},
                    ),
                ),
                ("GET", f"{base}/deposit/depositions/5", response(200, {"id": 5})),
            ]
        )
        client = zenodo.ZenodoClient("secret", base, transport)
        self.assertEqual(client.new_version(latest)["id"], 5)

    def test_hostile_api_link_is_rejected_before_token_egress(self):
        transport = FakeTransport([])
        client = zenodo.ZenodoClient("secret", zenodo.PRODUCTION_API, transport)
        with self.assertRaisesRegex(zenodo.ReleaseError, "untrusted"):
            client.get("https://attacker.example/api/deposit/depositions/1")
        self.assertEqual(transport.calls, [])

    def test_sandbox_doi_fallback_uses_sandbox_prefix(self):
        client = zenodo.ZenodoClient("secret", zenodo.SANDBOX_API, FakeTransport([]))
        self.assertEqual(client.concept_doi({"conceptrecid": 42}), "10.5072/zenodo.42")
        self.assertNotIn("10.5281", client.concept_doi({"conceptrecid": 42}))

    def test_replace_files_removes_all_inherited_files_then_uploads_exact_set(self):
        base = zenodo.SANDBOX_API
        bucket = "https://sandbox.zenodo.org/api/files/bucket-1"
        draft = {
            "id": 8,
            "links": {"bucket": bucket},
            "files": [
                {"links": {"self": f"{base}/deposit/depositions/8/files/1"}},
                {"links": {"self": f"{base}/deposit/depositions/8/files/2"}},
            ],
        }
        final = {
            "id": 8,
            "files": [
                {"filename": "CITATION.cff", "checksum": "md5:a"},
                {"filename": "source.tar.gz", "checksum": "md5:b"},
            ],
        }
        transport = FakeTransport(
            [
                ("DELETE", f"{base}/deposit/depositions/8/files/1", response(204)),
                ("DELETE", f"{base}/deposit/depositions/8/files/2", response(204)),
                ("PUT", f"{bucket}/CITATION.cff", response(200, {})),
                ("PUT", f"{bucket}/source.tar.gz", response(200, {})),
                ("GET", f"{base}/deposit/depositions/8", response(200, final)),
            ]
        )
        with tempfile.TemporaryDirectory() as temporary:
            citation = Path(temporary) / "CITATION.cff"
            source = Path(temporary) / "source.tar.gz"
            citation.write_bytes(b"citation")
            source.write_bytes(b"source")
            client = zenodo.ZenodoClient("secret", base, transport)
            self.assertEqual(
                client.replace_files(
                    draft, {source.name: source, citation.name: citation}
                ),
                final,
            )
        self.assertTrue(
            all(call[2]["Authorization"] == "Bearer secret" for call in transport.calls)
        )


if __name__ == "__main__":
    unittest.main()
