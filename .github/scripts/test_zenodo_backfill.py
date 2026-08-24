from __future__ import annotations

import copy
import importlib.util
import io
import json
import os
import subprocess
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).parents[2]

ZENODO_SPEC = importlib.util.spec_from_file_location(
    "zenodo", Path(__file__).with_name("zenodo.py")
)
assert ZENODO_SPEC and ZENODO_SPEC.loader
zenodo = importlib.util.module_from_spec(ZENODO_SPEC)
sys.modules[ZENODO_SPEC.name] = zenodo
ZENODO_SPEC.loader.exec_module(zenodo)

BACKFILL_SPEC = importlib.util.spec_from_file_location(
    "gammaloop_zenodo_backfill", Path(__file__).with_name("zenodo_backfill.py")
)
assert BACKFILL_SPEC and BACKFILL_SPEC.loader
backfill = importlib.util.module_from_spec(BACKFILL_SPEC)
sys.modules[BACKFILL_SPEC.name] = backfill
BACKFILL_SPEC.loader.exec_module(backfill)


def tar_gz(entries: dict[str, bytes]) -> bytes:
    output = io.BytesIO()
    with tarfile.open(fileobj=output, mode="w:gz") as archive:
        for name, data in entries.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            archive.addfile(info, io.BytesIO(data))
    return output.getvalue()


class ManifestFixture(unittest.TestCase):
    def setUp(self):
        self.registry = zenodo.load_registry(ROOT / ".zenodo/records.json")
        self.document = json.loads((ROOT / ".zenodo/backfill.json").read_text())

    def load(self, document):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "backfill.json"
            path.write_text(json.dumps(document))
            return backfill.BackfillManifest.load(path, self.registry)


class ManifestTests(ManifestFixture):
    def test_real_manifest_is_exactly_twenty_events(self):
        manifest = backfill.BackfillManifest.load(
            ROOT / ".zenodo/backfill.json", self.registry
        )
        self.assertEqual(len(tuple(manifest.events())), 20)
        self.assertEqual(
            tuple(family.name for family in manifest.families),
            backfill.EXPECTED_FAMILY_ORDER,
        )

    def test_rejects_deleted_historical_version(self):
        del self.document["families"]["gammaloop"]["versions"][3]
        with self.assertRaisesRegex(zenodo.ReleaseError, "version list"):
            self.load(self.document)

    def test_rejects_allowlist_extension(self):
        self.document["redirect_host_allowlists"]["pypi"].append("example.test")
        with self.assertRaisesRegex(zenodo.ReleaseError, "allowlist differs"):
            self.load(self.document)

    def test_rejects_component_set_change(self):
        self.document["families"]["spenso"]["versions"][0]["component_versions"].pop(
            "spynso3"
        )
        with self.assertRaisesRegex(zenodo.ReleaseError, "component set"):
            self.load(self.document)

    def test_rejects_crate_not_bound_to_component_version(self):
        artifact = self.document["families"]["vakint"]["versions"][0]["artifacts"][0]
        artifact["version"] = "0.1.1"
        artifact["filename"] = "vakint-0.1.1.crate"
        artifact["url"] = "https://crates.io/api/v1/crates/vakint/0.1.1/download"
        with self.assertRaisesRegex(zenodo.ReleaseError, "component versions"):
            self.load(self.document)

    def test_rejects_wrong_git_scope_and_ref(self):
        artifact = self.document["families"]["gammaloop"]["versions"][3]["artifacts"][0]
        artifact["archive_paths"] = ["crates"]
        with self.assertRaisesRegex(zenodo.ReleaseError, "full tree"):
            self.load(self.document)

    def test_rejects_deleted_frozen_git_hash(self):
        artifact = self.document["families"]["gammaloop-monorepo"]["versions"][0][
            "artifacts"
        ][0]
        del artifact["sha256"]
        with self.assertRaisesRegex(zenodo.ReleaseError, "schema differs"):
            self.load(self.document)

    def test_rejects_wheel_filename(self):
        with self.assertRaisesRegex(zenodo.ReleaseError, "wheel"):
            backfill.BackfillManifest._safe_filename("release.whl", "fixture")


class FakeTransport:
    def __init__(self, responses):
        self.responses = list(responses)
        self.urls = []

    def request(self, method, url, headers, body):
        self.urls.append(url)
        return self.responses.pop(0)


class DownloaderTests(unittest.TestCase):
    def test_follows_only_allowlisted_redirects(self):
        transport = FakeTransport(
            [
                zenodo.Response(
                    302, {"Location": "https://static.crates.io/file"}, b""
                ),
                zenodo.Response(200, {}, b"crate"),
            ]
        )
        result = backfill.Downloader(transport).fetch(
            "https://crates.io/download", ("crates.io", "static.crates.io")
        )
        self.assertEqual(result, b"crate")
        self.assertEqual(transport.urls[-1], "https://static.crates.io/file")

    def test_rejects_redirect_to_untrusted_host(self):
        transport = FakeTransport(
            [zenodo.Response(302, {"Location": "https://attacker.test/file"}, b"")]
        )
        with self.assertRaisesRegex(zenodo.ReleaseError, "untrusted"):
            backfill.Downloader(transport).fetch(
                "https://crates.io/download", ("crates.io",)
            )

    def test_rejects_redirect_with_userinfo_or_port(self):
        for url in ("https://user@crates.io/file", "https://crates.io:443/file"):
            with self.subTest(url=url), self.assertRaises(zenodo.ReleaseError):
                backfill.Downloader._trusted(url, ("crates.io",))


class CanonicalArchiveTests(unittest.TestCase):
    def test_crate_vcs_info_is_verified(self):
        artifact = {
            "filename": "fixture-1.2.3.crate",
            "crate": "fixture",
            "version": "1.2.3",
            "cargo_vcs_info": {"sha1": "a" * 40, "path_in_vcs": "crates/fixture"},
        }
        data = tar_gz(
            {
                "fixture-1.2.3/.cargo_vcs_info.json": json.dumps(
                    {"git": {"sha1": "a" * 40}, "path_in_vcs": "crates/fixture"}
                ).encode(),
                "fixture-1.2.3/Cargo.toml": b'[package]\nname = "fixture"\nversion = "1.2.3"\n',
            }
        )
        backfill.BackfillPreparer._verify_crate(data, artifact)
        artifact["cargo_vcs_info"]["sha1"] = "b" * 40
        with self.assertRaisesRegex(zenodo.ReleaseError, "differs"):
            backfill.BackfillPreparer._verify_crate(data, artifact)

    def test_sdist_package_metadata_is_verified(self):
        artifact = {
            "filename": "gammaloop-0.3.2.tar.gz",
            "project": "gammaloop",
        }
        data = tar_gz(
            {
                "gammaloop-0.3.2/PKG-INFO": b"Metadata-Version: 2.1\nName: GammaLoop\nVersion: 0.3.2\n"
            }
        )
        backfill.BackfillPreparer._verify_sdist(data, artifact)
        artifact["filename"] = "gammaloop-0.3.1.tar.gz"
        with self.assertRaisesRegex(zenodo.ReleaseError, "differs"):
            backfill.BackfillPreparer._verify_sdist(data, artifact)

    def test_generated_git_archive_is_deterministic_and_ref_pinned(self):
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
            (repository / "source.txt").write_text("historical\n")
            (repository / "dangling").symlink_to("target/release/cli")
            subprocess.run(["git", "add", "."], cwd=repository, check=True)
            environment = {
                **os.environ,
                "GIT_AUTHOR_DATE": "2024-01-02T03:04:05+00:00",
                "GIT_COMMITTER_DATE": "2024-01-02T03:04:05+00:00",
            }
            subprocess.run(
                ["git", "commit", "-qm", "fixture"],
                cwd=repository,
                env=environment,
                check=True,
            )
            subprocess.run(["git", "tag", "v1.2.3"], cwd=repository, check=True)
            commit = zenodo._commit(repository, "refs/tags/v1.2.3")
            timestamp = int(
                zenodo._git(repository, "show", "-s", "--format=%ct", commit)
            )
            commit_date = (
                zenodo._git(repository, "show", "-s", "--format=%cI", commit)
                .decode()
                .strip()
            )
            artifact = {
                "ref": "refs/tags/v1.2.3",
                "tag_object_type": "commit",
                "commit": commit,
                "source_date_epoch": timestamp,
                "commit_date": commit_date,
                "archive_paths": ["."],
                "archive_prefix": "fixture-1.2.3/",
            }
            preparer = backfill.BackfillPreparer(repository, None, None)
            first = repository / "first.tar.gz"
            second = repository / "second.tar.gz"
            self.assertEqual(
                preparer._git_archive(artifact, first),
                preparer._git_archive(artifact, second),
            )
            self.assertEqual(first.read_bytes(), second.read_bytes())
            with tarfile.open(first, "r:gz") as archive:
                dangling = archive.getmember("fixture-1.2.3/dangling")
                self.assertTrue(dangling.issym())
                self.assertEqual(dangling.linkname, "target/release/cli")
            artifact["commit"] = "0" * 40
            with self.assertRaisesRegex(zenodo.ReleaseError, "changed commit"):
                preparer._git_archive(artifact, first)

    def test_git_archive_rejects_symlink_escaping_prefix(self):
        output = io.BytesIO()
        with tarfile.open(fileobj=output, mode="w:") as archive:
            link = tarfile.TarInfo("fixture/link")
            link.type = tarfile.SYMTYPE
            link.linkname = "../outside"
            archive.addfile(link)
        with self.assertRaisesRegex(zenodo.ReleaseError, "symlink escapes"):
            backfill.BackfillPreparer._verify_git_tar(output.getvalue(), "fixture/")


class PlanSecurityTests(ManifestFixture):
    def test_prepare_requires_plan_directly_in_output(self):
        manifest = backfill.BackfillManifest.load(
            ROOT / ".zenodo/backfill.json", self.registry
        )
        preparer = backfill.BackfillPreparer(ROOT, manifest, self.registry)
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "output"
            with self.assertRaisesRegex(zenodo.ReleaseError, "directly inside"):
                preparer.prepare(output, output / "nested" / "plan.json")

    def test_resolved_path_rejects_symlink_escape(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            trusted = root / "trusted"
            trusted.mkdir()
            outside = root / "outside"
            outside.write_text("untrusted")
            (trusted / "link").symlink_to(outside)
            with self.assertRaisesRegex(zenodo.ReleaseError, "escapes"):
                backfill.PreparedBackfill._resolved(trusted, "link", "asset")

    def test_frozen_git_hash_is_enforced_during_prepare(self):
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
            (repository / "file").write_text("content")
            subprocess.run(["git", "add", "."], cwd=repository, check=True)
            subprocess.run(
                ["git", "commit", "-qm", "fixture"], cwd=repository, check=True
            )
            subprocess.run(["git", "tag", "v1.0.0"], cwd=repository, check=True)
            commit = zenodo._commit(repository, "refs/tags/v1.0.0")
            artifact = {
                "kind": "git-archive",
                "filename": "source.tar.gz",
                "ref": "refs/tags/v1.0.0",
                "tag_object_type": "commit",
                "commit": commit,
                "commit_date": zenodo._git(
                    repository, "show", "-s", "--format=%cI", commit
                )
                .decode()
                .strip(),
                "source_date_epoch": int(
                    zenodo._git(repository, "show", "-s", "--format=%ct", commit)
                ),
                "archive_paths": ["."],
                "archive_prefix": "source/",
                "sha256": "0" * 64,
            }
            preparer = backfill.BackfillPreparer(
                repository, SimpleNamespace(hosts={}), None
            )
            version = backfill.HistoricalVersion(
                "fixture", "1.0.0", "2024-01-01", {}, (artifact,), None
            )
            output = repository / "output"
            output.mkdir()
            with self.assertRaisesRegex(
                zenodo.ReleaseError, "generated Git archive differs"
            ):
                preparer._prepare_artifact(version, artifact, output, repository)


class PublicationInvariantTests(ManifestFixture):
    def test_historical_components_link_only_to_umbrella_concept(self):
        manifest = backfill.BackfillManifest.load(
            ROOT / ".zenodo/backfill.json", self.registry
        )
        prepared = SimpleNamespace(root=ROOT, registry=self.registry, manifest=manifest)
        client = SimpleNamespace(production=True)
        publisher = backfill.BackfillPublisher(prepared, client)
        historical = manifest.family("spenso").versions[0]
        event = {
            "metadata_version": "spenso-v0.6.0",
            "publication_date": "2026-01-14",
        }
        metadata = publisher._metadata(event, historical, "10.5281/zenodo.42")
        relations = {
            (item["identifier"], item["relation"])
            for item in metadata["related_identifiers"]
        }
        self.assertIn(("10.5281/zenodo.42", "isPartOf"), relations)
        self.assertFalse(any(relation == "hasPart" for _, relation in relations))

    def test_retry_tolerates_next_version_draft_but_rejects_missing_history(self):
        published = {"id": 1, "submitted": True, "metadata": {"version": "v0.2.0"}}
        cloned_draft = {"id": 2, "submitted": False, "metadata": {"version": "v0.2.0"}}
        existing, draft = backfill.BackfillPublisher._lineage_state(
            [published, cloned_draft], "0.2.0"
        )
        self.assertEqual(existing["id"], 1)
        self.assertEqual(draft["id"], 2)
        with self.assertRaisesRegex(zenodo.CollisionError, "missing before"):
            backfill.BackfillPublisher._lineage_state(
                [{"id": 3, "submitted": True, "metadata": {"version": "v0.3.0"}}],
                "0.2.0",
            )

    def test_authenticated_linnet_metadata_is_fully_rechecked(self):
        manifest = backfill.BackfillManifest.load(
            ROOT / ".zenodo/backfill.json", self.registry
        )
        adoption = manifest.family("linnet").versions[0].adoption
        assert adoption is not None
        record = {
            "id": int(adoption["record_id"]),
            "conceptrecid": adoption["concept_recid"],
            "doi": adoption["record_doi"],
            "conceptdoi": adoption["concept_doi"],
            "created": adoption["published_at"],
            "metadata": {
                "title": adoption["title"],
                "version": adoption["metadata_version"],
                "license": adoption["license"],
                "creators": [{"name": name} for name in adoption["creators"]],
                "related_identifiers": [{"identifier": adoption["related_tag_url"]}],
            },
            "files": [
                {
                    "key": item["filename"],
                    "size": item["size"],
                    "checksum": f"md5:{item['md5']}",
                }
                for item in adoption["files"]
            ],
        }
        backfill.BackfillPreparer._verify_adoption_record(record, adoption)
        changed = copy.deepcopy(record)
        changed["metadata"]["title"] = "changed"
        with self.assertRaisesRegex(zenodo.ReleaseError, "metadata differs"):
            backfill.BackfillPreparer._verify_adoption_record(changed, adoption)

        ref_record = {
            "ref": adoption["tag_ref"],
            "object": {
                "type": adoption["tag_object_type"],
                "sha": adoption["tag_object"],
            },
        }
        tag_record = {
            "sha": adoption["tag_object"],
            "object": {"type": "commit", "sha": adoption["commit"]},
        }
        backfill.BackfillPreparer._verify_adoption_git(ref_record, tag_record, adoption)
        tag_record["object"]["sha"] = "0" * 40
        with self.assertRaisesRegex(zenodo.ReleaseError, "Git tag differs"):
            backfill.BackfillPreparer._verify_adoption_git(
                ref_record, tag_record, adoption
            )

    def test_global_preflight_detects_later_lineage_before_any_mutation(self):
        umbrella_version = backfill.HistoricalVersion(
            "gammaloop-monorepo", "0.2.0", "2024-10-08", {}, (), None
        )
        gamma_version = backfill.HistoricalVersion(
            "gammaloop", "0.0.1", "2023-09-13", {}, (), None
        )
        families = (
            backfill.HistoricalFamily(
                "gammaloop-monorepo",
                "gammaloop-monorepo",
                "new",
                (umbrella_version,),
                False,
            ),
            backfill.HistoricalFamily(
                "gammaloop", "gammaloop", "new", (gamma_version,), False
            ),
        )
        registry = SimpleNamespace(
            families={
                "gammaloop-monorepo": SimpleNamespace(
                    name="gammaloop-monorepo",
                    marker="https://example.test#umbrella",
                    concept_recid=None,
                ),
                "gammaloop": SimpleNamespace(
                    name="gammaloop",
                    marker="https://example.test#gamma",
                    concept_recid=None,
                ),
            }
        )
        prepared = SimpleNamespace(
            manifest=SimpleNamespace(families=families),
            registry=registry,
            events=lambda: iter(()),
        )

        class Client:
            production = True

            def __init__(self):
                self.mutations = 0

            def list_depositions(self):
                return [
                    {
                        "id": 2,
                        "conceptrecid": "20",
                        "submitted": True,
                        "metadata": {
                            "version": "v0.0.1",
                            "related_identifiers": [
                                {
                                    "identifier": "https://example.test#gamma",
                                    "relation": "isAlternateIdentifier",
                                }
                            ],
                        },
                    }
                ]

            def __getattr__(self, _name):
                def mutate(*_args, **_kwargs):
                    self.mutations += 1
                    raise AssertionError("mutation must not be reached")

                return mutate

        client = Client()
        publisher = backfill.BackfillPublisher(prepared, client)
        with self.assertRaisesRegex(zenodo.CollisionError, "after incomplete"):
            publisher.publish()
        self.assertEqual(client.mutations, 0)


if __name__ == "__main__":
    unittest.main()
