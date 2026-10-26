"""Append authenticated 5L cache, suite and paired-control evidence to a NEW report.

This is a report adapter, not an authenticator or a reduction runner. The native
drivers own proof validation. Functional records never establish performance.
Nothing is written unless every supplied receipt passes its scoped checks.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import re
import statistics
from datetime import datetime, timedelta
from pathlib import Path

PRIMARY = {
    "quadratic-four-rung-uv": "quadratic-transverse",
    "zero-momentum-diagnostic": "zero-momentum-metric",
}
METRICS = (
    "data_loading_seconds",
    "cache_loading_seconds",
    "owner_proof_validation_seconds",
    "owner_preparation_seconds",
    "saved_application_including_owner_seconds",
    "ibp_excluding_owner_seconds",
    "evaluation_seconds",
    "end_to_end_seconds",
)
STAGES = (
    "frontend_canonicalization",
    "frontend_tensor_prepass",
    "frontend_admitted_witness_search",
    "parent_transport_and_scalar_lowering",
    "saved_rule_application",
    "terminal_geometry_and_lower_factors",
)
LIGHT_CACHE_LIMIT_BYTES = 95_000_000
SOURCE_PROOF_LIMIT_BYTES = 100 * 1024 * 1024
CUMULATIVE_SCHEMA = "vakint.cumulative-native-cache-report.v1"
LIGHT_SCHEMA = "vakint.strict-light-cache-export-validation.v1"


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def same_json_value(before, after):
    """Preserve JSON types: Python equality otherwise accepts true == 1 == 1.0."""
    return json.dumps(before, sort_keys=True, allow_nan=False) == json.dumps(
        after, sort_keys=True, allow_nan=False
    )


def evidence(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": digest(path)}


def hashed_file(reference, base):
    require(isinstance(reference, dict), "missing hashed evidence reference")
    path = Path(reference["path"])
    path = path if path.is_absolute() else base / path
    require(digest(path) == reference["sha256"], f"changed evidence: {path}")
    return path.resolve()


def hashed_record(reference, base):
    path = hashed_file(reference, base)
    return path, read(path)


def validate_physical_scope(scope):
    """A scope partitions known examples; it never substitutes for their receipts."""
    require(
        isinstance(scope, dict) and set(scope) == {"required", "deferred"},
        "physical scope needs exactly required and deferred examples",
    )
    for kind in ("required", "deferred"):
        values = scope[kind]
        require(
            isinstance(values, list)
            and all(isinstance(value, str) and value in PRIMARY for value in values)
            and len(values) == len(set(values)),
            "unknown or duplicate physical scope example",
        )
    require(
        scope["required"]
        and not set(scope["required"]) & set(scope["deferred"])
        and set(scope["required"]) | set(scope["deferred"]) == set(PRIMARY),
        "physical scope must partition every example with a nonempty required set",
    )
    return copy.deepcopy(scope)


def load_physical_scope(path):
    if path is None:
        return None, None
    document = read(path)
    require(
        isinstance(document, dict)
        and document.get("schema") == "vakint.benchmark-physical-scope.v1"
        and isinstance(document.get("reason"), str)
        and document["reason"].strip(),
        "physical scope needs an explicit reason and known schema",
    )
    scope = validate_physical_scope(document.get("physical_scope"))
    return scope, {"evidence": evidence(path), "raw_receipt": document}


def physical_completion(completed, scope=None):
    """Keep historical scope-less reports bound to both physical examples."""
    require(set(completed) <= set(PRIMARY), "unknown completed physical example")
    effective = (
        validate_physical_scope(scope)
        if scope is not None
        else {"required": list(PRIMARY), "deferred": []}
    )
    result = {
        "completed_physical_examples": sorted(set(completed)),
        "physical_examples_complete": set(effective["required"]).issubset(completed),
    }
    if scope is not None:
        result["physical_scope"] = effective
    return result


def physical_scope_note(authorization):
    reference = authorization["evidence"]
    return (
        "Explicit physical scope: "
        + authorization["raw_receipt"]["reason"]
        + f" Scope artifact {reference['path']} (SHA-256 {reference['sha256']}). "
        + "Deferred examples remain unmeasured for this run."
    )


def native_cache_identity(directory):
    """Match the collector's opaque index identity, without decoding native data."""
    directory = Path(directory)
    require(
        directory.is_dir() and not directory.is_symlink(),
        "native cache needs a regular directory",
    )
    bundle, pack = directory / "bundle.rrindex", directory / "cache.rrindex"

    def present(path):
        return path.exists() or path.is_symlink()

    require(
        present(bundle) != present(pack),
        "native cache needs exactly one bundle or pack index",
    )
    names = ["cache.rrindex"]
    if present(bundle):
        bindings = directory / "bindings"
        require(
            bindings.is_dir() and not bindings.is_symlink(),
            "invalid native bindings directory",
        )
        names = ["bundle.rrindex"]
        for child in bindings.iterdir():
            require(
                re.fullmatch(r"[0-9a-f]{64}", child.name)
                and child.is_dir()
                and not child.is_symlink(),
                "invalid native binding namespace",
            )
            names.append(f"bindings/{child.name}/cache.rrindex")
        require(len(names) > 1, "native bundle has no child indexes")
    indexes = {}
    for name in sorted(names):
        path = directory / name
        require(
            path.is_file() and not path.is_symlink(),
            "native index needs a regular file",
        )
        indexes[name] = {"bytes": path.stat().st_size, "sha256": digest(path)}
    identity = {
        "schema": "vakint.native-cache-index-identity.v1",
        "layout": "bundle" if present(bundle) else "pack",
        "indexes": indexes,
    }
    identity["sha256"] = hashlib.sha256(
        json.dumps(identity, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return validated_native_cache_identity(identity)


def validated_native_cache_identity(identity):
    """Check the same opaque index-metadata envelope for live and archived sources."""
    require(
        isinstance(identity, dict)
        and set(identity) == {"schema", "layout", "indexes", "sha256"}
        and identity["schema"] == "vakint.native-cache-index-identity.v1"
        and identity["layout"] in ("pack", "bundle")
        and isinstance(identity["indexes"], dict),
        "invalid native index identity envelope",
    )
    indexes = identity["indexes"]
    require(
        set(indexes) == {"cache.rrindex"}
        if identity["layout"] == "pack"
        else (
            "bundle.rrindex" in indexes
            and len(indexes) > 1
            and all(
                isinstance(name, str)
                and re.fullmatch(r"bindings/[0-9a-f]{64}/cache\.rrindex", name)
                for name in indexes
                if name != "bundle.rrindex"
            )
        ),
        "native index identity has an invalid layout/path inventory",
    )
    require(
        all(
            isinstance(item, dict)
            and set(item) == {"bytes", "sha256"}
            and type(item["bytes"]) is int
            and item["bytes"] > 0
            and isinstance(item["sha256"], str)
            and re.fullmatch(r"[0-9a-f]{64}", item["sha256"])
            for item in indexes.values()
        ),
        "native index identity has invalid file lengths/digests",
    )
    envelope = {key: value for key, value in identity.items() if key != "sha256"}
    require(
        identity["sha256"]
        == hashlib.sha256(
            json.dumps(envelope, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest(),
        "native index identity aggregate differs",
    )
    return identity


def cache_binding(producer):
    report = producer["raw_report"]
    if report.get("schema") in (CUMULATIVE_SCHEMA, LIGHT_SCHEMA):
        return report["semantic_binding"]
    return report["preparation"]["options"]["single_file"]["binding"]


def cache_identity(record):
    if record.get("schema") == LIGHT_SCHEMA:
        require(
            record.get("cache_native_index_identity") is None
            and record.get("cache_sha256") is None,
            "light producer must retain its selected-file identity separately from its source",
        )
        return record["output"]["cache"]["sha256"]
    native = record.get("cache_native_index_identity")
    require(
        not (native is not None and record.get("cache_sha256") is not None),
        "ambiguous file/native cache identity",
    )
    return native if native is not None else record.get("cache_sha256")


def light_attachment(statistics, producer):
    """Current cached workers validate every selected file record before attaching it."""
    if producer["raw_report"].get("schema") != LIGHT_SCHEMA:
        return
    require(
        isinstance(statistics, dict)
        and statistics.get("read_backing") is True
        and isinstance(statistics.get("file"), dict)
        and type(statistics["file"].get("records")) is int
        and statistics["file"]["records"]
        == producer["export_validation"]["unique_selected_records"]
        and all(
            type(statistics["file"].get(key)) is int
            and statistics["file"][key]
            == producer["raw_report"]["output"]["identity"][key]
            for key in ("compressed_bytes", "decoded_bytes")
        )
        and same_json_value(statistics.get("binding"), cache_binding(producer))
        and statistics.get("pack") is None
        and statistics.get("running") is None,
        "LIGHT consumer did not attach the exact strict selected file/count/binding",
    )


def artifact_binaries(artifact):
    """Use the compiler receipt's actual binaries, never a runtime source label."""
    binary = artifact.get("binary")
    if isinstance(binary, dict):
        return [binary]
    binaries = artifact.get("artifacts")
    require(isinstance(binaries, list) and binaries, "missing compiler artifacts")
    return binaries


COHORT_STAMPS = (
    ("main_revision", "VAKINT_CACHE_BUILD_GIT_REVISION"),
    ("native_revision", "VAKINT_CACHE_BUILD_RUSTRED_REVISION"),
    ("kernel_revision", "VAKINT_CACHE_BUILD_SYMBOLICA_REVISION"),
    ("main_source_snapshot_sha256", "VAKINT_CACHE_BUILD_MAIN_SOURCE_SNAPSHOT_SHA256"),
    (
        "native_source_snapshot_sha256",
        "VAKINT_CACHE_BUILD_NATIVE_SOURCE_SNAPSHOT_SHA256",
    ),
)


def validate_compiled_artifact(observed, artifact):
    """Bind a Native reader to its own build, independently of its producer."""
    cohort = observed["compiled_cohort"]
    require(
        all(
            artifact["stamps"].get(stamp) == cohort[field]
            and isinstance(cohort[field], str)
            and re.fullmatch(
                r"[0-9a-f]{64}" if field.endswith("_sha256") else r"[0-9a-f]{40}",
                cohort[field],
            )
            for field, stamp in COHORT_STAMPS
        )
        and observed["executable_sha256"]
        in {item["sha256"] for item in artifact_binaries(artifact)},
        "Native reader artifact or own compiled cohort differs",
    )
    for item in artifact_binaries(artifact):
        if item.get("immutable_path") and Path(item["immutable_path"]).exists():
            require(
                digest(item["immutable_path"]) == item["sha256"],
                "Native reader executable changed",
            )


def selected_inputs(workload):
    return {
        "source_census_sha256": workload["source_census_sha256"],
        **{
            kind: [
                {
                    "id": row["id"],
                    "powers": row["powers"],
                    "public_input_source_sha256": hashlib.sha256(
                        row["public_input"].encode()
                    ).hexdigest()
                    if row.get("public_input") is not None
                    else None,
                }
                for row in workload[kind]
            ]
            for kind in ("training", "tuning", "held_out")
        },
    }


class ReportEvidence:
    """Bind raw reports to existing launch/build and completed window receipts.

    A window receipt is an external acceptance decision, never manufactured by
    this adapter. Its explicit observations must agree with the referenced
    launchers and artifacts. Merely hashing a functional run does not qualify it.
    """

    def __init__(self, path):
        self.path = Path(path).resolve()
        self.index = read(path)
        require(
            self.index.get("schema") == "vakint.five-loop-report-evidence.v1",
            "unknown report evidence schema",
        )
        self.records = {}
        self.used = set()
        for item in self.index["records"]:
            record_path, record = hashed_record(item["record"], self.path.parent)
            key = str(record_path)
            require(key not in self.records, "duplicate report evidence")
            launch_path, launcher = hashed_record(item["launcher"], self.path.parent)
            artifact_path, artifact = hashed_record(item["artifact"], self.path.parent)
            self.records[key] = {
                "record": record,
                "launcher": launcher,
                "artifact": artifact,
                "references": {
                    "record": evidence(record_path),
                    "launcher": evidence(launch_path),
                    "artifact": evidence(artifact_path),
                },
                "window": item.get("window"),
                "case_id": item["case_id"],
            }

    def bind(self, path):
        record = self.records.get(str(Path(path).resolve()))
        require(record is not None, "report has no launcher/artifact evidence")
        self.used.add(str(Path(path).resolve()))
        return copy.deepcopy(record)

    @staticmethod
    def strict_record(reference, directory):
        def unique(pairs):
            result = {}
            for key, value in pairs:
                require(key not in result, "duplicate metadata JSON key")
                result[key] = value
            return result

        require(isinstance(reference, dict), "missing hashed metadata reference")
        candidate = Path(reference["path"])
        candidate = candidate if candidate.is_absolute() else directory / candidate
        require(
            candidate.is_file() and not candidate.is_symlink(),
            "metadata reference needs a regular file",
        )
        source, document = hashed_record(reference, directory)
        require(
            json.loads(source.read_text(), object_pairs_hook=unique) == document,
            "metadata JSON differs",
        )
        return source, document

    def light_cache_record(self, path, cache, physical_scope):
        """A selected file is produced by strict export, never by a finite workload."""
        path, report = self.strict_record(evidence(path), self.path.parent)
        acceptance_path, acceptance = self.strict_record(
            self.index.get("cache_requalification"), self.path.parent
        )
        relation_path, relation = self.strict_record(
            acceptance.get("data_metadata_requalification"), acceptance_path.parent
        )
        source_manifest, _ = self.strict_record(
            relation["source_manifest"], relation_path.parent
        )
        consumer_manifest, _ = self.strict_record(
            relation["consumer_manifest"], relation_path.parent
        )
        root = consumer_manifest.parent
        replacement = relation.get("prepared_cache_replacement", {})
        require(
            replacement.get("loop_count") == 5
            and replacement.get("export_validation", {}).get("sha256") == digest(path)
            and hashed_file(replacement["export_validation"], relation_path.parent)
            == path,
            "LIGHT producer differs from the dataset's selected export",
        )
        source_path, _ = self.strict_record(report["source"]["producer"], root)
        _, native = self.strict_record(report["validation"]["receipt"], root)
        source = cumulative_cache_record(
            source_path,
            None,
            physical_scope,
            source_identity=native["source_cache_native_index_identity"],
        )
        require(
            source["status"] == "validated-committed-snapshot",
            "LIGHT needs a validated source snapshot, not prefix preservation alone",
        )
        validation = self.light_cache_export(
            path, report, source_manifest, root, source, self.strict_record
        )
        phases = validation["phase_seconds"]
        export_phases = ("encode_selected_completed_cache", "write_selected_cache")
        require(
            Path(cache).is_file()
            and not Path(cache).is_symlink()
            and Path(cache).stat().st_size == validation["compressed_bytes"],
            "LIGHT runtime file differs from its strict export",
        )
        require(
            digest(cache) == validation["cache_sha256"],
            "LIGHT runtime file differs from its strict export",
        )
        _, metadata = self.strict_record(report["output"]["metadata"], root)
        self.light_cache_metadata(metadata, report, root)
        return {
            "evidence": evidence(path),
            "cache_file": evidence(cache),
            "raw_report": report,
            "source_producer": source,
            "export_validation": validation,
            "status": "validated-strict-selected-file",
            "production_seconds": {
                "selection": None,
                "preparation": phases.get("load_saved_application"),
                "export": number(sum(phases[key] for key in export_phases))
                if all(key in phases for key in export_phases)
                else None,
                "reimport": phases.get("reimport_selected_cache"),
                "comparison": phases.get("compare_selected_native_visitor"),
            },
            "scope": "Selected records and file identity are validated; source HEAVY receipts remain separate. Export-tool timing and runtime performance are not inferred.",
        }

    def cache_qualification(self, identity, producer, *, primary=False):
        """Check an accepted provenance bridge; never redo native authentication."""
        build = producer["raw_report"]
        if build.get("schema") in (CUMULATIVE_SCHEMA, LIGHT_SCHEMA):
            return self.cumulative_cache_qualification(
                identity, producer, primary=primary
            )
        require(
            all(
                identity[key] == build[key]
                for key in (
                    "cache_sha256",
                    "data_manifest_sha256",
                    "source_census_sha256",
                )
            ),
            "build/suite cache, data or census mismatch",
        )
        if (
            build.get("phase") != "pilot"
            and identity["workload_sha256"] == build["workload_sha256"]
            and identity["owner_proof"] == build["owner_proof"]
            and identity["build_git_revision"]
            == build["build_identity"]["build"]["git_revision"]
        ):
            return None
        path, accepted = hashed_record(
            self.index.get("cache_requalification"), self.path.parent
        )
        require(
            accepted.get("schema") == "vakint.cache-consumer-requalification.v1"
            and accepted.get("accepted") is True
            and all(
                accepted.get("checks", {}).get(check) is True
                for check in (
                    "canonical_consumer_proof",
                    "twelve_distinct_held_out_sets",
                    "training_raw_and_nonzero_owner_disjoint",
                    "cache_import_validated",
                )
            ),
            "cross-cohort cache reuse has no accepted requalification",
        )
        projected = accepted.get("qualification_projection") is not None
        require(
            accepted.get("checks", {}).get("unchanged_selection") is (not projected),
            "qualification selection difference must be recorded explicitly",
        )
        documents = {}
        for role, names in (
            (
                "producer",
                ("cache", "build_report", "workload", "owner_proof", "compiled_cohort"),
            ),
            (
                "consumer",
                (
                    "collector_artifacts",
                    "compiled_cohort",
                    "qualification_workload",
                    "final_workload",
                    "owner_proof",
                    "protected_union",
                    "training_route_qualification",
                    "alias_exclusion",
                ),
            ),
        ):
            documents[role] = {}
            for name in names:
                source = hashed_file(accepted[role][name], path.parent)
                if name != "cache":  # Native cache bytes are not a JSON document.
                    documents[role][name] = read(source)
        old, new = documents["producer"], documents["consumer"]
        require(
            accepted["producer"]["build_report"]["sha256"]
            == producer["evidence"]["sha256"]
            and old["build_report"] == build
            and accepted["producer"]["cache"]["sha256"]
            == producer["cache_file"]["sha256"]
            == identity["cache_sha256"],
            "requalification names another producer/cache",
        )
        binding = build["preparation"]["options"]["single_file"]["binding"]
        require(
            isinstance(binding, dict)
            and binding
            and accepted["semantic_binding"] == binding,
            "requalification native binding differs",
        )
        for role, doc, workload_name, expected in (
            ("producer", old, "workload", build),
            ("consumer", new, "final_workload", identity),
        ):
            workload = doc[workload_name]
            proof = doc["owner_proof"]
            revision = (
                expected["build_identity"]["build"]["git_revision"]
                if role == "producer"
                else expected["build_git_revision"]
            )
            require(
                doc["compiled_cohort"]["main_revision"] == revision,
                "requalification recorded revision differs from its own cohort",
            )
            require(
                accepted[role][workload_name]["sha256"] == expected["workload_sha256"]
                and accepted[role]["owner_proof"]["sha256"]
                == expected["owner_proof"]["sha256"]
                == workload["owner_proof"]["sha256"],
                "requalification workload/proof identity differs",
            )
            require(
                proof["compiled_cohort"]
                == doc["compiled_cohort"]
                == expected["owner_proof"]["compiled_cohort"]
                and proof["semantic_binding"]
                == doc["compiled_cohort"]["semantic_binding"]
                == binding,
                "requalification proof/cohort/binding differs",
            )
            require(
                workload["data_manifest_sha256"]
                == proof["data_manifest_sha256"]
                == expected["data_manifest_sha256"]
                and workload["source_census_sha256"]
                == expected["source_census_sha256"],
                "requalification changed frozen data/census",
            )
        selection = selected_inputs(old["workload"])
        require(
            selection
            == accepted["immutable_selection"]
            == selected_inputs(new["final_workload"]),
            "requalification changed frozen selection",
        )
        projection = None
        if projected:
            projection = self.qualification_projection(accepted, documents, path.parent)
        else:
            require(
                selection == selected_inputs(new["qualification_workload"]),
                "requalification changed frozen selection",
            )
        artifact = new["collector_artifacts"]
        require(
            all(
                artifact["stamps"].get(stamp) == new["compiled_cohort"][field]
                for field, stamp in COHORT_STAMPS
            )
            and new["compiled_cohort"]["main_revision"]
            == identity["build_git_revision"]
            and identity["executable_sha256"]
            == new["owner_proof"]["measurement_executable_sha256"]
            and identity["executable_sha256"]
            in {row["sha256"] for row in artifact_binaries(artifact)},
            "requalification consumer artifact/stamps differ",
        )
        result = {"evidence": evidence(path), "raw_receipt": accepted}
        if projection is not None:
            result["qualification_projection"] = projection
        return result

    def cumulative_cache_qualification(self, identity, producer, *, primary=False):
        """Require external native-import and consumer/exclusion acceptance."""
        source = producer.get("source_producer", producer)
        build = source["raw_report"]
        light = producer["raw_report"].get("schema") == LIGHT_SCHEMA
        require(
            cache_identity(identity) == cache_identity(producer["raw_report"]),
            "cumulative cache/data identity differs",
        )
        path, accepted = hashed_record(
            self.index.get("cache_requalification"), self.path.parent
        )
        require(
            accepted.get("schema")
            == "vakint.cumulative-cache-consumer-requalification.v1"
            and accepted.get("accepted") is True
            and all(
                accepted.get("checks", {}).get(key) is True
                for key in (
                    "canonical_consumer_proof",
                    "cache_import_validated",
                    "source_capture_inventory_authenticated",
                    "unchanged_selection",
                    "training_raw_and_nonzero_owner_disjoint",
                    "twelve_distinct_held_out_sets",
                )
            ),
            "cumulative cache lacks accepted native consumer qualification",
        )
        require(
            hashed_file(accepted["producer"], path.parent)
            == Path(producer["evidence"]["path"])
            and accepted["producer"]["sha256"] == producer["evidence"]["sha256"]
            and cache_identity(accepted) == cache_identity(producer["raw_report"])
            and accepted["semantic_binding"] == build["semantic_binding"]
            and accepted["data_manifest_sha256"] == build["data_manifest_sha256"],
            "cumulative qualification references another producer/snapshot/binding",
        )
        docs = {}
        for name in (
            "collector_artifacts",
            "compiled_cohort",
            "final_workload",
            "owner_proof",
            "protected_union",
            "training_route_qualification",
            "alias_exclusion",
        ):
            _, docs[name] = hashed_record(accepted["consumer"][name], path.parent)
        workload, proof, cohort = (
            docs[name] for name in ("final_workload", "owner_proof", "compiled_cohort")
        )
        refs = accepted["consumer"]
        require(
            refs["final_workload"]["sha256"] == identity["workload_sha256"]
            and refs["owner_proof"]["sha256"]
            == identity["owner_proof"]["sha256"]
            == workload["owner_proof"]["sha256"]
            and selected_inputs(workload) == accepted["immutable_selection"]
            and workload["source_census_sha256"] == identity["source_census_sha256"]
            and workload["data_manifest_sha256"]
            == proof["data_manifest_sha256"]
            == identity["data_manifest_sha256"],
            "cumulative qualification changed frozen workload/proof/data",
        )
        require(
            cohort
            == proof["compiled_cohort"]
            == identity["owner_proof"]["compiled_cohort"]
            and cohort["main_revision"] == identity["build_git_revision"]
            and cohort["semantic_binding"]
            == proof["semantic_binding"]
            == build["semantic_binding"]
            and all(
                docs["collector_artifacts"]["stamps"].get(stamp) == cohort[field]
                for field, stamp in COHORT_STAMPS
            )
            and (
                primary
                or identity["executable_sha256"]
                == proof["measurement_executable_sha256"]
            )
            and proof["measurement_executable_sha256"]
            in {row["sha256"] for row in artifact_binaries(docs["collector_artifacts"])}
            and identity["executable_sha256"]
            in {
                row["sha256"] for row in artifact_binaries(docs["collector_artifacts"])
            },
            "cumulative consumer proof/cohort/artifact differs",
        )
        exposure = accepted["case_exposure"]
        require(isinstance(exposure, dict), "missing explicit case exposure")
        for kind, label in (
            ("training", "seen-starting-input"),
            ("held_out", "held-out"),
        ):
            require(
                all(exposure.get(row["id"]) == label for row in workload[kind]),
                "seen starting input cannot become held-out evidence",
            )
        training = {tuple(row["powers"]) for row in workload["training"]}
        held_out = {tuple(row["powers"]) for row in workload["held_out"]}
        require(
            not training.intersection(held_out),
            "explicit starting training key is not held out",
        )
        self.cumulative_source_exclusion(accepted, docs, source, path.parent)
        if light:
            imported = None
            require(
                accepted.get("cache_import", {}).get("sha256")
                == producer["export_validation"]["native_validation"]["sha256"]
                and hashed_file(accepted["cache_import"], path.parent)
                == Path(producer["export_validation"]["native_validation"]["path"]),
                "LIGHT qualification must bind its own strict export/reimport receipt",
            )
        else:
            imported = self.cumulative_cache_import(
                accepted.get("cache_import"), producer, cohort, path.parent
            )
        metadata = None
        relation = accepted.get("data_metadata_requalification")
        require(
            not light or relation is not None,
            "LIGHT requires source-bound dataset requalification",
        )
        if relation is not None:
            metadata = self.data_metadata_requalification(
                relation, identity, source, accepted, docs, imported, path.parent
            )
            if light:
                require(
                    metadata["prepared_cache_replacement"]["export_validation"]
                    == producer["evidence"]
                    and metadata["prepared_cache_replacement"]["cache_sha256"]
                    == cache_identity(identity),
                    "LIGHT measurement uses another selected dataset export",
                )
        else:
            require(
                identity["data_manifest_sha256"] == build["data_manifest_sha256"],
                "changed data identity requires explicit metadata requalification",
            )
        result = {
            "evidence": evidence(path),
            "raw_receipt": accepted,
            "data_metadata_requalification": metadata,
        }
        result["strict_export_reimport" if light else "public_import"] = (
            producer["export_validation"] if light else imported
        )
        return result

    def data_metadata_requalification(
        self, reference, identity, producer, accepted, docs, imported, base
    ):
        """Compute the narrow release-metadata relation; preserve source identities."""

        strict_record = self.strict_record

        path, relation = strict_record(reference, base)
        require(
            relation.get("schema")
            == "vakint.external-data-metadata-requalification.v1",
            "unknown data metadata relation",
        )
        old_path, old = strict_record(relation["source_manifest"], path.parent)
        new_path, new = strict_record(relation["consumer_manifest"], path.parent)
        require(
            digest(old_path) == producer["raw_report"]["data_manifest_sha256"]
            and digest(new_path) == identity["data_manifest_sha256"],
            "metadata relation names another source/consumer dataset",
        )
        require(
            old["schema"] == new["schema"] == "vakint.external-data.v1"
            and all(
                isinstance(item["version"], str)
                and re.fullmatch(
                    r"v[0-9]+\.[0-9]+\.[0-9]+(?:[-+][A-Za-z0-9.-]+)?", item["version"]
                )
                for item in (old, new)
            ),
            "invalid release schema/version in metadata relation",
        )
        require(
            len(
                {"runtime_cache_snapshot", "prepared_cache_replacement"}.intersection(
                    relation
                )
            )
            == 1,
            "metadata relation needs exactly one cache representation",
        )
        selected = relation.get("runtime_cache_snapshot")
        light = None
        light_evidence = None
        if "prepared_cache_replacement" in relation:
            replacement = relation["prepared_cache_replacement"]
            require(
                isinstance(replacement, dict)
                and set(replacement) == {"loop_count", "export_validation"}
                and type(replacement["loop_count"]) is int
                and replacement["loop_count"] == 5
                and old.get("runtime_cache_snapshots", []) == []
                and new.get("runtime_cache_snapshots", []) == [],
                "light replacement needs only 5L and no selected native snapshot",
            )
            light_path, light = strict_record(
                replacement["export_validation"], path.parent
            )
            light_evidence = self.light_cache_export(
                light_path, light, old_path, new_path.parent, producer, strict_record
            )
        else:
            require(
                isinstance(selected, dict)
                and set(selected)
                == {"loop_count", "manifest", "manifest_sha256", "manifest_bytes"}
                and type(selected["loop_count"]) is int
                and selected["loop_count"] == 5
                and re.fullmatch(
                    r"5L/cache/[A-Za-z0-9][A-Za-z0-9._-]*\.json", selected["manifest"]
                )
                and old.get("runtime_cache_snapshots", []) == []
                and new.get("runtime_cache_snapshots") == [selected],
                "metadata relation must add exactly the selected 5L cache reference",
            )
            descriptor_path, descriptor = strict_record(
                {"path": selected["manifest"], "sha256": selected["manifest_sha256"]},
                new_path.parent,
            )
            _, packet = hashed_record(imported["packet"], path.parent)
            imported_descriptor = hashed_file(
                packet["descriptor"], Path(imported["packet"]["path"]).parent
            )
            require(
                type(selected["manifest_bytes"]) is int
                and selected["manifest_bytes"] == descriptor_path.stat().st_size
                and descriptor_path.read_bytes() == imported_descriptor.read_bytes()
                and descriptor["binding"] == producer["raw_report"]["semantic_binding"]
                and {
                    name.removeprefix("cache/"): item
                    for name, item in descriptor["files"].items()
                    if name.endswith(".rrindex")
                }
                == producer["raw_report"]["cache_native_index_identity"]["indexes"],
                "metadata cache descriptor differs from authenticated import/binding/index identity",
            )
        old_packs, new_packs = old["loop_packs"], new["loop_packs"]
        require(
            [item["loop_count"] for item in old_packs]
            == [item["loop_count"] for item in new_packs]
            == list(range(1, 6)),
            "metadata relation changed loop inventory/order",
        )
        changed = []
        verified_files = 0

        def record(pointer, before, after):
            if before != after:
                changed.append(
                    {"pointer": pointer, "source": before, "consumer": after}
                )

        def payloads(directory):
            result = {}
            for payload in directory.rglob("*"):
                require(
                    not payload.is_symlink(), "metadata payload tree contains a symlink"
                )
                require(
                    payload.is_file() or payload.is_dir(),
                    "metadata payload tree contains a special file",
                )
                if payload.is_file() and payload != directory / "manifest.json":
                    result[payload.relative_to(directory).as_posix()] = {
                        "sha256": digest(payload),
                        "bytes": payload.stat().st_size,
                    }
            return result

        for count, (before_pack, after_pack) in enumerate(zip(old_packs, new_packs), 1):
            manifest_name = f"{count}L/manifest.json"
            require(
                before_pack["manifest"] == after_pack["manifest"] == manifest_name,
                "metadata relation changed loop manifest path",
            )
            require(
                not (old_path.parent / f"{count}L").is_symlink()
                and not (new_path.parent / f"{count}L").is_symlink(),
                "metadata loop directory is a symlink",
            )
            source_path, source = strict_record(
                {"path": manifest_name, "sha256": before_pack["manifest_sha256"]},
                old_path.parent,
            )
            target_path, target = strict_record(
                {"path": manifest_name, "sha256": after_pack["manifest_sha256"]},
                new_path.parent,
            )
            require(
                source["schema"] == target["schema"] == "vakint.external-data.loop.v1"
                and source["loop_count"] == target["loop_count"] == count
                and source["version"] == old["version"]
                and target["version"] == new["version"],
                "metadata loop schema/version differs",
            )
            source_files, target_files = source["files"], target["files"]
            addition = (
                {
                    selected["manifest"].removeprefix("5L/"): {
                        "sha256": selected["manifest_sha256"],
                        "bytes": selected["manifest_bytes"],
                    }
                }
                if count == 5 and selected is not None
                else {}
            )
            replacements = {}
            if count == 5 and light is not None:
                before_cache, after_cache = source["cache"], target["cache"]
                cache_keys = {
                    "available",
                    "binding",
                    "decoded_sha256",
                    "format",
                    "identity",
                    "metadata",
                    "path",
                }
                require(
                    set(before_cache) == set(after_cache) == cache_keys
                    and before_cache["available"] is True
                    and after_cache["available"] is True
                    and before_cache["format"]
                    == "rustred.saved-application-cache-file.v1"
                    and same_json_value(
                        before_cache["binding"], light["semantic_binding"]
                    )
                    and before_cache["metadata"] == "cache.json"
                    and re.fullmatch(
                        r"cache/[A-Za-z0-9][A-Za-z0-9._-]*\.gz", before_cache["path"]
                    )
                    and same_json_value(
                        {
                            key: value
                            for key, value in before_cache.items()
                            if key not in ("identity", "decoded_sha256")
                        },
                        {
                            key: value
                            for key, value in after_cache.items()
                            if key not in ("identity", "decoded_sha256")
                        },
                    )
                    and same_json_value(
                        after_cache["identity"], light["output"]["identity"]
                    )
                    and after_cache["decoded_sha256"]
                    == light["output"]["decoded_sha256"],
                    "light replacement changed prepared-cache binding/format/paths/identity",
                )
                for kind, name in (
                    ("cache", before_cache["path"]),
                    ("metadata", before_cache["metadata"]),
                ):
                    output = light["output"][kind]
                    replacements[name] = {
                        key: output[key] for key in ("sha256", "bytes")
                    }
                    require(
                        name in source_files, "light replacement added a cache payload"
                    )
                    record(f"/5L/files/{name}", source_files[name], replacements[name])
                require(
                    same_json_value(target_files, {**source_files, **replacements}),
                    "light replacement changed a scientific file inventory",
                )
                metadata_path, metadata = strict_record(
                    {
                        "path": after_cache["metadata"],
                        **replacements[after_cache["metadata"]],
                    },
                    target_path.parent,
                )
                self.light_cache_metadata(metadata, light, new_path.parent)
                require(
                    metadata_path.stat().st_size
                    == light["output"]["metadata"]["bytes"],
                    "light cache metadata size differs",
                )
                for key in ("identity", "decoded_sha256"):
                    record(f"/5L/cache/{key}", before_cache[key], after_cache[key])
            else:
                require(
                    not set(source_files).intersection(addition)
                    and same_json_value(target_files, {**source_files, **addition}),
                    "metadata relation added/changed a scientific file inventory",
                )
            before_tree, after_tree = (
                payloads(source_path.parent),
                payloads(target_path.parent),
            )
            require(
                same_json_value(before_tree, source_files)
                and same_json_value(after_tree, target_files),
                "metadata relation changed original payload bytes or added hidden payload",
            )
            verified_files += len(before_tree) - len(replacements)
            target_compare = copy.deepcopy(target)
            target_compare.update(version=source["version"], files=source_files)
            if replacements:
                target_compare["cache"] = source["cache"]
            require(
                same_json_value(target_compare, source),
                "metadata relation changed loop scientific semantics",
            )
            for document, version in (
                (before_pack, old["version"]),
                (after_pack, new["version"]),
            ):
                require(
                    document["archive"]
                    == f"vakint-external-data-{version}-{count}L.tar.gz",
                    "invalid release archive name",
                )
                archive_sha, archive_bytes = (
                    document["archive_sha256"],
                    document["archive_bytes"],
                )
                require(
                    (archive_sha is None and archive_bytes is None)
                    or (
                        isinstance(archive_sha, str)
                        and re.fullmatch(r"[0-9a-f]{64}", archive_sha)
                        and type(archive_bytes) is int
                        and archive_bytes > 0
                    ),
                    "invalid release archive identity",
                )
            pack_compare = copy.deepcopy(after_pack)
            for key in (
                "manifest_sha256",
                "archive",
                "archive_sha256",
                "archive_bytes",
            ):
                record(
                    f"/loop_packs/{count - 1}/{key}", before_pack[key], after_pack[key]
                )
                pack_compare[key] = before_pack[key]
            require(
                same_json_value(pack_compare, before_pack),
                "metadata relation changed loop-pack semantics",
            )
            record(f"/{count}L/version", source["version"], target["version"])
            if addition:
                record(
                    "/5L/files/" + next(iter(addition)),
                    None,
                    next(iter(addition.values())),
                )
        root_compare = copy.deepcopy(new)
        root_compare.update(version=old["version"], loop_packs=old_packs)
        if selected is not None:
            if "runtime_cache_snapshots" in old:
                root_compare["runtime_cache_snapshots"] = old["runtime_cache_snapshots"]
            else:
                root_compare.pop("runtime_cache_snapshots")
        require(
            same_json_value(root_compare, old),
            "metadata relation changed other release semantics/provenance",
        )
        record("/version", old["version"], new["version"])
        if selected is not None:
            record(
                "/runtime_cache_snapshots",
                old.get("runtime_cache_snapshots"),
                new["runtime_cache_snapshots"],
            )
        _, frozen_workload = strict_record(relation["frozen_workload"], path.parent)
        frozen_proof_path, frozen_proof = strict_record(
            relation["frozen_owner_proof"], path.parent
        )
        require(
            frozen_workload["owner_proof"]["sha256"]
            == relation["frozen_owner_proof"]["sha256"],
            "metadata frozen workload/proof mismatch",
        )
        workload_compare = copy.deepcopy(docs["final_workload"])
        workload_compare.update(
            data_manifest_sha256=frozen_workload["data_manifest_sha256"],
            owner_proof=frozen_workload["owner_proof"],
        )
        require(
            same_json_value(workload_compare, frozen_workload)
            and same_json_value(
                selected_inputs(frozen_workload), accepted["immutable_selection"]
            ),
            "metadata relation transformed frozen workload/owner keys",
        )
        require(
            frozen_workload["data_manifest_sha256"]
            == frozen_proof["data_manifest_sha256"]
            == digest(old_path)
            and docs["final_workload"]["data_manifest_sha256"]
            == docs["owner_proof"]["data_manifest_sha256"]
            == digest(new_path),
            "metadata frozen/consumer proof names another dataset",
        )
        require(
            frozen_proof["loop_manifest_sha256"] == old_packs[4]["manifest_sha256"]
            and docs["owner_proof"]["loop_manifest_sha256"]
            == new_packs[4]["manifest_sha256"],
            "metadata proof names another loop manifest",
        )
        old_cohort = frozen_proof["compiled_cohort"]
        frozen_artifact_path, frozen_artifact = hashed_record(
            relation["frozen_collector_artifacts"], path.parent
        )
        require(
            all(
                frozen_artifact["stamps"].get(stamp) == old_cohort[field]
                and isinstance(old_cohort[field], str)
                and re.fullmatch(
                    r"[0-9a-fA-F]{64}"
                    if field.endswith("_sha256")
                    else r"[0-9a-fA-F]{40}",
                    old_cohort[field],
                )
                for field, stamp in COHORT_STAMPS
            )
            and frozen_proof["kernel_revision"] == old_cohort["kernel_revision"]
            and frozen_proof["rustred_revision"] == old_cohort["native_revision"]
            and frozen_proof["mapper_git_revision"] == old_cohort["main_revision"]
            and frozen_proof["measurement_executable_sha256"]
            in {item["sha256"] for item in artifact_binaries(frozen_artifact)}
            and frozen_proof["mapper_executable_sha256"]
            in {item["sha256"] for item in artifact_binaries(frozen_artifact)},
            "metadata frozen proof/artifact/cohort differs",
        )
        for item in artifact_binaries(frozen_artifact):
            if item.get("immutable_path") and Path(item["immutable_path"]).exists():
                require(
                    digest(item["immutable_path"]) == item["sha256"],
                    "metadata frozen executable changed",
                )
        for field, expected in (
            ("data_manifest_sha256", digest(old_path)),
            ("loop_manifest_sha256", old_packs[4]["manifest_sha256"]),
        ):
            if field in old_cohort:
                require(
                    old_cohort[field] == expected,
                    "metadata frozen cohort names another dataset/loop manifest",
                )
        consumer_proof_path = hashed_file(accepted["consumer"]["owner_proof"], base)
        controls = [
            row["id"]
            for group in ("tuning", "held_out")
            for row in frozen_workload[group]
        ]
        for proof, artifact, proof_path in (
            (frozen_proof, frozen_artifact, frozen_proof_path),
            (docs["owner_proof"], docs["collector_artifacts"], consumer_proof_path),
        ):
            own_cohort = proof["compiled_cohort"]
            binaries = {item["sha256"] for item in artifact_binaries(artifact)}
            require(
                proof["kernel_revision"] == own_cohort["kernel_revision"]
                and proof["rustred_revision"] == own_cohort["native_revision"]
                and proof["mapper_git_revision"] == own_cohort["main_revision"]
                and proof["mapper_executable_sha256"] in binaries
                and proof["measurement_executable_sha256"] in binaries
                and [row["id"] for row in proof["rows"]] == controls,
                "metadata proof is not bound to its own cohort/control inventory",
            )
            for row in proof["rows"]:
                _, worker = strict_record(
                    {"path": row["receipt_path"], "sha256": row["receipt_sha256"]},
                    proof_path.parent,
                )
                require(
                    worker.get("schema") == "vakint.five-loop-cache-measurement.v1"
                    and worker.get("success") is True
                    and worker.get("worker_success") is True
                    and worker.get("fresh_process") is True
                    and worker.get("variant") == "uncached"
                    and type(worker.get("sample")) is int
                    and worker["sample"] == 0
                    and worker["supervision"]["success"] is True
                    and type(worker["supervision"]["exit_code"]) is int
                    and worker["supervision"]["exit_code"] == 0
                    and worker["supervision"]["deadline_exhausted"] is False
                    and worker["case_id"] == row["id"]
                    and worker["input_sha256"] == row["input_sha256"]
                    and worker["actual_lowered_roots"] == row["actual_lowered_roots"]
                    and worker["identity"]["data_manifest_sha256"]
                    == proof["data_manifest_sha256"]
                    and worker["identity"]["build_git_revision"]
                    == own_cohort["main_revision"]
                    and worker["identity"]["executable_sha256"]
                    == proof["measurement_executable_sha256"],
                    "metadata proof worker receipt belongs to another dataset/cohort/input",
                )
        fields = (
            "id",
            "logical_cohort",
            "public_input_source_sha256",
            "input_sha256",
            "actual_lowered_roots",
            "protected_owner_roots",
            "protected_root_set_sha256",
            "application_binding",
            "complete_owner_transport",
            "reduction_succeeded",
            "mapping_error",
            "native_transports",
        )
        require(
            frozen_proof["schema"]
            == docs["owner_proof"]["schema"]
            == "vakint.public-lowering-owner-map-proof.v1"
            and frozen_proof["semantic_binding"]
            == docs["owner_proof"]["semantic_binding"]
            == producer["raw_report"]["semantic_binding"]
            and all(
                row["reduction_succeeded"] is True
                for proof in (frozen_proof, docs["owner_proof"])
                for row in proof["rows"]
            )
            and same_json_value(
                [{key: row[key] for key in fields} for row in frozen_proof["rows"]],
                [
                    {key: row[key] for key in fields}
                    for row in docs["owner_proof"]["rows"]
                ],
            ),
            "metadata relation transformed frozen native owner transports/guards",
        )
        return {
            "evidence": evidence(path),
            "source_manifest": evidence(old_path),
            "consumer_manifest": evidence(new_path),
            "verified_original_files": verified_files,
            "computed_metadata_changes": changed,
            "frozen_workload": relation["frozen_workload"],
            "frozen_owner_proof": relation["frozen_owner_proof"],
            "frozen_collector_artifacts": evidence(frozen_artifact_path),
            "prepared_cache_replacement": light_evidence,
            "scope": "Dataset metadata parity only; original source/import/physics provenance and all other consumer qualifications remain distinct.",
        }

    def light_cache_export(self, path, report, source_manifest, root, producer, record):
        """Bind an existing Native export/reimport acceptance, without decoding it."""
        require(
            report.get("schema") == "vakint.strict-light-cache-export-validation.v1"
            and report.get("accepted") is True,
            "light replacement lacks accepted native export validation",
        )
        source, selection, output, validation = (
            report[key] for key in ("source", "selection", "output", "validation")
        )
        source_producer, _ = record(source["producer"], root)
        require(
            digest(source_producer) == producer["evidence"]["sha256"]
            and source["data_manifest_sha256"] == digest(source_manifest)
            and same_json_value(
                source["cache_native_index_identity"],
                producer["raw_report"]["cache_native_index_identity"],
            )
            and same_json_value(
                report["semantic_binding"], producer["raw_report"]["semantic_binding"]
            ),
            "light export names another source producer/index/data/binding",
        )
        for ref in (source["inventory"], selection["keys"], selection["usage_profile"]):
            self.light_source_file(ref, root)
        count = selection["unique_keys"]
        require(
            type(count) is int
            and count > 0
            and selection["final_held_out_used"] is False
            and isinstance(selection["allowed_input_classes"], list)
            and len(selection["allowed_input_classes"]) == 2
            and set(selection["allowed_input_classes"])
            == {"tuning", "physical-example"},
            "light selection must use deduplicated eligible keys, not final held-out data",
        )
        identity = output["identity"]
        require(
            set(identity)
            == {
                "compressed_bytes",
                "compressed_blake3",
                "decoded_bytes",
                "decoded_blake3",
            }
            and all(
                type(identity[key]) is int and identity[key] > 0
                for key in ("compressed_bytes", "decoded_bytes")
            )
            and all(
                isinstance(identity[key], list)
                and len(identity[key]) == 32
                and all(
                    type(value) is int and 0 <= value <= 255 for value in identity[key]
                )
                for key in ("compressed_blake3", "decoded_blake3")
            )
            and re.fullmatch(r"[0-9a-f]{64}", output["decoded_sha256"]),
            "light export has an invalid native cache identity",
        )
        for kind in ("cache", "metadata"):
            item = output[kind]
            candidate = Path(item["path"])
            candidate = candidate if candidate.is_absolute() else root / candidate
            require(
                candidate.is_file()
                and not candidate.is_symlink()
                and type(item["bytes"]) is int
                and item["bytes"] > 0
                and candidate.stat().st_size == item["bytes"],
                "light export output needs matching regular file bytes",
            )
            if kind == "cache":
                require(
                    item["bytes"] <= LIGHT_CACHE_LIMIT_BYTES
                    and item["bytes"] == identity["compressed_bytes"],
                    "light cache exceeds 95,000,000 compressed bytes or differs from identity",
                )
            hashed_file(item, root)
        native_path, native = record(validation["receipt"], root)
        phases = native.get("phase_seconds", {})
        require(
            isinstance(phases, dict)
            and set(phases)
            <= {
                "mandatory_admission",
                "input_documents_and_indexes",
                "load_saved_application",
                "open_native_pack",
                "compiled_artifact_identity",
                "open_legacy_cache_file",
                "compare_full_native_visitor",
                "final_inputs_indexes_artifact_reference",
                "encode_selected_completed_cache",
                "write_selected_cache",
                "reimport_selected_cache",
                "compare_selected_native_visitor",
                "selected_cache_output_identity",
            },
            "invalid Native producer phase timings",
        )
        for seconds in phases.values():
            number(seconds)
        counts = ("selected_records", "reimported_records", "complete_records")
        zeroes = ("missing_keys", "partial_records")
        checks = (
            "exact_record_and_ordered_guard_equality",
            "declared_terminal_validation",
            "source_unchanged",
        )
        require(
            native.get("schema") == "vakint.strict-light-cache-native-validation.v1"
            and native.get("success") is True
            and all(
                type(validation[key]) is int and validation[key] == count
                for key in counts
            )
            and all(
                type(validation[key]) is int and validation[key] == 0 for key in zeroes
            )
            and all(validation[key] is True for key in checks)
            and all(
                type(native.get(key)) is int and native[key] == count for key in counts
            )
            and all(type(native.get(key)) is int and native[key] == 0 for key in zeroes)
            and all(native.get(key) is True for key in checks)
            and type(validation["codec_state_format"]) is int
            and validation["codec_state_format"] == 6
            and type(native.get("codec_state_format")) is int
            and all(
                native[key] == validation[key]
                for key in (*counts, *zeroes, *checks, "codec_state_format")
            )
            and same_json_value(native["semantic_binding"], report["semantic_binding"])
            and native["cache_sha256"] == output["cache"]["sha256"]
            and same_json_value(native["cache_identity"], identity)
            and native["decoded_sha256"] == output["decoded_sha256"]
            and native["selected_keys_sha256"] == selection["keys"]["sha256"]
            and native["source_inventory_sha256"] == source["inventory"]["sha256"]
            and same_json_value(
                native["source_cache_native_index_identity"],
                source["cache_native_index_identity"],
            )
            and native["usage_profile_sha256"] == selection["usage_profile"]["sha256"],
            "light native strict reimport/counts/guards/source proof differs",
        )
        require(
            isinstance(report.get("selection_admission"), dict),
            "light selection lacks external eligibility/frequency admission",
        )
        admission_path, admission = record(report["selection_admission"], root)
        require(
            admission.get("schema") == "vakint.light-cache-selection-admission.v1"
            and admission.get("accepted") is True
            and all(
                admission.get("checks", {}).get(key) is True
                for key in (
                    "authenticated_eligible_sources",
                    "dispatch_counts_reconciled",
                    "positive_frequency_selected_subset",
                    "final_held_out_excluded",
                )
            )
            and type(admission["selected_unique_keys"]) is int
            and admission["selected_unique_keys"] == count
            and type(admission["eligible_profile_unique_keys"]) is int
            and admission["eligible_profile_unique_keys"] >= count
            and admission["cache_metadata_sha256"] == output["metadata"]["sha256"]
            and admission["selected_keys_sha256"] == selection["keys"]["sha256"]
            and admission["usage_profile_sha256"]
            == selection["usage_profile"]["sha256"]
            and admission["source_inventory_sha256"] == source["inventory"]["sha256"]
            and admission["native_validation_sha256"] == validation["receipt"]["sha256"]
            and same_json_value(
                admission["semantic_binding"], report["semantic_binding"]
            )
            and same_json_value(
                admission["source_cache_native_index_identity"],
                source["cache_native_index_identity"],
            ),
            "light selection lacks matching external eligibility/frequency admission",
        )
        for observed in (report, native):
            _, artifact = record(observed["collector_artifacts"], root)
            validate_compiled_artifact(observed, artifact)
        return {
            "export_validation": evidence(path),
            "native_validation": evidence(native_path),
            "selection_admission": evidence(admission_path),
            "cache_sha256": output["cache"]["sha256"],
            "compressed_bytes": output["cache"]["bytes"],
            "unique_selected_records": count,
            "phase_seconds": phases,
            "scope": "Strict selected-cache dataset replacement only; no family coverage or light performance qualification.",
        }

    def light_source_file(self, reference, root, *, portable=False):
        """Hash source-only evidence; Native owns key/profile decoding, including compression."""
        name = reference["path"]
        if portable:
            require(
                isinstance(name, str)
                and name.startswith("provenance/")
                and "\\" not in name
                and ":" not in name
                and "\0" not in name
                and not Path(name).is_absolute()
                and all(part not in ("", ".", "..") for part in name.split("/")),
                "light cache metadata proof reference is not root-relative provenance",
            )
            parts = name.split("/")
            require(
                all(
                    not root.joinpath(*parts[:i]).is_symlink()
                    for i in range(1, len(parts) + 1)
                ),
                "light cache metadata proof reference contains a symlink",
            )
        candidate = Path(name)
        candidate = candidate if candidate.is_absolute() else root / candidate
        require(
            candidate.is_file()
            and not candidate.is_symlink()
            and candidate.stat().st_size < SOURCE_PROOF_LIMIT_BYTES,
            "light source proof needs a regular file smaller than 100 MiB",
        )
        return hashed_file(reference, root)

    def light_cache_metadata(self, metadata, report, root):
        """Source-only proof references use the repository root, outside loop packs."""
        selection = metadata.get("selection", {})
        require(
            set(metadata) == {"schema", "strict", "coverage_claim", "selection"}
            and metadata["schema"] == "vakint.external-data.saved-cache.v1"
            and metadata["strict"] is True
            and metadata["coverage_claim"] is False
            and set(selection)
            == {
                "kind",
                "source",
                "keys",
                "unique_keys",
                "usage_profile",
                "allowed_input_classes",
                "final_held_out_used",
                "validation",
            }
            and selection["kind"] == "usage-selected-strict-snapshot",
            "light cache metadata is not the strict usage-selected source projection",
        )

        def portable(actual, expected):
            self.light_source_file(actual, root, portable=True)
            expected_path = self.light_source_file(expected, root)
            require(
                actual["sha256"] == digest(expected_path),
                "light cache metadata proof identity differs",
            )

        source, expected_source = selection["source"], report["source"]
        require(
            set(source) == set(expected_source),
            "light cache metadata source fields differ",
        )
        for key in source:
            if key in ("producer", "inventory"):
                portable(source[key], expected_source[key])
            else:
                require(
                    same_json_value(source[key], expected_source[key]),
                    "light cache metadata source differs",
                )
        for key, value in report["selection"].items():
            if key in ("keys", "usage_profile"):
                portable(selection[key], value)
            else:
                require(
                    same_json_value(selection[key], value),
                    "light cache metadata selection differs",
                )
        portable(selection["validation"], report["validation"]["receipt"])

    def cumulative_cache_import(self, reference, producer, cohort, base):
        """Bind the actual public-import job, not a declared import check flag."""
        path, observed = hashed_record(reference, base)
        require(
            observed.get("schema")
            == "vakint.public-native-import.actual-terminal-observation.v1"
            and all(
                observed.get(key) is True
                for key in (
                    "pidfd_completion_observed",
                    "saved_launch_and_journal_identity_joined",
                    "public_import_terminal_accepted",
                    "frozen_inputs_unchanged",
                )
            )
            and observed["current_service"]["ActiveState"] == "inactive"
            and observed["current_service"]["Result"] == "success"
            and observed["current_service"]["ExecMainStatus"] == "0",
            "public cache import has no successful actual process completion",
        )
        pins = observed["artifact_sha256"]

        def pinned(relative):
            target = (path.parent / relative).resolve()
            require(
                str(target) in pins, "public import observation omits pinned evidence"
            )
            return hashed_record(
                {"path": str(target), "sha256": pins[str(target)]}, path.parent
            )

        packet_path, packet = pinned("packet.json")
        go_path, go = pinned("root-go.json")
        _, launch = pinned("actual-launch-receipt.json")
        _, journal = pinned("actual-invocation-journal.json")
        terminal_path, terminal = pinned("run1/terminal-outcome.json")
        worker_path, worker = hashed_record(
            terminal["worker_result"], terminal_path.parent
        )
        require(
            pins.get(str(worker_path)) == terminal["worker_result"]["sha256"],
            "public import observation omits its worker receipt",
        )
        require(
            terminal.get("schema") == "vakint.public-native-import-controller.v1"
            and terminal.get("state") == "PASS"
            and terminal.get("worker_returncode") == 0
            and terminal.get("error") is None
            and terminal["packet_sha256"] == digest(packet_path)
            and terminal["root_go_sha256"] == digest(go_path)
            and terminal == observed["controller_terminal"]
            and number(terminal["whole_elapsed_seconds"])
            < number(terminal["whole_seconds"]),
            "public import controller/worker/packet completion differs",
        )
        require(
            packet.get("schema")
            == "vakint.public-native-large-prefix-import-source-packet.v1"
            and go["packet_sha256"] == digest(packet_path)
            and go["decision"] == "GO_ONCE_EXACT_PUBLIC_IMPORT"
            and terminal["whole_seconds"]
            == packet["whole_seconds"]
            == go["whole_seconds"]
            and go["unit"] == packet["unit"] == observed["unit"] == launch["unit"]
            and launch.get("schema")
            == "vakint.public-large-native-import.actual-launch.v1"
            and launch["service"]["InvocationID"] == observed["invocation"]
            and observed["current_service"].get("InvocationID", "")
            in ("", observed["invocation"])
            and journal
            and all(
                row.get("USER_UNIT") == observed["unit"]
                and row.get("USER_INVOCATION_ID") == observed["invocation"]
                for row in journal
            ),
            "public import actual launch/journal identity differs",
        )
        # These are the actual admission and child invocation, not a second
        # eligibility decision. Native bytes remain authenticated by the importer.
        launch_records = {}
        for name in ("admission.json", "invocation-start.json", "actual-host.json"):
            target = (path.parent / "run1" / name).resolve()
            require(
                str(target) in launch["artifact_sha256"],
                "public import launch omits process/host evidence",
            )
            _, launch_records[name] = hashed_record(
                {"path": str(target), "sha256": launch["artifact_sha256"][str(target)]},
                path.parent,
            )
        admission, invocation = (
            launch_records["admission.json"],
            launch_records["invocation-start.json"],
        )
        require(
            admission == launch["admission"]
            and admission["packet_sha256"] == digest(packet_path)
            and admission.get("frozen_inputs_and_host_hashes_verified") is True
            and invocation["packet_sha256"] == digest(packet_path)
            and invocation["command"] == packet["command"]
            and invocation["controller_pid"] == int(launch["service"]["MainPID"])
            and all(
                any(
                    row["pid"] == pid
                    and type(row.get("start_ticks")) is int
                    and row["start_ticks"] > 0
                    for row in launch["processes"]
                )
                for pid in (invocation["pid"], invocation["controller_pid"])
            ),
            "public import actual admitted process differs",
        )
        descriptor_path, descriptor = hashed_record(
            packet["descriptor"], packet_path.parent
        )
        build = producer["raw_report"]
        require(
            descriptor.get("schema") == "vakint.nativepack-transport.v1"
            and descriptor["loop_count"] == 5
            and descriptor["binding"] == build["semantic_binding"]
            and {
                name.removeprefix("cache/"): item
                for name, item in descriptor["files"].items()
                if name.endswith(".rrindex")
            }
            == build["cache_native_index_identity"]["indexes"],
            "public import descriptor names another native snapshot/binding",
        )
        files = descriptor["files"]
        chunks = {
            name: item for name, item in files.items() if name.endswith(".rrchunk.gz")
        }
        indexes = {
            name: item for name, item in files.items() if name.endswith(".rrindex")
        }
        namespaces = sum(name.startswith("cache/bindings/") for name in indexes)
        expected = {
            "format": "multipart_manifest",
            "output": packet["fresh_destination"],
            "sha256": digest(descriptor_path),
            "archive_bytes": descriptor_path.stat().st_size
            + sum(part["archive_bytes"] for part in descriptor["parts"]),
            "extracted_bytes": sum(item["bytes"] for item in files.values()),
            "entries": len(files),
            "namespaces": namespaces,
            "copied_bases": namespaces,
            "copied_chunks": len(chunks),
            "preserved_record_references": descriptor["stored_record_references"],
            "compressed_bytes": sum(item["bytes"] for item in chunks.values()),
            "index_bytes": sum(item["bytes"] for item in indexes.values()),
            "transfer_buffer_bytes": 65536,
        }
        require(
            worker.get("schema") == "vakint.public-large-native-import.v1"
            and worker.get("state") == "PASS"
            and worker["import_report"] == expected
            and worker["exact_native_files_verified"] == len(files)
            and worker.get("refused_repeat_preserved_full_output_and_cache_state")
            is True
            and worker.get("cache_attachment_performed") is False
            and worker.get("integral_reduction_or_numerical_master_evaluation") is False
            and worker.get("cache_coefficient_lookup_or_capture") is False
            and worker["repeat_refusal"]["type"] == "ValueError"
            and "destination already exists" in worker["repeat_refusal"]["message"]
            and terminal["output"] == expected["output"]
            and native_cache_identity(expected["output"])
            == build["cache_native_index_identity"],
            "actual public import report/inventory or scope differs",
        )
        require(
            expected["preserved_record_references"]
            == build["inventory"]["stored_record_references"]
            and expected["compressed_bytes"] == build["inventory"]["compressed_bytes"]
            and expected["copied_chunks"] == build["inventory"]["physical_chunks"],
            "public import inventory differs from cumulative snapshot",
        )
        sources = packet["public_sources"]
        host_path = Path(worker["actual_host"]).resolve()
        require(
            str(host_path) in launch["artifact_sha256"],
            "public import launch omits actual host artifact",
        )
        host = read(
            hashed_file(
                {
                    "path": str(host_path),
                    "sha256": launch["artifact_sha256"][str(host_path)],
                },
                path.parent,
            )
        )
        require(
            worker["public_sources"]
            == sources
            == host["sources"]
            == launch["actual_host"]["sources"]
            and sources["rustred"] == cohort["native_revision"]
            and sources["official_kernel"] == cohort["kernel_revision"]
            and host == launch["actual_host"] == launch_records["actual-host.json"]
            and host["core_sha256"] == packet["core"]["sha256"]
            and host["wheel_sha256"] == packet["wheel"]["sha256"]
            and host["mapped_core_paths"] == [packet["core"]["path"]],
            "public import host artifact/native/kernel cohort differs",
        )
        manifest = Path(packet["frozen_data"]) / "manifest.json"
        require(
            digest(manifest) == build["data_manifest_sha256"],
            "public import used another data manifest",
        )
        loaded = worker["data_report"]
        state = worker["final_attachment_state"]
        require(
            read(worker_path.parent / "public-import-report.json")
            == worker["import_report"]
            and read(worker_path.parent / "data-load-report.json") == loaded
            and read(worker_path.parent / "cache-state-before-repeat.json") == state
            and Path(loaded["path"]).resolve() == Path(packet["frozen_data"]).resolve()
            and loaded["version"] == read(manifest)["version"]
            and loaded["load_prepared_caches"] is False
            and {row["loop_count"] for row in loaded["loops"]} == set(range(1, 6))
            and all(row["cache_loaded"] is False for row in loaded["loops"])
            and set(state) == {str(loop) for loop in range(1, 6)}
            and all(row["attached"] is False for row in state.values())
            and state["5"]["saved_application"]["binding"] == build["semantic_binding"]
            and state["5"]["saved_application"]["read_backing"] is False
            and state["5"]["saved_application"]["native"]["rule_applications"] == 0
            and state["5"]["saved_application"]["persistent"]["lookups"]
            == state["5"]["saved_application"]["persistent"]["captured_records"]
            == 0,
            "public import data/binding/attachment state differs",
        )
        return {
            "observation": evidence(path),
            "terminal": evidence(terminal_path),
            "worker": evidence(worker_path),
            "packet": evidence(packet_path),
            "raw_import_report": worker["import_report"],
            "importer_sources": sources,
            "scope": "Public import and native-byte preservation; attachment, physics and performance remain separate.",
        }

    def cumulative_source_exclusion(self, accepted, docs, producer, base):
        """Check all explicit source roots against native route/protected facts."""
        inventory_path, inventory = hashed_record(
            accepted["source_starting_keys"], base
        )
        require(
            inventory.get("schema") == "vakint.cumulative-source-starting-keys.v1",
            "missing cumulative starting-key inventory",
        )
        expected = {row["id"]: row for row in producer["raw_report"]["sources"]}
        rows = inventory["sources"]
        require(
            len(rows) == len(expected) and {row["id"] for row in rows} == set(expected),
            "incomplete cumulative starting-key inventory",
        )
        starting = set()
        for row in rows:
            source = expected[row["id"]]
            require(
                row["producer"]["sha256"] == source["producer"]["sha256"],
                "starting keys belong to another producer",
            )
            origin_path, origin = hashed_record(row["producer"], inventory_path.parent)
            input_path, inputs = hashed_record(row["inputs"], inventory_path.parent)
            if row["kind"] == "selected-workload":
                require(
                    inputs.get("schema") == "vakint.five-loop-cache-workload.v1"
                    and origin.get("workload_sha256") == digest(input_path),
                    "starting workload differs from actual producer",
                )
                starting.update(tuple(item["powers"]) for item in inputs["training"])
            elif row["kind"] == "captured-tuning-workload":
                require(
                    inputs.get("schema") == "vakint.five-loop-cache-workload.v1"
                    and isinstance(inputs.get("tuning"), list)
                    and all(isinstance(item, dict) for item in inputs["tuning"]),
                    "invalid captured tuning workload",
                )
                cases = [
                    item for item in inputs["tuning"] if item["id"] == row["case_id"]
                ]
                require(len(cases) == 1, "captured tuning case is absent or ambiguous")
                case = cases[0]
                require(
                    isinstance(origin.get("identity"), dict),
                    "captured tuning measurement lacks its original identity",
                )
                identity = origin["identity"]
                require(
                    origin.get("schema") == "vakint.five-loop-cache-measurement.v1"
                    and origin.get("success") is True
                    and (
                        "worker_success" not in origin
                        or origin["worker_success"] is True
                    )
                    and (
                        "supervision" not in origin
                        or isinstance(origin["supervision"], dict)
                        and origin["supervision"].get("success") is True
                        and origin["supervision"].get("deadline_exhausted") is False
                        and type(origin["supervision"].get("exit_code")) is int
                        and origin["supervision"]["exit_code"] == 0
                    )
                    and all(
                        origin.get(key) is None
                        for key in ("error", "failure", "resource_failure")
                    )
                    and origin.get("provisional") is True
                    and origin.get("native_keys_authenticated") is False
                    and identity.get("provisional") is True
                    and identity["suite"] == "tuning"
                    and identity["only_case"] == origin["case_id"] == row["case_id"]
                    and identity["workload_sha256"] == digest(input_path)
                    and identity["source_census_sha256"]
                    == inputs["source_census_sha256"]
                    and isinstance(identity.get("data_manifest_sha256"), str)
                    and re.fullmatch(r"[0-9a-f]{64}", identity["data_manifest_sha256"])
                    and identity["data_manifest_sha256"]
                    == inputs["data_manifest_sha256"]
                    and same_json_value(
                        origin["cache_statistics"]["binding"],
                        producer["raw_report"]["semantic_binding"],
                    ),
                    "captured tuning measurement differs from its legacy workload/case/binding",
                )
                proof = docs["owner_proof"]
                require(
                    proof.get("schema") == "vakint.public-lowering-owner-map-proof.v1"
                    and same_json_value(
                        proof["compiled_cohort"], docs["compiled_cohort"]
                    )
                    and same_json_value(
                        proof["semantic_binding"],
                        producer["raw_report"]["semantic_binding"],
                    ),
                    "captured tuning owner proof/cohort differs",
                )
                matches = [
                    item for item in proof["rows"] if item["id"] == row["case_id"]
                ]
                require(len(matches) == 1, "captured tuning owner proof omits its case")
                matched = matches[0]
                require(
                    isinstance(case.get("public_input"), str)
                    and matched["logical_cohort"] == "tuning"
                    and matched["public_input_source_sha256"]
                    == hashlib.sha256(case["public_input"].encode()).hexdigest()
                    and matched["input_sha256"] == origin["input_sha256"]
                    and matched.get("complete_owner_transport") is True
                    and matched.get("reduction_succeeded") is True
                    and matched.get("mapping_error") is None
                    and same_json_value(
                        matched["application_binding"],
                        producer["raw_report"]["semantic_binding"],
                    ),
                    "captured tuning input or complete owner transport differs",
                )
                roots = case["lowered_roots"]
                require(
                    isinstance(roots, list)
                    and roots
                    and all(
                        isinstance(key, list)
                        and len(key) == 15
                        and all(type(power) is int for power in key)
                        for key in roots
                    )
                    and len({tuple(key) for key in roots}) == len(roots)
                    and same_json_value(roots, origin["actual_lowered_roots"])
                    and same_json_value(roots, matched["actual_lowered_roots"]),
                    "captured tuning raw starting keys differ",
                )
                # Historical authentication flags remain false. This current
                # proof binds the raw input; postcombine seeds cannot replace it.
                starting.update(tuple(key) for key in roots)
            else:
                require(
                    row["kind"] == "physical-census"
                    and inputs.get("schema") == "vakint.qcd-ladder.campaign-census.v1",
                    "unknown starting-key source",
                )
                config, config_path = origin, origin_path
                if source["scope"] == "completed-example":
                    config_path, config = hashed_record(source["config"], base)
                prepared = config["prepared_input"]
                require(
                    config.get("schema")
                    == "vakint.inactive-capture-physical-user-unit.v1"
                    and inputs["source_sha256"] == prepared["input"]["sha256"]
                    and inputs["projection"] == prepared["manifest_projection"]
                    and type(inputs["source_bytes"]) is int
                    and inputs["source_bytes"] == prepared["input"]["bytes"]
                    and inputs["family_fingerprint"]
                    == producer["raw_report"]["semantic_binding"]["family_fingerprint"],
                    "physical census differs from producer input",
                )
                if "path" in prepared["input"]:
                    scalar = hashed_file(prepared["input"], config_path.parent)
                else:
                    manifest_path, manifest = hashed_record(
                        prepared["manifest"], config_path.parent
                    )
                    require(
                        manifest.get("schema") == "vakint.qcd-ladder.prepared-input.v1"
                        and manifest["definition"]
                        == "massless-qcd-four-rung-native-uv-v1"
                        and manifest["projection"] == prepared["manifest_projection"]
                        and same_json_value(manifest["input"], prepared["input"])
                        and Path(prepared["directory"]).is_absolute()
                        and Path(prepared["directory"]).resolve()
                        == manifest_path.parent,
                        "physical prepared manifest identity differs",
                    )
                    # Reuse the prepared format's established payload name,
                    # without adding a path to the immutable configuration.
                    scalar = manifest_path.parent / "vakint-input.symbolica"
                require(
                    type(prepared["input"]["bytes"]) is int
                    and prepared["input"]["bytes"] > 0
                    and scalar.stat().st_size == prepared["input"]["bytes"]
                    and digest(scalar) == inputs["source_sha256"],
                    "physical starting input changed",
                )
                starting.update(
                    tuple(item["powers"]) for item in inputs["combined_roots"]
                )
        require(
            starting
            and all(
                len(key) == 15 and all(type(power) is int for power in key)
                for key in starting
            ),
            "invalid or empty source starting keys",
        )
        routes, protected, aliases = (
            docs[name]
            for name in (
                "training_route_qualification",
                "protected_union",
                "alias_exclusion",
            )
        )
        require(
            routes.get("schema") == "vakint.source-candidate-owner-qualification.v2"
            and routes.get("all_routes_succeeded") is True
            and routes.get("selected_rows_complete") is True
            and routes["compiled_cohort"] == docs["compiled_cohort"]
            and routes["semantic_binding"]
            == producer["raw_report"]["semantic_binding"],
            "cumulative native route proof/cohort differs",
        )
        route_rows = routes["rows"]
        require(
            len(route_rows) == len(starting)
            and {tuple(row["raw_powers"]) for row in route_rows} == starting,
            "native route proof omits source starting keys",
        )
        require(
            protected.get("schema") == "vakint.protected-raw-owner-union.v1"
            and protected["owner_proof"]["sha256"]
            == accepted["consumer"]["owner_proof"]["sha256"]
            and protected["compiled_cohort"]["sha256"]
            == accepted["consumer"]["compiled_cohort"]["sha256"],
            "protected owner union belongs to another proof/cohort",
        )
        forbidden = {
            tuple(key) for key in protected["forbidden_direct_target_or_owner_keys"]
        }
        require(
            forbidden
            == {
                tuple(key)
                for name in (
                    "raw_lowered_keys",
                    "direct_owner_keys",
                    "declared_case_powers",
                )
                for key in protected[name]
            },
            "protected union omits declared/raw/owner keys",
        )
        required = set()
        for row in docs["final_workload"]["held_out"]:
            required.add(tuple(row["powers"]))
            required.update(
                tuple(key)
                for name in ("protected_roots", "lowered_roots")
                for key in row[name]
            )
        require(
            required.issubset(forbidden),
            "protected union omits frozen held-out input/owner keys",
        )
        require(
            aliases.get("schema")
            == "vakint.source-candidate-external-protected-alias-exclusion.v1"
            and aliases["qualification"]["sha256"]
            == accepted["consumer"]["training_route_qualification"]["sha256"]
            and aliases["protected_union"]["sha256"]
            == accepted["consumer"]["protected_union"]["sha256"],
            "alias exclusion belongs to another route/protected proof",
        )
        require(
            len(aliases["rows"]) == len(starting), "alias proof omits starting keys"
        )
        alias_rows = {tuple(row["powers"]): row for row in aliases["rows"]}
        for row in route_rows:
            raw = tuple(row["raw_powers"])
            owners = {
                tuple(endpoint["powers"])
                for endpoint in row["owner_endpoints"]
                if endpoint["coefficient"] != "0"
            }
            require(
                row.get("route_success") is True
                and row.get("route_error") is None
                and owners
                and raw not in forbidden
                and not owners.intersection(forbidden),
                "source starting/raw owner key overlaps held-out protection",
            )
            alias = alias_rows.get(raw, {})
            require(
                alias.get("route_complete_nonempty") is True
                and alias.get("raw_protected_intersection") is False
                and alias.get("owner_protected_intersections") == []
                and alias.get("eligible_after_external_alias_exclusion") is True,
                "native alias exclusion rejects source key",
            )

    def qualification_projection(self, accepted, documents, base):
        """Authenticate unused training-prefix metadata, never project timings."""
        path, relation = hashed_record(accepted["qualification_projection"], base)
        require(
            relation.get("schema") == "vakint.qualification-workload-projection.v1"
            and relation.get("relation") == "unexecuted-training-prefix-pruning",
            "unknown qualification projection",
        )
        detail_path, detail = hashed_record(relation.get("evidence"), path.parent)
        review_path, review = hashed_record(
            relation.get("independent_review"), path.parent
        )
        scope_path, scope = hashed_record(relation.get("validator_scope"), path.parent)
        failed_path, failed = hashed_record(
            relation.get("incomplete_attempt"), path.parent
        )
        packet_path, packet = hashed_record(
            review.get("evidence_packet"), review_path.parent
        )
        require(
            review.get("schema")
            == "vakint.independent-fallback2-semantic-projection-review.v1"
            and review.get("verdict")
            == "SOUND_FOR_NATIVE_COMPATIBILITY_WITH_EXPLICIT_AUTHENTICATED_PROJECTION_NOT_A_NEW_F022_SUITE_PASS"
            and scope.get("schema")
            == "vakint.independent-projection-validator-scope-addendum.v1"
            and scope["review"]["sha256"] == digest(review_path)
            and packet["exact_semantic_projection_evidence"]["sha256"]
            == digest(detail_path)
            and review["new_attempt_preserved"]["postflight"]["sha256"]
            == digest(failed_path),
            "projection lacks its reviewed source/validator scope",
        )
        _, reviewed = hashed_record(
            packet["native_compatibility_draft"], packet_path.parent
        )
        require(
            hashed_file(scope["review"], scope_path.parent) == review_path,
            "validator scope references another independent review",
        )
        for role in ("producer", "consumer"):
            for key, reference in accepted[role].items():
                require(
                    reviewed[role][key]["sha256"] == reference["sha256"],
                    "projection review belongs to another cohort or role",
                )
        require(
            detail.get("schema")
            == "vakint.qualification-workload-projection-evidence.v1"
            and detail["source_fresh_qualification_workload"]["sha256"]
            == accepted["consumer"]["qualification_workload"]["sha256"]
            and detail["compatibility_target_workload"]["sha256"]
            == accepted["consumer"]["final_workload"]["sha256"]
            and detail["producer_workload"]["sha256"]
            == accepted["producer"]["workload"]["sha256"]
            and detail["canonical_proof"]["sha256"]
            == accepted["consumer"]["owner_proof"]["sha256"]
            and detail["compiled_cohort"]["sha256"]
            == accepted["consumer"]["compiled_cohort"]["sha256"]
            and type(detail["fresh_training_reductions_executed"]) is int
            and detail["fresh_training_reductions_executed"] == 0,
            "projection source/target identities differ",
        )
        old, new = documents["producer"], documents["consumer"]
        source, target = new["qualification_workload"], new["final_workload"]
        require(
            len(source["training"]) == 3
            and len(target["training"]) == 2
            and target["training"]
            == old["workload"]["training"]
            == source["training"][:2]
            and source["data_manifest_sha256"] == target["data_manifest_sha256"]
            and source["source_census_sha256"] == target["source_census_sha256"]
            and source.get("owner_proof") is None,
            "projection must retain the original two-root training prefix",
        )
        expected = []
        for role, count in (("tuning", 1), ("held_out", 12)):
            require(
                len(source[role]) == len(target[role]) == count,
                "projection omits public cases",
            )
            for original, current in zip(source[role], target[role]):
                require(
                    all(
                        original[key] == current[key]
                        for key in ("id", "powers", "public_input", "rationale")
                    ),
                    "projection changed public case identity, role or order",
                )
                expected.append((role, current))
        proof, cohort = new["owner_proof"], new["compiled_cohort"]
        rows = {row["id"]: row for row in proof["rows"]}
        details = detail["full_case_key_value_and_guard_evidence"]
        require(
            len(rows) == len(proof["rows"]) == 13
            and set(rows) == {row["id"] for _, row in expected}
            and [row["id"] for row in details] == [row["id"] for _, row in expected],
            "projection proof omits, duplicates or reorders public cases",
        )
        artifact = new["collector_artifacts"]
        snapshot_path, snapshot = hashed_record(artifact["main_source_snapshot"], base)
        require(
            digest(snapshot_path) == cohort["main_source_snapshot_sha256"]
            and snapshot["head"] == cohort["main_revision"],
            "projection compiled source snapshot differs",
        )
        source_hashes = {}
        for reference in review["source_files"]:
            source_file = hashed_file(reference, review_path.parent)
            source_hashes[str(source_file)] = reference["sha256"]

        def key_set(values):
            return {tuple(value) for value in values}

        forbidden, signatures, complete = set(), set(), {}
        for (role, root), checked in zip(expected, details):
            row = rows[root["id"]]
            worker_path, worker = hashed_record(
                checked["actual_qualification_receipt"], detail_path.parent
            )
            _, historical = hashed_record(
                checked["original_producer_receipt"], detail_path.parent
            )
            supervision, identity = worker["supervision"], worker["identity"]
            require(
                worker["success"] is True
                and worker["worker_success"] is True
                and worker["fresh_process"] is True
                and worker["variant"] == "uncached"
                and type(worker["sample"]) is int
                and worker["sample"] == 0
                and supervision["success"] is True
                and type(supervision["exit_code"]) is int
                and supervision["exit_code"] == 0
                and supervision["deadline_exhausted"] is False
                and not any(worker["persistent_cache_attachments_before"].values())
                and not any(
                    worker["persistent_cache_attachments_during_evaluation"].values()
                ),
                "projection requires the original complete uncached public success",
            )
            require(
                worker["case_id"] == root["id"] == checked["id"]
                and checked["role"] == role
                and worker["declared_common_powers"] == root["powers"]
                and row["logical_cohort"] == role.replace("_", "-")
                and row["public_input_source_sha256"]
                == hashlib.sha256(root["public_input"].encode()).hexdigest()
                and row["receipt_sha256"] == digest(worker_path)
                and identity["workload_sha256"]
                == row["workload_sha256"]
                == accepted["consumer"]["qualification_workload"]["sha256"]
                and identity["data_manifest_sha256"] == target["data_manifest_sha256"]
                and identity["source_census_sha256"] == target["source_census_sha256"]
                and identity["build_git_revision"] == cohort["main_revision"]
                and identity["executable_sha256"]
                == proof["measurement_executable_sha256"]
                and identity["cache_sha256"] is None
                and worker["cache_statistics"]["binding"]
                == row["application_binding"]
                == accepted["semantic_binding"],
                "projection source worker belongs to another case/data/cohort",
            )
            for field, filename in (
                ("compiled_driver_sha256", "measure_cache.rs"),
                ("compiled_workload_schema_sha256", "build_cache/workload.rs"),
                ("compiled_owner_proof_sha256", "build_cache/owner_proof.rs"),
            ):
                name = "crates/vakint/examples/" + filename
                matches = [
                    sha for p, sha in source_hashes.items() if p.endswith("/" + name)
                ]
                require(
                    matches == [snapshot["source_sha256"][name]] == [identity[field]],
                    "projection uses a different reviewed worker",
                )
            require(
                row["reduction_succeeded"] is True
                and row["complete_owner_transport"] is True
                and row["mapping_error"] is None
                and key_set(worker["actual_lowered_roots"])
                == key_set(row["actual_lowered_roots"])
                == key_set(root["lowered_roots"])
                and key_set(row["protected_owner_roots"])
                == key_set(root["protected_roots"])
                and root["lowered_roots"]
                and root["protected_roots"],
                "projection changed the authenticated raw/owner sets",
            )
            for filename, field in (
                ("value.symbolica", "value_sha256"),
                ("guards.json", "guards_sha256"),
            ):
                require(
                    worker[field] == historical[field],
                    "projection changed an original value or complete guard set",
                )
                hashed_file(
                    {
                        "path": str(worker_path.with_name(filename)),
                        "sha256": worker[field],
                    },
                    base,
                )
                original_path = hashed_file(
                    checked["original_producer_receipt"], detail_path.parent
                )
                hashed_file(
                    {
                        "path": str(original_path.with_name(filename)),
                        "sha256": historical[field],
                    },
                    base,
                )
            closure = worker["catalogue_closure"]
            unknown = hashed_file(
                {
                    "path": str(worker_path.with_name("unknown-masters.json")),
                    "sha256": closure["unknown_sha256"],
                },
                base,
            )
            require(
                closure["closed"] is True
                and closure["unknown_count"] == 0
                and read(unknown) == [],
                "projection source closure is unavailable",
            )
            forbidden.update(
                key_set(root["lowered_roots"])
                | key_set(root["protected_roots"])
                | {tuple(root["powers"])}
            )
            if role == "held_out":
                signatures.add(tuple(sorted(key_set(root["protected_roots"]))))
            complete[root["id"]] = worker
        union = new["protected_union"]
        require(
            len(signatures) == 12
            and len(forbidden) == 51
            and forbidden == key_set(union["forbidden_direct_target_or_owner_keys"])
            and union["owner_proof"]["sha256"]
            == accepted["consumer"]["owner_proof"]["sha256"]
            and union["compiled_cohort"]["sha256"]
            == accepted["consumer"]["compiled_cohort"]["sha256"],
            "projection changed the complete protected union",
        )
        routes, aliases = new["training_route_qualification"], new["alias_exclusion"]
        require(
            routes["semantic_binding"] == accepted["semantic_binding"]
            and routes["compiled_cohort"] == cohort
            and routes["source_census_sha256"] == target["source_census_sha256"]
            and aliases["qualification"]["sha256"]
            == accepted["consumer"]["training_route_qualification"]["sha256"]
            and aliases["protected_union"]["sha256"]
            == accepted["consumer"]["protected_union"]["sha256"]
            and aliases["canonical_cohort_sha256"]
            == accepted["consumer"]["compiled_cohort"]["sha256"],
            "projection routes belong to another runtime or protected union",
        )
        for root in target["training"]:
            route = [r for r in routes["rows"] if r["raw_powers"] == root["powers"]]
            alias = [r for r in aliases["rows"] if r["powers"] == root["powers"]]
            require(
                len(route) == len(alias) == 1
                and route[0]["route_success"] is True
                and route[0]["route_error"] is None,
                "projection lacks complete selected-root routing",
            )
            endpoints = route[0]["owner_endpoints"]
            require(
                endpoints
                and all(
                    isinstance(e["coefficient"], str)
                    and e["coefficient"].strip()
                    and e["coefficient"].strip() != "0"
                    for e in endpoints
                ),
                "projection lost a nonzero direct owner",
            )
            require(
                tuple(root["powers"]) not in forbidden
                and key_set([e["powers"] for e in endpoints]).isdisjoint(forbidden)
                and alias[0]["eligible_after_external_alias_exclusion"] is True
                and alias[0]["raw_protected_intersection"] is False
                and alias[0]["owner_protected_intersections"] == [],
                "projection training leaks a protected raw or owner root",
            )
        _, validated = hashed_record(
            detail["actual_existing_validator_two_root_result"], detail_path.parent
        )
        require(
            validated["valid"] is True
            and validated["reduction_or_cache_import_performed"] is False
            and validated["data_manifest_sha256"] == target["data_manifest_sha256"]
            and validated["public_cases_checked"] == 13
            and validated["workload_sha256"]
            == accepted["consumer"]["final_workload"]["sha256"]
            and validated["identity"]["sha256"]
            == accepted["consumer"]["owner_proof"]["sha256"]
            and validated["identity"]["compiled_cohort"] == cohort,
            "projection lacks actual current-target OwnerProof validation",
        )
        import_path, imported = hashed_record(
            detail["strict_import_equality_validation"], detail_path.parent
        )
        cached_path, cached = hashed_record(
            detail["two_root_cached_public_case"], detail_path.parent
        )
        _, observed = hashed_record(imported["supervision"], import_path.parent)
        require(
            imported["actual_cached_worker"]["sha256"] == digest(cached_path)
            and all(
                imported["checks"].get(k) is True
                for k in (
                    "compressed_and_decoded_identity",
                    "single_member_and_framing",
                    "native_binding_and_every_record_validation",
                    "fresh_cached_case_exact_value_and_all_guards",
                    "catalogue_closure",
                )
            )
            and type(observed["exit_code"]) is int
            and observed["exit_code"] == 0
            and observed["deadline_exhausted"] is False
            and cached["success"] is True
            and cached["native_keys_authenticated"] is True
            and cached["variant"] == "cached"
            and cached["identity"]["workload_sha256"]
            == accepted["consumer"]["final_workload"]["sha256"]
            and cached["identity"]["build_git_revision"] == cohort["main_revision"]
            and cached["identity"]["data_manifest_sha256"]
            == target["data_manifest_sha256"]
            and cached["identity"]["source_census_sha256"]
            == target["source_census_sha256"]
            and cached["owner_proof"]["sha256"]
            == accepted["consumer"]["owner_proof"]["sha256"]
            and cached["owner_proof"]["compiled_cohort"] == cohort
            and set(cached["persistent_cache_attachments_before"])
            == {str(n) for n in range(1, 6)}
            and all(
                v is False
                for v in cached["persistent_cache_attachments_before"].values()
            )
            and cached["persistent_cache_attachments_during_evaluation"]
            == {str(n): n == 5 for n in range(1, 6)}
            and cached["identity"]["executable_sha256"]
            == proof["measurement_executable_sha256"]
            and cached["identity"]["cache_sha256"]
            == accepted["producer"]["cache"]["sha256"]
            and cached["cache_statistics"]["binding"] == accepted["semantic_binding"]
            and cached["cache_statistics"]["file"]["records"]
            == imported["actual_file_statistics"]["records"]
            == imported["producer_records"]
            > 0,
            "projection lacks actual strict native import/evaluation",
        )
        baseline = complete.get(cached["case_id"])
        require(
            baseline is not None,
            "cached import case is outside the unchanged public suite",
        )
        for filename, field in (
            ("value.symbolica", "value_sha256"),
            ("guards.json", "guards_sha256"),
        ):
            require(
                cached[field] == baseline[field],
                "cached import changed full value/guards",
            )
            hashed_file(
                {"path": str(cached_path.with_name(filename)), "sha256": cached[field]},
                base,
            )
        require(
            failed["complete"] is False
            and failed["accepted"] is False
            and all(type(value) is int for value in failed["counts"].values())
            and failed["counts"]
            == {
                "functional_success": 10,
                "deadline_failed": 1,
                "not_started_after_failure": 2,
                "started_native_workers": 11,
                "planned_public_cases": 13,
            }
            and failed["failed_case"]["exit_code"] == 124
            and failed["failed_case"]["worker_json_exists"] is False
            and failed["mapper_launched"] is False,
            "the incomplete target attempt must not be promoted or omitted",
        )
        return {
            "evidence": evidence(path),
            "raw_relation": relation,
            "scope": "Native compatibility only; original13 source identities retained, no fresh target-suite or timing claim.",
        }

    def final_selection(self, producer):
        """Keep a selected pilot's original build separate from the later decision."""
        path, selection = hashed_record(
            self.index.get("final_selection"), self.path.parent
        )
        require(
            selection.get("schema") == "vakint.cache-final-selection.v1"
            and selection.get("accepted") is True
            and selection.get("decision") == "retain_complete_candidate"
            and selection.get("producer_phase") == "pilot"
            and isinstance(selection.get("reason"), str)
            and selection["reason"].strip(),
            "pilot cache has no accepted final selection",
        )
        require(
            isinstance(selection.get("decided_utc"), str)
            and datetime.fromisoformat(selection["decided_utc"]).utcoffset()
            == timedelta(0),
            "final selection needs its actual UTC decision time",
        )
        refs = {
            key: evidence(hashed_file(selection.get(key), path.parent))
            for key in (
                "cache",
                "build_report",
                "producer_workload",
                "consumer_workload",
                "candidate_acceptance",
                "cache_requalification",
            )
        }
        bridge_path, bridge = hashed_record(
            self.index.get("cache_requalification"), self.path.parent
        )
        build = producer["raw_report"]
        require(
            refs["cache"] == producer["cache_file"]
            and refs["build_report"] == producer["evidence"]
            and refs["producer_workload"]["sha256"] == build["workload_sha256"]
            and refs["consumer_workload"]["sha256"]
            == bridge["consumer"]["final_workload"]["sha256"]
            and refs["cache_requalification"] == evidence(bridge_path),
            "final selection names another cache, build, workload or bridge",
        )
        candidate = read(refs["candidate_acceptance"]["path"])
        require(
            candidate.get("schema")
            == "vakint.strict-complete-small-candidate-acceptance.v1"
            and candidate.get("state")
            == "STRICT_COMPLETE_CANDIDATE_NOT_FINAL_SELECTION"
            and candidate["reports"]["build"]["sha256"]
            == producer["evidence"]["sha256"]
            and candidate["cache_file"]["sha256"] == producer["cache_file"]["sha256"]
            and candidate["workload_sha256"] == build["workload_sha256"]
            and candidate["selected_training_roots"]
            == candidate["completed_strict_roots"]
            == build["preparation"]["total_roots"],
            "final selection names another completed candidate",
        )
        _, workload = hashed_record(selection["consumer_workload"], path.parent)
        _, cohort = hashed_record(
            bridge["consumer"]["compiled_cohort"], bridge_path.parent
        )
        _, proof = hashed_record(bridge["consumer"]["owner_proof"], bridge_path.parent)
        identity = {
            "cache_sha256": refs["cache"]["sha256"],
            "data_manifest_sha256": workload["data_manifest_sha256"],
            "source_census_sha256": workload["source_census_sha256"],
            "workload_sha256": refs["consumer_workload"]["sha256"],
            "build_git_revision": cohort["main_revision"],
            "executable_sha256": proof["measurement_executable_sha256"],
            "owner_proof": {
                **bridge["consumer"]["owner_proof"],
                "compiled_cohort": cohort,
            },
        }
        # Reuse the full provenance/semantic checks, even for an unchanged cohort.
        qualification = self.cache_qualification(identity, producer)
        return {
            "evidence": evidence(path),
            "raw_receipt": selection,
            "cache_requalification": qualification,
        }

    def qualify(self, bundle, binary, case_id, variant, sample, supervision=None):
        launcher, artifact = bundle["launcher"], bundle["artifact"]
        require(bundle["case_id"] == case_id, "evidence belongs to another case")
        require(
            launcher.get("schema")
            not in (
                "vakint.primary-full-functional-run.v1",
                "vakint.q0-full-weighted-launch.v1",
            )
            and not launcher.get(
                "functional_only_not_isolated_three_run_performance", False
            )
            and launcher.get("final_performance_eligible") is True,
            "functional launcher is not final performance evidence",
        )
        require(
            launcher.get("fresh_process") is True
            and launcher.get("binary_sha256") == binary,
            "launcher does not identify this fresh executable",
        )
        require(
            launcher.get("case_id") == case_id
            and launcher.get("variant") == variant
            and launcher.get("sample") == sample,
            "launcher sample identity differs",
        )
        # Support the existing single-primary and multi-driver artifact receipts.
        binary_records = artifact_binaries(artifact)
        if isinstance(artifact.get("binary"), dict):
            binaries = [artifact["binary"]["sha256"]]
            codegen = artifact["actual_codegen"]
            require(
                str(codegen.get("optimization", codegen.get("opt_level"))) == "3"
                and codegen.get("fat_lto") is True
                and codegen.get("test") is False,
                "artifact is not the qualified release executable",
            )
        else:
            candidates = binary_records
            binaries = [
                item["sha256"]
                for item in candidates
                if item.get("normal_cli") is True
                and item["profile"]["opt_level"] == "3"
                and item["profile"]["test"] is False
                and item["actual_rustc"]["lto"] == "fat"
            ]
            binary_records = candidates
        for record in binary_records:
            path = record.get("immutable_path")
            if path and Path(path).exists():
                require(
                    digest(path) == record["sha256"], "available executable changed"
                )
        require(
            binary in binaries
            and artifact.get("source_drift") == []
            and artifact.get("source_unchanged", True) is True,
            "artifact binary/source mismatch",
        )
        stamps = artifact["stamps"]
        required = {
            "VAKINT_CACHE_BUILD_GIT_REVISION": 40,
            "VAKINT_CACHE_BUILD_RUSTRED_REVISION": 40,
            "VAKINT_CACHE_BUILD_SYMBOLICA_REVISION": 40,
            "VAKINT_CACHE_BUILD_MAIN_SOURCE_SNAPSHOT_SHA256": 64,
            "VAKINT_CACHE_BUILD_NATIVE_SOURCE_SNAPSHOT_SHA256": 64,
        }
        for field, length in required.items():
            value = stamps.get(field)
            require(
                isinstance(value, str)
                and len(value) == length
                and all(c in "0123456789abcdefABCDEF" for c in value),
                f"missing actual artifact provenance: {field}",
            )
        require(
            stamps.get("VAKINT_CACHE_BUILD_PROFILE") == "release",
            "artifact profile is not release",
        )
        window_path, window = hashed_record(bundle["window"], self.path.parent)
        require(
            window.get("schema") == "vakint.serialized-performance-window.v1"
            and window.get("final_performance_eligible") is True
            and window.get("accepted") is True
            and window.get("serialized") is True
            and window.get("all_processes_drained") is True
            and window.get("source_unchanged") is True,
            "no completed serialized performance window",
        )
        intervals = sorted(
            window["samples"], key=lambda row: row["started_unix_seconds"]
        )
        require(intervals, "empty performance window")
        require(
            len({row["report_sha256"] for row in intervals}) == len(intervals),
            "window duplicates a report",
        )
        window_start = number(window["started_unix_seconds"])
        window_end = number(window["finished_unix_seconds"])
        last_end = None
        matches = []
        for row in intervals:
            start, end = (
                number(row["started_unix_seconds"]),
                number(row["finished_unix_seconds"]),
            )
            require(
                window_start <= start < end <= window_end
                and (last_end is None or start >= last_end),
                "performance samples overlap or have invalid intervals",
            )
            require(
                row.get("exit_code") == 0 and row.get("fresh_process") is True,
                "window includes an unsuccessful or reused process",
            )
            last_end = end
            if row["report_sha256"] == bundle["references"]["record"]["sha256"]:
                matches.append(row)
        require(len(matches) == 1, "window omits or duplicates this report")
        row = matches[0]
        if supervision is not None:
            require(
                bundle["record"]["supervision"] == supervision,
                "sample supervision differs from its hashed report",
            )
            require(
                supervision.get("interval_scope")
                == "supervisor pre-spawn setup through supervisor-observed worker completion",
                "missing actual supervisor observation scope",
            )
            for field in ("child_pid", "supervisor_pid"):
                require(
                    type(supervision.get(field)) is int and supervision[field] > 0,
                    "missing actual supervisor process identity",
                )
            for field in ("started_unix_seconds", "finished_unix_seconds"):
                number(supervision.get(field))
            for field in (
                "child_pid",
                "supervisor_pid",
                "started_unix_seconds",
                "finished_unix_seconds",
                "interval_scope",
            ):
                require(
                    row.get(field) == launcher.get(field) == supervision[field],
                    "window/launcher differs from actual supervisor observations",
                )
        require(
            row["launcher_sha256"] == bundle["references"]["launcher"]["sha256"]
            and row["artifact_sha256"] == bundle["references"]["artifact"]["sha256"]
            and row["case_id"] == case_id
            and row["variant"] == variant
            and row["sample"] == sample,
            "window evidence binding differs",
        )
        affinity = window.get("cpu_affinity")
        threads = window.get("thread_environment")
        require(
            isinstance(affinity, str)
            and affinity
            and launcher.get("cpu_affinity") == affinity,
            "affinity not fixed by launcher/window",
        )
        require(
            isinstance(threads, dict)
            and all(
                threads.get(k) == "1"
                for k in ("RAYON_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
            )
            and launcher.get("thread_environment") == threads,
            "worker environment not fixed by launcher/window",
        )
        return {
            "window": evidence(window_path),
            "stamps": stamps,
            "binary_sha256": binary,
            "cpu_affinity": affinity,
            "thread_environment": threads,
        }


def number(value):
    require(
        type(value) in (int, float) and math.isfinite(value) and value >= 0,
        "missing, negative or nonfinite timing",
    )
    return value


def summarize(values):
    values = [number(value) for value in values]
    require(len(values) >= 3, "at least three actual samples required")
    return {
        "samples_seconds": values,
        "median_seconds": statistics.median(values),
        "min_seconds": min(values),
        "max_seconds": max(values),
    }


def preserved_native_source(source, base):
    """Read the actual opaque-union proof without promoting it to semantic acceptance."""
    proof_path, proof = hashed_record(source["barrier"], base)
    producer_path, producer = hashed_record(source["producer"], base)
    metadata_path, metadata = hashed_record(source["capture"], base)
    index = hashed_file(source["snapshot"], base)
    require(
        proof.get("schema") == "vakint.postflight-real-union-prefix-proof-qualified.v1"
        and proof.get("data_preservation_state") == "PASS"
        and proof.get("original_runtime_receipts_unchanged") is True,
        "source lacks actual qualified union-prefix proof",
    )
    artifacts = {
        str(hashed_file(ref, proof_path.parent)): ref for ref in proof["artifacts"]
    }
    require(
        str(index) in artifacts and str(metadata_path) in artifacts,
        "prefix snapshot/metadata not bound by union proof",
    )
    terminal = [
        read(Path(name))
        for name in artifacts
        if Path(name).name == "terminal-outcome.json"
    ]
    require(len(terminal) == 1, "union proof lacks terminal receipt")
    terminal = terminal[0]
    verification = terminal["independent_verification"]
    require(
        terminal.get("schema") == "vakint.real-opaque-union-prefix-diagnostic.v1"
        and terminal.get("state") == "PASS"
        and verification.get("schema")
        == "vakint.independent-opaque-union-prefix-preservation.v1"
        and verification.get("state") == "PASS"
        and all(
            verification.get(key) is True
            for key in (
                "every_captured_source_and_output_chunk_sha256_and_length_equal",
                "all_chunk_files_stable_during_reads",
                "native_exact_run_concatenation_verified",
            )
        )
        and verification.get("completed_q0_or_final_cache_claimed") is False
        and verification.get("semantic_qualification_or_cache_speedup_claimed")
        is False,
        "prefix preservation must retain its actual limited scope",
    )
    rows = [
        row
        for row in metadata["immutable_indices"]
        if row["role"] == source["prefix_role"]
    ]
    require(
        len(rows) == 1
        and Path(rows[0]["path"]).resolve() == index
        and rows[0]["bytes"] == index.stat().st_size
        and metadata.get("metadata_only") is True
        and metadata.get("run_concat_verified") is True
        and metadata["generation"] == proof["generation"]
        and metadata["compressed_bytes"] == proof["compressed_bytes"]
        and metadata["preserved_record_references"]
        == proof["stored_record_references"],
        "captured prefix metadata differs from actual union proof",
    )
    require(
        source["scope"] == "committed-prefix",
        "union-prefix proof cannot complete an example",
    )
    require(
        source["prefix_role"] in ("old", "q0"), "union output is not a source prefix"
    )
    if producer.get("schema") == "vakint.five-loop-cache-build.v2":
        require(source["prefix_role"] == "old", "old producer cannot claim q0 prefix")
        proposal_path, proposal = hashed_record(proof["proposal"], proof_path.parent)
        binding_path, binding = hashed_record(proposal["binding"], proposal_path.parent)
        require(
            hashed_file(binding["original_build_receipt"], binding_path.parent)
            == producer_path
            and producer.get("complete") is False
            and producer["preparation"]["status"] != "complete"
            and producer["preparation"]["native_profile"]["strict"] is True
            and producer["preparation"]["native_profile"]["unresolved_policy"]
            == "reject",
            "prefix source must preserve its authenticated failed/incomplete producer",
        )
    else:
        proposal_path, proposal = hashed_record(proof["proposal"], proof_path.parent)
        require(
            source["prefix_role"] == "q0"
            and producer.get("schema")
            == "vakint.inactive-capture-physical-user-unit.v1"
            and any(
                row["role"] == "q0"
                and Path(row["path"]).resolve()
                == Path(producer["running_cache_directory"]).resolve()
                for row in proposal["inputs_newest_first"]
            ),
            "prefix producer is not a captured union input",
        )
    return {
        **source,
        "receipts": {
            "producer": {"evidence": evidence(producer_path), "raw_receipt": producer},
            "capture": {"evidence": evidence(metadata_path), "raw_receipt": metadata},
            "snapshot": {"evidence": evidence(index)},
            "barrier": {"evidence": evidence(proof_path), "raw_receipt": proof},
        },
    }


def strict_preparation_source(source, base):
    """Join an original finite preparation to its later full Native check.

    Preserve both cohorts and the producer's original data identity. A
    selected-key reimport cannot replace the full source-parity receipt.
    Completing this finite preparation never completes a physical example.
    """
    require(
        source.get("kind") == "strict-preparation"
        and source.get("scope") == "completed-preparation",
        "finite preparation must retain its own source kind and scope",
    )
    documents = {}
    for name in (
        "producer",
        "acceptance",
        "preservation",
        "native_validation",
        "native_supervision",
        "native_artifacts",
        "inventory",
    ):
        path, document = hashed_record(source[name], base)
        documents[name] = {"evidence": evidence(path), "raw_receipt": document}
    producer, accepted, preserved, native, supervised, artifact, inventory = (
        documents[name]["raw_receipt"] for name in documents
    )
    preparation = producer["preparation"]
    roots = preparation["total_roots"]
    single_file = preparation["options"]["single_file"]
    count = single_file["records"]
    original_cohort = producer["owner_proof"]["compiled_cohort"]
    require(
        producer.get("schema") == "vakint.five-loop-cache-build.v2"
        and producer.get("complete") is True
        and producer.get("publishable") is True
        and isinstance(producer.get("data_manifest_sha256"), str)
        and re.fullmatch(r"[0-9a-f]{64}", producer["data_manifest_sha256"])
        and producer["data_manifest_sha256"]
        == original_cohort["data_manifest_sha256"]
        == producer["build_identity"]["data_manifest_sha256"]
        and producer["build_identity"]["build"]["git_revision"]
        == original_cohort["main_revision"]
        and preparation["status"] == "complete"
        and type(roots) is int
        and roots > 0
        and same_json_value(preparation["completed_roots"], roots)
        and all(
            preparation.get(key) is None
            for key in ("failure", "resource_failure", "checkpoint_failure")
        )
        and preparation["native_profile"]["strict"] is True
        and preparation["native_profile"]["unresolved_policy"] == "reject"
        and producer["supervision"]["success"] is True
        and producer["supervision"]["hard_deadline_exhausted"] is False
        and type(count) is int
        and count > 0,
        "finite Native source lacks an actually completed strict preparation",
    )
    require(
        accepted.get("schema") == "vakint.strict-complete-small-candidate-acceptance.v1"
        and accepted.get("state") == "STRICT_COMPLETE_CANDIDATE_NOT_FINAL_SELECTION"
        and accepted["reports"]["build"]["sha256"] == source["producer"]["sha256"]
        and same_json_value(
            [
                accepted["selected_training_roots"],
                accepted["completed_strict_roots"],
                accepted["record_count"],
            ],
            [roots, roots, count],
        )
        and accepted["cache_file"]["sha256"] == producer["cache_sha256"]
        and accepted["cache_file"]["decoded_sha256"] == single_file["decoded_sha256"]
        and same_json_value(
            accepted["cache_file"]["native_identity"], single_file["identity"]
        ),
        "finite source acceptance differs from its actual producer",
    )
    require(
        preserved.get("schema") == "vakint.accepted-small-native-cache-preservation.v1"
        and preserved.get("state") == "PASS_STORAGE_PRESERVATION"
        and all(
            preserved.get(key) is True
            for key in (
                "all_source_and_destination_compressed_bytes_equal",
                "all_gzip_members_crc_lengths_hashes_and_header_bindings_valid",
                "all_original_sources_stable",
                "all_output_inodes_independent",
                "read_only_preserved_bundle",
            )
        )
        and single_file["checkpoint_path"] == preserved["original_bundle"]
        and preserved["accepted_gzip"]["sha256"] == producer["cache_sha256"]
        and same_json_value(preserved["stored_record_references"], count),
        "finite source lost its original checkpoint/storage lineage",
    )
    identity = inventory["cache_native_index_identity"]
    recorded_indexes = {
        row["relative_path"]: {
            "sha256": row["sha256"],
            "bytes": row["destination"]["bytes"],
        }
        for row in preserved["indices"]
    }
    require(
        same_json_value(recorded_indexes, identity["indexes"])
        and native.get("schema")
        == "vakint.strict-cache-source-parity-native-validation.v1"
        and native.get("success") is True
        and same_json_value(native["source_cache_native_index_identity"], identity)
        and native["source_inventory_sha256"] == source["inventory"]["sha256"]
        and same_json_value(native["semantic_binding"], single_file["binding"])
        and native["application_sha256"]
        == producer["owner_proof"]["compiled_cohort"]["application_sha256"]
        and native["legacy_cache_sha256"] == producer["cache_sha256"]
        and same_json_value(
            [
                native[key]
                for key in (
                    "native_unique_records",
                    "legacy_validated_records",
                    "source_native_stored_record_references",
                )
            ],
            [count, count, count],
        )
        and same_json_value(native["source_native_generation"], preserved["generation"])
        and same_json_value([native["missing_keys"], native["partial_records"]], [0, 0])
        and all(
            native.get(key) is True
            for key in (
                "exact_record_and_ordered_guard_equality",
                "declared_terminal_validation",
                "source_unchanged",
            )
        )
        and native.get("functional_IBP_comparison") is False
        and native.get("new_reduction_or_capture") is False
        and same_json_value(native["codec_state_format"], 6)
        and native["collector_artifacts"]["sha256"]
        == source["native_artifacts"]["sha256"],
        "finite source lacks full bound current-owner Native/ordered-guard parity",
    )
    validate_compiled_artifact(native, artifact)
    admission = native["resource_guard_admission"]
    require(
        supervised.get("schema") == "vakint.guarded-small-native-gzip-source-parity.v1"
        and supervised.get("accepted") is True
        and same_json_value(supervised["worker_exit_code"], 0)
        and supervised["errors"] == []
        and supervised["raw_validation"]["sha256"]
        == source["native_validation"]["sha256"]
        and admission.get("admitted") is True
        and admission.get("schema") == "vakint.resource-guard-admission.v1"
        and same_json_value(admission, supervised["admission_receipt"])
        and same_json_value(
            admission["worker"],
            {
                **{
                    key: supervised["worker"][key]
                    for key in ("pid", "start_ticks", "uid")
                },
                "executable_sha256": native["executable_sha256"],
            },
        )
        and supervised["registration"]["signals_sent"] == []
        and supervised["guard_after_worker_exit"]["ActiveState"] == "inactive",
        "finite Native source lacks its actual successful guarded validation run",
    )
    return {
        **source,
        "receipts": documents,
        "original_producer_cohort": original_cohort,
        "current_validation_cohort": native["compiled_cohort"],
        "inventory_counts": {
            "compressed_bytes": preserved["compressed_chunk_bytes"],
            "stored_record_references": count,
            "physical_chunks": preserved["chunks"],
        },
        "scope_note": "Completed finite strict preparation; no completed physical example or final selection.",
    }


def cumulative_cache_record(path, cache, physical_scope=None, *, source_identity=None):
    """Consume accepted snapshot evidence; a preserved prefix is not a finished run."""
    path = Path(path).resolve()
    report = read(path)
    require(
        report.get("schema") == CUMULATIVE_SCHEMA, "unknown cumulative source schema"
    )
    if source_identity is None:
        identity = native_cache_identity(cache)
    else:
        require(
            cache is None and isinstance(source_identity, dict) and source_identity,
            "archived source identity must remain distinct from a runtime cache directory",
        )
        identity = validated_native_cache_identity(source_identity)
    require(
        report["cache_sha256"] is None
        and same_json_value(report["cache_native_index_identity"], identity),
        "cumulative native indexes changed",
    )
    require(
        isinstance(report["semantic_binding"], dict) and report["semantic_binding"],
        "missing cumulative native binding",
    )
    validation_path, validation = hashed_record(report["validation"], path.parent)
    preservation_only = (
        validation.get("schema")
        == "vakint.postflight-real-union-prefix-proof-qualified.v1"
    )
    if preservation_only:
        require(
            cache is not None,
            "archived LIGHT source needs native validation, not prefix preservation alone",
        )
        require(
            validation.get("data_preservation_state") == "PASS"
            and validation.get("original_runtime_receipts_unchanged") is True,
            "invalid prefix preservation proof",
        )
        require(
            Path(validation["output"]).resolve() == Path(cache).resolve(),
            "prefix proof names another cache",
        )
        require(
            report["inventory"]
            == {
                key: validation[key]
                for key in (
                    "compressed_bytes",
                    "stored_record_references",
                    "physical_chunks",
                )
            },
            "prefix inventory differs from actual preservation proof",
        )
        artifacts = {
            str(hashed_file(ref, validation_path.parent)): ref
            for ref in validation["artifacts"]
        }
        for name, expected in identity["indexes"].items():
            require(
                artifacts.get(str((Path(cache) / name).resolve()), {}).get("sha256")
                == expected["sha256"],
                "prefix proof names another native index",
            )
    else:
        require(
            validation.get("schema") == "vakint.cumulative-native-cache-validation.v1"
            and validation.get("accepted") is True
            and all(
                validation.get("checks", {}).get(key) is True
                for key in (
                    "immutable_index_snapshot",
                    "authenticated_chunks",
                    "strict_completed_records_only",
                    "authenticated_bindings",
                    "committed_source_snapshots",
                )
            ),
            "cumulative snapshot lacks accepted native validation",
        )
        for key in (
            "cache_native_index_identity",
            "semantic_binding",
            "data_manifest_sha256",
            "sources",
            "inventory",
        ):
            require(
                validation[key] == report[key],
                "cumulative validation snapshot/binding/source inventory differs",
            )
    inventory = report["inventory"]
    for key in ("compressed_bytes", "stored_record_references", "physical_chunks"):
        require(
            type(inventory.get(key)) is int and inventory[key] > 0,
            "invalid cumulative stored inventory",
        )
    sources, ids, completed_physical = [], set(), set()
    require(report["sources"], "empty cumulative source inventory")
    for source in report["sources"]:
        require(
            source["id"] not in ids
            and (
                source.get("kind") == "strict-preparation"
                and source["scope"] == "completed-preparation"
                or source.get("kind") != "strict-preparation"
                and source["scope"] in ("committed-prefix", "completed-example")
            ),
            "duplicate/invalid cumulative source",
        )
        ids.add(source["id"])
        if source.get("kind") == "strict-preparation":
            prepared = strict_preparation_source(source, path.parent)
            require(
                same_json_value(
                    prepared["receipts"]["native_validation"]["raw_receipt"][
                        "semantic_binding"
                    ],
                    report["semantic_binding"],
                ),
                "finite preparation differs from the cumulative binding",
            )
            if len(report["sources"]) == 1:
                require(
                    same_json_value(
                        prepared["receipts"]["inventory"]["raw_receipt"][
                            "cache_native_index_identity"
                        ],
                        report["cache_native_index_identity"],
                    )
                    and same_json_value(prepared["inventory_counts"], inventory)
                    and prepared["receipts"]["producer"]["raw_receipt"][
                        "data_manifest_sha256"
                    ]
                    == report["data_manifest_sha256"],
                    "one-source finite report substituted its snapshot, inventory or original data",
                )
            sources.append(prepared)
            continue
        if source.get("kind") == "native-index-prefix":
            require(
                source["barrier"]["sha256"] == report["validation"]["sha256"]
                if preservation_only
                else True,
                "prefix source belongs to another preservation proof",
            )
            sources.append(preserved_native_source(source, path.parent))
            continue
        require(
            source.get("kind") == "running-capture", "unknown cumulative source kind"
        )
        docs = {}
        for name in ("producer", "capture", "snapshot", "barrier"):
            refpath, document = hashed_record(source[name], path.parent)
            docs[name] = {"evidence": evidence(refpath), "raw_receipt": document}
        barrier = docs["barrier"]["raw_receipt"]
        require(
            barrier.get("schema") == "vakint.native-cache-snapshot-barrier.v1"
            and barrier.get("accepted") is True
            and barrier.get("committed_snapshot") is True
            and barrier.get("writer_fault") is None
            and type(barrier.get("generation")) is int
            and barrier["generation"] > 0
            and isinstance(barrier.get("manifest_digest"), str)
            and barrier["manifest_digest"],
            "source lacks accepted committed snapshot barrier",
        )
        for name in ("producer", "capture", "snapshot"):
            require(
                hashed_file(
                    barrier[name], Path(docs["barrier"]["evidence"]["path"]).parent
                )
                == Path(docs[name]["evidence"]["path"]),
                "snapshot barrier names another source",
            )
        snapshot = docs["snapshot"]["raw_receipt"]
        index_path = hashed_file(source["index"], path.parent)
        require(
            Path(snapshot["delta"]["manifest_path"]).resolve() == index_path
            and snapshot["delta"]["manifest_bytes"] == index_path.stat().st_size,
            "source snapshot index differs",
        )
        docs["index"] = {"evidence": evidence(index_path)}
        capture = docs["capture"]["raw_receipt"]
        native = capture.get("report", capture.get("native", capture))
        require(
            snapshot.get("schema") == "rustred.running-application-cache-snapshot.v1"
            and snapshot["delta"]["generation"] == barrier["generation"]
            and snapshot["delta"]["manifest_digest"] == barrier["manifest_digest"]
            and native.get("schema") == "rustred.running-application-cache.v1"
            and native.get("writer_error") is None
            and native.get("last_error") is None
            and native.get("generation") == barrier["generation"]
            and all(
                type(native.get(key)) is int and native[key] >= 0
                for key in (
                    "committed_records",
                    "accepted_records",
                    "pending_records",
                    "pending_bytes",
                )
            )
            and native["committed_records"] <= native["accepted_records"]
            and all(
                native.get(key) == barrier.get(key)
                for key in (
                    "accepted_records",
                    "committed_records",
                    "pending_records",
                    "pending_bytes",
                )
            ),
            "snapshot barrier differs from actual native capture/snapshot",
        )
        require(
            barrier["producer_completed"] is (source["scope"] == "completed-example"),
            "partial source cannot become completed example",
        )
        if source["scope"] == "completed-example":
            require(
                docs["producer"]["raw_receipt"].get("schema")
                != "vakint.five-loop-cache-build.v2",
                "failed finite workload cannot become completed physical example",
            )
            # Reuse the complete physical validator, including its actual input,
            # ELF/source, output, full guards, close and supervisor bindings.
            from merge_physical_capture_results import completed_run

            config_path = hashed_file(source["config"], path.parent)
            config = read(config_path)
            jobs = {"zero-momentum-metric": "q0", "quadratic-transverse": "quadratic"}
            projection = config["prepared_input"]["manifest_projection"]
            require(projection in jobs, "unknown completed physical projection")
            completed = completed_run(config_path, expected_job_id=jobs[projection])
            require(
                completed["completion"] == docs["producer"]["raw_receipt"]
                and completed["cache_report"] == capture
                and completed["cache_report"]["snapshot"] == snapshot,
                "completed source differs from accepted physical run",
            )
            docs["config"] = {
                "evidence": evidence(config_path),
                "raw_receipt": completed["config"],
            }
            completed_physical.update(
                case
                for case, projection in PRIMARY.items()
                if completed["report"]["projection"] == projection
            )
        sources.append({**source, "receipts": docs})
    unions = {}
    for source in sources:
        if source.get("kind") == "native-index-prefix":
            unions.setdefault(source["barrier"]["sha256"], []).append(source)
    for group in unions.values():
        expected_roles = {
            row["role"]
            for row in group[0]["receipts"]["capture"]["raw_receipt"][
                "immutable_indices"
            ]
        } - {"output"}
        roles = [row["prefix_role"] for row in group]
        require(
            len(roles) == len(set(roles)) and set(roles) == expected_roles,
            "union must retain every captured source prefix exactly once",
        )
    return {
        "evidence": evidence(path),
        **(
            {"cache_directory": str(Path(cache).resolve())}
            if cache is not None
            else {
                "source_index_evidence_scope": "Source receipts and the strict export's authenticated HEAVY index identity; the HEAVY payload was not re-read."
            }
        ),
        "raw_report": report,
        "validation": {
            "evidence": evidence(validation_path),
            "raw_receipt": validation,
        },
        "sources": sources,
        "status": "preserved-committed-prefix"
        if preservation_only
        else "validated-committed-snapshot",
        "all_source_runs_completed": all(
            row["scope"] in ("completed-example", "completed-preparation")
            for row in report["sources"]
        ),
        **physical_completion(completed_physical, physical_scope),
        "performance_eligible": False,
        "scope": "Stored references may overlap; snapshot preservation does not establish producer completion, unique keys, or performance.",
    }


def build_record(path, cache, provenance=None, physical_scope=None):
    build = read(path)
    if build.get("schema") == LIGHT_SCHEMA:
        require(
            provenance is not None,
            "LIGHT producer requires its source-bound consumer evidence",
        )
        return provenance.light_cache_record(path, cache, physical_scope)
    if build.get("schema") == "vakint.cumulative-native-cache-report.v1":
        return cumulative_cache_record(path, cache, physical_scope)
    require(
        build.get("schema") == "vakint.five-loop-cache-build.v2", "unknown build schema"
    )
    prep = build["preparation"]
    require(
        build["phase"] in ("final", "pilot")
        and build["complete"] is True
        and build["publishable"] is True
        and prep["status"] == "complete",
        "only completed cache preparation may enter a measured report",
    )
    require(
        build["supervision"]["success"] is True
        and build["supervision"]["hard_deadline_exhausted"] is False,
        "cache preparation lacks successful supervision",
    )
    require(
        prep["total_roots"] > 0 and prep["completed_roots"] == prep["total_roots"],
        "selected training roots are incomplete",
    )
    require(
        all(
            prep.get(field) is None
            for field in ("failure", "resource_failure", "checkpoint_failure")
        ),
        "cache preparation carries a failure",
    )
    require(
        prep["native_profile"]["strict"] is True
        and prep["native_profile"]["unresolved_policy"] == "reject",
        "cache records must be strict complete reductions",
    )
    require(
        build["cache_sha256"] == digest(cache)
        and build["compressed_bytes"] == Path(cache).stat().st_size,
        "cache file differs from completed build",
    )
    # Historical size-target fields remain provenance, not an aggregate cap.
    producer = {
        "evidence": evidence(path),
        "cache_file": evidence(cache),
        "raw_report": build,
    }
    if build["phase"] == "pilot":
        require(provenance is not None, "pilot cache requires final-selection evidence")
        producer["final_selection"] = provenance.final_selection(producer)
    return producer


def suite_record(path, kind, producer, provenance):
    build = producer["raw_report"]
    path = Path(path).resolve()
    suite = read(path)
    require(
        suite.get("schema") == "vakint.five-loop-cache-suite.v1", "unknown suite schema"
    )
    identity = suite["identity"]
    require(
        identity["suite"] == kind
        and suite["complete"] is True
        and suite["provisional"] is False
        and suite["native_keys_authenticated"] is True
        and suite["execution_provenance_complete"] is True,
        "suite is not qualified",
    )
    require(
        suite[
            "final_held_out_comparison" if kind == "held-out" else "selection_allowed"
        ]
        is True,
        "driver did not qualify this suite",
    )
    count = identity["repetitions"]
    require(type(count) is int and count >= 3, "three repetitions required")
    requalification = provenance.cache_qualification(identity, producer)
    require(
        identity["build_profile"] == "release"
        and identity["rayon_num_threads"] == "1"
        and identity["cpu_affinity"]
        and identity["load_prepared_caches"] is False,
        "final measurement profile/worker/cache policy differs",
    )
    raw = suite["measurements"]
    ids = [case["id"] for case in suite["cases"]]
    require(
        len(ids) == len(set(ids)) and len(raw) == len(ids) * count * 2,
        "missing or duplicate case/variant samples",
    )
    require(
        kind != "held-out" or len(ids) >= 12, "at least twelve held-out cases required"
    )
    rows, seen, root_sets, contexts = [], set(), set(), set()
    for index, case in enumerate(suite["cases"]):
        case_id = case["id"]
        require(
            case["successful"] is True
            and case["exact_output_and_guard_parity"] is True
            and case["native_keys_authenticated"] is True,
            "case failed exact validation",
        )
        samples = [row for row in raw if row["case_id"] == case_id]
        require(len(samples) == 2 * count, "case sample count differs")
        require(
            len({row["input_sha256"] for row in samples}) == 1
            and len({row["value_sha256"] for row in samples}) == 1
            and len({row["guards_sha256"] for row in samples}) == 1,
            "case input/output/guards differ across samples or modes",
        )
        protected = samples[0]["frozen_protected_roots"]
        require(
            protected and all(len(key) == 15 for key in protected), "empty owner set"
        )
        root_sets.add(tuple(sorted(tuple(key) for key in protected)))
        require(
            isinstance(case.get("catalogue_closure"), dict),
            "final collector must report catalogue closure separately from Retain success",
        )
        exposure = None
        if build.get("schema") in (CUMULATIVE_SCHEMA, LIGHT_SCHEMA):
            exposure = requalification["raw_receipt"]["case_exposure"].get(case_id)
            require(
                exposure in ("held-out", "tuning", "seen-starting-input")
                and (kind != "held-out" or exposure == "held-out"),
                "seen starting input cannot become held-out evidence",
            )
        output = {
            "id": case_id,
            "status": "measured",
            "catalogue_closure": case["catalogue_closure"],
            "requested_root_hits": case["requested_root_hits"],
            "intermediate_hits": case["intermediate_hits"],
            "backing_reuse": case["backing_reuse"],
            "variants": {},
        }
        if exposure is not None:
            output["cache_exposure"] = exposure
        for sample in samples:
            key = (case_id, sample["variant"], sample["sample"])
            require(
                key not in seen
                and sample["variant"] in ("uncached", "cached")
                and sample["sample"] in range(count),
                "duplicate or invalid sample",
            )
            seen.add(key)
            require(
                sample["identity"] == identity
                and sample["success"] is True
                and sample["worker_success"] is True
                and sample["fresh_process"] is True
                and sample["native_keys_authenticated"] is True
                and sample["supervision"]["success"] is True
                and sample["supervision"]["exit_code"] == 0
                and sample["supervision"]["deadline_exhausted"] is False,
                "unsuccessful or mismatched sample",
            )
            require(
                sample["owner_proof"] == identity["owner_proof"], "owner proof mismatch"
            )
            require(
                sample["frozen_protected_roots"] == protected
                and sample["actual_lowered_roots"] == sample["frozen_lowered_roots"],
                "actual root set differs",
            )
            before = sample["persistent_cache_attachments_before"]
            during = sample["persistent_cache_attachments_during_evaluation"]
            require(
                before == {str(loop): False for loop in range(1, 6)},
                "baseline cache present",
            )
            require(
                during
                == {
                    str(loop): sample["variant"] == "cached" and loop == 5
                    for loop in range(1, 6)
                },
                "wrong persistent backing",
            )
            require(
                sample["cache_statistics"]["binding"] == cache_binding(producer),
                "cache semantic binding mismatch",
            )
            if sample["variant"] == "cached":
                light_attachment(sample["cache_statistics"], producer)
            directory = (
                path.parent
                / f"case-{index:03}-{sample['sample']:02}-{sample['variant']}"
            )
            require(
                read(directory / "measurement.json") == sample, "saved sample differs"
            )
            bundle = provenance.bind(directory / "measurement.json")
            if build.get("schema") in (CUMULATIVE_SCHEMA, LIGHT_SCHEMA):
                expected_cache = (
                    cache_identity(build) if sample["variant"] == "cached" else None
                )
                require(
                    cache_identity(bundle["launcher"]) == expected_cache,
                    "suite launcher native cache identity differs",
                )
            if requalification is not None:
                require(
                    bundle["references"]["artifact"]["sha256"]
                    == requalification["raw_receipt"]["consumer"][
                        "collector_artifacts"
                    ]["sha256"],
                    "sample artifact differs from accepted consumer",
                )
            qualified = provenance.qualify(
                bundle,
                identity["executable_sha256"],
                case_id,
                sample["variant"],
                sample["sample"],
                supervision=sample["supervision"],
            )
            require(
                qualified["stamps"]["VAKINT_CACHE_BUILD_GIT_REVISION"]
                == identity["build_git_revision"]
                and qualified["cpu_affinity"] == identity["cpu_affinity"],
                "suite/artifact/launcher context differs",
            )
            cohort = identity["owner_proof"]["compiled_cohort"]
            for field, stamp in COHORT_STAMPS:
                require(
                    cohort[field] == qualified["stamps"][stamp],
                    "suite compiled cohort differs",
                )
            contexts.add(json.dumps(qualified, sort_keys=True))
            output.setdefault("performance_evidence", []).append(qualified)
            for filename, field in (
                ("value.symbolica", "value_sha256"),
                ("guards.json", "guards_sha256"),
            ):
                require(
                    digest(directory / filename) == sample[field],
                    "result/guard sidecar changed",
                )
            closure = sample["catalogue_closure"]
            unknown = closure["unknown_masters"]
            require(
                unknown == sorted(set(unknown))
                and closure["unknown_count"] == len(unknown)
                and closure["closed"] is (not unknown)
                and read(directory / "unknown-masters.json") == unknown
                and digest(directory / "unknown-masters.json")
                == closure["unknown_sha256"],
                "unknown-master classification or sidecar changed",
            )
        unknown_union = sorted(
            {
                master
                for sample in samples
                for master in sample["catalogue_closure"]["unknown_masters"]
            }
        )
        require(
            output["catalogue_closure"]["unknown_masters"] == unknown_union
            and output["catalogue_closure"]["unknown_count"] == len(unknown_union)
            and output["catalogue_closure"]["closed"] is (not unknown_union),
            "case catalogue closure differs from samples",
        )
        for variant in ("uncached", "cached"):
            selected = sorted(
                (row for row in samples if row["variant"] == variant),
                key=lambda row: row["sample"],
            )
            metrics = {
                field: summarize(row[field] for row in selected) for field in METRICS
            }
            stages = {
                field: summarize(
                    sum(number(v) for v in row["stages"][field]) for row in selected
                )
                if all(field in row["stages"] for row in selected)
                else None
                for field in STAGES
            }
            expected = metrics["ibp_excluding_owner_seconds"]
            reported = case[variant]
            require(
                all(
                    expected[key] == reported[key]
                    for key in ("median_seconds", "min_seconds", "max_seconds")
                ),
                "driver median/range mismatch",
            )
            metrics.update(
                graph_preparation_seconds=None,
                pure_ibp_seconds=None,
                terminal_substitution_seconds=None,
            )
            output["variants"][variant] = {"metrics": metrics, "stages": stages}
        cached_samples = [row for row in samples if row["variant"] == "cached"]
        for field in ("requested_root_hits", "intermediate_hits"):
            require(
                type(case[field]) is int
                and case[field] >= 0
                and all(
                    type(row["cache_statistics"]["persistent"][field]) is int
                    and row["cache_statistics"]["persistent"][field] >= 0
                    for row in cached_samples
                )
                and sum(
                    row["cache_statistics"]["persistent"][field]
                    for row in cached_samples
                )
                == case[field],
                "driver cache-hit aggregate differs",
            )
        # Preserve the collector's counters and label above as reported evidence.
        # In Native419, requested_root_hits counts Route records only: a requested
        # owner's Apply hit is classified as intermediate. Neither a revision
        # string nor that legacy label proves proper-descendant-only reuse.
        output["cache_hit_scope"] = {
            "some_exact_requested_root_reuse_proven": case["requested_root_hits"] > 0,
            "other_hits": (
                "unclassified" if case["intermediate_hits"] > 0 else "none-observed"
            ),
        }
        u = output["variants"]["uncached"]["metrics"]["ibp_excluding_owner_seconds"][
            "median_seconds"
        ]
        c = output["variants"]["cached"]["metrics"]["ibp_excluding_owner_seconds"][
            "median_seconds"
        ]
        output["uncached_over_cached_ratio"] = u / c if c > 0 else None
        require(
            output["uncached_over_cached_ratio"]
            == case.get("uncached_over_cached_ratio"),
            "driver ratio mismatch",
        )
        rows.append(output)
    require(
        len(contexts) == 1,
        "suite samples require one cohort and adopted performance window",
    )
    require(
        kind != "held-out" or len(root_sets) >= 12,
        "aliased owner sets cannot count as twelve",
    )
    return {
        "evidence": evidence(path),
        "identity": identity,
        "rows": rows,
        "raw_suite": suite,
        "cache_requalification": requalification,
        "metric": "saved application minus owner preparation",
        "not_separately_timed": [
            "pure IBP",
            "terminal substitution",
            "graph preparation",
        ],
    }


def primary_input(report, directory):
    """Bind persisted bytes/provenance; scientific validation belongs to the native loader."""
    mode = report.get("input_mode", "fresh")
    require(
        mode in ("fresh", "prepared")
        and report["fresh_native_preparation"] is (mode == "fresh"),
        "primary input/preparation mode differs",
    )
    persisted = report.get("prepared_input")
    if mode == "prepared":
        require(isinstance(persisted, dict), "prepared input requires its manifest")
        require(
            report.get("preparation") is None
            and report.get("preparation_seconds") is None,
            "historical preparation cannot become current-run work",
        )
        number(report.get("input_loading_seconds"))
    else:
        number(report["preparation_seconds"])
        require(
            report.get("input_loading_seconds") is None,
            "fresh preparation cannot report prepared-input loading",
        )
    input_path = Path(directory) / "vakint-input.symbolica"
    if persisted is not None:
        path, manifest = hashed_record(
            {
                "path": persisted["manifest_path"],
                "sha256": persisted["manifest_sha256"],
            },
            Path(directory),
        )
        require(
            path.name == "prepared-input.json"
            and manifest == persisted["manifest"]
            and manifest["schema"] == "vakint.qcd-ladder.prepared-input.v1"
            and manifest["definition"] == "massless-qcd-four-rung-native-uv-v1"
            and manifest["projection"] == report["projection"],
            "prepared manifest identity differs",
        )
        assembly = manifest["scientific_assembly_sha256"]
        require(
            isinstance(assembly, str)
            and len(assembly) == 64
            and all(c in "0123456789abcdefABCDEF" for c in assembly),
            "missing scientific assembly identity",
        )
        for name, entry in (
            ("vakint-input.symbolica", manifest["input"]),
            ("projector-conditions.symbolica", manifest["nonzero_conditions_file"]),
        ):
            payload = path.parent / name
            require(
                payload.stat().st_size == entry["bytes"]
                and digest(payload) == entry["sha256"],
                "prepared input/guard payload changed",
            )
        require(
            (path.parent / "projector-conditions.symbolica").read_text()
            == "\n".join(manifest["nonzero_conditions"])
            and report["prepared_nonzero_conditions"] == manifest["nonzero_conditions"],
            "prepared projector guards differ",
        )
        producer, consumer = (
            manifest["producer"]["fixture_sources_sha256"],
            report["consumer_fixture_sources_sha256"],
        )
        require(
            report["fixture_sources_sha256"] == producer
            and set(manifest["scientific_sources_sha256"])
            == {
                "physical.rs",
                "projector.rs",
                "vacuum.rs",
                "ladder.dot",
                "qcd.json",
                "phi3.json",
            }
            and all(
                producer.get(name) == consumer.get(name) == value
                for name, value in manifest["scientific_sources_sha256"].items()
            )
            and report["contributions"] == manifest["contributions"]
            and report["scalar_terms"] == manifest["scalar_terms"],
            "prepared producer/consumer definition differs",
        )
        number(manifest["preparation_seconds"])
        for value in manifest["timings"].values():
            number(value)
        if mode == "fresh":
            require(
                report["preparation_seconds"] == manifest["preparation_seconds"]
                and report["preparation"] == manifest["timings"],
                "fresh preparation receipt differs",
            )
        input_path = path.parent / "vakint-input.symbolica"
        require(
            report["input_sha256"] == manifest["input"]["sha256"],
            "prepared scalar identity differs",
        )
    require(digest(input_path) == report["input_sha256"], "primary input changed")
    return {"mode": mode, "prepared": persisted if mode == "prepared" else None}


def accepted_generator(reference, projection, base):
    """Read existing scientific acceptance; this does not issue a new acceptance."""
    path, accepted = hashed_record(reference, Path(base))
    expected = {
        "quadratic-transverse": (
            "vakint.current-quadratic-preparation.scientific-parity.v1",
            "PASS_ALL21_FRESH_CURRENT_DEPENDENCY_RAW_INPUT_PARITY",
            21,
        ),
        "zero-momentum-metric": (
            "vakint.current-q0-preparation.scientific-parity.v1",
            "PASS_Q0_FRESH_CURRENT_DEPENDENCY_RAW_INPUT_PARITY",
            1,
        ),
    }[projection]
    require(
        (accepted.get("schema"), accepted.get("status"), accepted.get("contributions"))
        == expected
        and accepted.get("collector_exit_code") == 0
        and accepted.get("source_drift") == []
        and accepted["input"]["exact_canonical_bytes_equal_accepted_input"] is True
        and accepted["projector_guards"]["bytes_equal_accepted"] is True,
        "primary generator acceptance is incomplete",
    )
    verified = evidence(path)
    if reference.get("provenance_addendum") is not None:
        extra_path, extra = hashed_record(reference["provenance_addendum"], Path(base))
        require(
            extra["schema"] == "vakint.current-q0-preparation.provenance-addendum.v1"
            and extra["status"]
            == "CORRECTS_FEYNKIT_OWNER_REVISION_ONLY_RAW_PARITY_UNCHANGED"
            and extra["original_scientific_parity"]["sha256"] == verified["sha256"]
            and extra["actual_compiled_workspace_revision"] == accepted["source_commit"]
            and extra["immutable_elf"]["sha256"] == accepted["immutable_elf_sha256"]
            and extra["input_sha256"] == accepted["input"]["sha256"]
            and extra["scientific_definition_and_raw_input_changes"] is False,
            "generator provenance addendum changed the scientific identity",
        )
        verified["provenance_addendum"] = evidence(extra_path)
    return verified, accepted


def primary_generator(report, launcher, base):
    """Bind the accepted physical input, not merely mutually equal samples."""
    reference, accepted = accepted_generator(
        launcher.get("generator_acceptance"), report["projection"], base
    )
    require(
        accepted["projection"] == report["projection"]
        and accepted["input"]["sha256"] == report["input_sha256"]
        and accepted["input"]["bytes"] == report["input_bytes"]
        and accepted["contributions"] == report["contributions"]
        and accepted["projector_guards"]["values"]
        == report["prepared_nonzero_conditions"]
        and accepted["projector_guards"]["sha256"]
        == hashlib.sha256(
            "\n".join(report["prepared_nonzero_conditions"]).encode()
        ).hexdigest(),
        "primary differs from the accepted complete physical fixture",
    )
    consumer = report.get(
        "consumer_fixture_sources_sha256", report["fixture_sources_sha256"]
    )
    require(
        all(
            consumer.get(name) == accepted["compiled_fixture_source_hashes"].get(name)
            for name in (
                "physical.rs",
                "projector.rs",
                "vacuum.rs",
                "ladder.dot",
                "qcd.json",
                "phi3.json",
            )
        ),
        "primary scientific sources differ from accepted generator",
    )
    if report.get("input_mode", "fresh") == "prepared":
        manifest = report["prepared_input"]["manifest"]
        require(
            manifest["producer"]["fixture_sources_sha256"]
            == accepted["compiled_fixture_source_hashes"]
            and manifest["producer"]["source_commit"] == accepted["source_commit"]
            and manifest["producer"]["executable_sha256"]
            == accepted["immutable_elf_sha256"],
            "prepared producer provenance differs from accepted generator",
        )
    return reference


def primary_record(bundle):
    """Retain original functional provenance; verify full outputs when present."""
    report = bundle["record"]
    case_id = bundle["case_id"]
    require(case_id in PRIMARY, "unknown primary case")
    schema = report.get("schema")
    require(
        schema
        in (
            "vakint.qcd-ladder.benchmark.v1",
            "vakint.q0-full-weighted-native-diagnostic.v2",
        ),
        "unknown primary report schema",
    )
    result = {
        "evidence": bundle["references"],
        "raw_report": report,
        "launcher": bundle["launcher"],
        "artifact": bundle["artifact"],
        "performance_eligible": False,
        "catalogue_closure": None,
        "status": "functional-incomplete",
    }
    if schema == "vakint.q0-full-weighted-native-diagnostic.v2":
        require(
            case_id == "zero-momentum-diagnostic" and bundle["window"] is None,
            "exploratory native q0 cannot become final primary performance",
        )
        artifact = bundle["artifact"]
        cohort = report["compiled_cohort"]
        require(
            all(
                cohort.get(key) == artifact.get(key) and cohort.get(key) is not None
                for key in (
                    "native",
                    "symbolica_revision",
                    "driver_source_sha256",
                    "input_sha256",
                    "data",
                    "native_profile",
                )
            )
            and report["input_sha256"] == artifact["input_sha256"],
            "exploratory q0 cohort differs",
        )
        bindings = bundle["launcher"]["bindings"]
        require(
            any(
                row["sha256"] == bundle["references"]["artifact"]["sha256"]
                for row in bindings
            )
            and any(row["sha256"] == artifact["binary_sha256"] for row in bindings),
            "exploratory launcher does not bind this artifact and executable",
        )
        result["status"] = "exploratory-functional"
        return result
    require(report["projection"] == PRIMARY[case_id], "primary projection differs")
    require(
        report["executable_sha256"]
        in {row["sha256"] for row in artifact_binaries(bundle["artifact"])}
        and report["executable_sha256"] == bundle["launcher"]["binary_sha256"],
        "primary report/launcher/artifact executable mismatch",
    )
    if report.get("retained_result_available") is not True:
        require(
            bundle["window"] is None,
            "incomplete primary cannot enter a performance window",
        )
        return result
    directory = Path(bundle["references"]["record"]["path"]).parent
    input_identity = primary_input(report, directory)
    result["input_mode"] = input_identity["mode"]
    result["prepared_input"] = input_identity["prepared"]
    require(
        digest(directory / "result.symbolica") == report["result_sha256"],
        "primary exact result changed",
    )
    closure = report["catalogue_closure"]
    unknown_path = directory / closure["unknown_file"]
    unknown = read(unknown_path)
    require(
        isinstance(unknown, list)
        and all(isinstance(item, str) for item in unknown)
        and unknown == sorted(set(unknown))
        and digest(unknown_path) == closure["unknown_sha256"]
        and closure["unknown_count"] == len(unknown)
        and closure["master_count"] >= len(unknown)
        and closure["closed"] is (not unknown),
        "invalid primary unknown inventory",
    )
    guards = report["guards"]
    guard_path = directory / guards["file"]
    require(
        guards["includes_preparation_conditions"] is True
        and digest(guard_path) == guards["sha256"]
        and len(guard_path.read_text().splitlines()) == guards["count"],
        "missing or changed complete primary guards",
    )
    result["catalogue_closure"] = {**closure, "unknown_masters": unknown}
    result["status"] = "functional-closed" if not unknown else "retained-incomplete"
    require(
        report["reduction_completed"] is (not unknown),
        "primary completion and declared catalogue closure disagree",
    )
    return result


def primary_records(provenance, physical_scope=None):
    deferred = (
        validate_physical_scope(physical_scope)["deferred"]
        if physical_scope is not None
        else []
    )
    output = {}
    for case_id in PRIMARY:
        bundles = [
            bundle
            for bundle in provenance.records.values()
            if bundle["case_id"] == case_id
        ]
        records = [primary_record(bundle) for bundle in bundles]
        provenance.used.update(
            bundle["references"]["record"]["path"] for bundle in bundles
        )
        item = {
            "status": "pending",
            "functional_records": records,
            "uncached": None,
            "cached": None,
            "catalogue_closure": None,
            "exact_output_and_guard_parity": None,
        }
        if records:
            item["status"] = "functional-only"
        if case_id in deferred:
            item["status"] = "deferred"
            output[case_id] = item
            continue
        candidates = [
            (bundle, record)
            for bundle, record in zip(bundles, records)
            if bundle["window"] is not None
        ]
        if not candidates:
            output[case_id] = item
            continue
        samples, signatures, contexts = {"uncached": {}, "cached": {}}, set(), set()
        for bundle, record in candidates:
            report, launcher = bundle["record"], bundle["launcher"]
            require(
                record["status"] == "functional-closed",
                "primary catalogue is incomplete",
            )
            variant, index = launcher["variant"], launcher["sample"]
            require(
                variant in samples
                and type(index) is int
                and index >= 0
                and index not in samples[variant],
                "duplicate or invalid primary sample",
            )
            qualified = provenance.qualify(
                bundle, report["executable_sha256"], case_id, variant, index
            )
            generator = primary_generator(
                report, launcher, Path(bundle["references"]["launcher"]["path"]).parent
            )
            record["generator_acceptance"] = generator
            require(
                report["source_commit"]
                == qualified["stamps"]["VAKINT_CACHE_BUILD_GIT_REVISION"]
                and report["source_dirty"] is False,
                "primary source/preparation context differs",
            )
            require(
                report["persistent_cache_attachments_before"]
                == {str(loop): False for loop in range(1, 6)}
                and report["persistent_cache_attachments_during_evaluation"]
                == {
                    str(loop): variant == "cached" and loop == 5 for loop in range(1, 6)
                },
                "primary persistent cache policy differs",
            )
            cache = cache_identity(launcher)
            require(
                (variant == "uncached" and cache is None)
                or (variant == "cached" and cache is not None),
                "missing actual manual-cache identity",
            )
            if variant == "cached":
                if isinstance(cache, dict):
                    require(
                        report.get("cache_native_index_identity") == cache
                        and native_cache_identity(report["persistent_cache"]) == cache,
                        "primary native cache changed",
                    )
                else:
                    require(
                        isinstance(cache, str)
                        and len(cache) == 64
                        and digest(report["persistent_cache"]) == cache,
                        "primary cache changed",
                    )
            mode = record["input_mode"]
            require(
                launcher.get("input_mode", "fresh") == mode,
                "primary launcher input mode differs",
            )
            manifest_sha = (
                record["prepared_input"]["manifest_sha256"]
                if mode == "prepared"
                else None
            )
            if mode == "prepared":
                require(
                    launcher.get("prepared_input_manifest_sha256") == manifest_sha,
                    "launcher did not bind the prepared input",
                )
            representation = report.get("coefficient_representation", "sparse")
            require(
                representation in ("sparse", "factorized"),
                "unknown coefficient representation",
            )
            signatures.add(
                json.dumps(
                    {
                        key: report[key]
                        for key in (
                            "projection",
                            "input_sha256",
                            "fixture_sources_sha256",
                            "data_manifest_sha256",
                            "result_sha256",
                            "guards",
                        )
                    },
                    sort_keys=True,
                )
            )
            contexts.add(
                json.dumps(
                    {
                        "stamps": qualified["stamps"],
                        "affinity": qualified["cpu_affinity"],
                        "threads": qualified["thread_environment"],
                        "representation": representation,
                        "input_mode": mode,
                        "generator_acceptance": generator,
                        "prepared_input_manifest_sha256": manifest_sha,
                        "consumer_fixture_sources_sha256": report.get(
                            "consumer_fixture_sources_sha256"
                        ),
                        "artifact": bundle["references"]["artifact"]["sha256"],
                        "window": qualified["window"]["sha256"],
                    },
                    sort_keys=True,
                )
            )
            record["performance_eligible"] = True
            record["performance_evidence"] = qualified
            samples[variant][index] = (report, cache)
        require(
            len(signatures) == 1 and len(contexts) == 1,
            "primary input/guards/result/representation/cohort parity differs",
        )
        require(
            set(samples["uncached"]) == set(samples["cached"])
            and len(samples["uncached"]) >= 3,
            "three matched fresh samples per mode required",
        )
        require(
            len(
                {
                    json.dumps(cache, sort_keys=True)
                    for _, cache in samples["cached"].values()
                }
            )
            == 1,
            "primary samples used different caches",
        )
        first_report = next(iter(samples["uncached"].values()))[0]
        item["generator_acceptance"] = generator
        item["input_mode"] = first_report.get("input_mode", "fresh")
        item["prepared_input"] = (
            first_report["prepared_input"] if item["input_mode"] == "prepared" else None
        )
        for variant, indexed in samples.items():
            reports = [indexed[index][0] for index in sorted(indexed)]
            for row in reports:
                for field in (
                    "data_loading_seconds",
                    "reduction_seconds",
                    "end_to_end_seconds",
                ):
                    number(row[field])
            metric_fields = (
                "preparation_seconds",
                "input_loading_seconds",
                "data_loading_seconds",
                "cache_loading_seconds",
                "reduction_seconds",
                "end_to_end_seconds",
            )
            metrics = {
                key: summarize(row[key] for row in reports)
                if all(row.get(key) is not None for row in reports)
                else None
                for key in metric_fields
            }
            metrics["owner_preparation_seconds"] = (
                summarize(
                    row["owner_preparation"]["preparation_during_reduction_seconds"]
                    for row in reports
                )
                if all(
                    row.get("owner_preparation", {}).get(
                        "preparation_during_reduction_seconds"
                    )
                    is not None
                    for row in reports
                )
                else None
            )
            metrics.update(
                pure_ibp_seconds=None,
                terminal_substitution_seconds=None,
                ibp_excluding_owner_seconds=None,
            )
            stages = {
                field: summarize(
                    sum(number(v) for v in row["stages"][field]) for row in reports
                )
                if all(field in row.get("stages", {}) for row in reports)
                else None
                for field in STAGES
            }
            preparation = {
                field: summarize(row["preparation"][field] for row in reports)
                if all(
                    (row.get("preparation") or {}).get(field) is not None
                    for row in reports
                )
                else None
                for field in (
                    "graph_generation_seconds",
                    "ingress_seconds",
                    "scalar_contraction_seconds",
                    "taylor_seconds",
                    "angular_seconds",
                )
            }
            item[variant] = {
                "metrics": metrics,
                "stages": stages,
                "preparation": preparation,
            }
        selected_cache = next(iter(samples["cached"].values()))[1]
        item["cache_sha256"] = (
            selected_cache if isinstance(selected_cache, str) else None
        )
        if isinstance(selected_cache, dict):
            item["cache_native_index_identity"] = selected_cache
        item["data_manifest_sha256"] = reports[0]["data_manifest_sha256"]
        item.update(
            status="measured",
            catalogue_closure=True,
            exact_output_and_guard_parity=True,
        )
        output[case_id] = item
    return output


WARM_REUSE_CASE = "connected-five-loop-warm-reuse-control"


def warm_reuse_records(provenance):
    """One fresh process contains the first call and its same-engine repeat."""
    result = {
        "status": "pending",
        "functional_records": [],
        "first_call": None,
        "warm_call": None,
        "warm_over_first": None,
        "catalogue_closure": None,
        "pure_ibp_seconds": None,
        "owner_preparation_seconds": None,
        "terminal_substitution_seconds": None,
    }
    samples, contexts, invocations = {}, set(), set()
    for bundle in provenance.records.values():
        if bundle["case_id"] != WARM_REUSE_CASE:
            continue
        report, launcher = bundle["record"], bundle["launcher"]
        directory = Path(bundle["references"]["record"]["path"]).parent
        launch_base = Path(bundle["references"]["launcher"]["path"]).parent

        def payload(reference, name=None, base=directory):
            if name is not None:
                require(reference["file"] == name, "unexpected timed-pair payload")
                reference = {**reference, "path": name}
            path = hashed_file(reference, base)
            require(
                type(reference["bytes"]) is int
                and reference["bytes"] == path.stat().st_size,
                "timed-pair payload length differs",
            )
            return path

        require(
            report.get("schema") == "vakint.five-loop-connected-timed-pair.v1"
            and all(
                report.get(key) is True
                for key in (
                    "completed",
                    "exact_value_parity",
                    "ordered_full_guard_parity",
                    "family_identity_parity",
                )
            )
            and report.get("additional_demo_checks_run") is False,
            "timed pair is incomplete or lacks exact parity",
        )
        input_path, inputs = hashed_record(launcher.get("input_manifest"), launch_base)
        require(input_path == directory / "input.json", "timed pair pins another input")
        require(
            inputs["schema"] == report["schema"]
            and inputs["powers_native_order"] == [1, 2, 1, 1, 1, 1, 1, 1, 1]
            and inputs["mass_squared"] == inputs["numerator"] == "1",
            "connected-control input convention differs",
        )
        dot = payload(inputs["dot"], "input.dot")
        payload(inputs["input"], "input.symbolica")
        require(
            payload(inputs["dot_source"]).read_bytes() == dot.read_bytes(),
            "DOT source differs",
        )
        fixture_path, fixture = hashed_record(
            launcher.get("scientific_fixture"), launch_base
        )
        _, old = hashed_record(
            {
                "path": fixture["source_report"],
                "sha256": fixture["source_report_sha256"],
            },
            fixture_path.parent,
        )
        require(
            fixture["dot_sha256"] == digest(dot)
            and fixture["dot_bytes"] == dot.stat().st_size
            and old[fixture["dot_field"]].encode() == dot.read_bytes()
            and fixture["historical_pinched_loop_blocks"]
            == old["pinched_loop_blocks"]
            == [5]
            and fixture["powers_native_order"]
            == old["pinched_powers"]
            == inputs["powers_native_order"]
            and old["exact_equality"] is True,
            "input is not the accepted genuinely connected five-loop control",
        )
        manifests = inputs["data_manifests"]
        require(
            set(manifests)
            == {"manifest.json", *(f"{n}L/manifest.json" for n in range(1, 6))},
            "timed pair lacks the complete loaded data identity",
        )
        for reference in manifests.values():
            payload(reference)
        payload(inputs["application"])
        require(
            manifests["manifest.json"]["sha256"] == old["data_manifest_sha256"],
            "connected-control data differs from its accepted fixture",
        )
        snapshot_path, snapshot = hashed_record(
            launcher.get("source_snapshot"), launch_base
        )
        stamps = bundle["artifact"]["stamps"]
        bindings = {
            "build_rustc": "VAKINT_CACHE_BUILD_RUSTC",
            "build_profile": "VAKINT_CACHE_BUILD_PROFILE",
            **{
                "build_" + field: stamp
                for field, stamp in (
                    ("git_revision", "VAKINT_CACHE_BUILD_GIT_REVISION"),
                    ("rustred_revision", "VAKINT_CACHE_BUILD_RUSTRED_REVISION"),
                    ("symbolica_revision", "VAKINT_CACHE_BUILD_SYMBOLICA_REVISION"),
                    (
                        "main_source_snapshot_sha256",
                        "VAKINT_CACHE_BUILD_MAIN_SOURCE_SNAPSHOT_SHA256",
                    ),
                    (
                        "native_source_snapshot_sha256",
                        "VAKINT_CACHE_BUILD_NATIVE_SOURCE_SNAPSHOT_SHA256",
                    ),
                )
            },
        }
        require(
            all(
                inputs.get(field) is not None and inputs[field] == stamps.get(stamp)
                for field, stamp in bindings.items()
            )
            and snapshot["head"] == inputs["build_git_revision"]
            and digest(snapshot_path) == inputs["build_main_source_snapshot_sha256"],
            "timed pair belongs to another compiled source/cohort",
        )
        for field, name in (
            (
                "compiled_example_sha256",
                "crates/vakint/examples/five_loop_saved_application.rs",
            ),
            (
                "compiled_dot_adapter_sha256",
                "crates/vakint/examples/saved_application/dot.rs",
            ),
            ("model_sha256", "crates/vakint/tests/fixtures/phi_3_4.json"),
            ("compiled_cargo_lock_sha256", "Cargo.lock"),
        ):
            require(
                inputs[field] == snapshot["source_sha256"][name],
                "compiled pair ingredient differs",
            )
        binary = inputs["executable"]["sha256"]
        payload(inputs["executable"])
        require(
            binary == launcher["binary_sha256"]
            and binary
            in {entry["sha256"] for entry in artifact_binaries(bundle["artifact"])},
            "timed-pair executable differs",
        )
        phases, outputs = {}, {}
        snapshots = [inputs["cache_statistics_before_cold"]]
        for phase in ("cold", "warm"):
            path = payload(report[phase + "_record"], phase + ".json")
            record = read(path)
            require(
                record["phase"] == phase and record["completed"] is True,
                "incomplete timed call",
            )
            number(record["whole_public_application_seconds"])
            value = payload(record["value"], phase + ".symbolica")
            guards = payload(
                record["ordered_nonzero_conditions"], phase + "-guards.json"
            )
            families = payload(record["families"], phase + "-families.json")
            ordered, family_rows = read(guards), read(families)
            require(
                value.stat().st_size > 0
                and isinstance(ordered, list)
                and all(isinstance(item, str) for item in ordered)
                and type(record["nonzero_condition_count"]) is int
                and len(ordered) == record["nonzero_condition_count"]
                and isinstance(family_rows, list)
                and all(
                    isinstance(row, list) and len(row) == 2 and isinstance(row[0], str)
                    for row in family_rows
                )
                and type(record["family_count"]) is int
                and len(family_rows) == record["family_count"],
                "invalid complete guard/family payload",
            )
            outputs[phase] = [digest(path) for path in (value, guards, families)]
            snapshots.append(record["cache_statistics_after"])
            phases[phase] = record
        require(
            outputs["cold"] == outputs["warm"],
            "timed-pair full output/ordered guard/family parity failed",
        )
        for snapshot in snapshots:
            attached = snapshot["persistent"]
            require(
                set(attached) == {str(n) for n in range(1, 6)}
                and all(attached[str(n)]["attached"] is (n < 5) for n in range(1, 6)),
                "connected-control lower-loop/5L cache attachments differ",
            )
        record = {
            "evidence": bundle["references"],
            "input_manifest": evidence(input_path),
            "scientific_fixture": evidence(fixture_path),
            "raw_report": report,
            "input": inputs,
            "calls": phases,
            "performance_eligible": False,
        }
        result["functional_records"].append(record)
        provenance.used.add(bundle["references"]["record"]["path"])
        if bundle["window"] is None:
            continue
        sample = launcher["sample"]
        require(
            type(sample) is int and sample in range(3) and sample not in samples,
            "three distinct pair indices required",
        )
        qualification = provenance.qualify(
            bundle, binary, WARM_REUSE_CASE, "cold-warm-pair", sample
        )
        window = read(qualification["window"]["path"])
        row = next(
            row
            for row in window["samples"]
            if row["report_sha256"] == bundle["references"]["record"]["sha256"]
        )
        for key in ("child_pid", "supervisor_pid", "child_start_ticks"):
            require(
                type(launcher.get(key)) is int and launcher[key] > 0,
                "missing actual pair process identity",
            )
        for key in ("started_unix_seconds", "finished_unix_seconds"):
            number(launcher.get(key))
        require(
            launcher.get("interval_scope")
            == "supervisor pre-spawn setup through supervisor-observed worker completion"
            and all(
                row.get(key) == launcher.get(key)
                for key in (
                    "child_pid",
                    "supervisor_pid",
                    "child_start_ticks",
                    "started_unix_seconds",
                    "finished_unix_seconds",
                    "interval_scope",
                )
            ),
            "pair launcher/window process observations differ",
        )
        invocation = (launcher["child_pid"], launcher["child_start_ticks"])
        require(invocation not in invocations, "warm pair reuses a process invocation")
        invocations.add(invocation)
        contexts.add(
            json.dumps(
                {
                    "artifact": bundle["references"]["artifact"]["sha256"],
                    "window": qualification["window"]["sha256"],
                    "fixture": digest(fixture_path),
                    "snapshot": digest(snapshot_path),
                    "input": inputs["input"]["sha256"],
                    "dot": inputs["dot"]["sha256"],
                    "data": {k: v["sha256"] for k, v in manifests.items()},
                    "application": inputs["application"]["sha256"],
                    "settings": inputs["settings_debug"],
                    "output": outputs["cold"],
                },
                sort_keys=True,
            )
        )
        record["performance_eligible"] = True
        record["performance_evidence"] = qualification
        samples[sample] = phases
    if result["functional_records"]:
        result["status"] = "functional-only"
    if samples:
        require(
            set(samples) == {0, 1, 2} and len(contexts) == 1,
            "three identical-context fresh pairs in one window required",
        )
        result["first_call"] = summarize(
            samples[n]["cold"]["whole_public_application_seconds"] for n in range(3)
        )
        result["warm_call"] = summarize(
            samples[n]["warm"]["whole_public_application_seconds"] for n in range(3)
        )
        first = result["first_call"]["median_seconds"]
        result["warm_over_first"] = (
            result["warm_call"]["median_seconds"] / first if first else None
        )
        result["status"] = "measured-paired-reuse"
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument(
        "--build",
        type=Path,
        help="Completed strict-build or authenticated cumulative-snapshot receipt",
    )
    parser.add_argument(
        "--cache",
        type=Path,
        help="Single-file cache or immutable native cache directory",
    )
    parser.add_argument("--held-out", type=Path)
    parser.add_argument("--tuning", type=Path)
    parser.add_argument(
        "--physical-scope",
        type=Path,
        help="Explicit scope authorization JSON; omitted historical scope requires both examples",
    )
    parser.add_argument(
        "--evidence",
        type=Path,
        required=True,
        help="Hashed report/launcher/artifact/window index; functional records need no window",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--select-current", action="store_true")
    parser.add_argument(
        "--select-tier",
        choices=("heavy", "light"),
        help="Select this authenticated cache tier alongside the other tier; retain the legacy singular selection",
    )
    args = parser.parse_args()
    output = args.output.resolve()
    require(
        not output.exists(),
        "output must be a new path; never overwrite canonical/history files",
    )
    base = read(args.base)
    merged = copy.deepcopy(base)
    require(
        bool(args.build) == bool(args.cache),
        "build and cache must be supplied together",
    )
    require(
        not (args.held_out or args.tuning) or args.build,
        "suite measurements require their authenticated cache receipt",
    )
    provenance = ReportEvidence(args.evidence)
    physical_scope, scope_authorization = load_physical_scope(args.physical_scope)
    require(args.build or provenance.records, "no report evidence supplied")
    build = (
        build_record(args.build, args.cache, provenance, physical_scope)
        if args.build
        else None
    )
    tier = (
        {CUMULATIVE_SCHEMA: "heavy", LIGHT_SCHEMA: "light"}.get(
            build["raw_report"].get("schema")
        )
        if build
        else None
    )
    require(
        args.select_tier is None or args.select_tier == tier,
        "selected tier does not match the actual cache producer",
    )
    run = {
        "schema": "vakint.five-loop-report-run.v1",
        "cache": build,
        "held_out": suite_record(args.held_out, "held-out", build, provenance)
        if args.held_out
        else None,
        "tuning": suite_record(args.tuning, "tuning", build, provenance)
        if args.tuning
        else None,
        "primary": primary_records(provenance, physical_scope),
        "warm_reuse": warm_reuse_records(provenance),
        "evidence_index": evidence(args.evidence),
        "scope": "Functional evidence is retained separately from eligible matched performance samples.",
    }
    if tier is not None:
        run["cache_tier"] = tier
    if scope_authorization is not None:
        run.update(
            physical_scope=physical_scope,
            physical_scope_authorization=scope_authorization,
            notes=[physical_scope_note(scope_authorization)],
        )
    require(
        provenance.used == set(provenance.records), "unused or unknown report evidence"
    )
    for primary in run["primary"].values():
        if primary["status"] == "measured":
            require(
                build is not None,
                "primary performance requires its completed strict cache build",
            )
            if build["raw_report"].get("schema") in (CUMULATIVE_SCHEMA, LIGHT_SCHEMA):
                require(
                    build["status"]
                    in (
                        "validated-committed-snapshot",
                        "validated-strict-selected-file",
                    ),
                    "primary performance needs semantic/import qualification, not prefix preservation",
                )
                acceptance_path, acceptance = hashed_record(
                    provenance.index.get("cache_requalification"),
                    provenance.path.parent,
                )
                _, workload = hashed_record(
                    acceptance["consumer"]["final_workload"], acceptance_path.parent
                )
                _, cohort = hashed_record(
                    acceptance["consumer"]["compiled_cohort"], acceptance_path.parent
                )
                identity = {
                    "cache_native_index_identity": primary.get(
                        "cache_native_index_identity"
                    ),
                    "cache_sha256": primary.get("cache_sha256"),
                    "data_manifest_sha256": primary["data_manifest_sha256"],
                    "workload_sha256": acceptance["consumer"]["final_workload"][
                        "sha256"
                    ],
                    "source_census_sha256": workload["source_census_sha256"],
                    "owner_proof": {
                        **acceptance["consumer"]["owner_proof"],
                        "compiled_cohort": cohort,
                    },
                    "build_git_revision": cohort["main_revision"],
                    "executable_sha256": next(
                        row["raw_report"]["executable_sha256"]
                        for row in primary["functional_records"]
                        if row["performance_eligible"]
                    ),
                }
                primary["cache_requalification"] = provenance.cache_qualification(
                    identity, build, primary=True
                )
                for record in primary["functional_records"]:
                    if (
                        record["performance_eligible"]
                        and record["launcher"]["variant"] == "cached"
                    ):
                        light_attachment(
                            record["raw_report"].get("saved_application_cache"), build
                        )
                require(
                    all(
                        row["evidence"]["artifact"]["sha256"]
                        == acceptance["consumer"]["collector_artifacts"]["sha256"]
                        for row in primary["functional_records"]
                        if row["performance_eligible"]
                    ),
                    "primary artifact differs from accepted native consumer",
                )
            require(
                cache_identity(primary) == cache_identity(build["raw_report"])
                and (
                    primary["data_manifest_sha256"]
                    == build.get("source_producer", build)["raw_report"][
                        "data_manifest_sha256"
                    ]
                    or primary.get("cache_requalification", {}).get(
                        "data_metadata_requalification"
                    )
                    is not None
                ),
                "primary performance cache/data differs from completed build",
            )
    identity_parts = [
        run["evidence_index"],
        build["evidence"] if build else None,
        run["held_out"]["evidence"] if run["held_out"] else None,
        run["tuning"]["evidence"] if run["tuning"] else None,
    ]
    if scope_authorization is not None:
        identity_parts.append(scope_authorization["evidence"])
    identity = json.dumps(identity_parts, sort_keys=True)
    run_id = "five-loop-" + hashlib.sha256(identity.encode()).hexdigest()[:24]
    runs = merged.setdefault("five_loop_runs", {})
    require(run_id not in runs, "this evidence is already appended")
    runs[run_id] = run
    if args.select_current:
        merged["selected_five_loop_run"] = run_id
    if args.select_tier is not None:
        selections = merged.setdefault("selected_five_loop_runs", {})
        require(
            isinstance(selections, dict) and set(selections) <= {"heavy", "light"},
            "invalid existing cache-tier selection",
        )
        for selected_tier, selected_id in selections.items():
            require(
                selected_id in runs
                and runs[selected_id].get("cache_tier") == selected_tier
                and runs[selected_id]
                .get("cache", {})
                .get("raw_report", {})
                .get("schema")
                == {"heavy": CUMULATIVE_SCHEMA, "light": LIGHT_SCHEMA}[selected_tier],
                "existing cache-tier selection names another producer tier",
            )
        selections[args.select_tier] = run_id
        shared_uncached = {}
        for selected_id in selections.values():
            selected_run = runs[selected_id]
            scope = selected_run.get(
                "physical_scope", {"required": list(PRIMARY), "deferred": []}
            )
            for primary in selected_run["primary"].values():
                for record in primary["functional_records"]:
                    if (
                        not record["performance_eligible"]
                        or record["launcher"]["variant"] != "uncached"
                    ):
                        continue
                    refs = record["evidence"]
                    signature = {
                        "launcher": refs["launcher"]["sha256"],
                        "artifact": refs["artifact"]["sha256"],
                        "window": record["performance_evidence"]["window"]["sha256"],
                        "physical_scope": {
                            kind: sorted(cases)
                            for kind, cases in validate_physical_scope(scope).items()
                        },
                    }
                    previous = shared_uncached.setdefault(
                        refs["record"]["sha256"], signature
                    )
                    require(
                        same_json_value(previous, signature),
                        "shared uncached primary must retain identical scope, launcher/cohort and accepted window",
                    )
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as stream:
        json.dump(merged, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(
        json.dumps({"run_id": run_id, "output": str(output), "sha256": digest(output)})
    )


if __name__ == "__main__":
    main()
