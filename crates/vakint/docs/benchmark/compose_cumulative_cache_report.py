"""Compose existing preparation or opaque-union evidence without Native imports.

The finite preparation or Native union and its original source receipts must
already exist. This records their authenticated relationship using the report
adapter's source validators. It creates no new scientific authority.
"""

import argparse
import json
import os
import re
import tempfile
from pathlib import Path

import merge_five_loop_results as reports


def compose(source_path, cache, output_directory):
    source_path, cache = Path(source_path).resolve(), Path(cache).resolve()
    output = Path(output_directory).absolute()
    reports.require(not output.exists(), "report output already exists")
    source = reports.read(source_path)
    verified = reports.strict_preparation_source(source, source_path.parent)
    producer = verified["receipts"]["producer"]["raw_receipt"]
    inventory = verified["receipts"]["inventory"]["raw_receipt"]
    native = verified["receipts"]["native_validation"]["raw_receipt"]
    identity = reports.native_cache_identity(cache)
    reports.require(
        reports.same_json_value(identity, inventory["cache_native_index_identity"]),
        "report source differs from the actual preserved Native indexes",
    )
    reports.require(
        all(
            not ((cache / name).stat().st_mode & 0o222) for name in identity["indexes"]
        ),
        "finite source indexes must be frozen read-only copies",
    )
    # Rebase explicit evidence references, never the historical paths embedded
    # inside immutable runtime receipts. No payload or coefficient is read.
    for name, document in verified["receipts"].items():
        source[name] = {
            "path": os.path.relpath(document["evidence"]["path"], output),
            "sha256": document["evidence"]["sha256"],
        }
    report = {
        "schema": reports.CUMULATIVE_SCHEMA,
        "cache_sha256": None,
        "cache_native_index_identity": identity,
        "semantic_binding": native["semantic_binding"],
        "data_manifest_sha256": producer["data_manifest_sha256"],
        "sources": [source],
        "inventory": verified["inventory_counts"],
        "scope": "One completed finite preparation, not final HEAVY/LIGHT or a completed physical example.",
    }
    validation = {
        "schema": "vakint.cumulative-native-cache-validation.v1",
        "accepted": True,
        **{
            key: report[key]
            for key in (
                "cache_native_index_identity",
                "semantic_binding",
                "data_manifest_sha256",
                "sources",
                "inventory",
            )
        },
        "checks": {
            key: True
            for key in (
                "immutable_index_snapshot",
                "authenticated_chunks",
                "strict_completed_records_only",
                "authenticated_bindings",
                "committed_source_snapshots",
            )
        },
        "authority": "Composition of the referenced actual preservation and full Native-validation receipts; no new Native run.",
    }
    return write_report(report, validation, cache, output)


def write_report(report, validation, cache, output):
    # Publish only after the ordinary consumer accepts every original source.
    with tempfile.TemporaryDirectory(
        dir=output.parent, prefix=f".{output.name}-"
    ) as temporary:
        staging = Path(temporary)
        with (staging / "native-validation.json").open("x") as stream:
            json.dump(validation, stream, indent=2)
            stream.write("\n")
        report["validation"] = {
            "path": "native-validation.json",
            "sha256": reports.digest(staging / "native-validation.json"),
        }
        with (staging / "cumulative-report.json").open("x") as stream:
            json.dump(report, stream, indent=2)
            stream.write("\n")
        result = reports.cumulative_cache_record(
            staging / "cumulative-report.json", cache
        )
        reports.require(not output.exists(), "report output already exists")
        staging.rename(output)
        result["evidence"] = reports.evidence(output / "cumulative-report.json")
        result["validation"]["evidence"] = reports.evidence(
            output / "native-validation.json"
        )
        return result


def compose_union(
    source_paths,
    union_validation_path,
    cache,
    output_directory,
    *,
    data_manifest,
    proof_root,
):
    """Record an actual ordered Native union; completion remains source-owned."""
    output = Path(output_directory).absolute()
    reports.require(not output.exists(), "report output already exists")
    raw_path = Path(union_validation_path).resolve()
    raw_reference = reports.evidence(raw_path)
    _, native = reports.ReportEvidence.strict_record(raw_reference, raw_path.parent)
    proof_root = Path(proof_root).resolve()
    reports.require(proof_root.is_dir(), "union needs an explicit proof root")
    cache = Path(cache)
    reports.require(
        not cache.is_symlink()
        and not (raw_path.parent / "cache").is_symlink()
        and cache.resolve() == raw_path.parent / "cache",
        "union cache differs from the Native output directory",
    )
    cache = cache.resolve()

    def unsigned(value):
        return type(value) is int and 0 <= value < 2**64

    def sha256(value):
        return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value)

    fields = {
        "generation",
        "runs",
        "chunk_references",
        "stored_record_references",
        "referenced_compressed_bytes",
        "manifest_bytes",
    }
    summed = fields - {"generation", "manifest_bytes"}

    def census(value):
        return (
            isinstance(value, dict)
            and set(value) == fields
            and all(unsigned(count) for count in value.values())
            and value["generation"] > 0
            and value["manifest_bytes"] > 0
        )

    reports.require(
        native.get("schema") == "vakint.native-opaque-union-validation.v1"
        and native.get("success") is True
        and sha256(native.get("request_sha256"))
        and native.get("source_indexes_unchanged") is True
        and native.get("frozen_readers_used") is True
        and "distinct_keys" in native
        and native["distinct_keys"] is None
        and type(native.get("compiled_codec_state_format")) is int
        and native["compiled_codec_state_format"] == 6
        and native.get("coefficient_context_created") is False
        and all(
            type(native.get(key)) is int and native[key] == 0
            for key in ("state_or_coefficient_imports", "IBP_applications")
        ),
        "union lacks an actual successful opaque Native validation",
    )
    binding = native["semantic_binding"]
    reports.require(
        isinstance(binding, dict)
        and set(binding)
        == {"artifact_digest", "family_fingerprint", "physical_arity", "layer"}
        and isinstance(binding["artifact_digest"], list)
        and len(binding["artifact_digest"]) == 32
        and all(
            type(byte) is int and 0 <= byte <= 255
            for byte in binding["artifact_digest"]
        )
        and isinstance(binding["family_fingerprint"], str)
        and binding["family_fingerprint"]
        and type(binding["physical_arity"]) is int
        and binding["physical_arity"] == 15
        and binding["layer"] == "saved-application-route-apply-v1",
        "union semantic binding is malformed",
    )
    _, artifact = reports.ReportEvidence.strict_record(
        native["collector_artifacts"], proof_root
    )
    reports.validate_compiled_artifact(native, artifact)
    reports.require(
        reports.same_json_value(artifact["compiled_cohort"], native["compiled_cohort"])
        and reports.same_json_value(
            artifact["stamps"], native["compiled_cohort"]["stamps"]
        )
        and native["executable_sha256"]
        == native["compiled_cohort"]["measurement_executable_sha256"],
        "union producer differs from its own compiled artifact",
    )
    admission = native["resource_guard_admission"]
    reports.require(
        admission.get("schema") == "vakint.resource-guard-admission.v1"
        and admission.get("admitted") is True
        and sha256(admission.get("launch_spec_sha256"))
        and sha256(admission.get("guard_config_sha256"))
        and admission["worker"]["executable_sha256"] == native["executable_sha256"]
        and all(
            unsigned(admission[who][key]) and admission[who][key] > 0
            for who in ("worker", "guard")
            for key in ("pid", "start_ticks")
        )
        and unsigned(admission["worker"]["uid"])
        and reports.same_json_value(
            admission["guard"]["uid"], admission["worker"]["uid"]
        )
        and isinstance(admission.get("slot"), str)
        and 0 < len(admission["slot"]) <= 256
        and isinstance(admission.get("admission_utc"), str)
        and 0 < len(admission["admission_utc"]) <= 128
        and unsigned(admission.get("MemAvailable_bytes")),
        "union lacks its own successful resource admission",
    )
    source_paths = [Path(path).resolve() for path in source_paths]
    sources, source_references = [], []
    for path in source_paths:
        reference = reports.evidence(path)
        _, source = reports.ReportEvidence.strict_record(reference, path.parent)
        # Only transport explicit outer references; historical receipt paths stay exact.
        for key, value in source.items():
            if isinstance(value, dict) and "path" in value and "sha256" in value:
                source[key] = reports.evidence(reports.hashed_file(value, path.parent))
        sources.append(source)
        source_references.append(reference)
    raw_sources = native["sources_newest_first"]
    snapshots, censuses = native["captured_snapshots"], native["source_censuses"]
    reports.require(
        sources
        and len(sources) == len(raw_sources) == len(snapshots) == len(censuses)
        and [source["id"] for source in sources]
        == [source["id"] for source in raw_sources]
        and len({source["id"] for source in sources}) == len(sources)
        and all(
            isinstance(source["id"], str)
            and re.fullmatch(r"[A-Za-z0-9_-]+", source["id"])
            for source in sources
        ),
        "union source IDs/order differ from the logical source inventory",
    )
    directories, identities = [], set()
    for position, (source, raw_source, frozen, actual) in enumerate(
        zip(sources, raw_sources, snapshots, censuses)
    ):
        identity = reports.validated_native_cache_identity(
            raw_source["cache_native_index_identity"]
        )
        directory = Path(raw_source["directory"])
        reports.require(
            not directory.is_relative_to(raw_path.parent)
            and not raw_path.parent.is_relative_to(directory),
            "union source overlaps the Native output directory",
        )
        reports.require(
            directory.is_absolute()
            and directory.resolve() == directory
            and not any(
                directory.is_relative_to(previous) or previous.is_relative_to(directory)
                for previous in directories
            )
            and identity["sha256"] not in identities
            and raw_source["kind"] == identity["layout"]
            and reports.same_json_value(
                reports.native_cache_identity(directory), identity
            )
            and (identity["layout"] == "pack" or len(identity["indexes"]) == 2),
            "union source directory or frozen Native index identity differs",
        )
        directories.append(directory)
        identities.add(identity["sha256"])
        child = next(name for name in identity["indexes"] if name != "bundle.rrindex")
        expected = raw_source["expected"]
        nullable = {"runs", "chunk_references", "referenced_compressed_bytes"}
        reports.require(
            census(actual)
            and isinstance(expected, dict)
            and set(expected) == fields
            and all(
                unsigned(count) or key in nullable and count is None
                for key, count in expected.items()
            )
            and all(
                count is None or count == actual[key] for key, count in expected.items()
            )
            and frozen["manifest"] == f"snapshots/source-{position}.rrindex"
            and frozen["refresh_allowed"] is False
            and type(frozen["generation"]) is int
            and frozen["generation"] == actual["generation"]
            and type(frozen["manifest_bytes"]) is int
            and frozen["manifest_bytes"] == actual["manifest_bytes"]
            and frozen["manifest_bytes"] == identity["indexes"][child]["bytes"]
            and frozen["manifest_sha256"] == identity["indexes"][child]["sha256"]
            and isinstance(frozen["manifest_blake3"], list)
            and len(frozen["manifest_blake3"]) == 32
            and all(
                type(byte) is int and 0 <= byte <= 255
                for byte in frozen["manifest_blake3"]
            ),
            "union frozen source snapshot/census differs",
        )
        snapshot_path = raw_path.parent / frozen["manifest"]
        reports.require(
            not snapshot_path.parent.is_symlink()
            and snapshot_path.is_file()
            and not snapshot_path.is_symlink()
            and snapshot_path.stat().st_size == frozen["manifest_bytes"]
            and reports.digest(snapshot_path) == frozen["manifest_sha256"],
            "union captured snapshot bytes changed",
        )
        if source["kind"] == "strict-preparation":
            _, inventory = reports.hashed_record(source["inventory"], output)
            _, preserved = reports.hashed_record(source["preservation"], output)
            reports.require(
                reports.same_json_value(
                    inventory["cache_native_index_identity"], identity
                )
                and reports.same_json_value(
                    preserved["generation"], actual["generation"]
                )
                and reports.same_json_value(
                    preserved["stored_record_references"],
                    actual["stored_record_references"],
                ),
                "union finite source differs from its preserved original snapshot",
            )
        elif source["kind"] == "native-index-prefix":
            preserved = reports.preserved_native_source(source, output)
            summary = preserved["receipts"].get("native_summary")
            reports.require(
                summary is not None,
                "union prefix lacks its authenticated Native per-role census",
            )
            summary = summary["raw_receipt"]
            role = source["prefix_role"]
            index = reports.hashed_file(source["snapshot"], output)
            reports.require(
                reports.digest(index) == frozen["manifest_sha256"]
                and index.stat().st_size == frozen["manifest_bytes"]
                and reports.same_json_value(
                    {key: summary[f"{role}_census"][key] for key in fields}, actual
                )
                and reports.same_json_value(
                    summary[f"{role}_snapshot"]["manifest_digest"],
                    frozen["manifest_blake3"],
                )
                and reports.same_json_value(
                    summary["binding"], native["semantic_binding"]
                ),
                "union prefix differs from its authenticated per-role census/snapshot/binding",
            )
        else:
            reports.require(
                source["kind"] == "running-capture", "unknown union logical source kind"
            )
            index = reports.hashed_file(source["index"], output)
            _, descriptor = reports.hashed_record(source["snapshot"], output)
            _, capture = reports.hashed_record(source["capture"], output)
            captured_native = capture.get("report", capture.get("native", capture))
            reports.require(
                reports.digest(index) == frozen["manifest_sha256"]
                and index.stat().st_size == frozen["manifest_bytes"]
                and reports.same_json_value(
                    descriptor["delta"]["generation"], actual["generation"]
                )
                and descriptor["delta"]["manifest_digest"]
                == bytes(frozen["manifest_blake3"]).hex()
                and reports.same_json_value(
                    captured_native["committed_records"],
                    actual["stored_record_references"],
                ),
                "union captured source differs from its original snapshot",
            )
    totals = {key: sum(census_row[key] for census_row in censuses) for key in summed}
    destination, copied = native["output_census"], native["copy_statistics"]
    copy_fields = {
        "namespaces",
        "copied_bases",
        "unique_physical_chunks",
        "physical_compressed_bytes",
        "index_bytes",
        "preserved_record_references",
        "transfer_buffer_bytes",
        "charged_disk_bytes",
    }
    identity = reports.native_cache_identity(cache)
    reports.require(
        census(destination)
        and destination["generation"] == 1
        and all(unsigned(value) for value in totals.values())
        and reports.same_json_value(native["source_reference_totals"], totals)
        and all(destination[key] == value for key, value in totals.items())
        and isinstance(copied, dict)
        and set(copied) == copy_fields
        and all(unsigned(value) for value in copied.values())
        and copied["namespaces"] == 1
        and copied["copied_bases"] == len(sources)
        and copied["preserved_record_references"] == totals["stored_record_references"]
        and copied["unique_physical_chunks"] <= totals["chunk_references"]
        and copied["physical_compressed_bytes"] <= totals["referenced_compressed_bytes"]
        and copied["charged_disk_bytes"]
        == copied["physical_compressed_bytes"] + copied["index_bytes"]
        and copied["index_bytes"]
        == sum(row["bytes"] for row in identity["indexes"].values())
        and reports.same_json_value(identity, native["cache_native_index_identity"])
        and identity["layout"] == "bundle"
        and len(identity["indexes"]) == 2
        and destination["manifest_bytes"]
        == next(
            row["bytes"]
            for name, row in identity["indexes"].items()
            if name != "bundle.rrindex"
        ),
        "union output census/copy arithmetic or actual indexes differ",
    )
    limits = native["limits"]
    reports.require(
        set(limits)
        == {"max_disk_bytes", "max_manifest_bytes", "max_directory_owned_bytes"}
        and all(unsigned(value) and value > 0 for value in limits.values())
        and sum(
            source["expected"]["referenced_compressed_bytes"] or 0
            for source in raw_sources
        )
        <= limits["max_disk_bytes"]
        and copied["charged_disk_bytes"] <= limits["max_disk_bytes"]
        and all(
            row["manifest_bytes"] <= limits["max_manifest_bytes"]
            for row in [*censuses, destination]
        ),
        "union copy exceeds its declared Native grants",
    )
    manifest_path = Path(data_manifest).resolve()
    manifest_reference = reports.evidence(manifest_path)
    _, manifest = reports.ReportEvidence.strict_record(
        manifest_reference, manifest_path.parent
    )
    reports.require(
        manifest.get("schema") == "vakint.external-data.v1",
        "unknown union report data manifest",
    )
    report = {
        "schema": reports.CUMULATIVE_SCHEMA,
        "cache_sha256": None,
        "cache_native_index_identity": identity,
        "semantic_binding": binding,
        "data_manifest_sha256": manifest_reference["sha256"],
        "sources": sources,
        "inventory": {
            "compressed_bytes": copied["physical_compressed_bytes"],
            "stored_record_references": totals["stored_record_references"],
            "physical_chunks": copied["unique_physical_chunks"],
        },
        "scope": "Opaque committed-source union; no distinct-key, physical-completion, selection or performance qualification.",
    }
    validation = {
        "schema": "vakint.cumulative-native-cache-validation.v1",
        "accepted": True,
        **{
            key: report[key]
            for key in (
                "cache_native_index_identity",
                "semantic_binding",
                "data_manifest_sha256",
                "sources",
                "inventory",
            )
        },
        "native_union_validation": raw_reference,
        "source_definitions": source_references,
        "report_data_manifest": manifest_reference,
        "proof_root": str(proof_root),
        "checks": {
            key: True
            for key in (
                "immutable_index_snapshot",
                "authenticated_chunks",
                "strict_completed_records_only",
                "authenticated_bindings",
                "committed_source_snapshots",
            )
        },
        "authority": "Actual Native opaque-union validation plus the cumulative adapter's original-source validators; no new Native run or source-exposure qualification.",
    }
    return write_report(report, validation, cache, output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, action="append", required=True)
    parser.add_argument("--union-validation", type=Path)
    parser.add_argument("--data-manifest", type=Path)
    parser.add_argument("--proof-root", type=Path)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    options = parser.parse_args()
    if options.union_validation is not None:
        if options.data_manifest is None or options.proof_root is None:
            parser.error("union composition requires --data-manifest and --proof-root")
        compose_union(
            options.source,
            options.union_validation,
            options.cache,
            options.output_directory,
            data_manifest=options.data_manifest,
            proof_root=options.proof_root,
        )
    else:
        if (
            len(options.source) != 1
            or options.data_manifest is not None
            or options.proof_root is not None
        ):
            parser.error(
                "multiple sources and union context require --union-validation"
            )
        compose(options.source[0], options.cache, options.output_directory)
