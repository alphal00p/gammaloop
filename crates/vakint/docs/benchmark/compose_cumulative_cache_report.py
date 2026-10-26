"""Compose existing strict-preparation evidence; never copy or load Native data.

The finite producer, preserved checkpoint and full guarded Native source-parity
run must already exist. This records their authenticated relationship using the
report adapter's source validator. It creates no new scientific authority.
"""

import argparse
import json
import os
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
    output.mkdir(mode=0o700)
    with (output / "native-validation.json").open("x") as stream:
        json.dump(validation, stream, indent=2)
        stream.write("\n")
    report["validation"] = {
        "path": "native-validation.json",
        "sha256": reports.digest(output / "native-validation.json"),
    }
    with (output / "cumulative-report.json").open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    # Exercise the same cumulative consumer used by benchmark reporting.
    return reports.cumulative_cache_record(output / "cumulative-report.json", cache)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    options = parser.parse_args()
    compose(options.source, options.cache, options.output_directory)
