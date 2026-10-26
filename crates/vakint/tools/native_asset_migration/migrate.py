#!/usr/bin/env python3
"""Transport Vakint's 15 trusted native assets from Symbolica state 6 to 11."""

import argparse
import gzip
import hashlib
import json
import shutil
import struct
import subprocess
from pathlib import Path

OLD_REVISION = "a19c760dd567c239f30d87e4e924ca2f8b8457ab"
NEW_REVISION = "f514a7fe1b19e851207f9421ef9a4a4157d5f0e5"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def sections(data):
    assert data[:8] == b"RRPBIN\r\n"
    assert struct.unpack_from("<I", data, 8)[0] == 1
    offset = 20
    result = {}
    for _ in range(struct.unpack_from("<I", data, 16)[0]):
        tag, flags, length = struct.unpack_from("<HHQ", data, offset)
        assert flags == 0 and tag not in result
        offset += 12
        result[tag] = data[offset : offset + length]
        assert len(result[tag]) == length
        offset += length
    assert offset == len(data)
    return result


def record_digest(path):
    checksum = hashlib.sha256()
    count = 0
    with gzip.open(path, "rb") as stream:
        for line in stream:
            checksum.update(line)
            count += 1
    return checksum.hexdigest(), count


def build_bridge(source, label, work, bridge):
    revision = subprocess.check_output(
        ["git", "-C", str(source), "rev-parse", "HEAD"], text=True
    ).strip()
    assert revision == {"old": OLD_REVISION, "new": NEW_REVISION}[label], revision
    package = work / label
    package.mkdir(parents=True, exist_ok=True)
    symbolica = json.dumps(str(source.resolve()))
    numerica = json.dumps(str((source / "lib/numerica").resolve()))
    graphica = json.dumps(str((source / "lib/graphica").resolve()))
    rust = json.dumps(str(bridge.resolve()))
    (package / "Cargo.toml").write_text(f"""[package]
name = "vakint-state-{label}"
version = "0.1.0"
edition = "2024"

[[bin]]
name = "vakint-state-{label}"
path = {rust}

[dependencies]
symbolica = {{ path = {symbolica}, default-features = false, features = ["integer-gmp", "float-mpfr"] }}
serde = {{ version = "1", features = ["derive"] }}
serde_json = "1"
flate2 = "1"

[patch.crates-io]
numerica = {{ path = {numerica} }}
graphica = {{ path = {graphica} }}

[profile.dev]
opt-level = 1
debug = 0
""")
    subprocess.run(
        [
            "cargo",
            "build",
            "--manifest-path",
            str(package / "Cargo.toml"),
            "--target-dir",
            str(package / "target"),
            "-j",
            "4",
        ],
        check=True,
    )
    return package / "target/debug" / f"vakint-state-{label}"


def migrate(source, destination, work, old, new):
    previous_receipt = (source / "codec-migration-2026-10-05.json").read_bytes()
    previous_assets = {
        entry["path"]: entry for entry in json.loads(previous_receipt)["assets"]
    }
    files = [source / f"unit_mass_vacuum_k{k}.rrbin" for k in [1, 3, 6]]
    files += [
        source / "four_loop" / f"{family}.{suffix}"
        for family in ["h", "x", "bmw", "fg"]
        for suffix in ["candidates.rrbin.gz", "rrcat.bin.gz", "rrnorm.bin"]
    ]
    destination.mkdir(parents=True, exist_ok=False)
    report = {
        "schema": "vakint-native-asset-codec-migration.v2",
        "date": "2026-10-26",
        "symbolica_revision": {"decoder": OLD_REVISION, "encoder": NEW_REVISION},
        "symbolica_export_format": {"old": 6, "new": 11},
        "source_receipt_sha256": digest(previous_receipt),
        "method": "native_decoder_exact_sparse_polynomials_and_canonical_text_native_encoder",
        "rule_generation": False,
        "master_value_recomputation": False,
        "candidate_closure_claim": False,
        "numerical_acceptance_claim": False,
        "assets": [],
    }
    for original in files:
        relative = original.relative_to(source)
        assert (
            digest(original.read_bytes())
            == previous_assets[str(relative)]["new_sha256"]
        ), relative
        prefix = work / str(relative).replace("/", "_")
        before_records = prefix.with_name(prefix.name + ".before.jsonl.gz")
        after_records = prefix.with_name(prefix.name + ".after.jsonl.gz")
        native_output = prefix.with_name(prefix.name + ".native")
        subprocess.run([old, "decode", original, before_records], check=True)
        subprocess.run(
            [new, "encode", original, before_records, native_output], check=True
        )
        subprocess.run([new, "decode", native_output, after_records], check=True)
        records_sha, count = record_digest(before_records)
        assert (records_sha, count) == record_digest(after_records), relative
        stored_before = original.read_bytes()
        native_before = (
            gzip.decompress(stored_before)
            if original.suffix == ".gz"
            else stored_before
        )
        native_after = native_output.read_bytes()
        old_sections, new_sections = sections(native_before), sections(native_after)
        assert native_before[:20] == native_after[:20]
        assert old_sections.keys() == new_sections.keys()
        assert struct.unpack_from("<H", old_sections[1], 4)[0] == 6
        assert struct.unpack_from("<H", new_sections[1], 4)[0] == 11
        structural = [tag for tag in old_sections if tag not in [1, 2, 7]]
        assert all(old_sections[tag] == new_sections[tag] for tag in structural), (
            relative
        )
        output = destination / relative
        output.parent.mkdir(parents=True, exist_ok=True)
        if original.suffix == ".gz":
            with output.open("wb") as stream:
                subprocess.run(
                    ["gzip", "-n", "-9", "-c", native_output], stdout=stream, check=True
                )
        else:
            shutil.copyfile(native_output, output)
        stored_after = output.read_bytes()
        report["assets"].append(
            {
                "path": str(relative),
                "frames": count,
                "exact_portable_records_sha256": records_sha,
                "structural_sections_unchanged": structural,
                "old_sha256": digest(stored_before),
                "new_sha256": digest(stored_after),
                "old_stored_bytes": len(stored_before),
                "new_stored_bytes": len(stored_after),
                "old_native_sha256": digest(native_before),
                "new_native_sha256": digest(native_after),
                "old_native_bytes": len(native_before),
                "new_native_bytes": len(native_after),
            }
        )
    (destination / "codec-migration-2026-10-26.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("old_symbolica", type=Path)
    parser.add_argument("new_symbolica", type=Path)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("work", type=Path)
    args = parser.parse_args()
    assert not args.destination.exists(), "destination must be new"
    args.work.mkdir(parents=True, exist_ok=True)
    bridge = Path(__file__).with_name("bridge.rs")
    old = build_bridge(args.old_symbolica, "old", args.work, bridge)
    new = build_bridge(args.new_symbolica, "new", args.work, bridge)
    report = migrate(args.source, args.destination, args.work, old, new)
    print(f"Verified {len(report['assets'])} migrated native assets")


if __name__ == "__main__":
    main()
