#!/usr/bin/env python3
"""Refresh the offline 1L–3L snapshot from a validated external-data checkout."""

import argparse
import hashlib
import json
import shutil
import tempfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument(
        "--output", type=Path, default=Path(__file__).resolve().parents[1] / "data"
    )
    args = parser.parse_args()
    release = json.loads((args.source / "manifest.json").read_text())
    if release["schema"] != "vakint.external-data.v1":
        raise ValueError("unsupported external-data schema")
    packs = {pack["loop_count"]: pack for pack in release["loop_packs"]}
    release["loop_packs"] = [packs[n] for n in (1, 2, 3)]
    with tempfile.TemporaryDirectory(prefix="vakint-vendored-") as scratch:
        staging = Path(scratch)
        for pack in release["loop_packs"]:
            n = pack["loop_count"]
            payload = (args.source / pack["manifest"]).read_bytes()
            if hashlib.sha256(payload).hexdigest() != pack["manifest_sha256"]:
                raise ValueError(f"{n}L manifest checksum mismatch")
            manifest = json.loads(payload)
            manifest["cache"] = {"available": False}
            manifest["files"] = {
                name: digest
                for name, digest in manifest["files"].items()
                if name != "cache.json" and not name.startswith("cache/")
            }
            for name, digest in manifest["files"].items():
                relative = Path(name)
                if relative.is_absolute() or ".." in relative.parts:
                    raise ValueError(f"invalid relative path {name}")
                source = args.source / f"{n}L" / relative
                contents = source.read_bytes()
                if (
                    len(contents) != digest["bytes"]
                    or hashlib.sha256(contents).hexdigest() != digest["sha256"]
                ):
                    raise ValueError(f"{n}L/{name} checksum mismatch")
                destination = staging / f"{n}L" / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(contents)
            payload = (
                json.dumps(manifest, ensure_ascii=False, indent=2) + "\n"
            ).encode()
            (staging / pack["manifest"]).write_bytes(payload)
            pack.update(
                manifest_sha256=hashlib.sha256(payload).hexdigest(),
                archive_sha256=None,
                archive_bytes=None,
            )
        (staging / "manifest.json").write_text(
            json.dumps(release, ensure_ascii=False, indent=2) + "\n"
        )
        args.output.mkdir(parents=True, exist_ok=True)
        for n in (1, 2, 3):
            destination = args.output / f"{n}L"
            if destination.exists():
                shutil.rmtree(destination)
            shutil.copytree(staging / f"{n}L", destination)
        shutil.copyfile(staging / "manifest.json", args.output / "manifest.json")


if __name__ == "__main__":
    main()
