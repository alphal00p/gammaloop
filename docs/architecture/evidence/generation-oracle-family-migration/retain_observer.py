"""Retain the emitted observer ELF, source overlay and real fingerprint closure."""

import json
import shutil
from pathlib import Path

import common as M
from fingerprints import collect

matches = list((M.OUT / "commands").glob("*-observer-build-bounded-shell/stdout.txt"))
M.require(len(matches) == 1, "observer", "actual build transcript")
messages = [
    json.loads(line)
    for line in matches[0].read_text().splitlines()
    if line.startswith("{")
]
items = [
    r
    for r in messages
    if r.get("reason") == "compiler-artifact"
    and r["target"]["name"] == "oracle_family_migration_capture"
    and r.get("executable")
]
M.require(len(items) == 1, "observer", "emitted ELF")
item = items[0]
elf = Path(item["executable"])
destination = M.OUT / "bin/oracle_family_migration_capture"
shutil.copy2(elf, destination)
found = []
for candidate in Path("/common/build-cache/lcnbr").glob(
    "*/*/dev-optim/.fingerprint/feynkit-generator-*/example-oracle_family_migration_capture.json"
):
    build_dir = candidate.parent.parent.parent
    meta = candidate.parent.name.rsplit("-", 1)[1]
    for binary in (
        build_dir / "examples/oracle_family_migration_capture",
        build_dir / ("examples/oracle_family_migration_capture-" + meta),
    ):
        if binary.is_file() and M.sha(binary) == M.sha(destination):
            found.append(candidate)
            break
M.require(
    len(found) == 1,
    "observer",
    "unique normal-cache fingerprint with byte-identical emitted ELF",
)
fp = found[0]
data = collect(fp.parent.parent, fp)
source = (
    M.OUT
    / "observer-source/crates/feynkit-generator/examples/oracle_family_migration_capture.rs"
)
shutil.copy2(source, M.HERE / "observer-source.rs.txt")
data.update(
    artifact=str(destination),
    artifact_sha256=M.sha(destination),
    emitted_artifact=item,
    build_transcript=str(matches[0]),
    build_transcript_sha256=M.sha(matches[0]),
    source_sha256=M.sha(source),
    source=str(source),
    source_note="private approved-base copy + pool argument observer overlay only; main.rs not compiled by observer; linked through complete actual test dumps",
)
M.write(M.OUT / "provenance/observer.json", data)
print(data["artifact_sha256"], data["closure_units"])
