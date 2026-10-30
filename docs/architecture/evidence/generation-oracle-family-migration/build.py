"""Fresh source-bound test compilation; no cache deletion/target override."""

import json
import os
import shutil
from pathlib import Path

import common as M
from fingerprints import collect

ROOT = M.OUT / "build-qualified"


def main():
    proc, _ = M.run(
        ["jj", "file", "list", "--ignore-working-copy", "-r", M.APPROVED],
        ROOT / "source-list",
    )
    M.require(proc.returncode == 0, "source", "approved file list")
    listed = [M.WORKSPACE / p for p in proc.stdout.splitlines()]
    files = [p for p in listed if p.is_file()]
    M.write(
        ROOT / "non-file-entries.json",
        {
            "entries": [str(p) for p in listed if not p.is_file()],
            "qualification": "immutable tracked gitlink entries, not compiled regular files; no invented content hashes",
        },
    )
    before = {str(p): M.sha(p) for p in files}
    M.write(ROOT / "sources-before.json", before)
    # Force a genuine compiler invocation with byte-identical source, rather
    # than calling an already fresh cache artifact a newly recorded build.
    os.utime(M.WORKSPACE / M.MAIN, None)
    argv = [
        "cargo",
        "test",
        "--locked",
        "--profile",
        "dev-optim",
        "-p",
        "feynkit-generator",
        "--test",
        "generation_oracle",
        "--no-run",
        "--message-format",
        "json-render-diagnostics",
    ]
    proc, command = M.run(argv, ROOT / "cargo", inputs=files)
    M.require(proc.returncode == 0, "build", "candidate compile")
    messages = []
    for line in proc.stdout.splitlines():
        if line.startswith("{"):
            messages.append(json.loads(line))
    candidates = [
        r
        for r in messages
        if r.get("reason") == "compiler-artifact"
        and r["target"]["name"] == "generation_oracle"
        and r.get("executable")
    ]
    M.require(
        len(candidates) == 1 and not candidates[0]["fresh"],
        "build",
        "non-fresh exact test executable",
    )
    item = candidates[0]
    elf = Path(item["executable"])
    destination = M.OUT / "bin/generation_oracle"
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(elf, destination)
    fp = (
        elf.parent.parent
        / ".fingerprint"
        / ("feynkit-generator-" + elf.name.rsplit("-", 1)[1])
        / "test-integration-test-generation_oracle.json"
    )
    data = collect(fp.parent.parent, fp)
    data.update(
        artifact=str(destination),
        artifact_sha256=M.sha(destination),
        emitted_artifact=item,
        build_command=command,
    )
    M.write(M.OUT / "provenance/new.json", data)
    old = M.load(M.OUT / "provenance/old.json")
    projection = lambda records: sorted(
        (
            r["unit"],
            r["kind"],
            r["features"],
            r["profile"],
            r["rustc"],
            str(r["rustflags"]),
        )
        for r in records
    )
    M.require(
        projection(old["closure"]) == projection(data["closure"]),
        "features",
        "complete dependency feature/profile/compiler closure",
    )
    for key in (
        "rustc",
        "profile",
        "features",
        "rustflags",
        "compile_kind",
        "target",
        "config",
    ):
        M.require(old["root"][key] == data["root"][key], "features", "root " + key)
    after = {str(p): M.sha(p) for p in files}
    M.require(before == after, "source", "source changed during compile")
    M.write(ROOT / "sources-after.json", after)
    for label, args in (("rustc", ["rustc", "-vV"]), ("cargo", ["cargo", "--version"])):
        p, _ = M.run(args, ROOT / "toolchain" / label)
        M.require(p.returncode == 0, "toolchain", label)
    print(
        "retained",
        destination,
        data["artifact_sha256"],
        "dependency units",
        data["closure_units"],
    )


if __name__ == "__main__":
    main()
