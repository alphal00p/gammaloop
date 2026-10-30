"""New pristine assertion and bounded observer builds in private source copies."""

import argparse
import json
import shutil
from pathlib import Path

import common as M
from fingerprints import collect


def fingerprints(elf, name, kind):
    hits = []
    for fp in Path("/common/build-cache/lcnbr").glob(
        f"*/*/dev-optim/.fingerprint/feynkit-generator-*/{kind}-{name}.json"
    ):
        meta = fp.parent.name.rsplit("-", 1)[1]
        base = fp.parent.parent.parent
        candidate = (
            base
            / ("deps" if kind.startswith("test") else "examples")
            / (name + "-" + meta)
        )
        if candidate.is_file() and M.sha(candidate) == M.sha(elf):
            hits.append(fp)
    M.require(len(hits) == 1, "fingerprints", f"{name}: unique emitted/cache artifact")
    return collect(hits[0].parent.parent, hits[0])


def compile_one(route, name, kind):
    root = M.OUT / "source" / route
    test = kind == "test-integration-test"
    argv = [
        "cargo",
        "test" if test else "build",
        "--locked",
        "--profile",
        "dev-optim",
        "-p",
        "feynkit-generator",
        "--test" if test else "--example",
        name,
        *(["--no-run"] if test else []),
        "--message-format",
        "json-render-diagnostics",
    ]
    sources = [
        p
        for p in root.rglob("*")
        if p.is_file() and "target" not in p.relative_to(root).parts
    ]
    before = {str(p): M.sha(p) for p in sources}
    out = M.OUT / "build" / route / name
    M.write(out / "sources-before.json", before)
    proc, command = M.run(argv, out / "execution", cwd=root, inputs=sources)
    M.require(proc.returncode == 0, "build", route + "/" + name)
    messages = [
        json.loads(line) for line in proc.stdout.splitlines() if line.startswith("{")
    ]
    rows = [
        r
        for r in messages
        if r.get("reason") == "compiler-artifact"
        and r["target"]["name"] == name
        and r.get("executable")
    ]
    M.require(len(rows) == 1, "build", "actual executable")
    elf = Path(rows[0]["executable"])
    dest = M.OUT / "bin" / route / name
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(elf, dest)
    data = fingerprints(dest, name, kind)
    data.update(
        artifact=str(dest),
        artifact_sha256=M.sha(dest),
        emitted=rows[0],
        command=command,
        source_revision=M.REVISIONS[route],
        sources_before_after=before,
    )
    M.write(out / "qualification.json", data)
    M.write(out / "sources-after.json", {str(p): M.sha(p) for p in sources})
    print(route, name, data["artifact_sha256"], data["closure_units"], flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("route", choices=M.REVISIONS)
    args = ap.parse_args()
    compile_one(args.route, "generation_oracle", "test-integration-test")
    source = M.HERE / "observer.rs.txt"
    target = (
        M.OUT
        / "source"
        / args.route
        / "crates/feynkit-generator/examples/oracle_xs_qed_l2_capture.rs"
    )
    shutil.copy2(source, target)
    compile_one(args.route, "oracle_xs_qed_l2_capture", "example")
    for name, argv in (("rustc", ["rustc", "-vV"]), ("cargo", ["cargo", "--version"])):
        p, _ = M.run(argv, M.OUT / "build" / args.route / name)
        M.require(p.returncode == 0, "toolchain", name)


if __name__ == "__main__":
    main()
