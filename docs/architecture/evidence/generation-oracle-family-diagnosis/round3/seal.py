"""Detached freeze receipt; no source edits, ref moves, or external actions."""

import argparse
import os
from pathlib import Path

import common as C


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("revision")
    args = ap.parse_args()
    out = C.OUT / "freeze"
    expected = {
        "generation-oracle-family-round1": "e814b6b4c3463d9aa1c018939365b624619af24b",
        "generation-oracle-family-round2": "78c0754c8fb732878123c02999ced3fc506c5a1d",
        "generation-oracle-family-historical-diagnostic": "5e61ef89f4df4e818884e40b8d0c46108eab75a8",
        "generation-oracle-family-round3": args.revision,
    }
    refs = {}
    for name, revision in expected.items():
        proc, command = C.run(
            [
                "jj",
                "log",
                "--ignore-working-copy",
                "--no-graph",
                "-r",
                name,
                "-T",
                "commit_id",
            ],
            out / name,
        )
        C.require(proc.returncode == 0 and proc.stdout == revision, "freeze", name)
        refs[name] = command
    proc, _ = C.run(
        [
            "jj",
            "diff",
            "--ignore-working-copy",
            "--from",
            expected["generation-oracle-family-round2"],
            "--to",
            args.revision,
            "--no-pager",
        ],
        out / "diff",
    )
    C.require(proc.returncode == 0, "freeze", "diff")
    proc, _ = C.run(
        [
            "jj",
            "--ignore-working-copy",
            "-R",
            str(C.WORKSPACE.parent / "generation-oracle-family-historical"),
            "status",
            "--no-pager",
        ],
        out / "historical-status",
    )
    C.require(
        proc.returncode == 0
        and "has no changes" in proc.stdout
        and "6d476227" in proc.stdout,
        "freeze",
        "historical holder",
    )
    processes = []
    for path in Path("/proc").iterdir():
        if not path.name.isdecimal() or int(path.name) == os.getpid():
            continue
        try:
            argv = (path / "cmdline").read_bytes().decode().split("\0")
        except (OSError, UnicodeDecodeError):
            continue
        if not argv:
            continue
        if (
            argv[0] == str(C.B.ELFS["current"])
            or argv[0] == str(C.B.ELFS["historical"])
            or any(
                a == str(C.HERE / n)
                for a in argv
                for n in ("suite.py", "validate.py", "report_checks.py")
            )
        ):
            processes.append({"pid": int(path.name), "argv": argv})
    C.require(not processes, "freeze", "owned proof processes remain")
    commands = C.files(C.OUT / "commands")
    open_commands = [
        str(p.parent)
        for p in (C.OUT / "commands").glob("*/command.json")
        if "ended_utc" not in C.load(p)
    ]
    commands = {
        p: h
        for p, h in commands.items()
        if not any(p.startswith(root + "/") for root in open_commands)
    }
    C.write(
        out / "completed-top-level-command-manifest.json",
        {
            "files": commands,
            "excluded_open_self_commands": open_commands,
            "reason": "A command record cannot contain its own final hash; this detached seal is not a mathematical execution.",
        },
    )
    C.write(
        out / "attempt-manifest.json",
        {
            "files": C.files(C.OUT / "attempts"),
            "credit": "none; exploratory and failed setup history",
        },
    )
    inputs = [
        C.OUT / "inventory/closure-final.json",
        C.OUT / "inventory/frozen-before.json",
        C.OUT / "inventory/frozen-after.json",
        C.OUT / "docs/render-validation.json",
        C.OUT / "docs/registry.json",
        C.HERE / "receipt.json",
        C.HERE / "dispositions.json",
        out / "completed-top-level-command-manifest.json",
        out / "attempt-manifest.json",
    ]
    C.write(
        C.OUT / "freeze.json",
        {
            "candidate_bookmark": "generation-oracle-family-round3",
            "revision": args.revision,
            "independent_review": "pending; root commissions; no review launched by this lane",
            "refs": refs,
            "owned_processes": processes,
            "source_files": C.files(C.HERE),
            "artifacts": {str(p): C.sha(p) for p in inputs},
            "source_inventory_limit": "Inventories count files, not a general runtime/physics invariance theorem.",
            "scope": "four cases independently, None/Auto, Grouping None, pool1; no migration or expansion",
        },
    )
    print(args.revision, C.sha(C.OUT / "freeze.json"))


if __name__ == "__main__":
    main()
