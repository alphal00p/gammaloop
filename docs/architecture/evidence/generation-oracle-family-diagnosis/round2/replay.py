"""Fresh baseline, actual B1 regression, parser and aggregate predicate replays."""

import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

SPEC = importlib.util.spec_from_file_location(
    "certificate", Path(__file__).with_name("validate.py")
)
V = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V)
OUT = (
    Path(
        os.environ.get(
            "FAMILY_ROUND2_OUT",
            str(V.WORKSPACE / "target/oracle-numerator-family-round2"),
        )
    )
    / "replay-final"
)
REVIEW = V.WORKSPACE / "target/oracle-family-independent-review-r1"


def command(argv, root, env=None):
    root.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        [str(x) for x in argv], env=env, capture_output=True, text=True
    )
    (root / "process.log").write_text(proc.stdout + proc.stderr)
    row = {"argv": [str(x) for x in argv], "exit": proc.returncode}
    V.write(root / "command.json", row)
    return row


def algebra(pairs, root, orders=(False, True)):
    root.mkdir(parents=True, exist_ok=True)
    path = root / "pairs.json"
    V.write(path, pairs)
    commands = []
    for reverse in orders:
        try:
            commands += V.run_exact(path, root, reverse)
        except V.Reject as error:
            return {"disposition": "reject", "stage": error.stage, "detail": str(error)}
    return {"disposition": "accept", "commands": commands}


def main():
    baseline = []
    for case in V.CASES:
        for mode in V.MODES:
            root = OUT / "baseline" / case / mode
            argv = [
                sys.executable,
                "-I",
                "-B",
                Path(__file__).with_name("validate.py"),
                "--case",
                case,
                "--mode",
                mode,
                "--out",
                root,
            ]
            row = command(argv, root)
            row.update(case=case, mode=mode, verdict=V.load(root / "verdict.json"))
            baseline.append(row)
            V.write(OUT / "baseline.json", baseline)
            if row["exit"] != 0:
                raise RuntimeError("baseline failed " + case + mode)
            print("baseline", case, mode, "accepted", flush=True)

    # The actual three reviewer mutants are consumed unchanged by both old and
    # new validators. Output directories are round-two supplements only.
    regressions = []
    for mutation, case in (
        ("f-attribute-dropped-current", "ym_g6_L0"),
        ("g-symmetric-dropped-current", "qed_eemm_L1"),
        ("gamma-linear-dropped-current", "qed_eemm_L1"),
    ):
        source = REVIEW / "work/mutations" / mutation
        root = OUT / "B1" / mutation
        env = dict(
            os.environ,
            REVIEW_HIST=str(source / "hist"),
            REVIEW_CUR=str(source / "cur"),
            REVIEW_OUT=str(root / "old-pairs"),
            PYTHONDONTWRITEBYTECODE="1",
        )
        argv = [
            sys.executable,
            "-I",
            "-B",
            REVIEW / "probes/reviewed_validator_redirectable.py",
            case,
        ]
        old = command(argv, root / "old-validator", env)
        old_exact = []
        for mode in V.MODES:
            for kind in ("local", "graph"):
                path = (
                    root
                    / "old-pairs"
                    / case
                    / f"classes-{mode}-t1"
                    / f"{kind}-pairs.json"
                )
                # Deliberately bypass new symbol gate to reproduce old acceptance.
                old_exact.append(
                    algebra(V.load(path), root / "old-checkers" / mode / kind, (False,))
                )
        new = []
        for mode in V.MODES:
            dest = root / "new" / mode
            argv = [
                sys.executable,
                "-I",
                "-B",
                Path(__file__).with_name("validate.py"),
                "--case",
                case,
                "--mode",
                mode,
                "--historical",
                source / "hist",
                "--current",
                source / "cur",
                "--out",
                dest,
            ]
            row = command(argv, dest)
            row["verdict"] = V.load(dest / "verdict.json")
            new.append(row)
        regressions.append(
            {
                "mutation": mutation,
                "old_validator": old,
                "old_exact": old_exact,
                "new": new,
            }
        )
        V.write(OUT / "B1.json", regressions)

    parser_rows = []
    probes = V.load(REVIEW / "work/parser-probes/manifest.json")
    for probe in probes:
        path = Path(probe["file"])
        if not path.is_absolute():
            path = REVIEW / path
        pairs = V.load(path)
        root = OUT / "parser" / probe["name"]
        try:
            checked = [
                V.pair(p["historical"], p["current"], probe["name"]) for p in pairs
            ]
            observed = algebra(checked, root)
        except V.Reject as error:
            observed = {
                "disposition": "reject",
                "stage": error.stage,
                "detail": str(error),
            }
        # Also invoke each kernel on the original strings so string-level
        # corruptions have real parser-path evidence, distinct from our gate.
        raw = algebra(pairs, root / "raw-parser")
        row = {
            **probe,
            "observed": observed,
            "raw_parser": raw,
            "matched": observed["disposition"] == probe["expect"],
        }
        parser_rows.append(row)
        V.write(OUT / "parser.json", parser_rows)

    aggregate_controls = []
    for case in V.CASES:
        for mode in V.MODES:
            rows = V.load(OUT / "baseline" / case / mode / "pairs.json")
            chosen = next(
                p for p in rows if p["label"].startswith("aggregate-product/current/")
            )
            corrupt = {**chosen, "current": "2*(" + chosen["current"] + ")"}
            result = algebra([corrupt], OUT / "aggregate-predicate" / case / mode)
            aggregate_controls.append(
                {
                    "case": case,
                    "mode": mode,
                    "label": chosen["label"],
                    "observed": result,
                    "matched": result["disposition"] == "reject",
                }
            )
    V.write(OUT / "aggregate-predicate.json", aggregate_controls)
    if not all(r["matched"] for r in parser_rows + aggregate_controls):
        raise RuntimeError("unmatched replay intent")


if __name__ == "__main__":
    main()
