"""Measure the example Marimo notebooks in alternating, fresh processes."""

import argparse
import ast
import hashlib
import html
import json
import math
import os
import subprocess
import sys
import tempfile
from itertools import chain
from pathlib import Path

from gammabench.gate_policy import PairedGate

COLLECTOR = Path(__file__).resolve().parents[2] / "scripts/measure-notebook-cells.py"
NOTEBOOKS = Path("examples/notebooks")
SIDES = ("baseline", "candidate")
PAIRS = 5
RELATIVE_LIMIT = 0.20
ABSOLUTE_LIMIT = 0.050  # seconds per cell
TIMEOUT_SECONDS = 600
VERDICTS = {0: "passed", 1: "regression", 2: "inconclusive"}


def notebooks(root):
    """Recognize Marimo constructors using the imports in each Python file."""
    paths = set()
    for path in chain(
        root.glob("examples/**/*.py"), root.glob("crates/*/examples/**/*.py")
    ):
        source = path.read_bytes()
        if b"marimo" not in source:
            continue
        tree = ast.parse(source, filename=str(path))
        modules, constructors = set(), set()
        for stmt in tree.body:
            if isinstance(stmt, ast.Import):
                modules.update(
                    alias.asname or alias.name
                    for alias in stmt.names
                    if alias.name == "marimo"
                )
            elif isinstance(stmt, ast.ImportFrom) and stmt.module == "marimo":
                constructors.update(
                    alias.asname or alias.name
                    for alias in stmt.names
                    if alias.name == "App"
                )
        for stmt in tree.body:
            if not isinstance(stmt, (ast.Assign, ast.AnnAssign)) or not isinstance(
                stmt.value, ast.Call
            ):
                continue
            func = stmt.value.func
            if (isinstance(func, ast.Name) and func.id in constructors) or (
                isinstance(func, ast.Attribute)
                and func.attr == "App"
                and isinstance(func.value, ast.Name)
                and func.value.id in modules
            ):
                paths.add(path.relative_to(root))
                break
    return paths


def collect(root, python, notebook):
    env = os.environ.copy()
    for name in ("PYTHONPATH", "PYTHONHOME", "LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH"):
        env.pop(name, None)
    env.update(PYTHONNOUSERSITE="1", PYTHONDONTWRITEBYTECODE="1")
    with tempfile.TemporaryDirectory(prefix="marimo-cell-") as temporary:
        output = Path(temporary) / "measurement.json"
        subprocess.run(
            [
                python,
                str(COLLECTOR),
                "--notebook",
                str(root / notebook),
                "--output",
                str(output),
            ],
            cwd=root,
            env=env,
            check=True,
            timeout=TIMEOUT_SECONDS,
        )
        return json.loads(output.read_text())


def validate(sample):
    cells = sample["cells"]
    if not isinstance(cells, dict) or not cells:
        raise ValueError("collector returned no timed cells")
    for name, elapsed in cells.items():
        if (
            not isinstance(name, str)
            or not name
            or type(elapsed) not in (float, int)
            or not math.isfinite(elapsed)
            or elapsed <= 0
        ):
            raise ValueError(f"invalid cell timing for {name!r}")
    if "labels" in sample and (
        not isinstance(sample["labels"], dict)
        or cells.keys() != sample["labels"].keys()
        or any(not isinstance(label, str) for label in sample["labels"].values())
    ):
        raise ValueError("invalid cell labels")
    return cells


def notebook_gate(path, sources, entry, *, bootstrap=False):
    entry.update(phases={}, cells={})
    inventories = {}
    scenarios = {}
    labels = {}

    def run(side):
        root, python = sources[side]
        sample = collect(root, python, path)
        cells = validate(sample)
        if set(cells) != inventories.setdefault(side, set(cells)):
            raise ValueError(f"{side} cell coverage changed between repetitions")
        name = sample["scenario"]
        settings = sample["settings"]
        if not isinstance(name, str) or not name or not isinstance(settings, dict):
            raise ValueError(f"{side} returned invalid scenario metadata")
        scenario = {"name": name, "settings": settings}
        if scenarios.setdefault(side, scenario) != scenario:
            raise ValueError(f"{side} notebook scenario changed between repetitions")
        labels[side] = sample.get("labels", {})
        entry["scenarios"] = scenarios
        return cells

    if bootstrap or len(sources) == 1:
        side = "candidate" if "candidate" in sources else "baseline"
        status = (
            "baseline-unavailable"
            if bootstrap
            else {"baseline": "removed", "candidate": "added"}[side]
        )
        if not bootstrap or side == "candidate":
            cells = run(side)
            entry["phases"]["observation"] = {side: cells}
            entry["cells"] = {
                cell: {"status": status, "label": labels[side].get(cell, cell)}
                for cell in cells
            }
        entry["status"] = status
        return

    policy = PairedGate(
        None, None, {"cell_seconds": (RELATIVE_LIMIT, ABSOLUTE_LIMIT)}, ()
    )
    for phase in ("discovery", "confirmation"):
        pairs = []
        entry["phases"][phase] = pairs
        for i in range(PAIRS):
            order = SIDES if (i + (phase == "confirmation")) % 2 == 0 else SIDES[::-1]
            pair = {"order": list(order)}
            pairs.append(pair)
            for side in order:
                pair[side] = run(side)
            if scenarios["baseline"] != scenarios["candidate"]:
                raise ValueError(
                    "baseline and candidate executed different scenarios or settings"
                )

        for cell in sorted(inventories["baseline"] | inventories["candidate"]):
            record = entry["cells"].setdefault(
                cell,
                {
                    "status": "unmeasured",
                    "label": labels["candidate"].get(
                        cell, labels["baseline"].get(cell, cell)
                    ),
                },
            )
            if cell not in inventories["baseline"]:
                record["status"] = "added"
            elif cell not in inventories["candidate"]:
                record["status"] = "removed"
            else:
                record[phase] = policy.summarize(
                    [(p["baseline"][cell], p["candidate"][cell]) for p in pairs],
                    "cell_seconds",
                )
                first = record["discovery"]["status"]
                second = record[phase]["status"]
                if first == second:
                    record["status"] = VERDICTS[first]
                elif {first, second} == {0, 2}:
                    # One clear batch can resolve a noisy batch only if the
                    # combined samples also rule out a material slowdown.
                    record["combined"] = policy.summarize(
                        [
                            (pair["baseline"][cell], pair["candidate"][cell])
                            for batch in entry["phases"].values()
                            for pair in batch
                        ],
                        "cell_seconds",
                    )
                    record["status"] = (
                        "passed"
                        if record["combined"]["status"] == 0
                        else "inconclusive"
                    )
                else:
                    record["status"] = "inconclusive"
        if phase == "discovery" and not any(
            record.get("discovery", {}).get("status")
            for record in entry["cells"].values()
        ):
            break
    entry["status"] = outcome(entry["cells"].values())


def outcome(records):
    statuses = {record["status"] for record in records}
    for status in ("error", "regression", "inconclusive", "passed"):
        if status in statuses:
            return status
    return "uncompared"


def checkout_revision(root):
    if (root / ".jj").is_dir():
        command = [
            "jj",
            "log",
            "--ignore-working-copy",
            "-R",
            str(root),
            "--no-graph",
            "-r",
            "@",
            "-T",
            "commit_id",
        ]
    elif (root / ".git").exists():
        command = ["git", "-C", str(root), "rev-parse", "--verify", "HEAD"]
    else:
        # An exported tree has no revision metadata of its own.
        return None
    try:
        return subprocess.check_output(
            command,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def write_report(path, report):
    # A job timeout should retain all completed notebook comparisons.
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as file:
        json.dump(report, file, indent=2, allow_nan=False)
        file.write("\n")
    Path(file.name).replace(path)


def report_summary(report):
    rows = [f"### Marimo cell timings — {report['status']}", ""]
    for path, entry in report["notebooks"].items():
        counts = {}
        for record in entry.get("cells", {}).values():
            counts[record["status"]] = counts.get(record["status"], 0) + 1
        detail = ", ".join(
            f"{count} {status}" for status, count in sorted(counts.items())
        )
        rows.append(
            f"- <code>{html.escape(path)}</code> — {entry['status']}"
            + (f" ({detail})" if detail else "")
        )
        if "error" in entry:
            rows.append("  - " + html.escape(entry["error"]))
        for record in entry.get("cells", {}).values():
            if record["status"] not in ("regression", "inconclusive"):
                continue
            phase = record.get("confirmation", record["discovery"])
            rows.append(
                f"  - <code>{html.escape(record['label'])}</code> — {record['status']}, "
                f"{phase['baseline_median']:.4g} → {phase['candidate_median']:.4g} s "
                f"({phase['median_relative_delta']:+.1%})"
            )
    if "error" in report:
        rows.append(html.escape(report["error"]))
    elif "baseline_unavailable_reason" in report:
        rows.append(html.escape(report["baseline_unavailable_reason"]))
    return "\n".join(rows) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for side in SIDES:
        parser.add_argument(f"--{side}-root", type=Path, required=True)
        parser.add_argument(f"--{side}-python", type=Path, required=side == "candidate")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path)
    parser.add_argument(
        "--notebook",
        type=Path,
        action="append",
        help="select a repository-relative notebook",
    )
    args = parser.parse_args()
    sources = {}
    for side in SIDES:
        python = getattr(args, side + "_python")
        sources[side] = (
            getattr(args, side + "_root").resolve(),
            python.absolute() if python else None,
        )
    report = {
        "schema_version": 1,
        "status": "incomplete",
        "clock_unit": "seconds",
        "relative_limit": RELATIVE_LIMIT,
        "absolute_limit_seconds": ABSOLUTE_LIMIT,
        "phase_pairs": PAIRS,
        "interpreters": {
            s: str(python.resolve()) if python else None
            for s, (_, python) in sources.items()
        },
        "revisions": {s: checkout_revision(root) for s, (root, _) in sources.items()},
        "notebooks": {},
    }
    try:
        for root, python in sources.values():
            if not root.is_dir():
                raise ValueError(f"missing checkout: {root}")
            if python is not None and not (
                python.is_file() and os.access(python, os.X_OK)
            ):
                raise ValueError(f"missing or non-executable interpreter: {python}")
        paths = {s: notebooks(root) for s, (root, _) in sources.items()}
        if not args.notebook and not paths["candidate"]:
            raise ValueError("Full suite found no candidate Marimo notebooks")
        bootstrap = args.baseline_python is None
        if bootstrap:
            host_manifest = NOTEBOOKS / "symbolica-host/Cargo.toml"
            if (sources["baseline"][0] / host_manifest).exists() or not (
                sources["candidate"][0] / host_manifest
            ).is_file():
                raise ValueError("--baseline-python required for this comparison")
            report["baseline_unavailable_reason"] = (
                "The baseline predates the native notebook host. Only candidate execution was measured."
            )
            report["baseline_inventory"] = sorted(
                p.as_posix() for p in paths["baseline"]
            )
        found = paths["baseline"] | paths["candidate"]
        selected = set(args.notebook) if args.notebook else found
        if not selected or selected - found:
            raise ValueError(
                f"No notebook found for selection: {sorted(selected - found)}"
            )
        if bootstrap and not selected & paths["candidate"]:
            raise ValueError("bootstrap requires at least one candidate notebook")
        for path in sorted(selected):
            sides = {s: pair for s, pair in sources.items() if path in paths[s]}
            entry = {
                "status": "incomplete",
                "source_sha256": {
                    s: hashlib.sha256((root / path).read_bytes()).hexdigest()
                    for s, (root, _) in sides.items()
                },
            }
            report["notebooks"][path.as_posix()] = entry
            write_report(args.output, report)
            try:
                notebook_gate(path, sides, entry, bootstrap=bootstrap)
            except (
                ValueError,
                OSError,
                KeyError,
                TypeError,
                subprocess.SubprocessError,
            ) as error:
                entry["status"] = "error"
                entry["error"] = str(error)
                print(f"{path}: {error}", file=sys.stderr)
        report["status"] = outcome(report["notebooks"].values())
        if bootstrap and report["status"] == "uncompared":
            report["status"] = "baseline-unavailable"
    except (OSError, ValueError, SyntaxError) as error:
        report["status"] = "error"
        report["error"] = str(error)
    write_report(args.output, report)
    summary = report_summary(report)
    print(summary, end="")
    if args.summary:
        with args.summary.open("a") as handle:
            handle.write(summary)
    return {"error": 2, "regression": 1, "inconclusive": 2}.get(report["status"], 0)


if __name__ == "__main__":
    sys.exit(main())
