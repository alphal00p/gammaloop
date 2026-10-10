#!/usr/bin/env python3
"""Time drawing fixtures; adapted from the supplied drawing-bench/run.py."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import shutil
import statistics
import subprocess
import sys
import tarfile
import time
from datetime import UTC, datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ARCHIVE_SHA256 = "324fe89dfb3950960b48455f4d02cee228e91883ccd55b1fb08b5f31af6adea6"
# Stock Typst Universe releases the drawing sources import.
PACKAGES = ("preview/cetz/0.5.2", "preview/mitex/0.2.6", "preview/oxifmt/1.0.0")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_revision() -> str | None:
    """Name the checked-out commit; source_hashes record working-copy edits.

    An exported tree outside any checkout has no commit to name.
    """
    if (REPO / ".jj").is_dir():
        command = [
            "jj",
            "--ignore-working-copy",
            "log",
            "--no-graph",
            "-r",
            "@",
            "-T",
            "commit_id",
        ]
    elif (REPO / ".git").exists():
        command = ["git", "rev-parse", "HEAD"]
    else:
        return None
    return subprocess.check_output(command, cwd=REPO, text=True).strip()


def inventory(root: Path) -> dict[str, str]:
    return {
        str(path.relative_to(root)): sha256(path)
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def extract(destination: Path) -> None:
    """Extract only regular immutable data, never archive links or absolute paths."""
    archive = HERE / "fixtures.tar"
    if sha256(archive) != ARCHIVE_SHA256:
        sys.exit("Fixture archive SHA256 mismatch.")
    with tarfile.open(archive) as handle:
        for member in handle:
            path = Path(member.name)
            if not member.isfile() or path.is_absolute() or ".." in path.parts:
                sys.exit(f"Unsafe fixture member: {member.name}")
            data = handle.extractfile(member)
            assert data is not None
            target = destination / path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data.read())


def build_root(build: Path, fixtures: Path) -> Path:
    """Assemble the current drawing sources; never apply old overlays."""
    root = build / "current"
    root.mkdir(parents=True)
    for name in ("linnest", "kurvst"):
        source = REPO / "crates" / name / "typst"
        target = root / "crates" / name / "typst"
        shutil.copytree(source / "src", target / "src")
        for path in sorted(source.iterdir()):
            if path.is_file() and (
                path.suffix == ".wasm"
                or path.name == "typst.toml"
                or path.name.startswith("LICENSE")
            ):
                shutil.copy2(path, target / path.name)
    templates = Path("assets/embedded/drawing/templates")
    shutil.copytree(REPO / templates, root / templates)
    shutil.copytree(fixtures / "cases", root / "cases")
    return root


def typst_env() -> dict[str, str]:
    """TYPST_FONT_PATHS could otherwise add non-embedded fonts, and the package
    cache could supply different copies of the checked packages."""
    env = os.environ.copy()
    env.pop("TYPST_FONT_PATHS", None)
    env.pop("TYPST_PACKAGE_CACHE_PATH", None)
    return env


def compile_case(
    typst: str,
    root: Path,
    packages: Path,
    case: str,
    output: Path,
    ppi: int | None = None,
    inputs: tuple[str, ...] = (),
) -> float:
    """Fresh process; wall time includes launch, compile, SVG export and exit."""
    command = [
        typst,
        "compile",
        "--root",
        str(root),
        "--package-path",
        str(packages),
        "--ignore-system-fonts",
        "--format",
        "png" if ppi else "svg",
    ]
    if ppi:
        command += ["--ppi", str(ppi)]
    for value in inputs:
        command += ["--input", value]
    command += [str(root / "cases" / f"{case}.typ"), str(output)]
    start = time.perf_counter()
    result = subprocess.run(
        command, capture_output=True, text=True, env=typst_env(), check=False
    )
    elapsed = (time.perf_counter() - start) * 1000
    if result.returncode:
        sys.exit(f"{case} failed:\n{result.stderr}")
    return elapsed


def label_report(typst: str, root: Path, packages: Path, case: str) -> list:
    """Each drawing's chosen labels and what they really meet: painted edges
    (`crossings`, summed length per edge), other labels, and the nodes,
    self-loop bends and momentum arrows whose boxes they overlap."""
    command = [
        typst,
        "eval",
        "query(<linnest-label-report>).map(record => record.value)",
        "--in",
        str(root / "cases" / f"{case}.typ"),
        "--root",
        str(root),
        "--package-path",
        str(packages),
        "--ignore-system-fonts",
        "--input",
        "linnest-label-report=true",
    ]
    result = subprocess.run(
        command, capture_output=True, text=True, env=typst_env(), check=False
    )
    if result.returncode:
        sys.exit(f"{case} label report failed:\n{result.stderr}")
    return json.loads(result.stdout)


def obstacle_name(obstacle: dict) -> str:
    """A node as `n<index>`, a self-loop bend as `e<edge>-loop`."""
    if "node" in obstacle:
        return f"n{obstacle['node']}"
    return f"e{obstacle['self-loop']}-loop"


def contacts(drawings: list) -> set[tuple[str, object, object]]:
    """(kind, label edge, what it meets) of every drawing in a case: crossed
    painted edges, obstacles and the edges of overlapped momentum arrows."""
    return {
        (kind, label["edge"], met)
        for drawing in drawings
        for label in drawing["labels"]
        for kind, met in [
            *(("crossing", crossing["edge"]) for crossing in label["crossings"]),
            *(("obstacle", obstacle_name(item)) for item in label["obstacles"]),
            *(("arrow", edge) for edge in label["arrows"]),
        ]
    }


def label_summary(drawings: list) -> dict[str, object]:
    labels = [label for drawing in drawings for label in drawing["labels"]]
    return {
        "labels": len(labels),
        "crossing_labels": sum(1 for label in labels if label["crossings"]),
        "label_pairs": sum(len(label["labels"]) for label in labels) // 2,
        "crossings": " ".join(
            f"e{label['edge']}>e{crossing['edge']}:{crossing['length']:.3f}"
            for label in labels
            for crossing in label["crossings"]
        ),
        "obstacles": " ".join(
            f"e{label['edge']}>{obstacle_name(item)}"
            for label in labels
            for item in label["obstacles"]
        ),
        "arrows": " ".join(
            f"e{label['edge']}>e{edge}" for label in labels for edge in label["arrows"]
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--typst", default=os.environ.get("TYPST", "typst"))
    parser.add_argument(
        "--package-path",
        type=Path,
        default=os.environ.get("TYPST_PACKAGE_CACHE_PATH"),
        help="Typst package directory with stock "
        + ", ".join(PACKAGES)
        + " (default: TYPST_PACKAGE_CACHE_PATH, as the dev shell sets it)",
    )
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--cases", nargs="+", default=[])
    parser.add_argument("--png", type=int, metavar="PPI")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--inspection-targets",
        action="store_true",
        help="render the interactive viewer's link targets for SVG/browser review "
        "(the corpus documents are static by default)",
    )
    parser.add_argument(
        "--label-report",
        action="store_true",
        help="after timing, record which painted edges, labels, nodes and "
        "momentum arrows every chosen label meets (label-report.csv/.json, untimed)",
    )
    parser.add_argument(
        "--label-baseline",
        type=Path,
        metavar="LABEL_REPORT_JSON",
        help="with --label-report, fail when a case gains a label meeting a "
        "painted edge, another label, a node or a momentum arrow relative to this "
        "earlier label-report.json",
    )
    parser.add_argument(
        "--out", type=Path, default=REPO / "target/drawing-bench/current"
    )
    args = parser.parse_args()
    if args.label_baseline and not args.label_report:
        parser.error("--label-baseline requires --label-report")
    if args.rounds < 3:
        parser.error("--rounds must be at least 3")
    if args.png is not None and args.png <= 0:
        parser.error("--png must be positive")
    if args.package_path is None:
        parser.error(
            "--package-path is required when TYPST_PACKAGE_CACHE_PATH is unset"
        )
    packages = args.package_path.resolve()
    for relative in PACKAGES:
        if not (packages / relative / "typst.toml").is_file():
            parser.error(f"missing package {relative} in {packages}")
    out = args.out.resolve()
    if not out.is_relative_to(REPO / "target/drawing-bench"):
        parser.error("--out must be inside target/drawing-bench")
    if out.exists():
        parser.error(
            "--out must be a new directory; existing results are never overwritten"
        )
    executable = shutil.which(args.typst)
    if executable is None:
        parser.error(f"Typst executable not found: {args.typst}")
    version = subprocess.check_output([executable, "--version"], text=True).strip()
    out.mkdir(parents=True)
    fixtures = out / "fixtures"
    extract(fixtures)
    cases = sorted(path.stem for path in (fixtures / "cases").glob("*.typ"))
    if len(cases) != 24:
        sys.exit("Expected exactly 24 archived cases.")
    cases = [
        case for case in cases if not args.cases or any(s in case for s in args.cases)
    ]
    if not cases:
        parser.error("No cases match --cases")
    root = build_root(out / "build", fixtures)
    metadata = {
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "source_revision": source_revision(),
        "source_identity": "CURRENT repository, including working-tree sources",
        "harness_sha256": sha256(HERE / "run.py"),
        "typst_executable": executable,
        "typst_sha256": sha256(Path(executable)),
        "typst_version": version,
        "host": platform.node(),
        "os": platform.platform(),
        "machine": platform.machine(),
        "python": platform.python_version(),
        "rounds": args.rounds,
        "samples": len(cases) * args.rounds,
        "warmup": "one untimed first-case compile",
        "fonts": "embedded only; ignore-system-fonts; TYPST_FONT_PATHS removed",
        "timing": "perf_counter wall ms, process launch through exit; SVG export included",
        "order": "round, sorted case",
        "fixture_archive_sha256": ARCHIVE_SHA256,
        "case_groups": {case: case.rsplit("-", 1)[0] for case in cases},
        "source_hashes": inventory(root),
        "package_path": str(packages),
        "package_hashes": {
            relative: inventory(packages / relative) for relative in PACKAGES
        },
        "policy": (
            "bounded arrow/label differences + explicit image review; "
            "links and hit coverage reviewed in --inspection-targets renders"
        ),
        "inspection_targets": args.inspection_targets,
        "label_report": args.label_report,
        "visual_review": "pending (not automated by this harness)",
        "dry_run": args.dry_run,
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    if args.dry_run:
        print(f"Prepared {len(cases)} cases: {out}")
        return
    renders = out / "current"
    renders.mkdir()
    inputs = ("linnest-inspection-targets=true",) if args.inspection_targets else ()
    compile_case(
        executable, root, packages, cases[0], renders / f"{cases[0]}.svg", None, inputs
    )
    times: dict[str, list[float]] = {}
    with (out / "timings.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["case", "round", "ms"])
        for round_index in range(args.rounds):
            for case in cases:
                ms = compile_case(
                    executable,
                    root,
                    packages,
                    case,
                    renders / f"{case}.svg",
                    None,
                    inputs,
                )
                times.setdefault(case, []).append(ms)
                writer.writerow([case, round_index + 1, f"{ms:.6f}"])
                handle.flush()
    medians = {case: statistics.median(times[case]) for case in cases}
    report = ["case,median_ms"]
    report += [f"{case},{median:.6f}" for case, median in medians.items()]
    report.append(f"TOTAL_MS,{sum(medians.values()):.6f}")
    (out / "medians.csv").write_text("\n".join(report) + "\n")
    print("\n".join(report))
    metadata["svg_hashes"] = {case: sha256(renders / f"{case}.svg") for case in cases}
    if args.png:
        for case in cases:
            compile_case(
                executable,
                root,
                packages,
                case,
                renders / f"{case}.png",
                args.png,
                inputs,
            )
        metadata["png_ppi"] = args.png
        metadata["png_hashes"] = {
            case: sha256(renders / f"{case}.png") for case in cases
        }
    if args.label_report:
        reports = {
            case: label_report(executable, root, packages, case) for case in cases
        }
        (out / "label-report.json").write_text(json.dumps(reports, indent=2) + "\n")
        summaries = {case: label_summary(reports[case]) for case in cases}
        with (out / "label-report.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            fields = [
                "labels",
                "crossing_labels",
                "label_pairs",
                "crossings",
                "obstacles",
                "arrows",
            ]
            writer.writerow(["case", *fields])
            for case, summary in summaries.items():
                writer.writerow([case, *(summary[field] for field in fields)])
        overlapping = {
            case: summary
            for case, summary in summaries.items()
            if contacts(reports[case]) or summary["label_pairs"]
        }
        print(f"Label overlaps in {len(overlapping)} of {len(cases)} cases:")
        for case, summary in overlapping.items():
            print(
                f"  {case}: crossings {summary['crossings'] or '-'}, "
                f"obstacles {summary['obstacles'] or '-'}, "
                f"arrows {summary['arrows'] or '-'}, "
                f"{summary['label_pairs']} label pairs"
            )
    regressions = []
    if args.label_baseline:
        baseline = json.loads(args.label_baseline.read_text())
        for case in cases:
            before = baseline.get(case, [])
            gained = contacts(reports[case]) - contacts(before)
            pairs = (
                summaries[case]["label_pairs"] - label_summary(before)["label_pairs"]
            )
            if gained or pairs > 0:
                regressions.append(
                    f"  {case}: new (kind, label edge, met) {sorted(gained, key=str)}, "
                    f"{max(pairs, 0)} new label pairs"
                )
        metadata["label_baseline"] = str(args.label_baseline.resolve())
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    if regressions:
        sys.exit("New label overlaps:\n" + "\n".join(regressions))


if __name__ == "__main__":
    main()
