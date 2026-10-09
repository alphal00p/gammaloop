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
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ARCHIVE_SHA256 = "8f5cf9c5ba70d466f69b800882552aa0b41171efd96cd8fe73b05fb585a99e1f"
CETZ_VERSION = "0.5.2"
CETZ_ARCHIVE_SHA256 = "77cf8490114ae04c6e665a11efa691d284a0cadb9719771b5708c1197292f23f"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inventory(root: Path) -> dict[str, str]:
    return {
        str(path.relative_to(root)): sha256(path)
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def extract(archive: Path, expected_hash: str, destination: Path) -> None:
    """Extract immutable regular files/directories, never links or escaping paths."""
    if sha256(archive) != expected_hash:
        sys.exit(f"Archive SHA256 mismatch: {archive}")
    with tarfile.open(archive) as handle:
        for member in handle:
            path = Path(member.name)
            if (
                not (member.isfile() or member.isdir())
                or path.is_absolute()
                or ".." in path.parts
            ):
                sys.exit(f"Unsafe archive member: {member.name}")
            target = destination / path
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            data = handle.extractfile(member)
            assert data is not None
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data.read())


def build_root(variant: str, build: Path, fixtures: Path) -> Path:
    """Use identical current sources for both variants; never apply old overlays."""
    root = build / variant
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


def compile_case(
    typst: str,
    root: Path,
    packages: Path,
    case: str,
    output: Path,
    ppi: int | None = None,
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
    command += [str(root / "cases" / f"{case}.typ"), str(output)]
    # TYPST_FONT_PATHS could otherwise add non-embedded fonts.
    env = os.environ.copy()
    env.pop("TYPST_FONT_PATHS", None)
    start = time.perf_counter()
    result = subprocess.run(command, capture_output=True, text=True, env=env)
    elapsed = (time.perf_counter() - start) * 1000
    if result.returncode:
        sys.exit(f"{root.name}/{case} failed:\n{result.stderr}")
    return elapsed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--typst", default=os.environ.get("TYPST", "typst"))
    parser.add_argument(
        "--variants",
        nargs="+",
        choices=["current", "current-stock"],
        default=["current"],
    )
    parser.add_argument("--stock-package-path", type=Path)
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--cases", nargs="+", default=[])
    parser.add_argument("--png", type=int, metavar="PPI")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--out", type=Path, default=REPO / "target/drawing-bench/current"
    )
    args = parser.parse_args()
    if args.rounds < 3:
        parser.error("--rounds must be at least 3")
    if args.png is not None and args.png <= 0:
        parser.error("--png must be positive")
    if len(set(args.variants)) != len(args.variants):
        parser.error("variants must be unique")
    if "current-stock" in args.variants and args.stock_package_path is None:
        parser.error(
            "current-stock requires --stock-package-path (no historical overlay)"
        )
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
    bundled = REPO / "crates/linnet-py/vendor/typst-packages"
    packages = {"current": out / "packages"}
    shutil.copytree(bundled / "preview", packages["current"] / "preview")
    extract(
        bundled / "archives" / f"cetz-{CETZ_VERSION}.tar.gz",
        CETZ_ARCHIVE_SHA256,
        packages["current"] / "preview" / "cetz" / CETZ_VERSION,
    )
    if args.stock_package_path:
        packages["current-stock"] = args.stock_package_path.resolve()
    for variant in args.variants:
        for relative in (
            f"preview/cetz/{CETZ_VERSION}",
            "preview/mitex/0.2.6",
            "preview/oxifmt/1.0.0",
        ):
            if not (packages[variant] / relative / "typst.toml").is_file():
                parser.error(f"{variant}: missing package {relative}")
    fixtures = out / "fixtures"
    extract(HERE / "fixtures.tar", ARCHIVE_SHA256, fixtures)
    cases = sorted(path.stem for path in (fixtures / "cases").glob("*.typ"))
    if len(cases) != 24:
        sys.exit("Expected exactly 24 archived cases.")
    cases = [
        case for case in cases if not args.cases or any(s in case for s in args.cases)
    ]
    if not cases:
        parser.error("No cases match --cases")
    roots = {v: build_root(v, out / "build", fixtures) for v in args.variants}
    metadata = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "source_revision": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
        ).strip(),
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
        "samples": len(cases) * len(args.variants) * args.rounds,
        "warmup": "one untimed first-case compile per variant",
        "fonts": "embedded only; ignore-system-fonts; TYPST_FONT_PATHS removed",
        "timing": "perf_counter wall ms, process launch through exit; SVG export included",
        "order": "round, sorted case, variant in requested order",
        "variants": args.variants,
        "fixture_archive_sha256": ARCHIVE_SHA256,
        "cetz_version": CETZ_VERSION,
        "cetz_archive_sha256": CETZ_ARCHIVE_SHA256,
        "case_groups": {case: case.rsplit("-", 1)[0] for case in cases},
        "source_hashes": {v: inventory(roots[v]) for v in args.variants},
        "package_paths": {v: str(packages[v]) for v in args.variants},
        "package_hashes": {v: inventory(packages[v]) for v in args.variants},
        "policy": (
            "bounded arrow/label differences + explicit image review; "
            "links and hit coverage required"
        ),
        "visual_review": "pending (not automated by this harness)",
        "dry_run": args.dry_run,
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    if args.dry_run:
        print(f"Prepared {len(cases)} cases: {out}")
        return
    for variant in args.variants:
        (out / variant).mkdir()
        compile_case(
            executable,
            roots[variant],
            packages[variant],
            cases[0],
            out / variant / f"{cases[0]}.svg",
        )
    times: dict[tuple[str, str], list[float]] = {}
    with (out / "timings.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["variant", "case", "round", "ms"])
        for round_index in range(args.rounds):
            for case in cases:
                for variant in args.variants:
                    ms = compile_case(
                        executable,
                        roots[variant],
                        packages[variant],
                        case,
                        out / variant / f"{case}.svg",
                    )
                    times.setdefault((variant, case), []).append(ms)
                    writer.writerow([variant, case, round_index + 1, f"{ms:.6f}"])
                    handle.flush()
    report = ["case," + ",".join(f"{v}_ms" for v in args.variants) + ",last/first"]
    totals = dict.fromkeys(args.variants, 0.0)
    for case in cases:
        medians = [statistics.median(times[v, case]) for v in args.variants]
        for variant, median in zip(args.variants, medians):
            totals[variant] += median
        report.append(
            case
            + ","
            + ",".join(f"{m:.6f}" for m in medians)
            + f",{medians[-1] / medians[0]:.6f}"
        )
    report.append(
        "TOTAL_MS,"
        + ",".join(f"{totals[v]:.6f}" for v in args.variants)
        + f",{totals[args.variants[-1]] / totals[args.variants[0]]:.6f}"
    )
    (out / "medians.csv").write_text("\n".join(report) + "\n")
    print("\n".join(report))
    metadata["svg_hashes"] = {
        v: {case: sha256(out / v / f"{case}.svg") for case in cases}
        for v in args.variants
    }
    if args.png:
        for variant in args.variants:
            for case in cases:
                compile_case(
                    executable,
                    roots[variant],
                    packages[variant],
                    case,
                    out / variant / f"{case}.png",
                    args.png,
                )
        metadata["png_ppi"] = args.png
        metadata["png_hashes"] = {
            v: {case: sha256(out / v / f"{case}.png") for case in cases}
            for v in args.variants
        }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
