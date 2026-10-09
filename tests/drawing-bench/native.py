"""Focused installed FeynKit Python native render acceptance (no Typst corpus)."""

import argparse
import hashlib
import json
import math
import platform
import re
import statistics
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "crates/feynkit-py/tests/fixtures/scalars_2p_3p.json"
DOT = """digraph native_cube_minus_edge {
  ext [style=invis];
  edge [particle="scalar_0"];
  ext -> a; b -> ext;
  b -> c; c -> d; d -> a;
  e -> f; f -> g; g -> h; h -> e;
  a -> e; b -> f; c -> g; d -> h;
}
"""
# FeynKit's public native API takes a JSON dictionary, not ln.RenderConfig.
# This is the equivalent of LayoutOptions(steps=100), with no label pages.
CONFIG = {
    "layouts": {"steps": 100},
    "template_options": {
        "momentum-arrows": True,
        "show-momentum": False,
        "show-particle": False,
    },
}
NS = "{http://www.w3.org/2000/svg}"
NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def command(*args):
    return subprocess.check_output(args, cwd=ROOT, text=True).strip()


def inventory():
    """Hash actual inputs, including dirty sources, without copying private code."""
    roots = [
        "crates/feynkit-py",
        "crates/feynkit-graph",
        "crates/feynkit-model",
        "crates/linnet-py",
        "crates/linnest/src",
        "crates/kurvst/src",
        "crates/typst-renderer",
        "examples/notebooks/symbolica-host",
    ]
    paths = [ROOT / "Cargo.lock", ROOT / "Cargo.toml", Path(__file__).resolve()]
    for root in roots:
        paths.extend(
            p
            for p in (ROOT / root).rglob("*")
            if p.is_file()
            and p.suffix in {".rs", ".py", ".json", ".toml", ".lock"}
            and not {"target", "__pycache__"}.intersection(p.parts)
        )
    return {
        str(p.relative_to(ROOT)): digest(p.read_bytes()) for p in sorted(set(paths))
    }


def inspect_svg(svg):
    """Check painted open-V arrow geometry and actual interactive hit targets."""
    root = ET.fromstring(svg)
    assert root.tag == NS + "svg", "not an SVG document"
    viewbox = [float(v) for v in root.attrib["viewBox"].split()]
    assert len(viewbox) == 4 and all(math.isfinite(v) for v in viewbox)
    assert viewbox[2] > 0 and viewbox[3] > 0
    heads = []
    for path in root.iter(NS + "path"):
        data = path.get("d", "")
        # Native mark paths carry explicit join/miter paint attributes. Scalar
        # edge shafts do not. The default momentum mark is an open stroked V,
        # not a filled triangle: closing it would test the wrong geometry.
        if "stroke-miterlimit" not in path.attrib:
            continue
        assert path.get("stroke") not in {None, "none", "transparent"}
        assert float(path.attrib["stroke-width"]) > 0
        coords = [float(v) for v in re.findall(NUMBER, data)]
        assert coords and len(coords) % 2 == 0
        assert all(math.isfinite(v) for v in coords)
        points = list(zip(coords[::2], coords[1::2]))
        assert len(points) >= 3
        # Noncollinear vertices/control points certify an actual arrow outline,
        # not merely a path tag, a repeated point or a straight shaft.
        area = (
            abs(
                sum(
                    x * ny - nx * y
                    for (x, y), (nx, ny) in zip(points, points[1:] + points[:1])
                )
            )
            / 2
        )
        assert area > 1e-8, "degenerate momentum arrow"
        heads.append(data)
    targets = []
    for link in root.iter(NS + "a"):
        kind = link.get("data-linnet-kind")
        if kind is None:
            continue
        assert kind in {"node", "edge", "halfedge"}
        assert link.attrib["data-linnet-id"].isdigit()
        assert isinstance(json.loads(link.attrib["data-linnet-detail"]), dict)
        assert link.get("role") == "button" and link.get("aria-label")
        boxes = list(link.iter(NS + "rect"))
        assert boxes, "interactive target has no hit geometry"
        assert all(
            math.isfinite(float(box.attrib[key])) and float(box.attrib[key]) > 0
            for box in boxes
            for key in ("width", "height")
        )
        targets.append((kind, link.attrib["data-linnet-id"]))
    assert targets, "missing interactive links"
    for element in root.iter():
        href = element.get("href") or element.get("{http://www.w3.org/1999/xlink}href")
        if href and href.startswith("#"):
            assert any(n.get("id") == href[1:] for n in root.iter()), href
    return {"painted_heads": heads, "interactive_targets": sorted(set(targets))}


def run(out, metadata):
    import linnet as ln
    import symbolica.community.hepkit as fk

    # Record the extensions actually loaded, not just Python wrapper identities.
    modules = {}
    for name, module in sorted(sys.modules.copy().items()):
        if name == "linnet" or name.startswith(("linnet.", "symbolica")):
            path = getattr(module, "__file__", None)
            if path and Path(path).is_file():
                modules[name] = {
                    "path": str(Path(path).resolve()),
                    "sha256": digest(Path(path).read_bytes()),
                }
    metadata["modules"] = modules
    metadata["linnet_module"] = ln.__name__
    metadata["feynkit_module"] = fk.__name__
    model = fk.Model(str(FIXTURE))
    diagram = fk.FeynmanDiagram.from_dot(model, DOT)
    diagram.validate()
    assert diagram.loop_count == 4
    assert len(diagram.internal_edges) == 11
    assert len(diagram.vertices) == 8
    assert len(diagram.external_edges) == 2
    metadata["diagram"] = {
        "loops": 4,
        "internal_edges": 11,
        "vertices": 8,
        "external_edges": 2,
    }
    without_arrows = {
        **CONFIG,
        "template_options": {
            **CONFIG["template_options"],
            "momentum-arrows": False,
        },
    }
    plain = diagram.render(config=without_arrows)
    arrows = diagram.render(config=CONFIG)
    (out / "arrows-off.svg").write_text(plain)
    (out / "arrows-on.svg").write_text(arrows)
    plain_check, arrow_check = inspect_svg(plain), inspect_svg(arrows)
    # No scalar flow arrows: the on/off differential isolates momentum heads.
    assert not plain_check["painted_heads"]
    assert len(arrow_check["painted_heads"]) == 13
    assert plain_check["interactive_targets"] == arrow_check["interactive_targets"]
    metadata["geometry"] = arrow_check
    for _ in range(10):
        diagram.render(config=CONFIG)
    samples = []
    hashes = []
    for index in range(21):
        start = time.perf_counter_ns()
        svg = diagram.render(config=CONFIG)
        elapsed = time.perf_counter_ns() - start
        samples.append(elapsed)
        # Parsing, correctness checks, hashing and disk writes are untimed.
        assert inspect_svg(svg) == arrow_check
        data = svg.encode()
        hashes.append(digest(data))
        (out / f"sample-{index:02d}.svg").write_bytes(data)
    median_ms = statistics.median(samples) / 1_000_000
    metadata.update(
        status="passed" if median_ms <= 100 else "target-exceeded",
        samples_ns=samples,
        median_ms=median_ms,
        svg_sha256=hashes,
        target_ms=100,
        target_met=median_ms <= 100,
        nonregression="not established: no genuine old same-host baseline",
        visual_review="pending: geometry/link checks are not image approval",
    )
    print(f"native Python median: {median_ms:.3f} ms; target <= 100 ms")
    return 0 if median_ms <= 100 else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out", type=Path, default=ROOT / "target/drawing-bench/native-acceptance"
    )
    args = parser.parse_args()
    out = args.out.resolve()
    base = (ROOT / "target/drawing-bench").resolve()
    if base not in out.parents:
        parser.error("--out must be a new directory under target/drawing-bench")
    out.mkdir(parents=True, exist_ok=False)
    metadata = {
        "status": "started",
        "python": sys.version,
        "executable": sys.executable,
        "executable_sha256": digest(Path(sys.executable).read_bytes()),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "head": command("git", "rev-parse", "HEAD"),
        "working_diff_sha256": digest(
            subprocess.check_output(["git", "diff", "HEAD", "--binary"], cwd=ROOT)
        ),
        "source_sha256": inventory(),
        "toolchain": {
            tool: command(tool, "--version") for tool in ("rustc", "cargo", "maturin")
        },
        "dot": DOT,
        "config": CONFIG,
        "fixture_sha256": digest(FIXTURE.read_bytes()),
        "warmups": 10,
        "sample_count": 21,
        "clock": vars(time.get_clock_info("perf_counter")),
        "timing_boundary": "installed diagram.render(config=CONFIG) only",
        "build_source_match": "unverified; module hashes identify installed artifacts",
        "native_route": "no labels/title; Scene.pages empty, Scene::render branch",
    }
    (out / "input.dot").write_text(DOT)
    (out / "config.json").write_text(json.dumps(CONFIG, indent=2) + "\n")
    try:
        result = run(out, metadata)
    except Exception as error:  # noqa: BLE001 -- redact arbitrary licensed host errors
        # License errors can contain private user data. Never serialize arbitrary
        # exception text or environment values into shared artifacts.
        metadata.update(status="blocked", error_type=type(error).__name__)
        if isinstance(error, ModuleNotFoundError):
            metadata["missing_module"] = error.name
            print(f"blocked: missing installed module {error.name}", file=sys.stderr)
        else:
            print(
                f"blocked: {type(error).__name__}; no performance result",
                file=sys.stderr,
            )
        result = 2
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return result


if __name__ == "__main__":
    sys.exit(main())
