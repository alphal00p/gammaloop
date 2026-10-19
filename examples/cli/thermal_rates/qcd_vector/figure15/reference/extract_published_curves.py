"""Extract Figure 15 vector vertices from the pinned paper, without interpolation."""

import argparse
from collections import Counter
import hashlib
import json
import re
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("paper", type=Path)
parser.add_argument("output", type=Path)
parser.add_argument("--mutool", default="mutool")
args = parser.parse_args()
paper_sha256 = hashlib.file_digest(args.paper.open("rb"), "sha256").hexdigest()
# Tick positions below belong to this PDF, not arbitrary revisions of the paper.
expected_paper = "b92fd775f7b5ebb15ea382925acc6e89b6ff1f1f9eec1898432a319c259a7d65"
if paper_sha256 != expected_paper:
    parser.error("paper checksum differs from the calibrated arXiv:1910.07552v2 PDF")

with tempfile.TemporaryDirectory(prefix="figure15-vectors-") as directory:
    svg = Path(directory) / "page.svg"
    subprocess.run(
        [args.mutool, "draw", "-F", "svg", "-o", str(svg), str(args.paper), "33"],
        check=True,
    )
    svg_bytes = svg.read_bytes()
document = ET.fromstring(svg_bytes)

# All relevant paths share the same affine page transform. Work before that
# transform, so tick calibration and curve vertices use identical coordinates.
anchors = {
    "E_over_T": [[0, 71.8672], [0.5, 128.559], [1, 185.254], [1.5, 241.945], [2, 298.637]],
    "T": [[0, 265.566], [0.5, 312.813], [1, 360.059]],
    "L": [[-1, 71.8672], [-0.5, 119.109], [0, 166.355], [0.5, 213.598]],
}
calibration = {}
for axis, ticks in anchors.items():
    mean_value = sum(value for value, _ in ticks) / len(ticks)
    mean_coordinate = sum(coordinate for _, coordinate in ticks) / len(ticks)
    scale = sum((value - mean_value) * (coordinate - mean_coordinate) for value, coordinate in ticks) / sum(
        (value - mean_value) ** 2 for value, _ in ticks
    )
    origin = mean_coordinate - scale * mean_value
    calibration[axis] = {
        "ticks": ticks,
        "origin": origin,
        "units_per_value": scale,
        "maximum_tick_residual": max(abs(coordinate - origin - scale * value) for value, coordinate in ticks),
    }

colors = {"#ff0000": 0.5, "#800080": 1.0, "#0000ff": 1.5, "#ffa500": 0.0}
gray_index = Counter()
curves = {}
token_pattern = re.compile(r"[A-Za-z]|[-+]?(?:\d*\.\d+|\d+\.?)(?:[eE][-+]?\d+)?")
for path_index, path in enumerate(document.iter("{http://www.w3.org/2000/svg}path")):
    stroke = path.get("stroke")
    data = path.get("d", "")
    if stroke not in {*colors, "#d3d3d3"} or len(data) < 200:
        continue
    if path.get("transform") != "matrix(.9,0,0,-.9,163.855,443.832)":
        raise ValueError("unexpected curve coordinate transform")
    tokens = token_pattern.findall(data)
    if "".join(tokens) != re.sub(r"[\s,]", "", data):
        raise ValueError("unparsed SVG path syntax")
    points = []
    index = 0
    command = None
    x = y = 0.0
    while index < len(tokens):
        if tokens[index].isalpha():
            command = tokens[index]
            index += 1
        if command in ("M", "L"):
            if command == "M" and points:
                raise ValueError("unexpected interior move; retain it as a separate segment")
            x, y = float(tokens[index]), float(tokens[index + 1])
            index += 2
            command = "L"  # SVG's implicit lineto after the first moveto pair.
        elif command == "H":
            x = float(tokens[index])
            index += 1
        elif command == "V":
            y = float(tokens[index])
            index += 1
        else:
            raise ValueError(f"unsupported SVG command {command!r}")
        points.append([x, y])
    component = "T" if points[0][1] > 227.773 else "L"
    order = "lo" if stroke == "#d3d3d3" else "lo_plus_nlo_g3"
    if order == "lo":
        k = (0.5, 1.0, 1.5)[gray_index[component]]
        gray_index[component] += 1
    else:
        k = colors[stroke]
    energy_axis, ordinate_axis = calibration["E_over_T"], calibration[component]
    vertices = [
        [(x - energy_axis["origin"]) / energy_axis["units_per_value"],
         (y - ordinate_axis["origin"]) / ordinate_axis["units_per_value"]]
        for x, y in points
    ]
    key = (component, k, order)
    curve = curves.setdefault(key, {
        "component": component, "k_over_T": k, "order": order,
        "source_stroke": stroke, "source_dasharray": path.get("stroke-dasharray"),
        "segments": [], "source_paths": [],
    })
    curve["segments"].append(vertices)
    curve["source_paths"].append({
        "svg_path_index": path_index,
        "path_sha256": hashlib.sha256(data.encode()).hexdigest(),
        "vertex_count": len(vertices),
    })

assert len(curves) == 14 and sum(len(curve["segments"]) for curve in curves.values()) == 15
assert gray_index == {"T": 3, "L": 3}
assert len(curves[("T", 1.5, "lo_plus_nlo_g3")]["segments"]) == 2
assert all(len(segment) >= 75 for curve in curves.values() for segment in curve["segments"])
record = {
    "source": {
        "url": "https://arxiv.org/pdf/1910.07552v2",
        "title": "G. Jackson, Two-loop thermal spectral functions with general kinematics",
        "figure": 15, "page": 33, "paper_sha256": paper_sha256,
        "svg_sha256": hashlib.sha256(svg_bytes).hexdigest(),
        "mutool": args.mutool,
    },
    "calibration": calibration,
    "coordinate_definition": "E/T and rho_{T,L}/(eT)^2; full published LO+NLO uses g=3",
    "method": "Affine mapping of original vector vertices; no smoothing, interpolation, new quadrature, or reconnection of separate PDF paths",
    "limitation": "Published graphical precision and clipping; not an independent numerical integration target",
    "curves": list(curves.values()),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(record, indent=2) + "\n")
print(f"Extracted {len(curves)} curves, 15 paths, {sum(len(segment) for curve in curves.values() for segment in curve['segments'])} vertices")
