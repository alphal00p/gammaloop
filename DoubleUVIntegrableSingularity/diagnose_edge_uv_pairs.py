#!/usr/bin/env python3
"""Diagnose UV-degenerate edge pairs at aa_aa 2L max-weight points.

The diagnostic intentionally works from the generated dot files rather than the
binary GammaLoop state.  It evaluates the spatial edge momenta encoded in
``lmb_rep`` for selected max-weight points and reports edge pairs that become
identical or nearly identical up to fixed external shifts in the hard limit.
"""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GRAPH_DIR = REPO / "examples" / "cli" / "aa_aa" / "2L" / "graphs"
OUT = ROOT / "data" / "current_maxweight_diagnostics" / "edge_uv_pair_diagnostics.json"

E_CM = 300.0
EXTERNAL_SPATIAL = {
    0: (0.0, 0.0, 86.5),
    1: (0.0, 0.0, -86.5),
    2: (0.0, -79.27855952273603, -34.6),
}


@dataclass(frozen=True)
class Point:
    label: str
    graph: str
    graph_id: int
    xs: tuple[float, float, float, float, float, float]


POINTS = [
    Point(
        "GL06 short-run real-minus max",
        "GL06",
        4,
        (
            0.9995215227215548,
            0.5787928959151282,
            0.22059957686624845,
            0.9999299728315765,
            0.8989299227500787,
            0.42030419319463375,
        ),
    ),
    Point(
        "GL06 short-run prob1000 real/imag-plus max",
        "GL06",
        4,
        (
            0.99988492079262081,
            0.44226073563398893,
            0.29422178156572731,
            0.99995583397035226,
            0.25056458851596020,
            0.069123913846451557,
        ),
    ),
    Point(
        "GL14 short-run bins128 imag-minus max",
        "GL14",
        7,
        (
            0.99928560786268694,
            0.69209294959349577,
            0.28659450247430496,
            0.99985395518723508,
            0.11911699688086141,
            0.050153678028659998,
        ),
    ),
    Point(
        "GL18 short-run bins128 real-minus/imag-plus max",
        "GL18",
        10,
        (
            0.99972605511944601,
            0.59527066218796598,
            0.028121655386487698,
            0.99929042279752922,
            0.99432932400174268,
            0.86220910818747232,
        ),
    ),
    Point(
        "GL18 80M real-plus max",
        "GL18",
        10,
        (
            0.9999738238928433,
            0.6761869937119277,
            0.479166862842833,
            0.9999638043698932,
            0.017088467296021152,
            0.8935175085829773,
        ),
    ),
]


def add(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return tuple(x + y for x, y in zip(a, b))


def sub(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return tuple(x - y for x, y in zip(a, b))


def mul(c: float, a: tuple[float, float, float]) -> tuple[float, float, float]:
    return tuple(c * x for x in a)


def dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def norm(a: tuple[float, float, float]) -> float:
    return math.sqrt(dot(a, a))


def spherical_to_vec(xs: tuple[float, float, float]) -> tuple[float, float, float]:
    radius = E_CM * xs[0] / (1.0 - xs[0])
    phi = 2.0 * math.pi * xs[1]
    cos_theta = -1.0 + 2.0 * xs[2]
    sin_theta = math.sqrt(max(0.0, 1.0 - cos_theta * cos_theta))
    return (
        radius * sin_theta * math.cos(phi),
        radius * sin_theta * math.sin(phi),
        radius * cos_theta,
    )


def parse_lmb_reps(graph: str) -> dict[int, str]:
    text = (GRAPH_DIR / f"{graph}.dot").read_text()
    return {
        int(match.group(1)): match.group(2)
        for match in re.finditer(r'\[id=(\d+)[^\]]*lmb_rep="([^"]+)"', text)
    }


def evaluate_rep(
    rep: str, loop_momenta: dict[int, tuple[float, float, float]]
) -> tuple[float, float, float]:
    total = (0.0, 0.0, 0.0)
    for coeff, kind, idx in re.findall(r"(?:(-?\d+(?:\.\d+)?)\*)?([KP])\((\d+),", rep):
        factor = float(coeff) if coeff else 1.0
        vector = loop_momenta[int(idx)] if kind == "K" else EXTERNAL_SPATIAL[int(idx)]
        total = add(total, mul(factor, vector))
    return total


def diagnose(point: Point) -> dict[str, object]:
    loop_momenta = {
        0: spherical_to_vec(point.xs[:3]),
        1: spherical_to_vec(point.xs[3:]),
    }
    reps = parse_lmb_reps(point.graph)
    edges = {
        edge_id: evaluate_rep(rep, loop_momenta)
        for edge_id, rep in reps.items()
        if edge_id >= 4
    }

    pair_rows = []
    edge_ids = sorted(edges)
    for i, first_id in enumerate(edge_ids):
        for second_id in edge_ids[i + 1 :]:
            first = edges[first_id]
            second = edges[second_id]
            first_norm = norm(first)
            second_norm = norm(second)
            scale = max(first_norm, second_norm)
            if scale == 0.0:
                continue
            cosine = dot(first, second) / (first_norm * second_norm)
            pair_rows.append(
                {
                    "edges": [first_id, second_id],
                    "difference_over_max_norm": norm(sub(first, second)) / scale,
                    "sum_over_max_norm": norm(add(first, second)) / scale,
                    "cosine": cosine,
                    "norm_ratio": min(first_norm, second_norm) / scale,
                    "first_norm": first_norm,
                    "second_norm": second_norm,
                    "first_rep": reps[first_id],
                    "second_rep": reps[second_id],
                }
            )

    return {
        "label": point.label,
        "graph": point.graph,
        "graph_id": point.graph_id,
        "basis_norms": [norm(loop_momenta[0]), norm(loop_momenta[1])],
        "basis_cosine": dot(loop_momenta[0], loop_momenta[1])
        / (norm(loop_momenta[0]) * norm(loop_momenta[1])),
        "closest_pairs": sorted(pair_rows, key=lambda row: row["difference_over_max_norm"])[:8],
    }


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    diagnostics = [diagnose(point) for point in POINTS]
    OUT.write_text(json.dumps(diagnostics, indent=2))
    print(json.dumps(diagnostics, indent=2))


if __name__ == "__main__":
    main()
