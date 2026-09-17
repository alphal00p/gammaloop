"""Reconstruct the retained GL638 spike from native logs, without GammaLoop calls.

The affine-star reconstruction follows the existing archived
analyze_max_native_trace.py. The CM graph routing is written explicitly here;
the retained registry confirms native cut-1 parent [3, 6, 7, 10].
Run with any Python 3 interpreter; only standard-library modules are used.
"""

import gzip
import json
import re
from decimal import Decimal as D
from decimal import getcontext
from pathlib import Path

getcontext().prec = 330
HERE = Path(__file__).resolve().parent
NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"


def table(text):
    return [
        [D(value.strip()) for value in line.split("│")[2:-1]]
        for line in text.splitlines()
        if line.startswith("│")
    ]


def root(text):
    tail = text[text.rfind("NewtonIterationResult") :]
    return list(map(D, re.findall("float: (" + NUMBER + ")", tail)[:3]))


def norm(vector):
    return sum(value * value for value in vector).sqrt()


def geometry(native):
    p, u, s, v = native  # q3, q6, q7, q10
    q = {3: p, 6: u, 7: s, 10: v, 12: p, 8: s, 5: v}
    q.update(
        {
            4: [u[a] + p[a] - v[a] for a in range(3)],
            2: [u[a] - v[a] for a in range(3)],
            13: [p[a] - v[a] for a in range(3)],
            14: [s[a] - p[a] for a in range(3)],
        }
    )
    energy = {
        e: (
            sum(a * a for a in vector)
            + D(125 if e == 2 else 0 if e in (13, 14) else 173) ** 2
        ).sqrt()
        for e, vector in q.items()
    }
    h = sum(energy[e] for e in (2, 4, 12)) - 600
    z = sum(energy[e] for e in (3, 10, 13)) - 600
    host = sum(energy[e] for e in (2, 6, 10)) - 600
    p_surface = 2 * energy[3] - 600
    radius = norm([h - host, p_surface * z / 600])
    return {
        "H_GeV": h,
        "Z_GeV": z,
        "host_residual_GeV": host,
        "R_HZ_GeV": norm([h, z]),
        "P_GeV": p_surface,
        "multiplier_denominator_radius_GeV": radius,
        "W_H": (h - host) ** 2 / radius**2,
        "E13_GeV": energy[13],
        "E14_GeV": energy[14],
        "energies_GeV": energy,
    }


records = [
    json.loads(line)
    for line in gzip.decompress((HERE / "native.jsonl.gz").read_bytes())
    .decode()
    .splitlines()
]
starts = [
    i for i, r in enumerate(records) if r.get("message", "").startswith("loop moms:")
]
assert len(starts) == 4, starts  # ordinary and forced-Arb, each identity and rotation
start, end = starts[2:]
raw = table(records[start]["message"])
roots = {}
pending = None
for row in records[start:end]:
    message = row.get("message", "")
    if message.startswith("representative esurface:"):
        pending = tuple(
            map(
                int,
                re.findall(
                    r"EdgeIndex\(\s*(\d+),?\s*\)", message.split("external_shift")[0]
                ),
            )
        )
    if message.startswith("solution:"):
        roots[pending] = root(message)
assert len(roots) == 6
t = roots[2, 6, 10][0]
p, r, s, v = [[value * t for value in vector] for vector in raw]
base = [p, [r[a] - p[a] + v[a] for a in range(3)], s, v]
centers = []
stars = []
for index in range(start, end):
    row = records[index]
    message = row.get("message", "")
    if row.get("stage") == "lu_threshold_center_values":
        assert row["rotation_id"] == "Identity rotation"
        centers.append((index, row, root(records[index + 2]["message"])))
    if not message.startswith("LU left evaluator input: cut_group_id=2"):
        continue
    local = int(re.search(r"left_threshold_id=(\d+)", message)[1])
    ci, center, alpha = next(
        c for c in reversed(centers) if c[1]["selected_esurface_id"] == local
    )
    active = list(
        map(int, re.findall(r"LoopIndex\((\d+)\)", center["subspace_loop_indices"]))
    )
    c = table(center["active_center"])
    fixed = table(center["center_with_fixed_complement"])
    join = max(
        abs(fixed[n][a] - (c[n][a] if n in active else base[n][a]))
        for n in range(4)
        for a in range(3)
    )
    assert join < D("1e-290"), join
    radius = norm([base[n][a] - c[n][a] for n in active for a in range(3)])
    measured_radius = D(re.search(r", r=(" + NUMBER + ")", message)[1])
    measured_star = D(re.search(r", rstar=(" + NUMBER + ")", message)[1])
    assert abs(radius - measured_radius) < D("1e-290")
    assert abs(radius * alpha[0] - measured_star) < D("1e-290")
    native_star = [
        [
            c[n][a] + alpha[0] * (base[n][a] - c[n][a]) if n in active else base[n][a]
            for a in range(3)
        ]
        for n in range(4)
    ]
    preceding = records[index - 10 : index]
    edges_row = next(
        row
        for row in preceding
        if row.get("message", "").startswith("edges in esurface:")
    )
    edges = list(map(int, re.findall(r'value: "e(\d+)"', edges_row["message"])))
    geo = geometry(native_star)
    residual = sum(geo["energies_GeV"][e] for e in edges) - 600
    assert abs(residual) < D("1e-290"), (edges, residual)
    stars.append(
        {
            "local_id": local,
            "edges": edges,
            "center_trace_index": ci,
            "evaluator_trace_index": index,
            "active_slots": active,
            "alpha": alpha[0],
            "radius_GeV": radius,
            "rstar_GeV": measured_star,
            "d_eta_dr": alpha[1] / radius,
            "join_error_GeV": join,
            "selected_surface_residual_GeV": residual,
            "geometry": geo,
        }
    )

normal = json.loads((HERE / "inspect.json").read_text())["evaluation"]
arb = json.loads((HERE / "inspect_arb.json").read_text())["evaluation"]
snapshot = json.loads((HERE / "snapshot_iter_0231.json").read_text())["slots"][0]
extrema = {
    row["component"]: row for row in snapshot["max_weight_info"] if row["sign"] == "-"
}
factors = {
    phase: extrema[phase]["max_eval"] / normal["integrand_result"][phase]
    for phase in ("re", "im")
}
assert abs(factors["re"] / factors["im"] - 1) < 2e-13
relative = {
    phase: abs(normal["integrand_result"][phase] / arb["integrand_result"][phase] - 1)
    for phase in ("re", "im")
}
assert max(relative.values()) < 2e-12
weight = factors[
    "re"
]  # Inferred from replay: the JSON snapshot omits the original outer weight.
cuts = [
    {
        "cut_id": event["cut_info"]["cut_id"],
        "unit_outer_weight": event["weight"],
        "reconstructed_MC_weight_pb": {
            k: v * weight for k, v in event["weight"].items()
        },
    }
    for event in normal["event_groups"][0]
]
terms = []
cut = None
for line in (
    gzip.decompress((HERE / "detailed_inspect.log.gz").read_bytes())
    .decode()
    .splitlines()
):
    cells = [cell.strip() for cell in line.split("│")[1:-1]]
    if len(cells) == 2 and cells[0] == "cut":
        cut = int(cells[1])
    if len(cells) == 7 and cells[0].isdigit():
        values = list(map(float, re.findall(NUMBER, cells[5])))
        terms.append(
            {
                "component_id": int(cells[0]),
                "cut_id": cut,
                "effective_multiplier": cells[3],
                "unit_outer_weight": values,
                "reconstructed_MC_weight_pb": [value * weight for value in values],
            }
        )
report = {
    "scope": "One retained MC maximum; no asymptotic or global-variance claim. Native star coordinates use forced-Arb identity trace; CT weights use ordinary detailed replay. Outer MC weight is inferred independently from Re and Im, not read from a retained Sample.",
    "extrema": extrema,
    "normal_result": normal["integrand_result"],
    "forced_arb_result": arb["integrand_result"],
    "relative_difference": relative,
    "outer_weight_inferred": factors,
    "physical_cut1": geometry(base),
    "cut1_stars": stars,
    "cuts": cuts,
    "components": terms,
}
(HERE / "analysis.json").write_text(json.dumps(report, default=str, indent=2) + "\n")
print(
    "Verified six physical roots, eight native cut-1 stars, and normal/Arb agreement."
)
print("Wrote", HERE / "analysis.json")
