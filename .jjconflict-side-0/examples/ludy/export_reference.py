"""Export physical inputs and separate oracles in the pinned LUDY environment."""

import argparse
import ast
import json
import subprocess
import sys
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("reference", type=Path)
parser.add_argument(
    "output", type=Path, help="destination directory for both JSON files"
)
args = parser.parse_args()
root = args.reference.resolve()
revision = "e5554019e398bbbc114ae3dafedc88bf09cb6eeb"
assert (
    subprocess.check_output(
        ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
    ).strip()
    == revision
)
subprocess.run(
    [
        "git",
        "-C",
        str(root),
        "diff",
        "--exit-code",
        "HEAD",
        "--",
        "calculations",
        "notebooks",
    ],
    check=True,
    stdout=subprocess.DEVNULL,
)
args.output.mkdir(parents=True, exist_ok=True)
sys.path[:0] = [str(root), str(root / "notebooks")]
from calculations.dis.python import physical as dis
from calculations.dy.python.physical import CutMembership, DYContext, DYGraphSpec
from notebooks.ludy import (
    gamma_trace_to_dots,
    momentum_printer,
    normalize_rational,
    replace_all,
    scalar_printer,
)
from symbolica import E, Expression, S
from symbolica.community.spenso import Representation, TensorName, dot

# These are the explicit dependencies of the four pinned input cells.
ns = {
    "E": E,
    "S": S,
    "Expression": Expression,
    "Representation": Representation,
    "TensorName": TensorName,
    "dot": dot,
    "gamma_trace_to_dots": gamma_trace_to_dots,
    "normalize_rational": normalize_rational,
    "replace_all": replace_all,
    "momentum_printer": momentum_printer,
    "scalar_printer": scalar_printer,
    "DYContext": DYContext,
    "DYGraphSpec": DYGraphSpec,
    "CutMembership": CutMembership,
}
tree = ast.parse((root / "calculations/dy/python/notebook.py").read_text())
for target in (
    "dy_context",
    "dy_graphs",
    "dy_trace_specs",
    "dy_qg_ward_trace_specs",
    "_expanded_numerator_oracles",
):
    cell = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef)
        and any(
            isinstance(child, ast.Name)
            and isinstance(child.ctx, ast.Store)
            and child.id == target
            for child in ast.walk(node)
        )
    )
    statements = []
    for stmt in cell.body:
        if isinstance(stmt, ast.Return):
            break
        statements.append(stmt)
        if target == "_expanded_numerator_oracles" and any(
            isinstance(child, ast.Name)
            and isinstance(child.ctx, ast.Store)
            and child.id == target
            for child in ast.walk(stmt)
        ):
            break
    exec(  # noqa: S102 -- explicitly requested, pinned local reference source
        compile(
            ast.Module(body=statements, type_ignores=[]),
            str(root / "calculations/dy/python/notebook.py"),
            "exec",
        ),
        ns,
    )

aliases = {
    name: S("ludy_port::" + name)
    for name in (
        "d",
        "kk",
        "kp1",
        "kp2",
        "p1_squared",
        "p2_squared",
        "p1_dot_p2",
        "mass_squared",
        "x",
        "Q_squared",
    )
}
# Scalar products contain d, so replace them before changing its namespace.
rules = [
    (ns[name], aliases[name])
    for name in (
        "kk",
        "kp1",
        "kp2",
        "p1_squared",
        "p2_squared",
        "p1_dot_p2",
        "mass_squared",
        "d",
    )
]
oracles = {"source_commit": revision, "dy": {}, "dy_ward": {}, "dis": {}}
for name, value in ns["_expanded_numerator_oracles"].items():
    assert normalize_rational(value - ns["dy_numerators"][name]) == 0, name
    oracles["dy"][name] = (
        replace_all(value, rules).together().cancel().factor().format_plain()
    )
for name, value in ns["dy_qg_ward_numerators"].items():
    oracles["dy_ward"][name] = (
        replace_all(value, rules).together().cancel().factor().format_plain()
    )

inputs = {"source_commit": oracles["source_commit"], "dy": {}, "dis": {}}
# Serialize only physical input; oracle expressions are kept in a separate file.
for name, graph in ns["dy_graphs"].items():
    spec = ns["dy_trace_specs"][name]
    vectors = []
    for vector in spec["vectors"]:
        if vector is None:
            vectors.append(None)
        else:
            match = next(
                (a, b, c)
                for a in (-1, 0, 1)
                for b in (-1, 0, 1)
                for c in (-1, 0, 1)
                if (a, b, c) != (0, 0, 0)
                and (a * ns["k"] + b * ns["p1"] + c * ns["p2"] - vector).expand() == 0
            )
            vectors.append(match)
    props = []
    for vector, mass in graph.propagators:
        match = next(
            (a, b, c)
            for a in (-1, 0, 1)
            for b in (-1, 0, 1)
            for c in (-1, 0, 1)
            if (a, b, c) != (0, 0, 0)
            and (a * ns["k"] + b * ns["p1"] + c * ns["p2"] - vector).expand() == 0
        )
        props.append({"momentum": match, "mass_squared": bool(mass != 0)})
    inputs["dy"][name] = {
        "vectors": vectors,
        "labels": spec["lorentz labels"],
        "propagators": props,
        "powers": graph.propagator_powers,
        "divisor_power": {"qg triangle": 1, "qg bubble": 2}.get(name, 0),
        "cuts": {
            key: {
                "mask": mask,
                "weight": str(graph.cut_weights[key]),
                "pm": graph.cut_memberships[key].pm_bare,
                "mp": graph.cut_memberships[key].mp,
                "analytic": key in graph.analytic_cuts,
            }
            for key, mask in graph.cutsets.items()
        },
        "observable_phase": str(graph.scheme_observable_phase),
    }
for name, graph in dis.DIS_GRAPH_SPECS.items():
    inputs["dis"][name] = {
        "vectors": [
            None if v is None else [v.k_coefficient, v.p_coefficient, v.q_coefficient]
            for v in graph.trace.vectors
        ],
        "labels": graph.trace.lorentz_labels,
        "photon_slots": graph.trace.photon_slots,
        "powers": graph.powers,
        "channel": graph.channel,
        "relative_weight": str(graph.relative_weight),
        "cuts": [
            {
                "name": c.name,
                "initial": c.initial_cut_lines,
                "final": c.final_cut_lines,
                "orders": c.cut_orders,
                "relative_sign": str(c.relative_sign),
            }
            for c in dis.KNOWN_DIS_CUTS
            if c.graph_name == name
        ],
    }
    nums = dis.derive_projected_numerators(graph)
    dis_rules = [
        (dis.d, aliases["d"]),
        (dis.kk, aliases["kk"]),
        (dis.kp, aliases["kp1"]),
        (dis.kq, aliases["kp2"]),
        (dis.p_squared, aliases["p1_squared"]),
        (dis.x, aliases["x"]),
        (dis.Q_squared, aliases["Q_squared"]),
    ]
    oracles["dis"][name] = {
        field: replace_all(getattr(nums, field), dis_rules)
        .together()
        .cancel()
        .factor()
        .format_plain()
        for field in ("metric", "longitudinal", "F1", "F2", "FL", "ward")
    }
(args.output / "inputs.json").write_text(json.dumps(inputs, indent=2) + "\n")
(args.output / "oracles.json").write_text(json.dumps(oracles, indent=2) + "\n")
print(
    "exported",
    len(inputs["dy"]),
    len(inputs["dis"]),
    sum(len(g["cuts"]) for g in inputs["dis"].values()),
)
