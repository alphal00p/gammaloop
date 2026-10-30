"""Bounded four-domain certificate. No production generation or graph expansion.

Exit 2 is a checked rejection; exit 3 is a harness/setup error. Retained ELFs run
in new processes with a private runtime license. Nix is never invoked here.
"""

from __future__ import annotations

import argparse
import collections
from functools import cache
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import subprocess
import sys
import time

WORKSPACE = Path(__file__).resolve().parents[5]
BASE = WORKSPACE / "target/oracle-numerator-family-diagnosis"
HIST_WORKSPACE = WORKSPACE.parent / "generation-oracle-family-historical"
ROOTS = {
    "historical": HIST_WORKSPACE
    / "target/oracle-numerator-family-diagnosis/adapted-historical/captures",
    "current": BASE / "adapted-current/captures",
}
ELFS = {
    "historical": HIST_WORKSPACE
    / "target/oracle-numerator-family-diagnosis/adapted-historical/bin/oracle_family_exact",
    "current": BASE / "adapted-current/bin/oracle_family_exact",
}
CASES = ("qed_eemm_L1", "ctorder_qcd_uugg_L1", "sm_eemm_L0", "ym_g6_L0")
MODES = ("None", "Auto")
F = "spenso::{antisymmetric,real,spenso::tensor}::f"
SPELLINGS = {
    **{
        f"gammalooprs::{{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:{n.lower()}}}::{n}": f"gammalooprs::{{spenso::rank1,spenso::tensor,spenso::component-labels:subscript}}::{n}"
        for n in ("Q", "K", "P")
    },
    "gammalooprs::{spenso::index,spenso::index-label:e}::edge": "gammalooprs::{spenso::index}::edge",
    "gammalooprs::{spenso::index,spenso::index-label:v}::vertex": "gammalooprs::{spenso::index}::vertex",
}
TOKEN = re.compile(r"(?:[^\W\d]\w*::)+(?:\{[^{}]*\}::)?[^\W\d]\w*|𝑖|[0-9]+|[()+*/^,\-]")


class Reject(Exception):
    def __init__(self, stage, detail):
        super().__init__(detail)
        self.stage = stage


def require(condition, stage, detail):
    if not condition:
        raise Reject(stage, detail)


def load(path):
    def unique(items):
        result = {}
        for key, value in items:
            require(key not in result, "json", f"duplicate key {key}")
            result[key] = value
        return result

    try:
        return json.loads(Path(path).read_text(), object_pairs_hook=unique)
    except json.JSONDecodeError as error:
        raise Reject("json", str(error)) from error


def write(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


@cache
def parse(text, mapped=False):
    """Consume the complete captured syntax; retain every symbol definition token."""
    require(
        isinstance(text, str) and bool(text.strip()),
        "syntax",
        "empty/non-string expression",
    )
    tokens, position = [], 0
    while position < len(text):
        if text[position].isspace():
            position += 1
            continue
        match = TOKEN.match(text, position)
        require(
            match is not None,
            "syntax",
            f"unsupported syntax at {position}: {text[position : position + 30]}",
        )
        token = match.group()
        tokens.append(SPELLINGS.get(token, token) if mapped else token)
        position = match.end()
    cursor = 0

    def expression(minimum=0):
        nonlocal cursor
        require(cursor < len(tokens), "syntax", "truncated expression")
        token = tokens[cursor]
        cursor += 1
        if token == "-":
            left = ("neg", expression(30))
        elif token == "(":
            left = expression()
            require(
                cursor < len(tokens) and tokens[cursor] == ")",
                "syntax",
                "missing closing parenthesis",
            )
            cursor += 1
        elif token.isdecimal():
            left = ("integer", token)
        else:
            require("::" in token or token == "𝑖", "syntax", f"invalid primary {token}")
            if cursor < len(tokens) and tokens[cursor] == "(":
                cursor += 1
                args = []
                if cursor < len(tokens) and tokens[cursor] != ")":
                    while True:
                        args.append(expression())
                        if cursor >= len(tokens) or tokens[cursor] != ",":
                            break
                        cursor += 1
                require(
                    cursor < len(tokens) and tokens[cursor] == ")",
                    "syntax",
                    "invalid call termination",
                )
                cursor += 1
                left = ("call", token, tuple(args))
            else:
                left = ("symbol", token)
        while cursor < len(tokens):
            operator = tokens[cursor]
            precedence = {"+": 10, "-": 10, "*": 20, "/": 20, "^": 40}.get(operator, -1)
            if precedence < minimum:
                break
            cursor += 1
            right = expression(precedence if operator == "^" else precedence + 1)
            left = (operator, left, right)
        return left

    tree = expression()
    require(
        cursor == len(tokens),
        "syntax",
        f"trailing tokens {tokens[cursor : cursor + 3]}",
    )
    return tree


@cache
def render(tree):
    kind = tree[0]
    if kind in ("integer", "symbol"):
        return tree[1]
    if kind == "call":
        return tree[1] + "(" + ",".join(render(x) for x in tree[2]) + ")"
    if kind == "neg":
        return "-(" + render(tree[1]) + ")"
    return "(" + render(tree[1]) + kind + render(tree[2]) + ")"


def definitions(tree, result=None):
    if result is None:
        result = {}
    if tree[0] in ("symbol", "call"):
        spelling = tree[1]
        if spelling != "𝑖":
            name = re.sub(r"::\{[^{}]*\}::", "::", spelling)
            require(
                name not in result or result[name] == spelling,
                "symbol-definition",
                f"inconsistent definition {name}",
            )
            result[name] = spelling
    children = (
        tree[2]
        if tree[0] == "call"
        else tree[1:]
        if tree[0] not in ("symbol", "integer")
        else ()
    )
    for child in children:
        definitions(child, result)
    return result


@cache
def cyclic(tree):
    if tree[0] == "call":
        args = tuple(cyclic(x) for x in tree[2])
        if tree[1] == F:
            require(len(args) == 3, "syntax", "f arity")
            args = min(args[i:] + args[:i] for i in range(3))
        return ("call", tree[1], args)
    if tree[0] in ("integer", "symbol"):
        return tree
    return (tree[0], *(cyclic(x) for x in tree[1:]))


def pair(left, right, label, mapped=False):
    a, b = parse(left, mapped), parse(right)
    require(definitions(a) == definitions(b), "symbol-definition", label)
    return {
        "historical": render(a),
        "current": render(b),
        "label": label,
        "cyclic_f_literal": cyclic(a) == cyclic(b),
    }


def owners(payload):
    result = {}
    for i, vertex in enumerate(payload["vertices"]):
        result[f"vertex:{i}"] = (vertex["numerator"], vertex.get("num_interface"))
    for i, edge in enumerate(payload["edges"]):
        require(isinstance(edge, list) and len(edge) == 2, "edge-tuple", str(i))
        require(
            isinstance(edge[0], dict) and set(edge[0]) == {"source", "target"},
            "edge-tuple",
            str(i),
        )
        result[f"edge:{i}"] = (edge[1]["numerator"], edge[1].get("num_interface"))
    for name, interface in (
        ("numerator_prefactor", "numerator_prefactor_interface"),
        ("projector", "projector_interface"),
        ("numerator", "num_interface"),
    ):
        result[name] = (payload[name], payload.get(interface))
    return result


def permutation(left, right):
    require(
        len(left) == len(set(left)) and len(right) == len(set(right)),
        "port",
        "repeated port",
    )
    require(
        collections.Counter(left) == collections.Counter(right),
        "port",
        "typed port occurrence/dimension/variance mismatch",
    )
    return [right.index(x) for x in left]


def same_shape(actual, anchor, path=""):
    require(type(actual) is type(anchor), "schema", path)
    if isinstance(anchor, dict):
        require(actual.keys() == anchor.keys(), "schema", path)
        for key in anchor:
            same_shape(actual[key], anchor[key], path + "/" + key)
    elif isinstance(anchor, list):
        require(len(actual) == len(anchor), "schema", path)
        for i, (x, y) in enumerate(zip(actual, anchor)):
            same_shape(x, y, path + f"/{i}")


def ct_provenance(model, payload, request):
    def degree(tree):
        if tree[0] == "symbol":
            return int(tree[1].startswith("UFO::") and tree[1].endswith("_CT"))
        if tree[0] in ("integer", "call"):
            return 0
        if tree[0] == "neg":
            return degree(tree[1])
        left, right = degree(tree[1]), degree(tree[2])
        if tree[0] in ("+", "-"):
            require(left == right, "ct-order", "inhomogeneous CT sum")
            return left
        if tree[0] == "*":
            return left + right
        if tree[0] == "/":
            return left - right
        if left == 0:
            return 0
        require(
            tree[0] == "^" and tree[2][0] == "integer",
            "ct-order",
            "unsupported CT power",
        )
        return left * int(tree[2][1])

    anti = {p["name"]: p.get("antiname") or p["name"] for p in model["particles"]}
    base = lambda name: min(name, anti[name])
    rules, counts, rows = model["vertex_rules"], collections.Counter(), []
    coupling_orders = {c["name"]: dict(c["orders"]) for c in model["couplings"]}
    require(
        request["orders"]["plain"]["CT"] == {"inclusive": {"min": 1, "max": 1}},
        "ct-order",
        "resolved order",
    )
    for graph in payload["diagrams"]:
        p = graph["payload"]
        incidence = [collections.Counter() for _ in p["vertices"]]
        for ends, data in p["edges"]:
            particle = model["particles"][data["particle"]]["name"]
            for side in ("source", "target"):
                if ends[side] is not None:
                    incidence[ends[side]][base(particle)] += 1
        selected = []
        for i, vertex in enumerate(p["vertices"]):
            rule = rules[vertex["interaction"]]
            require(
                incidence[i] == collections.Counter(base(x) for x in rule["particles"]),
                "ct-incidence",
                graph["id"],
            )
            allowed = {c for row in rule["couplings"] for c in row if c is not None}
            found = {
                name.rsplit("::", 1)[-1]
                for name in definitions(parse(vertex["numerator"]))
                if name.startswith("UFO::") and name.rsplit("::", 1)[-1].endswith("_CT")
            }
            if rule["name"].endswith("_CT"):
                require(found and found <= allowed, "ct-coupling", graph["id"])
                require(
                    all(coupling_orders[c].get("CT") == 1 for c in found),
                    "ct-order",
                    graph["id"],
                )
                require(
                    degree(parse(vertex["numerator"])) == 1,
                    "ct-order",
                    "vertex CT degree",
                )
                counts[rule["name"]] += 1
                selected.append(
                    {
                        "vertex": i,
                        "interaction": vertex["interaction"],
                        "rule": rule["name"],
                        "couplings": sorted(found),
                    }
                )
            else:
                require(not found, "ct-coupling", "CT coupling at ordinary vertex")
        require(len(selected) == 1, "ct-order", graph["id"])
        require(degree(parse(p["numerator"])) == 1, "ct-order", "aggregate CT degree")
        rows.append({"id": graph["id"], "selected": selected})
    return {"rule_counts": dict(counts), "occurrences": rows}


def validate(case, mode, roots, output):
    data, pairs, port_records, ct = {}, [], [], {}
    # Frozen anchors bind observed graph/source data and the actual raw producer
    # decoding. Mathematical expressions are checked separately by the ELF.
    for route in ROOTS:
        folder = roots[route] / case / f"classes-{mode}-t1"
        anchor_folder = ROOTS[route] / case / f"classes-{mode}-t1"
        for name in ("model.json", "request.json", "options.json", "grouping.json"):
            same_shape(load(folder / name), load(anchor_folder / name), name)
            require(
                load(folder / name) == load(anchor_folder / name),
                "input-identity",
                route + "/" + name,
            )
        p, t = load(folder / "payload.json"), load(folder / "typed-owners.json")
        ap, at = (
            load(anchor_folder / "payload.json"),
            load(anchor_folder / "typed-owners.json"),
        )
        same_shape(p, ap, "payload")
        same_shape(t, at, "typed-owners")
        require(set(p) == {"diagrams", "groups"}, "schema", "top-level")
        require(
            len(t) == len(p["diagrams"]) == len(ap["diagrams"]), "owner-domain", route
        )
        ids = [g["id"] for g in p["diagrams"]]
        require(len(set(ids)) == len(ids), "owner-domain", "duplicate payload id")
        tm = {g["id"]: g for g in t}
        require(len(tm) == len(t) and set(tm) == set(ids), "owner-domain", route)
        require(p["groups"] == ap["groups"], "groups", route)
        require(len(p["groups"]) == len(ids), "groups", "None singleton count")
        defs = {}
        for index, (g, ag) in enumerate(zip(p["diagrams"], ap["diagrams"])):
            require(set(g) == set(ag), "schema", "diagram")
            for field in set(g) - {"payload"}:
                require(
                    g[field] == ag[field],
                    "physical-capture",
                    f"{route}/{index}/{field}",
                )
            payload, anchor_payload = g["payload"], ag["payload"]
            # Interface arrays can be transposed; compare their entries below.
            for field in set(payload) | set(anchor_payload):
                if field not in {
                    "vertices",
                    "edges",
                    "numerator",
                    "numerator_prefactor",
                    "projector",
                    "num_interface",
                    "numerator_prefactor_interface",
                    "projector_interface",
                }:
                    require(
                        payload.get(field) == anchor_payload.get(field),
                        "physical-capture",
                        field,
                    )
            require(set(payload) == set(anchor_payload), "schema", "payload")
            po, ao = owners(payload), owners(anchor_payload)
            for collection in ("vertices", "edges"):
                require(
                    len(payload[collection]) == len(anchor_payload[collection]),
                    "owner-domain",
                    collection,
                )
                for a, b in zip(payload[collection], anchor_payload[collection]):
                    if collection == "edges":
                        require(a[0] == b[0], "physical-capture", "incidence")
                        a, b = a[1], b[1]
                    require(set(a) == set(b), "schema", collection)
                    for field in set(a) - {"numerator", "num_interface"}:
                        require(a[field] == b[field], "physical-capture", field)
            typed = tm[g["id"]]
            require(set(typed) == {"id", "owners"}, "schema", "typed diagram")
            om = {o["owner"]: o for o in typed["owners"]}
            anchor_om = {o["owner"]: o for o in at[index]["owners"]}
            require(
                len(om) == len(typed["owners"]) and om.keys() == po.keys(),
                "owner-domain",
                g["id"],
            )
            for owner, (expression, interface) in po.items():
                o, anchor = om[owner], anchor_om[owner]
                require(set(o) == set(anchor), "schema", owner)
                require(
                    o["is_metric"] == anchor["is_metric"]
                    and o["is_composite"] == anchor["is_composite"],
                    "owner-flags",
                    owner,
                )
                require(
                    expression == o["producer_decoded_expression"],
                    "payload-live-link",
                    owner,
                )
                if route == "historical":
                    require(
                        o["expression"] == anchor["expression"],
                        "raw-producer-binding",
                        owner,
                    )
                else:
                    require(
                        o["expression"] == expression, "raw-producer-binding", owner
                    )
                tree = parse(expression, route == "historical")
                for name, definition in definitions(tree).items():
                    require(
                        name not in defs or defs[name] == definition,
                        "symbol-definition",
                        name,
                    )
                    defs[name] = definition
                ports = [
                    render(parse(x, route == "historical")) for x in o["ordered_ports"]
                ]
                anchor_ports = [
                    render(parse(x, route == "historical"))
                    for x in anchor["ordered_ports"]
                ]
                perm = permutation(anchor_ports, ports)
                if interface is not None:
                    require(
                        interface == o["ordered_ports"], "interface-live-link", owner
                    )
                port_records.append(
                    {"route": route, "id": g["id"], "owner": owner, "permutation": perm}
                )
                pairs.append(
                    pair(
                        render(
                            parse(
                                anchor["producer_decoded_expression"],
                                route == "historical",
                            )
                        ),
                        render(tree),
                        f"anchor/{route}/{g['id']}/{owner}",
                    )
                )
            group = p["groups"][index]
            require(
                group
                == {
                    "master": index,
                    "members": [
                        {
                            "id": g["id"],
                            "name": g["name"],
                            "overall_factor": g["overall_factor"],
                            "ratio": "1",
                            "source_index": index,
                        }
                    ],
                },
                "group-chain",
                g["id"],
            )
            factors = (
                [payload["numerator_prefactor"]]
                + [v["numerator"] for v in payload["vertices"]]
                + [e[1]["numerator"] for e in payload["edges"]]
            )
            pairs.append(
                pair(
                    render(
                        parse(
                            "*".join("(" + f + ")" for f in factors),
                            route == "historical",
                        )
                    ),
                    render(parse(payload["numerator"], route == "historical")),
                    f"aggregate-product/{route}/{g['id']}",
                )
            )
        if case == "ctorder_qcd_uugg_L1":
            ct[route] = ct_provenance(
                load(folder / "model.json"), p, load(folder / "request.json")
            )
        data[route] = (p, tm, defs)
    hp, ht, hd = data["historical"]
    cp, ct_owners, cd = data["current"]
    require(hd == cd, "symbol-definition", "complete cross-route definition table")
    hkeys = [json.dumps(g["key"], sort_keys=True) for g in hp["diagrams"]]
    ckeys = [json.dumps(g["key"], sort_keys=True) for g in cp["diagrams"]]
    require(
        collections.Counter(hkeys) == collections.Counter(ckeys)
        and len(set(hkeys)) == len(hkeys),
        "graph-domain",
        case,
    )
    require(hp["groups"] == cp["groups"], "groups", "cross-route")
    for h, c in zip(hp["diagrams"], cp["diagrams"]):
        require(h["id"] == c["id"], "graph-domain", "source order")
        hom = {o["owner"]: o for o in ht[h["id"]]["owners"]}
        com = {o["owner"]: o for o in ct_owners[c["id"]]["owners"]}
        for owner in hom:
            a, b = hom[owner], com[owner]
            permutation(
                [render(parse(x, True)) for x in a["ordered_ports"]],
                [render(parse(x)) for x in b["ordered_ports"]],
            )
            pairs.append(
                pair(
                    a["producer_decoded_expression"],
                    b["producer_decoded_expression"],
                    f"cross/{h['id']}/{owner}",
                    True,
                )
            )
    write(output / "pairs.json", pairs)
    all_definitions = {}
    for p in pairs:
        definitions(parse(p["historical"]), all_definitions)
        definitions(parse(p["current"]), all_definitions)
    result = {
        "case": case,
        "mode": mode,
        "graphs": len(hkeys),
        "pairs": len(pairs),
        "definitions": cd,
        "ct": ct,
        "ports": port_records,
        "cyclic_f_literal_pairs": sum(p["cyclic_f_literal"] for p in pairs),
        "raw_historical_decoding": "capture-bound to producer raw expression; decoded value anchored through semantic pairs",
    }
    write(output / "structural.json", result)
    return result


def run_exact(pairs_file, output, reverse=False):
    pairs = load(pairs_file)
    # Literal AST identity includes unit and zero owners. The round-one ELF's
    # post-parse power control cannot distinguish 1 from 1^2; only nonliteral
    # identities need its algebra normalization. Keep the complete table bound
    # in pairs.json and record the exact selected subset in a separate file.
    pairs = [p for p in pairs if p["historical"] != p["current"]]
    pairs_file = output / "pairs-nonliteral.json"
    write(pairs_file, pairs)
    if reverse:
        pairs = [
            {**p, "historical": p["current"], "current": p["historical"]} for p in pairs
        ]
        pairs_file = output / "pairs-reversed.json"
        write(pairs_file, pairs)
    env = dict(os.environ)
    for key in list(env):
        if key.startswith("FEYNKIT_ORACLE_") or key == "CARGO_TARGET_DIR":
            env.pop(key)
    if not env.get("SYMBOLICA_LICENSE"):
        config = Path(env.get("XDG_CONFIG_HOME", str(Path.home() / ".config")))
        env["SYMBOLICA_LICENSE"] = (config / "symbolica/license").read_text().strip()
    rows = []
    for kernel, elf in ELFS.items():
        result = output / f"results-{kernel}-{'reverse' if reverse else 'forward'}.json"
        log = result.with_suffix(".log")
        argv = [str(elf), str(pairs_file.resolve()), str(result.resolve())]
        proc = subprocess.run(argv, env=env, capture_output=True, text=True)
        log.write_text(proc.stdout + proc.stderr)
        row = {
            "kernel": kernel,
            "argv": argv,
            "exit": proc.returncode,
            "log": str(log),
            "result": str(result),
            "selected_pairs": len(pairs),
        }
        rows.append(row)
        write(output / f"commands-{'reverse' if reverse else 'forward'}.json", rows)
    require(all(row["exit"] == 0 for row in rows), "exact-checker", json.dumps(rows))
    for row in rows:
        records = load(row["result"])
        require(
            len(records) == len(pairs) and all(all(r.values()) for r in records),
            "exact-result",
            row["kernel"],
        )
        row["result_sha256"] = hashlib.sha256(
            Path(row["result"]).read_bytes()
        ).hexdigest()
    write(output / f"commands-{'reverse' if reverse else 'forward'}.json", rows)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=CASES, required=True)
    parser.add_argument("--mode", choices=MODES, required=True)
    parser.add_argument("--historical", type=Path, default=ROOTS["historical"])
    parser.add_argument("--current", type=Path, default=ROOTS["current"])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--structural-only", action="store_true")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    execution = {
        "argv": sys.argv,
        "validator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "python": sys.executable,
    }
    write(args.out / "execution.json", execution)
    try:
        validate(
            args.case,
            args.mode,
            {"historical": args.historical, "current": args.current},
            args.out,
        )
        if not args.structural_only:
            run_exact(args.out / "pairs.json", args.out)
            run_exact(args.out / "pairs.json", args.out, True)
        write(
            args.out / "verdict.json",
            {
                "disposition": "accept",
                "stage": "complete",
                "case": args.case,
                "mode": args.mode,
            },
        )
    except Reject as error:
        write(
            args.out / "verdict.json",
            {
                "disposition": "reject",
                "stage": error.stage,
                "detail": str(error),
                "case": args.case,
                "mode": args.mode,
            },
        )
        sys.exit(2)
    except Exception as error:
        write(
            args.out / "verdict.json",
            {
                "disposition": "harness-error",
                "exception": type(error).__name__,
                "detail": str(error),
            },
        )
        raise
    finally:
        execution.update(
            elapsed_seconds=time.monotonic() - started,
            max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        )
        write(args.out / "execution.json", execution)


if __name__ == "__main__":
    main()
