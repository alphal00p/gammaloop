"""Mutate fresh capture copies and invoke the certificate as a subprocess.

Mutation exceptions are harness failures, never validation rejections. No digest
or FNV controls contribute to the semantic verdict.
"""

import argparse
import copy
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


def mutate(name, payload, typed, model, request):
    p, t = payload, typed
    g = p["diagrams"][0]
    owner = t[0]["owners"][0]
    vertex = g["payload"]["vertices"][0]

    def change_expression(fn, aggregate=False):
        field = "numerator"
        obj = g["payload"] if aggregate else vertex
        path = "numerator" if aggregate else "vertex:0"
        target = next(o for o in t[0]["owners"] if o["owner"] == path)
        value = fn(obj[field])
        obj[field] = value
        target["producer_decoded_expression"] = value
        target["expression"] = value

    if name == "missing-owner":
        t[0]["owners"].pop()
    elif name in ("extra-owner", "duplicate-owner"):
        extra = copy.deepcopy(owner)
        if name == "extra-owner":
            extra["owner"] = "vertex:999"
        t[0]["owners"].append(extra)
    elif name == "owner-swap":
        a, b = t[0]["owners"][:2]
        a["owner"], b["owner"] = b["owner"], a["owner"]
    elif name == "extra-typed-domain":
        t.append({"id": "ffffffffffffffffffffffffffffffff", "owners": []})
    elif name == "live-link":
        owner["producer_decoded_expression"] = "0"
    elif name == "raw-link":
        owner["expression"] = "0"
    elif name in ("aggregate", "mirrored-aggregate"):
        change_expression(lambda s: "2*(" + s + ")", True)
    elif name == "wrong-path":
        g["payload"]["edges"][0][0]["numerator"] = "0"
    elif name == "edge-tuple":
        g["payload"]["edges"][0].append({"numerator": "0"})
    elif name == "rule":
        vertex["interaction"] += 1
    elif name == "source-factor":
        g["payload"]["overall_factor"] = "2*(" + g["payload"]["overall_factor"] + ")"
    elif name == "routing":
        signatures = g["payload"]["loop_momentum_basis"]["edge_signatures"]
        changed = False
        for signature in signatures.values():
            for field in ("loops", "external"):
                if any(signature[field]):
                    signature[field] = [-x for x in signature[field]]
                    changed = True
                    break
            if changed:
                break
        if not changed:
            raise RuntimeError("routing mutation requires nonzero signature")
    elif name == "model":
        model["particles"][0]["mass"] = "review_changed_mass"
    elif name == "request":
        request["classes"] = "none" if request["classes"] == "auto" else "auto"
    elif name in ("group-ratio", "group-source"):
        p["groups"][0]["members"][0][
            "ratio" if name == "group-ratio" else "source_index"
        ] = "2" if name == "group-ratio" else 999
    elif name == "scalar-type":
        p["groups"][1]["members"][0]["source_index"] = True
    elif name in (
        "port-transpose",
        "port-dimension",
        "port-variance",
        "port-occurrence",
        "interface-trailing",
    ):
        target = next(o for o in t[0]["owners"] if len(o["ordered_ports"]) >= 2)
        ports = target["ordered_ports"]
        if name == "port-transpose":
            ports.reverse()
        elif name == "port-dimension":
            ports[0] = (
                ports[0]
                .replace("(4,", "(5,", 1)
                .replace("(8,", "(9,", 1)
                .replace("(3,", "(6,", 1)
            )
        elif name == "port-variance":
            ports[0] = "spenso::{spenso::lower}::dind(" + ports[0] + ")"
        elif name == "port-occurrence":
            ports.append(ports[0])
        else:
            ports[0] += ";"
        kind, _, index = target["owner"].partition(":")
        if kind == "vertex":
            obj = g["payload"]["vertices"][int(index)]
        elif kind == "edge":
            obj = g["payload"]["edges"][int(index)][1]
        else:
            obj = g["payload"]
        field = {
            "numerator_prefactor": "numerator_prefactor_interface",
            "projector": "projector_interface",
        }.get(kind, "num_interface")
        if field in obj:
            obj[field] = list(ports)
    elif name in (
        "string-sign",
        "string-coefficient",
        "string-power",
        "trailing-term",
        "trailing-word",
        "trailing-paren",
        "trailing-semicolon",
        "truncated",
        "empty",
    ):
        functions = {
            "string-sign": lambda s: "-(" + s + ")",
            "string-coefficient": lambda s: "2*(" + s + ")",
            "string-power": lambda s: "(" + s + ")^2",
            "trailing-term": lambda s: s + "+review::{}::extra",
            "trailing-word": lambda s: s + " junk",
            "trailing-paren": lambda s: s + ")",
            "trailing-semicolon": lambda s: s + ";",
            "truncated": lambda s: s[: len(s) // 2],
            "empty": lambda s: "",
        }
        change_expression(functions[name])
    elif name in ("even-f", "odd-f", "metric-swap"):
        head = (
            V.F
            if name != "metric-swap"
            else "spenso::{symmetric,linear,real,spenso::tensor}::g"
        )
        found = False
        for gi, graph in enumerate(p["diagrams"]):
            for vi, vert in enumerate(graph["payload"]["vertices"]):
                if head + "(" not in vert["numerator"]:
                    continue

                def transform(tree):
                    nonlocal found
                    if tree[0] == "call":
                        args = tree[2]
                        if tree[1] == head and not found:
                            found = True
                            order = (
                                (1, 2, 0)
                                if name == "even-f"
                                else (1, 0, 2)
                                if name == "odd-f"
                                else (1, 0)
                            )
                            return ("call", head, tuple(args[i] for i in order))
                        return ("call", tree[1], tuple(transform(x) for x in args))
                    if tree[0] in ("symbol", "integer"):
                        return tree
                    return (tree[0], *(transform(x) for x in tree[1:]))

                value = V.render(transform(V.parse(vert["numerator"])))
                vert["numerator"] = value
                target = next(
                    o
                    for x in t
                    if x["id"] == graph["id"]
                    for o in x["owners"]
                    if o["owner"] == f"vertex:{vi}"
                )
                target["expression"] = target["producer_decoded_expression"] = value
                return
        raise RuntimeError("meaningful tensor permutation unavailable")
    elif name in (
        "f-attribute",
        "g-attribute",
        "gamma-attribute",
        "hedge-tag",
        "unmapped-old-Q",
    ):
        old, new = {
            "f-attribute": (V.F, "spenso::{real,spenso::tensor}::f"),
            "g-attribute": (
                "spenso::{symmetric,linear,real,spenso::tensor}::g",
                "spenso::{linear,real,spenso::tensor}::g",
            ),
            "gamma-attribute": (
                "spenso::{linear,spenso::tensor}::gamma",
                "spenso::{spenso::tensor}::gamma",
            ),
            "hedge-tag": (
                "gammalooprs::{spenso::index}::hedge",
                "gammalooprs::{}::hedge",
            ),
            "unmapped-old-Q": (
                next(v for k, v in V.SPELLINGS.items() if k.endswith("::Q")),
                next(k for k in V.SPELLINGS if k.endswith("::Q")),
            ),
        }[name]
        count = 0

        def replace(value):
            nonlocal count
            if isinstance(value, str):
                count += value.count(old)
                return value.replace(old, new)
            if isinstance(value, list):
                return [replace(x) for x in value]
            if isinstance(value, dict):
                return {k: replace(v) for k, v in value.items()}
            return value

        np, nt = replace(p), replace(t)
        p.clear()
        p.update(np)
        t[:] = nt
        if not count:
            raise RuntimeError("symbol mutation was a no-op")
    elif name not in ("whole-key-reorder", "malformed-json"):
        raise RuntimeError("unknown mutation " + name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--case", choices=V.CASES)
    ap.add_argument("--mode", choices=V.MODES)
    args = ap.parse_args()
    common = [
        "missing-owner",
        "extra-owner",
        "duplicate-owner",
        "owner-swap",
        "extra-typed-domain",
        "live-link",
        "raw-link",
        "aggregate",
        "mirrored-aggregate",
        "wrong-path",
        "edge-tuple",
        "rule",
        "source-factor",
        "routing",
        "model",
        "request",
        "group-ratio",
        "group-source",
        "scalar-type",
        "port-transpose",
        "port-dimension",
        "port-variance",
        "port-occurrence",
        "interface-trailing",
        "string-sign",
        "string-coefficient",
        "string-power",
        "trailing-term",
        "trailing-word",
        "trailing-paren",
        "trailing-semicolon",
        "truncated",
        "empty",
        "hedge-tag",
        "g-attribute",
        "metric-swap",
        "whole-key-reorder",
        "malformed-json",
    ]
    rows = []
    for case in (args.case,) if args.case else V.CASES:
        special = (
            ["even-f", "odd-f", "f-attribute"]
            if case in ("ym_g6_L0", "ctorder_qcd_uugg_L1")
            else ["gamma-attribute"]
        )
        if case in ("qed_eemm_L1", "ctorder_qcd_uugg_L1"):
            special += ["unmapped-old-Q"]
        for mode in (args.mode,) if args.mode else V.MODES:
            for name in common + special:
                root = args.out / case / mode / name
                root.mkdir(parents=True, exist_ok=False)
                originals = {}
                for route, source in V.ROOTS.items():
                    d = root / route / case / f"classes-{mode}-t1"
                    d.mkdir(parents=True)
                    for file in (source / case / f"classes-{mode}-t1").iterdir():
                        (d / file.name).symlink_to(file)
                    originals[route] = [
                        V.load(d / f)
                        for f in (
                            "payload.json",
                            "typed-owners.json",
                            "model.json",
                            "request.json",
                        )
                    ]
                before = json.dumps(originals["current"], sort_keys=True)
                # Exceptions in setup escape and fail the harness immediately.
                mutate(name, *originals["current"])
                if name == "mirrored-aggregate":
                    hp, ht, _, _ = originals["historical"]
                    value = "2*(" + hp["diagrams"][0]["payload"]["numerator"] + ")"
                    hp["diagrams"][0]["payload"]["numerator"] = value
                    next(o for o in ht[0]["owners"] if o["owner"] == "numerator")[
                        "producer_decoded_expression"
                    ] = value
                for route in V.ROOTS:
                    d = root / route / case / f"classes-{mode}-t1"
                    for file, value in zip(
                        (
                            "payload.json",
                            "typed-owners.json",
                            "model.json",
                            "request.json",
                        ),
                        originals[route],
                    ):
                        path = d / file
                        path.unlink()
                        if name == "whole-key-reorder":

                            def reverse(x):
                                if isinstance(x, dict):
                                    return {
                                        k: reverse(v)
                                        for k, v in reversed(list(x.items()))
                                    }
                                if isinstance(x, list):
                                    return [reverse(v) for v in x]
                                return x

                            value = reverse(value)
                        V.write(path, value)
                    if route == "current" and name == "malformed-json":
                        path = d / "payload.json"
                        path.write_text(path.read_text()[:-2])
                after = json.dumps(originals["current"], sort_keys=True)
                if (
                    name not in ("whole-key-reorder", "malformed-json")
                    and before == after
                ):
                    raise RuntimeError("mutation no-op: " + name)
                output = root / "result"
                argv = [
                    sys.executable,
                    "-I",
                    "-B",
                    str(Path(__file__).with_name("validate.py")),
                    "--case",
                    case,
                    "--mode",
                    mode,
                    "--historical",
                    str(root / "historical"),
                    "--current",
                    str(root / "current"),
                    "--out",
                    str(output),
                ]
                proc = subprocess.run(argv, capture_output=True, text=True)
                (root / "validator.log").write_text(proc.stdout + proc.stderr)
                verdict = V.load(output / "verdict.json")
                expected = (
                    "accept"
                    if name
                    in ("port-transpose", "even-f", "metric-swap", "whole-key-reorder")
                    else "reject"
                )
                row = {
                    "case": case,
                    "mode": mode,
                    "mutation": name,
                    "intended": expected,
                    "observed": verdict,
                    "argv": argv,
                    "exit": proc.returncode,
                    "matched": verdict["disposition"] == expected
                    and proc.returncode == (0 if expected == "accept" else 2),
                }
                V.write(root / "control.json", row)
                rows.append(row)
                V.write(args.out / "results.json", rows)
                print(
                    case, mode, name, row["matched"], verdict.get("stage"), flush=True
                )
    if not all(r["matched"] for r in rows):
        sys.exit(1)


if __name__ == "__main__":
    main()
