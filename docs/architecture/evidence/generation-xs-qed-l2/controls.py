"""Small new-case consumed-input controls; matched anchors isolate predicates."""

import copy
import importlib.util
import sys

import boundary as B
import common as M

SPEC = importlib.util.spec_from_file_location(
    "audited_mutations",
    M.HERE.parent / "generation-oracle-family-diagnosis/round2/controls.py",
)
R = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(R)

STAGES = {
    "missing-owner": "owner-domain",
    "extra-owner": "owner-domain",
    "duplicate-owner": "owner-domain",
    "owner-swap": "payload-live-link",
    "live-link": "payload-live-link",
    "raw-link": "raw-producer-binding",
    "rule": "cross-physical",
    "routing": "cross-physical",
    "source-factor": "cross-physical",
    "model": "cross-input",
    "request": "selected-mode",
    "group-ratio": "group-chain",
    "group-source": "group-chain",
    "scalar-type": "schema",
    "edge-tuple": "edge-tuple",
    "wrong-path": "schema",
    "port-dimension": "cross-ports",
    "port-variance": "cross-ports",
    "port-occurrence": "cross-ports",
    "port-transpose": "complete",
    "whole-key-reorder": "complete",
    "string-sign": "cross-equality",
    "string-coefficient": "cross-equality",
    "string-power": "cross-equality",
    "trailing-word": "syntax",
    "trailing-semicolon": "syntax",
    "gamma-attribute": "cross-definitions",
    "g-attribute": "cross-definitions",
    "hedge-tag": "cross-ports",
    "mirrored-aggregate": "aggregate-product",
    "projector-real": "native-parse",
    "internal-definition": "symbol-definition",
    "interface-link": "interface-live-link",
    "typed-flags": "cross-flags",
    "capture-interface": "capture-interface",
    "grouping": "selected-grouping",
    "extra-edge": "cross-owner-domain",
    "missing-edge": "cross-owner-domain",
    "edge-particle": "cross-physical",
    "bool-master": "schema",
    "anchor-model": "anchor-input",
}


def main():
    rows = []
    for mode in ("None", "Auto"):
        for name, stage in STAGES.items():
            root = M.OUT / "controls" / mode / name
            roots, docs = {}, {}
            for route, original in B.ROOTS.items():
                roots[route] = root / route
                folder = roots[route] / M.CASES[0] / f"classes-{mode}-t1"
                folder.mkdir(parents=True, exist_ok=False)
                for path in (original / M.CASES[0] / f"classes-{mode}-t1").iterdir():
                    (folder / path.name).symlink_to(path)
                docs[route] = {f: M.load(folder / f) for f in B.FILES[:6]}
            p, t, model, req = (
                docs["current"][f]
                for f in (
                    "payload.json",
                    "typed-owners.json",
                    "model.json",
                    "request.json",
                )
            )
            g = p["diagrams"][0]["payload"]
            if name in (
                "projector-real",
                "internal-definition",
                "interface-link",
                "typed-flags",
                "capture-interface",
                "grouping",
                "extra-edge",
                "missing-edge",
                "edge-particle",
                "bool-master",
                "anchor-model",
                "whole-key-reorder",
            ):
                if name == "projector-real":
                    for route, d in docs.items():
                        value = (
                            "("
                            + d["payload.json"]["diagrams"][0]["payload"]["projector"]
                            + ")*xs_probe::{ real}::bad"
                        )
                        d["payload.json"]["diagrams"][0]["payload"]["projector"] = value
                        o = next(
                            o
                            for o in d["typed-owners.json"][0]["owners"]
                            if o["owner"] == "projector"
                        )
                        o["producer_decoded_expression"] = value
                        if route == "current":
                            o["expression"] = value
                elif name == "internal-definition":
                    vertex = next(
                        (i, v)
                        for i, v in enumerate(g["vertices"])
                        if "linear" in v["numerator"]
                    )
                    i, v = vertex
                    value = v["numerator"].replace("linear,", "", 1)
                    M.require(value != v["numerator"], "control-setup", name)
                    v["numerator"] = value
                    o = next(o for o in t[0]["owners"] if o["owner"] == f"vertex:{i}")
                    o["expression"] = o["producer_decoded_expression"] = value
                elif name in ("interface-link", "capture-interface"):
                    obj = next(
                        v
                        for graph in p["diagrams"]
                        for v in graph["payload"]["vertices"]
                        if "num_interface" in v
                    )
                    if name == "interface-link":
                        obj["num_interface"].reverse()
                    else:
                        del obj["num_interface"]
                elif name == "typed-flags":
                    t[0]["owners"][0]["is_composite"] = not t[0]["owners"][0][
                        "is_composite"
                    ]
                elif name == "grouping":
                    for d in docs.values():
                        d["grouping.json"] = "identical"
                elif name in ("extra-edge", "missing-edge"):
                    if name == "extra-edge":
                        edge = copy.deepcopy(g["edges"][0])
                        edge[1]["numerator"] = "1"
                        edge[1].pop("num_interface", None)
                        index = len(g["edges"])
                        g["edges"].append(edge)
                        o = copy.deepcopy(
                            next(o for o in t[0]["owners"] if o["owner"] == "edge:0")
                        )
                        o.update(
                            owner=f"edge:{index}",
                            expression="1",
                            producer_decoded_expression="1",
                            ordered_ports=[],
                        )
                        t[0]["owners"].append(o)
                    else:
                        index = len(g["edges"]) - 1
                        g["edges"].pop()
                        t[0]["owners"] = [
                            o for o in t[0]["owners"] if o["owner"] != f"edge:{index}"
                        ]
                elif name == "edge-particle":
                    g["edges"][0][1]["particle"] += 1
                elif name == "bool-master":
                    for d in docs.values():
                        d["payload.json"]["groups"][0]["master"] = False
                elif name == "anchor-model":
                    for d in docs.values():
                        d["model.json"]["particles"][0]["mass"] = "xs_anchor_mass"
                elif name == "whole-key-reorder":

                    def reverse(x):
                        if isinstance(x, dict):
                            return {k: reverse(v) for k, v in reversed(list(x.items()))}
                        if isinstance(x, list):
                            return [reverse(v) for v in x]
                        return x

                    docs = reverse(docs)
            else:
                R.mutate(name, p, t, model, req)
            if name == "mirrored-aggregate":
                d = docs["historical"]
                pg = d["payload.json"]["diagrams"][0]["payload"]
                value = "2*(" + pg["numerator"] + ")"
                pg["numerator"] = value
                next(
                    o
                    for o in d["typed-owners.json"][0]["owners"]
                    if o["owner"] == "numerator"
                )["producer_decoded_expression"] = value
            for route, d in docs.items():
                for f, value in d.items():
                    path = roots[route] / M.CASES[0] / f"classes-{mode}-t1" / f
                    path.unlink()
                    M.write(path, value)
            # Capture-integrity tests intentionally keep frozen anchors;
            # mathematical controls mirror changed consumed inputs and anchors.
            anchors = (
                B.ROOTS if name in ("capture-interface", "anchor-model") else roots
            )
            argv = [
                sys.executable,
                "-B",
                M.HERE / "validate.py",
                "--mode",
                mode,
                "--out",
                root / "result",
            ]
            for route, value in roots.items():
                argv += ["--" + route, value, "--anchors-" + route, anchors[route]]
            proc, command = M.run(argv, root / "execution")
            verdict = M.load(root / "result/verdict.json")
            expected = "accept" if stage == "complete" else "reject"
            matched = (
                verdict["stage"] == stage
                and verdict["disposition"] == expected
                and proc.returncode == (0 if expected == "accept" else 2)
            )
            rows.append(
                {
                    "name": name,
                    "mode": mode,
                    "expected_stage": stage,
                    "expected_disposition": expected,
                    "observed": verdict,
                    "matched": matched,
                    "command": command,
                }
            )
            M.write(M.OUT / "controls/results.json", rows)
            print(mode, name, matched, verdict["stage"], flush=True)
    M.require(
        all(r["matched"] for r in rows),
        "controls",
        "unexpected predicate/disposition; no credit",
    )


if __name__ == "__main__":
    main()
