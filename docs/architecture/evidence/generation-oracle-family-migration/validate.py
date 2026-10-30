"""Bind actual runtime dumps/observers to the approved finite capture domain."""

import collections
import json
import re

import common as M
import derive as D


def same(a, b, label):
    M.require(M.C.strict_equal(a, b), "runtime-binding", label)


def modes(a, b, ta, tb, case):
    key = lambda g: json.dumps(g["key"], sort_keys=True, separators=(",", ":"))
    ak, bk = [key(g) for g in a["diagrams"]], [key(g) for g in b["diagrams"]]
    M.require(
        collections.Counter(ak) == collections.Counter(bk) and len(set(ak)) == len(ak),
        "mode-domain",
        case,
    )
    mapping = []
    for i, g in enumerate(a["diagrams"]):
        j = bk.index(key(g))
        other = b["diagrams"][j]
        same(
            {k: v for k, v in g.items() if k != "name"},
            {k: v for k, v in other.items() if k != "name"},
            case + ": mode physical",
        )
        for p, graph, index in ((a, g, i), (b, other, j)):
            same(
                p["groups"][index],
                {
                    "master": index,
                    "members": [
                        {
                            "id": graph["id"],
                            "name": graph["name"],
                            "ratio": "1",
                            "overall_factor": graph["overall_factor"],
                            "source_index": index,
                        }
                    ],
                },
                case + ": source/group",
            )
        mapping.append(
            {
                "none": i,
                "auto": j,
                "id": g["id"],
                "key": g["key"],
                "none_name": g["name"],
                "auto_name": other["name"],
                "groups": [a["groups"][i], b["groups"][j]],
            }
        )
    same(
        {g["id"]: g["owners"] for g in ta},
        {g["id"]: g["owners"] for g in tb},
        case + ": typed modes",
    )
    return mapping


def debug_row(text):
    body = re.sub(
        r"Fraction \{ numerator: (-?\d+), denominator: (\d+) \}", r"(\1, \2)", text
    )
    body = body.replace("symmetry_factors: [", "symmetry_factors: &[")
    return D.row("debug: " + body.replace("Summary {", "Row {", 1), "debug")


def main():
    derivation = M.load(M.OUT / "derivation.json")
    approved_order = M.load(
        M.WORKSPACE
        / "target/oracle-numerator-family-round3/qualification/mode-order.json"
    )
    base = M.C.B.ROOTS["current"]
    rows, bindings = [], []
    for row in derivation["rows"]:
        case = row["case"]
        for pool in M.POOLS:
            observations = {}
            for route in ("old", "new"):
                test = M.OUT / "runtime-final" / route / case / f"pool-{pool}"
                caps = M.OUT / "companions" / route / case / f"pool-{pool}/dumps" / case
                payloads, typed = {}, {}
                for mode in ("None", "Auto"):
                    folder = caps / f"classes-{mode}-t{pool}"
                    frozen = base / case / f"classes-{mode}-t1"
                    for file in (
                        "model.json",
                        "request.json",
                        "options.json",
                        "grouping.json",
                        "payload.json",
                        "typed-owners.json",
                    ):
                        same(
                            M.load(folder / file),
                            M.load(frozen / file),
                            f"{case}/{route}/{pool}/{mode}/{file}: actual frozen-domain binding including ordered interfaces/ports",
                        )
                    meta = M.load(folder / "capture.json")
                    M.require(
                        meta["case"] == case
                        and meta["classes"] == mode
                        and meta["grouping"] == "None",
                        "metadata",
                        case,
                    )
                    M.require(
                        meta["run"]
                        == f"Run {{ threads: Some({pool}), parallel: true, check_heredity: false, cancelled: false }}",
                        "run",
                        case,
                    )
                    payload = M.load(folder / "payload.json")
                    same(
                        meta["prepared_len"],
                        len(payload["diagrams"]),
                        "prepared length",
                    )
                    same(meta["built_len"], len(payload["diagrams"]), "built length")
                    same(
                        M.load(folder / "request.json")["classes"],
                        mode,
                        "resolved request mode",
                    )
                    same(M.load(folder / "grouping.json"), "none", "resolved grouping")
                    actual = test / "dumps" / f"{case}.classes-{mode}-t{pool}.json"
                    sources = (
                        test / "dumps" / f"{case}.classes-{mode}-t{pool}.sources.json"
                    )
                    same(
                        M.load(actual),
                        payload,
                        "complete actual test -> observer payload",
                    )
                    same(M.load(sources), payload, "complete actual source companion")
                    typed[mode] = M.load(folder / "typed-owners.json")
                    payloads[mode] = payload
                    owners = {g["id"]: g["owners"] for g in typed[mode]}
                    for graph in payload["diagrams"]:
                        live = {o["owner"]: o for o in owners[graph["id"]]}
                        local = M.C.B.owners(graph["payload"])
                        same(set(live), set(local), "owner domain")
                        for name, (expression, interface) in local.items():
                            same(
                                live[name]["expression"], expression, "live/raw binding"
                            )
                            same(
                                live[name]["producer_decoded_expression"],
                                expression,
                                "live/decoded binding",
                            )
                            if interface is not None:
                                same(
                                    interface,
                                    live[name]["ordered_ports"],
                                    "actual ordered interface",
                                )
                    bindings.append(
                        {
                            "case": case,
                            "pool": pool,
                            "route": route,
                            "mode": mode,
                            "run": meta["run"],
                            "test_observer_byte_equal": actual.read_bytes()
                            == (folder / "payload.json").read_bytes(),
                            "source_byte_equal": sources.read_bytes()
                            == actual.read_bytes(),
                            "approved_payload_byte_equal": (
                                folder / "payload.json"
                            ).read_bytes()
                            == (frozen / "payload.json").read_bytes(),
                            "files": M.C.files(folder),
                            "actual_test_sha256": M.sha(actual),
                            "source_companion_sha256": M.sha(sources),
                        }
                    )
                mapping = modes(
                    payloads["None"],
                    payloads["Auto"],
                    typed["None"],
                    typed["Auto"],
                    case,
                )
                expected_mapping = next(
                    r["mapping"]
                    for r in approved_order
                    if r["route"] == "current" and r["case"] == case
                )
                same(
                    mapping,
                    expected_mapping,
                    "certified complete mode/source/output correspondence",
                )
                compact = (test / "dumps" / (case + ".census")).read_text()
                computed = D.summary(payloads["None"], compact)
                same(computed, row["current"], "full independently recomputed Row")
                observations[route] = (payloads, compact)
                if route == "old":
                    stderr = (test / "execution/stderr.txt").read_text()
                    left = re.search(
                        r"^  left: (Summary \{.*\})$", stderr, re.MULTILINE
                    )
                    right = re.search(
                        r"^ right: (Summary \{.*\})$", stderr, re.MULTILINE
                    )
                    M.require(left and right, "assertion-summary", case)
                    same(debug_row(left[1]), computed, "complete old assertion left")
                    same(
                        debug_row(right[1]),
                        row["expected"],
                        "complete old assertion right",
                    )
                    M.require(
                        [k for k in computed if computed[k] != row["expected"][k]]
                        == ["numerators"],
                        "failure-only",
                        case,
                    )
            same(
                observations["old"],
                observations["new"],
                "old/new complete census and payloads",
            )
            rows.append(
                {
                    "case": case,
                    "pool": pool,
                    "old_failed_only_numerator": True,
                    "new_full_row_pass": True,
                    "current_row": row["current"],
                    "mode_permutation": [
                        (p["none"], p["auto"])
                        for p in mapping
                        if p["none"] != p["auto"]
                    ],
                }
            )
    M.write(
        M.OUT / "validation.json",
        {
            "rows": rows,
            "capture_bindings": bindings,
            "limits": "same-source observer joined to each actual assertion ELF by complete payloads; not the assertion executable's in-memory owner instrumentation; interface presence/order and Run explicitly observed, not derived by approved mathematical validator",
        },
    )
    print(
        "validated",
        len(rows),
        "old/new pairs and",
        len(bindings),
        "complete observer-mode bindings",
    )


if __name__ == "__main__":
    main()
