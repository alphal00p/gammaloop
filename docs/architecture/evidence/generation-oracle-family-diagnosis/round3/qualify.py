"""Derive exact None/Auto correspondence without normalizing unmatched fields."""

import collections
import json

import common as C


def main():
    rows = []
    for route, root in C.B.ROOTS.items():
        for case in C.B.CASES:
            folders = {m: root / case / f"classes-{m}-t1" for m in C.B.MODES}
            before = {
                str(p): C.sha(p)
                for d in folders.values()
                for p in d.iterdir()
                if p.is_file()
            }
            payloads = {m: C.load(folders[m] / "payload.json") for m in C.B.MODES}
            typed = {m: C.load(folders[m] / "typed-owners.json") for m in C.B.MODES}

            def key(g):
                return json.dumps(g["key"], sort_keys=True, separators=(",", ":"))

            a, b = payloads["None"], payloads["Auto"]
            ak, bk = [key(g) for g in a["diagrams"]], [key(g) for g in b["diagrams"]]
            C.require(
                collections.Counter(ak) == collections.Counter(bk)
                and len(set(bk)) == len(bk),
                "mode-domain",
                case,
            )
            mapping = []
            for i, g in enumerate(a["diagrams"]):
                j = bk.index(key(g))
                other = b["diagrams"][j]
                C.require(
                    C.strict_equal(
                        {k: v for k, v in g.items() if k != "name"},
                        {k: v for k, v in other.items() if k != "name"},
                    ),
                    "mode-payload",
                    g["id"],
                )
                groups = []
                for p, graph, index in ((a, g, i), (b, other, j)):
                    expected = {
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
                    }
                    C.require(
                        C.strict_equal(p["groups"][index], expected),
                        "mode-group",
                        g["id"],
                    )
                    groups.append(expected)
                mapping.append(
                    {
                        "none": i,
                        "auto": j,
                        "id": g["id"],
                        "key": g["key"],
                        "none_name": g["name"],
                        "auto_name": other["name"],
                        "groups": groups,
                    }
                )
            C.require(
                C.strict_equal(
                    {g["id"]: g["owners"] for g in typed["None"]},
                    {g["id"]: g["owners"] for g in typed["Auto"]},
                ),
                "mode-typed",
                case,
            )
            C.require(before == {p: C.sha(p) for p in before}, "input-mutated", case)
            rows.append(
                {
                    "case": case,
                    "route": route,
                    "mapping": mapping,
                    "moved": [p for p in mapping if p["none"] != p["auto"]],
                    "inputs_before_after": before,
                }
            )
    C.write(C.OUT / "qualification/mode-order.json", rows)


if __name__ == "__main__":
    main()
