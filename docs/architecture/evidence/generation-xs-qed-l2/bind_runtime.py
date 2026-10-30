"""Actual assertion/observer/Row binding before symbolic proof; no digest shortcut."""

import collections
import importlib.util
import json
import re
from fractions import Fraction

import boundary as B
import common as M

SPEC = importlib.util.spec_from_file_location(
    "audited_digest_encoding",
    M.HERE.parent / "generation-oracle-family-migration/derive.py",
)
D = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(D)


def weight(text):
    def evaluate(tree):
        kind = tree[0]
        if kind == "integer":
            return Fraction(int(tree[1]))
        if kind == "neg":
            return -evaluate(tree[1])
        if kind == "call":
            head = tree[1].split("::")[-1]
            M.require(
                tree[1].startswith("feynkit_generator_factor::")
                and head
                in (
                    "AutG",
                    "NumeratorIndependentSymmetryGrouping",
                    "InternalFermionLoopSign",
                    "ExternalFermionOrderingSign",
                    "AntiFermionSpinSumSign",
                )
                and len(tree[2]) == 1,
                "bookkeeping-factor",
                text,
            )
            value = evaluate(tree[2][0])
            M.require(
                value.denominator == 1
                and (
                    value > 0
                    if head in ("AutG", "NumeratorIndependentSymmetryGrouping")
                    else value in (-1, 1)
                ),
                "bookkeeping-factor",
                text,
            )
            return value
        M.require(kind in ("*", "/", "^"), "bookkeeping-factor", text)
        a, b = evaluate(tree[1]), evaluate(tree[2])
        if kind == "*":
            return a * b
        if kind == "/":
            return a / b
        M.require(b.denominator == 1, "bookkeeping-factor", "integer power")
        return a ** int(b)

    return evaluate(M.C.B.parse(text))


def main():
    D.weight = weight
    summaries, records = {}, []
    before = B.hashes()
    for route, root in B.ROOTS.items():
        runtime = M.OUT / "runtime" / route
        payloads, typed = {}, {}
        for mode in ("None", "Auto"):
            folder = root / M.CASES[0] / f"classes-{mode}-t1"
            test = runtime / "dumps" / f"{M.CASES[0]}.classes-{mode}-t1.json"
            sources = runtime / "dumps" / f"{M.CASES[0]}.classes-{mode}-t1.sources.json"
            p = M.load(folder / "payload.json")
            M.require(
                M.C.strict_equal(M.load(test), p)
                and M.C.strict_equal(M.load(sources), p),
                "assertion-observer-binding",
                route + mode,
            )
            t = M.load(folder / "typed-owners.json")
            payloads[mode], typed[mode] = p, t
            meta = M.load(folder / "capture.json")
            M.require(
                meta["case"] == M.CASES[0]
                and meta["classes"] == mode
                and meta["run"]
                == "Run { threads: Some(1), parallel: true, check_heredity: false, cancelled: false }"
                and meta["grouping"] == "None",
                "Run-provenance",
                route,
            )
            records.append(
                {
                    "route": route,
                    "mode": mode,
                    "graphs": len(p["diagrams"]),
                    "owners": sum(len(g["owners"]) for g in t),
                    "axes": sum(
                        len(o["ordered_ports"]) for g in t for o in g["owners"]
                    ),
                    "vertex_count_histogram": dict(
                        collections.Counter(
                            len(g["payload"]["vertices"]) for g in p["diagrams"]
                        )
                    ),
                    "edge_count_histogram": dict(
                        collections.Counter(
                            len(g["payload"]["edges"]) for g in p["diagrams"]
                        )
                    ),
                    "test_payload_sha256": M.sha(test),
                    "source_companion_sha256": M.sha(sources),
                    "observer_sha256": M.sha(folder / "payload.json"),
                    "byte_identity": test.read_bytes()
                    == (folder / "payload.json").read_bytes()
                    == sources.read_bytes(),
                    "Run": meta,
                    "capture_bytes": sum(f.stat().st_size for f in folder.iterdir()),
                }
            )
        # Full-key join and complete singleton source-chain; names are not
        # globally stripped, each occurrence retains positional output identity.
        key = lambda g: json.dumps(g["key"], sort_keys=True, separators=(",", ":"))
        a, b = payloads["None"], payloads["Auto"]
        bk = [key(g) for g in b["diagrams"]]
        ak = [key(g) for g in a["diagrams"]]
        M.require(
            collections.Counter(ak) == collections.Counter(bk)
            and len(set(ak)) == len(ak),
            "mode-domain",
            route,
        )
        mapping = []
        for i, g in enumerate(a["diagrams"]):
            j = bk.index(key(g))
            ag = b["diagrams"][j]
            M.require(
                M.C.strict_equal(
                    {k: v for k, v in g.items() if k != "name"},
                    {k: v for k, v in ag.items() if k != "name"},
                ),
                "mode-payload",
                route,
            )
            for p, graph, index in ((a, g, i), (b, ag, j)):
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
                M.require(
                    M.C.strict_equal(p["groups"][index], expected),
                    "mode-source-chain",
                    route,
                )
            mapping.append(
                {
                    "none": i,
                    "auto": j,
                    "id": g["id"],
                    "key": g["key"],
                    "names": [g["name"], ag["name"]],
                    "groups": [a["groups"][i], b["groups"][j]],
                }
            )
        M.require(
            M.C.strict_equal(
                {g["id"]: g["owners"] for g in typed["None"]},
                {g["id"]: g["owners"] for g in typed["Auto"]},
            ),
            "mode-typed-owners",
            route,
        )
        M.write(
            M.OUT / "binding" / (route + "-modes.json"),
            {"mapping": mapping, "byte_equal": json.dumps(a) == json.dumps(b)},
        )
        compact = (runtime / "dumps" / (M.CASES[0] + ".census")).read_text()
        summaries[route] = D.summary(a, compact)
        source = (M.OUT / "source" / route / M.MAIN).read_text()
        expected = D.row(source, M.CASES[0])
        differences = [k for k in expected if expected[k] != summaries[route][k]]
        M.require(
            differences == ([] if route == "historical" else ["numerators"]),
            "Row",
            route + str(differences),
        )
        if route == "current":
            stderr = (runtime / "generation_oracle/stderr.txt").read_text()
            for side, target in (("  left", summaries[route]), (" right", expected)):
                text = re.search(
                    "^" + side + r": (Summary \{.*\})$", stderr, re.MULTILINE
                )[1]
                text = (
                    re.sub(
                        r"Fraction \{ numerator: (-?\d+), denominator: (\d+) \}",
                        r"(\1, \2)",
                        text,
                    )
                    .replace("symmetry_factors: [", "symmetry_factors: &[")
                    .replace("Summary {", "Row {", 1)
                )
                M.require(
                    D.row("observed: " + text, "observed") == target,
                    "complete-assertion-summary",
                    side,
                )
    M.require(
        {k: v for k, v in summaries["historical"].items() if k != "numerators"}
        == {k: v for k, v in summaries["current"].items() if k != "numerators"},
        "non-numerator-Row",
        "all fields",
    )
    M.require(
        before == {p: M.sha(p) for p in before}, "input-change", "runtime binding"
    )
    M.write(
        M.OUT / "binding/runtime.json",
        {
            "summaries": summaries,
            "domains": records,
            "inputs_before_after": before,
            "limit": "observer metadata joins actual assertion producers through complete payload, not in-memory assertion instrumentation; exact rational bookkeeping heads per immutable evaluate_overall_factor, no physical parameter substitution",
        },
    )
    print(records)


if __name__ == "__main__":
    main()
