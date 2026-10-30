"""Four fixed domains: semantic alignment first; capture anchors separately.

Exit 0 accept, 2 checked semantic/contract rejection, 3 infrastructure/setup.
"""

import argparse
import collections
import datetime
import json
import re
import sys
import time
from pathlib import Path

import common as C

B = C.B
FILES = (
    "model.json",
    "request.json",
    "options.json",
    "grouping.json",
    "payload.json",
    "typed-owners.json",
)
MARKER = "family_certificate::{}::round3_parse_witness"
EXPECTED_ELFS = {
    "historical": "1ae86dcea5b15ff04274a91557a28dcf75d57154ee0882ed51d5c86c53be81b6",
    "current": "41c68004f03d25e4a47c58cc5d2bb75e04af39cebe525313fab67873127108ae",
}
SCHEMA = C.load(C.HERE / "schema.json")


def fixed_schema(value, path):
    C.require(
        type(value).__name__ in SCHEMA["types"].get(path, []), "schema", path + ": type"
    )
    if isinstance(value, dict):
        signatures = SCHEMA["object_keys"].get(path, [])
        numeric = "numeric-map" in signatures
        C.require(
            (numeric and all(re.fullmatch(r"0|[1-9][0-9]*", k) for k in value))
            or sorted(value) in signatures,
            "schema",
            path + ": fields",
        )
        for k, v in value.items():
            fixed_schema(v, path + "/" + ("#" if numeric else k))
    elif isinstance(value, list):
        edge_tuple = path.endswith("/payload/edges/*")
        if edge_tuple:
            C.require(len(value) == 2, "edge-tuple", path)
        for i, v in enumerate(value):
            fixed_schema(v, path + "/" + (str(i) if edge_tuple else "*"))


def parse(text, mapped=False):
    try:
        return B.parse(text, mapped)
    except B.Reject as error:
        raise C.Rejection(error.stage, str(error)) from error


def mapped_text(text, historical):
    parse(text)
    return B.TOKEN.sub(
        lambda m: B.SPELLINGS.get(m.group(), m.group()) if historical else m.group(),
        text,
    )


def eq(a, b, stage, detail):
    C.require(C.strict_equal(a, b), stage, detail)


def owner_table(p, typed, route):
    ids = [g["id"] for g in p["diagrams"]]
    C.require(
        len(ids) == len(set(ids)), "graph-domain", route + ": duplicate payload ID"
    )
    tm = {g["id"]: g for g in typed}
    C.require(
        len(tm) == len(typed) == len(ids) and set(tm) == set(ids),
        "owner-domain",
        route + ": graph IDs",
    )
    tables = {}
    for graph in p["diagrams"]:
        g = tm[graph["id"]]
        po = B.owners(graph["payload"])
        om = {o["owner"]: o for o in g["owners"]}
        C.require(
            len(om) == len(g["owners"]) and set(om) == set(po),
            "owner-domain",
            route + "/" + graph["id"],
        )
        for name, (expression, interface) in po.items():
            owner = om[name]
            eq(
                expression,
                owner["producer_decoded_expression"],
                "payload-live-link",
                route + "/" + graph["id"] + "/" + name,
            )
            if route == "current":
                eq(
                    expression,
                    owner["expression"],
                    "raw-producer-binding",
                    route + "/" + name,
                )
            if interface is not None:
                eq(
                    interface,
                    owner["ordered_ports"],
                    "interface-live-link",
                    route + "/" + name,
                )
        tables[graph["id"]] = om
    return tables


def physical_alignment(h, c, ht, ct):
    C.require(
        len(h["diagrams"]) == len(c["diagrams"]), "cross-graph-domain", "diagram count"
    )
    keys = lambda p: [
        json.dumps(g["key"], sort_keys=True, separators=(",", ":"))
        for g in p["diagrams"]
    ]
    hk, ck = keys(h), keys(c)
    C.require(
        collections.Counter(hk) == collections.Counter(ck) and len(set(hk)) == len(hk),
        "cross-graph-domain",
        "full canonical keys with multiplicity",
    )
    for gi, (a, b) in enumerate(zip(h["diagrams"], c["diagrams"])):
        for field in set(a) - {"payload"}:
            eq(a[field], b[field], "cross-physical", f"diagram/{gi}/{field}")
        x, y = a["payload"], b["payload"]
        # Only these exact expression/interface fields are transformed below.
        transformed = {
            "vertices",
            "edges",
            "numerator",
            "numerator_prefactor",
            "projector",
            "num_interface",
            "numerator_prefactor_interface",
            "projector_interface",
        }
        for field in (set(x) | set(y)) - transformed:
            eq(x.get(field), y.get(field), "cross-physical", f"payload/{gi}/{field}")
        C.require(
            len(x["vertices"]) == len(y["vertices"])
            and len(x["edges"]) == len(y["edges"]),
            "cross-owner-domain",
            f"graph/{gi}/positional counts",
        )
        C.require(
            set(ht[a["id"]]) == set(ct[b["id"]]),
            "cross-owner-domain",
            f"graph/{gi}/owner paths",
        )
        for vi, (v, w) in enumerate(zip(x["vertices"], y["vertices"])):
            for field in (set(v) | set(w)) - {"numerator", "num_interface"}:
                eq(
                    v.get(field),
                    w.get(field),
                    "cross-physical",
                    f"vertex/{gi}/{vi}/{field}",
                )
        for ei, (v, w) in enumerate(zip(x["edges"], y["edges"])):
            eq(v[0], w[0], "cross-physical", f"edge/{gi}/{ei}/incidence")
            for field in (set(v[1]) | set(w[1])) - {"numerator", "num_interface"}:
                eq(
                    v[1].get(field),
                    w[1].get(field),
                    "cross-physical",
                    f"edge/{gi}/{ei}/{field}",
                )
        for owner in ht[a["id"]]:
            for flag in ("is_metric", "is_composite"):
                eq(
                    ht[a["id"]][owner][flag],
                    ct[b["id"]][owner][flag],
                    "cross-flags",
                    f"{gi}/{owner}/{flag}",
                )


def group_chain(p, route):
    C.require(len(p["groups"]) == len(p["diagrams"]), "group-chain", route + ": count")
    for i, (g, graph) in enumerate(zip(p["groups"], p["diagrams"])):
        expected = {
            "master": i,
            "members": [
                {
                    "id": graph["id"],
                    "name": graph["name"],
                    "overall_factor": graph["overall_factor"],
                    "ratio": "1",
                    "source_index": i,
                }
            ],
        }
        eq(g, expected, "group-chain", route + "/" + graph["id"])


def pair(a, b, label, historical=False):
    left, right = parse(a, historical), parse(b)
    try:
        C.require(
            B.definitions(left) == B.definitions(right), "symbol-definition", label
        )
    except B.Reject as error:
        raise C.Rejection("symbol-definition", str(error)) from error
    return {
        "historical": B.render(left),
        "current": B.render(right),
        "label": label,
        "token_tree_equal": left == right,
        "even_f_literal": B.cyclic(left) == B.cyclic(right),
    }


def native(rows, root, stage, both_orders=False):
    """Real parser/equality via retained ELF; explicit failure classification."""
    root.mkdir(parents=True, exist_ok=True)
    if not rows:
        C.write(root / "commands.json", [])
        return []
    # Marker avoids the old helper's post-parse 1^2=1 safeguard. A+z=B+z
    # is equivalent to A=B by exact additive cancellation, no distribution.
    for row in rows:
        C.require(
            "family_certificate::round3_parse_witness"
            not in B.definitions(parse(row["historical"]))
            and "family_certificate::round3_parse_witness"
            not in B.definitions(parse(row["current"])),
            "syntax",
            "reserved proof marker",
        )
    commands = []
    for reverse in (False, True) if both_orders else (False,):
        wrapped = [
            {
                "historical": "("
                + p["current" if reverse else "historical"]
                + ")+"
                + MARKER,
                "current": "("
                + p["historical" if reverse else "current"]
                + ")+"
                + MARKER,
            }
            for p in rows
        ]
        path = root / ("reverse-pairs.json" if reverse else "forward-pairs.json")
        C.write(path, wrapped)
        failures = []
        for kernel, elf in B.ELFS.items():
            digest = C.sha(elf)
            if digest != EXPECTED_ELFS[kernel]:
                raise RuntimeError("unqualified ELF " + str(elf))
            dest = root / (kernel + ("-reverse" if reverse else "-forward"))
            result = dest / "result.json"
            pair_digest = C.sha(path)
            proc, record = C.run(
                [elf, path.absolute(), result.absolute()],
                dest,
                runtime=True,
                bindings={
                    "pair_file_sha256": pair_digest,
                    "elf_sha256": digest,
                    "stage": stage,
                    "selected_pairs": len(rows),
                },
            )
            if C.sha(path) != pair_digest or C.sha(elf) != digest:
                raise RuntimeError("native input or ELF changed during execution")
            record.update(
                stage=stage,
                pair_file_sha256=C.sha(path),
                elf_sha256=digest,
                selected_pairs=len(rows),
                labels=[r["label"] for r in rows],
            )
            text = proc.stderr + proc.stdout
            if proc.returncode == 0:
                data = C.load(result)
                C.require(
                    len(data) == len(rows) and all(all(r.values()) for r in data),
                    "checker-result",
                    stage,
                )
                record.update(reason="equal", result_sha256=C.sha(result))
            elif (
                proc.returncode == 101
                and "assertion failed: row.parsed_equal && row.canonical_equal" in text
            ):
                record.update(reason="inequality")
                failures.append(("equality", kernel))
            elif (
                proc.returncode == 101
                and re.search(r"oracle_family_exact.rs:(35|36):", text)
                and "Err" in text
            ):
                record.update(reason="native-parse-error")
                failures.append(("parse", kernel))
            else:
                record.update(reason="infrastructure-or-unclassified-crash")
                C.write(dest / "command.json", record)
                raise RuntimeError(
                    f"unclassified {kernel} ELF failure {proc.returncode}; {dest}"
                )
            C.write(dest / "command.json", record)
            commands.append(record)
        C.write(root / "commands.json", commands)
        if failures:
            actual_stage = (
                "native-parse" if any(x[0] == "parse" for x in failures) else stage
            )
            raise C.Rejection(actual_stage, f"{stage}: {failures}; {root}")
    return commands


def validate(case, mode, roots, anchors, output):
    documents, anchor_docs, tables = {}, {}, {}
    for route in B.ROOTS:
        folder, anchor = (
            roots[route] / case / f"classes-{mode}-t1",
            anchors[route] / case / f"classes-{mode}-t1",
        )
        docs, ads = {}, {}
        for file in FILES:
            docs[file], ads[file] = C.load(folder / file), C.load(anchor / file)
            fixed_schema(docs[file], file)
        C.require(docs["request.json"]["classes"] == mode, "selected-mode", route)
        C.require(docs["grouping.json"] == "none", "selected-grouping", route)
        documents[route], anchor_docs[route] = docs, ads
        tables[route] = owner_table(
            docs["payload.json"], docs["typed-owners.json"], route
        )
        group_chain(docs["payload.json"], route)
    h, c = documents["historical"], documents["current"]
    for file in ("model.json", "request.json", "options.json", "grouping.json"):
        eq(h[file], c[file], "cross-input", file)
    physical_alignment(
        h["payload.json"], c["payload.json"], tables["historical"], tables["current"]
    )
    eq(
        h["payload.json"]["groups"],
        c["payload.json"]["groups"],
        "cross-groups",
        "ordered source/group records",
    )

    definitions, parsed_texts, cross, aggregate, anchor_pairs, ports = (
        {},
        {},
        [],
        [],
        [],
        [],
    )
    for route, docs in documents.items():
        defs = {}
        for graph in docs["payload.json"]["diagrams"]:
            for name, owner in tables[route][graph["id"]].items():
                raw = mapped_text(
                    owner["producer_decoded_expression"], route == "historical"
                )
                tree = parse(raw)
                try:
                    B.definitions(tree, defs)
                except B.Reject as error:
                    raise C.Rejection("symbol-definition", str(error)) from error
                parsed_texts.setdefault(raw, []).append([route, graph["id"], name])
                for port in owner["ordered_ports"]:
                    mapped = mapped_text(port, route == "historical")
                    try:
                        B.definitions(parse(mapped), defs)
                    except B.Reject as error:
                        raise C.Rejection("symbol-definition", str(error)) from error
                    parsed_texts.setdefault(mapped, []).append(
                        [route, graph["id"], name, "port"]
                    )
        definitions[route] = defs
    pending_cross = []
    for hg, cg in zip(h["payload.json"]["diagrams"], c["payload.json"]["diagrams"]):
        ht, ct = tables["historical"][hg["id"]], tables["current"][cg["id"]]
        C.require(set(ht) == set(ct), "cross-owner-domain", hg["id"])
        for name in ht:
            hp = [B.render(parse(x, True)) for x in ht[name]["ordered_ports"]]
            cp = [B.render(parse(x)) for x in ct[name]["ordered_ports"]]
            try:
                perm = B.permutation(hp, cp)
            except B.Reject as error:
                raise C.Rejection("cross-ports", str(error)) from error
            ports.append(
                {
                    "id": hg["id"],
                    "owner": name,
                    "historical_ports": hp,
                    "current_ports": cp,
                    "permutation": perm,
                }
            )
            pending_cross.append(
                (
                    ht[name]["producer_decoded_expression"],
                    ct[name]["producer_decoded_expression"],
                    "cross/" + hg["id"] + "/" + name,
                )
            )
    eq(
        definitions["historical"],
        definitions["current"],
        "cross-definitions",
        "complete namespace/name definitions",
    )
    cross = [pair(a, b, label, True) for a, b, label in pending_cross]
    C.write(output / "cross-ports.json", ports)
    C.write(output / "definitions.json", definitions)
    faithfulness = [
        {
            "historical": text,
            "current": B.render(parse(text)),
            "label": "faithful/" + str(i),
            "owners": refs,
        }
        for i, (text, refs) in enumerate(parsed_texts.items())
    ]
    C.write(output / "faithfulness.json", faithfulness)
    # Every distinct text reaches each native parser, including literal owners.
    native(faithfulness, output / "native-faithfulness", "faithfulness")
    C.write(output / "cross-pairs.json", cross)
    native(
        [p for p in cross if not p["token_tree_equal"]],
        output / "native-cross",
        "cross-equality",
        True,
    )
    for route, docs in documents.items():
        for graph in docs["payload.json"]["diagrams"]:
            p = graph["payload"]
            factors = [
                p["numerator_prefactor"],
                *[v["numerator"] for v in p["vertices"]],
                *[e[1]["numerator"] for e in p["edges"]],
            ]
            product = "*".join("(" + x + ")" for x in factors)
            aggregate.append(
                pair(
                    mapped_text(product, route == "historical"),
                    mapped_text(p["numerator"], route == "historical"),
                    "aggregate-product/" + route + "/" + graph["id"],
                )
            )
    C.write(output / "aggregate-pairs.json", aggregate)
    native(aggregate, output / "native-aggregate", "aggregate-product", True)

    # Additional provenance guards. Baseline anchor comparisons are explicitly
    # self-comparisons, not independent algebra or cross-route evidence.
    for route, docs in documents.items():
        ads = anchor_docs[route]
        for file in ("model.json", "request.json", "options.json", "grouping.json"):
            eq(docs[file], ads[file], "anchor-input", route + "/" + file)
        p, ap = docs["payload.json"], ads["payload.json"]
        at = owner_table(ap, ads["typed-owners.json"], route)
        try:
            physical_alignment(p, ap, tables[route], at)
        except C.Rejection as error:
            raise C.Rejection("anchor-physical", str(error)) from error
        eq(p["groups"], ap["groups"], "anchor-groups", route)
        am = {g["id"]: g for g in ads["typed-owners.json"]}
        eq(
            [g["id"] for g in p["diagrams"]],
            [g["id"] for g in ap["diagrams"]],
            "anchor-order",
            route,
        )
        for g, ag in zip(p["diagrams"], ap["diagrams"]):
            for field in set(g) - {"payload"}:
                eq(g[field], ag[field], "anchor-physical", field)
            amo = {o["owner"]: o for o in am[g["id"]]["owners"]}
            C.require(
                set(amo) == set(tables[route][g["id"]]), "anchor-owner-domain", route
            )
            for name, o in tables[route][g["id"]].items():
                ao = amo[name]
                if route == "historical":
                    eq(o["expression"], ao["expression"], "anchor-raw-cook", name)
                anchor_pairs.append(
                    pair(
                        mapped_text(
                            ao["producer_decoded_expression"], route == "historical"
                        ),
                        mapped_text(
                            o["producer_decoded_expression"], route == "historical"
                        ),
                        "anchor/" + route + "/" + g["id"] + "/" + name,
                    )
                )
    C.write(output / "anchor-pairs.json", anchor_pairs)
    native(
        [p for p in anchor_pairs if not p["token_tree_equal"]],
        output / "native-anchor",
        "anchor-equality",
        True,
    )
    ct = {}
    if case == "ctorder_qcd_uugg_L1":
        for route, docs in documents.items():
            try:
                ct[route] = B.ct_provenance(
                    docs["model.json"], docs["payload.json"], docs["request.json"]
                )
            except B.Reject as error:
                raise C.Rejection(error.stage, str(error)) from error
    C.write(
        output / "summary.json",
        {
            "case": case,
            "mode": mode,
            "graphs": len(h["payload.json"]["diagrams"]),
            "owners": len(cross),
            "nontrivial_cross_permutations": sum(
                p["permutation"] != list(range(len(p["permutation"]))) for p in ports
            ),
            "distinct_native_texts": len(faithfulness),
            "anchor_self_comparisons": sum(p["token_tree_equal"] for p in anchor_pairs),
            "nonliteral_cross": sum(not p["token_tree_equal"] for p in cross),
            "ct": ct,
        },
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", choices=B.CASES, required=True)
    ap.add_argument("--mode", choices=B.MODES, required=True)
    ap.add_argument("--historical", type=Path, default=B.ROOTS["historical"])
    ap.add_argument("--current", type=Path, default=B.ROOTS["current"])
    ap.add_argument("--anchors-historical", type=Path, default=B.ROOTS["historical"])
    ap.add_argument("--anchors-current", type=Path, default=B.ROOTS["current"])
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    roots = {"historical": args.historical, "current": args.current}
    anchors = {"historical": args.anchors_historical, "current": args.anchors_current}
    try:
        before = {
            str(p): C.sha(p)
            for root in [*roots.values(), *anchors.values()]
            for p in (root / args.case / f"classes-{args.mode}-t1").iterdir()
            if p.is_file()
        }
    except Exception as error:  # noqa: BLE001 — CLI boundary: any setup fault is exit 3, never proof credit.
        C.write(
            args.out / "verdict.json",
            {
                "disposition": "harness-error",
                "stage": "infrastructure",
                "exception": type(error).__name__,
                "detail": str(error),
                "exit": 3,
            },
        )
        sys.exit(3)
    try:
        elf_before = {k: C.sha(p) for k, p in B.ELFS.items()}
    except Exception as error:  # noqa: BLE001 — CLI boundary: retain infrastructure faults without semantic credit.
        C.write(
            args.out / "verdict.json",
            {
                "disposition": "harness-error",
                "stage": "infrastructure",
                "exception": type(error).__name__,
                "detail": str(error),
                "exit": 3,
            },
        )
        sys.exit(3)
    feature_evidence = (
        C.WORKSPACE
        / ".local/review-family-r2-a-20261030T010834Z/results/a01_fingerprints.json"
    )
    execution = {
        "argv": sys.argv,
        "started_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "source_hashes": C.source_hashes(),
        "inputs_before": before,
        "elf_hashes_before": elf_before,
        "feature_evidence": {
            "path": str(feature_evidence),
            "sha256": C.sha(feature_evidence),
            "attribution": "independent scope A, not original build logs",
        },
        "source_revisions": {
            "historical": "2b1ef4f07a44bc061bcc97c378bbdd9c41cb0e94",
            "current": "b9f788922341e7bde2e801a054423ae70ed7e9ae",
        },
        "compiler_features": "A a01_fingerprints.json: rustc 48a229cea; root default/native; 336 historical and 341 current dependency units, per-unit profiles qualified within each route; cache corroboration, not original build logs",
    }
    C.write(args.out / "execution.json", execution)
    tick, code = time.monotonic(), 0
    try:
        validate(args.case, args.mode, roots, anchors, args.out)
        verdict = {"disposition": "accept", "stage": "complete"}
    except C.Rejection as error:
        code, verdict = (
            2,
            {"disposition": "reject", "stage": error.stage, "detail": str(error)},
        )
    except Exception as error:  # noqa: BLE001 — CLI boundary: unexpected exceptions are harness failures.
        code, verdict = (
            3,
            {
                "disposition": "harness-error",
                "stage": "infrastructure",
                "exception": type(error).__name__,
                "detail": str(error),
            },
        )
    try:
        execution.update(
            inputs_after={p: C.sha(p) for p in before},
            elf_hashes_after={k: C.sha(p) for k, p in B.ELFS.items()},
        )
    except Exception as error:  # noqa: BLE001 — CLI boundary: failed final rehash invalidates the execution.
        execution["rehash_error"] = str(error)
    if execution.get("inputs_after") != before or execution[
        "elf_hashes_before"
    ] != execution.get("elf_hashes_after"):
        code, verdict = (
            3,
            {"disposition": "harness-error", "stage": "input-mutated-during-run"},
        )
    execution.update(
        elapsed_seconds=time.monotonic() - tick,
        ended_utc=datetime.datetime.now(datetime.UTC).isoformat(),
        exit=code,
    )
    C.write(args.out / "execution.json", execution)
    C.write(
        args.out / "verdict.json",
        {**verdict, "case": args.case, "mode": args.mode, "exit": code},
    )
    sys.exit(code)


if __name__ == "__main__":
    main()
