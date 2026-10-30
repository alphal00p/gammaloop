"""Fixed new-domain validator: approved checker mechanisms, independently bound shape."""

import argparse
import sys
import time

import boundary as B
import common as M


def validate(mode, roots, anchors, out):
    v = B.checker()
    # Complete registered fixed shape and cross-route input/physical alignment
    # are enforced by v.validate before any low-level native comparisons.
    v.validate(M.CASES[0], mode, roots, anchors, out)
    # Capture integrity is distinct from named-axis mathematical transport.
    # Unlike prior partial anchor layer, bind actual presence AND port ORDER.
    for route in roots:
        d = roots[route] / M.CASES[0] / f"classes-{mode}-t1"
        a = anchors[route] / M.CASES[0] / f"classes-{mode}-t1"
        typed, old = M.load(d / "typed-owners.json"), M.load(a / "typed-owners.json")
        am = {g["id"]: g for g in old}
        for g in typed:
            ao = {o["owner"]: o for o in am[g["id"]]["owners"]}
            for owner in g["owners"]:
                M.require(
                    M.C.strict_equal(
                        owner["ordered_ports"], ao[owner["owner"]]["ordered_ports"]
                    ),
                    "capture-port-order",
                    route + "/" + owner["owner"],
                )
        p, ap = M.load(d / "payload.json"), M.load(a / "payload.json")
        for g, ag in zip(p["diagrams"], ap["diagrams"]):
            for obj, ref in [
                (g["payload"], ag["payload"]),
                *zip(g["payload"]["vertices"], ag["payload"]["vertices"]),
                *[
                    (e[1], ae[1])
                    for e, ae in zip(g["payload"]["edges"], ag["payload"]["edges"])
                ],
            ]:
                for key in (
                    "num_interface",
                    "numerator_prefactor_interface",
                    "projector_interface",
                ):
                    M.require(
                        (key in obj) == (key in ref)
                        and M.C.strict_equal(obj.get(key), ref.get(key)),
                        "capture-interface",
                        route + "/" + key,
                    )
        meta = M.load(d / "capture.json")
        M.require(
            meta["case"] == M.CASES[0]
            and meta["classes"] == mode
            and meta["grouping"] == "None"
            and meta["run"]
            == "Run { threads: Some(1), parallel: true, check_heredity: false, cancelled: false }",
            "run-provenance",
            route,
        )
        M.require(
            meta["prepared_len"] == meta["built_len"] == len(p["diagrams"]),
            "run-provenance",
            route + ": counts",
        )
    # Characterization only: no novel transformation is added if not explained
    # by literal token identity or the already justified even-f rotation.
    pairs = M.load(out / "cross-pairs.json")
    unexplained = [
        p for p in pairs if not p["token_tree_equal"] and not p["even_f_literal"]
    ]
    M.require(
        not unexplained,
        "new-nonliteral-obstruction",
        str([p["label"] for p in unexplained]),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=("None", "Auto"), required=True)
    ap.add_argument("--out", type=M.Path, required=True)
    for r in B.ROOTS:
        ap.add_argument("--" + r, type=M.Path, default=B.ROOTS[r])
        ap.add_argument("--anchors-" + r, type=M.Path, default=B.ROOTS[r])
    args = ap.parse_args()
    roots = {r: getattr(args, r) for r in B.ROOTS}
    anchors = {r: getattr(args, "anchors_" + r) for r in B.ROOTS}
    before = B.hashes(
        {r + "-input": root for r, root in roots.items()}
        | {r + "-anchor": root for r, root in anchors.items()}
    )
    elf_before = {r: M.sha(p) for r, p in M.C.B.ELFS.items()}
    begin = time.monotonic()
    code = 0
    try:
        validate(args.mode, roots, anchors, args.out)
        verdict = {"disposition": "accept", "stage": "complete"}
    except M.C.Rejection as error:
        code, verdict = (
            2,
            {"disposition": "reject", "stage": error.stage, "detail": str(error)},
        )
    except Exception as error:  # noqa: BLE001 — CLI boundary never counts infrastructure as semantic rejection.
        code, verdict = (
            3,
            {
                "disposition": "harness-error",
                "stage": "infrastructure",
                "exception": type(error).__name__,
                "detail": str(error),
            },
        )
    after = {p: M.sha(p) for p in before}
    elf_after = {r: M.sha(p) for r, p in M.C.B.ELFS.items()}
    M.require(before == after and elf_before == elf_after, "input-mutated", "no credit")
    M.write(
        args.out / "execution.json",
        {
            "argv": sys.argv,
            "inputs_before": before,
            "inputs_after": after,
            "elf_before": elf_before,
            "elf_after": elf_after,
            "source_hashes": M.source_hashes(),
            "reused_mechanism_hashes": M.C.source_hashes(),
            "elapsed_seconds": time.monotonic() - begin,
            "exit": code,
        },
    )
    M.write(
        args.out / "verdict.json",
        {**verdict, "exit": code, "mode": args.mode, "case": M.CASES[0]},
    )
    sys.exit(code)


if __name__ == "__main__":
    main()
