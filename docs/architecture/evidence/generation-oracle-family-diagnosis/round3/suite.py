"""Real subprocess baselines and predicate-isolating staged-capture controls."""

import argparse
import importlib.util
import sys
from pathlib import Path

import common as C

R2_SPEC = importlib.util.spec_from_file_location(
    "round2_controls", C.HERE.parent / "round2/controls.py"
)
R2 = importlib.util.module_from_spec(R2_SPEC)
R2_SPEC.loader.exec_module(R2)
PYTHON = sys.executable
VALIDATOR = C.HERE / "validate.py"
REVIEW_B = C.WORKSPACE / ".local/review-family-r2-b-20261030T010814"


def launch(
    case,
    mode,
    root,
    roots=None,
    anchors=None,
    expected="accept",
    stage="complete",
    label="baseline",
    extra_env=None,
):
    argv = [
        PYTHON,
        "-B",
        VALIDATOR,
        "--case",
        case,
        "--mode",
        mode,
        "--out",
        root / "result",
    ]
    for route in C.B.ROOTS:
        if roots:
            argv += ["--" + route, roots[route]]
        if anchors:
            argv += ["--anchors-" + route, anchors[route]]
    proc, record = C.run(argv, root / "invocation", extra_env=extra_env)
    verdict = C.load(root / "result/verdict.json")
    row = {
        "label": label,
        "case": case,
        "mode": mode,
        "expected": expected,
        "expected_stage": stage,
        "observed": verdict,
        "command": record,
        "matched": verdict["disposition"] == expected
        and verdict["stage"] == stage
        and proc.returncode
        == (0 if expected == "accept" else 2 if expected == "reject" else 3),
    }
    C.write(root / "control.json", row)
    print(case, mode, label, row["matched"], verdict["stage"], flush=True)
    return row


def stage_world(root, case, mode):
    roots = {}
    docs = {}
    for route, source in C.B.ROOTS.items():
        roots[route] = root / route
        folder = roots[route] / case / f"classes-{mode}-t1"
        folder.mkdir(parents=True, exist_ok=False)
        for file in (source / case / f"classes-{mode}-t1").iterdir():
            (folder / file.name).symlink_to(file)
        docs[route] = {
            name: C.load(folder / name)
            for name in (
                "payload.json",
                "typed-owners.json",
                "model.json",
                "request.json",
                "options.json",
                "grouping.json",
            )
        }
    return roots, docs


def save_world(roots, docs, case, mode):
    for route, files in docs.items():
        folder = roots[route] / case / f"classes-{mode}-t1"
        for name, value in files.items():
            (folder / name).unlink()
            C.write(folder / name, value)


def legacy_stage(name):
    if name in (
        "missing-owner",
        "extra-owner",
        "duplicate-owner",
        "extra-typed-domain",
    ):
        return "owner-domain"
    if name in ("owner-swap", "live-link"):
        return "payload-live-link"
    if name == "raw-link":
        return "raw-producer-binding"
    if name == "aggregate":
        return "cross-equality"
    if name == "mirrored-aggregate":
        return "aggregate-product"
    if name in ("wrong-path", "scalar-type"):
        return "schema"
    if name == "edge-tuple":
        return "edge-tuple"
    if name in ("rule", "source-factor", "routing"):
        return "cross-physical"
    if name == "model":
        return "cross-input"
    if name == "request":
        return "selected-mode"
    if name in ("group-ratio", "group-source"):
        return "group-chain"
    if name in ("port-transpose", "even-f", "metric-swap", "whole-key-reorder"):
        return "complete"
    if name in ("port-dimension", "port-variance", "port-occurrence", "hedge-tag"):
        return "cross-ports"
    if name in (
        "interface-trailing",
        "trailing-word",
        "trailing-paren",
        "trailing-semicolon",
        "truncated",
        "empty",
    ):
        return "syntax"
    if name in ("string-sign", "string-coefficient", "string-power", "odd-f"):
        return "cross-equality"
    if name in (
        "trailing-term",
        "g-attribute",
        "gamma-attribute",
        "f-attribute",
        "unmapped-old-Q",
    ):
        return "cross-definitions"
    if name == "malformed-json":
        return "json"
    raise RuntimeError("unclassified legacy control " + name)


def legacy_controls(out, case=None, mode=None):
    previous = C.load(
        C.WORKSPACE
        / "target/oracle-numerator-family-round2-candidate/controls-certified/results.json"
    )
    results = []
    for old in previous:
        if case and case != old["case"] or mode and mode != old["mode"]:
            continue
        name, c, m = old["mutation"], old["case"], old["mode"]
        root = out / c / m / name
        roots, docs = stage_world(root, c, m)
        R2.mutate(
            name,
            *(
                docs["current"][n]
                for n in (
                    "payload.json",
                    "typed-owners.json",
                    "model.json",
                    "request.json",
                )
            ),
        )
        if name == "mirrored-aggregate":
            hp, ht = (
                docs["historical"]["payload.json"],
                docs["historical"]["typed-owners.json"],
            )
            new = "2*(" + hp["diagrams"][0]["payload"]["numerator"] + ")"
            hp["diagrams"][0]["payload"]["numerator"] = new
            next(o for o in ht[0]["owners"] if o["owner"] == "numerator")[
                "producer_decoded_expression"
            ] = new
        save_world(roots, docs, c, m)
        if name == "whole-key-reorder":

            def reverse(x):
                if isinstance(x, dict):
                    return {k: reverse(v) for k, v in reversed(list(x.items()))}
                if isinstance(x, list):
                    return [reverse(v) for v in x]
                return x

            for route in roots:
                for file in docs[route]:
                    C.write(
                        roots[route] / c / f"classes-{m}-t1" / file,
                        reverse(docs[route][file]),
                    )
        if name == "malformed-json":
            path = roots["current"] / c / f"classes-{m}-t1/payload.json"
            path.write_text(path.read_text()[:-2])
        expected = old["intended"]
        row = launch(
            c, m, root, roots, expected=expected, stage=legacy_stage(name), label=name
        )
        results.append(row)
        C.write(out / "results.json", results)
    return results


# Explicit dispositions for the exact reviewer-held worlds. These expectations
# are fixed before execution; an earlier anchor check is never sufficient.
HELD_STAGES = {
    "B1m-f-antisym-dropped-current": "cross-definitions",
    "B1m-g-sym-dropped-current": "cross-definitions",
    "B1m-gamma-linear-dropped-current": "cross-definitions",
    "B1m-f-antisym-dropped-historical": "cross-definitions",
    "B1m-internal-conflict-one-vertex": "symbol-definition",
    "B1m-old-Q-constructor-in-current": "cross-definitions",
    "B1m-extra-tag-current-only": "cross-definitions",
    "B1m-unknown-tag-both-routes": "complete",
    "B1m-attr-reorder-current": "cross-definitions",
    "B1m-unparseable-spelling-projector-only": "native-parse",
    "N3m-aggregate-doubled-both-routes": "aggregate-product",
    "N3m-local-factor-negated-both-routes": "aggregate-product",
    "N3m-aggregate-even-f-both-routes": "complete",
    "N3m-aggregate-odd-f-both-routes": "aggregate-product",
    "X-odd-f-current-only": "cross-equality",
    "X-vertex-numerators-swapped-current": "cross-equality",
    "P-port-dimension-current": "cross-ports",
    "P-port-transpose-current": "complete",
    "P-port-repeated-current": "cross-ports",
    "P-interface-mismatch-current": "interface-live-link",
    "O-typed-owner-removed-current": "owner-domain",
    "O-extra-typed-owner-current": "owner-domain",
    "O-decoded-differs-from-payload-current": "payload-live-link",
    "O-current-raw-differs": "raw-producer-binding",
    "O-historical-raw-garbage": "complete",
    "O-hidden-vertex-expression-current": "schema",
    "O-hidden-typed-field-current": "schema",
    "O-hidden-payload-owner-current": "schema",
    "G-ratio-two-both-routes": "group-chain",
    "G-bool-source-index-both-routes": "schema",
    "G-bool-master-both-routes": "schema",
    "X-vertex-interaction-current": "cross-physical",
    "X-edge-particle-current": "cross-physical",
    "X-edge-incidence-current": "cross-physical",
    "X-routing-current": "cross-physical",
    "X-model-mass-current": "cross-input",
    "X-request-mode-both-routes": "selected-mode",
    "X-grouping-both-routes": "selected-grouping",
    "X-payload-overall-factor-current": "cross-physical",
    "X-ct-vertex-rule-current": "cross-physical",
    "X-extra-unit-edge-current": "cross-owner-domain",
    "X-missing-edge-current": "cross-owner-domain",
}


def reviewer_controls(out):
    old = C.load(REVIEW_B / "results/mirrored_heldouts.json")
    results = []
    for row in old:
        name, case, mode = row["name"], row["case"], row["mode"]
        world = REVIEW_B / "scratch/worlds" / name
        roots = {
            r: world
            / (
                "generation-oracle-numerator-families"
                if r == "current"
                else "generation-oracle-family-historical"
            )
            / "target/oracle-numerator-family-diagnosis"
            / ("adapted-current" if r == "current" else "adapted-historical")
            / "captures"
            for r in C.B.ROOTS
        }
        expected = "accept" if HELD_STAGES[name] == "complete" else "reject"
        result = launch(
            case, mode, out / name, roots, roots, expected, HELD_STAGES[name], name
        )
        result["round2_reviewer_observation"] = row["verdict"]
        results.append(result)
        C.write(out / "results.json", results)
        if name.startswith("X-") and name not in (
            "X-odd-f-current-only",
            "X-vertex-numerators-swapped-current",
        ):
            # Same immutable old validator/world, new output path: reproduce the
            # exact old acceptance or its real KeyError/exit-1 harness error.
            oldval = (
                world
                / "generation-oracle-numerator-families/docs/architecture/evidence/generation-oracle-family-diagnosis/round2/validate.py"
            )
            _proc, command = C.run(
                [
                    PYTHON,
                    "-I",
                    "-B",
                    oldval,
                    "--case",
                    case,
                    "--mode",
                    mode,
                    "--out",
                    out / name / "old-result",
                ],
                out / name / "old-invocation",
            )
            result["old_command"] = command
            result["old_verdict"] = C.load(out / name / "old-result/verdict.json")
            C.write(out / "results.json", results)
    return results


def anchor_controls(out):
    rows = []
    for case in C.B.CASES:
        for mode in C.B.MODES:
            root = out / case / mode
            roots, docs = stage_world(root, case, mode)
            for route in docs:
                docs[route]["model.json"]["particles"][0]["mass"] = (
                    "round3_anchor_only_mass"
                )
            save_world(roots, docs, case, mode)
            rows.append(
                launch(
                    case,
                    mode,
                    root,
                    roots,
                    expected="reject",
                    stage="anchor-input",
                    label="anchor-only-both-models",
                )
            )
            C.write(out / "results.json", rows)
    return rows


def domain_controls(out):
    """Every domain exercises semantic predicates with matching mutated anchors."""
    rows = []
    mutations = (
        "missing-owner",
        "extra-owner",
        "duplicate-owner",
        "owner-swap",
        "live-link",
        "rule",
        "model",
        "request",
        "group-ratio",
        "scalar-type",
        "edge-tuple",
        "port-dimension",
        "port-transpose",
        "interface-trailing",
        "mirrored-aggregate",
        "wrong-path",
        "whole-key-reorder",
        "projector-real",
    )
    for case in C.B.CASES:
        for mode in C.B.MODES:
            for name in mutations:
                root = out / case / mode / name
                roots, docs = stage_world(root, case, mode)
                if name not in ("whole-key-reorder", "projector-real"):
                    R2.mutate(
                        name,
                        *(
                            docs["current"][n]
                            for n in (
                                "payload.json",
                                "typed-owners.json",
                                "model.json",
                                "request.json",
                            )
                        ),
                    )
                if name == "mirrored-aggregate":
                    p, t = (
                        docs["historical"]["payload.json"],
                        docs["historical"]["typed-owners.json"],
                    )
                    value = "2*(" + p["diagrams"][0]["payload"]["numerator"] + ")"
                    p["diagrams"][0]["payload"]["numerator"] = value
                    next(o for o in t[0]["owners"] if o["owner"] == "numerator")[
                        "producer_decoded_expression"
                    ] = value
                if name == "projector-real":
                    for d in docs.values():
                        p, t = d["payload.json"], d["typed-owners.json"]
                        value = (
                            "("
                            + p["diagrams"][0]["payload"]["projector"]
                            + ")*review::{ real}::x"
                        )
                        p["diagrams"][0]["payload"]["projector"] = value
                        o = next(o for o in t[0]["owners"] if o["owner"] == "projector")
                        o["producer_decoded_expression"] = value
                        o["expression"] = value
                if name == "request":
                    for d in docs.values():
                        d["request.json"]["classes"] = (
                            "Auto" if mode == "None" else "None"
                        )
                if name == "whole-key-reorder":

                    def reverse(value):
                        if isinstance(value, dict):
                            return {
                                k: reverse(v) for k, v in reversed(list(value.items()))
                            }
                        if isinstance(value, list):
                            return [reverse(v) for v in value]
                        return value

                    docs = reverse(docs)
                save_world(roots, docs, case, mode)
                stage = (
                    "native-parse" if name == "projector-real" else legacy_stage(name)
                )
                rows.append(
                    launch(
                        case,
                        mode,
                        root,
                        roots,
                        roots,
                        "accept" if stage == "complete" else "reject",
                        stage,
                        name,
                    )
                )
                C.write(out / "results.json", rows)
    return rows


def staged_controls(out):
    rows = []
    for old in C.load(REVIEW_B / "results/staged_json_heldouts.json"):
        name, case, mode = old["name"], old["case"], old["mode"]
        roots = {r: REVIEW_B / "scratch/staged" / name / r for r in C.B.ROOTS}
        # Same rejection intents; stronger real predicates than the old
        # anchor-relative schema. Both old and smoke observations are retained.
        stage = {
            "json-nan-symmetry-factor": "json",
            "route-swap": "raw-producer-binding",
            "wrong-case-model": "cross-input",
        }.get(name, old["verdict"]["stage"])
        rows.append(
            launch(
                case,
                mode,
                out / name,
                roots,
                expected=old["expected"],
                stage=stage,
                label=name,
            )
        )
        C.write(out / "results.json", rows)
    return rows


def infrastructure_controls(out):
    case, mode = C.B.CASES[0], C.B.MODES[0]
    rows = [
        launch(
            case,
            mode,
            out / "missing-input",
            roots={r: out / "absent-input" / r for r in C.B.ROOTS},
            expected="harness-error",
            stage="infrastructure",
            label="missing-input",
        )
    ]
    empty = out / "empty-config"
    empty.mkdir(parents=True)
    rows.append(
        launch(
            case,
            mode,
            out / "missing-license",
            expected="harness-error",
            stage="infrastructure",
            label="missing-license",
            extra_env={"XDG_CONFIG_HOME": str(empty.absolute())},
        )
    )
    C.write(out / "results.json", rows)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "kind",
        choices=(
            "baseline",
            "legacy",
            "reviewer",
            "anchor",
            "domain",
            "staged",
            "infrastructure",
        ),
    )
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--case", choices=C.B.CASES)
    ap.add_argument("--mode", choices=C.B.MODES)
    args = ap.parse_args()
    if args.kind == "baseline":
        rows = []
        for case in (args.case,) if args.case else C.B.CASES:
            for mode in (args.mode,) if args.mode else C.B.MODES:
                rows.append(launch(case, mode, args.out / case / mode))
                C.write(args.out / "results.json", rows)
    elif args.kind == "legacy":
        rows = legacy_controls(args.out, args.case, args.mode)
    elif args.kind == "reviewer":
        rows = reviewer_controls(args.out)
    elif args.kind == "domain":
        rows = domain_controls(args.out)
    elif args.kind == "staged":
        rows = staged_controls(args.out)
    elif args.kind == "infrastructure":
        rows = infrastructure_controls(args.out)
    else:
        rows = anchor_controls(args.out)
    C.require(
        rows and all(r["matched"] for r in rows),
        "suite",
        "unmatched intended disposition or predicate",
    )


if __name__ == "__main__":
    main()
