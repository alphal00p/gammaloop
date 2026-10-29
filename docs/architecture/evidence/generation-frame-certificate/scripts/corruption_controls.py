"""Adversarial semantic replay with injected reads; no SHA guard is credited."""

import contextlib
import io
import json
import traceback
from pathlib import Path
from unittest.mock import patch

root = Path("target/frame-certificate")
scripts = Path("docs/architecture/evidence/generation-frame-certificate/scripts")
original_read = Path.read_text
original_write = Path.write_text
out = root / "round2-controls"
out.mkdir(exist_ok=True)


def path(name):
    return (root / name).resolve()


def document(name):
    return json.loads(original_read(path(name)))


def pair(side, case="comparison-qed", index=65, stage="comparison"):
    return f"pair-captures-{side}-graphproof/{case}/frame-pairs/source-{index}/{stage}.json"


def frame(side, case="comparison-qed", index=65):
    return (
        f"pair-captures-{side}-graphproof/{case}/frame-pairs/source-{index}/frame.txt"
    )


def add_json(changes, name, change):
    value = document(name)
    change(value)
    changes[path(name)] = json.dumps(value, ensure_ascii=False)


def add_frame(changes, name, key, change):
    lines = original_read(path(name)).splitlines()
    seen = 0
    for index, line in enumerate(lines):
        fields = line.split("\t")
        if fields[0] == key:
            fields[1] = change(fields[1])
            lines[index] = "\t".join(fields)
            seen += 1
    assert seen == 1, (name, key, seen)
    changes[path(name)] = "\n".join(lines) + "\n"


def run(script, changes, label):
    script_path = scripts / f"{script}.py"
    code = original_read(script_path)
    log = io.StringIO()

    def read(self, *args, **kwargs):
        return changes.get(self.resolve(), original_read(self, *args, **kwargs))

    def write(self, data, *args, **kwargs):
        assert self.parent == root, self
        return original_write(out / f"{label}-{self.name}", data, *args, **kwargs)

    with (
        patch.object(Path, "read_text", read),
        patch.object(Path, "write_text", write),
        contextlib.redirect_stdout(log),
    ):
        try:
            exec(compile(code, str(script_path), "exec"), {"__name__": "__main__"})  # noqa: S102
            status = {"passed": True}
        except (AssertionError, ValueError, KeyError) as error:
            status = {
                "passed": False,
                "error_type": type(error).__name__,
                "error": str(error)[:500],
                "site": traceback.extract_tb(error.__traceback__)[-1].name,
            }
    original_write(out / f"{label}.log", log.getvalue())
    return status


def mutate(case):
    changes = {}
    if case in ("model_mass", "model_denominator", "model_parameter"):
        for side in ("historical", "current"):
            for name in ("comparison-qed", "representative-qed", "representative-lr"):
                file = f"model-hooks-{side}/{name}/effective-model.json"

                def change(model):
                    if case == "model_mass":
                        next(p for p in model["particles"] if p["name"] == "d")[
                            "mass"
                        ] = "MD"
                    elif case == "model_denominator":
                        next(p for p in model["propagators"] if p["particle"] == "s")[
                            "denominator"
                        ] += "-MS^2"
                    else:
                        model["parameters"][0]["name"] += "_corrupt"

                add_json(changes, file, change)
        # The detached auxiliary is changed too: rejection must bind live graph owners.
        file = "frame-model.json"

        def change(model):
            if case == "model_mass":
                next(p for p in model["particles"] if p["name"] == "d")["mass"] = "MD"
            elif case == "model_denominator":
                next(p for p in model["propagators"] if p["particle"] == "s")[
                    "denominator"
                ] += "-MS^2"
            else:
                model["parameters"][0]["name"] += "_corrupt"

        add_json(changes, file, change)
    elif case.startswith("prefactor_"):
        if case == "prefactor_trailing_coefficient":
            change = lambda s: s + "*2"
        elif case == "prefactor_trailing_symbol":
            change = lambda s: s + "*review::z"
        elif case == "prefactor_dimension":
            change = lambda s: s.replace("::mink(4,", "::mink(5,", 1)
        elif case == "prefactor_representation":
            change = lambda s: s.replace("::mink(", "::euc(", 1)
        elif case == "prefactor_dummy":
            change = lambda s: s.replace("::vertex(0,9)", "::vertex(0,8)", 1)
        elif case == "prefactor_owner":
            change = lambda s: s.replace("::vertex(0,9)", "::vertex(1,9)", 1)
        elif case == "prefactor_effective_only":
            change = lambda s: s.replace("-1*", "2*", 1)
        else:
            raise ValueError(case)
        for side in ("historical", "current"):
            if case != "prefactor_effective_only":
                add_frame(changes, frame(side), "comparison_prefactor", change)
            add_json(
                changes,
                pair(side),
                lambda d: d.__setitem__(
                    "numerator_prefactor", change(d["numerator_prefactor"])
                ),
            )
    elif case == "routing_both_producers":
        for side in ("historical", "current"):
            add_json(
                changes,
                pair(side),
                lambda d: d["loop_momentum_basis"]["edge_signatures"]["5"][
                    "loops"
                ].__setitem__(
                    0, d["loop_momentum_basis"]["edge_signatures"]["5"]["loops"][0] + 1
                ),
            )
    elif case.startswith("final_") or case == "master_swap_both":
        sides = (
            ("historical", "current") if case == "master_swap_both" else ("current",)
        )
        for side in sides:
            name = f"{side}-frozen/comparison-qed.json"

            def change(d):
                if case == "final_edge":
                    d["diagrams"][0]["payload"]["edges"][0][1]["numerator"] = "2"
                elif case == "final_aggregate":
                    d["diagrams"][0]["payload"]["numerator"] = (
                        "2*(" + d["diagrams"][0]["payload"]["numerator"] + ")"
                    )
                elif case == "final_interface":
                    d["diagrams"][0]["payload"]["edges"][0][1]["num_interface"] = [
                        "corrupt"
                    ]
                elif case == "master_swap_both":
                    d["groups"][0]["master"], d["groups"][1]["master"] = (
                        d["groups"][1]["master"],
                        d["groups"][0]["master"],
                    )
                else:
                    raise ValueError(case)

            add_json(changes, name, change)
    elif case == "member_ratio":
        for side in ("historical", "current"):
            add_json(
                changes,
                f"{side}-frozen/comparison-qed.json",
                lambda d: next(g for g in d["groups"] if len(g["members"]) == 2)[
                    "members"
                ][1].__setitem__("ratio", "-1"),
            )
    elif case == "signed_edge_map":
        for side in ("historical", "current"):
            name = frame(side)
            lines = original_read(path(name)).splitlines()
            for i, line in enumerate(lines):
                fields = line.split("\t")
                if fields[:3] == ["edge_map", "comparison", "0"]:
                    fields[4] = str(1 - int(fields[4]))
                    lines[i] = "\t".join(fields)
                    break
            else:
                raise AssertionError("edge map absent")
            changes[path(name)] = "\n".join(lines) + "\n"
    elif case == "combined_prefactor_routing":
        for side in ("historical", "current"):
            add_json(
                changes,
                pair(side),
                lambda d: (
                    d.__setitem__(
                        "numerator_prefactor", d["numerator_prefactor"] + "*2"
                    ),
                    d["loop_momentum_basis"]["edge_signatures"]["5"][
                        "loops"
                    ].__setitem__(
                        0,
                        d["loop_momentum_basis"]["edge_signatures"]["5"]["loops"][0]
                        + 1,
                    ),
                ),
            )
    else:
        raise ValueError(case)
    return changes


controls = [
    ("model_mass", "model_binding"),
    ("model_denominator", "model_binding"),
    ("model_parameter", "model_binding"),
    ("prefactor_trailing_coefficient", "prefactor_proof"),
    ("prefactor_trailing_symbol", "prefactor_proof"),
    ("prefactor_dimension", "prefactor_proof"),
    ("prefactor_representation", "prefactor_proof"),
    ("prefactor_dummy", "prefactor_proof"),
    ("prefactor_owner", "prefactor_proof"),
    ("prefactor_effective_only", "prefactor_proof"),
    ("routing_both_producers", "routing_proof"),
    ("final_edge", "full_key_join"),
    ("final_aggregate", "full_key_join"),
    ("final_interface", "full_key_join"),
    ("master_swap_both", "full_key_join"),
    ("member_ratio", "member_proof"),
    ("signed_edge_map", "generated_map_proof"),
    ("combined_prefactor_routing", "prefactor_proof"),
    ("combined_prefactor_routing", "routing_proof"),
]
results = []
for script in sorted({script for _, script in controls}):
    status = run(script, {}, "positive-" + script)
    assert status["passed"], (script, status)
for number, (case, script) in enumerate(controls):
    changes = mutate(case)
    status = run(script, changes, f"{number:02d}-{case}-{script}")
    assert not status["passed"], (case, script, status)
    results.append(
        {
            "case": case,
            "semantic_script": script,
            "changed_reads": [str(p) for p in changes],
            **status,
        }
    )
summary = {
    "positive_replays": len({script for _, script in controls}),
    "semantic_mutation_rejections": len(results),
    "sha_integrity_guards_used": False,
    "controls": results,
}
(root / "corruption-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(
    json.dumps(
        {
            "positive_replays": summary["positive_replays"],
            "semantic_mutation_rejections": len(results),
        },
        indent=2,
    )
)
