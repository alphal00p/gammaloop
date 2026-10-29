import collections
import json
from pathlib import Path

root = Path("target/frame-certificate")
auxiliary_model = json.loads((root / "frame-model.json").read_text())
results = []
for side in ("historical", "current"):
    model = json.loads(
        (root / f"model-hooks-{side}/comparison-qed/effective-model.json").read_text()
    )
    receipt = json.loads(
        (
            root / f"model-hooks-{side}/comparison-qed/effective-bindings.json"
        ).read_text()
    )
    assert model == auxiliary_model
    particles = model["particles"]
    propagators = {record["particle"]: record for record in model["propagators"]}
    assert receipt["model_fingerprint_bytes"] == list(
        bytes.fromhex(receipt["model_fingerprint"])
    )
    for name in ("u", "u~", "d", "d~", "s", "s~"):
        assert next(p for p in particles if p["name"] == name)["mass"] == "ZERO"
        assert name in propagators
        assert all(
            mass not in propagators[name]["denominator"] for mass in ("MU", "MD", "MS")
        )
    for name in ("d", "d~", "s", "s~"):
        assert propagators[name]["denominator"] == "(P(idx(1,1)))^2"
    final = json.loads((root / f"{side}-frozen/comparison-qed.json").read_text())
    pair_root = root / f"pair-captures-{side}-graphproof/comparison-qed/frame-pairs"
    counts = collections.Counter()
    witnesses = []
    for group in final["groups"]:
        if len(group["members"]) == 1:
            continue
        assert len(group["members"]) == 2
        first, second = group["members"]
        a = json.loads(
            (pair_root / f"source-{first['source_index']}/comparison.json").read_text()
        )
        b = json.loads(
            (pair_root / f"source-{second['source_index']}/comparison.json").read_text()
        )
        assert a["model"] == b["model"] == receipt["model_fingerprint_bytes"]
        assert a["loop_momentum_basis"] == b["loop_momentum_basis"]
        assert a["cuts"] == b["cuts"]
        assert len(a["edges"]) == len(b["edges"])
        switches = []
        for i, (old, new) in enumerate(zip(a["edges"], b["edges"])):
            assert old[0] == new[0]
            old_data, new_data = old[1], new[1]
            assert {
                k: v
                for k, v in old_data.items()
                if k not in ("particle", "numerator", "num_interface")
            } == {
                k: v
                for k, v in new_data.items()
                if k not in ("particle", "numerator", "num_interface")
            }
            old_particle = particles[old_data["particle"]]
            new_particle = particles[new_data["particle"]]
            old_propagator = propagators[old_particle["name"]]
            new_propagator = propagators[new_particle["name"]]
            assert old_particle["mass"] == new_particle["mass"]
            assert old_propagator["denominator"] == new_propagator["denominator"]
            counts["edges"] += 1
            if old_data["particle"] != new_data["particle"]:
                counts["particle_switches"] += 1
                switches.append(
                    {
                        "edge": i,
                        "old": old_particle["name"],
                        "new": new_particle["name"],
                    }
                )
        counts["merged_pairs"] += 1
        witnesses.append(
            {"members": [first["name"], second["name"]], "switches": switches}
        )
    results.append({"side": side, "counts": dict(counts), "witnesses": witnesses})
(root / "member-denominator-summary.json").write_text(
    json.dumps(results, indent=2) + "\n"
)
print(json.dumps(results, indent=2))
