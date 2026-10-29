import hashlib
import json
import os
import re
import subprocess
import tempfile
from pathlib import Path

root = Path("target/frame-certificate")
cases = ("comparison-qed", "representative-qed", "representative-lr")
probes = {
    "current": (
        root / "binaries/current-fingerprint-probe",
        "frame_certificate::model_fingerprint_from_file",
    ),
    "historical": (
        root / "binaries/historical-fingerprint-probe",
        "frames::certificate_model_fingerprint_from_file",
    ),
}
light = {"u", "u~", "d", "d~", "s", "s~"}
frozen_effective_fingerprint = (
    "9d826a6ef2be57828b443ffcee31307e75f8c51d3de7c87b17817a21b5914876"
)


def fingerprint(side, source):
    binary, test = probes[side]
    env = os.environ.copy()
    if not env.get("SYMBOLICA_LICENSE"):
        configuration = Path(env.get("XDG_CONFIG_HOME", str(Path.home() / ".config")))
        env["SYMBOLICA_LICENSE"] = (
            (configuration / "symbolica/license").read_text().strip()
        )
    with tempfile.NamedTemporaryFile(prefix="frame-model-", suffix=".json") as file:
        file.write(source.encode())
        file.flush()
        env["FEYNKIT_FRAME_MODEL_INPUT"] = file.name
        result = subprocess.run(
            [
                str(binary),
                "--test-threads",
                "1",
                "--nocapture",
                "--ignored",
                "--exact",
                test,
            ],
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )
    assert result.returncode == 0, (side, result.stdout, result.stderr)
    found = re.findall(r"FRAME_MODEL_FINGERPRINT=([0-9a-f]{64})", result.stdout)
    assert len(found) == 1, (side, result.stdout)
    return found[0]


def validate_light_model(model):
    particles = {particle["name"]: particle for particle in model["particles"]}
    propagators = {
        propagator["particle"]: propagator for propagator in model["propagators"]
    }
    assert light <= particles.keys() and light <= propagators.keys()
    for name in light:
        assert particles[name]["mass"] == "ZERO", name
        for field in ("numerator", "denominator"):
            assert all(
                mass not in propagators[name][field] for mass in ("MU", "MD", "MS")
            )
    for name in ("d", "d~", "s", "s~"):
        assert propagators[name]["denominator"] == "(P(idx(1,1)))^2"
    assert propagators["d"]["numerator"] == propagators["s"]["numerator"]
    assert propagators["d~"]["numerator"] == propagators["s~"]["numerator"]
    return len(model["parameters"]), len(model["propagators"])


auxiliary = (root / "frame-model.json").read_text()
rows = []
model_bytes = []
for side in ("historical", "current"):
    for case in cases:
        directory = root / f"model-hooks-{side}" / case
        source = (directory / "effective-model.json").read_text()
        receipt = json.loads((directory / "effective-bindings.json").read_text())
        assert source == auxiliary
        model_bytes.append(source.encode())
        model = json.loads(source)
        parameters, propagators = validate_light_model(model)
        computed = fingerprint(side, source)
        expected = receipt["model_fingerprint"]
        assert (
            computed
            == expected
            == receipt["roundtrip_fingerprint"]
            == frozen_effective_fingerprint
        )
        expected_bytes = list(bytes.fromhex(computed))
        assert receipt["model_fingerprint_bytes"] == expected_bytes
        first = receipt["first_build"]
        second = receipt["second_build"]
        grouped = receipt["grouped"]
        names = {record["name"] for record in first}
        assert len(names) == len(first) == len(second)
        assert names == {record["name"] for record in second}
        assert names == set(receipt["named_loop_momentum_bases"])
        assert receipt["by_id_loop_momentum_bases_count"] == 0
        for record in first + second + grouped:
            assert record["model_fingerprint"] == computed
        for record in first:
            assert (
                record["loop_edges"]
                == receipt["named_loop_momentum_bases"][record["name"]]
            )
        for record in second:
            assert (
                record["loop_edges"]
                == receipt["named_loop_momentum_bases"][record["name"]]
            )
        pair_root = root / f"pair-captures-{side}-graphproof" / case / "frame-pairs"
        pair_names = {folder.name for folder in pair_root.iterdir() if folder.is_dir()}
        assert len(pair_names) == len(second)
        checked_pairs = 0
        for folder in pair_root.iterdir():
            if not folder.is_dir():
                continue
            header = dict(
                line.split("\t", 1)
                for line in (folder / "frame.txt").read_text().splitlines()
                if line.startswith(
                    ("source_index\t", "source_name\t", "requested_prefactor\t")
                )
            )
            assert header["source_index"] == folder.name.split("-")[1]
            assert header["requested_prefactor"] == receipt["requested_prefactor"]
            for stage in ("source", "comparison", "representative"):
                diagram = json.loads((folder / f"{stage}.json").read_text())
                assert diagram["model"] == expected_bytes
                assert diagram["name"] == header["source_name"]
                checked_pairs += 1
        final = json.loads((root / f"{side}-frozen/{case}.json").read_text())
        for diagram in final["diagrams"]:
            assert diagram["payload"]["model"] == expected_bytes
        assert {record["name"] for record in grouped} == {
            diagram["name"] for diagram in final["diagrams"]
        }
        rows.append(
            {
                "side": side,
                "case": case,
                "fingerprint": computed,
                "full_model_sha256": hashlib.sha256(source.encode()).hexdigest(),
                "parameters": parameters,
                "propagators": propagators,
                "first_build_names_and_bases": len(first),
                "second_build_names_and_bases": len(second),
                "grouped_owner_fingerprints": len(grouped),
                "pair_stage_graph_fingerprints": checked_pairs,
                "final_graph_fingerprints": len(final["diagrams"]),
                "kernel_fingerprint_probe_executions": 1,
            }
        )
assert all(source == model_bytes[0] for source in model_bytes)
summary = {
    "rows": rows,
    "all_six_effective_model_definitions_literal_equal": True,
    "authoritative_contract": "Model::fingerprint hashes the canonical ModelDefinition with BLAKE3; DiagramSerde.model stores its 32 bytes. Each producer probe reparses the full effective JSON with its own locked kernel and recomputes that fingerprint.",
}
(root / "model-binding-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
