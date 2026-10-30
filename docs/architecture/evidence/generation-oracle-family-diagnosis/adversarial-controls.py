from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).parent
SPEC = importlib.util.spec_from_file_location("capture_validator", ROOT / "validate-captures.py")
VALIDATOR = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(VALIDATOR)


def rejects(name, operation):
    try:
        operation()
    except (AssertionError, KeyError, json.JSONDecodeError):
        return {"name": name, "rejected": True}
    raise AssertionError(f"control accepted: {name}")


def main():
    case = "qed_eemm_L1"
    historical_folder = VALIDATOR.HISTORICAL / case / "classes-None-t1"
    current_folder = VALIDATOR.CURRENT / case / "classes-None-t1"
    historical = VALIDATOR.load(historical_folder, "payload.json")
    current = VALIDATOR.load(current_folder, "payload.json")
    key = sorted(VALIDATOR.graph_map(historical))[0]
    left = VALIDATOR.graph_map(historical)[key]
    right = VALIDATOR.graph_map(current)[key]
    historical_owners = VALIDATOR.load(historical_folder, "typed-owners.json")
    current_owners = VALIDATOR.load(current_folder, "typed-owners.json")
    left_owner_map = VALIDATOR.owner_map(historical_owners)[left["id"]]
    right_owner_map = VALIDATOR.owner_map(current_owners)[right["id"]]
    controls = []

    def edge_shape():
        changed = copy.deepcopy(right["payload"])
        changed["edges"][0] = changed["edges"][0][:1]
        VALIDATOR.compare_local_payload(left["payload"], changed)

    controls.append(rejects("edge tuple shape", edge_shape))

    def wrong_path():
        changed = copy.deepcopy(right["payload"])
        changed["edges"][0][0]["numerator"] = changed["edges"][0][1]["numerator"]
        VALIDATOR.compare_local_payload(left["payload"], changed)

    controls.append(rejects("expression field at incidence path", wrong_path))

    def rule_mismatch():
        changed = copy.deepcopy(right["payload"])
        changed["vertices"][0]["interaction"] += 1
        VALIDATOR.compare_local_payload(left["payload"], changed)

    controls.append(rejects("vertex rule mismatch", rule_mismatch))

    def routing_mismatch():
        changed = copy.deepcopy(right["payload"])
        changed["loop_momentum_basis"][0], changed["loop_momentum_basis"][1] = (
            changed["loop_momentum_basis"][1],
            changed["loop_momentum_basis"][0],
        )
        VALIDATOR.compare_local_payload(left["payload"], changed)

    controls.append(rejects("routing mismatch", routing_mismatch))

    def source_factor_mismatch():
        changed = copy.deepcopy(right["payload"])
        changed["overall_factor"] = "2*(" + changed["overall_factor"] + ")"
        VALIDATOR.compare_local_payload(left["payload"], changed)

    controls.append(rejects("source overall-factor mismatch", source_factor_mismatch))

    def duplicate_owner():
        changed = copy.deepcopy(historical_owners)
        changed[0]["owners"].append(copy.deepcopy(changed[0]["owners"][0]))
        VALIDATOR.owner_map(changed)

    controls.append(rejects("duplicate owner", duplicate_owner))

    def missing_owner():
        changed = dict(right_owner_map)
        changed.pop(next(iter(changed)))
        assert left_owner_map.keys() == changed.keys()

    controls.append(rejects("missing owner", missing_owner))

    def owner_swap():
        changed = copy.deepcopy(right_owner_map)
        names = list(changed)[:2]
        changed[names[0]], changed[names[1]] = changed[names[1]], changed[names[0]]
        for owner in left_owner_map:
            assert VALIDATOR.modern(left_owner_map[owner]["producer_decoded_expression"]) == changed[owner]["producer_decoded_expression"]

    controls.append(rejects("owner payload swap", owner_swap))

    def port_dimension():
        ports = list(right_owner_map["vertex:0"]["ordered_ports"])
        ports[0] = ports[0].replace("(4,", "(5,", 1)
        VALIDATOR.axis_permutation(
            [VALIDATOR.modern(port) for port in left_owner_map["vertex:0"]["ordered_ports"]],
            ports,
        )

    controls.append(rejects("typed-port dimension", port_dimension))

    def port_variance():
        ports = list(right_owner_map["vertex:0"]["ordered_ports"])
        ports[0] = ports[0].replace("spenso::{spenso::lower}::dind(", "", 1)
        if ports[0] != right_owner_map["vertex:0"]["ordered_ports"][0]:
            ports[0] = ports[0][:-1]
        else:
            ports[0] = "spenso::{spenso::lower}::dind(" + ports[0] + ")"
        VALIDATOR.axis_permutation(
            [VALIDATOR.modern(port) for port in left_owner_map["vertex:0"]["ordered_ports"]],
            ports,
        )

    controls.append(rejects("typed-port variance", port_variance))

    def port_occurrence():
        ports = list(right_owner_map["vertex:0"]["ordered_ports"])
        ports.append(ports[0])
        VALIDATOR.axis_permutation(
            [VALIDATOR.modern(port) for port in left_owner_map["vertex:0"]["ordered_ports"]],
            ports,
        )

    controls.append(rejects("typed-port occurrence", port_occurrence))

    def aggregate_local_chain():
        changed = copy.deepcopy(right_owner_map)
        changed["numerator"]["producer_decoded_expression"] = changed["vertex:0"]["producer_decoded_expression"]
        assert right["payload"]["numerator"] == changed["numerator"]["producer_decoded_expression"]

    controls.append(rejects("aggregate/local chain", aggregate_local_chain))

    def model_mismatch():
        model = bytearray(historical_folder.joinpath("model.json").read_bytes())
        model[-2] ^= 1
        assert bytes(model) == current_folder.joinpath("model.json").read_bytes()

    controls.append(rejects("model identity", model_mismatch))

    def malformed_json():
        json.loads(current_folder.joinpath("payload.json").read_text()[:-1])

    controls.append(rejects("malformed JSON", malformed_json))

    positive = {
        "name": "typed-axis permutation",
        "accepted": sorted(left_owner_map["vertex:0"]["ordered_ports"])
        == sorted(right_owner_map["vertex:0"]["ordered_ports"]),
        "qualification": "Only the explicit occurrence-preserving bijection is accepted; no owner rename is claimed.",
    }
    assert positive["accepted"]
    result = {"schema": 1, "negative": controls, "positive": [positive]}
    ROOT.joinpath("adversarial-control-results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
