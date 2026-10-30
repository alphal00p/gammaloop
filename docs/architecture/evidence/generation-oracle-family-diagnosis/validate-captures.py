from __future__ import annotations

import hashlib
import json
import sys
from collections import defaultdict, deque
from pathlib import Path

HISTORICAL = Path(
    "/shared/lcnbr/dev/gammaloop/generation-oracle-family-historical/target/"
    "oracle-numerator-family-diagnosis/adapted-historical/captures"
)
CURRENT = Path(
    "/shared/lcnbr/dev/gammaloop/generation-oracle-numerator-families/target/"
    "oracle-numerator-family-diagnosis/adapted-current/captures"
)
OUTPUT = Path(
    "/shared/lcnbr/dev/gammaloop/generation-oracle-numerator-families/target/"
    "oracle-numerator-family-diagnosis/proof"
)
ALLOWED_CASES = {"qed_eemm_L1", "ctorder_qcd_uugg_L1", "sm_eemm_L0", "ym_g6_L0"}
MODES = ("None", "Auto")

# These are exact old/current constructor spellings verified in the immutable
# feynkit-graph symbols.rs definitions. They are not a general normalization.
SPELLINGS = {
    "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:q}::Q":
        "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript}::Q",
    "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:k}::K":
        "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript}::K",
    "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:p}::P":
        "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript}::P",
    "gammalooprs::{spenso::index,spenso::index-label:e}::edge":
        "gammalooprs::{spenso::index}::edge",
    "gammalooprs::{spenso::index,spenso::index-label:v}::vertex":
        "gammalooprs::{spenso::index}::vertex",
}


def encoded(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def modern(expression: str) -> str:
    for historical, current in SPELLINGS.items():
        expression = expression.replace(historical, current)
    return expression


def load(folder: Path, name: str) -> object:
    with (folder / name).open() as stream:
        return json.load(stream)


def graph_map(payload: dict) -> dict[str, dict]:
    result = {}
    for diagram in payload["diagrams"]:
        key = encoded(diagram["key"])
        assert key not in result, "duplicate full canonical key"
        result[key] = diagram
    assert len(result) == len(payload["diagrams"])
    return result


def owner_map(records: list[dict]) -> dict[str, dict[str, dict]]:
    result = {}
    for diagram in records:
        assert diagram["id"] not in result, "duplicate typed-owner diagram ID"
        owners = {}
        for owner in diagram["owners"]:
            assert owner["owner"] not in owners, "duplicate owner path"
            owners[owner["owner"]] = owner
        result[diagram["id"]] = owners
    return result


def axis_permutation(historical: list[str], current: list[str]) -> list[int]:
    positions = defaultdict(deque)
    for index, port in enumerate(current):
        positions[port].append(index)
    permutation = []
    for port in historical:
        assert positions[port], f"missing typed port {port}"
        permutation.append(positions[port].popleft())
    assert all(not remaining for remaining in positions.values()), "extra typed port"
    assert sorted(permutation) == list(range(len(current))), "non-bijective typed ports"
    return permutation


def compare_local_payload(historical: dict, current: dict) -> list[tuple[str, str, str, object, object]]:
    exact_fields = [
        "cuts",
        "half_edge_order",
        "loop_momentum_basis",
        "model",
        "numerator_prefactor",
        "overall_factor",
        "projector",
        "symmetry_factor",
        "topology_threshold_candidates",
    ]
    for field in exact_fields:
        assert historical[field] == current[field], f"payload field {field}"
    assert len(historical["vertices"]) == len(current["vertices"])
    assert len(historical["edges"]) == len(current["edges"])
    pairs = [
        ("numerator_prefactor", historical["numerator_prefactor"], current["numerator_prefactor"], historical.get("numerator_prefactor_interface"), current.get("numerator_prefactor_interface")),
        ("projector", historical["projector"], current["projector"], historical.get("projector_interface"), current.get("projector_interface")),
        ("numerator", historical["numerator"], current["numerator"], historical.get("num_interface"), current.get("num_interface")),
    ]
    for index, (left, right) in enumerate(zip(historical["vertices"], current["vertices"])):
        left_fixed = {key: value for key, value in left.items() if key not in {"numerator", "num_interface"}}
        right_fixed = {key: value for key, value in right.items() if key not in {"numerator", "num_interface"}}
        assert left_fixed == right_fixed, f"vertex:{index} physical/rule payload"
        pairs.append((f"vertex:{index}", left["numerator"], right["numerator"], left.get("num_interface"), right.get("num_interface")))
    for index, (left, right) in enumerate(zip(historical["edges"], current["edges"])):
        assert isinstance(left, list) and isinstance(right, list)
        assert len(left) == len(right) == 2, f"edge:{index} positional tuple"
        assert left[0] == right[0], f"edge:{index} incidence"
        left_fixed = {key: value for key, value in left[1].items() if key not in {"numerator", "num_interface"}}
        right_fixed = {key: value for key, value in right[1].items() if key not in {"numerator", "num_interface"}}
        assert left_fixed == right_fixed, f"edge:{index} physical payload"
        pairs.append((f"edge:{index}", left[1]["numerator"], right[1]["numerator"], left[1].get("num_interface"), right[1].get("num_interface")))
    return pairs


def validate(case: str, mode: str) -> dict:
    historical = HISTORICAL / case / f"classes-{mode}-t1"
    current = CURRENT / case / f"classes-{mode}-t1"
    for name in ["model.json", "request.json", "options.json", "grouping.json"]:
        assert historical.joinpath(name).read_bytes() == current.joinpath(name).read_bytes(), name
    assert load(historical, "grouping.json") == "none"
    request = load(historical, "request.json")
    assert request["classes"].lower() == mode.lower()

    historical_payload = load(historical, "payload.json")
    current_payload = load(current, "payload.json")
    historical_graphs = graph_map(historical_payload)
    current_graphs = graph_map(current_payload)
    assert historical_graphs.keys() == current_graphs.keys(), "complete canonical-key domain"

    historical_owners = owner_map(load(historical, "typed-owners.json"))
    current_owners = owner_map(load(current, "typed-owners.json"))
    assert historical_owners.keys() == current_owners.keys(), "typed-owner diagram domain"

    local_unique = {}
    graph_pairs = []
    axis_records = []
    local_records = []
    for key in sorted(historical_graphs):
        left = historical_graphs[key]
        right = current_graphs[key]
        for field in ["aut", "basis", "id", "key", "name", "overall_factor"]:
            assert left[field] == right[field], f"graph {left['id']} field {field}"
        payload_pairs = compare_local_payload(left["payload"], right["payload"])
        left_owners = historical_owners[left["id"]]
        right_owners = current_owners[right["id"]]
        assert left_owners.keys() == right_owners.keys(), f"owner domain {left['id']}"
        for owner in sorted(left_owners):
            old = left_owners[owner]
            new = right_owners[owner]
            assert old["is_metric"] == new["is_metric"]
            assert old["is_composite"] == new["is_composite"]
            old_ports = [modern(port) for port in old["ordered_ports"]]
            permutation = axis_permutation(old_ports, new["ordered_ports"])
            old_expression = modern(old["producer_decoded_expression"])
            new_expression = new["producer_decoded_expression"]
            pair_key = (old_expression, new_expression)
            proof = "literal" if old_expression == new_expression else local_unique.setdefault(
                pair_key, len(local_unique)
            )
            axis_records.append(
                {"id": left["id"], "owner": owner, "permutation": permutation, "ports": old_ports}
            )
            local_records.append({"id": left["id"], "owner": owner, "proof": proof})
        payload_by_owner = {owner: (old, new, old_interface, new_interface) for owner, old, new, old_interface, new_interface in payload_pairs}
        assert payload_by_owner.keys() == left_owners.keys(), f"serialized/typed owner paths {left['id']}"
        for owner, (old, new, old_interface, new_interface) in payload_by_owner.items():
            assert modern(old) == modern(left_owners[owner]["producer_decoded_expression"])
            assert new == right_owners[owner]["producer_decoded_expression"]
            if old_interface is not None:
                assert [modern(port) for port in old_interface] == [modern(port) for port in left_owners[owner]["ordered_ports"]]
            if new_interface is not None:
                assert new_interface == right_owners[owner]["ordered_ports"]
        graph_pairs.append(
            {
                "id": left["id"],
                "historical": modern(left["payload"]["numerator"]),
                "current": right["payload"]["numerator"],
            }
        )

    unique_pairs = [None] * len(local_unique)
    for (historical_expression, current_expression), index in local_unique.items():
        unique_pairs[index] = {
            "historical": historical_expression,
            "current": current_expression,
        }
    destination = OUTPUT / case / f"classes-{mode}-t1"
    destination.mkdir(parents=True, exist_ok=True)
    destination.joinpath("local-pairs.json").write_text(encoded(unique_pairs) + "\n")
    destination.joinpath("graph-pairs.json").write_text(encoded(graph_pairs) + "\n")
    destination.joinpath("axis-bijections.json").write_text(encoded(axis_records) + "\n")
    destination.joinpath("local-inventory.json").write_text(encoded(local_records) + "\n")

    model = load(historical, "model.json")
    ct_detail = None
    if case == "ctorder_qcd_uugg_L1":
        rules = [rule for rule in model["vertex_rules"] if rule["name"].endswith("_CT")]
        couplings = [coupling for coupling in model["couplings"] if coupling["name"].endswith("_CT")]
        assert rules and couplings
        assert all(["CT", 1] in coupling["orders"] for coupling in couplings)
        ct_detail = {
            "ct_rules": sorted(rule["name"] for rule in rules),
            "ct_couplings": sorted(coupling["name"] for coupling in couplings),
            "selected_vertex_interactions": sorted(
                {vertex["interaction"] for graph in historical_graphs.values() for vertex in graph["payload"]["vertices"]}
            ),
        }
    return {
        "case": case,
        "mode": mode,
        "graphs": len(historical_graphs),
        "owners": sum(len(owners) for owners in historical_owners.values()),
        "axes": sum(len(record["ports"]) for record in axis_records),
        "permuted_owners": sum(record["permutation"] != list(range(len(record["permutation"]))) for record in axis_records),
        "nonliteral_unique_local_pairs": len(unique_pairs),
        "input_sha256": {
            name: digest(historical / name)
            for name in ["model.json", "request.json", "options.json", "grouping.json"]
        },
        "ct_order_detail": ct_detail,
    }


def main() -> None:
    cases = sys.argv[1:]
    assert cases and set(cases) <= ALLOWED_CASES
    results = [validate(case, mode) for case in cases for mode in MODES]
    OUTPUT.mkdir(parents=True, exist_ok=True)
    OUTPUT.joinpath("capture-validation.json").write_text(
        json.dumps({"schema": 1, "results": results, "spellings": SPELLINGS}, indent=2) + "\n"
    )
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
