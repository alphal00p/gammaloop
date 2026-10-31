"""Factual types of this pinned Auto ledger, before its frozen math predicates.

This is the observed serializer shape, not a schema for arbitrary future captures.
It deliberately leaves the original sealed helpers unchanged.
"""

from collections import Counter


class Reject(ValueError):
    def __init__(self, predicate, detail):
        self.predicate = predicate
        super().__init__(f"{predicate}: {detail}")


def need(condition, predicate, detail):
    if not condition:
        raise Reject(predicate, detail)


def check(rows):
    stages = {}
    for row in rows:
        need(type(row) is dict, "strict-event-shape", "event object")
        need("stage" in row and "value" in row, "strict-event-shape", "stage/value")
        stages.setdefault(row["stage"], []).append(row["value"])
    for stage, count in [
        ("sample-zero", 30),
        ("prepare-zero", 30),
        ("group-cancellation", 29),
        ("comparison-decision", 3),
    ]:
        need(len(stages.get(stage, [])) == count, "strict-stage-count", stage)
        need(all(type(x) is dict for x in stages[stage]), "strict-event-shape", stage)
    for value in stages["sample-zero"]:
        need(
            "sample_count" in value
            and type(value["sample_count"]) is int
            and value["sample_count"] == 5,
            "strict-sample-count",
            "all 30 sample-zero.sample_count are actual integer 5",
        )
    for stage in ["prepare-zero", "sample-zero"]:
        for value in stages[stage]:
            need(
                "zero" in value and type(value["zero"]) is bool and value["zero"] is False,
                "strict-zero-boolean",
                stage,
            )
    for value in stages["group-cancellation"]:
        need(
            "cancelled" in value
            and type(value["cancelled"]) is bool
            and value["cancelled"] is False,
            "strict-cancellation-boolean",
            "all 29 cancellation flags are actual false",
        )
    members = []

    def walk(value, path):
        if type(value) is list:
            for i, child in enumerate(value):
                walk(child, f"{path}/{i}")
        elif type(value) is dict:
            # The two actual member serializers use diagram or id/name. The
            # separate collapse-plan.members is a list of names, not records.
            if {"source_index", "overall_factor"} <= set(value):
                need(
                    "ratio" in value
                    and ("diagram" in value or {"id", "name"} <= set(value)),
                    "strict-member-shape",
                    path,
                )
                need(
                    type(value["ratio"]) is str and value["ratio"] == "1",
                    "strict-member-ratio",
                    "the observed member serializer supplies a nonnull string 1",
                )
                members.append(path)
            for key, child in value.items():
                walk(child, f"{path}/{key}")

    walk(rows, "ledger")
    need(len(members) == 210, "strict-member-shape", "all seven actual 30-member record populations")
    ratios = []
    for value in stages["comparison-decision"]:
        need("ratio" in value, "strict-comparison-ratio", "explicit nullable ratio")
        ratio = value["ratio"]
        need(
            ratio is None or (type(ratio) is str and ratio == "1"),
            "strict-comparison-ratio",
            "actual domain is null or string 1; integers/booleans are not ratios",
        )
        ratios.append(ratio)
    need(
        Counter(ratios) == Counter([None, None, "1"]),
        "strict-comparison-ratio",
        "exactly two sampled nonmerges and one ratio-one decision",
    )
    return {
        "sample_counts_integer_5": 30,
        "zero_flags_boolean_false": 60,
        "cancellation_flags_boolean_false": 29,
        "nonnull_member_ratio_string_1_occurrences": len(members),
        "comparison_ratios": ratios,
        "scope": "Factual guard for this source-defined shape; no general capture validation approval.",
    }
