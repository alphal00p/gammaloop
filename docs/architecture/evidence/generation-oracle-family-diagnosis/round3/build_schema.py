"""Declare the fixed supported producer shapes, independent of mutable anchors.

This is an observed finite-domain contract, not a universal JSON schema. Actual
producer types are EdgePayload(EdgeEndpoints, EdgeDataPayload), PhysicalPayload,
Entry, Group/Member in generation_oracle/census.rs and the retained adapter.
"""

import collections

import common as C

types = collections.defaultdict(set)
objects = collections.defaultdict(set)
inputs = {}


def walk(value, path):
    types[path].add(type(value).__name__)
    if isinstance(value, dict):
        numeric_map = path.endswith("/edge_signatures")
        objects[path].add("numeric-map" if numeric_map else tuple(sorted(value)))
        for key, child in value.items():
            walk(child, path + "/" + ("#" if numeric_map else key))
    elif isinstance(value, list):
        for i, child in enumerate(value):
            walk(
                child,
                path + "/" + (str(i) if path.endswith("/payload/edges/*") else "*"),
            )


for root in C.B.ROOTS.values():
    for case in C.B.CASES:
        for mode in C.B.MODES:
            folder = root / case / f"classes-{mode}-t1"
            for name in (
                "payload.json",
                "typed-owners.json",
                "model.json",
                "request.json",
                "options.json",
                "grouping.json",
            ):
                path = folder / name
                inputs[str(path)] = C.sha(path)
                walk(C.load(path), name)
C.write(
    C.HERE / "schema.json",
    {
        "schema": 1,
        "scope": list(C.B.CASES),
        "types": {k: sorted(v) for k, v in sorted(types.items())},
        "object_keys": {
            k: sorted([list(x) if isinstance(x, tuple) else x for x in v], key=str)
            for k, v in sorted(objects.items())
        },
        "derivation_inputs": inputs,
        "rules": [
            "exact field kinds; Python bool is never int",
            "only declared object signatures",
            "numeric edge-signature maps only",
            "positional edge tuples have exactly two elements",
            "variable graph/owner cardinality is checked symmetrically across routes",
        ],
    },
)
print("declared paths", len(types), "object paths", len(objects), "inputs", len(inputs))
