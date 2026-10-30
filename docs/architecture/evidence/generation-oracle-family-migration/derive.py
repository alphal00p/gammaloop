"""Independent, source-matched Census::Summary byte encoding (no algebra)."""

import collections
import copy
import json
import re
from fractions import Fraction

import common as M

OFFSET, PRIME, MASK = 0xCBF29CE484222325, 0x100000001B3, (1 << 64) - 1


def digest(data, state=OFFSET):
    for byte in data:
        state = ((state ^ byte) * PRIME) & MASK
    return state


def number(n, state=OFFSET):
    return digest(n.to_bytes(8, "little"), state)


def text(value, state):
    return digest(value.encode() + b"\xff", state)


def json_bytes(value):
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode()


def multiset(values):
    state = OFFSET
    for value in sorted(values):
        state = number(value, state)
    return state


def fragment(payload):
    p = copy.deepcopy(payload)
    n = OFFSET
    for field in ("numerator", "numerator_prefactor", "projector"):
        n = digest(json_bytes(p.pop(field)), n)
    for field in (
        "num_interface",
        "numerator_prefactor_interface",
        "projector_interface",
    ):
        if field in p:
            n = digest(json_bytes(p.pop(field)), n)
    for local in [*p["vertices"], *[edge[1] for edge in p["edges"]]]:
        n = digest(json_bytes(local["numerator"]), n)
        local["numerator"] = None
        if "num_interface" in local:
            n = digest(json_bytes(local.pop("num_interface")), n)
    for field in ("model", "symmetry_factor", "overall_factor"):
        p.pop(field, None)
    return digest(json_bytes(p)), n


def weight(factor):
    # These four fixed ungrouped domains contain only exact AutG rational and
    # fermion +/-1 factors. Do not substitute arbitrary parameters.
    aut = re.search(r"AutG\((\d+)\)", factor)
    M.require(aut is not None, "weight", factor)
    signs = re.findall(
        r"(?:ExternalFermionOrderingSign|InternalFermionLoopSign)\((-?1)\)", factor
    )
    residual = re.sub(
        r"feynkit_generator_factor::\{\}::(?:AutG\(\d+\)|(?:ExternalFermionOrderingSign|InternalFermionLoopSign)\(-?1\))",
        "x",
        factor,
    )
    M.require(re.fullmatch(r"\(x\)\^\(-1\)(?:\*x)*", residual), "weight", residual)
    value = Fraction(1, int(aut[1]))
    for sign in signs:
        value *= int(sign)
    return value


def summary(payload, report):
    content, nums, ids, named, inverse, total, hist = (
        [],
        [],
        [],
        OFFSET,
        Fraction(0),
        Fraction(0),
        collections.Counter(),
    )
    for g in payload["diagrams"]:
        key = digest(json_bytes(g["key"]))
        frame, numerator = fragment(g["payload"])
        aut = int(g["aut"])
        inverse += Fraction(1, aut)
        total += weight(g["overall_factor"])
        hist[aut] += 1
        content.append(
            number(frame, text(g["overall_factor"], number(aut, number(key))))
        )
        nums.append(number(numerator, number(key)))
        named = number(
            key,
            text(g["name"], named) if g["name"] is not None else digest(b"\xfe", named),
        )
        ids.append(digest(int(g["id"], 16).to_bytes(16, "little")))
    for group in payload["groups"]:
        named = number(group["master"], named)
        for member in group["members"]:
            named = (
                text(member["name"], named)
                if member["name"] is not None
                else digest(b"\xfe", named)
            )
            named = text(member["overall_factor"], text(member["ratio"], named))
    counters = {
        k: int(v)
        for k, v in re.findall(
            r"(topology_count|interaction_assignment_count): (\d+)", report
        )
    }
    return {
        "diagrams": len(payload["diagrams"]),
        "topologies": counters["topology_count"],
        "assignments": counters["interaction_assignment_count"],
        "inverse_symmetry": [inverse.numerator, inverse.denominator],
        "weights": [total.numerator, total.denominator],
        "symmetry_factors": [list(p) for p in sorted(hist.items())],
        "content": content and multiset(content),
        "numerators": multiset(nums),
        "named": named,
        "ids": multiset(ids),
    }


def row(source, case):
    line = next(
        line
        for line in source.splitlines()
        if re.search(r"\b" + case + r": Row \{", line)
    )
    result = {
        k: int(v, 0) for k, v in re.findall(r"(\w+): (0x[0-9a-f]+|\d+)(?=[, }])", line)
    }
    for field in ("inverse_symmetry", "weights"):
        result[field] = [
            int(v) for v in re.search(field + r": \((-?\d+), (\d+)\)", line).groups()
        ]
    result["symmetry_factors"] = [
        [int(a), int(b)]
        for a, b in re.findall(
            r"\((\d+), (\d+)\)", re.search(r"symmetry_factors: &\[(.*?)\]", line)[1]
        )
    ]
    return result


def main():
    source, record = M.immutable(M.APPROVED, M.MAIN)
    historical_source, historical_record = M.immutable(
        "2b1ef4f07a44bc061bcc97c378bbdd9c41cb0e94", M.MAIN
    )
    records = []
    for case in M.CASES:
        roots = {
            "historical": M.C.B.ROOTS["historical"],
            "current": M.C.B.ROOTS["current"],
        }
        rows = {}
        for route, root in roots.items():
            d = root / case / "classes-None-t1"
            rows[route] = summary(
                M.load(d / "payload.json"), (d / "built-report.txt").read_text()
            )
        expected = row(source, case)
        M.require(row(historical_source, case) == expected, "source-row", case)
        M.require(rows["historical"] == expected, "derive", case + ": historical Row")
        M.require(
            {k: v for k, v in rows["current"].items() if k != "numerators"}
            == {k: v for k, v in expected.items() if k != "numerators"},
            "derive",
            case + ": non-numerator Row",
        )
        records.append(
            {
                "case": case,
                "expected": expected,
                **rows,
                "old_literal": f"0x{expected['numerators']:016x}",
                "new_literal": f"0x{rows['current']['numerators']:016x}",
            }
        )
    M.write(
        M.OUT / "derivation.json",
        {
            "rows": records,
            "source": record,
            "historical_source": historical_record,
            "algorithm": "Census::Summary FNV-1a; exact sorted JSON and source field order; no semantic equality inferred from hashes",
        },
    )
    print([(r["case"], r["old_literal"], r["new_literal"]) for r in records])


if __name__ == "__main__":
    main()
