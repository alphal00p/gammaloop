#!/usr/bin/env python3
"""Write the frozen inputs of the p-integral benchmark grid.

    python3 tools/p_benchmark/generate_inputs.py [--check]
    python3 tools/p_benchmark/generate_inputs.py --five-loop [TOP,...] [--check]

from `crates/vakint` (or any directory: paths are relative to this file). The
output is `tests/fixtures/benchmark_p_integrals.json`; `--check` only verifies
that the committed file is what this script writes.

Grid (M2 plan section 9, user Q1 defaults): L = 1..4 loops x dots {0, 1, 2} x
numerator rank r in {0, 2, 4, 6, 8}. One input per cell: the sum over the top
families of that L (1/1/3/11, the registry `P_FAMILIES` of `src/p_integral.rs`,
in registry order) of `(N_r)*top(dots)`, every top in its canonical long form
(the family lines with `pow(i)` replaced by the line power, external `p(1)`).

- dots 0: every line has power 1; dots 1: power 2 on line 1; dots 2: power 2 on
  line 1 and on the last line.
- N_r is the product of r/2 rank-2 blocks. Block b (0-based) uses the loop
  pair pairs[b % len(pairs)] with pairs = [(i, j), 1 <= i < j <= L] in
  lexicographic order, and is `dot(k(i),p(1))*dot(k(j),p(1))` for even b and
  `dot(k(i),k(j))` for odd b. At one loop the only pair is (1, 1) and every
  block is `dot(k(1),p(1))*dot(k(1),p(1))`: `dot(k(1),k(1))` would cancel
  line 1 and make the one-loop cells scaleless.
- N_0 = 1.

The strings are frozen: the benchmark, the Forcer grid fixtures and the
exact-agreement marks all refer to them. Never edit the JSON by hand.

`--five-loop` (t162-export, V7 preparation) writes the five-loop cells to
`tests/fixtures/benchmark_p_5l.json` instead (the frozen grid is untouched):
L = 5 x the same dots and ranks, one input per cell and top (default the s1
package's top `P3dc41bbf2115e501_e15`; a comma-separated list of registry
names otherwise), in the SHORT form `topo(NAME(p(1), a_1, ..., a_14))`: a
five-loop long form is not matched onto the registry yet (plan V1b). RustRed
only (Forcer has no five-loop topologies); the values end on PResidual
terminals, so the harness must run them with master_mode=retain.
"""
import json
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent
CRATE = HERE.parent.parent
REGISTRY = CRATE / "src" / "p_integral.rs"
OUTPUT = CRATE / "tests" / "fixtures" / "benchmark_p_integrals.json"
OUTPUT_5L = CRATE / "tests" / "fixtures" / "benchmark_p_5l.json"
DEFAULT_5L_TOPS = ["P3dc41bbf2115e501_e15"]

LOOPS = [1, 2, 3, 4]
DOTS = [0, 1, 2]
RANKS = [0, 2, 4, 6, 8]


def families():
    """(name, loops, long form) of every registry family, in registry order."""
    source = REGISTRY.read_text()
    entries = re.findall(
        r'PFamilyData \{\s*name: "(\w+)",\s*loops: (\d+),.*?lines: "([^"]+)"',
        source,
        re.S,
    )
    # 16 families of one to four loops, then 64 of five loops (t162-export).
    if len(entries) != 80 or sum(1 for _, loops, _ in entries if int(loops) <= 4) != 16:
        raise SystemExit(f"expected 16 + 64 registry families in {REGISTRY}, found {len(entries)}")
    return [(name, int(loops), lines) for name, loops, lines in entries]


def powers(loops, dots):
    lines = 3 * loops - 1
    result = [1] * lines
    if dots >= 1:
        result[0] = 2
    if dots >= 2:
        result[-1] = 2
    return result


def long_form(lines, power_list):
    def substitute(match):
        return str(power_list[int(match.group(1)) - 1])

    text, count = re.subn(r"pow\((\d+)\)", substitute, lines)
    if count != len(power_list):
        raise SystemExit(f"{lines}: {count} pow slots for {len(power_list)} powers")
    return text


def numerator(loops, rank):
    if rank == 0:
        return "1"
    pairs = [(i, j) for i in range(1, loops + 1) for j in range(i + 1, loops + 1)] or [(1, 1)]
    blocks = []
    for b in range(rank // 2):
        i, j = pairs[b % len(pairs)]
        if b % 2 == 0 or loops == 1:
            blocks.append(f"dot(k({i}),p(1))*dot(k({j}),p(1))")
        else:
            blocks.append(f"dot(k({i}),k({j}))")
    return "*".join(blocks)


def cells():
    registry = families()
    result = []
    for loops in LOOPS:
        tops = [(name, lines) for name, l, lines in registry if l == loops]
        for dots in DOTS:
            for rank in RANKS:
                p = powers(loops, dots)
                n = numerator(loops, rank)
                result.append(
                    {
                        "name": f"L{loops}_d{dots}_r{rank}",
                        "loops": loops,
                        "dots": dots,
                        "rank": rank,
                        "numerator": n,
                        "tops": [{"family": name, "powers": p} for name, _ in tops],
                        "input": "+".join(f"({n})*{long_form(lines, p)}" for _, lines in tops),
                    }
                )
    return result


def cells_5l(tops):
    registry = {name: int(loops) for name, loops, _ in families()}
    for top in tops:
        if registry.get(top) != 5:
            raise SystemExit(f"{top} is no five-loop registry family")
    result = []
    for top in tops:
        for dots in DOTS:
            for rank in RANKS:
                p = powers(5, dots)
                n = numerator(5, rank)
                result.append(
                    {
                        "name": f"L5_{top}_d{dots}_r{rank}",
                        "loops": 5,
                        "dots": dots,
                        "rank": rank,
                        "numerator": n,
                        "tops": [{"family": top, "powers": p}],
                        "input": f"({n})*topo({top}(p(1),{','.join(map(str, p))}))",
                    }
                )
    return result


def document_5l(tops):
    return {
        "schema": "vakint.p-benchmark-inputs.v1",
        "generator": "crates/vakint/tools/p_benchmark/generate_inputs.py --five-loop",
        "grid": {"loops": [5], "dots": DOTS, "ranks": RANKS, "tops": tops},
        "master_mode": "retain",
        "cells": cells_5l(tops),
    }


def document():
    return {
        "schema": "vakint.p-benchmark-inputs.v1",
        "generator": "crates/vakint/tools/p_benchmark/generate_inputs.py",
        "grid": {"loops": LOOPS, "dots": DOTS, "ranks": RANKS},
        "cells": cells(),
    }


def main():
    arguments = sys.argv[1:]
    output, document_text = OUTPUT, None
    if "--five-loop" in arguments:
        index = arguments.index("--five-loop")
        following = arguments[index + 1] if index + 1 < len(arguments) else ""
        tops = following.split(",") if following and not following.startswith("--") else DEFAULT_5L_TOPS
        output, document_text = OUTPUT_5L, json.dumps(document_5l(tops), indent=2) + "\n"
    text = document_text or json.dumps(document(), indent=2) + "\n"
    if "--check" in arguments:
        if output.read_text() != text:
            raise SystemExit(f"{output} differs from the generator output")
        print(f"{output}: up to date")
        return
    output.write_text(text)
    print(f"wrote {output} ({len(text)} bytes)")


if __name__ == "__main__":
    main()
