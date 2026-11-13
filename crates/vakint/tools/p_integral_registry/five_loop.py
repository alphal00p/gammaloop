#!/usr/bin/env python3
"""five_loop.py P_FAMILIES_5L_JSON: the 64 five-loop entries of vakint's P_FAMILIES (Rust source on stdout).

Input: RustRed's `crates/rustred-core/data/p/p-families-5L.v1.json` (schema `rustred.p-families.v1`). The
library's top presentation of each family has its line slots in the JSON order, so the entries keep the JSON's
line order and momenta exactly (vakint's `check_top` compares them up to an overall sign); only the vertices are
relabelled (q enters at vertex 1, leaves at vertex 2, the others 3.. in ascending order) and each line is oriented
so that momentum is conserved with that flow (vakint's `prop(i, edge(a, b), m)`: `m` leaves `a`). Unlike the 1-4
loop entries, line i <= L does NOT carry k(i) alone: the routing is RustRed's."""
import itertools, json, sys

def fmt(row, loops):
    out = ""
    for i, c in enumerate(row):
        if c == 0:
            continue
        name = f"k({i + 1})" if i < loops else "p(1)"
        out += ("+" if c > 0 and out else "") + ("-" if c < 0 else "") + (f"{abs(c)}*" if abs(c) != 1 else "") + name
    return out

def oriented(family):
    loops, lines = family["loops"], family["lines"]
    ext_in, ext_out = family["external"]["in"], family["external"]["out"]
    others = sorted({v for line in lines for v in line["endpoints"]} - {ext_in, ext_out})
    label = {ext_in: 1, ext_out: 2, **{v: i + 3 for i, v in enumerate(others)}}
    q = [0] * loops + [1]
    for signs in itertools.product((1, -1), repeat=len(lines)):
        balance = {}
        edges = []
        for sign, line in zip(signs, lines):
            a, b = (label[v] for v in line["endpoints"])
            if sign < 0:
                a, b = b, a
            edges.append((a, b, line["momentum"]))
            for vertex, factor in ((a, 1), (b, -1)):
                row = balance.setdefault(vertex, [0] * (loops + 1))
                for k, c in enumerate(line["momentum"]):
                    row[k] += factor * c
        expected = lambda v: q if v == 1 else [-c for c in q] if v == 2 else [0] * (loops + 1)
        if all(balance[v] == expected(v) for v in balance):
            return edges
    raise SystemExit(f"{family['name']}: no orientation conserves momentum")

data = json.load(open(sys.argv[1]))
assert data["schema"] == "rustred.p-families.v1", data["schema"]
families = data["families"]
assert len(families) == 64 and all(f["loops"] == 5 for f in families)
for family in families:
    edges = oriented(family)
    props = "*".join(f"prop({i + 1},edge({a},{b}),{fmt(m, 5)},0,pow({i + 1}))" for i, (a, b, m) in enumerate(edges))
    isps = ", ".join(f'"{fmt(row, 5)}"' for row in family["isps"])
    print("    PFamilyData {")
    print(f'        name: "{family["name"]}",')
    print("        loops: 5,")
    print("        forcer_topology: None,")
    print('        vacuum_parent: "6L, RustRed p-families-5L.v1",')
    print(f'        lines: "topo({props})",')
    print(f"        isps: &[{isps}],")
    print("    },")
