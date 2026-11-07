#!/usr/bin/env python3
"""Reduction graphs of an external-data root (benchmark.typ, Figures 1L-5L).

usage: reduction_graphs.py DATA_ROOT OUT.json

For every loop count L = 1..5 the graph is built only from the data root:
  * parent topologies: the family graphs of L/graphs (1L-4L) and the five-loop
    maximal parents L/graphs/parent-*.dot (5L);
  * nodes: vacuum topologies (graphs up to isomorphism) reachable from the
    parents by removing (pinching) lines, as in the reduction graph of
    B. Ruijl's thesis (Fig. 21); a line on a loop of its own (a one-loop
    tadpole factor) cannot be removed, its removal gives a scaleless zero;
  * arrows: "removing one line of A gives B", with the number of such lines;
  * node data: the shipped RustRed terminals (masters) whose graph.dot falls in
    the topology, factorisation into lower-loop blocks, the rule program
    that owns the topology (1L-3L closed artifact, 4L family program, 5L
    saved-application owner sector with its payload size), and a layered
    layout (layers = number of lines) for the Typst figure.
Standard library only; no rule table is decoded.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

EDGE = re.compile(r"^\s*v(\d+)\s*->\s*v(\d+)\s*(\[(.*)\])?\s*;?\s*$")
NODE = re.compile(r"^\s*v(\d+)\s*\[")
ATTR = re.compile(r'(\w+)=("(?:[^"\\]|\\.)*"|[^,\]\s]+)')


def parse_dot(path: Path):
    """Vertices and undirected edges (with attributes) of a DOT graph."""
    vertices, edges = set(), []
    for line in path.read_text().splitlines():
        m = EDGE.match(line)
        if m:
            a, b = int(m[1]), int(m[2])
            attrs = {k: v.strip('"') for k, v in ATTR.findall(m[4] or "")}
            vertices |= {a, b}
            edges.append((a, b, attrs))
            continue
        m = NODE.match(line)
        if m:
            vertices.add(int(m[1]))
    index = {v: i for i, v in enumerate(sorted(vertices))}
    return len(index), [(index[a], index[b], attrs) for a, b, attrs in edges]


# ----------------------------------------------------------- canonical form --

def adjacency(n, edges):
    a = [[0] * n for _ in range(n)]
    for u, v in edges:
        if u == v:
            a[u][u] += 1
        else:
            a[u][v] += 1
            a[v][u] += 1
    return a


def _refine(a, colors):
    n = len(a)
    while True:
        sig = [
            (colors[i], a[i][i],
             tuple(sorted((colors[j], a[i][j]) for j in range(n) if j != i and a[i][j])))
            for i in range(n)
        ]
        keys = sorted(set(sig))
        new = [keys.index(s) for s in sig]
        if len(keys) == len(set(colors)):
            return new
        colors = new


def canonical(n, edges):
    """Canonical certificate of an undirected multigraph with self-loops.

    Individualisation-refinement: refine by (colour, loops, neighbour colours
    with multiplicities); individualise every vertex of the first non-trivial
    cell; the lexicographically smallest relabelled adjacency is the form."""
    a = adjacency(n, edges)
    best = None

    def search(colors):
        nonlocal best
        colors = _refine(a, colors)
        cells = Counter(colors)
        if len(cells) == n:
            perm = colors  # vertex i -> position colors[i]
            inv = [0] * n
            for i, p in enumerate(perm):
                inv[p] = i
            cert = tuple(a[inv[i]][inv[j]] for i in range(n) for j in range(i, n))
            if best is None or cert < best:
                best = cert
            return
        target = min(c for c, k in cells.items() if k > 1)
        for v in [i for i in range(n) if colors[i] == target]:
            # individualise v: it keeps `target`, the rest of its cell moves up
            new = [2 * c + (1 if (c == target and i != v) else 0) for i, c in enumerate(colors)]
            search(new)

    search([0] * n)
    return (n, best)


def from_certificate(cert):
    n, flat = cert
    edges, k = [], 0
    for i in range(n):
        for j in range(i, n):
            edges += [(i, j)] * flat[k]
            k += 1
    return n, edges


def contract(n, edges, e):
    """Contract edge index e (not a self-loop); returns (n-1, edges)."""
    u, v = edges[e]
    keep, drop = min(u, v), max(u, v)
    out = []
    for i, (a, b) in enumerate(edges):
        if i == e:
            continue
        a = keep if a == drop else a
        b = keep if b == drop else b
        a = a - 1 if a > drop else a
        b = b - 1 if b > drop else b
        out.append((min(a, b), max(a, b)))
    return n - 1, out


def blocks(n, edges):
    """Biconnected blocks (as edge-index lists); every self-loop is a block."""
    simple = defaultdict(set)
    for i, (a, b) in enumerate(edges):
        if a != b:
            simple[a].add(b)
            simple[b].add(a)
    disc, low, stack, comps, t = {}, {}, [], [], [0]

    def dfs(u, parent):
        disc[u] = low[u] = t[0]
        t[0] += 1
        for w in simple[u]:
            if w == parent:
                continue
            if w not in disc:
                stack.append((u, w))
                dfs(w, u)
                low[u] = min(low[u], low[w])
                if low[w] >= disc[u]:
                    comp = set()
                    while True:
                        x = stack.pop()
                        comp |= set(x)
                        if x == (u, w):
                            break
                    comps.append(comp)
            elif disc[w] < disc[u]:
                stack.append((u, w))
                low[u] = min(low[u], disc[w])

    sys.setrecursionlimit(10000)
    for v in range(n):
        if v not in disc and simple[v]:
            dfs(v, None)
    result = []
    for comp in comps:
        result.append([i for i, (a, b) in enumerate(edges) if a != b and a in comp and b in comp])
    result += [[i] for i, (a, b) in enumerate(edges) if a == b]
    return result


def sub_certificate(n, edges, idx):
    verts = sorted({x for i in idx for x in edges[i]})
    m = {v: k for k, v in enumerate(verts)}
    return canonical(len(verts), [(m[edges[i][0]], m[edges[i][1]]) for i in idx])


def loops(n, edges):
    return len(edges) - n + 1


# --------------------------------------------------------------- the graph --

def topology_key(n, edges):
    """Node identity: irreducible graphs by their certificate; one-vertex-
    reducible graphs (vacuum products) by the multiset of their blocks, since a
    product integral does not depend on where its blocks are glued."""
    edges = [(min(a, b), max(a, b)) for a, b in edges]
    bl = blocks(n, edges)
    if len(bl) > 1:
        return ("product", tuple(sorted(sub_certificate(n, edges, b) for b in bl)))
    return ("graph", canonical(n, edges))


def build(parents, family_names):
    """BFS over topologies from the parents; returns nodes and arrows."""
    nodes, arrows, queue = {}, Counter(), []

    def visit(key, graph):
        if key not in nodes:
            nodes[key] = {"graph": graph, "parent_of": [], "family_names": []}
            queue.append(key)

    for label, (n, edges) in parents:
        edges = [(min(a, b), max(a, b)) for a, b in edges]
        key = topology_key(n, edges)
        visit(key, (n, edges))
        nodes[key]["parent_of"].append(label)
        if label in family_names:
            nodes[key]["family_names"].append(family_names[label])
    while queue:
        key = queue.pop()
        n, edges = nodes[key]["graph"]
        for e, (a, b) in enumerate(edges):
            if a == b:
                continue  # removing a tadpole line leaves a scaleless loop
            child = contract(n, edges, e)
            ckey = topology_key(*child)
            arrows[(key, ckey)] += 1
            visit(ckey, child)
    return nodes, arrows


SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")


def product_name(factor_names):
    counts = Counter(factor_names)
    order = sorted(counts, key=lambda s: (s != "T1", s))
    return "·".join(f + (str(counts[f]).translate(SUP) if counts[f] > 1 else "") for f in order)


def certificate_key(key):
    return hashlib.sha256(json.dumps(key).encode()).hexdigest()[:12]


BASIS_SYMBOL = re.compile(r"\b(PR\d+d?|miD6|miBN|miT111|miD5|miDM|miE3)\b")


# ------------------------------------------------------------------ layout --

def crossings(order, keys, arrows_by_parent):
    """Edge crossings between consecutive layers (edges join adjacent layers; weight 1 per edge)."""
    pos = {i: k for ln in order for k, i in enumerate(order[ln])}
    total = 0
    for ln in keys:
        edges = sorted((pos[u], pos[v]) for u in order[ln] for v in arrows_by_parent.get(u, ()))
        # count inversions of the child positions (pairs with u1 < u2 and v1 > v2)
        for a in range(len(edges)):
            ua, va = edges[a]
            for b in range(a + 1, len(edges)):
                ub, vb = edges[b]
                if ub != ua and vb < va:
                    total += 1
    return total


def layout(layers, arrows_by_child, arrows_by_parent, sweeps=400):
    """Layer ordering: barycentre sweeps (keeping the ordering with the fewest crossings), then adjacent
    transpositions while they remove crossings; then coordinate relaxation (unit spacing)."""
    order = {ln: sorted(ids) for ln, ids in layers.items()}
    keys = sorted(layers, reverse=True)
    pos = {i: k for ln in keys for k, i in enumerate(order[ln])}
    best, best_order = crossings(order, keys, arrows_by_parent), {k: list(v) for k, v in order.items()}
    for sweep in range(sweeps):
        seq = keys[1:] if sweep % 2 == 0 else list(reversed(keys[:-1]))
        for ln in seq:
            nbrs = arrows_by_child if sweep % 2 == 0 else arrows_by_parent

            def bary(i):
                ps = [pos[p] for p in nbrs.get(i, ())]
                return (sum(ps) / len(ps)) if ps else pos[i]

            order[ln].sort(key=lambda i: (bary(i), pos[i]))
            for k, i in enumerate(order[ln]):
                pos[i] = k
        c = crossings(order, keys, arrows_by_parent)
        if c < best:
            best, best_order = c, {k: list(v) for k, v in order.items()}
    order = best_order

    def local(ln):
        """Crossings of the two layer pairs touching layer ln."""
        sub = [k for k in keys if abs(k - ln) <= 1]
        return crossings({k: order[k] for k in sub}, [k for k in sub if k != min(sub)], arrows_by_parent)

    # sifting: move every node to the position of its layer with the fewest local crossings
    for _ in range(4):
        moved = False
        for ln in keys:
            for node in list(order[ln]):
                row = order[ln]
                row.remove(node)
                scores = []
                for k in range(len(row) + 1):
                    row.insert(k, node)
                    scores.append((local(ln), abs(k - len(row) / 2), k))
                    row.pop(k)
                k = min(scores)[2]
                row.insert(k, node)
                moved |= row.index(node) != k
        c = crossings(order, keys, arrows_by_parent)
        if c >= best and not moved:
            break
        best = min(best, c)
    layout.crossings = crossings(order, keys, arrows_by_parent)
    # x coordinates: centred, then relaxed towards neighbours keeping order
    x = {}
    width = max(len(v) for v in order.values())
    for ln in keys:
        m = len(order[ln])
        for k, i in enumerate(order[ln]):
            x[i] = k - (m - 1) / 2
    for it in range(60):
        for ln in (keys if it % 2 == 0 else list(reversed(keys))):
            ids = order[ln]
            want = []
            for i in ids:
                nb = list(arrows_by_child.get(i, ())) + list(arrows_by_parent.get(i, ()))
                want.append(sum(x[j] for j in nb) / len(nb) if nb else x[i])
            # place in order with unit gaps, as close as possible to `want`
            cur = [0.0] * len(ids)
            for k in range(len(ids)):
                cur[k] = want[k] if k == 0 else max(want[k], cur[k - 1] + 1.0)
            for k in range(len(ids) - 2, -1, -1):
                cur[k] = min(cur[k], cur[k + 1] - 1.0)
            shift = (sum(want) - sum(cur)) / len(ids)
            for k, i in enumerate(ids):
                x[i] = cur[k] + shift
    lo = min(x.values())
    hi = max(x.values())
    return order, {i: (x[i] - lo) for i in x}, max(hi - lo, 1.0), width


# -------------------------------------------------------------------- main --

def masters_of(root: Path, level: int):
    """(name, graph.dot path) of every shipped terminal of a loop count."""
    manifest = json.loads((root / f"{level}L/manifest.json").read_text())
    out = []
    for m in manifest["masters"]:
        name = m.get("name") or Path(m.get("directory", "")).name
        directory = m.get("directory") or f"masters/{name}"
        out.append((name, root / f"{level}L" / directory / "graph.dot", m))
    return out


def parents_of(root: Path, level: int):
    d = root / f"{level}L/graphs"
    if level == 5:
        files = sorted(d.glob("parent-*.dot"))
    else:
        files = sorted(d.glob("*.dot"))
    out = []
    for f in files:
        n, edges = parse_dot(f)
        out.append((f.stem, (n, [(a, b) for a, b, _ in edges])))
    return out


def family_names(root: Path, level: int):
    """Family graph file stem -> family alias (4L: x, fg, bmw, h)."""
    if level != 4:
        return {}
    manifest = json.loads((root / "4L/manifest.json").read_text())
    names = {}
    for m in manifest["masters"]:
        alias = (m.get("alias") or "").replace("4L_", "").upper()
        names.setdefault(m["family_sha256"], alias)
    return names


def owners_5l(root: Path):
    """5L saved-application owner sectors -> (parent stem, contracted parent slots, bytes)."""
    sel = json.loads((root / "5L/ibps/selection.json").read_text())
    src = (root / "provenance/five_loop_data.rs")
    census = []
    if src.exists():
        text = src.read_text()
        for block in re.findall(r"    Witness \{(.*?)\n    \}", text, re.DOTALL):
            census.append({
                "sector": int(re.search(r"sector: (\d+)", block)[1]),
                "parent": int(re.search(r"parent: (\d+)", block)[1]),
                "deleted": json.loads(re.search(r"deleted_slots: &(\[[^]]*\])", block)[1]),
                "source": json.loads(re.search(r"source_slots: &(\[[^]]*\])", block)[1]),
            })
        parent_sectors = [int(s) for s in re.findall(r"    Parent \{.*?sector: (\d+)", text, re.DOTALL)]
    else:
        parent_sectors = []
    out = []
    for o in sel["owners"]:
        # the owner's census representative (census sectors: slot i = bit 14 - i)
        w = next((w for w in census if w["sector"] == o["representative"]), None)
        out.append({
            "mask": o["mask"], "bytes": o["bytes"], "parent": o["parent"],
            "deleted": None if w is None else w["deleted"],
            "parent_stem": None if w is None or not parent_sectors else f"parent-{parent_sectors[w['parent']]}",
        })
    return out


def extract(root: Path):
    result = {"schema": "vakint.benchmark.reduction-graphs.v1", "loops": {}}
    lower_names = {}  # certificate of an irreducible graph -> name (lower loops)
    std = {1: {1: "T1"}, 2: {3: "T111"}, 3: {4: "BN", 5: "D5", 6: "D6"}}
    for level in range(1, 6):
        parents = parents_of(root, level)
        fam = family_names(root, level)
        nodes, arrows = build(parents, fam)
        names, info = {}, {}
        # irreducible topologies: families, standard names, else L·lines + letter
        per_layer = defaultdict(list)
        for key in nodes:
            n, edges = nodes[key]["graph"]
            info[key] = {"lines": len(edges), "vertices": n, "product": key[0] == "product"}
            if key[0] == "graph":
                per_layer[len(edges)].append(key)
        for ln, keys in per_layer.items():
            keys.sort()
            for k, key in enumerate(keys):
                if nodes[key]["family_names"]:
                    names[key] = "/".join(sorted(set(nodes[key]["family_names"])))
                elif level in std and ln in std[level] and len(keys) == 1:
                    names[key] = std[level][ln]
                elif info[key]["vertices"] == 2:
                    names[key] = "BN" if level == 3 else f"BN{level}"
                else:
                    names[key] = f"{level}L·{ln}" + ("abcdefghijklmnopqrstuvwxyz"[k] if len(keys) > 1 else "")
        # four loops: an irreducible topology whose shipped masters include
        # exactly one plain FMFT master PRk (besides its dotted PRkd) takes its name
        if level == 4:
            pure = defaultdict(set)
            for nm, path, entry in masters_of(root, level):
                expr = (entry.get("exact_basis_expression") or "").strip()
                m = re.fullmatch(r"(?:vakint::\{\}::)?(PR\d+)", expr)
                if m and path.exists():
                    n, edges = parse_dot(path)
                    pure[topology_key(n, [(a, b) for a, b, _ in edges])].add(m[1])
            for key, syms in pure.items():
                if key in nodes and key[0] == "graph" and not nodes[key]["family_names"] and len(syms) == 1:
                    names[key] = next(iter(syms))
        for key in nodes:
            if key[0] == "product":
                fn = []
                for f in key[1]:
                    fl = loops(*from_certificate(f))
                    fn.append(lower_names.get(f) or f"{fl}L·{len(from_certificate(f)[1])}?")
                names[key] = product_name(fn)
                info[key]["factors"] = fn
        # shipped masters (terminals) per topology, with their basis symbols
        masters, symbols, unmatched = defaultdict(list), defaultdict(set), []
        for nm, path, entry in masters_of(root, level):
            if not path.exists():
                unmatched.append(nm)
                continue
            n, edges = parse_dot(path)
            key = topology_key(n, [(a, b) for a, b, _ in edges])
            if key in nodes:
                masters[key].append(nm)
                symbols[key] |= set(BASIS_SYMBOL.findall(entry.get("exact_basis_expression") or ""))
            else:
                unmatched.append(nm)
        # rule programs that own a topology
        owner_bytes, owner_count = defaultdict(int), Counter()
        if level == 5:
            pgraphs = dict(parents)
            for o in owners_5l(root):
                if o["parent_stem"] not in pgraphs or o["deleted"] is None:
                    continue
                pdot = parse_dot(root / "5L/graphs" / f"{o['parent_stem']}.dot")[1]
                labels = [int(at.get("parent_slot", -1)) for _, _, at in pdot]
                cur_n, cur = pgraphs[o["parent_stem"]]
                cur = [(min(a, b), max(a, b)) for a, b in cur]
                ok = True
                for s in o["deleted"]:
                    e = labels.index(s)
                    if cur[e][0] == cur[e][1]:
                        ok = False
                        break
                    cur_n, cur = contract(cur_n, cur, e)
                    labels.pop(e)
                if ok:
                    key = topology_key(cur_n, cur)
                    owner_bytes[key] += o["bytes"]
                    owner_count[key] += 1
        else:
            programs = rule_programs(root, level)
            for key in nodes:
                for label in nodes[key]["parent_of"]:
                    if label in programs:
                        owner_bytes[key] += programs[label]
                        owner_count[key] += 1
        # layout and assembly
        ids = {key: certificate_key(key) for key in nodes}
        layers = defaultdict(list)
        for key in nodes:
            layers[info[key]["lines"]].append(ids[key])
        by_child, by_parent = defaultdict(list), defaultdict(list)
        for (p, ch), k in arrows.items():
            by_child[ids[ch]].append(ids[p])
            by_parent[ids[p]].append(ids[ch])
        order, x, span, width = layout(layers, by_child, by_parent)
        top = max(layers)
        node_list = []
        for key in nodes:
            d = info[key]
            node_list.append({
                "id": ids[key],
                "name": names[key],
                "lines": d["lines"],
                "vertices": d["vertices"],
                "product": d["product"],
                "factors": d.get("factors", []),
                "parent": bool(nodes[key]["parent_of"]),
                "families": sorted(set(nodes[key]["family_names"])),
                "masters": sorted(masters.get(key, []), key=lambda s: (len(s), s)),
                "basis_symbols": sorted(symbols.get(key, set()), key=lambda s: (len(s), s)),
                "owner_programs": owner_count.get(key, 0),
                "owner_bytes": owner_bytes.get(key, 0),
                "layer": top - d["lines"],
                "x": round(x[ids[key]], 4),
                "order": order[d["lines"]].index(ids[key]),
            })
        node_list.sort(key=lambda v: (v["layer"], v["order"]))
        arrow_list = sorted(
            ({"from": ids[p], "to": ids[ch], "lines": k} for (p, ch), k in arrows.items()),
            key=lambda e: (e["from"], e["to"]),
        )
        result["loops"][str(level)] = {
            "loop_count": level,
            "parents": [p for p, _ in parents],
            "nodes": node_list,
            "arrows": arrow_list,
            "span": round(span, 4),
            "crossings": layout.crossings,
            "max_layer_width": width,
            "layers": len(layers),
            "stats": {
                "topologies": len(nodes),
                "irreducible": sum(1 for key in nodes if key[0] == "graph"),
                "products": sum(1 for key in nodes if key[0] == "product"),
                "transitions": len(arrows),
                "removable_lines": sum(arrows.values()),
                "masters": sum(len(v) for v in masters.values()),
                "master_topologies": len(masters),
                "unmatched_masters": unmatched,
                "owner_programs": sum(owner_count.values()),
                "owner_bytes": sum(owner_bytes.values()),
            },
        }
        for key in nodes:
            if key[0] == "graph":
                lower_names.setdefault(key[1], names[key])
    return result


def rule_programs(root: Path, level: int):
    """Family graph stem -> bytes of the RustRed rule program it owns (1L-4L)."""
    manifest = json.loads((root / f"{level}L/manifest.json").read_text())
    out = {}
    ibps = root / f"{level}L/ibps"
    if level <= 3:
        size = (ibps / "rules.rrbin").stat().st_size if (ibps / "rules.rrbin").exists() else 0
        for f in (root / f"{level}L/graphs").glob("*.dot"):
            out[f.stem] = size
        return out
    alias = {}
    for m in manifest["masters"]:
        alias.setdefault(m["family_sha256"], (m.get("alias") or "").replace("4L_", ""))
    for sha, a in alias.items():
        out[sha] = sum(p.stat().st_size for p in ibps.glob(f"{a}.*") if p.is_file())
    return out


def main(argv):
    if len(argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    root, out = Path(argv[1]), Path(argv[2])
    result = extract(root)
    manifest = root / "manifest.json"
    result["data_root"] = str(root)
    result["root_manifest_sha256"] = hashlib.sha256(manifest.read_bytes()).hexdigest() if manifest.exists() else None
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=1, ensure_ascii=False) + "\n")
    for level, g in result["loops"].items():
        s = g["stats"]
        print(f"{level}L: {s['topologies']} topologies ({s['irreducible']} irreducible, {s['products']} products), "
              f"{s['transitions']} arrows, {s['masters']} masters on {s['master_topologies']} topologies, "
              f"unmatched {len(s['unmatched_masters'])}, owner programs {s['owner_programs']}, layers {g['layers']}, widest {g['max_layer_width']}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
