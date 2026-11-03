"""Usage: python3 -I registry.py <this directory> <forcer/notation.h> <out.json>
(prints the canonical routings; the Rust table in src/p_integral.rs is rendered from out.json).
"""
"""Canonical 1L-4L p-integral registry: q enters at vertex 1, leaves at vertex 2,
loop lines first (line i carries k(i)), ISPs completing the scalar-product basis."""
import itertools, sys, json
from fractions import Fraction
from collections import defaultdict
sys.path.insert(0, sys.argv[1])
import p_families as pf

FORCER = {}

def families(L, notation):
    if L == 1:
        return [("P1L", None, [(0, 1), (0, 1)], (0, 1), 2)]
    n = 2 * L
    parents = [g for g in pf.cubic_graphs(n) if pf.connected(n, g) and pf.three_edge_connected(n, g)]
    forcer = pf.forcer_topologies(notation, 3 * L - 1) if notation else {}
    fc = defaultdict(list)
    for name, (nv, edges, ext) in forcer.items():
        if nv == n: fc[pf.canon(n, edges, ext)].append(name)
    out, seen = [], set()
    for pi, g in enumerate(parents):
        orbits, aut = pf.edge_orbits(n, g)
        for o in orbits:
            e = o[0]
            rest = [x for x in g if tuple(sorted(x)) != e]
            c = pf.canon(n, rest, e)
            if c in seen: continue
            seen.add(c)
            out.append((fc.get(c, [None])[0], pi, rest, e, n))
    return out

def route(L, edges, ext, n):
    a, b = ext
    order = [a, b] + [v for v in range(n) if v not in ext]
    # BFS order of the remaining vertices from a for readable labels
    adj = defaultdict(list)
    for i, (u, v) in enumerate(edges): adj[u].append((v, i)); adj[v].append((u, i))
    lab = {a: 1, b: 2}; nxt = 3; fr = [a, b]
    while fr:
        nf = []
        for v in fr:
            for w, _ in sorted(adj[v]):
                if w not in lab: lab[w] = nxt; nxt += 1; nf.append(w)
        fr = nf
    E = [(lab[u], lab[v]) for u, v in edges]
    V = n
    # spanning tree by BFS from vertex 1 over relabelled graph
    adj = defaultdict(list)
    for i, (u, v) in enumerate(E): adj[u].append((v, i)); adj[v].append((u, i))
    parent = {1: None}; tree = set(); fr = [1]
    while fr:
        nf = []
        for v in fr:
            for w, i in sorted(adj[v]):
                if w not in parent: parent[w] = (v, i); tree.add(i); nf.append(w)
        fr = nf
    chords = [i for i in range(len(E)) if i not in tree]
    assert len(chords) == L
    # orient chords low->high; momentum vectors over (k1..kL, q)
    mom = {}
    for j, i in enumerate(chords):
        u, v = E[i]
        E[i] = (min(u, v), max(u, v))
        vec = [0] * (L + 1); vec[j] = 1; mom[i] = vec
    # tree momenta from conservation: process leaves upward
    # balance(v) = sum_out - sum_in must equal +q at 1, -q at 2, 0 else
    target = {v: [0] * (L + 1) for v in range(1, V + 1)}
    target[1][L] = 1; target[2][L] = -1
    depth = {1: 0}
    for v in sorted(parent, key=lambda x: 0):
        pass
    order_v = sorted(parent, key=lambda v: -bfs_depth(v, parent))
    for v in order_v:
        if parent[v] is None: continue
        pv, i = parent[v]
        # current balance at v from known lines
        bal = [0] * (L + 1)
        for w, k in adj[v]:
            if k == i or k not in mom: continue
            s = 1 if E[k][0] == v else -1
            if E[k][0] == E[k][1]: s = 0
            bal = [x + s * y for x, y in zip(bal, mom[k])]
        need = [t - x for t, x in zip(target[v], bal)]
        # line i oriented v -> pv carries `need` (outgoing from v)
        if any(need) and next(x for x in need if x) < 0:
            E[i] = (pv, v); mom[i] = [-x for x in need]
        else:
            E[i] = (v, pv); mom[i] = need
    # verify conservation
    for v in range(1, V + 1):
        bal = [0] * (L + 1)
        for k, (x, y) in enumerate(E):
            if x == v: bal = [p + q for p, q in zip(bal, mom[k])]
            if y == v: bal = [p - q for p, q in zip(bal, mom[k])]
        assert bal == target[v], (v, bal)
    lines = [(E[i], mom[i]) for i in chords] + sorted(((E[i], mom[i]) for i in tree), key=lambda t: (sum(map(abs, t[1])), [-x for x in t[1]]))
    return lines

def bfs_depth(v, parent):
    d = 0
    while parent[v] is not None: v = parent[v][0]; d += 1
    return d

def quad(vec, L):
    """coefficients of (sum a_i k_i + b q)^2 on k_i.k_j (i<=j) and k_i.q"""
    out = []
    for i in range(L):
        for j in range(i, L):
            out.append(vec[i] * vec[j] * (1 if i == j else 2))
    for i in range(L):
        out.append(2 * vec[i] * vec[L])
    return out

def rank(rows):
    rows = [[Fraction(x) for x in r] for r in rows]; r = 0
    if not rows: return 0
    cols = len(rows[0])
    for c in range(cols):
        p = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c] / rows[r][c]
                rows[i] = [x - f * y for x, y in zip(rows[i], rows[r])]
        r += 1
    return r

def isps(L, lines):
    rows = [quad(m, L) for _, m in lines]
    need = L * (L + 1) // 2 + L
    out = []
    cands = []
    for coeffs in itertools.product([0, 1, -1], repeat=L):
        for b in (0, -1):
            v = list(coeffs) + [b]
            if not any(coeffs): continue
            if next(x for x in coeffs if x) < 0: continue
            cands.append(v)
    cands.sort(key=lambda v: (sum(map(abs, v)), [abs(x) for x in v][::-1], v[::-1]))
    r = rank(rows)
    for v in cands:
        if r == need: break
        if rank(rows + [quad(v, L)]) > r:
            rows.append(quad(v, L)); out.append(v); r += 1
    assert r == need, (r, need)
    return out

def fmt(vec, L):
    s = ""
    for i, c in enumerate(vec):
        name = f"k({i+1})" if i < L else "p(1)"
        if c == 0: continue
        s += ("+" if c > 0 and s else "") + ("-" if c < 0 else "") + (str(abs(c)) + "*" if abs(c) != 1 else "") + name
    return s

notation = sys.argv[2]
names3 = {}
result = []
for L in (1, 2, 3, 4):
    for fam in families(L, notation):
        fname, pi, edges, ext, n = fam
        if L == 1:
            lines = [((1, 2), [1, 0]), ((2, 1), [1, -1])]
        else:
            lines = route(L, edges, ext, n)
        name = {1: "P1L", 2: "P2L"}.get(L)
        if L == 3:
            # prism cut rung -> LA, prism cut triangle edge -> BE, K33 -> NO
            name = None
        result.append(dict(L=L, forcer=fname, parent=pi, cut=ext, lines=lines, isps=isps(L, lines) if L > 1 else []))
for r in result:
    L = r["L"]
    print(f"L={L} parent={r['parent']} cut={r['cut']} forcer={r['forcer']}")
    for i, ((a, b), m) in enumerate(r["lines"]):
        print(f"   prop({i+1},edge({a},{b}),{fmt(m, L)},0,pow({i+1}))")
    for m in r["isps"]:
        print(f"   ISP ({fmt(m, L)})^2")
json.dump(result, open(sys.argv[3], "w"), indent=1)
