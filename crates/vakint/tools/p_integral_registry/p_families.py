"""Enumerate L-loop massless p-integral top families (cut one edge per orbit of the
3-edge-connected cubic (L+1)-loop vacuum graphs) and map them to Forcer topology names."""
import itertools, re, sys
from collections import defaultdict

def cubic_graphs(n):
    """All simple cubic graphs on n vertices up to isomorphism (brute force, n <= 8)."""
    verts = list(range(n))
    pairs = list(itertools.combinations(verts, 2))
    seen, out = {}, []
    def rec(i, deg, edges):
        if len(edges) == 3 * n // 2:
            if all(d == 3 for d in deg):
                inv = invariant(n, edges)
                if inv not in seen:
                    seen[inv] = True; out.append(list(edges))

            return
        # pick lowest-index vertex with deg < 3 and give it its next neighbour
        v = next(x for x in verts if deg[x] < 3)
        for w in range(v + 1, n):
            if deg[w] < 3 and (v, w) not in edges and (not edges or (v, w) > edges[-1] or edges[-1][0] != v):
                if edges and edges[-1][0] == v and (v, w) <= edges[-1]:
                    continue
                deg[v] += 1; deg[w] += 1; edges.append((v, w))
                rec(i, deg, edges)
                edges.pop(); deg[v] -= 1; deg[w] -= 1
    rec(0, [0] * n, [])
    return out

def invariant(n, edges):
    adj = defaultdict(set)
    for a, b in edges: adj[a].add(b); adj[b].add(a)
    def dist(s):
        d = {s: 0}; fr = [s]
        while fr:
            nf = []
            for v in fr:
                for w in adj[v]:
                    if w not in d: d[w] = d[v] + 1; nf.append(w)
            fr = nf
        return tuple(sorted(d.values()))
    tri = [sum(1 for a, b in itertools.combinations(adj[v], 2) if b in adj[a]) for v in range(n)]
    return tuple(sorted((tri[v], dist(v)) for v in range(n)))

def canon(n, edges, marked):
    """Canonical form; `marked` vertices (unordered set) must map onto positions 0..k-1."""
    best = None
    k = len(marked)
    rest = [v for v in range(n) if v not in marked]
    for head in itertools.permutations(marked):
        for tail in itertools.permutations(rest):
            perm = {v: i for i, v in enumerate(list(head) + list(tail))}
            form = tuple(sorted(tuple(sorted((perm[a], perm[b]))) for a, b in edges))
            if best is None or form < best:
                best = form
    return best

def connected(n, edges):
    adj = defaultdict(set)
    for a, b in edges: adj[a].add(b); adj[b].add(a)
    seen, stack = {0}, [0]
    while stack:
        v = stack.pop()
        for w in adj[v]:
            if w not in seen: seen.add(w); stack.append(w)
    return len(seen) == n

def three_edge_connected(n, edges):
    for k in (1, 2):
        for cut in itertools.combinations(range(len(edges)), k):
            rest = [e for i, e in enumerate(edges) if i not in cut]
            if not connected(n, rest): return False
    return True

def edge_orbits(n, edges):
    es = [tuple(sorted(e)) for e in edges]
    eset = sorted(es)
    autos = []
    for p in itertools.permutations(range(n)):
        if sorted(tuple(sorted((p[a], p[b]))) for a, b in es) == eset:
            autos.append(p)
    orbits, assigned = [], set()
    for e in es:
        if e in assigned: continue
        orb = {tuple(sorted((p[e[0]], p[e[1]]))) for p in autos}
        assigned |= orb; orbits.append(sorted(orb))
    return orbits, len(autos)

def forcer_topologies(path, lines):
    text = open(path).read()
    out = {}
    for e in re.split(r'\n\+', '\n' + text.strip().lstrip('+')):
        m = re.search(r'TOPO\((\w+)\)', e)
        if not m or 'ex(' in e: continue
        vxs = re.findall(r'vx\(([^)]*)\)', e)
        ends = defaultdict(list); ext = []
        for i, v in enumerate(vxs):
            for a in v.split(','):
                a = a.strip()
                if a.lstrip('-') == 'Q': ext.append(i)
                else: ends[a.lstrip('-')].append(i)
        if len(ends) != lines: continue
        edges = [tuple(v) for v in ends.values()]
        if any(len(v) != 2 for v in edges): continue
        out[m.group(1)] = (len(vxs), edges, tuple(ext))
    return out

if __name__ == '__main__':
    L = int(sys.argv[1]); n = 2 * L
    parents = [g for g in cubic_graphs(n) if connected(n, g) and three_edge_connected(n, g)]
    assert len({canon(n, g, ()) for g in parents}) == len(parents)
    forcer = forcer_topologies(sys.argv[2], 3 * L - 1) if len(sys.argv) > 2 else {}
    fcanon = defaultdict(list)
    for name, (nv, edges, ext) in forcer.items():
        if nv == n:
            fcanon[canon(n, edges, ext)].append(name)
    print(f"{L}L: {len(parents)} 3-edge-connected cubic vacuum parents on {n} vertices")
    fams = {}
    for pi, g in enumerate(parents):
        orbits, aut = edge_orbits(n, g)
        print(f" parent {pi}: edges={g} |Aut|={aut} orbit sizes={[len(o) for o in orbits]}")
        for o in orbits:
            e = o[0]
            rest = [x for x in g if tuple(sorted(x)) != e]
            c = canon(n, rest, e)
            fams.setdefault(c, (pi, e, len(o)))
    print(f" families: {len(fams)}")
    for c, (pi, e, size) in fams.items():
        print(f"  parent {pi} cut {e} (orbit {size}): forcer={fcanon.get(c, ['-'])}")
