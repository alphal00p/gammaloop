"""Canonization on the folded node encoding agrees with brute-force half-edge symmetries."""

import random
import unittest
from itertools import permutations

from symbolica import S
from symbolica.community import graph as g

UNDIRECTED = g.Orientation.Undirected
DEFAULT = g.Orientation.Default
REVERSED = g.Orientation.Reversed


def key(data):
    return "0" if data is None else str(data)


def leg(node, data=0, *, incoming=False, orientation=DEFAULT, port=None):
    end = g.sink if incoming else g.source
    return g.edge(end(node, data=port), data=data, orientation=orientation)


def link(a, b, data=0, *, orientation=UNDIRECTED, ports=(None, None)):
    return g.edge(
        g.source(a, data=ports[0]),
        second=g.sink(b, data=ports[1]),
        data=data,
        orientation=orientation,
    )


def edges(graph):
    """Edges as (edge key, directed, ends) with ends (node, half-edge key) tail first."""
    out = []
    for edge in graph.edges():
        ends = [
            (h.node.index, key(h.data))
            for h in (edge.source, edge.sink)
            if h is not None
        ]
        if edge.source is None or edge.orientation == REVERSED:
            ends.reverse()
        directed = edge.orientation != UNDIRECTED
        if len(ends) == 1 and directed:
            # A dangling half-edge is directed into or out of its node.
            outgoing = (edge.source is not None) == (edge.orientation == DEFAULT)
            out.append((key(edge.data), "out" if outgoing else "in", ends))
        else:
            out.append((key(edge.data), directed, ends))
    return out


def brute_force_count(a, b):
    """Count the half-edge bijections from `a` to `b` that respect nodes and keys."""
    nodes_a = [key(n.data) for n in a.nodes()]
    nodes_b = [key(n.data) for n in b.nodes()]
    edges_a, edges_b = edges(a), edges(b)
    if sorted(nodes_a) != sorted(nodes_b) or len(edges_a) != len(edges_b):
        return 0

    def count(sigma, i, used):
        if i == len(edges_a):
            return 1
        data, directed, ends = edges_a[i]
        image = [(sigma[n], k) for n, k in ends]
        found = 0
        for j, (data_b, directed_b, ends_b) in enumerate(edges_b):
            if j in used or (data, directed) != (data_b, directed_b):
                continue
            maps = [image]
            if directed is False and len(image) == 2:
                maps.append(image[::-1])
            found += sum(m == ends_b for m in maps) * count(sigma, i + 1, used | {j})
        return found

    return sum(
        count(sigma, 0, frozenset())
        for sigma in permutations(range(len(nodes_b)))
        if all(nodes_a[i] == nodes_b[sigma[i]] for i in range(len(nodes_a)))
    )


def shape(graph):
    """Everything canonization promises to fix: element order, keys, roles and incidence."""
    return (
        [key(n.data) for n in graph.nodes()],
        [
            (
                key(e.data),
                e.orientation,
                None if e.source is None else (e.source.node.index, key(e.source.data)),
                None if e.sink is None else (e.sink.node.index, key(e.sink.data)),
            )
            for e in graph.edges()
        ],
        [(h.node.index, h.edge.index, h.flow, key(h.data)) for h in graph.half_edges()],
    )


def cells(orbit):
    return sorted(sorted(i for i, o in enumerate(orbit) if o == r) for r in set(orbit))


def random_graph(rng, n_nodes, n_edges, n_legs):
    nodes = [g.node(str(i), data=rng.randrange(2)) for i in range(n_nodes)]
    items = []
    for _ in range(n_edges):
        a, b = rng.randrange(n_nodes), rng.randrange(n_nodes)
        items.append(
            link(
                str(a),
                str(b),
                rng.randrange(2),
                orientation=rng.choice([UNDIRECTED, UNDIRECTED, DEFAULT, REVERSED]),
                ports=(rng.randrange(2), rng.randrange(2)),
            )
        )
    for _ in range(n_legs):
        items.append(
            leg(
                str(rng.randrange(n_nodes)),
                rng.randrange(2),
                incoming=rng.random() < 0.5,
                orientation=rng.choice([UNDIRECTED, DEFAULT, REVERSED]),
                port=rng.randrange(2),
            )
        )
    return g.build(*nodes, *items)


def relabelled(graph, rng):
    """The same graph with shuffled nodes, edges and stored undirected source/sink roles."""
    order = list(range(graph.n_nodes))
    rng.shuffle(order)
    items = [g.node(str(i), data=graph.node(i).data) for i in order]
    shuffled = list(graph.edges())
    rng.shuffle(shuffled)
    for edge in shuffled:
        ends = [edge.source, edge.sink]
        orientation = edge.orientation
        if None in ends:
            (h,) = [h for h in ends if h is not None]
            flow_source = edge.source is not None
            if orientation == UNDIRECTED:
                flow_source = rng.random() < 0.5
            elif rng.random() < 0.5:
                flow_source = not flow_source
                orientation = REVERSED if orientation == DEFAULT else DEFAULT
            end = g.source if flow_source else g.sink
            items.append(
                g.edge(
                    end(str(h.node.index), data=h.data),
                    data=edge.data,
                    orientation=orientation,
                )
            )
            continue
        if orientation == UNDIRECTED and rng.random() < 0.5:
            ends.reverse()
        elif orientation != UNDIRECTED and rng.random() < 0.5:
            ends.reverse()
            orientation = REVERSED if orientation == DEFAULT else DEFAULT
        items.append(
            g.edge(
                g.source(str(ends[0].node.index), data=ends[0].data),
                second=g.sink(str(ends[1].node.index), data=ends[1].data),
                data=edge.data,
                orientation=orientation,
            )
        )
    return g.build(*items)


class CanonizeTests(unittest.TestCase):
    def assert_group(self, graph, expected):
        self.assertEqual(graph.canonize()[2], expected)
        self.assertEqual(brute_force_count(graph, graph), expected)

    def test_identical_legs_at_one_node_permute(self):
        self.assert_group(g.build(g.node("a"), *(leg("a") for _ in range(4))), 24)
        self.assert_group(
            g.build(g.node("a"), leg("a", 1), leg("a", 1), leg("a", 2)), 2
        )
        # The stored flow of an undirected leg is no invariant.
        self.assert_group(
            g.build(
                g.node("a"),
                leg("a", orientation=UNDIRECTED),
                leg("a", incoming=True, orientation=UNDIRECTED),
            ),
            2,
        )
        # Directed legs differ by their physical direction only.
        self.assert_group(g.build(g.node("a"), leg("a"), leg("a", incoming=True)), 1)
        self.assert_group(
            g.build(
                g.node("a"), leg("a"), leg("a", incoming=True, orientation=REVERSED)
            ),
            2,
        )

    def test_distinct_legs_and_symmetric_nodes(self):
        self.assert_group(g.build(g.node("a"), *(leg("a", i) for i in range(4))), 1)
        self.assert_group(g.build(g.node("a"), leg("a", port=1), leg("a", port=2)), 1)
        both = g.build(
            g.node("a"), g.node("b"), link("a", "b"), leg("a", 1), leg("b", 1)
        )
        self.assert_group(both, 2)
        self.assert_group(
            g.build(g.node("a"), g.node("b"), link("a", "b"), leg("a", 1), leg("b", 2)),
            1,
        )

    def test_undirected_self_loop_with_distinct_ends_loses_its_flip(self):
        self.assert_group(g.build(g.node("a"), link("a", "a")), 2)
        self.assert_group(g.build(g.node("a"), link("a", "a", ports=(1, 2))), 1)
        self.assert_group(
            g.build(
                g.node("a"), link("a", "a", ports=(1, 1)), link("a", "a", ports=(1, 1))
            ),
            8,
        )
        self.assert_group(
            g.build(
                g.node("a"), link("a", "a", ports=(1, 2)), link("a", "a", ports=(2, 1))
            ),
            2,
        )

    def test_parallel_bundles(self):
        self.assert_group(
            g.build(g.node("a"), g.node("b"), *(link("a", "b") for _ in range(3))), 12
        )
        self.assert_group(
            g.build(
                g.node("a", data=1),
                g.node("b", data=2),
                *(link("a", "b") for _ in range(3)),
            ),
            6,
        )
        self.assert_group(
            g.build(
                g.node("a"),
                g.node("b"),
                link("a", "b", 1),
                link("a", "b", 1),
                link("a", "b", 2),
            ),
            4,
        )
        self.assert_group(
            g.build(
                g.node("a"),
                g.node("b"),
                link("a", "b", orientation=DEFAULT),
                link("a", "b", orientation=DEFAULT),
                link("b", "a", orientation=DEFAULT),
            ),
            2,
        )
        # Distinct end keys orient every edge of the bundle the same way.
        self.assert_group(
            g.build(
                g.node("a"),
                g.node("b"),
                link("a", "b", ports=(1, 2)),
                link("b", "a", ports=(2, 1)),
            ),
            2,
        )

    def test_directed_and_undirected_edges_differ(self):
        undirected = g.build(g.node("a"), g.node("b"), link("a", "b", ports=(1, 2)))
        directed = g.build(
            g.node("a"), g.node("b"), link("a", "b", ports=(1, 2), orientation=DEFAULT)
        )
        reversed_ = g.build(
            g.node("a"), g.node("b"), link("b", "a", ports=(2, 1), orientation=REVERSED)
        )
        self.assertFalse(undirected.is_isomorphic(directed))
        self.assertTrue(directed.is_isomorphic(reversed_))
        self.assertEqual(shape(directed.canonize()[0]), shape(reversed_.canonize()[0]))
        self.assertNotEqual(
            shape(directed.canonize()[0]), shape(undirected.canonize()[0])
        )
        self.assert_group(
            g.build(g.node("a"), g.node("b"), link("a", "b", orientation=DEFAULT)), 1
        )
        self.assert_group(g.build(g.node("a"), g.node("b"), link("a", "b")), 2)

    def test_random_graphs_match_brute_force(self):
        rng = random.Random(7)
        graphs = [
            random_graph(
                rng, rng.randrange(1, 5), rng.randrange(0, 5), rng.randrange(0, 4)
            )
            for _ in range(150)
        ]
        for graph in graphs:
            with self.subTest(graph=shape(graph)):
                self.assertEqual(graph.canonize()[2], brute_force_count(graph, graph))
        outcomes = set()
        for graph in graphs:
            # Relabelled copies, most with one half-edge key toggled: often a near miss.
            other = relabelled(graph, rng)
            if other.n_half_edges and rng.random() < 0.7:
                half_edge = rng.choice(other.half_edges())
                half_edge.data = 1 - int(key(half_edge.data))
            expected = brute_force_count(graph, other) > 0
            outcomes.add(expected)
            self.assertEqual(graph.is_isomorphic(other), expected)
            same = shape(graph.canonize()[0]) == shape(other.canonize()[0])
            self.assertEqual(same, expected)
        self.assertEqual(outcomes, {True, False})

    def test_relabelling_gives_identical_canonical_graphs(self):
        rng = random.Random(11)
        for _ in range(60):
            graph = random_graph(
                rng, rng.randrange(1, 6), rng.randrange(0, 7), rng.randrange(0, 4)
            )
            canonical, mapping, size, orbit = graph.canonize()
            for old, new in enumerate(mapping):
                self.assertEqual(
                    key(canonical.node(new).data), key(graph.node(old).data)
                )
            self.assertEqual(brute_force_count(graph, canonical), size)
            for _ in range(4):
                other = relabelled(graph, rng)
                self.assertTrue(graph.is_isomorphic(other))
                result = other.canonize()
                self.assertEqual(shape(result[0]), shape(canonical))
                self.assertEqual(result[2], size)
                self.assertEqual(cells(result[3]), cells(orbit))

    def test_generated_graphs_keep_their_symmetry_factors(self):
        gluon = g.EdgeSignature(S("g"))
        quark = g.EdgeSignature(S("q"), direction=False)
        rng = random.Random(3)
        generated = g.Graph.generate(
            [(1, gluon), (2, gluon)],
            [[gluon] * 3, [quark.flip(), quark, gluon], [gluon] * 4],
            max_loops=3,
            max_bridges=0,
            allow_self_loops=True,
        )
        self.assertEqual(len(generated), 210)
        for graph, size in generated:
            canonical = graph.canonize()[0]
            self.assertEqual(graph.canonize()[2], size)
            other = relabelled(graph, rng)
            self.assertEqual(other.canonize()[2], size)
            self.assertEqual(shape(other.canonize()[0]), shape(canonical))

    def test_dangling_legs_count_like_identical_external_vertices(self):
        # Graphica sees external legs as degree-one vertices; with one shared colour their
        # permutations are graph automorphisms, which the folded legs must reproduce.
        phi = g.EdgeSignature(S("phi"))
        generated = g.Graph.generate(
            [(i, phi) for i in range(1, 5)], [[phi] * 3, [phi] * 4], max_loops=1
        )
        for graph, _ in generated:
            degree = [0] * graph.n_nodes
            for h in graph.half_edges():
                degree[h.node.index] += 1
            external = [i for i, d in enumerate(degree) if d == 1]
            vertices = g.build(
                *(
                    g.node(str(i), data=0 if i in external else 1)
                    for i in range(graph.n_nodes)
                ),
                *(
                    link(str(e.source.node.index), str(e.sink.node.index), e.data)
                    for e in graph.edges()
                ),
            )
            legs = g.build(
                *(
                    g.node(str(i), data=1)
                    for i in range(graph.n_nodes)
                    if i not in external
                ),
                *(
                    link(str(a), str(b), e.data)
                    if external.count(a) + external.count(b) == 0
                    else leg(
                        str(b if a in external else a), e.data, orientation=UNDIRECTED
                    )
                    for e in graph.edges()
                    for a, b in [(e.source.node.index, e.sink.node.index)]
                ),
            )
            self.assertEqual(legs.canonize()[2], vertices.canonize()[2])

    def test_edge_canonicalization_sorts_without_moving_vertices(self):
        graph = g.build(
            g.node("a", data=object()),
            g.node("b", data=object()),
            link("b", "a", 2),
            leg("b", 1),
            link("a", "b", 1),
            link("a", "a", 3, ports=(2, 1)),
        )
        nodes = [n.data for n in graph.nodes()]
        graph.canonize_edges()
        self.assertEqual([n.data for n in graph.nodes()], nodes)
        self.assertEqual([e.data for e in graph.edges()], [3, 1, 2, 1])
        swapped = g.build(
            g.node("a", data=object()),
            g.node("b", data=object()),
            link("a", "a", 3, ports=(1, 2)),
            link("b", "a", 1),
            leg("b", 1),
            link("a", "b", 2),
        )
        swapped.canonize_edges()
        self.assertEqual([e.data for e in swapped.edges()], [3, 1, 2, 1])
        self.assertEqual(
            [h.edge.index for h in swapped.half_edges()], [0, 0, 1, 1, 2, 2, 3]
        )


if __name__ == "__main__":
    unittest.main()
