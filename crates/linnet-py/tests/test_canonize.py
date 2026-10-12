"""Native canonization agrees with brute-force half-edge symmetries and ranks keys."""

import math
import numbers
import random
import unittest
from decimal import Decimal
from fractions import Fraction
from itertools import permutations

from linnet_under_test import COMMUNITY, community_only, field
from linnet_under_test import graph as g

UNDIRECTED = g.Orientation.Undirected
DEFAULT = g.Orientation.Default
REVERSED = g.Orientation.Reversed


def key(data):
    return "0" if data is None else str(data)


def nodes(*names):
    return [g.Node(name) for name in names]


def leg(node, data=0, *, incoming=False, orientation=DEFAULT, port=None):
    end = g.HalfEdge(node, data=port)
    if incoming:
        return g.Edge(sink=end, data=data, orientation=orientation)
    return g.Edge(end, data=data, orientation=orientation)


def link(a, b, data=0, *, orientation=UNDIRECTED, ports=(None, None)):
    return g.Edge(
        g.HalfEdge(a, data=ports[0]),
        g.HalfEdge(b, data=ports[1]),
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
    return g.Graph(
        nodes=[g.Node(str(i), data=rng.randrange(2)) for i in range(n_nodes)],
        edges=items,
    )


def relabelled(graph, rng):
    """The same graph with shuffled nodes, edges and stored undirected source/sink roles."""
    order = list(range(graph.n_nodes))
    rng.shuffle(order)
    shuffled = list(graph.edges())
    rng.shuffle(shuffled)
    items = []
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
            end = g.HalfEdge(str(h.node.index), data=h.data)
            items.append(
                g.Edge(end, data=edge.data, orientation=orientation)
                if flow_source
                else g.Edge(sink=end, data=edge.data, orientation=orientation)
            )
            continue
        if orientation == UNDIRECTED and rng.random() < 0.5:
            ends.reverse()
        elif orientation != UNDIRECTED and rng.random() < 0.5:
            ends.reverse()
            orientation = REVERSED if orientation == DEFAULT else DEFAULT
        items.append(
            g.Edge(
                g.HalfEdge(str(ends[0].node.index), data=ends[0].data),
                g.HalfEdge(str(ends[1].node.index), data=ends[1].data),
                data=edge.data,
                orientation=orientation,
            )
        )
    return g.Graph(
        nodes=[g.Node(str(i), data=graph.node(i).data) for i in order], edges=items
    )


class CanonizeTests(unittest.TestCase):
    def assert_group(self, graph, expected):
        self.assertEqual(graph.canonize()[2], expected)
        self.assertEqual(brute_force_count(graph, graph), expected)

    def test_identical_legs_at_one_node_permute(self):
        self.assert_group(
            g.Graph(nodes=nodes("a"), edges=[leg("a") for _ in range(4)]), 24
        )
        self.assert_group(
            g.Graph(nodes=nodes("a"), edges=[leg("a", 1), leg("a", 1), leg("a", 2)]), 2
        )
        # The stored flow of an undirected leg is no invariant.
        self.assert_group(
            g.Graph(
                nodes=nodes("a"),
                edges=[
                    leg("a", orientation=UNDIRECTED),
                    leg("a", incoming=True, orientation=UNDIRECTED),
                ],
            ),
            2,
        )
        # Directed legs differ by their physical direction only.
        self.assert_group(
            g.Graph(nodes=nodes("a"), edges=[leg("a"), leg("a", incoming=True)]), 1
        )
        self.assert_group(
            g.Graph(
                nodes=nodes("a"),
                edges=[leg("a"), leg("a", incoming=True, orientation=REVERSED)],
            ),
            2,
        )

    def test_distinct_legs_and_symmetric_nodes(self):
        self.assert_group(
            g.Graph(nodes=nodes("a"), edges=[leg("a", i) for i in range(4)]), 1
        )
        self.assert_group(
            g.Graph(nodes=nodes("a"), edges=[leg("a", port=1), leg("a", port=2)]), 1
        )
        both = g.Graph(
            nodes=nodes("a", "b"), edges=[link("a", "b"), leg("a", 1), leg("b", 1)]
        )
        self.assert_group(both, 2)
        self.assert_group(
            g.Graph(
                nodes=nodes("a", "b"), edges=[link("a", "b"), leg("a", 1), leg("b", 2)]
            ),
            1,
        )

    def test_undirected_self_loop_with_distinct_ends_loses_its_flip(self):
        self.assert_group(g.Graph(nodes=nodes("a"), edges=[link("a", "a")]), 2)
        self.assert_group(
            g.Graph(nodes=nodes("a"), edges=[link("a", "a", ports=(1, 2))]), 1
        )
        self.assert_group(
            g.Graph(
                nodes=nodes("a"),
                edges=[link("a", "a", ports=(1, 1)), link("a", "a", ports=(1, 1))],
            ),
            8,
        )
        self.assert_group(
            g.Graph(
                nodes=nodes("a"),
                edges=[link("a", "a", ports=(1, 2)), link("a", "a", ports=(2, 1))],
            ),
            2,
        )

    def test_parallel_bundles(self):
        self.assert_group(
            g.Graph(nodes=nodes("a", "b"), edges=[link("a", "b") for _ in range(3)]),
            12,
        )
        self.assert_group(
            g.Graph(
                nodes=[g.Node("a", data=1), g.Node("b", data=2)],
                edges=[link("a", "b") for _ in range(3)],
            ),
            6,
        )
        self.assert_group(
            g.Graph(
                nodes=nodes("a", "b"),
                edges=[link("a", "b", 1), link("a", "b", 1), link("a", "b", 2)],
            ),
            4,
        )
        self.assert_group(
            g.Graph(
                nodes=nodes("a", "b"),
                edges=[
                    link("a", "b", orientation=DEFAULT),
                    link("a", "b", orientation=DEFAULT),
                    link("b", "a", orientation=DEFAULT),
                ],
            ),
            2,
        )
        # Distinct end keys orient every edge of the bundle the same way.
        self.assert_group(
            g.Graph(
                nodes=nodes("a", "b"),
                edges=[link("a", "b", ports=(1, 2)), link("b", "a", ports=(2, 1))],
            ),
            2,
        )

    def test_directed_and_undirected_edges_differ(self):
        undirected = g.Graph(
            nodes=nodes("a", "b"), edges=[link("a", "b", ports=(1, 2))]
        )
        directed = g.Graph(
            nodes=nodes("a", "b"),
            edges=[link("a", "b", ports=(1, 2), orientation=DEFAULT)],
        )
        reversed_ = g.Graph(
            nodes=nodes("a", "b"),
            edges=[link("b", "a", ports=(2, 1), orientation=REVERSED)],
        )
        self.assertFalse(undirected.is_isomorphic(directed))
        self.assertTrue(directed.is_isomorphic(reversed_))
        self.assertEqual(shape(directed.canonize()[0]), shape(reversed_.canonize()[0]))
        self.assertNotEqual(
            shape(directed.canonize()[0]), shape(undirected.canonize()[0])
        )
        self.assert_group(
            g.Graph(nodes=nodes("a", "b"), edges=[link("a", "b", orientation=DEFAULT)]),
            1,
        )
        self.assert_group(g.Graph(nodes=nodes("a", "b"), edges=[link("a", "b")]), 2)

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
            # Orbits are indexed by canonical vertex and named by their smallest one.
            for vertex, first in enumerate(orbit):
                self.assertLessEqual(first, vertex)
                self.assertEqual(orbit[first], first)
            self.assertEqual(brute_force_count(graph, canonical), size)
            for _ in range(4):
                other = relabelled(graph, rng)
                self.assertTrue(graph.is_isomorphic(other))
                result = other.canonize()
                self.assertEqual(shape(result[0]), shape(canonical))
                self.assertEqual(result[2], size)
                self.assertEqual(cells(result[3]), cells(orbit))
                self.assertEqual(result[3], orbit)

    def test_generated_graphs_keep_their_symmetry_factors(self):
        gluon = g.EdgeSignature(field("g"))
        quark = g.EdgeSignature(field("q"), direction=False)
        rng = random.Random(3)
        generated = list(
            g.Theory([[gluon] * 3, [quark.flip(), quark, gluon], [gluon] * 4])
            .query(
                g.Legs([(1, gluon), (2, gluon)]),
                loops=range(4),
                structure=g.Structure(bridges=g.Bridges.one_pi(), self_loops=True),
            )
            .diagrams()
        )
        self.assertEqual(len(generated), 210)
        for diagram in generated:
            graph, size = diagram.graph(), diagram.symmetry_factor()
            canonical = graph.canonize()[0]
            self.assertEqual(graph.canonize()[2], size)
            other = relabelled(graph, rng)
            self.assertEqual(other.canonize()[2], size)
            self.assertEqual(shape(other.canonize()[0]), shape(canonical))

    def test_dangling_legs_count_like_identical_external_vertices(self):
        # With one shared colour, permutations of dangling legs at a node are
        # automorphisms, exactly as for legs ending in identical degree-one vertices.
        phi = g.EdgeSignature(field("phi"))
        plan = g.Theory([[phi] * 3, [phi] * 4]).query(
            g.Legs([(i, phi) for i in range(1, 5)]), loops=range(2)
        )
        for diagram in plan.diagrams():
            graph = diagram.graph()
            names = [str(i) for i in range(graph.n_nodes)]
            vertices = [g.Node(name, data=1) for name in names]
            links = []
            for e in graph.edges():
                ends = [str(h.node.index) for h in (e.source, e.sink) if h is not None]
                if len(ends) == 1:
                    ends.append(f"leg{len(vertices)}")
                    vertices.append(g.Node(ends[1], data=0))
                links.append(link(*ends, e.data))
            external = g.Graph(nodes=vertices, edges=links)
            dangling = graph.canonize(
                node_key=lambda _: 1, half_edge_key=lambda _: None
            )
            self.assertEqual(dangling[2], external.canonize()[2])

    def test_edge_canonicalization_sorts_without_moving_vertices(self):
        graph = g.Graph(
            nodes=[g.Node("a", data=object()), g.Node("b", data=object())],
            edges=[
                link("b", "a", 2),
                leg("b", 1),
                link("a", "b", 1),
                link("a", "a", 3, ports=(2, 1)),
            ],
        )
        payloads = [n.data for n in graph.nodes()]
        graph.canonize_edges()
        self.assertEqual([n.data for n in graph.nodes()], payloads)
        self.assertEqual([e.data for e in graph.edges()], [3, 1, 2, 1])
        swapped = g.Graph(
            nodes=[g.Node("a", data=object()), g.Node("b", data=object())],
            edges=[
                link("a", "a", 3, ports=(1, 2)),
                link("b", "a", 1),
                leg("b", 1),
                link("a", "b", 2),
            ],
        )
        swapped.canonize_edges()
        self.assertEqual([e.data for e in swapped.edges()], [3, 1, 2, 1])
        self.assertEqual(
            [h.edge.index for h in swapped.half_edges()], [0, 0, 1, 1, 2, 2, 3]
        )


class Opaque:
    """Equal by identity, ordered among other values by type name and repr."""

    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return f"Opaque({self.name})"


class Named:
    """Equal by name, with the default address-based repr."""

    def __init__(self, name):
        self.name = name

    def __eq__(self, other):
        return isinstance(other, Named) and self.name == other.name

    def __hash__(self):
        return hash(self.name)


class Shown(Named):
    """Equal by name, shown as given."""

    def __init__(self, name, shown):
        super().__init__(name)
        self.shown = shown

    def __repr__(self):
        return self.shown


class Cyclic(Named):
    """A real number whose `<` is not a total order."""

    def __lt__(self, other):
        return (self.name * 7 + other.name * 13) % 3 == 0

    def __repr__(self):
        return f"Cyclic({self.name})"


class Unordered(Cyclic):
    def __lt__(self, other):
        raise ArithmeticError("no order")


numbers.Real.register(Cyclic)


class BadHash:
    def __hash__(self):
        raise TypeError("bad hash")


def isolated(*payloads):
    return g.Graph(nodes=[g.Node(data=p) for p in payloads])


class KeyTests(unittest.TestCase):
    def test_mixed_kinds_sort_by_kind_then_value(self):
        late, early = Opaque("b"), Opaque("a")
        payloads = [late, (2, "x"), b"b", "s", 2, None, 1.5, (1, "y"), early]
        expected = [None, 1.5, 2, "s", b"b", (1, "y"), (2, "x"), early, late]
        if COMMUNITY:
            from symbolica import E

            payloads.insert(4, E("x"))
            expected.insert(expected.index(early), E("x"))
        rng = random.Random(5)
        for _ in range(5):
            rng.shuffle(payloads)
            canonical = isolated(*payloads).canonize()[0]
            self.assertEqual(
                [repr(n.data) for n in canonical.nodes()], [repr(p) for p in expected]
            )

    def test_python_equality_identifies_keys(self):
        self.assertTrue(isolated(1, "a").is_isomorphic(isolated("a", 1.0)))
        self.assertTrue(isolated(True, (1, 2)).is_isomorphic(isolated((1.0, 2), 1)))
        self.assertFalse(isolated(1).is_isomorphic(isolated("1")))
        self.assertFalse(isolated(None).is_isomorphic(isolated(0)))
        self.assertFalse(isolated(Opaque("a")).is_isomorphic(isolated(Opaque("a"))))
        shared = Opaque("a")
        self.assertTrue(isolated(shared).is_isomorphic(isolated(shared)))
        # Two identical nodes swap; two distinct keys do not.
        self.assertEqual(isolated(1, 1.0).canonize()[2], 2)
        self.assertEqual(isolated(1, "1").canonize()[2], 1)

    @community_only
    def test_expressions_are_their_own_kind(self):
        from symbolica import E

        self.assertTrue(isolated(E("x+1")).is_isomorphic(isolated(E("1+x"))))
        self.assertFalse(isolated(E("1")).is_isomorphic(isolated(1)))
        self.assertFalse(isolated((E("1"), 2)).is_isomorphic(isolated((1, 2))))
        self.assertEqual(isolated(E("1"), 1).canonize()[2], 1)
        self.assertEqual(isolated(E("1"), E("1")).canonize()[2], 2)
        graph = g.Graph(
            nodes=nodes("a", "b"), edges=[link("a", "b", E("p")), link("a", "b", 1)]
        )
        # Numbers sort before expressions.
        self.assertEqual(graph.canonize()[0].edge(0).data, 1)

    def test_key_callbacks_project_payloads(self):
        a = isolated({"kind": 1}, {"kind": 2})
        b = isolated({"kind": 2}, {"kind": 1})
        with self.assertRaisesRegex(TypeError, "node_key="):
            a.canonize()
        self.assertTrue(a.is_isomorphic(b, node_key=lambda d: d["kind"]))
        self.assertFalse(a.is_isomorphic(b, node_key=lambda d: id(d)))
        canonical = b.canonize(node_key=lambda d: -d["kind"])[0]
        self.assertEqual([n.data["kind"] for n in canonical.nodes()], [2, 1])

    def test_unhashable_keys_name_their_callback(self):
        edge = g.Graph(nodes=nodes("a"), edges=[leg("a", [1])])
        with self.assertRaisesRegex(TypeError, "edge_key="):
            edge.canonize()
        with self.assertRaisesRegex(TypeError, "edge_key="):
            edge.canonize_edges()
        self.assertEqual(edge.canonize(edge_key=tuple)[2], 1)
        half_edge = g.Graph(nodes=nodes("a"), edges=[leg("a", port={"p": 1})])
        with self.assertRaisesRegex(TypeError, "half_edge_key="):
            half_edge.is_isomorphic(half_edge)
        with self.assertRaisesRegex(TypeError, "node_key="):
            isolated(1).canonize(node_key=lambda d: [d])
        with self.assertRaisesRegex(TypeError, "node_key="):
            isolated((1, [2])).canonize()

    def test_nan_keys_are_rejected(self):
        for payload in (math.nan, (1, math.nan)):
            with self.assertRaisesRegex(ValueError, "NaN"):
                isolated(payload).canonize()
        with self.assertRaisesRegex(ValueError, "NaN"):
            isolated(1).is_isomorphic(isolated(1), node_key=lambda d: math.nan)

    def test_kinds_never_compare_equal(self):
        self.assertTrue(isolated(Fraction(2)).is_isomorphic(isolated(2)))
        canonical = isolated(1, Fraction(1, 2), 0.25).canonize()[0]
        self.assertEqual([n.data for n in canonical.nodes()], [0.25, Fraction(1, 2), 1])
        for other in (Decimal(1), 1 + 0j):
            self.assertFalse(isolated(other).is_isomorphic(isolated(1)))
            for payloads in ((1, other), (other, 1)):
                canonical, _, size, _ = isolated(*payloads).canonize()
                kinds = [type(n.data) for n in canonical.nodes()]
                self.assertEqual(kinds, [int, type(other)])
                self.assertEqual(size, 1)

    def test_other_keys_need_reprs_that_identify_them(self):
        a, b = isolated(Named("g"), Named("q")), isolated(Named("q"), Named("g"))
        self.assertTrue(a.is_isomorphic(b))
        with self.assertRaisesRegex(TypeError, "default repr.*node_key="):
            a.canonize()
        canonical = b.canonize(node_key=lambda d: d.name)[0]
        self.assertEqual([n.data.name for n in canonical.nodes()], ["g", "q"])
        # Equal keys must share their repr, and distinct keys must not.
        equal = isolated(Shown(0, "a"), Shown(1, "b"), Shown(0, "c"))
        self.assertEqual(equal.canonize(node_key=lambda d: d.name)[2], 2)
        with self.assertRaisesRegex(TypeError, "equal node keys a and c differ"):
            equal.canonize()
        with self.assertRaisesRegex(TypeError, r"share the type and repr Opaque\(a\)"):
            isolated(Opaque("a"), Opaque("a")).canonize()
        legs = g.Graph(nodes=nodes("a"), edges=[leg("a", Named(1)), leg("a", Named(2))])
        self.assertTrue(legs.is_isomorphic(legs))
        with self.assertRaisesRegex(TypeError, "edge_key="):
            legs.canonize_edges()

    def test_inconsistent_orders_raise_instead_of_panicking(self):
        canonical, _, size, _ = isolated(*[Cyclic(i) for i in range(200)]).canonize()
        self.assertEqual(len(canonical.nodes()), 200)
        self.assertEqual(size, 1)
        with self.assertRaisesRegex(ArithmeticError, "no order"):
            isolated(Unordered(1), Unordered(2)).canonize()

    def test_hash_errors_keep_their_cause(self):
        with self.assertRaisesRegex(TypeError, "bad hash") as raised:
            isolated(BadHash()).canonize()
        self.assertNotIn("node_key=", str(raised.exception))
        with self.assertRaisesRegex(TypeError, "unhashable node key") as raised:
            isolated([1]).canonize()
        self.assertIsInstance(raised.exception.__cause__, TypeError)

    def test_group_orders_beyond_128_bits_overflow(self):
        def legs(n):
            return g.Graph(nodes=nodes("a"), edges=[leg("a") for _ in range(n)])

        self.assertEqual(legs(34).canonize()[2], math.factorial(34))
        with self.assertRaises(OverflowError):
            legs(35).canonize()


if __name__ == "__main__":
    unittest.main()
