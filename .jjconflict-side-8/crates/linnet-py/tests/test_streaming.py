"""Public position-stream API tests, runnable with native and Pyodide wheels."""

import math
import unittest

import linnet_py as lp

DOT = """digraph {
    incoming [style=invis]
    outgoing [style=invis]
    a [id=0 pos="1,2!"]
    b [id=1]
    incoming -> a [id=0]
    a -> b [id=1 "spring-length"="1.5"]
    b -> outgoing [id=2]
}"""


class LayoutStreamTests(unittest.TestCase):
    def test_lazy_iteration_keeps_topology_and_batches_one_run(self):
        stream = lp.LayoutStream.from_dot(
            DOT, steps=5, epochs=3, every=4, early_tolerance=0
        )
        self.assertIs(iter(stream), stream)
        self.assertEqual(stream.node_names, ["a", "b"])
        self.assertEqual(stream.endpoints, [(None, 0), (0, 1), (1, None)])
        initial = next(stream)
        self.assertEqual(initial.iteration, 0)
        self.assertFalse(initial.done)
        second = next(stream)
        self.assertEqual(second.iteration, 4)
        self.assertFalse(second.done)
        frames = [initial, second, *stream]
        self.assertEqual([frame.iteration for frame in frames], [0, 4, 8, 12, 15])
        self.assertTrue(frames[-1].done)
        self.assertEqual(list(stream), [])
        for frame in frames:
            self.assertEqual(len(frame.nodes), 2)
            self.assertEqual(len(frame.edges), 3)
            self.assertEqual(frame.nodes[0], (1.0, 2.0))
            self.assertTrue(math.isfinite(frame.max_movement))
            self.assertTrue(
                all(
                    math.isfinite(x)
                    for point in frame.nodes + frame.edges
                    for x in point
                )
            )
        whole = list(
            lp.LayoutStream.from_dot(
                DOT, steps=5, epochs=3, every=100, early_tolerance=0
            )
        )[-1]
        self.assertEqual(whole.nodes, frames[-1].nodes)
        self.assertEqual(whole.edges, frames[-1].edges)
        self.assertEqual(stream.endpoints, [(None, 0), (0, 1), (1, None)])

    def test_frames_and_topology_are_independent_copies(self):
        stream = lp.LayoutStream.from_dot(DOT, every=1, steps=6, epochs=1)
        frame = next(stream)
        original_nodes, original_edges = frame.nodes, frame.edges
        frame.nodes[0] = (999.0, 999.0)
        frame.edges.clear()
        stream.node_names.clear()
        stream.endpoints.clear()
        with self.assertRaises(AttributeError):
            frame.iteration = 99
        list(stream)
        self.assertEqual(frame.nodes, original_nodes)
        self.assertEqual(frame.edges, original_edges)
        self.assertEqual(frame.iteration, 0)
        self.assertEqual(stream.node_names, ["a", "b"])
        self.assertEqual(len(stream.endpoints), 3)

    def test_zero_budget_yields_only_initial_final_frame(self):
        for kwargs in [{"steps": 0}, {"epochs": 0}]:
            with self.subTest(**kwargs):
                frames = list(lp.LayoutStream.from_dot(DOT, **kwargs))
                self.assertEqual(len(frames), 1)
                self.assertEqual(frames[0].iteration, 0)
                self.assertTrue(frames[0].done)
        empty = list(lp.LayoutStream.from_dot("digraph {}", steps=2, epochs=1))
        self.assertTrue(empty[-1].done)
        self.assertEqual(empty[-1].nodes, [])
        self.assertEqual(empty[-1].edges, [])

    def test_invalid_scalars_are_rejected_before_iteration(self):
        for kwargs in [
            {"every": 0},
            {"steps": -1},
            {"epochs": -1},
            {"seed": -1},
            {"step": float("nan")},
            {"cool": 1.1},
            {"spring_strength": -1},
            {"repulsion": float("inf")},
            {"length_scale": 0},
            {"depth_scale": -1},
            {"flattening_end": 1.1},
            {"delta": -1},
            {"early_tolerance": -1},
        ]:
            with self.subTest(**kwargs), self.assertRaises((ValueError, OverflowError)):
                lp.LayoutStream.from_dot(DOT, **kwargs)

    def test_invalid_dot_is_rejected_before_iteration(self):
        for dot in [
            "this is not DOT",
            "digraph { a -> }",
            "digraph {} digraph {}",
            'digraph { a -> b ["spring-length"="-1"] }',
        ]:
            with self.subTest(dot=dot), self.assertRaises(ValueError):
                lp.LayoutStream.from_dot(dot)


if __name__ == "__main__":
    unittest.main()
