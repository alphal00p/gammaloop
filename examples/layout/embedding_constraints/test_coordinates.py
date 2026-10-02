"""Geometric crossing reconstruction is independent of physical graph topology."""

import copy
import unittest

from coordinates import open_crossings


class CrossingCoordinatesTest(unittest.TestCase):
    def test_alternating_routes_get_distinct_knots_and_one_proper_crossing(self):
        positions = {
            "a": [-3, 0],
            "b": [2, 1],
            "c": [0, -4],
            "d": [-1, 3],
            "crossing": [0, 0],
        }
        routes = {"first": ["a", "crossing", "b"], "second": ["c", "crossing", "d"]}
        original = copy.deepcopy((positions, routes))
        points, paths, crossings = open_crossings(positions, routes, {"crossing"})
        self.assertEqual((positions, routes), original)
        self.assertNotIn("crossing", points)
        self.assertEqual(set(paths), set(routes))
        self.assertEqual(len(crossings), 1)
        self.assertEqual(len(set(paths["first"]) & set(paths["second"])), 0)
        for edge, route in paths.items():
            self.assertEqual((route[0], route[-1]), (routes[edge][0], routes[edge][-1]))
            self.assertEqual(len(route), 4)
        self.assertGreater(crossings[0]["clearance_radius"], 0)

    def test_two_crossing_disks_stay_separate(self):
        positions = {
            "left": [-3, 0],
            "right": [3, 0],
            "a": [-1, -2],
            "b": [-1, 2],
            "c": [1, -2],
            "d": [1, 2],
            "first": [-1, 0],
            "second": [1, 0],
        }
        routes = {
            "horizontal": ["left", "first", "second", "right"],
            "vertical-a": ["a", "first", "b"],
            "vertical-b": ["c", "second", "d"],
        }
        points, paths, crossings = open_crossings(
            positions, routes, {"first", "second"}
        )
        self.assertEqual(len(crossings), 2)
        self.assertEqual(len(paths["horizontal"]), 6)
        self.assertLess(sum(record["clearance_radius"] for record in crossings), 2)
        self.assertTrue(all(key in points for route in paths.values() for key in route))

    def test_non_alternating_routes_are_rejected(self):
        positions = {"a": [-1, 0], "b": [0, 1], "c": [1, 0], "d": [0, -1], "x": [0, 0]}
        with self.assertRaisesRegex(ValueError, "do not alternate"):
            open_crossings(
                positions, {"one": ["a", "x", "b"], "two": ["c", "x", "d"]}, {"x"}
            )


if __name__ == "__main__":
    unittest.main()
