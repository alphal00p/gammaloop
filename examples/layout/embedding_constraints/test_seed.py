"""Physical identities and metric constraints survive EC planarization."""

import json
import os
import subprocess
import sys
import unittest
from pathlib import Path

from seed import initialize


class SeedTest(unittest.TestCase):
    def check_seed(self, diagram, side_mode=True):
        seed = initialize(diagram)
        self.assertEqual(len(seed["node_ids"]), len(diagram["vertices"]))
        self.assertEqual(
            set(seed["routes"]), {str(edge["id"]) for edge in diagram["edges"]}
        )
        self.assertTrue(
            all(len(seed["routes"][edge]) == 2 for edge in seed["external_ids"])
        )
        for name in ("contacts", "overlaps", "degenerate"):
            self.assertFalse(seed["report"]["geometry"][name])
        if side_mode:
            for name, sign in (("incoming", -1), ("outgoing", 1)):
                values = {seed["positions"][key][0] for key in seed[name + "_ids"]}
                self.assertEqual(len(values), 1)
                self.assertGreater(sign * values.pop(), 0)
        for record in seed["crossings"]:
            self.assertNotIn(record["id"], seed["positions"])
            self.assertTrue(
                all(
                    key.startswith("b:")
                    for pair in record["route_points"]
                    for key in pair
                )
            )
        return seed

    def test_single_interaction_two_external_legs(self):
        self.check_seed(
            {
                "vertices": [{}],
                "edges": [
                    {"id": 0, "source": None, "target": 0},
                    {"id": 1, "source": 0, "target": None},
                ],
            }
        )

    def test_four_external_legs_share_one_owner(self):
        self.check_seed(
            {
                "vertices": [{}],
                "edges": [
                    {"id": 0, "source": None, "target": 0},
                    {"id": 1, "source": None, "target": 0},
                    {"id": 2, "source": 0, "target": None},
                    {"id": 3, "source": 0, "target": None},
                ],
            }
        )

    def test_self_loop_is_a_drawing_triangle(self):
        seed = self.check_seed(
            {
                "vertices": [{}],
                "edges": [
                    {"id": 0, "source": None, "target": 0},
                    {"id": 1, "source": 0, "target": None},
                    {"id": 2, "source": 0, "target": 0},
                ],
            }
        )
        self.assertEqual(len(seed["routes"]["2"]), 4)

    def test_single_flow_keeps_straight_cofacial_legs(self):
        seed = self.check_seed(
            {
                "vertices": [{}, {}],
                "edges": [
                    {"id": 0, "source": None, "target": 0},
                    {"id": 1, "source": None, "target": 1},
                    {"id": 2, "source": 0, "target": 1},
                ],
            },
            side_mode=False,
        )
        self.assertFalse(seed["report"]["external_side_constraints"]["realized"])

    def test_eight_fixtures_keep_physical_topology(self):
        fixtures = json.loads(Path(__file__).with_name("diagrams.json").read_text())
        for case in fixtures:
            with self.subTest(case=case["id"]):
                seed = self.check_seed(case)
                self.assertEqual(
                    len(seed["crossings"]), len(seed["report"]["geometry"]["crossings"])
                )
                if case["id"] in {"nonplanar-higgs-k33", "nonplanar-gluon-k5"}:
                    self.assertGreater(len(seed["crossings"]), 0)

    def test_coordinates_do_not_depend_on_process_hash_seed(self):
        # NetworkX's canonical-ordering ready set uses pop(). String vertex
        # hashes otherwise choose different integer-grid drawings per process.
        script = """
import json
from pathlib import Path
from seed import initialize
fixtures = json.loads(Path('diagrams.json').read_text())
print(json.dumps([initialize(case) for case in fixtures], sort_keys=True, allow_nan=False))
"""
        results = [
            subprocess.check_output(
                [sys.executable, "-c", script],
                cwd=Path(__file__).parent,
                env={**os.environ, "PYTHONHASHSEED": value},
                text=True,
                timeout=30,
            )
            for value in ("0", "1", "42")
        ]
        self.assertEqual(results[0], results[1])
        self.assertEqual(results[0], results[2])


if __name__ == "__main__":
    unittest.main()
