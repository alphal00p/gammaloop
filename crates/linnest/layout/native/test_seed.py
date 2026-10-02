"""Native EC seeds against frozen Python reference inputs and physical routes."""

import json
import os
import subprocess
import unittest
from pathlib import Path


class NativeSeedTests(unittest.TestCase):
    binary = Path(
        os.environ.get(
            "EC_SPQR_BINARY",
            Path(__file__).resolve().parents[4] / "target/ec-planarity-native/ec-spqr",
        )
    )

    def test_accepted_gallery_routes_crossings_and_coordinates(self):
        fixtures = json.loads(
            (Path(__file__).parent / "fixtures/seeds.json").read_text()
        )
        result = subprocess.run(
            [str(self.binary)],
            input="".join(
                json.dumps({"diagram": f["diagram"]}) + "\n" for f in fixtures
            ),
            capture_output=True,
            text=True,
            check=True,
            timeout=60,
        )
        lines = result.stdout.splitlines()
        self.assertEqual(len(lines), len(fixtures))
        for fixture, line in zip(fixtures, lines, strict=True):
            with self.subTest(diagram=fixture["name"]):
                got = json.loads(line)
                self.assertNotIn("error", got)
                expected = fixture["expected"]
                for key in ("routes", "node_ids", "external_ids", "edge_endpoints"):
                    self.assertEqual(got[key], expected[key], key)
                self.assertEqual(set(got["positions"]), set(expected["positions"]))
                for key, point in expected["positions"].items():
                    for actual, reference in zip(
                        got["positions"][key], point, strict=True
                    ):
                        self.assertAlmostEqual(actual, reference, delta=1e-12)
                self.assertEqual(
                    [
                        (c["id"], c["edges"], c["route_points"])
                        for c in got["crossings"]
                    ],
                    [
                        (c["id"], c["edges"], c["route_points"])
                        for c in expected["crossings"]
                    ],
                )
                for key in ("contacts", "overlaps", "degenerate"):
                    self.assertEqual(got["report"]["geometry"][key], [])


if __name__ == "__main__":
    unittest.main()
