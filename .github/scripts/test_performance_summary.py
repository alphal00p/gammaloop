import json
import tempfile
import unittest
from pathlib import Path

from performance_summary import summary


class PerformanceSummaryTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.manifest = Path(directory.name) / "run.json"
        self.ci = {
            "candidate": "a" * 40,
            "baseline": "b" * 40,
            "event": "pull_request",
            "pr": 110,
            "run_id": 17,
            "run_attempt": 3,
        }

    def capture(self, cases):
        self.manifest.write_text(
            json.dumps({"ci": self.ci, "cases": cases, "status": 2})
        )
        return summary(self.manifest)

    def test_mixed_outcome_from_older_capture_and_exact_retry_link(self):
        case = self.manifest.parent / "gl1"
        case.mkdir()
        report = {
            "status": 1,
            "metrics": {
                "generation_seconds": {
                    "status": 1,
                    "discovery": {"median_relative_delta": 0.43},
                    "confirmation": {"median_relative_delta": 0.4075},
                }
            },
            "correctness_mismatches": 0,
        }
        (case / "gate.json").write_text(json.dumps({"report": report}))
        text = self.capture(["missing", "gl1"])
        self.assertIn("1 failed · 1 inconclusive · 0 passed", text)
        self.assertIn("+40.75% in confirmation", text)
        self.assertLess(text.index("**gl1**"), text.index("**missing**"))
        self.assertIn(f"previews/pr-110/{self.ci['candidate']}/17-3/index.html", text)

    def test_unmeasured_and_manual_runs(self):
        self.assertIn("before it recorded a capture", summary(self.manifest))
        self.ci.update(event="push", pr=None)
        text = self.capture([])
        self.assertIn("did not record any measured cases", text)
        self.assertIn(
            f"docs/results/ci-runs/gammaloop-{self.ci['candidate']}-17-3/index.html",
            text,
        )
        self.ci["event"] = "workflow_dispatch"
        self.assertNotIn("github.io", self.capture([]))
