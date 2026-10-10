"""Check the metadata-only boundary of JUnit timing exports."""

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class ExportTestTimingsTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.report = self.directory / "junit.xml"
        self.output = self.directory / "metrics.csv"

    def export(self, xml):
        self.report.write_text(xml)
        return subprocess.run(
            [
                sys.executable,
                str(Path(__file__).with_name("export-test-timings.py")),
                "--check",
                "public-check",
                "--output",
                str(self.output),
                str(self.report),
            ],
            capture_output=True,
            text=True,
            check=False,
        )

    def test_exports_metadata_and_statuses_without_private_output(self):
        private = "PRIVATE_CAPTURED_PAYLOAD"
        self.output.write_text("stale CSV")
        result = self.export(f"""<testsuite>
          <testcase classname="public-binary" name="pass" time="1.230">
            <system-out>{private}</system-out></testcase>
          <testcase classname="public-binary" name="fail" time="2.5">
            <failure message="{private}">{private}</failure></testcase>
          <testcase classname="public-binary" name="error" time="3">
            <error>{private}</error><system-err>{private}</system-err></testcase>
          <testcase classname="public-binary" name="skip">
            <skipped message="{private}"/></testcase>
        </testsuite>""")
        self.assertEqual(result.returncode, 0, result.stderr)
        with self.output.open() as output:
            self.assertEqual(
                list(csv.reader(output)),
                [
                    ["check", "test_binary", "test_name", "status", "seconds"],
                    ["public-check", "public-binary", "pass", "passed", "1.230"],
                    ["public-check", "public-binary", "fail", "failed", "2.5"],
                    ["public-check", "public-binary", "error", "error", "3"],
                    ["public-check", "public-binary", "skip", "skipped", ""],
                ],
            )
        self.assertNotIn(private, self.output.read_text())
        self.assertNotIn(private, result.stdout + result.stderr)

    def test_invalid_reports_remove_stale_csv_without_partial_rows(self):
        reports = ["<testsuite>"] + [
            f'<testsuite><testcase name="valid" time="1"/>'
            f'<testcase name="invalid" time="{seconds}"/></testsuite>'
            for seconds in ["not-a-number", "NaN", "Infinity", "-1"]
        ]
        for report in reports:
            with self.subTest(report=report):
                self.output.write_text("stale CSV")
                result = self.export(report)
                self.assertEqual(result.returncode, 2)
                self.assertFalse(self.output.exists())


if __name__ == "__main__":
    unittest.main()
