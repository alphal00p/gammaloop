import unittest

from performance_changes import requires_measurement


class PerformanceChangesTests(unittest.TestCase):
    def test_documentation_and_test_only_changes_skip_measurements(self):
        self.assertFalse(
            requires_measurement(
                [
                    "docs/architecture/ci.typ",
                    "CONTRIBUTING.typ",
                    "README.md",
                    "tests/tests/test_runs.rs",
                    ".github/scripts/test_performance_gate.py",
                    ".github/workflows/docs-pages.yml",
                    "nix-ci.nix",
                ]
            )
        )

    def test_production_build_and_benchmark_changes_require_measurements(self):
        for path in [
            "Cargo.lock",
            "Cargo.toml",
            "build.rs",
            "rust-toolchain.toml",
            "crates/gammalooprs/src/integrands/evaluation.rs",
            "crates/spenso/src/network/graph.rs",
            "crates/spenso/Cargo.toml",
            "tests/Cargo.toml",
            "crates/vakint/form_src/vakint.frm",
            "crates/vakint/templates/integral.py",
            "crates/gammalooprs/assets/input.json",
            "assets/models/sm-default.json",
            ".cargo/config.toml",
            "flake.lock",
            "nix/rust-workspace.nix",
            "nix/ci-workspace-graph.json",
            "benchmarks/performance/suite.toml",
            "benchmarks/performance/scalar.toml",
            "benchmarks/performance/gg-hhh-gl15.toml",
            "benchmarks/performance/scalar-3l-gl04-q1-squared.dot",
            ".github/scripts/performance_changes.py",
            ".github/workflows/nixci-readiness.yml",
        ]:
            with self.subTest(path=path):
                self.assertTrue(requires_measurement([path]))

    def test_mixed_changes_and_empty_diff(self):
        self.assertTrue(requires_measurement(["README.md", "crates/idenso/src/lib.rs"]))
        self.assertFalse(requires_measurement([]))


if __name__ == "__main__":
    unittest.main()
