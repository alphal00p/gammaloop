"""Execution metadata and immutable standalone-runtime publication regressions."""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location(
    "cells", Path(__file__).with_name("prepare-documentation-cells.py")
)
cells = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cells)


class DocumentationCellsTest(unittest.TestCase):
    def test_fenced_api_example_retains_imports_and_multiline_statements(self):
        examples = list(
            cells.python_examples("""A docstring without an Examples heading.

```python
>>> from package import Graph
>>> graph = Graph(
...     name="diagram",
... )
>>> assert graph.name == "diagram"
```
""")
        )
        self.assertEqual(
            examples,
            [
                'from package import Graph\ngraph = Graph(\n    name="diagram",\n)\nassert graph.name == "diagram"'
            ],
        )

    def test_api_method_uses_the_displayed_class_prerequisite(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stubs = root / "docs/api/python"
            stubs.mkdir(parents=True)
            for name in (
                "linnet-graph",
                "feynkit-community",
                "spynso3",
                "vakint-community",
                "gammaloop-python",
            ):
                (stubs / (name + ".pyi")).write_text("")
            (stubs / "feynkit-community.pyi").write_text('''
class Matrix:
    """Examples
    --------
    >>> value = 2
    >>> assert value == 2
    """
    def squared(self):
        """Examples
        --------
        Using the setup in the ``Matrix`` class example:

        >>> assert value ** 2 == 4
        """
''')
            entries = cells.python_reference_inventory(root)
            parent, method = entries
            self.assertEqual(method["source"], "assert value ** 2 == 4")
            self.assertEqual(method["depends"], [parent["id"]])
            self.assertEqual(method["setup"], "")
            cells.validate_contexts(entries, {entry["id"]: entry for entry in entries})

    def test_tidy_keeps_hidden_execution_and_display_only_lines_distinct(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "crates/linnest/typst/src/graph.typ"
            source.parent.mkdir(parents=True)
            source.write_text("""/// ```example
/// >>>#let expected = 4
/// #let value = 2 * 2
/// <<<#value
/// #assert.eq(value, expected)
/// ```
""")
            (entry,) = cells.typst_reference_inventory(root)
            self.assertNotIn("expected = 4", entry["source"])
            self.assertIn("#let expected = 4", entry["execution"])
            self.assertIn("#value", entry["source"])
            self.assertNotIn("\n#value\n", entry["execution"])
            self.assertIn("#assert.eq(value, expected)", entry["execution"])
            self.assertIn("example-scope.typ", entry["prefix"])

    def test_context_rejects_cycles_missing_and_excluded_dependencies(self):
        entries = [{"id": "a"}, {"id": "b"}]
        for contexts in (
            {"a": {"depends": ["b"]}, "b": {"depends": ["a"]}},
            {"a": {"depends": ["missing"]}, "b": {}},
            {"a": {"depends": ["b"]}, "b": {"reason": "native only"}},
            {"a": {}},
        ):
            with self.subTest(contexts=contexts), self.assertRaises(ValueError):
                cells.validate_contexts(entries, contexts)
        cells.validate_contexts(entries, {"a": {"depends": ["b"]}, "b": {}})

    def test_manifest_and_publication_reject_content_mutation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bundle, runtime, site = (
                root / name for name in ("bundle", "runtime", "site")
            )
            bundle.mkdir()
            runtime.mkdir()
            (bundle / "python-files.json").write_text("{}")
            (bundle / "fonts.json").write_text("[]")
            (bundle / "cache-worker.js").write_text("// cache")
            (bundle / "worker.wasm").write_bytes(b"runtime bytes")
            manifest = cells.content_manifest(bundle)
            (bundle / "manifest.json").write_text(json.dumps(manifest))
            (runtime / "runtime-manifest.json").write_text(
                json.dumps(
                    {
                        "identity": "python",
                        "shared_files": {
                            "wheels/symbolica-test.whl": "wheel-hash",
                            "python_stdlib.zip": "stdlib-hash",
                        },
                    }
                )
            )
            (runtime / "wheels").mkdir()
            (runtime / "wheels/symbolica-test.whl").write_bytes(b"wheel")
            target_runtime = site / "assets/notebook-runtime/python"
            target_runtime.mkdir(parents=True)
            (target_runtime / "native-loader.js").write_text("export {}")
            (site / "assets/site.js").write_text("// site")
            cells.publish(site, runtime, bundle)
            config = (site / "assets/documentation-cells-config.json").read_bytes()
            cells.publish(site, runtime, bundle)
            self.assertEqual(
                config, (site / "assets/documentation-cells-config.json").read_bytes()
            )
            target = site / "assets/documentation-runtime" / manifest["identity"]
            (target / "worker.wasm").write_bytes(b"collision")
            with self.assertRaisesRegex(ValueError, "collision"):
                cells.publish(site, runtime, bundle)
            self.assertEqual(
                config, (site / "assets/documentation-cells-config.json").read_bytes()
            )
            (bundle / "worker.wasm").write_bytes(b"tampered")
            with self.assertRaisesRegex(ValueError, "manifest"):
                cells.publish(site, runtime, bundle)


if __name__ == "__main__":
    unittest.main()
