"""Validate generated user documentation against an installed native extension.

Run with the notebook dependencies (NumPy, Linnet, Typst) and a C++ compiler.
The examples exercise the same entry points users copy from the Python stubs.
"""

import ast
import doctest
import inspect
import io
import unittest
from pathlib import Path

from symbolica.community import tensor

STUB = Path(__file__).resolve().parents[3] / "docs/api/python/spynso3.pyi"


def declarations():
    """Yield each documented API entry once, combining identical overload docs."""
    module = ast.parse(STUB.read_text())
    yield "tensor", module
    seen = set()
    for declaration in module.body:
        if not isinstance(declaration, (ast.ClassDef, ast.FunctionDef)):
            continue
        entries = [(declaration.name, declaration)]
        if isinstance(declaration, ast.ClassDef):
            entries.extend(
                (f"{declaration.name}.{member.name}", member)
                for member in declaration.body
                if isinstance(member, ast.FunctionDef)
            )
        for name, node in entries:
            identity = name, ast.get_docstring(node)
            if identity not in seen:
                seen.add(identity)
                yield name, node


class DocumentationTests(unittest.TestCase):
    def test_public_entries_explain_results_and_have_examples(self):
        for name, node in declarations():
            member = name.rsplit(".", 1)[-1]
            # Python supplies the standard repr/comparison/slot help. Constructors,
            # indexing, arithmetic, and named tensor operations have authored help.
            if member.startswith("_") and member not in {
                "_AutoIndex",
                "__new__",
                "__call__",
                "__getitem__",
                "__setitem__",
                "__iter__",
                "__len__",
                "__add__",
                "__sub__",
                "__mul__",
                "__truediv__",
                "__pow__",
            }:
                continue
            with self.subTest(entry=name):
                doc = ast.get_docstring(node) or ""
                self.assertIn("Examples\n--------", doc)
                if isinstance(node, ast.FunctionDef):
                    self.assertIn("Returns\n-------", doc)
                    args = node.args
                    named = args.posonlyargs + args.args + args.kwonlyargs
                    if any(a.arg not in {"self", "cls"} for a in named) or args.vararg:
                        self.assertIn("Parameters\n----------", doc)

    def test_named_runtime_methods_and_stubs_share_documentation(self):
        for name, node in declarations():
            if not isinstance(node, ast.FunctionDef):
                continue
            parts = name.split(".")
            if parts[-1].startswith("_"):
                continue  # CPython replaces some special-method documentation.
            owner = tensor
            for part in parts:
                owner = getattr(owner, part)
            with self.subTest(entry=name):
                self.assertEqual(inspect.getdoc(owner), ast.get_docstring(node))

    def test_examples_run_as_written(self):
        parser = doctest.DocTestParser()
        for name, node in declarations():
            doc = ast.get_docstring(node) or ""
            if ">>>" not in doc:
                continue
            with self.subTest(entry=name):
                runner = doctest.DocTestRunner(optionflags=doctest.ELLIPSIS)
                output = io.StringIO()
                example = parser.get_doctest(
                    doc,
                    {},
                    name,
                    str(STUB),
                    node.lineno if hasattr(node, "lineno") else 0,
                )
                result = runner.run(example, out=output.write)
                self.assertEqual(result.failed, 0, output.getvalue())


if __name__ == "__main__":
    unittest.main()
