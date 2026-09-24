"""Run with the installed extension, marimo, Typst, and a C++ compiler on PATH."""

import ast
import importlib.util
import inspect
import unittest
from pathlib import Path

from symbolica import E
from symbolica.community import spenso as sp

ROOT = Path(__file__).resolve().parents[3]
STUB = ROOT / "docs/api/python/spynso3.pyi"


class ApiSurfaceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = ROOT / "examples/notebooks/spenso_api_tour.py"
        spec = importlib.util.spec_from_file_location("spenso_api_tour", path)
        tour = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(tour)
        cls.examples, cls.values = tour.run_examples()
        cls.declarations = ast.parse(STUB.read_text()).body

    @classmethod
    def tearDownClass(cls):
        cls.values["_build"].cleanup()

    def test_every_declared_type_has_an_executable_display_example(self):
        declared = {n.name for n in self.declarations if isinstance(n, ast.ClassDef)}
        self.assertEqual({name for name, *_ in self.examples}, declared)
        for name, _, _, value in self.examples:
            with self.subTest(type=name):
                expected = type(sp.AUTO) if name == "_AutoIndex" else getattr(sp, name)
                self.assertIsInstance(value, expected)
                self.assertTrue(repr(value))
                self.assertNotRegex(repr(value), r"object at 0x[0-9a-f]+")

    def test_declared_members_and_parameter_names_match_runtime(self):
        for declaration in self.declarations:
            if isinstance(declaration, ast.FunctionDef):
                self.assertTrue(callable(getattr(sp, declaration.name)))
            if not isinstance(declaration, ast.ClassDef):
                continue
            owner = (
                type(sp.AUTO)
                if declaration.name == "_AutoIndex"
                else getattr(sp, declaration.name)
            )
            declared = set()
            for member in declaration.body:
                if isinstance(member, ast.AnnAssign):
                    self.assertTrue(hasattr(owner, member.target.id))
                    declared.add(member.target.id)
                if not isinstance(member, ast.FunctionDef):
                    continue
                declared.add(member.name)
                with self.subTest(type=declaration.name, member=member.name):
                    # Tensor implements Python's sequence iteration via __getitem__.
                    if owner is sp.Tensor and member.name == "__iter__":
                        self.assertEqual(len(list(self.values["tensor"])), 4)
                        continue
                    actual = getattr(owner, member.name)
                    if member.name.startswith("__") or not callable(actual):
                        continue
                    try:
                        signature = inspect.signature(actual)
                    except (ValueError, TypeError):
                        continue  # Built-in exception methods have no text signature.
                    arguments = member.args
                    names = {
                        a.arg
                        for a in arguments.posonlyargs
                        + arguments.args
                        + arguments.kwonlyargs
                    }
                    names.update(
                        a.arg for a in (arguments.vararg, arguments.kwarg) if a
                    )
                    self.assertEqual(
                        names - {"self", "cls"},
                        set(signature.parameters) - {"self", "cls"},
                    )
            for name, value in vars(owner).items():
                if not name.startswith("_") and (
                    callable(value) or isinstance(value, property)
                ):
                    self.assertIn(
                        name, declared, f"{declaration.name}.{name} has no stub"
                    )

    def test_settings_repr_preserves_every_property(self):
        for name, _, _, value in self.examples:
            if not name.endswith(("Settings", "Filter")):
                continue
            with self.subTest(type=name):
                restored = eval(repr(value), vars(sp))
                self.assertIs(type(restored), type(value))
                self.assertEqual(repr(restored), repr(value))
                declaration = next(
                    n
                    for n in self.declarations
                    if isinstance(n, ast.ClassDef) and n.name == name
                )
                for member in declaration.body:
                    if isinstance(member, ast.FunctionDef) and any(
                        isinstance(d, ast.Name) and d.id == "property"
                        for d in member.decorator_list
                    ):
                        self.assertEqual(
                            repr(getattr(restored, member.name)),
                            repr(getattr(value, member.name)),
                        )

    def test_mathematical_values_have_rich_displays_without_mutation(self):
        for name in (
            "rep",
            "rep_name",
            "slot",
            "structure",
            "A",
            "indexed",
            "tensor",
            "network",
        ):
            value = self.values[name]
            with self.subTest(value=name):
                before = repr(value)
                html = value._repr_html_()
                if isinstance(
                    value,
                    (
                        sp.Representation,
                        sp.RepresentationName,
                        sp.Slot,
                        sp.TensorStructure,
                    ),
                ):
                    self.assertIn("data-spenso-metadata", html)
                    self.assertEqual(repr(value), before)
                    continue
                if isinstance(value, sp.Tensor):
                    self.assertIn("data-spenso-explorer", html)
                    html = value.to_html(
                        settings=sp.DisplaySettings(tensor_view="matrix")
                    )
                if isinstance(value, sp.TensorNetwork):
                    self.assertIn("data-linnet-interactive", html)
                    self.assertIn("TensorNetwork", html)
                    html = value.expression().to_html()
                self.assertIn("<math", html)
                self.assertIn("data-spenso-math", html)
                self.assertTrue(value._repr_latex_())
                self.assertEqual(repr(value), before)
        network = self.values["network"]
        self.assertIn("digraph", network.to_dot())
        self.assertNotIn("digraph", str(network))

    def test_real_and_compiled_evaluators_agree(self):
        expected = [2.0, 1.0, 0.0, 3.0]
        self.assertEqual(list(self.values["evaluated"]), expected)
        self.assertEqual(list(self.values["compiled_result"]), expected)
        self.assertIn("real=True", repr(self.values["evaluator"]))

    def test_constant_dense_and_sparse_evaluators_need_no_parameters(self):
        structure = sp.TensorName("api_surface_tests::constant")(
            sp.Representation.euc(2)
        )
        for sparse in (False, True):
            for complex_coefficients in (False, True):
                with self.subTest(
                    sparse=sparse, complex_coefficients=complex_coefficients
                ):
                    first = E("1") + (1j if complex_coefficients else 0)
                    tensor = sp.Tensor.dense(structure, [first, E("0")])
                    if sparse:
                        tensor.to_sparse()
                    evaluator = tensor.evaluator({}, {}, [], iterations=1, n_cores=1)
                    self.assertEqual(
                        list(evaluator.evaluate_complex([[]])[0]),
                        [1 + (1j if complex_coefficients else 0), 0],
                    )
                    if complex_coefficients:
                        self.assertIn("real=False", repr(evaluator))
                        with self.assertRaisesRegex(ValueError, "complex coefficients"):
                            evaluator.evaluate([[]])
                    else:
                        self.assertEqual(list(evaluator.evaluate([[]])[0]), [1.0, 0.0])


if __name__ == "__main__":
    unittest.main()
