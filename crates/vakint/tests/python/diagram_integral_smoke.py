"""Exercise native FeynKit graph/family objects in one Symbolica/Vakint host."""

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

if os.environ.get("VAKINT_TEST_HOST") == "standalone":
    import vakint_hepkit_standalone as fk
    from vakint_hepkit_standalone import (
        E,
        S,
        Vakint,
        VakintEvaluationMethod,
        integral_from_diagram,
    )
else:
    from symbolica import E, S
    from symbolica.community import hepkit as fk
    from symbolica.community.hepkit.vakint import (
        Vakint,
        VakintEvaluationMethod,
        integral_from_diagram,
    )


def _model_path():
    return Path(__file__).resolve().parents[1] / "fixtures/phi_3_4.json"


def _produce_tensor_graph(directory):
    model = fk.Model(_model_path())
    diagram = fk.FeynmanDiagram.from_dot(
        model,
        'digraph tensor { a -> b [particle="phi"]; a -> b [particle="phi"]; a -> b [particle="phi"]; }',
    )
    payload = json.loads(diagram.to_json())
    payload["numerator"] = "spenso::g(spenso::mink(4,mu),spenso::mink(4,nu))"
    payload["vertices"][0]["numerator"] = payload["numerator"]
    tensor = fk.FeynmanDiagram.from_json(model, json.dumps(payload))
    assert (
        "spenso::{spenso::representation,spenso::self_dual}::mink" in tensor.to_json()
    )
    (directory / "tensor.json").write_text(tensor.to_json())
    (directory / "tensor.dot").write_text(tensor.to_dot())


def _consume_tensor_graph(directory, encoding):
    # Load before parsing a tensor expression: representation metadata must come
    # from the native graph initializer, not an earlier expression in this process.
    model = fk.Model(_model_path())
    json_text = (directory / "tensor.json").read_text()
    dot_text = (directory / "tensor.dot").read_text()
    loader = (
        fk.FeynmanDiagram.from_json
        if encoding == "json"
        else fk.FeynmanDiagram.from_dot
    )
    encoded = json_text if encoding == "json" else dot_text
    tensor = loader(model, encoded)
    assert json.loads(tensor.to_json()) == json.loads(json_text)
    assert tensor.to_dot() == dot_text
    numerator = tensor.numerator_expression()
    assert fk.TensorReducer.feynkit(E("4")).reduce(numerator) == numerator

    canonical = "spenso::{spenso::representation,spenso::self_dual}::mink"
    assert canonical in encoded
    # Registered representations restore their full metadata from abbreviated tags.
    abbreviated = loader(
        model, encoded.replace(canonical, "spenso::{spenso::self_dual}::mink")
    )
    assert abbreviated.to_json() == tensor.to_json()
    abbreviated_numerator = abbreviated.numerator_expression()
    assert len(abbreviated_numerator.axes) == len(numerator.axes) == 2
    assert abbreviated_numerator.axes == numerator.axes
    assert abbreviated_numerator.to_expression() == numerator.to_expression()
    for malformed in [
        "spenso::{spenso::representation,spenso::dualizable}::mink",
    ]:
        assert canonical in encoded
        try:
            loader(model, encoded.replace(canonical, malformed))
        except fk.DiagramError as error:
            assert "redefined with new attributes" in str(error), error
        else:
            raise AssertionError("contradictory representation tags were accepted")


class DiagramIntegralSmoke(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = fk.Model(_model_path())

    def test_tensor_graph_roundtrip_preserves_registered_representation_metadata(self):
        with tempfile.TemporaryDirectory(
            prefix="vakint-tensor-roundtrip-"
        ) as directory:
            for arguments in [
                ["--produce-tensor-graph", directory],
                ["--consume-tensor-graph", directory, "json"],
                ["--consume-tensor-graph", directory, "dot"],
            ]:
                result = subprocess.run(
                    [sys.executable, __file__, *arguments],
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(
                    result.returncode,
                    0,
                    f"{arguments[0]} failed:\n{result.stdout}\n{result.stderr}",
                )

    def test_exact_family_algebra_preserves_declared_coordinates(self):
        k, p, s, d1, d2 = S(
            "py_family::k",
            "py_family::p",
            "py_family::s",
            "py_family::d1",
            "py_family::d2",
        )
        kin = fk.Kinematics(momenta=[k, p]).with_scalar_product(p, p, s)
        family = fk.IntegralFamily(
            [k],
            [p],
            [kin.scalar_product(k, k), kin.scalar_product(k - p, k - p)],
            kinematics=kin,
        )
        self.assertEqual(family.rank, 2)
        self.assertTrue(family.is_independent)
        self.assertTrue(family.is_complete)
        self.assertEqual(family.loop_momenta, [k])
        self.assertEqual(family.external_momenta, [p])
        expected = (d1 - d2 + s) / 2
        rewritten = family.rewrite_numerator(kin.scalar_product(k, p), [d1, d2])
        self.assertEqual((rewritten - expected).expand(), E("0"))
        partial = fk.IntegralFamily([k], [p], family.denominators[:1], kinematics=kin)
        self.assertFalse(partial.is_complete)
        completed = partial.complete(candidates=[kin.scalar_product(k, p)])
        self.assertTrue(completed.is_complete)
        self.assertEqual(completed.denominators[0], partial.denominators[0])

    def test_native_diagram_reduction_agrees_with_canonical_tadpole(self):
        topology = 'digraph tadpole { a -> a [particle="phi"]; }'
        with self.assertRaises(fk.DiagramError):
            fk.FeynmanDiagram.from_dot(self.model, topology)
        engine = Vakint(
            evaluation_order=[
                VakintEvaluationMethod.new_rustred_method(master_mode="retain")
            ],
            form_exe_path="/native-python-ingress-must-not-use-form",
            integral_normalization_factor="FMFTandMATAD",
        )
        diagram = next(
            master.graph
            for symbol in engine.master_symbols()
            if (master := engine.get_master(symbol)).loop_count == 1
        )
        self.assertEqual(
            diagram.propagator_family().kinematics.dimension, E("gammalooprs::dim")
        )
        family = diagram.propagator_family(kinematics=fk.Kinematics())
        complete = fk.IntegralFamily.from_diagram(diagram, kinematics=family.kinematics)
        self.assertEqual(complete.denominators, family.denominators)
        self.assertEqual(
            diagram.integral_family(kinematics=family.kinematics).denominators,
            family.denominators,
        )
        integral = integral_from_diagram(diagram, family, E("1")).to_expression()
        k = family.loop_momenta[0]
        mass_squared = (
            family.kinematics.scalar_product(k, k) - family.denominators[0]
        ).expand()
        canonical = E("topo(I1L(m2,1))", default_namespace="vakint").replace(
            E("m2", default_namespace="vakint"), mass_squared
        )
        actual = engine.evaluate_integral(integral)
        expected = engine.evaluate_integral(canonical)
        self.assertEqual((actual - expected).expand().together(), E("0"))
        with self.assertRaises(ValueError):
            integral_from_diagram(diagram, family, E("1"), powers=[True])

    def test_independent_two_loop_closure_preserves_compound_common_mass(self):
        model = fk.Model.phi4()
        topology = 'digraph figure_eight { a -> a [particle="phi", lmb_id=0]; a -> a [particle="phi", lmb_id=1]; }'
        engine = Vakint(
            evaluation_order=[
                VakintEvaluationMethod.new_rustred_method(master_mode="retain")
            ],
            form_exe_path="/native-python-closure-must-not-use-form",
            integral_normalization_factor="FMFTandMATAD",
        )
        for mass, factors, remaining, expanded in (
            (E("2"), 1, 1, False),
            (E("3/2"), 2, 2, True),
            (E("2*py_closure::mass"), 2, 1, False),
            (E("2*py_closure::mass"), 2, 1, True),
            (E("py_closure::mass+1"), 2, 1, False),
            (E("py_closure::mass+1"), 2, 1, True),
        ):
            with self.subTest(mass=str(mass), expanded=expanded):
                source = fk.FeynmanDiagram.from_dot(model, topology)
                target = fk.FeynmanDiagram.from_dot(model, topology)
                self.assertIsNot(source, target)
                families = []
                for diagram in (source, target):
                    family = diagram.propagator_family(kinematics=fk.Kinematics())
                    loop = family.loop_momenta[0]
                    mass_squared = (
                        family.kinematics.scalar_product(loop, loop)
                        - family.denominators[0]
                    ).expand()
                    symbols = mass_squared.get_all_symbols()
                    self.assertEqual(len(symbols), 1)
                    parameter = symbols[0]
                    self.assertEqual(mass_squared, parameter**2)
                    families.append(
                        fk.IntegralFamily(
                            family.loop_momenta,
                            [],
                            [d.replace(parameter, mass) for d in family.denominators],
                            kinematics=family.kinematics,
                        )
                    )
                numerator = families[0].denominators[0] ** factors
                if expanded:
                    numerator = numerator.expand()
                left = integral_from_diagram(
                    source,
                    families[0],
                    numerator,
                    powers=[remaining + factors, 1],
                    parameter_substitutions={parameter: mass},
                ).to_expression()
                right = integral_from_diagram(
                    target,
                    families[1],
                    E("1"),
                    powers=[remaining, 1],
                    parameter_substitutions={parameter: mass},
                ).to_expression()
                difference = (
                    (engine.evaluate_integral(left) - engine.evaluate_integral(right))
                    .expand()
                    .together()
                )
                self.assertEqual(difference, E("0"))

    def test_auxiliary_numerator_keeps_native_scalar_product_and_rejects_positive_power(
        self,
    ):
        diagram = fk.FeynmanDiagram.from_dot(
            self.model,
            'digraph sunset { a -> b [particle="phi"]; a -> b [particle="phi"]; a -> b [particle="phi"]; }',
        )
        family = diagram.propagator_family(kinematics=fk.Kinematics())
        product = family.kinematics.scalar_product(*family.loop_momenta)
        extended = fk.IntegralFamily(
            family.loop_momenta,
            [],
            family.denominators + [product],
            kinematics=family.kinematics,
        )
        direct = integral_from_diagram(diagram, family, product).to_expression()
        auxiliary = integral_from_diagram(
            diagram, extended, E("1"), powers=[1, 1, 1, -1]
        ).to_expression()
        self.assertEqual(direct, auxiliary)
        with self.assertRaises(Exception) as error:
            integral_from_diagram(diagram, extended, E("1"), powers=[1, 1, 1, 1])
        self.assertIn("auxiliary", str(error.exception))


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--produce-tensor-graph":
        _produce_tensor_graph(Path(sys.argv[2]))
    elif len(sys.argv) == 4 and sys.argv[1] == "--consume-tensor-graph":
        _consume_tensor_graph(Path(sys.argv[2]), sys.argv[3])
    else:
        unittest.main(verbosity=2)
