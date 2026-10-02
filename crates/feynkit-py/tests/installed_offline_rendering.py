"""Community rendering uses no Python graph or Typst packages, including in WASM."""

import builtins
import sys
import xml.etree.ElementTree as ET

from symbolica.community import hepkit as hep
from symbolica.community import tensor as spenso

original_import = builtins.__import__


def import_without_renderers(name, globals=None, locals=None, fromlist=(), level=0):
    if name.split(".")[0] in {"typst", "linnet", "gammaloop"}:
        raise AssertionError(f"rendering tried to import {name}")
    return original_import(name, globals, locals, fromlist, level)


builtins.__import__ = import_without_renderers
try:
    model = hep.Model.qcd()
    process = model.process(["g"], ["g"])
    assert "<svg" in process._repr_html_()
    amplitude = process.generate_amplitude(loops=2, progress=None)
    assert len(amplitude.diagrams) == 48
    assert "<svg" in amplitude._repr_html_()
    diagram = amplitude.diagrams[0]
    ET.fromstring(diagram.render())
    ET.fromstring(
        diagram.render(
            config={
                "layouts": {"impred_steps": 2},
                "template_options": {"show-particle": False},
            }
        )
    )
    expression = diagram.numerator_expression()
    assert isinstance(expression, spenso.TensorExpression)
    assert "<math" in expression.to_html()
    ET.fromstring(expression.to_svg())
    assert not hasattr(sys.modules["symbolica.community"], "linnet")
    assert "linnet" not in sys.modules
    assert "typst" not in sys.modules
    print("Offline QCD process, 48 two-loop diagrams, and tensor math rendered")
finally:
    builtins.__import__ = original_import
