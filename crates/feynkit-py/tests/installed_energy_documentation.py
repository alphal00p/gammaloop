"""Run the energy API examples and check the docs used by notebook hovers.

Requires the installed community package and notebook development dependencies.
"""

import ast
import doctest
import inspect
import io
import re
from contextlib import redirect_stdout
from pathlib import Path

import jedi
from marimo._runtime.complete import _convert_docstring_to_markdown
from symbolica.community import hepkit as hep

CLASSES = {
    "CffRepresentation",
    "CffOrientation",
    "CrossFreeFamily",
    "CffReport",
    "LtdRepresentation",
    "LtdResidue",
    "LtdReport",
    "SurfacePair",
    "OnShellEnergy",
    "EnergySurface",
    "SurfaceFactor",
}
STUB = Path(hep.__file__).with_suffix(".pyi")
PARSER = doctest.DocTestParser()


def check_doc(name, doc, setup=""):
    assert "Examples\n--------" in doc, name
    assert PARSER.get_examples(doc), name
    namespace = {}
    for source in (setup, doc):
        for example in PARSER.get_examples(source):
            with redirect_stdout(io.StringIO()) as output:
                mode = "single" if example.want else "exec"
                # Execute trusted examples from the installed, generated stub.
                exec(compile(example.source, name, mode), namespace)  # noqa: S102
            if example.want:
                assert doctest.OutputChecker().check_output(
                    example.want, output.getvalue(), doctest.ELLIPSIS
                ), name
    rendered = _convert_docstring_to_markdown(doc)
    assert "<blockquote>" not in rendered, name
    assert '<div class="language-python codehilite">' in rendered, name
    if ":math:" in doc or "$$" in doc:
        assert "<marimo-tex" in rendered, name
        assert ":math:" not in rendered and "$$" not in rendered, name
    assert not re.search(r"[\x00-\x08\x0b\x0c\x0e-\x1f]", doc), name


checked = 0
for cls in ast.parse(STUB.read_text()).body:
    if not isinstance(cls, ast.ClassDef) or cls.name not in CLASSES:
        continue
    native = getattr(hep, cls.name)
    setup = ast.get_docstring(cls)
    assert inspect.cleandoc(native.__doc__) == setup, cls.name
    check_doc(cls.name, setup)
    checked += 1
    for method in cls.body:
        if not isinstance(method, ast.FunctionDef):
            continue
        if method.name.startswith("_repr_") or method.name == "_display_":
            continue
        name = f"{cls.name}.{method.name}"
        doc = ast.get_docstring(method)
        # CPython slot wrappers replace __repr__/__len__ documentation with
        # generic text; their authored stub examples are still checked below.
        if not method.name.startswith("__"):
            assert inspect.cleandoc(getattr(native, method.name).__doc__) == doc, name
        args = [*method.args.args, *method.args.kwonlyargs]
        expected = {arg.arg for arg in args} - {"self", "cls"}
        if expected:
            assert doc.index("Examples\n--------") < doc.index(
                "Parameters\n----------"
            ), name
            for parameter in expected:
                assert re.search(rf"^{parameter} :", doc, re.MULTILINE), name
        check_doc(name, doc, setup)
        checked += 1

# Notebook hover must retain both overloads; qualified @typing.overload caused
# Jedi to show only LTD, even when the cell called integrate_energy(method="cff").
source = (
    "from symbolica.community import hepkit as hep\n"
    "diagram: hep.FeynmanDiagram\n"
    'diagram.integrate_energy(method="cff")'
)
method = jedi.Script(source).infer(3, 20)[0]
signatures = [signature.to_string() for signature in method.get_signatures()]
assert len(signatures) == 2, signatures
assert any("CffRepresentation" in signature for signature in signatures)
assert any("LtdRepresentation" in signature for signature in signatures)
assert method.docstring(raw=True) == inspect.cleandoc(
    hep.FeynmanDiagram.integrate_energy.__doc__
)
check_doc(
    "FeynmanDiagram.integrate_energy",
    method.docstring(raw=True),
    hep.CffRepresentation.__doc__,
)
print(
    f"{checked} energy API docstrings: examples, native/stub parity and hover conversion passed; both overloads visible"
)
