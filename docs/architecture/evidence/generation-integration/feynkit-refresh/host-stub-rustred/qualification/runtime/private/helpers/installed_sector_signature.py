"""PROPOSED future installed-host signature check; not executed by this preflight."""
import ast
import hashlib
import inspect
import json
import sys
from pathlib import Path


def main():
    from symbolica.community import hepkit as hep
    identities = []
    for name in ("symbolica", "symbolica.core", "symbolica.community.hepkit"):
        module = sys.modules.get(name)
        path = getattr(module, "__file__", None)
        identities.append({"module": name, "file": path,
                           "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest()
                           if path and Path(path).is_file() else None})
    Path(sys.argv[1]).with_name("native-identities.json").write_text(
        json.dumps({"executable": sys.executable, "modules": identities,
                    "owners": {name: getattr(hep, name).__module__
                               for name in ("FeynmanDiagram", "IntegralFamily")}},
                   indent=2) + "\n")
    expected = ["regulator", "kinematics", "dimension", "powers", "numerator",
                "scalar_values", "auxiliary_momenta", "measure_multiplier",
                "max_order", "coefficient_expansion", "observer", "progress"]
    defaults = [inspect.Parameter.empty, None, None, None, None, None, None,
                None, 0, "physical", None, "auto"]
    records = []
    for name in ("FeynmanDiagram", "IntegralFamily"):
        signature = inspect.signature(getattr(hep, name).sector_decompose)
        parameters = [p for p in signature.parameters.values() if p.name != "self"]
        assert [p.name for p in parameters] == expected, signature
        assert all(p.kind == inspect.Parameter.KEYWORD_ONLY for p in parameters), signature
        assert [p.default for p in parameters] == defaults, signature
        records.append({"owner": name, "signature": str(signature)})
    paths = ["docs/api/python/feynkit-community.pyi",
             "crates/feynkit-py/python/symbolica/community/hepkit/__init__.pyi",
             "examples/notebooks/symbolica-host/python/symbolica/community/hepkit/__init__.pyi"]
    checks = []
    for path in paths:
        tree = ast.parse((Path("/build/source") / path).read_text())
        methods = {}
        for owner in tree.body:
            if isinstance(owner, ast.ClassDef) and owner.name in ("FeynmanDiagram", "IntegralFamily"):
                selected = [n for n in owner.body if isinstance(n, ast.FunctionDef)
                            and n.name == "sector_decompose"]
                assert len(selected) == 1, (path, owner.name)
                method = selected[0]
                assert [a.arg for a in method.args.kwonlyargs] == expected
                methods[owner.name] = ast.unparse(method.args.kwonlyargs[3].annotation)
        assert len(methods) == 2, (path, methods)
        assert "Dict[int, int]" in methods["FeynmanDiagram"], methods
        assert "Sequence[int]" in methods["IntegralFamily"], methods
        checks.append({"path": path, "methods": methods})
    assert len(records) == 2 and len(checks) == 3
    Path(sys.argv[1]).write_text(json.dumps({"native_methods": records,
        "stub_checks": checks, "backend_execution_claimed": False}, indent=2) + "\n")


if __name__ == "__main__":
    main()
