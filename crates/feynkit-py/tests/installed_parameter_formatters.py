"""Parameter texnames format native expressions as well as tensor expressions."""

import json
import xml.etree.ElementTree as ET

from symbolica import E
from symbolica.community import hepkit as hep
from symbolica.community import tensor as spenso

model = hep.Model.standard_model()
mass = model.particle("e-").mass
muon_mass = model.particle("mu-").mass
assert mass == model.parameter("Me").symbol
assert "m_e" in mass.to_latex()
assert r"m_{\mu}" in muon_mass.to_latex()
assert "m_e" in mass.to_typst()
for parameter in model.parameters:
    assert parameter.texname in parameter.symbol.to_latex(), parameter.name

expression = mass**2 + muon_mass**2
plain = expression.format_plain()
assert E(plain) == expression
assert "Me" in plain and "MM" in plain
assert "m_e" in expression.to_latex() and r"m_{\mu}" in expression.to_latex()
ET.fromstring(spenso.TensorExpression(expression).to_svg())
assert "m_e" in spenso.TensorExpression(expression).to_latex()

phi4 = hep.Model.phi4()
assert "m" in phi4.parameter("mass").symbol.to_latex()
assert r"\lambda" in phi4.parameter("lam").symbol.to_latex()
ET.fromstring(spenso.TensorExpression(phi4.parameter("lam").symbol).to_svg())

# A forward reference must receive its printer before any expressions are parsed.
definition = {
    "name": "parameter_printer_test",
    "orders": [],
    "parameters": [
        {
            "name": "printForward",
            "nature": "internal",
            "parameter_type": "real",
            "expression": "printLater^2",
        },
        {
            "name": "printLater",
            "nature": "external",
            "parameter_type": "real",
            "value": [1.0, 0.0],
            "texname": r"\alpha_s",
        },
    ],
    "particles": [],
    "propagators": [],
    "lorentz_structures": [],
    "couplings": [],
    "vertex_rules": [],
}
custom = hep.Model.from_json(json.dumps(definition))
reference = custom.parameter("printLater").symbol
forward = custom.parameter("printForward").expression
assert r"\alpha_s" in forward.to_latex()
plain = reference.format_plain()

# Reloading updates the label on existing expressions; removing it restores the
# ordinary symbol spelling without altering its identity or plain serialization.
definition["parameters"][1]["texname"] = r"\beta"
reloaded = hep.Model.from_json(json.dumps(definition))
assert reloaded.parameter("printLater").symbol == reference
assert r"\beta" in reference.to_latex() and r"\beta" in forward.to_latex()
definition["parameters"][1].pop("texname")
hep.Model.from_json(json.dumps(definition))
assert "printLater" in reference.to_latex()
assert "mitex" not in reference.to_typst()
assert reference.format_plain() == plain

# Labels can also be added to parameters initially declared without a texname.
definition["parameters"][0]["texname"] = r"\gamma"
definition["parameters"][0]["typstname"] = "gamma"
definition["parameters"][1]["texname"] = r"\gamma"
definition["parameters"][1]["typstname"] = "gamma"
custom = hep.Model.from_json(json.dumps(definition))
distinct = custom.parameter("printForward").symbol + reference
assert distinct.to_latex().count(r"\gamma") == 2
assert distinct.to_typst().count("gamma") == 2
assert len(distinct.get_all_symbols()) == 2

print(
    "Native parameter formatters: TeX/Typst, forward references, reloads and exact algebra passed"
)
