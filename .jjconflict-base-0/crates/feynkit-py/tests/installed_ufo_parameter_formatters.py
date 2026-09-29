"""Run in a fresh interpreter: UFO parsing must not create unformatted symbols."""

import json
import tempfile
from pathlib import Path

from symbolica.community import feynkit as hep

root = Path(__file__).resolve().parents[3]
model = hep.UfoLoader(simplify_model=False).load(root / "assets/models/ufo/sm").model
assert model.particle("e-").electric_charge == -model.parameter("ee").symbol
for parameter in model.parameters:
    if parameter.texname:
        assert parameter.texname in parameter.symbol.to_latex(), parameter.name
        assert "mi(" in parameter.symbol.to_typst(), parameter.name

# The loader also accepts JSON. Preserve its labels across the intermediate
# Python model, which does not include texnames in its own JSON export.
definition = json.loads(hep.Model.scalar_qed().to_json())
for parameter in definition["parameters"]:
    if parameter["name"] == "mass":
        parameter["name"] = "importMass"
        parameter["texname"] = "m_{test}"
for particle in definition["particles"]:
    if particle["mass"] == "mass":
        particle["mass"] = "importMass"
for propagator in definition["propagators"]:
    propagator["denominator"] = propagator["denominator"].replace("mass", "importMass")
with tempfile.TemporaryDirectory() as directory:
    path = Path(directory) / "labelled.json"
    path.write_text(json.dumps(definition))
    imported = hep.UfoLoader(simplify_model=False).load(path).model
    mass = imported.parameter("importMass")
    assert mass.texname == "m_{test}"
    assert "m_{test}" in mass.symbol.to_latex()
    assert "m_{test}" in mass.symbol.to_typst()
    assert imported.particle("phi+").electric_charge == imported.parameter("e").symbol
    assert (
        imported.particle("phi-").electric_charge
        == -imported.particle("phi+").electric_charge
    )

print("UFO-first and JSON loader parameter formatters passed")
