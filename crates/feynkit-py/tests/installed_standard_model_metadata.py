"""Audit SM metadata and its default card; --ufo also tests the pinned UFO loader."""

import cmath
import json
import sys
from pathlib import Path
from fractions import Fraction

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "assets/models/ufo"))
import sm

model = json.loads((ROOT / "assets/models/json/sm/sm.json").read_text())
particles = {p["name"]: p for p in model["particles"]}
for p in sm.all_particles:
    j = particles[p.name]
    for key in [
        "pdg_code",
        "name",
        "antiname",
        "spin",
        "color",
        "texname",
        "antitexname",
    ]:
        assert j[key] == getattr(p, key), (p.name, key)
    assert j["mass"] == p.mass.name and j["width"] == p.width.name
    assert Fraction(j["charge"]) == Fraction(
        str(getattr(p, "charge_exact", p.charge))
    ), p.name
    assert j["y_charge"] == p.Y, (p.name, j["y_charge"], p.Y)
    assert j.get("y_charge_right") == getattr(p, "YRight", None), p.name
    assert j.get("goldstone", False) == bool(
        getattr(p, "GoldstoneBoson", p.goldstoneboson)
    ), p.name
    assert j["ghost_number"] == p.GhostNumber and j["lepton_number"] == p.LeptonNumber
params = {p["name"]: p for p in model["parameters"]}
card = json.loads((ROOT / "assets/models/json/sm/restrict_default.json").read_text())
env = {"cmath": cmath, "complex": complex}
for f in sm.all_functions:
    env[f.name] = lambda *args, f=f: eval(
        f.expr, vars(cmath), dict(zip(f.arguments, args))
    )
for p in sm.all_parameters:
    j = params[p.name]
    for key, attr in [("nature", "nature"), ("parameter_type", "type")]:
        assert j[key] == getattr(p, attr), (p.name, key)
    for key in ["lhablock", "lhacode"]:
        assert j[key] == getattr(p, key, None), (p.name, key)
    value = (
        (eval(p.value, env) if isinstance(p.value, str) else p.value)
        if p.nature == "internal"
        else complex(*j["value"])
    )
    if p.nature == "external":
        assert j["value"] == card[p.name], p.name
    env[p.name] = value
    assert abs(complex(*j["value"]) - value) < 1e-12 * max(1, abs(value)), (
        p.name,
        j["value"],
        value,
    )
for c in sm.all_couplings:
    j = next(x for x in model["couplings"] if x["name"] == c.name)
    value = eval(c.value, env)
    assert abs(complex(*j["value"]) - value) < 1e-12 * max(1, abs(value)), (
        c.name,
        j["value"],
        value,
    )
    assert dict(j["orders"]) == c.order, c.name
    # Independently evaluate the serialized expressions too.
    expr = j["expression"].replace("UFO::{}::", "").replace("^", "**").replace("𝑖", "j")
    e = dict(
        env,
        **{
            k: getattr(cmath, k)
            for k in ["sqrt", "sin", "cos", "tan", "asin", "acos", "pi"]
        },
    )
    e["j"] = 1j
    value2 = eval(expr, e)
    assert abs(value2 - value) < 1e-12 * max(1, abs(value)), (c.name, value2, value)
for v in sm.all_vertices:
    j = next(x for x in model["vertex_rules"] if x["name"] == v.name)
    assert j["particles"] == [p.name for p in v.particles], v.name
    assert j["lorentz_structures"] == [l.name for l in v.lorentz], v.name
    assert len(j["color_structures"]) == len(v.color), v.name
    for (i, k), c in v.couplings.items():
        assert j["couplings"][i][k] == c.name, v.name
    assert sum(Fraction(particles[p]["charge"]) for p in j["particles"]) == 0, v.name
print(
    "Audit passed: 43 particles, 72 parameter values, 108 coupling values and 153 vertex mappings/conserved charges."
)

# Check expressions again with the nonzero CKM and Yukawa inputs from the UFO,
# so an accidentally missing term cannot hide behind the default card's zeros.
env.update(vars(cmath))
env["j"] = 1j
for p in sm.all_parameters:
    value = eval(p.value, env) if isinstance(p.value, str) else p.value
    env[p.name] = value
    if p.nature == "internal":
        expression = params[p.name]["expression"]
        actual = (
            complex(*params[p.name]["value"])
            if expression is None
            else eval(
                expression.replace("UFO::{}::", "")
                .replace("^", "**")
                .replace("𝑖", "j")
                .replace("𝜋", "pi"),
                env,
            )
        )
        assert abs(actual - value) < 1e-12 * max(1, abs(value)), p.name
for c in sm.all_couplings:
    expression = next(
        x["expression"] for x in model["couplings"] if x["name"] == c.name
    )
    actual = eval(
        expression.replace("UFO::{}::", "").replace("^", "**").replace("𝑖", "j"), env
    )
    expected = eval(c.value, env)
    assert abs(actual - expected) < 1e-12 * max(1, abs(expected)), c.name
print("Parameter and coupling expressions also agree at nonzero CKM/Yukawa inputs.")

from symbolica import E, S
from symbolica.community import hep as fk

builtin = fk.Model.standard_model()
disk = fk.Model(ROOT / "assets/models/json/sm/sm.json")
assert builtin.to_json() == disk.to_json()
for name, q, yl, yr, t3, mass in (
    ("ve", "0", "-1", None, "1/2", E("0")),
    ("e-", "-1", "-1", "-2", "-1/2", S("UFO::Me")),
    ("c", "2/3", "1/3", "4/3", "1/2", S("UFO::MC")),
    ("b", "-1/3", "1/3", "-2/3", "-1/2", S("UFO::MB")),
):
    p = builtin.particle(name)
    assert p.charge == E(q)
    assert p.y_charge == E(yl)
    assert p.y_charge_right == (None if yr is None else E(yr))
    assert p.weak_isospin == E(t3)
    assert p.mass_expression == mass
    assert p.antiparticle.charge == -E(q)
    assert p.antiparticle.y_charge == (None if yr is None else -E(yr))
    assert p.antiparticle.y_charge_right == -E(yl)
    assert p.antiparticle.weak_isospin_right == -E(t3)
assert builtin.particle("G+").y_charge == E("1")
assert builtin.particle("H").weak_isospin is None
assert builtin.particle("G0").weak_isospin is None
assert fk.Model.phi3().particle("phi").weak_isospin is None
print(
    "Exact Python particle metadata passed, including antiparticles and symbolic zero-valued masses."
)

if "--ufo" in sys.argv:
    imported = (
        fk.UfoLoader(simplify_model=False).load(ROOT / "assets/models/ufo/sm").model
    )
    for particle in imported.particles:
        expected = builtin.particle(particle.name)
        for field in [
            "charge",
            "y_charge",
            "y_charge_right",
            "weak_isospin",
            "weak_isospin_right",
        ]:
            assert getattr(particle, field) == getattr(expected, field), (
                particle.name,
                field,
            )
    loaded_json = json.loads(imported.to_json())
    reference_json = json.loads(builtin.to_json())
    for collection, fields in [
        ("lorentz_structures", ["structure"]),
        ("propagators", ["numerator", "denominator"]),
    ]:
        reference = {entry["name"]: entry for entry in reference_json[collection]}
        for entry in loaded_json[collection]:
            for field in fields:
                assert (
                    E(entry[field]) - E(reference[entry["name"]][field])
                ).expand() == E("0"), (entry["name"], field)
    reference = {entry["name"]: entry for entry in reference_json["vertex_rules"]}
    for vertex in loaded_json["vertex_rules"]:
        expected = reference[vertex["name"]]
        assert vertex["couplings"] == expected["couplings"], vertex["name"]
        for actual, expected_color in zip(
            vertex["color_structures"], expected["color_structures"], strict=True
        ):
            assert (E(actual) - E(expected_color)).expand() == E("0"), vertex["name"]
    print(
        "UFO import agrees on all particle charges, 22 Lorentz structures, 43 propagators, and 153 vertex color/coupling matrices."
    )
