import copy
import tomllib
from fractions import Fraction
from pathlib import Path

import pytest
from symbolica import Expression

from hepkit_renorm import model_from_dict
from hepkit_renorm import tensors as tn

SPEC = Path(__file__).parents[1] / "hepkit_renorm" / "models" / "sm_one_family.toml"


def sm_doc() -> dict:
    with open(SPEC, "rb") as fh:
        return tomllib.load(fh)


def test_field_content_is_derived(sm):
    F = sm.fields
    gauge_sector = [f"{p}{r}" for r in "GWB" for p in ("", "gh_", "agh_", "b_")]
    matter = ["Q", "Qbar", "U", "Ubar", "D", "Dbar", "L", "Lbar", "E", "Ebar", "H", "Hdag"]
    sources = [f"src_{x}" for x in ["G", "gh_G", "W", "gh_W", "B", "gh_B"] + matter]
    assert list(F) == gauge_sector + matter + sources
    assert sm.factors["B"].abelian and not F["B"].slots[1:]  # U(1): no adjoint slot
    assert [s.space.key for s in F["Q"].slots] == ["spinor", "fund:G", "fund:W"]  # one generation: no flavour slot
    assert [s.rep for s in F["Qbar"].slots[1:]] == ["Fbar", "Fbar"] and F["Qbar"].charges == {"B": Fraction(-1, 6)}
    assert [s.rep for s in F["src_Q"].slots[1:]] == ["Fbar", "Fbar"]  # the source transforms like Qbar


@pytest.mark.parametrize("name, odd, gh, dim", [
    # report tab:quantum-numbers
    ("G", False, 0, 1), ("Q", True, 0, Fraction(3, 2)), ("H", False, 0, 1), ("gh_G", True, 1, 1),
    ("agh_G", True, -1, 1), ("b_G", False, 0, 2), ("src_G", True, -1, 2), ("src_gh_G", False, -2, 2),
    ("src_Qbar", False, -1, Fraction(3, 2)), ("src_Hdag", True, -1, 2), ("src_H", True, -1, 2),
])
def test_quantum_numbers(sm, name, odd, gh, dim):
    fi = sm.fields[name]
    assert (fi.odd, fi.ghost_number, fi.dimension) == (odd, gh, dim)


def test_conjugation_conventions(sm):
    F = sm.fields
    # report 03b sec:Qbasis: gh^dag = gh, agh^dag = -agh, rho^dag = -rho, zeta^dag = -zeta, Rbar = R^dag gamma0
    assert (F["gh_W"].conjugate, F["gh_W"].conj_sign) == ("gh_W", 1)
    assert (F["agh_W"].conjugate, F["agh_W"].conj_sign) == ("agh_W", -1)
    assert (F["src_W"].conjugate, F["src_W"].conj_sign) == ("src_W", -1)
    assert (F["src_gh_G"].conjugate, F["src_gh_G"].conj_sign) == ("src_gh_G", -1)
    assert (F["src_Q"].conjugate, F["src_Q"].conj_sign) == ("src_Qbar", 1)
    assert (F["src_H"].conjugate, F["src_Hdag"].conjugate) == ("src_Hdag", "src_H")
    assert F["src_Q"].source_of == "Q" and F["Q"].source == "src_Q" and F["agh_G"].source is None


def test_yukawa_contractions_are_derived(sm):
    by = {y.coupling: sm.yukawa_structure(y) for y in sm.yukawas}
    assert by["y_u"] == ("PR", [("G", "delta", [(0, 1), (2, 1)]), ("W", "epsf", [(0, 2), (1, 0)])])
    assert by["y_d"] == ("PR", [("G", "delta", [(0, 1), (2, 1)]), ("W", "delta", [(0, 2), (1, 0)])])
    assert by["y_e"] == ("PR", [("W", "delta", [(0, 1), (1, 0)])])


def test_invariants_table(sm, field):
    """report tab:invariants and eq:SLR."""
    inv = sm.group_invariants(field)
    gs, g, gp = (sm.param(n) for n in ("g_s", "g", "gp"))
    q = lambda a, b: Expression.num(a) / b
    expected_S = {
        "G": {"Q": gs**2, "U": gs**2 / 2, "D": gs**2 / 2, "L": 0, "E": 0, "H": 0},
        "W": {"Q": 3 * g**2 / 2, "U": 0, "D": 0, "L": g**2 / 2, "E": 0, "H": g**2 / 2},
        "B": {"Q": gp**2 / 6, "U": 4 * gp**2 / 3, "D": gp**2 / 3, "L": gp**2 / 2, "E": gp**2, "H": gp**2 / 2},
    }
    for r, row in expected_S.items():
        for f, v in row.items():
            assert (inv["S"][(f, r)] - v).expand() == 0, (f, r)
    expected_C = {
        "Q": q(4, 3) * gs**2 + q(3, 4) * g**2 + q(1, 36) * gp**2, "U": q(4, 3) * gs**2 + q(4, 9) * gp**2,
        "D": q(4, 3) * gs**2 + q(1, 9) * gp**2, "L": q(3, 4) * g**2 + q(1, 4) * gp**2, "E": gp**2,
        "H": q(3, 4) * g**2 + q(1, 4) * gp**2,
    }
    for f, v in expected_C.items():
        assert (inv["C"][f] - v).expand() == 0, f
    assert [str(inv["C_ad"][r]) for r in "GWB"] == [str(3 * gs**2), str(2 * g**2), "0"]
    total = {r: sum((inv["S"][(f, r)] for f in sm.fermions), Expression.num(0)).expand() for r in "GWB"}
    assert [(total[r] - v).expand() for r, v in zip("GWB", (2 * gs**2, 2 * g**2, q(10, 3) * gp**2))] == [0, 0, 0]


def test_invariants_match_golden_evanescent_coefficients(sm, field, golden_counterterms):
    """E_A^(r) = -1/3 S^(r) and Etilde_F^(r) = -2/3 S^(r) (final_counterterms.json: delta_S.evanescent.SM)."""
    inv = sm.group_invariants(field)
    golden = golden_counterterms["delta_S"]["evanescent"]["SM"]
    names = {"g_s": sm.param("g_s"), "g": sm.param("g"), "gp": sm.param("gp")}
    for r in "GWB":
        s_r = sum((inv["S"][(f, r)] for f in sm.fermions), Expression.num(0))
        for key, factor in ((f"E_A^({r})", Expression.num(-1) / 3), (f"Etilde_F^({r})", Expression.num(-2) / 3)):
            gold = Expression.parse(golden[key].replace("**", "^"))
            for n, sym in names.items():
                gold = gold.replace(Expression.parse(n), sym)
            assert (factor * s_r - gold).expand() == 0, key


def test_anomaly(sm, field):
    assert field.is_zero(sm.anomaly_tensor(field))  # report eq:Dabc: D^{abc} = 0 for one family
    doc = sm_doc()
    doc["fermions"] = {k: v for k, v in doc["fermions"].items() if k not in ("L", "E")}
    doc["yukawa"] = [y for y in doc["yukawa"] if y["coupling"] != "y_e"]
    quarks_only = model_from_dict(doc)
    assert not field.is_zero(quarks_only.anomaly_tensor(field))  # V7 control: no leptons, D^{BBB} != 0


def test_three_generations_use_a_flavour_index():
    doc = sm_doc()
    doc["generations"] = 3
    m3 = model_from_dict(doc)
    assert [s.space.key for s in m3.fields["Q"].slots] == ["spinor", "fund:G", "fund:W", "flavour"]
    assert len(m3.fields["H"].slots) == 1
    i, j = tn.flav(1), tn.flav(2)
    assert m3.flavour_tensor("y_u", i, j) == m3.param("y_u")(i, j)


@pytest.mark.parametrize("mutate, error", [
    (lambda d: d["fermions"]["U"].update(chirality="L"), ValueError),  # Qbar H U with equal chiralities
    (lambda d: d["fermions"]["D"].update(B="1/3"), ValueError),  # hypercharge not conserved in y_d
    (lambda d: d["fermions"]["E"].update(X="F"), ValueError),  # unknown factor
    (lambda d: d["regularisation"].update(gauge_vertex="full"), NotImplementedError),
    (lambda d: d["fermions"]["E"].update(chirality="D"), NotImplementedError),
])
def test_invalid_specifications_are_rejected(mutate, error):
    doc = copy.deepcopy(sm_doc())
    mutate(doc)
    with pytest.raises(error):
        model_from_dict(doc)


def test_identities_modulo_group_algebra(sm, field):
    G = tn.factor_symbol("G")
    a, b, c = tn.adj(G, 1), tn.adj(G, 2), tn.adj(G, 3)
    i, j, k = tn.fund(G, 1), tn.fund(G, 2), tn.fund(G, 3)
    X = sm.param("lambda")  # any non-group factor
    comm = X * (tn.T(a, i, k) * tn.T(b, k, j) - tn.T(b, i, k) * tn.T(a, k, j) - tn.I * tn.f(a, b, c) * tn.T(c, i, j))
    assert sm.vanishes_modulo_group(comm.expand(), [a, b, i, j], field)
    assert not sm.vanishes_modulo_group((X * tn.T(a, i, k) * tn.T(b, k, j)).expand(), [a, b, i, j], field)
