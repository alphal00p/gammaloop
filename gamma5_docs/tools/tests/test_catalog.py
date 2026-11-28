"""Checks of the operator catalogue of the unbroken phase (``site/data/unbroken.json``, built by ``gdocs.catalog``).

(a) For each role of the order store (tree, ct_singular, ct_finite, breaking) the sum over the published terms of
    coefficient x operator equals the store functional, class of external fields by class: exact plane-wave
    fingerprints over a prime field other than the one of the build, other samples, in D = 4 and 6 (tree, finite)
    and D = 4, 6, 8 (singular, breaking). The coefficients are read from the written file (``data.values``), the
    operators are rebuilt by ``catalog.operator_specs`` and are checked to be the ones displayed.
(b) The report's numbers (final_counterterms.json): every non-zero r_alpha (39) appears with its name and value, the
    20 renormalisation constants and the 8 evanescent operators appear with the report's values, and every singular
    coefficient of a tree-level operator equals its decomposition through the constants of the run report.
(c) The displayed vertex functions are generated: each part is ``mathview.vertex`` of ``rules.vertex_rule`` of the
    term's operator times its coefficient, and the vertex functions of the terms of a class add up to the vertex
    function of the store functional on that class (fingerprints).

Needs a licensed Symbolica. Run from ``gamma5_docs/tools``::

    SYMBOLICA_LICENSE=... PYTHONPATH=../../python/hepkit_renorm python -m pytest tests/test_catalog.py

``GDOCS_SITE=<dir>`` checks another site directory (default ``gamma5_docs/site``).
"""

from __future__ import annotations

import json
import re
import os
import sys
from pathlib import Path

import numpy as np
import pytest

TOOLS = Path(__file__).resolve().parent.parent
REPO = TOOLS.parent.parent
SITE = Path(os.environ.get("GDOCS_SITE", TOOLS.parent / "site"))
sys.path.insert(0, str(TOOLS))
sys.path.insert(0, str(REPO / "python" / "hepkit_renorm"))

if not os.environ.get("SYMBOLICA_LICENSE"):
    pytest.skip("SYMBOLICA_LICENSE is not set", allow_module_level=True)

from symbolica import Expression  # noqa: E402

from gdocs import catalog  # noqa: E402
from gdocs import mathview as mv  # noqa: E402
from hepkit_renorm import LocalFunctional  # noqa: E402
from hepkit_renorm.exact import second_field  # noqa: E402
from hepkit_renorm.fingerprint import Fingerprinter  # noqa: E402
from hepkit_renorm.rules import vertex_rule  # noqa: E402

ROLES = {"tree": ("tree",), "ct_singular": ("pole", "evanescent"), "ct_finite": ("finite",),
         "breaking": ("breaking",)}
DIMENSIONS = {"tree": (4, 6), "ct_singular": (4, 6, 8), "ct_finite": (4, 6), "breaking": (4, 6, 8)}


@pytest.fixture(scope="module")
def ctx():
    return catalog.Context(SITE, log=lambda *a: None)


@pytest.fixture(scope="module")
def doc():
    return json.loads((SITE / "data" / "unbroken.json").read_text())


@pytest.fixture(scope="module")
def specs(ctx):
    return {s.id: s for s in catalog.operator_specs(ctx)}


@pytest.fixture(scope="module")
def coeffs(ctx, doc):
    return {t["id"]: {p: catalog.parse_value(ctx.model, v) for p, v in t["data"]["values"].items()}
            for t in doc["terms"]}


def _store(ctx, role):
    return {"tree": ctx.S0, "ct_singular": ctx.dS, "ct_finite": ctx.DS, "breaking": ctx.Delta1}[role]


def _assembled(ctx, doc, specs, coeffs, role):
    total = LocalFunctional(ctx.fields)
    for t in doc["terms"]:
        s = specs[t["id"]]
        for p in ROLES[role]:
            if p in coeffs[t["id"]]:
                op = s.op4 if p == "finite" else s.op
                total = LocalFunctional.integrate(total + op.scale(coeffs[t["id"]][p]))
    return total


# ---------------------------------------------------------------------------------------------- the file
def test_file_shape(doc, specs):
    assert "mock" not in doc and doc["phase"] == "unbroken"
    ids = [t["id"] for t in doc["terms"]]
    assert len(ids) == len(set(ids))
    for t in doc["terms"]:
        assert t["id"] in specs, t["id"]
        assert set(t["types"]) == set(t["coefficients"]) == set(t["data"]["values"]), t["id"]
        assert set(t["rules"][0]["parts"]) == set(t["types"]), t["id"]
        assert tuple(t["fields"]) == specs[t["id"]].cls
        assert t["evanescent"] == specs[t["id"]].evanescent


def test_displayed_operators_are_the_tested_ones(doc, specs):
    for t in doc["terms"]:
        s = specs[t["id"]]
        # an operator with a finite coefficient only is shown in its four-dimensional (barred) form, the operator
        # that the restoring coefficient multiplies (report Section 5.1)
        finite_only = set(t["types"]) == {"finite"} and s.op4 is not None and not catalog._has_spinor(s.op)
        shown = s.op4 if finite_only else s.op
        assert t["operator"]["src"] == mv.M(mv.operator(shown))["src"], t["id"]
        if finite_only:
            assert "pdbar" in t["operator"]["src"] or "bar(" in t["operator"]["src"], t["id"]


def test_finite_only_terms_are_barred(doc):
    """Physics review: the twelve terms with a finite coefficient only show the barred operator in 'operator' and
    'total', and say so in the finite note."""
    ids = [t["id"] for t in doc["terms"] if set(t["types"]) == {"finite"}]
    assert len(ids) == 12
    for t in doc["terms"]:
        if t["id"] in ids:
            assert "bar(" in t["operator"]["src"] and "bar(" in t["total"]["src"], t["id"]
            assert "barred" in t["coefficients"]["finite"]["note"], t["id"]


# ---------------------------------------------------------------------------------------------- (a)
@pytest.mark.parametrize("role", list(ROLES))
def test_terms_add_up_to_the_store(role, ctx, doc, specs, coeffs):
    target = _store(ctx, role)
    built = _assembled(ctx, doc, specs, coeffs, role)
    diff = LocalFunctional.integrate(target - built)
    classes = sorted({tuple(c) for c in target.leg_classes()} | {tuple(c) for c in built.leg_classes()})
    assert classes, role
    field = second_field()
    checked = 0
    for D in DIMENSIONS[role]:
        fp = Fingerprinter(ctx.model, field, D, seed=7)
        for c in classes:
            smp = fp.samples(list(c), 4, salt=29)
            values = fp.functional(catalog.restrict(diff, c), smp) % field.p
            assert not np.any(values), f"{role}: class {c} differs from the store in D = {D}"
            checked += 1
    print(f"{role}: {len(classes)} classes, {checked} class-dimension checks")


def test_every_store_class_has_terms(ctx, doc):
    covered = {tuple(t["fields"]) for t in doc["terms"]}
    for role in ROLES:
        for c in ctx.nonzero_classes(_store(ctx, role)):
            assert c in covered, f"{role}: class {c} has no term"


# ---------------------------------------------------------------------------------------------- (b)
def test_restoring_coefficients_are_the_reports(ctx, doc, coeffs):
    final = catalog._final_r(ctx)                       # {engine name: (typst name, value)}, 39 entries
    assert len(final) == 39
    by_name = {}
    for t in doc["terms"]:
        if "finite" in t["types"]:
            mono = t["data"]["monomial"]
            assert mono in final, t["id"]
            name, val = final[mono]
            assert t["data"]["report_name"] == name, t["id"]
            assert (coeffs[t["id"]]["finite"] - val).expand() == Expression.num(0), t["id"]
            assert t["coefficients"]["finite"].get("symbol"), t["id"]
            by_name.setdefault(name, []).append(t["id"])
    assert set(by_name) == {v[0] for v in final.values()}


def test_constants_are_the_reports(ctx, doc):
    z = ctx.final["delta_S"]["Z_factors"]["values_in_kappa_units"]
    consts = doc["constants"]
    assert len(consts) == 20
    for name, _src, key, path in catalog.CONSTANTS:
        ref = z
        for k in path:
            ref = ref[k]
        val = catalog.parse_value(ctx.model, consts[name]["value"])
        assert (val - catalog.parse_value(ctx.model, ref)).expand() == Expression.num(0), name
        run = ctx.constant(key)
        if key.startswith("d[g"):
            run = (run / ctx.model.param(key[2:-1])).expand()
        assert (val - run).expand() == Expression.num(0), name


def test_evanescent_operators_are_the_reports(ctx, doc, coeffs):
    ev = ctx.final["delta_S"]["evanescent"]["SM"]
    seen = set()
    for t in doc["terms"]:
        if "evanescent" not in t["types"]:
            continue
        name = t["data"]["report_name"]
        assert name in ev, t["id"]
        assert (coeffs[t["id"]]["evanescent"] - catalog.parse_value(ctx.model, ev[name])).expand() == \
            Expression.num(0), t["id"]
        seen.add(name)
    assert seen == set(ev) and len(seen) == 8


def test_singular_coefficients_through_the_constants(ctx, doc, specs, coeffs):
    n = 0
    for t in doc["terms"]:
        s = specs[t["id"]]
        c = coeffs[t["id"]]
        if "tree" not in c or s.family not in ("Q", "extra"):
            continue
        if not s.renormalised:
            assert "pole" not in c, t["id"]
            continue
        _, val = catalog.decomposition(ctx, s, c["tree"])
        assert (val - c.get("pole", Expression.num(0))).expand() == Expression.num(0), t["id"]
        n += 1
    assert n >= 70


# ---------------------------------------------------------------------------------------------- (c)
def _rule(ctx, s, part):
    op = s.op4 if part == "finite" else s.op
    legs = catalog.rule_legs(ctx, s)
    return vertex_rule(op, list(legs), ctx.model.simplifier(expand_projectors=False), eliminate=None)


def test_vertex_functions_are_generated(ctx, doc, specs, coeffs):
    for t in doc["terms"]:
        s = specs[t["id"]]
        legs = catalog.rule_legs(ctx, s)
        assert [lg["field"] for lg in t["rules"][0]["legs"]] == list(legs), t["id"]
        for p, M in t["rules"][0]["parts"].items():
            V = mv.vertex(_rule(ctx, s, p))
            expected = catalog._M(catalog._with_loop(mv.coefficient(coeffs[t["id"]][p]) * V, p))
            # the order of terms inside sums is Symbolica's and depends on the history of the process: compare
            # the parsed (canonical) expressions when the strings differ
            assert M["src"] == expected["src"] or mv.parse(M["src"]) == mv.parse(expected["src"]), (t["id"], p)


@pytest.mark.parametrize("role", list(ROLES))
def test_vertex_functions_add_up_to_the_store(role, ctx, doc, specs, coeffs):
    """Sum over the terms of a class of coefficient x vertex function = vertex function of the store functional."""
    by_class = {}
    for t in doc["terms"]:
        for p in ROLES[role]:
            if p in coeffs[t["id"]]:
                by_class.setdefault(tuple(t["fields"]), []).append((specs[t["id"]], p, coeffs[t["id"]][p]))
    field = second_field()
    target = _store(ctx, role)
    for D in DIMENSIONS[role][:2]:
        fp = Fingerprinter(ctx.model, field, D, seed=11)
        for c, items in by_class.items():
            legs = catalog.rule_legs(ctx, items[0][0])
            store_rule = vertex_rule(catalog.restrict(target, c), list(legs), ctx.model.simplifier(), eliminate=None)
            total = Expression.num(0)
            for s, p, k in items:
                total = total + k * _rule(ctx, s, p).value
            smp = fp.samples(list(legs), 3, salt=31)
            diff = (fp.evaluate(total.expand(), smp) - fp.evaluate(store_rule.value, smp)) % field.p
            assert not np.any(diff), f"{role}: vertex functions of {c} differ from the store in D = {D}"


# ---------------------------------------------------------------------------------------------- classes
def test_term_classes(doc):
    data = json.loads((SITE / "data" / "classes.json").read_text())
    classes = {c["id"]: c for c in data["classes"]}
    for t in doc["terms"]:
        cl = t["classes"]
        for cid in cl["pole"]:
            assert classes[cid]["ghost_number"] == 0, (t["id"], cid)
        for cid in cl["breaking"]:
            assert classes[cid]["ghost_number"] == 1, (t["id"], cid)
        if set(t["types"]) == {"tree"}:
            assert not cl["pole"] and not cl["breaking"], t["id"]
        if "pole" in t["types"] or "evanescent" in t["types"]:
            assert len(cl["pole"]) == 1, t["id"]
        if "finite" in t["types"]:
            assert cl["breaking"], t["id"]


def test_classes_json_decorated():
    data = json.loads((SITE / "data" / "classes.json").read_text())
    for c in data["classes"]:
        assert c["label"] and c["label"]["html"].startswith("<math"), c["id"]
        assert c["result"] and c["result"]["pole"]["html"].startswith("<math"), c["id"]
        if c["ghost_number"] == 1:
            assert c["result"]["breaking"]["html"].startswith("<math"), c["id"]
        assert c["summary"], c["id"]


def test_review_fixes(doc):
    """Interface and physics review: names without source spellings, Yukawa names, the classes that fix r_V,
    subscript chip symbols, equation citations by number."""
    by = {t["id"]: t for t in doc["terms"]}
    for t in doc["terms"]:
        assert t.get("name_html") and "_" not in re.sub(r"<[^>]+>", "", t["name_html"]), t["id"]
    assert by["yukawa-u"]["name"].startswith("Up-quark Yukawa term −O_u")
    for f, r in (("Q", "G"), ("Q", "W"), ("U", "G"), ("D", "G"), ("L", "W")):
        b = by[f"vertex-{f}-{r}"]["classes"]["breaking"]
        assert f"g1:gh_{r},gh_{r},{f},src_{f}" in b and f"g1:gh_{r},gh_{r},{f}bar,src_{f}bar" in b
    sym = {f["id"]: f["symbol"]["text"] for f in doc["fields"]}
    assert sym["gh_G"] == "ω_G" and sym["src_W"] == "ρ_W" and sym["src_gh_W"] == "ζ_W"
    text = json.dumps(doc, ensure_ascii=False)
    assert not re.search(r"eq\. [a-z]", text)
    assert "four-dimensional operator, " not in doc["overview"][0]["text"]


def test_ghost_number_one_poles_are_evanescent_forms():
    """Physics review: the poles of the omega Psibar Psi classes are shown as multiples of Delta_1 (report
    Section 4.2: lfe 1/2 C_ad Delta_1), the others integrated by parts with a numerical check of evanescence."""
    data = json.loads((SITE / "data" / "classes.json").read_text())
    by = {c["id"]: c for c in data["classes"]}
    assert "times the breaking" in by["g1:gh_G,Q,Qbar"]["summary"]
    assert "γ̂" in by["g1:gh_G,Q,Qbar"]["result"]["pole"]["text"]
    assert "P_L P_L" not in json.dumps(data, ensure_ascii=False)
    assert "evanescent" in by["g1:G,gh_G"]["summary"]
