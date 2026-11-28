"""Comparison of an end-to-end run of the driver (``hepkit_renorm.driver``) for the one-family Standard Model with
the report, coefficient by coefficient (gates V5 to V8 end to end). Test helper (golden data from
``tests/data/final_counterterms.json`` through ``report_golden.py``); used by ``test_driver_sm.py`` and by the
comparison script of a run directory.

The run is read from its output files (the deliverables themselves): the order store (``deltaS1``, ``DeltaS1``),
the order data (``C`` and ``X_1``), the run report (constants, beta functions, counts, checks) and the
renormalised model file.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from report_golden import REMOVED, Golden, parse, report_basis
from symbolica import Expression
from test_matching import equal_external

from hepkit_renorm import LocalFunctional
from hepkit_renorm.brst import BRST
from hepkit_renorm.class_results import physical_samples
from hepkit_renorm.classes import leg_classes
from hepkit_renorm.exact import default_field
from hepkit_renorm.fingerprint import Fingerprinter
from hepkit_renorm.matching import FunctionalData, OrderContext
from hepkit_renorm.orders import counterterm_model, loop_factor
from hepkit_renorm.st_solve import ElementCondition, Normalisation, solve_restoring
from hepkit_renorm.store import OrderStore


@dataclass
class Comparison:
    name: str
    passed: bool
    detail: str = ""


@dataclass
class Run:
    """The outputs of a run directory."""

    directory: Path
    report: dict
    store: OrderStore
    C: LocalFunctional
    X: LocalFunctional
    model_file: dict


def load_run(directory: str | Path, model, name: str = "sm_one_family") -> Run:
    d = Path(directory)
    report = json.loads((d / "run_report.json").read_text())
    store = OrderStore.from_json(model, (d / f"{name}_store.json").read_text())
    od = json.loads((d / f"{name}_order_data.json").read_text())["orders"][0]

    def fun(s):
        return LocalFunctional(model.fields, Expression.parse(s), _canonical=True)

    return Run(d, report, store, fun(od["C"]), fun(od["X"]),
               json.loads((d / f"{name}_renormalised.json").read_text()))


def _g1_classes(model):
    return [c.fields for c in leg_classes(model, 1)]


def _coordinates(model, functionals, target, classes, field, seed_salt=61):
    """Coordinates of ``target`` on ``functionals`` (D = 4, physical spinors) at the fingerprinter's random
    parameter point, with the rank of the fingerprint matrix."""
    from hepkit_renorm import linalg

    fp = Fingerprinter(model, field, 4, values={"xi": 1})
    rows_a, rows_b = [], []
    for c in classes:
        n = sum(1 for f in functionals if c in f.leg_classes()) + 4
        smp = physical_samples(fp, fp.samples(list(c), n, salt=seed_salt))
        rows_a.append(np.stack([fp.functional(f, smp) for f in functionals], axis=1))
        rows_b.append(fp.functional(target, smp))
    a, b = np.concatenate(rows_a), np.concatenate(rows_b)
    return linalg.solve(a, b[:, None], field, require_unique=True)[:, 0], linalg.rank(a, field)


def _canonical_coupling(c: dict) -> dict:
    """A coupling with its expressions in canonical form (the printed order of a product depends on the order in
    which a process registered its symbols)."""
    def canon(text):
        return Expression.parse(text).expand().to_canonical_string()

    out = dict(c)
    out["expression"] = canon(c["expression"])
    if c.get("laurent"):
        out["laurent"] = dict(c["laurent"], terms={k: canon(v) for k, v in c["laurent"]["terms"].items()})
    return out


def _same_model(a: dict, b: dict) -> bool:
    if {k: v for k, v in a.items() if k != "couplings"} != {k: v for k, v in b.items() if k != "couplings"}:
        return False
    return [_canonical_coupling(c) for c in a["couplings"]] == [_canonical_coupling(c) for c in b["couplings"]]


def compare(model, run: Run, doc: dict) -> list[Comparison]:
    """Every comparison of the module docstring (the report's numbers: ``doc`` = final_counterterms.json)."""
    out: list[Comparison] = []
    res = run.report["results"]["1"]
    gold = Golden(model)
    q = gold.q_basis()
    P = lambda s: parse(model, s)

    # ------------------------------------------------------------ renormalisation constants and beta functions
    z = doc["delta_S"]["Z_factors"]["values_in_kappa_units"]
    const = {k: P(v) for k, v in res["constants"].items()}
    beta = {k: P(v) for k, v in res["beta"].items()}
    expected = {}
    for r in "GWB":
        c = model.factors[r].coupling
        expected[f"Z[{r}]"] = P(z["dZ_A"][r])
        expected[f"Z[gh_{r}]"] = P(z["dZ_omega"][r])
        expected[f"d[{c}]"] = (P(z["dZ_g"][r]) * model.param(c)).expand()
    for f in "QUDLE":
        expected[f"Z[{f}]"] = P(z["dZ_f"][f])
    expected["Z[H]"] = P(z["dZ_H"])
    for k in "ude":
        expected[f"d[y_{k}]"] = P(z["dy"][k])
    expected["d[lambda]"] = P(z["dlambda"])
    expected["d[m2]"] = P(z["dm2"])
    bad = [k for k, v in expected.items() if const.get(k) != v]
    out.append(Comparison("renormalisation constants", not bad and set(const) == set(expected),
                          f"{len(expected) - len(bad)} of {len(expected)} equal; differing: {bad}"))
    badb = [k[2:-1] for k, v in expected.items() if k.startswith("d[") and
            beta.get(k[2:-1]) != (Expression.num(2) * v).expand()]
    out.append(Comparison("beta functions = 2 delta p", not badb,
                          f"{sum(k.startswith('d[') for k in expected) - len(badb)} parameters; differing: {badb}"))

    # ------------------------------------------------------------ the singular counterterm and its evanescent part
    Cg = gold.counterterm(doc, q)
    out.append(Comparison("deltaS^(1) (closed directions continued + evanescent) equals the report",
                          equal_external(model, run.C, Cg, (4, 6, 8)), "fingerprints in D = 4, 6, 8, every class"))
    evan = res["evanescent"]
    closed = sorted("[closed]" in v for v in evan.values())
    out.append(Comparison("evanescent operators: 8, four of them b_d-closed",
                          len(evan) == 8 and closed == [False] * 4 + [True] * 4, f"{evan}"))
    out.append(Comparison("evanescent part equals the report (E_A, Etilde_F, E_H, Etilde_H)",
                          equal_external(model, LocalFunctional.integrate(run.C.evanescent(model.simplifier())),
                                         LocalFunctional.integrate(Cg.evanescent(model.simplifier())), (6, 8)),
                          "Evan[C] in D = 6, 8"))

    # ------------------------------------------------------------ X_1
    classes1 = _g1_classes(model)
    Xg = gold.breaking(doc)
    out.append(Comparison("X_1 equals the report", equal_external(model, run.X, Xg, (4,), physical_classes=classes1),
                          "fingerprints in D = 4 on physical spinors, all ghost-number-1 classes"))
    names = list(gold.n_all())
    funcs = list(gold.n_all().values())
    field = default_field()
    xo, rank = _coordinates(model, funcs, run.X, classes1, field)
    xg, _ = _coordinates(model, funcs, Xg, classes1, field)
    nz = [names[k] for k in range(len(names)) if xo[k] % field.p]
    out.append(Comparison("X_1 on the report's N monomials, coefficient by coefficient",
                          rank == len(funcs) and np.array_equal(xo % field.p, xg % field.p) and len(nz) == 71,
                          f"{len(nz)} non-zero of {len(funcs)} (rank {rank}); all equal at a random point"))

    # ------------------------------------------------------------ the restoring counterterm
    lf = loop_factor(1)
    DS = LocalFunctional.integrate(run.store.entry("DeltaS1").coefficient(0).scale(1 / lf))
    DSg = gold.restoring(doc, q)
    out.append(Comparison("DeltaS^(1) equals the report", equal_external(model, DS, DSg, (4,)),
                          f"{len(res['restoring'])} non-zero r_alpha on the generated basis"))
    out.append(Comparison("39 non-zero r_alpha on the generated basis", len(res["restoring"]) == 39,
                          f"{len(res['restoring'])}"))
    brst = BRST(model)
    ctx = OrderContext(model, brst=brst)
    norm = Normalisation([ElementCondition(k, k, v) for k, v in REMOVED.items()])
    sol = solve_restoring(ctx, FunctionalData(brst, breaking=run.X), norm, basis=report_basis(model, q))
    golden_r = doc["r"]["all_97_engine_names"]
    badr = [n for n, v in golden_r.items() if ctx.space.expression(sol.coefficients[n]) != P(v)]
    nzr = sum(bool(sol.coefficients[n]) for n in golden_r)
    out.append(Comparison("r_alpha on the report's 97 monomials (its removal lists), by name",
                          not badr and nzr == 39 and sol.unique and (sol.allowed, sol.existence_total) == (64, 8),
                          f"{97 - len(badr)} of 97 equal, {nzr} non-zero, {sol.allowed} allowed, rank "
                          f"{sol.rank_allowed}, {sol.existence_total} existence conditions; differing: {badr}"))
    out.append(Comparison("existence conditions on the generated basis: 14", len(res["existence"]) == 14,
                          f"{len(res['existence'])}"))

    # ------------------------------------------------------------ checks of the run
    checks = {c["name"]: c for c in run.report["checks"]}
    for key in ("pole identity (ghost-number-1 poles = -V_{b_d deltaS})", "b X = 0", "no anomaly (X in b(span Q))",
                "four-dimensional limit of deltaS is b-closed"):
        c = checks.get(key)
        out.append(Comparison(f"run check: {key}", bool(c and c["passed"]), c["detail"] if c else "missing"))

    # ------------------------------------------------------------ V8: the counterterm Feynman rules (05a) and the file
    import test_sm_gates as gates

    for title, fn, args in (("deltaS^(1) rules equal 05a (cr-tab-inv, cr-tab-full, cr-ev-rules)",
                             gates.test_v8_singular_rules, (model, brst, run.store)),
                            ("DeltaS^(1) rules equal 05a (eq:sts-result, cr-fin)", gates.test_v8_finite_rules,
                             (model, brst, run.store))):
        try:
            fn(*args)
            out.append(Comparison(title, True, "class by class"))
        except AssertionError as exc:
            out.append(Comparison(title, False, str(exc)[:300]))
    written = json.loads(counterterm_model(model, run.store).json_string())
    out.append(Comparison("model file = the writer on the store", _same_model(written, run.model_file),
                          f"{len(run.model_file['vertex_rules'])} vertex rules, {len(run.model_file['couplings'])} "
                          "couplings (coupling expressions compared as expressions: the printed order of a product "
                          "depends on the order in which a process registered its symbols)"))
    try:
        gates.test_v8_written_rules_equal_entries(model, run.store, field)
        out.append(Comparison("written counterterm vertices (UFO strings) equal i V of the store entries", True,
                              "252 vertex groups"))
    except AssertionError as exc:
        out.append(Comparison("written counterterm vertices (UFO strings) equal i V of the store entries", False,
                              str(exc)[:300]))
    return out
