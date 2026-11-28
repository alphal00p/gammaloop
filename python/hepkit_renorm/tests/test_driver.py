"""The order loop and the command line (hepkit_renorm.driver, python -m hepkit_renorm) on chiral QED (fast, about
half a minute with 4 worker processes; ``HEPKIT_RENORM_PROCESSES`` to change): from the specification alone to the
renormalised model file and the run report. The one-family Standard Model end to end is ``test_driver_sm.py``
(slow). Slow here: the pole identity on every ghost-number-1 class of chiral QED, and an anomalous model."""

from __future__ import annotations

import json
import os

import pytest
from test_loopdata import QED

from hepkit_renorm import model_from_dict
from hepkit_renorm.__main__ import main
from hepkit_renorm.driver import OrderNotSupported, RunOptions, missing_for_order, renormalise
from hepkit_renorm.st_solve import AnomalyError
from hepkit_renorm.store import OrderStore

pytestmark = pytest.mark.filterwarnings("ignore:This process .* is multi-threaded:DeprecationWarning")
PROCESSES = int(os.environ.get("HEPKIT_RENORM_PROCESSES", "4"))

SPEC = """name = "chiral_qed"
[gauge]
X = { group = "U(1)", coupling = "e" }
[fermions]
P = { chirality = "L", X = "1" }
N = { chirality = "R", X = "1" }
[gauge_fixing]
kind = "linear_covariant"
xi = "xi"
"""


@pytest.fixture(scope="module")
def cli_run(tmp_path_factory):
    d = tmp_path_factory.mktemp("qed")
    spec = d / "chiral_qed.toml"
    spec.write_text(SPEC)
    out = d / "out"
    status = main([str(spec), "-o", str(out), "-j", str(PROCESSES), "--no-pole-identity", "-q"])
    return status, out


def test_cli_writes_the_outputs(cli_run):
    status, out = cli_run
    assert status == 0
    names = {p.name for p in out.iterdir()}
    assert {"chiral_qed_tree.json", "chiral_qed_renormalised.json", "chiral_qed_store.json",
            "chiral_qed_order_data.json", "run_report.json", "run_report.txt", "graph_cache"} <= names
    report = json.loads((out / "run_report.json").read_text())
    assert report["status"] == "done" and report["all_checks_passed"]
    steps = [s["name"] for s in report["steps"]]
    assert steps == ["action", "tree_model", "classes", "bases", "graphs", "poles", "pole_identity", "breaking",
                     "restoring", "store", "model"]
    res = report["results"]["1"]
    # chiral QED = QED with one Dirac fermion: known one-loop constants (units 1/(16 pi^2 eps), Feynman gauge)
    assert res["constants"] == {"Z[X]": "-4/3*e^2", "Z[gh_X]": "-4/3*e^2", "Z[P]": "-e^2", "Z[N]": "-e^2",
                                "d[e]": "2/3*e^3"}
    assert res["beta"] == {"e": "4/3*e^3"}           # beta_e = e^3 / (12 pi^2)
    assert len(res["evanescent"]) == 2
    assert res["breaking"]["nonzero"] > 0 and res["restoring"]
    model = model_from_dict(QED)
    store = OrderStore.from_json(model, (out / "chiral_qed_store.json").read_text())
    assert [e.name for e in store] == ["S0", "Delta1", "deltaS1", "DeltaS1", "bd_deltaS1", "bd_DeltaS1"]
    written = json.loads((out / "chiral_qed_renormalised.json").read_text())
    assert {v["operator"] for v in written["vertex_rules"]} == {e.name for e in store}


def test_api_rerun_is_served_by_the_cache(cli_run):
    _, out = cli_run
    res = renormalise(model_from_dict(QED), RunOptions(output=out, processes=PROCESSES, pole_identity=False))
    assert res.report.passed
    stats = res.report.results["1"]["loop_statistics"]
    assert stats["evaluations"] == 0 and stats["cached"] > 0
    assert [o.order for o in res.orders] == [1] and res.store.loop_orders() == [0, 1]


def test_higher_orders_are_refused_with_what_is_missing(tmp_path):
    assert missing_for_order(1) == []
    missing = missing_for_order(2)
    assert any("back-end" in m for m in missing) and any("S_q" in m for m in missing)
    with pytest.raises(OrderNotSupported, match="order 2 is not supported"):
        renormalise(model_from_dict(QED), RunOptions(output=tmp_path, order=2))
    assert not any(tmp_path.iterdir())          # refused before any computation
    assert main(["sm_one_family", "-o", str(tmp_path / "x"), "--order", "2", "-q"]) == 2


@pytest.mark.slow
def test_pole_identity_of_chiral_qed(tmp_path):
    res = renormalise(model_from_dict(QED), RunOptions(output=tmp_path, processes=PROCESSES))
    check = next(c for c in res.report.checks if c.name.startswith("pole identity"))
    assert check.passed and "12 of 12" in check.detail


ANOMALOUS = {"name": "anomalous_u1", "gauge": {"X": {"group": "U(1)", "coupling": "e"}},
             "fermions": {"P": {"chirality": "L", "X": "1"}, "N": {"chirality": "R", "X": "2"}},
             "gauge_fixing": {"kind": "linear_covariant", "xi": "xi"}}
"""Charges 1 (left) and 2 (right): the U(1)^3 anomaly 1 - 8 does not cancel."""


@pytest.mark.slow
def test_anomalous_model_is_reported(tmp_path):
    """An anomalous U(1): the breaking from the graphs is not b-exact, the run stops at the restoring counterterm
    with the anomaly, and the report says so."""
    with pytest.raises(AnomalyError):
        renormalise(model_from_dict(ANOMALOUS), RunOptions(output=tmp_path, processes=PROCESSES,
                                                           pole_identity=False))
    report = json.loads((tmp_path / "run_report.json").read_text())
    assert report["status"] == "failed" and "AnomalyError" in report["error"]
    checks = {c["name"]: c for c in report["checks"]}
    assert not checks["no anomaly (X in b(span Q))"]["passed"]
    assert not checks["restoring counterterm exists"]["passed"]
