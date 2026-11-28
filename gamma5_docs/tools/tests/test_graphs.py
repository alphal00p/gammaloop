"""Fast checks of the built class and graph data (``site/data/classes.json``, ``site/data/graphs``).

Reads only the built JSON files and the run report; no Symbolica needed. Run from ``gamma5_docs/tools``:

    python -m pytest tests/test_graphs.py
"""

from __future__ import annotations

import json
import os
import re
import sys
import xml.etree.ElementTree as ET
from fractions import Fraction
from pathlib import Path

import pytest

TOOLS = Path(__file__).resolve().parent.parent
SITE = Path(os.environ.get("GDOCS_SITE", TOOLS.parent / "site"))
sys.path.insert(0, str(TOOLS))

from gdocs.graphs import FIELD_ORDER, GROUPS, class_id, graph_file  # noqa: E402

SIZE_BUDGET = 25_000_000
ID_CHARS = re.compile(r"^[A-Za-z0-9_.:,-]+$")
SVG_NS = "http://www.w3.org/2000/svg"


@pytest.fixture(scope="module")
def classes():
    return json.loads((SITE / "data" / "classes.json").read_text())


@pytest.fixture(scope="module")
def graph_files(classes):
    return {c["id"]: json.loads((SITE / c["graphs"]).read_text()) for c in classes["classes"]}


@pytest.fixture(scope="module")
def run_counts():
    doc = json.loads((SITE / "resources/run_sm_one_family/run_report.json").read_text())
    return doc["results"]["1"]["graphs_per_class"]


def test_class_ids(classes):
    ids = [c["id"] for c in classes["classes"]]
    assert len(ids) == len(set(ids)) == 253
    for c in classes["classes"]:
        assert ID_CHARS.match(c["id"])
        assert c["id"] == class_id(c["fields"], c["ghost_number"], c["insertion"])
        assert c["id"] == class_id(list(reversed(c["fields"])), c["ghost_number"])        # order-independent
        assert c["fields"] == sorted(c["fields"], key=FIELD_ORDER.index)
        assert c["graphs"] == graph_file(c["id"])
        assert c["insertion"] == ("Delta1" if c["ghost_number"] == 1 else None)
        # None after gdocs.graphs, filled by gdocs.catalog.decorate_classes
        assert {"label", "result", "summary"} <= set(c)
    assert sum(c["ghost_number"] == 0 for c in classes["classes"]) == 85
    assert sum(c["ghost_number"] == 1 for c in classes["classes"]) == 168


def test_counts_equal_run_report(classes, graph_files, run_counts):
    assert len(run_counts) == len(classes["classes"])
    for c in classes["classes"]:
        key = " ".join(c["fields"])
        assert run_counts[key] == c["n_graphs"] == len(graph_files[c["id"]]["graphs"]), key
    assert sum(c["n_graphs"] for c in classes["classes"]) == 2873


def test_groups(classes):
    keys = [g[0] for g in GROUPS]
    members = [cid for g in classes["groups"] for cid in g["classes"]]
    assert sorted(members) == sorted(c["id"] for c in classes["classes"])                  # each class in one group
    by_id = {c["id"]: c for c in classes["classes"]}
    for g in classes["groups"]:
        assert g["key"] in keys
        assert all(by_id[cid]["group"] == g["key"] and by_id[cid]["ghost_number"] == g["ghost_number"]
                   for cid in g["classes"])
        assert g["report"] and all(1 <= r["page"] <= 186 for r in g["report"])
        assert g["n_graphs"] == sum(by_id[cid]["n_graphs"] for cid in g["classes"])
    assert [g["key"] for g in classes["groups"]][:18] == keys[:18]


def test_graph_records(classes, graph_files):
    ids = set()
    for cid, doc in graph_files.items():
        assert doc["class"] == cid
        for n, g in enumerate(doc["graphs"], start=1):
            assert g["id"] == f"{cid}-{n:03d}"
            ids.add(g["id"])
            assert g["topology"] in ("tadpole", "bubble", "triangle", "box", "pentagon", "hexagon")
            assert len(g["loop"]) == {"tadpole": 1, "bubble": 2, "triangle": 3, "box": 4, "pentagon": 5,
                                      "hexagon": 6}[g["topology"]]
            assert g["sign"] in ("1", "-1") and g["wick_sign"] in ("1", "-1")
            # the report's sign column: -1 exactly for a closed fermion loop or a closed ghost loop
            loop = g["loop"]
            ghost = [x.startswith("gh_") for x in loop]
            fermion = [x in ("Q", "U", "D", "L", "E", "Qbar", "Ubar", "Dbar", "Lbar", "Ebar") for x in loop]
            assert g["sign"] == ("-1" if all(ghost) or all(fermion) else "1"), g["id"]
            assert 0 < Fraction(g["symmetry_factor"]) <= 1
            assert g["insertion"] == ("Delta1" if cid.startswith("g1:") else None)
            assert sum(v["family"] == "Delta1" for v in g["vertices"]) == (1 if cid.startswith("g1:") else 0)
            assert g["caption"].startswith("<p>") and g["caption"].endswith("</p>")
    assert len(ids) == 2873


def _check_svg(svg: str) -> ET.Element:
    root = ET.fromstring(svg)
    assert root.tag == f"{{{SVG_NS}}}svg" and root.get("viewBox")
    low = svg.lower()
    assert "<script" not in low and "javascript:" not in low
    assert not re.search(r"\son[a-z]+=", low)                                            # no event handlers
    assert "href" not in low and "url(" not in low and "@import" not in low and "<image" not in low
    assert "base64" not in low and "@font-face" not in low                              # no embedded fonts
    assert not re.search(r'\sid="', svg)                                                 # no ids
    assert "#3d2645" not in svg and "data-linnet" not in svg                            # cleaned
    return root


def test_svgs(graph_files, classes):
    for doc in graph_files.values():
        for g in doc["graphs"]:
            _check_svg(g["svg"])
            assert 'stroke="currentColor"' in g["svg"]
    for e in classes["legend"]:
        _check_svg(e["svg"])
    kinds = {e["kind"] for e in classes["legend"]}
    assert {"gauge", "gluon", "fermion", "scalar", "ghost", "source", "vertex", "delta1"} <= kinds


def test_size_budget():
    size = sum(p.stat().st_size for p in (SITE / "data" / "graphs").glob("*.json"))
    assert size <= SIZE_BUDGET, size
    names = {p.name for p in (SITE / "data" / "graphs").glob("*.json")}
    assert len(names) == 253


def test_class_id_function():
    assert class_id(["G", "G"], 0) == "g0:G,G"
    assert class_id(["Qbar", "Q", "gh_G", "G"], 1, "Delta1") == "g1:G,gh_G,Q,Qbar"
    with pytest.raises(ValueError):
        class_id(["G"], 0, "Delta1")
    with pytest.raises(ValueError):
        class_id(["X"], 0)
