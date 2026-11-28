"""Checks of ``gdocs.theory`` (theory cards, group summaries, manifest texts).

Needs a licensed Symbolica (``SYMBOLICA_LICENSE`` in the environment). Run from ``gamma5_docs/tools``::

    PYTHONPATH=../../python/hepkit_renorm python -m pytest tests/test_theory.py
"""

from __future__ import annotations

import os
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

TOOLS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TOOLS))

if not os.environ.get("SYMBOLICA_LICENSE"):
    pytest.skip("SYMBOLICA_LICENSE is not set", allow_module_level=True)

from gdocs import theory  # noqa: E402

MATH = re.compile(r"<math.*?</math>", re.S)


def _valid_rich(html: str) -> None:
    for m in MATH.findall(html):
        ET.fromstring(m)
    assert "[[" not in html and "{sec:" not in html and "{card:" not in html
    assert "gdocs::" not in html


def test_sources_are_consistent():
    assert theory.check() == []


def test_every_card_key_and_reference():
    cards = theory.cards()
    assert set(theory.CARD_KEYS + theory.METHOD_KEYS) <= set(cards)
    for key, c in cards.items():
        assert c["title"] and c["report"], key
        for ref in c["report"]:
            assert ref["label"].startswith("Section ") and ref["page"] >= 1
        _valid_rich(c["body"])
        for target in re.findall(r'data-theory="([^"]+)"', c["body"]):
            assert target in cards, (key, target)


def test_groups_and_texts():
    g = theory.groups()
    assert set(theory.GROUP_KEYS) <= set(g)
    for e in g.values():
        _valid_rich(e["body"])
    m = theory.manifest_texts()
    assert m["title"] and m["subtitle"] and "<" not in m["title"]
    _valid_rich(m["status"])
    for phase in ("unbroken", "broken"):
        _valid_rich(theory.phase_intro(phase))
    assert "linear gauge" in theory.phase_intro("broken")


def test_report_index_matches_headings():
    assert theory.report_ref("sec:g0_VV") == {"label": "Section 4.1.1", "page": 43}
    assert theory.report_ref("cr-fin") == {"label": "Section 5.2.3", "page": 153}


def _source(key: str) -> str:
    return theory._raw_blocks()[("card", key)]["body"]


def test_review_statements():
    """Physics review findings on the cards: the five ghost-number-one source classes (report sec:normalisation,
    existence conditions), the vanishing pole of omega omega Rbar Psi as a consistency check, and the two kinds of
    Delta_1 graphs (report 04b, first and second kind)."""
    st = _source("st-source")
    for cls in ("ord(omega, omega, rhotilde)", "ord(omega, omega, rhotilde, A)", "ord(omega, omega, omega, zsrc)",
                "ord(omega, omega, sigmadag, H)", "ord(omega, omega, Hdag, sigma)"):
        assert cls in st
    assert "ord(omega, omega, zsrc) ]]" not in st and "Without the" not in st
    assert "consistent with" in st
    br = _source("breaking")
    assert "always sits next to" not in br and "external fermion leg" in br
    status = theory.texts()["status"]["body"]
    assert "generated with HepKIT" not in status and "cross-checked against the HepKIT" in status
    assert "hand check" in theory.texts()["intro-broken"]["body"]


def test_equation_references():
    """Equations are cited by the number printed in the PDF (report_equations.json), with their page."""
    r = theory.eq_ref("eq:Delta1")
    assert r == {"label": "Eq. (76)", "page": 26}
    assert theory.eq_anchor("bl-def") == '<a data-report-page="157" title="Report, page 157">Eq. (521)</a>'
    eqs = theory.report_equations()
    assert len(eqs) > 500 and all(1 <= e["page"] <= 186 for e in eqs.values())
