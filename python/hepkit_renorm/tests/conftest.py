"""Shared fixtures. The Symbolica licence is read from the environment only (never from a file).

The tests of this package run within Symbolica's unlicensed limits; a test that needs a licence is marked
``requires_licence`` and is skipped with a message when ``SYMBOLICA_LICENSE`` is not set.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data"


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line("markers", "requires_licence: needs SYMBOLICA_LICENSE in the environment")


def pytest_report_header(config: pytest.Config) -> str:
    return "SYMBOLICA_LICENSE: " + ("set" if os.environ.get("SYMBOLICA_LICENSE") else "not set (unlicensed mode)")


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    if os.environ.get("SYMBOLICA_LICENSE"):
        return
    skip = pytest.mark.skip(reason="SYMBOLICA_LICENSE is not set in the environment")
    for item in items:
        if "requires_licence" in item.keywords:
            item.add_marker(skip)


@pytest.fixture(scope="session")
def sm():
    """The one-family Standard Model specification shipped with the package."""
    from hepkit_renorm import builtin_model

    return builtin_model("sm_one_family")


@pytest.fixture(scope="session")
def s0_classical(sm):
    from hepkit_renorm import classical_action

    return classical_action(sm)


@pytest.fixture(scope="session")
def field():
    from hepkit_renorm.exact import default_field

    return default_field()


@pytest.fixture(scope="session")
def golden_counterterms() -> dict:
    with open(DATA / "final_counterterms.json") as fh:
        return json.load(fh)


@pytest.fixture(scope="session")
def golden_rules() -> dict:
    with open(DATA / "sm_unbroken_bmhv_rules_digest.json") as fh:
        return json.load(fh)
