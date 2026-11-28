"""The one-family Standard Model end to end from its specification (slow: about an hour on 128 worker processes
without a cache, minutes with a warm cache): ``renormalise("sm_one_family")`` and the comparison of its output files
with the report, coefficient by coefficient (gates V5 to V8, ``driver_golden.py``).

Run with ``pytest --runslow tests/test_driver_sm.py``; ``HEPKIT_RENORM_PROCESSES`` (default: all cores) and
``HEPKIT_RENORM_CACHE`` (a per-graph cache directory outside the repository, default: a temporary directory) apply;
``HEPKIT_RENORM_RUN`` names the output directory (default: a temporary directory).
"""

from __future__ import annotations

import os

import pytest
from driver_golden import compare, load_run

from hepkit_renorm.driver import RunOptions, renormalise

pytestmark = [pytest.mark.slow,
              pytest.mark.filterwarnings("ignore:This process .* is multi-threaded:DeprecationWarning")]


def test_standard_model_end_to_end(sm, golden_counterterms, tmp_path):
    out = os.environ.get("HEPKIT_RENORM_RUN") or tmp_path / "run"
    opts = RunOptions(output=out, cache=os.environ.get("HEPKIT_RENORM_CACHE"),
                      processes=int(os.environ.get("HEPKIT_RENORM_PROCESSES", str(os.cpu_count() or 1))))
    result = renormalise(sm, opts)
    assert result.report.passed, [c for c in result.report.checks if not c.passed]
    failed = [c for c in compare(sm, load_run(out, sm), golden_counterterms) if not c.passed]
    assert not failed, failed
