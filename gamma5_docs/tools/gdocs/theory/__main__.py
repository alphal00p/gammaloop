"""``python -m gdocs.theory [--site DIR] [--no-classes] [--check] [--index-report REPORT_DIR]``

Writes ``data/theory.json`` and merges the group summaries into ``data/classes.json``. ``--check`` only reports
problems of the sources. ``--index-report`` regenerates ``report_sections.json`` from the report's Typst sources.
Needs ``SYMBOLICA_LICENSE`` in the environment.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from gdocs import theory


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m gdocs.theory")
    ap.add_argument("--site", type=Path, default=theory.SITE)
    ap.add_argument("--no-classes", action="store_true", help="do not touch data/classes.json")
    ap.add_argument("--check", action="store_true", help="only check the sources")
    ap.add_argument("--index-report", type=Path, metavar="REPORT_DIR",
                    help="regenerate report_sections.json from the report's Typst sources")
    a = ap.parse_args(argv)
    if a.index_report:
        doc = theory.index_report(a.index_report, a.site / "resources" / "report_headings.txt")
        print(f"report_sections.json: {len(doc['sections'])} labelled headings")
    problems = theory.check(a.site)
    for p in problems:
        print("problem:", p)
    if a.check:
        return 1 if problems else 0
    data = theory.build_theory(a.site)
    print(f"theory.json: {len(data['cards'])} cards")
    if not a.no_classes:
        doc = theory.merge_class_groups(a.site)
        n = sum(1 for g in doc["groups"] if g.get("summary"))
        print(f"classes.json: {n} of {len(doc['groups'])} groups with a summary")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
