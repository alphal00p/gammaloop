"""Regressions for native graph inspection links in generated manual pages."""

from __future__ import annotations

import importlib.util
import json
import unittest
from html import escape
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "check_docs_html", Path(__file__).with_name("check-docs-html.py")
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class InspectionActionTests(unittest.TestCase):
    @staticmethod
    def action(kind: str, identity: int | str, fields: object) -> str:
        href = f"#linnet-{kind}-{identity}?{json.dumps(fields)}"
        return f'<a href="{escape(href, quote=True)}">inspect</a>'

    @staticmethod
    def errors(content: str) -> list[str]:
        audit = MODULE.DocumentAudit(Path("fixture.html"))
        audit.feed(
            '<!doctype html><html lang="en"><head><title>Graph</title>'
            '<meta name="viewport" content="width=device-width"></head><body>'
            '<a class="skip-link" href="#main-content">Skip</a>'
            '<nav class="portal-nav" aria-label="Primary">'
            '<a href="#main-content">One</a><a href="#main-content">Two</a>'
            '<a href="#main-content">Three</a></nav>'
            '<main id="main-content"><h1>Graph</h1>' + content + "</main></body></html>"
        )
        audit.close()
        return audit.finish(require_sidebar=False)

    def test_valid_actions_keep_boundary_and_inspection_id_conventions(self) -> None:
        edge = {"edge": 40, "source": 10, "sink": None, "orientation": "default"}
        half = {
            **edge,
            "half-edge": 81,
            "node": 10,
            "flow": "source",
            "pair": None,
        }
        actions = (
            self.action("edge", 40, edge)
            + self.action("edge", 40, {**edge, "half-edge": 81})
            + self.action("halfedge", 81, half)
            + self.action(
                "node", 10, {"edges": [40], "name": "literal &quot; & <text>"}
            )
        )
        self.assertEqual(self.errors(f"<svg\n>{actions}</svg><svg>{actions}</svg>"), [])

    def test_identity_and_json_numerical_types_match_native_contract(self) -> None:
        invalid = [
            ("node", 0, {"edges": [True]}),
            ("node", 0, {"edges": [0.0]}),
            ("node", 0, {"edges": [-1]}),
            ("node", 0, {"edges": [2**64]}),
            ("node", 0, {"edges": [], "node": True}),
            ("node", 0, {"edges": [], "node": 1}),
            ("node", 2**64, {"edges": []}),
            ("edge", 1, {"edge": 0, "source": 0, "sink": 1, "orientation": "default"}),
            (
                "edge",
                0,
                {"edge": False, "source": 0, "sink": 1, "orientation": "default"},
            ),
            (
                "edge",
                0,
                {"edge": 0, "source": 0.0, "sink": 1, "orientation": "default"},
            ),
            ("edge", 0, {"edge": 0, "source": 0, "sink": 1, "orientation": []}),
            (
                "halfedge",
                0,
                {"edge": 0, "source": 0, "sink": 1, "orientation": "default"},
            ),
            ("other", 0, {"edges": []}),
            ("node", 0, []),
        ]
        for kind, identity, fields in invalid:
            with self.subTest(kind=kind, identity=identity, fields=fields):
                errors = self.errors(
                    f"<svg>{self.action(kind, identity, fields)}</svg>"
                )
                self.assertTrue(
                    any("invalid SVG inspection action" in error for error in errors)
                )
        for payload in (
            '{"edges":[],"x":NaN}',
            '{"edges":[],"x":1e400}',
            '{"edges":[-0]}',
            '{"edges":[],"node":-0}',
            '{"edges":[0e0]}',
        ):
            with self.subTest(payload=payload):
                self.assertFalse(
                    MODULE.DocumentAudit._valid_inspection_action(
                        f"#linnet-node-0?{payload}"
                    )
                )
        self.assertTrue(
            MODULE.DocumentAudit._valid_inspection_action(
                "#linnet-node-" + "0" * 5000 + '?{"edges":[]}'
            )
        )

    def test_actions_require_actual_closed_svg_scope(self) -> None:
        action = self.action("node", 0, {"edges": []})
        for content in (
            action,
            f"<svg></svg>{action}",
            f"<svg-fake>{action}</svg-fake>",
            f"<svg:fake>{action}</svg:fake>",
            f'<div title="<svg>"></div>{action}<div title="</svg>"></div>',
            f"<svg>{action}",
        ):
            with self.subTest(content=content):
                self.assertTrue(
                    any(
                        "SVG inspection action" in error
                        for error in self.errors(content)
                    )
                )

    def test_malformed_actions_cannot_fall_back_to_dom_ids(self) -> None:
        content = (
            '<span id="linnet-node-0?invalid-json"></span>'
            '<svg><a href="#linnet-node-0?invalid-json">Bad</a></svg>'
        )
        self.assertTrue(
            any(
                "invalid SVG inspection action" in error
                for error in self.errors(content)
            )
        )

    def test_ordinary_fragments_and_aria_references_remain_checked(self) -> None:
        self.assertEqual(
            self.errors('<span id="local"></span><a href="#local">Local</a>'), []
        )
        errors = self.errors(
            '<a href="#missing">Missing</a><button aria-controls="absent">Open</button>'
        )
        self.assertTrue(
            any('href references missing id="missing"' in error for error in errors)
        )
        self.assertTrue(
            any(
                'aria-controls references missing id="absent"' in error
                for error in errors
            )
        )


if __name__ == "__main__":
    unittest.main()
