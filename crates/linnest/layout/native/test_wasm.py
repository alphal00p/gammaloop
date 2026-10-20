"""Run Linnest's layout inside Typst's actual WebAssembly sandbox."""

import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


class TypstSeedTests(unittest.TestCase):
    def test_direct_layout_preserves_pins_and_supplies_physical_carriers(self):
        typst = os.environ.get("TYPST") or shutil.which("typst")
        if not typst:
            self.skipTest("Typst is required for the WASM runtime test")
        root = Path(__file__).resolve().parents[4]
        package = root / "crates/linnest/typst/src"
        source_text = f'''#import "{package}/graph.typ" as graph
#import "{package}/layout.typ": layout
#import "{package}/render/layout.typ" as renderer
#set page(width: auto, height: auto)
#context {{
  let g = graph.build({{
    graph.node(<a>, pos: graph.pos(x: graph.pin(0), y: graph.pin(0)))
    graph.node(<b>)
    graph.edge(graph.sink(<a>), id: 0)
    graph.edge(graph.source(<a>), graph.sink(<b>), id: 1)
    graph.edge(graph.source(<a>), graph.sink(<b>), id: 2)
    graph.edge(graph.source(<a>), graph.sink(<b>), id: 3)
    graph.edge(graph.source(<b>), id: 4)
  }})
  let placed = layout(g, solver: (algorithm: "impred"), impred-steps: 5,
                      impred-parallel-balance: 1)
  assert(graph.nodes(placed).at(0).pos.x == 0)
  assert(graph.nodes(placed).at(0).pos.y == 0)
  assert(graph.edges(placed).all(edge => edge.data.at("layout-carrier").len() >= 2))
  renderer.layout-graph((layouts: ((impred-steps: 5,),),), g)
}}
'''
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "direct.typ"
            destination = Path(temporary) / "direct.svg"
            source.write_text(source_text)
            result = subprocess.run(
                [typst, "compile", "--root", "/", str(source), str(destination)],
                check=False,
                text=True,
                capture_output=True,
                timeout=90,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("<svg", destination.read_text())

    def test_external_carriers_honor_explicit_node_anchors(self):
        typst = os.environ.get("TYPST") or shutil.which("typst")
        if not typst:
            self.skipTest("Typst is required for the WASM runtime test")
        root = Path(__file__).resolve().parents[4]
        source_text = r"""#import "@PACKAGE@/src/impl/draw.typ" as drawing
#import "@PACKAGE@/src/curve.typ" as curve
#let outgoing = ((0, 0), (2, 1), (4, 0))
#let anchor = (1, 0)
#let actual = drawing._dangling-path(anchor, outgoing.last(), none, (:), carrier: outgoing)
#assert(curve.segments(actual) == curve.segments(curve.split-through((anchor, (2, 1), (4, 0))).curve))
#let incoming = outgoing.rev()
#let actual = drawing._dangling-path(incoming.first(), anchor, none, (:), dangling-at-start: true, carrier: incoming)
#assert(curve.segments(actual) == curve.segments(curve.split-through(((4, 0), (2, 1), anchor)).curve))
#let ordinary = drawing._dangling-path(outgoing.first(), outgoing.last(), none, (:), carrier: outgoing)
#assert(curve.segments(ordinary) == curve.segments(curve.split-through(outgoing).curve))
Anchor checks passed.
""".replace("@PACKAGE@", str(root / "crates/linnest/typst"))
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "anchors.typ"
            source.write_text(source_text)
            result = subprocess.run(
                [
                    typst,
                    "compile",
                    "--root",
                    "/",
                    str(source),
                    str(Path(temporary) / "anchors.svg"),
                ],
                check=False,
                text=True,
                capture_output=True,
                timeout=30,
            )
            self.assertEqual(result.returncode, 0, result.stderr)

    def test_geometry_edits_invalidate_carriers_but_styling_preserves_them(self):
        typst = os.environ.get("TYPST") or shutil.which("typst")
        if not typst:
            self.skipTest("Typst is required for the WASM runtime test")
        root = Path(__file__).resolve().parents[4]
        source_text = r"""#import "@PACKAGE@/src/graph.typ" as graph
#import "@PACKAGE@/src/layout.typ": layout
#import "@PACKAGE@/src/subgraph.typ" as subgraph
#let g = graph.build({
  graph.node(<a>); graph.node(<b>)
  graph.edge(<e>, graph.source(<a>), graph.sink(<b>))
})
#let g = layout(g, layout-algo: "impred", impred-steps: 0)
#let carrier(g) = graph.edges(g).first().data.at("layout-carrier", default: none)
#let original = carrier(g)
#assert(original != none)
#let styled = graph.map(g, edge: (e: (label: [hello], label-pos: (1, 2))))
#assert(carrier(styled) == original)
#let moved = graph.map(g, node: (a: (pos: graph.pos(x: graph.pin(-10), y: graph.pin(0)))))
#assert(carrier(moved) == none)
#assert(graph.nodes(moved).first().pos.x == -10)
#let routed = graph.map(g, source: h => (route-points: ((-2, 3),)))
#assert(carrier(routed) == none)
#assert(graph.edges(routed).first().source.route-points == ((x: -2.0, y: 3.0),))
#let opened = graph.cut(g, left: subgraph.select(g, source: (<e>,)), right: subgraph.select(g, sink: (<e>,)))
#assert(graph.edges(opened).all(e => e.data.at("layout-carrier", default: none) == none))
#let left = layout(graph.build({graph.node(<l>); graph.edge(graph.sink(<l>, statement: "j"))}), layout-algo: "impred", impred-steps: 0)
#let right = layout(graph.build({graph.node(<r>); graph.edge(graph.source(<r>, statement: "j"))}), layout-algo: "impred", impred-steps: 0)
#let joined = graph.join(left, right)
#assert(graph.edges(joined).len() == 1)
#assert(carrier(joined) == none)
Cache ownership checks passed.
""".replace("@PACKAGE@", str(root / "crates/linnest/typst"))
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "carriers.typ"
            source.write_text(source_text)
            result = subprocess.run(
                [
                    typst,
                    "compile",
                    "--root",
                    "/",
                    str(source),
                    str(Path(temporary) / "carriers.svg"),
                ],
                check=False,
                text=True,
                capture_output=True,
                timeout=30,
            )
            self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
