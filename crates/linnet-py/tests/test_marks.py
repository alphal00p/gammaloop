"""Shared mark data, strict contracts, and the actual Typst rendering bridge."""

import json
import unittest
import xml.etree.ElementTree as ET

import linnet as lp


class MarkTests(unittest.TestCase):
    def test_catalogue_and_both_fits(self):
        options = {
            "triangle": {"length": "0.21cm", "width": "0.1575cm", "inset": "0%"},
            "straight": {"length": "0.16cm", "width": "0.12cm", "rev": True},
            "stealth": {"inset": "40%", "fill": lp.Color("red")},
            "round": {"stroke": lp.Stroke(paint=lp.Color("blue"))},
            "tikz": {"width": "3pt + 450%"},
            "barb": {"arc": lp.Angle.degrees(180), "rev": True},
            "hooks": {"arc": lp.Angle.degrees(120)},
            "bar": {"align": "end"},
            "bracket": {"length": "2pt", "width": "3pt"},
            "circle": {"length": lp.Ratio.from_fraction(4), "fill": False},
            "square": {"width": lp.Length.mm(2), "align": "center"},
            "diamond": {"length": "565.69%"},
            "rays": {"n": 5, "phase": lp.Angle.degrees(30)},
        }
        for name, fields in options.items():
            for fit in ("chord", "bend"):
                with self.subTest(name=name, fit=fit):
                    data = lp.Mark(name, fit=fit, shorten="100%", **fields).to_dict()
                    self.assertEqual(data["kind"], "kurvst-mark")
                    self.assertEqual(data["shape"], name)
                    self.assertEqual(data["fit"], fit)
                    self.assertIsInstance(data["shorten"], lp.Ratio)
                    # Defaults remain engine-owned, not materialized in Python.
                    self.assertEqual(
                        lp.Mark(name).to_dict(),
                        {"kind": "kurvst-mark", "shape": name},
                    )

    def test_nested_composites_fit_once(self):
        for fit in ("chord", "bend"):
            nested = lp.Mark.combine(lp.Mark("bar"), "2pt", lp.Mark("circle"))
            mark = lp.Mark.combine(nested, "50%", lp.Mark("triangle"), fit=fit)
            self.assertEqual(mark.to_dict()["shape"], "combine")
            self.assertEqual(len(mark.to_dict()["parts"]), 3)
            self.assertIsInstance(mark.to_dict()["parts"][0], dict)
        for child in (
            lp.Mark("triangle", fit="chord"),
            lp.Mark("bar", shorten="100%"),
        ):
            with self.assertRaisesRegex(ValueError, "composite|children"):
                lp.Mark.combine(child)
        with self.assertRaises(ValueError):
            lp.Mark.combine()

    def test_legacy_and_shape_specific_fields_are_rejected(self):
        for field in ("start", "end", "symbol", "scale", "anchor", "shorten_to"):
            with self.subTest(field=field):
                with self.assertRaisesRegex(ValueError, "CeTZ.*Mark"):
                    lp.Mark("triangle", **{field: 1})
                with self.assertRaisesRegex(TypeError, "replace CeTZ"):
                    lp.Mark(**{field: 1})
        for name, fields in (
            ("tikz", {"arc": lp.Angle.degrees(180)}),
            ("straight", {"fill": lp.AUTO}),
            ("bar", {"length": "2pt"}),
            ("triangle", {"n": 3}),
        ):
            with self.subTest(name=name):
                with self.assertRaises(ValueError):
                    lp.Mark(name, **fields)
        for field, value in (
            ("thickness", lp.Length.pt(2)),
            ("cap", lp.StrokeCap.Round),
            ("join", lp.StrokeJoin.Bevel),
            ("miter_limit", 5),
        ):
            with self.subTest(field=field):
                with self.assertRaisesRegex(ValueError, "paint-only"):
                    lp.Mark("triangle", stroke=lp.Stroke(**{field: value}))

    def test_sizes_round_trip_through_native_codec(self):
        for size in (
            "3pt",
            "450%",
            "3pt + 450%",
            "0.21cm",
            "2mm",
            "1e-3in",
            "1e+1pt",
            "0.5pt + -100%",
            lp.Length.pt(3),
            lp.Ratio.from_fraction(4.5),
            lp.RelativeLength(lp.Ratio.from_fraction(4.5), lp.Length.pt(3)),
        ):
            with self.subTest(size=size):
                mark = lp.Mark("triangle", length=size)
                graph = lp.build(
                    lp.node("a"),
                    lp.node("b"),
                    lp.edge(lp.source("a"), "ab", lp.sink("b")),
                )
                edge = graph.edge("ab")
                edge.drawing.style = {"flow-arrow": mark}
                codec = lp.DotCodec.linnest()
                restored = lp.Graph.from_dot(graph.to_dot(codec), codec)
                decoded = restored.edge("ab").drawing.style["flow-arrow"]
                self.assertEqual(repr(decoded.to_dict()), repr(mark.to_dict()))
        for size in ("3", "nanpt", "1em", "3pt + ", "2pt; panic()", lp.Length.em(1)):
            with self.subTest(size=size):
                with self.assertRaises((TypeError, ValueError)):
                    lp.Mark("triangle", length=size)

    def test_public_data_keeps_typed_ratios_and_paint(self):
        mark = lp.Mark(
            "stealth",
            inset="40%",
            fill="red",
            stroke=lp.Stroke(paint=lp.Color("blue")),
        )
        data = mark.to_dict()
        self.assertIsInstance(data["inset"], lp.Ratio)
        self.assertEqual(data["inset"].percent, 40)
        self.assertIsInstance(data["fill"], lp.Color)
        self.assertIsInstance(data["stroke"], lp.Color)
        self.assertEqual(
            lp.Mark("triangle", width=lp.AUTO).to_dict(),
            {
                "kind": "kurvst-mark",
                "shape": "triangle",
            },
        )

    def test_portable_native_data_and_public_dictionary_normalization(self):
        mark = lp.Mark.combine(
            lp.Mark.combine(
                lp.Mark("bar", stroke=lp.Color("blue")),
                lp.Length.pt(2),
                lp.Mark("bar"),
            ),
            lp.Ratio.from_fraction(0.5),
            lp.Mark(
                "stealth",
                length=lp.RelativeLength(lp.Ratio.from_fraction(4.5), lp.Length.pt(3)),
                inset=lp.Ratio.from_fraction(0.4),
                fill=lp.Color.rgb(1, 2, 3, 4),
            ),
            fit="bend",
            shorten="80%",
        )
        expected = {
            "mark": {
                "shape": "combine",
                "fit": "bend",
                "shorten": 0.8,
                "parts": [
                    {
                        "shape": "combine",
                        "parts": [
                            {"shape": "bar", "stroke": True},
                            {"gap": {"points": 2.0, "ratio": 0.0}},
                            {"shape": "bar"},
                        ],
                    },
                    {"gap": {"points": 0.0, "ratio": 0.5}},
                    {
                        "shape": "stealth",
                        "length": {"points": 3.0, "ratio": 4.5},
                        "inset": 0.4,
                        "fill": True,
                    },
                ],
            },
            "paints": [{"stroke": "#0074d9"}, {}, {"fill": "#01020304"}],
        }
        self.assertEqual(json.loads(json.dumps(mark.to_native())), expected)
        data = mark.to_dict()
        before = repr(data)
        self.assertEqual(lp.Mark.from_dict(data).to_native(), expected)
        self.assertEqual(repr(data), before)
        primitive = {
            "kind": "kurvst-mark",
            "shape": "triangle",
            "length": "3pt + 450%",
            "fill": "red",
        }
        before = json.dumps(primitive)
        self.assertEqual(
            lp.Mark.from_dict(primitive).to_native(),
            {
                "mark": {
                    "shape": "triangle",
                    "length": {"points": 3.0, "ratio": 4.5},
                    "fill": True,
                },
                "paints": [{"fill": "#ff4136"}],
            },
        )
        self.assertEqual(json.dumps(primitive), before)
        self.assertEqual(
            lp.Mark("triangle").to_native(),
            {"mark": {"shape": "triangle"}, "paints": [{}]},
        )
        # Returned containers are independent, including nested public data.
        data["parts"][0]["parts"][0]["stroke"] = False
        native = mark.to_native()
        native["paints"][0]["stroke"] = "red"
        self.assertEqual(mark.to_native(), expected)

    def test_public_dictionary_codec_is_strict_and_data_only(self):
        class Pretender:
            def to_native(self):
                raise AssertionError("arbitrary serialization callbacks must not run")

            def to_dict(self):
                raise AssertionError("arbitrary serialization callbacks must not run")

        for data in (
            {"end": "straight", "scale": 0.8},
            {"kind": "kurvst-mark", "shape": "triangle", "scale": 0.8},
        ):
            before = repr(data)
            with self.assertRaisesRegex(ValueError, "CeTZ"):
                lp.Mark.from_dict(data)
            self.assertEqual(repr(data), before)
        for value in (Pretender(), lambda: lp.Mark("triangle")):
            with self.assertRaises(TypeError):
                lp.Mark.from_dict(value)
            data = {"kind": "kurvst-mark", "shape": "triangle", "fill": value}
            before = repr(data)
            with self.assertRaises(TypeError):
                lp.Mark.from_dict(data)
            self.assertEqual(repr(data), before)

    def test_catalogue_and_composites_render_through_typst(self):
        shapes = (
            "triangle",
            "straight",
            "stealth",
            "round",
            "tikz",
            "barb",
            "hooks",
            "bar",
            "bracket",
            "circle",
            "square",
            "diamond",
            "rays",
        )
        for fit in ("chord", "bend"):
            marks = [
                (name, lp.Mark(name, fit=fit, stroke=lp.Color("red")))
                for name in shapes
            ]
            marks.append(
                (
                    "combine",
                    lp.Mark.combine(
                        lp.Mark("bar", stroke=lp.Color("red")),
                        "2pt",
                        lp.Mark("triangle", stroke=lp.Color("red")),
                        fit=fit,
                    ),
                )
            )
            for name, mark in marks:
                with self.subTest(name=name, fit=fit):
                    graph = lp.build(
                        lp.node("a"),
                        lp.node("b"),
                        lp.edge(lp.source("a"), "ab", lp.sink("b")),
                        render_config=lp.RenderConfig(
                            drawing=lp.DrawOptions(
                                sink_style={
                                    "stroke": lp.Stroke(
                                        paint=lp.Color("black"),
                                        thickness=lp.Length.pt(0.7),
                                    ),
                                    "mark": mark,
                                },
                            ),
                        ),
                    )
                    drawing = ET.fromstring(graph.to_svg())
                    self.assertTrue(
                        any(
                            element.get("stroke") == "#ff4136"
                            for element in drawing.iter()
                        ),
                        "the original Typst red paint must reach the rendered head",
                    )


if __name__ == "__main__":
    unittest.main()
