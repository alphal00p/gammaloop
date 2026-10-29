"""Squared amplitudes retain their tensor's live and static notebook displays."""

import sys
from types import SimpleNamespace
from unittest.mock import patch

from symbolica.community import hepkit as hep
from symbolica.community.tensor import DisplaySettings

model = hep.Model.standard_model()
generated = model.process(
    ["e-", "e+"], ["a", "a"], vertex_allow=["V_98"]
).generate_diagrams(
    max_vertices=2, maximum_bridges=None, numerator_grouping=None, progress=None
)
squared = hep.Amplitude(generated.diagrams).squared()
fermions = squared.sum_spins([0, 1], average_initial=True)
assert squared.expression().rank == 8
assert fermions.expression().rank == 4
assert fermions.spin_summed == [0, 1]


def assert_tensor_notation(page):
    for layout in [page, *page.get("layouts", [])]:
        for internal in ("index_scope", ">bis<", ">hedge<", ">gamma<"):
            assert internal not in layout["html"], (internal, layout["html"])


# Electron–positron annihilation: the two spin projectors use distinct primed
# copies of their fermion indices. Allocating each leaf independently used to
# discard the second primed alias and expose bis(4,index_scope(bra,...)).
for value in (squared, fermions):
    exact = value.expression().to_expression()
    viewer = value.paged()
    try:
        page = viewer.page()
        assert_tensor_notation(page)
        for layout in [page, *page["layouts"]]:
            for path, _ in layout["markers"]:
                for popup in (False, True):
                    _, preview = viewer.handle_message(
                        {
                            "action": "preview",
                            "value": path,
                            "layout": layout["layout"],
                            "popup": popup,
                        }
                    )
                    assert_tensor_notation(preview)
        compact = page["layouts"][-1]
        if compact["markers"]:
            _, expanded = viewer.handle_message(
                {
                    "action": "expand",
                    "value": compact["markers"][0][0],
                    "layout": compact["layout"],
                    "size": 25,
                }
            )
            assert_tensor_notation(expanded)
        assert value.expression().to_expression() == exact
    finally:
        viewer.close()

# SquaredAmplitude must retain the tensor's live notebook path, including after
# summing spins. Static HTML alone bypasses inspection and keyboard zoom.
for value in (squared, fermions):
    tensor = value.expression()
    exact = tensor.to_expression()
    for options in ({}, {"page_size": None}, {"page_size": 25}):
        viewer = value.paged(**options)
        expected_viewer = tensor.paged(**options)
        try:
            page = viewer.page()
            expected = expected_viewer.page()
            for key in (
                "html",
                "start",
                "end",
                "total",
                "next",
                "page_size",
                "marker_ids",
            ):
                assert page[key] == expected[key], (options, key)
            if options.get("page_size") == 25:
                assert page["page_size"] == 25
        finally:
            viewer.close()
            expected_viewer.close()
    default_pages = None
    for settings in (None, DisplaySettings(index_style="raw", show_dimensions=True)):
        viewer = value.paged(page_size=1, settings=settings)
        expected_viewer = tensor.paged(page_size=1, settings=settings)
        try:
            page = viewer.page()
            assert page["page_size"] == 1
            assert page["html"] == expected_viewer.page()["html"]
            assert "<math" in page["html"]
            field, next_page = viewer.handle_message({"action": "next"})
            assert field == "page" and next_page["start"] > page["start"]
            _, expected_next = expected_viewer.handle_message({"action": "next"})
            assert next_page["html"] == expected_next["html"]
            displayed_pages = (page["html"], next_page["html"])
            if settings is None:
                default_pages = displayed_pages
            else:
                assert displayed_pages != default_pages
            field, previous_page = viewer.handle_message({"action": "previous"})
            assert field == "page" and previous_page["html"] == page["html"]
        finally:
            viewer.close()
            expected_viewer.close()
    for options in ({"page_size": 0}, {"page_size": 501}, {"notation_source": ""}):
        try:
            value.paged(**options)
        except ValueError:
            pass
        else:
            raise AssertionError("squared-amplitude pager accepted invalid options")
    with patch.dict(sys.modules, {"anywidget": None}):
        bundle = value._repr_mimebundle_()
        assert {"text/plain", "text/html"} <= bundle.keys()
        assert "application/vnd.jupyter.widget-view+json" not in bundle
        assert bundle == tensor._repr_mimebundle_()
        assert set(value._repr_mimebundle_(include=["text/plain"])) == {"text/plain"}
        assert (
            value._repr_mimebundle_(include=["text/html"], exclude=["text/html"]) == {}
        )
        # Exercise the notebook protocol without requiring an optional frontend.
        with patch.dict(
            sys.modules, {"marimo": SimpleNamespace(Html=lambda html: html)}
        ):
            assert value._display_() == tensor._display_()
    assert value.expression().to_expression() == exact

print(
    "Squared-amplitude pager, display settings, navigation and static fallback passed"
)
