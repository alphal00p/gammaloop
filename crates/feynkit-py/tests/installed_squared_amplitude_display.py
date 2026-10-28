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

# SquaredAmplitude must retain the tensor's live notebook path, including after
# summing spins. Static HTML alone bypasses inspection and keyboard zoom.
for value in (squared, fermions):
    tensor = value.expression()
    exact = tensor.to_expression()
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
