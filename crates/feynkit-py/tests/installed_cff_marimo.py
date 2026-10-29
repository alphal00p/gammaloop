"""Marimo dispatch uses the shared CFF measurement widget and native payload."""

import json
import re
from unittest.mock import patch

import marimo as mo
from symbolica.community import hepkit as hep


def external_metrics_snapshot(view):
    # EXTERNAL METRICS FIXTURE in points, not browser measurements.
    request = view.measurement_request()
    return view.measured(
        measurements={
            "id": request["id"],
            "metrics": [
                {"id": label["id"], "width": 24.0, "height": 12.0, "baseline": 9.0}
                for label in request["labels"]
            ],
        }
    )


model = hep.Model.phi3()
diagram = (
    model.process(["phi"], ["phi", "phi"])
    .generate_diagrams(loops=1, max_vertices=3, maximum_bridges=0, progress=None)
    .diagrams[0]
)
cff = diagram.integrate_energy(method="cff")
assert type(cff) is hep.CffRepresentation

# Exercise native-owner and prepared-view dispatch through the same widget adapter.
view = cff.render()
request = view.measurement_request()
snapshot = external_metrics_snapshot(view)
assert callable(cff._mime_) and callable(cff._repr_mimebundle_)
assert not hasattr(cff, "_repr_html_")
assert not hasattr(cff, "_repr_svg_")
for value in (cff, view):
    with patch.object(mo.ui, "anywidget", wraps=mo.ui.anywidget) as adapter:
        markup = mo.as_html(value).text
    adapter.assert_called_once()
    widget = adapter.call_args.args[0]
    try:
        assert "marimo-ui-element" in markup
        assert widget.request["labels"] == request["labels"]
        assert "measure" in widget._esm
        if value is view:
            html = widget.rendered["html"]
            # Marimo and static exports retain exactly the same native algebra.
            payloads = [
                json.loads(
                    re.search(
                        r'<script type="application/json" data-cff>(.*?)</script>',
                        source,
                    ).group(1)
                )
                for source in (html, snapshot.to_html())
            ]
            assert payloads[0] == payloads[1]
            assert len(payloads[0]["orientations"]) == len(cff.orientations)
    finally:
        widget.close()
print(
    "CffRepresentation: Marimo dispatch retains the shared widget and native explorer"
)
