"""External metrics fixtures for geometry tests, not browser measurements.

The browser roundtrip has its own harness. These fixed point dimensions keep
identity, styling, and snapshot tests independent of a browser or compiler.
"""


def external_metrics_snapshot(view):
    request = view.measurement_request()
    receipt = {
        "id": request["id"],
        "metrics": [
            {"id": label["id"], "width": 24.0, "height": 12.0, "baseline": 9.0}
            for label in request["labels"]
        ],
    }
    return view.measured(measurements=receipt)
