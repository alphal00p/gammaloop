"""Collection receipts complete their original children in source order."""
# ruff: noqa: F821 -- scenes, first, second and collection are injected by Rust.

import copy


def receipt(view):
    request = view.measurement_request()
    return {
        "id": request["id"],
        "metrics": [
            {
                "id": label["id"],
                "width": 30.0 + label["id"],
                "height": 18.0,
                "baseline": 12.0,
            }
            for label in request["labels"]
        ],
    }


request = collection.measurement_request()
assert len(request["labels"]) == 12
assert len(scenes.measurement_request()["labels"]) == 12
assert len(collection.diagrams) == 2
assert not scenes.diagrams
assert "child-a" in request["labels"][0]["mathml"]
assert "child-b" in request["labels"][6]["mathml"]
assert (
    len(
        {
            view.measurement_request()["id"]
            for view in (collection, first, second, scenes)
        }
    )
    == 4
)
measurements = receipt(collection)
snapshot = collection.measured(measurements=measurements)
assert len(snapshot.to_svg_pages()) == len(snapshot.diagrams) == 2
assert first.to_svg() == snapshot.diagrams[0].to_svg()
assert second.to_svg() == snapshot.diagrams[1].to_svg()
html = snapshot.to_html()
assert "<section data-collection='native'>" in html
assert html.count("data:font/woff2;base64,") == 1
assert html.index("Native child-b") < html.index("Native child-a")
assert "<!--linnet-svg-" not in html
assert collection.to_html() == html
assert first.measurement_request()["id"] != second.measurement_request()["id"]
bad = copy.deepcopy(measurements)
bad["id"] = first.measurement_request()["id"]
try:
    collection.measured(measurements=bad)
except ValueError as error:
    assert "stale" in str(error)
else:
    raise AssertionError("a child receipt completed its parent")
channels = scenes.measured(measurements=receipt(scenes))
assert len(channels.to_svg_pages()) == 1
assert "Native channel-a" in channels.to_svg()
assert "Native channel-b" in channels.to_svg()
assert channels.to_svg().count('data-linnet-font="true"') == 1
