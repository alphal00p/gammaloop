"""Marimo must embed CFF scripts in an iframe when formatting the native object."""

import json
import re
from html.parser import HTMLParser

import marimo as mo
from symbolica.community import hepkit as hep


class FrameParser(HTMLParser):
    def handle_starttag(self, tag, attrs):
        if tag == "iframe":
            self.frame = dict(attrs)


model = hep.Model.phi3()
diagram = (
    model.process(["phi"], ["phi", "phi"])
    .generate_diagrams(loops=1, max_vertices=3, maximum_bridges=0, progress=None)
    .diagrams[0]
)
cff = diagram.integrate_energy(method="cff")
assert type(cff) is hep.CffRepresentation

# Exercise Marimo's actual dispatch: a raw text/html MIME hook bypasses its
# script isolation and leaves the CFF controls visible but inactive.
parser = FrameParser()
parser.feed(mo.as_html(cff).text)
html = parser.frame["srcdoc"]
assert "allow-scripts" in parser.frame.get("sandbox", "allow-scripts")
# Marimo normalizes HTML whitespace; compare the native data after parsing it.
payloads = [
    json.loads(
        re.search(
            r'<script type="application/json" data-cff>(.*?)</script>', source
        ).group(1)
    )
    for source in (html, cff._repr_html_())
]
assert payloads[0] == payloads[1]
assert len(payloads[0]["orientations"]) == len(cff.orientations)
print("CffRepresentation: Marimo dispatch embeds the native explorer in an iframe")
