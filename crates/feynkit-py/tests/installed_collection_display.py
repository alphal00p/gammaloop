"""Collection displays retain real diagrams and their aligned tensor terms."""

import importlib
import sys
from html.parser import HTMLParser
from pathlib import Path

fk = importlib.import_module(
    f"symbolica.community.{sys.argv[1] if len(sys.argv) > 1 else 'feynkit'}"
)


class Collection(HTMLParser):
    def __init__(self, source):
        super().__init__()
        self.templates = 0
        self.rows = 0
        self.buttons = 0
        self.feed(source)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        self.templates += tag == "template"
        self.rows += tag == "details" and attrs.get("class") == "fk-row"
        self.buttons += tag == "button" and "aria-pressed" in attrs


model = fk.Model(Path(__file__).parent / "fixtures/scalars_2p_3p.json")
process = model.process(
    ["scalar_0"], ["scalar_0", "scalar_0"], vertex_allow=["V_3_SCALAR_000"]
)
amplitude = process.generate_amplitude(loops=0, max_vertices=1, progress=None)
before = amplitude.expression()
html = amplitude._repr_html_()
parsed = Collection(html)
assert parsed.rows == parsed.templates == len(amplitude.diagrams) == 1
assert parsed.buttons == 0
assert amplitude.terms[0]._repr_html_() in html
assert amplitude.expression() == before
assert "Conjugate amplitude" in amplitude.conjugate()._repr_html_()

large = fk.Amplitude(amplitude.diagrams * 7)
large_html = large._repr_html_()
assert Collection(large_html).rows == 6
assert "Showing 6 of 7 diagrams" in large_html

cross_section = process.generate_cross_section(
    loops=1, max_vertices=2, allow_self_loops=True, progress=None
)
snapshots = [diagram.to_json() for diagram in cross_section.diagrams]
html = cross_section._repr_html_()
parsed = Collection(html)
assert parsed.buttons == parsed.templates == len(cross_section.diagrams) > 0
assert parsed.rows == 0
assert "<strong>Cross section</strong>" in html
assert all(diagram.cuts for diagram in cross_section.diagrams)
assert [diagram.to_json() for diagram in cross_section.diagrams] == snapshots

empty = process.generate_diagrams(
    loops=0,
    max_vertices=1,
    filter=lambda diagram, completed_vertices: False,
    progress=None,
)
html = empty._repr_html_()
assert Collection(html).templates == 0
assert "No diagrams retained." in html
assert "<strong>Generation result</strong>" in html
print(
    "Collection displays: terms, conjugation, preview bounds, cuts, empty results passed"
)
