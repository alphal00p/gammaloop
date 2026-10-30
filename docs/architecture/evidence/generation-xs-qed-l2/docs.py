"""One new canonical note; actual semantic PDF/HTML and inherited docs gate."""

import re
import tomllib
from html.parser import HTMLParser

import common as M

TYPST = "/nix/store/2szzba1xi60y8i3qcyxwvgq2i6w9iy5s-typst-0.15.1/bin/typst"


class SemanticHTML(HTMLParser):
    def __init__(self):
        super().__init__()
        self.title = []
        self.level = None
        self.buffer = []
        self.headers = 0
        self.active = []

    def handle_starttag(self, tag, attrs):
        if tag in ("h1", "h2", "h3", "h4", "h5", "h6"):
            self.level = tag
            self.buffer = []
        self.headers += tag == "thead"
        if tag == "script" or any(k.startswith("on") for k, _ in attrs):
            self.active.append(tag)

    def handle_data(self, data):
        if self.level:
            self.buffer.append(data)

    def handle_endtag(self, tag):
        if tag == self.level:
            self.title.append("".join(self.buffer))
            self.level = None


def main():
    registry = M.WORKSPACE / "docs/developers.toml"
    table = tomllib.loads(registry.read_text())
    notes = [n for section in table["section"] for n in section["note"]]
    M.require(table["schema"] == 3, "docs", "actualschema")
    for key in ("id", "source", "title"):
        values = [n[key] for n in notes]
        M.require(len(values) == len(set(values)), "docs", "unique" + key)
    entry = next(n for n in notes if n["id"] == "generation-xs-qed-l2")
    note = M.WORKSPACE / entry["source"]
    M.require(
        re.findall(r"^= (.+)$", note.read_text(), re.MULTILINE) == [entry["title"]],
        "docs",
        "canonicaltitle",
    )
    M.require(
        entry["owner"] == "unassigned"
        and entry["lifecycle"] == "investigation"
        and entry["products"] == ["feynkit", "gammaloop"],
        "docs",
        "actualregistryfields",
    )
    all_sources = {
        str(p.relative_to(M.WORKSPACE))
        for p in (M.WORKSPACE / "docs/architecture").glob("*.typ")
    }
    M.require(
        all_sources - {n["source"] for n in notes}
        == {"docs/architecture/generation-integration-source46.typ"},
        "docs",
        "no unregisterednewnote",
    )
    inputs = [
        note,
        registry,
        M.HERE / "receipt.json",
        M.HERE / "schema.json",
        M.HERE / "docs.py",
        M.OUT / "docs/main.typ",
    ]
    commands = []
    for label, argv in (
        (
            "pdf",
            [
                TYPST,
                "compile",
                "--root",
                ".",
                "--creation-timestamp",
                "0",
                note,
                M.OUT / "docs/report.pdf",
            ],
        ),
        (
            "html",
            [
                TYPST,
                "compile",
                "--root",
                ".",
                "--features",
                "html,bundle",
                "--format",
                "bundle",
                "--creation-timestamp",
                "0",
                M.OUT / "docs/main.typ",
                M.OUT / "docs/bundle",
            ],
        ),
        ("builder", [M.WORKSPACE / "target/debug/alphal00p-docs-builder", "check"]),
    ):
        p, record = M.run(argv, M.OUT / "docs" / label, inputs=inputs)
        if label == "builder":
            M.require(
                p.returncode == 1
                and "shrinking legacy inventory" in p.stderr
                and p.stderr.count(
                    '"docs/architecture/evidence/generation-integration/mark-geometry/'
                )
                == 8
                and '".local/review-family-r2-b-20261030T010814/commands.md"'
                in p.stderr,
                "docs",
                "inheritedfirstgate",
            )
            record["qualification"] = (
                "firstgate8inheritedHTML+retainedreviewcommands.md; masked unknownsunsetowners then source46; not fullgreen"
            )
        else:
            M.require(p.returncode == 0, "docs", label)
        commands.append(record)
    html = M.OUT / "docs/bundle/note.html"
    parser = SemanticHTML()
    parser.feed(html.read_text())
    M.require(
        parser.title[0] == entry["title"]
        and parser.title.count(entry["title"]) == 1
        and parser.headers == 3
        and not parser.active,
        "HTML",
        "exacttitle/semantictables",
    )
    M.write(
        M.OUT / "docs/validation.json",
        {
            "entry": entry,
            "commands": commands,
            "html_sha256": M.sha(html),
            "pdf_sha256": M.sha(M.OUT / "docs/report.pdf"),
            "semantic_headers": parser.headers,
            "actual_render_times": "in commandrecords;creationtimestamp0 deterministic metadata",
            "full_docs": "notgreen;inheritedfirstgateandmaskedlimits",
        },
    )


if __name__ == "__main__":
    main()
