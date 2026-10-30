"""Scoped canonical classification and pinned rendered-artifact bindings."""

import re
import tomllib
from html.parser import HTMLParser

import common as M

TYPST = "/nix/store/2szzba1xi60y8i3qcyxwvgq2i6w9iy5s-typst-0.15.1/bin/typst"


class HTML(HTMLParser):
    def __init__(self):
        super().__init__()
        self.headings, self.heading, self.text, self.headers, self.active = (
            [],
            None,
            [],
            0,
            [],
        )

    def handle_starttag(self, tag, attrs):
        if tag in ("h1", "h2", "h3", "h4", "h5", "h6"):
            self.heading, self.text = tag, []
        self.headers += tag == "thead"
        if tag == "script" or any(k.startswith("on") for k, _ in attrs):
            self.active.append(tag)

    def handle_data(self, data):
        if self.heading:
            self.text.append(data)

    def handle_endtag(self, tag):
        if self.heading == tag:
            self.headings.append("".join(self.text))
            self.heading = None


def main():
    registry = M.WORKSPACE / "docs/developers.toml"
    data = tomllib.loads(registry.read_text())
    notes = [n for section in data["section"] for n in section["note"]]
    M.require(data["schema"] == 3, "docs", "actual registry schema")
    for field in ("id", "source", "title"):
        items = [n[field] for n in notes]
        M.require(len(set(items)) == len(items), "docs", "unique " + field)
    selected = [n for n in notes if n["id"].startswith("generation-oracle-family-")]
    M.require(
        len(selected) == 4,
        "docs",
        "three certificate classifications plus separate migration",
    )
    for n in selected:
        source = (M.WORKSPACE / n["source"]).read_text()
        heads = re.findall(r"^= (.+)$", source, re.MULTILINE)
        M.require(heads == [n["title"]], "docs", "exact single title")
        M.require(
            n["owner"] == "unassigned"
            and n["lifecycle"] in ("investigation", "superseded")
            and n["products"] == ["feynkit", "gammaloop"],
            "docs",
            n["id"],
        )
    sources = {
        str(p.relative_to(M.WORKSPACE))
        for p in (M.WORKSPACE / "docs/architecture").glob("*.typ")
    }
    M.require(
        sources - {n["source"] for n in notes}
        == {"docs/architecture/generation-integration-source46.typ"},
        "docs",
        "inherited classification only",
    )
    inputs = [
        registry,
        M.HERE / "receipt.json",
        M.HERE / "docs.py",
        M.OUT / "docs/main.typ",
        *[M.WORKSPACE / n["source"] for n in selected],
    ]
    records = []
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
                "docs/architecture/generation-oracle-family-migration.typ",
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
        proc, command = M.run(argv, M.OUT / "docs" / label, inputs=inputs)
        if label == "builder":
            M.require(
                proc.returncode == 1
                and "shrinking legacy inventory" in proc.stderr
                and proc.stderr.count(
                    '"docs/architecture/evidence/generation-integration/mark-geometry/'
                )
                == 8
                and '".local/review-family-r2-b-20261030T010814/commands.md"'
                in proc.stderr,
                "docs",
                "actual inherited first gate",
            )
            command["qualification"] = (
                "first gate only; masked two sunset unknown owners then source46; no full docs pass"
            )
        else:
            M.require(proc.returncode == 0, "docs", label)
        records.append(command)
    outputs = []
    for n in selected:
        label = {
            "generation-oracle-family-diagnosis": "round3",
            "generation-oracle-family-diagnosis-round1": "round1",
            "generation-oracle-family-diagnosis-round2": "round2",
            "generation-oracle-family-migration": "migration",
        }[n["id"]]
        path = M.OUT / "docs/bundle" / (label + ".html")
        parser = HTML()
        parser.feed(path.read_text())
        M.require(
            parser.headings[0] == n["title"]
            and parser.headings.count(n["title"]) == 1
            and not parser.active,
            "html",
            label,
        )
        if label == "migration":
            M.require(parser.headers == 3, "html", "semantic headers")
        outputs.append(
            {"note": n["id"], "sha256": M.sha(path), "headers": parser.headers}
        )
    M.write(
        M.OUT / "docs/validation.json",
        {
            "registry_notes": selected,
            "outputs": outputs,
            "commands": records,
            "report_pdf_sha256": M.sha(M.OUT / "docs/report.pdf"),
            "masked_inherited": ["unknown sunset owners", "unregistered source46"],
        },
    )


if __name__ == "__main__":
    main()
