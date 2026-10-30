"""Pinned report execution records with actual inputs/tool hashes before/after."""

from html.parser import HTMLParser

import common as C

TYPST = "/nix/store/2szzba1xi60y8i3qcyxwvgq2i6w9iy5s-typst-0.15.1/bin/typst"
RUFF = "/nix/store/0zh3g7fz7h3nh34cnzda1dbvg0cgwfsy-ruff-0.16.6/bin/ruff"


class SemanticHTML(HTMLParser):
    def __init__(self):
        super().__init__()
        self.headings, self.heading, self.text = [], None, []
        self.headers, self.active = 0, []

    def handle_starttag(self, tag, attrs):
        if tag in ("h1", "h2", "h3", "h4", "h5", "h6"):
            self.heading, self.text = tag, []
        self.headers += tag == "thead"
        if tag == "script" or any(k.startswith("on") for k, _ in attrs):
            self.active.append((tag, attrs))

    def handle_data(self, data):
        if self.heading:
            self.text.append(data)

    def handle_endtag(self, tag):
        if self.heading == tag:
            self.headings.append("".join(self.text))
            self.heading = None


def main():
    notes = [
        ("round1", "generation-oracle-family-diagnosis.typ"),
        ("round2", "generation-oracle-family-diagnosis-round2.typ"),
        ("round3", "generation-oracle-family-diagnosis-round3.typ"),
    ]
    inputs = [C.WORKSPACE / "docs/architecture" / file for _, file in notes]
    inputs += [
        C.HERE / "report_checks.py",
        C.HERE / "receipt.json",
        C.HERE / "dispositions.json",
        C.WORKSPACE / "docs/developers.toml",
        C.OUT / "docs/main.typ",
    ]
    before = {str(p): C.sha(p) for p in inputs}
    commands = []
    for label, file in notes:
        commands.append(
            (
                "pdf-" + label,
                [
                    TYPST,
                    "compile",
                    "--root",
                    ".",
                    "--creation-timestamp",
                    "0",
                    "docs/architecture/" + file,
                    C.OUT / "docs" / (label + ".pdf"),
                ],
            )
        )
    commands += [
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
                C.OUT / "docs/main.typ",
                C.OUT / "docs/bundle",
            ],
        ),
        ("ruff", [RUFF, "check", "--no-cache", C.HERE]),
        ("format", [RUFF, "format", "--no-cache", "--check", C.HERE]),
        (
            "full-builder",
            [C.WORKSPACE / "target/debug/alphal00p-docs-builder", "check"],
        ),
    ]
    rows = []
    for label, argv in commands:
        binary_before = C.sha(argv[0])
        proc, record = C.run(
            argv,
            C.OUT / "docs/commands" / label,
            bindings={"tool_sha256": binary_before, "inputs": before},
        )
        C.require(
            C.sha(argv[0]) == binary_before
            and before == {str(p): C.sha(p) for p in inputs},
            "report-input-changed",
            label,
        )
        if label == "full-builder":
            C.require(
                proc.returncode == 1
                and "shrinking legacy inventory" in proc.stderr
                and proc.stderr.count(
                    '"docs/architecture/evidence/generation-integration/mark-geometry/'
                )
                == 8
                and '".local/review-family-r2-b-20261030T010814/commands.md"'
                in proc.stderr,
                "docs-first-gate",
                proc.stderr,
            )
            record["qualification"] = (
                "8 inherited tracked HTML paths + preserved reviewer-owned commands.md; first gate only; source46/unknown-owner gates masked"
            )
        else:
            C.require(proc.returncode == 0, "report-check", label)
        record["tool_sha256_after"] = C.sha(argv[0])
        rows.append(record)
    semantic = []
    for label, file in notes:
        path = C.OUT / "docs/bundle" / (label + ".html")
        parser = SemanticHTML()
        parser.feed(path.read_text())
        title = (
            (C.WORKSPACE / "docs/architecture" / file).read_text().splitlines()[0][2:]
        )
        C.require(
            parser.headings
            and parser.headings[0] == title
            and parser.headings.count(title) == 1,
            "semantic-html",
            label + ": title",
        )
        C.require(not parser.active, "semantic-html", label + ": active content")
        if label == "round3":
            C.require(
                parser.headers >= 2, "semantic-html", "new tables need semantic headers"
            )
        semantic.append(
            {
                "note": label,
                "html_sha256": C.sha(path),
                "pdf_sha256": C.sha(C.OUT / "docs" / (label + ".pdf")),
                "first_heading": parser.headings[0],
                "semantic_table_headers": parser.headers,
            }
        )
    C.write(
        C.OUT / "docs/render-validation.json",
        {"inputs_before_after": before, "commands": rows, "semantic": semantic},
    )


if __name__ == "__main__":
    main()
