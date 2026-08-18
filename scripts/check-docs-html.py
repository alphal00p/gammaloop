#!/usr/bin/env python3
"""Check deterministic accessibility invariants in every generated docs shell page."""

from __future__ import annotations

import re
import sys
from collections import Counter
from html.parser import HTMLParser
from pathlib import Path


VOID_ELEMENTS = {
    "area",
    "base",
    "br",
    "col",
    "embed",
    "hr",
    "img",
    "input",
    "link",
    "meta",
    "param",
    "source",
    "track",
    "wbr",
}


class DocumentAudit(HTMLParser):
    """Collect semantic facts while retaining enough ancestry for labels and navigation."""

    def __init__(self, path: Path) -> None:
        super().__init__(convert_charrefs=True)
        self.path = path
        self.errors: list[str] = []
        self.stack: list[dict[str, object]] = []
        self.ids: Counter[str] = Counter()
        self.label_targets: set[str] = set()
        self.labelled_controls: list[tuple[str, dict[str, str], bool]] = []
        self.references: list[tuple[str, str]] = []
        self.doctype = False
        self.html_language = ""
        self.html_classes: set[str] = set()
        self.main_count = 0
        self.h1_count = 0
        self.heading_levels: list[int] = []
        self.viewport = False
        self.skip_link = False
        self.sidebar_count = 0
        self.sidebar_links = 0
        self.primary_navigation_links = 0
        self.title = ""

    def handle_decl(self, declaration: str) -> None:
        self.doctype |= declaration.lower().strip() == "doctype html"

    def handle_starttag(self, tag: str, raw_attributes: list[tuple[str, str | None]]) -> None:
        attributes = {name.lower(): value or "" for name, value in raw_attributes}
        tag = tag.lower()
        ancestors = tuple(self.stack)

        element_id = attributes.get("id")
        if element_id:
            self.ids[element_id] += 1
        if tag == "html":
            self.html_language = attributes.get("lang", "")
            self.html_classes = set(attributes.get("class", "").split())
        elif tag == "meta" and attributes.get("name", "").lower() == "viewport":
            self.viewport = bool(attributes.get("content", "").strip())
        elif tag == "main":
            self.main_count += 1
            if attributes.get("id") != "main-content":
                self.errors.append("the main landmark must use id=\"main-content\"")
        elif tag in {"h1", "h2", "h3", "h4", "h5", "h6"}:
            level = int(tag[1])
            self.heading_levels.append(level)
            self.h1_count += level == 1
        elif tag == "label" and attributes.get("for"):
            self.label_targets.add(attributes["for"])
        elif tag in {"input", "select", "textarea"}:
            if tag != "input" or attributes.get("type", "").lower() != "hidden":
                inside_label = any(ancestor["tag"] == "label" for ancestor in ancestors)
                self.labelled_controls.append((tag, attributes, inside_label))
        elif tag == "img" and "alt" not in attributes:
            self.errors.append("an img element is missing its alt attribute")
        elif attributes.get("role") == "img" and not self._has_accessible_name(attributes):
            self.errors.append("an element with role=\"img\" has no accessible name")
        elif tag == "iframe" and not attributes.get("title", "").strip():
            self.errors.append("an iframe has no title")
        elif tag == "th":
            table = next(
                (ancestor for ancestor in reversed(ancestors) if ancestor["tag"] == "table"),
                None,
            )
            if table is not None:
                table["has_header"] = True

        if tag in {"nav", "aside"} and not self._has_accessible_name(attributes):
            article = "an" if tag == "aside" else "a"
            self.errors.append(
                f"{article} {tag} landmark has no aria-label or aria-labelledby"
            )

        href = attributes.get("href", "") if tag == "a" else ""
        if href.startswith("#") and len(href) > 1:
            self.references.append(("href", href[1:]))
        for attribute in ("aria-controls", "aria-labelledby"):
            for target in attributes.get(attribute, "").split():
                self.references.append((attribute, target))

        classes = set(attributes.get("class", "").split())
        in_sidebar = any(
            ancestor["attributes"].get("id") == "docs-sidebar" for ancestor in ancestors
        )
        in_primary_navigation = any(
            "portal-nav" in ancestor["attributes"].get("class", "").split()
            for ancestor in ancestors
        )
        if tag == "aside" and element_id == "docs-sidebar":
            self.sidebar_count += 1
            for forbidden in ("hidden", "inert"):
                if forbidden in attributes:
                    self.errors.append(
                        f"the documentation sidebar is {forbidden} without JavaScript"
                    )
            if attributes.get("aria-hidden", "").lower() == "true":
                self.errors.append("the documentation sidebar is aria-hidden without JavaScript")
        elif tag == "a" and in_sidebar and href:
            self.sidebar_links += 1
        elif tag == "a" and in_primary_navigation and href:
            self.primary_navigation_links += 1
        if tag == "a" and "skip-link" in classes and href == "#main-content":
            self.skip_link = True

        node: dict[str, object] = {
            "tag": tag,
            "attributes": attributes,
            "text": [],
            "has_header": False,
        }
        if tag not in VOID_ELEMENTS:
            self.stack.append(node)

    def handle_startendtag(
        self, tag: str, raw_attributes: list[tuple[str, str | None]]
    ) -> None:
        self.handle_starttag(tag, raw_attributes)
        if tag.lower() not in VOID_ELEMENTS:
            self.handle_endtag(tag)

    def handle_data(self, data: str) -> None:
        for ancestor in self.stack:
            if ancestor["tag"] in {"button", "title"}:
                ancestor["text"].append(data)

    def handle_endtag(self, tag: str) -> None:
        tag = tag.lower()
        matching = next(
            (
                index
                for index in range(len(self.stack) - 1, -1, -1)
                if self.stack[index]["tag"] == tag
            ),
            None,
        )
        if matching is None:
            return
        node = self.stack[matching]
        del self.stack[matching:]
        text = " ".join("".join(node["text"]).split())
        attributes = node["attributes"]
        if tag == "title":
            self.title = text
        elif tag == "button" and not text and not self._has_accessible_name(attributes):
            self.errors.append("a button has no text or accessible name")
        elif tag == "table" and not node["has_header"]:
            self.errors.append("a data table has no header cells")

    def finish(self, require_sidebar: bool) -> list[str]:
        if not self.doctype:
            self.errors.append("the HTML5 doctype is missing")
        if self.html_language != "en":
            self.errors.append("the root html element must declare lang=\"en\"")
        if "js" in self.html_classes:
            self.errors.append("the static html element must not enable JavaScript-only navigation")
        if not self.title:
            self.errors.append("the document title is empty")
        if not self.viewport:
            self.errors.append("the viewport meta element is missing or empty")
        if self.main_count != 1:
            self.errors.append(f"expected one main landmark, found {self.main_count}")
        if self.h1_count != 1:
            self.errors.append(f"expected one h1, found {self.h1_count}")
        for previous, current in zip(self.heading_levels, self.heading_levels[1:]):
            if current > previous + 1:
                self.errors.append(f"heading level jumps from h{previous} to h{current}")
        if not self.skip_link:
            self.errors.append("the skip link to #main-content is missing")

        for element_id, count in self.ids.items():
            if count > 1:
                self.errors.append(f'id="{element_id}" occurs {count} times')
        for attribute, target in self.references:
            if target not in self.ids:
                self.errors.append(f'{attribute} references missing id="{target}"')
        for tag, attributes, inside_label in self.labelled_controls:
            labelled_by = attributes.get("aria-labelledby", "").split()
            labelled = (
                inside_label
                or bool(attributes.get("aria-label", "").strip())
                or bool(labelled_by)
                or attributes.get("id") in self.label_targets
            )
            if not labelled:
                self.errors.append(f"a {tag} control has no associated accessible label")

        if require_sidebar:
            if self.sidebar_count != 1:
                self.errors.append(
                    f"expected one documentation sidebar, found {self.sidebar_count}"
                )
            if not self.sidebar_links:
                self.errors.append(
                    "the documentation sidebar has no ordinary links for no-JS navigation"
                )
        elif self.primary_navigation_links < 3:
            self.errors.append("the portal primary navigation has fewer than three ordinary links")
        return self.errors

    @staticmethod
    def _has_accessible_name(attributes: dict[str, str]) -> bool:
        return (
            attributes.get("aria-hidden", "").lower() == "true"
            or bool(attributes.get("aria-label", "").strip())
            or bool(attributes.get("aria-labelledby", "").strip())
        )


def generated_shell_pages(root: Path) -> list[tuple[str, Path, bool]]:
    """Discover non-redirect pages owned by the shared portal/manual shell."""

    pages: list[tuple[str, Path, bool]] = []
    stylesheet = re.compile(r'<link\b[^>]*\bhref="[^"]*assets/site\.css"')
    for path in sorted(root.rglob("*.html")):
        document = path.read_text(encoding="utf-8")
        if not stylesheet.search(document):
            continue
        relative = path.relative_to(root)
        require_sidebar = bool(relative.parts) and relative.parts[0] in {
            "developers",
            "products",
        }
        pages.append((relative.as_posix(), path, require_sidebar))
    return pages


def progressive_navigation_errors(root: Path) -> list[str]:
    """Confirm that off-canvas mobile navigation is an enhancement, not the default."""

    assets = {name: root / "assets" / name for name in ("site.css", "site.js")}
    missing = [path for path in assets.values() if not path.is_file()]
    if missing:
        return [f"missing navigation asset {path.relative_to(root)}" for path in missing]
    css = assets["site.css"].read_text(encoding="utf-8")
    script = assets["site.js"].read_text(encoding="utf-8")
    css_rules = {
        "the no-JS mobile sidebar is not in normal flow": (
            r"\.docs-sidebar\s*\{[^}]*position:\s*static"
        ),
        "the off-canvas sidebar is not scoped to the JavaScript marker": (
            r"\.js\s+\.docs-sidebar\s*\{[^}]*position:\s*fixed[^}]*"
            r"transform:\s*translateX\(-105%\)"
        ),
        "the mobile menu button is not scoped to the JavaScript marker": (
            r"\.js\s+\.menu-button\s*\{[^}]*display:\s*inline-block"
        ),
        "the enhanced mobile sidebar has no visible open state": (
            r"\.js\s+body\.sidebar-open\s+\.docs-sidebar\s*\{[^}]*"
            r"transform:\s*translateX\(0\)"
        ),
    }
    errors = [message for message, pattern in css_rules.items() if not re.search(pattern, css)]
    if not re.search(r"\.classList\.add\([\"']js[\"']\)", script):
        errors.append("site.js does not add the progressive-enhancement marker")
    return errors


def main() -> int:
    if len(sys.argv) != 2:
        print(f"usage: {Path(sys.argv[0]).name} BUILT_DOCS_ROOT", file=sys.stderr)
        return 2
    root = Path(sys.argv[1])
    failures: list[str] = []
    try:
        pages = generated_shell_pages(root)
    except (OSError, UnicodeError) as error:
        print(f"could not discover generated documentation pages: {error}", file=sys.stderr)
        return 1
    if not pages:
        print("no generated documentation shell pages found", file=sys.stderr)
        return 1
    try:
        failures.extend(
            f"navigation assets: {error}" for error in progressive_navigation_errors(root)
        )
    except (OSError, UnicodeError) as error:
        failures.append(f"navigation assets could not be read: {error}")
    for name, path, require_sidebar in pages:
        if not path.is_file():
            failures.append(f"{name}: missing {path.relative_to(root)}")
            continue
        audit = DocumentAudit(path)
        try:
            audit.feed(path.read_text(encoding="utf-8"))
            audit.close()
        except (OSError, UnicodeError) as error:
            failures.append(f"{name}: could not parse {path.relative_to(root)}: {error}")
            continue
        failures.extend(
            f"{name} ({path.relative_to(root)}): {error}"
            for error in audit.finish(require_sidebar)
        )

    if failures:
        print("documentation HTML accessibility checks failed:", file=sys.stderr)
        for failure in failures:
            print(f"- {failure}", file=sys.stderr)
        return 1
    print(
        f"checked all {len(pages)} generated HTML shell pages "
        "with deterministic semantic/no-JS rules"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
