#!/usr/bin/env python3
"""Check application test-source isolation and build-script change tracking."""

import argparse
import json
import re
import tempfile
from pathlib import Path


class RustSource:
    """Keep source offsets while hiding comments and literal contents."""

    raw_string = re.compile(r'r(\#*)"')
    character = re.compile(
        r"'(?:\\(?:x[0-9a-fA-F]{2}|u\{[0-9a-fA-F_]+\}|[^\n])|[^'\\\n])'"
    )

    def __init__(self, text):
        self.text = text
        self.literals = []
        code = list(text)
        index = 0
        while index < len(text):
            start = index
            if text.startswith("//", index):
                end = text.find("\n", index)
                index = len(text) if end < 0 else end
            elif text.startswith("/*", index):
                depth = 1
                index += 2
                while index < len(text) and depth:
                    if text.startswith("/*", index):
                        depth += 1
                        index += 2
                    elif text.startswith("*/", index):
                        depth -= 1
                        index += 2
                    else:
                        index += 1
            elif match := self.raw_string.match(text, index):
                closing = '"' + match[1]
                end = text.find(closing, index + len(match[0]))
                if end < 0:
                    raise ValueError("unterminated raw string")
                index = end + len(closing)
                self.literals.append(text[start:index])
            elif text[index] == '"':
                index += 1
                while index < len(text):
                    if text[index] == "\\":
                        index += 2
                    elif text[index] == '"':
                        index += 1
                        break
                    else:
                        index += 1
                self.literals.append(text[start:index])
            elif match := self.character.match(text, index):
                index += len(match[0])
            else:
                index += 1
                continue
            for offset in range(start, index):
                if code[offset] != "\n":
                    code[offset] = " "
        self.code = "".join(code)

    def block_end(self, opening):
        depth = 0
        for index in range(opening, len(self.code)):
            if self.code[index] == "{":
                depth += 1
            elif self.code[index] == "}":
                depth -= 1
                if depth == 0:
                    return index
        raise ValueError("unterminated Rust block")

    def external_tests(self, path, in_test_source=False):
        declarations = re.compile(
            r"(?P<attrs>(?:#\s*\[[^\]]*\]\s*)*)"
            r"(?:pub(?:\s*\([^)]*\))?\s+)?mod\s+(?P<name>\w+)\s*;"
        )
        for match in declarations.finditer(self.code):
            masked_attributes = self.code[match.start("attrs") : match.end("attrs")]
            if not in_test_source and not re.search(
                r"#\s*\[\s*cfg\s*\(\s*test\s*\)\s*\]", masked_attributes
            ):
                continue
            attributes = self.text[match.start("attrs") : match.end("attrs")]
            explicit = re.search(r'#\s*\[\s*path\s*=\s*"([^"\n]+)"\s*\]', attributes)
            if explicit:
                yield path.parent / explicit[1]
            else:
                parent = (
                    path.parent
                    if path.stem in ("lib", "main", "mod")
                    else path.with_suffix("")
                )
                flat = parent / f"{match['name']}.rs"
                directory = parent / match["name"] / "mod.rs"
                yield flat if flat.exists() else directory


def violations(root, configuration):
    declared = [
        root / entry
        for entries in configuration["ownTest"].values()
        for entry in entries
    ]

    def covered(path):
        return any(path == entry or entry in path.parents for entry in declared)

    errors = []
    for package in ("gammalooprs", "gammaloop-api"):
        for path in (root / "crates" / package / "src").rglob("*.rs"):
            for test in RustSource(path.read_text()).external_tests(
                path, covered(path)
            ):
                test = test.resolve()
                if not test.is_relative_to(root):
                    errors.append(
                        f"{path.relative_to(root)}: test source is outside the repository"
                    )
                    continue
                if not test.exists():
                    errors.append(
                        f"{path.relative_to(root)}: test module does not resolve: {test}"
                    )
                elif not covered(test):
                    errors.append(
                        f"{test.relative_to(root)}: add test-only source to ownTest"
                    )

    normalized = set(configuration["normalizedBuildScripts"])
    for relative in configuration["buildScripts"]:
        source = RustSource((root / relative).read_text())
        # Asset::generate emits file/env watches for its declared data inputs.
        tracked = any(
            "cargo:rerun-if-" in literal or "cargo::rerun-if-" in literal
            for literal in source.literals
        )
        tracked |= "symbolica_data::Asset::generate" in source.code
        tracked |= bool(
            re.search(r"\buse\s+symbolica_data::Asset\s*;", source.code)
            and re.search(r"\bAsset::generate\s*\(", source.code)
        )
        if not tracked and relative not in normalized:
            errors.append(
                f"{relative}: declare build-script rerun inputs or normalize its fingerprint"
            )
    return errors


def self_test():
    source = RustSource("""
// #[cfg(test)] mod absent;
/* outer /* nested */ { } */
fn sample<'a>(input: &'a str) { let _ = '{'; let _ = r###" } \" { "###; }
#[cfg(test)]
#[path = "other_tests.rs"]
mod tests;
""")
    opening = source.code.index("{")
    assert source.text[source.block_end(opening)] == "}"
    assert (
        source.text[source.block_end(opening) + 1 :].lstrip().startswith("#[cfg(test)]")
    )
    assert list(source.external_tests(Path("src/example.rs"))) == [
        Path("src/other_tests.rs")
    ]
    assert list(
        RustSource('#[path = "other_tests.rs"] #[cfg(test)] mod tests;').external_tests(
            Path("src/example.rs")
        )
    ) == [Path("src/other_tests.rs")]
    assert not RustSource('// println!("cargo:rerun-if-changed=build.rs");').literals
    assert any(
        "cargo:rerun-if-changed=" in value
        for value in RustSource('println!("cargo:rerun-if-changed=build.rs");').literals
    )
    with tempfile.TemporaryDirectory(prefix="gammaloop-source-lint-") as temporary:
        root = Path(temporary).resolve()
        package = root / "crates/gammalooprs/src"
        package.mkdir(parents=True)
        (package / "lib.rs").write_text("#[cfg(test)] mod tests;")
        (package / "tests.rs").write_text("#[test] fn example() {}")
        assert list(
            RustSource("#[cfg(test)] mod tests;").external_tests(package / "lib.rs")
        ) == [package / "tests.rs"]
        build = root / "build.rs"
        build.write_text('// println!("cargo:rerun-if-changed=build.rs");')
        configuration = {
            "ownTest": {},
            "buildScripts": ["build.rs"],
            "normalizedBuildScripts": [],
        }
        assert len(violations(root, configuration)) == 2
        configuration["ownTest"] = {"gammalooprs": ["crates/gammalooprs/src/tests.rs"]}
        build.write_text('fn main() { println!("cargo:rerun-if-changed=build.rs"); }')
        assert not violations(root, configuration)
        build.write_text(
            "use symbolica_data::Asset;\n"
            "fn main() { Asset::generate(todo!(), todo!(), todo!()); }"
        )
        assert not violations(root, configuration)
        build.write_text(
            "// use symbolica_data::Asset;\n"
            "fn main() { Asset::generate(todo!(), todo!(), todo!()); }"
        )
        assert len(violations(root, configuration)) == 1
        build.write_text('fn main() { println!("cargo:rerun-if-changed=build.rs"); }')
        (package / "tests.rs").write_text('#[path = "support.rs"] mod support;')
        (package / "support.rs").write_text("fn fixture() {}")
        assert len(violations(root, configuration)) == 1
        configuration["ownTest"]["gammalooprs"].append(
            "crates/gammalooprs/src/support.rs"
        )
        assert not violations(root, configuration)
        build.write_text("fn main() {}")
        configuration["normalizedBuildScripts"] = ["build.rs"]
        assert not violations(root, configuration)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("configuration", type=Path)
    parser.add_argument("root", type=Path)
    arguments = parser.parse_args()
    self_test()
    errors = violations(
        arguments.root.resolve(), json.loads(arguments.configuration.read_text())
    )
    if errors:
        parser.exit(1, "\n".join(errors) + "\n")
    print("Rust test-source boundaries and build-script declarations match")


if __name__ == "__main__":
    main()
