"""Exact four-literal patch and immutable donor/outer-guard qualification."""

import difflib
import re

import common as M


def main():
    source, _ = M.immutable(M.APPROVED, M.MAIN)
    current = (M.WORKSPACE / M.MAIN).read_text()
    expected = source
    rows = M.load(M.OUT / "derivation.json")["rows"]
    for row in rows:
        pattern = (
            r"(\b" + row["case"] + r": Row \{[^\n]*?numerators: )" + row["old_literal"]
        )
        expected, count = re.subn(
            pattern, lambda match, value=row["new_literal"]: match[1] + value, expected
        )
        M.require(count == 1, "diff", row["case"])
    M.require(current == expected, "diff", "only four numerator literals")
    patch = "".join(
        difflib.unified_diff(
            source.splitlines(True),
            current.splitlines(True),
            fromfile="a/" + M.MAIN,
            tofile="b/" + M.MAIN,
        )
    )
    (M.OUT / "four-literal.patch").write_text(patch)
    changed = [
        i + 1
        for i, (a, b) in enumerate(zip(source.splitlines(), current.splitlines()))
        if a != b
    ]
    M.require(
        len(changed) == 4 and len(source.splitlines()) == len(current.splitlines()),
        "diff",
        "four changed lines",
    )
    old_registry, _ = M.immutable(M.APPROVED, "docs/developers.toml")
    registry = (M.WORKSPACE / "docs/developers.toml").read_text()
    marker = '[[section.note]]\nid = "generation-oracle-family-migration"'
    start = registry.index(marker)
    end = registry.index("[[section.note]]", start + len(marker))
    M.require(
        registry[:start] + registry[end:] == old_registry,
        "registry",
        "insert one migration entry only; three original classifications untouched",
    )
    (M.OUT / "migration-registry-entry.toml").write_text(registry[start:end])
    manifest = M.load(
        M.WORKSPACE
        / "target/oracle-numerator-family-round3/inventory/closure-final.json"
    )["files"]
    differences = []
    for path, data in manifest.items():
        if M.sha(path) != data["sha256"]:
            differences.append(path)
    expected_difference = str(M.WORKSPACE / "docs/developers.toml")
    M.require(differences == [expected_difference], "outer-guard", str(differences))
    M.require(
        M.C.hashlib.sha256(old_registry.encode()).hexdigest()
        == manifest[expected_difference]["sha256"],
        "immutable-donor",
        "approved registry blob",
    )
    M.write(
        M.OUT / "audit.json",
        {
            "approved": M.APPROVED,
            "changed_code_paths": [M.MAIN],
            "changed_lines": changed,
            "patch_sha256": M.sha(M.OUT / "four-literal.patch"),
            "new_main_sha256": M.sha(M.WORKSPACE / M.MAIN),
            "outer_manifest": {
                "verdict": "EXPECTED_FAILURE",
                "mismatches": differences,
                "credit": 0,
                "qualification": "original guard not disabled; approved donor registry verified by immutable commit; portable semantic/capture/control replay is separate",
            },
            "old_main_sha256": M.C.hashlib.sha256(source.encode()).hexdigest(),
            "three_xs_rows": "byte-preserved at THIS approved preparation base; MAIN values deliberately not imported",
            "reconciliation": "apply only four-literal.patch to freshly qualified MAIN; no main.rs transplant",
        },
    )
    print("four-line patch", changed, "outer guard expected failure", differences)


if __name__ == "__main__":
    main()
