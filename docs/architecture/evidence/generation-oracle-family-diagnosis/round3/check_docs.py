"""Scoped registry checks from the actual DeveloperConfig/Note schema.

This does not bypass or claim success of the full docs-builder, whose earlier
source-format gate masks unrelated architecture classification failures.
"""

import re
import tomllib

import common as C


def main():
    path = C.WORKSPACE / "docs/developers.toml"
    data = tomllib.loads(path.read_text())
    C.require(type(data["schema"]) is int and data["schema"] == 3, "docs", "schema")
    notes = [n for s in data["section"] for n in s["note"]]
    for field in ("id", "source", "title"):
        values = [n[field] for n in notes]
        C.require(len(values) == len(set(values)), "docs", "duplicate " + field)
    ids = {n["id"] for n in notes}
    task = [
        n for n in notes if n["id"].startswith("generation-oracle-family-diagnosis")
    ]
    C.require(
        len(task) == 3 and sum(n["lifecycle"] != "superseded" for n in task) == 1,
        "docs",
        "one current task owner; two archived revisions",
    )
    owners = {o["id"] for o in data["owner"]}
    for n in task:
        for field in (
            "id",
            "title",
            "summary",
            "source",
            "kind",
            "status",
            "lifecycle",
            "owner",
        ):
            C.require(type(n[field]) is str and bool(n[field].strip()), "docs", field)
        C.require(
            n["lifecycle"]
            in ("current", "proposal", "investigation", "archived", "superseded"),
            "docs",
            n["id"],
        )
        C.require(
            n["owner"] in owners and bool(n["products"]) and bool(n["topics"]),
            "docs",
            n["id"],
        )
        C.require(("reviewed_at" in n) == ("review_ref" in n), "docs", "review pair")
        C.require((C.WORKSPACE / n["source"]).is_file(), "docs", n["source"])
        if n["lifecycle"] == "superseded":
            C.require(n["superseded_by"] in ids, "docs", "replacement")
    for n in task:
        title = re.search(
            r"^= (.+)$", (C.WORKSPACE / n["source"]).read_text(), re.MULTILINE
        )
        C.require(
            title and title[1] == n["title"], "docs", "exact task title " + n["id"]
        )
        C.require(
            n["lifecycle"] != "superseded"
            or n["superseded_by"] == "generation-oracle-family-diagnosis",
            "docs",
            "task replacement",
        )
    sources = {
        str(p.relative_to(C.WORKSPACE))
        for p in (C.WORKSPACE / "docs/architecture").glob("*.typ")
    }
    unclassified = sorted(sources - {n["source"] for n in notes})
    C.require(
        unclassified == ["docs/architecture/generation-integration-source46.typ"],
        "docs",
        "unexpected classification drift",
    )
    C.write(
        C.OUT / "docs/registry.json",
        {
            "task": task,
            "all_note_count": len(notes),
            "checks": "actual schema field kinds, unique IDs/sources/titles, exact task titles and supersession; scoped not full builder",
            "inherited_masked_unclassified": unclassified,
            "inherited_unknown_owner_references": [
                {"id": n["id"], "owner": n["owner"]}
                for n in notes
                if n["owner"] not in owners
            ],
            "registry_sha256": C.sha(path),
            "builder_source_sha256": C.sha(
                C.WORKSPACE / "crates/alphal00p-docs-builder/src/lib.rs"
            ),
        },
    )


if __name__ == "__main__":
    main()
