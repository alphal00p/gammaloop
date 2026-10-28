"""Scope real generated metadata to its unchanged staged CLI owner."""

import hashlib
import json
from pathlib import Path
import subprocess

BASE = "d48f9d1aedc766eea5cb185ab6029d1f04462043"
ROOT = Path("docs/architecture/evidence/generation-projector-flow-staged")
SOURCE = "crates/gammaloop-api/src/commands/generate.rs"


def baseline(path):
    return subprocess.check_output(["jj", "file", "show", "-r", BASE, path]).decode()


def spec_args(source):
    return source[
        source.index("pub struct SpecArgs {") : source.index("impl SpecArgs {")
    ]


before = spec_args(baseline(SOURCE))
after = spec_args(Path(SOURCE).read_text())
new_doc = "/// Global prefactor projector (Symbolica atom string). The k-th model index\n    /// of half-edge h is `hedge(h, k)`, counted from one"
old_doc = "/// Global prefactor projector (Symbolica atom string)"
assert new_doc in after
assert after.replace(new_doc, old_doc) == before


def differences(old, new, path=()):
    if isinstance(old, dict) and isinstance(new, dict):
        result = []
        for key in sorted(old.keys() | new.keys()):
            if key not in old or key not in new:
                result.append(
                    {"path": [*path, key], "old": old.get(key), "new": new.get(key)}
                )
            else:
                result.extend(differences(old[key], new[key], (*path, key)))
        return result
    if isinstance(old, list) and isinstance(new, list) and len(old) == len(new):
        return [
            item
            for i, (a, b) in enumerate(zip(old, new))
            for item in differences(a, b, (*path, i))
        ]
    return [] if old == new else [{"path": list(path), "old": old, "new": new}]


records = []
for path in [
    "assets/schemas/runhistory.json",
    "docs/api/generated/gammaloop-reference.json",
]:
    old = json.loads(baseline(path))
    new = json.loads(Path(path).read_text())
    changes = differences(old, new)
    assert changes
    for change in changes:
        owned = isinstance(change["new"], str) and change["new"].startswith(
            "Global prefactor projector"
        )
        if path.startswith("docs/api/"):
            assert owned and change["path"][-1] == "help"
        else:
            assert "SpecArgs" in change["path"] and "properties" in change["path"]
            if owned:
                assert change["path"][-2:] == [
                    "global_prefactor_projector",
                    "description",
                ]
        change["scope"] = (
            "source24 projector description"
            if owned
            else "unchanged approved d48 SpecArgs baseline metadata alignment"
        )
    records.append(
        {
            "path": path,
            "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
            "changes": changes,
        }
    )
(ROOT / "generated-metadata.json").write_text(
    json.dumps(
        {
            "base": BASE,
            "staged_spec_args_source_unchanged_except_projector_doc": True,
            "baseline_spec_args_sha256": hashlib.sha256(before.encode()).hexdigest(),
            "outputs": records,
        },
        indent=2,
        ensure_ascii=False,
    )
    + "\n"
)
print(
    "real generated output: projector help metadata plus schema alignment to unchanged approved SpecArgs"
)
