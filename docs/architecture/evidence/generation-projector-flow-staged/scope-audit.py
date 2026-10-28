"""Audit source22–26 scopes through the existing canonical-rule checker.

The original f2e32c44 presence-only checker and false breakdown remain historical
at that immutable revision. Default output is ignored integration evidence.
"""

import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import tomllib
from collections import Counter

OWNER = Path("docs/architecture/evidence/generation-integration/wp3/oracle-registry-audit.py")
spec = importlib.util.spec_from_file_location("canonical_scope_owner", OWNER)
owner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(owner)


def read(revision, path):
    if revision is None:
        p = Path(path)
        return p.read_text() if p.exists() else None
    result = subprocess.run(
        ["jj", "--ignore-working-copy", "file", "show", "-r", revision, path],
        capture_output=True,
        check=False,
    )
    return result.stdout.decode() if result.returncode == 0 else None


def key(note, scope):
    return (note["id"], scope["path"], scope.get("symbol"), scope["anchor"], scope.get("end_anchor"))


def audit(revision, base):
    manifest = tomllib.loads(read(revision, "docs/developers.toml"))
    previous = tomllib.loads(read(base, "docs/developers.toml"))
    prior = {key(n, s): s for n, s in owner.scopes(previous)}
    rows = []
    for note, scope in owner.scopes(manifest):
        result = owner.check(scope, read(revision, scope["path"]))
        old = prior.get(key(note, scope))
        before = owner.check(old, read(base, scope["path"])) if old else None
        if result["status"] == "valid":
            origin = "valid"
        elif old and old["digest"] == scope["digest"] and before == result:
            origin = "inherited_exact"
        elif old and old["digest"] == scope["digest"] and before["status"] != "valid":
            origin = "preexisting_failure_changed_source"
        else:
            origin = "new_mismatch"
        rows.append({"note": note["id"], **scope, **result, "origin": origin, "base_check": before})
    return {
        "revision": revision or "current working snapshot",
        "base": base,
        "by_type": dict(Counter(r["status"] for r in rows)),
        "by_origin": dict(Counter(r["origin"] for r in rows)),
        "rows": rows,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("target/wp3-units22-26/scope-audit.json"))
    args = parser.parse_args()
    result = {
        "canonical_owner": str(OWNER),
        "rules": "Exactly one nonblank start even without end; unique nonblank end when present; ordered nonoverlapping range; Git-blob hash. Whole source hashed when end absent.",
        "preparation": audit("f2e32c44c986058ede1d7e585b9e2e4fefe1110c", "d48f9d1aedc766eea5cb185ab6029d1f04462043"),
        "integration": audit(None, "0a38ee4e3e7d16af3978f2c6ccde04dc176dfea9"),
        "claim": "Per-scope status and origin, not counts of files/notes. No full documentation pass; no registry digest refreshed.",
    }
    assert result["preparation"]["by_type"].get("stale-digest") == 23
    assert sum(v for k, v in result["preparation"]["by_type"].items() if k not in ["valid", "stale-digest"]) == 4
    assert result["preparation"]["by_origin"].get("inherited_exact") == 27
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    for name in ["preparation", "integration"]:
        print(name, result[name]["by_type"], result[name]["by_origin"])


if __name__ == "__main__":
    main()
