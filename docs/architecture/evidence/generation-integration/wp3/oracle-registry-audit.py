# Corrected ORACLE-P2-REGISTRY audit. The original presence-only audit is retained as history.
# Mirrors SiteBuilder::check uniqueness, range and git-blob digest rules; inventory is not a global docs pass.
from pathlib import Path
import json, tomllib, subprocess, hashlib

ROOT = Path.cwd()
OUT = ROOT / "target/oracle-p2-fix"
BASE = "438b61500da4b32841c64211103e5da6fc4c4a72"


def old(path):
    r = subprocess.run(
        ["jj", "--ignore-working-copy", "file", "show", "-r", BASE, path],
        capture_output=True,
    )
    return r.stdout.decode() if r.returncode == 0 else None


def scopes(reg):
    for section in reg["section"]:
        for note in section.get("note", []):
            for scope in note.get("scope", []):
                yield note, scope


def check(s, text):
    if text is None:
        return {"status": "missing-source"}
    if not s["anchor"].strip() or text.count(s["anchor"]) != 1:
        return {"status": "missing-or-ambiguous-anchor"}
    if "end_anchor" in s:
        if not s["end_anchor"].strip() or text.count(s["end_anchor"]) != 1:
            return {"status": "missing-or-ambiguous-end"}
        start = text.index(s["anchor"])
        end = text.index(s["end_anchor"])
        if end < start + len(s["anchor"]):
            return {"status": "invalid-range"}
        text = text[start:end]
    data = text.encode()
    actual = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
    return {
        "status": "valid" if actual == s["digest"] else "stale-digest",
        "actual": actual,
    }


def main():
    reg = tomllib.loads((ROOT / "docs/developers.toml").read_text())
    baseline = tomllib.loads(old("docs/developers.toml"))
    cache = {}
    prior = {}
    for note, s in scopes(baseline):
        cache.setdefault(s["path"], old(s["path"]))
        prior[
            (note["id"], s["path"], s.get("symbol"), s.get("anchor"), s.get("end_anchor"))
        ] = (s, check(s, cache[s["path"]]))
    rows = []
    for note, s in scopes(reg):
        p = ROOT / s["path"]
        result = check(s, p.read_text() if p.exists() else None)
        key = (note["id"], s["path"], s.get("symbol"), s.get("anchor"), s.get("end_anchor"))
        before = prior.get(key)
        inherited = bool(
            before and before[0]["digest"] == s["digest"] and before[1]["status"] != "valid"
        )
        rows.append(
            {
                "note": note["id"],
                **s,
                **result,
                "inherited_at_438": inherited,
                "baseline_check": before[1] if before else None,
            }
        )
    result = {
        "base": BASE,
        "scopes": rows,
        "valid": sum(r["status"] == "valid" for r in rows),
        "inherited_mismatches": [
            r for r in rows if r["status"] != "valid" and r["inherited_at_438"]
        ],
        "new_mismatches": [
            r for r in rows if r["status"] != "valid" and not r["inherited_at_438"]
        ],
        "claim": "Inventory only; not a full documentation pass. Global source31/final scope reconciliation remains pending.",
    }
    (OUT / "registry-audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(
        "Valid",
        result["valid"],
        "inherited mismatches",
        len(result["inherited_mismatches"]),
        "new mismatches",
        len(result["new_mismatches"]),
    )
    for r in result["new_mismatches"]:
        print(r["note"], r["path"], r["status"])


if __name__ == "__main__":
    main()
