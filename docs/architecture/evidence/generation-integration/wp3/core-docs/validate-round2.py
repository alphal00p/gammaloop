import copy
import hashlib
import json
import pathlib
import subprocess
import sys
import tomllib

ROOT = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
CORRECTION = ROOT / "docs/architecture/evidence/generation-integration/wp3/core-docs/round2-evidence-correction.json"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def text(path):
    return path.read_text(encoding="utf-8").strip()


def exit_value(path):
    return int(text(path))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def manifest_ok(path, expected_entries, revision=None):
    lines = path.read_text(encoding="utf-8").splitlines()
    require(len(lines) == expected_entries, f"manifest entries {path}")
    for line in lines:
        digest, relative = line.split("  ", 1)
        if revision is not None and not relative.startswith("target/"):
            blob = subprocess.run(
                ["jj", "--ignore-working-copy", "file", "show", "-r", revision, relative],
                cwd=ROOT,
                check=True,
                capture_output=True,
            ).stdout
            actual = hashlib.sha256(blob).hexdigest()
        else:
            actual = sha(ROOT / relative)
        require(actual == digest, f"manifest hash {relative}")


def diff_paths(old, new):
    result = subprocess.run(
        ["jj", "--ignore-working-copy", "diff", "--from", old, "--to", new, "--summary"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return {line[2:] for line in result.stdout.splitlines() if len(line) > 2}


def validate(model, require_current=True):
    frozen = model["frozen_round1_evidence"]
    for value in frozen.values():
        path = ROOT / value["path"]
        require(sha(path) == value["sha256"], f"frozen hash {path}")
    manifest_ok(ROOT / frozen["runtime_manifest"]["path"], frozen["runtime_manifest"]["entries"])
    manifest_ok(
        ROOT / frozen["final_supplement"]["path"],
        frozen["final_supplement"]["entries"],
        model["round1_candidate"],
    )

    stages = model["stage_table"]
    expected = {
        "code": ("target/wp3-core-docs/pages-source-path.txt", "/nix/store/n3z69s8mgc59ndwg5nm2cdc2lc368jm9-source"),
        "intermediate_round1": ("target/wp3-core-docs/pages-source-final-path.txt", "/nix/store/rq2sf8iv7aknnpcc8lxhpssivs3wlyvv-source"),
        "round1_final": ("target/wp3-core-docs/pages-source-final-round2-path.txt", "/nix/store/vvlqfaqqpgwl312cgwcc0mqnxyrg275f-source"),
    }
    for stage, (record, expected_path) in expected.items():
        page = stages[stage]["pages_source"]
        require(page["path_record"] == record, f"stage pointer {stage}")
        require(text(ROOT / record) == expected_path == page["path"], f"stage source {stage}")
        require(exit_value(ROOT / page["exit_record"]) == page["exit"] == 0, f"stage exit {stage}")
    require("rev=9b18f329e992ccc642c2cdaa62e934318f55fa76" in text(ROOT / stages["code"]["docs_check"]["argv"]), "code docs argv")
    require("rev=3f14fbe02c3b91a59243f1a16abfed58235d465a" in text(ROOT / stages["intermediate_round1"]["docs_check"]["argv"]), "intermediate docs argv")
    require("rev=8f874aa727f3a5e702d3a4b1c30f7798bb0d33d3" in text(ROOT / stages["round1_final"]["docs_check"]["argv"]), "round1 docs argv")
    for stage in ["code", "intermediate_round1", "round1_final"]:
        check = stages[stage]["docs_check"]
        require(exit_value(ROOT / check["exit_record"]) == check["exit"] == 0, f"docs exit {stage}")
    require(exit_value(ROOT / stages["failed_local_check"]["exit_record"]) == stages["failed_local_check"]["exit"] == 1, "failed local check")
    metadata = json.loads((ROOT / stages["round1_final"]["full_source"]["metadata_record"]).read_text())
    require(metadata["revision"] == stages["round1_final"]["revision"], "round1 full revision")
    require(metadata["path"] == stages["round1_final"]["full_source"]["path"], "round1 full source")

    manifest = tomllib.loads((ROOT / model["generated_evidence_inventory"]["manifest"]).read_text())
    artifacts = manifest["artifact"]
    derived = {
        "artifact": sorted({item["path"] for item in artifacts}),
        "source": sorted({item["source"] for item in artifacts}),
        "artifact_owner": sorted({item["artifact_owner"] for item in artifacts}),
        "source_owner": sorted({item["source_owner"] for item in artifacts}),
        "tool_owner": sorted({item["tool_owner"] for item in artifacts}),
    }
    inventory = model["generated_evidence_inventory"]
    require(len(artifacts) == inventory["artifacts"] == 8, "artifact count")
    require(derived == inventory["roles"], "role inventory")
    role_counts = {key: len(value) for key, value in derived.items()}
    distinct = sorted(set().union(*(set(value) for value in derived.values())))
    require(role_counts == inventory["role_counts"], "role counts")
    require(sum(role_counts.values()) == inventory["role_slots"] == 24, "role slots")
    require(len(distinct) == inventory["distinct_binding_files"] == 22, "distinct bindings")

    base = "aa841fde3195270cc00328a9356dfe523c0d7a85"
    code = "9b18f329e992ccc642c2cdaa62e934318f55fa76"
    round1 = "8f874aa727f3a5e702d3a4b1c30f7798bb0d33d3"
    code_paths = diff_paths(base, code)
    evidence_paths = diff_paths(code, round1)
    round1_paths = diff_paths(base, round1)
    expected_inventory = model["path_inventory"]["round1_from_base"]
    require(len(code_paths) == expected_inventory["code_paths"] == 9, "code path count")
    require(len(evidence_paths) == expected_inventory["evidence_paths"] == 7, "evidence path count")
    require(sorted(code_paths & evidence_paths) == expected_inventory["overlap"] == ["docs/developers.toml"], "path overlap")
    require(len(round1_paths) == expected_inventory["distinct_paths"] == 15, "round1 distinct paths")

    current = None
    if require_current:
        current_record = ROOT / stages["current_final"]["record"]
        current = json.loads(current_record.read_text())
        current_revision = subprocess.run(
            ["jj", "--ignore-working-copy", "log", "-r", "@", "--no-graph", "-T", "commit_id"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        require(current["revision"] == current_revision, "current revision")
        require(current["change"] == model["evidence_change"], "current change")
        full = current["full_source"]
        require(sha(ROOT / full["argv_record"]) == full["argv_sha256"], "current full argv hash")
        require(sha(ROOT / full["metadata_record"]) == full["metadata_sha256"], "current full metadata hash")
        require(exit_value(ROOT / full["exit_record"]) == full["exit"] == 0, "current full exit")
        require(text(ROOT / full["nar_hash_record"]) == full["nar_hash"], "current full nar")
        full_metadata = json.loads((ROOT / full["metadata_record"]).read_text())
        require(full_metadata["revision"] == current["revision"], "current full metadata revision")
        require(full_metadata["path"] == full["path"], "current full metadata path")
        pages = current["pages_source"]
        require(sha(ROOT / pages["argv_record"]) == pages["argv_sha256"], "current pages argv hash")
        require(sha(ROOT / pages["path_record"]) == pages["path_record_sha256"], "current pages path hash")
        require(text(ROOT / pages["nar_hash_record"]) == pages["nar_hash"], "current pages nar")
        require(text(ROOT / current["pages_source"]["path_record"]) == current["pages_source"]["path"], "current pages path")
        require(exit_value(ROOT / current["pages_source"]["exit_record"]) == current["pages_source"]["exit"] == 0, "current pages exit")
        require(current["revision"] in text(ROOT / current["pages_source"]["argv_record"]), "current pages argv revision")
        gate = current["metadata_free_gate"]
        require(sha(pathlib.Path(gate["probe_path"])) == gate["probe_sha256"], "current probe hash")
        require(sha(ROOT / gate["positive_argv"]) == gate["positive_argv_sha256"], "current positive argv hash")
        require(sha(ROOT / gate["negative_argv"]) == gate["negative_argv_sha256"], "current negative argv hash")
        require(current["pages_source"]["path"] in text(ROOT / gate["positive_argv"]), "current positive argv source")
        require(current["pages_source"]["path"] in text(ROOT / gate["negative_argv"]), "current negative argv source")
        require(exit_value(ROOT / gate["positive_exit_record"]) == gate["positive_exit"] == 0, "current positive gate")
        require(exit_value(ROOT / gate["negative_exit_record"]) == gate["negative_exit"] == 2, "current negative gate")
        docs = current["docs_check"]
        require(sha(ROOT / docs["argv_record"]) == docs["argv_sha256"], "current docs argv hash")
        require(sha(ROOT / docs["source_binding"]) == docs["source_binding_sha256"], "current docs source binding hash")
        require(sha(ROOT / docs["post_run_diff"]) == docs["post_run_diff_sha256"], "current docs diff hash")
        require(sha(ROOT / docs["exit_record"].replace(".exit", ".log")) == docs["log_sha256"], "current docs log hash")
        require(current["revision"] in text(ROOT / docs["argv_record"]), "current docs argv revision")
        require(full["path"] in text(ROOT / docs["argv_record"]), "current docs argv origin")
        require(exit_value(ROOT / docs["exit_record"]) == docs["exit"] == 0, "current docs exit")
        registry = current["registry"]
        require(sha(ROOT / registry["scope_record"]) == registry["scope_record_sha256"], "current scope record hash")
        require(sha(ROOT / registry["summary_record"]) == registry["summary_record_sha256"], "current registry record hash")
        producer = current["producer_invariance"]
        require(text(ROOT / producer["core_record"]) == producer["core_derivation"], "current core derivation")
        require(text(ROOT / producer["docs_record"]) == producer["docs_derivation"], "current docs derivation")
        require(current["revision"] in text(ROOT / producer["core_argv_record"]), "current core eval argv")
        require(current["revision"] in text(ROOT / producer["docs_argv_record"]), "current docs eval argv")
        require(exit_value(ROOT / producer["core_exit_record"]) == 0, "current core eval exit")
        require(exit_value(ROOT / producer["docs_exit_record"]) == 0, "current docs eval exit")
        correction_paths = diff_paths(round1, current_revision)
        require(len(correction_paths) == model["path_inventory"]["round2_correction"]["paths"], "round2 path count")
        require(len(round1_paths | correction_paths) == model["path_inventory"]["cumulative_from_base"]["distinct_paths"], "cumulative path count")

    return {
        "stages": {key: value.get("classification", "bound") for key, value in stages.items()},
        "binding_inventory": {
            "artifacts": len(artifacts),
            "role_counts": role_counts,
            "role_slots": sum(role_counts.values()),
            "distinct_binding_files": len(distinct),
            "distinct_files": [{"path": path, "sha256": sha(ROOT / path)} for path in distinct],
        },
        "path_inventory": {
            "code": len(code_paths),
            "evidence": len(evidence_paths),
            "overlap": sorted(code_paths & evidence_paths),
            "round1_distinct": len(round1_paths),
            "current": current,
        },
    }


model = json.loads(CORRECTION.read_text())
result = validate(model, require_current="--pre-freeze" not in sys.argv)

negative_results = {}
bad_count = copy.deepcopy(model)
bad_count["generated_evidence_inventory"]["distinct_binding_files"] = 25
try:
    validate(bad_count, require_current=False)
    negative_results["wrong_25_count"] = "UNEXPECTED_ACCEPT"
except ValueError:
    negative_results["wrong_25_count"] = "REJECT"

swapped = copy.deepcopy(model)
swapped["stage_table"]["intermediate_round1"]["pages_source"], swapped["stage_table"]["round1_final"]["pages_source"] = swapped["stage_table"]["round1_final"]["pages_source"], swapped["stage_table"]["intermediate_round1"]["pages_source"]
try:
    validate(swapped, require_current=False)
    negative_results["swapped_intermediate_final"] = "UNEXPECTED_ACCEPT"
except ValueError:
    negative_results["swapped_intermediate_final"] = "REJECT"

bad_exit = copy.deepcopy(model)
bad_exit["stage_table"]["code"]["pages_source"]["exit"] = 2
try:
    validate(bad_exit, require_current=False)
    negative_results["wrong_code_exit"] = "UNEXPECTED_ACCEPT"
except ValueError:
    negative_results["wrong_code_exit"] = "REJECT"

require(set(negative_results.values()) == {"REJECT"}, "negative controls")
result["negative_controls"] = negative_results
print(json.dumps(result, indent=2, sort_keys=True))
