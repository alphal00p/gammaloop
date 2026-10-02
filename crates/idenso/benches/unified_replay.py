"""Replay existing benchmark children against an independently checked saved cohort.

This does not implement tensor operations or an alternative timing path. The
existing source-bound child retains admission, clocks, warmups, phase captures,
completion assertions and result serialization. Each result must equal the
previously independently checked expression exactly, including its input hash.
"""

import argparse
import hashlib
import json
import os
import subprocess
import time
from pathlib import Path
from statistics import median

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--reference", type=Path, required=True)
parser.add_argument("--candidate-frontend", type=Path, required=True)
parser.add_argument("--candidate-python", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
reference = args.reference.resolve()
old = json.loads(reference.read_text())
assert old["completed"]
frontend = json.loads(args.candidate_frontend.read_text())
core = frontend["core_sha256"]
for name, digest in frontend["sources"].items():
    assert hashlib.sha256(Path(name).read_bytes()).hexdigest() == digest
out = args.output.resolve()
out.mkdir(exist_ok=False)
report_path = out.with_suffix(".json")
report = {
    k: old[k]
    for k in [
        "cases",
        "rounds",
        "timing",
        "diagnostic_scope",
        "scheduling",
        "host_scope",
        "ladder_routes",
        "contraction_order_scope",
    ]
}
report.update(
    schema=1,
    protocol="Exact replay of measured child commands; same clocks/warmups and rotating labels; new checked source/core and fresh baseline samples.",
    reference={
        "path": str(reference),
        "sha256": hashlib.sha256(reference.read_bytes()).hexdigest(),
        "core_identities": old["core_identities"],
        "checker_scope": "The saved cohort independently checked baseline+host11 against FORM and exact HEP components. New input hashes and raw serialized expressions must equal those checked outputs, so equality transfers that certificate. FORM clocks retain their original timestamps.",
    },
    source_manifest=frontend,
    records=[],
    summary={},
    completed=False,
)
start = time.monotonic()
for prior in old["records"]:
    command = prior["command"].copy()
    label = prior["label"]
    name = Path(prior["output"]).name
    output = out / name
    command[command.index("--output") + 1] = str(output)
    if label == "candidate":
        # Preserve the venv interpreter symlink and its pyvenv.cfg lookup.
        command[0] = str(args.candidate_python.absolute())
        command[1] = frontend["entrypoint"]
        command[command.index("--expected-core") + 1] = core
    call_start = time.monotonic()
    with (
        output.with_suffix(".stdout").open("w") as stdout,
        output.with_suffix(".stderr").open("w") as stderr,
    ):
        result = subprocess.run(
            command,
            stdout=stdout,
            stderr=stderr,
            timeout=900,
            check=False,
            env=os.environ | {"PYTHONNOUSERSITE": "1"},
        )
    row = {
        "label": label,
        "case": prior["case"],
        "route": prior["route"],
        "command": command,
        "exit_code": result.returncode,
        "process_wall_seconds": time.monotonic() - call_start,
        "output": str(output),
        "reference_output_sha256": hashlib.sha256(
            Path(prior["record"]["expression"]).read_bytes()
        ).hexdigest(),
    }
    report["records"].append(row)
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    assert result.returncode == 0, (command, output.with_suffix(".stderr"))
    value = json.loads(output.read_text())
    row["record"] = value
    row["same_input_sha256"] = (
        value["input"]["sha256"] == prior["record"]["input"]["sha256"]
    )
    row["exact_serialized_output_equal"] = (
        Path(value["expression"]).read_bytes()
        == Path(prior["record"]["expression"]).read_bytes()
    )
    assert row["same_input_sha256"] and row["exact_serialized_output_equal"], row
    for case in old["cases"]:
        summary = {}
        for variant in ["baseline", "candidate"]:
            records = [
                r["record"]
                for r in report["records"]
                if r["case"] == case
                and r["label"] == variant
                and r["route"] != "scalar"
            ]
            if not records:
                continue
            item = {
                key: median(r[key] / r["calls"] for r in records)
                for key in ["wall_ns", "process_ns", "thread_ns"]
            }
            if all("complete_algebra" in r for r in records):
                item["complete_algebra"] = {
                    key: median(
                        r["complete_algebra"][key] / r["complete_algebra"]["calls"]
                        for r in records
                    )
                    for key in ["wall_ns", "process_ns", "thread_ns"]
                }
            summary[variant] = item
        if "FORM" in old["summary"][case]:
            summary["FORM"] = old["summary"][case]["FORM"] | {
                "retained_from_checked_cohort": str(reference)
            }
        report["summary"][case] = summary
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(label, prior["case"], prior["route"], "exact", flush=True)
for name, digest in frontend["sources"].items():
    assert hashlib.sha256(Path(name).read_bytes()).hexdigest() == digest
report["completed"] = True
report["process_wall_seconds"] = time.monotonic() - start
report_path.write_text(json.dumps(report, indent=2) + "\n")
print("completed", report["process_wall_seconds"], flush=True)
