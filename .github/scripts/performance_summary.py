"""Print the Actions job summary from the same capture sent to gammabench."""

import json
import sys
from pathlib import Path


def summary(manifest):
    if not manifest.is_file():
        return "### Performance\n\nThe job stopped before it recorded a capture."

    run = json.loads(manifest.read_text())
    ci = run["ci"]
    candidate = ci["candidate"]
    attempt = f"{ci['run_id']}-{ci['run_attempt']}"
    page = None
    if ci["event"] == "pull_request" and ci.get("pr"):
        page = f"previews/pr-{ci['pr']}/{candidate}/{attempt}"
    elif ci["event"] == "push":
        page = f"docs/results/ci-runs/gammaloop-{candidate}-{attempt}"

    counts = [0, 0, 0]
    findings = []
    for case in run["cases"]:
        gate = manifest.parent / case / "gate.json"
        if not gate.is_file():
            counts[2] += 1
            findings.append((2, case, "No complete measurement"))
            continue
        report = json.loads(gate.read_text())["report"]
        status = report["status"]
        counts[status] += 1
        if status == 0:
            continue
        details = []
        for metric, result in report["metrics"].items():
            if result["status"] == 0:
                continue
            phase = "confirmation" if "confirmation" in result else "discovery"
            change = result[phase]["median_relative_delta"]
            outcome = "regression" if result["status"] == 1 else "inconclusive"
            details.append(f"`{metric}` {change:+.2%} in {phase} ({outcome})")
        if report.get("invalid_input") or report.get("identity_mismatch"):
            details.append("Incomplete or invalid evidence")
        else:
            mismatches = report["correctness_mismatches"]
            details.append(
                f"{mismatches} correctness mismatch{'es' if mismatches != 1 else ''}"
            )
        findings.append((status, case, "; ".join(details)))

    lines = ["### Performance", ""]
    if run["cases"]:
        lines.append(
            f"**{counts[1]} failed · {counts[2]} inconclusive · {counts[0]} passed**"
        )
    else:
        lines.append("**Inconclusive** · The job did not record any measured cases.")
    if findings:
        lines.append("")
        for _, case, detail in sorted(findings):
            lines.append(f"- **{case}** — {detail}")
    lines += [
        "",
        f"`{candidate[:8]}` against `{ci['baseline'][:8]}` · attempt {ci['run_attempt']}",
    ]
    if page:
        url = f"https://alphal00p.github.io/gammabench/{page}/index.html"
        lines += [
            "",
            f"[Performance report for this run]({url}) (available after publication).",
        ]
    return "\n".join(lines)


if __name__ == "__main__":
    print(summary(Path(sys.argv[1])))
