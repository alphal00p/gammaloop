#!/usr/bin/env python3
"""Collect GammaLoop generation and integration diagnostics for one manual graph."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

ANSI = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")
SAMPLE_CORE = re.compile(
    r"([0-9]+(?:\.[0-9]+)?\s*(?:ns|µs|us|ms|s|min))\s*/sample/core"
)
GRAPH_IDS = [
    "GL076",
    "GL084",
    "GL086",
    "GL090",
    "GL098",
    "GL120",
    "GL132",
    "GL134",
    "GL136",
    "GL144",
    "GL148",
    "GL160",
    "GL164",
    "GL168",
    "GL170",
    "GL190",
    "GL194",
    "GL280",
    "GL286",
    "GL290",
    "GL292",
    "GL293",
    "GL297",
    "GL304",
    "GL310",
    "GL330",
    "GL334",
    "GL358",
    "GL362",
    "GL364",
    "GL368",
    "GL382",
    "GL438",
    "GL443",
    "GL448",
    "GL460",
    "GL488",
    "GL490",
    "GL494",
    "GL495",
    "GL500",
    "GL516",
    "GL524",
    "GL534",
    "GL554",
    "GL558",
    "GL562",
    "GL564",
    "GL568",
    "GL572",
    "GL588",
    "GL590",
    "GL592",
    "GL598",
    "GL612",
    "GL616",
    "GL636",
    "GL638",
    "GL644",
    "GL648",
    "GL652",
    "GL654",
    "GL662",
    "GL664",
    "GL672",
    "GL673",
    "GL674",
    "GL677",
    "GL679",
    "GL701",
    "GL703",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph-id")
    parser.add_argument("--state", type=Path)
    parser.add_argument("--integration-log", type=Path)
    parser.add_argument("--threshold-edited", choices=("yes", "no"))
    parser.add_argument("--all-limits", choices=("yes", "no", "inconclusive"))
    parser.add_argument("--uv-mode", choices=("full", "local-only", "blocked"))
    parser.add_argument("--observation", default="")
    parser.add_argument("--output-record", type=Path)
    parser.add_argument("--output-md", type=Path)
    parser.add_argument("--summary-states-root", type=Path)
    parser.add_argument("--summary-output", type=Path)
    args = parser.parse_args()
    summary_mode = (
        args.summary_states_root is not None or args.summary_output is not None
    )
    if summary_mode:
        if args.summary_states_root is None or args.summary_output is None:
            parser.error(
                "summary mode requires both --summary-states-root and --summary-output"
            )
        if args.graph_id is not None:
            parser.error("summary mode cannot be combined with per-graph collection")
        return args
    required = {
        "--graph-id": args.graph_id,
        "--state": args.state,
        "--threshold-edited": args.threshold_edited,
        "--all-limits": args.all_limits,
        "--uv-mode": args.uv_mode,
    }
    missing = [name for name, value in required.items() if value is None]
    if missing:
        parser.error("per-graph mode requires " + ", ".join(missing))
    return args


def duration_seconds(raw: dict) -> float:
    return float(raw.get("secs", 0)) + float(raw.get("nanos", 0)) * 1.0e-9


def generation_data(state: Path) -> dict:
    candidates = sorted(state.rglob("generation_summary.json"))
    if not candidates:
        return {"status": "missing", "time_seconds": None, "peak_ram_bytes": None}
    reports = []
    peak_ram = 0
    for path in candidates:
        raw = json.loads(path.read_text())
        peak_ram = max(peak_ram, int(raw.get("peak_ram_bytes", 0)))
        reports.extend(raw.get("reports", []))
    total = sum(
        duration_seconds(report.get("stats", {}).get("total_time", {}))
        for report in reports
    )
    return {
        "status": "available",
        "time_seconds": total,
        "peak_ram_bytes": peak_ram,
        "reports": len(reports),
        "sources": [str(path) for path in candidates],
    }


def complex_value(raw: dict | None) -> dict | None:
    if raw is None:
        return None
    return {"re": raw.get("re"), "im": raw.get("im")}


def integration_data(state: Path, log_path: Path | None) -> dict:
    candidates = sorted(state.rglob("integration_result.json"))
    if not candidates:
        return {"status": "missing"}
    path = candidates[-1]
    raw = json.loads(path.read_text())
    slots = raw.get("slots", [])
    if len(slots) != 1:
        return {
            "status": "invalid",
            "source": str(path),
            "error": f"expected one integration slot, found {len(slots)}",
        }
    slot = slots[0]
    integral = slot.get("integral", {})
    absolute = slot.get("absolute", {})
    absolute_integral = absolute.get("integral", {})
    statistics = slot.get("integration_statistics", {})
    absolute_relative_errors = {}
    active_relative_errors = []
    for row in absolute.get("table_results", []):
        component = row.get("component")
        relative = row.get("relative_error_percent")
        value = row.get("value")
        absolute_relative_errors[component] = relative
        if isinstance(value, (int, float)) and value != 0.0:
            if not isinstance(relative, (int, float)) or not math.isfinite(relative):
                active_relative_errors.append(math.inf)
            else:
                active_relative_errors.append(float(relative))
    nan_percentage = statistics.get("nan_percentage", 0.0)
    finite_values = [
        value
        for container in (
            integral.get("result", {}),
            integral.get("error", {}),
            absolute_integral.get("result", {}),
            absolute_integral.get("error", {}),
        )
        for value in container.values()
        if isinstance(value, (int, float))
    ]
    converged = (
        len(finite_values) == 8
        and bool(integral.get("neval", 0))
        and all(math.isfinite(value) for value in finite_values)
    )
    converged = converged and nan_percentage == 0.0 and bool(active_relative_errors)
    converged = converged and all(
        relative <= 75.0 for relative in active_relative_errors
    )

    sample_core = None
    if log_path is not None and log_path.exists():
        matches = SAMPLE_CORE.findall(
            ANSI.sub("", log_path.read_text(errors="replace"))
        )
        if matches:
            sample_core = matches[-1].replace(" ", "")
    return {
        "status": "available",
        "source": str(path),
        "neval": integral.get("neval"),
        "result": complex_value(integral.get("result")),
        "error": complex_value(integral.get("error")),
        "absolute_result": complex_value(absolute_integral.get("result")),
        "absolute_error": complex_value(absolute_integral.get("error")),
        "absolute_relative_error_percent": absolute_relative_errors,
        "converged": converged,
        "runtime_per_sample_per_core": sample_core,
        "statistics": statistics,
        "max_weight_info": slot.get("max_weight_info", []),
        "absolute_max_weight_info": absolute.get("max_weight_info", []),
    }


def format_complex(value: dict | None) -> str:
    if value is None:
        return "—"
    return f"({value.get('re')!s}, {value.get('im')!s})"


def markdown_cell(value: object) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def summary_markdown(states_root: Path) -> str:
    header = (
        "| graph_id | threshold struct. edited? | converged? | all limits OK? | "
        "central value | error | absolute central value | absolute error | n_points | "
        "generation time | runtime per sample per core | short observation |\n"
        "|---|---|---|---|---|---|---|---|---:|---|---|---|\n"
    )
    rows = []
    completed = 0
    for graph_id in GRAPH_IDS:
        record_path = states_root / f"state_{graph_id}" / "result_record.json"
        if not record_path.exists():
            values = [
                graph_id,
                "pending",
                "pending",
                "pending",
                "—",
                "—",
                "—",
                "—",
                "—",
                "—",
                "—",
                "Awaiting audit.",
            ]
        else:
            record = json.loads(record_path.read_text())
            generation = record.get("generation", {})
            integration = record.get("integration", {})
            generation_time = generation.get("time_seconds")
            completed += 1
            values = [
                graph_id,
                record.get("threshold_edited", "inconclusive"),
                "yes" if integration.get("converged") else "no",
                record.get("all_limits", "inconclusive"),
                format_complex(integration.get("result")),
                format_complex(integration.get("error")),
                format_complex(integration.get("absolute_result")),
                format_complex(integration.get("absolute_error")),
                integration.get("neval") or "—",
                "—" if generation_time is None else f"{generation_time:.6g} s",
                integration.get("runtime_per_sample_per_core") or "—",
                record.get("observation") or "Audit recorded.",
            ]
        rows.append("| " + " | ".join(markdown_cell(value) for value in values) + " |")
    return (
        "# Manual NNLO ttH IR-safety audit\n\n"
        f"Progress: **{completed}/{len(GRAPH_IDS)}** graph audits recorded.\n\n"
        + header
        + "\n".join(rows)
        + "\n"
    )


def markdown(record: dict) -> str:
    generation = record["generation"]
    integration = record["integration"]
    generation_time = generation.get("time_seconds")
    generation_display = "—" if generation_time is None else f"{generation_time:.6g} s"
    return f"""# {record["graph_id"]}

## Status

- Threshold structure edited: {record["threshold_edited"]}
- All limits OK: {record["all_limits"]}
- UV mode: {record["uv_mode"]}
- Converged: {integration.get("converged", False)}
- Observation: {record["observation"] or "Pending detailed graph audit."}

## Generation

- Internal generation time: {generation_display}
- Peak generation RAM: {generation.get("peak_ram_bytes")}
- Status: {generation.get("status")}

## IR limits

Record every limit's edges, midpoint, axis, rank, target-distance fits, signed-branch p envelope, and classification here.

## Threshold directives

Record the original and final compact specification, subspace rationale, and threshold-only checks here.

## Integration

- Signed central value `(re, im)`: {format_complex(integration.get("result"))}
- Signed error `(re, im)`: {format_complex(integration.get("error"))}
- Absolute central value `(|re|, |im|)`: {format_complex(integration.get("absolute_result"))}
- Absolute error `(|re|, |im|)`: {format_complex(integration.get("absolute_error"))}
- Samples: {integration.get("neval", "—")}
- Runtime per sample per core: {integration.get("runtime_per_sample_per_core") or "—"}
- Absolute relative errors: {integration.get("absolute_relative_error_percent", {})}
- Stability statistics: {integration.get("statistics", {})}

## Max-weight inspection

Signed:

```json
{json.dumps(integration.get("max_weight_info", []), indent=2)}
```

Componentwise absolute:

```json
{json.dumps(integration.get("absolute_max_weight_info", []), indent=2)}
```
"""


def main() -> int:
    args = parse_args()
    if args.summary_output is not None:
        args.summary_output.parent.mkdir(parents=True, exist_ok=True)
        args.summary_output.write_text(summary_markdown(args.summary_states_root))
        return 0
    assert args.graph_id is not None
    assert args.state is not None
    assert args.threshold_edited is not None
    assert args.all_limits is not None
    assert args.uv_mode is not None
    record = {
        "schema_version": 1,
        "graph_id": args.graph_id,
        "threshold_edited": args.threshold_edited,
        "all_limits": args.all_limits,
        "uv_mode": args.uv_mode,
        "observation": args.observation,
        "generation": generation_data(args.state),
        "integration": integration_data(args.state, args.integration_log),
    }
    if args.output_record is not None:
        args.output_record.parent.mkdir(parents=True, exist_ok=True)
        args.output_record.write_text(
            json.dumps(record, indent=2, sort_keys=True) + "\n"
        )
    if args.output_md is not None:
        args.output_md.parent.mkdir(parents=True, exist_ok=True)
        args.output_md.write_text(markdown(record))
    if args.output_record is None and args.output_md is None:
        print(json.dumps(record, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
