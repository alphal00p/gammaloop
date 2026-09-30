"""Attribute perf-script samples to disjoint reduction categories, leaf first."""

import collections
import json
import pathlib
import re
import sys

HEADER = re.compile(r"^\S+\s+\d+/\d+\s+\d+\.\d+:\s+(\d+)\s+cpu-clock:u:")
FRAME = re.compile(r"^\s+[0-9a-f]+\s+(.*)\s+\([^)]*\)\s*$")
MARKERS = {
    "observations": ("::DomainObservations", "::RegionObservation"),
    "trace_kernel": ("idenso::dirac::simplify::trace_kernel::",),
    "dot_notation": ("schoonschip::normalize_dots::", "::to_dots"),
}


def categories(frame):
    return [
        category
        for category, markers in MARKERS.items()
        if any(marker in frame for marker in markers)
    ]


def samples(path):
    period = None
    frames = []
    with path.open() as source:
        for line in source:
            match = HEADER.match(line)
            if match:
                if period is not None:
                    yield period, frames
                period = int(match[1])
                frames = []
            else:
                match = FRAME.match(line)
                if match:
                    frames.append(match[1])
    if period is not None:
        yield period, frames


def analyze(path):
    counts = collections.Counter()
    weights = collections.Counter()
    inclusive = collections.Counter()
    inclusiven = collections.Counter()
    overlaps = collections.Counter()
    periods = collections.Counter()
    depths = collections.Counter()
    reached = collections.Counter()
    reached_weight = collections.Counter()
    self_symbols = collections.defaultdict(collections.Counter)
    other_native = collections.Counter()
    other_native_n = collections.Counter()
    for period, frames in samples(path):
        labels = [
            (depth, category)
            for depth, frame in enumerate(frames)
            for category in categories(frame)
        ]
        bucket = min(labels)[1] if labels else "other"
        present = {category for _, category in labels}
        for category in present:
            inclusive[category] += period
            inclusiven[category] += 1
        overlaps[",".join(sorted(present))] += 1
        counts[bucket] += 1
        weights[bucket] += period
        periods[period] += 1
        depths[len(frames)] += 1
        if frames:
            self_symbols[bucket][frames[0]] += period
        for marker, label in [
            (
                "spynso3::expression::TensorExpression>::simplify_algebra",
                "python_binding",
            ),
            (">::plan_reduction", "planner"),
            ("observation::SettledRegions", "settled_regions_wrapper"),
        ]:
            if any(marker in frame for frame in frames):
                reached[label] += 1
                reached_weight[label] += period
        if bucket == "other":
            native = next(
                (
                    frame
                    for frame in frames
                    if frame.startswith(("idenso::", "<idenso::"))
                ),
                "(no principal idenso frame)",
            )
            other_native[native] += period
            other_native_n[native] += 1
    total = sum(weights.values())
    n = sum(counts.values())
    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "method": "Leaf-to-root nearest marker; same-frame tie alphabetical; observations only DomainObservations/RegionObservation; trace kernel exact namespace; dots normalize_dots namespace or ::to_dots; residual other.",
        "samples": n,
        "total_period": total,
        "buckets": {
            key: {
                "samples": counts[key],
                "period": weights[key],
                "percent": 100 * weights[key] / total,
            }
            for key in sorted(counts)
        },
        "inclusive": {
            key: {
                "samples": inclusiven[key],
                "period": inclusive[key],
                "percent": 100 * inclusive[key] / total,
            }
            for key in sorted(inclusive)
        },
        "reached": {
            key: {
                "samples": reached[key],
                "period": reached_weight[key],
                "percent": 100 * reached_weight[key] / total,
            }
            for key in sorted(reached)
        },
        "overlaps": dict(overlaps),
        "period_distribution": dict(periods),
        "depth": {
            "min": min(depths),
            "max": max(depths),
            "mean": sum(key * value for key, value in depths.items()) / n,
        },
        "other_native": [
            {"symbol": key, "samples": other_native_n[key], "period": value}
            for key, value in other_native.most_common(15)
        ],
        "self_symbols": {
            key: dict(value.most_common(12)) for key, value in self_symbols.items()
        },
    }


for name in sys.argv[1:]:
    path = pathlib.Path(name)
    result = analyze(path)
    path.with_name(path.stem + ".exclusive-independent.json").write_text(
        json.dumps(result, indent=2)
    )
    print(
        json.dumps(
            {
                key: value
                for key, value in result.items()
                if key not in ("other_native", "self_symbols")
            },
            indent=2,
        )
    )
    print("OTHER PRINCIPAL FRAMES", json.dumps(result["other_native"][:6], indent=2))
