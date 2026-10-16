"""Live ImPrEd reruns from certified EC seeds with native constraint projection.

Each request starts from a saved, certified seed. Temporary inputs isolate all
native projection inputs; a bounded cache avoids repeating completed layouts.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
import tempfile
import time
from collections import OrderedDict
from concurrent.futures import CancelledError
from pathlib import Path

from package_stages import OUT, STAGES, load_stage, seed_metadata
from run_constrained import NativeProjector
from solver import solve

# The metadata is also the validation contract consumed by the browser.
_PARAMETER_ROWS = [
    (
        "impred",
        "target",
        "Spacing",
        0.5,
        6.0,
        0.1,
        2.4,
        "Characteristic length for ImPrEd forces and movement limits; the certified seed stays fixed.",
    ),
    (
        "impred",
        "repulsion",
        "Point repulsion",
        0.0,
        4.0,
        0.1,
        2.5,
        "Multiplier on the existing point–point repulsion.",
    ),
    (
        "impred",
        "attraction",
        "Segment attraction",
        0.0,
        4.0,
        0.1,
        2.5,
        "Multiplier on each segment's attraction. Added route points let a leg lengthen and bend; unequal segment lengths pull the points toward more even spacing.",
    ),
    (
        "impred",
        "parallel_attraction_balance",
        "Parallel-edge attraction balancing",
        0.0,
        2.0,
        0.05,
        1.0,
        "Reduce each parallel internal edge's attraction by multiplicity raised to this exponent. Zero preserves unbalanced attraction; one compensates the edge count. Bridges, external legs and self-loops keep their original attraction.",
    ),
    (
        "impred",
        "pull",
        "External pull",
        0.0,
        32.0,
        0.05,
        0.45,
        "Constant left/right pull for mixed flow; radial pull when all legs have the same flow.",
    ),
    (
        "impred",
        "pull_balance",
        "Pull balancing",
        0.0,
        4.0,
        0.1,
        1.0,
        "Exponent on the topology multiplier. Zero gives uniform pull; values above one reduce multipliers below one and amplify multipliers above one.",
    ),
    (
        "impred",
        "pull_attachment",
        "Distributed attachment pull",
        0.0,
        32.0,
        0.5,
        4.0,
        "Extra pull for endpoints sharing an X coordinate whose attachment nodes are distributed through the network. Zero keeps only bottleneck balancing; terminal attachments and radial flow are unaffected.",
    ),
    (
        "impred",
        "external_max_points",
        "External route-point limit",
        0,
        3,
        1,
        2,
        "Maximum intermediate points per uncrossed external leg. Zero disables flexibility; subdivision and contraction thresholds control when points change.",
    ),
    (
        "impred",
        "split_length_ratio",
        "Subdivision threshold / spacing",
        0.25,
        6.0,
        "any",
        1.5,
        "Split a flexible segment above this multiple of spacing, subject to the point cap. Must exceed the contraction threshold.",
    ),
    (
        "impred",
        "contract_chord_ratio",
        "Contraction threshold / spacing",
        0.1,
        5.0,
        "any",
        1.25,
        "Remove a free route point when its neighbours are closer than this multiple of spacing and topology permits. Must be below the subdivision threshold.",
    ),
    (
        "impred",
        "edge_clearance",
        "Edge clearance",
        0.05,
        1.0,
        0.05,
        0.4,
        "Point–edge repulsion range as a fraction of spacing.",
    ),
    (
        "impred",
        "node_edge_strength",
        "Point–edge repulsion",
        0.0,
        12.0,
        0.25,
        4.0,
        "Repulsion where a point projects inside an edge segment, with balanced reactions on its endpoints.",
    ),
    (
        "impred",
        "steps",
        "Iterations",
        100,
        3000,
        100,
        1500,
        "Iterations over the complete exponent and cooling schedules.",
    ),
]
PARAMETERS = [
    dict(
        zip(
            ("stage", "key", "label", "min", "max", "step", "default", "description"),
            (stage, f"{stage}.{key}", *values),
        )
    )
    for stage, key, *values in _PARAMETER_ROWS
]
DEFAULTS = {p["key"]: p["default"] for p in PARAMETERS}
_STAGE_INDEX = {spec[0]: index for index, spec in enumerate(STAGES)}


class PlaygroundRunner:
    """Serialized runner with immutable source snapshots and a bounded ImPrEd cache."""

    def __init__(self, root=OUT):
        self.root = Path(root)
        self._sources = {}
        cases = {}
        for manifest in (
            "cases.json",
            "cases-multileg-extra.json",
            "cases-nonplanar-extra.json",
            "cases-xs-extra.json",
        ):
            path = self.root / manifest
            if path.exists():
                for case in json.loads(path.read_text()):
                    cases.setdefault(case["id"], case)
        for case in cases.values():
            folder = self.root / case["id"]
            seed_path = folder / "constrained-initial.json"
            if not seed_path.exists():
                continue
            self._sources[case["id"]] = {
                "case": {key: case[key] for key in ("id", "letter", "title")},
                "seed": json.loads(seed_path.read_text()),
                "graph": seed_path.with_suffix(".projected.graph.json").read_bytes(),
                "fingerprint": hashlib.sha256(seed_path.read_bytes()).hexdigest(),
            }
        self._cache = OrderedDict()
        self._config = None

    def config(self):
        if self._config is None:
            cases = []
            with tempfile.TemporaryDirectory(prefix="impred-playground-config-") as tmp:
                for identifier, source in self._sources.items():
                    item = self._case(source)
                    item["stages"] = [
                        load_stage(
                            self.root / identifier, spec, Path(tmp) / "audit.json"
                        )
                        for spec in STAGES
                    ]
                    cases.append(item)
            self._config = {
                "cases": cases,
                "parameters": PARAMETERS,
                "defaults": DEFAULTS,
            }
        return copy.deepcopy(self._config)

    def validate(self, diagram, stage, parameters):
        if not isinstance(diagram, str) or diagram not in self._sources:
            raise ValueError("Unknown diagram")
        if not isinstance(stage, str) or stage not in _STAGE_INDEX:
            raise ValueError("Unknown layout stage")
        if not isinstance(parameters, dict):
            raise ValueError("Parameters must be an object")  # noqa: TRY004 - shared HTTP validation contract
        unknown = set(parameters) - DEFAULTS.keys()
        if unknown:
            raise ValueError(
                f"Unknown parameter: {', '.join(sorted(map(str, unknown)))}"
            )
        normalized = dict(DEFAULTS)
        for spec in PARAMETERS:
            key = spec["key"]
            value = parameters.get(key, spec["default"])
            try:
                finite = isinstance(value, (int, float)) and math.isfinite(value)
            except OverflowError:
                finite = False
            if isinstance(value, bool) or not finite:
                raise ValueError(f"{key} must be a finite number")
            if not spec["min"] <= value <= spec["max"]:
                raise ValueError(
                    f"{key} must be between {spec['min']} and {spec['max']}"
                )
            if isinstance(spec["default"], int):
                if value != int(value):
                    raise ValueError(f"{key} must be an integer")
                value = int(value)
            else:
                value = float(value)
            normalized[key] = value
        if (
            normalized["impred.contract_chord_ratio"]
            >= normalized["impred.split_length_ratio"]
        ):
            raise ValueError(
                "Contraction threshold must be below subdivision threshold"
            )
        return normalized

    @staticmethod
    def _case(source):
        return {
            **source["case"],
            **seed_metadata(source["seed"]),
            "seed_fingerprint": source["fingerprint"],
            "constraints_summary": (
                "Incoming pull left · Outgoing pull right"
                if source["seed"]["incoming_ids"] and source["seed"]["outgoing_ids"]
                else "External legs pull radially outward"
            ),
            "constraints_note": "The graph’s coordinate and alignment constraints remain active. External legs may gain route points.",
        }

    @staticmethod
    def _check_cancelled(cancelled):
        if cancelled is not None and cancelled():
            raise CancelledError()

    def _key(self, diagram, parameters):
        return (
            diagram,
            self._sources[diagram]["fingerprint"],
            tuple(parameters.items()),
        )

    def _remember(self, key, value):
        self._cache[key] = value
        self._cache.move_to_end(key)
        if len(self._cache) > 32:
            self._cache.popitem(last=False)

    def run(self, diagram, stage, parameters, cancelled=None):
        parameters = self.validate(diagram, stage, parameters)
        started = time.monotonic()
        self._check_cancelled(cancelled)
        source = self._sources[diagram]
        item = self._case(source)
        with tempfile.TemporaryDirectory(prefix="impred-playground-run-") as tmp:
            folder = Path(tmp)
            seed_path = folder / STAGES[0][2]
            seed_path.write_text(json.dumps(source["seed"]))
            graph = json.loads(source["graph"])
            graph["layout_config"]["external-pull-attachment"] = parameters[
                "impred.pull_attachment"
            ]
            seed_path.with_suffix(".projected.graph.json").write_text(json.dumps(graph))
            stages = [load_stage(folder, STAGES[0], folder / "audit.json")]
            for spec in STAGES[1:]:
                current = spec[0]
                self._check_cancelled(cancelled)
                if _STAGE_INDEX[current] > _STAGE_INDEX[stage]:
                    stages.append(
                        {
                            "id": current,
                            "label": spec[1],
                            "source": spec[2],
                            "available": False,
                            "reason": "Choose this stage to run it with the current parameters.",
                        }
                    )
                    continue
                key = self._key(diagram, parameters)
                cached = self._cache.get(key)
                output = folder / spec[2]
                if cached is not None:
                    payload = cached
                    self._cache.move_to_end(key)
                else:
                    session = NativeProjector(seed_path)
                    try:
                        knobs = {
                            key.split(".", 1)[1]: value
                            for key, value in parameters.items()
                            if key != "impred.pull_attachment"
                        }
                        state = solve(
                            source["seed"],
                            flexible=False,
                            project=session,
                            pull_scales=session.ready["external_pull_scales"],
                            cancelled=cancelled,
                            **knobs,
                        )
                    finally:
                        session.close()
                    output.write_text(json.dumps(state))
                    payload = load_stage(folder, spec, folder / "audit.json")
                    self._check_cancelled(cancelled)
                    self._remember(key, payload)
                stages.append(copy.deepcopy(payload))
            item["stages"] = stages
        self._check_cancelled(cancelled)
        return {
            "case": item,
            "parameters": parameters,
            "elapsed": time.monotonic() - started,
        }
