#!/usr/bin/env python3
"""Single-CPU author reference; no grid generation or native integration."""

import argparse
import hashlib
import json
import math
import os
import pathlib
import signal
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent
SOURCE_SHA = "207a2cb3d39933a57a5c1a8592d664c7dc9d5750eb8cd10801e98ac9560446fb"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--energy", type=float, required=True)
parser.add_argument("--momentum", type=float, required=True)
parser.add_argument("--tolerance", type=float, default=1e-5)
parser.add_argument("--basis", choices=["both", "vector", "density"], default="both")
parser.add_argument(
    "--timeout",
    type=float,
    default=1800,
    help="total seconds across all upstream calls",
)
parser.add_argument(
    "--rest-fraction",
    type=float,
    default=0.05,
    help="largest rest-control k/E; others are half and quarter",
)
parser.add_argument(
    "--cache-directory",
    type=pathlib.Path,
    default=pathlib.Path(
        os.environ.get(
            "QCD_FIG15_REFERENCE_CACHE",
            os.path.join(os.environ.get("TMPDIR", "/tmp"), "qcd-fig15-reference-cache"),
        )
    ),
)
parser.add_argument(
    "--binary",
    type=pathlib.Path,
    default=pathlib.Path(
        os.environ.get("QCD_FIG15_REFERENCE_BINARY", str(ROOT / "reference"))
    ),
)
parser.add_argument("--cpu", type=int)
parser.add_argument("--output", type=pathlib.Path)
args = parser.parse_args()
if args.cpu is not None:
    os.sched_setaffinity(0, {args.cpu})
for name in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS",
):
    os.environ[name] = "1"
if not all(
    math.isfinite(x)
    for x in (
        args.energy,
        args.momentum,
        args.tolerance,
        args.timeout,
        args.rest_fraction,
    )
):
    parser.error("all numeric arguments must be finite")
if args.energy <= 0 or args.momentum < 0 or args.tolerance <= 0 or args.timeout <= 0:
    parser.error("require energy>0, momentum>=0, tolerance>0, timeout>0")
if not 0 < args.rest_fraction < 0.25:
    parser.error("rest-fraction must be between0 and.25")
if args.energy == args.momentum:
    parser.error(
        "strict fixed-order author integrals are singular on the light cone; select a flank"
    )
args.cache_directory.mkdir(parents=True, exist_ok=True)
started = time.monotonic()
used = []


def interrupted(signum, frame):
    # Raising lets subprocess.run kill and reap its current upstream child.
    raise RuntimeError(f"author reference interrupted by signal{signum}")


signal.signal(signal.SIGTERM, interrupted)
signal.signal(signal.SIGHUP, interrupted)


def raw(E, k, component):
    parameters = [SOURCE_SHA, E.hex(), k.hex(), args.tolerance.hex(), component]
    key = hashlib.sha256(json.dumps(parameters).encode()).hexdigest()
    path = args.cache_directory / f"{key}.json"
    if path.exists():
        item = json.loads(path.read_text())
        if item.get("parameters") != parameters:
            raise RuntimeError(f"cache identity mismatch: {path}")
        item = dict(item, cache_hit=True)
    else:
        remaining = args.timeout - (time.monotonic() - started)
        if remaining <= 0:
            raise TimeoutError("total author-reference wall-time exhausted")
        command = [str(args.binary), str(E), str(k), str(args.tolerance), component]
        begin = time.monotonic()
        result = subprocess.run(
            command, capture_output=True, text=True, timeout=remaining, check=True
        )
        fields = list(map(float, result.stdout.strip().splitlines()[-1].split(",")))
        lo, nlo = fields[3:5]
        if not math.isfinite(lo) or not math.isfinite(nlo):
            raise RuntimeError(f"nonfinite author result: {command}")
        item = {
            "parameters": parameters,
            "E": E,
            "k": k,
            "tolerance": args.tolerance,
            "component": component,
            "lo": lo,
            "nlo": nlo,
            "command": command,
            "driver_binary_sha256": hashlib.sha256(
                args.binary.read_bytes()
            ).hexdigest(),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "wall_seconds": time.monotonic() - begin,
            "cache_hit": False,
        }
        temporary = path.with_suffix(f".{os.getpid()}.tmp")
        temporary.write_text(json.dumps(item, indent=2) + "\n")
        temporary.replace(path)
    used.append(item)
    return item


try:
    E, k = args.energy, args.momentum
    result = {
        "status": "ok",
        "E": E,
        "k": k,
        "tolerance": args.tolerance,
        "basis": args.basis,
        "source_archive_sha256": SOURCE_SHA,
        "units": "T=1, Nc=3, unit charge; NLO is g^2 coefficient",
        "numerical_caveat": "upstream GSL statuses are not exposed; tolerance/extrapolation differences are empirical checks",
    }
    if k == 0:
        controls = [raw(E, E * args.rest_fraction / 2**i, "vector") for i in range(3)]
        large = (4 * controls[1]["nlo"] - controls[0]["nlo"]) / 3
        small = (4 * controls[2]["nlo"] - controls[1]["nlo"]) / 3
        result.update(
            rhoV_lo=3 * E * E * math.tanh(E / 4) / (4 * math.pi),
            rhoV_nlo=small,
            method="vector k^2 extrapolation then exact rest isotropy",
            rest_extrapolation_difference=abs(small - large),
            rest_quartic_extrapolation=(16 * small - large) / 15,
            rest_controls=controls,
        )
        for component in ("rhoT", "rhoL"):
            result[component + "_lo"] = result["rhoV_lo"] / 3
            result[component + "_nlo"] = result["rhoV_nlo"] / 3
    else:
        if args.basis in ("vector", "both"):
            vector = raw(E, k, "vector")
            result.update(rhoV_lo=vector["lo"], rhoV_nlo=vector["nlo"])
        if args.basis in ("density", "both"):
            density = raw(E, k, "density")
            factor = -(E * E - k * k) / (k * k)
            result.update(
                minus_ImPi00_lo=density["lo"],
                minus_ImPi00_nlo=density["nlo"],
                rhoL_lo=factor * density["lo"],
                rhoL_nlo=factor * density["nlo"],
            )
        if args.basis == "both":
            for order in ("lo", "nlo"):
                result["rhoT_" + order] = (
                    result["rhoV_" + order] - result["rhoL_" + order]
                ) / 2
        result["method"] = (
            "unchanged upstream assembly and paper Eq2.5 polarization conversion"
        )
    for component in ("rhoV", "rhoT", "rhoL"):
        if component + "_lo" in result and component + "_nlo" in result:
            result[component + "_g3"] = (
                result[component + "_lo"] + 9 * result[component + "_nlo"]
            )
    result.update(
        raw=used,
        wall_seconds=time.monotonic() - started,
        cache_hits=sum(item["cache_hit"] for item in used),
        cpu_workers=1,
        cpu=args.cpu,
    )
    code = 0
except (
    subprocess.CalledProcessError,
    subprocess.TimeoutExpired,
    TimeoutError,
    InterruptedError,
    RuntimeError,
    ValueError,
    IndexError,
    KeyboardInterrupt,
) as error:
    result = {
        "status": "incomplete",
        "E": args.energy,
        "k": args.momentum,
        "tolerance": args.tolerance,
        "error": str(error),
        "raw": used,
        "wall_seconds": time.monotonic() - started,
    }
    code = 3
rendered = json.dumps(result, indent=2) + "\n"
if args.output:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
print(rendered, end="")
sys.exit(code)
