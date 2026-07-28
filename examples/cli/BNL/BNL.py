#!/usr/bin/env python3
"""Compute BNL complex-result ratios with first-order error propagation."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ComplexResult:
    re: float
    im: float
    re_err: float
    im_err: float

    @property
    def value(self) -> complex:
        return complex(self.re, self.im)


def read_integration_result(path: Path) -> ComplexResult:
    with path.open() as f:
        data = json.load(f)
    slot = data["slots"][0]
    integral = slot["integral"]
    return ComplexResult(
        re=float(integral["result"]["re"]),
        im=float(integral["result"]["im"]),
        re_err=float(integral["error"]["re"]),
        im_err=float(integral["error"]["im"]),
    )


def ratio_with_uncertainty(num: ComplexResult, den: ComplexResult) -> ComplexResult:
    x = num.value
    y = den.value
    if y == 0:
        raise ZeroDivisionError("Cannot form a ratio with a zero denominator")

    value = x / y
    inv_y = 1.0 / y
    d_x_re = inv_y
    d_x_im = 1j * inv_y
    d_y_re = -x / (y * y)
    d_y_im = -1j * x / (y * y)

    derivatives = [
        (d_x_re, num.re_err),
        (d_x_im, num.im_err),
        (d_y_re, den.re_err),
        (d_y_im, den.im_err),
    ]
    re_var = sum((deriv.real * sigma) ** 2 for deriv, sigma in derivatives)
    im_var = sum((deriv.imag * sigma) ** 2 for deriv, sigma in derivatives)
    return ComplexResult(value.real, value.imag, re_var**0.5, im_var**0.5)


def parse_result_arg(raw: str) -> tuple[str, ComplexResult]:
    try:
        name, re, im, re_err, im_err = raw.split(":")
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "manual results must be NAME:RE:IM:RE_ERR:IM_ERR"
        ) from exc
    return name, ComplexResult(float(re), float(im), float(re_err), float(im_err))


def _ndec(value: float, offset: int) -> int:
    ans = int(offset - math.log10(value))
    thresholds = [0.5, 9.5, 99.5]
    if ans > 0 and value * 10.0**ans >= thresholds[offset]:
        ans -= 1
    return max(ans, 0)


def _normalize_exponent(exponent: str) -> str:
    return str(int(exponent))


def _format_uncertainty(mean: float, error: float) -> str:
    value = mean
    delta = abs(error)

    if math.isnan(value) or math.isnan(delta):
        return f"{value:e} +/- {delta:e}"
    if math.isinf(delta):
        return f"{value:e} +/- inf"
    if value == 0.0 and not (1e-4 <= delta < 1e5):
        if delta == 0.0:
            return "0(0)"
        mantissa, exponent = f"{delta:.1e}".split("e")
        return f"0.0({mantissa})e{_normalize_exponent(exponent)}"
    if value == 0.0:
        if delta >= 9.95:
            return f"0({delta:.0f})"
        if delta >= 0.995:
            return f"0.0({delta:.1f})"
        decimals = _ndec(delta, 2)
        return f"{value:.{decimals}f}({delta * 10.0**decimals:.0f})"
    if delta == 0.0:
        mantissa, exponent = f"{value:e}".split("e")
        exponent = _normalize_exponent(exponent)
        return f"{mantissa}(0)e{exponent}" if exponent != "0" else f"{mantissa}(0)"
    if delta > 1e4 * abs(value):
        return f"{value:.1e} +/- {delta:.2e}"
    if abs(value) >= 1e6 or abs(value) < 1e-5:
        exponent = math.floor(math.log10(abs(value)))
        scale = 10.0**exponent
        mantissa = _format_uncertainty(value / scale, delta / scale)
        return f"{mantissa}e{exponent}"
    if delta >= 9.95:
        if abs(value) >= 9.5:
            return f"{value:.0f}({delta:.0f})"
        decimals = _ndec(abs(value), 1)
        return f"{value:.{decimals}f}({delta:.{decimals}f})"
    if delta >= 0.995:
        if abs(value) >= 0.95:
            return f"{value:.1f}({delta:.1f})"
        decimals = _ndec(abs(value), 1)
        return f"{value:.{decimals}f}({delta:.{decimals}f})"

    decimals = max(_ndec(abs(value), 1), _ndec(delta, 2))
    return f"{value:.{decimals}f}({delta * 10.0**decimals:.0f})"


def format_uncertainty(mean: float, error: float) -> str:
    formatted = _format_uncertainty(mean, error)
    return formatted if math.copysign(1.0, mean) < 0.0 else f"+{formatted}"


def format_result(result: ComplexResult) -> str:
    return f"{format_uncertainty(result.re, result.re_err)}  {format_uncertainty(result.im, result.im_err)} i"


def default_result_paths(root: Path) -> dict[str, Path]:
    workspaces = root / "workspaces"
    return {
        "R1": workspaces
        / "current_bnl_R1_mt400_threshold_10M_5c"
        / "integration_result.json",
        "R2": workspaces
        / "final_bnl_R2_mt450_threshold_2M_5c"
        / "integration_result.json",
        "R3": workspaces
        / "final_bnl_R3_mt510_no_threshold_3M_5c"
        / "integration_result.json",
        "R4": workspaces
        / "final_bnl_R4_mt550_no_threshold_3M_5c"
        / "integration_result.json",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--result",
        action="append",
        type=parse_result_arg,
        default=[],
        metavar="NAME:RE:IM:RE_ERR:IM_ERR",
        help="Provide one manual complex result. Repeat for R1, R2, R3 and R4.",
    )
    parser.add_argument(
        "--integration-result",
        action="append",
        default=[],
        metavar="NAME=PATH",
        help="Read one GammaLoop integration_result.json. Repeat for R1..R4.",
    )
    parser.add_argument(
        "--json-output",
        type=Path,
        help="Optional path where the ratio summary is written as JSON.",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
    results: dict[str, ComplexResult] = {}

    if not args.result and not args.integration_result:
        for name, path in default_result_paths(root).items():
            results[name] = read_integration_result(path)

    for name, result in args.result:
        results[name] = result

    for raw in args.integration_result:
        if "=" not in raw:
            raise SystemExit("--integration-result entries must be NAME=PATH")
        name, raw_path = raw.split("=", 1)
        results[name] = read_integration_result(Path(raw_path))

    missing = [name for name in ["R1", "R2", "R3", "R4"] if name not in results]
    if missing:
        raise SystemExit(f"Missing required results: {', '.join(missing)}")

    ratios = {
        "R2/R1": ratio_with_uncertainty(results["R2"], results["R1"]),
        "R3/R1": ratio_with_uncertainty(results["R3"], results["R1"]),
        "R4/R1": ratio_with_uncertainty(results["R4"], results["R1"]),
    }

    print("Input results:")
    for name in ["R1", "R2", "R3", "R4"]:
        print(f"  {name}: {format_result(results[name])}")

    print("\nRatios:")
    for name, ratio in ratios.items():
        print(f"  {name}: {format_result(ratio)}")

    if args.json_output:
        payload = {
            "inputs": {name: result.__dict__ for name, result in results.items()},
            "ratios": {name: ratio.__dict__ for name, ratio in ratios.items()},
        }
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
