"""Refresh the FeynCalc gallery requirement inventory from its published pages.

Only URLs and function names are saved, not copied notebook implementations.
End-to-end validation remains a separate, manually reviewed claim.
"""

import argparse
import json
import re
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime
from html import unescape
from pathlib import Path
from urllib.parse import urljoin
from urllib.request import urlopen

SOURCE = "https://feyncalc.github.io/examples"
REQUIREMENTS = {
    "spin_sums": {"FermionSpinSum", "DoPolarizationSums"},
    "symbolic_kinematics": {"SetMandelstam", "TrickMandelstam", "ExpandScalarProduct"},
    "dirac_algebra": {
        "DiracSimplify",
        "DiracTrace",
        "DiracSubstitute67",
        "DiracEquation",
    },
    "color_algebra": {"SUNSimplify"},
    "loop_families": {
        "FCLoopFindTopologies",
        "FCLoopFindTopologyMappings",
        "FCLoopFindIntegralMappings",
        "FCLoopApplyTopologyMappings",
    },
    "loop_family_completion": {"FCLoopBasisFindCompletion"},
    "sector_analysis": {
        "FCLoopFindSubtopologies",
        "FCLoopScalelessQ",
        "FCLoopPakScalelessQ",
    },
    "dependent_propagators": {
        "ApartFF",
        "FCApart",
        "FCLoopBasisOverdeterminedQ",
        "FCLoopRewriteOverdeterminedTopologies",
    },
    "tensor_reduction_external_momenta": {"FCLoopTensorReduce", "TID"},
    "ibp_kira": {"KiraRunReduction", "KiraCreateJobFile"},
    "ibp_fire": {"FIRERunReduction", "FIRECreateConfigFile"},
    "analytic_loop_integrals": {
        "PaXEvaluate",
        "PaXEvaluateUV",
        "PaVeUVPart",
        "PaXEvaluateUVIRSplit",
    },
    "phase_space_integration": {"Integrate"},
    "gamma5_scheme": {"FCSetDiracGammaScheme"},
    "light_cone": {"LightConePerpendicularComponent"},
    "eikonal_propagators": {"FCLoopReplaceQuadraticEikonalPropagators"},
    "feynman_parameters": {"FCFeynmanParametrize"},
}


def inspect(url):
    row = {"url": url, "retrieval": "inspected"}
    try:
        with urlopen(url, timeout=60) as response:
            page = response.read().decode()
        code = "\n".join(
            unescape(re.sub("<[^>]+>", "", block))
            for block in re.findall(r"<code[^>]*>(.*?)</code>", page, re.DOTALL)
        )
        functions = set(re.findall(r"\b([A-Z][A-Za-z0-9]+)\s*(?:\[|@|/@)", code))
        # Mathematica's postfix form is common in the gallery: expr // DiracSimplify.
        functions.update(re.findall(r"//\s*([A-Z][A-Za-z0-9]+)\b", code))
        row.update(
            requirements=[
                key for key, names in REQUIREMENTS.items() if functions & names
            ],
            functions=sorted(functions),
        )
    except OSError as error:
        row.update(
            retrieval="unavailable", requirements=[], functions=[], error=str(error)
        )
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--date", default=datetime.now(UTC).date().isoformat())
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1]
        / "docs/products/feynkit/feyncalc-gallery.json",
    )
    args = parser.parse_args()
    previous = (
        json.loads(args.output.read_text())
        if args.output.exists()
        else {"examples": []}
    )
    validations = {
        row["url"]: row.get("end_to_end_validation", "pending")
        for row in previous["examples"]
    }
    components = {
        row["url"]: row["validated_components"]
        for row in previous["examples"]
        if "validated_components" in row
    }
    with urlopen(SOURCE, timeout=60) as response:
        gallery = response.read().decode()
    links = list(
        dict.fromkeys(
            urljoin(SOURCE, href)
            for href in re.findall(r'href="([^"]*FeynCalcExamples[^"]*)"', gallery)
        )
    )
    if not links:
        raise RuntimeError(
            "the gallery contains no recognized example links; inventory retained"
        )
    with ThreadPoolExecutor(max_workers=4) as pool:
        rows = list(pool.map(inspect, links))
    for row in rows:
        row["end_to_end_validation"] = validations.get(row["url"], "pending")
        if row["url"] in components:
            row["validated_components"] = components[row["url"]]
    args.output.write_text(
        json.dumps(
            {"source": SOURCE, "inspected_on": args.date, "examples": rows}, indent=2
        )
        + "\n"
    )
    print(
        f"Inspected {len(rows)} links; {sum(r['retrieval'] == 'unavailable' for r in rows)} unavailable"
    )


if __name__ == "__main__":
    main()
