#!/usr/bin/env python3
"""Extract the first Arb identity/Euler CT pair from an existing display log.

This reads a completed trace only; it never loads a GammaLoop state. The center
tokens are round-trippable binary64 output, whereas r/rstar and alpha residuals
are native Arb decimal output. Original line numbers are retained as provenance.
"""
import gzip
import hashlib
import json
import re
import sys
from decimal import Decimal, localcontext
from pathlib import Path

path = Path(sys.argv[1])
raw = gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes()
lines = raw.decode().splitlines()
start = next(i for i, line in enumerate(lines) if "ArbPrec parameterization succeeded" in line)
end = next(i for i in range(start + 1, len(lines)) if "parameterization succeeded" in lines[i])
rotations = [i for i in range(start, end) if "Evaluating rotation:" in lines[i]]
assert len(rotations) == 2
records = []
for begin, stop in zip(rotations, rotations[1:] + [end]):
    center_start = next(i for i in range(begin, stop) if "left overlap structure:" in lines[i])
    right_start = next(i for i in range(center_start + 1, stop) if "right overlap structure:" in lines[i])
    center_rows = []
    for i in range(center_start, right_start):
        cells = [cell.strip() for cell in lines[i].split("│")]
        if len(cells) == 6 and cells[1].isdigit():
            center_rows.append({"index": int(cells[1]), "tokens": cells[2:5], "line": i + 1})
    membership = next(line.split(":", 1)[1].strip() for line in lines[center_start:right_start]
                      if line.startswith("existing esurfaces in group:"))
    subspace_index = next(i for i in range(right_start, stop) if "subspace: SubspaceData" in lines[i])
    parent = re.search(r"lmb: LmbIndex\((\d+)\)", lines[subspace_index]).group(1)
    active = re.search(r"lmb_indices: (\[.*?\])", lines[subspace_index]).group(1)
    alpha = []
    for i in range(subspace_index, stop):
        if "LU left evaluator input:" in lines[i]:
            break
        if "alpha solution:" in lines[i]:
            residual = re.search(r"error_of_function: F\(VarFloat \{ float: ([^ }]+)", lines[i]).group(1)
            alpha.append({"residual": residual, "line": i + 1})
    thresholds = {}
    for i in range(begin, stop):
        match = re.search(r"LU (left|right) evaluator input: cut_group_id=0, (?:left|right)_threshold_id=(\d+).*?, r=([^,]+), rstar=(\S+)", lines[i])
        if match:
            side, index, radius, star = match.groups()
            thresholds[f"{side}{index}"] = {"r": radius, "rstar": star, "line": i + 1}
    assert len(center_rows) == 4 and len(alpha) == 4 and len(thresholds) == 4
    with localcontext() as context:
        context.prec = 350
        norm = sum(Decimal.from_float(float(token)) ** 2
                   for row in center_rows for token in row["tokens"]).sqrt()
    records.append({"rotation": lines[begin].split("Evaluating rotation:", 1)[1].strip(),
                    "rotation_line": begin + 1, "center_rows": center_rows,
                    "center_norm": str(norm), "members": membership,
                    "parent_lmb_index": parent, "active_indices": active,
                    "alpha_solve_residuals_left0_left1_right0_right1": alpha,
                    "thresholds": thresholds})

with localcontext() as context:
    context.prec = 350
    first, second = records
    differences = {"center_norm_GeV": str(Decimal(second["center_norm"]) - Decimal(first["center_norm"]))}
    differences["center_norm_relative"] = str(Decimal(differences["center_norm_GeV"]) / Decimal(first["center_norm"]))
    for key in first["thresholds"]:
        differences[key] = {field: str(Decimal(second["thresholds"][key][field]) - Decimal(first["thresholds"][key][field]))
                            for field in ["r", "rstar"]}
result = {"source_sha256": hashlib.sha256(raw).hexdigest(), "source_lines": len(lines),
          "scope": "first explicit Arb physical call, CT-on identity/Euler, physical cut3/group0; trace timing is not performance evidence",
          "records": records, "Euler_minus_identity": differences,
          "same_members_parent_active": all(first[key] == second[key] for key in ["members", "parent_lmb_index", "active_indices"])}
Path(sys.argv[2]).write_text(json.dumps(result, indent=2) + "\n")
