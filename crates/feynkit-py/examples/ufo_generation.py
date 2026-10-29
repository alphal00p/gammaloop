"""Generate tree-level e- e+ -> mu- mu+ diagrams from a UFO model.

Usage: python ufo_generation.py [UFO_MODEL_DIR]

The model directory defaults to $FEYNKIT_UFO_MODEL, then to the Standard Model
tracked at assets/models/ufo/sm in this repository.
"""

import os
import sys
from pathlib import Path

from symbolica.community.hepkit import UfoLoader

REPOSITORY_SM = Path(__file__).resolve().parents[3] / "assets/models/ufo/sm"

model_dir = Path(
    sys.argv[1]
    if len(sys.argv) > 1
    else os.environ.get("FEYNKIT_UFO_MODEL", REPOSITORY_SM)
)
loaded = UfoLoader().load(model_dir)
process = loaded.model.process(["e-", "e+"], ["mu-", "mu+"])
result = process.generate_diagrams()

for diagram in result.diagrams:
    print(diagram.to_dot())
