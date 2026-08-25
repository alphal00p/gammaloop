"""Smoke test for an installed combined Symbolica community wheel.

This is intentionally not a workspace test: `feynkit-py` is an rlib, while the
external Symbolica community host builds and installs the single native wheel.
"""

import sys

from symbolica import Expression
from symbolica.community import feynkit


assert feynkit.__name__ == "symbolica.community.feynkit"
assert "symbolica.community.feynkit_native" in sys.modules
assert "_gammaloop" not in sys.modules

for exported_type in (
    feynkit.Model,
    feynkit.Generator,
    feynkit.Process,
    feynkit.FeynmanDiagram,
    feynkit.CffGenerator,
    feynkit.FourMomentum,
):
    assert exported_type.__module__ == "symbolica.community.feynkit"

momentum = feynkit.ThreeMomentum(3.0, 4.0, 0.0).on_shell()
assert momentum.components() == (5.0, 3.0, 4.0, 0.0)
assert Expression.__module__ == "symbolica.core"
