"""Smoke test for an installed combined Symbolica community wheel.

This is intentionally not a workspace test: `feynkit-py` is an rlib, while the
external Symbolica community host builds and installs the single native wheel.
"""

import sys

from symbolica import Expression, S
from symbolica.community import feynkit


assert feynkit.__name__ == "symbolica.community.feynkit"
assert "symbolica.community.feynkit_native" in sys.modules
assert "_gammaloop" not in sys.modules

for exported_type in (
    feynkit.Model,
    feynkit.Generator,
    feynkit.Process,
    feynkit.FeynmanDiagram,
    feynkit.TensorReducer,
    feynkit.CffGenerator,
    feynkit.FourMomentum,
):
    assert exported_type.__module__ == "symbolica.community.feynkit"

momentum = feynkit.ThreeMomentum(3.0, 4.0, 0.0).on_shell()
assert momentum.components() == (5.0, 3.0, 4.0, 0.0)
assert Expression.__module__ == "symbolica.core"

dimension = S("feynkit_install_test::D")
mu = S("feynkit_install_test::mu")
nu = S("feynkit_install_test::nu")
mink = S("spenso::mink")
dot = S("spenso::dot")
loop_momentum = S("feynkit_install_test::k")
projector_momentum = S("feynkit_install_test::p")
tensor_input = (
    loop_momentum(mink(dimension, mu))
    * loop_momentum(mink(dimension, nu))
    * projector_momentum(mink(dimension, mu))
    * projector_momentum(mink(dimension, nu))
)
tensor_reduced = (
    feynkit.TensorReducer(dimension)
    .with_integrated_head("feynkit_install_test::k")
    .reduce(tensor_input)
)
tensor_expected = (
    dot(loop_momentum(mink(dimension)), loop_momentum(mink(dimension)))
    * dot(
        projector_momentum(mink(dimension)),
        projector_momentum(mink(dimension)),
    )
    / dimension
)
assert tensor_reduced == tensor_expected
