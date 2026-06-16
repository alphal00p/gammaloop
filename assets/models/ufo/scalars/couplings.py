# type: ignore
from __future__ import absolute_import
from .object_library import all_couplings, Coupling


# L_int = -lam * product(phi_a**n_a) / product(n_a!) gives this vertex rule.
SCALAR_COUPLING = Coupling(
    name="SCALAR_COUPLING", value="-lam*complex(0,1)", order={"QCD": 1}
)
