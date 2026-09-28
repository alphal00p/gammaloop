"""The electroweak fields and interactions use the xi=1 Feynman gauge."""

from .object_library import Propagator, all_propagators

Z = Propagator("Z_propFeynman", "-complex(0,1)*Metric(1,2)", "P(1)**2-MZ**2")
Wplus = Propagator("W+_propFeynman", "-complex(0,1)*Metric(1,2)", "P(1)**2-MW**2")
Wminus = Propagator("W-_propFeynman", "-complex(0,1)*Metric(1,2)", "P(1)**2-MW**2")
