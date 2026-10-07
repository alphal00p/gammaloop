#import "../../shared.typ": source-link

#let showcase-numerical-integration = [
= Numerical integration with CFF and LTD

The #source-link("examples/notebooks/feynkit/08_numerical_cff_ltd_marimo.py", label: "numerical integration notebook")
follows one generated scalar triangle from its momentum routing to a numerical
integral. Its three LTD residues and generated CFF expression are checked for
exact symbolic equality before comparing floating-point behavior.

== What to explore

- Compare the individual LTD cuts and their sum against CFF along a momentum
  ray. A high-precision evaluation, including the on-shell square roots,
  separates algebraic equality from numerical stability.
- Change the internal mass and external kinematics to inspect cancelling
  H-surfaces and physical E-surfaces in a momentum slice.
- Compare Cartesian and spherical coordinate maps, uniform sampling, adaptive
  continuous grids, and adaptive discrete channels at equal evaluation budgets.
- Inspect convergence, production-sample uncertainty, and weight distributions
  against an independent OneLOop scalar-triangle value. Repeat with another
  seed and with massless spacelike kinematics.

The numerical observable has numerator one and omits the model's couplings and
graph weights. Each representation includes the spatial loop measure once;
the OneLOop comparison uses $I=-C_0/(16 pi^2)$.
Adaptive runs train for two batches, freeze their sampling grids, and estimate
the integral from independent production samples. Training evaluations count
toward the common budget. Physical-threshold kinematics remain available for
the geometry plots; their real-axis integration requires a separate treatment
and is disabled in this notebook.

== Run locally

This showcase requires a full native Symbolica community host with
`symbolica.community.hepkit`, `hepkit.oneloop`, and Symbolica's numerical
integrator, plus Marimo, NumPy, and Matplotlib. Use that host's interpreter:

// docs-example: syntax
```sh
just notebook feynkit/08_numerical_cff_ltd_marimo /path/to/python
```

The smaller browser notebook host bundled with this repository does not
include OneLOop; this example is offered as a native notebook rather than a
browser embed. Select *Run integration comparison* to begin sampling.

Start with the #link("guides/showcases/cff/")[CFF structure showcase] for an
introduction to orientations and surface metadata, or return to the
#link("guides/showcases/")[notebook gallery].
]
