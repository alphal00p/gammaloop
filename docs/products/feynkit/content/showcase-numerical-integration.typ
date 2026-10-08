#import "../../shared.typ": source-link

#let showcase-numerical-integration = [
= Numerical integration with CFF and LTD

The #source-link("examples/notebooks/feynkit/08_numerical_cff_ltd_marimo.py", label: "numerical integration notebook")
follows one generated scalar triangle from its momentum routing to a numerical
integral. Its three LTD residues and generated CFF expression are checked for
exact symbolic equality before comparing floating-point behavior.

The *Two- and three-loop surfaces* gallery adds four symbolic examples: a
two-loop double box and crossed double box, a three-loop triple box, and a
three-loop vacuum tetrahedron. Its diagram picker updates both native
representations through `integrate_energy(method="cff")` and
`integrate_energy(method="ltd")`. The compact DOT definitions in the notebook
are editable; `lmb_id` sets the order of the reference loop edges.
Shift-click CFF factors to compare multiple surface regions, then inspect the
corresponding LTD cuts and E/H surface pairs. Numerical sampling and the OneLOop
comparison continue to use the triangle.

== What to explore

`diagram.integrate_energy(method="ltd")` returns a native `symbolica.community.hepkit.LtdRepresentation`.
`diagram.integrate_energy(method="cff")` returns `symbolica.community.hepkit.CffRepresentation`.
These are the two energy representations described together in the
#link("reference/interfaces/")[interface guide], including their conversion defaults,
routing conventions and supported inputs.
Its notebook display starts with the residue sum and uses Linnet's graph layout,
half-edge shading, and viewport. Drag to pan, use Ctrl/Command-scroll or `+`/`-`
to zoom, and `0` to fit; the camera stays fixed while inspecting other residues.
Click a cut gap to visit the next residue cutting that edge; click a tree edge
to switch its surface pair. Hover a cut to trace the fundamental cycle relative
to the ordered reference momentum basis. The shaded component stays fixed
within a pair. A gap pointing out of it denotes a positive pole, and a gap
pointing into it denotes a negative pole; solid arrows retain the momentum
orientation. Details contains the routing, cut vector, and boundary-sign product.

`residue.to_expression()` and `ltd.to_expression()` return rank-zero
`TensorExpression` values for the native scalar terms and their sum, including
the on-shell energy factors and the clockwise
$d q^0/(2 pi i)$ contour sign. They omit numerators, graph weights and the spatial
measure. The numerical cells apply the notebook's measure convention explicitly.

- Compare the individual LTD cuts and their sum against CFF along a momentum
  ray. Symbolica's `FunctionDefinition` supplies the on-shell square roots to
  the same evaluator at double and arbitrary precision, separating algebraic
  equality from numerical stability.
- Change the internal mass and external kinematics to inspect cancelling
  H-surfaces and physical E-surfaces in a momentum slice. Their equations are
  extracted from the generated LTD denominators and normalized to expanded
  affine forms before classifying their energy signs. The default spacelike
  example has visible H-contours. For an E-contour, select on-shell external
  legs and $m=0.25$; on-shell legs with $m>0.5$ have no real zero contours.
- Compare Cartesian and spherical coordinate maps, uniform sampling, adaptive
  continuous grids, and adaptive discrete channels at equal evaluation budgets.
- Inspect convergence, statistical uncertainty, and weight distributions
  against an independent OneLOop scalar-triangle value. Repeat with another
  seed and with massless spacelike kinematics.

The numerical observable has numerator one and omits the model's couplings and
graph weights. Each representation includes the spatial loop measure once;
the OneLOop comparison uses $I=-C_0/(16 pi^2)$.
All batches contribute to Symbolica's native mean and error estimates. Adaptive
grids update between batches; uniform grids keep their initial density. Each
method receives the same total evaluation budget, and the retained sample
weights are used only for the distribution plots. Coordinate maps and their
Jacobian factors are evaluated by Symbolica as well. Physical-threshold
kinematics remain available for the geometry plots; their real-axis integration
requires a separate treatment and is disabled in this notebook.

== Run locally

This showcase requires a full native Symbolica community host with
`symbolica.community.hepkit`, `hepkit.oneloop`, and Symbolica's numerical
integrator, plus Marimo, NumPy, and Matplotlib. Use that host's interpreter:

// docs-example: syntax
```sh
just notebook feynkit/08_numerical_cff_ltd_marimo /path/to/python
```

Use an optimized native build, such as the community host's
`release-performance` profile, when comparing timings.

The smaller browser notebook host bundled with this repository does not
include OneLOop; this example is offered as a native notebook rather than a
browser embed. Select *Run integration comparison* to begin sampling.

Start with the #link("guides/showcases/cff/")[CFF structure showcase] for an
introduction to orientations and surface metadata, or return to the
#link("guides/showcases/")[notebook gallery].
]
