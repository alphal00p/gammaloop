#import "../../shared.typ": catalog-contract

#let api = [
= Rust and Python interfaces

#catalog-contract(
  rust-scope: "feynkit, feynkit-model, feynkit-ufo, feynkit-kinematics, feynkit-graph, feynkit-amplitude, feynkit-generator, feynkit-cff, feynkit-tensor",
  python-scope: "symbolica.community.hepkit",
)

== Rust ownership map

- #link("reference/rust/feynkit/")[`feynkit`] re-exports the focused crates through feature gates.
- #link("reference/rust/feynkit_model/")[`feynkit-model`] owns validated models, parameter cards,
  typed IDs, and the explicit recomputation protocol.
- #link("reference/rust/feynkit_ufo/")[`feynkit-ufo`] loads UFO models through a caller-owned Python
  interpreter and `ufo_model_loader`.
- #link("reference/rust/feynkit_kinematics/")[`feynkit-kinematics`] owns momenta, transformations,
  signatures, helicities, symbolic scalar products, Mandelstam substitutions, and generalized-kt clustering.
- #link("reference/rust/feynkit_graph/")[`feynkit-graph`] owns finalized diagrams, cuts, routing,
  serialization, integral-family bases, and Linnest source output.
- #link("reference/rust/feynkit_amplitude/")[`feynkit-amplitude`] owns coherent symbolic amplitudes,
  physical conjugation, squared amplitudes, and shared external-state spin/color sums.
- #link("reference/rust/feynkit_generator/")[`feynkit-generator`] owns process selectors,
  generation options, progress/cancellation, and diagram grouping.
- #link("reference/rust/feynkit_cff/")[`feynkit-cff`] owns CFF expressions, surface arenas,
  orientations, residues, and topology conversion.
- #link("reference/rust/feynkit_tensor/")[`feynkit-tensor`] owns covariant tensor reduction, contraction
  orbits, and exact coefficient tables.

The facade's default features are `amplitude`, `cff`, `generator`, `graph`, `kinematics`, `model`, and
`tensor`. `ufo` is opt-in. Native Rustdoc for this product includes that optional interface;
applications choose their features explicitly. When using `feynkit-model` in isolation, enable
its `native` feature to select the GMP integer and MPFR floating-point backends. The kinematics crate selects those backends through its default `native` feature;
use `default-features = false` and `wasm` for the portable backends.
The APIs shown here correspond to the source
revision of this manual, so pin matching dependencies when reproducing a calculation.

== Python ownership map

The #link("reference/python/feynkit-community/")[generated Python reference] covers native
classes, properties, signatures, examples, and error types. The primary owners are `Model`,
`Process`, `SnailFilterOptions`, `NumeratorGrouping`, `FeynmanDiagram`, `Amplitude`, `SquaredAmplitude`, `AmplitudeLeg`, `Subgraph`, `TensorReducer`,
`CffRepresentation`, `LtdRepresentation`, `Kinematics`, `IntegralFamily`, `JetDefinition`,
and the optional `UfoLoader`. Expressions cross the boundary as Symbolica values.

Use the #link("quickstart/python/")[quickstart] for an installed host or the
#link("guides/community-host/")[community-host guide] when packaging the module. The Python
adapter does not publish a standalone wheel. Exception classes keep model, generation, graph,
amplitude, CFF, kinematics, and tensor-reduction failures distinguishable.

== Three-dimensional energy representations

`diagram.integrate_energy(method="cff")` returns a `CffRepresentation`;
`diagram.integrate_energy(method="ltd")` returns a `LtdRepresentation`.
Both organize the result of symbolic loop-energy integration. The remaining
spatial loop momenta are still variables: constructing either result does not
perform the numerical integral or insert a diagram numerator, couplings or
graph weight.

`integrate_energy` is the only construction entry point. Its `method` argument
is required. The Python stub uses `Literal` overloads: a known method gives a
single concrete return type; a runtime choice annotated `Literal["cff", "ltd"]`
gives `CffRepresentation | LtdRepresentation`. Unknown methods raise
`ValueError`. CFF-specific `max_orientations`, `fixed_orientations`,
`contracted_edges` and `initial_state_edges` are keyword options on this same
method; supplying them for LTD raises `ValueError`.

```python
from symbolica.community import hepkit as hep

diagram = (
    hep.Model.phi3().process(["phi"], ["phi", "phi"])
    .generate_diagrams(loops=1, max_vertices=3, maximum_bridges=0, progress=None)
    .diagrams[0]
)
cff = diagram.integrate_energy(method="cff")
ltd = diagram.integrate_energy(method="ltd")

orientation = cff.orientations[0]
family = orientation.families[0]  # CrossFreeFamily: one denominator family
residue = ltd.residues[0]        # LtdResidue: one retained signed assignment
assert sum(r.to_expression() for r in ltd.residues) == ltd.to_expression()
```

Both results expose `diagram`, `surfaces`, `on_shell_energies`, `report`, and
`to_expression()`. Displaying a result opens its interactive explorer. Displaying
an orientation shows only its family sum and fixed energy-flow graph; a family
shows only its contribution and circlings. A residue shows its signed expression
and fixed cut assignment, with selectable surface pairs. These smaller displays
keep their graphs collapsed initially and do not navigate to sibling objects.
Surfaces, factors, pairs and on-shell energies have compact mathematical displays;
generation reports have matching statistics tables.
All inspection objects are read-only. CFF reports candidate/acyclic orientations,
unfolded terms and interned surfaces; LTD reports residues, distinct trees,
unfolded terms and interned surfaces. `len(cff)` counts unfolded denominator
products and `len(ltd)` counts retained signed residues.

=== Expressions and on-shell energies

Both `to_expression()` methods return the scalar energy integral as a rank-zero
`TensorExpression`, including
contour coefficients and on-shell factors for one $d q^0/(2 pi i)$ contour
closed below per loop. The spatial measure, numerator, couplings and graph
weights remain separate. CFF multiplies its denominator sum by
$product_(e in "uncontracted internal") (-1)/(2 E_e)$; LTD retains its native
residue coefficients and energy factors. There is no separate normalization flag.
Use the same explicit remaining measure for either result.

Orientations, families, residues, surfaces, factors and on-shell energies also
convert to `TensorExpression`. This subclass of Symbolica's `Expression` provides
the shared mathematical renderer and paged notebook output. Its own
`to_expression()` method returns an ordinary Symbolica expression when needed.
Tensor equality includes the ordered tensor interface: use `is_zero()` for a
zero test, or convert to an ordinary expression before comparing with raw
Symbolica expressions. Algebraic transformations that preserve the tensor
interface keep the richer display; inherited transformations may return a plain
`Expression`.

Surface expansion defaults to `True` and replaces result-local surface
placeholders with affine energy combinations. It leaves
`gammalooprs::OSE(edge_id)` symbolic. CFF and LTD use the same OSE symbol for the
same physical diagram edge. `rep.on_shell_energies` maps physical edge IDs to
`OnShellEnergy` objects: `symbol` returns OSE, and `to_expression()` returns the
positive square root of the routed spatial momentum squared plus the symbolic
mass squared. Routing uses `gammalooprs::K` and `gammalooprs::P`, indexed by
positions in `loop_momentum_basis.loop_edges` and `.external_edges`, with
`spenso::cind(1..3)` Cartesian components. OSE itself remains keyed by physical
edge ID. These definitions work directly with Symbolica replacement and
custom function evaluation:

```python
from symbolica import FunctionDefinition
energy_functions = [
    FunctionDefinition(e.symbol, [], e.to_expression())
    for e in cff.on_shell_energies.values()
]
```

Use those same definitions for either representation; choose coordinate and
parameter values through Symbolica. This keeps surface expansion separate from
square-root expansion. External energy components remain
`gammalooprs::Q(edge_id, spenso::cind(0))`. CFF retains external-edge energies,
while LTD applies the stored routing; impose the same external
momentum-conservation relations when comparing their expanded expressions.
The #link("guides/showcases/numerical-integration/")[numerical showcase] uses
these shared energies and verifies the scalar triangle identity exactly.

=== Structured inspection

`EnergySurface` exposes `kind` (`"E"` or `"H"`), `symbol`, `index`, exact signed
`energy_coefficients`, `external_shift`, `constant`, and `to_expression()`.
Coefficients are keyed by physical edge ID. Surface indices and symbols are local
to a representation, so equal indices in different results need not mean the same
surface. `vertices` records a CFF circling; LTD regions belong to occurrences on
a particular tree rather than a globally stored surface.

`orientation.families` contains `CrossFreeFamily` objects. Each family's
`factors` are `SurfaceFactor` objects carrying a stored `surface`, a relative
`sign` and a positive `power`. A factor's expression is the denominator
$("sign" dot "surface")^"power"$, before inversion. Family and orientation
expression conversions include their common energy prefactor, so summing the
families gives the orientation and summing orientations gives the representation.
Raised-surface grouping, pole coefficients and analytic residues remain available;
residues now use that same energy normalization by default.

A CFF orientation's `edge_signs` maps physical edge IDs to `+1` when the flow
agrees with the stored diagram arrow, `-1` when reversed, and `None` for an
undirected edge. These signs describe the orientation; a surface coefficient
also depends on the boundary crossing. `edge_orientations` supplies the same
directions as `"default"`, `"reversed"` or `"undirected"` strings.

An LTD residue exposes `tree_edges`, `cut_edges`, `pole_signs`, `energy_map` and
`surface_pairs`. Pole signs mean $q_e^0 = sigma_e E_e$ relative to the stored
momentum arrows. `energy_map` gives the resulting $q_e^0$ for each internal edge.
`surface_pairs[edge].minus` and `.plus` retain the factors $q_e^0-E_e$ and
$q_e^0+E_e$ respectively. Their occurrence signs are separate from the cut pole
signs. Pairs are reported for surviving uncut propagators; repeated propagators
can contribute several algebraic terms to one residue. Its `to_expression()`
retains the exact coefficients and multiplicities.

=== Routing, pole signs and surface pairs

For LTD, the ordered `diagram.loop_momentum_basis.loop_edges` supplies the
reference energy coordinates and integration order. Use
`diagram.with_loop_momentum_edges([...])` to construct a diagram with another
ordered basis before calling `integrate_energy(method="ltd")`. Reordering the reference
basis can change individual retained residues and their pole signs.

For a retained cut set $C$, write the routed internal energies as
$q^0 = S k^0 + p^0$, with $p^0$ the known external shifts. The native residue
algorithm selects a pole assignment $sigma_c$ from the successive lower
contours, then solves
$ k^0 = S_C^(-1) (sigma_C E_C - p_C^0). $
Here $sigma_C E_C$ means componentwise multiplication. Each uncut tree edge $e$
therefore has one routed energy
$ a_e = S_e S_C^(-1) (sigma_C E_C - p_C^0) + p_e^0 $
and the two denominator factors $a_e - E_e$ and $a_e + E_e$. The cut-energy
coefficients are fixed by the same pole assignment throughout that residue;
switching within a pair changes the tree-edge factor. For repeated propagators,
higher-order residues can produce several terms in one signed cut assignment.

After factoring out an overall sign, a surface with on-shell energy coefficients
of one sign is an E-surface; mixed signs give an H-surface. This is an algebraic
classification, not a guarantee of a real zero for the selected kinematics.
A shaded graph component helps read the boundary incidence. Choosing its
complement reverses all boundary signs and does not select the other factor of
the pair. Pole signs come from the residue calculation, not the shading choice.

=== Inspection and supported inputs

Both classes display automatically through `_display_()` in Marimo and
`_repr_html_()` in other HTML notebooks. The expression comes first; expand its
explorer to inspect the original graph with shared Linnet layout, shading and
zoom. Browser selections do not change the symbolic result.

- CFF: click an orientation arrowhead, select a family, and click a surface to
  inspect its circling. Shift-click combines surface regions.
- LTD: click a cut gap to navigate residues containing that cut, and click an
  uncut tree edge to switch its surface pair. Hover a cut to trace its sign
  through the reference routing. Parallel arrows show momentum orientation;
  gap chevrons show the pole direction relative to the shaded component.

CFF accepts orientation constraints, contractions and selected subgraphs. Its
`max_orientations` limit bounds candidate enumeration; constrained output need
not represent the full unconstrained sum. LTD currently accepts complete
ordinary diagrams, with scalar numerator one, and supports repeated propagators.
Partial subgraphs and auxiliary cut edges are rejected. Neither constructor
accepts an arbitrary energy-dependent numerator.

=== CFF coefficients and generalized cuts

`diagram.subgraph(selection).integrate_energy(method="cff")` constructs the selected region's CFF;
`integrate_energy(method="cff")` accepts the diagram or its `Subgraph` view.
`result.to_expression()` includes on-shell factors and expands surfaces by
default, using the same contour convention as LTD. Set `expand_surfaces=False`
to retain canonical surface placeholders. Numerators, projectors, spatial
measure and global weights remain separate.

`result.raised_surface_groups(edge_representatives)` identifies equivalent
surfaces after mapping repeated propagators to their representative edge.
`result.pole_coefficients(group)` returns coefficients in increasing pole order;
these are not yet analytic residues. `result.residue` additionally requires the
integration variable, its root, the surface in that variable, and the complete
remaining coefficient. Raised-surface groups belong to the result that
created them; passing a group to another result is rejected. EnergySurface
objects carry their own affine definitions and can be inspected independently.

```python
q0 = Expression.symbol("q0")
cut = CutPropagator(q0, 2, power=3, normalization=1)
assert cut.apply(q0**2, q0) == -Expression.num(1) / 128
```

`CutPropagator` follows the generalized distributional cutting rules in
#link("https://arxiv.org/abs/2203.11038")[Local Unitarity, section 2.2].
Its default factor is $-2 pi i$ for a reciprocal propagator. The orientation
chooses the positive or negative energy root; the prescription specifies the
sign of the infinitesimal imaginary part. A raised cut differentiates both
the remaining integrand and the uncut energy factor before imposing the shell.
The energy-space symbol `Delta(n,x)` denotes the normalized divided-difference
distribution whose action is the derivative of order `n-1` divided by `(n-1)!`.
Route to an independent energy coordinate before calling `apply`.

]
