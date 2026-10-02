= CFF surface-cache ownership
<cff-surface-cache-ownership-proposal>
== Status
<status>
Implemented with one shared recursion engine in `three-dimensional-reps`.
`feynkit-cff` exposes that generalized engine and adapts it to FeynKit's
standalone topology and surface APIs. GammaLoop uses the generalized
engine for production expressions, including repeated poles and
numerator energy dependence.

== Shared generation
<generation-modes>
`CffGenerationGraph::try_surface_chains` owns source/sink selection,
contraction, cycle rejection, and surface-chain enumeration. Its surface
callback lets each consumer retain its own surface representation while
using the same recursion. Fixed boundary edges retain their direction
when contracted subgraphs are enumerated.

`feynkit_cff::generate_3d_expression` and `feynkit_cff::generalized` expose
the generalized engine directly. The topology-only `CffGenerator`
converts graph boundaries to that engine and reconstructs FeynKit trees
from the returned chains. It does not implement a second recursion.

== Standalone surface invariant
<invariant>
A FeynKit surface identifier is meaningful only with the `SurfaceCache`
supplied with its expression. Consumers discover the referenced surface
set from that expression's trees; a shared arena can contain additional
entries from other expressions.

`CffGenerator::generate` returns an expression, fresh arena, and report in
`CffResult`. `generate_into` accepts a caller-owned arena for related
standalone expressions. IDs remain append-only and stable. Raised-surface
analysis, residues, lowering, and diagnostics resolve only the IDs used
by the source expression.

Combining independent arenas requires the explicit `SurfaceIdMap`
returned by `SurfaceCache::merge` and applied by
`CffExpression::remap_surfaces`. Index coincidence is not an identity map.

== Runtime surfaces
<downstream-extensions>
GammaLoop supplies the generalized engine with its runtime E-surface and
H-surface types. These carry the numerical and threshold information
needed for production integrands and UV subtraction. Runtime caches and
FeynKit topology caches are separate representations; their identifiers
must never be interchanged implicitly.

The shared engine owns combinatorics, while each adapter owns conversion,
surface interning, and symbolic lowering. This retains main's generalized
production functionality while making it available through FeynKit.

== Required checks
<required-checks>
Tests cover stable standalone arena IDs, exact referenced surface sets,
explicit remapping, fixed boundaries, selected orientations, and
contracted/subgraph generation. Generalized and runtime regressions
cover repeated poles, signed contours, numerator energy maps, and UV
sources. Both API families must continue to use the shared recursion.
