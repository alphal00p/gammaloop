# Native sampling precision and rescue

Status: native components and geometry pass 122 isolated numerical tests;
production binding and original-source stability rescue remain to implement.
This note complements [ADVANCED_SAMPLING_PLAN.md](../../../ADVANCED_SAMPLING_PLAN.md)
and [LU_H_MATCHED_SAMPLING.md](LU_H_MATCHED_SAMPLING.md).

## First divergence and required invariant

For an X-space evaluation both routes start with the same original cube sample,
canonical `SamplingChannelId`, graph and orientation. Ordinary `F<T>` mapping
preserves T. The compiled route in `gammaloop_sample.rs::parameterize` currently
converts coordinates and improved externals to f64, calls the f64 bridge, then
`SamplingChannelBridgeEvaluation::to_momentum_sample<T>` promotes its answer.
Increasing the physical evaluator precision cannot recover information lost
there. Summed channels and direct-momentum inverse partitions have the same
boundary. `evaluate_from_source` also propagates a map error through `?` before
trying another stability level.

The invariant is one mathematical proposal evaluated throughout in the active
precision: source -> native kinematics -> selected map -> every foreign inverse
and score -> partition -> physical/reference value -> final weighted result.
Only original f64 settings/cube inputs and final public f64 reporting are valid
conversion boundaries. Recomputing after promotion of a rounded point, root,
rotation, improved external vector or `t*` does not implement rescue.

Preserve the exact binary value of a generated Monte Carlo cube coordinate when
reconstructing it at higher precision. The general `from_f64` conversion uses
the shortest decimal spelling for physical input conventions; that is a
different source point. Reuse the existing exact-binary constructors at the
retained-draw boundary, while keeping the established decimal interpretation of
user model/kinematic settings. Test a value such as binary64 `0.1` whose exact
binary embedding differs from decimal `0.1`.

GammaLoop's `QuadFloat` wraps Symbolica's two-binary64 `DoubleFloat`: it adds
mantissa precision but retains binary64's exponent range. Native tests must
distinguish these two capabilities. Sub-f64 coordinate differences can recover
in Quad; range failures such as `exp(-1000)` or `10^400` require Arb. A Quad
retry that remains unrepresentable is an expected precision outcome, not proof
of a broken evaluator. Test that it continues to Arb without dropping support.

## One owner and object-safe native components

Parameterize the existing trait, rather than adding generic methods to a trait
object: `SamplingMapComponent<T: FloatLike = f64>`. Its `forward`/`inverse`
take `&[T]` and return `SamplingMapEvaluation<T>` with native scalar fields.
Trait-level T is object-safe; `Box<dyn SamplingMapComponent<T>>` retains existing
conditional compositions. Parameterize the compiled map, affine/embedding/
composition owners, bridge, contexts and partition consistently, defaulting to
f64 only for existing caller ergonomics. Math uses GammaLoop's `F<T>` wrappers:
`F(value)` and `.0` preserve every native digit and active Arb precision; they
are not numerical conversions. Keep existing kernel/radial algorithms, remove
actual f64 conversions, and call the same implementation for every T. Defaults
permit native bridge gates before the host migration without a parallel engine.

Keep the existing compiled-map enum as its current dispatch owner. Replacing
all dynamic callbacks with a new closed enum is unnecessary for this migration;
three methods named `forward_f64`/`forward_quad`/`forward_arb` would duplicate the
precision dispatch. Do not introduce a second map interpreter or channel axis.

Resolve metadata, canonical IDs, routing and Symbolica expressions once. Extend
`LmbMultiChannelingSetup.sampling_bridge: RuntimeCache<_>` with lazy native
numeric instances under that same resolved catalogue. An instance contains
native bindings/kernels indexed by the immutable catalogue, never a reselected
or compacted list. The Double/Quad/Arb branches already present in the stability
owner select the corresponding typed cache; follow the existing native external
cache access pattern if generic lookup is needed, without `Any` downcasts.

`warm_up_sampling` validates/builds the immutable plan transactionally. Reuse
compiled expression programs and routing across native instances and workers;
worker-local mutable evaluator storage follows the current clone policy.
Invalidate bindings on model, external/improvement, routing, selection and
profile changes. Precision and native external-cache identity belong to binding
identity. A rotated or rescued point must never reuse a stale f64 binding.

## Native equations, contexts and eager evaluation

Make `ImplicitSurfaceRadialMap<T>` callbacks, centers and context vectors native.
Generalize `Esurface::sampling_radial_map` from the actual native graph/model
mass values and native external cache. Reuse `compute_self_and_r_derivative`,
including the existing massless-origin right derivative. Exact routing integers
stay exact; affine shifts, inverse matrices and determinants are computed in T.
The existing generic graph mass accessor avoids promoting a derived f64 mass.

Generalize `PreparedCutSamplingContext`, `PreparedSurfaceStatus`, prepared
surfaces and `DeferredCrossSectionSamplingState`: their current f64 `t*`, roots,
loop/external momenta cannot be inputs to native rescue. Retain graph/cut/side,
complete parent frame, branch and root identity with native values. Build each
cut's LU rescaling and threshold external data anew at T before conditional
left/right maps. `cast_sample` is not rescue for derived contexts/Jacobians.
Rebuild or rotate typed native contexts consistently with the mapped point;
copying opaque vectors cannot rotate their hidden momentum/basis entries.
Serialized diagnostic views may convert after computation.

Generalize `SamplingExpressionEvaluator` values, inputs, dual seeds, real-part
checks and determinant calculation. `FloatLike` already includes
`GenericEvaluatorFloat`; `<T as GenericEvaluatorFloat>::get_evaluator` selects
existing native eager evaluators. Parse/optimize once and reuse those evaluators;
instantiate map bindings lazily and preserve exact constants. Use native imaginary
roundoff tests, not fixed `1e-13`. Stable hyperbolic radial primitives remain
native until equivalent eager primitives and dual derivatives are registered.

Generalize callback/proxy scores and `SamplingPartition::from_log_scores` to
native log-sum-exp. Evaluate every foreign inverse at the SAME native raw point,
in its own native context. Keep support absence distinct from evaluation errors;
no failed inverse may silently acquire zero score. Preserve positive denominator
coverage and test extreme native log-score separations. The final selected
`J_c*w_c` product stays native even when either f64 factor is unrepresentable.

## Root localization and numerical errors

Let h be half an ulp of a positive radius R, and s=R/(R+beta). Rounding to R
occurs approximately inside `|u-s| < s*(h/R)^(1/p)` on the inner branch and
`|u-s| < (1-s)*z/(1+z)`, z=(h/beta)^(1/p), on the outer branch. At R=3,
beta=2 the combined f64 width is about 9.38e-9 for p=2 and 4.44e-6 for p=3.
Translations/angular reconstruction can add cancellation. This is an arithmetic
limit, not permission to clip a band.

Root tolerance must also resolve the intended signed distance delta=r-R. A
fixed `1e-11` energy residual can misalign the density by much more than delta.
Reuse `utils/newton_solver.rs` native safeguarded solver, `RadialRootIdentity`
and precision-carrying `RadialRootDiagnostics`; migrate the sampling solver
rather than copying its iteration loop. Add a sampling-local accuracy target:
require a root uncertainty bound below, for example, `|delta|/16`, in addition
to positive slope, native residual and bracket checks. That bound preserves
localization, but density accuracy needs the stronger budget
`alpha*root_uncertainty/|delta| << requested_density_error`, alpha=1-1/p,
plus independent inverse/Jacobian checks. Refine or escalate if
that target is below attainable precision; do not relax it to the previous
precision's absolute residual target merely because the residual improved.

A Newton correction alone is not a certified uncertainty for arbitrary
callbacks. Use reliable signed brackets, or a residual error bound together
with a positive derivative lower bound on the bracket. Physical energy sums can
supply their roundoff/conditioning bound; its derivation must be checked. The
existing Newton diagnostics explicitly do not promise a forward-error bound
for arbitrary callbacks. Keep that limitation explicit until supplied evidence
justifies the stronger claim. A stored root bracket must remain valid after
any angular/frame reconstruction; reconstruct/refine using the same equation.

Preserve selected root/frame identity through forward, inverse and rescue;
foreign channels obtain independent evidence for their own roots. A computed
direction-dependent root still gives the radial determinant (angular root
terms cancel in the wedge product), but exact proposal bookkeeping alone does
not establish localization on the true physical surface. The mathematical root
must depend only on direction/context: u may request tighter numerical accuracy,
not redefine the chart through an uncontrolled u-dependent stopping error.

Introduce one typed numerical sampling-error classification at the existing map
boundary, carrying channel, operation, precision and root/coordinate evidence.
Rounded seams/tails and numerically unresolved interiors/roots are retryable.
Invalid metadata/dimensions, unsupported geometry, structural missing interior
centers and demonstrably invalid equations remain hard errors. Uncertain signs
near the numerical floor must not masquerade as proven absent surfaces.

Catch only the typed retryable errors in `evaluate_from_source`, record the
failed level, and rebuild from its ORIGINAL `EvaluationSource` at the next
configured precision. Recompute selected mapping, all inverses, cut contexts
and the final weight; discard partial events/results from the failed attempt.
Apply this to MC, summed channels, direct-momentum scoring and precise API
reevaluation. A deterministic root identity includes graph, canonical channel,
map block/cut and frame; occurrence bookkeeping must survive early exits.
After final precision exhaustion return a diagnostic error, never a zero sample.

## Staged owners and acceptance gates

1. `sampling_maps.rs`, `sampling_evaluator.rs`: generic component/results,
   affine/composition/implicit callbacks, native eager duals and constants.
   Test f64/Quad/Arb roundtrips, independent 3D/6D determinants, tails and seams.
2. `sampling_selection.rs`, `sampling_partition.rs`, `sampling_context.rs`:
   native bridge/foreign densities/contexts and one-catalogue lazy bindings.
   Test mixed channels, conditional contexts, extreme log scores, unchanged IDs,
   warmup reuse/invalidation and worker-local evaluator independence.
3. `cff/esurface.rs`, graph warmup owners, `gammaloop_sample.rs`, `process/mod.rs`:
   native physical bindings, original-source reconstruction and typed retry.
   Extend the existing Newton owner with focused-distance evidence; never alter
   LU acceptance semantics globally merely to accommodate sampling.
4. API state acceptance: a chosen cube point fails the f64 seam check but passes
   Quad/Arb and agrees with direct native evaluation, including all foreign
   scores, events and second moments. Test summed=explicit-MC, rotations,
   saved/reloaded states and conditional cut contexts. Structural errors must
   not retry; an unresolvable exact seam must fail explicitly at the final level.
5. Use a curved physical surface whose fixed-tolerance root error exceeds delta:
   confirm native refinement restores actual signed-distance scaling and inverse
   reciprocity, rather than merely integrating a displaced spherical proposal.
   Repeat generated kite/cut and Gaussian normalization/moment gates, then a
   long-run error/weight audit before making bounded-weight or GL638 claims.

Parallel ownership after the cache milestone: kernel/bridge agent owns
`sampling_maps.rs` and `sampling_selection.rs`; expression/partition agent owns
`sampling_evaluator.rs`, `sampling_partition.rs` and reference/API acceptance;
physics/stability agent owns `sampling_context.rs`, `cff/esurface.rs`, graph
warmup/cache factories, `gammaloop_sample.rs`, stability dispatch and the shared
Newton accuracy handoff. Agree the generic signatures, cache lookup and typed
error first; complete the first two interfaces before activating native physics.
Preserve borrowed per-point map/proxy evaluation; only worker ownership clones
stateful eager storage, never the partition loop.
