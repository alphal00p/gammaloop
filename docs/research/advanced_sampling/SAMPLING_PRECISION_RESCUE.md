# Native sampling precision and rescue

Status: native components and geometry pass 122 isolated numerical tests.
Production native binding and original-source stability rescue now pass eight
focused tests, including actual kite/cut reconstruction, event factors and
precise API range. The broader 169-test core suite, both saved-state/summed API
regressions and six precise/event/histogram integration API tests also pass.
Formatting, core/API and integration-test checking pass; clippy reports no
warnings on changed lines.
Root-distance
density certification, native reference-harness retry, derived-expression mass
precision and the final adaptive-grid range boundary remain open.
This note complements [ADVANCED_SAMPLING_PLAN.md](../../../ADVANCED_SAMPLING_PLAN.md)
and [LU_H_MATCHED_SAMPLING.md](LU_H_MATCHED_SAMPLING.md).

## First divergence and required invariant

The first divergence before the native host migration was the compiled route in
`gammaloop_sample.rs::parameterize`: it converted coordinates and improved
externals to f64, called the f64 bridge, then promoted the answer. Ordinary
`F<T>` mapping preserved T. Both started from the same cube sample, canonical
channel, graph and orientation; increasing only physical evaluator precision
could not recover information lost at the map boundary. Summed channels and
direct-momentum inverse partitions shared that loss, while map errors escaped
before a later stability level. These production paths now use native bindings
and typed retry from the retained source in one shared precise-result loop.

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
Direct model masses retain the repository's original f64 input convention.
The 13 September owner audit found that `EdgeMass::value<T>` also evaluates
`EdgeMass::Evaluator` expressions in f64 before converting to T; merely calling
the generic graph accessor does not recover those derived digits. The current
sampling binding agrees with the physical evaluator's mass cache. Generalizing
that existing expression evaluator is a remaining precision requirement for
derived-expression masses, separate from native kinematics and map rescue.

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

The selected-channel route must also combine `J_c*w_c*f` before narrowing the
physical/reference result. For example, a native `J_c=10^400` and
`w_c*f=10^-400` have a finite product but become `Inf*0` after separate f64
conversion. Reuse `GraphEvaluationResult::apply_sampling_factor` for values,
events, full multiplicative factors and reference moments; report a unit
unapplied outer Jacobian, as summed channels already do. Preserve the actual
map Jacobian from the bridge map evaluation in replay diagnostics, separately
from its partition weight; the selected momentum sample carries their product.
Apply ordinary
X-space Jacobians at this same boundary; direct momentum input retains unit J.

Outer adaptive-grid weighting is a distinct remaining range boundary: a native
map/partition/physics product `10^-400` times a representable grid weight
`10^100` should contribute `10^-300`. Narrowing before the integration
accumulator applies that weight loses the contribution. Audit this through the
existing estimator/training owner in the trained-grid acceptance work; native
map rescue alone does not establish range-safe final grid weighting. Likewise,
Rust precise-result APIs must not require representability of a preliminary
f64 reporting pass before returning their native result.

## Root localization and numerical errors

Let h be half an ulp of a positive radius R, and s=R/(R+beta). Rounding to R
occurs approximately inside `|u-s| < s*(h/R)^(1/p)` on the inner branch and
`|u-s| < (1-s)*z/(1+z)`, z=(h/beta)^(1/p), on the outer branch. At R=3,
beta=2 the combined f64 width is about 9.38e-9 for p=2 and 4.44e-6 for p=3.
Translations/angular reconstruction can add cancellation. This is an arithmetic
limit, not permission to clip a band.

For faithful focusing on the true physical surface, root tolerance must resolve
the intended signed distance delta=r-R. A fixed `1e-11` energy residual can
misalign the density by much more than delta. This is a stronger requirement
than consistency of the actual numerical proposal, distinguished below.
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

### Minimum practical proposal-consistency gate

A deterministic numerical root `Rhat(n, context, T)` can itself define a valid
proposal. Its angular derivatives cancel from the radial determinant whether
or not it is the exact physical root. A small physical-root displacement then
changes efficiency and asymptotic alignment, rather than inherently biasing the
proposal. Keep this root independent of the sampled radial coordinate `u`;
adapting its numerical stopping rule to `u` without corresponding inverse
semantics would undermine that simple map definition. Raising the precision
rebuilds the whole proposal from the original draw as already implemented.

The 13 September audit found an inverse reconstructing coordinates, calling
forward on them and returning that nearby point's inverse determinant, with
only a finite coordinate-residual check. Small coordinate errors near a shell
can still change the density substantially. The bounded correction is:

1. Extend the existing radial inverse owner to return its derivative at the
   **supplied radius**, using its fixed recovered direction/root/context. Reuse
   this owner in analytic and implicit charts. Forward reconstruction becomes a
   checked diagnostic, not the source of an unchecked nearby-point density.
2. In the canonical bridge forward path compare `J_forward*q_selected` with one
   at the actual raw point. Map-density partitions already compute the selected
   inverse log score; reuse it. In proxy mode check only the selected channel's
   exact inverse, without requiring unrelated foreign inverses or replacing
   the supplied proxy scores. A failed finite/positive or relative-accuracy check
   raises the existing typed numerical error and retries the original sample.
3. Thread a relative-density budget from the configured stability requirement,
   retaining the stricter budget when the same precision occurs more than once.
   A native `sqrt(epsilon)` default is suitable for fresh standalone constructors;
   it is not a replacement for stricter user requests. Compare in log space to
   avoid overflow of individually representable Jacobian/density factors.
   The implementation uses `log(1 +/- tau)=+/-2*atanh(tau/(2 +/- tau))`
   for small budgets. This preserves the budget, but native logarithm rounding
   still limits what the numerical comparison establishes; it is not an
   interval proof of arbitrarily sub-epsilon accuracy.

For the signed-distance profile, let `alpha=1-1/p`, `beta>0`, and fix the actual
inverse radius `r` and root `Rhat`. Its positive radial density is

```
r < Rhat: q_r = ((Rhat-r)/Rhat)^(-alpha) / [p (Rhat+beta)],
r > Rhat: z = ((r-Rhat)/beta)^(1/p),
          q_r = ((r-Rhat)/beta)^(-alpha) / [p (Rhat+beta) (1+z)^2].
```

The same outer formula at `Rhat=0` preserves the existing rootless power law.
For `p=1` its ordinary limit is `beta/(beta+r)^2`. Evaluate from the supplied
signed distance, avoiding a second subtraction of recovered `u` and its split.
The Cartesian density includes `1/[J_Omega*r^(D-1)]` in the same angular chart.
A true pole/angle seam keeps its explicit measure-zero convention; a finite band
of failed arithmetic may never be clipped away or assigned zero support.

Decisive tests are: a reversible wrong-inverse fixture rejected in both density
and proxy modes; a proxy fixture whose foreign inverse is unavailable but whose
selected inverse is correct; native reconstruction where Double fails and Quad
or Arb passes; and the generated kite C surface plus the physical bubble cut at
representable distances that are small enough to expose normal-density errors.
For the massive rest kite with k=-l, the C equation is
`2*sqrt(1+|k|^2)+1-5=0`, giving a six-dimensional radius `sqrt(6)` when the shared
spatial direction is unit-normalized. Use the actual routed equation, and
account for the represented direction norm in a high-precision oracle. This
provides an independent localization check rather than another copy of the map.
Repeat direct Jacobians, all foreign densities, rotations, and shifted Gaussian
normalization/moments. After these bounded gates, the power-2 amplitude matrix
can measure efficiency while reporting near-surface/high-precision comparisons;
it need not wait for a formal global bounded-weight theorem.

The generated rest-kite gate supplied a concrete numerical check on 13 September:
with angular cube coordinates `[0.27,0.61,0.39,0.72,0.58]` and `u=s+1e-5`,
a fresh f64 C bridge measured `log(J*q)=1.553046278e-7`, exceeding its native
`sqrt(epsilon)=1.490116119e-8` budget although the mapped radius was finite.
The fixture now requires that typed failure and checks the same original draw
in Quad. Production bridges use their configured stability-derived budget;
this fresh-constructor result is not a claim that every production budget fails.

A stronger **rigorous physical-root certificate** still requires reliable
function enclosures. The smallest addition can use the current safeguarded
solver for a candidate and prove opposite signs at `Rhat +/- e` using the
existing physical E-surface owner. No second root solver is needed. Floating
energy evaluation must supply a justified roundoff enclosure; an arbitrary
multiple of epsilon, a Newton correction, or agreement between two precisions
is not by itself a proof. For a certified root interval, bound the radial
density over that interval at the actual raw point; an interval crossing r
cannot certify the shell side. Arbitrary callbacks lacking an enclosure retain
that explicit limitation. This stronger evidence is needed for rigorous
true-surface asymptotic claims, not to label the preceding numerical measurements
as measurements of the implemented proposal.

### Smooth LU-h profile cross-check

For X1's full-parent zero-centered cut map, `t=R/r` and
`q_r=p(t|R)*R/r^2`; equivalently `|dr/du|=r/G_y` for `y=log(t)` and
`G_y=partial_y G`. The broad scale `a=R/beta` induces the normalized raw floor
`epsilon*beta/(beta+r)^2`, independent of the root. Positive epsilon therefore
protects the smooth Gaussian reference's radial variance at the origin and
infinity; it does not bound physical angular, soft or threshold weights.

With the documented `u=G(t)` convention, broad-only sampling is
`r=beta*(1-u)/u`. Recognize coincident component quantiles before invoking a
solver requiring strict endpoint signs. Upper-tail inversion should solve
`(1-u)-S(t)=0`, whose derivative is positive; solve in log(t) and avoid evaluating
overflowing inactive formulas. These choices preserve one monotone map and the
full mixture density, with no latent branch/channel axis.

For log-logistic shape kappa, its log-density derivative in y lies between
`-kappa-1` and `kappa-1`; the broad component lies between -2 and 0. The mixture
log derivative is their positive weighted average. Consequently
`|partial_y log J_K| <= D+max(kappa,1)`. A log-quantile enclosure of width e_y
bounds determinant variation by `exp((D+max(kappa,1))*e_y)-1`. A scaled CDF
residual alone is not a uniform quantile/density certificate for arbitrarily
small user-supplied shape. The direct inverse-density and bridge consistency
checks apply to this smooth profile too, without inheriting the threshold-seam
veto. Existing absent/pinched classification remains separate from failed
numerical root evaluation.

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
