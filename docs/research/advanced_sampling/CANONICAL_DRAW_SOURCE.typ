= A fixed sampling proposal across physical precision rescue
<a-fixed-sampling-proposal-across-physical-precision-rescue>
14 September 2026. Implemented; the acceptance gates below pass. This corrects the discrete-policy-only transport discussed in #link("AFFINE_STAR_IMPLEMENTATION_AUDIT.typ")[the preparation audit] and #link("TWO_NORMAL_PROPOSAL.typ")[the joint-normal derivation];. Their component geometry and normalization results remain useful, but do not prove invariance under a physical evaluation that selects different numerical maps.

== Why matching discrete decisions is insufficient
<why-matching-discrete-decisions-is-insufficient>
On a unit circle, let `T0(x) = x` and `T1(x) = x + delta (mod 1)`. Each is normalized, invertible, and has `J = q = 1`. Choose `T1` only when `x` lies in `[0, delta)`, and `T0` elsewhere. The resulting proposal misses `[0, delta)` and covers `[delta, 2 delta)` twice. A normalized target equal to `1/delta` on the missing interval integrates to zero under this switched rule. This holds for arbitrarily small positive delta, despite exact agreement of the two Jacobians and densities. The joint chart has a periodic circle coordinate, so continuous map changes are relevant to the production design.

Comparing complete compact/fallback decisions and dyadic radius exponents does not fix this issue. Neither does a small displacement or a same-point density comparison. Such comparisons can diagnose numerical accuracy; they cannot certify one proposal law under active-point-dependent physical rescue.

== One source, multiple physical approximations
<one-source-multiple-physical-approximations>
Prepare all production unit-cube maps once at the existing fixed 1000-bit Arb budget, before any physical target or stability-driven lane choice. Retain the complete graph-parent point and combined `J*w`, including every foreign inverse or proxy needed for the partition. A selected direct-momentum input retains its original point and the partition `w`, without a forward Jacobian. A direct input without a selected partition needs no sampling-map phase. Direct momentum input is already represented binary64 data and is embedded exactly in every physical lane, including values such as `0.1`. The former decimal-settings conversion at that boundary changed the point in higher precision. User-authored runtime settings retain their existing decimal policy.

Each physical Double, Quad or Arb attempt materializes directly from that immutable anchor. It must not redraw coordinates, repeat sampling root solves, recompute partitions, or reinterpret the selected LMB. Never promote a previous Double or Quad payload. Fixed-budget canonical preparation failure is a source error before the physical body; rescue cannot choose a different proposal.

This replaces the former prohibition on transporting numerical payloads with its actual purpose: prevent promotion of stale lower-precision geometry. Canonical point, factor and retained host data always remain the authority. Warmup must require a usable fixed-Arb binding for each participating graph. Failure of an optional native map binding cannot reject that canonical proposal solely because a physical lane uses a different numeric type. Native component and diagnostic APIs remain available with their own binding diagnostics. Lower precision is an approximation to physical evaluation of this fixed draw, subject to the existing stability controls. This is not a proof of accuracy for arbitrary discontinuous numerical targets.

== Existing ownership and factorization
<existing-ownership-and-factorization>
The existing Gamma/Discrete sample owners become prepared graph rows with their event-group boundaries. There is one channel catalogue and one map traversal; no additional numerical lookup table or unprepared summed route.

#figure(
  align(center)[#table(
    columns: 2,
    align: (auto,auto,),
    table.header([Source mode], [Preparation before physical evaluation],),
    table.hline(),
    [Selected channel], [Retain master forward map, partition and initial selected-host authority; foreign grouped graphs do not adopt the master\'s hosts.],
    [Explicit channel sum], [Move the existing per-graph bridge enumeration and `J_c*w_c` construction out of each physical attempt.],
    [Default or selected LMB], [Complete each graph\'s existing affine reinterpretation before freezing the row.],
    [Tropical], [Freeze the existing map and separate energy prefactor, preserving its current integrand-only application semantics.],
    [Direct selected channel], [Evaluate the partition once at the supplied raw point; apply `w`, without `J`.],
  )]
  , kind: table
  )

Apply the combined sampling factor exactly once to physical results and events. Preserve graph-group event concatenation, graph/orientation normalization and separate grid selection probabilities. Do not duplicate the same factor between the physical sample Jacobian and a second row field. Preserve tropical behavior without extending this change into an unrelated event-normalization correction.

The existing `EvaluationSource` borrows the prepared anchor outside the stability loop. Original input stays available for diagnostics. Norm selection and debug output inspect the anchor, without replaying maps or consuming root diagnostic occurrences. For an explicit sum, the norm check uses actual mapped rows rather than a dummy spherical point. The former no-op physical traversal used only to collect policies is removed. Distinct prepared rows receive distinct loop-cache identities. Reserve rotation identity ranges before evaluating bodies so a failed probe cannot leave used identities available to a subsequent draw.

Checked high-to-native conversion belongs to existing numeric, sample and ray owners. It rejects nonfinite values and nonzero underflow. Combined row factors also respect the allocated relative sampling-error budget. Exact source values remain available for physical alignment checks. Native host materialization is not a new Newton solve and must not invent diagnostic observations.

== What this does not settle
<what-this-does-not-settle>
Hosted H/Z binding still requires original equation constants, exact shared spatial routing, and physical adoption of the same host. In GL638, `H_global = H_cut + delta_host`; keep this residual explicit. Replacing `C1 = Q - Eh(a)` by an exact-on-cut identity at finite host residual silently moves the target. Physical CFF residue selection preserves surviving original surface equations; it does not universally replace them by cut-eliminated representatives. A generic hosted-joint check can bound the completed point\'s cut residual relative to its joint normal radius, alongside materialization of both original normal equations. It must include the physical rescaling\'s roundoff. This does not certify an arbitrary user multiplier: GL\'s rationalized WH needs its separate relation to the original equations established. Full conditional-density sensitivity and CT-star alignment remain separate from fixing the proposal law. Uniform relative accuracy over an entire disk reaching radius zero cannot be required at finite precision; pointwise unresolved cases must retry or error, not alter the sealed support.

== Acceptance and timing
<acceptance-and-timing>
Required gates include real original-source body retries preserving the anchor point and factor; independent Arb-to-Double and Arb-to-Quad materialization; selected, summed, direct and default/LMB routes; unchanged event factors; source failure before any physical evaluation; and the analytic circle control. Retain generated amplitude Gaussian/moment tests and physical conditional/raised fixtures. The generic hosted joint witness precedes a bounded GL638 replay, then final six-cut/all-936-orientation physical comparisons.

The correction passes 214 unique core tests, including actual source/body retries, direct materialization, range and factor boundaries, route matching, and conditional/raised physics. The complete all-18-orientation kite passes in 937.885 s with unchanged 8192-draw 6% normalization and 8% moment bounds, native physical rescue and exact forced-Quad agreement. Three API regressions pass; the reloaded bubble/kite/cut acceptance took 502.151 s with unchanged criteria. Its former five-minute timeout was insufficient and is now fifteen minutes. Routine fixtures and error/status assertions were migrated to the fixed-draw contract; no numerical tolerance was relaxed. Core/API checking, formatting and diff checks pass. All-target clippy reports 50 existing warnings and none on changed lines (56.86 s).

Charge canonical preparation, conversion, validation, inverse/proxy work and host adoption once to sampling time. Exclude that shared work from physical time. Removing native redraws may offset canonical cost; measure it rather than assuming an improvement. The production gate remains warmed optimized `T_sampling / T_physical <= 0.10` on representative and hard GL638 samples, using 20 cores and the actual all-orientation physical settings. Stop speed optimization once that gate is met. No current unoptimized component or amplitude-fixture duration establishes the GL638 timing result. Initially also report the conservative bound `(wrapper_total - T_physical) / T_physical`, which includes unclassified overhead such as row rotation and factor application. If this bound passes, further timing subdivisions are unnecessary. Time the already-warmed core process call: the batch API performs warmup per request. Exhausted errors fail the benchmark gate; they cannot be dropped for lacking successful-result timing metadata.
