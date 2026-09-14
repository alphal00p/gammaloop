# GL638 threshold-weight bias controls

Recorded 14 September 2026. The common positive weight-bias family passed 256
selected-point and local-limit evaluations at `sqrt(s)=600 GeV`,
`mu_r=91.188 GeV`, `m_uv=50 GeV`, with all 936 orientations summed, six cuts
and nineteen threshold variants. These controls establish identity at the
original prescription and preserve the measured limiting powers. They do not
rank Monte Carlo variance or establish a global weight bound.

## Prescription and cancellation conditions

On physical Cut 1 `(2,6,10)`, the existing threshold labels are
`A=(7,8)`, `U=(8,12,14)`, `F=(8,10,13,14)`, `V=(3,7,14)`,
`P=(3,12)` and `Z=(3,10,13)`. The parent LMB is `[3,6,7,10]`;
native two-loop variants use `[3,7]`, and shared one-loop variants use `[7]`.
All signed cycles, groups, centers and threshold multichanneling are retained.

With `Q=sqrt(s)`, the original native weights are

\[
W_H=\frac{H^2}{H^2+(PZ/Q)^2},\qquad
W_F=\frac{(G_0G_4)^2}{(G_0G_4)^2+(PZ)^2}.
\]

Here `H=eta(2,4,12)` is physical Cut 3 evaluated on the Cut-1 shell;
`G0=eta(2,6,12,13)` and `G4=eta(2,4,10,13)` are the corresponding Cut-0
and Cut-4 equations. H is not an additional normal-threshold CT of Cut 1.
The coordinate definitions and original mechanism are in
[the surface-sampling study](gl638.md) and
[the threshold-subtraction report](../gl638/README.md).

For finite `c>0`, apply

\[
T_c(w)=\frac{w}{w+c(1-w)},\qquad
1-T_c(w)=\frac{c(1-w)}{w+c(1-w)}.
\]

Native A/U receive `T_c(WH)`, native F receives `T_c(WF)`, and their shared
one-loop variants receive the respective exact complements. Native V/P/Z
remain at one. Exactly six multiplier expressions change. Each expression
still uses its own variant's `star`; local and integrated pieces retain the
same multiplier and the radial PV identity. The localization h-function is
unchanged. Complementarity holds at a common physical point, including a
regular pole, rather than requiring different off-shell star points to have
complementary numerical values.

| c | Effect |
|---:|---|
| 1 | Original prescription |
| 0.25 | Wider region weighted toward native two-loop subtraction |
| 4 | Wider region weighted toward shared one-loop subtraction |

The implementation combines the existing paired expressions `b=w`, `a=1-w`
as `b/(b+c*a)` and `c*a/(b+c*a)`, preserving the lazy common-zero convention
of WH. For `0<=w<=1`, the denominator is at least `min(1,c)`;
`T_c(w)=w/c+O(w^2)` near zero and
`1-T_c(w)=c*(1-w)+O((1-w)^2)` near one. Its derivative is at most four for
the tested interval. Endpoint orders and existing smooth common-weight
cancellations therefore survive, with changed finite coefficients.

The relevant sufficient conditions from the
[dual-weight audit](/common/dev/gl638_threshold_research/dual_weight_audit.md)
are preserved:

| Locus | Retained condition |
|---|---|
| Regular A/P and A/Z, away from H/Z contact | Complement vanishes as `P^2` or `Z^2`, matching the unit-weight partner |
| V/Z | Both native weights remain one |
| A/U/P, soft gluon 14 | A/U share one function, equal to one on P |
| A/V/P, soft gluon 14 | Additional mismatch remains `O(abs(q14)^2)` |
| U/F, soft gluon 13 | Native WH remains quadratic and WF quartic in `abs(q13)`; shared weights tend to one |
| F/Z, including H=0 | `1-WF=O(Z^2)` with the existing positive G0/G4 gaps on the relevant roots |

There is no general requirement that WH equal WF at every intersection.
A/F, P/F and P/Z have no relevant regular on-shell intersection here; U/F
requires the stated soft orders. The common transform is a conservative
one-parameter family. It changes the H/Z crossover width but does not remove
its inverse-radius body, CT image poles or coalescing radial roots.

## Executed controls

The original state and the regenerated states at `c=1,0.25,4` were evaluated
with Eager evaluators, native Fixed256 sampling, ordinary stability/rescue and
independently forced Arb evaluation:

- Eleven retained Samples × four prescriptions × two precision lanes:
  **88 evaluations**.
- Twelve H/Z points and nine single/double-soft points × four prescriptions
  × two lanes: **168 evaluations**.

All were finite, valid and stable. At c=1, real totals, complete Sample
estimators, retained cut values, CT weighted components and effective
multipliers matched the original within the declared scalar tolerances.
The largest relative discrepancies were `2.257e-12` against `1e-6` in the
ordinary lane and `5.662e-286` against `1e-10` in forced Arb. Real components
were checked independently of the complex norm.

Native points, forward/inverse Jacobians, complete partition scores/weights
and selected `J*w` were exactly identical for all prescriptions. Shared
parameter blocks were unchanged; the additional-parameter block alone grew
from `[200,200)` to `[200,201)`. Both saved-state trees remained unchanged.

The following are the last consecutive-interval exponents
`abs(quantity) proportional to radius^(-p)`. H/Z uses the complete real
sampling weight on two retained angles over six decades. Soft limits use the
raw real integrand at `lambda=1,1e-3,1e-6 GeV`, before the measure.

| Prescription | H/Z angle 0 p | H/Z angle 1 p | Gluon 13 p | Gluon 14 p | Both p |
|---|---:|---:|---:|---:|---:|
| Original | -1.000000 | approximately 0 | 0.999914 | 0.999979 | 0.999995 |
| c=1 | -1.000000 | approximately 0 | 0.999914 | 0.999979 | 0.999995 |
| c=0.25 | -1.000000 | approximately 0 | 0.999907 | 0.999980 | 0.999995 |
| c=4 | -1.000000 | approximately 0 | 0.999957 | 0.999977 | 0.999996 |

The H/Z weights vanish or plateau on these rays. The approximately
`lambda^-1` soft behavior is integrable with `lambda^2 d(lambda)` for one
soft gluon or `lambda^5 d(lambda)` for both scaled together. Actual radii and
binary64 input bit patterns are retained for the public raw-momentum controls.
These rays do not cover every correlated soft/threshold approach or tangency.

The retained right-Cut-3 control weight, `3.24161e-5`, is effectively unchanged.
The deepest angle-1 H/Z weight changes from `-3.67047e-12` at c=1 to
`-3.80438e-12` at c=0.25 and `-3.21705e-12` at c=4. This is no variance
ranking. The signed integral is invariant under the subtraction change by the
radial identity; its numerical compatibility still requires MC comparison.
The absolute integral `integral(abs(Re f))` can legitimately change with the
subtraction prescription.

## Reproduction and evidence

All scratch paths below are relative to
`/tmp/gl638-final-physics-screen/weight_variants/`. The actual generated DOT
declares the existing generic parameter interface:

```dot
params = "threshold_weight_bias";
```

After generating that DOT, a runtime overlay selects c in declaration order:

```toml
[general]
additional_param_values = [0.25]
```

This requires one isolated regeneration; an additional runtime value cannot
replace multiplier expressions in the original saved state. The successful
generation used build9 (`generation-build9-run.json`, exit 0, 346.50 s),
followed by build10 controls. `generation-card.toml` preserves the original
1000 GeV generation catalogue, all orientations, process-selected physical
cuts, factorized numerator, full direct-3D local UV and integrated UV, and
Eager generation. The complete `full-runtime-{identity,prefer_2l,prefer_1l}.toml`
cards select the 600 GeV nine-channel setup. `metadata.toml` and
`expressions.json` retain the exact changes.

`manifest.json` records inputs and 3003 exact-rational algebraic checks;
its preparation-status string predates generation. Completed acceptance is
recorded separately in `RESULTS_BUILD10.md`, `controls-audit-build10.json`
and `local-audit-build10.json`. Raw captures and before/after state hashes are
in `controls-{original,parameterized}-build10/summary.json` and the matching
`local-*` directories. Requests, run records and
`weight-controls-build10.build.json` bind inputs and the control executable.
The optional c-specific `native-boundary-requests/` were prepared, not executed
as part of these 256 controls; see the distinct
[native-precision acceptance](GL638_NATIVE_PRECISION_ACCEPTANCE.md).

| Artifact | SHA-256 |
|---|---|
| GL638_weight_bias.dot | `4e1ed29cec1144d9678eae4055081b41606f744e706f2f505268616d6959e820` |
| generation-card.toml | `11ff7d7db0d3268e92057a7ba883a32db18a99b516d153afefcfb44904f64177` |
| controls-audit-build10.json | `0e47bfc7df8d3421102e7c12c9d2616f4faddfd231b77e3f21b2c6e3f7b04b20` |
| local-audit-build10.json | `edb0c7e01e8d07dc933a3b704ebe8a09d93ee399dba34eeff8f8a02c78193375` |

The preparation reasoning and predeclared gates remain in
`WEIGHT_VARIANT_REVIEW.md`. No long-run estimate, performance comparison or
Monte Carlo variance conclusion follows from this control study.
