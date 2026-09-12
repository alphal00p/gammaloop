# GL638 runtime stability and sampling follow-up

The current generic Euler stability probe rescues all four tested nearby A/P
cancellation controls. The previously reported nearby false acceptance used
the historical pi/2 z probe. Exact root coincidence still needs a combined
evaluation. Separately, the hard H/Z corner retains its measured inverse-radius
growth. An additional existing LMB improves a short integration comparison;
it does not change that local power.

These experiments used commit `d957817a7`, Symbolica
`a277a8dadbaae5a8d6605f2f381991a601e1f345`, and release SHA256
`029d4ecee30c4e038889ef62725ccfcdbf6e88f5d9e7dbfea08b627d238712a3`.
The freshly generated state includes all 936 orientations, the WH/WF metadata,
direct three-dimensional local UV with orientation localization, integrated UV,
and the complete six-cut sum. Q=1000 GeV, mt=173 GeV and mh=125 GeV.
Every state load was read-only and its files remained unchanged. The sampling
pilots used 20 pinned CPUs; precision replays used one. No failed samples were
removed from the reported statistics.

## Precision: nearby cancellation versus exact coincidence

The 48-case runtime matrix compares generic Euler, the historical z probe,
componentwise checks, extra rotations, stricter Quad tolerances, direct Arb,
and escalation controls. A six-case follow-up checks component accuracy at
three easier inputs. All points use raw generation-LMB momenta, without a
sampling Jacobian or adaptive density. The offsets below are threshold energies,
not gluon soft radii or measurements of the native radial-root gap.

| ηP [GeV], ηA/ηP=1.0001 | Historical z Im relative error | Euler norm: accepted precision / Im error | Euler components: precision / Im error |
|---:|---:|---|---|
| 1 | 29.3 | Quad / 2.48e−15 | Quad / 6.53e−16 |
| 0.01 | 2.06e7 | Quad / 5.93e−10 | Arb / 0 exported |
| 0.0001 | 2.08e13 | Arb / 0 exported | Arb / 0 exported |
| 0.000001 | 1119 | Arb / 0 exported | Arb / 0 exported |

Errors are `abs(normal−reference)/abs(reference)`. “0 exported” means equality
of the recorded values, including the f64 exports of Arb results. References
use matching norm/component return conventions. At the highlighted 1e−6-GeV
point, the Euler Quad discrepancy is 0.932519 and correctly triggers Arb.
Every primary cut then agrees with direct Arb. Component checking alone with
the historical z rotation still accepts the wrong value: that exact coordinate
permutation is a poor independent roundoff probe.

For routine integration, the current card's Euler(0.1,0.2,0.3) ladder is a
useful mitigation. Its Double/Quad/Arb relative tolerances are 1e−6/1e−10/1e−12.
For a smaller component that matters, the existing
`-D stability_check_on_norm=false` is worth testing. At ηP=0.01 it promotes
Quad to Arb and restores the reference, increasing evaluation time from
0.219 to 2.630 seconds. It retains Quad at ηP=1 (0.227 versus 0.261 seconds).
This is a targeted cost/accuracy observation, not a global prescription.

The implementation compares magnitudes when `check_on_norm=true`, using the
Re tolerance, and returns the unrotated value. With `false`, it tests the two
components separately and returns their rotation average; retained events
still describe the primary evaluation. A Stable label is a heuristic,
especially for a small component. Three generic Euler rotations also recover
the difficult point, but cost about 5.71 seconds versus 2.66 for one. Tightening
Quad to 1e−20 forces recovery even with z at the highlighted nearby point,
while unnecessarily promoting the two easier Euler controls. It cannot rescue
an incorrect value already accepted at Double. Skipping Quad likewise loses those inexpensive
rescues. These results do not motivate imposing those costs globally.

Arb is compiled as `VarFloat<1000>`: approximately 301 decimal digits. The
card's `arb_precision` placeholder sets the acceptance tolerance, not the bit
count. Direct Arb variants agree on the nearby point; their reported internal
rotation discrepancies are about 1e−282 to 1e−266. There is no evidence here
that additional arithmetic precision is needed for these distinct-root inputs.

At exact coincidence, every one of the sixteen configurations reaching Arb
still reports `is_nan=true`, final Unstable status and an invalid cut-1 value.
The two Quad-only controls return finite but Unstable values. The relevant
algebraic model is

```text
g(a)/(a-b) + g(b)/(b-a) = [g(a)-g(b)]/(a-b) → g′(a).
```

Arb can recover the nearby cancellation; it cannot make the separate divisions
defined at a=b. A combined divided-difference/confluent evaluator must preserve
the full smooth kernel, including projection Jacobians, damping, multipliers
and group factors. Denominator co-occurrence within an amplitude matters;
membership of a shared SOCP group alone does not establish the residue identity.
This is not a request for same-amplitude iterated threshold counterterms. The
exact equality is measure zero and does not invalidate the useful nearby rescue.
Nor do these controls settle every GL638 stability flag: center construction
includes f64 geometry, and failed rotation checks need point-specific diagnosis.

The large-weight escalation setting was not exercised by raw `inspect`, which
supplies zero as the historical maximum. Loop-norm escalation tests total input
momentum size, not a threshold-root gap. In component mode, the present
large-weight predicate escalates only when both phase thresholds are exceeded;
it should not be advertised as an either-component safeguard.

## Existing sampling: three matched pilots

Each pilot completed three 10,000-sample iterations with seed 1337, linear
mapping, b=0.3, inverse-Jacobian LMB partition weights, real-component training,
and the same Euler norm ladder. Saved runtime snapshots agree outside sampling.
Basis IDs refer to the actual generated catalogue, not local grid positions.

| Proposal | Signed error Re / Im [pb] | Largest absolute weight Re / Im | Absolute-integral maximum impact Re / Im | Seconds | Final Unstable / NaN |
|---|---:|---:|---:|---:|---:|
| Spherical, six LMBs | 8.18e−5 / 5.20e−4 | 1.785 / 15.136 | 9.60% / 32.69% | 370.58 | 0 / 0 |
| Spherical, six plus LMB84 | 7.17e−5 / 1.34e−4 | 0.873 / 2.535 | 3.68% / 6.69% | 367.63 | 0 / 0 |
| Common radial, six LMBs | 9.49e−5 / 1.73e−4 | 1.773 / 3.227 | 9.77% / 9.97% | 389.61 | 3 / 0 |

Maximum impact is `abs(wmax)/(N*absolute_integral_estimate)`. It is descriptive,
not the exact contribution under the combination of adaptive iterations.
Final Double/Quad/Arb counts were 29659/341/0, 29531/469/0 and 29172/825/3.
The three common-radial failures were retained; their input points were not
recorded, so the aggregate flags cannot identify their cause.

The five distinct saved Re/Im maxima were replayed in normal and forced Arb
precision (the common-radial Re/Im maxima share the same point). All ten
evaluations were finite and Stable, with normal runs accepting Double.
The largest normal/Arb differences were 7.87e−13 in the total complex value,
1.13e−12 in an individual total phase, and 1.47e−12 in any cut. The original
adaptive PDFs are absent from these inspections, so their returned values are
not reconstructions of the saved Monte Carlo weights. These replays validate
precision at the extrema, without identifying each extremum as an H/Z point.

The six original IDs are `[1,2,13,14,21,22]`, with edge sets
`[6,12,13,14]`, `[4,12,13,14]`, `[6,7,13,14]`, `[4,7,13,14]`,
`[6,10,13,14]`, `[4,10,13,14]`. Basis84 is `[2,10,13,14]`, including
the Higgs and both gluon edges. It is already in the saved state's 127-basis
catalogue. The existing runtime setting is

```toml
[sampling]
lmb_basis_ids = { GL638 = [1, 2, 13, 14, 21, 22, 84] }
```

The seven-channel error-squared times wall-time ratios are 0.760 (Re) and
0.066 (Im), relative to the six-channel run. The baseline Im error was driven
by a second-iteration outlier. These are single-seed finite-budget statistics:
changing channels or coordinates changes physical samples and learned grids,
so this is not a paired-sample experiment or a measured asymptotic variance
ratio. The additional basis is a promising optional choice; longer independent
seeds are needed before treating it as a new default. The common-radius map
does not improve both components in this comparison.

## Density needed at the remaining corner

The current-binary replay evaluated eight H/Z ray points with both forced-Arb
rotation choices, plus four normal-Euler controls. All twenty were finite and
Stable. The sixteen Arb totals and 96 per-cut values are bit-identical to the
earlier hard-corner scan. Both components retain approximately `f=C/R`, with
last-decade fitted powers between −1.00027 and −0.99999. The normal/Arb
complex-total discrepancy is at most 2.45e−10 on the four controls.

Here H=η(2,4,12), P=η(3,12), Z=η(3,10,13), and locally
`u=H`, `v=abs(P0)*Z/Q`, `R=sqrt(u²+v²)`, with P0 nonzero at the corner.
R is an energy distance in the two normal directions, not a gluon soft radius
or the twelve-dimensional loop hyperradius. The sampled range is 0.2 to
0.0002 GeV. All eligible LMB edge radii remain nonzero there.

If these are regular independent normals with a bounded coefficient that is
nonzero on an angular/tangential neighborhood, the local measure is a smooth
factor times `R dR dθ`. Consequently `abs(f)` has a finite local integral, while a
smooth nonzero proposal leaves the second moment proportional to `dR/R`.
Two rays support the power; they do not prove its angular or tangential bounds.

A direct solution for that transverse degree is to sample R uniformly on
`(0,Rmax)` and θ uniformly on `(0,2π)`. In the normal plane this gives

```text
q(u,v) = 1/(2π Rmax R).
```

It cancels the leading 1/R weight growth. The physical momentum density also
needs the inverse Jacobian of the normal/tangent coordinate chart. Add this
as a normalized channel mixed with ordinary full-support sampling. Every LU
rescaling, threshold projection, center dependence, inverse branch and support
used to define that channel must enter its actual density. A channel targeting
one cut or one star is legitimate.

A two-normal map is not mathematically necessary merely for finite variance.
Sampling one correctly aligned normal δ=H with

```text
qδ(δ) = (1-β)/(2Δ) * (abs(δ)/Δ)^(-β),  |δ|<Δ,  0<β<1
δ = random_sign * Δ * U^(1/(1-β))
```

already gives a local second moment proportional to
`∫ R^(β−1) dR ∫ abs(cos θ)^β dθ`, which is finite under the same regularity
assumptions. β=1/2 is a simple trial. Its typical leading weights still grow
as `R^(-1/2)`; the two-normal `1/R` density is better targeted at bounding
maximum weights from this corner.

Alignment must use the actual offending geometry. Earlier integrated-A maxima
had H close to zero at their A star while sampled H was not close. A channel
targeting only base H can miss that counterterm's pulled-back singularity.
An approximate star or fixed center may help a finite window but does not
establish the asymptotic result. Piecewise overlap/center choices also require
consistent branch accounting. Existing root/projection and multichannel
infrastructure can be reused; a threshold-centered singular density and its
inverse probability are additional functionality, not a current TOML knob.
No integrand or localization multiplier is altered by importance sampling, so
the required local-CT symmetry about r-star is preserved.

Existing LMB maps, radial scales, channel probabilities and adaptive bins can
improve finite-budget behavior. At this particular hard point their checked
unadapted densities remain smooth: adding basis84 raises the equal-mixture
density by about 1.43, while indiscriminately using all 127 lowers it to 0.66
of the original six. A fixed finite histogram does not change the strict
local power. Its trained density must not be confused with the unadapted
inverse-Jacobian representation partition. For branch probabilities πc,
densities qc and representation weights Wc, the second moment is
`∫ sum_c f² Wc²/(πc qc)`, reducing to `∫f²/qmix` only for matching balance
weights. A new channel must respect that distinction.

## Hyperspherical Jacobian correction

For Cartesian ordering `ki=r cos(θi) product_{j<i} sin(θj)` with uniformly
sampled cosines, the angular powers are `D−3−i`. Both forward and inverse
owners incorrectly used `i`, so a round-trip check alone concealed the error.
The correction changes those two powers without changing the coordinate map.
The new adjacent regression independently differentiates the actual Cartesian
map and compares its determinant with the reported Jacobian, then checks the
inverse. It covers D=3,6,9,12 and linear, log and power-1.7 radial maps;
the D12 case uses `higher_loops`. The fix and regression are published in
commit `608e61353`; all validation commands below passed.

Validation commands (six parameterization tests, including the new regression):

```bash
cargo fmt --all --check
cargo check --locked -p gammalooprs --tests --features higher_loops
cargo nextest run --locked --profile test_gammaloop -p gammalooprs --lib --features higher_loops -E 'test(utils::test_utils::test_)' --retries 0
cargo clippy --locked -p gammalooprs --all-targets --features higher_loops -- -D warnings
```

The preceding experiments deliberately used the frozen pre-correction binary
and spherical/common-radial maps, which are unaffected by this change. The
two-loop restrictions of `relative_spherical` and
`spherical_product_common_radial` were not expanded. No tropical experiment
was used.

## Reproduce and inspect the records

The portable overrides are [six-channel spherical](runtime-followup-baseline6.toml),
[seven-channel spherical](runtime-followup-union7.toml), and
[six-channel common radial](runtime-followup-common6.toml). They parse to the
same settings as the executed pilot inputs. With a compatible complete state
generated using [the current card](../../../examples/cli/epem_a_ttxh/NNLO/epem_a_tth_NNLO_test_GL638.toml),
run from the repository root, substituting the state folder:

```bash
./target/release/gammaloop --state-folder ./GL638_full/state --read-only-state --no-save-state \
  examples/cli/epem_a_ttxh/NNLO/epem_a_tth_NNLO_test_GL638.toml run -c '
  run set_model_parameters configure_runtime;
  set process -p epem_a_tth -i NNLO defaults;
  set process -p epem_a_tth -i NNLO file docs/research/gl638/runtime-followup-union7.toml;
  integrate -p epem_a_tth -i NNLO --n-cores 20 --workspace-path ./GL638_union7_30k --renderer tabled --show-phase both --show-max-weight-info --write-results-for-each-iteration'
```

Choose a fresh workspace for each comparison.
Use the same reset/override sequence for the other two profiles. Generation is
not required to select basis84 in the recorded state. Regenerated catalogues
should be checked against the edge sets, rather than assuming stable numeric IDs.

[Per-iteration records](runtime-followup-pilots.csv) contain signed/absolute
estimates, errors, extrema, precision counts and snapshot provenance.
[The pilot provenance](runtime-followup-pilot-provenance.json) preserves actual
runtime differences, budgets and state checks;
[the artifact manifest](runtime-followup-artifacts.json) records hashes and
the exact local-channel/basis/edge mappings.
[Maximum-point inputs](runtime-followup-maxima-points.json) include complete
runtime/model reset data and the five exact saved samples;
[the ten replay rows](runtime-followup-maxima-values.csv) include all cut values
and flags. Archive paths in these records identify the local executed evidence;
the embedded point and setting data are available without those archives.

[The H/Z replay inputs](runtime-followup-hz-inputs.json) contain the eight
original points and three executed mode combinations;
[the twenty rows](runtime-followup-hz-values.csv) retain totals, all six cuts,
flags and precision histories. The H/Z normal controls cover four endpoints;
the two Arb rotation choices cover all eight points.

[The A/P inputs](runtime-followup-ap-inputs.json) define five unique points,
eighteen runtime modes and the exact model/reset commands;
[the 54 rows](runtime-followup-ap-values.csv) preserve raw totals, cuts,
precision histories, errors, timings and source hashes. The mode and point IDs
remove repeated input data without changing any recorded numerical fields.
