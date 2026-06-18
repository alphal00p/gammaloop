# qq_hhh_2L Current Handoff

This directory is the active handoff bundle for the subtracted
`d d~ -> H H H` two-loop top-pentagon amplitude in GammaLoop.

Workspace root:

```text
/Users/vjhirsch/Documents/Work/gammaloop_lucien
```

Use the repository venv for Python-side GammaLoop work:

```bash
source /Users/vjhirsch/Documents/Work/gammaloop_lucien/.venv/bin/activate
```

## Current Goal

The production steering card is:

```text
examples/cli/qq_hhh_2L/qq_hhh_2L.toml
```

It is configured so that the usual pipeline:

```bash
./target/dev-optim/gammaloop --clean-state examples/cli/qq_hhh_2L/qq_hhh_2L.toml quit
./target/dev-optim/gammaloop -o -s ./examples/cli/qq_hhh_2L/state run generate_diagrams generate_integrands integrate
```

runs the optimized setup found so far:

- threshold subtraction enabled;
- threshold E-surface generation check disabled;
- generated integrands compiled with `summedFunctionMap`;
- all LMBs exposed by setting `override_lmb_heuristics = true`;
- graph Monte Carlo enabled;
- LMB-channel Monte Carlo enabled with inverse-Jacobian channel prefactors;
- `spherical_common_radial` coordinates;
- orientations explicitly summed;
- training/adaptation on the imaginary component via `integrator.integrated_phase = "imag"`;
- event and additional-weight storage disabled in the production integration block for speed;
- 3 integration cores;
- `n_start = 1000000`, `n_increase = 1000000`, `n_max = 10000000000`.

Do not enable Monte Carlo over orientations. The local IR subtraction terms only
work after explicitly summing all orientations.

## DOT Provenance

The DOT file in this directory was copied from pygloop:

```text
/Users/vjhirsch/Documents/Work/pygloop/outputs/dot_files/qqbar_nX/qqbar_nX_d_dbar_h_h_h_2L_top_pentagon_isr_subtracted.dot
```

The DOT has 12 graph blocks:

```text
GL05
GL05_isr_p1_ct
GL05_isr_p2_ct
GL06
GL06_isr_p1_ct
GL06_isr_p2_ct
GL12
GL12_isr_p1_ct
GL12_isr_p2_ct
GL13
GL13_isr_p1_ct
GL13_isr_p2_ct
```

Groups:

```text
group_id = 0: GL05, GL05_isr_p1_ct, GL05_isr_p2_ct, GL06, GL06_isr_p1_ct, GL06_isr_p2_ct
group_id = 1: GL12, GL12_isr_p1_ct, GL12_isr_p2_ct, GL13, GL13_isr_p1_ct, GL13_isr_p2_ct
```

The four original graphs are the final-state-symmetrized representatives
`GL05`, `GL06`, `GL12`, and `GL13`.

No UV subtraction is wanted here. The steering card disables GammaLoop UV
subtraction globally, and CT graph local pieces use very negative `dod` values
to avoid accidental local UV CT generation.

## Exact-xi Convention

External ordering:

```text
Q(0) Q(1) Q(2) Q(3) -> Q(4) Q(5) Q(6) Q(7) Q(8)
```

Meaning:

- `Q(0)`: incoming physical `d`;
- `Q(1)`: incoming physical `d~`;
- `Q(2)` and `Q(4)`: helper pair for the `p1` ISR CT;
- `Q(3)` and `Q(5)`: helper pair for the `p2` ISR CT;
- `Q(6)`, `Q(7)`: explicit physical Higgs momenta;
- `Q(8)`: dependent physical Higgs momentum.

The current production card uses:

```text
Q(2) = Q(3) = Q(4) = Q(5) = xi = (3000, 0, 0, 0)
```

Earlier `xi = (1000,0,0,0)` is bad for the usual `sqrt(s)=1000` setup because
it equals the total incoming momentum and creates artificial degeneracies
between helper and physical threshold surfaces.

The corrected exact-xi auxiliary denominator is:

```text
k_gi^2 -> (k_gi - xi)^2 - xi^2
```

where `k_gi` is the adjacent gluon momentum. In arXiv:2603.03169, eqs. 12 and
13 have a denominator typo: it should be `k_{g_1}-xi` and `k_{g_2}-xi`, not
`k_1-xi`.

## Physical Setup

The standard point used in the card is:

```text
sqrt(s) = 1000 GeV
MT = ymt = 173 GeV
MH = 125 GeV
WT = WH = 0
helicities = [1, -1, 0, 0, 0, 0, 0, 0, 0]
```

Higgs momenta:

```text
Q(6) = (438.55556622469447, 155.33220018353779, 348.01603965135871, -177.37736157184119)
Q(7) = (356.36963749219223, -16.802389008511, -318.72911024360047, 97.487191636880979)
Q(8) = dependent
```

## Threshold-NaN Investigation

The initial NaNs were threshold-subtraction related, not a malformed DOT:

- With `check_esurface_at_generation = true`, threshold CT construction saw
  fake exact-xi helper topology in a dangerous way.
- The event weights became huge or NaN, with additional weights dominated by
  entries such as `threshold_counterterm:*`.
- Setting `check_esurface_at_generation = false` is part of the current setup.
- Switching away from `xi = (1000,0,0,0)` removed the raw pass-on CT NaNs seen
  in focused diagnostics.

Focused xi scan at a known bad point:

```text
xi = (1000,0,0,0): raw_nonfinite = 16, coincidence cluster {8,10,16,33}
xi = (2000,0,0,0): raw_nonfinite = 0, clusters {5,8,16,31} and {10,33}
xi = (3000,0,0,0): raw_nonfinite = 0, clusters {5,31}, {8,16}, and {10,33}
xi = (2000,0,0,1000): partial/inconclusive because the run was too slow
```

The remaining `{10,33}` coincidence did not produce observed `1/0` or NaN in
the focused raw CT diagnostics. `xi = (3000,0,0,0)` gave the cleanest inspected
surface structure among the practical choices tried, so it is the production
default.

## Stability Finding

A separate instability was found in the base/original CFF evaluator at one
xi2000 point, independent of raw threshold CT generation. The rescue stack was
patched in:

```text
crates/gammalooprs/src/integrands/process/mod.rs
```

so nonfinal stability levels escalate when a rotated/orientation result is
nonfinite.

Known diagnostic point:

```text
graph group = 1
x = [0.04243809892765449, 0.6501739569562547, 0.6345802747849795,
     0.9999282161172643, 0.3850124895740257, 0.3432911025875355]
```

After the rescue-stack fix, default evaluation escalates through f64/f128 to
arb and returns:

```text
-1.9680617192257578e-9 + 1.7812184200189852e-9 i
```

This confirmed that not every NaN-looking event was a threshold CT bug; some
were normal precision-rescue failures that needed escalation.

## Integration Setup Observations

Uncompiled `SingleParametric` diagnostics:

```text
xi2000, spherical, LMB summed, 2048 points:
  -4.268903936148298e-05 - 1.0356238568615117e-03 i
  error 1.5903700965315063e-03 + 1.793993148690794e-03 i
  nonfinite fraction before rescue was about 0.195 percent

xi2000, spherical, LMB summed, 1024 probes after rescue:
  re 0.00447 +/- 0.01765
  im 0.00263 +/- 0.02438
  arb fraction about 0.195 percent

xi3000, spherical, LMB summed, 1024 probes:
  re 0.00870 +/- 0.00522
  im -0.00295 +/- 0.00352
  arb fraction about 0.391 percent

xi3000, relative-spherical, LMB summed, 1024 probes:
  re -0.000800 +/- 0.001413
  im -0.000163 +/- 0.001596
  imaginary relative error remained bad

xi3000, common-radial, LMB summed, 512 probes:
  re -0.002632 +/- 0.001604
  im 0.002562 +/- 0.002045
  arb fraction 0

xi3000, no LMB multichanneling, 512 probes:
  re -0.003053 +/- 0.002592
  im 0.001277 +/- 0.001353
```

The product common-radial map was worse. `spherical_common_radial` was the most
useful coordinate choice from the small scans.

Longer uncompiled LMB-summed xi3000/common-radial run:

```text
16384 points:
  re 0.01463 +/- 0.01215, about 83 percent relative error
  im -0.02569 +/- 0.02685, about 105 percent relative error
  NaN fraction 0
  arb fraction about 0.012 percent
  chi2 real about 8.0, imaginary about 2.59
  average total about 0.194 s/sample
```

Compiled `summedFunctionMap` observations:

```text
clean xi3000 all-LMB SFM state:
  1080 generated C++ sources
  peak RAM about 11.12 GiB

xi3000, common-radial, LMB summed, 2048 points:
  re 0.0064336 +/- 0.00464
  im -0.004373 +/- 0.00363
  NaN fraction 0
  arb fraction 0
  average total about 0.1416 s/sample

xi3000, spherical, LMB-channel MC, 4096 points:
  re 0.002346 +/- 0.01198
  im 0.007888 +/- 0.01609
  average total about 0.00753 s/sample
  variance poor

xi3000, common-radial, LMB-channel MC, 4096 points:
  re -0.004131 +/- 0.004326
  im -0.010739 +/- 0.010986
  average total about 0.00763 s/sample

xi3000, common-radial, LMB-channel MC, 32768 points:
  re 0.0004741 +/- 0.0024297
  im 0.014239 +/- 0.012576
  NaN fraction 0
  arb fraction about 0.0061 percent
  average total about 0.00705 s/sample
```

Compiled LMB-channel MC is much faster per sample than LMB-summed evaluation,
but the variance is still large. The main remaining issue is variance from
stable high-weight structures, not NaNs or precision instability.

## Current Rebuild And xi Scan

After the threshold-CT/LMB ordering fixes, both interfaces must be rebuilt
before trusting new diagnostics:

```bash
cargo clean -p gammalooprs -p gammaloop-api
just build-cli
just build-api-wheel
source /Users/vjhirsch/Documents/Work/gammaloop_lucien/.venv/bin/activate
python -m pip install --force-reinstall --no-deps target/wheels/gammaloop-0.3.3-cp312-abi3-macosx_11_0_arm64.whl
```

The local Python binding currently exposes `run`, `evaluate_sample`, and
`evaluate_samples`, but not `GammaLoopAPI.integrate`; integrations are therefore
launched through the rebuilt CLI and inspected through the Python API.

The active monitoring state was regenerated from the rebuilt CLI with
non-compiled `SingleParametric` evaluators:

```bash
./target/dev-optim/gammaloop --clean-state examples/cli/qq_hhh_2L/qq_hhh_2L.toml quit
./target/dev-optim/gammaloop -o -s ./examples/cli/qq_hhh_2L/state \
  examples/cli/qq_hhh_2L/qq_hhh_2L.toml \
  run generate_diagrams generate_integrands_single_parametric
```

Low-stat xi scan driver:

```text
examples/cli/qq_hhh_2L/scan_xi_variants_lowstat.py
```

It varies the helper vector as `xi = t * (p1 + p2) = (1000 t, 0, 0, 0)` while
keeping graph/LMB Monte Carlo, `spherical_common_radial`, inverse-Jacobian LMB
weights, threshold subtraction on, and orientations summed.

1024-point non-compiled scan:

```text
t = 0.5:
  re = -9.5560e-03 +/- 1.21e-02  (126.2%)
  im = -1.2544e-02 +/- 1.38e-02  (109.7%)
  im max-weight impact = 1.09
  im max range = [-14.0, +0.890]
  arb = 0.0977%, nans+unstable = 0

t = 1.0:
  pathological for this setup; killed after more than two minutes with no
  first-iteration result, while other xi values completed promptly.
  Avoid exact xi = p1 + p2.

t = 1.1:
  re = +1.5913e-01 +/- 1.29e-01  (81.3%)
  im = +5.3450e-02 +/- 4.58e-02  (85.7%)
  im max-weight impact = 0.829
  im max range = [-6.28, +45.4]
  arb = 0.0977%, nans+unstable = 0

t = 2.0:
  re = +1.3883e-02 +/- 1.24e-02  (89.5%)
  im = -1.2168e-02 +/- 2.47e-02  (203.1%)
  im max-weight impact = 1.70
  im max range = [-21.2, +13.4]
  arb = 0.0977%, nans+unstable = 0

t = 3.0:
  re = -5.5502e-03 +/- 6.93e-03  (124.8%)
  im = -4.8240e-03 +/- 5.71e-03  (118.3%)
  im max-weight impact = 1.17
  im max range = [-5.76, +0.824]
  arb = 0.0977%, nans+unstable = 0
```

4096-point non-compiled refinement for the most relevant choices:

```text
t = 1.1:
  re = +6.0642e-02 +/- 3.72e-02  (61.4%)
  im = +4.7175e-03 +/- 1.32e-02  (279.4%)
  im max-weight impact = 2.35
  im max range = [-23.5, +45.4]
  arb = 0.0244%, nans+unstable = 0

t = 3.0:
  re = -2.9373e-02 +/- 2.71e-02  (92.2%)
  im = -6.2011e-02 +/- 5.37e-02  (86.5%)
  im max-weight impact = 0.861
  im max range = [-219, +4.90]
  arb = 0.0488%, nans+unstable = 0
```

Interpretation:

- Exact `t = 1` remains a bad helper-vector choice.
- `t = 2` was finite but had poor imaginary variance in the 1024-point scan.
- `t = 1.1` is finite and useful as a near-degenerate stress test, but the
  4096-point refinement did not improve the imaginary variance.
- `t = 3` remains the best practical default among these checks.
- No xi value that completed showed NaNs or unstable samples.

The 4096-point `t = 3` negative-imaginary max point was approach-scanned:

```text
examples/cli/qq_hhh_2L/maxweight_approach_scans/xi_t3_4096_approach/
```

The integration max table reported `im[-] = -218.78020604861075`, but the
fixed-bin unit-sampling value at the same coordinates was only:

```text
-2.894425863126332 - 5.757373843384493 i
```

The largest sampled fixed-bin absolute values along the three approach lines
were:

```text
axis_x1          ~= 39.93
axis_x4          ~= 69.06
center_diagonal  ~= 32.14
```

The arb relative delta at the peak was about `6.2e-12`, and no scan row was
NaN. This does not look like the old incorrect-threshold-subtraction failure;
the larger integration-table number is again driven by sampling/adaptive
weights rather than an unresolved fixed-bin singularity.

## Max-Weight Inspects

Uncompiled xi3000/common-radial/LMB-summed 16384 max points were finite and
agreed with arb:

```text
graph group 1
x = [0.5340123207153026, 0.18929422980978153, 0.27293481721650031,
     0.96463329204710657, 0.27286422732872495, 0.076919046071880992]
default = 76.163074005184839 + 45.880456999624741 i

graph group 0
x = [0.20214898059753, 0.4873218279157891, 0.79935255628431778,
     0.99904552609362041, 0.51217914486600502, 0.76630827600601747]
default = 3.8839380452326684 - 18.265505503870678 i
```

Compiled xi3000/common-radial/LMB-channel MC 32768 max point:

```text
graph group = 0
LMB channel = 14
x = [0.30784578180460542, 0.90170572955766448, 0.71240971451878965,
     0.75722185209864434, 0.36452634459630484, 0.265760901346222]
default = -1.3627481340754974 + 10.782175613902147 i
arb     = -1.3627481340613690 + 10.782175613790022 i
```

Cluster run max-weight points with imaginary training phase:

```text
Global maxima:
  re[+]: +3.2657068339264096e5
    graph = 0
    LMB channel = 3
    x = [0.37277353160468185, 0.39959369059694888, 0.65556145913795516,
         0.62413302507246593, 0.65872296153231757, 0.77736662361123232]

  re[-]: -5.8344662283583556e5
    graph = 0
    LMB channel = 11
    x = [0.27543703182272550, 0.025975374287740258, 0.50782398548203456,
         0.68942356972268137, 0.80460957897792340, 0.52443570508794690]

  im[+]: +3.5795025591653502e6
    graph = 0
    LMB channel = 11
    x = [0.27543703182272550, 0.025975374287740258, 0.50782398548203456,
         0.68942356972268137, 0.80460957897792340, 0.52443570508794690]

  im[-]: -7.7414299360872700e5
    graph = 1
    LMB channel = 9
    x = [0.23262693392408407, 0.00000034231896575351528, 0.49485992364488252,
         0.99442214604999712, 0.33719614224884459, 0.78136672187917255]

Imaginary maxima by graph group:
  group #0 im[+]: +1.0696548987873162e6
    graph = 0
    LMB channel = 11
    x = [0.27543703182272550, 0.025975374287740258, 0.50782398548203456,
         0.68942356972268137, 0.80460957897792340, 0.52443570508794690]

  group #0 im[-]: -1.7707154593472226e5
    graph = 0
    LMB channel = 18
    x = [0.26378766476195209, 0.87106269195833919, 0.30334636296707945,
         0.32983414216533991, 0.26734416329469662, 0.48596502944428138]

  group #1 im[+]: +1.5144411558846670e4
    graph = 1
    LMB channel = 2
    x = [0.34335711109144562, 0.0000017072811633777772, 0.88839091396443770,
         0.99801107280213031, 0.68188596968339055, 0.057247742844868046]

  group #1 im[-]: -1.9216672554917634e5
    graph = 1
    LMB channel = 9
    x = [0.23262693392408407, 0.00000034231896575351528, 0.49485992364488252,
         0.99442214604999712, 0.33719614224884459, 0.78136672187917255]
```

The max points were numerically stable. Additional weights were finite. One
representative large threshold CT contribution was of order:

```text
threshold_counterterm:11:0 ~= -1.93e-17 + 2.27e-16 i
full factor ~= 4.5066e16
```

This points to an integrable high-variance region rather than a remaining NaN
bug.

## Max-Weight Approach Scans

The script:

```text
examples/cli/qq_hhh_2L/scan_imag_maxweight_approaches.py
```

scans three 201-point x-space lines through the two global imaginary cluster
max-weight points listed above. It writes:

```text
examples/cli/qq_hhh_2L/maxweight_approach_scans/imag_maxweight_approach_scan.json
examples/cli/qq_hhh_2L/maxweight_approach_scans/imag_maxweight_approach_scan.tsv
examples/cli/qq_hhh_2L/maxweight_approach_scans/imag_maxweight_approach_inspect_spotchecks.log
examples/cli/qq_hhh_2L/maxweight_approach_scans/imag_plus_global_approaches.png
examples/cli/qq_hhh_2L/maxweight_approach_scans/imag_minus_global_approaches.png
```

Important convention: these scans are fixed-bin diagnostics with explicit
`graph, LMB channel` selection and unit sampling weight. The cluster maximum
table also included the adaptive sampling/grid weight, so its normalization is
not expected to match the Python API fixed-bin scan directly. In the TSV/JSON:

```text
raw_re/raw_im/raw_abs = raw evaluator output before the parameterization Jacobian
re/im/abs             = raw output times parameterization_jacobian and integrator_weight
weight_convention     = parameterization_jacobian_x_unit_fixed_bin_sampling
```

The scanner uses non-compiled `SingleParametric` evaluation and disables the
polarization cache for these repeated point probes. With the cache enabled, a
debug `ParamCache::checked_push` assertion can be hit during the dense fixed-bin
scan; disabling it avoids a diagnostic-only cache bookkeeping issue without
changing the fixed-bin numerical value. This does not change the production
integration card, which remains configured for compiled `SummedFunctionMap`.

Peak arb spot checks agree well:

```text
imag_minus_global fixed-bin peak:
  default = +1.28243235e6 - 1.79737247e7 i
  arb relative delta ~= 8.8e-9

imag_plus_global fixed-bin peak:
  default = -1.43925284e4 + 8.82995810e4 i
  arb relative delta ~= 1.0e-7
```

Shape observations from the line scans:

```text
imag_plus_global:
  The fixed-bin weight is localized near the cluster point.
  The largest sampled |weight| on the three lines was:
    axis x1          ~= 1.36e6 near signed fraction +6.6e-9
    axis x4          ~= 1.45e6 near signed fraction -4.3e-9
    center diagonal  ~= 6.61e6 near signed fraction -1.2e-8

imag_minus_global:
  The fixed-bin peak itself is finite and stable, but the lines are dominated
  by much larger boundary growth away from the cluster point:
    peak |weight|    ~= 1.80e7
    axis x1 boundary ~= 8.63e12
    axis x3 boundary ~= 1.00e20
    diagonal boundary ~= 8.64e12
```

This means the cluster max-weight point is not necessarily the largest
fixed-bin structure once the adaptive sampling weight is removed. For variance
work, inspect both the cluster max points and the boundary approaches exposed by
these scans.

## Raw Threshold CT Diagnosis

The script:

```text
examples/cli/qq_hhh_2L/inspect_imag_maxweight_thresholds.py
```

uses the Python API to evaluate the cluster imaginary max-weight points and the
largest points found in the approach scans with event generation and
`additional_weights` enabled. The output is:

```text
examples/cli/qq_hhh_2L/maxweight_approach_scans/threshold_decomposition/threshold_decomposition.md
examples/cli/qq_hhh_2L/maxweight_approach_scans/threshold_decomposition/threshold_decomposition.json
examples/cli/qq_hhh_2L/maxweight_approach_scans/threshold_decomposition/threshold_decomposition.tsv
```

Both cluster imaginary max-weight points are dominated by one threshold CT, not
by the original integrand:

```text
imag_plus_global, graph group 0, LMB channel 11:
  dominant graph term = GL05_isr_p1_ct
  dominant additional weight = threshold_counterterm:11:0
  threshold edges = [9, 11, 14]
  fixed-bin weighted value ~= -1.43925271e4 + 8.82995520e4 i
  original weighted contribution ~= 2.24e-2 in absolute value

imag_minus_global, graph group 1, LMB channel 9:
  dominant graph term = GL12
  dominant additional weight = threshold_counterterm:23:0
  threshold edges = [11, 14, 15]
  fixed-bin weighted value ~= +1.28110283e6 - 1.79725935e7 i
  original weighted contribution ~= 5.58e-2 in absolute value
```

Selected raw pass-on CT evaluator dumps were written with:

```text
GAMMALOOP_DUMP_RAW_THRESHOLD_CT="GL12:23,GL05_isr_p1_ct:11"
```

to:

```text
examples/cli/qq_hhh_2L/maxweight_approach_scans/raw_threshold_ct_dumps/
```

At the exact cluster points, the selected raw pass-on evaluators are finite:

```text
GL05_isr_p1_ct, threshold 11:
  raw_nonfinite = false
  raw_result ~= +1.459157566e9 - 2.290372760e9 i
  h_left_th ~= 2.508900088e-3
  r_left ~= 455.9476158307
  r*_left ~= 449.3709013577

GL12, threshold 23:
  raw_nonfinite = false
  raw_result ~= -2.234829712e21 + 3.875507162e21 i
  h_left_th ~= 2.405299977e-3
  r_left ~= 469.8723087812
  r*_left ~= 469.1171583616
```

The raw expressions are products of regular factors and inverse pass-on energy
denominators. Ranking the `(... )^(-1)` factors numerically shows the concrete
origin of the large weights.

For `GL12` threshold 23 `[11, 14, 15]`:

```text
E9  ~= 499.999898
E10 ~= 1.0360623e-4
E11 ~= 1.0360623e-4

E9 + E11 - Q0^0 ~= 1.1558011e-6
E9 + E10 - Q1^0 ~= 1.1558011e-6
E10 ~= 1.0360623e-4
E11 ~= 1.0360623e-4
```

These four small factors alone give a very large finite amplification. This
region is close to the selected threshold 23, but the worst raw amplification is
caused by overlap with other near-on-shell structures involving the soft
massless lines 10 and 11. The corresponding active threshold edge sets are also
present in the group, for example `[9, 11]` and `[9, 10]`.

For `GL05_isr_p1_ct` threshold 11 `[9, 11, 14]`:

```text
E9 ~= 406.309799
E10 ~= 9.00698195
E11 ~= 403.857397
E12 ~= 3000.01352
E15 ~= 237.307047

E9 + E12 + E15 + Q7^0 - Q0^0 - Q1^0 - Q2^0 ~= 4.0549328e-6
```

This is again a different near-on-shell pass-on denominator inside the selected
threshold CT evaluator. The related active threshold edge set `[9, 12, 15]`
appears in the group.

The active-threshold metadata also shows duplicate edge sets for the problematic
surfaces:

```text
group 0: [9, 11, 14] appears as threshold IDs 11 and 24
group 1: [11, 14, 15] appears as threshold IDs 23 and 37
```

So duplicate/signature-degenerate E-surfaces are present. However, a selected
dump of the duplicate partners:

```text
GAMMALOOP_DUMP_RAW_THRESHOLD_CT="GL12:37,GL05_isr_p1_ct:24"
```

at the same two spot-checks produced no raw dump files, while IDs 11 and 23 do
dump and dominate the event-level additional weights. The inspected large
finite weights are therefore better explained by overlapping pass-on
denominators inside one selected threshold CT. At these points the problem is
not a raw `0/0`, and increasing precision is not expected to remove the peak.
It is a threshold-subtraction variance problem near overlapping E-surface
regions.

Useful next checks:

- compare the paired duplicate surfaces, especially IDs 11/24 and 23/37, by
  dumping both members at the same x-space points;
- extend the approach scans along the small-denominator directions above;
- consider damping or overlap handling for threshold CTs whose pass-on
  expression is close to other E-surfaces;
- inspect whether the selected threshold CT should be suppressed or combined
  differently when `E10`/`E11` become soft in the graph 1/LMB 9 region.

Follow-up cancellation check:

```text
imag_minus_global, dominant event GL12:
  threshold_counterterm:23:0, [11,14,15]
    weighted ~= +1.28110283e6 - 1.79725935e7 i

  candidate dual partners for the small pass-on denominators:
    threshold IDs 0 and 12, [9,11]: missing in this event
    threshold IDs 6 and 9, [9,10]: missing in this event
    duplicate edge-set partner ID 37, [11,14,15]: missing in this event

imag_plus_global, dominant event GL05_isr_p1_ct:
  threshold_counterterm:11:0, [9,11,14]
    weighted ~= -1.43925271e4 + 8.82995520e4 i

  candidate dual partners for the small pass-on denominator:
    threshold ID 4, [9,12,15]: missing in this event
    threshold ID 33, [9,12,15]:
      weighted ~= +1.00e-10 - 5.31e-12 i
    duplicate edge-set partner ID 24, [9,11,14]: missing in this event
```

This is the important failure pattern: the CT that blows up near the overlapping
surface is present, but the CT associated with the overlapping surface is either
not present in the selected graph-term overlap structure or is many orders of
magnitude too small to cancel it.

`profile bulk` selector checks with `SingleParametric` runtime:

```text
GL12 T(t23): accepted
GL12 T(t0): rejected, threshold does not exist in selected overlap structure
GL12 T(t6): rejected, threshold does not exist in selected overlap structure

GL05_isr_p1_ct T(t11): accepted
GL05_isr_p1_ct T(t4): rejected, threshold does not exist in selected overlap structure
GL05_isr_p1_ct T(t33): accepted
```

The rejected selectors are precisely the natural partners for the small
pass-on denominators identified above. This strongly suggests an incomplete
threshold-overlap subtraction/cancellation setup for these graph terms, rather
than a precision issue.

Correction from a deeper current-code trace:

```text
The dangerous graph-0 pass-on denominator is local to GL05_isr_p1_ct and is:

  D_B = -Q0^0 - Q1^0 - Q2^0 + Q7^0 + E9 + E12 + E15

The matching local GL05_isr_p1_ct E-surfaces are:

  t13, group_esurface_id 162, edges [9,12,15]
  t32, group_esurface_id 180, edges [9,12,15]

At the imag_plus_global point, with selected CT T(t11), the overlap/existence
check rejects both:

  t13: candidate_exists=false, shift_part=-1000, mass_sum=3346,
       existence_margin=-1.0195716e7
  t32: candidate_exists=false, shift_part=+1000, mass_sum=3346,
       existence_margin=-1.0195716e7

The rstar check at r = r*_t11 also only marks t11 coincident:

  T(t11): value = 0
  T(t13): value ~= +2.8101807165639723e3
  T(t32): value ~= +4.8101807165639723e3

However the raw pass-on evaluator contains the raw_atom with +Q7^0, while the
CFF existence atom for t13/t32 uses a helper external slot P5 in the
corresponding shift. Numerically, the raw denominator is near zero at the same
rstar point:

  -Q0^0 - Q1^0 - Q2^0 + Q7^0 + E9 + E12 + E15 ~= 4e-6

Using the raw shifted surface gives a positive approximate existence margin,
while the current CFF existence object gives a negative one. This means the
partner surface is not physically absent; the existence criterion is being
applied to an object whose external-shift mapping is not synchronized with the
raw CT pass-on evaluator.
```

Current-code fix applied in this workspace:

```text
Loop-momentum bases are now canonicalized so that external slots are ordered by
the external half-edge ids. This keeps the numerical external momenta aligned
between DOT parsing, diagram generation, threshold-existence checks, and local
threshold CT evaluation.

Touched areas:
  crates/gammalooprs/src/graph/lmb.rs
  crates/gammalooprs/src/graph/mod.rs
  crates/gammalooprs/src/graph/feynman_graph.rs
  crates/gammalooprs/src/graph/parse/mod.rs
  crates/gammalooprs/src/feyngen/diagram_generator.rs
  crates/gammalooprs/src/processes/cross_section.rs
```

Validation after the LMB external-order fix:

```text
cargo fmt --check: passed
git diff --check: passed
cargo test --no-run --profile dev-optim -p gammalooprs: passed
just build-cli: passed

just test_gammaloop --no-fail-fast:
  still blocked before running tests by the pre-existing spenso dead-code
  warning under -Dwarnings.

just test_gammaloop --no-warnings-as-errors --no-fail-fast:
  same baseline result before and after the fix:
  1331 tests run, 1279 passed, 52 failed, 129 skipped.
```

Focused old-cluster max-point checks after the fix:

```text
imag_minus_global old point, graph group 1, LMB channel 9:
  fixed-bin evaluation = +1.2824323595646839e6 - 1.7973724813762158e7 i
  raw integrand        = +1.5940391908503898e1 - 2.2341000322560794e2 i
  jacobian             = +8.0451745912252663e4
  arb-stable, no huge 1e269 threshold CTs, no NaN.

imag_plus_global old point, graph group 0, LMB channel 11:
  finite and stable; no huge 1e269 threshold CTs, no NaN.
```

This addresses the original NaN/missing-partner failure mode. The old cluster
max-weight points can still be large fixed-bin points, but they are no longer
catastrophic nonfinite threshold-subtraction events.

Post-fix compiled reintegration diagnostic:

```text
Command block:
  integrate_reintegration_200k

State:
  clean compiled xi3000 summedFunctionMap state generated with the production
  generate_diagrams/generate_integrands setup.

Runtime:
  3 cores
  imaginary training phase
  graph MC
  LMB-channel MC with inverse-Jacobian channel prefactors
  spherical_common_radial
  orientations summed
  n_start = 200000
  n_increase = 200000
  interrupted cleanly after two local iterations

Iteration 1, 200000 cumulative evaluations:
  re = -3.1083970011e-03 +/- 4.711e-03  (151.6 percent)
  im = -4.2744536294e-03 +/- 7.390e-03  (172.9 percent)
  f128 = 0.0045 percent
  arb  = 0.0055 percent
  NaN/unstable = 0

Iteration 2, 600000 cumulative evaluations:
  re = +3.1887943101e-02 +/- 3.595e-02  (112.7 percent)
  im = +5.9074670078e-03 +/- 1.738e-02  (294.1 percent)
  f128 = 0.0025 percent
  arb  = 0.0051666667 percent
  NaN/unstable = 0
```

Post-fix integration max weights after two iterations:

```text
re[+]:
  value = +2.116318e4
  graph group = 0
  LMB channel = 14
  x = [0.22215871102735293, 0.95151497131053109,
       0.10251622320444055, 0.17567342269336178,
       0.34693037158790663, 0.51463740505783795]

im[+]:
  value = +9.769326e3
  same graph group, LMB channel, and x as re[+]

re[-]:
  value = -3.601238e3
  graph group = 0
  LMB channel = 17
  x = [0.33749134101414835, 0.52728630237801777,
       0.36580145953442500, 0.87532168016194667,
       0.59155973698535502, 0.74659416774463649]

im[-]:
  value = -3.075815e3
  graph group = 1
  LMB channel = 16
  x = [0.57661611487341058, 0.045616448806269319,
       0.77070029584810219, 0.99532604805865467,
       0.17460243684638546, 0.10034458098675393]
```

The largest observed post-fix weight is orders of magnitude below the previous
cluster-scale values and is finite. The corresponding focused inspect with
event generation and additional weights enabled was:

```text
Fixed bin:
  graph group = 0
  LMB channel = 14

Fixed-bin evaluation:
  +1.8866605925766407e2 + 8.7091833076629058e1 i

Raw integrand before the parameterization jacobian:
  +1.7870952341180687e-13 + 8.2495707195160286e-14 i

Parameterization jacobian:
  +1.0557135157420459e15

Dominant event:
  GL06_isr_p1_ct
  LMB channel = 14
  basis = (9,12)
  event weight = +1.8791086792864749e2 + 8.5348817911559877e1 i

Dominant additional weight in that event:
  threshold_counterterm:17:0
  +1.8839120169205068e-13 + 8.5574689204001705e-14 i
```

Debug trace for this post-fix max point identifies the dominant local CT as:

```text
graph = GL06_isr_p1_ct
group_esurface_id = 35
local esurface_id = 17
raised_esurface_id = 17
edges = [12,14]
raw_atom = -Q(0,cind(0)) - Q(2,cind(0)) + OSE(12) + OSE(14)
candidate_exists = true
existence_margin ~= +3.0e6
mass_sum ~= +3.0e3
shift_part ~= -3.5e3
```

Only the self-coincident threshold is selected in the rstar coincidence list
for this CT. The raw pass-one evaluator is large, but the pass-two/local CT
result is finite and arb-stable. In this post-fix run the observed max weight
looks like an integrable finite threshold-CT variance spike amplified by the
common-radial jacobian and adaptive sampling/discrete weights, not the old
wrongly shifted missing-E-surface NaN.

## Useful Logs And Workspaces

Most diagnostic work was written under:

```text
examples/cli/qq_hhh_2L_experiment_xi2000/
```

Useful files:

```text
debug_logs/xi_choice_scan_summary.json
debug_logs/integration_compiled_xi3000_lmb_mc/summary.tsv
debug_logs/max_weight_inspects_compiled_lmb_mc_32768/summary.tsv
workspaces/compiled_xi3000_lmb_mc/compiled_xi3000_common_radial_lmb_mc_32768/integration_result.json
```

Post-fix reintegration logs and workspace:

```text
/tmp/qq_hhh_reintegration_generate.log
/tmp/qq_hhh_reintegration_200k.log
/tmp/qq_hhh_reintegration_maxweight_g0_c14_inspect.log
/tmp/qq_hhh_reintegration_maxweight_g0_c14_debug_inspect.log
/tmp/qq_hhh_display_integrands_after_patch.log

examples/cli/qq_hhh_2L/workspaces/reintegration_patch_xi3000_lmb_mc_200k/
examples/cli/qq_hhh_2L/state/logs/gammalog-reintegration_g0c14.jsonl.jsonl
```

The old `examples/cli/qq_hhh_2L/state` may be stale if it was generated before
the xi3000/check-disabled setup. Regenerate cleanly before trusting it.

## Next Useful Work

The pointwise integrand is now usable. The remaining path to percent-level
accuracy is variance reduction.

Recommended next diagnostics:

- continue longer post-fix integrations; the short local 600k cumulative run is
  not statistically meaningful for a central value, but it shows no
  NaN/unstable evaluations;
- keep monitoring max weight and top discrete bins during production runs;
- inspect any new max-weight point with event generation and additional weights
  enabled before treating it as a threshold bug;
- if the graph 0/LMB channel 14 structure recurs, study the
  `GL06_isr_p1_ct` local threshold `edges=[12,14]` and whether its finite
  pass-two CT variance can be improved by mapping/channel changes;
- keep the old cluster imaginary max points as regression checks for the
  external-order/LMB canonicalization fix.

Do not disable threshold subtraction in production. It is useful only as an
isolation test.

## Latest t=4 Sliver30 Run

After the t-scan indicated that `t = 4` was competitive with or better than the
previous `t = 3` default, a longer local run was performed with the rebuilt
non-compiled `SingleParametric` setup:

```text
xi = 4 * (p1 + p2) = (4000, 0, 0, 0)
sliver_width = 30
gaussian_width = 1
threshold h_function_settings.sigma = 1
neval = 81920
cores = 3
```

The run took about 40 minutes locally and wrote:

```text
examples/cli/qq_hhh_2L/workspaces/xi_t4_sliver30_hsigma1_81920/xi_t4_81920/results/integration_result_iter_0001.json
examples/cli/qq_hhh_2L/maxweight_approach_scans/xi_t4_sliver30_hsigma1_81920/xi_lowstat_summary.tsv
```

Central value and stability:

```text
re = +9.04527952218285e-3 +/- 9.285267341099707e-3  (102.65%)
im = +2.7017558565915895e-3 +/- 2.3609704868351104e-3 (87.39%)

f64 = 99.8876953125%
f128 = 0.001220703125%
arb = 0.111083984375%
nan_or_unstable = 0%
```

The max weights from this run are not large compared to the earlier broken
threshold-subtraction behavior:

```text
re[+] = +742.8988425547569
  graph = 0, LMB channel = 10
  x = [0.38598909479699239, 0.70666856531181010, 0.82233679407730331,
       0.89665017350639320, 0.52597618726715933, 0.20763469492626785]

re[-] = -73.36142376299647
  graph = 0, LMB channel = 16
  x = [0.46933350530084328, 0.71053972047096614, 0.19689946281369208,
       0.94155137245600085, 0.11454738293328448, 0.73494706196444592]

im[+] = +89.10839395271918
  graph = 0, LMB channel = 14
  x = [0.64490788982327474, 0.84548394414228900, 0.32619750670511971,
       0.60929269967474309, 0.10565752732380618, 0.65627442275460857]

im[-] = -62.98723405864441
  graph = 0, LMB channel = 9
  x = [0.53875892494388933, 0.24174091659991892, 0.94411904877480457,
       0.75278429956824044, 0.55851127303951287, 0.25621170988207465]
```

Fixed-bin approach scans through the two imaginary max-weight points were run
with the same `xi`, sliver, gaussian width, and threshold `h_sigma` settings:

```text
examples/cli/qq_hhh_2L/maxweight_approach_scans/xi_t4_sliver30_hsigma1_81920_approach/
```

Artifacts:

```text
highstat_iter_0001_imag_plus_approaches.png
highstat_iter_0001_imag_minus_approaches.png
imag_maxweight_approach_scan.tsv
threshold_decomposition/threshold_decomposition.md
```

The fixed-bin scan found no NaNs. Arb checks at the exact max points agree well:

```text
im[+] peak fixed-bin weight = -1.5880861774984467 + 2.3449577355978732 i
  |w| = 2.8321095474755094
  arb relative delta ~= 1.0e-12

im[-] peak fixed-bin weight = -0.6478055190454123 - 1.6575587910169582 i
  |w| = 1.7796497228902366
  arb relative delta ~= 1.4e-9
```

Largest fixed-bin absolute values found along the approach lines:

```text
im[+]:
  axis_x1          max |w| ~= 16.91
  axis_x4          max |w| ~= 95.99
  center_diagonal  max |w| ~= 16.42

im[-]:
  axis_x1          max |w| ~= 10.83
  axis_x2          max |w| ~= 15.64
  center_diagonal  max |w| ~= 59.97
```

Threshold decomposition at the exact fixed-bin max points shows modest
threshold CT dominance but no runaway:

```text
im[+] fixed bin:
  dominant graph term = GL06_isr_p2_ct
  largest threshold CT = threshold_counterterm:36:0, edges [9,14,16]
  weighted CT ~= -3.21e-1 + 4.16e-1 i, |w| ~= 5.25e-1

im[-] fixed bin:
  dominant graph term = GL05_isr_p2_ct
  largest threshold CT = threshold_counterterm:4:0, edges [9,12,15]
  weighted CT ~= -3.20e-1 - 7.59e-1 i, |w| ~= 8.24e-1
```

Interpretation: this t=4/sliver30 run does not show the old missing/incorrect
threshold-subtraction failure. The fixed-bin approach behavior is bounded and
arb-stable. The integration-table max weights are mostly sampling/adaptive-grid
amplification, not a large unresolved local threshold CT spike.

### axis_x2 CT-by-CT monitor

The sharp-looking feature in the `im[-]` `axis_x2` approach scan was monitored
with individual event additional weights:

```text
examples/cli/qq_hhh_2L/monitor_axis_x2_ct_contributions.py
examples/cli/qq_hhh_2L/maxweight_approach_scans/xi_t4_sliver30_hsigma1_81920_approach/axis_x2_ct_monitor.md
examples/cli/qq_hhh_2L/maxweight_approach_scans/xi_t4_sliver30_hsigma1_81920_approach/axis_x2_ct_monitor_points.tsv
examples/cli/qq_hhh_2L/maxweight_approach_scans/xi_t4_sliver30_hsigma1_81920_approach/axis_x2_ct_monitor_terms.tsv
examples/cli/qq_hhh_2L/maxweight_approach_scans/xi_t4_sliver30_hsigma1_81920_approach/axis_x2_ct_monitor_abs.png
examples/cli/qq_hhh_2L/maxweight_approach_scans/xi_t4_sliver30_hsigma1_81920_approach/axis_x2_ct_monitor_components.png
```

The reference root from the rstar probe is:

```text
x2_root = 0.9441181235200773
```

The monitored term set is identical at every sampled point; no threshold CT
appears or disappears across the drop. The feature is instead saturated by four
present threshold CTs:

```text
GL05_isr_p2_ct threshold_counterterm:4:0   edges [9,12,15]
GL05_isr_p1_ct threshold_counterterm:34:0  edges [15,16]
GL05           threshold_counterterm:31:0  edges [9,10,15]
GL05           threshold_counterterm:5:0   edges [9,12,16]
```

These four terms flip sign at `x2_root` and decay coherently like
`1 / |x2 - x2_root|`. The sum of all other event terms is about `8.7e-5` in
absolute value around the feature, while the dominant CT sum ranges from order
one on the saved scan to order `1e4` very close to the root. Thus this specific
drop is not an existence-condition switch; it is a smooth near-pole crossing
where several CT contributions add with the same phase instead of cancelling.
