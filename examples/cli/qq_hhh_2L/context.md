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

Superseded implementation hypothesis from the earlier trace:

```text
The first proposed fix was to canonicalize all generated loop-momentum bases so
external slots were ordered by external half-edge ids.

This was not the final fix retained in this workspace. The actual source of the
bad shifted symbolic surfaces was narrower: the symbolic CFF surface atoms
generated in an LMB were using external carrier edge ids instead of user
external slots.

Superseded touched-area list from that hypothesis:
  crates/gammalooprs/src/graph/lmb.rs
  crates/gammalooprs/src/graph/mod.rs
  crates/gammalooprs/src/graph/feynman_graph.rs
  crates/gammalooprs/src/graph/parse/mod.rs
  crates/gammalooprs/src/feyngen/diagram_generator.rs
  crates/gammalooprs/src/processes/cross_section.rs
```

Current-code fix applied in this workspace:

```text
E- and H-surface `to_atom_in_lmb` now translate external shifts through the
LMB edge signature into compact external slots Q0, Q1, ... rather than through
`lmb.ext_edges[external_index]`. This keeps the symbolic pass-one CFF
denominators aligned with the numerical external momenta supplied by the user.

Amplitude threshold CT generation must use the original/non-generated-basis
residue selector for this target:

  threshold_residue_in_generated_basis = false

The generated-basis selector was useful as a diagnostic, but in the current
code path it produces larger raw pass-one residues at the engineered root.

Touched areas for the retained fix:
  crates/gammalooprs/src/cff/esurface.rs
  crates/gammalooprs/src/cff/generation.rs
  crates/gammalooprs/src/cff/hsurface.rs
  crates/gammalooprs/src/cff/mod.rs
  crates/gammalooprs/src/cff/surface.rs
  crates/gammalooprs/src/integrands/process/amplitude/mod.rs
  crates/gammalooprs/src/processes/amplitude.rs
  crates/gammalooprs/src/subtraction/overlap.rs
```

Validation after the LMB external-order/local-existence fix:

```text
All commands below were run through
  python3 run_with_memory_watch.py --limit-gb 30 -- ...

cargo fmt --check: passed, peak RSS 0.10 GiB
git diff --check: passed
cargo check -p gammalooprs: passed, peak RSS 1.18 GiB
cargo test -p gammalooprs to_atom_in_lmb_uses_external_slots_not_carrier_edges -- --nocapture:
  passed, 2 tests, peak RSS 5.59 GiB
just build-cli:
  passed, peak RSS 7.35 GiB
just build-api-wheel plus .venv reinstall:
  passed, peak RSS 7.59 GiB
qq_hhh_2L SingleParametric scratch regeneration in state_watchdog_probe:
  passed, peak RSS 0.82 GiB

just test_gammaloop --no-fail-fast:
  launched under the watchdog and reached test 1339/1339, but stop.order
  interrupted the last UV test and the harness exited with code 130. The
  partial run showed broad unrelated failures, including generated example
  card discovery under qq_hhh workspaces, missing large-spenso fixture files,
  UV/integration tolerance failures, and one symjit assertion. Treat this as
  an incomplete baseline audit, not a clean regression signal for this fix.
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
GL05_isr_p1_ct threshold_counterterm:34:0  edges [9,10,15]
GL05           threshold_counterterm:31:0  edges [9,10,15]
GL05           threshold_counterterm:5:0   edges [9,12,16]
```

These four terms flip sign at `x2_root` and decay coherently like
`1 / |x2 - x2_root|`. The sum of all other event terms is about `8.7e-5` in
absolute value around the feature, while the dominant CT sum ranges from order
one on the saved scan to order `1e4` very close to the root. Thus this specific
drop is not an existence-condition switch; it is a smooth near-pole crossing
where several CT contributions add with the same phase instead of cancelling.

### T4/T34 dual-cancellation deep dive

The engineered intersection probe for the graph-0 group uses:

```text
x_root = [0.5387589249438893, 0.24174091659991892,
          0.9441181235200773, 0.7527842995682404,
          0.5585112730395129, 0.25621170988207465]
graph group = 0
LMB channel = 9
```

At this point, the two dangerous E-surfaces have the same local radial
projection:

```text
GL05_isr_p2_ct threshold_counterterm:4:0
  group E-surface = 196
  local raised/local E-surface = 4/4
  edges = [9,12,16]
  r = 1061.8603586938355
  rstar = 1061.8603586937081
  r - rstar = 1.2733e-10
  deta/dr = +2.3502742820730385

GL05_isr_p1_ct threshold_counterterm:34:0
  group E-surface = 192
  local raised/local E-surface = 34/34
  edges = [9,10,15]
  r = 1061.8603586938355
  rstar = 1061.8603586937081
  r - rstar = 1.2733e-10
  deta/dr = +2.3502742820730385
```

The raw pass-one residues at the same point have the same phase:

```text
GL05_isr_p2_ct T4  raw ~= -1.2920e-3 + 3.2616e-3 i
GL05_isr_p1_ct T34 raw ~= -1.1832e-3 + 2.9237e-3 i
```

Consequently the local helper multiplies both by the same near-pole
`1 / (r - rstar)` structure and the two threshold CTs add instead of
dual-cancelling.

The relevant E-surface equations at the engineered point are:

```text
T4  / group 196 / GL05_isr_p2_ct:
  -Q5 - Q6 + OSE12 + OSE16 + OSE9 + Q0 + Q1 + Q3 = 0

T34 / group 192 / GL05_isr_p1_ct:
  -Q4 - Q7 + OSE10 + OSE15 + OSE9 + Q0 + Q1 + Q2 = 0
```

The overlap machinery sees both group E-surfaces, but the actual CT evaluation
maps a group E-surface through the current graph position:

```text
GL05_isr_p2_ct:
  group 196 -> local 4, present and on shell
  group 192 -> absent in current graph

GL05_isr_p1_ct:
  group 192 -> local 34, present and on shell
  group 196 -> absent in current graph
```

There are same-edge-pattern local surfaces in the opposite graphs, but they are
not the dangerous partner at the intersection:

```text
GL05_isr_p2_ct group 97,  edges [9,10,15], eta(rstar) ~= -2.0316e3
GL05_isr_p1_ct group 95,  edges [9,12,16], eta(rstar) ~= -2.0316e3
```

The code path responsible for this is:

```text
AmplitudeGraph::preprocess
  -> build_threshold_counterterm_parametric_integrand
     builds threshold CT atoms from one graph's CFF expression only

AmplitudeCountertermData::evaluate
  -> OverlapBuilder::new_esurface_builder
     skips group E-surfaces absent in the current graph position

AmplitudeGraphTerm::evaluate_impl
  evaluates original graph - graph-local threshold CTs

generic evaluate_graph_group
  sums the already-subtracted graph terms
```

This is why a local sign flip is not the right fix. In the toy cancellation,
the relative minus comes from evaluating the other denominator in the selected
residue:

```text
R_A ~ N(r*) / (g_B (r_A* - r_B*))
R_B ~ N(r*) / (g_A (r_B* - r_A*))
```

The current graph-local CFF branch selection never produces the opposite
graph's denominator in the pass-one residue, so the source of that relative
sign is absent before the helper is evaluated.

The robust fix should move amplitude threshold-subtraction generation to a
graph-group object, or otherwise construct a group-level rational residue in
group E-surface variables before compiling the pass-one CT evaluators. Merely
delaying subtraction until after `evaluate_graph_group` is not sufficient,
because `(G1 - CT1) + (G2 - CT2)` is algebraically identical to
`(G1 + G2) - (CT1 + CT2)` if the CTs themselves remain graph-local.

The raw selected pass-one CFF factors confirm this diagnosis. Extracting the
linear inverse factors from the selected residues gives, for
`GL05_isr_p1_ct`/T34, factors such as:

```text
-Q4 - Q6 - Q7 + OSE15 + OSE16 + Q0 + Q1 + Q2
-Q4 - Q7 + OSE11 + OSE15 + Q0 + Q1 + Q2
-Q4 + OSE10 + OSE12
-Q4 + OSE10 + OSE13 + OSE9 + Q0 + Q1 + Q2
-Q7 + OSE12 + OSE15 + OSE9 + Q0 + Q1 + Q2
-Q7 + OSE13 + OSE15
-Q7 + OSE14 + OSE15 + OSE9 + Q1
```

and, for `GL05_isr_p2_ct`/T4:

```text
-Q5 - Q6 - Q7 + OSE15 + OSE16 + Q0 + Q1 + Q3
-Q5 - Q6 + OSE13 + OSE16 + Q0 + Q1 + Q3
-Q5 + OSE10 + OSE12
-Q5 + OSE11 + OSE12 + OSE9 + Q0 + Q1 + Q3
-Q6 + OSE10 + OSE16 + OSE9 + Q0 + Q1 + Q3
-Q6 + OSE11 + OSE16
-Q6 + OSE14 + OSE16 + OSE9 + Q0
```

The selected T34 residue no longer contains a T34 inverse factor by
construction, but it also contains no T4-like factor

```text
-Q5 - Q6 + OSE12 + OSE16 + OSE9 + Q0 + Q1 + Q3
```

Likewise, the selected T4 residue contains no T34-like factor

```text
-Q4 - Q7 + OSE10 + OSE15 + OSE9 + Q0 + Q1 + Q2
```

Thus the toy-model relative sign

```text
1 / (r_A* - r_B*)  versus  1 / (r_B* - r_A*)
```

cannot arise in the current selected residues. The overlap/rstar solver knows
both group E-surfaces exist and uses the same center, but the compiled
pass-one numerator for each graph has already been formed from a graph-local
CFF expression in which the other group surface is absent.

Additional implementation checks:

- `AmplitudeCountertermData::evaluate` applies each selected E-surface CT
  independently and sums the results. The final graph term uses
  `original - sum_of_cts`, so the final sign is common to all local threshold
  CTs and cannot generate the missing relative sign.
- The amplitude overlap prefactor is not the culprit. It is built from the
  complement of each overlap group and only partitions different overlap
  groups. For the dangerous point all relevant surfaces are in the same
  overlap group, and the raw dumps show the prefactor is `1` for the suspect
  CTs.
- The cross-section LU threshold path has the missing structural ingredient:
  it explicitly generates left, right, and left/right cartesian-product
  threshold residues, then combines them with inclusion-exclusion in
  `LUCounterTerm::evaluate`. The amplitude threshold path has no analogous
  same-overlap or graph-group iterated threshold residue.

Consequently, the next real fix should not be another local sign flip in
`threshold_counterterm_helper_atom`. The required change is architectural:
amplitude threshold subtraction needs a grouped/overlap-aware pass-one
construction that can form the residue of the graph-group object in common
group E-surface variables, or an equivalent same-overlap inclusion-exclusion
term, before the local radial helper is applied.

### Current status after external-slot fix and retained selector

Fresh validation was run on 2026-06-19 after rebuilding the CLI, cleaning the
state, and regenerating `qq_hhh_2L` with the non-compiled SingleParametric
generation block:

```text
./run_with_memory_watch.py --limit-gb 30 --interval 2 -- just build-cli
./run_with_memory_watch.py --limit-gb 30 --interval 2 -- \
  ./target/dev-optim/gammaloop --clean-state examples/cli/qq_hhh_2L/qq_hhh_2L.toml quit
./run_with_memory_watch.py --limit-gb 30 --interval 2 -- \
  ./target/dev-optim/gammaloop -o -s ./examples/cli/qq_hhh_2L/state \
  run generate_diagrams generate_integrands_single_parametric
```

The retained source-level settings/fixes are:

```text
crates/gammalooprs/src/cff/esurface.rs:
  Esurface::to_atom_in_lmb uses external slot ids from the LMB edge signature.

crates/gammalooprs/src/cff/hsurface.rs:
  Hsurface::to_atom_in_lmb uses the same external slot convention.

crates/gammalooprs/src/processes/amplitude.rs:
  amplitude threshold CT generation sets
  threshold_residue_in_generated_basis = false.
```

The focused regression tests and package check passed:

```text
./run_with_memory_watch.py --limit-gb 30 --interval 2 -- \
  cargo test -p gammalooprs to_atom_in_lmb_uses_external_slots_not_carrier_edges
  -> 2 passed

./run_with_memory_watch.py --limit-gb 30 --interval 2 -- cargo check -p gammalooprs
  -> passed
```

The engineered root probe was rerun with:

```text
x_root = [0.5387589249438893, 0.24174091659991892,
          0.9441181235200773, 0.7527842995682404,
          0.5585112730395129, 0.25621170988207465]
graph group = 0
LMB channel = 9
```

Fresh raw pass-one threshold CT residues at this root are now small and finite:

```text
GL05 T5:
  raw_result = +7.22605196011759e-7 + 2.12240111969347e-6 i

GL05 T31:
  raw_result = -2.57103769264393e-7 + 2.59605085238467e-6 i

GL05_isr_p2_ct T4:
  raw_result = +8.72856100696575e-8 - 4.72733492950589e-6 i

GL05_isr_p1_ct T34:
  raw_result = +1.33108172701426e-6 - 5.46209673633603e-6 i
```

This is the important comparison against the diagnostic generated-basis
selector:

```text
Using threshold_residue_in_generated_basis = true at the same engineered root:
  GL05 T5          ~= -2.41e-4 + 6.85e-4 i
  GL05 T31         ~= -3.07e-4 + 7.30e-4 i
  GL05_isr_p2 T4   ~= -2.64e-4 + 4.48e-3 i
  GL05_isr_p1 T34  ~= -1.10e-5 + 3.76e-3 i

So the generated-basis selector reintroduces inflated pass-one residues and
should not be used for this debugging target.
```

Latest watchdog probe at the same engineered root after rebuilding both CLI and
Python API:

```text
state:
  examples/cli/qq_hhh_2L/state_watchdog_probe
command block:
  inspect_threshold_overlap_root

CLI inspect:
  integrand_result          = -4.2848606320878566e-16 +1.0478274563803044e-15 i
  parameterization_jacobian = +5.4229032752001450e19
  integrator_weight         = +1.0
  threshold_counterterm:10:0 = -1.7514504877544214e-16 +7.4837130458710868e-17 i
  threshold_counterterm:33:0 = -2.7407962379201643e-16 -2.1202011538495107e-17 i

Python API evaluate_sample with the same state, runtime settings, and
discrete_dim=[0,9]:
  integrand_result          = -4.2848606320878566e-16 +1.0478274563803044e-15 i
  parameterization_jacobian = +5.422903275200145e19
  threshold_counterterm:10:0 abs = 1.904636033625483e-16
  threshold_counterterm:33:0 abs = 2.7489846393032416e-16

This confirms that the CLI and Python API are using the same fixed threshold
counterterm code path.
```

Scoped replacement validation after the UV regression check:

```text
The LMB external-slot replacement is now scoped to amplitude threshold CT
cutsets only:

  CutSet::canonicalize_external_shifts = true
    only in AmplitudeGraph::build_threshold_counterterm_parametric_integrand

  CutSet::canonicalize_external_shifts = false
    for ordinary UV/CFF and cross-section cutsets

Reason:
  globally replacing every CFF surface in LMB-slot convention leaked
  Q(0,cind(0)) into UV/local evaluator expressions whose parameter maps still
  use the edge-indexed external-energy symbols Q(4,cind(0)), Q(5,cind(0)), ...

Validation:
  cargo nextest run -p gammalooprs --cargo-profile dev-optim \
    -P test_gammaloop --run-ignored all \
    -E 'test(uv::tests::scalars_profile)' --no-fail-fast
    -> 2 tests run, 2 passed

  cargo test -p gammalooprs to_atom_in_lmb_uses_external_slots_not_carrier_edges
    -> 2 passed

  cargo check -p gammalooprs
    -> passed
```

Fresh root probe after rebuilding both CLI and Python API and regenerating
`examples/cli/qq_hhh_2L/state_clean_probe` with SingleParametric output:

```text
inspect_threshold_overlap_root:
  integrand_result          = -4.2848606419508805e-16 +1.0478274571169974e-15 i
  parameterization_jacobian = +5.4229032752001450e19
  integrator_weight         = +1.0
  generated events          = 6
  accepted events           = 6
  arb stability             = Stable(2 samples), rel. accuracy 7.5290827595602410e-26

  threshold_counterterm:10:0 =
    -1.7514504877544214e-16 +7.4837130458710868e-17 i
  threshold_counterterm:33:0 =
    -2.7407962379201643e-16 -2.1202011538495107e-17 i
```

Fresh 4096-point imaginary smoke integration after the scoped replacement:

```text
command block:
  integrate_debug_imag_4096

state:
  examples/cli/qq_hhh_2L/state_clean_probe

result:
  All re = -3.2(9.1)e-3   285%
  All im = -1.19(85)e-2   70.9%

global max weights:
  re[+] = +2.1768988663429763e1
    graph 0, LMB channel 11,
    x = [0.23844051619928841, 0.41480916417720726,
         0.84657084021138962, 0.33683279500664520,
         0.51358515737196997, 0.078273790141412225]

  re[-] = -1.3851544772165392e1
    graph 0, LMB channel 5,
    x = [0.19383529716105175, 0.31833693697749332,
         0.23544541436522859, 0.45637354424474219,
         0.53736706089196995, 0.44996245133302959]

  im[+] = +1.2478804222096382e1
    graph 1, LMB channel 11,
    x = [0.71906345171460317, 0.26040351173928078,
         0.80739917696325181, 0.93515792609772008,
         0.57309382535847941, 0.76652808410758222]

  im[-] = -1.5177455697466318e1
    graph 1, LMB channel 1,
    x = [0.61777759691350409, 0.30698347744194177,
         0.89614285895604617, 0.44171512552571912,
         0.47700201088427285, 0.18926861664962313]

statistics:
  f64             = 100.00%
  f128            = 0.00%
  arb             = 0.00%
  nans+unstable   = 0.00%
  peak RSS        = 2.91 GiB under run_with_memory_watch.py --limit-gb 30
```

Fresh fixed-bin event decomposition at the low-stat imaginary max points:

```text
global_im_plus_4096:
  weighted_abs = 3.76020481e-1, nan = false
  dominant event = GL12_isr_p2_ct, |event| = 3.76167550e-1
  largest CT = threshold_counterterm:48:0,
    -2.299715e-1 +2.968640e-1 i, |.| = 3.75519770e-1

global_im_minus_4096:
  weighted_abs = 4.44734039e-1, nan = false
  dominant event = GL12_isr_p2_ct, |event| = 4.43617356e-1
  largest CT = threshold_counterterm:45:0 edges [9,13,16],
    -2.761941e-1 -3.750059e-1 i, |.| = 4.65738722e-1

group0_im_plus_4096:
  weighted_abs = 1.80475169e-1, nan = false
  dominant event = GL06, |event| = 1.76244869e-1
  largest CT = threshold_counterterm:16:0 edges [9,11,12],
    -8.205254e-2 +8.723207e-2 i, |.| = 1.19758313e-1

group0_im_minus_4096:
  weighted_abs = 3.18161921e-2, nan = false
  dominant event = GL06_isr_p1_ct, |event| = 1.29901248e-1
  largest CT = threshold_counterterm:17:0 edges [9,13,14],
    +9.341700e-2 +8.139351e-2 i, |.| = 1.23901730e-1

Conclusion:
  these fixed-bin probes do not show the previous runaway missed-dual-
  cancellation pattern. The individual CT pieces are O(1), finite, and in the
  same scale as the corresponding event weights.
```

Graph-by-graph dual-cancellation conclusion:

```text
The original GL05 graph-local dangerous pair is fixed:

  T5  group 140, edges [9,12,16]
  T31 group 143, edges [9,10,15]

At the engineered root both are present in GL05 and have the same rstar and
deta/dr, so the graph-local residue construction has the required ingredients
for the usual dual-cancellation mechanism.

The ISR T4/T34 pair is different:

  GL05_isr_p2_ct T4   group 196, edges [9,12,16]
  GL05_isr_p1_ct T34  group 192, edges [9,10,15]

At the same root:
  group 196 is present in GL05_isr_p2_ct and absent in GL05_isr_p1_ct.
  group 192 is present in GL05_isr_p1_ct and absent in GL05_isr_p2_ct.

The common graph-group overlap center is therefore doing what it can, but the
usual toy-model sign flip from evaluating the "other" denominator in the same
selected graph residue cannot occur across these two ISR graph terms. Any
remaining large finite weights from this structure should be studied as a
graph-group threshold-subtraction issue, not as a missing graph-local
dual-cancellation partner inside GL05 itself.
```

### Graph-local overlap localization update

A later code pass made the overlap treatment match the latest diagnosis more
closely:

```text
Common graph-group behavior retained:
  - group-level E-surface existence is true if any graph-local incarnation
    exists;
  - center solving still uses all locally existing incarnations in the group;
  - this keeps one common center across the graph group, which is required for
    comparable radial variables.

Graph-local behavior added:
  - each graph term records local_esurface_exists[group_esurface_id];
  - threshold CT evaluation skips group surfaces that do not exist in that
    graph term;
  - the group overlap structure is localized before building that graph term's
    overlap prefactor evaluators;
  - complements are recomputed after localization, so graph-local
    dual-cancellation patterns are no longer polluted by surfaces that exist
    only in sibling graph terms.
```

Relevant source locations:

```text
crates/gammalooprs/src/integrands/process/amplitude/mod.rs
  get_existing_esurfaces now populates per-graph local_esurface_exists.
  overlap.localized_to_existing_surfaces(...) is applied before assigning the
  threshold_counterterm.overlap for each graph term.

crates/gammalooprs/src/subtraction/overlap.rs
  OverlapInput includes local_esurface_exists.
  center checks require every locally existing incarnation of a group surface
  to be inside the shared center.
  OverlapStructure::localized_to_existing_surfaces remaps the common group
  overlap to the graph-local existing-surface set and recomputes complements.

crates/gammalooprs/src/subtraction/amplitude_counterterm.rs
  evaluate skips selected E-surfaces whose group surface does not exist in the
  current graph term.
```

Focused validation after this localization:

```text
cargo fmt --check
  -> passed

./run_with_memory_watch.py --limit-gb 30 --interval 2 -- cargo check -p gammalooprs
  -> passed, peak RSS ~= 0.98 GiB

./run_with_memory_watch.py --limit-gb 30 --interval 2 -- \
  cargo test -p gammalooprs overlap_structure_localizes_to_graph_existing_surfaces
  -> passed

./run_with_memory_watch.py --limit-gb 30 --interval 2 -- \
  cargo test -p gammalooprs to_atom_in_lmb_uses_external_slots_not_carrier_edges
  -> 2 passed

./run_with_memory_watch.py --limit-gb 30 --interval 2 -- just build-cli
  -> passed
```

The engineered root was reprobed after cleaning the state and regenerating with
the non-compiled SingleParametric block. For the selected dangerous surfaces:

```text
x_root = [0.5387589249438893, 0.24174091659991892,
          0.9441181235200773, 0.7527842995682404,
          0.5585112730395129, 0.25621170988207465]
graph group = 0
LMB channel = 9

GL05 T5:
  raw_result ~= +7.23e-7 +2.12e-6 i
  prefactor  = 1

GL05 T31:
  raw_result ~= -2.57e-7 +2.60e-6 i
  prefactor  = 1

GL05_isr_p2_ct T4:
  raw_result ~= +8.73e-8 -4.73e-6 i
  prefactor  = 1

GL05_isr_p1_ct T34:
  raw_result ~= +1.33e-6 -5.46e-6 i
  prefactor  = 1
```

The same probe reports a common radial geometry for the T5/T31/T4/T34 root
structures:

```text
radius      ~= 1.0618603586916703e3
radius_star ~= 1.0618603586915431e3
deta/dr     ~= 2.3502742820729676
```

The center dumps also confirm that localization is active while the center is
shared:

```text
GL05 selected surfaces:
  localized overlap_group_size = 15

GL05_isr_p2_ct selected surfaces:
  localized overlap_group_size = 14

GL05_isr_p1_ct selected surfaces:
  localized overlap_group_size = 14
```

A short post-localization SingleParametric t=4 sanity integration with
`sliver_width=30`, threshold `h_sigma=1`, `gaussian_width=1`, and 4096 points
was finite:

```text
re = +1.296154082749754e-2 +/- 2.1563212648705775e-2
im = -5.585394162946225e-3 +/- 1.1526410120066537e-2

max_re_plus  = +8.016765903727686e1
max_re_minus = -3.372170192412362e1
max_im_plus  = +1.7304666135846336e1
max_im_minus = -4.238957848210082e1

f128 = 0.0244%, arb = 0.0244%, nans+unstable = 0%
```

This is not a precision run, but it is a useful regression signal: after
graph-local overlap localization, the quick scan does not reproduce the old
catastrophic threshold CT spikes or NaNs. A production compiled/SummedFunctionMap
state should still be regenerated cleanly before a serious run; one compiled
generation attempt in this debugging session was interrupted while compiling,
so do not treat the current `state` directory as a confirmed fresh production
artifact unless `run generate_diagrams generate_integrands` has completed
normally in the same checkout.
