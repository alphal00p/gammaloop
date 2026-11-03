#show table: set text(size: 8pt)

= Causal LU centers for thermal two-point functions

This note records the causal-center search for the spatial Local Unitarity flow, its
validation on the massive scalar double triangle (DT) and triangle-box-triangle (TBT), and
the remaining limitations. Both graphs use the
#link("../massive_double_triangle/benchmark.typ")[double-triangle model and normalization]:
`mH=0.6`, `mL=0.25`, `T=0.5`, all cubic couplings one, graph multiplier `1/2`, and the
observable $Gamma_rho=-op("Im") Pi_R(P)/E$ at probe momentum $P=(E,0,0,p)$. The TBT
topology follows `tests/resources/graphs/triangle_box_triangle_phys.dot` with the scalar
assignment in #link("triangle_box_triangle.dot")[`triangle_box_triangle.dot`].

== What a center has to do

LU dilates the spatial loop momenta about a fixed point, $bold(y)=bold(c)+t(bold(x)-bold(c))$,
with Jacobian $t^(3L)$, and localizes every cut $eta_i(bold(y))=0$ at its radial roots. Two
cuts that cross at an overlap cancel locally only if their slopes
$nabla eta_i dot (bold(y)-bold(c))$ share a sign there; otherwise the pair keeps a
$1\/Delta$ pole and the Monte Carlo variance is infinite. A ray tangent to a cut surface (a
fold) gives an integrable $abs(Delta t)^(-1\/2)$ weight for a simple cut, whose variance
still diverges logarithmically, and a nonintegrable one for raised cuts. Each physical cut
equation is oriented as emitted minus absorbed energies minus $E$, so positive slopes are
causal. E-surfaces are convex and positive about every interior center; thermal H-surfaces
are not.

The proposed argument holds where it applies. At a finite smooth point a common causal
direction exists unless the active gradients have a nonzero nonnegative linear dependence
(Gordan's alternative, as in Collins' Theorem 3). The obstruction set is semialgebraic, and
along each smooth stratum the dependence forces $dif omega=0$, so only finitely many external
energies are exceptional. It says nothing at infinity, at massless kinks, or about radial
centers. Three facts decide the radial problem instead.

+ For one surface point the valid centers form an open half-space, so the centers valid on a
  sampled set of cut points solve a linear program, and the exact valid set is convex.
+ The complete set of cut points is invariant under the rotations that fix the external
  momenta. Averaging a valid center over them keeps it valid, so a single center can be
  sought on the probe axis, and at rest the origin is the only candidate. Nonzero centers at
  rest are useful only as patches.
+ Thermal H-surfaces with one absorbed particle have unbounded lobes around collinear
  configurations, and the TBT at rest has rotation-invariant H-surfaces that do not contain
  the origin. No finite set of radial centers sees all of them star-shaped. At spacelike
  momenta the two sheets of a two-particle absorption pair are oriented by a common center
  only at their hyperboloid center $bold(p)\/2$, where the margin decays at infinity.

== Implementation

`lu_flow.centers = "search"` prepares a plan at warm-up in
#link("../../../../crates/gammalooprs/src/integrands/process/cross_section/lu_flow.rs")[`cross_section/lu_flow.rs`].
It shoots random rays through the origin and through random base points and samples every
active cut surface with the certified signed radial root solver. It keeps points whose
absorbed bath energies are within the `relevance` Boltzmann factor of the most relevant point.
Every point, and every other cut surface within `slope` of it, gives a causal half-space with
the slope threshold `slope` times the energy scale, decaying quadratically beyond that scale.
A weighted L1 relaxation of these half-spaces (Clarabel) gives the first center in the
invariant subspace, and a pattern search on the covered Boltzmann weight refines it. Further
centers cover the remaining points while each adds at least `min_patch_gain`. A second sampling
round shoots rays through the first centers, and an independent sample verifies the plan.
The warm-up log reports the causal, partial and uncovered shares and the overlap points that
keep opposite slopes.

Several centers split the integrand with a smooth partition of unity in the physical point.
The weight of patch $k$ is $phi_k\/(sum_l phi_l)$, where $phi_k$ is the product over cut
surfaces of an infinitely differentiable gate. The gate closes where a surface is close to the
point and center $k$ is not causal on it. A small floor on the first patch keeps the partition
exact everywhere. Each patch solves the roots of all cut groups about its own center, and every
root of every cut carries the same weight. The weight multiplies the `h` profile and its jets
along the patch ray, so raised residues differentiate it. Momenta, thermal occupations,
events and residue jets use the patch map, whose Jacobian remains $t^(3L)$. Centers are
binary64 points in the warm-up frame, embedded exactly in every precision and rotated with the
stability rotations. Positive-E sampling hosts belong to the first patch. Threshold subtraction
and zero-temperature Fermi sectors still require origin flow.

`integrator.training_iterations` freezes the sampling grid after that many iterations. The
#link("run.toml")[run card] generates both graphs, and #link("campaign.py")[`campaign.py`] runs
independent seeds, each adapting for 10 iterations of 10,000 samples and then sampling the
frozen grid. #link("analyze.py")[`analyze.py`] reconstructs frozen-sample sums from the
per-iteration archives. It compares the spread of disjoint block estimates with their internal
errors, fits both against $N^b$, and bootstraps the slope over seeds.

== Independent references

The #link("reference/derivation.typ")[reference engine] performs the Matsubara sums by iterated
residues at fixed spatial momenta, continues to $E+i epsilon$, integrates over spatial momenta,
and extrapolates $epsilon -> 0$. It uses no GammaLoop expression. It reproduces the published
DT targets at rest (pulls between $-1.05$ and $1.22$, including $T=0.25$ and vacuum), analytic
scalar bubbles at rest, timelike and spacelike momenta, and vacuum Lorentz invariance. It then
supplies the spacelike DT and all TBT targets. An independent review reran points with new
seeds and $epsilon$ ladders and checked the TBT vacuum against a Feynman-parameter calculation
(0.1%). The TBT $(1.6,0)$ value combines both runs with a PDG scale factor, and the TBT $(1,0)$
error is inflated 1.5 times for its low final $chi^2$ probability.
#link("reference.json")[`reference.json`] holds the values used by `analyze.py`, and
#link("reference/results.json")[`reference/results.json`] the runs behind them.

#table(
  columns: 4,
  table.header([Graph, $P$], [Reference], [Uncertainty], [Source]),
  [DT $(1,0)$], [`3.057651940618655e-4`], [quadrature], [Blaizot–Reinosa Eq. (58)],
  [DT $(1.6,0)$], [`5.282087736915636e-6`], [quadrature], [Blaizot–Reinosa, Bauberger–Böhm],
  [DT $(0.5,1.5)$], [`6.67521e-5`], [`1.3e-7`], [reference engine],
  [TBT $(1,0)$], [`5.22013e-6`], [`1.4e-7`], [reference engine],
  [TBT $(1.6,0)$], [`-1.91543e-6`], [`5.7e-9`], [reference engine, two runs],
  [TBT $(0.5,1.5)$], [`1.13479e-6`], [`2.7e-8`], [reference engine, two runs],
)

== Validation

Unit tests check the gate limits, that patch weights sum to one at a common physical point,
that their dual jets agree with central differences, and that the search orients the scalar
vertex overlap where origin flow fails. The integration test probes the DT overlap of the LL
and HH two-particle absorption cuts at $P=(0.5,0,0,1.5)$, off the probe axis, at relative
distances $plus.minus 10^(-3)$ and $plus.minus 10^(-6)$. It requires the origin cut sum to grow
by more than 300 and to change sign across the overlap, and the searched-center sums on both
sides to agree within $10^(-2)$ of their scale. The `approach_dt_overlap` block of the
#link("run.toml")[run card] scans the same overlap with the `approach` command, moving the LL
root by a relative distance $t$ from $10^(-2)$ down to $10^(-10)$ on both sides and recording every cut's
event weight. With origin flow the LL and HH cuts both grow like $1.15 times 10^(-7)\/t$ with
the same sign, so their sum is a $1\/t$ pole that reaches $plus.minus 2.3 times 10^3$. With
the searched center and with the apex center the two cuts carry opposite poles,
$plus.minus 1.07 times 10^(-7)\/t$ and $plus.minus 1.56 times 10^(-7)\/t$, and their sums stay
at `1.308635e-5` and `1.810939e-5` from $10^(-4)$ to $10^(-10)$. Pointwise values depend on
the flow; only the integral is center independent. Two numerical caveats appear at the
extremes. With forced arbitrary precision (`approach -f`), the searched plan returns NaN at
$t=plus.minus 10^(-4)$. A patch-1 ray meets an H-surface lobe at radius 134, where the
absolute residual tolerance cannot be met without a lower-precision root to compare with.
The default precision ladder evaluates these points correctly. At $t=-10^(-10)$, the apex
flow accepts a binary64 value of `2.68e-4` instead of `1.81e-5`: both stability rotations
reproduce the same root-tolerance error, so the rotation check does not escalate.

#figure(
  image("results/dt_overlap_approach.png", width: 100%),
  caption: [Complete DT cut sum approaching the LL/HH overlap at $P=(0.5,0,0,1.5)$ with
    origin, searched and apex centers (`approach_dt_overlap`, default precision ladder);
    dashed lines are negative values.],
)

#table(
  columns: 4,
  table.header([Configuration], [Plan], [Causal / uncovered weight], [Overlap points with opposite slopes]),
  [DT $(1,0)$], [origin], [1.000 / 0.000], [0 of 4317],
  [DT $(1.6,0)$], [origin], [1.000 / 0.000], [0 of 5854],
  [DT $(0.5,1.5)$], [$(0.708, 0.775) hat(bold(p))$ and 3 patches], [0.949 / 0.046], [237 of 1908],
  [TBT $(1,0)$], [origin and 1 patch], [0.964 / 0.025], [261 of 21616],
  [TBT $(1.6,0)$], [origin], [1.000 / 0.000], [0 of 23064],
  [TBT $(0.5,1.5)$], [$(0.731, 0.788, 0.894) hat(bold(p))$ and 5 patches], [0.813 / 0.175], [3653 of 20190],
)

The table lists the plans of the final build for the default `lu_flow.search` settings. The
first center of a plan is written in units of $bold(p)$ per loop momentum of the generation
basis.

Every configuration uses eight independent seeds. The frozen samples per seed are $10^6$ for the
DT and $4 times 10^5$ for the TBT. The ratio compares the spread of disjoint block estimates
with their median internal error at the smallest and largest block. A ratio near one with
slope $-1\/2$ is a finite-variance $N^(-1\/2)$ estimate. A large ratio means rare blocks carry
outliers that internal errors miss, which is a heavy-tailed or infinite-variance weight
distribution. In that regime the block-spread slope alone can still look like $-1\/2$.

#table(
  columns: 7,
  table.header(
    [Graph, $P$], [Center], [Weighted mean (pull)], [Seed mean (pull)], [Spread / error], [Slope (bootstrap)], [ms/sample],
  ),
  [DT $(1,0)$], [origin = search], [`3.0588e−4 ± 1.4e−7` (+0.84)], [`3.0588e−4 ± 2.1e−7` (+0.56)], [1.02 – 1.16], [−0.49 (−0.53, −0.46)], [0.87],
  [DT $(1.6,0)$], [origin = search], [`5.34e−6 ± 1.3e−7` (+0.40)], [`5.34e−6 ± 1.6e−7` (+0.33)], [1.00 – 1.25], [−0.48 (−0.51, −0.46)], [1.19],
  [DT $(0.5,1.5)$], [origin#super[\*]], [`6.58e−5 ± 1.0e−6` (−0.94)], [`5.7e−5 ± 1.2e−5` (−0.84)], [84 – 14], [heavy tail], [1.04],
  [DT $(0.5,1.5)$], [L1 center, unrefined#super[\*]], [`6.703e−5 ± 4.5e−7` (+0.59)], [`7.47e−5 ± 6.5e−6` (+1.22)], [36 – 19], [heavy tail], [0.83],
  [DT $(0.5,1.5)$], [apex, fixed#super[\*]], [`6.688e−5 ± 1.9e−7` (+0.59)], [`6.689e−5 ± 2.3e−7` (+0.53)], [1.59 – 1.34], [−0.44 (−0.47, −0.43)], [1.23],
  [DT $(0.5,1.5)$], [search, 4 patches], [`6.690e−5 ± 1.0e−7` (+0.93)], [`6.6907e−5 ± 6.6e−8` (+1.09)], [1.24 – 1.00], [−0.52 (−0.54, −0.51)], [3.05],
  [TBT $(1,0)$], [origin#super[\*]], [`5.17e−6 ± 3.5e−7` (−0.13)], [`7.0e−6 ± 1.3e−6` (+1.38)], [16 – 6.1], [heavy tail], [6.44],
  [TBT $(1,0)$], [search, 2 patches], [`4.78e−6 ± 2.5e−7` (−1.54)], [`5.24e−6 ± 3.4e−7` (+0.06)], [3.24 – 1.23], [−0.53 (−0.56, −0.52)], [13.7],
  [TBT $(1.6,0)$], [origin = search], [`−1.952e−6 ± 2.4e−8` (−1.47)], [`−1.942e−6 ± 4.2e−8` (−0.62)], [1.10 – 1.43], [−0.38 (−0.44, −0.35)], [5.65],
  [TBT $(0.5,1.5)$], [origin#super[\*]], [`1.30e−6 ± 4.2e−7` (+0.40)], [`−1.9e−5 ± 1.2e−5` (−1.58)], [119 – 35], [heavy tail], [7.73],
  [TBT $(0.5,1.5)$], [L1 centers, 2 patches, unrefined#super[\*]], [`7.1e−7 ± 6.8e−7` (−0.62)], [`−9e−7 ± 1.7e−6` (−1.18)], [15 – 3.8], [heavy tail], [9.87],
  [TBT $(0.5,1.5)$], [apex, fixed#super[\*]], [`1.11e−6 ± 8.4e−7` (−0.03)], [`−1.0e−3 ± 1.0e−3` (−1.00)], [9601 – 2181], [heavy tail], [5.37],
  [TBT $(0.5,1.5)$], [search, 6 patches], [not completed], [], [], [], [172],
)

#super[\*] Earlier builds of this branch. Their origin and fixed-center paths are unchanged:
the first ten iterations of the spacelike DT origin and apex runs repeat bit for bit with the
final build, and the final TBT $(1.6,0)$ campaign repeats the earlier one bit for bit in every
seed. The unrefined L1 rows document the superseded search without the coverage refinement.
Costs are per sample
and core on a host running at a load of about 420 on 384 cores, so compare them only within
the table. Every `campaign.py` output directory records the binary and state hashes and the
tree hash of the working directory.

The six-patch TBT plan at $(0.5,1.5)$ was not run to completion. Every patch solves the roots of
all thermal cut groups. Rays that graze the uncovered H-surfaces escalate 36% of the
evaluations to quadruple and 1.2% to arbitrary precision, so a sample costs 22 times the origin
flow. The plan still leaves 18% of its overlap points misoriented. Moreover, 0.02–0.06% of its
samples stay nonfinite even in arbitrary precision and are dropped, from the far-root tolerance
described above. Every other configuration reports 0.00% nonfinite or unstable final samples.

Pulls include the reference uncertainty. The weighted mean combines the seeds with
inverse-variance weights; the seed mean is their plain average, with its error from their
scatter. With finite variance the two agree. "Heavy tail" marks configurations whose block
spread exceeds the internal errors by at least five times. Their fitted block-spread slope is
near $-1\/2$, since a few outliers dominate every block size, so it is no evidence of finite
variance. In those configurations a seed that hits a large weight also reports a large error,
so inverse-variance weighting suppresses exactly the outliers that carry the tail. The
weighted means then agree with the references while the seed means scatter far, down to
$-1.0 times 10^(-3)$ for the spacelike TBT apex. The two-patch TBT at $(1,0)$ shows the same
effect mildly: its weighted mean is 1.5 standard deviations low, while its seed mean agrees.

#figure(
  grid(columns: 2, gutter: 4pt,
    image("results/dt_space_origin.png"), image("results/dt_space_search.png"),
    image("results/tbt_rest1_origin.png"), image("results/tbt_rest1_search.png"),
  ),
  caption: [Block spread and internal errors against samples per block, and cumulative
    estimates of every seed against the independent reference, for origin flow and the searched
    plan at the spacelike DT and the TBT at rest. Every configuration of the table has its figure
    in `results/`.],
)

== Investigation of large weights

The largest spacelike DT weights with origin flow sit at the LL/HH overlap, where cut groups
13 and 16 have roots $2 times 10^(-6)$ apart with slopes $+1.67$ and $-0.78$. Their
contributions add instead of cancelling. Three of eight seeds then carry weights up to 90
times the typical ones, and the estimates scatter far beyond their errors. The apex center
$(bold(p)\/2, bold(p)\/2)$ orients both two-particle pairs; its eight seeds agree with each
other and with the reference. The first searched center, before the coverage refinement, put
the HH center beyond the apex and produced a fold of the HH sheet, with two roots of slopes
$plus.minus 6 times 10^(-4)$. One of eight seeds had a weight 52 times the typical one. The L1
relaxation minimizes total violation rather than covered weight; refining the center on the
covered weight moves the first center to $(0.708, 0.775) hat(bold(p))$, next to the apex,
and three further patches raise the causal weight to 94.9%. Its eight seeds are consistent
and its block spread follows the internal errors.

At rest the TBT origin is misoriented on the rotation-invariant H-surfaces of the four-particle
cuts. At its largest weight an E-surface and such an H-surface have roots $1.2 times 10^(-5)$
apart with slopes $+0.26$ and $-0.039$, contributing $82$ and $41$ with the same sign. Since
the invariant single center is the origin, only asymmetric patches could help. The greedy
cover adds one patch that covers part of these surfaces. It lowers the block-spread ratio from
16 to 3.2 and brings the seed mean onto the reference, but 261 of 21616 sampled overlap points
still keep opposite slopes, and a sample costs twice as much.

== Reproduction

The card generates both graphs into its state folder (`/tmp/gammaloop-thermal-lu-centers-state`
unless `-s` overrides it). The campaigns, analyses and overlap scan of this note then run on
that state in read-only mode. The two TBT kinematics sample 40 frozen iterations per seed
instead of 100.

```sh
gammaloop examples/cli/thermal_rates/lu_centers/run.toml run generate generate_tbt save
python3 examples/cli/thermal_rates/lu_centers/campaign.py --state <state> --output <dt> \
  --process dt --kinematics rest1 rest16 space --centers origin search apex --cores 4
python3 examples/cli/thermal_rates/lu_centers/campaign.py --state <state> --output <tbt> \
  --process tbt --centers origin search apex --frozen 40 --cores 4
python3 examples/cli/thermal_rates/lu_centers/analyze.py <dt>
gammaloop -n --read-only-state -s <state> examples/cli/thermal_rates/lu_centers/run.toml \
  run kinematics approach_dt_overlap -D energy=0.5 -D momentum=1.5 \
  -D e_cm=1.4142135623730951 -D centers=search -D output=<scan>.json
python3 assets/plot_approach_result.py <scan>.json --component real --x-log-scale \
  --branch-layout split --include-contribution '^total_weight$' \
  --include-contribution '^event_weight .* cut(0|2)$' --output <scan>.pdf
```

For the apex scan, add the `apex_dt` block after `kinematics` and set `centers=origin`.

== Supported cases and limitations

- Supported: fully cut massive thermal two-point functions with origin, fixed or searched
  centers, including smooth patches and raised cuts, with threshold subtraction disabled and
  without zero-temperature Fermi sectors.
- Validated with finite-variance $N^(-1\/2)$ scaling and reference agreement: the DT at rest
  below and above threshold, and the spacelike DT with the searched plan (block spread over
  internal error 1.24 to 1.00, slope $-0.52$) or the apex center (1.59 to 1.34). Origin flow
  fails there with an uncancelled overlap pole.
- Consistent with the reference, with a mild tail: the TBT at $(1.6,0)$. Every sampled cut
  point is causal for the origin, yet the block spread grows to 1.4 times the internal error
  and the slope is $-0.39$. The tail does not come from a misoriented overlap; its source was
  not isolated.
- Improved, not resolved: the TBT at $(1,0)$, below the HH threshold. Origin flow is heavy
  tailed (16 to 6) on the rotation-invariant H-surfaces. The two-patch plan reduces the ratio
  to 3.2 at the smallest and 1.2 at the largest block (slope $-0.53$), and its seed mean
  agrees with the reference (pull $+0.06$). The invariant center is the origin, so only
  asymmetric patches can help, and 261 of 21616 sampled overlap points keep opposite slopes.
- Not resolved: the spacelike TBT. Origin, apex, unrefined L1 and searched plans all leave
  misoriented H-surface overlaps, and the six-patch plan is too expensive to run. The
  collinear and rotation-invariant H-surfaces need a structured patch family. Only the
  inverse-variance-weighted means agree with the reference; the seed means do not.
- The search samples surfaces and certifies coverage only statistically. Its objective is the
  Boltzmann-weighted share of causal cut points, a proxy for the integrand's variance.
  Thermal suppression bounds the region considered; the slope threshold decays at large
  momenta.
- Massless particles, soft kinks and the light cone are outside this validation. The
  finiteness argument for exceptional energies does not cover them.
- Patches multiply the root work by their number; weights skip integrand evaluations where
  they vanish.
- Patch rays reach far roots on unbounded H-surface lobes. Without a lower-precision root to
  compare with, the radial solver cannot meet its absolute residual tolerance there. Forced
  arbitrary precision (`inspect -f`, `approach -f`) can then return NaN, and the six-patch
  spacelike TBT plan drops a few samples in $10^4$. Within $10^(-10)$ of an overlap, a binary64
  root-tolerance error can also pass the rotation stability check.
