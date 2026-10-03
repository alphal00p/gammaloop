= Exact-zero-temperature Fermi-surface sampling

#show table: it => block(breakable: false, it)

== Findings

*At exact zero temperature, these Fermi channels show no substantial time
advantage at the tested budgets.* In a matched Double-to-Arb comparison, the
sunrise Fermi pair is 0.30× as efficient as ordinary optimized LMB sampling;
uncrossed tennis ball is 1.07× with three shells and 1.04× with its repeated
shell; Mercedes is 1.00×. The three-loop gains are small compared with their
variation across seeds. All seven control means are within 1.53 reported
standard errors of the supplied references.

*The original Double/Quad/Arb ladder fails numerical validation with Fermi
sampling in all three graphs.* It accepts spurious large tail contributions
in Quad while reporting stability. Those failures are retained and reproduced
at fixed coordinates. Skipping Quad avoids the observed failures, and all
70 final control extrema agree with forced Arb within $4 times 10^(-10)$
relative difference. This is an experimental precision control, not a
production fix or a proof of accuracy for all possible samples.

== Scope and conventions

This repeats the #link("../report.pdf")[sunrise] and
#link("../three_loop/report.pdf")[three-loop] sampling measurements at exact
$T = 0$ on rebased revision `8ec568b6c736da6ccb7900569923ab3f6e1a0578`.
The two three-loop graphs are Mercedes (GL1) and uncrossed tennis ball (GL3).
The crossed diagram is not part of this comparison.

Every integrand is generated freshly with
`generation.medium.mode = "zero_temperature_equilibrium"` and vacuum
subtraction enabled. This is the exact-zero-temperature implementation;
`general.inverse_temperature` does not set its temperature. The driver records
`beta: null` and `zero_temperature: true` rather than using a large finite beta.
The old finite-temperature states and measurements are retained separately.

User-supplied reference imaginary integrals are:

#table(
  columns: 2,
  table.header([*Graph*], [*Exact $T = 0$ reference*]),
  [Sunrise], [−0.0161257672165997446],
  [Mercedes], [−0.030797946764053006],
  [Uncrossed tennis ball], [+0.011121480775908030],
)

The physical setup matches the preceding measurements: one massless quark,
`aS = 1`, `muB = 3`, hence quark chemical potential and Fermi radius both 1.
The sunrise uses `m_uv = 2`; the supplied three-loop setups use `m_uv = 1`.
All use `mu_r = renormalization_localization_scale = 2`, local and integrated UV
counterterms, disabled threshold subtraction, and summed physical orientations.
Assembly evaluators and the Double/Quad/Arb stability checks are retained.

Exact-zero-temperature distribution terms use the rebased implementation's
localization machinery. Its fixed positive-scale profile is
`h_function = { function = "poly_exponential", sigma = 1, power = 0 }`.
That physical distribution localizer is shared by all sampling modes. It is
distinct from the optional Fermi importance-sampling channel. A finite-beta
peak can become a localized contribution whose variance need not respond to the
same sampling map. Changes from the previous high-beta results also span a
code rebase, so they are not a controlled temperature-only comparison.

== Sampling comparison

Ordinary sampling uses `auto:optimized_lmb`. Augmented sampling retains all
ordinary channels and adds the product of the two independent sunrise Fermi
shells, or of all three independent fermion momenta `[0,1,3]` for either
three-loop graph. Tennis ball also tests `fermi_repeated`: only edge 0's shell,
with edges 1 and 3 as ordinary complement coordinates. Edge 5 shares edge 0's
momentum. This is the channel that improved the earlier finite-beta benchmark.

Both catalogues use Monte Carlo channel selection, density partitioning,
`power = 2`, and ordinary spherical linear mapping with `b = e_cm = 1`.
Each mode trains its own fresh adaptive grid. All other physics and stability
settings are shared. The first differing boundary is the channel catalogue
and its momentum map and partition weight; physical evaluation remains shared.

Efficiency is ordinary mean $sigma^2 t$ divided by candidate mean $sigma^2 t$,
using reported integration errors and wall time. Values above one favor the
candidate. This estimates equal-error time efficiency at the measured sample
budget, not a direct fixed-precision speedup. Means across seeds use equal
weights, with propagated error $sqrt(sum_i sigma_i^2)/n$. The observed spread
of per-seed efficiency ratios is reported separately and is not a confidence
interval. Agreement with the supplied integral is checked in addition to
reported variance.

Wall time includes state loading, warm-up, integration and persistence.
Generation and compilation are separate shared costs. CPU time, mapping and
physical-evaluation costs, precision escalation, channel populations, and
rejected samples are retained in the raw records. All integrations use one
worker pinned to a distinct CPU on the shared AMD EPYC 9754 host.
The replicated controls run five concurrent seed batches, with modes sequential
within each batch. This is a shared-host measurement, not an isolated-machine
benchmark; small timing differences warrant cautious interpretation.

== Sunrise: accuracy before efficiency

Five seeds (1337, 2843, 4517, 6007, 7919) each use 100,000 draws over ten
equal-size adaptive iterations. The initial comparison uses the original
Double/Quad/Arb ladder. Ordinary sampling gives
$-0.0159613 plus.minus 0.0001839$ for the five-seed mean, within 0.90 standard
errors of the reference, at 10.10 s per run. All five ordinary runs agree within
1.46 reported standard errors.

*The default-ladder Fermi result fails validation.* Seeds 4517 and 7919 give
−2,743,480.40 ± 66,756.56 and −12.9227 ± 0.4615, respectively, instead of
−0.0161258. The other three seeds are within two reported standard errors.
All runs nevertheless report zero rejected samples. Thus a clean rejection
counter is insufficient evidence of numerical correctness here. Every seed,
including both failures, is retained in `sunrise/results.json`; its numerical
efficiency summary is not a valid sampling-performance estimate.

The smallest reproducer replays two stored extreme points from seed 4517 in
the `fermi_pair` channel (graph/channel indices `0 3`). Both use the same
generated integrand, physical settings, subtraction, channel, coordinates and
partitioning. Only precision selection changes:

#table(
  columns: 3,
  table.header([*Point*], [*Default inspection*], [*Forced Arb inspection*]),
  [Positive extreme], [$9.2420876 times 10^15$], [$-4.1463525 times 10^(-254)$],
  [Negative extreme], [$-6.0108215 times 10^7$], [$2.0682123 times 10^(-262)$],
)

Default inspection accepts Quad as stable for both points, reporting relative
accuracy around $10^(-31)$. Arb instead returns negligible contributions.
These are directly inspected integrand weights, before the integration grid's
importance weights; they need not equal the recorded integration extrema.
This identical-coordinate precision test excludes ordinary Monte Carlo
fluctuation and reference-normalization differences as explanations. It
demonstrates a false finite-precision stability acceptance, without assigning
the defect to a particular algebraic term. No production implementation was
changed in this investigation.

Both points lie far in the Fermi map's outer tail: their first radial
coordinates are 0.9999979939221901 and 0.9999836379973899. With shell radius,
map scale and power equal to 1, 1 and 2, the radial law
$r = 1 + ((u - 1/2)/(1-u))^2$ gives approximately
$6.21 times 10^10$ and $9.34 times 10^8$. The coordinate vectors, runtime
overlay, executable inspection card and detailed precision metadata are
preserved in `sunrise/failure-inspect/`.

=== Matched Arb-only control

To measure sampling without those erroneous contributions, both modes were
rerun with Arb as the sole precision level. The budget and all five seeds are
unchanged. Each seed batch ran on a separate pinned CPU (110–114); modes within
each batch ran sequentially. Reported timings are per run, not the concurrent
batch's elapsed time. CPU time is within 1% of wall time.

#table(
  columns: 5,
  table.header([*Mode*], [*Five-seed mean ± error*], [*RMS error/run*],
    [*Wall/run*], [*Time efficiency*]),
  [Ordinary], [−0.0161290 ± 0.0002164], [3.00%], [115.9 s], [1.00×],
  [\+ Fermi pair], [−0.0158916 ± 0.0002034], [2.82%], [161.0 s], [0.815×],
)

All ten Arb-only runs agree with the reference within 0.96 reported standard
errors, and none rejects samples. The Fermi mean is 1.15 standard errors above
the reference. The observed matched-seed efficiency range is 0.60–1.02×.
The variance gain is only 1.13×, while runtime rises by 39%. Mapping time rises
from 30 to 439 microseconds per draw; physical evaluation changes much less,
from 1054 to 1084 microseconds. Approximately 31.9% of augmented-mode draws
select the Fermi channel. *Fermi sampling does not pay off for sunrise at this
budget*, even in the controlled precision comparison. These timings should
not be conflated with those of the much faster default ladder.

#figure(image("sunrise/arb_only/comparison.png", width: 100%),
  caption: [Exact-zero-temperature sunrise, with Arb for every sample.])

=== Practical Double-to-Arb control

The same five-seed, 100,000-draw comparison was also run with the Quad level
removed, preserving the Double and Arb tolerances and every other setting.
This tests a less costly way of avoiding the observed false Quad acceptance.
The full Arb-only results above remain a separate numerical cross-check.

#table(
  columns: 5,
  table.header([*Mode*], [*Five-seed mean ± error*], [*RMS error/run*],
    [*Wall/run*], [*Time efficiency*]),
  [Ordinary], [−0.0159352 ± 0.0001833], [2.54%], [18.92 s], [1.00×],
  [\+ Fermi pair], [−0.0158273 ± 0.0001957], [2.71%], [55.11 s], [0.301×],
)

All ten runs are within 1.97 reported standard errors of the reference;
the two means have pulls of 1.04 and 1.53. No sample is rejected. About 11% of
draws use Arb in either mode. The matched-seed efficiency range is 0.16–0.44×.
Mapping takes 28 versus 376 microseconds per draw, while physical evaluation
takes 133 versus 136 microseconds. This reinforces the sunrise conclusion:
the Fermi map's cost is substantial compared with this simple integrand.
The 20 final positive/negative extrema across all ten runs agree with forced
Arb to better than $6 times 10^(-11)$ relative difference.
The small variance differences between precision controls also reflect
different adaptive-grid histories and finite-sample fluctuations. No timing
ratio mixes different precision policies.

== Three-loop precision checks

The uncrossed tennis-ball graph also exhibits false Quad acceptance. In the
three-shell channel, seed 2843 visits two points whose default inspections give
+15.6043241 and −38.2691467. Forced Arb gives approximately
$1.0290614 times 10^(-94)$ and $-1.0084201 times 10^(-294)$ at the same coordinates.
The Double-to-Arb ladder returns those same Arb values without forcing precision.
Eight other early extrema, from the other four seeds, agree with Arb.
The complete ten-point diagnostic and a runnable card are saved in
`tennis_ball/failure-inspect/`.

The 200-draw tennis-ball cost pilot gives matching default and Arb-only
estimates in all three modes. Arb-only evaluation costs approximately 69–71 ms
per draw. The main precision control therefore uses Double-to-Arb for all
sampling modes, with the same 50,000 draws, ten iterations and five seeds as
the default-ladder scan. Default and control results are stored separately;
failed seeds are retained. Agreement with Arb at selected extrema is a
diagnostic, not a proof that every sampled value is numerically accurate.

Mercedes exhibits the same problem. Two extreme points from the three-shell
seed 2843 inspect as +26.2962 and −1,342,425.08 under the default ladder,
accepted as stable in Quad. Forced Arb gives approximately
$-1.5201656 times 10^(-300)$ and $-3.3400669 times 10^(-264)$.
The Double-to-Arb ladder reproduces those Arb values exactly. Its diagnostic
cards and detailed results are in `mercedes/failure-inspect/`.
Thus the observed precision failure affects all three tested diagrams.

== Uncrossed tennis-ball results

The primary usable comparison uses Double-to-Arb, 50,000 draws per run,
ten equal-size iterations, and the same five seeds. All modes use the same
precision policy. The table reports per-run RMS errors and mean wall times:

#table(
  columns: 5,
  table.header([*Mode*], [*Five-seed mean ± error*], [*Error/run*],
    [*Wall/run*], [*Efficiency*]),
  [Ordinary], [+0.0114506 ± 0.0003307], [6.65%], [350.3 s], [1.00×],
  [\+ Three shells], [+0.0114259 ± 0.0002859], [5.75%], [439.1 s], [1.074×],
  [\+ Repeated shell], [+0.0115024 ± 0.0003007], [6.05%], [410.1 s], [1.037×],
)

The three mean pulls are 1.00, 1.06 and 1.27. Fourteen of fifteen runs are
within two reported standard errors; the remaining ordinary run is 2.25
standard errors high. No sample is rejected. All 30 final positive/negative
extrema agree with forced Arb to better than $4 times 10^(-10)$ relative
difference. These checks support agreement at the measured Monte Carlo
precision, rather than certifying the implementation at arbitrary accuracy.

*Both Fermi choices are near break-even.* The observed matched-seed efficiency
ranges are 0.85–1.50× for three shells and 0.86–1.47× for the repeated shell.
The mean gains of 7% and 4% are small compared with this variation. The variance
reductions are mostly consumed by the extra runtime. Mapping takes 0.12, 1.45
and 1.40 ms per draw for ordinary, three-shell and repeated-shell sampling;
physical evaluation takes 6.74, 7.15 and 6.65 ms. Roughly 9–10% of draws use Arb.
The Fermi channel receives 12.2% and 9.0% of draws in the two augmented modes.

With the default ladder, the repeated-shell variant measures only 0.79×
efficiency: 6.65% → 6.05% error but 124.9 → 190.5 s wall time. Its ten extrema
agree with Arb to better than $4 times 10^(-12)$. The three-shell variant's
seed 2843 instead finishes at −0.287629 ± 0.041368, a 7.22-standard-error
failure, so its variance-based efficiency is invalid. This is why timing
conclusions must identify the precision policy.

#figure(image("tennis_ball/double_arb/comparison.png", width: 100%),
  caption: [Exact-zero-temperature uncrossed tennis ball, Double-to-Arb.])

== Mercedes results

Both sampling modes use 50,000 draws, ten equal-size iterations and the same
five seeds. With Double-to-Arb:

#table(
  columns: 5,
  table.header([*Mode*], [*Five-seed mean ± error*], [*Error/run*],
    [*Wall/run*], [*Efficiency*]),
  [Ordinary], [−0.0310531 ± 0.0005281], [3.83%], [600.4 s], [1.00×],
  [\+ Three shells], [−0.0308384 ± 0.0004791], [3.48%], [733.3 s], [0.996×],
)

The two mean pulls are −0.48 and −0.08. All ten runs agree with the reference
within 1.49 reported standard errors, with no rejected samples. The 20 final
extrema agree with forced Arb to better than $8 times 10^(-11)$ relative
difference. The matched-seed efficiency range is 0.80–1.17×.

*The variance reduction is consumed by the added runtime.* Mapping rises
from 0.155 to 1.880 ms per draw, and physical evaluation from 11.636 to
12.483 ms. Approximately 15% of draws use Arb in either mode, and the Fermi
channel receives 10.2% of augmented-mode draws.

The default-ladder ordinary mean is −0.0310523 ± 0.0005280 at 256.8 s/run.
The augmented mean is −0.00225394 ± 0.01884146, including seed 2843's
+0.111124 ± 0.094183. Despite the confirmed numerical failure, every individual
default run lies within 1.60 reported standard errors of the target, and the
augmented mean is within 1.52. The inflated errors mask the defect. Agreement
with the reference and a zero rejection counter do not, by themselves,
validate these runs. The invalid augmented efficiency is not used in the
conclusion.

#figure(image("mercedes/double_arb/comparison.png", width: 100%),
  caption: [Exact-zero-temperature Mercedes, Double-to-Arb.])

== Setup cost and interpretation

Fresh integrand generation and compilation took 3.58 s for sunrise,
1490.04 s (24 min 50 s) for tennis ball, and 2657.30 s (44 min 17 s) for
Mercedes. These are shared setup costs, paid once before reusing the same
state for every mode and seed. They are excluded from the per-integration
efficiencies above. Building the CLI is also excluded. The CLI was built
with `nix build .#gammaloop` at the stated revision; generated assembly
evaluators use `O3` with fast math disabled. Binary hashes, generation-card
hashes, CPU affinities and individual CPU/wall times are retained with the
measurements. Generation logs are preserved as `generate.log.gz`.

The experiment contains 80 full integration runs, totaling 5.5 million draws,
plus separate 200-draw cost pilots and point diagnostics. It compares the
specific channels and fixed localization profile described above. The exact
$T = 0$ implementation analytically localizes distribution terms, so the
sampling problem differs from resolving a narrow finite-temperature peak.
The earlier finite-beta improvements need not persist in this representation.
The code rebase also prevents attributing differences from the earlier report
solely to temperature. These measurements do not rule out benefits for other
integrands, channel choices or localization profiles.

== Reproduction

Use a binary built from the measured revision, inside `nix develop` with its
configured Symbolica runtime license. For each of `sunrise`, `tennis_ball` and
`mercedes`, run the existing driver with `--zero-temperature --generate` into a
fresh output directory. Subsequent runs reuse that state, omitting `--generate`.
The driver writes the exact generated card and its hash into the output and
checks the saved binary and source identities before reuse. The generation
card changes only the medium mode and fixed localization profile relative to
the corresponding finite-temperature card.

```sh
python examples/cli/eos/fermi_surface_sampling/benchmark.py \
  --example sunrise --zero-temperature --generate \
  --binary /path/to/gammaloop --output /tmp/fermi-T0-sunrise --cpu 96 \
  --samples 100000 --iterations 10 --seeds 1337 2843 4517 6007 7919
python examples/cli/eos/fermi_surface_sampling/summarize.py \
  /tmp/fermi-T0-sunrise/results.json --output /tmp/fermi-T0-sunrise/summary
```

Exact $T = 0$ defaults to ordinary and augmented catalogues, adding the
single-shell control for tennis ball. `--betas` and `--zero-temperature` are
mutually exclusive. The existing finite-temperature benchmark remains available.
Use `--samples 50000` for the two three-loop comparisons and retain the same
five seeds and ten iterations.

For the precision controls, add `--arb-only` or `--skip-quad` and use a separate
output directory with the same generated state and its `generation.json`.
These mutually exclusive switches copy the retained precision-level settings
from the generation card. The short sunrise failure reproducer can be run from
the repository root without doing an integration scan:

```sh
benchmark_dir=examples/cli/eos/fermi_surface_sampling
sunrise_repro_dir="$benchmark_dir/zero_temperature/sunrise"
gammaloop -s /tmp/fermi-T0-sunrise-reproducer \
  -n "$sunrise_repro_dir/generate.toml" quit -n
gammaloop -s /tmp/fermi-T0-sunrise-reproducer --read-only-state \
  -n "$sunrise_repro_dir/failure-inspect/reproduce.toml" quit -n
```
