= Fermi-surface sampling on the thermal sunrise

== Measured result

Adding the independent Fermi-pair channel to ordinary optimized LMB sampling
does not improve this sunrise benchmark, including at $beta = 1000$. Its
integration uncertainty is nearly unchanged while the map adds substantial cost.

Measured on 3 October 2026 (Europe/Zurich), at revision
`2a236eaab703fe7e9b5505522bf6884efa4b9ef7`, using the Nix `ci-optim` binary
(GammaLoop optimization level 2, dependency optimization level 3, no debug
information) and assembly-compiled integrand evaluators. Each entry uses three
seeds, each with 100,000 samples in ten iterations. Error percentages below are
RMS reported errors of individual runs, relative to the analytic target; times
are arithmetic means. Generation took 3.22 seconds, shared by every run.

#table(
  columns: 6,
  table.header([*$beta$*], [*Ordinary error*], [*+ Fermi error*],
    [*Ordinary time*], [*+ Fermi time*], [*Efficiency*]),
  [1], [0.655%], [0.693%], [6.05 s], [40.18 s], [0.135],
  [10], [2.092%], [2.066%], [8.27 s], [43.46 s], [0.195],
  [100], [2.458%], [2.537%], [26.65 s], [64.37 s], [0.389],
  [1000], [2.529%], [2.666%], [44.00 s], [79.54 s], [0.498],
)

Efficiency is the ordinary catalogue's mean $sigma^2 t$ divided by the Fermi
catalogue's mean $sigma^2 t$; larger is better. The corresponding estimated
equal-error time penalties are 7.42, 5.13, 2.57 and 2.01. They describe convergence
at this sample budget, not an independently timed run to a fixed target precision.

#image("comparison.png", width: 100%)

All 24 primary-comparison runs agree with the analytic value within 1.61 reported
standard errors. All eight three-seed means agree within 1.57 standard errors.
The real component and its uncertainty are exactly zero. Across all 60 runs
(six million draws, including diagnostic modes), no NaN or unstable sample was
rejected. Fermi channels received 16.2%, 31.0%, 32.4% and 32.3% of draws as beta
increased; the map was actively exercised.

The average mapping cost rises from 26–33 microseconds per draw to 354–390
microseconds. The colder runs have more arbitrary-precision evaluation under the
unchanged stability policy: at $beta = 1000$, the ordinary and augmented
catalogues use Arb for 21.5% and 19.8% of samples, respectively. This makes the
relative map overhead smaller at high beta; it does not produce a variance gain.
CPU times are within 1% of wall times, and matching preliminary sequential runs
reproduce identical numerical values with timings within about 2% of these
concurrent batches.

A more expensive integrand can amortize the roughly 330–360 microsecond extra
mapping cost, but a useful speedup also requires variance reduction. This
sunrise shows no substantial reduction. Its only nominal improvement is a 2.5%
variance reduction at $beta = 10$; holding that noisy estimate fixed would require
about 13.8 milliseconds of extra common work per draw to break even. This is a
cost-model extrapolation, not evidence for a more complicated integrand.

The separate #link("three_loop/report.pdf")[three-loop follow-up] measures the
user's tennis-ball and Mercedes examples, including a targeted single-shell
channel that improves cold tennis-ball efficiency.

=== Convergence diagnostic

The single-LMB, Fermi-only and same-parent mixed modes omit explicit ordinary
sampling of the soft gluon momentum. At $beta = 10$, individual single-LMB and
Fermi-only runs miss the reference by 5.31 and 4.90 reported errors, respectively;
the three-seed single-LMB and mixed estimates miss by 4.47 and 3.82 errors.
Their reported errors also vary greatly across seeds. These failures demonstrate
that a small nominal error from a single such run is insufficient validation.
The optimized catalogue includes gluon-momentum LMBs and passes the independent
reference checks, including a separate 200,000-sample control. This behavior is
consistent with undersampling the soft-gluon region in the single-parent maps.
Keep ordinary soft-momentum coverage when adding Fermi channels.

The complete measurements are in #link("results.json")[results.json]; all-mode
metrics are in #link("summary.csv")[summary.csv] and
#link("summary.json")[summary.json]. The independent control is in
#link("validation-results.json")[validation-results.json]. Per-run command cards,
logs and resumable workspaces remain under
`/tmp/fermi-sunrise-benchmark/scan/`. No production implementation was changed.

== Scope and reference

This benchmark compares changes of variables for the physical massless quark–gluon
sunrise, using `tests/resources/graphs/sunrise_qcd_vacuum.dot`. It tests correctness,
variance reduction and sampling overhead on a simple two-loop integral. It does
not measure the performance of more complicated integrands.

The model is SU(3), with one massless quark, `aS = 1`, `muB = 3`, hence quark
chemical potential $mu = 1$ and Fermi radius $p_F = 1$. Generation includes vacuum
subtraction and local plus integrated UV counterterms, with
`m_uv = mu_r = renormalization_localization_scale = 2`. Threshold subtraction is
disabled. All physical orientations are summed.

For inverse temperature $beta$, the imaginary integral is compared with
$
  I_"ref" = - (5 pi / (18 beta^4) + 1 / (pi beta^2) + 1 / (2 pi^3)).
$
This is the one-flavor part of the order-$g^2$ pressure in
#link("https://arxiv.org/pdf/hep-ph/0305183")[Vuorinen, Eq. (3.15)], with
$g^2 = 4 pi$, omitting the pure-gluon contribution. The real component must vanish.

== Controlled comparison

`ordinary` uses the two fermion momenta as its ordinary spherical LMB.
`fermi` replaces both radial maps with the product
`product(block(lmb(0),fermi(0)),block(lmb(1),fermi(1)))`.
`mixed` samples one of these two channels per draw. All use `map_density`
partitioning, `power = 2`, `b = e_cm = 1`, and the same physics and stability settings.
The primary comparison is `optimized` versus `optimized_fermi`: automatic
optimized LMB coverage, respectively without and with the additional Fermi
channel. The ordinary catalogue contains `(1,2)`, `(0,2)` and `(0,1)`; the first
two explicitly sample the massless gluon momentum. The single-parent modes serve
as convergence diagnostics.

Both pipelines share model loading, generation, evaluator compilation, stability
checks, thermal physics, subtraction, and accumulation. The first difference is
the sampling catalogue and its map/partition; subsequent adaptive grids train
independently. Every run starts a fresh integration workspace. Jobs are shuffled
deterministically; one integration worker is used throughout. The recorded scan
runs one sequential batch per temperature on CPUs 96, 97, 98 and 99, respectively,
with the four temperature batches concurrent. The host is a shared AMD EPYC 9754
machine. Matching preliminary sequential runs provide a check on timing repeatability.

For each mode and temperature, the analysis reports RMS integration uncertainty
over seeds and the arithmetic mean runtime. The combined estimate is the
unweighted mean over independent seeds, with error
$sqrt(sum_i sigma_i^2) / n$. Efficiency is the optimized ordinary catalogue's mean
$sigma^2 t$ divided by the candidate's mean $sigma^2 t$; values above one favor the
candidate. Both the single-LMB and optimized-LMB comparison columns are retained
in the machine-readable summary. This is an estimate of equal-error time efficiency at the measured
sample budget, not a measured long-run speedup. Startup, loading, warm-up and
result persistence are included in run wall time. Generation is timed separately.
CPU time, internal mapping/evaluation time, precision escalation and rejected
sample fractions are retained in the raw results.

When variance improves by a factor $G > 1$, an extra common evaluation cost $c$
per sample would give the approximate efficiency
$G (t_o + c) / (t_f + c)$, where $t_o$ and $t_f$ are measured per-sample costs.
The break-even extra cost is $max(0, (t_f - G t_o)/(G - 1))$.
This holds the measured variance and adaptation behavior fixed; it does not
establish a gain for a different physical integrand.

== Reproduction

Run from the repository root inside `nix develop`, with a binary built from the
measured revision and its configured Symbolica runtime license:

```sh
python examples/cli/eos/fermi_surface_sampling/benchmark.py \
  --binary /path/to/gammaloop --output /tmp/fermi-sunrise \
  --generate --cpu 96
python examples/cli/eos/fermi_surface_sampling/summarize.py \
  /tmp/fermi-sunrise/results.json --output /tmp/fermi-sunrise/summary
```

The driver uses only the Python 3.11+ standard library; plotting additionally requires
Matplotlib. Each run retains its overlay, command card, full log, integration
workspace and exact result JSON. `results.json` records source revision, binary
identity and generation-card digest. The default scan is $beta = 1, 10, 100, 1000$,
three seeds (1337, 2843, 4517), five sampling modes, and 100,000 draws in ten equal
iterations.
