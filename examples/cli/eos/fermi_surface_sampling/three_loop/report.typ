= Fermi-surface sampling on two three-loop thermal graphs

== Tennis-ball result

The three-shell Fermi channel improves estimated time efficiency at $beta = 100$
by a factor of 1.66. It costs more at $beta = 1$ and 10, and is approximately
break-even at $beta = 1000$. A benefit therefore does not increase monotonically
with beta in this measured setup.

#table(
  columns: 6,
  table.header([*$beta$*], [*Ordinary error*], [*+ Fermi error*],
    [*Ordinary time*], [*+ Fermi time*], [*Efficiency*]),
  [1], [6.755%], [6.694%], [46.49 s], [112.91 s], [0.42],
  [10], [5.482%], [4.464%], [37.53 s], [106.53 s], [0.54],
  [100], [15.036%], [9.051%], [109.72 s], [182.45 s], [1.66],
  [1000], [34.840%], [34.449%], [474.19 s], [503.92 s], [0.95],
)

Errors are RMS reported errors per run, relative to the reference. Times are
means across three seeds. Efficiency is ordinary mean $sigma^2 t$ divided by
augmented mean $sigma^2 t$. At $beta = 100$, the individual seed ratios range
from 1.21 to 2.52; at $beta = 1000$, they span 0.39 to 1.92. These observed ranges
are not confidence intervals. The coldest estimates have large uncertainty, so
the approximately unit efficiency does not rule out a modest benefit or penalty.

All 24 primary runs agree with the supplied target within 1.74 reported standard
errors; all eight combined three-seed means agree within 1.53. No sample is
rejected as NaN or unstable. The Fermi channel receives 6.9%, 11.3%, 17.0% and
30.8% of draws as beta increases, confirming that it remains active at the
coldest temperature.

Mapping costs about 114–128 microseconds per draw for ordinary channels versus
1398–1440 microseconds with the Fermi channel. At $beta = 1000$, physical
evaluation costs 9.19 versus 8.46 milliseconds per draw, and the Arb fractions
are 16.6% versus 14.7%. The expensive physical evaluation largely amortizes the
mapping overhead, but the three-shell map provides little variance reduction
there at this sample budget.

#image("tennis_ball/comparison.png", width: 100%)

Primary measurements and all counters are in
#link("tennis_ball/results.json")[tennis_ball/results.json], with tabulated
metrics in #link("tennis_ball/summary.csv")[tennis_ball/summary.csv].

=== Repeated-momentum diagnostic

After observing the primary result, a separate exploratory comparison focused
only the repeated fermion momentum (edges 0 and 5), using
`product(block(lmb(0),fermi(0)),complement(1,3))` with parent `[0,1,3]`.
It retains the same ordinary optimized catalogue and uses the same three seeds,
50,000 draws and ten iterations at each of $beta = 100$ and 1000. The original
primary results remain unchanged.

#table(
  columns: 5,
  table.header([*$beta$*], [*Single-shell error*], [*Single-shell time*],
    [*Efficiency*], [*Observed seed range*]),
  [100], [9.454%], [184.84 s], [1.50], [1.11–1.97],
  [1000], [17.864%], [495.64 s], [3.63], [3.07–4.37],
)

At $beta = 1000$, this reduces variance by a factor of 3.80 while increasing
wall time by only about 4.5%. Mapping still costs 1.38 milliseconds per draw,
close to the three-shell map's 1.44 milliseconds. The gain comes primarily from
variance reduction, not cheaper mapping. All three seed comparisons favor the
single-shell channel. This establishes a useful measured low-temperature case
for the feature and shows why channel choice matters.

All six diagnostic runs agree with the target within 2.54 reported errors;
the combined means at $beta = 100$ and 1000 are 0.32 and −1.76 errors from the
reference. These remain Monte Carlo checks at finite precision, not a general
proof that the map is unbiased. The roughly 18% individual-run error at the
coldest point is still substantial.

#image("tennis_ball/repeated_momentum/comparison.png", width: 100%)

The exploratory results and matched ordinary baselines are preserved in
#link("tennis_ball/repeated_momentum/results.json")[repeated_momentum/results.json].

== Mercedes result

The three-shell channel does not repay its overhead at $beta = 1$, 10 or 100.
At $beta = 1000$, it gives a nominal efficiency gain of 1.52, with all three
seed ratios between 1.11 and 2.23. The reference check below qualifies that cold
result: variance-based efficiency alone is insufficient to certify accuracy.

#table(
  columns: 6,
  table.header([*$beta$*], [*Ordinary error*], [*+ Fermi error*],
    [*Ordinary time*], [*+ Fermi time*], [*Efficiency*]),
  [1], [7.792%], [7.596%], [86.23 s], [181.97 s], [0.50],
  [10], [2.434%], [2.322%], [77.89 s], [175.18 s], [0.49],
  [100], [3.377%], [3.213%], [472.28 s], [592.59 s], [0.88],
  [1000], [8.423%], [6.519%], [563.29 s], [613.21 s], [1.52],
)

The last row uses 20,000 draws per seed; the other rows use 50,000. At the coldest
point, physical evaluation costs 27.42 versus 28.15 milliseconds per draw and
mapping costs 0.178 versus 1.911 milliseconds. Arb handles 28.5% versus 26.9%
of draws. The Fermi channel receives 10.2% of samples. There are no rejected
NaN or unstable samples in any primary run; real components stay at rounding
scale (below $3 times 10^(-17)$ in magnitude).

#image("mercedes/comparison.png", width: 100%)

=== Reference-validation limitation

Every individual primary run is within 2.33 reported errors of its supplied
target. However, the three-seed Fermi mean at $beta = 1000$ is
$-0.03389618 plus.minus 0.00115927$, or 2.67 errors below the reference
$-0.030799358173484585$. An independent seed 9991, with 20,000 draws in each
mode, gives:

#table(
  columns: 3,
  table.header([*Mode*], [*Independent estimate*], [*Reference pull*]),
  [Ordinary], [−0.03112346 ± 0.00186873], [−0.17],
  [Three Fermi shells], [−0.03338256 ± 0.00176779], [−1.46],
)

Pooling all four equally weighted seeds gives ordinary
$-0.03240350 plus.minus 0.00121656$ (−1.32 errors from reference) and augmented
$-0.03376778 plus.minus 0.00097533$ (−3.04 errors). The two sampling modes remain
compatible at their present precision, so this does not isolate bias in the
Fermi map. It does leave a convergence concern against the supplied target.
The cold Mercedes case is therefore *not fully validated* by this measurement;
its 1.52 efficiency should be read as a nominal reported-variance result,
not an established gain at verified accuracy. The independent validation is excluded from the primary
performance table and retained in
#link("mercedes/validation-results.json")[validation-results.json].

Primary data and metrics are in
#link("mercedes/results.json")[mercedes/results.json] and
#link("mercedes/summary.csv")[mercedes/summary.csv].

== Setup

This extends the #link("../report.pdf")[sunrise investigation] to the user's
thermal tennis-ball (GL3) and Mercedes (GL1) examples. The original cards remain
unchanged under `ignore/my_runs/thermal_tennis_ball/` and
`ignore/my_runs/thermal_mercedes/`. The local `generate.toml` copies use the
current `sampling_multichanneling` and `sampling_channels` setting names, select
massless bottom quarks as in the supplied integration blocks, and write to
isolated states. The exported DOT files record the exact generated graphs.

Both examples have `aS = 1`, baryon chemical potential `muB = 3`, quark chemical
potential $mu = 1$, and Fermi radius $p_F = 1$. Vacuum subtraction, local and
integrated UV counterterms, summed physical orientations, cached evaluation,
assembly compilation, and the three-level Double/Quad/Arb stability policy are
retained from the user cards. `m_uv = 1` and
`mu_r = renormalization_localization_scale = 2`. Threshold subtraction is
disabled. No production code is changed.

The ordinary catalogue uses `auto:optimized_lmb`. The augmented catalogue adds
one `fermi_shells` channel with parent basis `[0,1,3]` and composition
`product(block(lmb(0),fermi(0)),block(lmb(1),fermi(1)),block(lmb(3),fermi(3)))`.
These three fermion momenta are independent in both graphs. The tennis-ball
edge 5 shares edge 0's momentum and is not an additional independent block.
All ordinary optimized LMBs remain present to cover soft gluon momenta.
The Fermi proposal has the implementation's fixed radial law, rather than a
temperature-dependent width. Both compared catalogues explicitly use
`power = 2` and `sampling_channel_weight = "map_density"`, as in the sunrise
benchmark; the original cards leave power at its default of 1. With power 2,
the Fermi density grows as the inverse square root of distance to the shell.
The ordinary spherical mapping remains linear, with `b = e_cm = 1`.

The supplied reference values are:

#table(
  columns: 3,
  table.header([*$beta$*], [*Tennis ball, Im I*], [*Mercedes, Im I*]),
  [1], [0.7113504235805803], [−1.739558840398423],
  [10], [0.01398321747318332], [−0.04062795126001039],
  [100], [0.011149085893291288], [−0.030915837358603345],
  [1000], [0.011121756598238203], [−0.030799358173484585],
)

These are user-supplied references; their analytic provenance is not established
by this investigation. The two colder targets were supplied during measurement
as temperatures $T = 0.01$ and 0.001. Tennis ball is the user's “uncrossed” graph;
the supplied crossed-graph targets are not used. Some tennis-ball runs therefore
did not have a CLI display target at execution time; their pulls are calculated
from the supplied values afterwards, without changing the sampled integral.

== Measurement convention

Measured on 2 October 2026 (UTC). The production revision is
`2a236eaab703fe7e9b5505522bf6884efa4b9ef7`, using the same
Nix `ci-optim` binary as the sunrise investigation. Each temperature and sampling
mode uses seeds 1337, 2843 and 4517, with ten adaptation iterations. Tennis ball
uses 50,000 draws per run at every temperature. Mercedes uses 50,000 at
$beta = 1, 10, 100$ and 20,000 at $beta = 1000$. Both sampling modes always have
the same budget at a given temperature. The cold Mercedes budget was selected
for its measured evaluation cost, before using any completed estimate to judge
convergence.

Each integration starts a fresh workspace with one worker on a pinned CPU.
Sampling modes share all physics and stability settings; their adaptive grids
train independently. Jobs within each batch are shuffled deterministically.
The driver records wall and child-process CPU time, reported integration error,
precision escalation, rejected samples, and the actual channel population.
Wall time includes state loading, warm-up and result persistence. Graph
generation is timed separately and shared across sampling modes.

Generation took 993.01 seconds for tennis ball and 1291.111 seconds for Mercedes.
The two graphs were generated concurrently; generation used a broader thread
pool than the integration measurements. These setup times describe this run,
not a separately optimized generation benchmark. The generated states are
reused by every temperature, seed and sampling mode.

The integration batches use distinct CPUs on the shared AMD EPYC 9754 host.
Two tennis-ball baseline timings at $beta = 10$ overlapped an early Mercedes
pilot on the same CPU. Both were rerun on dedicated cores and reproduced the
original estimates and errors exactly; only the replacement timings enter the
comparison. Superseded records are retained in the data. Incomplete exploratory
runs and the 200-draw timing pilots are excluded from all primary estimates.

Efficiency is the ratio of ordinary to augmented mean $sigma^2 t$ across seeds.
Values greater than one favor Fermi sampling. This is an estimate of time to
equal error at the measured budget, not a direct run to fixed precision.
Combined estimates use equal seed weights; their propagated error is
$sqrt(sum_i sigma_i^2)/n$. Relative errors use the supplied target when available,
otherwise the measured mean. Raw absolute errors remain available so that a
small or poorly converged mean cannot silently imply a variance gain.

== Reproduction

From the repository root, in `nix develop`, using the measured binary and its
configured runtime Symbolica license:

```sh
python examples/cli/eos/fermi_surface_sampling/benchmark.py \
  --example tennis_ball --binary /path/to/gammaloop \
  --output /tmp/fermi-tennis-ball --generate --cpu 96 --samples 50000
python examples/cli/eos/fermi_surface_sampling/benchmark.py \
  --example mercedes --binary /path/to/gammaloop \
  --output /tmp/fermi-mercedes --generate --cpu 97 --samples 50000 --betas 1 10 100
python examples/cli/eos/fermi_surface_sampling/benchmark.py \
  --example mercedes --binary /path/to/gammaloop \
  --output /tmp/fermi-mercedes --cpu 97 --samples 20000 --betas 1000
python examples/cli/eos/fermi_surface_sampling/summarize.py \
  /tmp/fermi-tennis-ball/results.json --output /tmp/fermi-tennis-ball/summary
```

The driver records binary and generation-card digests and keeps each run's
command card, settings overlay, log and integration workspace.
Separate invocations overwrite the aggregate `results.json` but retain all
per-run `record.json` files. Combine both Mercedes budgets with:

```sh
python - <<'PY'
import json
from pathlib import Path
root = Path("/tmp/fermi-mercedes")
data = json.loads((root / "results.json").read_text())
data["runs"] = [json.loads(p.read_text()) for p in root.glob("beta_*/record.json")]
data["metadata"]["samples"] = None
(root / "all-results.json").write_text(json.dumps(data, indent=2))
PY
python examples/cli/eos/fermi_surface_sampling/summarize.py \
  /tmp/fermi-mercedes/all-results.json --output /tmp/fermi-mercedes/summary
```

For the exploratory tennis-ball channel, use a generated tennis-ball output
directory with `--modes optimized_fermi_repeated --betas 100 1000 --samples 50000`.
The independent Mercedes check uses `--seeds 9991 --betas 1000 --samples 20000`.
The saved measured aggregate records each original batch and CPU; per-run logs
and resumable workspaces remain under `/tmp/fermi-three-loop/`.
