# Small positive virtuality without resummation

The strict fixed-order reference remains finite at each tested positive
`K^2`, but grows towards the light cone. GammaLoop completes these runs with
large Monte Carlo errors and occasional nonfinite samples. Compatibility
within these large errors does not establish a precision validation.

The scan holds `T=1`, `k0=4`, `e=g=1`, `Nc=3`, and the UV scales fixed, and sets
`kz=sqrt(16-K^2)`. It uses the same massless down-quark model, current projector,
and normalization as [README.md](README.md): `rho_V = 36*Gamma_rho`. Values
below are `rho_V/T^2`; NLO denotes the correction alone. Local and integrated
UV subtraction remain enabled, with `medium.vacuum_subtraction=false`.
No resummation, soft counterterms, or threshold subtraction is used.

## Physical fixed-order behaviour

The author's [reference code](https://doi.org/10.5281/zenodo.3478144) for
[arXiv:1910.07552v2](https://arxiv.org/pdf/1910.07552v2) gives:

| K²/T² | kz/T | LO | NLO correction | NLO/LO at g=1 |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 3.872983346207417 | 0.1595160 | 0.1080440 | 0.677 |
| 0.1 | 3.987480407475377 | 0.01582953 | 0.2405138 | 15.19 |
| 0.01 | 3.998749804626441 | 0.001581740 | 0.4057995 | 256.55 |

LO vanishes proportionally to `K^2` at fixed positive `k0`. The NLO correction
instead has logarithmic light-cone growth, of order
`g^2*T^2*log(T^2/K^2)`. The paper discusses this in Eq. (5.10), and identifies
`|K^2|` of order `g^2*T^2` or smaller as a regime requiring resummation. Thus
even an exact integration of this unresummed function would not provide a
controlled physical rate in the two smaller-mass cases at `g=1`.

The reference was evaluated with absolute quadrature tolerances `1e-4` and
`1e-5`. The changes in NLO are respectively `2.72e-8`, `1.85e-7`, and
`9.28e-8`, far below the native errors. As in the main benchmark, this is a
convergence check of the author's code, not a rigorous combined error bound.

## Native numerical behaviour

The scan keeps `e_cm=sqrt(15)` as a fixed parameterization scale so that only
the external spatial momentum changes between points. Graphs, orientations,
and OSE-weighted LMB channels are summed, with eager `SummedFunctionMap`
evaluation and the existing Double/Quad/Arb stability ladder.

| K²/T² | Samples | Seed | Native NLO | Nonfinite samples discarded |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 100,000 | 191007552 | 0.0379 ± 0.0479 | 0 |
| 0.1 | 100,000 | 191007552 | 0.0687 ± 0.1271 | 0 |
| 0.01 | 100,000 | 191007552 | 0.7299 ± 0.3960 | 0 |
| 0.1 | 500,000 | 271828 | 0.2040 ± 0.0725 | 4 |
| 0.01 | 500,000 | 271828 | 0.4738 ± 0.0890 | 8 |

All five integrations completed. The last two differ from the reference by
-0.50 and +0.76 quoted Monte Carlo errors, but their relative errors remain
36% and 19%. The largest individual weighted sample is 31% and 11% of the
absolute mean times sample count, respectively. In the first 100,000-sample
runs those diagnostics are 78%, 143%, and 52%. Large sample weights limit
confidence in error estimates from these run sizes.

Ordinary precision escalation remains relatively rare: the two longer runs
finish 99.89% and 99.88% of samples at Double. Nevertheless, 4 and 8 points
remain nonfinite at the final precision and are replaced by zero. Their small
count alone does not bound the bias from discarding them. These runs did not
encounter the fatal radial-root isolation error seen in the spacelike test;
that does not rule it out at other samples or smaller virtualities.

LO was also integrated with 100,000 samples per point, giving respectively
`0.152609 ± 0.007584`, `0.014304 ± 0.001862`, and
`0.001715 ± 0.000507`. All agree with the analytic LO reference, while their
sampling errors also worsen as the final-state momenta become collinear.

A control at `K^2=0.01` changed only the parameterization scale to
`e_cm=sqrt(K^2)=0.1`. Its 100,000-sample NLO estimate is
`-0.0042 ± 0.6593`, with two discarded nonfinite samples. This particular scale
change did not improve convergence. No optimization of the sampling was
attempted, and no conclusion about a divergent mathematical variance follows
from these finite runs.

## Reproducible nonfinite point

At `K^2=0.01`, the following unit-hypercube point reproduces a nonfinite final
evaluation under both the default stability ladder and forced Arb precision:

```text
0.38284976323786246 0.42928160419108263 0.010235373341083381
0.9998525100417789  0.43510709710570933 0.010667453207883319
```

Its parameterization Jacobian is `1.1376501886039257e21`. Both graphs also
produce the warning when selected separately with the existing one-graph
Monte Carlo selector. It has therefore not been isolated to the raised
self-energy cut. The generated expression and cut tables at all three scan
points are exactly identical to the hard timelike benchmark's standalone
export. The first changed inputs are the external momentum and the resulting
radial roots/parameters during numerical evaluation. A different generated
graph set or residue table cannot explain this scan's deterioration.

The diagnostic warning is `process_final_nonfinite_sample`, with
`result="(NaN+NaNi)"`. **The subsequent inspect output is zero because the
sample was discarded, not because the integrand is zero.** The raw inspect
JSON has no evaluation metadata exposing that distinction. The warning
records and sample coordinates are retained in `lightcone_results.json`.
The cause inside the individual graph evaluations remains unresolved; no
core fix was attempted.

## Reproduction

From the repository root, run the pilot for each listed `kz` (using a distinct
output prefix and state directory for each point):

```sh
./gammaloop --dev-optim -n -L warn \
  -s /tmp/gammaloop-qcd-lightcone/01-state \
  examples/cli/thermal_production/qcd_vector/run.toml \
  run model near_lightcone lo nlo export integrate \
  -D kz=3.987480407475377 -D samples=100000 -D batch=25000 -D cores=2 \
  -D output=/tmp/gammaloop-qcd-lightcone/01
```

For the independent repeats at `K^2=0.1,0.01`, run `model near_lightcone nlo
integrate`, with `-D samples=500000 -D batch=50000 -D cores=4 -D seed=271828`
and distinct state/output prefixes. For the `e_cm=0.1` control at `K^2=0.01`,
add `-D sampling_scale=0.1` to the 100,000-sample pilot command.

The nonfinite point can be replayed without an integration:

```sh
./gammaloop --dev-optim -n -L warn \
  -s /tmp/gammaloop-qcd-lightcone/probe-state \
  examples/cli/thermal_production/qcd_vector/run.toml \
  run model near_lightcone nlo lightcone_nonfinite \
  -D kz=3.998749804626441 \
  -D output=/tmp/gammaloop-qcd-lightcone/probe
```

`-n` prevents saving persistent process state; this diagnostic deliberately
allows structured warning logs in the isolated state directory. Do not add
`--read-only-state` when the warning evidence is needed, because it disables
file logging. Build `reference.cpp` as described in the main README, then run:

```sh
/tmp/gammaloop-gauge-reference/reference 4 3.872983346207417 0.00001
/tmp/gammaloop-gauge-reference/reference 4 3.987480407475377 0.00001
/tmp/gammaloop-gauge-reference/reference 4 3.998749804626441 0.00001
```

[lightcone_results.json](lightcone_results.json) contains all ten native
integration records, both reference tolerances, the scale control, and the
nonfinite replay evidence. No core implementation or tests were changed.
