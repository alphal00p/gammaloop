# QCD vector-current spectral function

This example compares GammaLoop's native equilibrium cross-section integrand
with the QCD vector channel in G. Jackson,
[arXiv:1910.07552v2](https://arxiv.org/abs/1910.07552v2), specifically Eq. (2.6)
and Section 7. The reference is evaluated using the author's
[public implementation, v0.9-spf](https://doi.org/10.5281/zenodo.3478144), rather
than digitizing Figure 15.

The benchmark has one massless quark flavour, `Nc=3`, `CF=4/3`, zero chemical
potentials, `T=1`, and `K=(4,0,0,1)`, hence `K^2=15`. It tests a timelike virtual
photon, not real-photon emission on the light cone. The paper's discussion below
Eq. (1.1) identifies `|K^2|` of order `g^2 T^2` as a region requiring
resummation; this benchmark stays away from it.

## Native calculation and normalization

The existing `sm-default.json` model supplies a down quark and gluon. Setting
`aEWM1=4*pi` and `aS=1/(4*pi)` gives `e=g=1`, while the down-quark charge is
`Qd=-1/3`. No custom model or bath-averaging script is needed.

The external photon uses the `summed` Feynman-gauge projector `-g^{mu nu}` and
global prefactor 1, with no polarization average. Generation selects the
inclusive `a > d d~` channel at orders `QED^2=2` and `QCD^2=0,2`. At NLO it
contains the quark self-energy insertions (`GL0`, including both placements)
and the vertex graph (`GL2`). Equilibrium LU residues include the thermal
crossings as well as the vacuum real and virtual cuts.

Define the unit-charge, one-flavour current spectral function by

```
rho_V = -Im Pi^mu_mu = rho_V_LO + g^2 rho_V_NLO + O(g^4).
```

The author's `Rho_V` includes `Nc=3`, but removes the electromagnetic charge
factor. GammaLoop's one-incoming-particle flux is `1/(2*k0)`; its native
cross-section integral is therefore

```
Gamma_rho = e^2 Qd^2 rho_V / k0,
rho_V = 36 * Gamma_rho                     (this benchmark).
```

Here NLO always means the order-`g^2` correction alone. The conversion is fixed
before comparison, and is checked against the analytic leading-order result:

```
rho_V_LO = Nc*K^2/(4*pi) * [1 + 2*T/k * log(
    (1 + exp(-(k0+k)/(2*T))) / (1 + exp(-(k0-k)/(2*T))))].
```

At this point it is `2.7034759765132232`. In vacuum it is
`Nc*K^2/(4*pi)`, with NLO/LO `3*CF*g^2/(16*pi^2) = 1/(4*pi^2)`.
For the corresponding lesser/production kernel, multiply `Gamma_rho` by
`n_B(k0)`. Conversion to a dilepton yield also requires the electromagnetic
and lepton phase-space factors; those are not part of this benchmark.

The card uses local three-dimensional UV subtraction plus integrated MSbar
counterterms, `mu_r=4`, `m_uv=4`, and localization scale 4. It keeps
`medium.vacuum_subtraction=false`. Both soft counterterms and threshold
subtraction are disabled. Graphs, orientations, and OSE-weighted loop-momentum
basis channels are summed. Evaluators use eager `SummedFunctionMap` and the
existing temporary hyperbolic-function workaround for raised residues.

Structural exports contain:

| Mode | Graph | Orientations | Cut groups | Order-1 / order-2 residue entries |
| --- | --- | ---: | ---: | ---: |
| vacuum | GL0 | 14 | 2 | 2 / 1 |
| vacuum | GL2 | 18 | 4 | 4 / 0 |
| equilibrium | GL0 | 32 | 10 | 10 / 3 |
| equilibrium | GL2 | 32 | 20 | 20 / 0 |

These are internal group/entry counts, not counts of independent scattering
processes. There are no generated threshold counterterms.

## Results

All entries are `rho_V/T^2`, for `g=1`. Uncertainties are the native one-sigma
Monte Carlo estimates, and NLO denotes only the correction.

| Calculation | Samples | GammaLoop | Reference | Difference / MC error |
| --- | ---: | ---: | ---: | ---: |
| Thermal LO | 100,000 | 2.69652944 ± 0.01369827 | 2.70347598 | -0.51 |
| Thermal NLO | 1,000,000 | 0.10785547 ± 0.00101861 | 0.10648332 | +1.35 |
| Vacuum LO | 100,000 | 3.57184470 ± 0.01816606 | 3.58098622 | -0.50 |
| Vacuum NLO, seed 191007552 | 1,000,000 | 0.09241921 ± 0.00062052 | 0.09070744 | +2.76 |
| Vacuum NLO, seed 271828 | 2,000,000 | 0.09068381 ± 0.00043810 | 0.09070744 | -0.05 |
| Thermal NLO, auxiliary scales 2 | 300,000 | 0.10918683 ± 0.00239492 | 0.10648332 | +1.13 |

The thermal NLO result agrees with the literature within 1.35 quoted errors.
The auxiliary UV-scale control changes `m_uv` and the localization scale from
4 to 2 while keeping `mu_r=4`; it is consistent with both the baseline and the
reference. The author-code tolerance change from `1e-4` to `1e-5` changes the
NLO reference by `4.46e-8`, well below the native integration errors.

The first vacuum NLO seed lies 2.76 quoted errors above its analytic target;
the independent two-million-sample run differs by only 0.05 errors. Both runs
are retained here rather than treating the second as a replacement. The
baseline thermal NLO relative error is 0.94%.

One sample out of the million-sample thermal NLO run remained nonfinite after
Double/Quad/Arb precision escalation and was discarded. The other five runs
in this table had no nonfinite or unstable samples. Thus this is numerical
agreement with a rare unresolved stability failure, not a demonstration that
all soft/collinear configurations are numerically safe. The saved terminal
warning identifies the final nonfinite evaluation boundary but does not
retain the sample coordinates; no particular physical origin is established.
A replay with structured warning logs is needed to isolate that sample.

The strict massless gauge-theory sum went through using the existing core
implementation and the previously added hyperbolic-function workaround. No
new threshold, soft, or collinear counterterm or resummation was implemented.
This conclusion is restricted to the chosen hard timelike kinematics.

[`results.json`](results.json) retains the complete native result records,
reference values, normalization, seeds, sample counts, and stability statistics.

## Reproduction

Run from the repository root using the existing dev-optim CLI:

```sh
./gammaloop --dev-optim --read-only-state \
  -s /tmp/gammaloop-qcd-vector-thermal-final-state -n -L off \
  examples/cli/thermal_production/qcd_vector/run.toml \
  run model lo nlo export inspect integrate \
  -D output=/tmp/gammaloop-qcd-vector-thermal-final

./gammaloop --dev-optim --read-only-state \
  -s /tmp/gammaloop-qcd-vector-vacuum-final-state -n -L off \
  examples/cli/thermal_production/qcd_vector/run.toml \
  run model lo nlo integrate -D mode=vacuum \
  -D output=/tmp/gammaloop-qcd-vector-vacuum-final

./gammaloop --dev-optim --read-only-state \
  -s /tmp/gammaloop-qcd-vector-vacuum-seed-state -n -L off \
  examples/cli/thermal_production/qcd_vector/run.toml \
  run model nlo integrate -D mode=vacuum -D seed=271828 \
  -D samples=2000000 -D batch=200000 \
  -D output=/tmp/gammaloop-qcd-vector-vacuum-seed

./gammaloop --dev-optim --read-only-state \
  -s /tmp/gammaloop-qcd-vector-thermal-uv-state -n -L off \
  examples/cli/thermal_production/qcd_vector/run.toml \
  run model nlo integrate -D m_uv=2.0 -D localization=2.0 \
  -D seed=314159 -D samples=300000 \
  -D output=/tmp/gammaloop-qcd-vector-thermal-uv
```

`run.toml` also accepts `samples`, `batch`, `lo_samples`, `cores`, and `seed`
macros. Runtime scale changes leave the physical MSbar scale `mu_r=4` fixed.
The `export` block saves standalone JSON and DOT files under the output prefix.
Integration results are written to `<output>-{lo,nlo}-integral/integration_result.json`.
No persistent GammaLoop state is modified.

The reference requires a C++11 compiler and GSL (2.8 was used here). Download
and build the public source, using `gsl-config` from the installed GSL prefix:

```sh
mkdir -p /tmp/gammaloop-gauge-reference
curl -fL \
  https://zenodo.org/api/records/3478144/files/gw3g/spectral-v0.9-spf.zip/content \
  -o /tmp/gammaloop-gauge-reference/spectral.zip
unzip /tmp/gammaloop-gauge-reference/spectral.zip \
  -d /tmp/gammaloop-gauge-reference
spectral_src=/tmp/gammaloop-gauge-reference/gw3g-spectral-dca510d
c++ -std=c++11 -O2 $(gsl-config --cflags) \
  -I "$spectral_src" -I "$spectral_src/inc" "$spectral_src"/src/*.cpp \
  examples/cli/thermal_production/qcd_vector/reference.cpp \
  $(gsl-config --libs) -o /tmp/gammaloop-gauge-reference/reference
/tmp/gammaloop-gauge-reference/reference 4 1 0.0001
/tmp/gammaloop-gauge-reference/reference 4 1 0.00001
```

`reference.cpp` reuses the author's `Rho_V` without changing its formulas,
master integrals, or integration domains. It selects the vector trace, sets
its existing absolute quadrature tolerance, and prints more digits. The
original example's two unrelated compiler warnings (renamed entry point
without return; deprecated `sprintf`) do not affect the invoked `Rho_V` path.
The author's code suppresses GSL errors and does not return a combined error
estimate, so tolerance variation is a convergence check, not a rigorous error
bound. The reference archive has SHA256
`207a2cb3d39933a57a5c1a8592d664c7dc9d5750eb8cd10801e98ac9560446fb`
and source commit `dca510d2d2c8bee2288345e972ee3a54b37cdbd5`.
