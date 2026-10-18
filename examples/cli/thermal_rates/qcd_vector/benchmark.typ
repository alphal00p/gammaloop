= Hard timelike thermal QCD vector current

This example restores the infrared-finite hard timelike comparison with
#link("https://arxiv.org/abs/1910.07552v2")[G. Jackson, _Two-loop thermal
spectral functions with general kinematics_], Eq. (2.6) and Section 7.
The independent NLO reference is the author's
#link("https://doi.org/10.5281/zenodo.3478144")[public implementation,
version v0.9-spf]. This is a virtual-current spectral function at
$T=1$ and $K=(4,0,0,1)$, hence $K^2=15$. Light-cone resummation and
thermal threshold/infrared subtraction are outside this example.

== Shared native calculation

The Standard Model supplies one massless down-quark flavor and its QCD
interactions, with $N_c=3$ and $C_F=4/3$. Setting `aEWM1=4*pi` and
`aS=1/(4*pi)` gives $e=g=1$, while $Q_d=-1/3$.
The incoming photon uses the native Feynman-gauge `summed` projector
$-g^(mu nu)$, global prefactor one, and no polarization average.
The inclusive process is `a > d d~` at `QED^2==2`, with `QCD^2==0`
for LO and `QCD^2==2` for NLO. The latter contains the quark
self-energy supergraph `GL0` and vertex supergraph `GL2`, including
the real, virtual, raised, and crossed thermal cuts selected by the
shared CFF/LU implementation.

The card explicitly selects local three-dimensional UV subtraction with
integrated MSbar counterterms. The renormalization scale is `mu_r=4`;
the auxiliary `m_uv` and renormalization-localization scales are initially
four. The `auxiliary_scales` block changes the latter two to two while
leaving the physical renormalization scale fixed. Such an auxiliary-scale
comparison is meaningful after integration, not pointwise.
Threshold subtraction and soft counterterms remain disabled. Graphs,
orientations, and OSE-weighted sampling channels are summed.
The restored example uses the current evaluator functions; it contains no
historical replacement of hyperbolic functions or alternative thermal
residue engine.

== Observable and independent reference

Define the unit-charge one-flavor current spectral trace as
$rho_V=-op("Im") Pi^mu_mu=rho_V^"LO"+g^2 rho_V^"NLO"+O(g^4)$.
The author's implementation includes $N_c$ and removes the electromagnetic
charge factor. With the native one-particle flux $1/(2k_0)$,

$ Gamma_rho=e^2 Q_d^2 rho_V/k_0, quad rho_V=36 Gamma_rho $

at this benchmark point. NLO denotes the correction alone. A lesser
production kernel is $n_B(k_0) Gamma_rho$; dilepton yields require further
electromagnetic and lepton phase-space factors.
The analytic leading-order control is

$ rho_V^"LO" = (N_c K^2)/(4pi) [1+(2T)/k ln(
  (1+exp(-(k_0+k)/(2T)))/(1+exp(-(k_0-k)/(2T))))]. $

At the specified point it is $2.7034759765132232$.
The vacuum value is $N_c K^2/(4pi)=3.580986219567645$ and the vacuum
NLO/LO ratio at $g=1$ is $3C_F/(16pi^2)=1/(4pi^2)$, giving
$rho_(V,"vac")^"NLO"=0.09070744059337357$.

`reference.cpp` only renames the upstream example entry point, invokes
its existing `Rho_V`, selects the requested absolute tolerance, and prints
the result. It does not change the master's integrands or domains.
The historical reference at tolerances $10^(-4)$ and $10^(-5)$ gives
$rho_V^"NLO"=0.106483279519403$ and $0.10648332415676968$.
Their difference is $4.4637366678 times 10^(-8)$; tolerance variation is
a convergence check, not a rigorous combined error bound.

== Historical evidence and current validation scope

#link("historical_validation.json")[The historical record] is explicitly
separate from validation on the current checkout. Its source implementation
was commit `866387db39eda30bb1bc2ed75c8923586daeeda7`; the preserved
investigation branch tip is `a0b2d4191cad17d8dc2f20891bb0f245f0df18b7`.

#table(
  columns: (1fr, auto, 1fr, auto),
  [Historical calculation], [Samples], [$rho_V$], [Reference pull],
  [Thermal LO], [100,000], [$2.69652944 plus.minus 0.01369827$], [$-0.51$],
  [Thermal NLO], [1,000,000], [$0.10785547 plus.minus 0.00101861$], [$1.35$],
  [Vacuum NLO, original seed], [1,000,000], [$0.09241921 plus.minus 0.00062052$], [$2.76$],
  [Vacuum NLO, independent seed], [2,000,000], [$0.09068381 plus.minus 0.00043810$], [$-0.05$],
  [Thermal NLO, auxiliary scales two], [300,000], [$0.10918683 plus.minus 0.00239492$], [$1.13$],
)

The historical thermal NLO run discarded one nonfinite sample after its
precision ladder; its coordinates were not retained. Consequently these
numbers establish agreement at the chosen hard point with a documented
stability limitation. They do not establish every soft/collinear limit.
The old spacelike NLO comparison was discrepant and a separate repeat
aborted; near-light-cone runs had large errors and reproducible nonfinite
points. Those investigations are not successful precision benchmarks and
are not restored as acceptance targets here.

The current fast regressions in `tests/tests/test_runs/thermal_inclusive_qcd.rs`
compare LO with an independent pointwise Fermi phase-space density and
its integrated analytic value. The NLO regression checks the inclusive
graph inventory, real/virtual and raised/crossed cut support, finite
native values, retained-event closure and conservation, auxiliary-scale
evaluation, and saved-state reload. A passing point or structural test is
not an integrated thermal NLO literature comparison. Record the current
card's completed integrations and their numerical-stability counters
before claiming that comparison has been reproduced.

#link("validation.json")[The current validation record] includes the final
116-test focused run, all passing in 64.499 seconds:
`07460a80-6bc8-4df6-b20a-d00fc0f77ac0`.
The LO extrapolation is $2.7034759917903606$, agreeing with the analytic
reference to $5.66 times 10^(-9)$ relative. The NLO point is finite and
real for both masters and their inclusive sum at auxiliary scales four
and two; its 15 retained events conserve momentum and close to the
returned weight. Reloading and returning the scales to four reproduces
the baseline exactly. The initial 44-of-45 test run and its subsequently
resolved large-root materialization failure remain in the record.

The rebuilt CLI generated, exported, inspected, and saved the example
successfully. Its SHA256 is
`a7752162aa0b917e42fd41522d4752d676d8f20be5a0b3440eb3eafa05083dfc`.
A fresh pilot and three control runs then produced the following results.
Each value is retained separately; the independent repeat does not replace
the pilot.

#table(
  columns: (1fr, auto, 1fr, auto),
  [Contribution / seed], [Samples], [Native $rho_V$], [Reference pull],
  [Thermal LO / 191007552], [50,000], [$2.69117422 plus.minus 0.02009022$], [$-0.61$],
  [Thermal NLO / 191007552], [25,000], [$0.10379516 plus.minus 0.00647305$], [$-0.42$],
  [Thermal NLO / 271828], [50,000], [$0.10896215 plus.minus 0.00393767$], [$0.63$],
  [Thermal NLO, auxiliary scales two / 314159], [25,000], [$0.12158622 plus.minus 0.00642792$], [$2.35$],
  [Vacuum LO / 191007552], [25,000], [$3.56792187 plus.minus 0.03761782$], [$-0.35$],
  [Vacuum NLO / 191007552], [25,000], [$0.09274639 plus.minus 0.00410994$], [$0.50$],
)

The independent thermal NLO result agrees within 0.63 quoted errors,
with a 3.61% Monte Carlo uncertainty. The auxiliary-scale result lies
2.35 errors above the reference and 1.67 combined errors above that
independent run. It remains statistically compatible at this precision;
the visible scatter is not a precise scale-invariance certificate.
The current runs do not reproduce the historical million-sample precision,
and their limits are retained without treating the run-card accuracy goal
as a prerequisite for further development.

All 200,000 samples finished without final nonfinite or unstable discards.
Intermediate thermal NLO nonfinite attempts were recovered by the
precision ladder. The pilot, independent repeat, and auxiliary-scale
control respectively finished 11.692%, 10.308%, and 11.396% at Quad;
the rest finished at Double, with none requiring Arb. The two LO runs
finished entirely at Double. Vacuum NLO required Quad for 0.012%.

The four-worker pilot took 400 seconds including LO. The independent
four-worker NLO run took 841 seconds, the two-worker auxiliary-scale
control 834 seconds, and the two-worker vacuum LO/NLO control 24 seconds.
Thermal NLO sample costs averaged 58.5--64.3 ms, exceeding the 10.6 ms
measured at the fixed inspection point. Random configurations and their
precision escalations must therefore determine larger-run budgets.
The validation JSON retains all raw native results, commands, executable
and artifact hashes, standalone cut inventories, and stability statistics.

== Reproduction

Build the current CLI using the repository workflow, then run from the root:

```sh
./gammaloop -n examples/cli/thermal_rates/qcd_vector/run.toml run model lo nlo export inspect integrate_lo integrate -D output=/tmp/thermal-rates-qcd-thermal
./gammaloop -n examples/cli/thermal_rates/qcd_vector/run.toml run model lo nlo integrate_lo integrate -D mode=vacuum -D output=/tmp/thermal-rates-qcd-vacuum
./gammaloop -n examples/cli/thermal_rates/qcd_vector/run.toml run model nlo auxiliary_scales integrate -D seed=314159 -D samples=300000 -D output=/tmp/thermal-rates-qcd-auxiliary
```

The `lo` and `nlo` blocks generate without integrating. Options `samples`,
`batch`, `lo_samples`, `cores`, and `seed` set the respective run budgets.
Use distinct output prefixes: the integration commands use `--restart`.
Preserve the executable/source identity, complete result JSON, and
discarded-sample counters for every comparison. The target accuracy in
the illustrative run card is not a prerequisite for further development.

To reuse generation across integrations, run `model lo nlo save` once with
an explicit `-s /tmp/thermal-rates-qcd-saved` state. Subsequent invocations
with that state and this card can run `integrate` or
`auxiliary_scales integrate` directly. The explicit `save` block persists
the generated evaluators even when automatic saving is disabled by `-n`.

The independent driver requires a C++11 compiler and GSL. In a scratch
directory, download the exact source archive and verify its hash:

```sh
mkdir -p /tmp/thermal-rates-qcd-reference
curl -fL https://zenodo.org/api/records/3478144/files/gw3g/spectral-v0.9-spf.zip/content -o /tmp/thermal-rates-qcd-reference/spectral.zip
sha256sum /tmp/thermal-rates-qcd-reference/spectral.zip
unzip /tmp/thermal-rates-qcd-reference/spectral.zip -d /tmp/thermal-rates-qcd-reference
spectral_src=/tmp/thermal-rates-qcd-reference/gw3g-spectral-dca510d
c++ -std=c++11 -O2 $(gsl-config --cflags) -I "$spectral_src" -I "$spectral_src/inc" "$spectral_src"/src/*.cpp examples/cli/thermal_rates/qcd_vector/reference.cpp $(gsl-config --libs) -o /tmp/thermal-rates-qcd-reference/reference
/tmp/thermal-rates-qcd-reference/reference 4 1 0.0001
/tmp/thermal-rates-qcd-reference/reference 4 1 0.00001
```

The archive SHA256 is
`207a2cb3d39933a57a5c1a8592d664c7dc9d5750eb8cd10801e98ac9560446fb`;
its source commit is `dca510d2d2c8bee2288345e972ee3a54b37cdbd5`.
The historical evaluation used GSL 2.8. The upstream example's renamed
entry point emits a missing-return warning, but that function is not
invoked by this driver.
