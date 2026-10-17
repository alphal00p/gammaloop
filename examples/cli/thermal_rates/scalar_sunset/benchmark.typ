= Scalar sunset acceptance example

This example defines four distinct real scalars with interaction
$cal(L)_("int") = -lambda X B_1 B_2 B_3$.
The probe mass is $M=1$, all three bath masses and chemical potentials vanish,
and $lambda=T=1$. The portable JSON model is imported directly; no UFO Python
package is required. Generating `xs X > B1 B2 B3 FI^2==2 [{{2}}]` gives the
two-loop forward sunset and its three-particle physical cut support.

The native observable is the spectral relaxation kernel
$Gamma_rho(P)=-op("Im") Pi_R(P)/E$, with $E=sqrt(M^2+p^2)$.
The native one-particle flux already supplies $1/(2E)$.
For an initially empty bosonic probe, KMS gives
$S(P)=n_B(E) Gamma_rho(P)$ and
$dif gamma / dif^3 p = S(P)/(2pi)^3$.
The Bose factor in this conversion does not assign an equilibrium occupation
to the probe. The forward graph has only bath fields on its internal lines.

The `thermal_excess` process enables vacuum subtraction in the existing thermal weights
and returns $Gamma_rho - Gamma_("vac")$. The production-source formula above requires
the full `equilibrium` kernel; add the vacuum contribution back before applying KMS
if starting from `thermal_excess`.

== Independent reference

`reference.py` derives the answer directly from Lorentz-invariant phase space.
It factorizes one bath momentum $k$ and a massless two-particle pair $Q$, then
integrates the latter pair's angle analytically. Write $w=Q^0$, $q=|bold(Q)|$,
and $n(x)=(exp(x/T)-1)^(-1)$. Its angular occupation average is

$ A(w,q) = n(w) [1 + (2T)/q ln(
  (1-exp(-(w+q)/(2T))) / (1-exp(-(w-q)/(2T))))]. $

The continuous $q=0$ limit is $n(w/2)^2$. With $c$ the angle between
$bold(p)$ and the singled-out bath momentum, the source contributions are

$ S_("inverse decay") =
  lambda^2/(128pi^3 E) integral_(-1)^1 dif c
  integral_0^(M^2/(2(E-p c))) dif k k n(k)
  A(E-k,sqrt(p^2+k^2-2p k c)), $

$ S_("scattering") =
  (3lambda^2)/(128pi^3 E) integral_(-1)^1 dif c
  integral_0^infinity dif k k (1+n(k))
  A(E+k,sqrt(p^2+k^2+2p k c)). $

The factor three counts the distinct choices of the crossed bath species.
There is no identical-particle factorial. The vacuum control is exactly
$Gamma_("vac")=lambda^2 M^2/(512pi^3 E)$.
The Maxwell--Boltzmann controls printed by the script are separate dilute
limits, not substitutes for the Bose reference.

The script uses Gauss--Legendre quadrature in the bath frame. A fourth-power
coordinate smooths the inverse-decay endpoint; a rational coordinate covers
the scattering half-line. The output retains all requested quadrature orders,
both channel contributions, normalization, and physical parameters. Differences
between orders estimate convergence; they are not certified error bounds.

== Native validation

From the repository root, generate and inspect at rest:

```sh
./gammaloop -n examples/cli/thermal_rates/scalar_sunset/run.toml run generate rest inspect
python3 examples/cli/thermal_rates/scalar_sunset/reference.py --orders 64 128 256 512 --output /tmp/thermal-rates-sunset-reference.json
```

Replace `rest` by `boosted` or `fast` for $p=1$ or $p=3$.
The `integrate` block integrates `equilibrium`, `vacuum`, and `thermal_excess` together.
Compare the first two with `spectral_rate` and `vacuum_spectral_rate` in the reference
output, and the third with their difference. The vacuum-subtracted result must also
reproduce that difference pointwise and after event selection. Use a separate integration
workspace for each momentum; the example's `--restart` deliberately replaces
its scratch integration. Save the card, output, sample count, seed, and numerical
stability statistics for every comparison.

The integration regression in `tests/tests/test_runs/thermal_rates.rs` exercises
the same card. It checks rest-frame phase-space densities against explicit
occupation factors and the normalized radial test function, all three scattering
channels plus decay, positive-energy events, four-momentum conservation,
event-weight closure, boosted kinematics, and the cold finite-temperature limit.
An outgoing-B1 selector compares the selected result with the corresponding retained-event
subset in both equilibrium and vacuum, including evaluation with event retention disabled.
It does not establish support for massive bath lines, raised cuts, nonzero
chemical potentials, UV or threshold counterterms, or resummed rates.

== Initial milestone verification record

The independent oracle was evaluated through quadrature order 512. Its sources
at $p=0,1,3$ are $0.004229085051719782$, $0.001501240527174083$, and
$0.00006989233735691003$. Orders 256 and 512 differ by at most
$4.2 times 10^(-18)$ for these sources.

Native point checks use loop momenta $(0.6,0,0)$ and $(0,0.8,0)$. Every boost
retains four equilibrium events (decay and the three crossed bath species) and
one vacuum event. Their weights reproduce the integrand, with maximum
four-momentum residual $1.34 times 10^(-15)$. The cold finite-temperature
control, $beta=1000$, equals the vacuum result.

Every signed event was also checked against direct phase space. For these rays,
independently solve
$F_sigma(t)=(0.6sigma_1+0.8sigma_2)t
 +sigma_3 sqrt(t^2+p^2)-E=0$ by bisection.
The energy magnitudes are $(0.6t,0.8t,sqrt(t^2+p^2))$ and the cut weight is
the forward-minus-reverse Bose product. The resulting density is
$(2pi t^6 h(t) W_sigma)/
 (abs(F'_sigma(t)) (product_i 2E_i) (2pi)^6 2E)$,
where $h(t)=2/sqrt(pi) exp(2-t^2-t^(-2))$.
This comparison uses neither native roots nor generated residues.
All twelve signed-event weights agree within $3.4 times 10^(-13)$ relative.

The integrated results below are $Gamma_rho$; errors are the quoted Monte Carlo
estimates. Sample counts apply to each of the two correlated slots.

#table(
  columns: (auto, auto, 1fr, 1fr),
  [Momentum], [Samples / seed], [Equilibrium], [Vacuum],
  [$p=0$], [50,000 / 1701026],
  [$0.00718607 plus.minus 0.00016511$], [$6.22805 times 10^(-5) plus.minus 1.48444 times 10^(-6)$],
  [$p=1$], [100,000 / 1701026],
  [$0.00456812 plus.minus 0.00013077$], [$4.50695 times 10^(-5) plus.minus 1.88866 times 10^(-6)$],
  [$p=3$ pilot], [100,000 / 1701026],
  [$0.00129701 plus.minus 0.00005061$], [$1.44446 times 10^(-5) plus.minus 1.72838 times 10^(-6)$],
  [$p=3$ repeat], [1,000,000 / 1701027],
  [$0.00158991 plus.minus 0.00004526$], [$1.95674 times 10^(-5) plus.minus 5.59874 times 10^(-7)$],
)

The rest and $p=1$ results agree within 0.81 quoted errors. The first $p=3$
pilot was low by 5.62 errors thermally and 3.17 errors in vacuum; it is retained
in the record. An independent-seed adaptive repeat agrees within 0.19 and
0.63 errors, respectively. Its thermal and vacuum uncertainties are 2.85%
and 2.86%: the historical 0.5% target was not reached before the sample cap.
That target is not a prerequisite for validating additional infrared-finite examples.

The repeat's iteration record shows the central values moving towards the
references as angular coverage improves. Its thermal uncertainty increases
when larger contributions appear late. The successful signed point checks and
this convergence evidence support a sampling explanation for the pilot, without
establishing general sampling efficiency or rigorous error coverage.

All four integrations, totaling 2.5 million integrand evaluations, reported
zero final discarded nonfinite or unstable samples. Intermediate nonfinite
evaluations did occur and were recovered by the precision ladder. The
million-sample repeat used double precision for 98.5105%, quadruple for 1.4894%,
and arbitrary precision for 0.0001% of evaluations (two evaluations).
Real-zero counts include ordinary vanishing or underflowed weights and are
recorded separately from final failures.

#link("validation.json")[The compact validation record] preserves commands,
physical parameters, settings hashes and paths, native/reference values,
precision diagnostics, and the repeat's iteration evolution. Raw inspection
JSON, event data, integration results, and logs remain at the recorded scratch
paths. The executable SHA256 is
`4036e004f2ed2522b4d59a09bdeb999e1c9714e885e8db634ef68d3501598fca`.
It predates the final structural preflight guards; these integrations are
evidence for the numerical implementation in that executable.

At the initial milestone, the focused source-test selection passed all 14 tests across
four binaries, including the thermal saved-state round trip and existing
raised-E regressions. The nextest run was
`6f6c60c1-6c35-4214-b5bd-ac55fa148f6d`, using `ci-optim` and
`ci_gammaloop`; its log is `/tmp/thermal-rates-tests.log`.
This was a focused selection, not the full repository test suite.

The initial milestone's final CLI build (SHA256 `9bfb0fe2a556ad14ca6ba58ec419727af34990d5ccc3bb32b8795e3a222a5916`) also passed the rest-frame generation/inspection smoke and rejected enabled threshold subtraction and a spacelike probe at their explicit support guards; commands and logs are recorded separately in `validation.json`, preserving the earlier numerical executable's provenance.

== Shared-LU refactor verification

The subsequent shared E/H surface-group and radial-equation refactor passed 195 tests
across five binaries (nextest `e8fd0805-f408-427a-bc7d-6c883d997ae4`,
`ci-optim` / `ci_gammaloop`; 2142 tests skipped). This selection includes
the saved-state round trip and outgoing-B1 selector with retained and hidden events.
Clippy passed without warnings, and formatting checks passed for the changed Rust sources.
Logs and exact commands are recorded separately in `validation.json`.
The rebuilt CLI passed fresh rest-frame generation and inspection, reproducing the original
thermal and vacuum point values exactly, with four and one events, exact event-weight closure,
and four-momentum residuals below $2.23 times 10^(-16)$.
Its SHA256 is `cba9ab0c65ae1e22da98b4c96f35fb2eeac6b76f4453d49a2e090cde1700fa97`.
The CLI smoke did not persist a state. The regression exercised the saved-state round trip
with the temporary manifest number 17; the unchanged layout was subsequently assigned
version 16, consolidating the uncommitted changes into one version bump from 15.
The historical integrations above were not rerun during that refactor; its supported
thermal physics scope was unchanged. These records predate the vacuum-subtracted extension.

== Vacuum-subtraction verification

The continuation started with 195 passing baseline tests on clean commit `3e73e6b07`.
After enabling vacuum subtraction through the existing shared thermal weights, all 12
targeted tests passed. The scalar sunset regression checks
$Gamma_("excess") = Gamma_rho - Gamma_("vac")$ both in total and per physical event
at $p=0,1,3$, including saved-state reload and the existing selector checks.
An independent occupation polynomial also verifies a positive excess of approximately
$1.7071821495 times 10^(-30)$ at $beta=200$. The vacuum constant is kept separate in
this oracle, so the comparison remains meaningful below binary64's relative resolution
for the full kernel.

The fresh CLI rest smoke gives `thermal_excess = 5.248854262631212e-05` with four
events. Total and per-event equilibrium-minus-vacuum identities and event closure hold
exactly; conservation residuals are at most $2.23 times 10^(-16)$.
Clippy, Rust formatting, and the edited documentation renders passed. Commands,
source/executable hashes, run IDs, and local log paths are retained under
`vacuum_subtraction_continuation` in #link("validation.json")[the validation record].
These are focused checks; the historical integrations were not rerun.
