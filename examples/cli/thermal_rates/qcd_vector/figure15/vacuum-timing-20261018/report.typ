= Matched vacuum and thermal QCD NLO timing

The small comparison uses source commit
`116767bd3ed30491a2b15559d228171d20082abe` and the same immutable CLI
as the Figure 15 scan. Four jobs evaluate transverse and longitudinal
components at $(E,k)=(1,0.5)$ and $(1.05,1)$, using 1,000 common random
samples per component and seed `2026101819`. Each job also benchmarks
three fixed six-dimensional sampling points for both components, with
ten warmups and eight timing batches per point. All four jobs finish
successfully with zero unresolved or nonfinite samples.

The model files are byte-identical. Effective runtime settings after the
explicit overrides are identical between thermal and vacuum jobs, including
the projectors, kinematics, sampling, local three-dimensional UV subtraction,
integrated MSbar counterterms, scales, evaluator method, and default
stability checks. The vacuum state is regenerated with medium mode `vacuum`;
the thermal state uses `thermodynamic_equilibrium`. Artifact paths differ.
No Rust code or production-scan physics settings are changed.

Measured mean worker time per integrand evaluation:

#table(
  columns: 5,
  [$(E,k)$], [Category], [Thermal], [Vacuum], [Ratio],
  [$(1,0.5)$], [Cut/kinematic core], [60.39 ms], [1.368 ms], [44.2],
  [$(1,0.5)$], [Total], [63.78 ms], [2.123 ms], [30.0],
  [$(1.05,1)$], [Cut/kinematic core], [110.74 ms], [1.349 ms], [82.1],
  [$(1.05,1)$], [Total], [118.98 ms], [2.108 ms], [56.4],
)

Here core means the saved inclusive integrand time minus expression-evaluator
and observable time. T and L share the random-batch timing counter, so their
identical timing rows are not independent measurements. One joint T/L sample
requires two evaluations. The fixed-point benchmarks corroborate the core
ratios: 42.1 at the ordinary point and 85.9 near the light cone, using an
equal-weight mean over the three coordinates and two components. Near the
light cone, thermal fixed-point total times range from 57.7 to 145.5 ms;
the corresponding vacuum times are 1.86 to 1.97 ms.

The thermal quad-precision fractions are 12.6% and 74.2%, compared with
0.10% and 0.15% in vacuum. These percentages count samples, not time.
The measured overhead includes natural precision escalation and stability
probes. This test does not isolate the cost of occupation functions at
fixed precision and fixed residue content.

Both modes contain master graphs `GL0` and `GL2` and the same ordinary cut
counts, three and four. Thermal orientation counts are 32 and 32, versus
14 and 18 in vacuum. Generated evaluator counts are 15 and 21 thermally,
versus five and five in vacuum. Medium generation also retains additional
H-surface bath-crossing residue groups. Thus the first pipeline difference
arises during medium-specific generation, before physical evaluation.
The broad core timer cannot distinguish root solving, parameter construction,
residue assembly, and their precision retries.

#link("comparison.json")[The full report] preserves timing batches and the
configuration audit; #link("summary.csv")[the CSV] contains the main ratios.
#link("manifest.json")[The original manifest] retains all commands and point
coordinates. The command driver, per-job integration results, fixed-point
benchmark batches, and cut display logs are preserved in the verified
#link("../archive.json")[archive], under the original checkout's
`vacuum-timing-20261018/` and `tmp/qcd-vacuum-timing-20261018T191055Z/` paths.
Generation time is excluded from the reported evaluation timings. These are
historical costs for `116767bd3`; the later `75a0ced4e` optimization was not
measured here.

The user authorized up to six additional workers for this test. Those jobs
have ended. The main scan subsequently finished within its 60-worker ceiling
and original hard cutoff of 19 October 2026 at 07:19:29 Europe/Zurich.
