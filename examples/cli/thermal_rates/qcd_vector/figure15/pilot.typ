= Figure 15 pilot, 18 October 2026

The source baseline is `116767bd3ed30491a2b15559d228171d20082abe`.
The native executable, seeds, commands, sample histories, reference values,
and measured timing records are preserved in
#link("pilot-20261018/pilot-report.json")[the pilot report] and
#link("pilot-20261018/pilot-analysis.json")[the detailed analysis].
The worker ceiling was 60 throughout this investigation. These pilots and
controller checks precede the separately timed twelve-hour production scan.

#image("pilot-20261018/pilot-comparison.png", width: 100%)

The plot compares the isolated unit-$g^2$ NLO coefficients, with one quoted
Monte Carlo error per bar. Figure 15 itself uses the sum
$rho^"LO"+9rho^"NLO"$ at $g=3$. The conversion, current projections,
and reference sign convention are documented in
#link("comparison.typ")[the comparison run card].
The two rest entries use the same trace estimate and are fully correlated.

== Pilot results

Seven nonzero-momentum points used correlated direct transverse and
longitudinal integrations, with a cap of 10,000 samples per component and
900 seconds per point. A separate 25,000-sample rest run used an independent
seed to check the low estimate from the original 10,000-sample trace pilot.
The extra hard control $(E/T,k/T)=(4,1)$ lies outside the plotted energy
range of Figure 15. Default stability rotations and the precision ladder
were retained, and no completed sample was unresolved or nonfinite.

The direct LO longitudinal projections passed three independent phase-space
controls within 0.64 quoted errors. At NLO, the completed timelike points
are compatible with the author implementation at their current precision.
The longitudinal point $(1.05,1)$ is 2.37 quoted errors high and remains
a convergence watch point. No persistent discrepancy has been established:
these short pilots have few adaptive iterations and cannot exclude rare
tails or a smaller systematic difference.

None of the pilot polarizations reaches 0.5% relative statistical error.
The independent rest run has 2.36% error. Ordinary timelike direct-component
errors are about 6–18%; near-light-cone and spacelike components are much
noisier. The earlier rest offset decreases from 2.66 to 1.25 quoted errors
with the independent seed; both estimates remain in the preserved records.

== Cost and production policy

The planning estimate for 35 points is approximately 170 allocated worker
hours for 100,000 samples per component, 420 hours for 250,000, and 1,670
hours for one million. At 60 workers these correspond to about 3, 7, and
28 wall hours. These are conditional estimates obtained from nearby pilot
timings and measured worker-time efficiency. The initial adaptive grids
are not trained; a factor-two workload contingency is appropriate, and
the host is shared with other work.

At the pilot variance, the ordinary point $(1,0.5)$ alone needs roughly
1.8 million samples and 70 allocated worker hours for both components to
reach 0.5%. The spacelike and light-cone pilots give much larger,
unreliable extrapolations. Uniform 0.5% precision over the full grid is
therefore not expected within the authorized 720 worker-hour ceiling.

Production first gives every point time, then continues toward 100,000,
250,000, and larger sample checkpoints. It reassigns six-worker slots
after a completed native checkpoint. Statistical precision and residual
against the independent target are recorded separately; agreement never
controls stopping. The relative goal is 0.5% on each isolated NLO
coefficient, with no absolute-error floor. Exact light-cone points and
$E=0$ are excluded. Each point needs at least 100,000 samples before
precision acceptance.

The production limit is 60 workers including author-reference calculations,
with a hard twelve-hour cutoff from the recorded production start.
An independent watchdog guards the entire compute process group.
Final completed-iteration results and logs are retained if the deadline
arrives before the precision goal. No automatic chat completion
notification is configured; the controller writes local JSON and CSV reports.
