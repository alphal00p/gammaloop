#set page(paper: "a4", margin: 20mm)
#set text(size: 10pt)

= Figure 15 leading-order scan

All 35 points and all 70 physical polarization rows reached the declared
0.5% statistical precision goal. The scan finished in 35.965646 seconds,
inside its three-hour budget. All final components agree with their analytic
targets within 2.421 quoted Monte Carlo errors. Central values are within
0.5% of target for 57 of 70 rows; statistical precision does not guarantee
that every central value falls inside that band.

This new LO measurement uses source
`31e8c8dcc92284cc52d3a36e4a713783ee3dcb4a` and binary SHA256
`fa17c6c38bcc590e456860726a1f0fbbb818d03494bc313497902bcb7da761c9`.
The #link("../baseline.typ")[historical baseline] and earlier NLO, profiling,
and vacuum result directories remain measurements of `116767bd3`.

== Completion, precision, and cost

The 35-point grid covers $k/T=0,0.5,1,1.5$ and
$E/T=0.25,0.5,...,2$, excluding $E=k$, with $E=k plus.minus 0.05$
for each nonzero momentum. The run began on 20 October 2026 at
12:25:16.266038 UTC and finished at 12:25:52.231684. Its hard deadline
was 15:25:16.266038. All 35 jobs ended on precision. The watchdog disarmed
after compute exited; no compute group members remain. Saved-state hashes
and metadata fingerprints were unchanged.

The audit confirms a peak of 60 numerical workers on disjoint CPU sets,
six workers per active point. Controller, signal handlers, gates, and watchdog
are additional control tasks; this is not a bound on all operating-system
tasks. Samples per point range from 100,000 to 256,000, totaling 4.226 million
workspace sampling points and 7.568 million native process evaluations.
The allocation ledger records 1,539.697 worker-seconds (25.662 worker-minutes),
including reserved-session overhead. This is allocation cost, not measured
operating-system CPU time. Build and later rendering are separate; campaign
duration includes its own hashing and launch setup.

The declared rule is
$sigma_(op("MC"))/abs(rho_(op("LO,target"))) <= 0.005$, with at least
100,000 completed samples, zero absolute-error floor, and zero final nonfinite
or unstable evaluations. Both components must pass at a moving point.
Acceptance references are frozen; residuals never control stopping.
Final target-normalized errors range from 0.09554% to 0.49957%; all rows
have zero final nonfinite or unstable samples. Precision-complete CLI shutdown
uses a deliberate signal, so saved exit code −15 is expected.

Rest T/L share one trace estimate and are fully correlated. There are 70
physical rows and 62 distinct native slots: eight rest traces and two direct
components at each of 27 moving points. Moving components also share samples.
Samples are counted once per workspace and process evaluations once per native
slot, avoiding duplicated rest rows.

== Agreement and near-light-cone points

The largest offset is T at $(E/T,k/T)=(0.25,1)$:
$0.19572370 plus.minus 0.00078915$ versus 0.19381362, or +0.9855%
and +2.4204 errors. Its L is
$-0.34842722 plus.minus 0.00143750$ versus −0.34536149, or −0.8877%
and −2.1327 errors. Signed relative differences divide by the absolute target.
These are the only final rows beyond two errors; no row reaches three errors.
This run establishes no persistent target discrepancy. Precision, a central
value within 0.5% of target, and absence of bias are separate claims.

All twelve light-cone-flank components pass precision and lie within
1.474 errors. The largest above-cone offsets are at $(1.55,1.5)$:
T +0.5863% (+1.4733 errors), L +0.6822% (+1.3832 errors).
At $(1.05,1)$, L is −0.5863% (−1.1860 errors).
At $(0.95,1)$, T/L are −0.0907%/+0.0829% (−0.1892/+0.1675 errors).
#link("agreement-audit.json")[The audit] retains every residual and flank comparison.

== Targets, normalization, and figure

Values use $T=1$, $N_c=3$, one massless flavor, and factored-out $e^2 Q_d^2$.
The native model uses $e=1$, $Q_d=-1/3$; direct results use $9E$ times the
native integral, while each physical rest component uses $3E$ times the trace.
No analytic term is added and no NLO $g=3$ multiplier is applied.

The 70 #link("reference/targets.json")[targets] use independent smooth
one-dimensional Fermi phase-space moments equivalent to the LO formulas in
Jackson's #link("https://arxiv.org/abs/1910.07552v2")[paper] and
#link("https://doi.org/10.5281/zenodo.3478144")[pinned author source].
Direct quadrature avoids its truncated trilogarithm; no NLO master is evaluated.
Rest values use $rho_T=rho_L=N_c E^2 tanh(E/4)/(12 pi)$; spacelike T/L are
positive/negative. The earlier NLO $E=3k$ anomaly does not apply to LO.
64/96-node refinement shifts targets by at most $2.73 times 10^(-16)$ relative;
closed-trace and rest-limit checks agree within $3.67 times 10^(-15)$ and
$8.89 times 10^(-16)$. The 64 saved author LO values differ by at most
$4.38 times 10^(-6)$ relative. #link("reference/provenance.json")[Provenance],
#link("reference/validation.json")[validation], and
#link("reference/reference.typ")[the derivation] record sources and checks.

The #link("figure15-lo-overlay.pdf")[full-range overlay] and
#link("figure15-lo-overlay-zoom.pdf")[fixed-range view] use only the six
published vector-extracted LO paths: $k/T=0.5,1,1.5$, in red, purple, and blue.
The extraction has no rest LO path; orange rest markers have no added line.
Both views contain 70 native markers, 1,794 published vertices, and no clipped
native error bars. The eight published LO-plus-NLO curves are excluded.
These graphical paths are distinct from targets used for numerical residuals.

#link("lo-residuals.pdf")[Residuals],
#link("figure15-lo-render.json")[rendering provenance],
#link("summary.json")[unchanged summary], #link("summary.csv")[CSV],
#link("lo-comparison.csv")[analytic comparison],
#link("native-results.json")[final raw objects], and
#link("manifest.json")[compact provenance] preserve measurements.
Complete raw records, executable, states, checkpoints, logs, and executed
controller sources remain in scratch and a local snapshot identified by
#link("raw-run.json")[the raw-run record]. No compressed archive or remote
publication is asserted. Earlier historical numerical results are unchanged.
