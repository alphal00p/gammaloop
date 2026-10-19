#set page(paper: "a4", margin: 20mm)
#set text(size: 10pt)

= Figure 15 full-scan outcome

This is a historical result for source `116767bd3`, before the radial-root
optimization in `75a0ced4e`. See #link("../baseline.typ")[the baseline and
archive record]; the commit containing this report is not the measured revision.

The scan finished within its budget, with results for all 35 kinematic points.
It did not reach 0.5% statistical precision across the full grid. Seven
rest-frame points passed the declared reference-normalized precision goal
for both polarizations. Moving-frame and near-light-cone points generally
remain limited by Monte Carlo variance. The earlier near-light-cone
longitudinal offset did not persist.

There is one major comparison anomaly at $(E/T,k/T)=(1.5,0.5)$, exactly
$E=3k$. Its saved author reference is suspect because strict domain branches
omit equality. This point is flagged for reference validation and excluded
from the ordinary agreement summary below. A corrected reference has not
been computed, so this investigation does not establish the causal size of
that omission or resolve the point's physics.

== Run completion and coverage

The immutable executable uses commit
`116767bd3ed30491a2b15559d228171d20082abe`. The scan began on 18 October
2026 at 19:19:29 Europe/Zurich and all compute exited on 19 October at
07:17:30, after 11 h 58 min 1 s. The original hard deadline was
07:19:29. Graceful interruption started two minutes before that deadline;
the watchdog disarmed only after all compute had exited. No scan jobs remain.

Main allocations, including the author-reference workers, peaked at 60
workers. The six reference workers initially left 54 native workers; their
completion released those CPUs to native integration. Separately authorized
vacuum/profiling probes are outside this campaign accounting. No source
changes or saved-state fingerprint changes were made to the production run.

The grid covers $k/T=0,0.5,1,1.5$ and $E/T=0.25,0.5,...,2$, excluding
the exact light cone, with additional $E=k plus.minus 0.05$ points. Every
point has valid final checkpoints. Samples per workspace range from 166,000
to 1,584,000; the total is 15.716 million workspace sampling points and
25.912 million native process evaluations. Rest-frame isotropy derives
$T=L$ from the same trace integral; those are not independent measurements.
All 70 physical component rows have zero final unresolved/nonfinite samples.
This numerical-stability check does not establish Monte Carlo convergence.

== Precision reached

The controller's acceptance rule is
$sigma_(op("MC"))/abs(rho_(op("reference"))) <= 0.005$ with at least
100,000 samples and zero final unresolved evaluations. Agreement of the
central value with the reference is a separate check. A 0.5% statistical
error does not guarantee a central value within 0.5% of the target.

#table(
  columns: 4,
  [$E/T$, $k=0$], [Samples], [MC error / target], [Result],
  [0.25], [200,000], [0.4998%], [Passed, T=L],
  [0.50], [244,000], [0.4988%], [Passed, T=L],
  [0.75], [334,000], [0.4989%], [Passed, T=L],
  [1.00], [462,000], [0.4994%], [Passed, T=L],
  [1.25], [642,000], [0.4998%], [Passed, T=L],
  [1.50], [882,000], [0.4997%], [Passed, T=L],
  [1.75], [1,172,000], [0.4996%], [Passed, T=L],
  [2.00], [1,584,000], [0.4797% of estimate], [Reference unavailable],
)

The seven passed rest-frame results agree with their targets within 1.18
quoted Monte Carlo errors. The eighth rest point has sufficient provisional
precision relative to its own estimate, but no independent target was
obtained within the reference timeout.

The unmodified controller reports 15 accepted component flags: fourteen
rest-frame rows plus transverse at $(1.5,0.5)$. That transverse flag uses
the suspect reference as its denominator: its MC error is 0.413% of that
reference but 0.781% of its own estimate. It is not counted here as a
successful 0.5% validation. The job at that point remained budget-stopped
because its longitudinal component did not pass.

Some rest-frame errors are marginally above 0.5% when divided by their
current estimate rather than the predeclared reference. Both ratios are
retained in the original summary. Exactly ten component rows have strict
error/estimate at or below 0.5%; eight are four passed rest integrations
and two are the targetless $E/T=2$ rest result. This distinction reflects
the denominator definition, not a change to results or stopping policy.

Away from rest, even ordinary timelike points mostly remain above target.
For example, $(1,0.5)$ ends at 498,000 samples with T/L errors of
0.727%/0.937% of their targets. At $(1.25,1)$, 408,000 samples give
1.33%/2.21%. No moving-frame point passes for both polarizations.

== Near-light-cone longitudinal result

Values below are isolated NLO coefficients including vacuum and thermal
contributions, with unit couplings and the electromagnetic charge factored
out; errors are one Monte Carlo standard
error, without a propagated author-code quadrature uncertainty.

#table(
  columns: 6,
  [$(E/T,k/T)$], [Samples], [Native L], [Reference L], [Pull], [Error / target],
  [$(1.05,1)$], [250,000], [$0.00949367 plus.minus 0.00191959$], [0.01211950], [−1.37], [15.84%],
  [$(0.95,1)$], [170,000], [$0.13832947 plus.minus 0.10321419$], [−0.01059693], [+1.44], [974%],
)

The above-cone transverse result is
$0.04115141 plus.minus 0.00879869$ versus 0.05320425, also −1.37 errors
from target, with 16.54% target-normalized uncertainty. The original
longitudinal pilot offset of +2.37 errors is absent in this scan. Across
125 above-cone saved checkpoints, neither component exceeds two errors
from target; longitudinal stays between roughly −1.72 and −0.76 after
100,000 samples. Neither component below the cone exceeds 1.45 errors.

These are compatibility statements under broad errors, not proof of
convergence. Large individual weights still matter: longitudinal max-weight
impact is 14.9% above and 68.8% below the cone. Rare contributions increased
estimated uncertainty during the ladder, so smooth $1/sqrt(N)$ forecasts
are not reliable for these points. A chi-square near one does not remove
this limitation.

== Reference anomaly and remaining agreement

At $(1.5,0.5)$, the raw saved values are:

#table(
  columns: 4,
  [Component], [Native NLO, 528,000 samples], [Author reference], [Pull],
  [T], [$0.02130332 plus.minus 0.00016631$], [0.04030709], [−114.3],
  [L], [$0.02044940 plus.minus 0.00018202$], [−0.00399609], [+134.3],
)

The author values were freshly computed with exact inputs and the recorded
binary hash; this is not a stale-cache comparison. In the pinned upstream
`src/6_.cpp`, lines 598 and 685 use `k0 > 3*k` and `k0 < 3*k`, with neither
branch handling equality. The limiting domain is not generally empty at
this equality. Neighboring grid points at $E/T=1.25$ and 1.75 have native
and author agreement within 2.15 errors. This is strong evidence for a
reference-boundary problem, but a flank/limit evaluation or corrected
reference is still required before attributing the entire discrepancy to it.

Excluding those two flagged component rows, all 62 available reference
comparisons lie within 2.79 quoted Monte Carlo errors. Many of those errors
are large, especially in spacelike regions, so this does not establish
sub-percent accuracy there.

Independent T/L targets are missing for $(0.25,1)$, $(0.25,1.5)$, and
$(2,0)$. Each reference job hit its 1,800-second cap. The spacelike cases
retained a vector result but lacked density completion, so no T/L
projection is accepted. Their raw partial outputs remain in the archive.

== Figure 15 comparison and saved artifacts

The scan integrates the isolated NLO coefficient only. Its Figure 15
overlay's native markers use the author's analytic LO plus nine times native
NLO, matching $g=3$. Smooth comparison lines now use the original published
Figure 15 vector paths, calibrated to the printed axes; gray dashed lines
are the published LO curves. The fourteen curves contain fifteen original
path segments and 4,187 vertices. The blue transverse curve retains its
original clipping gap. No smoothing, new reference integration, or changes
to the 64 displayed native values and error bars were applied. Three points
without complete references remain unplotted, as in the prior overlay.

Published paths have graphical precision and retain the original figure's
clipping, including in the full-range overlay. They are not substituted for
the independently computed targets used for residuals or acceptance.

#figure(
  image("convergence-and-agreement.png", width: 100%),
  caption: [Final MC precision and reference agreement. Circles are T,
  open squares L. The suspected $E=3k$ reference and missing targets are
  excluded explicitly; rest T/L overlap. The lower panel's ordinary
  residuals do not demonstrate small errors in the upper panel.],
)

#link("figure15-overlay.pdf")[Full-range Figure 15 overlay] and
#link("figure15-overlay-zoom.pdf")[the fixed-range zoom] show those published
paths with the native markers. The zoom's transverse range begins at −0.4.
The published curves do not show the suspect $E=3k$ reference spike; the
raw computed-reference discrepancy remains in the numerical diagnostics.
The prior sparse-grid overlays are in the archive's checkout copy under
`full-scan-20261019/rendering/recomputed-reference-grid/`. The zoom clips some
native error bars; the full-range plot and CSV preserve those error bars.
`../reference/published-curves.json` contains the calibrated paths;
`figure15-render.json` records their provenance and all displayed markers.
#link("summary.csv")[Physical NLO results],
#link("figure15-comparison.csv")[the LO-plus-NLO comparison],
#link("audit.json")[the interpreted audit], and
#link("manifest.json")[the compact run provenance] record the numbers and
acceptance definitions. This manifest is explicitly curated; the original
controller manifest is archived with its hash recorded here. Original
summary flags are preserved even where this report qualifies their interpretation.

#link("native-results.json")[Consolidated native results] retain all 35 final
raw result objects. #link("reference/point-results.json")[Author point results]
retain all 35 complete or partial source outputs, and the compact manifest
retains runtime settings and seeds. The original per-slot stability records,
checkpoint ladders, executed drivers, and complete generated states are in
the verified archive described by #link("../archive.json")[the archive record],
under the original checkout and `tmp/qcd-fig15-scan-20261018T171925Z/` paths.
Archival and plotting perform no new integrations. No physics changes were made.

The next comparison that can resolve the strongest discrepancy is a
one-sided author-reference check around $E=3k$. Further native convergence
work should focus on the high-variance signed-cut regions and their large
weights, alongside the independently identified interval-root cost documented
in #link("../profiling-20261018/report.typ")[the profiling report].
