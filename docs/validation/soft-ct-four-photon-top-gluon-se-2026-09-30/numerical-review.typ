= Four-photon numerical validation contract

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.

Read-only review of the prepared native evaluation commands at source commit
`032969abd45c53e0a5bb10d883d7b1d3a4bb149b`. This note describes what the
measurements can establish; it does not claim that pending measurements
have completed. Paths below are relative to `crates/gammalooprs/src/`
unless another root is given.

== Point comparisons

Compare signed complex values in a fixed graph, subtraction prescription,
basis, external configuration, and physical momentum point. The prepared
analysis compares Arb-evaluated results between routes at relative
tolerance `1e-9`; inspect exports those results as f64. It separately
records Double-versus-Arb relative differences. A completed finite point
phase is not a claim that those precision differences pass a tolerance.
All three selected points being nonzero is an analysis requirement, not
a physical theorem.

`inspect --momentum-space --graph-id 0` without an orientation selector
uses the summed evaluation for these states. The forced-precision flag
is `--use_arb_prec`; its CLI declaration explicitly documents f64 output
(`crates/gammaloop-api/src/commands/inspect.rs:52–53`). Pass the nine
momentum coordinates as separate argument tokens. The first ghost point
attempt failed to parse the initial single comma-separated token; its
failed receipt was retained and the corrected `points-spaced` phase rerun.
The runner preserves the frozen binary and saved state and supplies
`--read-only-state --no-save-state`.

== Selected central soft limit

`profile bulk --select 'GL262 S(e5)'` selects one independent soft momentum.
The corresponding GL256 command changes only the graph selector. Since
`Q7 = -Q5`, selecting both edges as separate soft directions would be
incorrect. Absence of `--per-orientation` gives one orientation-summed
report (`integrands/process/ir.rs:678–711`). It does not validate each
orientation separately; local cancellation is tested in their full sum.

The profiler selects a physical basis containing edge 5 and holds its
other loop coordinates fixed. It may choose companions different from
the generation basis, then routes samples into the canonical basis
(`ir.rs:1000–1020,1037–1089`). Therefore compare the saved ray fingerprints
between routes; matching seeds alone do not certify identical physical
rays. A gluon/ghost sum additionally needs the signed bubble-coordinate
mapping in `routing-notes.typ`.

For the saved states' `e_cm = 173`, the prepared range
`lambda = 1e-4 ... 1e-8` gives `|Q5| = 0.0173 ... 0.00000173`;
companion loop vectors are random unit vectors,
not rescaled by `e_cm` (`ir.rs:1017–1020,2065–2067`). These are spatial soft
samples with `lambda > 0`. Bulk profiling forces arbitrary precision,
independently of the card's Double ladder (`ir.rs:2081–2097`).

The fitted observable is the magnitude of the complete complex sum,
`|F(lambda)|` (`ir.rs:2143–2146`). The fitting routine takes adjacent
differences on a geometric grid, removing a constant before log-log
regression (`utils/fitting.rs:82–91,160–188`). The reported soft score is
that fitted exponent plus three per soft momentum (`ir.rs:700–711`).
For a singular `F ~ lambda^(-1)`, the score is therefore `+2`. If the
integrand tends to a nonzero constant, however, the fit describes the
first varying correction rather than the leading constant; a purely
constant magnitude has no such fit.

Power counting motivates a generic FullH score near `+2` or stronger
after removing the first two self-energy soft coefficients. The completed
native `S(e5)` profiles measured GL262 score `2.0003899636` with
R-squared `0.99999968275`, and GL256 score `2.0000075054` with R-squared
`0.999999999883`. These establish numerical behavior on the native
selected rays, not a graph-pair cancellation or an algebraic identity.
Additional gauge, ghost, or orientation cancellations can change the
leading power on other observables.
The native verdict tests only `score > 0`; it has no R-squared gate.
The prepared analyzer reproduces that verdict and requires finite
fit fields, so report the actual score and R-squared separately. A
poor-quality positive fit cannot establish the expected asymptotic power.
The saved native bulk JSON contains the fit and ray fingerprint, not the
raw signed samples, and cannot by itself demonstrate a combined
gluon-plus-ghost cancellation.

== Explicit matched gluon-plus-ghost ray

A separate completed diagnostic evaluated 17 shared values of
`lambda = 1e-4 ... 1e-8` at canonical coordinates
`q = lambda*(30,30,114.5)`, `h = (47,53,145.5)`,
`b_gluon = (37,41,43)`, and `b_ghost = -b_gluon`.
Both graphs used forced-Arb inspect, their complete orientation sums,
and their original included graph factors. No extra sign or symmetry
factor was applied when adding the exported complex values.

The direct log-log fit of `lambda^3 * |F_gluon + F_ghost|` gives slope
`1.9969923056` and R-squared `0.9999976745` over all 17 points. Its final
nine points (`1e-6 ... 1e-8`) give slope `1.9999182746` and R-squared
`0.999999999415`; the corresponding individual slopes are
`1.9999110414` for GL262 and `2.0000660983` for GL256. These fits do not
drop a constant. The result is consistent with raw combined behavior
`|F| ~ lambda^(-1)` on this explicit ray.

The minimum ratio `|F_gluon + F_ghost| / (|F_gluon| + |F_ghost|)` is
`0.99999208`, so the two graphs do not produce an appreciable extra
cancellation here. Although the Arb outputs are rounded to f64 before
the graph sum, this ratio shows that this particular addition does not
suffer strong cancellation. It remains a numerical check, not an exact
pole-cancellation certificate.

Both bounded phases completed with unchanged binary and saved state:
GL262 in 18.688 s, GL256 in 6.924 s. The campaign files
`paired-soft-plan.json` and `paired-soft-analysis.json` retain the exact
commands, per-point raw signed complex values, fit data, hashes, and
receipts. This ray is distinct from the native profiler's selected
`[5,10,13]` ray.

== Ultraviolet coverage and verdict

`profile ultra-violet --selected-limits all` means every nonempty loop
subset in every generated physical basis (`uv/profile.rs:69–132`).
Generation stores all spanning-tree bases
(`processes/amplitude.rs:1468–1475`, `graph/lmb.rs:1443–1455`); the
amplitude profiler uses that full saved list over the physical graph
with no Cutkosky restriction (`uv/profile.rs:2345–2355`).

Exact enumeration of the ten internal edges on eight vertices gives
40 spanning trees for each graph. GL262 and GL256 have the same undirected
topology; reversing ghost edge 9 does not alter the count. Three loops
give seven nonempty subsets per basis, hence 280 rays and 4,760 summed
evaluations at 17 points per graph. A completed all-limit report should
contain 40 distinct basis labels, seven distinct nonempty free subsets
for each, 280 rows in total, and 17 valid samples per row. These are
expected coverage counts derived from the source and graph topology;
they are not a claim that all rays have been evaluated. The initial
analysis helper requires nonempty rows but does not itself enforce the
280-row count.

At each UV ray the native profiler multiplies the integrand magnitude
by `scale^(3*n_scaled_loops)` before fitting
(`uv/profile.rs:2851–2874,3140–3148`). The displayed UV slope is therefore
a measure-weighted degree, not the raw integrand power. The prepared
range is `scale = 1e8 ... 1e12`. The CLI flag `--use_f128` currently flows
into the arbitrary-precision evaluation switch
(`uv/profile.rs:2821–2868,3094–3110`); the historical flag spelling is
not evidence that this path used binary128.

The native pass criterion is finite slope at most `-0.9` and finite
R-squared at least `0.99`; the analysis helper reproduces it
(`uv/profile.rs:59–62,2219–2237`). Amplitudes have no exemption for a
missing fit. The reported slope may use a high-scale suffix rather than
the entire range: the fitter chooses the longest suffix reaching
R-squared `0.99`, retaining at least nine of these 17 valid points
(`uv/profile.rs:3048–3082`). Report the full-range slope and R-squared
alongside the accepted fit when they differ materially. The analysis
requires all 17 samples to be present, which is stricter than accepting
a native fit after dropping nonpositive or nonfinite samples.

The minimum evaluation budget implied by all-limit coverage is substantial:
4,760 times the measured per-point cost, plus loading and routing. Use the
completed point timings to judge feasibility under the existing cap;
an interrupted sweep is not a completed all-limit validation.
