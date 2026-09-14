# GL638 matched-count Monte Carlo screen

All seven candidate catalogues complete the predeclared screen: 21 runs,
43,008 physical draws, three independent seeds and 6,144 samples per method.
All reported samples are stable and finite. The frozen confirmation methods
are **optimized LMB**, **six cuts**, and **six cuts + composed Cut1 LU-h→H/Z**.
The last choice has a score only 0.03947% below direct joint; this is a selection
result, not evidence that the composed chart has lower true variance.

The [seven-candidate preflight](GL638_HOSTED_JOINT_GATE.md#completed-cut-matched-acceptance)
and [local H/Z comparison](LOCAL_CUT_MATCHED_HZ.md) precede this screen.
Rust source is unchanged from `4d192b5be`, using optimized build5 and the frozen
`b53df2b4` driver source. Every method evaluates all 936 orientations and all
six physical cuts with the same numerator, UV terms, threshold counterterms,
WH/WF factors and precision settings. All 35 saved-file hashes are unchanged.
The seven cards and exact direct/composed/soft definitions remain linked from
the local comparison; this study changes neither the physical calculation nor
its acceptance tolerances.

## Fixed-count results

Each run is one fresh, untrained iteration of 2,048 samples, with seeds 11113,
22229 and 33331 and 20 workers. The complete sampling probabilities enter the
reported estimator once. Reusing a seed across catalogues does not mean that
they evaluate identical physical points. The table pools the original workspace
statistics using the [predeclared formula](gl638_hosted_joint_gate/cut_matched_acceptance/mc/MC_ANALYSIS_PROTOCOL.md),
including within-run and between-run contributions; it does not average three
reported errors.

Means and empirical standard errors are in **10⁻⁴ pb**. `Abs Re` and `Abs Im`
are separate absolute-component integrals, not the integral of the complex
norm. Maxima are the largest absolute **signed-component** extrema among the
included samples, in pb; they are not a bound or the largest sampled complex
norm. Score is the worst of the four empirical sample variances divided by the
corresponding common baseline absolute-component mean squared.

| Method | Re ± SE | Abs Re ± SE | Im ± SE | Abs Im ± SE | Max Re / Im (pb) | Score |
|---|---:|---:|---:|---:|---:|---:|
| optimized_lmb | 1.4317 ± 1.927 | 6.1834 ± 1.925 | 1.0504 ± 2.465 | 8.9363 ± 2.463 | 1.04892 / 1.34675 | 596.6986 |
| cuts | 1.8994 ± 1.471 | 2.1214 ± 1.471 | 1.915 ± 1.699 | 2.3827 ± 1.699 | 0.875679 / 1.0342 | 347.7284 |
| cuts_joint | 0.19781 ± 0.5581 | 1.1609 ± 0.5579 | 0.35659 ± 0.2857 | 0.8743 ± 0.2855 | 0.224439 / 0.112123 | 50.04617 |
| cuts_soft | -4.5453 ± 6.468 | 10.479 ± 6.467 | 20.417 ± 19.47 | 23.437 ± 19.47 | 3.76756 / 11.9369 | 29160.24 |
| cuts_joint_soft | -5.4561 ± 7.389 | 11.7 ± 7.388 | 23.353 ± 22.25 | 26.571 ± 22.25 | 4.30578 / 13.6421 | 38085.04 |
| cuts_combined_joint | 0.18392 ± 0.558 | 1.1552 ± 0.5578 | 0.35645 ± 0.2853 | 0.88482 ± 0.2851 | 0.224439 / 0.112123 | 50.02642 |
| cuts_combined_joint_soft | -5.4704 ± 7.389 | 11.711 ± 7.388 | 23.36 ± 22.25 | 26.582 ± 22.25 | 4.30578 / 13.6421 | 38085.04 |

The common absolute scales are 0.0006183385979 pb (real) and
0.0008936315639 pb (imaginary). Both joint variants without the soft channel
have smaller observed errors and component maxima than their cut-only control
in this screen. The soft variants retain large outliers and errors; none are
excluded after observing them. The local composed-chart reduction of about
7.9 on the two H/Z rays does not translate into a similarly separated screen
score: other sampled regions still contribute. Physical attribution of the
retained maxima is a separate follow-up.

These are empirical errors from a small screen of a potentially heavy-tailed
integrand. Stability acceptance and finite observed maxima do not establish
finite variance, a global weight bound, or a final cross-section estimate.
The full report retains individual seed results, signed and absolute estimates,
component extrema and matched comparisons, including unfavorable outcomes.

## Worker choice and frozen confirmation

The declared worker control precedes advanced-method statistics. It uses
2,048 baseline samples and the existing Integrate interval, including warmup;
peak RSS is the process measurement. Thirty workers are chosen only if their
throughput is at least 1.35 times the 20-worker throughput and RSS is at most
300 GB. Both completed controls are retained.

| Workers | Integrate time (s) | Samples/s | Peak RSS (GB, decimal) | Included in screen |
| --- | ---: | ---: | ---: | --- |
| 20 | 82.14215 | 24.93239 | 70.07079 | Yes, baseline seed 11113 once |
| 30 | 110.54657 | 18.52613 | 100.48298 | No |

The measured throughput ratio is 0.74305, so the screen and confirmation use
20 workers. This is the prescribed resource choice for this run, not a general
scaling claim. Performance optimization has stopped.

The fixed score rule chooses `c = cuts_combined_joint`; its matched non-joint
control has the same absent soft policy, `b = cuts`, and
`a = optimized_lmb`. The [frozen confirmation definition](gl638_hosted_joint_gate/mc_screen/confirmation/manifest.json)
uses these exact cards, source, preflight and common scales with independent
seeds 10007, 20011, 30013, 40009 and 50021: 32,768 samples per run,
163,840 per method, 491,520 in total. It keeps `c` as the primary estimate even
if confirmation does not show an improvement. Confirmation results are not
included in the screen archive; no reselection or additional favorable runs
are authorized by this screen.

The first confirmation attempt subsequently stopped on `optimized_lmb` seed
20011: canonical physical overlap preparation for cut group 3 could not certify
its Arb1000 root before body evaluation. The three seed-10007 runs completed;
they are not pooled as confirmation. The saved-state hashes remain unchanged.
The [pending failure note](/tmp/gl638-hosted-joint-gate/CONFIRMATION_FAILURE_PENDING.md)
and [attempt record](/tmp/gl638-hosted-joint-gate/confirmation-build5-run.json)
retain this interruption separately while exact-source/root-owner investigation
proceeds. No cause or repair is established by the screen results.

## Retained evidence

The [independent audit](gl638_hosted_joint_gate/mc_screen/independent/summary.json)
passes 167 checks, including Decimal70 recomputation of all 28 pooled
observables, exact mode/seed/count/worker and selection checks, all 93 recorded
analysis input hashes, and the unchanged saved-state hashes. The
[archive index](gl638_hosted_joint_gate/mc_screen/artifact_hashes.json) retains
all four original frozen directory trees, the excluded worker control, exact
run/log/resource records, frozen selection, independent audit source/results,
and confirmation declaration. Twenty-two original integration checkpoints
preserve complete Samples: 21 included runs and the excluded control.
No maximum-replay or confirmation outcomes are included.

Identical settings, reports and state-hash files share stored bytes; each
original path and its uncompressed hash remain explicit. These are the actual
separate workspaces, not a constructed aggregate checkpoint. The
[93-input reference map](gl638_hosted_joint_gate/mc_screen/analysis_input_references.json)
and [provenance](gl638_hosted_joint_gate/mc_screen/provenance.json) link the
unchanged cards, source/build records, analyzer and pre-results protocol.
Unmodified audit source retains its original absolute paths; restore the
indexed files in an isolated layout to repeat the artifact-only analysis.
No generated-process saved-state payloads, binaries or new physics runs are
part of this archive.
