# GL638 native sampling acceptance: build10

This milestone validates the faster native sampling source and the corrected
absolute-integral observable. It does **not** yet provide a new GL638 Monte
Carlo estimate. Evidence was collected on 14 September 2026 from build10,
based on `61e37ce2089c6dc73083d9a1fc3dbe8f8c6a1e25` plus its frozen source
snapshot.

## Source precision and acceptance

Mixed GL638 catalogues now prepare every draw at **Fixed256**, rather than
1000-bit arbitrary precision. This remains MPFR arithmetic on every mixed-map
draw; it is not a double-precision fast path. The precision policy is fixed
before sampling, using the requirements of all physical stability lanes. The
completed point, Jacobian/partition factor, and retained cut hosts are embedded
exactly in the canonical Arb representation. Physical rotations and rescue
levels consume that frozen source; they do not regenerate it.

Foreign-channel density evaluation avoids reconstructing coordinates needed
only for a full diagnostic inverse. The selected channel retains its full
inverse/Jacobian check. The LU-h proposal keeps the existing log-logistic
approximation and safeguarded inverse CDF described in
[LU_H_MATCHED_SAMPLING.md](LU_H_MATCHED_SAMPLING.md).

The private saved-state gate independently resolves original cut equations at
the actual promoted source, checks hosted H/Z alignment against the original
surfaces, compares physical cuts and counterterms, and exercises frozen Double,
Quad, and Arb rotation probes. All runs use the complete 936-orientation graph,
six cuts, and 19 threshold variants.

| Gate | Complete Samples | Result | Hosted H/Z cases |
|---|---:|---|---:|
| 1000 GeV, previous seven-channel timing catalogue | 9 | 9 passed | 1 |
| 600 GeV, leading nine channels, m_uv=50 GeV | 11 | 11 passed | 3 |

The 600 GeV deck covers all six standalone cut channels, the right-Cut3
threshold channel, the explicit soft channel, an ordinary joint anchor, and
the two closest retained H/Z approaches. It keeps mu_r=91.188 GeV and uses
componentwise physical stability checks. The three joint real totals agree
with independently solved hosts to relative differences between approximately
1.6e-77 and 4.3e-66. These are agreement checks, not rigorous error bounds.

The earlier 1000 GeV gate failed only its individual near-zero real Original/CT
pieces: differences around 1e-308 were compared purely relatively. Its complete
real estimator already agreed to 1.64e-78. With explicit user approval, only
the individual decomposition comparisons now use each cut's sum of absolute
**same-component** contributions. Real and imaginary budgets remain separate.
Complete estimators, cut totals, geometry, and final Arb probes retain their
strict criteria. The original failed report is preserved; the new report records
every scale, difference, budget, and outcome. Perturbation regressions reject
meaningful cancelling term errors and changes to a tiny complete total.

Ten focused tests passed, covering this criterion, native source/amplitude
external-data ownership, cut maps, absolute-observable covariance and stability,
reporting range, workspace versions, and SymJIT constant slots including O3.
The two private executions independently authenticated the test binary, all
28 frozen source paths, inputs, and all 35 saved-state files before and after.

## Absolute observable and compiler support

Both channel execution modes map the graph-group master once per channel and
reuse that physical point across its members. Physical graph/cut/CT cancellation
occurs before taking the absolute value. For explicit channel sums the observable
is `sum_c |J_c w_c F(k_c)|`, not `|sum_c J_c w_c F(k_c)|`: channel points differ.
One outer draw supplies one statistics update, preserving cross-channel covariance.
Native, Rust, and Python outputs retain this separate absolute payload. Workspace
version 3 requires restarting older checkpoints whose absolute moments used the
previous convention.

GL638 sums all orientations. If orientations are themselves sampled, the
absolute monitor instead describes that sampled-orientation domain; cancellation
between unevaluated orientations cannot be recovered. Separate graph groups
also remain separate integration domains. See the implemented conventions in
[architecture-current.md](../../architecture/architecture-current.md).

Shared process SymJIT activation now passes the requested O0–O3 level
through, including genuine O3 instead of clamping it to O2. Compression remains
separately configurable. This establishes support, not an O3 speed advantage.
The matched build10 evaluator/whole-sample benchmark is recorded in
[GL638_BUILD10_BACKEND_TIMING.md](GL638_BUILD10_BACKEND_TIMING.md). The
previous build9 map-only measurement was **7.569 ms per call**, one worker,
384 calls with native checks passing; it excludes the physical prepass and
evaluators. It must not be substituted for a production timing breakdown.

## First screen and remaining work

The first requested 4096-sample, 600 GeV optimized-LMB screen aborted before
producing an integral estimate. Ordinary output conversion rejected a native
signed component near `-6.2504430634584e-336` because it rounds to binary64 zero.
The remaining initial-grid factor was seven, so the complete component still
rounds to zero. The message does not identify real versus imaginary. This is
a reporting-policy failure, not evidence of a large weight or failed SOCP solve.

A separately reviewed, unapplied proposal permits ordinary rounding to zero
only after checking the remaining factors in native arithmetic. Overflow,
amplified underflow, unrepresentable separately reported factors, and raw
factorized event entries remain guarded. Its nontrivial test-expectation change
awaits approval at this milestone. The same frozen screen seeds should be
retried after its focused checks; no failed sample is replaced or discarded.

There is no final 600 GeV central value/error, new completed variance comparison,
global bounded-weight claim, or convergence claim. The optional m_uv=91.188 GeV
private request is prepared but was not executed in this milestone.

## Evidence locations

These are exact local artifact locations; hashes identify the retained files.

| Artifact | Path | SHA256 |
|---|---|---|
| Optimized build, successful and source unchanged | `/tmp/hosted-joint-optimized-build10-manifest.json` | `19e474d7ece1e53c7f0c180c9379c5ef9e206b964fbba14fbacb2c324cae4a1e` |
| Ten-test execution record | `/tmp/gl638-final-physics-screen/tests-build10-attempt1.json` | `264bccb3068cb38df2fe572a395d8d4f746c5fecce5aae002aefc2fd2a6fb597` |
| New 1000 GeV gate | `/tmp/gl638-hosted-joint-gate/native-state-boundary-build10/report1000.json` | `a53a123b1c36fae3ecf69f6d450d25aa134e969bf38db4f8444474e643453802` |
| New 600 GeV gate | `/tmp/gl638-hosted-joint-gate/native-state-boundary-build10/report600_leading9_muv50.json` | `4c3ed9d349cc71ddf7b8380afc01ebe5f95629d83f4543e4fc046419fac067a9` |
| Preserved original gate failure | `/tmp/gl638-hosted-joint-gate/native-state-boundary-build9/report1000.json` | `988cd64dbb11842a1e26d39fe59da5e6ee432157a69720a49e635e637d38e7e9` |
| Failed first screen | `/tmp/gl638-final-physics-screen/results-screen-component-build10/summary.json` | `4d09201cf0a50d5dccfd55ea53db7833e39b0a979ea88089db1edb63568eab69` |

The build10 gate directory also contains `authentication1000.json`,
`authentication600_leading9_muv50.json`, and `completed_gate_summary.json`.
The map-only measurement is `/tmp/gl638-map-profile/result-build9-fixed256.json`.
The pending conversion proposal, failure triage, and independent audit are in
`/tmp/gl638-final-physics-screen/underflow-proposal/`. Prior compiler timing and
physics context remain in [GL638_COMPRESSED_SYMJIT_TIMING.md](GL638_COMPRESSED_SYMJIT_TIMING.md),
[GL638_HOSTED_JOINT_GATE.md](GL638_HOSTED_JOINT_GATE.md), and
[GL638_SOFT_AND_CUT3_LOCAL_REAL.md](GL638_SOFT_AND_CUT3_LOCAL_REAL.md).

## Build11 follow-up: suppressed real-component instability

The reporting correction described above was subsequently approved and applied.
The identical 600 GeV, m_uv=50 GeV, mu_r=91.188 GeV screen then completed all
4096 draws but retained 16 native instability flags and no NaNs. It used seed
60101, 20 workers, compressed SymJIT O2, all 936 orientations, six cuts and
19 threshold variants. This is diagnostic evidence, not an accepted final
cross-section estimate.

An independent replay reconstructed the original initial Grid/RNG stream,
including complete Sample weights and contiguous worker ranges. The native
slot seed was `16511191977557783548`; all saved checkpoint maxima matched
exact Samples at indices 2290 and 3991. The replay reproduced all 16 failures,
the precision counts (3778 Double, 155 Quad, 163 Arb), and the signed mean/error
exactly. Imaginary and absolute-real aggregates agreed within `4e-16` relative.
Native follow-ups selected the existing rotation objects individually; their
averages reproduced the full Arb calls. Input/state hashes stayed unchanged.
The first follow-up's missing clone warmup was retained as a failed diagnostic;
the corrected continuation evaluated only the authenticated 16 failing Samples.

Every flag came from real or absolute-real components. Imaginary probes agreed
to roughly `1e-294` relative or better, so loosening imaginary precision would
not address the failure. Nine cases had every complete real probe round to
binary64 zero; fifteen had every probe below its minimum normal value.
Draw 3559 exposed why requiring *every* probe to be subnormal was unnecessarily
restrictive: its complete real probes were `+1.26070e-308` and `-3.03194e-308`,
but their mean absolute magnitude was `2.14632e-308`, below
`f64::MIN_POSITIVE = 2^-1022 ≈ 2.22507e-308`.

The implemented componentwise criterion is

```text
mean_j |remaining_weight × native_probe[j].component| < f64::MIN_POSITIVE
```

The triangle inequality bounds the returned rotation average without relying
on cancellation. Every remaining factor is multiplied in native arithmetic
before this test; map Jacobians and channel partitions are already included
at this boundary. Real/imaginary and signed/absolute observables are checked
independently. No imaginary scale or arbitrary absolute floor is borrowed.
All 16 recorded failures satisfy this fixed normal-range bound. Nonfinite
probes/products cannot qualify, and amplifying outer weights can remove the
waiver. Norm checking also bounds its returned primary probe, because that
owner returns the primary rather than the rotation average.

Only the relative-error rejection is waived. Native values, counterterms,
events, and measured relative discrepancies remain intact. Consequently a
stable suppressed component may still report a large native relative error;
that status does not claim accurate relative digits for its tiny value.

Check, clippy, and six focused tests passed. The optimized build12 replay also
passed: all 4096 original draws completed with zero unstable evaluations and
zero NaNs, using the identical card, seed, backend and 20-worker partition.
The new tests are
`stability_checks_bound_complete_underflow_without_cancellation` and
`stability_underflow_waivers_are_componentwise`; they cover the straddling probe,
outer-weight amplification, normal cancellation, mixed components, nonfinite
weights and norm-primary ordering. Existing finite/nonfinite, reporting-range
and independent absolute-observable checks remain applicable.

An independent artifact audit found exact equality of the build11/build12
signed and absolute Re/Im means, errors, zero counts, channel breakdowns and
maxima, including their complete stored coordinates. Precision counts remain
3778 Double, 155 Quad and 163 Arb. The real signed estimate is
`-6.381323993097482e-5 ± 4.2697958507087165e-5`; the absolute-real estimate is
`2.0483568597702613e-4 ± 4.2589482332175714e-5`. These are the small regression
stream's unchanged estimates, not the final production result. No changed
subnormal aggregate is observed; the integration report does not retain every
per-draw value, so this is not a new per-draw equality claim. Only instability
counts, timings and workspace paths differ in the native result. All four
recorded before/after state hash maps agree on the same 35 files, and the
optimized build completed with unchanged sources. Timing differences are not
used to claim a performance improvement.

All paths below are under `/tmp/gl638-final-physics-screen/unstable-replay/`:

| Evidence | Relative artifact path | SHA256 |
|---|---|---|
| Exact 4096-draw replay, including retained failed follow-up | `results-build11/summary.json` | `75d07b3329e3dde923463ccfa322b3384eaca1f856fb0e2be0e1a56e2639463a` |
| Corrected native capture of all 16 failures | `results-native-only-build11/summary.json` | `71fb1c0671708fae9af5efe3c1cca45f516130fe73e64167ea13fc78f15962dc` |
| Complete probe/aggregate classification | `analysis-native-only-build11-normal-range.json` | `7e67195d34c98f49ec6bf289f4b3fc5a81f4c0ef7c4ca367cf52d48360dbbe95` |
| Mean absolute bound for every failing component | `mean-absolute-bound-audit-build11.json` | `bf2e926d047c1ef4f927a9e144a8f585f43ad378bc4ad7a3dcfe0378fa48c2e8` |
| Independent post-fix regression audit and full provenance | `postfix-regression-audit-build12.json` | `62bbb44ebfc91f8b5df353e90711c9b146bbd17c65d306dd03200d68bae50e3d` |

The successful integration report is
`/tmp/gl638-final-physics-screen/results-underflow-replay-o2-build12/summary.json`
(SHA256 `96690391cddb76d34771f78ef8d8692954524824bded2220688c2cc87c8382ed`).
The audit records hashes for its request, execution, native result, driver link,
source build and six-test log; the source-build manifest is
`/tmp/hosted-joint-optimized-build12-manifest.json`
(SHA256 `df8b1f614730003b61628f66ad94c0fd5a66e790cd0388f7c0212b19401277b3`).

`RESULTS_BUILD11.md` contains the 16-row magnitude table. The authenticated
source hashes are `553b202946e5f2a3fa11893d0ce63f57af66bfaf7ab7afa3076b2ff3916b5951`
for `unstable_replay.rs` and
`ecd1b0a156aea8a4c5367364a6664cd9b437d567e89a34c6ae221fef1d5718aa`
for `unstable_native_replay.rs`; their link and execution records are retained
beside the captures. This correction changes no sampling map or threshold
metadata and establishes no new variance or bounded-weight claim.

## E_cm-relative agreement for a vanishing component

A subsequent 50-worker MC calibration retained fifteen unresolved evaluations.
The first iteration's eight failures were reproduced as exact Samples: all
probes were finite, with complete real magnitudes between approximately
`2.75e-308` and `7.39e-306`, above the machine-underflow rule's bound. Their
imaginary components agreed accurately. This is a distinct case: increasing a
physical tolerance must not be disguised as changing floating-point range.
The other seven flags and the later c=4 SUM validation's single flag retain
their failed status until independently replayed or reproduced.

The optional runtime agreement test is now dimensionless and defaults to zero:

```toml
[stability]
integrated_energy_dimension = -2

# Add these fields to the final Arb entry of the existing stability stack.
[[stability.levels]]
precision = "Arb"
required_precision_for_re = 1e-12
required_precision_for_im = 1e-12
ecm_relative_tolerance_for_re = 1e-100
ecm_relative_tolerance_for_im = 0.0
escalate_for_large_weight_threshold = -1.0
```

For every component, compare the native fully weighted probe disagreement with
the quantity the checker actually returns, divided by `U * E_cm^d`. Here `d`
is the declared integrated energy dimension after flux normalization and `U`
is the conversion actually applied to reported units. The component checker
returns its rotation average; norm checking returns the primary probe. Direct
momentum densities subtract their missing `3L` powers from `d`. Physical
evaluations require an explicit dimension when this allowance is enabled;
arbitrary numerator dimensions are not inferred from UV power counting.
Gaussian reference values and raw second moments instead have known dimensions
zero and two and use no cross-section unit conversion.

The optional test waives only relative-disagreement rejection. It preserves
native values, measured relative discrepancies, large-weight escalation,
nonfinite rejection and the separately configured exact-zero policy. Signed
and absolute estimators, and real and imaginary components, remain independent.
Reference second moments use the stricter of the two component allowances.
No dimensional cutoff or unscaled absolute-tolerance field is introduced.

Compilation, Clippy and twelve focused tests passed, including the existing
production cut fixture extended to check actual cube/raw source handling,
reference admission and physical-dimension errors. Energy/unit rescaling,
outer-weight amplification, returned-value bounds and independent absolute
payloads are covered. The Python settings test is added but awaits a rebuilt
extension. The test run is
`/tmp/gl638-final-physics-screen/tests-build13-attempt2.log`; earlier compile
fixture mistakes and a nested serialization-guard deadlock are retained in the
preceding logs and were corrected without relaxing assertions.

The optimized build13 replay now passes all eleven exact saved Samples: the
eight captured failures and three stable controls. With the allowance disabled,
the original eight unstable statuses are reproduced; enabling only the final
Arb real allowance accepts them. All four returned components (signed and
absolute, real and imaginary) remain exactly unchanged, including the controls.
Independent native probes satisfy the configured component bounds. The replay
took 150.0 seconds and left the state and inputs unchanged; it is a saved-point
regression, not a new integration or validation of the other unresolved flags.
The complete capture is
`/tmp/gl638-final-physics-screen/unstable-replay/results-ecm-regression-build13/summary.json`
(SHA256 `ebf30236b45bc7300549036b17c567981b288a3e76518e12b46452b971c2efcc`);
`analysis-ecm-regression-build13.json` in the parent directory contains the
sixteen aggregate checks and all eleven component-by-component audits.
