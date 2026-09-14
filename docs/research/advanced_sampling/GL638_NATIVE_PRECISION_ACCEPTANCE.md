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
