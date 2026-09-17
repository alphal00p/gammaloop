= Physical expression-generation performance
<physical-expression-generation-performance>
The first large regression occurs at C3. GL1 rises from 38.904 s / 2.137 GB at C2 to 166.063 s / 17.700 GB, and the attachment first exceeds 30 GB with the common sparse contraction order. Fresh repeats confirm both observations. The earlier order does not recover GL1 performance; most cost remains with integrated UV disabled. Complete upstream captures and fixed-executor controls now establish material costs from lost factor sharing and larger tensor sums. The separate C6 wait is localized to relative epsilon expansion of an exact zero. The #link("raised-energy-cff-structural-causes.typ")[structural-cause report] records these causal checks; no production remedy or whole-generation speedup has been validated. The smaller sunrise cases remain inexpensive.

All 46 primary attempts are recorded: 37 passes, six memory caps, two program failures and one monitoring failure. All 13 repetitions/one-setting controls are terminal: 12 passes and one memory cap. Separate diagnostic results and limits are reported below. The reviewed branch is published at `01982f688b6810d5152f0317713db412ae7ea794`. The comparison fixes current remote main at `bc8bfbad717b80cb83ca1bf7ac2aec704917d424`, the original stack base at `395610143576507503fd2c785db3ba62340f4277`, and every reviewed C1--C8. No performance fix has been applied.

=== Workloads and comparability
<workloads-and-comparability>
The supplied inputs define four workloads: the two-loop projected g→g GL1 graph with integrated UV; vacuum QCD sunrise with integrated UV; sunrise with local UV and a computed-forest export; and the attached three-loop photon graph with a central quark bubble and selected orientation58. The two sunrise cards also have different sampling settings, so their difference is not solely an integrated-UV toggle.

The attached `intermediate_cost` contraction strategy is unavailable on main and C1--C2. Those exact combinations are unsupported. The same physical workload is measured there with an explicitly labeled `sparse_atom_aware` control, also available on all later commits. The default for the other cards changes at C3; a separate GL1 control holds the older order fixed. Each revision retains its own native locked dependencies and model files; no dependency override is applied.

All runs use fresh state directories, one generation core as requested, eight Rayon workers, 128 MiB Rust worker stacks, and a 30 decimal GB process-tree RSS limit. Binaries are freshly checked/built and copied with hashes before each source change. Timings exclude builds. Sampled process-tree RSS and the application's own RAM report are distinct measurements.

=== Source pins and build environment
<source-pins-and-build-environment>
#figure(
  align(center)[#table(
    columns: (33.33%, 33.33%, 33.33%),
    align: (auto,auto,auto,),
    table.header([Revision], [Commit], [Change relevant to this investigation],),
    table.hline(),
    [main-current], [`bc8bfbad717b80cb83ca1bf7ac2aec704917d424`], [Remote main at the recorded fetch],
    [main-base], [`395610143576507503fd2c785db3ba62340f4277`], [Original stack base],
    [C1], [`665658b168981c5c508e6cd51d3b57bbb55a27fb`], [Symbolic/tensor foundations and dependency migration],
    [C2], [`a9397ee0f86b115fdfab50e81c48dde2137692e7`], [Shared CFF engine, not yet used by these CLI paths],
    [C3], [`3ea313789a1a5290093579febd7d1c601b92594e`], [Exact branch maps and local UV integration; first large regression],
    [C4], [`04637884f24fb28a247574f04e03308cf4c04afb`], [Command/evaluation workflows; large-workload costs persist],
    [C5], [`1f2cf6d8236de91138951af8380bf1fd1eb4e783`], [Physical phases/model sewing; export normalization changes],
    [C6], [`8e4e4664f3f9303760eb18d7972b91ff867b12ea`], [D-dimensional integrated algebra; native GL1 observation incomplete],
    [C7], [`6d7a9220129c0fac37bdc0e04ebefccb1688c7be`], [Projected local UV changes; immediate GL1 parser failure and larger attachment scalar],
    [C8], [`01982f688b6810d5152f0317713db412ae7ea794`], [Documentation/evidence only relative to C7],
  )]
  , kind: table
  )

Both main pins use their locked Symbolica dev revision `0441bd7a511209dce2ca99925fe87f8b18e4bf03`; C1--C8 use the locked alphal00p revision `4d0a833eb8e059d1f95bdae5abed2559830b235f`. No moving dependency branch was updated for this comparison. C1 therefore bundles source and dependency effects. C1 and C2 produce byte-identical CLI binaries (`b88334faa510c6c9a7783e5dbb4238eeb5d81e59815c7a9a0aca34cfef748e9b`).

The host uses AMD EPYC 9754 processors, Linux 6.18.45, Rust/Cargo 1.97.0 and jj 0.44.0. Runs share this host with unrelated processes; the measurement controller excludes competing investigation builds/workloads, not all host activity. Repeats characterize observed variability and are not confidence intervals. Builds use `dev-optim`, four Cargo workers and each revision's native locked dependencies. Native `cargo check` and build commands and the source, binary, card and model hashes are retained in individual receipts. The separate scratch probes also run formatting checks and Clippy.

=== Complete primary matrix
<complete-primary-matrix>
These are individual native observations. Full-card wall time includes imports, generation, explicit exports and saving. The original attachment setting is unavailable on four older revisions; each has an explicitly labeled common-strategy observation instead.

#figure(
  align(center)[#table(
    columns: (16.67%, 16.67%, 16.67%, 16.67%, 16.67%, 16.67%),
    align: (auto,auto,auto,auto,auto,auto,),
    table.header([Revision], [g GL1 integrated], [Sunrise integrated], [Sunrise local export], [Orientation58 exact intermediate], [Orientation58 sparse control],),
    table.hline(),
    [main-current], [PASS 73.639s / 3.180GB], [PASS 5.501s / 0.053GB], [PASS 2.538s / 0.047GB], [EXACT\_UNSUPPORTED], [PASS 109.066s / 8.443GB],
    [main-base], [PASS 75.909s / 3.459GB], [PASS 5.698s / 0.055GB], [PASS 2.663s / 0.049GB], [EXACT\_UNSUPPORTED], [PASS 114.278s / 8.383GB],
    [c1], [PASS 39.082s / 1.860GB], [PASS 4.687s / 0.050GB], [PASS 2.151s / 0.046GB], [EXACT\_UNSUPPORTED], [PASS 37.497s / 3.264GB],
    [c2], [PASS 38.904s / 2.137GB], [PASS 4.708s / 0.048GB], [PASS 2.135s / 0.043GB], [EXACT\_UNSUPPORTED], [PASS 36.660s / 3.102GB],
    [c3], [PASS 166.063s / 17.700GB], [PASS 5.941s / 0.052GB], [PASS 3.050s / 0.054GB], [PASS 59.884s / 7.015GB], [MEMORY\_CAP 49.706s / 30.618GB],
    [c4], [PASS 164.877s / 17.767GB], [PASS 5.683s / 0.054GB], [PASS 2.640s / 0.045GB], [PASS 59.493s / 8.354GB], [MEMORY\_CAP 46.912s / 30.146GB],
    [c5], [PASS 161.936s / 17.696GB], [PASS 5.550s / 0.058GB], [PASS 2.583s / 0.048GB], [PASS 60.117s / 6.885GB], [MEMORY\_CAP 48.355s / 30.242GB],
    [c6], [MONITOR\_FAILURE 875.207s / 0.330GB], [PASS 4.307s / 0.053GB], [PASS 2.602s / 0.048GB], [PASS 59.523s / 7.033GB], [MEMORY\_CAP 48.703s / 30.484GB],
    [c7], [CHILD\_FAILURE 0.931s / 0.041GB], [PASS 4.075s / 0.052GB], [PASS 2.299s / 0.049GB], [PASS 67.629s / 7.263GB], [MEMORY\_CAP 43.015s / 30.724GB],
    [c8], [CHILD\_FAILURE 0.859s / 0.042GB], [PASS 4.216s / 0.054GB], [PASS 2.620s / 0.049GB], [PASS 62.144s / 7.734GB], [MEMORY\_CAP 39.235s / 30.680GB],
  )]
  , kind: table
  )

Failed and capped observations are not completed generation timings. The cap is sampled, so the recorded peak can exceed its threshold. The first three cards retain each revision's native default contraction strategy; this changes at C3. Their initial differences cannot isolate that setting from other changes.

On main-current, local sunrise generation itself took 0.925 s: 0.411 s Spenso preparation/execution and 0.032 s Symbolica evaluator construction. The 0.482 s non-evaluator remainder also includes preprocessing and other generation setup; it is not a separate expression-construction timer. Its exported bare node reproduces the user's quoted main expression exactly. The GL1 generation took 72.490 s, of which Spenso took 15.685 s and Symbolica evaluator construction 31.528 s.

The attached case already has substantial cost on main. Its trace reports 69.294 s in tensor preparation/execution and 33.114 s in Symbolica evaluator construction. Within the first interval, contraction execution took 27.450 s and the remaining post-execution interval requires further attribution. No stack regression follows from these baseline numbers alone.

C8 with the requested intermediate strategy produces a 394,421,795-byte scalar Atom. Tensor preparation/execution takes 14.581 s and Symbolica evaluator construction 44.114 s. Changing only the contraction strategy to the common sparse control exceeds the cap during tensor execution. This establishes a material strategy-dependent memory effect on C8; the prefix sweep locates its onset at C3.

C8 GL1 rejects a color sum because `trace(cof(3), sym(t(a), t(b), t(c)))` is inferred as scalar while the companion structure-constant term exposes adjoint indices. Source inspection locates the newly exposed analytic UV parsing boundary at C6; C5 completes this card. C6 runs for more than 14 minutes before a monitoring failure; its native attempt has no completed timing. C7 reproduces C8's immediate parse error; C6 remains inconclusive as a completed workload. Relevant model rules agree between main and C8. The failure must not be reported as faster generation.

=== First two commits
<first-two-commits>
Both C1 and C2 complete all four supported physical workloads. The common sparse attachment control takes 37.497 s / 3.264 GB on C1 and 36.660 s / 3.102 GB on C2. GL1 takes 39.082 s and 38.904 s respectively. These are first-pass observations; they exclude C1 or C2 alone as the onset of the later 30 GB failure.

C1 and C2 CLI binaries are byte-identical. C2 adds the shared engine, which GammaLoop starts using at C3. Main/C1/C2 saved GL1 graphs are byte-identical, complete saved models and parameter cards agree, and all eight sunrise node DOTs and numerator values agree exactly. C1's attachment contraction falls from main's 27.450 s to 2.560 s; the following scalar retrieval and alias-restoration interval falls from 41.229 s to 3.481 s. Matching alias summaries support comparable preparation, but complete main/C1 attachment input equality has not been certified. C1 bundles source and Symbolica dependency changes, so these measurements do not isolate one optimization.

=== Measured onset at C3
<measured-onset-at-c3>
C3 is the first substantial regression in the larger workloads. GL1 grows from 38.904 s / 2.137 GB at C2 to 166.063 s / 17.700 GB at C3. C4 and C5 reproduce that scale: 164.877 s / 17.767 GB and 161.936 s / 17.696 GB. The common sparse attachment control first exceeds the 30 GB cap at C3 and also caps at C4/C5. With the requested intermediate strategy, the attachment completes in 59.884 s, 59.493 s and 60.117 s at C3/C4/C5 respectively.

GL1's recorded Spenso stage grows from 6.728 s at C2 to 71.859 s at C3, and Symbolica evaluator construction from 15.954 s to 64.271 s. The Spenso interval includes color/gamma simplification and network preparation, so it must not be attributed entirely to contraction. Fresh GL1 repetitions preserve the onset: main-current takes 75.744 s, C1 40.988 s, and C3 168.300 s. Completed one-setting controls are detailed below. The small sunrise cases show changed expression sizes and modest absolute generation costs, without a large process-memory regression.

C7 changes the attached scalar result from approximately 334 MB to 394 MB. Its tensor stage takes 14.730 s and Symbolica construction 49.465 s, compared with C6's 18.724 s and 35.152 s. C8 takes 14.581 s and 44.114 s for those stages. The complete C7→C8 diff contains documentation/evidence only, so that timing difference cannot be credited to a production optimization. Native repetitions are needed to characterize variability; these are individual observations. An isolated C7 build removing only the two new direct-route early color calls restores the complete C6 pre-network input byte-for-byte after stripping outer whitespace. This identifies the cause of the input regrouping; it is not a measurement of the resulting scalar size or full generation time.

=== Repetitions and one-setting controls
<repetitions-and-one-setting-controls>
Fresh unchanged GL1 runs take 73.639--75.744 s on main-current, 39.082--40.988 s on C1, and 166.063--168.300 s on C3 (two observations each). The attachment takes 109.066--115.130 s on main-current with sparse order, 37.328--37.497 s on C1 with sparse order, and 59.884--60.709 s on C3 with intermediate order. Both C3 sparse attempts exceed 30 GB. These ranges confirm the material onset without assigning statistical confidence from two samples.

One-setting diagnostics on C3 distinguish the source of cost:

#figure(
  align(center)[#table(
    columns: (21.43%, 28.57%, 28.57%, 21.43%),
    align: (auto,right,right,auto,),
    table.header([Change from native card], [Full-card wall], [Sampled RSS], [Interpretation],),
    table.hline(),
    [GL1: sparse contraction order], [201.664 s], [17.920 GB], [The earlier order does not recover C1 performance],
    [GL1: integrated UV disabled], [138.212 s], [17.630 GB], [Large cost remains without integrated counterterms],
    [GL1: UV subtraction disabled], [1.315 s], [0.062 GB], [Large cost depends on UV contributions],
    [Attachment: UV subtraction disabled], [0.887 s], [0.039 GB], [The selected bare graph alone does not cause the large cost],
    [Attachment: algebra enabled], [113.020 s], [14.001 GB], [Enabling this preprocessing does not remove the regression],
    [C8 GL1: integrated UV disabled], [144.602 s], [18.929 GB], [The integrated-path error disappears, while the local cost remains],
  )]
  , kind: table
  )

The UV-toggle rows intentionally change physics content to locate the cost; they are not replacements for the native workload benchmarks. In particular, `generate_integrated=false` still constructs local four-dimensional Taylor expressions. Disabling all UV subtraction excludes the forest contributions, so it does not identify one local-UV operation by itself. The attachment algebra control changes preprocessing while retaining physics content. C8 local-only GL1 completes, separating its integrated parser failure from the large local contribution cost. Its Spenso stage takes 67.448 s and Symbolica construction 69.307 s. The corresponding C3 local-only stages take 68.891 s and 61.095 s, respectively.

The attachment without UV returns a 32,634-byte scalar Atom, versus 334,489,168 bytes in the native C3 workload. Enabling algebra returns a smaller 256,421,594-byte scalar, yet its Spenso stage takes 82.490 s, including 74.909 s across the six sequential execution intervals. The preprocessing changes the incoming expression; this is not an identical-input executor comparison. Final scalar size alone does not characterize intermediate cost.

=== Upstream expression structure
<upstream-expression-structure>
The user requested that actual expression structure before tensor execution take priority over contraction heuristics. Four separate diagnostic runs captured the complete C1/C3 amplitude before evaluator algebra, after gamma simplification when enabled, and immediately before tensor parsing. All four completed; their logging overhead is excluded from the primary performance table.

GL1's raw expression grows from 874,634 to 2,162,342 UTF-8 bytes. After gamma simplification, packed Atom storage grows from 1,960,490 to 4,160,869 bytes; immediately before tensor parsing it is 1,960,105 versus 4,160,869 bytes. Top-level terms decrease from 31 to 18, so term counts alone cannot describe the relevant structure.

The attachment's pre-network input grows from 1,027,567 to 1,766,072 packed bytes. C3 contains only the requested exact selector `sigma(58)`, ruling out additional selected production orientations as a sufficient explanation for this case. Its complete bare factor audit retains the same coefficient, denominator, six vertex-gamma factors and product of six fermion numerator factors. Compact slashed-vector sums become explicit indexed gamma factors multiplied by two-term vector factors. The product remains factorized; the location of tensor sums changes. The separate proper-UV witness below identifies additional recoverable sharing; its isolated replay shows little timing benefit from extracting those gamma factors alone.

Raw-amplitude logging hides namespaces, while the post-gamma and pre-network payloads preserve full plain syntax. No diagnostic distributed a graph numerator. Serialized lengths, packed Atom bytes and process RSS remain separate quantities.

==== Complete factor-sharing witness
<complete-factor-sharing-witness>
A captured C3 proper-UV subexpression has three summands with exactly the same eight fully qualified gamma factors. Its complete enclosing contribution is `-2 * (G*R0 + G*R1 + G*R2) * metric * sigma(58)`. Extracting `G*(R0+R1+R2)` preserves every coefficient, sign, power, denominator and index. The 559,179-byte original sum, all complete remainders, common product and a factor-preserving alternative are saved. Independent text checks confirm the full payload and a reversible replacement of exactly this subexpression. The grouping matches the earlier raw-amplitude hook and reaches tensor parsing unchanged. This is an exact within-C3 identity, not a matched C1/C3 UV-sector certificate or a measured fix for the memory cap. Those original aggregate captures lacked per-forest provenance. Later owner captures join every C1/C3 pre-network contribution byte-exactly to its original payload; term order or gamma subsets alone still do not certify cross-version scalar equality.

The physical GL1 bare contribution shows a different form of lost sharing: C1 keeps its intact six-term three-gluon factor and quark numerator outside a sum of 18 scalar selector weights. C3 repeats sign-specialized complete tensor products in 18 branches. All five-edge sign-support rows and all bare tensor factors are accounted for. The three-gluon sum itself and the quark vector sums remain factorized. Full scalar CFF-weight equivalence and independent runtime catalog decoding are separate, outstanding checks.

==== Fixed-executor factor replay
<fixed-executor-factor-replay>
Both complete enclosing contributions parse and contract successfully using one immutable C3 executable, `Sequential / MinIntermediateCost`, the native completion pass and scalar-alias resolution. Formatting, Cargo check, optimized build and Clippy passed; both full source audits passed. The only input change is the certified extraction of the eight common gamma factors.

#figure(
  align(center)[#table(
    columns: 3,
    align: (auto,right,right,),
    table.header([Observation], [Original], [Gamma factors extracted],),
    table.hline(),
    [Input UTF-8 bytes], [559,305], [557,291],
    [Network parse], [0.082 s], [0.076 s],
    [Tensor execution], [3.140 s], [3.088 s],
    [Network parse/prepare/contract/resolve interval], [3.320 s], [3.237 s],
    [Packed scalar Atom bytes], [136,523,268], [139,627,496],
    [Canonical output writing], [12.523 s], [13.053 s],
    [Full guarded wall time], [16.944 s], [17.247 s],
    [Sampled process-tree RSS], [4.984 GB], [4.659 GB],
  )]
  , kind: table
  )

This one pair shows no material contraction-time improvement or total-time win. The RSS includes output serialization, and a single pair does not establish a robust memory improvement. Removing this small repeated gamma product alone is therefore insufficient evidence for a performance fix. The later metric and recursive factor controls address different boundaries and show material improvements, as detailed in the structural-cause report.

The exact input identity is proven by factor-preserving reassembly. Complete canonical scalar outputs have different representations and sizes. A separate checker passes formatting, Cargo check, optimized build and Clippy, then runs in 97.037 s / 1.014 GB. Direct Atom equality is false. All three predetermined finite, nonzero scalar evaluations fail the unchanged `2^-128` threshold: the scaled residuals are `1.317e-35`, `3.994e-30` and `5.226e-30`, about `1e-16` relative to the actual output magnitudes. All 18 scalar leaves and six distinct negative-power bases are checked; `sigma(58)` and `tree_denoms(1)` are fixed to one. This is a formal scalar comparison, not physical phase space.

Requesting 256 bits does not provide that precision throughout this check. The saved canonical outputs contain decimal coefficients without precision annotations. The native parser reads the observed `8.000…` and `16.000…` forms at 53 bits; evaluation preserves that lower coefficient precision, and subsequent multiplication/division propagate it. This limits the reparsed output check; it does not measure the original in-memory Atoms' precision. The receipt therefore records successful checker execution with semantic `NUMERICAL_DISAGREEMENT`. The result does not establish a contraction defect or certify exact returned-value equality. No tolerance, input or oracle was changed, and no graph numerator was expanded. This table remains an isolated contraction replay of one proper-UV contribution, not a full generation benchmark or proof of a matched C1/C3 contribution.

=== Established source boundaries
<established-source-boundaries>
C3 first replaces a shared parametric numerator times theta-selected CFF with a sum of exact-residue selectors times independently mapped complete numerators. C1/C2 retain the earlier assembly. Complete owner/selector captures now certify this duplication at proper-UV Taylor entry. C3 also changes the contraction default and processes existing tensor summands separately.

Factor-preserving controls replay complete C3 inputs with one fixed executor. Moving two shared metrics outside the attachment's hot sum lowers its rank 10 → 6 and its parse/prepare/contract/restore interval 11.300 → 0.738 seconds. Recursive extraction of intact common factors lowers the two hot GL1 intervals 10.656 → 4.320 and 49.460 → 20.148 seconds. These establish a material cost of expression structure, not a full-generation remedy; one GL1 scalar output grows despite faster contraction.

Two earlier representation boundaries are concrete. C1 collects common factors after multiplying the cograph numerator (`crates/gammalooprs/src/uv/approx/final_integrand.rs`, lines 130--143 at C1). C3 omits that collection and preserves separate denominator sums (the same file, lines 283--297 at C3). C1 reaches rank-one compaction before the remaining aggregate momentum split, whereas C3 maps to `Q3 + E*delta` first. The unchanged Idenso rule cannot compact an additive vector as an intact rank-one function. This ordering is not a general explanation for the hot contributions: attachment owner 14, forest `{11vs} · {168C} · {16C4}`, has no compact gamma arguments in either input, and GL1 hot terms 015/016 are gamma-free. The matched owner 14 instead has 12 shared outer gamma factors and maximum tensor-Add rank three at C1, versus 108 nested gamma occurrences and rank six at C3. Complete scalar equality between those versions remains unproved.

Restoring `collect_factors()` globally would have a tradeoff. The actual locked Symbolica implementation extracts common factors recursively but also extracts negative powers missing from other summands. For example, `G*A/D1 + G*B/D2` becomes `G*(A*D2+B*D1)/(D1*D2)`. This is algebraically valid, but combines denominator topologies that C3 deliberately keeps separate. The saved gamma-only witness preserves each denominator inside its original remainder. No production change is proposed as already validated.

The default HedgePoset pipeline still performs 4D Taylor work when integrated UV is disabled. Computed UV exports recompute forests. C5 separately removes a duplicate spatial-measure normalization from amplitude exports, explaining the quoted pi^12 versus pi^6 difference. These boundaries must be distinguished from expression growth.

=== Minimal reproduction of the GL1 parser failure
<minimal-reproduction-of-the-gl1-parser-failure>
The independent fixed-C8 public-API probe reaches all four cases. Formatting, Cargo check, optimized build and Clippy pass. The executable exits 101 because the unchanged contract assertions identify a bug; this is recorded as `CHILD_FAILURE` / `OBSERVED_CONTRACT_FAILURE`, not a passing test run.

All four mathematical inputs should expose three adjoint indices before closure. The ordered trace and its two complete decompositions have the independent SU(3) closed oracle `Tr(Ta Tb Tc) * f_abc = 6i` in the native normalization. The symmetric trace alone contracts to zero. The independent color algebra control passes before any parser cases run.

#figure(
  align(center)[#table(
    columns: (25%, 25%, 25%, 25%),
    align: (auto,auto,auto,auto,),
    table.header([Public parser input], [Exposed indices], [Complete closed value], [Result],),
    table.hline(),
    [Symmetric three-generator trace], [Incorrectly empty], [Correctly zero], [Contract mismatch],
    [Ordered three-generator trace], [Correct three indices], [`6i`], [Pass],
    [Small projector-expanded color sum], [Correct three indices], [`6i`], [Pass],
    [Native symmetric-plus-antisymmetric sum], [Rejected as scalar plus tensor], [Not reached], [Parse error],
  )]
  , kind: table
  )

The last error matches the physical GL1 failure. The scalar zero in the first row alone would miss the public tensor-structure bug. The native trace helper unwraps `cyclic` but does not expose the generators inside `sym`; subsequent structure inference treats that trace as scalar.

C6 introduces the analytic UV parser call. C7 changes the input before that unchanged call: early color simplification introduces the symmetric-trace form, and local Taylor processing defers Dirac reduction to the integrated path. Idenso and Spenso are unchanged across C6--C8. The immediate failure first occurs in the measured C7/C8 native runs; the C6 native workload remains incomplete. This is a pre-existing parsing limitation exposed by changed upstream input, not evidence of a Vakint integration-kernel bug. The standalone color-projector expansion is small and does not distribute a graph numerator.

=== C6 integrated expression boundaries
<c6-integrated-expression-boundaries>
A separate C6 diagnostic reaches its 1200-second limit without completing; the guarded observation lasts 1200.742 s and peaks at 629,071,872 bytes RSS. Its final log contains 43 complete expression payloads spanning only the first 4.921 seconds. One local result holds the same three color generators and one structure constant outside two intact sums. A later integrated expression repeats that exact color product in all 248 summands: 744 generator occurrences and 248 structure constants. The full expressions are retained; this is an observation of lost sharing, not an equality certificate between those two records. The log omits forest/owner identities, so adjacency cannot certify the local input corresponding to that integrated call.

The local expression contains no gamma, chain, bispinor or compact trace function. A residual color trace therefore cannot trigger the conditional analytic spin expansion for that input. The remaining source boundary before the integrated log collects Lorentz tensors, simplifies metrics, collects bispinors and collects gamma chains. These calls can reorganize tensor-bearing sums, but the current capture does not isolate which call creates the repeated color factors.

The 2,357,013-byte integrated expression is followed by ten further simplifier entries within those first five seconds. Their scalar Vakint coefficients multiply an external Lorentz/color basis; they contain no remaining gamma chains. Their syntax is consistent with the projected-numerator caller of the same simplifier. That filter did not record their completion or backend phases. The optional stack attachment was denied by the host; no stack frames were captured.

Later phase logging completes all ten projected-numerator simplifications and localizes the wait to term 1 in Vakint's `epsilon_pole_order`. Its complete 565-packed-byte coefficient is an unchanged tensor product times an exact scalar zero. The locked Symbolica relative-series loop keeps doubling depth while searching for a nonzero leading term. Denominator-only algebra proves the cancellation; no graph numerator is distributed.

The captured ten-term batch times out at 60 seconds with the original relative query. Changing only that pole-order query to absolute order zero completes all ten terms in 0.834 seconds, including 47 coefficient pole checks and ten backend evaluations. Settings, inputs and native dependencies are preserved; final epsilon expansions are unchanged. This diagnostic localizes the wait but does not replace the incomplete native full-card timing or validate all integrated results. The structural-cause report retains the complete proof.

=== Sunrise representation comparison
<sunrise-representation-comparison>
The actual saved main and C8 bare roots match branch by branch after compact-vector materialization and the known `(2*pi)^6` export-measure conversion. Complete nonselector factors agree, including coefficients, denominators, gamma incidence and vector sums. The graph DOTs are byte-identical and graph-used model rules agree. This comparison preserves factorization; it does not distribute numerator products.

The two selected bodies analytically contract to the same scalar `8*(E0*E1 + dot_E(vecQ0,vecQ1))`. That gives a possible sharing opportunity after contraction. The complete C8 runtime catalog independently confirms `sigma(0)` selects `(-,+,+)` and `sigma(1)` selects `(+,-,-)` on edges 0, 1 and 2. All eight sign rows match the bare-expression proof; the other six select neither branch. Production and export code preserve the ID/row pairing. Bare-root correspondence is not a full-forest or numerical validation.

=== Reproduction and interpretation
<reproduction-and-interpretation>
The recorded native build commands are:

```sh
cargo check --locked -p gammaloop-api --bin gammaloop --profile dev-optim
cargo build --locked -p gammaloop-api --bin gammaloop --profile dev-optim
```

The controller runs each preserved executable as `gammaloop -n -s FRESH_STATE EFFECTIVE_CARD run generate` through the existing licensed launcher and RAM watchdog. Each actual argv, working directory, source audit, native lock/model hashes, effective card and immutable executable hash is retained in its receipt. The card inventory binds every adapted card back to the supplied input. Existing attempt directories are immutable; reproduction needs fresh directories and the matching audited source checkout.

`measure.py` is the shared command/measurement implementation. `controls/post-sweep-plan.json` records all seven repeats and six one-setting deltas. `build-scratch-probes.py`, `run-factor-pair.py`, `run-color-parser-probe.py`, `run-c6-boundary-diagnostic-v2.py` and `diagnostics/common-factor-output-comparison/run-checker.py` record the separate diagnostic protocols. The earlier C6 launcher is retained as an explicitly invalid-filter setup attempt. None updates snapshots, automatically retries a failed workload or expands a graph numerator for comparison.

The matrix measures generation and requested export/save operations. Native traces do not independently time every export/recomputation operation; full-card wall time minus generation also includes imports, saving and other setup. It is not an integration-accuracy test, a full physical numerical equivalence certificate, or a benchmark of every contraction heuristic. The smaller sunrise workloads complete throughout with low absolute cost. The attachment remains faster than main with the newer intermediate strategy, but regresses relative to C1/C2 and fails under the common sparse strategy. Those are distinct comparisons.

=== Coverage and remaining limits
<coverage-and-remaining-limits>
The complete primary matrix and 13 controls establish the C3 onset. The corrected C6 diagnostic ends at its separate 1200-second limit; the original native observation remains classified as a monitoring failure. Its process scan exceeded the watchdog's 10-second deadline at 856.121 s and cleanup finished at 875.207 s. Neither observation supplies a completed C6 GL1 time. The first diagnostic's invalid filter, explicit stop and child exit 143 remain recorded as a setup failure, separate from both observations.

The complete GL1 DOT and all 15 graph-used saved model records match main; complete parameter cards agree. Unrelated whole-model records differ. The sunrise runtime catalog and all eight selector truth-table rows pass. These certificates do not establish full GL1 scalar CFF-weight equality. Later captures establish matching proper-UV owners and exact within-C3 input identities, without assuming complete scalar equality across revisions.

The original gamma-only replay has an exact input-reassembly certificate; its returned outputs still fail the strict numerical check, with reparsed coefficient precision as a documented limitation. The separate metric and GL1 factor pairs pass the unchanged three-point threshold after independently certified exact coefficient parsing. Neither is an exact full-output identity or a physical phase-space certificate.

The later controls establish causal costs from lost factor sharing and larger tensor sums, localize the C6 zero-series wait, and reverse C7's input regrouping. They do not establish a production remedy or a whole-generation speedup. No production source, existing test or published stack commit has changed; the diagnosed color parser and relative-series bugs remain unfixed.

=== Evidence
<evidence>
The original matrix evidence directory is `/tmp/raised-energy-expression-perf-20260914`. The later causal evidence is under `/tmp/raised-energy-causal-20260914` and is indexed by the #link("raised-energy-cff-structural-causes.typ")[structural-cause report];. The archive below retains the original evidence selection, including the revision/card inventory, measurement and control drivers, terminal receipts, full expression captures, factor-preserving proofs, public probe sources, graph/model certificates and source-boundary audits. The #link("raised-energy-cff-expression-performance-evidence.tar.gz")[evidence archive] and #link("raised-energy-cff-expression-performance-evidence.json")[per-file manifest] retain explicitly selected text evidence, including both complete scalar outputs. The archive writer verifies terminal receipts, input hashes, UTF-8 content and every archived member, and scans the selected content for listed credential formats. The licensed launcher, process environments and compiled binaries are excluded. Reproduction requires a separately configured license and rebuilding the recorded source revisions.
