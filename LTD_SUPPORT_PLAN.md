# LTD support for GammaLoop

## 1. Goal and implementation startup

**Goal:** Add native LTD generation and evaluation alongside CFF, preserving factorized numerators and complete LTD residue sums, with representation selection at each stability level. Establish local equivalence against current main’s CFF implementations, including individual physical cut weights and threshold-counterterm contributions.

Upon starting implementation:

1. Create `ltd_support` from the current clean main checkout.
2. Write this complete plan, including the corrected request below, verbatim into `./LTD_SUPPORT_PLAN.md`.
3. Register the goal statement above with `create_goal`, referring to that file as the implementation and acceptance specification.
4. Set the top-level `enable = false` in `nix-ci.nix`.

Current main does not implement LTD. Adapt the core implementation from `origin/ltd_in_gammaloop_symbolica_update_monolith`, inspected at commit `1b73243b4b3201eacde9c187b09302a8482bb631`. Build its GammaLoop integration against current main.

The explicitly requested Markdown plan overrides the repository’s documentation-format default. Add a narrowly scoped exception for this exact root file to documentation validation; keep other documentation in Typst and do not create a parallel plan.

## 2. Public behavior and interfaces

- Add `global.generation.3dreps`, accepting `"cff"`, `"ltd"`, or a nonempty list. Default to `["cff"]`. Preserve list order, remove duplicate entries while retaining first occurrence, and serialize canonically as a list.
- Reuse the existing `RepresentationMode` enum across generation, settings, persistence, and runtime selection. Add the necessary schema, serialization, and Python support rather than introducing a competing representation enum.
- Add optional `3drep = "cff" | "ltd"` to each stability level. **An omitted selection uses the first generated representation**, as agreed.
- Reject an explicitly requested representation that was not generated. Validate before expensive generation and before runtime evaluation.
- Any generation request containing LTD requires `local_uv_cts_from_expanded_4d_integrands = true`; never silently change the UV prescription.
- Requesting LTD automatically implies complete residue summation. Reject generation/runtime orientation filters, orientation Monte Carlo, and explicit individual-orientation inspection. Apply this consistently to mixed CFF/LTD builds.
- Keep ordinary CFF-only defaults and orientation-local operation unchanged.
- Remove the rejected `runtime.general.use_ltd` placeholder.
- Support ordered stacks such as Double/LTD → Double/CFF → Quad/CFF. Repeated precision values remain separate attempts.
- Forced arbitrary precision retains the configured Arb level’s representation. If an Arb level must be synthesized, inherit the first configured level’s resolved representation.
- Expose the representation in attempt diagnostics and recorded stability results. Existing precision-wide summaries may remain aggregates.
- Update settings completion, generated schemas, saved defaults, examples, and API documentation, respecting `SHOWDEFAULTS`.

## 3. Implementation changes

### Native LTD generation

Reuse current `energy_residues`, exact routing solvers, `LinearEnergyExpr`, repeated-denominator detection, surface storage, and expression assembly.

For ordinary LTD:

- Solve complete affine energy maps at each residue basis.
- Retain basis-local on-shell-energy factors and the contour prefactor.
- Represent each uncut quadratic denominator as its two signed linear factors, preserving explicit E/H surfaces.
- Do not perform CFF numerator-degree analysis or capacity optimization.

For raised propagators:

- Adapt the reference’s confluent residue construction and exact rational finite-difference sampling.
- Reuse numerator-degree analysis **only where raised poles require it**, following the corrected requirement.
- Identify degeneracies using complete signed routing, external shifts, masses, powers, and existing denominator provenance.
- For residue order \(\alpha\), combine numerator sampling of order \(\beta\leq\alpha\) with scalar-denominator derivatives of order \(\alpha-\beta\), including exact factorial coefficients.
- Obtain numerator derivative contributions through exact polynomial sampling weights, without differentiating or expanding the graph numerator.
- Preserve one consistent affine loop-energy map and its signed edge bindings for physical sampling.
- Do not import the reference’s restricted no-bound affine fallback as a general algorithm. Do not import its contact-completion machinery or compensating signs without establishing their necessity against current conventions.

Adapt core algorithms selectively; do not transplant the old GammaLoop integration or restore its compatibility paths.

### Factorized numerators and UV composition

- Extend the existing numerator-family machinery. It already represents exact coefficients of internal on-shell energies, external energies, \(M\), and constants.
- Retain one factorized numerator definition per family and small scalar argument rows for its residue maps. Preserve complete map identity; equal orientation labels do not establish equal maps.
- Keep numerator definitions shared through tensor preprocessing, UV composition, threshold operations, and evaluator construction.
- Carry representation selection through physical graph generation, local-4D UV components, and the outer cograph.
- Preserve current exact numerator/denominator reconstruction certificates. LTD bypasses CFF-specific capacity assignment and proposal scoring, not reconstruction correctness.
- Share representation-independent graph preparation, local-4D Taylor data, and integrated Vakint results. Dispatch into CFF/LTD projection below that shared boundary, for both supported forest backends.
- Include representation in caches whose contents depend on it.

### Cuts, thresholds, and normalization

- Share canonical physical cut and threshold identities across representations.
- Discover candidates from topology and match them against LTD’s explicit denominator factors. Do not require fully directed CFF orientations for LTD.
- Preserve exact factor signs, proportionality coefficients, support, and multiplicities. Normalize selected poles into current positive-energy cut/threshold coordinates once; retain spectator signs.
- After threshold localization, transform affine numerator arguments and affected denominator factors consistently. Preserve physical cut identity even when a transformed factor becomes an H-surface.
- Handle signs and multiplicities per denominator chain, including numerator cancellations and raised poles.
- Extend the existing residue-selection abstraction with representation-aware behavior while preserving CFF semantics.
- Retain current physical LU normalization. Avoid importing graph-dependent parity corrections, origin-string contact deletion, or whole-numerator expansion from the reference.
- Use current family-aware symbolic operations for required threshold derivatives and series.

### Evaluators, runtime, and persistence

- Store representation-specific generated expressions, numerator families, and evaluator programs within the existing graph/evaluator owners. Share physical topology, sampling geometry, cut identities, and event grouping.
- Separate “complete physical residue sum required” from the evaluator’s internal execution layout.
- Preserve source-local coefficient rows for `SingleParametric`; sum every required row internally. Independent UV-child and outer-source sums must remain independent, so counterterms are neither omitted nor multiplied by the production residue count.
- Support existing evaluator methods and compiled backends without duplicating numerator definitions per residue.
- Pass resolved representation through evaluation contexts rather than mutating settings during retries.
- Reuse the same prepared sample, graph/channel selection, Jacobian, and overlap geometry across representation changes. Prevent representation-dependent cache reuse and retain events only from the accepted attempt.
- Persist the ordered generated list and representation-indexed payloads. Update state loading, standalone export, generated standalone evaluation, and compiled artifact naming. Existing saved integrands may require regeneration; no legacy migration layer is required.

## 4. Validation and completion

**Independent algebraic checks**

- Use current main’s CFF generation as the oracle. Some reference comparisons reuse the LTD raised-pole builder through their CFF path and are therefore insufficient.
- Cover simple poles; powers two through four; equal and unequal masses; reversed routing; repeated spectators; multiple raised axes; initial-state cuts; trees; and disconnected components.
- Test polynomial numerators through relevant convergence boundaries, mixed energy dependence, and sampling-scale independence.
- Extend existing repeated-channel routing tests and add simple analytic residue checks.
- Verify signed cut normalization, raised cuts, left/right and iterated thresholds, and localization that changes E/H classification.

**Factorization and runtime checks**

- Verify that numerator tensor preprocessing remains shared as residue count increases.
- Verify exact contribution counting for independent UV-source sums.
- Test scalar/list settings, empty/unknown rejection, ordered defaults, unavailable representations, and inferred complete summation.
- Test Double/LTD → Double/CFF fallback and later precision escalation, including cached evaluation, forced Arb, and accepted-event isolation.
- Exercise eager/compiled evaluation, relevant evaluator methods, state round trips, standalone export, and Python/settings round trips.

**Existing integration tests**

- Extend the existing scalar cross-section matrix with LTD/from4D, retaining orientation-local CFF/from3D, explicit-sum CFF/from3D, and CFF/from4D.
- Preserve all existing scalar graph/numerator cases, subtraction settings, tolerances, and UV-profile checks.
- Compare each graph and sample by physical cut and sampling channel after complete residue summation. Preserve orientation-local CFF assertions separately.
- Extend `assert_evaluation_outputs_match` to compare channel identity and detailed threshold-counterterm decomposition: component identities, occurrences, original/bare/weighted values, multipliers, and skip flags. Require nonempty decompositions in targeted fixtures.
- Retain Double totals, Quad component comparisons, and precise Arb checks. Use adequate precision for LTD cancellation; do not weaken CFF expectations.
- Extend existing raised self-energy, numerator-cancellation, normal/reversed/dotted threshold, `gg → hhh`, LMB-normalization, and iterated-threshold inspect tests.
- Run the complete scalar suite with `just test_LU_scalar_xs`, including slow cases.

For discrepancies, follow the repository’s pipeline-first triage: identify the first differing boundary, then reduce through graph, forest, cut, residue order, and term while preserving numerator factorization.

Before final review, review the diff for duplication, run formatting and Clippy, update dependency/CI metadata where required, enable CI, and complete `just ci-checks-and-upload` on `itphlies` before pushing CI-enabled work. Preserve logs and commit identity; apply the repository’s non-draft/`final-review` sequence when preparing the PR. Claim completion notifications only if actually registered.

## 5. Original request, with b.1 corrected as subsequently instructed

The following preserves the original wording except for the corrected numerator-degree requirement in b.1. The decisions above incorporate the subsequent clarifications.

> Ok the goal here is to add LTD support (the particular 3-dimensional representation) to gammaloop.
> The existing crate `three_dimensional_rep` already implements the LTD mode I believe (there was an upstream abandonned PR that already had LTD support but we want to redo it from scratch here using the updated main).
> Because the LTD expression is identical (in terms of a numerical evaluation; faster to evaluate but numerically less stable) to the CFF one, from the user point of view, there is only one change:
>
> a) An additional gammaloop global generation option called `3dreps` allowing to specify which 3-dimensional representations to build the integrand for. By default, this list would contain only `cff` but the user could set it to `ltd` or `[cff,ltd]` (empty not allowed). Moreover, whenever `ltd` is requested, the UV local CTs generation *MUST* be specified as being obtained from 4D (the other option is ineligible for LTD). This is because one key motivation for generating the local UV CTs from 3D is that the UV subtraction becomes orientation local (i.e. each individual key of the residue map becomes separately UV finite and integrable).
> However, in the LTD expression, each residue in the residue map can never be integrated separately anyway because they the H-surface dual cancellation mechanism only works when summing over all residues, meaning that there is no point in keep any residue localization for LTD.
> Moreover, now in the stability stack, the user can add one specification for each level which is what 3d representation to use (i.e. `cff` or `ltd`, of course provided that those have been selected at generation time for this integrand). The reason being that one may often want to try the faster `ltd` expression in double-precision (as it is faster, but less stable), and only fallback to `cff` (still in double precision) if unstable, as it's more stable but slower, and only if still unstable escalate to higher precision in either of the modes.
>
> Then in terms of implementation, everything is analoguous to the current implementation except for the couple of caveats to keep in mind:
>
> b.1) Unlike the `cff` expression, the `ltd` expression without raised propagators does not care about the maximal power of each `emr` energy in the numerator: its residue map remains valid independently of this information, so long as the whole energy integral in the LMB energy variables is UV finite, which it is in the cases considered here. For raised propagators, however, numerator EMR energy powers do matter for the derivative-free LTD residue construction, and the existing numerator-degree analysis must be reused to provide the required bounds.
>
> b.2) However, unlike the `cff` expression, the `ltd` expression requires special care in the presence of raised propagators, as it would normally require numerator derivatives to handle the higher-order residue. Thankfully, the residue map of the generalised_ltd 3drep crate offers a way around this and a special handling of such cases that avoids the need for derivatives. However, it does require then to either pass the information about such raised degenerate propagators (or make sure the generalised_3drep crate can reconstruct this information fine.
>
> b.3) The `ltd` expression has keys of the residue maps that are bit more complicated, i.e. written not only as multiples of the sampling scale M, but also as a linear combination of all external momenta energies and LMB OSEs. It is very important however that you keep the numerator factorizing the sum of all the residue map values, so that (for generation efficiency point of view). so you must right the energy argument of the numerator function in terms of parametric coefficients of that linear decomposition, and then only concretize them at the end when summing over all entries of the residue map (or not even, when using the `SingleParametric` evalautor approach). Basically follow how it's already done for the CFF expression, except that the generic decomposition supporting all possible keys of the LTD residue map are a bit more complicated.
>
> b.4) Cut / threshold identification is more difficult with LTD (especially getting the correct sign) but it *does* remain possible, because the values of the residue map have the quadratic denominators written as products of E-surfaces and H-surfaces, so that the requires E-surfaces can be found explicitly.
>
> b.5) Slow scalar cross-section test suite and local equivalence with CFF from either local UV CTs from 3D or and from 4D (both should already match in current main and they are authoritative). Make sure to compare individual cut weights for each graph and sample, as well as threshold CT weights (from the additional weights of the event decomposition) and they must be equal. Bake this into the existing slow scalar cross-section tests but also add such comparison to other exising local inspect tests of diverse integrands.
>
> Thoroughly review the existing code and implementation in view of the plan above include in the plan a goal statement, this prompt verbatim, and that upon starting its implementation you must write it verbatim into `./LTD_SUPPORT_PLAN.md` and add it to yourself as a goal statement.
>
> Work your implementation within a new feature branch called `ltd_support`.
>
> Ask for anything unclear to you.

## 6. Accepted clarification: physical higher-pole sums and large one-loop benchmarks

The subsequent user clarification supersedes the per-derivative-order equality requirement above without changing the original plan text:

- LTD generation, physical cut identification, and threshold subtraction must use only the LTD expression. Do not construct CFF expressions or enumerate CFF residue keys to produce any LTD payload.
- For simple poles, compare the individual physical Cutkosky cut/residue and threshold contributions locally, retaining exact signed normalization.
- For a degenerate higher-order pole, sum all partial weights belonging to that one physical cut/residue and the same observable function. Compare that sum locally between CFF and LTD; its individual derivative-order pieces need not agree. Keep different physical cuts, threshold occurrences, sampling channels, and observable components separate.
- Retain detailed comparisons within CFF and preserve existing numerical tolerances. Cross-representation tests must reject compensating errors between distinct physical contributions.
- Reproduce cases (a) and (b), the decagon and triacontagon, from https://arxiv.org/pdf/1906.06138 using LTD exclusively. Add the decagon to the test suite and the triacontagon as a scalar CLI example with its TOML command card and DOT graph input, following the repository's existing scalar-example directory layout.
- Verify that these one-loop LTD constructions contain respectively 10 and 30 physical residues and retain factorized numerators, without a hidden CFF construction that would defeat the intended scaling.

This clarification resolves the formal off-shell threshold-component blocker: complete physical higher-pole contributions are the acceptance boundary, rather than CFF's separate formal derivative-order components.

## 7. Accepted clarification: dimensioned comparisons for certified zeros

The user approved absolute comparisons for analytically certified-zero totals and physical cut/threshold contributions, including CFF comparisons, with this condition:

> It's ok to add absolute zero comparison, but the threshold you compare to must be a dimensionless parameter multiplying some "characteristic" scale that carries the dimension, like E com or something like this, this way the dimensionless threshold can be assessed as being expectedly small given arithmetic precision considered. Of course keeping in ming that it's expected for LTD to be numerically less stable than CFF.

Every such comparison requires an exact source/cut certificate and a dimensionless precision-dependent tolerance multiplied by a characteristic scale with the tested quantity's dimensions. Preserve nonzero relative tolerances, physical identities, multipliers, skip flags, and the adaptive/Double and precise Arb checks. Transport raw additional weights through their recorded normalization before applying a physical scale; detailed event decomposition values already include event normalization.

## 8. Accepted clarification: disclose upstream dependency issues

The user requested testing an upgrade to SymJIT 2.26.0 and instructed:

> Continue as planned, but do not hide any bug or issues found in Symbolica or symjit behind a workaround you implement. I have to escalate those upstream for proper fixes, but for now it's good that you have a local patch to be able to continue.

Keep confirmed upstream failures, exact dependency versions, independent reproducers, and local patch scope explicit so they can be escalated upstream. Distinguish upstream defects or unsupported capabilities from GammaLoop integration errors. A local patch permits continued implementation; it does not establish that upstream has fixed the issue.

The user subsequently directed:

> Yes these performance issues are fixed in some other upstream branch, I would not focus on these now if you can still push through with current implementation on the setup as-is now.

Continue correctness validation on the current setup. Retain only the local adaptations needed to complete it; do not broaden the task into unrelated performance investigations or another dependency migration.

## 9. Additional request: configurable overlap-center objectives and heuristics

The user requested parallel investigation of additional SOCP center objectives:

> Continue as planned, but is there an objective function set for the solving of the center in the SOCP problem? Because this would make that center unique then.
> We should actually expose option in the runtime parameters for which convex objective function to pick. Something like the relaxed chebychev center objective function would be a good option and another would be to minimize the sum of the absolute values (or square if absolute values are inlegible) of the E-surfaces involved. Can you run subagents to investigate the possibility of adding these two modes (and the tests for them, all the while keeping the current mode default).

The subsequent corrections determine the requested signed-sum objective and heuristic behavior:

> Sorry I meant MAXIMIZE the sum of... before

> Actually they must all be negative within the domain, so indeed minimize the signed summed value of all the e-surfaces the center must be within,. Send that correction tot he agent.

> (but of course only when the heuristic centers are not valid, and it should be possible to choose to enable or not the use of testing heuristic centers)

Preserve the existing objective and heuristic behavior by default. Investigate a conservative relaxed Chebyshev objective and minimization of the signed E-surface sum, with runtime selection and a switch controlling heuristic-center tests. Enabled valid heuristics take precedence over optimization; explicitly forced centers retain their validation. An objective alone does not guarantee uniqueness. The signed-sum objective also requires an explicit interior-margin policy because its constrained minimum can lie on a boundary despite a nonempty strict overlap; do not silently change objectives or misclassify that case as no overlap.

The proposed implementation uses `max_min_depth` (the default), `relaxed_chebyshev`, and `min_sum`, plus `enable_heuristics = true`. The final strict-interior check remains authoritative, including when a solver returns a valid point without certifying objective optimality. The relaxed Chebyshev mode uses a conservative Lipschitz bound in the active solve LMB's Euclidean metric; it does not promise an exact nonlinear insphere or a basis-independent radius.

The initial proposed signed-sum safeguard retained 1% of a certified common energy depth. The user subsequently emphasized that the objective must favor geometric separation from all thresholds for numerical stability and explicitly selected **“Preserve best clearance (recommended)”**. The accepted signed-sum mode therefore preserves the best relaxed-Chebyshev clearance within solver accuracy while minimizing the signed sum as a secondary objective. Remove the proposed 1% tradeoff and its configuration field. The distinction between energy margin and geometric clearance must remain explicit in the documentation and tests.

The user further required:

> (keep in mind that the result of the SOCP problem with these different choices of objective functions must remain deterministic, and also it should never fail when it would not fail in the absence of the objective function, and also not drastically slow down center finding. Of course during the maximal overlap structure building, not such objective is needed, only for the center finding part)

Keep overlap discovery independent of the new optional objectives. Refine only the final overlap-group centers, reusing certified feasibility witnesses. Preserve deterministic constraint ordering and bound refinement work. An optional optimization failure must not discard a valid overlap or fail an otherwise valid sample: retain the certified witness and expose the refinement outcome in diagnostics. This explicitly required nonfailure behavior takes precedence over rejecting a valid witness solely because an optional optimization did not certify optimality.

After the geometry tests and initial CLI measurements, the user narrowed further work:

> it's ok, don't overdo optimization of the new center finding modes for now, and focus on finalizing LTD support.

Retain the implemented bounded center modes and their validation. Prioritize the complete LTD scalar suite and final checks over further center optimization or benchmarking.

## 10. Additional acceptance test: fast integrated center comparison

The user subsequently required:

> (just make sure there are enough tests for the new center finding modes, including some *FAST* numerical integration ones, e.g. a scalar one-loop box with non-trivial overlap structure and centers, where we demonstrate that all three center setting strategies agree at the integrated level within MC error, (aiming for about a 3 minutes total testing time). Put a subagent on implementing this test, and make it explicitly test that the centers used and found differ in the three integration setups compared)

Add a fast physical scalar integration regression covering all three objectives,
with a target of about three minutes total execution time. Compare the centers
actually used, require meaningful differences beyond solver roundoff, and check
integrated agreement within Monte Carlo uncertainties at useful statistical
precision. Keep the physical fixture and integration controls common across
the three runs. A small two-loop fixture is eligible when the ordinary one-loop
box's proportional depth and relaxed-Chebyshev objectives cannot provide a
genuine three-way center comparison. Preserve the selected clearance policy.

The user suggested the published Box4E fixture:

> For the test I mentioned, consider the kinematics of the `Box4E` example given in sect. 3.1 of https://arxiv.org/pdf/1912.09291.
> It has a non-trivial E-surface structure with 4 centers needed and there the various center definition will not match.

Retain a Box4E geometry regression with its four maximal overlap groups. For the
accepted Lipschitz relaxation all its surfaces have the same bound, so the
default and relaxed-Chebyshev objectives are proportional. Document that
limitation and use a mixed-mass two-loop fixture for the genuine three-way
integrated comparison; do not change objectives to force different centers.

## 11. Additional acceptance: native LTD CLI diagnostics

The user required:

> (continue as planned, but of course make sure the diagnostic CLI command `3drep` also works well with the new LTD mode, and displays the residue map elegantly, with some tests, etc.. etc.. Use inspiration from how it was handled in the upstream abandonned branch you took LTD from. We still don't want a separate "evaluation" track in this `3drep` command though, so you must not port that in here).

Use the existing symbolic validation/build command and display owner. Identify
the selected representation correctly, render full affine residue-map keys,
retain signed surface factors and raised-pole sampling maps, and test CLI
selection, exported artifacts, and readable summary/detail output. Do not add
a separate numerical evaluation subcommand or pathway.

## 12. Additional acceptance: generated-content display and representation timings

The user required:

> Continue as planned, but also review the "display" command, and make sure it clearly identifies what was generated, and also for the generation statistics it would be good to show separately the broken down timings for CFF and LTD when both were requested.

Extend the existing display and generation-report owners. Display the actual
generated representations, including after state loading, rather than inferring
them from mutable generation settings. Distinguish imported, ungenerated data
from generated integrands. Report separately measured CFF and LTD generation
stages and evaluator counts while retaining aggregate totals. Account for shared
graph preparation, four-dimensional Taylor data, and integrated UV work once;
do not charge that shared work to whichever representation is generated first.
Cover single-representation and mixed builds, report aggregation and persistence,
and the rendered display with focused tests.

## 13. Follow-up: explicit setting names and precision/performance studies

The user requested the following after the initial implementation acceptance:

> a) Would 1024 anyway be enough for our present usecase here.
>
> b) Do a longer integration run for the three sets of center, until the 1% barrier is crossed on both real and imaginary part (using best possible integration setup and sampling, with some exploration runs to see which one give smallest variance) to verify that they match at that level. Also tell me how many sample points are needed to reach that level for each of the three sets of centers.
>
> c) In:
> ```ini
> [cli_settings.global.generation]
> "3dreps" = ["ltd", "cff"]
> ```
> the name "3dreps" is bad as it forces quotes. Use the explicit "three_dimensional_representations" name for that option instead. Update this option name everywhere. Same in the name "3drep" of the "levels". Use "three_dimensional_representation" there instead too.
>
> d )Give me a comparison of generation and runtime for the GL638 example using advanced sampling in the cli epem_tth/NNLO example, between cff and LTD, and also adjust this card to use the setup above of generating both CFF and LTD and the new stability rescue method.

The active generation setting is now
`global.generation.three_dimensional_representations`; each stability level uses
`three_dimensional_representation`. Update active settings, schemas, Python,
completion, diagnostics, export selectors, documentation and cards consistently.
The symbolic diagnostic CLI command `3drep` retains its name. Historical prompt
quotations above remain verbatim and are superseded by these setting names.

Preserve the fast center acceptance test. Run separate reproducible sampling
pilots and extended numerical studies for the same triangle-box fixture, checking
both real and imaginary relative MC errors below 1% for each center strategy and
recording sample counts and agreement. Do not equate the integrator's stopping
criterion for one trained component with this two-component acceptance rule.

Update `examples/cli/epem_a_ttxh/NNLO/gl638_advanced_sampling.toml` to generate
both representations with local UV counterterms from 4D and ordered
Double/LTD, Double/CFF and higher-precision rescue attempts. Preserve its physical
input, complete subtraction and advanced sampling channels. Measure generation
and runtime under matching conditions, retaining accuracy and stability evidence
alongside timings. Determine the scope in which a 1024-argument dependency cap
would cover actual retained calls, separately from the unresolved serializer
limit in stock SymJIT 2.26.1.
