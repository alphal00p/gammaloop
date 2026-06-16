> Historical phase-stack plan and execution record. All results below apply to
> the recorded source checkpoints, including superseded intermediate steps.
> They do not certify the reconstructed reviewed stack; see its
> [review](docs/architecture/raised-energy-cff-stack-review.md) and
> [fresh validation](docs/architecture/raised-energy-cff-stack-review-validation.md).
> Current work follows [CONTRIBUTING.md](CONTRIBUTING.md) and
> [current architecture](docs/architecture/architecture-current.md).

> All final-source local gates pass: 2,084 exact `just test_gammaloop` tests, 114 focused tests, and all 166 scalar LU cases.
> The four signed ddx/ttx LO/NLO acceptances, fresh API/schema/parity, check, Clippy, formatting, and Hakari pass.
> Four fast numerical PySecDec tests were run manually and pass; all 24 physical PySecDec tests remain excluded from automatic local/CI coverage.
> CP symmetrization remains an explicit opt-in with warnings.
> The originally accepted plan is preserved below.
> The accepted research-driven revision after the verbatim addenda supersedes
> the proposed numerator-conjugation modes and their implementation steps.

# GammaLoop phase conventions fix — proposed implementation plan

This is an investigation-backed draft for acceptance. Implementation has not started. At implementation start, copy the accepted plan verbatim to PHASE_CONVENTIONS_FIX.md at the repository root, including these addenda, and set the goal statement below as the active task goal. Keep later implementation evidence in an appended execution log so the accepted plan remains intact.

## Goal statement

Make GammaLoop return the complete textbook Minkowski Feynman-rule amplitude contribution A_Gamma = -i M_Gamma, with a consistent scalar UFO/JSON model, and physical LU cross-sections obtained by correctly conjugating the right-hand amplitude. Establish this through independent scalar contour and phase-space oracles across loop counts, exact and scalar_only numerator-conjugation tests, and signed ddx/ttx LO and NLO acceptances matching the published values without compensating phases or absolute-value comparisons. Only after those conventions are certified, migrate all affected tests, snapshots, examples, and documentation, complete the required checks, and push the finished changes to the current branch name with suffix _phase_fix.

## Findings that determine the plan

1. The shared energy representation is J = integral product[dq0/(2*pi*i)] N/product(D^a), with D = q^2 - m^2 + i0. GammaLoop's existing source-frame and inverse-energy signs correctly recover J. Production then multiplies by (-i)^L/(2*pi)^(3L), whereas the physical Minkowski measure requires i^L/(2*pi)^(3L). The current uncut amplitude therefore differs from the Feynman-rule result by (-1)^L. This statement holds for the audited scalar contour boundary; each subtraction branch must be verified consistently during implementation.
2. The independent tadpole and tadpole-times-bubble contour check in crates/gammalooprs/src/cff/mod.rs:6879 explicitly strips the production (-i)^L factor at line 6946 before comparison with J. This distinguishes loop parity from denominator parity. The three production normalization sites are cff/mod.rs:392, :699, and :828. The real source-frame adapter at :224 and the 1/product(-2E) convention in graph/mod.rs:724 must not be independently flipped.
3. Integrated UV terms have the same historical convention: uv/settings.rs:425 defaults additional_normalization to "-1", and uv/approx/integrated.rs:977 applies it once per integrated loop. The amplitude correction must change this default coherently to "1", including its serialization-default predicate. Preserve forest subtraction signs and actual causal prescriptions.
4. The vendored scalar UFO uses +i*lambda, while the standard interaction L_int = -lambda*phi^4/4! requires -i*lambda. Its full JSON export has +i scalar propagator numerators, as expected. The actively used scalars_2p_3p.json instead has propagator numerators 1; its coupling expression and cached numerical value also use +i*lambda. Both JSON variants and the UFO must be corrected consistently.
5. LU currently retains the full forward-graph CFF loop phase, keeps the final-cut propagator numerators, and multiplies by a real 2*pi cut factor. The right integrated threshold coefficient is conjugated, but the RHS numerator and the remaining RHS virtual-loop measure are not treated as a complete conjugated amplitude.
6. An isolated diagnostic with the existing, unrebuilt target/dev-optim/gammaloop binary selected one symmetric scalar contact graph and one Born cut for each multiplicity n = 2,3,4,5. With no UV or remaining threshold subtraction, the phases were respectively -i,+i,-i,+i at forward loop counts L = 1,2,3,4. This corroborates the source-derived phase (-1)^(n+1)*i. It is diagnostic evidence, not fresh-build acceptance or an MC integral. Artifacts: /var/folders/xx/c9yqv8gs07b7n6njwdjymg4w0000gn/T/gammaloop-phase-born-55povae7/.
7. The present scalar LU matrix's "no_numerator" cases retain complete UFO Feynman rules; the label means no additional momentum numerator. A graph-level DOT num="1" also does not clear edge and vertex factors. Reuse Graph::with_global_numerator_only for a genuinely complete unit numerator.
8. RHS ownership first becomes insufficient when CrossSectionCut is reduced to LuCutSelection: physical cut alternatives survive, but the side-dependent numerator information does not. UV forest values are also computed and shared before the current cut loop. Applying conjugation only to the final numerator would therefore leave UV counterterms incorrect.
9. This checkout has no method literally named complex_conjugate. Existing spenso_conj, scalar/tensor network mapping, and Dirac simplifiers provide the required operations. dirac_adjoint additionally inserts gamma0 and transposes open spin-line endpoints; it must not be substituted blindly for elementwise conjugation of an amputated RHS tensor.
10. Current ddx/ttx acceptances contain phase-selecting probes, norm-based LO/NLO comparisons, and manually supplied +/-i generation prefactors. Some individual graph comparisons already retain signs, but their inputs still contain these compensating conventions.

## Convention to impose

Use metric (+,-,-,-), positive-energy physical cuts, and the ordinary measure product[d^4k/(2*pi)^4]. Return the complete amputated Feynman graph, including its vertex and propagator numerator factors, as A_Gamma = -i M_Gamma in the user's notation. The operational anchor is the rule itself: a phi^4 contact graph is -i*lambda and an internal scalar propagator is i/(q^2-m^2+i0). Some textbooks name the complete graph i*M instead; no numerical ambiguity remains once this anchor is fixed.

For mixed scalar interactions with occupation numbers n_a, use L_int = -lambda*product(phi_a^n_a)/product(n_a!). Their vertex is -i*lambda without another factorial; feyngen independently owns diagram and identical-state symmetry factors.

For a literal complete unit numerator and P = sum of denominator powers, a Wick-rotatable integral obeys I_M = i^L*(-1)^P*I_E. I_E is real and is positive for a convergent positive Euclidean scalar integrand. A no-threshold condition alone does not guarantee UV convergence or the sign of a renormalized finite remainder. The phase axis still alternates between real and imaginary with loop parity. In an ordinary four-point phi^4 topology P = 2L, so the unit-numerator sequence at L = 0,1,2,3,4 is 1,i,-1,-i,1. Physical UFO amplitudes include additional vertex and propagator phases; they must be tested separately.

For LU, impose d sigma_cut = positive physical normalization * A_left * conjugate(A_right) * dPhi_n. For identical scalar Born halves this is nonnegative, and strictly positive at nonzero interior points. There is no residual phase depending on the forward loop count. Do not impose positivity on individual interference, UV-counterterm, threshold-counterterm, or gauge-projector graph contributions.

The source-derived LU conversion to certify is: after fixing the shared contour normalization to i^L, retain each cut-line propagator numerator exactly once, conjugate RHS numerator factors and the RHS virtual contour, and replace the current 2*pi LU residue factor by 2*pi*i. In this representation the additional RHS contour factor is (-1)^L_R, where cut legs are external and L_R counts RHS virtual loops. Thus the combined LU conversion is i*(-1)^L_R relative to the physical forward contour. This is an analytic residue-to-cut conversion, not a process-dependent fitted prefactor. Verify it with an independently sewn amplitude and phase-space measure before making it the production contract.

References:
- David Tong's Feynman rules and phi^4 example: https://www.damtp.cam.ac.uk/user/tong/qft/qfthtml/S3.html
- Local Unitarity: cutting raised propagators and localising renormalisation, arXiv:2203.11038, especially Eq. (7.1) and the signed MSbar results table on PDF page 55: https://arxiv.org/pdf/2203.11038

## Implementation sequence

1. Record the accepted plan and establish the branch and goal.
   - Confirm the current clean/dirty state and preserve any user work.
   - The inspected current branch is codex/raised_energy_cff_wip_optimized. Create codex/raised_energy_cff_wip_optimized_phase_fix from its current tip.
   - Write the accepted plan verbatim into repository-root PHASE_CONVENTIONS_FIX.md, including the goal statement, original request, and follow-up instructions.
   - Set the goal statement above as the active task goal.
   - Record baseline revision, relevant settings, signed source diagnostics, and reproducible commands. Keep the existing test expectations until the independent convention checks establish the intended replacements.

2. Make the scalar model a textbook oracle.
   - Change assets/models/ufo/scalars/couplings.py to -i*lambda.
   - Regenerate/update assets/models/json/scalars/scalars.json and scalars_2p_3p.json, including symbolic expressions and cached coupling values [0,-1] at lambda=1.
   - Restore scalar propagator numerators +i in scalars_2p_3p.json while preserving its restricted interaction inventory.
   - Keep real lambda=+1 restrictions and zero-width benchmark settings; do not compensate by negating lambda.
   - Explicitly import this repository's UFO when regenerating, since the installed loader has its own independent bundled scalar model.
   - Rebuild compile-embedded JSON assets, then compare direct UFO and JSON imports, contact amplitudes, and an exchange tree with an internal propagator. Check mixed-species factorials and graph symmetry factors.
   - Do not apply a universal potential-coupling sign change to the separate scalar_gravity model.

3. Correct the amplitude measure at its existing owner.
   - Change all three physical CFF normalization paths from (-i)^L to i^L.
   - Change GammaLoop's matching Vakint additional normalization from -1 to 1, including default serialization and persisted-settings expectations.
   - Leave the validated real source-frame conversion and energy factors intact.
   - Verify ordinary, raised, exact-source, and stored-production paths against analytic energy residues.
   - Check numerator-one and full-UFO graphs independently, including zero-loop contact/exchange anchors and one through four loop scalar contours. Use convergent examples for actual MC targets.
   - Check local UV and integrated addback terms together, including contracted subgraphs whose remaining CFF loop count differs from the original graph.
   - Preserve and verify amplitude threshold discontinuities and causal signs with an independent contour/discontinuity example.

4. Introduce cut-aware numerator conjugation before UV construction.
   - Add global.generation.complex_conjugation = "exact" | "scalar_only", with exact as the default, on the existing GenerationSettings owner.
   - Reuse CrossSectionCut::amplitude_side_subjects for side membership. Cut spin sums/projectors and globally supplied complete-forward prefactors remain outside the amputated RHS conjugation.
   - Preserve original edge/vertex ownership through the transformation, numerator maps, Taylor terms, and UV reconstruction.
   - Use immutable generation inputs per distinct RHS conjugation signature, with caches scoped/keyed by that input. Feed the transformed factors into both direct-3D and projected-4D UV construction, including integrated UV terms.
   - Simplify the assembled RHS network as well as local factors, so identities crossing factor boundaries can cancel.
   - Carry RHS information through physical cut alternatives and raised-residue construction. Share a representative only after certifying equivalent transformed payloads; do not split a combined raised pole in a way that duplicates its Laurent residue.
   - Rebuild or certify energy-rank bounds and numerator maps after transformation rather than assuming they remain valid.
   - A manually supplied complete global forward numerator has no inferable left/right factorization; retain its documented full-forward meaning. Local factors or explicit side metadata are required for RHS-specific treatment.

5. Implement and certify the two modes using existing tensor APIs.
   - Exact: faithfully conjugate the entire RHS tensor, including couplings, explicit i, color, and Dirac structures. Select the existing whole-tensor operation according to the actual endpoint/sewing convention; verify it against direct componentwise conjugation. Apply simplify_gamma_conj, chain-aware simplify_gamma, gamma0 and metric/Schoonschip simplification. Respect the existing four-dimensional limits of gamma0 identities. Avoid applying Dirac-adjoint endpoint transpositions twice.
   - Scalar_only: parse the factorized symbolic network with its existing scalar/tensor separation, conjugate scalar coefficients, and leave tensor leaves/index structure unchanged. Remove conjugation wrappers only for couplings whose resolved numerical imaginary part is exactly zero at generation. Conjugate genuinely complex and unresolved couplings; do not give UFO symbols global Real attributes.
   - Ensure scalar conjugates lower to supported evaluator operations; spenso::conj must not survive as an unsupported opaque scalar function.
   - Test real, imaginary, and general complex couplings, nested scalar factors, sums, powers, compact chains, gamma5/projectors, open spinor lines, color generators, and epsilon structures.
   - Exact must remain correct under subsequent numerical model updates. Scalar_only is specialized to generation-time coupling reality; record that assumption and require regeneration when an update invalidates it.
   - Check defaults, overrides, schema, serialization/state round trips, SHOWDEFAULTS, completion, and compiled versus interpreted evaluation. Regenerate existing integrands/states for the new convention.

6. Complete the LU phase and causal-conjugation conversion.
   - Implement the certified i*(-1)^L_R residue/contour conversion at the existing cut boundary, retaining cut-propagator numerator factors once.
   - Count RHS virtual loops with cut legs external, and account consistently for loops moved into integrated UV terms.
   - Apply the same cut context to bare, local UV, integrated UV, and left/right/iterated threshold branches.
   - Retain the existing opposite integrated threshold prescriptions and prove that they are not conjugated a second time.
   - Add a one-virtual-loop RHS example and its left/right mirror: Born-only cuts have L_R=0 and cannot validate this part of the convention.
   - Check that conjugation does not change cut eligibility, event association, group multiplicities, or the physical state-sum normalization.

7. Establish independent scalar LU acceptance across multiplicities.
   - Generate symmetric single-contact Born halves with n = 2,3,4,5,6 final scalars, giving forward loop counts 1 through 5. The existing scalar model already supports the required interaction degrees.
   - Restrict the process/vertex inventory so exactly one graph and one physical Born cut are retained; assert these facts instead of depending on graph names.
   - At deterministic interior points, compare the complete cut weight to an independently evaluated |A_Gamma|^2 times the LU phase-space Jacobian. Check a vanishing imaginary component and positive real value, in both conjugation modes.
   - Run MC integrations for each multiplicity against independent massless n-body phase-space normalization:
     Phi_n(s) = s^(n-2) / [2*(4*pi)^(2*n-3)*Gamma(n)*Gamma(n-1)].
     For one scalar of mass M decaying to n identical massless scalars, the rate oracle is |g|^2*Phi_n(M^2)/(2*M*n!) when the physical flux is enabled and unit conversion is disabled.
   - Verify flux and symmetry conventions separately and include a two-incoming-scalar example so a decay-only normalization cannot conceal a scattering error.
   - Add symmetric multi-vertex Born trees and a massive smooth case to exercise RHS internal propagators and momentum factors.
   - These fixtures have no remaining threshold after the Born cut. Their positive-weight acceptance is independent of local UV/threshold subtraction and is followed by the separate all-subtractions-enabled tests.

8. Require signed ddx and ttx physical acceptances.
   - Run independent acceptance tests for all four channels: e+ e- > a > ddx, e+ e- > a > ttx, a^star > ddx, and a^star > ttx, after rebuilding and regenerating their integrands.
   - Each channel must explicitly check both components of its own LO and NLO result: positive real LO with zero imaginary part, and the published signed real NLO correction with zero imaginary part, within the appropriate numerical/statistical tolerances. A ratio check or a phase inferred from another channel does not satisfy this requirement.
   - Require positive real LO and compare the signed real NLO correction, including its published negative and positive graph components where the current tests provide those oracles.
   - Remove manual +/-i generation prefactors, extraction of -Im(LO), norm/hypot-based physical-result comparisons, and dominant-phase selection. Keep mathematically justified real spin averages, flux factors, symmetry multiplicities, Eq. (7.1) conversion, and units.
   - A phase fix must not absorb a normalization or graph-multiplicity discrepancy. Record and derive each retained real normalization.
   - In the paper's direct-current normalization and MSbar scheme, check:

     | Channel | LO | NLO vertex | NLO self-energy | Total NLO correction |
     | --- | ---: | ---: | ---: | ---: |
     | ddx, Ecm=400 GeV | 0.5031049 | +0.0503926 | -0.0314956 | +0.0188970 |
     | ttx, Ecm=600 GeV | 2.876302 | +0.348276 | -0.146756 | +0.201520 |

   - The entries are the published rows with their stated multiplicities; do not add another self-energy factor blindly.
   - For the lepton processes, apply the established Eq. (7.1) real conversion and picobarn units; compare signed absolute LO/NLO values as well as ddx NLO/LO = alpha_s/pi and the ttx published ratio.
   - Keep the fully-MSbar ttx scheme used by the existing acceptance.
   - Compare complete complex values across orientation-local direct 3D, explicit-sum direct 3D, and projected local 4D with UV and threshold terms enabled. Require imaginary integrated results to be consistent with zero with precision/statistical criteria; no post-hoc real projection or absolute value may hide a phase.
   - Preserve ddx UV-mass and localization-scale independence, cancellation of total renormalization-scale dependence, and opposite graphwise scale slopes.
   - Use exact as the physical acceptance default. Compare scalar_only where its algebraic assumptions are justified; it must not become a workaround for a failure of exact.

9. Migrate the affected tests only after independent certification.
   - First establish the scalar contour, scalar Born phase-space, exact-conjugation, and signed ddx/ttx checks.
   - Then update every affected phase-sensitive expectation, snapshot, run card, benchmark helper, and architecture statement to the new convention. The user's follow-up explicitly authorizes these changes without repeated permission requests.
   - Separate changes caused by the physical loop measure, scalar vertex sign, restricted-model propagator repair, and RHS conjugation; provide derivations for new expected signs.
   - Correct the misleading historical cubic-contour comment and retain/move existing comments rather than deleting them.
   - Replace obsolete magnitude acceptances with signed comparisons while retaining meaningful tolerances and coverage. Do not weaken tests merely to make them pass.
   - Update docs/architecture/architecture-current.md and related convention/state documentation; update differential_lu.md if event or cut-group behavior changes.

10. Complete validation and push the finished branch.
    - Format and run cargo check before compilation, then meaningful focused unit and integration tests via the repository's just/nextest workflow.
    - Run the full 166-case just test_LU_scalar_xs matrix, including slow cases, under its existing 30 GB watchdog.
    - Run the four signed ddx/ttx acceptances and relevant amplitude, numerator, CFF, UV, threshold, event, settings, and model-import regressions.
    - Run clippy and required broader repository checks. Report any resource stop as incomplete, not a pass.
    - Review the diff for concise idiomatic changes, existing helper reuse, preserved comments, and any unaccounted process-specific phase adjustment.
    - Append actual commands, signed numerical results and errors, source references, and any remaining limitations to the execution log.
    - Commit with the complete user request in the description as required by CONTRIBUTING.md, push codex/raised_energy_cff_wip_optimized_phase_fix to origin, and report the branch/commit and test evidence.
    - Mark the active goal complete only when the required work and publication are finished.

## Addendum A — original user request, verbatim

```text
Alright, time to address and fix once and for all the phase issues conventions we have in gammaloop.

First we need to plan this by carefully and diligiently identifying what is currently the phase convention of our amplitude type of calculations in gammaloop and of our cross-section calculations in gammaloop.

Then we need to decide what is the actual phase convention we want to *impose*. For this the guidelines are as follows:

a) For amplitudes, the phase is of course not observable. Yet, I would like our convention to be such that if feyngen were to build an amplitude graph from a phi^4 UFO model, then our result for the MC integral of that graph is directly what enters the textbook "scattering amplitude", i.e. the quantity typically denoted `-i M`, ready to be squared to build the squared matrix element forming the cross-section.
Also reflect on what this imply for resulting phase returned by gammaloop when computing a topology loaded from a long graph where the user forced the numerator to be exactly one (i.e. a "scalar" loop), and supposing it has no threshold, so that the CFF integrand has a uniform phase and no integrated threshold CTs, so no source of imaginary part, For example would the phase in such an example keep switching at various loop count with the phase choice yielding `-i M` like requested above?

b) For the cross-section mode when handling forward scattering graphs, the desired phase requirement is much simpler and clearer: the cross-section must be real and positive when considering a case with no threshold once each Cutkosky cut is imposed with a symmetric left and right graph. (e.g. when reproducing the Born cross-section of x > y (even with y having many particles, so the forward-cross-section is multiloop)). Indeed, each of these symmetric Born-like cut corresponds to | A\_\Gamma |^2 which must corresponds to the norm of an amplitude graph and must be positive and real so the cross-section is positive and real. This works only processes with external scalars (otherwise for e.g. vectors, the unphysical spin-sum of -g^{\mu\nu} means that this property may not be realised per graph), but it's enough to fix the phase-choice.&#x20;

c) To get the property b) on must of course consider properly the complex-conjugation on the right of the cutkosky cut. This is currently properly applied only for the imaginary prefactor of the integrated threshold CT (i.e. the causal prescription of the propagators to the right of the cut). But we need to apply it to the numerator too. When using the currently implemented "complex\_conjugate()" on a tensor network, it is a bit over-eager as it really applies a CP transformation to the tensor network. It is not really wrong but it needlessly complexifies the algebra. However, now that we support the schoonship notation, it may be easier to do. I want you to add a new generation mode "complex\_conjugation" which has two modes "exact", "scalar\_only". In the "exact\_mode" (default) you faithfully call "complex\_conjugate" on the whole tensor numerator on the r.h.s of the cut (i.e. CP) and try to simplify the complications that appears before proceeding (so any (gamma^0)^2 should be set to 1 for instance) and then proceed.  In "scalar\_only" mode, you only complexify the scalar coefficients of that r.h.s tensor, leaving all couplings whose current numerical value is real-values (at the time of the generation) untouched.
You must then verify that the phase requirement I mentioned for the LU cross-section result is now exactly met for various loop counts of a scalar LU cross-section (similar to the LU check of the test matrix but engineering only symmetric graphs with a single "born-like" cut of progressively higher multiplicity, where we know it must yield a positive and real cross-section no matter the loop count of the forward-scattering graph.

d) You'll be working and pushing to the branch "raised\_energy\_cff\_optimization\_phase\_fix" (not exact name, essentially I want the same as current name, but with suffix "\_phase\_fix")

Upon starting implementation of the plan, specify to write it verbatim in the file PHASE\_CONVENTIONS\_FIX.md, including a goal statement, and setting it to yourself as the goal.
In that plan, had an addendum that will include this user prompt here verbatim.
```

## Addendum B — follow-up user instructions, verbatim

```text
Ah and you can override the repo  requirement to ask for permission to make changes to the tests and things like this because there may be many here.
Also ask questions if anything is unclear in the plan after investigation of course.
```

```text
Continue elaborating the plan, but also make sure the conventions of the UFO model for the scalar theory (sign/phase of the prefactor of interaction feynman rules and sign/phase of the scalar propagator numerator is also matching text-book expectations).
```

```text
Ah and also add to the plan than then you must run the ddx and ttx tests and verify that you get a positive real LO cross-section and also an NLO correction matching our published numbers, but this time not in absolute value (like the test currently tests) but with the correct sign without any corrective ad-hoc measure.
```

```text
Ah and then of course add to the plan that when you're certain the phases are fixed according to my guideline, then all tests sensitive to that must be updated. I see you want to change the scalar model vendored in to match textbook conventions, this is fine indeed, but make sure then to update both the UFO and the json model there too.
```

```text
Add to the plan that you should check both e+ e- > a > ddx / ttx and a^star > ddx / ttx results (incl. phase) in the tests.
```

## Addendum C — implementation instructions, verbatim

```text
Implement the plan.
```

```text
(continue as planned, but you can safely assume that a model parameter declared "real" in the UFO will always remain real, the output would be invalidated otherwise, so at generation time it should be clear whether complex conjugation hits them or not (i.e. giving those symbols the attribute "real" should simply do the trick already actually))
```

```text
But actually spin out a research agent too in order to really make a deep study of whether doing the complex-conjugation on the spin tensors on the r.h.s. of the cutkosky cut is really the correct things to do (i.e. it's not sort of already done somehow due to how the feynman rules and edge orientation works)
```

```text
Wait this is important explain why it is wrong
```

```text
OK then remove the request for the two "complex conjugation modes" then just keep the code doing one thing, but the *right thing* as you described. Do another deep research to really pin down what is the right thing to do there (also if we had complex-values coupling like off-diagonal CKM or in the complex mass scheme).
```

```text
Ok but if it comes out that complex-valued couplings are difficult to support, then simply make sure that the code errors nicely if any are found in a cross-section and implement your solution for real-valued couplings.
```

```text
Continue as planned, but also pay attention to the cut edges in cross-section, indeed one must be careful to overall obtain the correct interplay between the origin `i` in the propagator numerator and the corresponding scalar phase factors used within gammaloop when taking the cut of that propagator.
Anyway, in the end it must all play out to give the requirements we asked for, with principled handling of every source of phase.
```

## Accepted research-driven revision

The revised goal retains the full physical amplitude, scalar-model, signed cross-section, regression and publication requirements, but removes both proposed complex-conjugation modes. Use one marked-forward prescription on the already inverse-oriented Hermitian UFO numerator; do not apply an additional scalar, tensor or Dirac adjoint to that raw RHS.

The independent second audit covers intrinsic complex CKM phases, the complex-mass scheme, cut-line numerator factors, and the complete marking/contour ledger. Its implementation contract and primary references are documented in [phase-conventions.md](docs/architecture/phase-conventions.md):

- Preserve existing oriented Hermitian-partner vertices, cut spin sums, fermion statistics and factorized numerator ownership. The unnecessary per-cut numerator/UV cloning and conjugation APIs are removed.
- Right-side vertex and virtual-propagator markings together with virtual-contour reversal yield `(-1)^(V_R+I_R+L_R)=(-1)^C_R`. Count the physical RHS components with cut hairs retained, and require equal parity among physical cuts sharing a raised residue.
- The complete CFF/LU conversion is `2*pi*i*(-1)^C_R`, hence `-2*pi*i` for connected RHS. Cut propagator numerator `i` factors remain in the common source and must cancel exactly once against the appropriate cut scalar factor. The physical whole-propagator replacement is `i/(q^2-m^2+i0) -> 2*pi*delta_+(q^2-m^2)`.
- Multiply bare, local/integrated UV and threshold terms by the same group factor; retain opposite left/right threshold prescriptions. No process-dependent compensating phase is allowed.
- Intrinsic complex CKM couplings in a Hermitian model are already paired as `V` and `conj(V)` by the inverse vertex assignment. Add generated-graph component certificates for a charged chiral scalar and complex CKM current; do not reject conventional imaginary UFO Feynman-rule factors.
- CMS requires complex poles, derived parameters and counterterms beyond these stable-particle cutting rules. GammaLoop's current real-energy implementation must report a clear error for actual complex masses at generation and after model updates. Ordinary width metadata does not activate CMS. The user's fallback authorizes explicit errors for genuinely unsupported complex conventions.
- Replace two-mode duplicates in the new scalar acceptance fixtures by one physical path. Keep all independent contour, cut-line, multiplicity, virtual-mirror, signed ddx/ttx and full-regression requirements.

## Execution log

- Implementation started from clean revision `ac99d32d194254f0c8931873d0a32ebd7d27f82f` on `codex/raised_energy_cff_wip_optimized_phase_fix`. The accepted plan above was copied verbatim and its goal was activated before implementation.
- Direct Python import of this repository's scalar UFO and both JSON inventories confirms the vertex `-i*lambda` and scalar propagator numerator `+i`; the full and restricted interaction inventories are preserved. The installed canonical UFO loader stalls in `symbolica_processing.py:56` while repeatedly replacing square roots in parameter expressions, before loading couplings or propagators, including for a phi4-only inventory. Its diagnostic was stopped; no successful full canonical-loader round trip is claimed.
- Initial targeted `cargo check` passed. Independent scalar phase-space and contour acceptance tests are being compiled; numerical certification and broader regression results will be appended when complete.
- A separate read-only research agent is auditing RHS spin tensors, forward edge orientation, crossing and cut sewing against independent amplitude conjugation. The proposed exact spin operation remains provisional until this audit and component certificates settle whether any conjugation is already implicit.
- The independent spin audit falsified the proposed raw-RHS adjoint: actual SM charged Goldstone inverse vertices already contain the opposite chiral projector. A direct explicit Weyl-matrix sewing gives 21, whereas adjointing the inverse vertex again gives zero. An actual complex charged-current vertex pair gives 105 with correct sewing, but scalar conjugation or the extra adjoint of the inverse vertex gives -63+84i. These are independent component calculations, not yet generated-integrand acceptances. The raw forward representation must instead obey marked-vertex and anti-causal-propagator rules; see Anselmi, arXiv:1606.06348, sections 2.2–3. The originally proposed tensor-conjugation implementation will not be shipped as physically exact. A clarification about replacing the now-redundant generation modes is pending.
- The first scalar test execution stopped before physics evaluation because the new fixture used the Rust field name `lu_h_function` rather than its serialized TOML name `h_function`. Both affected fixtures were corrected. A subsequent `cargo check -p gammalooprs -p gammaloop-api -p gammaloop-integration-tests --tests` passed; independent contour and component unit tests are compiling.
- The first focused unit run executed 18 tests: 15 passed, including the independent amplitude contour sequence, powered-pole physical measures and Vakint default. Three failures were in malformed gamma shorthand fixtures of the discarded conjugation-mode prototype; those prototype APIs and tests have been removed, rather than retained as physical acceptance. Actual generated-UFO sewing tests replace them. Log: `/tmp/gammaloop-phase-unit-tests.log`.
- The second literature audit was performed with the deep-research workflow. `update_plan` is not exposed by this environment, so discovery, follow-up and synthesis were tracked in the task's structured plan and this execution record. Primary claims were independently checked against Anselmi (2016), Denner–Lang (2015), Denner–Dittmaier (2006), and Bauer et al. (2012). The repository Markdown artifact was structurally verified against its canonical source; no visual rendering is claimed.

## Addendum D — derivative-coupling audit request, verbatim

```text
Also launch a separate agent to do a deep investigation of derivative couplings in the SM. For SMQCD the only one is the three-gluon vertex, and my concern here is that I believe (but to be checked) that our convention is "all-incoming" for the meaning of a $$p^\mu$$ in the feynman rules, but I am not convinced that the way the three-gluon feynman rules, as well as the QCD ghost one (and others beyond QCD) is currently encoded in our SM model actually fits this convention in its choice of sign (e.g. vs literature).
This is delicate as we're talking about signs here, so the online literature authoritative result must be triple-checked, and then really dig into how our current momenta in the feynman rules is being interpreted in gammaloop and confirms that it is consistent. All such derivative couplings must be checked this way and if there is any mismatch requiring a change in the SM model, make sure to escalate those changes to me so that I can decide them to greenlight them or not.
```

### Added audit requirement

A separate research agent will catalogue every momentum-dependent SM UFO and JSON Lorentz structure, trace its momentum labels through GammaLoop's oriented graph assignment and numerator generation, and triple-check the aligned signs against authoritative primary sources and independent action-level derivations. The three-gluon, QCD ghost, electroweak gauge, scalar/Goldstone and ghost derivative rules are in scope. Any proposed SM-model correction requires the user's explicit approval after presenting the precise discrepancy, convention ledger and minimal reproducer; research and independent diagnostics may proceed without that approval.

```text
Ah and also ask this agent to also double-check the sign of the momentum in the numerator of the fermion propagator (and anti-fermion). However there, I'm more confident that it's also correct (but it's tricky to get right!).
```

The independent derivative audit also covers fermion and antifermion propagator numerator momentum signs, distinguishing fermion flow, graph-edge orientation and physical momentum, and checking their consistency with crossing and cut spin sums. The same approval boundary applies to any proposed SM-model correction.

```text
Wait what are the two SM JSON variants? There should only be one no?
```

```text
Ah ok, well those should match exactly (only difference is indices wrapping). Actually I would just remove the sm\_wrapped version.
```

The user explicitly authorized removing `assets/models/json/sm/sm_wrapped_indices.json`; that obsolete duplicate is removed. The standard model loader selects `sm/sm.json`, and source searches found no explicit consumers of the wrapped file. The independent physics audit proceeds against the canonical SM UFO and `sm.json`; approval remains required for any proposed correction to their physical rules.

```text
Yes this ALOHA mismatch is exactly what I'm afraid about, but look, in the end we don't use/care about aloha, so the only thing that you need to do is to really compare the textbook correct rules vs the ones in our SM model and its interpretation within gammaloop, to make sure it lines up. ALOHA can do whatever it wants, we don't care I think.
```

The derivative audit's acceptance criterion is exclusively the independently established textbook/action rule versus the actual SM model interpreted through GammaLoop's momentum, particle, color and spinor assignments. ALOHA behavior can explain historical provenance but cannot establish correctness or justify a model correction.

### Validation progress

- `cargo fmt --all -- --check` and `git diff --check` pass; all changed JSON and TOML files parse.
- `CARGO_BUILD_JOBS=2 cargo check -p gammalooprs -p gammaloop-api -p gammaloop-integration-tests --tests` passes after applying the four signed acceptance changes and correcting the new tree fixture. Log: `/tmp/gammaloop-phase-signed-check.log`.
- `CARGO_BUILD_JOBS=2 just --set ci_cargo_profile dev-optim clippy` passes across the workspace and all targets, with no warnings, in 1m42s. Log: `/tmp/gammaloop-phase-clippy.log`.
- The rebuilt 23-test scalar/contour/sewing suite is in progress. Its newly added tree fixture initially referenced nonexistent `mass_scalar_0`; that particle is intrinsically massless. The assignment was removed without changing the physical oracle. The failed fixture must be rerun. Log: `/tmp/gammaloop-phase-physical-tests.log`.
- The four signed ddx/ttx acceptance executables are being rebuilt. Their physical-result comparisons use the signed real component, separately constrain the imaginary component, and retain only independently derived real averages, flux and unit factors.

```text
Continue as planned, but realise that for this phase-test you can generate a single orientation (with that cut contained), it does the job equally well to test the phase.
```

For phase-only cut diagnostics, a single orientation containing the physical cut is sufficient. Absolute-integral comparisons must use the independently derived weight of the selected contribution. The uncut unit-numerator MC fixture selects one acyclic orientation at every loop count, including four loops: each equal-mass vacuum triangle contributes one sixth of its full integral, giving `[-i/(192*pi^2*m^2)]^L`. Its pointwise residue and MC value are both checked against that selected-contribution oracle.

- The first rebuilt suite completed with 19 passes and four outstanding cases: the corrected tree fixture, the intentionally stopped expensive uncut MC fixture, a DOT mass-symbol namespace fixture, and an unprojected color identity in the complex-CKM component fixture. All scalar Born multiplicities two through six passed pointwise and MC acceptance; massive decay, two-incoming scattering, exchange-tree sewing, both virtual-loop mirrors, cut-propagator residue accounting, independent scalar contour tests and the Vakint default also passed. The outstanding fixture corrections and MC rerun remain required.

```text
but even for the four-loop uncut topology you can select a single orientation.
Also make sure that you propagate any change you do to all generation modes and pieces (UV CTs, local and integrated threshold CTs, etc..) but ok this is likely automatic.
```

The implementation uses the common physical contour measure and LU group factor, so their changes apply to bare terms, local/integrated UV terms and local/integrated/iterated threshold counterterms. The remaining validation explicitly checks these owners and opposite left/right threshold prescriptions across generation modes.

```text
I think the diagram generator *does* include a negative sign for closed ghost loops! I implemented it generically on any anti-commutating field, and this includes both ghosts and fermions I believe!
```

The user correctly recalled commit `6e92007a1d3e1a1202de590f01482224e1c21df6`, which introduces a generic anticommutating-field predicate and applies it to loop statistics and external ordering. That commit is reachable on `origin/alpha_output_v2`, but is not an ancestor of this checkout. The pre-change executable generates a two-ghost-edge gluon self-energy with `overall_factor_evaluated = 1`, confirming that the historical behavior is missing here. The correction must restore the generic statistics semantics without adding a second sign; a production generation regression will supplement the topology counter tests.

```text
Ok then indeed pull that upstream change here already.
```

The user authorized porting the upstream anticommutating-field change to the current repository layout, including its loop-filter naming and external ordering semantics. Obsolete cache files and historical layout artifacts are excluded from the port.

```text
Ok but also simplify all combinations of `P_L-P_R` (or `P_R-P_L`) so as to explicitly write them int the `gamma_5`-free form as it would otherwise cause problems in the d-dimensional regularization when doing d-dimensional gamma algebra to build vaking input.
```

The user approved the four local SM corrections in both the UFO and canonical JSON. A clarification is pending for the additional algebra request: projector sums equal the identity, while differences equal a genuine pseudoscalar `+/-gamma_5` and cannot be replaced by a gamma-five-free single-vertex expression. Any cancellation after contraction must respect the existing dimensional prescription.

```text
Yes sorry, simplify `P_L+P_R` without gamma_5 but then for the differences, write them explicitly as gamma_5!
```

```text
So essentially use the Vector vs axial decomposition, and not the chiral weyl decomposition.
Also does that mean then that you deemed our three-gluon current derivative coupling sign to be valid?
```

The clarification resolves the pending algebra question: use Identity/Gamma5 for scalar/pseudoscalar vertices and Gamma/Gamma*Gamma5 for vector/axial vertices. The canonical UFO and JSON now contain no chiral projectors in their Lorentz structures. Independent exact checks confirm all 22 structures agree between representations, 768 matrix entries match the intended coefficients, and all 26 Goldstone action identities hold with symbolic complex CKM. Existing dimension restrictions on gamma-five algebra remain in force. The current three-gluon sign is valid with GammaLoop's incoming momenta and actual transposed color-generator ordering: the independent quark–antiquark to two-gluon Ward check vanishes for the current sign and fails for its reversal.

- The rebuilt signed acceptance run passed all four e+e-/photon ddx/ttx LO and NLO tests, as well as the single-orientation uncut MC integrals through four loops, runtime complex-mass check and complex-CKM sewing certificate. Its two fixture failures were an empty-tree LMB for a contact graph (now repaired at the existing spanning-tree owner) and a provisional ghost-sign test replaced by the upstream generic anticommutating-field regression. Log: `/tmp/gammaloop-phase-signed-tests.log`.
- The common LU factor now also owns computed UV-forest export normalization. Exports consume already-finalized spatial measures, so their extra spatial divisions were removed. Exact production/export comparisons and independent Born normalization tests cover both UV orchestrators.
- Above-threshold virtual-mirror diagnostics isolate a possible absolute causal-sign error in the single integrated LU threshold coefficient. Local subtraction and shared CFF factors are excluded by CT-on/off comparisons. A fresh analytic triangle regression will be run before changing that coefficient; opposite left/right signs alone do not certify the correct boundary value.

- The fresh 14-test run completed in 4.25 seconds after a 15m41s optimized rebuild: 11 passed. This certifies the textbook contact/exchange amplitudes, scalar MC phases through four loops, negative-real absorptive part of the above-threshold amplitude bubble, generic ghost loop/exchange statistics and API filter, singleton LMB conservation, explicit SM vector/axial model basis, complex-CKM sewing and existing nonempty computed exports. Log: `/tmp/gammaloop-phase-causal-tests.log`.
- Three new diagnostics need another run. The triangle fixture incorrectly enabled runtime threshold subtraction below threshold while omitting its generated threshold data; its runtime flag now matches generation. The charged Ward fixture mixed radical representations, so it now uses homogeneous rational on-shell spinors and checks all generated matrix entries before contracting. The new amplitude/Born computed-export comparison reached a Symbolica polynomial-division failure; its first differing stage is under investigation. None of these failures is recorded as a physical acceptance pass.
- The CP optimization guard is implemented at the model's coupling-dependency owner, with one persisted generation flag revalidated by direct integrand warm-up. It distinguishes intrinsic complex parameters from ordinary imaginary Feynman-rule factors and allows ordinary Hermitian sewing. Fresh generation/runtime tests remain required.
- The separate W/Z gauge mismatch and physical-cut projector boundary are documented in `docs/architecture/sm-conventions-audit.md` and the review proposal `/tmp/gammaloop-phase-plan/SM_GAUGE_CONSISTENCY_PROPOSAL.md`. The user has been asked whether to include that separate coordinated repair; no W/Z model propagator change has been applied. The requested scalar and photon-mediated QCD acceptances are unaffected.
- A native debugger localized the export diagnostic failure precisely to `export.rs`'s test-only `.together()` comparison: production generation, export recomputation and DOT reparse all succeed. The comparison now uses exact unfactorized rational polynomials over the existing `Q_I` field, retaining the independent Born normalization oracle. Trace: `/tmp/gammaloop-phase-plan/computed_export_lldb_crash.log`.
- Formatting and the complete three-package test-target `cargo check` pass after the CP guard and fixture fixes (7.81 seconds; `/tmp/gammaloop-phase-causal-rerun-check.log`). The focused rerun is compiling; the absolute LU integrated-threshold coefficient is still unchanged, so the corrected triangle fixture can establish the expected failing absorptive-sign case before a production correction.

```text
I'm not sure what you mean with: "The diagnostic now uses an exact unfactorized rational-polynomial comparison over \\(\mathbb Q(i)\\), preserving the phase check without that failing factorization step." but make sure you *NEVER EXPAND* the numerator.
```

The user forbids numerator expansion, including temporary diagnostic copies. The rational-polynomial comparison proposed above is superseded and must be removed: it can expand a test expression even though it avoids factorization. The in-progress rerun was stopped during compilation before executing that diagnostic. All new comparisons must preserve the graph numerator's factors, using structural equality or direct signed evaluation without polynomial conversion, `together()` or distributive expansion of the numerator. Production phase ownership is unchanged; the LU causal-sign assertion remains pending.

- The export diagnostic now compares stored factorized expressions structurally, or evaluates their leaves at three exact nonzero assignments if their stored forms differ. The latter checks signed values at those points, not a symbolic identity. Its independent Born normalization oracle remains structural. The new SM component checks and the scalar matrix's two existing edge-numerator checks likewise no longer expand expressions; the edge certificate requires exact factor cancellation and structural reconstruction.
- Formatting and the complete three-package test-target `cargo check` pass with these changes (`/tmp/gammaloop-phase-factorized-check.log`). The focused causal/export/sewing/CP-guard rerun is compiling under the 30 GB watchdog (`/tmp/gammaloop-phase-factorized-tests.log`). No numerical result from this rerun is claimed yet.

- The fresh seven-test rerun completed after a 15m01s build: the generated charged-scalar sewing and charged Ward/ghost-momentum checks pass without expansion. The triangle now reaches the discriminating above-threshold point and confirms the physical sign error: its left contribution is `-3.0906335800611633e-6 - 1.0168165900800143e-7 i`, whereas its absorptive part must be positive. The single shared integrated-threshold helper is corrected to left `-i*pi`, right `+i*pi`, because runtime subtracts the helpers. Local poles are unchanged; iterated mixed terms and standalone exports inherit the correction. This correction still requires a fresh passing run.
- The export's factorized production comparison succeeds; the independent Born comparison encounters the transparent runtime identity `tree_denoms(1)`. The fixture now removes exactly this empty tree product before its structural normalization check. The three CP guard fixtures also need explicit complex CKM setup and full test initialization; their failed setup does not yet certify or invalidate the guard. Full log: `/tmp/gammaloop-phase-factorized-tests.log`.

```text
(continue as planned, but of course also make sure that int he vendored SM, no input parameter that is meant to be real is declared as something else than real-valued).
```

The SM audit also covers all input-parameter type declarations in the vendored UFO and canonical JSON. Inputs intended to be real must be declared real in both representations; genuinely complex derived quantities, including general CKM elements, retain complex declarations. This metadata audit must not silently change physical parameter values or introduce the separately pending W/Z gauge repair.

The current-file audit finds all 26 external inputs declared `real` with zero imaginary defaults in both representations: three SMINPUTS, four Wolfenstein inputs, six YUKAWA inputs, eight MASS inputs and five DECAY inputs. All 72 parameter names and nature/type declarations agree between UFO and JSON. No declaration correction is needed; derived CKM entries retain their complex type.

- The corrected 37-test run completed: 34 passed and three failed (`/tmp/gammaloop-phase-causal-fixed-tests.log`). All four signed photon/lepton ddx/ttx LO/NLO acceptances pass across their UV routes. All scalar Born multiplicities two through six, massive decay/scattering, exchange, uncut one-through-four-loop MC phases, causal bubble and both below/above-threshold triangle mirrors pass. The exact higher-pole/iterated distribution test, all three CP guards, SM vector/axial and sewing checks, generic ghost statistics and nonempty UV exports pass. Developed three-loop ghost UV checks pass in 49.2 seconds and the two-loop checks in 6.3 seconds.
- A separate exact-selector run confirms the factor-preserving computed-export normalization test passes (`/tmp/gammaloop-phase-export-factorized-tests.log`). A narrower package selection initially chose a different dependency-feature build; that unnecessary rebuild was stopped before testing, and the successful run retained the full workspace feature selection.
- Two snapshot failures have source-derived migrations applied and await rerun. Kaapo's historical `F>=1` control still reproduces 52 graphs and `-44/3`; corrected statistics on those graphs give `-5/6`, and 39 additional ghost-only graphs contribute `-15/2`, yielding `91 | -25/3 = -25/3`. GL20 changes only nine imaginary values: cut 2's threshold and event weights in four LMB channels, plus the summed imaginary result. All real values, original contributions, multiplicative factors, kinematics and event counts remain unchanged. These follow the independently certified threshold sign correction.
- The additional `slow::epem_ttxh_gl00_uv` check fails in the preexisting UV-orchestrator comparison's tensor-index canonicalizer. A read-only walk of its existing factorized expression finds all 62 top-level terms closed, with no inconsistent nested-sum spinor signatures. The used FFV1 and FFS4 model structures are unchanged. Investigation is at the canonicalizer's nested-sum contraction bookkeeping; no numerator-expansion workaround is permitted. Symbolica's existing `expand_num()` in that comparison distributes numerical coefficients only, not products of momentum/spin factors.
- The required full 166-case release scalar matrix is now running through `just test_LU_scalar_xs` and its existing 30 GB watchdog (`/tmp/gammaloop-phase-full-scalar-matrix.log`). Its result, the two snapshot reruns, final lint checks and branch commit/push remain pending.

```text
Continue as planned, but can you now give a summary of what you fixed in the end, and also does your amplitude phase targets `M` or `-i M`? And finally did you check that now the LO and NLO contribution to bother ``gamma^\star > d dx / t tx` and `e+ e- > \gamma^\star > d dx / t tx`` come out with the right sign?
```

The requested status summary confirms the amplitude target is the complete Feynman-rule graph `A_Gamma = -i M_Gamma`, operationally `-i*lambda` for the quartic scalar contact. All four signed LO/NLO acceptances pass, including positive real LO and the published signed NLO graph contributions. Final matrix and migrated-reference validation continue; no claim of a completely green broader suite is made while the separate ttH tensor canonicalizer failure remains.

- Final review confirms all physical phase owners and the three scalar propagator/coupling definitions agree. The old bubble addback-host diagnostic's whole-integrand `.expand()` is replaced by common-factor collection, preserving its exact host assertions. Its integration now trains on the imaginary component, matching the signed analytic targets. The one-loop quark and ghost UV pole references and both bubble checks join the final snapshot rerun.
- The full scalar release build completed in 22m06s under its existing watchdog and started all 166 selected cases. No build profile or memory limit was changed. The bubble diagnostic edits were applied only after compilation, outside the executing matrix's test paths.
- The optional CP guard's documented limitation is explicit: a physical phase written directly in a coupling, such as `i*g*exp(i*theta)` with real parameter declarations, is not diagnosed. Such models must disable the CP optimization and provide Hermitian partner vertices. The canonical SM instead exposes its CKM phases through named complex parameters covered by the guard.
- The required complete release scalar matrix PASSED: 166/166 cases, including slow cases, with no failures or retries. Test execution took 1945.287 seconds after the 22m06s build; the existing 30 GB watchdog recorded a peak process-tree footprint of 19,920,104,704 bytes. Log: `/tmp/gammaloop-phase-full-scalar-matrix.log`; watchdog: `target/test_LU_scalar_xs/run.ZnTSq7/watchdog.jsonl`. GL16's three-route quartic case passed in 477.454 seconds; a brief read-only stack sample located its finite cost in final UV metric/dot simplification. No profiling settings or evaluator paths were changed.
- Formatting and the three-package test-target check passed after the bubble diagnostic migration (1.69 seconds). The final seven-reference run is now active under its own 30 GB watchdog: the two derived snapshot migrations, one-loop quark/ghost UV pole references, both integrated-bubble checks and the active historical scalar-bubble MC case. Log: `/tmp/gammaloop-phase-final-reference-tests.log`.
- The seven-reference run completed in 32.436 seconds after an 81-second build: five passed, including GL20, Kaapo, the ghost pole and both bubble MC checks. Two assertions need further work. Common-factor collection alone reports two bubble addback hosts, so the diagnostic now records its factorized host differences to locate the uncancelled representation. The quark pole reaches the test's RQFT adapter with the correct external SU(3) trace, but that adapter maps the fundamental index to a closed-loop flavor count. The test now restores `T_F = 3*C_F/8` for its external projection before alignment, with an exact invariant check; the physical `-C_F` pole target and finite-slope oracle are unchanged. The two-test rerun is compiling (`/tmp/gammaloop-phase-reference-triage-tests.log`).
- The preserved commented-out threshold example now states the same left/right signs as the live helper. Ruff formatting passes for both edited UFO files; it only wrapped two Lorentz declarations. The latest formatting/type check passes (6.16 seconds). These cleanup changes do not alter the production algebra tested by the full scalar matrix.
- The two-test diagnostic confirms the quark pole and finite MUV slope PASS with the unchanged physical oracles (0.356 seconds). The bubble's second apparent host is exactly a factorized `-(A+B)+A+B` multiplied by common factors. The fixture now factors each selected route separately before subtraction; this permits exact structural cancellation without distributing any numerator. Its host-count and independent high-precision contour assertions remain unchanged, and its rerun is active (`/tmp/gammaloop-phase-bubble-factor-test.log`).
- The final bubble diagnostic PASSED (0.380 seconds), including exact single-host cancellation and the independent high-precision contour checks. All seven final reference cases now pass across their focused reruns. The temporary full-expression failure payload was reduced to host IDs after resolving the representation mismatch; no assertion or numerical oracle was weakened.
- Final `cargo fmt --all -- --check`, Ruff formatting for both edited UFO files, changed JSON/TOML/Python parsing, and `git diff --check` pass. The three-package test-target `cargo check` passes (1.13 seconds). Workspace/all-target Clippy passes without warnings after removing two redundant copies in a scalar fixture (1.15 seconds on the final cached run; `/tmp/gammaloop-phase-final-clippy-clean.log`).

## Core-commit checkpoint and then-retained limitations

This records commit `dd134b710` before the user's follow-up request. The continuation below supersedes its ttH, W/Z and optional-CP limitations.

The implemented amplitude is the complete Feynman-rule graph `A_Gamma = -i M_Gamma`. Scalar Born LU results are positive real through five forward loops, and unit-numerator amplitude MC phases are certified through four loops. All four signed photon/lepton ddx/ttx LO/NLO acceptances pass across the tested UV routes, including their published signed graph contributions. The complete 166-case release scalar matrix passes, as do the focused causal, sewing, export, complex-mass/CP, ghost-statistics and developed two-/three-loop ghost UV checks. No graph-numerator expansion or compensating generation phase is used to obtain these acceptances.

The additional `uv::slow::epem_ttxh_gl00_uv` diagnostic still fails in the `UVOrchestrator::Compare` tensor-index canonicalizer with `ExternalIndexMismatch`. Its factorized spinor-index incidence was checked, but the proposed upstream bookkeeping repair was not validated or applied. This is an explicitly retained failure, not a claim that the whole historical test suite is green. Historical unfinished scalar references and three-loop ghost placeholders remain documented as such.

The separately identified W/Z virtual-propagator/cut-projector gauge contract awaits user approval and was not changed. Complex masses are rejected; ordinary Hermitian sewing supports complex CKM, while the optional CP optimization has the documented restrictions. Regenerate integrands and saved compiled `gammaloop_state/` artifacts: their old numerical phases are not migrated. The generic statistics port also replaces the old fermion-only loop-filter names with anticommutating-loop names. The implementation is delivered on `codex/raised_energy_cff_wip_optimized_phase_fix`.

## Broader verification and pull request

```text
Summarize what you changed and confirm that `just test_gammaloop` pass. Also make sure to push a PR to merge into \`[...]\_optimization\`
```

The exact remote base is `raised_energy_cff_wip_optimized`, at the original implementation base `ac99d32d194254f0c8931873d0a32ebd7d27f82f`. Commit `dd134b710c148da0196bdb08f9a333bfc32dd93d` is pushed; [PR #103](https://github.com/alphal00p/gammaloop/pull/103) is open as a draft while the newly requested broader command is being verified. Its default selection includes eleven packages, excludes `slow::` and `failing::`, and retains the recipe's existing ARM/SymJIT exclusion. No new test exclusion is added for this work.

The pre-run audit found older full-numerator expansions in canonicalization, diagnostics and the Vakint handoff that the previously selected phase checks did not exhaustively cover. Graph/CFF diagnostic comparisons and the explicit Vakint handoff/conversion expansions are being migrated to preserve factorization. Ordinary vacuum integration and already-integrated coefficient algebra remain distinct from distributing an input graph numerator for a diagnostic. A successful exact `just test_gammaloop` run is not yet claimed.

```text
Continue
```

```text
Summarize what you changed and confirm that `just test_gammaloop` pass. Also make sure to push a PR to merge into \`[...]\_optimization\`
```

```text
Also you said in the PR:
"""
Retained limitations: `uv::slow::epem_ttxh_gl00_uv` fails in the `UVOrchestrator::Compare` tensor canonicalizer with `ExternalIndexMismatch`; the proposed upstream repair is unverified and not included. The separate W/Z virtual-propagator/cut-projector gauge change awaits approval. Complex masses are rejected; Hermitian complex CKM sewing is supported, with documented restrictions on the optional CP optimization.
"""
This needs to be fixed.
```

The continuation therefore includes fixing the ttH tensor canonicalizer failure, resolving the W/Z virtual-propagator and cut-projector gauge contract with matching UFO/JSON changes, and closing the optional CP optimization's direct-coupling phase-validation gap. The user's instruction authorizes the previously escalated W/Z correction. Whether full complex-mass-scheme support replaces the previously authorized clear rejection is being clarified separately. The complete requested suite and the additional ttH diagnostic must be rerun on the final changes; prior phase-check results do not establish those new fixes.

- The broad inventory build was intentionally stopped before test execution when the user requested the additional source fixes. It is not a full-suite result. The watchdog's recorded `signal 15` is this deliberate stop, not an observed memory-limit failure.
- A native pinned-Symbolica reproducer confirmed the ttH first-divergence boundary: `(A(i)+B(i))*(C(i)+D(i))` wrongly reports its contracted index as external; adding an independent scalar summand raises `ExternalIndexMismatch`. The repository-owned Idenso correction canonicalizes only already-existing outer summands, preserves factorization, and uses validated Spenso incidence for final dummy renaming. No dependency source/cache mutation or numerator expansion is involved.
- The combined `cargo check -p idenso -p gammalooprs -p gammaloop-api -p gammaloop-integration-tests --tests` passes after correcting two fixture pattern conversions (8.82 seconds; `/tmp/gammaloop-phase-followup-check-fixed.log`). The 17-test Idenso canonicalization run has 15 passes, including exact renamed-index, idempotence, factorization and genuine-external-slot controls. Two invalid-input controls reach the intended errors but abort while initiating a caught panic (`error 5`) on the current macOS linker setup; that runtime issue is being isolated before the full rerun. Log: `/tmp/gammaloop-phase-idenso-canonical-tests.log`.
- The uncertified initial/final CP canonicalization option and its guard metadata are removed. Ordinary Hermitian inverse-vertex sewing, same-side permutations and numerator-aware grouping remain. Named and directly inlined complex-CKM generation/model-update regressions replace the old incomplete-guard tests. Commands and persisted schema migrate without changing physical numerical oracles; graph labels/counts still require validation.
- For W/Z, the chosen consistent contract uses Feynman virtual and covariant cut numerators with complete model-declared BRST cut multiplets. This preserves the common source numerator required by raised-cut residues. Generation must include all mixed vector/Goldstone/ghost channels and the correct existing statistics/symmetry factors; event labels must use the requested physical representative so observables preserve cancellations. Matching UFO/JSON metadata and actual-model virtual-current and `H -> WW/ZZ` quartet certificates are in progress. A single imported gauge-dependent graph is a contribution to this sum, not a complete physical vector observable.
- Independent source review confirms the physical cut-state closure reaches early/late generation filters and unresolved classes, and that bare/UV/local/integrated/iterated threshold contributions share the same event-label mapping. The mapping is persisted with graph terms. Mixed explicit-unphysical and physical-vector requests are rejected when their event labels would be ambiguous; intentional graph/vertex subsets remain gauge-dependent contributions and do not synthesize missing graphs.
- The macOS caught-panic abort is independently reproduced with GCC's Darwin `-no_compact_unwind` driver behavior. The repository now selects Apple clang for macOS Rust links, retaining the required GNU runtime paths. A standalone catch-unwind control succeeds with those libraries; the complete Idenso and ttH reruns are pending.
- A fresh SM UFO import exposed an infinite square-root normalization in `ufo_model_loader`: Symbolica canonicalizes the replacement to the same expression, so a repeat-until-no-match rewrite never terminates. The source-owner fixed-point repair passes its bounded regression and the full wrapped SM import. A separate compatibility run exposed the loader's obsolete `symbolica>=0.18.0` declaration: its existing replacement/evaluation API requires Symbolica 2. The final [loader PR #1](https://github.com/alphal00p/ufo_model_loader/pull/1), commit `0bdd60af9927114a38df73d4823da07aafb1d158`, passes the bounded regression and full wrapped SM import on public Symbolica 2.0.0 and 2.2.0 (0.031 seconds per import). Both GammaLoop dependency declarations pin that revision and require Symbolica 2.0 or newer; `uv.lock` selects tested public 2.2.0. No installed-package or dependency-cache source patch is used.
- The combined five-package test-target check, including `vakint` and the API's `ufo_support` feature, passes with warnings treated as errors. The final cached check after preserving direct PySecDec indexed-sum support took 7.47 seconds (`/tmp/gammaloop-phase-pysecdec-scalar-check.log`). The focused eleven-package dev-optim build/run now includes the complete Idenso canonicalization controls, scalar-dot round trips, GL06 cancellation, direct scalar-contraction boundary, UV approximation diagnostics, W/Z sewing and gauge sums, CP-removal regressions, and the additional slow ttH UV test. Its results remain pending.

- The focused run completed: 94 tests selected, 90 passed and four failed (`/tmp/gammaloop-phase-followup-focused-tests.log`; 46m46s build, 3.363s execution, watchdog peak 9,774,948,736 bytes). All 17 Idenso canonicalization controls now pass, including caught invalid-input errors. The earlier filter mistakenly included the binary name in the ttH test name; it did not select ttH. Its required exact rerun is `binary(=uv) and test(=slow::epem_ttxh_gl00_uv)`. No ttH success is claimed from this run.
- Two gauge fixtures first disagree in numerical representation (decimal atoms versus exact rational atoms); their exact input boundary is being corrected before certifying the complete channel sums. The GL24 diagnostic rebuilt positive denominator factors from typed kinematics that intentionally discard owner tags, producing false ranks 12/2/2 instead of 6/4/4. It now retains the original tagged wrappers, with all original physical/rank assertions unchanged. The raised-LU oracle needs explicit finite-component contraction for preserved indexed sums; it now uses the existing tensor-network executor before its independently supplied exact coordinates and Arb jets. No whole numerator is distributed. These fixes await the next run.
- Cross-section mass validation now also rejects used nonstandard UFO propagator denominators, including hidden finite-width terms whose particle mass is declared real. It validates the current model at generation and warm-up, while retaining ordinary width metadata, unused custom propagators and amputated external lines. Only denominator expressions are normalized. Three focused regressions are queued.
- Full-suite Python subprocess tests require a fresh extension; the host default installation predates this branch. An isolated environment at `/tmp/gammaloop-phase-plan/python-api-env` now contains the locked public Symbolica 2.2.0 and tested loader revision. Its branch extension build and the exact broad suite remain pending; no global Python installation was changed.

- The nine-case rerun completed after a 14m13s build: four passed and five failed (35.121s execution, peak 8,852,958,208 bytes; `/tmp/gammaloop-phase-four-fixes-all-targets.log`). The complete raised-LU value/derivative and production-boundary comparisons pass with finite component contraction, as does the expanded open-index rejection coverage. Standard denominator normalization and current-model validation pass. The amputated custom-denominator fixture used an uncontracted fermion external leg; it now uses the scalar Higgs and its width expression to reach the intended guard without disabling tensor validation.
- The gauge rerun stopped in the attempted zero-error rationalization API, which requires a nonempty approximation interval. That attempted fixture conversion is removed. Comparing the original source snapshot with its reported failure line establishes the actual first precision boundary: the raw Goldstone propagator was exact, while the component library produced floating-point spatial Minkowski metric entries. Both fixtures now override that library's existing generic Minkowski key with its existing exact parametric metric constructor. A bounded native test using the pinned libraries proves exact diagonal entries, including a factorized complex coupling; no model coefficient, graph numerator or computed result is rationalized. The scalar-result extraction is simply a shorter exact API, not a fix to a conversion in result_tensor.
- GL24 now passes its original owner/rank assertions and reaches a mass-only coefficient form difference. The diagnostic normalizes numeric coefficients with expand_num(), whose implementation leaves all momentum factors and products of sums intact. A pinned native probe certifies the exact equality and preserves a separate product-of-sums control. No general expand() or polynomial conversion is used.
- The first correctly selected ttH execution no longer stopped during generation, but its old label selected a different graph after removing CP canonicalization: current GL00 has no UV subgraph, so the mandatory scale-dependence assertions correctly failed. The original preserved failure numerator encodes a top self-energy bubble. An exact directed graph isomorphism, including particles, interaction IDs, fermion arrows and external assignments, uniquely maps that original topology to current GL14. The existing GL00-named fixture now selects GL14; UVOrchestrator::Compare and every UV-profile, scale-dependence and integrated-invariance assertion remain unchanged. Mapping evidence: `/tmp/gammaloop-phase-plan/TTH_GL00_TO_GL14_MAPPING.md`. The mapped regression remains to be rerun.
- Three additional selected Idenso gluon diagnostics now preserve their numerator products, including settings that otherwise force sum contraction expansion despite a false const argument. Index-closure assertions validate every sum branch and require the independently known exact external slots on both original and simplified expressions. Bound internal names may remain textually inside a factorized product; their literal disappearance is not required. Historical timing comments and useful nonexpanding comparisons are retained. Contributor guidance now explicitly applies the no-expansion requirement to diagnostic copies.

- The mapped follow-up run PASSED all 18 cases (`/tmp/gammaloop-phase-mapped-tth-tests.log`): the intended ttH self-energy case, all sewing and gauge checks, GL24, custom-denominator rejection and all selected gluon diagnostics. Build: 5m58s; execution: 188.317s; watchdog peak: 6,551,865,024 bytes. The mapped ttH case passes UVOrchestrator::Compare, its UV profile, required mUV/muR inspect dependence and integrated scale invariance. The exact generated H→WW/ZZ covariant sums equal the independent three-physical-polarization norms at both benchmark energies and boosts. Named/inline complex-current and charged-scalar sewing remain correct.
- A fresh Python extension is now building in the isolated locked environment with the required compile-time NO_SYMBOLICA_OEM_LICENSE=1 setting. It will be selected explicitly for the exact full suite; this switch is scoped to the extension build. No full-suite success is claimed yet.

- The fresh isolated Python extension built and installed successfully (25m23s; `/tmp/gammaloop-phase-fresh-python-api-build.log`). Final read-only review then found two narrow follow-ups: direct PySecDec parameter discovery must inspect the reduced scalar numerator so canceled coefficients are not requested, and an imported covariant partner collection needs an explicit physical `--process-spec` because raw graph states cannot establish observable intent. The existing scalar-reduction and exact Higgs fixtures now cover these boundaries; no numerical oracle is changed or numerator expanded. The imported Higgs fixture uses the same six/four graphs through `Process::from_graph_list`.
- Formatting and the five-package test-target check with UFO support pass after those edits (13.60s; `/tmp/gammaloop-phase-final-review-dev-check.log`). An accidentally selected dev-optim check was stopped before completion to retain the established dev check profile; it is not a test failure. The incremental final extension rebuild is running before public-API model parity and the exact broad suite.

- The final incremental Python extension rebuild passes (1m04s; `/tmp/gammaloop-phase-final-python-api-build.log`). Public `api.run`/`get_model` import parity passes for all 43 particles and propagators, 72 parameters, 108 couplings, 153 vertices and 22 Lorentz structures, including all twelve W/Z-sector propagators, Goldstone/mass metadata, cut multiplets and 26 real input declarations. Signed real/imaginary coefficients agree at a common real input point with nonzero CKM phase. Six local Lorentz expressions differ only by an external versus distributed numerical 1/2; scalar-only `expand_num()` certifies equality with tensor functions and products opaque. No graph numerator is constructed or expanded. This uses the loader's existing default restriction without simplification, then explicitly resets all input parameters; it does not claim an unrestricted-loader test.
- Ruff was applied to the new propagator module, UFO initialization and added object-library class; Python AST equality is verified. Existing generated formatting outside the edited class and legacy process-test command strings is retained. The post-format public API parity rerun passes in 0.22s (`/tmp/gammaloop-phase-public-api-sm-parity-final.log`). Rust formatting, changed-file Python/JSON/TOML parsing and `git diff --check` pass.
- The exact `just test_gammaloop` command is now running under the 30 GB watchdog with the recipe's five test threads and explicit fresh Python interpreter (`/tmp/gammaloop-phase-exact-test-gammaloop.log`; `target/phase_fix_first_checks/exact-test-gammaloop-watchdog.jsonl`). No arguments or new exclusions were added. Its result remains pending.

- The first exact broad run completed after a 7m11s build: 2,050 tests executed across 46 binaries, 2,031 passed, 19 failed, and 296 remained outside the unchanged default selection. Execution took 211.825s; watchdog peak 6,097,816,768 bytes. All four signed photon/lepton ddx/ttx LO/NLO acceptances and the five fresh Python subprocess tests pass in this run. This is not a green full-suite result.
- The failures reduce to twelve cases at the same Vakint momentum-system boundary, three exact identity diagnostics needing numeric-sign normalization while retaining numerator products, and four graph-label/grouping fixtures. The shared Vakint failure precedes integration: fusion contracts the denominators outside a nested denominator sum before its branches are separated, leaving a fictitious bridge carrying nonzero loop momentum. The former whole-expression expansion accidentally separated these branches. The repair is denominator-only collection before fusion, with numerator coefficients opaque; neither disabling the momentum solver nor expanding numerators is acceptable.
- Exact graph recovery from preserved prior artifacts maps old ttH GL18 uniquely to current GL22 with identical node/edge IDs, and old GL20 uniquely to current GL38. The latter requires an explicit edge/channel mapping before its event snapshot can be migrated. A subsequent read-only export from the original CLI confirms that all physical edge directions and default loop coordinates are preserved; an apparent gluon reversal in the preliminary raw-tensor reconstruction was a node-label ambiguity. Numerical physical references are retained. The CP grouping snapshot differs only by provenance IDs; its graph count, coefficient multiset and signed sum are unchanged. The full-SM grouping change still needs attribution to the corrected model rules before updating its snapshot.

- The denominator-incidence repair now collects only complete denominator atoms before graph contraction, retaining each numerator as a factorized coefficient. The new actual-to-Vakint regression checks the two independently expected one-loop topologies, raised powers and unchanged `(a+b)*(c+d)` coefficient, including powers above the signed-eight-bit range. The collection uses the existing propagator power range (`i32`). The native signed-denominator collection proof passed; the rebuilt regression and complete suite remain pending.
- The three identity diagnostics now use factor-preserving exact identities: numerical sign normalization for affine momenta, evenness only of the squared GL04 contraction with its linear occurrences preserved, and common-energy denominator clearing on scalar CFF kernels before inserting any mapped numerator. Native exact controls pass, including a retained product-of-sums spectator. No whole-numerator rational-polynomial conversion or numerical sampling substitutes for an exact identity.
- The ttH selector/channel migrations are applied with all numerical targets, sample coordinates and cut/channel order retained. A model-clone 2×2 control attributes the full-SM grouping change solely to the approved FFS1/FFS3 projector corrections: all 47 raw directed graphs remain identical, while eight formerly spurious Lorentz zeros become the required nonzero mass-insertion traces. Restoring only the old W/Z propagators leaves the result unchanged. Full-SM and CP grouping references are updated from these controls; their signed sums are independently checked.

- Final review confirms denominator collection matches the existing signed-32-bit propagator power range and keeps numerator coefficients opaque. The five-package test-target check with UFO support passes in 8.49s, followed by a fresh isolated Python API rebuild in 6m13s (`/tmp/gammaloop-phase-broad-repair-i32-python-build.log`). All edited Python/JSON/TOML files parse. Independent read-only reviews found no remaining whole-numerator expansion or lost useful comments in the final changes.
- The repaired exact `just test_gammaloop` rerun is active without added arguments or exclusions, using the freshly built Python extension and the recipe's five test threads (`/tmp/gammaloop-phase-exact-test-gammaloop-repaired.log`; `target/phase_fix_first_checks/exact-test-gammaloop-repaired-watchdog.jsonl`). The earlier Linux CI run on the core commit has only the same explained full-SM grouping failure (163/164 integration tests pass); it does not establish the final follow-up result.

```text
Yes, continue as planned and obviously no need to concern yourself with or support at all backward compatibility!
```

The continuation uses the new schema and physical conventions directly. Old states and integrands must be regenerated; no compatibility layer is introduced.

```text
Continue as planned but what was your rational for the new schema?
```

The data-model changes have three owners: remove the uncertified CP option and frozen assumptions; declare covariant cut multiplets in the model; persist the physical representatives used for cut-event labels across bare, UV and threshold evaluation. These changes make gauge-state closure and observable semantics explicit. They do not require a separate file-format redesign.

```text
But I want to keep:
"""
**Remove&#x20;****`symmetrize_left_right_states`****.** That optimization exchanges the forward graph’s sides using CP, which is not generally valid with complex couplings. Its configuration and saved assumptions are removed.
"""
because the user may know it's valid for his own theory, but you can emit a warning at runtime when the user uses it if you want.
```

This instruction supersedes the earlier removal of `symmetrize_left_right_states`. Retain its CLI/API/serialized configuration and original CP/crossing behavior, disabled by default. Enabling CP-based forward-side canonicalization is a user assertion about its validity for the theory, process and current coupling point. Generation and cross-section runtime warm-up warn about that assertion; complex inputs alone are not rejected by a CP-parameter guard. Preserve the separate real-pole/CMS checks. Retain the flag in generated integrands so the warning also applies after model updates.

- The repaired exact-suite build was intentionally stopped before test execution when this instruction arrived (watchdog exit 125 after SIGTERM). It is not a test result or memory-limit failure. No build or native probe remains active while the option is restored.
- Original example opt-ins and the four signed QCD acceptance opt-ins are restored, with their numerical targets unchanged. The first exact run already passed those four acceptances with CP disabled. The mapped ttH UV/event fixtures now select CP=false explicitly to preserve their certified directed-topology/edge mappings. Existing CP-specific grouping fixtures and API/runtime checks are being restored; the final exact suite will run after this update.

- The CP option is restored in production, configuration, runtime metadata and original example/test commands. Focused source coverage now includes CLI/serialization defaults and opt-in, actual mirrored scalar topology canonicalization, and named/inline complex-CKM generation and model-update warm-up with the option both disabled and enabled. Complex-input acceptance tests certify the user-controlled option contract, not CP validity of a violating coupling point.
- Formatting and the restored five-package test-target check with UFO support pass (8.95s; final cached check after the one Clippy cleanup: 5.92s). Workspace/all-target Clippy with UFO support passes without warnings (20.43s; `/tmp/gammaloop-phase-cp-opt-in-workspace-clippy-clean.log`). The preliminary check caught three explicit fixture initializers missing the restored bool; their original choices were restored. The serialization regression retains the required-field policy, with no compatibility default.
- A fresh Python API rebuild is active for the restored CP option (`/tmp/gammaloop-phase-cp-opt-in-python-build.log`). The full-SM CP-enabled grouping reference will be recertified through this public API before the final exact suite. No full-suite success is claimed from the deliberately interrupted build.

- The restored-CP Python API build completed in 6m49s. Fresh full-SM generation with CP enabled gives 40 groups and signed coefficient sum −47. Independent model-clone controls restore either the old Goldstone Yukawa structures, the old W/Z propagators, or both: only restoring the old Yukawa structures returns 36 groups and eight Lorentz zeros. All 47 raw directed graph structures remain identical in all four cases. The eight restored mass-insertion graphs form four CP pairs, explaining the updated snapshot without changing a physical target (`/tmp/gammaloop-phase-plan/full-sm-cp-opt-in-snapshot-6cqz2lfl/report.json`, `/tmp/gammaloop-phase-plan/full-sm-cp-opt-in-controls-m7aldx_t/report.json`).
- The mirrored scalar canonicalization control uses a cubic box on one side and a quartic Born vertex on the other. Their different external-neighbor valences exclude a side-preserving graph automorphism, so the disabled-option case is a valid negative control. Formatting and the final five-package test-target check pass (5.55s; `/tmp/gammaloop-phase-cp-final-test-check.log`). The exact `just test_gammaloop` rerun, with no added arguments or exclusions, is active using the final fresh Python extension (`/tmp/gammaloop-phase-exact-test-gammaloop-cp-retained.log`).

```text
Continue, and also:

Launch a subagent to study the slowdown occurring during generation and comparison with numerators of 50k+ diagrams. It's not due to the evaluation time of a single numerator but it's during the lookup with all previous accepted entries (of course it's unavoidable to have this growth in lookup time but maybe an agent can find a non-invasive easy improvement/optimization in feyngen for this).
```

A separate agent is tracing the accepted-numerator lookup in Feyngen and looking for a small optimization that preserves grouping semantics. This investigation proceeds independently of the active phase-validation suite; no generation-source changes are made during its build.

The generation slowdown investigation traced the accepted-numerator lookup to `diagram_generator.rs` around the `pooled_bare_graphs` occupied branch. A single global mutex protects the topology map, and its `find_map` linearly invokes symbolic/sample numerator comparison for every accepted class in a topology bucket. A bucket with k classes therefore has O(k) comparisons per candidate, while the mutex also serializes unrelated worker lookups. This establishes quadratic worst-case comparison growth from source, not a measured 50,000-diagram speedup. The agent is evaluating a conservative index over already computed sample results, with the existing comparator and a full-scan fallback retained. Raw polynomial-variable hashes or naive specialization are not certified keys: nonpolynomial variable slots may encode algebraic relationships that the existing comparison restores. The study does not change generation source during phase validation.

- The exact restored-CP suite completed: 2,053 tests ran, 2,050 passed and three failed, with 296 existing skips (22m12s build, 253.533s execution, watchdog peak 10,390,829,248 bytes; `/tmp/gammaloop-phase-exact-test-gammaloop-cp-retained.log`). All four signed QCD acceptances and the restored CP regressions passed. The three failures are isolated below; this is not yet a successful full-suite result.
- The GL20/GL38 event snapshot differs only in CFF discovery order: within every channel, cut 0 and cut 2 exchange their event enumeration positions. Comparing all sixteen events by physical cut/channel IDs with JSON numbers retained as literal strings proves every momentum, complex weight, counterterm, normalization and top-level result unchanged. The snapshot migration reorders whole records and their event indices; no numerical reference changes.
- The intermittent differential JSON failure is a test-directory race: two nextest subprocesses used the same fixture root, and one cleaned it between its sibling creating a workspace and writing an observable. Timestamp and surviving-directory evidence identify this boundary directly. Six JSON/HWU fixtures now use their own test names with the existing setup helper and derive workspaces from that state folder; assertions and runtime settings are unchanged.
- The affine identity's separate normalized expressions can be identical while their unevaluated subtraction remains `a+b-(a+b)` under factor collection. A bounded native scalar control demonstrates this case. The test now compares the two numerically normalized factorized expressions directly with exact equality, preserving its physical mapping oracle and improving failure diagnostics. No graph numerator is distributed.

- The bounded Feyngen study is complete at `/tmp/gammaloop-phase-plan/FEYNGEN_LOOKUP_PERFORMANCE.md`. The smallest recommended changes are early rejection once sampled ratios fail or disagree, and separate locks per topology. The current comparator evaluates all five default samples for every candidate even after acceptance is impossible. A more ambitious candidate index needs conservative algebra-preserving keys and explicit fallback coverage; no 50,000-diagram speedup is claimed and no generation optimization is applied in this phase checkpoint.

- The final focused run passes all six selected checks: the three repaired tests, the two sibling differential-output fixtures and `uv::slow::epem_ttxh_gl00_uv`, including its `UVOrchestrator::Compare`, UV profile, scale-dependence and integrated-invariance assertions. Build 4m22s, execution 218.517s, watchdog peak 5,284,904,640 bytes (`/tmp/gammaloop-phase-final-three-and-slow-tth.log`). The exact full-suite command is now rerunning on these unchanged binaries, without additional arguments or exclusions (`/tmp/gammaloop-phase-exact-test-gammaloop-final.log`).

- Final exact `just test_gammaloop` **passes**: 2,053 tests run, 2,053 passed, zero failed, 296 existing skips; nextest run `a1bcc3fc-06ae-4f10-809e-d07df8b48a10`. Cached build 0.60s, execution 274.024s, watchdog peak 1,607,286,784 bytes (`/tmp/gammaloop-phase-exact-test-gammaloop-final.log`, `target/phase_fix_first_checks/exact-test-gammaloop-final-watchdog.jsonl`). No arguments or exclusions were added to the requested command. All four signed photon-current/lepton-annihilation ddx/ttx LO/NLO acceptances pass on the final implementation with their original CP opt-ins restored. The separate final slow ttH UV test also passes, as recorded above.

- Final formatting, whitespace checks and `cargo clippy --workspace --all-targets --locked --features gammaloop-api/ufo_support` pass with warnings treated as errors (21.75s; `/tmp/gammaloop-phase-final-workspace-clippy.log`). This is the verified phase checkpoint before the separately requested Feyngen lookup optimizations.

```text
Implement all not-too-invasive recommendations for the improvement in feyngen.
```

Phase checkpoint `5c8b46f39` is committed and pushed. Implement early sample rejection, cached zero-mask filtering after the symbolic comparison route, and per-topology locks. Preserve the existing accepted-class ordering within each bucket and all exact grouping/ratio semantics. The broader algebraic-signature index remains outside this small change because nonpolynomial symbolic slots need additional specialization and fallback handling. Verify comparator edge cases and normalized graph memberships with one and several generation workers, then rerun the exact full suite before delivery.

- Implemented the three bounded lookup changes and removed redundant graph clones. Scalar/sign sample comparisons return immediately on a missing or inconsistent ratio; the cached zero-mask guard runs after the canonical symbolic route; the global map lock ends at bucket retrieval and each topology retains atomic first-match search/insertion. The existing scalar-ratio validation and representative normalization remain. New exact coefficient-cache tests cover nineteen cases plus canonical precedence/open-index rejection, and a small SM regression compares complete grouping provenance, directed DOTs and signed factors across both grouping modes with one/four workers.
- Formatting, whitespace checks and the five-package test-target check pass (11.39s; `/tmp/gammaloop-phase-feyngen-lookup-check.log`). A pre-optimization public-API baseline independently certifies the four small SM cases at checkpoint `5c8b46f39`, each with two graphs and signed factor sum −3 (`/tmp/gammaloop-phase-plan/feyngen-lookup-baseline-kv5uh7vf/report.json`). The fresh optimized Python API build is active; no post-optimization full-suite success is claimed yet.

```text
Stop yourself temporarily as soon as you are at a good stopping point, getting out of quota.
```

Paused at the user's request after the active Python API build completed successfully (6m12s; `/tmp/gammaloop-phase-feyngen-lookup-python-build.log`, watchdog peak 6,217,828,864 bytes). No build, test or native probe remains running. The verified phase checkpoint `5c8b46f39c485bf2f7162319df6c2de27f8ac031` is pushed to PR #103; its exact full suite passed 2,053 tests and its slow ttH check and final Clippy passed. The subsequent Feyngen optimizations and regression tests remain uncommitted, with formatting and cargo check passing. Their fresh Python extension is installed, but optimized old/new API comparison, new regression execution, exact full-suite rerun, final Clippy, commit/push and PR readiness remain pending. Resume details: `/tmp/gammaloop-phase-plan/PAUSE_CHECKPOINT.md`.

The locking optimization preserves a serialized search/insertion order within each topology, equivalent to an allowed scheduling order of the previous global mutex. The existing comparator is not universally transitive for noncanonical zero samples; no new claim of arrival-independent grouping is made. The fixed physical fixture tests compare normalized outputs, while redesigning that pre-existing comparison behavior remains outside this optimization.


```text
Run multiple agents (like 4) to cross-audit the code and the changes in depth
```

Four independent source audits completed on 2026-09-07 (phase/CFF/LU, model/sewing, tensor/UV, Feyngen), with root cross-checks. Reports are under `/tmp/gammaloop-phase-plan/CROSS_AUDIT_*.md`. They found three Vakint follow-ups: the pure scalar fallback still distributed the complete numerator inside FORM, odd open loop tensors could be projected to zero before scalar validation, and partial/full numerical consumers were not migrated when dot conversion stopped contracting across sums. A fourth finding concerns the missing covariant-sector diagnostic for anticommutating-loop and related diagram filters. No additional phase-accounting or introduced Feyngen concurrency/comparison defect was established. Source inspection is not runtime reproduction; the reports distinguish test gaps and inherited limitations.

Implementation resumed under the active goal after the audit. Replace the pure scalar FORM fallback with a factor-preserving contraction/validation boundary and migrate all numerical consumers. Preserve deliberate diagram filtering while extending the existing gauge-subset warning. Recheck the explicit integrated vacuum-projection boundary separately: removing a local expansion call must not merely relocate whole-numerator distribution downstream. Then run the final required pipeline and update PR #103; checkpoint test results do not certify the changed tree.

- The prebuilt optimized Python extension passes exact public-API comparison with checkpoint `5c8b46f39` in all four small SM cases (sign/scalar grouping, one/four workers), including complete sorted DOT blocks, provenance and signed factors: `/tmp/gammaloop-phase-plan/feyngen-lookup-optimized-e4htikx0/report.json`, `/tmp/gammaloop-phase-feyngen-lookup-public-api-optimized.log`.
- The broader full-SM CP-enabled generation also reproduces all 40 complete DOT blocks from the pre-optimization checkpoint, with signed factor sum -47: `/tmp/gammaloop-phase-plan/full-sm-cp-opt-in-snapshot-90tyr0rf/graphs.dot`, `/tmp/gammaloop-phase-feyngen-full-sm-optimized.log`. These checks certify the existing compiled Feyngen optimization, not the subsequent Vakint/diagnostic source edits. Both ran under the 30 GB watchdog and completed successfully.

- The pure Vakint contraction fix is implemented in its private `lorentz` module. The converter validates sum/product index domains, preserves independent scalar factors and disconnected contractions, handles symbolic metric traces, and does not use FORM. Its Result-returning API and the partial numerical evaluator were migrated across Vakint sources, tests, examples and documentation. The first compile check caught one metric-constructor typo; after that correction, the five-package test-target check passes (11.93s; `/tmp/gammaloop-phase-audit-contraction-check-2.log`). The coordinated focused nextest pipeline is still active at this record; no test success is claimed yet.
- The covariant-filter diagnostic now covers restrictive diagram/topology filters. The existing generated Higgs gauge oracle checks the actual WW six-to-four graph reduction under the anticommutating-loop filter and the warning predicate; complete WW/ZZ and explicit Goldstone diagnostic requests remain quiet. Low-level automatic zero-flow cleanup was excluded from this warning because no missing Born partner was established for it.
- The integrated vacuum path needs a further change: universal tensor kernels and loop-dot monomials may enter FORM, while retained graph coefficients remain factorized outside it. A source-only prototype and owner patch are in `/tmp/gammaloop-phase-plan/vacuum-projection-prototype/`; they are not yet live or verified. They account for projector epsilon poles and combine numerical terms into one estimator. The existing PySecDec adapter also discards phases by taking norms of complex parameters (and signs of real-valued complex external components); `/tmp/gammaloop-phase-plan/PYSECDEC_PARAMETERS.patch` prepares its correction, again not yet applied or tested. These remaining items prevent a completion claim.

- The current four-reviewer cross-audit is recorded in `/tmp/gammaloop-phase-plan/CROSS_AUDIT_CURRENT_SUMMARY.md`. It confirms seven remaining issues: the inherited PySecDec norm conversion and full-numerator FORM paths; temporary machine-local test instrumentation; exponent validation; malformed Lorentz domains entering explicit vacuum projection; two new test initialization failures; and weakened contraction-progress assertions in two Idenso tests. Valid odd-rank vacuum integrals reducing to zero are intentional and are not a defect. No additional core CFF/LU phase or Feyngen synchronization defect was established.
- The focused contraction run completed with 85 passes and two failures, both symbol-registry initialization panics before the new tests' assertions (`/tmp/gammaloop-phase-audit-contraction-targeted.log`). The 50,000-record lookup diagnostic passed: one same-topology last-match scan took 82.237 ms in sign mode and 322.749 ms in scalar mode. These are one-sided lookup timings, not an end-to-end generation speedup. Its temporary absolute include is removed from the repository test source after the completed run. The historical 2,053-test full-suite pass does not certify this changed worktree.
- Implementation resumes with separate ownership for numerical-parameter transport, factor-preserving vacuum projection, and Lorentz/test validation. Apply the prepared proposals only after reviewing their current owners; verify actual backend calls, not just source applicability. Keep all graph coefficients outside FORM, preserve a single numerical estimator for combined scalar sectors, and validate Lorentz domains before odd-rank projection. Follow with affected tests, the exact full suite, the separate slow ttH UV check, fresh model parity and final local gates before commit/push.

```text
And once you've finished the audits and included all the changes you deemed necessary (make sure you do the quick fast pySecDec tests in vakint since you changed the pysecdec eval method and those corresponding tests are disabled by default (keep them disabled but run them once at least here, only those fast enough). Then make sure you commit+push a merge-ready PR (with all local gates and test ready) and then babysit the PR at a rate of one check every 30min.
```

Run the fast existing PySecDec tests explicitly, preserving their default disabled state, together with a fresh-output numerical regression of complex coefficients and negative real values. After all required local gates pass, commit and push the complete change, mark PR #103 ready against `raised_energy_cff_wip_optimized`, and start a thread follow-up every thirty minutes to monitor its checks and actionable review changes. Remain quiet while the remote state is unchanged; report failures or required user decisions and completion.

- The audit implementation now preserves complex PySecDec parameters and signed real momentum components, keeps retained graph coefficients outside FORM using universal tensor/scalar kernels, and accounts for coefficient epsilon poles before requesting backend orders. Lorentz validation includes exponents and occurs before odd-rank projection. The optional Vakint Python interface compiles with the new complex parameter map. The five-package test-target check with UFO and community-module features passes (`/tmp/gammaloop-phase-audit-fixes-check-4.log`).
- The first current-tree focused run completed on 2026-09-08: 106 tests ran, 100 passed and six failed (25m55s build, 3.362s execution, watchdog peak 11,023,350,016 bytes; `/tmp/gammaloop-phase-audit-fixes-targeted.log`). The eleven projector module tests, parameter transport and exponent/domain regressions, Feyngen comparison/grouping tests and generated covariant-cut oracle passed. PySecDec numerical execution was disabled in this run. The six failures are fixture/oracle boundary issues: two component oracles left implicit scalar products unassigned, the new analytic fixture used a literal mass instead of a recognized mass symbol, and three decorated-index fixtures mixed scalar and rank-two summands. Fix those boundaries while preserving numerical targets, then rerun before broader acceptance.
- A final source review confirms the epsilon order handoff for supported vacuum integrals, including poles from high-rank projectors. An initially suspected nonvacuum double-pole serialization issue was withdrawn as a supported-path failure: current topology matching rejects external denominator momenta. No nonvacuum support is added in this review. Source inspection does not replace the pending enabled PySecDec run or final local gates.

- The corrected current-tree focused run passes all 106 selected tests (21.43s rebuild, 5.693s execution; `/tmp/gammaloop-phase-audit-fixtures-targeted.log`). The fixture-only edits preserve every numerical reference. The Idenso oracle now supplies the implicit Minkowski representation on bare vectors so the existing component materializer evaluates compact scalar products; exact before/after assignments and nonzero/progress assertions remain enabled.
- Explicitly enabled fast PySecDec acceptance passes all four selected one-loop tests, generating fresh packages: simple tadpole, cross-product numerator, additional numerator symbols, and the new complex/signed-coefficient regression. The latter checks `1+2i`, `-3`, and `(1+2i)/epsilon`, including the finite coefficient requiring an extra integral order. Runtime 61.498s, watchdog peak 1,076,841,152 bytes (`/tmp/gammaloop-phase-fast-pysecdec-explicit.log`, nextest run `8e5a5b27-8acb-44f5-9014-2d2943d2ec59`). `RUN_PYSECDEC_TESTS=1` was set only for this filtered invocation; all default disabling remains unchanged. The full suite and final broad gates are still pending on this tree.

- Workspace/all-target Clippy with UFO support and warnings treated as errors passes after removing unnecessary clone/borrow operations and retaining the existing local convention for a private tuple return (`/tmp/gammaloop-phase-review-final-clippy-2.log`, 29.66s). The five-package test-target check, including the optional Vakint Python module, also passes after this cleanup (9.97s). No physics or numerical algorithm changed. Refresh the Python extension after these source-level cleanups before final native API and exact-suite acceptance.

- The exact current-tree `just test_gammaloop` run completed with 2,052 passes and 23 failures out of 2,075 tests, with 296 existing skips (179.697s execution; `/tmp/gammaloop-phase-review-final-test-gammaloop.log`). This is not a green final suite. The four signed NLO acceptances fail during integrated UV generation before their numerical assertions. Fresh native UFO/JSON parity on the final rebuilt Python extension passed separately (`/tmp/gammaloop-phase-plan/public-api-sm-parity-dknr4ob9/report.json`).
- Pipeline triage isolates two representation boundaries. Retained external-dot coefficients were restored after the analytic backend had applied `use_dot_product_notation=false`, so bare Vakint dots leaked into GammaLoop's evaluator; apply output notation after coefficient restoration at the owning boundary. Separately, the new Lorentz validator knows momentum/metric slots but not the gamma tensors restored by GammaLoop's pre-Vakint conversion. The canonical ddx trace contains contracted `gamma(mu)*k(mu)` terms with different dummy names, which the validator incorrectly classifies as different free-index sets. Carry explicit opaque tensor-slot metadata from the existing GammaLoop tensor parser into Vakint; keep the tensor bodies and graph coefficients factorized outside FORM. Preserve all signed acceptance targets and rerun the affected paths before the exact suite. Trace: `/tmp/gammaloop-phase-plan/integrated-domain-trace/ddx-canonical.json`.

```text
Continue as planned, and get this PR in a merge-ready stat with green CI (once all local gates are green, babysit the PR at a 30min interval)
```

The final acceptance remains the complete local gates followed by a pushed, ready PR and green CI. Start thirty-minute monitoring only after the local gates pass; keep the monitor quiet while the remote state is unchanged.

- The explicit tensor-slot handoff and analytic output-notation repair pass the affected physics paths: all four signed photon/lepton ddx/ttx LO/NLO acceptances, scalar UV cases, ghost renormalization and raised-cut selection. The focused run completed with 154 passes and two test-oracle failures out of 156 (204.967s execution, 15m52s rebuild; `/tmp/gammaloop-phase-review-tensor-boundary-targeted.log`, run `8afd977a-0621-4d55-b21a-1f6b477c8df2`). One oracle hard-coded a former UV dummy name; the new closed-gamma oracle incorrectly required a fresh index when an original dummy is equally valid. Correct those assumptions using existing contraction normalization and exact alpha-renaming, preserving all numerical targets and spectator factors, then rerun. Compilation checking and all-target Clippy passed before this run.

- Final-tree enabled fast PySecDec passes 4/4 (52.984s, run `00e88e52-e3dd-460d-bc60-5c22bc9e313b`, `/tmp/gammaloop-phase-review-fast-pysecdec-final.log`). Fresh Python rebuild and SM UFO/JSON parity pass (`/tmp/gammaloop-phase-plan/public-api-sm-parity-yxnwcs62/report.json`). Exact `just test_gammaloop` passes all 2,079 tests, with 296 existing skips and no added arguments (258.858s; `/tmp/gammaloop-phase-review-exact-suite-final.log`). All four signed photon/lepton ddx/ttx LO/NLO acceptances are included and pass at their existing MC statistics.
- The separate slow ttH UV test fails during `UVOrchestrator::Compare` with `ContractedMoreThanOnce` in Idenso tensor canonicalization, so the PR is not ready. Both raw Vakint outputs are identical and the two integrated expressions are exactly equal after renaming their complete dummy slots (`/tmp/gammaloop-phase-plan/tth-integrated-return-alpha-comparison.json`). This excludes a phase or integrated-algebra discrepancy. Investigate generic capture-free dummy allocation in Spenso parsing and exact canonical reconstruction; do not weaken comparison or special-case ttH. Trace: `/tmp/gammaloop-phase-review-tth-index-trace.log`.

```text
OK, and make sure the fix to `Idenso reports a Lorentz dummy contracted more than once` is generic. We really need a robust and generic gammaloop.
```

- The generic capture regression reproduces the original duplicate-index error with one explicit dummy and a compact scalar product, independent of ttH. Spenso now reserves serialized input-index names before materializing compact dots/chains/traces; Idenso's final relabeling also reserves genuine external slots. Three focused regressions pass, including exact alpha-equivalence, idempotence, external-slot preservation and retained scalar factors. The UV comparator no longer calls `expand_num`, which distributed numeric coefficients over sums. Check and all-target Clippy pass (`/tmp/gammaloop-phase-review-generic-indices-check.log`, `/tmp/gammaloop-phase-review-generic-indices-clippy.log`).
- Slow ttH passes the original capture boundary but fails later in Symbolica's nested-sum reconstruction (`tensors.rs:469`), so the previous 2,079-test pass does not certify this changed tree. Two independent source audits identify a completed contraction counter being retained as open. The minimal factorized expression `s*((A(i)+B(i))*(C(i)+D(i))+(E(i)+F(i))*(G(i)+H(i)))` reproduces the identical panic in a native regression (`/tmp/gammaloop-phase-review-nested-sum-baseline.log`). Repair the counter at its Symbolica owner, based on the exact locked revision; do not add spare index names, expand products, or bypass validation. Current upstream `dev` still has this defect and has diverged from the locked revision, so an unrelated upstream update is inappropriate.
- The minimal Symbolica repair is committed and pushed to the existing `alphal00p/symbolica` fork at `4d0a833eb8e059d1f95bdae5abed2559830b235f`, directly atop the previously locked revision. Nested sum branches merge occurrence counts by maximum, preserving completed contractions and the existing open-index checks. All nine tensor canonicalizer tests pass, including nested closed sums, squared sums, shared external slots, exact alpha-equivalence and factor preservation (run `1c0ed517-8972-4ef7-905f-f19b5456b6ef`; 5m31s build, 0.026s execution, watchdog peak 9,335,084,800 bytes; `/tmp/gammaloop-phase-review-symbolica-counters-tests.log`).
- GammaLoop pins all three related crates to that immutable repair. `cargo hakari generate` and `cargo hakari verify` pass using official cargo-hakari 0.9.38; the generated feature-unification manifest retains its feature lists and the lockfile changes only the three source entries. Idenso now validates the complete factorized expression through Symbolica, removing the earlier top-level-summand workaround. Final relabeling uses validated dummy metadata, reserves genuine external slots and ignores names belonging to canceled terms. New dual-slot and mixed-representation cancellation regressions complement the generic capture and nested-sum cases. Formatting, test-target checking and all-target Clippy pass before the pending GammaLoop native acceptance (`/tmp/gammaloop-phase-review-canonical-unified-check.log`, `/tmp/gammaloop-phase-review-canonical-unified-clippy.log`).
- The rebuilt GammaLoop focused acceptance passes all eight selected tests, including the complete slow ttH UV test and genuine malformed-index rejection checks (run `ffe36421-65af-4998-93b5-8b6e0358ae28`; 40m15s build, 344.513s execution, watchdog peak 9,659,919,040 bytes; `/tmp/gammaloop-phase-review-canonical-owner-targeted.log`). Nextest noted one passed test with a delayed output-handle cleanup; no worker remained after completion. The subsequent eight-test generic/UV-comparator recheck passes cleanly in 0.160s, including the numeric spectator-factor assertion (run `c3645c71-afca-419a-b7a5-d9ba73bf0cea`; `/tmp/gammaloop-phase-review-canonical-cleanup-check.log`). Final full-suite, Python/parity and scalar-matrix gates remain pending on this dependency pin.

- The first full acceptance on the generic index repair passed fresh Python/API model parity, exact `just test_gammaloop` (2,084/2,084; run `010d820d-f6db-4582-91b8-14b77c7b3344`) and all four explicitly enabled fast PySecDec tests (run `839e0276-0d26-475b-a35a-15d0901ea515`). The full scalar matrix then failed at GL29 after 146 passes, leaving 19 cases unrun (run `8d8fbcb6-eb50-44f1-8852-c344155a9f43`). This was a generation failure, before numerical assertions: an exactly-zero external coefficient survived angular projection in a factorized form and induced spurious CFF energy bounds. The first differing boundary is Vakint scalar-coefficient accumulation, before integrated-localizer CFF construction. Graph selection and Feyngen lookup are excluded; no CFF rank limit or physical reference is changed.
- `ScalarTerms::insert` now uses the existing equal-numerical-coefficient collector on a temporary candidate solely to certify zero, with a cycle guard. A nonzero stored coefficient retains its original AST. This avoids numerator expansion, polynomial conversion and numerical-coefficient GCDs, including their problematic complex-coefficient path. The regression reproduces the mixed angular sector while preserving a nonzero product-of-sums spectator; another covers cycles and accumulated cancellation with rational, complex and floating coefficients. The two tests supplement the generic index regressions rather than specializing the implementation to GL29 or ttH.
- Current-source test-target checking (9.04s), workspace/all-target Clippy with UFO support and warnings denied (15.81s), and the focused 28-test native run all pass. The latter includes GL29, complete slow ttH, Vakint Lorentz validation and tensor reduction (run `4ecf08c9-2923-481e-bcac-47f6e5592984`; 346.363s execution; watchdog peak 5,405,081,728 bytes; `/tmp/gammaloop-phase-review-zero-coefficient-safe-v2-targeted.log`). Fresh API/parity, the exact suite, enabled fast PySecDec and the complete 166-case scalar matrix must pass again on this final zero-coefficient repair before publishing the root commit.

- The refreshed Python/API build and SM UFO/JSON parity pass (`public-api-sm-parity-f1wjcf6w/report.json`), exact `just test_gammaloop` passes 2,086/2,086 with 296 existing skips (run `c6082736-0127-42af-9082-3d491e08f3a5`, 248.150s), and explicitly enabled fast PySecDec passes 4/4 (run `b9d9b6ad-cbdb-4534-9e29-d5d0631c2468`, 49.627s). Nextest reported one successful test with delayed output-handle cleanup; no worker survived the run, and the standard output settings did not retain its identity. This was not an assertion failure.
- The complete scalar release gate passed 152 cases, including GL29 and every quadratic/quartic probe, before GL35 failed with the same outer-CFF energy bounds (run `b33668d4-ebc4-49ef-8ba9-e2b3a0c81ad4`, 28m22s build, 1,810.805s execution, peak 20,215,427,072 bytes). A separate run of the 13 previously unrun cases passed ten and exposed the same generation failure in GL38, GL40 and GL46. All 166 cases have now been exercised on this checkpoint; four fail, so publication remains blocked on a code repair, not approval.
- Structured traces locate all four failures at the same scalar-coefficient boundary. With A=p3², B=p3·p4, C=p4², GL35 adds the zero `2*(-A-B)+2*(-B-C)+2*A+4*B+2*C` inside complex epsilon/logarithmic prefactors. Equal-coefficient collection reaches a nonzero-looking fixed point. GL38 and GL46 have exactly the same tensor-reduced numerator; GL40 combines the same zeros with surviving mass and loop terms. The exact numerical-content collector proves this real-rational zero without expansion. No additional phase, tensor, graph-selection or CFF-sector defect is implicated.
- Scalar-coefficient accumulation now certifies additive subexpressions from the inside out, committing only exact-zero replacements and retaining surviving factors. Its scratch candidates use numerical-content collection only when arithmetic numeric leaves are exact real rationals and arithmetic powers are positive integers; functions remain opaque. Other domains retain equal-coefficient grouping, and every subtree has a cycle guard. This avoids unsafe complex GCDs and powers that could introduce complex/infinite coefficients. A new regression covers both a zero nested inside a surviving complex coefficient and the complete mixed angular product; nonzero complex, floating, inverse and fractional-power expressions retain their original ASTs. Checking, linting and native acceptance must be refreshed on this repair.

- The guarded bottom-up zero repair passes test-target checking (9.10s), all-target Clippy with UFO support and warnings denied (14.73s), and all 33 focused native tests (run `61226ff3-e000-484d-a9ed-80e3a55e9f42`, 345.256s execution, peak 5,439,783,552 bytes; `/tmp/gammaloop-phase-review-zero-subtree-safe-targeted.log`). This includes GL29, GL35, GL38, GL40, GL46, all relevant Vakint projection/Lorentz/tensor tests and complete slow ttH. The source is frozen for a fresh full API/parity, exact suite, enabled fast PySecDec and complete scalar-matrix acceptance.

- Final acceptance of the guarded zero-subtree and generic index repairs passes on the same frozen source. The rebuilt Python API passes SM UFO/JSON parity, including the real-input declarations and complex CKM values (`/tmp/gammaloop-phase-plan/public-api-sm-parity-4gadzay7/report.json`). Exact `just test_gammaloop`, without added arguments or exclusions, passes 2,087/2,087 tests with 296 existing skips (run `acd44189-9934-4810-8877-6016d7f566f5`, 237.100s). All four signed photon/lepton ddx/ttx LO/NLO acceptances pass at their existing MC statistics. Nextest reports two passed-test output-handle cleanup notes, with no assertion failures or surviving workers.
- The explicitly enabled fast PySecDec gate passes 4/4 using fresh packages (run `077046ee-e24e-4817-b4de-66c61b209421`, 47.855s), preserving the default disabled state. Complete `just test_LU_scalar_xs` passes 166/166, with 105 cases outside that recipe skipped (run `3a409ab2-ad07-4e4e-a26d-d8ded29b7125`, 23m08s release build, 1,966.289s execution, watchdog peak 20,190,474,112 bytes). This includes every quadratic, quartic and mixed-angular case, including GL29, GL35, GL38, GL40 and GL46. The driver exits successfully and preserves all JUnit reports under `/tmp/gammaloop-phase-plan/zero-subtree-{focused,final-exact,final-fast-pysecdec,final-scalar}-junit.xml`; final logs use `/tmp/gammaloop-phase-review-zero-subtree-final-*.log`.
- The repairs remain at their shared tensor and coefficient owners: no graph-specific branch, weakened comparison, changed numerical reference, numerator expansion, or CFF fallback is introduced. Independent publication review confirms the amplitude target remains `-i M`, CP side symmetrization remains opt-in, and the PR describes the supported complex-coupling contract accurately. All required local test gates are green; publish the reviewed change to PR #103 and monitor checks on its new head every thirty minutes before declaring CI readiness.


### Complete d-dimensional algebra of analytically integrated UV subgraphs

The goal of this follow-up is to finish every Lorentz and ordinary Dirac contraction of the complete analytically integrated subgraph at symbolic d=4-2epsilon before Laurent coefficients are discarded, including contractions introduced by vacuum tensor projection. Expansion is explicitly permitted and desired inside that analytic subgraph; retained cograph and numerical-production numerators remain factorized. Reject gamma5 anywhere in the analytic UV subgraph before simplification can erase it, with a contextual error. Gamma5 prescriptions are out of scope.

The shared input already reconstructs symbolic Feynman-rule tensors and maps their Minkowski slots to d before GammaLoop evaluates ordinary closed fermion traces. The first missing boundary is after Vakint angular projection: previously opaque open spin tensors can acquire contracted slots, but the current return path merely collects chains and later specializes remaining slots to four dimensions. A second GammaLoop algebra pass must therefore occur before Vakint scalar integration performs coefficient-order selection and Laurent truncation. Preserve separate denominator topologies and complete loop-momentum mappings throughout this round trip.

Implement and vet three coordinated pieces: (1) restore projected tensors, complete their d-dimensional gamma/metric algebra and re-encode the numerator before scalar backend evaluation; (2) audit coefficient and integral epsilon valuations, nested-forest composition and backend depth limits so pole multiplication cannot discard evanescent contributions, with explicit errors where sufficient orders cannot be supplied; (3) reject gamma5 and implicit gamma5-bearing chiral tensors at the earliest shared analytic-subgraph boundary, without rejecting an unrelated cograph. Use independent cross-audits and exact single/double-pole regressions, ordinary closed traces, new projector contractions and gamma5 scope tests. Run formatting/checking and Clippy, focused native tests, then the required full gates on the final tree. Previous 26b57b786 results do not certify these changes; PR #103 stays draft during implementation and is republished only after verification.

User clarifications (verbatim):

```text
Also note tha tin general my statement about the numerator not being expanded always has the one exception of that of the subgraph being integrated over analytically, there the expansion is necessary and desired precisely because it needs to be done in d-dimension, not 4, so it's necessary.
```

```text
the tensor reduction in vaking is done in d-dimension, but any trace of a closed fermion loop for instance, would need to be done by gammaloop before supplying the input to vakint, and that trace will need to be computed in d-dimension!
```

```text
Ok, use multiple subagent to thoroughly edit and vet the d-dimensional algebra needed within the complete analytically integrated subgraph. Numerator expansion is ok there, but we must really never afford to drop a d-dimensional term, i.e. the complete epsilon-expansion must be fully resolved before any truncation happen.
If a gamma5 is found anywhere in the UV subgraph being analytically integrated over for now, error cleanly (this is out of scope for now)
```

```text
The complete numerator of the UV subgraph being integrated over by vaking can be expanded, it's fine, it needs to because it must be done in d-dimension. Any gamma5 within that UV subgraph should error cleanly for now (out of scope).
I believe there is no guarantee for now that he whole algebra is done in (d=4-2eps)-dimension within that subgraph, in particular after injecting vakint's result in there.
Fix and review this thoroughly with subagents, and also tell me why in the first place in this PR did you change the way you supplied the inputs to vakint (no longer one monolithic input).
```

The earlier backend split was motivated by an overly broad interpretation of the numerator-expansion prohibition: FORM distributed complete coefficients even when no Rust-side expansion was called. GammaLoop still supplied a complete Vakint expression; Vakint introduced universal tensor kernels and scalar coefficient aliases to protect those factors. The user has clarified that this restriction does not apply inside the analytically integrated subgraph. Dimensional algebra completion takes precedence there. The retained cograph remains outside the analytic input, and its factorization remains required.

- The independent epsilon audit found an additional real order-budget defect: poles in a custom loop normalization were applied after the raw backend had already truncated. The repair evaluates the normalization at the actual topology loop count, extracts its epsilon pole, requests the necessary additional raw orders, and restores the pole before the final Laurent projection. Numerator and normalization pole costs add. Analytic and explicitly enabled fast PySecDec regressions cover both sources, their combination, and quadratic evanescence. Supported mathematical coefficient functions use the existing series algebra; unresolved epsilon-dependent coefficients fail instead of passing through truncation.
- The suspected disconnected-union depth defect was excluded by tracing actual consumers: later analytic UV operations use the factorized local value built from each component’s full global-depth recursion input, not the aggregate union summary. The shorter aggregate only feeds numerical 3D replay/final assembly, and terminal renormalization likewise uses individual component projections. No forest-order change was warranted.
- CI on the previously published commit 26b57b786 is fully green (83 checks, including the deployment marker, checked 2026-09-08 07:59 UTC). This does not certify the uncommitted dimensional-algebra follow-up. The current PR remains draft until the new source passes its gates.


### Whole-numerator Vakint comparison mode

Keep the pre-existing conservative forest epsilon bound; positive epsilon prefactors never reduce requested orders. Add global.generation.uv.project_integrated_uv_cts_onto_tensor_integrals, true by default. True keeps the kernel/opaque-coefficient route. False restores the base branch's complete-numerator FORM tensor reduction and complete scalar-backend input per denominator topology. The base already kept distinct denominator topologies separate. Both routes share raw gamma5 rejection, pre-reduction and post-reduction d-dimensional algebra, sufficient conservative epsilon orders, unresolved tensor errors, and final truncation/restoration. Add exact pole/finite parity tests, including newly projected gamma contractions and scalar/closed-trace examples. Test configuration defaults, false round-trip and SHOWDEFAULTS exposure; regenerate schemas.

User instructions (verbatim):

```text
I don't want you to tune the depth of vakint epsilon order asked for, just always set it to the more conservative value; even if there is an order-epsilon prefactor, it's ok to overshoot the target; it's ok the truncation will eventually take care of this. So there was a conservative bound before which you should keep as-is now. Confirm you understand.
```

```text
Also I actually fear that you do something wrong here, so I want you to add a new generation mode in the uv settings, called "project_onto_tensor_integrals" (set to true by default, in which case you do what you're doing now), but then when set to false, you should do essentially the equivalent of what was done before (look into the branch you're based on to verify how it was done), i.e. pass the whole integrand at once to vakint at once.
In both cases I want you to make sure d-dimensional algebra is done consistently, and also add a couple tests to verify that both modes agree.
```

The user-facing UV option is now named `project_integrated_uv_cts_onto_tensor_integrals`, with no compatibility alias. Vakint's generic internal tensor-reduction strategy setting retains its library-level name. All local gates must pass on the final renamed implementation before committing, pushing and resuming the existing thirty-minute PR monitor.

```text
Continue as planned, but rename the option "`project_onto_tensor_integrals`" to "`project_integrated_uv_cts_onto_tensor_integrals`"
```

```text
Also make sure you have all local gates passing before commiting+pushing and monitoring the PR at a 30min interval rate.
```

- The first both-mode focused run executes 76 tests: 70 pass and six fail (26m28s build, 195.561s execution, peak 9,786,634,816 bytes; `/tmp/gammaloop-phase-d-dimensional-targeted.log`, saved `d-dimensional-mode-baseline-junit.xml`). Complete slow ttH passes. The failing boundaries are compact Lorentz materialization with explicit spectator spin slots, an exact rational-input comparison lacking denominator cancellation, repeated inline snapshot use, and monolithic normalization/complex-number serialization. These are not green local gates.
- Spenso's shared materializer now identifies exactly one unindexed compact axis while preserving explicit spectator slots. Analytic UV validation rejects remaining ambiguous Lorentz products. The input oracle uses unfactorized exact rational cancellation; the quark test preserves its signed snapshot under the test framework's repeated-assertion guard. The scalar parity oracle keeps transcendental master constants symbolic, so it compares exact algebra without floating accumulation order.
- A retained-binary FORM diagnostic locates the monolithic tensor failure before reduction: `(1+2i)*(a+b)*tensor*loop_momenta` was serialized as `(a+b)*1+2*i_*tensor*loop_momenta`, changing precedence and creating scalar summands with no tensor slots. The repair belongs to canonical backend serialization, preserving full symbol attributes and bypassing custom display printers. Standard backend normalizations are built as exponentials of their epsilon-weighted formal real logarithms before user coefficients are combined; arbitrary complex coefficient functions and custom normalizations retain their branches. No production reference or phase is adjusted to obtain parity.

- The second focused checkpoint passes 106/109 tests, including complete slow ttH (190.010s; `/tmp/gammaloop-phase-d-dimensional-v2-targeted.log`, saved `d-dimensional-v2-baseline-junit.xml`). Three remaining failures are isolated at representation boundaries: custom display callbacks corrupt canonical backend identifiers, the whole-input FORM route cancels rational denominator factors that simple expansion does not cancel, and the quark oracle assumes an explicit `log(μ_R²)` term. Captured FORM inputs and outputs prove that both tensor routes retain the same complex coefficient and projector: `2/(2ε−4) = 1/(ε−2)`.
- The backend serializer now follows the existing expression tree with raw symbol names, attributes and tags, preserving complex parentheses and antisymmetric signs without display callbacks or temporary wrapper symbols. Both tensor modes are compared against the same exact `k²/d` projector over complex rationals. The quark scale-log coefficient is extracted as `μ_R² ∂/∂μ_R²`, retaining the original signed pole residue and normalization without imposing logarithm branch rewrites. These changes require a fresh native run and the final full gates; no commit or push has occurred.

- The completed final focused rerun passes all 109 tests, including complete slow ttH, both-mode pole/finite comparisons, projected Dirac evanescence, callback/antisymmetric serialization and compact-tensor guards (run `04707236-e43e-43c9-bcce-80e6bfe6b100`, 1m18s build, 190.664s execution; `/tmp/gammaloop-phase-d-dimensional-v4-targeted.log`, saved `d-dimensional-v4-focused-junit.xml`). Formatting, Hakari verification, test-target checking (9.64s) and workspace/all-target Clippy with warnings denied (16.51s) pass on this source. The preceding 108/109 checkpoint differed only by a malformed namespace in the new test fixture, now corrected; no production adjustment was required. Fresh API/schema/parity, exact full suite, explicitly enabled fast PySecDec and complete scalar LU validation follow before commit/push.


### Manual-only numerical PySecDec validation

All 24 physical PySecDec tests are ignored and explicitly excluded from automatic nextest profiles and the `just test_gammaloop` recipe, including its slow/failing selections. Running these tests manually requires both an explicit ignored-test selection and `RUN_PYSECDEC_TESTS=1`; missing opt-in or a missing backend now fails clearly instead of reporting a skipped computation as a pass. Automatic analytic tests cannot select PySecDec as an installed fallback. Pure Rust adapter tests remain automatic because they do not launch PySecDec. Run only the four fast numerical validations once for this backend change; do not add them to CI or the PR heartbeat.

User instruction (verbatim):

```text
Continue as planned, but make sure that the pySecDec tests are really just run now because you modified something related to them, but they will not appear in CI or in any automatic test coverage in the future.
```

- The full d-dimensional checkpoint completed with 2,104 passes and two failures (run `cf9dc267-c6b9-4a74-96f7-d901759f7561`, 256.852s execution). All four signed photon/lepton ddx/ttx LO/NLO acceptances pass. The analytic one-loop reference differed only by equivalent mass/scale logarithm notation; its unchanged Laurent coefficients are now written in the dimensionless-ratio basis and independently verified through epsilon squared. The scalar-sunrise UV orchestrator comparison still requires first-difference diagnosis. These results do not certify the final tree.

- Sunrise triage locates the first unequal expression at the direct two-loop analytically integrated node; the other six integrated outputs match. An exact symbolic check proves that full analytic expansion makes the two expressions identical, whereas numerical-coefficient expansion and factor collection do not. The common master constants were treated as independent symbols, without relying on rounding. Normalize the completed analytic subgraph at `Integrated::integrate` before Laurent projection and cograph reinsertion. The existing factorized comparator remains unchanged and now records a structured expression pair on mismatch. The scalar finite/pole regression compares both forest orchestrators in both Vakint input modes. All three reviewers approve this bounded repair; native validation follows.

- Final focused validation of the analytic-return repair passes 111/111 tests, including unchanged sunrise and slow ttH acceptance, both forest owners, and both Vakint modes (run `60625d82-21a1-4000-8868-3d7ee0dcca57`, 11m33s build, 205.959s execution, peak 8,109,173,760 bytes). Formatting, Hakari, test-target checking (9.81s), and all-target Clippy (16.80s) pass. One output-handle cleanup note is reported; no test failure or surviving worker is found. The list-only selection audit also passes: all seven automatic profiles and the slow-test bypass exclude the 24 physical PySecDec tests even with the opt-in set; manual inventories are exactly 24 ignored tests/four fast checks, and three pure Rust adapter tests remain automatic. Final full gates continue on the same frozen sources.

- The refreshed API and SM UFO/JSON parity pass (`public-api-sm-parity-la_86efs/report.json`), including persisted defaults and the explicit false mode. The next exact suite runs 2,082 tests: 2,081 pass and `dod0_bubble_uv` fails only its pointwise RLS diagnostic (run `1161550c-d4d2-4e5a-b883-3718163c62bc`, 255.466s). All four signed ddx/ttx acceptances, the integrated RLS-invariance check, and the bubble physical target pass. At this fixture's equal vacuum/renormalization scales the finite counterterm vanishes, so pointwise RLS independence is allowed; the test incorrectly demanded exact floating-point equality. Its best relative difference is 2.365e-16, far below the reported 1e-10 numerical accuracy. Correct the shared numerical test classification to recognize accuracy-compatible absence, retain the strict 1000-times-accuracy requirement for detected dependence and the unresolved band between them, and reject nonfinite probes before ranking. No production algebra or physical reference is changed.

- The accuracy-aware test-only repair passes all 113 focused tests, including the unchanged bubble and sunrise physical targets, both-mode UV comparisons, the new numerical-classification regression, and slow ttH (run `e9a4b435-7c52-48c5-826a-c48783a796bb`, 26.10s build, 211.571s execution). Formatting, Hakari, test-target checking (1.75s), and all-target Clippy (1.38s) pass. One output-handle cleanup note accompanies the successful run. No production code changed; the fresh API/schema/parity evidence remains current. Exact full-suite, manual numerical PySecDec, and scalar LU gates follow on the frozen tree.

- Exact `just test_gammaloop` now passes 2,083/2,083 tests on the frozen source (run `1139c626-f800-4919-99d6-e6a5c6776112`, 260.470s), including all four signed photon/lepton ddx/ttx LO/NLO acceptances. The separately invoked manual PySecDec gate passes three tests and fails the complex-parameter test on its first whole-input case, after all five projected-mode cases pass (run `c4cd9aae-e9ad-4fb3-8d09-516528fe2b59`, 87.380s). The complete numeric coefficient is serialized as a decimal complex literal; downstream FORM rejects the decimal tokens before numerical integration. Preserve the whole polynomial and its represented coefficient exactly using the existing lossless finite-float-to-rational conversion at that backend boundary. Validate a copied retained package and the manual cases before refreshing final gates. Numerical PySecDec remains excluded from automatic coverage; no scalar release run or publication has started.

- Controlled copied-package tests isolate two PySecDec input boundaries. Exact complex rationals fix FORM parsing; with no complex parameter aliases, a literal imaginary coefficient additionally requires `enforce_complex` on the parsed numerator. Both copied generation/compilation/integration runs pass for `1+2i` and the exact IEEE representation of `-0.1+3i/8`, preserving their signed real/imaginary pole and finite parts. The implementation uses existing `Rational::try_from(f64)` without additional rounding, and derives the complex-return flag from PySecDec's parsed numerator. A pure Rust test checks the whole polynomial, multiple monomials, exact coefficients, epsilon powers, and absence of aliases. Both reviewers approve the two-file repair; physical references, epsilon depths, and manual-only selection remain unchanged. Run the manual fast checks first after rebuilding, then refresh the broader final gates.

### Final manual PySecDec validation on the reviewed d-dimensional source

The four explicitly selected fast numerical tests pass (run `a7435d78-90b4-448f-b753-eefd844aba41`, 111.527 s), including the complex/signed coefficient cases in projected and whole-numerator modes. The repaired whole-numerator boundary preserves finite IEEE coefficients as exact rationals for FORM and requests a complex generated return type when a literal imaginary coefficient remains. Their numerical references and conservative epsilon requests are unchanged. All 24 physical numerical tests remain ignored and excluded from automatic profiles and recipes; the list-only selection audit passes. Current-source check, Clippy, formatting, and Hakari pass. Focused UV, refreshed API/schema/parity, exact full suite, and scalar LU validation must finish before publication.

- On the same frozen source, all 114 focused tests pass (`cd81102f-33de-419f-a823-c3d4b9daf529`, 203.812 s), followed by a fresh API build (67.623 s), schema/default/explicit-false persistence checks and SM UFO/JSON parity (`public-api-sm-parity-rckwksou/report.json`). Exact `just test_gammaloop` passes 2,084/2,084 tests (`7cf84abe-3538-4987-af8a-60912af801ec`, 258.224 s), including all four signed photon/lepton ddx/ttx LO/NLO acceptances. None of the physical numerical PySecDec tests enters that suite. The separate 166-case scalar LU release gate follows before publication.

### Completed final d-dimensional local gates

All 166 scalar LU cases pass on the reviewed source (`ebc0888f-b80a-490c-a255-3d1a5caa20c3`, 25m17s release build, 1971.052 s test execution). This includes the full quadratic/quartic numerator matrix and GL35/38/40/46; peak process-tree memory including compilation is 20,008,611,328 bytes, below the 30 GB watchdog cap. No reference, phase, evaluator path or assertion was adjusted during validation.

The final evidence collector verifies unchanged source/gate hashes, zero JUnit failures, 114 focused tests, 2,084 exact `just test_gammaloop` tests, all four signed photon/lepton ddx/ttx LO/NLO acceptances, the four explicitly invoked numerical PySecDec tests, and all 166 scalar LU cases. Fresh API/schema/persistence and UFO/JSON parity, formatting, Hakari, test-target checking and all-target Clippy also pass. Numerical PySecDec remains manual-only and is excluded from CI, automatic profiles and the PR heartbeat. Results are retained in `/tmp/gammaloop-phase-plan/d-dimensional-final-evidence.json` and the referenced JUnit/log files. All local gates are green before the follow-up commit and push; remote CI must certify the resulting commit separately.
