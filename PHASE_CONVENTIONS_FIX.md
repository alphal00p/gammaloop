> Implementation and the requested phase, full-suite and slow ttH validation are complete.
> Exact `just test_gammaloop`: 2,053 passed, zero failed, 296 existing skips.
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
