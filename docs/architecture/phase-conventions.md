# Phase conventions and forward-cut conjugation

This note fixes the convention for GammaLoop amplitudes and LU cross-sections and records the independent audit of right-side conjugation. It is intended for GammaLoop contributors. The supported physical setting is a Hermitian UFO theory with real propagator masses. The complex-mass scheme requires additional machinery and is outside this implementation. The companion [SM convention audit](sm-conventions-audit.md) covers derivative vertices, scalar/pseudoscalar and vector/axial structures, fermion momentum signs and the separate W/Z gauge limitation.

## Amplitudes

GammaLoop returns the complete amputated Feynman graph, including vertex and propagator numerator factors, with the measure

\[
\prod_{\ell=1}^{L}\frac{d^4k_\ell}{(2\pi)^4}.
\]

We call this quantity \(A_\Gamma=-i\mathcal M_\Gamma\), following the requested convention. Its operational definition is unambiguous: for \(\mathcal L_{\rm int}=-\lambda\phi^4/4!\), the contact graph is \(-i\lambda\), and the scalar propagator is \(i/(q^2-m^2+i0)\). Some texts call the same complete Feynman-rule graph \(i\mathcal M\). [Tong, Feynman rules](https://www.damtp.cam.ac.uk/user/tong/qft/qfthtml/S3.html).

The native CFF source-sign adapter yields the energy representation with measure \(dq^0/(2\pi i)\). Its physical conversion is therefore \(i^L/(2\pi)^{3L}\). The existing real source-frame signs and inverse-energy products remain in place. Integrated UV addbacks must use the same conversion for the loops integrated into their coefficients; GammaLoop's additional Vakint normalization is consequently `1`.

For a complete literal numerator of one and total denominator power \(P\), a convergent Wick-rotatable integral satisfies

\[
I_M=i^L(-1)^P I_E.
\]

Thus absence of thresholds does not make all scalar topology integrals real. Their phase can alternate with loop count. In four-point quartic topologies \(P=2L\), the sequence is \(1,i,-1,-i,1,\ldots\). A full UFO amplitude also contains its vertex and propagator factors. Ultraviolet subtraction can change the sign of a finite real coefficient; absence of a threshold alone is not a positivity theorem for a renormalized amplitude.

## The right side is an inverse-process graph

Conjugating an independently constructed amplitude is necessary. Conjugating the raw right-hand numerator of a forward graph again generally acts on the wrong object.

Consider a vertex \(V=-ig\Gamma\), and define \(\bar\Gamma=\gamma^0\Gamma^\dagger\gamma^0\). In a Hermitian model, the inverse transition has vertex

\[
V_{\rm inv}=-ig^*\bar\Gamma,
\qquad
\bar V=+ig^*\bar\Gamma=-V_{\rm inv}.
\]

The inverse vertex already carries the conjugate physical coupling and appropriate spin structure. Another adjoint would give \(\bar V_{\rm inv}=+ig\Gamma\), undoing both transformations.

In GammaLoop this distinction occurs before numerator algebra: oriented particle extraction in `graph/parse/mod.rs` selects particle/antiparticle fields; `fix_cp_vertex_rules` selects the corresponding UFO vertex; numerator generation uses that vertex's actual coupling and its ordered spin/color slots. Cut spin sums and existing fermion statistics supply their own sewing factors. A generic matrix-adjoint identity does not certify that applying it at this point is physically correct.

The independent audit used explicit Weyl matrices and actual SM vertex pairs:

| Example | Independently sewn norm | Adjoint applied again to raw RHS |
| --- | ---: | ---: |
| Charged scalar partners \(+y_\tau P_R\), \(-y_\tau P_L\), normalized test kinematics | \(21\) | \(0\) |
| Charged current with \(V=1+2i\), common real weak factor stripped | \(105\) | \(-63+84i\) |

The first failure is a chirality error and cannot be repaired with an overall phase. Real scalar and ordinary real-coupling vector vertices hide this problem because their tensor adjoints leave their spin structures unchanged.

## Marking, contours and cut lines

The largest-time construction marks the right-side vertices with a minus sign while retaining their inverse-process values. Marked virtual propagators have the anti-causal prescription. Hermiticity then identifies the marked region with the conjugated amplitude; this statement does not assume CP invariance. [Anselmi, sections 2.2–3](https://arxiv.org/pdf/1606.06348).

For standard scalar, Dirac and covariant vector propagators, the complete phase accounting at a physical source graph is:

| Source of phase | Right-side conversion |
| --- | --- |
| Each inverse-process interaction vertex | \(-1\) |
| Each uncut propagator's explicit numerator \(i\) | \(-1\), retaining its physical tensor and reversing \(i0\) |
| Each virtual energy contour, expressed through the common causal CFF representation | \(-1\) |
| Physical cut propagator | Replace the complete Feynman propagator by the real positive-energy on-shell measure |

Since \(L_R=I_R-V_R+C_R\), the first three rows combine to

\[
(-1)^{V_R+I_R+L_R}=(-1)^{C_R}.
\]

Here \(C_R\) counts the connected interacting components on the right, with cut boundary hairs retained. Removing those hairs before counting would incorrectly erase a contact vertex. A contact, an exchange tree and a connected virtual-loop correction all have \(C_R=1\). Derivative-created copies of raised denominators are not extra action-level propagators for this count.

For a scalar line, the distribution identity is

\[
\frac{i}{x+i0}-\frac{i}{x-i0}=2\pi\delta(x),
\qquad
\frac{i}{q^2-m^2+i0}\ \longrightarrow\ 2\pi\theta(q^0)\delta(q^2-m^2).
\]

At the positive-energy pole this gives \(2\pi/(2E)\). Equivalently, the clockwise cut factor \(-2\pi i\) multiplies the retained pole numerator \(i/(2E)\). There is no residual cut-line \(i\).

GammaLoop has already integrated the forward energies and retained the cut-line numerator factors in its CFF residue. Its LU conversion must account for that representation exactly once: the common residue conversion is \(2\pi i\), and the right-side marking/contour factor gives \(2\pi i(-1)^{C_R}\), hence \(-2\pi i\) for connected RHS. This is independent of forward loop count and cut multiplicity. It is not an additional replacement applied independently to each cut edge after CFF reduction.

The same cut-group factor multiplies bare, local UV, integrated UV, and threshold terms. The left CFF threshold has the prescription `eta-i0`, with integrated contribution `+i*pi*delta(eta)`; the right has its complex conjugate. Runtime subtracts the single-threshold helpers, so their integrated coefficients are respectively `-i*pi` and `+i*pi`. Iterated helpers use the same factors with the existing inclusion-exclusion signs; no additional conjugation is applied. Grouped cuts must have consistent marking parity. The phase formula does not by itself implement disconnected amplitudes requiring several independent energy constraints; existing CFF surface support still applies.

Computed UV-forest exports contain the same finalized physical expressions as production. They reuse the cut-group conversion and must not divide the spatial measure a second time. Comparisons for both UV orchestrators preserve the factorized expression trees: they use structural equality where possible and exact signed values at sampled nonsingular rational points otherwise. These sampled checks are not a general symbolic identity proof. The independent scalar Born residue normalization is checked structurally with \(\pi\) left symbolic.

## Complex couplings and the complex-mass scheme

Intrinsic complex couplings in a Hermitian interaction are compatible with marking. In the vendored SM, the unique conjugate-field partners `V_125` and `V_95` contain \(iV_{ub}\) and \(iV_{ub}^*\), respectively. The parameter \(V_{ub}=A\lambda^3(\rho-i\eta)\) is complex even though its Wolfenstein inputs are real. Leaving the inverse vertex value intact preserves \(|V_{ub}|^2\). Testing whether a numerical UFO coupling has a nonzero imaginary part would wrongly reject conventional Feynman-rule factors of \(i\). Model Hermiticity and correct partner assignment, rather than that numerical test, are the relevant assumptions.

The optional `symmetrize_left_right_states` setting performs an additional CP-based graph transformation. It defaults to false. A guard rejects this optimization when a used vertex coupling depends on a complex-declared parameter with a nonreal or unresolved value. It follows coupling aliases, ignores unrelated parameters, and does not treat an explicit Feynman-rule `i` as an intrinsic phase. The assumption is retained in the generated integrand and rechecked after model updates. Disable the optimization and regenerate to use ordinary Hermitian sewing with complex CKM. Passing this guard is not a proof of CP invariance for every possible UFO interaction: a phase written directly in a coupling, such as `i*g*exp(i*theta)` with real `g` and `theta`, is not diagnosed. Such models must leave this optional optimization disabled; ordinary sewing still requires the Hermitian partner vertex.

CMS analytically continues masses and derived parameters. Its renormalized action alone is non-Hermitian; for example, the same complex field-renormalization constant applies to both charged W fields. Charge reversal therefore does not supply all conjugations required by amplitude sewing. [Denner–Dittmaier, equations (1)–(4)](https://arxiv.org/pdf/hep-ph/0605312).

An exact squared resonant amplitude contains conjugate denominators \((s-\mu^2)^{-1}(s-\mu^{2*})^{-1}\). [Bauer et al., equation (30)](https://arxiv.org/pdf/1211.1684). If a coefficient is \(Vf(\mu^2)\), the inverse vertex can contain \(V^*f(\mu^2)\), whereas its conjugated amplitude needs \(V^*f(\mu^2)^*\). Conjugating the entire inverse coefficient would incorrectly undo the CKM conjugation.

Ordinary stable-particle Cutkosky delta functions are insufficient for finite-width propagators. CMS unitarity requires the appropriate mass counterterms and cuts through stable states; its fixed-order relation holds up to omitted perturbative orders. [Denner–Lang, equations (4.4), (4.11)–(4.12), and section 4.3](https://arxiv.org/pdf/1406.6280).

GammaLoop currently constructs real on-shell energies and does not implement those complex poles or CMS cutting rules. Actual complex masses must produce a clear error before cross-section generation or evaluation, including after model updates. A UFO's width metadata alone does not activate CMS in this implementation. No tensor-conjugation switch can supply the missing analytic continuation and counterterm prescription.

The mass check does not prove that an arbitrary UFO is Hermitian. Unpaired complex interaction coefficients, including CMS-derived coefficients in a graph whose propagator masses happen to be real, violate the supported model contract and are not comprehensively diagnosed by this check.

## Acceptance and user-supplied numerators

For symmetric scalar Born halves, the full cut result must equal \(|A_\Gamma|^2d\Phi_n\), with physical real flux and symmetry factors. The independent massless phase-space oracle is

\[
\Phi_n(s)=\frac{s^{n-2}}{2(4\pi)^{2n-3}\Gamma(n)\Gamma(n-1)}.
\]

The acceptance suite covers multiplicities two through six, massive decay, two incoming particles, an exchange tree and left/right virtual-loop mirrors. Separate signed acceptances compare both components of the four photon/lepton ddx/ttx channels, including negative and positive NLO graph contributions. Analytic scalar threshold examples test discontinuity signs separately from Born positivity. Execution status and actual numerical results belong in `PHASE_CONVENTIONS_FIX.md`.

An explicit complete global numerator retains its full-forward meaning; it does not specify two separately supplied amplitude factors. A literal numerator-one topology need not have the phase or positivity of a squared physical UFO amplitude. No new `exact`/`scalar_only` modes are introduced. Regenerate integrands and saved compiled states after the convention change; old numerical artifacts contain the old phases.
