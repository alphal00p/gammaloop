= Five-loop quark-ladder transverse coefficient

The compact seed #raw("ladder.dot") specifies one closed massless down-quark loop,
four noncrossing gluon rungs, and two external gluons. Native GammaLoop/FeynKit
generation selects its exact colored graph and instantiates the existing QCD
model rules; importing the bare seed alone does not dress its numerator.
The conventions are Feynman gauge, fixed SU(3), and Dirac trace of unity four.
Every Lorentz representation is continued to generic #raw("rustred::d"); color
and spinor dimensions are not changed to that dimension.

The calculated quantity is the coefficient of the transverse quadratic
two-point structure for this single graph. It is not a complete beta function,
a gauge-complete graph sum, or an R-subtracted result. The contraction is
#raw("delta_ab * (g_munu - q_mu*q_nu/q^2) / [8*(d-1)*q^2]").
Both quadratic tensor structures are retained, so transversality of the
individual graph is never assumed. Native controls retain the coefficient of
#raw("A*q^2*g_munu") and remove #raw("B*q_mu*q_nu").

The homogeneous quadratic coefficient is extracted from the complete numerator
and propagator product with Symbolica's native series API. This is the logarithmic
Taylor term of eq. (4.16) of
#link("https://arxiv.org/abs/2203.11038")[Local Unitarity: cutting raised propagators and localising renormalisation].
As in eqs. (4.17–4.18), a common UV mass enters only the vacuum denominators;
the original fermion slash numerators remain massless. All raised-denominator
terms are retained. An auxiliary fixed-radius angular average of the external
direction uses the existing FeynKit rank-two/rank-four tensor reducer. It adds
no physical integration loop and no extra loop-measure normalization.

The five native upper-rail chords keep #raw("p0,...,p4") compact during Dirac
algebra, before expansion into a vacuum reduction basis. Removing the two
external carriers merges two equal fermion pairs at zero external momentum:
there are five loops, twelve vacuum incidence edges, two insertion dots, and initial
total denominator power fourteen. These twelve edges are not twelve distinct
quadratic forms: the three other upper/lower fermion pairs also coincide at zero
external momentum. Native exact denominator comparison and graph ingress verify
nine distinct forms, with five doubled fermion lines and four gluon lines, while
preserving all five loops and the numerator. Quadratic Taylor terms have total powers
fourteen through sixteen and physical #raw("A-R=10"). After exact scalar lowering,
the report distinguishes these physical degrees from the actual shifted-basis
keys and explicit UV-mass powers; no higher key is declared an admitted
descendant without a recorded scope witness.

The full preparation runs inside the ignored #raw("slow::") test and is excluded
from the ordinary test selection. The preparation command shares the same helper:

#raw("cargo run --release -p vakint --no-default-features --features feynkit-ingress,rustred-app/capacity-dispatch --example five_loop_qcd_ladder -- TMP/qcd-ladder", block: true)

#raw("--dress-only") performs only the native graph/rule preparation. The slow
test requires the explicitly configured #raw("VAKINT_SAVED_APPLICATION") path;
#raw("VAKINT_QCD_LADDER_EVIDENCE") optionally selects its evidence directory.
The test computes the numerator anew, reports preparation time separately from
application time, and checks every surviving lowered key against the actual
saved entry union. A refusal is a failing acceptance test with a retained census,
not a successful reduction. No FORM fallback, new IBP generation, default
artifact download, or entry-domain widening is permitted.

Only the compact source fixture and reproducible code belong in version control.
Large dressed/projected numerators and diagnostic DOT output remain generated
workspace-TMP evidence, never a precomputed substitute for the slow test.

The separately ignored prepared-input diagnostic can set
#raw("VAKINT_QCD_DIAGNOSTIC_MAX_WORK") to an explicit operational work budget;
without it, the existing default of sixteen million charged work units remains.
The diagnostic records that policy and still checks all original visible input
poles, mass homogeneity, and the actual saved entry union without applying IBPs.
A raised-budget replay is not a default-limit success or fresh-preparation test.
Contribution ordinals are local to each run: compare the recorded family, key,
and exact numerator rather than assuming that equal ordinals identify equal
inputs. Its per-contribution histogram is also not a globally coalesced sum
or proof that every equivalent momentum presentation is unsupported.
