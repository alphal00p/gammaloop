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
This bare quadratic building block enters the R-operation formula; it is not
the complete recursive counterterm, which also requires the forest terms of
eq. (3.3) and the local subtraction prescription around eq. (4.16) below.

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
fourteen through sixteen and unlowered #raw("A-r_mom=10"). After exact scalar lowering,
the report distinguishes original momentum rank from the actual shifted-basis
keys and explicit UV-mass powers; no higher key is declared an admitted
descendant without a recorded scope witness.

Here #raw("r_mom=4,5,6") counts momentum scalar products before the denominator-basis
rewrite. The native signed-index rank #raw("r_idx") counts negative denominator powers.
For #raw("D_i=Q_i^2-M_UV^2"), exact substitution
#raw("Q_i^2=D_i+M_UV^2") can produce a rank-zero key with a mass prefactor.
A branch with coefficient #raw("(M_UV^2)^j") satisfies engineering mass-dimension
balance #raw("D_key-j=10"). This is not that branch's hard UV degree when the
UV mass is fixed. Its superficial degree at four dimensions is
#raw("omega_UV=20-2*D_key"); an #raw("A_key=16,r_idx=0,j=6") branch therefore
has #raw("omega_UV=-12"), whereas the original #raw("A=16,r_mom=6") source has
#raw("omega_UV=0"). Engineering balance does not establish saved starting-domain
membership. The fixture performs no optional numerator/denominator
cancellation; integer-index lowering is required by the saved application.
No mass branch is discarded solely because its overall UV degree is finite;
subdivergent poles must be assessed together with the complete forest formula.

The full preparation runs inside the ignored #raw("slow::") test and is excluded
from the ordinary test selection. The preparation command shares the same helper:

#raw("cargo run --release -p vakint --features feynkit-ingress --example five_loop_qcd_ladder -- --work-dir ../workspace/5l-push/quadratic --artifacts ../vakint-external-data", block: true)

#raw("--dress-only") performs native graph/rule preparation and validates the
vacuum incidence without gamma algebra. #raw("--prepare-only") additionally
recomputes the complete scalar input without applying IBPs. The default
#raw("--projection quadratic") preserves the full transverse coefficient above.
The additional #raw("--projection zero-momentum") diagnostic retains both current
vertices and contracts with unaveraged #raw("g_munu*delta_ab") at zero external
momentum. Its unlowered source has #raw("A=14, r_mom=5, D_source=9") and superficial
#raw("omega_UV=2") at four dimensions. It does not use a transverse division,
angular projector, or polarization average.
It is an additional diagnostic of the bare graph, rather than the full soft
mass term of the recursively subtracted counterterm.

The collector records native preparation stages, data loading, geometric routing,
weighted scalar lowering, rule application, terminal substitution, and total
elapsed time. Its input and result hashes, executable hash, source revision,
data identity, exact conditions, and cache counters are written alongside the
expressions in #raw("benchmark.json"). #raw("--cache PATH") explicitly attaches a
prepared five-loop cache; #raw("--warm") performs a second evaluation and requires
exact result, guard, and family parity. Repeated uncached/cached measurements use
separate fresh processes. Preparation is rerun for each complete measurement.
Uncached runs explicitly detach any backing attached by the data loader;
cached runs clear resident results while retaining the selected backing.
The receipt separates native owner preparation during loading and reduction,
and records persistent root and intermediate hits. An explicit
#raw("--scalar-work-limit N") changes only the finite scalar arithmetic budget;
the receipt records the chosen limits, so a raised-budget result is never
reported as a default-limit success.
The native angular reducer has an explicit one-million-term output allowance.
A preserved fresh run of the indexed adapter reached the native default
100,000-term limit at the second angular group; that operational failure
does not establish anything about IBP coverage or surviving terminals.

The slow test requires the explicit #raw("VAKINT_SAVED_APPLICATION") and
#raw("VAKINT_QCD_LADDER_EVIDENCE") paths. Its preparation census distinguishes a
verified geometric witness from the recorded starting predicates. Out-of-domain
inputs receive an advisory warning and may still reduce safely. A successful
primary benchmark must report any unknown residuals explicitly; it never changes
the scientific input or widens the campaign. No FORM fallback or new IBP
generation is permitted.

The separately ignored #raw("prepared_ladder_selected_campaign_census_without_ibps")
may replay an explicitly selected, hashed prepared input for early coverage and
cache-workload inspection. It uses the same selected routing as public evaluation,
retains every exact mass factor and globally combines roots before classifying
coverage. Its JSON clearly labels diagnostic replay and does not claim fresh
preparation or successful IBP reduction. The native saved predicates, rather than
rectangular rank/degree maxima, remain the admission authority.
The census can set #raw("VAKINT_QCD_CENSUS_MAX_TERM_OPERATIONS") to a finite work
ceiling. Completed checkpoints are addressed by the exact source identity,
bound to the input and application hashes and family fingerprint, and checked
again against the active admission predicates. Imported checkpoints retain
their original budget and policy receipt. Their local contribution ordinals
are presentation metadata and are reassigned for the current enumeration.
The default census requires all 21 quadratic groups and checks engineering
mass-dimension balance #raw("D_key-j=10"). For the additional diagnostic, explicitly selecting
#raw("VAKINT_QCD_CENSUS_PROJECTION=zero-momentum") requires one group and
#raw("D_key-j=9"); other or mixed projection values are rejected. These checks
do not certify individual-term UV degree or campaign coverage.

Only the compact source fixture and reproducible code belong in version control.
Large dressed/projected numerators and diagnostic DOT output remain generated
workspace evidence, never a precomputed substitute for the slow test.

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
