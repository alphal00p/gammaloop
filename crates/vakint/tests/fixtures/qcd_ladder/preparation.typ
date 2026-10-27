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

The ignored #raw("slow::") acceptance test and the collector share the same native
preparation helper. Prepare the complete scalar input once, without applying IBPs:

```sh
cargo run --locked --release -p vakint \
  --features "feynkit-ingress experimental-rustred \
weighted-lowering-diagnostics symbolica_community_module" \
  --example five_loop_qcd_ladder -- \
  --work-dir /path/to/workspace/quadratic-input --fresh-prepare --prepare-only
```

This saves #raw("prepared-input.json"), the exact #raw("vakint-input.symbolica") and
all projector conditions. Subsequent invocations reuse a valid bundle in their
#raw("--work-dir") automatically. An explicit #raw("--prepared-input PATH") selects
that directory for a separate fresh-process reduction, without graph generation,
Dirac/color contraction, Taylor expansion or angular projection:

```sh
five_loop_qcd_ladder --work-dir /path/to/workspace/uncached-run \
  --prepared-input /path/to/workspace/quadratic-input \
  --artifacts /path/to/vakint-external-data
```

The loader checks the input and condition hashes, projection, complete contribution
inventory, scientific source definitions and producer provenance. A corrupt or
incompatible existing bundle, including payloads with a missing manifest, is an
error before native preparation, never a reason to regenerate silently.
#raw("--fresh-prepare") explicitly repeats native validation in a new work directory;
it cannot overwrite a completed bundle or combine with #raw("--prepared-input").
#raw("--prepare-only") creates or validates the selected bundle and exits before
IBPs. #raw("--dress-only") performs native graph/rule preparation and validates the
vacuum incidence without gamma algebra; it cannot consume a prepared bundle. The default
#raw("--projection quadratic") preserves the full transverse coefficient above.
The additional #raw("--projection zero-momentum") diagnostic retains both current
vertices and contracts with unaveraged #raw("g_munu*delta_ab") at zero external
momentum. Its unlowered source has #raw("A=14, r_mom=5, D_source=9") and superficial
#raw("omega_UV=2") at four dimensions. It does not use a transverse division,
angular projector, or polarization average.
It is an additional diagnostic of the bare graph, rather than the full soft
mass term of the recursively subtracted counterterm.

The collector records native preparation stages, data loading and cache import,
geometric routing, weighted scalar lowering, the combined native application,
subsequent terminal geometry and lower-loop factors, and total elapsed time. Its input and result hashes, executable hash, source revision,
data identity, exact conditions, and cache counters are written alongside the
expressions in #raw("benchmark.json"). #raw("--cache PATH") explicitly attaches a
prepared five-loop cache; #raw("--warm") performs a second evaluation and requires
exact result, guard, and family parity. Repeated uncached/cached measurements use
separate fresh processes consuming the same verified scalar-input bundle.
Native preparation is measured separately. For reused input, #raw("input_mode") is
#raw("prepared"), current #raw("preparation") and #raw("preparation_seconds") are null,
and #raw("input_loading_seconds") measures verification and parsing. The manifest
retains the original producer's preparation timings and source hashes as historical
provenance; they are excluded from the current end-to-end duration. A fresh
#raw("--fresh-prepare") run reports its actual native preparation separately.
Uncached runs explicitly detach any backing attached by the data loader;
cached runs clear resident results while retaining the selected backing.
The receipt measures deferred owner preparation from actual cumulative counters
at the data-loading, cache-setup, and reduction boundaries, and records persistent
root and intermediate hits. The combined native application includes entry
transport, lazy owner preparation, IBP application, exact summation, and native
terminal substitution. The current native API does not expose disjoint IBP and
terminal-substitution durations; these fields remain null, with no estimate by
subtracting overlapping wall times. Execution work is governed by time and RAM
availability, rather than arbitrary term, operation, routing-call, or rule counts.
The quadratic auxiliary-rank check is a scientific invariant and remains in place.
Historical sealed runs retain their original work ceilings and failure receipts:
an earlier indexed adapter reached its then-default 100,000-term allowance,
and later validation used an explicit one-million-term allowance. Those receipts
remain provenance; neither an execution ceiling nor its failure establishes IBP
coverage or surviving terminals. The scalar-work override interface is removed.

The slow test requires the explicit #raw("VAKINT_SAVED_APPLICATION") and
#raw("VAKINT_QCD_LADDER_EVIDENCE") paths. It automatically reuses the verified bundle
in the evidence directory after the first successful preparation; optionally set
#raw("VAKINT_QCD_LADDER_PREPARED_INPUT") to a separate saved-input directory. Set
#raw("VAKINT_QCD_LADDER_FRESH_PREPARE=1") and choose a new evidence directory for
explicit fresh native validation. The complete census and exact cold/warm result,
family and guard checks run for either input mode. Its preparation census distinguishes a
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
Completed checkpoints are addressed by the exact source identity,
bound to the input and application hashes and family fingerprint, and checked
again against the active admission predicates. Imported historical checkpoints
retain their original execution-policy receipt; they do not impose a new work
ceiling or change the current time-and-RAM policy. Their local contribution ordinals
are presentation metadata and are reassigned for the current enumeration.
The default census requires all 21 quadratic groups and checks engineering
mass-dimension balance #raw("D_key-j=10"). For the additional diagnostic, explicitly selecting
#raw("VAKINT_QCD_CENSUS_PROJECTION=zero-momentum") requires one group and
#raw("D_key-j=9"); other or mixed projection values are rejected. These checks
do not certify individual-term UV degree or campaign coverage.

The compact source fixture and reproducible code remain the scientific definition.
Large dressed/projected numerators and diagnostic DOT output remain generated
workspace evidence. A verified prepared scalar-input bundle is a reusable input to
the IBP acceptance test; it does not replace the separate fresh-preparation check
or contain a precomputed reduction result.

The separately ignored prepared-input diagnostic records the native execution
policy and still checks all original visible input poles, mass homogeneity, and
the actual saved entry union without applying IBPs. Historical budget failures
remain sealed evidence, rather than physical-coverage failures or successful
fresh-preparation tests. No operation-count environment override remains.
Contribution ordinals are local to each run: compare the recorded family, key,
and exact numerator rather than assuming that equal ordinals identify equal
inputs. Its per-contribution histogram is also not a globally coalesced sum
or proof that every equivalent momentum presentation is unsupported.

== Saved zero-momentum coefficients

#link("../../../docs/five_loop_q0_coefficients.typ")[The saved-result companion]
loads the complete retained zero-momentum reduction and its 300-digit scalar
master-product coefficients without repeating IBPs. The derived scalar jets
include $epsilon^5$ with coefficientwise remainder $cal(O)(epsilon^6)$;
all master integrals remain formal.
