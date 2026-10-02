= Triage: collective three-component forest cancellation

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.
<triage-collective-three-component-forest-cancellation>
#strong[Current assessment: the captured fixture has an exact factor-preserving cancellation certificate. The failing `collect_factors().is_zero()` predicate does not recognize that cancellation; this is a test normalization regression, not evidence of incorrect forest selection or a nonzero physical soft-counterterm residual. The unchanged repository test remains failed.]

No source or test files were edited for this investigation. The existing compiled test was rerun once with existing logging filters; no rebuild occurred.

== Observed failure
<observed-failure>
Test: `uv::hedge_poset::tests::collective_three_component_region_keeps_only_its_divergent_join`.

Original fresh campaign: `/tmp/soft-ct-validation-2026-09-14/core.log:2238`, 0.131s, ordinary failure. Focused reproduction: `/tmp/soft-ct-forest-trace.log`, 0.14s.

Assertion: `crates/gammalooprs/src/uv/hedge_poset.rs:3633`:

```rust
assert!(
    (actual - expected).collect_factors().is_zero(),
    "the complete collective Taylor value must equal -T_U(1-T_A)(1-T_B)I"
);
```

The failure reports no numerical residual. It says the chosen symbolic simplification did not yield literal zero. It is an exact symbolic comparison, not a numerical tolerance failure.

== Complete-pipeline comparison (before low-level attribution)
<complete-pipeline-comparison-before-low-level-attribution>
Shared inputs: a scalar graph containing two disconnected bubble cycle components linked to a scalar triangle by bridge edges. Bubble A=\[e0,e1\], bubble B=\[e3,e4\], triangle C=\[e6,e7,e8\]. Each bubble has edge-local momentum numerator factors and superficial degree2; C has degree−2. The union U=A∪B∪C has degree2. Both pipelines use the same Graph, model, settings, common physical chart, and ordinary local4D Taylor implementation where called. Integrated counterterms are disabled. This fixture exercises ordinary UV forest closure, not the H=U+S−US soft scheme itself.

#figure(
  align(center)[#table(
    columns: (33.33%, 33.33%, 33.33%),
    align: (auto,auto,auto,),
    table.header([Stage], [Production side (`actual`)], [Explicit oracle (`expected`)],),
    table.hline(),
    [Cycle inventory], [Shared raw graph/spinneys], [Shared raw graph/spinneys],
    [Region selection: first differing boundary], [`classified_spinneys` keeps disconnected regions only when every connected component was independently retained. Triangle C is rejected, so U is rejected.], [Direct `Spinney::new` explicitly constructs U and mixed A+C/B+C without using production classification.],
    [Prefixes], [Retain A, B, their factorized union A∪B, and root.], [Explicitly construct root, signed −T\_A input, signed −T\_B input, and their factorized product (+T\_AT\_B).],
    [Complete collective operations], [No generated operation has cover U. Filtering for cover U sums the empty set, giving0.], [Invoke ordinary local4D `uv_limit` for U with all four prefixes, then sum: −T\_U(1−T\_A)(1−T\_B)I.],
    [Algebra comparison], [Erase momentum provenance, replace denominator wrappers by full polynomials, substitute graph momentum routing, normalize dots.], [Identical normalization.],
    [Final predicate], [Compare difference after `collect_factors`, require literal zero.], [Same predicate.],
  )]
  , kind: table
  )

The production selection is explicit in `uv/uv_graph.rs:165–205`; `Wood::new` calls it at `hedge_poset.rs:174–197`. The test itself verifies all retained components are independently divergent before reaching the failure. The trace confirms production replay only computes its two bubble operations, while manual replay performs the excluded collective Taylor calls.

CFF source reconstruction, residue generation, contour signs, orientation selectors, evaluator precision, and Vakint analytic integration are downstream or unused here; they cannot explain this failure. A classification bug that accidentally retains C has also been excluded by earlier passing assertions. The remaining boundary is the manually constructed collective expression and its normalization/zero test.

== Smallest result matrix already obtained
<smallest-result-matrix-already-obtained>
All these checks ran before the failing assertion or are separately passing in the same fresh core run:

#figure(
  align(center)[#table(
    columns: (50%, 50%),
    align: (auto,auto,),
    table.header([Check], [Outcome],),
    table.hline(),
    [Component count3 and degrees\[-2,2,2\], total degree2], [Pass],
    [Convergent C absent from all retained disconnected UV factors], [Pass],
    [Exactly two independently divergent factors retained], [Pass],
    [A∪B is retained with join factors exactly{A,B}], [Pass],
    [Mixed A+C atomic term nonzero, atomic+nested remainder cancels], [Pass],
    [Mixed B+C atomic term nonzero, atomic+nested remainder cancels], [Pass],
    [Full-U atomic Taylor term nonzero], [Pass],
    [Four-prefix full-U sum compared with omitted production U using collect\_factors], [Fail],
    [Separate two-component bubble+triangle regression], [Pass],
    [Separate disconnected completed soft union regression], [Pass],
  )]
  , kind: table
  )

Assertions after the failing line (the final no-C-cover and surviving-join checks) were not executed; earlier assertions and production source/trace support the same inventory, but the report must not count those later assertions as tested.

== Rebase history
<rebase-history>
`jj diff -r ywptrwuz crates/gammalooprs/src/uv/hedge_poset.rs` shows that change `ywptrwuz` (f61b5f90, \"validate both local UV routes and integrate soft counterterms\") changed two relevant things:

+ Renamed the former `collective_three_component_region_matches_complete_taylor_cancellation` test to the current name and updated it to require omission of U under independent-component selection. Previously production retained collective operations and the test compared their complete value with the explicit four-prefix value.
+ Replaced final `(actual - expected).expand().is_zero()` with `(actual - expected).collect_factors().is_zero()`.

The later `wxtyqpvs` rebase patch only fixes moved-wood usages around this region; it does not repair the final zero comparison. Thus treating the failing expectation as merely an untouched stale API assertion would be inaccurate. The changed symbolic reduction is a concrete suspect, but the historical change also strengthened the omitted-collective cancellation expectation.

== Exact symbolic certificate and residual formula
<exact-symbolic-certificate-and-residual-formula>
A second diagnostic proves the four captured expressions cancel for arbitrary symbolic scalar masses and denominators, without expanding the graph numerator. It uses existing installed SymPy/mpmath libraries and the exact post-Taylor expressions from the unchanged executable.

Artifacts:

- `/tmp/soft-ct-forest-trace/prove-cancellation.py`: reproducible certificate implementation.
- `/tmp/soft-ct-forest-trace/symbolic-certificate.json`: coefficient formula, denominator witnesses, and all four zero certificates.

Define

$ N = Q_0 \( alpha_0 \)^2 Q_3 \( alpha_1 \)^2 \, quad D_A = Q_0^2 - u^2 \, med D_B = Q_3^2 - u^2 \, med D_C = Q_6^2 - u^2 \, $

where $u = m_(U V)$, $e = m_(U V e x p)$, and $m$ is the physical scalar mass. The four logged raw full-U prefixes, ordered root, A, B, A∪B, are exactly

$ \( F \, - F \, - F \, F \) \, #h(2em) F = frac(i N [- D_A D_B D_C + \( e^2 - m^2 \) \( 3 D_A D_B + 2 D_A D_C + 2 D_B D_C \)], D_A^3 D_B^3 D_C^4) . $

The subsequent common subtraction minus in `uv_limit` changes their signs to $\( - F \, F \, F \, - F \)$, leaving their sum zero. Thus the exact physical residual is

$ upright(a c t u a l) - upright(e x p e c t e d) = 0 - \( - F + F + F - F \) = 0 . $

Proof boundaries:

+ Erase each momentum-provenance wrapper to its stored literal carrier, preserving signs.
+ Apply this vacuum graph\'s exact common chart: Q1=−Q0, Q4=−Q3, Q7=Q8=Q6, Q2=Q5=0. Every captured full-U call uses compatible coordinates\[e0,e3,e6\].
+ Check each physical denominator owner against its complete polynomial: owners0/1 match D\_A,3/4 match D\_B,6/7/8 match D\_C. Keep all denominator powers.
+ Represent graph numerator factors by a separate exponent witness. Arithmetic adds terms only when their numerator powers are identical; zero terms are discarded, products add exponent witnesses, and numerator powers are forbidden in denominators. Every complete prefix has exactly powers\[2,2\], so N remains intact throughout.
+ Apply `together`/`cancel` only to the remaining scalar rational coefficient, with D\_A,D\_B,D\_C,m,u,e symbolic. The script certifies root+A=0, root+B=0, root−AB=0, and all-four-sum=0 exactly.

This is a symbolic certificate for this captured fixture, not a claim about arbitrary unrelated forests. It is stronger than the coordinate samples below. It confirms the physical zero expected by the test, while showing the repository\'s current normalization predicate is insufficient. The test was neither edited nor marked passing.

== Additional exact numerical evidence from unchanged production expressions
<additional-exact-numerical-evidence-from-unchanged-production-expressions>
Focused command (existing binary):

```bash
GL_DISPLAY_FILTER=off \
GL_LOGFILE_FILTER='off,[{#uv,#integrated,#series}]=debug,[{#generation,#uv,#four_d,#trace}]=debug' \
GL_TEST_LOG_DIR=/tmp/soft-ct-forest-trace \
target/dev-optim/deps/gammalooprs-00fead27b5c7d936 \
  --exact uv::hedge_poset::tests::collective_three_component_region_keeps_only_its_divergent_join \
  --nocapture
```

Durable trace: `/tmp/soft-ct-forest-trace/gammaloop-test-3622015.jsonl`. Extracted post-Taylor, t=1 expressions:

- `collective-prefix-0.txt` (root)
- `collective-prefix-FU.txt` (B)
- `collective-prefix-F.txt` (A)
- `collective-prefix-Fj.txt` (A∪B)

Every full-U call uses loop coordinates\[e0,e3,e6\], with crown edges\[e2,e5\]. A rational evaluator applied the existing test\'s physical-chart normalization directly to the logged scalar expressions: provenance resolves to its stored momentum payload; Q1=−Q0, Q4=−Q3, Q7=Q8=Q6; bridge momenta Q2=Q5=0; Minkowski signature(+−−−); denominator wrappers evaluate to their fourth argument. The five function types appearing in the expressions are Q, momentum\_provenance, mink, dot, and denom. Arithmetic used Python Fraction throughout. The expressions are purely imaginary, so the stored data report their imaginary coefficients; uv\_limit\'s subsequent common overall minus does not affect cancellation.

Three distinct integer four-momentum triples were tested, each with two independent mass choices (including mUVexp≠mUV):

#figure(
  align(center)[#table(
    columns: (25%, 25%, 25%, 25%),
    align: (auto,auto,auto,auto,),
    table.header([(Q0,Q3,Q6)], [(physical mass,mUV,mUVexp)], [Four-prefix imaginary coefficients], [Sum],),
    table.hline(),
    [((11,2,3,5),(13,7,2,1),(17,2,5,7))], [(2,3,5)], [(+X,−X,−X,+X), X=647452729/200889964234905364736], [Exactly0],
    [Same], [(5,7,7)], [(+X,−X,−X,+X), X=130585/297775864704384], [Exactly0],
    [((23,3,7,9),(19,2,3,5),(29,7,4,3))], [Both mass choices], [(X,−X,−X,X), X nonzero], [Exactly0],
    [((31,5,2,1),(37,3,4,2),(41,6,9,7))], [Both mass choices], [(X,−X,−X,X), X nonzero], [Exactly0],
  )]
  , kind: table
  )

Full exact values: `/tmp/soft-ct-forest-trace/rational-evaluation.json`.

These six exact-coordinate checks support the independent symbolic certificate above. On their own they would be additional checks, #strong[not a universal symbolic proof];. Neither the samples nor the standalone certificate change the outcome of the unchanged repository test.

== Remaining work
<remaining-work>
The next required action is to restore a factor-preserving exact zero comparison in the existing test while retaining the same expected physical result and assertions. The adjacent disconnected replay test already uses `together().cancel()` at `hedge_poset.rs:2318`, which is a candidate existing approach; its suitability must be checked against the graph-numerator-factorization requirement. The standalone certificate demonstrates the intended mathematics without changing production behavior.

No expected value, test, production code, or comment was changed. Repository guidance requires user confirmation before changing a failing test. Until an authorized correction is made and the repository test reruns successfully, report this as #strong[failed test, independently certified physical cancellation, normalization repair pending];.
