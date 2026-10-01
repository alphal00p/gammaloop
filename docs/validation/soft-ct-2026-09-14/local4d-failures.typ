= Local 4D soft-counterterm failure triage

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.
<local-4d-soft-counterterm-failure-triage>
Input revision: `efaf2da158c5b01bddbc5435592e7cfea7a45a91`. Read-only triage of the four local-4D assertion failures in core.log and core.junit.xml;. No repository files, tests or assertions were changed; no Rust tests were rerun by this triage.

== Assessment
<assessment>
The four observed failures separate into two comparison-boundary problems:

+ #strong[Three production-versus-paper comparisons normalize color differently.] Production calls `.color_simplify()` before the local Taylor operation, while each independent oracle keeps raw color generators. The comparison helpers simplify metrics and scalar rational structure, but never perform color reduction. The first unequal expressions retain explicit color-identity residuals. The free-fermion case fails before any UV rescaling, which directly excludes mass grading as the cause of its current failure. This shared boundary is established from source and logs; later assertions in those tests have not been exercised successfully.
+ #strong[The top-bubble test mistakes a factorized zero for a nonzero soft coefficient.] Both full printed coefficient residuals cancel exactly using common-factor extraction and regrouping, preserving every tensor factor and original denominator owner. This is an exact identity check on the captured expressions, not just a numeric substitution. It establishes that the reported constant and linear residuals are zero, but does not convert the still-failing Rust test into a pass.

These findings do not justify changing production UV algebra, CFF recursion, residue signs, or the independent physics expectations. The tests remain failing until their actual comparison boundaries are resolved and rerun. Any proposed edits to failing tests require the user\'s confirmation under repository guidance.

== Pipeline comparison before deeper diagnosis
<pipeline-comparison-before-deeper-diagnosis>
#figure(
  align(center)[#table(
    columns: 5,
    align: (auto,auto,auto,auto,auto,),
    table.header([Cases], [Last shared input], [Production route], [Oracle route], [First differing boundary],),
    table.hline(),
    [B.7 massive outer H2; B.11 massless child; free-fermion mass grading], [Same parsed graph, component, numerator and d-dimensional conversion], [numerator → color simplification → denominator division → metric simplification; then U/S as applicable], [numerator → raw atom without color simplification → denominator/independent direct derivatives], [Color representation immediately after d-dimensional numerator conversion],
    [Actual top-bubble H2], [Same grown child atom and compatible LMB, local-only IR, d=2], [Build U, S, US and H; production H agrees with explicit U+S−US], [Erase momentum provenance, route, scale external crown, extract constant/linear coefficients, finalize metrics, collect factors, structural is\_zero], [Final algebraic zero detection on an undistributed factorized coefficient],
  )]
  , kind: table
  )

All four cases are purely local 4D at their failure points. No CFF generation, selector assignment, contour integration, evaluator precision or profile fitting is involved. Neither analytic integration nor backend topology conversion can cause these specific assertions.

== Smallest result matrix
<smallest-result-matrix>
#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,auto,auto,auto,),
    table.header([Test], [Mode / first assertion reached], [Result / elapsed], [What it narrows],),
    table.hline(),
    [production\_outer\_u\_grades\_a\_free\_fermion\_mass\_but\_keeps\_muv\_fixed], [Massive one-loop fermion self-energy, MUV, integrated=false; compare grow with raw N/D before U], [Fail, 0.066 s], [Color mismatch exists before UV deformation or mass grading],
    [paper\_appendix\_b\_10\_b\_12\_massless\_child\_and\_local\_b\_16\_nested\_branch], [Massless child H1, integrated=false; compare physical S0 with B.11], [Fail, 0.069 s], [Same mismatch persists with zero physical fermion mass; nested B.16 assertions are not reached],
    [paper\_appendix\_b\_6\_b\_7\_production\_outer\_h2\_matches\_the\_routed\_rhs], [Massive outer H2, integrated=false; compare physical S1 with B.7], [Fail, 0.089 s], [Failure is already in the physical branch; UV-quadratic and completed-H comparisons are not reached],
    [top\_bubble\_actual\_fermion\_child\_h2\_removes\_both\_soft\_jet\_coefficients], [Actual massive fermion bubble H2, integrated=false; independent q\_g^0/q\_g^1 coefficients], [Fail, 0.087 s; both captured residuals exactly cancel under factor-preserving normalization], [Final zero comparison, with no remaining evidence of a nonzero soft coefficient in this captured fixture],
  )]
  , kind: table
  )

Already-passing neighboring checks include the scalar H/refinement identities, paper Eqs. 4.15/4.16/4.19, early color simplification across open/nested boundaries, denominator/provenance classes, and immutable momentum-provenance erasure. Those passes support examining the shallow representation difference first; they are not a substitute for the failed physical fixtures.

== Three color-normalization mismatches
<three-color-normalization-mismatches>
Production #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:1805")[grow\_factor] obtains the numerator with

```text
numerator(...).to_d_dim(GS.dim).color_simplify().get_single_atom()
```

The paper B.7 oracle at #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:6237")[line 6237];, B.11 oracle at #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:6577")[line 6577];, and fermion oracle at #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:7194")[line 7194] omit that color-simplification step.

#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:2321")[finalize] only calls `simplify_metrics()`. Neither #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:5749")[paper\_assert\_normalized\_zero] nor #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:6009")[paper\_assert\_atomic\_rhs] closes color generators. Thus raw generator products and their already-reduced Casimir/trace values survive as different Symbolica atoms.

=== Massless B.11 child
<massless-b11-child>
The #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:6710")[first assertion] reports `Appendix B.11 massless physical tensor/gamma term`, at the shared helper line 6029. The displayed difference is, suppressing unchanged explicit tensor slots,

```text
GC_11^2 * (-t^a_ij t^a_jk + delta_ik C_F) * Q gamma gamma gamma / D_k^2
```

See core.log:586;. The exact color factor is the raw fundamental generator product minus its quadratic-Casimir reduction. It multiplies a common momentum/gamma/denominator factor; no Dirac contraction is necessary to locate the mismatch.

Before this assertion, direct checks of physical and UV denominator routing, their p=0 values, and their linear shifts completed successfully (#link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:6652")[lines 6652 onward];). The test never reaches its outer nested B.16 acceptance, so the failed result must not be described as a demonstrated nested-B.16 algebra error.

=== Massive B.7 outer H2
<massive-b7-outer-h2>
The #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:6456")[first failed assertion] reports `Appendix B.7 physical tensor/gamma lines`, again at helper line 6029. Its residual contains the unreduced four-generator color chain on one side and `g(color0,color1) C_F T_R` on the other, attached to the same factorized fermion/gamma structure.

Independent scalar-denominator/routing checks precede and pass this assertion. The mismatch occurs in physical S1, before the UV-quadratic comparison at line 6462 or full-H comparison at line 6468. Early color closure is the first required comparison; changing Taylor signs or UV mass grading would not address that boundary.

=== Free fermion mass
<free-fermion-mass>
The #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:7210")[precondition] compares production grow with raw numerator/denominator and reports:

```text
production grow uses the numerator whose fermion-mass term is tested
```

The difference contains the same `-t^a t^a + delta C_F` factor multiplying `(MT*g + Q*gamma)`, the outer gamma factors and propagators. The actual UV rescaling begins only at #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:7218")[line 7218];. This failure therefore says nothing yet about where the free fermion mass enters the Laurent coefficients.

#strong[Next comparison for this cluster:] put both sides at the same color-normalized numerator boundary using the existing color simplifier, preserving open tensor boundaries and factorization, then compare the actual full physical branch. This is a proposed diagnostic comparison, not an assertion that the later unexecuted tests already pass.

== Actual top-bubble residual: exact zero after factoring
<actual-top-bubble-residual-exact-zero-after-factoring>
The fixture compares the two independent soft coefficients at #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3385")[lines 3385--3388];. The earlier explicit H=U+S−US identity and the nonzero ordinary-U constant-remainder control have passed. The coefficient extractor ends with

```text
finalize(coefficient.replace(dim).with(4)).collect_factors()
```

at #link("../../../crates/gammalooprs/src/uv/approx/local_4d.rs:3375")[line 3375];, then uses structural `Atom::is_zero`.

The failure message contains the complete q\_g^0 and q\_g^1 expressions, respectively 8,330 and 8,012 characters. I parsed those captured expressions into an exact arithmetic AST, replacing each full indexed tensor/function atom with a unique opaque commuting symbol. Original `denom(7,...)` and `denom(8,...)` stayed distinct; no denominator classes, momentum slots, masses or physical color factors were identified.

Using the installed Sympy 1.14 `factor_terms(clear=False, fraction=False)` routine, whose implementation was inspected and performs common-factor extraction without expansion:

- The linear coefficient reduces to exactly zero under repeated factor extraction.
- The constant coefficient reduces to a common prefactor times

```text
-A*N + B*N + A*(N - B*G) - B*(N - A*G).
```

Here N and G remain complete product factors; A and B are the distinct denominator atoms. Grouping the first/third and second/fourth terms produces `-A*B*G + A*B*G = 0`. The diagnostic implements exactly this common-factor regrouping, with a fixed iteration bound and no numerator expansion.

Exact diagnostic output ends with:

```text
coefficient 0: exact result after factor-preserving pair grouping = 0
coefficient 1: exact result after factor-preserving pair grouping = 0
```

The reproducible script runs with the existing Python executable `/nix/store/iizfzmqjdfvbyx82gqw92glb5zl6m4my-python3-3.12.12/bin/python3` and cached Sympy/mpmath package paths. It does not import repository code or rerun Rust tests. Its preliminary exact-rational substitutions are supplementary only; the final factorization identity is the evidence for cancellation.

This excludes a nonzero constant/linear coefficient, missing gamma-trace contraction, and denominator-owner identification as explanations of this specific captured failure. The next comparison belongs in the existing factor-preserving symbolic normalization/zero-check boundary, not in soft routing or CFF reconstruction.

== Reporting limits
<reporting-limits>
The actual runner still records all four cases as failures. No pass count should be adjusted based on triage. The three color cases establish a first representation mismatch but have unexecuted downstream assertions; the top-bubble result proves cancellation only of the captured fixture coefficients. The current physical CLI acceptance and full rebase validation must still be reported from their own executed results.
