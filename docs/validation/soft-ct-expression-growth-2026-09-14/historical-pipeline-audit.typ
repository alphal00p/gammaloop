= Historical direct Taylor comparison: c08b0d9f → efaf2da1

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.
<historical-direct-taylor-comparison-c08b0d9f--efaf2da1>
#strong[Ordinary U already used hard-loop rescaling and a Laurent series before this stack.] The removed external-momentum Taylor path and `OSE_FOR_LOCAL_3D_SERIES` derivative helper belonged to the old #strong[soft S (`t_tilde`) path];, not ordinary U. Their removal cannot directly explain a slowdown in a pure-MUV run. Current U nevertheless has changed routing and OSE/mass representation, so a same-revision U/H comparison does not exclude a historical regression affecting both.

Exact read-only snapshots:

- Old kernel;, commit `c08b0d9f1605f7e3ec7445508c69edc07a4b599f`.
- Measured kernel;, commit `efaf2da158c5b01bddbc5435592e7cfea7a45a91`.

No rerun, generation, code edits or numerator expansion was performed for this comparison.

#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,auto,auto,auto,),
    table.header([Boundary], [Old ordinary U], [Measured ordinary U], [Potential size implication],),
    table.hline(),
    [Coordinate chart], [Starts from `current.lmb()`, retains loop carriers by edge incidence; fallback searches compatible bases. Old lines 75--121;], [Starts from the graph-canonical induced component LMB; retains whole prior frames and rejects affine retained carriers; returns canonical and selected operation charts. New lines 63--119;], [Earlier than Taylor algebra: a different routed polynomial or chart may change expression size. Actual selected LMBs must be compared before blaming series internals.],
    [Common start], [Multiplies CFF by mapped numerator, opens active OSE arguments, splits temporal/spatial momentum. Old start;], [Same sequence, with physical masses tagged by component scope. New start;], [Scope prevents premature identification of masses from different completed components; may retain distinctions longer. Tags are unit-valued and removed in final assembly, not before nested Taylor.],
    [Loop rescaling], [Rescales every selected loop; fixed external shift comes from `current.lmb().ext_atom(...)`. Old rescale;], [Same hard-loop formula, but shift comes from the selected passed `lmb.ext_atom(...)`. New rescale;], [This is a real U input change when charts differ, not replacement of external Taylor by hard Taylor.],
    [OSE rearrangement], [Every active OSE gets `(m_exp² λ² + P − m_exp²)/λ²`, with outside λ². Old OSE;], [Fresh OSE gets `(m_vac² scope² λ² + P − m_exp² scope²)/λ²`. An existing vacuum OSE instead co-scales its actual vacuum mass and divides by λ², avoiding a second fresh rearrangement. New OSE;], [U\'s actual symbolic input and opportunities for cancellation changed. Distinct vacuum/expansion heads and owner tags can affect intermediate size; no size factor was measured. The existing-vacuum guard is intended to reduce redundant rearrangement in nested U.],
    [Series/materialization], [Multiply by λ^(3L), invert λ, `series(λ,0,0)`, `to_atom`, set λ=1. Old series;], [Same algorithm and endpoint for ordinary U. New series;], [`series/to_atom` was not newly introduced. Any historical difference here must be explained by different inputs, Symbolica behavior, or surrounding retained work.],
  )]
  , kind: table
  )

The substantive #strong[soft] rewrite is separate. Old `t_tilde` scales external momenta in the newly added numerator, supplies routed external shifts, replaces OSE by the special helper whose derivative is one for argument index two and zero otherwise, takes `series(...,-1)`, then rewrites helper derivatives and restores OSE. Old helper;, old S path;, helper use/derivative cleanup;.

Measured S uses the complete common-start CFF, co-scales component physical masses with hard loops, factors its physical OSE without inserting MUV, and takes the same hard-dual Laurent endpoint −1. Measured H uses `S(X)+U(X−S(X))`; old H separately computed U(X), S(X), and U(S(X)). Thus the stack removed one explicit overlap materialization while changing the S input representation substantially. New S mass scaling;, new S OSE;, new H;, old H;.

Neither old nor measured ordinary U explicitly expands graph numerators in this code. Both materialize Taylor coefficients, so that source-level observation does not certify factorization or bounded size after the series. Historical speed alone also does not certify equivalent subtraction: routing, mass ownership and the S derivative contract changed.

The smallest historical comparison adds #strong[old U] to the matched #strong[measured U / measured H] matrix on the exact same graph, retained forest, cut and one production residue key. Compare selected charts first; then byte sizes and factorized OSE/numerator structure immediately before and after series. If old/current U already diverge in size, changing only the soft prescription cannot isolate that regression. If U stays bounded while H grows, the first S output or U-on-remainder input is the next useful reproducer. These are measurement proposals, not completed historical timing results.

The historical evaluator diff contains only three test API spelling updates (`emr_vec_index`→`emr_vec`), with #strong[no production constructor or preprocessing changes];. Selector scanning and generic evaluator alias restoration therefore predate this stack; the observed scan hotspot can be triggered by newly large inputs without being a newly introduced evaluator loop. The Spenso network diff changes odd tensor-power execution and comments; it does not alter scalar alias restoration. This rules out attributing those specific pre-existing passes to a new implementation change in these files, but does not exclude changed inputs reaching them or changes in other preprocessing components. Evidence: exact evaluator diff;, exact Spenso network diff;.
