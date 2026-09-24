= GL638: physical attribution of three retained screen maxima
<gl638-physical-attribution-of-three-retained-screen-maxima>
The three scoped native traces reproduce the previous Arb totals, factors and six cut events exactly. They sample the largest complex norms among retained signed extrema for the selected methods, rather than unrecorded norm maxima over every draw. All physical comparisons below refer to the recorded #strong[identity probe];.

#figure(
  align(center)[#table(
    columns: 3,
    align: (auto,right,auto,),
    table.header([Retained method / seed / channel], [Complete sample magnitude (pb)], [Largest contributions and their geometry],),
    table.hline(),
    [cuts / 22229 / lu\_cut\_5], [1.35512983844], [Cut 1 has norm 1.2323 times the total. The largest individual CT is cut-3 left \[8,12,14\] local (2.1842 times total), with substantial cancellation against other terms. Cut-1 native-2l P/U locals also contribute materially.],
    [cuts\_combined\_joint / 33331 / lu\_cut\_1], [0.250887363126], [Cut 3 has norm 1.1310 times total. Its ordinary right \[5,10\] local and integrated terms have norms 1.5188 and 1.4998 times total. The right star remains far from H/Z: R=284.862 GeV, alpha=0.635661.],
    [cuts\_combined\_joint\_soft / 11113 / soft\_6\_12], [14.3055223644], [Cut 1 has norm 1.0199 times total; its native-2l U local contributes 0.9165 times total. The larger cut-3 U local (0.9952) largely cancels, leaving cut 3 norm 0.03366 times total.],
  )]
  , kind: table
  )

These component norm ratios are #strong[not signed rate fractions] and do not add to one. Returned event and CT weights already include the complete outer Sample factor.

For the large soft-channel maximum, the physical cut-1 point has original `R=sqrt(H²+Z²)=81.9219 GeV`. Its native-2l U star moves to `H=−1.21485`, `Z=−8.37989`, hence `R=8.46749 GeV`. The actual WH denominator radius falls from 63.5949 to 2.04380 GeV; the reconstructed WH multiplier agrees with the measured 0.35332155645116065. Thus this maximum probes materially sharper #strong[star geometry] than its physical base point.

This is not an actual soft-gluon endpoint: the dominant star has soft-13/14 energies 209.664/147.582 GeV, and every physical, single-star and merged-star point reconstructed for these three identity probes has both energies above 43.9 GeV. LU scales span 0.489--1.206 and the threshold radial derivatives remain above 0.443. The retained points do not indicate an extreme LU radial tail or nearly vanishing radial derivative. The composed maximum instead has an ordinary right-threshold origin with substantial local/integrated cancellation.

The joins use measured LU root, native parent/active slots, full centers, alpha/radius, evaluator cut-group/side/local IDs and the semantic variant registry. All 54 single-star joins are unique, with 12 merged residue points retained; maximum affine consistency error is 5.59e−299 GeV. For every reconstructed selected threshold, the sum of its listed positive energies minus Q vanishes within diagnostic tolerance; and 12 actual cut-1 A/U multipliers agree within 9.24e−300. Decimal350 reconstruction preserves native text but supplies no new directed certificate or root solve.

The added joint partition is zero at the composed and composed-soft maxima. This fact alone does not identify a support failure or refute the independently tested supported H/Z approach. A cut-3 projected star can also have H=Z=0 while failing the cut-1 host equation; those are different geometries. The three-point attribution establishes moved-star relevance, not an asymptotic exponent, global boundedness, or a variance cure. The interrupted confirmation run remains separate unresolved evidence.

/* Links use the verified mc_maximum_replay archive layout. */
Reproducible records: #link("gl638_hosted_joint_gate/mc_maximum_replay/traces/analysis/analysis.json")[analysis];, #link("gl638_hosted_joint_gate/mc_maximum_replay/traces/analysis/analyze_max_native_trace.py")[native trace postprocessor];, and #link("gl638_hosted_joint_gate/mc_maximum_replay/traces/analysis/ATTRIBUTION.md")[detailed scope];. The same existing maximum replayer is reserved for the final complete confirmation, after its statistics and program provenance pass their existing gates.
