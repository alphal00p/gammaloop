= Three retained maximum samples: native CT attribution
<three-retained-maximum-samples-native-ct-attribution>
All three traced Arb evaluations exactly reproduce the previous native totals, factors and six cut events. The saved-state and source-input hashes remain unchanged. These are maxima among retained signed-extremum Samples, not unrecorded global complex-norm maxima.

The reconstruction uses the actual first (identity) probe\'s six measured LU roots, original routing and native parent bases. Each center is joined through the explicit evaluator\'s cut-group/side/local ID, the solve-group variant registry, matching parent/active slots, and the measured fixed complement/radius/alpha. All 54 single-star joins are unique. For every reconstructed selected threshold, the sum of its listed positive energies minus Q vanishes within diagnostic tolerance. The six physical cut shifts are also structurally serialized and checked; the other selected threshold shifts are not a separate structural certificate. Twelve merged left/right residue points are also retained. Maximum fixed-complement/radius/alpha discrepancies are 5.59e−299 GeV. Twelve independently reconstructed cut-1 native-2l A/U WH multipliers agree with the measured multipliers within 9.24e−300.

The computations use Decimal350 on the preserved native decimal text. They do not invoke Gamma, solve another root, or supply a new directed numerical certificate. `analysis.json` retains full-precision coordinates, root/alpha values, semantic component records, trace indices and hashes. `analyze_max_native_trace.py` reproduces it from the immutable trace directory.

#figure(
  align(center)[#table(
    columns: 3,
    align: (auto,right,auto,),
    table.header([Method / seed / selected channel], [Complete sample magnitude (pb)], [Measured attribution],),
    table.hline(),
    [cuts / 22229 / lu\_cut\_5], [1.35512983844], [Cut 1 supplies 1.2323 times the total norm. Several CTs cancel: the largest is cut-3 left \[8,12,14\] local (2.1842 times total), followed by its left/right local-local combination (1.6116), and cut-1 native-2l \[3,12\] and \[8,12,14\] locals (1.4478 and 1.4097).],
    [cuts\_combined\_joint / 33331 / lu\_cut\_1], [0.250887363126], [Cut 3 supplies 1.1310 times total. Its right \[5,10\] local and integrated terms have norms 1.5188 and 1.4998 times total. This is an ordinary projected-threshold contribution, with substantial cancellation, at a point far from the direct H/Z corner.],
    [cuts\_combined\_joint\_soft / 11113 / soft\_6\_12], [14.3055223644], [Cut 1 supplies 1.0199 times total. Its native-2l \[8,12,14\] local term supplies 0.9165 times total. Cut 3 has a larger individual \[8,12,14\] local term (0.9952), but its net contribution is only 0.03366 times total after cancellation.],
  )]
  , kind: table
  )

Component norms are compared with the complex total norm and do not add as signed fractions. The native decomposition already contains the full outer Sample weight; no grid, Jacobian or partition factor was applied again.

== Actual physical and star distances
<actual-physical-and-star-distances>
All distances and soft energies below are in GeV. H and Z are the original routed rows \[2,4,12\] and \[3,10,13\]. The WH radius is `sqrt((H−eta_host)^2+(P Z/Q)^2)` and is identified with the actual multiplier only for the checked cut-1 native-2l A/U consumers.

#figure(
  align(center)[#table(
    columns: 7,
    align: (auto,right,auto,right,right,right,right,),
    table.header([Case], [Cut-1 physical R(H,Z)], [Relevant CT star], [Alpha], [Star R(H,Z)], [Star WH radius], [Star soft-13 / soft-14],),
    table.hline(),
    [cuts / 22229], [302.4515], [cut-1 native-2l \[8,12,14\]], [0.7956642], [39.0003], [19.7776], [139.973 / 237.647],
    [composed / 33331], [273.7441], [cut-3 right \[5,10\]], [0.6356612], [284.8623], [not that consumer], [330.804 / 347.071],
    [composed-soft / 11113], [81.92188], [cut-1 native-2l \[8,12,14\]], [1.1140280], [8.467493], [2.043800], [209.664 / 147.582],
  )]
  , kind: table
  )

For the large soft-channel sample, the physical cut-1 WH radius is 63.5949 GeV. Its A/U star moves to H=−1.21485 and Z=−8.37989 GeV, and the measured WH multiplier is 0.35332155645116065. Thus the base point does not describe the sharper normal geometry sampled by that dominant CT. This is concrete evidence that moved-star geometry matters for this maximum; it does not establish an asymptotic divergence or prove that a particular new proposal would cure its variance.

The soft-channel label does not identify an actual soft-gluon limit here. Across all six physical LU points and all reconstructed single and merged stars, the soft energies stay at ordinary nonzero scales; the single-star minimum soft-14 energy for the large soft-channel sample is 43.9258 GeV. Its dominant cut-1 star has soft energies about 210 and 148 GeV. Physical LU scales lie between 0.489 and 1.206 across these three cases, with derivatives at least 636.06 GeV; threshold radial derivatives are at least 0.443. These observations do not indicate an extreme radial tail or nearly vanishing radial derivative at the retained points.

Some cut-3 right stars satisfy both the original H row and Z row to reconstruction precision: H vanishes because cut 3 is on H, while the right projection sets Z. Their cut-1 host residual is not zero. They must not be labeled the same hosted cut-1 H/Z corner merely from those two zeros.

The added joint channel has exactly zero partition at the composed and composed-soft maxima. That fact alone proves neither failed compact support nor a boundedness defect. The retained canonical map supplies the complete existing density/partition factors; this trace did not separately query the joint\'s directed domain certificate. The narrower supported H/Z approach tests remain separate evidence.

== Reproduction and limits
<reproduction-and-limits>
Run only the postprocessor, after the trace\'s exact controls have passed:

```bash
source /tmp/gammaloop-rebase-dev-env.sh
python /tmp/gl638-hosted-joint-gate/analyze_max_native_trace.py \
  /tmp/gl638-hosted-joint-gate/results-max-native-trace-build5 \
  /tmp/gl638-hosted-joint-gate/NEW-max-native-trace-analysis
```

Single-star roots, the common merged residue point and LI/IL/II multiplier contexts are distinct. The report retains the common merged point and each single root; it does not relabel every integrated component as an LL star. The three-point attribution cannot certify global boundedness, the tail distribution or MC convergence. The separate interrupted confirmation run remains unresolved by these successful maximum replays.
