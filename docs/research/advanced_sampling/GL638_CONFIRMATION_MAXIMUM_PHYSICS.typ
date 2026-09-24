= GL638 confirmation maximum physics
<gl638-confirmation-maximum-physics>
The retained confirmation maxima have two sources of large weights: counterterms evaluated near their own threshold intersections, and large opposing terms on other physical cuts. None of the three investigated unique b/c extrema approaches a soft-gluon endpoint. This is pointwise evidence; it does not establish boundedness or finite variance elsewhere.

The fixed confirmation remains a mixed result: adding composed Cut 1/HZ to the six cuts gives real variance 1.10754 times cuts-only and imaginary variance 0.64931 times cuts-only. Its real maximum is larger. The preselected method and all reported samples remain unchanged; there is no post-hoc winner selection.

== Exact explanation of the shared real maximum
<exact-explanation-of-the-shared-real-maximum>
Cuts-only Re maximum 10.4734102782pb and joint-added Re maximum 12.2189786579pb, both seed 40009 / lu\_cut\_1, have exactly the same canonical 12-coordinate point, six shared partition weights, raw density sum, selected J\*w and pre-outer native integrand. The added joint partition is exactly zero. Their actual complete Sample weights are 6 and 7, so the final complex estimator changes by 7/6, with relative decimal discrepancy below 7.27e-302. The joint real maximum supplies 48.17% of its empirical second moment, versus 39.20% for cuts-only.

This explains that single maximum\'s increase through actual sampling probabilities. It neither explains the whole variance ratio nor proves an invalid proposal. The event and CT weights already contain the outer factor.

== Measured physical and counterterm geometry
<measured-physical-and-counterterm-geometry>
Here b denotes six cut channels, and c adds the composed joint channel. H and Z denote the signed energy surfaces `[2,4,12]` and `[3,10,13]`; the displayed radius is `R = sqrt(H^2 + Z^2)`. The actual multiplier is `W_H = H^2 / (H^2 + (P Z/Q)^2)`, so its denominator radius, quoted separately below, is `sqrt(H^2 + (P Z/Q)^2)`. The existing metadata defines P and Q.

#figure(
  align(center)[#table(
    columns: 5,
    align: (auto,right,auto,right,right,),
    table.header([Actual maximum], [Physical cut1 R \[GeV\]], [Main measured CT configuration], [Star R \[GeV\]], [Star gluon13/14 energies \[GeV\]],),
    table.hline(),
    [Shared b/c Re, seed40009], [125.907], [Cut3 right\[5,10\], parent\[4,5,7,12\], active\[5\]], [51.1309], [89.217/349.104],
    [b Im, seed40009/lu\_cut\_5], [29.0311], [Cut1 native2l A=\[7,8\], parent\[3,6,7,10\], active\[3,7\]], [1.57488], [185.803/478.966],
    [c absolute-Im, seed20011/lu\_cut\_0], [23.9737], [Cut1 native2l U=\[8,12,14\], same native parent/active slots], [12.1308], [187.177/271.100],
  )]
  , kind: table
  )

For the shared Re point, the c estimator is 12.21898+2.63825i pb. The largest terms are cut3 right\[5,10\] local−42.16297i and integrated+30.18901, with substantial opposing other counterterm contributions. Net cut3 is 12.17651−1.35361i. The dominant right-star has cut1-host defect90.1152GeV: its H=0 follows from physical cut3, not the hosted cut1 condition. This is another cut/threshold configuration with large cancellations; it must not be relabeled as failure of the cut1-hosted corner chart.

For b Im, total 1.76318+11.51988i is almost entirely net cut1. The dominant native2l A local term is+19.26394i, opposed by native2l Z=\[3,10,13\] local−8.91995i. At the A-star, the actual WH radius decreases21.6219→1.19117GeV. At the opposing Z-star, Z vanishes within diagnostic precision and H=−0.0155391GeV, with hard gluons186.080/480.046GeV. This is measured enhancement near a moved counterterm-star intersection. Here Z names\[3,10,13\]; the historical P surface is\[3,12\]. The Z-star radius is geometry only; no A/U multiplier is assigned to it.

For c Im, total 0.791036−6.226383i has net cut1=0.171433−5.795819i. Native2l U local−4.595399i dominates, with native2l\[3,7,14\]+1.136323i and cut3 cancellations. The U-star has alpha1.0219255, radial derivative1.80614, and actual WH radius7.29315→2.62232GeV. This is finite-gluon near-threshold counterterm enhancement, less sharp in H/Z than b Im. A cut3 right-star also has H=Z≈0 but cut1 defect20.71794GeV; it is not the same hosted intersection.

Across all reconstructed physical, single-star and merged points of these three cases, both gluon energies stay above 53.84GeV. LU scales lie 0.456--1.665; physical radial derivatives exceed 461.94 and threshold radial derivatives exceed 0.1748. These are observed values at three points, not uniform bounds. There is no sampled soft endpoint, extreme LU scale or unresolved root explaining these extrema.

== Verification and retained evidence
<verification-and-retained-evidence>
The original top2 replay had 60 bit-exact signed-extremum aliases, 48 unique Samples and 6 Arb controls, with 551 metric checks. The true c absolute-Im maximum ranks 4 by retained complex norm. A separately named existing replayer changed only take(2) to take(4) and descriptive scope, producing 12 real Arb controls and 605 passing checks. An additional 78 artifact checks confirm unchanged original Samples/reported values and exactly identical old six native controls. These replay calls do not add MC samples or alter any estimate.

Three new traces reproduce their retained Arb totals, factors and six events exactly. The unchanged numerical analysis reconstructs 56 single stars and 12 merged points from native occurrence/parent/center/root/alpha fields; maximum affine join discrepancy is below 4.87e-299GeV. Twelve actual native2l A/U multiplier checks agree within 4.53e-300. Native texts are retained. Only identity probes are reconstructed; selected threshold energy sums are checked numerically, and only the six physical-cut shifts additionally have serialized structural checks. No new root solve or directed channel-support certificate is claimed.

The artifact extension is #link("gl638_hosted_joint_gate/mc_confirmation/maximum_extension_build7/provenance.json")[maximum\_extension\_build7];, with exact source/hash index, top4 replays, immutable JSONL traces, native analyses and the same-point factor audit. The original #link("GL638_MC_CONFIRMATION.typ")[confirmation result] and #link("GL638_MAXIMUM_PHYSICS_ATTRIBUTION.typ")[screen maximum study] remain separate. Component norm ratios are not signed rate fractions. Finite maxima and two regular H/Z approaches cannot establish global boundedness, soft integrability or finite variance; those require the separately scoped soft/replacement study.
