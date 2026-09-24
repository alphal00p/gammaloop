= GL638 local H/Z scaling with cut channels
<gl638-local-hz-scaling-with-cut-channels>
Joint H/Z focusing keeps the complete sampled weight nearly constant on both tested hard directions when combined with the six cut channels. Applying Cut1 LU-h focusing before the joint block lowers these local plateaus by about 7.9 times. Cut-only mixtures retain approximately `1/R` growth. These are measured finite-range directional results, not global bounds or Monte Carlo rankings.

== Matched calculation and channel laws
<matched-calculation-and-channel-laws>
All seven catalogues use optimized build5, whose Rust source is exactly commit `4d192b5be`: the #link("gl638_hosted_joint_gate/local_cut_matched_hz/source_commit_proof.json")[source proof] matches its committed diff to the frozen build patch. The calculation retains all 936 orientations, all six cuts, the physical numerator, local and integrated UV, WH/WF and canonical physical centers. Both runs preserve all 35 saved-file hashes. Physics settings and stability tolerances are unchanged across modes.

The six cut channels target `(2,6,12,13)`, `(2,6,10)`, `(2,6,7,13,14)`, `(2,4,12)`, `(2,4,10,13)` and `(2,4,7,14)`. Each uses the same LU-h log-logistic radial profile with broad fraction 0.02. Soft coverage adds the LMB `[6,12,13,14]`. The original joint uses ordinary prior sampling, followed by the original H=`surface(2,4,12)` and Z=`surface(3,10,13)` under cut1. The genuinely composed channel instead performs:

```text
phase_space(cut(2,6,10)) on (6,10)
  -> complement(7)
  -> joint H/Z on (3), hosted by cut1
```

Its ordered parent is `[3,6,7,10]`; the first block is six-dimensional. It retains the auxiliary raw LU radial coordinate. The actual preceding outputs prepare the same host and joint geometry; the determinant includes every block and the physical-to-raw host transformation. The six full-volume cut siblings supply support outside a compact joint chart.

== Measured slopes and endpoint weights
<measured-slopes-and-endpoint-weights>
There are 84 map/point rows: seven methods on twelve anchors, comprising six radii on each of two hard directions. Every row retains ordinary and forced-Arb physics at both the supplied binary64 raw point and its actual sampled source, for 336 native evaluations. All 4,952 independent checks pass.

The radius is the actual native joint diagnostic `R=sqrt(H_global^2+Z_global^2)`, approximately `0.251624` to `2.51624e-6 GeV`. WH\'s finite-residual H-cut differs from H-global by the host residual; the original target equation is not rewritten. The fits use each method\'s last three actual forwarded radii, with `beta=-d log|Y|/d log R`. Y is the complete sampled weight in pb, including the actual grid and channel probabilities once. The full raw physical norm grows approximately as `1/R` on both directions.

#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,right,right,right,),
    table.header([Catalogue], [Channels], [beta Y, −Z / +Z], [Smallest-R |Y| \[pb\], −Z / +Z],),
    table.hline(),
    [Optimized LMB], [6], [1.000003 / 1.000002], [1,455,081.9 / 1,677,810.1],
    [Six cuts], [6], [1.000002 / 1.000002], [10,353,404 / 11,938,192],
    [Six cuts + soft LMB], [7], [1.000003 / 1.000002], [6,370,791.9 / 7,345,963.7],
    [Six cuts + direct joint], [7], [2.561e−6 / 2.192e−6], [0.0072421491 / 0.0083506996],
    [Six cuts + direct joint + soft LMB], [8], [2.573e−6 / 2.203e−6], [0.0082767419 / 0.0095436566],
    [Six cuts + composed Cut1 LU-h→joint], [7], [2.550e−6 / 2.181e−6], [0.00091684014 / 0.0010571802],
    [Six cuts + composed Cut1 LU-h→joint + soft LMB], [8], [2.551e−6 / 2.182e−6], [0.0010478173 / 0.0012082059],
  )]
  , kind: table
  )

The direct-joint component is the same law as in the #link("LOCAL_HZ_BOUNDING.typ")[earlier local study];; only its full-volume siblings change. Its positive density contributes `C/R` in this regular patch, so the other positive partition scores cannot remove that damping. The composed channel changes the prior density coefficient and improves these two local plateaus. It does not establish that this channel is best over the whole domain.

Adding the soft channel leaves the selected joint point, cube and full Jacobian exactly unchanged in all 24 paired rows. Here its density contribution is small, while the uniform channel-selection factor changes from 7 to 8, increasing the local joint weights by approximately `8/7`. The soft LMB may still be useful in soft regions or for global variance. Conversely, the poor cut-only coefficients on these hard rays do not rank their performance elsewhere.

Different maps recover slightly different binary64 cubes. Their actual source coordinates differ from the intended raw point by at most `1.481737e-13 GeV`; the radius differs from the common native anchor by up to `4.415459e-8` relatively. Each physical value is combined with the density and factors at its own source. The diagnostic joint inverse is not a separate original-energy equation oracle or a certificate for arbitrary CT-star normals.

The largest ordinary/Arb norm discrepancy is `2.034931e-7`, within the unchanged Double budget `1e-6`. Event-sum and CT-decomposition checks pass, with maxima `3.385079e-15` relative to the complete norm and `4.755161e-16` relative to the sum of term norms. The largest `|J*w*q_sum−1|` is `1.567125e-291`. The postprocessor uses the actual 106-bit Quad rounding model; physical Double, Quad and Arb tolerances remain `1e-6`, `1e-10` and `1e-12`.

Moved CT-star H/Z loci, soft13/14, tangent or small-Gram geometry, shrinking support, radial-center limits, tails and coincident roots remain separate. Neither this finite scan nor the compact kernel proves a global maximum or finite variance. Full-state normalization, matched-count MC errors, retained maxima attribution and a supported central estimate remain separate gates.

== Evidence and reproduction
<evidence-and-reproduction>
The #link("gl638_hosted_joint_gate/local_cut_matched_hz/artifact_hashes.json")[archive index] records uncompressed hashes. All 84 point files are unchanged gzip copies; inventories and accepted settings are shared once, while compact mode/summary records preserve the other fields and original raw hashes. The prior local archive is untouched. The two old extension seed files are linked explicitly. The retained hash manifest matches all 56 base point files; its creation after the extension is not presented as an earlier snapshot.

The driver source is the same `b53df2b4…`, linked against build5. Both processes exit zero: 422.513 s for eight anchors and 259.518 s for four extension anchors. They use one worker and `RAYON_NUM_THREADS=1`, while a separate 20-worker normalization preflight runs concurrently. These elapsed times establish no sampling-cost or parallel-scaling result. Exact commands, build/source links, state hashes and the cumulative-extension-RSS caveat are in #link("gl638_hosted_joint_gate/local_cut_matched_hz/provenance.json")[provenance.json];.

For postprocessing only, decompress the archived script and point files into scratch directories, including the two linked old seed points, then run:

```sh
python analyze_local_approach_extension.py BASE8_DIRECTORY EXTENSION4_DIRECTORY \
  OUTPUT_DIRECTORY base8-frozen.json OLD_SEED_DIRECTORY
```

The #link("gl638_hosted_joint_gate/local_cut_matched_hz/analysis.typ")[independent report] contains signed endpoint weights and the full slope table. Replaying from the archive reproduces that Markdown byte-for-byte and the full numerical JSON except for relocated input-path keys; all 4,952 checks pass again. No saved state is loaded by this analysis; complete event weights are not multiplied again.
