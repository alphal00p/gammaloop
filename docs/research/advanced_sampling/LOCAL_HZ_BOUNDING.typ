= GL638 local H/Z weight scaling
<gl638-local-hz-weight-scaling>
The direct H/Z joint channel keeps the complete sampled weight nearly constant on the two tested hard approach directions, through five decades in the actual normal radius. Optimized LMB grows approximately as `1/R`; direct-H power2 grows as `1/sqrt(R)`. This is a finite-range directional result, not a global bound, a Monte Carlo variance measurement or an integral estimate.

== Calculation and exact comparison boundary
<calculation-and-exact-comparison-boundary>
Both completed runs use optimized build4, whose Rust source is exactly commit `8006662289c17f177469c67a10a274b190aa2ccf`: base `6e9bf401d` plus patch `1d78e98b73e6a8dafaf6e446353c2a1607dc0724f8ea7d1c2d88c4962b81bba7`. The #link("gl638_hosted_joint_gate/local_hz/source_commit_proof.json")[source proof] compares the committed Rust diff with the original build patch. The physical calculation retains all 936 orientations, all six cuts, the physical numerator, local and integrated UV, current WH/WF metadata and the canonical center repair. All 35 saved-state hashes remain unchanged in both runs.

The first run uses eight saved points: four radii on each of two mirrored hard H-positive/Z-positive-or-negative directions. The extension adds two decades per direction. There are 36 completed map/point rows across three methods, with 144 native physical results: ordinary and forced-Arb evaluations at the binary64 raw point and at each actual sampled source. The original 24 point files remain hash-identical.

The fitted radius is the actual joint diagnostic `R=sqrt(H_global^2+Z_global^2)`, approximately `0.251624` to `2.51624e-6 GeV`. It is not the historical WH-scaled radius. The map retains the original H=`surface(2,4,12)` and Z=`surface(3,10,13)` under host cut1 `(2,6,10)`; WH\'s finite-residual H-cut differs from H-global by the host residual.

Each method reconstructs its own binary64 cube, then forwards that cube through the canonical Arb owner. Those actual Cartesian points are close but not identical. The largest coordinate displacement from the intended binary64 point is `1.481737e-13 GeV`; in the extension the actual source radius differs from the common native comparison anchor by up to `4.415459e-8` relatively. Slopes and weights use each method\'s actual forwarded point, its own inverse density and physical result. No quotient combines raw binary64 physics with a density at a different native anchor.

== Measured scaling
<measured-scaling>
Let F be the full raw physical result at the actual source, `q_sum=sum_i q_i` the map-density sum, and Y the complete sampled weight in pb, including the actual discrete/continuous grid probabilities exactly once. The scan cards have six optimized LMBs, one direct-H power2 channel, or joint plus six LMBs. Their different channel probabilities are included in Y. These are the original local-scan cards, not the later matched six-cut candidate matrix. In particular, the tested hosted joint uses ordinary prior sampling; it does not apply LU-h radial focusing inside that same channel.

The table uses `beta=-d log|value|/d log R`, fitted to the last three actual radii. F has beta `1.000003` / `1.000002` on the negative/positive Z directions.

#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,right,right,right,),
    table.header([Method], [beta q\_sum, −Z / +Z], [beta Y, −Z / +Z], [Smallest-R |Y| \[pb\], −Z / +Z],),
    table.hline(),
    [Optimized LMB], [approximately 0 / 0], [1.000003 / 1.000002], [1,455,081.9 / 1,677,810.1],
    [Direct-H power2], [0.500247 / 0.500247], [0.499755 / 0.499755], [2,295.0718 / 2,646.3768],
    [Direct joint H/Z + LMB], [1.000000 / 1.000000], [0.000002640 / 0.000002271], [0.0072421491 / 0.0083506995],
  )]
  , kind: table
  )

The joint endpoint weights are approximately `(-0.0029120368,+0.0066308947)` pb for −Z and `(-0.00017059052,+0.0083489569)` pb for +Z. Their norms change by only `1.10520e-6` / `9.50639e-7` over the final decade. These are individual local sample weights, not cross-section estimates or measured global efficiency gains.

Independent Decimal380 postprocessing passes 2040/2040 checks. These include native validity, raw-point physics/cut invariance across methods, actual-source grid factors, selected `J*w*q_sum=1`, six-event sums and CT decompositions. The largest reciprocity residual is `5.435883e-293`. The maximum ordinary/Arb complex-norm discrepancy is `1.591925e-7`, within the ordinary Double budget `1e-6`; this does not promise relative accuracy for every negligible component. The largest event-sum discrepancy is `2.831740e-15` on the complete-result norm, and decomposition error is at most `3.229585e-16` relative to the sum of term norms. No physical tolerance changed.

== Scope of the result
<scope-of-the-result>
All actual forwarded points have regular compact joint diagnostics. In this patch the complete density supplies the expected `1/R` factor and damps the observed completed physical remainder. Shared s-only A/U projections preserve the original p-normal geometry, but native two-loop A/U, right-side and iterated CT stars can move it. The scan does not supply a proposal for those star loci.

Soft13/soft14, shrinking-F or radial-center limits, small Gram/circle margins, shrinking support, tails and angular hierarchies remain separate. The center repair restores a consistent physical prescription; it does not add CT-star sampling or define coincident A/P residues. A finite-variance conclusion would also require control of a neighborhood and all remaining coefficients. The next physics comparisons use actual normalized candidate laws and equal sample counts, including the requested cut/joint/soft combinations.

== Archive and reproduction
<archive-and-reproduction>
The #link("gl638_hosted_joint_gate/local_hz/")[archive] contains exact compressed point files, manifests, cards, logs, unchanged state hashes, frozen driver sources and build records. Large duplicate mode/summary reports are represented by hashes and compact execution records. No binary or physics state is copied. The archive\'s `artifact_hashes.json` uses uncompressed hashes. Both analysis scripts and their original results are retained unchanged; the extension JSON records absolute input-path keys, which naturally differ when replayed elsewhere.

The original driver source is `ab24de37…`; the extended source is `b53df2b4…`. They link the same build4 libraries. Both processes exit zero. The extension completes in 132.897 seconds including one state load; this is not a production timing, parallel-scaling or sampling-cost claim. The exact root invocations are:

```sh
RAYON_NUM_THREADS=20 GL_DISPLAY_FILTER=off \
  /tmp/gl638-hosted-joint-gate/drivers/gate-physics \
  /tmp/gl638-hosted-joint-gate/manifest-local-approach.json \
  /tmp/gl638-hosted-joint-gate/results-local-approach1 \
  optimized_lmb,direct_h_p2,joint_hz_plus_lmb local-approach 1 1 1 1337

RAYON_NUM_THREADS=20 GL_DISPLAY_FILTER=off \
  /tmp/gl638-hosted-joint-gate/drivers/gate-physics-extended \
  /tmp/gl638-hosted-joint-gate/manifest-local-approach-extension.json \
  /tmp/gl638-hosted-joint-gate/results-local-approach-extension1 \
  optimized_lmb,direct_h_p2,joint_hz_plus_lmb local-approach 1 1 1 1337
```

For postprocessing only, decompress each directory\'s `*.approach.json.gz` into separate scratch directories, then run:

```sh
python analyze_local_approach.py SCAN8_DIRECTORY OUTPUT_DIRECTORY
python analyze_local_approach_extension.py SCAN8_DIRECTORY EXTENSION_DIRECTORY \
  OUTPUT_DIRECTORY local-eight-frozen.json
```

The original analysis reproduces byte-for-byte. The extended Markdown and all numerical JSON content reproduce exactly; only relocated input-hash path keys change. Neither script loads the saved state or evaluates GammaLoop.
