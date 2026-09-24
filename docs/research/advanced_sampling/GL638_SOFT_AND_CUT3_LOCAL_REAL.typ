= GL638: real soft limits, Cut1 replacement, and a Cut3 channel
<gl638-real-soft-limits-cut1-replacement-and-a-cut3-channel>
The cut-only and joint catalogues do not keep real weights bounded on the tested single-soft approaches. The soft LMB channel strongly suppresses their growth. Replacing the ordinary Cut1 channel with its composed joint channel is worse here than keeping both. Separately, adding the supported Cut3 right `[5,10]` channel reduces the recorded real maximum locally from #strong[12.21898 to 5.78174 pb];. These are bounded local diagnostics, not new Monte Carlo selection or a global finite-variance proof.

== Scope and acceptance
<scope-and-acceptance>
The frozen build7 client evaluated the same full GL638 physical settings, six cuts, subtraction terms and multipliers. The soft scan contains two hard anchors, single soft edge13, single soft edge14 and simultaneous double-soft approaches, five nominal radii from `1` to `1e-8 GeV`, and four catalogues: six cuts; six cuts plus the composed Cut1 joint; five other cuts plus that joint replacing Cut1; and six cuts plus `lmb(6,12,13,14)`. All use the same generating channel `lu_cut_0`; the actual native forwarded points are exactly equal across catalogues.

All #strong[120 soft rows / 480 native evaluations] completed with no missing or failed branch. The final existing-owner audit passes #strong[15,050 checks];, including current complex total/six-cut precision criteria, event identities, signed CT decomposition, actual Grid probabilities and full factors. Ordinary calls selected Double 88, Quad 104 and Arb 48 times; each has a forced-Arb control. The maximum total complex-norm discrepancy is `4.21e-7`, within its existing precision-specific budget. No physical tolerance was changed.

Two inherited auxiliary roundoff assumptions were corrected with their original failures retained: Quad is a 106-bit DoubleFloat (`epsilon=2^-105`), and its default display writes only 31 significant decimal digits. The factor audit retains its 16-epsilon arithmetic allowance and separately propagates half-last-decimal-place uncertainty of the displayed operands/result. The one arithmetic-only miss remains visible in the final record. This is a serialization correction, not a looser physical comparison.

The raw signed routing is `q13=q3-q10+Qsp`, `q14=q7-q3-Qsp`, with exact CM `Qsp=0`. The largest actual-source/input coordinate change is `1.30e-13 GeV`; the smallest soft scales differ from their nominal values by at most `1.67e-5` relatively. Every slope uses the actual forwarded soft scale. Returned native cut energies confirm the soft momenta remain soft after physical rescaling: their energy/source-scale ratios lie between `.7621` and `1.0031`. No second LU solve is used.

== Single-soft real weights
<single-soft-real-weights>
Without the soft LMB, the measured density tends to a finite value and resolved real weights grow roughly as `1/lambda`. With the soft LMB, the density grows as `lambda^-2`, and real weights decrease overall. The table uses the smallest actual scales near `1e-8 GeV`; signs are retained.

#figure(
  align(center)[#table(
    columns: 6,
    align: (auto,right,right,right,right,right,),
    table.header([Hard anchor], [Soft edge], [Cuts Re Y (pb)], [Additive joint Re Y (pb)], [Replacement Re Y (pb)], [Cuts + soft LMB Re Y (pb)],),
    table.hline(),
    [real\_composed], [13], [5.406815e+13], [6.30795e+13], [6.463766e+13], [6.662206e-08],
    [real\_composed], [14], [-3.134226e+13], [-3.656597e+13], [-3.806645e+13], [-6.891519e-10],
    [real\_ordinary], [13], [-1.290969e+14], [-1.506131e+14], [-1.544718e+14], [-3.698944e-08],
    [real\_ordinary], [14], [3.789865e+13], [4.421509e+13], [4.517445e+13], [4.667849e-10],
  )]
  , kind: table
  )

The composed joint component has #strong[zero inverse density at all 30 tested soft points];. Adding it leaves the unweighted density sum exactly unchanged while the actual flat-grid channel count grows 6→7; consequently each complete real weight is multiplied by exactly 7/6. Replacement removes a positive Cut1 density contribution. On the resolved single-soft points it gives weights 1.192--1.215 times the cuts control, and 1.022--1.041 times the additive catalogue. Fewer channels therefore do not improve this sector.

Write `|Re F|~lambda^-p` and `q_eff~lambda^-s`. In a d-dimensional radial soft measure, the absolute integral has margin `d-p`, and the second moment has margin `d-2p+s`; positive margins are the radial convergence condition if the scaling holds uniformly over angles. For the d3 single-soft rays, descriptive last-three-point slopes give p=.893--1.177. Cut/joint q has s≈0, giving positive absolute margins 1.82--2.11 and second-moment margins .646--1.214. Soft-LMB q has s≈2, improving the latter margins to 2.646--3.214.

This evidence is consistent with finite absolute integral and finite variance along these single-soft directions, despite unbounded cut/joint weights. It is not a global proof: only two angular directions were examined, and the real coefficient is not uniformly asymptotic. One edge13 series changes sign and its final adjacent slopes vary from .0085 to 2.0677; the full signed values and adjacent slopes are retained. The soft LMB is the only tested coverage that suppresses the observed growth, making it useful for the bounded-maximum objective even where the directional variance margin is already positive.

== Double-soft real resolution
<double-soft-real-resolution>
The simultaneous double-soft real values cannot support a power fit. At forced Arb they are only about .0007--.008 epsilon times the sum of complex bare/CT term norms. All 40 double-soft methodpoints are below the existing 32\*n epsilon summation scale; the eight real curve fits are explicitly suppressed. Their changing tiny signs and apparent lambda^-5 slope are numerical remnants, not a measured real divergence. The soft-channel density itself scales as lambda^-4 on these directions, but this gives no d6 real integrability or variance conclusion without a resolved real signal.

The existing complex accuracy criterion can accept a tiny real leakage beside a large imaginary result. For example, one smallest double-soft source has ordinary Quad Re Y≈3.43e20 and Im Y≈8.32e47, while forced Arb retains the same imaginary result but Re Y≈1.47e-248. Thus “valid” under the current complex norm does not imply relative accuracy of an almost-zero real component. The present real-only analysis uses resolved forced-Arb values and leaves these double-soft cases unresolved; no runtime precision criterion was changed.

== Supported Cut3 right channel
<supported-cut3-right-channel>
The native expression is `then(block(lmb(4,12),phase_space(cut(2,4,12))),complement(7),block(lmb(5),at_cut(cut(2,4,12),surface(5,10))))`, parent `[4,5,7,12]`, on cut 3. Its LU-h profile acts on the phase-space prefix. The side surface uses explicitly changed global power 2. A separate global-power control was retained; it is exactly unchanged from power 1 in q and complete Y at all seven probed points.

The Cut3 scan completed #strong[21 rows / 84 native calls];, with #strong[2,681 checks passing];. At the recorded maximum, adding the channel increases the unweighted density sum 2.41528 times. Including the actual 7→8 channel probability gives a 2.11337-fold full proposal-density gain and reduces Re Y by 52.68%. All three catalogues again give exactly the same actual native source point.

#figure(
  align(center)[#table(
    columns: 4,
    align: (right,right,right,right,),
    table.header([Actual Cut3 sphere distance (GeV)], [Control Re Y (pb)], [Added-channel Re Y (pb)], [Full density gain],),
    table.hline(),
    [27.99381709], [12.218979], [5.7817416], [2.1133733],
    [-10], [9.1158933], [2.4357134], [3.7425969],
    [-1], [9.7881092], [0.89754552], [10.905418],
    [-0.1], [9.8681885], [0.29979082], [32.916914],
    [0.1], [9.8861684], [0.38552196], [25.643593],
    [1], [9.9678068], [1.2144391], [8.2077454],
    [10], [10.810282], [3.690723], [2.9290418],
  )]
  , kind: table
  )

The sphere distances are recovered algebraically from returned native Cut3 leg data and the actual common spatial rescaling. The original CM surface is `2sqrt(|q10|^2+173^2)-1000=0`. This does not identify a CT star. The pure right `[5,10]` local CT has a real contribution below `1.01e-336` in raw-density units; it is numerically negligible rather than an exact serialized zero. Its integrated real contribution remains broad, about 1.19--1.24e-36 over the six sphere probes versus total Re F about .47--.54e-36. These signed raw-density terms undergo cancellations and are not rate fractions. The new map improves local coverage without establishing a real divergence on the sphere or a global variance improvement.

== What remains established and open
<what-remains-established-and-open>
The user's leading criterion is to minimize known regions with unbounded real weights, ideally leaving none. Explicit soft LMB coverage is therefore retained in the next leading configuration; it is not removed merely because a finite Monte Carlo sample happens to have a smaller error without it. Keep the ordinary Cut1 coverage when assessing the joint map; the tested replacement is worse in the soft sector. Retaining explicit soft coverage addresses the observed single-soft weight growth. The new Cut3 channel helps the recorded real maximum locally. A future post-optimization comparison may test these coverages with learned channel probabilities and explicit channel summation, as already requested; this note selects no new Monte Carlo winner and changes no historical confirmation result. New replacement/Cut3 catalogue normalization and any longer sampling comparison remain pending. Neither the two-direction soft scan nor the finite Cut3 profile proves a global bound or finite variance.

The prepared leading nine-channel card contains all six standalone LU cut channels, the composed Cut1 joint, Cut3 right `[5,10]`, and `lmb(6,12,13,14)`. Its matched eight-channel control omits only the Cut3 channel. Both use global power 2, and their non-sampling settings are identical. Static card checks passed; native inventory, normalized catalogue acceptance and actual local HZ/soft/recorded-real-maximum gates remain required before any statistical comparison. The cards are under `leading-real-soft-cut3/` in the retained protocol. They are a coverage proposal, not a measured Monte Carlo winner.

== Retained evidence
<retained-evidence>
Lossless inputs, native outputs, exact driver/card/build hashes, prior diagnostic failures and final analyses are indexed in `gl638_hosted_joint_gate/local_soft_cut3_build7/artifact_hashes.json`. No saved-state payload, compiled executable or library is copied. Original runtime/pilot directories remain unchanged.

Final soft analysis SHA256: `33650951043376970164fdbda398a53669fa3af6d97cb303c7e4d9834e0e6083`. Soft analyzer SHA256: `4da70a97d7a91db0622f117f7b8b57199f8333d0d5599d8fad95edfa098644de`. Cut3 analysis SHA256: `8cb217f0848cad806a38e629d4ea9e7432e6fa3c64419313990363306fa60e36`.

An independent Decimal420 audit passes 113 additional checks and authenticates all 120 native point files. It reproduces the single-soft slopes, probability factors and radial margins, checks the explicit Quad serialization allowance, and confirms that all double-soft real fits and ratios are suppressed. Its script and results are retained under `local_soft_cut3_build7/independent_audit/`.
