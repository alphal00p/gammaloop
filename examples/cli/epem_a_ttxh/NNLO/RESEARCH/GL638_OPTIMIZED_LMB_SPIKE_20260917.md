# GL638 optimized-LMB error spike, 17 September 2026

The new maximum is a genuine, reproducible large threshold-counterterm weight.
It is associated with Cut 1 and approaches the H/Z region more closely after
threshold projection than at the sampled physical point. It is not a soft-gluon
endpoint or a loss of numerical precision. One point does not establish an
asymptotic scaling law or prove that a particular advanced channel bounds it.

The point first enters `integration_result_iter_0231.json` in
`GL638_optimized_lmbs_workspace`. Its real weight is **−310.7628500264454 pb**;
the same point has imaginary weight **−3397.274568766552 pb**. The sampled channel
is #3, LMB **[4,7,13,14]**, at **sqrt(s) = 600 GeV**. The retained coordinates are:

```text
0.65326967347543075 0.99071624997667884 0.34598647260057075
0.73874307614468926 0.90145414024659398 0.087144592922320935
0.58803498486189254 0.099566804759442573 0.80480600515392220
0.67392249894159018 0.66470428853177543 0.060665899033751615
```

## Numerical control and normalization

Ordinary inspection gives `−1.949730998643977 − 21.31455396636133 i` with unit
outer integrator weight; forced 1000-bit Arb evaluation gives
`−1.9497309986438203 − 21.314553966359853 i`. The relative discrepancies are
8.04e−14 and 6.93e−14, respectively. The Arb rotation comparison reports
1.92e−299. Thus the outlier remains after high-precision rescue.

Inspection includes the sampling-map weight but sets the outer adaptive-grid
and discrete-selection weight to one. Dividing the stored maximum by this
inspection independently in Re and Im reconstructs that outer factor as
159.38755153535465 and 159.3875515353564. This factor is inferred, not read from
a retained native Sample. It is applied exactly once in the tables below.

## Physical cuts and counterterms

| Physical cut | Edges | Real contribution to stored sample [pb] |
|---|---|---:|
| 0 | [2,6,12,13] | +0.000018071 |
| 1 | [2,6,10] | **−329.679856** |
| 2 | [2,6,7,13,14] | negligible |
| 3 | [2,4,12] | +18.916967 |
| 4 | [2,4,10,13] | +0.000020908 |
| 5 | [2,4,7,14] | +0.000000000040 |
| **Sum** | | **−310.762850** |

The two main real Cut-1 terms are:

| Threshold / component | Variant | Multiplier | Real contribution [pb] |
|---|---|---:|---:|
| A = [7,8], integrated CT, component 35 | native_2l, subspace [3,7] | 0.755138093 | **−239.631944** |
| Z = [3,10,13], integrated CT, component 29 | native_2l, subspace [3,7] | 1 | **−88.600466** |

Their native parent is [3,6,7,10]. Cut 3 also has large individual terms, for
example components 7 and 13 contribute −221.102904 and +220.093369 pb, but they
largely cancel. This is why attribution to the largest individual component
alone would be misleading; the Cut-1 total controls the outlier.

## H/Z geometry

Define `H = E2 + E4 + E12 − Q`, `Z = E3 + E10 + E13 − Q`,
`P = E3 + E12 − Q`, and `Q = 600 GeV`. The multiplier uses
`W_H = H² / (H² + (P Z/Q)²)` on the Cut-1 host. All rows below satisfy the
Cut-1 host equation `E2 + E6 + E10 − Q = 0` to below 1e−290 GeV, so they
cannot be confused with an H=0 point that belongs only to Cut 3.

| Configuration | H [GeV] | Z [GeV] | sqrt(H²+Z²) [GeV] | sqrt(H²+(PZ/Q)²) [GeV] | E13 / E14 [GeV] |
|---|---:|---:|---:|---:|---:|
| Physical Cut-1 point | +5.959606 | +2.724014 | 6.552645 | 5.984287 | 126.686 / 183.432 |
| A counterterm star | **+0.580884** | **−1.588175** | 1.691072 | **0.668461** | 125.062 / 179.141 |
| Z counterterm star | +2.573497 | ~0 | 2.573497 | 2.573497 | 125.654 / 180.737 |

The physical LU root is `t* = 0.4930766295125468`. For the dominant A star,
the radius changes from 301.251976656 to 294.204829643 GeV, with
`alpha = 0.9766071343634826` and `d eta / d r* = 1.3612058436253613`.
Neither the root nor its radial derivative is ill-conditioned. The actual
multiplier denominator radius is 0.0011141 Q: much smaller than the sampled
6.55 GeV H/Z radius, but still a finite displacement, not an asymptotic scan.

The A integrated residue is already −1.990966409 at unit outer weight before
its multiplier; `W_H` reduces it to −1.503454577. Its projected point is also
near the other threshold Z, whose energy defect is only −1.588 GeV. The trace
has first-order LU and threshold residues and `t_variable=None`; no
higher-power LU derivative of the multiplier is involved. The evidence supports
enhancement of a projected threshold residue near intersecting surfaces,
amplified by the adaptive proposal weight. It does not isolate a universal
power law from this one event.

The pointwise support and density of the advanced H/Z channel at this exact
canonical point have not been replayed here. Physical proximity alone therefore
does not establish the amount by which that channel would reduce this weight.

## Retained verification

The [artifact directory](gl638_spike_20260917/) contains the frozen iteration
snapshot, ordinary/Arb inspection JSON, runtime settings, exact command arrays,
compressed native trace, detailed component text and registry display.
[analysis.json](gl638_spike_20260917/analysis.json) retains all cut totals,
component weights and eight reconstructed Cut-1 stars.

The native-star reconstruction uses the measured identity-frame center,
complement and alpha without solving another threshold. It checks native-frame
joins, both radii and every selected threshold residual below 1e−290 GeV. It
also checks normal/Arb agreement and independently reconstructs the outer
factor from both complex components. Reproduce those offline checks with:

```bash
cd /common/dev/gammaloop_ir_safe_thresholds
examples/cli/epem_a_ttxh/NNLO/.venv/bin/python \
  examples/cli/epem_a_ttxh/NNLO/RESEARCH/gl638_spike_20260917/analyze.py
```

Numerical calls used an isolated copy of the saved state under the ignored
`workspaces/GL638_max_spike_diagnostic_20260917`; the production state and live
integration workspace were not modified. The initial detailed JSON route exposed
a separate enum-map-key serialization limitation. The forced-Arb detailed
report then rejected a separately stored, exponentially suppressed Original
real component of order 1e−341 at the f64 reporting boundary. Disabling only
additional-weight retention permits the ordinary and forced-Arb totals and
event results to serialize; the detailed ordinary CT output and already
completed native Arb geometry remain retained. Neither reporting limitation
causes the measured large weight, and no production source was patched here.
