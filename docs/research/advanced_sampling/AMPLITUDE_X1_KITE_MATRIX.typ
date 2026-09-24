= X1 massive-kite sampling pilot
<x1-massive-kite-sampling-pilot>
The mixed C-surface plus LMB proposal improves the measured real-component error and maximum weight in this small pilot. The C-only proposal performs worse. This is amplitude evidence for the generic channel machinery, not a GL638 result or a bounded-weight claim.

Exact source: `3f9b1a28473f313d6242a556ada2e8dbcd6f6477`; optimized binary SHA256 `f700ee575153067ccc77894ea776a96f139105812a15014d094dd67ba7db6d61`. All 18 physical orientations are summed. Three seeds (101, 202, 303) each use 10,000 draws in five 2,000-draw batches, one core, power 2, real-component training, and exact map-density partition weights. Threshold subtraction stays enabled. All twelve runs have zero reported NaN/unstable samples.

#figure(
  align(center)[#table(
    columns: 5,
    align: (auto,right,right,right,right,),
    table.header([Channels], [Mean Re / Im], [RMS reported error Re / Im], [Max absolute weight Re / Im], [Mean evaluation µs],),
    table.hline(),
    [optimized\_lmb], [-1.3773e-05 / 6.6078e-06], [3.276e-07 / 1.015e-06], [2.297e-03 / 7.694e-04], [76.53],
    [all\_lmb], [-1.3773e-05 / 6.6078e-06], [3.276e-07 / 1.015e-06], [2.297e-03 / 7.694e-04], [76.26],
    [surface], [-1.4235e-05 / 7.6714e-06], [8.685e-07 / 1.752e-06], [6.352e-03 / 6.705e-03], [69.61],
    [mixed], [-1.3794e-05 / 5.9924e-06], [2.794e-07 / 9.728e-07], [1.205e-03 / 8.324e-04], [79.16],
  )]
  , kind: table
  )

Here C is `surface(2,4,6)` in parent/subspace `[4,6]`. Its partner D coincides geometrically at rest, so the performance channel list does not duplicate it. Both LMB selectors produce identical values, errors and maxima for each common seed; the reference catalogue has eight LMB channels. Mixed adds C as a ninth channel.

Relative to LMB, mixed has reported-variance ratios 0.727 (real) and 0.918 (imaginary), with an approximately 3.4% larger evaluation cost. The real maximum falls by 48%; the imaginary maximum rises by 8%. C alone has variance ratios 7.03 and 2.98 and larger maxima. The profile means agree within the uncertainties reported by these short runs; the seeds are paired, so a formal significance claim would require their covariance. Evaluation cost is measured inside the integrand; end-to-end times and every maximum point are retained in the JSON.

#figure(
  align(center)[#table(
    columns: 3,
    align: (auto,right,right,),
    table.header([Reference channels], [Normalization (target 1)], [Raw second moment (target 9.1025)],),
    table.hline(),
    [optimized\_lmb], [0.98341131], [8.97610299],
    [all\_lmb], [0.98341131], [8.97610299],
    [surface], [0.98705210], [8.98678016],
    [mixed], [0.98740466], [9.01105246],
  )]
  , kind: table
  )

Reference checks use 8,192 shared Halton points per channel and a normalized Gaussian of width 1.2 centered at `[0.4,-0.3,0.2,-0.2,0.1,0.35]`. All pass the existing absolute tolerances 0.02 on normalization and 0.2 on the second moment. Per-channel normalization dispersions are retained in the JSON; they are not combined as independent errors because the draws are shared. This original driver did not collect moment dispersion. It links the exact frozen X1 libraries and uses that checkpoint's f64 reference path; the newer native reference-rescue changes are outside this evidence.

The portable #link("AMPLITUDE_X1_KITE_MATRIX.json")[JSON artifact] contains binary/library/source hashes, graph and generation inputs, cards, resolved settings, actual orientation keys, original CLI argv, full sample counts, timing and stability statistics, reference contributions, and retained maximum points. Scratch paths are provenance: reproduce by extracting the embedded inputs into a new directory and substituting those paths. `tests/resources/graphs/massive_kite.dot` is the same graph up to whitespace.

The same-binary #link("AMPLITUDE_X1_DOUBLE_BOX_MATRIX.typ")[boosted-kite and rest-frame double-box matrices] are now complete with their 18/98 orientations and reference gates. Next use paired wall-time budgets to test whether the mixed gain survives equal computational cost. Proper-subspace channels and the pending conditional maps need their own gated binary before entering this comparison; keep that comparison bounded to the same existing graphs.
