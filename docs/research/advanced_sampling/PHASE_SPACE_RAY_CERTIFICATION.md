# Phase-space ray: retained foreign-inverse failure

The seven-candidate GL638 preflight found a finite-precision routing failure in
the existing full-frame Cut1 inverse. The exact retained Arb1000 point reproduces
it in a small graph-based test. Routing each fixed velocity before radial
scaling resolves that point with the **same solver and residual budget** and
satisfies a directed check of the original graph equation.

The production change is confined to `phase_space(cut(...))` callbacks. The
rebuilt GL638 point replay now passes; complete seven-candidate acceptance is
still pending. The original failed preflight remains preserved. No physical
tolerance, root policy, LU-h profile, channel definition or proposal precision
changed.

## Failed preflight and exact isolation

The source is the same build4 used for the
[local H/Z study](LOCAL_HZ_BOUNDING.md): base `6e9bf401db75cb13ec20e18377a444bc2fcd4f33`
plus patch `1d78e98b73e6a8dafaf6e446353c2a1607dc0724f8ea7d1c2d88c4962b81bba7`,
equivalent to committed core source `8006662289c17f177469c67a10a274b190aa2ccf`.
The physical state retains all 936 orientations, six cuts and the same CT/UV
metadata. Both completed diagnostic runs preserve all 35 saved-file hashes.

| Catalogue | Last reference result | Outcome |
| --- | --- | --- |
| Optimized LMB | 32768 points; normalization 1.006498; raw moment 1053863.74 | Preflight passed |
| Six LU-h cuts | 32768 points; normalization 1.031806; raw moment 1126940.37 | Preflight passed |
| Six cuts + original H/Z joint | Halton index 3363 in the 8192-point batch | Numerical failure |
| Remaining four candidates | Not reached | No result |

The reference targets are normalization 1 and raw second moment 1085600 GeV².
The unchanged acceptance limits are 0.06 absolute and 0.08 relative, respectively.
The deterministic-Halton dispersion is diagnostic, not a confidence interval.

The one-point replay retains exact binary64 cube bits and all native forward
points. All six ordinary generated rows and their Cut1 inverses are identical
between `cuts` and `cuts_joint`. The additional joint row forwards successfully;
its original Cut1 inverse alone fails. The same error is retained by both the
selected-joint and summed reference evaluations. Thus the matrix contains 13
successful forwards and 12 successful foreign inverses, with one failed inverse.

This point uses the joint component's **normalized ordinary fallback**, because
the compact disk was declined. It is not a compact H/Z corner or an R=0 seam.
Its master momenta include two vectors of order 10⁶ GeV whose difference is of
order 10³ GeV. The failure precedes native materialization, the Gaussian or
physical body, rotations, CT evaluation and grid accumulation.

## First differing arithmetic boundary

The Cut1 equation is the original `surface(2,6,10)`, with masses 125, 173 and
173 GeV and temporal shift −1000 GeV. In generation parent `[3,4,7,10]`, its
spatial routes are `q2=q4−q3`, `q6=q4−q3+q10`, `q10=q10`. Its inverse normalizes
all twelve master components before solving a directional radial root; dropping
host-null components would define a different radial coordinate.

The old callback formed each scaled master momentum before subtracting the
large common components. The existing `Esurface::routed_ray` instead forms the
fixed signed velocities first, then evaluates their affine energies at the
radius. Both retain every original energy, mass and external shift. The latter
is already the explicit LU preparation convention. Unqualified amplitude,
static and fiber callbacks retain their prior arithmetic convention.

The fixture parses all twelve 303-significant-digit native point strings
directly into Arb1000 and requires exact format roundtrip. It uses the actual
GL638 graph parser, verifies the parent and signed routes, and supplies the two
recorded external momenta `(500,0,0,±500)`. No generated state is loaded and no
joint map is regenerated for this comparison.

| Callback / root owner | Iterations | Residual | Result |
| --- | ---: | ---: | --- |
| Scale master point, then route / strict | 17 | −1.1467943344×10⁻²⁹⁶ | Failed |
| Same callback / existing bracket-resolution wrapper | 17 | Same | Failed |
| Route fixed velocity, then scale / strict | 8 | +9.5566194535×10⁻²⁹⁹ | Passed |
| Same prepared ray / existing wrapper | 8 | Same | Passed |

Every call uses guess 300, origin scale 529, tolerance multiplier 64, limits
2048/96 and epsilon 2⁻⁹⁹⁹. The old printed bracket width is approximately
4.9×10⁻²⁹⁶; adding derivative times that width does not rescue its residual.
No wrapper or tolerance adjustment is part of the fix.

The independent directed original-routing evaluation at the new root gives
an interval near **+3.0645044730×10⁻²⁹⁹**, of width 1.26744×10⁻⁶¹³. Both
endpoints lie within the unchanged `64·2⁻⁹⁹⁹·529 = 6.3193146136×10⁻²⁹⁷`
budget, with about 206-fold margin. The interval is positive and does not contain
zero: this certifies the requested approximate-root accuracy, not an exact root.
The fixture now asserts this original-equation bound in addition to the native
callback checks.

The changed production path and its existing controls pass all 15 selected
core tests (451.960 s), including implicit roots, LU-h profiles, standalone
cut sampling, conditional raised cuts and hosted-joint normalization. Core
checking, formatting and clippy also pass; none of the 52 clippy diagnostics
touches a changed line. The optimized core/API/CLI build also passes in
426.392 s with the frozen patch unchanged.

## Rebuilt GL point replay

Build5 reuses the same diagnostic driver source `254a7465…`, linked against
the new core/API artifacts. At Halton3363, all **13 native forward rows remain
exactly identical** to the old build. All 13 foreign Cut1 inverses now succeed,
as do both summed references and the selected joint reference. The original
cube bits, selected Sample, physics settings, six-cut/936-orientation/CT
inventory and all 35 saved-state hashes remain unchanged.

The [independent matrix audit](gl638_hosted_joint_gate/phase_space_ray/replay_build5/independent_matrix_analysis.json)
confirms these comparisons. The run completed in 81.800 s with one worker. The selected joint report is
finite and returns zero Gaussian weight at this far-tail ordinary-fallback
point. This establishes the repaired pointwise pipeline; it supplies neither
a one-point normalization check nor a physical-integrand acceptance result.
The full seven-candidate preflight is a separate, still-pending gate.

## Evidence and reproduction

The [archive](gl638_hosted_joint_gate/phase_space_ray/artifact_hashes.json)
records uncompressed hashes and storage paths. It retains the failed preflight,
exact minimal matrix, cards, manifests, run records, native fixture log and
compiler/link provenance. Repeated giant summary payloads are replaced by
compact summaries retaining the original summary hash and size; the individual
mode reports remain byte-preserving gzip copies. No binary or state payload is
duplicated. The initial fixture port-count failure is retained and explicitly
identified as a correction made before the callback comparison.

Run the read-only numerical audit from any checkout:

```bash
python docs/research/advanced_sampling/gl638_hosted_joint_gate/phase_space_ray/analyze.py
```

Exact client commands and environment are in
[provenance.json](gl638_hosted_joint_gate/phase_space_ray/provenance.json).
The minimal old-binary replay exits successfully because it records all
independent outcomes; its selected and summed joint rows still report failure.
The core reproduction is
`cff::esurface::tests::phase_space_cut_root_replays_canonical_foreign_inverse_cancellation`.
The rebuilt binder replay passes. The full seven-candidate preflight must also
pass before the equal-N physics comparison proceeds. This correctness repair supplies no
new cost, variance, integral or global boundedness result.
