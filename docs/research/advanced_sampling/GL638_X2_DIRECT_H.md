# GL638 direct-H sampling gate

The conditional H chart passes an actual saved-GL638 sampling acceptance gate at
commit `f2f64fb17`. This is a map/measure result, not a physical variance result.
Both the selected-orientation state and the full state with 936 production
orientation keys pass. The latter exposes one summed orientation selector; its
936 actual orientation vectors exactly equal the archived generation set.
Reference substitution bypasses physical amplitude evaluation, so these results
do not establish the required all-orientation physical cut/counterterm
cancellation gate.

The exact numerical records, eight original raw points, source/library/binary
hashes, state hashes and orientation key are in
[GL638_X2_DIRECT_H.json](GL638_X2_DIRECT_H.json). The original threshold metadata
and the full local-3D/integrated-UV generated payload are retained.

## Channel and independent checks

The [runtime card](gl638_x2_direct_h.toml) selects one canonical channel:

```toml
[sampling.channel_definitions.GL638.direct_H]
around = "then(complement(6,7,10),block(lmb(3),at_cut(cut(2,6,10),surface(2,4,12))))"
parent_lmb = [3,6,7,10]
subspace_lmb = [3]
on_cut = [1]
```

It uses `map_density`, surface power 2, and ordinary spherical/linear complement
coordinates with scale `b * e_cm = 300 GeV`. The graph's raw generation frame is
`[3,4,7,10]`. The host is cut `(2,6,10)` and H is the direct global surface
`(2,4,12)`, independently of whether H has a host threshold-CT association.
The stored cut-surface registry has six entries; this is distinct from the
selected orientation's active-cut mask.

Writing the prepared momenta as `p`, `a`, `t`, the independent check uses

\[
H=E_t(a+p)+E_t(p)-B,\qquad B=E_t(a+t)+E_t(t).
\]

Its center is `-a/2`; along unit direction `n`,

\[
R_H^2=\frac{B^2-|a|^2-4m_t^2}{4[1-(a\cdot n/B)^2]}.
\]

The host rescaling is solved independently from
`E_h(tau*a)+E_t(tau*(a+t))+E_t(tau*t)=1000 GeV`, with `mt=173 GeV`,
`mh=125 GeV`. At the chosen cube point, the production map and independent
ellipsoid predict the same prepared radius, `162.39857492083064 GeV`.
This includes the native-to-generation affine map and the conditional `tau^-3`
volume factor.

A full 12-dimensional central finite-difference determinant gives:

| Cube step | Relative error against complete map Jacobian |
|---:|---:|
| `1e-5` | `1.46e-9` |
| `2e-6` | `1.46e-10` |

The complete Jacobian is `3.7545074795e31 GeV^12`. Cut-to-active off-diagonal
entries are nonzero (largest magnitude about `476.28`); this check includes
those conditional derivatives, rather than multiplying independent chart
Jacobians by assumption.

## Normalized-reference evidence

The target is a normalized 12-dimensional Gaussian of width `300 GeV`, centered
at `[30,-20,10]` repeated four times. Its raw second moment is
`1,085,600 GeV^2`. One graph and one channel are explicitly asserted; orientations
are summed, so the normalized target remains 1 without an orientation-count
factor.

| Production orientations | Halton points | Normalization | Raw second moment | Relative moment deviation |
|---:|---:|---:|---:|---:|
| 1 | 8,192 | 0.998664 | 1,129,609 | +4.05% |
| 1 | 32,768 | 1.010322 | 1,047,957 | -3.47% |
| 936 | 32,768 | 1.010322 | 1,047,957 | -3.47% |

The full-state driver initially confused the exposed orientation selector count
with the production orientation catalogue: `explicit_orientation_sum_only`
correctly exposes one summed selector. Its inventory assertion is now corrected
to check `selected_production_orientation_keys().len()` and report the selector
count separately. No production change was needed; the failed 527-second debug
run is preserved in the JSON. The corrected driver passed with optimized
libraries in 48.77 seconds, including full-state loading, geometry, rays and
reference. Its estimates exactly equal the selected-orientation values at the
same cube points. The JSON records exact orientation-set equality and hashes
of all 936 production keys.

All completed reference draws were finite. The two deterministic estimates pass
the stated 6% normalization and 8% moment acceptance tolerances. The reported
normalization dispersions are `0.14124` and `0.09902`: these apply the ordinary
per-draw standard-error formula to deterministic Halton values, and are **not
randomized-QMC confidence intervals**. They signal a broad weight distribution;
these short checks do not provide a precise normalization measurement or a
physical variance comparison.

## Eight H/Z rays

All eight archived raw points pass native inverse/forward checks independently
in f64, Quad and Arb, with both the selected and full production catalogues. Original binary64 inputs are promoted with
`T::from_f64_exact_binary`; all derived map arithmetic is rebuilt natively.
The JSON's H/Z values are independent binary64 oracle/archive checks, not
claims of direct native Esurface evaluation. Maximum roundtrip errors are
recorded per point. These direct native calls do not exercise the process's
precision-rescue stack.

The complete raw density, including ordinary complement coordinates, behaves
as follows (`GeV^-12`; two fixed H/Z directions):

| Normal-plane radius R [GeV] | q, Hplus/Zminus | q, Hplus/Zplus |
|---:|---:|---:|
| 0.2 | 3.661071e-36 | 3.661014e-36 |
| 0.02 | 1.190252e-35 | 1.190250e-35 |
| 0.002 | 3.796754e-35 | 3.796753e-35 |
| 0.0002 | 1.203934e-34 | 1.203934e-34 |

Using `q(R) ~ R^(-beta)`, the last decade gives `beta=0.50119026` and
`0.50119032`, consistent with the signed quadratic H profile's
`|H|^(-1/2)` focus. This does not imply bounded weights. If the actual physical
residual scales as `F~1/R`, then `F/q~R^(-1/2)` on these rays; local finite
variance follows only under that physical-scaling assumption. Direct H also
does not cover affine-star images. Physical raw-point replay, the actual rescue
stack, and complete all-orientation comparisons remain required.

## Reproduction

[gl638_x2_direct_h_driver.rs](gl638_x2_direct_h_driver.rs) uses existing public
state, map and reference owners; no production code is duplicated. It includes
the additional full determinant check. The final portable linker and corrected
driver were copied to fresh scratch directories and compiled against matching
current libraries. They passed the selected-orientation geometry-only gate
(11.85 seconds including compilation), then the full-catalogue 32,768-point
reference gate with optimized libraries. Both gates include native ray inverses
and the full determinant.

Copy the driver, [linker script](compile_gl638_x2_direct_h.sh), runtime card and
JSON to a scratch directory. Build the current API library normally, then set
`X2_TARGET` to its profile directory and `X2_API_FINGERPRINT` to the matching
built `lib-gammaloop_api.json`. The linker resolves exact dependency fingerprints,
not the newest filename among incompatible cached Rust artifacts. Then run:

```sh
./compile_gl638_x2_direct_h.sh
./gl638_x2_direct_h_driver STATE gl638_x2_direct_h.toml GL638_X2_DIRECT_H.json 32768 result 1
```

Use the all-orientation state and final argument `936` for the full-state
reference gate. A point count of `0` runs only geometry/native inverse and
finite-difference checks. The loaded state is read-only. The supplied Symbolica
license and ordinary gammaLoop build/runtime environment are required. If the
API was built with `ufo_support`, the matching Python shared-library directory
must be on `LD_LIBRARY_PATH`. The frozen optimized build linked Python 3.14;
prepending its `sysconfig.get_config_var("LIBDIR")` resolved an initial loader
error before any state or sample was processed.
