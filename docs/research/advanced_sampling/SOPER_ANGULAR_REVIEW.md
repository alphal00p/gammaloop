# Soper comparison and generic angular control

Research note, 12 September 2026. No GammaLoop source changes.

## What Soper supplies beyond an isotropic shell map

[Soper (2001)](https://arxiv.org/pdf/hep-ph/0103262) supplies several complementary improvements:

| Location | Improvement |
|---|---|
| IV, V; (19)–(24) | Cut-conditioned virtual-loop sampling together with enhanced soft/collinear final-state sampling; channel labels remain independent of the physical cut sum. |
| VI; (35), (39)–(42) | Ellipsoidal Jacobians provide soft factors; angular endpoint enhancement strengthens near ellipsoid/sphere tangencies, while shell width follows deformation suppression. |
| VII; (49)–(54) | A correlated soft/threshold density resolves an anisotropic narrow region using sequentially normalized azimuth, tangent and normal distributions. |
| VIII; (62)–(65) | Sphere sampling couples angular separation and radial widths, including finite separation from the other surface. |
| IX; (66)–(70) | Separate soft-endpoint and collinear-region components cover exceptional self-energy configurations. |

The soft/threshold kernel contains

`1 / [(|delta| + a omega^2)(|delta| + b omega)]`,

where `delta` is a shell normal and `omega` measures distance within the shell to the soft point. Its powers reflect the quadratic vanishing of contour deformation. They are not universal powers for a real, threshold-subtracted integrand.

## Generic extension proposed for GammaLoop

The prior proposal already includes conditional contexts, joint surface normals, multiple patches and inverse density evaluation. It should additionally make **tangential targets and normal widths conditioned on tangential distance** explicit first-class ingredients. Uniform angular sampling, or an ordinary adaptive grid on unmodified angles, cannot guarantee good weights in a narrowing anisotropic region.

A generic local construction is:

1. Resolve a target surface and a distinguished soft, collinear, intersection or near-tangency locus from the actual graph energy equations and routing. These targets are geometric, so they apply to amplitudes and cross sections alike.
2. Choose a regular chart `(delta, z, y)`, with surface-normal coordinate `delta`, tangential coordinates `z`, and complement `y`. More than one normal is permitted. Rank loss requires another patch or a dedicated degenerate chart.
3. Sample the complement, then tangential coordinates, then the normal. The exact conditional density is `q(y) q(z|y) q(delta|z,y)` with the chart determinant. Angularly dependent widths do not generate an omitted Jacobian: the ordered conditional map is triangular, provided that this dependency order is obeyed.
4. Supply normalized normal profiles whose widths can depend on distances measured in `z`, finite gaps, and relevant scales. Exponents and width laws are configurable and can be tuned from the actual subtracted integrand.
5. Add full-support components and a declared fallback where a regular surface fiber disappears. Retain explicit forward/inverse branch accounting and the same raw-space density in the multichannel partition.

For a concrete generic mathematical family, let `omega=|z|`, take a finite normal interval `-D<delta<D`, and set

`A(omega)=a omega^p`, `B(omega)=b omega^q`,

`K(delta,omega)=1/[(|delta|+A)(|delta|+B)]`.

For `A != B`, its exact conditional normalization is

`I(omega)=2/(B-A) log[B(D+A)/(A(D+B))]`.

For `A=B`, use the continuous limit

`I(omega)=2[1/A-1/(D+A)]`.

Then `q(delta|z)=K/I` is normalized for every `omega>0`. Choosing the tangential density to include `I(omega)` recovers the combined correlated kernel rather than accidentally removing its enhancement by normalizing only the normal variable. The tangential measure must be included when deciding whether that marginal is normalizable. For the illustrative `p=2`, `q=1` case, `I(omega)` behaves as `log(1/omega)/omega`; in two tangential dimensions, the radial marginal is proportional to `log(1/omega) d omega`, which is locally integrable. This is an example of a topology-independent profile family, not an assumption that this particular profile is optimal for GL638.

Axes, scales and angular metrics can be derived from momentum routing and projected gradients; numerical one-dimensional CDFs or normalized analytic cap mixtures are possible generic building blocks. More complicated singular sets may need several local patches. Full automatic classification of every overlapping, pinched or hierarchical multiloop limit is a substantial separate problem. The generic API can expose these capabilities without promising that an initial automatic mode discovers the optimal one for every graph.

## Soft coverage independent of regular E-surfaces

Every loop-dependent massless edge defines a legitimate candidate soft target `q_e=0`, even if no nonpinched E-surface exists. A suitable routing can use that edge momentum as an active three-vector, with a complete complementary basis. A normalized Cartesian proposal with `q ~ 1/|q_e|` near zero compensates an isolated LTD `1/E_e` factor; in spherical coordinates this means the radial PDF behaves as `r`, not `1/r`.

Coverage of each isolated edge does not automatically cover all products when several edges become soft together. Compatible simultaneous soft sets need joint densities, respecting linear dependencies of their routed momenta. Their powers must be checked against the summed/subtracted integrand and the sampling measure. Mere inclusion of a particular LMB name is not a proof of the required density.

Thus `auto:surfaces` should either visibly expand to surface-plus-soft coverage, or be paired with an explicitly documented automatic coverage policy. The earlier literal meaning “surface channels only” should be revised if this expansion is adopted. Manual exact selection and generated-catalogue inspection remain valuable.
