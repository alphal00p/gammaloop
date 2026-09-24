= Surface-adapted sampling: literature and mathematical checks
<surface-adapted-sampling-literature-and-mathematical-checks>
A surface-adapted extension is well motivated and has a direct predecessor in Soper's numerical QCD program. A regular, compact, convex E-surface admits an explicit radial map to a sphere. The proposed #emph[exactly linear] vanishing Jacobian across an entire threshold shell, however, is not a normalizable importance-sampling map. Integrable fractional powers, finite-width profiles, or a density adapted to the joint transverse distance from an intersection resolve this distinction. Different active subspaces do not by themselves make several surface maps independent.

== Primary literature
<primary-literature>
Page numbers below are the numbers printed in the papers, which equal the one-based PDF page index for Soper 1999/2001. Kleiss--Pittau's cover is unnumbered, so printed page 1 is PDF page 2.

=== Soper, 1999
<soper-1999>
#link("https://arxiv.org/pdf/hep-ph/9910292")[Techniques for QCD calculations by numerical integration];, Sections VI.A--C, printed pp. 18--21, Eqs. (71)--(89), establishes normalized importance sampling of the cut sum, its second moment, and additive mixtures of normalized elementary densities. Equation (89) couples soft radius and angle relative to the scattering surface's tangent plane; purely radial enhancement misses the angular concentration. Section VI.B explicitly distinguishes generic soft directions from an increasingly thin angular region, with different powers in those limits. These results justify inspecting anisotropic limits, rather than inferring global integrability from a few fixed-angle rays.

=== Soper, 2001: the closest predecessor
<soper-2001-the-closest-predecessor>
#link("https://arxiv.org/pdf/hep-ph/0103262")[Choosing integration points for QCD calculations by numerical integration];:

- Printed pp. 7--11, Eq. (18): channels select complete edge-based loop routings; physical cut sums and sampling-channel cut sums remain independent.
- Printed pp. 16--17, Eqs. (25)--(35): massless one-loop ellipsoidal coordinates obey \[ A\_\\pm=(|\\ell|\\pm|\\ell-p|)/(2\\kappa),\\quad 2\\kappa=|p|, \\quad d^3\\ell=\\kappa^3(A\_+^2-A\_-^2)dA\_+dA\_-d\\phi. \] The target is (A\_+=S\_+).
- Printed p. 18, Eqs. (39)--(42): a normalized shell profile is \[ \\rho\'=\\frac{N}{(1-A\_-^2+\\tau)A\_+(|A\_+-S\_+|+\\lambda)}, \\quad\\lambda=(S\_+-1)^2/S\_+\>0. \]
- Printed pp. 19--21, Eqs. (49)--(54): the correlated soft-point component uses \[ \\rho\_t=\\frac{A}{(|\\delta|+\\beta S\_+\\omega^2)(|\\delta|+\\beta\\gamma S\_+\\omega)}, \\quad\\delta=A\_+-S\_+,\\quad\\omega=|A\_--S\_-|+|\\phi|/\\pi. \] Ordered conditional sampling supplies normalization. The paper explicitly improves the older tangent-plane ridge by placing it on the ellipsoid itself.
- Printed pp. 21--22, Eqs. (62)--(65): spherical-shell enhancement also depends on tangential approach to an intersecting surface.

These are one-loop massless constructions, not a theorem about arbitrary massive multi-loop or projected surfaces.

=== Other primary references
<other-primary-references>
#link("https://arxiv.org/pdf/hep-ph/9405257")[Kleiss and Pittau, Weight optimization in multichannel Monte Carlo];, printed pp. 1--3, Eqs. (1)--(7), defines nonnegative normalized channel densities, their normalized mixture, the common weight (f/\\sum\_i\\alpha\_iq\_i), and adaptation of channel probabilities using second-moment estimates. All component densities must be evaluated at a sampled physical point. The derivation supplies a useful implementation contract; it does not say that arbitrary unnormalized peak functions may be inserted as channel densities.

#link("https://arxiv.org/pdf/hep-ph/9806432")[Ohl, Vegas Revisited: Adaptive Monte Carlo Integration Beyond Factorization];, printed pp. 4--6, Eqs. (9)--(15), combines inverse coordinate maps, their Jacobians, and adaptive cube densities. Section 4.1 discusses the cost of cross-channel maps and density evaluation. Eqs. (25)--(26), printed pp. 9--10, provide an analytically normalized example combining spherical, cylindrical, and Cartesian peaks. Section 4.3 gives a cheaper construction based on fixed analytic partition functions with separate adaptation. This distinction is relevant if the existing GammaLoop estimator uses partitioned channel integrals rather than the fully adaptive mixture density.

#link("https://arxiv.org/pdf/hep-ph/9804454")[Soper, QCD calculations by numerical integration];, printed pp. 4--5, Eq. (9), keeps all cuts under one momentum-space integral to preserve pointwise cancellation. A density named after one cut must therefore remain a coordinate choice for the full cut-summed integrand.

== Independent derivation: when an E-surface can be made spherical
<independent-derivation-when-an-e-surface-can-be-made-spherical>
Fix every complement momentum and external datum. Suppose the active Cartesian variables are (x\\in\\mathbb R^d), and

\[ E(x)=\\sum\_e\\sqrt{|A\_ex+p\_e|^2+m\_e^2}-Q. \]

This is convex. If the stacked active routing map (x\\mapsto(A\_ex)#emph[e) has trivial kernel, the positive energy sum is coercive. If (\\min E\<0), its nonempty strict interior is bounded. Choose a center (c) with (E(c)\<0). Each ray meets (E=0) at one positive radius (r];\*(n)), for (n\\in S^{d-1}). Convexity and strict interior exclude a radial segment of boundary; differentiable regular points have positive outward radial derivative.

Define

\[ x=c+R,r\_\*(n)n,\\qquad R\\ge0. \]

The target becomes the unit sphere (R=1). The exact Cartesian volume element is

\[ d^dx=r\_\*(n)^dR^{d-1},dR,d\\Omega. \]

Derivatives of (r\_\*) with respect to angles do not add factors: their contributions point radially and disappear in the determinant after the radial column is included. Any unit-cube angular parametrization contributes its usual sphere Jacobian. This is a nonlinear map, not generally a linear transformation turning the original equation into a quadratic ellipsoid.

A chosen full routing and active cycle are essential. A surface can be cylindrical in a larger loop space because it is independent of spectator directions; factor those directions out. Empty surfaces, pinched zero-volume sets, and surfaces with nonsmooth massless points need explicit cases. A compact arbitrary set need not have spherical topology; here convexity is doing the work. A residual evaluated after a nonlinear LU or CT-star projection need not retain this convex structure. The construction must be checked for the actual pulled-back target.

== Independent derivation: the normalization obstruction
<independent-derivation-the-normalization-obstruction>
At a regular nonpinched point of the shell, write (\\delta=r-r\_\*). Radial offset, Euclidean normal distance, and the E-residual differ by finite nonzero local factors. If a proposed finite-interval map obeys

\[ \\left|\\frac{dr}{du}\\right|\\asymp|\\delta|, \]

its generated conditional density is (p\_r\\asymp1/|\\delta|). The integral over either side of the shell diverges logarithmically. Equivalently, integrating (du\\propto d\\delta/|\\delta|) requires an infinite coordinate range to reach the shell. Finite enclosed volume does not remove a divergence transverse to its boundary.

A valid singular shell profile uses (0\<\\beta\<1). On (\\delta\\in\[-a,b\]),

\[ p(\\delta)=\\frac{1-\\beta}{a^{1-\\beta}+b^{1-\\beta}}|\\delta|^{-\\beta}. \]

Choose the negative branch with probability (a^{1-\\beta}/(a^{1-\\beta}+b^{1-\\beta})), then set (\\delta=-aU^{1/(1-\\beta)}); use the analogous positive branch. For (\\beta=1/2), the radial map is quadratic and its Jacobian vanishes as (\\sqrt{|\\delta|}). Include the angular/radial-volume factors above in the physical density. A finite shell component can be mixed with an ordinary full-support channel to retain both tails.

A finite-width reciprocal profile (1/(|\\delta|+\\epsilon)) is also normalizable on a finite interval, but becomes smooth below (\\epsilon); it cannot prove an asymptotic variance cure. Log-corrected borderline kernels are mathematically possible, for example (1/\[|\\delta|\\log^{1+\\eta}(a/|\\delta|)\]) in an appropriately restricted interval, (\\eta\>0), but create extreme proximity to roots and are an unattractive first implementation.

== Independent derivation: one surface versus an intersection
<independent-derivation-one-surface-versus-an-intersection>
Let (u,v) be two regular transverse coordinates and let the remaining tangential variables be (y). Assume a bounded nonsingular change-of-variable Jacobian and

\[ f(u,v,y)\\sim C(\\theta,y)/R,\\qquad R=\\sqrt{u^2+v^2}. \]

With bounded (C), a smooth proposal gives a finite local first moment but second moment proportional to (\\int dR/R). A two-normal density (q\\propto1/R) is normalizable because the normal measure is (R,dR,d\\theta); its leading weights are bounded. The exact power forbidden across a codimension-one shell is therefore allowed around a codimension-two intersection.

A single shell profile (q\\propto|u|^{-\\beta}), (0\<\\beta\<1), already gives local second moment

\[ \\int R^{\\beta-1}dR\\int|\\cos\\theta|^\\beta d\\theta\<\\infty, \]

but weights still grow as (R^{\\beta-1}). Two valid shell profiles used as a #strong[product in genuine joint coordinates] are stronger:

\[ q\\propto|u|^{-\\beta\_1}|v|^{-\\beta\_2}. \]

Each exponent must be below one; their sum controls radial behavior. A positive sum makes the local second moment finite. A sum at least one bounds the leading radial weights. In particular, (\\beta\_1=\\beta\_2=1/2) gives (q\\propto1/\\sqrt{|uv|}), hence (q\\propto R^{-1}/\\sqrt{|\\cos\\theta\\sin\\theta|}), with bounded leading weights. An #strong[additive mixture] of two half-power shells does not have this product strength.

These are local sufficiency statements, not a proof that all angular, tangential, boundary, or other GL638 limits are controlled.

== Multiple-surface composition
<multiple-surface-composition>
A product of maps on disjoint Cartesian blocks has a product Jacobian when each block map depends only on that block and an already fixed complement. A triangular conditional construction can allow dependence on previously generated variables. Merely naming different subspaces does not establish either property: an energy in the first surface can depend on the second block's momentum through external data.

For overlapping targets, a local joint chart ((E\_1,E\_2,y)) exists where the target gradients have rank two, but its Jacobian and possible branches must be included. Tangencies, coincident surfaces, moving centers, and rank loss are separate cases. A map that spheres one target can displace the other target or leave its residual angular rather than radial. Consequently, a user-facing product notation should compile to a validated routing and dependency structure, and reject unsupported cycles of conditional dependencies rather than silently multiplying guessed Jacobians.

== Validation implications
<validation-implications>
Normalization of a generated density is necessary but not sufficient. Checking the mean of the implementation's own (qJ) can be tautological if both sides contain the same Jacobian mistake. Useful independent targets include a Cartesian cube of known volume, Gaussian moments, analytically integrable spheres/ellipsoids, finite-difference determinant checks, forward/inverse round trips, cross-channel density reciprocity, and a known normal-intersection integrand.

An especially transparent two-normal test on (\[-1,1\]^2) is

\[ f(u,v)=\\frac1{\\sqrt{u^2+v^2}},\\qquad \\int f,du,dv=8\\operatorname{asinh}(1). \]

The product proposal (q=1/(16\\sqrt{|uv|})) is exactly normalized and gives

\[ \\frac f q=\\frac{16\\sqrt{|uv|}}{\\sqrt{u^2+v^2}}\\le8\\sqrt2. \]

Uniform sampling has divergent second moment. The one-shell proposal (q=1/(8\\sqrt{|u|})) has finite second moment but unbounded weights. On a unit disk, (q=1/(2\\pi R)) makes the weight of (f=1/R) identically (2\\pi). These distinguish normalization, variance improvement, and bounded weights without production physics or numerical cancellation obscuring the outcome.

For GL638, the decisive question is whether the proposed targets vanish on the actual hard H/Z corner #emph[after] the relevant LU and CT-star maps, with a nonsingular joint-normal chart. If so, either one fractional shell or a joint two-normal construction can address the measured local power. Raw cut-edge lists alone cannot establish that premise. Importance sampling also does not repair undefined exact-coincident-root expressions or replace the shared CT subspace/center required for dual cancellation.
