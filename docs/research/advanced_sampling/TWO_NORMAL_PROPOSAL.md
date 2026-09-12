# Exact two-normal H/Z proposal and its A-star pullback

12 September 2026. Independent analytic and arithmetic audit. This is a proposal,
not an implemented or integrated sampler. It strengthens the one-normal map:
on a regular compact patch, a density proportional to `1/R` bounds a leading
`1/R` weight, whereas `|H|^(-beta)`, `0<beta<1`, gives finite local variance but
leaves the weight unbounded on generic rays.

## Exact conditional momentum chart

Work first in prepared physical cut1 coordinates, with fixed a,t, top mass m,
`Et=E(t)`, `B=E(a+t)+Et`, and

```
h = H = E(a+p)+E(p)-B,
z = Z = E(p)+Et+|p-t|-Q,
u = E(p), S=B+h, D=Q+z-Et.
```

The original unsquared energy signs require
`u>=m`, `S-u>=m`, and `D-u>=0`. Squaring with those conditions gives

```
2 a.p = S^2 - 2 S u - |a|^2,
2 t.p = 2 D u - D^2 - m^2 + |t|^2.
```

For `a x t != 0`, solve these two linear equations in their plane:
`p_parallel=d+e*u`. An orthonormal frame built from a and the part of t
perpendicular to a avoids an unnecessary ill-conditioned Gram-matrix inversion.
Let n0 be the unit normal to that plane. Then

```
p = d + e*u + w*n0,
w^2 = u^2-m^2-|d+e*u|^2.
```

The inverse Jacobian follows directly from the triple product of the gradients:

```
|d^3p / (dh dz du)| = u(S-u)(D-u) / (|a x t| |w|).
```

Both signs of w must be included. A simpler complete parameterization removes
the apparent endpoint singularity. For a physical nonempty H ellipsoid,
`S>sqrt(|a|^2+4m^2)>|a|`; hence `kappa=|e|^2-1>0`. Set

```
uc = -(d.e)/kappa,
delta_u^2 = uc^2 - (m^2+|d|^2)/kappa,
u = uc + delta_u*cos(phi),
w = sqrt(kappa)*delta_u*sin(phi),  0<=phi<2*pi.
```

A full circle includes both inverse branches, including their joining points.
The remaining Jacobian is smooth at `w=0`:

```
Jp = |d^3p / (dh dz dphi)|
   = u(S-u)(D-u) / (|a x t| sqrt(kappa)).
```

Require `delta_u^2>0` and the physical energy signs throughout the sampled
circle. If the interval `[uc-delta_u,uc+delta_u]` must be clipped to
`[m,min(S-m,D)]`, its endpoints determine the allowed phi arcs explicitly via
arccos. Uniform sampling on those arcs uses their exact total length in its
density; it does not reject points or retain only one square-root branch.
The initial compact-patch design should use a complete physical circle with
strict margins, so this additional bookkeeping is unnecessary.

## Normal-plane sampling and normalization

Choose a fixed positive scale `alpha=|P0|/Q` within a patch, with P0 separated
from zero, and use

```
h=R*cos(theta), z=R*sin(theta)/alpha,
|d(h,z)/d(R,theta)|=R/alpha.
```

Uniform R in `(0,Rmax)`, theta in `[0,2*pi)`, and phi in `[0,2*pi)` induce

```
q_p(p | a,t) = alpha*|a x t|*sqrt(kappa)
              / [4*pi^2*R*Rmax*u(S-u)(D-u)].
```

All factors, including Rmax and any restricted angular-domain lengths, belong
in the inverse density. This proposal has `q_p~1/R` on a regular patch. If the
complete integrand is `C(theta,phi,y)/R+O(1)` with bounded coefficient and
smooth nonzero remaining densities/Jacobians, its leading importance weight
is bounded there. This does not prove bounded weights or finite variance in
other regions or near tangencies.

Not every complement a,t admits H=Z=0. A normalized implementation can avoid
both unknown rejection probabilities and missing support in either of two ways:

1. Choose a normalized compact complement patch around a known regular
   intersection. Certify a disk `|h|<=Rmax`, `|z|<=Rmax/alpha` throughout that
   complement support, with strict margins on `|a x t|`, kappa, delta_u squared,
   and the original energy signs. Mix this local proposal with a globally
   supported ordinary proposal. This only addresses the certified patch.
2. Draw complements y from any normalized g(y). Decide from y alone whether
   an admissible disk can be certified. If so, use its explicitly calculated
   Rmax(y); otherwise use a normalized ordinary conditional momentum law.
   The conditional distribution integrates to one in either case. Recompute
   the same deterministic decision when evaluating its inverse density, and
   include the full outer mixture density in every weight.

The explicit polynomial/rational coefficients permit interval bounds for such
a disk. Adaptive reduction of Rmax is acceptable only when the final domain
has a certificate; testing a finite grid does not establish domain coverage.
Computing the entire global feasible H/Z domain is unnecessary. An outer
mixture supplies positive density outside a compact focused chart.

## Targeting the actual A-star image

The same conditional p chart can be used for p-star rather than sampled p.
For the A=(7,8) counterterm, let `RA=sqrt((Q/2)^2-m^2)`, sample
`s_star=RA*n` uniformly on its sphere, and use the same center as its native
two-loop subtraction group:

```
(p_cut,s_cut) = c + ell*((p_star,s_star)-c),  ell>0.
```

For fixed complements and `|cs|<RA`, this covers every `s_cut!=cs` exactly once.
Writing `v=s_cut-cs` and `D_A=RA^2-|cs|^2`, its inverse is

```
ell = [cs.v + sqrt((cs.v)^2+D_A*|v|^2)]/D_A,
(p_star,s_star) = c + ((p_cut,s_cut)-c)/ell.
```

For `cs.v<0`, rationalize the ell expression as
`|v|^2/[sqrt((cs.v)^2+D_A*|v|^2)-cs.v]`. Its six-dimensional Jacobian is

```
J_A = ell^5 * RA^2 * (RA-cs.n) * d^3p_star d(ell) dOmega.
```

Substituting the two-normal p-star chart gives the positive Jacobian factor

```
ell^5 * RA^2 * (RA-cs.n) * (R/alpha) * Jp.
```

With original, unrescaled a,t drawn from a normalized six-dimensional density,
cut1's tau depends only on those complements. Returning the two active
momenta to original coordinates adds **tau^(-6)** to this Jacobian. The density
is the product of the normalized complement, ell, sphere, R, theta and phi laws,
divided by the complete Jacobian. Dependence of tau and c on the complements
is block triangular and adds no other determinant. Drawing prepared a,t freely
instead would require a separate shell/radial parameterization and its Jacobian.

The fixed-center premise is supported for this full GL638 native `[3,7]` group:
its six surfaces are from cut1, and the SOCP's inactive energies/spatial shifts
depend only on the cut-preserving complements. Its objective is unweighted
minimax. If the origin shortcut does not apply, Z and F fix the ideal unique
center `(p,s)=(t,t)`; the actual recorded numerical solver center differs by
about 1e-6 GeV at the examined maximum. Subsets, other metadata or additional
cuts in a group can change this conclusion.

Any chosen fixed approximate center still defines a normalized sampling chart
when its own inverse and Jacobian are used. It does **not** place the singular
density on the exact runtime CT image, and therefore cannot establish its
asymptotic variance regularization. Use the actual common center for that claim,
cover every relevant A-containing overlap center with a normalized mixture,
and maintain positive conditioning margins. The existing full graph and its
threshold-group operations must remain intact.

## Independent checks performed

No GammaLoop execution was needed for this audit. On the existing prepared
hard-max complements, h=0.1 GeV and z=0.2 GeV give
`u in [263.0012,414.4574] GeV`, below the physical upper bound 555.0715 GeV,
with kappa=12.9427. Four circle angles, including the joining point phi=pi,
reproduce the unsquared H/Z equations within 9e-13 GeV. A direct central
finite-difference determinant agrees with Jp within 3.8e-9 relative.

Separately, the six-dimensional A-star cone map with nonzero cp and cs gives
the claimed determinant within 1.5e-9 relative; its inverse recovers ell within
2.3e-16 and p-star within 2.9e-14 GeV. These checks substantiate the formulas
at regular points, not a global domain certificate or an implementation claim.
