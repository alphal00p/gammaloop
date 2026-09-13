# Matching Cutkosky sampling to the LU localization

## Result and scope

For the current multiplicative LU flow, a useful Cutkosky channel samples the
dimensionless root rescaling **t with density h(t)**, then sets the raw radial
coordinate to `r=R*(n)/t`. It should not sample t with h(t)/t. The latter factor
arises only when expressing a different coordinate measure. Matching h exactly
removes the auxiliary radial variation of a simple cut residue at fixed shape;
it does not remove physical angular, soft, threshold, or CT-star singularities.

An inexpensive, exactly invertible log-logistic proposal already approximates
the currently supported polynomial h profiles well in controlled one-dimensional
tests. An exact or tabulated CDF is possible too. Raised cut residues require
attention to derivatives of h; a positive derivative envelope or sufficiently
broad tails are preferable to assuming h alone bounds their weights.

This is a research proposal. The checks below integrate scalar radial models,
not GL638 or any generated GammaLoop integrand. It develops the observation in
[API_AND_GENERICITY_ADDENDUM.md](API_AND_GENERICITY_ADDENDUM.md), under
physical-cut/left/right composition. The distinct GL638 shape problem remains
as described in [gl638.md](gl638.md) and
[TWO_NORMAL_PROPOSAL.md](TWO_NORMAL_PROPOSAL.md).

## Implementation and physics boundary

The relevant implementation owners are:

| Owner | Fact used in this derivation |
|---|---|
| [utils/mod.rs](../../../crates/gammalooprs/src/utils/mod.rs), `h` and `h_dual` | LU profiles and their ordinary/hyperdual evaluations use the same family, sigma and power. |
| [settings/runtime.rs](../../../crates/gammalooprs/src/settings/runtime.rs), `HFunctionSettings` | The default is `poly_exponential`, sigma=1, power=None. |
| [settings/mod.rs](../../../crates/gammalooprs/src/settings/mod.rs), `lu_h_function` | The public TOML block is `[h_function]`; this is the auxiliary LU h. |
| [processes/cross_section.rs](../../../crates/gammalooprs/src/processes/cross_section.rs), `CutGroup` LU prefactor and `build_derivative_structure_atom` | The numerator includes `t*^(3L) h(t*)`; raised residues differentiate the complete residue structure. |
| [integrands/process/cross_section/mod.rs](../../../crates/gammalooprs/src/integrands/process/cross_section/mod.rs), cut-group evaluation | LU solves the cut at mapped raw momenta, evaluates h/h_dual, and keeps the sampling partition outside the residue derivatives. |
| [cff/esurface.rs](../../../crates/gammalooprs/src/cff/esurface.rs), `sampling_radial_map` | The current standalone physical-cut chart uses the graph energy equation and the zero center, with external data fixed. |

The normalized auxiliary h and the LU-flow Jacobian are part of the residue
representation, including for higher-order residues; see Eq. (2.13) of
Capatti, Hirschi and Ruijl, *Local Unitarity: cutting raised propagators and
localising renormalisation*.[^lu] The formulas below are an independent
specialization to the multiplicative flow implemented here. They do not assume
massless particles or homogeneous energy equations: external shifts and masses
remain fixed inside the energy equation.

Do not confuse the global LU `[h_function]` with a local threshold CT's
localization around its radial pole. Sampling changes a measure/partition,
not either integrand multiplier. In particular it must not modify the local
CT localization's required symmetry around r-star or its zero real PV integral.

## Correct measures and the exact map

Let `K` be all D=3L raw loop coordinates in the graph's full parent frame,
`K=r n`, `r>0`, `n` on the unit `(D-1)`-sphere. Assume the selected physical cut
has a unique regular positive root

```
eta_c(R_c(n) n) = 0,
a_c(n) = partial_R eta_c(R n)|_(R=R_c) != 0.
```

The LU solution is `t_c*=R_c/r`, and the physical point `Q_c=t_c*K=R_c n`
does not depend on r. Write the simple-cut contribution, omitting fixed phase
and normalization constants, as

```
F_c(K) = t_c^D h(t_c) A_c(Q_c) / |partial_t eta_c(t K)|
       = (R_c/r)^D h(R_c/r) A_c(Q_c) / (r |a_c|).
```

Since `d^D K=r^(D-1) dr dOmega`, changing from r to t gives

```
F_c(K) d^D K = [R_c^(D-1) A_c(Q_c)/|a_c|] h(t) dt dOmega.
```

Thus p_t=h cancels the entire radial dependence for this simple cut at a fixed
physical shape. This claim presumes that A_c is evaluated at Q_c and has no
additional dependence on the unscaled raw radius; any auxiliary runtime
dependence outside that factorization must be audited independently.

For a normalized positive proposal p_t, generate `t=G^-1(u)`, `r=R_c(n)/t`.
The exact conditional radial density and Cartesian determinant are

```
q_r(r|n) = p_t(R_c/r) R_c/r^2,
|dr/du| = R_c/[t^2 p_t(t)],
J_K = J_Omega R_c^D/[t^(D+1) p_t(t)]
    = J_Omega r^(D+1)/[R_c p_t(R_c/r)].
```

`J_Omega` is the actual determinant of the sphere coordinate map. Angular
derivatives of R_c add components along the radial column and cancel from the
determinant. The inverse is `n=K/|K|`, `t=R_c(n)/|K|`, `u=G(t)` plus the inverse
angular map. These identities supply both exact forward and inverse density.

The equivalent variables clarify several easy measure mistakes:

| Variable | Density induced by normalized h(t)dt |
|---|---|
| t, with measure dt | h(t) |
| y=log(t), with measure dy | exp(y) h(exp(y)) |
| ell=r/R_c=1/t, with measure d ell | h(1/ell)/ell^2 |
| signed raw distance delta=r-R_c | `R_c h(R_c/(R_c+delta))/(R_c+delta)^2`, delta>-R_c |

There is no singular density at delta=0 unless h itself has one. This proposal
matches the LU flow's auxiliary localization; it is a different objective from
an integrable fractional-power density across a surviving threshold shell.

For translated sampling centers, arbitrary nonlinear causal flows, multiple
radial roots, or a chart that rescales only part of the LU variables, the
identity `t*=R_c/r` cannot simply be assumed. Such a map needs its own exact
flow coordinate and determinant. Existing/pinched/absent root handling must
also stay explicit. A normalized ordinary fallback handles a fiber without a
regular root; it must be the same deterministic branch in forward and inverse.

## Profiles available in the current code

Write `z=t/sigma`. The polynomial families have

```
poly_exponential:            h(t) proportional z^(-p) exp(2-z^2-z^(-2)),
poly_left_right_exponential: h(t) proportional z^(-p) exp(2-z-z^(-1)).
```

Both suppress both endpoints faster than any power. The code currently
supports p=None/0,1,3,4,6,7,9,10,12,13,15,16; a proposal must mirror that
support rather than claim that every integer p is accepted by the integrand.
The normalization constants are, respectively,

```
N_poly = sigma exp(2) K_((1-p)/2)(2),
N_left_right = 2 sigma exp(2) K_(1-p)(2).
```

These follow by `y=z^2` or `y=z` from the standard integral representation of
the modified Bessel K function.[^bessel] Their associated transformed
densities are generalized inverse-Gaussian forms; no new special-function
dependency is necessary for an approximate sampling proposal.

For default `poly_exponential`, p=0, the normalized CDF has a useful exact form:

```
H(t) = 1/2 [erfc(1/z-z) - exp(4) erfc(z+1/z)].
1-H(t) = 1/2 [erfc(z-1/z) + exp(4) erfc(z+1/z)].
```

Differentiation gives `2 exp(2-z^2-z^-2)/(sqrt(pi) sigma)`. The inverse is
one-dimensional and monotone, so bracketed Newton/bisection is sufficient.
Use survival probabilities in the upper tail. In the small-z lower tail the
difference of erfc terms loses relative precision, so a tail expansion, scaled
erfc/log-CDF implementation or direct positive quadrature is needed. Clamping
CDF values to arbitrary interior endpoints would truncate support and is not
an acceptable exact implementation.

For `exponential`, `h(t)=2 exp(-(t/sigma)^2)/(sqrt(pi) sigma)` and
`H(t)=erf(t/sigma)`. Here h(0)>0: a proposal with
`p_t~t^(kappa-1)` has unbounded h/p_t for kappa>1 and infinite second moment
near zero for kappa>=2. A logistic exponent kappa<=1, a half-Cauchy proposal,
or a positive broad mixture floor avoids that problem.

`exponential_ct` takes an additional pole location in h/h_dual and is not a
normalized auxiliary density on positive t. The LU call passes no such
location, so automatic h matching must reject that family with an explicit
diagnostic. It must not invent a location or reuse a threshold-CT profile.

## A fast analytic approximation

An algebraically invertible log-logistic proposal is

```
t = s [u/(1-u)]^(1/kappa),
G(t) = (t/s)^kappa/[1+(t/s)^kappa],
p_t(t) = kappa (t/s)^kappa/[t (1+(t/s)^kappa)^2],
dt/du = t/[kappa u(1-u)].
```

It is strictly positive and normalized on t>0, has heavier tails than either
polynomial h family, and has a simple exact inverse and determinant. Its
formula can be built in the same eager Symbolica evaluator as the map and
proxy. Evaluate its tails in logarithmic form to avoid intermediate overflow.

A generic parameter choice matches the mode and curvature of the *log-t*
density, which includes the extra factor t. For
`h(t) proportional z^-p exp(2-z^a-z^-a)` with a=2 or a=1, set

```
y0 = asinh((1-p)/(2a))/a,
s = sigma exp(y0),
kappa = 2a sqrt(cosh(a y0)).
```

The default profile gives `s/sigma=1.1317139...`, kappa=4.06108637. A
lognormal approximation with mean log(s) and log-variance
`1/[2a^2 cosh(a y0)]` is also valid and has heavy enough tails for the
polynomial h families. The log-logistic construction has a simpler exact
algebraic inverse, which is attractive for the current evaluator design.

At fixed shape, for a normalized simple-cut radial integrand h, the mean is
one, the second moment is `integral h(t)^2/p_t(t) dt`, and the weight is h/p_t.
The following checks use sigma=1 and p=0:

| Proposal | Second moment | Maximum h/p_t found |
|---|---:|---:|
| Exact h | 1.000000 | 1.000000 |
| Matched log-logistic | 1.098670 | 1.183053 |
| Rational raw radius, R_c/beta=0.1 | 11.788964 | 16.794966 |
| Rational raw radius, R_c/beta=1 | 3.530859 | 4.822672 |
| Rational raw radius, R_c/beta=10 | 9.814200 | 13.682222 |
| Rational raw radius, R_c/beta=100 | 81.534049 | 115.108787 |

The rational comparator is exactly `r=beta*u/(1-u)`, which induces
`p_t(t)=a/(a+t)^2`, a=R_c/beta. It is a controlled ordinary radial proposal,
not a measurement of the current advanced surface-power chart or VEGAS grid.
Across all polynomial powers supported in the code, the matched log-logistic
second moment ranged from 1.060442 to 1.100891; the largest observed h/p_t
was 1.188258. Heavier analytic tails ensure bounded ratios at both endpoints;
the reported interior maxima were searched numerically, not certified bounds.

The portable [script](lu_h_profile_checks.py) uses only Python's standard
library. Its [JSON output](lu_h_profile_checks.json) records all families,
parameters, maxima and 8000/16000-interval Simpson refinement differences.
Both grids integrate in log(t) on [-7,7]. This interval safely resolves these
particular exponentially suppressed profiles; it is not a general tail-bound
algorithm. Refinement differences were below 1.2e-13. A separate symmetric
finite-difference check of the p=0 exact CDF at five points gave a maximum
relative derivative discrepancy 4.3e-9.

```
python3 docs/research/advanced_sampling/lu_h_profile_checks.py
```

## Numerical CDFs with an exact sampling weight

There are two distinct valid approaches:

1. Evaluate the exact normalized CDF and solve its inverse to a controlled
   tolerance. Differentiate the defining relation, `dt/du=1/p_t(t)`, and bound
   the inversion residual. For p=0, the erfc expression supplies the CDF;
   general profiles can use positive quadrature and a cached lookup bracket.
2. Define an approximate density whose normalization, CDF and inverse are exact
   for the approximation. For example, approximate the *log-t density* by
   positive piecewise exponentials, with integrable exponential tails on the
   entire real log-t axis. Every segment has an elementary integral and a
   log1p inverse. Normalize the sum of those exact segment masses, and use the
   resulting piecewise density in the Jacobian and partition.

The second approach affects variance, not the integral, provided its actual
map determinant is used. Never sample from an interpolated quantile while
claiming the determinant is `1/h` if that quantile's derivative differs from
`1/h`. Monotonicity, positive segment masses, consistent branch lookup and
full tail coverage are correctness requirements. Build/certify the profile at
warmup per h setting and derivative order, not with a new quadrature per sample.

## Raised residues and cancellation

For an m-th-order cut pole, use `y=t r` in the residue before performing the
radial change of variables. The change of residue variable contributes 1/r,
the LU prefactor contributes r^-D, and each derivative hitting h(y/r)
contributes another 1/r. Consequently the radial measure has the form

```
F_c(K) d^D K = sum_(j=0)^(m-1) C_j(n) t^j h^(j)(t) dt dOmega,
```

where the C_j include all physical R_c factors, energy derivatives and physical
amplitude derivatives; they do not depend on the auxiliary raw radius under
the factorization assumptions stated earlier. Matching h alone then leaves
`t^j h^(j)/h`. For the polynomial profiles these are growing Laurent
polynomials in t/sigma: moments may be finite while pointwise tail weights are
unbounded.

A generic nonnegative envelope is

```
p_t(t) proportional sum_(j=0)^(m-1) c_j |t^j h^(j)(t)|,
c_0>0, c_j>0 for every derivative to be covered.
```

Use a normalized positive approximation or exact piecewise CDF of this
envelope. Derivative zeros cause cusps in the envelope but no zero of the whole
density since h>0. A smooth sum of squares can replace the absolute-value
envelope if its tails dominate the same terms. Automatically inspect the
generated maximum pole order, but allow proposal coefficients to be steered;
no universal choice of c_j captures unknown angular amplitude coefficients
optimally. The heavier log-logistic proposal also dominates every finite-order
derivative at both endpoints of the polynomial families, though it can have
larger interior weights than a tuned derivative envelope.

Provided endpoint terms vanish,
`integral t^j h^(j)(t) dt=(-1)^j j!`. This is a useful independent analytic
raised-residue acceptance oracle. The sampling determinant and partition stay
outside the residue/hyperdual differentiation. Moving them inside changes the
physical derivative and is incorrect. Positivity of a proposal does not remove
signed cancellation among the C_j, between cuts, or between CTs.

## Multiple cuts and GL638

One chosen host cut defines the coordinate map; it does not select the physical
integrand contribution. Evaluate the complete cut/CT sum at that one raw point.
For another cut d at the same direction,

```
t_d = (R_d/R_c) t_c = a_d t_c,
F_d d^D K proportional a_d h(a_d t_c) dt_c dOmega
```

for a simple residue, with its own shape coefficient. A profile matched only
to c can be poor for d when the root ratio is large. A canonical channel per
useful distinct physical surface, using the common raw-frame exact density
sum, gives the appropriate mixture. Degenerate cut IDs with identical geometry
need not duplicate channels. Profiles for raised groups must cover their
largest relevant derivative order. A more elaborate conditional radial law
can mix the scaled h functions at fixed direction, but is unnecessary before
measuring whether the existing canonical multichannel mixture suffices.

For angular density p_Omega with respect to dOmega, the actual Cartesian score
is `q_c(K)=p_Omega(n) p_tc(t_c) t_c/r^D`; uniform sphere sampling has
`p_Omega=1/area(S^(D-1))`. Comparing only the radial h or p_t values omits both
the geometric factor and any channel-dependent angular density. Forward and
inverse evaluation must reproduce all root, angular and fallback branches.

GL638's six cuts can benefit from reduced auxiliary radial variation and fewer
extreme raw-to-physical scale ratios. This may reduce wasted evaluations and
stability-rescue cost. It is not a cure for the measured H/Z corner:
`t_c K=R_c n` is invariant as the raw global radius changes. The hard H/Z and
CT-star loci vary in shape/conditional momentum coordinates. The conditional
H, exact A-star H pullback and joint H/Z constructions in the existing studies
remain necessary for their established asymptotic density improvements.
No GL638 maximum-weight or integration improvement follows from the 1D table.

## Full support is not enough for a useful Gaussian acceptance test

Exact h matching has a subtle limitation for the proposed universal reference
harness. Both polynomial families and the `exponential` profile give a raw
density exponentially tiny as r tends to zero. A Gaussian reference has a
nonzero value there, so
`integral g(K)^2/q(K) dK` diverges near zero under that single exact-h channel.
The map is normalized and invertible, but a naive Gaussian Monte Carlo
acceptance run has infinite variance. Some tail directions can also need
coverage independently of the near-zero issue.

For a log-logistic exponent kappa, `q_r~r^(kappa-1)` near zero, hence a smooth
nonzero reference has finite near-origin second moment only if kappa<2D.
Default p=0 satisfies this in D>=3; some high supported powers in D=3 do not.
The lognormal approximation also has infinite Gaussian-reference variance near
the raw origin, for every finite log-variance: its density decays faster than
any power there. It needs the same broad coverage for this acceptance test.

A small normalized broad floor resolves this and protects the physical
integrand's less-focused regions. For example, mix with
`g_t(t)=s/(s+t)^2`. Its raw radial density has a nonzero finite limit at zero
and a 1/r^2 tail. Define the mixture CDF and invert it monotonically, or include
an ordinary canonical channel in the common exact density sum with a positive
sampling probability. If a profile mixture is used, include its entire density
in the determinant; do not add a weight floor only in the partition denominator.
The floor is also a convenient solution to the `exponential` family's lower
endpoint mismatch. The reference harness should report/enforce its moment
assumptions; it must not mistake a slow infinite-variance reference estimate
for a biased map.

## Proposed runtime controls and acceptance gates

Extend the existing named channel's radial-profile specification rather than
add another channel axis. A possible minimal public form is:

```toml
[sampling.channel_definitions.GL638.born_cut]
around = "phase_space(cut(2,6,10))"
parent_lmb = [3,4,7,10]
subspace_lmb = [3,4,7,10]
radial_profile = { kind = "lu_h", approximation = "log_logistic", broad_fraction = 0.02 }
```

The shorthand `radial_profile="lu_h"` can use the same defaults; both forms
resolve to one profile owner. Amplitude-only surface channels cannot inherit
an LU h and must reject this profile with a clear error.
The full parent and subspace must equal the actual graph's ordered parent LMB
for the current standalone physical-cut implementation. These example edge
lists do not enable an alternative native-parent chart automatically.

This is proposed syntax, not an implemented setting. `kind="lu_h"` reads the
actual integrand `[h_function]` settings by default; changing sigma/power there
must invalidate the warmup proposal cache. Expert overrides change only the
proposal, for example an explicit `scale`, `shape`, `derivative_order` or
positive derivative-envelope coefficients. An `exact_cdf` or
`piecewise_exponential` approximation is an alternative behind the same
profile owner. Broad-fraction changes require the corresponding exact mixture
CDF/inverse. Auto mode may select a profile from pole order and dimension after
the moment checks above; users retain normal named-channel selection and can
choose only these channels or combine them with existing LMB channels.

The decisive tests are:

1. Forward/inverse round trips across both tails; positive Jacobian and actual
   independent Cartesian finite-difference determinant in several dimensions,
   masses, external boosts and direction-dependent roots.
2. At fixed n, compare the production LU solve to `R_c/r`. Use an actual generated
   simple-cut graph and show its full weight is independent of t under exact
   h matching; compare the analytic physical shell measure, not a hand-copied
   implementation of the same Jacobian.
3. For every supported h setting, compare h and h_dual values/derivatives,
   CDF normalization and inverse derivatives. Test raised residues against the
   integration-by-parts oracle and the generated residue engine.
4. Run Gaussian and asymmetric-moment acceptance through saved generated
   states, with a broad density component where required for finite variance.
   Exercise exact and proxy partitions at a common raw point, all canonical
   channels, rotations and precision rescues.
5. Compare GL638 with equal work budgets, all orientations and all physical
   cuts/UV/CT terms: signed integral, absolute first and second moments, maximum
   weights, root ratios, rescue rate and wall time. Repeat saved H/Z and CT-star
   approaches to verify that unchanged shape singularities are diagnosed as
   unchanged, then combine with their dedicated conditional channels.

## Sources

[^lu]: Z. Capatti, V. Hirschi and B. Ruijl, *Local Unitarity: cutting raised
    propagators and localising renormalisation*, JHEP 10 (2022) 120,
    [arXiv:2203.11038](https://arxiv.org/pdf/2203.11038), Sec. 2.1.1,
    Eq. (2.13). The independent specialization and proposal calculations in
    this note are not claims taken from that paper.
[^bessel]: NIST Digital Library of Mathematical Functions,
    [Eq. 10.32.10](https://dlmf.nist.gov/10.32.E10), integral representation of
    the modified Bessel K function.

## Bounded X1 implementation decision (2026-09-13)

This section specifies the next implementation slice; these controls are not
implemented yet. X1 changes the auxiliary radial proposal for the existing
full-parent, zero-centered `phase_space(cut(...))` map. Conditional side maps,
threshold-distance focusing certificates, exact h CDFs, and automatic channel
discovery remain separate work.

Use the simple `radial_profile="lu_h"` spelling and an equivalent table with
`kind="lu_h"`, `broad_fraction` (default 0.02), and optional positive `scale`
and `shape` overrides. The scale and shape affect the proposal alone. Omission
of `radial_profile` retains the existing signed-distance power profile. A LU-h
profile replaces that power transformation entirely; it must not inherit its
`power != 1 && r == R` seam rejection, since r=R is a regular point here.

The focused component is the fitted log-logistic law above for the two
polynomial families. For `exponential`, start with the endpoint-safe rational
law (log-logistic shape 1) and scale sigma/2. Among rational laws this minimizes
the simple half-Gaussian h's radial second moment: it is
`sqrt(2/pi)+2/pi`, approximately 1.4345. Keep the supported polynomial powers
identical to h/h_dual; `exponential` ignores `power`, as the actual h does.
Reject `exponential_ct`, invalid sigma, or unsupported polynomial powers at
warmup with a diagnostic rather than reaching h's panic. Arbitrary algebraic
h families are not part of this claim.

### A broad floor independent of the physical root

Improve the earlier fixed-t broad example by using the existing raw scale
`beta=e_cm*b`. For each direction and physical root R, define a=R/beta and

```
p(t|R) = (1-epsilon) p_loglogistic(t;s,kappa) + epsilon a/(a+t)^2,
G(t|R) = (1-epsilon) (t/s)^kappa/(1+(t/s)^kappa) + epsilon t/(a+t).
```

The broad component then induces exactly
`q_r,broad(r)=beta/(beta+r)^2`, independently of R. For epsilon>0 the raw
proposal has an ordinary radial floor even when a cut root is very large or
small. This gives stronger reference-function coverage than fixing the broad
scale in t. There is no additional user scale knob in X1. A certified absent
or pinched surface uses this same normalized raw law with power 1, in both
forward and inverse, rather than inheriting a focused power. Numerical failure
to establish a root remains an error, not an absent classification.

Invert G monotonically using the existing safeguarded scalar owner. The two
analytic component quantiles give an immediate bracket:

```
Q_f(u)=s (u/(1-u))^(1/kappa), Q_b(u)=a u/(1-u),
min(Q_f,Q_b) <= G^-1(u) <= max(Q_f,Q_b).
```

The epsilon=0/1 cases use their analytic quantiles directly; epsilon=1 is
exactly the ordinary raw law and need not solve R. Evaluate the lower CDF for u<=1/2 and the survival probability for u>1/2, using positive reciprocal
forms in the tails. Scale residuals by u or 1-u; a fixed absolute CDF tolerance
must not erase a small tail. Inverse evaluation at a raw point is direct:
`t=R/r`, `u=G(t|R)`. Preserve native arithmetic and typed numerical retries.

Do not split the input interval into a focused branch and a broad branch and
return only that branch's determinant. Both branches cover the entire positive
t axis, so that would introduce two preimages and the wrong estimator under
the existing single-map contract. The monotone mixture CDF avoids new branch
accounting and avoids another canonical channel domain.

### Existing-owner implementation boundary

Compile the CDF, its stable tail forms, and r=R/t through the existing
`SamplingExpressionEvaluator` during the same warmup epoch as proxy programs.
Use its dual derivative `partial_t G=p(t|R)` and implicit differentiation:
`dt/du=1/p`, `|dr/du|=R/(t^2 p)`. The Cartesian determinant is still
`J_Omega R^D/(t^(D+1) p)`. Do not differentiate the scalar root iterations.
Angular derivatives of R remain radial-column contributions and cancel from
this determinant. Test any eager conditional/tail expressions in all native
precisions; do not evaluate an overflowing inactive expression eagerly.

Extend the existing setup's compiled-program cache to hold these neutral
profile programs alongside proxies. Native bindings clone evaluator buffers;
they do not reparse strings, refit parameters in binary64, or compile programs
again. Build the default fit algebraically from the h input parameters in the
symbolic expression, then evaluate it natively. One catalogue remains the
source of IDs. A typed radial-profile option containing floats will require
removing now-inapplicable Eq derives from the settings/catalogue records.

Pass the actual runtime `lu_h_function` to the fresh compile/warmup owner;
never substitute defaults because its current argument list lacks h. The
cross-section binder must aggregate `max_occurence` over every active,
geometrically equivalent cut group before deduplicating its root map. Preserve
its existing incompatible-external-shift diagnostic. Store the actual maximum
order with the bound profile for inspection and tests. The initial fitted
proposal does not invent order-dependent tuning: its tails cover every finite
h derivative packet for the supported families, but large raised-order interior
weights remain possible.

Keep the cut geometry registration independent of a particular named profile.
After cloning a registered implicit map, the canonical channel compiler attaches
that named channel's profile. Thus two names may target the same cut with
narrower/wider proposals without overwriting each other or introducing a second
target resolver. Physical LU evaluation still sums every active cut and CT;
all sampling factors remain outside h_dual/residue differentiation.

### Required X1 gates

- Native profile normalization, both-tail round trips, positive J and selected
  factor, analytic-quantile endpoints, and mixture inversion residuals.
- Independent Cartesian finite differences for the actual bubble cut and a
  direction-dependent massive graph; compare eager/dual derivatives directly.
- On a generated simple cut, verify at fixed shape that the complete weighted
  result is proportional to h(t)/p(t|R). Approximate matching does not predict
  a constant weight. Vary inherited h family, sigma, and every supported power.
- Use `extract_t_derivatives` to undo HyperDual factorial normalization when
  checking `integral t^j h^(j)(t) dt = (-1)^j j! integral h(t) dt`; separately
  check h normalization to the accuracy of its current tabulated constants.
  Several powers use binary64 normalization constants even in Arb, so exact
  Arb-epsilon agreement with one would be a false fixture expectation. Retain
  a generated raised-cut comparison using the existing triple-dotted bubble.
- Saved-state Gaussian/moment acceptance with the broad floor, common-frame
  multiple-cut density checks, warmup invalidation, and two names sharing one
  cut geometry but carrying different profiles.
- Only then run a fixed-budget GL638 A/B pilot. This can improve auxiliary radial
  variance or rescue cost; it does not change the H/Z shape singularity and
  does not establish bounded GL638 weights.
