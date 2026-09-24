= GL638 assessment of E-surface sampling channels
<gl638-assessment-of-e-surface-sampling-channels>
Research only, 12 September 2026. No GammaLoop source, graph metadata, runtime card, state or integration workspace was modified, and no integrand integration was run. This note combines the existing measured results with new analytic work and lightweight arithmetic checks. The latter are in `gl638_geometry_checks.json` and their standalone script alongside this note.

== Conclusions for the requested design
<conclusions-for-the-requested-design>
+ A single radial map around #strong[H in cut-1 kinematics] can give finite local variance at the measured hard H/Z corner. A signed quadratic distance map, giving density proportional to `|H|^(-1/2)`, is sufficient under the measured regular codimension-two model. It does not bound the leading maximum weights.
+ H and Z occupy #strong[the same three-dimensional p block];, not two independent loop subspaces. To compensate their combined inverse-radius growth, use one joint chart with coordinates `(H,Z,tangent)`, with either a polar normal radius sampled uniformly or two signed-quadratic normal coordinates. Calling this an independent product of two three-dimensional radial E-surface maps would double-count that block.
+ A raw H channel alone does not cover an integrated CT singularity at `H(A-star)=Z(A-star)=0`. However, for native A in this GL638 graph, there is a particularly simple exact conditional map: #strong[at fixed s, the A-star projection is affine in p];, so H(A-star) is another ellipsoid in p. The API need only express the target frame and exact native variant/center to expose this useful special case. A full six-dimensional star-first chart remains a valid more general alternative.
+ Focusing the global radial distance to a Cutkosky shell is generally not a solution to LU-projected singularities: the physical cut-prepared point is invariant under global scaling of the original input. Targeting a cut\'s own equation after its own LU projection is identically zero. Similarly, A evaluated at its own A-star is identically zero. The target and its evaluation frame are indispensable parts of the channel.
+ This is useful new sampler functionality, but a channel must include the graph partition/energy equation, signed-cycle coordinate embedding, evaluation frame, conditional dependencies, normalized radial law, support/inverse branches and exact physical density. An edge-list-only shorthand is appropriate only when these can be resolved unambiguously.

== Physics and measured evidence being targeted
<physics-and-measured-evidence-being-targeted>
The published current prescription uses generation LMB `[3,4,7,10]`, where `K=(p0,p0+a0,s0,t0)` and the gluons are `q13=p0-t0`, `q14=s0-p0`. The native cut-1 parent is `[3,6,7,10]`; its two-loop active block `[3,7]` varies p and s while preserving a and t. Native one-loop `[7]` varies s. These statements and the current WH/WF metadata are recorded in #link("../gl638/README.md")[the tracked report];.

The six physical cut equations and the full threshold catalogue are recorded in #link("/common/dev/gl638_integration_validation/GL638_CURRENT_PHYSICS_ASSESSMENT.md:11")[the detailed assessment];. For physical cut 1, edges `(2,6,10)`, the LU rescaling is

```
Eh(tau*a0) + E(tau*(a0+t0)) + E(tau*t0) = Q,
a=tau*a0, t=tau*t0, p=tau*p0, s=tau*s0,
E(v)=sqrt(m_t^2+|v|^2), Eh(v)=sqrt(m_h^2+|v|^2).
```

Thus tau depends only on original a0,t0, not on p0 or s0. Q=1000 GeV, mt=173 GeV and mh=125 GeV. Write `B=E(a+t)+E(t)` in prepared coordinates. The relevant equations are

```
H = eta(2,4,12) = E(a+p)+E(p)-B,
P = eta(3,12)   = 2E(p)-Q,
Z = eta(3,10,13)= E(p)+E(t)+|p-t|-Q,
A = eta(7,8)   = 2E(s)-Q.
```

H is physical cut 3 evaluated on the cut-1 shell. It is not a normal threshold CT of cut 1: this Higgs-crossing singularity is removed through LU. Therefore a sampling target resolver must allow a graph-cut/E-surface geometry #strong[even when no corresponding CT is generated];. Restricting targets to the generated normal-threshold metadata registry would wrongly exclude H.

The multiplier `WH=H^2/[H^2+(P Z/Q)^2]` is evaluated at each relevant CT\'s own star. A and U have shared-one-loop weight `1-WH` and native-two-loop weight WH. F uses the analogous WF splitting; P, Z and V have unit native-two-loop weight. Sampling must leave that integrand prescription, common-center grouping, dual cancellation and local CT PV symmetry unchanged. The supplied #link("/common/dev/gammaloop_ir_safe_thresholds/TMP_TO_IGNORE/ttH_defo.pdf")[ttH\_defo.pdf];, especially its last page, motivates the complementary sector construction; it does not supply a proof that every Higgs-crossing surface is separated from every two-top surface. The independently established H/Z exception is exactly the target here.

The current complete sum over all 936 orientations, all six physical cuts and both local/integrated UV retains approximately `f=C/R`, where `R=sqrt(H^2+(P0 Z/Q)^2)` and P0 is separated from zero. Eight points on two rays range from 0.2 to 0.0002 GeV; current forced-Arb replays reproduce the previous values. These are hard momenta, not gluon-soft radii. The regular two-normal measure is `R dR dtheta`, giving a finite local absolute first moment and a logarithmic second moment under a smooth proposal, conditional on bounded coefficients on an angular/tangential neighborhood. See #link("../gl638/runtime-followup.typ")[measured powers and scope];.

A previously diagnosed integrated-A maximum had `Hstar=-0.90172` GeV and `Zstar=-12.78656` GeV, while sampled H and Z were `+13.46493` and `+27.06685` GeV. Thus a base H shell alone cannot be assumed to cover the CT image. See #link("../gl638/README.md")[maximum decomposition];. The more recent saved integration extrema pass their precision replays, but they have not all been decomposed into a proof that they are H/Z-dominated; avoid that inference.

== What “sphere” means in this actual graph
<what-sphere-means-in-this-actual-graph>
At fixed a,t, H=0 is an ellipsoid centered at `b=-a/2`. For a unit vector n, its positive radial root is exactly

```
rH(n)^2 = (B^2-|a|^2-4m_t^2) / [4(1-(a.n/B)^2)],
p=b+rH(n)*n.
```

The numerator is positive except at the degenerate pair `t=-a/2`, when the ellipsoid collapses. The denominator is positive since `B>|a|`. Rescaling each ray by rH(n) makes H=0 a unit sphere. This is a nonlinear coordinate change with a nontrivial Jacobian; it does not make the full momentum measure flat or make the E-surface energy equal the radial displacement everywhere. Locally, at a regular root, H is a smooth nonzero multiple of `r-rH`.

The same claim is #strong[false in an unnecessarily large active block];: native two-loop A(p,s)=2E(s)-Q is the noncompact cylinder `R_p^3 x S_s^2`, not a five-sphere. Its null directions must remain spectators, or the sampler must use a star-manifold chart that keeps those directions. This provides a direct GL638 reason to validate the active rank instead of equating an existing CT subspace with a compact radial sampling domain.

A new sampling center can be chosen for its own target (for H, b=-a/2). It need not equal the common CT center unless the target explicitly means a #strong[runtime CT-star pullback];. The two center roles should remain distinct in the API.

== Normalizable density and the one-H channel
<normalizable-density-and-the-one-h-channel>
For one regular codimension-one shell, a Jacobian proportional exactly to its distance would induce `q_delta~1/|delta|`; its normal integral diverges, even if the shell bounds a finite volume. Instead choose `0<beta<1`, for example beta=1/2:

```
q_delta(delta) = (1-beta)/(2 Delta) * (|delta|/Delta)^(-beta),
delta = random_sign * Delta * U^(1/(1-beta)),  |delta|<Delta.
```

The signed quadratic case has a Jacobian proportional to `sqrt(|delta|)`. Use it around the compactified H radial root, or within a finite shell neighborhood mixed with a normalized ordinary radial law. A complete map must cover both sides with their actual available lengths; near the center or domain boundary the two lengths need not agree.

For `f=C/R`, taking delta=H gives

```
f^2/q * dH dZ ~ R^(beta-1) |cos(theta)|^beta dR dtheta.
```

This is locally integrable; a typical leading weight still grows as `R^(beta-1)`. Thus a single correct H channel with beta=1/2 is an attractive first implementation and a finite-variance remedy for this regular corner, not a universal bound on maximum weights.

A useful declarative description is

```
parent=[3,6,7,10]
condition on complement=[6,7,10] in original coordinates
host_cut=[2,6,10], at=lu_cut
radial_block=[3], surface=[2,4,12], center=interior_auto,
distance_density_power=1/2
```

The actual syntax should be supplied by the API-design study. Here the complement includes original s0 and the pair (q6\_0,t0) determining a0; the p block is generated conditionally. Computing the map in prepared p and returning to original p adds tau^(-3) to its Jacobian. Equivalently, use the raw center and raw radial root divided by tau, and compute the ordinary raw spherical Jacobian directly. Complement dependence is block triangular.

The earlier scratch H proposal already passed forward/inverse and independent finite-difference checks. It changed the measured hard-ray second-moment exponent from about -1 to about -0.54 but its 256-draw full-graph pilot showed no global gain. This separates demonstrated local density behavior from unestablished global efficiency; see #link("../gl638/README.md")[the earlier proposal evidence];.

== Exact A-star targeting with another one-block ellipsoid
<exact-a-star-targeting-with-another-one-block-ellipsoid>
Let c=(cp,cs) be the #strong[actual selected native-two-loop common center];, and `RA=sqrt((Q/2)^2-m_t^2)`. At fixed original a0,t0,s0, tau, c and prepared s are fixed. Define v=s-cs and D=RA^2-|cs|^2\>0. The A-star cone coordinate is

```
ell = [cs.v + sqrt((cs.v)^2+D|v|^2)]/D,
sstar = cs+(s-cs)/ell,
pstar = cp+(p-cp)/ell.
```

The cancellation-safe expression for negative `cs.v` is `ell=|v|^2/[sqrt((cs.v)^2+D|v|^2)-cs.v]`. It obeys A(sstar)=0 and `ell>0` except the cone center.

#strong[Since ell is independent of p, H(Astar) is an affine pullback of the H ellipsoid.] Its prepared-p center is

```
b_eff=cp+ell*(-a/2-cp),
r_eff(n)=ell*rH(n).
```

A signed power map around this exact p-only radial root, with fixed s complement, targets the real runtime star image. In raw coordinates divide both b\_eff and r\_eff by tau. A change through pstar has physical Jacobian `(ell/tau)^3`; direct raw spherical coordinates contain that factor automatically. This avoids needing a separate six-dimensional star-first engine for this particular A target.

The input must identify host cut `(2,6,10)`, A threshold `(7,8)`, native variant active `[3,7]`, and the exact participating overlap center. Proposed semantic descriptor:

```
radial_block=[3], complement=[6,7,10], parent=[3,6,7,10]
surface=[2,4,12]
at=threshold_star(host=[2,6,10], threshold=[7,8],
                  variant_subspace=[3,7], center=each_containing_overlap)
distance_density_power=1/2
```

Each center becomes a normalized mixture branch. The current native group has only cut-1 surfaces; its inactive data and therefore its geometrical center depend only on a,t. This dependency is supported by the #link("TWO_NORMAL_PROPOSAL.typ")[existing center audit];. General graphs/groups must certify the analogous dependency or use the complete implicit Jacobian. An approximate center defines a normalized sampler if its own inverse is used, but does not place the singular density exactly on the runtime pullback. For an exact-alignment claim, the sampler and integrand should reuse the same prepared group-center context, or reproduce it identically from the complement. Independent SOCP recomputations that differ by numerical noise can misalign the density below that difference; reusing geometry also avoids duplicated solver cost.

New arithmetic checks using the recorded actual center tested 45 affine-map cases: component errors \<=1.43e-13 GeV, Hstar-shell residuals \<=2.28e-13 GeV. On those exact Hstar shells, base H ranges from -139 to +9927 GeV. No new integrand value was calculated. These checks substantiate the geometry, not global variance.

Along the #strong[whole native six-dimensional cone] `(p,s)=c+ell*(pstar-cp,sstar-cs)`, the projected star is independent of ell. Thus changing only that cone\'s overall radial sampling cannot concentrate on Hstar. A-star\'s own E-surface equation is zero by definition. This is why selecting a surface, its frame and the block to vary matters.

== Two-normal channel in the same p block
<two-normal-channel-in-the-same-p-block>
At the recorded hard intersection, the new arithmetic check gives `|grad H x grad Z|=1.73380575` and cosine 0.777009: these are regular independent normals in p-space. They are not independent #emph[three-dimensional subspaces];.

An exact conditional GL638 chart already has an analytic derivation. Let h=H, z=Z, u=E(p), S=B+h, D=Q+z-E(t). Then

```
2 a.p = S^2-2Su-|a|^2,
2 t.p = 2Du-D^2-m_t^2+|t|^2.
```

For `a x t !=0`, write `p_parallel=d+e*u`, and `p=p_parallel+w*n0`. With kappa=|e|^2-1, `uc=-(d.e)/kappa`, and `du^2=uc^2-(m_t^2+|d|^2)/kappa`, set

```
u=uc+du*cos(phi), w=sqrt(kappa)*du*sin(phi),
Jp=|d^3p/(dh dz dphi)|=u(S-u)(D-u)/(|a x t|sqrt(kappa)).
```

The full phi circle includes both inverse branches. The unsquared energy signs and du^2\>0 must hold; if arcs are restricted, their total length belongs in the PDF. The prior independent finite-difference determinant agreed within 3.8e-9 at a regular point. See #link("TWO_NORMAL_PROPOSAL.typ")[complete derivation, support and validation];.

Two appropriate normal laws are:

- Polar normal radius: `h=R cos(theta)`, `z=R sin(theta)/alpha`, with fixed positive alpha, uniform R and theta. It gives `q(h,z)~1/R` and compensates leading `C/R` weights.
- Product of two #strong[scalar normal-coordinate laws];: `q(h,z)~|h|^(-1/2)|z|^(-1/2)` on a certified finite rectangle. It also scales as 1/R and gives leading `f/q~C sqrt(|cos(theta)sin(theta)|)` up to smooth factors. Both marginals are normalizable. This is a product within one rank-two joint chart, not the user\'s product of separate loop blocks.

The polar law is isotropic in the physically scaled normal plane; the product law may be easier to compose with reusable one-dimensional density laws. Both require a certified feasible normal domain for each complement, or a deterministic ordinary-map fallback chosen from the complement alone and reproduced by inverse-density evaluation. Sampling then rejecting infeasible points without the acceptance normalization is incorrect.

The same (h,z,phi) chart can be used for pstar in the A pullback. With s held fixed, returning to p is affine and adds ell^3. Alternatively, a complete star-first A cone chart has Jacobian `ell^5 RA^2(RA-cs.n)`, followed by tau^(-6) for the two active momenta. #link("TWO_NORMAL_PROPOSAL.typ")[Exact star-cone formula];.

Native U has a genuinely coupled star equation `U=E(p)+E(s)+|s-p|-Q`. Its A-style p-only affine simplification is not automatic. There is nevertheless a known exact star-first chart: choose `|pstar|<RA`, direction n and

```
rho=Q*(Q-2E(pstar))/(2*(Q-E(pstar)+pstar.n)),
sstar=pstar+rho*n.
```

The full H ellipsoid fits strictly inside this pstar ball because of the Higgs mass gap. The coarea/cone Jacobian and boundary cusp need to be retained; #link("/common/dev/gl638_integration_validation/twenty-core/a-star-pullback/U_STAR_PROPOSAL.md:5")[the prior U-star audit] gives details. U-star coverage is motivated by WH appearing there; current maxima do not prove U dominance.

== LU and product-channel traps
<lu-and-product-channel-traps>
For this e+e- CM graph, globally scaling the original input by lambda changes tau to tau/lambda; prepared `tau*K` is unchanged. New arithmetic checks for lambda=.07,.3,2,17 reproduced the same physical cut-1 point within 1.14e-13 GeV. The actual owner #link("../../../crates/gammalooprs/src/integrands/process/cross_section/mod.rs")[rescaling all loop momenta] and #link("../../../crates/gammalooprs/src/integrands/process/cross_section/mod.rs")[preparing all cut points before group geometry] are consistent with this reasoning.

Consequently, placing more raw hyperradial samples at a cut\'s own shell may help match the auxiliary LU localization, but does not by itself change the sampled physical event shape or tackle the H/Z shape corner. A sampler targeting a threshold #emph[after] LU must either vary a cut-preserving block with fixed tau, as above, or explicitly pull back the threshold equation and use its real shape derivatives. It must not replace the physical projected residual by an unprojected one for convenience.

For a product `cut(...) x cut(...) x complement(...)`, validate exact signed-cycle span and conditional dependency, not merely nonoverlapping edge labels. Conditional sequential blocks can have a triangular Jacobian when every target\'s needed complement is already sampled. Mutual dependencies require a joint solve/chart, an explicit ordering with recomputed density, or rejection at definition time. Two constraints within the same p-space consume two scalar degrees plus a tangent, not two three-vector blocks.

== Validation matrix for an implementation proposal
<validation-matrix-for-an-implementation-proposal>
#figure(
  align(center)[#table(
    columns: 2,
    align: (auto,auto,),
    table.header([Check], [Why it matters for GL638],),
    table.hline(),
    [Independent finite-difference Cartesian determinant, including tau and ell factors], [Forward/inverse consistency alone can share the same wrong Jacobian, as the hyperspherical bug demonstrated.],
    [Analytic ellipsoid volume at fixed complements], [For H, transverse semiaxis is \`a\_perp=sqrt(B^2-],
    [Independent integration of q against a known reference proposal; test smooth Gaussian and asymmetric moments], [Tests normalization, both root sides, branch accounting and offsets beyond the tautology `q(F)J=1`. The whole momentum domain has infinite volume, so a constant full-space integrand is not a unit-volume check.],
    [HZ chart signs, both inverse branches, phi joining points, angular and tangent domains], [Protects against a squaring-induced extra root or silently omitted branch.],
    [Empty/pinched H fibers; a parallel to t; normal-rank degeneration; center boundaries], [Must enter declared deterministic ordinary-map fallback or a different chart with exactly evaluated PDF; never drop/resample unnormalized.],
    [Round-trip target residuals in raw, LU and each CT-star frame], [Detects using base H to target Hstar or treating own-cut/own-star zero as a varying normal.],
    [Many H/Z angles, including axes; multiple tangent positions and complements], [Existing two rays are evidence of power, not a neighborhood theorem.],
    [Hierarchical H/Z approaches and tangent-contact limits], [Control singular coefficients beyond fixed-ratio rays; the tangent and shrinking-F corners can have different measures.],
    [Soft13, soft14, double soft, shrinking-F and native-origin controls], [Retain the ordinary soft-aligned channels and test that focusing hard thresholds has not sacrificed those regions.],
    [Actual A and U overlap centers, center changes, star-projection ratios ell small/large], [Captures propagated CT singularities and center-dependent conditioning.],
    [Saved extrema plus independent seeds, equal sample and walltime budgets, Re/Im and absolute first/second moments], [A local density improvement need not improve global cost or maxima. Identify dominant terms rather than assuming every extreme is H/Z.],
    [Full six-cut sum, all 936 orientations, physical numerator and complete UV], [Preserve the previous final-prescription scope; per-cut or selected-orientation results are useful diagnostics only.],
    [Arb checks on nearby cancellation controls and exact-degeneracy reporting], [Importance sampling does not define the separate A/P coincident-residue limit and can deliberately visit more delicate points.],
  )]
  , kind: table
  )

The proposal has a solid mathematical route to the measured hard corner and explicit reusable GL638 charts. It should be developed in stages: single H normal in cut-preserving p, exact A-star H in the same conditional block, then a rank-two H/Z chart and U-star coverage. Normalization/correctness should be settled independently of finite-run performance. No claim of complete GL638 finite variance follows before angular, hierarchy, star and global-tail tests.
