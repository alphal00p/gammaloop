= Massive planar double-box geometry check
<massive-planar-double-box-geometry-check>
Bounded read-only study, 13 September 2026. No production code changed. The existing debug CLI was used for graph generation/export only; this is not validation of the latest sampling implementation or a variance measurement.

== Exact parent frame and surfaces
<exact-parent-frame-and-surfaces>
The scratch graph/card are `/tmp/advanced-sampling-amplitude-research/double_box.dot` and `double_box.toml`. All seven internal masses are one, all numerators are one, and the complete native parent is `[5,8]`. Incoming p0,p1 and outgoing p2,p3 obey Q=p0+p1=p2+p3. Momentum conservation at the graph vertices gives

```
k=q5, l=q8,
q4=k-p1, q6=k-l, q7=k-Q, q9=l-p2, q10=l-Q.
```

Write E(v)=sqrt(1+|v|^2), using spatial vectors in E and physical external energies in the constant shifts. The four generated existing thresholds are

#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,auto,auto,right,),
    table.header([Name], [Edges], [Equation], [Loop-signature rank],),
    table.hline(),
    [A], [(5,7)], [E(k)+E(k-Qsp)-Q0], [1],
    [B], [(8,10)], [E(l)+E(l-Qsp)-Q0], [1],
    [C], [(5,6,10)], [E(k)+E(k-l)+E(l-Qsp)-Q0], [2],
    [D], [(6,7,8)], [E(k-l)+E(k-Qsp)+E(l)-Q0], [2],
  )]
  , kind: table
  )

For C the signature rows are (1,0),(1,-1),(0,1); D has the same independent rows in a different order. C and D therefore define bounded six-dimensional bodies. A and B need proper three-dimensional active blocks and complements; they must not be registered as full-rank six-dimensional charts.

The saved point is p0=(3,0,0,3), p1=(3,0,0,-3), p2=(3,3,0,0), with dependent p3=(3,-3,0,0). Thus Q=(6,0,0,0), s=36 and t=u=-18. At k=l=0, A=B=-4 and C=D=-3: zero is strictly inside both full-rank surfaces. Their actual six-dimensional radius is direction-dependent: C\'s positive root is sqrt(21)/2 along (k,l)=(r ex,0), and sqrt(21/2) along (k,l)=(r ex/sqrt(2),r ex/sqrt(2)). A constant sphere would be an incorrect replacement even though the zero-center condition is satisfied.

At these rest-frame equal-mass kinematics C and D coincide as scalar equations. Their edge sets remain different. They do not provide two independent normal coordinates or a regular C/D intersection. Their physical -Q0 cut orientations also require opposite signs of edge 6: C needs (+5,-6,-10), while D needs (+6,-7,+8). Consequently they cannot both be physical threshold factors of one fixed orientation. This is not a reason to drop either channel globally.

== Regular intersection in an actual generated orientation
<regular-intersection-in-an-actual-generated-orientation>
Take k=(sqrt(8),0,0), l=(0,sqrt(8),0). Then A=B=0 and C=D=sqrt(17)\>0. The two gradients in the complete (k,l) space occupy independent blocks; their Gram determinant is (32/9)^2=1024/81\>0. All internal energies are strictly positive, so this is a regular threshold intersection without a massless soft singularity or a coincident third threshold.

The original state used `store_atom=false`, so standalone export correctly refused to export expressions. A new scratch card `double_box_atoms.toml`, identical except for `store_atom=true`, generated an independent scratch state and retained all 98 orientations. The exported original CFF expression for orientation \#56, `[0,0,0,0,+,+,+,-,+,+,-]`, contains

```
(1/(B X) + 1/(X Y)) / (A D U product(E4,...,E10))
```

up to the common finite numerator/normalization. Here U=E4+E7-p0^0, X=E9+E10-(Q0-p2^0), and Y=E6+E7+E9-(Q0-p2^0). At the witness, U=X=sqrt(18), Y=sqrt(17)+sqrt(18), D=sqrt(17). The first term therefore contains both A and B in the same real orientation; the remaining factors are nonzero. This is stronger than merely intersecting two independently listed graph thresholds.

Exact scratch evidence is in `double_box_atoms_standalone/processes/amplitudes/` `massive_double_box/default/standalone_evaluators.json` and `double_box_orientation_56.txt`, under the same `/tmp` research directory. The expression comes from the `residue_map_id=56` branch of `graph_terms[0].original_integrand.single_parametric.exprs[0]`.

== Bounded promotion
<bounded-promotion>
This graph can now be promoted as the second UV-finite amplitude benchmark: first C/D full-rank zero-center channels plus ordinary LMB coverage, with the same production graph-derived map and Gaussian/moment gates as the kite. The A/B witness is reserved for the generic proper-subspace/complement and joint-normal implementation; it does not justify an amplitude-only shortcut. A boost can separate the coincident C/D equations, but no additional boosted benchmark is needed before the shared cross-section work progresses.
