#let integral-families = [
= Loop-integral families

`IntegralFamily` organizes inverse propagators as affine forms in loop scalar
products. It detects dependence, completes a family with irreducible scalar
products, partial-fractions dependent propagators, finds verified momentum shifts,
and rewrites scalar numerators
in propagator variables. The Rust
implementation belongs to `feynkit-graph`; Python exposes the same implementation
through `symbolica.community.feynkit`. Symbolica supplies the rank calculation
and exact linear solve.

Leave a family as the final expression in a notebook cell to see its ordered
inverse propagators, loop and external momenta, dimension, rank, completeness,
and independence. Expressions use Symbolica's native printer. `repr(family)`
gives a compact summary; IPython's plain-text display also lists the inverse
propagators in the same order.

== Conventions

Supply distinct unindexed loop and external momentum names. The external list
must be an independent basis; eliminate dependent momenta first. With `L` loop
momenta and `E` external momenta, the scalar-product space has
`L*(L+1)/2 + L*E` entries. On-shell conditions may constrain external products,
but must not constrain an integrated loop momentum.

The `denominators` argument contains inverse propagators, such as `(k+p)^2-m^2`,
not their reciprocals. Construct these expressions with
`Kinematics.scalar_product` so their Spenso notation and dimension agree.
Linear eikonal forms such as `k.p + delta` are also affine in this basis.
Prescriptions and integration contours are not inferred from these algebraic
expressions.

== Start from a generated diagram

`diagram.integral_family()` uses the diagram's stored loop routing, internal
propagators and model masses, then completes the basis with auxiliary scalar
products. Physical denominators come first, following `diagram.internal_edges`
in ascending edge-ID order. External carriers and dummy edges are excluded.
The shared annotated denominator builder preserves model mass conventions.

The family exposes its loop momenta, independent external momenta and scoped
`kinematics`. Use those returned names when setting external invariants, then
rebuild the family from the diagram. This recomputes propagators under the new
assumptions. The default dimension matches the diagram's symbolic dimension;
pass `kinematics=fk.Kinematics(D)` to choose a dimension explicitly.

Run this example from a repository checkout; substitute your own model path in
an installed workflow.

// docs-example: compile feynkit-diagram-integral-family
```python
from symbolica import E, S
import symbolica.community.feynkit as fk

model = fk.Model("crates/feynkit-py/tests/fixtures/scalars_2p_3p.json")
diagrams = model.generate_diagrams(
    ["scalar_0"], ["scalar_0"], loops=1, max_vertices=2,
    vertex_allow=["V_3_SCALAR_000"], allow_self_loops=False,
).diagrams
diagram = diagrams[0]
family = diagram.integral_family()
p = family.external_momenta[0]
s, x, y = S("s", "x", "y")
kin = family.kinematics.with_scalar_product(p, p, s)
family = diagram.integral_family(kinematics=kin)
U, F = family.symanzik([x, y])
assert U == x + y
assert (F + s*x*y).expand() == E("0")
```

Tree diagrams raise `DiagramError`, since `IntegralFamily` requires at least one
loop. Widths, imaginary prescriptions and custom UFO denominator formulas are
excluded, matching the existing diagram denominator API. For just the physical
propagators, use `diagram.propagator_family()`: it preserves repeated propagators
and bridges without completing or partial-fractioning them. To use another routing, apply
`diagram.with_loop_momentum_edges(...)` before extracting the family.

For a complete basis, use `diagram.integral_family()` or the equivalent
`fk.IntegralFamily.from_diagram(diagram)`. Both retain
the graph propagators first, then append independent scalar products until the
basis spans every loop-loop and loop-external product. Auxiliary denominators
have zero powers in the original scalar integral, or negative powers for
numerator factors.

Pass `independent_dot_products=[...]` to prefer particular
auxiliary denominators, expressed in the diagram's routed momentum names.
Candidates are tried in order; redundant entries are skipped and automatic
choices fill any remaining directions. The optional `kinematics` keyword has
the same meaning as for `diagram.propagator_family()`. For example, a two-loop
sunrise can use the two loop-external products:

```python
raw = sunrise.propagator_family()
p = raw.external_momenta[0]
products = [raw.kinematics.scalar_product(k, p) for k in raw.loop_momenta]
family = sunrise.integral_family(independent_dot_products=products)
assert family.is_complete and family.is_independent
```

Dependent graph propagators raise `DiagramError` in this constructor. Extract
them with `diagram.propagator_family()`, partial-fraction them, then call
`complete()` on each resulting independent family.

== A massive one-loop bubble

// docs-example: compile feynkit-integral-family-bubble
```python
from symbolica import E, S
import symbolica.community.feynkit as fk

D, k, p, s, m1sq, m2sq, d1, d2 = S(
    "family_docs::D", "family_docs::k", "family_docs::p", "s",
    "m1sq", "m2sq", "d1", "d2"
)
kin = fk.Kinematics(D, momenta=[k, p]).with_scalar_product(p, p, s)
family = fk.IntegralFamily(
    [k], [p],
    [kin.scalar_product(k, k) - m1sq,
     kin.scalar_product(k+p, k+p) - m2sq],
    kinematics=kin,
)
assert family.is_complete and family.is_independent
numerator = kin.scalar_product(k, p)**2
rewritten = family.rewrite_numerator(numerator, [d1, d2])
expected = (d2 - d1 - s + m2sq - m1sq)**2 / 4
assert (rewritten - expected).expand() == E("0")
```

`scalar_product_rules([d1, d2])` exposes the simultaneous substitution pairs
behind the rewrite. Labels must be distinct symbols or labeled calls and must
not collide with existing invariant names. Unknown numerator coefficients are
preserved; this operation does not integrate the numerator.

== Completion and dependent propagators

`family.complete()` appends scalar products until the inverse propagators form
a basis. Existing entries retain their positions. For the four denominators
`k^2`, `q^2`, `(k-p)^2` and `(q-p)^2`, the missing fifth direction is `k.q`.
The returned family includes this auxiliary inverse propagator; nonpositive
powers represent irreducible numerator factors in an IBP integral list.

To reuse inverse propagators from another family, pass them in priority order:
`family.complete(candidates=other_family.denominators)`. Express them in the
current family's momentum coordinates first. The same scoped kinematics are
applied before testing independence. For the example above,
`family.complete(candidates=[kin.scalar_product(k-q, k-q)])` appends `(k-q)^2`
instead of `k.q`. Existing denominators keep their positions; dependent
candidates are skipped, and scalar products fill any directions the pool does
not span. Every supplied candidate is validated as affine, even if an earlier
entry already completed the basis. Repeating completion preserves the result.

Completeness and independence are different. The pair `k^2` and `k^2-m^2`
spans the one-dimensional scalar-product space, but its two entries are
dependent. Completion and numerator solving reject this family until it has
been partial-fractioned.

== Partial fractions in family coordinates

`family.partial_fraction(powers)` returns `(coefficient, powers)` pairs in the
original family order. Positive powers denote denominators and negative powers
denote numerator factors. Every output term has an independent set of
positive-power propagators. Both affine relations with mass shifts and homogeneous
relations are handled, including repeated poles.

// docs-example: compile feynkit-integral-family-partial-fractions
```python
from symbolica import E, S
import symbolica.community.feynkit as fk

k, m2 = S("apart_docs::k", "m2")
kin = fk.Kinematics()
kk = kin.scalar_product(k, k)
family = fk.IntegralFamily([k], [], [kk, kk-m2], kinematics=kin)
terms = family.partial_fraction([1, 1])
reconstructed = E("0")
for coefficient, powers in terms:
    term = coefficient
    for denominator, power in zip(family.denominators, powers):
        term *= denominator**(-power)
    reconstructed += term
assert (reconstructed - 1/(kk*(kk-m2))).together() == E("0")
```

The algorithm retains propagator labels and gives an exact rational-function
identity. It does not shift loop momenta or discard scaleless terms. Apply
exceptional kinematics before constructing the family: generic external
invariants can occur in coefficient denominators. The optional `max_states`
budget limits intermediate exponent vectors; exhaustion raises an error without
returning a truncated decomposition.

== Verified momentum shifts and subtopologies

`source.find_mapping(target)` derives affine loop-shift candidates from independent
quadratic propagators and verifies every denominator. The search includes loop
permutations, reversals, mixtures, and external shifts. A target may have extra
propagators, allowing the source to embed as a subtopology. Every accepted map has
real coefficients and loop-variable determinant `+1` or `-1`, so the loop measure
is unchanged. Dimensions, external names, and external scalar-product assumptions
must agree.

Fractional shifts and reciprocal loop rescalings are supported. The shared
coefficient extraction preserves rational row factors before constructing
Symbolica matrices. Clearing these factors is valid when solving equations,
but would change a Jacobian, a Symanzik determinant, or a relation between
the original inverse propagators. Regression cases include the light-cone
shift `k-q+n-nbar/2`, the unit-Jacobian rescaling `(k,q) -> (2*l,r/2)`,
and exact partial fractions with fractional quadratic coefficients.

// docs-example: compile feynkit-integral-family-mapping
```python
from symbolica import E, S
import symbolica.community.feynkit as fk

k, l, p, s, m1sq, m2sq = S("map_docs::k", "map_docs::l", "map_docs::p", "s", "m1sq", "m2sq")
kin = fk.Kinematics(momenta=[k, l, p]).with_scalar_product(p, p, s)
source = fk.IntegralFamily(
    [k], [p], [kin.scalar_product(k, k)-m1sq, kin.scalar_product(k+p, k+p)-m2sq],
    kinematics=kin,
)
target = fk.IntegralFamily(
    [l], [p], [kin.scalar_product(l, l)-m2sq, kin.scalar_product(l-p, l-p)-m1sq],
    kinematics=kin,
)
mapping = source.find_mapping(target)
assert mapping is not None
assert mapping.denominator_map == [1, 0]
assert mapping.map_powers([2, -1]) == [-1, 2]
for i, j in enumerate(mapping.denominator_map):
    assert (mapping.apply(source.denominators[i])-target.denominators[j]).together() == E("0")
```

`mapping.momentum_rules` exposes the loop images. `mapping.apply(numerator)` uses
the verified simultaneous scalar-product substitutions; contract tensor indices
to Spenso dots first. `map_powers` reorders signed powers and inserts zeros for
unused target propagators.

Supply a known transformation with `source.mapping_to(target, loop_images)`.
This also supports purely eikonal families, which cannot supply a quadratic basis
for automatic search. Unknown or complex coefficients do not pass the real-map
check; known-real Symbolica assumptions can establish realness. A failed
equivalence check returns `None`; incompatible inputs and exhausted
`max_candidates` budgets raise an error.

The Rust regression suite checks the published three-loop `topo4 -> topo1`
example from #link("https://feyncalc.github.io/FeynCalcBookDev/FCLoopFindTopologyMappings.html")[FeynCalc's mapping documentation],
including automatic discovery. The search does not implement Pak's parametric
equivalence algorithm, external-momentum exchanges, or identities without an
affine loop-momentum map. The separate polynomial comparison below handles
parameter permutations, and the scaling certificate below detects scaleless
sectors. Kira/FIRE execution and master-integral evaluation remain separate work.

== Symanzik polynomials

`family.symanzik(parameters)` returns `(U, F)` with one Feynman parameter per
inverse propagator, in the original order. With
`sum_i x_i D_i = k.M.k + 2 k.Q + J`, the convention is
`U = det(M)` and `F = Q.adj(M).Q - U*J`. It matches Minkowski denominators
`k^2-m^2`; for a single massive tadpole, `U=x` and `F=m^2*x^2`.
Symbolica computes the determinant, cofactors and rational cancellation.
For nonsingular matrices this equals `U*(Q.M^-1.Q - J)`.

// docs-example: compile feynkit-integral-family-symanzik
```python
from symbolica import E, S
import symbolica.community.feynkit as fk

k, p, s, m1sq, m2sq, x1, x2 = S("uf_docs::k", "uf_docs::p", "s", "m1sq", "m2sq", "x1", "x2")
kin = fk.Kinematics(momenta=[k, p]).with_scalar_product(p, p, s)
family = fk.IntegralFamily(
    [k], [p], [kin.scalar_product(k, k)-m1sq, kin.scalar_product(k-p, k-p)-m2sq],
    kinematics=kin,
)
U, F = family.symanzik([x1, x2])
assert U == x1 + x2
assert (F - ((x1+x2)*(m1sq*x1+m2sq*x2)-s*x1*x2)).expand() == E("0")
```

Parameters must be distinct new symbols or labeled calls. The quadratic loop
matrix may be singular: the adjugate formula still defines algebraic polynomials
without an inverse. For example, the two-loop denominators
`(k+l)^2-m^2` and `(k-l).p+delta`, with `p^2=s`, give `U=0`, `F=s*x*y^2`.
This does not define a Gaussian integral or prove the sector scaleless.
`parametric_mapping` also rejects singular forms because their polynomials can
lose physical information, such as the masses in this example.
Eikonal denominators may accompany quadratic propagators. Choose a
family containing the propagators of the desired sector before preparation;
auxiliary numerator denominators are not automatically omitted.

The returned polynomials do not include a parameter measure, Gamma prefactor,
propagator powers, numerator derivatives or an integration prescription.
They are a building block for Feynman parametrization and topology comparison.
Regression tests compare the bubble and massless two-loop self-energy with
#link("https://feyncalc.github.io/FeynCalcBookDev/FCFeynmanPrepare.html")[FeynCalc's preparation examples].

== Compare families in parameter space

`source.parametric_mapping(target, parameters)` finds a parameter permutation
that makes both `U` and `F` identical. It returns a `PropagatorMapping`, or `None`
when the polynomials differ. Source and target must have equal propagator counts
and compatible external kinematics. The supplied parameters must be absent from
both families' expressions.

Each polynomial becomes an incidence graph. Parameter vertices carry no names;
term vertices retain their exact coefficients and whether they belong to `U` or
`F`; edge labels retain exponents. Symbolica supplies the graph canonicalization.
This addresses the parameter-permutation comparison performed by
#link("https://feyncalc.github.io/FeynCalcBookDev/FCLoopToPakForm.html")[FeynCalc's Pak-form backend]
without duplicating a canonicalization algorithm in FeynKit.

The following published example has equal parameter polynomials even though the
fixed-external-momentum search finds no shift. It appears as `topos4` in
#link("https://feyncalc.github.io/FeynCalcBookDev/FCLoopFindTopologyMappings.html")[FeynCalc's mapping examples].

// docs-example: compile feynkit-integral-family-parametric-mapping
```python
from symbolica import E, S
import symbolica.community.feynkit as fk

k, l, p, q, s, m2 = S("uf_map_docs::k", "uf_map_docs::l", "uf_map_docs::p", "uf_map_docs::q", "s", "m2")
xs = list(S("x1", "x2", "x3", "x4", "x5"))
kin = (fk.Kinematics(momenta=[k, l, p, q])
       .with_scalar_product(p, p, E("0"))
       .with_scalar_product(q, q, E("0"))
       .with_scalar_product(p, q, s/2))
def square(v):
    return kin.scalar_product(v, v)
source = fk.IntegralFamily(
    [k, l], [p, q],
    [square(k+p)-m2, square(k-l), square(l+p)-m2, square(l-q)-m2, square(l)],
    kinematics=kin,
)
target = fk.IntegralFamily(
    [k, l], [p, q],
    [square(k-l)-m2, square(k-q), square(l-q)-m2, square(l+p)-m2, square(l)],
    kinematics=kin,
)
assert source.find_mapping(target) is None
mapping = source.parametric_mapping(target, xs)
assert mapping is not None
mapped_parameters = [xs[j] for j in mapping.denominator_map]
source_U, source_F = source.symanzik(mapped_parameters)
target_U, target_F = target.symanzik(xs)
assert (source_U-target_U).expand() == E("0")
assert (source_F-target_F).expand() == E("0")
target_powers = mapping.map_powers([1, 2, 1, 1, -1])
```

Polynomial equality is a formal parameter-space result. Prescriptions and
contours are additional input, and no external-momentum or tensor-numerator
substitution is implied. The result intentionally exposes only the propagator
permutation and power reordering. Express scalar numerator factors in family
coordinates first; use a verified momentum mapping to transform a general
contracted numerator. Sector discovery and subtopology minimization remain
separate work.

== Family collections

`IntegralFamily.find_mappings(families)` groups a collection and returns one
`(target_index, mapping)` per input family, in input order. Target indices refer
to the original collection. A retained representative maps to itself; every
other map goes directly to a representative, so there are no mapping chains.
Earlier entries are preferred. Apply `mapping.map_powers` to retain signed
propagator powers, and `mapping.apply` for scalar numerator substitutions.

The shared implementation canonizes each Symanzik pair once with Symbolica and
compares affine momentum shifts only inside matching groups. Every merge still
requires a real shift with unit absolute Jacobian and exact denominator
identities. Equal Symanzik polynomials without a verified momentum map remain
separate. All inputs need compatible external kinematics and nonsingular
quadratic forms; unsupported shift searches and per-pair candidate-budget
exhaustion raise errors rather than certifying inequivalence.

Select positive-power sectors and remove certified scaleless contributions
before grouping. Complete the retained families afterwards, preserving their
original denominator order and adding zero auxiliary powers. Different
propagator counts remain separate in this operation. Subtopology embeddings
remain available through `source.find_mapping(target)`; external-momentum
exchanges and integration prescriptions are not inferred.

// docs-example: compile feynkit-integral-family-grouping
```python
from symbolica import E, S
import symbolica.community.feynkit as fk

k, l, p, s = S("group_docs::k", "group_docs::l", "group_docs::p", "s")
kin = fk.Kinematics(momenta=[k, l, p]).with_scalar_product(p, p, s)
square = lambda v: kin.scalar_product(v, v)
families = [
    fk.IntegralFamily([k], [p], [square(k)-1, square(k+p/2)-2], kinematics=kin),
    fk.IntegralFamily([l], [p], [square(l)-2, square(l-p/2)-1], kinematics=kin),
]
mappings = fk.IntegralFamily.find_mappings(families)
assert [target for target, mapping in mappings] == [0, 0]
assert mappings[1][1].map_powers([2, -1]) == [-1, 2]
```

== Sector support and scalelessness

`family.sector(powers)` retains the positive-power propagators in their original
order. It preserves loop variables and external kinematics. Zero-power entries
are absent, and negative-power entries are numerator factors rather than sector
denominators. This selects support for analysis; it does not discard a numerator
from the original integral.

`sector.scaleless_scaling(parameters)` returns exact weights satisfying
`sum_i w_i*x_i*dG/dx_i = G`, where `G=U+F`. Symbolica solves the linear system
formed from monomial exponents. A solution is a scalelessness certificate in
dimensional regularization at generic dimension. This uses
#link("https://arxiv.org/html/1310.1145")[Lee's parametric zero-sector criterion,
equations (13)--(16)]. `None` means the criterion found no certificate; it does
not assert that the integral is nonzero. A singular quadratic loop matrix raises
an error rather than treating a degenerate polynomial identity as a scaling
certificate. This restriction is stronger than the algebraic `symanzik` API.

// docs-example: compile feynkit-integral-family-scaleless-sector
```python
from symbolica import E, S
import symbolica.community.feynkit as fk

k, l, m2, M2, x, y, z = S("scale_docs::k", "scale_docs::l", "m2", "M2", "x", "y", "z")
kin = fk.Kinematics(momenta=[k, l])
kk = kin.scalar_product(k, k)
family = fk.IntegralFamily(
    [k, l], [], [kk, kin.scalar_product(l, l)-m2, kk-M2], kinematics=kin,
)
assert family.scaleless_scaling([x, y, z]) is None
sector = family.sector([1, 1, -2])
assert len(sector.denominators) == 2
weights = sector.scaleless_scaling([x, y])
assert weights is not None
U, F = sector.symanzik([x, y])
G = U + F
euler = weights[0]*x*G.derivative(x) + weights[1]*y*G.derivative(y)
assert (euler-G).expand() == E("0")
```

Here the selected sector contains a massless integration over `k`, despite the
massive propagator in the independent `l` integration. The certificate also
covers positive powers and polynomial numerators in this sector. Apply any
exceptional kinematics before constructing the family so vanished scales are
visible to the test.

`sector.scaleless_transverse_direction()` covers a different source of zero
sectors, including degenerate quadratic forms with eikonal denominators. It
returns a nonzero real vector `w` in loop order such that every denominator is
invariant under `k_i -> k_i + w_i*r_perp`, with `r_perp` orthogonal to the
external span. The implementation stacks the individual quadratic loop
matrices and uses Symbolica to find a common null direction. It also requires
a nonsingular external Gram matrix and a nonempty transverse space. Symbolic
dimension is generic; concrete dimension must exceed the external basis size.
`None` is inconclusive, and a vanishing Symanzik `U` alone is insufficient.

A real change of loop coordinates isolates an unrestricted transverse integral
whose integrand is polynomial, even when the remaining denominators carry
masses or eikonal scales. Its vanishing follows from the polynomial-integral
corollary in the
#link("https://arxiv.org/html/2203.13014v3")[SAGEX review, section 1.1]. The
common-null-direction test is our application of that result. Polynomial
numerators do not obstruct the proof; exclude their entries with `sector`
before testing. Complex directions are rejected because they would require
an additional contour argument.

// docs-example: compile feynkit-integral-family-transverse-sector
```python
from symbolica import E, S
import symbolica.community.feynkit as fk

k, q, p, s, m2, delta = S("transverse_docs::k", "transverse_docs::q", "transverse_docs::p", "s", "m2", "delta")
kin = fk.Kinematics(momenta=[k, q, p]).with_scalar_product(p, p, s)
family = fk.IntegralFamily(
    [k, q], [p],
    [kin.scalar_product(k+q, k+q)-m2, kin.scalar_product(k-q, p)+delta],
    kinematics=kin,
)
assert family.scaleless_transverse_direction() == [E("1"), E("-1")]
```

Regressions additionally check an on-shell massless bubble and the massless and
massive eikonal examples from
#link("https://feyncalc.github.io/FeynCalcBookDev/FCLoopPakScalelessQ.html")[FeynCalc's scalelessness documentation].
The API returns a certificate; it does not automatically erase terms, separate
UV and IR poles, or enable the generator's `only_scaleless` filter.

The conventions match the algebraic tasks described by FeynCalc's
#link("https://feyncalc.github.io/FeynCalcBookDev/FCLoopBasisFindCompletion.html")[basis completion]
and #link("https://feyncalc.github.io/FeynCalcBookDev/FCLoopBasisOverdeterminedQ.html")[dependence checks].
The propagator decomposition addresses the algebraic part of
#link("https://feyncalc.github.io/FeynCalcBookDev/ApartFF.html")[ApartFF]; its optional
momentum shifts and scaleless-integral removal are separate tasks.
The #link("guides/feyncalc-coverage/")[gallery audit] tracks the remaining
end-to-end loop calculations separately.
]
