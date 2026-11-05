# Massless propagator integrals (p-integrals) in vakint: binding interface

Status: binding M1 interface (2026-11-03). Implementation agents follow it; a deviation needs a recorded decision.
Code vocabulary: `src/p_integral.rs` (registry data and shared types, compiled and tested).
Scope: 1L-4L massless single-scale propagator integrals with Forcer coverage; 5L is next and must not be blocked.
Out of scope: R*, massive p-integrals, tensor p-integrals (open loop indices or loop momenta contracted with
spectator vectors), non-integer user powers.

Ownership in M1:
- frontend agent: sections 1-3 (input, inference, normalization, registry, catalogue, matching, routing, DOT);
- Forcer agent: sections 4-5 (vendoring, FORM 5 detection, `EvaluationMethod::Forcer`, the shared p-finalizer,
  the p-master table, `LoopNormalizationFactor::GScheme`);
- RustRed agents: section 6 (M2).
- Both vakint agents share `src/p_integral.rs`: extend it, never renumber `P_FAMILIES`.

## 1. Atom syntax

A p-integral is an ordinary vakint integral:

```text
numerator * topo(prop(id, edge(a, b), momentum, 0, power) * ...)
```

- Momentum: a linear combination of loop momenta `k(j)` and exactly one external head `p(i)`, for any
  `i >= 1`.
  - `edge(a, b)` means the momentum flows from `a` to `b`, as for vacuum integrals.
  - At every vertex `v`, `balance(v) = sum(out) - sum(in)`. The integral is valid only if
    `balance = +p(i)` at one vertex (where q enters), `-p(i)` at one other vertex (where q leaves), and `0` everywhere
    else.
  - Prop ids are distinct integers (positive by convention; the DOT ingress uses FeynKit edge ids) and need not be
    contiguous; vertex ids are nonnegative integers.
- Mass: the literal `0` in every prop.
- Power: an integer. `power <= 0` pinches the line, and a negative power becomes a numerator factor.
- Numerator: a polynomial in `dot(k(j),k(l))`, `dot(k(j),p(i))`, `dot(p(i),p(i))` (or the contracted index form
  `k(j,mu)*p(i,mu)`), with any ε- or parameter-dependent scalar coefficients.
  - Other external vectors `p(m)`, `m != i`, may appear only in loop-momentum-free factors.
- Metric: vakint's Minkowski metric (+,-,-,-). Propagators are `1/(k_e^2)^a`. The Euclidean scale is
  `Q^2 = -dot(p(i),p(i))`, which is positive for spacelike `p(i)`.
- Short form: `topo(NAME(±p(i), a_1, ..., a_N))`.
  - `NAME` is a family of section 2, the `a_n` are the integer powers of its N lines in canonical order (`0` pinches
    a line), and the canonical `p(1)` is the first argument.
  - ISP powers never appear in the short form: ISPs go in the numerator, as dot products.

Example, the 2-loop kite with one dot and a rank-2 numerator in the user's own routing, with q written as `p(3)`:

```text
dot(k(1),p(3))^2 * topo(prop(1,edge(1,3),k(1),0,2)*prop(2,edge(1,4),p(3)-k(1),0,1)
  *prop(3,edge(3,4),k(2),0,1)*prop(4,edge(3,2),k(1)-k(2),0,1)*prop(5,edge(4,2),p(3)-k(1)+k(2),0,1))
```

Kind inference (no flag):
- The integral is a p-integral (`IntegralKind::Propagator`) iff some prop momentum contains a head `p(_)`, or its
  short-form head is a p family name. Otherwise it is a vacuum integral (`IntegralKind::Vacuum`), exactly as today.
- This is unambiguous: no vacuum topology has `p` in its prop momenta, and the vacuum `UNKNOWN` matcher refuses them.
  The Atom therefore stays self-describing without a flag argument.
- Out of domain (M1 review F1): a `topo` with `p` in its momenta that is well formed but not a p-integral vakint covers is
  an unknown integral under `allow_unknown_integrals` (the default), exactly as before p-integrals existed, and an error
  (items 1, 2, 4 and 8 below) otherwise. Out of domain means: a massive line (item 1), several external heads (item 2),
  a non-integer power (item 4), more than four loops or a graph outside the registry (item 8). The unknown integral is
  `topo(UNKNOWN(...))` of the long form (after the exact normalization of items 5-6 when only item 8 applies); its
  momenta keep `p(i)` and its loop count is the number of `k(j)`. No backend evaluates it (pySecDec takes external
  momenta only in numerators), so evaluating it is an error; `to_canonical` passes it through. A short form out of domain
  (a non-integer power) is no unknown integral and keeps its error. Malformed inputs (not a `prop`, a bad edge or
  vertex id, repeated prop ids, no external head, broken momentum conservation, a loop-count mismatch, a momentum that is
  not an integer combination of heads, a non-unimodular routing, a tensor numerator) are errors under any setting.
- In code: `p_integral::input::PRejection { OutOfDomain(String), Malformed(VakintError) }`; the fallback is
  `Topologies::match_propagator(input, allow_unknown_integrals)`.
- Inference lives in `p_integral::infer_kind(integral: AtomView) -> Result<IntegralKind, VakintError>`.
  `Vakint::match_topologies_for_settings` dispatches on it before any vacuum or five-loop matcher runs.

Validation, in this order. Every error is `VakintError::InvalidIntegralFormat("p-integral: ...")` unless stated:
1. A nonzero mass in any prop: "p-integrals are massless; prop {id} has mass {m}".
2. More than one external head, e.g. `p(1)` and `p(2)`: "a p-integral has one external momentum; found {heads}".
   A sum `p(1)+p(2)` injected at one vertex is rejected the same way: define `q` as a single head.
3. Vertex balances not of the form above: "momentum conservation: {vertex} -> {balance}".
4. Non-integer power: "powers must be integers". A cycle rank different from the number of distinct `k(j)` heads:
   "the routing does not match the graph's loop count".
5. Normalization (exact, before matching):
   - a. A line whose momentum has no loop momentum (it carries only `±p(i)`) is a tree factor. It becomes the
     numerator factor `dot(p(i),p(i))^(-a)`, and the line is contracted.
   - b. Lines with identical quadrics (`q_e^2` equal after expansion) merge into the lowest id, adding their powers.
   - c. Lines with `power <= 0` are contracted (their vertices merge). A negative power multiplies the numerator by
     `dot(q_e,q_e)^(-a)` in the user's momenta.
6. Scaleless check, which gives an exact `0` term and is not an error. After normalization the integral vanishes iff:
   - the two external vertices coincide; or
   - a self-loop remains; or
   - a vertex `v` separates the graph into a component that contains neither external vertex.

   The same rule also makes a vacuum `topo` with every mass literally `0` and no `p` exactly zero.
7. Numerator: an open loop index, or a loop momentum contracted with `p(m)`, `m != i`, gives
   `VakintError::InvalidNumerator("p-integral: tensor p-integrals are not supported yet: ...")`.
   The tensor prepass (`tensor_reduce_for_evaluation`) must not project p-integral terms; it only puts them in
   scalar-product form.
8. More than 4 loops, or a graph that is in no catalogue (section 2): `VakintError::UnreckognizedIntegral`, with the
   note "p-integral not covered by the 1L-4L registry".

## 2. Topology registry and matching

### Families

Families: cut one edge orbit of every (L+1)-loop cubic vacuum graph that has no 1- or 2-edge cut. This gives
1 / 1 / 3 / 11 families at 1 / 2 / 3 / 4 loops, frozen in `p_integral::P_FAMILIES`. The generator lives in
`tools/p_integral_registry/`; its 4L tops match Forcer's named top-level topologies one to one.

| name | L | Forcer | vacuum parent, cut edge orbit | lines | ISPs |
|---|---|---|---|---|---|
| `P1L` | 1 | - | 2L sunset | 2 | - |
| `P2L` | 2 | - | 3L tetrahedron K4 (orbit 6) | 5 | - |
| `P3L_LA` | 3 | (Mincer LA) | 4L prism, a rung (3) | 8 | (k(1)-k(2))^2 |
| `P3L_BE` | 3 | (Mincer BE) | 4L prism, a triangle edge (6) | 8 | (k(2)-p(1))^2 |
| `P3L_NO` | 3 | (Mincer NO) | 4L K3,3 (9) | 8 | (k(1)-k(2))^2 |
| `P4L_LA4` | 4 | la4 | 5L graph A, edge 25 (1) | 11 | (k(1)-k(2))^2, (k(1)-k(3))^2, (k(2)-k(3))^2 |
| `P4L_J1` | 4 | j1 | graph A, edge 01 (2) | 11 | (k(1)-k(3))^2, (k(2)-p(1))^2, (k(3)-p(1))^2 |
| `P4L_J2` | 4 | j2 | graph A, edge 02 (4) | 11 | (k(1)-k(2))^2, (k(1)-p(1))^2, (k(2)-p(1))^2 |
| `P4L_J4` | 4 | j4 | graph B, edge 36 (6) | 11 | (k(1)-k(4))^2, (k(2)-k(4))^2, (k(1)-p(1))^2 |
| `P4L_J5` | 4 | j5 | graph A, edge 03 (4) | 11 | (k(1)-k(2))^2, (k(1)-k(3))^2, (k(3)-p(1))^2 |
| `P4L_J6` | 4 | j6 | graph A, edge 34 (1) | 11 | (k(1)-k(3))^2, (k(1)-k(4))^2, (k(4)-p(1))^2 |
| `P4L_J11` | 4 | j11 | graph B, edge 03 (3) | 11 | (k(1)-k(2))^2, (k(1)-k(3))^2, (k(2)-k(3))^2 |
| `P4L_HAHA` | 4 | haha | cube Q3 (12) | 11 | (k(1)-k(2))^2, (k(1)-k(4))^2, (k(2)-p(1))^2 |
| `P4L_NO1` | 4 | no1 | Wagner M8, a rim edge (8) | 11 | (k(1)-k(2))^2, (k(1)-k(3))^2, (k(2)-p(1))^2 |
| `P4L_NO2` | 4 | no2 | Wagner M8, a spoke (4) | 11 | (k(1)-k(2))^2, (k(2)-k(3))^2, (k(3)-p(1))^2 |
| `P4L_NO3` | 4 | no3 | graph B, edge 01 (3) | 11 | (k(2)-k(3))^2, (k(1)-p(1))^2, (k(3)-p(1))^2 |

Graph A has edges {01,02,03,12,14,25,34,36,47,56,57,67}. Graph B has edges {01,02,03,12,14,25,36,37,46,47,56,57}.

Canonical routing (the `lines` long form of each entry):
- q (`p(1)`) enters at vertex 1 and leaves at vertex 2. The lines carrying `p(1)` form one path between them.
- Line `i <= L` carries `k(i)` alone; the other lines are BFS-tree lines. Coefficients are in {-1,0,1}.
- ISPs are greedy, in the order (differences of two loop momenta, then k(i)-p(1)). The family denominators are the
  lines followed by the ISPs. This order is the RustRed slot order.
- Three `src/p_integral.rs` tests pin the counts, the q flow and the completeness of the scalar products.

### Registry code

Extend the vacuum code; do not duplicate it:
- Add `Topology::Propagator(Integral)`, a distinct variant, so that no vacuum backend can claim a p-family by
  loop count.
- Add the constructor `Integral::new_propagator(family: &PFamilyData)`.
  - It admits `p(1)` in momenta and mass `0`, and its short pattern is `NAME(p1_, pow1_, ..., powN_)`.
  - `applicable_evaluation_methods` is `[RustRed, Forcer]`.
- `Topologies::generate_topologies()` appends the 16 families with no contraction entries: the catalogue covers
  pinches.
- Admit external heads in `Graph::new_from_atom`, `canonicalize_lmb` (the loop-id check counts `k` heads only) and
  `get_individual_momenta`.

### Catalogue

Add `PCatalogue`, modelled on `FiveLoopCatalogue` / `momentum_mapping`, built lazily per loop count and shared
behind `OnceLock`:
- For each family and each line subset S, contract S. Drop scaleless results (rule of section 1, item 6).
- Key: the Symbolica `CanonicalForm` of the contracted multigraph, with the two external vertices coloured alike
  (q reversal is a symmetry).
- Value: `(family index, sorted S, vertex map)`. The first family in table order wins, then the lexicographically
  smallest S. This is deterministic.
- 4L size: 11 x 2^11 = 22,528 contractions.
- At 5L (64 x 2^14 ≈ 1 M) eager enumeration is not affordable. The catalogue API must allow a shipped index or an
  on-demand search, and nothing may assume eager enumeration.
- Implemented (frontend, M1): the catalogue is split into strata of one loop count and one line count, each built
  under its own `OnceLock` on first lookup. A lookup builds only its stratum (at 5L, inputs with few pinches touch
  small strata), and a shipped index can replace one stratum without changing callers.

### Matching

Matching of a normalized input:
1. Look up the canonical form of the input graph.
   - Implemented: the input is first labelled canonically with its powers as edge colours, and that labelled graph
     is then keyed without powers. The line correspondence thus depends only on the (graph, powers) isomorphism
     class: reroutings, edge flips, q reversal and head renaming give identical canonical output, and vertex
     relabelings do too unless an automorphism preserves the powers (then the outputs are equivalent).
2. Read off the line correspondence (try each assignment of parallel lines) and the orientation `σ`.
   `σ = +1` iff the user's in-vertex maps to canonical vertex 1. Parallel lines are assigned in order of power.
3. Solve the affine routing (below).
4. Build `ReplacementRules`:
   - `canonical_topology = Topology::Propagator(family)`;
   - `canonical_expression_substitutions`: `pow(n)` -> the matched power (`0` for pinched lines),
     `p(1) -> σ·p(i)`;
   - `numerator_substitutions`: `k(j,idx___)` -> its image;
   - `edge_ids_canonical_to_input_map`;
   - a new field `propagator: Option<PropagatorWitness>`. The witness also carries `powers` (the matched power of
     every canonical line, 0 for pinched ones), so the glue reads the backend `powers` and `external`
     (`canonical_external() = σ·p(i)`) from it directly.
5. No graph match: `UnreckognizedIntegral`. There is no quadric-search fallback in M1. Normalization (merging
   identical quadrics, contracting pinched lines) turns self-energy-insertion graphs into contractions of the tops.
   The census test of section 7 proves coverage, and any uncovered graph is added as a new data family, never as
   special code.

### Affine routing

Extend `routing.rs`; the vacuum call is the `None` case:

```rust
pub(crate) struct AffineRouting {
    /// (input k(j), image Σ M_jl k(l) + c_j·p_user), M integer with |det M| = 1.
    pub images: Vec<(Atom, Atom)>,
    /// The canonical p(1) is external_sign·p_user.
    pub external_sign: i8,
}
pub(crate) fn input_loop_momenta_in_canonical_basis_affine(
    equations: &[(Atom, Atom)],          // (input line momentum, ± canonical line momentum)
    external: Option<(&Atom, &Atom, i8)>, // (user head p(i), canonical head p(1), σ from the vertex map)
) -> Result<AffineRouting, VakintError>;
```

Algorithm:
1. Substitute `p(1) -> σ·p_user` in the canonical momenta, with σ taken from the vertex map.
2. Solve for the input loop momenta on an independent set of lines, treating `p_user` as known.
3. Check every line up to its sign and `|det M| = 1`.
4. Search the line signs as today. If no sign choice works, retry with `-σ` before failing.
   Implemented order: written orientations with `σ`, then with `-σ`, then the line-sign search with `σ` and `-σ`.
   An input already written in the family routing (every line equal to its oriented image) keeps its loop momenta.

The canonical numerator is the user numerator with these images; it keeps the user's head `p(i)`. It is returned
as an expanded polynomial in `dot(k(j),k(l))`, `dot(k(j),p(i))` and `dot(p(i),p(i))`.

`Vakint::to_canonical` returns the family long form with matched powers (zero-power lines kept) and `σ·p(i)`, or the
short form `NAME(σ·p(i), a_1..a_N)`.

## 3. DOT ingress (HEPKit / FeynKit conventions)

Input: a compact FeynKit DOT with exactly two dangling legs, written as invisible endpoint nodes:

```text
digraph p2l_kite {
  exte0 [style=invis]; exte1 [style=invis];
  exte0 -> v1 [id=0 particle="phi"];
  v2 -> exte1 [id=1 particle="phi"];
  v2 -> v4 [id=2 particle="phi" lmb_id=0];
  v3 -> v4 [id=3 particle="phi" lmb_id=1];
  v4 -> v1 [id=4 particle="phi"];
  v3 -> v2 [id=5 particle="phi"];
  v1 -> v3 [id=6 particle="phi" power=2];
}
```

`VakintExpression::from_diagram` extension, behind the `feynkit-ingress` feature:
- It accepts uncut diagrams with exactly two external legs. Every internal denominator must be massless after
  `parameter_substitutions`.
- `p(1)` is the momentum flowing into the graph through the external leg with the smallest edge id. It is derived
  from that leg's stored signature, whichever external FeynKit chose as dependent.
- Leg particle and orientation (both legs incoming is allowed) only fix that sign.
- The new option `DiagramIntegralOptions::external_momentum: Option<Atom>` renames `p(1)` (it must be a bare
  `p(i)`).
- `external_momenta` (spectators) must be empty for p-graphs in M1.
- Powers, physical slots, contraction and numerator handling are unchanged. FeynKit `P(0)` scalar products in the
  numerator map to `p(1)` products. Implemented: every graph line is written as a `prop` with its power, zero and
  negative ones included, and the ordinary normalization of section 1 (tree lines, merged quadrics, contractions,
  scaleless zero) applies to the resulting Atom exactly as to a written one.

Export: `p_integral::to_dot(integral: AtomView, particle: &str) -> Result<String, VakintError>`.
- Input: a p-integral `topo` with nonnegative powers and no tree line: zero-power lines are contracted and repeated
  quadrics merged first (section 1, item 5); anything that needs a numerator factor, or is scaleless, is rejected.
- It writes the DOT above: `exte0 -> v{in}`, `v{out} -> exte1`, one edge per line with `power=a` when `a != 1`,
  and `lmb_id` on the lines that carry a single `k(j)`. Numerators are not exported.
- Test model: one scalar `phi` with vertices of valence 3 to 10 (fixture `tests/fixtures/phi_p_graphs.json`).

Round trip:
- `atom -> to_dot -> FeynmanDiagram::from_dot -> from_diagram` gives an atom with the same family, sector, witness
  sign convention and identical Retain-mode results.
- `dot -> atom -> dot -> atom` is a fixed point after the first export.

## 4. Normalization and output

Definition (Minkowski, as for vacuum):

```text
I = [F]^L ∫ ∏_j d^d k_j  N(k, q) / ∏_e (k_e^2)^(a_e),    d = 4 - 2ε
```

- `F` is `VakintSettings::integral_normalization_factor`. One setting serves vacuum and p-integrals: an expression
  may mix both, and the measure must be common.
- Wick rotation gives `I = [F·i·π^(d/2)]^L (-1)^(Σa + r) I_E`.
  - `I_E` uses the measure `d^d k_E/π^(d/2)` and Euclidean props.
  - `r` is the degree of a numerator monomial in Minkowski scalar products, `q·q` included.
- `I_E = (Q^2)^(w - Lε) (εG(1,1))^L C_G(ε)`, where:
  - `Q^2 = -dot(p(i),p(i))`;
  - `w = 2L - Σa + r` (an integer);
  - `G(1,1) = Γ(ε)Γ(1-ε)^2/Γ(2-2ε)`;
  - `C_G` is the G-scheme value at `Q^2 = 1`. This is exactly Forcer's normalization, checked by probe: the bubble
    is `GschemeConstants(0,0) = 1/ε`.

Backend boundary (`p_integral::ReducedPropagatorValue`): every backend returns, per homogeneous part,
`scale_power = w` and `reduced = (-1)^(Σa+r) C_G`:
- `reduced` is exact in ε, in the common p-master basis:
  - `vakint::PMaster(name)`: Forcer's `Master(name)`, e.g. `no1`, `haha`, `t105`, `fastar2`;
  - `vakint::PG(x1,x2)`, with `x1 >= x2 >= 0`: Forcer's `GschemeConstants(x1,x2)`, which is exactly
    `(1+x1+x2)·G(1+x1ε, 1+x2ε)/(ε G(1,1))`, with
    `G(a,b) = Γ(a+b-d/2)Γ(d/2-a)Γ(d/2-b)/(Γ(a)Γ(b)Γ(d-a-b))`. This was checked numerically against the
    `G10` and `G20` expansions of Masters.prc.
- Coefficients are rational in ε.

Shared p-finalizer (`p_integral::finalize::finalize(values, loops, external, settings, master_mode, backend)`, owned
by the Forcer agent, used by both backends; `backend` names the backend in `MissingMasterOrders` reports):

```text
I = Σ_w [F·i·π^(d/2)·εG(1,1)]^L · (-dot(p,p))^w · (-dot(p,p)/mursq)^(-Lε) · mursq^(-Lε) · reduced_w
```

- The explicit `(Q^2/μ^2)^(-Lε)` factor is kept separate from `(μ^2)^(-Lε)`. This mirrors the vacuum's
  `logmUVmu`/`log_mu_sq` split: normalizations that carry `exp(log_mu_sq)^(ε L)` (MSbar, GScheme) cancel
  `mursq^(-Lε)` symbolically.
- Master modes are identical to RustRed's (`RustRedMasterMode`):
  - `Retain`: exact, masters symbolic, no expansion;
  - `Symbolic`: expanded, `PG` expanded exactly, `PMaster` replaced by its exact Masters.prc coefficients with the
    constants symbolic (`zeta(n)`, `MZV(...)`), as RustRed's `Symbolic` mode substitutes exact master formulas;
  - `Numerical`: `PG` from exact Γ series; `PMaster` from the table, with `Oep` tails.
- `Symbolic` and `Numerical` expand in ε through `number_of_terms_in_epsilon_expansion - L - 1` and finish with
  `truncation_guard::ensure_known_orders`. The numerical constants are evaluated by the caller's numerical boundary
  (`EvaluationMethod::numerical_output`), as for the other backends.
- The expansion (`p_integral/laurent.rs`) writes every master monomial as
  `rational(ε) × exp(Σ log Γ(1+cε) + scale logarithms) × PMaster(ε)` and multiplies the three Laurent series with
  expanded coefficients, each to the depth the others' valuations require. Symbolica's generic series of the Γ
  products took minutes for one 4L ε^0 coefficient; this takes milliseconds.
- Logs appear as `log(-dot(p(i),p(i))/mursq)`, always in dot notation. The polynomial factors follow
  `use_dot_product_notation`.

New normalization `LoopNormalizationFactor::GScheme`:
- The per-loop factor is `(i·π^(d/2))^(-1) · exp(log_mu_sq)^ε / (ε G(1,1))`, built in `to_atom` from exact Γ
  functions. Custom strings still exclude Γ.
- With it, `I = (-1)^(Σa+r)·(Q^2)^w·(Q^2/μ^2)^(-Lε)·C_G`: the Forcer G-scheme value at `Q^2 = μ^2`.
- The default normalization stays `pySecDec`, as for vacuum.

p-master table (`PMaster` values):
- It is parsed once (`OnceLock`) from the vendored `form_src/forcer/forcer/Masters.prc` (the verbatim upstream file
  is the single source).
- `id Master(name) = Σ ep^k*(...)`, where `+ep^K*(ALARMname)` becomes `ε^K·Oep(K, PMaster(name))`.
- Constant map:
  - `z2..z12` -> `zeta(n)`;
  - `z5z3` -> `MZV(5,3)`, `z7z3` -> `MZV(7,3)`, `z5z3z3` -> `MZV(5,3,3)`, `z9z3` -> `MZV(9,3)`,
    `z6z4z1z1` -> `MZV(6,4,1,1)`.
- `vakint::MZV(n_1,...,n_k)` is a new built-in constant, `Σ_{m_1>...>m_k>=1} ∏ m_i^(-n_i)` (Forcer/datamine
  convention). Its exact definition (inspection) is the GPL `(-1)^k G(0^(n_1-1),1,...,0^(n_k-1),1; 1)`. Symbolica's
  GPL evaluator stops at weight six, so numerically `MZV(...)` becomes the evaluable terminal
  `MZVvalue(mzv_indices(...))` of `p_integral/mzv.rs` (Hölder convolution at 1/2, any weight and precision).
- Pinned by tests: `MZV(5,3) = 0.03770767298484754401130478229...`, and `ζ(5)ζ(3) = MZV(5,3) + MZV(3,5) + ζ(8)`
  numerically.
- The table is built in (the Forcer backend must work without external data). Masters.prc is under the vendoring
  permission.

Exact backend comparison:
- Retain-mode results are compared coefficient by coefficient over the master monomials (products of
  `PMaster`/`PG`): `(c_A - c_B).together() == 0`.
- Numerical results are compared at `p(i) = (0,0,0,1)` (`Q^2 = 1`) and `mursq = 1`, to `1e-30` relative.
- Both backends share the same finalizer, so a disagreement is a reduction disagreement.

## 5. Backend contract

Each backend implements:

```rust
fn reduce_p(family: &PFamilyData, powers: &[i64], canonical_numerator: AtomView,
            external: &Atom, settings: &VakintSettings) -> Result<Vec<ReducedPropagatorValue>, VakintError>
```

Inputs:
- `powers` is in canonical line order (0 = pinched, all >= 0 after normalization).
- `canonical_numerator` is a scalar polynomial in canonical `k(j)` and the user head.
- `external` is the canonical `p(1)` written in the user head, `σ·p(i)`.
- These are independent of `ReplacementRules`, so the backends can be built in parallel with the frontend. The
  glue reads them from the witness: `ReplacementRules::propagator_request(numerator)` returns the
  `PropagatorRequest { family, powers, canonical_numerator, external }` that `Vakint::forcer_evaluate` passes to
  `reduce_p` and the finalizer.

Results are exact in ε. Expansion, normalization, master substitution and the `Oep` guard happen only in the shared
finalizer.

### Forcer backend (`EvaluationMethod::Forcer(ForcerOptions)`; `ForcerOptions { master_mode }`)

Wiring:
- `supports`: `Topology::Propagator` with `1 <= L <= 4`.
- `dependencies() = [VakintDependency::Forcer]`: FORM at `form_exe_path` with version `>= 5.0.0`
  (`MINIMAL_FORCER_FORM_VERSION`), probed once per process like FORM today.
- `EvaluationOrder::all()` = [RustRed, Forcer, AlphaLoop, MATAD, FMFT, PySecDec]. `available()` drops Forcer without
  FORM 5. Add `EvaluationOrder::forcer_only()`. Python accepts `"forcer"`.
- `CitationSource::Forcer` cites arXiv:1704.06650 and the Zenodo DOI.
- An explicit Forcer request without FORM 5 is the new
  `VakintError::ForcerUnavailable { form_exe_path, reason }`:
  > "the Forcer backend (massless propagator integrals) needs FORM >= 5.0.0, but '{form_exe_path}' {reason}.
  > Install FORM 5 (https://github.com/form-dev/form) and expose it through PATH, FORM_PATH or form_exe_path."
  The same error replaces `NoEvaluationMethodFound` when no method of the order evaluates a p-integral and FORM 5 is
  missing (Forcer is the only p-integral backend until RustRed ships p data).

Vendoring:
- `crates/vakint/form_src/forcer/`: verbatim upstream `forcer.h` and the 19 files of `forcer/` from
  benruijl/forcer `68a30bf` (`68a30bfb2aec3672a2f674ca8666bfc970105e24`), unmodified. vakint never calls
  `gschemeexpand` or `msbarexpand`, so the depth patch of the study copy is not needed.
- `build.rs` embeds the files gzip-compressed (flate2 build-dependency; ~1.1 MB instead of 9.6 MB).
- At runtime they are materialized once per process into `<temporary_directory|TMPDIR>/vakint-forcer-68a30bf-<uid>/`
  (atomic rename) and symlinked into each FORM run directory. FORM includes `forcer.h` from there, so the directory
  is reused only when it is a private (mode 0700, not a link) directory of the current user holding byte for byte
  the embedded sources; anything else (a directory planted by another user of a shared /tmp, a damaged copy) is
  ignored and the sources go to a fresh private directory of the process.

`form_src/forcer/NOTICE` (exact content):

```text
Forcer: a FORM program for massless four-loop propagator integrals.
Authors: Ben Ruijl, Takahiro Ueda, Jos Vermaseren (ERC Advanced Grant 320651, HEPGAME).
Upstream: https://github.com/benruijl/forcer, commit 68a30bfb2aec3672a2f674ca8666bfc970105e24 (vendored unmodified).
Reference: B. Ruijl, T. Ueda, J.A.M. Vermaseren, arXiv:1704.06650; Zenodo DOI 10.5281/zenodo.1043944.
Masters.prc expansions follow R.N. Lee, A.V. Smirnov, V.A. Smirnov, arXiv:1108.0732.
The upstream repository carries no LICENSE file. These sources are redistributed inside vakint with the
permission of the Forcer authors (Ben Ruijl, collaborator), granted on 2026-11-03.
```

Input generation:
- The input is the uncontracted (cubic) family graph, a top-level Forcer topology: Forcer momenta `p1..pN` are all
  the family lines in canonical order, oriented as the canonical edges. `Q` is the canonical `p(1)`.
- One `vx(...)` per family vertex: `+p` for incoming lines, `-p` for outgoing lines, `Q` at vertex 1 and `-Q` at
  vertex 2.
- Propagators are written `1/(pj.pj)^a` for the positive-power lines only. A pinched line keeps its place in the
  vertices without a propagator, and Forcer's own reduction hands the sector over to its subtopology.
- Never contract the pinched lines: Forcer matches the input only against its own topologies, whose vertices carry
  at most four momenta, and silently drops any other input (a wrong zero). Contracted sectors have vertices with
  five or six momenta for 3700 of the 4703 non-scaleless sectors of the registry (3L: 97, 4L: 3606).
- Fail closed: Forcer's warning `the following diagram did not match a topology` on standard output is the error
  `VakintError::ForcerUnrecognizedTopology`, never a zero.
- The canonical `k(j)` are solved in line momenta and `Q` (a pinched loop line is solved through the remaining
  lines).
- The numerator becomes FORM dot products (`p1.p3`, `Q.p2`, `Q.Q`). It is split by degree `r` into separate
  `Local F{r}` expressions, because Forcer sets `Q.Q = 1`.
- Template `templates/run_forcer.txt`: `#include- forcer.h`, the locals, `#call Forcer()` (exact, no expansion),
  then `#write` per local.
- Run `form` when `n_cores = 1`, otherwise the sibling `tform -w{n_cores}` (fall back to `form`).

Output parsing:
- `rat(a,b)` -> `a/b`, `ep` -> ε, `GschemeConstants(x1,x2)` -> `PG(x1,x2)`, `Master(n)` -> `PMaster(n)`.
- Any other symbol is a `FormOutputParsingError`.
- Each part becomes `ReducedPropagatorValue { scale_power: 2L-Σa+r, reduced: (-1)^(Σa+r)·parsed }`.
- Measured (FORM 5.0.0, one core): a 4L `haha` with a rank-2 numerator, exact, takes 5.1 s and 31 MB.

### RustRed backend

`EvaluationMethod::RustRed`, unchanged options:
- `supports` for `Topology::Propagator` is true when the p data of that loop count is loaded: 1L-3L built in, 4L from
  external data. Otherwise it yields to Forcer, exactly as 4L vacuum yields to FMFT.
- It returns the same `ReducedPropagatorValue`s; see section 6.

## 6. RustRed backend outline (M2; no implementation in M1)

vakint side:
- New `src/rustred_evaluation/p_integral.rs`.
- Numerator lowering generalizes `MatchedScalarInput::prepare` to families with coordinates `k_i·k_j` and `k_i·q`.
  It produces an `IntegralKey` over the family slots (lines, then ISPs; ISP numerators become negative powers).
- One saved application per loop count ("p<L>L"). All families of that loop count share its owner classes.
- Calls into `rustred_app` per request: `PreparedSavedApplication::<N>::route_entry_to_owner`, then
  `reduce_routed_linear_combination` under the Ordered scheduler (`ordered_scheduler(settings)`), with
  `OrderedWeightArithmetic::ModularPolynomial` (Z_p(d) replay) by default.
- All p contributions of one loop count in a request share one frontier, as for vacuum.
- Arity: `N = 2, 5, 9, 14` for L = 1..4, dispatched over RustRed's compiled arities by a const-generic match. Never
  hard-code 16; 5L needs 20.

Family identity:
- vakint `P_FAMILIES[i]` maps to the RustRed family compiled from the same denominators: Minkowski, massless,
  `externals(q)`, `gram(q,q,-1)` (Euclidean `Q^2 = 1`, field Q(d)), slots in vakint order.
- The data manifest lists `name -> fingerprint`. vakint validates it at load, and a test recomputes it.
- vakint never handles owner ids. RustRed routes a family sector to its owner class (marked-matroid class, shared
  across families) through signed affine maps with q -> ±q.

Shareable cache:
- One Ordered step cache per p program. Its records are keyed by (owner class, owner-local `IntegralKey`, RustRed
  step-semantics version), so they are independent of the family and the route the request came from.
- Directory `p<L>L-<program>.<rule tables>.<stamp>`, with version constants and no digests of sources. Shipped
  caches use the same layout.

Data layout:
- Built-in `crates/vakint/data/p1L`, `p2L`, `p3L`; `p4L` in vakint-external-data. The release manifest gains
  `propagator_packs`, with the same `LoopPack` fields.
- Loop manifest `vakint.external-data.loop.v1` with:
  - `"integral_kind": "propagator"`, `"external_momenta": ["q"]`, `"gram": [["-1"]]`;
  - coordinates of kind `loop_loop` and `loop_external`;
  - `"normalization": "Forcer-GScheme-Euclidean-Q2=1"`;
  - terminal records whose `exact_basis_expression` is written in the common basis (`PMaster`, `PG`, exact in d);
  - `ibps/` packages and `ordered_cache/` as for vacuum.
- `load_official_external_data(4)` also installs `p4L` when present.

What RustRed must provide (for the RustRed agents):
1. Families with one external momentum (`externals(q)`, numeric `gram(q,q,-1)`) through every vacuum-only gate of
   plan section 3.3: source port, app admission, scalar-numerator service with `k·q` coordinates and q^2
   homogeneity, terminal normalization, relations and products. Fingerprints stay stable (`E1;1:q`).
2. Sector symmetry that colours the two external vertices alike. q reversal is a symmetry with q -> -q.
3. A zero analyzer that agrees with vakint's scaleless rule on every sector of the 16 families (test).
4. Saved application packages for p owners: closure by construction, routes with q -> ±q and Whitney twists, the
   Ordered engine, Z_p(d) replay and an `OrderedStepCacheBinding` per program.
5. An analytic insertion lane:
   - recursively one-loop sectors return exact `PG` products in Forcer's canonical `GschemeConstants` form, with
     the same product identities Forcer applies;
   - insertion masters with ε-power lines map to `PMaster`.
6. Exact-in-d terminal -> master relations onto the common basis: built-in for 1L-3L, the M4 Kira basis change for
   4L.
7. Runtime arities 2, 5, 9 and 14 now, and 20 for 5L.
8. `SavedApplicationReduction` terms with exact Q(d) coefficients, plus the unresolved set. No FORM.

## 7. Test plan

Rules:
- Every test is one short call of a shared helper in the new `tests/support/p_integrals.rs`. Helpers:
  `forcer_settings()`, `rustred_p_settings()`, `assert_p_backends_agree(input, laurent_terms)` (exact Retain
  comparison plus Numerical to 1e-30), `assert_p_routing_invariant(input, seeds)`, `assert_p_dot_round_trip(atom)`,
  `assert_p_scaleless(input)`, `assert_p_error(input, fragment)` (with `allow_unknown_integrals: false`),
  `assert_p_malformed(input, fragment)` (an error under either setting), `assert_p_unknown(input, fragment)` (section 1,
  "Out of domain": `topo(UNKNOWN(input))` by default, the error otherwise, never a value),
  `assert_p_identity(lhs, rhs)` (exact live Forcer equality, symbolic masters, G-scheme) and `fake_form(version)` (a
  FORM that only reports its version).
- Forcer battery data are frozen in `tests/fixtures/forcer/p{L}l_battery.json`: the canonical input (family, powers,
  numerator, external), the raw exact Forcer output of each homogeneous part, FORM version and Forcer commit.
  Without FORM, they test the input generation, the parser, the finalizer and RustRed.
- With FORM 5 (`VAKINT_TEST_FORM5=/path/to/form`, or the default `form` when it is FORM 5, as in the dev shell and
  CI) every case also runs live and must reproduce the fixture; `VAKINT_FORCER_FIXTURES=write` regenerates them.
- Fast gate (`just test_vakint`): everything at 1L-3L plus a few frozen 4L cases. Slow 4L RustRed cases go in the
  `test_vakint` profile's default-filter (never in `slow::` modules); `test_vakint_full` runs everything.

Per loop count:
- 1L (`P1L`):
  - powers `a, b` in -2..3 and numerators `dot(k,q)^n`, `n <= 4`: exact against the closed G-function form and
    against Forcer;
  - scaleless cases (`a <= 0` or `b <= 0`) are exactly 0.
- 2L (`P2L`): every sector with dots 0-2 and rank 0-4, symbolic and numerical against Forcer; the kite gives
  `6ζ(3)` at ε^0 (G-scheme).
- 3L (`LA`, `BE`, `NO`): tops and sectors with rank 0/2/4 and dots 0-2. The 3L reductions produce
  `PMaster(no)` and `PMaster(t105)` (measured), which need their Masters.prc orders.
- 4L (11 tops): scalar, rank 2 and one dot per top. Benchmark-grid subset in the full suite; `haha` rank 2 in the
  gate (frozen).
- Every loop count:
  - Routing invariance: random unimodular `k -> Mk + c·q` maps, vertex relabelings, edge flips, q reversal
    (`p(1) -> -p(1)` in every prop) and head renaming (`p(1) -> p(3)`) give identical Retain results.
  - DOT round trip of every top and of three pinched sectors per loop count.
  - `Oep` guard, at 3L and 4L (1L-2L have only exact `PG`): a case where a `1/ε^k` coefficient pulls a `PMaster`
    tail into the requested order gives `MissingMasterOrders`, and the same request one order lower passes.
  - Error cases of section 1: massive line, two heads, broken conservation, non-integer power, tensor numerator,
    5L input, non-unimodular routing; the same out-of-domain inputs are `UNKNOWN` under `allow_unknown_integrals`, the
    malformed ones stay errors, and unknown vacuum integrals are unchanged.
  - Normalization: tree line, merged identical quadrics, negative power; evaluated (live Forcer): tree lines of power 1
    and 2, a negative short-form power, a reversed short-form head, and a bubble chain through a cut vertex (`1/ε^2`,
    not scaleless).
  - Analytic identities (live Forcer, `tests/p_integral_forcer_tests.rs::identities`): odd and even rank and the
    reversed external momentum (`σ = -1`) at 1L-4L (`k·q = (k² + q² - (k-q)²)/2` and its partial fractions on the
    bubble, the kite, `NO` and `LA4`); index notation and linear `dot` arguments (also FORM-free on canonical forms).
  - The `Oep` guard through `Vakint::evaluate` at 4L (`P4L_NO1`, symbolic and numerical masters).
  - G-scheme on vacuum integrals (FORM-free, RustRed: the massive tadpole is `1/ε - 2 + ζ2 ε`) and on a vacuum plus a
    p-integral (FORM 5).
  - FORM older than 5 (fake FORM 4.3.1, FORM-free): `EvaluationOrder::all()` validates, a vacuum integral evaluates, a
    p-integral reports `ForcerUnavailable`, and an order of Forcer alone fails validation.
  - The MZV evaluator (`src/p_integral/mzv.rs`): stuffle relations at weights 8-12 and duality for every `ζ(n)` and
    `MZV(...)` of the master table, directly and through the numerical evaluator.
  - Forcer sources: a foreign or modified materialized directory is never used
    (`forcer::tests::untrusted_source_directories_are_never_used`).
- Registry and catalogue:
  - the `src/p_integral.rs` data tests;
  - catalogue determinism;
  - census: every 1PI cubic p-multigraph with <= 4 loops is covered after normalization;
  - the scaleless rule equals RustRed's zero analyzer on all sectors (M2).
