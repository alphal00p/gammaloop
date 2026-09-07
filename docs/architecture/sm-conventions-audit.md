# Standard Model conventions audit

This records the local SM conventions used by GammaLoop after the phase-convention corrections. The canonical UFO and JSON now agree on the corrected ghost derivative and Goldstone Yukawa rules, using an explicit scalar/pseudoscalar and vector/axial basis. The three-gluon, electroweak triple-vector and derivative scalar rules were retained after independent action and component checks. Fermion and antifermion propagator signs were retained. General electroweak gauge consistency remains a separate limitation, described below.

The full amplitude and LU marking/contour contract is in [phase-conventions.md](phase-conventions.md). The supported setting is a Hermitian theory with real propagator masses. The local audit does not establish a complex-mass prescription or complete electroweak LU correctness.

## Momentum, spin and color conventions

We use metric `(+,-,-,-)`, fields Fourier transformed with `exp(-ip.x)`, and a vertex equal to `i` times the interaction action. Numbered UFO `P(mu,leg)` is the incoming momentum of that vertex leg. GammaLoop assigns the underlying source-to-sink momentum `q` with `+q` at a sink and `-q` at a source. This conversion is owned by [UFO reindexing](../../crates/gammalooprs/src/numerator/ufo.rs), after vertex legs, tensor slots and momenta have received the same permutation in [graph parsing](../../crates/gammalooprs/src/graph/parse/mod.rs).

These choices agree with [UFO 2.0, Table 2](https://arxiv.org/pdf/2304.09883). The action checks used [Romão–Silva, sections 2, 4 and 5](https://arxiv.org/pdf/1209.6213), whose sign parameters expose the convention choices; [Denner, Appendix A](https://arxiv.org/pdf/0709.1075), which fixes incoming momenta and propagator factors; and [Schwinn, equations 4.81–4.87](https://www.tep.physik.uni-freiburg.de/lectures/archive/QCD-WS-14/qcd), which derives the QCD and Faddeev–Popov interactions. ALOHA compatibility was not the criterion for accepting a sign.

The electroweak bosonic rules fix

```text
D = partial - i K,
K = [[ e A + (g/cw)(1/2-sw^2) Z,   g W+/sqrt(2) ],
     [ g W-/sqrt(2),              -g Z/(2cw)    ]],
Phi = ( i G+, (v+H-i G0)/sqrt(2) )^T,
g = e/sw.
```

This Goldstone phase follows from the actual nonderivative bosonic vertices. Differentiating the action checked all 31 scalar/Higgs VVS, VVSS, SSS and SSSS entries and then all seven VSS currents. An isolated reversal of the Z–G0–H or W–G–H coupling would break this agreement.

The QCD spin and fundamental-color index orders must be distinguished. With the canonical `Gamma(3,1,2)` and `T(3,2,1)`, the physical quark color matrix is transposed relative to its spin matrix. Equivalently `R^a=-T^{aT}` satisfies `[R^a,R^b]=i f^{abc}R^c`, with a matching `D=partial+igR` and `F=dA-gfAA` convention. A sign comparison that reverses only one color or momentum ordering is not comparing the same rule.

## Complete derivative inventory

Write `S=P_2-P_3` for `VSS1`, `U=p_antighost` for the corrected `UUV1`, and

```text
K123 = g12 (p1-p2)_3 + g13 (p3-p1)_2 + g23 (p2-p3)_1
```

for `VVV1`. All momenta and field lists below are incoming; color and Lorentz leg numbers refer to the same ordered list. The table contains all 23 momentum-dependent vertices in the canonical SM, using its three derivative Lorentz structures. `c` denotes a ghost and `barc` its antighost; the W label identifies the model species.

| Vertex | Ordered fields | Current coefficient and structure |
|---|---|---|
| V_11 | A,G−,G+ | `−ie S` |
| V_12 | cA,barcW−,W− | `−ie U` |
| V_13 | cA,barcW+,W+ | `+ie U` |
| V_15 | cW−,barcA,W+ | `−ie U` |
| V_18 | cW−,barcW−,A | `+ie U` |
| V_19 | cW−,barcW−,Z | `+ig cw U` |
| V_21 | cW−,barcZ,W+ | `−ig cw U` |
| V_23 | cW+,barcA,W− | `+ie U` |
| V_26 | cW+,barcW+,A | `−ie U` |
| V_27 | cW+,barcW+,Z | `−ig cw U` |
| V_29 | cW+,barcZ,W− | `+ig cw U` |
| V_31 | cZ,barcW−,W− | `−ig cw U` |
| V_33 | cZ,barcW+,W+ | `+ig cw U` |
| V_35 | cG,barcG,g | `−gs f(3,1,2) U` |
| V_36 | g,g,g | `−gs f(1,2,3) K123` |
| V_41 | W−,G0,G+ | `+ig/2 S` |
| V_42 | W−,G+,H | `−g/2 S` |
| V_43 | A,W−,W+ | `+ie K123` |
| V_47 | W+,G0,G− | `−ig/2 S` |
| V_48 | W+,G−,H | `−g/2 S` |
| V_54 | W−,W+,Z | `+ig cw K123` |
| V_57 | Z,G0,H | `−g/(2cw) S` |
| V_58 | Z,G−,G+ | `−ig(cw^2-sw^2)/(2cw) S` |

The final export certificate reads the current UFO `lorentz.py`, `vertices.py` and `couplings.py` directly through their Python syntax trees and compares them with the current JSON. All 23 field orders, color factors, Lorentz and coupling slots, symbolic coupling expressions and three derivative tensors agree. Historical pre-correction inventory files are not the source for this comparison.

The correction changes `UUV1` from `−P(3,2)` to `+P(3,2)` in both [the UFO](../../assets/models/ufo/sm/lorentz.py) and [the JSON](../../assets/models/json/sm/sm.json). No derivative coupling expression or cached coupling value changes. For example, the charged ghost kinetic term and scalar electromagnetic current fix `D=partial-ieA`; `−barc partial D c` then gives V_26 equal to `−ie pbar`. The previous `UUV1` reversed that result. The same action comparison fixes all twelve electroweak ghost currents and the QCD ghost current.

The three-gluon sign is consequently valid **with the conventions above**. An independent cubic Yang–Mills differentiation checked all 64 Lorentz components for each ordered color basis of V_36, V_43 and V_54. A separate exact `q qbar -> gg` Ward calculation used the SU(2) subgroup of SU(3), `f123=+1`, fixed rational on-shell momenta, one transverse spectator polarization, four independent external-spin choices and all entries of each resulting color matrix. It vanishes for the actual transposed quark color and current three-gluon sign; changing either independently gives four nonzero entries. This is a discriminating relative-sign certificate, not a claim to have sampled every SU(3) color assignment or every kinematic point. Changing both conventions together is an equivalent representation. No three-gluon correction is needed.

## Scalar/pseudoscalar and vector/axial basis

Let `I` be the spin identity, `G5=gamma5`, `V^mu=gamma^mu` and `A^mu=gamma^mu gamma5`, in that order. The final canonical Lorentz definitions are:

| Structure | Final expression |
|---|---|
| FFS1 | `(I+G5)/2` |
| FFS2 | `G5` |
| FFS3 | `(I-G5)/2` |
| FFS4 | `I` |
| FFV1 | `V` |
| FFV2 | `(V-A)/2` |
| FFV3 | `(-V-3A)/2` |
| FFV4 | `(3V+A)/2` |
| FFV5 | `(5V+3A)/2` |

The physical corrections are the exchange of the pure-scalar FFS1/FFS3 chiralities and the reversal of the FFS2 pseudoscalar sign. FFS4 remains the CP-even Higgs identity. The vector entries only change representation; their original left/right weights and couplings are preserved. The canonical SM Lorentz structures now contain no `ProjM` or `ProjP` calls.

For the bosonic phase fixed above, the Yukawa action requires

```text
G+ ubar d : V_ij (yd PR - yu PL),
G- dbar u : V_ij* (yu PR - yd PL),
G0 ubar u: +yu/sqrt(2) (PR-PL),
G0 dbar d: -yd/sqrt(2) (PR-PL),
PL=(I-G5)/2, PR=(I+G5)/2.
```

Here `V_ij` is defined from the model's existing W+ flavor coefficient; no CKM naming convention or value is changed. This action check covers all 26 Goldstone–fermion vertices symbolically, including complex CKM and both Yukawa terms where present. Six ordinary Higgs fermion couplings remain unchanged. The affected neutral vertices are V_77,101,102,103,138,139; the charged ones are V_82–88,110–112,116–122,145–147.

The sum `PL+PR=I` is gamma5-free; the difference `PR-PL=G5` is a genuine pseudoscalar. Writing this basis explicitly lets scalar/vector cancellations occur before dimensional Dirac algebra. It does not authorize discarding an uncanceled axial term or imposing a new d-dimensional gamma5 prescription. [Idenso's existing simplifier](../../crates/idenso/src/dirac/simplify.rs) restricts gamma5 anticommutation and special-factor reductions to their supported dimensions. The tests `gamma5_does_not_move_in_dimension_generic_chain` and `gamma5_square_is_four_dimensional` preserve that boundary. Normal numerator processing, analytical amplitude evaluation and integrated UV preparation already invoke the same gamma simplification owner before Vakint.

The model regression `sm_lorentz_structures_use_the_scalar_and_vector_axial_basis` in [model/mod.rs](../../crates/gammalooprs/src/model/mod.rs) checks the exact basis coefficients, scalar sum/difference and absence of chiral-projector symbols. The audit-time independent UFO/JSON comparison checked all 22 Lorentz structures, 768 explicit matrix entries and all 26 Goldstone action identities.

## Generated local certificates and RHS sewing

[sewing_tests.rs](../../crates/gammalooprs/src/processes/cross_section/sewing_tests.rs) uses actual `Graph::from_string` vertex assignment and tensor contraction, with independent explicit matrices/spinors as the oracle.

`generated_sm_charged_ward_and_ghost_momentum_follow_the_action` fixes `MW=5`, `m_tau=2`, `g=1`, `y_tau=sqrt(2)/5` and real on-shell momenta

```text
k_W=(5,0,0,0),
p_nu=(21/10,0,0,21/10),
p_tau=(29/10,0,0,-21/10).
```

With normalized spinors, the fixed left-handed neutrino gives `k.A_W=-i sqrt(42)` and corrected `A_G=sqrt(42)/5`, so `k.A_W+i MW A_G=0`. The Rust regression divides the spinors by their common normalization factors, using `ubar=(0,0,0,1)` and `v=(0,-2/5,0,1)`, and writes the Yukawa as `(2/5)/sqrt(2)`. This homogeneous identity then reduces exactly without a radical-simplification assumption. It first checks every generated scalar and temporal-vector matrix entry against explicit Weyl matrices. The previous scalar chirality gave `A_G=0`, violating the identity. This tests chirality rather than only a spin-summed norm. No internal vector propagator is involved.

The same test checks the all-incoming charged ghost vertex with `e=1` and `pbar^0=3`, obtaining `−3i`. Actual graph exports before the correction instead gave `+3i`, with `overall_factor=1` and the antighost momentum attached to the expected edge. There is no hidden per-vertex graph-ordering minus that could compensate the old `UUV1`.

`generated_charged_scalar_forward_vertex_is_already_the_hermitian_partner` compares the inverse-process vertex with the explicit matrix adjoint and independently sewn spin sums. `generated_complex_charged_current_preserves_the_ckm_norm` retains a genuinely complex CKM coefficient and tests the sewn norm. These establish why inverse UFO vertices on the RHS must not be subjected to a second whole-tensor conjugation. The LU marking and contour factors remain separately required by [the phase contract](phase-conventions.md).

The neutral `f fbar -> ZH/G0H` audit also evaluated all four tree diagrams with exact rational components: consistent field phases satisfy the Ward identity, whereas the previous scalar-current/Goldstone-Yukawa pair does not. That calculation explicitly uses common Feynman-gauge propagators. It is an action diagnostic, not a literal complete amplitude from the still-hybrid canonical JSON.

## Input-parameter reality

The vendored UFO and canonical JSON agree on the names, nature and type of all 72 SM parameters. All 26 external inputs are declared `real` and have zero imaginary defaults: three SMINPUTS, four Wolfenstein inputs, six YUKAWA inputs, eight MASS inputs and five DECAY inputs. Derived CKM entries retain their complex declarations, allowing their physical phases when the real Wolfenstein inputs are nonzero. No parameter-declaration correction was needed.

The loader preserves these declarations as `ParameterType`, which the API and CP-optimization guard use. It does not automatically attach a Symbolica reality attribute to the underlying symbol. The marking prescription uses the model's inverse-process vertices and does not depend on such an attribute to undo a second numerator adjoint.

## Fermion flow and Grassmann statistics

For an underlying source-to-sink momentum `q`, the actual propagator ledger is:

| Particle | Numerator | Spin slots |
|---|---|---|
| Positive-PDG fermion | `i(slash(q)+m)` | `(sink,source)` |
| Negative-PDG antifermion | `i(-slash(q)+m)` | `(source,sink)` |

Joint reversal changes the endpoints, momentum and PDG together, preserving the physical propagator. For future-directed `p`, `u ubar=slash(p)+m` and `v vbar=slash(p)-m`, so the negative-PDG numerator is `−i v vbar`, not `+i v vbar`. Its positive mass term is correct. All 24 spinor propagator exports were checked; explicit massive Weyl components and inverse-Dirac-operator identities independently confirm both signs. A scalar-current sewing trace gives the same positive norm as the explicit spin sum after the graph's fermionic statistics sign is included.

The generic `Particle::is_anticommutating` predicate and corresponding loop/external-ordering logic have been restored from the previously implemented upstream change. Ghosts as well as fermions receive one minus per closed Grassmann loop. This correction is independent of UUV1: a two-vertex ghost loop contains two UUV1 factors, so reversing that vertex cannot supply its missing loop minus. The public anticommutating-loop filter uses the same predicate. [Feyngen tests](../../crates/gammalooprs/src/feyngen/test.rs) cover closed topologies, identical external ghost exchange and an actual generated ghost loop (`closed_anticommutating_loop_counts`, `open_ghost_chain_exchange_has_grassmann_sign`, `generated_ghost_loop_has_one_statistics_minus`).

## Remaining electroweak gauge limitation

The canonical JSON still stores

`D_U^{mu nu}=-i(g^{mu nu}-q^mu q^nu/M^2)/(q^2-M^2+i0)`

for W/Z while retaining finite-mass Goldstone and electroweak ghost interactions. This is a mixture of unitary-form vector propagators and the xi=1 unphysical spectrum. In a consistent Feynman gauge the virtual vector numerator is `−ig^{mu nu}` and the Goldstone/ghost masses equal M; [Romão–Silva, equations 51–57](https://arxiv.org/pdf/1209.6213) and [Denner, Appendix A, equation A.2](https://arxiv.org/pdf/0709.1075) specify these together.

Generation does not automatically exclude Goldstone or electroweak ghost vertices. The `vector_polarization_sum_gauge` setting controls external initial-state sums, not virtual propagators. Moreover, final LU cut edges currently retain their model propagator numerator. A physical massive-vector cut requires the three-polarization sum `−g+qq/M^2` even when virtual vectors use Feynman gauge. Consequently a future virtual-propagator correction must also preserve the physical cut projector or provide a certified complete covariant cut-state sum; replacing three JSON numerators alone does not establish general W/Z LU correctness.

Fresh UFO import must be included in that future correction: the current loader's default massive-vector branch also creates the unitary numerator when called with `Feynman`. The canonical JSON additionally lacks the Goldstone metadata for G0/G±, so a proposed filter based on `is_goldstone()` must first restore that metadata. A coherent alternative is an explicitly selected electroweak unitary model that excludes Goldstone and electroweak ghost contributions while retaining QCD ghosts. Neither choice is silently imposed by the local sign corrections.

This limitation does not affect the scalar phase fixtures or the four requested photon-mediated ddbar/ttbar LO/NLO QCD acceptances: their allowed-particle sets exclude virtual and cut W/Z and all electroweak unphysical fields. The local charged Ward test likewise has no internal gauge propagator. General electroweak amplitudes and physical W/Z LU final states require the separate gauge contract to be resolved and tested.
