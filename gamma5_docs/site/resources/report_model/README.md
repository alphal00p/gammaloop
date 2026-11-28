# HepKIT model of the regularised unbroken SM (BMHV, sterile partners) and one-loop 1PI enumeration

Conventions: `report/sections/03_setup.typ` (labels quoted as `@eq:...`). Tool: `symbolica.community.hepkit`
(Symbolica 3.0.x, venv `<scratch>/venv`).

## Files

| file | content |
|---|---|
| `build_model.py` | builds `sm_unbroken_bmhv.json` (commuting ghost encoding) and, with `--native-ghosts`, `sm_unbroken_bmhv_native.json` |
| `sm_unbroken_bmhv.json` | the model: 61 particles, 174 vertex rules, 54 couplings, orders QCD GW GY YUK LAM DELTA SRC |
| `classes.py` | Green-function classes, multiplet-level assignments allowed by rigid symmetry, component tuples, process construction |
| `enumerate_diagrams.py` | generation of all one-loop 1PI diagrams (writes `parts/<class>.json`, `dot/<class>/<assignment>/<tuple>_D<k>.dot`) |
| `merge_and_report.py` | builds `enumeration.json` and `ENUMERATION.md` (multiplet-level topologies, sanity checks) |
| `tests/t_model_consistency.py` | checks the JSON couplings against the group theory and Feynman rules of 03_setup (passes) |
| `tests/t_st_signs.py` | zero-dimensional Grassmann check of `@eq:ST`, `@eq:gf-exact` against `@eq:Lext` (finds two sign errors in the draft) |
| `tests/t_speed.py` | timing reproduction for the `max_vertices` issue |

Run (licence variable set inline):
```
python build_model.py && python build_model.py --native-ghosts
SYMBOLICA_LICENSE=... python enumerate_diagrams.py --jobs=12      # about 20 s wall time
python merge_and_report.py
```

## Fields

| HepKIT particle | field of 03_setup | spin (UFO) | colour | notes |
|---|---|---|---|---|
| `Q1`,`Q2` | $\Psi_Q$, SU(2) components $i=1,2$ ($u_L$, $d_L$ partners) | 2 | 3 | Dirac package $\psi+\chi$, propagator $i\,\slashed k/k^2$ |
| `U`,`D` | $\Psi_U$, $\Psi_D$ | 2 | 3 | |
| `L1`,`L2`,`E` | $\Psi_L$ components, $\Psi_E$ | 2 | 1 | |
| `X~` | $\overline\Psi_X$ (leg of the barred field) | | | |
| `G` | $G^A_\mu$ | 3 | 8 | Feynman gauge $-i g_{\mu\nu}/k^2$ |
| `W1`,`W2`,`W3`,`B` | $W^I_\mu$, $B_\mu$ (real basis of @eq:Tg) | 3 | 1 | self-conjugate |
| `H1`,`H2`; `Hk~` | $H_k$; $H^\dagger_k$ | 1 | 1 | propagator $i/(k^2-M_H^2)$ with $M_H^2 \equiv m^2$ |
| `ghG`,`ghW1..3`,`ghB`; `~` | $\omega^a$; leg `~` = $\bar\omega^a$ | -1 | 8/1 | propagator $i/k^2$ |
| `rhoG`,`rhoW1..3`,`rhoB` | $\rho^\mu_a$ | 3 | 8/1 | external only |
| `zetaG`,`zetaW1..3`,`zetaB` | $\zeta_a$ | 1 | 8/1 | external only; `zetaB` has no vertex ($\mathfrak f^{Bbc}=0$) |
| `RQ1`..`RE`; `R*~` | $R_f$; leg `~` = $\bar R_f$ | 2 | 3/1 | external only |
| `sigma1`,`sigma2`; `~` | $\sigma_k$; leg `~` = $\sigma^\dagger_k$ | 1 | 1 | external only |

The antighost is internal only (antighost functions follow from the ghost equation `@eq:linear-ids`); $b^a$ is eliminated.

## Vertices (all of 03_setup; factor = coupling x colour x Lorentz, coupling = $i\times$ coefficient in $\mathcal L$)

| family | rules | Lorentz structure (UFO, legs as listed) | coupling |
|---|---|---|---|
| $\overline\Psi\Psi A$ `@eq:V-ffv` | 23 | `Gamma(idx3,idx1,d1) ProjM/ProjP(d1,idx2) * bmhv::bar(idx3)` ($\bar\gamma^\mu P_{c_f}$) | $i(\mathcal T^a_f)_{ij}$ |
| $AAA$ `@eq:V-vvv` | 2 (`GGG`, `W1W2W3`) | $g^{12}(p_1-p_2)^3+g^{23}(p_2-p_3)^1+g^{31}(p_3-p_1)^2$ | $\mathfrak f^{abc}$ |
| $AAAA$ `@eq:V-vvvv` | 1 + 3 | three metric pairings | $-i\,\mathfrak f\mathfrak f$ |
| $\omega\bar\omega A$ `@eq:V-ghost` | 1 + 6 | legs $(\omega^c,\bar\omega^a,A^b)$: `P(idx3, 2)` $=-p'$ | $\mathfrak f^{abc}$ |
| $H^\dagger H A$, $H^\dagger HAA$, $H^4$ `@eq:V-ssv`, `@eq:V-ssss` | 8, 14, 3 | $P_2-P_1$; $g^{\mu\nu}$; 1 | $i\mathcal T$; $i\{\mathcal T^a,\mathcal T^b\}$; $-2i\lambda(\dots)$ |
| Yukawa `@eq:V-yuk` | 12 | `ProjP`/`ProjM` | $\mp i y\,\epsilon_{ij}$, $-iy\,\delta_{ij}$ |
| sources `@eq:V-src-gauge`, `@eq:V-src-matter` | 78 (order SRC) | $\rho\omega$: `P(idx1,2)`; $\rho A\omega$: metric; $\zeta\omega\omega$: 1; $\bar R\omega\Psi$: $P_c$; $\overline\Psi\omega R$: $P_{\bar c}$; $\sigma^\dagger\omega H$, $H^\dagger\omega\sigma$: 1 | $\delta$, $i\mathfrak f$, $-i\mathfrak f$, $-\mathcal T$, $-\mathcal T$, $-\mathcal T_H$, $+\mathcal T_H$ |
| $\Delta_1$ `@eq:V-delta` | 23 (order DELTA) | legs $(\omega,\overline\Psi,\Psi)$: `Gamma(d2,idx2,d1) P(d2,3) Proj_cbar(d1,idx3) hat(d2) + Gamma(d2,idx2,d1) P(d2,2) Proj_c(d1,idx3) hat(d2)` | $(\mathcal T^a_f)_{ij}$ |

$P(\mu,k)$ is the incoming momentum of leg $k$, so $p=P_3$ and $p'=-P_2$ in `@eq:V-delta`; the instantiated numerator
reads $+\hat{\slashed Q}_{\rm in}P_{\bar c} - \hat{\slashed Q}_{\rm out}P_c$ with edge momenta (checked in a DOT file).
`bmhv::hat(i)` marks the evanescent dummy index and `bmhv::bar(i)` the four-dimensional gauge index; both pass the lowering
with their index renamed to the diagram labels, so a post-processor can move them to the hat/bar spaces.
Colour: native `T`, `f`, `Identity`. SU(2): explicit components with numerical Pauli matrices and $\epsilon^{IJK}$.

## Encoding choices and deviations from a faithful model

1. **Ghosts and sources are commuting (`ghost_number: 0`)** in `sm_unbroken_bmhv.json`. HepKIT treats every field with
   non-zero ghost number as Grassmann odd and rejects vertices with three odd legs ($\Delta_1$, $\bar R\omega\Psi$), and
   rejects external $\rho\,\omega$ chains. Consequences: (a) closed ghost loops lack the factor $-1$; it is recorded as
   `hand_closed_ghost_loop_sign`; (b) the ordering signs of ghost and source legs are not computed; HepKIT's
   `ExternalFermionOrderingSign` refers to spinor legs only and treats $R_f$ (bosonic in 03_setup) as a fermion and
   $\sigma$ (fermionic) as a boson. Signs of diagrams with external ghosts or sources must be fixed by hand. The
   `_native` model (true ghost numbers) is used as a cross-check where HepKIT accepts it (ghost-number-0 classes without
   sources): 149 component tuples agree in diagram count and in the ghost-loop sign.
2. **Real $W^I$ basis.** Matches `@eq:Tg`; it gives no charge selection rule, so component processes such as
   $\langle W^1 W^2\rangle$ produce diagrams that cancel only after summing over internal components.
3. **Sterile partners** are not separate particles: $\Psi_f=\psi_f+\chi_f$ is one Dirac field with the $d$-dimensional
   propagator, and the chirality sits in the vertices (`@eq:Lferm`). Exact, since $\chi_f$ is never external.
4. **Sources external only.** `propagating: false` is set but ignored by the generator. Internal source lines are
   excluded by `coupling_orders={'SRC': n}` with $n$ the number of external sources (every source vertex is linear in
   sources, so $n$ source vertices are all used by the external sources), plus the filter callback as a safeguard.
5. **Exactly one $\Delta_1$**: `coupling_orders={'DELTA': 1}` for ghost number 1, `0` otherwise (ghost-number
   conservation would exclude it anyway).
6. **Generation options**: `loops=1`, default `maximum_bridges` (1PI), `self_energy=tadpoles=zero_snails=None`, and
   `max_vertices = number of external legs`. The last one is exact for one-loop 1PI graphs (every vertex lies on the
   loop and carries at least one external leg; the only two-point vertex, $\rho\omega$, cannot lie on the loop) and
   removes a slowdown of up to ~500x. Verified: all 860 runs of the first 17 classes give identical diagram sets with
   and without it.
7. **Higgs mass**: parameter `MH` with $M_H^2 = m^2$ in the propagator denominator; complex $m^2$ is allowed by the
   parameter type but no numerical value is used.
8. **Coupling values** are symbolic expressions in `gs, gw (=g), gp (=g'), yu, yd, ye, lam, MH`; numerical `value` fields
   are null (unit placeholders in the parameter list).
9. Hermiticity of the Lagrangian is built in: both Yukawa terms and both matter source terms are explicit vertices.

## Component tuples and multiplet-level grouping

For every class, `classes.assignments` lists the multiplet-level external assignments with total ghost number 0 or 1,
zero total hypercharge, and an SU(3) and SU(2) singlet in the product of representations (a single adjoint index with no
other charged field is forbidden; two equal Grassmann ghosts need an antisymmetric invariant, so $\omega_B\omega_B$ and
$\delta^{ab}\omega^a\omega^b$ cases are dropped). Each assignment is expanded into component tuples: all doublet
components, and SU(2)-adjoint indices up to relabelling of the values 1,2,3 (set-partition patterns). Signed permutations
of the adjoint axes are SU(2) rotations, so every multiplet-level diagram with a non-vanishing invariant group factor
appears in at least one tuple. The multiplet-level diagram is the cyclic word of the loop (vertices labelled by their
external legs, lines by multiplet and particle-flow direction), canonical under rotation and reflection.
