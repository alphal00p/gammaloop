# Handoff: remaining GL04 Taylor-`T2` discrepancy

> **Temporary handoff material. This file is intentionally pushed under
> `HANDOFF_REMOVE_FOR_PR/`; remove that directory before the merge-ready PR.**
>
> This document records the state of the investigation at the handoff point on
> 2026-09-04. It deliberately distinguishes verified evidence from hypotheses
> and retains the exact reproduction commands. Read this file before changing
> the direct-local-3D Taylor kernel, projected-local-4D construction, EMR
> dispatcher, CFF adapter, or `three-dimensional-reps`.

## 1. Executive summary

The smallest authoritative production discrepancy currently known is the
single DOD-two UV-forest operation `{1zs}` of scalar LU graph `GL04`, with the
edge-local numerator

```text
(Q_5^0)^2
```

and with both integrated UV counterterms and threshold counterterms disabled.
The comparison is between:

1. **explicit-sum local UV from 3D**: generate the complete generalized CFF of
   the original graph, keep each complete residue-map key as one branch, and
   apply the UV Taylor operator directly to that post-energy-integration branch;
2. **projected local UV from 4D**: apply the Taylor operator to the factorized
   four-dimensional integrand, reconstruct the UV skeleton from retained
   physical-edge provenance, rewrite the coefficient in the resulting graph's
   EMR variables, generate generalized CFFs for the factorized pieces, and sum
   all their residues explicitly.

`T0` and `T1` agree at the 1000-bit Arb floor on every physical cut. `T2`
agrees on cut `[3,4]`, but disagrees exactly on cuts `[1,2]` and `[1,7,8]`.
The distinction is not that one cut has a DOD-two term and another does not:
all three have the same `{1zs}`/`T2` forest term. The distinction is that edge
3 carries the retained soft/outer momentum `p`. Cut `[3,4]` fixes `Q_3^0` as a
cut energy and produces no `Q_3^0=0` generalized outer-contact map. The other
two cuts leave `Q_3^0` unresolved and produce respectively six and four such
contact maps. At `mUV=0`, the complete discrepancy on each failing cut is
exactly `-(100/7)` times its complete projected `Q_3^0=0` contact family.

The most important remaining experiment is an **identical-input, one-call CFF
oracle for the complete actual GL04 factorized `T2` rational source**, compared
only after the complete physical residue is summed. A temporary scaffold exists
but has not produced a result yet; see section 15.

The local-3D construction is currently the **more likely oracle** because it is
conceptually simpler, applies `T` directly to the complete keyed CFF, and can be
checked for UV cancellation per residue-map key. This is a strong prior, not a
proof: it could still contain an implementation error. Conversely, a passing
UV profile does not prove equality of finite terms. If both routes subtract all
UV limits while their complete residue sums differ, the defect is confined to a
UV-finite contribution and must still be resolved by the local finite
comparison; neither route may be declared correct from UV cancellation alone.

## 2. Workspace, branch, and mandatory setup

At the time this handoff was written:

```text
workspace: /common/dev/gammaloop/higher-power-energies
branch:    raised_energy_cff_wip
base HEAD: c8e7631734f65b7fae00ea28219fb7c871b6ed1e
remote:    origin/raised_energy_cff_wip at the same commit before handoff cleanup
```

The handoff/cleanup commit made by the parent agent may advance the branch. Use
`git log` rather than assuming the base hash remains HEAD.

Use:

```bash
cd /common/dev/gammaloop/higher-power-energies
export SYMBOLICA_LICENSE='dc3cebc3#6c76981c#e3f8341d-402d-522d-87fa-6415319915a3'
export RUST_MIN_STACK=67108864
```

The user explicitly stated that the license value is not sensitive.

Run Cargo through the Nix development shell. For a focused diagnostic, prefer
an isolated target directory so a stale or concurrent shared-target build does
not confound the result:

```bash
export CARGO_TARGET_DIR=/tmp/gammaloop-t2-handoff-target
```

Before any git operation, use the identity required by the user:

```bash
git -c user.name=ValentinHirschi \
    -c user.email=valentin.hirschi@gmail.com <operation>
```

Do not run or use `tests::failing` as a gate. Follow `CONTRIBUTING.md`, including
its bird's-eye debugging protocol and local-4D reconstruction certificates.

## 3. Non-negotiable conceptual contracts

### 3.1 Local UV from 3D is the simple lane

For each complete generalized residue-map key `k`, write the full CFF as

```text
I = sigma_k * CFF_k + sigma_l * CFF_l + ...
```

where `sigma_k` is the opaque selector for the **complete residue-map key**, not
merely a coarse physical edge orientation. Apply the UV Taylor operator to each
complete branch body. The selector is opaque and remains outside the Taylor
operation. This route must perform no UV graph reconstruction.

The orientation-local and explicit-orientation-sum versions of this same route
must differ only by whether the residue-key selectors are materialized. After
summing all keys, their expressions must agree. Per-key UV finiteness is a
required property only of the orientation-local local-3D route.

### 3.2 Local UV from 4D is the reconstructed/projected lane

This lane first Taylor-expands the factorized 4D integrand. It then:

1. retains original physical edge IDs as provenance through the Taylor
   operation;
2. uses that provenance and the original graph to construct the UV skeleton,
   avoiding generic algebraic graph reconstruction from an incidence/Kirchhoff
   matrix;
3. represents denominator powers by distinct serial/degenerate edges;
4. rewrites the Taylor coefficient in the EMR variables of the corresponding
   factorized UV graph and proves equality after both expressions are written in
   the same LMB;
5. calls generalized CFF for each factorized component and combines the maps by
   a Cartesian product;
6. sums residues explicitly. This mode does **not** require per-orientation UV
   cancellation and only supports `explicit_orientation_sum_only=true`.

Future on-shell schemes may yield non-vacuum reconstructed children. The
provenance-based skeleton construction must remain applicable to them.

### 3.3 Compare complete residues, not internal pieces

The two lanes are required to agree only after summing the **complete physical
residue for one Cutkosky cut**, including every generalized ordinary/contact
piece belonging to it. Individual contact pieces, interpolation nodes, source
hosts, or synthetic LU moments are not invariant between the two constructions
and must not be compared as if they were physical terms.

All three live GL04 cuts have `max_occurrence=1`; this is an ordinary physical
residue containing generalized numerator-contact pieces, not a raised physical
LU2 residue. Earlier experiments that compared separate LU1/LU2 representatives
of a synthetic raised fixture were not production oracles.

### 3.4 Keep the numerator factorized

The numerator must stay factorized in production. Energy-degree analysis walks
the factorized product and composes degrees; it must not expand the full
numerator. No upstream algebraic cancellation of numerator and denominator is
performed. Generalized CFF realizes lower/contact sectors internally.

The original/root numerator remains attached to its original physical owner.
The EMR dispatcher may act only on **new temporal factors generated when a
Taylor derivative hits a denominator**, and only among the new serial copies of
that same physical line. It must not use global energy conservation to move an
original factor to another physical line.

## 4. Important clarification: why `T2` contains `(Q_3^0)^2`

There are two mechanisms which must not be conflated.

For the active `{1zs}` child, use

```text
q = Q5
p = Q3
Q6 = -q + p
```

The original `(Q5^0)^2` remains a fixed hard-child numerator and is never
dispatched to a new edge. When a Taylor derivative hits a denominator, every
new hard temporal factor is assigned to the new serial edge created by that
derivative; no quadratic hard load is thereby accumulated on an occupied edge.

However, the second Taylor coefficient of the shifted denominator
`D(-q + t p)` contains terms proportional to

```text
p^2 / D(q)^2
(q.p)^2 / D(q)^3
```

so the *coefficient itself* has genuine quadratic dependence on the retained
outer energy `p^0=Q3^0`. These `p^0` factors are not hard-child derivative
factors and cannot be assigned to child serial copies: `p` belongs to the outer
factorized graph. Thus the projected outer CFF correctly needs `(3,2)`. This
does not violate the dispatcher contract and does not mean `(5,2)` was moved by
the dispatcher.

## 5. Exact production MRE

The focused integration test constructs GL04 through FeynGen, attaches the two
temporal factors locally to edge 5, enables UV subtraction, and compares fresh
states for explicit-sum direct local-3D and projected local-4D.

Test filter:

```text
scalar_3l_cross_section_inspects::default_scalar_3l_cross_section_inspects::scalar_3l_cross_section_gl04_q5_temporal_square_inspects_match
```

### 5.1 Clean-branch production command

After diagnostic hooks are removed, the broader discrepancy is still exposed
by the production acceptance test itself, with all local/integrated UV and
threshold counterterms enabled and all three CFF routes built in fresh states:

```bash
export SYMBOLICA_LICENSE='dc3cebc3#6c76981c#e3f8341d-402d-522d-87fa-6415319915a3'
export RUST_MIN_STACK=67108864
export CARGO_TARGET_DIR=/tmp/gammaloop-t2-handoff-clean-target

nix develop --command cargo test \
  -p gammaloop-integration-tests \
  --test test_runs \
  scalar_3l_cross_section_inspects::default_scalar_3l_cross_section_inspects::scalar_3l_cross_section_gl04_q5_temporal_square_inspects_match \
  -- --exact --nocapture
```

This is the first command to run on the clean handoff branch. At the last
broader failing checkpoint the Arb totals were approximately

```text
orientation-local direct 3D = -9.1734755803268345365e-6 i
explicit-sum direct 3D      = -9.1734755803268345365e-6 i
projected local 4D          = -1.2003542408873861243e-5 i
```

and the test failed with a projected/direct relative difference around
`2.36e-1`. Artifact: `/tmp/gl04-coherent-wrapper-ab-20260903.log`. The exact
clean-handoff values and pass/fail status must be recorded afresh because the
parent agent is removing diagnostics and may retain later fixes not present at
that historical checkpoint.

### 5.2 Historical production-path isolation of `{1zs}`/`T2`

The following was the exact pre-cleanup production-path diagnostic command. It disables
integrated and threshold counterterms, isolates forest node 1 (`{1zs}`), sets
`mUV=0`, and selects Laurent power 0, which is `T2` for a DOD-two operation:

```bash
export SYMBOLICA_LICENSE='dc3cebc3#6c76981c#e3f8341d-402d-522d-87fa-6415319915a3'
export RUST_MIN_STACK=67108864
export CARGO_TARGET_DIR=/tmp/gammaloop-t2-handoff-target
export GL_SCALAR_DIAGNOSTIC_SKIP_PROFILE=1
export GL_SCALAR_DIAGNOSTIC_SKIP_ORIENTATION_LOCAL=1
export GL_SCALAR_DIAGNOSTIC_SUBTRACT_UV=true
export GL_SCALAR_DIAGNOSTIC_INTEGRATED_UV=false
export GL_SCALAR_DIAGNOSTIC_THRESHOLDS=false
export GL_SCALAR_DIAGNOSTIC_M_UV=0
export GL_UV_DIAGNOSTIC_FOREST_NODE=1
export GL_UV_DIAGNOSTIC_TAYLOR_LAURENT_POWER=0

nix develop --command cargo test \
  -p gammaloop-integration-tests \
  --test test_runs \
  scalar_3l_cross_section_inspects::default_scalar_3l_cross_section_inspects::scalar_3l_cross_section_gl04_q5_temporal_square_inspects_match \
  -- --exact --nocapture 2>&1 | tee /tmp/gl04-current-node1-t2-muv0.out
```

The test may print `ok` because its current final total comparison has an f64
input-floor escape hatch and different physical-cut differences partially
cancel. The authoritative evidence is the per-cut `ARB NODE DIAG` output. A
future regression must assert complete residues **per physical cut**, so this
known discrepancy cannot hide in the graph total.

The handoff cleanup intentionally removes all `GL_UV_DIAGNOSTIC_*` hooks. These
variables are therefore historical documentation, not surviving branch APIs.
Section 23 contains the exact temporary code needed to restore the isolation in
an expendable checkout. Do not restore broad debug printing to production.

### 5.3 Isolate `T0` and `T1`

Use the same command with:

```bash
# T0, coefficient t^-2 for DOD two
export GL_UV_DIAGNOSTIC_TAYLOR_LAURENT_POWER=-2

# T1, coefficient t^-1 for DOD two
export GL_UV_DIAGNOSTIC_TAYLOR_LAURENT_POWER=-1
```

Known outputs are:

```text
/tmp/gl04-current-node1-t0-muv0.out
/tmp/gl04-current-node1-t1-muv0.out
/tmp/gl04-current-node1-t2-muv0.out
```

### 5.4 What the harness actually builds

Relevant defaults in
`tests/tests/test_runs/scalar_3L_cross_section_inspects.rs` are:

```text
external Q0 = (1,0,0,0)
loop x-space sample = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
mu_r = 3
default mUV = 20 (overridden above to 0)
M = numerator_sampling_scale = 1
orientations = summed
LMB channels = summed, OSE weights
compile = false
evaluator = SingleParametric for both explicit routes
```

## 6. GL04 graph identity, LMB, and cuts

The in-memory graph has nine physical edges `0..8`; split DOT half-edges for
the initial cut are sewn before cut processing and do not survive as extra
physical edges. Exact-source occurrence IDs therefore start at 9.

The production LMB used in the audit is:

```text
e1 =  K0
e2 = -K0 + P0
e3 = -K2 - P0 + K0          = p
e4 = -K0 + K2
e5 =  K1                    = q
e6 = -K1 - K2 - P0 + K0    = -q + p
e7 =  K2
e8 = -K2 - P0 + K0         = p
```

The physical CutGroup/event-ID permutation is easy to misread:

| CutGroupId | physical cut edges | generated-event `cut_id` |
| ---: | --- | ---: |
| 0 | `[1,2]` | 2 |
| 1 | `[3,4]` | 0 |
| 2 | `[1,7,8]` | 1 |

Do not label an `ARB NODE DIAG` line by its event ID without applying this
table. An earlier investigation accidentally treated event cut 0 as CutGroup 0
and briefly claimed the wrong cut had been certified.

The exact discriminator is:

| physical cut | is edge 3 cut? | projected `Q3^0=0` contact maps | `T2` result |
| --- | --- | ---: | --- |
| `[1,2]` | no | 6 | exact failure |
| `[3,4]` | yes | 0 | Arb-floor agreement |
| `[1,7,8]` | no; edge 8 is cut | 4 | exact failure |

This is why `[3,4]` can pass while `[1,2]` fails even though both contain the
same DOD-two Taylor term.

## 7. Authoritative `mUV=0` numerical results

All values below are complete generated event weights for the isolated
`{1zs}` forest node and selected Taylor coefficient, evaluated with native
1000-bit Arb inputs. Values are purely imaginary; `P-D` means projected local
4D minus explicit-sum direct local 3D.

### 7.1 `T0`

Both routes are exactly identical in the recorded Arb values on all cuts:

| physical cut | direct = projected |
| --- | ---: |
| `[1,2]` | `-7.74568222542355281117531158015e-15 i` |
| `[3,4]` | `+2.98360199155556879528881350381e-14 i` |
| `[1,7,8]` | `+6.72657428539934271594808808301e-15 i` |

The total in both routes is
`+1.27573903229764644365034331194e-17 i`.

### 7.2 `T1`

| physical cut | direct/projected leading value | `P-D` |
| --- | ---: | ---: |
| `[1,2]` | `-6.29017661624079957615467602210e-17 i` | `+2.69394e-316 i` |
| `[3,4]` | `+2.42294777054658630229822699251e-16 i` | `-1.12399e-313 i` |
| `[1,7,8]` | `+5.46257115203456783212330878244e-17 i` | `-1.579996e-312 i` |

The complete totals agree to approximately `7.5e-319 i`. The projected outer
bound report is `(3,1)`.

### 7.3 `T2`: live discrepancy

| physical cut | direct local 3D | projected local 4D | `P-D` |
| --- | ---: | ---: | ---: |
| `[1,2]` | `+2.05373267604338749125286594161e-14 i` | `+8.14723844168360258663814015022e-14 i` | `+6.09350576564021509538527420860e-14 i` |
| `[3,4]` | `-1.96468412057704763771151284163e-13 i` | same | `+6.36599e-313 i` |
| `[1,7,8]` | `-4.52074606449374150469316409625e-14 i` | `-4.33807204985919856881055209736e-14 i` | `+1.82674014634542935882611998898e-15 i` |

The route totals are:

```text
direct    = -9.78991346621616381867126491468e-17 i
projected = -7.01141744755450283582366521063e-17 i
```

The projected outer bound report is exactly `(3,2)`.

## 8. Exact symbolic on-cut evidence

### 8.1 Passing cut `[3,4]`

Imposing the physical cut equation exactly,

```text
Q0 = E3 + E4,
```

at an exact rational/radical point makes all 12 complete outer host groups
agree separately. The full 48-branch direct sum and 12-branch projected sum are

```text
direct = projected = -5.43331072437902428657...e-8 i
relative difference = 1.43e-300
```

This certifies only `[3,4]`, which was already the passing cut.

### 8.2 Failing cut `[1,2]`

At `mUV=0` and an exact Pythagorean/radical point satisfying the `[1,2]` cut
exactly, the residual is symbolically nonzero:

```text
-i*(378559 + 157288*sqrt(29/5)) /
  (2189426688*pi^8*(-4/5+sqrt(29/5)/2)*(6/5+sqrt(29/5)/2)*
   (2+sqrt(29/5)/2)*(16/5+sqrt(29/5)/2))
```

The corresponding artifact is
`/tmp/gl04-cut0-muv0-pythag-symbolic.out` (the historical `cut0` filename refers
to CutGroupId 0, i.e. physical `[1,2]`). Exact expressions are also retained in
`/tmp/gl04-expr-direct-cut0-muv0-oncut-pythag.txt` and
`/tmp/gl04-expr-projected-cut0-muv0-oncut-pythag.txt`.

An exact on-cut test of `[1,7,8]` is likewise nonzero. Therefore the live
failure is not inherited f64 off-shellness, an approximate Cutkosky root, or
ordinary Arb rounding.

### 8.3 Exact contact relation

For `[1,2]`, the projected sum of the six `Q3^0=0` outer-contact maps was found
symbolically as

```text
C_[1,2] = i * (7/819200 + 3577*b/62914560
               + 343*b^2/4718592 + 245*b^3/9437184)
          / (b*(b+1/10)*(b+3/10)*(b+11/10)*(b+13/10))
```

in the audit's exact on-cut energy variable `b`. The complete physical-cut
residual obeys

```text
projected - direct = -(100/7) * C_[1,2].
```

For `[1,7,8]`, the four-map contact sum is

```text
C_[1,7,8] = +7 i / (327680 * denominator),
projected - direct = -5 i / (16384 * denominator)
                 = -(100/7) * C_[1,7,8].
```

The same exact ratio on both failing cuts, and the absence of such maps on the
passing cut, is the sharpest current localization of the defect. It does not by
itself prove whether the contact family is wrong in generalized CFF or is
correct but mis-composed by GammaLoop.

## 9. Projected-4D `T2` structure and certificates

The isolated coefficient uses two denominator copies from owner 5 and three
from owner 6:

```text
source owners: 5, 5, 6, 6, 6
energy IDs:    9,12,10,11,13
```

They remain five distinct topology edges and five distinct public EMR/residue
slots. Generalized CFF may correlate them internally as one repeated algebraic
energy channel for confluent residues, but that must not merge their topology
or public residue-map entries.

In common `q,p` coordinates, up to the graph's common scalar and `i` factor, the
certified coefficient has the form

```text
-q5_fixed_0^2 *
  [-D5*D6^2 + (mUV^2+p^2)*D5*D6 - 4*(p.(-q))^2*D5
   + mUV^2*D6^2 + 2*(p.(-q))*D5*D6]
  / (D5^2*D6^3)
```

with `D5=D(q)` and, in the Taylor-vacuum child, `D6=D(-q)=D(q)`.

Verified certificates:

1. The exact 4D Taylor coefficient is algebraically correct.
2. The provenance skeleton has the expected five serial occurrences.
3. The pre-dispatch and post-dispatch numerators become symbolically identical
   after both are substituted into the same child LMB.
4. The owner-built `D5^2 D6^3` denominator equals the reconstructed
   five-occurrence denominator in the common loop chart.
5. Original edge-5 factors are `TaylorFixed` and use base occurrence 9.
6. An owner-5 denominator-derived factor uses new occurrence 12.
7. Owner-6 denominator-derived factors use new occurrences 11 and 13, never
   base occurrence 10.
8. Owner-5 raw signatures are `+q`; owner-6 raw signatures are `-q`. All five
   parsed serial signatures use one canonical `+q` routing, and each odd
   numerator conversion retains `H^0 = h*r*P^0` with the hard-to-raw and
   raw-to-parsed signs. The earlier phrase `(Q,Q,-Q)` mixed distinct physical
   owners and did not show mixed routing within one serial chain.
9. The child energy assignment plan is sufficient and includes
   `[(9,2),(11,2),(12,2),(13,2)]`.
10. The outer source requests and receives `(3,2)`.
11. Every projected `Q(i,...)` factor is attached to the CFF component whose
    graph contains physical edge `i`; no wrong-component owner was found.
12. Child and outer residue maps are combined as a Cartesian product, and the
    source-energy map survives the final mapped multiplication once.

The structural regression is:

```text
graph::three_d_source::gl04_t2_planned_lift_matches_post_t_numerator_in_common_loop_coordinates
```

## 10. Direct-3D branch certificate and caveat

The direct root generalized map for `(Q5^0)^2` samples the original numerator at

```text
Q5^0 = -E5, 0, +E5.
```

For `N(q0)=q0^2`, the contact finite difference is exactly

```text
[2*N(-E5) - 4*N(0) + 2*N(+E5)] / (2*E5)^2 = 1.
```

The mapped branch numerator is attached before the OSE rewrite, rescaling, and
Taylor series. The complete residue key is retained outside as an opaque
selector. Tracing found no missing quadratic root key, misplaced mapped
numerator, or wrong finite-difference ordering.

One caveat must be remembered when reading the expression trace: physical tree
denominators are multiplied into the root atom under the opaque
`tree_denoms(...)` wrapper. They are *syntactically* inside the atom passed to
the Taylor kernel, but Symbolica treats that wrapper as independent of `t`, so
they are operationally frozen and factorized outside the Taylor action. Do not
claim that the Taylor series acts on those denominators merely from their
printed syntactic position. This is intended to retain untouched outer/tree
factors, but the identical-input one-call oracle must respect and test the same
factorization boundary.

Relevant code is in:

```text
crates/gammalooprs/src/uv/approx/direct_3d/forest.rs
crates/gammalooprs/src/uv/approx/direct_3d/branches.rs
crates/gammalooprs/src/uv/approx/direct_3d/kernel.rs
```

## 11. Rank-envelope experiments

The default physical/root numerator analysis reports `(5,2)` for the original
`(Q5^0)^2`. The projected `T2` outer coefficient correctly reports `(3,2)`.
Because Taylor expansion exposes a quadratic retained outer energy, a serious
hypothesis was that the global root CFF had not been provisioned with enough
capacity to remain closed under the Taylor action.

This was tested in an isolated checkout by forcing supported root bounds:

1. `[(3,2),(5,2)]`;
2. `[(3,2),(5,2),(8,2)]`, covering the original hard owner and both physical
   carriers of the retained soft routing.

Both experiments leave every per-cut Arb value and route delta unchanged to all
printed digits. Compare:

```text
baseline:
  /tmp/gl04-current-node1-t2-muv0.out
forced (3,2),(5,2):
  /tmp/gl04-root-capacity-node1-diagnostic.out
forced (3,2),(5,2),(8,2):
  /tmp/gl04-root-capacity-e3-e5-e8-diagnostic.out
```

The user also requested an intentionally excessive `(edge,2)` bound for every
internal edge. That experiment reached no evaluation because the current Rust
generalized-CFF port rejects this multi-edge higher-rank sector:

```text
Some([(1,2),(2,2),(3,2),(4,2),(5,2),(6,2),(7,2),(8,2)])
-> this generalized CFF higher energy-numerator sector is not supported by the
   current Rust port
```

Artifact: `/tmp/gl04-root-capacity-all-edge2-diagnostic.out`.

The supported, evidence-directed overprovisioning is already enough to rule out
missing `(3,2)` root capacity as the cause. Do not generalize the unsupported
all-edge experiment into a production fix.

Separately, exact-CFF cache capacity pooling over certified copies of one
repeated algebraic channel was tested with a componentwise-max envelope and with
the cache bypassed entirely. Both were inert. Rank pooling is therefore not the
known missing contact.

## 12. Scale, precision, and parameter experiments

### 12.1 Arb precision

For the stripped temporal-square comparison at default `mUV=20`, the total
route difference is stable for roughly 70 decimal digits:

```text
256-bit  +2.262532013744616378263364907633471491633795616772479646570887905061893401622075e-12 i
1000-bit +2.262532013744616378263364907633471491633795616772479646570887905061893572314082...e-12 i
```

At 1000 bits:

```text
direct    = -9.1734755805810912839071256602802601e-6 i
projected = -9.1734733180490775392907473969153524e-6 i
relative difference = 2.4663847348480072054e-7
```

This is structural, not Arb noise. Quad and ordinary precision show the same
leading behavior.

### 12.2 Cut-surface precision

Production native-Arb Cutkosky solving satisfies the three physical surface
equations at approximately `1e-301`:

```text
[1,2]:   2*|Q1| - Q0 = O(1e-301)
[3,4]:   2*|Q1-Q7| - Q0 = O(1e-301)
[1,7,8]: |Q1| + |Q7| + sqrt(|Q1-Q7|^2+0.1^2) - Q0 = O(1e-301)
```

The f64 records have `O(1e-17)` residuals, but the Arb path recomputes `t*` and
the rescaled sample natively rather than merely upcasting that defect.

A separate precision-promotion bug was fixed: static
`ParamBuilder<f64> -> ArbPrec` values now convert directly rather than passing
through Quad. That regression is useful but did not change the GL04 mismatch.

### 12.3 Parameter vectors and cache

The complete ordered `ParamBuilder` vectors and replacement maps are
byte-identical between direct and projected routes. Audited slots include:

```text
t*, all external and loop components, mUV, mUVexp, mu_r^2, M,
all residue/orientation slots, and all model parameters.
```

Disabling the exact-source cache does not change the discrepancy. This rules
out a stale route cache as the cause.

### 12.4 `mUV`

`mUV=0` is valid and makes the failure much larger, as section 7 records. Thus
the discrepancy cannot be solely a missing explicit `mUV^2` term in `T2`.

At fixed default kinematics, one smaller failing-cut delta was observed as:

```text
mUV=5:  approximately 3.6468503e-20 i
mUV=20: approximately 5.7041057e-22 i
mUV=80: approximately 8.9132420e-24 i
```

It scales as `mUV^-3`, while the complete `T2` contribution scales as
`mUV^-1`; the relative discrepancy scales as `mUV^-2`. Scaling the whole
nine-coordinate sample by `0.4` amplifies the delta to about `1.3001942e-10 i`
while retaining a relative discrepancy around `2.5e-5`.

### 12.5 Numerator sampling scale `M`

Running the actual production `ParamBuilder` with `M=0.5,1,2` leaves direct,
projected, and all per-cut deltas bit-for-bit unchanged. The auxiliary CFF
sampling scale cancels as intended. No `M`/`mUV` confusion was found.

### 12.6 `mu_r`

The standalone expressions use `mUV`/`mUVexp` for the Taylor deformation;
`mu_r^2` remains a separate runtime parameter. No substitution of `mu_r` for
`mUV`, or vice versa, was found in this MRE.

## 13. Verified UV scaling law

The 4D implementation uses the inverse-scaling convention `Q=H/t+S`. After the
common `1/t^2` is exposed, `Graph::uv_rescaled` produces

```text
[(H+t*S)^2 - t^2*(m^2-mUVexp^2) - mUVvac^2] / t^2.
```

When `mUVexp=mUVvac=mUV`, this is precisely the user's required deformation,
up to the inverse variable convention:

```text
1 / [(k+lambda*p)^2 - lambda^2*(m^2-mUV^2) - mUV^2].
```

The direct-local-3D OSE rewrite is the positive-energy square root of the same
polynomial with the correct 3D measure scaling. `ParamBuilder` maps both
`GS.m_uv_vacuum` and `GS.m_uv_expansion` to the same runtime `general.m_uv` in
the current MSbar/MUV setup.

Relevant code:

```text
crates/gammalooprs/src/uv/approx/local_4d.rs        Graph::uv_rescaled
crates/gammalooprs/src/uv/approx/direct_3d/kernel.rs  OSE rescaling/Taylor series
crates/gammalooprs/src/integrands/process/param_builder.rs
```

Given exact `T0`/`T1`, exact common-LMB `T2` coefficient certificates, and a
failure that persists at `mUV=0`, the Taylor mass deformation is now a low-rank
hypothesis. Reopen it only if a new symbolic counterexample is produced.

## 14. Other experiments and what they proved

### 14.1 Spatial-square control

Replacing `(Q5^0)^2` by a localized spatial product `Q5^1 Q5^1` keeps the UV
degree but removes generalized temporal contact interpolation. At one sample:

```text
direct local 3D  = -6.6390634494155886293366743389748287e-10 i
projected local4D= -6.6390534378633180985233346720561283e-10 i
P-D              = +1.0011552270530813339666918700e-15 i
```

This ruled out a wholesale projected-`T2` sign reversal, but it is not a proof
that the spatial case is correct: a stripped projected UV profile for this
control failed, with an initially expected DOD-two limit fitting about `+2.000`
and a subsequently expected DOD-zero limit fitting about zero rather than a
negative subtracted power. Treat the tiny pointwise delta as accidental at that
sample.

### 14.2 Powered-denominator identities

The independent identity

```text
D(Q) / [D(Q) D(-Q)] = 1 / D(-Q)
```

now passes in the raw generalized engine and GammaLoop. It found a real adapter
bug: a generalized child was receiving an ordinary CFF core sign a second time
when retained positive local `1/(2E)` factors were converted to GammaLoop's
historical `1/prod(-2E)` convention. The corrected conversion is typed by
energy-factor ownership and is uniform over the complete Laurent functional.

This fix is principled and covered by unit-numerator/rank-capacity controls. It
cannot turn `F+C` into `F-C` only for one GL04 contact subfamily, so it is not a
fit for the remaining discrepancy.

Relevant tests include:

```text
cff::tests::exact_cff_uncancelled_powered_denominator_matches_lower_source
cff::tests::exact_cff_powered_rational_identities_match_at_arb_precision
cff::tests::production_rank_capacity_only_uses_typed_energy_factor_conversion
cff::tests::two_loop_rank_capacity_rejects_a_second_generalized_core_sign
three_dimensional_reps::eval::unused_rank_capacity_preserves_raw_scalar_cff
```

### 14.3 Historical child-only `T2` oracle (removed during cleanup)

The actual DOD-two temporal self-energy child, stripped of the outer graph,
agrees between direct Taylor-on-CFF and exact local-4D CFF at

```text
p0 = 0, +1, -1, +sqrt(13)/10, -sqrt(13)/10.
```

It is even under `p0 -> -p0`. This certifies the complete child function at
several outer-energy values, not only one interpolation node.

The focused ignored diagnostic test was:

```text
uv::approx::tests::quadratic_temporal_self_energy_t2_matches_exact_local_4d_cff
```

It was intentionally removed with the diagnostic cleanup and is **not present
in the clean handoff tree**. Its result above remains evidence, not a surviving
regression gate. If this exact child oracle is needed again, recover its source
from the pre-cleanup patch/history or reconstruct it in an isolated checkout;
do not expect `cargo test -- --list` to show it. Its historical invocation used
`GL_UV_DIAGNOSTIC_TAYLOR_LAURENT_POWER=0` and `--ignored`.

### 14.4 Component composition and hosting

Audits found:

- every pre-outer active projected coefficient is cloned identically over its
  production hosts;
- every authoritative non-`None` source-energy map from the outer CFF survives
  `zip_mul_mapped_factor` exactly once;
- bypassing `localized_orientation_terms` and rehosting all raw reduced-outer
  branches under one bookkeeping ID leaves the discrepancy unchanged;
- `CffGenerationContext::Standalone` versus `EmbeddedCffFactor` is inert here,
  because this nonterminal generalized source invokes the same bounded builder;
- the projected production MRE contains no deferred `projected_cff_sum` call;
  it is already materialized into ordinary integrands.

These observations weaken simple host-ID/filtering explanations, but they do
not yet replace the decisive one-call complete-source comparison.

### 14.5 Generalized contact identities

The complete exact child satisfies

```text
(q0)^2 = D(q) + |q_vec|^2 + m^2
```

with a **positive** lower-topology contact. The opposite sign fails. The raw
generalized engine also has extensive repeated-pole, routing-reversal, lower
sector, disconnected-component, and high-rank contact tests. Therefore do not
change `three-dimensional-reps` merely because GammaLoop's two production
routes disagree. Change it only if the identical-input raw one-call oracle
fails against a trustworthy lower-topology identity.

### 14.6 `3Drep` CLI status

Earlier simple temporal/powered examples were reproduced through the diagnostic
`3Drep` CLI and supported the generalized CFF residue maps. The CLI is
diagnostic only: its input preparation and expression display are not a
GammaLoop production contract. It can help test an identical rational source,
but a CLI-only difference must not drive a production-specific patch.

### 14.7 Cache and stale-state controls

Every UV route is built in a fresh graph/evaluator state in the scalar harness.
Disabling exact-source cache reuse leaves the GL04 discrepancy unchanged.
State contamination and cache-key reuse are therefore exonerated for this MRE.

## 15. One-call complete-source oracle: current status

This is the highest-priority next step. Existing green tests certify adjacent
boundaries but not their complete composition:

- `gl04_t2_planned_lift_matches_post_t_numerator_in_common_loop_coordinates`
  proves the exact coefficient, skeleton, occurrence plan, and common-LMB EMR
  rewrite, but stops before CFF generation;
- `exact_raised_lu_residue_factorizes_from_quadratic_cubic_spectator` and the
  child-only test prove generalized CFF factorizations for related inputs, but
  not the full actual GL04 coefficient embedded in its actual outer graph.

An optional temporary, uncommitted scaffold exists at:

```text
/tmp/gl04-onecall.KO3L9q
```

This directory is ephemeral and nonportable. Section 24 embeds the
authoritative source; use that embedded copy if the `/tmp` checkout is absent
or disagrees with this document.

The test is appended to:

```text
/tmp/gl04-onecall.KO3L9q/crates/gammalooprs/src/uv/approx/mod.rs
```

Resume it with:

```bash
export SYMBOLICA_LICENSE='dc3cebc3#6c76981c#e3f8341d-402d-522d-87fa-6415319915a3'
export RUST_MIN_STACK=67108864
export CARGO_TARGET_DIR=/tmp/gl04-onecall-target
export GL_UV_DIAGNOSTIC_TAYLOR_LAURENT_POWER=0

nix develop --command cargo test \
  --manifest-path /tmp/gl04-onecall.KO3L9q/Cargo.toml \
  -p gammalooprs \
  --profile dev-optim \
  gl04_complete_t2_one_call_oracle \
  -- --ignored --nocapture
```

There is **no authoritative result yet**. An earlier prototype in a different
temporary checkout first failed to compile, then a broader GL04 run failed
during source extraction with:

```text
source edge 6 has no unique routing for rewritten exact signature [0,0,-1]
modulo the factorized opposite-domain loop span
```

That was an oracle-scaffold/reconstruction problem, not evidence that either
production route was wrong. Do not record the one-call test as passing or
failing until it accepts the exact same complete physical source and sums the
complete physical residue.

The decisive interpretation is:

1. **one-call = projected, not direct:** investigate direct post-CFF Taylor
   action or its consumption of the complete Q5 root key;
2. **one-call = direct, not projected:** investigate sequential child-times-
   outer CFF composition and outer contact attachment in projected local 4D;
3. **one-call differs from both:** the scaffold is not identical-input yet, or
   a shared GammaLoop CFF adapter/source-map boundary is wrong;
4. **raw generalized one-call violates a lower-topology identity:** only then
   investigate `three-dimensional-reps` itself.

## 16. Ranked remaining hypotheses

The ranking below incorporates a separate clean-room review of all evidence.

### 1. Projected sequential child-times-outer composition/contact attachment

This is currently the strongest live hypothesis. The exact discrepancy occurs
only where an unresolved quadratic outer energy creates `Q3^0=0` contact maps,
and it is exactly proportional to their complete family on both failing cuts.
The child alone is certified. A subtle mismatch can therefore first appear when
that child coefficient is interpreted as the numerator of the reduced outer
CFF and its complete contact family is composed with the child map.

**Decisive test:** finish the one-call complete-source oracle, then compare its
complete physical-cut result to the sequential projected construction.

### 2. Direct Q5-key Taylor/CFF commutation or mapped-key consumption

Mathematically the direct route is the more likely oracle in design: apply `T`
to each complete post-energy-integration residue-key branch, and verify local
subtraction per key. That simplicity is evidence, not an assumption that it
cannot be wrong. If the sequential one-call
agrees with projected, the direct branch may be losing/misweighting a contact
when Taylor acts on the original quadratic Q5 functional. However, exact traces
already exonerate the finite-difference formula, mapped-numerator ordering, and
supported root rank capacity, so this hypothesis is narrower than it first
appeared.

**Test:** compare one complete failing-cut direct branch sum against the same
identical one-call rational source, keeping the opaque `tree_denoms` caveat in
view.

### 3. Edge-3/edge-8 retained-carrier or source-map alias mismatch

`e3` and `e8` carry the same `p` routing in the production LMB, while the two
failing cuts resolve them differently. The physical common-LMB rewrite is
correct, but a later source-energy map could still distinguish aliases when it
should not, especially in a zero-sample contact.

**Tests:** exchange the legal retained carrier (`e3` versus `e8`) without
changing the common-LMB polynomial; compare complete cut sums and emitted outer
maps. Require invariance, not equality of individual internal pieces.

### 4. Lower-sector normalization/parity at the GammaLoop adapter boundary

The previous typed normalization bug was real and fixed, and broad identities
now pass. A more local component-product parity mistake remains possible, but a
new change is justified only if an identical-input ordinary/generalized control
isolates it. The exact `-(100/7)` contact relation is a clue, not permission to
insert another sign or multiplicity.

**Test:** extract the complete lower/contact sector from the identical one-call
and sequential constructions before GammaLoop conversion, then after typed
conversion. The first boundary where their complete sums diverge owns the bug.

### 5. Direct rescaling or 4D `T2` algebra

Low probability. `T0`/`T1`, the scaling law, the exact 4D coefficient, and the
common-LMB numerator certificate all pass; `mUV=0` removes the explicit mass
piece without curing the discrepancy.

### 6. Core `three-dimensional-reps`

Lowest priority until the identical-input one-call fails a trusted algebraic
identity. The crate has extensive independent generalized-CFF validation, and
the child-only, powered-denominator, routing-reversal, and rank-capacity tests
pass. Do not modify core contour/contact signs to fit GL04.

### UV-profile adjudication caveat

Run both profile modes because a failure immediately identifies a wrong
subtraction:

- direct local 3D must be tested per complete residue-map key;
- projected local 4D is tested only after its complete residue sum.

If only one fails, that route is wrong. If both pass, the route discrepancy may
be a finite term invisible to the asymptotic profile. In that case the complete
residue-summed pointwise comparison and identical-source oracle remain the
adjudicators; do not use “both UV-subtract” as permission to ignore the mismatch.

## 17. Ruled-out or strongly exonerated causes

Do not reopen these without new evidence at an earlier boundary:

| Candidate cause | Evidence against it |
| --- | --- |
| integrated UV localization | MRE disables integrated UV completely |
| threshold CT or raised threshold | MRE disables thresholds completely |
| physical raised-Cutkosky derivative | all live physical cuts have `max_occurrence=1` |
| imprecise f64-to-Arb cut solution | native Arb surface residuals are `O(1e-301)`; exact radical cuts still fail where expected |
| mismatched runtime inputs | complete ordered parameter vectors and replacement maps are identical |
| `mUV` confused with `mu_r` or `M` | symbols/slots are distinct; `M` scan inert; exact `mUV=0` failure remains |
| wrong UV mass deformation | both lanes implement the same certified polynomial; T0/T1 pass |
| missing projected outer `(3,2)` | bound report is exactly `(3,2)` |
| missing direct root `(3,2)` capacity | supported forced root bounds are digit-identical to baseline |
| generic derivative EMR dispatcher | original numerator stays fixed; new hard factors use created copies; common-LMB symbolic certificate passes |
| two derivatives hard-coded | occurrence vectors, Laurent order, and assignments are dynamic; add a T3 regression after resolution |
| mixed same-owner `Q,Q,-Q` serial topology | same-owner chains use one parsed routing; raw minus belongs to owner 6 and sign is retained |
| collapsing serial edges | five distinct topology/EMR/residue slots are retained |
| upstream numerator/denominator cancellation | none is performed; generalized CFF owns contacts |
| exact-CFF cache | disabled-cache and overprovisioned-envelope results are unchanged |
| selector-host choice/filtering | bypassing localized host placement is inert |
| Standalone vs Embedded context | same bounded builder here; forced context is inert |
| global overall `T2` sign | spatial control has same leading sign; T0/T1 and `[3,4]` pass |
| changing `M` sampling nodes | `M=0.5,1,2` is bitwise invariant |
| comparing individual residue pieces | explicitly invalid; authoritative result sums complete physical residue |

## 18. Relevant production code paths

### Production CFF and bounds

```text
crates/gammalooprs/src/processes/cross_section.rs
  CrossSectionGraph::generate_cff / root CFF creation

crates/gammalooprs/src/cff/generation.rs
  Graph::production_cff_3d_expression_options
  exact-source generation/cache/capacity logic

crates/gammalooprs/src/cff/mod.rs
  exact 4D-denominator CFF generation
  typed GammaLoop/generalized-CFF prefactor conversion
```

### Provenance skeleton, serial occurrences, and EMR maps

```text
crates/gammalooprs/src/graph/three_d_source.rs
  physical_owner_skeleton_from_denominators
  serial_power_chain_incidences
  ExactSourceEnergyMapper
  common-LMB and owner/routing certificates

crates/gammalooprs/src/numerator/energy_degree.rs
  factorized degree analysis
  original-versus-DenominatorDerived occurrence eligibility
  deterministic minimax assignment
```

### The two UV lanes

```text
crates/gammalooprs/src/uv/approx/direct_3d/
  branches.rs
  forest.rs
  kernel.rs

crates/gammalooprs/src/uv/approx/projected_4d.rs
crates/gammalooprs/src/uv/approx/local_4d.rs
crates/gammalooprs/src/uv/approx/final_integrand.rs
crates/gammalooprs/src/uv/approx/local_3d/residue_localizer.rs
crates/gammalooprs/src/uv/hedge_poset.rs
```

### Evaluator inputs and precision

```text
crates/gammalooprs/src/integrands/process/param_builder.rs
tests/tests/test_runs/scalar_3L_cross_section_inspects.rs
```

### Generalized CFF engine

```text
crates/three-dimensional-reps/src/generation.rs
crates/three-dimensional-reps/src/eval.rs
```

## 19. Focused tests and their cleanup status

Except where explicitly marked historical, the following names were present in
the clean tree when this handoff was written. Check the actual handoff commit
with `rg` or `cargo test -- --list` before relying on a name; cleanup was still
ongoing in the parent task.

### Structural/provenance/dispatcher

```text
gl04_t2_planned_lift_matches_post_t_numerator_in_common_loop_coordinates
denominator_derived_dispatch_composes_each_occurrence_routing_sign
exact_energy_mapper_separates_fixed_owner_routing_from_derived_dispatch
exact_energy_bounds_use_the_explicit_base_independent_of_occurrence_order
source_energy_candidates_fix_originals_and_offer_all_derived_occurrences
partitioned_source_candidates_do_not_infer_the_base_from_energy_id_order
partitioned_source_candidates_reject_derived_factors_without_a_copy
positive_denominator_completion_obeys_embedded_emr_provenance
positive_denominator_completion_rejects_mixed_emr_provenance
```

### Generalized CFF and normalization

```text
exact_cff_uncancelled_powered_denominator_matches_lower_source
exact_cff_uncancelled_powered_denominator_matches_lower_lu_residues
exact_raised_lu_residue_factorizes_from_quadratic_cubic_spectator
exact_cff_powered_rational_identities_match_at_arb_precision
production_rank_capacity_only_uses_typed_energy_factor_conversion
two_loop_rank_capacity_rejects_a_second_generalized_core_sign
nonterminal_quadratic_contact_matches_independent_lower_cff_frame
factorized_lower_sector_restores_global_core_sign_and_orientation_sum
independent_bubbles_bounded_maps_match_component_cff_product
two_loop_mixed_carrier_denominator_numerator_pinches_locally
reversed_vacuum_triangle_incidence_preserves_generalized_cff
```

### Local UV composition

```text
complete_self_energy_taylor_sum_matches_direct_3d_for_raised_lu_jets
factorized_owned_dot_child_cff_matches_direct_3d_for_uncut_self_energy
nested_scalar_banana_direct_3d_matches_typed_local_4d_with_arb
```

The following was a temporary ignored diagnostic and was removed during
cleanup; it is not a clean-tree test:

```text
quadratic_temporal_self_energy_t2_matches_exact_local_4d_cff
```

### End-to-end scalar MRE

```text
scalar_3l_cross_section_gl04_q5_temporal_square_inspects_match
scalar_3l_cross_section_gl04_q5_spatial_square_inspects_match
scalar_3l_cross_section_gl04_q5_temporal_component_inspects_match
scalar_3l_cross_section_gl04_q6_temporal_component_inspects_match
scalar_3l_cross_section_gl04_inspects_match
```

The complete scalar matrix, UV profiles, DDx/TTx acceptances, and
`just test_gammaloop` are **not** green gates at this handoff until the parent
agent's cleanup/test report says otherwise. Do not infer a matrix pass from a
focused diagnostic that prints `ok` while retaining nonzero per-cut deltas.

## 20. Key artifacts

`/tmp` is ephemeral. Preserve/copy an artifact before relying on it across
machine restarts.

### Current production isolation

```text
/tmp/gl04-current-node1-t0-muv0.out
/tmp/gl04-current-node1-t1-muv0.out
/tmp/gl04-current-node1-t2-muv0.out
/tmp/gl04-muv-scan.out
/tmp/gl04-param-vector.log
```

### Exact complete expressions and on-cut checks

```text
/tmp/gl04-direct-g0.expr       # CutGroup 0 = [1,2]
/tmp/gl04-projected-g0.expr
/tmp/gl04-direct-g1.expr       # CutGroup 1 = [3,4]
/tmp/gl04-projected-g1.expr
/tmp/gl04-direct-g2.expr       # CutGroup 2 = [1,7,8]
/tmp/gl04-projected-g2.expr
/tmp/gl04-expr-direct-cut0-muv0-oncut-pythag.txt
/tmp/gl04-expr-projected-cut0-muv0-oncut-pythag.txt
/tmp/gl04-cut0-muv0-pythag-symbolic.out
```

### Rank-capacity controls

```text
/tmp/gl04-root-capacity-node1-diagnostic.out
/tmp/gl04-root-capacity-e3-e5-e8-diagnostic.out
/tmp/gl04-root-capacity-all-edge2-diagnostic.out
```

### Standalone production exports

```text
/tmp/gl04-standalone.cPQrzj/direct/processes/cross_sections/
  scalar_3l_gl04_q5_temporal_square_gl04_component_0_pair_numerator/
  numerator/standalone_cross_section.json

/tmp/gl04-standalone.cPQrzj/projected/processes/cross_sections/
  scalar_3l_gl04_q5_temporal_square_gl04_component_0_pair_numerator/
  numerator/standalone_cross_section.json
```

### One-call work

```text
/tmp/gl04-onecall.KO3L9q
/tmp/gl04-onecall-target
/tmp/gl04-complete-oneshot3.out
/tmp/gl04-complete-oneshot2.GkhJhk/run.out
```

### Historical detailed traces

There are many additional `/tmp/gl04-*` traces. The most useful families are:

```text
/tmp/gl04-family-audit.PYlVad/
/tmp/gl04-group0-families-regenerated.out
/tmp/gl04-group0-muv-1000.out
/tmp/gl04-outer-map-audit.KXDnko/
/tmp/gl04-pv-audit.hJYTb0/
/tmp/gl04-rank-audit.g8jMN9/
/tmp/gl04-localuv-exports.sFHbuC/
```

Treat historical filenames containing `cut0`, `cut1`, or `cut2` carefully:
some mean CutGroupId and some mean generated-event `cut_id`. Always inspect the
logged physical `CutSet` and apply section 6's permutation.

## 21. Dirty-tree cautions

At the start of handoff cleanup the worktree contained intended modifications
across the CFF, exact-source, numerator-degree, both UV lanes, ParamBuilder,
UV-profile, generalized-3D-rep, architecture docs, and scalar/DDx acceptance
tests. It also contained generated state files, target directories, `.snap.new`
files, and `.jj/` metadata.

Do not use destructive reset/checkout commands. Do not delete untracked files
indiscriminately. Compare the actual handoff commit to `c8e763173` to learn what
the parent retained. In particular:

- the entire `HANDOFF_REMOVE_FOR_PR/` directory is intentionally committed and
  pushed for handoff traceability, but must be removed before the merge-ready
  pull request;
- generated root files such as `model.json`, `global_settings.toml`,
  `default_runtime_settings.toml`, `state_manifest.toml`, `symbolica_state.bin`,
  `processes/`, and local target directories are not source;
- do not accept or delete `.snap.new` files without inspecting whether the
  corresponding existing expectation was explicitly approved by the user;
- if changing an existing test expectation, repository policy requires explicit
  user approval. Several specific expectation updates were already approved in
  the long thread, but do not infer blanket approval for new ones.

## 22. Recommended next sequence

1. Read `CONTRIBUTING.md`, this file, and
   `HANDOFF_REMOVE_FOR_PR/HANDOFF_GOAL.md`.
2. Re-run the exact `mUV=0` T0/T1/T2 MRE and confirm the cut-ID mapping.
3. Finish the temporary one-call exact-source oracle before changing signs,
   multiplicities, or generalized-CFF core logic.
4. Locate the first boundary at which the **complete failing-cut residue sum**
   diverges: raw exact source, generalized generation, GammaLoop conversion,
   sequential child/outer composition, or direct Taylor consumption.
5. Fix that boundary generically. Do not insert a GL04-, cut-, `T2`-, or
   `100/7`-specific patch.
6. Add a focused regression that asserts complete per-cut equality and cannot
   hide behind graph-total cancellation.
7. Run focused structural/unit tests after every meaningful change.
8. Then restore integrated UV and threshold counterterms and run the complete
   scalar LU three-route matrix plus both required UV-profile modes.
9. Only after an all-green scalar line should the broader goal proceed to
   refactoring, DDx/TTx, full `just test_gammaloop`, and merge-readiness.

The key discipline is to keep the two physical constructions simple and
separate: direct local 3D is `T` acting on each complete keyed CFF branch;
projected local 4D is a certified 4D Taylor coefficient rewritten into the
provenance-derived factorized UV graph and then residue-summed. The remaining
bug must be found at the first concrete implementation boundary where those two
complete sums cease to agree, not by adding another inferred sign bridge.

## 23. Recreating the temporary forest/Taylor isolation after cleanup

The handoff commit deliberately removes every `GL_UV_DIAGNOSTIC_*` production
hook. Those variables were useful while reducing the failure, but retaining
environment-controlled semantic branches in production would be the wrong
design. Consequently, commands earlier in this document that set
`GL_UV_DIAGNOSTIC_FOREST_NODE` or
`GL_UV_DIAGNOSTIC_TAYLOR_LAURENT_POWER` describe historical runs and will not
isolate anything in the clean handoff branch.

The clean branch still exposes the broad discrepancy through the normal test,
without any diagnostic environment variable:

```bash
export SYMBOLICA_LICENSE='dc3cebc3#6c76981c#e3f8341d-402d-522d-87fa-6415319915a3'
export RUST_MIN_STACK=67108864
export CARGO_TARGET_DIR=/tmp/gammaloop-handoff-target
nix develop --command cargo test \
  -p gammaloop-integration-tests \
  --test test_runs \
  --profile dev-optim \
  scalar_3l_cross_section_inspects::default_scalar_3l_cross_section_inspects::scalar_3l_cross_section_gl04_q5_temporal_square_inspects_match \
  -- --exact --nocapture
```

At the last full-configuration checkpoint captured in
`/tmp/gl04-coherent-wrapper-ab-20260903.log`, with integrated and threshold CTs
enabled, the direct orientation-local and direct explicit-sum lanes agreed at
`-9.1734755803268345365e-6 i`, while projected local 4D gave
`-1.2003542408873861243e-5 i` (relative discrepancy about `0.2358`). Reconfirm
the numbers from the clean handoff commit rather than treating that historical
snapshot as an acceptance target.

For exact forest-node and Taylor-coefficient isolation, make a throw-away
checkout or worktree and apply the following *temporary* edits. Do not apply
them in the shared dirty checkout, do not commit them, and revert the temporary
worktree after extracting the trace.

One reproducible way to start from the handoff commit is:

```bash
git -c user.name=ValentinHirschi \
    -c user.email=valentin.hirschi@gmail.com \
    worktree add --detach /tmp/gammaloop-t2-isolation HEAD
cd /tmp/gammaloop-t2-isolation
```

If that explicit target already exists, choose another narrow `/tmp` target;
do not remove or overwrite an existing worktree merely to reuse the name.

### 23.1 Select forest node 1

In
`crates/gammalooprs/src/uv/hedge_poset.rs`, inside
`Wood::orientation_parametric_exprs`, temporarily enumerate the compatible
topological order and keep only index 1:

```diff
-            for nidx in self.compatible_topological_order(compatible_subset)? {
+            for (node_index, nidx) in self
+                .compatible_topological_order(compatible_subset)?
+                .into_iter()
+                .enumerate()
+            {
+                // HANDOFF DIAGNOSTIC ONLY: GL04 child bubble R={e5,e6}.
+                if node_index != 1 {
+                    continue;
+                }
                 let operation = &self.graph[nidx];
```

The integer is a diagnostic ordering index, not a stable physical identifier.
Verify in the trace that index 1 is still the child bubble `R={e5,e6}` before
interpreting results. A production regression must select a physical forest
operation, not depend on this index.

### 23.2 Select one Laurent coefficient in both UV lanes

The GL04 node above has DOD 2. In the series conventions used by the isolated
run, Laurent powers `-2`, `-1`, and `0` isolate respectively `T0`, `T1`,
and `T2`. To isolate one coefficient, import
`symbolica::domains::rational::Rational` where needed and temporarily replace
the full truncated series in both lanes.

In `crates/gammalooprs/src/uv/approx/direct_3d/kernel.rs`:

```diff
         let series = atomarg.series(GS.rescale, Atom::Zero, 0).unwrap();
-        let series_atom = series.to_atom();
+        // HANDOFF DIAGNOSTIC ONLY. Use -2 for T0, -1 for T1, 0 for T2.
+        let series_atom = series.coefficient(Rational::from(0));
```

In `crates/gammalooprs/src/uv/approx/local_4d.rs`:

```diff
-    let series = rescaled
-        .series(GS.rescale, Atom::Zero, 0)
-        .unwrap()
-        .to_atom();
+    let raw_series = rescaled.series(GS.rescale, Atom::Zero, 0).unwrap();
+    // HANDOFF DIAGNOSTIC ONLY. Use -2 for T0, -1 for T1, 0 for T2.
+    let series = raw_series.coefficient(Rational::from(0));
```

Do not generalize the `-2/-1/0` mapping blindly to another DOD or forest node.
It is the observed convention for this exact DOD-2 MRE. The production
implementation correctly consumes the complete Taylor truncation.

When a selected coefficient vanishes for one cut, preserve the typed zero in
the test trace rather than dropping the cut and accidentally shifting the
event-index mapping. If the temporary harness filters zero terms, adjust that
harness only so it emits one value for every physical cut; do not add such a
special case to production evaluation.

### 23.3 Use the reduced runtime configuration

In the throw-away diagnostic version of
`tests/tests/test_runs/scalar_3L_cross_section_inspects.rs`, make the GL04
temporal-square test use:

- `generate_integrated = false`;
- threshold generation and runtime threshold subtraction disabled;
- `m_uv = 0`, `mu_r^2 = 9`, and numerator sampling scale `M = 1`;
- Arb evaluation of every cut, printed together with its physical `CutSet`;
- explicit orientation sum and explicit LMB sum/multichanneling, not Monte
  Carlo orientation sampling.

The relevant temporary replacements in
`setup_scalar_3l_cross_section_cli` are mechanical and should read:

```diff
-global.generation.uv.generate_integrated=true
+global.generation.uv.generate_integrated=false
-global.generation.threshold_subtraction.enable_thresholds=true
+global.generation.threshold_subtraction.enable_thresholds=false
```

and in the embedded runtime settings:

```diff
-m_uv = 20.0
+m_uv = 0.0
 numerator_sampling_scale = 1.0
@@
-disable_threshold_subtraction = false
+disable_threshold_subtraction = true
```

Temporarily set `profile_scalar_uv_routes = false`; UV profiling is a later
production gate, not part of coefficient isolation. To retain exact per-cut
evidence, change only the diagnostic Arb closure from minimal/no-events to:

```diff
-                minimal_output: true,
-                return_generated_events: Some(false),
+                minimal_output: false,
+                return_generated_events: Some(true),
```

and print the returned events before returning the total:

```rust
for (group_index, group) in result.event_groups.iter().enumerate() {
    for event in group {
        eprintln!(
            "ARB EVENT group={group_index} cut={} orientation={:?} weight={:e}",
            event.cut_info.cut_id,
            event.cut_info.orientation_id,
            event.weight,
        );
    }
}
```

Finally, an isolated coefficient may be identically zero on one physical cut.
In `crates/gammalooprs/src/uv/approx/final_integrand.rs`, temporarily retain
that typed cut entry while running the coefficient oracle:

```diff
         let localized_local = match localized_local {
             Some(localized) => localized,
+            None => localizer.localize(&Atom::Zero, graph, current)?.combine()?,
-            None => {
-                return Err(eyre::eyre!(
-                    "factorized local term has no active UV sectors"
-                ));
-            }
         };
```

That unconditional zero fallback is diagnostic-only. Do not leave it in
production; outside this deliberately isolated coefficient, absence of active
UV sectors can indicate malformed construction and must remain an error.

Then run the same focused command as above. Repeat with coefficient `-2`, then
`-1`, then `0`. The authoritative expected pattern is the table in section 7:
T0 and T1 agree on every cut; T2 agrees on `[3,4]` and differs on `[1,2]` and
`[1,7,8]`.

### 23.4 Why this patch is preferable to restoring environment hooks

The patch makes the semantic change visible in the diff, keeps it confined to
an isolated checkout, and cannot silently affect a normal process run because
someone happened to export a variable. It also makes the meaning of `T0/T1/T2`
explicit at the two call sites. Once the first divergent boundary is captured,
remove the patch and validate the eventual generic fix through complete
residue sums in normal production code.

## 24. Embedded complete-T2 one-call oracle source

The checkout `/tmp/gl04-onecall.KO3L9q` and target
`/tmp/gl04-onecall-target` are ephemeral and nonportable. They may disappear
at any time. The source below is therefore the authoritative copy of the
unfinished one-call experiment. Append it temporarily to the existing
`#[cfg(test)]` module in `crates/gammalooprs/src/uv/approx/mod.rs` in an isolated
checkout. It intentionally expects the T2 coefficient-isolation patch from
section 23 and is `#[ignore]` so it cannot enter the ordinary suite by accident.

The test has not yet produced an authoritative result. Its purpose is to take
the complete actual projected T2 source—outer denominators and the five child
occurrences—in one generalized-CFF call. The outcome determines whether the
first failing boundary is already the raw reconstructed source/CFF call or the
later sequential child-by-outer composition. Compile or routing fixes may
still be needed; preserve the experiment's physical construction if repairing
the scaffold.

```rust

    #[test]
    #[ignore = "temporary complete GL04 T2 one-call oracle"]
    fn gl04_complete_t2_one_call_oracle() -> Result<()> {
        test_initialise()?;
        assert_eq!(
            std::env::var("GL_UV_DIAGNOSTIC_TAYLOR_LAURENT_POWER").as_deref(),
            Ok("0")
        );
        let mut graph: Graph = dot!(digraph gl04_complete_t2_one_call {
            edge [num=1 particle="scalar_0"]
            node [num=1]

            v0 -> v1 [id=0 is_cut=0 particle="scalar_1"]
            v0 -> v3 [id=1 lmb_id=0]
            v0 -> v5 [id=2]
            v1 -> v2 [id=3]
            v1 -> v3 [id=4]
            v2 -> v4 [id=5 lmb_id=1 num="Q(5,spenso::cind(0))^2"]
            v2 -> v4 [id=6]
            v3 -> v5 [id=7 lmb_id=2]
            v4 -> v5 [id=8 particle="scalar_1"]
        }, "scalars")?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let numerator = graph.production_numerator_atom_for_full_3d_expression();
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production_contract = graph.tree_edges.subtract(&graph.initial_state_cut);
        let production = graph.generate_3d_expression_for_integrand(
            &graph.paired_edges(&production_contract),
            &canonization,
            &options,
            Some(&numerator),
        )?;
        let target_support = BTreeSet::from([EdgeIndex(1), EdgeIndex(2)]);
        let lu_surface = graph
            .determine_raised_esurfaces_from_expression(&production.expression)
            .raised_groups
            .into_iter()
            .find(|group| {
                group.esurface_ids.iter().any(|surface_id| {
                    production.expression.surfaces.esurface_cache[*surface_id]
                        .energies
                        .iter()
                        .copied()
                        .collect::<BTreeSet<_>>()
                        == target_support
                })
            })
            .expect("the GL04 e1/e2 physical cut exists");
        assert_eq!(lu_surface.max_occurence, 1);
        let mut cutset = CutSet::empty(graph.n_hedges());
        cutset.residue_selector.lu = Some(LuCutSelection {
            raised_group: lu_surface,
            cut_edge_alternatives: vec![vec![EdgeIndex(1), EdgeIndex(2)]],
        });
        let child_filter = graph
            .get_edge_subgraph(EdgeIndex(5))
            .union(&graph.get_edge_subgraph(EdgeIndex(6)));
        let child_subgraph =
            InternalSubGraph::cleaned_filter_optimist(child_filter, graph.as_ref());
        let child_spinney = Spinney::with_scheme(
            child_subgraph,
            &graph,
            &graph.loop_momentum_basis,
            ApproximationType::MUV,
            2,
        )
        .expect("the GL04 e5/e6 bubble has DOD two");
        let orientation_pattern = OrientationPattern::default();

        let build = |mut route_graph: Graph,
                     projected_4d: bool|
         -> Result<(Atom, Vec<local_4d::FourDTerm>)> {
            let settings = UVgenerationSettings {
                generate_integrated: false,
                local_uv_cts_from_expanded_4d_integrands: projected_4d,
                add_marker: false,
                ..Default::default()
            };
            let localizer = Localizer::new(
                &cutset,
                OrientationProjection::exact_expression(
                    &production,
                    &options,
                    &orientation_pattern,
                    true,
                ),
            );
            let mut root = Approximation::new(Spinney::empty(&route_graph));
            root.root(&mut route_graph, localizer, &settings)?;
            let mut child = Approximation::new(child_spinney.clone());
            child.simple_approx = Some(
                root.simple_approx
                    .as_ref()
                    .expect("the root approximation exists")
                    .dependent(child.spinney.subgraph.clone()),
            );
            child.compute_4d(
                &route_graph,
                (crate::utils::vakint()?, &vakint::VakintSettings::default()),
                &root,
                &settings,
            )?;
            let physical_terms = child
                .local(&route_graph)?
                .active_sectors()
                .iter()
                .map(|sector| sector.physical_terms())
                .collect::<Result<Vec<_>>>()?
                .into_iter()
                .flatten()
                .collect::<Vec<_>>();
            child.compute_3d(&root, &mut route_graph, localizer, &settings)?;
            let selected = child
                .final_integrand(&route_graph)?
                .iter()
                .filter(|(index, _)| index.lu_cut_order == Some(1))
                .map(|(_, atom)| atom)
                .fold(Atom::Zero, |sum, atom| sum + atom);
            Ok((selected, physical_terms))
        };

        let (direct, _) = build(graph.clone(), false)?;
        let (projected, physical_terms) = build(graph.clone(), true)?;
        assert_eq!(
            physical_terms.len(),
            1,
            "T2 must be one additive rational term"
        );
        let child_term = &physical_terms[0];
        assert_eq!(
            child_term
                .denominators
                .iter()
                .map(|denominator| usize::from(denominator.source_edge))
                .collect::<Vec<_>>(),
            vec![5, 5, 6, 6, 6]
        );
        let outer_denominators = [1usize, 2, 3, 4, 7, 8].map(|edge| {
            let edge = EdgeIndex(edge);
            FourDDenominator {
                source_edge: edge,
                momentum: symbolica::atom::FunctionBuilder::new(GS.emr_mom)
                    .add_arg(usize::from(edge))
                    .finish(),
                mass_squared: graph.underlying[edge].particle.mass_atom().pow(2),
                full_expr: Atom::one(),
            }
        });
        let complete_denominators = outer_denominators
            .into_iter()
            .chain(child_term.denominators.iter().cloned())
            .collect::<Vec<_>>();
        assert_eq!(complete_denominators.len(), 11);
        let (complete_cff, _) = graph.cff_from_4d_denominators_in_uv_edges(
            &complete_denominators,
            [EdgeIndex(5), EdgeIndex(6)],
            &cutset,
            &options,
            &child_term.numerator,
            None,
        )?;
        let one_call = complete_cff
            .terms
            .iter()
            .filter(|(index, _)| index.lu_cut_order == Some(1))
            .flat_map(|(_, term)| {
                term.orientations.iter().map(|orientation| {
                    Ok(orientation.expression.clone()
                        * term.map_exact_source_numerator(&orientation.orientation)?)
                })
            })
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .fold(Atom::Zero, |sum, atom| sum + atom)
            * Atom::num(complete_cff.production_prefactor_factor());

        let split_momentum_replacements = graph
            .iter_edges_of(
                &graph
                    .full_filter()
                    .subtract(&graph.initial_state_cut)
                    .subtract(&graph.tree_edges),
            )
            .filter_map(|(pair, edge, _)| {
                pair.is_paired().then(|| GS.split_mom_pattern_simple(edge))
            })
            .collect::<Vec<_>>();
        let expressions = [direct, projected, one_call].map(|atom| {
            atom.replace_multiple(&split_momentum_replacements)
                .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
                .with(W_.d_)
        });
        let rational =
            |numerator, denominator| F::<ArbPrec>::from(&Rational::from((numerator, denominator)));
        let loop_moms: LoopMomenta<F<ArbPrec>> = [
            ThreeMomentum::new(rational(3, 10), rational(4, 10), rational(0, 1)),
            ThreeMomentum::new(rational(1, 10), rational(2, 10), rational(3, 10)),
            ThreeMomentum::new(rational(-2, 10), rational(1, 10), rational(1, 10)),
        ]
        .into_iter()
        .collect();
        let external_moms: ExternalFourMomenta<F<ArbPrec>> = [[
            rational(1, 1),
            rational(0, 1),
            rational(0, 1),
            rational(0, 1),
        ]
        .into()]
        .into_iter()
        .collect();
        let sample = MomentumSample {
            sample: BareMomentumSample {
                loop_moms,
                dual_loop_moms: None,
                loop_mom_cache_id: 0,
                loop_mom_base_cache_id: 0,
                external_moms,
                external_mom_cache_id: 0,
                external_mom_base_cache_id: 0,
                jacobian: rational(1, 1),
                orientation: None,
                parameterization_branch: None,
            },
        };
        let orientations = TiVec::<OrientationID, EdgeVec<Orientation>>::new();
        let orientation_filter = SubSet::full(orientations.len());
        let mut param_builder = graph.param_builder.clone();
        param_builder.m_uv_value(Complex::new_re(F(0.0)));
        param_builder.mu_r_sq_value(Complex::new_re(F(9.0)));
        param_builder.numerator_sampling_scale_value(Complex::new_re(F(1.0)));
        let (mut evaluator, _) = EvaluatorStack::new_explicit_sum_with_timings(
            &expressions,
            &param_builder,
            None,
            &EvaluatorSettings::default(),
        )?;
        let input = <ArbPrec as GenericEvaluatorFloat>::get_parameters(
            &mut param_builder,
            (false, false),
            &graph,
            &sample,
            &[],
            &[],
            None,
            None,
            None,
        );
        let values = evaluator
            .evaluate(
                input,
                SingleOrAllOrientations::All {
                    all: &orientations,
                    filter: &orientation_filter,
                },
                &RuntimeSettings::default(),
                &mut EvaluationMetaData::new_empty(),
                false,
            )?
            .into_iter()
            .map(DualOrNot::unwrap_real)
            .collect::<Vec<_>>();
        eprintln!(
            "GL04 COMPLETE T2 ONE CALL direct={:e} projected={:e} one_call={:e} projected-direct={:e} one_call-direct={:e} one_call-projected={:e}",
            values[0],
            values[1],
            values[2],
            values[1].clone() - &values[0],
            values[2].clone() - &values[0],
            values[2].clone() - &values[1],
        );
        Ok(())
    }
```

## 25. Clean checkpoint gate

After removing the temporary diagnostics described above, the curated
warnings-as-errors gate was classified exhaustively across all 1,967 selected
tests: 1,937 pass, 29 fail, and one times out. This includes a separate run of
the 12 tests withheld by nextest's failure cap; 11 of those pass, including all
four DDx/TTx NLO acceptance tests, while `integrations::v_diag` is the sole
additional failure. See
`HANDOFF_REMOVE_FOR_PR/CHECKPOINT_TEST_STATUS.md` for the exact commands, every
failing test, the timeout, and the residual-set status table.
