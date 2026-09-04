# T2 sign-debug plan

> Temporary handoff material committed under `HANDOFF_REMOVE_FOR_PR/`.
> Remove this directory before the merge-ready PR is finalized.

This is the live debugging plan for the remaining local-UV discrepancy exposed
by the GL04 scalar LU test. It must be reread before changing either the
projected local-4D lane or generalized CFF sign conventions, and updated after
each decisive experiment.

## 2026-09-04 correction: numerator Taylor migration is not EMR dispatch

The current GL04 `T2` discrepancy does **not** come from the dispatcher for
denominator-derivative numerator factors. These are two separate mechanisms:

- When a Taylor derivative hits a denominator, that derivative creates both a
  new serial denominator occurrence and a new temporal numerator factor. The
  dispatcher assigns that factor to the new, initially unoccupied occurrence.
  This path does not create a quadratic power on an already occupied edge.
- In the isolated GL04 `{1zs}` frame the actual routing is `Q5=q` and
  `Q6=-q+p`, with the retained soft shift `p=Q3`. The original numerator
  `(Q5^0)^2` therefore remains a pure hard-child factor. The second Taylor
  coefficient of the shifted denominator `D(-q+t p)` contains both a
  `p^2/D(q)^2` term and a `(q.p)^2/D(q)^3` term. These create a genuine
  quadratic outer dependence `(Q3^0)^2`.
- This still does not pile two hard energies onto one derivative-created
  occurrence. Every new hard `q^0` factor from `q.p` is assigned to its new
  serial child edge. The `p^0=Q3^0` factors are external to the vacuum child;
  they remain on the factorized outer graph and are not candidates for the
  child EMR dispatcher.

Consequently, the projected local-4D outer CFF must—and does—request the bound
`(3,2)`. The direct lane starts from the complete generalized root CFF for the
original `(Q5^0)^2`, whose physical bound is `(5,2)`, and then applies the
Taylor operator to each complete residue-map branch. Whether that root must
also be provisioned for the newly exposed outer `(3,2)` dependence is a diagnostic
question about closure of the generalized functional under the Taylor
operation; it is not a reason to dispatch the original numerator or attach
its factors to derivative-created edges.

The first true-root experiment added `(3,2)` alongside `(5,2)` in an isolated
checkout. All three per-cut Arb values and route differences were unchanged
digit-for-digit. A second supported experiment used `(3,2),(5,2),(8,2)`, thus
covering the original hard owner and both physical carriers of the retained
soft routing; it was again digit-for-digit identical to baseline. Requesting
degree two on every internal edge was attempted, but the current Rust
generalized-CFF port rejects that multi-edge sector before evaluation. The
complete evidence-directed envelope is therefore supported and inert: missing
root capacity is not the explanation.

The cut discriminator is exact:

| physical cut | edge 3 cut? | projected `Q3^0=0` contact maps | `T2` match |
| --- | --- | ---: | --- |
| `[1,2]` | no | 6 | fails exactly |
| `[3,4]` | yes | 0 | matches at the 1000-bit floor |
| `[1,7,8]` | no (`e8`, not `e3`, is cut) | 4 | fails exactly |

`T0` and `T1` match on every cut. At `mUV=0`, the two failing complete-cut
differences are exactly proportional to their complete projected `Q3^0=0`
contact families. This localizes the unresolved boundary to the rank-two
outer-energy part of the original numerator's `T2` coefficient; it does not
implicate the generic derivative-factor dispatcher.

The direct root-key placement has also been traced exactly. For the quadratic
`Q5^0` bound its generalized samples are `-OSE5`, `0`, and `+OSE5` (this case
does not use the uniform `n*M` samples reserved for higher degree). The mapped
numerator is multiplied into the complete branch before the OSE rewrite,
rescaling, and Taylor series. In particular the contact finite difference is

```text
[2 N(-E5) - 4 N(0) + 2 N(+E5)] / (2 E5)^2 = 1
```

for `N(q0)=q0^2`. CFF variant prefactors, numerator surfaces, half-edge energy
factors, residual denominators, and the mapped numerator are all inside the
Taylor operation; only the opaque residue-key selector remains outside. No
incorrect ordering or missing quadratic root key has been found in this
direct family. The next authoritative test is one generalized-CFF call on the
complete actual factorized `T2` rational source, compared only after summing
the full physical residue.

## Non-negotiable construction contracts

1. **Local UV from 3D:** apply the UV Taylor operator directly to the complete,
   residue-key-local generalized CFF expression. This lane performs no UV graph
   reconstruction. Its explicit-orientation-sum form is the same expression
   with the residue-key selectors removed.
2. **Local UV from 4D:** apply the Taylor operator to the factorized 4D
   integrand, reconstruct the resulting UV skeleton from retained physical-edge
   provenance, express its numerator in that skeleton's EMR variables, and call
   generalized CFF. This lane is residue-summed and supports only explicit
   orientation sums.
3. Compare the two lanes only after summing the **complete residue** (all its
   generalized pieces) for a fixed Cutkosky cut. Individual contact/remainder
   pieces are not expected to agree.
4. The original numerator remains factorized and attached to its original
   physical owner. The EMR dispatcher may move only temporal factors created by
   denominator derivatives, and only among newly created degenerate copies of
   that same physical line.
5. A derivative-created denominator power is represented by a new serial edge.
   Repeated derivatives of one denominator create indistinguishable temporal
   factors; no derivative-event identity is needed. Deterministically distribute
   the degree-one factors over the new serial copies, reserving the canonical
   base occurrence for the original numerator.
6. Do not algebraically cancel numerator and denominator factors upstream of
   generalized CFF. The generalized residue map realizes contact/lower-sector
   terms.
7. Repeated denominator occurrences remain distinct topology edges, distinct
   EMR IDs, and distinct orientation/residue-key slots. Generalized CFF may—and
   for confluent residues must—group them internally as one logical energy
   channel and sample them on their common LMB diagonal. That internal
   correlation must never collapse the public edge/key structure.
8. Different UV-forest terms can legitimately have residue-map keys of
   different lengths because their reconstructed denominator multiplicities
   differ. The local-3D explicit-sum route ultimately sums those residues, so no
   artificial common key length is required beyond the Cutkosky/threshold
   construction.

## Current status: one cut resolved, two cuts remain

The direct local-3D and projected local-4D constructions agree for the complete
GL04 `{1zs}` Taylor-`T2` contribution on the physical cut `[3,4]`, once the
comparison is made on the exact Cutkosky surface and after summing every
generalized residue piece for that cut. Impose

\[
Q^0 = E_3 + E_4.
\]

At an exact rational/radical Arb point, all twelve complete outer host groups
agree separately, and the full 48-branch direct sum equals the full 12-branch
projected sum:

```text
direct = projected = -5.43331072437902428657...e-8 i
projected - direct = -7.79e-308 i
relative difference = 1.43e-300
```

The earlier interpretation then made a bookkeeping mistake: the event's
`cut_info.cut_id` is not the `CutGroupId`. The exact CutGroup/event permutation
is

```text
CutGroupId 0: [1,2]   -> event cut_id 2
CutGroupId 1: [3,4]   -> event cut_id 0
CutGroupId 2: [1,7,8] -> event cut_id 1
```

Thus the exact `[3,4]` proof certified the cut that was already green. It did
not explain the two live discrepancies. At 1000-bit Arb, before the common
ordinary-cut pass-two normalization, the complete pass-one values are

```text
[1,2]:
  direct    = -9.7358832009585117881199719792504e-20 i
  projected = -5.8222698776306930273704916168530e-20 i

[3,4]:
  direct = projected
         = +4.3997865139100253906601972010589e-19 i  (to ~1e-300)

[1,7,8]:
  direct    = +9.7220545326254829481117326753448e-20 i
  projected = +9.7222977344583806241029024530620e-20 i
```

The Arb Newton residuals, rescaled `MomentumSample`, and production FunctionMap
surface evaluation remain at the roughly `1e-301` rounding floor. Therefore
the live discrepancy is not inherited f64 energy-momentum nonconservation, an
inaccurate `t*`, rescaled-momentum construction, or pass-two contraction. The
current MRE is the complete `[1,2]` pass-one expression: 56 direct branches
versus 20 projected branches, always compared only after summing the complete
physical residue.

The child-only oracle independently agrees at 1000-bit Arb precision for
`p0={0,+1,-1,+sqrt(13)/10,-sqrt(13)/10}` and is even under `p0 -> -p0`.
This certifies the complete child function, rather than a single accidental
sample.

### Production parameter-vector and cut-surface audit

The native-Arb production evaluation does **not** inherit the small violation
of the cut equation present in the original double-precision evaluation.  A
focused trace recorded 36 `ParamBuilder` records, each containing the complete
ordered vector of 56 parameter slots.  The final native-Arb direct-local3D
records are byte-for-byte identical to the corresponding projected-local4D
records: this includes `t*`, the external four-momentum, all three rescaled loop
momenta, `mUV`, `mUVexp`, `mu_r^2`, `M`, the model parameters, and every
residue/orientation slot.  The preliminary native-Arb stability-pass records
are likewise identical between the two routes.  There is therefore no first
parameter-vector mismatch to explain the route difference.

Reconstructing the three physical surface equations directly from the logged
production values gives

```text
[1,2]:   2 |Q1| - Q0                                    = O(1e-301)
[3,4]:   2 |Q1-Q7| - Q0                                 = O(1e-301)
[1,7,8]: |Q1| + |Q7| + sqrt(|Q1-Q7|^2 + 0.1^2) - Q0   = O(1e-301)
```

For the two evaluated rotations, the six residuals reconstructed from the
decimal trace are respectively `4.62e-302`, `1.30e-301`, `1.34e-301`,
`-6.59e-302`, `-4.07e-302`, and `9.27e-302`.  These are at the 1000-bit Arb
rounding/logging floor.  The analogous residuals reconstructed from the f64
records are of order `1e-17`; the precision escalation therefore recomputes
`t*` and the rescaled momenta natively and removes that f64 defect rather than
upcasting it.  The exact-on-cut standalone oracle uses a different valid
rational/radical sample point, but its relevant contract is the same exact
surface equation, which production satisfies to Arb precision.  Consequently
neither an inaccurate native-Arb cut solution nor different runtime inputs in
the two UV routes is the current explanation.  The existing `mUV=5,20,80`
scan changes no kinematic or `t*` input and is consistent with this exclusion.

## Superseded off-surface discrepancy record

The section below is retained as investigation history.  Its numerical
differences are diagnostic values away from the physical cut surface and must
not be treated as failures of either UV route.

The production GL04 control uses the original numerator

\[
N=(Q_5^0)^2.
\]

With integrated UV counterterms and threshold counterterms disabled, the first
remaining mismatch is in the second-order Taylor contribution. Its projected
local-4D denominator multiset is

```text
owners:      5, 5, 6, 6, 6
energy IDs:  9,12,10,11,13
```

The reconstructed exact source therefore contains two serial copies owned by
edge 5 and three serial copies owned by edge 6. The pre-dispatch and
post-dispatch numerator expressions agree symbolically after both are written
in the same LMB. The current full rank request is sufficient (it contains rank
two where required), so the discrepancy is not explained by simply omitting a
rank-two contact.

At the diagnostic phase-space point, the two complete results differ only by a
small recursive-contact subfamily:

```text
projected local 4D: F + C
direct local 3D:    F - C
difference:         2 C
```

Forcing the two simple derivative-created temporal factors onto two distinct
new occurrences changes the numerical result only below roughly `1e-306`; it is
still the correct structural default, but it does not by itself cure this sign.
Disabling the exact-source cache likewise leaves the mismatch unchanged.

## Immediate controls

### A. Spatial-quadratic numerator

Replace `(Q_5^0)^2` by a product of two spatial components, for example
`Q_5^1 Q_5^1`. This preserves the overall UV degree but removes the raised EMR
energy bound and its generalized contact interpolation.

- If direct 3D and projected 4D disagree, the fault is in T2 topology/sign or
  source normalization and is not specific to temporal contacts.
- If they agree and both UV profiles pass, return to `(Q_5^0)^2` and isolate the
  generalized contact path.
- Compare complete sums per Cutkosky cut, not individual event/residue pieces.

Current result with integrated and threshold counterterms disabled:

```text
orientation-local direct 3D = explicit-sum direct 3D
  -6.6390634494155886293366743389748287e-10 i
projected local 4D
  -6.6390534378633180985233346720561283e-10 i
projected - direct
  +1.0011552270530813339666918700e-15 i
```

Both routes have the same sign. The remaining spatial difference is at the
known double-precision-input absolute floor, whereas the temporal-square T2
difference is about three orders of magnitude larger. This rules out an
overall projected-T2 sign reversal and sharpens the temporal/contact path as
the material GL04 discrepancy.

### B. Powered-denominator cancellation identity (resolved)

The smaller no-forest identity

\[
\frac{D(Q)}{D(Q)D(-Q)}=\frac{1}{D(-Q)}
\]

now agrees in the raw generalized engine and in every GammaLoop consumption
lane. The earlier opposite production sign was caused by applying an ordinary
CFF core sign a second time while converting retained positive `1/(2E)`
factors to GammaLoop's historical negative-energy convention. The corrected
typed conversion uses `(-1)^N S_den` for a generalized child and multiplies the
complete Laurent functional uniformly. This test contains no local-4D
reconstruction, Taylor forest, or threshold/integrated counterterm, and now
exonerates generalized contact generation and the adapter boundary. It cannot
be reused to fit the remaining GL04 contact-relative discrepancy.

### C. Algebraically equivalent temporal oracle (only if A and B are
inconclusive)

Use

\[
(Q^0)^2 = D(Q)+|\vec Q|^2+m^2
\]

and cancel the explicit `D(Q)` in a handcrafted diagnostic representation.
The resulting sum uses an ordinary lower topology plus a contact-free spatial
numerator. It must locally equal the uncancelled generalized-CFF result. This is
an independent direct test of generalized contact generation, but it is more
work and should not precede the simpler controls.

## Routing-sign audit

For each same-owner serial group, record separately:

1. the immutable hard/source momentum `H` stored in the provenance tag;
2. the raw rewritten denominator momentum `R` in the Taylor term;
3. the parsed serial-edge momentum `P` after topology canonicalization;
4. `hard_to_raw_sign` in `H = h R`;
5. `raw_to_parsed_sign` in `R = r P`;
6. the final odd-numerator conversion `H^0 = h r P^0`.

Because `D(Q)=D(-Q)`, canonicalizing a denominator from `Q` to `-Q` is valid,
but the sign cannot be discarded for odd temporal numerator factors. Every
serial copy created from one physical denominator should use one consistent
parsed signature. Prefer constructing the whole chain as `(Q,Q,Q)` (or all
`-Q` when anchored to the original edge) rather than permitting a mixed
`(Q,Q,-Q)` group. If a mixed spelling is merely the raw momentum of distinct
physical owners, do not incorrectly force those owners into one group; compose
their individual `h` and `r` signs instead.

The production experiment is:

1. identify whether the reported `Q,Q,-Q` belongs to one same-owner serial
   group or combines owner-5 and owner-6 groups;
2. assert all parsed signatures inside each same-owner group are identical;
3. assert every emitted orientation and edge-energy map retains a slot for each
   occurrence;
4. if same-owner parsed signatures are mixed, build the serial chain from one
   representative signature and retain per-occurrence raw-to-parsed signs only
   in numerator mapping;
5. rerun the exact common-LMB numerator certificate, the powered-denominator
   identity, the spatial GL04 control, and the temporal GL04 T2 comparison.

## GammaLoop/generalized-CFF normalization audit

Do not retain a sign because one numerical identity prefers it.  GammaLoop and
`three-dimensional-reps` use two explicit on-shell-energy conventions at their
adapter boundary:

- a generalized variant owns positive factors `1/(2 E_i)` locally;
- GammaLoop's ordinary CFF evaluator removes those factors and restores the
  historical global product `1/prod(-2 E_i)`.

For one rational component with `N` denominator-energy factors, let `S_core`
be its ordinary CFF core conversion and `S_den` the uniform sign of the same
scalar denominator source. The adapter conversion is determined by two typed
facts: whether the aggregate retained positive local energy factors and
whether this component is ordinary or generalized. It is

```text
all-ordinary aggregate:
    B_i = S_core_i

variant-local aggregate, ordinary child component:
    B_i = (-1)^N_i * S_core_i

variant-local aggregate, generalized child component:
    B_i = (-1)^N_i * S_den_i

S_GL<-3d = product_i B_i .
```

The `(-1)^N` factor appears exactly when GammaLoop retains the generator's
positive local `1/(2E_i)` factors instead of replacing them by its historical
global `1/prod(-2E_i)` product. An ordinary child retains its established core
conversion. A generalized child uses `S_den`, the numerator-bound-independent
scalar reference. Its generalized `S_core` already belongs to the raw Laurent
functional and must **not** be multiplied again.

This conversion must be tested without UV reconstruction or a nontrivial
numerator: generate the same physical graph once with ordinary scalar capacity
and once with unused rank-two capacity, evaluate numerator `1`, and require
exact equality. The native `three-dimensional-reps` two-loop control already
does so before GammaLoop conversion and agrees at `5e-17`, exonerating graph
reconstruction, residue-map evaluation, and generalized contact algebra. The
corresponding GammaLoop two-loop control proves that multiplying the
generalized core a second time flips the result. The conversion is therefore
ownership-dependent but route-independent: persisted production, `Graph::cff`,
and reconstructed exact-4D sources all use the same adapter rule.

For the GL04 T2 five-denominator parent and its four-denominator lower contact,
this ratio is `-1` in both cases.  Consequently it cannot change the relative
contact sign `F+C` into `F-C`; it acts on the complete Laurent functional.

## Precision and cut reduction

The stripped GL04 temporal-square comparison has been evaluated through the
production evaluator at both 256-bit and 1000-bit Arb precision. The total
route difference is

```text
256 bit  = +2.262532013744616378263364907633471491633795616772479646570887905061893401622075e-12 i
1000 bit = +2.262532013744616378263364907633471491633795616772479646570887905061893572314082...e-12 i
```

Roughly seventy significant digits are stable. At 1000 bits the route values
are

```text
explicit-sum direct 3D = -9.1734755805810912839071256602802601e-6 i
projected local 4D      = -9.1734733180490775392907473969153524e-6 i
projected - direct      = +2.2625320137446163782633649076334715e-12 i
relative difference    =  2.4663847348480072054e-7
```

The same leading values are seen in the ordinary evaluation, so this is not
roundoff introduced by the Arb evaluator. The events returned by the
cross-section evaluator are already complete fixed-Cutkosky residues for one
evaluator/LMB group when thresholds are disabled; they are not individual
generalized contact pieces. Aggregating all groups by cut gives

```text
cut 0 projected-direct: +1.513e-302 i                 (numerical zero)
cut 1 projected-direct: +5.11064462007056688928949e-9 i
cut 2 projected-direct: +5.46868617282760250236730e-14 i
```

Cut 1 is therefore the primary MRE; cut 2 is a stable secondary oracle about
`9.35e4` smaller in event-weight space. The next reduction must key complete
events by `(group, cut_id)`, start from cut 1, and only then isolate the first
unequal UV-forest node and its Taylor order. Positional event zipping is
diagnostic only and must not be mistaken for residue-piece equality.

## Acceptance gates after the sign is resolved

1. Focused Arb-to-Arb exact-source identities and common-LMB numerator
   certificates.
2. GL04 temporal and spatial controls with integrated and threshold terms first
   disabled, then enabled.
3. Complete scalar LU three-way matrix. Orientation-local local-3D must also
   pass per-residue-key UV profiling; projected local-4D is profiled after the
   complete residue sum.
4. DDx and TTx NLO acceptance tests, scale-independence checks, full
   `just test_gammaloop`, performance report, cleanup, commit/push, and PR.

## Investigation log

- 2026-09-04: Corrected a decisive cut-index bookkeeping error. `CutGroupId`
  and generated-event `cut_id` are permuted as `0->2`, `1->0`, `2->1` for
  GL04. The exact on-surface 48-versus-12 proof covered group 1 `[3,4]`, which
  is already green. The actual pass-one failures are group 0 `[1,2]` and group
  2 `[1,7,8]`; pass two only preserves their discrepancies.
- 2026-09-04: Exonerated the proposed Arb Cutkosky-root precision failure.
  Newton, the rescaled sample, and the production FunctionMap all satisfy the
  selected cut equations at about `1e-301`. A separate generic promotion bug
  was fixed: static `ParamBuilder<f64>` values now lift directly to Arb instead
  of passing through Quad and inheriting its roughly 32-decimal-digit ceiling.
  Its focused production-path regression passes, but it does not change the
  GL04 discrepancies.

- 2026-09-04: Confirmed exact-source topology retains five distinct serial
  edges for the GL04 T2 term. Internal logical-channel rank pooling correlates
  their samples but does not merge their EMR or orientation slots.
- 2026-09-04: Confirmed cache removal and a forced two-copy temporal dispatch do
  not change the GL04 T2 discrepancy.
- 2026-09-04: Adopted production candidate semantics that reserve the base
  occurrence for the original numerator and distribute derivative-created
  degree-one factors over the additional copies. Focused tests are being run.
- 2026-09-04: Reproduced the smaller `D/D^2 -> 1/D` sign discrepancy on current
  source. The exact powered source equals the exact lower source; only the
  graph-generated generalized production path has the opposite sign. This
  points to a production CFF normalization boundary independent of projected
  T2 graph reconstruction.
- 2026-09-04: Audited the actual serial-chain constructor. Occurrences are
  grouped by immutable `source_edge`; one parsed signature is chosen for that
  owner and cloned onto every member of its serial chain. In the GL04 T2 term,
  owner 5 is raw `(+Q,+Q)` and owner 6 is raw `(-Q,-Q,-Q)`. The previously
  quoted cubic `(Q,Q,-Q)` fixture contains three *different* physical owners
  (4, 5, and 7), so it is not evidence that a same-owner chain has mixed
  routing. The next check is at the generated residue-map layer, where each
  distinct occurrence can still acquire its own contour orientation while its
  momentum signature remains common.
- 2026-09-04: The first focused dispatcher build in the shared target failed at
  link time with stale/corrupted incremental LLVM objects, not a Rust/test
  failure. Re-run this test in an isolated target before drawing a conclusion;
  do not change physics code in response to the linker artifact.
- 2026-09-04: The spatial-square GL04 control gives the same negative sign in
  both direct-3D modes and projected local 4D. Its absolute residual is about
  `1.0e-15`; the material T2 discrepancy is therefore specific to temporal
  contact structure rather than an overall forest/source sign.
- 2026-09-04: The GL04 T2 trace confirms the reconstructed parsed topology is
  already `(Q,Q,Q,Q,Q)`: owner-5 raw occurrences 9 and 12 map with sign `+1`,
  while owner-6 raw occurrences 10, 11, and 13 start as `-Q` and each map with
  `raw_to_parsed=-1` into that same `+Q` frame. Thus the requested same-sign
  serial construction is already the production behavior. The independent
  plus/minus values in a residue-map key are contour choices, not denominator
  momentum signatures, and must not be forced equal.
- 2026-09-04: Made the same-routing invariant explicit in production code by
  removing the unused per-chain sign parameter from
  `serial_power_chain_incidences`. The helper can now construct only one
  consistently tail-to-head serial chain; a single group-level signature
  canonicalization handles `Q` versus `-Q` for all copies together. This is a
  semantics-preserving simplification of behavior already observed in GL04.
- 2026-09-04: Found a bookkeeping weakness in the first reserve-base dispatch
  patch: it inferred the base occurrence from the smallest canonical energy ID.
  Canonical edge ordering can interleave different owners, so the base role
  must instead be carried explicitly through occurrence canonicalization.
  Replace that inference with one explicit base plus explicit derivative-copy
  IDs per physical owner before treating the dispatcher patch as complete.
- 2026-09-04: A provisional route-specific normalization conclusion was
  falsified by a same-graph rank-capacity control. Persisting a generated CFF
  does not change its normalization, and an unchanged exact-denominator source
  can be literally the same generated object as a physical source. The
  conversion therefore cannot depend on object lifetime or the UV route.
- 2026-09-04: A two-loop unit-numerator control isolated the flaw in the first
  componentwise formula. With three denominator energies it gave ordinary
  `+1/(76800*pi^6)`, the formula with an extra generalized-core sign gave the
  negative, and `(-1)^N*S_den` gave the ordinary value exactly. A native
  `three-dimensional-reps` comparison of the same scalar source with ordinary
  and unused rank-two capacity agrees before GammaLoop conversion. The defect
  is therefore the GammaLoop adapter's asymmetric treatment of global versus
  variant-local energy factors, not generalized residue generation, EMR
  dispatch, or UV graph reconstruction.
- 2026-09-04: The complete exact GL04 R/T2 child satisfies the symbolic identity
  `q0^2 = D(Q) + |q|^2 + m^2` with a positive lower-topology contact:
  `original - spatial_mass_remainder - lower_contact == 0`, whereas the
  opposite contact sign fails. Thus generalized contact generation and the
  exact-child normalization bridge are correct at that boundary. The remaining
  projected mismatch is downstream of complete child generation/mapping.
- 2026-09-04: The stripped spatial-square GL04 projected lane also fails UV
  profiling (initial expected DOD 2 fits approximately `+2.000028`; the
  subsequently expected DOD 0 fits approximately `-0.000624`). Its tiny
  `1.0012e-15` pointwise difference is therefore accidental at the chosen
  sample, not a proof of correctness. The projected assembly defect is more
  general than temporal contact interpolation.
- 2026-09-04: Replaced the provisional lowest-ID base inference with explicit
  occurrence provenance. Each owner now carries exactly one `is_base` marker
  through topology canonicalization; unprovenanced/fixed factors may use only
  that base, while `DenominatorDerived` factors may use only the remaining
  serial copies. A derived factor with no created copy is an error. Canonical
  energy IDs may interleave across owners, so neither numeric order nor map
  order participates in this decision. This preserves the invariant
  `H^0 = h r P^0`: provenance role determines occurrence eligibility only,
  while every selected fixed or derived occurrence obeys the same routing-sign
  conversion.
- 2026-09-04: Corrected positive `GS.den(...)` assignment planning to read the
  EMR provenance embedded in its factorized value. An untagged original
  denominator completion and a fixed original factor stay on the explicit base
  occurrence; only a `DenominatorDerived` completion may use a
  derivative-created serial copy. Mixed provenance inside one coherent wrapper
  is rejected. The focused ownership tests and the unchanged Arb
  powered-rational identity pass.
- 2026-09-04: The stripped GL04 temporal-square residual survives the
  production 256-bit and 1000-bit Arb evaluators with the values recorded in
  the precision section above. Cut 0 is exonerated; cut 1 is the primary MRE
  target and cut 2 is secondary. No further normalization or generalized-CFF
  sign change is permitted unless an identical-input raw-engine comparison
  fails.
- 2026-09-04: Replaced the obsolete off-common-diagonal GL04 T2 formula in the
  exact-source certificate with the physical assignment ledger. The original
  `Q_5^0` factors use base occurrence 9; owner-5 derivative factors use only
  new occurrence 12; owner-6 derivative factors use both new occurrences 11
  and 13 and never base occurrence 10. The mapped numerator then agrees
  symbolically with the post-T numerator after every serial copy is expressed
  in the common loop coordinate. The independent denominator certificate and
  the complete raised-LU child comparison are also green. This removes EMR
  dispatch, rank capacity, and same-line signature conversion from the first
  suspect boundary.
- 2026-09-04: Audited projected outer-factor attachment in the production Arb
  trace. Every pre-outer active coefficient is cloned identically over its 216
  production hosts; each authoritative non-None source-energy map emitted by
  the outer CFF survives `zip_mul_mapped_factor` exactly once. Rehosting while
  retaining those maps cannot change the mapped value. The stable residual is
  therefore not explained by selector-host pairing or a dropped Cartesian
  branch.
- 2026-09-04: Isolated the single `{1zs}` DOD-two forest node and then selected
  its Laurent layers through both production routes. The leading `T0` layer is
  exactly equal on all three physical cuts. The `T1` layer agrees to at least
  about 29 relative decimal digits on cut 1 and about 300 digits on cuts 0 and
  2. The finalized `T2` layer alone reproduces the complete node discrepancy:
  `+5.704105674507342249085253385...e-22 i` on cut 1 and
  `+1.989051502902968035605313928...e-16 i` on cut 2, while cut 0 is zero to
  Arb noise. Therefore “T2 matches” must be stated narrowly: the exact 4D
  coefficient, provenance-preserving EMR lift, serial topology, rank request,
  and child generalized-CFF normalization all pass their symbolic
  certificates. The first still-unequal object is the *completed* T2
  contribution after that correct child is embedded in the untouched outer
  CFF and the complete physical-cut residue is summed. Do not reopen T0,
  dispatcher ownership, graph reconstruction, or the uniform normalization
  conversion without a new counterexample at one of those earlier boundaries.
- 2026-09-04: The finalized T2 outer CFF requests exactly the required
  quadratic soft-energy capacity: every `{1zs}` call reports
  `physical_parent_bounds=[(3,2)]` and
  `assigned_cff_source_bounds=[(3,2)]`. Missing or truncated outer rank is not
  the cause. Bypassing `localized_orientation_terms` completely and hosting
  every raw reduced-outer CFF branch under one bookkeeping ID leaves every
  per-cut Arb value and discrepancy unchanged. Valid-ID filtering,
  representative choice, and source-map rehosting are therefore also
  exonerated.
- 2026-09-04: Varied the final production-path diagnostic at fixed kinematics.
  On cut 1 the projected-minus-direct T2 delta is approximately
  `3.6468503e-20 i`, `5.7041057e-22 i`, and `8.9132420e-24 i` for
  `M_UV=5,20,80`, respectively: it scales as `M_UV^-3`, while the complete T2
  contribution scales as `M_UV^-1`. Its relative size therefore scales as
  `M_UV^-2`. Scaling the whole nine-coordinate sample by `0.4` amplifies the
  cut-1 delta to `1.3001942e-10 i` while retaining the roughly `2.5e-5`
  relative discrepancy. This is a structural lower/contact contribution, not
  floating-point noise. Both explicit `M_UV^2` Taylor terms are present in the
  exact coefficient and survive the symbolic common-LMB EMR-lift certificate;
  wholesale loss of the leading mass sector would instead contaminate the
  `M_UV^-1` behavior.
- 2026-09-04: Forced `CffGenerationContext::Standalone` for the finalized
  reduced outer T2 source. Every complete-cut value and route delta was
  unchanged. This is required by the implementation: the GL04 outer source is
  a nonterminal generalized source, for which `Standalone` and
  `EmbeddedCffFactor` both invoke the same bounded CFF builder. The context
  distinction is therefore excluded from this MRE.
- 2026-09-04: Varied the runtime numerator sampling scale through the actual
  scalar-LU `ParamBuilder` at `M=0.5,1,2`. Direct, projected, and their per-cut
  differences were bit-for-bit unchanged. The generalized sampling-scale
  interpolation cancels as it should, so neither the `M` parameter slot nor a
  missing uniform-scale normalization explains the T2 residual.
- 2026-09-04: Audited the batched exact-CFF capacity envelope. The cache may
  redistribute a total rank over occurrences belonging to one certified
  algebraic repeated-energy channel; generalized_3drep then detects that
  channel, sums its member bounds, and samples it as one common energy
  variable. Existing higher-rank tests cover different distributions,
  reversed routing, and lowering. A diagnostic componentwise-max envelope
  (which over-provisions each occurrence separately) leaves the GL04 result
  and the `2.262532e-12 i` total residual unchanged, as did a genuine uncached
  generation. Do not mistake this capacity metadata for a lost contact term;
  retain the lower-cost pooled envelope unless a counterexample involves
  genuinely independent algebraic channels.
- 2026-09-04: The production Taylor/EMR/CFF pipeline is structurally generic in
  order. Negative denominator powers create an arbitrary occurrence vector,
  the Laurent truncation retains every coefficient through `t^0`, derived
  temporal factors are assigned over an arbitrary candidate vector, and CFF
  capacities and residue orders are dynamic. After the present mismatch is
  closed, add an explicit production-shaped T3 regression with one base and
  three derivative-created serial copies to guard that contract.
- 2026-09-04: Reproduced the isolated GL04 `{1zs}` T2 discrepancy through the
  actual cross-section evaluator and compared the two routes' complete
  `ParamBuilder` inputs. Their ordered symbol/value vectors agree exactly,
  including `mUV=20`, `mUVexp=20`, `rls=1000`, `mu_r^2=9`, `M=1`, every
  residue/orientation slot, model parameter, external momentum, and loop
  spatial component. Their complete replacement maps agree exactly as well,
  so the `Q` and `OSE` definitions are the same. Both are explicit residue
  sums, hence the runtime residue-map ID and orientation slots are inert.
  Inspection of every projected evaluator function map found no
  `projected_cff_sum` call: this production MRE is already materialized by
  `Localizer::projected_cff`, which builds ordinary `Integrands`; deferred
  bodies currently belong to the test-only whole-projection path. Therefore
  neither a deferred-function boundary nor `ParamBuilder` accounts for the
  Arb-stable cut-1/cut-2 residual. The previous `M=0.5,1,2` scan independently
  excludes a misplaced numerator-sampling-scale value, but it alone would not
  have excluded the rest of the parameter schema; the exact full-vector and
  replacement-map comparison now does.
- 2026-09-04: Distinguished an invalid per-order outer-moment probe from the
  authoritative GL04 failure.  The synthetic outer-basis test originally
  forced DOD two on a unit-numerator bubble whose physical DOD is zero, then
  compared the separate LU1/LU2 representatives at a generic point.  Those
  representatives are not individually invariant: terms proportional to the
  selected surface can move between LU1 and LU2 without changing the raised
  residue.  After replacing the child by a genuine DOD-two temporal numerator,
  a one-call exact four-dimensional CFF agreed with the sequential projected
  construction for every tested representative; the remaining odd LU2-only
  difference therefore establishes neither a contraction defect nor a
  production mismatch.  It must be tested only through the complete raised
  residue.
- 2026-09-04: The existing isolated GL04 `{1zs}`/T2 result is already such a
  complete physical-residue comparison and must not be dismissed with the
  per-order probe.  `evaluate_sample_precise` calls
  `CrossSectionGraphTerm::evaluate`; for each physical cut group that method
  iterates all `num_esurfaces = 1..=max_occurrence`, evaluates the matching
  `CutCFFIndex::lu_cut_order`, passes the required HyperDual integrand and
  surface derivatives to `pass_two_evaluators[order - 1]`, and accumulates
  every returned residue contribution in `bare_cut_total`.  With threshold
  subtraction disabled, the generated Arb event weight is exactly that full
  sum; the graph result then sums all physical cuts.  Consequently the
  recorded cut-1 and cut-2 GL04 deltas remain the authoritative complete-T2
  blocker, whereas raw LU-order values and isolated basis moments are only
  diagnostics.

## Audit of the earlier "standalone/emulated T2 oracle"

The earlier passing generalized-CFF test closest to the intended GL04 shape is

```text
cff::tests::exact_raised_lu_residue_factorizes_from_quadratic_cubic_spectator
```

It is **not** an identical-input oracle for the live GL04 `{1zs}` Taylor-
`T2` term.  A second passing test,
`graph::three_d_source::gl04_t2_planned_lift_matches_post_t_numerator_in_common_loop_coordinates`,
does use the exact GL04 Taylor coefficient, but it stops before CFF generation.
Together they certify two adjacent boundaries; neither certifies their
composition.  The field-by-field distinction is:

The freshly extracted production LMB used in this audit is

```text
e1 =  K0
e2 = -K0 + P0
e3 = -K2 - P0 + K0 = p
e4 = -K0 + K2
e5 =  K1               = q
e6 = -K1 - K2 - P0 + K0 = -q + p
e7 =  K2
e8 = -K2 - P0 + K0     = p
```

Thus `p` is not an invented bridge: it is the literal soft part of edge 6 in
the child sub-LMB and the physical EMR momentum of outer edges 3 and 8.

| Boundary | Live GL04 `{1zs}` Laurent-0 / `T2` | Passing standalone/emulated CFF test | Same input? |
| --- | --- | --- | --- |
| 4D coefficient | With `q=Q5=K1`, `p=Q3=-K2-P0+K0`, `D5=D(q)` and `D6=D(-q)`, the common-denominator term is `-q5_fixed_0^2 * [-D5*D6^2 +(mUV^2+p^2)*D5*D6 -4*(p.(-q))^2*D5 +mUV^2*D6^2 +2*(p.(-q))*D5*D6] /(D5^2*D6^3)`, up to the graph's common scalar/`i` prefactor. | The initial comparison uses `Q1_0^2*Q5_0^2 /(D1*D2*D3*D5*D6^2)`. Later subchecks use a single `D5` numerator cancellation, odd UV components, or `Q5_0^2*Q1_0*O - Q5_0*Q5_1*Q1_1*O`, with `O=prod(e=1,2,3)(Qe_0+1)`. No one subcheck contains the five-term GL04 Taylor polynomial. | **No: this is the first actual difference.** |
| UV denominator occurrences | Five serial occurrences: owner 5 has energy IDs `9,12`; owner 6 has `10,11,13`. In the common child LMB they are `D(+q)^2 D(-q)^3=D(q)^5`. | Three UV occurrences: one owned by edge 5 and two by edge 6, alongside three cograph occurrences. The first combined call has exact IDs `7,8,9` for the cograph component and `10,11,12` for the UV component. | No. |
| Full topology | Three-loop GL04. The child is the one-loop `e5/e6` bubble; after shrinking it, the two-loop outer graph retains physical edges `1,2,3,4,7,8`. Exact child residues and the outer physical residue are generated sequentially. | A two-loop graph: one outer triangle `e1/e2/e3` and one UV bubble `e5/e6` sharing a vertex. The primary check sends all six denominator occurrences through one standalone exact-CFF call, then compares with a product of two separately generated components. | No. |
| Factorized numerator provenance | The edge-local production factor is `i*(Q5^0)^2` (the common scalar/`i` factor is irrelevant to the route ratio). Both original `Q5^0` factors are `TaylorFixed` and stay on base occurrence 9. Denominator derivatives create tagged factors only: owner 5 may use new occurrence 12; owner 6 may use new occurrences 11 and 13. The soft `p=Q3` remains external to the child CFF and is mapped only by the outer CFF. | Numerators are hand-written, untagged physical `Q1`/`Q5` factors. Denominator owner IDs are retained, but no Taylor operation creates the live mixture of fixed and `DenominatorDerived` provenance. | No. |
| Child energy assignment | Exact bound plan `[(9,2),(11,2),(12,2),(13,2)]`. Its common-LMB substitution is symbolically identical to the pre-dispatch Taylor numerator. | For the first combined quadratic test the plan is physical `[(1,2),(5,2)]` -> exact `[(7,2),(12,2)]`. The later GL0-like cross-component source requests physical `[(1,2),(2,1),(3,1),(5,2)]` -> exact `[(7,2),(8,1),(9,1),(12,2)]`. Each physical owner has only one eligible exact occurrence in those plans. | No. |
| Child CFF context and maps | `GraphThreeDSource::from_exact_denominators_in_uv_sub_lmb(..., TaylorVacuum)` in the child sub-LMB, empty child `CutSet`, `CffGenerationContext::EmbeddedCffFactor`. The child mapper owns five occurrence maps; owner-5 raw `+q` has parsed sign `+1`, owner-6 raw `-q` has parsed sign `-1`, and all map to the same parsed `+q` carrier. | `cff_from_4d_denominators_in_uv_edges` uses the full graph LMB and the default `Standalone` context. Its combined source maps the cograph and UV contours in one call; it does not exercise a Taylor-vacuum child map followed by an outer map. | No. |
| Outer numerator/capacity | The graph numerator is localized on edge 5, so the untouched cograph numerator `resnum` is exactly one. Dependence on `p=Q3` is carried by the already factorized child coefficient. The outer generation nevertheless correctly requests `physical_parent_bounds=[(3,2)]` and `assigned_cff_source_bounds=[(3,2)]`. | The test deliberately adds `Q1_0^2` or the cubic `O` to the cograph. Its representative bounds are therefore generated by a genuinely non-unit outer numerator rather than by a child coefficient depending on a soft outer momentum. | No. |
| Physical CutSet | All live GL04 cuts are ordinary: `[1,2]`, `[3,4]`, `[1,7,8]`, each with `max_occurrence=1`. The simplest nonzero discrepancy is `[3,4]`; `[1,7,8]` has the larger signal. | The selected outer triangle surface is explicitly raised with `max_occurrence=2`, yielding `lu_cut_order=1,2`. The test often compares each order's value and first radial derivative separately. | No. |
| Runtime parameters and derivative assembly | Production `ParamBuilder` supplies the full model/external/loop values (`mUV=mUVexp=20`, `mu_r^2=9`, `M=1`, etc.). Since `max_occurrence=1`, `CrossSectionGraphTerm::evaluate` uses a nondual integrand value, the ordinary `h(t*)`, and `eta'(t*)`; only `CutCFFIndex::lu_cut_order=1` enters `pass_two_evaluators[0]`. There is no raised-LU integrand derivative in this MRE. | A test-local Symbolica substitution sets its own masses and momenta, tunes the external energy so `eta(t=1)=0`, series-expands through the first radial derivative, and evaluates with an Arb evaluator whose only parameter is `pi`. It does not use production `ParamBuilder`, the LU `h` function, or `CrossSectionGraphTerm`'s pass-two assembly. | No. |

The exact structural GL04 test named above agrees through the end of the
"child energy assignment" row: it proves the 4D coefficient, five-occurrence
topology, provenance dispatch, denominator diagonal, and common-LMB numerator
identity.  Its first missing boundary is CFF generation itself.  Conversely,
the standalone/emulated CFF test proves several useful generalized-CFF
factorizations, but differs already in its 4D coefficient.  Therefore its green
result cannot be used to infer that the complete live `T2` production input is
green, and its separate raised-LU representatives must not be compared with
the ordinary live cuts.
- 2026-09-04: Corrected the residue-order interpretation for the actual GL04
  production graph.  Its three physical Cutkosky groups, supported on
  `[1,2]`, `[3,4]`, and `[1,7,8]`, all have `max_occurrence = 1`.  The live
  discrepancy is therefore a Taylor-`T2` discrepancy inside an *ordinary*
  physical residue, not a raised-LU2 contraction discrepancy.  Statements
  that two individual LU orders or derivative-generated "moments" agree came
  from a synthetic raised-residue fixture and are not a production oracle.
  The EMR localization of separate derivative-generated pieces may differ
  between the direct-3D and projected-4D constructions; only the complete
  Taylor coefficient, summed over every generalized residue-map piece that
  belongs to one physical Cutkosky cut, is required to agree.  The current
  production comparison and all future regressions must enforce that summed
  contract.  The smallest live boundary is cut `[3,4]`: cut `[1,2]` agrees to
  Arb noise, while `[3,4]` retains a nonzero projected-minus-direct delta.
- 2026-09-04: Repeated the authoritative comparison on the exact physical
  `[3,4]` cut surface instead of the earlier generic `Q0=1` point.  With
  `Q0=E3+E4` imposed exactly, each of the twelve complete outer host groups
  agrees at about 300 decimal digits, and the full 48-branch direct sum agrees
  with the full 12-branch projected sum to relative `1.43e-300`.  The previous
  residual was off-surface and does not diagnose CFF or UV construction.  The
  first remaining numerical boundary is production preservation of the cut
  equation during precision escalation.

## Exact GL04 graph identity and edge namespaces

The persisted split-cut DOT and the in-memory production graph use different
surface spellings for the same initial-state cut.  `Graph::from_parsed` calls
`HedgeGraph::sew` before processing cuts: two dangling half-edges carrying the
same `is_cut` tag and opposite flows are replaced by one paired edge.  The
merge closure retains the incoming (`Flow::Sink`) half-edge data, and
`CutProcessingResult::permute` subsequently moves the sewn initial-cut edges
to the front of the physical edge namespace.  The active parser snapshot with
two split initial-cut pairs likewise turns eleven input edge statements into
nine in-memory edges.  Thus a split input with incoming `id=0` and outgoing
`id=9` does **not** leave a physical edge 9: those statements become the one
physical cut edge 0.  The ordinary `debug_dot()` serialization shows that
already-sewn representation directly.

For an identical-input GL04 oracle, the least ambiguous fixture is therefore
the ordinary in-memory/debug-DOT form with exactly physical edges `0..8`:
edge 0 is the paired initial-state cut `v1 -> v0`, and edges 1 through 8 have
the production incidences and LMB routing listed above.  Multiply the existing
autogenerated scalar numerator of edge 5 by `(Q5^0)^2` (rather than replacing
its common `i` factor), set that edge's DOD to zero, form the genuine `e5/e6`
DOD-two child, and select the complete ordinary physical residue on cut
`[3,4]`.  Do not force a DOD and do not construct a raised-LU cut.

Exact-source occurrence IDs cannot collide with a physical edge.  Their
canonical namespace starts at `graph.underlying.n_edges()`; with the correct
nine-edge graph, the five Taylor-`T2` denominator occurrences are `9..13`.
Had an extra physical edge really survived, the allocator would instead start
at 10, so it would produce `10..14`, not collide with physical edge 9.  Such a
graph would nevertheless be a different topology and every hard-coded
occurrence expectation would be shifted.  Oracle code should derive these IDs
from the source plan/energy map, using `9..13` only as an assertion of the
known production graph identity.
