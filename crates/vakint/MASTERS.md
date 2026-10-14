# Master-integral evaluation handoff

Snapshot: 13 October 2026. This Markdown handoff is explicitly requested by
Valentin, as an exception to the repository's usual Typst documentation policy.

## 1. Repositories and immediate objective

Use **`alphal00p/gammaloop`, branch `vakint_rustred`**, for Vakint/FeynKit, and
**`alphal00p/rustred`, branch `main`**, for the reduction engine. The partial
five-loop integration pins RustRed
`1fca5bd0ae11b8d25c89658a9c78a476a2b76048` in Cargo. Respect the checked-in lockfile:
this is the Symbolica 3.0.1 integration, with a shared native Symbolica kernel.

The next agent's numerical task is to extend the existing one-, two-, three-
and four-loop master data by one epsilon order, and obtain five-loop master
Laurent coefficients through the finite term, epsilon^0. **Inventory the existing
orders per master first: they are not uniform.** Higher orders may additionally
be needed when reduction coefficients contain poles; see section 7.

This branch is a tested partial implementation, not a declaration that every
five-loop application test passes. It has an external-package symbolic backend,
native DOT input, factor products and exact coefficients. Five-loop numerical
substitution is deliberately unavailable. Keep the new numerical work separate
from ongoing application/eligibility work; do not regenerate IBPs, modify a
running campaign, or reinterpret an unknown reduction leaf as a master.

## 2. Definitions, values and precision are different objects

A RustRed terminal is identified by **(full family fingerprint, integer powers)**.
A displayed family hash or a friendly master name is not its mathematical
definition. The family supplies ordered denominators, loop momenta, dimension,
scalar-product coordinates and power shifts. Zero powers remove a denominator;
negative powers put that *denominator polynomial*, including its mass term, in
the numerator. Never replace it silently by a homogeneous dot product.

There are three layers to inspect independently:

1. Generated rules and their declared terminal keys.
2. Exact terminal-to-terminal, symmetry, factorization or MATAD/FMFT-basis maps.
3. Laurent coefficients and their available numerical precision.

Many terminal keys can have the same value or project onto a much smaller
conventional basis. Neither 74 four-loop output keys nor 905 five-loop terminal
keys means that many independent numerical constants must be computed.

All paths in the next sections are relative to `crates/vakint/` unless marked
**RustRed**. Do not dump the complete numerical source files into a chat: some
lines contain tens of thousands of digits.

## 3. One through three loops

| What | Definition / mapping | Existing values and extensions |
|---|---|---|
| One-loop tadpole | [artifact.rs](src/rustred_evaluation/artifact.rs), K1 terminal `[1]` | Exact source `-Gam(1,1)/(ep*(ep-1))`; extend the gamma expansion, not a numerical IBP computation |
| Two-loop sunset and pinch | Same file, K3 `[1,1,1] -> -miT111`, `[0,1,1] -> (Gam(1,1)/(ep*(ep-1)))^2` | `miT111` in [matad.rs](src/matad.rs) and [matad_numerics.rs](src/matad_numerics.rs) |
| Three-loop fully massive family | [terminal/k6.rs](src/rustred_evaluation/terminal/k6.rs), 38 explicit terminal keys | Exact projections onto tadpole products, sunset-times-tadpole, `miBN`, `miD6`, and a five-line expression involving `miD5` |
| Broader legacy MATAD basis | [matad-ng.hh](form_src/matad/matad-ng.hh), including the line/mass definitions and integral-to-`mi*` identities | Nine three-loop names: `miD6`, `miD5`, `miD4`, `miD3`, `miDN`, `miDM`, `miE3`, `miBN`, `miBN1` |

The K6 denominator order is
`k1²-1, k2²-1, k3²-1, (k3-k1)²-1, (k1-k2)²-1, (k2-k3)²-1`.
Its 38 terminals form six sector orbits but five distinct mapped expressions;
that is not 38 independent masters. In particular, **retain the full five-line
mapping and its outer minus sign**; it is not simply `miD5`.

The shipped binaries are `data/rustred/unit_mass_vacuum_k{1,3,6}.rrbin`.
The integers K=1,3,6 count basis denominators, not loops.

The main value entry points are:

- `matad::MASTERS_EXPANSION`: symbolic epsilon expansions and explicit unknown tails.
- `matad::MASTERS_SUBSTITUTION`: analytic constants, including harmonic polylogarithms.
- `matad_numerics::DIRECT_SUBSTITUTIONS`: high-precision direct numerical series.
- `matad_numerics::{HPL_SUBSTITUTIONS, ADDITIONAL_CONSTANTS}`: finite stored period data.
- Symbolica's native `gamma`, `polygamma`, and `euler_gamma`: exact normalization and
  numerical evaluation at the requested precision, without fixed polygamma tables.
- [alphaloop_numerics.rs](src/alphaloop_numerics.rs): the alternative backend's data;
  retain its conventions and comparisons when extending shared constants.

Current explicit unknown-tail boundaries in the MATAD symbolic/direct tables:

| Names | First unknown power | Next coefficient to supply |
|---|---:|---:|
| `miT111` | 5 | epsilon^5 |
| `miD6`, `miD5`, `miD4`, `miD3`, `miDM`, `miDN` | 3 | epsilon^3 |
| `miE3` | 4 | epsilon^4 |
| `miBN`, `miBN1` | 5 | epsilon^5 |
| Generic gamma helper expansions | 6 | epsilon^6, or replace the finite helper with an appropriate native series service |

These boundaries describe the explicit table tails, not a proof that every
composite result is available through the preceding power. The direct MATAD
tables contain approximately 20,000-digit numerical coefficients, but that does
not extend their epsilon depth. Check native precision metadata and every
dependent constant rather than assuming all paths have the same accuracy.

## 4. Four loops: evaluate the PR basis, not 1,155 terminals independently

The four installed parents are H, X, BMW and FG. In
[`data/rustred/four_loop/`](data/rustred/four_loop/):

- `*.toml` and `*.csv` define the actual denominator/momentum ordering.
- `*.candidates.rrbin.gz` contain the saved rules and raw terminal sets.
- `*.rrnorm.bin` provide exact weighted terminal normalization.
- `*.rrcat.bin.gz` map the retained output keys **exactly onto FMFT PR symbols**.
  They are not new high-precision numerical evaluations.

| Parent | Raw terminals | Retained catalog keys |
|---|---:|---:|
| H | 386 | 22 |
| X | 445 | 19 |
| BMW | 179 | 17 |
| FG | 145 | 16 |
| Total, family-local | 1,155 | 74 |

Use [four_loop.rs](src/rustred_evaluation/four_loop.rs) for loading/binding and
[`OfflineTerminalCatalog`](src/rustred_evaluation/experimental/catalog.rs) for
the typed value map. Its `terms()` returns integral keys and native Symbolica
expressions. RustRed's underlying `ExactTerminalCatalog` uses native Atom/state
serialization. Decompress/read through existing code; do not invent a text
artifact decoder or manually edit the compressed catalogs.

The numerical source is [fmft_numerics.rs](src/fmft_numerics.rs):
`MASTERS_EXPANSION`, `MASTERS_NUMERIC_SUBSTITUTIONS`, and dependent constants.
The names include `PR0` through `PR15` and dotted/auxiliary variants `PR4d`,
`PR9d`, `PR9x`, `PR11d`. The conventional count "16 masters" is **not** a rule
that all these extra symbols are independent. Integral definitions and exact
relations are in [fmft.frm](form_src/fmft/fmft.frm); search for assignments to
`PR*` and `id PR9x =`. The latter explicitly eliminates `PR9x` in favor of other
PR symbols. Keep physical momentum definitions and normalization with each name.

Current declared first-unknown orders:

- `PR1`, `PR2`: epsilon^3.
- `PR9d`, `PR11d`: epsilon^1.
- `PR9x`: epsilon^0 as a standalone table entry; use its exact relation first.
- The remaining `PR0`--`PR15` entries, and `PR4d`: epsilon^2.

Audit dependent tails: a coefficient symbol may itself contain `Oep(...)`.
For example `PR13ep1` retains a dependent `Oep(2,iD6)` sentinel. `PR3ep2` is
tabulated but not automatically exposed by the current `PR3` expansion row.
**Audit the existing `PR6` row before extending it:** it currently adds
`PR6ep1` without an explicit `ep` multiplier, although its numeric substitution
is a constant. This handoff records the discrepancy; it does not silently fix
or certify that table. Independent validation matters because both RustRed and
FMFT evaluation paths can reuse the same numerical tables.

Four-loop numeric source precision varies, roughly 26--50 digits for the
tabulated PR coefficients; some entries carry explicit "TO CHECK" comments.
Do not advertise these as 20,000-digit evaluations. See
[master_precision.rs](src/master_precision.rs): raising working precision emits
a warning when it exceeds a stored source's precision, but cannot add digits.

## 5. Five-loop definitions and external package

No five-loop numerical values are included in this partial implementation.
The fixed application package derives from THE_ONE **generation 54**, scope
`f7ee869008fccc6ed6cf2f4e3342c78d2b8607f59ac4d523dcae504595c1b0f8`:
67 owners, 9,982 rule records, 804 recorded requests, 939 terminal source keys,
and **905 remaining family-qualified keys** after exact finite substitutions.
These are a finite, potentially nonminimal evaluation basis, not an independence
claim. They are not the moving D<=10 campaign publication.

Local package directory on the shared machine:

```text
/shared/codex-2/rustred/TMP/five-loop-vakint-20261013/application-gen54-first
```

Size: 281,312,247 bytes. SHA256 of `application.json`:

```text
01425e1cecc71893d2e3a592a162600100326f285914926df61e43dd0c870e4b
```

This is the manifest checksum, not an archive checksum. Copy the **entire
directory**, including compressed native chunks and `terminal-maps.rrbin`.
It is portable and does not require the campaign traversal ledger. The user
will host it; **never commit/push this five-loop package or generated campaign
outputs**. If the receiving machine cannot access it, request the hosted copy.

From the RustRed checkout, export the master-key inventory without generating
rules or reducing any integral:

```bash
export RUSTRED_APPLICATION=/shared/codex-2/rustred/TMP/five-loop-vakint-20261013/application-gen54-first
mkdir -p TMP/master-evaluations
cargo run --release --locked -p rustred-app --features capacity-dispatch \
  --example saved_application -- terminals "$RUSTRED_APPLICATION" \
  > TMP/master-evaluations/terminals-gen54.json
jq '{count: (.terminals | length), families: (.families | length), minimality_claim}' \
  TMP/master-evaluations/terminals-gen54.json
```

Expected: 905 terminals, one family, `minimality_claim: false`. This command
loads the native package; its initialization is not instantaneous. The JSON
family field is an index into the accompanying full fingerprint list, not a
standalone master number. Do not replace this inventory with the 813 application
test witnesses or the 804 query regions; neither is a master list.

Use RustRed's `load_saved_application::<15>` and `PreparedSavedApplication`
accessors `remaining_terminals()`, `terminal_sources()`, `output_family()` and
`output_family_shared()` for typed access. Inspect the returned `IntegralFamily`
through `loop_momenta()`, `dimension()`, `denominators()`, `coordinates()` and
`power_shifts()`. A displayed hash alone cannot reconstruct these objects.

The readable common geometry is **RustRed**
`examples/input/tide_five_loop.toml`, with its manifest and exact routing-witness
JSON files next to it. All fifteen coordinates are shifted forms `q_i²-1`;
`D14=(k1+k2-k4)²-1` is intentional. Auxiliary coordinates are not a contiguous
trailing block, and the 15-coordinate basis is not a graph with 15 physical
propagators. Always bind a value to its actual family plus powers.

For the exact 939-to-905 substitution rows, the source publication is locally:

```text
/shared/codex-2/rustred/campaigns/five-loop-rank-ladder-THE-ONE/master-reduction/scopes/f7ee869008fccc6ed6cf2f4e3342c78d2b8607f59ac4d523dcae504595c1b0f8
```

Read it with RustRed's `load_master_reduction`, its `collection()` and
`apply_terminal(family, key)` interfaces. Do not parse `terminal-maps.rrbin`
as text. The runtime already applies these substitutions; they must not be
applied twice to its outputs.

Before scheduling expensive five-loop numerical work, use the existing native
terminal-product and product-chart services to identify exact lower-loop
products. `prepare_terminal_products` / `terminal_product` cover verified
scalar factorization; numerator-bearing products use product charts and the
existing FeynKit tensor projector. Vakint's `ShippedFactorProvider` and
`VacuumProductReducer` then reduce the factors. Do not infer factorization from
names, and do not assume the 19 factorized topology classes equal 19 masters.
Reuse known lower-loop values, with their required extra epsilon orders.

The `terminals` command lists the flat published 905-key convention. Vakint can
instead emit products containing lower-family `RustRedMaster` symbols. Resolve
those through the returned family map and the native factor/shipped-family
registry; they are not extra rows secretly missing from that flat inventory.
The existing 52 factor-product controls are not an independence count.

## 6. Fix the integration convention before computing numbers

Use `d=4-2*epsilon`, unit mass squared one, and the **stored denominator signs**.
The new five-loop `RustRedMaster(family_hash, powers...)` represents a unit-mass
Minkowski integral with the unnormalized loop measure, not a pre-normalized
MATAD number. [five_loop.rs](src/rustred_evaluation/five_loop.rs) returns a family
map alongside these symbols; the hash is SHA256 of the full family fingerprint.

The adapter restores `(m²)^(L*d/2-sum(powers))` and the configured loop
normalization once, including factorized answers. Extra numerator mass powers
are included by the weighted lowering. Do not apply the original five-loop
mass/measure factor a second time to each lower-loop factor.

The older one-through-four-loop value tables use the existing pre-normalization
MATAD/FMFT conventions. In particular,
`Gam(a,b)=exp(b*epsilon*EulerGamma)*Gamma(a+b*epsilon)` in
`MATAD::substitute_gam_functions`. The finalizers
`MATAD::finalize_reduced_masters` and `FMFT::finalize_master_expression`, plus
`LoopNormalizationFactor` in [lib.rs](src/lib.rs), are authoritative for conversion.
The default `pySecDec` factor is `(i*pi^(d/2))^(-L)`; `FMFTandMATAD` uses
`(i*pi^(d/2)*exp(-EulerGamma*epsilon))^(-L)`. In particular, the exponential
inside that parenthesis is also inverted. Custom and MSbar settings differ.
Record which quantity each new table evaluates.

Euclidean `q²+1`, Minkowski `q²-1`, TIDE's normalization by powers of its
tadpole J, and AMFlow's chosen measure are not interchangeable. Establish the
conversion on the exact one-loop tadpole, sunset, an irreducible three-loop
anchor and a factorized product before importing a five-loop number. Carry
the sign/phase, mass dimension, epsilon convention and normalization as data.

## 7. Epsilon depth is not determined by loop count alone

For `c(epsilon) * product_i M_i(epsilon)`, let `v_c` be the coefficient's
lowest epsilon power and `ell_j` the other masters' leading powers. To know
the product through epsilon^T, a safe required order for factor i is

```text
K_i >= T - v_c - sum_{j != i} ell_j.
```

Thus a coefficient proportional to `1/epsilon^p` generally needs its master
through epsilon^p to obtain a finite answer. Five-loop masters through epsilon^0
are a valuable first target, **not a guarantee that every reduced finite part
can already be evaluated**. Coalesce exact coefficients before measuring pole
orders, and exploit cancellations only when proved.

Even with regular coefficients, worst-case leading poles give requirements
such as one-loop through epsilon^4 times a four-loop factor through epsilon^1,
or two-loop through epsilon^3 times three-loop through epsilon^2. Actual leading
orders can be less singular. Do not multiply finite-only factor tables and
assume their finite coefficient is determined.

Existing finalizers use
`last_output_power = number_of_terms_in_epsilon_expansion - L - 1`.
Increasing that setting requests more output; it does not supply absent table
coefficients. RustRed's current numerical dispatch also has a through-epsilon^0
four-loop-era request ceiling (`number_of_terms_in_epsilon_expansion <= 5`);
audit/extend it only alongside real new data and tail tests. Preserve every
`Oep(...)` sentinel and unknown-order error until its coefficient is genuinely
known. Never substitute zero for an unavailable coefficient.

## 8. Suggested handoff workflow and acceptance

1. Export the pinned terminal inventory and record family definitions,
   normalization, exact relations and existing leading/last-known orders.
2. Reuse factor products and the exact K6/PR maps. Extend existing independent
   constants rather than evaluating every redundant terminal separately.
3. Compute new data using the approved master-evaluation toolchain. Store exact
   expressions where available, otherwise high-precision values with source,
   precision/error evidence and epsilon order. Working precision alone is not
   an accuracy certificate.
4. Compare old orders against genuinely independent references. Two reduction
   modes sharing the same value table are not independent numerical validation.
5. Extend the pure-Rust substitution/finalization services. Keep FORM available
   as an offline legacy oracle, never as a new runtime fallback.
6. Test nonunit mass, factor products, spurious coefficient poles, insufficient
   order/precision diagnostics, and all existing one-through-four-loop numerical
   regressions. New five-loop values must cover exactly their declared typed keys.

Useful tests: `rustred_scalar_evaluation_tests`, `k6_matad_oracle_tests`,
`rustred_four_loop_tests`, `experimental_rustred_4l_tests`,
`integral_alphaloop_vs_matad_tests`, and the `master_precision` unit tests.
The symbolic five-loop examples are
[`five_loop_saved_application.rs`](examples/five_loop_saved_application.rs) and
[`five_loop_saved_application.py`](examples/five_loop_saved_application.py).
Read [five_loop_application.typ](five_loop_application.typ) for API details,
known incomplete/resource-limited validations and the external artifact setting.

This partial push does not claim all 813 application witnesses have reduced,
nor numerical five-loop acceptance. The physical QCD ladder is also a separate
ongoing eligibility test, not a replacement for the master inventory. Coordinate
source edits with that lane, use separate workspaces, and never mutate the live
THE_ONE campaign to begin the numerical master task.

The latest native unit run passed 119 tests (seven external/manual tests ignored),
and the loaded bounded-cache and isolated Python DOT checks passed. Four existing
integration failures were independently reproduced unchanged at baseline
`6203c6cbba6ae5e90329ba5081fad55319e678db`: two tensor-presentation assertions,
the legacy repeated-index-power scalar control, and the nonunit three-loop
momentum-head declaration conflict. Their assertions remain intact; this partial
push does not claim a completely green legacy integration suite.
