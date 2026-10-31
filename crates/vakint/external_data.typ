= External reduction data and master capabilities

Vakint ships its 1L–3L RustRed programs, master coefficients and graphs without
prepared caches. The separate #link("https://github.com/alphal00p/vakint-external-data")[external-data repository]
contains independently packaged loop data and optional complete native caches.
FORM reduction rules and their hard-coded master formulas remain in Vakint.
See #link("external_data_validation.typ")[the migration validation report] for
release identities, numerical evidence, completed checks and remaining gates.

== Explicit loading

The official-download API and options below are available on native Unix targets.
They are not exposed on Windows or WebAssembly targets.

```rust
use vakint::{OfficialExternalDataOptions, Vakint};
let mut engine = Vakint::new()?;
let report = engine.load_official_external_data(4, &OfficialExternalDataOptions {
    local_artifacts_path: Some("./data".into()),
    ..Default::default()
})?;
```

The matching Python method is:

```python
report = engine.load_official_external_data(
    4, version="latest", local_artifacts_path="./data", refresh=False
)
```

The loop argument is always a maximum: four loads 1L, 2L, 3L and 4L. All packs
come from one immutable release. The default path is the current working
directory, with managed installations below `.vakint-external-data/<version>`.
The report identifies the version, per-loop paths, available caches, and whether
each pack was local, reused, or downloaded. Prepared caches are loaded by default:
`OfficialExternalDataOptions.load_prepared_caches` is `true`, and the Python
keyword `load_prepared_caches` defaults to `True`. Set the Rust option to `false`
or pass `load_prepared_caches=False` in Python to load the same authenticated data
without native cache decoding or import. Checksums still cover every listed
file, including optional cache bytes; this is an activation policy, not a change
to the dataset. Use this option with `local_artifacts_path` for a local checkout.

`cache_available` reports the manifest declaration. Per-loop `cache_loaded`
reports actual automatic attachment, and the top-level `load_prepared_caches`
field records the requested policy. Availability may therefore be true while
loading is false. Explicit cache load methods remain available after data loading.

If the supplied directory contains the repository manifest, Vakint validates and
loads it directly. It never downloads into or modifies this checkout. `latest`
means the checkout's current content; an explicit version must match its declared
identity. A missing file or incorrect checksum is an error.
Loop-pack manifests without a valid root release manifest also identify damaged
local data and cause an error without downloading into that checkout. An
unrelated valid `manifest.json` is preserved and does not block a managed install;
empty loop-named directories alone do not identify an external-data checkout.

Managed installations reuse validated files without downloading or extracting
again. `refresh=True` checks the `data-latest` release alias; pinned versions remain
immutable. A failed refresh reports an error and leaves the active engine intact.
Checksums, relative paths and schemas are checked before activation. Native
program/cache bindings are checked when their native payloads are activated. Loading a new snapshot replaces the prior external selection;
old clones keep their existing immutable data and cache sessions.

`load_external_data(path, maximum)` loads an existing checkout directly. The
`load_rustred_artifacts(loop_count, path)` method accepts a local pack/IBP path
or a sealed 1L–3L rules file; bare programs retain the active master values and
discard the preceding program's cache attachment. The
explicit native cache load/import/download methods and the saved-application
path interface remain available. No ordinary evaluation initiates a download.

== Four-loop acceptance tests

The four-loop numerical acceptance gates and native asset tests call the official
loader with maximum loop count four, the compatible immutable release declared
by the vendored manifest, and `refresh=false`. This avoids a stale cached `latest`
marker while preserving offline reuse.
Tests that inspect the five-loop application or graphs explicitly request five.
On first use it downloads
missing packs into the test process's current directory; subsequent processes
reuse the validated installation without network access. Cargo normally runs a
test binary from its package directory. Set `VAKINT_TEST_EXTERNAL_DATA` to an
existing checkout to override that location.
Managed `.vakint-external-data` directories are ignored by Git.

The scientific gates retain their explicit execution flags and FORM requirements.
For example, this existing gate uses the RustRed backend and an archived
independent reference without invoking FORM:

```sh
cargo test -p vakint --test rustred_four_loop_tests \
  h_rank_four_matches_reference_without_form -- --ignored --nocapture
```

Ordinary tests using only the vendored snapshot continue to run offline.

== Master representation and numerical capability

Rust uses `RustRedEvaluationOptions { master_mode: RustRedMasterMode::Symbolic }`.
Python uses `VakintEvaluationMethod.new_rustred_method(master_mode="symbolic")`.
The other values are `retain` and `numerical`; numerical substitution is default.

Retain mode returns family-qualified master integrals. Symbolic mode expands
available exact Laurent coefficients and leaves numerical-only coefficients as
snapshot-qualified symbols. Evaluate those through the same engine's numerical
evaluation method so it can resolve their data. Switching to a different data
snapshot does not silently reinterpret old coefficient symbols.

Numerical mode evaluates exact formulas at the requested precision. Numerical-only
coefficients enforce their validated significant-digit limit. Required epsilon
powers are computed from actual reduced coefficients, including reduction poles,
products and normalization, after exact cancellation. Missing entries are not zero.
Errors name the integral, coefficient power, requested capability, available data,
and package version. This source-data check is not a final-result error estimate.

Records of `constants.json` may carry `symbolic_definition`, the defining
integral printed by Symbolica with full namespaces in the inert heads of
#link("masters.typ")[the masters documentation]. The prose `definition` field is
unchanged. The loader parses it, requires every variable to be bound and every
head to be admitted, and compares it with the definition reconstructed from the
record's `definition_metadata` when that exists; a different loop-momentum
routing of the same propagators, related by a linear change of loop momenta
with unit Jacobian, agrees. Unparseable, unbound or contradicting definitions
reject the pack. A recorded definition must be explicit: it may not reference
another named constant (no `vakint` constant, HPL, GPL or master symbol, and
`zeta`/`polylog` only as the complete definition of that standard function), and
every normalization factor is written out. The data definition is authoritative
for its symbol. Names of built-in constants, such as `cl2`, `cl4`, `s6`, `D5bar1`,
`D5bar2`, `D6bar1`, `D6bar2`, `T111bar4`, `BNbar4`, `D6`, `T1ep`, `T1ep2` or `S2`,
cannot be redefined by a record. Packs without the field remain valid:
the four `K5Moment`/`K5LogMoment` records and the `PR7D2`, `PR10D2`, `PR15D2`,
`PR7D2Eps1` and `PR10D2Eps1` periods use that reconstruction, and any other
record is returned without a definition, with the reason, on inspection.

Built-in identities are reserved. This version adds the built-in constants
`log(q)` for every positive rational `q ≠ 1`, `cl4`, `D5bar1`, `D5bar2`, `D6bar1`,
`D6bar2`, `T111bar4` and `BNbar4` (the last one is a convenience consistent with
`matad.rs`, not a basis generator), and the alternating Euler sum
`s6 = Σ_(m≥k≥1) (-1)^(m+k)/(m^5 k) = ζ(-5,-1) + ζ(6)` (diagonal included), evaluated
as `π^6/945 + GPLReal(0,0,0,0,-1,1,1)` at any precision. Three-loop and four-loop
formulas both use the built-in `s6`; the four-loop record `s6weak` holds the same
number and stays loadable. Formulas should use these names directly:
`validate_parameter` rejects them as runtime parameter names, and a
`constants.json` record with one of these identities, with another built-in
identity, or with a definition head (`Integral`, `Interval`, `LoopIntegral`,
`SeriesCoefficient`, `EuclideanDot`) rejects the pack as an invalid finite
constant identity. The loader also rejects the pack when a `symbolic_definition`
is unparseable, has a free symbol or an unadmitted head, is an inexact number,
references another named constant or contradicts the reconstruction from
`definition_metadata`, and when that metadata is inconsistent with the symbol
(a `K5Moment`/`K5LogMoment` record whose metadata states another power, or a
`PR7D2`/`PR10D2`/`PR15D2`/`PR7D2Eps1`/`PR10D2Eps1` record whose Taylor order,
dimension, family or powers disagree with its name and the manifest). Packs that satisfy
these rules, including `ca49f328` and later, load unchanged.

The constant card states the normalization of each definition class:

#table(
  columns: (1fr, 1.3fr, 2.4fr),
  [*Class*], [*Examples*], [*Normalisation*],
  [Elementary and polylogarithmic numbers],
  [`𝜋`, `γ`, `zeta(n)`, `log(q)`, `polylog(n,q)`],
  [Plain real number given by the convergent integral shown, first interval
   outermost; no loop measure and no further factor.],
  [Alternating Euler sum],
  [`s6` (and the external-data record `s6weak`, same number)],
  [`s6 = Σ_(m≥k≥1) (-1)^(m+k)/(m^5 k) = ζ(-5,-1) + ζ(6)`, diagonal included:
   the double integral `∫_0^1 ∫_0^1 log(x_1)^4/24 · 1/(1+x_2) · (1/(1+x_1) + x_2/(1-x_1 x_2)) dx_2 dx_1`;
   a plain real number, no further factor.],
  [Clausen values],
  [`cl2`, `cl4`, `S2`],
  [`Cl_n(π/3) = Im Li_n(exp(i π/3))`, with `Cl_2(π/3) = -∫_0^(π/3) log(2 sin(t/2)) dt`
   and `S2 = 4/(9√3) Cl_2(π/3)`; plain real numbers.],
  [HPL projections],
  [`HR(w, z)`, `HI(w, z)`],
  [Real or imaginary part of `H(w; z)` at `z = (1±i√3)/2` on the straight path
   `t = z s`, `1 > s_1 > … > 0`; letters `1/t`, `1/(1-t)`, `1/(1+t)`; trailing zeros
   become `log(z s)^m/m!` on the principal branch.],
  [GPL and polylogarithm projections],
  [`GPLReal`, `GPLImag`, `PolylogReal`, `PolylogImag`],
  [Real or imaginary part of `G(a_1,…,a_n; 1) = ∫_(1>t_1>…>t_n>0) ∏ dt_i/(t_i-a_i)`,
   resp. `Li_n(z) = z/(n-1)! ∫_0^1 (-log t)^(n-1)/(1-z t) dt`.],
  [Bessel moments],
  [`K5Moment(k)`, `K5LogMoment(k)`, `K5MomentTaylor(k, j)`],
  [`∫_0^∞ t^k K_0(t)^5 dt`; `3 ∫_0^∞ t^k K_0(t)^5 (log t + γ - log 2) dt`; the
   coefficient of `ε^j` of the Taylor expansion (not the derivative) of
   `(e^γ/2)^(3ε) ∫_0^∞ t^(k+3ε) K_ε(t)^5 dt`.],
  [MATAD masters],
  [`D6`, `T1ep`, `T1ep2`, MATAD `NumericalSource` leaves],
  [Euclidean measure `d^d k/π^(d/2)` per loop, `d = 4-2ε`, unit masses, no sign,
   factor `e^(γε)` per loop; the coefficient of `ε^n` of the Laurent expansion,
   not the derivative.],
  [MATAD bars],
  [`D5bar1`, `D5bar2`, `D5bar3`, `D6bar1`, `D6bar2`, `T111bar4`, `T111bar5`,
   `BNbar4`, `BNbar5`],
  [The coefficient of `ε^n` (not the derivative) of the `matad.rs` prefactor
   times `LoopIntegral(4-2ε, …)` (measure `d^d k/π^(d/2)` per loop, no sign):
   `(1-2ε)(1-ε) e^(3γε)` for D5 and D6; `(1-ε)(1-2ε)(1-3ε)(2-3ε) Γ(1+ε)^(-3)`
   for BN and `√3 (1-2ε)(ε-1) Γ(1+ε)^(-2)` for T111, with no `e^(γε)` factor
   because the explicit `e^(-Lγε)` of `matad.rs` cancels the per-loop `e^(γε)`
   of the MATAD master.],
  [Public four-loop masters],
  [`PR4fin`, `PR4dfin`, `PR11dep0`, `PR4ep1`, `PR4dep1`, `PR14ep2`, `PR14ep3`,
   FMFT `NumericalSource` leaves],
  [Minkowski measure `(i π^(d/2) e^(-γε))^(-1)` per loop written as
   `(-1)^(Σa) e^(4γε) LoopIntegral(4-2ε, ∏ D_E^(-a))`, unit masses; the
   coefficient of `ε^n` of the Laurent expansion, not the derivative.],
  [`D = 2` periods],
  [`PR7D2`, `PR10D2`, `PR15D2`],
  [`LoopIntegral(2, ∏ 1/(q_E²+1))`: measure `d^2 k/π` per loop, unit masses,
   positive (hat) sign convention without `(-1)^(Σa)`, no exponential factor.],
  [`D = 2-2δ` jets],
  [`PR7D2Eps1`, …, `PR11dD2Eps4`, `PR0D2Eps2`],
  [The coefficient of `δ^n` of the Taylor expansion (not the derivative) of
   `e^(Lγδ) LoopIntegral(2-2δ, ∏ 1/(q_E²+1))`: measure `d^D k/π^(D/2)` per
   loop, unit masses, positive (hat) sign convention, `e^(4γδ)` for four loops.],
  [No definition],
  [runtime inputs, pySecDec estimates, AlphaLoop literals, `PR3ep2`, `PR13ep1`],
  [`definition` is `None`; the value is the recorded source value and
   `definition_unavailable_reason` states why.],
)

For a definition built from `LoopIntegral` the statement is not looked up by
name: it is derived from the definition itself. It names the coefficient taken,
the measure and dimension, the number of loops, the propagators with their
masses and `Σa`, and every factor in front of the integral (sign, `e^(Lγε)`,
`Γ` factors, polynomials, added terms), and ends with "No other factor". Its
label states the convention those factors realise (public Vakint/FMFT master,
MATAD master, MATAD bar with the matching `matad.rs` prefactor, `D = 2` period
or `D = 2-2δ` jet), or "Loop integral with explicit factors" when none applies;
a statement therefore never contradicts its definition.

HPLs at `exp(±i*pi/3)` use Symbolica's arbitrary-precision numerical algorithms.
The same evaluator supplies FORM-backed HPL substitution. Finite constants without
an evaluatable exact definition retain explicit numerical precision; resizing a
stored decimal does not provide more accurate digits.

== Partial campaigns

Recorded coverage predicates describe the actual campaign scope. Summary maxima
and completed cache workloads do not enlarge that scope. Outside certified
coverage, Vakint warns and applies only rules whose mathematical conditions hold.
Unsupported integrals remain exact, family-qualified residuals. Corruption, unsafe
guards, invalid transports and operational-limit failures remain errors.

After cancellation, warnings list remaining undeclared masters and parent-topology
DOT paths. Numerical evaluation reports these masters again, together with declared
terminals lacking coefficients. The initial 5L release has 905 declared campaign
terminals and no computed 5L coefficients; it does not certify a complete basis.

== Maintaining the vendored snapshot

Run `python3 tools/vendor_external_data.py /path/to/vakint-external-data` from the
Vakint crate. The tool validates input hashes and copies all available 1L–3L
coefficients, graphs and programs, excluding caches. Native cache compatibility is
determined by codec and authenticated program identities, not producer revision
strings alone. Rebuild only caches invalidated by actual program/format changes.

Mathematica may be used to develop independent numerical reference values. It is
not a build or runtime dependency of Vakint or the official-data loader.
