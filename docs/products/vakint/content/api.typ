#import "../../shared.typ": boundary, catalog-contract, source-link

#let api = [
= Rust and Python APIs

#catalog-contract(
  rust-scope: "vakint",
  python-scope: "symbolica.community.hepkit.vakint",
)

== Rust package

The main Rust types make the workflow state explicit:

- `Vakint` owns the topology library and matching/evaluation operations;
- `VakintSettings` owns symbols, precision, normalization, tool paths, validation, and backend
  order;
- `VakintExpression` and `VakintTerm` separate numerators from topology atoms;
- `EvaluationOrder` and `EvaluationMethod` select RustRed, AlphaLoop, MATAD, FMFT, or pySecDec;
- `NumericalEvaluationResult` stores Laurent coefficients in the dimensional regulator.

For a dependency-free canonicalization setup, select an empty backend order before validating:

```rust
use vakint::{EvaluationOrder, Vakint, VakintSettings};

let engine = Vakint::new()?;
let settings = VakintSettings {
    evaluation_order: EvaluationOrder::empty(),
    ..VakintSettings::default()
};
engine.validate_settings(&settings)?;
```

`vakint_parse!` and the related helpers parse Symbolica expressions in Vakint's namespace. Handle
parse failures before passing a user-supplied expression to matching or evaluation.

== External RustRed data

`Vakint::load_official_external_data(max_loop_count, &options)` installs and loads
all loop packs from one through the requested maximum. `OfficialExternalDataOptions`
defaults to version `latest`, the current working directory, and no refresh.
The matching Python method accepts `max_loop_count`, `version`,
`local_artifacts_path`, and `refresh`. A path containing a checked-out external-data
manifest is validated and loaded without network access or checkout modification.
Installed releases are reused offline; explicit refresh checks for updates.

`RustRedMasterMode` chooses `Retain`, `Symbolic`, or `Numerical` (default). Python
selects the same behavior with `new_rustred_method(master_mode="symbolic")`.
Symbolic mode leaves numerical-only Laurent coefficients as named symbols, which
the same engine can resolve later. Numerical mode reports finite-source precision limits,
and rejects missing epsilon coefficients or unresolved masters with their exact identities.
All three RustRed modes preserve inexact input coefficients for later tracked evaluation.
Passing retained uncertainty sources to pySecDec returns an explicit error before input generation.
Use #source-link("crates/vakint/external_data.typ", label: "the external-data guide")
for campaign coverage, cache bindings, and release conventions.

== Incremental five-loop caches and monitoring

`attach_running_cache(5, path, ...)` creates or resumes a standalone native delta;
`load_running_cache` requires an existing one. Native capture keeps strictly
completed roots and descendants with their exact nonzero conditions, including
useful completed work from a parent that ultimately remains partially reduced.
Partial nodes never become strict shortcuts. Immutable bases stay separate from
the delta; use `read_backing=False` for capture without persistent reads, or
`read_delta=False` to keep reads restricted to explicit immutable bases/snapshots.

The returned `RunningCacheControl` exposes `report()`, `flush()`, `snapshot(path)`
and `close()` independently of the reducer lock. Accepted records are queued;
committed records are durably acknowledged. `save_running_cache(5)` flushes;
`freeze_running_cache(5, path)` writes a new native index outside the delta and
returns its descriptor.
`load_cache_snapshot(5, snapshots=[descriptor])` reads a frozen generation without
capturing new work. Retain the descriptor, referenced data and explicit bases;
freezing does not copy payloads. Closing the control releases its writer lease.
Bundle-managed child directories cannot be used as standalone writable deltas.

Python's `evaluate_rustred_symbolic(integral_expression, ...)` returns the retained
value, nonzero conditions and family map. This method and `prepare_cache` accept
`monitor="off"`, `"json"` or `"tui"`, a `progress_snapshot` path and
`progress_interval_seconds`. JSON mode requires a snapshot path.
Monitoring defaults off; the TUI requires
`progress-tui`. Pressing `q` closes only the display: reduction and capture continue.

On supported native platforms, `import_cache` and `download_cache` transport
directory packs for one through five loops. Select `format="multipart_manifest"`
explicitly for a large native pack and authenticate its manifest with `sha256`;
`max_manifest_bytes` optionally adjusts its read limit. Import/download publishes
a new directory without attaching it. There is no 50/95 MiB aggregate native-pack
cap. See #source-link("crates/vakint/cache_preparation.typ", label: "the cache guide")
for the lifecycle example, explicit operational limits and snapshot semantics.

The #link("reference/rust/")[Rust orientation] identifies the compiled crates, versions, and
feature matrix and links directly to native Rustdoc. The generated
#link("reference/topologies/")[topology reference] records the runtime topology library and
supported external-tool versions for this release.

== Python community module

#boundary("Python availability", [
  Python users import `symbolica.community.hepkit.vakint`; Vakint is not distributed as a standalone
  Python package. Check that the installed Symbolica distribution includes the Vakint community
  module before relying on this interface.
])

Install and verify the published assembly with
`python -m pip install --upgrade symbolica` and
`python -c "import symbolica.community.hepkit.vakint"`. Custom Symbolica builds can add Vakint to the
#link("https://github.com/symbolica-dev/symbolica-community")[symbolica-community] assembly, enable the
`symbolica_community_module` feature, and register `VakintWrapper` while building the extension.

The structured #link("reference/python/")[Python reference] covers the engine, expression,
evaluation-method and numerical-result classes, plus typed constant/master identities,
inspection records, precision metadata and warnings. This example constructs a matching-only
engine without probing external backends:

```python
from symbolica import E
from symbolica.community.hepkit.vakint import Vakint

vakint = Vakint(evaluation_order=[])
expr = E(
    "topo(prop(18,edge(7,7),k(99),muvsq,1))",
    default_namespace="vakint",
)
canonical = vakint.to_canonical(expr, short_form=True)
assert "I1L" in str(canonical)
```

`Vakint` exposes canonicalization, tensor reduction, integral evaluation, complete evaluation,
and numerical conversion. `VakintEvaluationMethod` constructs configured backend choices.
`VakintNumericalResult.coefficients` maps epsilon powers to read-only `NumericalCoefficient`
objects. Their `numerical_value` is a native Symbolica `ComplexFloat`, and `precision.real`
and `precision.imaginary` expose achieved significant digits, absolute error, exact-zero status
and error quality. Absolute errors and real constant/master-coefficient values use Symbolica's
existing arbitrary-precision `Float` class. `requested_decimal_precision` and `warnings` belong
to the complete result.
`numerical_evaluation` continues to return `(central_result, optional_statistical_error)`.

Parameter values and external four-vector components accept native `Float` and `ComplexFloat`
objects without changing their precision. Python `float` and `complex` inputs carry their
binary64 precision, even when more working digits are requested. Python integers and numerical
Symbolica expressions are substituted symbolically before tracked evaluation, so exact rational
masses and exact cancellations in momentum contractions are preserved. Supply high-precision
data as, for example, `Float("1.234567890123456789", decimal_digits=200)`, rather than first
converting it to a Python float.
Floating coefficients embedded in an `Expression` retain separate source identities for
different parameters. Equal central values therefore do not turn independent uncertain inputs
into a certified exact cancellation.
Parameter mappings cannot override registered constants or reserved data-source identities;
qualified ordinary parameter names remain supported.

`to_expression()` preserves independently achieved component precision and raises an error for
uncertain zeros or values without reliably encodable digits. `to_list()` intentionally rounds
to Python floats and drops accuracy metadata. `compare_to()` retains its relative threshold
and optional statistical-error arguments.

```python
from symbolica.community.hepkit.vakint import Vakint, ConstantSymbol

engine = Vakint(evaluation_order=[])
constant = engine.get_constant(ConstantSymbol("z3"), prec=40)
print(constant.available_precision.kind)
print(constant.precision.real.significant_digits)
print(constant.definition)  # 1/2*vakint::Integral((1-vakint::t)^-1*log(vakint::t)^2, vakint::Interval(vakint::t,0,1))
print(constant.normalisation)  # Plain real number: the convergent integral shown, …

master = engine.get_master(engine.master_symbols()[0], prec=32)
for power, coefficient in master.coefficients.items():
    print(power, coefficient.numerical_value, coefficient.symbolic_expression.kind)
```

`get_constant` and `get_master` query loaded or vendored RustRed data without downloading packs.
Both accept typed identifiers and default to 16 decimal digits. Discovery returns qualified
identifiers; ambiguous friendly names report the matching alternatives. `VakintConstant`,
`VakintMaster`, `MasterCoefficient`, and the shared precision/provenance records have structured
terminal representations and `to_html()` / `_repr_html_()` notebook views. Native Symbolica
rendering preserves long formulas and numbers behind expandable sections, and visible warnings
explain precision shortfalls.
Loaded finite constants such as `K5Moment(1)` report their validated source cap.
`constant.definition` is the defining multiple integral as an exact `Expression` in inert
`vakint::Integral`, `vakint::LoopIntegral` and `vakint::SeriesCoefficient` heads; it is never
evaluated and never references another named constant. Constants without a defining integral,
such as runtime input literals, keep their value: `definition` is `None` and
`definition_unavailable_reason` states why. `constant.normalisation` is a plain-language
statement of the conventions, also shown as "Normalisation" on the constant card:

#table(
  columns: (1fr, 1.3fr, 2.4fr),
  [*Class*], [*Examples*], [*Normalisation*],
  [Elementary and polylogarithmic numbers],
  [`𝜋`, `γ`, `zeta(n)`, `log(q)`, `polylog(n,q)`, `s6weak`],
  [Plain real number given by the convergent integral shown, first interval
   outermost; no loop measure and no further factor.],
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

Master formulas retain their original snapshot internally; notebook and terminal views show the
mathematical names, while the expression getter preserves the binding needed for evaluation.
Symbolic reductions preserve the same finite-period bindings when their formulas are evaluated
later; loading another data snapshot does not silently replace those sources.

`master.graph` is the normal FeynKit `FeynmanDiagram`, with its existing DOT, SVG and HTML
interfaces. In Marimo, display the master directly or use `mo.Html(master.to_html())`.
Graph previews require FeynKit's optional `typst-py` renderer; metadata and coefficient displays
remain available without it. Rendering does not reevaluate numbers or change the data snapshot.

Construct one engine and reuse it. When enabling pySecDec, provide numerical masses and external
momenta where the topology requires them. Reuse existing pySecDec output only while debugging,
and remove it before evaluating changed inputs.

Source starting points are #source-link("crates/vakint/src/lib.rs", label: "the Rust API") and
#source-link("crates/vakint/src/symbolica_community_module/mod.rs", label: "the Python module").
]
