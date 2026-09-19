# Vakint

Vakint matches, canonicalizes, tensor-reduces, and evaluates single-scale vacuum integrals for
high-energy physics. It combines Symbolica expression handling with analytic FORM-backed methods
and numerical pySecDec evaluation.

Use the canonical Typst documentation for maintained workflows:

- [overview and backend boundaries](https://alphal00p.github.io/gammaloop/products/vakint/latest/)
- [matching-only tutorial](https://alphal00p.github.io/gammaloop/products/vakint/latest/tutorial/)
- [staged evaluation, numerical substitution, result interpretation, and citations](https://alphal00p.github.io/gammaloop/products/vakint/latest/guides/evaluation/)
- [Rust and Python interfaces](https://alphal00p.github.io/gammaloop/products/vakint/latest/reference/interfaces/)

The crate also contains source examples under [`examples/`](examples/); check them against the
versioned API reference before adapting them. Vakint is exposed to Python through
`symbolica.community.vakint` when that community module is included in the installed Symbolica
assembly, not as an independent Python distribution.

For supported one-, two- and three-loop equal-mass vacuum families, the opt-in
`EvaluationOrder::rustred_only()` scalar backend instead applies closing IBP
artifacts shipped with Vakint and reuses Vakint's pure-Rust master
substitutions. This scalar tail neither invokes nor falls back to FORM. The
default evaluation order is unchanged; tensor-bearing inputs continue to use
Vakint's existing tensor prepass before scalar evaluation.

The same opt-in backend also has a four-loop candidate-program lane for the
H, FG, BMW and X parents and their matched contractions. Its
[bundled programs and offline master catalogs](data/rustred/four_loop/README.md)
are loaded lazily once; ordinary evaluation never regenerates IBPs or invokes
FORM. RustRed applies the rules and checks applicability, descent and known
terminal leaves at each requested integer target. These programs are explicitly
**not** unrestricted closing artifacts: unresolved targets produce an error,
and passing a numerical test is not a proof of arbitrary-rank coverage. The
four-loop public-backend gates passed on 2026-09-19: all fifteen original
numerical references and sixteen expanded-numerator/pinch comparisons against
FMFT, with invalid FORM paths for both FeynKit and the RustRed scalar tail.
The [test instructions and evidence](tests/experimental_rustred_4l/README.md)
record the original precision/tolerances and corrected offline catalog audit.
The measured 34.42 s and 104.55 s are full comparison-process times, not
RustRed scalar-only performance. This finite acceptance does not establish
unrestricted family closure or arbitrary-rank coverage.

The embedded artifacts and Vakint's pinned RustRed revision form one atomic
build-time contract and must use RustRed's current artifact schema. Vakint does
not migrate or decode obsolete RustRed artifact schemas. This is separate from
Vakint's backward-compatibility contract for its public API, defaults, and
existing FORM-backed evaluation methods.

The one-, two- and three-loop certified programs now use the same native binary
envelope as the four-loop candidate programs, with separate certified/candidate
payload kinds. Certified schema V6 stores Symbolica's native coefficients and
shared symbol state and still replays the proof on cold loading. The embedded
files are trusted generated data, not an arbitrary-file import interface.
The K1/K3/K6 migration re-encoded the existing saved rules without a new search;
the exact families, guards, source traces and terminal sets remain unchanged.
Their native files are 2,790, 36,692 and 3,725,739 bytes respectively. That migration
did not change terminal catalogs or the scalar application algorithm.

The four-loop terminal-value catalogs now also use RustRed's native Atom/state
codec, under a separate value-only envelope kind. All 1,155 exact PR expressions,
typed keys, family identities and coverage declarations were compared unchanged
in fresh processes; dictionary encoding reduces their total from 83,020 to
30,307 bytes. This is storage deduplication, not a claim of fewer physical masters
or stronger closure. Vakint only steers the RustRed loader and retains the exact
raw terminal-set check. No FORM invocation or IBP generation occurs during
catalog loading. See the [catalog migration details](data/rustred/four_loop/README.md).
On RustRed `d51721b6`, the unchanged 83-test selection through three loops and
all fifteen four-loop references plus sixteen expanded-numerator/pinch pairs
pass again, with invalid FORM paths in the native lanes. Seven catalog/loader
checks and three fixture checks also pass. That I/O-only milestone did not
enable terminal aliases.

The subsequent `f91c47ab` pin enables RustRed's verified terminal-routing plan
once per loaded four-loop parent. The 1,155 raw catalog declarations remain
unchanged; the applier coalesces them onto 505 family-local representatives
before memoizing coefficients. These are not a minimal master basis or a new
closure claim. Vakint delegates preparation and application entirely to
RustRed; its options, defaults and existing FORM-backed modes stay unchanged.
All 83 lower-loop cases, all 31 four-loop comparisons and nine focused checks
pass on this pin. In the matched nine-input benchmark, first H/X cubed-line
calls decrease from 14.86/84.17 s to 7.16/30.99 s; warm calls remain faster than
FMFT on the tested matrix, but initial cubed-parent reductions are still slower than FMFT
and preparation is not free. See the [complete timings and boundaries](data/rustred/four_loop/README.md#public-scalar-timing-probe).

The [historical RustRed acceptance inventory](tests/RUSTRED_ACCEPTANCE.md) maps
all existing single-common-mass inputs through three loops to native peers.
Its unchanged 83-test selection passes again after the V6 migration on
2026-09-19, with Symbolica 3 and RustRed `d6718733`, including the offline MATAD
check of all 38 K6 terminals. Native FeynKit/RustRed stages use an invalid FORM
path; separate legacy oracles use FORM. The full 16-test K6 pipeline, 26 focused
native library checks and three four-loop fixture checks also pass; these
overlap the selection rather than adding distinct numerical inputs. The
historical inventory distinguishes the 83-test matrix from the 40-entry/46-input
legacy census; the four-loop acceptance results are recorded separately above.

The migration gate used release builds and the actual public backend, with no
regeneration of shipped artifacts or changed tolerances. Both RustRed dependencies are pinned to
the published revision together with the embedded native assets. Workspace
dependency verification and CI metadata regeneration pass; this focused gate
does not claim a complete GammaLoop workspace CI run.

The shared parent-routing validator is also tested against all 19 registered
four-loop classes and their 123 surviving propagator slots. It reuses the
matcher's signed simultaneous witness, without graph rematching or topology-name
dispatch. The new four-loop candidate-program lane uses this same witness.
The earlier routing-only follow-up passes 39
focused tests, including the existing 16-test K6 pipeline; it is not a fresh
rerun of the complete 83-test selection or a four-loop numerical-parity claim.

`run_time_decimal_precision` controls arithmetic precision, not the accuracy of
stored master data. If a required numerical master or period constant carries
less precision, Vakint emits a `WARNING` through its usual Rust logger and
continues at the requested working precision. Configure the application's
logger to display warnings (for example, `RUST_LOG=warn` with an env logger).
Resizing a stored number does not supply additional accurate digits; in
particular, requesting 20,000 digits does not establish a 20,000-digit result.
FMFT checks constants surviving Laurent truncation; direct MATAD/AlphaLoop
substitution and numerical RustRed catalogs conservatively check each active
source record before resizing or materialization. Unsupported Laurent orders
and missing master data remain errors.

```rust
use vakint::{EvaluationOrder, Vakint, VakintSettings, vakint_parse};

let vakint = Vakint::new()?;
let settings = VakintSettings {
    evaluation_order: EvaluationOrder::rustred_only(),
    ..VakintSettings::default()
};
let result = vakint.evaluate_integral(
    &settings,
    vakint_parse!("topo(I2L(1,2,1,1))")?.as_view(),
)?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

## Python usage

`Vakint` is also part of the [`symbolica-community`](https://github.com/benruijl/symbolica-community) `Python` module where vaccuum graphs can be evaluated directly using the `Python` API exposed in that module.

Numerator tensor reduction has two backends. `"feynkit"` is the default: it uses
the native FeynKit projector, supports ranks through 20, and does not need FORM
for the tensor-reduction step. `"alphaloop"` explicitly selects the historical
FORM implementation, whose bundled projector tables cover ranks through 10:

```python
from symbolica import E
from symbolica.community.vakint import Vakint

vakint = Vakint()
integral = E(
    "k(1,mu)*k(1,nu)*p(1,mu)*p(1,nu)"
    "*topo(prop(1,edge(1,1),k(1),muvsq,1))",
    default_namespace="vakint",
)
reduced = vakint.tensor_reduce(integral)
```

The Rust API uses the same default through `VakintSettings`. Set
`TensorReductionMethod::AlphaLoop` explicitly to request the legacy FORM
projector:

```rust
use vakint::{TensorReductionMethod, VakintSettings};

let settings = VakintSettings {
    tensor_reduction_method: TensorReductionMethod::AlphaLoop,
    ..VakintSettings::default()
};
```

## Rust usage

Checkout the examples directory for a few examples of how to use this library.
For example:

```rust
use vakint::{vakint_parse, Vakint, VakintExpression, VakintSettings};

fn main() {
    let settings = VakintSettings {
        allow_unknown_integrals: false,
        ..VakintSettings::default()
    };
    let vakint = Vakint::new().unwrap();
    vakint.validate_settings(&settings).unwrap();

    //println!("Supported topologies:\n{}", vakint.topologies);

    let input = vakint_parse!(
        "(k(11,2)*k(11,2)+k(11,77)*k(22,77)+k(22,33)*p(42,33))*topo(\
                    prop(9,edge(7,10),k(11),mUVsqA,1)*\
                    prop(33,edge(7,10),k(22),mUVsqA,2)*\
                    prop(55,edge(7,10),k(11)+k(22),mUVsqA,1)\
                )+\
                (2*k(33,2)*k(33,2)+17*k(44,33)*p(17,33))*topo(\
                    prop(7,edge(9,21),k(33),mUVsqB,1)*\
                    prop(13,edge(9,21),k(44),mUVsqB,2)*\
                    prop(17,edge(9,21),k(33)+k(44),mUVsqB,1)\
                )"
    )
    .unwrap();

    let output = vakint
        .to_canonical(&settings, input.as_view(), true)
        .unwrap();

    println!(
        "\nInput:\n\n{}\n\nhas been matched to\n\n{}\n",
        VakintExpression::try_from(input).unwrap(),
        VakintExpression::try_from(output).unwrap()
    );
}
```

yields:

![result_matching](result_readme.png)

When carrying a complete evaluation, the workflow is as follows:

```rust
use ahash::HashMap;
use vakint::{vakint_parse, Vakint, VakintExpression, VakintSettings};

fn main() {
    let settings = VakintSettings {
        allow_unknown_integrals: false,
        use_dot_product_notation: true,
        integral_normalization_factor: vakint::LoopNormalizationFactor::Custom("1".into()),
        run_time_decimal_precision: 16,
        ..VakintSettings::default()
    };
    let vakint = Vakint::new().unwrap();
    vakint.validate_settings(&settings).unwrap();

    let mut integral = vakint_parse!(
        "(k(3,11)*k(3,22)+k(3,77)*p(8,77))*topo(\
        prop(9,edge(66,66),k(3),MUVsq,1)\
    )"
    )
    .unwrap();
    println!(
        "\nInput integral:\n{}\n",
        VakintExpression::try_from(integral.clone()).unwrap()
    );

    integral = vakint
        .to_canonical(&settings, integral.as_view(), true)
        .unwrap();
    println!(
        "Matched integral:\n{}\n",
        VakintExpression::try_from(integral.clone()).unwrap()
    );

    integral = vakint.tensor_reduce(&settings, integral.as_view()).unwrap();
    println!(
        "Tensor reduced integral:\n{}\n",
        VakintExpression::try_from(integral.clone()).unwrap()
    );

    integral = vakint
        .evaluate_integral(&settings, integral.as_view())
        .unwrap();
    println!("Evaluated integral:\n{}\n", integral);

    let mut params = HashMap::default();
    params.insert("MUVsq".into(), settings.real_to_prec("1.0"));
    params.insert("mursq".into(), settings.real_to_prec("1.0"));

    let numerical_partial_eval = Vakint::partial_numerical_evaluation(
        &settings,
        integral.as_view(),
        &params,
        &HashMap::default(),
        None,
    );
    println!("Partial eval:\n{}\n", numerical_partial_eval);

    params.insert("g(11,22)".into(), settings.real_to_prec("1.0"));
    let numerical_full_eval = Vakint::full_numerical_evaluation_without_error(
        &settings,
        integral.as_view(),
        &params,
        &HashMap::default(),
        None,
    )
    .unwrap();
    println!(
        "Full eval (metric substituted with 1):\n{}\n",
        numerical_full_eval
    );
}
```

yielding:

![result_evaluation](result_2_readme.png)

## Symbolica license

This library uses `Symbolica` for some of its computations.
`Symbolica` is licensed under a commercial license, *however* its usage is free for Hobbyist and also free under single-core restriction for Academic use.
