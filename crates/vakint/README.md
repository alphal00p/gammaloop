# vakint

Library for the computation of (single-scale) vacuum integrals in High Energy Physics.

The code will perform the matching of user input onto known topologies, and then perform their computation either:

* analytically, if the topology is known, using tensor reduction and parametric integration by parts.
To this end, it uses a combination of [`Symbolica`](https://symbolica.io/) and [`FORM`](https://arxiv.org/abs/1203.6543) scripts, as well as [`MATAD`](https://arxiv.org/abs/hep-ph/0009029) and [`FMFT`](https://arxiv.org/abs/1707.01710) for the four-loop case.

* numerically, if the topology is not known, using the sector decomposition algorithm implemented in `pySecDec`[https://arxiv.org/pdf/1703.09692].

## Python usage

`Vakint` is also part of the [`symbolica-community`](https://github.com/benruijl/symbolica-community) `Python` module where vaccuum graphs can be evaluated directly using the `Python` API exposed in that module.

Tensor reduction validates Lorentz index domains before vacuum projection. It
passes universal loop-momentum tensor kernels to FORM while retaining graph
coefficients as factorized expressions. Analytic integration uses temporary
coefficient symbols and restores their values afterward; numerical integration
combines the scalar sectors into one estimator. Numerator and loop-normalization
poles in epsilon increase the required integral order. Positive epsilon prefactors never reduce the conservative
requested depth. Rational dependence and scalar log/exp/sin/cos series are
supported; retained Laurent coefficients must be epsilon-independent. Opaque
epsilon-dependent functions produce an error.

Callers with additional spin tensors declare their Lorentz slots using
`vakint::tensor(body, slot, ...)`. The loop-independent body stays opaque;
metrics rename complete slots without modifying spin or color indices.
Projection can introduce `vakint::tensor_index(family, loop, ordinal)` dummy
slots, which callers translate to their own index representation before
unwrapping tensor bodies. GammaLoop performs both conversions automatically
and completes any new d-dimensional Dirac contractions before scalar backend
evaluation truncates epsilon. Vakint does not evaluate opaque spin algebra.

`VakintSettings.project_onto_tensor_integrals` defaults to `true`. Setting it to
`false` sends each complete numerator through FORM tensor reduction and the
scalar backend, instead of using universal tensor kernels and coefficient
aliases. Distinct denominator topologies remain separate in both modes.
GammaLoop exposes this choice as
`global.generation.uv.project_integrated_uv_cts_onto_tensor_integrals`; both
choices share its d-dimensional algebra and gamma5 checks.

`PySecDecOptions.numerical_parameters` stores mass and numerator parameters as
complex values. Real external momentum components retain their signs. This
adapter requires real pole masses and real external momenta; complex numerator
coefficients are supported. Physical PySecDec integration tests are manual only:
they are ignored and excluded from automatic local/CI test selection. Missing
opt-in or missing PySecDec causes an explicit manual invocation to fail instead
of reporting an unevaluated test as passed. Pure Rust adapter tests remain in
the normal suite because they do not launch PySecDec.

Run the four fast one-loop numerical checks explicitly when changing that backend:

```bash
RUN_PYSECDEC_TESTS=1 cargo nextest run -p vakint \
  --cargo-profile dev-optim -P test_gammaloop --test integral_evaluation_pysecdec_tests \
  --run-ignored only --ignore-default-filter \
  -E 'test(/^test_integrate_1l_(simple|complex_and_signed_parameters|cross_product|cross_product_with_additional_symbols_numerator)$/)'
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
    ).unwrap();
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
