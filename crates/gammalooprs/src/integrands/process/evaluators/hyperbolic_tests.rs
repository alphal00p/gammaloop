//! Numerical contracts for the production complex evaluator paths.
//!
//! These tests intentionally expose numerical failures in upstream functions and
//! composed derivative expressions. Runtime arguments prevent constant folding
//! from replacing the backend under test. MPC supplies an independent oracle;
//! agreement between GammaLoop backends alone is insufficient.
//! ArbPrec keeps a fixed working precision, unlike Symbolica's precision-tracking
//! Float. Derivative accuracy is reported separately: cancellation in a composed
//! expression can lose significant bits even when its built-in functions are accurate.
use std::collections::BTreeMap;

use rand::{Rng, SeedableRng, rngs::StdRng};
use rug::{Complex as Mpc, Float, float::Constant};
use spenso::algebra::complex::Complex;
use symbolica::{
    atom::{Atom, AtomCore},
    domains::float::Float as SymbolicaFloat,
    evaluate::{FunctionMap, OptimizationSettings},
    symbol,
    transcendental::{cosh, coth, csch, sech, sinh, tanh},
};

use super::{ActiveF64Backend, GenericEvaluator};
use crate::{
    initialisation::test_initialise,
    processes::EvaluatorSettings,
    settings::global::{
        CompilationOptimizationLevel, FrozenCompilationMode, GammaloopCompileOptions,
    },
    utils::{ArbPrec, F, FloatLike, QuadFloat},
};

const REFERENCE_PRECISION: u32 = 2048;
const FUNCTIONS: [&str; 6] = ["sinh", "cosh", "tanh", "coth", "sech", "csch"];

/// Independent MPC values and analytic derivatives, in function-major order.
///
/// For the extreme cases |Re(z)| >= 10^12, even MPC overflows sinh/cosh.
/// Their reciprocals correctly underflow: sech/csch and their first two
/// derivatives have magnitude <= 16 exp(-|Re(z)|), while the tanh/coth
/// derivatives are bounded by 16 exp(-2|Re(z)|), below MPC's exponent range.
/// MPC also preserves the zero sinh/cosh components on the real axis.
fn reference_values(input: &Mpc, precision: u32) -> Vec<Mpc> {
    // Near zero, derivative components can cancel inverse powers of subnormal
    // inputs. Guard bits cover that cancellation through the second derivative.
    let input = Mpc::with_val(precision + 4096, input);
    let sinh = input.clone().sinh();
    let cosh = input.clone().cosh();
    let tanh = input.tanh();
    let coth = tanh.clone().recip();
    let sech = cosh.clone().recip();
    let csch = sinh.clone().recip();
    let sech_squared = sech.clone().square();
    let csch_squared = csch.clone().square();

    // Using the reciprocal squares preserves exponentially small derivatives
    // after tanh/coth have rounded to +/-1; subtracting their squares would not.
    let tanh_second = tanh.clone() * &sech_squared * -2;
    let coth_first = -csch_squared.clone();
    let coth_second = coth.clone() * &csch_squared * 2;
    let sech_first = -sech.clone() * &tanh;
    let sech_second = sech.clone() * (1 - sech_squared.clone() * 2);
    let csch_first = -csch.clone() * &coth;
    let csch_second = csch.clone() * (1 + csch_squared * 2);

    [
        sinh.clone(),
        cosh.clone(),
        sinh.clone(),
        cosh.clone(),
        sinh,
        cosh,
        tanh,
        sech_squared,
        tanh_second,
        coth,
        coth_first,
        coth_second,
        sech,
        sech_first,
        sech_second,
        csch,
        csch_first,
        csch_second,
    ]
    .into_iter()
    .map(|value| Mpc::with_val(precision, value))
    .collect()
}

fn cases() -> Vec<(&'static str, Mpc)> {
    let mut cases = Vec::new();
    let mut add = |regime, re, im| {
        // These are deliberately binary64 boundary inputs, lifted exactly to MPC.
        cases.push((regime, Mpc::with_val(REFERENCE_PRECISION, (re, im))));
    };
    for re in [-1.0, -0.75, -0.0, 0.0, 0.75, 1.0] {
        for im in [-1.0, -0.5, -0.0, 0.0, 0.5, 1.0] {
            add("ordinary/axes", re, im);
        }
    }
    for magnitude in [
        2.0_f64.powi(-30),
        2.0_f64.powi(-500),
        2.0_f64.powi(-600),
        1e-300,
        f64::MIN_POSITIVE,
        f64::MIN_POSITIVE.next_down(),
        f64::from_bits(1),
    ] {
        for x in [-magnitude, magnitude] {
            for (re, im) in [(x, 0.0), (0.0, x), (x, x), (x, -x)] {
                add("tiny", re, im);
            }
        }
    }
    let log_max = f64::MAX.ln();
    let log_two = 2.0_f64.ln();
    for boundary in [
        20.0_f64,
        355.0,
        356.0,
        709.0,
        710.0,
        711.0,
        740.0,
        745.0,
        746.0,
        (log_max + log_two) / 2.0,
        log_max / 2.0 + log_two,
        log_max,
        log_max + log_two,
        -f64::MIN_POSITIVE.ln() + log_two,
        -f64::from_bits(1).ln() + log_two,
        -f64::from_bits(1).ln() + 2.0 * log_two,
    ] {
        for x in [boundary.next_down(), boundary, boundary.next_up()] {
            for re in [-x, x] {
                for im in [0.0, 0.5, -1.0] {
                    add("range boundary", re, im);
                }
            }
        }
    }
    for x in [1000.0, 1e12, f64::MAX / 2.0, f64::MAX] {
        for re in [-x, x] {
            for im in [0.0, 1.0] {
                add("UV", re, im);
            }
        }
    }
    for y in [1e6, 1e12, 1e100, f64::MAX / 2.0, f64::MAX] {
        for im in [-y, y] {
            for re in [0.0, 0.75, -1.0] {
                add("argument reduction", re, im);
            }
        }
    }
    for x in [356.0, 711.0, 1000.0] {
        for y in [1e-300, f64::from_bits(1)] {
            for (re, im) in [(x, y), (-x, y), (y, x), (y, -x)] {
                add("unequal components", re, im);
            }
        }
    }
    let mut rng = StdRng::seed_from_u64(0x4859_5045_5242_4f4c);
    for _ in 0..64 {
        let re = rng.random_range(-1.0..1.0) * 2.0_f64.powi(rng.random_range(-1074..=10));
        let im = rng.random_range(-1.0..1.0) * 2.0_f64.powi(rng.random_range(-40..=40));
        add("seeded magnitudes", re, im);
    }

    // Construct pi and perturbations at reference precision, then round to each
    // actual input domain. Rounded multiples of pi are not exact numeric poles.
    let pi = Float::with_val(REFERENCE_PRECISION, Constant::Pi);
    for k in -3..=3 {
        let center: Float = Float::with_val(REFERENCE_PRECISION, &pi * k) / 2;
        for exponent in [10, 30, 100, 500, 1000] {
            let delta = Float::with_val(REFERENCE_PRECISION, 1) >> exponent;
            for sign in [-1, 1] {
                let offset = Float::with_val(REFERENCE_PRECISION, &delta * sign);
                cases.push((
                    "near poles",
                    Mpc::with_val(REFERENCE_PRECISION, (&offset, &center)),
                ));
                let im = Float::with_val(REFERENCE_PRECISION, &center + &offset);
                cases.push(("near poles", Mpc::with_val(REFERENCE_PRECISION, (0, im))));
            }
        }
        let binary_center = center.to_f64();
        for im in [
            binary_center.next_down(),
            binary_center,
            binary_center.next_up(),
        ] {
            cases.push((
                "pole neighbors",
                Mpc::with_val(REFERENCE_PRECISION, (0, im)),
            ));
        }
    }
    for exponent in [80, 200, 900] {
        let offset = Float::with_val(REFERENCE_PRECISION, 1) >> exponent;
        let re = Float::with_val(REFERENCE_PRECISION, 1) + &offset;
        cases.push((
            "precision",
            Mpc::with_val(REFERENCE_PRECISION, (re, offset)),
        ));
    }
    cases
}

fn to_mpc<T: FloatLike>(value: Complex<F<T>>) -> Mpc {
    let parts = [value.re, value.im].map(|component| {
        let symbolic: SymbolicaFloat = component.into();
        symbolic.into_inner()
    });
    Mpc::with_val(REFERENCE_PRECISION, (&parts[0], &parts[1]))
}

struct Checks {
    precision: u32,
    checked: usize,
    failures: BTreeMap<(bool, String), (usize, String)>,
}

impl Checks {
    fn compare(&mut self, regime: &str, index: usize, input: &Mpc, actual: &Mpc, expected: &Mpc) {
        for (component, actual, reference) in [
            ("re", actual.real(), expected.real()),
            ("im", actual.imag(), expected.imag()),
        ] {
            self.checked += 1;
            // Double-double has ~106 mantissa bits but binary64 exponent range.
            // Its relative precision decreases near the subnormal boundary.
            let expected = if self.precision == 53
                || (self.precision == 106 && reference.clone().abs() > f64::MAX)
            {
                Float::with_val(REFERENCE_PRECISION, reference.to_f64())
            } else {
                reference.clone()
            };
            let minimum = if self.precision <= 106 {
                Float::with_val(REFERENCE_PRECISION, 1) >> 1074
            } else {
                Float::with_val(REFERENCE_PRECISION, 1) << (rug::float::exp_min() - 1)
            };
            let error = Float::with_val(REFERENCE_PRECISION, actual - &expected);
            let mut tolerance =
                (expected.clone().abs() >> (self.precision - 8)) + minimum.clone() * 2;
            let (valid, reason) = if expected.is_nan() {
                (actual.is_nan(), "expected NaN")
            } else if expected.is_infinite() {
                (
                    actual.is_infinite()
                        && actual.is_sign_negative() == expected.is_sign_negative(),
                    "expected signed infinity",
                )
            } else if !actual.is_finite() {
                (false, "unexpected nonfinite result")
            } else if actual.is_zero() && !expected.is_zero() {
                // A small absolute tolerance must never hide a representable tail.
                tolerance = minimum / 2;
                (expected.clone().abs() <= tolerance, "premature zero")
            } else {
                (error.clone().abs() <= tolerance, "component accuracy")
            };
            if !valid {
                self.record(
                    regime,
                    index,
                    format!(
                        "z={input:.8e}, {component}: got {actual:.16e}, expected {expected:.16e}, error={error:.3e}, tolerance={tolerance:.3e}; {reason}",
                    ),
                );
            }
        }
    }

    fn record(&mut self, regime: &str, index: usize, detail: String) {
        let key = (
            !index.is_multiple_of(3),
            format!("{regime}: {} d{}", FUNCTIONS[index / 3], index % 3),
        );
        let entry = self.failures.entry(key).or_insert((0, detail));
        entry.0 += 1;
    }

    fn finish(self, backend: &str) {
        let sections = [
            (false, "Built-in function failures"),
            (
                true,
                "Symbolic derivative failures (composed expressions; cancellation may reduce accuracy)",
            ),
        ];
        let report = sections
            .into_iter()
            .filter_map(|(derivative, title)| {
                let rows = self
                    .failures
                    .iter()
                    .filter(|((is_derivative, _), _)| *is_derivative == derivative)
                    .map(|((_, key), (count, example))| {
                        format!("{key}: {count} failures; first {example}")
                    })
                    .collect::<Vec<_>>();
                (!rows.is_empty()).then(|| format!("{title}:\n{}", rows.join("\n")))
            })
            .collect::<Vec<_>>()
            .join("\n");
        assert!(
            self.failures.is_empty(),
            "{backend}: {} component checks\n{report}",
            self.checked
        );
    }
}

fn check_backend<T: FloatLike>(
    mode: FrozenCompilationMode,
    precision: u32,
    from_rug: impl Fn(&Float) -> F<T>,
) {
    test_initialise().unwrap();
    let parameter = symbol!("hyperbolic_robustness_z");
    let expressions = [sinh(), cosh(), tanh(), coth(), sech(), csch()]
        .into_iter()
        .flat_map(|function| {
            let value = function.call_args([Atom::var(parameter)]);
            let first = value.derivative(parameter);
            let second = first.derivative(parameter);
            [value, first, second]
        });
    let mut evaluator = GenericEvaluator::new_from_raw_params(
        expressions,
        &[Atom::var(parameter)],
        &FunctionMap::default(),
        vec![],
        OptimizationSettings::default(),
        // Symbolic derivatives work independently of the currently unsupported
        // native hyperbolic hyperdual-vectorization callbacks.
        None,
        &EvaluatorSettings::default(),
    )
    .unwrap();
    let backend = format!("{mode}, {precision} bits");
    // Generate binary64 boundary values before loading fast-math code, which
    // may change the process floating-point environment on some platforms.
    let cases = cases();
    let directory = std::env::temp_dir().join(format!(
        "gammaloop-hyperbolic-{}-{precision}-{}-{}",
        mode.active_backend_name(),
        mode.external_options()
            .is_some_and(|options| options.fast_math),
        std::process::id(),
    ));
    let expected_backend = match &mode {
        FrozenCompilationMode::Eager => ActiveF64Backend::Eager,
        FrozenCompilationMode::Symjit(level) => {
            evaluator.activate_symjit(*level).unwrap();
            ActiveF64Backend::Symjit
        }
        FrozenCompilationMode::Cpp(_) => {
            std::fs::create_dir_all(&directory).unwrap();
            evaluator
                .compile_external(
                    directory.join("hyperbolic.cpp"),
                    "hyperbolic",
                    directory.join("hyperbolic.so"),
                    &mode,
                )
                .unwrap();
            ActiveF64Backend::Cpp
        }
        FrozenCompilationMode::Assembly(_) => unreachable!("C++ is tested without inline assembly"),
    };
    assert_eq!(evaluator.active_f64_backend(), expected_backend);
    let mut evaluate = |input: &Mpc| {
        let argument = Complex::new(from_rug(input.real()), from_rug(input.imag()));
        let actual_input = to_mpc(argument.clone());
        let output = T::get_evaluator(&mut evaluator)(&[argument])
            .into_iter()
            .map(|v| to_mpc(v.unwrap_real()))
            .collect::<Vec<_>>();
        assert_eq!(output.len(), FUNCTIONS.len() * 3);
        (actual_input, output)
    };
    let mut checks = Checks {
        precision,
        checked: 0,
        failures: BTreeMap::new(),
    };
    let mut oracle_checks = Checks {
        precision: REFERENCE_PRECISION - 32,
        checked: 0,
        failures: BTreeMap::new(),
    };
    for (regime, input) in cases {
        let (input, actual) = evaluate(&input);
        let expected = reference_values(&input, REFERENCE_PRECISION);
        let refined = reference_values(&input, REFERENCE_PRECISION * 2);
        for index in 0..actual.len() {
            let function = index / 3;
            if input.real().is_zero() && input.imag().is_zero() && [3, 5].contains(&function) {
                // There is no direction-independent value or derivative at a pole.
                if index % 3 == 0
                    && actual[index].real().is_finite()
                    && actual[index].imag().is_finite()
                {
                    checks.record(
                        "exact pole",
                        index,
                        format!("z={input}, got {}", actual[index]),
                    );
                }
                continue;
            }
            if refined[index].real().is_nan() || refined[index].imag().is_nan() {
                oracle_checks.record(
                    "unexpected oracle NaN",
                    index,
                    format!("z={input:.8e}, got {:.8e}", refined[index]),
                );
            }
            oracle_checks.compare(
                "oracle convergence",
                index,
                &input,
                &expected[index],
                &refined[index],
            );
            checks.compare(regime, index, &input, &actual[index], &refined[index]);
        }
    }

    // Exact transformations avoid the conditioning problems of identities near
    // poles. A pi shift is tested only at these moderate, regular arguments.
    for (re, im) in [(1, 1), (-1, 1), (1, -1), (-1, -1)] {
        let input = Mpc::with_val(REFERENCE_PRECISION, (re, im));
        let (_, base) = evaluate(&input);
        for (regime, transformed) in [
            ("parity", -input.clone()),
            ("conjugation", input.clone().conj()),
            ("periodicity", {
                let mut shifted = input.clone();
                *shifted.mut_imag() += Float::with_val(REFERENCE_PRECISION, Constant::Pi);
                shifted
            }),
        ] {
            let (argument, result) = evaluate(&transformed);
            for (index, actual) in result.iter().enumerate() {
                let odd = [0, 2, 3, 5].contains(&(index / 3));
                let expected = match regime {
                    "conjugation" => base[index].clone().conj(),
                    "parity" if odd != (index % 3 == 1) => -base[index].clone(),
                    "periodicity" if ![2, 3].contains(&(index / 3)) => -base[index].clone(),
                    _ => base[index].clone(),
                };
                checks.compare(regime, index, &argument, actual, &expected);
            }
        }
    }

    // Fast math does not promise IEEE nonfinite or signed-zero semantics. Its
    // finite robustness cases above have exactly the same accuracy requirements.
    let strict = !matches!(&mode, FrozenCompilationMode::Cpp(options) if options.fast_math || options.unsafe_math);
    if strict {
        for (re, im) in [
            (f64::NAN, 0.0),
            (0.0, f64::NAN),
            (f64::NAN, 1.0),
            (1.0, f64::NAN),
            (f64::INFINITY, 0.0),
            (f64::NEG_INFINITY, 0.0),
            (f64::INFINITY, 1.0),
            (f64::NEG_INFINITY, -1.0),
        ] {
            let (input, result) = evaluate(&Mpc::with_val(REFERENCE_PRECISION, (re, im)));
            for (function, actual) in result.iter().step_by(3).enumerate() {
                if re.is_nan() || im.is_nan() {
                    if !actual.real().is_nan() && !actual.imag().is_nan() {
                        checks.record(
                            "NaN input",
                            function * 3,
                            format!("z={input}, got {actual}"),
                        );
                    }
                } else if function >= 2 {
                    let expected = Mpc::with_val(
                        REFERENCE_PRECISION,
                        (if function < 4 { re.signum() } else { 0.0 }, 0),
                    );
                    checks.compare(
                        "infinite real part",
                        function * 3,
                        &input,
                        actual,
                        &expected,
                    );
                }
            }
        }
        for (re, im) in [(0.0, f64::INFINITY), (1.0, f64::NEG_INFINITY)] {
            let (input, result) = evaluate(&Mpc::with_val(REFERENCE_PRECISION, (re, im)));
            for (function, actual) in result.iter().step_by(3).enumerate() {
                if !actual.real().is_nan() && !actual.imag().is_nan() {
                    checks.record(
                        "infinite imaginary part",
                        function * 3,
                        format!("z={input}, got {actual}"),
                    );
                }
            }
        }
        // Mixed infinities/NaNs have no single analytic limiting direction.
        // Exercise them for evaluation/ABI errors without prescribing a value.
        for re in [f64::INFINITY, f64::NEG_INFINITY, f64::NAN] {
            for im in [f64::INFINITY, f64::NEG_INFINITY, f64::NAN] {
                let _ = evaluate(&Mpc::with_val(REFERENCE_PRECISION, (re, im)));
            }
        }
        if precision == 53 {
            for re in [-0.0_f64, 0.0] {
                for im in [-0.0_f64, 0.0] {
                    let (_, result) = evaluate(&Mpc::with_val(REFERENCE_PRECISION, (re, im)));
                    for function in [0, 2] {
                        let actual = &result[function * 3];
                        if actual.real().is_sign_negative() != re.is_sign_negative()
                            || actual.imag().is_sign_negative() != im.is_sign_negative()
                        {
                            checks.record(
                                "signed zero",
                                function * 3,
                                format!("z=({re:?},{im:?}), got {actual}"),
                            );
                        }
                    }
                }
            }
        }
    }
    if mode.requires_external_compilation() {
        std::fs::remove_dir_all(directory).unwrap();
    }
    oracle_checks.finish("MPC reference convergence");
    checks.finish(&backend);
}

#[test]
fn hyperbolic_eager_f64() {
    check_backend(FrozenCompilationMode::Eager, 53, |x| F(x.to_f64()));
}

#[test]
fn hyperbolic_eager_quad() {
    check_backend(FrozenCompilationMode::Eager, 106, |x| {
        F(QuadFloat::from(x.clone()))
    });
}

#[test]
fn hyperbolic_eager_arb() {
    check_backend(FrozenCompilationMode::Eager, 1000, |x| {
        F(ArbPrec::from(x.clone()))
    });
}

#[test]
fn hyperbolic_cpp_strict() {
    let options = GammaloopCompileOptions {
        fast_math: false,
        unsafe_math: false,
        ..Default::default()
    };
    check_backend(
        FrozenCompilationMode::Cpp(options.external_options_snapshot()),
        53,
        |x| F(x.to_f64()),
    );
}

#[test]
fn hyperbolic_cpp_production() {
    let options = GammaloopCompileOptions::default();
    check_backend(
        FrozenCompilationMode::Cpp(options.external_options_snapshot()),
        53,
        |x| F(x.to_f64()),
    );
}

#[test]
fn hyperbolic_symjit_o0() {
    check_backend(
        FrozenCompilationMode::Symjit(CompilationOptimizationLevel::O0),
        53,
        |x| F(x.to_f64()),
    );
}

#[test]
fn hyperbolic_symjit_o2() {
    check_backend(
        FrozenCompilationMode::Symjit(CompilationOptimizationLevel::O2),
        53,
        |x| F(x.to_f64()),
    );
}
