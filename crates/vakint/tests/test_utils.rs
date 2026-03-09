use std::f64::consts::LOG2_10;

use colored::Colorize;
use log::{debug, info};
use symbolica::{
    atom::{Atom, AtomView},
    domains::float::{Complex, Float},
    printer::{AtomPrinter, PrintOptions},
    try_parse,
};

use std::collections::HashMap;
use std::process::{Command, Stdio};
use std::sync::{Once, OnceLock};
use vakint::{EvaluationMethod, NumericalEvaluationResult, Vakint, VakintError};
use vakint::{EvaluationOrder, LoopNormalizationFactor, Momentum, VakintSettings};

pub struct TestVakint {
    pub vakint: Vakint,
    pub settings: VakintSettings,
}

impl std::ops::Deref for TestVakint {
    type Target = Vakint;

    fn deref(&self) -> &Self::Target {
        &self.vakint
    }
}

impl std::ops::DerefMut for TestVakint {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.vakint
    }
}

// Helpers are shared across test binaries; not every binary uses every helper, so allow dead_code.
#[allow(dead_code)]
impl TestVakint {
    pub fn to_canonical(&self, input: AtomView, short_form: bool) -> Result<Atom, VakintError> {
        self.vakint.to_canonical(&self.settings, input, short_form)
    }

    pub fn tensor_reduce(&self, input: AtomView) -> Result<Atom, VakintError> {
        self.vakint.tensor_reduce(&self.settings, input)
    }

    pub fn evaluate_integral(&self, input: AtomView) -> Result<Atom, VakintError> {
        self.vakint.evaluate_integral(&self.settings, input)
    }

    pub fn evaluate(&self, input: AtomView) -> Result<Atom, VakintError> {
        self.vakint.evaluate(&self.settings, input)
    }

    pub fn params_from_f64(
        &self,
        params: &HashMap<String, f64>,
    ) -> HashMap<String, Float, ahash::RandomState> {
        self.vakint.params_from_f64(&self.settings, params)
    }

    pub fn params_from_complex_f64(
        &self,
        params: &HashMap<String, Complex<f64>>,
    ) -> HashMap<String, Complex<Float>, ahash::RandomState> {
        self.vakint.params_from_complex_f64(&self.settings, params)
    }

    pub fn externals_from_f64(
        &self,
        externals: &HashMap<usize, (f64, f64, f64, f64)>,
    ) -> HashMap<usize, Momentum, ahash::RandomState> {
        self.vakint.externals_from_f64(&self.settings, externals)
    }

    #[allow(clippy::type_complexity)]
    pub fn externals_from_complex_f64(
        &self,
        externals: &HashMap<usize, (Complex<f64>, Complex<f64>, Complex<f64>, Complex<f64>)>,
    ) -> HashMap<usize, Momentum, ahash::RandomState> {
        self.vakint
            .externals_from_complex_f64(&self.settings, externals)
    }

    pub fn numerical_evaluation(
        &self,
        expression: AtomView,
        params_real: &HashMap<String, Float, ahash::RandomState>,
        params_complex: &HashMap<String, Complex<Float>, ahash::RandomState>,
        externals: Option<&HashMap<usize, Momentum, ahash::RandomState>>,
    ) -> Result<(NumericalEvaluationResult, Option<NumericalEvaluationResult>), VakintError> {
        self.vakint.numerical_evaluation(
            &self.settings,
            expression,
            params_real,
            params_complex,
            externals,
        )
    }
}

pub fn get_vakint(mut vakint_settings: VakintSettings) -> TestVakint {
    if evaluation_requires_pysecdec(&vakint_settings.evaluation_order)
        && !pysecdec_available(&vakint_settings.python_exe_path)
    {
        vakint_settings
            .evaluation_order
            .0
            .retain(|method| !matches!(method, EvaluationMethod::PySecDec(_)));
    }
    static VAKINT_INSTANCE: OnceLock<Vakint> = OnceLock::new();
    let vakint = VAKINT_INSTANCE
        .get_or_init(|| Vakint::new().expect("Failed to initialize vakint"))
        .clone();
    TestVakint {
        vakint,
        settings: vakint_settings,
    }
}

#[allow(dead_code)]
pub fn should_skip_pysecdec_tests() -> bool {
    static WARNED: Once = Once::new();
    let skip = match std::env::var("VAKINT_SKIP_PYSECDEC_TESTS") {
        Ok(value) => {
            let trimmed = value.trim();
            !(trimmed.is_empty() || trimmed == "0" || trimmed.eq_ignore_ascii_case("false"))
        }
        Err(_) => false,
    };
    if skip {
        WARNED.call_once(|| {
            eprintln!("Skipping PySecDec tests because VAKINT_SKIP_PYSECDEC_TESTS is set.");
        });
    }
    skip
}

fn pysecdec_available(python_exe: &str) -> bool {
    for probe in ["import pySecDec", "import pysecdec"] {
        let ok = Command::new(python_exe)
            .arg("-c")
            .arg(probe)
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .map(|status| status.success())
            .unwrap_or(false);
        if ok {
            return true;
        }
    }
    false
}

fn evaluation_requires_pysecdec(evaluation_order: &EvaluationOrder) -> bool {
    evaluation_order
        .0
        .iter()
        .any(|method| matches!(method, EvaluationMethod::PySecDec(_)))
}

#[allow(unused)]
pub fn compare_output(output: Result<AtomView, &VakintError>, expected_output: Atom) -> Atom {
    match output {
        Ok(r) => {
            let r_processed = try_parse!(
                AtomPrinter::new_with_options(
                    r,
                    PrintOptions {
                        hide_namespace: Some("tests"),
                        ..PrintOptions::file()
                    }
                )
                .to_string()
                .as_str()
            )
            .unwrap();
            let r_processed_view = r_processed.as_view();
            if r != expected_output.as_view() {
                println!(
                    "Output does not match expected output:\n{}\n!=\n{}",
                    format!(
                        "{}",
                        AtomPrinter::new_with_options(r_processed_view, PrintOptions::file())
                    )
                    .red(),
                    format!(
                        "{}",
                        AtomPrinter::new_with_options(
                            expected_output.as_view(),
                            PrintOptions::file()
                        )
                    )
                    .green()
                );
            }
            assert_eq!(r, expected_output.as_view());
            r.to_owned()
        }
        Err(err) => panic!("Error: {}", err),
    }
}

#[allow(unused)]
pub fn compare_numerical_output(
    output: Result<&NumericalEvaluationResult, &VakintError>,
    expected_output: Vec<(i64, (String, String))>,
    prec: u32,
) {
    let binary_prec: u32 = ((prec as f64) * LOG2_10).floor() as u32;

    match output {
        Ok(numerical_result) => {
            let r = numerical_result.get_epsilon_coefficients();
            if r.len() != expected_output.len() {
                println!(
                    "Output does not match expected output: length mismatch: {} != {}",
                    r.len(),
                    expected_output.len()
                );
            }
            assert!(r.len() == expected_output.len());
            for ((o_pwr, o), (e_pwr, (e_real, e_cmplx))) in r.iter().zip(expected_output) {
                if *o_pwr != e_pwr {
                    println!("Power mismatch: {} != {}", o_pwr, e_pwr);
                    assert_eq!(*o_pwr, e_pwr);
                }
                let trgt = Complex::new(
                    Float::parse(e_real.as_str(), Some(binary_prec)).unwrap(),
                    Float::parse(e_cmplx.as_str(), Some(binary_prec)).unwrap(),
                );
                let scale = if trgt.norm_squared() > Float::with_val(binary_prec, 0.0) {
                    trgt.norm_squared()
                } else {
                    Float::with_val(binary_prec, 0.0)
                };
                let o_prec = (o.clone() - trgt.clone()).norm_squared() / scale;
                let trgt_prec = Float::with_val(binary_prec, (10.0_f64).powi(-(prec as i32)));
                if o_prec > trgt_prec {
                    println!(
                        "Output does not match expected output:\n{}\n!=\n{} (error: {} > target precision {})",
                        format!("{}", o).red(),
                        format!("{}", trgt).green(),
                        o_prec,
                        trgt_prec
                    );
                }
                assert!(o_prec < trgt_prec);
            }
        }
        Err(err) => panic!("Error: {}", err),
    }
}

// EvaluationOrder::analytic_only()

#[allow(unused, clippy::too_many_arguments)]
pub fn compare_two_evaluations(
    vakint_default_settings: VakintSettings,
    evaluation_orders: ((&EvaluationOrder, bool), (&EvaluationOrder, bool)),
    integra_view: AtomView,
    numerical_masses: HashMap<String, Float, ahash::RandomState>,
    numerical_external_momenta: HashMap<usize, Momentum, ahash::RandomState>,
    rel_threshold: f64,
    max_pull: f64,
    quiet: bool,
) {
    let python_exe = vakint_default_settings.python_exe_path.clone();
    let mut mod_evaluation_order_a = evaluation_orders.0.0.clone();
    let mut mod_evaluation_order_b = evaluation_orders.1.0.clone();
    for eval_order in [&mut mod_evaluation_order_a, &mut mod_evaluation_order_b] {
        eval_order.adjust(
            Some(quiet),
            rel_threshold * 1.0e-2,
            &numerical_masses,
            &HashMap::default(),
            &numerical_external_momenta,
        );
    }

    // First perform the first evaluation type

    let mut vakint_analytic_settings = VakintSettings {
        allow_unknown_integrals: false,
        use_dot_product_notation: true,
        mu_r_sq_symbol: "mursq".into(),
        integral_normalization_factor: LoopNormalizationFactor::pySecDec,
        evaluation_order: mod_evaluation_order_a.clone(),
        ..vakint_default_settings
    };
    if (evaluation_requires_pysecdec(&mod_evaluation_order_a)
        || evaluation_requires_pysecdec(&mod_evaluation_order_b))
        && !pysecdec_available(&python_exe)
    {
        eprintln!("Skipping test: PySecDec not available.");
        return;
    }
    let mut vakint = get_vakint(vakint_analytic_settings);

    let mut eval_params = HashMap::default();
    if numerical_masses.contains_key("user_space::muv") {
        eval_params.insert(
            "user_space::muv".into(),
            numerical_masses
                .get("user_space::muv")
                .unwrap_or_else(|| panic!("user_space::muv not found in numerical_masses"))
                .to_owned(),
        );
    } else {
        eval_params.insert(
            "muvsq".into(),
            numerical_masses
                .get("muvsq")
                .unwrap_or_else(|| panic!("muvsq not found in numerical_masses"))
                .to_owned(),
        );
    }
    eval_params.insert(
        "mursq".into(),
        numerical_masses
            .get("mursq")
            .unwrap_or_else(|| panic!("mursq not found in numerical_masses"))
            .to_owned(),
    );

    let integral = vakint.to_canonical(integra_view, true).unwrap();

    let integral_reduced = if [evaluation_orders.0, evaluation_orders.1]
        .iter()
        .any(|(_a, do_reduction)| *do_reduction)
    {
        vakint.tensor_reduce(integral.as_view()).unwrap()
    } else {
        integral.clone()
    };

    debug!(
        "Evaluating integral with evaluation order: {}",
        format!("{}", mod_evaluation_order_a).green(),
    );
    let benchmark_evaluated_integral = match vakint.evaluate_integral(if evaluation_orders.0.1 {
        integral_reduced.as_view()
    } else {
        integral.as_view()
    }) {
        Ok(eval) => eval,
        Err(e) => {
            panic!(
                "Error during benchmark integral evaluation with {} :: error :\n{}",
                format!("{}", mod_evaluation_order_a).red(),
                e
            );
        }
    };

    let (benchmark_central, benchmark_error) = match Vakint::full_numerical_evaluation(
        &vakint.settings,
        benchmark_evaluated_integral.as_view(),
        &eval_params,
        &HashMap::default(),
        Some(&numerical_external_momenta),
    ) {
        Ok(eval) => eval,
        Err(e) => {
            panic!(
                "Error during numerical evaluation of benchmark result with {} :: error:\n{}",
                format!("{}", mod_evaluation_order_a).red(),
                e
            );
        }
    };
    debug!(
        "Benchmark {} :: central :\n{}",
        format!("{}", mod_evaluation_order_a).green(),
        benchmark_central.clone()
    );
    if let Some(err) = &benchmark_error {
        debug!(
            "Benchmark {} :: error   :\n{}",
            format!("{}", mod_evaluation_order_a).green(),
            err.clone()
        );
    }

    // Now perform the second evaluation
    vakint.settings.evaluation_order = mod_evaluation_order_b.clone();
    debug!(
        "Evaluating integral with evaluation order: {}",
        format!("{}", mod_evaluation_order_b).green(),
    );
    let comparison_eval = match vakint.evaluate_integral(if evaluation_orders.1.1 {
        integral_reduced.as_view()
    } else {
        integral.as_view()
    }) {
        Ok(eval) => eval,
        Err(e) => {
            panic!(
                "Error during comparison valuation with: {} :: error :\n{}",
                format!("{}", mod_evaluation_order_b).red(),
                e
            );
        }
    };
    let (tested_central, tested_error) = match Vakint::full_numerical_evaluation(
        &vakint.settings,
        comparison_eval.as_view(),
        &eval_params,
        &HashMap::default(),
        Some(&numerical_external_momenta),
    ) {
        Ok(eval) => eval,
        Err(e) => {
            panic!(
                "Error during numerical evaluation of benchmark result with {} :: error:\n{}",
                format!("{}", mod_evaluation_order_b).red(),
                e
            );
        }
    };
    debug!(
        "Tested {} :: central :\n{}",
        format!("{}", mod_evaluation_order_b).green(),
        tested_central.clone()
    );
    if let Some(err) = &tested_error {
        debug!(
            "Tested {} :: error   :\n{}",
            format!("{}", mod_evaluation_order_b).green(),
            err.clone()
        );
    }

    let mut combined_error = match (&benchmark_error, &tested_error) {
        (Some(b), Some(t)) => Some(b.aggregate_errors(t)),
        (Some(b), None) => Some(b.clone()),
        (None, Some(t)) => Some(t.clone()),
        (None, None) => None,
    };
    let (matches, msg) = benchmark_central.does_approx_match(
        &tested_central,
        combined_error.as_ref(),
        rel_threshold,
        max_pull,
    );
    if !matches || !quiet {
        println!("\n{}\n", "<><><><><>".green());
        println!(
            "Benchmark {} :: central :\n{}",
            format!("{}", mod_evaluation_order_a).green(),
            benchmark_central
        );
        if let Some(err) = &benchmark_error {
            println!(
                "Benchmark {} :: error   :\n{}",
                format!("{}", mod_evaluation_order_a).green(),
                err
            );
        }
        println!(
            "Tested {} :: central :\n{}",
            format!("{}", mod_evaluation_order_b).green(),
            tested_central
        );
        if let Some(err) = &tested_error {
            println!(
                "Tested {} :: error   :\n{}",
                format!("{}", mod_evaluation_order_b).green(),
                err
            );
        }
        println!("{}", msg);
        println!("\n{}\n", "<><><><><>".green());
    } else {
        println!("\n{}\n", "<><><><><>".green());
        info!(
            "Benchmark {} :: central :\n{}",
            format!("{}", mod_evaluation_order_a).green(),
            benchmark_central
        );
        if let Some(err) = &benchmark_error {
            info!(
                "Benchmark {} :: error   :\n{}",
                format!("{}", mod_evaluation_order_a).green(),
                err
            );
        }
        info!(
            "Tested {} :: central :\n{}",
            format!("{}", mod_evaluation_order_b).green(),
            tested_central
        );
        if let Some(err) = &tested_error {
            info!(
                "Tested {} :: error   :\n{}",
                format!("{}", mod_evaluation_order_b).green(),
                err
            );
        }
        info!("{}", msg);
        println!("\n{}\n", "<><><><><>".green());
    }
    assert!(matches, "Benchmark and numerical result do not match.");
}

#[allow(unused, clippy::too_many_arguments)]
pub fn compare_vakint_evaluation_vs_reference(
    vakint_default_settings: VakintSettings,
    evaluation_order: EvaluationOrder,
    integra_view: AtomView,
    numerical_masses: HashMap<String, Float, ahash::RandomState>,
    numerical_external_momenta: HashMap<usize, Momentum, ahash::RandomState>,
    expected_output: Vec<(i64, (String, String))>,
    prec: u32,
    max_pull: f64,
) {
    let python_exe = vakint_default_settings.python_exe_path.clone();
    // Adjust evaluation method options
    let mut mod_evaluation_order = evaluation_order.clone();
    mod_evaluation_order.adjust(
        None,
        10.0_f64.powi(-(prec as i32)),
        &numerical_masses,
        &HashMap::default(),
        &numerical_external_momenta,
    );

    // First perform Vakint evaluation
    let mut vakint_settings = VakintSettings {
        allow_unknown_integrals: true,
        use_dot_product_notation: true,
        //        mu_r_sq_symbol: "vakint::{}::mursq".into(),
        mu_r_sq_symbol: vakint_default_settings.mu_r_sq_symbol.clone(),
        run_time_decimal_precision: prec,
        evaluation_order: mod_evaluation_order.clone(),
        ..vakint_default_settings
    };

    if evaluation_requires_pysecdec(&mod_evaluation_order) && !pysecdec_available(&python_exe) {
        eprintln!("Skipping test: PySecDec not available.");
        return;
    }

    let mut vakint = get_vakint(vakint_settings);

    let integral = vakint.to_canonical(integra_view, true).unwrap();

    let integral = vakint.tensor_reduce(integral.as_view()).unwrap();

    let integral = vakint
        .evaluate_integral(integral.as_view())
        .unwrap_or_else(|op| panic!("Failed to evaluate integral: {}", op));

    let (result, error) = Vakint::full_numerical_evaluation(
        &vakint.settings,
        integral.as_view(),
        &numerical_masses,
        &HashMap::default(),
        Some(&numerical_external_momenta),
    )
    .unwrap();

    let binary_prec: u32 = ((prec.max(16) as f64) * LOG2_10).floor() as u32;

    let reference = NumericalEvaluationResult(
        expected_output
            .iter()
            .map(|(eps_pwr, (re, im))| {
                (
                    *eps_pwr,
                    Complex::new(
                        Float::parse(re.as_str(), Some(binary_prec)).unwrap(),
                        Float::parse(im.as_str(), Some(binary_prec)).unwrap(),
                    ),
                )
            })
            .collect::<Vec<_>>(),
    );
    let (matches, msg) = result.does_approx_match(
        &reference,
        error.as_ref(),
        0.1_f64.powi((prec as i32) - 2),
        max_pull,
    );
    if !matches {
        println!(
            "Result from {}:\n{}",
            format!("{}", mod_evaluation_order).green(),
            result
        );
        if let Some(err) = &error {
            println!("Error:\n{}", err);
        }
        println!("Reference:\n{}", reference);
        println!("{}", msg)
    } else {
        debug!(
            "Result from {}:\n{}",
            format!("{}", mod_evaluation_order).green(),
            result
        );
        if let Some(err) = &error {
            debug!("Error:\n{}", err);
        }
        debug!("Reference:\n{}", reference);
        debug!("{}", msg)
    }
    assert!(matches, "Vakint result and reference result do not match");
}
