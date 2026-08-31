use std::f64::consts::LOG2_10;

use colored::Colorize;
use log::{debug, info};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    domains::float::{Complex, Float},
    printer::{AtomPrinter, PrintOptions},
    try_parse,
};

use std::collections::HashMap;
use std::process::{Command, Stdio};
use std::sync::{Once, OnceLock};
use vakint::{
    EvaluationMethod, MATADOptions, NumericalEvaluationResult, RustRedEvaluationOptions, Vakint,
    VakintError,
};
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
    let skip = match std::env::var("RUN_PYSECDEC_TESTS") {
        Ok(value) => {
            let trimmed = value.trim();
            !(trimmed.is_empty() || trimmed == "1" || trimmed.eq_ignore_ascii_case("true"))
        }
        Err(_) => true,
    };
    if skip {
        WARNED.call_once(|| {
            eprintln!("Skipping PySecDec tests because RUN_PYSECDEC_TESTS is not set.");
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
                        hide_namespace: Some("tests".into()),
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

/// Which shared input one evaluation lane consumes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EvaluationTestInput {
    Canonical,
    TensorReduced,
}

/// Whether the comparison prepares one shared tensor-reduced input with FORM.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TensorPrepass {
    Skip,
    Form,
}

/// Acceptance authority for one strict RustRed scalar lane.
///
/// A common symbolic master basis permits the stronger exact-coefficient
/// comparison.  A certified finite, nonminimal RustRed terminal basis instead
/// uses numerical parity after each backend has substituted its own masters.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub enum RustRedParityPolicy {
    ExactMatadBasis,
    NumericalOnly,
}

/// One owned lane in the acceptance-test multi-method evaluator.
#[derive(Clone, Debug)]
pub struct EvaluationTestLane {
    evaluation_order: EvaluationOrder,
    input: EvaluationTestInput,
    forbid_form_after_prepass: bool,
    rustred_parity: Option<RustRedParityPolicy>,
}

#[allow(dead_code)]
impl EvaluationTestLane {
    pub fn oracle(evaluation_order: EvaluationOrder, input: EvaluationTestInput) -> Self {
        Self {
            evaluation_order,
            input,
            forbid_form_after_prepass: false,
            rustred_parity: None,
        }
    }

    /// Strict RustRed scalar lane with the requested acceptance authority.
    ///
    /// [`RustRedParityPolicy::NumericalOnly`] still forbids FORM and fallback in
    /// the scalar tail; it only disables an invalid raw symbolic comparison
    /// when RustRed intentionally uses a different finite master basis.
    pub fn rustred_scalar(input: EvaluationTestInput, parity: RustRedParityPolicy) -> Self {
        Self {
            evaluation_order: EvaluationOrder::rustred_only(),
            input,
            forbid_form_after_prepass: true,
            rustred_parity: Some(parity),
        }
    }

    pub fn rustred_parity(&self) -> Option<RustRedParityPolicy> {
        self.rustred_parity
    }

    pub fn forbids_form_after_prepass(&self) -> bool {
        self.forbid_form_after_prepass
    }
}

/// The three strict scalar-reduction peers used by the Stage-1 acceptance suite.
#[allow(dead_code)]
pub fn alphaloop_matad_lanes() -> [EvaluationTestLane; 2] {
    [
        EvaluationTestLane::oracle(
            EvaluationOrder::alphaloop_only(),
            EvaluationTestInput::TensorReduced,
        ),
        EvaluationTestLane::oracle(
            EvaluationOrder::matad_only(None),
            EvaluationTestInput::TensorReduced,
        ),
    ]
}

/// The existing AlphaLoop/MATAD oracle pair, extended by a strict RustRed lane.
#[allow(dead_code)]
pub fn alphaloop_matad_rustred_lanes(
    rustred_parity: RustRedParityPolicy,
) -> [EvaluationTestLane; 3] {
    let [alphaloop, matad] = alphaloop_matad_lanes();
    [
        alphaloop,
        matad,
        EvaluationTestLane::rustred_scalar(EvaluationTestInput::TensorReduced, rustred_parity),
    ]
}

/// The existing analytic/MATAD oracle pair.
#[allow(dead_code)]
pub fn analytic_matad_lanes() -> [EvaluationTestLane; 2] {
    [
        EvaluationTestLane::oracle(
            EvaluationOrder::analytic_only(),
            EvaluationTestInput::TensorReduced,
        ),
        EvaluationTestLane::oracle(
            EvaluationOrder::matad_only(None),
            EvaluationTestInput::TensorReduced,
        ),
    ]
}

/// Existing analytic/MATAD oracles, extended by a strict RustRed lane.
#[allow(dead_code)]
pub fn analytic_matad_rustred_lanes(
    rustred_parity: RustRedParityPolicy,
) -> [EvaluationTestLane; 3] {
    let [analytic, matad] = analytic_matad_lanes();
    [
        analytic,
        matad,
        EvaluationTestLane::rustred_scalar(EvaluationTestInput::TensorReduced, rustred_parity),
    ]
}

const FORBIDDEN_SCALAR_FORM_PATH: &str = "/this/path/must/not/be/invoked/by-rustred-acceptance";

fn adjusted_lanes(
    lanes: &[EvaluationTestLane],
    quiet: Option<bool>,
    relative_precision: f64,
    numerical_masses: &HashMap<String, Float, ahash::RandomState>,
    numerical_external_momenta: &HashMap<usize, Momentum, ahash::RandomState>,
) -> Vec<EvaluationTestLane> {
    lanes
        .iter()
        .cloned()
        .map(|mut lane| {
            lane.evaluation_order.adjust(
                quiet,
                relative_precision,
                numerical_masses,
                &HashMap::default(),
                numerical_external_momenta,
            );
            lane
        })
        .collect()
}

fn prepare_comparison_inputs(
    vakint: &mut TestVakint,
    integra_view: AtomView,
    tensor_prepass: TensorPrepass,
) -> (Atom, Atom) {
    let canonical = vakint.to_canonical(integra_view, true).unwrap();
    let tensor_reduced = match tensor_prepass {
        TensorPrepass::Skip => canonical.clone(),
        TensorPrepass::Form => vakint.tensor_reduce(canonical.as_view()).unwrap(),
    };
    (canonical, tensor_reduced)
}

fn evaluate_lane(
    vakint: &mut TestVakint,
    lane: &EvaluationTestLane,
    canonical: &Atom,
    tensor_reduced: &Atom,
    original_form_path: &str,
) -> Atom {
    vakint.settings.evaluation_order = lane.evaluation_order.clone();
    vakint.settings.form_exe_path = if lane.forbid_form_after_prepass {
        FORBIDDEN_SCALAR_FORM_PATH.to_owned()
    } else {
        original_form_path.to_owned()
    };
    let input = match lane.input {
        EvaluationTestInput::Canonical => canonical.as_view(),
        EvaluationTestInput::TensorReduced => tensor_reduced.as_view(),
    };
    vakint.evaluate_integral(input).unwrap_or_else(|error| {
        panic!(
            "integral evaluation with {} failed: {error}",
            lane.evaluation_order
        )
    })
}

fn assert_exact_raw_rustred_matad_peer(
    vakint: &mut TestVakint,
    lanes: &[EvaluationTestLane],
    canonical: &Atom,
    tensor_reduced: &Atom,
    original_form_path: &str,
) {
    let Some(rustred_input) = lanes
        .iter()
        .find(|lane| lane.rustred_parity == Some(RustRedParityPolicy::ExactMatadBasis))
        .map(|lane| lane.input)
    else {
        return;
    };
    let matad = EvaluationTestLane::oracle(
        EvaluationOrder::matad_only(Some(MATADOptions {
            expand_masters: false,
            ..MATADOptions::default()
        })),
        rustred_input,
    );
    let rustred = EvaluationTestLane {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                substitute_masters: false,
            },
        )]),
        input: rustred_input,
        forbid_form_after_prepass: true,
        rustred_parity: Some(RustRedParityPolicy::ExactMatadBasis),
    };
    let matad_result = evaluate_lane(
        vakint,
        &matad,
        canonical,
        tensor_reduced,
        original_form_path,
    );
    let rustred_result = evaluate_lane(
        vakint,
        &rustred,
        canonical,
        tensor_reduced,
        original_form_path,
    );
    let difference = (rustred_result - matad_result).together().expand();
    assert!(
        difference.is_zero(),
        "RustRed and MATAD raw master coefficients differ exactly: {difference}"
    );
}

#[allow(unused, clippy::too_many_arguments)]
fn compare_evaluations_impl(
    vakint_default_settings: VakintSettings,
    lanes: &[EvaluationTestLane],
    tensor_prepass: TensorPrepass,
    integra_view: AtomView,
    numerical_masses: HashMap<String, Float, ahash::RandomState>,
    numerical_external_momenta: HashMap<usize, Momentum, ahash::RandomState>,
    rel_threshold: f64,
    max_pull: f64,
    quiet: bool,
) {
    assert!(lanes.len() >= 2, "a comparison needs at least two lanes");
    let python_exe = vakint_default_settings.python_exe_path.clone();
    let lanes = adjusted_lanes(
        lanes,
        Some(quiet),
        rel_threshold * 1.0e-2,
        &numerical_masses,
        &numerical_external_momenta,
    );
    if lanes
        .iter()
        .any(|lane| evaluation_requires_pysecdec(&lane.evaluation_order))
        && !pysecdec_available(&python_exe)
    {
        eprintln!("Skipping test: PySecDec not available.");
        return;
    }

    let settings = VakintSettings {
        allow_unknown_integrals: false,
        use_dot_product_notation: true,
        mu_r_sq_symbol: "mursq".into(),
        integral_normalization_factor: LoopNormalizationFactor::pySecDec,
        evaluation_order: lanes[0].evaluation_order.clone(),
        ..vakint_default_settings
    };
    let original_form_path = settings.form_exe_path.clone();
    let mut vakint = get_vakint(settings);
    let (canonical, tensor_reduced) =
        prepare_comparison_inputs(&mut vakint, integra_view, tensor_prepass);
    assert_exact_raw_rustred_matad_peer(
        &mut vakint,
        &lanes,
        &canonical,
        &tensor_reduced,
        &original_form_path,
    );
    let mut evaluation_parameters = HashMap::default();
    let mass_parameter = if numerical_masses.contains_key("user_space::muv") {
        "user_space::muv"
    } else {
        "muvsq"
    };
    evaluation_parameters.insert(
        mass_parameter.to_owned(),
        numerical_masses
            .get(mass_parameter)
            .unwrap_or_else(|| panic!("{mass_parameter} not found in numerical_masses"))
            .clone(),
    );
    evaluation_parameters.insert(
        "mursq".to_owned(),
        numerical_masses
            .get("mursq")
            .unwrap_or_else(|| panic!("mursq not found in numerical_masses"))
            .clone(),
    );

    let mut results = Vec::with_capacity(lanes.len());
    for lane in &lanes {
        debug!("Evaluating integral with {}", lane.evaluation_order);
        let evaluated = evaluate_lane(
            &mut vakint,
            lane,
            &canonical,
            &tensor_reduced,
            &original_form_path,
        );
        let numerical = Vakint::full_numerical_evaluation(
            &vakint.settings,
            evaluated.as_view(),
            &evaluation_parameters,
            &HashMap::default(),
            Some(&numerical_external_momenta),
        )
        .unwrap_or_else(|error| {
            panic!(
                "numerical evaluation of {} failed: {error}",
                lane.evaluation_order
            )
        });
        results.push((lane.evaluation_order.clone(), numerical));
    }

    let (benchmark_order, (benchmark_central, benchmark_error)) = &results[0];
    for (tested_order, (tested_central, tested_error)) in &results[1..] {
        let combined_error = match (benchmark_error, tested_error) {
            (Some(benchmark), Some(tested)) => Some(benchmark.aggregate_errors(tested)),
            (Some(benchmark), None) => Some(benchmark.clone()),
            (None, Some(tested)) => Some(tested.clone()),
            (None, None) => None,
        };
        let (matches, message) = benchmark_central.does_approx_match(
            tested_central,
            combined_error.as_ref(),
            rel_threshold,
            max_pull,
        );
        if !matches || !quiet {
            println!("\n{}\n", "<><><><><>".green());
            println!("Benchmark {benchmark_order}:\n{benchmark_central}");
            if let Some(error) = benchmark_error {
                println!("Benchmark error:\n{error}");
            }
            println!("Tested {tested_order}:\n{tested_central}");
            if let Some(error) = tested_error {
                println!("Tested error:\n{error}");
            }
            println!("{message}");
            println!("\n{}\n", "<><><><><>".green());
        } else {
            info!("Benchmark {benchmark_order}:\n{benchmark_central}");
            info!("Tested {tested_order}:\n{tested_central}");
            info!("{message}");
        }
        assert!(
            matches,
            "{tested_order} does not match benchmark {benchmark_order}: {message}"
        );
    }
}

const MULTI_LANE_TEST_STACK_BYTES: usize = 32 * 1024 * 1024;

pub fn run_multi_lane_acceptance(work: impl FnOnce() + Send + 'static) {
    let handle = std::thread::Builder::new()
        .name("vakint-multi-lane-acceptance".to_owned())
        .stack_size(MULTI_LANE_TEST_STACK_BYTES)
        .spawn(work)
        .unwrap_or_else(|error| panic!("failed to spawn acceptance-test worker: {error}"));
    if let Err(payload) = handle.join() {
        std::panic::resume_unwind(payload);
    }
}

/// Compare multiple scalar-reduction peers on a bounded, test-local stack.
///
/// Symbolica's exact raw-expression normalization can recurse more deeply than
/// the Rust test harness's default worker stack.  Keeping this stack on one
/// short-lived acceptance worker avoids requiring a process-wide
/// `RUST_MIN_STACK` override or inflating every parallel test worker.
#[allow(unused, clippy::too_many_arguments)]
pub fn compare_evaluations(
    vakint_default_settings: VakintSettings,
    lanes: &[EvaluationTestLane],
    tensor_prepass: TensorPrepass,
    integra_view: AtomView,
    numerical_masses: HashMap<String, Float, ahash::RandomState>,
    numerical_external_momenta: HashMap<usize, Momentum, ahash::RandomState>,
    rel_threshold: f64,
    max_pull: f64,
    quiet: bool,
) {
    let lanes = lanes.to_vec();
    let integral = integra_view.to_owned();
    run_multi_lane_acceptance(move || {
        compare_evaluations_impl(
            vakint_default_settings,
            &lanes,
            tensor_prepass,
            integral.as_view(),
            numerical_masses,
            numerical_external_momenta,
            rel_threshold,
            max_pull,
            quiet,
        )
    });
}

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
    let lanes = [
        EvaluationTestLane::oracle(
            evaluation_orders.0.0.clone(),
            if evaluation_orders.0.1 {
                EvaluationTestInput::TensorReduced
            } else {
                EvaluationTestInput::Canonical
            },
        ),
        EvaluationTestLane::oracle(
            evaluation_orders.1.0.clone(),
            if evaluation_orders.1.1 {
                EvaluationTestInput::TensorReduced
            } else {
                EvaluationTestInput::Canonical
            },
        ),
    ];
    let tensor_prepass = if evaluation_orders.0.1 || evaluation_orders.1.1 {
        TensorPrepass::Form
    } else {
        TensorPrepass::Skip
    };
    compare_evaluations(
        vakint_default_settings,
        &lanes,
        tensor_prepass,
        integra_view,
        numerical_masses,
        numerical_external_momenta,
        rel_threshold,
        max_pull,
        quiet,
    );
}

#[allow(unused, clippy::too_many_arguments)]
fn compare_vakint_evaluations_vs_reference_impl(
    vakint_default_settings: VakintSettings,
    lanes: &[EvaluationTestLane],
    tensor_prepass: TensorPrepass,
    integra_view: AtomView,
    numerical_masses: HashMap<String, Float, ahash::RandomState>,
    numerical_external_momenta: HashMap<usize, Momentum, ahash::RandomState>,
    expected_output: Vec<(i64, (String, String))>,
    prec: u32,
    max_pull: f64,
) {
    assert!(!lanes.is_empty(), "a reference comparison needs a lane");
    let python_exe = vakint_default_settings.python_exe_path.clone();
    let lanes = adjusted_lanes(
        lanes,
        None,
        10.0_f64.powi(-(prec as i32)),
        &numerical_masses,
        &numerical_external_momenta,
    );
    if lanes
        .iter()
        .any(|lane| evaluation_requires_pysecdec(&lane.evaluation_order))
        && !pysecdec_available(&python_exe)
    {
        eprintln!("Skipping test: PySecDec not available.");
        return;
    }

    let settings = VakintSettings {
        allow_unknown_integrals: true,
        use_dot_product_notation: true,
        mu_r_sq_symbol: vakint_default_settings.mu_r_sq_symbol.clone(),
        run_time_decimal_precision: prec,
        evaluation_order: lanes[0].evaluation_order.clone(),
        ..vakint_default_settings
    };
    let original_form_path = settings.form_exe_path.clone();
    let mut vakint = get_vakint(settings);
    let (canonical, tensor_reduced) =
        prepare_comparison_inputs(&mut vakint, integra_view, tensor_prepass);
    assert_exact_raw_rustred_matad_peer(
        &mut vakint,
        &lanes,
        &canonical,
        &tensor_reduced,
        &original_form_path,
    );

    let binary_prec: u32 = ((prec.max(16) as f64) * LOG2_10).floor() as u32;
    let reference = NumericalEvaluationResult(
        expected_output
            .iter()
            .map(|(epsilon_power, (real, imaginary))| {
                (
                    *epsilon_power,
                    Complex::new(
                        Float::parse(real.as_str(), Some(binary_prec)).unwrap(),
                        Float::parse(imaginary.as_str(), Some(binary_prec)).unwrap(),
                    ),
                )
            })
            .collect(),
    );

    for lane in &lanes {
        let evaluated = evaluate_lane(
            &mut vakint,
            lane,
            &canonical,
            &tensor_reduced,
            &original_form_path,
        );
        let (result, error) = Vakint::full_numerical_evaluation(
            &vakint.settings,
            evaluated.as_view(),
            &numerical_masses,
            &HashMap::default(),
            Some(&numerical_external_momenta),
        )
        .unwrap_or_else(|evaluation_error| {
            panic!(
                "numerical evaluation of {} failed: {evaluation_error}",
                lane.evaluation_order
            )
        });
        let (matches, message) = result.does_approx_match(
            &reference,
            error.as_ref(),
            0.1_f64.powi((prec as i32) - 2),
            max_pull,
        );
        if !matches {
            println!("Result from {}:\n{}", lane.evaluation_order, result);
            if let Some(error) = &error {
                println!("Error:\n{error}");
            }
            println!("Reference:\n{reference}");
            println!("{message}");
        } else {
            debug!("Result from {}:\n{}", lane.evaluation_order, result);
            debug!("Reference:\n{reference}");
        }
        assert!(
            matches,
            "{} does not match the Vakint reference: {message}",
            lane.evaluation_order
        );
    }
}

/// Compare multiple scalar-reduction peers to one fixed reference on the same
/// bounded, test-local stack as [`compare_evaluations`].
#[allow(unused, clippy::too_many_arguments)]
pub fn compare_vakint_evaluations_vs_reference(
    vakint_default_settings: VakintSettings,
    lanes: &[EvaluationTestLane],
    tensor_prepass: TensorPrepass,
    integra_view: AtomView,
    numerical_masses: HashMap<String, Float, ahash::RandomState>,
    numerical_external_momenta: HashMap<usize, Momentum, ahash::RandomState>,
    expected_output: Vec<(i64, (String, String))>,
    prec: u32,
    max_pull: f64,
) {
    let lanes = lanes.to_vec();
    let integral = integra_view.to_owned();
    run_multi_lane_acceptance(move || {
        compare_vakint_evaluations_vs_reference_impl(
            vakint_default_settings,
            &lanes,
            tensor_prepass,
            integral.as_view(),
            numerical_masses,
            numerical_external_momenta,
            expected_output,
            prec,
            max_pull,
        )
    });
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
    compare_vakint_evaluations_vs_reference(
        vakint_default_settings,
        &[EvaluationTestLane::oracle(
            evaluation_order,
            EvaluationTestInput::TensorReduced,
        )],
        TensorPrepass::Form,
        integra_view,
        numerical_masses,
        numerical_external_momenta,
        expected_output,
        prec,
        max_pull,
    );
}
