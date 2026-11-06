//! Byte comparison of vacuum outputs against the frozen report of the base.
//!
//! The p-integral frontend edits vakint's shared paths (canonicalization,
//! matching, tensor reduction, `prepare_integral_specs`, `validate_settings`),
//! and every one of them must leave a vacuum integral alone. This test renders
//! each public step on vacuum inputs, on the backends of [`lanes`], and must
//! equal `tests/fixtures/vacuum_byte_identity.txt` byte for byte. The report
//! also holds the symbol table, which fixes the order of every printed term.
//!
//! That file was written by this same test compiled against the base commit,
//! with `VAKINT_VACUUM_BYTES_OUT=<path>`, which writes the report and compares
//! nothing. Regenerate it only for an intended vacuum change, from a tree
//! without the change under review.

use std::panic::{AssertUnwindSafe, catch_unwind};

use symbolica::{
    atom::{Atom, AtomCore},
    state::State,
};
use vakint::{
    EvaluationMethod, EvaluationOrder, RustRedEvaluationOptions, RustRedMasterMode,
    TensorReductionMethod, Vakint, VakintDependency, VakintError, VakintSettings, vakint_parse,
};

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/vacuum_byte_identity.txt"
);
const OUTPUT_VARIABLE: &str = "VAKINT_VACUUM_BYTES_OUT";
/// The lanes that run FORM start with this; the others must never need it.
const FORM_LANE: &str = "form/";
/// Part of the first column of every line of the symbol table.
const SYMBOL_SECTION: &str = "/symbols";
/// `form_exe_path` of the FORM-free lanes: a FORM call fails loudly.
const NO_FORM: &str = "/vacuum-byte-identity-must-not-use-form";

const TADPOLE: &str = "topo(prop(1,edge(1,1),k(1),muvsq,3))";
const SUNSET: &str = "topo(prop(1,edge(1,2),k(1),muvsq,2)*prop(2,edge(1,2),k(2),muvsq,1)*prop(3,edge(2,1),k(1)+k(2),muvsq,1))";
const TWO_MASSES: &str = "topo(prop(1,edge(1,2),k(1),m1sq,1)*prop(2,edge(2,1),k(1),m2sq,1))";

/// Short and long forms, tensor numerators of even and odd rank, sums, float
/// coefficients, unknown and massless integrals, malformed inputs.
fn inputs() -> Vec<String> {
    let mut inputs: Vec<String> = [
        "0",
        "topo(I1L(muvsq,3))",
        "dot(k(1),k(1))*topo(I2L(muvsq,2,1,1))",
        "topo(I3L(muvsq,2,1,1,1,1,1))",
        "dot(k(1),k(2))*topo(I3L(muvsq,1,1,1,2,1,1))",
        "topo(I1L(muvsq,2))+2*topo(I2L(muvsq,1,1,1))",
        "1.5`30*topo(I1L(muvsq,2))",
        "topo(prop(1,edge(1,1),k(1),0,1))",
        "topo(foo(1))",
        "dot(k(1),k(1))",
        "topo(prop(1,edge(1,2),k(1),muvsq,1))",
        "topo(prop(1,edge(1,2),k(1),muvsq))",
        "topo(prop(1,edge(1,1),k(1),muvsq,1/2))",
        "topo(prop(1,edge(1,1),k(1),muvsq,a))",
    ]
    .map(str::to_owned)
    .to_vec();
    for numerator in [
        "1",
        "k(1,1)*p(1,1)",
        "k(1,1)*k(1,2)*p(1,1)*p(1,2)",
        "k(1,1)*k(1,1)",
        "dot(k(1),p(1))",
        "dot(k(1),p(1))^2",
        "2.5`30*dot(k(1),p(1))^2",
    ] {
        inputs.push(format!("{numerator}*{TADPOLE}"));
    }
    inputs.push(SUNSET.to_owned());
    inputs.push(format!("dot(k(1),k(2))*{SUNSET}"));
    inputs.push(format!("k(1,1)*k(2,1)*{SUNSET}"));
    inputs.push(TWO_MASSES.to_owned());
    inputs.push(format!("dot(k(1),p(1))*{TWO_MASSES}"));
    inputs.push(format!(
        "topo(UNKNOWN({}))",
        &TWO_MASSES["topo(".len()..TWO_MASSES.len() - 1]
    ));
    inputs
}

fn settings(evaluation_order: EvaluationOrder, form_exe_path: Option<&str>) -> VakintSettings {
    let defaults = VakintSettings::default();
    VakintSettings {
        evaluation_order,
        form_exe_path: form_exe_path.map_or(defaults.form_exe_path.clone(), str::to_owned),
        // The default order and core count depend on the host.
        n_cores: 1,
        ..defaults
    }
}

fn rustred(master_mode: RustRedMasterMode) -> EvaluationOrder {
    EvaluationOrder(vec![EvaluationMethod::RustRed(RustRedEvaluationOptions {
        master_mode,
    })])
}

/// `EvaluationOrder::all()` without pySecDec, which depends on the host.
fn all_installed() -> EvaluationOrder {
    let mut order = EvaluationOrder::all();
    order
        .0
        .retain(|method| !matches!(method, EvaluationMethod::PySecDec(_)));
    order
}

fn lanes() -> Vec<(&'static str, VakintSettings)> {
    let retain = settings(rustred(RustRedMasterMode::Retain), Some(NO_FORM));
    let no_unknown = VakintSettings {
        allow_unknown_integrals: false,
        ..retain.clone()
    };
    let dot_notation = VakintSettings {
        use_dot_product_notation: true,
        ..retain.clone()
    };
    let verified = VakintSettings {
        verify_numerator_identification: true,
        ..retain.clone()
    };
    let alphaloop_tensors = VakintSettings {
        tensor_reduction_method: TensorReductionMethod::AlphaLoop,
        ..settings(EvaluationOrder::alphaloop_only(), None)
    };
    vec![
        ("rustred_retain", retain),
        ("rustred_no_unknown", no_unknown),
        ("rustred_dot_notation", dot_notation),
        ("rustred_verified", verified),
        (
            "rustred_numerical",
            settings(rustred(RustRedMasterMode::Numerical), Some(NO_FORM)),
        ),
        (
            "form/alphaloop",
            settings(EvaluationOrder::alphaloop_only(), None),
        ),
        ("form/alphaloop_tensors", alphaloop_tensors),
        (
            "form/matad",
            settings(EvaluationOrder::matad_only(None), None),
        ),
        ("form/all_installed", settings(all_installed(), None)),
    ]
}

/// One line of the report: values never contain a newline or a tab.
fn line(key: &[&str], value: &str) -> String {
    format!(
        "{}\t{}\n",
        key.join("\t"),
        value.replace('\n', "\\n").replace('\t', "\\t")
    )
}

type Step = Result<Result<Atom, VakintError>, String>;

/// A panic is an output too, and it must not hide the rest of the report.
fn step(run: impl FnOnce() -> Result<Atom, VakintError>) -> Step {
    catch_unwind(AssertUnwindSafe(run)).map_err(|payload| {
        payload
            .downcast_ref::<String>()
            .cloned()
            .or_else(|| payload.downcast_ref::<&str>().map(|text| text.to_string()))
            .unwrap_or_default()
    })
}

fn render(step: &Step) -> String {
    match step {
        Ok(Ok(value)) => format!("Ok {}", value.to_canonical_string()),
        // FORM's log carries the time of the run, and the node a graph error
        // names comes from a hash map.
        Ok(Err(VakintError::FormError(..))) => "Err FormError".to_owned(),
        Ok(Err(VakintError::MalformedGraph(..))) => "Err MalformedGraph".to_owned(),
        Ok(Err(error)) => format!("Err {error:?}"),
        Err(panic) => format!("Panic {panic}"),
    }
}

fn symbol_table(section: &str) -> String {
    let names = State::symbol_iter().map(|(_, name)| name);
    names.map(|name| line(&[section], name)).collect()
}

fn lane_report(vakint: &Vakint, lane: &str, settings: &VakintSettings) -> String {
    let mut report = line(
        &[lane, "validate_settings"],
        &format!("{:?}", vakint.validate_settings(settings)),
    );
    for input in inputs() {
        let atom = vakint_parse!(&input).unwrap();
        let view = atom.as_view();
        let mut record = |name: &str, value: String| {
            report.push_str(&line(&[lane, &input, name], &value));
        };
        let canonical_short = step(|| vakint.to_canonical(settings, view, true));
        record("canonical_short", render(&canonical_short));
        let canonical = step(|| vakint.to_canonical(settings, view, false));
        record("canonical", render(&canonical));
        let reduced_raw = step(|| vakint.tensor_reduce(settings, view));
        record("tensor_reduce_raw", render(&reduced_raw));
        let integral_raw = step(|| vakint.evaluate_integral(settings, view));
        record("evaluate_integral_raw", render(&integral_raw));
        record(
            "evaluate",
            render(&step(|| vakint.evaluate(settings, view))),
        );
        // The three public steps by hand, while each one gives a value.
        if let Ok(Ok(canonical)) = canonical {
            let reduced = step(|| vakint.tensor_reduce(settings, canonical.as_view()));
            record("tensor_reduce", render(&reduced));
            if let Ok(Ok(reduced)) = reduced {
                let value = step(|| vakint.evaluate_integral(settings, reduced.as_view()));
                record("evaluate_integral", render(&value));
            }
        }
    }
    report
}

fn assert_same_bytes(report: &str, expected: &str, with_form: bool) {
    // A later symbol extends the table without moving a term of any output, and
    // the numerical source symbols of the FORM lanes change order between
    // processes: the outputs are what is compared, the table is to read by hand.
    let free = |text: &&str| {
        let lane = text.split('\t').next().unwrap_or_default();
        !lane.contains(SYMBOL_SECTION) && (with_form || !lane.starts_with(FORM_LANE))
    };
    let actual = report.lines().filter(free).collect::<Vec<_>>();
    let expected = expected.lines().filter(free).collect::<Vec<_>>();
    for (index, (actual, expected)) in actual.iter().zip(&expected).enumerate() {
        assert_eq!(actual, expected, "line {} of the report", index + 1);
    }
    assert_eq!(actual.len(), expected.len(), "report length in lines");
}

#[test]
fn vacuum_outputs_match_the_base_byte_for_byte() {
    let vakint = Vakint::new().unwrap();
    let default_form = VakintSettings::default().form_exe_path;
    let with_form = VakintDependency::Forcer
        .version(&default_form, "python3")
        .is_ok();
    let mut report = symbol_table("free/symbols_after_new");
    for (lane, settings) in lanes() {
        if !lane.starts_with(FORM_LANE) {
            report.push_str(&lane_report(&vakint, lane, &settings));
        }
    }
    report.push_str(&symbol_table("free/symbols"));
    if with_form {
        for (lane, settings) in lanes() {
            if lane.starts_with(FORM_LANE) {
                report.push_str(&lane_report(&vakint, lane, &settings));
            }
        }
        report.push_str(&symbol_table("form/symbols"));
    } else {
        eprintln!("FORM 5 is unavailable: compared the FORM-free lanes only");
    }
    if let Some(path) = std::env::var_os(OUTPUT_VARIABLE) {
        std::fs::write(path, report).unwrap();
        return;
    }
    let expected = std::fs::read_to_string(FIXTURE).unwrap();
    assert_same_bytes(&report, &expected, with_form);
}
