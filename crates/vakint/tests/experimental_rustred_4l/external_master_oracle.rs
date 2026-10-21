//! Independently reduce with FORM and evaluate using the official coefficient
//! dataset. Its precision can exceed the historical hard-coded FORM constants.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::OnceLock;

use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::domains::float::Float;
use symbolica::domains::rational::Rational;
use symbolica::id::Replacement;
use vakint::{
    EvaluationOrder, FMFTOptions, TensorReductionMethod, VakintSettings, vakint_parse,
    vakint_symbol,
};

use super::{acceptance_inputs, test_utils};

struct ReferenceData {
    substitutions: Vec<(Atom, Atom)>,
    tails: Vec<Symbol>,
}

fn reference_data() -> &'static ReferenceData {
    static DATA: OnceLock<ReferenceData> = OnceLock::new();
    DATA.get_or_init(|| {
        let root = test_utils::external_data_path().join("4L");
        let manifest: serde_json::Value =
            serde_json::from_slice(&std::fs::read(root.join("manifest.json")).unwrap()).unwrap();
        let mut primitives = BTreeSet::new();
        let mut substitutions = Vec::new();
        let mut tails = Vec::new();
        let epsilon = vakint_parse!("ep").unwrap();
        let bits = VakintSettings {
            run_time_decimal_precision: 80,
            ..Default::default()
        }
        .get_binary_precision();
        for master in manifest["masters"].as_array().unwrap() {
            let primitive =
                vakint_parse!(master["exact_basis_expression"].as_str().unwrap()).unwrap();
            if !matches!(primitive.as_view(), AtomView::Var(_))
                || !primitives.insert(primitive.clone())
            {
                continue;
            }
            let read = |name: &str| -> BTreeMap<i64, String> {
                let path = master["files"][name]["path"].as_str().unwrap();
                let rows: Vec<(i64, String)> =
                    serde_json::from_slice(&std::fs::read(root.join(path)).unwrap()).unwrap();
                rows.into_iter().collect()
            };
            let exact = read("symbolic.json");
            let numerical = read("numerical.json");
            let mut polynomial = Atom::Zero;
            let mut last = i64::MIN;
            for metadata in master["coefficients"].as_array().unwrap() {
                let power = metadata["epsilon_power"].as_i64().unwrap();
                let value = if metadata["symbolic"].as_bool().unwrap() {
                    vakint_parse!(&exact[&power]).unwrap()
                } else {
                    assert!(metadata["significant_digits"].as_u64().unwrap() >= 32);
                    Atom::num(Float::parse(&numerical[&power], Some(bits)).unwrap())
                };
                polynomial += value * epsilon.pow(Atom::num(power));
                last = last.max(power);
            }
            let tail = vakint_symbol!(format!(
                "external_oracle_tail_{}",
                master["name"].as_str().unwrap()
            ));
            polynomial += Atom::var(tail) * epsilon.pow(Atom::num(last + 1));
            tails.push(tail);
            substitutions.push((primitive, polynomial));
        }
        ReferenceData {
            substitutions,
            tails,
        }
    })
}

pub(super) fn compare(case: &acceptance_inputs::AcceptanceInput, form: &str) {
    let engine = test_utils::external_vakint();
    let settings = VakintSettings {
        tensor_reduction_method: TensorReductionMethod::FeynKit,
        use_dot_product_notation: true,
        form_exe_path: form.to_owned(),
        evaluation_order: EvaluationOrder::fmft_only(Some(FMFTOptions {
            expand_masters: false,
            ..Default::default()
        })),
        ..case.settings.clone()
    };
    let canonical = engine
        .to_canonical(&settings, case.input.as_view(), false)
        .unwrap();
    let scalar = engine
        .tensor_reduce(&settings, canonical.as_view())
        .unwrap();
    let raw_form = engine
        .evaluate_integral(&settings, scalar.as_view())
        .unwrap();
    let data = reference_data();
    let epsilon = vakint_parse!(&settings.epsilon_symbol).unwrap();
    let substitutions: Vec<_> = data
        .substitutions
        .iter()
        .map(|(primitive, value)| {
            Replacement::new(
                primitive.clone(),
                value
                    .replace(vakint_parse!("ep").unwrap().to_pattern())
                    .with(epsilon.clone()),
            )
        })
        .collect();
    let reference = raw_form
        .replace_multiple(&substitutions)
        .replace(vakint_parse!("Gamma(x_)").unwrap().to_pattern())
        .with(
            symbolica::transcendental::gamma()
                .call(vakint_parse!("x_").unwrap())
                .to_pattern(),
        )
        .replace(vakint_parse!("EulerGamma").unwrap().to_pattern())
        .with(Atom::var(symbolica::transcendental::euler_gamma()))
        .series(
            vakint_symbol!(&settings.epsilon_symbol),
            Atom::Zero.as_view(),
            Rational::from(settings.number_of_terms_in_epsilon_expansion - 5),
        )
        .unwrap()
        .to_atom()
        .expand();
    for tail in &data.tails {
        assert!(
            !reference.contains_symbol(*tail),
            "{} needs an unavailable external master order",
            case.name
        );
    }
    let working = VakintSettings {
        run_time_decimal_precision: 80,
        ..settings.clone()
    };
    // The reduction above is independently performed by FORM. Its exact
    // coefficient formulas use the same authenticated scientific inputs and
    // tracked HPL/K5 boundary as the native reduction, never parameter overrides
    // of registered mathematical constants.
    let (reference, _) = engine
        .numerical_evaluation(
            &working,
            reference.as_view(),
            &case.parameters_at_reference_point(80),
            &Default::default(),
            Some(&case.externals_at_reference_point(80)),
        )
        .unwrap();
    let native_settings = VakintSettings {
        form_exe_path: "/external-master-oracle-must-not-use-form".into(),
        evaluation_order: EvaluationOrder::rustred_only(),
        ..settings
    };
    let actual = engine
        .evaluate_integral(&native_settings, scalar.as_view())
        .unwrap();
    let (actual, _) = engine
        .numerical_evaluation(
            &native_settings,
            actual.as_view(),
            &case.parameters_at_reference_point(32),
            &Default::default(),
            Some(&case.externals_at_reference_point(32)),
        )
        .unwrap();
    let (matches, detail) = actual.does_approx_match(&reference, None, 1e-30, 0.0);
    assert!(
        matches,
        "{}: independent FORM reduction with external coefficients: {detail}\nFORM: {reference}\nRustRed: {actual}",
        case.name
    );
}
