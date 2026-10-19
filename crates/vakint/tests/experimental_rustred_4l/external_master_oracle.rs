//! Independently reduce with FORM and evaluate using the official coefficient
//! dataset. Its precision can exceed the historical hard-coded FORM constants.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::sync::OnceLock;

use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::domains::float::{Complex, Float};
use symbolica::domains::rational::Rational;
use symbolica::id::Replacement;
use vakint::{
    EvaluationOrder, FMFTOptions, TensorReductionMethod, Vakint, VakintSettings, vakint_parse,
    vakint_symbol,
};

use super::{acceptance_inputs, test_utils};

struct ReferenceData {
    substitutions: Vec<(Atom, Atom)>,
    tails: Vec<Symbol>,
    hpls: BTreeMap<String, Float>,
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
        let mut hpls = BTreeMap::new();
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
                value.visitor(&mut |part| {
                    let AtomView::Fun(function) = part else {
                        return true;
                    };
                    let imaginary = function.get_symbol() == vakint_symbol!("HI");
                    if !imaginary && function.get_symbol() != vakint_symbol!("HR") {
                        return true;
                    }
                    let arguments: Vec<_> = function.iter().collect();
                    let (point, indices) = arguments.split_last().unwrap();
                    let omega = vakint_parse!("(1+3^(1/2)*𝑖)/2").unwrap();
                    assert!(
                        *point == vakint_parse!("om1").unwrap().as_view()
                            || (point.to_owned() - omega).expand().is_zero()
                    );
                    let word: Vec<_> = indices
                        .iter()
                        .map(|index| index.to_canonical_string().parse::<i8>().unwrap())
                        .collect();
                    let value = symbolica::transcendental::sixth_root_hpl(&word, 1)
                        .unwrap()
                        .evaluate_with_prec(&HashMap::<Atom, Complex<Float>>::new(), bits)
                        .unwrap();
                    hpls.insert(
                        part.to_canonical_string(),
                        if imaginary { value.im } else { value.re },
                    );
                    false
                });
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
            hpls,
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
    let parameters: ahash::HashMap<_, _> = case
        .parameters
        .iter()
        .map(|(name, value)| ((*name).to_owned(), *value))
        .collect();
    let external: ahash::HashMap<_, _> = case.external_momenta.iter().copied().collect();
    let working = VakintSettings {
        run_time_decimal_precision: 80,
        ..settings.clone()
    };
    let mut reference_parameters = vakint::params_from_f64(&parameters, 80);
    reference_parameters.extend(
        data.hpls
            .iter()
            .map(|(name, value)| (name.clone(), value.clone())),
    );
    let (reference, _) = Vakint::full_numerical_evaluation(
        &working,
        reference.as_view(),
        &reference_parameters,
        &Default::default(),
        Some(&vakint::externals_from_f64(&external, 80)),
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
    let (actual, _) = Vakint::full_numerical_evaluation(
        &native_settings,
        actual.as_view(),
        &vakint::params_from_f64(&parameters, 32),
        &Default::default(),
        Some(&vakint::externals_from_f64(&external, 32)),
    )
    .unwrap();
    let (matches, detail) = actual.does_approx_match(&reference, None, 1e-30, 0.0);
    assert!(
        matches,
        "{}: independent FORM reduction with external coefficients: {detail}\nFORM: {reference}\nRustRed: {actual}",
        case.name
    );
}
