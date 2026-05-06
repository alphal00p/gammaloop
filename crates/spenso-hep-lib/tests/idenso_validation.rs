use ahash::{HashMap, HashMapExt};
use idenso::{
    color::CS,
    reference_cases::{ReferenceCase, ReferenceDomain, ReferenceValidation, reference_cases},
};
use spenso::{
    algebra::upgrading_arithmetic::FallibleSub,
    iterators::IteratableTensor,
    network::parsing::ParseSettings,
    shadowing::Concretize,
    structure::{
        TensorStructure,
        abstract_index::AbstractIndex,
        representation::{Minkowski, RepName},
    },
    symbolica_atom::ProjectorExpander,
    tensors::{
        data::{DenseTensor, SparseOrDense},
        parametric::atomcore::TensorAtomMaps,
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    coefficient::CoefficientView,
    domains::float::RealLike,
};

use crate::common::{HepAtomExt, NetExt, test_initialize};

mod common;

type Complex64 = symbolica::domains::float::Complex<f64>;

#[test]
fn enabled_idenso_reference_cases_validate_against_explicit_tensors() {
    test_initialize();

    let mut validated = Vec::new();
    for case in reference_cases() {
        if !case.validation.is_enabled() {
            continue;
        }

        assert_case_valid(case);
        validated.push(case.name);
    }

    assert!(
        !validated.is_empty(),
        "the shared idenso reference table has no enabled validation cases"
    );
}

#[test]
fn blocked_idenso_reference_cases_stay_documented() {
    let blocked = reference_cases()
        .iter()
        .filter_map(|case| match case.validation {
            ReferenceValidation::Enabled => None,
            ReferenceValidation::Blocked(reason) => Some((case.name, reason)),
        })
        .collect::<Vec<_>>();

    assert!(
        !blocked.is_empty(),
        "remove this test once every shared idenso reference case validates explicitly"
    );
}

#[track_caller]
fn assert_case_valid(case: &ReferenceCase) {
    let original = case.expression();
    let simplified = case.simplify(&original);
    let constants = constants_for(case.domain);

    assert_valid(original, simplified, &constants);
}

fn constants_for(domain: ReferenceDomain) -> HashMap<Atom, Complex64> {
    let mut constants = scalar_constants();
    if domain == ReferenceDomain::Dirac4 {
        insert_vector_values(&mut constants, "p", 1.);
        insert_vector_values(&mut constants, "q", 5.);
    }
    constants
}

fn scalar_constants() -> HashMap<Atom, Complex64> {
    let mut constants = HashMap::new();
    insert_scalar(&mut constants, CS.ca, 3.);
    insert_scalar(&mut constants, CS.cf, 4. / 3.);
    insert_scalar(&mut constants, CS.tr, 0.5);
    insert_scalar(&mut constants, CS.nc, 3.);
    insert_scalar(&mut constants, CS.na, 8.);
    insert_scalar_aliases(&mut constants, "m", 5.);
    insert_scalar_aliases(&mut constants, "c1", 7.);
    insert_scalar_aliases(&mut constants, "c2", 11.);
    constants
}

fn insert_scalar_aliases(constants: &mut HashMap<Atom, Complex64>, name: &str, value: f64) {
    insert_scalar(constants, symbolica::symbol!(name), value);
    insert_scalar(
        constants,
        symbolica::symbol!(format!("idenso::{name}")),
        value,
    );
}

fn insert_scalar(constants: &mut HashMap<Atom, Complex64>, symbol: Symbol, value: f64) {
    constants.insert(Atom::var(symbol), Complex64::new(value, 0.));
}

fn insert_vector_values(constants: &mut HashMap<Atom, Complex64>, name: &str, offset: f64) {
    let tensor: DenseTensor<Atom, _> =
        spenso::network::parsing::ShadowedStructure::<AbstractIndex>::from_iter(
            [Minkowski {}.new_slot(4, 1)],
            symbolica::symbol!(format!("spenso::{name}")),
            None,
        )
        .structure
        .to_shell()
        .concretize(None);

    for (index, component) in tensor.iter_flat() {
        constants.insert(
            component.clone(),
            Complex64::new(offset + usize::from(index) as f64, 0.),
        );
    }
}

#[track_caller]
fn assert_valid(original: Atom, simplified: Atom, constants: &HashMap<Atom, Complex64>) {
    if simplified.is_zero() {
        assert_hep_tensor_zero(evaluate(original, constants));
        return;
    }

    if original.is_zero() {
        assert_hep_tensor_zero(evaluate(simplified, constants));
        return;
    }

    let original_result = evaluate(original, constants);
    let simplified_result = evaluate(simplified, constants);
    let difference = original_result.sub_fallible(&simplified_result).unwrap();

    assert_hep_tensor_zero(difference);
}

#[track_caller]
fn assert_hep_tensor_zero(mut difference: spenso_hep_lib::HepTensor<AbstractIndex>) {
    difference.to_param();
    let difference = difference.try_into_parametric().unwrap();
    let zero = difference
        .zero_test(10, 0.01)
        .iter_flat()
        .fold(symbolica::id::ConditionResult::True, |acc, (_, value)| {
            acc & *value
        });

    match zero {
        symbolica::id::ConditionResult::True => {}
        symbolica::id::ConditionResult::False
            if difference
                .iter_flat()
                .all(|(_, value)| is_numerically_small(value, 1e-12)) => {}
        symbolica::id::ConditionResult::Inconclusive => {
            panic!("validation inconclusive")
        }
        symbolica::id::ConditionResult::False => {
            panic!("validation difference should be zero:\n{difference}")
        }
    }
}

fn is_numerically_small(value: AtomView, tolerance: f64) -> bool {
    let AtomView::Num(number) = value else {
        return false;
    };

    match number.get_coeff_view() {
        CoefficientView::Natural(re_num, re_den, im_num, im_den) => {
            let re = re_num as f64 / re_den as f64;
            let im = im_num as f64 / im_den as f64;
            re.hypot(im) <= tolerance
        }
        CoefficientView::Float(re, im) => {
            let re = re.to_float().to_f64();
            let im = im.to_float().to_f64();
            re.hypot(im) <= tolerance
        }
        _ => false,
    }
}

#[track_caller]
fn evaluate(
    expression: Atom,
    constants: &HashMap<Atom, Complex64>,
) -> spenso_hep_lib::HepTensor<AbstractIndex> {
    let expression = expression.expand_projectors().expand();
    let net = expression
        .parse_to_hep_net(&ParseSettings::default())
        .unwrap_or_else(|err| {
            panic!("failed to parse validation expression `{expression}`: {err}")
        });

    let mut result = net.execute_and_res().unwrap_or_else(|err| {
        panic!("failed to execute validation expression `{expression}`: {err}")
    });

    let function_map = HashMap::new();
    result.evaluate_complex(|c| c.into(), constants, &function_map);
    result.to_dense()
}
