use std::{
    collections::{BTreeMap, BTreeSet},
    sync::{Arc, Mutex},
};

use ahash::{HashMap, HashMapExt};
use idenso::{
    color::CS,
    epsilon::EPSILON_SYMBOL,
    reference_cases::{ReferenceCase, reference_cases},
};
use spenso::{
    algebra::complex::RealOrComplexRef,
    algebra::upgrading_arithmetic::FallibleSub,
    iterators::IteratableTensor,
    network::{parsing::ParseSettings, tags::SPENSO_TAG},
    shadowing::ProjectorExpander,
    structure::{
        abstract_index::{AIND_SYMBOLS, AbstractIndex},
        representation::RepName,
        slot::{DualSlotTo, IsAbstractSlot},
    },
    sym,
    tensors::{
        data::SparseOrDense,
        parametric::{AtomViewOrConcrete, atomcore::TensorAtomMaps},
    },
    trace,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    coefficient::CoefficientView,
    domains::float::RealLike,
    evaluate::EvaluationError,
    function,
};
use tabled::{Table, Tabled, settings::Style};

use crate::common::{HepAtomExt, NetExt, test_initialize};
use idenso::representations::{ColorAdjoint, ColorFundamental};

mod common;

type Complex64 = symbolica::domains::float::Complex<f64>;

#[derive(Tabled)]
struct ReferenceValidationRow {
    #[tabled(rename = "case")]
    case_name: &'static str,
    #[tabled(rename = "domain")]
    domain: String,
    #[tabled(rename = "evaluations")]
    evaluations: usize,
    #[tabled(rename = "rank1 functions")]
    rank_one_functions: usize,
    #[tabled(rename = "evidence")]
    evidence: String,
    #[tabled(rename = "status")]
    status: &'static str,
}

#[test]
fn enabled_idenso_reference_cases_validate_against_explicit_tensors() {
    test_initialize();

    let mut rows = Vec::new();
    for case in reference_cases() {
        if !case.validation.is_enabled() {
            continue;
        }

        rows.push(assert_case_valid(case));
    }

    assert!(
        !rows.is_empty(),
        "the shared idenso reference table has no enabled validation cases"
    );

    let mut table = Table::new(rows);
    table.with(Style::rounded());
    println!("validated enabled idenso reference cases\n{table}");
}

#[test]
fn color_trace_structure_contraction_matches_commutator() {
    test_initialize();

    let cof = ColorFundamental {}.new_rep(3);
    let coad = ColorAdjoint {}.new_rep(8);
    let [a, b, c, d, x] = ["a", "b", "c", "d", "x"]
        .map(|name| coad.slot::<AbstractIndex, _>(symbolica::symbol!(name)));

    let lhs = Atom::i()
        * idenso::color_f!(c, d, x)
        * trace!(
            &cof,
            idenso::color_t!(a),
            idenso::color_t!(b),
            idenso::color_t!(x)
        );
    let rhs = trace!(
        &cof,
        idenso::color_t!(a),
        idenso::color_t!(b),
        idenso::color_t!(c),
        idenso::color_t!(d)
    ) - trace!(
        &cof,
        idenso::color_t!(a),
        idenso::color_t!(b),
        idenso::color_t!(d),
        idenso::color_t!(c)
    );

    assert_valid(lhs, rhs, &scalar_constants());
}

#[test]
fn color_middle_trace_structure_contraction_matches_commutator() {
    test_initialize();

    let cof = ColorFundamental {}.new_rep(3);
    let coad = ColorAdjoint {}.new_rep(8);
    let [a, b, c, d, x] = ["a", "b", "c", "d", "x"]
        .map(|name| coad.slot::<AbstractIndex, _>(symbolica::symbol!(name)));

    let lhs = Atom::i()
        * idenso::color_f!(c, d, x)
        * trace!(
            &cof,
            idenso::color_t!(a),
            idenso::color_t!(x),
            idenso::color_t!(b)
        );
    let rhs = trace!(
        &cof,
        idenso::color_t!(a),
        idenso::color_t!(c),
        idenso::color_t!(d),
        idenso::color_t!(b)
    ) - trace!(
        &cof,
        idenso::color_t!(a),
        idenso::color_t!(d),
        idenso::color_t!(c),
        idenso::color_t!(b)
    );

    assert_valid(lhs, rhs, &scalar_constants());
}

#[test]
fn color_explicit_middle_generator_contraction_matches_commutator() {
    test_initialize();

    let cof = ColorFundamental {}.new_rep(3);
    let coad = ColorAdjoint {}.new_rep(8);
    let [a, b, c, d, x] = ["a", "b", "c", "d", "x"]
        .map(|name| coad.slot::<AbstractIndex, _>(symbolica::symbol!(name)));
    let [i, j, k] =
        ["i", "j", "k"].map(|name| cof.slot::<AbstractIndex, _>(symbolica::symbol!(name)));

    let lhs = Atom::i()
        * idenso::color_f!(c, d, x)
        * function!(CS.t, a.to_atom(), i.to_atom(), j.dual().to_atom())
        * function!(CS.t, x.to_atom(), j.to_atom(), k.dual().to_atom())
        * function!(CS.t, b.to_atom(), k.to_atom(), i.dual().to_atom());
    let rhs = trace!(
        &cof,
        idenso::color_t!(a),
        idenso::color_t!(c),
        idenso::color_t!(d),
        idenso::color_t!(b)
    ) - trace!(
        &cof,
        idenso::color_t!(a),
        idenso::color_t!(d),
        idenso::color_t!(c),
        idenso::color_t!(b)
    );

    assert_valid(lhs, rhs, &scalar_constants());
}

#[test]
fn color_trace_materialization_is_cyclic() {
    test_initialize();

    let cof = ColorFundamental {}.new_rep(3);
    let coad = ColorAdjoint {}.new_rep(8);
    let [a, b, c] =
        ["a", "b", "c"].map(|name| coad.slot::<AbstractIndex, _>(symbolica::symbol!(name)));

    let lhs = trace!(
        &cof,
        idenso::color_t!(a),
        idenso::color_t!(b),
        idenso::color_t!(c)
    );
    let rhs = trace!(
        &cof,
        idenso::color_t!(b),
        idenso::color_t!(c),
        idenso::color_t!(a)
    );

    assert_valid(lhs, rhs, &scalar_constants());
}

#[test]
fn color_three_trace_antisymmetric_part_matches_structure_constant() {
    test_initialize();

    let cof = ColorFundamental {}.new_rep(3);
    let coad = ColorAdjoint {}.new_rep(8);
    let [a, b, c] =
        ["a", "b", "c"].map(|name| coad.slot::<AbstractIndex, _>(symbolica::symbol!(name)));

    let lhs = trace!(
        &cof,
        idenso::color_t!(a),
        idenso::color_t!(b),
        idenso::color_t!(c)
    ) - trace!(
        &cof,
        idenso::color_t!(a),
        idenso::color_t!(c),
        idenso::color_t!(b)
    );
    let rhs = Atom::i() * Atom::var(CS.tr) * idenso::color_f!(a, b, c);

    assert_valid(lhs, rhs, &scalar_constants());
}

#[test]
fn color_explicit_three_trace_antisymmetric_part_matches_structure_constant() {
    test_initialize();

    let cof = ColorFundamental {}.new_rep(3);
    let coad = ColorAdjoint {}.new_rep(8);
    let [a, b, c] =
        ["a", "b", "c"].map(|name| coad.slot::<AbstractIndex, _>(symbolica::symbol!(name)));
    let [i, j, k] =
        ["i", "j", "k"].map(|name| cof.slot::<AbstractIndex, _>(symbolica::symbol!(name)));

    let abc = function!(CS.t, a.to_atom(), i.to_atom(), j.dual().to_atom())
        * function!(CS.t, b.to_atom(), j.to_atom(), k.dual().to_atom())
        * function!(CS.t, c.to_atom(), k.to_atom(), i.dual().to_atom());
    let acb = function!(CS.t, a.to_atom(), i.to_atom(), j.dual().to_atom())
        * function!(CS.t, c.to_atom(), j.to_atom(), k.dual().to_atom())
        * function!(CS.t, b.to_atom(), k.to_atom(), i.dual().to_atom());
    let rhs = Atom::i() * Atom::var(CS.tr) * idenso::color_f!(a, b, c);

    assert_valid(abc - acb, rhs, &scalar_constants());
}

#[test]
fn color_symmetric_trace_structure_contraction_matches_commutator_average() {
    test_initialize();

    let cof = ColorFundamental {}.new_rep(3);
    let coad = ColorAdjoint {}.new_rep(8);
    let [a, b, c, d, x] = ["a", "b", "c", "d", "x"]
        .map(|name| coad.slot::<AbstractIndex, _>(symbolica::symbol!(name)));

    let lhs = Atom::i()
        * idenso::color_f!(c, d, x)
        * trace!(
            &cof,
            sym!(
                idenso::color_t!(a),
                idenso::color_t!(b),
                idenso::color_t!(x)
            )
        );
    let rhs = (trace!(
        &cof,
        idenso::color_t!(a),
        idenso::color_t!(b),
        idenso::color_t!(c),
        idenso::color_t!(d)
    ) - trace!(
        &cof,
        idenso::color_t!(a),
        idenso::color_t!(b),
        idenso::color_t!(d),
        idenso::color_t!(c)
    ) + trace!(
        &cof,
        idenso::color_t!(a),
        idenso::color_t!(c),
        idenso::color_t!(d),
        idenso::color_t!(b)
    ) - trace!(
        &cof,
        idenso::color_t!(a),
        idenso::color_t!(d),
        idenso::color_t!(c),
        idenso::color_t!(b)
    )) / 2;

    assert_valid(lhs, rhs, &scalar_constants());
}

#[test]
fn color_trace_materialization_matches_explicit_generator_product() {
    test_initialize();

    let cof = ColorFundamental {}.new_rep(3);
    let coad = ColorAdjoint {}.new_rep(8);
    let [a, b, c, d] =
        ["a", "b", "c", "d"].map(|name| coad.slot::<AbstractIndex, _>(symbolica::symbol!(name)));
    let [i, j, k, l] =
        ["i", "j", "k", "l"].map(|name| cof.slot::<AbstractIndex, _>(symbolica::symbol!(name)));

    let lhs = trace!(
        &cof,
        idenso::color_t!(a),
        idenso::color_t!(b),
        idenso::color_t!(c),
        idenso::color_t!(d)
    );
    let rhs = function!(CS.t, a.to_atom(), i.to_atom(), j.dual().to_atom())
        * function!(CS.t, b.to_atom(), j.to_atom(), k.dual().to_atom())
        * function!(CS.t, c.to_atom(), k.to_atom(), l.dual().to_atom())
        * function!(CS.t, d.to_atom(), l.to_atom(), i.dual().to_atom());

    assert_valid(lhs, rhs, &scalar_constants());
}

#[test]
fn color_permuted_trace_materialization_matches_explicit_generator_product() {
    test_initialize();

    let cof = ColorFundamental {}.new_rep(3);
    let coad = ColorAdjoint {}.new_rep(8);
    let [a, b, c, d] =
        ["a", "b", "c", "d"].map(|name| coad.slot::<AbstractIndex, _>(symbolica::symbol!(name)));
    let [i, j, k, l] =
        ["i", "j", "k", "l"].map(|name| cof.slot::<AbstractIndex, _>(symbolica::symbol!(name)));

    let lhs = trace!(
        &cof,
        idenso::color_t!(a),
        idenso::color_t!(c),
        idenso::color_t!(d),
        idenso::color_t!(b)
    );
    let rhs = function!(CS.t, a.to_atom(), i.to_atom(), j.dual().to_atom())
        * function!(CS.t, c.to_atom(), j.to_atom(), k.dual().to_atom())
        * function!(CS.t, d.to_atom(), k.to_atom(), l.dual().to_atom())
        * function!(CS.t, b.to_atom(), l.to_atom(), i.dual().to_atom());

    assert_valid(lhs, rhs, &scalar_constants());
}

#[test]
fn color_middle_trace_materialization_matches_explicit_generator_product() {
    test_initialize();

    let cof = ColorFundamental {}.new_rep(3);
    let coad = ColorAdjoint {}.new_rep(8);
    let [a, b, x] =
        ["a", "b", "x"].map(|name| coad.slot::<AbstractIndex, _>(symbolica::symbol!(name)));
    let [i, j, k] =
        ["i", "j", "k"].map(|name| cof.slot::<AbstractIndex, _>(symbolica::symbol!(name)));

    let lhs = trace!(
        &cof,
        idenso::color_t!(a),
        idenso::color_t!(x),
        idenso::color_t!(b)
    );
    let rhs = function!(CS.t, a.to_atom(), i.to_atom(), j.dual().to_atom())
        * function!(CS.t, x.to_atom(), j.to_atom(), k.dual().to_atom())
        * function!(CS.t, b.to_atom(), k.to_atom(), i.dual().to_atom());

    assert_valid(lhs, rhs, &scalar_constants());
}

#[test]
fn color_structure_constant_matches_generator_commutator() {
    test_initialize();

    let cof = ColorFundamental {}.new_rep(3);
    let coad = ColorAdjoint {}.new_rep(8);
    let [a, b, c] =
        ["a", "b", "c"].map(|name| coad.slot::<AbstractIndex, _>(symbolica::symbol!(name)));
    let [i, j, k] =
        ["i", "j", "k"].map(|name| cof.slot::<AbstractIndex, _>(symbolica::symbol!(name)));

    let lhs = Atom::i()
        * idenso::color_f!(a, b, c)
        * function!(CS.t, c.to_atom(), i.to_atom(), j.dual().to_atom());
    let rhs = function!(CS.t, a.to_atom(), i.to_atom(), k.dual().to_atom())
        * function!(CS.t, b.to_atom(), k.to_atom(), j.dual().to_atom())
        - function!(CS.t, b.to_atom(), i.to_atom(), k.dual().to_atom())
            * function!(CS.t, a.to_atom(), k.to_atom(), j.dual().to_atom());
    assert_valid(lhs, rhs, &scalar_constants());
}

#[test]
fn color_four_generator_symmetric_trace_matches_permutation_average() {
    test_initialize();

    let cof = ColorFundamental {}.new_rep(3);
    let coad = ColorAdjoint {}.new_rep(8);
    let [a, b, c, d] =
        ["a", "b", "c", "d"].map(|name| coad.slot::<AbstractIndex, _>(symbolica::symbol!(name)));
    let factors = [a, b, c, d].map(|slot| idenso::color_t!(slot));
    let lhs = trace!(&cof, sym!(; factors.clone()));
    let mut rhs = Atom::Zero;
    for i in 0..4 {
        for j in 0..4 {
            for k in 0..4 {
                for l in 0..4 {
                    if [j, k, l].contains(&i) || k == j || l == j || l == k {
                        continue;
                    }
                    rhs += trace!(
                        &cof,
                        factors[i].clone(),
                        factors[j].clone(),
                        factors[k].clone(),
                        factors[l].clone()
                    );
                }
            }
        }
    }
    rhs /= Atom::num(24);

    assert_valid(lhs, rhs, &scalar_constants());
}

#[test]
fn color_four_generator_trace_decomposition_matches_direct_trace() {
    test_initialize();

    let cof = ColorFundamental {}.new_rep(3);
    let coad = ColorAdjoint {}.new_rep(8);
    let [a, b, c, d, x] = ["a", "b", "c", "d", "x"]
        .map(|name| coad.slot::<AbstractIndex, _>(symbolica::symbol!(name)));

    let lhs = trace!(
        &cof,
        idenso::color_t!(a),
        idenso::color_t!(b),
        idenso::color_t!(c),
        idenso::color_t!(d)
    );
    let rhs = trace!(
        &cof,
        sym!(
            idenso::color_t!(a),
            idenso::color_t!(b),
            idenso::color_t!(c),
            idenso::color_t!(d)
        )
    ) + Atom::i() / Atom::num(2)
        * idenso::color_f!(a, b, x)
        * trace!(
            &cof,
            sym!(
                idenso::color_t!(c),
                idenso::color_t!(d),
                idenso::color_t!(x)
            )
        )
        + Atom::i() / Atom::num(2)
            * idenso::color_f!(c, d, x)
            * trace!(
                &cof,
                sym!(
                    idenso::color_t!(a),
                    idenso::color_t!(b),
                    idenso::color_t!(x)
                )
            )
        - Atom::var(CS.tr) / Atom::num(6) * idenso::color_f!(a, c, x) * idenso::color_f!(b, d, x)
        + Atom::var(CS.tr) / Atom::num(3) * idenso::color_f!(a, d, x) * idenso::color_f!(b, c, x);

    assert_valid(lhs, rhs, &scalar_constants());
}

#[test]
fn color_single_contracted_structure_constants_satisfy_jacobi() {
    test_initialize();

    let coad = ColorAdjoint {}.new_rep(8);
    let [a, b, c, d, x] = ["a", "b", "c", "d", "x"]
        .map(|name| coad.slot::<AbstractIndex, _>(symbolica::symbol!(name)));
    let expr = idenso::color_f!(a, b, x) * idenso::color_f!(c, d, x)
        - idenso::color_f!(a, c, x) * idenso::color_f!(b, d, x)
        + idenso::color_f!(a, d, x) * idenso::color_f!(b, c, x);

    assert_valid(expr, Atom::Zero, &scalar_constants());
}

#[track_caller]
fn assert_case_valid(case: &ReferenceCase) -> ReferenceValidationRow {
    let original = case.expression();
    let simplified = case.simplify(&original);
    let constants = scalar_constants();

    let validation = assert_valid_with_label(original, simplified, &constants, case.name);
    assert!(
        validation.evaluations > 0,
        "validation case `{}` did not evaluate any expression",
        case.name
    );

    ReferenceValidationRow {
        case_name: case.name,
        domain: format!("{:?}", case.domain),
        evaluations: validation.evaluations,
        rank_one_functions: validation.rank_one_functions,
        evidence: validation.evidence_for_table(),
        status: "ok",
    }
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

#[track_caller]
fn assert_valid(
    original: Atom,
    simplified: Atom,
    constants: &HashMap<Atom, Complex64>,
) -> ValidationOutcome {
    let label = format!("`{original}` should equal `{simplified}`");
    assert_valid_with_label(original, simplified, constants, &label)
}

#[track_caller]
fn assert_valid_with_label(
    original: Atom,
    simplified: Atom,
    constants: &HashMap<Atom, Complex64>,
    label: &str,
) -> ValidationOutcome {
    if simplified.is_zero() {
        let evaluation = evaluate(original, constants);
        let outcome = ValidationOutcome::from_evaluation("original", &evaluation);
        assert_hep_tensor_zero(evaluation.tensor, label);
        return outcome;
    }

    if original.is_zero() {
        let evaluation = evaluate(simplified, constants);
        let outcome = ValidationOutcome::from_evaluation("simplified", &evaluation);
        assert_hep_tensor_zero(evaluation.tensor, label);
        return outcome;
    }

    let original_result = evaluate(original, constants);
    let simplified_result = evaluate(simplified, constants);
    original_result.assert_nontrivial(label, "original");
    simplified_result.assert_nontrivial(label, "simplified");
    let mut outcome = ValidationOutcome::default();
    outcome.add_evaluation("original", &original_result);
    outcome.add_evaluation("simplified", &simplified_result);
    let difference = original_result
        .tensor
        .sub_fallible(&simplified_result.tensor)
        .unwrap();

    assert_hep_tensor_zero(difference, label);
    outcome
}

#[track_caller]
fn assert_hep_tensor_zero(mut difference: spenso_hep_lib::HepTensor<AbstractIndex>, label: &str) {
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
            panic!(
                "validation difference should be zero for `{label}`:\n{}",
                validation_difference_summary(&difference)
            )
        }
    }
}

fn validation_difference_summary<T>(difference: &T) -> String
where
    T: IteratableTensor,
    for<'a> T::Data<'a>: Into<AtomView<'a>>,
{
    let mut count = 0usize;
    let mut shown = Vec::new();

    for (index, value) in difference.iter_flat() {
        let value = value.into();
        if is_numerically_small(value, 1e-12) {
            continue;
        }

        count += 1;
        if shown.len() < 16 {
            let expanded = difference
                .expanded_index(index)
                .map(|index| format!("{index:?}"))
                .unwrap_or_else(|_| format!("flat {}", usize::from(index)));
            shown.push(format!("{expanded}: {value}"));
        }
    }

    format!(
        "{count} entries exceed tolerance 1e-12; first {}:\n{}",
        shown.len(),
        shown.join("\n")
    )
}

fn is_numerically_small(value: AtomView, tolerance: f64) -> bool {
    atom_magnitude(value).is_some_and(|magnitude| magnitude <= tolerance)
}

#[track_caller]
fn evaluate(expression: Atom, constants: &HashMap<Atom, Complex64>) -> ValidationEvaluation {
    let expression = expression.expand_projectors().expand();

    evaluate_term(expression, constants)
}

fn evaluate_term(expression: Atom, constants: &HashMap<Atom, Complex64>) -> ValidationEvaluation {
    let validation_functions = ValidationFunctions::for_expression(expression.as_view());
    let rank_one_functions = validation_functions.rank_one_functions;
    let samples = validation_functions.samples.clone();
    let net = expression
        .parse_to_hep_net(&ParseSettings::default())
        .unwrap_or_else(|err| {
            panic!("failed to parse validation expression `{expression}`: {err}")
        });

    let mut result = net.execute_and_res().unwrap_or_else(|err| {
        panic!("failed to execute validation expression `{expression}`: {err}")
    });

    let mut evaluation_constants = constants.clone();
    validation_functions.insert_values_for_tensor(&result, &mut evaluation_constants);
    result.evaluate_complex(&evaluation_constants);
    let tensor = result.to_dense();
    let rank_one_samples = {
        let mut samples = samples.lock().unwrap();
        RankOneSampleSummary::from_samples(std::mem::take(&mut *samples))
    };
    rank_one_samples.assert_nontrivial(rank_one_functions, &expression);
    ValidationEvaluation {
        value_summary: TensorValueSummary::from_tensor(&tensor),
        rank_one_samples,
        tensor,
        rank_one_functions,
    }
}

#[derive(Default)]
struct ValidationOutcome {
    evaluations: usize,
    rank_one_functions: usize,
    evidence: Vec<String>,
}

impl ValidationOutcome {
    const MAX_TABLE_LINES: usize = 8;

    fn from_evaluation(role: &str, evaluation: &ValidationEvaluation) -> Self {
        let mut outcome = Self::default();
        outcome.add_evaluation(role, evaluation);
        outcome
    }

    fn add_evaluation(&mut self, role: &str, evaluation: &ValidationEvaluation) {
        self.evaluations += 1;
        self.rank_one_functions += evaluation.rank_one_functions;
        self.evidence.push(format!(
            "{role} tensor = {}",
            evaluation.value_summary.format()
        ));
        self.evidence
            .extend(evaluation.rank_one_samples.table_lines(
                role,
                Self::MAX_TABLE_LINES.saturating_sub(self.evidence.len()),
            ));
    }

    fn evidence_for_table(&self) -> String {
        if self.evidence.is_empty() {
            return "-".to_owned();
        }

        let mut evidence = self
            .evidence
            .iter()
            .take(Self::MAX_TABLE_LINES)
            .cloned()
            .collect::<Vec<_>>();
        if self.evidence.len() > Self::MAX_TABLE_LINES {
            evidence.push(format!(
                "... {} more evidence lines",
                self.evidence.len() - Self::MAX_TABLE_LINES
            ));
        }
        evidence.join("\n")
    }
}

struct ValidationEvaluation {
    tensor: spenso_hep_lib::HepTensor<AbstractIndex>,
    value_summary: TensorValueSummary,
    rank_one_samples: RankOneSampleSummary,
    rank_one_functions: usize,
}

impl ValidationEvaluation {
    #[track_caller]
    fn assert_nontrivial(&self, label: &str, role: &str) {
        self.value_summary.assert_nontrivial(label, role);
    }
}

struct EvaluationSample {
    symbol: String,
    args: Vec<ComplexSampleKey>,
    value: ComplexSampleKey,
    magnitude: f64,
    is_complex: bool,
}

impl EvaluationSample {
    fn new(symbol: Symbol, args: &[Complex64], value: Complex64) -> Self {
        Self {
            symbol: symbol.to_string(),
            args: args.iter().map(Self::complex_key).collect(),
            value: Self::complex_key(&value),
            magnitude: value.re.hypot(value.im),
            is_complex: value.im.abs() > TensorValueSummary::TOLERANCE,
        }
    }

    fn complex_key(value: &Complex64) -> ComplexSampleKey {
        (value.re.to_bits(), value.im.to_bits())
    }
}

type ComplexSampleKey = (u64, u64);

struct RankOneSampleSummary {
    symbols: BTreeMap<String, RankOneSymbolSummary>,
}

impl RankOneSampleSummary {
    fn from_samples(samples: Vec<EvaluationSample>) -> Self {
        let mut summary = Self {
            symbols: BTreeMap::new(),
        };
        for sample in samples {
            summary
                .symbols
                .entry(sample.symbol.clone())
                .or_default()
                .record(sample);
        }
        summary
    }

    #[track_caller]
    fn assert_nontrivial(&self, expected_symbols: usize, expression: &Atom) {
        if expected_symbols == 0 {
            return;
        }

        assert_eq!(
            self.symbols.len(),
            expected_symbols,
            "rank-1 validation for `{expression}` found {expected_symbols} rank-1 function(s), but only evaluated {}: {}",
            self.symbols.len(),
            self.symbols.keys().cloned().collect::<Vec<_>>().join(", ")
        );

        for (symbol, summary) in &self.symbols {
            summary.assert_nontrivial(symbol, expression);
        }
    }

    fn table_lines(&self, role: &str, limit: usize) -> Vec<String> {
        self.symbols
            .iter()
            .take(limit)
            .map(|(symbol, summary)| format!("{role} {symbol}: {}", summary.format_for_table()))
            .collect()
    }
}

#[derive(Default)]
struct RankOneSymbolSummary {
    calls: usize,
    nonzero_values: usize,
    complex_values: usize,
    unique_args: BTreeSet<Vec<ComplexSampleKey>>,
    unique_values: BTreeSet<ComplexSampleKey>,
}

impl RankOneSymbolSummary {
    fn record(&mut self, sample: EvaluationSample) {
        self.calls += 1;
        if sample.magnitude > TensorValueSummary::TOLERANCE {
            self.nonzero_values += 1;
        }
        if sample.is_complex {
            self.complex_values += 1;
        }
        self.unique_args.insert(sample.args);
        self.unique_values.insert(sample.value);
    }

    #[track_caller]
    fn assert_nontrivial(&self, symbol: &str, expression: &Atom) {
        assert!(
            self.calls > 0,
            "rank-1 function `{symbol}` was discovered in `{expression}` but never evaluated"
        );
        assert!(
            self.nonzero_values > 0,
            "rank-1 function `{symbol}` only produced zero values in `{expression}`"
        );
        assert!(
            self.complex_values > 0,
            "rank-1 function `{symbol}` only produced real values in `{expression}`"
        );
        if self.unique_args.len() > 1 {
            assert!(
                self.unique_values.len() > 1,
                "rank-1 function `{symbol}` produced a constant value across {} sampled components in `{expression}`",
                self.unique_args.len()
            );
        }
    }

    fn format_for_table(&self) -> String {
        format!(
            "calls={}, nonzero={}, complex={}, distinct_args={}, distinct_values={}",
            self.calls,
            self.nonzero_values,
            self.complex_values,
            self.unique_args.len(),
            self.unique_values.len()
        )
    }
}

type ValidationFunction = Box<dyn Fn(&[Complex64]) -> Complex64 + Send + Sync>;

struct ValidationFunctions {
    functions: HashMap<Symbol, ValidationFunction>,
    samples: Arc<Mutex<Vec<EvaluationSample>>>,
    rank_one_functions: usize,
}

impl ValidationFunctions {
    const SEED: u64 = 0x9e37_79b9_7f4a_7c15;

    fn for_expression(expression: AtomView<'_>) -> Self {
        let mut validation = Self {
            functions: HashMap::new(),
            samples: Arc::new(Mutex::new(Vec::new())),
            rank_one_functions: 0,
        };
        validation.insert_builtin_functions();
        validation.insert_rank_one_tensor_functions(expression);
        validation
    }

    fn insert_builtin_functions(&mut self) {
        self.functions.insert(
            AIND_SYMBOLS.cind,
            Box::new(|args| Complex64::new(encode_cind(args) as f64, 0.)),
        );
        self.functions.insert(
            *EPSILON_SYMBOL,
            Box::new(|args| Complex64::new(0., -levi_civita(args))),
        );
    }

    fn insert_rank_one_tensor_functions(&mut self, expression: AtomView<'_>) {
        match expression {
            AtomView::Fun(fun) => {
                let symbol = fun.get_symbol();
                if symbol.has_tag(&SPENSO_TAG.rank1) && !self.functions.contains_key(&symbol) {
                    self.rank_one_functions += 1;
                    let samples = self.samples.clone();
                    self.functions.insert(
                        symbol,
                        Box::new(move |args| {
                            let value = Self::seeded_tensor_component(symbol, args);
                            samples
                                .lock()
                                .unwrap()
                                .push(EvaluationSample::new(symbol, args, value));
                            value
                        }),
                    );
                }
                for arg in fun.iter() {
                    self.insert_rank_one_tensor_functions(arg);
                }
            }
            AtomView::Add(add) => {
                for arg in add.iter() {
                    self.insert_rank_one_tensor_functions(arg);
                }
            }
            AtomView::Mul(mul) => {
                for arg in mul.iter() {
                    self.insert_rank_one_tensor_functions(arg);
                }
            }
            AtomView::Pow(pow) => {
                let (base, exponent) = pow.get_base_exp();
                self.insert_rank_one_tensor_functions(base);
                self.insert_rank_one_tensor_functions(exponent);
            }
            _ => {}
        }
    }

    fn insert_values_for_tensor(
        &self,
        tensor: &spenso_hep_lib::HepTensor<AbstractIndex>,
        constants: &mut HashMap<Atom, Complex64>,
    ) {
        for (_, value) in tensor.iter_flat() {
            if let AtomViewOrConcrete::Atom(atom) = value {
                self.insert_function_values(atom, constants)
                    .unwrap_or_else(|err| {
                        panic!("failed to evaluate function applications in `{atom}`: {err}")
                    });
            }
        }
    }

    fn insert_function_values<'a>(
        &self,
        expression: AtomView<'a>,
        constants: &mut HashMap<Atom, Complex64>,
    ) -> Result<(), EvaluationError> {
        match expression {
            AtomView::Fun(fun) => {
                for arg in fun.iter() {
                    self.insert_function_values(arg, constants)?;
                }

                if let Some(function) = self.functions.get(&fun.get_symbol()) {
                    let args = fun
                        .iter()
                        .map(|arg| arg.evaluate(constants))
                        .collect::<Result<Vec<_>, _>>()?;
                    constants.insert(expression.to_owned(), function(&args));
                }
            }
            AtomView::Add(add) => {
                for arg in add.iter() {
                    self.insert_function_values(arg, constants)?;
                }
            }
            AtomView::Mul(mul) => {
                for arg in mul.iter() {
                    self.insert_function_values(arg, constants)?;
                }
            }
            AtomView::Pow(pow) => {
                let (base, exponent) = pow.get_base_exp();
                self.insert_function_values(base, constants)?;
                self.insert_function_values(exponent, constants)?;
            }
            _ => {}
        }

        Ok(())
    }

    fn seeded_tensor_component(symbol: Symbol, args: &[Complex64]) -> Complex64 {
        let mut state = Self::SEED;
        for byte in symbol.to_string().bytes() {
            state = Self::mix(state ^ u64::from(byte));
        }
        for arg in args {
            state = Self::mix(state ^ arg.re.to_bits());
            state = Self::mix(state ^ arg.im.to_bits().rotate_left(1));
        }

        Complex64::new(Self::sample(&mut state), Self::sample(&mut state))
    }

    fn sample(state: &mut u64) -> f64 {
        *state = Self::mix(*state);
        let unit = ((*state >> 11) as f64) * (1. / ((1_u64 << 53) as f64));
        2. * unit - 1.
    }

    fn mix(mut value: u64) -> u64 {
        value = value.wrapping_add(0x9e37_79b9_7f4a_7c15);
        value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        value ^ (value >> 31)
    }
}

struct TensorValueSummary {
    entries: usize,
    nonzero_entries: usize,
    max_abs: f64,
}

impl TensorValueSummary {
    const TOLERANCE: f64 = 1e-12;

    fn from_tensor(tensor: &spenso_hep_lib::HepTensor<AbstractIndex>) -> Self {
        let mut summary = Self {
            entries: 0,
            nonzero_entries: 0,
            max_abs: 0.,
        };

        for (_, value) in tensor.iter_flat() {
            summary.entries += 1;
            let magnitude = Self::value_magnitude(value);
            if let Some(magnitude) = magnitude {
                summary.max_abs = summary.max_abs.max(magnitude);
            }
            if !Self::is_nonzero(magnitude) {
                continue;
            }

            summary.nonzero_entries += 1;
        }

        summary
    }

    #[track_caller]
    fn assert_nontrivial(&self, label: &str, role: &str) {
        assert!(
            self.nonzero_entries > 0,
            "{role} side of `{label}` evaluated to an all-zero tensor; this would make the validation trivial"
        );
    }

    fn format(&self) -> String {
        if self.entries == 0 {
            return "empty".to_owned();
        }

        format!(
            "nonzero={}/{}, max_abs={:.6}",
            self.nonzero_entries, self.entries, self.max_abs
        )
    }

    fn value_magnitude(value: AtomViewOrConcrete<'_, RealOrComplexRef<'_, f64>>) -> Option<f64> {
        match value {
            AtomViewOrConcrete::Atom(value) => atom_magnitude(value),
            AtomViewOrConcrete::Concrete(value) => match value {
                RealOrComplexRef::Real(value) => Some(value.abs()),
                RealOrComplexRef::Complex(value) => Some(value.re.hypot(value.im)),
            },
        }
    }

    fn is_nonzero(magnitude: Option<f64>) -> bool {
        magnitude
            .map(|magnitude| magnitude > TensorValueSummary::TOLERANCE)
            .unwrap_or(true)
    }
}

fn atom_magnitude(value: AtomView<'_>) -> Option<f64> {
    let AtomView::Num(number) = value else {
        return None;
    };

    match number.get_coeff_view() {
        CoefficientView::Natural(re_num, re_den, im_num, im_den) => {
            let re = re_num as f64 / re_den as f64;
            let im = im_num as f64 / im_den as f64;
            Some(re.hypot(im))
        }
        CoefficientView::Float(re, im) => {
            let re = re.to_float().to_f64();
            let im = im.to_float().to_f64();
            Some(re.hypot(im))
        }
        _ => None,
    }
}

fn levi_civita(args: &[Complex64]) -> f64 {
    let indices = if args.len() == 1 {
        decode_cind(args[0].re.round() as usize)
    } else {
        args.iter()
            .map(|arg| arg.re.round() as usize)
            .collect::<Vec<_>>()
    };

    if indices
        .iter()
        .enumerate()
        .any(|(position, index)| indices[position + 1..].contains(index))
    {
        return 0.;
    }

    let inversions = indices
        .iter()
        .enumerate()
        .map(|(i, left)| {
            indices[i + 1..]
                .iter()
                .filter(|right| left > *right)
                .count()
        })
        .sum::<usize>();

    if inversions % 2 == 0 { 1. } else { -1. }
}

fn encode_cind(args: &[Complex64]) -> usize {
    args.iter()
        .map(|arg| arg.re.round() as usize + 1)
        .fold((0, 1), |(encoded, place), index| {
            (encoded + index * place, place * CIND_BASE)
        })
        .0
}

fn decode_cind(mut encoded: usize) -> Vec<usize> {
    let mut indices = Vec::new();
    while encoded != 0 {
        indices.push(encoded % CIND_BASE - 1);
        encoded /= CIND_BASE;
    }
    indices
}

const CIND_BASE: usize = 32;
