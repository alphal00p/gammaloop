use ahash::{HashMap, HashMapExt};
use idenso::{
    color::CS,
    epsilon::EPSILON_SYMBOL,
    reference_cases::{ReferenceCase, ReferenceDomain, reference_cases},
};
use spenso::{
    algebra::upgrading_arithmetic::FallibleSub,
    iterators::IteratableTensor,
    network::parsing::ParseSettings,
    shadowing::Concretize,
    structure::{
        TensorStructure,
        abstract_index::{AIND_SYMBOLS, AbstractIndex},
        representation::{Minkowski, RepName},
        slot::{DualSlotTo, IsAbstractSlot},
    },
    sym,
    symbolica_atom::ProjectorExpander,
    tensors::{
        data::{DenseTensor, SparseOrDense},
        parametric::atomcore::TensorAtomMaps,
    },
    trace,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    coefficient::CoefficientView,
    domains::float::RealLike,
    evaluate::EvaluationFn,
    function,
};

use crate::common::{HepAtomExt, NetExt, test_initialize};
use idenso::representations::{ColorAdjoint, ColorFundamental};

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
fn color_trace_structure_contraction_matches_commutator() {
    test_initialize();

    let cof = ColorFundamental {}.new_rep(3);
    let coad = ColorAdjoint {}.new_rep(8);
    let [a, b, c, d, x] = ["a", "b", "c", "d", "x"]
        .map(|name| coad.slot::<AbstractIndex, _>(symbolica::symbol!(name)));

    let lhs = Atom::i()
        * idenso::color_f!(c.clone(), d.clone(), x.clone())
        * trace!(
            &cof,
            idenso::color_t!(a.clone()),
            idenso::color_t!(b.clone()),
            idenso::color_t!(x)
        );
    let rhs = trace!(
        &cof,
        idenso::color_t!(a.clone()),
        idenso::color_t!(b.clone()),
        idenso::color_t!(c.clone()),
        idenso::color_t!(d.clone())
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
        * idenso::color_f!(c.clone(), d.clone(), x.clone())
        * trace!(
            &cof,
            idenso::color_t!(a.clone()),
            idenso::color_t!(x),
            idenso::color_t!(b.clone())
        );
    let rhs = trace!(
        &cof,
        idenso::color_t!(a.clone()),
        idenso::color_t!(c.clone()),
        idenso::color_t!(d.clone()),
        idenso::color_t!(b.clone())
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
        * idenso::color_f!(c.clone(), d.clone(), x.clone())
        * function!(
            CS.t,
            a.clone().to_atom(),
            i.clone().to_atom(),
            j.dual().to_atom()
        )
        * function!(CS.t, x.to_atom(), j.to_atom(), k.dual().to_atom())
        * function!(CS.t, b.clone().to_atom(), k.to_atom(), i.dual().to_atom());
    let rhs = trace!(
        &cof,
        idenso::color_t!(a.clone()),
        idenso::color_t!(c.clone()),
        idenso::color_t!(d.clone()),
        idenso::color_t!(b.clone())
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
        idenso::color_t!(a.clone()),
        idenso::color_t!(b.clone()),
        idenso::color_t!(c.clone())
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
        idenso::color_t!(a.clone()),
        idenso::color_t!(b.clone()),
        idenso::color_t!(c.clone())
    ) - trace!(
        &cof,
        idenso::color_t!(a.clone()),
        idenso::color_t!(c.clone()),
        idenso::color_t!(b.clone())
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

    let abc = function!(
        CS.t,
        a.clone().to_atom(),
        i.clone().to_atom(),
        j.dual().to_atom()
    ) * function!(
        CS.t,
        b.clone().to_atom(),
        j.clone().to_atom(),
        k.dual().to_atom()
    ) * function!(
        CS.t,
        c.clone().to_atom(),
        k.clone().to_atom(),
        i.dual().to_atom()
    );
    let acb = function!(
        CS.t,
        a.clone().to_atom(),
        i.clone().to_atom(),
        j.dual().to_atom()
    ) * function!(CS.t, c.clone().to_atom(), j.to_atom(), k.dual().to_atom())
        * function!(CS.t, b.clone().to_atom(), k.to_atom(), i.dual().to_atom());
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
        * idenso::color_f!(c.clone(), d.clone(), x.clone())
        * trace!(
            &cof,
            sym!(
                idenso::color_t!(a.clone()),
                idenso::color_t!(b.clone()),
                idenso::color_t!(x)
            )
        );
    let rhs = (trace!(
        &cof,
        idenso::color_t!(a.clone()),
        idenso::color_t!(b.clone()),
        idenso::color_t!(c.clone()),
        idenso::color_t!(d.clone())
    ) - trace!(
        &cof,
        idenso::color_t!(a.clone()),
        idenso::color_t!(b.clone()),
        idenso::color_t!(d.clone()),
        idenso::color_t!(c.clone())
    ) + trace!(
        &cof,
        idenso::color_t!(a.clone()),
        idenso::color_t!(c.clone()),
        idenso::color_t!(d.clone()),
        idenso::color_t!(b.clone())
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
        idenso::color_t!(a.clone()),
        idenso::color_t!(b.clone()),
        idenso::color_t!(c.clone()),
        idenso::color_t!(d.clone())
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
        idenso::color_t!(a.clone()),
        idenso::color_t!(c.clone()),
        idenso::color_t!(d.clone()),
        idenso::color_t!(b.clone())
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
        idenso::color_t!(a.clone()),
        idenso::color_t!(x.clone()),
        idenso::color_t!(b.clone())
    );
    let rhs = function!(CS.t, a.to_atom(), i.clone().to_atom(), j.dual().to_atom())
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
        * idenso::color_f!(a.clone(), b.clone(), c.clone())
        * function!(CS.t, c.to_atom(), i.clone().to_atom(), j.dual().to_atom());
    let rhs = function!(
        CS.t,
        a.clone().to_atom(),
        i.clone().to_atom(),
        k.dual().to_atom()
    ) * function!(
        CS.t,
        b.clone().to_atom(),
        k.clone().to_atom(),
        j.dual().to_atom()
    ) - function!(CS.t, b.to_atom(), i.to_atom(), k.dual().to_atom())
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
        idenso::color_t!(a.clone()),
        idenso::color_t!(b.clone()),
        idenso::color_t!(c.clone()),
        idenso::color_t!(d.clone())
    );
    let rhs = trace!(
        &cof,
        sym!(
            idenso::color_t!(a.clone()),
            idenso::color_t!(b.clone()),
            idenso::color_t!(c.clone()),
            idenso::color_t!(d.clone())
        )
    ) + Atom::i() / Atom::num(2)
        * idenso::color_f!(a.clone(), b.clone(), x.clone())
        * trace!(
            &cof,
            sym!(
                idenso::color_t!(c.clone()),
                idenso::color_t!(d.clone()),
                idenso::color_t!(x.clone())
            )
        )
        + Atom::i() / Atom::num(2)
            * idenso::color_f!(c.clone(), d.clone(), x.clone())
            * trace!(
                &cof,
                sym!(
                    idenso::color_t!(a.clone()),
                    idenso::color_t!(b.clone()),
                    idenso::color_t!(x.clone())
                )
            )
        - Atom::var(CS.tr) / Atom::num(6)
            * idenso::color_f!(a.clone(), c.clone(), x.clone())
            * idenso::color_f!(b.clone(), d.clone(), x.clone())
        + Atom::var(CS.tr) / Atom::num(3)
            * idenso::color_f!(a, d, x.clone())
            * idenso::color_f!(b, c, x);

    assert_valid(lhs, rhs, &scalar_constants());
}

#[test]
fn color_single_contracted_structure_constants_satisfy_jacobi() {
    test_initialize();

    let coad = ColorAdjoint {}.new_rep(8);
    let [a, b, c, d, x] = ["a", "b", "c", "d", "x"]
        .map(|name| coad.slot::<AbstractIndex, _>(symbolica::symbol!(name)));
    let expr = idenso::color_f!(a.clone(), b.clone(), x.clone())
        * idenso::color_f!(c.clone(), d.clone(), x.clone())
        - idenso::color_f!(a.clone(), c.clone(), x.clone())
            * idenso::color_f!(b.clone(), d.clone(), x.clone())
        + idenso::color_f!(a, d, x.clone()) * idenso::color_f!(b, c, x);

    assert_valid(expr, Atom::Zero, &scalar_constants());
}

#[track_caller]
fn assert_case_valid(case: &ReferenceCase) {
    let original = case.expression();
    let simplified = case.simplify(&original);
    let constants = constants_for(case.domain);

    assert_valid_with_label(original, simplified, &constants, case.name);
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
    let label = format!("`{original}` should equal `{simplified}`");
    assert_valid_with_label(original, simplified, constants, &label);
}

#[track_caller]
fn assert_valid_with_label(
    original: Atom,
    simplified: Atom,
    constants: &HashMap<Atom, Complex64>,
    label: &str,
) {
    if simplified.is_zero() {
        assert_hep_tensor_zero(evaluate(original, constants), label);
        return;
    }

    if original.is_zero() {
        assert_hep_tensor_zero(evaluate(simplified, constants), label);
        return;
    }

    let original_result = evaluate(original, constants);
    let simplified_result = evaluate(simplified, constants);
    let difference = original_result.sub_fallible(&simplified_result).unwrap();

    assert_hep_tensor_zero(difference, label);
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

    evaluate_term(expression, constants)
}

fn evaluate_term(
    expression: Atom,
    constants: &HashMap<Atom, Complex64>,
) -> spenso_hep_lib::HepTensor<AbstractIndex> {
    let net = expression
        .parse_to_hep_net(&ParseSettings::default())
        .unwrap_or_else(|err| {
            panic!("failed to parse validation expression `{expression}`: {err}")
        });

    let mut result = net.execute_and_res().unwrap_or_else(|err| {
        panic!("failed to execute validation expression `{expression}`: {err}")
    });

    let function_map = validation_functions();
    result.evaluate_complex(|c| c.into(), constants, &function_map);
    result.to_dense()
}

fn validation_functions() -> HashMap<Symbol, EvaluationFn<Atom, Complex64>> {
    let mut functions = HashMap::new();
    functions.insert(
        AIND_SYMBOLS.cind,
        EvaluationFn::new(Box::new(|args, _, _, _| {
            Complex64::new(encode_cind(args) as f64, 0.)
        })),
    );
    functions.insert(
        *EPSILON_SYMBOL,
        EvaluationFn::new(Box::new(|args, _, _, _| {
            Complex64::new(0., -levi_civita(args))
        })),
    );
    functions
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
