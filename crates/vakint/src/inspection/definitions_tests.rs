use std::collections::BTreeSet;

use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::domains::float::{Complex, Float, Real};
use symbolica::domains::rational::Rational;

use super::*;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{Vakint, VakintSettings};

fn arguments(view: AtomView) -> Vec<AtomView> {
    match view {
        AtomView::Fun(function) => function.iter().collect(),
        _ => panic!("{view} is not a function"),
    }
}

fn head(view: AtomView) -> Option<Symbol> {
    match view {
        AtomView::Fun(function) => Some(function.get_symbol()),
        _ => None,
    }
}

fn numerical(expression: &Atom, digits: u32) -> Complex<Float> {
    let settings = VakintSettings {
        run_time_decimal_precision: digits,
        ..Default::default()
    };
    crate::hpl::substitute_numerically(expression.as_view(), &settings)
        .unwrap()
        .evaluate_with_prec(&ahash::HashMap::<Atom, Complex<Float>>::default(), 300)
        .unwrap()
}

fn small(value: &Complex<Float>, bound: f64) -> bool {
    let bound = Float::with_val(300, bound);
    value.re.norm() < bound && value.im.norm() < bound
}

#[test]
fn every_builtin_constant_has_a_valid_round_tripping_definition() {
    Vakint::initialize_vakint_symbols();
    let mut symbols = super::super::constants::builtin_symbols();
    symbols.extend(
        [
            "Zeta(7)",
            "symbolica::zeta(5)",
            "Sqrt(5)",
            "symbolica::polylog(4,1/2)",
            "symbolica::polylog(2,-1/3)",
            "PolylogReal(2,(1+3^(1/2)*𝑖)/2)",
            "PolylogImag(3,3^(1/2)*𝑖/3)",
            "HR(0,0,om1)",
            "HI(1,-1,0,1,0,0,(1-3^(1/2)*𝑖)/2)",
            "GPLImag(0,1,(1+3^(1/2)*𝑖)/2,1)",
        ]
        .map(|text| vk_parse!(text).unwrap()),
    );
    assert!(symbols.len() > 40);
    for symbol in symbols {
        let definition = builtin(symbol.as_view())
            .unwrap()
            .unwrap_or_else(|| panic!("{symbol} has no definition"));
        validate(definition.as_view()).unwrap_or_else(|error| panic!("{symbol}: {error}"));
        let text = definition.to_canonical_string();
        assert_eq!(vk_parse!(&text).unwrap(), definition, "{symbol}: {text}");
        assert!(
            !definition.contains_symbol(vk_symbol!("om1")),
            "{symbol}: the HPL endpoint must be explicit"
        );
    }
    for invalid in ["Zeta(1)", "Sqrt(-1)", "HR(2,om1)", "HI(1,0,1/2)"] {
        assert!(
            builtin(vk_parse!(invalid).unwrap().as_view()).is_err(),
            "{invalid}"
        );
    }
    for foreign in [
        "user_space::HR(1,0,vakint::om1)",
        "user_space::D6",
        "PR4fin",
    ] {
        assert!(
            builtin(vk_parse!(foreign).unwrap().as_view())
                .unwrap()
                .is_none(),
            "{foreign}"
        );
    }
}

#[test]
fn validation_rejects_free_misplaced_and_unknown_symbols() {
    Vakint::initialize_vakint_symbols();
    for valid in [
        "vakint::Integral(vakint::x(1)*vakint::x(2),vakint::Interval(vakint::x(1),0,1),vakint::Interval(vakint::x(2),0,vakint::x(1)))",
        "vakint::SeriesCoefficient(vakint::LoopIntegral(2-2*vakint::delta,(vakint::EuclideanDot(vakint::k(1),vakint::k(1))+1)^-2,vakint::k(1)),vakint::delta,0,1)",
        "symbolica::real_part(vakint::Integral(symbolica::exp(-vakint::t)*𝜋*symbolica::euler_gamma,vakint::Interval(vakint::t,0,∞)))",
    ] {
        validate(vk_parse!(valid).unwrap().as_view())
            .unwrap_or_else(|error| panic!("{valid}: {error}"));
    }
    for (invalid, reason) in [
        (
            "vakint::Integral(vakint::t*vakint::y,vakint::Interval(vakint::t,0,1))",
            "free symbol",
        ),
        (
            "vakint::Integral(vakint::x(2),vakint::Interval(vakint::x(1),0,1))",
            "free variable",
        ),
        (
            "vakint::Integral(vakint::t,vakint::Interval(vakint::t,0,vakint::t))",
            "free symbol",
        ),
        (
            "vakint::Integral(vakint::t,vakint::Interval(vakint::t,0,1),vakint::Interval(vakint::t,0,1))",
            "bound twice",
        ),
        (
            "vakint::Integral(vakint::foo(vakint::t),vakint::Interval(vakint::t,0,1))",
            "not admitted",
        ),
        (
            "vakint::Integral(vakint::K5Moment(1)*vakint::t,vakint::Interval(vakint::t,0,1))",
            "not admitted",
        ),
        ("vakint::Integral(vakint::t)", "without an Interval"),
        (
            "vakint::Integral(1.5`20*vakint::t,vakint::Interval(vakint::t,0,1))",
            "exact",
        ),
        (
            "vakint::LoopIntegral(4-2*vakint::ep,(vakint::EuclideanDot(vakint::k(1),vakint::k(1))+1)^-1,vakint::k(1))",
            "free symbol",
        ),
        (
            "vakint::LoopIntegral(2,(vakint::EuclideanDot(vakint::k(1),vakint::k(2))+1)^-1,vakint::k(1))",
            "bound loop momentum",
        ),
        (
            "vakint::LoopIntegral(2,vakint::k(1),vakint::k(1))",
            "not admitted",
        ),
        ("vakint::Interval(vakint::t,0,1)", "not admitted"),
        (
            "vakint::SeriesCoefficient(vakint::ep,vakint::ep,0,1/2)",
            "integer",
        ),
    ] {
        let error = validate(vk_parse!(invalid).unwrap().as_view()).unwrap_err();
        assert!(error.contains(reason), "{invalid}: {error}");
    }
}

/// Count the nested intervals and find the innermost logarithm power.
fn hpl_structure(definition: &Atom) -> (Atom, usize, i64) {
    let [kernel] = arguments(definition.as_view())[..] else {
        panic!("projection of one kernel")
    };
    if head(kernel) != Some(vk_symbol!("Integral")) {
        return (kernel.to_owned(), 0, 0);
    }
    let parts = arguments(kernel);
    let mut log_power = 0;
    parts[0].to_owned().visitor(&mut |part| {
        if let AtomView::Pow(power) = part {
            let (base, exponent) = power.get_base_exp();
            if head(base) == Some(Symbol::LOG) {
                log_power = get_integer_from_atom(exponent).unwrap();
            }
        }
        if head(part) == Some(Symbol::LOG) && log_power == 0 {
            log_power = 1;
        }
        true
    });
    (parts[0].to_owned(), parts.len() - 1, log_power)
}

#[test]
fn hpl_definitions_follow_the_radial_path_convention() {
    Vakint::initialize_vakint_symbols();
    // The endpoint (1±𝑖√3)/2 keeps its numerator sum inside products.
    let positive = vk_parse!("1+𝑖*3^(1/2)").unwrap();
    let negative = vk_parse!("1-𝑖*3^(1/2)").unwrap();
    for (word, sign, depth, logarithms) in [
        (vec![1, 0], 1, 1, 1),
        (vec![1, 0, 0, 0], -1, 1, 3),
        (vec![1, -1, 0, 1], 1, 4, 0),
        (vec![-1, 1, 1, 0, 0], 1, 3, 2),
        (vec![0, 0, 1], -1, 3, 0),
    ] {
        let definition = hpl(&word, sign, false);
        assert_eq!(head(definition.as_view()), Some(real_part()));
        let (integrand, intervals, power) = hpl_structure(&definition);
        assert_eq!(intervals, depth, "{word:?}");
        assert_eq!(power, logarithms, "{word:?}");
        let endpoint = if sign == 1 { &positive } else { &negative };
        let conjugate = if sign == 1 { &negative } else { &positive };
        assert!(
            integrand.contains(endpoint.as_view())
                || word[..depth].iter().all(|letter| *letter == 0),
            "{word:?}: {integrand}"
        );
        assert!(!integrand.contains(conjugate.as_view()));
    }
    let all_zero = hpl(&[0, 0], 1, true);
    assert_eq!(head(all_zero.as_view()), Some(imaginary_part()));
    assert!(!all_zero.contains_symbol(vk_symbol!("Integral")));
}

/// Gamma-function value of a loop integral that factorizes, after an
/// unimodular change of loop momenta, into one-loop massive tadpoles:
/// `∫ d^Dk/π^(D/2) (k²+1)^(-a) = Γ(a-D/2)/Γ(a)`. `Γ(n+ep)` is written as
/// `Γ(1+ep)` times rational factors so the Laurent series stays regular.
fn tadpole_product(loop_integral: AtomView) -> Atom {
    let parts = arguments(loop_integral);
    let (dimension, integrand, momenta) = (parts[0], parts[1], &parts[2..]);
    let loops = momenta.len();
    let index = |momentum: AtomView| {
        momenta
            .iter()
            .position(|candidate| *candidate == momentum)
            .unwrap()
    };
    let factors = match integrand {
        AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
        _ => vec![integrand],
    };
    let mut vectors = Vec::new();
    let mut value = Atom::num(1);
    let epsilon = Atom::var(S.form_epsilon);
    assert_eq!(
        (dimension.to_owned() - Atom::num(4) + Atom::num(2) * &epsilon).expand(),
        Atom::Zero
    );
    for factor in factors {
        let AtomView::Pow(power) = factor else {
            panic!("{factor} is not a propagator power")
        };
        let (base, exponent) = power.get_base_exp();
        let power = -get_integer_from_atom(exponent).unwrap();
        let mut quadratic = vec![vec![Rational::zero(); loops]; loops];
        let mut constant = Rational::zero();
        let AtomView::Add(terms) = base else {
            panic!("{base} is not a massive propagator")
        };
        for term in terms.iter() {
            let (coefficient, dot) = match term {
                AtomView::Num(_) => {
                    constant += &Rational::try_from(term).unwrap();
                    continue;
                }
                AtomView::Mul(product) => {
                    let parts = product.iter().collect::<Vec<_>>();
                    (Rational::try_from(parts[0]).unwrap(), parts[1])
                }
                _ => (Rational::one(), term),
            };
            assert_eq!(head(dot), Some(S.euclidean_dot));
            let pair = arguments(dot);
            let (i, j) = (index(pair[0]), index(pair[1]));
            if i == j {
                quadratic[i][i] += &coefficient;
            } else {
                let half = &coefficient / &Rational::from(2);
                quadratic[i][j] += &half;
                quadratic[j][i] += &half;
            }
        }
        assert_eq!(constant, Rational::one());
        // A rank-one form v vᵀ with an integer vector v.
        let pivot = (0..loops).find(|i| !quadratic[*i][*i].is_zero()).unwrap();
        let root = (1..=4)
            .map(Rational::from)
            .find(|root| (root * root) == quadratic[pivot][pivot])
            .unwrap();
        let vector = (0..loops)
            .map(|j| &quadratic[pivot][j] / &root)
            .collect::<Vec<_>>();
        for i in 0..loops {
            for j in 0..loops {
                assert_eq!(quadratic[i][j], &vector[i] * &vector[j]);
            }
        }
        vectors.push(vector);
        // Γ(n+ep) = Γ(1+ep)/∏_{m=n}^{0}(m+ep) for n < 1, Γ(1+ep)∏_{m=1}^{n-1}(m+ep) otherwise.
        let mut gamma_value = symbolica::function!(gamma(), Atom::num(1) + &epsilon);
        let mut shift = power - 2;
        while shift < 1 {
            gamma_value /= Atom::num(shift) + &epsilon;
            shift += 1;
        }
        while shift > 1 {
            shift -= 1;
            gamma_value *= Atom::num(shift) + &epsilon;
        }
        value = value * gamma_value / factorial(power as u32 - 1);
    }
    assert_eq!(vectors.len(), loops, "one tadpole per loop");
    // Unimodular change of variables: |det V| = 1 by fraction-free elimination.
    let mut matrix = vectors;
    let mut determinant = Rational::one();
    for column in 0..loops {
        let row = (column..loops)
            .find(|row| !matrix[*row][column].is_zero())
            .unwrap();
        matrix.swap(row, column);
        let pivot = matrix[column][column].clone();
        determinant = &determinant * &pivot;
        let pivot_row = matrix[column].clone();
        for row in matrix.iter_mut().skip(column + 1) {
            let factor = &row[column] / &pivot;
            for (entry, above) in row.iter_mut().zip(&pivot_row).skip(column) {
                *entry = &*entry - &(&factor * above);
            }
        }
    }
    assert!(determinant == 1 || determinant == -1);
    value
}

/// Evaluate `SeriesCoefficient(prefactor * LoopIntegral(…), ep, 0, n)` for tadpole products.
fn evaluate_tadpole_coefficient(definition: &Atom) -> Atom {
    let parts = arguments(definition.as_view());
    assert_eq!(
        head(definition.as_view()),
        Some(vk_symbol!("SeriesCoefficient"))
    );
    assert_eq!(parts[1], Atom::var(S.form_epsilon).as_view());
    assert_eq!(parts[2], Atom::Zero.as_view());
    let order = get_integer_from_atom(parts[3]).unwrap();
    let mut loop_integral = None;
    parts[0].to_owned().visitor(&mut |part| {
        if head(part) == Some(vk_symbol!("LoopIntegral")) {
            loop_integral = Some(part.to_owned());
            return false;
        }
        true
    });
    let loop_integral = loop_integral.unwrap();
    let value = tadpole_product(loop_integral.as_view());
    let function = parts[0]
        .to_owned()
        .replace(loop_integral.to_pattern())
        .with(value.to_pattern());
    let series = function
        .series(
            S.form_epsilon,
            Atom::Zero.as_view(),
            Rational::from(order + 1),
        )
        .unwrap();
    series
        .coefficient(Rational::from(order))
        .unwrap_or(Atom::Zero)
}

#[test]
fn tadpole_products_fix_the_public_master_normalization_exactly() {
    Vakint::initialize_vakint_symbols();
    // The one-loop FMFT/MATAD tadpole is -Gam(1,1)/(ep (ep-1)) = -e^(γ ep) Γ(-1+ep).
    let one_loop = EuclideanFamily::from_momenta(1, &[(&[1], true)]);
    let reference = -symbolica::function!(gamma(), Atom::num(1) + Atom::var(S.form_epsilon))
        * exp(Atom::var(euler_gamma()) * Atom::var(S.form_epsilon))
        / (Atom::var(S.form_epsilon) * (Atom::var(S.form_epsilon) - Atom::num(1)));
    let reference = reference
        .series(S.form_epsilon, Atom::Zero.as_view(), Rational::from(4))
        .unwrap();
    for order in -1..=3 {
        let definition = one_loop.public_master(&[1], order).unwrap();
        let residual = evaluate_tadpole_coefficient(&definition)
            - reference
                .coefficient(Rational::from(order))
                .unwrap_or(Atom::Zero);
        let residual = numerical(&residual, 60);
        assert!(small(&residual, 1e-55), "order {order}");
    }

    // PR1 is the product of four tadpoles; compare the shipped FMFT table.
    let expansion = &crate::fmft_numerics::MASTERS_EXPANSION[&vk_parse!("PR1_exact").unwrap()].0;
    let exact = expansion
        .series(S.form_epsilon, Atom::Zero.as_view(), Rational::from(1))
        .unwrap();
    for order in -4..=0 {
        let definition = fmft_master("PR1", order).unwrap();
        let expected = crate::hpl::normalize_constants(
            exact
                .coefficient(Rational::from(order))
                .unwrap_or(Atom::Zero)
                .as_view(),
        )
        .unwrap();
        let residual = numerical(&(evaluate_tadpole_coefficient(&definition) - expected), 60);
        assert!(small(&residual, 1e-55), "PR1 order {order}");
    }
    for (order, name) in [(1, "PR1ep1"), (2, "PR1ep2")] {
        let table =
            &crate::fmft_numerics::MASTERS_NUMERIC_SUBSTITUTIONS[&vk_parse!(name).unwrap()].0;
        let residual = numerical(
            &(evaluate_tadpole_coefficient(&fmft_master("PR1", order).unwrap()) - table.clone()),
            60,
        );
        assert!(small(&residual, 1e-44), "{name}");
    }
    // An odd total power flips the sign of the Euclidean integral.
    let definition = fmft_master("PR4", 0).unwrap();
    let sign = arguments(definition.as_view())[0];
    assert!(matches!(sign, AtomView::Mul(product)
        if product.iter().any(|factor| factor == Atom::num(-1).as_view())));
}

#[test]
fn legacy_tables_name_complete_master_coefficients() {
    Vakint::initialize_vakint_symbols();
    let epsilon = Atom::var(S.form_epsilon);
    let expansions = &crate::fmft_numerics::MASTERS_EXPANSION;
    let mut covered = BTreeSet::new();
    for (name, master, order) in FMFT_TABLE {
        let full = &expansions[&vk_parse!(master).unwrap()].0;
        let coefficient = full
            .coefficient_list::<i32>(std::slice::from_ref(&epsilon))
            .into_iter()
            .find(|(monomial, _)| *monomial == epsilon.pow(Atom::num(order)))
            .unwrap()
            .1;
        let exact = vk_parse!(&format!("{master}_exact")).unwrap();
        let remainder = coefficient.clone() - Atom::var(vk_symbol!(name));
        assert!(
            remainder.is_zero() || remainder == exact,
            "{name}: {coefficient}"
        );
        if let Some((exact_expansion, _)) = expansions.get(&exact) {
            assert!(
                exact_expansion
                    .coefficient_list::<i32>(std::slice::from_ref(&epsilon))
                    .iter()
                    .all(|(monomial, _)| *monomial != epsilon.pow(Atom::num(order))),
                "{name}: the exact part overlaps"
            );
        }
        assert!(fmft_master(master, order).is_some());
        covered.insert(name);
    }
    for name in crate::fmft_numerics::MASTERS_NUMERIC_SUBSTITUTIONS.keys() {
        let name = name.to_canonical_string();
        let name = name.rsplit("::").next().unwrap();
        assert!(
            covered.contains(name) || ["PR3ep2", "PR13ep1", "D6", "T1ep", "T1ep2"].contains(&name),
            "{name} is neither defined nor documented as undefined"
        );
    }
    for name in crate::matad_numerics::DIRECT_SUBSTITUTIONS.keys() {
        let name = name.to_canonical_string();
        let name = name.rsplit("::").next().unwrap();
        assert!(matad_master(name, 0).is_some(), "{name}");
    }
}

fn staged_leaf(backend: &str, original: &str, table: &str) -> Atom {
    let staged = crate::precision::legacy::stage_table(
        vk_parse!(original).unwrap().as_view(),
        vk_parse!(table).unwrap().as_view(),
        S.form_epsilon,
        backend,
    )
    .unwrap();
    crate::precision::legacy::collect_sources(staged.as_view())
        .unwrap()
        .pop_first()
        .unwrap()
        .0
}

#[test]
fn legacy_leaves_are_defined_or_rejected_with_a_typed_error() {
    Vakint::initialize_vakint_symbols();
    let leaf = staged_leaf("FMFT", "PR4fin", "2.5`50");
    assert_eq!(
        legacy_source(leaf.as_view()).unwrap().unwrap(),
        fmft_master("PR4", 0).unwrap()
    );
    let leaf = staged_leaf("MATAD", "miD6", "2.25`40/ep");
    assert_eq!(
        legacy_source(leaf.as_view()).unwrap().unwrap(),
        matad_master("miD6", -1).unwrap()
    );
    assert_eq!(
        matad_master("miD6", 0).unwrap(),
        builtin(vk_parse!("D6").unwrap().as_view())
            .unwrap()
            .unwrap()
    );
    for (backend, original, table) in [
        ("FMFT", "PR3ep2", "717.31`50"),
        ("FMFT", "PR13ep1", "27.0`50"),
        ("MATAD", "miD6", "2.25`40*ep^-3"),
        ("MATAD", "unknown_master", "2.25`40"),
        ("Input", "input_coefficient", "0.5`20"),
        ("PySecDec", "pysecdec_result", "0.5`20"),
        ("AlphaLoop", "MASTER2m22EP", "0.67`30"),
    ] {
        let leaf = staged_leaf(backend, original, table);
        assert!(
            matches!(
                legacy_source(leaf.as_view()),
                Err(VakintError::Inspection(
                    InspectionError::UndefinedConstant { .. }
                ))
            ),
            "{backend} {original}"
        );
    }
    assert!(
        legacy_source(vk_parse!("PR4fin").unwrap().as_view())
            .unwrap()
            .is_none()
    );
}

#[test]
fn external_fallbacks_follow_the_recorded_metadata() {
    Vakint::initialize_vakint_symbols();
    let null = serde_json::Value::Null;
    let moment = external_fallback(&vk_parse!("K5Moment(3)").unwrap(), &null, &[])
        .unwrap()
        .unwrap();
    assert_eq!(
        moment,
        vk_parse!("vakint::Integral(vakint::t^3*symbolica::bessel_k(0,vakint::t)^5,vakint::Interval(vakint::t,0,∞))").unwrap()
    );
    let logarithmic = external_fallback(&vk_parse!("K5LogMoment(1)").unwrap(), &null, &[])
        .unwrap()
        .unwrap();
    assert_eq!(
        logarithmic,
        vk_parse!("3*vakint::Integral(vakint::t*symbolica::bessel_k(0,vakint::t)^5*(symbolica::log(vakint::t)+symbolica::euler_gamma-symbolica::log(2)),vakint::Interval(vakint::t,0,∞))").unwrap()
    );
    let wrong = serde_json::json!({"definition_metadata": {"power": 3, "logarithmic": false}});
    assert!(external_fallback(&vk_parse!("K5Moment(1)").unwrap(), &wrong, &[]).is_err());
    let wrong = serde_json::json!({"definition_metadata": {"power": 1, "logarithmic": false}});
    assert!(external_fallback(&vk_parse!("K5LogMoment(1)").unwrap(), &wrong, &[]).is_err());
    assert!(
        external_fallback(&vk_parse!("unrelated_period").unwrap(), &null, &[])
            .unwrap()
            .is_none()
    );

    // A two-loop sunset family: rows k1², k2², (k1+k2)² in the external-data schema.
    let row = |coefficients: [&str; 3]| serde_json::json!({"coefficients": coefficients, "constant": "-1"});
    let family = serde_json::json!({
        "sha256": "test-family",
        "coordinate_index_base": "0",
        "denominator_convention": "constant + sum_j coefficients[j] * coordinates[j]",
        "external_momenta": [],
        "coordinates": [
            {"kind": "loop_loop", "left": 0, "right": 0},
            {"kind": "loop_loop", "left": 0, "right": 1},
            {"kind": "loop_loop", "left": 1, "right": 1},
        ],
        "denominators": [row(["1", "0", "0"]), row(["0", "0", "1"]), row(["1", "2", "1"])],
    });
    let provenance = serde_json::json!({"definition_metadata": {
        "dimension": 2, "family_sha256": "test-family", "powers": [1, 1, 2]
    }});
    let period = external_fallback(
        &vk_parse!("PR7D2").unwrap(),
        &provenance,
        std::slice::from_ref(&family),
    )
    .unwrap()
    .unwrap();
    validate(period.as_view()).unwrap();
    let expected = vk_parse!(
        "vakint::LoopIntegral(2,(vakint::EuclideanDot(vakint::k(1),vakint::k(1))+1)^-1*(vakint::EuclideanDot(vakint::k(2),vakint::k(2))+1)^-1*(vakint::EuclideanDot(vakint::k(1)+vakint::k(2),vakint::k(1)+vakint::k(2))+1)^-2,vakint::k(1),vakint::k(2))"
    )
    .unwrap();
    assert_eq!(period, expected);
    let jet = serde_json::json!({"definition_metadata": {
        "Taylor_order": 1,
        "Taylor_coefficient_not_factorial_scaled_derivative": true,
        "portable_integral_definition": {
            "dimension": "2-2*delta", "ordered_family": family.clone(), "signed_index_tuple": [1, 1, 2]
        }
    }});
    let jet = external_fallback(
        &vk_parse!("PR10D2Eps1").unwrap(),
        &jet,
        std::slice::from_ref(&family),
    )
    .unwrap()
    .unwrap();
    validate(jet.as_view()).unwrap();
    assert_eq!(head(jet.as_view()), Some(vk_symbol!("SeriesCoefficient")));
    let mismatched = serde_json::json!({"definition_metadata": {"Taylor_order": 2,
        "portable_integral_definition": {"ordered_family": family.clone(), "signed_index_tuple": [1, 1, 2]}}});
    assert!(external_fallback(&vk_parse!("PR10D2Eps1").unwrap(), &mismatched, &[]).is_err());
    let mut contradicting = family.clone();
    contradicting["denominators"][2] = row(["1", "-2", "1"]);
    let embedded = serde_json::json!({"definition_metadata": {
        "dimension": 2, "ordered_family": family.clone(), "signed_index_tuple": [1, 1, 1]
    }});
    assert!(external_fallback(&vk_parse!("PR15D2").unwrap(), &embedded, &[contradicting]).is_err());
    assert!(agree(&period, &expected));
    assert!(!agree(&period, &jet));

    // The same integral in another loop-momentum routing (k2 -> -k1-k2) agrees;
    // a rescaled momentum or a massless line does not.
    let relabeled = vk_parse!(
        "vakint::LoopIntegral(2,(vakint::EuclideanDot(vakint::k(1),vakint::k(1))+1)^-1*(vakint::EuclideanDot(vakint::k(1)+vakint::k(2),vakint::k(1)+vakint::k(2))+1)^-1*(vakint::EuclideanDot(vakint::k(2),vakint::k(2))+1)^-2,vakint::k(1),vakint::k(2))"
    )
    .unwrap();
    assert!(agree(&period, &relabeled));
    for different in [
        "vakint::LoopIntegral(2,(vakint::EuclideanDot(2*vakint::k(1),2*vakint::k(1))+1)^-1*(vakint::EuclideanDot(vakint::k(2),vakint::k(2))+1)^-1*(vakint::EuclideanDot(vakint::k(1)+vakint::k(2),vakint::k(1)+vakint::k(2))+1)^-2,vakint::k(1),vakint::k(2))",
        "vakint::LoopIntegral(2,vakint::EuclideanDot(vakint::k(1),vakint::k(1))^-1*(vakint::EuclideanDot(vakint::k(2),vakint::k(2))+1)^-1*(vakint::EuclideanDot(vakint::k(1)+vakint::k(2),vakint::k(1)+vakint::k(2))+1)^-2,vakint::k(1),vakint::k(2))",
        "vakint::LoopIntegral(4,(vakint::EuclideanDot(vakint::k(1),vakint::k(1))+1)^-1*(vakint::EuclideanDot(vakint::k(1)+vakint::k(2),vakint::k(1)+vakint::k(2))+1)^-1*(vakint::EuclideanDot(vakint::k(2),vakint::k(2))+1)^-2,vakint::k(1),vakint::k(2))",
    ] {
        assert!(
            !agree(&period, &vk_parse!(different).unwrap()),
            "{different}"
        );
    }
    // PR10 is recorded on two public families with the same powers.
    let family = |rows: &[[i64; 4]; 10]| {
        EuclideanFamily::from_momenta(
            4,
            &rows
                .iter()
                .map(|momentum| (momentum.as_slice(), true))
                .collect::<Vec<_>>(),
        )
    };
    let powers = [0, 1, 0, 1, 1, 1, 1, 1, 1, 0];
    let on_h = family(&FAMILY_4L_H).d2_period(&powers).unwrap();
    let on_x = family(&FAMILY_4L_X).d2_period(&powers).unwrap();
    assert_ne!(on_h, on_x);
    assert!(agree(&on_h, &on_x));
    let pr9 = family(&FAMILY_4L_H)
        .d2_period(&[0, 1, 1, 0, 1, 1, 1, 1, 1, 0])
        .unwrap();
    assert!(!agree(&on_h, &pr9));
}

#[test]
fn hpl_and_elementary_definitions_have_the_registered_values() {
    // Cheap exact anchors; full numerical quadrature is done offline.
    Vakint::initialize_vakint_symbols();
    let definition = builtin(vk_parse!("S2").unwrap().as_view())
        .unwrap()
        .unwrap();
    assert!(definition.contains_symbol(vk_symbol!("Integral")));
    let sqrt = builtin(vk_parse!("sqrt3").unwrap().as_view())
        .unwrap()
        .unwrap();
    let value = numerical(&sqrt, 30);
    assert!((value.re - Float::with_val(300, 3).sqrt()).norm() < Float::with_val(300, 1e-28));
}

/// Input for the offline quadrature check (not CI):
/// `VAKINT_DEFINITION_DUMP=out.json [VAKINT_DEFINITION_SNAPSHOT=dir] cargo test … -- --ignored`.
/// The snapshot directory holds `4L-constants.json` and `4L-manifest.json`.
#[test]
#[ignore = "writes offline validation input"]
fn dump_definitions_for_offline_validation() {
    let Ok(path) = std::env::var("VAKINT_DEFINITION_DUMP") else {
        return;
    };
    let engine = Vakint::new().unwrap();
    let mut symbols = super::super::constants::builtin_symbols()
        .into_iter()
        .collect::<BTreeSet<_>>();
    let mut words: Vec<Vec<i8>> = vec![vec![]];
    for _ in 0..4 {
        words = words
            .iter()
            .flat_map(|word| {
                [-1, 0, 1].map(|letter| {
                    let mut longer = word.clone();
                    longer.push(letter);
                    longer
                })
            })
            .collect::<Vec<_>>();
        for word in &words {
            for head in ["HR", "HI"] {
                let letters = word.iter().map(i8::to_string).collect::<Vec<_>>().join(",");
                symbols.insert(vk_parse!(&format!("{head}({letters},om1)")).unwrap());
            }
        }
    }
    for text in [
        "HR(1,0,1,0,1,om1)",
        "HI(-1,0,0,1,1,om1)",
        "HR(1,1,-1,0,0,om1)",
        "HI(0,1,0,-1,1,om1)",
        "HR(1,0,0,0,0,0,om1)",
        "HI(1,0,1,0,1,0,om1)",
        "HR(-1,1,0,0,1,-1,om1)",
        "HI(1,1,1,0,-1,0,(1-3^(1/2)*𝑖)/2)",
        "Zeta(2)",
        "Zeta(5)",
        "symbolica::polylog(2,1/2)",
        "symbolica::polylog(3,-1/2)",
        "symbolica::polylog(4,1/3)",
        "PolylogReal(2,(1+3^(1/2)*𝑖)/2)",
        "PolylogImag(4,(1+3^(1/2)*𝑖)/2)",
        "PolylogReal(3,3^(1/2)*𝑖/3)",
        "PolylogImag(3,3^(1/2)*𝑖/3)",
        "GPLReal(0,1,(1+3^(1/2)*𝑖)/2,1)",
        "GPLImag(-1,0,(1-3^(1/2)*𝑖)/2,1)",
    ] {
        symbols.insert(vk_parse!(text).unwrap());
    }
    let mut records = Vec::new();
    for symbol in symbols {
        let constant = engine
            .get_constant(
                &crate::inspection::ConstantSymbol::try_from(symbol.clone()).unwrap(),
                45,
            )
            .unwrap();
        records.push(serde_json::json!({
            "symbol": symbol.to_canonical_string(),
            "definition": constant.definition().to_canonical_string(),
            "value": constant.numerical_value().to_string(),
            "source": "get_constant(prec=45)",
        }));
    }
    if let Ok(directory) = std::env::var("VAKINT_DEFINITION_SNAPSHOT") {
        let read = |name: &str| -> serde_json::Value {
            serde_json::from_slice(&std::fs::read(format!("{directory}/{name}")).unwrap()).unwrap()
        };
        let constants = read("4L-constants.json");
        let manifest = read("4L-manifest.json");
        let families = manifest["families"].as_array().unwrap();
        for record in constants.as_array().unwrap() {
            let symbol = vk_parse!(record["symbol"].as_str().unwrap()).unwrap();
            let definition = external_fallback(&symbol, &record["provenance"], families)
                .unwrap()
                .unwrap();
            records.push(serde_json::json!({
                "symbol": symbol.to_canonical_string(),
                "definition": definition.to_canonical_string(),
                "value": record["value"],
                "significant_digits": record["significant_digits"],
                "source": "external-data snapshot ca49f328 4L/constants.json (fallback definition)",
                "metadata": record["provenance"]["definition_metadata"],
            }));
        }
    }
    std::fs::write(path, serde_json::to_vec_pretty(&records).unwrap()).unwrap();
}
