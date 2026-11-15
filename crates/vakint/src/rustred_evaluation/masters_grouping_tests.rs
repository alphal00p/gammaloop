//! Exact normal-form gates for denominator collection, including map/field cases
//! in which algebraic equality alone would miss a changed canonical result.
use super::{ExactDenominatorGroups, MasterCatalog, combine_exact, proxy_tests};
use crate::utils::vakint_macros::vk_parse;
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::coefficient::CoefficientView;

fn check_grouping(input: Atom) -> bool {
    let expanded = input.expand();
    let reference = proxy_tests::reference_combine_exact(input.clone());
    assert_eq!(combine_exact(input), reference);
    if MasterCatalog::has_floating_coefficients(&expanded) {
        assert_eq!(reference, expanded);
        return false;
    }
    let Ok(grouped) = ExactDenominatorGroups::new(&expanded) else {
        return false;
    };
    assert_eq!(
        *grouped.variables,
        ExactDenominatorGroups::raw_variables(expanded.as_view(), grouped.complex).unwrap(),
        "original raw variable order changed: {expanded}"
    );
    assert_eq!(grouped.combine().unwrap(), reference, "{expanded}");
    true
}

#[test]
fn denominator_groups_preserve_variable_order_coefficient_field_and_root_basis() {
    crate::Vakint::initialize_vakint_symbols();
    let cases = [
        // Regrouping without the original map changes even this univariate case.
        "(1-ep)^(-3)*x+(1-ep)^(-2)*y+(1-ep)^(-3)*ep*z",
        "x/(1-ep)+2*y/(2-2*ep)+z/(3-3*ep)",
        "𝑖*x/(1-ep)+2*𝑖*x/(2*ep-2)+1/(ep^2+1)",
        "x/(1-ep)+2*x/(2*ep-2)+y/(1+ep)",
        "x/(1-ep)+2*x/(2*ep-2)",
        "x/(1-ep)+2*x/(2*ep-2)+y/(1-ep)+2*y/(2*ep-2)",
        "x^(-ep)/(1-ep)+y*x^(-ep)/(1-ep)",
        "x^(a+ep)/(1-ep)+y*x^(a+ep)/(1-ep)",
        "x^(1/2)/(1-ep)+x^(1/3)/(1-ep)+x^(5/6)/(2-2*ep)",
        "𝑖*x^(1/2)/(1-ep)+𝑖*x^(1/3)/(1-ep)+𝑖*x^(5/6)/(2-2*ep)",
        "x^(-1/2)/(1-ep)+x^(-1/3)/(1-ep)",
        "x^(ep/2)/(1-ep)+x^(ep/3)/(1-ep)+x^(5*ep/6)/(2-2*ep)",
        "x^(a+ep/2)/(1-ep)+x^(a+ep/3)/(1-ep)+x^(a+5*ep/6)/(2-2*ep)",
        "y/(x^(1/2)+x^(1/3))+z/(x^(1/2)+x^(1/3))",
        "y/(x^(ep/2)+x^(ep/3))+z/(x^(ep/2)+x^(ep/3))",
        "𝑖*y/(x^(ep/2)+x^(ep/3))+𝑖*z/(x^(ep/2)+x^(ep/3))",
        "(2+3*𝑖)*(x/(1-ep)+y/(2-2*ep))",
        "(1+𝑖)*((1-ep)^(-3)*x+(1-ep)^(-2)*y+(1-ep)^(-3)*ep*z)",
        "𝑖*x/(ep^2+1)+𝑖*y/(2*ep^2+2)",
        "𝑖*x/(ep^2+a^2)+𝑖*y/(2*ep^2+2*a^2)",
        "(1+𝑖)*(x^(1/2)/(ep^2+1)+x^(1/3)/(2*ep^2+2))",
        "𝑖*((-x)^(1/2)/(1-ep)+(-x)^(1/3)/(2-2*ep)+(-x)^(5/6)/(1-ep))",
        "𝑖*((-2)^(1/2)*x/(1-ep)+(-2)^(1/2)*y/(2-2*ep))",
        "𝑖*f(1+𝑖)*x/(1-ep)+𝑖*y/(1-ep)",
        "𝑖*m^((1+𝑖)*ep)*x/(1-ep)+𝑖*y/(1-ep)",
        "𝑖*(𝑖^(1/2)*x/(1-ep)+𝑖^(1/2)*y/(2-2*ep))",
        "𝑖*((1+𝑖)^ep*x/(1-ep)+(1+𝑖)^ep*y/(2-2*ep))",
        "RustRedMaster(f,1)*g(a,b)/(1-ep)+RustRedMaster(f,2)*g(b,a)/(1-ep)",
        "(x+y)/(a+b)+(x-y)/(a-b)+x/(a+b)",
        "(𝑖*x+y)/(a+𝑖*b)+(x-𝑖*y)/(a-𝑖*b)+z/(a+𝑖*b)",
        "0.125*x/(1-ep)+0.25*y/(1-ep)",
        "0",
    ];
    let mut admitted = 0;
    for (index, input) in cases.into_iter().enumerate() {
        // Disjoint namespaces allow different first-registration orders in the
        // same process without depending on concurrent test registration.
        for (order_index, order) in [
            ["x", "y", "z", "ep", "a", "b"],
            ["ep", "z", "y", "x", "b", "a"],
            ["a", "b", "ep", "x", "y", "z"],
        ]
        .into_iter()
        .enumerate()
        {
            let namespace: &'static str =
                Box::leak(format!("grouping_case_{index}_{order_index}").into_boxed_str());
            for symbol in order {
                let _ = Atom::parse(symbol, namespace, Default::default()).unwrap();
            }
            let input = Atom::parse(input, namespace, Default::default()).unwrap();
            admitted += usize::from(check_grouping(input));
        }
    }
    assert!(admitted >= 40, "only {admitted} cases exercised grouping");
    for (input, admitted) in [
        ("𝑖*x/(ep^2+1)+𝑖*y/(2*ep^2+2)", true),
        ("(2+3*𝑖)*(x/(1-ep)+y/(2-2*ep))", true),
        ("𝑖*x/(1-ep)+y/(1-ep)", false),
        ("𝑖*f(1+𝑖)*x/(1-ep)+𝑖*y/(1-ep)", false),
        ("𝑖*(𝑖^(1/2)*x/(1-ep)+𝑖^(1/2)*y/(2-2*ep))", false),
        ("0.5*𝑖*x/(1-ep)+𝑖*y/(1-ep)", false),
    ] {
        let expanded = vk_parse!(input).unwrap().expand();
        let phase = ExactDenominatorGroups::common_phase(&expanded);
        assert_eq!(!phase.is_one(), admitted, "{input}");
        if admitted {
            assert!(
                !(&expanded / phase)
                    .expand()
                    .as_view()
                    .has_complex_coefficients()
            );
        }
    }
}

#[test]
fn denominator_groups_match_deterministic_mixed_rational_sums() {
    crate::Vakint::initialize_vakint_symbols();
    let denominators = ["1-ep", "1+ep", "2-2*ep", "ep^2+1", "a+b", "a-b"];
    let numerators = ["x", "y", "z*ep", "x*y", "g(a,b)", "RustRedMaster(f,1)"];
    let mut state = 20_261_107_u64;
    let mut next = || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        (state >> 32) as usize
    };
    let mut admitted = 0;
    for index in 0..256 {
        let mut terms = Vec::new();
        for _ in 0..(2 + next() % 7) {
            let numerator = numerators[next() % numerators.len()];
            let denominator = denominators[next() % denominators.len()];
            let coefficient = format!("{}/{}", 1 + next() % 13, 1 + next() % 7);
            let phase = if next() % 3 == 0 { "𝑖" } else { "1" };
            terms.push(format!(
                "({coefficient})*{phase}*{numerator}/({denominator})"
            ));
        }
        if index % 5 == 0 {
            terms.extend(["x/(1-ep)".into(), "2*x/(2*ep-2)".into()]);
        }
        if index % 32 == 0 {
            terms.push("0.125*z/(1-ep)".into());
        }
        admitted += usize::from(check_grouping(vk_parse!(&terms.join("+")).unwrap()));
    }
    assert!(admitted >= 100, "only {admitted} cases exercised grouping");
}

#[test]
fn denominator_groups_are_used_by_large_master_combinations() {
    crate::Vakint::initialize_vakint_symbols();
    let expression = Atom::add_many((0..24).map(|index| {
        vk_parse!(&format!(
            "RustRedMaster(f,{index})*(x+ep*y)/(1-ep)^2+RustRedMaster(g,{index})*z/(1-ep)"
        ))
        .unwrap()
    }));
    assert!(expression.expand().as_view().nterms() >= 16);
    assert!(super::together_over_shared_remainder(&expression.expand()).is_err());
    assert_eq!(
        super::together_over_denominator_groups(&expression.expand()).unwrap(),
        proxy_tests::reference_combine_exact(expression.clone())
    );
    assert!(check_grouping(expression));
    let coefficient = vk_parse!(
        "123456789012345678901234567890123456789/7+987654321098765432109876543210987654321*𝑖/11"
    )
    .unwrap();
    assert!(matches!(
        coefficient.as_view(),
        AtomView::Num(number) if matches!(number.get_coeff_view(), CoefficientView::Large(..))
    ));
    assert!(
        ExactDenominatorGroups::raw_variables(coefficient.as_view(), true)
            .unwrap()
            .is_empty()
    );
    assert!(ExactDenominatorGroups::raw_variables(coefficient.as_view(), false).is_err());
    let expression = Atom::add_many((0..24).map(|index| {
        &coefficient
            * vk_parse!(&format!(
                "RustRedMaster(f,{index})/(1-ep)+g({index})/(2-2*ep)"
            ))
            .unwrap()
    }));
    assert!(
        !ExactDenominatorGroups::new(&expression.expand())
            .unwrap()
            .phase
            .is_one()
    );
    assert!(check_grouping(expression));
    assert!(ExactDenominatorGroups::new(&vk_parse!("x+y+z").unwrap()).is_err());
}

#[test]
fn unsupported_embedded_coefficients_keep_the_original_combination() {
    use symbolica::{
        coefficient::Coefficient,
        domains::{
            finite_field::{FiniteFieldCore, Zp64},
            integer::Z,
            rational::Q,
        },
    };

    crate::Vakint::initialize_vakint_symbols();
    let coefficient = Atom::num(Coefficient::RationalPolynomial(
        vk_parse!("1+hidden_parameter")
            .unwrap()
            .to_rational_polynomial::<_, _, u16>(&Q, &Z, None),
    ));
    let expression = Atom::add_many(
        (0..20).map(|index| &coefficient * vk_parse!(&format!("f({index})/(1-ep)")).unwrap()),
    );
    assert!(matches!(
        ExactDenominatorGroups::new(&expression.expand()),
        Err("unsupported_group_conversion")
    ));
    assert_eq!(
        combine_exact(expression.clone()),
        proxy_tests::reference_combine_exact(expression)
    );
    let field = Zp64::new(17);
    let field_element = field.to_element(5);
    for coefficient in [
        coefficient,
        Atom::num(Coefficient::from_finite_field(field, field_element)),
    ] {
        // Multiplying unsupported coefficients by i can yield Indeterminate;
        // phase admission must reject them before performing that arithmetic.
        let mixed =
            &coefficient * vk_parse!("x/(1-ep)").unwrap() + vk_parse!("𝑖*y/(1-ep)").unwrap();
        assert!(ExactDenominatorGroups::common_phase(&mixed.expand()).is_one());
        for complex in [false, true] {
            assert!(ExactDenominatorGroups::raw_variables(coefficient.as_view(), complex).is_err());
        }
        let denominator = (&coefficient + vk_parse!("ep").unwrap()).pow(-2);
        // Test the denominator-key boundary directly: expanding a mixed-field
        // finite-field fixture performs unsupported arithmetic before grouping.
        assert!(matches!(
            ExactDenominatorGroups::canonicalize(&denominator, None, false),
            Err("unsupported_group_conversion")
        ));
    }
}

#[test]
#[ignore = "replay locally captured complete exact-combination inputs without reductions"]
fn denominator_groups_replay_complete_captured_inputs() {
    crate::Vakint::initialize_vakint_symbols();
    let directory = std::env::var("VAKINT_COMBINE_CAPTURE_DIR")
        .expect("set VAKINT_COMBINE_CAPTURE_DIR to a directory of canonical .expr payloads");
    let mut paths = std::fs::read_dir(directory)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| {
            path.extension()
                .is_some_and(|extension| extension == "expr")
        })
        .collect::<Vec<_>>();
    paths.sort();
    assert!(!paths.is_empty());
    let mut admitted = 0;
    let mut phase_counts = [0usize; 3];
    let mut phase_seconds = [0f64; 3];
    for path in &paths {
        let input = std::fs::read_to_string(path).unwrap();
        let input = Atom::parse(&input, "combine_probe", Default::default()).unwrap();
        let expanded = input.expand();
        let started = std::time::Instant::now();
        let class = if !expanded.as_view().has_complex_coefficients() {
            0
        } else if ExactDenominatorGroups::common_phase(&expanded).is_one() {
            2
        } else {
            1
        };
        phase_seconds[class] += started.elapsed().as_secs_f64();
        phase_counts[class] += 1;
        admitted += usize::from(check_grouping(input));
    }
    assert!(admitted > 0);
    eprintln!(
        "Exactly replayed {} captured inputs; grouped {admitted}",
        paths.len()
    );
    eprintln!(
        "Phase admission (real/admitted/rejected): counts={phase_counts:?}, seconds={phase_seconds:?}"
    );
}
