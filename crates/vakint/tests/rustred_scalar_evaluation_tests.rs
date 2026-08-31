mod test_utils;

use std::collections::{BTreeSet, HashMap};
use symbolica::atom::{Atom, AtomCore};
use test_utils::{
    RustRedParityPolicy, TensorPrepass, alphaloop_matad_rustred_lanes,
    analytic_matad_rustred_lanes, compare_evaluations,
};
use vakint::{
    EvaluationMethod, EvaluationOrder, MATADOptions, RustRedEvaluationOptions, Vakint, VakintError,
    VakintSettings, params_from_f64, vakint_parse,
};

fn rustred_settings(substitute_masters: bool) -> VakintSettings {
    VakintSettings {
        form_exe_path: "/this/path/must/not/be/invoked/by/rustred-scalar".to_owned(),
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions { substitute_masters },
        )]),
        run_time_decimal_precision: 32,
        use_dot_product_notation: true,
        ..VakintSettings::default()
    }
}

fn matad_raw_settings() -> VakintSettings {
    VakintSettings {
        evaluation_order: EvaluationOrder::matad_only(Some(MATADOptions {
            expand_masters: false,
            ..MATADOptions::default()
        })),
        run_time_decimal_precision: 32,
        use_dot_product_notation: true,
        ..VakintSettings::default()
    }
}

fn matad_substituted_settings() -> VakintSettings {
    VakintSettings {
        evaluation_order: EvaluationOrder::matad_only(None),
        run_time_decimal_precision: 32,
        use_dot_product_notation: true,
        ..VakintSettings::default()
    }
}

fn assert_exactly_equal(left: Atom, right: Atom) {
    let difference = (left - right).together().expand();
    assert!(
        difference.is_zero(),
        "nonzero exact difference: {difference}"
    );
}

#[test]
fn acceptance_lane_policy_separates_exact_basis_from_numerical_parity() {
    for make_lanes in [alphaloop_matad_rustred_lanes, analytic_matad_rustred_lanes] {
        let exact = make_lanes(RustRedParityPolicy::ExactMatadBasis);
        let numerical = make_lanes(RustRedParityPolicy::NumericalOnly);

        assert_eq!(exact[0].rustred_parity(), None);
        assert_eq!(exact[1].rustred_parity(), None);
        assert_eq!(
            exact[2].rustred_parity(),
            Some(RustRedParityPolicy::ExactMatadBasis)
        );
        assert_eq!(
            numerical[2].rustred_parity(),
            Some(RustRedParityPolicy::NumericalOnly)
        );
        assert!(exact[2].forbids_form_after_prepass());
        assert!(numerical[2].forbids_form_after_prepass());
    }
}

#[test]
fn scalar_rustred_reduces_without_a_form_executable() {
    let vakint = Vakint::new().unwrap();
    let result = vakint
        .evaluate_integral(
            &rustred_settings(false),
            vakint_parse!("topo(I1L(muvsq,4))").unwrap().as_view(),
        )
        .unwrap();

    assert!(result.to_canonical_string().contains("Gam(1,1)"));

    let settings = rustred_settings(false);
    let scalar_identities = [
        (
            vakint_parse!("dot(k(1),k(1))*topo(I1L(muvsq,2))").unwrap(),
            vakint_parse!("topo(I1L(muvsq,1))+muvsq*topo(I1L(muvsq,2))").unwrap(),
        ),
        (
            vakint_parse!("dot(k(1),k(2))*topo(I2L(muvsq,1,1,1))").unwrap(),
            vakint_parse!(
                "1/2*(\
                    topo(I2L(muvsq,1,1,0))\
                   -topo(I2L(muvsq,0,1,1))\
                   -topo(I2L(muvsq,1,0,1))\
                   -muvsq*topo(I2L(muvsq,1,1,1))\
                )"
            )
            .unwrap(),
        ),
    ];
    for (scalar_numerator, affine_expansion) in scalar_identities {
        let lowered = vakint
            .evaluate_integral(&settings, scalar_numerator.as_view())
            .unwrap();
        let expanded = vakint
            .evaluate_integral(&settings, affine_expansion.as_view())
            .unwrap();
        assert_exactly_equal(lowered, expanded);
    }
}

#[test]
fn scalar_rustred_matches_raw_matad_through_two_loops() {
    let vakint = Vakint::new().unwrap();
    let rustred = rustred_settings(false);
    let matad = matad_raw_settings();
    let mut inputs = (1..=6)
        .map(|power| vakint_parse!(format!("topo(I1L(muvsq,{power}))").as_str()).unwrap())
        .collect::<Vec<_>>();
    inputs.extend([
        vakint_parse!("topo(I2L(muvsq,1,1,1))").unwrap(),
        vakint_parse!("topo(I2L(muvsq,2,1,1))").unwrap(),
        vakint_parse!("topo(I2L(muvsq,1,2,1))").unwrap(),
        vakint_parse!("topo(I2L(muvsq,1,1,2))").unwrap(),
        vakint_parse!("topo(I2L(muvsq,2,2,1))").unwrap(),
        vakint_parse!("topo(I2L(muvsq,3,1,1))").unwrap(),
        vakint_parse!("topo(I2L(muvsq,-1,2,2))").unwrap(),
        vakint_parse!("topo(I2L(muvsq,0,1,1))").unwrap(),
        vakint_parse!("topo(I2L(muvsq,0,2,1))").unwrap(),
        vakint_parse!("topo(I2L(muvsq,2,0,1))").unwrap(),
        vakint_parse!("topo(I2L(muvsq,2,1,0))").unwrap(),
        vakint_parse!("topo(I2L_pinch_3(muvsq,1,1,0))").unwrap(),
        vakint_parse!("topo(I2L_pinch_3(muvsq,2,1,0))").unwrap(),
    ]);

    for input in inputs {
        let rustred_result = vakint.evaluate_integral(&rustred, input.as_view()).unwrap();
        let matad_result = vakint.evaluate_integral(&matad, input.as_view()).unwrap();
        assert_exactly_equal(rustred_result, matad_result);
    }
}

#[test]
fn scalar_rustred_reuses_matcher_routing_and_preserves_spectators() {
    let vakint = Vakint::new().unwrap();
    let settings = rustred_settings(false);
    let equivalent_inputs = [
        (
            vakint_parse!(
                "7*user_space::c*topo(\
                    prop(77,edge(42,42),k(8),user_space::mass_squared,4)\
                )"
            )
            .unwrap(),
            vakint_parse!("7*user_space::c*topo(I1L(user_space::mass_squared,4))").unwrap(),
        ),
        (
            vakint_parse!(
                "user_space::c*topo(\
                    prop(55,edge(7,10),k(22),user_space::mass_squared,1)*\
                    prop(9,edge(7,10),k(11),user_space::mass_squared,2)*\
                    prop(33,edge(10,7),k(11)+k(22),user_space::mass_squared,1)\
                )"
            )
            .unwrap(),
            vakint_parse!("user_space::c*topo(I2L(user_space::mass_squared,2,1,1))").unwrap(),
        ),
        (
            vakint_parse!(
                "user_space::c*topo(\
                    prop(33,edge(10,10),k(22),user_space::mass_squared,2)*\
                    prop(55,edge(10,10),k(11),user_space::mass_squared,1)\
                )"
            )
            .unwrap(),
            vakint_parse!("user_space::c*topo(I2L_pinch_3(user_space::mass_squared,2,1,0))")
                .unwrap(),
        ),
    ];

    for (graph_input, short_input) in equivalent_inputs {
        let graph_result = vakint
            .evaluate_integral(&settings, graph_input.as_view())
            .unwrap();
        let short_result = vakint
            .evaluate_integral(&settings, short_input.as_view())
            .unwrap();
        assert_exactly_equal(graph_result, short_result);
    }
}

#[test]
fn scalar_rustred_accepts_literal_unit_mass_without_a_mass_symbol() {
    let vakint = Vakint::new().unwrap();
    let settings = rustred_settings(false);
    let cases = [
        (
            vakint_parse!("topo(I1L(1,4))").unwrap(),
            vakint_parse!("topo(I1L(user_space::mass_squared,4))").unwrap(),
        ),
        (
            vakint_parse!("topo(I2L(1,2,1,1))").unwrap(),
            vakint_parse!("topo(I2L(user_space::mass_squared,2,1,1))").unwrap(),
        ),
        (
            vakint_parse!("topo(I2L_pinch_3(1,2,1,0))").unwrap(),
            vakint_parse!("topo(I2L_pinch_3(user_space::mass_squared,2,1,0))").unwrap(),
        ),
    ];

    for (unit_mass, symbolic_mass) in cases {
        let unit_result = vakint
            .evaluate_integral(&settings, unit_mass.as_view())
            .unwrap();
        let symbolic_result = vakint
            .evaluate_integral(&settings, symbolic_mass.as_view())
            .unwrap()
            .replace(
                vakint_parse!("user_space::mass_squared")
                    .unwrap()
                    .to_pattern(),
            )
            .with(Atom::num(1).to_pattern());
        assert_exactly_equal(unit_result, symbolic_result);
    }

    let caller_m_result = vakint
        .evaluate_integral(
            &settings,
            vakint_parse!("topo(I1L(M,4))").unwrap().as_view(),
        )
        .unwrap();
    let reference = vakint
        .evaluate_integral(
            &settings,
            vakint_parse!("topo(I1L(user_space::mass_squared,4))")
                .unwrap()
                .as_view(),
        )
        .unwrap()
        .replace(
            vakint_parse!("user_space::mass_squared")
                .unwrap()
                .to_pattern(),
        )
        .with(vakint_parse!("M").unwrap().to_pattern());
    assert_exactly_equal(caller_m_result, reference);
}

#[test]
fn scalar_rustred_matches_substituted_matad_through_two_loops() {
    let vakint = Vakint::new().unwrap();
    let rustred = rustred_settings(true);
    let matad = matad_substituted_settings();
    let parameter_values = HashMap::from([("muvsq".to_owned(), 1.0), ("mursq".to_owned(), 1.0)]);
    let rustred_parameters = vakint.params_from_f64(&rustred, &parameter_values);
    let matad_parameters = vakint.params_from_f64(&matad, &parameter_values);
    let inputs = [
        vakint_parse!("topo(I1L(muvsq,4))").unwrap(),
        vakint_parse!("topo(I2L(muvsq,2,2,1))").unwrap(),
        vakint_parse!("topo(I2L_pinch_3(muvsq,2,1,0))").unwrap(),
    ];

    for input in inputs {
        let rustred_result = vakint.evaluate_integral(&rustred, input.as_view()).unwrap();
        let matad_result = vakint.evaluate_integral(&matad, input.as_view()).unwrap();
        let (rustred_numerical, _) = Vakint::full_numerical_evaluation(
            &rustred,
            rustred_result.as_view(),
            &rustred_parameters,
            &HashMap::default(),
            None,
        )
        .unwrap();
        let (matad_numerical, _) = Vakint::full_numerical_evaluation(
            &matad,
            matad_result.as_view(),
            &matad_parameters,
            &HashMap::default(),
            None,
        )
        .unwrap();
        let (matches, detail) =
            rustred_numerical.does_approx_match(&matad_numerical, None, 1.0e-25, 0.0);
        assert!(matches, "{detail}");
    }
}

#[test]
fn form_tensor_prepass_then_rustred_scalar_matches_matad() {
    let vakint = Vakint::new().unwrap();
    let mut rustred = rustred_settings(false);
    rustred.form_exe_path = "form".to_owned();
    let matad = matad_raw_settings();
    let inputs = [
        vakint_parse!(
            "(k(1,1)*k(1,2)+k(1,3)*p(1,3))*topo(\
                prop(1,edge(1,1),k(1),muvsq,1)\
            )"
        )
        .unwrap(),
        vakint_parse!(
            "(k(1,11)*p(1,11)*k(1,12)*p(1,12))*topo(\
                prop(1,edge(1,1),k(1),muvsq,2)\
            )"
        )
        .unwrap(),
        vakint_parse!(
            "(user_space::A*k(1,11)*p(1,11)*k(1,12)*p(1,12)+user_space::B)*topo(\
                prop(1,edge(1,1),k(1),muvsq,2)\
            )"
        )
        .unwrap(),
        vakint_parse!(
            "(k(1,1)*p(1,1)*k(1,2)*p(2,2))*topo(\
                prop(1,edge(1,1),k(1),muvsq,1)\
            )"
        )
        .unwrap(),
        vakint_parse!(
            "(\
                k(1,11)*k(2,22)*k(1,11)*k(2,22)\
                +p(1,11)*k(1,11)*k(1,22)*p(1,22)\
                +p(1,11)*p(2,11)*k(2,22)*k(2,22)\
            )*topo(\
                prop(1,edge(1,2),k(1),muvsq,1)*\
                prop(2,edge(1,2),k(2),muvsq,1)*\
                prop(3,edge(2,1),k(1)+k(2),muvsq,1)\
            )"
        )
        .unwrap(),
        vakint_parse!(
            "k(11,mu)*k(22,nu)*topo(\
                prop(33,edge(10,10),k(22),muvsq,1)*\
                prop(55,edge(10,10),k(11),muvsq,1)\
            )"
        )
        .unwrap(),
    ];

    for input in inputs {
        let rustred_result = vakint.evaluate(&rustred, input.as_view()).unwrap();
        let matad_result = vakint.evaluate(&matad, input.as_view()).unwrap();
        assert_exactly_equal(rustred_result, matad_result);
    }
}

#[derive(Clone, Copy)]
enum K6PeerSuite {
    AlphaLoopMatad,
    AnalyticMatad,
}

impl K6PeerSuite {
    fn lanes(self, parity: RustRedParityPolicy) -> [test_utils::EvaluationTestLane; 3] {
        match self {
            Self::AlphaLoopMatad => alphaloop_matad_rustred_lanes(parity),
            Self::AnalyticMatad => analytic_matad_rustred_lanes(parity),
        }
    }
}

#[derive(Clone, Copy)]
struct PendingK6Acceptance {
    legacy_test: &'static str,
    tensor_prepass: bool,
    parity: RustRedParityPolicy,
    peer_suite: K6PeerSuite,
}

const PENDING_K6_CORE_ACCEPTANCE: [PendingK6Acceptance; 11] = [
    PendingK6Acceptance {
        legacy_test: "integral_alphaloop_vs_matad_tests::test_integrate_3l_basketball_a",
        tensor_prepass: true,
        parity: RustRedParityPolicy::NumericalOnly,
        peer_suite: K6PeerSuite::AlphaLoopMatad,
    },
    PendingK6Acceptance {
        legacy_test: "integral_alphaloop_vs_matad_tests::test_integrate_3l_basketball_b",
        tensor_prepass: true,
        parity: RustRedParityPolicy::NumericalOnly,
        peer_suite: K6PeerSuite::AlphaLoopMatad,
    },
    PendingK6Acceptance {
        legacy_test: "integral_alphaloop_vs_matad_tests::test_integrate_3l_no_numerator",
        tensor_prepass: false,
        parity: RustRedParityPolicy::NumericalOnly,
        peer_suite: K6PeerSuite::AlphaLoopMatad,
    },
    PendingK6Acceptance {
        legacy_test: "integral_alphaloop_vs_matad_tests::test_integrate_3l_rank_4",
        tensor_prepass: true,
        parity: RustRedParityPolicy::NumericalOnly,
        peer_suite: K6PeerSuite::AlphaLoopMatad,
    },
    PendingK6Acceptance {
        legacy_test: "integral_alphaloop_vs_matad_tests::test_integrate_3l_rank_4_different_scales",
        tensor_prepass: true,
        parity: RustRedParityPolicy::NumericalOnly,
        peer_suite: K6PeerSuite::AlphaLoopMatad,
    },
    PendingK6Acceptance {
        legacy_test: "integral_evaluation_analytic_tests::test_integrate_3l",
        tensor_prepass: false,
        parity: RustRedParityPolicy::NumericalOnly,
        peer_suite: K6PeerSuite::AnalyticMatad,
    },
    PendingK6Acceptance {
        legacy_test: "integral_evaluation_analytic_tests::test_integrate_3l_rank_4",
        tensor_prepass: true,
        parity: RustRedParityPolicy::NumericalOnly,
        peer_suite: K6PeerSuite::AnalyticMatad,
    },
    PendingK6Acceptance {
        legacy_test: "integral_evaluation_analytic_tests::test_integrate_3l_rank_4_additional_symbols_numerator",
        tensor_prepass: true,
        parity: RustRedParityPolicy::NumericalOnly,
        peer_suite: K6PeerSuite::AnalyticMatad,
    },
    PendingK6Acceptance {
        legacy_test: "integral_evaluation_analytic_tests::test_integrate_3l_rank_4_matad",
        tensor_prepass: true,
        parity: RustRedParityPolicy::NumericalOnly,
        peer_suite: K6PeerSuite::AnalyticMatad,
    },
    PendingK6Acceptance {
        legacy_test: "integral_evaluation_analytic_tests::test_integrate_3l_rank_4_matad_additional_symbols_numerator",
        tensor_prepass: true,
        parity: RustRedParityPolicy::NumericalOnly,
        peer_suite: K6PeerSuite::AnalyticMatad,
    },
    PendingK6Acceptance {
        legacy_test: "integral_evaluation_analytic_tests::test_integrate_3l_matad",
        tensor_prepass: false,
        parity: RustRedParityPolicy::NumericalOnly,
        peer_suite: K6PeerSuite::AnalyticMatad,
    },
];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SupplementalThreeLoopRole {
    OptionalPySecDec,
    TopologyMatching,
    FrozenExperimentalTensorMode,
    MatadK6Oracle,
}

#[derive(Clone, Copy)]
struct SupplementalThreeLoopTest {
    test: &'static str,
    role: SupplementalThreeLoopRole,
}

// These are useful coverage, but they are not part of the 11-case gating
// AlphaLoop/MATAD numerical-parity matrix above. PySecDec remains non-gating;
// matching and tensor-mode tests exercise different layers; the MATAD records
// are offline closure diagnostics rather than RustRed acceptance.
const SUPPLEMENTAL_THREE_LOOP_TESTS: [SupplementalThreeLoopTest; 8] = [
    SupplementalThreeLoopTest {
        test: "integral_comparison_vs_pysecdec_tests::test_integrate_3l_pysecdec",
        role: SupplementalThreeLoopRole::OptionalPySecDec,
    },
    SupplementalThreeLoopTest {
        test: "integral_comparison_vs_pysecdec_tests::test_integrate_3l_rank_4",
        role: SupplementalThreeLoopRole::OptionalPySecDec,
    },
    SupplementalThreeLoopTest {
        test: "integral_comparison_vs_pysecdec_tests::test_integrate_3l_rank_4_matad",
        role: SupplementalThreeLoopRole::OptionalPySecDec,
    },
    SupplementalThreeLoopTest {
        test: "integral_evaluation_pysecdec_tests::test_integrate_3l_o_eps",
        role: SupplementalThreeLoopRole::OptionalPySecDec,
    },
    SupplementalThreeLoopTest {
        test: "input_matching_tests::test_3l_matching_with_zero_powers_in_short_form",
        role: SupplementalThreeLoopRole::TopologyMatching,
    },
    SupplementalThreeLoopTest {
        test: "tensor_reduction_mode_tests::rustred_projects_an_isp_completed_three_loop_pinch_without_form",
        role: SupplementalThreeLoopRole::FrozenExperimentalTensorMode,
    },
    SupplementalThreeLoopTest {
        test: "k6_matad_oracle_tests::matad_k6_exact_raw_master_records",
        role: SupplementalThreeLoopRole::MatadK6Oracle,
    },
    SupplementalThreeLoopTest {
        test: "k6_matad_oracle_tests::matad_k6_representative_laurent_records",
        role: SupplementalThreeLoopRole::MatadK6Oracle,
    },
];

/// The five matcher-owned graph classes that one sector-complete K=6 artifact
/// must cover. This is deliberately typed separately from legacy test names so
/// a renamed test cannot silently erase a topology-class obligation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum K6MatcherClass {
    Parent,
    Pinch6,
    Pinch3And6,
    Pinch1And6,
    Pinch1And3And6,
}

impl K6MatcherClass {
    const ALL: [Self; 5] = [
        Self::Parent,
        Self::Pinch6,
        Self::Pinch3And6,
        Self::Pinch1And6,
        Self::Pinch1And3And6,
    ];

    fn short_head(self) -> &'static str {
        match self {
            Self::Parent => "I3L(",
            Self::Pinch6 => "I3L_pinch_6(",
            Self::Pinch3And6 => "I3L_pinch_3_6(",
            Self::Pinch1And6 => "I3L_pinch_1_6(",
            Self::Pinch1And3And6 => "I3L_pinch_1_3_6(",
        }
    }

    fn scalar_fixture(self) -> Atom {
        match self {
            Self::Parent => vakint_parse!("topo(I3L(muvsq,1,1,1,1,1,1))").unwrap(),
            Self::Pinch6 => vakint_parse!("topo(I3L_pinch_6(muvsq,1,1,1,1,1,0))").unwrap(),
            Self::Pinch3And6 => vakint_parse!("topo(I3L_pinch_3_6(muvsq,1,1,0,1,1,0))").unwrap(),
            Self::Pinch1And6 => vakint_parse!("topo(I3L_pinch_1_6(muvsq,0,1,1,1,1,0))").unwrap(),
            Self::Pinch1And3And6 => {
                vakint_parse!("topo(I3L_pinch_1_3_6(muvsq,0,1,0,1,1,0))").unwrap()
            }
        }
    }
}

#[test]
fn k6_pending_fixtures_cover_all_registered_matcher_classes() {
    let vakint = Vakint::new().unwrap();
    let mut heads = BTreeSet::new();
    for matcher_class in K6MatcherClass::ALL {
        assert!(heads.insert(matcher_class.short_head()));
        let canonical = vakint
            .to_canonical(
                &VakintSettings::default(),
                matcher_class.scalar_fixture().as_view(),
                true,
            )
            .unwrap_or_else(|error| {
                panic!(
                    "registered matcher class {:?} rejected its fixture: {error}",
                    matcher_class
                )
            });
        assert!(
            canonical
                .to_canonical_string()
                .contains(matcher_class.short_head()),
            "fixture for {:?} routed to the wrong short class: {canonical}",
            matcher_class
        );
    }
    assert_eq!(heads.len(), 5);
}

fn run_pending_k6_numerical_parity(matcher_class: K6MatcherClass) {
    let raw_masses = [("muvsq".to_owned(), 1.0), ("mursq".to_owned(), 1.0)]
        .into_iter()
        .collect::<ahash::HashMap<_, _>>();
    let numerical_masses = params_from_f64(&raw_masses, 32);
    compare_evaluations(
        VakintSettings {
            number_of_terms_in_epsilon_expansion: 4,
            run_time_decimal_precision: 32,
            ..VakintSettings::default()
        },
        &alphaloop_matad_rustred_lanes(RustRedParityPolicy::NumericalOnly),
        TensorPrepass::Skip,
        matcher_class.scalar_fixture().as_view(),
        numerical_masses,
        ahash::HashMap::default(),
        1.0e-25,
        0.0,
        true,
    );
}

macro_rules! pending_k6_numerical_parity_test {
    ($name:ident, $matcher_class:expr) => {
        #[test]
        #[ignore = "pending certified sector-complete K=6 artifact"]
        fn $name() {
            run_pending_k6_numerical_parity($matcher_class);
        }
    };
}

pending_k6_numerical_parity_test!(pending_k6_numerical_parity_parent, K6MatcherClass::Parent);
pending_k6_numerical_parity_test!(pending_k6_numerical_parity_pinch_6, K6MatcherClass::Pinch6);
pending_k6_numerical_parity_test!(
    pending_k6_numerical_parity_pinch_3_6,
    K6MatcherClass::Pinch3And6
);
pending_k6_numerical_parity_test!(
    pending_k6_numerical_parity_pinch_1_6,
    K6MatcherClass::Pinch1And6
);
pending_k6_numerical_parity_test!(
    pending_k6_numerical_parity_pinch_1_3_6,
    K6MatcherClass::Pinch1And3And6
);

#[test]
fn k6_core_acceptance_inventory_is_explicitly_pending() {
    let names = PENDING_K6_CORE_ACCEPTANCE
        .iter()
        .map(|case| case.legacy_test)
        .collect::<BTreeSet<_>>();
    assert_eq!(names.len(), 11, "pending K6 cases must remain unique");
    assert!(
        PENDING_K6_CORE_ACCEPTANCE
            .iter()
            .all(|case| case.parity == RustRedParityPolicy::NumericalOnly)
    );
    for case in PENDING_K6_CORE_ACCEPTANCE {
        let lanes = case.peer_suite.lanes(case.parity);
        assert_eq!(
            lanes[2].rustred_parity(),
            Some(RustRedParityPolicy::NumericalOnly),
            "{} must be promoted as a numerical-only RustRed lane",
            case.legacy_test
        );
        assert!(lanes[2].forbids_form_after_prepass());
    }
    assert_eq!(
        PENDING_K6_CORE_ACCEPTANCE
            .iter()
            .filter(|case| case.tensor_prepass)
            .count(),
        8
    );
    assert_eq!(
        names
            .iter()
            .filter(|name| name.starts_with("integral_alphaloop_vs_matad_tests::"))
            .count(),
        5
    );
    assert_eq!(
        names
            .iter()
            .filter(|name| name.starts_with("integral_evaluation_analytic_tests::"))
            .count(),
        6
    );

    let supplemental_names = SUPPLEMENTAL_THREE_LOOP_TESTS
        .iter()
        .map(|case| case.test)
        .collect::<BTreeSet<_>>();
    assert_eq!(supplemental_names.len(), 8);
    for (role, expected) in [
        (SupplementalThreeLoopRole::OptionalPySecDec, 4),
        (SupplementalThreeLoopRole::TopologyMatching, 1),
        (SupplementalThreeLoopRole::FrozenExperimentalTensorMode, 1),
        (SupplementalThreeLoopRole::MatadK6Oracle, 2),
    ] {
        assert_eq!(
            SUPPLEMENTAL_THREE_LOOP_TESTS
                .iter()
                .filter(|case| case.role == role)
                .count(),
            expected
        );
    }

    // This assertion deliberately flips when the shipped K=6 artifact and
    // adapter arrive, forcing this inventory to become live peer-lane tests.
    let vakint = Vakint::new().unwrap();
    let result = vakint.evaluate_integral(
        &rustred_settings(false),
        vakint_parse!("topo(I3L(muvsq,1,1,1,1,1,1))")
            .unwrap()
            .as_view(),
    );
    assert!(
        matches!(result, Err(VakintError::NoEvaluationMethodFound(_, _))),
        "K6 support changed; promote or debug every pending numerical-parity case: {result:?}"
    );
}
