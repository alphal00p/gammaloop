mod test_utils;

use test_utils::compare_two_evaluations;
use vakint::{EvaluationMethod, EvaluationOrder, PySecDecOptions, VakintSettings};

use vakint::{externals_from_f64, params_from_f64, vakint_parse};

const N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC: u32 = 16;
const COMPARISON_WITH_PYSECDEC_REL_THRESHOLD: f64 = 1.0e-7;
// PySecDec QMC is often very optimistic
const MAX_PULL: f64 = 1.0e5;

#[test_log::test]
fn test_integrate_1l_pysecdec() {
    if test_utils::should_skip_pysecdec_tests() {
        return;
    }
    #[rustfmt::skip]
    compare_two_evaluations(
        VakintSettings::default(),
        ((&EvaluationOrder::analytic_only() ,true),
        (&EvaluationOrder::pysecdec_only(Some(PySecDecOptions { reuse_existing_output: Some("./tests_workspace/pysecdec_comparison_1l_pysecdec".into()),..PySecDecOptions::default() })) ,false)),
        vakint_parse!(
            "(1)*topo(\
            prop(1,edge(1,1),k(1),muvsq,1)\
        )"
        ).unwrap().as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
        COMPARISON_WITH_PYSECDEC_REL_THRESHOLD, MAX_PULL,
        true,
    );
}

#[test_log::test]
fn test_integrate_1l_pysecdec_non_unit_mass() {
    if test_utils::should_skip_pysecdec_tests() {
        return;
    }
    #[rustfmt::skip]
    compare_two_evaluations(
        VakintSettings::default(),
        ((&EvaluationOrder::analytic_only() ,true),
        (&EvaluationOrder::pysecdec_only(Some(PySecDecOptions { reuse_existing_output: Some("./tests_workspace/pysecdec_comparison_1l_pysecdec_non_unit_mass".into()),..PySecDecOptions::default() })) ,false)),
        vakint_parse!(
            "(1)*topo(\
            prop(1,edge(1,1),k(1),muvsq,1)\
        )"
        ).unwrap().as_view(),
        params_from_f64(&[("muvsq".into(), 2.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
        COMPARISON_WITH_PYSECDEC_REL_THRESHOLD, MAX_PULL,
        true,
    );
}

#[test_log::test]
fn test_integrate_1l_pysecdec_non_unit_scale() {
    if test_utils::should_skip_pysecdec_tests() {
        return;
    }
    #[rustfmt::skip]
    compare_two_evaluations(
        VakintSettings::default(),
        ((&EvaluationOrder::analytic_only() ,true),
        (&EvaluationOrder::pysecdec_only(Some(PySecDecOptions { reuse_existing_output: Some("./tests_workspace/pysecdec_comparison_1l_pysecdec_non_unit_scale".into()),..PySecDecOptions::default() })) ,false)),
        vakint_parse!(
            "(1)*topo(\
            prop(1,edge(1,1),k(1),muvsq,1)\
        )"
        ).unwrap().as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 2.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
        COMPARISON_WITH_PYSECDEC_REL_THRESHOLD, MAX_PULL,
        true,
    );
}

#[test_log::test]
fn test_integrate_1l_pysecdec_num_rank_two() {
    if test_utils::should_skip_pysecdec_tests() {
        return;
    }
    #[rustfmt::skip]
    compare_two_evaluations(
        VakintSettings::default(),
        ((&EvaluationOrder::analytic_only() ,true),
        (&EvaluationOrder::pysecdec_only(Some(PySecDecOptions { reuse_existing_output: Some("./tests_workspace/pysecdec_comparison_1l_pysecdec_num_rank_two".into()),..PySecDecOptions::default() })) ,false)),
        vakint_parse!(
            "((k(1,33)*k(1,33))^2+k(1,55)*p(1,55))*topo(\
            prop(1,edge(1,1),k(1),muvsq,1)\
        )"
        ).unwrap().as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 2.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
        COMPARISON_WITH_PYSECDEC_REL_THRESHOLD, MAX_PULL,
        true,
    );
}

#[test_log::test]
fn test_integrate_1l_pysecdec_dot_product_external() {
    if test_utils::should_skip_pysecdec_tests() {
        return;
    }
    #[rustfmt::skip]
    compare_two_evaluations(
        VakintSettings::default(),
        ((&EvaluationOrder::analytic_only() ,true),
        (&EvaluationOrder::pysecdec_only(Some(PySecDecOptions { reuse_existing_output: Some("./tests_workspace/pysecdec_comparison_1l_pysecdec_dot_product_external".into()),..PySecDecOptions::default() })) ,false)),
        vakint_parse!(
            "(k(1,1)*p(1,1)*k(1,2)*p(2,2))*topo(\
                prop(1,edge(1,1),k(1),muvsq,1)\
            )"
        ).unwrap().as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
        COMPARISON_WITH_PYSECDEC_REL_THRESHOLD, MAX_PULL,
        true,
    );
}

#[test_log::test]
fn test_integrate_2l_pysecdec() {
    if test_utils::should_skip_pysecdec_tests() {
        return;
    }
    #[rustfmt::skip]
    compare_two_evaluations(
        VakintSettings::default(),
        ((&EvaluationOrder::analytic_only() ,true),
        (&EvaluationOrder::pysecdec_only(Some(PySecDecOptions { reuse_existing_output: Some("./tests_workspace/pysecdec_comparison_2l_pysecdec".into()),..PySecDecOptions::default() })) ,false)),
        vakint_parse!(
            "(1)*topo(\
                prop(1,edge(1,2),k(1),muvsq,1)\
                *prop(2,edge(1,2),k(2),muvsq,1)\
                *prop(3,edge(2,1),k(1)+k(2),muvsq,1)\
            )"
        ).unwrap().as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 2.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
        COMPARISON_WITH_PYSECDEC_REL_THRESHOLD, MAX_PULL,
        true,
    );
}

#[test_log::test]
fn test_integrate_2l_pysecdec_pinched() {
    if test_utils::should_skip_pysecdec_tests() {
        return;
    }
    #[rustfmt::skip]
    compare_two_evaluations(
        VakintSettings::default(),
        ((&EvaluationOrder::analytic_only() ,true),
        (&EvaluationOrder::pysecdec_only(Some(PySecDecOptions { reuse_existing_output: Some("./tests_workspace/pysecdec_comparison_2l_pysecdec_pinched".into()),..PySecDecOptions::default() })) ,false)),
        vakint_parse!(
            "(1)*topo(\
                prop(1,edge(1,1),k(1),muvsq,1)\
                *prop(2,edge(1,1),k(2),muvsq,1)
            )"
        ).unwrap().as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
        COMPARISON_WITH_PYSECDEC_REL_THRESHOLD, MAX_PULL,
        true,
    );
}

#[test_log::test]
fn test_integrate_2l_pysecdec_pinched_other_lmb() {
    if test_utils::should_skip_pysecdec_tests() {
        return;
    }
    #[rustfmt::skip]
    compare_two_evaluations(
        VakintSettings::default(),
        ((&EvaluationOrder::analytic_only() ,true),
        (&EvaluationOrder::pysecdec_only(Some(PySecDecOptions { reuse_existing_output: Some("./tests_workspace/pysecdec_comparison_2l_pysecdec_pinched_other_lmb".into()),..PySecDecOptions::default() })) ,false)),
        vakint_parse!(
            "(1)*topo(\
                prop(1,edge(1,1),k(1),muvsq,1)\
                *prop(3,edge(1,1),k(1)+k(2),muvsq,1)
            )"
        ).unwrap().as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
        COMPARISON_WITH_PYSECDEC_REL_THRESHOLD, MAX_PULL,
        true,
    );
}

#[test_log::test]
fn test_integrate_2l_pysecdec_rank_four_num() {
    if test_utils::should_skip_pysecdec_tests() {
        return;
    }
    #[rustfmt::skip]
    compare_two_evaluations(
        VakintSettings::default(),
        ((&EvaluationOrder::analytic_only() ,true),
        (&EvaluationOrder::pysecdec_only(Some(PySecDecOptions { reuse_existing_output: Some("./tests_workspace/pysecdec_comparison_2l_pysecdec_rank_four_num".into()),..PySecDecOptions::default() })) ,false)),
        vakint_parse!(
            "(
                  k(1,11)*k(2,22)*k(1,11)*k(2,22)
                + p(1,11)*k(1,11)*k(1,22)*p(1,22)
                + p(1,11)*p(2,11)*k(2,22)*k(2,22)
            )
            *topo(\
                  prop(1,edge(1,2),k(1),muvsq,1)\
                * prop(2,edge(1,2),k(2),muvsq,1)\
                * prop(3,edge(2,1),k(1)+k(2),muvsq,1)\
            )"
        ).unwrap().as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 2.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
        COMPARISON_WITH_PYSECDEC_REL_THRESHOLD, MAX_PULL,
        true,
    );
}

#[test_log::test]
fn test_integrate_3l_pysecdec() {
    if test_utils::should_skip_pysecdec_tests() {
        return;
    }
    #[rustfmt::skip]
    compare_two_evaluations(
        VakintSettings::default(),
        ((&EvaluationOrder::analytic_only() ,true),
        (&EvaluationOrder::pysecdec_only(Some(PySecDecOptions { reuse_existing_output: Some("./tests_workspace/pysecdec_comparison_3l_pysecdec".into()),..PySecDecOptions::default() })) ,false)),
        vakint_parse!(
            "( 1 )
            *topo(\
                 prop(1,edge(1,2),k(1),muvsq,1)\
                *prop(2,edge(2,3),k(2),muvsq,1)\
                *prop(3,edge(3,1),k(3),muvsq,1)\
                *prop(4,edge(1,4),k(3)-k(1),muvsq,1)\
                *prop(5,edge(2,4),k(1)-k(2),muvsq,1)\
                *prop(6,edge(3,4),k(2)-k(3),muvsq,1)\
            )"
        ).unwrap().as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 2.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
        COMPARISON_WITH_PYSECDEC_REL_THRESHOLD, MAX_PULL,
        true,
    );
}

#[test_log::test]
fn test_integrate_3l_rank_4() {
    if test_utils::should_skip_pysecdec_tests() {
        return;
    }
    // pySecDec is not so great for such higher rank cases, so we need to set a very high threshold
    const ADJUSTED_THRESHOLD: f64 = 1.0e-2;
    #[rustfmt::skip]
    compare_two_evaluations(
        VakintSettings::default(),
        ((&EvaluationOrder::analytic_only() ,true),
         (&EvaluationOrder(vec![EvaluationMethod::PySecDec(
            PySecDecOptions{ min_n_evals: 100_000, max_n_evals: 10_000_000, reuse_existing_output: Some("./tests_workspace/pysecdec_comparison_integrate_3l_rank_4".into()),..PySecDecOptions::default()} )]) ,true)),
        vakint_parse!(
            "(
                  k(1,11)*k(2,11)*k(1,22)*k(2,22)
                + p(1,11)*k(3,11)*k(3,22)*p(2,22)
                + p(1,11)*p(2,11)*(k(2,22)+k(1,22))*k(2,22)
             )
            *topo(\
                 prop(1,edge(1,2),k(1),muvsq,1)\
                *prop(2,edge(2,3),k(2),muvsq,1)\
                *prop(3,edge(3,1),k(3),muvsq,1)\
                *prop(4,edge(1,4),k(3)-k(1),muvsq,1)\
                *prop(5,edge(2,4),k(1)-k(2),muvsq,1)\
                *prop(6,edge(3,4),k(2)-k(3),muvsq,1)\
            )"
        ).unwrap().as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 2.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
            ADJUSTED_THRESHOLD, MAX_PULL,
        true,
    );
}

#[test_log::test]
fn test_integrate_3l_rank_4_matad() {
    if test_utils::should_skip_pysecdec_tests() {
        return;
    }
    // pySecDec is not so great for such higher rank cases, so we need to set a very high threshold
    const ADJUSTED_THRESHOLD: f64 = 1.0e-2;
    #[rustfmt::skip]
    compare_two_evaluations(
        VakintSettings {number_of_terms_in_epsilon_expansion: 5,..VakintSettings::default()},
        ((&EvaluationOrder::matad_only(None) ,true),
         (&EvaluationOrder(vec![EvaluationMethod::PySecDec(
            PySecDecOptions{
                // reuse_existing_output: Some("./SAVED_PYSECDEC_OUTPUTS/test_integrate_3l_rank_4_matad".into()),
                min_n_evals: 100_000, max_n_evals: 10_000_000, reuse_existing_output: Some("./tests_workspace/pysecdec_comparison_3l_rank_4_matad".into()), ..PySecDecOptions::default()} )]) ,true)),
        vakint_parse!(
            "(
                  k(1,11)*k(2,11)*k(1,22)*k(2,22)
                + p(1,11)*k(3,11)*k(3,22)*p(2,22)
                + p(1,11)*p(2,11)*(k(2,22)+k(1,22))*k(2,22)
             )
            *topo(\
                 prop(1,edge(1,2),k(1),muvsq,1)\
                *prop(2,edge(2,3),k(2),muvsq,1)\
                *prop(3,edge(3,1),k(3),muvsq,1)\
                *prop(4,edge(1,4),k(3)-k(1),muvsq,1)\
                *prop(5,edge(2,4),k(1)-k(2),muvsq,1)\
                *prop(6,edge(3,4),k(2)-k(3),muvsq,1)\
            )"
        ).unwrap().as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 2.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_COMPARISON_WITH_PYSECDEC),
            ADJUSTED_THRESHOLD, MAX_PULL,
        true,
    );
}

#[allow(dead_code)]
fn run_integral_comparison_vs_pysecdec_tests() {
    // Convenience runner to execute all tests in this module.
    test_integrate_1l_pysecdec();
    test_integrate_1l_pysecdec_non_unit_mass();
    test_integrate_1l_pysecdec_non_unit_scale();
    test_integrate_1l_pysecdec_num_rank_two();
    test_integrate_1l_pysecdec_dot_product_external();
    test_integrate_2l_pysecdec();
    test_integrate_2l_pysecdec_pinched();
    test_integrate_2l_pysecdec_pinched_other_lmb();
    test_integrate_2l_pysecdec_rank_four_num();
    test_integrate_3l_pysecdec();
    test_integrate_3l_rank_4();
    test_integrate_3l_rank_4_matad();
}
