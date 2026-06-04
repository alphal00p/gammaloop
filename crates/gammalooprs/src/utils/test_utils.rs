use std::{env, path::PathBuf};

use insta::assert_snapshot;
use itertools::Itertools;
use linnet::half_edge::{
    HedgeGraph, NoData,
    builder::HedgeGraphBuilder,
    involution::{EdgeIndex, Orientation},
};
use momtrop::assert_approx_eq;
use spenso::structure::abstract_index::AIND_SYMBOLS;
use symbolica::{atom::AtomCore, parse_lit};

use crate::{
    momentum::ThreeMomentum,
    settings::runtime::{ParameterizationMapping, ParameterizationMode, ParameterizationSettings},
    utils::{F, GS, global_inv_parameterize, global_parameterize},
};

pub fn output_dir() -> PathBuf {
    if let Ok(pytest_output_path) = env::var("PYTEST_OUTPUT_PATH_FOR_RUST") {
        pytest_output_path.into()
    } else {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("src/test_resources")
    }
}

pub(crate) fn dummy_hedge_graph(n_edges: usize) -> HedgeGraph<NoData, NoData, NoData> {
    let mut builder = HedgeGraphBuilder::<NoData, NoData>::new();
    let left = builder.add_node(NoData {});
    let right = builder.add_node(NoData {});
    for _ in 0..n_edges {
        builder.add_edge(left, right, NoData {}, Orientation::Default);
    }
    builder.build()
}

#[test]
fn normalization() {
    let a = GS.emr_vec_index(EdgeIndex(1), AIND_SYMBOLS.cind.call_args([0]));
    let b = GS.emr_vec_index(EdgeIndex(1), AIND_SYMBOLS.cind.call_args([1]));

    assert_snapshot!(a.to_canonical_string(),@"0");
    assert_snapshot!(b.to_canonical_string(),@"gammalooprs::{spenso::rank1,spenso::tensor}::Q(1,spenso::{}::cind(1))");

    let c = GS.energy_delta(GS.cind(1));
    assert_snapshot!(c.to_canonical_string(),@"0");
    let c = GS.energy_delta(GS.cind(0));
    assert_snapshot!(c.to_canonical_string(),@"1");

    let expr = parse_lit!(f(a + p + r));
    assert_snapshot!(GS.linearize.call_args([expr]).to_canonical_string(),@"gammalooprs::{}::f(gammalooprs::{}::a)+gammalooprs::{}::f(gammalooprs::{}::p)+gammalooprs::{}::f(gammalooprs::{}::r)");
}

#[test]
fn test_inv_param() {
    let x = [
        F(0.1),
        F(0.2),
        F(0.3),
        F(0.4),
        F(0.5),
        F(0.6),
        F(0.7),
        F(0.8),
        F(0.9),
    ];
    let e_cm = F(42.2);
    let param_settings = ParameterizationSettings::default();

    let (momenta, jac_1) = global_parameterize(&x, e_cm, &param_settings);
    let actual_momenta = momenta
        .iter()
        .map(|p| ThreeMomentum {
            px: p[0],
            py: p[1],
            pz: p[2],
        })
        .collect_vec();

    let (xs, jac_2) = global_inv_parameterize(&actual_momenta, e_cm, &param_settings);

    for i in 0..9 {
        assert_approx_eq(&x[i], &xs[i], &F(1.0e-14));
    }

    let prod = jac_1 * jac_2;
    assert_approx_eq(&prod, &F(1.0), &F(1.0e-14));
}

#[test]
fn test_relative_spherical_inv_param() {
    let x = [F(0.17), F(0.23), F(0.31), F(0.43), F(0.59), F(0.71)];
    let e_cm = F(173.0);
    let param_settings = ParameterizationSettings {
        mode: ParameterizationMode::RelativeSpherical,
        ..Default::default()
    };

    let (momenta, jac_1) = global_parameterize(&x, e_cm, &param_settings);
    let actual_momenta = momenta
        .iter()
        .map(|p| ThreeMomentum {
            px: p[0],
            py: p[1],
            pz: p[2],
        })
        .collect_vec();

    let (xs, jac_2) = global_inv_parameterize(&actual_momenta, e_cm, &param_settings);

    for i in 0..6 {
        assert_approx_eq(&x[i], &xs[i], &F(1.0e-14));
    }

    let prod = jac_1 * jac_2;
    assert_approx_eq(&prod, &F(1.0), &F(1.0e-14));
}

#[test]
fn test_spherical_common_radial_inv_param() {
    let x = [F(0.17), F(0.23), F(0.31), F(0.43), F(0.59), F(0.71)];
    let e_cm = F(173.0);
    let param_settings = ParameterizationSettings {
        mode: ParameterizationMode::SphericalCommonRadial,
        ..Default::default()
    };

    let (momenta, jac_1) = global_parameterize(&x, e_cm, &param_settings);
    let actual_momenta = momenta
        .iter()
        .map(|p| ThreeMomentum {
            px: p[0],
            py: p[1],
            pz: p[2],
        })
        .collect_vec();

    let (xs, jac_2) = global_inv_parameterize(&actual_momenta, e_cm, &param_settings);

    for i in 0..6 {
        assert_approx_eq(&x[i], &xs[i], &F(1.0e-14));
    }

    let prod = jac_1 * jac_2;
    assert_approx_eq(&prod, &F(1.0), &F(1.0e-14));
}

#[test]
fn test_spherical_common_radial_inv_param_general_loop_counts() {
    let e_cm = F(173.0);
    let xs_by_loop_count = [
        vec![F(0.17), F(0.31), F(0.43)],
        vec![F(0.17), F(0.23), F(0.31), F(0.43), F(0.59), F(0.71)],
        vec![
            F(0.17),
            F(0.23),
            F(0.41),
            F(0.31),
            F(0.43),
            F(0.59),
            F(0.71),
            F(0.37),
            F(0.83),
        ],
        vec![
            F(0.17),
            F(0.23),
            F(0.41),
            F(0.29),
            F(0.31),
            F(0.43),
            F(0.59),
            F(0.71),
            F(0.37),
            F(0.83),
            F(0.19),
            F(0.67),
        ],
    ];
    let mappings = [
        ParameterizationMapping::Linear,
        ParameterizationMapping::Log,
        ParameterizationMapping::Power,
    ];

    for x in xs_by_loop_count {
        for mapping in &mappings {
            let param_settings = ParameterizationSettings {
                mode: ParameterizationMode::SphericalCommonRadial,
                mapping: mapping.clone(),
                power: 1.7,
                ..Default::default()
            };

            let (momenta, jac_1) = global_parameterize(&x, e_cm, &param_settings);
            let actual_momenta = momenta
                .iter()
                .map(|p| ThreeMomentum {
                    px: p[0],
                    py: p[1],
                    pz: p[2],
                })
                .collect_vec();

            let (xs, jac_2) = global_inv_parameterize(&actual_momenta, e_cm, &param_settings);

            for (expected, actual) in x.iter().zip(xs.iter()) {
                assert_approx_eq(expected, actual, &F(1.0e-13));
            }

            let prod = jac_1 * jac_2;
            assert_approx_eq(&prod, &F(1.0), &F(1.0e-12));
        }
    }
}

#[test]
fn test_spherical_common_radial_uses_lmb_ordered_stick_breaking() {
    let x = [
        F(0.17),
        F(0.23),
        F(0.41),
        F(0.31),
        F(0.43),
        F(0.59),
        F(0.71),
        F(0.37),
        F(0.83),
    ];
    let e_cm = F(173.0);
    let param_settings = ParameterizationSettings {
        mode: ParameterizationMode::SphericalCommonRadial,
        ..Default::default()
    };

    let (momenta, _) = global_parameterize(&x, e_cm, &param_settings);
    let radii = momenta
        .iter()
        .map(|p| (p[0].square() + p[1].square() + p[2].square()).sqrt())
        .collect_vec();
    let radius = radii.iter().fold(F(0.0), |acc, radius_i| acc + radius_i);

    let actual_first_fraction = radii[0] / radius;
    let actual_second_fraction = radii[1] / radius;
    let actual_third_fraction = radii[2] / radius;
    let expected_first_fraction = x[1];
    let expected_second_fraction = (F(1.0) - x[1]) * x[2];
    let expected_third_fraction = (F(1.0) - x[1]) * (F(1.0) - x[2]);

    assert_approx_eq(
        &actual_first_fraction,
        &expected_first_fraction,
        &F(1.0e-14),
    );
    assert_approx_eq(
        &actual_second_fraction,
        &expected_second_fraction,
        &F(1.0e-14),
    );
    assert_approx_eq(
        &actual_third_fraction,
        &expected_third_fraction,
        &F(1.0e-14),
    );
}
