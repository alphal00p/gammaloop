use std::collections::BTreeMap;

use momtrop::float::MomTropFloat;
use symbolica::parse;

use super::{
    ThresholdMultiplierEvaluator, ThresholdMultiplierInput, ThresholdMultiplierInputValues,
    ThresholdMultiplierLayout, ThresholdMultiplierPoint,
};
use crate::{
    graph::{Graph, parse::from_dot::IntoGraph, threshold_counterterms::ThresholdCountertermSpec},
    initialisation::test_initialise,
    integrands::evaluation::EvaluationMetaData,
    processes::EvaluatorSettings,
    utils::{ArbPrec, F, FloatLike, load_generic_model},
};

fn compare_multipliers<T: FloatLike>(
    layout: &ThresholdMultiplierLayout,
    evaluators: &mut [(ThresholdMultiplierEvaluator, ThresholdMultiplierEvaluator)],
    zero: F<T>,
) {
    let mt = zero.from_usize(173);
    let tolerance = zero.epsilon() * zero.from_usize(256);
    // All vectors are in the same CM routing as GL638. The last point deliberately lies
    // off the host shell and exercises the WH common-zero convention, not physical WF.
    for (label, soft_power) in [
        ("host", None),
        ("off_host", None),
        ("soft_10", Some(10)),
        ("soft_30", Some(30)),
        ("soft_45", Some(45)),
        ("WH_common_zero", None),
    ] {
        let common_zero = label == "WH_common_zero";
        let a = if common_zero {
            [0, 0, 0]
        } else {
            [16, -28, 21]
        }
        .map(|x| zero.from_isize(x));
        let t = if common_zero {
            [0, 0, 0]
        } else {
            [-31, 57, 91]
        }
        .map(|x| zero.from_isize(x));
        let p = if let Some(power) = soft_power {
            let distance = zero.from_usize(2).powi(-power);
            std::array::from_fn(|i| &t[i] + &distance * zero.from_isize([1, 2, -1][i]))
        } else {
            if common_zero {
                [0, 0, 0]
            } else {
                [85, -42, 130]
            }
            .map(|x| zero.from_isize(x))
        };
        let mut momenta = vec![[zero.clone(), zero.clone(), zero.clone()]; 17];
        momenta[2] = a.clone();
        momenta[3] = p.clone();
        momenta[4] = std::array::from_fn(|i| &a[i] + &p[i]);
        momenta[5] = t.clone();
        momenta[6] = std::array::from_fn(|i| &a[i] + &t[i]);
        momenta[7] = [101, 32, -9].map(|x| zero.from_isize(x));
        momenta[8] = momenta[7].clone();
        momenta[10] = t.clone();
        momenta[12] = p.clone();
        momenta[13] = std::array::from_fn(|i| &p[i] - &t[i]);
        momenta[14] = std::array::from_fn(|i| &momenta[7][i] - &p[i]);
        let energy = |edge: usize, factor: &F<T>| {
            let mass = match edge {
                2 => zero.from_usize(125),
                13 | 14 => zero.clone(),
                _ => mt.clone(),
            };
            (momenta[edge].iter().fold(mass.square(), |sum, component| {
                sum + (component * factor).square()
            }))
            .sqrt()
        };
        let q = match label {
            "off_host" => zero.from_usize(600),
            "WH_common_zero" => &mt * zero.from_usize(2),
            _ => energy(2, &zero.one()) + energy(6, &zero.one()) + energy(10, &zero.one()),
        };
        let mut values = ThresholdMultiplierInputValues::new(layout, zero.clone());
        for (index, input) in layout.inputs().iter().enumerate() {
            // Distinguish the two frames so an accidental effective/star alias is detected.
            let frame_factor = |point: &ThresholdMultiplierPoint| match point {
                ThresholdMultiplierPoint::Star => zero.one(),
                ThresholdMultiplierPoint::Effective => zero.from_usize(2),
            };
            let value = match input {
                ThresholdMultiplierInput::ModelParameter { index: 0 } => mt.clone(),
                ThresholdMultiplierInput::ExternalMomentum { component: 0, .. } => {
                    &q / zero.from_usize(2)
                }
                ThresholdMultiplierInput::EdgeMomentum {
                    point,
                    edge,
                    component,
                } => {
                    if *component == 0 {
                        zero.clone()
                    } else {
                        &momenta[*edge][*component - 1] * frame_factor(point)
                    }
                }
                ThresholdMultiplierInput::EdgeEnergy { point, edge } => {
                    energy(*edge, &frame_factor(point))
                }
                _ => zero.clone(),
            };
            values.set_real(index, value).unwrap();
        }
        for (index, (old, compact)) in evaluators.iter_mut().enumerate() {
            if common_zero && index >= 4 {
                continue;
            }
            let old = old
                .evaluate(&values, &mut EvaluationMetaData::new_empty())
                .unwrap();
            let compact = compact
                .evaluate(&values, &mut EvaluationMetaData::new_empty())
                .unwrap();
            let scale = old.abs().max(compact.abs());
            assert!(
                (&old - &compact).abs() <= &tolerance * &scale,
                "{label}, multiplier {index}: expanded={old}, compact={compact}"
            );
            if common_zero {
                assert_eq!(compact, zero.from_usize(usize::from(index % 2 == 0)));
            }
        }
    }
}

#[test]
fn gl638_compact_function_map_matches_expanded_multipliers() {
    test_initialise().unwrap();
    let graph: Graph = include_str!("../../../../../../tests/resources/graphs/GL638.dot")
        .into_graph(&load_generic_model("sm"))
        .unwrap();
    let original = ThresholdCountertermSpec::parse_toml(include_str!(
        "../../../../../../tests/resources/graphs/ir_safe_thresholds/GL638_legacy_multipliers.toml"
    ))
    .unwrap();
    let layout = ThresholdMultiplierLayout::new(
        vec![parse!("UFO::MT")],
        Vec::new(),
        2,
        (0..17).collect(),
        Vec::new(),
    )
    .unwrap();
    let mut evaluators = Vec::new();
    for (old_cut, cut) in original.cuts.iter().zip(&graph.threshold_counterterms.cuts) {
        assert_eq!(old_cut.edges, cut.edges);
        for (old_threshold, threshold) in old_cut.thresholds.iter().zip(&cut.thresholds) {
            assert_eq!(old_threshold.edges, threshold.edges);
            for (old, current) in old_threshold
                .counterterms
                .iter()
                .zip(&threshold.counterterms)
            {
                assert_eq!(old.subspace, current.subspace);
                assert_eq!(old.parent_lmb, current.parent_lmb);
                if let Some(old) = &old.multiplier {
                    let current = current.multiplier.as_ref().unwrap();
                    let old = layout
                        .parse_expression(&old.expression, &BTreeMap::new())
                        .unwrap();
                    let current = layout
                        .parse_expression(
                            &current.expression,
                            &graph.threshold_counterterms.function_map,
                        )
                        .unwrap();
                    evaluators.push((
                        layout
                            .build_evaluator(&old, &EvaluatorSettings::default())
                            .unwrap(),
                        layout
                            .build_evaluator(&current, &EvaluatorSettings::default())
                            .unwrap(),
                    ));
                }
            }
        }
    }
    assert_eq!(evaluators.len(), 6);
    compare_multipliers(&layout, &mut evaluators, F(0.0));
    compare_multipliers(&layout, &mut evaluators, F::<ArbPrec>::default());
}
