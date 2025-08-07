use std::fmt::Display;

use color_eyre::Result;
use eyre::eyre;
use itertools::Itertools;
use linnet::half_edge::involution::{EdgeIndex, Orientation};
use linreg::linear_regression_of;
use rand::Rng;
use symbolica::{
    domains::float::NumericalFloatLike, domains::float::Real, numerical_integration::MonteCarloRng,
};

use crate::{
    cff::surface::Surface,
    model::Model,
    momentum::ThreeMomentum,
    momentum_sample::{LoopIndex, LoopMomenta, MomentumSample},
    new_gammaloop_integrand::{cross_section_integrand::CrossSectionGraphTerm, GraphTerm},
    new_graph::FeynmanGraph,
    utils::{box_muller, f128, FloatLike, F},
    DependentMomentaConstructor, Settings,
};

const SLOPE_STABILITY_points: usize = 50;

/// The range is from 10^start to 10^end.
pub struct ApproachSettings {
    pub lambda_exp_start: F<f64>,
    pub lambda_exp_end: F<f64>,
    pub steps: usize,
}

impl CrossSectionGraphTerm {
    pub(crate) fn test_ir(
        &self,
        settings: &Settings,
        approach_settings: &ApproachSettings,
    ) -> Result<()> {
        let rng = MonteCarloRng::new(settings.integrator.seed, 0);

        for (cut_id, esurface) in self.cut_esurface.iter_enumerated() {
            let supergraph_loop_count = self.graph.underlying.get_loop_number();
            let cut_cardinality = esurface.get_positive_energies().count();

            let massless_edges_in_cut = esurface
                .get_positive_energies()
                .filter(|edge_id| self.graph.underlying[*edge_id].0.particle.is_massless())
                .copied()
                .sorted()
                .collect_vec();

            if massless_edges_in_cut.is_empty() {
                println!("No IR limits for cut {}", cut_id);
                continue;
            }

            if massless_edges_in_cut.len() == cut_cardinality {
                println!("todo: implement loop over cmbs")
            } else {
                // find a loop momentum basis to easily do the ir limit in
                let (lmb_id, lmb) = self
                    .lmbs
                    .iter_enumerated()
                    .find(|(_, lmb)| {
                        lmb.loop_edges
                            .iter()
                            .all(|edge| massless_edges_in_cut.contains(edge))
                    })
                    .ok_or(eyre!(
                        "could not find cut momentum basis for cut: {}",
                        cut_id
                    ))?;

                let loop_indices = massless_edges_in_cut.iter().map(|edge| {
                    LoopIndex(
                        lmb.loop_edges
                            .iter()
                            .position(|loop_edge| loop_edge == edge)
                            .unwrap_or_else(|| unreachable!("corrupted lmb and cut")),
                    )
                });

                // enumeration of all ir limits
                todo!()
            }
        }

        Ok(())
    }

    pub(crate) fn test_single_ir_limit(
        &self,
        limit: &str,
        rng: &mut MonteCarloRng,
        settings: &Settings,
        approach_settings: &ApproachSettings,
        model: &Model,
    ) -> Result<()> {
        let ir_limit = IrLimit::parse_limit(limit)?;
        self.test_single_ir_limit_impl(&ir_limit, rng, settings, approach_settings, model)
    }

    fn test_single_ir_limit_impl(
        &self,
        ir_limit: &IrLimit,
        rng: &mut MonteCarloRng,
        settings: &Settings,
        approach_settings: &ApproachSettings,
        model: &Model,
    ) -> Result<()> {
        let edges_in_limit = ir_limit.get_all_edges()?;

        // find cut that for that as all edges of the limit
        let (cut_id, _esurface) = self
            .cut_esurface
            .iter_enumerated()
            .find(|(_cut_id, esurface)| {
                edges_in_limit
                    .iter()
                    .all(|edge| esurface.energies.contains(edge))
            })
            .ok_or(eyre!(
                "could not find cut with all edges of the limit: {}",
                ir_limit
            ))?;

        let cs_cut = &self.cuts[cut_id];

        let edges_to_flip = cs_cut
            .cut
            .iter_edges(&self.graph.underlying)
            .map(|(or, _)| or)
            .zip(
                self.graph
                    .underlying
                    .iter_edges_of(&cs_cut.cut)
                    .map(|x| x.1),
            )
            .filter_map(|(orientation, edge_id)| {
                if edges_in_limit.contains(&edge_id) && matches!(orientation, Orientation::Reversed)
                {
                    Some(edge_id)
                } else {
                    None
                }
            })
            .collect_vec();

        // find lmb
        let (_, lmb) = self
            .lmbs
            .iter_enumerated()
            .find(|(_, lmb)| {
                edges_in_limit
                    .iter()
                    .all(|edge| lmb.loop_edges.contains(edge))
            })
            .ok_or(eyre!(
                "could not find lmb to approach limit in: {}",
                ir_limit
            ))?;

        let model_parameter_cache = model.generate_values::<f128>();

        let momenta = ir_limit.get_momenta(rng, settings, approach_settings);
        let non_limit_loops = lmb
            .loop_edges
            .iter_enumerated()
            .filter_map(|(loop_id, edge_id)| {
                if !edges_in_limit.contains(edge_id) {
                    Some(loop_id)
                } else {
                    None
                }
            })
            .collect_vec();

        let non_limit_momenta = non_limit_loops
            .iter()
            .map(|loop_id| (*loop_id, sample_random_unit_vector(rng)))
            .collect_vec();

        let dependent_momenta_constructor = DependentMomentaConstructor::CrossSection {
            external_connections: todo!(), //;self.graph.external_connections.as_ref().unwrap(),
        };

        let loop_number = lmb.loop_edges.len();
        let external_particles = self.graph.underlying.get_external_partcles();
        let polarizations = settings
            .kinematics
            .externals
            .generate_polarizations(&external_particles, dependent_momenta_constructor);

        let limit_data = LimitData {
            data: momenta
                .into_iter()
                .map(|lambda_point| {
                    let mut loop_moms: LoopMomenta<F<_>> = (0..loop_number)
                        .map(|_| {
                            ThreeMomentum::new(F::from_f64(0.0), F::from_f64(0.0), F::from_f64(0.0))
                        })
                        .collect();

                    for (loop_id, momentum) in non_limit_momenta.iter() {
                        loop_moms[*loop_id] = momentum.clone();
                    }

                    for tagged_momenta in &lambda_point.momenta {
                        let edge_id = tagged_momenta.tag;
                        let loop_id = lmb
                            .loop_edges
                            .iter()
                            .position(|loop_edge| loop_edge == &edge_id)
                            .unwrap_or_else(|| {
                                unreachable!("corrupted lmb and ir limit: {}", ir_limit);
                            });

                        loop_moms[LoopIndex(loop_id)] = tagged_momenta.momentum.clone();
                    }

                    let sample = MomentumSample::new(
                        loop_moms,
                        &settings.kinematics.externals,
                        F::from_f64(1.0),
                        &polarizations,
                        dependent_momenta_constructor,
                        None,
                    );

                    LambdaPointEval {
                        lambda_point,
                        value: self
                            .evaluate(&sample, settings, todo!())
                            .norm_squared()
                            .sqrt(),
                    }
                })
                .collect(),
        };

        let slope = limit_data.extract_power();

        Ok(println!("slope: {}", slope))
    }
}

struct IrLimit {
    colinear: Vec<Vec<HardOrSoft>>,
    soft: Vec<EdgeIndex>,
}

enum MomentumBuilder<T: FloatLike> {
    Colinear {
        edge_id: EdgeIndex,
        x: F<T>,
        colinear_direction: ThreeMomentum<F<T>>,
        perpendicular_direction: ThreeMomentum<F<T>>,
        is_soft: bool,
    },
    Soft {
        edge_id: EdgeIndex,
        direction: ThreeMomentum<F<T>>,
    },
}

impl Display for IrLimit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for colinear_set in &self.colinear {
            write!(f, "C[")?;

            for i in 0..(colinear_set.len() - 1) {
                write!(f, "{},", colinear_set[i])?;
            }

            write!(f, "{}]", colinear_set.last().unwrap())?;
        }

        for soft in self.soft.iter() {
            write!(f, "S({})", soft)?;
        }

        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum HardOrSoft {
    Hard(EdgeIndex),
    Soft(EdgeIndex),
}

impl HardOrSoft {
    fn index(&self) -> EdgeIndex {
        match self {
            HardOrSoft::Hard(index) => *index,
            HardOrSoft::Soft(index) => *index,
        }
    }
}

impl Display for HardOrSoft {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            HardOrSoft::Hard(index) => write!(f, "{}", index),
            HardOrSoft::Soft(index) => write!(f, "S({})", index),
        }
    }
}

impl IrLimit {
    fn check_min_colinear_size(&self) -> bool {
        self.colinear
            .iter()
            .all(|colinear_set| colinear_set.len() >= 2)
    }

    fn is_valid(&self, loop_number: usize) -> Result<()> {
        if !self.check_min_colinear_size() {
            return Err(eyre!("colinear sets must have at least two edges"));
        }

        let all_edges = self.get_all_edges()?;

        if all_edges.len() > loop_number {
            return Err(eyre!("not enough degrees of freedom to setup IR limit"));
        }

        Ok(())
    }

    fn get_all_edges(&self) -> Result<Vec<EdgeIndex>> {
        let colinear_edges = self
            .colinear
            .iter()
            .flatten()
            .map(HardOrSoft::index)
            .collect_vec();
        let soft_edges = self.soft.iter().copied().collect_vec();

        let all_edges: Vec<EdgeIndex> = colinear_edges
            .into_iter()
            .chain(soft_edges)
            .sorted()
            .collect();

        // check for duplicates
        let mut unique_edges = all_edges.clone();
        unique_edges.dedup();

        if unique_edges.len() != all_edges.len() {
            return Err(eyre!("Edges specified in ir limit must be unique")); // duplicates found
        }

        Ok(all_edges)
    }

    fn parse_limit(limit: &str) -> Result<IrLimit> {
        let mut colinear_sets = Vec::new();
        let mut soft_edges = Vec::new();

        let mut char_iter = limit.chars().enumerate();

        while let Some((char_position, char)) = char_iter.next() {
            match char {
                'C' => {
                    let (_opening_bracket_position, opening_bracket) =
                        char_iter.next().ok_or_else(|| {
                            eyre!(
                                "Expected opening bracket after 'C' at position {}",
                                char_position
                            )
                        })?;

                    if opening_bracket != '[' {
                        return Err(eyre!(
                            "Expected '[' after 'C' at position {}, found '{}'",
                            char_position,
                            opening_bracket
                        ));
                    }

                    let mut colinear_set_str = String::new();

                    let mut closing_bracket_found = false;
                    while let Some((_next_char_position, next_char)) = char_iter.next() {
                        if next_char == ']' {
                            closing_bracket_found = true;
                            break;
                        }
                        colinear_set_str.push(next_char);
                    }

                    if !closing_bracket_found {
                        return Err(eyre!(
                            "Expected closing bracket ']' for colinear set at position {}",
                            char_position
                        ));
                    }

                    let edges = colinear_set_str.trim().split(',');

                    let mut colinear_set = Vec::new();

                    for edge in edges {
                        let trimmed_edge = edge.trim();
                        if trimmed_edge.is_empty() {
                            return Err(eyre!(
                                "Empty edge found in colinear set at position {}",
                                char_position
                            ));
                        }

                        if trimmed_edge.starts_with('S') {
                            let mut trimmed_edge_iter = trimmed_edge.chars().skip(1);
                            let opening_bracket = trimmed_edge_iter
                                .next()
                                .ok_or(eyre!("Expected '(' after 'S' in soft edge at position"))?;

                            if opening_bracket != '(' {
                                return Err(eyre!(
                                    "Expected '(' after 'S' in soft edge at position {}, found '{}'",));
                            }

                            let mut edge_str = String::new();
                            let mut closing_bracket_found = false;

                            while let Some(next_char) = trimmed_edge_iter.next() {
                                if next_char == ')' {
                                    closing_bracket_found = true;
                                    break;
                                }
                                edge_str.push(next_char);
                            }

                            if !closing_bracket_found {
                                return Err(eyre!(
                                    "Expected closing bracket ')' for soft edge at position {}",
                                    char_position
                                ));
                            }

                            let edge_index = Self::parse_edge(&edge_str)?;
                            colinear_set.push(HardOrSoft::Soft(edge_index));
                        } else {
                            let edge_index = Self::parse_edge(trimmed_edge)?;
                            colinear_set.push(HardOrSoft::Hard(edge_index));
                        }
                    }
                    colinear_sets.push(colinear_set);
                }
                'S' => {
                    let (_opening_bracket_position, opening_bracket) =
                        char_iter.next().ok_or_else(|| {
                            eyre!(
                                "Expected opening bracket '(' after 'S' at position {}",
                                char_position
                            )
                        })?;

                    if opening_bracket != '(' {
                        return Err(eyre!(
                            "Expected '(' after 'S' at position {}, found '{}'",
                            char_position,
                            opening_bracket
                        ));
                    }

                    let mut edge_str = String::new();
                    let mut closing_bracket_found = false;

                    while let Some((_next_char_position, next_char)) = char_iter.next() {
                        if next_char == ')' {
                            closing_bracket_found = true;
                            break;
                        }
                        edge_str.push(next_char);
                    }

                    if !closing_bracket_found {
                        return Err(eyre!(
                            "Expected closing bracket ')' for soft edge at position {}",
                            char_position
                        ));
                    }

                    let edge_index = Self::parse_edge(&edge_str)?;
                    soft_edges.push(edge_index);
                }
                _ => {
                    return Err(eyre!(
                        "Unexpected character '{}' at position {}",
                        char,
                        char_position
                    ))
                }
            }
        }

        let ir_limit = IrLimit {
            colinear: colinear_sets,
            soft: soft_edges,
        };

        Ok(ir_limit)
    }

    fn parse_edge(edge: &str) -> Result<EdgeIndex> {
        let mut edge = String::from(edge);
        edge = String::from(edge.trim());

        if edge.len() < 2 {
            return Err(eyre!("Edge must be at least two characters long"));
        }

        if edge.remove(0) != 'e' {
            return Err(eyre!("Edge must start with 'e'"));
        }

        let edge_id: usize = edge
            .parse()
            .map_err(|_| eyre!("Edge must be a valid integer, got: {}", edge))?;

        Ok(EdgeIndex::from(edge_id))
    }

    fn get_momentum_builders(&self, rng: &mut MonteCarloRng) -> Vec<MomentumBuilder<f128>> {
        let mut momentum_builder = Vec::new();

        for colinear_set in &self.colinear {
            let direction_for_set: ThreeMomentum<F<f128>> = sample_random_unit_vector(rng);

            let x_variables: Vec<F<f128>> = (0..colinear_set.len())
                .map(|_| F::from_f64(rng.random::<f64>() * 0.8 + 0.1))
                .sorted_by(|a, b| a.partial_cmp(b).unwrap())
                .collect_vec();

            for (edge, x) in colinear_set.iter().zip(x_variables) {
                let edge_id = edge.index();
                let direction: ThreeMomentum<F<f128>> = sample_random_unit_vector(rng);

                let perpendicular = direction.clone()
                    - direction.clone() * (direction.clone() * direction_for_set.clone());

                let perpendicular_norm = perpendicular.norm();
                let perpendicular = perpendicular * perpendicular_norm.inv();

                let is_soft = matches!(edge, HardOrSoft::Soft(_));

                momentum_builder.push(MomentumBuilder::Colinear {
                    edge_id,
                    x,
                    colinear_direction: direction_for_set.clone(),
                    perpendicular_direction: perpendicular,
                    is_soft,
                });
            }
        }

        for soft_edge in &self.soft {
            let direction = sample_random_unit_vector(rng);
            momentum_builder.push(MomentumBuilder::Soft {
                edge_id: *soft_edge,
                direction,
            });
        }

        momentum_builder
    }

    fn get_momenta(
        &self,
        rng: &mut MonteCarloRng,
        settings: &Settings,
        approach_settings: &ApproachSettings,
    ) -> Vec<LambdaPoint<f128>> {
        let momentum_builders = self.get_momentum_builders(rng);

        let lambda_exp_delta = (approach_settings.lambda_exp_start
            - approach_settings.lambda_exp_end)
            / F(approach_settings.steps as f64);

        let lambda_values: Vec<F<f128>> =
            (0..=approach_settings.steps)
                .map(|i| {
                    F::from_f64(10.0f64.powf(
                        approach_settings.lambda_exp_start.0 - &lambda_exp_delta.0 * (i as f64),
                    ))
                })
                .collect();

        lambda_values
            .into_iter()
            .map(|lambda| LambdaPoint {
                momenta: momentum_builders
                    .iter()
                    .map(|builder| match builder {
                        MomentumBuilder::Colinear {
                            edge_id,
                            x,
                            colinear_direction,
                            perpendicular_direction,
                            is_soft,
                        } => {
                            let momentum = if *is_soft {
                                (colinear_direction * x + perpendicular_direction * &lambda)
                                    * F::from_ff64(settings.kinematics.e_cm)
                                    * &lambda
                            } else {
                                (colinear_direction * x + perpendicular_direction * &lambda)
                                    * F::from_ff64(settings.kinematics.e_cm)
                            };
                            TaggedMomenta {
                                momentum,
                                tag: *edge_id,
                            }
                        }
                        MomentumBuilder::Soft { edge_id, direction } => {
                            let momentum =
                                direction * &lambda * F::from_ff64(settings.kinematics.e_cm);
                            TaggedMomenta {
                                momentum,
                                tag: *edge_id,
                            }
                        }
                    })
                    .collect(),
                lambda,
            })
            .collect()
    }
}

fn sample_random_unit_vector<T: FloatLike>(rng: &mut MonteCarloRng) -> ThreeMomentum<F<T>> {
    let x_1 = F::from_f64(rng.random::<f64>());
    let x_2 = F::from_f64(rng.random::<f64>());
    let x_3 = F::from_f64(rng.random::<f64>());
    let x_4 = F::from_f64(rng.random::<f64>());

    let (k_x, k_y) = box_muller(x_1, x_2);
    let (k_z, _) = box_muller(x_3, x_4);

    let unnormalized_momentum = ThreeMomentum::new(k_x, k_y, k_z);
    let norm = unnormalized_momentum.norm();
    unnormalized_momentum * norm.inv()
}

struct TaggedMomenta<T> {
    momentum: ThreeMomentum<T>,
    tag: EdgeIndex,
}

struct LambdaPoint<T: FloatLike> {
    lambda: F<T>,
    momenta: Vec<TaggedMomenta<F<T>>>,
}

struct LambdaPointEval<T: FloatLike> {
    lambda_point: LambdaPoint<T>,
    value: F<T>,
}

struct LimitData<T: FloatLike> {
    data: Vec<LambdaPointEval<T>>,
}

impl<T: FloatLike> LimitData<T> {
    fn extract_power(&self) -> f64 {
        let log_lambda = self
            .data
            .iter()
            .map(|eval| eval.lambda_point.lambda.0.to_f64().ln());

        let log_res = self.data.iter().map(|eval| eval.value.0.to_f64().ln());

        let log_data = log_lambda.zip(log_res).collect_vec();

        let mut stable_slope = false;
        let mut current_slope = 0.0;
        let mut previous_slope: f64 = linear_regression_of(&log_data).unwrap().0;
        let mut current_index_start = SLOPE_STABILITY_points;

        while !stable_slope {
            current_slope = linear_regression_of(&log_data[current_index_start..])
                .unwrap()
                .0;

            if (current_slope - previous_slope).abs()
                / ((current_slope + previous_slope) / 2.0).abs()
                < 0.01
            {
                stable_slope = true;
            } else {
                previous_slope = current_slope;
                current_index_start += SLOPE_STABILITY_points;
            }
        }

        if !stable_slope {
            println!("Warning: accuracy of xi not guaranteed")
        }

        current_slope
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_display() {
        let ir_limit = IrLimit {
            colinear: vec![
                vec![
                    HardOrSoft::Hard(EdgeIndex::from(1)),
                    HardOrSoft::Hard(EdgeIndex::from(2)),
                    HardOrSoft::Hard(EdgeIndex::from(3)),
                ],
                vec![
                    HardOrSoft::Hard(EdgeIndex::from(4)),
                    HardOrSoft::Soft(EdgeIndex::from(5)),
                ],
            ],
            soft: vec![EdgeIndex::from(6), EdgeIndex::from(7)],
        };

        let display = ir_limit.to_string();
        let expected = "C[e1,e2,e3]C[e4,S(e5)]S(e6)S(e7)";

        assert_eq!(display, expected);
    }

    #[test]
    fn parse_edge() {
        let edge_str = "e5";
        let edge_index = IrLimit::parse_edge(edge_str).unwrap();
        assert_eq!(edge_index, EdgeIndex::from(5));

        let invalid_edge_str = "5"; // missing 'e'
        assert!(IrLimit::parse_edge(invalid_edge_str).is_err());

        let invalid_edge_str2 = "e"; // too short
        assert!(IrLimit::parse_edge(invalid_edge_str2).is_err());

        let invalid_edge_str3 = "e5a"; // not a valid integer
        assert!(IrLimit::parse_edge(invalid_edge_str3).is_err());
    }

    #[test]
    fn parse_limit() {
        let limit_str = "C[e1,e2,e3]C[e4,S(e5)]S(e6)S(e7)";
        let ir_limit = IrLimit::parse_limit(limit_str).unwrap();

        assert_eq!(ir_limit.colinear.len(), 2, "Expected two colinear sets");
        assert_eq!(ir_limit.soft.len(), 2, "Expected two soft edges");

        assert_eq!(
            ir_limit.colinear[0],
            vec![
                HardOrSoft::Hard(EdgeIndex::from(1)),
                HardOrSoft::Hard(EdgeIndex::from(2)),
                HardOrSoft::Hard(EdgeIndex::from(3))
            ],
            "First colinear set does not match"
        );
        assert_eq!(
            ir_limit.colinear[1],
            vec![
                HardOrSoft::Hard(EdgeIndex::from(4)),
                HardOrSoft::Soft(EdgeIndex::from(5))
            ],
            "Second colinear set does not match"
        );

        assert_eq!(
            ir_limit.soft,
            vec![EdgeIndex::from(6), EdgeIndex::from(7)],
            "Soft edges do not match"
        );

        let invalid_limit_str = "C[e1,e2,e3]C[e4,S(e5)]S(e6)S(e7, e8)";
        assert!(
            IrLimit::parse_limit(invalid_limit_str).is_err(),
            "Expected error"
        );

        let invalid_limit_str2 = "C[e1,e2,e3C[e4,S(e5)]S(e6)";
        assert!(
            IrLimit::parse_limit(invalid_limit_str2).is_err(),
            "Expected error for unmatched brackets"
        );
    }
}
