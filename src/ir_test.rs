use std::fmt::Display;

use bincode::enc::write;
use color_eyre::Result;
use eyre::eyre;
use itertools::Itertools;
use linnet::half_edge::involution::EdgeIndex;
use rand::Rng;
use symbolica::{domains::float::NumericalFloatLike, numerical_integration::MonteCarloRng};

use crate::{
    cff::surface::Surface,
    momentum::ThreeMomentum,
    momentum_sample::LoopIndex,
    new_gammaloop_integrand::cross_section_integrand::CrossSectionGraphTerm,
    new_graph::FeynmanGraph,
    utils::{box_muller, f128, FloatLike, F},
    Settings,
};

impl CrossSectionGraphTerm {
    pub fn test_ir(&self, settings: &Settings) -> Result<()> {
        let mut rng = MonteCarloRng::new(settings.integrator.seed, 0);

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
            }
        }

        Ok(())
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
    fn is_valid(&self, loop_number: usize) -> Result<()> {
        if !self
            .colinear
            .iter()
            .all(|colinear_set| colinear_set.len() >= 2)
        {
            return Err(eyre!("colinear sets must have at least two edges"));
        }

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

        if all_edges.len() > loop_number {
            return Err(eyre!("not enough degrees of freedom to setup IR limit"));
        }

        Ok(())
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

        Ok(IrLimit {
            colinear: colinear_sets,
            soft: soft_edges,
        })
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

    fn get_momentum_builders(&self, rng: &mut MonteCarloRng) {
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
