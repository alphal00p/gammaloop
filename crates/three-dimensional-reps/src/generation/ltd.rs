//! Signed loop-energy residues with factorized, derivative-free numerator maps.
//!
//! Only repeated poles need a polynomial capacity: exact finite-difference
//! stencils replace their numerator derivatives. Every stencil shifts a common
//! set of loop energies, so all original edge arguments remain physically
//! consistent, including oppositely routed copies of one denominator.

use super::{
    CffEnergyFactorComponent, CffEnergyFactorOwnership, CffGlobalPrefactorSign,
    ExpressionAssembler, Generate3DExpressionOptions, GeneratedThreeDExpression, GenerationError,
    LinearEnergyExpr, MomentumSignature, ParsedGraph, Rational, RationalExt, RepresentationMode,
    Result, apply_initial_state_cut_edge_energy_exprs, denominator_tree_from_chain,
    edge_q0_from_loop_exprs, energy_residues, normalize_energy_degree_bounds, repeated_groups,
    solve_loop_energy_from_target_edge_exprs, solve_rational_system,
};
use crate::{
    ThreeDExpression,
    cut_structure::{ContourClosure, EnergyResidue},
    expression::CFFVariant,
    surface::{HybridSurfaceID, SurfaceOrigin},
    utils::determinant_i32,
};
use itertools::Itertools;
use linnet::half_edge::involution::EdgeIndex;
use std::collections::{BTreeMap, HashMap, hash_map::Entry};

struct Channel {
    members: Vec<usize>,
}

impl Channel {
    fn representative(&self) -> usize {
        self.members[0]
    }

    fn power(&self) -> usize {
        self.members.len()
    }
}

pub(super) struct LtdBuilder<'a> {
    parsed: &'a ParsedGraph,
    options: &'a Generate3DExpressionOptions,
    signatures: Vec<MomentumSignature>,
    channels: Vec<Channel>,
    bounds: Vec<usize>,
    stencils: HashMap<(usize, usize), Vec<(i64, Rational)>>,
    assembly: ExpressionAssembler,
}

impl<'a> LtdBuilder<'a> {
    pub(super) fn new(
        parsed: &'a ParsedGraph,
        options: &'a Generate3DExpressionOptions,
    ) -> Result<Self> {
        let repeated = repeated_groups(parsed);
        let bounds = if repeated.is_empty() {
            // Ordinary residues evaluate the complete numerator on their pole
            // maps; their validity has no dependence on an EMR degree bound.
            vec![0; parsed.internal_edges.len()]
        } else {
            normalize_energy_degree_bounds(
                options
                    .energy_degree_bounds
                    .as_deref()
                    .ok_or(GenerationError::LtdRepeatedPropagatorsRequireEnergyBounds)?,
                parsed.internal_edges.len(),
            )?
        };
        let repeated_by_member = repeated
            .iter()
            .flat_map(|group| {
                group
                    .edge_ids
                    .iter()
                    .map(move |edge| (*edge, &group.edge_ids))
            })
            .collect::<BTreeMap<_, _>>();
        let channels = parsed
            .denominator_internal_edge_ids()
            .into_iter()
            .filter_map(|edge| match repeated_by_member.get(&edge) {
                Some(members) if members[0] != edge => None,
                Some(members) => Some(Channel {
                    members: (*members).clone(),
                }),
                None => Some(Channel {
                    members: vec![edge],
                }),
            })
            .collect();
        Ok(Self {
            parsed,
            options,
            signatures: parsed
                .internal_edges
                .iter()
                .map(|edge| edge.signature.clone())
                .collect(),
            channels,
            bounds,
            stencils: HashMap::new(),
            assembly: ExpressionAssembler::new(ThreeDExpression::new_empty()),
        })
    }

    pub(super) fn build(mut self) -> Result<GeneratedThreeDExpression> {
        let loop_count = self.parsed.loop_names.len();
        let signatures = self
            .channels
            .iter()
            .map(|channel| {
                self.signatures[channel.representative()]
                    .loop_signature
                    .clone()
            })
            .collect::<Vec<_>>();
        let residues = if loop_count == 0 {
            vec![EnergyResidue {
                basis: Vec::new(),
                sigmas: Vec::new(),
                sign: 1,
                spanning_tree_index: 0,
            }]
        } else {
            if signatures.is_empty() {
                return Err(GenerationError::SingularBasis);
            }
            energy_residues(&signatures, &vec![ContourClosure::Below; loop_count])?
        };
        if residues.is_empty() {
            return Err(GenerationError::SingularBasis);
        }
        for residue in residues {
            self.add_residue(&residue)?;
        }
        let mut expression = self.assembly.expression.fuse_compatible_variants();
        super::assign_numerator_map_labels(&mut expression.orientations);
        let denominator_edges = self.parsed.denominator_internal_edge_ids();
        Ok(GeneratedThreeDExpression {
            representation: RepresentationMode::Ltd,
            expression,
            energy_factor_ownership: CffEnergyFactorOwnership::VariantLocal,
            energy_factor_components: (!denominator_edges.is_empty())
                .then_some(CffEnergyFactorComponent {
                    internal_edge_ids: denominator_edges,
                    ownership: CffEnergyFactorOwnership::VariantLocal,
                    denominator_only_global_prefactor_sign: CffGlobalPrefactorSign::default(),
                    core_global_prefactor_sign: CffGlobalPrefactorSign::default(),
                })
                .into_iter()
                .collect(),
            source_energy_degree_bounds: Vec::new(),
            // LTD is already the signed dq0/(2*pi*i) contour. Consumers must
            // not apply the CFF source-frame conversion to this representation.
            denominator_only_global_prefactor_sign: CffGlobalPrefactorSign::default(),
            core_global_prefactor_sign: CffGlobalPrefactorSign::default(),
        })
    }

    fn add_residue(&mut self, residue: &EnergyResidue) -> Result<()> {
        let basis = residue
            .basis
            .iter()
            .map(|channel| self.channels[*channel].representative())
            .collect::<Vec<_>>();
        let alpha = residue
            .basis
            .iter()
            .map(|channel| self.channels[*channel].power() - 1)
            .collect::<Vec<_>>();
        let mut targets = vec![LinearEnergyExpr::zero(); self.signatures.len()];
        for (edge, sigma) in basis.iter().zip(&residue.sigmas) {
            targets[*edge] = LinearEnergyExpr::ose(EdgeIndex(*edge), i64::from(*sigma));
        }
        let loop_energies =
            solve_loop_energy_from_target_edge_exprs(&self.signatures, &basis, &targets)?;
        let mut edge_energies = edge_q0_from_loop_exprs(&self.signatures, &loop_energies);
        apply_initial_state_cut_edge_energy_exprs(self.parsed, &mut edge_energies);

        // Differentiate only the scalar affine routing, never the graph
        // numerator. The same inverse defines each physical stencil map.
        let mut edge_derivatives = vec![vec![Rational::zero(); basis.len()]; self.signatures.len()];
        if alpha.iter().any(|order| *order != 0) {
            let zero_shift_signatures = self
                .signatures
                .iter()
                .cloned()
                .map(|mut signature| {
                    signature.external_signature.fill(0);
                    signature
                })
                .collect::<Vec<_>>();
            for (axis, edge) in basis.iter().enumerate() {
                let mut unit = vec![LinearEnergyExpr::zero(); self.signatures.len()];
                unit[*edge].constant = Rational::one();
                let direction = solve_loop_energy_from_target_edge_exprs(
                    &zero_shift_signatures,
                    &basis,
                    &unit,
                )?;
                for (row, derivative) in edge_derivatives
                    .iter_mut()
                    .zip(edge_q0_from_loop_exprs(&zero_shift_signatures, &direction))
                {
                    row[axis] = derivative.constant;
                }
            }
        }
        let degrees = (0..basis.len())
            .map(|axis| {
                edge_derivatives
                    .iter()
                    .zip(&self.bounds)
                    .filter_map(|(row, degree)| (!row[axis].is_zero()).then_some(*degree))
                    .try_fold(0usize, |sum, degree| sum.checked_add(degree))
                    .ok_or(GenerationError::CoefficientOutOfRange)
            })
            .collect::<Result<Vec<_>>>()?;
        let determinant = determinant_i32(
            &basis
                .iter()
                .map(|edge| self.signatures[*edge].loop_signature.clone())
                .collect::<Vec<_>>(),
        );
        let mut contour = Rational::from(residue.sign) / determinant.abs();
        for (sigma, power) in residue.sigmas.iter().zip(&alpha) {
            contour *= Rational::from(*sigma).pow_usize(*power);
        }
        let mut factors = Vec::new();
        for (channel_index, channel) in self.channels.iter().enumerate() {
            let representative = channel.representative();
            if let Some(axis) = residue
                .basis
                .iter()
                .position(|index| *index == channel_index)
            {
                let mut derivatives = vec![Rational::zero(); basis.len()];
                derivatives[axis] = Rational::from(residue.sigmas[axis]);
                factors.push(DenominatorFactor {
                    surface: FactorSurface::Cut(representative),
                    edges: channel.members.clone(),
                    power: channel.power(),
                    derivatives,
                });
            } else {
                for sign in [-1, 1] {
                    let surface = self.assembly.intern_surface(
                        edge_energies[representative].clone()
                            + LinearEnergyExpr::ose(EdgeIndex(representative), sign),
                        SurfaceOrigin::Physical,
                        false,
                    );
                    factors.push(DenominatorFactor {
                        surface: FactorSurface::Linear(surface),
                        edges: channel.members.clone(),
                        power: channel.power(),
                        derivatives: edge_derivatives[representative].clone(),
                    });
                }
            }
        }

        for beta in alpha
            .iter()
            .map(|order| 0..=*order)
            .multi_cartesian_product()
        {
            if beta
                .iter()
                .zip(&degrees)
                .any(|(order, degree)| order > degree)
            {
                continue;
            }
            let gamma = alpha
                .iter()
                .zip(&beta)
                .map(|(a, b)| a - b)
                .collect::<Vec<_>>();
            // binomial(alpha,beta)/alpha! = 1/(beta! gamma!). The
            // denominator expansion below supplies its own gamma!.
            let coefficient = beta
                .iter()
                .chain(&gamma)
                .fold(contour.clone(), |value, order| {
                    (1..=*order).fold(value, |value, factor| value / Rational::from(factor))
                });
            let denominator_terms = DenominatorFactor::differentiate_product(&factors, &gamma);
            let mut axes = Vec::new();
            for (axis, order) in beta
                .iter()
                .copied()
                .enumerate()
                .filter(|(_, order)| *order > 0)
            {
                let degree = degrees[axis];
                let stencil = match self.stencils.entry((degree, order)) {
                    Entry::Occupied(entry) => entry.into_mut(),
                    Entry::Vacant(entry) => entry.insert(Self::derivative_stencil(degree, order)?),
                };
                axes.push((
                    axis,
                    self.options
                        .numerator_sampling_scale
                        .is_active_for_degree(degree),
                    stencil.clone(),
                ));
            }
            for samples in axes
                .iter()
                .map(|(_, _, samples)| samples)
                .multi_cartesian_product()
            {
                let mut sample_targets = targets.clone();
                let mut sample_coefficient = coefficient.clone();
                let mut extra_half_edges = Vec::new();
                let mut uniform_scale_power = 0;
                for ((axis, uniform, _), (offset, weight)) in axes.iter().zip(samples) {
                    let edge = basis[*axis];
                    sample_coefficient *= weight;
                    sample_targets[edge] = sample_targets[edge].clone()
                        + if *uniform {
                            uniform_scale_power += beta[*axis];
                            LinearEnergyExpr::uniform_scale(*offset)
                        } else {
                            extra_half_edges
                                .extend(std::iter::repeat_n(EdgeIndex(edge), beta[*axis]));
                            sample_coefficient *= Rational::from(2).pow_usize(beta[*axis]);
                            LinearEnergyExpr::ose(EdgeIndex(edge), *offset)
                        };
                }
                let loop_map = solve_loop_energy_from_target_edge_exprs(
                    &self.signatures,
                    &basis,
                    &sample_targets,
                )?;
                let mut edge_map = edge_q0_from_loop_exprs(&self.signatures, &loop_map);
                apply_initial_state_cut_edge_energy_exprs(self.parsed, &mut edge_map);
                for term in &denominator_terms {
                    let mut half_edges = term.half_edges.clone();
                    half_edges.extend(extra_half_edges.iter().copied());
                    self.assembly.push_variant_for_maps(
                        loop_map.clone(),
                        edge_map.clone(),
                        CFFVariant {
                            origin: Some("ltd".to_string()),
                            prefactor: &sample_coefficient * &term.coefficient,
                            half_edges,
                            denominator_edges: term.edges.clone(),
                            denominator_surface_signs: BTreeMap::new(),
                            denominator_edge_support_signs: BTreeMap::new(),
                            uniform_scale_power,
                            numerator_surfaces: Vec::new(),
                            denominator: denominator_tree_from_chain(&term.chain),
                        },
                    );
                }
            }
        }
        Ok(())
    }

    fn derivative_stencil(degree: usize, order: usize) -> Result<Vec<(i64, Rational)>> {
        let nodes = (0..=degree)
            .map(|index| {
                let magnitude = i64::try_from(index.div_ceil(2))
                    .map_err(|_| GenerationError::CoefficientOutOfRange)?;
                Ok(if index % 2 == 0 {
                    -magnitude
                } else {
                    magnitude
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let matrix = (0..=degree)
            .map(|power| {
                nodes
                    .iter()
                    .map(|node| Rational::from(*node).pow_usize(power))
                    .collect()
            })
            .collect();
        let rhs = (0..=degree)
            .map(|power| {
                if power == order {
                    (1..=order).fold(Rational::one(), |value, factor| {
                        value * Rational::from(factor)
                    })
                } else {
                    Rational::zero()
                }
            })
            .collect();
        let weights = solve_rational_system(matrix, rhs).ok_or(GenerationError::SingularBasis)?;
        Ok(nodes
            .into_iter()
            .zip(weights)
            .filter(|(_, weight)| !weight.is_zero())
            .collect())
    }
}

#[derive(Clone, Copy)]
enum FactorSurface {
    Cut(usize),
    Linear(HybridSurfaceID),
}

struct DenominatorFactor {
    surface: FactorSurface,
    edges: Vec<usize>,
    power: usize,
    derivatives: Vec<Rational>,
}

#[derive(Clone)]
struct DenominatorTerm {
    coefficient: Rational,
    half_edges: Vec<EdgeIndex>,
    edges: Vec<EdgeIndex>,
    chain: Vec<HybridSurfaceID>,
}

impl DenominatorFactor {
    fn differentiate_product(factors: &[Self], gamma: &[usize]) -> Vec<DenominatorTerm> {
        let factorial = gamma.iter().fold(Rational::one(), |value, order| {
            (1..=*order).fold(value, |value, factor| value * Rational::from(factor))
        });
        let mut terms = vec![(
            gamma.to_vec(),
            DenominatorTerm {
                coefficient: factorial,
                half_edges: Vec::new(),
                edges: Vec::new(),
                chain: Vec::new(),
            },
        )];
        for factor in factors {
            terms = terms
                .into_iter()
                .flat_map(|(remaining, term)| {
                    remaining
                        .iter()
                        .zip(&factor.derivatives)
                        .map(|(order, derivative)| {
                            0..=if derivative.is_zero() { 0 } else { *order }
                        })
                        .multi_cartesian_product()
                        .map(move |delta| {
                            let order = delta.iter().sum::<usize>();
                            let mut derived = term.clone();
                            for index in 0..order {
                                derived.coefficient *= -Rational::from(factor.power + index);
                            }
                            for (order, derivative) in delta.iter().zip(&factor.derivatives) {
                                derived.coefficient *= derivative.pow_usize(*order);
                                for divisor in 1..=*order {
                                    derived.coefficient /= Rational::from(divisor);
                                }
                            }
                            let power = factor.power + order;
                            derived
                                .edges
                                .extend(factor.edges.iter().copied().map(EdgeIndex));
                            match factor.surface {
                                FactorSurface::Cut(edge) => derived
                                    .half_edges
                                    .extend(std::iter::repeat_n(EdgeIndex(edge), power)),
                                FactorSurface::Linear(surface) => {
                                    derived.chain.extend(std::iter::repeat_n(surface, power))
                                }
                            }
                            let remaining = remaining
                                .iter()
                                .zip(delta)
                                .map(|(left, right)| left - right)
                                .collect::<Vec<_>>();
                            (remaining, derived)
                        })
                        .collect::<Vec<_>>()
                })
                .collect();
        }
        terms
            .into_iter()
            .filter_map(|(remaining, mut term)| {
                if remaining.iter().any(|order| *order != 0) {
                    return None;
                }
                term.edges.sort_unstable();
                term.edges.dedup();
                Some(term)
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(feature = "eval")]
    use crate::NumeratorSamplingScaleMode;
    use crate::{
        generation::generate_3d_expression,
        graph_io::{ParsedGraphInternalEdge, test_graphs},
    };

    fn repeated_cycle(power: usize, signs: &[i32]) -> ParsedGraph {
        ParsedGraph {
            internal_edges: (0..power)
                .map(|edge_id| {
                    let (tail, head) = if signs[edge_id] > 0 {
                        (edge_id, (edge_id + 1) % power)
                    } else {
                        ((edge_id + 1) % power, edge_id)
                    };
                    ParsedGraphInternalEdge {
                        edge_id,
                        tail,
                        head,
                        label: format!("q{edge_id}"),
                        mass_key: Some("m".to_string()),
                        signature: MomentumSignature {
                            loop_signature: vec![signs[edge_id]],
                            external_signature: Vec::new(),
                        },
                        had_pow: false,
                    }
                })
                .collect(),
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["k".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: (0..power).map(|node| (format!("v{node}"), node)).collect(),
        }
    }

    #[test]
    fn repeated_ltd_requires_explicit_numerator_capacity() {
        let parsed = repeated_cycle(2, &[1, -1]);
        let options = Generate3DExpressionOptions {
            representation: RepresentationMode::Ltd,
            ..Default::default()
        };
        assert!(matches!(
            generate_3d_expression(&parsed, &options),
            Err(GenerationError::LtdRepeatedPropagatorsRequireEnergyBounds)
        ));
        assert!(
            generate_3d_expression(
                &parsed,
                &Generate3DExpressionOptions {
                    energy_degree_bounds: Some(Vec::new()),
                    ..options.clone()
                }
            )
            .is_ok()
        );
        assert!(matches!(
            generate_3d_expression(
                &parsed,
                &Generate3DExpressionOptions {
                    energy_degree_bounds: Some(vec![(0, usize::MAX), (1, usize::MAX)]),
                    ..options
                }
            ),
            Err(GenerationError::CoefficientOutOfRange)
        ));
    }

    #[test]
    fn ordinary_ltd_does_not_use_numerator_capacity() {
        let parsed = test_graphs::box_graph();
        let options = Generate3DExpressionOptions {
            representation: RepresentationMode::Ltd,
            ..Default::default()
        };
        let ordinary = generate_3d_expression(&parsed, &options).unwrap();
        let bounded = generate_3d_expression(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(0, 50)]),
                ..options
            },
        )
        .unwrap();
        assert_eq!(
            ordinary
                .expression
                .to_atom(crate::expression::AllOrientations),
            bounded
                .expression
                .to_atom(crate::expression::AllOrientations)
        );
        assert!(ordinary.expression.orientations.iter().all(|orientation| {
            edge_q0_from_loop_exprs(
                &parsed
                    .internal_edges
                    .iter()
                    .map(|edge| edge.signature.clone())
                    .collect::<Vec<_>>(),
                &orientation.loop_energy_map,
            ) == orientation.edge_energy_map
        }));
    }

    #[cfg(feature = "eval")]
    fn contour_value(
        parsed: &ParsedGraph,
        generated: &GeneratedThreeDExpression,
        numerator: &str,
        input: &crate::eval::EvaluationInput,
    ) -> f64 {
        let frame = if generated.representation == RepresentationMode::Ltd {
            1
        } else {
            generated
                .energy_factor_components
                .iter()
                .map(|component| {
                    let source_frame = match component.ownership {
                        CffEnergyFactorOwnership::GlobalSourceProduct => {
                            component.core_global_prefactor_sign
                        }
                        CffEnergyFactorOwnership::VariantLocal => {
                            component.denominator_only_global_prefactor_sign
                        }
                    };
                    CffGlobalPrefactorSign::from_exponent(component.internal_edge_ids.len())
                        .product(source_frame)
                        .factor()
                })
                .product()
        };
        frame as f64
            * crate::eval::evaluate_expression(parsed, &generated.expression, numerator, input)
                .unwrap()
                .value
    }

    #[cfg(feature = "eval")]
    #[test]
    fn raised_ltd_matches_analytic_contours_for_all_convergent_monomials() {
        for power in 2..=4 {
            for signs in [
                vec![1; power],
                (0..power)
                    .map(|index| if index % 2 == 0 { -1 } else { 1 })
                    .collect(),
            ] {
                let parsed = repeated_cycle(power, &signs);
                let mass = 0.87_f64;
                for degree in 0..=2 * power - 2 {
                    // Coefficient of y^(power-1) in (E+y)^degree/(2E+y)^power,
                    // with a minus for the clockwise positive-energy contour.
                    let expected = -(0..=degree.min(power - 1))
                        .map(|index| {
                            let remaining = power - 1 - index;
                            let numerator_choose = super::super::binomial(degree, index)
                                .to_i64_pair()
                                .unwrap()
                                .0 as f64;
                            let denominator_choose =
                                super::super::binomial(power + remaining - 1, remaining)
                                    .to_i64_pair()
                                    .unwrap()
                                    .0 as f64;
                            numerator_choose
                                * denominator_choose
                                * if remaining % 2 == 0 { 1.0 } else { -1.0 }
                                * mass.powi((degree - index) as i32)
                                / (2.0 * mass).powi((power + remaining) as i32)
                        })
                        .sum::<f64>();
                    let numerator = format!("({}*edges[0][0])**{degree}", signs[0]);
                    for sampling in [
                        NumeratorSamplingScaleMode::None,
                        NumeratorSamplingScaleMode::All,
                    ] {
                        let generated = generate_3d_expression(
                            &parsed,
                            &Generate3DExpressionOptions {
                                representation: RepresentationMode::Ltd,
                                energy_degree_bounds: Some(vec![(0, degree)]),
                                numerator_sampling_scale: sampling,
                                ..Default::default()
                            },
                        )
                        .unwrap();
                        for scale in [0.43, 1.37] {
                            let input = crate::eval::EvaluationInput {
                                external_momenta: Vec::new(),
                                loop_spatial_momenta: vec![[0.0; 3]],
                                masses: vec![mass; power],
                                uniform_scale: Some(scale),
                            };
                            let actual = contour_value(&parsed, &generated, &numerator, &input);
                            assert!(
                                (actual - expected).abs() < 2.0e-10 * expected.abs().max(1.0),
                                "power={power}, signs={signs:?}, degree={degree}, sampling={sampling:?}, M={scale}: {actual} != {expected}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[cfg(feature = "eval")]
    #[test]
    fn ltd_matches_current_cff_for_simple_and_repeated_sources() {
        let mut forward_bubble = test_graphs::initial_state_cut_line_graph(1);
        forward_bubble.internal_edges[1].tail = 1;
        forward_bubble.internal_edges[1].head = 0;
        let mut second = forward_bubble.internal_edges[1].clone();
        second.edge_id = 2;
        second.label = "p_minus_k".to_string();
        second.signature.loop_signature = vec![-1];
        second.signature.external_signature = vec![1];
        second.mass_key = Some("m_other".to_string());
        forward_bubble.internal_edges.push(second);
        assert!(crate::validate_parsed_graph(&forward_bubble).ok);
        for (parsed, bounds, numerator) in [
            (test_graphs::box_graph(), vec![(0, 2)], "edges[0][0]**2"),
            (
                test_graphs::box_pow3_graph(),
                vec![(3, 2)],
                "edges[3][0]**2",
            ),
            (
                test_graphs::box_pow3_graph(),
                vec![(0, 2)],
                "edges[0][0]**2",
            ),
            (
                test_graphs::sunrise_pow4_graph(),
                vec![(2, 2), (3, 1)],
                "edges[2][0]**2*edges[3][0]",
            ),
            (forward_bubble, Vec::new(), "1"),
            (test_graphs::pure_tree_graph(), Vec::new(), "1"),
        ] {
            let options = Generate3DExpressionOptions {
                energy_degree_bounds: Some(bounds.clone()),
                numerator_sampling_scale: NumeratorSamplingScaleMode::All,
                preserve_internal_edges_as_four_d_denominators: if parsed.loop_names.is_empty() {
                    parsed.denominator_internal_edge_ids()
                } else {
                    Vec::new()
                },
                ..Default::default()
            };
            let cff = generate_3d_expression(&parsed, &options).unwrap_or_else(|error| {
                panic!(
                    "CFF source loops={}, edges={}, bounds={bounds:?}: {error}",
                    parsed.loop_names.len(),
                    parsed.internal_edges.len()
                )
            });
            let ltd = generate_3d_expression(
                &parsed,
                &Generate3DExpressionOptions {
                    representation: RepresentationMode::Ltd,
                    ..options
                },
            )
            .unwrap();
            for seed in [5, 37, 291] {
                let input = crate::eval::EvaluationInput::deterministic(
                    &parsed,
                    seed,
                    &BTreeMap::new(),
                    Some(0.91),
                )
                .unwrap();
                let expected = contour_value(&parsed, &cff, numerator, &input);
                let actual = contour_value(&parsed, &ltd, numerator, &input);
                assert!(
                    (actual - expected).abs() < 2.0e-8 * actual.abs().max(expected.abs()).max(1.0),
                    "loops={}, edges={}, numerator={numerator}, bounds={bounds:?}, seed={seed}: LTD={actual}, CFF={expected}",
                    parsed.loop_names.len(),
                    parsed.internal_edges.len()
                );
            }
        }
    }

    #[cfg(feature = "eval")]
    #[test]
    fn raised_ltd_preserves_products_and_mixed_derivative_axes() {
        for connected in [false, true] {
            let mut parsed = repeated_cycle(2, &[1, -1]);
            for edge in &mut parsed.internal_edges {
                edge.signature.loop_signature.push(0);
            }
            let offset_node = |node| {
                if connected && node == 0 {
                    0
                } else {
                    node + 2 - usize::from(connected)
                }
            };
            for mut edge in repeated_cycle(3, &[-1, 1, 1]).internal_edges {
                edge.edge_id += 2;
                edge.tail = offset_node(edge.tail);
                edge.head = offset_node(edge.head);
                edge.mass_key = Some("n".to_string());
                let sign = edge.signature.loop_signature[0];
                edge.signature.loop_signature = vec![if connected { sign } else { 0 }, sign];
                parsed.internal_edges.push(edge);
            }
            parsed.loop_names.push("l".to_string());
            parsed.node_name_to_internal = parsed
                .internal_edges
                .iter()
                .flat_map(|edge| [edge.tail, edge.head])
                .map(|node| (format!("v{node}"), node))
                .collect();
            let input = crate::eval::EvaluationInput {
                external_momenta: Vec::new(),
                loop_spatial_momenta: vec![[0.0; 3]; 2],
                masses: vec![0.83, 0.83, 1.17, 1.17, 1.17],
                uniform_scale: Some(0.71),
            };
            let expected = -1.0 / (64.0 * 0.83 * 1.17_f64.powi(3));
            for sampling in [
                NumeratorSamplingScaleMode::None,
                NumeratorSamplingScaleMode::All,
            ] {
                let generated = generate_3d_expression(
                    &parsed,
                    &Generate3DExpressionOptions {
                        representation: RepresentationMode::Ltd,
                        energy_degree_bounds: Some(vec![(0, 2), (2, 2)]),
                        numerator_sampling_scale: sampling,
                        ..Default::default()
                    },
                )
                .unwrap();
                let actual =
                    contour_value(&parsed, &generated, "edges[0][0]**2*edges[2][0]**2", &input);
                assert!(
                    (actual - expected).abs() < 1.0e-11 * expected.abs(),
                    "connected={connected}, sampling={sampling:?}: {actual} != {expected}"
                );
            }
        }
    }
}
