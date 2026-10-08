//! Laurent coefficients of LTD poles, with the graph numerator kept opaque.

use std::collections::BTreeSet;

use super::*;
use crate::{
    generation::{GenerationError, classify_surface_kind, lagrange_basis},
    surface::{LinearSurface, LinearSurfaceID, SurfaceOrigin},
    utils::{RationalExt, solve_rational_system},
};

/// Independent affine energy used while extracting a physical pole.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum ResidueEnergyCoordinate {
    External(EdgeIndex),
    OnShell(EdgeIndex),
}

impl ResidueEnergyCoordinate {
    fn coefficient(self, expression: &LinearEnergyExpr) -> Rational {
        let (terms, edge) = match self {
            Self::External(edge) => (&expression.external_terms, edge),
            Self::OnShell(edge) => (&expression.internal_terms, edge),
        };
        terms
            .iter()
            .filter(|(id, _)| *id == edge)
            .fold(Rational::zero(), |sum, (_, coefficient)| sum + coefficient)
    }
}

impl<E: Clone, H: Clone> ThreeDExpression<OrientationID, E, H> {
    /// Extract the principal part in the canonical selected energy surface.
    ///
    /// Scalar factors are expanded exactly. Polynomial numerator coefficients
    /// are finite linear combinations of evaluations at affine energy maps;
    /// neither the numerator body nor its derivatives enter this operation.
    /// The affine direction keeps the other physical residue coordinates
    /// fixed, so mixed pole orders use the existing radial derivative contract.
    /// Bounds index edge-map slots first, followed by loop-map slots; `None`
    /// means unknown, whereas `Some(&[])` certifies energy independence.
    pub fn select_ltd_residue(
        &self,
        surface_expressions: &BTreeMap<HybridSurfaceID, LinearEnergyExpr>,
        selected_surface: &LinearEnergyExpr,
        other_selected_surfaces: &[LinearEnergyExpr],
        energy_degree_bounds: Option<&[(usize, usize)]>,
    ) -> Result<Vec<Self>, GenerationError> {
        if !self.residual_denominators.is_empty() {
            return Err(GenerationError::InvalidResidueCoordinate(
                "residual four-dimensional denominators must be resolved before pole selection"
                    .to_owned(),
            ));
        }
        let selected_surface = selected_surface.clone().canonical();
        let axes = std::iter::once(&selected_surface)
            .chain(other_selected_surfaces)
            .collect_vec();
        let coordinates = axes
            .iter()
            .flat_map(|surface| {
                surface
                    .external_terms
                    .iter()
                    .map(|(edge, _)| ResidueEnergyCoordinate::External(*edge))
                    .chain(
                        surface
                            .internal_terms
                            .iter()
                            .map(|(edge, _)| ResidueEnergyCoordinate::OnShell(*edge)),
                    )
            })
            .collect::<BTreeSet<_>>();
        let direction = coordinates
            .into_iter()
            .combinations(axes.len())
            .find_map(|coordinates| {
                let matrix = axes
                    .iter()
                    .map(|surface| {
                        coordinates
                            .iter()
                            .map(|coordinate| coordinate.coefficient(surface))
                            .collect()
                    })
                    .collect();
                let rhs = (0..axes.len())
                    .map(|index| Rational::from(i32::from(index == 0)))
                    .collect();
                solve_rational_system(matrix, rhs)
                    .map(|direction| coordinates.into_iter().zip(direction).collect_vec())
            })
            .ok_or_else(|| {
                GenerationError::InvalidResidueCoordinate(
                    "selected physical surfaces must be linearly independent energy coordinates"
                        .to_owned(),
                )
            })?;
        let slope = |expression: &LinearEnergyExpr| {
            direction
                .iter()
                .fold(Rational::zero(), |sum, (coordinate, direction)| {
                    sum + coordinate.coefficient(expression) * direction
                })
        };
        let project = |expression: LinearEnergyExpr| {
            let derivative = slope(&expression);
            expression - selected_surface.clone().scale_rational(derivative)
        };
        let mut surfaces = self.surfaces.clone();
        let mut surface_ids = surfaces
            .linear_surface_cache
            .iter_enumerated()
            .map(|(id, surface)| (surface.expression.clone(), HybridSurfaceID::Linear(id)))
            .collect::<HashMap<_, _>>();
        let mut intern = |expression: LinearEnergyExpr| {
            *surface_ids.entry(expression.clone()).or_insert_with(|| {
                let id = LinearSurfaceID(surfaces.linear_surface_cache.len());
                surfaces.linear_surface_cache.push(LinearSurface {
                    kind: classify_surface_kind(&expression),
                    expression,
                    origin: SurfaceOrigin::Physical,
                    numerator_only: false,
                });
                HybridSurfaceID::Linear(id)
            })
        };
        let mut coefficients = Vec::<Vec<OrientationExpression>>::new();

        for orientation in &self.orientations {
            let map_slopes = orientation.edge_energy_map.iter().map(&slope).collect_vec();
            let loop_slopes = orientation.loop_energy_map.iter().map(&slope).collect_vec();
            let numerator_slopes = map_slopes.iter().chain(&loop_slopes).collect_vec();
            let degree = energy_degree_bounds
                .unwrap_or_default()
                .iter()
                .filter(|(edge, _)| {
                    numerator_slopes
                        .get(*edge)
                        .is_some_and(|value| !value.is_zero())
                })
                .try_fold(0usize, |sum, (_, degree)| sum.checked_add(*degree))
                .ok_or(GenerationError::CoefficientOutOfRange)?;
            let zero_edge_maps = orientation
                .edge_energy_map
                .iter()
                .cloned()
                .map(&project)
                .collect_vec();
            let zero_loop_maps = orientation
                .loop_energy_map
                .iter()
                .cloned()
                .map(&project)
                .collect_vec();
            let nodes = (0..=degree)
                .map(|index| {
                    let magnitude = i32::try_from(index.div_ceil(2))
                        .map_err(|_| GenerationError::CoefficientOutOfRange)?;
                    Ok(if index % 2 == 0 {
                        -magnitude
                    } else {
                        magnitude
                    })
                })
                .collect::<Result<Vec<_>, GenerationError>>()?;
            let interpolation = nodes
                .iter()
                .enumerate()
                .map(|(index, node)| (*node, lagrange_basis(&nodes, index)))
                .collect_vec();

            for variant in &orientation.variants {
                for chain in denominator_tree_chains(&variant.denominator) {
                    let mut pole_order = 0isize;
                    let mut prefactor = variant.prefactor.clone();
                    let mut regular_factors = Vec::new();
                    let mut retained_half_edges = Vec::new();
                    for edge in &variant.half_edges {
                        let expression = LinearEnergyExpr::ose(*edge, 1);
                        let derivative = slope(&expression);
                        if derivative.is_zero() {
                            retained_half_edges.push(*edge);
                        } else {
                            let constant = project(expression);
                            prefactor /= Rational::from(2);
                            if constant.is_zero() {
                                pole_order += 1;
                                prefactor /= derivative;
                            } else {
                                regular_factors.push((intern(constant), derivative, true));
                            }
                        }
                    }
                    for (id, denominator) in chain.into_iter().map(|id| (id, true)).chain(
                        variant
                            .numerator_surfaces
                            .iter()
                            .copied()
                            .map(|id| (id, false)),
                    ) {
                        let expression = surface_expressions.get(&id).ok_or_else(|| {
                            GenerationError::InvalidResidueCoordinate(format!(
                                "missing exact expression for surface {id:?}"
                            ))
                        })?;
                        let derivative = slope(expression);
                        let constant = project(expression.clone());
                        if constant.is_zero() {
                            if derivative.is_zero() {
                                return Err(GenerationError::InvalidResidueCoordinate(format!(
                                    "surface {id:?} is identically zero"
                                )));
                            }
                            if denominator {
                                pole_order += 1;
                                prefactor /= derivative;
                            } else {
                                pole_order -= 1;
                                prefactor *= derivative;
                            }
                        } else {
                            regular_factors.push((intern(constant), derivative, denominator));
                        }
                    }
                    if pole_order <= 0 {
                        continue;
                    }
                    let pole_order = pole_order as usize;
                    if pole_order > 1
                        && energy_degree_bounds.is_none()
                        && numerator_slopes.iter().any(|slope| !slope.is_zero())
                    {
                        return Err(GenerationError::LtdResidueRequiresEnergyBounds);
                    }
                    coefficients.resize_with(coefficients.len().max(pole_order), Vec::new);
                    let mut base = variant.clone();
                    base.prefactor = prefactor;
                    base.half_edges = retained_half_edges;
                    base.denominator = Tree::from_root(HybridSurfaceID::Unit);
                    base.numerator_surfaces.clear();
                    // The exact scalar prefactor already contains these
                    // signs. Newly expanded linear factors have no CFF
                    // selected-denominator sign convention to consume.
                    base.denominator_surface_signs.clear();
                    let mut jets = vec![Vec::new(); pole_order];
                    jets[0].push(base);
                    for (id, derivative, denominator) in regular_factors {
                        let mut next = vec![Vec::new(); pole_order];
                        for (order, terms) in jets.into_iter().enumerate() {
                            let max_extra = if derivative.is_zero() {
                                0
                            } else if denominator {
                                pole_order - 1 - order
                            } else {
                                usize::from(order + 1 < pole_order)
                            };
                            for extra in 0..=max_extra {
                                for term in &terms {
                                    let mut term = term.clone();
                                    if denominator {
                                        term.prefactor *= (-derivative.clone()).pow_usize(extra);
                                        for _ in 0..=extra {
                                            let parent = term.denominator.get_bottom_layer()[0];
                                            term.denominator.insert_node(parent, id);
                                        }
                                    } else if extra == 0 {
                                        term.numerator_surfaces.push(id);
                                    } else {
                                        term.prefactor *= &derivative;
                                    }
                                    next[order + extra].push(term);
                                }
                            }
                        }
                        jets = next;
                    }
                    for (scalar_order, terms) in jets.into_iter().enumerate() {
                        if terms.is_empty() {
                            continue;
                        }
                        for numerator_order in 0..=degree.min(pole_order - scalar_order - 1) {
                            let samples = if numerator_order == 0 {
                                vec![(0, Rational::one())]
                            } else {
                                interpolation
                                    .iter()
                                    .filter_map(|(node, weights)| {
                                        weights
                                            .get(numerator_order)
                                            .filter(|weight| !weight.is_zero())
                                            .map(|weight| (*node, weight.clone()))
                                    })
                                    .collect_vec()
                            };
                            for (node, weight) in samples {
                                let sample_maps =
                                    |maps: &[LinearEnergyExpr], slopes: &[Rational]| {
                                        maps.iter()
                                            .zip(slopes)
                                            .map(|(map, slope)| {
                                                map.clone()
                                                    + LinearEnergyExpr::uniform_scale_with_coeff(
                                                        slope * &Rational::from(node),
                                                    )
                                            })
                                            .collect_vec()
                                    };
                                let mut data = orientation.data.clone();
                                data.numerator_map_index = None;
                                coefficients[pole_order - scalar_order - numerator_order - 1].push(
                                    OrientationExpression {
                                        data,
                                        edge_energy_map: sample_maps(&zero_edge_maps, &map_slopes),
                                        loop_energy_map: sample_maps(&zero_loop_maps, &loop_slopes),
                                        variants: terms
                                            .iter()
                                            .map(|term| {
                                                let mut term = term.clone();
                                                term.prefactor *= &weight;
                                                term.uniform_scale_power += numerator_order;
                                                term
                                            })
                                            .collect(),
                                    },
                                );
                            }
                        }
                    }
                }
            }
        }
        let coefficients = coefficients
            .into_iter()
            .map(|orientations| {
                let mut expression = Self {
                    orientations: orientations.into_iter().collect(),
                    surfaces: surfaces.clone(),
                    residual_denominators: self.residual_denominators.clone(),
                }
                .fuse_compatible_variants();
                assign_numerator_map_labels(&mut expression.orientations);
                expression
            })
            .collect();
        Ok(coefficients)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{surface::LinearSurfaceKind, symbols::numerator_sampling_scale};

    #[test]
    fn laurent_selection_retains_numerator_and_denominator_jets_exactly() {
        let a =
            LinearEnergyExpr::ose(EdgeIndex(0), 1) + LinearEnergyExpr::external(EdgeIndex(0), 1);
        let b =
            LinearEnergyExpr::ose(EdgeIndex(1), 1) + LinearEnergyExpr::external(EdgeIndex(1), 1);
        let numerator_map = LinearEnergyExpr::ose(EdgeIndex(2), 1)
            + LinearEnergyExpr::external(EdgeIndex(0), 2)
            + LinearEnergyExpr::external(EdgeIndex(1), 3);
        let n0 = LinearEnergyExpr::ose(EdgeIndex(2), 1) - LinearEnergyExpr::ose(EdgeIndex(0), 2)
            + LinearEnergyExpr::external(EdgeIndex(1), 3);
        let evaluate = |expression: &ThreeDExpression<OrientationID>| {
            expression
                .orientations
                .iter()
                .fold(Atom::Zero, |sum, orientation| {
                    sum + orientation.edge_energy_map[0].to_atom(&[]).pow(2)
                        * orientation
                            .to_atom()
                            .replace_multiple(expression.surfaces.get_all_replacements(&[]))
                })
        };
        for sign in [-1, 1] {
            for power in 1..=3 {
                let first = HybridSurfaceID::Linear(LinearSurfaceID(0));
                let second = HybridSurfaceID::Linear(LinearSurfaceID(1));
                let factors = BTreeMap::from([
                    (first, a.clone().scale(sign)),
                    (second, b.clone() + a.clone().scale(5)),
                ]);
                let mut source = ThreeDExpression::<OrientationID>::new_empty();
                for expression in factors.values() {
                    source.surfaces.linear_surface_cache.push(LinearSurface {
                        kind: LinearSurfaceKind::Esurface,
                        expression: expression.clone(),
                        origin: SurfaceOrigin::Physical,
                        numerator_only: false,
                    });
                }
                let mut chain = vec![first; power];
                chain.push(second);
                source.orientations.push(OrientationExpression {
                    data: OrientationData::new(EdgeVec::from_iter([Orientation::Default])),
                    loop_energy_map: vec![numerator_map.clone()],
                    edge_energy_map: vec![numerator_map.clone()],
                    variants: vec![CFFVariant {
                        origin: None,
                        prefactor: Rational::one(),
                        half_edges: Vec::new(),
                        denominator_edges: vec![EdgeIndex(0)],
                        denominator_surface_signs: BTreeMap::new(),
                        denominator_edge_support_signs: BTreeMap::new(),
                        uniform_scale_power: 0,
                        numerator_surfaces: Vec::new(),
                        denominator: denominator_tree_from_chains(&[chain]),
                    }],
                });
                if power == 2 && sign == 1 {
                    // Selecting A must not manufacture an A^-1 B^-2 term
                    // from A^-2 B^-1 when A and B share an external energy.
                    let shared_b = LinearEnergyExpr::ose(EdgeIndex(1), 1)
                        + LinearEnergyExpr::external(EdgeIndex(0), 1);
                    let mut independent_axes = source.clone();
                    independent_axes.surfaces.linear_surface_cache[LinearSurfaceID(1)].expression =
                        shared_b.clone();
                    let factors = BTreeMap::from([(first, a.clone()), (second, shared_b.clone())]);
                    let coefficients = independent_axes
                        .select_ltd_residue(&factors, &a, &[shared_b], Some(&[(0, 2)]))
                        .unwrap();
                    assert!(coefficients[0].orientations.is_empty());
                    assert!(!coefficients[1].orientations.is_empty());
                }
                let residues = source
                    .select_ltd_residue(&factors, &a, &[], Some(&[(0, 2)]))
                    .unwrap();
                let b_atom = b.to_atom(&[]);
                let n0_atom = n0.to_atom(&[]);
                let jet = [
                    n0_atom.clone().pow(2) / &b_atom,
                    Atom::num(4) * &n0_atom / &b_atom
                        - Atom::num(5) * n0_atom.clone().pow(2) / b_atom.clone().pow(2),
                    Atom::num(4) / &b_atom - Atom::num(20) * &n0_atom / b_atom.clone().pow(2)
                        + Atom::num(25) * n0_atom.clone().pow(2) / b_atom.clone().pow(3),
                ];
                assert_eq!(residues.len(), power);
                for (index, residue) in residues.iter().enumerate() {
                    let actual = evaluate(residue);
                    let expected = Atom::num(sign).pow(power as i64) * &jet[power - index - 1];
                    assert!(
                        (actual.clone() - &expected).together().is_zero(),
                        "signed order {power}, coefficient {}: {actual} != {expected}",
                        index + 1,
                    );
                    for scale in [1, 3] {
                        assert!(
                            (actual
                                .clone()
                                .replace(numerator_sampling_scale())
                                .with(Atom::num(scale))
                                - &expected)
                                .together()
                                .is_zero()
                        );
                    }
                }
                if power == 2 {
                    let first_residue = residues[0].clone();
                    let next_factors = first_residue
                        .surfaces
                        .linear_surface_cache
                        .iter_enumerated()
                        .map(|(id, surface)| {
                            (HybridSurfaceID::Linear(id), surface.expression.clone())
                        })
                        .collect();
                    let nested = first_residue
                        .select_ltd_residue(&next_factors, &b, &[], Some(&[(0, 2)]))
                        .unwrap();
                    let n00 = n0
                        .clone()
                        .substitute_external_energy(
                            EdgeIndex(1),
                            &LinearEnergyExpr::ose(EdgeIndex(1), -1),
                        )
                        .to_atom(&[]);
                    // Res_A Res_B N(A,B)^2/[A²(B+5A)] = (4−30)N(0,0).
                    assert!(
                        (evaluate(&nested[0]) + Atom::num(26) * n00)
                            .together()
                            .is_zero()
                    );
                }
            }
        }
    }

    #[test]
    fn merged_physical_poles_request_unknown_numerator_bounds_lazily() {
        let a =
            LinearEnergyExpr::ose(EdgeIndex(0), 1) + LinearEnergyExpr::external(EdgeIndex(0), 1);
        let b =
            LinearEnergyExpr::ose(EdgeIndex(1), 1) + LinearEnergyExpr::external(EdgeIndex(1), 1);
        let numerator_map = LinearEnergyExpr::ose(EdgeIndex(2), 1)
            + LinearEnergyExpr::external(EdgeIndex(0), 2)
            + LinearEnergyExpr::external(EdgeIndex(1), 3);
        let factors = [a.clone(), b.clone() + a.clone(), b.clone() - a.clone()];
        let table = factors
            .iter()
            .enumerate()
            .map(|(index, expression)| {
                (
                    HybridSurfaceID::Linear(LinearSurfaceID(index)),
                    expression.clone(),
                )
            })
            .collect();
        let mut source = ThreeDExpression::<OrientationID>::new_empty();
        for (index, expression) in factors.into_iter().enumerate() {
            source.surfaces.linear_surface_cache.push(LinearSurface {
                kind: if index == 2 {
                    LinearSurfaceKind::Hsurface
                } else {
                    LinearSurfaceKind::Esurface
                },
                expression,
                origin: SurfaceOrigin::Physical,
                numerator_only: false,
            });
        }
        source.orientations.push(OrientationExpression {
            data: OrientationData::new(EdgeVec::from_iter([Orientation::Default])),
            edge_energy_map: vec![numerator_map.clone()],
            loop_energy_map: vec![numerator_map],
            variants: vec![CFFVariant {
                origin: None,
                prefactor: Rational::one(),
                half_edges: Vec::new(),
                denominator_edges: vec![EdgeIndex(0)],
                denominator_surface_signs: BTreeMap::new(),
                denominator_edge_support_signs: BTreeMap::new(),
                uniform_scale_power: 0,
                numerator_surfaces: Vec::new(),
                denominator: denominator_tree_from_chains(&[(0..3)
                    .map(|index| HybridSurfaceID::Linear(LinearSurfaceID(index)))
                    .collect()]),
            }],
        });
        let first = source
            .select_ltd_residue(&table, &a, std::slice::from_ref(&b), None)
            .unwrap()
            .remove(0);
        let table = first
            .surfaces
            .linear_surface_cache
            .iter_enumerated()
            .map(|(id, surface)| (HybridSurfaceID::Linear(id), surface.expression.clone()))
            .collect();
        assert!(matches!(
            first.select_ltd_residue(&table, &b, std::slice::from_ref(&a), None),
            Err(GenerationError::LtdResidueRequiresEnergyBounds)
        ));
        let constant = first
            .select_ltd_residue(&table, &b, std::slice::from_ref(&a), Some(&[]))
            .unwrap();
        assert!(constant[0].orientations.is_empty());
        let n00 = (LinearEnergyExpr::ose(EdgeIndex(2), 1)
            - LinearEnergyExpr::ose(EdgeIndex(0), 2)
            - LinearEnergyExpr::ose(EdgeIndex(1), 3))
        .to_atom(&[]);
        for bound_slot in [0, 1] {
            let mut first = first.clone();
            if bound_slot == 1 {
                for orientation in &mut first.orientations {
                    orientation.edge_energy_map[0] = LinearEnergyExpr::zero();
                }
                assert!(matches!(
                    first.select_ltd_residue(&table, &b, std::slice::from_ref(&a), None),
                    Err(GenerationError::LtdResidueRequiresEnergyBounds)
                ));
            }
            let residues = first
                .select_ltd_residue(
                    &table,
                    &b,
                    std::slice::from_ref(&a),
                    Some(&[(bound_slot, 2)]),
                )
                .unwrap();
            for (expression, expected) in residues
                .iter()
                .zip([Atom::num(6) * &n00, n00.clone().pow(2)])
            {
                let actual = expression
                    .orientations
                    .iter()
                    .fold(Atom::Zero, |sum, orientation| {
                        let map = if bound_slot == 0 {
                            &orientation.edge_energy_map[0]
                        } else {
                            &orientation.loop_energy_map[0]
                        };
                        sum + map.to_atom(&[]).pow(2)
                            * orientation
                                .to_atom()
                                .replace_multiple(expression.surfaces.get_all_replacements(&[]))
                    });
                assert!((actual - expected).together().is_zero());
            }
        }
    }

    #[test]
    fn nested_shared_external_coordinate_includes_on_shell_energy_factor_jets() {
        let a =
            LinearEnergyExpr::ose(EdgeIndex(0), 1) + LinearEnergyExpr::external(EdgeIndex(0), 1);
        let b =
            LinearEnergyExpr::ose(EdgeIndex(1), 1) + LinearEnergyExpr::external(EdgeIndex(0), 1);
        let numerator_map =
            LinearEnergyExpr::ose(EdgeIndex(2), 1) + LinearEnergyExpr::external(EdgeIndex(0), 2);
        let first = HybridSurfaceID::Linear(LinearSurfaceID(0));
        let second = HybridSurfaceID::Linear(LinearSurfaceID(1));
        let factors =
            BTreeMap::from([(first, a.clone()), (second, b.clone() + a.clone().scale(5))]);
        let mut source = ThreeDExpression::<OrientationID>::new_empty();
        for expression in factors.values() {
            source.surfaces.linear_surface_cache.push(LinearSurface {
                kind: LinearSurfaceKind::Esurface,
                expression: expression.clone(),
                origin: SurfaceOrigin::Physical,
                numerator_only: false,
            });
        }
        source.orientations.push(OrientationExpression {
            data: OrientationData::new(EdgeVec::from_iter([Orientation::Default])),
            loop_energy_map: vec![numerator_map.clone()],
            edge_energy_map: vec![numerator_map],
            variants: vec![CFFVariant {
                origin: None,
                prefactor: Rational::one(),
                half_edges: vec![EdgeIndex(0)],
                denominator_edges: vec![EdgeIndex(0)],
                denominator_surface_signs: BTreeMap::new(),
                denominator_edge_support_signs: BTreeMap::new(),
                uniform_scale_power: 0,
                numerator_surfaces: Vec::new(),
                denominator: denominator_tree_from_chains(&[vec![first, first, second]]),
            }],
        });
        let mut first_residues = source
            .select_ltd_residue(&factors, &a, std::slice::from_ref(&b), Some(&[(0, 2)]))
            .unwrap();
        let first_residue = first_residues.remove(0);
        let next_factors = first_residue
            .surfaces
            .linear_surface_cache
            .iter_enumerated()
            .map(|(id, surface)| (HybridSurfaceID::Linear(id), surface.expression.clone()))
            .collect();
        let residues = first_residue
            .select_ltd_residue(&next_factors, &b, &[a], Some(&[(0, 2)]))
            .unwrap();
        let actual = residues[0]
            .orientations
            .iter()
            .fold(Atom::Zero, |sum, orientation| {
                sum + orientation.edge_energy_map[0].to_atom(&[]).pow(2)
                    * orientation
                        .to_atom()
                        .replace_multiple(residues[0].surfaces.get_all_replacements(&[]))
            });
        let e = LinearEnergyExpr::ose(EdgeIndex(1), 1).to_atom(&[]);
        let n = (LinearEnergyExpr::ose(EdgeIndex(2), 1) - LinearEnergyExpr::ose(EdgeIndex(1), 2))
            .to_atom(&[]);
        let expected = -Atom::num(10) * &n / &e - Atom::num(3) * n.pow(2) / e.pow(2);
        assert!(
            (actual.clone() - &expected).together().is_zero(),
            "{actual} != {expected}"
        );
    }
}
