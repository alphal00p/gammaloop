//! Normal flows tangent to massless spectator strata where they are transverse.
//!
//! A chart freezes an independent set of spectator momenta. For the remaining
//! coordinates, J is the gradient of h_a = q_a²/2 and G = J Jᵀ. The polynomial
//! field Jᵀ adj(G) e_i satisfies J field = det(G) e_i, even at a singular chart.
//! Multiplying by det(G), then normalizing the sum by Σ t_C det(G)², avoids
//! inverse-Gram poles from individual charts. E_i converts the h_i normal to
//! the energy normal: V_i(E_a) = δ_ia.
//!
//! A chart's t_C vanishes quadratically when an unprotected spectator becomes
//! soft. At a transverse soft stratum, surviving charts preserve that momentum
//! and the complete flow is tangent there. This construction assumes that at
//! least one chart remains transverse; a nontransverse stratum needs its own
//! distributional limit.

use super::*;

impl FermiSurfaceLocalizer<'_> {
    pub(super) fn normal_flow(
        &self,
        support: &[FermiDistribution],
        target: EdgeIndex,
        builder: &ParamBuilder,
        coefficient: &Atom,
    ) -> Result<FermiNormalFlow> {
        use crate::graph::lmb::LMBext;
        use itertools::Itertools;
        use symbolica::{domains::atom::AtomField, tensors::matrix::Matrix};

        let target_index = support
            .iter()
            .position(|factor| factor.edge == target)
            .ok_or_else(|| eyre!("normal-flow target is outside the Fermi support"))?;
        let original_edges = &self.graph.loop_momentum_basis.loop_edges;
        let active_momenta = support
            .iter()
            .map(|factor| self.momentum(factor.edge, builder))
            .collect::<Result<Vec<_>>>()?;
        let target_energy = self.energy(target, builder)?;
        let target_energy2 = target_energy.pow(2);
        let mut candidate_squares = BTreeMap::new();
        for edge in self.graph.iter_edge_ids() {
            let mass = self.graph[edge].mass_atom();
            if (mass.is_zero() || !matches!(mass.as_view(), AtomView::Num(_)))
                && !support.iter().any(|factor| factor.edge == edge)
            {
                candidate_squares.insert(edge, self.energy(edge, builder)?.pow(2));
            }
        }
        let mut candidates = BTreeSet::new();
        coefficient.as_view().visitor(&mut |part| {
            if let AtomView::Pow(power) = part {
                let (base, exponent) = power.get_base_exp();
                if !i64::try_from(exponent).is_ok_and(|power| power >= 0) {
                    for (edge, square) in &candidate_squares {
                        if base == square.as_view() {
                            candidates.insert(*edge);
                        }
                    }
                }
            }
            true
        });
        if candidates.is_empty() {
            let basis = self.basis(support)?;
            let target_loop = LoopIndex(
                basis
                    .loop_edges
                    .iter()
                    .position(|edge| *edge == target)
                    .unwrap(),
            );
            let radial = target_energy / Self::norm_squared(&active_momenta[target_index]);
            let flow = original_edges
                .iter()
                .map(|edge| {
                    let routing =
                        Self::sign_atom(basis.edge_signatures[*edge].internal[target_loop]);
                    std::array::from_fn(|axis| {
                        &routing * &radial * &active_momenta[target_index][axis]
                    })
                })
                .collect();
            return Ok(FermiNormalFlow {
                velocity: flow,
                energy_square_rates: BTreeMap::new(),
                domains: Vec::new(),
            });
        }

        // A closure determines the protected subspace. Different graph bases
        // spanning the same protected momenta need only contribute one chart.
        let mut charts = BTreeMap::new();
        for basis in self.graph.generate_loop_momentum_bases() {
            let eligible = basis
                .loop_edges
                .iter_enumerated()
                .filter_map(|(index, edge)| candidates.contains(edge).then_some(index))
                .collect::<Vec<_>>();
            for frozen in eligible.iter().copied().powerset() {
                let free = basis
                    .loop_edges
                    .iter_enumerated()
                    .filter_map(|(index, _)| (!frozen.contains(&index)).then_some(index))
                    .collect::<Vec<_>>();
                let closure = candidates
                    .iter()
                    .copied()
                    .filter(|edge| {
                        free.iter().all(|index| {
                            basis.edge_signatures[*edge].internal[*index] == SignOrZero::Zero
                        })
                    })
                    .collect::<Vec<_>>();
                let frozen_edges = frozen
                    .iter()
                    .map(|index| basis.loop_edges[*index])
                    .sorted()
                    .collect::<Vec<_>>();
                let (_, _, soft_bases) = charts
                    .entry(closure)
                    .or_insert_with(|| (basis.clone(), free, BTreeSet::new()));
                soft_bases.insert(frozen_edges);
            }
        }

        let determinant = |entries: Vec<Vec<Atom>>| -> Result<Atom> {
            if entries.is_empty() {
                return Ok(Atom::one());
            }
            let field = AtomField {
                statistical_zero_test: false,
                cancel_check_on_division: true,
                custom_normalization: None,
            };
            let matrix = Matrix::from_nested_vec(entries, field).map_err(|error| eyre!(error))?;
            // Determinants contain kinematics only; exact polynomial expansion
            // exposes rank-deficient charts without touching thermal weights.
            Ok(matrix
                .det()
                .map_err(|error| eyre!("{error}"))?
                .together()
                .expand())
        };
        let mut domains = Vec::new();
        let mut weight = Atom::Zero;
        let mut numerator = vec![[Atom::Zero, Atom::Zero, Atom::Zero]; original_edges.len()];
        let mut soft_velocities = candidates
            .iter()
            .map(|edge| (*edge, [Atom::Zero, Atom::Zero, Atom::Zero]))
            .collect::<BTreeMap<_, _>>();
        for (closure, (basis, free, soft_bases)) in charts {
            // Even charts with no free directions describe a possible soft
            // stratum; certify it before discarding a vanishing Gram matrix.
            if let Some(domain) =
                self.soft_domain(support, &basis, &free, soft_bases.into_iter().collect())?
            {
                domains.push(domain);
            }
            if 3 * free.len() < support.len() {
                continue;
            }
            let gram = active_momenta
                .iter()
                .enumerate()
                .map(|(a, qa)| {
                    active_momenta
                        .iter()
                        .enumerate()
                        .map(|(b, qb)| {
                            let routing_dot = Atom::add_many(free.iter().map(|index| {
                                Self::sign_atom(
                                    basis.edge_signatures[support[a].edge].internal[*index],
                                ) * Self::sign_atom(
                                    basis.edge_signatures[support[b].edge].internal[*index],
                                )
                            }));
                            routing_dot * Atom::add_many(qa.iter().zip(qb).map(|(x, y)| x * y))
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            let det = determinant(gram.clone())?;
            if det.is_zero() {
                continue;
            }
            let damping = candidates
                .iter()
                .filter(|edge| !closure.contains(edge))
                .fold(Atom::one(), |product, edge| {
                    let square = &candidate_squares[edge];
                    product * square / (square + &target_energy2)
                });
            weight += &damping * det.pow(2);
            let adjugate_column = (0..support.len())
                .map(|column| {
                    // adj(G)[column,target] is cofactor(G)[target,column].
                    let minor = gram
                        .iter()
                        .enumerate()
                        .filter(|(row, _)| *row != target_index)
                        .map(|(_, entries)| {
                            entries
                                .iter()
                                .enumerate()
                                .filter_map(|(index, entry)| {
                                    (index != column).then_some(entry.clone())
                                })
                                .collect::<Vec<_>>()
                        })
                        .collect::<Vec<_>>();
                    Ok(Atom::num(if (target_index + column) % 2 == 0 {
                        1
                    } else {
                        -1
                    }) * determinant(minor)?)
                })
                .collect::<Result<Vec<_>>>()?;
            for free_index in free {
                let chart_flow: [Atom; 3] = std::array::from_fn(|axis| {
                    Atom::add_many(support.iter().enumerate().map(|(a, factor)| {
                        Self::sign_atom(basis.edge_signatures[factor.edge].internal[free_index])
                            * &active_momenta[a][axis]
                            * &adjugate_column[a]
                    }))
                });
                // Protected charts have exactly zero spectator velocity.
                // Every remaining chart carries this edge's E² in its damping;
                // cancel it before summing, so the soft factor stays explicit
                // when the complete flow differentiates an energy pole.
                for edge in candidates.iter().filter(|edge| !closure.contains(edge)) {
                    let routing =
                        Self::sign_atom(basis.edge_signatures[*edge].internal[free_index]);
                    if routing.is_zero() {
                        continue;
                    }
                    let damping_without_square = &damping / &candidate_squares[edge];
                    for (axis, component) in chart_flow.iter().enumerate() {
                        soft_velocities.get_mut(edge).unwrap()[axis] +=
                            &damping_without_square * &det * &routing * component;
                    }
                }
                for (old_index, old_edge) in original_edges.iter().enumerate() {
                    let routing =
                        Self::sign_atom(basis.edge_signatures[*old_edge].internal[free_index]);
                    if routing.is_zero() {
                        continue;
                    }
                    for (axis, component) in chart_flow.iter().enumerate() {
                        numerator[old_index][axis] += &damping * &det * &routing * component;
                    }
                }
            }
        }
        if weight.is_zero() {
            return Err(eyre!(
                "no transverse normal-flow chart for Fermi support {support:?}"
            ));
        }
        let flow = numerator
            .into_iter()
            .map(|components| components.map(|component| &target_energy * component / &weight))
            .collect();
        let energy_square_rates = soft_velocities
            .into_iter()
            .map(|(edge, velocity)| {
                let momentum = self.momentum(edge, builder)?;
                Ok((
                    edge,
                    2 * &target_energy
                        * Atom::add_many(velocity.iter().zip(&momentum).map(|(v, q)| v * q))
                        / &weight,
                ))
            })
            .collect::<Result<BTreeMap<_, _>>>()?;
        Ok(FermiNormalFlow {
            velocity: flow,
            energy_square_rates,
            domains,
        })
    }
}
