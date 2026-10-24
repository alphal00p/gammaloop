//! Native graph automorphism proposals, with exact momentum-map replay.
mod edge_orbits;
#[cfg(test)]
mod tests;

use super::*;
use std::collections::BTreeSet;

/// Optional work-count diagnostics for scope-aware native witness selection.
/// Ordinary requests use the shared RAM/time budget instead of these counts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MomentumMatchLimits {
    pub diagnostics_enabled: bool,
    pub max_graph_automorphisms: usize,
    pub max_orientation_candidates: usize,
    pub max_verified_witnesses: usize,
}
impl Default for MomentumMatchLimits {
    fn default() -> Self {
        Self {
            diagnostics_enabled: false,
            max_graph_automorphisms: 100_000,
            max_orientation_candidates: 1_000_000,
            max_verified_witnesses: 100_000,
        }
    }
}
impl MomentumMatchLimits {
    pub fn with_diagnostics(mut self) -> Self {
        self.diagnostics_enabled = true;
        self
    }
}

fn charge(
    count: &mut usize,
    limit: Option<usize>,
    resource: &'static str,
) -> Result<(), VakintError> {
    symbolica::utils::resources::checkpoint()?;
    let requested = count.saturating_add(1);
    if let Some(limit) = limit
        && requested > limit
    {
        return Err(VakintError::MomentumMatchResourceLimit {
            resource,
            requested,
            limit,
        });
    }
    *count = requested;
    Ok(())
}

impl MomentumInput {
    pub fn match_topologies_with_acceptance<'a, K: Ord>(
        &self,
        topologies: impl IntoIterator<Item = &'a Topology>,
        limits: MomentumMatchLimits,
        mut accept: impl FnMut(&ReplacementRules) -> Result<Option<K>, VakintError>,
    ) -> Result<Option<ReplacementRules>, VakintError> {
        let topologies = topologies.into_iter().collect::<Vec<_>>();
        // Retain the ordinary exact witness for partial reduction if every
        // verified alternative is outside the saved entry union. Errors from
        // acceptance are never confused with ordinary out-of-scope rejection.
        let mut fallback = None;
        let mut best = None;
        let source_graph = edge_orbits::EdgeGraph::new(&self.incidence)?;
        let loops = self.family.loop_momenta();
        let variables = (1..=loops.len())
            .map(|axis| Atom::var(vk_symbol!(format!("momentum_match_source_{axis}"))))
            .collect::<Vec<_>>();
        let substitute = loops
            .iter()
            .zip(&variables)
            .map(|(source, variable)| Replacement::new(source.clone(), variable.clone()))
            .collect::<Vec<_>>();
        let sources = self
            .basis
            .iter()
            .map(|&slot| {
                self.propagators[slot]
                    .momentum
                    .replace_multiple(&substitute)
            })
            .collect::<Vec<_>>();
        let orientations = 1usize
            .checked_shl(loops.len().try_into().map_err(invalid)?)
            .ok_or_else(|| invalid("momentum orientation count overflow"))?;
        let mut tried = 0;
        let mut verified = 0;
        let mut graph_maps = 0;
        for topology in topologies {
            let integral = topology.get_integral();
            if integral.n_loops != loops.len()
                || integral.graph.edges.len() != self.propagators.len()
            {
                continue;
            }
            let mut unit_graph = integral.graph.clone();
            for edge in unit_graph.edges.values_mut() {
                edge.mass = Atom::one();
            }
            let target_graph = edge_orbits::EdgeGraph::new(&unit_graph.to_symbolica_graph(false))?;
            if source_graph.canonical.graph != target_graph.canonical.graph {
                continue;
            }
            let expression = integral
                .canonical_expression
                .as_ref()
                .ok_or_else(|| invalid("missing canonical topology"))?;
            let physical = (1..=integral.n_props)
                .filter_map(|id| {
                    get_prop_with_id(expression.as_view(), id)
                        .map(|propagator| (id, propagator[&vk_symbol!("q_")].clone()))
                })
                .collect::<Vec<_>>();
            if physical.len() != self.propagators.len() {
                return Err(invalid("canonical graph and propagators disagree"));
            }
            let kinematics = self.family.kinematics();
            let denominators = physical
                .iter()
                .map(|(_, momentum)| {
                    kinematics
                        .scalar_product(momentum, momentum)
                        .map(|s| s - Atom::one())
                })
                .collect::<Result<Vec<_>, _>>()
                .map_err(invalid)?;
            let target = IntegralFamily::new(loops.to_vec(), Vec::new(), denominators, kinematics)
                .map_err(invalid)?;
            let mut seen_coordinates = BTreeSet::new();
            let mut seen_permutations = BTreeSet::new();
            let mut evaluate = |permutation: &[usize]| {
                if !seen_permutations.insert(permutation.to_vec()) {
                    return Ok(());
                }
                let images = self
                    .basis
                    .iter()
                    .map(|&slot| &physical[permutation[slot]].1)
                    .collect::<Vec<_>>();
                for signs in 0..orientations {
                    charge(
                        &mut tried,
                        limits
                            .diagnostics_enabled
                            .then_some(limits.max_orientation_candidates),
                        "momentum orientation candidates",
                    )?;
                    let system = sources
                        .iter()
                        .zip(&images)
                        .enumerate()
                        .map(|(axis, (source, image))| {
                            source
                                - Atom::num(if signs & (1 << axis) == 0 { 1 } else { -1 })
                                    * (*image)
                        })
                        .collect::<Vec<_>>();
                    let coordinates = Topology::solve_loop_coordinates(&system, &variables)?;
                    if !seen_coordinates.insert(coordinates.clone()) {
                        continue;
                    }
                    let Some(mapping) = self
                        .family
                        .mapping_to(&target, &coordinates)
                        .map_err(invalid)?
                    else {
                        continue;
                    };
                    charge(
                        &mut verified,
                        limits
                            .diagnostics_enabled
                            .then_some(limits.max_verified_witnesses),
                        "verified momentum witnesses",
                    )?;
                    let candidate = self.replacement(topology, &physical, &mapping);
                    if verified <= 8 {
                        log::debug!(target: "vakint::five_loop_match",
                            "phase=verified_alternative ordinal={verified} topology={} source_rows={} target_rows={} denominator_map={:?} coordinates={:?}",
                            integral.name,self.propagators.len(),physical.len(),mapping.denominator_map(),coordinates);
                    }
                    if fallback.is_none() {
                        fallback = Some(candidate.clone());
                    }
                    let accepted = accept(&candidate)?;
                    if verified <= 8 {
                        log::debug!(target: "vakint::five_loop_match", "phase=alternative_admission ordinal={verified} admitted={}",accepted.is_some());
                    }
                    if let Some(key) = accepted
                        && best.as_ref().is_none_or(|(old, _)| key < *old)
                    {
                        best = Some((key, candidate));
                    }
                }
                Ok(())
            };
            // Try precisely the ordinary incidence witness first, under the
            // same budget, so a complete scope miss has a deterministic fallback.
            if let Some(first) = self.incidence_permutation(&integral.graph) {
                charge(
                    &mut graph_maps,
                    limits
                        .diagnostics_enabled
                        .then_some(limits.max_graph_automorphisms),
                    "native edge automorphism witnesses",
                )?;
                evaluate(&first)?;
            }
            edge_orbits::visit(
                &source_graph,
                &target_graph,
                limits
                    .diagnostics_enabled
                    .then_some(limits.max_graph_automorphisms),
                &mut graph_maps,
                &mut evaluate,
            )?;
        }
        log::debug!(target:"vakint::five_loop_match", "phase=admitted_witness_search graph_maps={graph_maps} orientation_candidates={tried} verified_witnesses={verified} admitted={}",best.is_some());
        Ok(best.map(|(_, candidate)| candidate).or(fallback))
    }
}
