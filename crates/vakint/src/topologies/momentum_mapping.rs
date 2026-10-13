//! Native exact momentum equivalence for the five-loop long-form registry.
//! Distinct attachments of vacuum factors can have identical integral families.
//! FeynKit verifies every denominator and the measure-preserving Jacobian;
//! Vakint only translates that witness into its existing replacement rules.

mod alternatives;
pub use alternatives::MomentumMatchLimits;

use std::{
    collections::{BTreeMap, BTreeSet},
    time::Instant,
};

use feynkit_graph::IntegralFamily;
use feynkit_kinematics::Kinematics;
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::graph::{CanonicalForm, Empty, Graph as NativeGraph};
use symbolica::id::Replacement;

use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{
    ReplacementRules, Topology, VakintError, get_individual_momenta_from_atom,
    get_integer_from_atom, get_prop_with_id, symbols::S,
};

struct Propagator {
    id: usize,
    nodes: (i64, i64),
    momentum: Atom,
    mass: Atom,
    power: Atom,
}

pub(super) struct MomentumInput {
    family: IntegralFamily,
    propagators: Vec<Propagator>,
    basis: Vec<usize>,
    incidence: NativeGraph<Empty, Atom>,
    canonical_incidence: CanonicalForm<Empty, Atom>,
}

#[derive(Clone, Copy, Debug)]
enum Proposal {
    Incidence,
    Parametric,
}

fn invalid(error: impl std::fmt::Display) -> VakintError {
    VakintError::InvalidIntegralFormat(format!("native momentum matching: {error}"))
}

impl MomentumInput {
    pub fn from_expression(
        expression: AtomView,
        loop_count: usize,
    ) -> Result<Option<Self>, VakintError> {
        let pattern = vk_parse!("prop(id_,edge(nl_,nr_),q_,m_,power_)")
            .unwrap()
            .to_pattern();
        let mut propagators = Vec::new();
        let mut momenta = BTreeSet::new();
        for matched in expression.pattern_match(&pattern, None, None) {
            let momentum = matched[&vk_symbol!("q_")].clone();
            for (symbol, (_, axis)) in get_individual_momenta_from_atom(momentum.as_view())? {
                if symbol != S.k || axis <= 0 {
                    return Ok(None);
                }
                momenta.insert(axis as usize);
            }
            let id = get_integer_from_atom(matched[&vk_symbol!("id_")].as_view())
                .and_then(|id| usize::try_from(id).ok())
                .ok_or_else(|| invalid("noninteger propagator ID"))?;
            propagators.push(Propagator {
                id,
                nodes: (
                    get_integer_from_atom(matched[&vk_symbol!("nl_")].as_view())
                        .filter(|id| *id >= 0)
                        .ok_or_else(|| invalid("invalid source vertex"))?,
                    get_integer_from_atom(matched[&vk_symbol!("nr_")].as_view())
                        .filter(|id| *id >= 0)
                        .ok_or_else(|| invalid("invalid target vertex"))?,
                ),
                momentum,
                mass: matched[&vk_symbol!("m_")].clone(),
                power: matched[&vk_symbol!("power_")].clone(),
            });
        }
        if momenta != (1..=loop_count).collect::<BTreeSet<_>>() || propagators.is_empty() {
            return Ok(None);
        }
        propagators.sort_by_key(|propagator| propagator.id);
        if propagators.windows(2).any(|pair| pair[0].id == pair[1].id) {
            return Err(invalid("duplicate propagator IDs"));
        }
        let mass = &propagators[0].mass;
        if mass.is_zero()
            || propagators
                .iter()
                .any(|propagator| &propagator.mass != mass)
        {
            return Ok(None);
        }
        let mut balances = BTreeMap::<i64, Atom>::new();
        let mut dense = Atom::one();
        let mut original = Atom::one();
        for (index, propagator) in propagators.iter().enumerate() {
            *balances.entry(propagator.nodes.0).or_default() += &propagator.momentum;
            *balances.entry(propagator.nodes.1).or_default() -= &propagator.momentum;
            dense *= S.prop.call((
                index + 1,
                S.edge.call(propagator.nodes),
                &propagator.momentum,
                &propagator.mass,
                &propagator.power,
            ));
            original *= S.prop.call((
                propagator.id,
                S.edge.call(propagator.nodes),
                &propagator.momentum,
                &propagator.mass,
                &propagator.power,
            ));
        }
        if expression != S.topo.call(original).as_view() {
            return Ok(None);
        }
        if balances.values().any(|balance| !balance.expand().is_zero()) {
            return Err(invalid(
                "the input routing violates vacuum momentum conservation",
            ));
        }
        let graph =
            crate::graph::Graph::new_from_atom(S.topo.call(dense).as_view(), propagators.len())?;
        let basis = graph
            .get_one_lmb()?
            .into_iter()
            .map(|edge| edge as usize - 1)
            .collect::<Vec<_>>();
        if basis.len() != loop_count {
            return Err(invalid(
                "graph cycle count differs from its declared loop momenta",
            ));
        }
        let mut unit_graph = graph;
        for edge in unit_graph.edges.values_mut() {
            edge.mass = Atom::one();
        }
        let incidence = unit_graph.to_symbolica_graph(false);
        let canonical_incidence = incidence.canonize();
        let loops = (1..=loop_count)
            .map(|axis| S.k.call(axis))
            .collect::<Vec<_>>();
        let kinematics = Kinematics::new()
            .with_momenta(loops.iter().cloned())
            .map_err(invalid)?;
        let denominators = propagators
            .iter()
            .map(|propagator| {
                kinematics
                    .scalar_product(&propagator.momentum, &propagator.momentum)
                    .map(|square| square - Atom::one())
            })
            .collect::<Result<Vec<_>, _>>()
            .map_err(invalid)?;
        let family =
            IntegralFamily::new(loops, Vec::new(), denominators, &kinematics).map_err(invalid)?;
        log::debug!(target: "vakint::five_loop_match",
            "phase=native_input_authenticated rows={} zero_power_rows={} basis={:?} powers={:?}",
            propagators.len(), propagators.iter().filter(|p|p.power.is_zero()).count(), basis,
            propagators.iter().map(|p|(p.id,&p.power)).collect::<Vec<_>>());
        Ok(Some(Self {
            family,
            propagators,
            basis,
            incidence,
            canonical_incidence,
        }))
    }

    pub fn match_topologies<'a>(
        &self,
        topologies: impl IntoIterator<Item = &'a Topology>,
    ) -> Result<Option<ReplacementRules>, VakintError> {
        let started = Instant::now();
        let topologies = topologies.into_iter().collect::<Vec<_>>();
        for proposal in [Proposal::Incidence, Proposal::Parametric] {
            log::debug!(target: "vakint::five_loop_match", "phase=proposal_begin proposal={proposal:?}");
            for topology in &topologies {
                if let Some(matched) = self.match_topology(topology, proposal)? {
                    log::debug!(target: "vakint::five_loop_match", "phase=matched proposal={proposal:?} elapsed_seconds={:.6}", started.elapsed().as_secs_f64());
                    return Ok(Some(matched));
                }
            }
        }
        Ok(None)
    }

    fn match_topology(
        &self,
        topology: &Topology,
        proposal: Proposal,
    ) -> Result<Option<ReplacementRules>, VakintError> {
        let integral = topology.get_integral();
        if integral.n_loops != self.family.loop_momenta().len() {
            return Ok(None);
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
            return Ok(None);
        }
        let started = Instant::now();
        let incidence = match proposal {
            Proposal::Incidence => {
                let Some(permutation) = self.incidence_permutation(&integral.graph) else {
                    return Ok(None);
                };
                Some(permutation)
            }
            Proposal::Parametric => None,
        };
        log::debug!(target: "vakint::five_loop_match", "phase=proposal_candidate proposal={proposal:?} topology={} elapsed_seconds={:.6}", integral.name, started.elapsed().as_secs_f64());
        let kinematics = self.family.kinematics();
        let denominators = physical
            .iter()
            .map(|(_, momentum)| {
                kinematics
                    .scalar_product(momentum, momentum)
                    .map(|square| square - Atom::one())
            })
            .collect::<Result<Vec<_>, _>>()
            .map_err(invalid)?;
        let target = IntegralFamily::new(
            self.family.loop_momenta().to_vec(),
            Vec::new(),
            denominators,
            kinematics,
        )
        .map_err(invalid)?;
        log::debug!(target: "vakint::five_loop_match", "phase=target_family proposal={proposal:?} topology={} elapsed_seconds={:.6}", integral.name, started.elapsed().as_secs_f64());
        let permutation = if let Some(permutation) = incidence {
            permutation
        } else {
            let parameters = (0..physical.len())
                .map(|axis| vk_symbol!("momentum_match_parameter").call(axis))
                .collect::<Vec<_>>();
            let Some(permutation) = self
                .family
                .parametric_mapping(&target, &parameters)
                .map_err(invalid)?
            else {
                return Ok(None);
            };
            permutation.denominator_map().to_vec()
        };
        log::debug!(target: "vakint::five_loop_match", "phase=permutation proposal={proposal:?} topology={} elapsed_seconds={:.6}", integral.name, started.elapsed().as_secs_f64());
        // The native parameter map proposes line images, not an integral map.
        // Only the orientations of an independent basis remain undetermined.
        // Reuse Vakint's existing native exact basis solve, then ask FeynKit to
        // verify every denominator and the Jacobian for the full candidate.
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
            .map(|slot| {
                self.propagators[*slot]
                    .momentum
                    .replace_multiple(&substitute)
            })
            .collect::<Vec<_>>();
        let images = self
            .basis
            .iter()
            .map(|slot| physical[permutation[*slot]].1.clone())
            .collect::<Vec<_>>();
        let attempts = 1usize
            .checked_shl(loops.len().try_into().map_err(invalid)?)
            .filter(|attempts| *attempts <= 1_000_000)
            .ok_or_else(|| invalid("momentum-map orientation budget exceeded"))?;
        let mut verified = None;
        for signs in 0..attempts {
            let system = sources
                .iter()
                .zip(&images)
                .enumerate()
                .map(|(axis, (source, image))| {
                    source - Atom::num(if signs & (1 << axis) == 0 { 1 } else { -1 }) * image
                })
                .collect::<Vec<_>>();
            let coordinates = Topology::solve_loop_coordinates(&system, &variables)?;
            log::debug!(target: "vakint::five_loop_match", "phase=solved signs={signs} topology={} elapsed_seconds={:.6}", integral.name, started.elapsed().as_secs_f64());
            if let Some(mapping) = self
                .family
                .mapping_to(&target, &coordinates)
                .map_err(invalid)?
            {
                verified = Some(mapping);
                log::debug!(target: "vakint::five_loop_match", "phase=verified signs={signs} topology={} elapsed_seconds={:.6}", integral.name, started.elapsed().as_secs_f64());
                break;
            }
        }
        let Some(mapping) = verified else {
            if matches!(proposal, Proposal::Incidence) {
                log::debug!(target: "vakint::five_loop_match",
                    "incidence proposal has no verified unit-Jacobian momentum map: target={}, basis={:?}, source={:?}, target_momenta={:?}, permutation={permutation:?}",
                    integral.name,
                    self.basis,
                    self.propagators
                        .iter()
                        .map(|propagator| &propagator.momentum)
                        .collect::<Vec<_>>(),
                    physical
                        .iter()
                        .map(|(_, momentum)| momentum)
                        .collect::<Vec<_>>()
                );
            }
            return Ok(None);
        };
        Ok(Some(self.replacement(topology, &physical, &mapping)))
    }

    fn replacement(
        &self,
        topology: &Topology,
        physical: &[(usize, Atom)],
        mapping: &feynkit_graph::IntegralMapping,
    ) -> ReplacementRules {
        let integral = topology.get_integral();
        let mut replacement = ReplacementRules {
            canonical_topology: topology.clone(),
            ..Default::default()
        };
        for slot in 1..=integral.n_props {
            replacement
                .canonical_expression_substitutions
                .insert(S.pow.call(slot), Atom::Zero);
        }
        replacement
            .canonical_expression_substitutions
            .insert(vk_symbol!("msq").call(1), self.propagators[0].mass.clone());
        for (source, target_slot) in mapping.denominator_map().iter().enumerate() {
            let target_id = physical[*target_slot].0;
            replacement.canonical_expression_substitutions.insert(
                S.pow.call(target_id),
                self.propagators[source].power.clone(),
            );
            replacement
                .edge_ids_canonical_to_input_map
                .insert(target_id, self.propagators[source].id);
        }
        for (axis, (_, image)) in mapping.momentum_rules().iter().enumerate() {
            let image = image
                .replace(vk_parse!("k(axis_)").unwrap().to_pattern())
                .allow_new_wildcards_on_rhs(true)
                .with(vk_parse!("k(axis_,idx___)").unwrap().to_pattern());
            replacement
                .numerator_substitutions
                .insert(S.k.call((axis + 1, Atom::var(vk_symbol!("idx___")))), image);
        }
        replacement
    }

    fn incidence_permutation(&self, graph: &crate::graph::Graph) -> Option<Vec<usize>> {
        let mut unit_graph = graph.clone();
        for edge in unit_graph.edges.values_mut() {
            edge.mass = Atom::one();
        }
        let target = unit_graph.to_symbolica_graph(false);
        let canonical = target.canonize();
        if self.canonical_incidence.graph != canonical.graph {
            return None;
        }
        let endpoints = |vertices: (usize, usize), mapping: &[usize]| {
            let (left, right) = (mapping[vertices.0], mapping[vertices.1]);
            (left.min(right), left.max(right))
        };
        let mut targets = BTreeMap::<(usize, usize), std::collections::VecDeque<usize>>::new();
        for (index, edge) in target.edges().iter().enumerate() {
            targets
                .entry(endpoints(edge.vertices, &canonical.vertex_map))
                .or_default()
                .push_back(index);
        }
        self.incidence
            .edges()
            .iter()
            .map(|edge| {
                targets
                    .get_mut(&endpoints(
                        edge.vertices,
                        &self.canonical_incidence.vertex_map,
                    ))?
                    .pop_front()
            })
            .collect()
    }
}
