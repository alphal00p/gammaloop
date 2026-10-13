//! Model-local interaction stars: ordered particle ports and their UFO-slot maps.
//!
//! A plan refers only to positions in a vertex's ordered incidence. Each use
//! binds the current edge/half-edge identities and flow; graph frames and
//! symbolic expression ownership remain with their existing owners.
use std::{
    collections::HashMap,
    sync::{OnceLock, RwLock},
};

use super::{GeneratedGraph, GenerationError, NumeratorHalfEdge, ParticleStub, incoming_stub};
use feynkit_model::{Model, ParticleId, VertexRule, VertexRuleId};
use linnet::half_edge::involution::Flow;

pub(super) struct InteractionStars {
    incoming: Vec<ParticleStub>,
    rules: Vec<RulePorts>,
}

struct Port {
    particle: ParticleId,
    stub: ParticleStub,
    spin: i64,
    color: i64,
}

pub(super) struct RulePorts {
    ports: Vec<Port>,
    signature: Vec<ParticleStub>,
    conjugate: Vec<ParticleStub>,
    cp_candidates: OnceLock<Vec<VertexRuleId>>,
    slots: RwLock<HashMap<Vec<ParticleStub>, Box<[usize]>>>,
}

/// Scratch incidence, reused across a diagram's vertices. Both ends of a loop
/// retain their separate positions even when they carry the same particle.
#[derive(Default)]
pub(super) struct InteractionStar {
    legs: Vec<NumeratorHalfEdge>,
    stubs: Vec<ParticleStub>,
}

impl InteractionStars {
    pub(super) fn new(model: &Model) -> Self {
        // Model construction has already resolved every particle reference.
        // This prepares metadata only; rule expressions are still lowered on demand.
        let incoming: Vec<_> = model
            .particles()
            .iter()
            .map(|p| incoming_stub(model, p))
            .collect();
        let rules = model
            .vertex_rules()
            .iter()
            .map(|rule| {
                let ports: Vec<_> = rule
                    .particles
                    .iter()
                    .map(|&id| {
                        let particle = &model.particles()[id.index()];
                        Port {
                            particle: id,
                            stub: incoming[id.index()].clone(),
                            spin: particle.spin,
                            color: particle.color,
                        }
                    })
                    .collect();
                let mut signature: Vec<_> = ports.iter().map(|p| p.stub.clone()).collect();
                signature.sort();
                let mut conjugate: Vec<_> = signature
                    .iter()
                    .cloned()
                    .map(ParticleStub::flipped)
                    .collect();
                conjugate.sort();
                RulePorts {
                    ports,
                    signature,
                    conjugate,
                    cp_candidates: OnceLock::new(),
                    slots: RwLock::default(),
                }
            })
            .collect();
        Self { incoming, rules }
    }

    pub(super) fn rule(&self, id: VertexRuleId) -> &RulePorts {
        &self.rules[id.index()]
    }

    pub(super) fn read(
        &self,
        graph: &GeneratedGraph,
        vertex: usize,
        model: &Model,
        star: &mut InteractionStar,
    ) -> Result<(), GenerationError> {
        star.legs.clear();
        star.stubs.clear();
        for (edge_id, flow) in graph.incident(vertex) {
            let edge = graph.edge(edge_id);
            let particle = model.particle_by_id(edge.particle)?;
            let entering = if flow == Flow::Sink || !edge.directed {
                edge.particle
            } else {
                particle.antiparticle
            };
            star.stubs.push(self.incoming[entering.index()].clone());
            star.legs.push(NumeratorHalfEdge {
                edge: edge_id,
                flow,
                hedge: 2 * edge_id + usize::from(flow == Flow::Sink),
                spin: 0,
                color: 1,
            });
        }
        Ok(())
    }

    pub(super) fn repair_rule(
        &self,
        rule_id: VertexRuleId,
        vertex: usize,
        star: &InteractionStar,
        model: &Model,
    ) -> Result<VertexRuleId, GenerationError> {
        let rule = model.vertex_rule_by_id(rule_id)?;
        let ports = self.rule(rule_id);
        if star.stubs.len() == ports.ports.len() && ports.with_slots(&star.stubs, |_| ()).is_ok() {
            return Ok(rule_id);
        }
        let mut actual = star.stubs.clone();
        actual.sort();
        let candidates = if actual == ports.conjugate {
            ports
                .cp_candidates
                .get_or_init(|| {
                    let orders = rule.coupling_orders(model);
                    model
                        .vertex_rules()
                        .iter()
                        .enumerate()
                        .filter_map(|(i, candidate)| {
                            (self.rules[i].signature == ports.conjugate
                                && candidate.coupling_orders(model) == orders)
                                .then_some(VertexRuleId::from_index(i))
                        })
                        .collect()
                })
                .as_slice()
        } else {
            &[]
        };
        match candidates {
            [candidate] => Ok(*candidate),
            _ => Err(GenerationError::InvalidCpInteractionRule {
                vertex,
                interaction: rule.name.clone(),
                candidates: candidates
                    .iter()
                    .map(|id| model.vertex_rules()[id.index()].name.clone())
                    .collect(),
            }),
        }
    }
}

impl RulePorts {
    /// Borrow a successful port permutation on a hit; allocate its key only on a
    /// miss. Failed matches retain no cached owner-dependent diagnostic.
    fn with_slots<R>(
        &self,
        stubs: &[ParticleStub],
        apply: impl FnOnce(&[usize]) -> R,
    ) -> Result<R, usize> {
        if let Some(slots) = self.slots.read().unwrap().get(stubs) {
            return Ok(apply(slots));
        }
        let mut slots = Vec::with_capacity(self.ports.len());
        for (slot, port) in self.ports.iter().enumerate() {
            let position = stubs
                .iter()
                .enumerate()
                .position(|(position, stub)| *stub == port.stub && !slots.contains(&position))
                .ok_or(slot)?;
            slots.push(position);
        }
        let result = apply(&slots);
        self.slots
            .write()
            .unwrap()
            .insert(stubs.to_vec(), slots.into_boxed_slice());
        Ok(result)
    }

    pub(super) fn bind(
        &self,
        star: &InteractionStar,
        vertex: usize,
        rule: &VertexRule,
        model: &Model,
        mut apply: impl FnMut(usize, NumeratorHalfEdge),
    ) -> Result<(), GenerationError> {
        self.with_slots(&star.stubs, |slots| {
            for (slot, &position) in slots.iter().enumerate() {
                let mut leg = star.legs[position];
                leg.spin = self.ports[slot].spin;
                leg.color = self.ports[slot].color;
                apply(slot, leg);
            }
        })
        .map_err(|slot| GenerationError::MissingInteractionLeg {
            vertex,
            interaction: rule.name.clone(),
            leg: slot + 1,
            particle: model.particles()[self.ports[slot].particle.index()]
                .name
                .clone(),
        })
    }

    pub(super) fn spin_color(&self, slot: usize) -> (i64, i64) {
        (self.ports[slot].spin, self.ports[slot].color)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::generation::graph::{Endpoint, ExternalNode, GeneratedEdge};
    use feynkit_graph::ExternalState;

    #[test]
    fn star_plan_rebinds_global_ids_and_distinguishes_self_loop_ports() {
        let model = Model::phi3();
        let rule = VertexRuleId::from_index(0);
        let particle = model.vertex_rules()[0].particles[0];
        let plans = InteractionStars::new(&model);
        let external = |index| {
            Endpoint::External(ExternalNode {
                index,
                state: ExternalState::Incoming,
                particle,
            })
        };
        let edge = |source, sink| GeneratedEdge {
            source,
            sink,
            particle,
            directed: false,
        };
        let first = GeneratedGraph::new(
            [rule],
            [
                edge(external(0), Endpoint::Vertex(0)),
                edge(Endpoint::Vertex(0), Endpoint::Vertex(0)),
            ],
        );
        let second = GeneratedGraph::new(
            [rule, rule],
            [
                edge(external(0), Endpoint::Vertex(0)),
                edge(Endpoint::Vertex(1), Endpoint::Vertex(1)),
                edge(external(1), Endpoint::Vertex(1)),
            ],
        );
        let mut star = InteractionStar::default();
        let mut bound = Vec::new();
        for (graph, vertex) in [(&first, 0), (&second, 1)] {
            plans.read(graph, vertex, &model, &mut star).unwrap();
            let mut legs = Vec::new();
            plans
                .rule(rule)
                .bind(&star, vertex, &model.vertex_rules()[0], &model, |_, leg| {
                    legs.push((leg.edge, leg.hedge, leg.flow, leg.spin, leg.color))
                })
                .unwrap();
            bound.push(legs);
        }
        assert_eq!(plans.rule(rule).slots.read().unwrap().len(), 1);
        assert_eq!(
            bound[0],
            [
                (0, 1, Flow::Sink, 1, 1),
                (1, 2, Flow::Source, 1, 1),
                (1, 3, Flow::Sink, 1, 1)
            ]
        );
        assert_eq!(
            bound[1],
            [
                (1, 2, Flow::Source, 1, 1),
                (1, 3, Flow::Sink, 1, 1),
                (2, 5, Flow::Sink, 1, 1)
            ]
        );
    }

    #[test]
    fn star_plan_orders_charged_ports_and_caches_only_successes() {
        let model = Model::qed();
        let rule_id = model.vertex_rule_id("V_98").unwrap();
        let rule = model.vertex_rule_by_id(rule_id).unwrap();
        let plans = InteractionStars::new(&model);
        let edges = [2, 0, 1].into_iter().enumerate().map(|(edge, slot)| {
            let particle = rule.particles[slot];
            GeneratedEdge {
                source: Endpoint::External(ExternalNode {
                    index: edge,
                    state: ExternalState::Incoming,
                    particle,
                }),
                sink: Endpoint::Vertex(0),
                particle,
                directed: model.particle_by_id(particle).unwrap().antiparticle != particle,
            }
        });
        let graph = GeneratedGraph::new([rule_id], edges);
        let mut star = InteractionStar::default();
        plans.read(&graph, 0, &model, &mut star).unwrap();
        assert_eq!(
            plans.repair_rule(rule_id, 0, &star, &model).unwrap(),
            rule_id
        );
        let mut legs = Vec::new();
        plans
            .rule(rule_id)
            .bind(&star, 7, rule, &model, |_, leg| legs.push(leg))
            .unwrap();
        assert_eq!(
            legs.iter().map(|leg| leg.edge).collect::<Vec<_>>(),
            [1, 2, 0]
        );
        for (leg, &particle) in legs.iter().zip(&rule.particles) {
            let particle = model.particle_by_id(particle).unwrap();
            assert_eq!((leg.spin, leg.color), (particle.spin, particle.color));
        }
        assert_eq!(plans.rule(rule_id).slots.read().unwrap().len(), 1);
        star.stubs.clear();
        for vertex in [7, 19] {
            assert!(
                matches!(plans.rule(rule_id).bind(&star, vertex, rule, &model, |_, _| {}),
                Err(GenerationError::MissingInteractionLeg { vertex: actual, leg: 1, .. }) if actual == vertex)
            );
        }
        assert_eq!(plans.rule(rule_id).slots.read().unwrap().len(), 1);
    }

    #[test]
    fn star_cp_candidates_preserve_coupling_uniqueness_and_error_owner() {
        let model = Model::standard_model();
        let rule = model.vertex_rule_id("V_82").unwrap();
        let conjugate = model.vertex_rule_id("V_117").unwrap();
        let plans = InteractionStars::new(&model);
        let star = InteractionStar {
            stubs: plans.rule(rule).conjugate.clone(),
            ..Default::default()
        };
        assert_eq!(
            plans.repair_rule(rule, 3, &star, &model).unwrap(),
            conjugate
        );
        assert_eq!(
            plans.repair_rule(rule, 4, &star, &model).unwrap(),
            conjugate
        );

        let mut definition: serde_json::Value =
            serde_json::from_str(&model.to_json().unwrap()).unwrap();
        let rules = definition["vertex_rules"].as_array_mut().unwrap();
        let mut duplicate = rules
            .iter()
            .find(|rule| rule["name"] == "V_117")
            .unwrap()
            .clone();
        duplicate["name"] = serde_json::json!("V_117_duplicate");
        rules.push(duplicate);
        let ambiguous = Model::from_json(&definition.to_string()).unwrap();
        let rule = ambiguous.vertex_rule_id("V_82").unwrap();
        let plans = InteractionStars::new(&ambiguous);
        for vertex in [3, 4] {
            assert!(matches!(plans.repair_rule(rule, vertex, &star, &ambiguous),
                Err(GenerationError::InvalidCpInteractionRule { vertex: actual, candidates, .. })
                if actual == vertex && candidates == ["V_117", "V_117_duplicate"]));
        }
    }
}
