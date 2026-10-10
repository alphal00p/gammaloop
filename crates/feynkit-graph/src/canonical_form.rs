//! Canonical forms of diagrams.
//!
//! A diagram's canonical form is Linnet's canonical form of a keyed graph: the diagram graph with
//! vertices keyed by their interaction, edges by their particle, external leg (without its name)
//! and dummy flag, and half-edges by their interaction slot, plus a pair of side vertices per cut
//! and threshold candidate, joined to each other and to the diagram vertices they contain. A
//! directed edge points from its source to its target vertex; stored orientations are not part of
//! the key. The structural form keeps particle orientations as stored and identifies the diagram
//! ([`DiagramId`]). The normalized form ([`FeynmanDiagram::canonical_key`]) turns every
//! antiparticle edge into a particle edge against its direction, so that it identifies the
//! physical diagram whatever edge representation it uses.

use std::collections::BTreeMap;

use feynkit_model::{Model, ParticleId, VertexRuleId};
use linnet::half_edge::{
    HedgeGraph, HedgeGraphError, NodeIndex,
    canon::CanonicalForm,
    involution::{EdgeData, Flow, Orientation},
};

use crate::{
    CanonicalCutSide, CanonicalDiagramEdge, CanonicalDiagramKey, CanonicalDiagramVertex,
    DiagramEdge, DiagramError, DiagramHalfEdge, DiagramId, DiagramVertex, ExternalLeg,
    ExternalState, FeynmanDiagram, VertexSlot,
};

/// Vertex key of the keyed graph, ordered like the [`CanonicalDiagramVertex`] it stands for.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum VertexKey<'a> {
    Diagram(Option<VertexRuleId>),
    /// A cut side with its coupling orders and loop count.
    CutSide(CanonicalCutSide, &'a BTreeMap<String, usize>, usize),
    ThresholdSide(CanonicalCutSide),
}

/// Edge key of the keyed graph, ordered like the [`CanonicalDiagramEdge`] it stands for. A leg is
/// its index, state and connection.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum EdgeKey {
    Diagram {
        particle: ParticleId,
        external: Option<(usize, ExternalState, usize)>,
        is_dummy: bool,
    },
    CutPair,
    CutMembership,
}

impl DiagramVertex {
    fn structural_key(&self) -> VertexKey<'static> {
        VertexKey::Diagram(self.interaction)
    }

    pub(super) fn same_structure(&self, other: &Self) -> bool {
        self.structural_key() == other.structural_key()
    }
}

impl DiagramEdge {
    fn structural_key(&self, particle: ParticleId) -> EdgeKey {
        EdgeKey::Diagram {
            particle,
            external: self
                .external
                .as_ref()
                .map(|leg| (leg.index, leg.state, leg.connection)),
            is_dummy: self.is_dummy,
        }
    }

    pub(super) fn same_structure(&self, other: &Self) -> bool {
        self.structural_key(self.particle) == other.structural_key(other.particle)
            && self.directed == other.directed
            && self.source_slot == other.source_slot
            && self.target_slot == other.target_slot
    }
}

type Form<'a> = CanonicalForm<VertexKey<'a>, EdgeKey, Option<VertexSlot>>;

impl FeynmanDiagram {
    /// Build a deterministic key for the model-resolved colored topology.
    ///
    /// Canonical labeling is delegated to Linnet's half-edge canonizer.
    /// Particle/antiparticle edge representations are normalized so reversing
    /// a directed edge with [`Self::reverse_edge`] does not change the key.
    pub fn canonical_key(&self) -> Result<CanonicalDiagramKey, DiagramError> {
        self.validate()?;
        let (words, vertices, edges, half_edges) = self.canonical_form(true)?.into_parts();
        Ok(CanonicalDiagramKey {
            words,
            vertices: vertices.into_iter().map(Into::into).collect(),
            edges: edges.into_iter().map(Into::into).collect(),
            half_edges,
        })
    }

    /// The id of the structural canonical form.
    pub(crate) fn structural_id(&self) -> Result<DiagramId, DiagramError> {
        DiagramId::of_form(&self.model, &self.canonical_form(false)?)
    }

    /// The canonical form of the keyed graph, with antiparticle edges turned into particle edges
    /// against their direction if `normalize`.
    fn canonical_form(&self, normalize: bool) -> Result<Form<'_>, DiagramError> {
        // If normalizing, per edge in index order: the particle of its key and whether it points
        // against the edge.
        let normalized = normalize
            .then(|| {
                self.edges()
                    .map(|(edge_id, _, edge)| {
                        let (resolved, base) = self.resolve_particle(edge_id, edge)?;
                        Ok((base, edge.directed && resolved.is_antiparticle()))
                    })
                    .collect::<Result<Vec<_>, DiagramError>>()
            })
            .transpose()?;
        let mut graph = self.graph.map_data_ref(
            |_, _, vertex| vertex.structural_key(),
            |_, index, _, edge| {
                let edge = edge.data;
                let (particle, reversed) = match &normalized {
                    Some(normalized) => normalized[index.0],
                    None => (edge.particle, false),
                };
                EdgeData {
                    orientation: match (edge.directed, reversed) {
                        (false, _) => Orientation::Undirected,
                        (true, false) => Orientation::Default,
                        (true, true) => Orientation::Reversed,
                    },
                    data: edge.structural_key(particle),
                }
            },
            |hedge, _| {
                let edge = &self.graph[self.graph[&hedge]];
                Some(match self.graph.flow(hedge) {
                    Flow::Source => edge.source_slot,
                    Flow::Sink => edge.target_slot,
                })
            },
        );
        for cut in &self.cuts {
            let sides = [
                (CanonicalCutSide::Left, &cut.left),
                (CanonicalCutSide::Right, &cut.right),
            ]
            .map(|(side, summary)| {
                VertexKey::CutSide(side, &summary.coupling_orders, summary.loop_count)
            });
            self.add_sides(
                &mut graph,
                sides,
                [&cut.left.half_edges, &cut.right.half_edges],
            )?;
        }
        for candidate in &self.topology_threshold_candidates {
            let sides =
                [CanonicalCutSide::Left, CanonicalCutSide::Right].map(VertexKey::ThresholdSide);
            self.add_sides(&mut graph, sides, [&candidate.left, &candidate.right])?;
        }
        Ok(graph.canonical_form_by(|&vertex| vertex, |&edge| edge, |&slot| slot))
    }

    /// Add a pair of side vertices to the keyed graph, joined to each other and to the diagram
    /// vertices of their half-edges.
    fn add_sides<'a>(
        &self,
        graph: &mut HedgeGraph<EdgeKey, VertexKey<'a>, Option<VertexSlot>>,
        sides: [VertexKey<'a>; 2],
        half_edges: [&[DiagramHalfEdge]; 2],
    ) -> Result<(), DiagramError> {
        let invariant = |error: HedgeGraphError| DiagramError::Invariant {
            operation: "canonicalizing the diagram",
            message: error.to_string(),
        };
        let sides = sides.map(|side| graph.add_node(side));
        graph
            .add_pair_mut(
                sides[0],
                sides[1],
                EdgeKey::CutPair,
                Orientation::Undirected,
            )
            .map_err(invariant)?;
        for (side, half_edges) in sides.into_iter().zip(half_edges) {
            let mut members: Vec<NodeIndex> = half_edges
                .iter()
                .filter_map(|&half_edge| self.half_edge_id(half_edge))
                .map(|hedge| self.graph.node_id(hedge))
                .collect();
            members.sort_unstable();
            members.dedup();
            for vertex in members {
                graph
                    .add_pair_mut(
                        side,
                        vertex,
                        EdgeKey::CutMembership,
                        Orientation::Undirected,
                    )
                    .map_err(invariant)?;
            }
        }
        Ok(())
    }
}

impl From<VertexKey<'_>> for CanonicalDiagramVertex {
    fn from(key: VertexKey<'_>) -> Self {
        match key {
            VertexKey::Diagram(interaction) => Self::Diagram { interaction },
            VertexKey::CutSide(side, coupling_orders, loop_count) => Self::CutSide {
                side,
                coupling_orders: coupling_orders.clone(),
                loop_count,
            },
            VertexKey::ThresholdSide(side) => Self::TopologyThresholdSide { side },
        }
    }
}

impl From<EdgeKey> for CanonicalDiagramEdge {
    fn from(key: EdgeKey) -> Self {
        match key {
            EdgeKey::Diagram {
                particle,
                external,
                is_dummy,
            } => Self::Diagram {
                particle,
                external: external.map(|(index, state, connection)| ExternalLeg {
                    name: String::new(),
                    index,
                    state,
                    connection,
                }),
                is_dummy,
            },
            EdgeKey::CutPair => Self::CutPair,
            EdgeKey::CutMembership => Self::CutMembership,
        }
    }
}

impl VertexKey<'_> {
    /// Append the id bytes of this key: its variant tag, then its fields.
    fn encode(self, model: &Model, bytes: &mut IdBytes) -> Result<(), DiagramError> {
        match self {
            VertexKey::Diagram(interaction) => {
                bytes.tag(0);
                match interaction {
                    None => bytes.tag(0),
                    Some(rule) => {
                        bytes.tag(1);
                        bytes.uint(rule.index());
                        bytes.text(&model.vertex_rule_by_id(rule)?.name);
                    }
                }
            }
            VertexKey::CutSide(side, coupling_orders, loop_count) => {
                bytes.tag(1);
                bytes.tag(side as u8);
                bytes.uint(coupling_orders.len());
                for (order, &power) in coupling_orders {
                    bytes.text(order);
                    bytes.uint(power);
                }
                bytes.uint(loop_count);
            }
            VertexKey::ThresholdSide(side) => {
                bytes.tag(2);
                bytes.tag(side as u8);
            }
        }
        Ok(())
    }
}

impl EdgeKey {
    /// Append the id bytes of this key: its variant tag, then its fields.
    fn encode(self, model: &Model, bytes: &mut IdBytes) -> Result<(), DiagramError> {
        match self {
            EdgeKey::Diagram {
                particle,
                external,
                is_dummy,
            } => {
                bytes.tag(0);
                bytes.uint(particle.index());
                bytes.text(&model.particle_by_id(particle)?.name);
                match external {
                    None => bytes.tag(0),
                    Some((index, state, connection)) => {
                        bytes.tag(1);
                        bytes.uint(index);
                        bytes.tag(state as u8);
                        bytes.uint(connection);
                    }
                }
                bytes.tag(is_dummy.into());
            }
            EdgeKey::CutPair => bytes.tag(1),
            EdgeKey::CutMembership => bytes.tag(2),
        }
        Ok(())
    }
}

impl DiagramId {
    /// The first 128 bits, little endian, of the BLAKE3 hash of the bytes of `form`: its words,
    /// then its node, edge and half-edge key tables, each after its length.
    fn of_form(model: &Model, form: &Form<'_>) -> Result<Self, DiagramError> {
        let mut bytes = IdBytes(Vec::with_capacity(256));
        bytes.uint(form.words().len());
        for &word in form.words() {
            bytes.uint(word as usize);
        }
        bytes.uint(form.node_keys().len());
        for key in form.node_keys() {
            key.encode(model, &mut bytes)?;
        }
        bytes.uint(form.edge_keys().len());
        for key in form.edge_keys() {
            key.encode(model, &mut bytes)?;
        }
        bytes.uint(form.hedge_keys().len());
        for slot in form.hedge_keys() {
            match slot {
                None => bytes.tag(0),
                Some(slot) => {
                    bytes.tag(1);
                    bytes.uint(slot.0);
                }
            }
        }
        let hash = blake3::hash(&bytes.0);
        let (id, _) = hash.as_bytes().split_first_chunk().expect("a 32-byte hash");
        Ok(Self(u128::from_le_bytes(*id)))
    }
}

/// The bytes a [`DiagramId`] hashes: unsigned LEB128 integers, one byte per enum tag, and
/// sequences and strings after their length. The encoding is prefix free, so distinct canonical
/// forms have distinct bytes.
struct IdBytes(Vec<u8>);

impl IdBytes {
    fn tag(&mut self, tag: u8) {
        self.0.push(tag);
    }

    fn uint(&mut self, mut value: usize) {
        while value >= 0x80 {
            self.0.push(value as u8 | 0x80);
            value >>= 7;
        }
        self.0.push(value as u8);
    }

    fn text(&mut self, text: &str) {
        self.uint(text.len());
        self.0.extend_from_slice(text.as_bytes());
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use feynkit_model::{ComplexValue, ParameterCard};
    use symbolica::atom::Atom;

    use super::*;
    use crate::{
        DiagramEdge, DiagramVertex, EdgeId, VertexId,
        tests::{
            cut_scalar_bubble, directed_fermion_line, fermion_model, one_loop, one_loop_in,
            one_loop_reordered, scalar_model,
        },
    };

    fn encoded(write: impl FnOnce(&mut IdBytes)) -> Vec<u8> {
        let mut bytes = IdBytes(Vec::new());
        write(&mut bytes);
        bytes.0
    }

    #[test]
    fn id_bytes_are_leb128_integers_and_length_prefixed_strings() {
        assert_eq!(encoded(|b| b.uint(0)), [0]);
        assert_eq!(encoded(|b| b.uint(127)), [0x7f]);
        assert_eq!(encoded(|b| b.uint(128)), [0x80, 0x01]);
        assert_eq!(encoded(|b| b.uint(300)), [0xac, 0x02]);
        assert_eq!(
            encoded(|b| b.uint(u32::MAX as usize)),
            [0xff, 0xff, 0xff, 0xff, 0x0f]
        );
        assert_eq!(encoded(|b| b.text("V_f")), [3, b'V', b'_', b'f']);
        assert_eq!(encoded(|b| b.text("")), [0]);
    }

    /// Every key variant encodes as its tag, then its fields; particles and rules as model index
    /// and name, which must exist.
    #[test]
    fn keys_encode_their_variant_tag_then_their_fields() {
        let model = fermion_model();
        let vertex = |key: VertexKey| encoded(|b| key.encode(&model, b).unwrap());
        let edge = |key: EdgeKey| encoded(|b| key.encode(&model, b).unwrap());
        let [f, anti_f] = ["f", "f~"].map(|name| model.particle_id(name).unwrap());

        assert_eq!(vertex(VertexKey::Diagram(None)), [0, 0]);
        assert_eq!(
            vertex(VertexKey::Diagram(Some(
                model.vertex_rule_id("V_af").unwrap()
            ))),
            [0, 1, 1, 4, b'V', b'_', b'a', b'f']
        );
        let orders = BTreeMap::from([("QCD".to_owned(), 2), ("QED".to_owned(), 130)]);
        assert_eq!(
            vertex(VertexKey::CutSide(CanonicalCutSide::Right, &orders, 1)),
            [
                1, 1, 2, 3, b'Q', b'C', b'D', 2, 3, b'Q', b'E', b'D', 0x82, 0x01, 1
            ]
        );
        assert_eq!(
            vertex(VertexKey::ThresholdSide(CanonicalCutSide::Left)),
            [2, 0]
        );

        let internal = EdgeKey::Diagram {
            particle: f,
            external: None,
            is_dummy: false,
        };
        assert_eq!(edge(internal), [0, 0, 1, b'f', 0, 0]);
        let dummy_leg = EdgeKey::Diagram {
            particle: anti_f,
            external: Some((2, ExternalState::Outgoing, 3)),
            is_dummy: true,
        };
        assert_eq!(edge(dummy_leg), [0, 1, 2, b'f', b'~', 1, 2, 1, 3, 1]);
        assert_eq!(edge(EdgeKey::CutPair), [1]);
        assert_eq!(edge(EdgeKey::CutMembership), [2]);

        let mut bytes = IdBytes(Vec::new());
        assert!(matches!(
            VertexKey::Diagram(Some(VertexRuleId::from_index(2))).encode(&model, &mut bytes),
            Err(DiagramError::Model(_))
        ));
        assert!(matches!(
            EdgeKey::Diagram {
                particle: ParticleId::from_index(2),
                external: None,
                is_dummy: false,
            }
            .encode(&model, &mut bytes),
            Err(DiagramError::Model(_))
        ));
    }

    /// The id hashes the documented encoding of the structural canonical form: the words, then
    /// the node, edge and half-edge key tables (sorted distinct keys), spelled out for
    /// - a fermion line: the rules `V_f` (index 0) and `V_af` (index 1), one directed `f` edge
    ///   (index 0) and the slot 0 at both ends;
    /// - the cut scalar bubble: the rule `V_1` (index 0), cut sides without coupling orders or
    ///   loops and threshold sides; an internal and an incoming `phi` edge (index 0), the cut
    ///   pair and membership edges; half-edges without a slot (at the sides) and with slots 0
    ///   to 2.
    ///
    /// The ids are pinned, so that a change of the canonizer or of the encoding, which changes
    /// every persisted id, is seen.
    #[test]
    fn diagram_ids_hash_the_encoded_structural_canonical_form() {
        let fermion_line: [&[u8]; 3] = [
            &[
                2, 0, 1, 0, 3, b'V', b'_', b'f', 0, 1, 1, 4, b'V', b'_', b'a', b'f',
            ],
            &[1, 0, 0, 1, b'f', 0, 0],
            &[1, 1, 0],
        ];
        let cut_bubble: [&[u8]; 3] = [
            &[
                5, 0, 1, 0, 3, b'V', b'_', b'1', 1, 0, 0, 0, 1, 1, 0, 0, 2, 0, 2, 1,
            ],
            &[
                4, 0, 0, 3, b'p', b'h', b'i', 0, 0, 0, 0, 3, b'p', b'h', b'i', 1, 0, 0, 0, 0, 1, 2,
            ],
            &[4, 0, 1, 0, 1, 1, 1, 2],
        ];
        for (diagram, tables, id) in [
            (
                directed_fermion_line(),
                fermion_line,
                "b44a9fc65c1401b45c6d1e51529823ff",
            ),
            (
                cut_scalar_bubble(),
                cut_bubble,
                "5cf5e5c03bb0f52a99ac16a5f34f653c",
            ),
        ] {
            let form = diagram.canonical_form(false).unwrap();
            let mut expected = encoded(|b| {
                b.uint(form.words().len());
                for &word in form.words() {
                    b.uint(word as usize);
                }
            });
            expected.extend(tables.concat());
            let hash = blake3::hash(&expected);
            let hashed = u128::from_le_bytes(hash.as_bytes()[..16].try_into().unwrap());
            assert_eq!(
                diagram.id().unwrap(),
                DiagramId(hashed),
                "{:?}",
                diagram.name()
            );
            assert_eq!(
                diagram.id().unwrap().to_string(),
                id,
                "{:?}",
                diagram.name()
            );
        }
    }

    /// The keys order like the public key types they convert to, so the canonical key is the
    /// canonical form of the keyed graph under the public keys.
    #[test]
    fn keys_order_like_the_canonical_key_types() {
        let (none, qcd) = (BTreeMap::new(), BTreeMap::from([("QCD".to_owned(), 2)]));
        let mut vertices = vec![
            VertexKey::Diagram(None),
            VertexKey::Diagram(Some(VertexRuleId::from_index(0))),
            VertexKey::Diagram(Some(VertexRuleId::from_index(1))),
        ];
        for side in [CanonicalCutSide::Left, CanonicalCutSide::Right] {
            vertices.push(VertexKey::ThresholdSide(side));
            for (orders, loops) in [(&none, 0), (&none, 1), (&qcd, 0)] {
                vertices.push(VertexKey::CutSide(side, orders, loops));
            }
        }
        let mut edges = vec![EdgeKey::CutPair, EdgeKey::CutMembership];
        for particle in [0, 1].map(ParticleId::from_index) {
            for external in [
                None,
                Some((0, ExternalState::Outgoing, 1)),
                Some((1, ExternalState::Incoming, 0)),
                Some((1, ExternalState::Outgoing, 0)),
            ] {
                for is_dummy in [false, true] {
                    edges.push(EdgeKey::Diagram {
                        particle,
                        external,
                        is_dummy,
                    });
                }
            }
        }
        for (&a, &b) in vertices
            .iter()
            .flat_map(|a| vertices.iter().map(move |b| (a, b)))
        {
            let public = CanonicalDiagramVertex::from(a).cmp(&CanonicalDiagramVertex::from(b));
            assert_eq!(a.cmp(&b), public, "{a:?} vs {b:?}");
        }
        for (&a, &b) in edges.iter().flat_map(|a| edges.iter().map(move |b| (a, b))) {
            let public = CanonicalDiagramEdge::from(a).cmp(&CanonicalDiagramEdge::from(b));
            assert_eq!(a.cmp(&b), public, "{a:?} vs {b:?}");
        }

        // The legs of the bubble, crossed so that their index and connection orders differ, are
        // keyed by their own index, state and connection.
        let crossed = one_loop()
            .map_data(
                |_, vertex| vertex.clone(),
                |_, _, edge| {
                    let mut edge = edge.clone();
                    if let Some(external) = &mut edge.external {
                        external.connection = 1 - external.index;
                    }
                    edge
                },
            )
            .unwrap();
        let mut legs: Vec<_> = crossed
            .edges()
            .filter_map(|(_, _, edge)| edge.external.clone())
            .map(|leg| ExternalLeg {
                name: String::new(),
                ..leg
            })
            .collect();
        legs.sort();
        let keyed: Vec<_> = crossed
            .canonical_key()
            .unwrap()
            .edges
            .into_iter()
            .filter_map(|edge| match edge {
                CanonicalDiagramEdge::Diagram { external, .. } => external,
                _ => None,
            })
            .collect();
        assert_eq!(keyed, legs);
    }

    /// Ids identify the colored topology: not names, insertion order or model parameters, but
    /// interactions, external legs and the names of particles and rules.
    #[test]
    fn diagram_ids_identify_the_structural_colored_topology() {
        let diagram = one_loop();
        assert_eq!(diagram.id().unwrap(), one_loop_reordered().id().unwrap());

        let mut heavier = scalar_model().as_ref().clone();
        heavier
            .apply_parameter_card(&ParameterCard::from_iter([(
                "M".to_owned(),
                ComplexValue::new(2.0, 0.0),
            )]))
            .unwrap();
        let heavier = one_loop_in(Arc::new(heavier));
        assert_ne!(heavier.model().fingerprint(), diagram.model().fingerprint());
        assert_eq!(heavier.id().unwrap(), diagram.id().unwrap());

        let renamed = |from: &str, to: &str| {
            let json = scalar_model().to_json().unwrap().replace(from, to);
            one_loop_in(Arc::new(Model::from_json(&json).unwrap()))
                .id()
                .unwrap()
        };
        assert_ne!(renamed("\"phi\"", "\"chi\""), diagram.id().unwrap());
        assert_ne!(renamed("\"V_1\"", "\"V_2\""), diagram.id().unwrap());

        let changed_connection = diagram
            .map_data(
                |_, vertex| vertex.clone(),
                |_, _, edge| {
                    let mut edge = edge.clone();
                    if let Some(external) = &mut edge.external
                        && external.index == 0
                    {
                        external.connection = 2;
                    }
                    edge
                },
            )
            .unwrap();
        assert_ne!(changed_connection.id().unwrap(), diagram.id().unwrap());
        let changed_interaction = diagram
            .map_data(
                |id, vertex| {
                    let mut vertex = vertex.clone();
                    if id == VertexId(0) {
                        vertex.interaction = None;
                    }
                    vertex
                },
                |_, _, edge| edge.clone(),
            )
            .unwrap();
        assert_ne!(changed_interaction.id().unwrap(), diagram.id().unwrap());
    }

    /// The normalized form ([`FeynmanDiagram::canonical_key`]) keys every antiparticle edge by
    /// its particle, directed ones against their direction, so it does not depend on the edge
    /// representation; the structural form and the id keep the stored particle and orientation.
    #[test]
    fn normalization_turns_antiparticle_edges_into_particle_edges() {
        let line = directed_fermion_line();
        let reversed = line.reverse_edge(EdgeId(0)).unwrap();
        assert_eq!(
            reversed.canonical_key().unwrap(),
            line.canonical_key().unwrap()
        );
        assert_ne!(reversed.id().unwrap(), line.id().unwrap());

        // An undirected edge between vertices without interaction, stored as `f` or as `f~`.
        let model = fermion_model();
        let undirected = |particle: &str| {
            let mut builder = FeynmanDiagram::builder(Arc::clone(&model), "undirected");
            let [a, b] = ["a", "b"].map(|name| {
                builder.add_vertex(DiagramVertex {
                    name: name.into(),
                    interaction: None,
                    numerator: Atom::one(),
                })
            });
            let particle = model.particle_id(particle).unwrap();
            builder
                .add_edge(a, b, DiagramEdge::new(particle, false))
                .unwrap();
            builder.build().unwrap()
        };
        let (particle, antiparticle) = (undirected("f"), undirected("f~"));
        let key = particle.canonical_key().unwrap();
        assert_eq!(antiparticle.canonical_key().unwrap(), key);
        assert!(matches!(
            key.edges[..],
            [CanonicalDiagramEdge::Diagram { particle, .. }]
                if particle == model.particle_id("f").unwrap()
        ));
        assert_ne!(antiparticle.id().unwrap(), particle.id().unwrap());
    }
}
