//! Independent native graphs for exact massive numerator/pinch identities.
//!
//! Neither graph is obtained by reducing or canonicalizing the other. A forest
//! contraction specifies the right incidence, and native momentum maps prove
//! the identity before either side reaches Vakint's matcher or the IBP engine.

use std::collections::BTreeMap;

use feynkit_graph::{DiagramEdge, DiagramVertex, EdgeId, FeynmanDiagram, IntegralFamily};
use feynkit_kinematics::Kinematics;
use feynkit_model::Model;
use idenso::tensor::SymbolicTensor;
use symbolica::{
    atom::{Atom, AtomCore},
    id::Replacement,
};
use vakint::{Vakint, VakintExpression, diagram_integral::DiagramIntegralOptions};
use vakint::{vakint_parse, vakint_symbol};

#[allow(dead_code)]
#[path = "../../src/topologies/five_loop_data.rs"]
mod data;
use data as five_loop_data;

#[path = "../../examples/five_loop_cancellation_fixtures/mod.rs"]
pub mod factorized;

pub type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

#[derive(Clone)]
pub struct Case {
    pub name: &'static str,
    pub parent: Option<usize>,
    /// Original edge indices to contract before constructing the numerator side.
    pub prepinches: &'static [usize],
    pub dots: &'static [(usize, i64)],
    pub cancellations: &'static [(usize, i64)],
    pub mass: i64,
    pub reroute: bool,
}

pub fn cases() -> Vec<Case> {
    vec![
        Case {
            name: "parent_32745_single",
            parent: Some(0),
            prepinches: &[0, 1],
            dots: &[],
            cancellations: &[(2, 1)],
            mass: 1,
            reroute: false,
        },
        Case {
            name: "parent_31740_single",
            parent: Some(1),
            prepinches: &[0, 1],
            dots: &[],
            cancellations: &[(4, 1)],
            mass: 1,
            reroute: false,
        },
        Case {
            name: "parent_30699_multiple",
            parent: Some(2),
            prepinches: &[0, 1],
            dots: &[],
            cancellations: &[(2, 1), (4, 1)],
            mass: 1,
            reroute: false,
        },
        Case {
            name: "parent_30527_multiple",
            parent: Some(3),
            prepinches: &[0, 1],
            dots: &[],
            cancellations: &[(2, 1), (4, 1)],
            mass: 1,
            reroute: false,
        },
        Case {
            name: "parent_32745_partial_dot",
            parent: Some(0),
            prepinches: &[0, 1, 2],
            dots: &[(3, 2)],
            cancellations: &[(3, 1)],
            mass: 1,
            reroute: false,
        },
        Case {
            name: "parent_32745_complete_dot",
            parent: Some(0),
            prepinches: &[0, 1],
            dots: &[(2, 2)],
            cancellations: &[(2, 2)],
            mass: 1,
            reroute: false,
        },
        Case {
            name: "parent_31740_mixed_dots",
            parent: Some(1),
            prepinches: &[0, 1, 2],
            dots: &[(3, 2), (4, 2)],
            cancellations: &[(3, 1), (4, 2)],
            mass: 1,
            reroute: false,
        },
        Case {
            name: "parent_30699_nonunit_routing",
            parent: Some(2),
            prepinches: &[0, 1],
            dots: &[],
            cancellations: &[(2, 1), (4, 1)],
            mass: 2,
            reroute: true,
        },
        Case {
            name: "ladder_left_partial_attachment",
            parent: None,
            prepinches: &[],
            dots: &[],
            cancellations: &[(0, 1), (2, 1), (3, 1), (4, 1)],
            mass: 1,
            reroute: false,
        },
        Case {
            name: "ladder_right_complete_attachment",
            parent: None,
            prepinches: &[],
            dots: &[],
            cancellations: &[(1, 2), (2, 1), (3, 1)],
            mass: 1,
            reroute: false,
        },
        Case {
            name: "ladder_attachment_and_rung",
            parent: None,
            prepinches: &[],
            dots: &[],
            cancellations: &[(0, 2), (9, 1), (4, 1)],
            mass: 1,
            reroute: false,
        },
        Case {
            name: "ladder_repeated_quadratics_nonunit",
            parent: None,
            prepinches: &[],
            dots: &[],
            cancellations: &[(0, 1), (1, 1), (2, 1), (5, 1)],
            mass: 2,
            reroute: true,
        },
    ]
}

struct Incidence {
    edges: Vec<(usize, usize)>,
    powers: Vec<i64>,
}

impl Incidence {
    fn for_case(case: &Case) -> Self {
        let (edges, mut powers) = if let Some(index) = case.parent {
            let parent = &data::PARENTS[index];
            let edges = parent
                .physical_slots
                .iter()
                .map(|slot| {
                    let id = (*slot + 1) as i64;
                    (
                        parent
                            .vertices
                            .iter()
                            .position(|row| row.contains(&id))
                            .unwrap(),
                        parent
                            .vertices
                            .iter()
                            .position(|row| row.contains(&-id))
                            .unwrap(),
                    )
                })
                .collect();
            (edges, vec![1; 12])
        } else {
            // Four-rung vacuum ladder, independently written in rail coordinates.
            // Vertices 0..4 are t1..t4, and 4..8 are b1..b4. The first two
            // edges are the collapsed external-current attachments, each dotted.
            (
                vec![
                    (0, 4),
                    (3, 7),
                    (0, 1),
                    (1, 2),
                    (2, 3),
                    (5, 4),
                    (6, 5),
                    (7, 6),
                    (0, 4),
                    (1, 5),
                    (2, 6),
                    (3, 7),
                ],
                vec![2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            )
        };
        for &edge in case.prepinches {
            powers[edge] = 0;
        }
        for &(edge, power) in case.dots {
            powers[edge] = power;
        }
        Self { edges, powers }
    }

    fn diagram(&self, model: &Model, name: &str, powers: &[i64]) -> Result<NativeGraph> {
        let mut roots = (0..8).collect::<Vec<_>>();
        for ((left, right), &power) in self.edges.iter().zip(powers) {
            if power != 0 {
                continue;
            }
            let source = roots[*left];
            let target = roots[*right];
            if source == target {
                return Err(format!("{name}: contractions remove an independent loop").into());
            }
            for root in &mut roots {
                if *root == target {
                    *root = source;
                }
            }
        }
        let particle = model.particle_id("phi")?;
        let mut builder = FeynmanDiagram::builder(model.clone(), name);
        let vertices = roots
            .iter()
            .copied()
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .map(|root| {
                (
                    root,
                    builder.add_vertex(DiagramVertex {
                        name: format!("v{root}"),
                        interaction: None,
                        numerator: SymbolicTensor::infer(Atom::one())
                            .expect("literal one is scalar"),
                    }),
                )
            })
            .collect::<BTreeMap<_, _>>();
        let mut slots = Vec::new();
        for (slot, (&power, &(left, right))) in powers.iter().zip(&self.edges).enumerate() {
            if power == 0 {
                continue;
            }
            if power < 0 {
                return Err("closure graph has a negative physical power".into());
            }
            let mut edge = DiagramEdge::new(particle, false);
            edge.propagator_power = i32::try_from(power)?;
            let id = builder.add_edge(vertices[&roots[left]], vertices[&roots[right]], edge)?;
            if id.0 != slots.len() {
                return Err("native graph changed ordered edge identity".into());
            }
            slots.push(slot);
        }
        // Stable DOT is an actual tested boundary, not an interchangeable string.
        let dot = builder.build()?.to_dot()?;
        let diagram = FeynmanDiagram::from_dot(model.clone(), &dot)?;
        if diagram.edges().count() != slots.len() {
            return Err("stable DOT changed the physical edge count".into());
        }
        for (id, _, edge) in diagram.edges() {
            let slot = slots
                .get(id.0)
                .ok_or("stable DOT changed a physical edge identity")?;
            if i64::from(edge.propagator_power) != powers[*slot] {
                return Err("stable DOT changed a physical propagator power".into());
            }
        }
        diagram.loop_momentum_basis().validate(&diagram)?;
        if diagram.loop_count() != 5 {
            return Err("closure graph changed loop count".into());
        }
        Ok(NativeGraph { diagram, slots })
    }
}

struct NativeGraph {
    diagram: FeynmanDiagram,
    /// Original incidence edge for each stable native edge ID.
    slots: Vec<usize>,
}

impl NativeGraph {
    /// Cycle ranks of the vertex-biconnected vacuum factors. In particular,
    /// connected incidence with tadpole factors is not called an irreducible 5L
    /// control. Repeated undirected edges keep distinct incidence identities.
    fn loop_blocks(&self) -> Vec<usize> {
        fn split(edges: &[(usize, usize)], blocks: &mut Vec<usize>) {
            blocks.extend(
                edges
                    .iter()
                    .filter(|(left, right)| left == right)
                    .map(|_| 1),
            );
            let edges = edges
                .iter()
                .copied()
                .filter(|(left, right)| left != right)
                .collect::<Vec<_>>();
            let nodes = edges
                .iter()
                .flat_map(|&(left, right)| [left, right])
                .collect::<std::collections::BTreeSet<_>>();
            for cut in &nodes {
                let mut remaining = nodes.clone();
                remaining.remove(cut);
                let mut components = Vec::new();
                while let Some(start) = remaining.pop_first() {
                    let mut found = std::collections::BTreeSet::from([start]);
                    let mut queue = vec![start];
                    while let Some(node) = queue.pop() {
                        for &(left, right) in &edges {
                            let next = if left == node {
                                right
                            } else if right == node {
                                left
                            } else {
                                continue;
                            };
                            if next != *cut && remaining.remove(&next) {
                                found.insert(next);
                                queue.push(next);
                            }
                        }
                    }
                    components.push(found);
                }
                if components.len() > 1 {
                    for component in components {
                        let subgraph = edges
                            .iter()
                            .copied()
                            .filter(|(left, right)| {
                                (component.contains(left) || left == cut)
                                    && (component.contains(right) || right == cut)
                            })
                            .collect::<Vec<_>>();
                        split(&subgraph, blocks);
                    }
                    return;
                }
            }
            if !edges.is_empty() && edges.len() >= nodes.len() {
                blocks.push(edges.len() + 1 - nodes.len());
            }
        }
        let edges = self
            .diagram
            .edges()
            .map(|(_, ends, _)| {
                (
                    ends.source.expect("vacuum edge has an internal source").0,
                    ends.target.expect("vacuum edge has an internal target").0,
                )
            })
            .collect::<Vec<_>>();
        let mut blocks = Vec::new();
        split(&edges, &mut blocks);
        blocks.sort();
        blocks
    }

    fn family(&self, mass_symbol: &Atom, mass: i64) -> Result<IntegralFamily> {
        let family = self.diagram.propagator_family(&Kinematics::new())?;
        let replacement = [Replacement::new(mass_symbol.clone(), Atom::num(mass))];
        Ok(IntegralFamily::new(
            family.loop_momenta().to_vec(),
            vec![],
            family
                .denominators()
                .iter()
                .map(|d| d.replace_multiple(&replacement))
                .collect(),
            family.kinematics(),
        )?)
    }
}

pub struct Pair {
    pub case: Case,
    pub factored: Atom,
    pub expanded: Atom,
    pub pinched: Atom,
    pub source_dot: String,
    pub pinched_dot: String,
    pub source_numerator: Atom,
    pub source_powers: Vec<i64>,
    pub pinched_powers: Vec<i64>,
    pub source_distinct_quadratics: usize,
    pub pinched_loop_blocks: Vec<usize>,
}

impl Pair {
    pub fn new(case: Case) -> Result<Self> {
        Vakint::initialize_vakint_symbols();
        let model = Model::from_json(include_str!("../fixtures/phi_3_4.json"))?;
        let mass_symbol = model
            .particle_by_id(model.particle_id("phi")?)?
            .symbolic_mass(&model);
        let incidence = Incidence::for_case(&case);
        let mut reduced_powers = incidence.powers.clone();
        for &(slot, power) in case.cancellations {
            if power <= 0 || reduced_powers[slot] < power {
                return Err(format!("{}: invalid requested cancellation", case.name).into());
            }
            reduced_powers[slot] -= power;
        }
        let mut source = incidence.diagram(&model, case.name, &incidence.powers)?;
        let target =
            incidence.diagram(&model, &format!("{}_pinched", case.name), &reduced_powers)?;
        // Target chords survive the forest contraction and form a valid source
        // basis. The alternative order exercises independently stored routing.
        let mut source_chords = target
            .diagram
            .loop_momentum_basis()
            .loop_edges
            .iter()
            .map(|edge| {
                let slot = target.slots[edge.0];
                source
                    .slots
                    .iter()
                    .position(|candidate| *candidate == slot)
                    .map(EdgeId)
                    .ok_or_else(|| "target chord absent from numerator graph".into())
            })
            .collect::<Result<Vec<_>>>()?;
        if case.reroute {
            source_chords.rotate_left(2);
            source_chords.reverse();
        }
        source.diagram = source.diagram.with_loop_momentum_edges(&source_chords)?;
        let source_family = source.family(&mass_symbol, case.mass)?;
        let target_family = target.family(&mass_symbol, case.mass)?;

        let mut images = vec![Atom::Zero; 5];
        for chord in &source_chords {
            let source_signature = source.diagram.loop_momentum_basis().edge_signatures[chord]
                .loops
                .integer_coefficients();
            let entries = source_signature
                .iter()
                .enumerate()
                .filter(|(_, coefficient)| **coefficient != 0)
                .collect::<Vec<_>>();
            let [(axis, coefficient)] = entries.as_slice() else {
                return Err("source chord is not an independent loop coordinate".into());
            };
            if **coefficient != 1 && **coefficient != -1 {
                return Err("source chord lacks unit normalization".into());
            }
            let slot = source.slots[chord.0];
            let target_id = EdgeId(
                target
                    .slots
                    .iter()
                    .position(|candidate| *candidate == slot)
                    .ok_or("source chord was pinched")?,
            );
            images[*axis] = target.diagram.loop_momentum_basis().edge_signatures[&target_id]
                .loops
                .integer_coefficients()
                .iter()
                .zip(target_family.loop_momenta())
                .map(|(value, momentum)| Atom::num(*value as i64 * **coefficient as i64) * momentum)
                .sum();
        }
        let retained = source
            .slots
            .iter()
            .enumerate()
            .filter(|(_, slot)| reduced_powers[**slot] > 0)
            .map(|(edge, _)| source_family.denominators()[edge].clone())
            .collect();
        let retained_family = IntegralFamily::new(
            source_family.loop_momenta().to_vec(),
            vec![],
            retained,
            source_family.kinematics(),
        )?;
        let mapping = retained_family
            .mapping_to(&target_family, &images)?
            .ok_or("independent pinch lacks a simultaneous unit-Jacobian momentum map")?;
        let mut source_weights = BTreeMap::<Atom, i64>::new();
        for (edge, &slot) in source.slots.iter().enumerate() {
            if reduced_powers[slot] > 0 {
                *source_weights
                    .entry(mapping.apply(&source_family.denominators()[edge]).expand())
                    .or_default() += reduced_powers[slot];
            }
        }
        let mut target_weights = BTreeMap::<Atom, i64>::new();
        for (edge, &slot) in target.slots.iter().enumerate() {
            *target_weights
                .entry(target_family.denominators()[edge].expand())
                .or_default() += reduced_powers[slot];
        }
        if source_weights != target_weights {
            return Err("independent pinch changed denominator powers or mass".into());
        }

        let numerator = case
            .cancellations
            .iter()
            .map(|&(slot, power)| {
                let position = source
                    .slots
                    .iter()
                    .position(|candidate| *candidate == slot)
                    .ok_or("canceled line absent from source")?;
                Ok(source_family.denominators()[position]
                    .clone()
                    .pow(Atom::num(power)))
            })
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .product::<Atom>();
        let options = DiagramIntegralOptions {
            parameter_substitutions: vec![(mass_symbol, Atom::num(case.mass))],
            ..Default::default()
        };
        let factored = Atom::from(VakintExpression::from_diagram(
            &source.diagram,
            &source_family,
            &numerator,
            &options,
        )?);
        let expanded = Atom::from(VakintExpression::from_diagram(
            &source.diagram,
            &source_family,
            &numerator.expand(),
            &options,
        )?);
        let pinched = Atom::from(VakintExpression::from_diagram(
            &target.diagram,
            &target_family,
            &Atom::one(),
            &options,
        )?);
        let source_distinct_quadratics = source_family
            .denominators()
            .iter()
            .cloned()
            .collect::<std::collections::BTreeSet<_>>()
            .len();
        let pinched_loop_blocks = target.loop_blocks();
        if pinched_loop_blocks.iter().sum::<usize>() != 5 {
            return Err("pinched factor blocks changed total loop count".into());
        }
        Ok(Self {
            case,
            factored,
            expanded,
            pinched,
            source_dot: source
                .diagram
                .clone()
                .with_numerator(SymbolicTensor::infer(numerator.clone())?)?
                .to_dot()?,
            pinched_dot: target.diagram.to_dot()?,
            source_numerator: numerator,
            source_powers: source
                .slots
                .iter()
                .map(|slot| incidence.powers[*slot])
                .collect(),
            pinched_powers: target
                .slots
                .iter()
                .map(|slot| reduced_powers[*slot])
                .collect(),
            source_distinct_quadratics,
            pinched_loop_blocks,
        })
    }
}
