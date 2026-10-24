//! Complete model-generated covariant cuts, tested against three physical
//! polarizations. Component contraction never expands a physical numerator.

use feynkit_generator::{
    GenerationFilter, GenerationOptions, GenerationType, Process as GenerationProcess,
};
use std::collections::BTreeMap;

use itertools::Itertools;
use linnet::half_edge::{
    involution::{EdgeIndex, Flow, HedgePair},
    subgraph::Inclusion,
};
use spenso::{
    network::{
        ExecutionResult, Sequential, SmallestDegree,
        library::{
            LibraryTensor,
            symbolic::{ETS, ExplicitKey, TensorLibrary},
        },
        parsing::ParseSettings,
    },
    structure::{
        Canonicalized,
        representation::{Minkowski, RepName},
    },
    tensors::parametric::ParamTensor,
};
use symbolica::{
    atom::{Atom, AtomCore},
    function, parse,
};

use crate::{
    graph::{FinalizedCut, Graph},
    initialisation::test_initialise,
    model::{
        Model, ModelGammaLoopExt, ParticleGammaLoopExt, ParticleId, ParticleIdGammaLoopExt,
        UFOSymbol,
    },
    numerator::{ParsingNet, aind::Aind},
    processes::{Process, ProcessCollection, ProcessDefinition},
    settings::{GlobalSettings, global::VectorPolarizationSumGauge},
    utils::{FUN_LIB, GS, W_, load_generic_model},
};

type ExactMomentum = [Atom; 4];
/// An external leg: its PDG code and exact momentum.
type Leg = (i64, ExactMomentum);

fn minkowski_square(p: &ExactMomentum) -> Atom {
    (&p[0] * &p[0] - &p[1] * &p[1] - &p[2] * &p[2] - &p[3] * &p[3]).expand()
}

fn scaled(sign: i64, p: &ExactMomentum) -> ExactMomentum {
    std::array::from_fn(|mu| Atom::num(sign) * &p[mu])
}

fn substitute(expression: &Atom, parameters: &[(&str, Atom)]) -> Atom {
    parameters
        .iter()
        .fold(expression.clone(), |expression, (name, value)| {
            expression
                .replace(Atom::from(UFOSymbol::from(*name)).to_pattern())
                .with(value.clone())
        })
}

/// Contract a graph's complete production numerator at exact edge momenta and
/// parameter values. The generic metric has Atom entries so that its spatial
/// signs remain exact during component contraction.
fn production_numerator_at(
    graph: &Graph,
    model: &Model,
    momenta: &BTreeMap<EdgeIndex, ExactMomentum>,
    parameters: &[(&str, Atom)],
) -> eyre::Result<Atom> {
    let mut library = spenso_hep_lib::hep_lib_atom::<Aind, ParamTensor<ExplicitKey<Aind>>>();
    library.insert_generic(
        TensorLibrary::<ParamTensor<ExplicitKey<Aind>>, Aind>::id(Minkowski {}.into()),
        TensorLibrary::<ParamTensor<ExplicitKey<Aind>>, Aind>::diag_unimodular_metric,
    );
    for (edge, momentum) in momenta {
        let key = ExplicitKey::from_iter(
            [Minkowski {}.new_rep(4)],
            GS.emr_mom,
            Some(vec![Atom::num(edge.0)]),
        );
        let tensor = ParamTensor::from_dense(key.into_canonical(), momentum.to_vec())?;
        library.insert_explicit(Canonicalized::identity(tensor));
    }
    let expression = substitute(
        &model.expand_couplings(&graph.production_numerator_atom_for_full_3d_expression()),
        parameters,
    );
    let mut network =
        ParsingNet::try_from_view(expression.as_view(), &library, &ParseSettings::default())?;
    network.execute::<Sequential, SmallestDegree, _, _, _>(&library, &*FUN_LIB)?;
    // Read the scalar Atom directly instead of materializing a rank-zero
    // tensor and extracting its component again.
    Ok(match network.result_scalar()? {
        ExecutionResult::Zero => Atom::Zero,
        ExecutionResult::One => Atom::one(),
        ExecutionResult::Val(value) => value.into_owned(),
    })
}

/// The particle crossing a split cut edge from left to right, with the sign
/// relating its outgoing momentum to the edge momentum.
fn left_to_right(pair: HedgePair, particle: ParticleId, model: &Model) -> (ParticleId, i64) {
    match pair {
        HedgePair::Split { split: flow, .. } | HedgePair::Unpaired { flow, .. } => {
            if flow == Flow::Source {
                (particle, 1)
            } else {
                (particle.antiparticle(model), -1)
            }
        }
        HedgePair::Paired { .. } => panic!("cut edge must be split"),
    }
}

/// Initial-state momenta enter the left amplitude. An initial-state edge whose
/// source lies on the left carries the antiparticle and minus that momentum.
fn enters_against_edge(pair: HedgePair, cut: &FinalizedCut) -> bool {
    let source = match pair {
        HedgePair::Paired { source, .. } | HedgePair::Split { source, .. } => source,
        HedgePair::Unpaired { hedge, .. } => hedge,
    };
    cut.left.includes(&source)
}

/// Exact momenta of every edge of one finalized cut, with the width-free
/// uncut propagator denominators that the production numerator omits.
struct CutKinematics {
    momenta: BTreeMap<EdgeIndex, ExactMomentum>,
    /// PDG codes crossing the cut from left to right, sorted.
    channel: Vec<i64>,
    /// Cut vector edges with their symbolic masses.
    vectors: Vec<(EdgeIndex, Atom)>,
    denominator: Atom,
}

impl CutKinematics {
    /// Cut edges take the outgoing momentum of the first unused final-state
    /// leg with their physical label, initial-state edges their incoming
    /// momentum, and conservation fixes the tree edges on both sides.
    fn new(
        graph: &Graph,
        model: &Model,
        cut: &FinalizedCut,
        labels: &BTreeMap<isize, isize>,
        (initial, finals): (&[Leg], &[Leg]),
        parameters: &[(&str, Atom)],
    ) -> Self {
        let mut momenta = BTreeMap::new();
        let mut channel = Vec::new();
        let mut vectors = Vec::new();
        let mut unused = finals.iter().collect_vec();
        for (pair, edge, data) in graph.iter_edges_of(&cut.cut) {
            let (particle, sign) = left_to_right(pair, data.data.particle().unwrap(), model);
            let particle = particle.resolve(model);
            let label = labels
                .get(&(particle.pdg_code as isize))
                .map_or(particle.pdg_code, |label| *label as i64);
            let leg = unused
                .iter()
                .position(|(pdg, _)| *pdg == label)
                .expect("every cut particle labels an unused final-state leg");
            momenta.insert(edge, scaled(sign, &unused.remove(leg).1));
            channel.push(particle.pdg_code);
            if particle.spin == 3 {
                vectors.push((edge, particle.symbolic_mass(model)));
            }
        }
        let mut unused = initial.iter().collect_vec();
        for (pair, edge, data) in graph.iter_edges_of(&graph.initial_state_cut) {
            let against = enters_against_edge(pair, cut);
            let particle = data.data.particle().unwrap();
            let pdg = if against {
                particle.antiparticle(model)
            } else {
                particle
            }
            .resolve(model)
            .pdg_code;
            let leg = unused
                .iter()
                .position(|(incoming, _)| *incoming == pdg)
                .expect("every initial-state edge carries an unused incoming leg");
            let sign = if against { -1 } else { 1 };
            momenta.insert(edge, scaled(sign, &unused.remove(leg).1));
        }

        // Momentum entering minus leaving a node, and its unrouted edges.
        let balance = |momenta: &BTreeMap<EdgeIndex, ExactMomentum>, node| {
            let mut balance: ExactMomentum = std::array::from_fn(|_| Atom::Zero);
            let mut unrouted = Vec::new();
            for hedge in graph.underlying.iter_crown(node) {
                let edge = graph.underlying[&hedge];
                let sign = if graph.underlying.flow(hedge) == Flow::Sink {
                    1
                } else {
                    -1
                };
                match momenta.get(&edge) {
                    Some(p) => balance = std::array::from_fn(|mu| &balance[mu] + sign * &p[mu]),
                    None => unrouted.push((edge, sign)),
                }
            }
            (balance, unrouted)
        };
        let mut denominator = Atom::one();
        while let Some((edge, momentum)) = graph.underlying.iter_nodes().find_map(|(node, _, _)| {
            let (balance, unrouted) = balance(&momenta, node);
            let [(edge, sign)] = unrouted[..] else {
                return None;
            };
            Some((edge, scaled(-sign, &balance).map(|p| p.expand())))
        }) {
            let mass = graph.underlying[edge]
                .particle()
                .unwrap()
                .resolve(model)
                .symbolic_mass(model);
            denominator *= minkowski_square(&momentum) - substitute(&mass, parameters).pow(2);
            momenta.insert(edge, momentum);
        }
        assert_eq!(momenta.len(), graph.underlying.n_edges(), "{}", graph.name);
        for (node, _, _) in graph.underlying.iter_nodes() {
            let (balance, _) = balance(&momenta, node);
            assert!(
                balance.iter().all(|p| p.expand().is_zero()),
                "{}",
                graph.name
            );
        }
        let denominator = denominator.expand();
        assert!(
            !denominator.is_zero(),
            "on-shell propagator in {}",
            graph.name
        );
        channel.sort();
        Self {
            momenta,
            channel,
            vectors,
            denominator,
        }
    }

    /// The covariant cut, and for a vector-only channel the same cut with
    /// every vector summed over its physical polarizations, -g + k k / M^2.
    fn evaluate(
        &self,
        graph: &Graph,
        model: &Model,
        parameters: &[(&str, Atom)],
    ) -> eyre::Result<(Atom, Option<Atom>)> {
        let value = |graph: &Graph| -> eyre::Result<Atom> {
            let numerator = production_numerator_at(graph, model, &self.momenta, parameters)?;
            Ok((numerator / &self.denominator).expand())
        };
        let covariant = value(graph)?;
        if self.vectors.len() != self.channel.len() {
            return Ok((covariant, None));
        }
        let mut physical = graph.clone();
        for (edge, mass) in &self.vectors {
            let metric = function!(ETS.metric, W_.a_, W_.b_);
            let proca = &metric - GS.emr_mom(*edge, W_.a_) * GS.emr_mom(*edge, W_.b_) / mass.pow(2);
            let numerator = &physical.underlying[*edge].num.value;
            let replaced = numerator
                .expression()
                .replace(metric.to_pattern())
                .with(proca.to_pattern());
            assert_ne!(
                &replaced,
                numerator.expression(),
                "cut vector without metric"
            );
            physical.underlying[*edge].num.value = numerator
                .with_checked_expression(replaced)
                .map_err(|error| eyre::eyre!("{error}"))?;
        }
        Ok((covariant, Some(value(&physical)?)))
    }
}

#[test]
fn declared_covariant_cut_states_preserve_multisets_and_physical_labels() -> eyre::Result<()> {
    test_initialise()?;
    let model = load_generic_model("sm");
    for (state, count) in [(vec![24, -24], 16), (vec![23, 23], 10)] {
        let process = ProcessDefinition {
            generation_type: feynkit_generator::GenerationType::CrossSection,
            process: GenerationProcess::new(Vec::<i64>::new(), state.clone()),
            ..Default::default()
        };
        let closure = process.covariant_cut_states(&model)?;
        assert_eq!(closure.len(), count);
        assert_eq!(model.covariant_cut_states(&closure)?, closure);
        let labels = process.covariant_cut_representatives(&model);
        for members in closure {
            assert_eq!(
                members
                    .iter()
                    .map(|pdg| labels[&(*pdg as isize)] as i64)
                    .sorted()
                    .collect_vec(),
                state.iter().copied().sorted().collect_vec()
            );
        }
    }
    let diagnostic = ProcessDefinition {
        generation_type: feynkit_generator::GenerationType::CrossSection,
        process: GenerationProcess::new(Vec::<i64>::new(), [250_i64, 250]),
        ..Default::default()
    };
    assert!(diagnostic.covariant_cut_representatives(&model).is_empty());
    assert_eq!(
        diagnostic.covariant_cut_states(&model)?,
        diagnostic.process.outgoing_pdgs(&model)?
    );
    let mixed = ProcessDefinition {
        process: diagnostic
            .process
            .clone()
            .with_final_state_alternatives([vec![23_i64, 23], vec![250, 25]])?,
        ..diagnostic.clone()
    };
    assert!(
        mixed
            .covariant_cut_states(&model)
            .unwrap_err()
            .to_string()
            .contains("observable labels remain unambiguous")
    );
    let unresolved = ProcessDefinition {
        generation_options: GenerationOptions::default().with_graph_filter(
            GenerationFilter::PerturbativeOrders([("QED".into(), 1)].into_iter().collect()),
        ),
        ..diagnostic
    };
    let mut unresolved_model = model.clone();
    // Native unresolved sets derive from massless interactions. A massless
    // restriction of a declared massive-vector quartet is invalid as well.
    let mut card = feynkit_model::ParameterCard::new();
    card.insert("MZ".into(), feynkit_model::ComplexValue::new(0.0, 0.0));
    unresolved_model.apply_parameter_card(&card)?;
    assert!(unresolved.covariant_cut_states(&unresolved_model).is_err());

    let incomplete = ProcessDefinition {
        generation_type: feynkit_generator::GenerationType::CrossSection,
        process: GenerationProcess::new(Vec::<i64>::new(), [24_i64, -24]),
        generation_options: GenerationOptions::default()
            .with_graph_filter(GenerationFilter::ParticleVeto(vec![251_i64.into()])),
        ..Default::default()
    };
    assert!(
        incomplete
            .covariant_cut_states(&model)
            .unwrap_err()
            .to_string()
            .contains("complete vector/Goldstone/ghost sector")
    );
    let mut invalid = model.clone();
    invalid.covariant_cut_multiplets_mut().get_mut(&24).unwrap()[3] = -9000003;
    assert!(invalid.validate_covariant_cut_multiplets().is_err());
    let mut invalid = serde_json::to_value(&model)?;
    let propagator = invalid["propagators"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .find(|p| p["particle"] == "G0")
        .unwrap();
    let denominator = &model
        .propagators()
        .iter()
        .find(|p| p.particle.resolve(&model).name == "G0")
        .unwrap()
        .denominator;
    propagator["denominator"] =
        serde_json::to_value((denominator + Atom::one()).to_canonical_string())?;
    assert!(serde_json::from_value::<crate::model::Model>(invalid).is_err());
    Ok(())
}

#[test]
fn generated_higgs_covariant_cuts_equal_three_physical_vector_polarizations() -> eyre::Result<()> {
    test_initialise()?;
    let model = load_generic_model("sm");
    let mut settings = GlobalSettings::default();
    settings.n_cores.feyngen = 1;

    let dot = |a: &[Atom; 4], b: &[Atom; 4]| {
        &a[0] * &b[0] - &a[1] * &b[1] - &a[2] * &b[2] - &a[3] * &b[3]
    };
    for (physical, vector_mass, graph_count, symmetry) in [
        (vec![24, -24], Atom::num(3), 6, Atom::one()),
        (vec![23, 23], Atom::num((15, 4)), 4, Atom::num((1, 2))),
    ] {
        let process = ProcessDefinition {
            generation_type: feynkit_generator::GenerationType::CrossSection,
            process: GenerationProcess::new([25_i64], physical),
            generation_options: GenerationOptions::default()
                .with_loop_count(1, 1)?
                .with_graph_filter(GenerationFilter::CouplingOrders(
                    [("QED".into(), (2, Some(2)))].into_iter().collect(),
                )),
            ..Default::default()
        };
        let graphs = process.generate(&model, &settings)?;
        assert_eq!(
            graphs.len(),
            graph_count,
            "the complete covariant Higgs cut inventory"
        );
        assert!(!process.may_filter_covariant_partners(&model));
        if process.process.outgoing_pdgs(&model)?[0] == [24, -24] {
            let mut filtered = process.clone();
            // In H -> WW this retains the four bosonic Born graphs while
            // removing the two ghost graphs required by covariant completeness.
            filtered.generation_options =
                filtered
                    .generation_options
                    .with_graph_filter(GenerationFilter::VertexVeto(vec![
                        "V_17".into(),
                        "V_25".into(),
                    ]));
            assert!(filtered.may_filter_covariant_partners(&model));
            let bosonic_graphs = filtered.generate(&model, &settings)?;
            assert_eq!(bosonic_graphs.len(), 4);
            assert!(bosonic_graphs.iter().all(|graph| {
                graph.iter_edges().all(|(_, _, edge)| {
                    !edge
                        .data
                        .particle()
                        .is_some_and(|particle| particle.resolve(&model).is_ghost())
                })
            }));
            filtered.process = filtered
                .process
                .with_final_state_alternatives([vec![251_i64, -251]])?;
            assert!(!filtered.may_filter_covariant_partners(&model));
        }
        // Raw imported states cannot establish physical-vector intent. The
        // same collection with its physical process declaration retains every
        // partner through the import owner and the cut matcher below.
        let inferred =
            ProcessDefinition::from_graph_list(&graphs, GenerationType::CrossSection, &model)?;
        assert!(
            inferred
                .covariant_cut_states(&model)
                .unwrap_err()
                .to_string()
                .contains("--process-spec")
        );
        let imported = Process::from_graph_list(
            "higgs_gauge".into(),
            "born".into(),
            graphs,
            GenerationType::CrossSection,
            Some(process),
            None,
            &model,
        )?;
        let ProcessCollection::CrossSections(cross_sections) = &imported.collection else {
            panic!("forward graphs must import as a cross section");
        };
        let process = &imported.definition;
        for (energy_ratio, momentum_ratio) in [((5, 3), (4, 3)), ((13, 5), (12, 5))] {
            let energy = &vector_mass * Atom::num(energy_ratio);
            let momentum = &vector_mass * Atom::num(momentum_ratio);
            let higgs_mass = Atom::num(2) * &energy;
            // ee=sw=3/5, cw=4/5, MW=3 and vev=6 obey the same action
            // relations as the model. Set MH^2=2 lambda vev^2 on shell.
            let coupling = vector_mass.pow(2) / Atom::num(3);
            for boosted in [false, true] {
                let boost = |p: [Atom; 4]| -> [Atom; 4] {
                    if !boosted {
                        return p;
                    }
                    [
                        Atom::num((5, 4)) * &p[0] + Atom::num((3, 4)) * &p[1],
                        Atom::num((3, 4)) * &p[0] + Atom::num((5, 4)) * &p[1],
                        p[2].clone(),
                        p[3].clone(),
                    ]
                };
                let q = [
                    boost([energy.clone(), Atom::Zero, Atom::Zero, momentum.clone()]),
                    boost([energy.clone(), Atom::Zero, Atom::Zero, -&momentum]),
                ];
                let total: [Atom; 4] = std::array::from_fn(|mu| &q[0][mu] + &q[1][mu]);
                let polarizations = [1, -1].map(|direction| {
                    [
                        boost([Atom::Zero, Atom::one(), Atom::Zero, Atom::Zero]),
                        boost([Atom::Zero, Atom::Zero, Atom::one(), Atom::Zero]),
                        boost([
                            &momentum / &vector_mass,
                            Atom::Zero,
                            Atom::Zero,
                            Atom::num(direction) * &energy / &vector_mass,
                        ]),
                    ]
                });
                let physical_norm = polarizations[0]
                    .iter()
                    .flat_map(|left| polarizations[1].iter().map(|right| dot(left, right).pow(2)))
                    .fold(Atom::Zero, |sum, norm| sum + norm)
                    * coupling.pow(2)
                    * &symmetry;
                assert_eq!(
                    physical_norm,
                    coupling.pow(2)
                        * &symmetry
                        * (Atom::num(2) + dot(&q[0], &q[1]).pow(2) / vector_mass.pow(4))
                );

                let mut channels = BTreeMap::<Vec<i64>, Atom>::new();
                for forward in &cross_sections["born"].supergraphs {
                    let graph = &forward.graph;
                    let cuts = forward.process_valid_cuts(&model, process, &settings.generation)?;
                    assert_eq!(cuts.len(), 1, "each Higgs Born graph has one cut");
                    let cut = &cuts[super::CutId(0)];
                    let mut particles = Vec::new();
                    let mut momenta = BTreeMap::new();
                    for (position, (pair, eid, edge)) in graph.iter_edges_of(&cut.cut).enumerate() {
                        let (particle, sign) =
                            left_to_right(pair, edge.data.particle().unwrap(), &model);
                        particles.push(particle.resolve(&model).pdg_code);
                        momenta.insert(eid, scaled(sign, &q[position]));
                    }
                    for (pair, eid, _) in graph.iter_edges_of(&graph.initial_state_cut) {
                        let sign = if enters_against_edge(pair, cut) {
                            -1
                        } else {
                            1
                        };
                        momenta.insert(eid, scaled(sign, &total));
                    }
                    let value = production_numerator_at(
                        graph,
                        &model,
                        &momenta,
                        &[
                            ("ee", Atom::num((3, 5))),
                            ("sw", Atom::num((3, 5))),
                            ("cw", Atom::num((4, 5))),
                            ("MW", Atom::num(3)),
                            ("MZ", Atom::num((15, 4))),
                            ("vev", Atom::num(6)),
                            ("lam", higgs_mass.pow(2) / Atom::num(72)),
                            ("MH", higgs_mass.clone()),
                        ],
                    )?;
                    *channels
                        .entry(particles.into_iter().sorted().collect())
                        .or_default() += value;
                }
                assert_eq!(channels.len(), graph_count);
                for (particles, value) in &channels {
                    let spins = particles
                        .iter()
                        .map(|pdg| {
                            model
                                .get_particle_from_pdg(*pdg as isize)
                                .resolve(&model)
                                .spin
                        })
                        .sorted()
                        .collect_vec();
                    let expected = match spins.as_slice() {
                        [3, 3] => Atom::num(4) * coupling.pow(2) * &symmetry,
                        [1, 1] => higgs_mass.pow(4) / Atom::num(36) * &symmetry,
                        [-1, -1] => -coupling.pow(2) / Atom::num(4),
                        [1, 3] => {
                            -coupling.pow(2)
                                * (Atom::num(2) * higgs_mass.pow(2) + vector_mass.pow(2))
                                / (Atom::num(4) * vector_mass.pow(2))
                        }
                        _ => panic!("unexpected Higgs gauge channel {particles:?}"),
                    };
                    assert_eq!(*value, expected, "channel {particles:?}, boosted={boosted}");
                }
                let covariant_norm = channels.values().fold(Atom::Zero, |sum, value| sum + value);
                assert_eq!(
                    covariant_norm, physical_norm,
                    "{channels:?}, boosted={boosted}"
                );
                // The VV channel alone is the incomplete covariant metric sum.
                let vector_channel = process.process.outgoing_pdgs(&model)?[0]
                    .iter()
                    .copied()
                    .sorted()
                    .collect_vec();
                assert_eq!(
                    channels[&vector_channel],
                    Atom::num(4) * coupling.pow(2) * &symmetry
                );
                assert_ne!(channels[&vector_channel], physical_norm);
            }
        }
    }
    Ok(())
}

/// One exact on-shell point: MH, m_tau, incoming and outgoing momenta.
type ExactPoint = (
    &'static str,
    &'static str,
    [[&'static str; 4]; 2],
    [[&'static str; 4]; 2],
);

/// tau- tau+ -> W+ W-. The s-channel photon, Z, G0 and H combine the neutral
/// ghost-vector vertices V_18/V_19/V_26/V_27 with V_16/V_24 and V_17/V_25.
const TAU_PAIR_TO_W_PAIR: [ExactPoint; 3] = [
    (
        "7",
        "7/5",
        [["5", "0", "0", "24/5"], ["5", "0", "0", "-24/5"]],
        [["5", "8/3", "4/3", "8/3"], ["5", "-8/3", "-4/3", "-8/3"]],
    ),
    (
        "5",
        "9/4",
        [
            ["483/80", "81/16", "0", "12/5"],
            ["267/80", "9/16", "0", "-12/5"],
        ],
        [
            ["579/112", "405/112", "-27/14", "27/28"],
            ["471/112", "225/112", "27/14", "-27/28"],
        ],
    ),
    (
        "11",
        "3",
        [
            ["1363/140", "72/35", "923/140", "216/35"],
            ["1003/140", "-72/35", "-13/140", "-216/35"],
        ],
        [
            ["209/20", "-12/5", "169/20", "24/5"],
            ["129/20", "12/5", "-39/20", "-24/5"],
        ],
    ),
];

/// tau- nu_tau~ -> W- Z. The s-channel W- and G- combine the charged
/// ghost-vector vertices V_21/V_29/V_31/V_33 with V_20/V_28/V_30/V_32.
const TAU_NEUTRINO_TO_W_Z: [ExactPoint; 3] = [
    (
        "7",
        "2",
        [
            ["61565/8648", "26/21", "73/21", "-1044257/181608"],
            ["35725/8648", "0", "0", "35725/8648"],
        ],
        [
            ["5", "8/3", "4/3", "8/3"],
            ["25/4", "-10/7", "15/7", "-30/7"],
        ],
    ),
    (
        "5",
        "5/2",
        [
            ["6205/712", "187237/37380", "-41165/7476", "282607/74760"],
            ["1805/712", "1805/2492", "5415/2492", "5415/4984"],
        ],
        [["5", "12/5", "0", "16/5"], ["25/4", "10/3", "-10/3", "5/3"]],
    ),
    (
        "11",
        "3/2",
        [
            ["66815/10952", "450423/76664", "-1829/9583", "9/14"],
            ["42705/10952", "-25623/10952", "8541/2738", "0"],
        ],
        [
            ["15/4", "-3/4", "3/2", "-3/2"],
            ["25/4", "30/7", "10/7", "15/7"],
        ],
    ),
];

/// tau- tau+ -> Z Z. Only the Higgs ghost vertex V_34 enters, which every
/// ghost-slot convention shares; this is a control of the kinematics.
const TAU_PAIR_TO_Z_PAIR: [ExactPoint; 1] = [(
    "7",
    "7/4",
    [["25/4", "0", "0", "6"], ["25/4", "0", "0", "-6"]],
    [
        ["25/4", "10/3", "10/3", "5/3"],
        ["25/4", "-10/3", "-10/3", "-5/3"],
    ],
)];

/// The tree-level BRST quartet identity in Feynman gauge: summed over every
/// covariant cut channel (V V with the -g numerator, V G, G V, G G and both
/// ghost-pair orientations), the forward cuts equal the V V channel with each
/// cut vector summed over its physical polarizations, -g + k k / M^2.
///
/// A massive tau couples the s-channel vector current to its Goldstone partner
/// through the Ward identity, so the ghost channels combine ghost-vector
/// (UUV1) with ghost-Goldstone (UUS1) vertices. The UUS1 particle lists fix
/// which slot holds the antighost (see docs/architecture/sm-conventions-audit.typ):
/// UUV1 = P(3,2) + P(3,3) satisfies the identity exactly, while P(3,2) alone
/// leaves a residual of 1.3-6.7% of the physical value at every W W and W Z
/// point below.
#[test]
fn massive_tau_covariant_cuts_satisfy_the_brst_quartet_identity() -> eyre::Result<()> {
    test_initialise()?;
    let model = load_generic_model("sm");
    let mut settings = GlobalSettings::default();
    settings.n_cores.feyngen = 1;
    let legs = |pdgs: [i64; 2], momenta: &[[&str; 4]; 2]| -> Vec<Leg> {
        pdgs.into_iter()
            .zip(momenta.map(|p| p.map(|component| parse!(component))))
            .collect()
    };

    for (initial_pdgs, final_pdgs, graph_count, ghost_orientations, points) in [
        ([15, -15], [24, -24], 114, 2, &TAU_PAIR_TO_W_PAIR[..]),
        ([15, -16], [-24, 23], 48, 2, &TAU_NEUTRINO_TO_W_Z[..]),
        ([15, -15], [23, 23], 20, 1, &TAU_PAIR_TO_Z_PAIR[..]),
    ] {
        let process = ProcessDefinition {
            generation_type: GenerationType::CrossSection,
            process: GenerationProcess::new(initial_pdgs, final_pdgs),
            generation_options: GenerationOptions::default()
                .with_loop_count(1, 1)?
                .with_graph_filter(GenerationFilter::CouplingOrders(
                    [("QED".into(), (4, Some(4)))].into_iter().collect(),
                )),
            ..Default::default()
        };
        let graphs = process.generate(&model, &settings)?;
        assert_eq!(graphs.len(), graph_count, "{final_pdgs:?}");
        let labels = process.covariant_cut_representatives(&model);
        let imported = Process::from_graph_list(
            "brst_quartet".into(),
            "born".into(),
            graphs,
            GenerationType::CrossSection,
            Some(process),
            None,
            &model,
        )?;
        let ProcessCollection::CrossSections(cross_sections) = &imported.collection else {
            panic!("forward graphs must import as a cross section");
        };
        let mut forward = Vec::new();
        for supergraph in &cross_sections["born"].supergraphs {
            // Spin-sum the initial states with the model's completeness
            // relations, one external edge at a time as in production.
            let mut graph = supergraph.graph.clone();
            let source = &supergraph.graph;
            for (_, edge, data) in source.iter_edges_of(&source.initial_state_cut) {
                let particle = data.data.particle().unwrap().resolve(&model);
                let Some(sum) = particle.polarization_sum(
                    &model,
                    edge,
                    false,
                    VectorPolarizationSumGauge::Feynman,
                )?
                else {
                    continue;
                };
                let projector = &graph.global_prefactor.projector;
                graph.global_prefactor.projector = projector
                    .with_checked_expression(projector.expression().replace_multiple(&[sum]))
                    .map_err(|error| eyre::eyre!("{error}"))?;
            }
            let cuts = supergraph.process_valid_cuts(
                &model,
                &imported.definition,
                &settings.generation,
            )?;
            forward.push((graph, cuts));
        }

        for (higgs_mass, tau_mass, incoming, outgoing) in points {
            let parameters = quartet_parameters(parse!(higgs_mass), parse!(tau_mass));
            let initial = legs(initial_pdgs, incoming);
            let finals = legs(final_pdgs, outgoing);
            for (pdg, momentum) in initial.iter().chain(&finals) {
                let mass = model.get_particle_from_pdg(*pdg as isize).resolve(&model);
                let mass = substitute(&mass.symbolic_mass(&model), &parameters);
                assert_eq!(minkowski_square(momentum), mass.pow(2), "{pdg} on shell");
            }
            // Identical final states average both momentum assignments.
            let assignments = if final_pdgs[0] == final_pdgs[1] {
                vec![finals.clone(), finals.iter().rev().cloned().collect()]
            } else {
                vec![finals]
            };
            let weight = Atom::num((1, assignments.len() as i64));
            let mut channels = BTreeMap::<Vec<i64>, Atom>::new();
            let mut physical = Atom::Zero;
            for (graph, cuts) in &forward {
                for (cut, finals) in cuts.iter().cartesian_product(&assignments) {
                    let kinematics = CutKinematics::new(
                        graph,
                        &model,
                        cut,
                        &labels,
                        (&initial, finals),
                        &parameters,
                    );
                    let (covariant, polarized) = kinematics.evaluate(graph, &model, &parameters)?;
                    *channels.entry(kinematics.channel).or_default() += &weight * covariant;
                    if let Some(polarized) = polarized {
                        physical += &weight * polarized;
                    }
                }
            }
            let covariant = channels
                .values()
                .fold(Atom::Zero, |sum, value| sum + value)
                .expand();
            let physical = physical.expand();
            assert_eq!(
                covariant, physical,
                "{final_pdgs:?}, MH={higgs_mass}: {channels:?}"
            );
            // The V V channel alone is the incomplete covariant metric sum, and
            // both ghost-pair orientations contribute.
            let vector_channel = final_pdgs.into_iter().sorted().collect_vec();
            assert_ne!(channels[&vector_channel], physical);
            let ghosts = channels
                .iter()
                .filter(|(channel, _)| {
                    model
                        .get_particle_from_pdg(channel[0] as isize)
                        .resolve(&model)
                        .is_ghost()
                })
                .collect_vec();
            assert_eq!(ghosts.len(), ghost_orientations, "{channels:?}");
            assert!(ghosts.iter().all(|(_, value)| !value.is_zero()));
        }
    }
    Ok(())
}

/// ee=sw=3/5, cw=4/5 and vev=6 obey MW=ee vev/(2 sw)=3 and MZ=MW/cw=15/4, the
/// tree-level relations the quartet identity requires. The Higgs self-coupling
/// and the tau Yukawa coupling follow from MH and m_tau.
fn quartet_parameters(higgs_mass: Atom, tau_mass: Atom) -> Vec<(&'static str, Atom)> {
    let vev = Atom::num(6);
    vec![
        ("ee", Atom::num((3, 5))),
        ("sw", Atom::num((3, 5))),
        ("cw", Atom::num((4, 5))),
        ("MW", Atom::num(3)),
        ("MZ", Atom::num((15, 4))),
        ("lam", higgs_mass.pow(2) / (Atom::num(2) * vev.pow(2))),
        ("MH", higgs_mass),
        // sqrt(2) = 2 (1/2)^(1/2) shares the radical of the model couplings, so
        // exact products of both reduce to rationals.
        (
            "ytau",
            Atom::num(2) * Atom::num((1, 2)).pow(Atom::num((1, 2))) * &tau_mass / &vev,
        ),
        ("ymtau", tau_mass.clone()),
        ("MTA", tau_mass),
        ("vev", vev),
    ]
}
