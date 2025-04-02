use std::sync::Arc;

use ahash::HashMap;
use ahash::HashMapExt;
use ahash::HashSet;
use insta::assert_snapshot;
use itertools::Itertools;
use symbolica::graph::Graph;

use crate::feyngen::diagram_generator::EdgeColor;
use crate::graph::BareGraph;
use crate::model::ColorStructure;
use crate::model::Model;
use crate::model::VertexRule;
use crate::numerator::GlobalPrefactor;
use crate::tests_from_pytest::load_generic_model;

use super::diagram_generator::NodeColorWithVertexRule;
use super::GenerationType;
use super::NumeratorAwareGraphGroupingOption;
use super::SewedFilterOptions;
use super::SnailFilterOptions;
use super::TadpolesFilterOptions;
use super::{diagram_generator::FeynGen, FeynGenFilter, FeynGenFilters, FeynGenOptions};

#[test]
fn cut_content() {
    let model = load_generic_model("sm");

    let mut coupling = HashMap::new();
    coupling.insert("QED".into(), (6, Some(6)));
    let mut pert = HashMap::new();
    pert.insert("QCD".into(), 1);
    let filters = FeynGen::new(FeynGenOptions {
        generation_type: GenerationType::CrossSection,
        initial_pdgs: vec![-11, 11],
        final_pdgs_lists: vec![vec![5, -5, 25]],
        loop_count_range: (3, 3),
        symmetrize_initial_states: true,
        symmetrize_final_states: true,
        symmetrize_left_right_states: true,
        allow_symmetrization_of_external_fermions_in_amplitudes: false,
        max_multiplicity_for_fast_cut_filter: 6,
        amplitude_filters: FeynGenFilters(vec![
            FeynGenFilter::SelfEnergyFilter(Default::default()),
            FeynGenFilter::ParticleVeto(vec![
                23, 24, 9000001, 9000002, 9000003, 9000004, 12, 14, 16, 2, 4, 6, 3, 250, 251, 13,
                15,
            ]),
            FeynGenFilter::TadpolesFilter(Default::default()),
            FeynGenFilter::ZeroSnailsFilter(Default::default()),
        ]),
        cross_section_filters: FeynGenFilters(vec![
            FeynGenFilter::SelfEnergyFilter(Default::default()),
            FeynGenFilter::ParticleVeto(vec![
                23, 24, 9000001, 9000002, 9000003, 9000004, 12, 14, 16, 2, 4, 6, 3, 250, 251, 13,
                15,
            ]),
            FeynGenFilter::TadpolesFilter(Default::default()),
            FeynGenFilter::ZeroSnailsFilter(Default::default()),
            FeynGenFilter::PerturbativeOrders(pert),
            FeynGenFilter::CouplingOrders(coupling),
            FeynGenFilter::LoopCountRange((3, 3)),
            FeynGenFilter::BlobRange(1..=1),
            FeynGenFilter::SpectatorRange(0..=0),
        ]),
    });

    let mut graph = Graph::new();
    #[allow(non_snake_case)]
    let bbH = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule("V_78"),
    };
    let bba = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule("V_73"),
    };
    let bbg = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule("V_76"),
    };
    let epema = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule("V_98"),
    };

    let dummy_external_vertex_rule = Arc::new(VertexRule {
        name: "external".into(),
        couplings: vec![],
        lorentz_structures: vec![],
        particles: vec![],
        color_structures: ColorStructure::new(vec![]),
    });

    let v0 = graph.add_node(bbH.clone());
    let v1 = graph.add_node(bbH.clone());
    let v2 = graph.add_node(bbg.clone());
    let v3 = graph.add_node(bba.clone());
    let v4 = graph.add_node(bbg.clone());
    let v5 = graph.add_node(bba.clone());
    let v6 = graph.add_node(epema.clone());
    let v7 = graph.add_node(epema.clone());
    let e1 = NodeColorWithVertexRule {
        external_tag: 1,
        vertex_rule: dummy_external_vertex_rule.clone(),
    };
    let mut e2 = e1.clone();
    e2.external_tag = 2;
    let mut e3 = e1.clone();
    e3.external_tag = 3;
    let mut e4 = e1.clone();
    e4.external_tag = 4;

    let v8 = graph.add_node(e1.clone());
    let v9 = graph.add_node(e2.clone());
    let v10 = graph.add_node(e3.clone());
    let v11 = graph.add_node(e4.clone());

    let h = EdgeColor::from_particle(model.get_particle("H"));

    let b = EdgeColor::from_particle(model.get_particle("b"));
    let bbar = EdgeColor::from_particle(model.get_particle("b").get_anti_particle(&model));
    let g = EdgeColor::from_particle(model.get_particle("g"));
    let a = EdgeColor::from_particle(model.get_particle("a"));
    let eplus = EdgeColor::from_particle(model.get_particle("e+"));
    let eminus = EdgeColor::from_particle(model.get_particle("e-"));
    graph.add_edge(v0, v1, true, h).unwrap();
    graph.add_edge(v0, v5, true, b).unwrap();
    graph.add_edge(v1, v3, true, b).unwrap();
    graph.add_edge(v2, v4, true, g).unwrap();
    graph.add_edge(v2, v4, true, b).unwrap();
    graph.add_edge(v3, v2, true, b).unwrap();
    graph.add_edge(v3, v7, true, a).unwrap();
    graph.add_edge(v4, v1, true, b).unwrap();
    graph.add_edge(v5, v0, true, b).unwrap();
    graph.add_edge(v5, v6, true, a).unwrap();
    graph.add_edge(v6, v10, false, eplus).unwrap();
    graph.add_edge(v6, v11, false, eminus).unwrap();
    graph.add_edge(v8, v7, false, eplus).unwrap();
    graph.add_edge(v9, v7, false, eminus).unwrap();

    let (n_unresolved, unresolved_type) = filters.unresolved_cut_content(&model);
    assert!(!filters.half_edge_filters(&model, &graph, &[], n_unresolved, &unresolved_type));

    let mut double_double_triangle = Graph::new();
    let v0 = double_double_triangle.add_node(e1.clone());
    let v1 = double_double_triangle.add_node(e2.clone());
    let v2 = double_double_triangle.add_node(e3.clone());
    let v3 = double_double_triangle.add_node(e4.clone());

    let v4 = double_double_triangle.add_node(bba.clone());
    let v5 = double_double_triangle.add_node(bba.clone());
    let v6 = double_double_triangle.add_node(bbH.clone());
    let v7 = double_double_triangle.add_node(bbH.clone());
    let v8 = double_double_triangle.add_node(bbg.clone());
    let v9 = double_double_triangle.add_node(bbg.clone());
    let v10 = double_double_triangle.add_node(bbg.clone());
    let v11 = double_double_triangle.add_node(bbg.clone());
    let v12 = double_double_triangle.add_node(epema.clone());
    let v13 = double_double_triangle.add_node(epema.clone());

    double_double_triangle
        .add_edge(v0, v13, true, eplus)
        .unwrap();
    double_double_triangle
        .add_edge(v1, v13, true, eminus)
        .unwrap();
    double_double_triangle.add_edge(v4, v13, false, a).unwrap();
    double_double_triangle
        .add_edge(v4, v11, true, bbar)
        .unwrap();
    double_double_triangle.add_edge(v4, v10, true, b).unwrap();
    double_double_triangle.add_edge(v10, v11, false, g).unwrap();
    double_double_triangle.add_edge(v6, v11, true, b).unwrap();
    double_double_triangle
        .add_edge(v6, v10, true, bbar)
        .unwrap();
    double_double_triangle.add_edge(v6, v7, true, h).unwrap();
    double_double_triangle.add_edge(v7, v8, true, bbar).unwrap();
    double_double_triangle.add_edge(v7, v9, true, b).unwrap();
    double_double_triangle.add_edge(v5, v8, true, b).unwrap();
    double_double_triangle.add_edge(v5, v9, true, bbar).unwrap();
    double_double_triangle.add_edge(v8, v9, true, g).unwrap();
    double_double_triangle.add_edge(v5, v12, true, a).unwrap();
    double_double_triangle
        .add_edge(v12, v2, true, eminus)
        .unwrap();
    double_double_triangle
        .add_edge(v12, v3, true, eplus)
        .unwrap();

    let mut coupling = HashMap::new();
    coupling.insert("QED".into(), (6, Some(6)));
    let mut pert = HashMap::new();
    pert.insert("QCD".into(), 2);
    let filters = FeynGen::new(FeynGenOptions {
        generation_type: GenerationType::CrossSection,
        initial_pdgs: vec![-11, 11],
        final_pdgs_lists: vec![vec![5, -5, 25]],
        loop_count_range: (4, 4),
        symmetrize_initial_states: true,
        symmetrize_final_states: true,
        symmetrize_left_right_states: true,
        allow_symmetrization_of_external_fermions_in_amplitudes: false,
        max_multiplicity_for_fast_cut_filter: 6,
        amplitude_filters: FeynGenFilters(vec![
            FeynGenFilter::SelfEnergyFilter(Default::default()),
            FeynGenFilter::ParticleVeto(vec![
                23, 24, 9000001, 9000002, 9000003, 9000004, 12, 14, 16, 2, 4, 6, 3, 250, 251, 13,
                15,
            ]),
            FeynGenFilter::TadpolesFilter(Default::default()),
            FeynGenFilter::ZeroSnailsFilter(Default::default()),
        ]),
        cross_section_filters: FeynGenFilters(vec![
            FeynGenFilter::SelfEnergyFilter(Default::default()),
            FeynGenFilter::ParticleVeto(vec![
                23, 24, 9000001, 9000002, 9000003, 9000004, 12, 14, 16, 2, 4, 6, 3, 250, 251, 13,
                15,
            ]),
            FeynGenFilter::TadpolesFilter(Default::default()),
            FeynGenFilter::ZeroSnailsFilter(Default::default()),
            FeynGenFilter::PerturbativeOrders(pert),
            FeynGenFilter::CouplingOrders(coupling),
            FeynGenFilter::LoopCountRange((4, 4)),
            FeynGenFilter::BlobRange(1..=1),
            FeynGenFilter::SpectatorRange(0..=0),
        ]),
    });

    let (n_unresolved, unresolved_type) = filters.unresolved_cut_content(&model);
    assert!(!filters.half_edge_filters(
        &model,
        &double_double_triangle,
        &[],
        n_unresolved,
        &unresolved_type
    ));
}

// fn fs_generate(loop_count: usize) -> Vec<BareGraph> {
//     let model = load_generic_model("sm_dis");

//     let options = dis_cart_prod(&["d", "d~", "g"], loop_count, &model);
//     chain_dis_generate(&options, &model)
// }

pub fn dis_options_impl(
    init: &[isize],
    final_states: &[Vec<isize>],
    pert: usize,
    loop_count: usize,
    coupling: usize,
) -> FeynGen {
    let mut amp_coupling = HashMap::new();
    amp_coupling.insert("QED".into(), (2, Some(2)));
    amp_coupling.insert("LT".into(), (1, Some(1)));
    let mut xs_coupling = HashMap::new();
    xs_coupling.insert("QCD".into(), (coupling, Some(coupling)));

    let mut xs_pert = HashMap::new();
    xs_pert.insert("QCD".into(), pert);
    FeynGen {
        options: FeynGenOptions {
            generation_type: GenerationType::CrossSection,
            initial_pdgs: init.iter().map(|a| *a as i64).collect(),
            final_pdgs_lists: final_states
                .iter()
                .map(|f| f.iter().map(|a| *a as i64).collect())
                .collect(),
            loop_count_range: (loop_count, loop_count),
            symmetrize_initial_states: true,
            symmetrize_final_states: true,
            symmetrize_left_right_states: false,
            allow_symmetrization_of_external_fermions_in_amplitudes: false,
            max_multiplicity_for_fast_cut_filter: 0,
            amplitude_filters: FeynGenFilters(vec![
                FeynGenFilter::ParticleVeto(vec![
                    23, 24, 9000001, 9000002, 9000003, 9000004, 9000005, -9000005, 12, 14, 16, 2,
                    4, 6, 3, 5, 25, 250, 251, 13, 15,
                ]),
                FeynGenFilter::TadpolesFilter(TadpolesFilterOptions {
                    veto_tadpoles_attached_to_massive_lines: true,
                    veto_tadpoles_attached_to_massless_lines: true,
                    veto_only_scaleless_tadpoles: false,
                }),
                FeynGenFilter::ZeroSnailsFilter(SnailFilterOptions {
                    veto_snails_attached_to_massive_lines: false,
                    veto_snails_attached_to_massless_lines: true,
                    veto_only_scaleless_snails: false,
                }),
                FeynGenFilter::CouplingOrders(amp_coupling),
            ]),
            cross_section_filters: FeynGenFilters(vec![
                FeynGenFilter::SewedFilter(SewedFilterOptions {
                    filter_tadpoles: true,
                }),
                FeynGenFilter::ParticleVeto(vec![
                    23, 24, 9000001, 9000002, 9000003, 9000004, 9000005, -9000005, 12, 14, 16, 2,
                    4, 6, 3, 5, 25, 250, 251, 13, 15,
                ]),
                FeynGenFilter::TadpolesFilter(TadpolesFilterOptions {
                    veto_tadpoles_attached_to_massive_lines: true,
                    veto_tadpoles_attached_to_massless_lines: true,
                    veto_only_scaleless_tadpoles: false,
                }),
                FeynGenFilter::ZeroSnailsFilter(SnailFilterOptions {
                    veto_snails_attached_to_massive_lines: false,
                    veto_snails_attached_to_massless_lines: true,
                    veto_only_scaleless_snails: false,
                }),
                FeynGenFilter::PerturbativeOrders(xs_pert),
                FeynGenFilter::CouplingOrders(xs_coupling),
                FeynGenFilter::LoopCountRange((loop_count, loop_count)),
                FeynGenFilter::BlobRange(1..=1),
                FeynGenFilter::SpectatorRange(0..=1),
            ]),
        },
    }
}

pub fn dis_options(
    init: &[&'static str],
    final_states: &[Vec<&'static str>],
    pert: usize,
    model: &Model,
) -> FeynGen {
    let initial_pdgs: Vec<_> = init
        .iter()
        .map(|a| model.get_particle(a).pdg_code)
        .collect();

    let final_pdgs_lists: Vec<_> = final_states
        .iter()
        .map(|f| f.iter().map(|a| model.get_particle(a).pdg_code).collect())
        .collect();

    let initial_state_mult = init.len();
    let loop_count = pert + 3 - initial_state_mult;
    let coupling = 2 * pert;
    dis_options_impl(&initial_pdgs, &final_pdgs_lists, pert, loop_count, coupling)
}
struct CombinationsWithRepetition<T> {
    elements: Vec<T>,
    k: usize,
    current: Vec<usize>,
    done: bool,
}

impl<T: Clone + Ord> CombinationsWithRepetition<T> {
    fn new(set: HashSet<T>, k: usize) -> Self {
        let elements: Vec<T> = set.into_iter().sorted().collect();
        CombinationsWithRepetition {
            elements,
            k,
            current: vec![0; k],
            done: false,
        }
    }
}

impl<T: Clone> Iterator for CombinationsWithRepetition<T> {
    type Item = Vec<T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        let result = self
            .current
            .iter()
            .map(|&i| self.elements[i].clone())
            .collect();

        // Generate next combination (non-decreasing sequence)
        let n = self.elements.len();
        let mut i = self.k - 1;
        while i < self.k && self.current[i] == n - 1 {
            i = i.wrapping_sub(1);
        }
        if i >= self.k {
            self.done = true;
        } else {
            self.current[i] += 1;
            for j in (i + 1)..self.k {
                self.current[j] = self.current[i];
            }
        }

        Some(result)
    }
}

pub fn dis_cart_prod_impl(initial_states: &[isize], loop_count: usize) -> Vec<FeynGen> {
    let mut options = vec![];

    let initial_states: HashSet<isize> = initial_states.iter().cloned().collect();

    let initial_state_template = vec![11];
    let final_states = initial_states
        .iter()
        .map(|a| {
            let mut temp = initial_state_template.iter().copied().collect_vec();

            temp.extend([*a]);
            temp
        })
        .collect_vec();

    for initial_state_mult in 1..(loop_count + 2) {
        for mut init_states in
            CombinationsWithRepetition::new(initial_states.clone(), initial_state_mult)
        {
            init_states.extend(initial_state_template.clone());

            let init_states = init_states.into_iter().collect_vec();

            // info!("initial states: {:?}\nfinal states: {:?}\ncross_section_orders:{}\nloop_count:{}\nn_unresolved:{}", init_states, final_states,2*loop_count ,loop_count+ 2 - initial_state_mult,loop_count);

            let pert = loop_count;
            let loop_count = pert + 2 - initial_state_mult;
            let coupling = 2 * pert;

            options.push(dis_options_impl(
                &init_states,
                &final_states,
                pert,
                loop_count,
                coupling,
            ));
        }
    }

    options
}

pub fn dis_cart_prod(
    initial_states: &[&'static str],
    loop_count: usize,
    model: &Model,
) -> Vec<FeynGen> {
    let initial_states: Vec<_> = initial_states
        .iter()
        .map(|a| model.get_particle(a).pdg_code)
        .collect();
    dis_cart_prod_impl(&initial_states, loop_count)
}

pub fn chain_dis_generate(options: &[FeynGen], model: &Model) -> Vec<BareGraph> {
    options
        .iter()
        .flat_map(|a| {
            a.generate(
                model,
                &NumeratorAwareGraphGroupingOption::OnlyDetectZeroes,
                true,
                "DIS".into(),
                None,
                None,
                None,
                GlobalPrefactor::default(),
                Some(10),
            )
            .unwrap()
        })
        .collect()
}

#[test]
fn nlo_fs_dis() {
    // setup_logger().unwrap();

    let model = load_generic_model("sm_dis");
    let options = dis_cart_prod(&["d", "d~", "g"], 1, &model);

    for option in options {
        let diagrams = chain_dis_generate(&[option.clone()], &model);

        let process_name = option
            .options
            .initial_pdgs
            .iter()
            .map(|a| model.get_particle_from_pdg(*a as isize).name.clone())
            .join(",");

        assert_snapshot!(
            format!("number of diagrams for {}", process_name),
            diagrams.len()
        );
        let dots = diagrams
            .iter()
            .map(|a| a.dot())
            .collect::<Vec<_>>()
            .join("\n");
        assert_snapshot!(format!("dots for {}", process_name), dots);
    }

    // println!("Number of fs: {}", fs_diagrams.len());
}
