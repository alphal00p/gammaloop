//! Generation workloads shared by the oracle tests, the search census and projector equivalence
//! unit tests and the `generation_workloads` example that the instruction gate
//! (`bin/generation_ir.py`) runs.
//!
//! Names read `<model>_<process>_L<loops>`. The profile cases of the generation redesign keep
//! their names; `phi*` cases have independent oracles (OEIS A005967 and zero-dimensional field
//! theory), `xs_*` cases are cross sections with physical cuts, `ct_*` models duplicate vertex
//! rules with identical field content, as split UFO vertices and counterterms do,
//! `identical`, `sign` and `scalar` group numerators that are identical, equal up to a sign or
//! proportional, `sym` symmetrizes identical external states and `lr` identifies cross-section
//! mirror images.

use std::collections::BTreeMap;

use feynkit_generator::{
    GenerationError, GenerationFilter, GenerationOptions, GenerationResult, GenerationType,
    GraphGroupingOptions, NumeratorGrouping, ParticleSelector, Process,
};
use feynkit_model::{Model, Particle};

type Configure = fn(&Model, GenerationOptions) -> GenerationOptions;

pub struct Case {
    pub name: &'static str,
    pub kind: GenerationType,
    model: fn() -> Model,
    incoming: &'static [&'static str],
    outgoing: &'static [&'static str],
    loops: usize,
    configure: Configure,
}

/// The reference workloads of the instruction gate: gg->gg L3, gg->ggg L2, uu~->gg L3, the
/// photon self-energy at five loops, ee->4f L2, SM ee->mumu L2, SM uu~->dd~ at QED<=2 L2 and
/// gg->6g.
pub const REFERENCE: [&str; 8] = [
    "qcd_gggg_L3",
    "qcd_ggggg_L2",
    "qcd_uugg_L3",
    "qed_ase_L5",
    "qed_ee4f_L2",
    "sm_eemm_L2",
    "sm_uudd_qed2_L2",
    "ym_g8_L0",
];

impl Case {
    pub fn named(name: &str) -> &'static Case {
        CASES
            .iter()
            .find(|case| case.name == name)
            .unwrap_or_else(|| panic!("unknown generation case {name}"))
    }

    pub fn model(&self) -> Model {
        (self.model)()
    }

    pub fn process(&self) -> Process {
        Process::new(self.incoming.iter().copied(), self.outgoing.iter().copied())
    }

    /// The case's options with an explicit thread count.
    pub fn options(&self, model: &Model, threads: usize) -> GenerationOptions {
        let options = GenerationOptions::default()
            .with_loop_count(self.loops, self.loops)
            .unwrap()
            .threads(threads);
        (self.configure)(model, options)
    }

    pub fn generate_with(
        &self,
        model: Model,
        options: &GenerationOptions,
    ) -> Result<GenerationResult, GenerationError> {
        match self.kind {
            GenerationType::Amplitude => self.process().generate_diagrams(model, options),
            GenerationType::CrossSection => self.process().generate_cross_section(model, options),
        }
    }
}

const fn amplitude(
    name: &'static str,
    model: fn() -> Model,
    incoming: &'static [&'static str],
    outgoing: &'static [&'static str],
    loops: usize,
    configure: Configure,
) -> Case {
    Case {
        name,
        kind: GenerationType::Amplitude,
        model,
        incoming,
        outgoing,
        loops,
        configure,
    }
}

/// A cross section: the forward graphs with physical cuts. `loops` counts forward-graph loops.
const fn cross_section(
    name: &'static str,
    model: fn() -> Model,
    incoming: &'static [&'static str],
    outgoing: &'static [&'static str],
    loops: usize,
    configure: Configure,
) -> Case {
    Case {
        kind: GenerationType::CrossSection,
        ..amplitude(name, model, incoming, outgoing, loops, configure)
    }
}

fn plain(_: &Model, options: GenerationOptions) -> GenerationOptions {
    options
}

fn one_pi(_: &Model, options: GenerationOptions) -> GenerationOptions {
    options.with_graph_filter(GenerationFilter::MaxNumberOfBridges(0))
}

fn tadpoles(_: &Model, options: GenerationOptions) -> GenerationOptions {
    options.allow_self_loops(true).allow_zero_flow_edges(true)
}

fn one_pi_tadpoles(model: &Model, options: GenerationOptions) -> GenerationOptions {
    one_pi(model, tadpoles(model, options))
}

fn loop_one_pi(_: &Model, options: GenerationOptions) -> GenerationOptions {
    options.with_graph_filter(GenerationFilter::LoopOneParticleIrreducible)
}

fn order_at_most(options: GenerationOptions, order: &str, maximum: usize) -> GenerationOptions {
    options.with_graph_filter(GenerationFilter::CouplingOrders(BTreeMap::from([(
        order.to_owned(),
        (0, Some(maximum)),
    )])))
}

fn no_qcd(_: &Model, options: GenerationOptions) -> GenerationOptions {
    order_at_most(options, "QCD", 0)
}

fn qed_at_most_two(_: &Model, options: GenerationOptions) -> GenerationOptions {
    order_at_most(options, "QED", 2)
}

fn orders(
    options: GenerationOptions,
    required: &[(&str, usize, Option<usize>)],
) -> GenerationOptions {
    options.with_graph_filter(GenerationFilter::CouplingOrders(
        required
            .iter()
            .map(|&(order, minimum, maximum)| (order.to_owned(), (minimum, maximum)))
            .collect(),
    ))
}

/// Mixed QCD and electroweak corrections: both orders exactly two.
fn qcd2_qed2(_: &Model, options: GenerationOptions) -> GenerationOptions {
    orders(options, &[("QCD", 2, Some(2)), ("QED", 2, Some(2))])
}

/// A lower bound only.
fn qcd_at_least_three(_: &Model, options: GenerationOptions) -> GenerationOptions {
    orders(options, &[("QCD", 3, None)])
}

fn one_counterterm(_: &Model, options: GenerationOptions) -> GenerationOptions {
    orders(options, &[("CT", 1, Some(1))])
}

fn two_counterterms(model: &Model, options: GenerationOptions) -> GenerationOptions {
    orders(tadpoles(model, options), &[("CT", 2, Some(2))])
}

fn at_most_one_counterterm(model: &Model, options: GenerationOptions) -> GenerationOptions {
    order_at_most(tadpoles(model, options), "CT", 1)
}

fn fermion_loops(options: GenerationOptions, minimum: usize, maximum: usize) -> GenerationOptions {
    options.with_graph_filter(GenerationFilter::FermionLoopCountRange((minimum, maximum)))
}

fn no_fermion_loops(_: &Model, options: GenerationOptions) -> GenerationOptions {
    fermion_loops(options, 0, 0)
}

fn one_fermion_loop(_: &Model, options: GenerationOptions) -> GenerationOptions {
    fermion_loops(options, 1, 1)
}

fn two_fermion_loops_one_pi(model: &Model, options: GenerationOptions) -> GenerationOptions {
    fermion_loops(one_pi(model, options), 2, 2)
}

/// QED order four at tree level and one loop: trees have order two, so the tree level is proven
/// empty and its topologies are not enumerated.
fn qed4_up_to_one_loop(_: &Model, options: GenerationOptions) -> GenerationOptions {
    let options = options.with_loop_count(0, 1).unwrap();
    orders(options, &[("QED", 4, Some(4))])
}

fn symmetrized(_: &Model, options: GenerationOptions) -> GenerationOptions {
    options.symmetrize_initial(true).symmetrize_final(true)
}

/// The exchange of identical beams; cross sections reject final-state symmetrization.
fn beam_exchange(_: &Model, options: GenerationOptions) -> GenerationOptions {
    options.symmetrize_initial(true)
}

/// Cross-section forward graphs identified with their mirror images.
fn left_right(_: &Model, options: GenerationOptions) -> GenerationOptions {
    options.symmetrize_left_right(true)
}

fn beam_exchange_left_right(model: &Model, options: GenerationOptions) -> GenerationOptions {
    left_right(model, beam_exchange(model, options))
}

fn identical_numerators(_: &Model, options: GenerationOptions) -> GenerationOptions {
    options.numerator_grouping(NumeratorGrouping::Identical(GraphGroupingOptions::default()))
}

fn numerators_up_to_sign(_: &Model, options: GenerationOptions) -> GenerationOptions {
    options.numerator_grouping(NumeratorGrouping::UpToSign(GraphGroupingOptions::default()))
}

fn numerators_up_to_scalar(_: &Model, options: GenerationOptions) -> GenerationOptions {
    options.numerator_grouping(NumeratorGrouping::UpToScalar(
        GraphGroupingOptions::default(),
    ))
}

/// The `rejected` particles rejected in flight by a user filter, as the Python dressing callback
/// does.
fn filtered(
    model: &Model,
    options: GenerationOptions,
    rejected: impl Fn(&Particle) -> bool,
) -> GenerationOptions {
    let rejected: Vec<_> = model
        .particles()
        .iter()
        .filter(|particle| rejected(particle))
        .map(|particle| model.particle_id(&particle.name).unwrap())
        .collect();
    options.prune_dressings(move |view| {
        view.edges()
            .all(|edge| edge.stub.is_none_or(|stub| !rejected.contains(stub.field)))
    })
}

fn ghosts_filtered(model: &Model, options: GenerationOptions) -> GenerationOptions {
    filtered(model, options, Particle::is_ghost)
}

/// One fermion loop, with the quarks of PDG code `pdg` rejected in flight by a user filter: it
/// rejects field assignments of the fermions that the fermion-loop budget follows.
fn one_fermion_loop_without(
    model: &Model,
    options: GenerationOptions,
    pdg: i64,
) -> GenerationOptions {
    filtered(model, fermion_loops(options, 1, 1), |particle| {
        particle.pdg_code.abs() == pdg
    })
}

fn one_fermion_loop_without_down_quarks(
    model: &Model,
    options: GenerationOptions,
) -> GenerationOptions {
    one_fermion_loop_without(model, options, 1)
}

fn one_fermion_loop_without_top_quarks(
    model: &Model,
    options: GenerationOptions,
) -> GenerationOptions {
    one_fermion_loop_without(model, options, 6)
}

/// Ghost fields removed by the process's particle veto.
fn ghosts_vetoed(model: &Model, options: GenerationOptions) -> GenerationOptions {
    let ghosts = model
        .particles()
        .iter()
        .filter(|particle| particle.is_ghost())
        .map(|particle| ParticleSelector::from(particle.name.as_str()))
        .collect();
    options.with_graph_filter(GenerationFilter::ParticleVeto(ghosts))
}

/// `model` with a copy of every vertex rule `duplicate` selects by its field content, named with
/// a `_CT` suffix.
fn with_counterterms(model: Model, duplicate: impl Fn(&[String]) -> bool) -> Model {
    let mut definition = serde_json::to_value(model).unwrap();
    let rules = definition["vertex_rules"].as_array_mut().unwrap();
    let copies: Vec<_> = rules
        .iter()
        .filter(|rule| {
            let particles: Vec<String> = serde_json::from_value(rule["particles"].clone()).unwrap();
            duplicate(&particles)
        })
        .map(|rule| {
            let mut copy = rule.clone();
            copy["name"] = format!("{}_CT", rule["name"].as_str().unwrap()).into();
            copy
        })
        .collect();
    assert!(!copies.is_empty(), "no vertex rule to duplicate");
    rules.extend(copies);
    serde_json::from_value(definition).unwrap()
}

/// phi^3 with a counterterm copy of its vertex.
fn ct_phi3() -> Model {
    with_counterterms(Model::phi3(), |_| true)
}

/// QCD with counterterm copies of the three-gluon and the up-quark-gluon vertices.
fn ct_qcd() -> Model {
    with_counterterms(Model::qcd(), qcd_counterterm)
}

fn qcd_counterterm(particles: &[String]) -> bool {
    let mut particles = particles.to_vec();
    particles.sort();
    particles == ["g", "g", "g"] || particles == ["g", "u", "u~"]
}

/// `model` with the counterterm copies of [`with_counterterms`], whose couplings carry one power
/// of an extra `CT` order: a copy and its rule have the same fields but different orders.
fn with_counterterm_order(model: Model, duplicate: impl Fn(&[String]) -> bool) -> Model {
    let mut definition = serde_json::to_value(with_counterterms(model, duplicate)).unwrap();
    let couplings = definition["couplings"].as_array().unwrap().clone();
    let mut copies: Vec<serde_json::Value> = vec![];
    for rule in definition["vertex_rules"].as_array_mut().unwrap() {
        if !rule["name"].as_str().unwrap().ends_with("_CT") {
            continue;
        }
        for entry in rule["couplings"]
            .as_array_mut()
            .unwrap()
            .iter_mut()
            .flat_map(|row| row.as_array_mut().unwrap())
        {
            let Some(name) = entry.as_str().map(str::to_owned) else {
                continue;
            };
            let copy = format!("{name}_CT");
            if !copies
                .iter()
                .any(|coupling| coupling["name"] == copy.as_str())
            {
                let mut coupling = couplings
                    .iter()
                    .find(|coupling| coupling["name"] == name.as_str())
                    .unwrap()
                    .clone();
                coupling["name"] = copy.clone().into();
                let orders = coupling["orders"].as_array_mut().unwrap();
                orders.push(serde_json::json!(["CT", 1]));
                copies.push(coupling);
            }
            *entry = copy.into();
        }
    }
    definition["couplings"]
        .as_array_mut()
        .unwrap()
        .extend(copies);
    definition["orders"]
        .as_array_mut()
        .unwrap()
        .push(serde_json::json!({"name": "CT", "expansion_order": 99, "hierarchy": 3}));
    serde_json::from_value(definition).unwrap()
}

/// phi^3 whose counterterm copy carries a `CT` order.
fn ct_order_phi3() -> Model {
    with_counterterm_order(Model::phi3(), |_| true)
}

/// QCD whose counterterm copies of the three-gluon and the up-quark-gluon vertices carry a `CT`
/// order.
fn ct_order_qcd() -> Model {
    with_counterterm_order(Model::qcd(), qcd_counterterm)
}

const EE: &[&str] = &["e-", "e+"];
const MUMU: &[&str] = &["mu-", "mu+"];
const FOUR_LEPTONS: &[&str] = &["mu-", "mu+", "ta-", "ta+"];
const GG: &[&str] = &["g", "g"];
const UU: &[&str] = &["u", "u~"];
const DD: &[&str] = &["d", "d~"];
const A: &[&str] = &["a"];
const G: &[&str] = &["g"];
const PHI: &[&str] = &["phi"];
const NONE: &[&str] = &[];

pub const CASES: &[Case] = &[
    // Scalar theories with independent oracles.
    amplitude("phi3_vac_L2", Model::phi3, NONE, NONE, 2, tadpoles),
    amplitude("phi3_vac_L3", Model::phi3, NONE, NONE, 3, tadpoles),
    amplitude("phi3_vac_L4", Model::phi3, NONE, NONE, 4, tadpoles),
    amplitude("phi3_vac_L5", Model::phi3, NONE, NONE, 5, tadpoles),
    amplitude("phi3_vac_L6", Model::phi3, NONE, NONE, 6, tadpoles),
    amplitude("phi3_vac_L7", Model::phi3, NONE, NONE, 7, tadpoles),
    amplitude("phi4_vac_L2", Model::phi4, NONE, NONE, 2, tadpoles),
    amplitude("phi4_vac_L3", Model::phi4, NONE, NONE, 3, tadpoles),
    amplitude("phi4_vac_L4", Model::phi4, NONE, NONE, 4, tadpoles),
    amplitude("phi4_vac_L5", Model::phi4, NONE, NONE, 5, tadpoles),
    amplitude("phi4_vac_L6", Model::phi4, NONE, NONE, 6, tadpoles),
    amplitude("phi4_vac_L7", Model::phi4, NONE, NONE, 7, tadpoles),
    amplitude("phi34_2pt_L1", Model::phi_3_4, PHI, PHI, 1, tadpoles),
    amplitude("phi34_2pt_L2", Model::phi_3_4, PHI, PHI, 2, tadpoles),
    amplitude("phi34_2pt_L3", Model::phi_3_4, PHI, PHI, 3, tadpoles),
    amplitude("phi34_2pt_L4", Model::phi_3_4, PHI, PHI, 4, tadpoles),
    amplitude("phi34_2pt_L5", Model::phi_3_4, PHI, PHI, 5, tadpoles),
    amplitude("phi3_2pt_L1", Model::phi3, PHI, PHI, 1, plain),
    amplitude("phi3_2pt_tad_L2", Model::phi3, PHI, PHI, 2, tadpoles),
    amplitude(
        "phi34_2pt_identical_L3",
        Model::phi_3_4,
        PHI,
        PHI,
        3,
        identical_numerators,
    ),
    // QED: the photon and all charged fermions.
    amplitude("qed_eemm_L1", Model::qed, EE, MUMU, 1, plain),
    amplitude("qed_eemm_L2", Model::qed, EE, MUMU, 2, plain),
    amplitude("qed_eemm_L3", Model::qed, EE, MUMU, 3, plain),
    amplitude("qed_ee4f_L0", Model::qed, EE, FOUR_LEPTONS, 0, plain),
    amplitude("qed_ee4f_L1", Model::qed, EE, FOUR_LEPTONS, 1, plain),
    amplitude("qed_ee4f_L2", Model::qed, EE, FOUR_LEPTONS, 2, plain),
    amplitude("qed_eeaaa_L1", Model::qed, EE, &["a", "a", "a"], 1, plain),
    amplitude("qed_eeaaa_L2", Model::qed, EE, &["a", "a", "a"], 2, plain),
    amplitude(
        "qed_eeaa_identical_L1",
        Model::qed,
        EE,
        &["a", "a"],
        1,
        identical_numerators,
    ),
    amplitude(
        "qed_eeaa_sign_L1",
        Model::qed,
        EE,
        &["a", "a"],
        1,
        numerators_up_to_sign,
    ),
    amplitude(
        "qed_eeaa_scalar_L1",
        Model::qed,
        EE,
        &["a", "a"],
        1,
        numerators_up_to_scalar,
    ),
    amplitude(
        "qed_ee4f_identical_L1",
        Model::qed,
        EE,
        FOUR_LEPTONS,
        1,
        identical_numerators,
    ),
    amplitude(
        "qed_ee4f_sign_L1",
        Model::qed,
        EE,
        FOUR_LEPTONS,
        1,
        numerators_up_to_sign,
    ),
    amplitude(
        "qed_ee4f_scalar_L1",
        Model::qed,
        EE,
        FOUR_LEPTONS,
        1,
        numerators_up_to_scalar,
    ),
    amplitude(
        "qed_eemm_identical_L2",
        Model::qed,
        EE,
        MUMU,
        2,
        identical_numerators,
    ),
    amplitude(
        "qed_eemm_sign_L2",
        Model::qed,
        EE,
        MUMU,
        2,
        numerators_up_to_sign,
    ),
    amplitude(
        "qed_eemm_scalar_L2",
        Model::qed,
        EE,
        MUMU,
        2,
        numerators_up_to_scalar,
    ),
    amplitude(
        "qed_eeaa_identical_L2",
        Model::qed,
        EE,
        &["a", "a"],
        2,
        identical_numerators,
    ),
    amplitude(
        "qed_eeaa_sign_L2",
        Model::qed,
        EE,
        &["a", "a"],
        2,
        numerators_up_to_sign,
    ),
    amplitude(
        "qed_eeaa_scalar_L2",
        Model::qed,
        EE,
        &["a", "a"],
        2,
        numerators_up_to_scalar,
    ),
    amplitude("qed_ase_L3", Model::qed, A, A, 3, one_pi),
    amplitude("qed_ase_L4", Model::qed, A, A, 4, one_pi),
    amplitude("qed_ase_L5", Model::qed, A, A, 5, one_pi),
    amplitude("qed_ese_L3", Model::qed, &["e-"], &["e-"], 3, one_pi),
    amplitude("qed_ese_L4", Model::qed, &["e-"], &["e-"], 4, one_pi),
    // QCD: gluons, ghosts and six quark flavours.
    amplitude("qcd_gggg_L1", Model::qcd, GG, GG, 1, plain),
    amplitude("qcd_gggg_L2", Model::qcd, GG, GG, 2, plain),
    amplitude("qcd_gggg_L3", Model::qcd, GG, GG, 3, plain),
    amplitude("qcd_ggggg_L1", Model::qcd, GG, &["g", "g", "g"], 1, plain),
    amplitude("qcd_ggggg_L2", Model::qcd, GG, &["g", "g", "g"], 2, plain),
    amplitude("qcd_uugg_L1", Model::qcd, UU, GG, 1, plain),
    amplitude("qcd_uugg_L2", Model::qcd, UU, GG, 2, plain),
    amplitude("qcd_uugg_L3", Model::qcd, UU, GG, 3, plain),
    amplitude("qcd_uuuu_L2", Model::qcd, UU, UU, 2, plain),
    amplitude("qcd_gse_L3", Model::qcd, G, G, 3, one_pi),
    amplitude("qcd_gse_L4", Model::qcd, G, G, 4, one_pi),
    amplitude("qcd_gse_L5", Model::qcd, G, G, 5, one_pi),
    amplitude("qcd_gse_all_L3", Model::qcd, G, G, 3, plain),
    amplitude("qcd_tad_L2", Model::qcd, NONE, G, 2, tadpoles),
    amplitude("qcd_tad_L3", Model::qcd, NONE, G, 3, tadpoles),
    amplitude("qcd_vac_L2", Model::qcd, NONE, NONE, 2, tadpoles),
    amplitude("qcd_vac_L3", Model::qcd, NONE, NONE, 3, tadpoles),
    amplitude("qcd_vac_L4", Model::qcd, NONE, NONE, 4, tadpoles),
    amplitude("qcd_vac1pi_L4", Model::qcd, NONE, NONE, 4, one_pi_tadpoles),
    amplitude("qcd_vac1pi_L5", Model::qcd, NONE, NONE, 5, one_pi_tadpoles),
    amplitude("qcd_gggg_sym_L1", Model::qcd, GG, GG, 1, symmetrized),
    amplitude(
        "qcd_gggg_identical_L1",
        Model::qcd,
        GG,
        GG,
        1,
        identical_numerators,
    ),
    amplitude(
        "qcd_gggg_sign_L1",
        Model::qcd,
        GG,
        GG,
        1,
        numerators_up_to_sign,
    ),
    amplitude(
        "qcd_gggg_scalar_L1",
        Model::qcd,
        GG,
        GG,
        1,
        numerators_up_to_scalar,
    ),
    amplitude(
        "qcd_gggg_ghost_filter_L1",
        Model::qcd,
        GG,
        GG,
        1,
        ghosts_filtered,
    ),
    amplitude(
        "qcd_gggg_ghost_veto_L1",
        Model::qcd,
        GG,
        GG,
        1,
        ghosts_vetoed,
    ),
    // Yang-Mills: gluons and ghosts.
    amplitude("ym_g6_L0", Model::yang_mills, GG, &["g"; 4], 0, plain),
    amplitude("ym_g7_L0", Model::yang_mills, GG, &["g"; 5], 0, plain),
    amplitude("ym_g8_L0", Model::yang_mills, GG, &["g"; 6], 0, plain),
    amplitude("ym_gse_L4", Model::yang_mills, G, G, 4, one_pi),
    amplitude("ym_gse_L5", Model::yang_mills, G, G, 5, one_pi),
    amplitude("ym_vac_L5", Model::yang_mills, NONE, NONE, 5, plain),
    // The Standard Model with the full CKM matrix.
    amplitude("sm_eemm_L0", Model::standard_model, EE, MUMU, 0, plain),
    amplitude("sm_eemm_L1", Model::standard_model, EE, MUMU, 1, plain),
    amplitude("sm_eemm_L2", Model::standard_model, EE, MUMU, 2, plain),
    amplitude(
        "sm_eemm_identical_L1",
        Model::standard_model,
        EE,
        MUMU,
        1,
        identical_numerators,
    ),
    amplitude(
        "sm_eemm_sign_L1",
        Model::standard_model,
        EE,
        MUMU,
        1,
        numerators_up_to_sign,
    ),
    amplitude(
        "sm_eemm_scalar_L1",
        Model::standard_model,
        EE,
        MUMU,
        1,
        numerators_up_to_scalar,
    ),
    amplitude(
        "sm_ee4f_L0",
        Model::standard_model,
        EE,
        &["mu-", "vm~", "u", "d~"],
        0,
        plain,
    ),
    amplitude(
        "sm_ee4f_L1",
        Model::standard_model,
        EE,
        &["mu-", "vm~", "u", "d~"],
        1,
        plain,
    ),
    amplitude(
        "sm_eeww_L0",
        Model::standard_model,
        EE,
        &["W+", "W-"],
        0,
        plain,
    ),
    amplitude(
        "sm_eeww_L1",
        Model::standard_model,
        EE,
        &["W+", "W-"],
        1,
        plain,
    ),
    amplitude(
        "sm_ggtt_L1",
        Model::standard_model,
        GG,
        &["t", "t~"],
        1,
        plain,
    ),
    amplitude(
        "sm_ggttg_L0",
        Model::standard_model,
        GG,
        &["t", "t~", "g"],
        0,
        plain,
    ),
    amplitude(
        "sm_ggttg_L1",
        Model::standard_model,
        GG,
        &["t", "t~", "g"],
        1,
        plain,
    ),
    amplitude(
        "sm_hse_L1",
        Model::standard_model,
        &["H"],
        &["H"],
        1,
        one_pi,
    ),
    amplitude(
        "sm_hse_L2",
        Model::standard_model,
        &["H"],
        &["H"],
        2,
        one_pi,
    ),
    amplitude("sm_ase_L2", Model::standard_model, A, A, 2, one_pi),
    amplitude(
        "sm_eemm_qcd0_L1",
        Model::standard_model,
        EE,
        MUMU,
        1,
        no_qcd,
    ),
    amplitude(
        "sm_eemm_qcd0_L2",
        Model::standard_model,
        EE,
        MUMU,
        2,
        no_qcd,
    ),
    amplitude(
        "sm_uudd_qed2_L1",
        Model::standard_model,
        UU,
        DD,
        1,
        qed_at_most_two,
    ),
    amplitude(
        "sm_uudd_qed2_L2",
        Model::standard_model,
        UU,
        DD,
        2,
        qed_at_most_two,
    ),
    // Several orders, lower bounds, a tree level that the orders prove empty, fermion-loop ranges
    // and a user filter rejecting fermions that the range counts.
    amplitude(
        "sm_uudd_qcd2_qed2_L1",
        Model::standard_model,
        UU,
        DD,
        1,
        qcd2_qed2,
    ),
    amplitude(
        "sm_uudd_qcd3_L1",
        Model::standard_model,
        UU,
        DD,
        1,
        qcd_at_least_three,
    ),
    amplitude(
        "sm_eemm_qed4_L0to1",
        Model::standard_model,
        EE,
        MUMU,
        1,
        qed4_up_to_one_loop,
    ),
    amplitude(
        "sm_eemm_floop0_L1",
        Model::standard_model,
        EE,
        MUMU,
        1,
        no_fermion_loops,
    ),
    amplitude(
        "qcd_gggg_floop1_L1",
        Model::qcd,
        GG,
        GG,
        1,
        one_fermion_loop,
    ),
    amplitude(
        "qcd_gse_all_floop1_dfilter_L2",
        Model::qcd,
        G,
        G,
        2,
        one_fermion_loop_without_down_quarks,
    ),
    amplitude(
        "sm_eemm_floop1_tfilter_L1",
        Model::standard_model,
        EE,
        MUMU,
        1,
        one_fermion_loop_without_top_quarks,
    ),
    amplitude(
        "qed_ase_floop2_L3",
        Model::qed,
        A,
        A,
        3,
        two_fermion_loops_one_pi,
    ),
    amplitude("ew_eemm_L1", Model::electroweak, EE, MUMU, 1, plain),
    amplitude("ew_eemm_L2", Model::electroweak, EE, MUMU, 2, plain),
    // Vertex rules with identical field content.
    amplitude("ct_phi3_2pt_L1", ct_phi3, PHI, PHI, 1, plain),
    amplitude("ct_qcd_gggg_L1", ct_qcd, GG, GG, 1, plain),
    amplitude("ct_qcd_uugg_L1", ct_qcd, UU, GG, 1, plain),
    // Their vacua and tadpoles have automorphisms that move vertices with duplicated rules.
    amplitude("ct_phi3_2pt_tad_L2", ct_phi3, PHI, PHI, 2, tadpoles),
    amplitude("ct_phi3_vac_L2", ct_phi3, NONE, NONE, 2, tadpoles),
    amplitude("ct_phi3_vac_L3", ct_phi3, NONE, NONE, 3, tadpoles),
    amplitude("ct_phi3_vac_L4", ct_phi3, NONE, NONE, 4, tadpoles),
    amplitude("ct_qcd_vac_L2", ct_qcd, NONE, NONE, 2, tadpoles),
    amplitude("ct_qcd_vac_L3", ct_qcd, NONE, NONE, 3, tadpoles),
    // Counterterm copies with their own order: the rule choice decides the order.
    amplitude(
        "ctorder_qcd_gggg_L1",
        ct_order_qcd,
        GG,
        GG,
        1,
        one_counterterm,
    ),
    amplitude(
        "ctorder_qcd_uugg_L1",
        ct_order_qcd,
        UU,
        GG,
        1,
        one_counterterm,
    ),
    amplitude(
        "ctorder_qcd_vac_L3",
        ct_order_qcd,
        NONE,
        NONE,
        3,
        at_most_one_counterterm,
    ),
    amplitude(
        "ctorder_phi3_vac_L4",
        ct_order_phi3,
        NONE,
        NONE,
        4,
        two_counterterms,
    ),
    // Cross sections.
    cross_section("xs_qed_eemm_L1", Model::qed, EE, MUMU, 1, plain),
    cross_section("xs_qed_eemm_L2", Model::qed, EE, MUMU, 2, plain),
    cross_section("xs_qed_eemm_L3", Model::qed, EE, MUMU, 3, plain),
    cross_section(
        "xs_qed_eemm_identical_L2",
        Model::qed,
        EE,
        MUMU,
        2,
        identical_numerators,
    ),
    cross_section(
        "xs_qed_eemm_sign_L2",
        Model::qed,
        EE,
        MUMU,
        2,
        numerators_up_to_sign,
    ),
    cross_section(
        "xs_qed_eemm_scalar_L2",
        Model::qed,
        EE,
        MUMU,
        2,
        numerators_up_to_scalar,
    ),
    cross_section("xs_qed_eedd_L1", Model::qed, EE, DD, 1, plain),
    cross_section("xs_qed_eedd_L2", Model::qed, EE, DD, 2, plain),
    cross_section("xs_qed_eeaa_L2", Model::qed, EE, &["a", "a"], 2, plain),
    // Identical beams under the beam exchange, the mirror identification and both.
    cross_section("xs_qed_aaee_L2", Model::qed, &["a", "a"], EE, 2, plain),
    cross_section(
        "xs_qed_aaee_sym_L2",
        Model::qed,
        &["a", "a"],
        EE,
        2,
        beam_exchange,
    ),
    cross_section(
        "xs_qed_aaee_lr_L2",
        Model::qed,
        &["a", "a"],
        EE,
        2,
        left_right,
    ),
    cross_section(
        "xs_qed_aaee_sym_lr_L2",
        Model::qed,
        &["a", "a"],
        EE,
        2,
        beam_exchange_left_right,
    ),
    cross_section("xs_qed_eeee_L1", Model::qed, EE, EE, 1, plain),
    cross_section("xs_qed_eeee_L2", Model::qed, EE, EE, 2, plain),
    cross_section("xs_qed_eeee_1pi_L1", Model::qed, EE, EE, 1, loop_one_pi),
    cross_section("xs_qed_eeee_1pi_L2", Model::qed, EE, EE, 2, loop_one_pi),
    cross_section("xs_sm_eemm_L1", Model::standard_model, EE, MUMU, 1, plain),
    cross_section("xs_sm_eemm_L2", Model::standard_model, EE, MUMU, 2, plain),
    cross_section("xs_sm_uudd_L1", Model::standard_model, UU, DD, 1, plain),
    cross_section("xs_qcd_gggg_L1", Model::qcd, GG, GG, 1, plain),
    cross_section("xs_qcd_gggg_sym_L1", Model::qcd, GG, GG, 1, beam_exchange),
    cross_section("xs_qcd_gggg_L2", Model::qcd, GG, GG, 2, plain),
    cross_section("xs_qcd_gguu_L1", Model::qcd, GG, UU, 1, plain),
    cross_section("xs_qcd_gguu_L2", Model::qcd, GG, UU, 2, plain),
    cross_section("xs_qcd_gguu_sym_L1", Model::qcd, GG, UU, 1, beam_exchange),
    cross_section("xs_qcd_gguu_1pi_L1", Model::qcd, GG, UU, 1, loop_one_pi),
];
