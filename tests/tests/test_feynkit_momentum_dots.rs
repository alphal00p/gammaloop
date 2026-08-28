use std::collections::BTreeMap;

use feynkit_generator::{
    GenerationFilter, GenerationOptions, Generator, ParticleSelector, Process,
};
use feynkit_model::Model;
use idenso::shorthands::{metric::MetricSimplifier, schoonschip::Schoonschip};
use symbolica::{atom::AtomCore, symbol};

#[test]
fn generated_qcd_numerator_contractions_convert_to_dots() {
    let model = Model::from_json(include_str!(
        "../../crates/feynkit-model/tests/fixtures/sm.json"
    ))
    .unwrap();
    let process = Process::amplitude(["g"], ["g"])
        .with_loop_count(1, 1)
        .unwrap();
    let options = GenerationOptions::default()
        .max_vertices(2)
        .with_graph_filter(GenerationFilter::CouplingOrders(BTreeMap::from([
            ("QCD".to_owned(), (2, Some(2))),
            ("QED".to_owned(), (0, Some(0))),
        ])))
        .with_graph_filter(GenerationFilter::ParticleVeto(
            [1, -1, 2, -2, 3, -3, 4, -4, 6, -6]
                .into_iter()
                .map(ParticleSelector::from)
                .collect(),
        ));
    let generated = Generator::new(model).generate(&process, &options).unwrap();
    let gluon_loop = generated
        .diagrams
        .iter()
        .find(|diagram| {
            diagram.edges().all(|(_, _, edge)| {
                diagram
                    .model()
                    .particle_by_id(edge.particle)
                    .is_ok_and(|particle| particle.pdg_code == 21)
            })
        })
        .expect("the QCD self-energy sample must contain its gluon loop");
    let numerator = gluon_loop.numerator().clone();
    let momentum = symbol!("FeynKit::Momentum");
    assert!(momentum.has_tag("spenso::tensor"));
    assert!(momentum.has_tag("spenso::rank1"));

    let contracted = numerator
        .expand()
        .simplify_metrics()
        .to_dots()
        .to_plain_string();

    assert!(contracted.contains("spenso::dot("), "{contracted}");
    assert!(contracted.contains("FeynKit::Momentum"), "{contracted}");
}
