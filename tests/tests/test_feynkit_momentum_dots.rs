use std::collections::BTreeMap;

use feynkit_generator::{GenerationFilter, GenerationOptions, ParticleSelector, Process};
use feynkit_model::Model;
use idenso::tensor::ContractSettings;
use symbolica::atom::AtomCore;

#[test]
fn generated_qcd_numerator_contractions_convert_to_dots() {
    let model = Model::from_json(include_str!(
        "../../crates/feynkit-model/tests/fixtures/sm.json"
    ))
    .unwrap();
    let process = Process::new(["g"], ["g"]);
    let options = GenerationOptions::default()
        .with_loop_count(1, 1)
        .unwrap()
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
    let generated = process.generate_diagrams(model, &options).unwrap();
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
    let momentum = feynkit_graph::momentum_symbol();
    assert!(momentum.has_tag("spenso::tensor"));
    assert!(momentum.has_tag("spenso::rank1"));

    let contracted = numerator
        .contract(ContractSettings::default())
        .unwrap()
        .to_dots()
        .unwrap()
        .into_expression();

    assert!(
        contracted.contains_symbol(spenso::network::tags::SPENSO_TAG.dot),
        "{contracted}"
    );
    assert!(contracted.contains_symbol(momentum), "{contracted}");
}
