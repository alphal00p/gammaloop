use std::collections::BTreeMap;

use feynkit_generator::{GenerationFilter, GenerationOptions, Generator, Process};
use feynkit_model::Model;
use idenso::shorthands::{metric::MetricSimplifier, schoonschip::Schoonschip};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol},
    parser::ParseSettings,
    symbol,
};

fn is_ufo_symbol(symbol: Symbol, name: &str) -> bool {
    symbol.get_stripped_name() == name
        && (symbol.get_namespace() == "UFO" || symbol.get_namespace().starts_with("UFO::"))
}

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
        .with_graph_filter(GenerationFilter::ParticleVeto(vec![
            1, -1, 2, -2, 3, -3, 4, -4, 6, -6,
        ]));
    let generated = Generator::new(model).generate(&process, &options).unwrap();
    let gluon_loop = generated
        .diagrams
        .iter()
        .find(|diagram| diagram.edges().all(|(_, _, edge)| edge.particle.pdg == 21))
        .expect("the QCD self-energy sample must contain its gluon loop");
    let numerator = Atom::parse(
        gluon_loop.numerator().unwrap(),
        "feynkit_generator_qcd_numerator_test",
        ParseSettings::default(),
    )
    .unwrap();
    let momentum = symbol!("FeynKit::Momentum");
    assert!(momentum.has_tag("spenso::tensor"));
    assert!(momentum.has_tag("spenso::rank1"));

    let tensorial = numerator.replace_map(|term, _, out| {
        let AtomView::Fun(function) = term else {
            return;
        };
        let arguments = function.iter().collect::<Vec<_>>();
        if function.get_symbol() == momentum {
            if let [edge, index] = arguments.as_slice() {
                **out = FunctionBuilder::new(momentum)
                    .add_arg(*edge)
                    .add_arg(spenso::mink!(4, index.to_owned()))
                    .finish();
            }
        } else if is_ufo_symbol(function.get_symbol(), "Metric")
            && let [left, right] = arguments.as_slice()
        {
            **out = spenso::g!(
                spenso::mink!(4, left.to_owned()),
                spenso::mink!(4, right.to_owned())
            );
        }
    });
    let contracted = tensorial
        .expand()
        .simplify_metrics()
        .to_dots()
        .to_plain_string();

    assert!(contracted.contains("spenso::dot("), "{contracted}");
    assert!(contracted.contains("FeynKit::Momentum"), "{contracted}");
}
