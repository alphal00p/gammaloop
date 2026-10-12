use feynkit_generator::{
    Loops, Naming, NumeratorOptions, Orders, Process, Restrict, Run, TopologyRules,
};
use feynkit_model::Model;
use idenso::tensor::{ContractSettings, SymbolicTensor};
use symbolica::{atom::AtomCore, symbol};

#[test]
fn generated_qcd_numerator_contractions_convert_to_dots() {
    let model = Model::from_json(include_str!(
        "../../crates/feynkit-model/tests/fixtures/sm.json"
    ))
    .unwrap();
    let request = Process::new(["g"], ["g"])
        .amplitude(Loops::total(1..=1))
        .with_topology(
            TopologyRules::new()
                .with_self_loops(false)
                .with_max_vertices(2),
        )
        .with_orders(Orders::new().exact("QCD", 2).exact("QED", 0))
        .with_restrict(Restrict::new().exclude_particles([1, -1, 2, -2, 3, -3, 4, -4, 6, -6]));
    let run = Run::default();
    let generated = request
        .generate(model.into(), &run)
        .unwrap()
        .build(
            &NumeratorOptions::default().with_naming(Naming::prefix("D")),
            &run,
        )
        .unwrap();
    let gluon_loop = generated
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

    // Indexed vectors contract explicitly; to_dots only normalizes the resulting notation.
    let contracted = SymbolicTensor::infer(numerator)
        .unwrap()
        .contract(ContractSettings {
            collect_chains: false,
            collect_traces: false,
            ..Default::default()
        })
        .unwrap()
        .to_dots()
        .unwrap()
        .into_expression();

    assert!(
        contracted.contains_symbol(symbol!("spenso::dot")),
        "{contracted}"
    );
    assert!(contracted.contains_symbol(momentum), "{contracted}");
}
