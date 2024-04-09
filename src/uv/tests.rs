use symbolica::{printer::PrintOptions, representations::Atom};

use crate::tests_from_pytest::load_amplitude_output;

#[test]
#[allow(unused)]
fn lbl() {
    let (model, amplitude) = load_amplitude_output("./src/test_resources/lbl/");

    model
        .export_coupling_replacement_rules("./ignore/", PrintOptions::mathematica())
        .unwrap();
    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    graph.generate_uv();

    let uv_graph = graph.derived_data.uv.unwrap();

    println!("{}", uv_graph.0.base_dot());

    // if let AtomView::Mul(mul) = numerator.as_view() {
    //     let net = SymbolicTensor::mul_to_tracking_network(mul).unwrap();

    //     println!("{}", net.dot());
    //     for (_, t) in net.graph.nodes {
    //         // println!("{}", t.structure());
    //     }
    // }
}
