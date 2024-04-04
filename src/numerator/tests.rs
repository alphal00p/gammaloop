use symbolica::printer::PrintOptions;

use crate::tests_from_pytest::load_amplitude_output;

#[test]
#[allow(unused)]
fn lbl() {
    let (model, amplitude) = load_amplitude_output("./src/test_resources/lbl/");

    model
        .export_coupling_replacement_rules("./ignore/", PrintOptions::mathematica())
        .unwrap();
    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    graph.generate_numerator(&model);

    let mut numerator = graph.derived_data.numerator.clone().unwrap();
    // println!("{}", numerator);
    for (lhs, rhs) in graph.generate_lmb_replacement_rules() {
        numerator =
            lhs.into_pattern()
                .replace_all(numerator.as_view(), &rhs.into_pattern(), None, None);
    }

    println!("{}", numerator);

    // if let AtomView::Mul(mul) = numerator.as_view() {
    //     let net = SymbolicTensor::mul_to_tracking_network(mul).unwrap();

    //     println!("{}", net.dot());
    //     for (_, t) in net.graph.nodes {
    //         // println!("{}", t.structure());
    //     }
    // }
}
