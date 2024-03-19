use symbolica::representations::AtomView;

use crate::{
    tensor::{SymbolicTensor, TensorStructure},
    tests_from_pytest::load_amplitude_output,
};

#[test]

fn lbl() {
    let (model, amplitude) = load_amplitude_output("./src/numerator/lbl/");
    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    graph.generate_numerator(&model);
    let numerator = graph.derived_data.numerator.unwrap();

    println!("{}", numerator);

    if let AtomView::Mul(mul) = numerator.as_view() {
        let mut net = SymbolicTensor::mul_to_tracking_network(mul).unwrap();

        println!("{}", net.dot());
        for (id, t) in net.graph.nodes {
            println!("{}", t.structure());
        }
    }
}
