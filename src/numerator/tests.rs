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

    graph.generate_numerator(&model);

    let mut numerator = graph.derived_data.numerator.clone().unwrap();
    // println!("{}", numerator);
    for (lhs, rhs) in graph.generate_lmb_replacement_rules() {
        numerator =
            lhs.into_pattern()
                .replace_all(numerator.as_view(), &rhs.into_pattern(), None, None);
    }

    let should_be = Atom::parse("
        16/81*ee^4*(MT*id(euc(4,173),euc(4,224))+(P0(lor(4,236))+P1(lor(4,236))+P2(lor(4,236))+K0(lor(4,236)))
            *γ(lor(4,236),bis(4,173),bis(4,224)))
        *(MT*id(euc(4,212),euc(4,172))+(P1(lor(4,233))+P2(lor(4,233))+K0(lor(4,233)))
            *γ(lor(4,233),bis(4,212),bis(4,172)))
        *(MT*id(euc(4,219),euc(4,211))+(P2(lor(4,234))+K0(lor(4,234)))
            *γ(lor(4,234),bis(4,219),bis(4,211)))
        *(MT*id(euc(4,225),euc(4,218))+γ(lor(4,235),bis(4,225),bis(4,218))*K0(lor(4,235)))
        *complex(0,1)^4
        *id(coaf(3,172),cof(3,173))
        *id(coaf(3,211),cof(3,212))
        *id(coaf(3,218),cof(3,219))
        *id(coaf(3,224),cof(3,225))
        *ϵ0(lor(4,174))*ϵ̅1(lor(4,213))*ϵ̅2(lor(4,220))*ϵ̅3(lor(4,226))
        *γ(lor(4,174),bis(4,173),bis(4,172))*γ(lor(4,213),bis(4,212),bis(4,211))*γ(lor(4,220),bis(4,219),bis(4,218))*γ(lor(4,226),bis(4,225),bis(4,224))").unwrap();

    assert_eq!(numerator, should_be);

    // if let AtomView::Mul(mul) = numerator.as_view() {
    //     let net = SymbolicTensor::mul_to_tracking_network(mul).unwrap();

    //     println!("{}", net.dot());
    //     for (_, t) in net.graph.nodes {
    //         // println!("{}", t.structure());
    //     }
    // }
}
