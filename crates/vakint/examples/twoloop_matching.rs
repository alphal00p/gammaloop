use symbolica::try_parse;
use vakint::{EvaluationOrder, Vakint, VakintExpression, VakintSettings};

fn main() {
    let settings = VakintSettings {
        allow_unknown_integrals: false,
        evaluation_order: EvaluationOrder::alphaloop_only(),
        ..VakintSettings::default()
    };
    let vakint = Vakint::new().unwrap();
    vakint.validate_settings(&settings).unwrap();

    //println!("Supported topologies:\n{}", vakint.topologies);

    let input = try_parse!(
        "(vakint::k(11,2)*vakint::k(11,2)+vakint::k(11,77)*vakint::k(22,77)+vakint::k(22,33)*p(42,33))*vakint::topo(\
                    vakint::prop(9,vakint::edge(7,10),vakint::k(11),mUVsqA,1)*\
                    vakint::prop(33,vakint::edge(7,10),vakint::k(22),mUVsqA,2)*\
                    vakint::prop(55,vakint::edge(7,10),vakint::k(11)+vakint::k(22),mUVsqA,1)\
                )+\
                (2*vakint::k(33,2)*vakint::k(33,2)+17*vakint::k(44,33)*p(17,33))*vakint::topo(\
                    vakint::prop(7,vakint::edge(9,21),vakint::k(33),mUVsqB,1)*\
                    vakint::prop(13,vakint::edge(9,21),vakint::k(44),mUVsqB,2)*\
                    vakint::prop(17,vakint::edge(9,21),vakint::k(33)+vakint::k(44),mUVsqB,1)\
                )"
    )
    .unwrap();

    let output = vakint
        .to_canonical(&settings, input.as_view(), true)
        .unwrap();

    println!(
        "\nInput:\n\n{}\n\nhas been matched to\n\n{}\n",
        VakintExpression::try_from(input).unwrap(),
        VakintExpression::try_from(output).unwrap()
    );
}
