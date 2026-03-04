use gammalooprs::{
    dot,
    graph::{Graph, parse::IntoGraph},
    initialisation::test_initialise,
    model::Model,
    processes::Amplitude,
    utils::{W_, test_utils::load_generic_model},
    uv::{
        UVgenerationSettings,
        settings::{MATADSettings, VakintSettings},
    },
};
use idenso::color::{CS, ColorSimplifier};
use spenso::shadowing::symbolica_utils::AtomCoreExt;
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    parse, parse_lit,
};

#[test]
fn scalar_pole_part() {
    test_initialise().unwrap();
    let sunrise: Vec<Graph> = dot!( digraph sunrise{
        edge [particle=scalar_1]
        A -> B    [ id=0]
        A -> B     [ id=1]
        A -> B   [ id=2]

    },"scalars")
    .unwrap();

    let mut amp = Amplitude::from_graph_list("bub", sunrise).unwrap();

    let model = load_generic_model("scalars");

    let a = amp.graphs[0]
        .renormalization_part(&UVgenerationSettings {
            softct: false,
            only_integrated: true,
            vakint: VakintSettings {
                normalization: "MSbar".to_string(),
                additional_normalization: "1".to_string(),
                ..Default::default()
            },
            ..Default::default()
        })
        .unwrap();

    println!("ren part: {:>}", a);
    println!(
        "ren part: {:>}",
        model.apply_parameter_replacement_rules(
            &model.apply_coupling_replacement_rules(&a.simplify_color().expand())
        )
    );
}

pub fn align_to_rqft(atom: &Atom, model: &Model) -> Atom {
    (model
        .apply_parameter_replacement_rules(
            &model.apply_coupling_replacement_rules(&atom.simplify_color().expand()),
        )
        .replace(parse_lit!(spenso::mink(gammalooprs::dim)))
        .with(parse_lit!(spenso::mink(4)))
        .collect_factors()
        .replace(CS.tr)
        .with(Atom::num((1, 2)))
        .replace(CS.nc)
        .with(CS.ca)
        .replace(parse!("UFO::aS"))
        .with(parse!("gs").npow(2) / (Atom::var(Symbol::PI) * 4))
        / -8)
        .expand_num()
}
#[test]
fn finite_part_ghost_3loop() {
    test_initialise().unwrap();

    let model = load_generic_model("sm");
    let g: Vec<Graph> = Graph::from_path(
        "./integration-tests/resources/graphs/rqft_ghost_3l.dot",
        &model,
    )
    .unwrap();

    let mut amp = Amplitude::from_graph_list("bub", g).unwrap();

    let settings = UVgenerationSettings {
        softct: false,
        only_integrated: true,
        vakint: VakintSettings {
            normalization: "(
                𝑖*(𝜋^((4-2*eps)/2))
             * (exp(-EulerGamma))^(eps)
             * (exp(-logmUVmu-log_mu_sq))^(eps)
             )^(-n_loops)"
                .to_string(),
            // evaluation_methods: vec!["matad".to_string()],
            // normalization: "(exp(log_mu_sq+logmUVmu)/(4*𝜋*exp(-EulerGamma)))^(eps*n_loops)"
            //     .to_string(),
            // normalization: "FMFTandMATAD".to_string(),
            additional_normalization: "1".to_string(),
            matad: MATADSettings {
                expand_masters: true,
                susbstitute_masters: true,
                substitute_hpls: false,
                direct_numerical_substition: false,
            },
            ..Default::default()
        },
        ..Default::default()
    };

    let a = amp.graphs[0].renormalization_part(&settings).unwrap();
    //p1.p1*gs^6*ca^3*rat( - 3/8*ep^-2 + 29/32*ep^-1)
    insta::assert_snapshot!(amp.graphs[0].graph.name,@"d1");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-1/4*ca^2*ε^2*𝜋^2+-27/4+-29/48*ca^2*ε^2+-3/2*ca^2+-87/32*ε^2+-9/8*ε^2*𝜋^2+11/12*ca^2*ε+33/8*ε)*ca*dot(P(0),P(0),mink(4))*gs^6*ε^(-3)"
    );

    //p1.p1*gs^6*ca^3*rat( - 1/16*ep^-2 + 5/192*ep^-1)
    let a = amp.graphs[1].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[1].graph.name,@"d2");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-5𝑖/192*ε+1𝑖/16)*ca^3*dot(P(0),P(0),mink(4))*gs^6*ε^(-2)"
    );

    //p1.p1*gs^6*ca^3*rat(9/128*ep^-3 - 39/256*ep^-2 + 9/128*ep^-1)
    let a = amp.graphs[2].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[2].graph.name,@"d3");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-9𝑖/128+-9𝑖/128*ε^2+39𝑖/256*ε)*ca^3*dot(P(0),P(0),mink(4))*gs^6*ε^(-3)"
    );

    //p1.p1*gs^6*ca^3*rat(9/128*ep^-3 - 39/256*ep^-2 + 27/128*ep^-1)
    let a = amp.graphs[3].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[3].graph.name,@"d4");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-21𝑖/128*ε+9𝑖/64)*ca^3*dot(P(0),P(0),mink(4))*gs^6*ε^(-2)"
    );

    //p1.p1*i_*gs^4*ca^2*rat( - 5/8*ep^-2 + 35/48*ep^-1)
    let a = amp.graphs[4].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"d5");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-27𝑖/128*ε^2+-9𝑖/128+39𝑖/256*ε)*ca^3*dot(P(0),P(0),mink(4))*gs^6*ε^(-3)"
    ); //-1/2 * target

    //p1.p1*i_*gs^4*ca^2*rat(1/24*ep^-1)
    let a = amp.graphs[5].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[5].graph.name,@"d6");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-131𝑖/128*ε+159𝑖/64+1𝑖*ε^2+27𝑖/64*ε^2*𝜋^2)*ca^3*dot(P(0),P(0),mink(4))*gs^6*ε^(-3)"
    );

    let a = amp.graphs[6].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"d5");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-35𝑖/128*ε+-3𝑖/64+1𝑖*ε^2)*ca^3*dot(P(0),P(0),mink(4))*gs^6*ε^(-3)"
    );

    let a = amp.graphs[7].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"d5");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-27𝑖/32+63𝑖/64*ε+99𝑖/128*ε^2)*ca^3*dot(P(0),P(0),mink(4))*gs^6*ε^(-3)"
    );

    let a = amp.graphs[8].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"d5");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"0"
    );

    let a = amp.graphs[9].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"d5");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"0"
    );

    let a = amp.graphs[10].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[11].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[12].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[13].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[14].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[15].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[16].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[17].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[18].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[19].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[20].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[21].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[22].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[23].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[24].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[25].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[26].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[27].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[28].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[29].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[30].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[31].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[32].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[33].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[34].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[35].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[36].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[37].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[38].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[39].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[40].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[41].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[42].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[43].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[44].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[45].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[46].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[47].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[48].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[49].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[50].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[51].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[52].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[53].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[54].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[55].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[56].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[57].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[58].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[59].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[60].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[61].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[62].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[63].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[64].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[65].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[66].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[67].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[68].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[69].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[70].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[71].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[72].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[73].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[74].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[75].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[76].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );

    let a = amp.graphs[77].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );
}

#[test]
fn finite_part_ghost_2loop() {
    test_initialise().unwrap();
    let mut g: Vec<Graph> = dot!(digraph d1 {//0
      overall_factor= "+1"
      num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))"
      in1 [style=invis];
      in1 -> v1:0
      [particle="ghG" pin="x:@-left"];

      out1 [style=invis];
      v2:1 -> out1
      [particle="ghG" pin="x:@+right"];

      v1 -> v3 [particle= "g"];
      v1 -> v4 [particle= "ghG"];
      v2 -> v3 [particle= "g"];
      v4 -> v2 [particle= "ghG"];
      v3 -> v4 [particle= "g"];
    }
    digraph d2 {//1
      overall_factor= "+1"
      num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))"

      in1 [style=invis];
      in1 -> v1:0
      [particle="ghG" pin="x:@-left"];

      out1 [style=invis];
      v2:1 -> out1
      [particle="ghG" pin="x:@+right"];

      v1 -> v3 [particle= "g"];
      v1 -> v4 [particle= "ghG"];
      v3 -> v2 [particle= "ghG"];
      v2 -> v4 [particle= "g"];
      v4 -> v3 [particle= "ghG"];
    }
    digraph d3 {//2
      overall_factor= "+1"
      num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))"

      in1 [style=invis];
      in1 -> v1:0
      [particle="ghG" pin="x:@-left"];

      out1 [style=invis];
      v2:1 -> out1
      [particle="ghG" pin="x:@+right"];

      v1 -> v2 [particle= "g"];
      v1 -> v3 [particle= "ghG"];
      v4 -> v2 [particle= "ghG"];
      v3 -> v4 [particle= "g"];
      v3 -> v4 [particle= "ghG"];
    }
    digraph d4 {//3
      overall_factor= "-1"
      num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))"
      in1 [style=invis];
      in1 -> v1:0
      [particle="ghG" pin="x:@-left"];

      out1 [style=invis];
      v2:1 -> out1
      [particle="ghG" pin="x:@+right"];

      v1 -> v2 [particle= "ghG"];
      v1 -> v3 [particle= "g"];
      v2 -> v4 [particle= "g"];
      v3 -> v4 [particle= "d"];
      v4 -> v3 [particle= "d"];
    }
    digraph d5 {//4
      overall_factor= "+1/2"
      num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))"
      in1 [style=invis];
      in1 -> v1:0
      [particle="ghG" pin="x:@-left"];

      out1 [style=invis];
      v2:1 -> out1
      [particle="ghG" pin="x:@+right"];

      v1 -> v2 [particle= "ghG"];
      v1 -> v3 [particle= "g"];
      v2 -> v4 [particle= "g"];
      v3 -> v4 [particle= "g"];
      v3 -> v4 [particle= "g"];
    }
    digraph d6 {//5
      overall_factor= "-1"
      num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))"
      in1 [style=invis];
      in1 -> v1:0
      [particle="ghG" pin="x:@-left"];

      out1 [style=invis];
      v2:1 -> out1
      [particle="ghG" pin="x:@+right"];

      v1 -> v2 [particle= "ghG"];
      v1 -> v3 [particle= "g"];
      v2 -> v4 [particle= "g"];
      v3 -> v4 [particle= "ghG"];
      v4 -> v3 [particle= "ghG"];
    }

    digraph GL00{//6
        num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))"
        exte0	 [style=invis];
        exte0	-> 3:0	 [id=0 dir=none particle="ghG" pin="x:@-left"];
        exte1	 [style=invis];
        2:1	-> exte1	 [id=1 dir=none particle="ghG" pin="x:@+right"];
        1:2	-> 0:3	 [id=2 particle="c"];
        0:4	-> 1:5	 [id=3 particle="c"];
        0:6	-> 3:7	 [id=4 particle="g"];
        1:8	-> 2:9	 [id=5 particle="g"];
        3:10	-> 2:11	 [id=6 particle="ghG"];
    }

    digraph GL01{//7
        num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))"
        exte0	 [style=invis];
        exte0	-> 3:0	 [id=0 dir=none particle="ghG" pin="x:@-left"];
        exte1	 [style=invis];
        2:1	-> exte1	 [id=1 dir=none particle="ghG" pin="x:@+right"];
        1:2	-> 0:3	 [id=2  particle="t"];
        0:4	-> 1:5	 [id=3 particle="t"];
        0:6	-> 3:7	 [id=4 dir=none particle="g"];
        1:8	-> 2:9	 [id=5 dir=none particle="g"];
        3:10	-> 2:11	 [id=6 dir=none particle="ghG"];
    }

    digraph GL02{//8
        num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))"
        exte0	 [style=invis];
        exte0	-> 3:0	 [id=0 dir=none particle="ghG" pin="x:@-left"];
        exte1	 [style=invis];
        2:1	-> exte1	 [id=1 dir=none particle="ghG" pin="x:@+right"];
        0:2	-> 1:3	 [id=2 dir=none particle="g"];
        0:4	-> 1:5	 [id=3 dir=none particle="ghG"];
        3:6	-> 0:7	 [id=4 dir=none particle="ghG"];
        1:8	-> 2:9	 [id=5 dir=none particle="ghG"];
        2:10-> 3:11	 [id=6 dir=none source=2 particle="g"];
    }

    digraph GL03{//9
        num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))"
           exte0	 [style=invis];
        exte0	-> 3:0	 [id=0 dir=none particle="ghG" pin="x:@-left"];
        exte1	 [style=invis];
        2:1	-> exte1	 [id=1 dir=none particle="ghG" pin="x:@+right"];
        1:2	-> 0:3	 [id=2 dir=none particle="ghG"];
        0:4	-> 1:5	 [id=3 dir=none particle="ghG"];
        0:6	-> 3:7	 [id=4 dir=none particle="g"];
        1:8	-> 2:9	 [id=5 dir=none particle="g"];
        3:10	-> 2:11	 [id=6 dir=none particle="ghG"];
    }

    digraph GL04{//10
        num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))"
           exte0	 [style=invis];
        exte0	-> 3:0	 [id=0 dir=none particle="ghG" pin="x:@-left"];
        exte1	 [style=invis];
        2:1	-> exte1	 [id=1 dir=none  particle="ghG" pin="x:@+right"];
        1:2	-> 0:3	 [id=2 dir=none particle="ghG"];
        0:4	-> 2:5	 [id=3 dir=none particle="ghG"];
        0:6	-> 3:7	 [id=4 dir=none particle="g"];
        1:8	-> 2:9	 [id=5 dir=none particle="g"];
        3:10	-> 1:11	 [id=6 dir=none particle="ghG"];
    }


    digraph GL05{//11
        num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))"
           exte0	 [style=invis];
        exte0	-> 2:0	 [id=0 dir=none particle="ghG" pin="x:@-left"];
        exte1	 [style=invis];
        1:1	-> exte1	 [id=1 dir=none particle="ghG" pin="x:@+right"];
        0:2	-> 1:3	 [id=2 dir=none particle="ghG"];
        2:4	-> 0:5	 [id=3 dir=none particle="ghG"];
        0:6	-> 3:7	 [id=4 dir=none particle="g"];
        1:8	-> 3:9	 [id=5 dir=none particle="g"];
        2:10	-> 3:11	 [id=6 dir=none particle="g"];
    }

    digraph GL06{//12
        num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))"
          exte0	 [style=invis];
        exte0	-> 1:0	 [id=0 dir=none  particle="ghG" pin="x:@-left"];
        exte1	 [style=invis];
        0:1	-> exte1	 [id=1 dir=none particle="ghG" pin="x:@+right"];
        1:2	-> 0:3	 [id=2 dir=none  particle="ghG"];
        0:4	-> 3:5	 [id=3 dir=none  particle="g"];
        1:6	-> 2:7	 [id=4 dir=none  particle="g"];
        3:8	-> 2:9	 [id=5 particle="u"];
        2:10	-> 3:11	 [id=6  particle="u"];
    }

    digraph GL07{//13
        num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))"
        exte0	 [style=invis];
        exte0	-> 1:0	 [id=0 dir=none  particle="ghG" pin="x:@-left"];
           exte1	 [style=invis];
        0:1	-> exte1	 [id=1 dir=none particle="ghG" pin="x:@+right"];
        1:2	-> 0:3	 [id=2 dir=none   particle="ghG"];
        0:4	-> 3:5	 [id=3 dir=none    particle="g"];
        1:6	-> 2:7	 [id=4 dir=none    particle="g"];
        2:8	-> 3:9	 [id=5 dir=none   particle="g"];
        2:10	-> 3:11	 [id=6 dir=none  particle="g"];
    }
    digraph GL10{//14
        num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))"
           exte0	 [style=invis];
        exte0	-> 3:0	 [id=0 dir=none  particle="ghG" pin="x:@-left"];
        exte1	 [style=invis];
        2:1	-> exte1	 [id=1 dir=none particle="ghG" pin="x:@+right"];
        1:2	-> 0:3	 [id=2  particle="b"];
        0:4	-> 1:5	 [id=3  particle="b"];
        0:6	-> 3:7	 [id=4 dir=none particle="g"];
        1:8	-> 2:9	 [id=5 dir=none particle="g"];
        3:10	-> 2:11	 [id=6 dir=none   particle="ghG"];
    })
    .unwrap();

    let mut amp = Amplitude::from_graph_list("bub", g).unwrap();

    let settings = UVgenerationSettings {
        softct: false,
        only_integrated: true,
        vakint: VakintSettings {
            normalization: "(
            𝑖*(𝜋^((4-2*eps)/2))
         * (exp(-EulerGamma))^(eps)
         * (exp(-logmUVmu-log_mu_sq))^(eps)
         )^(-n_loops)"
                .to_string(),
            // normalization: "(exp(log_mu_sq+logmUVmu)/(4*𝜋*exp(-EulerGamma)))^(eps*n_loops)"
            //     .to_string(),
            // normalization: "FMFTandMATAD".to_string(),
            additional_normalization: "1".to_string(),
            ..Default::default()
        },
        ..Default::default()
    };
    let model = load_generic_model("sm");

    let a = amp.graphs[0].renormalization_part(&settings).unwrap();
    //p1.p1*i_*gs^4*ca^2*rat( - 3/16*ep^-2 + 5/32*ep^-1)
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-3/16+5/32*ε)*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-2)"
    );

    //p1.p1*i_*gs^4*ca^2*rat( - 1/16*ep^-2 + 1/32*ep^-1)
    let a = amp.graphs[1].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-1/16+1/32*ε)*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-2)"
    ); //-1 * target

    //p1.p1*i_*gs^4*ca^2*rat( - 1/8*ep^-2 + 1/16*ep^-1)
    let a = amp.graphs[2].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-1/8+1/16*ε)*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-2)"
    ); //-1 * target

    //p1.p1*i_*gs^4*ca*nf*rat(1/4*ep^-2 - 5/24*ep^-1)
    let a = amp.graphs[3].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-5/8*ε+3/4)*dot(P(0),P(0),mink(4))*gs^4*ε^(-2)"
    );

    //p1.p1*i_*gs^4*ca^2*rat( - 5/8*ep^-2 + 35/48*ep^-1)
    let a = amp.graphs[4].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-5/8+35/48*ε)*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-2)"
    ); //-1/2 * target

    //p1.p1*i_*gs^4*ca^2*rat(1/24*ep^-1)
    let a = amp.graphs[5].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
    );
}

#[test]
fn finit_part_ghlo() {
    test_initialise().unwrap();

    let g: Vec<Graph> = dot!(digraph d1 {

        projector = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))"

          in1 [style=invis];
          in1 -> v1:0
          [id=0 particle="ghG" pin="x:@-left"];

          out1 [style=invis];
          v2:1 -> out1
          [id=1 particle="ghG" pin="x:@+right"];

          v1:2 -> v2:3 [id=2 particle= "g"];
          v2:4 -> v1:5 [id=3 particle= "ghG~"];
        }
    )
    .unwrap();

    let mut amp = Amplitude::from_graph_list("bub", g).unwrap();

    let model = load_generic_model("sm");

    let a = amp.graphs[0]
        .renormalization_part(&UVgenerationSettings {
            softct: false,
            only_integrated: true,
            ..Default::default()
        })
        .unwrap();

    println!("ren part: {:>}", a);
    insta::assert_snapshot!(
        align_to_rqft(&a,&model)
        .to_bare_ordered_string(),@"-1𝑖/2*ca*dot(P(0),P(0),mink(4))*gs^2*ε^(-1)*𝜋^2"
    );
}
